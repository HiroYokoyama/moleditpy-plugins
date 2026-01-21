PLUGIN_NAME = "Compound Info Report"
PLUGIN_VERSION = "2026.01.21"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Generate a compound info report with properties, adducts, and structure. Useful for organic synthesis experiments."
PLUGIN_ID = "compound_info_report"

import sys
import json

try:
    import urllib.request
    import urllib.parse
    URLLIB_AVAILABLE = True
except ImportError:
    URLLIB_AVAILABLE = False


from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QTextBrowser, 
    QCheckBox, QLabel, QMessageBox, QProgressDialog, QApplication, QFileDialog
)
from PyQt6.QtGui import QPixmap, QImage, QPainter, QTextDocument
from PyQt6.QtPrintSupport import QPrinter, QPrintDialog
from PyQt6.QtCore import Qt, QSize, QRectF, QCoreApplication, QByteArray, QBuffer, QIODevice

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import Descriptors
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None

class PubChemFetcher:
    """Helper to fetch data from PubChem."""
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    @staticmethod
    def get_synonyms(inchikey):
        """Fetch synonyms (CAS, etc.) for an InChIKey."""
        if not URLLIB_AVAILABLE or not inchikey:
            return []
        
        url = f"{PubChemFetcher.BASE_URL}/compound/inchikey/{inchikey}/synonyms/JSON"
        try:
            with urllib.request.urlopen(url) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode('utf-8'))
                    information = data.get("InformationList", {}).get("Information", [])
                    if information:
                        return information[0].get("Synonym", [])
        except Exception:
            pass
        return []


    @staticmethod
    def get_cid(inchikey):
        """Fetch CID for InChIKey."""
        if not URLLIB_AVAILABLE or not inchikey: return None
        try:
            cid_url = f"{PubChemFetcher.BASE_URL}/compound/inchikey/{inchikey}/cids/JSON"
            with urllib.request.urlopen(cid_url) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode('utf-8'))
                    cids = data.get("IdentifierList", {}).get("CID", [])
                    if cids:
                        return cids[0]
        except:
            pass
        return None

    @staticmethod
    def fetch_experimental_properties(cid):
        """Fetch Density and PhysDesc from Experimental Properties View."""
        if not URLLIB_AVAILABLE or not cid: return None, None
        
        density = None
        phys_desc = None
        
        # Query "Chemical and Physical Properties" section
        # We can fetch the whole section which includes Experimental Properties
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading=Chemical+and+Physical+Properties"
        
        try:
            with urllib.request.urlopen(url) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode('utf-8'))
                    
                    # Traversal
                    # Data -> Record -> Section -> [Chemical and Physical Properties] -> Section -> [Experimental Properties] -> Section -> [Property Name]
                    record = data.get("Record", {})
                    sections = record.get("Section", [])
                    
                    # Since we filtered by heading, the top sections should be relevant or the one we asked for
                    for sec in sections:
                        if sec.get("TOCHeading") == "Chemical and Physical Properties":
                            for sub in sec.get("Section", []):
                                if sub.get("TOCHeading") == "Experimental Properties":
                                    for prop in sub.get("Section", []):
                                        heading = prop.get("TOCHeading", "")
                                        
                                        # Density
                                        if heading == "Density" and not density:
                                            info = prop.get("Information", [])
                                            if info:
                                                val = info[0].get("Value", {})
                                                if "StringWithMarkup" in val:
                                                    density = val["StringWithMarkup"][0].get("String")
                                                elif "Number" in val:
                                                    density = str(val["Number"][0])
                                        
                                        # Physical Description / Color / Form
                                        if ("Color" in heading or "Physical Description" in heading) and not phys_desc:
                                            info = prop.get("Information", [])
                                            for info_item in info:
                                                val = info_item.get("Value", {})
                                                text = ""
                                                if "StringWithMarkup" in val:
                                                    text = val["StringWithMarkup"][0].get("String")
                                                elif "Number" in val:
                                                    text = str(val["Number"][0])
                                                
                                                if text:
                                                     # Heuristic check
                                                     text_lower = text.lower()
                                                     if any(x in text_lower for x in ["white", "colorless", "yellow", "red", "blue", "green", "solid", "powder", "crystal", "liquid", "oil"]):
                                                         phys_desc = text
                                                         break
                                        
                                        if density and phys_desc: break
                                if density and phys_desc: break
                        if density and phys_desc: break

        except:
            pass
            
        return density, phys_desc

    @staticmethod
    def extract_cas(synonyms):
        """Extract potential CAS numbers from synonyms list."""
        cas_list = []
        for syn in synonyms:
            # Simple heuristic for CAS format: ddd-dd-d
            # Filter out long names, keep short numeric-ish strings
            if len(syn) < 15 and syn.replace('-', '').isdigit():
                 # Check hyphen placement roughly (2 hyphens)
                 if syn.count('-') == 2:
                     cas_list.append(syn)
        # Return top 5 unique ones
        return sorted(list(set(cas_list)))[:5]

class ReportDialog(QDialog):
    def __init__(self, mol, parent=None):
        super().__init__(parent)
        self.mol = mol
        self.setWindowTitle("Compound Info Report")
        self.resize(700, 800)
        self.fetched_data = {}
        self.setup_ui()
        self.generate_initial_report()

    def setup_ui(self):
        layout = QVBoxLayout(self)

        # Controls
        ctrl_layout = QHBoxLayout()
        self.chk_pubchem = QCheckBox("Fetch Online Data (CAS via PubChem)")
        
        if not URLLIB_AVAILABLE:
            self.chk_pubchem.setChecked(False)
            self.chk_pubchem.setEnabled(False)
            self.chk_pubchem.setText("Fetch Online Data (Disabled - urllib missing)")
        else:
            self.chk_pubchem.setChecked(False)
            self.chk_pubchem.stateChanged.connect(self.generate_report)
        
        ctrl_layout.addWidget(self.chk_pubchem)
        ctrl_layout.addStretch()
        layout.addLayout(ctrl_layout)

        # Preview
        self.preview = QTextBrowser()
        self.preview.setOpenExternalLinks(True)
        layout.addWidget(self.preview)

        # Buttons
        btn_layout = QHBoxLayout()
        
        self.btn_print = QPushButton("Print...")
        self.btn_print.clicked.connect(self.print_report)
        
        self.btn_pdf = QPushButton("Save PDF...")
        self.btn_pdf.clicked.connect(self.save_pdf)
        
        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.accept)

        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_print)
        btn_layout.addWidget(self.btn_pdf)
        btn_layout.addWidget(self.btn_close)
        
        layout.addLayout(btn_layout)


    def generate_initial_report(self):
        self.generate_report()

    def calculate_adducts(self, exact_mass):
        # Constants (IUPAC/NIST standard atomic weights for monoisotopic mass)
        MASS_ELECTRON = 0.00054858
        MASS_H = 1.00782503
        MASS_NA = 22.98976928
        MASS_K  = 38.96370668
        MASS_CL = 34.96885268

        # [M+H]+ : Add H, remove electron
        m_h = exact_mass + MASS_H - MASS_ELECTRON  # = exact_mass + 1.007276

        # [M+Na]+ : Add Na, remove electron
        m_na = exact_mass + MASS_NA - MASS_ELECTRON # = exact_mass + 22.989221

        # [M+K]+ : Add K, remove electron
        m_k = exact_mass + MASS_K - MASS_ELECTRON   # = exact_mass + 38.963158

        # [M+Cl]- : Add Cl, add electron (Negative mode)
        m_cl = exact_mass + MASS_CL + MASS_ELECTRON # = exact_mass + 34.969401

        # [M-H]- : Remove H, add electron (net: remove proton)
        # Calculation: M - H + e_ = M - (H - e_) = M - Proton
        m_minus_h = exact_mass - (MASS_H - MASS_ELECTRON) # = exact_mass - 1.007276

        return [
            ("[M+H]<sup>+</sup>", f"{m_h:.4f}"),
            ("[M+Na]<sup>+</sup>", f"{m_na:.4f}"),
            ("[M+K]<sup>+</sup>", f"{m_k:.4f}"),
            ("[M-H]<sup>-</sup>", f"{m_minus_h:.4f}"),
            ("[M+Cl]<sup>-</sup>", f"{m_cl:.4f}"),
        ]


    def capture_scene_image(self):
        """Capture the molecule image from the Main Window scene."""
        mol_b64 = None
        img_w, img_h = 0, 0
        
        # Try capturing from parent scene (User's 2D View)
        if self.parent() and hasattr(self.parent(), 'scene'):
            try:
                scene = self.parent().scene
                
                # Calculate tight bounding box
                molecule_bounds = QRectF()
                for item in scene.items():
                    item_type = type(item).__name__
                    if item_type in ["AtomItem", "BondItem"] and item.isVisible():
                        molecule_bounds = molecule_bounds.united(item.sceneBoundingRect())
                
                if molecule_bounds.isEmpty():
                    molecule_bounds = scene.itemsBoundingRect()
                
                if not molecule_bounds.isEmpty():
                    # Add padding
                    padding = 20
                    source_rect = molecule_bounds.adjusted(-padding, -padding, padding, padding)
                    
                    w = int(source_rect.width())
                    h = int(source_rect.height())
                    
                    if w > 0 and h > 0:
                        img = QImage(w, h, QImage.Format.Format_ARGB32)
                        img.fill(Qt.GlobalColor.white)
                        
                        painter = QPainter(img)
                        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
                        scene.render(painter, QRectF(0, 0, w, h), source_rect)
                        painter.end()
                        
                        ba = QByteArray()
                        buf = QBuffer(ba)
                        buf.open(QIODevice.OpenModeFlag.WriteOnly)
                        img.save(buf, "PNG")
                        mol_b64 = ba.toBase64().data().decode("utf-8")
                        img_w, img_h = w, h
                        
            except Exception as e:
                print(f"Scene Capture Error: {e}")
        
        # Fallback to RDKit if scene capture failed
        if not mol_b64 and Chem:
            try:
                 # Standardize fallback size
                 img = Draw.MolToImage(self.mol, size=(400, 300))
                 from io import BytesIO
                 import base64
                 buffered = BytesIO()
                 img.save(buffered, format="PNG")
                 mol_b64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
                 img_w, img_h = 400, 300
            except:
                 pass
                 
        return mol_b64, img_w, img_h


    def generate_html(self, for_pdf=False, img_b64=None, img_w=0, img_h=0):
        # Gather Data
        formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
        mw = Descriptors.MolWt(self.mol)
        exact_mass = Descriptors.ExactMolWt(self.mol)
        logp = Descriptors.MolLogP(self.mol)

        # PubChem Data
        cas_numbers = []
        common_name = ""
        density = None
        phys_desc = None
        



    def fetch_pubchem_data(self, progress_callback=None):
        data = {
            "common_name": "",
            "density": None,
            "phys_desc": None,
            "cas_numbers": []
        }
        
        if not self.chk_pubchem.isChecked():
            return data

        try:
            inchikey = Chem.MolToInchiKey(self.mol)
            
            # Step 1: Synonyms
            if progress_callback: progress_callback("Fetching Synonyms...", 1)
            synonyms = PubChemFetcher.get_synonyms(inchikey)
            data["cas_numbers"] = PubChemFetcher.extract_cas(synonyms)
            if synonyms:
                data["common_name"] = synonyms[0]
            
            # Step 2: Properties (Density & Phys Desc)
            if progress_callback: progress_callback("Fetching Properties...", 2)
            
            cid = PubChemFetcher.get_cid(inchikey)
            if cid:
                if progress_callback: progress_callback("Fetching Experimental Data...", 3)
                density, phys_desc = PubChemFetcher.fetch_experimental_properties(cid)
                data["density"] = density
                data["phys_desc"] = phys_desc
            
        except Exception:
            pass
            
        return data

    def build_html(self, pubchem_data, for_pdf=False, img_b64=None, img_w=0, img_h=0):
        # Gather Basic Data
        formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
        
        # Format Formula with Subscripts
        import re
        formula_html = re.sub(r'(\d+)', r'<sub>\1</sub>', formula)
        
        mw = Descriptors.MolWt(self.mol)
        exact_mass = Descriptors.ExactMolWt(self.mol)
        logp = Descriptors.MolLogP(self.mol)
        
        adducts = self.calculate_adducts(exact_mass)
        
        # --- HTML Components ---
        
        # 1. Molecule Info Table
        mol_info_rows = f"""
            <tr><td>Formula</td><td class="prop-val">{formula_html}</td></tr>
            <tr><td>Mol. Weight</td><td class="prop-val">{mw:.4f}</td></tr>
            <tr><td>Exact Mass</td><td class="prop-val">{exact_mass:.4f}</td></tr>
            <tr><td>LogP</td><td class="prop-val">{logp:.2f}</td></tr>
        """
        mol_info_table = f"""
            <h3>Molecule Information</h3>
            <table>
                <tr><th width="40%">Property</th><th>Value</th></tr>
                {mol_info_rows}
            </table>
        """

        # 2. Adducts Table
        adduct_rows = ''.join([f'<tr><td>{name}</td><td>{mass}</td></tr>' for name, mass in adducts])
        adduct_table = f"""
            <h3>Adduct Ions</h3>
            <table>
                <tr><th width="40%">Ion Type</th><th>m/z</th></tr>
                {adduct_rows}
            </table>
        """

        # 3. PubChem Table
        pubchem_html = ""
        pd = pubchem_data
        if pd["common_name"] or pd["density"] or pd.get("cas_numbers") or pd["phys_desc"]:
            rows = []
            if pd["common_name"]: rows.append(f'<tr><td>Common Name</td><td>{pd["common_name"]}</td></tr>')
            if pd["phys_desc"]: rows.append(f'<tr><td>Physical State/Color</td><td>{pd["phys_desc"]}</td></tr>')
            if pd["density"]: rows.append(f'<tr><td>Density</td><td>{pd["density"]}</td></tr>')
            if pd.get("cas_numbers"): rows.append(f'<tr><td>CAS Numbers</td><td>{", ".join(pd["cas_numbers"])}</td></tr>')
            
            pubchem_html = f"""
            <h3>Online Data (PubChem)</h3>
            <table>
                <tr><th width="30%">Property</th><th>Value</th></tr>
                {''.join(rows)}
            </table>
            """

        # --- Layout Assembly ---
        

        # PDF Image Insertion
        img_html = ""
        if for_pdf and img_b64:

            # Calculate text-document friendly scaling
            # Target Box (The "Frame")
            MAX_W = 400
            MAX_H = 250
            
            disp_w = img_w
            disp_h = img_h
            
            # Scale down if needed to fit MAX_W x MAX_H
            if disp_w > 0 and disp_h > 0:
                scale_w = MAX_W / disp_w
                scale_h = MAX_H / disp_h
                scale = min(scale_w, scale_h)

                # Limit max magnification
                scale = min(scale, 0.5)
                
                if scale < 1.0:
                    disp_w = int(disp_w * scale)
                    disp_h = int(disp_h * scale)
                elif scale > 1.0:
                     disp_w = int(disp_w * scale)
                     disp_h = int(disp_h * scale)

            # Ensure strict integer for HTML
            disp_w = int(disp_w)
            disp_h = int(disp_h)

            img_html = f"""
            <h3>Structure</h3>
            <div style="text-align: center; margin-bottom: 20px;">
                 <img src="data:image/png;base64,{img_b64}" width="{disp_w}" height="{disp_h}" />
            </div>
            """


        # Side-by-Side Table Layout (50% / 50%)
        # Side-by-Side Table Layout (50% / 50%)
        # Use HTML attributes for column widths width="50%".
        
        # Ensure tables inside take full available width of their cell
        mol_info_table = mol_info_table.replace("<table>", '<table width="100%">')
        adduct_table = adduct_table.replace("<table>", '<table width="100%">')

        side_by_side = f"""
        <table width="100%" border="0" cellpadding="0" cellspacing="0">
            <tr>
                <td width="50%" valign="top" style="padding-right: 15px;">
                    {mol_info_table}
                </td>
                <td width="50%" valign="top" style="padding-left: 15px;">
                    {adduct_table}
                </td>
            </tr>
        </table>
        """

        html = f"""
        <html>
        <head>
            <style>
                body {{ font-family: sans-serif; padding: 30px; }}
                h1 {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; text-align: center; margin-bottom: 30px; }}
                h3 {{ color: #2c3e50; margin-top: 0px; margin-bottom: 8px; border-bottom: 1px solid #eee; padding-bottom: 4px; font-size: 1.1em; }}
                table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
                th, td {{ text-align: left; padding: 8px; border-bottom: 1px solid #ddd; font-size: 0.95em; }}
                th {{ background-color: #f2f2f2; }}
                .prop-val {{ font-weight: bold; }}
            </style>
        </head>
        <body>
            <h1>Compound Info Report</h1>
            
            {img_html}
            
            {side_by_side}
            
            {pubchem_html}
            
            <div style="margin-top: 50px; border-top: 1px solid #ddd; padding-top: 10px; text-align: right; color: gray; font-size: 0.8em;">
                Generated by MoleditPy Compound Info Report Plugin (Version {PLUGIN_VERSION})
            </div>
        </body>
        </html>
        """
        return html

    def generate_report(self):
        # Progress Dialog for Fetching
        pubchem_data = {
            "common_name": "",
            "density": None,
            "phys_desc": None,
            "cas_numbers": []
        }

        if self.chk_pubchem.isChecked():
            progress = QProgressDialog("Fetching PubChem Data...", "Cancel", 0, 4, self)
            progress.setWindowModality(Qt.WindowModality.WindowModal)
            progress.setMinimumDuration(0) # Show immediately
            progress.setValue(0)
            
            def update_progress(msg, val):
                progress.setLabelText(msg)
                progress.setValue(val)
                QApplication.processEvents()
            
            pubchem_data = self.fetch_pubchem_data(update_progress)
            progress.setValue(4)
        
        html = self.build_html(pubchem_data, for_pdf=False)
        self.preview.setHtml(html)
        
        # Cache data
        self.last_pubchem_data = pubchem_data

    def print_report(self):
        printer = QPrinter(QPrinter.PrinterMode.HighResolution)
        dialog = QPrintDialog(printer, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.preview.print(printer)

    def save_pdf(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Save PDF", "synthesis_report.pdf", "PDF Files (*.pdf)")
        if filename:
            # Progress for PDF
            progress = QProgressDialog("Generating PDF...", "Cancel", 0, 3, self)
            progress.setWindowModality(Qt.WindowModality.WindowModal)
            progress.setMinimumDuration(0)
            progress.setValue(0)
            QApplication.processEvents()

            try:
                # 1. Capture Image
                progress.setLabelText("Capturing Molecule Image...")
                img_b64, w, h = self.capture_scene_image()
                progress.setValue(1)
                QApplication.processEvents()

                # 2. Build HTML
                if not hasattr(self, 'last_pubchem_data'):
                     self.last_pubchem_data = self.fetch_pubchem_data()
                
                # Check for cancel
                if progress.wasCanceled(): return

                progress.setLabelText("Building PDF Document...")
                html = self.build_html(self.last_pubchem_data, for_pdf=True, img_b64=img_b64, img_w=w, img_h=h)
                progress.setValue(2)
                QApplication.processEvents()
                
                # 3. Print
                printer = QPrinter(QPrinter.PrinterMode.HighResolution)
                printer.setOutputFormat(QPrinter.OutputFormat.PdfFormat)
                printer.setOutputFileName(filename)
                
                from PyQt6.QtGui import QTextDocument
                doc = QTextDocument()
                doc.setHtml(html)
                doc.print(printer)
                
                progress.setValue(3)
                QMessageBox.information(self, "Export", f"Saved to {filename}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save PDF: {e}")
            finally:
                progress.close()


def show_report(context):
    """Callback for Analysis Menu"""
    mw = context.get_main_window()
    
    if not Chem:
        QMessageBox.warning(mw, "Error", "RDKit is not installed.")
        return

    # Get Molecule
    mol = None
    if hasattr(mw, 'current_mol') and mw.current_mol:
        mol = mw.current_mol
    elif hasattr(mw, 'data') and hasattr(mw.data, 'to_rdkit_mol'):
        try:
            mol = mw.data.to_rdkit_mol()
        except:
            pass
            
    if not mol:
        QMessageBox.warning(mw, "Error", "Please load a molecule first.")
        return

    # Generate Structure from SMILES if 2D coords not present (optional safety)
    try:
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
    except:
        pass

    dialog = ReportDialog(mol, mw)
    dialog.exec()

def initialize(context):
    """Plugin Initialization Hook"""
    context.add_analysis_tool("Compound Info Report", lambda: show_report(context))
