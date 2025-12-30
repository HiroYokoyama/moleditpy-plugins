"""
PubChem Structure Identifier Plugin (Cleaned - No SMILES)
Resolve chemical names and fetch basic molecular properties via PubChem.
"""
import sys
import json
import urllib.request
import urllib.parse
import ssl
from PyQt6.QtWidgets import (
    QInputDialog, QMessageBox, QProgressDialog, QDialog, 
    QVBoxLayout, QPushButton, QTextBrowser, QApplication, QHBoxLayout
)
from PyQt6.QtCore import Qt

# --- 【追加】truststoreのインポート ---
try:
    import truststore
except ImportError:
    truststore = None

try:
    from rdkit import Chem
except ImportError:
    Chem = None

# --- Metadata ---
PLUGIN_NAME = "PubChem Structure Identifier (truststore)"
PLUGIN_VERSION = "2025.12.30"
PLUGIN_AUTHOR = "HiroYokoyama"
# ... (中略) ...

class PubChemResolver:
    """Helper class for PubChem API interactions"""
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    # --- 【修正】最強のSSLコンテキスト作成メソッド ---
    @staticmethod
    def _create_ssl_context():
        """
        Create an SSL context that:
        1. Uses the system trust store (fixes security software/proxy issues).
        2. Lowers security level (fixes Basic Constraints errors).
        """
        ctx = None
        
        # 1. truststoreが使えるなら使い、Windowsの証明書情報を読み込む
        if truststore:
            try:
                ctx = truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
            except Exception:
                # 念のため失敗したら標準に戻す
                ctx = ssl.create_default_context()
        else:
            # インストールされていない場合は標準のものを使う
            ctx = ssl.create_default_context()

        # 2. OpenSSL 3.0対策: セキュリティレベルを1に下げる
        # これにより "Basic Constraints" エラーを回避する
        try:
            ctx.set_ciphers('DEFAULT@SECLEVEL=1')
        except (ssl.SSLError, AttributeError):
            pass
            
        return ctx

    @staticmethod
    def resolve_name_to_smiles(name):
        """
        Resolve a chemical name to a SMILES string (Required for Import).
        Returns: (smiles, error_message)
        """
        if not name:
            return None, "Empty name provided."

        try:
            encoded_name = urllib.parse.quote(name)
            url = f"{PubChemResolver.BASE_URL}/compound/name/{encoded_name}/property/IsomericSMILES/JSON"
            
            # 【修正】context引数を追加
            with urllib.request.urlopen(url, context=PubChemResolver._create_ssl_context()) as response:
                if response.status != 200:
                   return None, f"HTTP Error: {response.status}"
                
                data = json.loads(response.read().decode('utf-8'))
                properties = data.get("PropertyTable", {}).get("Properties", [])
                if properties:
                    smiles = properties[0].get("IsomericSMILES")
                    if smiles:
                        return smiles, None
                
                return None, "No SMILES found for this name."

        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, "Molecule not found in PubChem."
            return None, f"HTTP Error: {e.code}"
        except Exception as e:
            return None, f"Network/Parsing Error: {str(e)}"

    @staticmethod
    def get_compound_details(inchikey):
        """
        Fetch detailed properties for an InChIKey (Excluding SMILES).
        Returns: (details_dict, error_message)
        """
        if not inchikey:
            return None, "Empty InChIKey provided."

        # Initialize dictionary with the query key
        details = {"InChIKey": inchikey}
        
        try:
            # 1. Get Common Name (Title)
            desc_url = f"{PubChemResolver.BASE_URL}/compound/inchikey/{inchikey}/description/JSON"
            try:
                # 【修正】context引数を追加
                with urllib.request.urlopen(desc_url, context=PubChemResolver._create_ssl_context()) as response:
                    if response.status == 200:
                        data = json.loads(response.read().decode('utf-8'))
                        info_list = data.get("InformationList", {}).get("Information", [])
                        for info in info_list:
                            title = info.get("Title")
                            if title:
                                details["Common Name"] = title
                                break
            except:
                details["Common Name"] = "Unknown"

            # 2. Get Physical Properties (No SMILES, No XLogP, No TPSA)
            props_to_fetch = "MolecularFormula,MolecularWeight,IUPACName"
            prop_url = f"{PubChemResolver.BASE_URL}/compound/inchikey/{inchikey}/property/{props_to_fetch}/JSON"

            # 【修正】context引数を追加
            with urllib.request.urlopen(prop_url, context=PubChemResolver._create_ssl_context()) as response:
                if response.status != 200:
                     return None, f"HTTP Error: {response.status}"
                
                data = json.loads(response.read().decode('utf-8'))
                props = data.get("PropertyTable", {}).get("Properties", [])
                
                if props:
                    p = props[0]
                    if p.get("MolecularFormula"): details["Formula"] = p.get("MolecularFormula")
                    if p.get("MolecularWeight"): details["Mol. Weight"] = p.get("MolecularWeight")
                    if p.get("IUPACName"): details["IUPAC Name"] = p.get("IUPACName")
                    
                    return details, None
                
                return None, "No properties found for this InChIKey."

        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, "InChIKey not found in PubChem."
            return None, f"HTTP Error: {e.code}"
        except Exception as e:
            return None, f"Network/Parsing Error: {str(e)}"

class MoleculeDetailsDialog(QDialog):
    """
    Dialog to display fetched molecule details and allow copying data.
    """
    def __init__(self, details, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Molecule Details - PubChem")
        self.resize(500, 300)
        self.details = details
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)

        self.info_browser = QTextBrowser()
        self.info_browser.setOpenExternalLinks(True)
        
        # Build HTML dynamically
        html_parts = [
            "<style>",
            "th { text-align: left; color: #555; padding-right: 15px; white-space: nowrap; }",
            "td { padding-bottom: 5px; width: 100%; }",
            "h3 { margin-bottom: 5px; color: #2c3e50; }",
            "</style>",
            f"<h3>{self.details.get('Common Name', 'Molecule Details')}</h3>",
            '<table width="100%">'
        ]

        # Table Rows
        if 'IUPAC Name' in self.details:
            html_parts.append(f"<tr><th>IUPAC Name:</th><td>{self.details['IUPAC Name']}</td></tr>")
        if 'Formula' in self.details:
            html_parts.append(f"<tr><th>Formula:</th><td><b>{self.details['Formula']}</b></td></tr>")
        if 'Mol. Weight' in self.details:
            html_parts.append(f"<tr><th>Mol. Weight:</th><td>{self.details['Mol. Weight']} g/mol</td></tr>")
        if 'InChIKey' in self.details:
            html_parts.append(f"<tr><th>InChIKey:</th><td><small>{self.details['InChIKey']}</small></td></tr>")
        
        html_parts.append("</table>")

        self.info_browser.setHtml("".join(html_parts))
        layout.addWidget(self.info_browser)

        # -- Buttons --
        btn_layout = QHBoxLayout()
        
        self.btn_copy = QPushButton("Copy Info")
        self.btn_copy.setToolTip("Copy properties to clipboard")
        self.btn_copy.clicked.connect(self.copy_to_clipboard)
        
        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.accept)

        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_copy)
        btn_layout.addWidget(self.btn_close)
        
        layout.addLayout(btn_layout)

    def copy_to_clipboard(self):
        """Format details as plain text and copy to clipboard"""
        order = [
            "Common Name", "IUPAC Name", "Formula", "Mol. Weight", "InChIKey"
        ]
        
        lines = []
        for key in order:
            if key in self.details:
                lines.append(f"{key}: {self.details[key]}")
        
        text_data = "\n".join(lines)
        
        clipboard = QApplication.clipboard()
        clipboard.setText(text_data)
        
        original_text = self.btn_copy.text()
        self.btn_copy.setText("Copied!")
        self.btn_copy.setDisabled(True)
        
        from PyQt6.QtCore import QTimer
        QTimer.singleShot(1500, lambda: self.reset_copy_button(original_text))

    def reset_copy_button(self, text):
        self.btn_copy.setText(text)
        self.btn_copy.setEnabled(True)

def resolve_and_load(plugin_context):
    """Callback for 'Import from PubChem'"""
    mw = plugin_context.main_window
    
    name, ok = QInputDialog.getText(
        mw, "Import from PubChem", "Enter Chemical Name (e.g., Aspirin, Benzene):"
    )
    
    if not ok or not name.strip():
        return

    progress = QProgressDialog("Searching PubChem...", "Cancel", 0, 0, mw)
    progress.setWindowModality(Qt.WindowModality.WindowModal)
    progress.show()
    
    smiles, error = PubChemResolver.resolve_name_to_smiles(name.strip())
    
    progress.close()
    
    if error:
        QMessageBox.warning(mw, "PubChem Error", f"Could not resolve '{name}':\n{error}")
        return
        
    if not smiles:
        QMessageBox.warning(mw, "PubChem Error", f"No structure found for '{name}'.")
        return

    if hasattr(mw, 'main_window_string_importers'):
        try:
            mw.main_window_string_importers.load_from_smiles(smiles)
            mw.statusBar().showMessage(f"Loaded '{name}' from PubChem.", 5000)
            
            if hasattr(mw, 'plotter_controller'):
                 mw.plotter_controller.reset_camera()
                 
        except Exception as e:
            QMessageBox.critical(mw, "Import Error", f"Failed to load SMILES:\n{e}")
    else:
        QMessageBox.critical(mw, "Error", "String Importer module not found in Main Window.")


def identify_current_molecule(plugin_context):
    """Callback for 'Identify Molecule' action"""
    mw = plugin_context.main_window
    
    if not Chem:
        QMessageBox.warning(mw, "Error", "RDKit is required for this feature.")
        return

    # Get Current Molecule
    mol = None
    if hasattr(mw, 'current_mol') and mw.current_mol:
        mol = mw.current_mol
    elif hasattr(mw, 'data') and hasattr(mw.data, 'to_rdkit_mol'):
        try:
            mol = mw.data.to_rdkit_mol()
        except:
            pass
            
    if not mol:
        QMessageBox.warning(mw, "Error", "No valid molecule loaded to identify.")
        return

    # Calculate InChIKey Locally
    try:
        inchikey = Chem.MolToInchiKey(mol)
    except Exception as e:
        QMessageBox.warning(mw, "Error", f"Failed to calculate InChIKey: {e}")
        return

    # Fetch Details
    progress = QProgressDialog("Querying PubChem details...", "Cancel", 0, 0, mw)
    progress.setWindowModality(Qt.WindowModality.WindowModal)
    progress.show()

    details, error = PubChemResolver.get_compound_details(inchikey)
    
    progress.close()

    if details:
        dlg = MoleculeDetailsDialog(details, mw)
        dlg.exec()
    else:
        QMessageBox.information(mw, "PubChem Result", f"Could not identify molecule.\nError: {error}")


# --- Plugin Entry Points ---

def run(main_window):
    """Entry point for the plugin when launched from the Plugins menu."""
    class Context:
        def __init__(self, mw):
            self.main_window = mw
            
    identify_current_molecule(Context(main_window))