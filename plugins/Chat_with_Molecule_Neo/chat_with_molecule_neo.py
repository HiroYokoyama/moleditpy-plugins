#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PLUGIN_NAME = "Chat with Molecule Neo"
PLUGIN_VERSION = "2025.12.28"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Chat with Google Gemini about the current molecule. Automatically injects SMILES context. (Neo Version)"
PLUGIN_ID = "chat_with_molecule_neo"
"""


import sys
import os
import json
import io
import base64
import re
import urllib.request
import urllib.parse
# Matplotlib for LaTeX rendering (OO interface to avoid thread issues)
try:
    import matplotlib
    matplotlib.use('Agg') # Non-interactive backend
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# RDKit imports (Top-level)
from rdkit import Chem
from rdkit.Chem import AllChem

from PyQt6.QtCore import (
    Qt, QThread, pyqtSignal, QTimer, QSize, QEvent, 
    QPointF, QRunnable, QThreadPool, QObject
)
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QLineEdit,
    QPushButton, QLabel, QFrame, QScrollArea, QSizePolicy,
    QProgressBar, QMessageBox, QApplication, QMainWindow, QMenu,
    QFileDialog, QTextBrowser, QPlainTextEdit, QComboBox, QDialog
)
from PyQt6.QtGui import (
    QTextCursor, QColor, QDesktopServices, QAction, QIcon,
    QFont, QTextBlockFormat, QTextCharFormat, QPainter, QGuiApplication
)


class PubChemResolver:
    """Helper class for PubChem API interactions (Embedded)"""
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    @staticmethod
    def resolve_inchikey_to_name(inchikey):
        """
        Resolve an InChIKey to a chemical name (Title).
        Returns: (name, error_message)
        """
        if not inchikey:
            return None, "Empty InChIKey provided."

        try:
            url = f"{PubChemResolver.BASE_URL}/compound/inchikey/{inchikey}/description/JSON"
            
            with urllib.request.urlopen(url) as response:
                if response.status != 200:
                    return None, f"HTTP Error: {response.status}"
                
                data = json.loads(response.read().decode('utf-8'))
                
                # Parse JSON response
                info_list = data.get("InformationList", {}).get("Information", [])
                
                # Prioritize entries with 'Title'
                for info in info_list:
                    title = info.get("Title")
                    if title:
                        return title, None
                
                return None, "No name found for this InChIKey."

        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, "InChIKey not found in PubChem."
            return None, f"HTTP Error: {e.code}"
        except Exception as e:
            return None, f"Network/Parsing Error: {str(e)}"

# --- Metadata ---
PLUGIN_NAME = "Chat with Molecule Neo"
PLUGIN_VERSION = "2025.12.28"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Chat with Google Gemini about the current molecule. Automatically injects SMILES context. (Neo Version)"
PLUGIN_ID = "chat_with_molecule_neo"

SYSTEM_PROMPT = """You are an expert computational and organic chemistry assistant embedded within the advanced molecular editor software "MoleditPy". 
Your users are researchers, students, or chemistry enthusiasts. Adhere to the following guidelines:

### 1. Protocol for Molecular Identification & Graph Consistency (CRITICAL)
When provided with a SMILES string as `Context` (and optionally a Name), you MUST NOT rely on intuition.
**Graph & Atom Mapping**: The SMILES provided will include **Atom Map Numbers** on HEAVY atoms (e.g., `[C:1]`, `[O:2]`). Hydrogens are implicit.
*   These numbers correspond to the skeletal atoms in the graph.
*   **Reaction SMARTS**: You can still target the implicit hydrogens attached to these mapped atoms.
    *   Example: To fluorinate `[C:1]`, you can write `[C:1][H]>>[C:1][F]`. 
    *   The execution engine will handle the explicit hydrogen math.

**User Selection**: If the Context includes "[User Selection: ...]", you MUST prioritize transformations on those specific atoms.
    *   Example: "Add a methyl group" + Selection `[5, 6]` -> Methylate atoms 5 and 6.

You must perform the following internal reasoning steps before generating a response:
1.  **Stoichiometry Check**: Mentally count the atoms (e.g., Carbon count, Hydrogen count) to verify the molecular formula.
2.  **Topology Analysis**: Analyze the number of rings and their fusion pattern.
3.  **Identification**: If a Name is provided in the context, use it as the primary identification but verify it against the SMILES.

### 2. Tool Use & Function Calling
You have access to specific tools to interact with the software. To use a tool, you **MUST** output **exactly ONE** JSON block at the very END of your response.
**Do NOT** output multiple separate JSON blocks. If you need multiple tools, use a JSON Array.

**Format**:
```json
{
  "tool": "tool_name",
  "params": { ... }
}
```

**Best Practice**: For any operation that modifies or highlights specific atoms (e.g., `highlight_substructure`, `set_electronic_state`), **use `atom_indices`** to specify targets precisely. This avoids SMARTS errors and ensures reliable targeting. 
> [!WARNING]
> Modifying the molecular structure (e.g., `apply_transformation` or `load_molecule`) may **change or reset atom indices**. When performing complex multi-step modifications, consider generating separate responses to ensure targeting stays accurate.

**Available Tools**:
1.  **`apply_transformation`**: Apply a chemical reaction to the current molecule.
    *   `reaction_smarts`: The Reaction SMARTS string defining the transformation.
    *   **Tip**: Do NOT specify hydrogen counts (e.g., `H2`, `H3`) in the PRODUCT side. Let RDKit calculate implicit Hs automatically.
    *   **Carbon Matching**: Use `[#6]` to match ANY Carbon (aliphatic or aromatic). Pattern `C` (uppercase) ONLY matches aliphatic carbons, while `c` (lowercase) ONLY matches aromatic ones. Using `[#6]` is recommended for robust targeting.
    *   Example: `{"tool": "apply_transformation", "params": {"reaction_smarts": "[#6:1][H]>>[#6:1][Cl]"}}` (Chlorinate any carbon)

2.  **`highlight_substructure`**: Visually highlight atoms matching a pattern or indices.
    *   `smarts`: (Optional) The SMARTS pattern to find and highlight.
        *   **Tip**: Use `[#6]` to match ANY Carbon (aliphatic or aromatic). Pattern `C` ONLY matches aliphatic carbons.
    *   `atom_indices`: (Optional) List of Atom Map Numbers (Integers) to highlight e.g. [1, 5]. **STRONGLY RECOMMENDED** to use this for precise targeting if User Selection is available.
    *   Example: `{"tool": "highlight_substructure", "params": {"smarts": "[#8]"}}` (Highlight Hydroxyls/Ethers)

3.  **`calculate_descriptors`**: Calculate molecular properties.
    *   `properties`: List of properties to calculate. Available options:
        - "MW" (Molecular Weight)
        - "LogP" (Partition Coefficient)
        - "TPSA" (Topological Polar Surface Area)
        - "HBondDonor" or "HBD" (Hydrogen Bond Donors)
        - "HBondAcceptor" or "HBA" (Hydrogen Bond Acceptors)
        - "RotatableBonds" or "RB" (Number of Rotatable Bonds)
        - "AromaticRings" (Number of Aromatic Rings)
        - "Rings" (Total Ring Count)
    *   Example: `{"tool": "calculate_descriptors", "params": {"properties": ["MW", "LogP", "TPSA", "HBD", "HBA"]}}`

4.  **`orca_input_generator`**: Create an ORCA calculation input file.
    *   `filename`: Suggested filename (e.g., "[molecule]-opt.inp").
    *   `header`: The content to go BEFORE the coordinates (e.g., `%maxcore... ! B3LYP...`). **Do NOT** include Title, Charge, or Multiplicity/Spin here.
    *   `footer`: The content to go AFTER the coordinates.
    *   Example: `{"tool": "orca_input_generator", "params": {"filename": "opt.inp", "header": "! B3LYP def2-SVP Opt"}}`

5.  **`gaussian_input_generator`**: Create a Gaussian calculation input file.
    *   `filename`: Suggested filename (e.g., "[molecule]-opt.gjf").
    *   `header`: The content to go BEFORE the coordinates (e.g., `%nproc... #P B3LYP...`). **Do NOT** include Title, Charge, or Multiplicity here.
    *   `footer`: The content to go AFTER the coordinates.
    *   Example: `{"tool": "gaussian_input_generator", "params": {"filename": "opt.gjf", "header": "#P B3LYP/6-31G* Opt"}}`

6.  **`load_molecule`**: Load a new molecule from SMILES string (replaces current molecule).
    *   `smiles`: The SMILES string of the molecule to load.
    *   `name`: Optional human-readable name for the molecule.
    *   Example: `{"tool": "load_molecule", "params": {"smiles": "c1ccccc1", "name": "Benzene"}}`

7.  **`convert_to_3d`**: Trigger 3D structure generation and visualization.
    *   **Note**: This tool is RARELY needed. Tools like `generate_input` already perform 3D conversion automatically. Only use this if the user explicitly asks to "show me the 3D structure" or "visualize in 3D".
    *   No parameters required.
    *   Example: `{"tool": "convert_to_3d", "params": {}}`

8.  **`clear_canvas`**: Clear all molecules from the editor (start fresh).
    *   No parameters required.
    *   Example: `{"tool": "clear_canvas", "params": {}}`

9.  **`set_electronic_state`**: Set formal charge or multiplicity on a **single atom**.
    *   `atom_indices`: **ABSOLUTELY REQUIRED**. The Atom Map Number (Integer) e.g. `5`. **WILL FAIL WITHOUT THIS.**
    *   `charge`: Integer formal charge (e.g., 1 for +1, -1 for -1). **Do NOT use for charge 0** - that is the default.
    *   `multiplicity`: Integer multiplicity.
    *   **CRITICAL**: You MUST ALWAYS include `atom_indices`. Without it, the tool does NOTHING.
    *   Example: `{"tool": "set_electronic_state", "params": {"atom_indices": 1, "charge": 1}}` (Make Atom 1 +1)

10. **`save_file`**: Save text content to a file with any extension.
    *   `filename`: The filename (e.g., "script.py", "data.csv", "config.ini").
    *   `content`: The text content to write.
    *   **Tags**:
        *   `[[atom]]`: Replaced with "ElementSymbol X Y Z" lines (e.g., "C 0.0 0.0 0.0").
        *   `[[atom_count]]`: Replaced with the total number of atoms.
    *   Example: `{"tool": "save_file", "params": {"filename": "mol.xyz", "content": "[[atom_count]]\nTitle\n[[atom]]"}}`

**Multiple Operations**: You can propose multiple tool calls in a single response by using a JSON array **inside the single JSON block**. Operations will execute in order.
Example: `[{"tool": "apply_transformation", "params": {...}}, {"tool": "convert_to_3d", "params": {}}]`

### 3. Response Style
* **Tone**: Professional, intellectual, and helpful.
* **Formatting**:
    * Use **SIMPLE** LaTeX format for chemical formulas and math (e.g., $C_{16}H_{10}$).
    * **Avoid** complex commands like `\\xrightarrow` or `\\stackrel`.
    * **CRITICAL**: Always format ANY molecular structure you mention *in plain text* as a clickable link like this: `[Name](smiles:SMILES_STRING)`.
    * **IMPORTANT**: Do **NOT** put SMILES strings inside LaTeX equations.

### 4. MoleditPy Software Context
* **Role:** You are the integrated AI for **MoleditPy**.
* **Native Features:** Suggest built-in tools for immediate modeling and export.
* **Plugin Ecosystem:** Recommend extensions from the **[Plugin Exp lorer](https://hiroyokoyama.github.io/moleditpy-plugins/explorer/)**.
"""

GENERATION_CONFIG = {
    "temperature": 0.1,
    "top_p": 0.95,
    "top_k": 40,
    "max_output_tokens": 16384,
    "response_mime_type": "text/plain",
}

# DEMO MODE: Set to True to test UI without API connection
DEMO_MODE = False
# History Limit: 10 Turns = 20 Messages (User + AI)
MAX_HISTORY = 20

# Button Colors
BTN_COLOR_ACCEPT_SINGLE = "#90ee90"  # Green - for single tool accept
BTN_COLOR_ACCEPT_ALL = "#87CEEB"     # Sky Blue - for accept all
BTN_COLOR_REJECT = "#ffcccb"          # Light red
BTN_COLOR_RETRY = "#ffcc00"           # Orange

SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "chat_with_molecule_neo.json")
LOG_FILE = os.path.join(os.path.dirname(__file__), "chat_with_molecule_neo.log")

# --- Google GenAI Import ---
# --- Markdown Import ---
try:
    import markdown
    HAS_MARKDOWN = True
except ImportError:
    HAS_MARKDOWN = False

try:
    import google.generativeai as genai
    HAS_GENAI = True
except ImportError:
    HAS_GENAI = False

def load_settings():
    if DEMO_MODE:
        return {}  # Don't load settings in demo mode
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, 'r') as f:
                return json.load(f)
        except:
            pass
    return {}

def save_settings(settings):
    if DEMO_MODE:
        return  # Don't overwrite settings in demo mode
    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(settings, f)
    except Exception as e:
        print(f"Error saving settings: {e}")

def append_log(sender, text):
    """Append message to local log file with timestamp"""
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        with open(LOG_FILE, 'a', encoding='utf-8') as f:
            f.write(f"[{timestamp}] {sender}: {text}\n")
    except Exception as e:
        print(f"Logging failed: {e}")

# --- LaTeX Cache ---
LATEX_CACHE = {}

def latex_to_html(latex_str):
    """
    Render LaTeX string to an HTML <img> tag using Matplotlib (OO Interface).
    Thread-safe implementation avoiding pyplot global state.
    """
    if not HAS_MATPLOTLIB:
        return f"<i>{latex_str}</i>" # Fallback if no matplotlib

    # Check Cache
    if latex_str in LATEX_CACHE:
        return LATEX_CACHE[latex_str]

    try:
        # Create a Figure object directly (OO Interface)
        # Scaled to fit text; start small (0.1x0.1) and resize automatically
        fig = Figure(figsize=(0.1, 0.1), dpi=120) 
        canvas = FigureCanvasAgg(fig)
        
        # Add text (math) to calculate size
        # Use simple dollar signs if not present?
        content = latex_str if latex_str.startswith('$') else f"${latex_str}$"
        text = fig.text(0, 0, content, fontsize=12)
        
        # Get dimensions (Tight Bounding Box)
        # We need a renderer to calculate text size
        renderer = canvas.get_renderer()
        bbox = text.get_window_extent(renderer)
        
        # Convert bbox width/height (pixels) to inches
        w_in = bbox.width / fig.dpi
        h_in = bbox.height / fig.dpi
        
        # Re-create figure with exact size (+padding)
        pad = 0.05
        fig.set_size_inches(w_in + pad, h_in + pad)
        
        # Clear and redraw centered
        fig.clear()
        text = fig.text(0.5, 0.5, content, fontsize=12, ha='center', va='center')
        
        # Render to buffer
        buf = io.BytesIO()
        canvas.print_png(buf)
        buf.seek(0)
        
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        
        # Cache Result
        html_tag = f'<img src="data:image/png;base64,{img_base64}" style="vertical-align: middle;">'
        LATEX_CACHE[latex_str] = html_tag
        return html_tag
        
    except Exception as e:
        print(f"LaTeX render error: {e}")
        return f"<code>{latex_str}</code>"

class ChatBrowser(QTextBrowser):
    """
    Custom QTextBrowser that disables internal navigation to prevent 'No document for...' errors.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        # CRITICAL: Disable internal navigation. 
        # If openLinks is True (default), clicking a link calls setSource/setDocument,
        # which clears the chat if it's a web URL or unknown scheme.
        self.setOpenLinks(False) 
        self.setOpenExternalLinks(False) # We handle everything in text_browser.anchorClicked
        
    def setSource(self, url):
        # OVERRIDE: Do absolutely nothing. 
        # Even with setOpenLinks(False), this acts as a failsafe.
        pass

class InitWorker(QThread):
    finished = pyqtSignal(list, str) # available_models, error_message

    def __init__(self, api_key):
        super().__init__()
        self.api_key = api_key

    def run(self):
        try:
            genai.configure(api_key=self.api_key)
            models = [m.name for m in genai.list_models() if 'generateContent' in m.supported_generation_methods]
            self.finished.emit(models, "")
        except Exception as e:
            self.finished.emit([], str(e))

class NameResolverWorker(QRunnable):
    """
    Worker thread to resolve InChIKey to Name asynchronously.
    """
    class Signals(QObject):
        finished = pyqtSignal(str, str) # inchikey, name

    def __init__(self, inchikey):
        super().__init__()
        self.inchikey = inchikey
        self.signals = self.Signals()

    def run(self):
        try:
            if PubChemResolver:
                name, error = PubChemResolver.resolve_inchikey_to_name(self.inchikey)
                if name:
                    self.signals.finished.emit(self.inchikey, name)
        except Exception:
            pass

class GenAIWorker(QThread):
    """Worker thread to handle API calls to avoid freezing UI"""
    response_received = pyqtSignal(object) # Changed to object to pass full response
    chunk_received = pyqtSignal(str)       # Signal for streaming updates
    error_occurred = pyqtSignal(str)

    def __init__(self, chat_session, user_message):
        super().__init__()
        self.chat_session = chat_session
        self.user_message = user_message

    def run(self):
        try:
            # Enable streaming
            response = self.chat_session.send_message(self.user_message, stream=True)
            
            # Iterate through chunks
            for chunk in response:
                if hasattr(chunk, 'text') and chunk.text:
                    self.chunk_received.emit(chunk.text)
            
            # Emit final complete response object for usage logging etc.
            self.response_received.emit(response)
        except Exception as e:
            self.error_occurred.emit(str(e))

class ChatMoleculeWindow(QDialog):
    def __init__(self, main_window):
        # Force parent=None to ensure it's a top-level window (avoids being hidden inside parent)
        # We process 'main_window' manually.
        super().__init__(None) 
        try:
            # QMessageBox.information(None, "Debug", "Window Init Start")
            self.main_window = main_window
            self.setWindowTitle("Chat with Molecule Neo (Gemini)")
            self.resize(500, 700)
            self.settings = load_settings()
            self.chat_session = None
            self.worker = None
    
            self.chat_history_log = [] # List of tuples (sender, text)
            self.last_smiles = None # Track previous state
            self.last_error = None # Track previous error state
            self.last_inchikey = None # To track molecule identity
            
            # Async Thread Pool
            self.thread_pool = QThreadPool()
            
            self.pending_context_msg = None # Pending context update to send with next user message
            self.pending_info_text = None # Pending visual notification to show in chat on next send
            self.first_check_done = False # Track if we've done the initial check
    
            # Streaming State
            self.stream_accumulated_text = ""
            self.stream_start_pos = 0
    
            self.init_ui()
            
        except Exception as e:
            print(f"ChatMoleculeWindow Init Error: {e}")
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self.main_window, "Plugin Init Error", f"Failed to initialize Chat Window:\n{e}")
            raise e
        
        # Defer initialization to allow window to show immediately
        QTimer.singleShot(100, self.initialize_session)

        # Polling Timer for Auto-Update
        # User Request: "remove 2-second polling. update only on input focus or send."
        # self.poll_timer = QTimer(self)
        # self.poll_timer.timeout.connect(self.check_molecule_change)
        # self.poll_timer.start(2000)

    def get_selected_atom_indices(self):
        """Helper: Get selected atom indices (1-based for MapNum)"""
        selected_ids = []
        try:
            if hasattr(self.main_window, 'data') and hasattr(self.main_window.data, 'atoms'):
                for atom_id, atom_data in self.main_window.data.atoms.items():
                    item = atom_data.get('item')
                    if item and item.isSelected():
                        selected_ids.append(str(atom_id + 1)) # Use 1-based (Matches MapNum)
        except Exception:
            pass
        return selected_ids

    def check_molecule_change(self):
        """Check if the main window's molecule OR selection has changed"""
        
        # 1. Get current SMILES
        current_smiles, error = self.get_current_molecule_smiles()
        
        # 2. Get current selection and sort (for comparison)
        current_selection = self.get_selected_atom_indices()
        # Sort as integer values to ensure consistency ("1", "10", "2" -> "1", "2", "10")
        current_selection.sort(key=int) 

        # Initialize last_selection if not exists
        if not hasattr(self, 'last_selection'):
            self.last_selection = []

        is_first_check = not self.first_check_done
        
        # 3. Detect changes
        smiles_changed = (current_smiles != self.last_smiles)
        error_changed = (error != self.last_error)
        selection_changed = (current_selection != self.last_selection)
        
        has_changed = smiles_changed or error_changed or selection_changed

        if is_first_check or has_changed:
             # --- Sync states ---
             self.last_smiles = current_smiles
             self.last_error = error
             self.last_selection = current_selection # Update selection state
             
             # Calculate InChIKey (only if molecule changed)
             if smiles_changed:
                 self.last_inchikey = None
                 if current_smiles and Chem:
                     try:
                         mol = Chem.MolFromSmiles(current_smiles)
                         if mol:
                             self.last_inchikey = Chem.MolToInchiKey(mol)
                     except:
                         pass

             self.first_check_done = True
             
             # --- UI Label Update ---
             if current_smiles:
                 # Standard context label (only SMILES/Name to keep it clean)
                 self.lbl_context.setText(f"Context: {self.get_molecule_name() or 'Unknown'} ({current_smiles[:20]}...)")
             else:
                 context_label_text = "Context: No valid molecule found."
                 if error:
                      context_label_text += f" ({error})"
                 self.lbl_context.setText(context_label_text)

             # --- Chat Logging & AI Context ---
             if not is_first_check:
                 if current_smiles:
                     # Get Name (Cache only)
                     mol_name = self.get_molecule_name(allow_fetch=False)
                     
                     # Re-build AI Context Message (Always updated with latest selection)
                     self.pending_context_msg = self._build_context_msg(current_smiles, mol_name, lazy=True)
                     
                     # --- Feedback Log Message ---
                     log_text = ""
                     
                     # Prepare selection info HTML
                     sel_info_html = ""
                     if current_selection:
                         sel_info_html = f" <span style='color:green'>[Selection: {', '.join(current_selection)}]</span>"
                     elif selection_changed and not current_selection:
                         sel_info_html = " <span style='color:gray'>[Selection cleared]</span>"

                     if smiles_changed:
                         # Structure changed: show SMILES + Selection
                         log_text = f"Context updated: {current_smiles}{sel_info_html}"
                         append_log("INFO", f"Context updated: {current_smiles}, Sel: {current_selection}")
                     
                     elif selection_changed:
                         # Only selection changed: simpler message
                         display_sel = ", ".join(current_selection) if current_selection else "None"
                         log_text = f"Selection updated: {display_sel}"
                         append_log("INFO", f"Selection updated: {current_selection}")

                     # Queue feedback for next message or display
                     if log_text:
                        self.pending_info_text = log_text
                     
                 else:
                     # Molecule unloaded
                     self.pending_context_msg = "System Update: The user has unloaded the molecule."
                     if error:
                         append_log("INFO", f"Molecule unloaded/scan failed: {error}")
                     else:
                         append_log("INFO", "Molecule unloaded.")
                     self.pending_info_text = "Molecule unloaded."

    def init_ui(self):
        # Set window icon (same as main app)
        try:
            # Try package import first (installed library)
            from moleditpy.modules import main_window_main_init
            script_dir = os.path.dirname(os.path.abspath(main_window_main_init.__file__))
            icon_path = os.path.join(script_dir, 'assets', 'icon.png')
        except ImportError:
            # Fallback for script execution
            try:
                from modules import main_window_main_init
                script_dir = os.path.dirname(os.path.abspath(main_window_main_init.__file__))
                icon_path = os.path.join(script_dir, 'assets', 'icon.png')
            except Exception:
                icon_path = None
        
        if icon_path and os.path.exists(icon_path):
            try:
                app_icon = QIcon(icon_path)
                self.setWindowIcon(app_icon)
            except Exception:
                pass
        
        layout = QVBoxLayout(self)

        # --- Settings (API Key & Model) ---
        settings_frame = QFrame()
        settings_frame.setFrameShape(QFrame.Shape.StyledPanel)
        settings_layout = QVBoxLayout(settings_frame) # Changed to Vertical
        settings_layout.setContentsMargins(5, 5, 5, 5)

        # API Key Row
        key_layout = QHBoxLayout()
        key_layout.addWidget(QLabel("API Key:"))
        self.txt_api_key = QLineEdit()
        self.txt_api_key.setEchoMode(QLineEdit.EchoMode.Password)
        self.txt_api_key.setText(self.settings.get("api_key", ""))
        self.txt_api_key.setPlaceholderText("Enter Google Gemini API Key")
        key_layout.addWidget(self.txt_api_key)

        # Link to get key
        lbl_link = QLabel('<a href="https://aistudio.google.com/app/api-keys">Get API Key</a>')
        lbl_link.setOpenExternalLinks(True)
        lbl_link.setToolTip("Click to open Google AI Studio and generate an API key")
        key_layout.addWidget(lbl_link)
        settings_layout.addLayout(key_layout)

        # Model Row
        model_layout = QHBoxLayout()
        model_layout.addWidget(QLabel("Model:"))
        self.combo_model = QComboBox()
        self.combo_model.setEditable(True) # Allow custom entry
        default_model = self.settings.get("model", "gemini-flash-latest")
        self.combo_model.addItem(default_model)
        self.combo_model.setCurrentText(default_model)
        model_layout.addWidget(self.combo_model, 1)

        btn_fetch_models = QPushButton("Fetch Models")
        btn_fetch_models.clicked.connect(self.fetch_models)
        model_layout.addWidget(btn_fetch_models)
        
        btn_save = QPushButton("Save & Reload")
        btn_save.clicked.connect(self.save_settings)
        model_layout.addWidget(btn_save)

        settings_layout.addLayout(model_layout)

        layout.addWidget(settings_frame)

        # --- Export Button ---
        btn_export = QPushButton("Export Chat History")
        btn_export.clicked.connect(self.export_history)
        layout.addWidget(btn_export)

        # --- Warning Label ---
        lbl_warning = QLabel(
            "<div style='color: orange; text-align: center;'>"
            "⚠️ Do not include confidential information.<br>"
            "LLM can make mistakes. Please verify important information."
            "</div>"
        )
        lbl_warning.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(lbl_warning)

        # --- Chat Display ---
        self.chat_display = ChatBrowser(self) # Use custom subclass
        self.chat_display.anchorClicked.connect(self.handle_link)
        self.chat_display.setReadOnly(True)
        # Increase base font size
        font = self.chat_display.font()
        font.setPointSize(12)
        self.chat_display.setFont(font)
        layout.addWidget(self.chat_display)

        # --- Context Status ---
        self.lbl_context = QLabel("Context: No molecule loaded.")
        self.lbl_context.setStyleSheet("color: gray; font-style: italic;")
        layout.addWidget(self.lbl_context)

        # --- Input Area ---
        input_layout = QHBoxLayout()
        self.txt_input = QTextEdit()
        self.txt_input.setFixedHeight(60)
        self.txt_input.setPlaceholderText("Type your message here... (Enter to send, Shift+Enter for new line)")
        self.txt_input.setEnabled(False) # Disable until initialized
        input_layout.addWidget(self.txt_input)

        self.btn_send = QPushButton("Send")
        self.btn_send.setFixedHeight(60)
        self.btn_send.clicked.connect(self.send_message)
        self.btn_send.setEnabled(False) # Disable until initialized
        input_layout.addWidget(self.btn_send)

        # --- Tool Confirmation Widget (Hidden by default) ---
        self.tool_confirm_frame = QFrame()
        self.tool_confirm_frame.setFrameShape(QFrame.Shape.StyledPanel)
        self.tool_confirm_frame.setStyleSheet("background-color: #f0f8ff; border: 1px solid #87ceeb; border-radius: 5px;")
        tool_layout = QVBoxLayout(self.tool_confirm_frame) # Changed to Vertical
        tool_layout.setContentsMargins(5, 5, 5, 5)
        
        self.lbl_tool_info = QLabel("Gemini wants to execute a command...")
        self.lbl_tool_info.setWordWrap(True) # Enable wrapping for long text
        tool_layout.addWidget(self.lbl_tool_info)
        
        # Buttons Row
        btn_layout = QHBoxLayout()
        
        # Add stretch to left to maintain alignment
        btn_layout.addStretch()
        
        self.btn_tool_accept = QPushButton("Accept")
        self.btn_tool_accept.setStyleSheet(f"background-color: {BTN_COLOR_ACCEPT_ALL}")
        self.btn_tool_accept.clicked.connect(self.on_tool_accept_click)
        btn_layout.addWidget(self.btn_tool_accept)
        
        self.btn_tool_step = QPushButton("Step-by-Step")
        self.btn_tool_step.setStyleSheet(f"background-color: {BTN_COLOR_ACCEPT_SINGLE}")
        self.btn_tool_step.clicked.connect(self.on_tool_step_click)
        self.btn_tool_step.setVisible(False)
        btn_layout.addWidget(self.btn_tool_step)
        
        self.btn_tool_reject = QPushButton("Reject")
        self.btn_tool_reject.setStyleSheet(f"background-color: {BTN_COLOR_REJECT}")
        self.btn_tool_reject.clicked.connect(self.on_tool_reject_click)
        btn_layout.addWidget(self.btn_tool_reject)
        
        tool_layout.addLayout(btn_layout)
        
        self.tool_confirm_frame.setVisible(False)
        layout.addWidget(self.tool_confirm_frame)

        layout.addLayout(input_layout)

        # Handle Enter to send
        self.txt_input.installEventFilter(self)

    def eventFilter(self, obj, event):
        if obj is self.txt_input:
            if event.type() == event.Type.FocusIn:
                # Check for updates as soon as user focuses on input
                self.check_molecule_change()
                
            elif event.type() == event.Type.KeyPress:
                if event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
                    if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                        # Allow Shift+Enter to insert new line (default behavior)
                        return False
                    else:
                        # Enter (without Shift) -> Send
                        self.send_message()
                        return True
        return super().eventFilter(obj, event)



    def closeEvent(self, event):
        """Handle window closure to clean up workers"""
        try:
            # Clean up ThreadPool
            if hasattr(self, 'thread_pool'):
                self.thread_pool.clear() 

            # Disconnect signals from active workers to prevent UI updates after close
            if self.worker and self.worker.isRunning():
                try: self.worker.chunk_received.disconnect()
                except: pass
                try: self.worker.response_received.disconnect()
                except: pass
                try: self.worker.finished.disconnect() 
                except: pass
                try: self.worker.error_occurred.disconnect() 
                except: pass
            
            if hasattr(self, 'init_worker') and self.init_worker and self.init_worker.isRunning():
                try: self.init_worker.finished.disconnect()
                except: pass

        except Exception as e:
            print(f"CloseEvent Error: {e}")
            
        event.accept()

    def fetch_models(self):
        """Fetch available Gemini models"""
        if DEMO_MODE:
             self.on_models_fetched(["demo-model"], None)
             return

        api_key = self.txt_api_key.text().strip()
        if not api_key:
             QMessageBox.warning(self, "Error", "Please enter an API Key first.")
             return
             
        self.combo_model.setEnabled(False) 
        self.append_message("System", "Fetching models...", "gray")
        
        worker = InitWorker(api_key)
        worker.finished.connect(self.on_models_fetched)
        worker.start()
        self.temp_worker = worker 

    def on_models_fetched(self, models, error):
        self.combo_model.setEnabled(True)
        if error:
             QMessageBox.critical(self, "Error", f"Failed to fetch models:\n{error}")
        else:
             self.combo_model.clear()
             self.combo_model.addItems(models)
             current = self.settings.get("model", "gemini-1.5-flash")
             if current in models:
                 self.combo_model.setCurrentText(current)
             # Non-intrusive feedback
             self.append_message("System", f"Found {len(models)} models.", "green")

    def save_settings(self):
        """Save settings from UI"""
        api_key = self.txt_api_key.text().strip()
        model = self.combo_model.currentText().strip()
        
        self.settings["api_key"] = api_key
        self.settings["model"] = model
        
        # Save to file
        save_settings(self.settings) 
        
        # Re-initialize
        self.initialize_session()
        # Non-blocking notification
        self.append_message("System", "Settings saved and session re-initialized.", "green")



    def handle_link(self, url):
        """Handle clickable links in chat window"""
        scheme = url.scheme()
        
        # Debugging: Print scheme and URL
        # print(f"Link Clicked: {url.toString()}, Scheme: {scheme}")

        if scheme == "smiles":
            # Extract SMILES string
            full_url = url.toString()
            if full_url.startswith("smiles:"):
                smiles = full_url[7:] # Remove 'smiles:'
            else:
                smiles = url.path() # Fallback
            
            if not smiles: return

            self.append_message("System", f"Loading molecule from SMILES: {smiles} ...", "blue")

            # Use MainWindow's importer if available
            if hasattr(self.main_window, 'main_window_string_importers'):
                try:
                    self.main_window.main_window_string_importers.load_from_smiles(smiles)
                    # self.append_message("System", "Molecule loaded successfully.", "green")
                    
                    # Force immediate context update
                    self.last_smiles = "FORCE_UPDATE"
                    self.check_molecule_change()
                    
                except Exception as e:
                    self.append_message("System", f"Failed to load molecule: {e}", "red")
            else:
                self.append_message("System", "Error: Molecule importer not found in main application.", "red")
        else:
            # Open other links (e.g. http) in external browser
            QDesktopServices.openUrl(url)



    def render_content(self, text):
        """Helper to render text (Markdown + LaTeX + SMILES) to HTML"""
        processed_text = text
        
        # --- LaTeX / Chemical Formula Formatting ---
        def replace_latex_block(match):
            return latex_to_html(match.group(1))
            
        def replace_latex_inline(match):
            return latex_to_html(match.group(1))

        # 1. Block Math: $$...$$
        processed_text = re.sub(r'\$\$([^$]+)\$\$', replace_latex_block, processed_text)
        
        # 2. Inline Math: $...$
        processed_text = re.sub(r'(?<!\$)\$([^$]+)\$(?!\$)', replace_latex_inline, processed_text)

        # 2.5. Simplify Tool JSON (New)
        def replace_tool_json(match):
            json_str = match.group(1)
            try:
                payload = json.loads(json_str)
                is_tool = False
                tools = []
                
                if isinstance(payload, list):
                    if all("tool" in item for item in payload):
                        is_tool = True
                        tools = payload
                elif isinstance(payload, dict) and "tool" in payload:
                    is_tool = True
                    tools = [payload]
                
                if is_tool:
                    # Improved HTML view with details
                    details_html = ""
                    for i, t in enumerate(tools):
                        t_name = t.get("tool", "Unknown")
                        t_params = t.get("params", {})
                        
                        # Format params nicely
                        param_str = ""
                        if t_params:
                            items = []
                            for k, v in t_params.items():
                                # Truncate long values for display (e.g. file content)
                                val_str = str(v)
                                if len(val_str) > 50: val_str = val_str[:47] + "..."
                                items.append(f"<b>{k}</b>: {val_str}")
                            param_str = f"<div style='margin-left: 15px; font-size: 0.9em; color: #555;'>{'<br>'.join(items)}</div>"
                        
                        prefix = f"{i+1}. " if len(tools) > 1 else ""
                        details_html += f"<div style='margin-bottom: 5px;'><b>{prefix}{t_name}</b>{param_str}</div>"

                    return (
                        f'<div style="background-color: #e7f5ff; padding: 10px; border-radius: 5px; '
                        f'border-left: 4px solid #0d6efd; margin: 5px 0;">'
                        f'<span style="font-weight: bold; color: #198754; font-size: 1.1em;">🛠️ Tool Request</span><br>'
                        f'{details_html}'
                        f'</div>'
                    )
                else:
                    return match.group(0) # Not a tool, keep raw
            except:
                return match.group(0) # Parse error, keep raw

        # Replace ```json ... ``` blocks if they contain tools
        processed_text = re.sub(r'```json\s*([\[{].*?[\]}])\s*```', replace_tool_json, processed_text, flags=re.DOTALL)

        # 3. Markdown
        if HAS_MARKDOWN:
            try:
                processed_text = markdown.markdown(processed_text, extensions=['fenced_code', 'tables'])
            except Exception as e:
                print(f"Markdown error: {e}")
                processed_text = processed_text.replace(chr(10), '<br>')
        else:
             processed_text = processed_text.replace(chr(10), '<br>')

        # 4. SMILES Links
        pattern = r"\[([^\]]+)\]\((smiles:[^\)]+)\)"
        processed_text = re.sub(pattern, r'<a href="\2">\1</a>', processed_text)
        
        return processed_text

    def append_message(self, sender, text, color="black"):
        self.chat_history_log.append({"sender": sender, "text": text})

        # Auto-log to file
        append_log(sender, text)

        # Render Content
        processed_text = self.render_content(text)

        # CSS Styling for Markdown elements
        style = """
        <style>
            body { font-size: 14px; line-height: 1.6; }
            code { background-color: #f0f0f0; padding: 2px; border-radius: 3px; font-family: monospace; }
            pre { background-color: #f0f0f0; padding: 10px; border-radius: 5px; border: 1px solid #ddd; }
            span.formula { font-family: 'Times New Roman', serif; }
            h1, h2, h3 { color: #2c3e50; }
            a { color: #3498db; text-decoration: none; font-weight: bold; }
        </style>
        """

        html = f"{style}<b style='color:{color}'>{sender}:</b><br>{processed_text}<br><br>"
        
        # Scroll Logic: Prevent System messages from scrolling user away
        v_bar = self.chat_display.verticalScrollBar()
        old_scroll = v_bar.value()
        
        self.chat_display.moveCursor(QTextCursor.MoveOperation.End)
        self.chat_display.insertHtml(html)
        
        if sender == "System":
             v_bar.setValue(old_scroll)
        else:
             self.chat_display.moveCursor(QTextCursor.MoveOperation.End)

    def export_history(self):
        if not self.chat_history_log:
             QMessageBox.information(self, "Export", "No history to export.")
             return

        file_path, filter = QFileDialog.getSaveFileName(
            self, "Export Chat History", "chat_history.md",
            "Markdown Files (*.md);;Text Files (*.txt);;JSON (*.json);;All Files (*)"
        )

        if not file_path:
            return

        try:
            import datetime
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            with open(file_path, 'w', encoding='utf-8') as f:
                if file_path.endswith('.json'):
                     json.dump(self.chat_history_log, f, indent=2)
                elif file_path.endswith('.md'):
                    f.write(f"# Chat with Molecule History\n")
                    f.write(f"Date: {current_time}\n")
                    f.write(f"Plugin Version: {PLUGIN_VERSION}\n\n")

                    for entry in self.chat_history_log:
                        sender = entry['sender']
                        text = entry['text']
                        role = "**You**" if sender == "You" else f"**{sender}**"
                        f.write(f"### {role}\n\n{text}\n\n---\n\n")
                else:
                    f.write(f"Chat History ({current_time})\n\n")
                    for entry in self.chat_history_log:
                        sender = entry['sender']
                        text = entry['text']
                        f.write(f"[{sender}]\n{text}\n\n{'='*40}\n\n")
                        
            QMessageBox.information(self, "Success", f"Exported to {os.path.basename(file_path)}")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to export: {e}")



    def initialize_session(self):
        """Initialize chat session (called shortly after init or on reload)"""
        
        # --- Demo Mode Handling ---
        if DEMO_MODE:
            self.settings["api_key"] = "DEMO_KEY_123" # Dummy key
            self.txt_api_key.setText("DEMO_KEY_123")
            self.combo_model.clear()
            self.combo_model.addItem("demo-model")
            self.chat_session = type('obj', (object,), {'history': []})() # Dummy object
            self.demo_step = 0 # Initialize step counter
            
            self.append_message("System", "<b>DEMO MODE ACTIVE</b><br>No API Key required. Responses are simulated.", "orange")
            # Enable UI
            self.txt_input.setEnabled(True)
            self.btn_send.setEnabled(True)
            self.txt_input.setFocus()
            
            # Trigger first check for molecule
            self.check_molecule_change()
            return

        if not HAS_GENAI:
            self.append_message("System", "Error: 'google-generativeai' library is not installed.", "red")
            self.txt_input.setEnabled(False)
            self.btn_send.setEnabled(False)
            return

        api_key = self.settings.get("api_key")
        if not api_key:
            self.append_message("System", "Please enter your Google API Key above.", "orange")
            QMessageBox.warning(self, "Welcome", "Please enter your Google Gemini API Key in the settings above to start chatting.")
            self.txt_api_key.setFocus()
            return

        self.append_message("System", "Checking available models...", "blue")
        
        # Start background worker
        self.init_worker = InitWorker(api_key)
        self.init_worker.finished.connect(self.on_init_finished)
        self.init_worker.start()

    def on_init_finished(self, available_models, error_msg):
        if error_msg:
             self.append_message("System", f"Failed to list models: {error_msg}", "red")
             return

        if not available_models:
            self.append_message("System", "No models found that support 'generateContent'.", "red")
            return

        # Update Dropdown
        self.combo_model.blockSignals(True)
        self.combo_model.clear()
        self.combo_model.addItems(available_models)
        self.combo_model.blockSignals(False)

        # Determine Model to Use
        saved_model = self.settings.get("model", "")
        target_model_name = None

        if saved_model in available_models:
            target_model_name = saved_model
        else:
            # Prioritize efficient 'flash' models as requested
            # We check for exact matches or partial matches for robustness
            preferred_defaults = ["models/gemini-flash-latest"]
            
            for pref in preferred_defaults:
                if pref in available_models:
                    target_model_name = pref
                    break
            
            # Fallback: any gemini, then just the first available
            if not target_model_name:
                geminis = [m for m in available_models if 'gemini' in m]
                if geminis:
                    target_model_name = geminis[0]
                else:
                    target_model_name = available_models[0]
            
            self.append_message("System", f"Default/Saved model '{saved_model}' not found. Auto-selecting: {target_model_name}", "gray")

        # Set Dropdown selection
        self.combo_model.setCurrentText(target_model_name)
        
        # Save the valid model we found so next time it's the default
        if target_model_name != saved_model:
            self.settings["model"] = target_model_name
            save_settings(self.settings)

        # --- Start Chat ---
        self.append_message("System", f"Initializing chat with: {target_model_name}", "gray")
        
        # Detailed Log
        append_log("INFO", f"Session Started. Plugin Version: {PLUGIN_VERSION}")
        append_log("INFO", f"Target Model: {target_model_name}")
        
        try:
            model = genai.GenerativeModel(
                target_model_name,
                system_instruction=SYSTEM_PROMPT,
                generation_config=GENERATION_CONFIG
            )
            self.chat_session = model.start_chat(history=[])
        except Exception as e:
            self.append_message("System", f"Error starting chat: {e}", "red")
            self.chat_session = None
            return

        # Context Injection
        smiles, error = self.get_current_molecule_smiles()
        self.last_smiles = smiles # Sync state
        
        # Calculate InChIKey immediately for name resolution
        self.last_inchikey = None
        if smiles and Chem:
             try:
                 mol = Chem.MolFromSmiles(smiles)
                 if mol:
                     self.last_inchikey = Chem.MolToInchiKey(mol)
             except:
                 pass

        context_msg = ""
        if smiles:
            mol_name = self.get_molecule_name()
            name_display = mol_name if mol_name else "Loaded Molecule"
            
            context_msg = self._build_context_msg(smiles, mol_name) # Call the class methodhelper
            self.pending_context_msg = context_msg # Store to send with first user message
            self.pending_context_msg = context_msg # Store to send with first user message
            self.pending_info_text = f"Context loaded: {name_display} ({smiles})"
        else:
            context_label_text = "Context: No valid molecule found."
            if error:
                 context_label_text += f" ({error})"
            self.lbl_context.setText(context_label_text)
            
            # Even with no context, teach formatting
            context_msg = (
                "No molecule is currently loaded. "
            )
            if error:
                 append_log("INFO", f"No molecule context found: {error}")
            else:
                 append_log("INFO", "No molecule context found.")
            self.pending_context_msg = context_msg

        # Enable inputs immediately (nothing is being sent yet)
        self.txt_input.setEnabled(True)
        self.btn_send.setEnabled(True)
        self.txt_input.setFocus()

    def get_current_molecule_smiles(self):
        """Returns (smiles, error_message) with Atom Map Numbers matching UI IDs (1-based)."""
        if not Chem:
            return None, "RDKit library not found/loaded."
            
        if not self.main_window:
            return None, "Internal Error: No reference to main application."
            
        mol = None
        
        # 1. Try reconstructing from 2D Data (Preserves internal IDs)
        try:
            if hasattr(self.main_window, 'data') and hasattr(self.main_window.data, 'to_rdkit_mol'):
                if hasattr(self.main_window.data, 'atoms') and self.main_window.data.atoms:
                    mol = self.main_window.data.to_rdkit_mol()
                    if mol:
                        # Map internal IDs to MapNums
                        for atom in mol.GetAtoms():
                            if atom.HasProp("_original_atom_id"):
                                aid = atom.GetIntProp("_original_atom_id")
                                atom.SetAtomMapNum(aid + 1) # Use 1-based for MapNum
        except Exception:
            mol = None
    
        # 2. Fallback: Try cached 3D mol
        if mol is None:
            try:
                if hasattr(self.main_window, 'current_mol'):
                    mol = self.main_window.current_mol
            except Exception:
                pass

        if mol is None:
             if hasattr(self.main_window, 'data') and hasattr(self.main_window.data, 'atoms') and not self.main_window.data.atoms:
                 return None, "No molecule loaded (Scene is empty)."
             return None, "No valid molecule found."
            
        try:
            # Remove hydrogens for Gemini prompt, but preserve MapNums on heavy atoms
            mol_clean = Chem.RemoveHs(mol)
            
            # Ensure every heavy atom has a MapNum
            for atom in mol_clean.GetAtoms():
                if atom.GetAtomMapNum() == 0:
                    atom.SetAtomMapNum(atom.GetIdx() + 1)
                    
            return Chem.MolToSmiles(mol_clean), None
        except Exception as e:
            return None, f"SMILES conversion failed: {str(e)}"

    def load_smiles_undo_safe(self, smiles_string):
        """
        Load SMILES into MainWindow but preserve Undo Stack (Push state instead of Reset).
        Replicates logic from MainWindowStringImporters.load_from_smiles but without reset_undo_stack().
        """
        mw = self.main_window
        try:
            cleaned_smiles = smiles_string.strip()
            mol = Chem.MolFromSmiles(cleaned_smiles)
            
            # Fallback: Relaxed loading if strict failed
            if mol is None:
                try:
                    mol = Chem.MolFromSmiles(cleaned_smiles, sanitize=False)
                    if mol:
                        mol.UpdatePropertyCache(strict=False)
                        self.append_message("System", "Warning: Molecule loaded with relaxed validation (check structure).", "orange")
                except:
                    mol = None

            if mol is None:
                self.append_message("System", f"Error: Generated SMILES is invalid.\nInput: {cleaned_smiles[:50]}...", "red")
                return

            AllChem.Compute2DCoords(mol)
            Chem.Kekulize(mol)
            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)

            # --- UNDO MAGIC ---
            # Push current state BEFORE clearing
            mw.push_undo_state()
            
            # Clear editor WITHOUT pushing another undo state
            mw.clear_2d_editor(push_to_undo=False)
            
            mw.current_mol = None
            mw.plotter.clear()
            mw.analysis_action.setEnabled(False)

            # --- RECONSTRUCT SCENE ---
            SCALE_FACTOR = 50.0
            view_center = mw.view_2d.mapToScene(mw.view_2d.viewport().rect().center())
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            mol_center_x = sum(p.x for p in positions) / len(positions) if positions else 0.0
            mol_center_y = sum(p.y for p in positions) / len(positions) if positions else 0.0

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()
                
                relative_x = pos.x - mol_center_x
                relative_y = pos.y - mol_center_y
                
                scene_x = (relative_x * SCALE_FACTOR) + view_center.x()
                scene_y = (-relative_y * SCALE_FACTOR) + view_center.y()
                
                # Access mw.scene directly
                atom_id = mw.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge)
                rdkit_idx_to_my_id[i] = atom_id

            for bond in mol.GetBonds():
                b_idx, e_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble()
                b_dir = bond.GetBondDir()
                stereo = 0
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4

                if b_idx in rdkit_idx_to_my_id and e_idx in rdkit_idx_to_my_id:
                    a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                    a1_item = mw.data.atoms[a1_id]['item']
                    a2_item = mw.data.atoms[a2_id]['item']
                    mw.scene.create_bond(a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo)

            # --- FINALIZE ---
            mw.has_unsaved_changes = True
            
            # Update UI
            if mw.data.atoms:
                mw.view_2d.fitInView(mw.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
            mw.update_realtime_info()
            mw.update_undo_redo_actions()
            mw.update_window_title()

            # Push Final State so Undo works (Current state is now on top of stack)
            mw.push_undo_state()
            
            # --- Check for Chemistry Problems ---
            try:
                if mw.check_chemistry_problems_fallback():
                     pass
            except:
                pass
            
            mw.update_undo_redo_actions()
            mw.update_window_title()
            QTimer.singleShot(0, mw.fit_to_view)
            
            # Force Scene Update
            if hasattr(mw, 'scene'): mw.scene.update()
            if hasattr(mw, 'view_2d'): mw.view_2d.viewport().update()
            
            # Force check to update context immediately
            self.check_molecule_change()
            
        except Exception as e:
            self.append_message("System", f"Error loading molecule: {e}", "red")

    def get_molecule_name(self, allow_fetch=False):
        """Try to resolve name from current InChIKey using PubChemResolver (Async if allow_fetch=True)"""
        if not self.last_inchikey or not PubChemResolver:
            return None
        
        # Check Cache first
        if not hasattr(self, '_name_cache'):
            self._name_cache = {}
            
        if self.last_inchikey in self._name_cache:
            return self._name_cache[self.last_inchikey]
            
        # Optimization: Start Async Worker if not already fetching
        if allow_fetch:
            if not hasattr(self, '_fetching_inchikeys'):
                self._fetching_inchikeys = set()
                
            if self.last_inchikey not in self._fetching_inchikeys:
                self._fetching_inchikeys.add(self.last_inchikey)
                worker = NameResolverWorker(self.last_inchikey)
                worker.signals.finished.connect(self.on_name_resolved)
                self.thread_pool.start(worker)
            
        return None # Return None initially (loading)

    def on_name_resolved(self, inchikey, name):
        """Callback when name is resolved successfully"""
        if hasattr(self, '_fetching_inchikeys'):
            self._fetching_inchikeys.discard(inchikey)
            
        if not hasattr(self, '_name_cache'):
            self._name_cache = {}
            
        self._name_cache[inchikey] = name
        
        # Update UI if this is still the current molecule
        if self.last_inchikey == inchikey:
            self.lbl_context.setText(f"Context: {name} (Updated)")
            
            # Since check_molecule_change might have finished with "None/Unknown",
            # we should update the valid pending context message if it hasn't been sent yet.
            if self.pending_context_msg and "Name: " not in self.pending_context_msg:
                 # Inject name into pending message
                 # Simple replacement for "Name: None" or just prepend?
                 # It's cleaner to just regenerate the context string.
                 # Call check_molecule_change again? 
                 # Safer: just set it.
                 pass
            
            # Re-trigger check to update pending messages cleanly
            self.check_molecule_change()

    def _get_descriptors_str(self, smiles):
        """Calculate basic descriptors for unknown molecules to help LLM"""
        try:
            from rdkit.Chem import Descriptors, rdMolDescriptors
            
            # User Request: "Use mw.draw_molecule_3D... avoid freeze" 
            # Update: User said "Only just before sending". So freezing/blocking is acceptable here.
            # Strategy: Generate 3D locally for ACCURATE descriptors.
            
            mol = None
            if hasattr(self.main_window, 'current_mol') and self.main_window.current_mol:
                mol = self.main_window.current_mol
            
            # If no current 3D mol (because we disabled auto-trigger), generate temp one.
            if not mol or mol.GetNumConformers() == 0:
                 mol = Chem.MolFromSmiles(smiles)
                 if mol:
                     try:
                         mol = Chem.AddHs(mol)
                         AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                     except:
                         pass

            if not mol: return ""

            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            num_rings = rdMolDescriptors.CalcNumRings(mol)
            
            return f"Properties: Formula={formula}, MW={mw:.2f}, LogP={logp:.2f}, TPSA={tpsa:.2f}, NumRings={num_rings}. "
        except:
            return ""


    def propose_tool_action(self, payload):
        """Show the confirmation UI for the proposed action"""
        # Reset State
        self.tool_state = "CONFIRM"
        self.btn_tool_accept.setText("Accept")
        self.btn_tool_reject.setText("Reject")
        self.btn_tool_accept.setStyleSheet(f"background-color: {BTN_COLOR_ACCEPT_ALL}")

        self.pending_tool_payload = payload
        
        # Check if multiple tools
        is_multiple = "tools" in payload and len(payload["tools"]) > 1

        if "tools" in payload:
            tools_list = payload["tools"]
            tool_names = [t.get("tool", "?") for t in tools_list]
            desc = f"<b>{len(tools_list)} Operations</b>: " + " → ".join(tool_names)
        else:
            tool_name = payload.get("tool")
            params = payload.get("params", {})
            
            # Update UI label
            desc = f"Action: <b>{tool_name}</b>"
            if tool_name == "apply_transformation":
                desc += f" (Reaction SMARTS: {params.get('reaction_smarts')})"
            elif tool_name == "highlight_substructure":
                desc += f" (SMARTS: {params.get('smarts')})"
            elif tool_name == "calculate_descriptors":
                desc += f" (Properties: {', '.join(params.get('properties', []))})"
            elif tool_name == "orca_input_generator":
                desc += f" (ORCA Input: {params.get('filename')})"
            elif tool_name == "gaussian_input_generator":
                desc += f" (Gaussian Input: {params.get('filename')})"
            elif tool_name == "save_file":
                desc += f" (Save File: {params.get('filename')})"
            elif tool_name == "load_molecule":
                desc += f" (Name: {params.get('name', params.get('smiles', ''))})"
            elif tool_name == "set_electronic_state":
                desc += f" (Charge: {params.get('charge', '-')}, Mult: {params.get('multiplicity', '-')})"
            elif tool_name == "convert_to_3d":
                desc += " (Show 3D structure)"
            elif tool_name == "clear_canvas":
                desc += " (Clear all)"
            
        self.lbl_tool_info.setText(f"Gemini suggests: {desc}")
        self.lbl_tool_info.setTextFormat(Qt.TextFormat.RichText)
        
        # Show Step button only if multiple tools are queued
        if is_multiple:
            next_tool_name = payload["tools"][0].get("tool", "Next Step")
            self.btn_tool_step.setText("Accept")
            self.btn_tool_step.setStyleSheet(f"background-color: {BTN_COLOR_ACCEPT_SINGLE}")
            self.btn_tool_step.setVisible(True)
            self.btn_tool_accept.setText("Accept All")
            self.btn_tool_accept.setStyleSheet(f"background-color: {BTN_COLOR_ACCEPT_ALL}")
            self.btn_tool_accept.setVisible(True)  # Show Accept All
        else:
            self.btn_tool_step.setVisible(False)
            self.btn_tool_accept.setText("Accept")
            self.btn_tool_accept.setStyleSheet(f"background-color: {BTN_COLOR_ACCEPT_SINGLE}")
            self.btn_tool_accept.setVisible(True)  # Always show for single tool

        self.tool_confirm_frame.setVisible(True)
        self.txt_input.setEnabled(False)  # Freeze input while tool confirmation is shown
        # Scroll to bottom to ensure user sees it? It's above input, so it should be visible.
        
    def on_tool_step_click(self):
        """Execute only the next tool in the chain"""
        if not self.pending_tool_payload: return
        
        payload = self.pending_tool_payload
        
        # Must be a list if we are stepping
        if "tools" not in payload: 
             # Fallback to normal accept if somehow single
             self.accept_tool_action()
             return

        tools_list = payload["tools"]
        if not tools_list: return

        # 1. Pop first tool
        current_tool = tools_list.pop(0)
        
        # 2. Execute it
        tool_name = current_tool.get("tool")
        params = current_tool.get("params", {})
        
        self.append_message("System", f"Step Executing: {tool_name}...", "blue")
        QApplication.processEvents()
        
        error = self._dispatch_tool(tool_name, params)
        
        if error:
            # Step Fail: Stop here
            self.last_tool_error = error
            self.tool_state = "RETRY"
            self.lbl_tool_info.setText(f"Tool '{tool_name}' failed during step: {error}")
            self.btn_tool_accept.setText("Retry Step (Ask Gemini)")
            self.btn_tool_reject.setText("Stop Chain")
            self.btn_tool_step.setVisible(False) 
            self.btn_tool_accept.setStyleSheet("background-color: #ffcc00;") # Orange
            return

        # 3. Check Remaining
        if not tools_list:
            # All done
            self.tool_confirm_frame.setVisible(False)
            self.pending_tool_payload = None
            self.append_message("System", "All steps completed.", "green")
        else:
            # Re-propose remaining items
            self.pending_tool_payload = {"tools": tools_list}
            
            # Refresh UI (This will update the count "2 Operations: ..." and show Step button again)
            # We call propose_tool_action again with remaining list
            self.propose_tool_action(self.pending_tool_payload)
            self.append_message("System", "Step completed. Waiting for next...", "green")

    def accept_tool_action(self):
        """Execute the pending tool action(s)"""
        if not self.pending_tool_payload: return
        
        self.tool_confirm_frame.setVisible(False)
        self.txt_input.setEnabled(True)  # Unfreeze input
        self.btn_send.setEnabled(True)   # Re-enable send button
        payload = self.pending_tool_payload
        self.pending_tool_payload = None
        
        error = None
        failed_tool = ""
        
        # Check if multiple tools (array)
        if "tools" in payload:
            tools_list = payload["tools"]
            self.append_message("System", f"Executing {len(tools_list)} operations in order...", "blue")
            
            # Execute each tool in order
            for i, tool_item in enumerate(tools_list):
                tool_name = tool_item.get("tool")
                params = tool_item.get("params", {})
                
                self.append_message("System", f"[{i+1}/{len(tools_list)}] {tool_name}...", "gray")
                QApplication.processEvents()
                
                error = self._dispatch_tool(tool_name, params)
                if error:
                    failed_tool = tool_name
                    break
        else:
            # Single tool
            tool_name = payload.get("tool")
            params = payload.get("params", {})
            
            self.append_message("System", f"Executing tool: {tool_name}...", "blue")
            error = self._dispatch_tool(tool_name, params)
            if error: failed_tool = tool_name
            
        if error:
            # Switch to Retry Mode
            self.last_tool_error = error
            self.tool_state = "RETRY"
            
            self.lbl_tool_info.setText(f"Tool '{failed_tool}' failed: {error}")
            self.btn_tool_accept.setText(" Retry (Ask Gemini) ")
            self.btn_tool_reject.setText(" Ignore ")
            self.btn_tool_step.setVisible(False) # Hide Step button on error
            self.btn_tool_accept.setStyleSheet("background-color: #ffcc00;") # Orange
            self.tool_confirm_frame.setVisible(True)

    def _dispatch_tool(self, tool_name, params):
        """Dispatch a single tool execution"""
        # Dispatch
        if tool_name == "apply_transformation":
            return self.execute_apply_transformation(params)
        elif tool_name == "highlight_substructure":
            return self.execute_highlight_substructure(params)
        elif tool_name == "calculate_descriptors":
            return self.execute_calculate_descriptors(params)
        elif tool_name == "orca_input_generator":
            return self.execute_orca_input_generator(params)
        elif tool_name == "gaussian_input_generator":
            return self.execute_gaussian_input_generator(params)
        elif tool_name == "save_file":
            return self.execute_save_file(params)
        elif tool_name == "load_molecule":
            return self.execute_load_molecule(params)
        elif tool_name == "convert_to_3d":
            return self.execute_convert_to_3d(params)
        elif tool_name == "set_electronic_state":
            return self.execute_set_electronic_state(params)
        elif tool_name == "clear_canvas":
            return self.execute_clear_canvas(params)
        else:
            self.append_message("System", f"Unknown tool: {tool_name}", "red")
            return f"Unknown tool: {tool_name}"


    def reject_tool_action(self):
        """Cancel the pending tool action"""
        self.tool_confirm_frame.setVisible(False)
        self.pending_tool_payload = None
        
        # UI Feedback
        self.append_message("System", "Tool execution rejected by user.", "orange")
        
        # Batching: Queue this info for the next message to save an API turn
        if self.pending_context_msg:
            self.pending_context_msg += "\n[System: User rejected the previous tool suggestion.]"
        else:
            self.pending_context_msg = "[System: User rejected the previous tool suggestion.]"

    def on_tool_accept_click(self):
        state = getattr(self, 'tool_state', 'CONFIRM')
        if state == 'RETRY':
            self.retry_tool_action()
        else:
            self.accept_tool_action()

    def on_tool_reject_click(self):
        state = getattr(self, 'tool_state', 'CONFIRM')
        self.tool_confirm_frame.setVisible(False)
        
        # Re-enable UI and restore focus
        self.txt_input.setEnabled(True)
        self.btn_send.setEnabled(True)
        self.txt_input.setFocus()
        
        if state == 'CONFIRM':
            self.reject_tool_action()

    def retry_tool_action(self):
        """Ask Gemini to retry based on last error"""
        self.tool_confirm_frame.setVisible(False)
        self.tool_state = "CONFIRM"
        
        error_msg = getattr(self, 'last_tool_error', 'Unknown Error')
        
        # Send message
        retry_msg = f"The previous tool execution failed with error: {error_msg}. Please try again with corrected parameters."
        
        self.txt_input.setEnabled(True)
        self.btn_send.setEnabled(True)
        self.txt_input.setText(retry_msg)
        self.send_message()

    # --- Tool Implementations ---




    def execute_apply_transformation(self, params):
        """Execute chemical transformation using Reaction SMARTS"""
        try:
             reaction_smarts = params.get("reaction_smarts")
             if not reaction_smarts:
                 self.append_message("System", "Error: No reaction SMARTS provided.", "red")
                 return

             # Use method to safely get SMILES/Mol
             current_smiles, error = self.get_current_molecule_smiles()
             if not current_smiles:
                 self.append_message("System", f"Error: {error}", "red")
                 return

             if not Chem or not AllChem:
                  self.append_message("System", "RDKit is required for this feature.", "red")
                  return

             mol = Chem.MolFromSmiles(current_smiles)
             if not mol:
                 self.append_message("System", "Error: Could not parse current molecule.", "red")
                 return
             
             # Reaction Logic
             rxn = AllChem.ReactionFromSmarts(reaction_smarts)
             
             # User Request: "omitted hydrogens should be shown before H replacement"
             # If the user wants to substitute H (e.g. H->F), RDKit reactions on implicit Hs usually fail 
             # or don't work as expected unless mapped carefully.
             # It is safer to make Hs explicit before running the reaction.
             mol_with_h = Chem.AddHs(mol)
             
             # Run reaction
             products = rxn.RunReactants((mol_with_h,))
             
             if not products:
                  # Retry without explicit Hs just in case the SMARTS was designed for implicit
                  # But usually SMARTS with [H] require explicit Hs.
                  self.append_message("System", "Reaction produced no products with explicit Hs. Checking implicit...", "orange")
                  products = rxn.RunReactants((mol,))
                  if not products:
                      self.append_message("System", "Reaction produced no products. Check compatibility.", "orange")
                      return
             
             # Take first product
             new_mol = products[0][0]
             
             # ROBUST FIX: Implicit Hydrogens & Valence
             # 1. Update Property Cache
             try: new_mol.UpdatePropertyCache(strict=False)
             except: pass
             
             # 2. Sanitize to fix valence
             try: Chem.SanitizeMol(new_mol)
             except: pass
             
             # 3. Explicitly Remove hydrogens (often spectators in reaction)
             # This forces RDKit to recalculate implicit H count based on new connectivity
             try: new_mol = Chem.RemoveHs(new_mol, implicitOnly=False, updateExplicitCount=True, sanitize=True)
             except: pass
             
             # 4. SMILES Round-Trip (The Nuclear Option)
             # This guarantees a clean graph with correct valences
             temp_smiles = Chem.MolToSmiles(new_mol)
             clean_mol = Chem.MolFromSmiles(temp_smiles)
             
             if not clean_mol:
                 # Fallback if round-trip fails (unlikely)
                 clean_mol = new_mol
             
             # 5. Optimize 2D Coordinates (Fresh Layout)
             AllChem.Compute2DCoords(clean_mol)
             
             # Load back
             final_smiles = Chem.MolToSmiles(clean_mol)
             
             # Use local undo-safe loader
             self.load_smiles_undo_safe(final_smiles)
             
             # User Request: "Optimize 2D" after conversion
             if hasattr(self.main_window, 'clean_up_2d_structure'):
                 self.main_window.clean_up_2d_structure()
                 
             self.append_message("System", f"Transformation Applied.\nRule: `{reaction_smarts}`", "green")
                 
        except Exception as e:
            self.append_message("System", f"Transformation Failed: {e}\n(Tip: The SMARTS pattern might be invalid. Ask Gemini for a more precise one.)", "red")
            return str(e)

    def execute_highlight_substructure(self, params):
        """Highlight atoms matching SMARTS or explicit indices"""
        try:
            smarts = params.get("smarts")
            atom_indices_param = params.get("atom_indices")
            
            if not smarts and not atom_indices_param:
                return
                
            current_smiles, error = self.get_current_molecule_smiles()
            if not current_smiles:
                return

            mol = Chem.MolFromSmiles(current_smiles)
            if not mol: return
            
            # Use SET (Unique atoms)
            atom_indices_set = set()

            # A) Explicit Indices provided (Map Numbers)
            if atom_indices_param:
                 # Map explicit IDs/MapNums to RDKit internal indices
                 needed_maps = [int(x) for x in atom_indices_param]
                 for atom in mol.GetAtoms():
                     if atom.GetAtomMapNum() in needed_maps:
                         atom_indices_set.add(atom.GetIdx())

            # B) SMARTS provided
            if smarts:
                query = Chem.MolFromSmarts(smarts)
                
                if not query:
                     self.append_message("System", f"Invalid SMARTS: {smarts}", "red")
                     # If indices were okay, continue? Or stop?
                     # Let's stop to be safe, or just warn.
                     if not atom_indices_set: return
                else:
                    matches = mol.GetSubstructMatches(query)
                    if matches:
                        for match in matches:
                            atom_indices_set.update(match)
                    elif not atom_indices_set:
                        self.append_message("System", f"No matches found for: {smarts}", "orange")
                        return

            # Proceed to highlighting
            # ... existing highlighting logic ...
            
            # Highlighting using Map Numbers (Robust)
            # The SMILES from get_current_molecule_smiles() includes Map Numbers [C:1].
            # These MapNums correspond to the keys in mw.data.atoms via the exporter.
            
            highlight_ids = []
            
            for rdf_idx in atom_indices_set:
                atom = mol.GetAtomWithIdx(rdf_idx)
                # MapNum is 1-based usually, matching our dict keys if generated by our exporter
                map_num = atom.GetAtomMapNum() 
                if map_num > 0:
                    highlight_ids.append(map_num - 1) # Back to 0-based internal ID
                else:
                    pass

            count = 0
            if hasattr(self.main_window, 'data') and hasattr(self.main_window.data, 'atoms'):
                if hasattr(self.main_window, 'scene'):
                     self.main_window.scene.clearSelection()

                for aid in highlight_ids:
                    # Try to find aid in atoms.keys() (int or str mismatch check)
                    # dict keys might be int, aid is int.
                    if aid in self.main_window.data.atoms:
                        item = self.main_window.data.atoms[aid].get('item')
                        if item:
                            item.setSelected(True)
                            count += 1
            
            highlight_msg = f"Selected {count} atoms."
            self.append_message("System", highlight_msg, "green")
            
        except Exception as e:
            self.append_message("System", f"Highlight Error: {e}", "red")
            return str(e)


    def execute_set_electronic_state(self, params):
        """Set Formal Charge / Multiplicity on specific atoms or global state"""
        try:
            current_smiles, error = self.get_current_molecule_smiles()
            if not current_smiles:
                return "No molecule loaded."

            mol = Chem.MolFromSmiles(current_smiles)
            if not mol: return "Could not parse current molecule."

            # Params
            # Params
            target_smarts = params.get("target_smarts")
            atom_indices_param = params.get("atom_indices")
            new_charge = params.get("charge")
            new_mult = params.get("multiplicity")
            mode = params.get("mode", "atom") # 'atom' or 'global'
            
            # 1. Atom-Specific Changes (Robust RDKit modification)
            changes_made = False
            
            # Target Selection from SMARTS or User Selection or Indices
            target_indices = set()
            
            # A) Explicit Indices provided (Map Numbers)
            if atom_indices_param:
                 # Handle single int or list
                 if isinstance(atom_indices_param, int):
                     needed_maps = [atom_indices_param]
                 else:
                     needed_maps = [int(x) for x in atom_indices_param]
                 for atom in mol.GetAtoms():
                     if atom.GetAtomMapNum() in needed_maps:
                         target_indices.add(atom.GetIdx())

            # B) SMARTS provided
            elif target_smarts:
                query = Chem.MolFromSmarts(target_smarts)
                if query:
                    matches = mol.GetSubstructMatches(query)
                    for m in matches: target_indices.update(m)
            
            # B) User Selection (if no specific SMARTS)
            # If AI didn't specify SMARTS, check if there's a UI Selection to apply to
            elif hasattr(self.main_window, 'data') and hasattr(self.main_window.data, 'atoms'):
                 # We need to map UI Selection -> RDKit Indices
                 # This requires map numbers again
                 selected_map_nums = []
                 for aid, adata in self.main_window.data.atoms.items():
                      item = adata.get('item')
                      if item and item.isSelected():
                          selected_map_nums.append(aid + 1) # MapNum is ID + 1
                 
                 if selected_map_nums:
                      for atom in mol.GetAtoms():
                          if atom.GetAtomMapNum() in selected_map_nums:
                              target_indices.add(atom.GetIdx())
            
            # Apply Changes
            if target_indices and new_charge is not None:
                for idx in target_indices:
                    atom = mol.GetAtomWithIdx(idx)
                    atom.SetFormalCharge(int(new_charge))
                    # Reset radical info if charge changes (assumption)
                    # or handle multiplicity if provided
                    if new_mult is not None:
                         # Multiplicity is harder on generic atoms.
                         # RDKit uses NumRadicalElectrons.
                         # Mult = N_rad + 1.  N_rad = Mult - 1
                         n_rad = int(new_mult) - 1
                         if n_rad >= 0:
                             atom.SetNumRadicalElectrons(n_rad)
                    changes_made = True
            
            # Global Metadata fallback (if no atoms targeted but global requested)
            if not changes_made and (new_charge is not None or new_mult is not None):
                 # Just explicitly update the global state variables without changing atoms
                 # This is for "Metadata" updates
                 msg = []
                 if new_charge is not None and hasattr(self.main_window, 'current_charge'):
                     self.main_window.current_charge = int(new_charge)
                     msg.append(f"Global Charge={new_charge}")
                 if new_mult is not None and hasattr(self.main_window, 'current_mult'):
                     self.main_window.current_mult = int(new_mult)
                     msg.append(f"Global Mult={new_mult}")
                 
                 if msg:
                      if hasattr(self.main_window, 'update_realtime_info'):
                           self.main_window.update_realtime_info()
                      self.append_message("System", f"Updated Global Metadata: {', '.join(msg)} (Atoms unchanged)", "orange")
                      return
            
            if changes_made:
                 # Sanitize & Reload
                 try:
                     mol.UpdatePropertyCache(strict=False)
                     Chem.SanitizeMol(mol)
                     # Determine new global charge from atoms
                     total_charge = Chem.GetFormalCharge(mol)
                     if hasattr(self.main_window, 'current_charge'):
                         self.main_window.current_charge = total_charge
                         
                     new_smiles = Chem.MolToSmiles(mol)
                     self.load_smiles_undo_safe(new_smiles)
                     self.append_message("System", f"Applied Electronic State Changes. (Net Charge: {total_charge})", "green")
                 except Exception as e:
                     self.append_message("System", f"Error applying state: {e}", "red")
            else:
                 # DEBUG INFO
                 all_maps = [a.GetAtomMapNum() for a in mol.GetAtoms()]
                 debug_msg = f"Available MapNums: {all_maps}. "
                 if atom_indices_param: debug_msg += f"Requested: {atom_indices_param}. "
                 if target_smarts: debug_msg += f"SMARTS: {target_smarts}. "
                 
                 self.append_message("System", f"No atoms targeted. {debug_msg} Select atoms or use SMARTS.", "orange")

        except Exception as e:
            self.append_message("System", f"Error setting state: {e}", "red")
            return str(e)


    def execute_calculate_descriptors(self, params):
        """Calculate properties"""
        try:
             req_props = params.get("properties", ["MW", "LogP"])
             current_smiles, error = self.get_current_molecule_smiles()
             if not current_smiles:
                 return
             
             mol = Chem.MolFromSmiles(current_smiles)
             from rdkit.Chem import Descriptors, Lipinski
             
             results = []
             # Molecular Weight
             if "MW" in req_props:
                 results.append(f"MW: {Descriptors.MolWt(mol):.2f}")
             
             # LogP (Partition Coefficient)
             if "LogP" in req_props:
                 results.append(f"LogP: {Descriptors.MolLogP(mol):.2f}")
             
             # Topological Polar Surface Area
             if "TPSA" in req_props:
                 results.append(f"TPSA: {Descriptors.TPSA(mol):.2f}")
             
             # Hydrogen Bond Donors
             if "HBondDonor" in req_props or "HBD" in req_props:
                 results.append(f"H-Bond Donors: {Lipinski.NumHDonors(mol)}")
             
             # Hydrogen Bond Acceptors
             if "HBondAcceptor" in req_props or "HBA" in req_props:
                 results.append(f"H-Bond Acceptors: {Lipinski.NumHAcceptors(mol)}")
             
             # Rotatable Bonds
             if "RotatableBonds" in req_props or "RB" in req_props:
                 results.append(f"Rotatable Bonds: {Lipinski.NumRotatableBonds(mol)}")
             
             # Aromatic Rings
             if "AromaticRings" in req_props:
                 results.append(f"Aromatic Rings: {Lipinski.NumAromaticRings(mol)}")
             
             # Number of Rings
             if "Rings" in req_props:
                 results.append(f"Rings: {Lipinski.RingCount(mol)}")
             
             result_text = "Properties:\n" + "\n".join(results)
             
             # Show in chat
             self.append_message("System", result_text, "blue")
             
             # Show in popup dialog
             msg_box = QMessageBox(self)
             msg_box.setWindowTitle("Molecular Properties")
             msg_box.setText(result_text)
             msg_box.setIcon(QMessageBox.Icon.Information)
             msg_box.exec()
             
        except Exception as e:
             self.append_message("System", f"Calculation Error: {e}\n(Tip: Ask Gemini to ensure the molecule is valid.)", "red")
             return str(e)

    def execute_load_molecule(self, params):
        """Load a new molecule from SMILES"""
        try:
            smiles = params.get("smiles")
            name = params.get("name", "")
            
            if not smiles:
                self.append_message("System", "Error: No SMILES provided.", "red")
                return
            
            # Use undo-safe loader
            self.load_smiles_undo_safe(smiles)
            
            # Optimize 2D layout
            if hasattr(self.main_window, 'clean_up_2d_structure'):
                self.main_window.clean_up_2d_structure()
            
            # Success message
            if name:
                self.append_message("System", f"{name} loaded successfully.", "green")
            else:
                self.append_message("System", f"Molecule loaded from SMILES: {smiles}", "green")
            
        except Exception as e:
            self.append_message("System", f"Load Error: {e}\n(Tip: The SMILES string might be invalid. Ask Gemini for a more precise one.)", "red")
            return str(e)

    def execute_convert_to_3d(self, params):
        """Trigger 3D conversion and visualization"""
        try:
            self.append_message("System", "Triggering 3D conversion...", "blue")
            QApplication.processEvents()
            
            # Use the internal conversion helper
            self._ensure_main_window_3d_conversion()
            
            self.append_message("System", "3D structure displayed in viewer.", "green")
            
        except Exception as e:
            self.append_message("System", f"3D Conversion Error: {e}\n(Tip: Ask Gemini to try again.)", "red")
            return str(e)

    def execute_clear_canvas(self, params):
        """Clear all molecules from the editor (with undo support)"""
        try:
            mw = self.main_window
            
            # Push undo state BEFORE clearing
            mw.push_undo_state()
            
            # --- Manual Clear 2D Logic ---
            # 1. Clear flags
            if hasattr(mw.scene, 'clear_all_problem_flags'):
                mw.scene.clear_all_problem_flags()
            
            # 2. Remove Items from Scene
            # Note: We must remove from scene before clearing data refs
            for atom_id in list(mw.data.atoms.keys()):
                atom_data = mw.data.atoms[atom_id]
                if atom_data.get('item'):
                    if atom_data['item'].scene() == mw.scene:
                        mw.scene.removeItem(atom_data['item'])
            
            for (id1, id2) in list(mw.data.bonds.keys()):
                bond_data = mw.data.bonds.get((id1, id2))
                if bond_data and bond_data.get('item'):
                   if bond_data['item'].scene() == mw.scene:
                        mw.scene.removeItem(bond_data['item'])

            # 3. Clear Data
            mw.data.atoms.clear()
            mw.data.bonds.clear()
            mw.data.next_atom_id = 1
            
            # 4. Reset helper flags
            mw.is_xyz_derived = False
            if hasattr(mw, 'clear_2d_measurement_labels'):
                mw.clear_2d_measurement_labels()
            
            # --- Manual Clear 3D Logic ---
            mw.plotter.clear()
            mw.current_mol = None
            if hasattr(mw, '_enable_3d_features'):
                mw._enable_3d_features(False)
            
            # Update UI
            mw.has_unsaved_changes = True
            mw.update_undo_redo_actions()
            mw.update_window_title()
            
            # Push undo state AFTER clearing (Saves the "Empty" state on top of stack)
            mw.push_undo_state()
            
            # Update context
            self.check_molecule_change()
            
            self.append_message("System", "Canvas cleared.", "green")
            
        except Exception as e:
            self.append_message("System", f"Clear Error: {e}\n(Tip: Ask Gemini to try again.)", "red")
            return str(e)

    def _get_molecule_for_export(self):
        """Helper to retrieve or generate a 3D molecule for export."""
        current_smiles, _ = self.get_current_molecule_smiles()
        mol_to_export = None
        
        # 1. Try existing Main Window 3D Mol
        if hasattr(self.main_window, 'current_mol') and self.main_window.current_mol:
             if self.main_window.current_mol.GetNumConformers() > 0:
                 mol_to_export = self.main_window.current_mol
        
        # 2. If missing/flat, Trigger Conversion
        if not mol_to_export:
             self.append_message("System", "Triggering 3D conversion for input...", "blue")
             self._ensure_main_window_3d_conversion()
             
             if hasattr(self.main_window, 'current_mol') and self.main_window.current_mol:
                 if self.main_window.current_mol.GetNumConformers() > 0:
                     mol_to_export = self.main_window.current_mol

        # 3. Fallback to Local Generation (Blocking)
        if not mol_to_export:
             status_msg = ""
             try: status_msg = self.main_window.statusBar().currentMessage()
             except: pass
             self.append_message("System", f"Main Window conversion failed (Status: {status_msg}). Trying local generation...", "orange")

             QApplication.processEvents()
             try:
                 mol_3d = Chem.MolFromSmiles(current_smiles)
                 if mol_3d:
                     mol_3d = Chem.AddHs(mol_3d)
                     res = AllChem.EmbedMolecule(mol_3d, AllChem.ETKDGv3())
                     if res == 0:
                         mol_to_export = mol_3d
                     else:
                         self.append_message("Error", f"Local 3D Embed failed (Result code: {res})", "red")
                 else:
                     self.append_message("Error", "Local MolFromSmiles failed (Invalid SMILES?)", "red")
             except Exception as e:
                 self.append_message("Error", f"Failed to generate 3D: {e}", "red")

        return mol_to_export

    def _save_input_file(self, filename, content, filter_str):
        """Helper to save content to a file with dialog."""
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Input File", filename, filter_str)
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write(content)
                self.append_message("System", f"Input file saved to {os.path.basename(file_path)}", "green")
                return None
            except Exception as e:
                self.append_message("Error", f"Failed to save: {e}\n(Tip: Check file permissions.)", "red")
                return str(e)
        else:
            self.append_message("System", "Save cancelled.", "orange")
            return "Cancelled"

    def execute_orca_input_generator(self, params):
        """Generate ORCA Input File"""
        filename = params.get("filename", "input.inp")
        header = params.get("header", "")
        footer = params.get("footer", "")
        
        mol = self._get_molecule_for_export()
        if not mol:
             self.append_message("Error", "Could not generate 3D structure. Input file cancelled.", "red")
             return "No valid molecule"

        try:
            from rdkit.Chem import Descriptors
            charge = Chem.GetFormalCharge(mol)
            num_radicals = Descriptors.NumRadicalElectrons(mol)
            mult = num_radicals + 1

            content = header + "\n\n"
            content += f"* xyz {charge} {mult}\n"
            conf = mol.GetConformer()
            for i, atom in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                content += f"{atom.GetSymbol():2s} {pos.x:10.5f} {pos.y:10.5f} {pos.z:10.5f}\n"
            content += "*\n"
            content += footer + "\n"

            return self._save_input_file(filename, content, "ORCA Input (*.inp *.orca);;All Files (*)")
        except Exception as e:
            self.append_message("Error", f"ORCA Gen Error: {e}", "red")
            return str(e)

    def execute_gaussian_input_generator(self, params):
        """Generate Gaussian Input File"""
        filename = params.get("filename", "input.gjf")
        header = params.get("header", "")
        footer = params.get("footer", "")
        
        mol = self._get_molecule_for_export()
        if not mol:
             self.append_message("Error", "Could not generate 3D structure. Input file cancelled.", "red")
             return "No valid molecule"

        try:
            from rdkit.Chem import Descriptors
            charge = Chem.GetFormalCharge(mol)
            num_radicals = Descriptors.NumRadicalElectrons(mol)
            mult = num_radicals + 1

            content = header + "\n"
            # Gaussian expects Title section + Empty line
            if not content.endswith("\n\n"): 
                if content.endswith("\n"): content += "\n"
                else: content += "\n\n"
                
            content += "Generated by Chat with Molecule Neo\n\n"
            content += f"{charge} {mult}\n"
            
            conf = mol.GetConformer()
            for i, atom in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                content += f"{atom.GetSymbol():2s} {pos.x:10.5f} {pos.y:10.5f} {pos.z:10.5f}\n"
            content += "\n" # End of coords
            content += footer + "\n\n" # Footer often needs blank lines at end in Gaussian

            return self._save_input_file(filename, content, "Gaussian Input (*.gjf *.com);;All Files (*)")
        except Exception as e:
             self.append_message("Error", f"Gaussian Gen Error: {e}", "red")
             return str(e)

    def execute_save_file(self, params):
        """Save arbitrary text content to a file (with [[atom]]/[[atom_count]] injection)"""
        filename = params.get("filename", "output.txt")
        content = params.get("content", "")
        
        # Get Molecule if tags present
        tags = ["[[atom]]", "[[atom_count]]"]
        if any(tag in content for tag in tags):
            mol = self._get_molecule_for_export()
            if not mol:
                self.append_message("Error", "Geometry tag used but no valid molecule found.", "red")
                return "No molecule for XYZ tag"
                
            try:
                conf = mol.GetConformer()
                
                # [[atom_count]] - Total number of atoms
                if "[[atom_count]]" in content:
                    content = content.replace("[[atom_count]]", str(mol.GetNumAtoms()))

                # [[atom]] - Symbol + Coordinates
                if "[[atom]]" in content:
                    lines = []
                    for i, atom in enumerate(mol.GetAtoms()):
                        pos = conf.GetAtomPosition(i)
                        lines.append(f"{atom.GetSymbol():2s} {pos.x:10.5f} {pos.y:10.5f} {pos.z:10.5f}")
                    block = "\n".join(lines)
                    content = content.replace("[[atom]]", block)

            except Exception as e:
                self.append_message("Error", f"Failed to generate geometry block: {e}", "red")
                return str(e)
                
        # Determine filter based on extension
        ext = os.path.splitext(filename)[1]
        filter_str = f"File (*{ext});;All Files (*)" if ext else "All Files (*)"
        
        return self._save_input_file(filename, content, filter_str)


    def _build_context_msg(self, smiles, name, lazy=False):
        """Construct the context message string."""
        name_str = f"Name: {name}. " if name else ""

        # User Request: "for unknown compound. add RDKit descriptors... (3D auto)"
        descriptor_str = ""
        if not lazy and not name: # Only if unknown name?
             descriptor_str = self._get_descriptors_str(smiles)

        # --- User Selection Injection ---
        selection_info = ""
        selected_ids = self.get_selected_atom_indices()
        if selected_ids:
            selection_info = f"\n[User Selection: Atoms {', '.join(selected_ids)} (Matches SMILES Map Numbers)]"
            
        # --- State Injection (Always include if no name) ---
        state_info = ""
        try:
            # Try MainWindow state first
            charge = getattr(self.main_window, 'current_charge', None)
            mult = getattr(self.main_window, 'current_mult', None)
            
            # If not set, calculate from RDKit
            if charge is None or mult is None:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    from rdkit.Chem import Descriptors
                    charge = Chem.GetFormalCharge(mol)
                    num_radicals = Descriptors.NumRadicalElectrons(mol)
                    mult = num_radicals + 1
            
            # Include if we have values (especially useful when name is unknown)
            if charge is not None and mult is not None:
                state_info = f"State: Charge {charge}, Multiplicity {mult}. "
        except Exception:
            pass

        return (
            f"I am currently looking at a molecule. {name_str}{descriptor_str}"
            f"{state_info}SMILES string: {smiles}.{selection_info} "
            f"Please use this as context for our conversation. "
            f"Please do not reply to this information."
        )



    def send_message(self):
        # Force a check of the molecule state right before sending
        # This catches cases where the user loads a molecule and immediately clicks Send
        # before the auto-poll timer triggers.
        self.check_molecule_change()


        text = self.txt_input.toPlainText().strip()
        if not text:
            return

        if not self.chat_session:
            # If session not ready, try init or just warn
            self.append_message("System", "Chat session not initialized. Please wait...", "orange")
            return

        # Show pending info message if any (Context updated, etc.)
        if self.pending_info_text:
            self.append_message("System", self.pending_info_text, "green")
            self.pending_info_text = None

        # UI Updates
        self.append_message("You", text, "blue")
        self.txt_input.clear()
        self.txt_input.setEnabled(False)
        self.btn_send.setEnabled(False)
        self.chat_display.setFocus() # Keep focus away from input while sending

        # Prepare message with pending context if any
        full_text_to_send = text
        
        # Re-verify context freshness (Blocking Update if needed)
        # Even if pending_context_msg exists, it might be lazy (missing descriptors).
        # We need to Ensure Full Context (Name + Descriptors)
        
        # Trigger explicit update check logic
        current_smiles, _ = self.get_current_molecule_smiles()
        if current_smiles:
             # NAME Resolution
             mol_name = self.get_molecule_name(allow_fetch=False) # Check cache
             if not mol_name and self.last_inchikey and PubChemResolver:
                 self.append_message("System", "Resolving molecule name...", "gray")
                 QApplication.processEvents()
                 try:
                     name_res, _ = PubChemResolver.resolve_inchikey_to_name(self.last_inchikey)
                     if name_res:
                         if not hasattr(self, '_name_cache'): self._name_cache = {}
                         self._name_cache[self.last_inchikey] = name_res
                         mol_name = name_res
                 except:
                     pass
             
             # Re-build Context (Lazy=False -> Forces Descriptors -> Forces 3D if needed)
             # This might block slightly for 3D generation, but user requested "Only just before sending".
             
             # User Request: "trigger conversion. non blocking ... but block"
             # Interpretation: Use Main Window's Async Conversion (start via trigger_conversion), but WAIT for it here.
             
             # 1. Trigger Main Window Conversion if needed (and wait)
             self._ensure_main_window_3d_conversion()
             
             # 2. Now mw.current_mol should be updated (if conversion succeeded).
             #    _build_context_msg (via _get_descriptors_str) will pick it up.
             self.pending_context_msg = self._build_context_msg(current_smiles, mol_name, lazy=False)
        
        if self.pending_context_msg:
             # [CONTEXT INJECTION POINT]
             full_text_to_send = f"{self.pending_context_msg}\n\n{text}"
             self.pending_context_msg = None # Clear pending

        # Initialize Streaming UI
        self.start_stream_message("Gemini")

        # LOGGING: Log the full prompt including injected context
        append_log("PROMPT_FULL", full_text_to_send)

        # DEMO MODE CHECK
        if DEMO_MODE:
            response_text = "I'm not sure how to help with that in Demo Mode."
            
            if self.demo_step == 0:
                response_text = (
                    "Welcome! Let's start by loading Benzene using the `load_molecule` tool.\n"
                    "```json\n"
                    '{"tool": "load_molecule", "params": {"smiles": "c1ccccc1", "name": "Benzene"}}\n'
                    "```"
                )
                self.demo_step += 1
            elif self.demo_step == 1:
                response_text = (
                    "Great! Now let's try modifying it. I can replace a Hydrogen with Chlorine using a reaction.\n"
                    "```json\n"
                    '{"tool": "apply_transformation", "params": {"reaction_smarts": "[c:1][H]>>[c:1][Cl]"}}\n'
                    "```"
                )
                self.demo_step += 1
            elif self.demo_step == 2:
                 response_text = (
                    "Now let's highlight the aromatic carbons using atom indices!\n"
                    "```json\n"
                    '{"tool": "highlight_substructure", "params": {"atom_indices": [2, 3, 4, 5, 6, 7]}}\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 3:
                 response_text = (
                    "Now let's set a formal charge on the chlorinated carbon (atom 1)!\n"
                    "```json\n"
                    '{"tool": "set_electronic_state", "params": {"atom_indices": 1, "charge": 1}}\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 4:
                 response_text = (
                    "Let's calculate molecular properties and show 3D structure!\n"
                    "```json\n"
                    '[{"tool": "calculate_descriptors", "params": {"properties": ["MW", "LogP"]}}, {"tool": "convert_to_3d", "params": {}}]\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 5:
                 response_text = (
                    "Now let's generate an **ORCA input file** with proper header order!\n"
                    "```json\n"
                    '{"tool": "orca_input_generator", "params": {"filename": "chlorobenzene-H-opt.inp", "header": "%maxcore 2000\\n%pal nprocs 4 end\\n! B3LYP def2-SVP Opt", "footer": ""}}\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 6:
                 response_text = (
                    "And a **Gaussian input file** too!\n"
                    "```json\n"
                    '{"tool": "gaussian_input_generator", "params": {"filename": "chlorobenzene-H-opt.gjf", "header": "%nprocshared=4\\n%mem=4GB\\n#P B3LYP/6-31G* Opt", "footer": ""}}\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 7:
                 response_text = (
                    "Use the versatile `save_file` tool with atom injection!\n"
                    "```json\n"
                    '{"tool": "save_file", "params": {"filename": "custom_structure.xyz", "content": "[[atom_count]]\\nMy Molecule\\n[[atom]]"}}\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 8:
                 response_text = (
                    "Let's load a new molecule - **Aspirin**!\n"
                    "```json\n"
                    '{"tool": "load_molecule", "params": {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "Aspirin"}}\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 9:
                 response_text = (
                    "Now a **triple combo**: Highlight, calculate ALL properties, AND convert to 3D!\n"
                    "```json\n"
                    '[{"tool": "highlight_substructure", "params": {"atom_indices": [10,11,12]}}, {"tool": "calculate_descriptors", "params": {"properties": ["MW", "LogP", "TPSA", "HBondDonor", "HBondAcceptor", "RotatableBonds"]}}, {"tool": "convert_to_3d", "params": {}}]\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 10:
                 response_text = (
                    "Let's highlight the acetyl group using atom indices!\n"
                    "```json\n"
                    '{"tool": "highlight_substructure", "params": {"atom_indices": [1, 2, 3]}}\n'
                    "```"
                 )
                 self.demo_step += 1
            elif self.demo_step == 11:
                 response_text = (
                    "Finally, let's clear the canvas to start fresh!\n"
                    "```json\n"
                    '{"tool": "clear_canvas", "params": {}}\n'
                    "```"
                 )
                 self.demo_step += 1
            else:
                 response_text = (
                     "That concludes the demo! You've seen:\n"
                     "- Loading molecules (Benzene & Aspirin)\n"
                     "- Chemical transformations (Chlorination)\n"
                     "- Atom highlighting (indices & SMARTS)\n"
                     "- Setting electronic states (charges)\n"
                     "- Multiple tool combinations (up to 3 tools!)\n"
                     "- All descriptor properties\n"
                     "- ORCA & Gaussian input generation\n"
                     "- Custom file saving with atom injection\n"
                     "- Canvas clearing\n\n"
                     "Feel free to try your own prompts now!"
                 )

            QTimer.singleShot(500, lambda: self.on_chunk_received(response_text))
            QTimer.singleShot(600, lambda: self.on_final_response(None))
            return 

        # Worker
        self.worker = GenAIWorker(self.chat_session, full_text_to_send)
        self.worker.chunk_received.connect(self.on_chunk_received)
        self.worker.response_received.connect(self.on_final_response)
        self.worker.error_occurred.connect(self.on_error)
        self.worker.finished.connect(self.on_worker_finished)
        
        self.worker.start()



    def _ensure_main_window_3d_conversion(self):
        """
        Triggers Main Window's 3D conversion (Async) and waits for it to complete (Block locally).
        This ensures 3D descriptors and Context are accurate without freezing the UI thread permanently (uses processEvents).
        Time Limit: 60s (User Request for complex molecules).
        """


        if not hasattr(self.main_window, 'trigger_conversion'):
            return

        # Check if we really need it? (If current_mol exists and has conformers, maybe skip?)
        # User said "trigger conversion", implying they want it FRESH.
        
        # Trigger it
        try:
             self.append_message("System", "Converting 2D structure to 3D... (Max 60s)", "blue")
             QApplication.processEvents()

             import time
             
             # Capture conversion ID state to detect if it actually starts
             prev_conversion_id = getattr(self.main_window, 'next_conversion_id', -1)
             
             try:
                 self.main_window.trigger_conversion()
             except Exception as e:
                 self.append_message("Error", f"Exception during trigger_conversion: {e}", "red")
                 return
             
             # Check if conversion started (ID should increment)
             # Wait a brief moment to ensure UI update propagated
             time.sleep(0.2)
             curr_conversion_id = getattr(self.main_window, 'next_conversion_id', -1)
             
             if prev_conversion_id != -1 and curr_conversion_id == prev_conversion_id:
                 status_msg = "Unknown Error"
                 try:
                     status_msg = self.main_window.statusBar().currentMessage()
                 except:
                     pass
                 self.append_message("Error", f"3D Conversion failed to start: {status_msg}", "red")
                 return
             
             # Wait for start (active_worker_ids populated)
             deadline = time.time() + 60.0 # 60s timeout
             
             # Wait loop
             while time.time() < deadline:
                 QApplication.processEvents() # Keep UI responsive
                 
                 active_ids = getattr(self.main_window, 'active_worker_ids', None)
                 if not active_ids:
                     # Either finished or never started.
                     if time.time() - (deadline - 60.0) < 0.5:
                         # Still in potential startup window
                         time.sleep(0.05)
                         continue
                     else:
                         # Finished
                         break
                 
                 time.sleep(0.05) # Yield CPU
             
             # Force 3D Update as requested
             if hasattr(self.main_window, 'draw_molecule_3d'):
                 self.main_window.draw_molecule_3d(self.main_window.current_mol)
                  
        except Exception as e:
            self.append_message("Error", f"3D Conversion Wait Failed: {e}", "red")

    def start_stream_message(self, sender, color="black"):
        """Prepare UI for streaming a new message"""
        self.stream_accumulated_text = ""
        
        # Insert Header (Gemini:)
        style = "<style>h1, h2, h3 { color: #2c3e50; } a { color: #3498db; text-decoration: none; font-weight: bold; }</style>"
        html_header = f"{style}<b style='color:{color}'>{sender}:</b><br>"
        
        self.chat_display.moveCursor(QTextCursor.MoveOperation.End)
        self.chat_display.insertHtml(html_header)
        
        # Store position where the content starts
        # We add a temporary placeholder or just track the cursor?
        # Tracking the cursor position (int) is safest.
        self.stream_start_pos = self.chat_display.textCursor().position()
        
        # Add a newline for spacing after the message
        self.chat_display.insertHtml("<br><br>")

    def update_stream_message(self):
        """Re-render the accumulated text and update the chat window in-place"""
        processed_html = self.render_content(self.stream_accumulated_text)
        
        # Create a cursor to edit the document
        cursor = self.chat_display.textCursor()
        cursor.setPosition(self.stream_start_pos)
        
        # Select everything from start_pos to End-2 (skipping the <br><br> we added?)
        # Actually, if we just overwrite from start_pos to the end minus the trailing breaks?
        # It's easier to just keep the trailing breaks at the very end.
        
        # Let's try: Select from start_pos to the current end of the block/document?
        # Since we are appending at the end, we can select to End.
        # But we want to preserve the <br><br> if possible, or just re-add them.
        # Let's just re-add them to be simple.
        
        cursor.movePosition(QTextCursor.MoveOperation.End, QTextCursor.MoveMode.KeepAnchor)
        cursor.removeSelectedText()
        
        # Insert new content + Spacing
        cursor.insertHtml(processed_html + "<br><br>")
        
        # Scroll to bottom
        self.chat_display.moveCursor(QTextCursor.MoveOperation.End)

    def on_chunk_received(self, text):
        self.stream_accumulated_text += text
        self.update_stream_message()

    def on_error(self, error_msg):
        """Handle worker errors"""
        self.append_message("System", f"Gemini Error: {error_msg}", "red")
        
        # Log error
        append_log("Error", error_msg)
        
        # Re-enable UI
        self.txt_input.setEnabled(True)
        self.btn_send.setEnabled(True)
        self.txt_input.setFocus()
        
        self.worker = None

    def on_worker_finished(self):
        """Called when worker thread finishes (success or failure)"""
        # Just ensure worker reference is cleared
        # UI re-enabling is handled in on_final_response or on_error
        if self.worker:
            self.worker = None

    def on_final_response(self, response):
        """Called when stream is fully complete"""
        # Ensure final state matches accumulated text
        # (Usually redundant but good for safety)
        # self.update_stream_message() 
        
        # Log to history now that it's complete
        self.chat_history_log.append({"sender": "Gemini", "text": self.stream_accumulated_text})
        append_log("Gemini", self.stream_accumulated_text)
        
        self.log_usage(response)
        
        # --- Check for Tool Calls ---
        # Look for JSON block at the end (supports single object or array)
        text = self.stream_accumulated_text
        json_match = re.search(r'```json\s*([\[{].*?[\]}])\s*```', text, re.DOTALL)
        
        tool_proposed = False
        if json_match:
            try:
                json_str = json_match.group(1)
                payload = json.loads(json_str)
                
                # Support both single tool and multiple tools (array)
                if isinstance(payload, list):
                    # Multiple operations: propose all at once, preserve order
                    self.propose_tool_action({"tools": payload})
                    tool_proposed = True
                elif "tool" in payload and "params" in payload:
                    # Single operation
                    self.propose_tool_action(payload)
                    tool_proposed = True
                    
            except json.JSONDecodeError:
                print("Failed to parse Tool JSON")

        # Re-enable UI only if no tool was proposed (tool prompt freezes input)
        if not tool_proposed:
            self.txt_input.setEnabled(True)
            self.btn_send.setEnabled(True)
            self.txt_input.setFocus()

    # Renamed from on_response to keep old method signature available if needed
    def on_response(self, response):
        pass # Deprecated in favor of streaming flow

    def on_initial_response(self, response):
        """Handle initial context response silently (log only)"""
        text = response.text
        append_log("Gemini (Hidden Init)", text)
        self.log_usage(response)
        # We don't show the initial "Okay I understand" message to the user
        self.append_message("System", "Session initialized. Ready to chat.", "green")

    def log_usage(self, response):
        """Helper to log usage metadata"""
        try:
            if hasattr(response, 'usage_metadata'):
                usage = response.usage_metadata
                append_log("usage", str(usage))
                
            # User Request: "History Management" - Check token/turn count
            # Simple heuristic: If history > N turns, prone.
            if self.chat_session and hasattr(self.chat_session, 'history'):
                history_len = len(self.chat_session.history)
                if history_len > MAX_HISTORY:
                    self._prune_history()
                    
        except Exception as e:
            print(f"Failed to log usage / manage history: {e}")

    def _prune_history(self):
        """Rolling history: Keep last N turns + System prompt context (handled by restart)"""
        try:
            old_history = self.chat_session.history
            # Keep last 10 messages (5 turns)
            new_history = old_history[-10:]
            
            # Restart session with trimmed history
            # Note: start_chat history argument takes list of Content objects.
            model_name = self.settings.get("model", "gemini-flash-latest")
            
            model = genai.GenerativeModel(
                model_name,
                system_instruction=SYSTEM_PROMPT,
                generation_config=GENERATION_CONFIG
            )
            self.chat_session = model.start_chat(history=new_history)
            
            self.append_message("System", "History pruned to save tokens (Rolling update).", "gray")
            append_log("System", "History pruned.")
        except Exception as e:
             append_log("Error", f"Failed to prune history: {e}")



def run(main_window):
    """
    Entry point for the plugin.
    Display the ChatMoleculeWindow.
    """
    try:
        # DEBUG: Confirm entry
        # QMessageBox.information(None, "Debug", "Plugin Starting... Run() called.")
        
        if not hasattr(main_window, 'chat_molecule_window_instance'):
            main_window.chat_molecule_window_instance = None
    
        if main_window.chat_molecule_window_instance:
            try:
                main_window.chat_molecule_window_instance.close()
                main_window.chat_molecule_window_instance = None
            except Exception as e:
                print(f"Error closing existing instance: {e}")
    
        # Create new instance
        # QMessageBox.information(None, "Debug", "Creating Instance...")
        dialog = ChatMoleculeWindow(main_window)
        main_window.chat_molecule_window_instance = dialog
    
        # Use show() for modeless
        # QMessageBox.information(None, "Debug", "Calling show()")
        dialog.show()
        # QMessageBox.information(None, "Debug", "Run() returning")
        
    except BaseException as e:
        import traceback
        traceback.print_exc()
        try:
            QMessageBox.critical(None, "Plugin Critical Error", f"Failed to launch plugin:\n{e}\n\n{traceback.format_exc()}")
        except:
             print(f"CRITICAL: {e}")
