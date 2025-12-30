#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PLUGIN_NAME = "Chat with Molecule Neo (Local)"
PLUGIN_VERSION = "2025.12.30"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Chat with Local LLM (OpenAI-Compatible) about the current molecule. Automatically injects SMILES context. (Neo Version) Note: InChIKey is sent to PubChem."
PLUGIN_ID = "chat_with_molecule_neo_local"
"""


import sys
import os
import json
import io
import base64
import re
import urllib.request
import urllib.parse
import time
import urllib.error
import unicodedata

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
    QFileDialog, QTextBrowser, QPlainTextEdit, QComboBox, QDialog,
    QCheckBox
)
from PyQt6.QtGui import (
    QTextCursor, QColor, QDesktopServices, QAction, QIcon,
    QFont, QTextBlockFormat, QTextCharFormat, QPainter, QGuiApplication
)


class PubChemResolver:
    """Helper class for PubChem API interactions (Embedded)"""
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    # User-Agentを設定してAPIブロックを回避
    HEADERS = {'User-Agent': 'MoleditPy/2025.12.29 (Research Tool)'}
    TIMEOUT = 10  # 秒

    @staticmethod
    def _make_request(url):
        """Helper to create a request with headers"""
        req = urllib.request.Request(url, headers=PubChemResolver.HEADERS)
        return urllib.request.urlopen(req, timeout=PubChemResolver.TIMEOUT)

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
            
            with PubChemResolver._make_request(url) as response:
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


    @staticmethod
    def resolve_name_to_smiles(name):
        """
        指定された化合物名からPubChemを検索し、SMILES文字列を返す（1ステップ検索版）。
        """
        if not name:
            return None, "Empty name provided."

        # 正規化とURLエンコード
        clean_name = unicodedata.normalize('NFKC', name).strip()
        encoded_name = urllib.parse.quote(clean_name)

        # ---------------------------------------------------------
        # 名前から直接SMILESプロパティを取得するURL
        # ---------------------------------------------------------
        # 構造: /compound/name/{name}/property/SMILES/JSON
        url = f"{PubChemResolver.BASE_URL}/compound/name/{encoded_name}/property/SMILES/JSON"

        try:
            # 連続アクセス対策（必要に応じて調整）
            time.sleep(0.3)

            with PubChemResolver._make_request(url) as response:
                # -----------------------------------------------------
                # Case 1: 成功 (HTTP 200)
                # -----------------------------------------------------
                if response.status == 200:
                    raw_json = response.read().decode('utf-8')
                    data = json.loads(raw_json)

                    # APIが混雑中で "Waiting" が返るケースのハンドリング
                    if "Waiting" in data:
                        return None, "PubChem API is busy. Please try again."

                    # プロパティリストの取得
                    props = data.get("PropertyTable", {}).get("Properties", [])
                    
                    if props:
                        # ログにあった通り、キーは単純な "SMILES" になっている
                        compound_data = props[0]
                        if "SMILES" in compound_data:
                            # 成功
                            return compound_data["SMILES"], None
                        elif "CanonicalSMILES" in compound_data:
                             # 念のためCanonicalもチェック
                            return compound_data["CanonicalSMILES"], None
                        elif "IsomericSMILES" in compound_data:
                            return compound_data["IsomericSMILES"], None
                        else:
                            # データはあるがSMILESキーがない
                            return None, f"SMILES data missing in response for '{clean_name}'"
                    else:
                        return None, f"No properties returned for '{clean_name}'"

                # -----------------------------------------------------
                # Case 2: 見つからない (HTTP 404)
                # -----------------------------------------------------
                elif response.status == 404:
                    return None, f"Molecule '{clean_name}' not found."
                
                # その他のステータス
                else:
                    return None, f"API Error: Status {response.status}"

        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, f"Molecule '{clean_name}' not found (404)."
            return None, f"HTTP Error: {e.code} {e.reason}"
        except Exception as e:
            return None, f"Lookup Error: {str(e)}"

# --- Metadata ---
PLUGIN_NAME = "Chat with Molecule Neo (Local)"
PLUGIN_VERSION = "2025.12.30"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Chat with Local LLM (OpenAI-Compatible) about the current molecule. Automatically injects SMILES context. (Neo Version) Note: InChIKey is sent to PubChem."
PLUGIN_ID = "chat_with_molecule_neo_local"

SYSTEM_PROMPT = """You are an expert computational and organic chemistry assistant embedded within the advanced molecular editor software "MoleditPy". 
Your users are researchers, students, or chemistry enthusiasts. Adhere to the following guidelines:

### 1. Protocol for Molecular Identification & Graph Consistency (CRITICAL)
When provided with a SMILES string as `Context` (and optionally a Name), you MUST NOT rely on intuition.
    *   **Graph & Atom Mapping**: The SMILES provided will include **Atom Indices** on HEAVY atoms (e.g., `[C:1]`, `[O:2]`). Hydrogens are implicit.
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
You have access to specific tools to interact with the software. To use a tool, you must output your commands in JSON code blocks.

**Format**:
```json
{
  "tool": "tool_name",
  "params": { ... }
}
```

**Multiple Operations**: You can propose multiple tool calls in a single response. You can use EITHER:
1. A JSON array inside a single code block: `[{"tool": "...", ...}, {"tool": "...", ...}]`
2. Multiple separate JSON code blocks (each block will be collected and executed in order):
```json
{"tool": "apply_transformation", "params": {...}}
```
[Any messages here]
```json
{"tool": "convert_to_3d", "params": {}}
```

**Best Practice**: For any operation that modifies or highlights specific atoms (e.g., `highlight_substructure`, `set_electronic_state`), **use `atom_index`** to specify targets precisely. This avoids SMARTS errors and ensures reliable targeting. 

**Available Tools**:
1.  **`apply_transformation`**: Apply a chemical reaction to the current molecule.
    *   `reaction_smarts`: The Reaction SMARTS string defining the transformation.
    *   **Crucial Warning**: Do **NOT** use `[cH]` (aromatic carbon with implicit H). It fails because `AddHs` makes hydrogens explicit, reducing implicit H count to 0.
    *   **Workaround**: Use `[c]` (matches aromatic carbon regardless of H) or `[#6]`.
    *   **Allowed**: `[c:1][H]` works because it matches the explicit H neighbor.
    *   **Recommended**: `[c:1][H]>>[c:1][Cl]` (Explicitly target H for replacement).
    *   `atom_index`: **REQUIRED**. The Atom Index (Integer) from the SMILES context to apply the transformation to.
    *   Example: `{"tool": "apply_transformation", "params": {"reaction_smarts": "[c:15][H]>>[c:15][Cl]", "atom_index": 15}}`

2.  **`highlight_substructure`**: Visually highlight atoms matching a pattern or indices.
    *   `smarts`: (Optional) The SMARTS pattern to find and highlight.
        *   **Tip**: Use `[#6]` to match ANY Carbon (aliphatic or aromatic). Pattern `C` ONLY matches aliphatic carbons.
    *   `atom_indices`: **REQUIRED**. List of Atom Indices (Integers) to highlight e.g. [1, 5].
    *   Example: `{"tool": "highlight_substructure", "params": {"atom_indices": [1, 2]}}` (Highlight Specific Atoms)

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
    *   `header`: The content to go BEFORE the coordinates. **MUST include** `%pal nprocs N end` for parallel execution and `%maxcore XXXX` for memory per core (in MB), followed by the simple input line (e.g., `! B3LYP...`). **Do NOT** include Charge, Multiplicity/Spin here.
    *   `footer`: The content to go AFTER the coordinates.
    *   Example: `{"tool": "orca_input_generator", "params": {"filename": "benzene-opt.inp", "header": "%pal nprocs 8 end\n%maxcore 2000\n! B3LYP def2-SVP Opt Freq"}}`

5.  **`gaussian_input_generator`**: Create a Gaussian calculation input file.
    *   `filename`: Suggested filename (e.g., "[molecule]-opt.gjf").
    *   `header`: The content to go BEFORE the coordinates. **MUST include** `%nprocshared=N` and `%mem=XGB` directives for parallel execution and memory allocation, followed by the route card (e.g., `#P B3LYP...`). **Do NOT** include Title, Charge, or Multiplicity here. The `%chk` file is auto-generated from the filename.
    *   `footer`: The content to go AFTER the coordinates.
    *   Example: `{"tool": "gaussian_input_generator", "params": {"filename": "benzene-opt.gjf", "header": "%nprocshared=8\n%mem=16GB\n#P B3LYP/6-31G* Opt Freq"}}`

6.  **`load_molecule`**: Load a new molecule from SMILES string (replaces current molecule).
    *   **IMPORTANT**: This tool requires an **exact SMILES string**. If the user asks to load a molecule by **name** (e.g., "Load Aspirin", "カフェインをロードして"), you **MUST** use `load_molecule_by_name` instead. Do NOT generate SMILES yourself - use the PubChem lookup.
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
    *   `atom_index`: **ABSOLUTELY REQUIRED**. The Atom Index (Integer) e.g. `5`. **WILL FAIL WITHOUT THIS.**
    *   `charge`: Integer formal charge (e.g., 1 for +1, -1 for -1).
    *   `multiplicity`: Integer multiplicity.
    *   **CRITICAL**: You MUST ALWAYS include `atom_index`. Without it, the tool does NOTHING.
    *   Example: `{"tool": "set_electronic_state", "params": {"atom_index": 1, "charge": 1}}` (Make Atom 1 +1)

10. **`save_file`**: Save text content to a file with any extension.
    *   `filename`: The filename (e.g., "script.py", "data.csv", "config.ini").
    *   `content`: The text content to write.
    *   **Tags**:
        *   `[[atom]]`: Replaced with "ElementSymbol X Y Z" lines (e.g., "C 0.0 0.0 0.0").
        *   `[[atom_count]]`: Replaced with the total number of atoms.
    *   Example: `{"tool": "save_file", "params": {"filename": "mol.xyz", "content": "[[atom_count]]\nTitle\n[[atom]]"}}`

11. **`load_molecule_by_name`**: **(PREFERRED)** Search PubChem for a molecule by name and load it.
    *   **Use this tool whenever the user mentions a molecule by name** (e.g., "Load Aspirin", "Caffeine を表示", "アデノシンをロードして"). This queries PubChem for the correct structure.
    *   `name`: The common name or IUPAC name of the molecule (e.g., "Aspirin", "Caffeine", "Adenosine").
    *   Example: `{"tool": "load_molecule_by_name", "params": {"name": "Caffeine"}}`

### 3. CRITICAL: Reaction SMARTS & Transformation Strategy
You must adhere to the following rules to prevent "Molecular Destruction" bugs:

1.  **PRESERVE ATOM MAP NUMBERS**:
    * You **MUST** include Atom Map Numbers (e.g., `[C:1]`) in your SMARTS for ALL atoms that persist from reactant to product.
    * **Forbidden**: `[c][H]>>[c][Cl]` (No IDs -> RDKit assumes new atoms -> Old atoms are DELETED).
    * **Required**: `[c:1][H]>>[c:1][Cl]` (IDs preserved -> Atoms are modified in-place).

2.  **ANCHOR TO THE ROOT (NOT JUST THE SELECTION)**:
    * When the user selects an atom to *replace* or *delete* (e.g., a terminal methyl group), do not solely focus on that selected atom as the anchor.
    * **Think about the "Unchanging Root Atom"**: Identify the stable neighbor atom that will remain after the reaction. Use this root atom to anchor your SMARTS pattern.
    * *Example*: Replacing a Methyl group (`[C:31]`) attached to a Ring Carbon (`[C:30]`) with Chlorine.
        * *Risky*: Targeting only index 31.
        * *Safe*: Anchor to the root `[C:30]`. Rule: `[C:30][C:31]>>[C:30][Cl]`.
    * If you must delete an atom, ensure the SMARTS explicitly matches the context (neighbors) to avoid ambiguous matching.

3.  **AVOID `[cH]`**:
    * Never use `[cH]`. Use `[c]` or `[#6]` instead. `[cH]` fails if hydrogens are explicit.

4.  **EXPLICIT HYDROGENS IN REPLACEMENTS**:
    *   When performing a substitution (Replace), **DO NOT OMIT HYDROGENS**.
    *   You must explicitly target the Hydrogen atom being replaced in the reactant.
    *   **Bad**: `[C:1]>>[C:1][Cl]` (relying on implicit valence change is dangerous).
    *   **Good**: `[C:1][H]>>[C:1][Cl]` (Explicity replacing [H] with [Cl] is safe).

### 5. **GENERAL RULE: Heteroatom Substitution & Hydrogen Management**
* **Scope**: This applies to ANY heteroatom substitution (O, N, S, P, etc.).
* **Instruction**: When replacing a skeletal atom (e.g., C) with a Heteroatom, you MUST calculate the valence difference and explicitly remove the excess Hydrogens.
* **Critical Validation**: Ensure your SMARTS explicitly consumes these disappearing Hydrogens.
* **Failure Mode**: Simple symbol swapping (e.g., `[C:8]>>[O:8]`) preserves reactant Hydrogens. This causes "Valence Violations" (e.g., an Ether Oxygen with 2 bonds + 2 inherited Hydrogens → CRASH).
* **Strategy (Anchoring)**: Although the user targets a specific atom (e.g., Atom 8), you MUST anchor the tool execution to a **stable neighbor** (e.g., Atom 7) to strictly define the chemical context.
* **Ether Formation Example**: `{"tool": "apply_transformation", "params": {"reaction_smarts": "[C:7]-[C:8](-[H])(-[H])-[C:9]>>[C:7]-[O:8]-[C:9]", "atom_index": 7}}`
(Scenario: The user selects **Atom 8** to become Oxygen. **Execution**: The `atom_index` is set to the neighbor (**7**) to anchor the reaction, while the SMARTS explicitly maps `[C:7]-[C:8](-[H])(-[H])-[C:9]>>[C:7]-[O:8]-[C:9]` to ensure **Atom 8** correctly transforms into Oxygen.)

### 4. Response Style
* **Tone**: Professional, intellectual, and helpful.
* **Formatting**:
    * Use **SIMPLE** LaTeX format for chemical formulas and math (e.g., $C_{16}H_{10}$).
    * **Avoid** complex commands like `\\xrightarrow` or `\\stackrel`.
    * **CRITICAL**: Always format ANY molecular structure you mention *in plain text* as a clickable link like this: `[Name](smiles:SMILES_STRING)`.
    * **IMPORTANT**: Do **NOT** put SMILES strings inside LaTeX equations.

### 5. MoleditPy Software Context
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

SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "chat_with_molecule_neo_local.json")
LOG_FILE = os.path.join(os.path.dirname(__file__), "chat_with_molecule_neo_local.log")

# --- OpenAI Import ---
# --- Markdown Import ---
try:
    import markdown
    HAS_MARKDOWN = True
except ImportError:
    HAS_MARKDOWN = False

try:
    import openai
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False

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

    def __init__(self, api_key, base_url):
        super().__init__()
        self.api_key = api_key
        self.base_url = base_url

    def run(self):
        try:
            client = openai.OpenAI(api_key=self.api_key, base_url=self.base_url)
            # Fetch models
            all_models = client.models.list()
            # Local APIs (Ollama/LM Studio) might return objects differently, but usually standard list
            # Return ALL models, no filtering for "gpt"
            models = [m.id for m in all_models]
            models.sort()
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

class OpenAIWorker(QThread):
    """Worker thread to handle API calls to avoid freezing UI"""
    response_received = pyqtSignal(str) # Return full text
    chunk_received = pyqtSignal(str)       # Signal for streaming updates
    error_occurred = pyqtSignal(str)

    def __init__(self, client, history, model_name="gpt-4o"):
        super().__init__()
        self.client = client
        self.history = history
        self.model_name = model_name
        self._is_interrupted = False

    def stop(self):
        """Signal the worker to stop"""
        self._is_interrupted = True

    def run(self):
        try:
            # Enable streaming
            # Use request_timeout if possible, but standard client doesn't support easy cancel without breaking connection
            # We will handle loop break.
            stream = self.client.chat.completions.create(
                model=self.model_name,
                messages=self.history,
                stream=True,
                temperature=0.1
            )
            
            full_response = ""
            
            # Iterate through chunks
            for chunk in stream:
                if self._is_interrupted:
                    # Try to close stream if possible (not always exposed in sync client)
                    try: stream.response.close()
                    except: pass
                    break

                if chunk.choices[0].delta.content:
                    content = chunk.choices[0].delta.content
                    full_response += content
                    self.chunk_received.emit(content)
            
            # Emit final complete response
            if not self._is_interrupted:
                self.response_received.emit(full_response)
            else:
                self.error_occurred.emit("Interrupted by User")
        except Exception as e:
            self.error_occurred.emit(str(e))

class ChatMoleculeWindow(QDialog):
    def __init__(self, main_window):
        # User Request: "On top of MoleditPy, but NOT other apps"
        # Solution: Set parent to main_window (Owned windows stay on top of parent in Qt)
        super().__init__(main_window) 
        # Note: Previous "Force parent=None" was preventing this hierarchy. 
        # If clipping occurs, we might need Qt.WindowType.Window
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.Dialog)
        try:
            # QMessageBox.information(None, "Debug", "Window Init Start")
            self.main_window = main_window
            self.setWindowTitle("Chat with Molecule Neo (Local)")
            self.resize(500, 700)
            self.settings = load_settings()
            self.client = None # OpenAI Client
            self.chat_history_state = [] # List of {"role":..., "content":...} for API
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
                 self.lbl_context.setText(f"Context: {self.get_molecule_name(allow_fetch=False) or 'Unknown'} ({current_smiles[:20]}...)")
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

        # API Base URL Row (New)
        base_layout = QHBoxLayout()
        base_layout.addWidget(QLabel("API Base URL:"))
        self.txt_api_base = QLineEdit()
        self.txt_api_base.setText(self.settings.get("api_base", ""))
        self.txt_api_base.setPlaceholderText("e.g. http://localhost:1234/v1")
        base_layout.addWidget(self.txt_api_base)
        settings_layout.addLayout(base_layout)

        # API Key Row
        key_layout = QHBoxLayout()
        key_layout.addWidget(QLabel("API Key:"))
        self.txt_api_key = QLineEdit()
        self.txt_api_key.setEchoMode(QLineEdit.EchoMode.Password)
        self.txt_api_key.setText(self.settings.get("api_key", ""))
        self.txt_api_key.setPlaceholderText("API Key (Optional for Local)")
        key_layout.addWidget(self.txt_api_key)
        
        # Removed "Get API Key" link (User request: No standard OpenAI redirect)
        settings_layout.addLayout(key_layout)

        # Model Row
        model_layout = QHBoxLayout()
        model_layout.addWidget(QLabel("Model:"))
        self.combo_model = QComboBox()
        self.combo_model.setEditable(True) # Allow custom entry
        default_model = self.settings.get("model", "local-model")
        self.combo_model.addItem(default_model)
        self.combo_model.setCurrentText(default_model)
        model_layout.addWidget(self.combo_model, 1)

        btn_fetch_models = QPushButton("Fetch model")
        btn_fetch_models.clicked.connect(self.fetch_models)
        model_layout.addWidget(btn_fetch_models)
        
        # Removed btn_save from here to move it below
        # btn_save = QPushButton("Save & Reload")
        # btn_save.clicked.connect(self.save_settings)
        # model_layout.addWidget(btn_save)

        settings_layout.addLayout(model_layout)
        
        # Enable PubChem Option & Save Button Row
        pubchem_layout = QHBoxLayout()
        
        self.chk_enable_pubchem = QCheckBox("Enable PubChem InChIKey resolve")
        self.chk_enable_pubchem.setChecked(not self.settings.get("disable_pubchem", True))
        #self.chk_enable_pubchem.stateChanged.connect(self.update_warning_label)
        pubchem_layout.addWidget(self.chk_enable_pubchem)
        
        # Spacer to separate checkbox and button? Or just keep them next to each other.
        # pubchem_layout.addStretch() 
        
        # Spacer to separate checkbox and button? Or just keep them next to each other.
        # pubchem_layout.addStretch() 
        
        self.btn_save = QPushButton("Save & Reload")
        self.btn_save.clicked.connect(self.save_settings)
        pubchem_layout.addWidget(self.btn_save)
        
        settings_layout.addLayout(pubchem_layout)

        layout.addWidget(settings_frame)
        
        # --- Export Button ---
        btn_export = QPushButton("Export Chat History")
        btn_export.clicked.connect(self.export_history)
        layout.addWidget(btn_export)

        # --- Warning Label ---
        self.lbl_warning = QLabel()
        self.lbl_warning.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.lbl_warning)
        


        #Connect dynamic update
        #self.txt_api_base.textChanged.connect(self.update_warning_label)
        self.update_warning_label() # Initial call

        # --- Chat Display ---
        self.chat_display = ChatBrowser(self) # Use custom subclass
        self.chat_display.anchorClicked.connect(self.handle_link)
        self.chat_display.setReadOnly(True)
        # Increase base font size
        font = self.chat_display.font()
        font.setPointSize(12)
        self.chat_display.setFont(font)
        layout.addWidget(self.chat_display)

        # --- Thinking / Loading Indicator (New) ---
        self.status_container = QFrame()
        status_layout = QHBoxLayout(self.status_container)
        status_layout.setContentsMargins(0, 0, 0, 0)
        
        self.loading_bar = QProgressBar()
        self.loading_bar.setRange(0, 0) # Indeterminate (Busy)
        self.loading_bar.setFixedHeight(15)
        self.loading_bar.setVisible(False)
        
        self.lbl_thinking = QLabel("Thinking...")
        self.lbl_thinking.setStyleSheet("color: gray; font-style: italic;")
        self.lbl_thinking.setVisible(False)
        
        status_layout.addWidget(self.loading_bar, 1) # Stretch
        status_layout.addWidget(self.lbl_thinking)
        
        layout.addWidget(self.status_container)

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
        
        self.lbl_tool_info = QLabel("The AI wants to execute a command...")
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



        # Connect signals for settings changes
        self.txt_api_base.textChanged.connect(self.check_settings_changed)
        self.txt_api_key.textChanged.connect(self.check_settings_changed)
        self.combo_model.currentTextChanged.connect(self.check_settings_changed)
        self.chk_enable_pubchem.stateChanged.connect(self.check_settings_changed)

        # Initial check to disable button
        self.check_settings_changed()

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



        return super().eventFilter(obj, event)

    def check_settings_changed(self):
        """Enable Save button only if settings differ from saved values"""
        current_api_base = self.txt_api_base.text().strip()
        saved_api_base = self.settings.get("api_base", "")

        current_api_key = self.txt_api_key.text().strip()
        saved_api_key = self.settings.get("api_key", "")

        current_model = self.combo_model.currentText().strip()
        saved_model = self.settings.get("model", "local-model") # default fallback match init

        # Checkbox: Checked means NOT disabled. Saved: disable_pubchem=True/False
        # If checked (True) -> disable_pubchem should be False.
        # Comparison: 
        current_disable_pubchem = not self.chk_enable_pubchem.isChecked()
        saved_disable_pubchem = self.settings.get("disable_pubchem", True)

        has_changes = (
            current_api_base != saved_api_base or
            current_api_key != saved_api_key or
            current_model != saved_model or
            current_disable_pubchem != saved_disable_pubchem
        )

        self.btn_save.setEnabled(has_changes)
        if has_changes:
             self.btn_save.setText("Save & Reload *")
             # User Request: Blue if enabled
             self.btn_save.setStyleSheet("background-color: #0d6efd; color: white; font-weight: bold;")
        else:
             self.btn_save.setText("Save & Reload")
             self.btn_save.setStyleSheet("") # Revert to default (native disabled look)


    def _get_privacy_details(self):
        """Helper: Get privacy text and color based on current settings"""
        url = self.txt_api_base.text().lower().strip()
        
        # 1. Loopback (Strict Local Machine)
        is_loopback = "localhost" in url or "127.0.0.1" in url or "0.0.0.0" in url
        
        # 2. Local Network (Private IP ranges or .local domain)
        is_local_net = ".local" in url or "192.168." in url or "10." in url or "172." in url
        
        # OLD: disable_pubchem = not self.chk_enable_pubchem.isChecked()
        # NEW: Check SAVED settings, not current checkbox state
        disable_pubchem = self.settings.get("disable_pubchem", False)
        
        notice_color = "gray"
        privacy_text = ""

        if is_loopback:
            # Green: Completely local to this machine
            if disable_pubchem:
                notice_color = "green"
                privacy_text = "Privacy Notice: Chat is local. InChIKey is NOT sent to PubChem."
            else:
                notice_color = "orange"
                privacy_text = "Privacy Notice: Chat is local. Only InChIKey is sent to PubChem."

        elif is_local_net:
            # Blue: Local Network
            if disable_pubchem:
                notice_color = "blue"
                privacy_text = "Privacy Notice: Chat is in local network. InChIKey NOT sent to PubChem."
            else:
                notice_color = "orange"
                privacy_text = "Privacy Notice: Chat is in local network. Only InChIKey is sent to PubChem."
                
        elif not url:
            # Empty
            if disable_pubchem:
                privacy_text = "Note: InChIKey is NOT sent to PubChem."
            else:
                notice_color = "orange"
                privacy_text = "Note: InChIKey is sent to PubChem."
        else:
            # Remote / External
            notice_color = "red"
            if disable_pubchem:
                privacy_text = "Privacy Notice: Chat sent to external server. InChIKey NOT sent to PubChem."
            else:
                privacy_text = "Privacy Notice: Chat sent to external server. Only InChIKey is sent to PubChem."
                
        return privacy_text, notice_color

    def update_warning_label(self):
        """Update warning label based on API Base URL privacy status"""
        privacy_text, notice_color = self._get_privacy_details()

        self.lbl_warning.setText(
            "<div style='color: orange; text-align: center;'>"
            "LLM can make mistakes. Please verify important information.<br>"
            f"<span style='font-size: 0.9em; color: {notice_color};'><b>{privacy_text}</b></span>"
            "</div>"
        )

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
        """Fetch available models"""
        if DEMO_MODE:
             self.on_models_fetched(["demo-model"], None)
             return

        api_key = self.txt_api_key.text().strip()
        if not api_key:
             api_key = "lm-studio" # Dummy key for local servers
             
        self.combo_model.setEnabled(False) 
        self.append_message("System", "Fetching models...", "gray")
        
        api_base = self.txt_api_base.text().strip()
        worker = InitWorker(api_key, api_base)
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
             current = self.settings.get("model", "gpt-4o")
             if current in models:
                 self.combo_model.setCurrentText(current)
             # Non-intrusive feedback
             self.append_message("System", f"Found {len(models)} models.", "green")

    def save_settings(self):
        """Save settings to file"""
        api_key = self.txt_api_key.text().strip()
        api_base = self.txt_api_base.text().strip()
        model = self.combo_model.currentText().strip()
        disable_pubchem = not self.chk_enable_pubchem.isChecked()
        
        self.settings["api_key"] = api_key
        self.settings["api_base"] = api_base
        self.settings["model"] = model
        self.settings["disable_pubchem"] = disable_pubchem
        
        # Save to file
        save_settings(self.settings) 
        
        # Re-initialize
        self.initialize_session()
        
        # Force update label just in case
        self.update_warning_label()
        privacy_text, _ = self._get_privacy_details()
        
        # Non-blocking notification
        self.append_message("System", f"Settings saved and session re-initialized.<br><span style='color:gray; font-size:0.9em;'>({privacy_text})</span>", "green")
        
        # Update button state (Should be disabled now)
        self.check_settings_changed()



    def stop_generation(self):
        """Stop current AI generation"""
        if self.worker:
            self.append_message("System", "Stopping generation...", "orange")
            
            # --- FIX: Disconnect signals to prevent zombie updates ---
            try: self.worker.chunk_received.disconnect()
            except: pass
            try: self.worker.response_received.disconnect()
            except: pass
            try: self.worker.error_occurred.disconnect()
            except: pass
            # try: self.worker.finished.disconnect() # Worker is QThread, finished might not be used here but good practice
            # except: pass
            
            self.worker.stop()
            # Worker will emit error "Interrupted by User" which handles UI reset
            # But we force UI reset here just in case to feel responsive
            self.btn_send.setText("Send")
            self.btn_send.setStyleSheet("")
            self.btn_send.clicked.disconnect()
            self.btn_send.clicked.connect(self.send_message)
            self.txt_input.setEnabled(True)
            self.loading_bar.setVisible(False)
            self.lbl_thinking.setVisible(False)

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
            self.client = None
            self.chat_history_state = []
            self.demo_step = 0 # Initialize step counter
            
            self.append_message("System", "<b>DEMO MODE ACTIVE</b><br>No API Key required. Responses are simulated.", "orange")
            # Enable UI
            self.txt_input.setEnabled(True)
            self.btn_send.setEnabled(True)
            self.txt_input.setFocus()
            
            # Trigger first check for molecule
            self.check_molecule_change()
            return

        if not HAS_OPENAI:
            self.append_message("System", "Error: 'openai' library is not installed.", "red")
            self.txt_input.setEnabled(False)
            self.btn_send.setEnabled(False)
            return

        api_base = self.settings.get("api_base", "").strip()
        if not api_base:
            # STRICT Local Mode: Must provide URL. No default fallback.
            self.append_message("System", "Please enter a Local API Base URL above to start.", "orange")
            self.txt_api_base.setFocus()
            return
            
        api_key = self.settings.get("api_key", "").strip()
        if not api_key:
             api_key = "lm-studio" # Dummy key for local servers that don't need one

        self.append_message("System", f"Connecting to {api_base}...", "blue")
        
        # Start background worker
        self.init_worker = InitWorker(api_key, api_base)
        self.init_worker.finished.connect(self.on_init_finished)
        self.init_worker.start()

    def on_init_finished(self, available_models, error_msg):
        if error_msg:
             self.append_message("System", f"Failed to list models: {error_msg}", "red")
             return

        if not available_models:
            self.append_message("System", "No models found at this endpoint.", "red")
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
            # Fallback: Just take the first one available
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
            api_key = self.txt_api_key.text().strip() or "lm-studio"
            api_base = self.txt_api_base.text().strip()
            self.client = openai.OpenAI(api_key=api_key, base_url=api_base)
            self.chat_history_state = [{"role": "system", "content": SYSTEM_PROMPT}]
        except Exception as e:
            self.append_message("System", f"Error initializing Client: {e}", "red")
            self.client = None
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
                    # --- ROBUST RECONSTRUCTION START ---
                    # Build RDKit mol manually to guarantee ID mapping matches UI perfectly
                    try:
                        rwvm = Chem.RWMol()
                        id_map = {} # UI ID -> RDKit Idx
                        
                        # 1. Add Atoms strictly from data keys
                        sorted_ids = sorted(self.main_window.data.atoms.keys())
                        for aid in sorted_ids:
                            atom_data = self.main_window.data.atoms[aid]
                            symbol = atom_data.get('symbol', 'C')
                            a = Chem.Atom(symbol)
                            
                            # Set Properties
                            a.SetIntProp("_original_atom_id", aid)
                            a.SetAtomMapNum(aid + 1) # KEY FIX: Set MapNum here directly!
                            
                            if 'charge' in atom_data: a.SetFormalCharge(atom_data['charge'])
                            if 'explicit_valence' in atom_data: a.SetNumExplicitHs(atom_data['explicit_valence'])
                                
                            idx = rwvm.AddAtom(a)
                            id_map[aid] = idx

                        # 2. Add Bonds
                        if hasattr(self.main_window.data, 'bonds') and self.main_window.data.bonds:
                            sorted_bids = sorted(self.main_window.data.bonds.keys(), key=str)
                            for bid in sorted_bids:
                                bond_data = self.main_window.data.bonds[bid]
                                
                                # PRIORITIZE bond_data values (Source->Target) to preserve direction for stereo
                                id1 = bond_data.get('atom1')
                                id2 = bond_data.get('atom2')
                                
                                # Fallback to key if data is missing (legacy compat)
                                if id1 is None or id2 is None:
                                     if isinstance(bid, tuple) and len(bid) == 2:
                                         id1, id2 = bid

                                if id1 in id_map and id2 in id_map:
                                    bond_data = self.main_window.data.bonds[bid]
                                    order = bond_data.get('order', 1)
                                    stereo_val = bond_data.get('stereo', 0)
                                    
                                    # Translate order int to RDKit BondType
                                    btype = Chem.BondType.SINGLE
                                    if order == 2: btype = Chem.BondType.DOUBLE
                                    elif order == 3: btype = Chem.BondType.TRIPLE
                                    
                                    bond_idx = rwvm.AddBond(id_map[id1], id_map[id2], btype)
                                    
                                    # Set Stereo for Single Bonds (Wedge/Dash)
                                    if order == 1 and stereo_val > 0:
                                        # Assuming mw.data.bonds uses 1=Wedge, 2=Dash (matches scene.create_bond stereo param)
                                        # However, RDKit requires the bond to be set on the specific bond object
                                        # We need to retrieve the bond we just added.
                                        new_bond = rwvm.GetBondWithIdx(bond_idx - 1)
                                        if stereo_val == 1:
                                            new_bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                                        elif stereo_val == 2:
                                            new_bond.SetBondDir(Chem.BondDir.BEGINDASH)
                        
                        mol = rwvm.GetMol()
                        # Sanitize but keep our MapNums
                        pass # Don't rely on external to_rdkit_mol
                        
                    except Exception as e:
                        print(f"Local Mol Build Failed: {e}")
                        # Fallback to existing method if local build fails
                        mol = self.main_window.data.to_rdkit_mol()
                    
                    if mol:
                        # Map internal IDs to MapNums if not already set (fallback)
                        for atom in mol.GetAtoms():
                             if atom.GetAtomMapNum() == 0:
                                 if atom.HasProp("_original_atom_id"):
                                     aid = atom.GetIntProp("_original_atom_id")
                                     atom.SetAtomMapNum(aid + 1)
                                 else:
                                     atom.SetAtomMapNum(atom.GetIdx() + 1)
                    # --- ROBUST RECONSTRUCTION END ---
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
            # Assign Stereochemistry from Bond Directions (Critical for @/@@ in SMILES)
            Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

            # Remove hydrogens for AI prompt, but preserve MapNums on heavy atoms
            mol_clean = Chem.RemoveHs(mol)
            
            # Ensure every heavy atom has a MapNum
            # for atom in mol_clean.GetAtoms():
            #     if atom.GetAtomMapNum() == 0:
            #         atom.SetAtomMapNum(atom.GetIdx() + 1)
                    
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
        
            # Determine strict ID mapping if Map Numbers are present
            map_num_to_rd_idx = {}
            for i in range(mol.GetNumAtoms()):
                 atom = mol.GetAtomWithIdx(i)
                 m = atom.GetAtomMapNum()
                 if m > 0: map_num_to_rd_idx[m] = i

            # Create Atoms
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()
                map_num = atom.GetAtomMapNum()
                
                relative_x = pos.x - mol_center_x
                relative_y = pos.y - mol_center_y
                
                scene_x = (relative_x * SCALE_FACTOR) + view_center.x()
                scene_y = (-relative_y * SCALE_FACTOR) + view_center.y()
                
                # Access mw.scene directly
                atom_id = mw.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge)
                
                # --- STRICT ID REMAPPING (FIX) ---
                # Ensure ID = MapNum - 1 if MapNum exists
                if map_num > 0:
                    target_id = map_num - 1
                    if atom_id != target_id:
                        if target_id not in mw.data.atoms:
                            # Remap
                            try:
                                adata = mw.data.atoms.pop(atom_id)
                                mw.data.atoms[target_id] = adata
                                atom_id = target_id
                                # Update item internal ID if accessible (optional/risky if private)
                                if 'item' in adata and hasattr(adata['item'], 'atom_id'):
                                    adata['item'].atom_id = target_id
                            except Exception as e:
                                print(f"ID Remap Failed: {e}")
                        else:
                            print(f"ID Collision: {target_id} already exists. Keeping {atom_id}")
                
                rdkit_idx_to_my_id[i] = atom_id
                
                # Update next_atom_id to avoid future collisions
                if hasattr(mw.data, 'next_atom_id'):
                     if atom_id >= mw.data.next_atom_id:
                          mw.data.next_atom_id = atom_id + 1

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
        # User Option: Disable PubChem
        if self.settings.get("disable_pubchem", False):
            return None
            
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
            # Get current SMILES for label
            current_smiles = self.last_smiles or ""
            smiles_preview = f"({current_smiles[:20]}...)" if current_smiles else ""
            display_name = name if name else "Unknown"
            self.lbl_context.setText(f"Context: {display_name} {smiles_preview}")
            
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
            elif tool_name == "load_molecule_by_name":
                desc += f" (Search & Load: {params.get('name', 'Unknown')})"
            elif tool_name == "set_electronic_state":
                desc += f" (Charge: {params.get('charge', '-')}, Mult: {params.get('multiplicity', '-')})"
            elif tool_name == "convert_to_3d":
                desc += " (Show 3D structure)"
            elif tool_name == "clear_canvas":
                desc += " (Clear all)"
            
        self.lbl_tool_info.setText(f"{self.combo_model.currentText()} suggests: {desc}")
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
            self.btn_tool_accept.setText(f"Retry Step (Ask {self.combo_model.currentText()})")
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
        self.txt_input.setFocus()  # Return focus to input (avoid API key selection)
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
            self.btn_tool_accept.setText(f" Retry (Ask {self.combo_model.currentText()}) ")
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
        elif tool_name == "load_molecule_by_name":
            return self.execute_load_molecule_by_name(params)
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
        """Ask AI to retry based on last error"""
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
             
             # --- [FIX 1] Track Max Map ID ---
             max_map_num = 0
             for atom in mol.GetAtoms():
                 m = atom.GetAtomMapNum()
                 if m > max_map_num:
                     max_map_num = m

             # Reaction Logic
             try:
                rxn = AllChem.ReactionFromSmarts(reaction_smarts)
             except Exception as e:
                 self.append_message("System", f"Invalid SMARTS: {e}", "red")
                 return
             


             # Explicit Hydrogens handling
             mol_with_h = Chem.AddHs(mol)
             
             # DEBUG: Check if AddHs preserves maps
             # maps_h = [a.GetAtomMapNum() for a in mol_with_h.GetAtoms()]
             # print(f"DEBUG: MapNums after AddHs: {maps_h}")
             
             # Run reaction
             products = rxn.RunReactants((mol_with_h,))
             
             if not products:
                  # Retry without explicit Hs
                  self.append_message("System", "Reaction produced no products with explicit Hs. Checking implicit...", "orange")
                  products = rxn.RunReactants((mol,))
                  if not products:
                      self.append_message("System", "Reaction produced no products. Check compatibility.", "orange")
                      return

        
             # --- FILTER PRODUCTS BY ATOM INDEX ---
             target_index = params.get("atom_index")
             # Fallback for old param
             if target_index is None and params.get("atom_indices"):
                  val = params.get("atom_indices")
                  if isinstance(val, list) and len(val) > 0: target_index = val[0]

             selected_product_idx = 0
             
             if target_index is not None:
                 try:
                     # Convert MapNums (string/int) -> Int
                     target_mapnum = int(target_index)
                     
                     # Find matches in Reactant (mol_with_h used for reaction)
                     # Reaction template is reactant 0
                     reactant_template = rxn.GetReactants()[0]
                     matches = mol_with_h.GetSubstructMatches(reactant_template, uniquify=False)
                     
                     
                     # Iterate matches to find ALL matches containing our target MapNum
                     # [FIX] Strict Map Enforcement: Filter out matches that violate template map numbers
                     # 1. Build Template Map Constraints
                     template_constraints = {}
                     for atom in reactant_template.GetAtoms():
                         m = atom.GetAtomMapNum()
                         if m > 0:
                             template_constraints[atom.GetIdx()] = m
                     
                     candidate_indices = []
                     for i, match_indices in enumerate(matches):
                         # match_indices contains RDKit Atom Indices
                         
                         # Check Template Map Constraints
                         is_valid_map_match = True
                         for tmpl_idx, mol_idx in enumerate(match_indices):
                             if tmpl_idx in template_constraints:
                                 required_map = template_constraints[tmpl_idx]
                                 mol_atom = mol_with_h.GetAtomWithIdx(mol_idx)
                                 actual_map = mol_atom.GetAtomMapNum()
                                 if actual_map != required_map:
                                     is_valid_map_match = False
                                     break
                         
                         if not is_valid_map_match:
                             continue # Skip this match
                         
                         
                         for idx in match_indices:
                             atom = mol_with_h.GetAtomWithIdx(idx)
                             if atom.GetAtomMapNum() == target_mapnum:
                                 # Found a match anchored to our atom.
                                 # Calculate atom count of the resulting product to determine "preservation"
                                 if i < len(products):
                                     prod_mol = products[i][0]
                                     atom_count = prod_mol.GetNumAtoms()
                                     candidate_indices.append((i, atom_count))
                                 break
                     
                     if candidate_indices:
                         # Sort candidates by Atom Count (Descending) -> Maximizes retention
                         candidate_indices.sort(key=lambda x: x[1], reverse=True)
                         
                         best_idx, best_count = candidate_indices[0]
                         selected_product_idx = best_idx
                         
                         self.append_message("System", f"Selected match #{best_idx} (Atoms: {best_count}) from {len(candidate_indices)} candidates.", "gray")
                     else:
                         self.append_message("System", "Warning: Selected atom not found in reaction matches. Applying to first match.", "orange")
                         
                 except Exception as e:
                     print(f"Index filtering failed: {e}")

             # Take selected product
             if selected_product_idx < len(products):
                 new_mol = products[selected_product_idx][0]
             else:
                 new_mol = products[0][0]

             # ============================================================
             # 【Correction】 Post-Check: Verify Transformation Integrity
             # Prevent corrupted IDs or massive atom loss.
             # ============================================================

             # Check 1: Massive Atom Loss (Structure collapsed?)
             orig_count = mol.GetNumAtoms()
             new_count = new_mol.GetNumAtoms()
             if orig_count > 5 and new_count < orig_count * 0.7:
                 # If >30% atoms lost, assume destruction
                 error_msg = f"Safety Guard: Transformation caused massive atom loss ({orig_count} -> {new_count}). Aborted."
                 self.append_message("System", error_msg, "red")
                 return error_msg

             # Check 2: Duplicate ID Check (Are there colliding IDs?)
             seen_ids = set()
             duplicates = []
             for atom in new_mol.GetAtoms():
                 mid = atom.GetAtomMapNum()
                 if mid > 0:
                     if mid in seen_ids:
                         duplicates.append(str(mid))
                     seen_ids.add(mid)
             
             if duplicates:
                 # Check strictness: Abort if duplicates found
                 error_msg = f"Safety Guard: Transformation created duplicate Atom IDs ({', '.join(duplicates[:3])}...). Aborted."
                 self.append_message("System", error_msg, "red")
                 return error_msg
             
             # ============================================================
             # Post-Check Passed
             # ============================================================
             
             # ROBUST FIX: Implicit Hydrogens & Valence
             try: new_mol.UpdatePropertyCache(strict=False)
             except: pass
             
             # Sanitize
             try: Chem.SanitizeMol(new_mol)
             except: pass
             
             # Remove hydrogens
             try: new_mol = Chem.RemoveHs(new_mol, implicitOnly=False, updateExplicitCount=True, sanitize=True)
             except: pass

             # --- [FIX 2] Assign new IDs to new atoms (MapNum=0) ---
             # Essential for update_structure_diff_based to create bonds
             for atom in new_mol.GetAtoms():
                 if atom.GetAtomMapNum() == 0:
                     max_map_num += 1
                     atom.SetAtomMapNum(max_map_num)

             # Enforce Stereo assignment on product before SMILES generation
             try: Chem.AssignStereochemistry(new_mol, force=True, cleanIt=True)
             except: pass

             # 4. SMILES Round-Trip
             temp_smiles = Chem.MolToSmiles(new_mol)
             clean_mol = Chem.MolFromSmiles(temp_smiles)
             
             if not clean_mol:
                 clean_mol = new_mol
             
             # 5. Optimize 2D Coordinates (Fresh Layout)
             # UX IMPROVEMENT: Try to match original coordinates to prevent jumping
             try:
                 ref_mol = Chem.MolFromSmiles(current_smiles)
                 if ref_mol:
                     AllChem.Compute2DCoords(ref_mol) # Ensure ref has coords
                     # Matches topology to ref
                     AllChem.GenerateDepictionMatching2DStructure(clean_mol, ref_mol)
                 else:
                     raise Exception("Ref mol failed")
             except:
                 # Fallback if matching fails (e.g. significant structural change)
                 AllChem.Compute2DCoords(clean_mol)
             
             # Load back
             final_smiles = Chem.MolToSmiles(clean_mol)
             
             # Use local differential updater for smoother UX
             self.update_structure_diff_based(final_smiles)
             
             # User Request: "Optimize 2D" after conversion
             if hasattr(self.main_window, 'clean_up_2d_structure'):
                 self.main_window.clean_up_2d_structure()
                 
             self.append_message("System", f"Transformation Applied.\nRule: `{reaction_smarts}`", "green")
                 
        except Exception as e:
            self.append_message("System", f"Transformation Failed: {e}\n(Tip: The SMARTS pattern might be invalid. Ask {self.combo_model.currentText()} for a more precise one.)", "red")
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
                 # Handle single int or list
                 if isinstance(atom_indices_param, int):
                     needed_maps = [atom_indices_param]
                 else:
                     needed_maps = [int(x) for x in atom_indices_param]
                 
                 for atom in mol.GetAtoms():
                     if atom.GetAtomMapNum() in needed_maps:
                         atom_indices_set.add(atom.GetIdx())
            
            # Fallback: Check if atom_index was passed (legacy/hallucination compat)
            elif params.get("atom_index") is not None:
                  needed_map = int(params.get("atom_index"))
                  for atom in mol.GetAtoms():
                     if atom.GetAtomMapNum() == needed_map:
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
            
            # A) Explicit Index provided (Map Number)
            atom_index_param = params.get("atom_index")
            # Fallback
            if atom_index_param is None and params.get("atom_indices"):
                 val = params.get("atom_indices")
                 if isinstance(val, list) and len(val) > 0: atom_index_param = val[0]

            if atom_index_param is not None:
                 needed_map = int(atom_index_param)
                 for atom in mol.GetAtoms():
                     if atom.GetAtomMapNum() == needed_map:
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
                          selected_map_nums.append(aid) # MapNum is ID + 1
                 
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
             self.append_message("System", f"Calculation Error: {e}\n(Tip: Ask {self.combo_model.currentText()} to ensure the molecule is valid.)", "red")
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
            self.append_message("System", f"Load Error: {e}\n(Tip: The SMILES string might be invalid. Ask {self.combo_model.currentText()} for a more precise one.)", "red")
            return str(e)

    def execute_load_molecule_by_name(self, params):
        """Search PubChem by name and load the molecule"""
        try:
            name = params.get("name")
            if not name:
                self.append_message("System", "Error: No name provided for search.", "red")
                return "No name provided"

            self.append_message("System", f"Searching PubChem for '{name}'...", "blue")
            QApplication.processEvents()  # UI更新

            # PubChem検索を実行
            smiles, error = PubChemResolver.resolve_name_to_smiles(name)

            if error:
                self.append_message("System", f"Search Failed: {error}", "red")
                return error
            
            if not smiles:
                self.append_message("System", "No SMILES returned from PubChem.", "red")
                return "No SMILES found"

            self.append_message("System", f"Found: {smiles}", "gray")

            # 既存の安全なローダーを再利用して読み込み
            self.load_smiles_undo_safe(smiles)
            
            # 2D構造の整形
            if hasattr(self.main_window, 'clean_up_2d_structure'):
                self.main_window.clean_up_2d_structure()

            self.append_message("System", f"Loaded '{name}' successfully.", "green")

        except Exception as e:
            self.append_message("System", f"Load by Name Error: {e}", "red")
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
            self.append_message("System", f"3D Conversion Error: {e}\n(Tip: Ask {self.combo_model.currentText()} to try again.)", "red")
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
            self.append_message("System", f"Clear Error: {e}\n(Tip: Ask {self.combo_model.currentText()} to try again.)", "red")
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

            # Auto-generate %chk line from filename
            chk_basename = os.path.splitext(os.path.basename(filename))[0]
            content = f"%chk={chk_basename}.chk\n"
            content += header + "\n"
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


    def get_selected_atom_indices(self):
        """Get list of selected atom IDs (as strings) matching SMILES MapNums."""
        indices = []
        if hasattr(self.main_window, 'data') and hasattr(self.main_window.data, 'atoms'):
             for aid, adata in self.main_window.data.atoms.items():
                 item = adata.get('item')
                 if item and item.isSelected():
                     # Match the new MapNum logic: MapNum = ID + 1
                     indices.append(str(aid + 1))
        return indices

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

        if not self.client:
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
             # Fix: respect disable_pubchem setting
             disable_pubchem = self.settings.get("disable_pubchem", False)
             if not mol_name and self.last_inchikey and PubChemResolver and not disable_pubchem:
                 self.append_message("System", "Resolving molecule name...", "gray")
                 QApplication.processEvents()
                 try:
                     name_res, _ = PubChemResolver.resolve_inchikey_to_name(self.last_inchikey)
                     if name_res:
                         if not hasattr(self, '_name_cache'): self._name_cache = {}
                         self._name_cache[self.last_inchikey] = name_res
                         mol_name = name_res
                         # Update context label immediately
                         smiles_preview = f"({current_smiles[:20]}...)" if current_smiles else ""
                         self.lbl_context.setText(f"Context: {mol_name} {smiles_preview}")
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
        self.start_stream_message(self.combo_model.currentText())

        # LOGGING: Log the full prompt including injected context
        append_log("PROMPT_FULL", full_text_to_send)
        
        # --- UI UPDATE FOR STOP/LOADING ---
        self.lbl_thinking.setVisible(True)
        self.loading_bar.setVisible(True)
        
        # Change Send -> Stop
        self.btn_send.setText("Stop")
        self.btn_send.setStyleSheet("background-color: #dc3545; color: white; font-weight: bold;") # Red
        self.btn_send.clicked.disconnect()
        self.btn_send.clicked.connect(self.stop_generation)
        
        self.btn_send.setEnabled(True) # Ensure enabled so it can be clicked
        self.txt_input.setEnabled(False) # Still disable input

        # SHORTCUTS handled above... but if shortcut returns, we need to reset UI? 
        # Actually shortcuts mock the stream, so they should probably respect the UI flow too or handle their own reset.
        # But for safety, I'll let them run. Interrupting them isn't prioritized since they are fast/fake.

        # SHORTCUT: "br" triggers bromination demo sequence for SELECTED atoms
        if DEMO_MODE and text.lower().startswith("br"):
             # Get selected atom IDs (strings)
             selected_ids = self.get_selected_atom_indices()
             
             tools_json_list = []
             if selected_ids:
                 self.append_message("System", f"Bromination Shortcut: Applying to selected atoms {', '.join(selected_ids)}", "green")
                 for aid in selected_ids:
                     # aid is already 1-based MapNum (string)
                     tools_json_list.append(
                         f'  {{"tool": "apply_transformation", "params": {{"reaction_smarts": "[#6:{aid}][H]>>[#6:{aid}][Br]", "atom_index": {aid}}}}}'
                         #f'  {{"tool": "apply_transformation", "params": {{"reaction_smarts": "[c:{aid}][H]>>[c:{aid}][Br]", "atom_index": {aid}}}}}'
                     )
                 json_body = "[\n" + ",\n".join(tools_json_list) + "\n]"
             else:
                 # No selection - abort
                 self.append_message("System", "Bromination Shortcut: No atoms selected. Please select atoms to brominate.", "orange")
                 # Reset UI since we are returning early
                 self.btn_send.setText("Send")
                 self.btn_send.setStyleSheet("")
                 self.btn_send.clicked.disconnect()
                 self.btn_send.clicked.connect(self.send_message)
                 self.txt_input.setEnabled(True)
                 self.loading_bar.setVisible(False)
                 self.lbl_thinking.setVisible(False)
                 return

             response_text = (
                 "Applying Bromination Demo Sequence...\n"
                 "```json\n"
                 f"{json_body}\n"
                 "```"
             )
             
             # UI Reset Wrapper for Demo
             def reset_ui_demo(response):
                 self.on_final_response(response)
                 self.btn_send.setText("Send")
                 self.btn_send.setStyleSheet("")
                 self.loading_bar.setVisible(False)
                 self.lbl_thinking.setVisible(False)
                 
             QTimer.singleShot(500, lambda: self.on_chunk_received(response_text))
             QTimer.singleShot(600, lambda: reset_ui_demo(None))
             return

        # SHORTCUT: "load" triggers bromination demo sequence for SELECTED atoms
        if DEMO_MODE and text.lower().startswith("load"):
            # Get selected atom IDs (strings)
            selected_ids = self.get_selected_atom_indices()
             
            tools_json_list = []
            self.append_message("System", f"p-Xylene Shortcut: Loading p-xylene to 2D editor", "green")
            # aid is already 1-based MapNum (string)
            tools_json_list.append(
                f'  {{"tool": "load_molecule_by_name", "params": {{"name": "p-xylene"}}}}'
            )
            json_body = "[\n" + ",\n".join(tools_json_list) + "\n]"

            response_text = (
                "Loading p-xylene molecule...\n"
                "```json\n"
                f"{json_body}\n"
                "```"
            )
            
            # UI Reset Wrapper
            def reset_ui_demo(response):
                 self.on_final_response(response)
                 self.btn_send.setText("Send")
                 self.btn_send.setStyleSheet("")
                 self.loading_bar.setVisible(False)
                 self.lbl_thinking.setVisible(False)

            QTimer.singleShot(500, lambda: self.on_chunk_received(response_text))
            QTimer.singleShot(600, lambda: reset_ui_demo(None))
            return

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
                    '{"tool": "set_electronic_state", "params": {"atom_index": 1, "charge": 1}}\n'
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

            # UI Reset Wrapper
            def reset_ui_demo(response):
                 self.on_final_response(response)
                 self.btn_send.setText("Send")
                 self.btn_send.setStyleSheet("")
                 self.loading_bar.setVisible(False)
                 self.lbl_thinking.setVisible(False)

            QTimer.singleShot(500, lambda: self.on_chunk_received(response_text))
            QTimer.singleShot(600, lambda: reset_ui_demo(None))
            return

        # Append to History State
        self.chat_history_state.append({"role": "user", "content": full_text_to_send})

        # Worker
        model = self.settings.get("model", "gpt-4o")
        self.worker = OpenAIWorker(self.client, self.chat_history_state, model_name=model)
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
        
        # Insert Header (Model Name:)
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
        # Hide thinking indicators on first chunk
        if self.lbl_thinking.isVisible():
            self.lbl_thinking.setVisible(False)
            # We keep the loading/busy bar active to show it's still streaming/working?
            # User said "loading icon while model responding".
            # Usually spinners stop when text starts arriving, OR they persist until done.
            # I will persist it until final_response.
            pass

    def on_error(self, error_msg):
        """Handle worker errors"""
        self.append_message("System", f"{self.combo_model.currentText()} Error: {error_msg}", "red")
        
        # Log error
        append_log("Error", error_msg)
        
        # Reset UI (Stop -> Send)
        self.btn_send.setText("Send")
        self.btn_send.setStyleSheet("")
        try: self.btn_send.clicked.disconnect() 
        except: pass
        self.btn_send.clicked.connect(self.send_message)
        
        self.loading_bar.setVisible(False)
        self.lbl_thinking.setVisible(False)
        
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
        self.chat_history_log.append({"sender": self.combo_model.currentText(), "text": self.stream_accumulated_text})
        append_log(self.combo_model.currentText(), self.stream_accumulated_text)

        # Append to context history
        self.chat_history_state.append({"role": "assistant", "content": self.stream_accumulated_text})
        
        # self.log_usage(response) # Usage logging might differ for OpenAI
        
        # --- Check for Tool Calls ---
        # Look for ALL JSON blocks (supports multiple separate JSON objects)
        text = self.stream_accumulated_text
        json_matches = re.findall(r'```json\s*([\[{].*?[\]}])\s*```', text, re.DOTALL)
        
        tool_proposed = False
        if json_matches:
            all_tools = []
            for json_str in json_matches:
                try:
                    payload = json.loads(json_str)
                    
                    # Collect all tools from all JSON blocks
                    if isinstance(payload, list):
                        # Array of tools
                        all_tools.extend(payload)
                    elif "tool" in payload and "params" in payload:
                        # Single tool
                        all_tools.append(payload)
                except json.JSONDecodeError:
                    print(f"Failed to parse Tool JSON: {json_str[:50]}...")
            
            # Propose all collected tools at once
            if all_tools:
                if len(all_tools) == 1:
                    self.propose_tool_action(all_tools[0])
                else:
                    self.propose_tool_action({"tools": all_tools})
                tool_proposed = True

        # Re-enable UI only if no tool was proposed (tool prompt freezes input)
        if not tool_proposed:
            self.txt_input.setEnabled(True)
            self.btn_send.setEnabled(True)
            self.txt_input.setFocus()

        # Reset UI (Stop -> Send)
        self.btn_send.setText("Send")
        self.btn_send.setStyleSheet("")
        try: self.btn_send.clicked.disconnect() 
        except: pass
        self.btn_send.clicked.connect(self.send_message)
        
        self.loading_bar.setVisible(False)
        self.lbl_thinking.setVisible(False)

    # Renamed from on_response to keep old method signature available if needed
    def on_response(self, response):
        pass # Deprecated in favor of streaming flow

    def on_initial_response(self, response):
        """Handle initial context response silently (log only)"""
        text = response.text
        append_log(f"{self.combo_model.currentText()} (Hidden Init)", text)
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
            model_name = self.settings.get("model", "gpt-4o")
            
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




    def update_structure_diff_based(self, smiles_string):
        """
        SMILESの内容に基づいて、現在のシーン上の原子を「更新」する。
        全消去せず、IDが一致するものは座標と属性を変更する。
        """
        mw = self.main_window
        try:
            cleaned_smiles = smiles_string.strip()
            mol = Chem.MolFromSmiles(cleaned_smiles)
            if not mol: return

            # --- [IMPORTANT] Kekulize ---
            # Explicitly clear aromatic flags to ensure aromatic bonds (1.5) become localized double bonds (2.0)
            # This prevents them from being cast to single bonds (1.0) during int conversion for visualization.
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            except Exception as e:
                print(f"Kekulize Warning: {e}")

            # 1. 座標生成 (2D)
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()

            # --- 座標合わせ (Alignment) ---
            # RDKitの座標系を、現在のキャンバスの座標系に合わせるためのオフセットを計算
            
            # A. 現在のキャンバス上の選択原子（または全原子）の重心を計算
            current_points = []
            if mw.data.atoms:
                for ad in mw.data.atoms.values():
                    if ad.get('item'):
                        current_points.append(ad['item'].scenePos())
            
            if current_points:
                avg_x = sum(p.x() for p in current_points) / len(current_points)
                avg_y = sum(p.y() for p in current_points) / len(current_points)
                target_center = QPointF(avg_x, avg_y)
            else:
                # 原子がない場合はビューの中心
                if hasattr(mw, 'view_2d') and mw.view_2d:
                     target_center = mw.view_2d.mapToScene(mw.view_2d.viewport().rect().center())
                else:
                     target_center = QPointF(0, 0)

            # --- PATCH: Ensure next_atom_id exists ---
            # Fix for AttributeError: 'MolecularData' object has no attribute 'next_atom_id'
            if not hasattr(mw.data, 'next_atom_id'):
                if mw.data.atoms:
                    # simplistic max+1 (safe for unique ID generation)
                    mw.data.next_atom_id = max(mw.data.atoms.keys()) + 1
                else:
                    mw.data.next_atom_id = 0

            # B. RDKit側の重心を計算
            rd_positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            if rd_positions:
                rd_avg_x = sum(p.x for p in rd_positions) / len(rd_positions)
                rd_avg_y = sum(p.y for p in rd_positions) / len(rd_positions)
            else:
                rd_avg_x, rd_avg_y = 0.0, 0.0

            SCALE = 50.0  # RDKit -> Sceneのスケール係数

            # --- Undo State Push ---
            mw.push_undo_state()

            # --- 2. 原子の更新 / 新規作成 ---
            
            # 今回のSMILESに含まれる MapNum (ID+1) を記録
            processed_ids = set()
            
            # MapNumの最大値を追跡（新規追加用）


            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                
                # 座標変換: (RDKit座標 - RDKit重心) * スケール + キャンバス重心
                scene_x = ((pos.x - rd_avg_x) * SCALE) + target_center.x()
                scene_y = (-(pos.y - rd_avg_y) * SCALE) + target_center.y() # Y軸反転
                
                map_num = atom.GetAtomMapNum()
                
                target_id = -1
                
                # A) 既存IDの更新 (MapNumがある場合)
                if map_num > 0:
                    target_id = map_num - 1 # MapNum is 1-based, internal ID is 0-based
                    processed_ids.add(target_id)
                    
                    if target_id in mw.data.atoms:
                        # === UPDATE (既存原子の変更) ===
                        atom_data = mw.data.atoms[target_id]
                        item = atom_data.get('item')
                        
                        # 属性更新
                        atom_data['symbol'] = atom.GetSymbol()
                        atom_data['charge'] = atom.GetFormalCharge()
                        # 必要なら他の属性も...
                        
                        # アイテムの再描画・移動
                        if item:
                            # 座標移動
                            item.setPos(scene_x, scene_y)
                            # シンボル更新
                            if hasattr(item, 'set_symbol'): 
                                item.set_symbol(atom.GetSymbol())
                            item.update()
                    else:
                        # ID指定があるが、キャンバスにそのIDがない -> 新規作成（ID指定）
                        actual_id = mw.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=atom.GetFormalCharge())
                        
                        # --- STRICT ID REMAPPING (FIX) ---
                        if target_id >= 0 and actual_id != target_id:
                             if target_id not in mw.data.atoms:
                                  try:
                                       adata = mw.data.atoms.pop(actual_id)
                                       mw.data.atoms[target_id] = adata
                                       # Remap internal item ID if possible
                                       if 'item' in adata and hasattr(adata['item'], 'atom_id'):
                                           adata['item'].atom_id = target_id
                                       # print(f"Remapped New Atom {actual_id} -> {target_id}")
                                  except:
                                       pass
                             else:
                                  # Should not happen (checked in else block above)
                                  pass
                        
                        # Update next_atom_id
                        if hasattr(mw.data, 'next_atom_id'):
                            current_max = max(mw.data.atoms.keys()) if mw.data.atoms else 0
                            mw.data.next_atom_id = current_max + 1
                
                else:
                    # B) 新規原子 (MapNumなし、または新しい原子)
                    mw.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=atom.GetFormalCharge())
                    # 新規IDは自動採番

            # --- 3. 削除 (SMILESに含まれなかった原子を消す) ---
            # 「反応で消えた原子」を処理
            existing_ids = list(mw.data.atoms.keys())
            for aid in existing_ids:
                # 今回の更新対象（processed_ids）に含まれず、
                if aid not in processed_ids:
                    # 削除処理
                    if aid in mw.data.atoms:
                        item = mw.data.atoms[aid].get('item')
                        if item:
                            mw.scene.delete_items([item])

            # --- 4. 結合の再構築 (Differential Update) ---
            # UX IMPROVEMENT: Preserve selection by updating in-place instead of recreate
            
            existing_bond_keys = set(mw.data.bonds.keys())
            current_step_bond_keys = set()
            
            for bond in mol.GetBonds():
                b_idx = bond.GetBeginAtomIdx()
                e_idx = bond.GetEndAtomIdx()
                
                # RDKit Index -> MapNum -> Internal ID
                a1_map = mol.GetAtomWithIdx(b_idx).GetAtomMapNum()
                a2_map = mol.GetAtomWithIdx(e_idx).GetAtomMapNum()
                
                if a1_map > 0 and a2_map > 0:
                    id1 = a1_map - 1
                    id2 = a2_map - 1
                    
                    if id1 in mw.data.atoms and id2 in mw.data.atoms:
                        # Bond Key Normalization (Smaller ID first) - Assuming Main Window uses this?
                        # Actually main_window_edit_actions usually uses whatever order came in?
                        # Let's check how create_bond stores it. 
                        # Usually it stores as (id1, id2). But checking keys:
                        # We must match the key format. If undirected, check both?
                        # For robustness, we check both directions or assume canonical.
                        # Let's check if (id1, id2) OR (id2, id1) is in bonds.
                        
                        found_key = None
                        if (id1, id2) in mw.data.bonds: found_key = (id1, id2)
                        elif (id2, id1) in mw.data.bonds: found_key = (id2, id1)
                        
                        # Target Key for CREATE (Canonical logic preferred?)
                        # For consistency let's use what create_bond returns.
                        # We'll just define target_key = (id1, id2) for creating if missing.
                        
                        order = int(bond.GetBondTypeAsDouble())
                        
                        # Stereo Logic
                        b_dir = bond.GetBondDir()
                        stereo = 0
                        if b_dir == Chem.BondDir.BEGINWEDGE: stereo = 1
                        elif b_dir == Chem.BondDir.BEGINDASH: stereo = 2
                        
                        # Double Bond Stereo
                        if order == 2:
                             if bond.GetStereo() == Chem.BondStereo.STEREOZ: stereo = 3
                             elif bond.GetStereo() == Chem.BondStereo.STEREOE: stereo = 4

                        if found_key:
                            # === UPDATE ===
                            current_step_bond_keys.add(found_key)
                            b_data = mw.data.bonds[found_key]
                            
                            # Only update if changed (Minimize drawing)
                            if b_data.get('order') != order or b_data.get('stereo', 0) != stereo:
                                b_data['order'] = order
                                b_data['stereo'] = stereo
                                # Force redraw
                                item = b_data.get('item')
                                if item:
                                     # Assuming create_bond updates if called again? 
                                     # Or strictly manually update?
                                     # Safer to delete/recreate ONLY this modified bond if item API is unknown
                                     # But user wants preservation.
                                     # Let's use create_bond to overwrite properties (if supported)
                                     # or manual item update. Assuming update_bond_visuals exists?
                                     # If simple, create_bond usually handles "update" if data references match?
                                     # No, create_bond makes NEW item usually.
                                     
                                     # Update item properties manually
                                     item1 = mw.data.atoms[id1]['item']
                                     item2 = mw.data.atoms[id2]['item']
                                     # Re-call create_bond helps?
                                     # If we want to strictly preserve selection, we should check setBondOrder on item.
                                     if hasattr(item, 'set_order'): item.set_order(order)
                                     if hasattr(item, 'set_stereo'): item.set_stereo(stereo)
                                     item.update()

                        else:
                            # === CREATE ===
                            item1 = mw.data.atoms[id1]['item']
                            item2 = mw.data.atoms[id2]['item']
                            # Note: create_bond returns key usually?
                            new_key = mw.scene.create_bond(item1, item2, bond_order=order, bond_stereo=stereo)
                            if new_key:
                                 # Depending on implementation, new_key is key or None?
                                 # We just assume it updated mw.data.bonds internally.
                                 # We need to find what key was added.
                                 if (id1, id2) in mw.data.bonds: current_step_bond_keys.add((id1, id2))
                                 elif (id2, id1) in mw.data.bonds: current_step_bond_keys.add((id2, id1))
            
            # Ensure next_atom_id is correct after all additions/remappings
            if mw.data.atoms:
                 mw.data.next_atom_id = max(mw.data.atoms.keys()) + 1

            # === DELETE (Bonds not present in new) ===
            for b_key in existing_bond_keys:
                if b_key not in current_step_bond_keys:
                    if b_key in mw.data.bonds:
                        item = mw.data.bonds[b_key].get('item')
                        if item:
                            mw.scene.delete_items([item])

            # --- Finalize ---
            mw.has_unsaved_changes = True
            if hasattr(mw, 'update_realtime_info'):
                mw.update_realtime_info()
            mw.update_undo_redo_actions()
            
            mw.scene.update()

            self.append_message("System", "Structure updated in-place.", "green")

        except Exception as e:
            self.append_message("System", f"Update Error: {e}", "red")
            import traceback
            traceback.print_exc()

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
