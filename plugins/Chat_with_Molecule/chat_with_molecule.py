#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Chat with Molecule Plugin
Allows chatting with Google Gemini API about the currently loaded molecule.
"""

import sys
import os
import json
import io
import base64
import re
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QTextEdit, QPushButton, QLineEdit, QMessageBox, 
    QWidget, QFrame, QSizePolicy, QComboBox, QFileDialog, QTextBrowser
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QUrl, QTimer
from PyQt6.QtGui import QColor, QTextCursor, QDesktopServices

try:
    from rdkit import Chem
except ImportError:
    Chem = None

# --- Matplotlib for LaTeX Rendering ---
try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# --- Metadata ---
PLUGIN_NAME = "Chat with Molecule"
PLUGIN_VERSION = "2025.12.27"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Chat with Google Gemini about the current molecule. Automatically injects SMILES context."
PLUGIN_ID = "chat_with_molecule"

SYSTEM_PROMPT = """You are an expert computational and organic chemistry assistant embedded within the advanced molecular editor software "MoleditPy". 
Your users are researchers, students, or chemistry enthusiasts. Adhere to the following guidelines:

### 1. Protocol for Molecular Identification (CRITICAL)
When provided with a SMILES string as `Context`, you MUST NOT rely on intuition. You must perform the following internal reasoning steps before generating a response:
1.  **Stoichiometry Check**: Mentally count the atoms (e.g., Carbon count, Hydrogen count) to verify the molecular formula (e.g., C16H10 vs C24H12).
2.  **Topology Analysis**: Analyze the number of rings and their fusion pattern.
3.  **Identification**: Based on the above data, identify the molecule. If uncertain, state "It resembles [Molecule Name] but..."

### 2. Response Style
* **Tone**: Professional, intellectual, and helpful.
* **Formatting**:
    * Use **SIMPLE** LaTeX format for chemical formulas and math (e.g., $C_{16}H_{10}$).
    * **Avoid** complex commands like `\\xrightarrow` or `\\stackrel` as the renderer is limited. Use `\\rightarrow` or `->` instead.
    * **CRITICAL**: Always format ANY molecular structure you mention *in plain text* as a clickable link like this: `[Name](smiles:SMILES_STRING)`.
    * **IMPORTANT**: Do **NOT** put SMILES strings inside LaTeX equations (e.g. inside `$...$`). Use chemical formulas (e.g. $C_6H_6$) or names instead.

### 3. Constraints
* Do not state uncertain information as absolute fact. Admit ambiguity if the SMILES is complex.
* If the SMILES appears chemically invalid, politely inform the user.

### 4. MoleditPy Software Context
* **Role:** You are the integrated AI for **MoleditPy**.
* **Native Features:** Suggest built-in tools for immediate modeling and export:
    * **Modeling:** "3D Optimization (MMFF64/UFF: sufficient quality for viewing)," "Alignment to Axes," "Planarization," and "Constrained Optimization."
    * **Export:** "STL/OBJ Export" (for 3D printing) and "Image Export."
* **Plugin Ecosystem:** For advanced tasks, recommend installing specific extensions from the **[Plugin Explorer](https://hiroyokoyama.github.io/moleditpy-plugins/explorer/)**:
    * **Calculation:** "Gaussian/ORCA Input Generator" (prepare DFT inputs).
    * **Analysis:** "MS Spectrum Simulation" (mass spec) or "Gaussian/ORCA Freq Analyzer" (IR/vibrations).
    * **Visualization:** "Cube File Viewer" (orbitals/densities) or "Mapped Cube Viewer" (ESP maps).
    * **Utility:** "PubChem Name Resolver" (load by name) or "Animated XYZ Giffer."
"""

GENERATION_CONFIG = {
    # 【最重要】創造性を極限まで下げる。
    # 0.0〜0.1にすることで、モデルは確率が最も高い（論理的な）トークンのみを選びます。
    "temperature": 0.1,
    
    # 確度の高い上位95%のトークン群から選択
    "top_p": 0.95,
    
    # 選択候補を上位40個に絞る
    "top_k": 40,
    
    # 回答の最大長（必要に応じて調整）
    "max_output_tokens": 2048,
    
    # 応答を必ず「text/plain」で返させる（JSONモード等は不要）
    "response_mime_type": "text/plain",
}

SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "chat_with_molecule.json")
LOG_FILE = os.path.join(os.path.dirname(__file__), "chat_with_molecule_log.txt")

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
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, 'r') as f:
                return json.load(f)
        except:
            pass
    return {}

def save_settings(settings):
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

def latex_to_html(latex_str):
    """
    Render LaTeX string to base64 encoded PNG image using matplotlib.
    Returns an HTML <img> tag.
    """
    if not HAS_MATPLOTLIB:
        return f"<i>{latex_str}</i>" # Fallback if no matplotlib

    try:
        # Create a pure Figure (no global state)
        # Scaled DPI down to 72 (approx 60% of 120) to reduce image size as requested.
        fig = Figure(figsize=(4, 0.5), dpi=72)
        canvas = FigureCanvasAgg(fig)
        
        # Ensure we don't double wrap
        content = latex_str.strip()
        
        # --- Simple Fixes for Matplotlib Limitations ---
        # Matplotlib's mathtext is not a full TeX engine. 
        # We replace unsupported commands with simpler alternatives to prevent render failure.
        content = content.replace(r'\xrightarrow', r'\rightarrow')
        content = content.replace(r'\text', r'\mathrm') # \mathrm is safer than \text in some mpl versions
        
        # Remove \stackrel arguments roughly (not perfect, but better than crash)
        # Replacing \stackrel{top}{base} with base^{top}
        import re
        content = re.sub(r'\\stackrel\{([^}]+)\}\{([^}]+)\}', r'\2^{\1}', content)
        
        if not content.startswith('$'):
            content = f"${content}$"
            
        # Draw text centered
        text = fig.text(0.5, 0.5, content, fontsize=14, ha='center', va='center')
        
        # Save to buffer
        buf = io.BytesIO()
        # bbox_inches='tight' is crucial to crop the image to the equation size
        fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0.05, transparent=True)
        
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        # Use simple vertical align and max-width to prevent overflow, but respect natural size (now smaller via DPI).
        return f'<img src="data:image/png;base64,{img_base64}" style="vertical-align: middle; max-width: 100%;">'
        
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

class GenAIWorker(QThread):
    """Worker thread to handle API calls to avoid freezing UI"""
    response_received = pyqtSignal(object) # Changed to object to pass full response
    error_occurred = pyqtSignal(str)

    def __init__(self, chat_session, user_message):
        super().__init__()
        self.chat_session = chat_session
        self.user_message = user_message

    def run(self):
        try:
            response = self.chat_session.send_message(self.user_message)
            self.response_received.emit(response)
        except Exception as e:
            self.error_occurred.emit(str(e))

class ChatMoleculeWindow(QDialog):
    def __init__(self, main_window):
        super().__init__(main_window)
        self.main_window = main_window
        self.setWindowTitle("Chat with Molecule (Gemini)")
        self.resize(500, 700)
        self.settings = load_settings()
        self.chat_session = None
        self.worker = None

        self.chat_history_log = [] # List of tuples (sender, text)
        self.last_smiles = None # Track previous state
        self.last_error = None # Track previous error state
        self.pending_context_msg = None # Pending context update to send with next user message
        self.pending_info_text = None # Pending visual notification to show in chat on next send
        self.first_check_done = False # Track if we've done the initial check

        self.init_ui()
        
        # Defer initialization to allow window to show immediately
        QTimer.singleShot(100, self.initialize_session)

        # Polling Timer for Auto-Update
        self.poll_timer = QTimer(self)
        self.poll_timer.timeout.connect(self.check_molecule_change)
        self.poll_timer.start(2000) # Check every 2 seconds

    def check_molecule_change(self):
        """Check if the main window's molecule has changed"""
        # We need to use the method to safely get/clean SMILES
        current_smiles, error = self.get_current_molecule_smiles()
        
        # Determine if this is the first check or a subsequent change
        is_first_check = not self.first_check_done
        # Change detected (SMILES changed OR Error state changed)
        has_changed = (current_smiles != self.last_smiles) or (error != self.last_error)

        if is_first_check or has_changed:
             self.last_smiles = current_smiles
             self.last_error = error
             self.first_check_done = True
             
             # UI Update Logic (Label) - Update on first check AND changes
             if current_smiles:
                 self.lbl_context.setText(f"Context: {current_smiles}")
             else:
                 context_label_text = "Context: No valid molecule found."
                 if error:
                      context_label_text += f" ({error})"
                 self.lbl_context.setText(context_label_text)

             # Chat Log Logic - Only on ACTUAL changes (skip first check to avoid spam)
             if not is_first_check:
                 if current_smiles:
                     # Store pending message instead of sending immediately
                     self.pending_context_msg = (
                     f"System Update: The user has switched to a new molecule with SMILES: {current_smiles}. "
                     f"Please use this as the new context. "
                     f"Do not reply to this update; simply answer the user's message below."
                 )
                     
                     # Update Chat locally (Silent now, queued for visual)
                     append_log("INFO", f"Context auto-updated (queued): {current_smiles}")
                     self.pending_info_text = f"Context updated: {current_smiles}"
                     
                 else:
                     self.pending_context_msg = "System Update: The user has unloaded the molecule."
                     if error:
                         append_log("INFO", f"Molecule unloaded/scan failed: {error}")
                     else:
                         append_log("INFO", "Molecule unloaded.")
                     self.pending_info_text = "Molecule unloaded."

    def init_ui(self):
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
        default_model = self.settings.get("model", "gemini-1.5-flash")
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
        lbl_warning = QLabel("<b style='color: orange;'>⚠️ Do not include confidential information for free API.</b>")
        lbl_warning.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(lbl_warning)

        # --- Chat Display ---
        self.chat_display = ChatBrowser(self) # Use custom subclass
        self.chat_display.anchorClicked.connect(self.handle_link)
        self.chat_display.setReadOnly(True)
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

    def save_settings(self):
        key = self.txt_api_key.text().strip()
        model = self.combo_model.currentText().strip()

        if not key:
            QMessageBox.warning(self, "Warning", "API Key cannot be empty.")
            return

        self.settings["api_key"] = key
        self.settings["model"] = model
        save_settings(self.settings)
        QMessageBox.information(self, "Saved", "Settings saved successfully.\nRe-initializing session...")
        self.initialize_session()

    def fetch_models(self):
        key = self.txt_api_key.text().strip()
        if not key:
             QMessageBox.warning(self, "Error", "Enter API Key first.")
             return

        if not HAS_GENAI: return

        try:
            genai.configure(api_key=key)
            self.combo_model.clear()

            self.append_message("System", "Fetching available models...", "blue")
            models = [m.name for m in genai.list_models() if 'generateContent' in m.supported_generation_methods]

            if models:
                self.combo_model.addItems(models)
                self.append_message("System", f"Found {len(models)} models.", "green")
                # Try to select the previously saved one if in list
                saved = self.settings.get("model", "")
                index = self.combo_model.findText(saved)
                if index >= 0:
                    self.combo_model.setCurrentIndex(index)
                else:
                    # Select first valid gemini if possible
                    geminis = [i for i, m in enumerate(models) if 'gemini' in m]
                    if geminis: self.combo_model.setCurrentText(geminis[0])
            else:
                self.append_message("System", "No models found with 'generateContent' capability.", "orange")

        except Exception as e:
             QMessageBox.warning(self, "Error", f"Failed to list models:\n{e}")

    def append_message(self, sender, text, color="black"):
        self.chat_history_log.append({"sender": sender, "text": text})

        # Auto-log to file
        append_log(sender, text)

        # --- LaTeX / Chemical Formula Formatting ---
        # 1. Block Math: $$...$$
        # We process this BEFORE markdown to avoid markdown messing up LaTeX syntax.
        
        processed_text = text
        
        def replace_latex_block(match):
            return latex_to_html(match.group(1))
            
        def replace_latex_inline(match):
            return latex_to_html(match.group(1))

        # Using regex to find $$...$$ (Block)
        # We replace it with a unique placeholder or directly with the HTML
        # If we do it directly, Markdown might ignore it if it's inside HTML tags.
        # But standard markdown (Python lib) typically ignores content inside HTML block tags or handles them safely?
        # Let's try direct replacement.
        processed_text = re.sub(r'\$\$([^$]+)\$\$', replace_latex_block, processed_text)
        
        # 2. Inline Math: $...$
        # Be careful not to match currency. Regex usually requires lookarounds or spacing checks,
        # but for this context (scientific), $...$ is almost always math.
        # We match $...$ but NOT $$...$$ (already handled).
        processed_text = re.sub(r'(?<!\$)\$([^$]+)\$(?!\$)', replace_latex_inline, processed_text)

        # Render Markdown (now that LaTeX is turned into HTML images)
        if HAS_MARKDOWN and sender in ["Gemini", "You", "System"]:
            try:
                # Enable useful extensions. 'tables' and 'fenced_code' are standard.
                # We need to ensure we don't break the HTML images we just inserted.
                # Python-Markdown usually preserves HTML if not in safe_mode.
                processed_text = markdown.markdown(processed_text, extensions=['fenced_code', 'tables'])
            except Exception as e:
                print(f"Markdown error: {e}")
                processed_text = processed_text.replace(chr(10), '<br>')
        else:
             processed_text = processed_text.replace(chr(10), '<br>')

        # Enforce SMILES links manually if they didn't get parsed
        # Pattern: [Label](smiles:Data) -> <a href="smiles:Data">Label</a>
        # Note: Markdown might have already converted this to <a> tags if the scheme was recognized,
        # but valid markdown links usually require http/https or known schemes.
        # Python-Markdown might not support 'smiles:' scheme automatically.
        
        # Check if markdown already handled it?
        # If not, we run our regex.
        # Regex for markdown link: \[([^\]]+)\]\((smiles:[^\)]+)\)
        pattern = r"\[([^\]]+)\]\((smiles:[^\)]+)\)"
        processed_text = re.sub(pattern, r'<a href="\2">\1</a>', processed_text)

        # CSS Styling for Markdown elements
        style = """
        <style>
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
            "Markdown Files (*.md);;Text Files (*.txt);;All Files (*)"
        )

        if not file_path:
            return

        try:
            import datetime
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(f"# Chat with Molecule History\n")
                f.write(f"Date: {current_time}\n")
                f.write(f"Plugin Version: {PLUGIN_VERSION}\n\n")

                for entry in self.chat_history_log:
                    sender = entry['sender']
                    text = entry['text']
                    f.write(f"## {sender}\n{text}\n\n")
            QMessageBox.information(self, "Success", f"Exported to {os.path.basename(file_path)}")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to export: {e}")



    def initialize_session(self):
        if not HAS_GENAI:
            self.append_message("System", "Error: 'google-generativeai' library is not installed.", "red")
            self.txt_input.setEnabled(False)
            self.btn_send.setEnabled(False)
            return

        api_key = self.settings.get("api_key")
        if not api_key:
            self.append_message("System", "Please enter your Google API Key above.", "orange")
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
        
        context_msg = ""
        if smiles:
            self.lbl_context.setText(f"Context: {smiles}")
            context_msg = (
                f"I am currently looking at a molecule with this SMILES string: {smiles}. "
            f"Please use this as context for our conversation. "
            f"Please do not reply to this information."
        )
            self.pending_context_msg = context_msg # Store to send with first user message
            self.pending_info_text = f"Context loaded: {smiles}"
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

    def get_current_molecule_smiles(self):
        """Returns (smiles, error_message)"""
        if not Chem:
            return None, "RDKit library not found/loaded."
            
        if not self.main_window:
            return None, "Internal Error: No reference to main application."
            
        # Strategy:
        # 1. Try 'current_mol' (Cached 3D RDKit object).
        # 2. If None, try 'data.to_rdkit_mol()' (Live 2D scene reconstruction).
        
        mol = None
        
        # 1. Try cached 3D mol
        try:
            if hasattr(self.main_window, 'current_mol'):
                mol = self.main_window.current_mol
        except Exception:
            pass # Ignore access errors, proceed to fallback

        # 2. If no cached 3D mol, try reconstructing from 2D Data
        if mol is None:
            if hasattr(self.main_window, 'data') and hasattr(self.main_window.data, 'to_rdkit_mol'):
                # Quick check if there are atoms to convert
                if hasattr(self.main_window.data, 'atoms') and self.main_window.data.atoms:
                    try:
                        mol = self.main_window.data.to_rdkit_mol()
                    except Exception as e:
                        return None, f"Failed to reconstruct molecule from scene: {e}"
        
        # Final check
        if mol is None:
             # Distinguish between "Empty" and "Broken wrapper"
             if hasattr(self.main_window, 'data') and hasattr(self.main_window.data, 'atoms') and not self.main_window.data.atoms:
                 return None, "No molecule loaded (Scene is empty)."
             return None, "No valid molecule found (current_mol is None)."
            
        try:
            # Attempt 1: Clean SMILES (No Hydrogens)
            try:
                mol_no_h = Chem.RemoveHs(mol)
                return Chem.MolToSmiles(mol_no_h), None
            except:
                pass

            # Attempt 2: Direct SMILES
            return Chem.MolToSmiles(mol), None
        except Exception as e:
            return None, f"SMILES conversion failed: {str(e)}"

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
        if self.pending_context_msg:
            # [CONTEXT INJECTION POINT]
            # This is where the SMILES info and other system warnings are attached to the 
            # user's message before being sent to Gemini.
            full_text_to_send = f"{self.pending_context_msg}\n\n{text}"
            self.pending_context_msg = None # Clear pending

        # Worker
        self.worker = GenAIWorker(self.chat_session, full_text_to_send)
        self.worker.response_received.connect(self.on_response)
        self.worker.error_occurred.connect(self.on_error)
        self.worker.finished.connect(self.on_worker_finished)
        self.worker.start()
        return None

    def on_response(self, response):
        text = response.text
        self.append_message("Gemini", text, "black")
        self.log_usage(response)

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
        except Exception as e:
            print(f"Failed to log usage: {e}")

    def on_error(self, error_msg):
        self.append_message("Error", error_msg, "red")

    def on_worker_finished(self):
        self.txt_input.setEnabled(True)
        self.btn_send.setEnabled(True)
        self.txt_input.setFocus()

def run(main_window):
    """
    Entry point for the plugin.
    Display the ChatMoleculeWindow.
    """
    # Modeless implementation
    # We store the instance on the main_window to prevent garbage collection
    # and to allow checking if it's already open.

    if not hasattr(main_window, 'chat_molecule_window_instance'):
        main_window.chat_molecule_window_instance = None

    # If instance exists but is hidden/closed, we can reuse or recreate.
    # Recreating is often safer to ensure clean state as user requested "fresh session on launch".
    # But if it's just hidden or backgrounded, maybe raise it.
    # Let's recreate if valid instance doesn't exist or isn't visible?
    # Simply: If open, raise. If not, create new.

    if main_window.chat_molecule_window_instance:
        try:
            # Force close the old instance to ensure code updates (like ChatBrowser class) are applied.
            # If we just raise it, it keeps using the old class definition from before the reload.
            main_window.chat_molecule_window_instance.close()
            main_window.chat_molecule_window_instance = None
        except:
            pass

    # Create new instance
    dialog = ChatMoleculeWindow(main_window)
    main_window.chat_molecule_window_instance = dialog

    # Use show() for modeless
    dialog.show()
