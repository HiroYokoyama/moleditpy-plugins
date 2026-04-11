import io
import code
import traceback
import keyword
from contextlib import redirect_stdout, redirect_stderr
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QTextEdit, QPlainTextEdit,
    QLabel, QSplitter
)
from PyQt6.QtGui import (
    QFont, QColor, QSyntaxHighlighter, QTextCharFormat,
    QTextCursor
)
from PyQt6.QtCore import Qt, QRegularExpression, pyqtSignal
import rdkit.Chem as Chem

PLUGIN_VERSION = "2026.04.12"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Embedded Python console for interactive scripting."
PLUGIN_NAME = "Python Console"

class PythonHighlighter(QSyntaxHighlighter):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.highlighting_rules = []

        # Keyword format
        keyword_format = QTextCharFormat()
        keyword_format.setForeground(QColor("#FFA500"))  # Orange
        keyword_format.setFontWeight(QFont.Weight.Bold)
        keywords = keyword.kwlist
        for word in keywords:
            pattern = QRegularExpression(f"\\b{word}\\b")
            self.highlighting_rules.append((pattern, keyword_format))

        # String format ("" or '')
        string_format = QTextCharFormat()
        string_format.setForeground(QColor("#6A8759"))  # Greenish
        self.highlighting_rules.append((QRegularExpression("\".*\""), string_format))
        self.highlighting_rules.append((QRegularExpression("\'.*\'"), string_format))

        # Comment format
        comment_format = QTextCharFormat()
        comment_format.setForeground(QColor("#808080"))  # Grey
        self.highlighting_rules.append((QRegularExpression("#[^\n]*"), comment_format))
        
        # Number format
        number_format = QTextCharFormat()
        number_format.setForeground(QColor("#6897BB")) # Blueish
        self.highlighting_rules.append((QRegularExpression("\\b[0-9]+\\b"), number_format))

        # Builtins / Magic (Optional, simple set)
        builtin_format = QTextCharFormat()
        builtin_format.setForeground(QColor("#8888C6"))
        builtins = ["print", "len", "range", "list", "dict", "set", "str", "int", "float", "help"]
        for word in builtins:
             pattern = QRegularExpression(f"\\b{word}\\b")
             self.highlighting_rules.append((pattern, builtin_format))


    def highlightBlock(self, text):
        for pattern, format in self.highlighting_rules:
            match_iterator = pattern.globalMatch(text)
            while match_iterator.hasNext():
                match = match_iterator.next()
                self.setFormat(match.capturedStart(), match.capturedLength(), format)


class ConsoleInput(QPlainTextEdit):
    execute_signal = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.history = []
        self.history_index = 0
        self.setPlaceholderText("Enter Python code... (Shift+Enter for new line)")
        # Make tab strictly 4 spaces
        self.setTabStopDistance(self.fontMetrics().horizontalAdvance(' ') * 4)

    def append_history(self, text):
        if text and (not self.history or self.history[-1] != text):
            self.history.append(text)
        self.history_index = len(self.history)

    def keyPressEvent(self, event):
        # Enter -> Execute
        # Shift+Enter -> Newline
        if event.key() == Qt.Key.Key_Return:
            modifiers = event.modifiers()
            if modifiers & Qt.KeyboardModifier.ShiftModifier:
                # Insert newline
                super().keyPressEvent(event)
            else:
                # Execute
                self.execute_signal.emit()
                return # Don't insert newline
        
        elif event.key() == Qt.Key.Key_Up:
            # Navigate history only if we are at top (or 1 line) OR if modifiers? 
            # Standard console behavior: Up always goes to history if cursor is at top line?
            # Let's simple it: ONLY if Control+Up, or if text is empty? 
            # Or standard shell behavior: Up goes up through lines, then through history.
            
            # Simple implementation: Up moves cursor unless at top line, then history.
            cursor = self.textCursor()
            block_number = cursor.blockNumber()
            
            if block_number == 0 and self.history and self.history_index > 0:
                 # Check if we are really at the very top (start of doc)?
                 # Actually, usually shells allow Up to access history easily.
                 # Let's say: if single line, Up -> History. If multi-line, Up -> Move up.
                 # If at top line of multiline -> History.
                 
                 self.history_index -= 1
                 self.setPlainText(self.history[self.history_index])
                 self.moveCursor(QTextCursor.MoveOperation.End)
                 return
            
            super().keyPressEvent(event)

        elif event.key() == Qt.Key.Key_Down:
             cursor = self.textCursor()
             block_count = self.blockCount()
             block_number = cursor.blockNumber()

             if block_number == block_count - 1:
                 # At last line
                 if self.history_index < len(self.history) - 1:
                    self.history_index += 1
                    self.setPlainText(self.history[self.history_index])
                    self.moveCursor(QTextCursor.MoveOperation.End)
                    return
                 elif self.history_index == len(self.history) - 1:
                     # Go to empty
                     self.history_index += 1
                     self.clear()
                     return
            
             super().keyPressEvent(event)
        
        else:
            super().keyPressEvent(event)


class PythonConsoleDialog(QDialog):
    def __init__(self, context):
        super().__init__(context.get_main_window())
        self.context = context
        # Register window for V3 lifecycle management
        self.context.register_window("main_panel", self)
        
        self.setWindowTitle("MoleditPy Python Console")
        self.resize(700, 500)
        
        # UI Setup
        layout = QVBoxLayout()
        
        # Splitter to allow resizing output/input area
        splitter = QSplitter(Qt.Orientation.Vertical)
        
        # Output Area (Log)
        self.output_area = QTextEdit()
        self.output_area.setReadOnly(True)
        self.output_area.setStyleSheet("background-color: #1e1e1e; color: #dcdcdc;")
        self.output_font = QFont("Consolas", 10)
        self.output_area.setFont(self.output_font)
        splitter.addWidget(self.output_area)
        
        # Input Area with History
        self.input_area = ConsoleInput()
        self.input_area.setStyleSheet("background-color: #2b2b2b; color: #ffffff; border: 1px solid #3e3e3e;")
        self.input_area.setFont(self.output_font)
        
        # Syntax Highlighter
        self.highlighter = PythonHighlighter(self.input_area.document())
        
        self.input_area.execute_signal.connect(self.run_code)
        splitter.addWidget(self.input_area)
        
        # Set initial sizes (Output 70%, Input 30%)
        splitter.setSizes([350, 150])
        
        layout.addWidget(splitter)
        
        # Help Label
        help_text = QLabel("Avalable vars: 'mw' (MainWindow), 'mol' (current_mol), 'Chem' (rdkit.Chem)\nUse Shift+Enter for new line, Enter to Run.")
        help_text.setStyleSheet("color: gray; font-size: 10px;")
        layout.addWidget(help_text)

        self.setLayout(layout)
        
        # Initialize execution environment (namespace)
        # # [DIRECT ACCESS] to core objects for scripting
        self.local_scope = {
            'mw': self.context.get_main_window(),
            'Chem': Chem,
            'mol': self._get_best_mol(),
        }

        # Initialize Interpreter
        self.interpreter = code.InteractiveInterpreter(self.local_scope)

        self.append_output("MoleditPy Console Ready.")
        self.append_output(">>> Type commands and press Enter.")

    def _get_best_mol(self):
        """Helper to get the most relevant RDKit molecule object."""
        return self.context.current_molecule

    def append_output(self, text, color=None):
        if color:
            self.output_area.append(f"<span style='color: {color};'>{text}</span>")
        else:
            self.output_area.append(text)

    def run_code(self):
        command = self.input_area.toPlainText()
        if not command.strip():
            return

        # Handle History
        self.input_area.append_history(command)
        self.input_area.clear()

        # Display Input
        lines = command.split('\n')
        self.append_output(f">>> {lines[0]}", color="#4CAF50")
        for line in lines[1:]:
             self.append_output(f"... {line}", color="#4CAF50")
        
        # Sync variables
        self.local_scope['mol'] = self._get_best_mol()
        self.local_scope['mw'] = self.context.get_main_window()
        
        if self.local_scope['mol'] is None and 'mol' in command:
             print("Warning: 'mol' is None (no valid molecule found).")

        # Capture Output
        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        try:
            with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
                more = self.interpreter.runsource(command, "<console>", "exec" if '\n' in command else "single")
                if more:
                    self.append_output("... (Incomplete input/block)", color="#FF9800")
        except Exception:
            traceback.print_exc(file=stderr_capture)

        # Process Captured Output
        out_str = stdout_capture.getvalue()
        err_str = stderr_capture.getvalue()

        if out_str:
            self.output_area.append(out_str.strip())
        if err_str:
            self.append_output(err_str.strip(), color="#FF5252")

        # Scroll to bottom
        self.output_area.verticalScrollBar().setValue(self.output_area.verticalScrollBar().maximum())

def initialize(context):
    def run_console():
        win = context.get_window("main_panel")
        if win is None:
            win = PythonConsoleDialog(context)
        
        if win.isVisible():
            win.hide()
        else:
            win.show()
            win.raise_()
            win.activateWindow()
            

def run(mw):
    if hasattr(mw, 'host'):
        mw = mw.host
    from moleditpy.plugins.plugin_interface import PluginContext
    context = PluginContext(mw.plugin_manager, PLUGIN_NAME)
    
    win = context.get_window("main_panel")
    if win is None:
        win = PythonConsoleDialog(context)
    win.show()
    win.raise_()
    win.activateWindow()
