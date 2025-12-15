import sys
import io
import code
import traceback
from contextlib import redirect_stdout, redirect_stderr
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QTextEdit, QLineEdit, 
    QPushButton, QLabel, QWidget
)
from PyQt6.QtGui import QFont, QColor
from PyQt6.QtCore import Qt
import rdkit.Chem as Chem

PLUGIN_NAME = "Python Console"

class HistoryLineEdit(QLineEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.history = []
        self.history_index = 0
    
    def append_history(self, text):
        if text and (not self.history or self.history[-1] != text):
            self.history.append(text)
        self.history_index = len(self.history)
    
    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Up:
            if self.history_index > 0:
                self.history_index -= 1
                self.setText(self.history[self.history_index])
        elif event.key() == Qt.Key.Key_Down:
            if self.history_index < len(self.history) - 1:
                self.history_index += 1
                self.setText(self.history[self.history_index])
            else:
                self.history_index = len(self.history)
                self.clear()
        else:
            super().keyPressEvent(event)

class PythonConsoleDialog(QDialog):
    def __init__(self, main_window):
        super().__init__(main_window)
        self.main_window = main_window
        self.setWindowTitle("MoleditPy Python Console")
        self.resize(600, 400)
        
        # UI Setup
        layout = QVBoxLayout()
        
        # Output Area (Log)
        self.output_area = QTextEdit()
        self.output_area.setReadOnly(True)
        self.output_area.setStyleSheet("background-color: #1e1e1e; color: #dcdcdc;")
        self.output_area.setFont(QFont("Consolas", 10))
        layout.addWidget(self.output_area)
        
        # Input Area with History
        self.input_line = HistoryLineEdit()
        self.input_line.setPlaceholderText("Enter Python code...")
        self.input_line.setStyleSheet("background-color: #2d2d2d; color: #ffffff; border: 1px solid #3e3e3e;")
        self.input_line.setFont(QFont("Consolas", 10))
        self.input_line.returnPressed.connect(self.run_code)
        layout.addWidget(self.input_line)
        
        # Help Label
        help_text = QLabel("Available vars: 'mw' (MainWindow), 'mol' (current_mol), 'Chem' (rdkit.Chem)")
        help_text.setStyleSheet("color: gray; font-size: 10px;")
        layout.addWidget(help_text)

        self.setLayout(layout)
        
        # Initialize execution environment (namespace)
        self.local_scope = {
            'mw': self.main_window,
            'Chem': Chem,
            'mol': self._get_best_mol(),
        }

        # Initialize Interpreter
        self.interpreter = code.InteractiveInterpreter(self.local_scope)

        self.append_output("MoleditPy Console Ready.")
        self.append_output(">>> Type commands and press Enter.")

    def _get_best_mol(self):
        """Helper to get the most relevant RDKit molecule object.
        """
        mol = getattr(self.main_window, 'current_mol', None)

        return mol

    def append_output(self, text, color=None):
        if color:
            self.output_area.append(f"<span style='color: {color};'>{text}</span>")
        else:
            self.output_area.append(text)

    def run_code(self):
        command = self.input_line.text()
        if not command:
            return

        # Handle History
        self.input_line.append_history(command)
        self.input_line.clear()

        # Display Input
        self.append_output(f">>> {command}", color="#4CAF50")
        
        # Sync variables
        self.local_scope['mol'] = self._get_best_mol()
        
        if self.local_scope['mol'] is None:
             # Optional: warn user if they try to use 'mol' and it's still None
             # but only if 'mol' appears in command to avoid spam
             if 'mol' in command:
                 print("Warning: 'mol' is None (no valid 2D or 3D structure found).")

        # Capture Output
        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        try:
            with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
                # runsource handles compilation and syntax errors.
                # It returns True if more input is needed (incomplete code), which we can handle or just report.
                more = self.interpreter.runsource(command, "<console>", "single")
                
                if more:
                    print("(Incomplete input - multiline not fully supported yet)")
        except Exception:
            # This catches errors OUTSIDE runsource's internal handling if any
            traceback.print_exc(file=stderr_capture)

        # Process Captured Output
        out_str = stdout_capture.getvalue()
        err_str = stderr_capture.getvalue()

        if out_str:
            self.output_area.append(out_str.strip())
        if err_str:
            self.append_output(err_str.strip(), color="#FF5252")

        # Refresh UI if needed
        # If the user modified 'mol', we might want to push it back?
        # For now, read-only assumption for simpler integration, 
        # but if they modify 'mw.data' directly, we might need a refresh.
        pass

# Plugin Entry Point
def run(main_window):
    if not hasattr(main_window, 'python_console_dialog'):
        main_window.python_console_dialog = PythonConsoleDialog(main_window)
    
    main_window.python_console_dialog.show()
    main_window.python_console_dialog.raise_()
    main_window.python_console_dialog.activateWindow()
