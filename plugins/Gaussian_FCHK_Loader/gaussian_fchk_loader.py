import os
import sys
import importlib
import importlib.util
import fnmatch
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QPushButton, QLabel, 
                             QMessageBox, QHBoxLayout, QWidget, QDockWidget)
from PyQt6.QtCore import Qt, QTimer

# --- Plugin Metadata ---
PLUGIN_NAME = "Gaussian FCHK Loader"
PLUGIN_VERSION = "2026.01.23"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Dispatches FCHK, FCH, and FCK files to appropriate analyzers (Freq vs MO) with priority handling."

def find_file_recursive(root_dir, filename_pattern):
    """
    Recursively find a file matching the pattern in root_dir.
    Returns absolute path of the first match or None.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if fnmatch.fnmatch(filename, filename_pattern):
                return os.path.join(dirpath, filename)
    return None

def find_mo_analyzer_module(root_dir):
    """
    Find the Gaussian MO Analyzer package or file.
    Strategy:
    1. Look for 'gaussian_fchk_mo_analyzer' folder containing '__init__.py'
    """
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if "gaussian_fchk_mo_analyzer" in dirnames:
            # Check if it has __init__.py
            pkg_path = os.path.join(dirpath, "gaussian_fchk_mo_analyzer")
            if os.path.exists(os.path.join(pkg_path, "__init__.py")):
                return pkg_path
    return None

def load_module_from_path(name, path):
    try:
        spec = importlib.util.spec_from_file_location(name, path)
        if spec and spec.loader:
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            return module
    except Exception as e:
        print(f"Failed to load module {name} from {path}: {e}")
    return None

class FCHKLoaderDialog(QDialog):
    def __init__(self, parent, context, fchk_path):
        super().__init__(parent)
        self.context = context
        self.fchk_path = fchk_path
        self.mw = context.get_main_window()
        self.setWindowTitle("Select FCHK Analysis Mode")
        self.resize(300, 150)
        
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        lbl = QLabel(f"Opening: {os.path.basename(self.fchk_path)}\nSelect Analysis Mode:")
        lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(lbl)
        
        # Determine available plugins
        self.freq_analyzer_path = self.find_freq_analyzer()
        self.mo_analyzer_pkg = self.find_mo_analyzer()
        
        # Buttons
        self.btn_freq = QPushButton("Frequency Analysis (IR/Vib)")
        self.btn_freq.clicked.connect(self.launch_freq)
        self.btn_freq.setEnabled(self.freq_analyzer_path is not None)
        if not self.btn_freq.isEnabled():
            self.btn_freq.setToolTip("Plugin not found")
        layout.addWidget(self.btn_freq)
        
        self.btn_mo = QPushButton("MO Analysis (Orbitals)")
        self.btn_mo.clicked.connect(self.launch_mo)
        self.btn_mo.setEnabled(self.mo_analyzer_pkg is not None)
        if not self.btn_mo.isEnabled():
            self.btn_mo.setToolTip("Plugin not found")
        layout.addWidget(self.btn_mo)
        
 

    def find_freq_analyzer(self):
        start_dir = os.path.dirname(os.path.abspath(__file__))
        root = start_dir
        if os.path.basename(os.path.dirname(start_dir)).lower() == 'plugins':
             root = os.path.dirname(start_dir)
        elif os.path.basename(start_dir).lower() == 'plugins':
             root = start_dir
        else:
             root = os.path.dirname(start_dir)

        # Pattern for Freq Analyzer (handling the space version too)
        # "gaussian_fchk_freq_analyzer*.py"
        path = find_file_recursive(root, "gaussian_fchk_freq_analyzer.py")
        return path

    def find_mo_analyzer(self):
        start_dir = os.path.dirname(os.path.abspath(__file__))
        root = os.path.dirname(start_dir)
        return find_mo_analyzer_module(root)

    def close_existing_analyzers(self):
        # Close Frequency Analyzers (DockWidgets)
        for dock in self.mw.findChildren(QDockWidget):
            if dock.windowTitle() == "Gaussian Freq Analyzer":
                dock.close()
                dock.deleteLater()
        
        # Close MO Analyzers (DockWidgets)
        for dock in self.mw.findChildren(QDockWidget):
            if dock.windowTitle() == "Gaussian MO Analyzer":
                dock.close()
                dock.deleteLater()

    def launch_freq(self):
        if not self.freq_analyzer_path: return
        self.close_existing_analyzers()
        try:
            mod = load_module_from_path("gaussian_fchk_freq_analyzer_dynamic", self.freq_analyzer_path)
            if mod and hasattr(mod, 'GaussianFCHKFreqAnalyzer'):
                 # Create Dock Wrapper
                 dock = QDockWidget("Gaussian Freq Analyzer", self.mw)
                 dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
                 dock.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

                 # Instantiate Analyzer (passing dock for close control)
                 self.analyzer = mod.GaussianFCHKFreqAnalyzer(self.mw, dock_widget=dock)
                 
                 dock.setWidget(self.analyzer)
                 self.mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
                 
                 self.analyzer.load_file(self.fchk_path)
                 dock.show()
                 self.accept()
        except Exception as e:
            QMessageBox.critical(self.mw, "Error", f"Failed to launch Frequency Analyzer:\n{e}")

    def launch_mo(self):
        if not self.mo_analyzer_pkg: return
        self.close_existing_analyzers()
        try:
            # Add plugin parent dir to sys.path to support package import
            pkg_parent = os.path.dirname(self.mo_analyzer_pkg)
            if pkg_parent not in sys.path:
                sys.path.insert(0, pkg_parent)

            # Import as proper package
            mod = importlib.import_module("gaussian_fchk_mo_analyzer.gui")
            importlib.reload(mod)

            if mod and hasattr(mod, 'OrbitalWidget'):
                # Create Dock Wrapper
                dock = QDockWidget("Gaussian MO Analyzer", self.mw)
                dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
                dock.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

                # Instantiate Widget
                self.mo_widget = mod.OrbitalWidget(self.mw, self.context, self.fchk_path)
                
                dock.setWidget(self.mo_widget)
                self.mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
                dock.show()
                self.accept()
        except Exception as e:
             QMessageBox.critical(self.mw, "Error", f"Failed to launch MO Analyzer:\n{e}")

def initialize(context):
    """
    Register the FCHK file opener with HIGH PRIORITY (100).
    """
    def open_fchk(path):
        mw = context.get_main_window()
        
        def run_dialog():
            # Pre-check logic for auto-launch
            dlg = FCHKLoaderDialog(mw, context, path)
            
            can_freq = dlg.btn_freq.isEnabled()
            can_mo = dlg.btn_mo.isEnabled()
            
            if can_freq and not can_mo:
                dlg.launch_freq()
            elif can_mo and not can_freq:
                dlg.launch_mo()
            elif not can_freq and not can_mo:
                QMessageBox.warning(mw, "FCHK Loader", "No compatible analysis plugins found (Freq or MO Analyzer).")
            else:
                # If both available, show dialog
                dlg.exec()

        # Run after main window is ready
        QTimer.singleShot(0, run_dialog)
    
    # Register for file open (Import menu / CLI)
    context.register_file_opener(".fchk", open_fchk, priority=100)
    context.register_file_opener(".fck", open_fchk, priority=100)
    context.register_file_opener(".fch", open_fchk, priority=100)

    # Register for Drag and Drop
    def handle_drop(path):
        if path.lower().endswith((".fchk", ".fck", ".fch")):
            open_fchk(path)
            return True
        return False
        
    context.register_drop_handler(handle_drop, priority=100)
