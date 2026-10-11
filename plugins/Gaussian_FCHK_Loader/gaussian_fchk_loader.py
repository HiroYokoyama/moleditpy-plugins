import os
import logging
import sys
import importlib
import importlib.util
import fnmatch
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QPushButton,
    QLabel,
    QMessageBox,
    QDockWidget,
    QFileDialog,
)
from PyQt6.QtCore import Qt, QTimer

# --- Plugin Metadata ---
PLUGIN_NAME = "Gaussian FCHK Loader"
PLUGIN_VERSION = "2026.10.13"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Smart loader for FCHK, FCH, and FCK files. Dispatches to Frequency or MO Analyzer with priority handling."
PLUGIN_CONTEXT = None


def find_file_recursive(root_dir, filename_pattern):
    """
    Recursively find a file matching the pattern in root_dir.
    Returns absolute path of the first match or None.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = sorted(name for name in dirnames if name not in {"_old", ".git", "__pycache__", ".venv"})
        for filename in sorted(filenames):
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
        dirnames[:] = sorted(name for name in dirnames if name not in {"_old", ".git", "__pycache__", ".venv"})
        if "gaussian_fchk_mo_analyzer" in dirnames:
            # Check if it has __init__.py
            pkg_path = os.path.join(dirpath, "gaussian_fchk_mo_analyzer")
            if os.path.exists(os.path.join(pkg_path, "__init__.py")):
                return pkg_path
    return None


def load_module_from_path(name, path):
    """Load the selected file/package without changing the interpreter search path."""
    previous = {key: value for key, value in sys.modules.items() if key == name or key.startswith(name + ".")}
    cached = previous.get(name)
    if cached is not None and getattr(cached, "__file__", None) == path:
        return cached
    try:
        for key in previous:
            del sys.modules[key]
        spec = importlib.util.spec_from_file_location(name, path)
        if spec is None or spec.loader is None:
            raise ImportError(f"No loader for {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module
    except Exception:
        for key in list(sys.modules):
            if key == name or key.startswith(name + "."):
                del sys.modules[key]
        sys.modules.update(previous)
        logging.exception("Failed to load module %s from %s", name, path)
        return None


class AnalyzerDock(QDockWidget):
    cleanup_callback = None

    def closeEvent(self, event):
        widget = self.widget()
        if widget is not None:
            worker = getattr(widget, "worker", None)
            if worker is not None and worker.isRunning():
                worker._is_cancelled = True
                event.ignore()
                return
            if hasattr(widget, "close_plugin"):
                widget.close_plugin()
            widget.close()
        if self.cleanup_callback is not None:
            try:
                self.cleanup_callback()
            except (RuntimeError, AttributeError, KeyError, TypeError, ValueError):
                logging.exception("Could not clean up analyzer surfaces")
        super().closeEvent(event)


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

        lbl = QLabel(
            f"Opening: {os.path.basename(self.fchk_path)}\nSelect Analysis Mode:"
        )
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
        if os.path.basename(os.path.dirname(start_dir)).lower() == "plugins":
            root = os.path.dirname(start_dir)
        elif os.path.basename(start_dir).lower() == "plugins":
            root = start_dir
        else:
            root = os.path.dirname(start_dir)

        # Pattern for Freq Analyzer (handling the space version too)
        # "gaussian_fchk_freq_analyzer*.py"
        path = find_file_recursive(root, "gaussian_fchk_freq_analyzer.py")
        return path

    def find_mo_analyzer(self):
        start_dir = os.path.dirname(os.path.abspath(__file__))
        root = start_dir if os.path.basename(start_dir).lower() == "plugins" else os.path.dirname(start_dir)
        return find_mo_analyzer_module(root)

    def close_existing_analyzers(self):
        docks = [dock for dock in self.mw.findChildren(QDockWidget)
                 if dock.windowTitle() in {"Gaussian Freq Analyzer", "Gaussian MO Analyzer", "ORCA Output Freq Analyzer"}]
        # Check all running calculations before closing any other analyzer.
        for dock in docks:
            worker = getattr(dock.widget(), "worker", None)
            if worker is not None and worker.isRunning():
                worker._is_cancelled = True
                raise RuntimeError("The previous calculation is stopping. Retry when it has finished.")
        for dock in docks:
            widget = dock.widget()
            if widget is not None and hasattr(widget, "close_plugin"):
                widget.close_plugin()
            if not dock.close():
                raise RuntimeError("The previous analyzer could not close.")
            self.mw.removeDockWidget(dock)
            dock.deleteLater()

    def _remove_mo_surfaces(self):
        plotter = self.context.plotter
        if plotter is not None:
            plotter.remove_actor("fchk_iso_p")
            plotter.remove_actor("fchk_iso_n")
            plotter.render()

    def launch_freq(self):
        if not self.freq_analyzer_path:
            return
        self.btn_freq.setEnabled(False)
        self.btn_mo.setEnabled(False)
        dock = None
        try:
            mod = load_module_from_path("gaussian_fchk_freq_analyzer_dynamic", self.freq_analyzer_path)
            if mod is None or not hasattr(mod, "GaussianFCHKFreqAnalyzer"):
                raise ImportError("Frequency Analyzer could not be loaded.")
            # Validate the input before closing the currently displayed analyzer.
            parser = mod.FCHKParser()
            parser.parse(self.fchk_path)
            if not parser.atoms:
                raise ValueError("The FCHK contains no molecule for frequency analysis.")
            self.close_existing_analyzers()
            dock = AnalyzerDock("Gaussian Freq Analyzer", self.mw)
            dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
            dock.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
            self.analyzer = mod.GaussianFCHKFreqAnalyzer(self.context, dock_widget=dock)
            dock.setWidget(self.analyzer)
            dock.visibilityChanged.connect(self.analyzer.on_dock_visibility_changed)
            self.analyzer.load_file(self.fchk_path)
            self.mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
            dock.show()
            self.accept()
        except Exception as exc:
            if dock is not None:
                dock.close()
                dock.deleteLater()
            QMessageBox.critical(self.mw, "Error", f"Failed to launch Frequency Analyzer:\n{exc}")
        finally:
            self.btn_freq.setEnabled(self.freq_analyzer_path is not None)
            self.btn_mo.setEnabled(self.mo_analyzer_pkg is not None)

    def launch_mo(self):
        if not self.mo_analyzer_pkg:
            return
        self.btn_freq.setEnabled(False)
        self.btn_mo.setEnabled(False)
        dock = None
        try:
            mod = load_module_from_path("_fchk_loader_mo_analyzer", os.path.join(self.mo_analyzer_pkg, "__init__.py"))
            if mod is None or not hasattr(mod, "OrbitalWidget"):
                raise ImportError("MO Analyzer could not be loaded.")
            reader = mod.gui.FCHKReader(self.fchk_path)
            mod.gui.BasisSetEngine(reader)
            if not reader.get_xyz_block():
                raise ValueError("The FCHK contains no molecule for MO analysis.")
            self.close_existing_analyzers()
            dock = AnalyzerDock("Gaussian MO Analyzer", self.mw)
            dock.cleanup_callback = self._remove_mo_surfaces
            dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
            dock.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
            self.mo_widget = mod.OrbitalWidget(self.mw, self.context, self.fchk_path)
            dock.setWidget(self.mo_widget)
            self.mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
            dock.show()
            self.accept()
        except Exception as exc:
            if dock is not None:
                dock.close()
                dock.deleteLater()
            QMessageBox.critical(self.mw, "Error", f"Failed to launch MO Analyzer:\n{exc}")
        finally:
            self.btn_freq.setEnabled(self.freq_analyzer_path is not None)
            self.btn_mo.setEnabled(self.mo_analyzer_pkg is not None)


def initialize(context):
    """
    Register the FCHK file opener with HIGH PRIORITY (100).
    """
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context

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
                QMessageBox.warning(
                    mw,
                    "FCHK Loader",
                    "No compatible analysis plugins found (Freq or MO Analyzer).",
                )
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


def run(mw):
    path, _ = QFileDialog.getOpenFileName(
        mw, "Open Gaussian FCHK", "", "Gaussian FCHK (*.fchk *.fck);;All Files (*)"
    )
    if not path:
        return

    context = PLUGIN_CONTEXT
    if not context:
        return

    def run_dialog():
        dlg = FCHKLoaderDialog(mw, context, path)
        can_freq = dlg.btn_freq.isEnabled()
        can_mo = dlg.btn_mo.isEnabled()
        if can_freq and not can_mo:
            dlg.launch_freq()
        elif can_mo and not can_freq:
            dlg.launch_mo()
        elif not can_freq and not can_mo:
            QMessageBox.warning(
                mw,
                "FCHK Loader",
                "No compatible analysis plugins found (Freq or MO Analyzer).",
            )
        else:
            dlg.exec()

    QTimer.singleShot(0, run_dialog)
