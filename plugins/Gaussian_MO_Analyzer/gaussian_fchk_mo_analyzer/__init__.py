from .gui import OrbitalWidget
from .analyzer import UnsupportedBasisError
from PySide6.QtCore import QTimer

# --- Plugin Metadata ---
PLUGIN_NAME = "Gaussian MO Analyzer"
PLUGIN_VERSION = "2026.10.31"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0, <6.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Visualizes Molecular Orbitals from Gaussian FCHK files by generating Cube "
    "files. Supports S, P, D, F and spherical G shells; files using basis "
    "functions it cannot render exactly are rejected rather than shown wrong."
)
PLUGIN_CONTEXT = None


def _warn_unsupported(parent, exc):
    """Show the refusal message for a basis set we cannot render exactly."""
    from PySide6.QtWidgets import QMessageBox

    QMessageBox.warning(parent, "Unsupported basis set", str(exc))


def _warn_invalid(parent, exc):
    from PySide6.QtWidgets import QMessageBox

    QMessageBox.warning(parent, "Invalid FCHK data", str(exc))


def initialize(context):
    """
    Register the FCHK file opener.
    """
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context

    def open_fchk(path):
        mw = context.get_main_window()
        # Create and show dialog (keep reference to avoid GC)
        # Attach to context to keep alive
        try:
            context._fchk_dialog = OrbitalWidget(mw, context, path)
        except UnsupportedBasisError as exc:
            _warn_unsupported(mw, exc)
            return
        except (OSError, ValueError, TypeError, IndexError, RuntimeError) as exc:
            _warn_invalid(mw, exc)
            return
        context._fchk_dialog.show()

    context.register_file_opener(".fchk", open_fchk, priority=10)
    context.register_file_opener(".fck", open_fchk, priority=10)
    context.register_file_opener(".fch", open_fchk, priority=10)

    def handle_drop(path):
        if path.lower().endswith((".fchk", ".fck", ".fch")):
            QTimer.singleShot(0, context.get_main_window(), lambda: open_fchk(path))
            return True
        return False

    context.register_drop_handler(handle_drop, priority=10)


def run(mw):
    from PySide6.QtWidgets import QFileDialog

    path, _ = QFileDialog.getOpenFileName(
        mw, "Open Gaussian FCHK", "", "Gaussian FCHK (*.fchk *.fck);;All Files (*)"
    )
    if not path:
        return

    context = PLUGIN_CONTEXT
    if not context:
        return
    try:
        dialog = OrbitalWidget(mw, context, path)
    except UnsupportedBasisError as exc:
        _warn_unsupported(mw, exc)
        return
    except (OSError, ValueError, TypeError, IndexError, RuntimeError) as exc:
        _warn_invalid(mw, exc)
        return
    dialog.show()
