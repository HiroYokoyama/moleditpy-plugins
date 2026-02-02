from .gui import OrbitalWidget

# --- Plugin Metadata ---
PLUGIN_NAME = "Gaussian MO Analyzer"
PLUGIN_VERSION = "2026.02.02"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Visualizes Molecular Orbitals from Gaussian FCHK files by generating Cube files."

def initialize(context):
    """
    Register the FCHK file opener.
    """
    def open_fchk(path):
        mw = context.get_main_window()
        # Create and show dialog (keep reference to avoid GC)
        # Attach to context to keep alive
        context._fchk_dialog = OrbitalWidget(mw, context, path)
        context._fchk_dialog.show() 
    
    context.register_file_opener(".fchk", open_fchk, priority=10)
    context.register_file_opener(".fck", open_fchk, priority=10)
    context.register_file_opener(".fch", open_fchk, priority=10)

    def handle_drop(path):
        if path.lower().endswith((".fchk", ".fck", ".fch")):
            open_fchk(path)
            return True
        return False
        
    context.register_drop_handler(handle_drop, priority=10)
