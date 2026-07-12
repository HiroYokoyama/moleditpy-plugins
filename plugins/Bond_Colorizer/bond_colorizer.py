from PyQt6.QtWidgets import (
    QVBoxLayout,
    QPushButton,
    QLabel,
    QColorDialog,
    QMessageBox,
    QLineEdit,
    QGroupBox,
    QDialog,
)
from PyQt6.QtGui import QColor, QCloseEvent
from PyQt6.QtCore import Qt, QTimer
import logging

# Plugin Metadata
PLUGIN_NAME = "Bond Colorizer"
PLUGIN_VERSION = "2026.07.12"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Select bonds by index or atom pair in the 3D viewer and apply color."
PLUGIN_CONTEXT = None


class BondColorizerWindow(QDialog):
    """
    Plugin window for applying custom colors to bonds in the 3D scene.
    Supports both bond-index input and atom-pair input.
    """

    def __init__(self, context):
        super().__init__(parent=context.get_main_window())
        self.context = context
        self._restore_measurement_mode = False
        self._restore_edit_mode = False
        self._forced_measurement_mode = False

        # Set window properties
        self.setModal(False)
        self.setWindowTitle(PLUGIN_NAME)
        self.setWindowFlags(Qt.WindowType.Window)

        self.current_color = QColor(255, 0, 0)  # Default red
        self.init_ui()

        # Register window for V3 lifecycle management
        self.context.register_window("main_panel", self)

        # Ensure the host "3D Select" (measurement) mode is active while this
        # window is open, so clicking atoms in the viewer populates a selection.
        self._enter_select_mode()

        self.context.show_status_message("Bond Colorizer: 3D picking enabled.")

    def _enter_select_mode(self):
        """Remember the current 3D interaction mode and switch to select mode."""
        mw = self.context.get_main_window()
        if not mw or not hasattr(mw, "edit_3d_manager"):
            return
        try:
            e3d = mw.edit_3d_manager
            self._restore_measurement_mode = bool(
                getattr(e3d, "measurement_mode", False)
            )
            self._restore_edit_mode = bool(getattr(e3d, "is_3d_edit_mode", False))
            if not self._restore_measurement_mode:
                if hasattr(mw, "init_manager") and hasattr(
                    mw.init_manager, "measurement_action"
                ):
                    mw.init_manager.measurement_action.setChecked(True)
                e3d.toggle_measurement_mode(True)
                self._forced_measurement_mode = True
        except Exception as _e:
            logging.warning("[bond_colorizer] enter select mode silenced: %s", _e)

    def _restore_select_mode(self):
        """Clear any selection this window drove and restore the prior 3D mode.

        Runs unconditionally on close so the viewer never gets stuck in a
        half-selected / non-interactive state (VTK interactor reset happens
        inside toggle_measurement_mode / toggle_3d_edit_mode).
        """
        mw = self.context.get_main_window()
        if not mw or not hasattr(mw, "edit_3d_manager"):
            return
        e3d = mw.edit_3d_manager
        try:
            # Always drop the picked-atom selection and its 3D highlight.
            if hasattr(e3d, "clear_measurement_selection"):
                e3d.clear_measurement_selection()
            if hasattr(e3d, "selected_atoms_3d"):
                e3d.selected_atoms_3d.clear()
                if hasattr(e3d, "update_3d_selection_display"):
                    e3d.update_3d_selection_display()
        except Exception as _e:
            logging.warning("[bond_colorizer] clear selection silenced: %s", _e)
        try:
            if hasattr(mw, "init_manager") and hasattr(
                mw.init_manager, "measurement_action"
            ):
                mw.init_manager.measurement_action.setChecked(
                    self._restore_measurement_mode
                )
            e3d.toggle_measurement_mode(self._restore_measurement_mode)
            if hasattr(mw, "ui_manager"):
                mw.ui_manager.toggle_3d_edit_mode(self._restore_edit_mode)
        except Exception as _e:
            logging.warning("[bond_colorizer] restore mode silenced: %s", _e)

    def init_ui(self):
        layout = QVBoxLayout()

        # Information Label
        info_label = QLabel("Select bonds in the 3D viewer and apply color.")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)

        # Selection Group
        sel_group = QGroupBox("Selection")
        sel_layout = QVBoxLayout()

        self.le_bond_ids = QLineEdit()
        self.le_bond_ids.setPlaceholderText("Bond IDs: e.g. 0, 3, 5")
        sel_layout.addWidget(self.le_bond_ids)

        self.le_atom_pairs = QLineEdit()
        self.le_atom_pairs.setPlaceholderText("Atom pairs: e.g. 0-1, 2-3")
        sel_layout.addWidget(self.le_atom_pairs)

        # Auto-update timer
        self.sel_timer = QTimer(self)
        self.sel_timer.timeout.connect(self._auto_update_selection)
        self.sel_timer.start(200)  # Check every 200ms

        sel_group.setLayout(sel_layout)
        layout.addWidget(sel_group)

        # Color Group
        col_group = QGroupBox("Color")
        col_layout = QVBoxLayout()

        self.btn_color = QPushButton("Choose Color")
        # Initial style based on self.current_color (QColor object)
        self.btn_color.setStyleSheet(
            f"background-color: {self.current_color.name()}; color: {'black' if self.current_color.lightness() > 128 else 'white'};"
        )
        self.btn_color.clicked.connect(self.choose_color)
        col_layout.addWidget(self.btn_color)

        btn_apply = QPushButton("Apply Color")
        btn_apply.clicked.connect(self.apply_color)
        col_layout.addWidget(btn_apply)

        col_group.setLayout(col_layout)
        layout.addWidget(col_group)

        # Reset Group
        reset_group = QGroupBox("Reset")
        reset_layout = QVBoxLayout()

        btn_reset = QPushButton("Reset to Default Bond Color")
        btn_reset.clicked.connect(self.reset_colors)
        reset_layout.addWidget(btn_reset)

        reset_group.setLayout(reset_layout)
        layout.addWidget(reset_group)

        layout.addStretch()

        # Close Button
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)

        self.setLayout(layout)
        self.resize(300, 400)

    def get_selection_from_viewer(self):
        """Get selected atom indices from context/Edit3DManager, then derive
        the bonds whose both endpoints lie in the selection.

        Wrapped defensively: a transient error here must never propagate out of
        the polling timer, otherwise selection sync would stop for good.
        """
        try:
            # Fetch 2D indices from core API
            indices = set(self.context.get_selected_atom_indices())

            # Fetch 3D indices directly from Edit3DManager (measurement + select).
            mw = self.context.get_main_window()
            if mw and hasattr(mw, "edit_3d_manager"):
                indices_3d = getattr(mw.edit_3d_manager, "selected_atoms_3d", set())
                if isinstance(indices_3d, (set, list)):
                    indices.update(indices_3d)
                picked_for_measurement = getattr(
                    mw.edit_3d_manager, "selected_atoms_for_measurement", []
                )
                if isinstance(picked_for_measurement, (set, list, tuple)):
                    indices.update(int(i) for i in picked_for_measurement)

            mol = self.context.current_molecule
            if not mol:
                return

            pairs = set()
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if a1 in indices and a2 in indices:
                    pairs.add((min(a1, a2), max(a1, a2)))

            sorted_pairs = sorted(pairs)
            new_text = ",".join(f"{a}-{b}" for a, b in sorted_pairs)
            if self.le_atom_pairs.text() != new_text:
                self.le_atom_pairs.setText(new_text)
        except Exception as _e:
            logging.warning("[bond_colorizer] selection sync silenced: %s", _e)

    def _auto_update_selection(self):
        """Timer slot to auto-update selection."""
        if self.le_bond_ids.hasFocus() or self.le_atom_pairs.hasFocus():
            return
        self.get_selection_from_viewer()

    def choose_color(self):
        c = QColorDialog.getColor(initial=self.current_color, title="Select Color")
        if c.isValid():
            self.current_color = c
            # Update button style
            self.btn_color.setStyleSheet(
                f"background-color: {c.name()}; color: {'black' if c.lightness() > 128 else 'white'};"
            )

    def apply_color(self):
        bond_ids_txt = self.le_bond_ids.text().strip()
        atom_pairs_txt = self.le_atom_pairs.text().strip()

        if not bond_ids_txt and not atom_pairs_txt:
            QMessageBox.warning(
                self,
                "Warning",
                "No bonds selected. Please select bonds in the 3D viewer first.",
            )
            return

        mol = self.context.current_molecule
        if not mol:
            QMessageBox.warning(self, "Error", "No molecule loaded.")
            return

        target_indices = set()
        bad_ids = []
        bad_pairs = []

        if bond_ids_txt:
            try:
                str_ids = [x.strip() for x in bond_ids_txt.split(",") if x.strip()]
                ids = [int(x) for x in str_ids]
            except ValueError:
                QMessageBox.warning(self, "Error", "Invalid bond ID format.")
                return

            for idx in ids:
                if 0 <= idx < mol.GetNumBonds():
                    target_indices.add(idx)
                else:
                    bad_ids.append(idx)

        if atom_pairs_txt:
            str_pairs = [x.strip() for x in atom_pairs_txt.split(",") if x.strip()]
            pairs = []
            for token in str_pairs:
                parts = token.split("-")
                if len(parts) != 2:
                    QMessageBox.warning(
                        self, "Error", "Invalid atom-pair format. Use a-b, c-d."
                    )
                    return
                try:
                    a1, a2 = int(parts[0]), int(parts[1])
                except ValueError:
                    QMessageBox.warning(
                        self, "Error", "Invalid atom-pair format. Use a-b, c-d."
                    )
                    return
                pairs.append((a1, a2))

            for a1, a2 in pairs:
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is None:
                    bad_pairs.append((a1, a2))
                else:
                    target_indices.add(bond.GetIdx())

        try:
            controller = self.context.get_3d_controller()
            hex_color = self.current_color.name()

            for idx in sorted(target_indices):
                controller.set_bond_color(idx, hex_color)

            self.context.refresh_3d_view()

        except Exception as e:
            logging.exception("Failed to apply color: %s", e)
            QMessageBox.critical(self, "Error", f"Failed to apply color: {e}")
            return

        if bad_ids or bad_pairs:
            QMessageBox.warning(
                self,
                "Warning",
                f"Skipped invalid: bonds {bad_ids}, pairs {bad_pairs}",
            )

    def reset_colors(self):
        mol = self.context.current_molecule
        if not mol:
            return

        try:
            controller = self.context.get_3d_controller()
            for i in range(mol.GetNumBonds()):
                controller.set_bond_color(i, None)

            self.context.refresh_3d_view()
            self.le_bond_ids.clear()
            self.le_atom_pairs.clear()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reset colors: {e}")

    def closeEvent(self, event: QCloseEvent):
        """Reset selection and restore the prior 3D interaction mode on close."""
        try:
            if self.sel_timer.isActive():
                self.sel_timer.stop()
        except Exception as _e:
            logging.warning("[bond_colorizer] stop timer silenced: %s", _e)

        # Always reset — even if this window did not force select mode — so the
        # viewer is never left stuck with a stale selection / interactor state.
        self._restore_select_mode()
        super().closeEvent(event)


def launch(context):
    win = context.get_window("main_panel")
    if win:
        win.show()
        win.raise_()
        win.activateWindow()
        return

    win = BondColorizerWindow(context)
    win.show()


def initialize(context):
    """MoleditPy Plugin Entry Point (V3.0)"""
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context

    def save_handler():
        """Save bond color overrides to project file."""
        mw = context.get_main_window()
        v3d = getattr(mw, "view_3d_manager", None)
        if not v3d or not hasattr(v3d, "_plugin_bond_color_overrides"):
            return {}

        return {
            "bond_colors": {
                str(k): v for k, v in v3d._plugin_bond_color_overrides.items()
            }
        }

    def load_handler(data):
        """Restore bond color overrides from project file."""
        if not data:
            return

        bond_colors = data.get("bond_colors", {})
        controller = context.get_3d_controller()
        for bond_idx_str, hex_color in bond_colors.items():
            try:
                controller.set_bond_color(int(bond_idx_str), hex_color)
            except Exception as _e:
                logging.warning("[bond_colorizer.py:268] silenced: %s", _e)
        context.refresh_3d_view()

    def on_document_reset():
        """Reset colors when a new document is created."""
        mw = context.get_main_window()
        v3d = getattr(mw, "view_3d_manager", None)
        if v3d and hasattr(v3d, "_plugin_bond_color_overrides"):
            v3d._plugin_bond_color_overrides.clear()

        win = context.get_window("main_panel")
        if win:
            win.le_bond_ids.clear()
            win.le_atom_pairs.clear()

    context.register_save_handler(save_handler)
    context.register_load_handler(load_handler)
    context.register_document_reset_handler(on_document_reset)
    context.show_status_message(f"{PLUGIN_NAME} Loaded.")


def run(mw):
    """Legacy entry point for older plugin loaders."""
    if hasattr(mw, "host"):
        mw = mw.host

    if not hasattr(mw, "plugin_manager"):
        return

    context = PLUGIN_CONTEXT
    if not context:
        return
    launch(context)
