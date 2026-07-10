"""
Tests for the High Resolution Imager plugin (initialize -> add_export_action; PLUGIN_CONTEXT stored).
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
HI_RES_PATH = PLUGINS_DIR / "High_Resolution_Imager" / "high_res_imager.py"


class TestHighResolutionImager:
    def test_initialize_registers_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_screenshot(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            label = ctx.add_export_action.call_args[0][0]
            assert "Screenshot" in label

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_initialize_does_not_register_menu_action(self):
        """High Resolution Imager uses add_export_action, not add_menu_action."""
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_not_called()




# ---------------------------------------------------------------------------
# run() guards / rejected dialog short-circuit
# ---------------------------------------------------------------------------

from unittest.mock import MagicMock

HIRES_PATH = HI_RES_PATH


class TestHighResImagerRun:
    def test_run_without_plugin_manager_returns(self):
        with mock_optional_imports():
            mod = load_plugin(HIRES_PATH)
            mod.PLUGIN_CONTEXT = MagicMock()
            mw = MagicMock(spec=[])  # no plugin_manager attribute
            mod.run(mw)  # must not raise, must not open dialog
        mod.PLUGIN_CONTEXT.get_main_window.assert_not_called()

    def test_run_without_context_returns(self):
        with mock_optional_imports():
            mod = load_plugin(HIRES_PATH)
            assert mod.PLUGIN_CONTEXT is None
            mod.run(MagicMock())  # must not raise

    def test_rejected_dialog_skips_file_dialog(self):
        with mock_optional_imports():
            mod = load_plugin(HIRES_PATH)
            ctx = make_context()
            # QDialog is a MagicMock: exec() returns a MagicMock which never
            # equals DialogCode.Accepted, i.e. the "user cancelled" path.
            mod.take_screenshot(ctx)
            assert mod.QFileDialog.getSaveFileName.call_count == 0
            ctx.plotter.screenshot.assert_not_called()


# ---------------------------------------------------------------------------
# Full "accepted" dialog flow — driven with controllable fake Qt widgets so
# the radio-button / spinbox / file-dialog interactions are deterministic
# (the default MagicMock-based widgets are indistinguishable singletons).
# ---------------------------------------------------------------------------


class _FakeSignal:
    def __init__(self):
        self._slots = []

    def connect(self, slot):
        self._slots.append(slot)

    def emit(self, *a):
        for s in self._slots:
            s(*a)


class _FakeRadioButton:
    def __init__(self, text=""):
        self.text = text
        self._checked = False

    def setChecked(self, v):
        self._checked = v

    def isChecked(self):
        return self._checked


class _FakeButtonGroup:
    def __init__(self, parent=None):
        self.buttons = []
        self.buttonClicked = _FakeSignal()

    def addButton(self, b):
        self.buttons.append(b)


class _FakeSpinBox:
    def __init__(self):
        self._v = 0

    def setRange(self, *a):
        pass

    def setValue(self, v):
        self._v = v

    def value(self):
        return self._v


class _FakeButton:
    def __init__(self, text=""):
        self.text = text
        self.clicked = _FakeSignal()

    def setEnabled(self, v):
        pass


class _FakeLayout:
    def addWidget(self, *a, **k):
        pass

    def addLayout(self, *a, **k):
        pass

    def addStretch(self, *a, **k):
        pass


class _FakeLabel:
    def __init__(self, *a, **k):
        pass

    def setFixedSize(self, *a):
        pass

    def setStyleSheet(self, *a):
        pass


class _FakeDialog:
    class DialogCode:
        Accepted = 1
        Rejected = 0

    result = DialogCode.Accepted  # overridable default

    def __init__(self, parent=None):
        self.accept = MagicMock()
        self.reject = MagicMock()

    def setWindowTitle(self, *a):
        pass

    def setLayout(self, *a):
        pass

    def exec(self):
        return _FakeDialog.result


def _patch_fake_widgets(mod):
    """Replace the module-level Qt widget classes with lightweight, fully
    controllable fakes so take_screenshot()'s accepted-dialog branch can be
    driven end to end."""
    mod.QDialog = _FakeDialog
    mod.QVBoxLayout = lambda *a, **k: _FakeLayout()
    mod.QHBoxLayout = lambda *a, **k: _FakeLayout()
    mod.QLabel = _FakeLabel
    mod.QSpinBox = _FakeSpinBox
    mod.QPushButton = _FakeButton

    import PyQt6.QtWidgets as qtw  # same shared mocked module as the plugin

    qtw.QRadioButton = _FakeRadioButton
    qtw.QButtonGroup = _FakeButtonGroup


class TestTakeScreenshotAcceptedFlow:
    def _run(self, tmp_path, filename, raise_on_screenshot=False):
        """Drive take_screenshot() through the accepted-dialog branch.

        The plugin itself calls ``rb_trans.setChecked(True)`` right after
        constructing the radio buttons, so the default (untouched) run
        exercises the "Transparent Background" path — matching the dialog's
        real default.
        """
        with mock_optional_imports():
            mod = load_plugin(HIRES_PATH)
            _patch_fake_widgets(mod)
            _FakeDialog.result = _FakeDialog.DialogCode.Accepted
            mod.QFileDialog.getSaveFileName = MagicMock(return_value=(filename, ""))
            mod.QMessageBox.information = MagicMock()
            mod.QMessageBox.critical = MagicMock()

            ctx = make_context()
            ctx.plotter = MagicMock()
            mw = ctx.get_main_window.return_value
            mw.init_manager.settings.get.return_value = "#4f4f4f"

            if raise_on_screenshot:
                ctx.plotter.screenshot.side_effect = RuntimeError("gpu oom")

            mod.take_screenshot(ctx)
        return mod, ctx, mw

    def test_empty_filename_skips_screenshot(self, tmp_path):
        mod, ctx, mw = self._run(tmp_path, filename="")
        ctx.plotter.screenshot.assert_not_called()
        ctx.plotter.save_graphic.assert_not_called()

    def test_raster_screenshot_called_with_window_size(self, tmp_path):
        out = str(tmp_path / "shot.png")
        mod, ctx, mw = self._run(tmp_path, filename=out)
        ctx.plotter.screenshot.assert_called_once()
        args, kwargs = ctx.plotter.screenshot.call_args
        assert args[0] == out
        assert "window_size" in kwargs

    def test_svg_export_uses_save_graphic_not_screenshot(self, tmp_path):
        out = str(tmp_path / "shot.svg")
        mod, ctx, mw = self._run(tmp_path, filename=out)
        ctx.plotter.save_graphic.assert_called_once_with(out)
        ctx.plotter.screenshot.assert_not_called()

    def test_screenshot_exception_shows_critical_not_raise(self, tmp_path):
        out = str(tmp_path / "shot.png")
        mod, ctx, mw = self._run(tmp_path, filename=out, raise_on_screenshot=True)
        mod.QMessageBox.critical.assert_called_once()
        mod.QMessageBox.information.assert_not_called()
