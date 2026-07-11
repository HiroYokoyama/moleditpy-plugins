"""
Tests for the Dummy Atom Mode plugin (initialize -> add_toolbar_action, no add_menu_action).
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
DUMMY_PATH = PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py"


class TestDummyAtomMode:
    def test_initialize_calls_add_toolbar_action(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_toolbar_action.assert_called_once()

    def test_initialize_toolbar_text_is_dummy_atom_star(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            kwargs = ctx.add_toolbar_action.call_args[1]
            assert kwargs.get("text") == "Dummy Atom *"

    def test_initialize_does_not_call_add_menu_action(self):
        """Dummy Atom Mode registers a toolbar button, NOT a menu item."""
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_not_called()

    def test_initialize_tooltip_mentions_dummy(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            kwargs = ctx.add_toolbar_action.call_args[1]
            assert "dummy" in kwargs.get("tooltip", "").lower()


class TestDummyAtomChemPatch:
    """initialize() wraps Chem.Atom so '*' maps to atomic number 0 while all
    other symbols pass through to the previous constructor."""

    def _initialized_mod(self):
        ctx = make_context()
        mod = load_plugin(PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py")
        original = MagicMock(name="previous_Atom")
        mod.Chem = MagicMock()
        mod.Chem.Atom = original
        del mod.Chem._original_Atom  # force the guard branch to store it
        mod.initialize(ctx)
        return mod, ctx, original

    def test_star_symbol_maps_to_atomic_number_zero(self):
        with mock_optional_imports():
            mod, ctx, original = self._initialized_mod()
            # initialize() stored the pre-patch constructor as _original_Atom,
            # so '*' must reach it as atomic number 0, not as a string.
            mod.Chem.Atom("*")
            original.assert_called_once_with(0)

    def test_regular_symbol_passes_through(self):
        with mock_optional_imports():
            mod, ctx, original = self._initialized_mod()
            mod.Chem.Atom("C")
            original.assert_called_once_with("C")

    def test_cpk_colors_registered_for_star(self):
        with mock_optional_imports():
            mod, ctx, _ = self._initialized_mod()
            mod.CPK_COLORS.__setitem__.assert_any_call("*", mod.QColor.return_value)
            mod.CPK_COLORS_PV.__setitem__.assert_any_call("*", [0.5, 0.5, 0.5])


class _FakeAction:
    def __init__(self, text):
        self._text = text
        self.checked = []
        self.checkable = []

    def text(self):
        return self._text

    def setChecked(self, flag):
        self.checked.append(flag)

    def setCheckable(self, flag):
        self.checkable.append(flag)


class TestDummyAtomSetModeSync:
    """The patched UIManager.set_mode keeps the toolbar button checked state
    in sync with the active scene mode."""

    def _patched_set_mode(self):
        ctx = make_context()
        mod = load_plugin(PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py")
        previous = MagicMock(name="previous_set_mode")
        mod.UIManager = MagicMock()
        mod.UIManager.set_mode = previous
        mod.initialize(ctx)
        return mod.UIManager.set_mode, previous

    def _fake_ui_manager(self, action):
        fake = MagicMock()
        fake.host.init_manager.plugin_toolbar.actions.return_value = [action]
        return fake

    def test_dummy_mode_checks_button_and_chains_previous(self):
        with mock_optional_imports():
            set_mode, previous = self._patched_set_mode()
            action = _FakeAction("Dummy Atom *")
            fake_self = self._fake_ui_manager(action)
            set_mode(fake_self, "atom_*")
            previous.assert_called_once_with(fake_self, "atom_*")
            assert action.checked == [True]

    def test_other_mode_unchecks_button(self):
        with mock_optional_imports():
            set_mode, previous = self._patched_set_mode()
            action = _FakeAction("Dummy Atom *")
            set_mode(self._fake_ui_manager(action), "select")
            assert action.checked == [False]

    def test_unrelated_actions_left_alone(self):
        with mock_optional_imports():
            set_mode, previous = self._patched_set_mode()
            action = _FakeAction("Some Other Tool")
            set_mode(self._fake_ui_manager(action), "atom_*")
            assert action.checked == []


class TestDummyAtomToggleCallback:
    def _toggle_and_mw(self):
        ctx = make_context()
        mod = load_plugin(PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py")
        mod.initialize(ctx)
        toggle = ctx.add_toolbar_action.call_args[1]["callback"]
        return toggle, ctx.get_main_window.return_value

    def test_enters_dummy_mode_when_inactive(self):
        with mock_optional_imports():
            toggle, mw = self._toggle_and_mw()
            mw.init_manager.scene.mode = "select"
            toggle()
            mw.ui_manager.set_mode_and_update_toolbar.assert_called_once_with("atom_*")

    def test_toggles_back_to_select_when_active(self):
        with mock_optional_imports():
            toggle, mw = self._toggle_and_mw()
            mw.init_manager.scene.mode = "atom_*"
            toggle()
            mw.ui_manager.set_mode_and_update_toolbar.assert_called_once_with("select")


class TestDummyAtomDeferredCheckable:
    def test_timer_callback_makes_action_checkable(self):
        with mock_optional_imports():
            ctx = make_context()
            mod = load_plugin(PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py")
            mod.initialize(ctx)

            import PyQt6.QtCore as qtcore  # same mocked module the plugin imported

            delay, callback = qtcore.QTimer.singleShot.call_args[0]
            assert callable(callback)

            action = _FakeAction("Dummy Atom *")
            mw = ctx.get_main_window.return_value
            mw.init_manager.plugin_toolbar.actions.return_value = [action]
            mw.init_manager.scene.mode = "atom_*"
            callback()
            assert action.checkable == [True]
            assert action.checked == [True]
