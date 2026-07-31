"""
Headless GUI tests for the Orbital Comparator plugin.

Covers: OrbitalComparator, CubeSlot.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_orbital_comparator.py

numpy and pyvista are MagicMock here, so these tests deliberately do NOT parse
cube files — a stubbed numpy makes every geometry assertion vacuous, and the
real parsing/grid maths is covered against real numpy in
tests/test_plugin_orbital_comparator.py. What only real Qt can show is tested
here instead: that the signal wiring fires, that setValue on a live spin box
actually reaches render_all, and that widgets hold the state the code reads
back.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

COMPARATOR_PATH = PLUGINS_DIR / "Orbital_Comparator" / "orbital_comparator.py"

with mock_chemistry_imports():
    _comparator = load_plugin_for_gui(COMPARATOR_PATH)


class _FakeMesh:
    n_points = 10


class _FakeGrid:
    """Stands in for the parsed pyvista grid, without touching VTK."""

    def __init__(self):
        self.levels = []

    def contour(self, levels, scalars=None):
        self.levels.append(levels[0])
        return _FakeMesh()


# ===========================================================================
# OrbitalComparator  (visible plugin: "Orbital Comparator")
# ===========================================================================


class TestOrbitalComparator:
    @pytest.fixture
    def comp(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        w = _comparator.OrbitalComparator(context=ctx)
        w.mw = SimpleNamespace(plotter=MagicMock())
        yield w
        w.destroy()

    def _load(self, comp, index):
        """Put a slot into the loaded state without parsing anything."""
        slot = comp.slots[index]
        slot.grid = _FakeGrid()
        slot.path = f"/cubes/mo{index}.cube"
        slot.settings = {"grid": 40, "margin": 4.0, "version": "3.13.2"}
        slot.check_on.setChecked(True)
        return slot

    def _drawn(self, comp):
        return [
            c.kwargs["name"]
            for c in comp.mw.plotter.add_mesh.call_args_list
            if "name" in c.kwargs
        ]

    # -- construction ----------------------------------------------------

    def test_the_window_builds_with_four_slots(self, comp):
        assert len(comp.slots) == _comparator.SLOT_COUNT

    def test_the_window_has_a_title(self, comp):
        assert comp.windowTitle() == "Orbital Comparator"

    def test_it_accepts_drops(self, comp):
        assert comp.acceptDrops()

    def test_every_slot_starts_empty_and_unticked(self, comp):
        for slot in comp.slots:
            assert slot.grid is None
            assert not slot.check_on.isChecked()

    def test_every_slot_offers_the_styles(self, comp):
        for slot in comp.slots:
            assert slot.combo_style.count() == len(_comparator.STYLES)

    def test_the_isovalue_box_holds_a_real_value(self, comp):
        assert comp.slots[0].spin_iso.value() == pytest.approx(0.02)

    def test_slot_colours_round_trip_through_the_stylesheet(self, comp):
        comp.slots[0].set_color("p", "#123456")
        assert comp.slots[0].color("p") == "#123456"

    def test_slot_colours_start_distinct(self, comp):
        assert len({s.color("p") for s in comp.slots}) == _comparator.SLOT_COUNT

    # -- live wiring, which only real signals can show -------------------

    def test_setting_the_isovalue_redraws(self, comp):
        slot = self._load(comp, 0)
        comp.mw.plotter.add_mesh.reset_mock()
        slot.spin_iso.setValue(0.05)  # a real spin box emits valueChanged
        assert "orb_cmp0_p" in self._drawn(comp)

    def test_the_new_isovalue_reaches_the_contour(self, comp):
        slot = self._load(comp, 0)
        slot.grid.levels.clear()
        slot.spin_iso.setValue(0.05)
        assert slot.grid.levels[0] == pytest.approx(0.05)

    def test_setting_the_opacity_redraws(self, comp):
        slot = self._load(comp, 0)
        comp.mw.plotter.add_mesh.reset_mock()
        slot.spin_opacity.setValue(0.9)
        assert "orb_cmp0_p" in self._drawn(comp)

    def test_changing_the_style_redraws(self, comp):
        slot = self._load(comp, 0)
        comp.mw.plotter.add_mesh.reset_mock()
        slot.combo_style.setCurrentText("Wireframe")
        assert "orb_cmp0_p" in self._drawn(comp)

    def test_ticking_show_draws_the_slot(self, comp):
        slot = comp.slots[0]
        slot.grid = _FakeGrid()
        comp.mw.plotter.add_mesh.reset_mock()
        slot.check_on.setChecked(True)
        assert "orb_cmp0_p" in self._drawn(comp)

    def test_unticking_show_removes_the_actors(self, comp):
        slot = self._load(comp, 0)
        comp.mw.plotter.remove_actor.reset_mock()
        slot.check_on.setChecked(False)
        removed = [c[0][0] for c in comp.mw.plotter.remove_actor.call_args_list]
        assert "orb_cmp0_p" in removed

    def test_four_slots_draw_under_distinct_names(self, comp):
        for i in range(4):
            self._load(comp, i)
        comp.render_all()
        names = self._drawn(comp)
        for i in range(4):
            assert f"orb_cmp{i}_p" in names

    def test_syncing_the_isovalue_reaches_every_slot(self, comp):
        for i in range(4):
            self._load(comp, i)
        comp.slots[0].spin_iso.setValue(0.066)
        comp.sync_isovalue()
        assert all(s.spin_iso.value() == pytest.approx(0.066) for s in comp.slots)

    def test_the_status_line_counts_what_is_shown(self, comp):
        self._load(comp, 0)
        comp.render_all()
        assert "1 shown" in comp.lbl_status.text()

    def test_the_refresh_button_redraws_the_orbitals(self, comp):
        """A real click, through the real signal: the host can rebuild its
        scene at any time and drop the actors."""
        self._load(comp, 0)
        comp.mw.plotter.add_mesh.reset_mock()
        comp.btn_refresh.click()
        assert "orb_cmp0_p" in self._drawn(comp)

    # -- teardown --------------------------------------------------------

    def test_clearing_empties_every_slot(self, comp):
        for i in range(4):
            self._load(comp, i)
        comp.clear_all()
        assert all(s.grid is None for s in comp.slots)
        assert "No cubes" in comp.lbl_status.text()

    def test_closing_removes_the_actors(self, comp):
        self._load(comp, 0)
        comp.mw.plotter.remove_actor.reset_mock()
        comp.close()
        removed = [c[0][0] for c in comp.mw.plotter.remove_actor.call_args_list]
        assert "orb_cmp0_p" in removed
        assert "orb_cmp0_n" in removed
