"""
Headless GUI tests for the (retired) Gaussian Input Generator Neo plugin.

Gaussian Input Generator Neo has been replaced by Gaussian Input Generator
Pro; this whole class is skipped and kept only for historical continuity of
the test count.
"""

from __future__ import annotations

import pytest


# ===========================================================================
# RouteBuilderDialog  (visible plugin: "Gaussian Input Generator Neo")
# ===========================================================================


@pytest.mark.skip(reason="Gaussian Input Generator Neo retired; replaced by Gaussian Input Generator Pro")
class TestRouteBuilderDialog:
    """RouteBuilderDialog — tabbed QDialog with combo-driven route preview."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _gaussian.RouteBuilderDialog(parent=None, current_route="")
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Route Builder"

    def test_has_four_tabs(self, dlg):
        assert dlg.tabs.count() == 4

    def test_tab_method_basis_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Method" in n or "Basis" in n for n in names)

    def test_tab_job_type_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Job" in n for n in names)

    def test_tab_properties_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Propert" in n for n in names)

    def test_preview_label_exists(self, dlg):
        assert hasattr(dlg, "preview_label")

    def test_preview_non_empty_after_init(self, dlg):
        assert dlg.preview_label.text() != ""

    def test_preview_contains_default_basis_set(self, dlg):
        assert "6-31G" in dlg.preview_label.text()

    def test_ok_button_label(self, dlg):
        assert dlg.btn_ok.text() == "Apply to Job"

    def test_cancel_button_label(self, dlg):
        assert dlg.btn_cancel.text() == "Cancel"

    def test_method_type_combo_has_dft(self, dlg):
        items = [dlg.method_type.itemText(i) for i in range(dlg.method_type.count())]
        assert any("DFT" in it for it in items)

    def test_dft_method_list_contains_b3lyp(self, dlg):
        dlg.method_type.setCurrentText("DFT")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "B3LYP" in items

    def test_mp2_method_list_contains_ccsd(self, dlg):
        dlg.method_type.setCurrentText("MP2")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "CCSD" in items

    def test_hf_method_list_contains_hf(self, dlg):
        dlg.method_type.setCurrentText("Hartree-Fock")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "HF" in items

    def test_switching_method_updates_preview(self, dlg):
        dlg.method_type.setCurrentText("Hartree-Fock")
        dlg.method_name.setCurrentText("HF")
        assert "HF" in dlg.preview_label.text()

    def test_basis_set_change_updates_preview(self, dlg):
        dlg.method_type.setCurrentText("DFT")
        dlg.basis_set.setCurrentText("def2TZVP")
        assert "def2TZVP" in dlg.preview_label.text()

    def test_job_opt_only_shows_opt_group_hides_freq(self, dlg):
        # "Optimization Only (Opt)" is index 1
        dlg.job_type.setCurrentIndex(1)
        # Use isHidden() — isVisible() requires the window to be .show()n first
        assert not dlg.opt_group.isHidden()
        assert dlg.freq_group.isHidden()

    def test_job_freq_only_shows_freq_group_hides_opt(self, dlg):
        # "Frequency Only (Freq)" is index 2
        dlg.job_type.setCurrentIndex(2)
        assert dlg.opt_group.isHidden()
        assert not dlg.freq_group.isHidden()

    def test_job_sp_hides_both_groups(self, dlg):
        # "Single Point Energy (SP)" is index 3
        dlg.job_type.setCurrentIndex(3)
        assert dlg.opt_group.isHidden()
        assert dlg.freq_group.isHidden()

    def test_opt_freq_job_shows_both_groups(self, dlg):
        # "Optimization + Freq (Opt Freq)" is index 0
        dlg.job_type.setCurrentIndex(0)
        assert not dlg.opt_group.isHidden()
        assert not dlg.freq_group.isHidden()
