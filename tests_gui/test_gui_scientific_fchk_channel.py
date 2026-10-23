"""An unrestricted file must identify the supported spin channel explicitly."""

from unittest.mock import MagicMock, patch

from test_gui_plugin_gaussian_mo_analyzer import (
    _S_SHELL_MO_FCHK,
    _mo,
    _real_numpy_in_mo_modules,
)


def test_unrestricted_result_is_explicitly_alpha_only(qapp, tmp_path):
    source = tmp_path / "uhf.fchk"
    source.write_text(
        _S_SHELL_MO_FCHK
        + "Beta MO coefficients                        R   N=           1\n"
        + " -2.000000000000E+00\n"
        + "Number of alpha electrons                   I            1\n"
        + "Number of beta electrons                    I            0\n",
        encoding="utf-8",
    )
    with _real_numpy_in_mo_modules() as np:
        window = _mo.OrbitalWidget(None, MagicMock(), str(source))
        try:
            assert "Alpha orbitals only" in window.windowTitle()
            assert window.list_widget.item(0).text().startswith("Alpha MO")
            assert "HOMO" in window.list_widget.item(0).text()
            assert "Beta orbitals are not supported" in window.list_widget.toolTip()
            # The displayed alpha row must route alpha coefficients to the
            # worker; the deliberately different beta array exposes confusion.
            worker = MagicMock()
            with (
                patch.object(_mo.gui, "CalcWorker", return_value=worker) as factory,
                patch.object(_mo.gui.os.path, "exists", return_value=False),
            ):
                window.generation_queue = [1]
                window.process_generation_queue()
            np.testing.assert_array_equal(factory.call_args.args[5], [1.0])
            window.worker = None
        finally:
            window.destroy()
