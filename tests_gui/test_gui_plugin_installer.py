"""
Headless GUI tests for the Plugin Installer plugin.

Covers: PluginDetailsDialog, PluginInstallerWindow, and pure-function unit
tests (is_app_version_compatible, parse_dependency, check_dependency_satisfied).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PLUGIN_INSTALLER_PATH = PLUGINS_DIR / "Plugin_Installer" / "plugin_installer.py"

with mock_chemistry_imports():
    _plugin_installer = load_plugin_for_gui(PLUGIN_INSTALLER_PATH)


# ===========================================================================
# PluginDetailsDialog  (inside Plugin Installer)
# ===========================================================================


def _make_details_dialog(*, dependencies=None, local_info=None) -> object:
    """Create a PluginDetailsDialog with safe defaults."""
    return _plugin_installer.PluginDetailsDialog(
        None,
        "Test Plugin",
        "Test Author",
        "2026.01.01",
        "A headless test plugin.",
        dependencies if dependencies is not None else [],
        local_info,
        None,
        ">=4.0.0",
    )


class TestPluginDetailsDialog:
    """PluginDetailsDialog — lightweight details popup, no network calls."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _make_details_dialog()
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title_includes_plugin_name(self, dlg):
        assert "Test Plugin" in dlg.windowTitle()

    def test_no_error_with_empty_dependencies(self, qapp):
        d = _make_details_dialog(dependencies=[])
        assert d is not None
        d.destroy()

    def test_no_error_with_missing_dependency(self, qapp):
        d = _make_details_dialog(dependencies=["fake-nonexistent-pkg-xyz"])
        assert d is not None
        d.destroy()

    def test_no_error_with_installed_dependency(self, qapp):
        d = _make_details_dialog(dependencies=["pytest"])
        assert d is not None
        d.destroy()

    def test_target_file_stored(self, qapp):
        d = _plugin_installer.PluginDetailsDialog(
            None, "X", "A", "1.0", "Desc", [], None, "/some/path.py"
        )
        assert d.target_file == "/some/path.py"
        d.destroy()


# ===========================================================================
# PluginInstallerWindow  (visible plugin: "Plugin Installer")
# ===========================================================================


class TestPluginInstallerWindow:
    """PluginInstallerWindow — main window of the Plugin Installer.

    check_updates() is patched to prevent any network calls during __init__.
    """

    @pytest.fixture
    def installer(self, qapp):
        from unittest.mock import patch

        with patch.object(
            _plugin_installer.PluginInstallerWindow,
            "check_updates",
            lambda self: None,
        ):
            d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
        yield d
        d.destroy()

    def test_creates_without_error(self, installer):
        assert installer is not None

    def test_window_title(self, installer):
        assert installer.windowTitle() == "Plugin Installer"

    def test_table_has_six_columns(self, installer):
        assert installer.table.columnCount() == 6

    def test_search_input_exists(self, installer):
        assert installer.search_input.placeholderText() != ""

    def test_startup_checkbox_exists(self, installer):
        assert installer.chk_startup is not None


# ===========================================================================
# Plugin Installer — pure-function unit tests (no Qt needed)
# ===========================================================================


class TestPluginInstallerVersionCompat:
    """is_app_version_compatible(): pure Python, no Qt required."""

    def test_within_range(self):
        assert _plugin_installer.is_app_version_compatible("4.2.0", ">=4.0.0,<5.0.0")

    def test_equal_lower_bound(self):
        assert _plugin_installer.is_app_version_compatible("4.0.0", ">=4.0.0")

    def test_below_lower_bound(self):
        assert not _plugin_installer.is_app_version_compatible("3.9.9", ">=4.0.0")

    def test_equal_upper_bound_exclusive(self):
        assert not _plugin_installer.is_app_version_compatible("5.0.0", "<5.0.0")

    def test_wildcard_always_passes(self):
        assert _plugin_installer.is_app_version_compatible("1.0.0", "*")

    def test_empty_specifier_always_passes(self):
        assert _plugin_installer.is_app_version_compatible("1.0.0", "")

    def test_not_equal(self):
        assert not _plugin_installer.is_app_version_compatible("4.0.0", "!=4.0.0")

    def test_compatible_release_same_major(self):
        assert _plugin_installer.is_app_version_compatible("4.2.0", "~=4.0")

    def test_compatible_release_different_major(self):
        assert not _plugin_installer.is_app_version_compatible("5.0.0", "~=4.0")


class TestPluginInstallerParseDependency:
    """parse_dependency(): pure Python, no Qt required."""

    def test_name_only(self):
        name, spec = _plugin_installer.parse_dependency("numpy")
        assert name == "numpy"
        assert spec == ""

    def test_name_with_specifier(self):
        name, spec = _plugin_installer.parse_dependency("rdkit>=2022.03.5")
        assert name == "rdkit"
        assert "2022" in spec

    def test_empty_string_returns_empty_pair(self):
        name, spec = _plugin_installer.parse_dependency("")
        assert name == ""
        assert spec == ""

    def test_url_with_colon_is_rejected(self):
        name, spec = _plugin_installer.parse_dependency("http://example.com/package")
        assert name == ""

    def test_dash_in_name(self):
        name, spec = _plugin_installer.parse_dependency("my-package>=1.0")
        assert name == "my-package"


class TestPluginInstallerCheckDependency:
    """check_dependency_satisfied(): queries importlib.metadata, pure Python."""

    def test_installed_package_returns_true(self):
        assert _plugin_installer.check_dependency_satisfied("pytest")

    def test_nonexistent_package_returns_false(self):
        assert not _plugin_installer.check_dependency_satisfied(
            "fake-nonexistent-dep-xyz123"
        )

    def test_empty_string_returns_false(self):
        assert not _plugin_installer.check_dependency_satisfied("")


# ===========================================================================
# Plugin Installer — PluginInstallerWindow more widget tests
# ===========================================================================


class TestPluginInstallerWindowWidgets:
    """Additional widget-state tests for PluginInstallerWindow."""

    @pytest.fixture
    def installer(self, qapp):
        from unittest.mock import patch

        with patch.object(
            _plugin_installer.PluginInstallerWindow,
            "check_updates",
            lambda self: None,
        ):
            d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
        yield d
        d.destroy()

    def test_table_header_first_column(self, installer):
        assert installer.table.horizontalHeaderItem(0).text() == "Plugin Name"

    def test_table_header_action_column(self, installer):
        assert installer.table.horizontalHeaderItem(5).text() == "Action"

    def test_progress_bar_initially_hidden(self, installer):
        assert installer._progress_bar.isHidden()

    def test_upgrade_button_initially_hidden(self, installer):
        assert installer.btn_upgrade_app.isHidden()

    def test_update_all_initially_disabled(self, installer):
        assert not installer.btn_update_all.isEnabled()

    def test_table_empty_before_fetch(self, installer):
        assert installer.table.rowCount() == 0

    def test_remote_data_initially_empty(self, installer):
        assert installer.remote_data == []


# ===========================================================================
# Plugin Installer — supported_python_version compatibility (GUI)
# ===========================================================================

PY_OK_SPEC = ">=3.9, <3.15"
PY_BAD_SPEC = ">=99.0"
APP_OK_SPEC = ">=4.0.0, <5.0.0"
APP_BAD_SPEC = ">=99.0"


def _remote_entry(name, *, moleditpy_spec=APP_OK_SPEC, python_spec=None):
    entry = {
        "name": name,
        "visible": True,
        "version": "1.0.0",
        "author": "Test Author",
        "description": "d",
        "supported_moleditpy_version": moleditpy_spec,
        "downloadUrl": "../plugins/Fake/fake.py",
    }
    if python_spec is not None:
        entry["supported_python_version"] = python_spec
    return entry


class TestPythonVersionCompatGui:
    """populate_table() + install path honour supported_python_version."""

    @pytest.fixture
    def installer(self, qapp):
        from unittest.mock import MagicMock, patch

        mw = MagicMock()
        mw.plugin_manager.plugins = []
        with patch.object(
            _plugin_installer.PluginInstallerWindow,
            "check_updates",
            lambda self: None,
        ):
            d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
        d.main_window = mw
        d.get_app_version = lambda: "4.2.0"
        d.latest_app_version = "Checking..."
        yield d
        d.destroy()

    @staticmethod
    def _row_of(installer, name):
        for row in range(installer.table.rowCount()):
            if installer.table.item(row, 0).text() == name:
                return row
        raise AssertionError(f"row for {name!r} not found")

    def _populate(self, installer, entries):
        installer.remote_data = entries
        installer.populate_table(silent=True)

    def test_python_incompatible_marks_incompatible(self, installer):
        self._populate(
            installer, [_remote_entry("Py Bad", python_spec=PY_BAD_SPEC)]
        )
        row = self._row_of(installer, "Py Bad")
        status = installer.table.item(row, 4)
        assert status.text() == "Incompatible"
        assert "Python" in status.toolTip()
        assert PY_BAD_SPEC in status.toolTip()

    def test_python_compatible_is_not_installed(self, installer):
        self._populate(
            installer, [_remote_entry("Py Ok", python_spec=PY_OK_SPEC)]
        )
        row = self._row_of(installer, "Py Ok")
        status = installer.table.item(row, 4)
        assert status.text() == "Not Installed"
        assert "Python" not in status.toolTip()

    def test_missing_python_spec_is_compatible(self, installer):
        self._populate(installer, [_remote_entry("No Spec")])
        row = self._row_of(installer, "No Spec")
        assert installer.table.item(row, 4).text() == "Not Installed"

    def test_both_incompatible_tooltip_has_both_lines(self, installer):
        self._populate(
            installer,
            [
                _remote_entry(
                    "Both Bad",
                    moleditpy_spec=APP_BAD_SPEC,
                    python_spec=PY_BAD_SPEC,
                )
            ],
        )
        row = self._row_of(installer, "Both Bad")
        tip = installer.table.item(row, 4).toolTip()
        assert "Requires MoleditPy" in tip
        assert "Requires Python" in tip

    def test_install_button_tooltip_mentions_python(self, installer):
        self._populate(
            installer, [_remote_entry("Py Bad", python_spec=PY_BAD_SPEC)]
        )
        row = self._row_of(installer, "Py Bad")
        btn = installer.table.cellWidget(row, 5)
        assert btn is not None
        assert "Python" in btn.toolTip()
        assert PY_BAD_SPEC in btn.toolTip()

    def test_install_button_no_python_tooltip_when_compatible(self, installer):
        self._populate(
            installer, [_remote_entry("Py Ok", python_spec=PY_OK_SPEC)]
        )
        row = self._row_of(installer, "Py Ok")
        btn = installer.table.cellWidget(row, 5)
        assert btn is not None
        assert "Python" not in btn.toolTip()

    # -- install-time warning path ------------------------------------------

    def _click_install(self, installer, name):
        """Fire on_update_clicked via a real button so self.sender() works.

        download_url is left unset: after the compatibility gates the handler
        aborts at 'Missing download URL', so no network is ever touched.
        """
        from PyQt6.QtWidgets import QPushButton

        btn = QPushButton(parent=installer)
        btn.setProperty("plugin_name", name)
        btn.setProperty("target_file", None)
        btn.setProperty("dependencies", [])
        btn.clicked.connect(installer.on_update_clicked)
        btn.click()

    def test_install_warning_no_aborts(self, installer, qapp):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        installer.remote_data = [
            _remote_entry("Py Bad", python_spec=PY_BAD_SPEC)
        ]
        with patch.object(
            _plugin_installer.QMessageBox,
            "warning",
            return_value=QMessageBox.StandardButton.No,
        ) as warn:
            self._click_install(installer, "Py Bad")
        assert warn.call_count == 1
        assert "Incompatible Python Version Warning" in warn.call_args[0]

    def test_install_warning_yes_proceeds_past_check(self, installer, qapp):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        installer.remote_data = [
            _remote_entry("Py Bad", python_spec=PY_BAD_SPEC)
        ]
        with patch.object(
            _plugin_installer.QMessageBox,
            "warning",
            return_value=QMessageBox.StandardButton.Yes,
        ) as warn:
            self._click_install(installer, "Py Bad")
        # Warning #1: python check (answered Yes) → proceeds; warning #2:
        # the missing-download-URL abort, proving the gate was passed.
        assert warn.call_count == 2
        assert "Incompatible Python Version Warning" in warn.call_args_list[0][0]
        assert any(
            "Missing download URL" in str(a) for a in warn.call_args_list[1][0]
        )

    def test_install_compatible_no_python_warning(self, installer, qapp):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        installer.remote_data = [
            _remote_entry("Py Ok", python_spec=PY_OK_SPEC)
        ]
        with patch.object(
            _plugin_installer.QMessageBox,
            "warning",
            return_value=QMessageBox.StandardButton.No,
        ) as warn:
            self._click_install(installer, "Py Ok")
        # Only the missing-download-URL abort fires — no compatibility warning.
        assert warn.call_count == 1
        assert any("Missing download URL" in str(a) for a in warn.call_args[0])


class TestDetailsDialogSupportedPython:
    """PluginDetailsDialog shows the Supported Python line."""

    def _find_labels(self, dlg):
        from PyQt6.QtWidgets import QLabel

        return [w.text() for w in dlg.findChildren(QLabel)]

    def test_shows_passed_spec(self, qapp):
        d = _plugin_installer.PluginDetailsDialog(
            None, "X", "A", "1.0", "Desc", [], None, None,
            ">=4.0.0", PY_OK_SPEC,
        )
        import html as _html

        texts = self._find_labels(d)
        escaped = _html.escape(PY_OK_SPEC)
        assert any(
            "Supported Python:" in t and escaped in t for t in texts
        )
        d.destroy()

    def test_defaults_to_unknown(self, qapp):
        d = _make_details_dialog()
        texts = self._find_labels(d)
        assert any(
            "Supported Python:" in t and "Unknown" in t for t in texts
        )
        d.destroy()


# ===========================================================================
# Fetch-result handling, downloads, and the install pipeline
# ===========================================================================


@pytest.fixture
def bare_installer(qapp):
    """PluginInstallerWindow with a mock main_window and no network."""
    from unittest.mock import MagicMock, patch

    mw = MagicMock()
    mw.plugin_manager.plugins = []
    with patch.object(
        _plugin_installer.PluginInstallerWindow,
        "check_updates",
        lambda self: None,
    ):
        d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
    d.main_window = mw
    d.get_app_version = lambda: "4.2.0"
    d.latest_app_version = "Checking..."
    yield d
    d.destroy()


class TestSetUiBusy:
    def test_busy_disables_refresh_and_update_all(self, bare_installer):
        bare_installer._set_ui_busy(True)
        assert not bare_installer._btn_refresh.isEnabled()
        assert not bare_installer.btn_update_all.isEnabled()

    def test_unbusy_reenables_refresh(self, bare_installer):
        bare_installer._set_ui_busy(True)
        bare_installer._set_ui_busy(False)
        assert bare_installer._btn_refresh.isEnabled()


class TestOnFetchDone:
    def test_update_available(self, bare_installer):
        bare_installer._on_fetch_done("9.9.9", [])
        assert "Update Available" in bare_installer.lbl_latest_version.text()
        assert bare_installer.latest_app_version == "9.9.9"

    def test_up_to_date(self, bare_installer):
        bare_installer._on_fetch_done("4.2.0", [])
        assert "Up to date" in bare_installer.lbl_latest_version.text()

    def test_dev_newer(self, bare_installer):
        bare_installer._on_fetch_done("4.0.0", [])
        assert "Dev/Newer" in bare_installer.lbl_latest_version.text()

    def test_unknown(self, bare_installer):
        bare_installer._on_fetch_done("Unknown", [])
        assert "Unknown" in bare_installer.lbl_latest_version.text()

    def test_remote_data_populates_table(self, bare_installer):
        bare_installer._on_fetch_done("Unknown", [_remote_entry("Fetched Plugin")])
        assert bare_installer.table.rowCount() == 1
        assert bare_installer.table.item(0, 0).text() == "Fetched Plugin"

    def test_upgrade_button_configured_when_newer_on_pypi(self, bare_installer):
        bare_installer._on_fetch_done("9.9.9", [])
        assert "9.9.9" in bare_installer.btn_upgrade_app.text()


class TestCheckUpdatesSilent:
    def test_newer_pypi_sets_updates_found_and_title(self, bare_installer):
        from unittest.mock import patch

        with patch.object(
            _plugin_installer._FetchWorker, "_fetch_pypi", lambda self: "9.9.9"
        ), patch.object(
            _plugin_installer._FetchWorker, "_fetch_remote", lambda self: []
        ):
            bare_installer.check_updates_silent()
        assert bare_installer.updates_found is True
        assert "Updates Available" in bare_installer.windowTitle()

    def test_same_version_no_updates(self, bare_installer):
        from unittest.mock import patch

        with patch.object(
            _plugin_installer._FetchWorker, "_fetch_pypi", lambda self: "4.2.0"
        ), patch.object(
            _plugin_installer._FetchWorker, "_fetch_remote", lambda self: []
        ):
            bare_installer.check_updates_silent()
        assert bare_installer.updates_found is False
        assert bare_installer.windowTitle() == "Plugin Installer"


class _FakeResponse:
    """Context-manager stand-in for urllib.request.urlopen()."""

    def __init__(self, payload):
        self._chunks = [payload[i : i + 5] for i in range(0, len(payload), 5)]
        self.headers = {"Content-Length": str(len(payload))}

    def read(self, _size):
        return self._chunks.pop(0) if self._chunks else b""

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class TestDownloadChunked:
    def test_success_writes_file_and_reports_progress(
        self, bare_installer, tmp_path, monkeypatch
    ):
        payload = b"0123456789abc"
        monkeypatch.setattr(
            _plugin_installer.urllib.request,
            "urlopen",
            lambda url, timeout=15: _FakeResponse(payload),
        )
        dest = tmp_path / "out.py"
        seen = []
        ok = bare_installer._download_chunked(
            "https://example.com/x.py",
            str(dest),
            lambda recv, total: seen.append((recv, total)) or True,
        )
        assert ok is True
        assert dest.read_bytes() == payload
        assert seen[-1] == (len(payload), len(payload))

    def test_progress_false_aborts(self, bare_installer, tmp_path, monkeypatch):
        monkeypatch.setattr(
            _plugin_installer.urllib.request,
            "urlopen",
            lambda url, timeout=15: _FakeResponse(b"0123456789"),
        )
        ok = bare_installer._download_chunked(
            "https://example.com/x.py",
            str(tmp_path / "out.py"),
            lambda recv, total: False,
        )
        assert ok is False

    def test_network_error_retries_then_fails(
        self, bare_installer, tmp_path, monkeypatch
    ):
        from unittest.mock import MagicMock, patch

        calls = []

        def _raise(url, timeout=15):
            calls.append(url)
            raise _plugin_installer.urllib.error.URLError("boom")

        monkeypatch.setattr(_plugin_installer.urllib.request, "urlopen", _raise)
        monkeypatch.setattr(_plugin_installer.time, "sleep", lambda s: None)
        with patch.object(
            _plugin_installer.QMessageBox, "warning", MagicMock()
        ) as warn:
            ok = bare_installer._download_chunked(
                "https://example.com/x.py", str(tmp_path / "out.py")
            )
        assert ok is False
        assert len(calls) == 3
        assert warn.call_count == 1
        assert "after 3 attempts" in str(warn.call_args)

    def test_unexpected_error_fails_without_retry(
        self, bare_installer, tmp_path, monkeypatch
    ):
        from unittest.mock import MagicMock, patch

        calls = []

        def _raise(url, timeout=15):
            calls.append(url)
            raise ValueError("bad")

        monkeypatch.setattr(_plugin_installer.urllib.request, "urlopen", _raise)
        with patch.object(
            _plugin_installer.QMessageBox, "warning", MagicMock()
        ) as warn:
            ok = bare_installer._download_chunked(
                "https://example.com/x.py", str(tmp_path / "out.py")
            )
        assert ok is False
        assert len(calls) == 1
        assert warn.call_count == 1


class TestOverwriteFolderPlugin:
    def test_wipe_preserves_settings_and_drops_orphans(self, tmp_path):
        source = tmp_path / "src"
        target = tmp_path / "dst"
        source.mkdir()
        target.mkdir()
        (source / "__init__.py").write_text("new code")
        (target / "__init__.py").write_text("old code")
        (target / "orphan.txt").write_text("stale")
        (target / "settings.json").write_text('{"user": 1}')

        _plugin_installer.PluginInstallerWindow._overwrite_folder_plugin(
            str(source), str(target)
        )

        assert (target / "__init__.py").read_text() == "new code"
        assert not (target / "orphan.txt").exists()
        assert (target / "settings.json").read_text() == '{"user": 1}'


def _install_button(installer, *, name, download_url=None, target_file=None):
    from PyQt6.QtWidgets import QPushButton

    btn = QPushButton(parent=installer)
    btn.setProperty("plugin_name", name)
    btn.setProperty("download_url", download_url)
    btn.setProperty("target_file", target_file)
    btn.setProperty("dependencies", [])
    btn.clicked.connect(installer.on_update_clicked)
    return btn


class TestInstallPipeline:
    """on_update_clicked end-to-end with the download stubbed out."""

    PAYLOAD = b"PLUGIN_VERSION = '2.0.0'\n"

    def _entry(self, name, sha=True):
        import hashlib

        entry = _remote_entry(name)
        entry["downloadUrl"] = "https://example.com/fake_plugin.py"
        if sha:
            entry["sha256"] = hashlib.sha256(self.PAYLOAD).hexdigest()
        return entry

    def _stub_download(self, installer):
        data = self.PAYLOAD

        def _fake(url, dest, cb=None):
            with open(dest, "wb") as f:
                f.write(data)
            return True

        installer._download_chunked = _fake

    def test_successful_update_overwrites_target(self, bare_installer, tmp_path):
        from unittest.mock import MagicMock, patch
        from PyQt6.QtWidgets import QMessageBox

        target = tmp_path / "fake_plugin.py"
        target.write_bytes(b"OLD")
        bare_installer.remote_data = [self._entry("Fake Plugin")]
        self._stub_download(bare_installer)

        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.Yes,
        ), patch.object(
            _plugin_installer.QMessageBox, "information", MagicMock()
        ) as info:
            _install_button(
                bare_installer,
                name="Fake Plugin",
                download_url="https://example.com/fake_plugin.py",
                target_file=str(target),
            ).click()

        assert target.read_bytes() == self.PAYLOAD
        assert bare_installer._last_install_succeeded is True
        assert "Fake Plugin" in bare_installer._pending_installs
        assert bare_installer._needs_plugin_reload is True
        assert info.call_count == 1
        assert "Successfully updated" in str(info.call_args)

    def test_sha_mismatch_blocks_install(self, bare_installer, tmp_path):
        from unittest.mock import MagicMock, patch
        from PyQt6.QtWidgets import QMessageBox

        target = tmp_path / "fake_plugin.py"
        target.write_bytes(b"OLD")
        entry = self._entry("Fake Plugin")
        entry["sha256"] = "0" * 64
        bare_installer.remote_data = [entry]
        self._stub_download(bare_installer)

        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.Yes,
        ), patch.object(
            _plugin_installer.QMessageBox, "critical", MagicMock()
        ) as crit:
            _install_button(
                bare_installer,
                name="Fake Plugin",
                download_url="https://example.com/fake_plugin.py",
                target_file=str(target),
            ).click()

        assert crit.call_count == 1
        assert "SHA256 Checksum Mismatch" in str(crit.call_args)
        assert target.read_bytes() == b"OLD"
        assert bare_installer._last_install_succeeded is False
        assert "Fake Plugin" not in bare_installer._pending_installs

    def test_missing_sha_prompt_declined_aborts(self, bare_installer, tmp_path):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        target = tmp_path / "fake_plugin.py"
        target.write_bytes(b"OLD")
        bare_installer.remote_data = [self._entry("Fake Plugin", sha=False)]
        self._stub_download(bare_installer)

        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.Yes,
        ), patch.object(
            _plugin_installer.QMessageBox,
            "warning",
            return_value=QMessageBox.StandardButton.No,
        ) as warn:
            _install_button(
                bare_installer,
                name="Fake Plugin",
                download_url="https://example.com/fake_plugin.py",
                target_file=str(target),
            ).click()

        assert warn.call_count == 1
        assert "SHA256 checksum is not available" in str(warn.call_args)
        assert target.read_bytes() == b"OLD"
        assert bare_installer._last_install_succeeded is False

    def test_unsafe_scheme_rejected(self, bare_installer):
        from unittest.mock import MagicMock, patch
        from PyQt6.QtWidgets import QMessageBox

        bare_installer.remote_data = [self._entry("Fake Plugin")]

        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.Yes,
        ), patch.object(
            _plugin_installer.QMessageBox, "warning", MagicMock()
        ) as warn:
            _install_button(
                bare_installer,
                name="Fake Plugin",
                download_url="ftp://example.com/fake_plugin.py",
            ).click()

        assert warn.call_count == 1
        assert "unsafe scheme" in str(warn.call_args)
        assert bare_installer._last_install_succeeded is False


class TestUpdateAllPlugins:
    def test_no_updates_shows_info(self, bare_installer):
        from unittest.mock import MagicMock, patch

        bare_installer.remote_data = [_remote_entry("Some Plugin")]
        bare_installer.populate_table(silent=True)
        with patch.object(
            _plugin_installer.QMessageBox, "information", MagicMock()
        ) as info:
            bare_installer.update_all_plugins()
        assert info.call_count == 1
        assert "No plugins with available updates" in str(info.call_args)


class TestOnFinished:
    def test_reload_runs_once_and_clears_state(self, bare_installer):
        from unittest.mock import MagicMock, patch

        bare_installer._needs_plugin_reload = True
        bare_installer._pending_installs["X"] = {"name": "X"}
        with patch.object(
            _plugin_installer, "_refresh_plugin_menus", MagicMock()
        ) as refresh:
            bare_installer._on_finished(0)
        bare_installer.main_window.plugin_manager.discover_plugins.assert_called_once()
        refresh.assert_called_once()
        assert bare_installer._pending_installs == {}
        assert bare_installer._needs_plugin_reload is False

    def test_no_reload_when_nothing_installed(self, bare_installer):
        bare_installer._needs_plugin_reload = False
        bare_installer._on_finished(0)
        bare_installer.main_window.plugin_manager.discover_plugins.assert_not_called()


class TestCopyUpgradeCommand:
    def test_clipboard_receives_pip_command(self, bare_installer):
        from unittest.mock import MagicMock, patch
        from PyQt6.QtWidgets import QApplication

        with patch.object(
            _plugin_installer.QMessageBox, "information", MagicMock()
        ):
            bare_installer.copy_upgrade_command()
        assert QApplication.clipboard().text().startswith("pip install --upgrade ")


class TestPerformStartupCheck:
    def _fake_checker(self, updates_found):
        from unittest.mock import MagicMock

        checker = MagicMock()
        checker.updates_found = updates_found
        return checker

    def test_no_updates_no_prompt(self, qapp, monkeypatch):
        from unittest.mock import MagicMock, patch

        checker = self._fake_checker(False)
        monkeypatch.setattr(
            _plugin_installer,
            "PluginInstallerWindow",
            lambda mw, auto_check: checker,
        )
        with patch.object(
            _plugin_installer.QMessageBox, "question", MagicMock()
        ) as q:
            _plugin_installer.perform_startup_check(None)
        q.assert_not_called()
        checker.deleteLater.assert_called_once()

    def test_updates_declined_deletes_checker(self, qapp, monkeypatch):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        checker = self._fake_checker(True)
        monkeypatch.setattr(
            _plugin_installer,
            "PluginInstallerWindow",
            lambda mw, auto_check: checker,
        )
        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.No,
        ):
            _plugin_installer.perform_startup_check(None)
        checker.show.assert_not_called()
        checker.deleteLater.assert_called_once()

    def test_updates_accepted_shows_window(self, qapp, monkeypatch):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        checker = self._fake_checker(True)
        monkeypatch.setattr(
            _plugin_installer,
            "PluginInstallerWindow",
            lambda mw, auto_check: checker,
        )
        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.Yes,
        ):
            _plugin_installer.perform_startup_check(None)
        checker.show.assert_called_once()
        checker.raise_.assert_called_once()
        checker.deleteLater.assert_not_called()


class TestShowPluginDetails:
    """show_plugin_details builds the dialog from remote + local info."""

    @pytest.fixture
    def captured_dialog(self, monkeypatch):
        calls = []

        class _FakeDialog:
            def __init__(self, *args):
                calls.append(args)

            def exec(self):
                return 0

        monkeypatch.setattr(
            _plugin_installer, "PluginDetailsDialog", _FakeDialog
        )
        return calls

    def test_remote_only_plugin(self, bare_installer, captured_dialog):
        entry = _remote_entry("Remote Only", python_spec=PY_OK_SPEC)
        entry["dependencies"] = ["numpy"]
        bare_installer.remote_data = [entry]
        bare_installer.populate_table(silent=True)

        bare_installer.show_plugin_details(0, 0)

        assert len(captured_dialog) == 1
        args = captured_dialog[0]
        # (parent, name, author, version, description, deps, local, target, app_spec, py_spec)
        assert args[1] == "Remote Only"
        assert args[3] == "1.0.0"
        assert args[5] == ["numpy"]
        assert args[8] == APP_OK_SPEC
        assert args[9] == PY_OK_SPEC

    def test_installed_plugin_merges_versions(
        self, bare_installer, captured_dialog, tmp_path
    ):
        target = tmp_path / "remote_only.py"
        target.write_text("PLUGIN_VERSION = '0.9.0'\n")
        bare_installer.remote_data = [_remote_entry("Remote Only")]
        bare_installer.main_window.plugin_manager.plugins = [
            {
                "name": "Remote Only",
                "version": "0.9.0",
                "filepath": str(target),
                "author": "Local Author",
            }
        ]
        bare_installer.populate_table(silent=True)

        bare_installer.show_plugin_details(0, 0)

        args = captured_dialog[0]
        assert "Installed: 0.9.0" in args[3]
        assert "Latest: 1.0.0" in args[3]
        assert args[7] == str(target)

    def test_not_in_registry_uses_unknown_specs(
        self, bare_installer, captured_dialog, tmp_path
    ):
        target = tmp_path / "local_only.py"
        target.write_text("PLUGIN_VERSION = '0.5.0'\n")
        bare_installer.remote_data = []
        bare_installer.main_window.plugin_manager.plugins = [
            {
                "name": "Local Only",
                "version": "0.5.0",
                "filepath": str(target),
                "author": "Local Author",
            }
        ]
        bare_installer.populate_table(silent=True)

        bare_installer.show_plugin_details(0, 0)

        args = captured_dialog[0]
        assert args[3] == "Installed: 0.5.0"
        assert args[8] == "Unknown"
        assert args[9] == "Unknown"

    def test_empty_row_is_noop(self, bare_installer, captured_dialog):
        bare_installer.show_plugin_details(99, 0)
        assert captured_dialog == []


class TestBatchUpdateAll:
    PAYLOAD = b"PLUGIN_VERSION = '1.0.0'\n"

    def test_updates_one_plugin_end_to_end(self, bare_installer, tmp_path):
        import hashlib
        from unittest.mock import MagicMock, patch
        from PyQt6.QtWidgets import QMessageBox

        target = tmp_path / "fake_plugin.py"
        target.write_text("PLUGIN_VERSION = '0.9.0'\n")
        entry = _remote_entry("Fake Plugin")
        entry["downloadUrl"] = "https://example.com/fake_plugin.py"
        entry["sha256"] = hashlib.sha256(self.PAYLOAD).hexdigest()
        bare_installer.remote_data = [entry]
        bare_installer.main_window.plugin_manager.plugins = [
            {
                "name": "Fake Plugin",
                "version": "0.9.0",
                "filepath": str(target),
                "author": "Test Author",
            }
        ]
        bare_installer.populate_table(silent=True)
        assert bare_installer.table.item(0, 4).text() == "Update Available"

        data = self.PAYLOAD

        def _fake_download(url, dest, cb=None):
            with open(dest, "wb") as f:
                f.write(data)
            return True

        bare_installer._download_chunked = _fake_download

        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.Yes,
        ), patch.object(
            _plugin_installer.QMessageBox, "information", MagicMock()
        ) as info:
            bare_installer.update_all_plugins()

        assert target.read_bytes() == self.PAYLOAD
        assert "Fake Plugin" in bare_installer._pending_installs
        assert bare_installer._batch_updating is False
        assert bare_installer._batch_progress is None
        # Only the final summary dialog — per-plugin success is suppressed in batch
        assert info.call_count == 1
        assert "Updated 1 of 1" in str(info.call_args)

    def test_batch_confirm_declined_is_noop(self, bare_installer, tmp_path):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        target = tmp_path / "fake_plugin.py"
        target.write_text("PLUGIN_VERSION = '0.9.0'\n")
        entry = _remote_entry("Fake Plugin")
        bare_installer.remote_data = [entry]
        bare_installer.main_window.plugin_manager.plugins = [
            {
                "name": "Fake Plugin",
                "version": "0.9.0",
                "filepath": str(target),
                "author": "Test Author",
            }
        ]
        bare_installer.populate_table(silent=True)

        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.No,
        ):
            bare_installer.update_all_plugins()

        assert bare_installer._pending_installs == {}
        assert target.read_text() == "PLUGIN_VERSION = '0.9.0'\n"


class TestFreshInstallViaManager:
    PAYLOAD = b"PLUGIN_VERSION = '1.0.0'\n"

    def test_install_delegates_to_plugin_manager(self, bare_installer):
        import hashlib
        from unittest.mock import MagicMock, patch
        from PyQt6.QtWidgets import QMessageBox

        entry = _remote_entry("Fresh Plugin")
        entry["downloadUrl"] = "https://example.com/fresh_plugin.py"
        entry["sha256"] = hashlib.sha256(self.PAYLOAD).hexdigest()
        bare_installer.remote_data = [entry]

        data = self.PAYLOAD

        def _fake_download(url, dest, cb=None):
            with open(dest, "wb") as f:
                f.write(data)
            return True

        bare_installer._download_chunked = _fake_download

        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.Yes,
        ), patch.object(
            _plugin_installer.QMessageBox, "information", MagicMock()
        ) as info:
            _install_button(
                bare_installer,
                name="Fresh Plugin",
                download_url="https://example.com/fresh_plugin.py",
            ).click()

        install = bare_installer.main_window.plugin_manager.install_plugin
        install.assert_called_once()
        assert install.call_args[0][0].endswith("fresh_plugin.py")
        assert bare_installer._last_install_succeeded is True
        assert "Fresh Plugin" in bare_installer._pending_installs
        assert "Successfully installed" in str(info.call_args)


class TestZipFolderOverwrite:
    def _zip_bytes(self):
        import io
        import zipfile as zf

        buf = io.BytesIO()
        with zf.ZipFile(buf, "w") as z:
            z.writestr("My_Plugin/__init__.py", "new folder code")
            z.writestr("My_Plugin/extra.py", "helper")
        return buf.getvalue()

    def test_folder_plugin_updated_from_zip(self, bare_installer, tmp_path):
        import hashlib
        from unittest.mock import MagicMock, patch
        from PyQt6.QtWidgets import QMessageBox

        plugin_dir = tmp_path / "My_Plugin"
        plugin_dir.mkdir()
        target = plugin_dir / "__init__.py"
        target.write_text("old folder code")
        (plugin_dir / "orphan.txt").write_text("stale")
        (plugin_dir / "settings.json").write_text('{"keep": true}')

        payload = self._zip_bytes()
        entry = _remote_entry("Folder Plugin")
        entry["downloadUrl"] = "https://example.com/My_Plugin.zip"
        entry["sha256"] = hashlib.sha256(payload).hexdigest()
        bare_installer.remote_data = [entry]

        def _fake_download(url, dest, cb=None):
            with open(dest, "wb") as f:
                f.write(payload)
            return True

        bare_installer._download_chunked = _fake_download

        with patch.object(
            _plugin_installer.QMessageBox,
            "question",
            return_value=QMessageBox.StandardButton.Yes,
        ), patch.object(
            _plugin_installer.QMessageBox, "information", MagicMock()
        ):
            _install_button(
                bare_installer,
                name="Folder Plugin",
                download_url="https://example.com/My_Plugin.zip",
                target_file=str(target),
            ).click()

        assert target.read_text() == "new folder code"
        assert (plugin_dir / "extra.py").read_text() == "helper"
        assert not (plugin_dir / "orphan.txt").exists()
        assert (plugin_dir / "settings.json").read_text() == '{"keep": true}'
        # Manual overwrite path — the plugin manager was never involved
        bare_installer.main_window.plugin_manager.install_plugin.assert_not_called()
        assert bare_installer._last_install_succeeded is True
