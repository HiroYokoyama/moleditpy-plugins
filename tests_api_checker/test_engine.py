import importlib.util
import json
import tempfile
import textwrap
import unittest
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_TESTS_DIR = Path(__file__).resolve().parent
_PLUGIN_ROOT = _TESTS_DIR.parent
_WORKSPACE_ROOT = _PLUGIN_ROOT.parent

_DEFAULT_APP_CANDIDATES = [
    _WORKSPACE_ROOT / "python_molecular_editor",
    _PLUGIN_ROOT / "python_molecular_editor",  # CI path
]
_APP_PATH = next(
    (p for p in _DEFAULT_APP_CANDIDATES if p and (p / "moleditpy").exists()), None
)
HAS_APP = _APP_PATH is not None


# ---------------------------------------------------------------------------
# Load plugin_api_checker (always present in the repo)
# ---------------------------------------------------------------------------


def _load_checker():
    checker_path = _TESTS_DIR.parent / "api-checker/plugin_api_checker.py"
    spec = importlib.util.spec_from_file_location("plugin_api_checker", checker_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    mod.__file__ = str(checker_path)
    spec.loader.exec_module(mod)
    return mod


_checker = _load_checker()
Issue = _checker.Issue
APIInfo = _checker.APIInfo
AppAPIExtractor = _checker.AppAPIExtractor
PluginFileChecker = _checker.PluginFileChecker
_merge_allowlists = _checker._merge_allowlists
_load_site_allowlist = _checker._load_site_allowlist


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_api(
    mw_members=None,
    manager_members=None,
    manager_class_names=None,
    context_members=None,
):
    api = APIInfo()
    if mw_members:
        api.mw_members.update(mw_members)
    if manager_members:
        api.manager_members.update(manager_members)
    if manager_class_names:
        api.manager_class_names.update(manager_class_names)
    if context_members:
        api.context_members.update(context_members)
    return api


def _check_source(
    source: str, api: APIInfo, check_context: bool = False, allowlist=None
) -> list:
    """Run PluginFileChecker on a string of source via a temp file."""
    with tempfile.NamedTemporaryFile(
        suffix=".py", mode="w", encoding="utf-8", delete=False
    ) as f:
        f.write(textwrap.dedent(source))
        tmp = Path(f.name)
    try:
        checker = PluginFileChecker(
            tmp, api, check_context=check_context, allowlist=allowlist or {}
        )
        return checker.check()
    finally:
        tmp.unlink(missing_ok=True)


def _make_synthetic_app(tmp_dir: Path) -> None:
    """Write a minimal main-app tree so AppAPIExtractor has something to parse."""
    (tmp_dir / "main_window.py").write_text(
        textwrap.dedent("""\
        class MainWindow:
            mol_loaded = None

            def __init__(self):
                self.state_manager = StateManager(self)

            def get_molecule(self):
                pass

            @property
            def current_mol(self):
                pass
    """),
        encoding="utf-8",
    )

    (tmp_dir / "state_manager.py").write_text(
        textwrap.dedent("""\
        class StateManager:
            def __init__(self, host):
                self.host = host
                self.mol_data = None

            def load_mol(self, mol):
                pass
    """),
        encoding="utf-8",
    )

    (tmp_dir / "plugin_interface.py").write_text(
        textwrap.dedent("""\
        class PluginContext:
            def get_main_window(self): pass
            def register_file_opener(self, ext, cb, priority=0): pass
            def register_drop_handler(self, cb, priority=0): pass
            def add_menu_action(self, path, cb): pass
            def register_window(self, key, win): pass
            def get_window(self, key): pass
            def show_status_message(self, msg, duration=0): pass
    """),
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# TestIssue
# ---------------------------------------------------------------------------


class TestIssue(unittest.TestCase):
    def test_str_in_try_shows_tag(self):
        i = Issue("f.py", 10, "CODE", "msg", in_try=True)
        self.assertIn("[try]", str(i))

    def test_str_not_in_try_no_tag(self):
        i = Issue("f.py", 10, "CODE", "msg", in_try=False)
        self.assertNotIn("[try]", str(i))

    def test_key_returns_four_tuple(self):
        i = Issue("f.py", 5, "UNKNOWN_MW_ATTR", "bad attr")
        self.assertEqual(i.key(), ("f.py", 5, "UNKNOWN_MW_ATTR", "bad attr"))

    def test_key_excludes_in_try_flag(self):
        a = Issue("f.py", 1, "X", "y", in_try=True)
        b = Issue("f.py", 1, "X", "y", in_try=False)
        self.assertEqual(a.key(), b.key())


# ---------------------------------------------------------------------------
# TestMergeAllowlists
# ---------------------------------------------------------------------------


class TestMergeAllowlists(unittest.TestCase):
    def test_mw_sets_are_unioned(self):
        result = _merge_allowlists({"mw": {"attr1"}}, {"mw": {"attr2"}})
        self.assertEqual(result["mw"], {"attr1", "attr2"})

    def test_manager_members_are_unioned(self):
        a = {"manager": {"state_manager": {"data"}}}
        b = {"manager": {"state_manager": {"mol"}}}
        result = _merge_allowlists(a, b)
        self.assertEqual(result["manager"]["state_manager"], {"data", "mol"})

    def test_context_sets_are_unioned(self):
        a = {"context": {"get_main_window"}}
        b = {"context": {"register_file_opener"}}
        result = _merge_allowlists(a, b)
        self.assertIn("get_main_window", result["context"])
        self.assertIn("register_file_opener", result["context"])

    def test_empty_dicts_return_empty(self):
        self.assertEqual(_merge_allowlists({}, {}), {})

    def test_single_dict_is_preserved(self):
        result = _merge_allowlists({"mw": {"foo", "bar"}})
        self.assertEqual(result["mw"], {"foo", "bar"})


# ---------------------------------------------------------------------------
# TestLoadSiteAllowlist
# ---------------------------------------------------------------------------


class TestLoadSiteAllowlist(unittest.TestCase):
    def _write(self, tmp_dir: Path, data: dict):
        (tmp_dir / ".moleditpy-api-allowlist").write_text(
            json.dumps(data), encoding="utf-8"
        )

    def test_no_file_returns_empty(self):
        with tempfile.TemporaryDirectory() as d:
            result = _load_site_allowlist(Path(d))
        self.assertEqual(result, {})

    def test_mw_list_form(self):
        with tempfile.TemporaryDirectory() as d:
            self._write(Path(d), {"mw": ["foo", "bar"]})
            result = _load_site_allowlist(Path(d))
        self.assertEqual(result.get("mw"), {"foo", "bar"})

    def test_mw_dict_form(self):
        with tempfile.TemporaryDirectory() as d:
            self._write(Path(d), {"mw": {"foo": "legacy compat bridge"}})
            result = _load_site_allowlist(Path(d))
        self.assertIn("foo", result.get("mw", set()))

    def test_manager_list_form(self):
        with tempfile.TemporaryDirectory() as d:
            self._write(Path(d), {"manager": {"state_manager": ["data"]}})
            result = _load_site_allowlist(Path(d))
        self.assertIn("data", result.get("manager", {}).get("state_manager", set()))

    def test_context_list_form(self):
        with tempfile.TemporaryDirectory() as d:
            self._write(Path(d), {"context": ["register_file_opener"]})
            result = _load_site_allowlist(Path(d))
        self.assertIn("register_file_opener", result.get("context", set()))

    def test_invalid_json_returns_empty(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / ".moleditpy-api-allowlist").write_text(
                "not json", encoding="utf-8"
            )
            result = _load_site_allowlist(Path(d))
        self.assertEqual(result, {})


# ---------------------------------------------------------------------------
# TestPluginFileChecker  (synthetic code — no main app required)
# ---------------------------------------------------------------------------


class TestPluginFileChecker(unittest.TestCase):
    def setUp(self):
        self.api = _make_api(
            mw_members={
                "get_molecule": "method",
                "current_mol": "property",
                "state_manager": "manager:StateManager",
            },
            manager_members={"state_manager": {"load_mol", "mol_data"}},
            manager_class_names={"state_manager": "StateManager"},
            context_members={
                "get_main_window",
                "register_file_opener",
                "register_drop_handler",
                "add_menu_action",
                "register_window",
                "get_window",
                "show_status_message",
            },
        )

    def _issues(self, source, check_context=False, allowlist=None):
        return _check_source(
            source, self.api, check_context=check_context, allowlist=allowlist
        )

    # --- basic clean code ---

    def test_no_issues_for_clean_code(self):
        self.assertEqual(self._issues("x = 1 + 2\nprint(x)\n"), [])

    # --- unknown MW attr ---

    def test_unknown_mw_attr_flagged(self):
        issues = self._issues("def f(mw):\n    mw.nonexistent()\n")
        self.assertTrue(any(i.code == "UNKNOWN_MW_ATTR" for i in issues))

    def test_known_mw_attr_not_flagged(self):
        self.assertEqual(self._issues("def f(mw):\n    mw.get_molecule()\n"), [])

    # --- private and Qt-inherited skips ---

    def test_private_attr_not_flagged(self):
        self.assertEqual(self._issues("def f(mw):\n    mw._internal()\n"), [])

    def test_qt_inherited_show_not_flagged(self):
        self.assertEqual(self._issues("def f(mw):\n    mw.show()\n"), [])

    def test_qt_inherited_raise_not_flagged(self):
        self.assertEqual(self._issues("def f(mw):\n    mw.raise_()\n"), [])

    # --- safe call suppression ---

    def test_hasattr_not_flagged(self):
        source = "def f(mw):\n    if hasattr(mw, 'x'):\n        pass\n"
        self.assertEqual(self._issues(source), [])

    # --- try block flag ---

    def test_in_try_flag_set(self):
        source = textwrap.dedent("""\
            def f(mw):
                try:
                    mw.nonexistent()
                except Exception:
                    pass
        """)
        issues = self._issues(source)
        self.assertTrue(any(i.code == "UNKNOWN_MW_ATTR" and i.in_try for i in issues))

    # --- alias tracking ---

    def test_get_main_window_alias_tracked(self):
        source = textwrap.dedent("""\
            def f(context):
                w = context.get_main_window()
                w.nonexistent()
        """)
        issues = self._issues(source)
        self.assertTrue(any(i.code == "UNKNOWN_MW_ATTR" for i in issues))

    def test_self_mw_attr_tracked(self):
        source = textwrap.dedent("""\
            class Plugin:
                def __init__(self, mw):
                    self.main_window = mw
                def run(self):
                    self.main_window.nonexistent()
        """)
        issues = self._issues(source)
        self.assertTrue(any(i.code == "UNKNOWN_MW_ATTR" for i in issues))

    # --- allowlist suppression ---

    def test_mw_allowlist_suppresses_issue(self):
        source = "def f(mw):\n    mw.legacy_attr()\n"
        self.assertEqual(self._issues(source, allowlist={"mw": {"legacy_attr"}}), [])

    def test_manager_allowlist_suppresses_issue(self):
        source = "def f(mw):\n    mw.state_manager.runtime_attr()\n"
        al = {"manager": {"state_manager": {"runtime_attr"}}}
        self.assertEqual(self._issues(source, allowlist=al), [])

    # --- manager attr checks ---

    def test_unknown_manager_attr_flagged(self):
        issues = self._issues("def f(mw):\n    mw.state_manager.bad_attr()\n")
        self.assertTrue(any(i.code == "UNKNOWN_MANAGER_ATTR" for i in issues))

    def test_known_manager_attr_not_flagged(self):
        self.assertEqual(
            self._issues("def f(mw):\n    mw.state_manager.load_mol(None)\n"), []
        )

    # --- context checks ---

    def test_unknown_context_attr_flagged_when_enabled(self):
        source = "def initialize(context):\n    context.nonexistent_ctx_method()\n"
        issues = self._issues(source, check_context=True)
        self.assertTrue(any(i.code == "UNKNOWN_CONTEXT_ATTR" for i in issues))

    def test_known_context_attr_not_flagged(self):
        source = "def initialize(context):\n    context.get_main_window()\n"
        self.assertEqual(self._issues(source, check_context=True), [])

    def test_context_not_checked_by_default(self):
        source = "def initialize(context):\n    context.nonexistent_ctx_method()\n"
        self.assertEqual(self._issues(source, check_context=False), [])

    # --- plugin-specific context surface ---

    def test_register_file_opener_is_known_context_attr(self):
        source = "def initialize(context):\n    context.register_file_opener('.out', lambda p: None, priority=100)\n"
        self.assertEqual(self._issues(source, check_context=True), [])

    def test_register_drop_handler_is_known_context_attr(self):
        source = "def initialize(context):\n    context.register_drop_handler(lambda p: False, priority=100)\n"
        self.assertEqual(self._issues(source, check_context=True), [])

    def test_register_window_is_known_context_attr(self):
        source = "def f(context, win):\n    context.register_window('analyzer', win)\n"
        self.assertEqual(self._issues(source, check_context=True), [])

    def test_get_window_is_known_context_attr(self):
        source = "def f(context):\n    context.get_window('analyzer')\n"
        self.assertEqual(self._issues(source, check_context=True), [])

    # --- edge cases ---

    def test_syntax_error_produces_syntax_issue(self):
        issues = self._issues("def broken(\n")
        self.assertTrue(any(i.code == "SYNTAX_ERROR" for i in issues))

    def test_dedup_same_access_not_triple_counted(self):
        source = "def f(mw):\n    mw.nonexistent()\n    mw.nonexistent()\n"
        issues = [i for i in self._issues(source) if i.code == "UNKNOWN_MW_ATTR"]
        self.assertLessEqual(len(issues), 2)


# ---------------------------------------------------------------------------
# TestAppAPIExtractor  (synthetic app tree in temp dir — no main app required)
# ---------------------------------------------------------------------------


class TestAppAPIExtractor(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp_dir = Path(self._tmp.name)
        _make_synthetic_app(self.tmp_dir)

    def tearDown(self):
        self._tmp.cleanup()

    def test_returns_api_info_instance(self):
        api = AppAPIExtractor(self.tmp_dir, verbose=False).extract()
        self.assertIsInstance(api, APIInfo)

    def test_extracts_mw_method(self):
        api = AppAPIExtractor(self.tmp_dir, verbose=False).extract()
        self.assertIn("get_molecule", api.mw_members)

    def test_extracts_mw_property(self):
        api = AppAPIExtractor(self.tmp_dir, verbose=False).extract()
        self.assertIn("current_mol", api.mw_members)

    def test_extracts_mw_class_attr(self):
        api = AppAPIExtractor(self.tmp_dir, verbose=False).extract()
        self.assertIn("mol_loaded", api.mw_members)

    def test_extracts_manager_as_mw_member(self):
        api = AppAPIExtractor(self.tmp_dir, verbose=False).extract()
        self.assertIn("state_manager", api.mw_members)

    def test_extracts_manager_members(self):
        api = AppAPIExtractor(self.tmp_dir, verbose=False).extract()
        self.assertIn("state_manager", api.manager_members)
        self.assertIn("load_mol", api.manager_members["state_manager"])

    def test_extracts_context_members(self):
        api = AppAPIExtractor(self.tmp_dir, verbose=False).extract()
        self.assertIn("get_main_window", api.context_members)
        self.assertIn("register_file_opener", api.context_members)

    def test_scans_host_assignments(self):
        (self.tmp_dir / "extra_manager.py").write_text(
            textwrap.dedent("""\
            class ExtraManager:
                def setup(self, host):
                    self.host = host
                    self.host.plugin_manager = self
        """),
            encoding="utf-8",
        )
        api = AppAPIExtractor(self.tmp_dir, verbose=False).extract()
        self.assertIn("plugin_manager", api.mw_members)


# ---------------------------------------------------------------------------
# TestAPIChecker — integration test (skipped when main app is not present)
# ---------------------------------------------------------------------------
