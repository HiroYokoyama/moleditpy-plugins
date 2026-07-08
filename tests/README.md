# moleditpy-plugins Test Suite

All tests run headlessly — no GUI, no chemistry libraries, no network required.
Heavy optional dependencies (PyQt6, rdkit, moleditpy, pyvista, …) are
intercepted by a custom `MetaPathFinder` in `conftest.py` and replaced with
`MagicMock` objects for the duration of each test.

---

## Running the tests

```bash
# Full suite
pytest tests/ -v

# Single file
pytest tests/test_initialize.py -v

# Single test by name
pytest tests/ -k "Atom Colorizer"
```

---

## Test files

### `conftest.py` — Shared infrastructure

Not a test file. Provides helpers used by all smoke/API test files:

| Symbol | Purpose |
|---|---|
| `BLOCKED_TOPS` | Set of top-level package names replaced with MagicMock |
| `mock_optional_imports()` | Context manager: installs the MetaPathFinder, cleans up on exit |
| `load_plugin(path)` | Load a plugin `.py` file as an isolated module (call inside `mock_optional_imports()`) |
| `make_context()` | Return a stub `PluginContext` (MagicMock with a non-None main window) |
| `visible_py_plugins(entry_point)` | Iterate visible single-file `.py` plugins from the registry, optionally filtered by entry-point function name |
| `mocks_with_real_numpy()` | Like `mock_optional_imports()` but keeps real numpy in `sys.modules` for generators/parsers doing real vector math |
| `P3`, `FakeAtom`, `FakeBond`, `FakeConf`, `FakeMol` | Shared fake rdkit bond-graph objects |
| `extract_function(path, class_name, fn_name, extra_globals=None)` | AST-based method/function extractor for methods on Qt-derived (mocked-base) classes |

---

### `test_registry.py` — Registry integrity (7 tests)

Validates `REGISTRY/plugins.json` without loading any plugin code.

| Test | What it checks |
|---|---|
| `test_registry_required_fields[field]` | Every visible plugin has all required fields: `name`, `version`, `description`, `downloadUrl`, `sha256` |
| `test_registry_names_unique` | No duplicate plugin names |
| `test_registry_sha256_match` | `sha256` in registry matches the actual file on disk |
| `test_local_files_exist` | Every local `downloadUrl` resolves to an existing file |
| `test_plugin_source_metadata` | Each `.py` plugin defines all required `PLUGIN_*` constants |
| `test_version_consistency` | `PLUGIN_VERSION` in source matches `version` in registry |

---

### `test_imports.py` — Static analysis (3 × N tests, N = visible plugins)

AST-based checks. No code is executed.

| Test | What it checks |
|---|---|
| `test_plugin_syntax[name]` | Plugin file parses without `SyntaxError` |
| `test_plugin_has_run_entrypoint[name]` | Plugin defines at least one of `run`, `autorun`, or `initialize` |
| `test_plugin_no_stdlib_import_errors[name]` | No top-level imports of missing stdlib modules (optional deps are skipped) |

---

### `test_initialize.py` — `initialize()` smoke tests (29 tests)

For every visible plugin that exposes `initialize(context)`, loads the module
with all deps mocked and calls `initialize()` with a stub context.

**Catches:** runtime errors invisible to AST — bad attribute access, wrong
argument count, unguarded `None` dereferences, broken lazy imports inside
`initialize()`.

**Context setup:** `get_main_window()` returns a `MagicMock` (not `None`)
because plugins are allowed to assume the host app is fully started when
`initialize()` is called.

---

### `test_run.py` — `run()` / `autorun()` smoke tests (13 tests)

For every visible plugin that uses the legacy entry-point API (`run(main_window)`
or `autorun(main_window)`) without an `initialize()`, calls the entry point
with a `MagicMock` main window.

**Catches:** same class of runtime errors as `test_initialize.py`, for plugins
that pre-date the `PluginContext` API.

---

### `test_save_load.py` — Save/load handler round-trips (6 tests)

For plugins that call `context.register_save_handler()` during `initialize()`,
captures the registered handler and exercises the full persistence round-trip:

1. `save_handler()` — must not raise
2. `load_handler(save_data or {})` — must not raise

`save_handler` may legitimately return `None` (nothing to save); in that case
`load_handler` is called with `{}`, which every well-behaved handler must
tolerate.

**Catches:** crashes in serialisation / deserialisation logic that only surface
when the user saves or loads a project file.

Currently covers: **Atom Colorizer**, **XYZ Editor**, **Settings Saver**.

---

### `test_plugin_<name>.py` — per-plugin test files

Every visible plugin with non-trivial logic (pure functions, save/load
handlers, dialog methods extracted via AST, etc.) gets its own flat file in
`tests/`, named after its plugin: `test_plugin_atom_colorizer.py`,
`test_plugin_paste_xyz.py`, `test_plugin_povray_export.py`, and so on. Each
file typically:

- Loads the plugin module once at module scope inside `mock_optional_imports()`.
- Defines any plugin-specific fakes it needs (widget stand-ins, fake rdkit
  atoms/mols) — shared ones live in `conftest.py` (see below).
- Groups tests into one `Test*` class per method/function under test.

There is no 1:1 mapping to a single "kind" of test — a per-plugin file may
mix `initialize()` registration checks, pure-function unit tests, and
Qt-method extraction tests for the same plugin, since they all belong to the
same plugin's coverage.

`test_plugin_installer.py` is the deepest example: version comparison,
settings round-trip, startup-check scheduling, dependency parsing, chunked
downloads, folder-plugin overwrite logic, `populate_table` / `filter_plugins`
UI logic, and app-version detection — all for the `Plugin_Installer` plugin.

Shared helpers used across many per-plugin files, defined once in
`conftest.py`:

| Symbol | Purpose |
|---|---|
| `mocks_with_real_numpy()` | Like `mock_optional_imports()` but keeps real numpy in `sys.modules` (for export/parser generators doing real vector math) |
| `P3`, `FakeAtom`, `FakeBond`, `FakeConf`, `FakeMol` | Minimal fake rdkit bond-graph objects |
| `extract_function(path, class_name, fn_name, extra_globals=None)` | AST-based method/function extractor for methods on Qt-derived (mocked-base) classes |

Some plugin-specific fake-widget sets (e.g. `FakeSpin`/`FakeCombo`/`FakeCheck`
for QM input-generator dialogs) are *not* centralised — they differ enough
between plugins that each per-plugin file keeps its own copy.

---

### `test_shared_<topic>.py` — cross-plugin contract/regression files

A small number of test classes are genuinely about more than one plugin at
once (parametrized across the three Chat-with-Molecule-Neo variants, or a
regression that must hold for every legacy QM input generator). These live
in `tests/test_shared_<topic>.py` instead of being duplicated per plugin:

| File | Covers |
|---|---|
| `test_shared_chat_variants.py` | Chat Neo ChatGPT + Local: settings round-trip, `latex_to_html` fallback, `PubChemResolver` guards, `run()` smoke |
| `test_shared_ai_tool_parsing.py` | Chat Neo Gemini/ChatGPT/Local: tool-call regex, `collect_tools`, SMILES links, module constants, `append_log` — parametrized across all three variants |
| `test_shared_chat_optimizer.py` | Chat Neo Gemini/ChatGPT/Local: history-pruning regressions, Local `log_usage` |
| `test_shared_input_generator_guards.py` | All five legacy QM input generators (MOPAC/GAMESS/PySCF/Psi4/NWChem): `run(mw)` no-molecule warning guard |

---

### `test_menu_registration.py` — Menu registration contract (≈43 tests)

Parametrized over every visible `initialize()`-plugin. Asserts that at least
one `PluginContext` registration method is called so the plugin is reachable.

Recognised registration methods: `add_menu_action`, `add_export_action`,
`add_analysis_tool`, `add_plugin_menu`, `register_file_opener`,
`register_3d_style`, `register_save_handler`, `register_load_handler`,
`register_document_reset_handler`, `add_toolbar_action`.

**Exempt** (intentional — no menu registration in `initialize()`):

| Plugin | Reason |
|---|---|
| Dark Mode Theme | `autorun()` applies a stylesheet at load; no menu entry |
| Plugin Installer | `run(mw)` is auto-registered by the host; `initialize()` only schedules checks |
| All-Trans Optimizer | `run(mw)` auto-registered; `initialize()` stores the launch closure |
| Complex Molecule Untangler | same |
| PubChem Name Resolver | same |
| Vector Viewer | same; `initialize()` calls `show_status_message()` but no menu |

---

### `test_api.py` — Static API compatibility check (1 test)

Uses `api-checker/check_api.py` to perform a two-phase AST scan:

1. **Phase 1** — Extract the `MainWindow` API surface from the main app source
   (methods, manager classes, `self.host.X = ...` assignments).
2. **Phase 2** — Scan all visible plugins for `mw.attr` and
   `mw.manager.attr` accesses; report any attribute not found in Phase 1.

Both allowlists are active (`--default-allowlist --mw-allowlist`), suppressing
known false positives:

- Manager attrs set via `self.host.manager.X = ...` (AST-invisible)
- V2 legacy bridge attrs (`mw.host`, `mw.view3d`) — always guarded by `hasattr`

Expected result: **0 issues**.

**Skipped automatically** when `../python_molecular_editor` is not present.
On CI the `test-api` job clones the main app at `--depth 1` before running.

**Catches:** new plugins (or changes to existing plugins) that call a
`MainWindow` or manager attribute that no longer exists in the main app.

---

## CI

Two jobs in `.github/workflows/test-plugins.yml`:

| Job | Python matrix | Main app cloned | Tests run |
|---|---|---|---|
| `test` | 3.11, 3.12, 3.13 | No | All except `test_api` (skipped) |
| `test-api` | 3.11 | Yes (`--depth 1`) | `test_api.py` only |
