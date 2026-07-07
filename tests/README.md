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

### `test_plugin_installer.py` — Plugin Installer unit tests (14 tests)

Deep unit tests for the `Plugin_Installer` plugin specifically.

| Class / test | What it covers |
|---|---|
| `TestCompareVersions` (6) | Version comparison logic: greater, less, equal, mixed-length, date-style versions |
| `TestSettings` (3) | `load_settings` / `save_settings` round-trip, missing file, corrupt JSON |
| `TestInitialize` (5) | Startup check scheduling: `check_at_startup=True/False`, first-run flow, deduplication |
| `test_plugin_version_constant_present` | `PLUGIN_VERSION` constant exists and is non-empty |
| `test_plugin_version_matches_registry` | Source version matches registry entry |

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

### `test_plugin_cube_parsers.py` — Cube file parser tests (29 tests)

Pure-function tests for the three Cube viewer plugins.

| Area | What it covers |
|---|---|
| `parse_cube_data` | Valid files, short-file `ValueError`, excess-line trim, zero-pad, Angstrom header flag, multi-atom |
| `build_grid_from_meta` | Grid shape from metadata dict |
| `read_cube` | End-to-end read returning atom list + flat data array |
| `initialize()` | File-opener registration for `.cube` for all three plugins |

Plugins covered: **Cube File Viewer**, **Cube File Viewer Advanced**, **Mapped Cube Viewer**.

---

### `test_plugin_parsers_exports.py` — Parser and export-script tests (49 tests)

| Area | Plugin | What it covers |
|---|---|---|
| `is_valid_orca_file` | ORCA Freq Analyzer | Keyword variants, 500-line boundary, empty/missing files, case sensitivity |
| `OrcaParser.parse` | ORCA Freq Analyzer | Atom symbols/coords, frequency extraction, multi-block geometry, edge cases |
| `FCHKParser.parse` | Gaussian Freq Analyzer | Atomic numbers, Bohr→Å conversion, Vib-E2 frequency/intensity |
| `find_file_recursive` | Gaussian FCHK Loader | Exact match, wildcard, nested dirs, not found |
| `find_mo_analyzer_module` | Gaussian FCHK Loader | Found, missing `__init__.py`, not found |
| `initialize()` | Blender Export | Export action registered with correct label |
| `initialize()` | POV-Ray Export | Same |
| `initialize()` | ORCA / Gaussian Freq | File-opener registration |

---

### `test_plugin_optimize_resolvers.py` — Optimizer / resolver / viewer tests (≈30 tests)

| Class | Plugin | What it covers |
|---|---|---|
| `TestAllTransOptimizer` | All-Trans Optimizer | `_launch_fn` set by `initialize()`; `run_plugin(mol=None)` shows warning |
| `TestComplexMoleculeUntangler` | Complex Molecule Untangler | `PLUGIN_CONTEXT` and `_launch_fn` set; window reuse; `register_window` called |
| `TestConformationalSearch` | Conformational Search | Menu action registered; `run_plugin(mol=None)` warns; dialog accept; window registered |
| `TestPubChemNameResolver` | PubChem Name Resolver | `initialize()` doesn't raise; `run_search("")` short-circuits; `run(mw)` doesn't raise |
| `TestVectorViewer` | Vector Viewer | `_launch_fn` set; `show_status_message` called |
| `TestMoleculeComparator` | Molecule Comparator | `register_document_reset_handler` called; luminance text-colour formula |

---

### `test_plugin_ui_misc.py` — UI helpers and settings round-trips (48 tests)

| Area | Plugin | What it covers |
|---|---|---|
| `ConsoleInput.append_history` | Python Console | Deduplication, empty-string skip, index tracking |
| `GUIHelp.__repr__` / `__call__` | Python Console | Help text content |
| `MSSpectrumDialog.parse_formula_str` | MS Spectrum Simulation Neo | H₂O, glucose, parentheses, invalid chars |
| `to_subscript` / `to_superscript` | MS Spectrum Simulation Neo | Unicode conversion |
| `get_adduct_delta` | MS Spectrum Simulation Neo | Positive/negative mode, multi-charge adducts |
| `load_settings` / `save_settings` | VDW Radii Overlay | Round-trip, missing file, corrupt JSON |
| `load_library` / `save_library` | Settings Saver | Round-trip, empty library |
| `load_settings` / `save_settings` | Chat with Molecule Neo | Round-trip |
| `latex_to_html` | Chat with Molecule Neo | LaTeX → HTML conversion |
| `reconstruct_from_flat_text` | Paste from ChemDraw | Edge cases: empty, single atom, malformed |

---

### `test_plugin_chat_variants.py` — Chat Neo ChatGPT + Local tests (27 tests)

| Class | What it covers |
|---|---|
| `TestChatGPTSettings` (6) | `load_settings`/`save_settings` round-trip, missing file, corrupt JSON, DEMO_MODE guard, valid JSON output |
| `TestChatGPTLatexToHtml` (5) | No-matplotlib fallback returns `<i>text</i>`; empty string; mocked-matplotlib path; cache hit |
| `TestChatGPTPubChemResolver` (2) | Empty string and `None` → immediate error, no network call |
| `test_chatgpt_run_does_not_raise` (1) | `run(mw)` smoke |
| `TestLocalSettings` (6) | Same battery for the Local variant |
| `TestLocalLatexToHtml` (4) | Fallback + mocked + cache |
| `TestLocalPubChemResolver` (2) | Empty/None guard |
| `test_local_run_does_not_raise` (1) | `run(mw)` smoke |

---

### `test_plugin_pubchem_and_paste.py` — PubChem utilities & Paste XYZ tests (33 tests)

| Class | Plugin | What it covers |
|---|---|---|
| `TestPubChemResolver` (9) | PubChem Structure Identifier | Early-exit on empty/None name+InChIKey; mocked HTTP 200 success; HTTP 404 "not found"; empty properties; network error; `run(mw)` smoke |
| `TestPubChemFetcher` (14) | Compound Info Report | Early-exit on empty/None InChIKey for `get_synonyms`/`get_cid`; mocked success paths; `extract_cas` valid/filtered/capped; `initialize()` registers action |
| `TestPasteXYZParser` (9) | Paste XYZ | Valid 3-atom block; header lines skipped; short lines skipped; non-float coords skipped; non-alpha symbol skipped; empty input; blank lines; fractional coords; mixed-case elements |
| `TestPasteXYZSmoke` (1) | Paste XYZ | `run(mw)` exits cleanly when dialog returns non-Accepted result |

---

### `test_plugin_ai_tool_parsing.py` — AI/LLM tool-call parsing tests (84 tests)

| Class | Tests | What it verifies |
|---|---|---|
| `TestToolRegex` | 7 | Regex extracts single/array/multi-block tool JSON; non-tool JSON found but lacks "tool" key; invalid JSON matched but unparseable; no blocks → empty |
| `TestCollectTools` | 8 | Collection logic: single tool, array of 3, two separate blocks merged, invalid block silently dropped, non-tool JSON and missing "params" excluded |
| `TestSmilesLinkPattern` | 6 | `[Name](smiles:...)` → `<a href=...>`; multiple links; `https://` links not converted; no-link text unchanged |
| `TestLatexToHtmlFallback` | 15 | Fallback (`HAS_MATPLOTLIB=False`) returns `<i>text</i>`; mocked-matplotlib path returns non-empty string — all 3 Chat Neo variants |
| `TestPubChemResolverEarlyReturn` | 12 | `resolve_inchikey_to_name("")/None` and `resolve_name_to_smiles("")/None` return `(None, msg)` without network — all 3 variants |
| `TestModuleConstants` | 24 | `MAX_HISTORY` is even int ≥ 2; `DEMO_MODE=False`; `SYSTEM_PROMPT` non-empty and contains "tool"; `GENERATION_CONFIG` has `temperature` — all 3 variants |
| `TestAppendLog` | 12 | Creates file, writes sender+text, appends across calls, bad path doesn't raise — all 3 variants |

---

### `test_plugin_settings_saver_extended.py` — Settings Saver extended tests (26 tests)

| Class | What it covers |
|---|---|
| `TestGetPluginDataPath` (3) | Returns a string ending in `"settings_saver.json"`; parent dir exists on disk |
| `TestGetLiveSettings` (3) | Returns `mw.init_manager.settings` when present; `None` when `init_manager` absent or has no `settings` |
| `TestSyncLegacySettingsAlias` (5) | Clears+updates distinct dict; same-object no-op; no `init_manager` no-op; non-dict `settings` no-op; exception silently swallowed |
| `TestSaveLoadProject` (6) | `on_save_project` returns `None` when disabled; returns dict when live settings present; `on_load_project` handles empty dict, `"settings"` key, preset storage, non-dict input |
| `TestProjectMode` (5) | `enable_project_mode` sets `EMBED_SETTINGS["enabled"]=True`; `disable_project_mode` sets it `False`; both don't raise; round-trip |
| `TestOnDocumentReset` (4) | No context; mock context; clears `PROJECT_PRESETS`; disables embed when `always_save_to_project=False` |

---

### `test_plugin_encrypted_and_structural.py` — Encrypted Project & Structural Updater tests (33 tests)

| Class | Plugin | What it covers |
|---|---|---|
| `TestEncryptedProjectOnDrop` (7) | Encrypted Project | `on_drop()` returns True for `.pmeenc`/`.PMEENC`, False for other extensions and empty string |
| `TestEncryptedProjectDocumentReset` (2) | Encrypted Project | `on_document_reset()` clears `current_key`/`current_salt`; no-op from None state |
| `TestEncryptedProjectPatchMethod` (6) | Encrypted Project | `_patch_method` replaces `mw` method and marks it with `._is_pmeenc_patch=True`; `_unpatch_method` restores original; double-patch guard; missing-method no-raise |
| `TestEncryptedProjectInitialize` (4) | Encrypted Project | `initialize()` registers file opener for `.pmeenc`, export action, document-reset handler, drop handler |
| `TestStructuralUpdaterInitialize` (4) | Structural Updater | `initialize()` sets `_PLUGIN_INSTANCE`, calls `add_menu_action` with "Structural Updater" in path; `finalize()` does not raise |
| `TestStructuralUpdaterSettings` (8) | Structural Updater | `save_settings`/`load_settings` round-trip; default `enabled=True`; missing file creates default; JSON content verified |
| `TestStructuralUpdaterCheckState` (2) | Structural Updater | `check_state()` exits early when disabled; runs without raise when enabled |

---

### `test_plugin_advanced_and_misc.py` — Advanced rendering, misc visible plugins (42 tests)

| Class | Plugin | What it covers |
|---|---|---|
| `TestParseMultiFrameXYZ` (9) | Animated XYZ Giffer | `parse_multi_frame_xyz`: single/multi-frame files, empty file, blank lines, incomplete frame at EOF, bad coord skipped, non-integer atom count skipped, `run(mw)` no-raise smoke |
| `TestAdvancedRendering` (5) | Advanced Rendering | `get_icon()→None`; `initialize()` registers exactly 4 × `register_3d_style` + 1 × `add_menu_action`; menu path and style key names verified |
| `TestDummyAtomMode` (4) | Dummy Atom Mode | `initialize()` calls `add_toolbar_action` with text `"Dummy Atom *"` and tooltip containing "dummy"; does NOT call `add_menu_action` |
| `TestOpenBabelConversionTool` (6) | OpenBabel Conversion Tool | `initialize()` with mocked openbabel (available): registers export + drop handler; with `OBABEL_AVAILABLE=False`: returns early, no registrations; `open_file_with_openbabel` extensionless path warns and returns |
| `TestHighResolutionImager` (4) | High Resolution Imager | `initialize()` registers export action with "Screenshot" in label; stores `PLUGIN_CONTEXT`; does not call `add_menu_action` |
| `TestXYZEditor` (8) | XYZ Editor | `initialize()` registers menu action, save/load/document-reset handlers; stores `PLUGIN_CONTEXT`; save handler returns `{}` when `current_molecule=None`; load handler tolerates `{}` |
| `TestSymmetryAnalyzer` (6) | Symmetry Analyzer | `initialize()` registers menu action with "Symmetr" in path; stores `PLUGIN_CONTEXT`; `SymmetryAnalysisWorker.run()` always emits `finished(dict, bool)`, `found_any=False` when pymatgen loop runs 0 iterations |

**Technique note:** `SymmetryAnalysisWorker` and `AnimatedXYZPlayer` inherit from Qt classes (mocked → MagicMock metaclass → class body becomes MagicMock). Pure methods are extracted via `ast.get_source_segment` + `exec` into an isolated namespace to avoid the `object.__new__` limitation.

---

### `test_plugin_input_generators_extended.py` — QM input generators (84 tests)

Deep content tests for the MOPAC, GAMESS, PySCF, Psi4 and NWChem input
generators: generated input text (keywords, charge/multiplicity, coordinate
blocks), preset save/load round-trips, save-file rewrite logic, the legacy
`run(mw)` no-molecule warning guard, and NWChem `_scf_spin_lines`
(rhf/uhf + multiplicity keyword selection, parametrized over mult 1–9).

---

### `test_plugin_visualization_extended.py` — Visualization plugins (~50 tests)

Dark Mode stylesheet content and `autorun()`, Atom Colorizer save/load/reset
handlers and color application, Vector Viewer center-of-mass and visualization
update, High Resolution Imager `run()`, VDW Radii Overlay settings edge cases
and `initialize()` detail.

---

### `test_plugin_analysis_extended.py` — Analysis plugins (~70 tests)

ORCA Freq Analyzer final-mode parsing and the `self.context` →
`PLUGIN_CONTEXT` regression; Gaussian Freq Analyzer FCHK parser edge cases
including the inline-scalar (`Charge  I  -2`) regression; Gaussian MO Analyzer
FCHK reader, normalization prefactor and cube writer; Symmetry Analyzer
sort-key/op-label/symbol helpers; Molecule Comparator reset handler;
MS Spectrum Neo gaussian broadening; Compound Info Report property fetch.

---

### `test_plugin_io_export_extended.py` — Export script generation (~40 tests)

Blender and POV-Ray script/scene content generation (materials, cylinders,
per-bond splitting), `_calculate_double_bond_offset` math, and the POV-Ray
triple-bond `off_dir` UnboundLocalError regression.

**Technique note:** uses `mocks_with_real_numpy()`, which restores real numpy
(including *all* `numpy.*` submodules — restoring only the top-level package
lets the mocking MetaPathFinder corrupt numpy process-wide) while keeping
rdkit/PyQt6 mocked.

---

### `test_plugin_chat_optimizer_extended.py` — Chat history management (~15 tests)

History-pruning regressions for all three Chat Neo variants (system prompt
must survive pruning; Local variant must prune `chat_history_state`, not a
nonexistent Gemini `chat_session`) and Local `log_usage`.

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
