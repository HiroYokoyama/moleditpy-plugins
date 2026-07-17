# tests_gui — Headless GUI Tests

Headless GUI tests for plugin dialog widgets using **real PyQt6** with `QT_QPA_PLATFORM=offscreen`.

## How this differs from `tests/`

| | `tests/` | `tests_gui/` |
|---|---|---|
| PyQt6 | Mocked (MagicMock) | **Real** |
| rdkit / numpy / … | Mocked | Mocked |
| What is tested | Plugin logic, settings, parsers | Widget construction, signals, UI state |
| Qt requirement | None | PyQt6 installed |

## Running locally

```bash
# Easiest — runner script sets everything up automatically
python tests_gui/run_gui_tests.py

# Filter to a specific class or test
python tests_gui/run_gui_tests.py -k PasteXYZ
python tests_gui/run_gui_tests.py -k RouteBuilder -v

# Or use pytest directly (you must set QT_QPA_PLATFORM yourself)
QT_QPA_PLATFORM=offscreen pytest tests_gui/ -v
```

> **Windows note:** If both PyQt6 and PySide6 are installed, use `run_gui_tests.py` — it
> automatically adds the correct Qt6/bin DLL directory to avoid the version conflict.

## What is tested

All registry-visible plugins and UI integrations are covered by headless GUI-tier tests.

Tests are laid out **one file per plugin**, mirroring `tests/test_plugin_<snake_name>.py`:
`tests_gui/test_gui_plugin_<snake_name>.py`. Fakes/helpers shared by more than one plugin's
tests (e.g. the rdkit-style `FakeAtom`/`FakeBond`/`FakeMol` used by both editor plugins) live
in [`gui_test_helpers.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/gui_test_helpers.py) (a plain module, not collected as tests).

| Plugin | Test File | Key GUI Classes & Aspects Tested |
|---|---|---|
| Advanced Rendering | `test_gui_plugin_advanced_rendering.py` | `HideOnCloseDialog`, `AdvancedGraphicsWidget` (presets, sliders, PBR toggles) |
| All-Trans Optimizer | `test_gui_plugin_all_trans_optimizer.py` | `_select_torsions`, `run_plugin` guards (no dialog class) |
| Animated XYZ Giffer | `test_gui_plugin_animated_xyz_giffer.py` | `AnimatedXYZPlayer` (playback controls), `parse_multi_frame_xyz` |
| Atom Colorizer | `test_gui_plugin_atom_colorizer.py` | `AtomColorizerWindow` |
| Blender Export | `test_gui_plugin_blender_export.py` | Module constants, `initialize()`/`run()`, `export_to_blender` flow (no-molecule/no-3D guards, cancel, save, generator error) |
| Bond Colorizer | `test_gui_plugin_bond_colorizer.py` | `BondColorizerWindow`, `initialize()` save/load handlers |
| Bond Editor | `test_gui_plugin_bond_editor.py` | `BondEditorWindow`, `_ClickFilter`, bond-type label maps |
| Charge Editor | `test_gui_plugin_charge_editor.py` | `ChargeEditorWindow` |
| Chat with Molecule Neo (Gemini) | `test_gui_plugin_chat_with_molecule_neo.py` | `ChatMoleculeWindow` construction (shared behavior in `test_gui_shared_chat_variants.py`) |
| Chat with Molecule Neo (ChatGPT) | `test_gui_plugin_chat_with_molecule_neo_chatgpt.py` | `ChatMoleculeWindow` construction |
| Chat with Molecule Neo (Local) | `test_gui_plugin_chat_with_molecule_neo_local.py` | `ChatMoleculeWindow` construction (real QWidget main window) |
| Complex Molecule Untangler | `test_gui_plugin_complex_molecule_untangler.py` | `UntanglerDialog`, `run_plugin` |
| Compound Info Report | `test_gui_plugin_compound_info_report.py` | `ReportDialog` |
| Conformational Search | `test_gui_plugin_conformational_search.py` | `ConformerSearchDialog` |
| Cube File Viewer | `test_gui_plugin_cube_file_viewer.py` | `ChargeDialog`, `CubeViewerWidget` (isovalue/opacity sync, complementary color, settings round-trip, structure-change auto-detach, close cleanup) |
| Cube File Viewer Advanced | `test_gui_plugin_cube_file_viewer_advanced.py` | `ChargeDialog`, `FlexibleDoubleSpinBox`, `CubeViewerWidget` (controls, texture path, settings persistence, watch timer, close) |
| Dark Mode Theme | `test_gui_plugin_dark_mode_theme.py` | module constants, `autorun()` on a real `QWidget`, `initialize()` |
| Dummy Atom Mode | `test_gui_plugin_dummy_atom_mode.py` | `initialize()` with a real `QColor` |
| Encrypted Project | `test_gui_plugin_encrypted_project.py` | `PmeencPlugin` (patch/unpatch with real `QAction`s) |
| GAMESS Input Generator | `test_gui_plugin_gamess_input_generator.py` | `GamessSetupDialog` |
| Gaussian FCHK Loader | `test_gui_plugin_gaussian_fchk_loader.py` | `FCHKLoaderDialog` |
| Gaussian Freq Analyzer | `test_gui_plugin_gaussian_freq_analyzer.py` | `GaussianFCHKFreqAnalyzer`, `SpectrumDialog`, `FCHKParser` |
| Gaussian Input Generator Neo (retired) | `test_gui_plugin_gaussian_input_generator_neo.py` | `RouteBuilderDialog` — whole class skipped |
| Gaussian MO Analyzer | `test_gui_plugin_gaussian_mo_analyzer.py` | `initialize()`, `OrbitalWidget` |
| High Resolution Imager | `test_gui_plugin_high_resolution_imager.py` | Module constants, `initialize()` (no dialog class) |
| MOPAC Input Generator | `test_gui_plugin_mopac_input_generator.py` | `MopacSetupDialog` |
| MS Spectrum Simulation Neo | `test_gui_plugin_ms_spectrum_simulation_neo.py` | `MSSpectrumDialog` (formula parsing, adduct deltas & combo, charge/zoom behavior, sync timer), `HistogramWidget` (real profile-peak local-maxima labels, view range) |
| Mapped Cube Viewer | `test_gui_plugin_mapped_cube_viewer.py` | `MappedCubeSetupDialog` |
| Metadata Saver | `test_gui_plugin_metadata_saver.py` | `MetadataSaverDialog` + save/load handlers |
| Molecule Comparator | `test_gui_plugin_molecule_comparator.py` | `MoleculeComparator` (list management, set-as-reference, selection state, scope, RMSD results/clipboard, cleanup, document reset) |
| NWChem Input Generator | `test_gui_plugin_nwchem_input_generator.py` | `NwchemSetupDialog` |
| OpenBabel Conversion Tool | `test_gui_plugin_openbabel_conversion_tool.py` | `MoleculeSelectionDialog` |
| ORCA Freq Analyzer | `test_gui_plugin_orca_freq_analyzer.py` | `OrcaOutFreqAnalyzer`, `SpectrumDialog`, `OrcaParser` |
| Paste from ChemDraw | `test_gui_plugin_paste_from_chemdraw.py` | `reconstruct_from_flat_text()` extra edge cases |
| Paste XYZ | `test_gui_plugin_paste_xyz.py` | `PasteXYZDialog` |
| Plugin Installer | `test_gui_plugin_installer.py` | `PluginDetailsDialog`, `PluginInstallerWindow`, version/dependency parsing |
| POV-Ray Export | `test_gui_plugin_povray_export.py` | Module constants, `initialize()`/`run()`, `export_to_povray` flow (no-molecule/no-3D guards, cancel, save, generator error) |
| Psi4 Input Generator | `test_gui_plugin_psi4_input_generator.py` | `Psi4SetupDialog` |
| PubChem Name Resolver | `test_gui_plugin_pubchem_name_resolver.py` | `MoleculeResolverDialog` |
| PubChem Structure Identifier | `test_gui_plugin_pubchem_structure_identifier.py` | `MoleculeDetailsDialog`, `PubChemResolver` |
| PySCF Input Generator | `test_gui_plugin_pyscf_input_generator.py` | `PyscfSetupDialog` |
| Python Console | `test_gui_plugin_python_console.py` | `PythonConsoleDialog`, entry points, `ConsoleInput`, `PythonHighlighter` |
| Settings Saver | `test_gui_plugin_settings_saver.py` | `SettingsSaverDialog` |
| Step Optimizer | `test_gui_plugin_step_optimizer.py` | `StepOptimizerDialog` |
| Structural Updater | `test_gui_plugin_structural_updater.py` | `SettingsDialog`, `StructuralUpdaterPlugin` |
| Symmetry Analyzer | `test_gui_plugin_symmetry_analyzer.py` | `SymmetryAnalysisPlugin` |
| VDW Radii Overlay | `test_gui_plugin_vdw_radii_overlay.py` | `VDWConfigWindow` |
| Vector Viewer | `test_gui_plugin_vector_viewer.py` | `VectorViewerPlugin` |
| xTB Optimizer | `test_gui_plugin_xtb_optimizer.py` | `XtbOptimizerDialog`, `XtbWorker` |
| XYZ Editor | `test_gui_plugin_xyz_editor.py` | `XYZEditorWindow` (table load/edit, duplicate & adjust-H guards, XYZ export/clipboard, close reopen-fix, `initialize()` handlers) |
| *(cross-plugin)* | `test_gui_shared_chat_variants.py` | Shared `ChatMoleculeWindow` behavior (`append_message`, `save_settings`, `export_history`) parametrized across all three Chat Neo variants — the logic is copy-shared, so it is tested once for all three modules rather than per file |
| *(cross-plugin)* | [`test_integration_real_context.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_integration_real_context.py) | Validates every visible plugin's `initialize()` against the real `PluginContext` API contract (menu paths, callbacks, file loaders, document states) |

## Adding new tests


1. Confirm the plugin is **visible in the registry** (`REGISTRY/plugins.json`).
2. Load the plugin at module level inside `mock_chemistry_imports()`:
   ```python
   with mock_chemistry_imports():
       _my_plugin = load_plugin_for_gui(PLUGINS_DIR / "My_Plugin" / "my_plugin.py")
   ```
3. Write a test class with a `dlg` fixture that creates and destroys the widget:
   ```python
   class TestMyDialog:
       @pytest.fixture
       def dlg(self, qapp):
           d = _my_plugin.MyDialog(parent=None)
           yield d
           d.destroy()

       def test_window_title(self, dlg):
           assert dlg.windowTitle() == "My Dialog"
   ```
4. Use `isHidden()` (not `isVisible()`) to test widget show/hide state — `isVisible()` only
   returns `True` when the window is actually shown on screen.

## CI

The `test-gui` job in `.github/workflows/test-plugins.yml` runs these tests on every push and PR:

```yaml
test-gui:
  runs-on: ubuntu-latest
  steps:
    - pip install pytest PyQt6
    - run: pytest tests_gui/ -v --tb=short
      env:
        QT_QPA_PLATFORM: offscreen
```

No chemistry stack (rdkit, numpy, etc.) is required — all are mocked by `mock_chemistry_imports()`.
