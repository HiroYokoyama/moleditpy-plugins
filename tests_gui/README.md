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

| Test File | Covered Plugins / Components | Key GUI Classes & Aspects Tested |
|---|---|---|
| [`test_gui_dialogs.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_dialogs.py) | Paste XYZ, Gaussian Input Generator Neo, MS Spectrum Simulation Neo, MOPAC / GAMESS / PSI4 input generators, Atom Colorizer, Plugin Installer, Python Console | `PasteXYZDialog`, `RouteBuilderDialog`, `MSSpectrumDialogNeo`, input/setup dialogs, `AtomColorizerWindow`, `PluginInstallerWindow`, `PythonConsoleDialog` (widgets, syntax highlighting, version checks) |
| [`test_gui_analysis.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_analysis.py) | XYZ Editor, Symmetry Analyzer, Molecule Comparator, Conformational Search | `XYZEditorWindow` (table rows, edits), `SymmetryAnalysisPlugin` (labels, layouts), `MoleculeComparator` (labels), `ConformerSearchDialog` (methods, settings) |
| [`test_gui_chat_variants.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_chat_variants.py) | Chat with Molecule Neo (ChatGPT, Gemini, Local) | `ChatMoleculeWindow` variants (API key field echo mode, placeholder text, model combo defaults, read-only displays, worker state initializations) |
| [`test_gui_exporters.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_exporters.py) | High-Res Imager, Blender Export, POV-Ray Export, Animated XYZ Player / Giffer | Module registrations, export callbacks, `AnimatedXYZPlayer` slider/playback controls, loop check state, multi-frame XYZ files parsing |
| [`test_gui_freq_analyzers.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_freq_analyzers.py) | ORCA & Gaussian Frequency Analyzers, ORCA & Gaussian Spectrum Dialogs | `OrcaOutFreqAnalyzer`, `GaussianFCHKFreqAnalyzer`, `SpectrumDialog` curve/range controls, vector overlays, and output parser integrity |
| [`test_gui_openbabel_pubchem.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_openbabel_pubchem.py) | OpenBabel Conversion Tool, PubChem Structure Identifier | `MoleculeSelectionDialog` (mols list, selection state), `MoleculeDetailsDialog` (html description, IUPAC copy-to-clipboard), and fallback resolvers |
| [`test_gui_qchem_extras.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_qchem_extras.py) | NWChem Input Generator, PySCF Input Generator | `NwchemSetupDialog`, `PyscfSetupDialog` (basis combos, input options, file save flow) |
| [`test_gui_report_settings_fchk.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_report_settings_fchk.py) | Compound Info Report, Settings Saver, Gaussian FCHK Loader | `ReportDialog` (calculations, HTML preview), `SettingsSaverDialog` (presets, export/import options), `FCHKLoaderDialog` (paths validation) |
| [`test_gui_visualizers.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_visualizers.py) | Cube File Viewer, Cube File Viewer Advanced, Mapped Cube Viewer, VDW Radii Overlay, Vector Viewer | `CubeViewerChargeDialog` (contours, limits), `CubeViewerAdvancedDialog` (isovalue, density maps), `VdwRadiiOverlayDialog`, `VectorViewerDialog` |
| [`test_gui_no_dialog_plugins.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_gui_no_dialog_plugins.py) | Dark Mode Theme, Dummy Atom Mode, Paste from ChemDraw | Non-dialog integrations: application-wide stylesheet applications, status bar signals, CPK color updates, Chem.Atom monkey patching, ChemDraw clipboard data parsing |
| [`test_integration_real_context.py`](file:///e:/Research/Calculation/moleditpy/DEV_MAIN/moleditpy-plugins/tests_gui/test_integration_real_context.py) | Plugin-to-Host Integration | Validates plugins against real `PluginContext` API contract, verifying menu paths, callbacks, file loaders, and document states |

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
