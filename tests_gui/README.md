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

Only **registry-visible** plugins with Qt dialog classes are covered here.

| Plugin | Dialog | Tests |
|---|---|---|
| **Paste XYZ** | `PasteXYZDialog` | Window title, button labels, `get_data()` round-trip, placeholder text |
| **Gaussian Input Generator Neo** | `RouteBuilderDialog` | 4 tabs present, preview label, method/basis combos, job-type group visibility |
| **MS Spectrum Simulation Neo** | `MSSpectrumDialog` | Window title, formula input, sync checkbox disabled with no main window, adduct combo |

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
