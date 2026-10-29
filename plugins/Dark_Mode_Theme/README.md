# Dark Mode Theme

Applies a dark stylesheet to the application UI and makes the built-in bond,
ring-template, and chain toolbar icons white while the plugin is loaded.
It does not modify the application's settings dictionary or write any settings
to JSON. The saved 2D and 3D background colors and `icon_foreground` are preserved.
Third-party toolbar icons are not recolored.

Disable the plugin and restart MoleditPy to return to the application's normal
appearance. Disabling it in an already running session does not undo its stylesheet
or icon colors because the host does not provide a plugin unload callback.

## 2026.10.30

- Fixed startup overriding the user's saved 3D background color.
- Replaced the persistent white `icon_foreground` override with in-memory tinting
  of the built-in toolbar icons, deferred until the host has created the actions.

Earlier versions may already have saved `background_color: "#2b2b2b"` and
`icon_foreground: "#FFFFFF"` in `~/.moleditpy/settings.json`. This update leaves
those existing values intact because they cannot be distinguished from intentional
user choices. **Updating or disabling the plugin does not clear these previously
saved overrides.**

After updating, use **Settings > Reset All Settings** from the application's main
menu and restart MoleditPy to clear the old overrides. This resets all application
preferences, including your chosen background color; set your preferred colors
again afterward. Use the main-menu reset action, which deletes the settings file
and reloads the defaults, so plugin-added keys are removed as well.

To keep your other preferences, you can instead close MoleditPy, back up
`~/.moleditpy/settings.json`, restore the desired `background_color`, and remove
`icon_foreground` from that JSON file to restore automatic icon color selection.
