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
user choices. Change the 3D background through Settings if needed. To restore
automatic icon color selection, remove `icon_foreground` from that JSON file while
MoleditPy is closed (keep a backup first).
