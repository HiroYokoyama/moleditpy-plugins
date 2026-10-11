# Gaussian FCHK Loader and Vector Viewer audit

Both plugins are version 2026.10.13. Archived copies are excluded.

## Loader

- Discovery skips `_old` and cache/environment directories, and searches the
  correct directory when the loader is installed directly in `plugins`.
- MO packages load from their selected `__init__.py`; relative imports work
  without modifying `sys.path` or reloading live widget classes.
- Failed imports restore the previous module and its submodules.
- Import and input validation precede closing existing analyzers. Missing
  entry points and failed preflight checks report an error and restore the
  buttons according to available plugins.
- All running calculations are checked before closing an existing dock. A
  running MO worker is asked to stop; its dock remains alive until it finishes,
  after which closing or retrying is safe.
- Frequency docks connect visibility cleanup; title-bar closure invokes plugin
  cleanup. MO dock closure removes its named isosurface actors.

## Vector Viewer

- COM uses atomic masses and finite conformer coordinates, including isotope
  masses supplied by RDKit. Missing/invalid conformers are reported rather than
  anchoring a molecular vector at a made-up origin.
- Exactly three finite components are required. Commas and Fortran D exponents
  are supported. Empty, malformed, zero and overflowing vectors clear stale
  arrows. Input errors appear beside the vector controls.
- Molecule replacement clears the old vector and arrow. Geometry changes on the
  same molecule recenter the arrow; reversal updates immediately.
- Actors are removed from the plotter that created them. Closed windows stop
  their timer and clear the actor; reopening redraws a still-valid vector.
- PNG export reports filesystem failures.

## Validation

Tests exercise import rollback, relative package imports, archived discovery,
failed preflight preservation, worker-aware dock closure, MO surface cleanup,
mass-weighted COM, real PyVista arrow construction, invalid inputs, geometry
changes, molecule replacement, reopening and export failures. Unit and GUI
suites are run separately; combined coverage exceeds 80% for both plugins.
