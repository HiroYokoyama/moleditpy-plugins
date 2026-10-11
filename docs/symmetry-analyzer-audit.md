# Symmetry Analyzer audit

Scope: active in-repository Symmetry Analyzer only. Version 2026.10.13.

- Closing, Escape and acceptance cancel cooperatively and defer closure until native QThread completion. Workers remain owned by the dialog; no forced termination or GUI-thread wait. Scan reentry is blocked until the native completion callback releases the worker.
- Document reset or changed coordinates invalidate queued results and visualizations. Empty documents cannot receive late results. Starting another scan clears old analyzer state and overlays.
- Nonfinite coordinates and dummy elements are refused before pymatgen construction. Output coordinates must be finite and match the complete atom count before any coordinate write or undo checkpoint.
- The tolerance scan includes the requested upper endpoint even when the range is shorter than the grid spacing.
- Rotation powers retain their finite order: 144 degrees is a fivefold rotation, and 135 degrees is an eightfold rotation. Labels and sorting use the rational angle fraction.
- Overlays are removed from the plotter that created them. Their center follows the selected pymatgen analyzer's center of mass, matching its elemental-mass frame even when RDKit atoms carry isotope masses.
- Plugin imports are at module head, retaining optional dependency fallbacks.

Scientific reference: [pymatgen PointGroupAnalyzer documentation](https://pymatgen.org/pymatgen.symmetry.html). The local pymatgen-core implementation was also inspected for its centered-molecule frame.

Regression tests cover real Qt workers through close/reject/accept, stale and empty documents, invalid input/output, tolerance endpoints, rotation powers, plotter ownership, isotope frame alignment, and a real pymatgen/scipy translated-water projection.

Validation: 102 targeted unit and 84 targeted GUI tests passed; combined plugin coverage 94%. Full suites: 2,920 unit tests passed (122 existing skips, 3 existing warnings), 3,294 GUI tests passed (57 existing skips). Registry metadata and hash regenerated with the repository script.
