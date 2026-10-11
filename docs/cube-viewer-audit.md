# Cube viewer audit — 2026-10-12

Affected plugins: Cube File Viewer, Advanced Cube File Viewer and Mapped Cube
Viewer, all bumped to **2026.10.13** at the user's request.

## Reference implementations

- `moleditpy_pyscf-calculator/pyscf_calculator/vis.py`: reuse its strict parser
  checks and removal of the old mapped actor before computing a new contour.
- `moleditpy_orca_result_analyzer_plugin/orca_result_analyzer/vis.py`: confirm
  dataset IDs follow the atom block, datasets interleave per voxel, and the
  Cube Z-fast order becomes VTK X-fast order by C-reshape/F-flatten.
- [PySCF Cube writer](https://pyscf.org/_modules/pyscf/tools/cubegen.html):
  scalar records can have fewer than six values; line wrapping is independent
  of the grid dimensions.
- [PyVista sample](https://docs.pyvista.org/api/core/_autosummary/pyvista.dataobjectfilters.sample):
  `vtkValidPointMask` identifies whether each point was successfully sampled.

The plugins install as individual files, so their parsers remain self-contained.
The shared regression suite checks each copy against the same cases.

## Confirmed fixes

1. Valid short scalar lines no longer get skipped. Truncated, excessive,
   nonnumeric or nonfinite payloads raise errors instead of getting padded
   with zeros, shifted or silently accepted.
2. Negative NATOMS dataset identifiers are read after the atoms. Positive
   NATOMS NVAL datasets are also deinterleaved correctly. Fortran D exponents
   are accepted in numeric records.
3. Malformed/nonfinite atom coordinates and invalid/singular grids are rejected.
   The PySCF parser's limits are retained: 256 MiB files and two million scalar
   values. This is an explicit loading limit rather than silent truncation.
4. Cube replacement parses and constructs the molecule before closing the old
   viewer. Cancelling bond perception or failing to read a new file preserves
   the old view. Fallback conformers have the correct atom count.
5. Mapped replacement closes the old viewer before installing the new molecule.
   The previous order could clear the newly installed molecule.
6. Empty/failed mapped updates remove the previous actor. Failed samples are
   reported rather than shown as valid zero-valued properties. Invalid color
   ranges are reported. Colorbar plotters close even when export fails.
7. Closing the mapped viewer no longer invokes the host's entire-document reset.
   Widget closure removes its actor; a molecule change detaches the stale view.
   Basic/advanced Close also protect a replacement molecule before the watch
   timer has fired.

## Validation

Tests cover parsing with real NumPy; grid coordinates and interpolation with
real PyVista/VTK; and real Qt controls, closure and molecule replacement. An
analytic linear property confirms colors come from the property grid rather
than the density contour. Distinct grid resolutions are allowed when the
property grid covers every contour point. This checks viewer behavior, not
whether a supplied scientific dataset is physically correct.

## Post-push rescan

Scalar conversion stops after the expected count plus one, bounding allocation
for excessive payloads. Dock title-bar closure now invokes plugin cleanup in all
three viewers; real Qt tests cover that path. Rendering CI explicitly includes
the real-data regression module, whose optional imports are gated on all
required scientific libraries.

## Related plugins

Orbital Comparator and Gaussian MO Analyzer use the same strict PySCF-derived
Cube parsing rules, including NVAL, wrapped DSET_IDS, finite records and bounded
conversion. Orbital Comparator preserves a loaded slot when replacement fails.
Gaussian MO Analyzer's cache includes source content, grid size and margin;
its grid control is capped at 125 points per axis to fit the two million scalar
loading limit. Its downloadable ZIP is rebuilt and checked against each source
file. Both plugin versions are 2026.10.13.
