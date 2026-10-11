# Frequency analyzer audit (2026-10-12)

Affected plugins: Gaussian Freq Analyzer and ORCA Freq Analyzer, version
2026.10.13. Archived plugins are excluded.

## Reference logic

The local ORCA Result Analyzer `parser_spectra.py::parse_frequencies` scopes
IR/Raman/mode tables to the latest frequency calculation and locates intensity
columns by their headers. These rules are adapted into the independent ORCA
frequency plugin. Its fixed normal-mode header offset is not reused: integer
column headers and explicit coordinate row indices are validated instead.

The ORCA manual documents the frequency and IR table formats:
https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/vibrations.html

The Gaussian MO Analyzer FCHK reader accepts packed Fortran numbers and D
exponents. The frequency parser now accepts these while additionally checking
declared counts, finite values, atom/coordinate counts, and mode dimensions.

## Fixed findings

- ORCA tables no longer depend on fixed header line counts; modes are matched
  by their explicit IDs. Earlier calculations cannot donate missing IR data.
- Malformed, incomplete or reordered displacement rows are rejected.
- Soft and imaginary vibrational frequencies are retained. The ORCA reported
  first-vibration index is used when present, rather than a 10 cm^-1 cutoff.
- Missing IR intensities remain missing in the mode list and are excluded from
  the spectrum; they cannot cause arithmetic with `None` or fabricated peaks.
- Gaussian validates FCHK records and atom/mode dimensions, uses the declared
  normal-mode count when available, and accounts for linear geometry in the
  fallback. Reusing a parser resets charge and multiplicity.
- Failed parsing preserves the loaded analyzer. Successful replacement stops
  the previous animation. Gaussian's duplicate UI/molecule reload is removed.
- Timers belong to their widgets. Molecule identity checks stop stale animation
  and prevent geometry reset from redrawing a replaced document. Hiding or
  closing the analyzer stops playback and removes displacement vectors.

## Validation

Regression fixtures cover variable ORCA headers, explicit row order, repeated
jobs, missing intensities, soft/imaginary frequencies, D exponents, FCHK record
counts, linear geometry, parser reuse, and replacement molecule ownership.
Existing numeric animation, spectra and GIF tests are also rerun. These checks
validate parser and viewer behavior, not the physical accuracy of calculations.
