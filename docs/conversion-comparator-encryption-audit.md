# Conversion, comparator and encrypted project audit

Versions: 2026.10.13. Archived plugins excluded.

## OpenBabel Conversion Tool
Native conversion runs in a separate Qt-managed Python process, with GUI event processing, cancellation and a 30-second timeout. Unsupported formats and input files larger than 32 MiB are rejected before parser launch. Imports cap 500 structures, 10,000 atoms per structure and 32 MiB of converted structure data. Unicode input paths are staged and temporary files cleaned up. Export also runs in the isolated process and replaces the destination only after conversion succeeds. Import reentry and a changed active molecule are guarded.

## Molecule Comparator
Closing requests interruption and refuses closure while alignment runs. A second worker and molecule-list edits are blocked during a search; cancelled or stale results cannot replace entries. Bond-free MCS fallback honors the hydrogen option and checks element order. XYZ parsing rejects nonfinite coordinates and honors an empty frame. CSV uses standard quoting. Imports moved to the module head.

## Encrypted Project
Encrypted saves use a flushed, synced temporary file and atomic replacement; write failures preserve the original and cached state. Truncated and oversized input, oversized plaintext and non-object JSON are rejected. Only authentication failures prompt password retry. QAction argument adaptation uses signature binding, so internal TypeError failures are not retried. Cryptography is declared as a dependency; existing encryption format and key derivation remain compatible.

## Validation
103 targeted unit and 203 targeted GUI tests passed. Combined coverage: OpenBabel 95%, Molecule Comparator 94%, Encrypted Project 96%. Regression tests include live disposable-process timeout/cancellation and GUI responsiveness, malformed-input limits, Unicode cleanup, failed atomic replacement, real cryptography round-trip, comparator worker ownership and CSV quoting.
