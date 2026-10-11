# Important plugin audit

Scope: Conformational Search, xTB Optimizer, Step Optimizer, Structural Updater, XYZ Editor, Paste XYZ, Charge Editor, Bond Editor and Plugin Installer. Archived plugins excluded. Changed plugins use 2026.10.13.

- Conformational Search: closing now requests cancellation and completes closure after the native worker finishes, without indefinite GUI-thread waits. Thread ownership lasts until QThread.finished. Completion, preview, accept and reject check the current molecule identity. New runs clear previous results; missing conformers and nonfinite energies are handled.
- xTB Optimizer: result delivery no longer shadows QThread.finished; worker references survive until native thread completion. Reentry is guarded, failed cancellation waits veto closure, and all returned coordinates are checked before any write.
- Step Optimizer: stop, accept and reject act only on the active target, so a replaced document cannot receive stale checkpoints or coordinate restoration.
- Structural Updater: unload stops its timer, reconnects the original conversion button and clears stored patch state. Reinitialization first unloads the previous instance.
- XYZ Editor: nonfinite input is refused. Coordinate edits retain atom charge, isotope, atom map, chirality and custom properties; surviving bonds retain aromaticity, conjugation, direction and stereo. Typed molecule and bond properties are copied.
- Paste XYZ: honors the first declared frame, excludes its comment line, accepts headerless coordinates, rejects negative counts and nonfinite coordinates.
- Bond Editor: direct and table-driven length edits reject nonfinite/nonpositive values and invalid starting distances.
- Plugin Installer: rejects oversized declared/streamed downloads, detects truncated transfers, and validates archive paths, member count, expanded size and symlinks before extraction. Fetch reentry is blocked; a slow fetch keeps its worker alive and defers closure until thread completion. Existing staged replacement and hash checks were reviewed.
- Charge Editor: reviewed molecule signatures, edits, sanitization, undo, picking and close/reset cleanup; no additional confirmed issue in this pass.

Validation: 581 targeted unit tests and 641 targeted GUI tests passed. Combined coverage is 91–100% for every plugin in scope. CI also installs cryptography explicitly so the genuine encrypted round-trip test executes across supported Python versions.

Full local suites: 2,920 unit tests passed (122 existing skips), 3,274 GUI tests passed (57 existing skips). The final installer lifecycle adjustment additionally passed its 195-test combined regression run with Conformational Search.
