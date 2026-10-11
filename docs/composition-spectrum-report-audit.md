# Composition, spectrum and report audit

Scope: Elemental Analysis, MS Spectrum Simulation Neo and Compound Info Report. Changed plugin versions: 2026.10.13. Archived plugins excluded.

## Elemental Analysis

- Reject unmatched/empty parentheses, misplaced charges, zero atom counts and oversized formula strings.
- Preserve explicit RDKit isotope labels in synchronized formulas; support [13C]-style labels and pure D/T masses. Unknown isotopes are rejected. Hill ordering recognizes labeled carbon, and publication percentages combine isotopes under the corresponding element.
- Combine fractional solvate counts directly instead of reparsing decimal atom counts as an integer formula. Preserve small fractional counts in tables and the entered solvate notation in reports.
- Respect disabled synchronization, clear results on an empty document and restart polling when the retained dialog reopens.

## MS Spectrum Simulation Neo

- Apply the same strict formula syntax and explicit-isotope support, preserving pure isotope distributions and their exact/average masses.
- Initialize molecular ions with their actual formal charge. Track charge synchronization separately from formula text, restore synchronization after manual edits and clear stale data on document reset.
- Bound simulations to 10,000 atoms per ion and Gaussian grids to 500,000 points; show the atom limit in the interface and retain a stick spectrum when a profile cannot be generated.
- Preserve probability-weighted masses within convolution bins, preventing repeated rounding from shifting large isotope-labeled ions. Regression: D1000 retains its expected exact mass within 1e-8 Da.
- Respect disabled synchronization and restart the polling timer on reopen.

## Compound Info Report

- Bound PubChem requests with a 10-second timeout and an 8 MiB response limit. Cancellation stops subsequent requests and report generation; report reentry is guarded. An in-flight synchronous request completes or times out before cancellation takes effect.
- Validate CAS number syntax and its check digit; retain units for numeric density values.
- Report existing charged molecules at their actual m/z instead of applying neutral-molecule adduct templates; omit deprotonation for hydrogen-free molecules.
- Preserve isotope labels and render isotope mass numbers as superscripts.
- Copy the active molecule before generating 2D coordinates; cache its first structure image so subsequent reports remain consistent with that molecule snapshot.
- Imports in all three plugin main files are at module head.

## Validation

140 targeted unit tests and 171 targeted GUI tests passed. Combined coverage: Elemental Analysis 89%, MS Spectrum Simulation Neo 90%, Compound Info Report 95%.

Regression tests include real RDKit labeled-isotope composition and ion masses, doubly charged ions, fractional solvates, malformed/large inputs, bounded Gaussian grids, sync lifecycle, density units, CAS checksum validation, request bounds, cancellation and preservation of the active molecule's conformers.

References: [RDKit formula API](https://rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html), [CAS check-digit verification](https://www.cas.org/ja/training/documentation/chemical-substances/checkdig).

Full local suites: 2,942 unit tests passed (122 existing skips, 3 existing warnings), 3,309 GUI tests passed (57 existing skips). Registry metadata and hashes were regenerated with the repository script.
