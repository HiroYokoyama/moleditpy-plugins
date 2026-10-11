# Input generator validation — 2026-10-12

Scope: the five active input generators in this repository. Archived generators
are excluded. Each fix is committed separately with version `2026.10.12`.

## GAMESS (US)

Sources: [official input manual](https://www.msg.chem.iastate.edu/GAMESS/GAMESS_Manual/input.pdf),
[$CONTRL and input format](https://www.msg.chem.iastate.edu/GAMESS/GAMESS_Manual/intro.pdf).

The existing functional names, including `PBE0` and `M06-2X`, are valid.
`NGAUSS` applies to STO and Pople basis sets, not TZV/MINI/MIDI. STO allows
2–6 primitives, N21 allows 3 or 6, and N31 allows 4–6. Polarization counts
`NDFUNC` and `NPFUNC` cannot exceed 3. The dialog now limits these fields,
rejects invalid N21 counts, selects ROHF when an open shell would otherwise
use RHF, and excludes GVB/MCSCF jobs that need configuration blocks this
generator cannot supply. Control records are wrapped to fit 80 columns.

Regression tests cover these constraints, keyword preservation on wrapped
records, and preview/save behavior. Tests validate generated text and real Qt
behavior; a GAMESS executable is not installed in the available environment.

## MOPAC

Sources: [official usage guide](https://openmopac.net/guides/basics/),
[spin keyword](https://openmopac.net/Manual/ms.html),
[force calculation](https://openmopac.net/Manual/force.html),
[thermochemistry](https://openmopac.net/Manual/thermo00.html).

Open-shell defaults now include `UHF`; explicitly requested RHF/CI is preserved.
Incompatible UHF/CI combinations are rejected. Charge and spin controls replace
conflicting copies in free-form keywords, and NOSYM is not duplicated. `1SCF`
is matched as a whole keyword. The thermochemistry template includes `FORCE`;
the supplied geometry must already be optimized, as required by MOPAC.
Logic and Qt regressions validate keyword composition and preview/save behavior.

## NWChem

Sources: [DFT](https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html),
[SCF](https://nwchemgit.github.io/Hartree-Fock-Theory-for-Molecules.html),
[CCSD](https://nwchemgit.github.io/CCSD.html),
[TCE](https://nwchemgit.github.io/TCE.html),
[basis conventions](https://nwchemgit.github.io/Basis.html),
[tasks](https://nwchemgit.github.io/TASK.html).

PBE96, BLYP, LDA and TPSS expand to their documented exchange/correlation
components. Open-shell DFT uses ODFT; correlation jobs include an SCF spin
reference, with NOPEN for multiplicities above eight. Open-shell CCSD uses TCE,
and TCE explicitly selects CCSD. Correlation-consistent and def2 bases use
spherical functions. CCSD derivatives and meta-GGA frequencies request numerical
differentiation. Invalid start names block saving; titles cannot inject lines.
43 logic tests and 39 Qt tests passed. No NWChem executable is installed.

## Psi4

Sources: [SCF references](https://psicode.org/psi4manual/master/scf.html),
[DFT methods](https://psicode.org/psi4manual/master/dft.html).

Existing functional names are preserved. Closed-shell references cannot be used
for open-shell jobs, and HF/correlation jobs reject KS references. Output names
are Python string literals; saving replaces the complete output-setting line,
including names with quotes or parentheses. Imports are at module scope and
files use LF. 37 logic tests and 25 Qt tests passed. Psi4 is not installed.
