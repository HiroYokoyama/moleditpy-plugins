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
