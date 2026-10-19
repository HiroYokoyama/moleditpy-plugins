# MoleditPy 5 plugin licensing audit — 2026-10-19

The project owner authorized LGPL-3.0-or-later for eligible first-party plugin
source. This audit records the scope of that migration, not a grant for upstream
libraries or a declaration that every possible combined binary is distributable.
Historical tags/releases and third-party terms remain unchanged. Independent
plugin authors remain free to choose their own license under the host API permission.

## Scope

- 32 repository migrations, including the plugin collection and shared library.
- NMRShiftDB2 predictor: GPLv3 retained. Both predictor JARs match upstream r2333
  byte for byte; all ten classes in each JAR match r2097. CDK 2.11 matches its
  official checksum. Commit `21082a1` corrects existing NOTICE filenames, records
  provenance, and includes the separate database license. Complete source/build
  correspondence for the upstream stripped-down predictor remains unverified;
  binary provenance alone is not a full redistribution clearance. JARs are retained.
- PySCF Calculator: GPLv3 retained with an explicit copyright-holder section 7
  permission for the intended geomeTRIC combination (`0295672`). geomeTRIC 1.1.1
  has BSD 3-clause Non-AI terms; those terms are preserved, not relicensed or
  described as ordinary BSD-3-Clause. The permission covers only first-party
  plugin rights. PySCF itself is Apache-2.0; pyberny is MPL-2.0.
- ORCA Rust result analyzer: excluded at the owner's request; PR #1 closed.
- Collection `plugins/OpenBabel_Conversion_Tool/`: excluded; existing GPL terms retained.
  Open Babel's GPL library terms cannot be changed by relabeling its Python wrapper.
- Collection `plugins/_old/`: unchanged. Remote registry entries retain their authors' terms.
- The new local-only Open Babel optimizer is not included in this migration or published.

Each migrated repository includes `LICENSE` (LGPL additional permissions), `COPYING`
(incorporated GPLv3 text), and `LICENSING.md` (scope and provenance). Release archives,
Python distributions and single-file release companions carry those legal notices.

## Pinned first-party copies

Job Manager's copied Job Server HTTP client retains the recorded upstream commit
and source hashes. The owner's additional LGPL permission covers newly distributed
first-party copies; upstream historical grants are retained in provenance.
Job Server's pinned Job Manager gitlink is unchanged; newly distributed first-party
copies include the owner's permission and the original `UPSTREAM_LICENSE`.
No third-party copyright is relicensed by those permissions.

Four consumers pin the updated canonical periodic shared library. Its runtime
source is byte-identical to the previous pin; the change adds licensing/packaging notices.

## Validation completed locally

- 32 canonical license pairs match the main application; Python/TOML/YAML changes parse.
- 42 focused tests pass: shared release 2, Job Server packaging/release 15,
  Job Manager client 6, registry 19 (plus Job Manager subtests).
- Four Python package families built as wheels and sdists; all eight artifacts
  include the three legal notices. The actual Job Server plugin ZIP also includes
  the root/execution notices and original upstream license.
- Registry update changes only Metadata Saver's content hash; historical plugin
  files and the excluded Open Babel tool are unchanged.
- NMRShiftDB2 and PySCF release ZIPs include their corrected notices and license
  information; all three NMR JARs are unchanged. See `license-fix-verification.json`.
- All 32 LGPL migration commit heads passed applicable GitHub checks, including
  the manually triggered template tests pinned to host `v5.0.0a2`.
- These checks complement existing runtime tests; latest-head GitHub results are
  recorded separately in `lgpl-ci-status.json`.

## Migrated repositories

| Repository | Licensing commit |
| --- | --- |
| [moleditpy-mcp_server](https://github.com/HiroYokoyama/moleditpy-mcp_server) | [86448912](https://github.com/HiroYokoyama/moleditpy-mcp_server/commit/86448912fb81512305ed0adb2b509bfd38a2583f) |
| [moleditpy-periodic-shared](https://github.com/HiroYokoyama/moleditpy-periodic-shared) | [e3b03d5a](https://github.com/HiroYokoyama/moleditpy-periodic-shared/commit/e3b03d5a376b13a2084373cd1a8f1126e1ee7c38) |
| [moleditpy-plugins](https://github.com/HiroYokoyama/moleditpy-plugins) | [4343a0be](https://github.com/HiroYokoyama/moleditpy-plugins/commit/4343a0be0deb7bb6b21b3a769db52be65e1406c9) |
| [moleditpy_3d-molecule-on-2d](https://github.com/HiroYokoyama/moleditpy_3d-molecule-on-2d) | [17a7ce51](https://github.com/HiroYokoyama/moleditpy_3d-molecule-on-2d/commit/17a7ce51790daf471c7784e036e6d0775348950d) |
| [moleditpy_3d_functional_group_replacer](https://github.com/HiroYokoyama/moleditpy_3d_functional_group_replacer) | [66dc0295](https://github.com/HiroYokoyama/moleditpy_3d_functional_group_replacer/commit/66dc0295eec65e44eccb9a05213a1a359e2da500) |
| [moleditpy_auto_rotator](https://github.com/HiroYokoyama/moleditpy_auto_rotator) | [a4b06ba7](https://github.com/HiroYokoyama/moleditpy_auto_rotator/commit/a4b06ba70b0a74a256d0351382c0da120d894f95) |
| [moleditpy_blender-export-pro](https://github.com/HiroYokoyama/moleditpy_blender-export-pro) | [d7f0eb34](https://github.com/HiroYokoyama/moleditpy_blender-export-pro/commit/d7f0eb34dab0c4791807a778b81cd4423f6002b4) |
| [moleditpy_camera_position_setter](https://github.com/HiroYokoyama/moleditpy_camera_position_setter) | [ef53d43c](https://github.com/HiroYokoyama/moleditpy_camera_position_setter/commit/ef53d43c4c269830c2521d94540bd2b38b40cf80) |
| [moleditpy_cif_viewer](https://github.com/HiroYokoyama/moleditpy_cif_viewer) | [b9da7dd1](https://github.com/HiroYokoyama/moleditpy_cif_viewer/commit/b9da7dd128d811584b63b81a1794d45221b71db7) |
| [moleditpy_DECIMER_plugin](https://github.com/HiroYokoyama/moleditpy_DECIMER_plugin) | [b21d72bb](https://github.com/HiroYokoyama/moleditpy_DECIMER_plugin/commit/b21d72bbfd13763d03308e76d2bc9f7072f2d3a5) |
| [moleditpy_gaussian_input_generator_pro](https://github.com/HiroYokoyama/moleditpy_gaussian_input_generator_pro) | [d27c36cf](https://github.com/HiroYokoyama/moleditpy_gaussian_input_generator_pro/commit/d27c36cf2dac80be09a53d96b201a838c47248b4) |
| [moleditpy_group_template](https://github.com/HiroYokoyama/moleditpy_group_template) | [77364018](https://github.com/HiroYokoyama/moleditpy_group_template/commit/7736401808e309320baf97e84be5a1c5f0a15718) |
| [moleditpy_job_manager](https://github.com/HiroYokoyama/moleditpy_job_manager) | [d6f7b502](https://github.com/HiroYokoyama/moleditpy_job_manager/commit/d6f7b5028cfb425e69cf80bd3ecad6cef5aff1a9) |
| [moleditpy_nics_placer](https://github.com/HiroYokoyama/moleditpy_nics_placer) | [e05e48ba](https://github.com/HiroYokoyama/moleditpy_nics_placer/commit/e05e48ba6424c4fc417471f5eacc3b61c7092999) |
| [moleditpy_nmr_predictor_cascade](https://github.com/HiroYokoyama/moleditpy_nmr_predictor_cascade) | [86c7385b](https://github.com/HiroYokoyama/moleditpy_nmr_predictor_cascade/commit/86c7385b21642cc7bc156692e68eb39d2e6e0cce) |
| [moleditpy_opsin_name_to_structure](https://github.com/HiroYokoyama/moleditpy_opsin_name_to_structure) | [6b02eb83](https://github.com/HiroYokoyama/moleditpy_opsin_name_to_structure/commit/6b02eb83efe7bd7c09c63657bc9b2ddee0228d87) |
| [moleditpy_orca_input_generator_pro](https://github.com/HiroYokoyama/moleditpy_orca_input_generator_pro) | [56f7e193](https://github.com/HiroYokoyama/moleditpy_orca_input_generator_pro/commit/56f7e193f1f13bb1a698f221eb8a5482643fc3e1) |
| [moleditpy_orca_result_analyzer_plugin](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_plugin) | [d99d5a77](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_plugin/commit/d99d5a77d6ca8b7c09c239dc446ce477f617e854) |
| [moleditpy_pmeff-plugin](https://github.com/HiroYokoyama/moleditpy_pmeff-plugin) | [edb1044c](https://github.com/HiroYokoyama/moleditpy_pmeff-plugin/commit/edb1044cb97ac235bcbb1ea45ac7ec01b0f6776d) |
| [moleditpy_reaction_sketcher_plugin](https://github.com/HiroYokoyama/moleditpy_reaction_sketcher_plugin) | [82a2e88d](https://github.com/HiroYokoyama/moleditpy_reaction_sketcher_plugin/commit/82a2e88d04df68f47516d8b81c18f8a8a8949bb0) |
| [moleditpy_rotation_giffer](https://github.com/HiroYokoyama/moleditpy_rotation_giffer) | [c66a19af](https://github.com/HiroYokoyama/moleditpy_rotation_giffer/commit/c66a19afacaed533a8997c92552e9cebe9e32de1) |
| [moleditpy-plugin-template](https://github.com/HiroYokoyama/moleditpy-plugin-template) | [901f7643](https://github.com/HiroYokoyama/moleditpy-plugin-template/commit/901f7643c7caa2ed5fedf7b70f0e64c755fd1a7b) |
| [moleditpy_cp2k_input_generator](https://github.com/HiroYokoyama/moleditpy_cp2k_input_generator) | [4c5351df](https://github.com/HiroYokoyama/moleditpy_cp2k_input_generator/commit/4c5351dfb14605497aa7ab1df0caa3d2b7a3d6ac) |
| [moleditpy_job_server](https://github.com/HiroYokoyama/moleditpy_job_server) | [fd881e20](https://github.com/HiroYokoyama/moleditpy_job_server/commit/fd881e20abaf39623b0b2912b7179106e095569d) |
| [moleditpy_molecule_styler](https://github.com/HiroYokoyama/moleditpy_molecule_styler) | [0add3933](https://github.com/HiroYokoyama/moleditpy_molecule_styler/commit/0add39332fcf6bcdeef66d355242006a37f43195) |
| [moleditpy_orca_nics_analyzer](https://github.com/HiroYokoyama/moleditpy_orca_nics_analyzer) | [4ae19821](https://github.com/HiroYokoyama/moleditpy_orca_nics_analyzer/commit/4ae1982184b061f4d902cd6b6940775718f223ee) |
| [moleditpy_pets_in_editor](https://github.com/HiroYokoyama/moleditpy_pets_in_editor) | [c787eb67](https://github.com/HiroYokoyama/moleditpy_pets_in_editor/commit/c787eb67182ae5ab62408c275160a3db8d3e4bef) |
| [moleditpy_quantum_espresso_input_generator](https://github.com/HiroYokoyama/moleditpy_quantum_espresso_input_generator) | [1c2e18df](https://github.com/HiroYokoyama/moleditpy_quantum_espresso_input_generator/commit/1c2e18dfc12d1bbd0b3b946aaf0e50930f9ffe3f) |
| [moleditpy_reshiki_importer](https://github.com/HiroYokoyama/moleditpy_reshiki_importer) | [17ddaaa7](https://github.com/HiroYokoyama/moleditpy_reshiki_importer/commit/17ddaaa7fbe813cb6bd83749435a171bf2ae6fab) |
| [moleditpy_slab_builder](https://github.com/HiroYokoyama/moleditpy_slab_builder) | [ab54e078](https://github.com/HiroYokoyama/moleditpy_slab_builder/commit/ab54e07862da31135ae3de4c71cabd9987220bd7) |
| [moleditpy_strain_homodesmotic_reaction_generator](https://github.com/HiroYokoyama/moleditpy_strain_homodesmotic_reaction_generator) | [a7e74d0c](https://github.com/HiroYokoyama/moleditpy_strain_homodesmotic_reaction_generator/commit/a7e74d0c7f5152076589134fd23694ef5e52ebc4) |
| [moleditpy_vasp_input_generator](https://github.com/HiroYokoyama/moleditpy_vasp_input_generator) | [6a38a0b4](https://github.com/HiroYokoyama/moleditpy_vasp_input_generator/commit/6a38a0b4bf1723fe2f4eed43cc1636d5a3964113) |

## Published dependency metadata inspected

Versions below are the inspected metadata snapshot, not universal build pins.
Actual shipped versions and native/transitive components need their own notices.

| Dependency | Inspected version | Reported license |
| --- | --- | --- |
| [PySide6-Essentials](https://pypi.org/pypi/PySide6-Essentials/json) | 6.12.0 | LGPL-3.0-only OR GPL-2.0-only OR GPL-3.0-only |
| [rdkit](https://pypi.org/pypi/rdkit/json) | 2026.9.1 | BSD-3-Clause |
| [pyvista](https://pypi.org/pypi/pyvista/json) | 0.49.1 | MIT |
| [pyvistaqt](https://pypi.org/pypi/pyvistaqt/json) | 0.13.1 | MIT |
| [vtk](https://pypi.org/pypi/vtk/json) | 9.7.1 | BSD |
| [numpy](https://pypi.org/pypi/numpy/json) | 2.5.4 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 |
| [scipy](https://pypi.org/pypi/scipy/json) | 1.18.1 | Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers. All rights reserved. Redistribution and use in source and binary forms, with or wit... (see snapshot) |
| [matplotlib](https://pypi.org/pypi/matplotlib/json) | 3.11.2 | License agreement for matplotlib versions 1.3.0 and later ========================================================= 1. This LICENSE AGREEMENT is be... (see snapshot) |
| [Pillow](https://pypi.org/pypi/Pillow/json) | 13.0.0 | MIT-CMU |
| [pymatgen](https://pypi.org/pypi/pymatgen/json) | 2026.9.24 | MIT |
| [ase](https://pypi.org/pypi/ase/json) | 3.29.0 | LGPL-2.1-or-later |
| [tblite](https://pypi.org/pypi/tblite/json) | 0.7.0 | LGPL-3.0-or-later |
| [cryptography](https://pypi.org/pypi/cryptography/json) | 50.0.2 | Apache-2.0 OR BSD-3-Clause |
| [paramiko](https://pypi.org/pypi/paramiko/json) | 5.1.0 | LGPL-2.1 |
| [DECIMER](https://pypi.org/pypi/DECIMER/json) | 2.8.0 | MIT |
| [py2opsin](https://pypi.org/pypi/py2opsin/json) | 1.2.0 | MIT |
| [pyscf](https://pypi.org/pypi/pyscf/json) | 2.14.0 | Apache-2.0 |
| [geometric](https://pypi.org/pypi/geometric/json) | 1.1.1 | See published source license |
| [pyscf-dispersion](https://pypi.org/pypi/pyscf-dispersion/json) | 1.5.0 | Apache-2.0 |
| [pyberny](https://pypi.org/pypi/pyberny/json) | 0.7.0 | MPL-2.0 |
| [nmrsim](https://pypi.org/pypi/nmrsim/json) | 0.7.1 | MIT |
| [fastapi](https://pypi.org/pypi/fastapi/json) | 0.146.0 | MIT |
| [uvicorn](https://pypi.org/pypi/uvicorn/json) | 0.54.0 | BSD-3-Clause |
| [starlette](https://pypi.org/pypi/starlette/json) | 1.9.0 | BSD-3-Clause |
| [anyio](https://pypi.org/pypi/anyio/json) | 4.15.1 | MIT |
| [packaging](https://pypi.org/pypi/packaging/json) | 26.3 | Apache-2.0 OR BSD-2-Clause |
| [psutil](https://pypi.org/pypi/psutil/json) | 7.2.2 | BSD-3-Clause |
| [openai](https://pypi.org/pypi/openai/json) | 3.35.1 | Apache-2.0 |
| [google-genai](https://pypi.org/pypi/google-genai/json) | 2.31.0 | Apache-2.0 |
| [requests](https://pypi.org/pypi/requests/json) | 2.34.2 | Apache-2.0 |
| [markdown](https://pypi.org/pypi/markdown/json) | 3.11.1 | BSD-3-Clause |
| [truststore](https://pypi.org/pypi/truststore/json) | 0.10.5 | MIT |
| [openbabel](https://pypi.org/pypi/openbabel/json) | 3.2.1 | GPL-2.0-only |

## Primary sources and reproducible evidence

- [GNU license compatibility](https://www.gnu.org/licenses/license-compatibility.en.html)
- [GNU GPL FAQ](https://www.gnu.org/licenses/gpl-faq.en.html)
- [GPLv3 with LGPL additional permissions](https://www.gnu.org/licenses/lgpl+gpl-3.0-standalone.html)
- [Open Babel library introduction/license](https://openbabel.org/docs/UseTheLibrary/intro.html)
- [NMRShiftDB2 predictor JAR documentation](https://sourceforge.net/p/nmrshiftdb2/wiki/PredictorJars/)
- [Published geomeTRIC 1.1.1 source archive](https://files.pythonhosted.org/packages/4d/57/f53eb30c38aa4eeabddc3bddf283f4c10616de65d95fcd1f43900721f425/geometric-1.1.1.tar.gz)

Local evidence: `plugin-license-inventory.json`, `plugin-dependency-license-details.json`,
`dependency-official-pypi-licenses.json`, `geometric-published-license.json`,
`lgpl-static-verification.json`, `lgpl-shared-pins.json`, `lgpl-committed.json`,
`lgpl-pushed.json` and the build logs/artifacts in `lgpl-packages/`.
