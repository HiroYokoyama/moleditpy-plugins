# Licensing and independently developed plugins

The first-party source code in **moleditpy-plugins** is licensed under the GNU Lesser
General Public License, version 3 or (at your option) any later version
(LGPL-3.0-or-later). `LICENSE` contains LGPLv3's additional permissions;
`COPYING` contains the incorporated GPLv3 terms. Earlier releases retain their
original licenses.

Copyright: Hiromichi Yokoyama / HiroYokoyama and the contributors identified in
the existing copyright notices. Retain those notices, `LICENSE`, `COPYING`, and
this document when redistributing the covered source or binaries.

Modifications to the covered code remain subject to the LGPL. Distribution
terms must allow users to modify the covered components and reverse engineer
for debugging their modifications as required by LGPLv3. Provide corresponding
source and allow modification or replacement of covered libraries when
distributing a combined binary, as required by the license.

Independently developed plugins that use MoleditPy's documented plugin API may
be distributed under any license, including proprietary licenses, under the
[main application's copyright-holder permission](https://github.com/HiroYokoyama/python_molecular_editor/blob/main/LICENSING.md).
Their authors choose their licenses and retain their copyright. This repository's
LGPL license does not change the licenses of other independently developed
plugins, copied MoleditPy code, or third-party components.

Third-party libraries, executables, models, and data retain their own licenses
and obligations. Preserve their notices and comply with the terms of the
versions actually redistributed, including native and transitive dependencies.
Running a program in a separate process does not by itself establish license
compatibility.

## Collection exceptions

The grant covers the maintained first-party plugins except
`plugins/OpenBabel_Conversion_Tool/`, which keeps its existing GPL terms.
`plugins/_old/` is excluded and remains byte-for-byte unchanged under its
existing license. Third-party remote entries in REGISTRY are pointers to their
authors' distributions and do not adopt this repository's license.

Open Babel's current official metadata/library documentation specifies GPLv2.
No claim is made that an in-process Open Babel combination is LGPLv3-compatible.
Resolve that combination before distributing it; changing a wrapper label does
not change the engine's license. Source: https://openbabel.org/docs/UseTheLibrary/intro.html

## Dependency audit

See [the MoleditPy 5 plugin licensing audit](docs/PLUGIN_LICENSE_AUDIT.md) for
repository scopes, inspected dependency versions, packaging checks and exceptions.
