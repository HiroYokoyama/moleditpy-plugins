# Plugin license metadata

The `license` field describes the plugin's first-party code. It records existing
terms; adding metadata does not relicense any plugin. Dependency, executable,
model and dataset terms remain separate. Read the project's notices before
redistributing a combined package.

Initial values were checked against the maintained collection's
[licensing scope](../LICENSING.md) and the license documents at each external
entry's registered release tag. Retired/hidden entries were left unchanged.
Existing plugin source files did not receive new license constants.

Most maintained collection plugins use `LGPL-3.0-or-later`. Exceptions and scope
details include:

| Plugin | Recorded license | Additional context |
| --- | --- | --- |
| OpenBabel Conversion Tool | `GPL-3.0-only` | Excluded from the collection's LGPL grant; Open Babel has its own GPL terms. |
| Open Babel Optimizer | `GPL-3.0-or-later` | Declared in the [v1.0.0 README](https://github.com/HiroYokoyama/moleditpy_openbabel_optimizer/blob/v1.0.0/README.md#license). |
| PySCF Calculator | `GPL-3.0-only` | Its [v5.0.0 licensing notes](https://github.com/HiroYokoyama/moleditpy_pyscf-calculator/blob/v5.0.0/LICENSING.md) also grant a geomeTRIC combination permission. The base-license field does not replace that permission or geomeTRIC's terms. |
| NMR Predictor nmrshiftdb2 | `GPL-3.0-only` | The [v3.0.0 README](https://github.com/HiroYokoyama/moleditpy_nmr_predicator_nmrshiftdb2/blob/v3.0.0/README.md#licenses--credits) distinguishes plugin code from predictor JARs and data. |
| Molibrary Browser | `LGPL-3.0-or-later` | The [plugin licensing document](https://github.com/HiroYokoyama/chem_db_web/blob/2.0.0/moleditpy_plugin/LICENSING.md) applies; the web application has a separate root license. |

When a project only grants GPLv3 without an explicit later-version option, the
field uses `GPL-3.0-only`; it does not infer additional permissions from examples
in the standard license text. Unconfirmed values use `NOASSERTION`, displayed as
**Not specified**. Metadata is informational; the actual license documents govern.

For synchronization and reviewed registry-only changes, see
[CONTRIBUTING.md](../CONTRIBUTING.md#license-metadata).
