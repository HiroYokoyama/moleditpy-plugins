# Security Policy

## Supported Versions

The latest maintained collection for MoleditPy 5 / PySide6 receives security updates. MoleditPy 4 and earlier snapshots, including `2026.10.19-MoleditPy_4`, are historical and are not security-maintained. Standalone plugins have their own release and security policies; registry listing does not extend their support period.

## Trust and scope

- Plugins execute Python with the same permissions as MoleditPy and are not sandboxed. Install only code and dependencies from publishers you trust.
- Registry SHA-256 verification detects a mismatch against catalogued bytes. It is not a code audit, publisher authentication or a guarantee that a plugin is safe.
- Reports may concern downloaded code, archive/path handling, molecular-file parsers, external calculation commands, credentials or optional network integrations. Include the plugin ID, version, download URL and enabled features, without sending real credentials or private calculation data.
- Network-capable plugins can expose local HTTP/MCP services or use remote SSH/API endpoints. Their access controls and binding defaults are part of their security scope; report them to the owning plugin repository as well as this collection when relevant.

## Reporting a Vulnerability

**Please do not open a public issue for security vulnerabilities.**

To report a vulnerability, please contact the maintainers at titech.yoko.hiro[at]gmail.com.

Include the following in your report:

- A description of the vulnerability and its potential impact
- Steps to reproduce or a proof-of-concept
- The version of the software affected
- Your operating system and environment details

If the issue is confirmed, a fix will be released as soon as possible.
