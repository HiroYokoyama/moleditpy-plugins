"""Generate pubchem_structure_identifier_truststore.py from the original."""
from pathlib import Path

SRC = Path("plugins/PubChem_Structure_Identifier/pubchem_structure_identifier.py")
DST = Path("plugins/PubChem_Structure_Identifier/pubchem_structure_identifier_truststore.py")

# Block inserted after "import urllib.parse"
TRUSTSTORE_IMPORT_BLOCK = (
    "import ssl\n"
    "\n"
    "# --- SSL Truststore for Corporate Environments ---\n"
    "try:\n"
    "    import truststore\n"
    "except ImportError:\n"
    "    truststore = None\n"
    "\n"
)

# Static method inserted inside PubChemResolver class, before resolve_name_to_smiles
SSL_CONTEXT_METHOD = (
    "    @staticmethod\n"
    "    def _create_ssl_context():\n"
    '        """\n'
    "        Create an SSL context that uses the system trust store if available.\n"
    "        Fixes SSL errors in corporate environments with security software/proxy.\n"
    '        """\n'
    "        ctx = None\n"
    "        if truststore:\n"
    "            try:\n"
    "                ctx = truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)\n"
    "            except Exception:\n"
    "                ctx = ssl.create_default_context()\n"
    "        else:\n"
    "            ctx = ssl.create_default_context()\n"
    "        try:\n"
    "            ctx.set_ciphers('DEFAULT@SECLEVEL=1')\n"
    "        except (ssl.SSLError, AttributeError):\n"
    "            pass\n"
    "        return ctx\n"
    "\n"
)

text = SRC.read_text(encoding="utf-8")
lines = text.splitlines(keepends=True)
out = []

for line in lines:
    # 1. Patch PLUGIN_NAME
    if line.startswith("PLUGIN_NAME = ") and "(truststore)" not in line:
        stripped = line.rstrip("\n").rstrip("\r")
        line = stripped[: stripped.rindex('"')] + ' (truststore)"\n'

    # 2. Patch PLUGIN_ID
    elif line.startswith("PLUGIN_ID = ") and "_truststore" not in line:
        stripped = line.rstrip("\n").rstrip("\r")
        line = stripped[: stripped.rindex('"')] + '_truststore"\n'

    # 3. Insert truststore import block after "import urllib.parse"
    elif line.strip() == "import urllib.parse":
        out.append(line)
        out.append(TRUSTSTORE_IMPORT_BLOCK)
        continue

    # 4. Insert _create_ssl_context before "    @staticmethod" of resolve_name_to_smiles
    elif line.strip() == "@staticmethod" and out and "BASE_URL" in "".join(out[-5:]):
        out.append(SSL_CONTEXT_METHOD)
        out.append(line)
        continue

    # 5. Add context= to urlopen calls
    elif "urllib.request.urlopen(" in line and "context=" not in line:
        line = line.replace(
            "urllib.request.urlopen(",
            "urllib.request.urlopen(",
        )
        # Insert context argument: urlopen(url) -> urlopen(url, context=...)
        # and urlopen(url, timeout=N) -> urlopen(url, context=..., timeout=N)
        stripped = line.rstrip("\n").rstrip("\r")
        indent = len(stripped) - len(stripped.lstrip())
        if stripped.rstrip().endswith(") as response:"):
            # Simple: urlopen(expr) as response:
            stripped = stripped.rstrip().rstrip(") as response:").rstrip()
            # Find last open paren content
            inner = stripped[stripped.index("urlopen(") + len("urlopen("):]
            prefix = stripped[: stripped.index("urlopen(") + len("urlopen(")]
            line = " " * indent + prefix.lstrip() + inner + ", context=PubChemResolver._create_ssl_context()) as response:\n"
        else:
            # Fallback: just append note
            pass

    out.append(line)

result = "".join(out)
DST.write_text(result, encoding="utf-8")
print(f"Written: {DST}")

# Quick check
import subprocess, sys
result2 = subprocess.run([sys.executable, "-c", f"import ast; ast.parse(open(r'{DST}').read()); print('Syntax OK')"],
                        capture_output=True, text=True)
print(result2.stdout.strip() or result2.stderr.strip())
