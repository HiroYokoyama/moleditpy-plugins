"""Generate truststore variants from the original plugin files."""
from pathlib import Path

TRUSTSTORE_BLOCK = (
    "# --- SSL Truststore for Corporate Environments ---\n"
    "# This must be executed as early as possible to patch SSLContext\n"
    "try:\n"
    "    import truststore\n"
    "    import ssl\n"
    "    # Inject truststore into the system SSL defaults\n"
    "    # This forces Python (urllib, requests, pip, openai) to use the Windows System Certificate Store\n"
    "    # Solution for: \"SSL: CERTIFICATE_VERIFY_FAILED\" in corporate networks\n"
    "    truststore.inject_into_ssl()\n"
    "except ImportError:\n"
    "    # truststore is optional; proceed if not installed\n"
    "    pass\n"
    "except Exception as e:\n"
    "    print(f\"Warning: Truststore injection failed: {e}\")\n"
    "\n"
)

LOCAL_DOCSTRING = (
    '"""\n'
    "=== SSL Certificate Handling (truststore) ===\n"
    "\n"
    "This version uses the `truststore` library to resolve SSL certificate verification errors\n"
    "that occur in corporate environments with security software (e.g., antivirus, proxy).\n"
    "\n"
    "\u3010Problem\u3011\n"
    "- Standard Python uses its own certificate list (certifi) and ignores Windows' trusted certificates\n"
    "- Corporate security software intercepts HTTPS connections with its own certificate\n"
    "- Python rejects this certificate as \"untrusted\" even though Windows trusts it\n"
    "- Result: SSL errors when accessing PubChem API, even though browsers work fine\n"
    "\n"
    "\u3010Solution\u3011\n"
    "- `truststore` makes Python use Windows' system certificate store instead of certifi\n"
    "- Python now recognizes certificates trusted by Windows (including corporate security software)\n"
    "- This is the SAFE solution (no security bypass, just proper integration with OS trust)\n"
    "\n"
    "\u3010Technical Details\u3011\n"
    "- Import `truststore` and call `truststore.inject_into_ssl()` before making HTTPS requests\n"
    "- This patches Python's SSL context to query Windows Certificate Store\n"
    "- Changes persist for the entire process lifetime\n"
    "- Requires: `pip install truststore`\n"
    "\n"
    "\u3010Why This Works\u3011\n"
    "- Windows already trusts the security software's certificate (admin-configured)\n"
    "- truststore bridges Python and Windows, enabling certificate sharing\n"
    "- Future certificate updates on Windows automatically propagate to Python\n"
    '"""\n'
    "\n"
)

INSERT_ANCHOR = "from rdkit.Chem import AllChem\n"


def make_truststore(src: Path, dst: Path, name_suffix: str, id_suffix: str, add_docstring: bool):
    text = src.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)
    out = []

    for line in lines:
        # Patch PLUGIN_NAME: PLUGIN_NAME = "Foo"  ->  PLUGIN_NAME = "Foo (truststore)"
        if line.startswith("PLUGIN_NAME = ") and name_suffix not in line:
            stripped = line.rstrip("\n").rstrip("\r")
            # Remove trailing quote, append suffix and re-close
            line = stripped[: stripped.rindex('"')] + name_suffix + '"\n'

        # Patch PLUGIN_ID
        elif line.startswith("PLUGIN_ID = ") and id_suffix not in line:
            stripped = line.rstrip("\n").rstrip("\r")
            line = stripped[: stripped.rindex('"')] + id_suffix + '"\n'

        out.append(line)

        # Insert docstring right after PLUGIN_ID line (local only)
        if add_docstring and line.startswith("PLUGIN_ID = "):
            out.append("\n" + LOCAL_DOCSTRING)

    text = "".join(out)

    # Insert truststore injection block after rdkit import
    if INSERT_ANCHOR in text:
        text = text.replace(INSERT_ANCHOR, INSERT_ANCHOR + "\n" + TRUSTSTORE_BLOCK, 1)
    else:
        print(f"  WARNING: anchor '{INSERT_ANCHOR.strip()}' not found in {src.name}")

    dst.write_text(text, encoding="utf-8")
    print(f"Written: {dst.relative_to(Path('plugins').parent)}")


base = Path("plugins")

make_truststore(
    src=base / "Chat_with_Molecule_Neo/chat_with_molecule_neo.py",
    dst=base / "Chat_with_Molecule_Neo/chat_with_molecule_neo_truststore.py",
    name_suffix=" (truststore)",
    id_suffix="_truststore",
    add_docstring=False,
)

make_truststore(
    src=base / "Chat_with_Molecule_Neo_Local/chat_with_molecule_neo_local.py",
    dst=base / "Chat_with_Molecule_Neo_Local/chat_with_molecule_neo_local_truststore.py",
    name_suffix=" (truststore)",
    id_suffix="_truststore",
    add_docstring=True,
)

print("Done.")
