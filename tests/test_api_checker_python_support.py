"""The developer API checker must start on the plugin suite's interpreter."""

import subprocess
import sys
from pathlib import Path


def test_api_checker_cli_can_start():
    checker = Path(__file__).resolve().parents[1] / "api-checker/plugin_api_checker.py"
    result = subprocess.run(
        [sys.executable, str(checker), "--help"],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "--plugin" in result.stdout
