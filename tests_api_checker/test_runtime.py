import ast
import importlib.util
import subprocess
import sys
from pathlib import Path

CHECKER = Path(__file__).resolve().parents[1] / "api-checker/plugin_api_checker.py"


def test_consumer_importlib_loader():
    spec = importlib.util.spec_from_file_location("checker_probe", CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.Issue("x.py", 1, "CODE", "message").line == 1
    assert module.APIInfo().mw_members == {}


def test_python39_and_cli():
    ast.parse(CHECKER.read_text(encoding="utf-8"), feature_version=(3, 9))
    result = subprocess.run(
        [sys.executable, str(CHECKER), "--help"],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert result.returncode == 0, result.stderr
