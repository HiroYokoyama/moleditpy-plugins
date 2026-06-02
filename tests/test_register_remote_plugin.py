import hashlib
import json
import re
import sys
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import pytest

# Add scripts directory to path to import register_remote_plugin
sys.path.append(str(Path(__file__).resolve().parents[1] / "scripts"))

import register_remote_plugin

# ---------------------------------------------------------------------------
# Unit tests for Helper Functions
# ---------------------------------------------------------------------------

def test_parse_github_release_url_valid():
    url = "https://github.com/HiroYokoyama/moleditpy_rotation_giffer/releases/download/1.2.0/rotation_giffer.py"
    info = register_remote_plugin.parse_github_release_url(url)
    assert info["owner"] == "HiroYokoyama"
    assert info["repo"] == "moleditpy_rotation_giffer"
    assert info["tag"] == "1.2.0"
    assert info["filename"] == "rotation_giffer.py"

    url_v = "https://github.com/owner-name/repo-name/releases/download/v2.4.5/plugin.zip"
    info_v = register_remote_plugin.parse_github_release_url(url_v)
    assert info_v["owner"] == "owner-name"
    assert info_v["repo"] == "repo-name"
    assert info_v["tag"] == "v2.4.5"
    assert info_v["filename"] == "plugin.zip"

def test_parse_github_release_url_invalid():
    invalid_urls = [
        "https://github.com/owner/repo/releases/tag/v1.0.0",
        "https://github.com/owner/repo/archive/refs/tags/v1.0.0.zip",
        "http://github.com/owner/repo/releases/download/v1.0.0/file.py",
        "https://gitlab.com/owner/repo/releases/download/v1.0.0/file.py"
    ]
    for url in invalid_urls:
        with pytest.raises(ValueError):
            register_remote_plugin.parse_github_release_url(url)

def test_parse_version():
    assert register_remote_plugin.parse_version("1.2.3") == (1, 2, 3)
    assert register_remote_plugin.parse_version("v1.2.3") == (1, 2, 3)
    assert register_remote_plugin.parse_version("V2.0.0") == (2, 0, 0)
    assert register_remote_plugin.parse_version("2.1") == (2, 1)
    assert register_remote_plugin.parse_version("1.2.3b") == (1, 2, "3b")
    assert register_remote_plugin.parse_version("2.3.0") > register_remote_plugin.parse_version("2.2.9")
    assert register_remote_plugin.parse_version("10.0.0") > register_remote_plugin.parse_version("2.0.0")

def test_extract_metadata_from_code():
    code = """
    # My Plugin
    PLUGIN_NAME = "Vibrational Analyzer"
    PLUGIN_VERSION = "2.5.1"
    PLUGIN_AUTHOR = "Jane Doe"
    PLUGIN_DESCRIPTION = "Analyzes vibrational normal modes."
    
    def initialize(context):
        pass
    """
    meta = register_remote_plugin.extract_metadata_from_code(code)
    assert meta["name"] == "Vibrational Analyzer"
    assert meta["version"] == "2.5.1"
    assert meta["author"] == "Jane Doe"
    assert meta["description"] == "Analyzes vibrational normal modes."

def test_extract_metadata_from_code_missing():
    code = """
    # No metadata here
    PLUGIN_NAME = "Half-defined Plugin"
    """
    meta = register_remote_plugin.extract_metadata_from_code(code)
    assert meta["name"] == "Half-defined Plugin"
    assert "version" not in meta
    assert "author" not in meta
    assert "description" not in meta

# ---------------------------------------------------------------------------
# Integration/Functional Tests for Registry Matching
# ---------------------------------------------------------------------------

def test_find_existing_plugin():
    plugins = [
        {
            "id": "rotation_giffer",
            "name": "Rotation Giffer",
            "version": "1.2.0",
            "downloadUrl": "https://github.com/HiroYokoyama/moleditpy_rotation_giffer/releases/download/1.2.0/rotation_giffer.py",
            "projectUrl": "https://github.com/HiroYokoyama/moleditpy_rotation_giffer"
        },
        {
            "id": "cif_viewer",
            "name": "CIF Viewer",
            "version": "0.9.0",
            "downloadUrl": "https://github.com/HiroYokoyama/moleditpy_cif_viewer/releases/download/0.9.0/cif_viewer.zip",
            "projectUrl": "https://github.com/HiroYokoyama/moleditpy_cif_viewer"
        }
    ]

    # Match by ID
    p_id = register_remote_plugin.find_existing_plugin(plugins, "owner", "repo", "file.py", "cif_viewer")
    assert p_id["id"] == "cif_viewer"

    # Match by repo URL and filename (updating same file)
    p_url = register_remote_plugin.find_existing_plugin(
        plugins, "HiroYokoyama", "moleditpy_rotation_giffer", "rotation_giffer.py"
    )
    assert p_url["id"] == "rotation_giffer"

    # Match by project URL
    p_proj = register_remote_plugin.find_existing_plugin(
        plugins, "HiroYokoyama", "moleditpy_cif_viewer", "cif_viewer_v2.zip"
    )
    assert p_proj["id"] == "cif_viewer"

    # No match
    p_none = register_remote_plugin.find_existing_plugin(
        plugins, "OtherOwner", "other_repo", "other_file.py"
    )
    assert p_none is None

# ---------------------------------------------------------------------------
# Unit tests for Security SHA-256 Checks
# ---------------------------------------------------------------------------

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data="[]")
def test_security_check_non_hiro_missing_sha(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Third Party Plugin",
        "version": "1.0.0",
        "author": "ThirdParty",
        "description": "Some description"
    }
    
    # Configure command-line args for non-HiroYokoyama repo, missing expected SHA
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/ThirdParty/some_plugin/releases/download/v1.0.0/plugin.py",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args):
        register_remote_plugin.main()
        
    # Expecting sys.exit(1) to be called due to missing SHA
    mock_exit.assert_called_with(1)

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data="[]")
def test_security_check_non_hiro_mismatched_sha(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Third Party Plugin",
        "version": "1.0.0",
        "author": "ThirdParty",
        "description": "Some description"
    }
    
    # Configure command-line args for non-HiroYokoyama repo, mismatched expected SHA
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/ThirdParty/some_plugin/releases/download/v1.0.0/plugin.py",
        "--expected-sha256", "wrongsha1234567890",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args):
        register_remote_plugin.main()
        
    # Expecting sys.exit(1) to be called due to mismatched SHA
    mock_exit.assert_called_with(1)

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data="[]")
def test_security_check_non_hiro_matching_sha(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_content = b"mock content"
    mock_response.read.return_value = mock_content
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    actual_sha = hashlib.sha256(mock_content).hexdigest()
    
    mock_extract.return_value = {
        "name": "Third Party Plugin",
        "version": "1.0.0",
        "author": "ThirdParty",
        "description": "Some description"
    }
    
    # Configure command-line args for non-HiroYokoyama repo, matching expected SHA
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/ThirdParty/some_plugin/releases/download/v1.0.0/plugin.py",
        "--expected-sha256", actual_sha,
        "--supported-version", "3.*",
        "--dry-run"
    ]
    
    # We expect sys.exit to NOT be called since verification succeeds
    with patch('sys.argv', test_args):
        register_remote_plugin.main()
        
    # Expecting sys.exit to NOT be called with 1 (success)
    calls = [c[0][0] for c in mock_exit.call_args_list if c[0]]
    assert 1 not in calls

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data='[{"id": "rotation_giffer", "name": "Rotation Giffer", "version": "1.2.0", "downloadUrl": "https://github.com/HiroYokoyama/moleditpy_rotation_giffer/releases/download/1.2.0/rotation_giffer.py", "projectUrl": "https://github.com/HiroYokoyama/moleditpy_rotation_giffer"}]')
def test_same_version_dry_run_warning(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Rotation Giffer",
        "version": "1.2.0",
        "author": "HiroYokoyama",
        "description": "Creates a rotating GIF around global or view axes by orbiting the camera."
    }
    
    # Configure command-line args for hiroyokoyama repo, same version, with dry-run
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/moleditpy_rotation_giffer/releases/download/1.2.0/rotation_giffer.py",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args):
        register_remote_plugin.main()
        
    # Check that sys.exit was NOT called with 1
    calls = [c[0][0] for c in mock_exit.call_args_list if c[0]]
    assert 1 not in calls

def test_parse_python_list_or_string():
    assert register_remote_plugin.parse_python_list_or_string('["A", "B"]') == ["A", "B"]
    assert register_remote_plugin.parse_python_list_or_string("('C', 'D')") == ["C", "D"]
    assert register_remote_plugin.parse_python_list_or_string('"E, F"') == ["E", "F"]
    assert register_remote_plugin.parse_python_list_or_string('"G"') == ["G"]
    assert register_remote_plugin.parse_python_list_or_string("G, H") == ["G", "H"]

def test_extract_metadata_with_tags_and_dependencies():
    code = """
    PLUGIN_NAME = "Advanced Tool"
    PLUGIN_VERSION = "1.0.0"
    PLUGIN_AUTHOR = "Author Name"
    PLUGIN_DESCRIPTION = "Some desc"
    PLUGIN_TAGS = ["Visualization", "Utility"]
    PLUGIN_DEPENDENCIES = "numpy, rdkit"
    """
    meta = register_remote_plugin.extract_metadata_from_code(code)
    assert meta["tags"] == ["Visualization", "Utility"]
    assert meta["dependencies"] == ["numpy", "rdkit"]

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data="[]")
def test_author_mismatch_fails(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Third Party Plugin",
        "version": "1.0.0",
        "author": "Real Name Instead of GitHub Username",
        "description": "Some description"
    }
    
    # Configure command-line args for non-HiroYokoyama repo, mismatched author
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/ThirdParty/some_plugin/releases/download/v1.0.0/plugin.py",
        "--expected-sha256", hashlib.sha256(b"mock content").hexdigest(),
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args):
        register_remote_plugin.main()
        
    # Expecting sys.exit(1) to be called due to mismatched author
    mock_exit.assert_called_with(1)

def test_extract_metadata_multiline_parenthesized():
    code = """
    PLUGIN_NAME = "PySCF Calculator"
    PLUGIN_VERSION = "2.3.0"
    PLUGIN_AUTHOR = "HiroYokoyama"
    PLUGIN_DESCRIPTION = (
        "Perform PySCF quantum chemistry calculations directly in MoleditPy. "
        "Features: Single Point Energy (RHF/UHF/DFT), Geometry Optimization (Geometric/Berny), "
        "Frequency Analysis, and interactive 3D visualization."
    )
    """
    meta = register_remote_plugin.extract_metadata_from_code(code)
    assert meta["name"] == "PySCF Calculator"
    assert meta["version"] == "2.3.0"
    assert meta["author"] == "HiroYokoyama"
    assert "Perform PySCF quantum chemistry" in meta["description"]
    assert "interactive 3D visualization" in meta["description"]

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data="[]")
@patch('builtins.print')
def test_dry_run_prints_new_entry(mock_print, mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Dry Run Test Plugin",
        "version": "1.0.0",
        "author": "HiroYokoyama",
        "description": "Prints entry in dry run mode"
    }
    
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/dry_run_test_plugin/releases/download/v1.0.0/dry_run_test_plugin.py",
        "--supported-version", "3.*",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args):
        register_remote_plugin.main()
        
    # Verify that mock_print was called to output the JSON structure of the new entry
    printed_calls = [c[0][0] for c in mock_print.call_args_list if c[0]]
    # Find if any call has the JSON content we expect
    found_json = False
    for p in printed_calls:
        if isinstance(p, str) and '"id": "dry_run_test_plugin"' in p and '"supported_moleditpy_version": "3.*"' in p:
            found_json = True
            break
    assert found_json, "Expected the dry run output to contain the proposed plugin entry JSON with supported_moleditpy_version"

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data='[]')
def test_date_override_valid(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Date Override Plugin",
        "version": "1.0.0",
        "author": "HiroYokoyama",
        "description": "Some description"
    }
    
    # Configure command-line args with valid date override
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.0.0/some_plugin.py",
        "--date", "2026-05-20",
        "--supported-version", "3.*",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args), patch('builtins.print') as mock_print:
        register_remote_plugin.main()
        
    # Verify that the printed output JSON contains the overriden date in lastUpdated/firstAppeared
    printed_calls = [c[0][0] for c in mock_print.call_args_list if c[0]]
    found_date = False
    for p in printed_calls:
        if isinstance(p, str) and '"lastUpdated": "2026-05-20"' in p and '"firstAppeared": "2026-05-20"' in p:
            found_date = True
            break
    assert found_date, "Expected the dry run output to contain the overridden date"

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data='[]')
def test_date_override_invalid(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Date Override Plugin",
        "version": "1.0.0",
        "author": "HiroYokoyama",
        "description": "Some description"
    }
    
    # Configure command-line args with invalid date override
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.0.0/some_plugin.py",
        "--date", "2026/05/20",
        "--supported-version", "3.*",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args):
        register_remote_plugin.main()
        
    mock_exit.assert_called_with(1)

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data='[{"id": "some_plugin", "visible": true, "supported_moleditpy_version": "3.5", "name": "Some Plugin", "version": "1.0.0", "downloadUrl": "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.0.0/some_plugin.py", "projectUrl": "https://github.com/HiroYokoyama/some_plugin"}]')
def test_update_preserves_previous_supported_version(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Some Plugin",
        "version": "1.1.0",
        "author": "HiroYokoyama",
        "description": "Some description"
    }
    
    # Configure command-line args for update without passing --supported-version
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.1.0/some_plugin.py",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args), patch('builtins.print') as mock_print:
        register_remote_plugin.main()
        
    # Verify that the output JSON still contains `"supported_moleditpy_version": "3.5"`
    printed_calls = [c[0][0] for c in mock_print.call_args_list if c[0]]
    found_preserved = False
    for p in printed_calls:
        if isinstance(p, str) and '"supported_moleditpy_version": "3.5"' in p:
            found_preserved = True
            break
    assert found_preserved, "Expected the update output to preserve the existing supported version"

@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data='[{"id": "some_plugin", "visible": true, "supported_moleditpy_version": "3.5", "name": "Some Plugin", "version": "1.0.0", "downloadUrl": "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.0.0/some_plugin.py", "projectUrl": "https://github.com/HiroYokoyama/some_plugin"}]')
def test_update_overwrites_supported_version(mock_file, mock_extract, mock_urlopen, mock_exit):
    # Setup mocks
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Some Plugin",
        "version": "1.1.0",
        "author": "HiroYokoyama",
        "description": "Some description"
    }
    
    # Configure command-line args for update passing a new --supported-version
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.1.0/some_plugin.py",
        "--supported-version", "3.6",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args), patch('builtins.print') as mock_print:
        register_remote_plugin.main()
        
    # Verify that the output JSON contains the new `"supported_moleditpy_version": "3.6"`
    printed_calls = [c[0][0] for c in mock_print.call_args_list if c[0]]
    found_overwritten = False
    for p in printed_calls:
        if isinstance(p, str) and '"supported_moleditpy_version": "3.6"' in p:
            found_overwritten = True
            break
    assert found_overwritten, "Expected the update output to contain the overwritten supported version"


def test_extract_metadata_with_supported_version():
    code = """
    PLUGIN_NAME = "Supported Plugin"
    PLUGIN_VERSION = "1.0.0"
    PLUGIN_AUTHOR = "Author Name"
    PLUGIN_DESCRIPTION = "Some desc"
    PLUGIN_SUPPORTED_MOLEDITPY_VERSION = "3.2"
    """
    meta = register_remote_plugin.extract_metadata_from_code(code)
    assert meta["supported_moleditpy_version"] == "3.2"


@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data="[]")
def test_add_prioritizes_cli_version_over_code(mock_file, mock_extract, mock_urlopen, mock_exit):
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "CLI Priority Plugin",
        "version": "1.0.0",
        "author": "HiroYokoyama",
        "description": "Desc",
        "supported_moleditpy_version": "3.1"
    }
    
    # Passing --supported-version "3.5" (CLI)
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/cli_priority/releases/download/v1.0.0/plugin.py",
        "--supported-version", "3.5",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args), patch('builtins.print') as mock_print:
        register_remote_plugin.main()
        
    # Verify that the printed output JSON contains `"supported_moleditpy_version": "3.5"` (CLI wins)
    printed_calls = [c[0][0] for c in mock_print.call_args_list if c[0]]
    found_cli = False
    for p in printed_calls:
        if isinstance(p, str) and '"supported_moleditpy_version": "3.5"' in p:
            found_cli = True
            break
    assert found_cli, "Expected CLI version to override code constant"


@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data="[]")
def test_add_uses_code_constant_if_cli_omitted(mock_file, mock_extract, mock_urlopen, mock_exit):
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    mock_extract.return_value = {
        "name": "Code Priority Plugin",
        "version": "1.0.0",
        "author": "HiroYokoyama",
        "description": "Desc",
        "supported_moleditpy_version": "3.1"
    }
    
    # Omit --supported-version
    test_args = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/code_priority/releases/download/v1.0.0/plugin.py",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args), patch('builtins.print') as mock_print:
        register_remote_plugin.main()
        
    # Verify that the printed output JSON contains `"supported_moleditpy_version": "3.1"` (Code wins)
    printed_calls = [c[0][0] for c in mock_print.call_args_list if c[0]]
    found_code = False
    for p in printed_calls:
        if isinstance(p, str) and '"supported_moleditpy_version": "3.1"' in p:
            found_code = True
            break
    assert found_code, "Expected code constant to be used when CLI is omitted"


@patch('sys.exit')
@patch('urllib.request.urlopen')
@patch('register_remote_plugin.extract_metadata_from_file')
@patch('builtins.open', new_callable=mock_open, read_data='[{"id": "some_plugin", "visible": true, "supported_moleditpy_version": "3.0", "name": "Some Plugin", "version": "1.0.0", "downloadUrl": "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.0.0/some_plugin.py", "projectUrl": "https://github.com/HiroYokoyama/some_plugin"}]')
def test_update_prioritizes_cli_then_code_then_registry(mock_file, mock_extract, mock_urlopen, mock_exit):
    mock_response = MagicMock()
    mock_response.read.return_value = b"mock content"
    mock_urlopen.return_value.__enter__.return_value = mock_response
    
    # Case A: CLI provided, code constant provided, registry has value. CLI should win (3.5)
    mock_extract.return_value = {
        "name": "Some Plugin",
        "version": "1.1.0",
        "author": "HiroYokoyama",
        "description": "Some description",
        "supported_moleditpy_version": "3.2"
    }
    
    test_args_cli = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.1.0/some_plugin.py",
        "--supported-version", "3.5",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args_cli), patch('builtins.print') as mock_print_cli:
        register_remote_plugin.main()
        
    printed_calls = [c[0][0] for c in mock_print_cli.call_args_list if c[0]]
    found_cli = False
    for p in printed_calls:
        if isinstance(p, str) and '"supported_moleditpy_version": "3.5"' in p:
            found_cli = True
            break
    assert found_cli, "Expected CLI version to win in UPDATE"

    # Case B: CLI omitted, code constant provided, registry has value. Code constant should win (3.2)
    test_args_code = [
        "scripts/register_remote_plugin.py",
        "https://github.com/HiroYokoyama/some_plugin/releases/download/v1.1.0/some_plugin.py",
        "--dry-run"
    ]
    
    with patch('sys.argv', test_args_code), patch('builtins.print') as mock_print_code:
        register_remote_plugin.main()
        
    printed_calls = [c[0][0] for c in mock_print_code.call_args_list if c[0]]
    found_code = False
    for p in printed_calls:
        if isinstance(p, str) and '"supported_moleditpy_version": "3.2"' in p:
            found_code = True
            break
    assert found_code, "Expected code constant to win when CLI omitted in UPDATE"





