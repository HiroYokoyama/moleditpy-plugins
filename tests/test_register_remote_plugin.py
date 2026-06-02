import json
import re
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

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
