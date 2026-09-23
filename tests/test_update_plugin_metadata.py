from __future__ import annotations

import json
import sys
from pathlib import Path

from unittest.mock import MagicMock, patch

import pytest

# Add scripts directory to path to import scripts
sys.path.append(str(Path(__file__).resolve().parents[1] / "scripts"))

import update_plugin_metadata


def test_update_metadata_tags_and_description():
    plugins = [
        {
            "id": "test_plugin",
            "visible": True,
            "supported_moleditpy_version": ">=4.0.0, <5.0.0",
            "name": "Old Plugin Name",
            "version": "1.0.0",
            "author": "HiroYokoyama",
            "authorUrl": "https://github.com/HiroYokoyama",
            "projectUrl": "https://github.com/HiroYokoyama/moleditpy_test",
            "description": "Old description",
            "tags": ["OldTag"],
            "dependencies": ["rdkit"],
            "downloadUrl": "https://github.com/HiroYokoyama/moleditpy_test/releases/download/v1.0.0/test_plugin.zip",
            "lastUpdated": "2026-01-01",
            "sha256": "fakehash",
            "firstAppeared": "2026-01-01",
            "supported_python_version": ">=3.9, <3.15",
            "supported_os": ["Windows", "macOS", "Linux", "WSL"],
        }
    ]

    mock_meta = {
        "name": "Updated Plugin Name",
        "description": "Fresh description from code",
        "tags": ["3D", "Editing", "Chemistry"],
        "dependencies": ["rdkit", "PyQt6"],
        "supported_moleditpy_version": ">=4.0.0, <5.0.0",
    }

    url = "https://github.com/HiroYokoyama/moleditpy_test/releases/download/v1.0.0/test_plugin.zip"

    with patch("update_plugin_metadata.extract_metadata_from_file", return_value=mock_meta), \
         patch("urllib.request.urlopen"), \
         patch("json.load", return_value=plugins), \
         patch("builtins.open"):

        with patch("pathlib.Path.exists", return_value=True):
            res = update_plugin_metadata.update_metadata_only(
                release_url=url,
                sync_all_from_code=True,
                dry_run=True,
            )

    assert res["plugin_id"] == "test_plugin"
    assert "tags (from code)" in res["changed_fields"]
    assert "description (from code)" in res["changed_fields"]
    assert "name (from code)" in res["changed_fields"]
    assert "dependencies (from code)" in res["changed_fields"]
    assert plugins[0]["tags"] == ["3D", "Editing", "Chemistry"]
    assert plugins[0]["description"] == "Fresh description from code"
    assert plugins[0]["dependencies"] == ["rdkit", "PyQt6"]
    # Version, URL, sha256 should remain untouched
    assert plugins[0]["version"] == "1.0.0"
    assert plugins[0]["sha256"] == "fakehash"


def test_update_metadata_code_wins_input_fills_gaps():
    plugins = [
        {
            "id": "test_plugin",
            "visible": True,
            "name": "Plugin",
            "version": "1.0.0",
            "tags": ["OldTag"],
            "downloadUrl": "https://github.com/HiroYokoyama/moleditpy_test/releases/download/v1.0.0/test_plugin.zip",
            "sha256": "fakehash",
        }
    ]

    mock_meta = {"tags": ["CodeTag"]}
    url = "https://github.com/HiroYokoyama/moleditpy_test/releases/download/v1.0.0/test_plugin.zip"

    with patch("update_plugin_metadata.extract_metadata_from_file", return_value=mock_meta), \
         patch("urllib.request.urlopen"), \
         patch("json.load", return_value=plugins), \
         patch("builtins.open"):

        with patch("pathlib.Path.exists", return_value=True):
            res = update_plugin_metadata.update_metadata_only(
                release_url=url,
                tags="ManualTag1, ManualTag2",
                description="Manual override desc",
                dry_run=True,
            )

    # The code declares PLUGIN_TAGS, so it wins over the input; the code has no
    # description, so the input fills that gap.
    assert "tags (from code)" in res["changed_fields"]
    assert plugins[0]["tags"] == ["CodeTag"]
    assert "description (from input)" in res["changed_fields"]
    assert plugins[0]["description"] == "Manual override desc"


def test_update_metadata_code_wins_for_every_declared_field():
    plugins = [
        {
            "id": "test_plugin",
            "visible": True,
            "supported_moleditpy_version": ">=4.0.0, <5.0.0",
            "name": "Plugin",
            "version": "1.0.0",
            "dependencies": ["numpy"],
            "downloadUrl": "https://github.com/HiroYokoyama/moleditpy_test/releases/download/v1.0.0/test_plugin.zip",
            "sha256": "fakehash",
            "supported_os": ["Windows", "macOS", "Linux", "WSL"],
        }
    ]
    mock_meta = {
        "dependencies": ["scipy"],
        "supported_moleditpy_version": ">=4.1.0, <5.0.0",
        "supported_os": ["macOS", "Linux"],
    }
    url = "https://github.com/HiroYokoyama/moleditpy_test/releases/download/v1.0.0/test_plugin.zip"

    with patch("update_plugin_metadata.extract_metadata_from_file", return_value=mock_meta), \
         patch("urllib.request.urlopen"), \
         patch("json.load", return_value=plugins), \
         patch("builtins.open"), \
         patch("pathlib.Path.exists", return_value=True):
        update_plugin_metadata.update_metadata_only(
            release_url=url,
            dependencies="pandas",
            supported_version=">=3.0.0",
            supported_os="Windows",
            dry_run=True,
        )

    assert plugins[0]["dependencies"] == ["scipy"]
    assert plugins[0]["supported_moleditpy_version"] == ">=4.1.0, <5.0.0"
    assert plugins[0]["supported_os"] == ["macOS", "Linux"]
