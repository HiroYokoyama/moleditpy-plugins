#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import base64
import tempfile
import inspect

from PyQt6.QtWidgets import QInputDialog, QLineEdit, QMessageBox, QFileDialog
from PyQt6.QtGui import QAction
from PyQt6.QtCore import QTimer
import logging

# Cryptography imports
try:
    from cryptography.fernet import Fernet, InvalidToken
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC

    CRYPTOGRAPHY_AVAILABLE = True
except ImportError:
    CRYPTOGRAPHY_AVAILABLE = False

PLUGIN_NAME = "Encrypted Project"
PLUGIN_VERSION = "2026.10.13"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DEPENDENCIES = ["cryptography"]
MAX_PROJECT_BYTES = 64 * 1024 * 1024


def _call_original(method, args, kwargs):
    """Adapt QAction checked arguments before calling, without retrying failures."""
    signature = inspect.signature(method)
    try:
        signature.bind(*args, **kwargs)
    except TypeError:
        if len(args) == 1 and isinstance(args[0], bool) and not kwargs:
            signature.bind()
            return method()
        raise
    return method(*args, **kwargs)


PLUGIN_DESCRIPTION = "Securely saves molecular data using AES-128 encryption with password protection (.pmeenc file)."


class PmeencPlugin:
    def __init__(self, context):
        self.context = context
        self.mw = context.get_main_window()
        self._original_save_project = None
        self._original_clear_all = None
        self.current_key = None  # Store the derived key
        self.current_salt = None  # Store the salt used for the key

    def initialize(self):
        if not CRYPTOGRAPHY_AVAILABLE:
            print("PmeencPlugin: cryptography library not found. Plugin disabled.")
            return

        # Register file opener for .pmeenc
        self.context.register_file_opener(".pmeenc", self.on_import)

        # Register export action
        self.context.add_export_action("Export Encrypted (.pmeenc)...", self.on_export)

        # Register document reset handler
        self.context.register_document_reset_handler(self.on_document_reset)

        # Register drop handler
        self.context.register_drop_handler(self.on_drop)

    def on_drop(self, file_path: str) -> bool:
        """Handle file drops. Returns True if handled."""
        if file_path.lower().endswith(".pmeenc"):
            # Use QTimer to avoid blocking the drop event's GUI thread with a dialog
            def safe_import():
                try:
                    self.on_import(file_path)
                except ValueError as _e:
                    # Silently handle cancelled or failed authentication from UI
                    logging.warning("[encrypted_project.py:60] silenced: %s", _e)
                except Exception as e:
                    logging.exception("Unexpected import error: %s", e)
                    QMessageBox.critical(
                        self.mw, "Import Error", f"Unexpected error: {e}"
                    )

            QTimer.singleShot(0, safe_import)
            return True
        return False

    def _apply_patches(self):
        """Apply monkey patches and reconnect signals."""
        self._patch_method(
            "save_project", self._patched_save_project, ["Save Project", "Ctrl+S"]
        )
        self._patch_method(
            "clear_all",
            self._patched_clear_all,
            ["Clear All", "New", "Ctrl+Shift+C", "Ctrl+N"],
        )

    def _patch_method(self, method_name, patch_func, search_terms):
        if not hasattr(self.mw, method_name):
            return

        original_method = getattr(self.mw, method_name)
        # Avoid patching our own patches
        if getattr(original_method, "_is_pmeenc_patch", False):
            return

        # Store original for the plugin instance if not already stored
        attr_name = f"_original_{method_name}"
        if getattr(self, attr_name, None) is None:
            setattr(self, attr_name, original_method)

        # Create a wrapper function that allows attributes and forwards args
        def new_method(*args, **kwargs):
            return patch_func(self.mw, *args, **kwargs)

        new_method._is_pmeenc_patch = True
        new_method._original = original_method

        # Replace on instance
        setattr(self.mw, method_name, new_method)

        # Fix signal connections for QActions
        for action in self.mw.findChildren(QAction):
            text = action.text().replace("&", "")
            shortcut = action.shortcut().toString()

            should_reconnect = False
            for term in search_terms:
                if term in text or term == shortcut:
                    should_reconnect = True
                    break

            if should_reconnect:
                try:
                    # Surgical disconnect: only disconnect the original method
                    action.triggered.disconnect(original_method)
                except (TypeError, RuntimeError) as _e:
                    # Fallback for complex connections or if already disconnected
                    logging.warning("[encrypted_project.py:115] silenced: %s", _e)
                action.triggered.connect(new_method)

    def _unpatch_all(self):
        """Revert all monkey patches."""
        self._unpatch_method("save_project", ["Save Project", "Ctrl+S"])
        self._unpatch_method(
            "clear_all", ["Clear All", "New", "Ctrl+Shift+C", "Ctrl+N"]
        )

    def _unpatch_method(self, method_name, search_terms):
        attr_name = f"_original_{method_name}"
        original_method = getattr(self, attr_name, None)
        current_method = getattr(self.mw, method_name, None)

        if (
            original_method
            and current_method
            and getattr(current_method, "_is_pmeenc_patch", False)
        ):
            # Restore on instance
            setattr(self.mw, method_name, original_method)
            setattr(self, attr_name, None)

            # Reconnect signals back to original
            for action in self.mw.findChildren(QAction):
                text = action.text().replace("&", "")
                shortcut = action.shortcut().toString()

                should_reconnect = False
                for term in search_terms:
                    if term in text or term == shortcut:
                        should_reconnect = True
                        break

                if should_reconnect:
                    try:
                        # Surgical disconnect: only disconnect our patch
                        action.triggered.disconnect(current_method)
                    except (TypeError, RuntimeError) as _e:
                        logging.warning("[encrypted_project.py:150] silenced: %s", _e)
                    action.triggered.connect(original_method)

    def _patched_save_project(self, mw_instance, *args, **kwargs):
        """Monkey-patched version of save_project."""
        if (
            mw_instance.current_file_path
            and mw_instance.current_file_path.lower().endswith(".pmeenc")
        ):
            # For encrypted saving, we ignore external args
            self.export_encrypted(mw_instance.current_file_path)
        else:
            original = getattr(self, "_original_save_project", None)
            if original:
                return _call_original(original, args, kwargs)

    def _patched_clear_all(self, mw_instance, *args, **kwargs):
        """Monkey-patched version of clear_all."""
        original = getattr(self, "_original_clear_all", None)
        if original:
            return _call_original(original, args, kwargs)

    def on_document_reset(self):
        """Called when document is cleared or a new one is created."""
        self.current_key = None
        self.current_salt = None

        # Reset colors from Atom Colorizer plugin if present
        if hasattr(self.mw, "_plugin_color_overrides"):
            self.mw._plugin_color_overrides = {}
            # Trigger a redraw if there is a molecule
            mol = self.context.current_molecule
            if mol:
                self.context.draw_molecule_3d(mol)

        self._unpatch_all()

    def derive_key(self, password: str, salt: bytes) -> bytes:
        kdf = PBKDF2HMAC(
            algorithm=hashes.SHA256(),
            length=32,
            salt=salt,
            iterations=100000,
        )
        return base64.urlsafe_b64encode(kdf.derive(password.encode()))

    def on_export(self):
        """Export as .pmeenc file."""
        if not self.mw.state_manager.data.atoms and not self.context.current_molecule:
            self.context.show_status_message("Error: Nothing to save.")
            return

        default_name = "untitled"
        if self.mw.init_manager.current_file_path:
            base = os.path.basename(self.mw.init_manager.current_file_path)
            default_name = os.path.splitext(base)[0]

        default_path = default_name
        if self.mw.init_manager.current_file_path:
            default_path = os.path.join(
                os.path.dirname(self.mw.init_manager.current_file_path), default_name
            )

        file_path, _ = QFileDialog.getSaveFileName(
            self.mw,
            "Export Encrypted",
            default_path,
            "Encrypted Files (*.pmeenc);;All Files (*)",
        )

        if not file_path:
            return

        if not file_path.lower().endswith(".pmeenc"):
            file_path += ".pmeenc"

        self.export_encrypted(file_path)

    def export_encrypted(self, file_path: str):
        """The core encryption and saving logic."""
        # Check if we can do a seamless overwrite with cached key/salt
        if (
            self.current_key
            and self.current_salt
            and file_path == self.mw.init_manager.current_file_path
        ):
            key = self.current_key
            salt = self.current_salt
        else:
            # New file or no cached key: prompt for password
            prompt_text = (
                "Enter encryption password:"
                if not self.current_key
                else "Enter password for NEW file:"
            )
            while True:
                password, ok = QInputDialog.getText(
                    self.mw,
                    "Encryption Password",
                    prompt_text,
                    QLineEdit.EchoMode.Password,
                )

                if not ok:
                    self.context.show_status_message("Export cancelled.")
                    return

                if not password:
                    QMessageBox.warning(self.mw, "Error", "Password cannot be empty.")
                    continue

                # Derive new key with fresh salt
                salt = os.urandom(16)
                key = self.derive_key(password, salt)
                break

        try:
            # Prepare data
            json_data = self.mw.state_manager.create_json_data()
            raw_data = json.dumps(json_data, allow_nan=False).encode("utf-8")
            if len(raw_data) > MAX_PROJECT_BYTES:
                raise ValueError("Project exceeds 64 MiB.")

            # Encrypt
            f = Fernet(key)
            encrypted_data = f.encrypt(raw_data)

            # Store: [Salt(16)] + [EncryptedData]
            fd, staged = tempfile.mkstemp(
                dir=os.path.dirname(os.path.abspath(file_path)), suffix='.pmeenc.tmp'
            )
            try:
                with os.fdopen(fd, 'wb') as file:
                    file.write(salt)
                    file.write(encrypted_data)
                    file.flush()
                    os.fsync(file.fileno())
                os.replace(staged, file_path)
            finally:
                if os.path.exists(staged):
                    os.remove(staged)

            # Update app state
            self.current_key = key
            self.current_salt = salt
            self.mw.state_manager.has_unsaved_changes = False
            self.mw.init_manager.current_file_path = file_path
            self.context.refresh_ui()

            # Patch only after successful save
            self._apply_patches()

            self.context.show_status_message(f"Encrypted project saved to {file_path}")

        except Exception as e:
            logging.exception("Encryption error: %s", e)
            QMessageBox.critical(self.mw, "Encryption Error", f"Failed to encrypt: {e}")

    def on_import(self, file_path: str):
        """Import logic for .pmeenc files."""
        while True:
            password, ok = QInputDialog.getText(
                self.mw,
                "Decryption Password",
                "Enter password for decryption:",
                QLineEdit.EchoMode.Password,
            )

            if not ok:
                # User cancelled explicitly
                raise ValueError("Import cancelled by user.")

            if not password:
                QMessageBox.warning(self.mw, "Error", "Password cannot be empty.")
                continue

            try:
                with open(file_path, "rb") as file:
                    salt = file.read(16)
                    encrypted_data = file.read(MAX_PROJECT_BYTES * 2 + 1)
                if len(salt) != 16 or not encrypted_data:
                    raise ValueError("Encrypted project is truncated.")
                if len(encrypted_data) > MAX_PROJECT_BYTES * 2:
                    raise ValueError("Encrypted project exceeds the size limit.")

                key = self.derive_key(password, salt)
                f = Fernet(key)

                try:
                    decrypted_data = f.decrypt(encrypted_data)
                except InvalidToken:
                    QMessageBox.warning(
                        self.mw,
                        "Decryption Failed",
                        "Incorrect password or corrupted file.",
                    )
                    # Loop back to ask again
                    continue

                if len(decrypted_data) > MAX_PROJECT_BYTES:
                    raise ValueError("Decrypted project exceeds 64 MiB.")
                json_data = json.loads(decrypted_data.decode("utf-8"))
                if not isinstance(json_data, dict):
                    raise ValueError("Project data must be a JSON object.")

                self.mw.ui_manager.restore_ui_for_editing()
                self.mw.state_manager.load_from_json_data(json_data)

                # Reset app state to loaded file
                self.current_key = key
                self.current_salt = salt
                self.mw.state_manager.reset_undo_stack()
                self.mw.state_manager.has_unsaved_changes = False
                self.mw.init_manager.current_file_path = file_path
                self.context.refresh_ui()

                # Patch only after successful load
                self._apply_patches()

                self.context.show_status_message(
                    f"Encrypted project loaded from {file_path}"
                )
                QTimer.singleShot(0, self.context.fit_2d_view)

                # Success: break the retry loop
                break

            except Exception as e:
                logging.exception("Decryption error: %s", e)
                QMessageBox.critical(
                    self.mw, "Decryption Error", f"Failed to decrypt: {e}"
                )
                raise  # Re-raise to signal failure to load_command_line_file


def initialize(context):
    plugin = PmeencPlugin(context)
    plugin.initialize()
