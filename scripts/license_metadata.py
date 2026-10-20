"""Shared handling for explicitly declared plugin license metadata."""

import re


def normalize_license(value):
    """Return a single-line SPDX-style value, or None for an absent value."""
    if value is None or value == "":
        return None
    if not isinstance(value, str):
        raise ValueError("Plugin license must be a string")
    value = value.strip()
    if not value:
        return None
    if len(value) > 256 or not re.fullmatch(r"[A-Za-z0-9(][A-Za-z0-9.+(): -]*", value):
        raise ValueError("Plugin license must be a single-line SPDX identifier or expression")
    return value


def resolve_license(code_value=None, input_value=None, existing_value=None):
    """Declared code wins; a maintainer input fills gaps; stored values survive."""
    for value in (code_value, input_value, existing_value):
        normalized = normalize_license(value)
        if normalized:
            return normalized
    return "NOASSERTION"
