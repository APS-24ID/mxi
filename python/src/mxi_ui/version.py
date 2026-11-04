"""The version: stamped in where CMake copies the package beside mxi; from the
wheel, mxi-python, the wheel's own."""

VERSION = "@PROJECT_VERSION@"
if VERSION.startswith("@"):
    try:
        from importlib.metadata import version as _version

        VERSION = _version("mxi-python")
    except Exception:  # noqa: BLE001 -- run from a checkout, not installed
        VERSION = "unknown"
