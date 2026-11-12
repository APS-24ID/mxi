"""The version: stamped in where CMake copies the package beside mxi; from the
wheel -- mxi, or mxi-cuda -- the wheel's own."""

VERSION = "@PROJECT_VERSION@"
if VERSION.startswith("@"):
    try:
        from importlib.metadata import version as _version

        try:
            VERSION = _version("mxi")
        except Exception:  # noqa: BLE001 -- the CUDA build's wheel
            VERSION = _version("mxi-cuda")
    except Exception:  # noqa: BLE001 -- run from a checkout, not installed
        VERSION = "unknown"
