try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:  # pragma: no cover - Python <3.8 fallback
    from importlib_metadata import PackageNotFoundError, version

try:
    __version__ = version("mosaic")
except PackageNotFoundError:
    __version__ = "0+local"
