try:
    from importlib import metadata
except ImportError as e:
    import importlib_metadata as metadata
#__version__ = '0.1-beta'
__version__ = metadata.version('mosaic')
