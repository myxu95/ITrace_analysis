"""
ImmunoScope CLI Module

Command-line interface for ImmunoScope tools.

Main entry point: ims (ImmunoScope)
"""

def imn_main(*args, **kwargs):
    """Lazy import main entry point to keep module execution warning-free."""
    from .main import main
    return main(*args, **kwargs)


def batch_pdb_main(*args, **kwargs):
    """Lazy import batch_pdb entry point to avoid unrelated import failures."""
    from .batch_pdb import main
    return main(*args, **kwargs)


__all__ = ['imn_main', 'batch_pdb_main']
