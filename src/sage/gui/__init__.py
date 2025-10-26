"""sage.gui

Top-level package for the Sage GUI extension.

This module exposes a single callable `main()` so the project can
register an entry point such as ``sage-gui = "sage.gui:main"``.
"""

from .app import main as _main


def main():
    """Start the GUI application.

    This function simply forwards to :func:`sage.gui.app.main`.
    It exists so the project's entry point can reference
    ``sage.gui:main`` directly.
    """
    return _main()