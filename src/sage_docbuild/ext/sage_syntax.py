"""
Sphinx extension to translate Sage-specific markers to Sphinx field list entries.

Example::

    INPUT:
        - ``FOO`` -- Short description.

    OUTPUT:
        Result description.

    → becomes::

    :param FOO: Short description.
    :returns: Result description.
"""

from __future__ import annotations

from typing import List

from sphinx.application import Sphinx

from sage_docbuild.ext.sage_transformer import DoctestTransformer


def _process_docstring(
    app: Sphinx, what: str, name: str, obj: object, options: object, lines: List[str]
) -> None:
    """
    Process a docstring after autodoc has collected it.

        Called when autodoc has read and processed a docstring.

    INPUT:
        - app -- Sphinx application instance.
        - what -- Object type (module, class, exception, function, method, attribute).
        - name -- Fully qualified object name.
        - obj -- The object whose docstring is being processed.
        - options -- Autodoc directive options.
        - lines -- Mutable list of docstring lines (modified in place).
    """
    if not lines:
        return

    transformer = DoctestTransformer(lines, what, name, obj)
    lines[:] = transformer.transform()


def setup(app: Sphinx):
    app.setup_extension("sage_docbuild.ext.sage_autodoc")
    app.connect("autodoc-process-docstring", _process_docstring)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
