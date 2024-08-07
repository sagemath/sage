# sage_setup: distribution = sagemath-repl
r"""
The ``pretty_print`` command

Works similar to the ``print`` function, except that it always tries
to use a rich output for an object, as specified via the text display
preference. If such a rich output is not available, it falls back on the
plain text.

EXAMPLES::

    sage: pretty_print(1, 2, 3)
    1 2 3

    sage: pretty_print(x^2 / (x + 1))                                                   # needs sage.symbolic
    x^2/(x + 1)

TESTS::

    sage: dm = get_display_manager()
    sage: dm.preferences.text = 'ascii_art'

EXAMPLES::

    sage: %display ascii_art  # not tested
    sage: pretty_print(x^2 / (x + 1))                                                   # needs sage.symbolic
       2
      x
    -----
    x + 1

TESTS:

After the previous example, we need to reset the text display preferences::

    sage: dm.preferences.text = None

EXAMPLES:

Printing a graphics object just prints a string, whereas
:func:`pretty_print` does not print anything and just shows the
graphics instead::

    sage: print(plot(sin))                                                              # needs sage.plot sage.symbolic
    Graphics object consisting of 1 graphics primitive
    sage: pretty_print(plot(sin))                                                       # needs sage.plot sage.symbolic
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.sage_object import SageObject
from sage.repl.rich_output import get_display_manager


class SequencePrettyPrinter(SageObject):

    def __init__(self, *args, **kwds):
        r"""
        Pretty Printer for Multiple Arguments.

        INPUT/OUTPUT:

        Same as :func:`pretty_print`, except that the number of
        arguments must be >= 2. Otherwise its not a sequence of things
        to print.

        EXAMPLES::

            sage: pretty_print(1, 2, 3)   # indirect doctest
            1 2 3
            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: SequencePrettyPrinter(1, 2, 3).pretty_print()
            1 2 3
        """
        self.args = args
        assert len(self.args) >= 2
        self.kwds = kwds

    def is_homogeneous(self, common_type):
        """
        Return whether the pretty print items are homogeneous.

        INPUT:

        - ``common_type`` -- a type

        OUTPUT:

        boolean; whether all items to be pretty printed are of said
        type.

        EXAMPLES::

            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: seq = SequencePrettyPrinter(1, 2, 3)
            sage: seq.is_homogeneous(Integer)
            True
            sage: seq.is_homogeneous(str)
            False
        """
        return all(isinstance(arg, common_type) for arg in self.args)

    def _concatenate_graphs(self):
        """
        Plot multiple graphs into a single plot.

        OUTPUT: a graphics object

        EXAMPLES::

            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: plt = SequencePrettyPrinter(*list(graphs(3)))._concatenate_graphs()   # needs sage.graphs sage.plot
            sage: type(plt)                                                             # needs sage.graphs sage.plot
            <class 'sage.plot.multigraphics.GraphicsArray'>
            sage: plt                                                                   # needs sage.graphs sage.plot
            Graphics Array of size 1 x 4
        """
        import sage.graphs.graph_list as graphs_list
        return graphs_list.to_graphics_array(self.args, **self.kwds)

    def _concatenate_graphics(self):
        """
        Combine multiple graphics objects into one graphics array.

        OUTPUT: a graphics array

        EXAMPLES::

            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: ga = SequencePrettyPrinter(*[Graphics()]*5)._concatenate_graphics()   # needs sage.plot
            sage: type(ga)                                                              # needs sage.plot
            <class 'sage.plot.multigraphics.GraphicsArray'>
            sage: ga.nrows(), ga.ncols()                                                # needs sage.plot
            (2, 4)
        """
        from sage.plot.plot import graphics_array
        return graphics_array(self.args, ncols=4, **self.kwds)

    def pretty_print(self):
        """
        Actually do the pretty print.

        EXAMPLES::

            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: SequencePrettyPrinter(1, 2, 3).pretty_print()
            1 2 3

        The keyword arguments are only used the first time graphics
        output is generated::

            sage: seq = SequencePrettyPrinter(Graph(), Graph(), edge_labels=True)       # needs sage.graphs sage.plot
            sage: seq.pretty_print()   # does not pass edge_labels to graphics object   # needs sage.graphs sage.plot
            sage: seq._concatenate_graphs().show(edge_labels=True)                      # needs sage.graphs sage.plot
            Traceback (most recent call last):
            ...
            TypeError: ...matplotlib() got an unexpected keyword argument 'edge_labels'
        """
        try:
            from sage.plot.plot import Graphics
        except ImportError:
            Graphics = ()
        try:
            from sage.graphs.graph import GenericGraph
        except ImportError:
            GenericGraph = ()
        if self.is_homogeneous(GenericGraph):
            args = self._concatenate_graphs()
            kwds = dict()
        elif self.is_homogeneous(Graphics):
            args = self._concatenate_graphics()
            kwds = dict()
        else:
            args = self.args
            kwds = dict(self.kwds)
            kwds['concatenate'] = True
        get_display_manager().display_immediately(args, **kwds)


def pretty_print(*args, **kwds):
    r"""
    Pretty print the arguments using rich output if available.

    This function is similar to ``print()``, except that a rich output
    representation such as ``ascii_art`` or Latex is printed instead of the
    string representation.

    Note that the output depends on the global display preferences specified
    via
    :meth:`~sage.repl.rich_output.display_manager.DisplayManager.preferences`.
    If the display preference for ``text`` is not specified, Latex output is
    preferred.

    For graphical objects, a graphical output is used.

    For certain homogeneous multiple positional arguments, a suitable
    combined graphical output is generated. In particular, graphs and
    plots are treated special. Otherwise this function concatenates the
    textual representations.

    INPUT:

    - ``*args`` -- any number of positional arguments; the objects to
      pretty print

    - ``**kwds`` -- optional keyword arguments that are passed to the
      rich representation. Examples include:

        - ``dpi`` -- dots per inch

        - ``figsize``- [width, height] (same for square aspect)

        - ``axes`` -- (default: ``True``)

        - ``fontsize`` -- positive integer

        - ``frame`` -- boolean (default: ``False``); draw a MATLAB-like frame around
          the image

    EXAMPLES::

        sage: pretty_print(ZZ)
        Integer Ring

        sage: pretty_print("Integers = ", ZZ) # trac 11775
        'Integers = ' Integer Ring

    To typeset LaTeX code as-is, use :class:`LatexExpr`::

        sage: pretty_print(LatexExpr(r"\frac{x^2 + 1}{x - 2}"))
        \frac{x^2 + 1}{x - 2}

    For text-based backends, the default text display preference is to output
    plain text which is usually the same as using ``print()``::

        sage: pretty_print(x^2 / (x + 1))                                               # needs sage.symbolic
        x^2/(x + 1)

        sage: t = BinaryTrees(3).first()                                                # needs sage.graphs
        sage: pretty_print(t)                                                           # needs sage.graphs
        [., [., [., .]]]
        sage: print(t)                                                                  # needs sage.graphs
        [., [., [., .]]]

    TESTS::

        sage: dm = get_display_manager()
        sage: dm.preferences.text = 'ascii_art'

    EXAMPLES:

    Changing the text display preference affects the output of this function.
    The following illustrates a possible use-case::

        sage: %display ascii_art  # not tested
        sage: for t in BinaryTrees(3)[:3]:                                              # needs sage.graphs
        ....:     pretty_print(t)
        o
         \
          o
           \
            o
        o
         \
          o
         /
        o
          o
         / \
        o   o

        sage: pretty_print(x^2 / (x + 1))                                               # needs sage.symbolic
           2
          x
        -----
        x + 1

    TESTS:

    After the previous example, we need to reset the text display preferences::

        sage: dm.preferences.text = None

    ::

        sage: # needs sage.plot sage.symbolic
        sage: plt = plot(sin)
        sage: pretty_print(plt)             # graphics output
        sage: pretty_print(plt, plt)        # graphics output
        sage: pretty_print(ZZ, 123, plt)
        Integer Ring 123 Graphics object consisting of 1 graphics primitive
    """
    dm = get_display_manager()
    old_preferences_text = dm.preferences.text
    try:
        if dm.preferences.text is None:
            dm.preferences.text = 'latex'
        if len(args) == 0:
            pass
        elif len(args) == 1:
            dm.display_immediately(*args, **kwds)
        else:
            SequencePrettyPrinter(*args, **kwds).pretty_print()
    finally:
        dm.preferences.text = old_preferences_text


def show(*args, **kwds):
    r"""
    Alias for ``pretty_print``.

    This function is an alias for :func:`pretty_print`.

    INPUT/OUTPUT:

    See :func:`pretty_print`. Except if the argument is a graph, in
    which case it is plotted instead.

    EXAMPLES::

        sage: show(1)
        1
    """
    try:
        from sage.graphs.generic_graph import GenericGraph
    except ImportError:
        GenericGraph = ()
    if len(args) == 1 and isinstance(args[0], GenericGraph):
        # Graphs are special, they ride the short bus...
        # Please, somebody help me get rid of this! #18289
        args[0].show()
        return
    pretty_print(*args, **kwds)
