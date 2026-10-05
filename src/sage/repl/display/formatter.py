r"""
IPython Displayhook Formatters

The classes in this module can be used as IPython displayhook
formatters. It has two main features, by default the displayhook
contains a new facility for displaying lists of matrices in an easier
to read format::

    sage: [identity_matrix(i) for i in range(2, 5)]
    [
                    [1 0 0 0]
           [1 0 0]  [0 1 0 0]
    [1 0]  [0 1 0]  [0 0 1 0]
    [0 1], [0 0 1], [0 0 0 1]
    ]

This facility uses ``_repr_`` (and a simple string) to try do a nice read
format (see :meth:`sage.structure.parent.Parent._repr_option` for details).

With this displayhook there exists an other way for displaying object and more
generally, all sage expression as an ASCII art object::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell('%display ascii_art')
    sage: shell.run_cell('integral(x^2/pi^x, x)')
       -x / 2    2                      \
    -pi  *\x *log (pi) + 2*x*log(pi) + 2/
    --------------------------------------
                     3
                   log (pi)
    sage: shell.run_cell("i = var('i')")
    sage: shell.run_cell('sum(i*x^i, i, 0, 10)')
        10      9      8      7      6      5      4      3      2
    10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x
    sage: shell.run_cell('StandardTableaux(4).list()')
    [
    [                                                                  1  4    1  3
    [                 1  3  4    1  2  4    1  2  3    1  3    1  2    2       2
    [   1  2  3  4,   2      ,   3      ,   4      ,   2  4,   3  4,   3   ,   4   ,
    <BLANKLINE>
               1 ]
       1  2    2 ]
       3       3 ]
       4   ,   4 ]
    sage: shell.run_cell('%display default')
    sage: shell.quit()

This other facility uses a simple
:class:`~sage.typeset.ascii_art.AsciiArt` object (see and
:meth:`sage.structure.sage_object.SageObject._ascii_art_`).  """

# ****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from io import StringIO

from IPython.core.display import DisplayObject
from IPython.core.formatters import DisplayFormatter, PlainTextFormatter
from ipywidgets import Widget

from sage.features import FeatureNotPresentError
from sage.misc.lazy_import import LazyImport, lazy_import
from sage.repl.display.pretty_print import SagePrettyPrinter

IPYTHON_NATIVE_TYPES = (DisplayObject, Widget)

PLAIN_TEXT = 'text/plain'
TEXT_LATEX = 'text/latex'
TEXT_HTML = 'text/html'

lazy_import('matplotlib.figure', 'Figure')


def _resolve_lazy_import(obj):
    """
    Return the object behind a lazy import.

    Objects that are not lazy imports are returned unchanged.  A lazy import
    may resolve to another one, so resolution is repeated; a cyclic chain of
    lazy imports is left alone rather than followed forever.

    The lazy import itself is kept whenever it does not resolve to something
    displayable, so that the caller still sees the proxy and its established
    behaviour: :meth:`~sage.misc.lazy_import.LazyImport.__repr__` reports a
    missing feature instead of raising, and
    :meth:`~sage.repl.rich_output.display_manager.DisplayManager.displayhook`
    reads ``None`` as "nothing to display".

    EXAMPLES::

        sage: from sage.misc.lazy_import import LazyImport
        sage: from sage.repl.display.formatter import _resolve_lazy_import
        sage: _resolve_lazy_import(42)
        42
        sage: _resolve_lazy_import(LazyImport('sage.rings.integer_ring', 'ZZ'))
        Integer Ring

    Chains of lazy imports are resolved all the way down::

        sage: import sys, types
        sage: module = types.ModuleType('sage_doctest_lazy_module')
        sage: sys.modules['sage_doctest_lazy_module'] = module
        sage: module.inner = LazyImport('sage.rings.integer_ring', 'ZZ')
        sage: _resolve_lazy_import(LazyImport('sage_doctest_lazy_module', 'inner'))
        Integer Ring

    A lazy import that does not resolve to a displayable object is returned
    as is::

        sage: module.nothing = None
        sage: proxy = LazyImport('sage_doctest_lazy_module', 'nothing')
        sage: _resolve_lazy_import(proxy) is proxy
        True
        sage: del sys.modules['sage_doctest_lazy_module']

        sage: from sage.features import PythonModule
        sage: feature = PythonModule('math')
        sage: feature.hide()
        sage: proxy = LazyImport('math', 'sqrt', feature=feature)
        sage: _resolve_lazy_import(proxy) is proxy
        True
        sage: feature.unhide()
    """
    seen = []
    while isinstance(obj, LazyImport) and not any(obj is s for s in seen):
        seen.append(obj)
        try:
            resolved = obj._get_object()
        except FeatureNotPresentError:
            break
        if resolved is None:
            break
        obj = resolved
    return obj


class SageDisplayFormatter(DisplayFormatter):

    def __init__(self, *args, **kwds):
        """
        This is where the Sage rich objects are translated to IPython.

        INPUT/OUTPUT: see the IPython documentation

        EXAMPLES:

        This is part of how Sage works with the IPython output
        system. It cannot be used in doctests::

            sage: from sage.repl.display.formatter import SageDisplayFormatter
            sage: fmt = SageDisplayFormatter()
            Traceback (most recent call last):
            ...
            RuntimeError: check failed: current backend is invalid
        """
        super().__init__(*args, **kwds)
        from sage.repl.rich_output.display_manager import get_display_manager
        self.dm = get_display_manager()
        from sage.repl.rich_output.backend_ipython import BackendIPython
        self.dm.check_backend_class(BackendIPython)

        pt_formatter = self.formatters[PLAIN_TEXT]
        pt_formatter.observe(self._ipython_float_precision_changed,
                             names=['float_precision'])

    def format(self, obj, include=None, exclude=None):
        r"""
        Use the Sage rich output instead of IPython.

        INPUT/OUTPUT: see the IPython documentation

        EXAMPLES::

            sage: [identity_matrix(i) for i in range(3,7)]
            [
                                             [1 0 0 0 0 0]
                                [1 0 0 0 0]  [0 1 0 0 0 0]
                     [1 0 0 0]  [0 1 0 0 0]  [0 0 1 0 0 0]
            [1 0 0]  [0 1 0 0]  [0 0 1 0 0]  [0 0 0 1 0 0]
            [0 1 0]  [0 0 1 0]  [0 0 0 1 0]  [0 0 0 0 1 0]
            [0 0 1], [0 0 0 1], [0 0 0 0 1], [0 0 0 0 0 1]
            ]
            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%display ascii_art')   # indirect doctest
            sage: shell.run_cell("i = var('i')")
            sage: shell.run_cell('sum(i*x^i, i, 0, 10)')
                10      9      8      7      6      5      4      3      2
            10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x
            sage: shell.run_cell('%display default')
            sage: shell.quit()

        TESTS::

            sage: import os
            sage: import importlib.resources
            sage: import warnings
            sage: warnings.filterwarnings('ignore', category=DeprecationWarning,
            ....:     message=r'path is deprecated\. Use files\(\) instead\.')
            sage: from sage.repl.rich_output.backend_ipython import BackendIPython
            sage: backend = BackendIPython()
            sage: shell = get_test_shell()
            sage: backend.install(shell=shell)
            sage: shell.run_cell('get_ipython().display_formatter')
            <sage.repl.display.formatter.SageDisplayFormatter object at 0x...>
            sage: shell.run_cell('from IPython.display import Image')
            sage: with importlib.resources.path(sage.repl.rich_output, 'example.png') as example_png:
            ....:     shell.run_cell('ipython_image = Image("{0}")'.format(example_png))
            sage: shell.run_cell('ipython_image')
            <IPython.core.display.Image object>
            sage: shell.run_cell('get_ipython().display_formatter.format(ipython_image)')
            ({'image/png': ...,
              'text/plain': '<IPython.core.display.Image object>'},
            {})

        Test that IPython images and widgets still work even in latex output mode::

            sage: shell.run_cell('%display latex')   # indirect doctest
            sage: shell.run_cell('set(get_ipython().display_formatter.format(ipython_image)[0].keys())'
            ....:                ' == set(["text/plain", "image/png"])')
            True

            sage: shell.run_cell('import ipywidgets')
            sage: shell.run_cell('slider = ipywidgets.IntSlider()')
            sage: shell.run_cell('get_ipython().display_formatter.format(slider)')
            ...IntSlider(value=0)..., {})

            sage: shell.run_cell('%display default')
            sage: shell.quit()

        Test that ``__repr__`` is only called once when generating text output::

            sage: class Repper():
            ....:    def __repr__(self):
            ....:        print('__repr__ called')
            ....:        return 'I am repper'
            sage: Repper()
            __repr__ called
            I am repper

        A lazily imported class is displayed like the class it stands for.  In
        particular, IPython does not call an ordinary instance method such as
        :meth:`_repr_svg_` without an instance (:issue:`41697`)::

            sage: shell = get_test_shell()
            sage: shell.run_cell('Tableau')
            <class 'sage.combinat.tableau.Tableau'>
            sage: shell.run_cell('DyckWord')
            <class 'sage.combinat.dyck_word.DyckWord'>
            sage: shell.run_cell('PlanePartition')
            <class 'sage.combinat.plane_partition.PlanePartition'>

        The class itself is unaffected and still hands out its methods::

            sage: shell.run_cell('Tableau._repr_svg_')
            <function Tableau._repr_svg_ at ...>
            sage: shell.quit()
        """
        obj = _resolve_lazy_import(obj)
        sage_format, sage_metadata = self.dm.displayhook(obj)
        assert PLAIN_TEXT in sage_format, 'plain text is always present'

        # use Sage rich output for any except those native to IPython, but only
        # if it is not plain and dull
        if (not isinstance(obj, IPYTHON_NATIVE_TYPES) and
            not set(sage_format.keys()).issubset([PLAIN_TEXT]) and
                not isinstance(obj, Figure)):
            return sage_format, sage_metadata

        if self.ipython_display_formatter(obj):
            # object handled itself, don't proceed
            return {}, {}

        # try IPython display formatter
        if exclude is not None:
            exclude = list(exclude) + [PLAIN_TEXT]
        else:
            exclude = [PLAIN_TEXT]
        ipy_format, ipy_metadata = super().format(obj, include=include, exclude=exclude)
        if not ipy_format:
            return sage_format, sage_metadata
        ipy_format[PLAIN_TEXT] = sage_format[PLAIN_TEXT]
        if PLAIN_TEXT in sage_metadata:
            ipy_metadata[PLAIN_TEXT] = sage_metadata[PLAIN_TEXT]
        return ipy_format, ipy_metadata

    @staticmethod
    def _ipython_float_precision_changed(change):
        """
        Update the current float precision for the display of matrices in Sage.

        This function is called when the IPython ``%precision`` magic is
        invoked.

        TESTS::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%precision 4')
            '%.4f'
            sage: shell.run_cell('matrix.options.precision')  # indirect doctest
            4
            sage: shell.run_cell('%precision')
            '%r'
            sage: shell.run_cell('matrix.options.precision')  # indirect doctest
            None
        """
        try:
            from sage.matrix.constructor import options
        except ImportError:
            pass
        else:
            s = change.new
            if not s:
                # unset the precision
                options.precision = None
            else:
                try:
                    prec = int(s)
                    if prec >= 0:
                        options.precision = prec
                    # otherwise ignore the change
                except ValueError:
                    pass


class SagePlainTextFormatter(PlainTextFormatter):

    def __init__(self, *args, **kwds):
        r"""
        Improved plain text IPython formatter.

        In particular, it correctly print lists of matrices or other
        objects (see
        :meth:`sage.structure.parent.Parent._repr_option`).

        .. warning::

            This IPython formatter is NOT used. You could use it to
            enable Sage formatting in IPython, but Sage uses its own
            rich output system that is more flexible and supports
            different backends.

        INPUT/OUTPUT: see the IPython documentation

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.display_formatter.formatters['text/plain']
            <IPython.core.formatters.PlainTextFormatter object at 0x...>
            sage: shell.quit()
        """
        super().__init__(*args, **kwds)

    def __call__(self, obj):
        r"""
        Compute the pretty representation of the object.

        Adapted from ``IPython.core.formatters.PlainTextPrettyPrint``.

        INPUT:

        - ``obj`` -- anything

        OUTPUT: string; the plain text representation

        EXAMPLES::

            sage: from sage.repl.display.formatter import SagePlainTextFormatter
            sage: fmt = SagePlainTextFormatter()
            sage: fmt
            <sage.repl.display.formatter.SagePlainTextFormatter object at 0x...>
            sage: fmt(2)
            ---- calling ipython formatter ----
            '2'
            sage: a = identity_matrix(ZZ, 2)
            sage: fmt([a, a])
            ---- calling ipython formatter ----
            '[\n[1 0]  [1 0]\n[0 1], [0 1]\n]'
        """
        from sage.doctest import DOCTEST_MODE
        if DOCTEST_MODE:
            # Just to show that this is never executed in any other doctests in the Sage library
            print('---- calling ipython formatter ----')
        stream = StringIO()
        printer = SagePrettyPrinter(
            stream, self.max_width, self.newline)
        printer.pretty(obj)
        printer.flush()
        return stream.getvalue()
