# sage_setup: distribution = sagemath-repl
"""
The Sage pretty printer

Any transformation to a string goes through here. In other words, the
:class:`~sage.repl.displayhook.formatter.SagePlainTextFormatter` is
entirely implemented via :class:`SagePrettyPrinter`. Other formatters
may or may not use :class:`SagePrettyPrinter` to generate text output.

AUTHORS:

- Bill Cauchois (2009): initial version
- Jean-Baptiste Priez <jbp@kerios.fr> (2013): ASCII art
- Volker Braun (2013): refactored into DisplayHookBase
"""

# ****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from IPython.lib.pretty import PrettyPrinter

from sage.repl.display.fancy_repr import (TallListRepr, PlainPythonRepr,
                                          LargeMatrixHelpRepr,
                                          SomeIPythonRepr)


class SagePrettyPrinter(PrettyPrinter):

    DEBUG = False

    # These object representers will be tried, in this order, until
    # one is found that is able to deal with the object.
    pretty_repr = (
        TallListRepr(),
        LargeMatrixHelpRepr(),
        SomeIPythonRepr(),
        PlainPythonRepr(),
    )

    def toplevel(self):
        r"""
        Return whether we are currently at the top level.

        OUTPUT:

        boolean; whether we are currently pretty-printing an object at
        the outermost level (``True``), or whether the object is
        inside a container (``False``).

        EXAMPLES::

            sage: from sage.repl.display.pretty_print import SagePrettyPrinter
            sage: from io import StringIO
            sage: stream = StringIO()
            sage: spp = SagePrettyPrinter(stream, 78, '\n')
            sage: spp.toplevel()
            True
        """
        return len(self.stack) <= 1   # only the object currently being represented

    def __init__(self, output, max_width, newline, max_seq_length=None):
        """
        Pretty print Sage objects for the commandline.

        INPUT:

        See IPython documentation.

        EXAMPLES::

            sage: 123
            123

        IPython pretty printers::

            sage: set({1, 2, 3})
            {1, 2, 3}
            sage: dict(zzz=123, aaa=99, xab=10)    # sorted by keys
            {'aaa': 99, 'xab': 10, 'zzz': 123}

        These are overridden in IPython in a way that we feel is somewhat
        confusing, and we prefer to print them like plain Python which is
        more informative. See :issue:`14466` ::

            sage: 'this is a string'
            'this is a string'
            sage: type(123)
            <class 'sage.rings.integer.Integer'>
            sage: type
            <... 'type'>
            sage: import types
            sage: type('name', (), {})
            <class '__main__.name'>
            sage: types.BuiltinFunctionType
            <class 'builtin_function_or_method'>

            sage: def foo(): pass
            sage: foo
            <function foo at 0x...>
        """
        super().__init__(output, max_width, newline,
                         max_seq_length=max_seq_length)
        self.stack = []

    def pretty(self, obj):
        r"""
        Pretty print ``obj``.

        This is the only method that outside code should invoke.

        INPUT:

        - ``obj`` -- anything

        OUTPUT: string representation for object

        EXAMPLES::

            sage: from sage.repl.display.pretty_print import SagePrettyPrinter
            sage: from io import StringIO
            sage: stream = StringIO()
            sage: SagePrettyPrinter(stream, 78, '\n').pretty([type, 123, 'foo'])
            sage: stream.getvalue()
            "[<... 'type'>,"
        """
        obj_id = id(obj)
        cycle = obj_id in self.stack
        self.stack.append(obj_id)
        self.begin_group()
        try:
            ok = False
            for representation in self.pretty_repr:
                if self.DEBUG:
                    print('Trying {0}'.format(representation))
                ok = representation(obj, self, cycle)
                if self.DEBUG:
                    print('ok = {0}'.format(ok))
                if ok not in [True, False]:
                    raise RuntimeError('printer failed to return boolean')
                if ok:
                    break
            if not ok:
                raise RuntimeError('no printer registered for object')
        finally:
            self.end_group()
            self.stack.pop()
