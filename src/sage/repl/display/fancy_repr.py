# sage_setup: distribution = sagemath-repl
"""
Representations of objects
"""
# ****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import types
from io import StringIO

from IPython.lib.pretty import (
    _safe_getattr,
    _type_pprinters,
)

from sage.repl.display.util import format_list

_baseclass_reprs = (object.__repr__,)


class ObjectReprABC:
    """
    The abstract base class of an object representer.

    .. automethod:: __call__
    """

    def __repr__(self):
        """
        Return string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.display.fancy_repr import ObjectReprABC
            sage: ObjectReprABC()
            ObjectReprABC pretty printer
        """
        return '{0} pretty printer'.format(self.__class__.__name__)

    def __call__(self, obj, p, cycle):
        r"""
        Format object.

        INPUT:

        - ``obj`` -- anything; object to format

        - ``p`` -- PrettyPrinter instance

        - ``cycle`` -- boolean; whether there is a cycle

        OUTPUT:

        boolean; whether the representer is applicable to ``obj``. If
        ``True``, the string representation is appended to ``p``.

        EXAMPLES::

            sage: from sage.repl.display.fancy_repr import ObjectReprABC
            sage: ObjectReprABC().format_string(123)   # indirect doctest
            'Error: ObjectReprABC.__call__ is abstract'
        """
        p.text('Error: ObjectReprABC.__call__ is abstract')
        return True

    def format_string(self, obj):
        """
        For doctesting only: Directly return string.

        INPUT:

        - ``obj`` -- anything; object to format

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.display.fancy_repr import ObjectReprABC
            sage: ObjectReprABC().format_string(123)
            'Error: ObjectReprABC.__call__ is abstract'
        """
        from sage.repl.display.pretty_print import SagePrettyPrinter
        stream = StringIO()
        p = SagePrettyPrinter(stream, 79, '\n')
        ok = self(obj, p, False)
        if ok:
            p.flush()
            return stream.getvalue()
        else:
            return '--- object not handled by representer ---'


class SomeIPythonRepr(ObjectReprABC):

    def __init__(self):
        """
        Some selected representers from IPython.

        EXAMPLES::

            sage: from sage.repl.display.fancy_repr import SomeIPythonRepr
            sage: SomeIPythonRepr()
            SomeIPythonRepr pretty printer

        .. automethod:: __call__
        """
        type_repr = _type_pprinters.copy()
        del type_repr[type]
        del type_repr[types.BuiltinFunctionType]
        del type_repr[types.FunctionType]
        del type_repr[str]
        self._type_repr = type_repr

    def __call__(self, obj, p, cycle):
        """
        Format object.

        INPUT:

        - ``obj`` -- anything; object to format

        - ``p`` -- PrettyPrinter instance

        - ``cycle`` -- boolean; whether there is a cycle

        OUTPUT:

        boolean; whether the representer is applicable to ``obj``. If
        ``True``, the string representation is appended to ``p``.

        EXAMPLES::

            sage: from sage.repl.display.fancy_repr import SomeIPythonRepr
            sage: pp = SomeIPythonRepr()
            sage: pp.format_string(set([1, 2, 3]))
            '{1, 2, 3}'
        """
        try:
            pretty_repr = self._type_repr[type(obj)]
        except KeyError:
            return False
        pretty_repr(obj, p, cycle)
        return True


class LargeMatrixHelpRepr(ObjectReprABC):
    """
    Representation including help for large Sage matrices.

    .. automethod:: __call__
    """

    def __call__(self, obj, p, cycle):
        r"""
        Format matrix.

        INPUT:

        - ``obj`` -- anything; object to format

        - ``p`` -- PrettyPrinter instance

        - ``cycle`` -- boolean; whether there is a cycle

        OUTPUT:

        boolean; whether the representer is applicable to ``obj``. If
        ``True``, the string representation is appended to ``p``.

        EXAMPLES::

            sage: # needs sage.modules
            sage: from sage.repl.display.fancy_repr import LargeMatrixHelpRepr
            sage: M = identity_matrix(40)
            sage: pp = LargeMatrixHelpRepr()
            sage: pp.format_string(M)
            "40 x 40 dense matrix over Integer Ring (use the '.str()' method...)"
            sage: pp.format_string([M, M])
            '--- object not handled by representer ---'

        Leads to::

            sage: M                                                                     # needs sage.modules
            40 x 40 dense matrix over Integer Ring (use the '.str()' method...)
            sage: [M, M]                                                                # needs sage.modules
            [40 x 40 dense matrix over Integer Ring,
             40 x 40 dense matrix over Integer Ring]
        """
        if not p.toplevel():
            # Do not print the help for matrices inside containers
            return False
        from sage.structure.element import Matrix
        if not isinstance(obj, Matrix):
            return False
        from sage.matrix.constructor import options
        if obj.nrows() <= options.max_rows() and obj.ncols() <= options.max_cols():
            return False
        p.text(
            repr(obj) + " (use the '.str()' method to see the entries)"
        )
        return True


class PlainPythonRepr(ObjectReprABC):
    """
    The ordinary Python representation.

    .. automethod:: __call__
    """

    def __call__(self, obj, p, cycle):
        r"""
        Format matrix.

        INPUT:

        - ``obj`` -- anything; object to format

        - ``p`` -- PrettyPrinter instance

        - ``cycle`` -- boolean; whether there is a cycle

        OUTPUT:

        boolean; whether the representer is applicable to ``obj``. If
        ``True``, the string representation is appended to ``p``.

        EXAMPLES::

            sage: from sage.repl.display.fancy_repr import PlainPythonRepr
            sage: pp = PlainPythonRepr()
            sage: pp.format_string(type(1))
            "<class 'sage.rings.integer.Integer'>"

        Do not swallow a trailing newline at the end of the output of
        a custom representer. Note that it is undesirable to have a
        trailing newline, and if we don't display it you can't fix
        it::

            sage: class Newline():
            ....:     def __repr__(self):
            ....:         return 'newline\n'
            sage: n = Newline()
            sage: pp.format_string(n)
            'newline\n'
            sage: pp.format_string([n, n, n])
            '[newline\n, newline\n, newline\n]'
            sage: [n, n, n]
            [newline
             , newline
             , newline
             ]
        """
        klass = _safe_getattr(obj, '__class__', None) or type(obj)
        klass_repr = _safe_getattr(klass, '__repr__', None)
        if klass_repr in _baseclass_reprs:
            p.text(klass_repr(obj))
        else:
            # A user-provided repr. Find newlines and replace them with p.break_()
            try:
                output = repr(obj)
            except Exception:
                import sys
                import traceback
                objrepr = object.__repr__(obj).replace("object at", "at")
                exc = traceback.format_exception_only(sys.exc_info()[0], sys.exc_info()[1])
                exc = (''.join(exc)).strip()
                output = "<repr({}) failed: {}>".format(objrepr, exc)
            for idx, output_line in enumerate(output.split('\n')):
                if idx:
                    p.break_()
                p.text(output_line)
        return True


class TallListRepr(ObjectReprABC):
    """
    Special representation for lists with tall entries (e.g. matrices).

    .. automethod:: __call__
    """

    def __call__(self, obj, p, cycle):
        r"""
        Format list/tuple.

        INPUT:

        - ``obj`` -- anything; object to format

        - ``p`` -- PrettyPrinter instance

        - ``cycle`` -- boolean; whether there is a cycle

        OUTPUT:

        boolean; whether the representer is applicable to ``obj``. If
        ``True``, the string representation is appended to ``p``.

        EXAMPLES::

            sage: from sage.repl.display.fancy_repr import TallListRepr
            sage: format_list = TallListRepr().format_string
            sage: format_list([1, 2, identity_matrix(2)])                               # needs sage.modules
            '[\n      [1 0]\n1, 2, [0 1]\n]'

        Check that :issue:`18743` is fixed::

            sage: class Foo():
            ....:     def __repr__(self):
            ....:         return '''BBB    AA   RRR
            ....: B  B  A  A  R  R
            ....: BBB   AAAA  RRR
            ....: B  B  A  A  R  R
            ....: BBB   A  A  R   R'''
            ....:     def _repr_option(self, key):
            ....:         return key == 'ascii_art'
            sage: F = Foo()
            sage: [F, F]
            [
            BBB    AA   RRR    BBB    AA   RRR
            B  B  A  A  R  R   B  B  A  A  R  R
            BBB   AAAA  RRR    BBB   AAAA  RRR
            B  B  A  A  R  R   B  B  A  A  R  R
            BBB   A  A  R   R, BBB   A  A  R   R
            ]
        """
        if not (isinstance(obj, (tuple, list)) and len(obj) > 0):
            return False
        ascii_art_repr = False
        for o in obj:
            try:
                ascii_art_repr = ascii_art_repr or o._repr_option('ascii_art')
            except (AttributeError, TypeError):
                pass
            try:
                ascii_art_repr = ascii_art_repr or o.parent()._repr_option('element_ascii_art')
            except (AttributeError, TypeError):
                pass
        if not ascii_art_repr:
            return False
        output = format_list.try_format(obj)
        if output is None:
            return False
        p.text(output)
        return True
