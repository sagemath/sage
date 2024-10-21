# sage_setup: distribution = sagemath-gap
"""
Operations for LibGAP Elements

GAP functions for which several methods can be available are called
operations, so GAP ``Size`` is an example of an operation. This module
is for inspecting GAP operations from Python. In particular, it can
list the operations that take a particular LibGAP element as first
argument. This is used in tab completion, where Python ``x.[TAB]``
lists all GAP operations for which ``Operation(x, ...)`` is defined.
"""

import re
import string

from sage.structure.sage_object import SageObject
from sage.libs.gap.libgap import libgap

Length = libgap.function_factory('Length')
FlagsType = libgap.function_factory('FlagsType')
TypeObj = libgap.function_factory('TypeObj')
IS_SUBSET_FLAGS = libgap.function_factory('IS_SUBSET_FLAGS')
GET_OPER_FLAGS = libgap.function_factory('GET_OPER_FLAGS')
OPERATIONS = libgap.get_global('OPERATIONS')
NameFunction = libgap.function_factory('NameFunction')


NAME_RE = re.compile(r'(Setter|Getter|Tester)\((.*)\)')


class OperationInspector(SageObject):

    def __init__(self, libgap_element):
        """
        Information about operations that can act on a given LibGAP element.

        INPUT:

        - ``libgap_element`` -- libgap element

        EXAMPLES::

            sage: from sage.libs.gap.operations import OperationInspector
            sage: OperationInspector(libgap(123))
            Operations on 123
        """
        self._obj = libgap_element
        self.flags = FlagsType(TypeObj(self.obj))

    def _repr_(self):
        """
        Return the string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.libs.gap.operations import OperationInspector
            sage: opr = OperationInspector(libgap(123))
            sage: opr._repr_()
            'Operations on 123'
        """
        return 'Operations on {0}'.format(repr(self._obj))

    @property
    def obj(self):
        """
        The first argument for the operations.

        OUTPUT: a Libgap object

        EXAMPLES::

            sage: from sage.libs.gap.operations import OperationInspector
            sage: x = OperationInspector(libgap(123))
            sage: print(x.obj)
            123
        """
        return self._obj

    def operations(self):
        """
        Return the GAP operations for :meth:`obj`.

        OUTPUT: list of GAP operations

        EXAMPLES::

            sage: from sage.libs.gap.operations import OperationInspector
            sage: x = OperationInspector(libgap(123))
            sage: Unknown = libgap.function_factory('Unknown')
            sage: Unknown in x.operations()
            True
        """
        def mfi(o):
            filts = GET_OPER_FLAGS(o)
            return any(all(IS_SUBSET_FLAGS(self.flags, fl) for fl in fls)
                       for fls in filts)

        return (op for op in OPERATIONS if mfi(op))

    def op_names(self):
        """
        Return the names of the operations.

        OUTPUT: list of strings

        EXAMPLES::

            sage: from sage.libs.gap.operations import OperationInspector
            sage: x = OperationInspector(libgap(123))
            sage: 'Sqrt' in x.op_names()
            True
        """
        result = set()
        for f in self.operations():
            name = NameFunction(f).sage()
            if name[0] not in string.ascii_letters:
                continue
            match = NAME_RE.match(name)
            if match:
                result.add(match.groups()[1])
            else:
                result.add(name)
        return sorted(result)
