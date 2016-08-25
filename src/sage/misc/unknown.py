"""
The Unknown truth value

AUTHORS:

- Florent Hivert (2010): initial version.
"""
from __future__ import print_function

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation


class UnknownClass(UniqueRepresentation, SageObject):
    """
    TESTS::

        sage: TestSuite(Unknown).run()
    """
    def __init__(self):
        """
        The ``Unknown`` truth value

        EXAMPLES::

            sage: l = [False, Unknown, True]
            sage: for a in l: print([a and b for b in l])
            [False, False, False]
            [Unknown, Unknown, Unknown]
            [False, Unknown, True]

            sage: for a in l: print([a or b  for b in l])
            [False, Unknown, True]
            [False, Unknown, True]
            [True, True, True]

        ..warning:: Unless PEP 335 is accepted, in the following cases,
        ``and``, ``not`` and ``or`` return a somewhat wrong value::

            sage: not Unknown         # should return Unknown
            True
            sage: Unknown and False   # should return False
            Unknown
            sage: Unknown or False    # should return Unknown
            False
        """

    def _repr_(self):
        """
        TESTS::

            sage: Unknown
            Unknown
        """
        return "Unknown"

    def __nonzero__(self):
        """
        When evaluated in a boolean context ``Unknown()`` is evalutated into
        ``False``.

        EXAMPLES::

            sage: bool(Unknown)
            False
            sage: not Unknown
            True
        """
        return False

    def __and__(self, other):
        """
        The ``and`` logical connector.

        ..warning:: This is not used by ``and`` unless PEP 335 is accepted.

        EXAMPLES::

            sage: Unknown & False
            False
            sage: Unknown & Unknown
            Unknown
            sage: Unknown & True
            Unknown

        Compare with::

            sage: Unknown and False    # should return False
            Unknown
        """
        if other is False: return False
        else:              return self

    def __or__(self, other):
        """
        The ``or`` logical connector.

        ..warning:: This is not used by ``or`` unless PEP 335 is accepted.

        EXAMPLES::

            sage: Unknown | False
            Unknown
            sage: Unknown | Unknown
            Unknown
            sage: Unknown | True
            True

        Compare with::

            sage: Unknown or False    # should return Unknown
            False
        """
        if other: return True
        else:     return self

    def __not__(self):
        """
        The ``not`` logical connector.

        ..warning:: This is not used by ``not`` unless PEP 335 is accepted.

        EXAMPLES::

            sage: Unknown.__not__()
            Unknown

        Compare with::

            sage: not Unknown  # should return Unknown
            True
        """
        return self

    def __cmp__(self, other):
        """
        Comparison of truth value.

        EXAMPLES::

            sage: l = [False, Unknown, True]
            sage: for a in l: print([a < b for b in l])
            [False, True, True]
            [False, False, True]
            [False, False, False]

            sage: for a in l: print([a <= b for b in l])
            [True, True, True]
            [False, True, True]
            [False, False, True]
        """
        if other is self:
            return 0
        if isinstance(other, bool):
            if other:
                return -1
            else:
                return +1
        else:
            raise ValueError("Unable to compare {} with {}".format(self, other))

Unknown = UnknownClass()
