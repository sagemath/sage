# sage_setup: distribution = sagemath-objects
"""
The Unknown truth value

The ``Unknown`` object is used in Sage in several places as return value
in addition to ``True`` and ``False``, in order to signal uncertainty
about or inability to compute the result. ``Unknown`` can be identified
using ``is``, or by catching :exc:`UnknownError` from a boolean operation.

.. WARNING::

    Calling ``bool()`` with ``Unknown`` as argument will throw an
    :exc:`UnknownError`. This also means that in the following cases,
    ``and``, ``not``, and ``or`` fail or return a somewhat wrong value::

        sage: not Unknown         # should return Unknown
        Traceback (most recent call last):
        ...
        UnknownError: Unknown does not evaluate in boolean context
        sage: Unknown and False   # should return False
        Traceback (most recent call last):
        ...
        UnknownError: Unknown does not evaluate in boolean context
        sage: Unknown or False    # should return Unknown
        Traceback (most recent call last):
        ...
        UnknownError: Unknown does not evaluate in boolean context

EXAMPLES::

    sage: def func(n):
    ....:     if n > 0:
    ....:         return True
    ....:     elif n < 0:
    ....:         return False
    ....:     else:
    ....:         return Unknown

Using direct identification::

    sage: for n in [-3, 0, 12]:
    ....:    res = func(n)
    ....:    if res is True:
    ....:        print("n={} is positive".format(n))
    ....:    elif res is False:
    ....:        print("n={} is negative".format(n))
    ....:    else:
    ....:        print("n={} is neither positive nor negative".format(n))
    n=-3 is negative
    n=0 is neither positive nor negative
    n=12 is positive

Using :exc:`UnknownError`::

    sage: for n in [-3, 0, 12]:
    ....:    try:
    ....:        if func(n):
    ....:            print("n={} is positive".format(n))
    ....:        else:
    ....:            print("n={} is negative".format(n))
    ....:    except UnknownError:
    ....:        print("n={} is neither positive nor negative".format(n))
    n=-3 is negative
    n=0 is neither positive nor negative
    n=12 is positive

AUTHORS:

- Florent Hivert (2010): initial version.
- Ralf Stephan, Vincent Delecroix (2018-2020): redesign
"""
# ****************************************************************************
#       Copyright (C) 2010 Florent Hivert <florent.hivert at lri.fr>
#                     2018-2020 Ralf Stefan
#                     2018-2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import richcmp_method, rich_to_bool


class UnknownError(TypeError):
    """
    Raised whenever :class:`Unknown` is used in a boolean operation.

    EXAMPLES::

        sage: not Unknown
        Traceback (most recent call last):
        ...
        UnknownError: Unknown does not evaluate in boolean context
    """
    pass


@richcmp_method
class UnknownClass(UniqueRepresentation):
    """
    The Unknown truth value.

    The ``Unknown`` object is used in Sage in several places as return value
    in addition to ``True`` and ``False``, in order to signal uncertainty
    about or inability to compute the result. ``Unknown`` can be identified
    using ``is``, or by catching :exc:`UnknownError` from a boolean
    operation.

    .. WARNING::

        Calling ``bool()`` with ``Unknown`` as argument will throw an
        :exc:`UnknownError`. This also means that applying ``and``, ``not``,
        and ``or`` to ``Unknown`` might fail.

    TESTS::

        sage: TestSuite(Unknown).run()
    """
    def __repr__(self):
        """
        TESTS::

            sage: Unknown
            Unknown
        """
        return "Unknown"

    def __bool__(self):
        """
        When evaluated in a boolean context ``Unknown`` raises a :exc:`UnknownError`.

        EXAMPLES::

            sage: bool(Unknown)
            Traceback (most recent call last):
            ...
            UnknownError: Unknown does not evaluate in boolean context
            sage: not Unknown
            Traceback (most recent call last):
            ...
            UnknownError: Unknown does not evaluate in boolean context
        """
        raise UnknownError('Unknown does not evaluate in boolean context')

    def __and__(self, other):
        """
        The ``&`` logical operation.

        EXAMPLES::

            sage: Unknown & False
            False
            sage: Unknown & Unknown
            Unknown
            sage: Unknown & True
            Unknown

            sage: Unknown.__or__(3)
            NotImplemented
        """
        if other is False:
            return False
        elif other is True or other is Unknown:
            return self
        else:
            return NotImplemented

    def __or__(self, other):
        """
        The ``|`` logical connector.

        EXAMPLES::

            sage: Unknown | False
            Unknown
            sage: Unknown | Unknown
            Unknown
            sage: Unknown | True
            True

            sage: Unknown.__or__(3)
            NotImplemented
        """
        if other is True:
            return True
        elif other is False or other is Unknown:
            return self
        else:
            return NotImplemented

    def __richcmp__(self, other, op):
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
            return rich_to_bool(op, 0)
        if not isinstance(other, bool):
            return NotImplemented
        if other:
            return rich_to_bool(op, -1)
        else:
            return rich_to_bool(op, +1)


Unknown = UnknownClass()
