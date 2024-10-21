# sage_setup: distribution = sagemath-categories
"""
Callable dictionaries
"""
# ****************************************************************************
#       Copyright (C) 2015 Nicolas M. Thiéry <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


cdef class CallableDict(dict):
    r"""
    Callable dictionary.

    This is a trivial subclass of :class:`dict` with an alternative
    view as a function.

    Typical use cases involve passing a dictionary `d` down to some
    tool that takes a function as input. The usual idiom in such use
    cases is to pass the ``d.__getitem__`` bound method. A pitfall is
    that this object is not picklable. When this feature is desired, a
    :class:`CallableDict` can be used instead. Note however that, with
    the current implementation, :class:`CallableDict` is slightly
    slower than ``d.__getitem__`` (see :issue:`6484` for benchmarks, and
    :issue:`18330` for potential for improvement).

    EXAMPLES::

        sage: from sage.misc.callable_dict import CallableDict
        sage: d = CallableDict({'one': 1, 'zwei': 2, 'trois': 3})
        sage: d['zwei']
        2
        sage: d('zwei')
        2

    In case the input is not in the dictionary, a :exc:`ValueError`
    is raised, for consistency with the function call syntax::

        sage: d[1]
        Traceback (most recent call last):
        ...
        KeyError: 1
        sage: d(1)
        Traceback (most recent call last):
        ...
        ValueError: 1 is not in dict
    """
    def __call__(self, key):
        r"""
        Return ``self[key]``.

        INPUT:

        - ``x`` -- any hashable object

        A :exc:`ValueError` is raised if ``x`` is not in ``self``.

        TESTS::

            sage: from sage.misc.callable_dict import CallableDict
            sage: d = CallableDict({'one': 1, 'zwei': 2, 'trois': 3})
            sage: d('one'), d('zwei'), d('trois')
            (1, 2, 3)
            sage: d('x')
            Traceback (most recent call last):
            ...
            ValueError: 'x' is not in dict
        """
        try:
            return self[key]
        except KeyError:
            raise ValueError(repr(key) + " is not in dict")

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.misc.callable_dict import CallableDict
            sage: d = CallableDict({1: 'a', 3: 'b', 2: 'c'}); d
            {1: 'a', 2: 'c', 3: 'b'}
        """
        from pprint import pformat
        return pformat(dict(self))
