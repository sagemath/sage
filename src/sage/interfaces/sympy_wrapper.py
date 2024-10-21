# sage_setup: distribution = sagemath-symbolics
# sage.doctest: optional - sympy
"""
Wrapper Class for Sage Sets as SymPy Sets
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sympy.core.basic import Basic
from sympy.core.decorators import sympify_method_args
from sympy.core.sympify import sympify
from sympy.sets.sets import Set


@sympify_method_args
class SageSet(Set):
    r"""
    Wrapper for a Sage set providing the SymPy Set API.

    Parents in the category :class:`sage.categories.sets_cat.Sets`, unless
    a more specific method is implemented, convert to SymPy by creating
    an instance of this class.

    EXAMPLES::

        sage: F = Family([2, 3, 5, 7]); F
        Family (2, 3, 5, 7)
        sage: sF = F._sympy_(); sF            # indirect doctest
        SageSet(Family (2, 3, 5, 7))
        sage: sF._sage_() is F
        True
        sage: bool(sF)
        True
        sage: len(sF)
        4
        sage: list(sF)
        [2, 3, 5, 7]
        sage: sF.is_finite_set
        True
    """

    def __new__(cls, sage_set):
        r"""
        Construct a wrapper for a Sage set.

        TESTS::

            sage: from sage.interfaces.sympy_wrapper import SageSet
            sage: F = Set([1, 2]); F
            {1, 2}
            sage: sF = SageSet(F); sF
            SageSet({1, 2})
        """
        return Basic.__new__(cls, sage_set)

    def _sage_(self):
        r"""
        Return the underlying Sage set of the wrapper ``self``.

        EXAMPLES::

            sage: F = Family([1, 2])
            sage: F is Family([1, 2])
            False
            sage: sF = F._sympy_(); sF
            SageSet(Family (1, 2))
            sage: sF._sage_() is F
            True
        """
        return self._args[0]

    @property
    def is_empty(self):
        r"""
        Return whether the set ``self`` is empty.

        EXAMPLES::

            sage: Empty = Family([])
            sage: sEmpty = Empty._sympy_()
            sage: sEmpty.is_empty
            True
        """
        return self._sage_().is_empty()

    @property
    def is_finite_set(self):
        r"""
        Return whether the set ``self`` is finite.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: W = WeylGroup(["A",1,1])
            sage: sW = W._sympy_(); sW
            SageSet(Weyl Group of type ['A', 1, 1] (as a matrix group acting on the root space))
            sage: sW.is_finite_set
            False
        """
        return self._sage_().is_finite()

    @property
    def is_iterable(self):
        r"""
        Return whether the set ``self`` is iterable.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: W = WeylGroup(["A",1,1])
            sage: sW = W._sympy_(); sW
            SageSet(Weyl Group of type ['A', 1, 1] (as a matrix group acting on the root space))
            sage: sW.is_iterable
            True
        """
        from sage.categories.enumerated_sets import EnumeratedSets
        return self._sage_() in EnumeratedSets()

    def __iter__(self):
        r"""
        Iterator for the set ``self``.

        EXAMPLES::

            sage: sPrimes = Primes()._sympy_(); sPrimes
            SageSet(Set of all prime numbers: 2, 3, 5, 7, ...)
            sage: iter_sPrimes = iter(sPrimes)
            sage: next(iter_sPrimes), next(iter_sPrimes), next(iter_sPrimes)
            (2, 3, 5)
        """
        for element in self._sage_():
            yield sympify(element)

    def _contains(self, element):
        """
        Return whether ``element`` is an element of the set ``self``.

        EXAMPLES::

            sage: sPrimes = Primes()._sympy_(); sPrimes
            SageSet(Set of all prime numbers: 2, 3, 5, 7, ...)
            sage: 91 in sPrimes
            False

            sage: from sympy.abc import p
            sage: sPrimes.contains(p)
            Contains(p, SageSet(Set of all prime numbers: 2, 3, 5, 7, ...))

            sage: p in sPrimes
            Traceback (most recent call last):
            ...
            TypeError: did not evaluate to a bool: None
        """
        if element.is_symbol:
            # keep symbolic
            return None
        return element in self._sage_()

    def __len__(self):
        """
        Return the cardinality of the finite set ``self``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: sB3 = WeylGroup(["B", 3])._sympy_(); sB3
            SageSet(Weyl Group of type ['B', 3] (as a matrix group acting on the ambient space))
            sage: len(sB3)
            48
        """
        return len(self._sage_())

    def __str__(self):
        """
        Return the print representation of ``self``.

        EXAMPLES::

            sage: sPrimes = Primes()._sympy_()
            sage: str(sPrimes)  # indirect doctest
            'SageSet(Set of all prime numbers: 2, 3, 5, 7, ...)'
            sage: repr(sPrimes)
            'SageSet(Set of all prime numbers: 2, 3, 5, 7, ...)'
        """
        # Provide this method so that sympy's printing code does not try to inspect
        # the Sage object.
        return f"SageSet({self._sage_()})"

    __repr__ = __str__
