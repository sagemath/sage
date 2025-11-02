"""
The set of prime numbers and its subsets defined by congruence conditions

AUTHORS:

 - William Stein (2005): original version
 - Florent Hivert (2009-11): adapted to the category framework
 - Xavier Caruso (2025-10): implement congruence conditions
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#                     2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#                     2025 Xavier Caruso <xavier@caruso.ovh>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from .set import Set_generic
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.arith.misc import euler_phi
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import EmptySetError


class Primes(Set_generic, UniqueRepresentation):
    r"""
    The set of prime numbers and some of its subsets.

    EXAMPLES:

    The set of all primes numbers::

        sage: P = Primes(); P
        Set of all prime numbers: 2, 3, 5, 7, ...

    The arguments ``modulus`` and ``classes`` allows for constructing
    subsets given by congruence conditions::

        sage: Primes(modulus=4)
        Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...

    By default the congruence class `1` is selected, but we can specify any
    subset of congruence classes::

        sage: Primes(modulus=4, classes=[3])
        Set of prime numbers congruent to 3 modulo 4: 3, 7, 11, 19, ...
        sage: Primes(modulus=8, classes=[1, 3])
        Set of prime numbers congruent to 1, 3 modulo 8: 3, 11, 17, 19, ...

    If possible, the congruence conditions are simplified::

        sage: Primes(modulus=8, classes=[1, 5])
        Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...

    We can create a finite set of primes by passing in ``modulus=0``::

        sage: Primes(modulus=0, classes=[2, 3, 5, 11])
        Finite set of prime numbers: 2, 3, 5, 11

    We show various operations that can be performed on these sets::

        sage: P = Primes(modulus=4); P
        Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...
        sage: P.cardinality()
        +Infinity
        sage: P[:10]
        [5, 13, 17, 29, 37, 41, 53, 61, 73, 89]
        sage: P.next(500)
        509

    ::

        sage: Q = Primes(modulus=4, classes=[3])
        sage: PQ = P.union(Q)
        sage: PQ
        Set of all prime numbers with 2 excluded: 3, 5, 7, 11, ...
        sage: PQ.complement_in_primes()
        Finite set of prime numbers: 2
        sage: PQ.complement_in_primes().cardinality()
        1
    """
    @staticmethod
    def __classcall__(cls, modulus=1, classes=None, exceptions=None):
        """
        Normalize the input.

        INPUT:

        - ``modulus`` -- an integer (default: ``1``)

        - ``classes`` -- a list of integers (default: ``[1]``), the
          congruence classes (modulo ``modulus``) included in this
          set

        - ``exceptions`` -- a dictionary with items of the form
          ``c: v`` where ``c`` is an integer and ``v`` is a boolean;
          if ``v`` is ``True`` (resp. ``False``) then ``c`` is added
          to (resp. removed from) this set

        TESTS::

            sage: Primes(modulus=10)
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: Primes(modulus=10) == Primes(modulus=5)
            True

            sage: Primes(modulus=9, classes=[1, 3, 4, 7])
            Set of prime numbers congruent to 1 modulo 3 with 3 included: 3, 7, 13, 19, ...

            sage: Primes(modulus=5, exceptions={7: True, 11: False})
            Set of prime numbers congruent to 1 modulo 5 with 7 included and 11 excluded: 7, 31, 41, 61, ...
        """
        modulus = ZZ(modulus)
        if modulus < 0:
            modulus = -modulus
        if classes is None:
            classes = [ZZ(1)]
        if exceptions is None:
            exceptions = {}
        if isinstance(exceptions, (tuple, list)):
            exceptions = {c: v for c, v in exceptions}

        if modulus == 0:
            for c in classes:
                exceptions[c] = True
            modulus = ZZ(1)
            classes = []

        # We replace congruences of the form
        #   p = a (mod n) with gcd(a, n) > 1
        # (which only includes a finite number of primes)
        # with exceptions
        indic = modulus * [False]
        for c in classes:
            indic[ZZ(c) % modulus] = True
        for c in range(modulus):
            if modulus.gcd(c) > 1:
                if indic[c]:
                    if c == 0:
                        if modulus not in exceptions:
                            exceptions[modulus] = True
                    else:
                        if c not in exceptions:
                            exceptions[ZZ(c)] = True
                indic[c] = None

        # We normalize the congruence conditions
        # by minimizing the modulus
        for p, mult in modulus.factor():
            while mult > 0:
                m = modulus // p
                add_true = []
                add_false = []
                add_excluded = []
                for c in range(m):
                    cs = [indic[c + m*i] for i in range(p) if indic[c + m*i] is not None]
                    if not cs:
                        pass
                    elif all(cs):
                        if m.gcd(c) == 1:
                            add_true.append(c)
                        for i in range(p):
                            j = c + m*i
                            if indic[j] is None:
                                if j == 0:
                                    add_excluded.append(modulus)
                                else:
                                    add_excluded.append(j)
                    elif not any(cs):
                        if m.gcd(c) == 1:
                            add_false.append(c)
                    else:
                        mult = 0
                        break
                else:
                    for c in add_true:
                        indic[c] = True
                    for c in add_false:
                        indic[c] = False
                    for c in add_excluded:
                        if c not in exceptions:
                            exceptions[c] = False
                    modulus = m
                    mult -= 1

        # We format the final result
        classes = tuple([c for c in range(modulus) if indic[c] is True])
        exceptions_list = [(c, v) for c, v in exceptions.items()
                           if c.is_prime() and (v != (indic[c % modulus] is True))]
        exceptions_list.sort()

        return super().__classcall__(cls, modulus, classes, tuple(exceptions_list))

    def __init__(self, modulus, classes, exceptions):
        r"""
        Initialize this set.

        TESTS::

            sage: P = Primes(modulus=4)
            sage: P.category()
            Category of facade infinite enumerated sets
            sage: TestSuite(P).run()

        ::

            sage: Q = Primes(modulus=0, classes=[2, 3, 5])
            sage: Q.category()
            Category of facade finite enumerated sets
            sage: TestSuite(Q).run()
        """
        if classes:
            category = InfiniteEnumeratedSets()
        else:
            category = FiniteEnumeratedSets()
        super().__init__(facade=ZZ, category=category)
        self._modulus = modulus
        self._classes = set(classes)
        self._exceptions = dict(exceptions)
        if classes:
            self._elements = []
        else:
            self._elements = [c for c, v in exceptions if v]
            self._elements.sort()

    def _repr_(self):
        r"""
        Return a string representation of this subset.

        TESTS::

            sage: Primes(modulus=4).include(2).exclude(5)  # indirect doctest
            Set of prime numbers congruent to 1 modulo 4 with 2 included and 5 excluded: 2, 13, 17, 29, ...

            sage: E = Primes(modulus=0)
            sage: E  # indirect doctest
            Empty set of prime numbers

            sage: E.include(range(50), check=False)  # indirect doctest
            Finite set of prime numbers: 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
        """
        classes = sorted(list(self._classes))
        sc = ", ".join([str(c) for c in classes])
        included = []
        excluded = []
        for c, v in self._exceptions.items():
            if v:
                included.append(c)
            else:
                excluded.append(c)
        si = ", ".join([str(i) for i in sorted(included)])
        se = ", ".join([str(e) for e in sorted(excluded)])
        if not classes:
            if not included:
                return "Empty set of prime numbers"
            else:
                return "Finite set of prime numbers: %s" % si
        if self._modulus == 1:
            s = "Set of all prime numbers"
        else:
            s = "Set of prime numbers congruent to %s modulo %s" % (sc, self._modulus)
        if included:
            s += " with %s included" % si
        if excluded:
            if not included:
                s += " with %s excluded" % se
            else:
                s += " and %s excluded" % se
        s += ": %s, ..." % (", ".join([str(c) for c in self[:4]]))
        return s

    def __contains__(self, x):
        r"""
        Return ``True`` if `x` is in this set; ``False`` otherwise.

        INPUT:

        - ``x`` -- an integer

        EXAMPLES::

            sage: P = Primes(modulus=4)
            sage: 3 in P
            False
            sage: 9 in P
            False
            sage: 13 in P
            True

        TESTS::

            sage: x in P
            False
        """
        try:
            if x not in ZZ:
                return False
        except TypeError:
            return False
        x = ZZ(x)
        if not x.is_prime():
            return False
        e = self._exceptions.get(x, None)
        return (e is True) or (e is None and x % self._modulus in self._classes)

    def cardinality(self):
        r"""
        Return the cardinality of this set.

        EXAMPLES::

            sage: P = Primes(modulus=4); P
            Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...
            sage: P.cardinality()
            +Infinity

        ::

            sage: P = Primes(modulus=4, classes=[2]); P
            Finite set of prime numbers: 2
            sage: P.cardinality()
            1
        """
        if self.is_finite():
            return ZZ(len(self._elements))
        return infinity

    def first(self, n=None):
        r"""
        Return the first element in this set.

        EXAMPLES::

            sage: P = Primes(modulus=5); P
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: P.first()
            11
        """
        if self.is_empty():
            raise EmptySetError
        return self.next(1)

    def next(self, x):
        r"""
        Return the smallest element in this set strictly
        greater than ``x``.

        INPUT:

        - ``x`` -- an integer

        EXAMPLES::

            sage: P = Primes(modulus=5); P
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: P.next(1000)
            1021

        If there is no element greater than the given bound, an
        error is raised::

            sage: P = Primes(modulus=0, classes=[2, 5]); P
            Finite set of prime numbers: 2, 5
            sage: P.next(10)
            Traceback (most recent call last):
            ...
            ValueError: no element greater that 10 in this set
        """
        x = ZZ(x)
        if self._classes:
            while True:
                x = x.next_prime()
                e = self._exceptions.get(x, None)
                if (e is True) or (e is None and x % self._modulus in self._classes):
                    return x

        if not self._elements or x >= self._elements[-1]:
            raise ValueError("no element greater that %s in this set" % x)
        min = 0
        max = len(self._elements)
        while min < max:
            i = (min + max) // 2
            if self._elements[i] <= x:
                min = i + 1
            if self._elements[i] > x:
                max = i
        return self._elements[min]

    def _an_element_(self):
        r"""
        Return an element in this set.

        EXAMPLES::

            sage: P = Primes()
            sage: P.an_element()  # indirect doctest
            43

        If the set is empty, an error is raised::

            sage: P = Primes(modulus=12, classes=[4])
            sage: P.an_element()  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: this set is empty
        """
        if self.is_finite():
            if self._elements:
                return self._elements[0]
            raise ValueError("this set is empty")
        return self.next(42)

    def unrank(self, n):
        r"""
        Return the ``n``-th element of this set.

        INPUT:

        - ``n`` -- an integer

        EXAMPLES::

            sage: P = Primes(modulus=5); P
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: P[0]   # indirect doctest
            11
            sage: P[10]  # indirect doctest
            211

        If there is less than `n` elements in this set, an error
        is raised::

            sage: P = Primes(modulus=0, classes=[2, 5]); P
            Finite set of prime numbers: 2, 5
            sage: P[1]
            5
            sage: P[2]  # indirect doctest
            Traceback (most recent call last):
            ...
            IndexError: this set has not enough elements

        TESTS::

            sage: P.unrank(-1)
            Traceback (most recent call last):
            ...
            IndexError: index must be nonnegative
        """
        if n < 0:
            raise IndexError("index must be nonnegative")
        if self.is_finite():
            if len(self._elements) <= n:
                raise IndexError("this set has not enough elements")
            return self._elements[n]
        if self._elements:
            x = self._elements[-1]
        else:
            x = 1
        while len(self._elements) <= n:
            x = self.next(x)
            self._elements.append(x)
        return self._elements[n]

    def is_empty(self):
        r"""
        Return ``True`` if this set is empty; ``False`` otherwise.

        EXAMPLES::

            sage: P = Primes(modulus=6); P
            Set of prime numbers congruent to 1 modulo 3: 7, 13, 19, 31, ...
            sage: P.is_empty()
            False

        ::

            sage: P = Primes(modulus=6, classes=[6]); P
            Empty set of prime numbers
            sage: P.is_empty()
            True

        .. SEEALSO::

            :meth:`is_finite`, :meth:`is_cofinite`
        """
        return not bool(self._classes) and not bool(self._exceptions)

    def is_finite(self):
        r"""
        Return ``True`` if this set is finite; ``False`` otherwise.

        EXAMPLES::

            sage: P = Primes(modulus=5); P
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: P.is_finite()
            False

        ::

            sage: P = Primes(modulus=0, classes=[2, 5]); P
            Finite set of prime numbers: 2, 5
            sage: P.is_finite()
            True

        .. SEEALSO::

            :meth:`is_empty`, :meth:`is_cofinite`
        """
        return not bool(self._classes)

    def is_cofinite(self):
        r"""
        Return ``True`` if this set is cofinite in the set
        of all prime numbers; ``False`` otherwise.

        EXAMPLES::

            sage: P = Primes(modulus=4); P
            Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...
            sage: P.is_cofinite()
            False

        ::

            sage: P = Primes(modulus=4, classes=[1, 3]); P
            Set of all prime numbers with 2 excluded: 3, 5, 7, 11, ...
            sage: P.is_cofinite()
            True

        .. SEEALSO::

            :meth:`is_empty`, :meth:`is_finite`
        """
        return self._modulus == 1 and bool(self._classes)

    def density(self):
        r"""
        Return the density of this set in the set of all
        prime numbers.

        EXAMPLES::

            sage: P = Primes(modulus=4); P
            Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...
            sage: P.density()
            1/2
        """
        return len(self._classes) / euler_phi(self._modulus)

    def include(self, elements, check=True):
        r"""
        Return this set with the integers in ``elements`` included.

        INPUT:

        - ``elements`` -- an integer, or a tuple/list of integers

        - ``check`` -- a boolean (default: ``True``); if ``False``,
          do not raise an error if we try to add composite numbers

        EXAMPLES::

            sage: P = Primes(modulus=5); P
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: P.include(2)
            Set of prime numbers congruent to 1 modulo 5 with 2 included: 2, 11, 31, 41, ...
            sage: P.include([2, 3])
            Set of prime numbers congruent to 1 modulo 5 with 2, 3 included: 2, 3, 11, 31, ...

        If we try to include an element which is already in the set,
        nothing changes::

            sage: P.include(11)
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...

        Trying to include a composite number results in an error::

            sage: P.include(10)
            Traceback (most recent call last):
            ...
            ValueError: 10 is not a prime number

        We can avoid this by passing in ``check=False``; in this case,
        composite numbers are however not added to the set.
        This behavior can be convenient if one wants to add all prime
        numbers in a range::

            sage: P.include(range(20, 30), check=False)
            Set of prime numbers congruent to 1 modulo 5 with 23, 29 included: 11, 23, 29, 31, ...

        .. SEEALSO::

            :meth:`exclude`
        """
        if elements in ZZ:
            elements = [elements]
        exceptions = self._exceptions.copy()
        for x in elements:
            x = ZZ(x)
            if check and not x.is_prime():
                raise ValueError("%s is not a prime number" % x)
            exceptions[x] = True
        return Primes(self._modulus, self._classes, exceptions)

    def exclude(self, elements):
        r"""
        Return this set with the integers in ``elements`` excluded.

        INPUT:

        - ``elements`` -- an integer, or a tuple/list of integers

        EXAMPLES::

            sage: P = Primes(modulus=5); P
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: P.exclude(11)
            Set of prime numbers congruent to 1 modulo 5 with 11 excluded: 31, 41, 61, 71, ...
            sage: P.exclude([11, 31])
            Set of prime numbers congruent to 1 modulo 5 with 11, 31 excluded: 41, 61, 71, 101, ...

        If we try to exclude an element which is not in the set,
        nothing changes::

            sage: P.exclude(2)
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...

        .. SEEALSO::

            :meth:`include`
        """
        if elements in ZZ:
            elements = [elements]
        exceptions = self._exceptions.copy()
        for x in elements:
            exceptions[x] = False
        return Primes(self._modulus, self._classes, exceptions)

    def complement_in_primes(self):
        r"""
        Return the complement of this set in the set of all prime
        numbers.

        EXAMPLES::

            sage: P = Primes(modulus=4); P
            Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...
            sage: Q = P.complement_in_primes(); Q
            Set of prime numbers congruent to 3 modulo 4 with 2 included: 2, 3, 7, 11, ...

        We check that the union of `P` and `Q` is the whole set of
        prime numbers::

            sage: P.union(Q)
            Set of all prime numbers: 2, 3, 5, 7, ...

        and that the intersection is empty::

            sage: P.intersection(Q)
            Empty set of prime numbers

        .. SEEALSO::

            :meth:`intersection`, :meth:`union`
        """
        modulus = self._modulus
        classes = [c for c in range(modulus)
                   if c % self._modulus not in self._classes]
        exceptions = {c: not v for c, v in self._exceptions.items()}
        return Primes(modulus, classes, exceptions)

    def intersection(self, other):
        r"""
        Return the intesection of this set with ``other``.

        INPUT:

        - ``other`` -- a subset of the set of prime numbers

        EXAMPLES::

            sage: P = Primes(modulus=5); P
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: Q = Primes(modulus=3, classes=[2]); Q
            Set of prime numbers congruent to 2 modulo 3: 2, 5, 11, 17, ...
            sage: P.intersection(Q)
            Set of prime numbers congruent to 11 modulo 15: 11, 41, 71, 101, ...

        TESTS::

            sage: P.intersection(ZZ) == P
            True
            sage: P.intersection(RR)
            Traceback (most recent call last):
            ...
            NotImplementedError: boolean operations are only implemented with other sets of prime numbers

        .. SEEALSO::

            :meth:`complement_in_primes`, :meth:`union`
        """
        if other is ZZ:
            return self
        if not isinstance(other, Primes):
            raise NotImplementedError("boolean operations are only implemented with other sets of prime numbers")
        modulus = self._modulus.lcm(other._modulus)
        classes = [c for c in range(modulus)
                   if (c % self._modulus in self._classes
                   and c % other._modulus in other._classes)]
        exceptions = {c: v for c, v in self._exceptions.items()
                      if not v or c in other}
        exceptions.update((c, v) for c, v in other._exceptions.items()
                          if not v or c in self)
        return Primes(modulus, classes, exceptions)

    def union(self, other):
        r"""
        Return the intesection of this set with ``other``.

        INPUT:

        - ``other`` -- a subset of the set of prime numbers

        EXAMPLES::

            sage: P = Primes(modulus=5); P
            Set of prime numbers congruent to 1 modulo 5: 11, 31, 41, 61, ...
            sage: Q = Primes(modulus=3, classes=[2]); Q
            Set of prime numbers congruent to 2 modulo 3: 2, 5, 11, 17, ...
            sage: P.union(Q)
            Set of prime numbers congruent to 1, 2, 8, 11, 14 modulo 15 with 5 included: 2, 5, 11, 17, ...

        TESTS::
            sage: P = Primes(modulus=5, exceptions={5: True, 11: False}); P
            Set of prime numbers congruent to 1 modulo 5 with 5 included and 11 excluded: 5, 31, 41, 61, ...
            sage: Q = Primes(modulus=4, exceptions={13: False, 11: True}); Q
            Set of prime numbers congruent to 1 modulo 4 with 11 included and 13 excluded: 5, 11, 17, 29, ...
            sage: P.union(Q)
            Set of prime numbers congruent to 1, 9, 11, 13, 17 modulo 20 with 5 included and 13 excluded: 5, 11, 17, 29, ...

            sage: P.union(ZZ) == ZZ
            True
            sage: P.union(RR)
            Traceback (most recent call last):
            ...
            NotImplementedError: boolean operations are only implemented with other sets of prime numbers

        .. SEEALSO::

            :meth:`complement_in_primes`, :meth:`intersection`
        """
        if other is ZZ:
            return ZZ
        if not isinstance(other, Primes):
            raise NotImplementedError("boolean operations are only implemented with other sets of prime numbers")
        modulus = self._modulus.lcm(other._modulus)
        classes = [c for c in range(modulus)
                   if (c % self._modulus in self._classes
                    or c % other._modulus in other._classes)]
        exceptions = {}
        for c, v in self._exceptions.items():
            if v:
                exceptions[c] = True
            if not v and c not in other:
                exceptions[c] = False
        for c, v in other._exceptions.items():
            if v:
                exceptions[c] = True
            if not v and c not in self:
                exceptions[c] = False
        return Primes(modulus, classes, exceptions)

    def is_subset(self, other):
        r"""
        Return ``True`` is this set of is subset of ``other``;
        ``False`` otherwise.

        EXAMPLES::

            sage: P = Primes(modulus=4); P
            Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...
            sage: Q = Primes(modulus=8); Q
            Set of prime numbers congruent to 1 modulo 8: 17, 41, 73, 89, ...
            sage: P.is_subset(Q)
            False
            sage: Q.is_subset(P)
            True

        .. SEEALSO::

            :meth:`is_supset`, :meth:`is_disjoint`
        """
        return self.intersection(other) == self

    def is_supset(self, other):
        r"""
        Return ``True`` is this set of is supset of ``other``;
        ``False`` otherwise.

        EXAMPLES::

            sage: P = Primes(modulus=4); P
            Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...
            sage: Q = Primes(modulus=8); Q
            Set of prime numbers congruent to 1 modulo 8: 17, 41, 73, 89, ...
            sage: P.is_supset(Q)
            True
            sage: Q.is_supset(P)
            False

        .. SEEALSO::

            :meth:`is_subset`, :meth:`is_disjoint`
        """
        return self.intersection(other) == other

    def is_disjoint(self, other):
        r"""
        Return ``True`` is this set of is disjoint from ``other``;
        ``False`` otherwise.

        EXAMPLES::

            sage: P = Primes(modulus=4); P
            Set of prime numbers congruent to 1 modulo 4: 5, 13, 17, 29, ...
            sage: Q = Primes(modulus=4, classes=[3]); Q
            Set of prime numbers congruent to 3 modulo 4: 3, 7, 11, 19, ...
            sage: P.is_disjoint(Q)
            True

        ::

            sage: R = Primes(modulus=5, classes=[3]); R
            Set of prime numbers congruent to 3 modulo 5: 3, 13, 23, 43, ...
            sage: P.is_disjoint(R)
            False
            sage: Q.is_disjoint(R)
            False

        .. SEEALSO::

            :meth:`is_subset`, :meth:`is_disjoint`
        """
        return self.intersection(other).is_empty()
