"""
Basis matroids

In a matroid, a basis is an inclusionwise maximal independent set.
The common cardinality of all bases is the rank of the matroid.
Matroids are uniquely determined by their set of bases.

This module defines the class
:class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`, which
internally represents a matroid as a set of bases. It is a subclass of
:mod:`BasisExchangeMatroid <sage.matroids.basis_exchange_matroid>`, and as
such it inherits all method from that class and from the class
:mod:`Matroid <sage.matroids.matroid>`. Additionally, it provides the
following methods:

    - :meth:`is_distinguished() <sage.matroids.basis_matroid.BasisMatroid.is_distinguished>`
    - :meth:`relabel() <sage.matroids.basis_matroid.BasisMatroid.relabel>`

Construction
============

A ``BasisMatroid`` can be created from another matroid, from a list of bases,
or from a list of nonbases. For a full description of allowed inputs, see
:class:`below <sage.matroids.basis_matroid.BasisMatroid>`. It is recommended
to use the :func:`Matroid() <sage.matroids.constructor.Matroid>` function for
easy construction of a ``BasisMatroid``. For direct access to the
``BasisMatroid`` constructor, run::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`.

EXAMPLES::

    sage: from sage.matroids.advanced import *
    sage: M1 = BasisMatroid(groundset='abcd', bases=['ab', 'ac', 'ad', 'bc', 'bd', 'cd'])
    sage: M2 = Matroid(['ab', 'ac', 'ad', 'bc', 'bd', 'cd'])
    sage: M1 == M2
    True

Implementation
==============

The set of bases is compactly stored in a bitset which takes
`O(binomial(N, R))` bits of space, where `N` is the cardinality of the
groundset and `R` is the rank. ``BasisMatroid`` inherits the matroid oracle
from its parent class ``BasisExchangeMatroid``, by providing the elementary
functions for exploring the base exchange graph. In addition, ``BasisMatroid``
has methods for constructing minors, duals, single-element extensions, for
testing matroid isomorphism and minor inclusion.

AUTHORS:

- Rudi Pendavingh, Stefan van Zwam (2013-04-01): initial version

TESTS::

    sage: F = matroids.named_matroids.Fano()
    sage: M = Matroid(bases=F.bases())
    sage: TestSuite(M).run()

Methods
=======
"""
#*****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
include 'sage/data_structures/bitset.pxi'
from matroid cimport Matroid
from basis_exchange_matroid cimport BasisExchangeMatroid
from itertools import permutations
from sage.arith.all import binomial
from set_system cimport SetSystem
from itertools import combinations

# class of general matroids, represented by their list of bases

cdef class BasisMatroid(BasisExchangeMatroid):
    """
    Create general matroid, stored as a set of bases.

    INPUT:

    - ``M`` (optional) -- a matroid.
    - ``groundset`` (optional) -- any iterable set.
    - ``bases`` (optional) -- a set of subsets of ``groundset``.
    - ``nonbases`` (optional) -- a set of subsets of ``groundset``.
    - ``rank`` (optional) -- a natural number

    EXAMPLES:

    The empty matroid::

        sage: from sage.matroids.advanced import *
        sage: M = BasisMatroid()
        sage: M.groundset()
        frozenset()
        sage: M.full_rank()
        0

    Create a BasisMatroid instance out of any other matroid::

        sage: from sage.matroids.advanced import *
        sage: F = matroids.named_matroids.Fano()
        sage: M = BasisMatroid(F)
        sage: F.groundset() == M.groundset()
        True
        sage: len(set(F.bases()).difference(M.bases()))
        0

    It is possible to provide either bases or nonbases::

        sage: from sage.matroids.advanced import *
        sage: M1 = BasisMatroid(groundset='abc', bases=['ab', 'ac'] )
        sage: M2 = BasisMatroid(groundset='abc', nonbases=['bc'])
        sage: M1 == M2
        True

    Providing only groundset and rank creates a uniform matroid::

        sage: from sage.matroids.advanced import *
        sage: M1 = BasisMatroid(matroids.Uniform(2, 5))
        sage: M2 = BasisMatroid(groundset=range(5), rank=2)
        sage: M1 == M2
        True

    We do not check if the provided input forms an actual matroid::

        sage: from sage.matroids.advanced import *
        sage: M1 = BasisMatroid(groundset='abcd', bases=['ab', 'cd'])
        sage: M1.full_rank()
        2
        sage: M1.is_valid()
        False

    """
    def __init__(self, M=None, groundset=None, bases=None, nonbases=None, rank=None):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: F = matroids.named_matroids.Fano()
            sage: M = BasisMatroid(F)
            sage: F.groundset() == M.groundset()
            True
            sage: len(set(F.bases()).difference(M.bases()))
            0
        """
        cdef SetSystem NB
        cdef long i

        if isinstance(M, BasisMatroid):
            BasisExchangeMatroid.__init__(self, groundset=(<BasisMatroid>M)._E, rank=(<BasisMatroid>M)._matroid_rank)
            bitset_init(self._bb, binom[(<BasisMatroid>M)._groundset_size][(<BasisMatroid>M)._matroid_rank])
            bitset_copy(self._bb, (<BasisMatroid>M)._bb)
            bitset_init(self._b, (<BasisMatroid>M)._bitset_size)
            bitset_copy(self._b, (<BasisMatroid>M)._b)
            bitset_copy(self._current_basis, (<BasisMatroid>M)._current_basis)
            self._bcount = (<BasisMatroid>M)._bcount
            return

        if isinstance(M, BasisExchangeMatroid):
            BasisExchangeMatroid.__init__(self, groundset=(<BasisExchangeMatroid>M)._E, rank=(<BasisExchangeMatroid>M)._matroid_rank)
            binom_init(len(M), M.full_rank())
            bc = binom[(<BasisExchangeMatroid>M)._groundset_size][(<BasisExchangeMatroid>M)._matroid_rank]
            bitset_init(self._bb, bc)
            bitset_set_first_n(self._bb, bc)
            NB = M.nonbases()
            for i in xrange(len(NB)):
                bitset_discard(self._bb, set_to_index(NB._subsets[i]))

            bitset_init(self._b, (<BasisExchangeMatroid>M)._bitset_size)
            self.reset_current_basis()
            self._bcount = bc - len(NB)
            return

        if M is not None:
            rank = M.full_rank()
            nonbases = M.nonbases()
            groundset = sorted(M.groundset())

        if groundset is None:
            groundset = frozenset()
        if rank is None:
            if bases is not None:
                rank = len(min(bases))
            elif nonbases is not None:
                rank = len(min(nonbases))
            else:
                rank = 0

        BasisExchangeMatroid.__init__(self, groundset=groundset, rank=rank)

        size = len(groundset)
        binom_init(size, rank)
        bitset_init(self._bb, binom[size][rank])
        bitset_init(self._b, max(size, 1))
        bitset_clear(self._bb)

        if bases is not None:
            if len(bases) == 0:
                raise ValueError("set of bases must be nonempty.")
            self._bcount = 0
            for B in bases:
                b = frozenset(B)
                if len(b) != self._matroid_rank:
                    raise ValueError("basis has wrong cardinality.")
                if not b.issubset(self._groundset):
                    raise ValueError("basis is not a subset of the groundset")
                self.__pack(self._b, b)
                i = set_to_index(self._b)
                if not bitset_in(self._bb, i):
                    self._bcount += 1
                bitset_add(self._bb, i)
        else:
            bitset_complement(self._bb, self._bb)
            self._bcount = binom[size][rank]
            if nonbases is not None:
                for B in nonbases:
                    b = frozenset(B)
                    if len(b) != self._matroid_rank:
                        raise ValueError("nonbasis has wrong cardinalty.")
                    if not b.issubset(self._groundset):
                        raise ValueError("nonbasis is not a subset of the groundset")
                    self.__pack(self._b, b)
                    i = set_to_index(self._b)
                    if bitset_in(self._bb, i):
                        self._bcount -= 1
                    bitset_discard(self._bb, i)

        self.reset_current_basis()

    def __dealloc__(self):
            bitset_free(self._b)
            bitset_free(self._bb)

    # Sage special functions
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Fano())
            sage: repr(M)  # indirect doctest
            'Matroid of rank 3 on 7 elements with 28 bases'

        """
        return Matroid._repr_(self) + " with " + str(self.bases_count()) + " bases"

    # support for parent BasisExchangeMatroid

    cdef bint __is_exchange_pair(self, long x, long y) except -1:      # test if current_basis-x + y is a basis
        """
        Test if `B-e + f` is a basis of the current matroid.

        Here ``B`` is the 'current' basis, i.e. the one returned by
        ``self.basis()``, and ``e=self._E[x]``, ``f=self._E[y]``.

        INPUT:

        - ``x`` -- an integer
        - ``y`` -- an integer

        OUTPUT:

        ``True`` if `B-e + f` is a basis, ``False`` otherwise. Here `e`, `f`
        are the groundset elements which are internally named by the integers
        ``x``, ``y``.

        NOTE: this is an internal function, supporting the parent class
        BasisExchangeMatroid of BasisMatroid.

        """
        bitset_copy(self._b, self._current_basis)
        bitset_discard(self._b, x)
        bitset_add(self._b, y)
        return bitset_in(self._bb, set_to_index(self._b))

    cdef reset_current_basis(self):
        """
        Set the current basis to the (lexicographically) first basis of the
        matroid.
        """
        index_to_set(self._current_basis, bitset_first(self._bb), self._matroid_rank, self._groundset_size)  # set current basis of parent BasisExchangeMatroid

    # a function that is very efficient for this class

    cpdef _is_basis(self, X):
        """
        Test if input is a basis.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        .. WARNING::

            This method assumes that ``X`` has the right size to be a basis,
            i.e. ``len(X) == self.full_rank()``. Otherwise its behavior is
            undefined.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = Matroid(bases=matroids.named_matroids.Vamos().bases())
            sage: M._is_basis(set(['a', 'b', 'c', 'e']))
            True
            sage: M._is_basis(set(['a', 'b', 'c', 'd']))
            False
        """
        self.__pack(self._b, X)
        return bitset_in(self._bb, set_to_index(self._b))

    # dual and minors

    cpdef dual(self):
        """
        Return the dual of the matroid.

        Let `M` be a matroid with ground set `E`. If `B` is the set of bases
        of `M`, then the set `\{E - b : b \in B\}` is the set of bases of
        another matroid, the *dual* of `M`.

        OUTPUT:

        The dual matroid.

        EXAMPLES::

            sage: M = Matroid(bases=matroids.named_matroids.Pappus().bases())
            sage: M.dual()
            Matroid of rank 6 on 9 elements with 75 bases

        ALGORITHM:

        A BasisMatroid on `n` elements and of rank `r` is stored as a
        bitvector of length `n \choose r`. The `i`-th bit in this vector
        indicates that the `i`-th `r`-set in the lexicographic enumeration of
        `r`-subsets of the groundset is a basis. Reversing this bitvector
        yields a bitvector that indicates whether the complement of an
        `(n - r)`-set is a basis, i.e. gives the bitvector of the bases of the
        dual.

        """
        cdef long i, N
        cdef BasisMatroid D
        D = BasisMatroid(groundset=self._E, rank=self.full_corank())
        N = binom[self._groundset_size][self._matroid_rank]
        for i in xrange(N):
            if not bitset_in(self._bb, i):
                bitset_discard(D._bb, N - i - 1)
        D.reset_current_basis()
        D._reset_invariants()
        D._bcount = self._bcount
        return D

    cpdef _minor(self, contractions, deletions):
        """
        Return a minor.

        INPUT:

        - ``contractions`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.
        - ``deletions`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above,
            we assume the following:

            - ``contractions`` is independent
            - ``deletions`` is coindependent
            - ``contractions`` and ``deletions`` are disjoint.

        OUTPUT:

        A matroid.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Vamos())
            sage: M._minor(contractions=set(['a']), deletions=set(['b', 'c']))
            Matroid of rank 3 on 5 elements with 10 bases

        """
        E = self.groundset() - (contractions | deletions)
        mr = self.full_rank() - len(contractions)
        NB = [frozenset(B) for B in combinations(E, mr) if not self._is_basis(contractions | frozenset(B))]
        return BasisMatroid(groundset=E, nonbases=NB, rank=mr)

    cpdef truncation(self):
        r"""
        Return a rank-1 truncation of the matroid.

        Let `M` be a matroid of rank `r`. The *truncation* of `M` is the
        matroid obtained by declaring all subsets of size `r` dependent. It
        can be obtained by adding an element freely to the span of the matroid
        and then contracting that element.

        OUTPUT:

        A matroid.

        .. SEEALSO::

            :meth:`M.extension() <sage.matroids.matroid.Matroid.extension>`,
            :meth:`M.contract() <sage.matroids.matroid.Matroid.contract>`

        EXAMPLES::

            sage: M = Matroid(bases=matroids.named_matroids.N2().bases())
            sage: M.truncation()
            Matroid of rank 5 on 12 elements with 702 bases
            sage: M.f_vector()
            [1, 12, 66, 190, 258, 99, 1]
            sage: M.truncation().f_vector()
            [1, 12, 66, 190, 258, 1]

        """
        if self.full_rank() == 0:
            return None
        return BasisMatroid(groundset=self._E, nonbases=self.dependent_r_sets(self.full_rank() - 1), rank=self.full_rank() - 1)

    cpdef _extension(self, e, H):
        r"""
        Extend the matroid by a new element.

        The result is a matroid on ``self.groundset() + {element}``, where
        ``element`` is contained in exactly the hyperplanes of ``self``
        specified by ``hyperplanes``.

        INPUT:

        - ``element`` -- a hashable object not in ``self.groundset()``.
        - ``hyperplanes`` -- the set of hyperplanes of a linear subclass of
          ``self``.

        OUTPUT:

        A matroid.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.Uniform(3, 5))
            sage: H = M.hyperplanes()
            sage: M._extension('x', [H[0]])
            Matroid of rank 3 on 6 elements with 19 bases
            sage: M._extension('x', [H[0], H[1]])
            Matroid of rank 3 on 6 elements with 18 bases
            sage: M._extension('x', H)
            Matroid of rank 3 on 6 elements with 10 bases
            sage: len([M._extension('x', mc) for mc in M.linear_subclasses()])
            32
        """
        cdef bint found
        cdef frozenset B
        if self.full_rank() == 0:
            return BasisMatroid(groundset=self._E + (e,), bases=[set()])

        BB = self.bases()
        BT = self.independent_r_sets(self.full_rank() - 1)
        se = set([e])
        BE = []
        for B in BT:
            found = False
            for hyp in H:
                if B.issubset(hyp):
                    found = True
                    break
            if not found:
                BE.append(B | se)
        BE += BB
        return BasisMatroid(groundset=self._E + (e,), bases=BE)

    cpdef _with_coloop(self, e):
        r"""
        Return the matroid that arises by adding an element `e` to the
        groundset, that is a coloop of the resulting matroid.

        INPUT:

        - ``e`` -- the label of the new element. Assumed to be outside the
          current groundset.

        OUTPUT:

        The extension of this matroid by a coloop.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Fano())
            sage: M
            Matroid of rank 3 on 7 elements with 28 bases
            sage: M._with_coloop('x')
            Matroid of rank 4 on 8 elements with 28 bases

        """
        cdef frozenset se = frozenset([e])
        return BasisMatroid(groundset=self._E + (e,), bases=[B | se for B in self.bases()])

    cpdef relabel(self, l):
        """
        Return an isomorphic matroid with relabeled groundset.

        The output is obtained by relabeling each element ``e`` by ``l[e]``,
        where ``l`` is a given injective map. If ``e not in l`` then the
        identity map is assumed.

        INPUT:

        - ``l`` -- a python object such that `l[e]` is the new label of `e`.

        OUTPUT:

        A matroid.

        .. TODO::

            Write abstract ``RelabeledMatroid`` class, and add ``relabel()``
            method to the main ``Matroid`` class, together with ``_relabel()``
            method that can be replaced by subclasses. Use the code from
            ``is_isomorphism()`` in ``relabel()`` to deal with a variety of
            input methods for the relabeling.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Fano())
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: N = M.relabel({'g':'x'})
            sage: sorted(N.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'x']

        """
        M = BasisMatroid(M=self)
        M.__relabel(l)
        return M

    # enumeration

    cpdef bases_count(self):
        r"""
        Return the number of bases of the matroid.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: M = Matroid(bases=matroids.named_matroids.Fano().bases())
            sage: M
            Matroid of rank 3 on 7 elements with 28 bases
            sage: M.bases_count()
            28
        """
        if self._bcount is None:
            self._bcount = bitset_len(self._bb)
        return self._bcount

    cpdef bases(self):
        r"""
        Return the list of bases of the matroid.

        A *basis* is a maximal independent set.

        OUTPUT:

        An iterable containing all bases of the matroid.

        EXAMPLES::

            sage: M = Matroid(bases=matroids.named_matroids.Fano().bases())
            sage: M
            Matroid of rank 3 on 7 elements with 28 bases
            sage: len(M.bases())
            28
        """
        cdef long r, n
        r = self.full_rank()
        n = len(self)
        cdef SetSystem BB
        BB = SetSystem(self._E, capacity=bitset_len(self._bb))
        cdef long b
        b = bitset_first(self._bb)
        while b >= 0:
            index_to_set(self._b, b, r, n)
            BB._append(self._b)
            b = bitset_next(self._bb, b + 1)
        return BB

    cpdef nonbases(self):
        r"""
        Return the list of nonbases of the matroid.

        A *nonbasis* is a set with cardinality ``self.full_rank()`` that is
        not a basis.

        OUTPUT:

        An iterable containing the nonbases of the matroid.

        .. SEEALSO::

            :meth:`Matroid.basis() <sage.matroids.matroid.Matroid.basis>`

        EXAMPLES::

            sage: M = Matroid(bases=matroids.named_matroids.Fano().bases())
            sage: M
            Matroid of rank 3 on 7 elements with 28 bases
            sage: len(M.nonbases())
            7
        """
        if self._nonbases is not None:
            return self._nonbases
        cdef long r, n
        r = self.full_rank()
        n = len(self)
        cdef bitset_t bb_comp
        bitset_init(bb_comp, binom[self._groundset_size][self._matroid_rank])
        bitset_complement(bb_comp, self._bb)
        cdef SetSystem NB
        NB = SetSystem(self._E, capacity=bitset_len(bb_comp))
        cdef long b
        b = bitset_first(bb_comp)
        while b >= 0:
            index_to_set(self._b, b, r, n)
            NB._append(self._b)
            b = bitset_next(bb_comp, b + 1)
        bitset_free(bb_comp)
        self._nonbases = NB
        return NB

    # isomorphism test

    cpdef _bases_invariant(self):
        """
        Return an isomorphism invariant based on the incidences of groundset
        elements with bases.

        INPUT:

        - Nothing

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Fano())
            sage: N = BasisMatroid(matroids.named_matroids.Fano())
            sage: M._bases_invariant() == N._bases_invariant()
            True
        """
        if self._bases_invariant_var is not None:
            return self._bases_invariant_var
        cdef long i, j
        cdef bitset_t bb_comp
        cdef list bc
        cdef dict bi
        bc = [0 for i in xrange(len(self))]
        for i in xrange(binom[self._groundset_size][self._matroid_rank]):
            if not bitset_in(self._bb, i):
                index_to_set(self._b, i, self._matroid_rank, self._groundset_size)
                j = bitset_first(self._b)
                while j >= 0:
                    bc[j] += 1
                    j = bitset_next(self._b, j + 1)
        bi = {}
        for e in xrange(len(self)):
            if bc[e] in bi:
                bi[bc[e]].append(e)
            else:
                bi[bc[e]] = [e]
        self._bases_invariant_var = hash(tuple([(c, len(bi[c])) for c in sorted(bi)]))
        self._bases_partition_var = SetSystem(self._E, [[self._E[e] for e in bi[c]] for c in sorted(bi)])
        return self._bases_invariant_var

    cpdef _bases_partition(self):
        """
        Return an ordered partition based on the incidences of groundset
        elements with bases.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Vamos())
            sage: [sorted(p) for p in M._bases_partition()]
            [['c', 'd', 'g', 'h'], ['a', 'b', 'e', 'f']]
        """
        self._bases_invariant()
        return self._bases_partition_var

    cpdef _bases_invariant2(self):
        """
        Return an isomophism invariant of the matroid.

        Compared to BasisMatroid._bases_invariant() this invariant
        distinguishes more frequently between nonisomorphic matroids but
        takes more time to compute.
        See also :meth:`<BasisMatroid.basis_partition2>`.

        OUTPUT:

        an integer isomorphism invariant.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Fano())
            sage: N = BasisMatroid(matroids.named_matroids.NonFano())
            sage: M._bases_invariant2() == N._bases_invariant2()
            False
        """
        if self._bases_invariant2_var is None:
            CP = self.nonbases()._equitable_partition(self._bases_partition())
            self._bases_partition2_var = CP[0]
            self._bases_invariant2_var = CP[2]
        return self._bases_invariant2_var

    cpdef _bases_partition2(self):
        """
        Return an equitable partition which refines
        :meth:`<BasisMatroid._bases_partition2>`.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Vamos())
            sage: [sorted(p) for p in M._bases_partition2()]
            [['c', 'd', 'g', 'h'], ['a', 'b', 'e', 'f']]
        """
        self._bases_invariant2()
        return self._bases_partition2_var

    cpdef _bases_invariant3(self):
        """
        Return a number characteristic for the construction of
        :meth:`<BasisMatroid._bases_partition3>`.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Vamos())
            sage: N = BasisMatroid(matroids.named_matroids.Vamos())
            sage: M._bases_invariant3() == N._bases_invariant3()
            True
        """
        if self._bases_invariant3_var is None:
            CP = self.nonbases()._heuristic_partition(self._bases_partition2())
            self._bases_partition3_var = CP[0]
            self._bases_invariant3_var = CP[2]
        return self._bases_invariant3_var

    cpdef _bases_partition3(self):
        """
        Return an ordered partition into singletons which refines an equitable
        partition of the matroid.

        The purpose of this partition is to heuristically find an isomorphism
        between two matroids, by lining up their respective
        heuristic_partitions.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Vamos())
            sage: N = BasisMatroid(matroids.named_matroids.Vamos())
            sage: PM = M._bases_partition3()
            sage: PN = N._bases_partition3()
            sage: morphism = {}
            sage: for i in xrange(len(M)): morphism[min(PM[i])]=min(PN[i])
            sage: M._is_isomorphism(N, morphism)
            True
        """
        self._bases_invariant3()
        return self._bases_partition3_var

    cdef _reset_invariants(self):
        """
        Remove all precomputed invariants.
        """
        self._bcount = None
        self._nonbases = None
        self._bases_invariant_var = None
        self._bases_partition_var = None
        self._bases_invariant2_var = None
        self._bases_partition2_var = None
        self._bases_invariant3_var = None
        self._bases_partition3_var = None
        self._flush()

    cpdef bint is_distinguished(self, e):
        """
        Return whether ``e`` is a 'distinguished' element of the groundset.

        The set of distinguished elements is an isomorphism invariant. Each
        matroid has at least one distinguished element. The typical
        application of this method is the execution of an orderly algorithm
        for generating all matroids up to isomorphism in a minor-closed class,
        by successively enumerating the single-element extensions and
        coextensions of the matroids generated so far.

        INPUT:

        - ``e`` -- an element of the ground set

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.extensions() <sage.matroids.matroid.Matroid.extensions>`,
            :meth:`M.linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`,
            :mod:`sage.matroids.extension <sage.matroids.extension>`

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.N1())
            sage: sorted([e for e in M.groundset() if M.is_distinguished(e)])
            ['c', 'g', 'h', 'j']

        """
        P = self._bases_partition()
        p = P[0]
        if e not in p:
            return False
        if len(p) == 1:
            return True

        SP = self._bases_partition2()
        q = p
        for q2 in SP:
            if q2.issubset(p) and len(q2) < len(q):
                q = q2
        return e in q

    cpdef _is_relaxation(self, other, morphism):
        """
        Return if the application of a groundset morphism to this matroid
        yields a relaxation of the given matroid.

        M is a relaxation of N if the set of bases of M is a subset of the
        bases of N.

        INPUT:

        - ``other`` -- a BasisMatroid
        - ``morphism`` -- a dictionary with sends each element of the
          groundset of this matroid to a distinct element of the groundset
          of ``other``

        OUTPUT:

        ``True`` if ``morphism[self]`` is a relaxation of ``other``;
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.NonFano())
            sage: N = BasisMatroid(matroids.named_matroids.Fano())
            sage: m = {e:e for e in M.groundset()}
            sage: M._is_relaxation(N, m)
            True
            sage: N._is_relaxation(M, m)
            False
        """
        cdef long i, j
        cdef bitset_t b2
        cdef bitset_t bb_comp

        bitset_init(bb_comp, binom[self._groundset_size][self._matroid_rank])
        bitset_complement(bb_comp, self._bb)

        bitset_init(b2, max(len(self), 1))
        morph = [(<BasisMatroid>other)._idx[morphism[self._E[i]]] for i in xrange(len(self))]
        i = bitset_first(bb_comp)
        while i != -1:
            index_to_set(self._b, i, self._matroid_rank, self._groundset_size)
            bitset_clear(b2)
            j = bitset_first(self._b)
            while j != -1:
                bitset_add(b2, morph[j])
                j = bitset_next(self._b, j + 1)
            if bitset_in((<BasisMatroid>other)._bb, set_to_index(b2)):
                return False
            i = bitset_next(bb_comp, i + 1)
        return True

    cpdef _is_isomorphism(self, other, morphism):
        """
        Version of is_isomorphism() that does no type checking.

        INPUT:

        - ``other`` -- A matroid instance.
        - ``morphism`` -- a dictionary mapping the groundset of ``self`` to
          the groundset of ``other``.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`<sage.matroids.matroid.Matroid.is_isomorphism>`.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.NonFano())
            sage: N = BasisMatroid(matroids.named_matroids.Fano())
            sage: m = {e:e for e in M.groundset()}
            sage: M._is_relaxation(N, m)
            True
            sage: M._is_isomorphism(N, m)
            False
        """
        if not isinstance(other, BasisMatroid):
            ot = BasisMatroid(other)
        else:
            ot = other
        return self.bases_count() == (<BasisMatroid>ot).bases_count() and self._is_relaxation(ot, morphism)

    cpdef _isomorphism(self, other):
        """
        Return isomorphism from ``self`` to ``other``, if one exists.

        INPUT:

        - ``other`` -- a matroid.

        OUTPUT:

        A dictionary, or ``None``

        .. NOTE::

            Internal version that does no input checking.

        EXAMPLES::
        
            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.Wheel(3))
            sage: N = BasisMatroid(matroids.CompleteGraphic(4))
            sage: morphism = M._isomorphism(N)
            sage: M._is_isomorphism(N, morphism)
            True
            sage: M = BasisMatroid(matroids.named_matroids.NonFano())
            sage: N = BasisMatroid(matroids.named_matroids.Fano())
            sage: M._isomorphism(N) is not None
            False
        """
        if not isinstance(other, BasisMatroid):
            return self.isomorphism(BasisMatroid(other))
        if self is other:
            return {e:e for e in self.groundset()}
        if len(self) != len(other):
            return None
        if self.full_rank() != other.full_rank():
            return None
        if self.full_rank() == 0:
            return {self.groundset_list()[i]: other.groundset_list()[i] for i in xrange(len(self))}        
        if self.bases_count() != other.bases_count():
            return None

        if self._bases_invariant() != other._bases_invariant():
            return None
        PS = self._bases_partition()
        PO = other._bases_partition()
        if len(PS) == len(self) and len(PO) == len(other):
            morphism = {}
            for i in xrange(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            if self._is_relaxation(other, morphism):
                return morphism
            else:
                return None

        if self._bases_invariant2() != other._bases_invariant2():
            return None
        PS = self._bases_partition2()
        PO = other._bases_partition2()
        if len(PS) == len(self) and len(PO) == len(other):
            morphism = {}
            for i in xrange(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            if self._is_relaxation(other, morphism):
                return morphism
            else:
                return None

        if self._bases_invariant3() == other._bases_invariant3():
            PHS = self._bases_partition3()
            PHO = other._bases_partition3()
            morphism = {}
            for i in xrange(len(self)):
                morphism[min(PHS[i])] = min(PHO[i])
            if self._is_relaxation(other, morphism):
                return morphism

        return self.nonbases()._isomorphism(other.nonbases(), PS, PO)
        
    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Return if this matroid is isomorphic to the given matroid.

        INPUT:

        - ``other`` -- A matroid,
        - optional parameter ``certificate`` -- Boolean.

        OUTPUT:

        Boolean,
        and, if certificate = True, a dictionary giving the isomophism or None

        .. NOTE::

            Internal version that does no input checking.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.NonFano())
            sage: N = BasisMatroid(matroids.named_matroids.Fano())
            sage: M._is_isomorphic(N)
            False
            sage: M._is_isomorphic(N, certificate=True)
            (False, None)
        """
        if certificate:
            return self._is_isomorphic(other), self._isomorphism(other)
        if not isinstance(other, BasisMatroid):
            return BasisExchangeMatroid._is_isomorphic(self, other)
        if self is other:
            return True
        if len(self) != len(other):
            return False
        if self.full_rank() != other.full_rank():
            return False
        if self.full_rank() == 0:
            return True
        if self.bases_count() != other.bases_count():
            return False
        if self.full_rank() < 2 or self.full_corank() < 2:
            return True  # number of bases then determines matroid up to isomorphism

        if self._bases_invariant() != other._bases_invariant():
            return False
        PS = self._bases_partition()
        PO = other._bases_partition()
        if len(PS) == len(self) and len(PO) == len(other):
            morphism = {}
            for i in xrange(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            return self._is_relaxation(other, morphism)

        if self._bases_invariant2() != other._bases_invariant2():
            return False
        PS = self._bases_partition2()
        PO = other._bases_partition2()
        if len(PS) == len(self) and len(PO) == len(other):
            morphism = {}
            for i in xrange(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            return self._is_relaxation(other, morphism)

        if self._bases_invariant3() == other._bases_invariant3():
            PHS = self._bases_partition3()
            PHO = other._bases_partition3()
            morphism = {}
            for i in xrange(len(self)):
                morphism[min(PHS[i])] = min(PHO[i])
            if self._is_relaxation(other, morphism):
                return True

        return self.nonbases()._isomorphism(other.nonbases(), PS, PO) is not None

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to __richcmp__ (in Cython) and __cmp__ or
            __eq__/__ne__ (in Python). If you override one, you should (and in
            Cython: MUST) override the other!

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Fano())
            sage: N = BasisMatroid(matroids.named_matroids.Fano().dual()).dual()
            sage: O = BasisMatroid(matroids.named_matroids.NonFano())
            sage: hash(M) == hash(N)
            True
            sage: hash(M) == hash(O)
            False
        """
        return hash((self.groundset(), self.bases_count(), self._weak_invariant()))

    def __richcmp__(left, right, int op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. For BasisMatroids, this means that the groundsets and the
        sets of bases of the two matroids are equal.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Pappus())
            sage: N = BasisMatroid(matroids.named_matroids.NonPappus())
            sage: M == N
            False
        """
        if op in [0, 1, 4, 5]:  # <, <=, >, >=
            return NotImplemented
        if not isinstance(left, BasisMatroid) or not isinstance(right, BasisMatroid):
            return NotImplemented
        if op == 2:  # ==
            res = True
        if op == 3:  # !=
            res = False
        # res gets inverted if matroids are deemed different.
        if left.equals(right):
            return res
        else:
            return not res

    def __copy__(self):
        """
        Create a shallow copy.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Vamos())
            sage: N = copy(M)  # indirect doctest
            sage: M == N
            True

        """
        N = BasisMatroid(M=self)
        N.rename(getattr(self, '__custom_name'))
        return N

    def __deepcopy__(self, memo={}):
        """
        Create a deep copy.

        .. NOTE::

            Identical to shallow copy for BasisMatroid class.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Vamos())
            sage: N = deepcopy(M)  # indirect doctest
            sage: M == N
            True
            sage: M.groundset() is N.groundset()
            False
        """
        N = BasisMatroid(M=self)
        N.rename(getattr(self, '__custom_name'))
        return N

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle, (version, data))``, where ``unpickle`` is the
        name of a function that, when called with ``(version, data)``,
        produces a matroid isomorphic to ``self``. ``version`` is an integer
        (currently 0) and ``data`` is a tuple ``(E, R, name, BB)`` where
        ``E`` is the groundset of the matroid, ``R`` is the rank, ``name`` is a
        custom name, and ``BB`` is the bitpacked list of bases, as pickled by
        Sage's ``bitset_pickle``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.named_matroids.Vamos())
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.rename('Vamos')
            sage: loads(dumps(M))
            Vamos
        """
        import sage.matroids.unpickling
        BB = bitset_pickle(self._bb)
        data = (self._E, self._matroid_rank, getattr(self, '__custom_name'), BB)
        version = 0
        return sage.matroids.unpickling.unpickle_basis_matroid, (version, data)


cdef long binom[2956][33]   # Cached binomial table


cdef  binom_init(long N, long K):
    """
    Fill up the cached binomial table.
    """
    cdef long bin
    if binom[0][0] != 1:
        binom[0][0] = 1
        binom[0][1] = 0
        for n in xrange(1, 2955):
            bin = 1
            k = 0
            while bin < 2 ** 32 and k <= 32 and k <= n:
                binom[n][k] = bin
                k += 1
                bin = binom[n - 1][k - 1] + binom[n - 1][k]
            while k < 33:
                binom[n][k] = 0
                k += 1

    if N > 2954:
        raise ValueError("BasisMatroid: size of groundset exceeds 2954")  # if n > 2954 and k > 2, then binomial(n, k) > 2^32
    if K > 32:
        raise ValueError("BasisMatroid: rank exceeds 32")  # if n > 2954 and k > 2, then binomial(n, k) > 2^32
    if binom[N][K] == 0:
        raise ValueError("BasisMatroid: number of potential bases would exceed 2^32")

cdef long set_to_index(bitset_t S):
    """
    Compute the rank of a set of integers amongst the sets of integers
    of the same cardinality.
    """
    cdef long index = 0
    cdef long count = 1
    cdef long s
    s = bitset_first(S)
    while s >= 0:
        index += binom[s][count]
        count += 1
        s = bitset_next(S, s + 1)
    return index

cdef  index_to_set(bitset_t S, long index, long k, long n):
    """
    Compute the k-subset of `\{0, ..., n-1\}` of rank index
    """
    bitset_clear(S)
    cdef long s = n
    while s > 0:
        s -= 1
        if binom[s][k] <= index:
            index = index - binom[s][k]
            k = k - 1
            bitset_add(S, s)
