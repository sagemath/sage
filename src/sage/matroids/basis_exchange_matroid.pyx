"""
Basis exchange matroids

:class:`BasisExchangeMatroid <sage.matroids.basis_exchange_matroid.BasisExchangeMatroid>`
is an abstract class implementing default matroid functionality in terms of
basis exchange. Several concrete matroid classes are subclasses of this. They
have the following methods in addition to the ones provided by the parent
class :mod:`Matroid <sage.matroids.matroid>`.

- :func:`bases_count() <sage.matroids.basis_exchange_matroid.BasisExchangeMatroid.bases_count>`
- :func:`groundset_list() <sage.matroids.basis_exchange_matroid.BasisExchangeMatroid.groundset_list>`

AUTHORS:

- Rudi Pendavingh, Stefan van Zwam (2013-04-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem
from sage.data_structures.bitset_base cimport *

cdef class BasisExchangeMatroid(Matroid):
    r"""
    Class BasisExchangeMatroid is a virtual class that derives from Matroid.
    It implements each of the elementary matroid methods
    (:meth:`rank() <sage.matroids.matroid.Matroid.rank>`,
    :meth:`max_independent() <sage.matroids.matroid.Matroid.max_independent>`,
    :meth:`circuit() <sage.matroids.matroid.Matroid.circuit>`,
    :meth:`closure() <sage.matroids.matroid.Matroid.closure>` etc.),
    essentially by crawling the base exchange graph of the matroid. This is
    the graph whose vertices are the bases of the matroid, two bases being
    adjacent in the graph if their symmetric difference has 2 members.

    This base exchange graph is not stored as such, but should be provided
    implicitly by the child class in the form of two methods
    ``_is_exchange_pair(x, y)`` and ``_exchange(x, y)``, as well as an
    initial basis. At any moment, BasisExchangeMatroid keeps a current basis
    `B`. The method ``_is_exchange_pair(x, y)`` should return a boolean
    indicating whether `B - x + y` is a basis. The method ``_exchange(x, y)``
    is called when the current basis `B` is replaced by said `B-x + y`. It is
    up to the child class to update its internal data structure to make
    information relative to the new basis more accessible. For instance, a
    linear matroid would perform a row reduction to make the column labeled by
    `y` a standard basis vector (and therefore the columns indexed by `B-x+y`
    would form an identity matrix).

    Each of the elementary matroid methods has a straightforward greedy-type
    implementation in terms of these two methods. For example, given a subset
    `F` of the groundset, one can step to a basis `B` over the edges of the
    base exchange graph which has maximal intersection with `F`, in each step
    increasing the intersection of the current `B` with `F`. Then one computes
    the rank of `F` as the cardinality of the intersection of `F` and `B`.

    The following matroid classes can/will implement their oracle efficiently
    by deriving from ``BasisExchangeMatroid``:

    - :class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`: keeps
      a list of all bases.

        - ``_is_exchange_pair(x, y)`` reduces to a query whether `B - x + y`
          is a basis.
        - ``_exchange(x, y)`` has no work to do.

    - :class:`LinearMatroid <sage.matroids.linear_matroid.LinearMatroid>`:
      keeps a matrix representation `A` of the matroid so that `A[B] = I`.

        - ``_is_exchange_pair(x, y)`` reduces to testing whether `A[r, y]`
          is nonzero, where `A[r, x]=1`.
        - ``_exchange(x, y)`` should modify the matrix so that `A[B - x + y]`
          becomes `I`, which means pivoting on `A[r, y]`.

    - ``TransversalMatroid`` (not yet implemented): If `A` is a set of subsets
      of `E`, then `I` is independent if it is a system of distinct
      representatives of `A`, i.e. if `I` is covered by a matching of an
      appropriate bipartite graph `G`, with color classes `A` and `E` and an
      edge `(A_i,e)` if `e` is in the subset `A_i`. At any time you keep a
      maximum matching `M` of `G` covering the current basis `B`.

        - ``_is_exchange_pair(x, y)`` checks for the existence of an
          `M`-alternating path `P` from `y` to `x`.
        - ``_exchange(x, y)`` replaces `M` by the symmetric difference of
          `M` and `E(P)`.

    - ``AlgebraicMatroid`` (not yet implemented): keeps a list of polynomials
      in variables `E - B + e` for each variable `e` in `B`.

        - ``_is_exchange_pair(x, y)`` checks whether the polynomial that
          relates `y` to `E-B` uses `x`.
        - ``_exchange(x, y)`` make new list of polynomials by computing
          resultants.

    All but the first of the above matroids are algebraic, and all
    implementations specializations of the algebraic one.

    BasisExchangeMatroid internally renders subsets of the groundset as
    bitsets. It provides optimized methods for enumerating bases, nonbases,
    flats, circuits, etc.
    """

    def __init__(self, groundset, basis=None, rank=None):
        """
        Construct a BasisExchangeMatroid.

        A BasisExchangeMatroid is a virtual class. It is unlikely that you
        want to create a BasisExchangeMatroid from the command line. See the
        class docstring for
        :class:`BasisExchangeMatroid <sage.matroids.basis_exchange_matroid.BasisExchangeMatroid>`.

        INPUT:

        - ``groundset`` -- set
        - ``basis`` -- subset of groundset (default: ``None``)
        - ``rank`` -- integer (default: ``None``)

        This initializer sets up a correspondence between elements of
        ``groundset`` and ``range(len(groundset))``. ``BasisExchangeMatroid``
        uses this correspondence for encoding of subsets of the groundset as
        bitpacked sets of integers --- see ``_pack()`` and ``__unpack()``. In
        general, methods of ``BasisExchangeMatroid`` having a name starting
        with two underscores deal with such encoded subsets.

        A second task of this initializer is to store the rank and initialize
        the 'current' basis.

        EXAMPLES::

            sage: from sage.matroids.basis_exchange_matroid import BasisExchangeMatroid
            sage: M = BasisExchangeMatroid(groundset='abcdef', basis='abc')
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f']
            sage: sorted(M.basis())
            ['a', 'b', 'c']


        TESTS::

            sage: from sage.matroids.advanced import *
            sage: M = BasisExchangeMatroid(groundset=[1, 2, 3], rank=2)
            sage: TestSuite(M).run(skip='_test_pickling')

        .. NOTE::

            This is an abstract base class, without a data structure, so no
            pickling mechanism was implemented.
        """
        self._groundset_size = len(groundset)
        self._bitset_size = max(self._groundset_size, 1)
        if basis is None:
            self._matroid_rank = rank
        else:
            self._matroid_rank = len(basis)
        bitset_init(self._current_basis, self._bitset_size)
        bitset_init(self._inside, self._bitset_size)
        bitset_init(self._outside, self._bitset_size)
        bitset_init(self._input, self._bitset_size)
        bitset_init(self._input2, self._bitset_size)
        bitset_init(self._output, self._bitset_size)
        bitset_init(self._temp, self._bitset_size)

        self._groundset = frozenset(groundset)
        if not isinstance(groundset, tuple):
            self._E = tuple(groundset)
        else:
            self._E = groundset
        self._idx = {}
        cdef long i
        for i in range(self._groundset_size):
            self._idx[self._E[i]] = i

        if basis is not None:
            self._pack(self._current_basis, frozenset(basis))

    def __dealloc__(self):
        bitset_free(self._current_basis)
        bitset_free(self._inside)
        bitset_free(self._outside)
        bitset_free(self._input)
        bitset_free(self._input2)
        bitset_free(self._output)
        bitset_free(self._temp)

    cdef _relabel(self, mapping):
        """
        Relabel each element `e` as ``mapping[e]``, where ``mapping`` is a
        given injective map.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: none

        .. NOTE::

            For internal use. Matroids are immutable but this method does
            modify the matroid. The use of this method will only be safe in
            very limited circumstances, such as perhaps on a fresh copy of a
            matroid.
        """
        cdef long i
        E = []
        for i in range(self._groundset_size):
            try:
                E.append(mapping[self._E[i]])
            except LookupError:
                E.append(self._E[i])
        self._E = tuple(E)
        self._groundset = frozenset(E)

        self._idx = {}
        for i in range(self._groundset_size):
            self._idx[self._E[i]] = i

        if self._weak_partition_var:
            self._weak_partition_var._relabel(mapping)

        if self._strong_partition_var:
            self._strong_partition_var._relabel(mapping)
        if self._heuristic_partition_var:
            self._heuristic_partition_var._relabel(mapping)

    # the engine

    cdef _pack(self, bitset_t I, F):
        """
        Encode a subset F of the groundset into a bitpacked set of integers.
        """
        bitset_clear(I)
        for f in F:
            bitset_add(I, <mp_bitcnt_t> self._idx[f])

    cdef __unpack(self, bitset_t I):
        """
        Unencode a bitpacked set of integers to a subset of the groundset.
        """
        cdef long i
        F = set()
        i = bitset_first(I)
        while i >= 0:
            F.add(self._E[i])
            i = bitset_next(I, i + 1)
        return frozenset(F)

    # this method needs to be overridden by child class
    cdef bint _is_exchange_pair(self, long x, long y) except -1:
        """
        Test if current_basis-x + y is a basis.
        """
        raise NotImplementedError

    # if this method is overridden by a child class, the child class needs to call this method
    cdef int _exchange(self, long x, long y) except -1:
        """
        Put ``current_basis <-- current_basis-x + y``.
        """
        bitset_discard(self._current_basis, x)
        bitset_add(self._current_basis, y)

    cdef int _move(self, bitset_t X, bitset_t Y) except -1:
        """
        Change ``current_basis`` to minimize intersection with ``X``, maximize
        intersection with ``Y``.
        """
        cdef long x, y
        x = bitset_first(X)
        while x >= 0:
            y = bitset_first(Y)
            while y >= 0:
                if self._is_exchange_pair(x, y):
                    self._exchange(x, y)
                    bitset_discard(Y, y)
                    bitset_discard(X, x)
                    if bitset_isempty(Y):
                        return 0
                    break
                else:
                    y = bitset_next(Y, y + 1)
            x = bitset_next(X, x + 1)

    cdef __fundamental_cocircuit(self, bitset_t C, long x):
        """
        Return the unique cocircuit that meets ``self._current_basis`` in
        exactly element ``x``.
        """
        cdef long y
        bitset_clear(C)
        bitset_complement(self._temp, self._current_basis)
        y = bitset_first(self._temp)
        while y >= 0:
            if self._is_exchange_pair(x, y):
                bitset_add(C, y)
            y = bitset_next(self._temp, y + 1)
        bitset_add(C, x)

    cdef __fundamental_circuit(self, bitset_t C, long y):
        """
        Return the unique circuit contained in ``self._current_basis`` union ``y``.
        """
        cdef long x
        bitset_clear(C)
        x = bitset_first(self._current_basis)
        while x >= 0:
            if self._is_exchange_pair(x, y):
                bitset_add(C, x)
            x = bitset_next(self._current_basis, x + 1)
        bitset_add(C, y)

    cdef __max_independent(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``max_independent``.
        """
        bitset_difference(self._inside, self._current_basis, F)
        bitset_difference(self._outside, F, self._current_basis)
        self._move(self._inside, self._outside)
        bitset_intersection(R, self._current_basis, F)

    cdef __circuit(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``circuit``.
        """
        bitset_difference(self._inside, self._current_basis, F)
        bitset_difference(self._outside, F, self._current_basis)
        cdef long x, y
        y = bitset_first(self._outside)
        if y < 0:
            raise ValueError('no circuit in independent set.')
        while y >= 0:
            x = bitset_first(self._inside)
            while x >= 0:
                if self._is_exchange_pair(x, y):
                    self._exchange(x, y)
                    bitset_discard(self._outside, y)
                    bitset_discard(self._inside, x)
                    if bitset_isempty(self._outside):
                        raise ValueError('no circuit in independent set.')
                    break
                else:
                    x = bitset_next(self._inside, x + 1)
            if x == -1:
                self.__fundamental_circuit(R, y)
                return
            y = bitset_next(self._outside, y + 1)

    cdef __closure(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``closure``.
        """
        bitset_difference(self._inside, self._current_basis, F)
        bitset_difference(self._outside, F, self._current_basis)
        self._move(self._inside, self._outside)
        bitset_set_first_n(R, self._groundset_size)
        cdef long x = bitset_first(self._inside)
        while x >= 0:
            self.__fundamental_cocircuit(F, x)
            bitset_difference(R, R, F)
            x = bitset_next(self._inside, x + 1)

    cdef __max_coindependent(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``max_coindependent``.
        """
        bitset_complement(R, F)
        bitset_difference(self._inside, self._current_basis, R)
        bitset_difference(self._outside, R, self._current_basis)
        self._move(self._inside, self._outside)
        bitset_difference(R, F, self._current_basis)

    cdef __cocircuit(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``cocircuit``.
        """
        bitset_complement(R, F)
        bitset_difference(self._inside, self._current_basis, R)
        bitset_difference(self._outside, R, self._current_basis)
        cdef long x, y
        x = bitset_first(self._inside)
        if x < 0:
            raise ValueError('no cocircuit in coindependent set.')
        while x >= 0:
            y = bitset_first(self._outside)
            while y >= 0:
                if self._is_exchange_pair(x, y):
                    self._exchange(x, y)
                    bitset_discard(self._outside, y)
                    bitset_discard(self._inside, x)
                    if bitset_isempty(self._inside):
                        raise ValueError('no cocircuit in coindependent set.')
                    break
                else:
                    y = bitset_next(self._outside, y + 1)
            if y == -1:
                self.__fundamental_cocircuit(R, x)
                return
            x = bitset_next(self._inside, x + 1)

    cdef _coclosure_internal(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``closure``.
        """
        bitset_complement(R, F)
        bitset_difference(self._inside, self._current_basis, R)
        bitset_difference(self._outside, R, self._current_basis)
        self._move(self._inside, self._outside)
        bitset_set_first_n(R, self._groundset_size)
        cdef long y = bitset_first(self._outside)
        while y >= 0:
            self.__fundamental_circuit(F, y)
            bitset_difference(R, R, F)
            y = bitset_next(self._outside, y + 1)

    cdef __augment(self, bitset_t R, bitset_t X, bitset_t Y):
        """
        Bitpacked version of ``augment``.
        """
        bitset_difference(self._inside, self._current_basis, X)
        bitset_difference(self._outside, X, self._current_basis)
        self._move(self._inside, self._outside)
        bitset_difference(self._inside, self._inside, Y)
        bitset_difference(self._outside, Y, self._current_basis)
        self._move(self._inside, self._outside)
        bitset_intersection(R, self._current_basis, Y)

    cdef bint __is_independent(self, bitset_t F) except -1:
        """
        Bitpacked version of ``is_independent``.
        """
        bitset_difference(self._inside, self._current_basis, F)
        bitset_difference(self._outside, F, self._current_basis)
        self._move(self._inside, self._outside)
        return bitset_isempty(self._outside)

    cdef __move_current_basis(self, bitset_t X, bitset_t Y):
        """
        Bitpacked version of ``_move_current_basis``.
        """
        bitset_difference(self._inside, self._current_basis, X)
        bitset_difference(self._outside, X, self._current_basis)
        self._move(self._inside, self._outside)
        bitset_intersection(self._inside, self._current_basis, Y)
        bitset_complement(self._outside, self._current_basis)
        bitset_difference(self._outside, self._outside, Y)
        self._move(self._inside, self._outside)

    # functions for derived classes and for parent class

    cdef bint _set_current_basis(self, F) noexcept:
        """
        Set ``_current_basis`` to subset of the groundset ``F``.
        """
        self._pack(self._input, F)
        bitset_difference(self._inside, self._current_basis, self._input)
        bitset_difference(self._outside, self._input, self._current_basis)
        self._move(self._inside, self._outside)
        return bitset_isempty(self._outside) and bitset_isempty(self._inside)

    # groundset and full_rank

    cpdef frozenset groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT: set

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        """
        return self._groundset

    cpdef list groundset_list(self):
        """
        Return a list of elements of the groundset of the matroid.

        The order of the list does not change between calls.

        OUTPUT: list

        .. SEEALSO::

            :meth:`M.groundset() <sage.matroids.basis_exchange_matroid.BasisExchangeMatroid.groundset>`

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: type(M.groundset())
            <... 'frozenset'>
            sage: type(M.groundset_list())
            <... 'list'>
            sage: sorted(M.groundset_list())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']

            sage: E = M.groundset_list()
            sage: E.remove('a')
            sage: sorted(M.groundset_list())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        """
        return list(self._E)

    def __len__(self):
        """
        Return the size of the groundset.

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: len(M)
            7
            sage: len(M.groundset())
            7
        """
        return self._groundset_size

    cpdef full_rank(self):
        r"""
        Return the rank of the matroid.

        The *rank* of the matroid is the size of the largest independent
        subset of the groundset.

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M.full_rank()
            3
            sage: M.dual().full_rank()
            4
        """
        return self._matroid_rank

    cpdef full_corank(self):
        r"""
        Return the corank of the matroid.

        The *corank* of the matroid equals the rank of the dual matroid. It is
        given by ``M.size() - M.full_rank()``.

        OUTPUT: integer

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.full_rank() <sage.matroids.basis_exchange_matroid.BasisExchangeMatroid.full_rank>`

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M.full_corank()
            4
            sage: M.dual().full_corank()
            3
        """
        return self._groundset_size - self._matroid_rank

    # matroid oracles

    cpdef basis(self):
        r"""
        Return an arbitrary basis of the matroid.

        A *basis* is an inclusionwise maximal independent set.

        .. NOTE::

            The output of this method can change in between calls. It reflects
            the internal state of the matroid. This state is updated by lots
            of methods, including the method ``M._move_current_basis()``.

        OUTPUT: set of elements

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: sorted(M.basis())
            ['a', 'b', 'c']
            sage: M.rank('cd')
            2
            sage: sorted(M.basis())
            ['a', 'c', 'd']
        """
        return self.__unpack(self._current_basis)

    cpdef _move_current_basis(self, X, Y):
        """
        Change current basis so that intersection with X is maximized,
        intersection with Y is minimized.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``
        - ``Y`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: nothing

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: sorted(M.basis())
            ['a', 'b', 'c', 'e']
            sage: M._move_current_basis('ef', 'a')
            sage: sorted(M.basis())
            ['b', 'c', 'e', 'f']
        """
        self._pack(self._input, X)
        self._pack(self._input2, Y)
        self.__move_current_basis(self._input, self._input2)

    cpdef frozenset _max_independent(self, frozenset F):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``F`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: subset of ``F``

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: sorted(M._max_independent(frozenset(['a', 'c', 'd', 'e', 'f'])))
            ['a', 'c', 'd', 'e']

        .. NOTE::

            This is an unguarded method. For the version that verifies if
            the input is indeed a subset of the groundset,
            see :meth:`<sage.matroids.matroid.Matroid.max_independent>`.
        """
        self._pack(self._input, F)
        self.__max_independent(self._output, self._input)
        return self.__unpack(self._output)

    cpdef int _rank(self, frozenset F) except? -1:
        """
        Compute the rank of a subset of the groundset.

        INPUT:

        - ``F`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: M._rank(frozenset(['a', 'c', 'd', 'e', 'f']))
            4

        .. NOTE::

            This is an unguarded method. For the version that verifies if
            the input is indeed a subset of the groundset,
            see :meth:`<sage.matroids.matroid.Matroid.rank>`.
        """
        self._pack(self._input, F)
        self.__max_independent(self._output, self._input)
        return bitset_len(self._output)

    cpdef frozenset _circuit(self, frozenset F):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``F`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: circuit contained in ``F``, if it exists; otherwise, an error
        is raised

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: sorted(sage.matroids.matroid.Matroid._circuit(M,
            ....:                             frozenset(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(sage.matroids.matroid.Matroid._circuit(M,
            ....:                                       frozenset(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set.

        .. NOTE::

            This is an unguarded method. For the version that verifies if
            the input is indeed a subset of the groundset,
            see :meth:`<sage.matroids.matroid.Matroid.circuit>`.
        """
        self._pack(self._input, F)
        self.__circuit(self._output, self._input)
        return self.__unpack(self._output)

    cpdef frozenset _fundamental_circuit(self, frozenset B, e):
        r"""
        Return the `B`-fundamental circuit using `e`.

        Internal version that does no input checking.

        INPUT:

        - ``B`` -- a basis of the matroid
        - ``e`` -- an element not in ``B``

        OUTPUT: set of elements

        EXAMPLES::

            sage: M = matroids.catalog.P8()
            sage: sorted(M._fundamental_circuit(frozenset('abcd'), 'e'))
            ['a', 'b', 'c', 'e']
        """
        self._pack(self._input, B)
        bitset_clear(self._input2)
        self.__move_current_basis(self._input, self._input2)
        self.__fundamental_circuit(self._output, self._idx[e])
        return self.__unpack(self._output)

    cpdef frozenset _closure(self, frozenset F):
        """
        Return the closure of a set.

        INPUT:

        - ``F`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: the smallest closed set containing ``F``

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: sorted(M._closure(frozenset(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']

        .. NOTE::

            This is an unguarded method. For the version that verifies if the
            input is indeed a subset of the groundset, see
            :meth:`<sage.matroids.matroid.Matroid.closure>`.
        """
        self._pack(self._input, F)
        self.__closure(self._output, self._input)
        return self.__unpack(self._output)

    cpdef frozenset _max_coindependent(self, frozenset F):
        """
        Compute a maximal coindependent subset.

        INPUT:

        - ``F`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a maximal coindependent subset of ``F``

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: sorted(M._max_coindependent(frozenset(['a', 'c', 'd', 'e', 'f'])))
            ['a', 'c', 'd', 'f']

        .. NOTE::

            This is an unguarded method. For the version that verifies if the
            input is indeed a subset of the groundset,
            see :meth:`<sage.matroids.matroid.Matroid.max_coindependent>`.
        """
        self._pack(self._input, F)
        self.__max_coindependent(self._output, self._input)
        return self.__unpack(self._output)

    cpdef int _corank(self, frozenset F) noexcept:
        """
        Return the corank of a set.

        INPUT:

        - ``F`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: integer; the corank of `F`

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: M._corank(frozenset(['a', 'e', 'g', 'd', 'h']))
            4

        .. NOTE::

            This is an unguarded method. For the version that verifies if the
            input is indeed a subset of the groundset,
            see :meth:`<sage.matroids.matroid.Matroid.corank>`.
        """
        self._pack(self._input, F)
        self.__max_coindependent(self._output, self._input)
        return bitset_len(self._output)

    cpdef frozenset _cocircuit(self, frozenset F):
        """
        Return a minimal codependent subset.

        INPUT:

        - ``F`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a cocircuit contained in ``F``, if it exists; otherwise, an
        error is raised

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: sorted(sage.matroids.matroid.Matroid._cocircuit(M,
            ....:                             frozenset(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(sage.matroids.matroid.Matroid._cocircuit(M,
            ....:                                       frozenset(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no cocircuit in coindependent set.

        .. NOTE::

            This is an unguarded method. For the version that verifies if the
            input is indeed a subset of the groundset,
            see :meth:`<sage.matroids.matroid.Matroid.cocircuit>`.
        """
        self._pack(self._input, F)
        self.__cocircuit(self._output, self._input)
        return self.__unpack(self._output)

    cpdef frozenset _fundamental_cocircuit(self, frozenset B, e):
        r"""
        Return the `B`-fundamental circuit using `e`.

        Internal version that does no input checking.

        INPUT:

        - ``B`` -- a basis of the matroid
        - ``e`` -- an element of ``B``

        OUTPUT: set of elements

        EXAMPLES::

            sage: M = matroids.catalog.P8()
            sage: sorted(M._fundamental_cocircuit(frozenset('efgh'), 'e'))
            ['b', 'c', 'd', 'e']
        """
        self._pack(self._input, B)
        bitset_clear(self._input2)
        self.__move_current_basis(self._input, self._input2)
        self.__fundamental_cocircuit(self._output, self._idx[e])
        return self.__unpack(self._output)

    cpdef frozenset _coclosure(self, frozenset F):
        """
        Return the coclosure of a set.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: the smallest coclosed set containing ``X``

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: sorted(M._coclosure(frozenset(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']

        .. NOTE::

            This is an unguarded method. For the version that verifies if the
            input is indeed a subset of the groundset,
            see :meth:`<sage.matroids.matroid.Matroid.coclosure>`.
        """
        self._pack(self._input, F)
        self._coclosure_internal(self._output, self._input)
        return self.__unpack(self._output)

    cpdef frozenset _augment(self, frozenset X, frozenset Y):
        r"""
        Return a maximal subset `I` of `Y` such that `r(X + I)=r(X) + r(I)`.

        This version of ``augment`` does no type checking. In particular,
        ``Y`` is assumed to be disjoint from ``X``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``
        - ``Y`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``, and disjoint from ``X``

        OUTPUT: a subset `I` of ``Y`` such that `r(X + I)=r(X) + r(I)`.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: sorted(M._augment(frozenset(['a']), frozenset(['e', 'f', 'g', 'h'])))
            ['e', 'f', 'g']
        """
        self._pack(self._input, X)
        self._pack(self._input2, Y)
        self.__augment(self._output, self._input, self._input2)
        return self.__unpack(self._output)

    cpdef bint _is_independent(self, frozenset F) noexcept:
        """
        Test if input is independent.

        INPUT:

        - ``F`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: M._is_independent(frozenset(['a', 'b', 'c']))
            True
            sage: M._is_independent(frozenset(['a', 'b', 'c', 'd']))
            False

        .. NOTE::

            This is an unguarded method. For the version that verifies if
            the input is indeed a subset of the groundset,
            see :meth:`<sage.matroids.matroid.Matroid.is_independent>`.
        """
        self._pack(self._input, F)
        return self.__is_independent(self._input)

    # connectivity

    cpdef components(self):
        """
        Return an iterable containing the components of the matroid.

        A *component* is an inclusionwise maximal connected subset of the
        matroid. A subset is *connected* if the matroid resulting from
        deleting the complement of that subset is
        :meth:`connected <sage.matroids.matroid.Matroid.is_connected>`.

        OUTPUT: list of subsets

        .. SEEALSO::

            :meth:`M.is_connected() <sage.matroids.matroid.Matroid.is_connected>`,
            :meth:`M.delete() <sage.matroids.matroid.Matroid.delete>`

        EXAMPLES::

            sage: from sage.matroids.advanced import setprint
            sage: M = Matroid(ring=QQ, matrix=[[1, 0, 0, 1, 1, 0],
            ....:                              [0, 1, 0, 1, 2, 0],
            ....:                              [0, 0, 1, 0, 0, 1]])
            sage: setprint(M.components())
            [{0, 1, 3, 4}, {2, 5}]
        """
        cdef long i,j,e
        if not self._E:
            return SetSystem(self._E)
        cdef bitset_t *comp
        comp = <bitset_t*>sig_malloc((self.full_rank()) * sizeof(bitset_t))
        e = bitset_first(self._current_basis)
        i=0
        while e>=0:
            bitset_init(comp[i], self._bitset_size)
            self.__fundamental_cocircuit(comp[i], e)
            e=bitset_next(self._current_basis, e+1)
            i=i+1

        cdef bitset_t active_rows
        bitset_init(active_rows, <mp_bitcnt_t> self.full_rank()+1)
        bitset_set_first_n(active_rows, <mp_bitcnt_t> self.full_rank())
        i=0
        while i>=0:
            j = bitset_first(active_rows)
            while j>=0 and j<i:
                if not bitset_are_disjoint(comp[i], comp[j]):
                    bitset_union(comp[i], comp[i], comp[j])
                    bitset_discard(active_rows, j)
                j = bitset_next(active_rows, j+1)
            i = bitset_next(active_rows, i+1)

        res = SetSystem(self._E)
        i = bitset_first(active_rows)
        while i>=0:
            res._append(comp[i])
            i = bitset_next(active_rows, i+1)

        cdef bitset_t loop, loops
        bitset_init(loops, self._bitset_size)
        bitset_set_first_n(loops, len(self))
        i = bitset_first(active_rows)
        while i>=0:
            bitset_difference(loops, loops, comp[i])
            i = bitset_next(active_rows, i+1)

        bitset_init(loop, self._bitset_size)
        bitset_clear(loop)
        e = bitset_first(loops)
        while e>=0:
            bitset_add(loop, e)
            res._append(loop)
            bitset_discard(loop, e)
            e = bitset_next(loops, e+1)

        bitset_free(loops)
        bitset_free(loop)
        bitset_free(active_rows)
        for i in range(self.full_rank()):
            bitset_free(comp[i])
        sig_free(comp)
        return res

    cpdef _link(self, S, T):
        r"""
        Given disjoint subsets `S` and `T`, return a connector `I` and a separation `X`,
        which are optimal dual solutions in Tutte's Linking Theorem:

        .. MATH::

            \max \{ r_N(S) + r_N(T) - r(N) \mid N = M/I\setminus J, E(N) = S\cup T\}=\\
            \min \{ r_M(X) + r_M(Y) - r_M(E) \mid X \subseteq S, Y \subseteq T,
            E = X\cup Y, X\cap Y = \emptyset \}.

        Here `M` denotes this matroid.

        Internal version that does not verify that ``S`` and ``T``
        are sets, are disjoint, are subsets of the groundset.

        INPUT:

        - ``S`` -- a subset of the groundset
        - ``T`` -- a subset of the groundset disjoint from ``S``

        OUTPUT: tuple ``(I, X)`` of two frozensets

        ALGORITHM:

        Compute a maximum-cardinality common independent set `I` of
        of `M / S \setminus T` and `M \setminus S / T`.

        EXAMPLES::

            sage: M = matroids.catalog.BetsyRoss()
            sage: S = frozenset('ab')
            sage: T = frozenset('cd')
            sage: I, X = M._link(S, T)
            sage: M.connectivity(X)
            2
            sage: J = M.groundset() - (S | T | I)
            sage: N = (M / I).delete(J)
            sage: N.connectivity(S)
            2
        """
        # compute maximal common independent set of self\S/T and self/T\S
        cdef bitset_t SS, TT
        bitset_init(SS, self._groundset_size)
        bitset_init(TT, self._groundset_size)
        self._pack(SS,S)
        self._pack(TT,T)
        #F = set(self.groundset()) - (S | T)
        cdef bitset_t F, I
        bitset_init(F, self._groundset_size)
        bitset_init(I, self._groundset_size)
        bitset_union(self._input, SS, TT)
        bitset_complement(F, self._input)
        #I = self._augment(S|T, F)
        self.__augment(I, self._input, F)
        cdef bitset_t X, X1, X2, next_layer, todo, out_neighbors, R
        bitset_init(X, self._groundset_size)
        bitset_init(X1, self._groundset_size)
        bitset_init(X2, self._groundset_size)
        bitset_init(next_layer, self._groundset_size)
        bitset_init(todo, self._groundset_size)
        bitset_init(out_neighbors, self._groundset_size)
        bitset_init(R, self._groundset_size)
        cdef long* predecessor
        predecessor = <long*>sig_malloc(self._groundset_size*sizeof(long))
        cdef long e, u, y
        cdef bint found_path = True
        while found_path:
            #X = F - I
            bitset_difference(X,F,I)
            #X1 = X - self._closure(T|I)
            bitset_union(self._input, TT, I)
            self.__closure(X1, self._input)
            bitset_difference(X1,X,X1)
            #X2 = X - self._closure(S|I)
            bitset_union(self._input, SS, I)
            self.__closure(X2, self._input)
            bitset_difference(X2,X,X2)
            bitset_intersection(R, X1, X2)
            e = bitset_first(R)
            if e >= 0:
                bitset_add(I, e)
                continue
            #predecessor = {x: None for x in X1}
            e = bitset_first(X1)
            while e>=0:
                predecessor[e] = -1
                e = bitset_next(X1, e+1)
            #next_layer = set(X1)
            bitset_copy(next_layer, X1)
            bitset_union(R, SS, X1)
            found_path = False
            while not bitset_isempty(next_layer) and not found_path:
                #todo = next_layer
                bitset_copy(todo,next_layer)
                #next_layer = {}
                bitset_clear(next_layer)
                u = bitset_first(todo)
                while u>=0 and not found_path:
                    if bitset_in(X,u):
                        #out_neighbors = self._circuit(I|S.union([u])) - S.union([u])
                        bitset_union(self._input, I, SS)
                        bitset_add(self._input, u)
                        self.__circuit(out_neighbors, self._input)
                        bitset_discard(out_neighbors, u)
                    else:
                        #out_neighbors = X - self._closure(I|T - set([u]))
                        bitset_union(self._input, I, TT)
                        bitset_discard(self._input, u)
                        self.__closure(out_neighbors, self._input)
                        bitset_difference(out_neighbors, X, out_neighbors)
                    bitset_difference(out_neighbors, out_neighbors, R)
                    y = bitset_first(out_neighbors)
                    while y>=0:
                        predecessor[y] = u
                        if bitset_in(X2, y):
                            found_path = True
                            while y>=0:
                                bitset_flip(I, y)
                                y = predecessor[y]
                            break
                        bitset_add(R, y)
                        bitset_add(next_layer, y)
                        y = bitset_next(out_neighbors, y+1)
                    u = bitset_next(todo, u+1)
        II = self.__unpack(I)
        RR = self.__unpack(R)

        bitset_free(SS)
        bitset_free(TT)
        bitset_free(F)
        bitset_free(I)
        bitset_free(X)
        bitset_free(X1)
        bitset_free(X2)
        bitset_free(next_layer)
        bitset_free(todo)
        bitset_free(out_neighbors)
        bitset_free(R)
        sig_free(predecessor)

        return II, RR

    # enumeration

    cpdef list whitney_numbers2(self):
        r"""
        Return the Whitney numbers of the second kind of the matroid.

        The Whitney numbers of the second kind are here encoded as a vector
        `(W_0, \ldots, W_r)`, where `W_i` is the number of flats of rank `i`,
        and `r` is the rank of the matroid.

        OUTPUT: list of integers

        EXAMPLES::

            sage: M = matroids.catalog.S8()
            sage: M.whitney_numbers2()
            [1, 8, 22, 14, 1]
        """
        cdef bitset_t *flats
        cdef bitset_t *todo
        if self._matroid_rank == 0:
            return [1]
        flats = <bitset_t*>sig_malloc((self.full_rank() + 1) * sizeof(bitset_t))
        todo = <bitset_t*>sig_malloc((self.full_rank() + 1) * sizeof(bitset_t))

        for i in range(self.full_rank() + 1):
            bitset_init(flats[i], self._bitset_size)
            bitset_init(todo[i], self._bitset_size)
        f_vec = [0 for i in range(self.full_rank() + 1)]
        i = 0
        bitset_clear(todo[0])
        self.__closure(flats[0], todo[0])
        bitset_complement(todo[0], flats[0])
        self._whitney_numbers2_rec(f_vec, flats, todo, 0, 0)
        for i in range(self.full_rank() + 1):
            bitset_free(flats[i])
            bitset_free(todo[i])
        sig_free(flats)
        sig_free(todo)
        return f_vec

    cdef _whitney_numbers2_rec(self, object f_vec, bitset_t* flats, bitset_t* todo, long elt, long i):
        """
        Recursion for the ``whitney_numbers2`` method.
        """
        cdef long e
        f_vec[i] += 1
        e = bitset_next(todo[i], elt)
        while e >= 0:
            bitset_copy(self._input, flats[i])
            bitset_add(self._input, e)
            self.__closure(flats[i + 1], self._input)
            bitset_difference(todo[i], todo[i], flats[i + 1])
            bitset_difference(todo[i + 1], flats[i + 1], flats[i])
            if bitset_first(todo[i + 1]) == e:
                bitset_copy(todo[i + 1], todo[i])
                self._whitney_numbers2_rec(f_vec, flats, todo, e + 1, i + 1)
            e = bitset_next(todo[i], e)

    cpdef SetSystem flats(self, long k):
        """
        Return the collection of flats of the matroid of specified rank.

        A *flat* is a closed set.

        INPUT:

        - ``k`` -- integer

        OUTPUT: :class:`SetSystem`

        .. SEEALSO::

            :meth:`Matroid.closure() <sage.matroids.matroid.Matroid.closure>`

        EXAMPLES::

            sage: M = matroids.catalog.S8()
            sage: M.whitney_numbers2()
            [1, 8, 22, 14, 1]
            sage: len(M.flats(2))
            22
            sage: len(M.flats(8))
            0
            sage: len(M.flats(4))
            1
        """
        cdef bitset_t *flats
        cdef bitset_t *todo
        if k < 0 or k > self.full_rank():
            return SetSystem(self._E)
        if k == self.full_rank():
            return SetSystem(self._E, subsets=[self.groundset()])
        flats = <bitset_t*>sig_malloc((k + 1) * sizeof(bitset_t))
        todo = <bitset_t*>sig_malloc((k + 1) * sizeof(bitset_t))

        for i in range(k + 1):
            bitset_init(flats[i], self._bitset_size)
            bitset_init(todo[i], self._bitset_size)
        Rflats = SetSystem(self._E)
        i = 0
        bitset_clear(todo[0])
        self.__closure(flats[0], todo[0])
        bitset_complement(todo[0], flats[0])
        self._flats_rec(Rflats, k, flats, todo, 0, 0)
        for i in range(k + 1):
            bitset_free(flats[i])
            bitset_free(todo[i])
        sig_free(flats)
        sig_free(todo)
        return Rflats

    cdef _flats_rec(self, SetSystem Rflats, long R, bitset_t* flats, bitset_t* todo, long elt, long i):
        """
        Recursion for the ``flats`` method.
        """
        cdef long e
        if i == R:
            Rflats._append(flats[i])
            return
        e = bitset_next(todo[i], elt)
        while e >= 0:
            bitset_copy(self._input, flats[i])  # I fear that self._input is dangerous in a parallel computing environment. --SvZ
            bitset_add(self._input, e)          # It absolutely is! --RP
            self.__closure(flats[i + 1], self._input)
            bitset_difference(todo[i], todo[i], flats[i + 1])
            bitset_difference(todo[i + 1], flats[i + 1], flats[i])
            if bitset_first(todo[i + 1]) == e:
                bitset_copy(todo[i + 1], todo[i])
                self._flats_rec(Rflats, R, flats, todo, e + 1, i + 1)
            e = bitset_next(todo[i], e)

    cpdef SetSystem coflats(self, long k):
        """
        Return the collection of coflats of the matroid of specified corank.

        A *coflat* is a coclosed set.

        INPUT:

        - ``k`` -- integer

        OUTPUT: iterable containing all coflats of corank `k`

        .. SEEALSO::

            :meth:`Matroid.coclosure() <sage.matroids.matroid.Matroid.coclosure>`

        EXAMPLES::

            sage: M = matroids.catalog.S8().dual()
            sage: M.whitney_numbers2()
            [1, 8, 22, 14, 1]
            sage: len(M.coflats(2))
            22
            sage: len(M.coflats(8))
            0
            sage: len(M.coflats(4))
            1
        """
        cdef bitset_t *coflats
        cdef bitset_t *todo
        if k < 0 or k > self.full_corank():
            return SetSystem(self._E)
        if k == self.full_corank():
            return SetSystem(self._E, subsets=[self.groundset()])
        coflats = <bitset_t*>sig_malloc((k + 1) * sizeof(bitset_t))
        todo = <bitset_t*>sig_malloc((k + 1) * sizeof(bitset_t))

        for i in range(k + 1):
            bitset_init(coflats[i], self._bitset_size)
            bitset_init(todo[i], self._bitset_size)
        Rcoflats = SetSystem(self._E)
        i = 0
        bitset_clear(todo[0])
        self._coclosure_internal(coflats[0], todo[0])
        bitset_complement(todo[0], coflats[0])
        self._coflats_rec(Rcoflats, k, coflats, todo, 0, 0)
        for i in range(k + 1):
            bitset_free(coflats[i])
            bitset_free(todo[i])
        sig_free(coflats)
        sig_free(todo)
        return Rcoflats

    cdef _coflats_rec(self, SetSystem Rcoflats, long R, bitset_t* coflats, bitset_t* todo, long elt, long i):
        """
        Recursion for the ``coflats`` method.
        """
        cdef long e
        if i == R:
            Rcoflats._append(coflats[i])
            return
        e = bitset_next(todo[i], elt)
        while e >= 0:
            bitset_copy(self._input, coflats[i])
            bitset_add(self._input, e)
            self._coclosure_internal(coflats[i + 1], self._input)
            bitset_difference(todo[i], todo[i], coflats[i + 1])
            bitset_difference(todo[i + 1], coflats[i + 1], coflats[i])
            if bitset_first(todo[i + 1]) == e:
                bitset_copy(todo[i + 1], todo[i])
                self._coflats_rec(Rcoflats, R, coflats, todo, e + 1, i + 1)
            e = bitset_next(todo[i], e)

    cdef _flat_element_inv(self, long k):
        """
        Compute a flat-element invariant of the matroid.
        """
        cdef bitset_t *flats
        cdef bitset_t *todo
        if self._groundset_size == 0:
            return {}, tuple()
        flats = <bitset_t*>sig_malloc((k + 1) * sizeof(bitset_t))
        todo = <bitset_t*>sig_malloc((k + 1) * sizeof(bitset_t))

        for i in range(k + 1):
            bitset_init(flats[i], self._bitset_size)
            bitset_init(todo[i], self._bitset_size)
        f_inc = [[0 for e in range(self._groundset_size + 1)] for i in range(k + 1)]
        i = 0
        bitset_clear(todo[0])
        self.__closure(flats[0], todo[0])
        bitset_complement(todo[0], flats[0])
        self._flat_element_inv_rec(f_inc, k, flats, todo, 0, 0)
        for i in range(k + 1):
            bitset_free(flats[i])
            bitset_free(todo[i])
        sig_free(flats)
        sig_free(todo)
        fie = {}
        for e in range(self._groundset_size):
            t = tuple([f_inc[i][e] for i in range(k + 1)])
            if t in fie:
                fie[t].add(self._E[e])
            else:
                fie[t] = set([self._E[e]])
        f_vec = tuple([f_inc[i][self._groundset_size] for i in range(k + 1)])
        return fie, f_vec

    cdef _flat_element_inv_rec(self, object f_inc, long R, bitset_t* flats, bitset_t* todo, long elt, long i):
        """
        Recursion for ``_flat_element_inv``.
        """
        cdef long e, k
        k = bitset_len(flats[i])
        k = k * k
        e = bitset_first(flats[i])
        inc = f_inc[i]
        while e >= 0:
            inc[e] += k
            e = bitset_next(flats[i], e + 1)
        inc[self._groundset_size] += k
        if i == R:
            return
        e = bitset_next(todo[i], elt)
        while e >= 0:
            bitset_copy(self._input, flats[i])
            bitset_add(self._input, e)
            self.__closure(flats[i + 1], self._input)
            bitset_difference(todo[i], todo[i], flats[i + 1])
            bitset_difference(todo[i + 1], flats[i + 1], flats[i])
            if bitset_first(todo[i + 1]) == e:
                bitset_copy(todo[i + 1], todo[i])
                self._flat_element_inv_rec(f_inc, R, flats, todo, e + 1, i + 1)
            e = bitset_next(todo[i], e)

    cpdef bases_count(self):
        """
        Return the number of bases of the matroid.

        A *basis* is an inclusionwise maximal independent set.

        .. SEEALSO::

            :meth:`M.basis() <sage.matroids.matroid.Matroid.basis>`.

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.N1()
            sage: M.bases_count()
            184
        """
        if self._bcount is not None:
            return self._bcount
        cdef long res = 0
        bitset_clear(self._input)
        bitset_set_first_n(self._input, self._matroid_rank)
        repeat = True
        while repeat:
            if self.__is_independent(self._input):
                res += 1
            repeat = nxksrd(self._input, self._groundset_size, self._matroid_rank, True)
        self._bcount = res
        return self._bcount

    cpdef SetSystem independent_sets(self, long k=-1):
        r"""
        Return the independent sets of the matroid.

        INPUT:

        - ``k`` -- integer (optional); if specified, return the size-`k`
          independent sets of the matroid

        OUTPUT: :class:`SetSystem`

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: I = M.independent_sets()
            sage: len(I)
            57
            sage: M = matroids.catalog.N1()
            sage: M.bases_count()
            184
            sage: [len(M.independent_sets(k)) for k in range(M.full_rank() + 1)]
            [1, 10, 45, 120, 201, 184]

        TESTS::

            sage: len([B for B in M.bases()])
            184
        """
        cdef bitset_t *I
        cdef bitset_t *T
        cdef long i, e, r
        if k == -1:  # all independent sets
            res = SetSystem(self._E)
            bitset_clear(self._input)
            res._append(self._input)
            if not self._E:
                return res

            r = self.full_rank()
            I = <bitset_t*>sig_malloc((r + 1) * sizeof(bitset_t))
            T = <bitset_t*>sig_malloc((r + 1) * sizeof(bitset_t))
            for i in range(r + 1):
                bitset_init(I[i], self._bitset_size)
                bitset_init(T[i], self._bitset_size)

            i = 0
            bitset_clear(I[0])
            bitset_copy(self._input, I[0])
            self.__closure(T[0], self._input)
            while i >= 0:
                e = bitset_first_in_complement(T[i])
                if e >= 0:
                    bitset_add(T[i], e)
                    bitset_copy(I[i+1], I[i])
                    bitset_add(I[i+1], e)
                    res._append(I[i+1])
                    bitset_copy(self._input, I[i+1])
                    self.__closure(T[i+1], self._input)
                    bitset_union(T[i+1],T[i+1],T[i])
                    i = i + 1
                else:
                    i = i - 1
            for i in range(r + 1):
                bitset_free(I[i])
                bitset_free(T[i])
            sig_free(I)
            sig_free(T)
            return res

        # independent k-sets
        cdef SetSystem BB = SetSystem(self._E)
        if k < 0 or k > self.full_rank():
            return BB
        bitset_clear(self._input)
        bitset_set_first_n(self._input, k)
        repeat = True
        while repeat:
            if self.__is_independent(self._input):
                BB._append(self._input)
            repeat = nxksrd(self._input, self._groundset_size, k, True)
        return BB

    cpdef SetSystem dependent_sets(self, long k):
        """
        Return the dependent sets of fixed size.

        INPUT:

        - ``k`` -- integer

        OUTPUT: iterable containing all dependent sets of size ``k``

        EXAMPLES::

            sage: M = matroids.catalog.N1()
            sage: len(M.nonbases())
            68
            sage: [len(M.dependent_sets(k)) for k in range(M.full_rank() + 1)]
            [0, 0, 0, 0, 9, 68]

        TESTS::

            sage: binomial(M.size(), M.full_rank())-M.bases_count()                     # needs sage.symbolic
            68
            sage: len([B for B in M.nonbases()])
            68
        """
        cdef SetSystem NB
        NB = SetSystem(self._E)
        if k < 0 or k > self.size():
            return NB
        bitset_clear(self._input)
        bitset_set_first_n(self._input, k)
        repeat = True
        if k > self.full_rank():
            while repeat:
                NB._append(self._input)
                repeat = nxksrd(self._input, self._groundset_size, k, True)
        else:
            while repeat:
                if not self.__is_independent(self._input):
                    NB._append(self._input)
                repeat = nxksrd(self._input, self._groundset_size, k, True)
        NB.resize()
        return NB

    cpdef SetSystem nonspanning_circuits(self):
        """
        Return the nonspanning circuits of the matroid.

        A *nonspanning circuit* is a circuit whose rank is strictly smaller
        than the rank of the matroid.

        OUTPUT: iterable containing all nonspanning circuits

        .. SEEALSO::

            :meth:`Matroid.circuit() <sage.matroids.matroid.Matroid.circuit>`,
            :meth:`Matroid.rank() <sage.matroids.matroid.Matroid.rank>`

        EXAMPLES::

            sage: M = matroids.catalog.N1()
            sage: len(M.nonspanning_circuits())
            23
        """
        cdef SetSystem NSC
        NSC = SetSystem(self._E)
        if self._groundset_size == 0:
            return NSC
        bitset_clear(self._input)
        bitset_set_first_n(self._input, self._matroid_rank)
        cdef long e, f
        repeat = True
        while repeat:
            if self.__is_independent(self._input):
                bitset_complement(self._input2, self._current_basis)
                e = bitset_first(self._current_basis)
                while e >= 0:
                    self.__fundamental_cocircuit(self._output, e)
                    if e > bitset_first(self._output):
                        bitset_intersection(self._input2, self._input2, self._output)
                    e = bitset_next(self._current_basis, e + 1)
                f = bitset_first(self._input2)
                while f >= 0:
                    self.__fundamental_circuit(self._output, f)
                    if f == bitset_first(self._output) and bitset_len(self._output) <= self._matroid_rank:
                        NSC._append(self._output)
                    f = bitset_next(self._input2, f + 1)
            repeat = nxksrd(self._input, self._groundset_size, self._matroid_rank, True)
        NSC.resize()
        return NSC

    cpdef SetSystem noncospanning_cocircuits(self):
        """
        Return the noncospanning cocircuits of the matroid.

        A *noncospanning cocircuit* is a cocircuit whose corank is strictly
        smaller than the corank of the matroid.

        OUTPUT: iterable containing all nonspanning circuits

        .. SEEALSO::

            :meth:`Matroid.cocircuit() <sage.matroids.matroid.Matroid.cocircuit>`,
            :meth:`Matroid.corank() <sage.matroids.matroid.Matroid.corank>`

        EXAMPLES::

            sage: M = matroids.catalog.N1()
            sage: len(M.noncospanning_cocircuits())
            23
        """
        cdef SetSystem NSC
        NSC = SetSystem(self._E)
        if self._groundset_size == 0:
            return NSC
        bitset_clear(self._input)
        bitset_set_first_n(self._input, self._matroid_rank)
        cdef long e, f, corank
        corank = self._groundset_size - self._matroid_rank
        repeat = True
        while repeat:
            if self.__is_independent(self._input):
                bitset_copy(self._input2, self._current_basis)
                e = bitset_first(self._current_basis)
                for e in range(self._groundset_size):
                    if not bitset_in(self._current_basis, e):
                        self.__fundamental_circuit(self._output, e)
                        if e > bitset_first(self._output):
                            bitset_intersection(self._input2, self._input2, self._output)
                f = bitset_first(self._input2)
                while f >= 0:
                    self.__fundamental_cocircuit(self._output, f)
                    if f == bitset_first(self._output) and bitset_len(self._output) <= corank:
                        NSC._append(self._output)
                    f = bitset_next(self._input2, f + 1)
            repeat = nxksrd(self._input, self._groundset_size, self._matroid_rank, True)
        NSC.resize()
        return NSC

    cpdef SetSystem cocircuits(self):
        """
        Return the cocircuits of the matroid.

        OUTPUT: iterable containing all cocircuits

        .. SEEALSO::

            :meth:`Matroid.cocircuit() <sage.matroids.matroid.Matroid.cocircuit>`

        EXAMPLES::

            sage: M = Matroid(bases=matroids.catalog.NonFano().bases())
            sage: sorted([sorted(C) for C in M.cocircuits()])
            [['a', 'b', 'c', 'd', 'g'], ['a', 'b', 'c', 'e', 'g'],
             ['a', 'b', 'c', 'f', 'g'], ['a', 'b', 'd', 'e'],
             ['a', 'c', 'd', 'f'], ['a', 'e', 'f', 'g'], ['b', 'c', 'e', 'f'],
             ['b', 'd', 'f', 'g'], ['c', 'd', 'e', 'g']]
        """
        cdef SetSystem NSC
        NSC = SetSystem(self._E)
        if self._groundset_size == 0:
            return NSC
        bitset_clear(self._input)
        bitset_set_first_n(self._input, self._matroid_rank)
        cdef long e, f
        repeat = True
        while repeat:
            if self.__is_independent(self._input):
                bitset_copy(self._input2, self._current_basis)
                e = bitset_first(self._current_basis)
                for e in range(self._groundset_size):
                    if not bitset_in(self._current_basis, e):
                        self.__fundamental_circuit(self._output, e)
                        if e > bitset_first(self._output):
                            bitset_intersection(self._input2, self._input2, self._output)
                f = bitset_first(self._input2)
                while f >= 0:
                    self.__fundamental_cocircuit(self._output, f)
                    if f == bitset_first(self._output):
                        NSC._append(self._output)
                    f = bitset_next(self._input2, f + 1)
            repeat = nxksrd(self._input, self._groundset_size, self._matroid_rank, True)
        NSC.resize()
        return NSC

    cpdef SetSystem circuits(self, k=None):
        """
        Return the circuits of the matroid.

        OUTPUT: iterable containing all circuits

        .. SEEALSO::

            :meth:`M.circuit() <sage.matroids.matroid.Matroid.circuit>`

        EXAMPLES::

            sage: M = Matroid(matroids.catalog.NonFano().bases())
            sage: sorted([sorted(C) for C in M.circuits()])
            [['a', 'b', 'c', 'g'], ['a', 'b', 'd', 'e'], ['a', 'b', 'f'],
             ['a', 'c', 'd', 'f'], ['a', 'c', 'e'], ['a', 'd', 'e', 'f'],
             ['a', 'd', 'g'], ['a', 'e', 'f', 'g'], ['b', 'c', 'd'],
             ['b', 'c', 'e', 'f'], ['b', 'd', 'e', 'f'], ['b', 'd', 'f', 'g'],
             ['b', 'e', 'g'], ['c', 'd', 'e', 'f'], ['c', 'd', 'e', 'g'],
             ['c', 'f', 'g'], ['d', 'e', 'f', 'g']]
        """
        cdef SetSystem NSC
        NSC = SetSystem(self._E)
        if self._groundset_size == 0:
            return NSC
        bitset_clear(self._input)
        bitset_set_first_n(self._input, self._matroid_rank)
        cdef long e, f
        repeat = True
        while repeat:
            if self.__is_independent(self._input):
                bitset_complement(self._input2, self._current_basis)
                e = bitset_first(self._current_basis)
                while e >= 0:
                    self.__fundamental_cocircuit(self._output, e)
                    if e > bitset_first(self._output):
                        bitset_intersection(self._input2, self._input2, self._output)
                    e = bitset_next(self._current_basis, e + 1)
                f = bitset_first(self._input2)
                while f >= 0:
                    self.__fundamental_circuit(self._output, f)
                    if f == bitset_first(self._output):
                        NSC._append(self._output)
                    f = bitset_next(self._input2, f + 1)
            repeat = nxksrd(self._input, self._groundset_size, self._matroid_rank, True)
        NSC.resize()
        if k:
            return SetSystem(self.groundset(), [C for C in NSC if len(C) == k])
        else:
            return NSC

    # isomorphism

    cpdef _characteristic_setsystem(self):
        r"""
        Return a characteristic set-system for this matroid, on the same
        groundset.

        OUTPUT: :class:`<sage.matroids.set_system.SetSystem>`

        EXAMPLES::

            sage: M = matroids.catalog.N1()
            sage: M._characteristic_setsystem()
            SetSystem of 23 sets over 10 elements
            sage: len(M._characteristic_setsystem())
            23
        """
        if 2 * self._matroid_rank > self._groundset_size:
            return self.nonspanning_circuits()
        else:
            return self.noncospanning_cocircuits()

    cpdef _weak_invariant(self):
        """
        Return an isomorphism invariant of the matroid.

        Compared to ``BasisExchangeMatroid._strong_invariant()`` this invariant
        distinguishes less frequently between nonisomorphic matroids but takes
        less time to compute. See also
        :meth:`<BasisExchangeMatroid._weak_partition>`.

        OUTPUT: integer isomorphism invariant

        EXAMPLES::

            sage: M = Matroid(bases=matroids.catalog.Fano().bases())
            sage: N = Matroid(matroids.catalog.NonFano().bases())
            sage: M._weak_invariant() == N._weak_invariant()
            False
        """
        from sage.matroids.utilities import cmp_elements_key
        if self._weak_invariant_var is None:
            if self.full_rank() == 0 or self.full_corank() == 0:
                self._weak_invariant_var = 0
                self._weak_partition_var = SetSystem(self._E, [self.groundset()])
            else:
                k = min(self.full_rank() - 1, 2)
                fie, f_vec = self._flat_element_inv(k)
                self._weak_invariant_var = hash(tuple([tuple([(f, len(fie[f])) for f in sorted(fie, key=cmp_elements_key)]), f_vec]))
                self._weak_partition_var = SetSystem(self._E, [fie[f] for f in sorted(fie, key=cmp_elements_key)])
        return self._weak_invariant_var

    cpdef _weak_partition(self):
        """
        Return an ordered partition based on the incidences of elements with
        low-dimensional flats.

        EXAMPLES::

            sage: M = Matroid(matroids.catalog.Vamos().bases())
            sage: [sorted(p) for p in M._weak_partition()]
            [['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']]
        """
        self._weak_invariant()
        return self._weak_partition_var

    cpdef _strong_invariant(self):
        """
        Return an isomorphism invariant of the matroid.

        Compared to ``BasisExchangeMatroid._weak_invariant()`` this invariant
        distinguishes more frequently between nonisomorphic matroids but takes
        more time to compute. See also
        :meth:`<BasisExchangeMatroid._strong_partition>`.

        OUTPUT: integer isomorphism invariant

        EXAMPLES::

            sage: M = Matroid(matroids.catalog.Fano().bases())
            sage: N = Matroid(matroids.catalog.NonFano().bases())
            sage: M._strong_invariant() == N._strong_invariant()
            False
        """
        if self._strong_invariant_var is None:
            CP = self._characteristic_setsystem()._equitable_partition(self._weak_partition())
            self._strong_partition_var = CP[0]
            self._strong_invariant_var = CP[2]
        return self._strong_invariant_var

    cpdef _strong_partition(self):
        """
        Return an equitable partition which refines ``_weak_partition()``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: [sorted(p) for p in M._strong_partition()]
            [['a', 'b', 'e', 'f'], ['c', 'd', 'g', 'h']]
        """
        self._strong_invariant()
        return self._strong_partition_var

    cpdef _heuristic_invariant(self):
        """
        Return a number characteristic for the construction of
        ``_heuristic_partition()``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: N = BasisMatroid(matroids.catalog.Vamos())
            sage: M._heuristic_invariant() == N._heuristic_invariant()
            True
        """
        if self._heuristic_invariant_var is None:
            CP = self._characteristic_setsystem()._heuristic_partition(self._strong_partition())
            self._heuristic_partition_var = CP[0]
            self._heuristic_invariant_var = CP[2]
        return self._heuristic_invariant_var

    cpdef _heuristic_partition(self):
        """
        Return an ordered partition into singletons which refines an equitable
        partition of the matroid.

        The purpose of this partition is to heuristically find an isomorphism
        between two matroids, by lining up their respective
        ``heuristic_partitions``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: N = BasisMatroid(matroids.catalog.Vamos())
            sage: PM = M._heuristic_partition()
            sage: PN = N._heuristic_partition()
            sage: morphism = {}
            sage: for i in range(len(M)): morphism[min(PM[i])] = min(PN[i])
            sage: M._is_isomorphism(N, morphism)
            True
        """
        self._heuristic_invariant()
        return self._heuristic_partition_var

    cdef _flush(self):
        """
        Delete all invariants.
        """
        self._weak_invariant_var = None
        self._strong_invariant_var = None
        self._heuristic_invariant_var = None

    cpdef _equitable_partition(self, P=None):
        """
        Return the equitable refinement of a given ordered partition.

        The coarsest equitable partition of the groundset of this matroid
        that refines P.

        INPUT:

        - ``P`` -- ordered partition of the groundset (default: ``None``);
          if ``None``, the trivial partition is used

        OUTPUT: SetSystem

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Vamos())
            sage: [sorted(p) for p in M._equitable_partition()]
            [['a', 'b', 'e', 'f'], ['c', 'd', 'g', 'h']]
            sage: [sorted(p) for p in M._equitable_partition(['a', 'bcdefgh'])]
            [['a'], ['b'], ['e', 'f'], ['c', 'd', 'g', 'h']]
        """
        if P is not None:
            EQ = self._characteristic_setsystem()._equitable_partition(SetSystem(self._E, P))
        else:
            EQ = self._characteristic_setsystem()._equitable_partition()
        return EQ[0]

    cpdef _is_isomorphism(self, other, morphism):
        r"""
        Version of ``is_isomorphism()`` that does no type checking.

        INPUT:

        - ``other`` -- matroid
        - ``morphism`` -- dictionary mapping the groundset of ``self`` to
          the groundset of ``other``

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.catalog.Pappus()
            sage: N = BasisMatroid(matroids.catalog.NonPappus())
            sage: N._is_isomorphism(M, {e:e for e in M.groundset()})
            False

            sage: M = matroids.catalog.Fano().delete(['g'])
            sage: N = matroids.Wheel(3)
            sage: morphism = {'a':0, 'b':1, 'c': 2, 'd':4, 'e':5, 'f':3}
            sage: M._is_isomorphism(N, morphism)
            True

        TESTS:

        Check that :issue:`23300` was fixed::

            sage: def f(X):
            ....:     return min(len(X), 2)
            ....:
            sage: M = Matroid(groundset='abcd', rank_function=f)
            sage: N = Matroid(field=GF(3), reduced_matrix=[[1,1],[1,-1]])
            sage: N._is_isomorphism(M, {0:'a', 1:'b', 2:'c', 3:'d'})
            True
        """
        if not isinstance(other, BasisExchangeMatroid):
            from sage.matroids.basis_matroid import BasisMatroid
            ot = BasisMatroid(other)
        else:
            ot = other
        return self.__is_isomorphism(ot, morphism)

    cdef bint __is_isomorphism(self, BasisExchangeMatroid other, morphism) noexcept:
        """
        Bitpacked version of ``is_isomorphism``.
        """
        cdef long i
        morph = [other._idx[morphism[self._E[i]]] for i in range(len(self))]
        bitset_clear(self._input)
        bitset_set_first_n(self._input, self._matroid_rank)
        repeat = True
        while repeat:
            bitset_clear(other._input)
            i = bitset_first(self._input)
            while i != -1:
                bitset_add(other._input, <mp_bitcnt_t> morph[i])
                i = bitset_next(self._input, i + 1)
            if self.__is_independent(self._input) != other.__is_independent(other._input):
                return False
            repeat = nxksrd(self._input, self._groundset_size, self._matroid_rank, True)
        return True

    cpdef _isomorphism(self, other):
        """
        Return an isomorphism form ``self`` to ``other``, if one exists.

        Internal version that performs no checks on input.

        INPUT:

        - ``other`` -- matroid

        OUTPUT: dictionary or ``None``

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = matroids.Wheel(3)
            sage: M2 = matroids.CompleteGraphic(4)                                      # needs sage.graphs
            sage: morphism = M1._isomorphism(M2)                                        # needs sage.graphs
            sage: M1._is_isomorphism(M2, morphism)                                      # needs sage.graphs
            True
            sage: M1 = matroids.catalog.Fano()
            sage: M2 = matroids.catalog.NonFano()
            sage: M1._isomorphism(M2) is None
            True

        TESTS:

        Check that :issue:`23300` was fixed::

            sage: def f(X):
            ....:     return min(len(X), 2)
            ....:
            sage: M = Matroid(groundset='abcd', rank_function=f)
            sage: N = Matroid(field=GF(3), reduced_matrix=[[1,1],[1,-1]])
            sage: N._isomorphism(M) is not None
            True
        """
        if not isinstance(other, BasisExchangeMatroid):
            from sage.matroids.basis_matroid import BasisMatroid
            other = BasisMatroid(other)
        if self is other:
            return {e:e for e in self.groundset()}
        if len(self) != len(other):
            return None
        if self.full_rank() != other.full_rank():
            return None
        if self.full_rank() == 0 or self.full_corank() == 0:
            return {self.groundset_list()[i]: other.groundset_list()[i] for i in range(len(self))}

        if self._weak_invariant() != other._weak_invariant():
            return None
        PS = self._weak_partition()
        PO = other._weak_partition()
        if len(PS) == len(self) and len(PO) == len(other):
            morphism = {}
            for i in range(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            if self.__is_isomorphism(other, morphism):
                return morphism
            else:
                return None

        if self._strong_invariant() != other._strong_invariant():
            return False
        PS = self._strong_partition()
        PO = other._strong_partition()
        if len(PS) == len(self) and len(PO) == len(other):
            morphism = {}
            for i in range(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            if self.__is_isomorphism(other, morphism):
                return morphism
            else:
                return None

        if self._heuristic_invariant() == other._heuristic_invariant():
            PHS = self._heuristic_partition()
            PHO = other._heuristic_partition()
            morphism = {}
            for i in range(len(self)):
                morphism[min(PHS[i])] = min(PHO[i])
            if self.__is_isomorphism(other, morphism):
                return morphism

        return self._characteristic_setsystem()._isomorphism(other._characteristic_setsystem(), PS, PO)

    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Test if ``self`` is isomorphic to ``other``.

        Internal version that performs no checks on input.

        INPUT:

        - ``other`` -- matroid
        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, and, if ``certificate = True``, a dictionary giving
        the isomorphism or ``None``

        .. NOTE::

            Internal version that does no input checking.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = BasisMatroid(matroids.Wheel(3))
            sage: M2 = matroids.CompleteGraphic(4)                                      # needs sage.graphs
            sage: M1._is_isomorphic(M2)                                                 # needs sage.graphs
            True
            sage: M1._is_isomorphic(M2, certificate=True)                               # needs sage.graphs
            (True, {0: 0, 1: 1, 2: 2, 3: 3, 4: 5, 5: 4})
            sage: M1 = BasisMatroid(matroids.catalog.Fano())
            sage: M2 = matroids.catalog.NonFano()
            sage: M1._is_isomorphic(M2)
            False
            sage: M1._is_isomorphic(M2, certificate=True)
            (False, None)
        """
        if certificate:
            return self._is_isomorphic(other), self._isomorphism(other)
        if not isinstance(other, BasisExchangeMatroid):
            return other._is_isomorphic(self)
            # Either generic test, which converts other to BasisMatroid,
            # or overridden method.
        if self is other:
            return True
        if len(self) != len(other):
            return False
        if self.full_rank() != other.full_rank():
            return False
        if self.full_rank() == 0 or self.full_corank() == 0:
            return True
        if self.full_rank() == 1:
            return len(self.loops()) == len(other.loops())
        if self.full_corank() == 1:
            return len(self.coloops()) == len(other.coloops())

        if self._weak_invariant() != other._weak_invariant():
            return False
        PS = self._weak_partition()
        PO = other._weak_partition()
        if len(PS) == len(self) and len(PO) == len(other):
            morphism = {}
            for i in range(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            return self.__is_isomorphism(other, morphism)

        if self._strong_invariant() != other._strong_invariant():
            return False
        PS = self._strong_partition()
        PO = other._strong_partition()
        if len(PS) == len(self) and len(PO) == len(other):
            morphism = {}
            for i in range(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            return self.__is_isomorphism(other, morphism)

        if self._heuristic_invariant() == other._heuristic_invariant():
            PHS = self._heuristic_partition()
            PHO = other._heuristic_partition()
            morphism = {}
            for i in range(len(self)):
                morphism[min(PHS[i])] = min(PHO[i])
            if self.__is_isomorphism(other, morphism):
                return True

        return self._characteristic_setsystem()._isomorphism(other._characteristic_setsystem(), PS, PO) is not None

    cpdef is_valid(self, certificate=False):
        r"""
        Test if the data obey the matroid axioms.

        This method checks the basis axioms for the class. If `B` is the set
        of bases of a matroid `M`, then

        * `B` is nonempty
        * if `X` and `Y` are in `B`, and `x` is in `X - Y`, then there is a
          `y` in `Y - X` such that `(X - x) + y` is again a member of `B`.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, or (boolean, dictionary)

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BasisMatroid(matroids.catalog.Fano())
            sage: M.is_valid()
            True
            sage: M = Matroid(groundset='abcd', bases=['ab', 'cd'])
            sage: M.is_valid(certificate=True)
            (False, {'error': 'exchange axiom failed'})

        TESTS:

        Verify that :issue:`20172` was fixed::

            sage: M = Matroid(groundset='1234', bases=['12','13','23','34'])
            sage: M.is_valid()
            False
        """
        cdef long pointerX, pointerY, x, y, ln
        cdef bint foundpair
        cdef SetSystem BB
        BB = self.bases()
        pointerX = 0
        ln = len(BB)
        while pointerX < ln:  # for X in BB
            pointerY = 0
            while pointerY < ln:  # for Y in BB
                # Set current basis to Y
                bitset_difference(self._inside, self._current_basis, BB._subsets[pointerY])
                bitset_difference(self._outside, BB._subsets[pointerY], self._current_basis)
                self._move(self._inside, self._outside)
                if not bitset_eq(self._current_basis, BB._subsets[pointerY]):
                    # We failed to set the current basis to Y through basis exchanges.
                    # Therefore, the exchange axioms are violated!
                    return False if not certificate else (False, {"error": "exchange axiom failed"})
                bitset_difference(self._input, BB._subsets[pointerX], BB._subsets[pointerY])
                bitset_difference(self._input2, BB._subsets[pointerY], BB._subsets[pointerX])
                x = bitset_first(self._input)
                while x >= 0:  # for x in X-Y
                    foundpair = False
                    y = bitset_first(self._input2)
                    while y >= 0:  # for y in Y-X
                        if self._is_exchange_pair(y, x):
                            foundpair = True
                            y = -1
                        else:
                            y = bitset_next(self._input2, y + 1)
                    if not foundpair:
                        return False if not certificate else (False, {"error": "exchange axiom failed"})
                    x = bitset_next(self._input, x + 1)
                pointerY += 1
            pointerX += 1
        return True if not certificate else (True, {})

cdef bint nxksrd(bitset_s* b, long n, long k, bint succ) noexcept:
    """
    Next size-k subset of a size-n set in a revolving-door sequence.

    It will cycle through all such sets, returning each set exactly
    once. Each successive set differs from the last in exactly one
    element.

    This returns ``True`` if there is a next set, ``False`` otherwise.
    """
    # next k-subset of n-set in a revolving-door sequence
    if n == k or k == 0:
        return False
    if bitset_in(b, n - 1):
        if nxksrd(b, n - 1, k - 1, not succ):
            return True
        else:
            if succ:
                return False
            else:
                if k == 1:
                    bitset_add(b, n - 2)
                else:
                    bitset_add(b, k - 2)
                bitset_discard(b, n - 1)
                return True
    else:
        if nxksrd(b, n - 1, k, succ):
            return True
        else:
            if succ:
                if k == 1:
                    bitset_discard(b, n - 2)
                else:
                    bitset_discard(b, k - 2)
                bitset_add(b, n - 1)
                return True
            else:
                return False
