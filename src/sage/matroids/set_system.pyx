# sage_setup: distribution = sagemath-modules
"""
Set systems

Many matroid methods return a collection of subsets. In this module a class
:class:`SetSystem <sage.matroids.set_system.SetSystem>` is defined to do
just this. The class is intended for internal use, so all you can do as a user
is iterate over its members.

The class is equipped with partition refinement methods to help with matroid
isomorphism testing.

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

from cysignals.memory cimport check_allocarray, check_reallocarray, sig_free
from sage.data_structures.bitset_base cimport *

# SetSystem

cdef class SetSystem:
    """
    A ``SetSystem`` is an enumerator of a collection of subsets of a given
    fixed and finite groundset. It offers the possibility to enumerate its
    contents. One is most likely to encounter these as output from some
    Matroid methods::

        sage: M = matroids.catalog.Fano()
        sage: M.circuits()
        SetSystem of 14 sets over 7 elements

    To access the sets in this structure, simply iterate over them. The
    simplest way must be::

        sage: from sage.matroids.set_system import SetSystem
        sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
        sage: T = list(S)

    Or immediately use it to iterate::

        sage: from sage.matroids.set_system import SetSystem
        sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
        sage: [min(X) for X in S]
        [1, 3, 1]

    Note that this class is intended for runtime, so no loads/dumps mechanism
    was implemented.

    .. WARNING::

        The only guaranteed behavior of this class is that it is iterable. It
        is expected that M.circuits(), M.bases(), and so on will in the near
        future return actual iterators. All other methods (which are already
        hidden by default) are only for internal use by the Sage matroid code.
    """
    def __cinit__(self, groundset, subsets=None, capacity=1):
        """
        Init internal data structures.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: S
            SetSystem of 3 sets over 4 elements
        """
        cdef long i
        if not isinstance(groundset, tuple):
            self._groundset = tuple(groundset)
        else:
            self._groundset = groundset
        self._idx = {}
        for i in range(len(self._groundset)):
            self._idx[self._groundset[i]] = i

        self._groundset_size = len(groundset)
        self._bitset_size = max(self._groundset_size, 1)
        self._capacity = capacity
        self._subsets = <bitset_t*>check_allocarray(self._capacity, sizeof(bitset_t))
        bitset_init(self._temp, self._bitset_size)
        self._len = 0

    def __init__(self, groundset, subsets=None, capacity=1):
        """
        Create a SetSystem.

        INPUT:

        - ``groundset`` -- list or tuple of finitely many elements
        - ``subsets`` -- (default: ``None``) enumerator for a set of subsets of
          ``groundset``
        - ``capacity`` -- (default: ``1``) initial maximal capacity of the set
          system

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: S
            SetSystem of 3 sets over 4 elements
            sage: sorted(S[1])
            [3, 4]
            sage: for s in S: print(sorted(s))
            [1, 2]
            [3, 4]
            [1, 2, 4]
        """
        if subsets is not None:
            for e in subsets:
                self.append(e)

    def __dealloc__(self):
        cdef long i
        for i in range(self._len):
            bitset_free(self._subsets[i])
        sig_free(self._subsets)
        bitset_free(self._temp)

    def __len__(self):
        """
        Return the number of subsets in this SetSystem.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: S
            SetSystem of 3 sets over 4 elements
            sage: len(S)
            3
        """
        return self._len

    def __iter__(self):
        """
        Return an iterator for the subsets in this SetSystem.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: for s in S: print(sorted(s))
            [1, 2]
            [3, 4]
            [1, 2, 4]
        """
        return SetSystemIterator(self)

    def __getitem__(self, k):
        """
        Return the `k`-th subset in this SetSystem.

        INPUT:

        - ``k`` -- integer; the index of the subset in the system

        OUTPUT: the subset at index `k`

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: sorted(S[0])
            [1, 2]
            sage: sorted(S[1])
            [3, 4]
            sage: sorted(S[2])
            [1, 2, 4]
        """
        if k < len(self):
            return self.subset(k)
        else:
            raise ValueError("out of range")

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: repr(S)  # indirect doctest
            'SetSystem of 3 sets over 4 elements'
        """
        return f'SetSystem of {self._len} sets over {self._groundset_size} elements'

    cdef copy(self):
        cdef SetSystem S
        S = SetSystem(self._groundset, capacity=len(self))
        for i in range(len(self)):
            S._append(self._subsets[i])
        return S

    cdef _relabel(self, mapping):
        """
        Relabel each element `e` of the ground set as ``mapping[e]``, where
        ``mapping`` is a given injective map.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: none
        """
        cdef long i
        E = []
        for i in range(self._groundset_size):
            if self._groundset[i] in mapping:
                E.append(mapping[self._E[i]])
            else:
                E.append(self._E[i])
        self._groundset = E
        self._idx = {}
        for i in range(self._groundset_size):
            self._idx[self._groundset[i]] = i

    cpdef _complements(self):
        """
        Return a SetSystem containing the complements of each element in the
        groundset.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: T = S._complements()
            sage: for t in T: print(sorted(t))
            [3, 4]
            [1, 2]
            [3]
        """
        cdef SetSystem S
        if self._groundset_size == 0:
            return self
        S = SetSystem(self._groundset, capacity=len(self))
        for i in range(len(self)):
            bitset_complement(self._temp, self._subsets[i])
            S._append(self._temp)
        return S

    cdef inline resize(self, k=None):
        """
        Change the capacity of the SetSystem.
        """
        if k is None:
            k = self._len
        for i in range(k, self._len):
            bitset_free(self._subsets[i])
        self._len = min(self._len, k)
        k2 = max(k, 1)
        self._subsets = <bitset_t*>check_reallocarray(self._subsets, k2, sizeof(bitset_t))
        self._capacity = k2

    cdef inline _append(self, bitset_t X):
        """
        Append subset in internal, bitset format.
        """
        if self._capacity == self._len:
            self.resize(self._capacity * 2)
        bitset_init(self._subsets[self._len], self._bitset_size)
        bitset_copy(self._subsets[self._len], X)
        self._len += 1

    cdef inline append(self, X):
        """
        Append subset.
        """
        if self._capacity == self._len:
            self.resize(self._capacity * 2)
        bitset_init(self._subsets[self._len], self._bitset_size)
        bitset_clear(self._subsets[self._len])
        for x in X:
            bitset_add(self._subsets[self._len], <mp_bitcnt_t> self._idx[x])
        self._len += 1

    cdef inline _subset(self, long k):
        """
        Return the `k`-th subset, in index format.
        """
        return bitset_list(self._subsets[k])

    cdef subset(self, k):
        """
        Return the `k`-th subset.
        """
        cdef long i
        F = set()
        i = bitset_first(self._subsets[k])
        while i >= 0:
            F.add(self._groundset[i])
            i = bitset_next(self._subsets[k], i + 1)
        return frozenset(F)

    cpdef _get_groundset(self):
        """
        Return the groundset of this SetSystem.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: sorted(S._get_groundset())
            [1, 2, 3, 4]
        """
        return frozenset(self._groundset)

    cpdef is_connected(self):
        """
        Test if the :class:`SetSystem` is connected.

        A :class:`SetSystem` is connected if there is no nonempty proper subset
        ``X`` of the groundset so the each subset is either contained in ``X``
        or disjoint from ``X``.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: S.is_connected()
            True
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4]])
            sage: S.is_connected()
            False
            sage: S = SetSystem([1], [])
            sage: S.is_connected()
            True
        """
        if self._groundset_size <= 1:
            return True
        cdef long i
        bitset_clear(self._temp)
        cdef bitset_t active
        bitset_init(active, self._len)
        bitset_complement(active, active)

        # We compute the union of all sets containing 0, and deactivate them.
        for i in range(self._len):
            if bitset_in(self._subsets[i], 0):
                bitset_union(self._temp, self._subsets[i], self._temp)
                bitset_discard(active, i)

        cdef bint closed = False
        while not closed:
            closed = True

            # We update _temp with all active sets that intersects it. If there
            # is no such set, then _temp is closed (i.e. a connected component).
            i = bitset_first(active)
            while i>=0:
                if not bitset_are_disjoint(self._temp, self._subsets[i]):
                    bitset_union(self._temp, self._subsets[i], self._temp)
                    bitset_discard(active, i)
                    closed = False
                i = bitset_next(active, i+1)
        bitset_free(active)
        bitset_complement(self._temp, self._temp)
        return bitset_isempty(self._temp)

    # isomorphism

    cdef list _incidence_count(self, E):
        """
        For the sub-collection indexed by ``E``, count how often each element
        occurs.
        """
        cdef long i, e
        cdef list cnt
        cnt = [0 for v in range(self._groundset_size)]
        for e in E:
            i = bitset_first(self._subsets[e])
            while i >= 0:
                cnt[i] += 1
                i = bitset_next(self._subsets[e], i + 1)
        return cnt

    cdef SetSystem _groundset_partition(self, SetSystem P, list cnt):
        """
        Helper method for partition methods below.
        """
        cdef dict C
        cdef long i, v, t0, t
        cdef bint split

        C = {}
        for i in range(len(P)):
            v = bitset_first(P._subsets[i])
            if v < 0:
                continue
            t0 = cnt[v]
            v = bitset_next(P._subsets[i], v + 1)
            split = False
            while v >= 0:
                t = cnt[v]
                if t != t0:
                    split = True
                    if t < t0:
                        t0 = t
                v = bitset_next(P._subsets[i], v + 1)
            if split:
                C.clear()
                v = bitset_first(P._subsets[i])
                while v >= 0:
                    t = cnt[v]
                    if t != t0:
                        if t in C:
                            C[t].add(v)
                        else:
                            C[t] = set([v])
                    v = bitset_next(P._subsets[i], v + 1)
                for t in sorted(C):
                    bitset_clear(self._temp)
                    for v in C[t]:
                        bitset_add(self._temp, v)
                        bitset_discard(P._subsets[i], v)
                    P._append(self._temp)

    cdef long subset_characteristic(self, SetSystem P, long e) noexcept:
        """
        Helper method for partition methods below.
        """
        cdef long c
        c = 0
        for p in range(len(P)):
            c <<= bitset_len(P._subsets[p])
            bitset_intersection(self._temp, P._subsets[p], self._subsets[e])
            c += bitset_len(self._temp)
        return c

    cdef subsets_partition(self, SetSystem P=None, E=None):
        """
        Helper method for partition methods below.
        """
        if P is None:
            P = self.groundset_partition()
        if E is None:
            E = range(self._len)
        if len(E) == 0:
            return [E]

        ED = [(self.subset_characteristic(P, e), e) for e in E]
        ED.sort()

        EP = []
        ep = []
        d = ED[0][0]
        eh = [d]
        for ed in ED:
            if ed[0] != d:
                EP.append(ep)
                ep = [ed[1]]
                d = ed[0]
            else:
                ep.append(ed[1])
            eh.append(ed[0])
        EP.append(ep)
        return EP, hash(tuple(eh))

    cdef _distinguish(self, Py_ssize_t v):
        """
        Helper method for partition methods below.
        """
        cdef SetSystem S
        S = SetSystem(self._groundset, capacity=len(self) + 1)
        bitset_clear(self._temp)
        bitset_add(self._temp, v)
        for i in range(len(self)):
            bitset_difference(S._temp, self._subsets[i], self._temp)
            S._append(S._temp)
        S._append(self._temp)
        return S

    # partition functions
    cdef initial_partition(self, SetSystem P=None, E=None):
        """
        Helper method for partition methods below.
        """
        if E is None:
            E = range(self._len)
        if P is None:
            if self._groundset:
                P = SetSystem(self._groundset, [self._groundset], capacity=self._groundset_size)
            else:
                P = SetSystem([], [])
        cnt = self._incidence_count(E)
        self._groundset_partition(P, cnt)
        return P

    cpdef _equitable_partition(self, SetSystem P=None, EP=None):
        r"""
        Return an equitable ordered partition of the groundset of the
        hypergraph whose edges are the subsets in this SetSystem.

        Given any ordered partition `P = (p_1, ..., p_k)` of the groundset of
        a hypergraph, any edge `e` of the hypergraph has a characteristic
        intersection number sequence `i(e)=(|p_1\cap e|, ... , |p_k\cap e|))`.
        There is an ordered partition `EP` of the edges that groups the edges
        according to this intersection number sequence. Given this an ordered
        partition of the edges, we may similarly refine `P` to a new ordered
        partition `P'`, by considering the incidence numbers of groundset
        elements with each partition element of `EP`.

        The ordered partition `P` is equitable when `P' = P`.

        INPUT:

        - ``P`` -- an equitable ordered partition of the groundset, stored as
          a SetSystem
        - ``EP`` -- the corresponding equitable partition of the edges, stored
          as a list of lists of indices of subsets of this SetSystem

        OUTPUT:

        - ``P`` -- an equitable ordered partition of the groundset, stored as a
          SetSystem
        - ``EP`` -- the corresponding equitable partition of the edges, stored
          as a list of lists of indices of subsets of this SetSystem
        - ``h`` -- integer invariant of the SetSystem

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: for p in S._equitable_partition()[0]: print(sorted(p))
            [3]
            [4]
            [1, 2]
            sage: T = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 3, 4]])
            sage: for p in T._equitable_partition()[0]: print(sorted(p))
            [2]
            [1]
            [3, 4]

        .. NOTE::

            We do not maintain any well-defined order when refining a
            partition. We do maintain that the resulting order of the
            partition elements is an invariant of the isomorphism class of the
            hypergraph.
        """
        cdef long h
        cdef list EP2, H

        if P is None:
            P = self.initial_partition()
        else:
            P = P.copy()
        if EP is None:
            EP = self.subsets_partition(P)[0]

        h = len(EP)
        pl = 0
        while len(P) > pl:
            H = [h]
            pl = len(P)
            EP2 = []
            for ep in EP:
                SP, h = self.subsets_partition(P, ep)
                H.append(h)
                for p in SP:
                    cnt = self._incidence_count(p)
                    self._groundset_partition(P, cnt)
                    if len(p) > 1:
                        EP2.append(p)
            EP = EP2
            h = hash(tuple(H))

        return P, EP, h

    cpdef _heuristic_partition(self, SetSystem P=None, EP=None):
        """
        Return a heuristic ordered partition into singletons of the ground
        set of the hypergraph whose edges are the subsets in this SetSystem.

        This partition obtained as follows: make an equitable
        partition ``P``, and while ``P`` has a partition element ``p`` with
        more than one element, select an arbitrary ``e`` from the first such
        ``p`` and split ``p`` into ``p-e``. Then replace ``P`` with
        the equitable refinement of this partition.

        INPUT:

        - ``P`` -- (default: ``None``) an ordered partition of the groundset
        - ``EP`` -- (default: ``None``) the corresponding partition of the
          edges, stored as a list of lists of indices of subsets of this
          SetSystem

        OUTPUT:

        - ``P`` -- an ordered partition of the groundset into singletons,
          stored as a SetSystem
        - ``EP`` -- the corresponding partition of the edges, stored as a list
          of lists of indices of subsets of this SetSystem
        - ``h`` -- integer invariant of the SetSystem

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: for p in S._heuristic_partition()[0]: print(sorted(p))
            [3]
            [4]
            [2]
            [1]
            sage: T = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 3, 4]])
            sage: for p in T._heuristic_partition()[0]: print(sorted(p))
            [2]
            [1]
            [4]
            [3]
        """
        P, EP, h = self._equitable_partition(P, EP)
        for i in range(len(P)):
            if bitset_len(P._subsets[i]) > 1:
                return self._heuristic_partition(P._distinguish(bitset_first(P._subsets[i])), EP)
        return P, EP, h

    cpdef _isomorphism(self, SetSystem other, SetSystem SP=None, SetSystem OP=None):
        """
        Return a groundset isomorphism between this SetSystem and an other.

        INPUT:

        - ``other`` -- SetSystem
        - ``SP`` -- (optional) SetSystem storing an ordered partition of the
          groundset of ``self``
        - ``OP`` -- (optional) SetSystem storing an ordered partition of the
          groundset of ``other``

        OUTPUT:

        ``morphism`` -- dictionary containing an isomorphism respecting the
        given ordered partitions, or ``None`` if no such isomorphism exists.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: T = SetSystem(['a', 'b', 'c', 'd'], [['a', 'b'], ['c', 'd'],
            ....:                                      ['a', 'c', 'd']])
            sage: S._isomorphism(T)
            {1: 'c', 2: 'd', 3: 'b', 4: 'a'}
            sage: S = SetSystem([], [])
            sage: S._isomorphism(S)
            {}
        """
        cdef long l, p, v
        if SP is None or OP is None:
            SP, SEP, sh = self._equitable_partition()
            OP, OEP, oh = other._equitable_partition()
            if sh != oh:
                return None
        if len(SP) != len(OP):
            return None
        p = SP._groundset_size + 1
        for i in range(len(SP)):
            l = bitset_len(SP._subsets[i])
            if l != bitset_len(OP._subsets[i]):
                return None
            if l != 1 and l < p:
                p = l
        for i in range(len(SP)):
            if bitset_len(SP._subsets[i]) == p:
                SP2, SEP, sh = self._equitable_partition(SP._distinguish(bitset_first(SP._subsets[i])))
                v = bitset_first(OP._subsets[i])
                while v >= 0:
                    OP2, OEP, oh = other._equitable_partition(OP._distinguish(v))
                    if sh == oh:
                        m = self._isomorphism(other, SP2, OP2)
                        if m is not None:
                            return m
                    v = bitset_next(OP._subsets[i], v + 1)
                return None
        if sorted([self.subset_characteristic(SP, i) for i in range(len(self))]) != sorted([other.subset_characteristic(OP, i) for i in range(len(other))]):
            return None
        return dict([(self._groundset[bitset_first(SP._subsets[i])], other._groundset[bitset_first(OP._subsets[i])]) for i in range(len(SP))])

    cpdef _equivalence(self, is_equiv, SetSystem other, SetSystem SP=None, SetSystem OP=None):
        """
        Return a groundset isomorphism that is an equivalence between this
        SetSystem and an other.

        INPUT:

        - ``is_equiv`` -- a function that determines if a given groundset
          isomorphism is a valid equivalence
        - ``other`` -- SetSystem
        - ``SP`` -- (optional) SetSystem storing an ordered partition of the
          groundset of ``self``
        - ``OP`` -- (optional) SetSystem storing an ordered partition of the
          groundset of ``other``

        OUTPUT:

        ``morphism``, a dictionary containing an isomorphism respecting the
        given ordered partitions, so that ``is_equiv(self, other, morphism)``
        is ``True``; or ``None`` if no such equivalence exists.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: T = SetSystem(['a', 'b', 'c', 'd'], [['a', 'b'], ['c', 'd'],
            ....:                                      ['a', 'c', 'd']])
            sage: S._equivalence(lambda self, other, morph:True, T)
            {1: 'c', 2: 'd', 3: 'b', 4: 'a'}

        Check that :issue:`15189` is fixed::

            sage: M = Matroid(ring=GF(5), reduced_matrix=[[1,0,3],[0,1,1],[1,1,0]])
            sage: N = Matroid(ring=GF(5), reduced_matrix=[[1,0,1],[0,1,1],[1,1,0]])
            sage: M.is_field_isomorphic(N)
            False
            sage: any(M.is_field_isomorphism(N, p) for p in Permutations(range(6)))
            False
        """
        cdef long v
        if SP is None or OP is None:
            SP, SEP, sh = self._equitable_partition()
            OP, OEP, oh = other._equitable_partition()
            if sh != oh:
                return None
        if len(SP) != len(OP):
            return None
        for i in range(len(SP)):
            if bitset_len(SP._subsets[i]) != bitset_len(OP._subsets[i]):
                return None
        for i in range(len(SP)):
            if bitset_len(SP._subsets[i]) > 1:
                SP2, SEP, sh = self._equitable_partition(SP._distinguish(bitset_first(SP._subsets[i])))
                v = bitset_first(OP._subsets[i])
                while v >= 0:
                    OP2, OEP, oh = other._equitable_partition(OP._distinguish(v))
                    if sh == oh:
                        m = self._equivalence(is_equiv, other, SP2, OP2)
                        if m is not None:
                            return m
                    v = bitset_next(OP._subsets[i], v + 1)
                return None
        morph = dict([(self._groundset[bitset_first(SP._subsets[i])], other._groundset[bitset_first(OP._subsets[i])]) for i in range(len(SP))])
        if is_equiv(self, other, morph):
            return morph


cdef class SetSystemIterator:
    def __init__(self, H):
        """
        Create an iterator for a SetSystem.

        Called internally when iterating over the contents of a SetSystem.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: type(S.__iter__())
            <... 'sage.matroids.set_system.SetSystemIterator'>
        """
        self._H = H
        self._pointer = -1
        self._len = len(H)

    def __next__(self):
        """
        Return the next subset of a SetSystem.

        EXAMPLES::

            sage: from sage.matroids.set_system import SetSystem
            sage: S = SetSystem([1, 2, 3, 4], [[1, 2], [3, 4], [1, 2, 4]])
            sage: I = S.__iter__()
            sage: sorted(I.__next__())
            [1, 2]
            sage: sorted(I.__next__())
            [3, 4]
            sage: sorted(I.__next__())
            [1, 2, 4]
        """
        self._pointer += 1
        if self._pointer == self._len:
            self._pointer = -1
            raise StopIteration
        else:
            return self._H.subset(self._pointer)
