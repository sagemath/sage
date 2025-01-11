# sage.doctest: needs sage.rings.finite_rings
r"""
Evenly distributed sets in finite fields

This module consists of a simple class :class:`EvenlyDistributedSetsBacktracker`. Its main
purpose is to iterate through the evenly distributed sets of a given finite
field.

The naive backtracker implemented here is not directly used to generate
difference family as even for small parameters it already takes time to run.
Instead, its output has been stored into a database
:mod:`sage.combinat.designs.database`. If the backtracker is improved, then one
might want to update this database with more values.

Classes and methods
-------------------
"""

cimport cython

from libc.limits cimport UINT_MAX
from libc.string cimport memset, memcpy
from memory_allocator cimport MemoryAllocator

from sage.rings.integer cimport smallInteger

from sage.categories.fields import Fields


cdef class EvenlyDistributedSetsBacktracker:
    r"""
    Set of evenly distributed subsets in finite fields.

        **Definition:** Let `K` be a finite field of cardinality `q` and `k` an
        integer so that `k(k-1)` divides `q-1`. Let `H = K^*` be the
        multiplicative group of invertible elements in `K`. A `k`-*evenly
        distributed set* in `K` is a set `B = \{b_1, b_2, \ldots, b_k\}` of `k`
        elements of `K` so that the `k(k-1)` differences `\Delta B = \{b_i -
        b_j; i \not= j\}` hit each coset modulo `H^{2(q-1)/(k(k-1))}` exactly
        twice.

    Evenly distributed sets were introduced by Wilson [Wi72]_ (see also
    [BJL99-1]_, Chapter VII.6). He proved that for any `k`, and for any prime power
    `q` large enough such that `k(k-1)` divides `k` there exists a `k`-evenly
    distributed set in the field of cardinality `q`. This existence result based
    on a counting argument (using Dirichlet theorem) does not provide a simple
    method to generate them.

    From a `k`-evenly distributed set, it is straightforward to build a
    `(q,k,1)`-difference family (see :meth:`to_difference_family`). Another
    approach to generate a difference family, somehow dual to this one, is via
    radical difference family (see in particular
    :func:`~sage.combinat.designs.difference_family.radical_difference_family`
    from the module :mod:`~sage.combinat.designs.difference_family`).

    By default, this backtracker only considers evenly distributed sets up to
    affine automorphisms, i.e. `B` is considered equivalent to `s B + t` for any
    invertible element `s` and any element `t` in the field `K`. Note that the
    set of differences is just multiplicatively translated by `s` as `\Delta (s
    B + t) = s (\Delta B)`, and so that `B` is an evenly distributed set if and
    only if `sB` is one too. This behaviour can be modified via the argument
    ``up_to_isomorphism`` (see the input section and the examples below).

    INPUT:

    - ``K`` -- a finite field of cardinality `q`

    - ``k`` -- positive integer such that `k(k-1)` divides `q-1`

    - ``up_to_isomorphism`` -- boolean (default: ``True``); whether only consider
      evenly distributed sets up to automorphisms of the field of the form
      `x \mapsto ax + b`. If set to ``False`` then the iteration is over all
      evenly distributed sets that contain ``0`` and ``1``.

    - ``check`` -- boolean (default: ``False``); whether you want to check
      intermediate steps of the iterator. This is mainly intended for debugging
      purpose. Setting it to ``True`` will considerably slow the iteration.

    EXAMPLES:

    The main part of the code is contained in the iterator. To get one evenly
    distributed set just do::

        sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker
        sage: E = EvenlyDistributedSetsBacktracker(Zmod(151),6)
        sage: B = E.an_element()
        sage: B
        [0, 1, 69, 36, 57, 89]

    The class has a method to convert it to a difference family::

        sage: E.to_difference_family(B)
        [[0, 1, 69, 36, 57, 89],
         [0, 132, 48, 71, 125, 121],
         [0, 59, 145, 10, 41, 117],
         [0, 87, 114, 112, 127, 42],
         [0, 8, 99, 137, 3, 108]]

    It is also possible to run over all evenly distributed sets::

        sage: E = EvenlyDistributedSetsBacktracker(Zmod(13), 4, up_to_isomorphism=False)
        sage: for B in E: print(B)
        [0, 1, 11, 5]
        [0, 1, 4, 6]
        [0, 1, 9, 3]
        [0, 1, 8, 10]

        sage: E = EvenlyDistributedSetsBacktracker(Zmod(13), 4, up_to_isomorphism=True)
        sage: for B in E: print(B)
        [0, 1, 11, 5]


    Or only count them::

        sage: for k in range(13, 200, 12):
        ....:     if is_prime_power(k):
        ....:         K = GF(k,'a')
        ....:         E1 = EvenlyDistributedSetsBacktracker(K, 4, False)
        ....:         E2 = EvenlyDistributedSetsBacktracker(K, 4, True)
        ....:         print("{:3} {:3} {:3}".format(k, E1.cardinality(), E2.cardinality()))
         13   4   1
         25  40   4
         37  12   1
         49  24   2
         61  12   1
         73  48   4
         97  64   6
        109  72   6
        121 240  20
        157  96   8
        169 240  20
        181 204  17
        193 336  28

    Note that by definition, the number of evenly distributed sets up to
    isomorphisms is at most `k(k-1)` times smaller than without isomorphisms.
    But it might not be exactly `k(k-1)` as some of them might have symmetries::

        sage: B = EvenlyDistributedSetsBacktracker(Zmod(13), 4).an_element()
        sage: B
        [0, 1, 11, 5]
        sage: [9*x + 5 for x in B]
        [5, 1, 0, 11]
        sage: [3*x + 11 for x in B]
        [11, 1, 5, 0]
    """
    # PYTHON DATA
    cdef K                       # the underlying field
    cdef list list_K             # the elements of K  (i -> x)
    cdef dict K_to_int           # inverse of list_K  (x -> i)

    # FLAGS
    cdef int count               # do we count or do we iterate
    cdef int check               # do we need to check (debug)
    cdef int up_to_isom          # do we care only about isomorphisms

    # STATIC DATA
    cdef unsigned int q          # cardinality of the field
    cdef unsigned int k          # size of the subsets
    cdef unsigned int e          # k(k-1)/2
    cdef unsigned int m          # (q-1) / e
    cdef unsigned int ** diff    # qxq array: diff[x][y]  = x - y
    cdef unsigned int ** ratio   # qxq array: ratio[x][y] = x / y
    cdef unsigned int * min_orb  # q array  : min_orb[x]  = min {x, 1-x, 1/x,
                                 #                   1/(1-x), (x-1)/x, x/(x-1)}

    # DYNAMIC DATA
    cdef unsigned int * B        # current stack of elements of {0,...,q-1}
    cdef unsigned int * cosets   # e array: cosets of differences of elts in B
    cdef unsigned int * t        # e array: temporary variable for updates

    # MANAGEMENT OF MEMORY
    cdef MemoryAllocator mem

    def __dealloc__(self):
        pass

    def __init__(self, K, k, up_to_isomorphism=True, check=False):
        r"""
        TESTS::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: EvenlyDistributedSetsBacktracker(Zmod(4),2)
            Traceback (most recent call last):
            ...
            ValueError: Ring of integers modulo 4 is not a field

            sage: EvenlyDistributedSetsBacktracker(Zmod(71),7)
            Traceback (most recent call last):
            ...
            ValueError: k(k-1)=42 does not divide q-1=70

        For `q=421` which is congruent to 1 modulo `12`, `20`, `30` and `42` we
        run backtracker with the ``check`` argument set to ``True``::

            sage: for _ in EvenlyDistributedSetsBacktracker(Zmod(421), 4, check=True):
            ....:     pass
            sage: for _ in EvenlyDistributedSetsBacktracker(Zmod(421), 5, check=True):
            ....:    pass
            sage: for _ in EvenlyDistributedSetsBacktracker(Zmod(421), 6, check=True):
            ....:    pass
            sage: for _ in EvenlyDistributedSetsBacktracker(Zmod(421), 7, check=True):
            ....:    pass
        """
        self.check      = bool(check)
        self.up_to_isom = bool(up_to_isomorphism)
        self.count      = 0

        cdef unsigned int i,j

        if K not in Fields():
            raise ValueError(f"{K} is not a field")
        cdef unsigned int q = K.cardinality()
        cdef unsigned int e = k*(k-1)/2
        if (q-1) % (2*e):
            raise ValueError("k(k-1)={} does not divide q-1={}".format(k*(k-1),q-1))
        cdef unsigned int m = (q - 1) // e

        self.q = q
        self.e = e
        self.k = k
        self.m = (q - 1) // e
        self.K = K

        self.mem = MemoryAllocator()
        self.diff = <unsigned int **> self.mem.calloc(q, sizeof(unsigned int *))
        self.diff[0] = <unsigned int *> self.mem.malloc(q*q*sizeof(unsigned int))
        for i in range(1, self.q):
            self.diff[i] = self.diff[i-1] + q

        self.ratio = <unsigned int **> self.mem.calloc(q, sizeof(unsigned int *))
        self.ratio[0] = <unsigned int *> self.mem.malloc(q*q*sizeof(unsigned int))
        for i in range(1, self.q):
            self.ratio[i] = self.ratio[i-1] + q

        self.B = <unsigned int *> self.mem.malloc(k*sizeof(unsigned int))
        self.min_orb = <unsigned int *> self.mem.malloc(q*sizeof(unsigned int))
        self.cosets = <unsigned int *> self.mem.malloc(e*sizeof(unsigned int))
        self.t = <unsigned int *> self.mem.malloc(e*sizeof(unsigned int))

        x = K.multiplicative_generator()
        list_K = []
        for i in range(e):
            list_K.extend(sorted(x**(j*e+i) for j in range(m)))
        list_K.append(K.zero())
        self.list_K = list_K
        K_to_int = self.K_to_int = {y:i for i,y in enumerate(list_K)}

        zero = K.zero()
        one = K.one()
        assert self.K_to_int[zero] == q-1
        assert self.K_to_int[one] == 0
        assert set(K) == set(list_K)

        self.min_orb[0] = self.min_orb[q-1] = 0
        for i,x in enumerate(self.list_K):
            if x != zero and x != one:
                self.min_orb[i] = min(K_to_int[z] for z in
                        [x, one/x, one-x, one/(one-x), (x-one)/x, x/(x-one)])
            for j,y in enumerate(self.list_K):
                self.diff[i][j]  = K_to_int[x-y]
                if y:
                    self.ratio[i][j] = K_to_int[x/y]
                else:
                    self.ratio[i][j] = UINT_MAX

    def to_difference_family(self, B, check=True):
        r"""
        Given an evenly distributed set ``B`` convert it to a difference family.

        As for any `x\in K^*=H` we have `|\Delta B \cap x
        H^{2(q-1)/(k(k-1))}|=2` (see :class:`EvenlyDistributedSetsBacktracker`),
        the difference family is produced as `\{xB:x\in H^{2(q-1)/(k(k-1))}\}`

        This method is useful if you want to obtain the difference family from
        the output of the iterator.

        INPUT:

        - ``B`` -- an evenly distributed set

        - ``check`` -- boolean (default: ``True``); whether to check the result

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker
            sage: E = EvenlyDistributedSetsBacktracker(Zmod(41),5)
            sage: B = E.an_element(); B
            [0, 1, 13, 38, 31]
            sage: D = E.to_difference_family(B); D
            [[0, 1, 13, 38, 31], [0, 32, 6, 27, 8]]

            sage: from sage.combinat.designs.difference_family import is_difference_family
            sage: is_difference_family(Zmod(41),D,41,5,1)
            True

        Setting ``check`` to ``False`` is much faster::

            sage: timeit("df = E.to_difference_family(B, check=True)") # random
            625 loops, best of 3: 117 µs per loop

            sage: timeit("df = E.to_difference_family(B, check=False)")  # random
            625 loops, best of 3: 1.83 µs per loop
        """
        xe = self.K.multiplicative_generator() ** (self.e)
        df = [[xe**j*b for b in B] for j in range((self.q-1)/(2*self.e))]
        if check:
            from sage.combinat.designs.difference_family import is_difference_family
            if not is_difference_family(self.K, df, self.q, self.k, 1):
                raise RuntimeError("a wrong evenly distributed set was "
                        "produced by the Sage library for the parameters:\n"
                        "  q={}  k={}\n"
                        "Please send an e-mail to "
                        "sage-devel@googlegroups.com".format(self.q, self.k))
        return df

    def an_element(self):
        r"""
        Return an evenly distributed set.

        If there is no such subset raise an
        :class:`~sage.categories.sets_cat.EmptySetError`.

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: E = EvenlyDistributedSetsBacktracker(Zmod(41),5)
            sage: E.an_element()
            [0, 1, 13, 38, 31]

            sage: E = EvenlyDistributedSetsBacktracker(Zmod(61),6)
            sage: E.an_element()
            Traceback (most recent call last):
            ...
            EmptySetError: no 6-evenly distributed set in Ring of integers modulo 61
        """
        from sage.categories.sets_cat import EmptySetError
        it = iter(self)
        try:
            B = next(it)
        except StopIteration:
            raise EmptySetError("no {}-evenly distributed set in {}".format(self.k,self.K))
        self.to_difference_family(B, check=True) # check the validity
        return B

    def __repr__(self):
        r"""
        A string representative.

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: EvenlyDistributedSetsBacktracker(GF(25,'a'), 4)
            4-evenly distributed sets (up to isomorphism) in Finite Field in a of size 5^2
            sage: EvenlyDistributedSetsBacktracker(GF(25,'a'), 4,
            ....:                                  up_to_isomorphism=False)
            4-evenly distributed sets in Finite Field in a of size 5^2
        """
        return "{}-evenly distributed sets {} in {}".format(
                self.k,
                '(up to isomorphism)' if self.up_to_isom else '',
                self.K)

    def cardinality(self):
        r"""
        Return the number of evenly distributed sets.

        Use with precaution as there can be a lot of such sets and this method
        might be very long to answer!

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: E = EvenlyDistributedSetsBacktracker(GF(25,'a'), 4); E
            4-evenly distributed sets (up to isomorphism)
             in Finite Field in a of size 5^2
            sage: E.cardinality()
            4

            sage: E = EvenlyDistributedSetsBacktracker(GF(25,'a'), 4,
            ....:                                      up_to_isomorphism=False)
            sage: E.cardinality()
            40
        """
        cdef n = 0
        self.count = 1
        for a in self:
            n += a
        self.count = 0
        return smallInteger(n)

    def _B_relabelled_copies(self):
        r"""
        Check whether ``self.B`` is minimal among its relabelization.

        If `B=\{b_1,...,b_k\}` is an evenly distributed set and contains `0` and
        `1`, then for any two distinct `i,j` we define `f_{ij} : x \mapsto
        (x-b_j)/(b_i-b_j)` which maps `B` on another evenly distributed set of
        size `k` containing `0` and `1`. For each pair `i,j` we consider check
        whether the set `f_{ij}(B)` is smaller than `B`.

        This is an internal function and should only be call by the backtracker
        implemented in the method ``__iter__``.

        OUTPUT: ``False`` if ``self.B`` is not minimal

        - the list of evenly distributed sets isomorphic to ``self.B``
          given as a list of tuples if ``self.up_to_isom=0`` or list
          containing only ``self.B`` as a tuple if ``self.up_to_isom=1``.

        TESTS::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker
            sage: E = EvenlyDistributedSetsBacktracker(Zmod(13), 4, up_to_isomorphism=True)
            sage: E.cardinality()   # indirect doctest
            1
            sage: E = EvenlyDistributedSetsBacktracker(Zmod(13), 4, up_to_isomorphism=False)
            sage: E.cardinality()   # indirect doctest
            4

        .. NOTE::

            this method is not seriously optimized. The main goal of this backtracker
            is to generate one evenly distributed set. In that case, this method
            will be called only once.
        """
        cdef unsigned int i,j,k,tmp1,tmp2,verify
        cdef list B = [self.B[i] for i in range(1, self.k)]
        B.append(self.q-1)
        cdef list BB = [None]*self.k
        cdef set relabs = set([tuple(B)])

        #  z -> (z - B[i]) / (B[j] - B[i])
        for i in range(self.k):
            for j in range(self.k):
                if i == j:
                    continue
                tmp1 = self.diff[self.B[j]][self.B[i]]

                verify = 0
                for k in range(self.k):
                    if k == i:
                        BB[k] = self.q-1
                    elif k == j:
                        BB[k] = 0
                    else:
                        tmp2 = self.ratio[self.diff[self.B[k]][self.B[i]]][tmp1]
                        if tmp2 == 0 or tmp2 == self.q-1 or tmp2 < self.B[2]:
                            # the backtracker should never build a set which by
                            # relabelling is strictly smaller than B[:3]
                            raise RuntimeError("there is a problem got tmp2={}".format(tmp2,self.B[2]))
                        elif tmp2 == self.B[2]:
                            verify = 1
                        BB[k] = tmp2

                if verify:
                    BB.sort()
                    if BB < B:
                        return False

                if not self.up_to_isom:
                    if not verify:
                        BB.sort()
                    relabs.add(tuple(BB))

        return sorted(relabs)

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __iter__(self):
        r"""
        Iterator through all evenly distributed sets that start with `[0,1]`.

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: E = EvenlyDistributedSetsBacktracker(Zmod(13),4)
            sage: for B in E:
            ....:     print(B)
            [0, 1, 11, 5]
        """
        cdef unsigned int k_m_1 = self.k - 1
        cdef unsigned int q_m_1 = self.q - 1
        cdef unsigned int m = self.m

        # in the list B we store the candidate for being an e.d.s.
        # we always have B[0] = 0 and B[1] = 1
        # because 0 is in B, the cosets of the elements of B must be
        # disjoint.
        cdef unsigned int kk = 2
        cdef unsigned int * B  = self.B
        B[0] = q_m_1  # the element 0 in K
        B[1] = 0      # the element 1 in K

        memset(self.cosets, 0, self.e * sizeof(unsigned int))

        self.cosets[0] = 1  # coset 0 is hit by the difference 1-0

        cdef unsigned int x = m
        while True:
            if self.check:
                self._check_cosets(kk)
                if x < m or x >= q_m_1:
                    raise RuntimeError("got x < m or x > q_m_1 (x={})".format(x))
                if self.cosets[x/m]:
                    raise RuntimeError("got x={} in an already occupied coset".format(x))

            # try to append x
            B[kk] = x
            if self._check_last_element(kk):
                if kk == k_m_1:
                    ans = self._B_relabelled_copies()

                    if self.check and ans:
                        for a in ans:
                            r = [self.list_K[q_m_1]] + [self.list_K[a[r]] for r in range(k_m_1)]
                            self.to_difference_family(r, check=True)

                    if ans is False:
                        pass
                    elif self.count:
                        yield len(ans)
                    else:
                        for a in ans:
                            yield [self.list_K[q_m_1]] + [self.list_K[a[r]] for r in range(k_m_1)]

                    # remove the differences created by x and increment
                    for j in range(kk):
                        self.cosets[ self.diff[x][B[j]] / m ] = 0
                    x += 1
                else:
                    kk += 1
                    x += m - x%m
            else:
                x += 1

            if self.check:
                self._check_cosets(kk)

            # now we determine the next element x to be tested
            while True:
                if kk == 1:
                    return
                elif x == q_m_1:
                    kk -= 1
                    x = self.B[kk]
                    # remove the differences created by x and increment
                    for j in range(kk):
                        self.cosets[ self.diff[x][B[j]] / m ] = 0
                    x += 1
                    if self.check:
                        self._check_cosets(kk)
                elif self.cosets[x / m]:
                    x += m - x%m
                elif kk == 2:
                    if self.min_orb[x] < x:
                        x += 1
                    else:
                        break
                else:
                    if self.min_orb[x] < B[2]:
                        x += 1
                    else:
                        break

    @cython.cdivision(True)
    cdef inline int _check_last_element(self, unsigned int kk) except -1:
        r"""
        Add the element ``x`` to ``B`` in position kk if the resulting set is
        still evenly distributed.

        OUTPUT:

            1 if the element was added, and 0 otherwise.
        """
        cdef unsigned int i, j, x_m_i, x_m_j
        cdef unsigned int m = self.m
        cdef unsigned int * B = self.B
        cdef unsigned int ** diff = self.diff
        cdef unsigned int x = B[kk]

        # We check two things:
        # 1. that the newly created differences x-B[i] will not be in a coset
        #    already occuppied
        #
        # 2. that by applying some automorphisms we will not get an
        #    element smaller than B[2].
        #
        #    We should test all linear functions that send a subset of the form
        #    {x, B[i], B[j]} to some {0, 1, ?}.
        #
        #    Note that if {x, B[i], B[j]} can be mapped to {0, 1, z} by some
        #    function, then it can also be mapped to all {0, 1, z'} where z'=
        #    1/z, 1-z, 1/(1-z), (z-1)/z and z/(z-1). The attribute
        #    'min_orbit[z]' is exactly the minimum among these values.
        #
        #    So, it is enough to test one of these functions. We choose t -> (x
        #    - t)/ (x - B[j]) (that maps x to 0 and B[j] to 1). Its value at
        #    B[i] is just z = (x - B[i]) / (x - B[j]).
        #
        #    In the special case when kk=2, or equivalently when we are testing if x
        #    fits as a new B[2], then we just check that x is the minimum among
        #    {x, 1/x, 1-x, 1/(1-x), (x-1)/x and x/(x-1)}.

        if self.cosets[diff[x][0] / m] == 1:
            return 0

        self.cosets[x / m] = 1
        for i in range(2,kk):
            x_m_i = diff[x][B[i]]

            # 1. check that the difference x-B[i] was not already in an
            #    occuppied coset
            if self.cosets[x_m_i / m]:
                self.cosets[x / m] = 0
                return 0

            # 2. check relabeling
            for j in range(i):
                x_m_j = diff[x][B[j]]
                if self.min_orb[self.ratio[x_m_i][x_m_j]] < B[2]:
                    self.cosets[x / m] = 0
                    return 0

        # Now check that the x-B[i] belongs to distinct cosets
        memcpy(self.t, self.cosets, self.e*sizeof(unsigned int))
        for i in range(1,kk):
            x_m_i = diff[x][B[i]] / m
            if self.t[x_m_i]:
                self.cosets[x / m] = 0
                return 0
            self.t[x_m_i] = 1
        self.t, self.cosets = self.cosets, self.t
        return 1

    @cython.cdivision(True)
    cdef int _check_cosets(self, unsigned int kk) except -1:
        r"""
        Sanity check (only for debug purposes).
        """
        cdef unsigned int i,j
        cdef unsigned int m = self.m
        cdef unsigned int c

        # count the number of elements in self.cosets
        c = 0
        for i in range(self.e):
            c += self.cosets[i]
        if c != (kk * (kk-1)) / 2:
            raise RuntimeError("the number of elements in cosets is wrong! Got {} instead of {}.".format(c, (kk*(kk-1))/2))

        for i in range(kk):
            for j in range(i):
                if self.cosets[ self.diff[self.B[i]][self.B[j]] / m ] != 1:
                    raise RuntimeError("self.cosets misses the difference B[{}]-B[{}]".format(i,j))

        return 0
