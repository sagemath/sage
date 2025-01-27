"""
Baxter permutations
"""
from sage.combinat.permutation import Permutations
from sage.rings.integer_ring import ZZ
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class BaxterPermutations(UniqueRepresentation, Parent):
    r"""
    The combinatorial class of Baxter permutations.

    A Baxter permutation is a permutation avoiding the generalized
    permutation patterns `2-41-3` and `3-14-2`. In other words, a
    permutation `\sigma` is a Baxter permutation if for any subword `u
    := u_1u_2u_3u_4` of `\sigma` such that the letters `u_2` and `u_3`
    are adjacent in `\sigma`, the standardized version of `u` is
    neither `2413` nor `3142`.

    See [Gir2012]_ for a study of Baxter permutations.

    INPUT:

    - ``n`` -- nonnegative integer (default: ``None``); the size of
      the permutations

    OUTPUT:

    Return the combinatorial class of the Baxter permutations of size ``n``
    if ``n`` is not ``None``. Otherwise, return the combinatorial class
    of all Baxter permutations.

    EXAMPLES::

        sage: BaxterPermutations(5)
        Baxter permutations of size 5
        sage: BaxterPermutations()
        Baxter permutations
    """
    @staticmethod
    def __classcall_private__(classe, n=None):
        """
        EXAMPLES::

            sage: BaxterPermutations(5)
            Baxter permutations of size 5
            sage: BaxterPermutations()
            Baxter permutations
        """
        if n is None:
            return BaxterPermutations_all()
        return BaxterPermutations_size(n)


class BaxterPermutations_size(BaxterPermutations):
    r"""
    The enumerated set of Baxter permutations of a given size.

    See :class:`BaxterPermutations` for the definition of Baxter
    permutations.

    EXAMPLES::

        sage: from sage.combinat.baxter_permutations import BaxterPermutations_size
        sage: BaxterPermutations_size(5)
        Baxter permutations of size 5
    """

    def __init__(self, n):
        """
        EXAMPLES::

            sage: from sage.combinat.baxter_permutations import BaxterPermutations_size
            sage: BaxterPermutations_size(5)
            Baxter permutations of size 5
        """
        self.element_class = Permutations(n).element_class
        self._n = ZZ(n)
        from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
        super().__init__(category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.baxter_permutations import BaxterPermutations_size
            sage: BaxterPermutations_size(5)
            Baxter permutations of size 5
        """
        return "Baxter permutations of size %s" % self._n

    def __contains__(self, x):
        r"""
        Return ``True`` if and only if ``x`` is a Baxter permutation of
        size ``self._n``.

        INPUT:

        - ``x`` -- a permutation

        EXAMPLES::

            sage: Permutation([2, 1, 4, 3]) in BaxterPermutations(4)
            True
            sage: Permutation([2, 1, 4, 3]) in BaxterPermutations(5)
            False
            sage: Permutation([3, 1, 4, 2]) in BaxterPermutations(4)
            False
            sage: [len([p for p in Permutations(n) if p in BaxterPermutations(n)]) for n in range(7)]
            [1, 1, 2, 6, 22, 92, 422]
            sage: sorted([p for p in Permutations(6) if p in BaxterPermutations(6)]) == sorted(BaxterPermutations(6).list())
            True
        """
        if x not in Permutations(self._n):
            return False
        for i in range(1, len(x) - 1):
            a = x[i]
            b = x[i + 1]
            if a < b:  # Hunting pattern 3-14-2.
                max_l = 0
                for x_j in x[:i]:
                    if x_j > a and x_j < b and x_j > max_l:
                        max_l = x_j
                min_r = len(x) + 1
                for x_j in x[i + 2:]:
                    if x_j > a and x_j < b and x_j < min_r:
                        min_r = x_j
                if max_l > min_r:
                    return False
            else:  # Hunting pattern 2-41-3.
                min_l = len(x) + 1
                for x_j in x[:i]:
                    if x_j < a and x_j > b and x_j < min_l:
                        min_l = x_j
                max_r = 0
                for x_j in x[i + 2:]:
                    if x_j < a and x_j > b and x_j > max_r:
                        max_r = x_j
                if min_l < max_r:
                    return False
        return True

    def __iter__(self):
        r"""
        Efficient generation of Baxter permutations.

        OUTPUT: an iterator over the Baxter permutations of size ``self._n``

        EXAMPLES::

            sage: BaxterPermutations(4).list()
            [[4, 3, 2, 1], [3, 4, 2, 1], [3, 2, 4, 1], [3, 2, 1, 4], [2, 4, 3, 1],
             [4, 2, 3, 1], [2, 3, 4, 1], [2, 3, 1, 4], [2, 1, 4, 3], [4, 2, 1, 3],
             [2, 1, 3, 4], [1, 4, 3, 2], [4, 1, 3, 2], [1, 3, 4, 2], [1, 3, 2, 4],
             [4, 3, 1, 2], [3, 4, 1, 2], [3, 1, 2, 4], [1, 2, 4, 3], [1, 4, 2, 3],
             [4, 1, 2, 3], [1, 2, 3, 4]]
            sage: [len(BaxterPermutations(n)) for n in range(9)]
            [1, 1, 2, 6, 22, 92, 422, 2074, 10754]

        TESTS::

            sage: all(a in BaxterPermutations(n) for n in range(7)
            ....:     for a in BaxterPermutations(n))
            True

        ALGORITHM:

        The algorithm using generating trees described in [BBMF2008]_ is used.
        The idea is that all Baxter permutations of size `n + 1` can be
        obtained by inserting the letter `n + 1` either just before a left
        to right maximum or just after a right to left maximum of a Baxter
        permutation of size `n`.
        """
        if self._n == 0:
            yield Permutations(0)([])
        elif self._n == 1:
            yield Permutations(1)([1])
        else:
            for b in BaxterPermutations(self._n - 1):
                # Left to right maxima.
                for j in b.reverse().saliances():
                    i = self._n - 2 - j
                    yield Permutations(self._n)(b[:i] + [self._n] + b[i:])
                # Right to left maxima.
                for i in b.saliances():
                    yield Permutations(self._n)(b[:i + 1] + [self._n] + b[i + 1:])

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: BaxterPermutations(4)._an_element_()
            [4, 3, 2, 1]
        """
        return self.first()

    def cardinality(self):
        r"""
        Return the number of Baxter permutations of size ``self._n``.

        For any positive integer `n`, the number of Baxter
        permutations of size `n` equals

        .. MATH::

            \sum_{k=1}^n \dfrac
            {\binom{n+1}{k-1} \binom{n+1}{k} \binom{n+1}{k+1}}
            {\binom{n+1}{1} \binom{n+1}{2}} .

        This is :oeis:`A001181`.

        EXAMPLES::

            sage: [BaxterPermutations(n).cardinality() for n in range(13)]
            [1, 1, 2, 6, 22, 92, 422, 2074, 10754, 58202, 326240, 1882960, 11140560]

            sage: BaxterPermutations(3r).cardinality()
            6
            sage: parent(_)
            Integer Ring
        """
        if self._n == 0:
            return ZZ.one()
        n = self._n + 1
        return sum((n.binomial(k) *
                    n.binomial(k + 1) *
                    n.binomial(k + 2)) // (n * n.binomial(2))
                   for k in range(self._n))


class BaxterPermutations_all(DisjointUnionEnumeratedSets, BaxterPermutations):
    r"""
    The enumerated set of all Baxter permutations.

    See :class:`BaxterPermutations` for the definition of Baxter
    permutations.

    EXAMPLES::

        sage: from sage.combinat.baxter_permutations import BaxterPermutations_all
        sage: BaxterPermutations_all()
        Baxter permutations
    """

    def __init__(self, n=None):
        r"""
        EXAMPLES::

            sage: from sage.combinat.baxter_permutations import BaxterPermutations_all
            sage: BaxterPermutations_all()
            Baxter permutations
        """
        self.element_class = Permutations().element_class
        from sage.sets.non_negative_integers import NonNegativeIntegers
        from sage.sets.family import Family
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(NonNegativeIntegers(),
                                                    BaxterPermutations_size),
                                             facade=False, keepkey=False)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.baxter_permutations import BaxterPermutations_all
            sage: BaxterPermutations_all()
            Baxter permutations
        """
        return "Baxter permutations"

    def __contains__(self, x):
        r"""
        Return ``True`` if and only if ``x`` is a Baxter permutation.

        INPUT:

        - ``x`` -- any object

        EXAMPLES::

            sage: Permutation([4, 2, 1, 7, 3, 8, 5, 6]) in BaxterPermutations()
            False
            sage: Permutation([4, 3, 6, 9, 7, 5, 1, 2, 8]) in BaxterPermutations()
            True

        TESTS::

            sage: 42 in BaxterPermutations()
            False
        """
        if x not in Permutations():
            return False
        return x in BaxterPermutations(len(x))

    def to_pair_of_twin_binary_trees(self, p):
        r"""
        Apply a bijection between Baxter permutations of size ``self._n``
        and the set of pairs of twin binary trees with ``self._n`` nodes.

        INPUT:

        - ``p`` -- a Baxter permutation

        OUTPUT:

        The pair of twin binary trees `(T_L, T_R)` where `T_L`
        (resp. `T_R`) is obtained by inserting the letters of ``p`` from
        left to right (resp. right to left) following the binary search
        tree insertion algorithm. This is called the *Baxter P-symbol*
        in [Gir2012]_ Definition 4.1.

        .. NOTE::

            This method only works when ``p`` is a permutation. For words
            with repeated letters, it would return two "right binary
            search trees" (in the terminology of [Gir2012]_), which conflicts
            with the definition in [Gir2012]_.

        EXAMPLES::

            sage: BP = BaxterPermutations()
            sage: BP.to_pair_of_twin_binary_trees(Permutation([]))                      # needs sage.graphs
            (., .)
            sage: BP.to_pair_of_twin_binary_trees(Permutation([1, 2, 3]))               # needs sage.graphs
            (1[., 2[., 3[., .]]], 3[2[1[., .], .], .])
            sage: BP.to_pair_of_twin_binary_trees(Permutation([3, 4, 1, 2]))            # needs sage.graphs
            (3[1[., 2[., .]], 4[., .]], 2[1[., .], 4[3[., .], .]])
        """
        from sage.combinat.binary_tree import LabelledBinaryTree
        left = LabelledBinaryTree(None)
        right = LabelledBinaryTree(None)
        for a in p:
            left = left.binary_search_insert(a)
        for a in reversed(p):
            right = right.binary_search_insert(a)
        return (left, right)
