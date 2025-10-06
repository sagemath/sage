r"""
Decorated permutations

AUTHORS:

- Martin Rubey (2020): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2020 Martin Rubey <martin.rubey at tuwien.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.arith.misc import factorial
from sage.rings.integer import Integer
from sage.combinat.permutation import Permutations
from sage.combinat.subset import Subsets
from sage.combinat.colored_permutations import SignedPermutations
from sage.structure.list_clone import ClonableArray


class DecoratedPermutation(ClonableArray,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A decorated permutation.

    A decorated permutation is a signed permutation where all
    non-fixed points have positive sign.
    """
    @staticmethod
    def __classcall_private__(cls, pi):
        """
        Create a decorated permutation.

        EXAMPLES::

            sage: DecoratedPermutation([2, 1, 3])
            [2, 1, 3]

            sage: DecoratedPermutation([2, 1, -3])
            [2, 1, -3]

        TESTS:

        Check that hashing and comparison works::

            sage: S = DecoratedPermutations(3)
            sage: elt1 = S([2, 1, -3])
            sage: elt2 = DecoratedPermutation([2, 1, -3])
            sage: elt1 == elt2
            True

            sage: elt1 == [2, 1, -3]
            False

            sage: elt2 = DecoratedPermutation([2, 1, 3])
            sage: elt1 != elt2
            True

            sage: hash(elt1)                                                    # random
            915443076393556996
        """
        pi = list(pi)
        return DecoratedPermutations(len(pi))(pi)

    def __init__(self, parent, pi, check=True):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = DecoratedPermutations(3)
            sage: elt = S([2, 1, -3])
            sage: TestSuite(elt).run()
        """
        ClonableArray.__init__(self, parent, pi, check=check)

    def check(self):
        """
        Check that ``self`` is a valid decorated permutation.

        EXAMPLES::

            sage: S = DecoratedPermutations(3)
            sage: elt = S([2, 1, -3])
            sage: elt.check()
            sage: elt = S([2, -1, 3])
            Traceback (most recent call last):
            ...
            ValueError: invalid decorated permutation
        """
        if self not in self.parent():
            raise ValueError("{} is not a decorated permutation".format(self))

    def size(self):
        """
        Return the size of the decorated permutation.

        EXAMPLES::

            sage: DecoratedPermutation([2, 1, -3]).size()
            3
        """
        return len(self)

    def to_signed_permutation(self):
        """
        Return ``self`` as a signed permutation.

        EXAMPLES::

            sage: DecoratedPermutation([2, 1, -3]).to_signed_permutation()
            [2, 1, -3]
        """
        return SignedPermutations(len(self))(list(self))


class DecoratedPermutations(UniqueRepresentation, Parent):
    r"""
    Class of all decorated permutations of `n`.

    A decorated permutation is a signed permutation where all
    non-fixed points have positive sign.

    INPUT:

    - ``n`` -- integer; the size of the decorated permutations

    EXAMPLES:

    This will create an instance to manipulate the decorated
    permutations of size 3::

        sage: S = DecoratedPermutations(3); S
        Decorated permutations of size 3
        sage: S.cardinality()
        16
    """

    def __init__(self, n):
        r"""
        Initialize ``self``.

        TESTS::

            sage: S = DecoratedPermutations(4)
            sage: TestSuite(S).run()
        """
        self._n = n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: DecoratedPermutations(4)
            Decorated permutations of size 4
        """
        return "Decorated permutations of size %s" % self._n

    def __contains__(self, pi):
        """
        Check if ``pi`` is in ``self``.

        TESTS::

            sage: S = DecoratedPermutations(3)
            sage: [2, 1, -3] in S
            True
            sage: [2, -1, 3] in S
            False
        """
        if isinstance(pi, DecoratedPermutation):
            return len(pi) == self._n

        values = list(pi)
        if len(values) != self._n:
            return False
        abs_values = [abs(v) for v in values]
        for i, (v, abs_v) in enumerate(zip(values, abs_values), 1):
            if i != abs_v and v < 0:
                return False
        return sorted(abs_values) == list(range(1, self._n + 1))

    def _element_constructor_(self, pi, check=True):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: S = DecoratedPermutations(3)
            sage: elt = S([2, 1, -3])
            sage: elt.parent() is S
            True
        """
        if isinstance(pi, DecoratedPermutation):
            if pi.parent() is self:
                return pi
            raise ValueError("cannot convert between decorated permutations of different sizes")

        pi = tuple(pi)
        if check and pi not in self:
            raise ValueError('invalid decorated permutation')
        return self.element_class(self, pi)

    Element = DecoratedPermutation

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: S = DecoratedPermutations(3)
            sage: S._an_element_()
            [1, 2, 3]
        """
        return self.element_class(self, list(range(1, self._n + 1)))

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of decorated permutations of size `n` is equal to

        .. MATH::

            \sum_{k=0^n} \frac{n!}{k!}

        EXAMPLES::

            sage: [DecoratedPermutations(n).cardinality() for n in range(11)]
            [1, 2, 5, 16, 65, 326, 1957, 13700, 109601, 986410, 9864101]
        """
        return Integer(sum(factorial(self._n) // factorial(k)
                           for k in range(self._n + 1)))

    def __iter__(self):
        r"""
        Iterator on the decorated permutations of size `n`.

        TESTS::

            sage: S = DecoratedPermutations(2); S.list()
            [[1, 2], [-1, 2], [1, -2], [-1, -2], [2, 1]]
            sage: sum(1 for a in S)
            5
        """
        for sigma in Permutations(self._n):
            F = sigma.fixed_points()
            for X in Subsets(F):
                tau = list(sigma)
                for i in X:
                    tau[i - 1] = -tau[i - 1]
                yield DecoratedPermutations(self._n)(tau)
