# sage.doctest: needs sage.combinat sage.modules
r"""
Similarity class types of matrices with entries in a finite field

The notion of a matrix conjugacy class type was introduced by J. A. Green in
[Green55]_, in the context of computing the irreducible characters of finite
general linear groups. The class types are equivalence classes of similarity
classes of square matrices with entries in a finite field which, roughly
speaking, have the same qualitative properties.

For example, all similarity classes of the same class type have centralizers of
the same cardinality and the same degrees of elementary divisors. Qualitative
properties of similarity classes such as semisimplicity and regularity descend
to class types.

The most important feature of similarity class types is that, for any `n`, the
number of similarity class types of `n\times n` matrices is independent of `q`.
This makes it possible to perform many combinatorial calculations treating `q`
as a formal variable.

In order to define similarity class types, recall that similarity classes of
`n\times n` matrices with entries in `\GF{q}` correspond to functions

.. MATH::

    c: \mathrm{Irr}\GF{q[t]} \to \Lambda

such that

.. MATH::

    \sum_{f\in \mathrm{Irr}\GF{q[t]}} |c(f)|\deg f = n,

where we denote the set of irreducible monic polynomials in `\GF{q[t]}`
by `\mathrm{Irr}\GF{q[t]}`, the set of all partitions by `\Lambda`, and
the size of `\lambda \in \Lambda` by `|\lambda|`.

Similarity classes indexed by functions `c_1` and `c_2` as above are said to be
of the same type if there exists a degree-preserving self-bijection `\sigma` of
`\mathrm{Irr}\GF{q[t]}` such that `c_2 = c_1\circ \sigma`. Thus, the type
of `c` remembers only the degrees of the polynomials (and not the polynomials
themselves) for which `c` takes a certain value `\lambda`. Replacing each
irreducible polynomial of degree `d` for which `c` takes a non-trivial value
`\lambda` by the pair `(d, \lambda)`, we obtain a multiset of such pairs.
Clearly, `c_1` and `c_2` have the same type if and only if these multisets are
equal. Thus a similarity class type may be viewed as a multiset of pairs of the
form `(d, \lambda)`.

For `2 \times 2` matrices there are four types::

    sage: for tau in SimilarityClassTypes(2):
    ....:    print(tau)
    [[1, [1]], [1, [1]]]
    [[1, [2]]]
    [[1, [1, 1]]]
    [[2, [1]]]

These four types correspond to the regular split semisimple matrices, the
non-semisimple matrices, the central matrices and the irreducible matrices
respectively.

For any matrix `A` in a given similarity class type, it is possible to calculate
the number elements in the similarity class of `A`, the dimension of the algebra
of matrices in `M_n(A)` that commute  with `A`, and the cardinality of the
subgroup of `GL_n(\GF{q})` that commute with `A`. For each similarity
class type, it is also possible to compute the number of classes of that type
(and hence, the total number of matrices of that type). All these calculations
treat the cardinality `q` of the finite field as a formal variable::

    sage: M = SimilarityClassType([[1, [1]], [1, [1]]])
    sage: M.class_card()
    q^2 + q
    sage: M.centralizer_algebra_dim()
    2
    sage: M.centralizer_group_card()
    q^2 - 2*q + 1
    sage: M.number_of_classes()
    1/2*q^2 - 1/2*q
    sage: M.number_of_matrices()
    1/2*q^4 - 1/2*q^2

We now describe two applications of similarity class types.

We say that an `n \times n` matrix has rational canonical form type `\lambda` for
some partition `\lambda` of `n` if the diagonal blocks in the rational canonical
form have sizes given by the parts of `\lambda`. Thus the matrices with rational
canonical type `(n)` are the regular ones, while the matrices with rational
canonical type `(1^n)` are the central ones.

Using similarity class types, it becomes easy to get a formula for the number of
matrices with a given rational canonical type::

    sage: def matrices_with_rcf(la):
    ....:    return sum([tau.number_of_matrices() for tau in filter(lambda tau:tau.rcf()==la, SimilarityClassTypes(la.size()))])
    sage: matrices_with_rcf(Partition([2,1]))
    q^6 + q^5 + q^4 - q^3 - q^2 - q

Similarity class types can also be used to calculate the number of simultaneous
similarity classes of `k`-tuples of `n\times n` matrices with entries in
`\GF{q}` by using Burnside's lemma::

    sage: from sage.combinat.similarity_class_type import order_of_general_linear_group, centralizer_algebra_dim
    sage: q = ZZ['q'].gen()
    sage: def simultaneous_similarity_classes(n, k):
    ....:     return SimilarityClassTypes(n).sum(lambda la: q**(k*centralizer_algebra_dim(la)), invertible = True)/order_of_general_linear_group(n)
    sage: simultaneous_similarity_classes(3, 2)
    q^10 + q^8 + 2*q^7 + 2*q^6 + 2*q^5 + q^4

Similarity class types can be used to compute the coefficients of generating
functions coming from the cycle index type techniques of Kung and Stong (see
Morrison [Morrison06]_).

They can also be used to compute the number of invariant subspaces for a matrix
over a finite field of any given dimension. For this we use the elegant recursive
formula of Ramaré [R17]_ (see also [PR22]_).

Along with the results of [PSS13]_, similarity class types can be used to
calculate the number of similarity classes of matrices of order `n` with entries
in a principal ideal local ring of length two with residue field of cardinality
`q` with centralizer of any given cardinality up to `n = 4`. Among these, the
classes which are selftranspose can also be counted::

    sage: from sage.combinat.similarity_class_type import matrix_centralizer_cardinalities_length_two
    sage: list(matrix_centralizer_cardinalities_length_two(3))
    [(q^6 - 3*q^5 + 3*q^4 - q^3, 1/6*q^6 - 1/2*q^5 + 1/3*q^4),
    (q^6 - 2*q^5 + q^4, q^5 - q^4),
    (q^8 - 3*q^7 + 3*q^6 - q^5, 1/2*q^5 - q^4 + 1/2*q^3),
    (q^8 - 2*q^7 + q^6, q^4 - q^3),
    (q^10 - 2*q^9 + 2*q^7 - q^6, q^4 - q^3),
    (q^8 - q^7 - q^6 + q^5, 1/2*q^5 - q^4 + 1/2*q^3),
    (q^6 - q^5 - q^4 + q^3, 1/2*q^6 - 1/2*q^5),
    (q^6 - q^5, q^4),
    (q^10 - 2*q^9 + q^8, q^3),
    (q^8 - 2*q^7 + q^6, q^4 - q^3),
    (q^8 - q^7, q^3 + q^2),
    (q^12 - 3*q^11 + 3*q^10 - q^9, 1/6*q^4 - 1/2*q^3 + 1/3*q^2),
    (q^12 - 2*q^11 + q^10, q^3 - q^2),
    (q^14 - 2*q^13 + 2*q^11 - q^10, q^3 - q^2),
    (q^12 - q^11 - q^10 + q^9, 1/2*q^4 - 1/2*q^3),
    (q^12 - q^11, q^2),
    (q^14 - 2*q^13 + q^12, q^2),
    (q^18 - q^17 - q^16 + q^14 + q^13 - q^12, q^2),
    (q^12 - q^9, 1/3*q^4 - 1/3*q^2),
    (q^6 - q^3, 1/3*q^6 - 1/3*q^4)]

REFERENCES:

.. [Green55] Green, J. A.  *The characters of the finite general linear groups*.
   Trans. Amer. Math. Soc.  80  (1955), 402--447.
   :doi:`10.1090/S0002-9947-1955-0072878-2`

.. [Morrison06] Morrison, Kent E.
   *Integer sequences and matrices over finite fields*.
   J. Integer Seq. 9 (2006), no. 2, Article 06.2.1, 28 pp.
   https://cs.uwaterloo.ca/journals/JIS/VOL9/Morrison/morrison37.html

.. [PSS13] Prasad, A., Singla, P., and Spallone, S., *Similarity of matrices
   over local rings of length two*. :arxiv:`1212.6157`

.. [PR22] Prasad, A., Ram, S., *Splitting subspaces and a finite field
   interpretation of the Touchard-Riordan formula*. :arxiv:`2205.11076`

.. [R17] Ramaré, O., *Rationality of the zeta function of the subgroups of
   abelian p-groups*. Publ. Math. Debrecen 90.1-2.
   :doi:`10.5486/PMD.2017.7466`

AUTHOR:

- Amritanshu Prasad (2013-07-18): initial implementation

- Amritanshu Prasad (2013-09-09): added functions for similarity classes over
  rings of length two

- Amritanshu Prasad (2022-07-31): added computation of similarity class type of
  a given matrix and invariant subspace generating function
"""
# ****************************************************************************
#       Copyright (C) 2013 Amritanshu Prasad <amri@imsc.res.in>
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

from itertools import chain, product

from sage.arith.misc import divisors, factorial, moebius
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.combinat import CombinatorialElement
from sage.combinat.misc import IterableFunctionCall
from sage.combinat.partition import Partitions, Partition
from sage.misc.cachefunc import cached_in_parent_method, cached_function
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.misc_c import prod
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element import Element, Matrix
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


@cached_function
def fq(n, q=None):
    r"""
    Return `(1-q^{-1}) (1-q^{-2}) \cdots (1-q^{-n})`.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``q`` -- integer or an indeterminate

    OUTPUT: a rational function in ``q``

    EXAMPLES::

        sage: from sage.combinat.similarity_class_type import fq
        sage: fq(0)
        1
        sage: fq(3)
        (q^6 - q^5 - q^4 + q^2 + q - 1)/q^6
    """
    if q is None:
        q = ZZ['q'].gen()
    return prod(1 - q**(-i - 1) for i in range(n))


@cached_function
def primitives(n, invertible=False, q=None):
    """
    Return the number of similarity classes of simple matrices
    of order ``n`` with entries in a finite field of order ``q``.
    This is the same as the number of irreducible polynomials
    of degree `d`.

    If ``invertible`` is ``True``, then only the number of
    similarity classes of invertible matrices is returned.

    .. NOTE::

        All primitive classes are invertible unless ``n`` is `1`.

    INPUT:

    - ``n`` -- positive integer

    - ``invertible`` -- boolean; if set, only number of nonzero classes is returned

    - ``q`` -- integer or an indeterminate

    OUTPUT: a rational function of the variable ``q``

    EXAMPLES::

        sage: from sage.combinat.similarity_class_type import primitives
        sage: primitives(1)
        q
        sage: primitives(1, invertible = True)
        q - 1
        sage: primitives(4)
        1/4*q^4 - 1/4*q^2
        sage: primitives(4, invertible = True)
        1/4*q^4 - 1/4*q^2
    """
    if q is None:
        q = QQ['q'].gen()
    p = sum(moebius(n // d) * q**d for d in divisors(n)) / n
    if invertible and n == 1:
        return p - 1
    else:
        return p


@cached_function
def order_of_general_linear_group(n, q=None):
    r"""
    Return the cardinality of the group of `n \times n` invertible matrices
    with entries in a field of order ``q``.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``q`` -- integer or an indeterminate

    EXAMPLES::

        sage: from sage.combinat.similarity_class_type import order_of_general_linear_group
        sage: order_of_general_linear_group(0)
        1
        sage: order_of_general_linear_group(2)
        q^4 - q^3 - q^2 + q
    """
    if q is None:
        q = ZZ['q'].gen()
    return prod([q**n - q**i for i in range(n)])


@cached_function
def centralizer_algebra_dim(la):
    r"""
    Return the dimension of the centralizer algebra in `M_n(\GF{q})`
    of a nilpotent matrix whose Jordan blocks are given by ``la``.

    EXAMPLES::

        sage: from sage.combinat.similarity_class_type import centralizer_algebra_dim
        sage: centralizer_algebra_dim(Partition([2, 1]))
        5

    .. NOTE::

        If it is a list, ``la`` is expected to be sorted in decreasing order.
    """
    return sum([(2 * i + 1) * la[i] for i in range(len(la))])


@cached_function
def centralizer_group_cardinality(la, q=None):
    r"""
    Return the cardinality of the centralizer group in `GL_n(\GF{q})`
    of a nilpotent matrix whose Jordan blocks are given by ``la``.

    INPUT:

    - ``lambda`` -- a partition

    - ``q`` -- an integer or an indeterminate

    OUTPUT: a polynomial function of ``q``

    EXAMPLES::

        sage: from sage.combinat.similarity_class_type import centralizer_group_cardinality
        sage: q = ZZ['q'].gen()
        sage: centralizer_group_cardinality(Partition([2, 1]))
        q^5 - 2*q^4 + q^3
    """
    if q is None:
        q = ZZ['q'].gen()
    return q**centralizer_algebra_dim(la)*prod([fq(m, q=q) for m in la.to_exp()])


def invariant_subspace_generating_function(la, q=None, t=None):
    """
    Return the invariant subspace generating function of a nilpotent matrix with
    Jordan block sizes given by ``la``.

    INPUT:

    - ``la`` -- a partition
    - ``q`` -- (optional) an integer or an inderminate
    - ``t`` -- (optional) an indeterminate

    OUTPUT: a polynomial in ``t`` whose coefficients are polynomials in ``q``

    EXAMPLES::

        sage: from sage.combinat.similarity_class_type import invariant_subspace_generating_function
        sage: invariant_subspace_generating_function([2,2])
        t^4 + (q + 1)*t^3 + (q^2 + q + 1)*t^2 + (q + 1)*t + 1
    """
    if q is None:
        q = PolynomialRing(QQ,'q').gen()
    S = q.parent()
    if t is None:
        t = PolynomialRing(S,'t').gen()
    R = t.parent()
    Rff = R.fraction_field()
    if not la:
        return Rff(1)
    u = invariant_subspace_generating_function(la[1:], q=q, t=t)
    return R((t**(la[0]+1) * q**(sum(la[1:])) * u.substitute(t=t/q) - u.substitute(t=t*q)) / (t - 1))


class PrimarySimilarityClassType(Element,
                                 metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A primary similarity class type is a pair consisting of a partition and a positive
    integer.

    For a partition `\lambda` and a positive integer `d`, the primary similarity
    class type `(d, \lambda)` represents similarity classes of square matrices
    of order `|\lambda| \cdot d` with entries in a finite field of order `q`
    which correspond to the `\GF{q[t]}`-module

    .. MATH::

        \frac{\GF{q[t]}}{p(t)^{\lambda_1} } \oplus
        \frac{\GF{q[t]}}{p(t)^{\lambda_2}} \oplus \dotsb

    for some irreducible polynomial `p(t)` of degree `d`.
    """
    @staticmethod
    def __classcall_private__(cls, deg, par):
        r"""
        Create a primary similarity class type.

        EXAMPLES::

            sage: PrimarySimilarityClassType(2, [3, 2, 1])
            [2, [3, 2, 1]]

        The parent class is the class of primary similarity class types of order
        `d |\lambda|`::

            sage: PT = PrimarySimilarityClassType(2, [3, 2, 1])
            sage: PT.parent().size()
            12
        """
        par = Partition(par)
        P = PrimarySimilarityClassTypes(par.size()*deg)
        return P(deg, par)

    def __init__(self, parent, deg, par):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: elt =  PrimarySimilarityClassType(2, [3, 2, 1])
            sage: TestSuite(elt).run()
        """
        self._deg = deg
        self._par = par
        Element.__init__(self, parent)

    def __repr__(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: PrimarySimilarityClassType(2, [3, 2, 1])
            [2, [3, 2, 1]]
        """
        return "%s" % ([self._deg, self._par],)

    def __hash__(self):
        r"""
        TESTS::

            sage: PT1 = PrimarySimilarityClassType(2, [3, 2, 1])
            sage: PT2 = PrimarySimilarityClassType(3, [3, 2, 1])
            sage: PT3 = PrimarySimilarityClassType(2, [4, 2, 1])
            sage: hash(PT1) == hash(PrimarySimilarityClassType(2, [3, 2, 1]))
            True
            sage: abs(hash(PT1) - hash(PT2)) == 1
            True
            sage: hash(PT1) == hash(PT3)
            False
            sage: hash(PT2) == hash(PT3)
            False
        """
        return hash(self._deg) ^ hash(tuple(self._par))

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: PT1 =  PrimarySimilarityClassType(2, [3, 2, 1])
            sage: PT2 =  PrimarySimilarityClassType(2, Partition([3, 2, 1]))
            sage: PT1 == PT2
            True
            sage: PT3 =  PrimarySimilarityClassType(3, [3, 2, 1])
            sage: PT1 == PT3
            False
            sage: PT4 =  PrimarySimilarityClassType(2, [3, 2, 1, 0])
            sage: PT1 == PT4
            True
            sage: PT5 = PrimarySimilarityClassType(2, [4, 2, 1])
            sage: PT1 == PT5
            False
        """
        return isinstance(other, PrimarySimilarityClassType) and \
            self.degree() == other.degree() and \
            self.partition() == other.partition()

    def __ne__(self, other):
        r"""
        TESTS::

            sage: PT1 =  PrimarySimilarityClassType(2, [3, 2, 1])
            sage: PT2 =  PrimarySimilarityClassType(2, Partition([3, 2, 1]))
            sage: PT1 != PT2
            False
            sage: PT3 =  PrimarySimilarityClassType(3, [3, 2, 1])
            sage: PT1 != PT3
            True
        """
        return not isinstance(other, PrimarySimilarityClassType) or \
            self.degree() != other.degree() or \
            self.partition() != other.partition()

    def size(self):
        """
        Return the size of ``self``.

        EXAMPLES::

            sage: PT = PrimarySimilarityClassType(2, [3, 2, 1])
            sage: PT.size()
            12
        """
        return self.parent().size()

    def degree(self):
        """
        Return degree of ``self``.

        EXAMPLES::

            sage: PT = PrimarySimilarityClassType(2, [3, 2, 1])
            sage: PT.degree()
            2
        """
        return self._deg

    def partition(self):
        """
        Return partition corresponding to ``self``.

        EXAMPLES::

            sage: PT = PrimarySimilarityClassType(2, [3, 2, 1])
            sage: PT.partition()
            [3, 2, 1]
        """
        return Partition(self._par)

    def centralizer_algebra_dim(self):
        r"""
        Return the dimension of the algebra of matrices which commute with a
        matrix of type ``self``.

        For a partition `(d, \lambda)` this dimension is given by
        `d(\lambda_1 + 3\lambda_2 + 5\lambda_3 + \cdots)`.

        EXAMPLES::

            sage: PT = PrimarySimilarityClassType(2, [3, 2, 1])
            sage: PT.centralizer_algebra_dim()
            28
        """
        return self.degree()*centralizer_algebra_dim(self.partition())

    @cached_in_parent_method
    def statistic(self, func, q=None):
        r"""
        Return `n_{\lambda}(q^d)` where `n_{\lambda}` is the value returned by
        ``func`` upon input `\lambda`, if ``self`` is `(d, \lambda)`.

        EXAMPLES::

            sage: PT = PrimarySimilarityClassType(2, [3, 1])
            sage: q = ZZ['q'].gen()
            sage: PT.statistic(lambda la:q**la.size(), q = q)
            q^8
        """
        if q is None:
            q = ZZ['q'].gen()
        return q.parent()(func(self.partition()).substitute(q=q**self.degree()))

    @cached_in_parent_method
    def centralizer_group_card(self, q=None):
        """
        Return the cardinality of the centralizer group of a matrix of type
        ``self`` in a field of order ``q``.

        INPUT:

        - ``q`` -- integer or an indeterminate

        EXAMPLES::

            sage: PT = PrimarySimilarityClassType(1, [])
            sage: PT.centralizer_group_card()
            1
            sage: PT = PrimarySimilarityClassType(2, [1, 1])
            sage: PT.centralizer_group_card()
            q^8 - q^6 - q^4 + q^2
        """
        if q is None:
            q = FractionField(ZZ['q']).gen()
        return self.statistic(centralizer_group_cardinality, q=q)

    def invariant_subspace_generating_function(self, q=None, t=None):
        """
        Return the invariant subspace generating function of ``self``.

        INPUT:

        - ``q`` -- (optional) an integer or an inderminate
        - ``t`` -- (optional) an indeterminate

        EXAMPLES::

            sage: PrimarySimilarityClassType(1, [2, 2]).invariant_subspace_generating_function()
            t^4 + (q + 1)*t^3 + (q^2 + q + 1)*t^2 + (q + 1)*t + 1
        """
        if q is None:
            q = PolynomialRing(QQ, 'q').gen()
        S = q.parent()
        if t is None:
            t = PolynomialRing(S, 't').gen()
        return invariant_subspace_generating_function(self.partition()).substitute(q=q**self.degree(), t=t**self.degree())


class PrimarySimilarityClassTypes(UniqueRepresentation, Parent):
    r"""
    All primary similarity class types of size ``n`` whose degree is greater
    than that of ``min`` or whose degree is that of ``min`` and  whose partition
    is less than of ``min`` in lexicographic order.

    A primary similarity class type of size `n` is a pair `(\lambda, d)`
    consisting of a partition `\lambda` and a positive integer `d` such that
    `|\lambda| d = n`.

    INPUT:

    - ``n`` -- positive integer
    - ``min`` -- a primary matrix type of size ``n``

    EXAMPLES:

    If ``min`` is not specified, then the class of all primary similarity class
    types of size ``n`` is created::

        sage: PTC = PrimarySimilarityClassTypes(2)
        sage: for PT in PTC:
        ....:     print(PT)
        [1, [2]]
        [1, [1, 1]]
        [2, [1]]

    If ``min`` is specified, then the class consists of only those primary
    similarity class types whose degree is greater than that of ``min`` or whose
    degree is that of ``min`` and  whose partition is less than of ``min`` in
    lexicographic order::

        sage: PTC = PrimarySimilarityClassTypes(2, min = PrimarySimilarityClassType(1, [1, 1]))
        sage: for PT in PTC:
        ....:     print(PT)
        [1, [1, 1]]
        [2, [1]]
    """
    @staticmethod
    def __classcall_private__(cls, n, min=None):
        r"""
        Create the class of vector partitions of ``vec`` where all parts
        are greater than or equal to the vector ``min``.

        EXAMPLES::

            sage: PTC1 = PrimarySimilarityClassTypes(2)
            sage: PTC2 = PrimarySimilarityClassTypes(2, min = PrimarySimilarityClassType(1, [2]))
            sage: PTC1 is PTC2
            True
        """
        if min is None:
            min = (ZZ.one(), Partition([n]))
        elif isinstance(min, PrimarySimilarityClassType):
            min = (min.degree(), min.partition())
        elif len(min) == 2:
            min = (min[0], Partition(min[1]))
        else:
            raise ValueError("min must be a PrimarySimilarityClassType")
        return super().__classcall__(cls, n, min)

    def __init__(self, n, min):
        r"""
        Initialize ``self``.

        TESTS::

            sage: PTC = PrimarySimilarityClassTypes(2)
            sage: TestSuite(PTC).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._n = n
        self._min = min

    def _element_constructor_(self, deg, par):
        """
        Construct an element of ``self``.

        INPUT:

        - ``deg`` -- positive integer

        - ``par`` -- a partition

        EXAMPLES::

            sage: PTC = PrimarySimilarityClassTypes(2)
            sage: elt = PTC(1, [1, 1]); elt
            [1, [1, 1]]
            sage: elt.parent() is PTC
            True
        """
        return self.element_class(self, deg, par)

    Element = PrimarySimilarityClassType

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: PTC = PrimarySimilarityClassTypes(2)
            sage: PTC.cardinality()
            3
        """
        n = self._n
        if self._min[0].divides(n):
            for par in Partitions(n // self._min[0], starting=self._min[1]):
                yield self.element_class(self, self._min[0], par)
        for d in (d for d in divisors(n) if d > self._min[0]):
            for par in Partitions(n // d):
                yield self.element_class(self, d, par)

    def size(self):
        r"""
        Return size of elements of ``self``.

        The size of a primary similarity class type `(d, \lambda)` is
        `d |\lambda|`.

        EXAMPLES::

            sage: PTC = PrimarySimilarityClassTypes(2)
            sage: PTC.size()
            2
        """
        return self._n

###############################################################################

###############################################################################


class SimilarityClassType(CombinatorialElement):
    r"""
    A similarity class type.

    A matrix type is a multiset of primary similarity class types.

    INPUT:

    - ``tau`` -- list of primary similarity class types or a square matrix
      over a finite field

    EXAMPLES::

        sage: tau1 = SimilarityClassType([[3, [3, 2, 1]], [2, [2, 1]]]); tau1
        [[2, [2, 1]], [3, [3, 2, 1]]]

        sage: SimilarityClassType(Matrix(GF(2), [[1,1],[0,1]]))
        [[1, [2]]]
    """
    @staticmethod
    def __classcall_private__(cls, tau):
        """
        Create a similarity class type.

        EXAMPLES:

        The input can be a list of lists or a list of primary similarity class
        types, and the order in which this list is given does not matter::

            sage: tau1 = SimilarityClassType([[3, [3, 2, 1]], [2, [2, 1]]]); tau1
            [[2, [2, 1]], [3, [3, 2, 1]]]
            sage: types = [PrimarySimilarityClassType(2, [2, 1]), PrimarySimilarityClassType(3, [3, 2, 1])]
            sage: tau2 = SimilarityClassType(types)
            sage: tau1 == tau2
            True

        The input can also be a matrix with entries in a finite field::

            sage: SimilarityClassType(Matrix(GF(2), [[1,1],[0,1]]))
            [[1, [2]]]

        The parent class is the class of similarity class types of the sum of
        the sizes of the primary matrix types in ``tau``::

            sage: tau = SimilarityClassType([[3, [3, 2, 1]], [2, [2, 1]]])
            sage: tau.parent().size()
            24
        """
        if isinstance(tau, Matrix):
            n = tau.nrows()
            F = tau.base_ring()
            R = PolynomialRing(F, 't')
            t = R.gen()
            S = (t - tau).smith_form(transformation=False)
            L = [S[i,i] for i in range(n-1, -1, -1) if S[i,i]]
            f = [dict(list(p.factor())) for p in L]
            d = {p: Partition([h[p] for h in f if p in h]) for p in f[0]}
            return SimilarityClassType([[p.degree(), d[p]] for p in d])
        else:
            ret = []
            for l in tau:
                if isinstance(l, PrimarySimilarityClassType):
                    ret.append(l)
                else:
                    ret.append(PrimarySimilarityClassType(*l))
            n = sum([PT.size() for PT in ret])
            T = SimilarityClassTypes(n)
            return T(tau)

    def __init__(self, parent, tau):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: elt =  SimilarityClassType([[3, [3, 2, 1]], [2, [2, 1]]])
            sage: TestSuite(elt).run()
        """
        tau = sorted(tau, key=lambda PT: (PT.degree(), PT.partition()))
        CombinatorialElement.__init__(self, parent, tau)

    def size(self):
        """
        Return the sum of the sizes of the primary parts of ``self``.

        EXAMPLES::

            sage: tau = SimilarityClassType([[3, [3, 2, 1]], [2, [2, 1]]])
            sage: tau.size()
            24
        """
        return self.parent().size()

    def centralizer_algebra_dim(self):
        """
        Return the dimension of the algebra of matrices which commute with a
        matrix of type ``self``.

        EXAMPLES::

            sage: tau = SimilarityClassType([[1, [1]], [1, [1]]])
            sage: tau.centralizer_algebra_dim()
            2
        """
        return sum([PT.centralizer_algebra_dim() for PT in self])

    def centralizer_group_card(self, q=None):
        r"""
        Return the cardinality of the group of matrices in `GL_n(\GF{q})`
        which commute with a matrix of type ``self``.

        INPUT:

        - ``q`` -- integer or an indeterminate

        EXAMPLES::

            sage: tau = SimilarityClassType([[1, [1]], [1, [1]]])
            sage: tau.centralizer_group_card()
            q^2 - 2*q + 1
        """
        return prod([PT.centralizer_group_card(q=q) for PT in self])

    def as_partition_dictionary(self):
        r"""
        Return a dictionary whose keys are the partitions of types occurring in
        ``self`` and the value at the key `\lambda` is the partition formed by
        sorting the degrees of primary types with partition `\lambda`.

        EXAMPLES::

            sage: tau = SimilarityClassType([[1, [1]], [1, [1]]])
            sage: tau.as_partition_dictionary()
            {[1]: [1, 1]}
        """
        D = {}
        for PT in self:
            if PT.partition() in D:
                D[PT.partition()] = Partition(sorted(D[PT.partition()] + [PT.degree()]))
            else:
                D[PT.partition()] = Partition([PT.degree()])
        return D

    def number_of_classes(self, invertible=False, q=None):
        """
        Return the number of similarity classes of matrices of type ``self``.

        INPUT:

        - ``invertible`` -- boolean; return number of invertible classes if set
          to ``True``

        - ``q`` -- integer or an indeterminate

        EXAMPLES::

            sage: tau = SimilarityClassType([[1, [1]], [1, [1]]])
            sage: tau.number_of_classes()
            1/2*q^2 - 1/2*q
        """
        if q is None:
            q = ZZ['q'].gen()
        if self.size() == 0:
            return q.parent().one()
        list_of_degrees = [PT.degree() for PT in self]
        maximum_degree = max(list_of_degrees)
        numerator = prod([prod([primitives(d+1, invertible=invertible, q=q)-i for i in range(list_of_degrees.count(d+1))]) for d in range(maximum_degree)])
        tau_list = list(self)
        D = dict((i, tau_list.count(i)) for i in tau_list)
        denominator = prod(factorial(D[primary_type]) for primary_type in D)
        return numerator / denominator

    def is_semisimple(self):
        """
        Return ``True`` if every primary similarity class type in ``self`` has
        all parts equal to ``1``.

        EXAMPLES::

            sage: tau = SimilarityClassType([[2, [1, 1]], [1, [1]]])
            sage: tau.is_semisimple()
            True
            sage: tau = SimilarityClassType([[2, [1, 1]], [1, [2]]])
            sage: tau.is_semisimple()
            False
        """
        return all(PT.partition().get_part(0) == 1 for PT in self)

    def is_regular(self):
        """
        Return ``True`` if every primary type in ``self`` has partition with one
        part.

        EXAMPLES::

            sage: tau = SimilarityClassType([[2, [1]], [1, [3]]])
            sage: tau.is_regular()
            True
            sage: tau = SimilarityClassType([[2, [1, 1]], [1, [3]]])
            sage: tau.is_regular()
            False
        """
        return all(len(PT.partition()) == 1 for PT in self)

    def rcf(self):
        """
        Return the partition corresponding to the rational canonical form of a
        matrix of type ``self``.

        EXAMPLES::

            sage: tau = SimilarityClassType([[2, [1, 1, 1]], [1, [3, 2]]])
            sage: tau.rcf()
            [5, 4, 2]
        """
        out_list = list()
        i = 0
        while True:
            new_part = sum([PT.partition().get_part(i)*PT.degree() for PT in self])
            if new_part:
                out_list.append(new_part)
            else:
                return Partition(out_list)
            i = i+1

    def class_card(self, q=None):
        """
        Return the number of matrices in each similarity class of type ``self``.

        INPUT:

        - ``q`` -- integer or an indeterminate

        EXAMPLES::

            sage: tau = SimilarityClassType([[1, [1, 1, 1, 1]]])
            sage: tau.class_card()
            1
            sage: tau = SimilarityClassType([[1, [1]], [1, [1]]])
            sage: tau.class_card()
            q^2 + q
        """
        if q is None:
            q = ZZ['q'].gen()
        return order_of_general_linear_group(self.size(), q=q) / self.centralizer_group_card(q=q)

    def number_of_matrices(self, invertible=False, q=None):
        """
        Return the number of matrices of type ``self``.

        INPUT:

        - ``invertible`` -- a boolean; return the number of invertible
          matrices if set

        EXAMPLES::

            sage: tau = SimilarityClassType([[1, [1]]])
            sage: tau.number_of_matrices()
            q
            sage: tau.number_of_matrices(invertible = True)
            q - 1
            sage: tau = SimilarityClassType([[1, [1]], [1, [1]]])
            sage: tau.number_of_matrices()
            1/2*q^4 - 1/2*q^2
        """
        if q is None:
            q = ZZ['q'].gen()
        return self.class_card(q=q)*self.number_of_classes(invertible=invertible, q=q)

    def statistic(self, func, q=None):
        r"""
        Return.

        .. MATH::

            \prod_{(d, \lambda)\in \tau} n_{\lambda}(q^d)

        where `n_{\lambda}(q)` is the value returned by ``func`` on the input
        `\lambda`.

        INPUT:

        - ``func`` -- a function that takes a partition to a polynomial in ``q``

        - ``q`` -- integer or an indeterminate

        EXAMPLES::

            sage: tau = SimilarityClassType([[1, [1]], [1, [2, 1]], [2, [1, 1]]])
            sage: from sage.combinat.similarity_class_type import fq
            sage: tau.statistic(lambda la: prod([fq(m) for m in la.to_exp()]))
            (q^9 - 3*q^8 + 2*q^7 + 2*q^6 - 4*q^5 + 4*q^4 - 2*q^3 - 2*q^2 + 3*q - 1)/q^9
            sage: q = ZZ['q'].gen()
            sage: tau.statistic(lambda la: q**la.size(), q = q)
            q^8
        """
        if q is None:
            q = FractionField(ZZ['q']).gen()
        return prod([PT.statistic(func, q=q) for PT in self])

    def invariant_subspace_generating_function(self, q=None, t=None):
        r"""
        Return the invariant subspace generating function of ``self``.

        The invariant subspace generating function is the function is the
        polynomial

        .. MATH::

            \sum_{j\geq 0} a_j(q) t^j,

        where `a_j(q)` denotes the number of `j`-dimensional invariant subspaces
        of dimensiona `j` for any matrix with the similarity class type ``self``
        with entries in a field of order `q`.

        EXAMPLES::

            sage: SimilarityClassType([[1, [2, 2]]]).invariant_subspace_generating_function()
            t^4 + (q + 1)*t^3 + (q^2 + q + 1)*t^2 + (q + 1)*t + 1
            sage: A = Matrix(GF(2),[(0, 1, 0, 0), (0, 1, 1, 1), (1, 0, 1, 0), (1, 1, 0, 0)])
            sage: SimilarityClassType(A).invariant_subspace_generating_function()
            t^4 + 1
        """
        if q is None:
            q = PolynomialRing(QQ, 'q').gen()
        S = q.parent()
        if t is None:
            t = PolynomialRing(S, 't').gen()
        return prod(p.invariant_subspace_generating_function(q=q, t=t) for p in self)


class SimilarityClassTypes(UniqueRepresentation, Parent):
    r"""
    Class of all similarity class types of size ``n`` with all primary matrix
    types greater than or equal to the primary matrix type ``min``.

    A similarity class type is a multiset of primary matrix types.

    INPUT:

    - ``n`` -- nonnegative integer
    - ``min`` -- a primary similarity class type

    EXAMPLES:

    If ``min`` is not specified, then the class of all matrix types of size
    ``n`` is constructed::

        sage: M = SimilarityClassTypes(2)
        sage: for tau in M:
        ....:     print(tau)
        [[1, [1]], [1, [1]]]
        [[1, [2]]]
        [[1, [1, 1]]]
        [[2, [1]]]

    If ``min`` is specified, then the class consists of only those similarity
    class types which are multisets of primary matrix types which either have
    size greater than that of ``min``, or if they have size equal to that of
    ``min``, then they occur after ``min`` in the iterator for
    ``PrimarySimilarityClassTypes(n)``, where ``n`` is the size of ``min``::

        sage: M = SimilarityClassTypes(2, min = [1, [1, 1]])
        sage: for tau in M:
        ....:     print(tau)
        [[1, [1, 1]]]
        [[2, [1]]]
    """
    @staticmethod
    def __classcall_private__(cls, n, min=None):
        r"""
        Create the class of similarity class types of size ``n`` consisting of
        primary similarity class types greater than or equal to ``min``.

        EXAMPLES::

            sage: M1 = SimilarityClassTypes(2, min = [1, [1]])
            sage: M2 = SimilarityClassTypes(2)
            sage: M1 is M2
            True
        """
        if min is None:
            min = PrimarySimilarityClassType(1, Partition([1]))
        if isinstance(min, list):
            min = PrimarySimilarityClassType(min[0], min[1])
        if not isinstance(min, PrimarySimilarityClassType):
            raise ValueError("min must be a PrimarySimilarityClassType")
        return super().__classcall__(cls, n, min)

    def __init__(self, n, min):
        r"""
        Initialize ``self``.

        TESTS::

            sage: M = SimilarityClassTypes(2)
            sage: TestSuite(M).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._n = n
        self._min = min

    def _element_constructor_(self, tau):
        """
        Construct an element of ``self``.

        INPUT:

        - ``tau`` -- list of primary similarity class types

        EXAMPLES::

            sage: M = SimilarityClassTypes(2)
            sage: elt = M([[1, [1]], [1, [1]]]); elt
            [[1, [1]], [1, [1]]]
            sage: elt.parent() is M
            True
        """
        ret = []
        for l in tau:
            if isinstance(l, PrimarySimilarityClassType):
                ret.append(l)
            else:
                ret.append(PrimarySimilarityClassType(*l))
        return self.element_class(self, ret)

    Element = SimilarityClassType

    def __iter__(self):
        r"""
        Iterator for vector partitions.

        EXAMPLES::

            sage: SimilarityClassTypes(3).cardinality()
            8

        A good test of the iterator is to see that all elements of
        `M_n(\GF{q})` or `GL_n(\GF{q})` are enumerated through
        types::

            sage: from sage.combinat.similarity_class_type import order_of_general_linear_group
            sage: q = QQ['q'].gen()
            sage: def test(n):
            ....:     M = SimilarityClassTypes(n)
            ....:     return M.sum(lambda la:1) == q**(n**2) and M.sum(lambda la:1, invertible = True)== order_of_general_linear_group(n)
            sage: all(test(n) for n in range(5))
            True
            sage: all(test(n) for n in range(5, 10)) # long time
            True
        """
        n = self._n
        min = self._min
        if n == 0:
            yield self.element_class(self, [])  # dimension zero has only empty type
        if min.size() > n:
            return
        else:
            # choose first part
            for PT in chain(PrimarySimilarityClassTypes(min.size(), min=min), *[PrimarySimilarityClassTypes(k) for k in range(min.size() + 1, n + 1)]):
                if PT.size() == n:
                    yield self.element_class(self, [PT])
                else:  # recursively find all possibilities for what remains of n
                    for smaller_type in SimilarityClassTypes(n - PT.size(), min=PT):
                        yield self.element_class(self, [PT] + list(smaller_type))

    def size(self):
        """
        Return size of ``self``.

        EXAMPLES::

            sage: tau = SimilarityClassType([[3, [3, 2, 1]], [2, [2, 1]]])
            sage: tau.parent().size()
            24
        """
        return self._n

    def sum(self, stat, sumover='matrices', invertible=False, q=None):
        r"""
        Return the sum of a local statistic over all types.

        Given a set of functions `n_{\lambda}(q)` (these could be polynomials or
        rational functions in `q`, for each similarity class type `\tau` define

        .. MATH::

            n_\tau(q) = \prod_{(d,\lambda)\in \tau} n_{\lambda}(q^d).

        This function returns

        .. MATH::

            \sum n_{\tau(g)}(q)

        where `\tau(g)` denotes the type of a matrix `g`, and the sum is over
        all `n \times n` matrices if ``sumover`` is set to ``'matrices'``, is
        over all `n \times n` similarity classes if ``sumover`` is set to
        ``'classes'``, and over all `n \times n` types if ``sumover`` is set
        to ``'types'``. If ``invertible`` is set to ``True``, then the sum is
        only over invertible matrices or classes.

        INPUT:

        - ``stat`` -- a function which takes partitions and returns a function
          of ``q``
        - ``sumover`` -- can be one of the following:

          * ``'matrices'``
          * ``'classes'``
          * ``'types'``

        - ``q`` -- integer or an indeterminate

        OUTPUT: a function of ``q``

        EXAMPLES::

            sage: M = SimilarityClassTypes(2)
            sage: M.sum(lambda la:1)
            q^4
            sage: M.sum(lambda la:1, invertible = True)
            q^4 - q^3 - q^2 + q
            sage: M.sum(lambda la:1, sumover = "classes")
            q^2 + q
            sage: M.sum(lambda la:1, sumover = "classes", invertible = True)
            q^2 - 1

        Burside's lemma can be used to calculate the number of similarity
        classes of matrices::

            sage: from sage.combinat.similarity_class_type import centralizer_algebra_dim, order_of_general_linear_group
            sage: q = ZZ['q'].gen()
            sage: M.sum(lambda la:q**centralizer_algebra_dim(la), invertible = True)/order_of_general_linear_group(2)
            q^2 + q
        """
        if sumover == "matrices":
            return sum([tau.statistic(stat, q=q)*tau.number_of_matrices(invertible=invertible, q=q) for tau in self])
        elif sumover == "classes":
            return sum([tau.statistic(stat, q=q)*tau.number_of_classes(invertible=invertible, q=q) for tau in self])
        elif sumover == "types":
            return sum([tau.statistic(stat, invertible=invertible, q=q) for tau in self])
        else:
            raise ValueError("invalid parameter %s" % (sumover))

################################################################################
#                 Similarity over rings of length two                          #
################################################################################


def dictionary_from_generator(gen):
    r"""
    Given a generator for a list of pairs `(c,f)`, construct a dictionary whose
    keys are the distinct values for `c` and whose value at `c` is the sum of
    `f` over all pairs of the form `(c',f)` such that `c=c'`.

    EXAMPLES::

        sage: from sage.combinat.similarity_class_type import dictionary_from_generator
        sage: dictionary_from_generator(((x // 2, x) for x in range(10)))
        {0: 1, 1: 5, 2: 9, 3: 13, 4: 17}

    It also works with lists::

        sage: dictionary_from_generator([(x // 2, x) for x in range(10)])
        {0: 1, 1: 5, 2: 9, 3: 13, 4: 17}

    .. NOTE::

        Since the generator is first converted to a list, memory usage could be
        high.
    """
    L = list(gen)
    setofkeys = list(set(item[0] for item in L))
    return dict((key, sum(entry[1] for entry in (pair for pair in L if pair[0] == key))) for key in setofkeys)


def matrix_similarity_classes(n, q=None, invertible=False):
    r"""
    Return the number of matrix similarity classes over a finite field of order
    ``q``.

    TESTS::

        sage: from sage.combinat.similarity_class_type import matrix_similarity_classes
        sage: matrix_similarity_classes(2)
        q^2 + q
        sage: matrix_similarity_classes(2, invertible = True)
        q^2 - 1
        sage: matrix_similarity_classes(2, invertible = True, q = 4)
        15
    """
    if q is None:
        q = ZZ['q'].gen()
    basering = q.parent()
    if n == 0:
        return basering.one()
    if invertible:
        tilde = 1 - ~q
        return sum(q**max(la) *
                   tilde ** len([x for x in la.to_exp() if x > 0])
                   for la in Partitions(n))
    return sum(q**max(la) for la in Partitions(n))


def matrix_centralizer_cardinalities(n, q=None, invertible=False):
    """
    Generate pairs consisting of centralizer cardinalities of matrices over a
    finite field and their frequencies.

    TESTS::

        sage: from sage.combinat.similarity_class_type import matrix_centralizer_cardinalities
        sage: list(matrix_centralizer_cardinalities(1))
        [(q - 1, q)]
        sage: list(matrix_centralizer_cardinalities(2))
        [(q^2 - 2*q + 1, 1/2*q^2 - 1/2*q),
        (q^2 - q, q),
        (q^4 - q^3 - q^2 + q, q),
        (q^2 - 1, 1/2*q^2 - 1/2*q)]
        sage: list(matrix_centralizer_cardinalities(2, invertible = True))
        [(q^2 - 2*q + 1, 1/2*q^2 - 3/2*q + 1),
        (q^2 - q, q - 1),
        (q^4 - q^3 - q^2 + q, q - 1),
        (q^2 - 1, 1/2*q^2 - 1/2*q)]
    """
    for tau in SimilarityClassTypes(n):
        yield (tau.centralizer_group_card(q=q), tau.number_of_classes(invertible=invertible, q=q))


def input_parsing(data):
    """
    Recognize and return the intended type of ``input``.

    TESTS::

        sage: from sage.combinat.similarity_class_type import input_parsing
        sage: input_parsing(Partition([2, 1]))
        ('par', [2, 1])
        sage: input_parsing(PrimarySimilarityClassType(2, [2, 1]))
        ('pri', [2, [2, 1]])
        sage: input_parsing(SimilarityClassType([[2, [2, 1]]]))
        ('sim', [[2, [2, 1]]])
        sage: input_parsing([2, 1])
        ('par', [2, 1])
        sage: input_parsing([2, [2, 1]])
        ('pri', [2, [2, 1]])
        sage: input_parsing([[2, [2, 1]]])
        ('sim', [[2, [2, 1]]])
    """
    if isinstance(data, SimilarityClassType):
        case = 'sim'
    elif isinstance(data, PrimarySimilarityClassType):
        case = 'pri'
    elif isinstance(data, Partition):
        case = 'par'
    else:
        try:
            data = Partition(data)
            case = 'par'
        except (TypeError, ValueError):
            try:
                data = SimilarityClassType(data)
                case = 'sim'
            except (TypeError, ValueError):
                try:
                    data = PrimarySimilarityClassType(*data)
                    case = 'pri'
                except (TypeError, ValueError):
                    raise ValueError("expected a Partition, a SimilarityClassType or a PrimarySimilarityClassType, got a %s" % type(data))
    return case, data


def ext_orbits(input_data, q=None, selftranspose=False):
    r"""
    Return the number of orbits in `\mathrm{Ext}^1(M, M)` for the action of
    `\mathrm{Aut}(M, M)`, where `M` is the `\GF{q[t]}`-module constructed
    from ``input_data``.

    INPUT:

    - ``input_data`` -- input for :func:`input_parsing()`
    - ``q`` -- (default: `q`) an integer or an indeterminate
    - ``selftranspose`` -- boolean (default: ``False``); stating if we only
      want selftranspose type

    TESTS::

        sage: from sage.combinat.similarity_class_type import ext_orbits
        sage: ext_orbits([6, 1])
        q^7 + q^6 + q^5
        sage: ext_orbits([6, 1], selftranspose = True)
        q^7 + q^6 - q^5
        sage: ext_orbits([6, 1, 1])
        q^8 + 2*q^7 + 2*q^6 + 2*q^5
        sage: ext_orbits ([6, 1, 1], selftranspose = True)
        q^8 + 2*q^7
        sage: ext_orbits([2, 2])
        q^4 + q^3 + q^2
        sage: ext_orbits([2, 2], selftranspose = True)
        q^4 + q^3 + q^2
        sage: ext_orbits([2, 2, 2])
        q^6 + q^5 + 2*q^4 + q^3 + 2*q^2
        sage: ext_orbits([2, 2, 2], selftranspose = True)
        q^6 + q^5 + 2*q^4 + q^3
        sage: ext_orbits([2, 2, 2, 2])
        q^8 + q^7 + 3*q^6 + 3*q^5 + 5*q^4 + 3*q^3 + 3*q^2
        sage: ext_orbits([2, 2, 2, 2], selftranspose = True)
        q^8 + q^7 + 3*q^6 + 3*q^5 + 3*q^4 + q^3 + q^2
        sage: ext_orbits([2, [6, 1]])
        q^14 + q^12 + q^10
        sage: ext_orbits([[2, [6, 1]]])
        q^14 + q^12 + q^10
    """
    # Comments cite items in the paper "Similarity over rings of length two" by
    # Prasad, Singla, and Spallone.
    if q is None:
        q = FractionField(QQ['q']).gen()
    case, data = input_parsing(input_data)
    if case == 'par':
        la = data
        if la.size() == 0:
            return q.parent()(1)
        if max(la) == 1:
            return matrix_similarity_classes(len(la), q=q)
        elif len(la) == 1:
            return q**la.size()
        elif len(la) == 2 and list(la).count(1) == 1:  # see Table 3
            m = max(la) - 1
            if selftranspose:
                return q**(m + 2) + q**(m + 1) - q**m
            else:
                return q**(m + 2) + q**(m + 1) + q**m
        elif len(la) == 3 and list(la).count(1) == 2:  # see Table 4
            m = max(la) - 1
            if not selftranspose:
                return q**m*(q**3 + 2*q**2 + 2*q + 2)
            else:
                return q**m*(q**3 + 2*q**2)
        elif min(la) == 2 and max(la) == 2:
            return matrix_similarity_classes_length_two(len(la), q=q, selftranspose=selftranspose)
        else:
            raise ValueError('partition %s not implemented for ExtOrbitClasses.orbits' % (la))
    elif case == 'pri':
        tau = data
        return ext_orbits(tau.partition(), q=q, selftranspose=selftranspose).substitute(q=q**tau.degree())
    elif case == 'sim':
        tau = data
        return prod([ext_orbits(PT, q=q, selftranspose=selftranspose) for PT in tau])


def matrix_similarity_classes_length_two(n, q=None, selftranspose=False, invertible=False):
    """
    Return the number of similarity classes of matrices of order ``n`` with
    entries in a principal ideal local ring of length two.

    INPUT:

    - ``n`` -- the order
    - ``q`` -- (default: `q`) an integer or an indeterminate
    - ``selftranspose`` -- boolean (default: ``False``); stating if we only want
      selftranspose type
    - ``invertible`` -- boolean (default: ``False``); stating if we only want
      invertible type

    EXAMPLES:

    We can generate Table 6 of [PSS13]_::

        sage: from sage.combinat.similarity_class_type import matrix_similarity_classes_length_two
        sage: matrix_similarity_classes_length_two(2)
        q^4 + q^3 + q^2
        sage: matrix_similarity_classes_length_two(2, invertible = True)
        q^4 - q
        sage: matrix_similarity_classes_length_two(3)
        q^6 + q^5 + 2*q^4 + q^3 + 2*q^2
        sage: matrix_similarity_classes_length_two(3, invertible = true)
        q^6 - q^3 + 2*q^2 - 2*q
        sage: matrix_similarity_classes_length_two(4)
        q^8 + q^7 + 3*q^6 + 3*q^5 + 5*q^4 + 3*q^3 + 3*q^2
        sage: matrix_similarity_classes_length_two(4, invertible = True)
        q^8 + q^6 - q^5 + 2*q^4 - 2*q^3 + 2*q^2 - 3*q

    And also Table 7::

        sage: matrix_similarity_classes_length_two(2, selftranspose = True)
        q^4 + q^3 + q^2
        sage: matrix_similarity_classes_length_two(2, selftranspose = True, invertible = True)
        q^4 - q
        sage: matrix_similarity_classes_length_two(3, selftranspose = True)
        q^6 + q^5 + 2*q^4 + q^3
        sage: matrix_similarity_classes_length_two(3, selftranspose = True, invertible = True)
        q^6 - q^3
        sage: matrix_similarity_classes_length_two(4, selftranspose = True)
        q^8 + q^7 + 3*q^6 + 3*q^5 + 3*q^4 + q^3 + q^2
        sage: matrix_similarity_classes_length_two(4, selftranspose = True, invertible = True)
        q^8 + q^6 - q^5 - q
    """
    if q is None:
        q = FractionField(QQ['q']).gen()
    return sum([tau.number_of_classes(invertible=invertible, q=q)*ext_orbits(tau, q=q, selftranspose=selftranspose) for tau in SimilarityClassTypes(n)])


def ext_orbit_centralizers(input_data, q=None, selftranspose=False):
    r"""
    Generate pairs consisting of centralizer cardinalities of orbits in
    `\mathrm{Ext}^1(M, M)` for the action of `\mathrm{Aut}(M, M)`, where `M` is
    the `\GF{q[t]}`-module constructed from ``input`` and their frequencies.

    INPUT:

    - ``input_data`` -- input for :func:`input_parsing()`
    - ``q`` -- (default: `q`) an integer or an indeterminate
    - ``selftranspose`` -- boolean (default: ``False``); stating if we only want
      selftranspose type

    TESTS::

        sage: from sage.combinat.similarity_class_type import ext_orbit_centralizers
        sage: list(ext_orbit_centralizers([6, 1]))
        [(q^9 - 2*q^8 + q^7, q^6),
         (q^7 - 2*q^6 + q^5, q^7 - q^6),
         (q^7 - q^6, q^6 + q^5)]
        sage: list(ext_orbit_centralizers([6, 1], selftranspose = True))
        [(q^9 - 2*q^8 + q^7, q^6),
         (q^7 - 2*q^6 + q^5, q^7 - q^6),
         (q^7 - q^6, q^6 - q^5)]
        sage: list(ext_orbit_centralizers([6, 1, 1]))
        [(q^12 - 3*q^11 + 3*q^10 - q^9, 1/2*q^7 - 1/2*q^6),
         (q^8 - 3*q^7 + 3*q^6 - q^5, 1/2*q^8 - q^7 + 1/2*q^6),
         (q^12 - 2*q^11 + q^10, q^6),
         (q^8 - 2*q^7 + q^6, q^7 - q^6),
         (q^14 - 2*q^13 + 2*q^11 - q^10, q^6),
         (q^10 - 2*q^9 + 2*q^7 - q^6, q^7 - q^6),
         (q^12 - q^11 - q^10 + q^9, 1/2*q^7 - 1/2*q^6),
         (q^8 - q^7 - q^6 + q^5, 1/2*q^8 - q^7 + 1/2*q^6),
         (q^8 - 2*q^7 + q^6, q^7 - q^6),
         (q^8 - q^7, q^6 + 2*q^5),
         (q^10 - 2*q^9 + q^8, 2*q^6)]
        sage: list(ext_orbit_centralizers([6, 1, 1], selftranspose = True))
        [(q^12 - 3*q^11 + 3*q^10 - q^9, 1/2*q^7 - 1/2*q^6),
         (q^8 - 3*q^7 + 3*q^6 - q^5, 1/2*q^8 - q^7 + 1/2*q^6),
         (q^12 - 2*q^11 + q^10, q^6),
         (q^8 - 2*q^7 + q^6, q^7 - q^6),
         (q^14 - 2*q^13 + 2*q^11 - q^10, q^6),
         (q^10 - 2*q^9 + 2*q^7 - q^6, q^7 - q^6),
         (q^12 - q^11 - q^10 + q^9, 1/2*q^7 - 1/2*q^6),
         (q^8 - q^7 - q^6 + q^5, 1/2*q^8 - q^7 + 1/2*q^6),
         (q^8 - 2*q^7 + q^6, q^7 - q^6),
         (q^8 - q^7, q^6)]
        sage: list(ext_orbit_centralizers([2, [6, 1, 1]], selftranspose = True))
        [(q^24 - 3*q^22 + 3*q^20 - q^18, 1/2*q^14 - 1/2*q^12),
         (q^16 - 3*q^14 + 3*q^12 - q^10, 1/2*q^16 - q^14 + 1/2*q^12),
         (q^24 - 2*q^22 + q^20, q^12),
         (q^16 - 2*q^14 + q^12, q^14 - q^12),
         (q^28 - 2*q^26 + 2*q^22 - q^20, q^12),
         (q^20 - 2*q^18 + 2*q^14 - q^12, q^14 - q^12),
         (q^24 - q^22 - q^20 + q^18, 1/2*q^14 - 1/2*q^12),
         (q^16 - q^14 - q^12 + q^10, 1/2*q^16 - q^14 + 1/2*q^12),
         (q^16 - 2*q^14 + q^12, q^14 - q^12),
         (q^16 - q^14, q^12)]
        sage: list(ext_orbit_centralizers([[2, [6, 1, 1]]], selftranspose = True))
        [(q^24 - 3*q^22 + 3*q^20 - q^18, 1/2*q^14 - 1/2*q^12),
         (q^16 - 3*q^14 + 3*q^12 - q^10, 1/2*q^16 - q^14 + 1/2*q^12),
         (q^24 - 2*q^22 + q^20, q^12),
         (q^16 - 2*q^14 + q^12, q^14 - q^12),
         (q^28 - 2*q^26 + 2*q^22 - q^20, q^12),
         (q^20 - 2*q^18 + 2*q^14 - q^12, q^14 - q^12),
         (q^24 - q^22 - q^20 + q^18, 1/2*q^14 - 1/2*q^12),
         (q^16 - q^14 - q^12 + q^10, 1/2*q^16 - q^14 + 1/2*q^12),
         (q^16 - 2*q^14 + q^12, q^14 - q^12),
         (q^16 - q^14, q^12)]
    """
    # Comments cite items in the paper "Similarity over rings of length two" by
    # Prasad, Singla, and Spallone.
    if q is None:
        q = FractionField(QQ['q']).gen()
    case, data = input_parsing(input_data)
    if case == 'par':
        la = data
        if len(la) == 0:
            yield (1, 1)
            return
        elif max(la) == 1:
            for item in matrix_centralizer_cardinalities(len(la), q=q):
                yield item
            return
        elif len(la) == 1:
            yield (q**la[0] - q**(la[0]-1), q**la[0])
            return
        elif len(la) == 2 and list(la).count(1) == 1:  # see Table 3
            m = max(la) - 1
            yield (q**(m + 4) - 2*q**(m + 3) + q**(m + 2), q**(m + 1))  # (8.5.1)
            yield (q**(m + 2) - 2*q**(m + 1) + q**m, q**(m + 2) - q**(m + 1))  # (8.5.2)
            if selftranspose:
                yield (q**(m + 2) - q**(m + 1), q**(m+1) - q**m)  # (8.5.3) and (8.5.4)
            else:
                yield (q**(m + 2) - q**(m + 1), q**(m + 1) + q**m)  # (8.5.3) and (8.5.4)
            return
        elif len(la) == 3 and list(la).count(1) == 2:  # see Table 4
            m = max(la) - 1
            for item in matrix_centralizer_cardinalities(2, q=q):
                yield (item[0]*(q**(m + 5) - q**(m + 4)), item[1]*q**m)  # (8.6.1)
                yield (item[0]*(q**(m + 1) - q**m), item[1]*(q**(m + 1) - q**m))  # (8.6.2)
            yield (q**(m + 3) - 2*q**(m + 2) + q**(m+1), q**(m + 2) - q**(m + 1))  # (8.6.3)
            if selftranspose:
                yield (q**(m + 3) - q**(m+2), q**(m+1))  # (8.6.4), (8.6.5) and (8.6.7)
            else:
                yield (q**(m + 3) - q**(m+2), q**(m + 1) + 2*q**m)  # (8.6.4), (8.6.5) and (8.6.7)
                yield (q**(m + 5) - 2*q**(m + 4) + q**(m + 3), 2*q**(m + 1))  # (8.6.6) and (8.6.8)
            return
        elif max(la) == 2 and min(la) == 2:
            for item in matrix_centralizer_cardinalities_length_two(len(la), q=q, selftranspose=selftranspose):
                yield item
        else:
            raise ValueError('partition %s not implemented for ExtOrbitClasses.orbit_centralizers' % (la))
    elif case == 'pri':
        tau = data
        for item in ext_orbit_centralizers(tau.partition(), selftranspose=selftranspose):
            yield (item[0].substitute(q=q**tau.degree()), item[1].substitute(q=q**tau.degree()))
    elif case == 'sim':
        tau = data
        for item in product(*[IterableFunctionCall(lambda x: ext_orbit_centralizers(x, q=q, selftranspose=selftranspose), PT) for PT in tau]):
            size = prod([list(entry)[0] for entry in item])
            freq = prod([list(entry)[1] for entry in item])
            yield (size, freq)


def matrix_centralizer_cardinalities_length_two(n, q=None, selftranspose=False, invertible=False):
    r"""
    Generate pairs consisting of centralizer cardinalities of matrices over a
    principal ideal local ring of length two with residue field of order ``q``
    and their frequencies.

    INPUT:

    - ``n`` -- the order
    - ``q`` -- (default: `q`) an integer or an indeterminate
    - ``selftranspose`` -- boolean (default: ``False``); stating if we only want
      selftranspose type
    - ``invertible`` -- boolean (default: ``False``); stating if we only want
      invertible type

    TESTS::

        sage: from sage.combinat.similarity_class_type import matrix_centralizer_cardinalities_length_two
        sage: list(matrix_centralizer_cardinalities_length_two(1))
        [(q^2 - q, q^2)]
        sage: list(matrix_centralizer_cardinalities_length_two(2))
        [(q^4 - 2*q^3 + q^2, 1/2*q^4 - 1/2*q^3),
        (q^4 - q^3, q^3),
        (q^6 - 2*q^5 + q^4, 1/2*q^3 - 1/2*q^2),
        (q^6 - q^5, q^2),
        (q^8 - q^7 - q^6 + q^5, q^2),
        (q^6 - q^4, 1/2*q^3 - 1/2*q^2),
        (q^4 - q^2, 1/2*q^4 - 1/2*q^3)]
        sage: from sage.combinat.similarity_class_type import dictionary_from_generator
        sage: dictionary_from_generator(matrix_centralizer_cardinalities_length_two(2, q = 2))
        {4: 4, 8: 8, 12: 4, 16: 2, 32: 4, 48: 2, 96: 4}
    """
    if q is None:
        q = FractionField(QQ['q']).gen()
    for tau in SimilarityClassTypes(n):
        for pair in ext_orbit_centralizers(tau, q=q, selftranspose=selftranspose):
            yield (q**tau.centralizer_algebra_dim()*pair[0], tau.number_of_classes(invertible=invertible, q=q)*pair[1])
