r"""
Symmetric functions of `r`-coloured cycles

This module provides the algebra `\Lambda(r)` of `r`-coloured symmetric
functions in the sense of Henderson [Henderson2004]_.  The case `r = 2`
is the *hyperoctahedral* case.

The corresponding species live in
:mod:`sage.rings.species_hyperoctahedral`; the combinatorial model of
the canonical permutation representation of the wreath product
`W(r,n) = C_r \wr S_n` lives in
:mod:`sage.groups.perm_gps.hyperoctahedral_group`.

Following Henderson, `\Lambda(r)` is the polynomial algebra

.. MATH::

    \Lambda(r) = \mathbb{C}[p_i(\zeta) \mid i \in \mathbb{Z}^+, \zeta \in \mu_r]

graded by `\deg(p_i(\zeta)) = i`.  Rather than representing `\zeta` by
an actual root of unity, we label the `r` possible values of `\zeta` by
the integers `0, 1, \ldots, r-1`.  The algebra `\Lambda(r)` is thus
identified with the tensor product `\Lambda(1)^{\otimes r}` of `r`
copies of the ordinary symmetric function algebra, the `j`-th tensor
factor carrying the power sums `p_i(\zeta^j)`.

Macdonald [Macdonald1995]_ (Appendix B) treats the same algebra with
the colours read either as conjugacy classes `\zeta^j` of the cyclic
group `C_r` or as its irreducible characters `\chi^j`, the two
conventions being related by the discrete Fourier transform (his
(7.1), (7.1')):

.. MATH::

    p_i(\chi^k) = \frac{1}{r} \sum_{j=0}^{r-1} \zeta^{jk}\, p_i(\zeta^j),
    \qquad
    p_i(\zeta^j) = \sum_{k=0}^{r-1} \zeta^{-jk}\, p_i(\chi^k).

Every basis `B` of the ordinary symmetric functions and either colour
convention therefore provides a basis of `\Lambda(r)`, indexed by
partition tuples of level `r`:  the basis element
`B_{\lambda^{(0)}, \ldots, \lambda^{(r-1)}}` is the product over colours
`j` of `B_{\lambda^{(j)}}` with each power sum `p_i` replaced by
`p_i(\zeta^j)` or `p_i(\chi^j)`, respectively.

EXAMPLES::

    sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
    sage: H = HyperoctahedralSymmetricFunctions(2)
    sage: S = SymmetricFunctions(QQ); p = S.p(); s = S.s()

The power sums indexed by classes `\zeta^j` are the algebra generators::

    sage: P = H(p)
    sage: P[[2, 1], [3]]
    p_{2,1}(ζ^0)*p_3(ζ^1)
    sage: P[[2], []] * P[[], [2, 1]]
    p_2(ζ^0)*p_{2,1}(ζ^1)

Every other classical basis gives a basis with the classical structure
constants, applied colourwise::

    sage: Y = H(s)
    sage: Y[[2], []] * Y[[2], []]
    s_{2,2}(ζ^0) + s_{3,1}(ζ^0) + s_4(ζ^0)

The colour may equally label the irreducible characters `\chi^j` of
`C_r` (Macdonald, Appendix B, (7.1)); over `\QQ` this is available for
`r \leq 2`::

    sage: Q = H(p, 'chi')
    sage: Q[[], [2]] == (P[[2], []] - P[[], [2]]) / 2
    True
    sage: Z = H(s, 'chi')
    sage: Z[[], [2]]
    s_2(χ^1)

All conversions between these bases go through the class-indexed power
sums and the discrete Fourier transform::

    sage: P(Z[[], [2]])
    1/8*p_{1,1}(ζ^1) - 1/4*p_2(ζ^1) - 1/4*p_1(ζ^0)*p_1(ζ^1) + 1/8*p_{1,1}(ζ^0) + 1/4*p_2(ζ^0)
    sage: Z(P(Z[[], [2]])) == Z[[], [2]]
    True
    sage: Y(Z[[], [2]])
    3/8*s_{1,1}(ζ^1) - 1/8*s_2(ζ^1) - 1/4*s_1(ζ^0)*s_1(ζ^1) - 1/8*s_{1,1}(ζ^0) + 3/8*s_2(ζ^0)

The two plethysms of Henderson are available through the public method
:meth:`~sage.rings.sf_hyperoctahedral.HyperoctahedralSymmetricFunctions.Element.plethysm`::

    sage: p[2].plethysm(P[[3], []])
    p_6(ζ^0)
    sage: P[[2], [1]].plethysm(p[3])
    p_6(ζ^0)*p_3(ζ^1)

The two generator rules are (with `p_i` an ordinary power sum and
`p_j(\zeta^a)` a generator of `\Lambda(r)`)::

    sage: H = HyperoctahedralSymmetricFunctions(3)
    sage: p = SymmetricFunctions(QQ).p()
    sage: p[2].plethysm(H(p)[[], [3], []])     # p_2 o p_3(ζ^1) = p_6(ζ^1)
    p_6(ζ^1)
    sage: H(p)[[2], [], []].plethysm(p[3])     # p_2(ζ^0) o p_3 = p_6(ζ^0)
    p_6(ζ^0)
    sage: H(p)[[], [2], []].plethysm(p[3])     # p_2(ζ^1) o p_3 = p_6(ζ^3) = p_6(ζ^0)
    p_6(ζ^0)
    sage: H(p)[[], [2], []].plethysm(p[2])     # p_2(ζ^1) o p_2 = p_4(ζ^2)
    p_4(ζ^2)

Both plethysms are linear in the left argument, and the second one is
additive in the right argument for a single generator::

    sage: f = p[2] + 3*p[1, 1]
    sage: g = p[3] - p[1]
    sage: (f + g).plethysm(P[[3], []]) == f.plethysm(P[[3], []]) + g.plethysm(P[[3], []])
    True
    sage: P[[2], []].plethysm(f + g) == P[[2], []].plethysm(f) + P[[2], []].plethysm(g)
    True
    sage: F = P[[], [2]] + P[[1], [1]]*P[[1], []]
    sage: (F + P[[], [3]]).plethysm(f) == F.plethysm(f) + P[[], [3]].plethysm(f)
    True

The plethysms are associative whenever both sides are defined::

    sage: h = P[[], [2]] + P[[1], [1]]*P[[1], []]
    sage: f.plethysm(g).plethysm(h) == f.plethysm(g.plethysm(h))
    True
    sage: F.plethysm(g).plethysm(f) == F.plethysm(g.plethysm(f))
    True

REFERENCES:

.. [Henderson2004] Anthony Henderson.
   *Representations of wreath products on cohomology of De Concini-Procesi
   compactifications*.
   International Mathematics Research Notices 2004, no. 20, 981-1021.

.. [Macdonald1995] I. G. Macdonald.
   *Symmetric functions and Hall polynomials*. 2nd ed., Oxford University
   Press, 1995, Appendix B.

AUTHORS:

- Martin Rubey (2026): initial version
"""


# ****************************************************************************
#       Copyright (C) 2026 Martin Rubey <martin.rubey@tuwien.ac.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import product as iproduct

from sage.arith.misc import multinomial
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.tensor import tensor
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.partition import Partition
from sage.combinat.partition_tuple import PartitionTuple, PartitionTuples_level
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.sf.sfa import SymmetricFunctionAlgebra_generic
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.structure.sage_object import SageObject

# aliases for the two colour conventions
_INDEXING_ALIASES = {'zeta': 'zeta', 'classes': 'zeta',
                     'chi': 'chi', 'characters': 'chi'}


class _HyperoctahedralSymmetricFunctionsFactory(SageObject):
    r"""
    Factory for the bases of the algebra `\Lambda(r)` of `r`-coloured
    symmetric functions.

    Created by calling :class:`HyperoctahedralSymmetricFunctions` with
    only the level `r`.  See the documentation there for examples.
    """
    def __init__(self, r):
        r"""
        TESTS::

            sage: from sage.rings.sf_hyperoctahedral import _HyperoctahedralSymmetricFunctionsFactory
            sage: _HyperoctahedralSymmetricFunctionsFactory(2)
            Hyperoctahedral symmetric functions of level 2
        """
        self._r = r

    def level(self):
        r"""
        Return the level `r` of ``self``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: HyperoctahedralSymmetricFunctions(3).level()
            3
        """
        return self._r

    def __call__(self, classical, indexing='zeta'):
        r"""
        Return the basis of `\Lambda(r)` built from the classical basis
        ``classical`` and the colour convention ``indexing``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: H(S.p(), 'chi')
            Hyperoctahedral symmetric functions of level 2 over Rational Field in the powersum basis with χ-colours
        """
        return HyperoctahedralSymmetricFunctions(self._r, classical, indexing)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: HyperoctahedralSymmetricFunctions(3)
            Hyperoctahedral symmetric functions of level 3
        """
        return f"Hyperoctahedral symmetric functions of level {self._r}"


class HyperoctahedralSymmetricFunctions(CombinatorialFreeModule):
    r"""
    A basis of the algebra `\Lambda(r)` of `r`-coloured symmetric
    functions.

    INPUT:

    - ``r`` -- positive integer; the order of the cyclic group
    - ``classical`` -- a basis of the ordinary symmetric functions, that
      is, a realisation of
      :class:`~sage.combinat.sf.sf.SymmetricFunctions`; its base ring
      is the base ring of the result
    - ``indexing`` -- ``'zeta'`` (default; aliases: ``'classes'``) to
      read the colour `j` as the conjugacy class `\zeta^j` of `C_r`, or
      ``'chi'`` (alias: ``'characters'``) to read it as the irreducible
      character `\chi^j` of `C_r`.  Character indexing requires a
      primitive `r`-th root of unity and `1/r` in the base ring, so over
      `\QQ` it is available only for `r \leq 2`; use
      ``CyclotomicField(r)`` or ``UniversalCyclotomicField()`` for
      larger `r`.

    Calling ``HyperoctahedralSymmetricFunctions(r)`` with the level
    alone returns a factory whose call takes ``classical`` and
    ``indexing``.

    The basis is indexed by :class:`~sage.combinat.partition_tuple.PartitionTuple`\ s
    of level `r`, the partition tuple `(\lambda^{(0)}, \ldots, \lambda^{(r-1)})`
    corresponding to the product over colours `j` of the classical
    basis element `B_{\lambda^{(j)}}` in the given colour convention.
    Multiplication is componentwise classical multiplication, and
    conversions between any two such bases (with the same level and
    base ring) go through the class-indexed power sums and the discrete
    Fourier transform (Macdonald, Appendix B, (7.1), (7.1')).

    EXAMPLES::

        sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
        sage: H = HyperoctahedralSymmetricFunctions(2)
        sage: S = SymmetricFunctions(QQ); s = S.s()
        sage: P = H(S.p()); Y = H(s, 'chi')
        sage: P
        Hyperoctahedral symmetric functions of level 2 over Rational Field in the powersum basis with ζ-colours
        sage: Y
        Hyperoctahedral symmetric functions of level 2 over Rational Field in the Schur basis with χ-colours
        sage: tau = P[[2, 1], [3]]
        sage: tau
        p_{2,1}(ζ^0)*p_3(ζ^1)
        sage: tau.degree()
        6

    The product of Schur basis elements is the classical Littlewood-Richardson
    product, applied colourwise::

        sage: Y[[2], [1]] * Y[[2], [1]]
        s_{2,2}(χ^0)*s_{1,1}(χ^1) + s_{2,2}(χ^0)*s_2(χ^1)
         + s_{3,1}(χ^0)*s_{1,1}(χ^1) + s_{3,1}(χ^0)*s_2(χ^1)
         + s_4(χ^0)*s_{1,1}(χ^1) + s_4(χ^0)*s_2(χ^1)

    Conversions between any two of these bases are registered with the
    coercion model, so mixed arithmetic works as well::

        sage: P(Y[[2], []])
        1/8*p_{1,1}(ζ^1) + 1/4*p_2(ζ^1) + 1/4*p_1(ζ^0)*p_1(ζ^1) + 1/8*p_{1,1}(ζ^0) + 1/4*p_2(ζ^0)
        sage: Y(P(Y[[2], [1]])) == Y[[2], [1]]
        True
        sage: P[[2], []] + Y[[1], []]
        1/2*p_1(ζ^1) + 1/2*p_1(ζ^0) + p_2(ζ^0)

    TESTS::

        sage: TestSuite(P).run()
        sage: TestSuite(Y).run()
    """
    @staticmethod
    def __classcall__(cls, r, classical=None, indexing='zeta'):
        r"""
        Normalize the arguments for unique representation.

        TESTS::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: S = SymmetricFunctions(QQ)
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: HyperoctahedralSymmetricFunctions(2, S.p()) is H(S.powersum())
            True
            sage: HyperoctahedralSymmetricFunctions(2, S.p(), 'classes') is H(S.p(), 'zeta')
            True
            sage: HyperoctahedralSymmetricFunctions(2, S.s(), 'characters') is H(S.schur(), 'chi')
            True
            sage: HyperoctahedralSymmetricFunctions(0)
            Traceback (most recent call last):
            ...
            ValueError: r must be a positive integer
            sage: HyperoctahedralSymmetricFunctions(QQ, 2)
            Traceback (most recent call last):
            ...
            ValueError: r must be a positive integer
            sage: HyperoctahedralSymmetricFunctions(2, S.p(), 'colors')
            Traceback (most recent call last):
            ...
            ValueError: indexing must be one of 'zeta', 'classes', 'chi', 'characters'
            sage: HyperoctahedralSymmetricFunctions(2, 'schur')
            Traceback (most recent call last):
            ...
            ValueError: classical must be a basis of ordinary symmetric functions, such as SymmetricFunctions(QQ).schur()

        Note that the previous signature
        ``HyperoctahedralSymmetricFunctions(base_ring, r)`` of this
        class has been replaced by the current one.
        """
        try:
            r = ZZ(r)
        except (TypeError, ValueError) as exc:
            raise ValueError("r must be a positive integer") from exc
        if r < 1:
            raise ValueError("r must be a positive integer")
        if classical is None:
            return _HyperoctahedralSymmetricFunctionsFactory(r)
        if not isinstance(classical, SymmetricFunctionAlgebra_generic):
            raise ValueError("classical must be a basis of ordinary symmetric "
                             "functions, such as SymmetricFunctions(QQ).schur()")
        try:
            indexing = _INDEXING_ALIASES[indexing]
        except (KeyError, TypeError):
            raise ValueError("indexing must be one of 'zeta', 'classes', "
                             "'chi', 'characters'") from None
        return super().__classcall__(cls, r, classical, indexing)

    def __init__(self, r, classical, indexing):
        r"""
        Initialize the basis of `\Lambda(r)`.

        TESTS::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: S = SymmetricFunctions(QQ)
            sage: HyperoctahedralSymmetricFunctions(3, S.p(), 'chi')
            Traceback (most recent call last):
            ...
            ValueError: character colours require a primitive root of unity of order 3 and 1/3 in the base ring; use CyclotomicField(3) or UniversalCyclotomicField()
        """
        from sage.rings.rational_field import QQ
        self._r = ZZ(r)
        self._classical = classical
        self._indexing = indexing
        R = classical.base_ring()
        if indexing == 'chi':
            if r > 1:
                try:
                    self._zeta = R.zeta(r)
                    R(QQ((1, r)))
                except (AttributeError, ValueError, TypeError,
                        ArithmeticError, ZeroDivisionError) as exc:
                    raise ValueError("character colours require a primitive "
                                     f"root of unity of order {r} and 1/{r} "
                                     f"in the base ring; use CyclotomicField({r}) "
                                     "or UniversalCyclotomicField()") from exc
            else:
                self._zeta = R.one()
        category = GradedAlgebrasWithBasis(R).Commutative()
        # PartitionTuples(1) collapses to Partitions; use the level-1
        # parent directly so that keys are always partition tuples.
        CombinatorialFreeModule.__init__(self, R,
                                         basis_keys=PartitionTuples_level(self._r),
                                         category=category,
                                         prefix=classical.prefix(),
                                         bracket=False)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: H(S.p())
            Hyperoctahedral symmetric functions of level 2 over Rational Field in the powersum basis with ζ-colours
            sage: H(S.schur(), 'chi')
            Hyperoctahedral symmetric functions of level 2 over Rational Field in the Schur basis with χ-colours
        """
        name = str(self._classical)
        if ' in the ' in name:
            name = name.rsplit(' in the ', 1)[1]
        else:
            name = self._classical.prefix() + " basis"
        glyph = 'χ' if self._indexing == 'chi' else 'ζ'
        return (f"Hyperoctahedral symmetric functions of level {self._r} "
                f"over {self.base_ring()} in the {name} with {glyph}-colours")

    def level(self):
        r"""
        Return the level `r` of ``self``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: HyperoctahedralSymmetricFunctions(2)(SymmetricFunctions(QQ).p()).level()
            2
        """
        return self._r

    def classical_basis(self):
        r"""
        Return the classical basis of ordinary symmetric functions used
        to construct ``self``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: H(S.m()).classical_basis()
            Symmetric Functions over Rational Field in the monomial basis
        """
        return self._classical

    def indexing(self):
        r"""
        Return the colour convention of ``self``: ``'zeta'`` (classes) or
        ``'chi'`` (characters).

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: H(SymmetricFunctions(QQ).p(), 'characters').indexing()
            'chi'
        """
        return self._indexing

    def __getitem__(self, x):
        r"""
        Return the basis element indexed by ``x``.

        The index may be given as a
        :class:`~sage.combinat.partition_tuple.PartitionTuple`, as a
        tuple or list of partitions (each partition itself given as a
        list, tuple or :class:`~sage.combinat.partition.Partition`), or,
        when the level is `1`, as a single partition or an integer.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: P[[2, 1], [3]]
            p_{2,1}(ζ^0)*p_3(ζ^1)
            sage: P[([2, 1], [3])]
            p_{2,1}(ζ^0)*p_3(ζ^1)
            sage: P[PartitionTuple([[2, 1], [3]])]
            p_{2,1}(ζ^0)*p_3(ζ^1)

        For level one, a single partition or an integer is accepted::

            sage: P1 = HyperoctahedralSymmetricFunctions(1)(SymmetricFunctions(QQ).s())
            sage: P1[[2]]
            s_2(ζ^0)
            sage: P1[2]
            s_2(ζ^0)

        TESTS::

            sage: P[[2, 1]]
            Traceback (most recent call last):
            ...
            ValueError: expected a tuple of 2 partitions
            sage: P[2]
            Traceback (most recent call last):
            ...
            ValueError: a single partition can only be used for level 1
            sage: P[PartitionTuple([[2], [1], []])]
            Traceback (most recent call last):
            ...
            ValueError: expected a partition tuple of level 2, got level 3
        """
        C = self._indices
        if isinstance(x, PartitionTuple):
            if x.level() == self._r:
                return self.monomial(x)
            raise ValueError(f"expected a partition tuple of level {self._r}, "
                             f"got level {x.level()}")
        if isinstance(x, Partition):
            if self._r == 1:
                return self.monomial(C([list(x)]))
            raise ValueError("a single partition can only be used for level 1")
        if isinstance(x, (list, tuple)):
            if x and all(isinstance(part, (list, tuple, Partition))
                         for part in x):
                if len(x) != self._r:
                    raise ValueError(f"expected a tuple of {self._r} "
                                     f"partitions, got {len(x)}")
                return self.monomial(C([list(part) for part in x]))
            if self._r == 1:
                return self.monomial(C([list(x)]))
            raise ValueError(f"expected a tuple of {self._r} partitions")
        if x in ZZ:
            if self._r == 1:
                return self.monomial(C([[ZZ(x)]]))
            raise ValueError("a single partition can only be used for level 1")
        raise TypeError(f"cannot convert {x!r} to a basis element of {self}")

    def one_basis(self):
        r"""
        Return the index of the multiplicative unit.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: H(SymmetricFunctions(QQ).p()).one_basis()
            ([], [])
        """
        return self._indices([[] for _ in range(self._r)])

    def degree_on_basis(self, tau):
        r"""
        Return the degree of the basis element ``tau``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: P.degree_on_basis(P._indices([[2, 1], [3]]))
            6
        """
        return tau.size()

    def product_on_basis(self, mu, nu):
        r"""
        Return the product of the two basis elements ``mu`` and ``nu``.

        Since the identification of `\Lambda(r)` with
        `\Lambda(1)^{\otimes r}` is an algebra map, the structure
        constants are the classical ones, applied colourwise.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: Y = H(SymmetricFunctions(QQ).s())
            sage: Y[[2], []] * Y[[], [1, 1]]
            s_2(ζ^0)*s_{1,1}(ζ^1)
            sage: Y[[2], []] * Y[[2], []]
            s_{2,2}(ζ^0) + s_{3,1}(ζ^0) + s_4(ζ^0)
        """
        B = self._classical
        d = {}
        for combo in iproduct(*[(B[mu[c]] * B[nu[c]])
                                .monomial_coefficients().items()
                                for c in range(self._r)]):
            tau = self._indices([list(lam) for lam, _ in combo])
            coeff = prod(c for _, c in combo)
            d[tau] = d.get(tau, self.base_ring().zero()) + coeff
        return self._from_dict(d)

    def _repr_term(self, tau):
        r"""
        Return a string representation of the basis element ``tau``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: P[[2], [1]]
            p_2(ζ^0)*p_1(ζ^1)
            sage: P[[2, 1], [3]]
            p_{2,1}(ζ^0)*p_3(ζ^1)
            sage: P.one()
            1
        """
        prefix = self._classical.prefix()
        glyph = 'χ' if self._indexing == 'chi' else 'ζ'
        factors = []
        for j in range(self._r):
            lam = tau[j]
            if not lam:
                continue
            if len(lam) == 1:
                factors.append(f"{prefix}_{lam[0]}({glyph}^{j})")
            else:
                factors.append(prefix + "_{" + ",".join(str(k) for k in lam) +
                               f"}}({glyph}^{j})")
        return "*".join(factors) if factors else "1"

    def _latex_term(self, tau):
        r"""
        Return a latex representation of the basis element ``tau``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: latex(P[[2, 1], [3]])
            p_{2,1}(\zeta^{0}) \cdot p_{3}(\zeta^{1})
            sage: latex(H(SymmetricFunctions(QQ).s(), 'chi')[[2], [1]])
            s_{2}(\chi^{0}) \cdot s_{1}(\chi^{1})
        """
        prefix = self._classical.prefix()
        glyph = r'\chi' if self._indexing == 'chi' else r'\zeta'
        factors = []
        for j in range(self._r):
            lam = tau[j]
            if not lam:
                continue
            factors.append(prefix + "_{" + ",".join(str(k) for k in lam) +
                           f"}}({glyph}^{{{j}}})")
        return r" \cdot ".join(factors) if factors else "1"

    @cached_method
    def _powersum(self):
        r"""
        Return the ordinary powersum basis over the same base ring.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: H(SymmetricFunctions(QQ).s())._powersum()
            Symmetric Functions over Rational Field in the powersum basis
        """
        return SymmetricFunctions(self.base_ring()).power()

    @cached_method
    def _hub(self):
        r"""
        Return the class-indexed powersum basis with the same level and
        base ring, through which all conversions pass.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: S = SymmetricFunctions(QQ)
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(S.p())
            sage: H(S.s(), 'chi')._hub() is P
            True
            sage: P._hub() is P
            True
        """
        p = self._powersum()
        if self._indexing == 'zeta' and self._classical is p:
            return self
        return HyperoctahedralSymmetricFunctions(self._r, p, 'zeta')

    @cached_method
    def _chi_power_sum(self, i, c):
        r"""
        Return the character-coloured power sum
        `\pi_i^{(c)} = p_i(\chi^c) = \frac{1}{r} \sum_j \zeta^{jc} p_i(\zeta^j)`
        (Macdonald, Appendix B, (7.1)) as a hub element.

        Only meaningful for character-indexed ``self``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: H(S.p(), 'chi')._chi_power_sum(2, 1)
            -1/2*p_2(ζ^1) + 1/2*p_2(ζ^0)
        """
        hub = self._hub()
        result = hub.zero()
        for j in range(self._r):
            parts = [[] for _ in range(self._r)]
            parts[j] = [ZZ(i)]
            result += (self._zeta ** (ZZ(j * c) % self._r)
                       * hub.monomial(hub._indices(parts)))
        return result / self._r

    def _to_hub_on_basis(self, tau):
        r"""
        Return the expansion of the basis element ``tau`` of ``self`` in
        the class-indexed powersum basis (the hub).

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: Yz = H(S.s())
            sage: Yz._to_hub_on_basis(Yz._indices([[2], [1]]))
            1/2*p_{1,1}(ζ^0)*p_1(ζ^1) + 1/2*p_2(ζ^0)*p_1(ζ^1)
        """
        p = self._powersum()
        hub = self._hub()
        result = hub.zero()
        for combo in iproduct(*[p(self._classical[lam])
                                .monomial_coefficients().items()
                                for lam in tau]):
            coeff = prod(c for _, c in combo)
            if self._indexing == 'zeta':
                result = result + coeff * hub.monomial(
                    hub._indices([list(lam) for lam, _ in combo]))
            else:
                term = hub.one()
                for c, (lam, _) in enumerate(combo):
                    for i in lam:
                        term = term * self._chi_power_sum(i, c)
                result = result + coeff * term
        return result

    def _chi_expand(self, rho):
        r"""
        Expand the hub basis element indexed by the partition tuple
        ``rho`` in character-coloured powersum monomials.

        This implements Macdonald's (7.1'), `p_i(\zeta^j) = \sum_k
        \zeta^{-jk} p_i(\chi^k)`, separately for each part size `i`:  if
        the monomial contains `m_j` parts of size `i` in class colour
        `j`, these are distributed among the character colours `k` with
        multinomial weights `\zeta^{-jk}` per part.

        OUTPUT: a dictionary mapping partition tuples of level `r` to
        coefficients.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: Q = H(S.p(), 'chi')
            sage: rho = Q._hub()._indices([[], [2]])
            sage: Q._chi_expand(rho)
            {([], [2]): -1, ([2], []): 1}
        """
        r = self._r
        zeta = self._zeta
        R = self.base_ring()
        counts = {}
        for j in range(r):
            for i in rho[j]:
                counts.setdefault(i, [0] * r)[j] += 1
        result = {tuple([()] * r): R.one()}
        for i in counts:
            m = counts[i]
            new = {}
            for key, w in result.items():
                for vecs in iproduct(*[IntegerVectors(m[j], r)
                                       for j in range(r)]):
                    coeff = w
                    parts = [list(key[c]) for c in range(r)]
                    for j in range(r):
                        if not m[j]:
                            continue
                        coeff = coeff * multinomial(*vecs[j])
                        for c in range(r):
                            n = vecs[j][c]
                            if n:
                                coeff = coeff * zeta ** (ZZ(-j * c * n) % r)
                                parts[c].extend([i] * n)
                    tau = tuple(tuple(sorted(part, reverse=True))
                                for part in parts)
                    if tau in new:
                        new[tau] = new[tau] + coeff
                    else:
                        new[tau] = coeff
            result = new
        return {self._indices([list(key[j]) for j in range(r)]): coeff
                for key, coeff in result.items()}

    def _from_hub_on_basis(self, rho):
        r"""
        Return the expansion of the hub basis element ``rho`` in the
        basis of ``self``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: Q = H(S.p(), 'chi')
            sage: rho = Q._hub()._indices([[], [2]])
            sage: Q._from_hub_on_basis(rho)
            -p_2(χ^1) + p_2(χ^0)
        """
        B = self._classical
        p = self._powersum()
        if self._indexing == 'zeta':
            monos = [(rho, self.base_ring().one())]
        else:
            monos = self._chi_expand(rho).items()
        d = {}
        for key, w in monos:
            for combo in iproduct(*[B(p[lam])
                                    .monomial_coefficients().items()
                                    for lam in key]):
                tau = self._indices([list(mu) for mu, _ in combo])
                coeff = w * prod(c for _, c in combo)
                d[tau] = d.get(tau, self.base_ring().zero()) + coeff
        return self._from_dict(d)

    @cached_method
    def _to_hub_map(self):
        r"""
        Return the module morphism from ``self`` to the hub.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: Q = H(S.p(), 'chi')
            sage: Q._to_hub_map()(Q[[], [2]])
            -1/2*p_2(ζ^1) + 1/2*p_2(ζ^0)
        """
        return self.module_morphism(on_basis=self._to_hub_on_basis,
                                    codomain=self._hub())

    @cached_method
    def _from_hub_map(self):
        r"""
        Return the module morphism from the hub to ``self``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: Q = H(S.p(), 'chi')
            sage: Q._from_hub_map()(Q._hub()[[], [2]])
            -p_2(χ^1) + p_2(χ^0)
        """
        return self._hub().module_morphism(on_basis=self._from_hub_on_basis,
                                           codomain=self)

    def _coerce_map_from_(self, other):
        r"""
        Return a coercion map from ``other`` to ``self``, if there is one.

        Any basis of the same level and base ring is connected to
        ``self`` through the class-indexed powersum basis (the hub), so
        the composite morphism through the hub is returned.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: S = SymmetricFunctions(QQ)
            sage: P = H(S.p()); Q = H(S.p(), 'chi')
            sage: P.has_coerce_map_from(Q)
            True
            sage: P(Q[[], [2]])
            -1/2*p_2(ζ^1) + 1/2*p_2(ζ^0)
            sage: P[[2], []] + Q[[2], []]
            1/2*p_2(ζ^1) + 3/2*p_2(ζ^0)
        """
        if isinstance(other, HyperoctahedralSymmetricFunctions):
            if other is self:
                return True
            if (other._r == self._r
                    and other.base_ring() == self.base_ring()):
                return self._from_hub_map() * other._to_hub_map()
        return super()._coerce_map_from_(other)

    @cached_method
    def _tensor(self):
        r"""
        Return the tensor product realizing `\Lambda(r) \cong \Lambda(1)^{\otimes r}`.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: H(SymmetricFunctions(QQ).p())._tensor()
            Symmetric Functions over Rational Field in the powersum basis # Symmetric Functions over Rational Field in the powersum basis
        """
        return tensor([self._powersum()] * self._r)

    def _to_tensor(self, g):
        r"""
        Convert the element ``g`` of ``self`` to the tensor realization.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: P._to_tensor(P[[2], []] + P[[], [1]])
            p[] # p[1] + p[2] # p[]
        """
        T = self._tensor()
        return T._from_dict({tuple(tau): c for tau, c in self._to_hub_map()(g)},
                            remove_zeros=False)

    def _from_tensor(self, x):
        r"""
        Convert the tensor element ``x`` to an element of ``self``.

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: T = P._tensor()
            sage: P._from_tensor(T._from_dict({((2,), (1,)): QQ(3)}, remove_zeros=False))
            3*p_2(ζ^0)*p_1(ζ^1)
        """
        hub = self._hub()
        return self._from_hub_map()(hub._from_dict(
            {hub._indices(list(key)): c for key, c in x},
            remove_zeros=False))

    def _plethysm_type1(self, f, g):
        r"""
        Return the plethysm `f \circ g`, where `f` is ordinary.

        This is Henderson's first plethysm
        `\Lambda(1) \times \Lambda(r) \to \Lambda(r)`, defined by

        .. MATH::

            p_i \circ p_j(\zeta^a) = p_{ij}(\zeta^a).

        It is computed by identifying `\Lambda(r)` with
        `\Lambda(1)^{\otimes r}` and using Sage's plethysm of ordinary
        symmetric functions on each tensor factor.

        INPUT:

        - ``f`` -- an ordinary symmetric function
        - ``g`` -- an element of ``self``

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: p = SymmetricFunctions(QQ).p()
            sage: P._plethysm_type1(p[2], P[[3], []])
            p_6(ζ^0)
            sage: P._plethysm_type1(p[2], P[[2], []]*P[[3], []])
            p_{6,4}(ζ^0)
        """
        p = f.parent().realization_of().power()
        return self._from_tensor(p(f)(self._to_tensor(g)))

    def _plethysm_type2(self, f, g):
        r"""
        Return the plethysm `f \circ g`, where ``g`` is ordinary.

        This is Henderson's second plethysm
        `\Lambda(r) \times \Lambda(1) \to \Lambda(r)`, defined by

        .. MATH::

            p_i(\zeta^a) \circ p_j = p_{ij}(\zeta^{aj}).

        INPUT:

        - ``f`` -- an element of ``self``
        - ``g`` -- an ordinary symmetric function

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(3)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: p = SymmetricFunctions(QQ).p()
            sage: P._plethysm_type2(P[[], [2], []], p[2])
            p_4(ζ^2)
            sage: P._plethysm_type2(P[[], [2], []], p[2] + p[1, 1])
            p_4(ζ^2) + p_{2,2}(ζ^1)
        """
        p = g.parent().realization_of().power()
        g = p(g)
        hub = self._hub()
        result = hub.zero()
        for nu, c in self._to_hub_map()(f):
            term = hub.one()
            for a in range(self._r):
                for k in nu[a]:
                    term = term * hub._type2_factor(a, k, g)
            result = result + c * term
        return self._from_hub_map()(result)

    def _type2_factor(self, a, k, g):
        r"""
        Return `p_k(\zeta^a) \circ g` for an ordinary power sum ``g``.

        This must only be called when ``self`` is a class-indexed
        powersum basis (the hub), since it constructs powersum monomials.

        INPUT:

        - ``a`` -- the colour of the left power sum
        - ``k`` -- the degree of the left power sum
        - ``g`` -- an ordinary symmetric function in the power-sum basis

        EXAMPLES::

            sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
            sage: H = HyperoctahedralSymmetricFunctions(2)
            sage: P = H(SymmetricFunctions(QQ).p())
            sage: p = SymmetricFunctions(QQ).p()
            sage: P._type2_factor(1, 2, p[3])
            p_6(ζ^1)
            sage: P._type2_factor(1, 2, p[2])
            p_4(ζ^0)
        """
        result = self.zero()
        for lam, d in g:
            parts = [[] for _ in range(self._r)]
            for j in lam:
                parts[(a * j) % self._r].append(k * j)
            result += d * self.monomial(self._indices(parts))
        return result

    class Element(CombinatorialFreeModule.Element):
        r"""
        An element of `\Lambda(r)`.
        """
        def plethysm(self, x):
            r"""
            Return the plethysm of ``self`` with ``x``.

            The dispatch follows Henderson: if ``x`` is an ordinary
            symmetric function, the second plethysm
            `\Lambda(r) \times \Lambda(1) \to \Lambda(r)` is used.  If
            ``x`` is an element of `\Lambda(r)`, the `r`-by-`r`
            plethysm is not implemented yet.

            EXAMPLES::

                sage: from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
                sage: H = HyperoctahedralSymmetricFunctions(3)
                sage: P = H(SymmetricFunctions(QQ).p())
                sage: p = SymmetricFunctions(QQ).p()
                sage: P[[], [2], []].plethysm(p[2])
                p_4(ζ^2)
                sage: P[[], [2], []].plethysm(p[2] + p[1, 1])
                p_4(ζ^2) + p_{2,2}(ζ^1)

            The `r`-by-`r` plethysm is not implemented yet::

                sage: P[[], [2], []].plethysm(P[[3], [], []])
                Traceback (most recent call last):
                ...
                NotImplementedError: plethysm of two r-species is not implemented yet
            """
            if isinstance(x, HyperoctahedralSymmetricFunctions.Element):
                raise NotImplementedError("plethysm of two r-species is not implemented yet")
            return self.parent()._plethysm_type2(self, x)
