r"""
Lazy polynomial `r`-species

The graded algebra of polynomial `r`-species, described in
:mod:`sage.rings.species_hyperoctahedral`, is completed here in the same
way as the ring of polynomial species is completed in
:mod:`sage.rings.lazy_species`: the coefficient of degree `n` of a lazy
`r`-species is a finite linear combination of molecular `r`-species of
degree `n`, where the degree is the number of `C_r`-blocks.  Addition
and multiplication are performed coefficientwise.

In addition, both kinds of substitutions introduced by Henderson are
available: an ordinary lazy species `F`, possibly multisort, can be
evaluated at lazy `r`-species `G_1, \dots, G_k`, where `k` is the
number of sorts of `F`.  Also, an `r`-species `F`, possibly multisort,
can be composed with ordinary species `G_1, \dots, G_k`, where `k` is
again the number of sorts of `F`. In both cases, the result is again a
lazy `r`-species.

EXAMPLES:

We construct the ring of lazy `2`-species over the rationals.  The
generator `X` corresponds to the trivial subgroup of the wreath product
`W(r, 1)`, the molecular species `X^\circ` corresponds to its full
cyclic stabilizer::

    sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
    sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
    sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
    sage: Xo = L(_wreath_group(2, 1))
    sage: X, Xo
    (X, X°)

Addition and multiplication are coefficientwise::

    sage: F = X + Xo
    sage: F
    (X+X°)
    sage: F**2
    (X^2+2*X*X°+X°^2)

The generating series of a lazy `r`-species `F` is
`\sum_{n \geq 0} |F[n]_r| x^n/(n! r^n)`, where `F[n]_r` denotes the
`F`-structures on the canonical free `W(r,n)`-set with `n` blocks.  For
molecular `r`-species, this amounts to replacing each molecular species
`M` with coefficient `c` by `c / |H_M|`, where `H_M` is the stabilizer
group of `M`::

    sage: F.generating_series()
    3/2*X + O(X^7)
    sage: Xo.generating_series()
    1/2*X + O(X^7)

The isotype generating series counts isomorphism classes of structures::

    sage: F.isotype_generating_series()
    2*X + O(X^7)

Substituting the molecular species `X^\circ` into the species of sets
provides the species `E^{(r)}` of `r`-sets, i.e., the unique
`r`-species with precisely one structure in each degree::

    sage: from sage.rings.lazy_species import LazyCombinatorialSpecies
    sage: E = LazyCombinatorialSpecies(QQ, "Z").Sets()
    sage: Er = E(Xo); Er
    1 + X° + E_2(X°) + E_3(X°) + E_4(X°) + E_5(X°) + E_6(X°) + O^7
    sage: Er.generating_series()
    1 + 1/2*X + 1/8*X^2 + 1/48*X^3 + 1/384*X^4 + 1/3840*X^5 + 1/46080*X^6 + O(X^7)
    sage: exp(Xo.generating_series())[:7] == Er.generating_series()[:7]
    True

The type 1 substitution distributes over addition and scalar
multiplication, and is associative::

    sage: from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    sage: L1 = LazyCombinatorialSpecies(QQ, "Z")
    sage: E_2 = L1(SymmetricGroup(2))
    sage: E_2(X + Xo)[2]
    E_2(X) + X*X° + E_2(X°)
    sage: E_2(2*X)[2]
    2*E_2(X) + X^2
    sage: E_2(Xo)[2]
    E_2(X°)
    sage: E_3 = L1(SymmetricGroup(3))
    sage: E_2(E_3(F))[4] == E_2(E_3)(F)[4]
    True

The generating series of a substitution is the substitution of the
generating series, see [Henderson2004]_, Equation (1.5), Corollary 4.7
and Corollary 4.10::

    sage: OC = L.OrientedCycles()
    sage: E1 = E.restrict(1)
    sage: OC(E1).generating_series() - OC.generating_series()(E1.generating_series())
    O(Z^7)
    sage: E(OC).generating_series() - E.generating_series()(OC.generating_series())
    O(X^7)

    sage: E(X + Xo).generating_series()
    1 + 3/2*X + 9/8*X^2 + 9/16*X^3 + 27/128*X^4 + 81/1280*X^5 + 81/5120*X^6 + O(X^7)
    sage: E(X + Xo).generating_series() - E.generating_series()((X + Xo).generating_series())
    O(X^7)

Multisort `r`-species are supported::

    sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
    sage: P2 = PolynomialHyperoctahedralSpecies(QQ, 2, "X2, Y2")
    sage: L2.<X2, Y2> = LazyHyperoctahedralSpecies(QQ, 2)
    sage: Yo = L2(P2(_wreath_group(2, 1), {1: [1, 2]}))
    sage: E(X2 + Yo)[:3]
    [1, X2 + Y2°, E_2(X2) + X2*Y2° + E_2(Y2°)]

For `r = 1` the type 1 substitution reduces to the usual substitution
of species::

    sage: L1r.<X1> = LazyHyperoctahedralSpecies(QQ, 1)
    sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
    sage: P1 = PolynomialHyperoctahedralSpecies(QQ, 1, "X1")
    sage: E_2(X1)[2] == L1r(P1(SymmetricGroup(2)))[2]
    True

Not every lazy `r`-species is a type 1 substitution of an ordinary
species.  For example, the species of signed graphs is available::

    sage: S = L.SignedGraphs()
    sage: S[2]
    E_2(X°) + X°(E_2)

Likewise, the species of Seidel graphs, i.e., of labelings of the complete
graph with elements of `Z_r` acted upon by switching, is not a type 1
substitution.  For `r = 2` its isomorphism types are the two-graphs::

    sage: L.SeidelGraphs()[2]
    X°(E_2)

Henderson's trees can be implemented as follows.  First, using
Proposition 5.3, we define the trees for the ordinary case::

    sage: L1 = LazyCombinatorialSpecies(QQ, "X")
    sage: X1 = L1._first_ngens(1)[0]
    sage: E = L1.Sets()
    sage: T1 = L1.undefined()
    sage: T1.define(X1 + E.restrict(2)(T1))
    sage: [T1.generating_series()[n]*factorial(n) for n in range(4)]
    [0, 1, 1, 4]

Then we use Proposition 5.5::

    sage: T = L.undefined()
    sage: T.define((1 + T) * E(Xo).restrict(1)(T1))
    sage: [T.generating_series()[n]*factorial(n)*2^n for n in range(4)]
    [0, 1, 5, 47]

REFERENCES:

.. [Henderson2004] Anthony Henderson.
   *Representations of wreath products on cohomology of De Concini-Procesi
   compactifications*.
   International Mathematics Research Notices 2004, no. 20, 981-1021.

AUTHORS:

- Martin Rubey (2026): initial version

"""
# ****************************************************************************
#       Copyright (C) 2026 Martin Rubey <martin.rubey@tuwien.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools
from collections import defaultdict

from sage.categories.tensor import tensor
from sage.combinat.sf.sf import SymmetricFunctions
from sage.graphs.graph import Graph
from sage.graphs.graph_generators import graphs
from sage.groups.perm_gps.constructor import PermutationGroupElement
from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_list import lazy_list
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.lazy_series import LazyCompletionGradedAlgebraElement
from sage.rings.lazy_series_ring import (LazyCompletionGradedAlgebra,
                                         LazyPowerSeriesRing,
                                         LazySymmetricFunctions)
from sage.rings.lazy_species import (LazyCombinatorialSpecies,
                                     LazyCombinatorialSpeciesElement,
                                     weighted_vector_compositions)
from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
from sage.rings.species import PolynomialSpecies
from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
from sage.data_structures.stream import (Stream_exact,
                                         Stream_function,
                                         Stream_zero)
from sage.structure.element import get_coercion_model, parent
from sage.structure.unique_representation import UniqueRepresentation


class LazyHyperoctahedralSpeciesElement(LazyCompletionGradedAlgebraElement):
    r"""
    A lazy `r`-species.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
        sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
        sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
        sage: Xo = L(_wreath_group(2, 1))
        sage: F = X + Xo
        sage: F
        (X+X°)
        sage: F**2
        (X^2+2*X*X°+X°^2)
    """

    def generating_series(self):
        r"""
        Return the generating series of ``self``.

        The coefficient of `x^n` is `|F[n]_r|/(n! r^n)`, where `F[n]_r`
        denotes the set of structures of ``self`` on the canonical free
        `W(r,n)`-set with `n` blocks.

        For multisort species the variables correspond to the sorts,
        the coefficient of `x_1^{n_1} \cdots x_k^{n_k}` is
        `|F[n_1, \ldots, n_k]_r| / |W(r; n_1, \ldots, n_k)|`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: (X + Xo).generating_series()
            3/2*X + O(X^7)

            sage: E = LazyCombinatorialSpecies(QQ, "Z").Sets()
            sage: E(Xo).generating_series()
            1 + 1/2*X + 1/8*X^2 + 1/48*X^3 + 1/384*X^4 + 1/3840*X^5 + 1/46080*X^6 + O(X^7)

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: P2 = PolynomialHyperoctahedralSpecies(QQ, 2, "X, Y")
            sage: L2.<X, Y> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Yo = L2(P2(_wreath_group(2, 1), {1: [1, 2]}))
            sage: E(X + Yo).generating_series()[:3]
            [1, X + 1/2*Y, 1/2*X^2 + 1/2*X*Y + 1/8*Y^2]
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        if P._arity == 1:
            def coefficient(n):
                return sum(c / M.permutation_group()[0].cardinality()
                           for M, c in self[n].monomial_coefficients().items())
        else:
            def coefficient(n):
                return sum(c / M.permutation_group()[0].cardinality()
                           * P.base_ring().prod(v ** d for v, d in zip(L.gens(), M.grade()))
                           for M, c in self[n].monomial_coefficients().items())
        return L(coefficient)

    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of ``self``.

        The coefficient of `x^n` is the number of isomorphism classes of
        structures of ``self`` on the canonical free `W(r,n)`-set with
        `n` blocks.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: (X + Xo).isotype_generating_series()
            2*X + O(X^7)

            sage: E = LazyCombinatorialSpecies(QQ, "Z").Sets()
            sage: E(Xo).isotype_generating_series()
            1 + X + X^2 + X^3 + X^4 + X^5 + X^6 + O(X^7)
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        if P._arity == 1:
            def coefficient(n):
                return sum(c for M, c in self[n].monomial_coefficients().items())
        else:
            def coefficient(n):
                return sum(c * P.base_ring().prod(v ** d for v, d in zip(L.gens(), M.grade()))
                           for M, c in self[n].monomial_coefficients().items())
        return L(coefficient)

    def cycle_index_series(self):
        r"""
        Return the cycle index series of ``self``.

        The coefficient of degree `n` is the sum of the cycle indices of
        the molecular `r`-species of ``self[n]``, in the power sum basis
        of the hyperoctahedral symmetric functions
        :mod:`sage.rings.sf_hyperoctahedral` over the fraction field of
        the base ring, see
        :meth:`~sage.rings.species_hyperoctahedral.MolecularHyperoctahedralSpecies.Element.cycle_index`.

        For multisort species, the result lives in the lazy completion of
        a tensor product of hyperoctahedral symmetric functions, one
        factor per sort, and is graded by the total degree.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.rings.lazy_species import LazyCombinatorialSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group, _wreath_young_subgroup
            sage: L = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: E2Xo = L(_wreath_young_subgroup(2, [2]))
            sage: Z = (Xo + E2Xo).cycle_index_series()
            sage: Z
            (1/2*p_1(ζ^1)+1/2*p_1(ζ^0)) + (1/8*p_{1,1}(ζ^1)+1/4*p_2(ζ^1)+1/4*p_1(ζ^0)*p_1(ζ^1)+1/8*p_{1,1}(ζ^0)+1/4*p_2(ζ^0)) + O^7
            sage: Z[1]
            1/2*p_1(ζ^1) + 1/2*p_1(ζ^0)
            sage: Z[2]
            1/8*p_{1,1}(ζ^1) + 1/4*p_2(ζ^1) + 1/4*p_1(ζ^0)*p_1(ζ^1) + 1/8*p_{1,1}(ζ^0) + 1/4*p_2(ζ^0)

        In particular, the cycle index series of the `r`-species of
        sets `E^{(r)} = E(X^\circ)` is Henderson's
        `exp(\sum_i \sum_a p_i(\zeta^a)/(r i))` [Henderson2004, (4.1)]_:

            sage: E = LazyCombinatorialSpecies(QQ, "Z").Sets()
            sage: ZE = E(Xo).cycle_index_series()
            sage: ZE[3]
            1/48*p_{1,1,1}(ζ^1) + 1/8*p_{2,1}(ζ^1) + 1/6*p_3(ζ^1) + 1/16*p_1(ζ^0)*p_{1,1}(ζ^1) + 1/8*p_1(ζ^0)*p_2(ζ^1) + 1/16*p_{1,1}(ζ^0)*p_1(ζ^1) + 1/48*p_{1,1,1}(ζ^0) + 1/8*p_2(ζ^0)*p_1(ζ^1) + 1/8*p_{2,1}(ζ^0) + 1/6*p_3(ζ^0)

        Infinite series are supported::

            sage: F = 1/(2 - Xo)
            sage: ZF = F.cycle_index_series()
            sage: ZF[1]
            1/8*p_1(ζ^1) + 1/8*p_1(ζ^0)
            sage: ZF[2]
            1/32*p_{1,1}(ζ^1) + 1/16*p_1(ζ^0)*p_1(ζ^1) + 1/32*p_{1,1}(ζ^0)

        For multisort species, the result is graded by the total
        degree::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: LXY = LazyHyperoctahedralSpecies(QQ, 2, "X, Y")
            sage: AXY = AtomicHyperoctahedralSpecies(2, "X, Y")
            sage: d0 = AXY(_wreath_group(2, 1), {0: [1, 2]})
            sage: d1 = AXY(_wreath_group(2, 1), {1: [1, 2]})
            sage: W = _wreath_young_subgroup(2, [1, 1])
            sage: dd = AXY(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
            sage: ZXY = (LXY(d0) + LXY(d1) + LXY(dd)).cycle_index_series()
            sage: ZXY[1]
            1/2*1 # p_1(ζ^1) + 1/2*1 # p_1(ζ^0) + 1/2*p_1(ζ^1) # 1 + 1/2*p_1(ζ^0) # 1
            sage: ZXY[2]
            1/2*p_1(ζ^1) # p_1(ζ^1) + 1/2*p_1(ζ^0) # p_1(ζ^0)

        The coefficient of the key `([1, \dots, 1], [], \dots, [])` of
        ``Z[n]`` is the coefficient of `x^n` of the generating series,
        and the sum of the coefficients of ``Z[n]`` is the coefficient
        of `x^n` of the isotype generating series::

            sage: from sage.combinat.partition_tuple import PartitionTuples_level
            sage: Pt = PartitionTuples_level(2)
            sage: Z[2].coefficient(Pt([[1, 1], []])) == (Xo + E2Xo).generating_series()[2]
            True
            sage: sum(Z[2].coefficients()) == (Xo + E2Xo).isotype_generating_series()[2]
            True
        """
        P = self.parent()
        r = P._internal_poly_ring.base_ring()._r
        H = HyperoctahedralSymmetricFunctions(
            r, SymmetricFunctions(P.base_ring().fraction_field()).powersum())
        if P._arity == 1:
            L = LazySymmetricFunctions(H)

            def coefficient(n):
                return sum(c * M.cycle_index(parent=H)
                           for M, c in self[n].monomial_coefficients().items())
        else:
            T = tensor([H for _ in range(P._arity)])
            L = LazySymmetricFunctions(T)

            def coefficient(n):
                return sum(c * M.cycle_index(parent=T)
                           for M, c in self[n].monomial_coefficients().items())
        return L(coefficient)

    def structures(self, *labels):
        r"""
        Iterate over the structures on the given free `C_r`-set of labels.

        The labels are given as one list of `C_r`-orbits per sort, each
        orbit an iterable of exactly `r` labels in cyclic phase order;
        the degree is the number of orbits.  See
        :meth:`~sage.rings.species_hyperoctahedral.PolynomialHyperoctahedralSpecies.Element.structures`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: X = L(_wreath_group(2, 1).subgroup([]))
            sage: list(Xo.structures([('a', 'b')]))
            [(X°, (('a', 'b'),))]
            sage: sorted((X + Xo).structures([('a', 'b')]), key=str)
            [(X, (('a', 'b'),)), (X, (('b', 'a'),)), (X°, (('a', 'b'),))]

            The number of structures agrees with `r^n n!` times the
            coefficient of the generating series::

            sage: E = LazyCombinatorialSpecies(QQ, "Z").Sets()
            sage: sorted(E(Xo).structures([('a', 'b'), ('c', 'd'), ('e', 'f')]))
            [(E_3(X°), (('a', 'b', 'c', 'd', 'e', 'f'),))]
            sage: E(Xo).generating_series()[3]
            1/48
        """
        yield from self[sum(map(len, labels))].structures(*labels)

    def _test_structures(self, tester=None, max_size=4, **options):
        r"""
        Check that structures and generating series are consistent.

        We check all structures with less than ``max_size`` orbits.

        TESTS::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group, _wreath_young_subgroup
            sage: L = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: X = L(_wreath_group(2, 1).subgroup([]))
            sage: (X + Xo)._test_structures()
            sage: (X*Xo)._test_structures()

            Multisort species are checked for every vector of orbit
            counts::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies, AtomicHyperoctahedralSpecies
            sage: LXY = LazyHyperoctahedralSpecies(QQ, 2, "X, Y")
            sage: MXY = MolecularHyperoctahedralSpecies(2, "X, Y")
            sage: AXY = AtomicHyperoctahedralSpecies(2, "X, Y")
            sage: W = _wreath_young_subgroup(2, [1, 1])
            sage: d = AXY(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
            sage: LXY(MXY({d: 1}))._test_structures()
        """
        if tester is None:
            tester = self._tester(**options)
        P = self.parent()
        r = int(P._internal_poly_ring.base_ring()._r)
        for n in range(max_size):
            if P._arity == 1:
                labels = [[(i, j) for j in range(r)] for i in range(n)]
                s = list(self.structures(labels))
                tester.assertEqual(len(s), len(set(s)),
                                   f"structures for {labels} are {s}, which is not a set")
                coeff = self.generating_series()[n]
                tester.assertEqual(len(s) / ZZ(r)**n / ZZ(n).factorial(), coeff,
                                   f"the number of structures for {labels} is {len(s)}, "
                                   f"but the generating series gives {coeff}")
            else:
                from sage.combinat.integer_vector import IntegerVectors
                for shape in IntegerVectors(n, length=P._arity):
                    labels = [[[(i, j, sort) for j in range(r)]
                               for i in range(k)]
                              for sort, k in enumerate(shape)]
                    s = list(self.structures(*labels))
                    tester.assertEqual(len(s), len(set(s)),
                                       f"structures for {labels} are {s}, which is not a set")
                    coeff = self.generating_series()[n].coefficient(list(shape))
                    scale = ZZ.prod(ZZ(r)**k * ZZ(k).factorial() for k in shape)
                    tester.assertEqual(len(s) / scale, coeff,
                                       f"the number of structures for {labels} is {len(s)}, "
                                       f"but the generating series gives {coeff}")


    def __call__(self, *args):
        r"""
        Return the substitution of ``args`` into ``self``.

        This implements the composition `F \circ_2 (G_1, \ldots, G_k)`
        of [Henderson2004]_, Equation (4.8), where ``self`` is a lazy
        `C_r`-equivariant species and each `G_i` is an ordinary lazy
        species, see
        :meth:`~sage.rings.species_hyperoctahedral.PolynomialHyperoctahedralSpecies.Element.__call__`.

        Substitution is associative, and linear in the outer species.

        The arguments may be multisort, in which case the composite is
        an `r`-species in the sorts of the arguments.  This is
        analogous to ordinary composition of species.

        EXAMPLES:

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: L1.<Z> = LazyCombinatorialSpecies(QQ)

        Substituting an ordinary species into an `r`-species yields
        the argument as an `r`-species::

            sage: X(Z + Z^2)
            Z + Z^2 + O^8
            sage: Xo(Z + Z^2)
            Z° + {((1,2)(3,4),)} + O^8

        Substitution is compatible with passing to the cycle index
        series and the generating series, see [Braunsteiner2010]_,
        Theorem 4.1.9 and [Henderson2004]_, Equation (1.5), Corollary
        4.7 and Corollary 4.10.  For example::

            sage: OC = L.OrientedCycles()
            sage: E = L1.Sets()
            sage: E1 = E.restrict(1)
            sage: OC(E1).generating_series() - OC.generating_series()(E1.generating_series())
            O(Z^7)
            sage: E(OC).generating_series() - E.generating_series()(OC.generating_series())
            O(X^7)

        Weighted species are supported.  The composite carries the
        weights of both ``self`` and the arguments, so the weights of
        ``self`` are re-based to their common base ring::

            sage: R.<q> = QQ[]
            sage: L1q.<Zq> = LazyCombinatorialSpecies(R)
            sage: OC((1+q)*Zq).truncate(4)
            ((q+1)*Zq°) + ((q^2+2*q+1)*Zq°^2) + ((q^3+q^2+q+1)*E_2(Zq°)*Zq°+(q^2+q)*Zq°^3)

        For `r = 1` the composition specializes to the ordinary
        composition of species::

            sage: L1r.<X1> = LazyHyperoctahedralSpecies(QQ, 1)
            sage: E_2 = L1(SymmetricGroup(2))
            sage: C_4 = L1r(CyclicPermutationGroup(4))
            sage: C_4(E_2)[8]
            {((7,8), (1,3,5,7)(2,4,6,8))}
            sage: L1(CyclicPermutationGroup(4))(E_2)[8]
            {((7,8), (1,3,5,7)(2,4,6,8))}

        TESTS:

        Check the case of arity zero::

            sage: L0 = LazyHyperoctahedralSpecies(QQ, 2, [])
            sage: L0.one()()
            1
            sage: (5*L0.one())()
            5

            sage: Xo()
            Traceback (most recent call last):
            ...
            ValueError: number of args must match arity of self
            sage: Xo(2)
            Traceback (most recent call last):
            ...
            ValueError: all args must be ordinary lazy species
            sage: Xo(E_2, E_2)
            Traceback (most recent call last):
            ...
            ValueError: number of args must match arity of self
            sage: Xo(Xo)
            Traceback (most recent call last):
            ...
            ValueError: all args must be ordinary lazy species
            sage: Xo(L1.one())
            Traceback (most recent call last):
            ...
            ValueError: can only compose with a positive valuation series

            sage: E_2(Xo(E_2))._test_structures()
            sage: E_2(Xo(2*E_2))._test_structures()
            sage: TestSuite(Xo(E_2)).run(skip=['_test_category', '_test_pickling'])

        """
        if not args and self.parent()._arity == 0:
            return self
        return Type2CompositionSpeciesElement(self, *args)


class Type1CompositionSpeciesElement(LazyHyperoctahedralSpeciesElement):
    r"""
    The type 1 substitution of lazy `r`-species into an ordinary lazy
    species.
    """

    def __init__(self, left, *args):
        r"""
        Initialize the type 1 substitution of ``args`` into ``left``.

        INPUT:

        - ``left`` -- an ordinary lazy species with `k` sorts

        - ``args`` -- `k` lazy `r`-species, all with the same `r` and
          sorts

        TESTS::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.rings.lazy_species import LazyCombinatorialSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: L1 = LazyCombinatorialSpecies(QQ, "Z")
            sage: E_2 = L1(SymmetricGroup(2))
            sage: F = E_2(X + Xo)
            sage: TestSuite(F).run(skip=['_test_category', '_test_pickling'])

            sage: L1.zero()(Xo)
            0
            sage: E_2(L.zero())
            0
            sage: (1 + E_2)(L.zero())
            1
            sage: (1 + E_2)(Xo)
            1 + E_2(X°) + O^7

        Substituting a constant series raises an error::

            sage: E_2(1 + Xo)
            Traceback (most recent call last):
            ...
            ValueError: can only compose with a positive valuation series

        All substituted species must be lazy `r`-species with common
        `r` and sorts::

            sage: L2.<A, B> = LazyCombinatorialSpecies(QQ)
            sage: L3.<Z3> = LazyHyperoctahedralSpecies(QQ, 3)
            sage: (A + B)(X, Z3)
            Traceback (most recent call last):
            ...
            ValueError: unable to find a common parent for the substituted species (X, Z3)
        """
        fP = left.parent()
        # Find a good parent for the result
        cm = get_coercion_model()
        try:
            P = cm.common_parent(*[parent(g) for g in args])
        except TypeError:
            raise ValueError(f"unable to find a common parent for the "
                             f"substituted species {args}")
        if not isinstance(P, LazyHyperoctahedralSpecies):
            raise ValueError(f"the substituted species {args} must be "
                             f"lazy r-species with the same r and sorts")
        try:
            BR = cm.common_parent(fP.base_ring(), P.base_ring())
        except TypeError:
            raise ValueError(f"unable to find a common base ring for {left} "
                            f"and the substituted species {args}")
        if P.base_ring() is BR:
            args = [P(g) for g in args]
        else:
            # the args stay in their own ring, whose weights coerce
            # into the common base ring of the result
            P = LazyHyperoctahedralSpecies(
                BR, P._laurent_poly_ring._r,
                P._laurent_poly_ring._indices._indices._names,
                sparse=P._sparse)

        R = P._internal_poly_ring.base_ring()
        molecules = R._indices

        # f = 0
        if isinstance(left._coeff_stream, Stream_zero):
            super().__init__(P, Stream_zero())
            self._left = left
            self._args = args
            return

        # all substituted species are zero, so the composition is constant
        if all(isinstance(g._coeff_stream, Stream_zero) for g in args):
            c = left[0]
            if c:
                c = list(c)[0][1]
            else:
                c = P.base_ring().zero()
            super().__init__(P, P(c)._coeff_stream)
            self._left = left
            self._args = args
            return

        # the outer species is a constant polynomial
        if (isinstance(left._coeff_stream, Stream_exact)
                and not left._coeff_stream._constant
                and left._coeff_stream._degree == 1):
            c = left._coeff_stream[0]
            B = c.parent()
            if not (B is ZZ or B is QQ or B == fP.base_ring()):
                c = c.coefficients()[0]
            super().__init__(P, P(c)._coeff_stream)
            self._left = left
            self._args = args
            return

        for g in args:
            if g._coeff_stream._approximate_order == 0:
                if not g._coeff_stream.is_uninitialized() and g[0]:
                    raise ValueError("can only compose with a positive valuation series")
                g._coeff_stream._approximate_order = 1

        sorder = left._coeff_stream._approximate_order
        gv = min(g._coeff_stream._approximate_order for g in args)
        L = fP._internal_poly_ring.base_ring()
        if L.base_ring() is R.base_ring():
            LR = L

            def lcoeff(c):
                return c
        else:
            # the homogeneous components of the outer species carry its
            # weights, so they have to be re-based to the common base
            # ring, in which also the weights of the args live
            LR = PolynomialSpecies(
                R.base_ring(),
                fP._laurent_poly_ring._indices._indices.variable_names())

            def lcoeff(c):
                return R.base_ring()(c)

        def coeff(g, i):
            c = g._coeff_stream[i]
            if not isinstance(c, PolynomialSpecies.Element):
                return R(c)
            return c

        # args_flat and weights contain one list for each substituted species
        weight_exp = [lazy_list(lambda j, g=g: len(coeff(g, j+1)))
                      for g in args]

        def flat(g):
            # function needed to work around python's scoping rules
            return itertools.chain.from_iterable(coeff(g, j) for j in itertools.count())

        args_flat1 = [lazy_list(flat(g)) for g in args]

        def coefficient(n):
            if not n:
                if left[0]:
                    return R(list(left[0])[0][1])
                return R.zero()
            result = R.zero()
            for i in range(1, n // gv + 1):
                # skip i=0 because it produces a term only for n=0

                # compute homogeneous components
                lF = defaultdict(LR)
                for M, c in left[i]:
                    lF[M.grade()] += LR._from_dict({M: lcoeff(c)})
                for mc, F in lF.items():
                    for degrees in weighted_vector_compositions(mc, n, weight_exp):
                        args_flat = [list(a[0:len(degrees[j])])
                                     for j, a in enumerate(args_flat1)]
                        multiplicities = [c for alpha, g_flat in zip(degrees, args_flat)
                                          for d, (_, c) in zip(alpha, g_flat) if d]
                        subs = [M for alpha, g_flat in zip(degrees, args_flat)
                                for d, (M, _) in zip(alpha, g_flat) if d]
                        non_zero_degrees = [[d for d in alpha if d] for alpha in degrees]
                        names = ["X%s" % i for i in range(len(subs))]
                        FX = F._compose_with_weighted_singletons(names,
                                                                 multiplicities,
                                                                 non_zero_degrees)
                        FG = [(molecules._type1_substitute_molecular(M, subs), c)
                              for M, c in FX]
                        result += R.sum_of_terms(FG)
            return result

        coeff_stream = Stream_function(coefficient, P._sparse, sorder * gv)
        super().__init__(P, coeff_stream)
        self._left = left
        self._args = args

    def generating_series(self):
        r"""
        Return the generating series of ``self``.

        The generating series of a type 1 substitution is the
        substitution of the generating series [Henderson2004, Section 4].

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.rings.lazy_species import LazyCombinatorialSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: E = LazyCombinatorialSpecies(QQ, "Z").Sets()
            sage: E(X + Xo).generating_series()
            1 + 3/2*X + 9/8*X^2 + 9/16*X^3 + 27/128*X^4 + 81/1280*X^5
             + 81/5120*X^6 + O(X^7)
            sage: exp(Xo.generating_series())[:7] == E(Xo).generating_series()[:7]
            True
        """
        return self._left.generating_series()(*[G.generating_series()
                                                for G in self._args])


class Type2CompositionSpeciesElement(LazyHyperoctahedralSpeciesElement):
    r"""
    The type 2 substitution of ordinary lazy species into a lazy
    `r`-species.
    """
    def __init__(self, left, *args):
        r"""
        Initialize the type 2 substitution of ``args`` into ``left``.

        INPUT:

        - ``left`` -- a lazy `r`-species with `k` sorts

        - ``args`` -- `k` ordinary lazy species with the same sorts

        TESTS::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.rings.lazy_species import LazyCombinatorialSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Xo = L(_wreath_group(2, 1))
            sage: L1.<Z> = LazyCombinatorialSpecies(QQ)
            sage: E2 = L1(SymmetricGroup(2))
            sage: F = Xo(E2)
            sage: TestSuite(F).run(skip=['_test_category', '_test_pickling'])

            sage: L.zero()(E2)
            0
            sage: Xo(L1.zero())
            0
            sage: (1 + Xo)(L1.zero())
            1
            sage: (1 + Xo)(E2)
            1 + Z°(E_2) + O^7

        Substituting a constant series raises an error::

            sage: Xo(1 + E2)
            Traceback (most recent call last):
            ...
            ValueError: can only compose with a positive valuation series

        All substituted species must be ordinary lazy species with the
        same sorts::

            sage: Xo(Xo)
            Traceback (most recent call last):
            ...
            ValueError: all args must be ordinary lazy species
            sage: Xo(E2, E2)
            Traceback (most recent call last):
            ...
            ValueError: number of args must match arity of self
        """
        fP = left.parent()
        if len(args) != fP._arity:
            raise ValueError("number of args must match arity of self")
        if not all(isinstance(g, LazyCombinatorialSpeciesElement)
                   for g in args):
            raise ValueError("all args must be ordinary lazy species")

        # Find a good parent for the result: an r-species with the
        # sorts of the args over the common base ring
        cm = get_coercion_model()
        try:
            P0 = cm.common_parent(*[parent(g) for g in args])
        except TypeError:
            raise ValueError(f"unable to find a common parent for the "
                             f"substituted species {args}")
        if not isinstance(P0, LazyCombinatorialSpecies):
            raise ValueError(f"the substituted species {args} must be "
                             f"ordinary lazy species")
        try:
            BR = cm.common_parent(fP.base_ring(), P0.base_ring())
        except TypeError:
            raise ValueError(f"unable to find a common base ring for {left} "
                            f"and the substituted species {args}")

        args = [P0(g) for g in args]
        P = LazyHyperoctahedralSpecies(
            BR, fP._laurent_poly_ring._r,
            P0._laurent_poly_ring._indices._indices._names,
            sparse=fP._sparse)

        R = P._internal_poly_ring.base_ring()

        # f = 0
        if isinstance(left._coeff_stream, Stream_zero):
            super().__init__(P, Stream_zero())
            self._left = left
            self._args = args
            return

        # all substituted species are zero, so the composition is constant
        if all(isinstance(g._coeff_stream, Stream_zero) for g in args):
            c = left[0]
            if c:
                c = list(c)[0][1]
            else:
                c = P.base_ring().zero()
            super().__init__(P, P(c)._coeff_stream)
            self._left = left
            self._args = args
            return

        # the outer species is a constant polynomial
        if (isinstance(left._coeff_stream, Stream_exact)
                and not left._coeff_stream._constant
                and left._coeff_stream._degree == 1):
            c = left._coeff_stream[0]
            B = c.parent()
            if not (B is ZZ or B is QQ or B == fP.base_ring()):
                c = c.coefficients()[0]
            super().__init__(P, P(c)._coeff_stream)
            self._left = left
            self._args = args
            return

        for g in args:
            if g._coeff_stream._approximate_order == 0:
                if not g._coeff_stream.is_uninitialized() and g[0]:
                    raise ValueError("can only compose with a positive valuation series")
                g._coeff_stream._approximate_order = 1

        sorder = left._coeff_stream._approximate_order
        gv = min(g._coeff_stream._approximate_order for g in args)
        L = fP._internal_poly_ring.base_ring()
        if L.base_ring() is BR:
            LR = L

            def lcoeff(c):
                return c
        else:
            # the homogeneous components of the outer species carry its
            # weights, so they have to be re-based to the common base
            # ring, in which also the weights of the args live
            LR = PolynomialHyperoctahedralSpecies(
                BR, fP._laurent_poly_ring._r,
                fP._laurent_poly_ring._indices._indices._names)

            def lcoeff(c):
                return BR(c)

        def coeff(g, i):
            c = g._coeff_stream[i]
            if not isinstance(c, PolynomialSpecies.Element):
                return R(c)
            return c

        # args_flat and weights contain one list for each substituted species
        weight_exp = [lazy_list(lambda j, g=g: len(coeff(g, j+1)))
                      for g in args]

        def flat(g):
            # function needed to work around python's scoping rules
            return itertools.chain.from_iterable(coeff(g, j) for j in itertools.count())

        args_flat1 = [lazy_list(flat(g)) for g in args]

        def coefficient(n):
            if not n:
                if left[0]:
                    return R(list(left[0])[0][1])
                return R.zero()
            result = R.zero()
            for i in range(1, n // gv + 1):
                # skip i=0 because it produces a term only for n=0

                # compute homogeneous components
                lF = defaultdict(LR)
                for M, c in left[i]:
                    lF[M.grade()] += LR._from_dict({M: lcoeff(c)})
                for mc, F in lF.items():
                    for degrees in weighted_vector_compositions(mc, n, weight_exp):
                        args_flat = [list(a[0:len(degrees[j])])
                                     for j, a in enumerate(args_flat1)]
                        multiplicities = [c for alpha, g_flat in zip(degrees, args_flat)
                                          for d, (_, c) in zip(alpha, g_flat) if d]
                        subs = [M for alpha, g_flat in zip(degrees, args_flat)
                                for d, (M, _) in zip(alpha, g_flat) if d]
                        non_zero_degrees = [[d for d in alpha if d] for alpha in degrees]
                        names = ["X%s" % i for i in range(len(subs))]
                        FX = F._compose_with_weighted_singletons(names,
                                                                 multiplicities,
                                                                 non_zero_degrees)
                        FG = [(M._type2_substitute_molecular(subs), c)
                              for M, c in FX]
                        result += R.sum_of_terms(FG)
            return result

        coeff_stream = Stream_function(coefficient, P._sparse, sorder * gv)
        super().__init__(P, coeff_stream)
        self._left = left
        self._args = args


class LazyHyperoctahedralSpecies(LazyCompletionGradedAlgebra):
    r"""
    The ring of lazy `r`-species.

    This is the completion of the graded algebra of polynomial
    `r`-species, where coefficients of degree `n` are finite linear
    combinations of molecular `r`-species of degree `n`, i.e., of
    `r`-species with `n` `C_r`-blocks.

    INPUT:

    - ``base_ring`` -- the base ring

    - ``r`` -- positive integer; the order of the cyclic group `C_r`

    - ``names`` -- names of the sorts

    - ``sparse`` -- boolean (default: ``True``); whether we use a sparse
      or a dense representation

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
        sage: LazyHyperoctahedralSpecies(QQ, 2)
        Lazy completion of Polynomial 2-species in X over Rational Field
        sage: LazyHyperoctahedralSpecies(QQ, 2, "X, Y")
        Lazy completion of Polynomial 2-species in X, Y over Rational Field

    TESTS::

        sage: LazyHyperoctahedralSpecies(QQ, 2, "X") is LazyHyperoctahedralSpecies(QQ, 2, "X")
        True
    """

    Element = LazyHyperoctahedralSpeciesElement

    @staticmethod
    def __classcall_private__(cls, base_ring, r, names="X", sparse=True):
        r"""
        Normalize input to ensure a unique representation.
        """
        from sage.structure.category_object import normalize_names
        names = normalize_names(-1, names)
        return super().__classcall__(cls, base_ring, ZZ(r), names, sparse)

    def _first_ngens(self, n):
        r"""
        Used by the preparser for ``F.<x> = ...``.

        We do not use the generic implementation of
        :class:`sage.combinat.CombinatorialFreeModule`, because we do
        not want to implement `gens`.

        Only the first ``n`` degree one species are returned::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L = LazyHyperoctahedralSpecies(QQ, 2, "X, Y")
            sage: L._first_ngens(2)
            (X, Y)
        """
        return tuple(self(g) for g in self._laurent_poly_ring._first_ngens(n))

    def __init__(self, base_ring, r, names, sparse):
        r"""
        Initialize the ring of lazy `r`-species.

        TESTS::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: LazyHyperoctahedralSpecies(QQ, 2)._arity
            1
            sage: LazyHyperoctahedralSpecies(QQ, 2, "X, Y")._arity
            2
        """
        super().__init__(PolynomialHyperoctahedralSpecies(base_ring, r, names),
                         sparse=sparse)
        self._arity = len(names)

    def SignedGraphs(self, connected=False):
        r"""
        Return the species of signed graphs.

        A signed graph is a simple graph in which every edge carries a
        label in `Zmod(r)`, thought of as a power of a primitive `r`-th
        root of unity.  A relabeling of a signed graph permutes the
        vertices, whereas a sign change in the free `C_r`-set of labels
        multiplies the labels of the adjacent edges.

        For `r = 1` this is the species of simple graphs.

        INPUT:

        - ``connected`` -- boolean; whether the graphs should be
          connected

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: L.SignedGraphs()[2]
            E_2(X°) + X°(E_2)
            sage: L.SignedGraphs(connected=True)[2]
            X°(E_2)

            sage: sorted(str(G.edges()) for G in L.SignedGraphs().isotypes(3))
            ['[(1, 2, 0), (1, 3, 0), (2, 3, 0)]', '[(1, 2, 0), (1, 3, 0), (2, 3, 1)]',
             '[(1, 3, 0), (2, 3, 0)]', '[(2, 3, 0)]', '[]']

        The species of signed graphs is not a type 1 substitution of an
        ordinary species, because the stabilizers of its structures are
        in general not wreath products.  We list the orders of the
        stabilizer groups of the molecular components, since the choice
        of their generators depends on the order of computation::

            sage: sorted((M.permutation_group()[0].order(), c) for M, c in L.SignedGraphs()[3])
            [(4, 1), (8, 1), (12, 2), (48, 1)]

        TESTS::

            sage: LazyHyperoctahedralSpecies(QQ, 2, "X, Y").SignedGraphs()
            Traceback (most recent call last):
            ...
            ValueError: the species of signed graphs is only implemented for a single sort
        """
        if self._arity != 1:
            raise ValueError("the species of signed graphs is only implemented for a single sort")
        return SignedGraphSpecies(self, connected=bool(connected))

    def SeidelGraphs(self):
        r"""
        Return the species of Seidel graphs.

        A Seidel graph is a complete graph whose edges carry labels in
        `Zmod(r)`.  A relabeling of a graph permutes the vertices,
        whereas a sign change in the free `C_r`-set of labels adds `1`
        to the labels of the adjacent edges.

        For `r = 2`, identifying the label `1` of a pair of vertices
        with the presence of an edge, this is the species of simple
        graphs together with the action of the hyperoctahedral group by
        relabelings and switchings.  Its isomorphism types are the
        two-graphs.

        Note that connected graphs do not form a subspecies, since
        switching does not preserve connectedness: switching at a
        vertex of the empty graph produces a star.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: L.SeidelGraphs()[2]
            X°(E_2)

            sage: sorted(str(G.edges()) for G in L.SeidelGraphs().isotypes(3))
            ['[(1, 2, 0), (1, 3, 0), (2, 3, 0)]', '[(1, 2, 0), (1, 3, 0), (2, 3, 1)]']

        TESTS::

            sage: LazyHyperoctahedralSpecies(QQ, 2, "X, Y").SeidelGraphs()
            Traceback (most recent call last):
            ...
            ValueError: the species of Seidel graphs is only implemented for a single sort
        """
        if self._arity != 1:
            raise ValueError("the species of Seidel graphs is only implemented for a single sort")
        return SeidelGraphSpecies(self)

    def _require_unisort_r2(self, name):
        r"""
        Raise a ``ValueError`` unless ``self`` is a ring of unisort
        lazy `C_2`-species.

        INPUT:

        - ``name`` -- string; the name of the species, used in the error
          message

        TESTS::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: LazyHyperoctahedralSpecies(QQ, 3).SignSpecies()
            Traceback (most recent call last):
            ...
            ValueError: the species of signs is only implemented for r = 2
            sage: LazyHyperoctahedralSpecies(QQ, 2, "X, Y").SignSpecies()
            Traceback (most recent call last):
            ...
            ValueError: the species of signs is only implemented for a single sort
        """
        if self._laurent_poly_ring._r != 2:
            raise ValueError("the species of %s is only implemented for r = 2" % name)
        if self._arity != 1:
            raise ValueError("the species of %s is only implemented for a single sort" % name)

    @cached_method
    def SignSpecies(self):
        r"""
        Return the species of signs, the sum `\sum_n \mathrm{sgn}_n`.

        A structure of degree `n` is a sign, which a signed permutation
        acts upon by multiplying with the parity of the length of a
        shortest reduced word, i.e., with the product of the sign of the
        induced permutation of the orbits and the number of orbits whose
        sign is flipped.  In particular, there are two structures on `n`
        orbits for every `n`, including the empty one.

        The homogeneous component of degree `n` is the molecular
        species whose structure stabilizer is the alternating group of
        signed permutations.  Its cycle index is `s_n(t) + s_{1^n}(u)`,
        where `t` and `u` are the virtual alphabets whose power sums
        are `p_k(t) = (p_k(ζ^0) + p_k(ζ^1))/2` and
        `p_k(u) = (p_k(ζ^0) - p_k(ζ^1))/2`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: S = L.SignSpecies()

        There are two structures at every degree, so that the
        generating series is `2 e^{z/2}`::

            sage: S.generating_series().truncate(5)
            2 + X + 1/4*X^2 + 1/24*X^3 + 1/192*X^4

        The sign of a signed permutation on a single orbit is the flip,
        so that the degree one component is the free singleton `X`; the
        even signed permutations on two orbits are the powers of a
        negative 2-cycle::

            sage: S[1] == X
            True
            sage: S[2] == L.NegativeCycles()[2]
            True

        The cycle index of the degree `n` component is `s_n(t) +
        s_{1^n}(u)`::

            sage: sum(c * M.cycle_index()
            ....:     for M, c in S[2].monomial_coefficients().items())
            1/4*p_{1,1}(ζ^1) + 1/2*p_2(ζ^1) + 1/4*p_{1,1}(ζ^0)

        TESTS::

            sage: TestSuite(S).run(skip=['_test_pickling'])
        """
        self._require_unisort_r2("signs")
        P = self._laurent_poly_ring

        def coefficient(n):
            if not n:
                return 2 * P.one()
            return P(_alternating_hyperoctahedral_group(n))
        return self(coefficient)

    @cached_method
    def ParitySpecies(self):
        r"""
        Return the species of parities, the sum `\sum_n \mathrm{prt}_n`.

        A structure of degree `n` is a parity, which a signed
        permutation acts upon by multiplying with the parity of the
        number of flipped orbits.  Equivalently, the structure `+`
        corresponds to the signed permutations with an even number of
        flipped orbits, and the structure `-` to those with an odd
        number.  In particular, there are two structures on `n` orbits
        for every `n`, including the empty one.

        The homogeneous component of degree `n` is the molecular
        species whose structure stabilizer is the group of signed
        permutations flipping an even number of orbits.  Its cycle
        index is `h_n(t) + h_n(u)`, where `t` and `u` are the virtual
        alphabets whose power sums are `p_k(t) = (p_k(ζ^0) +
        p_k(ζ^1))/2` and `p_k(u) = (p_k(ζ^0) - p_k(ζ^1))/2`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: Par = L.ParitySpecies()

            sage: Par.generating_series().truncate(5)
            2 + X + 1/4*X^2 + 1/24*X^3 + 1/192*X^4

        The parity of a signed permutation on a single orbit is the flip,
        so that the degree one component is the free singleton `X`; the
        signed permutations on two orbits flipping an even number of
        orbits stabilize a positive 2-cycle::

            sage: Par[1] == X
            True
            sage: Par[2] == L.PositiveCycles()[2]
            True

        The cycle index of the degree `n` component is `h_n(t) + h_n(u)`::

            sage: sum(c * M.cycle_index()
            ....:     for M, c in Par[2].monomial_coefficients().items())
            1/4*p_{1,1}(ζ^1) + 1/4*p_{1,1}(ζ^0) + 1/2*p_2(ζ^0)

        TESTS::

            sage: TestSuite(Par).run(skip=['_test_pickling'])
        """
        self._require_unisort_r2("parities")
        P = self._laurent_poly_ring

        def coefficient(n):
            if not n:
                return 2 * P.one()
            return P(_even_flip_group(n))
        return self(coefficient)

    @cached_method
    def PositiveCycles(self):
        r"""
        Return the species of positive cycles, the sum `\sum_n C^+_n`.

        A structure of degree `n` is a positive `n`-cycle, a cycle
        through all `n` orbits whose signs multiply to `+1`; the
        structures are in bijection with the signed permutations that
        consist of a single positive cycle.  The structure stabilizer
        is generated by the rotation of the cycle and the flip of all
        its orbits, and has order `2n`.

        The generating series is `-log(1-z)/2`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: C = L.PositiveCycles()

            sage: C.generating_series().truncate(5)
            1/2*X + 1/4*X^2 + 1/6*X^3 + 1/8*X^4

        A positive cycle on a single orbit has full cyclic stabilizer,
        so that the degree one component is the singleton `X^\circ`,
        and on two orbits the stabilizer coincides with the group of
        the parity species::

            sage: C[1] == L(_wreath_group(2, 1))
            True
            sage: C[2] == L.ParitySpecies()[2]
            True

        Positive and negative cycles on an odd number of orbits have
        conjugate stabilizers, so that the two species coincide there::

            sage: all(C[n] == L.NegativeCycles()[n] for n in range(1, 6) if n % 2)
            True
            sage: C[2] == L.NegativeCycles()[2]
            False

        TESTS::

            sage: TestSuite(C).run(skip=['_test_pickling'])
        """
        self._require_unisort_r2("positive cycles")
        P = self._laurent_poly_ring
        return self(lambda n: P(_positive_cycle_group(n)) if n else P.zero())

    @cached_method
    def NegativeCycles(self):
        r"""
        Return the species of negative cycles, the sum `\sum_n C^-_n`.

        A structure of degree `n` is a negative `n`-cycle, a cycle
        through all `n` orbits whose signs multiply to `-1`; the
        structures are in bijection with the signed permutations that
        consist of a single negative cycle.  The structure stabilizer
        is the cyclic group generated by a negative `n`-cycle, and has
        order `2n`.

        The generating series is `-log(1-z)/2`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: C = L.NegativeCycles()

            sage: C.generating_series().truncate(5)
            1/2*X + 1/4*X^2 + 1/6*X^3 + 1/8*X^4

        A negative cycle on a single orbit has full cyclic stabilizer,
        so that the degree one component is the singleton `X^\circ`;
        the even signed permutations on two orbits are the powers of a
        negative 2-cycle::

            sage: C[1] == L(_wreath_group(2, 1))
            True
            sage: C[2] == L.SignSpecies()[2]
            True

        Negative and positive cycles on an odd number of orbits have
        conjugate stabilizers, so that the two species coincide there::

            sage: all(C[n] == L.PositiveCycles()[n] for n in range(1, 6) if n % 2)
            True
            sage: C[2] == L.PositiveCycles()[2]
            False

        TESTS::

            sage: TestSuite(C).run(skip=['_test_pickling'])
        """
        self._require_unisort_r2("negative cycles")
        P = self._laurent_poly_ring
        return self(lambda n: P(_negative_cycle_group(n)) if n else P.zero())

    @cached_method
    def SignedCycles(self):
        r"""
        Return the species of signed cycles, the sum `\sum_n C^{\pm}_n`
        of the positive and the negative cycles.

        A structure of degree `n` is a signed permutation that consists
        of a single cycle through all `n` orbits.

        The cycle index of the homogeneous component of degree `n` is
        `\frac{1}{n} \sum_{k \mid n} \phi(k) (p_k^{n/k}(ζ^0) +
        p_k^{n/k}(ζ^1))`, and the generating series is `-log(1-z)`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.rings.lazy_species import LazyCombinatorialSpecies
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: C = L.SignedCycles()
            sage: C == L.PositiveCycles() + L.NegativeCycles()
            True

            sage: C.generating_series().truncate(5)
            X + 1/2*X^2 + 1/3*X^3 + 1/4*X^4

            sage: sum(c * M.cycle_index()
            ....:     for M, c in C[3].monomial_coefficients().items())
            1/3*p_{1,1,1}(ζ^1) + 2/3*p_3(ζ^1) + 1/3*p_{1,1,1}(ζ^0) + 2/3*p_3(ζ^0)

        A set of signed cycles is a signed permutation, so that the
        type 1 substitution of the signed cycles into the species of
        sets is the species of signed permutations, whose generating
        series is `1/(1-z)`::

            sage: E = LazyCombinatorialSpecies(QQ, "Z").Sets()
            sage: P = E(C)
            sage: [P.generating_series()[n] for n in range(5)]
            [1, 1, 1, 1, 1]

        The isomorphism types are the signed cycle types, so that their
        number is the number of pairs of partitions of total size `n`::

            sage: P.isotype_generating_series()[:5]
            [1, 2, 5, 10, 20]

        TESTS::

            sage: TestSuite(C).run(skip=['_test_pickling'])

        """
        return self.PositiveCycles() + self.NegativeCycles()

    @cached_method
    def OrientedCycles(self):
        r"""
        Return the species of oriented cycles, the sum `\sum_n C^o_n`.

        A structure of degree `n` is an oriented cycle on `n-1` of the
        `n` orbits, the remaining orbit being the missing label of the
        cycle.  A signed permutation acts by relabelling, except that
        changing the sign of the missing label reverses the cycle.

        There are `n (n-2)!` structures on `n` orbits, so that the
        generating series is `z/2 (1 - log(1-z/2))`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: O = L.OrientedCycles()

            sage: O.generating_series().truncate(5)
            1/2*X + 1/4*X^2 + 1/16*X^3 + 1/48*X^4

        An oriented cycle on a single orbit has full cyclic stabilizer,
        so that the degree one component is the singleton `X^\circ`,
        and an oriented cycle on two orbits is a pair of them::

            sage: O[1] == L(_wreath_group(2, 1))
            True
            sage: O[2] == L(_wreath_group(2, 1))^2
            True

        An oriented cycle on three orbits consists of the choice of the
        missing label together with an orientation of the cycle on the
        remaining two orbits::

            sage: len(list(O.structures([(1, 2), (3, 4), (5, 6)])))
            3

        TESTS::

            sage: TestSuite(O).run(skip=['_test_pickling'])
        """
        self._require_unisort_r2("oriented cycles")
        P = self._laurent_poly_ring
        return self(lambda n: P(_oriented_cycle_group(n)) if n else P.zero())


def _lift_permutation(sigma, r):
    r"""
    Return the lift of the permutation ``sigma`` to the wreath product.

    The permutation ``sigma`` of the blocks `1, \ldots, n` is lifted to
    the permutation of `\{1, \ldots, rn\}` mapping the `t`-th point of
    block `i` to the `t`-th point of block `sigma(i)`, for all `t` in
    `range(r)`.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _lift_permutation
        sage: from sage.groups.perm_gps.constructor import PermutationGroupElement
        sage: _lift_permutation(PermutationGroupElement("(1,3)"), 2)
        (1,5)(2,6)
    """
    cycles = []
    for cyc in sigma.cycle_tuples():
        for t in range(r):
            cycles.append(tuple((i - 1) * r + t + 1 for i in cyc))
    return PermutationGroupElement(cycles)


def _signed_permutation(eps, sigma):
    r"""
    Return the signed permutation with signs ``eps`` and underlying
    permutation ``sigma`` as a permutation of `\{1, \ldots, 2n\}`.

    The points `2i+1` and `2i+2` are the phases `0` and `1` of the
    orbit `i`, for `i` in `range(n)`.  The signed permutation
    `(\epsilon; \sigma)` maps the phase `t` of the orbit `i` to the
    phase `t + \epsilon_i` of the orbit `sigma(i)`.

    INPUT:

    - ``eps`` -- list of `n` integers; the signs of the orbits, taken
      modulo `2`

    - ``sigma`` -- list of `n` integers; the images `sigma[i]` of the
      orbits, i.e., a permutation of `range(n)` in one-line notation

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _signed_permutation
        sage: _signed_permutation([1, 0], [1, 0])   # a negative 2-cycle
        (1,4,2,3)
        sage: _signed_permutation([1, 1], [0, 1])   # the flip of all orbits
        (1,2)(3,4)
    """
    n = len(eps)
    images = []
    for i in range(n):
        j = sigma[i]
        for t in range(2):
            images.append(2*j + 1 + (t + eps[i]) % 2)
    return PermutationGroupElement(images)


def _alternating_hyperoctahedral_group(n):
    r"""
    Return the subgroup of even signed permutations of
    ``_wreath_group(2, n)``.

    This is the kernel of the sign character of `W(2, n)`, which maps
    `(\epsilon; \sigma)` to `sign(\sigma) \prod_i \epsilon_i`.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _alternating_hyperoctahedral_group
        sage: [_alternating_hyperoctahedral_group(n).order() for n in range(1, 5)]
        [1, 4, 24, 192]
    """
    W = _wreath_group(2, n)
    gens = []
    for k in range(2, n):
        sigma = list(range(n))
        sigma[0], sigma[1], sigma[k] = 1, k, 0
        gens.append(W(_signed_permutation([0] * n, sigma)))
    if n > 1:
        gens.append(W(_signed_permutation([1, 1] + [0] * (n - 2),
                                         list(range(n)))))
        sigma = list(range(n))
        sigma[0], sigma[1] = 1, 0
        gens.append(W(_signed_permutation([0, 1] + [0] * (n - 2), sigma)))
    return W.subgroup(gens)


def _even_flip_group(n):
    r"""
    Return the subgroup of ``_wreath_group(2, n)`` of signed
    permutations flipping an even number of orbits.

    This is the kernel of the character `(\epsilon; \sigma) \mapsto
    \prod_i \epsilon_i`, and also the intersection of
    ``_wreath_group(2, n)`` with the alternating group `A_{2n}`.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _even_flip_group
        sage: [_even_flip_group(n).order() for n in range(1, 5)]
        [1, 4, 24, 192]
    """
    W = _wreath_group(2, n)
    gens = []
    for i in range(n - 1):
        sigma = list(range(n))
        sigma[i], sigma[i + 1] = i + 1, i
        gens.append(W(_signed_permutation([0] * n, sigma)))
    if n > 1:
        gens.append(W(_signed_permutation([1, 1] + [0] * (n - 2),
                                         list(range(n)))))
    return W.subgroup(gens)


def _negative_cycle_group(n):
    r"""
    Return the cyclic subgroup of ``_wreath_group(2, n)`` generated by
    a negative `n`-cycle.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _negative_cycle_group
        sage: [_negative_cycle_group(n).order() for n in range(1, 5)]
        [2, 4, 6, 8]
    """
    sigma = [(i + 1) % n for i in range(n)]
    W = _wreath_group(2, n)
    return W.subgroup([W(_signed_permutation([1] + [0] * (n - 1), sigma))])


def _positive_cycle_group(n):
    r"""
    Return the subgroup of ``_wreath_group(2, n)`` stabilizing a
    positive `n`-cycle.

    This is generated by the positive `n`-cycle and the flip of all
    orbits.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _positive_cycle_group
        sage: [_positive_cycle_group(n).order() for n in range(1, 5)]
        [2, 4, 6, 8]
    """
    sigma = [(i + 1) % n for i in range(n)]
    W = _wreath_group(2, n)
    gens = [W(_signed_permutation([0] * n, sigma)),
            W(_signed_permutation([1] * n, list(range(n))))]
    return W.subgroup(gens)


def _oriented_cycle_group(n):
    r"""
    Return the stabilizer in ``_wreath_group(2, n)`` of an oriented
    cycle on `n-1` of the `n` orbits.

    The `n`-th orbit is the missing label of the cycle.  The stabilizer
    is generated by the rotations of the cycle, the flip of an orbit of
    the cycle, and the element flipping the missing label together with
    reversing the cycle.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _oriented_cycle_group
        sage: [_oriented_cycle_group(n).order() for n in range(1, 6)]
        [2, 4, 16, 48, 128]
    """
    m = n - 1
    rho = list(range(n))
    for i in range(m):
        rho[i] = (i + 1) % m
    rev = list(range(n))
    for i in range(m):
        rev[i] = m - 1 - i
    W = _wreath_group(2, n)
    gens = [W(_signed_permutation([0] * m + [1], rev))]
    if n > 1:
        gens.append(W(_signed_permutation([1] + [0] * (n - 1),
                                         list(range(n)))))
        gens.append(W(_signed_permutation([0] * n, rho)))
    return W.subgroup(gens)


def _labeling_orbits(WG, E, r):
    r"""
    Iterate over the orbits of the edge labelings of a graph.

    The labelings of the edges ``E`` of a graph on the vertices
    `1, \ldots, n` are acted upon by the subgroup ``WG`` of the wreath
    product `W(r, n) = C_r \wr S_n`: the induced permutation of the
    blocks relabels the vertices, whereas a sign change in a block `i`
    adds `\zeta` to the labels of the edges adjacent to vertex `i`.

    INPUT:

    - ``WG`` -- a subgroup of the wreath product `_wreath_group(r, n)`
      on the standard domain

    - ``E`` -- the sorted list of edges of a graph on the vertices
      `1, \ldots, n`, where `n` is the number of blocks of ``WG``

    - ``r`` -- positive integer; the order of the cyclic group `C_r`

    OUTPUT: pairs ``(ell, H)``, where ``ell`` is a tuple of edge labels
    in ``range(r)``, indexed by ``E``, and ``H`` is the stabilizer of
    the corresponding labeled graph, a subgroup of ``WG``.

    ALGORITHM:

    The orbits are enumerated using breadth first search with the
    generators of ``WG``, and the stabilizers are computed with GAP.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _labeling_orbits
        sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
        sage: W = _wreath_group(2, 3)
        sage: E = [(1, 2), (1, 3), (2, 3)]
        sage: [(ell, H.cardinality()) for ell, H in _labeling_orbits(W, E, 2)]
        [((0, 0, 0), 12), ((0, 0, 1), 12)]
    """
    n = WG.degree() // r
    index = {e: k for k, e in enumerate(E)}
    m = len(E)
    gens_WG = WG.gens()

    # for each generator precompute the induced permutation of the
    # edges together with the shifts of the edge labels: the label
    # of edge k is shifted by the sum of the signs at the endpoints
    # of the edge
    gen_maps = []
    for g in gens_WG:
        sigma = [0] * (n + 1)
        eps = [0] * (n + 1)
        for i in range(1, n + 1):
            p = g((i - 1) * r + 1)
            sigma[i] = (p - 1) // r + 1
            eps[i] = (p - 1) % r
        srcs = [0] * m
        shifts = [0] * m
        for k, (i, j) in enumerate(E):
            u, v = sigma[i], sigma[j]
            if u > v:
                u, v = v, u
            target = index[(u, v)]
            srcs[target] = k
            shifts[target] = (eps[i] + eps[j]) % r
        gen_maps.append((srcs, shifts))

    # iterate over the orbits of the edge labelings under the
    # subgroup, using breadth first search with its generators
    seen = set()
    for ell in itertools.product(range(r), repeat=m):
        if ell in seen:
            continue
        orbit = [ell]
        seen.add(ell)
        frontier = [ell]
        while frontier:
            x = frontier.pop()
            for srcs, shifts in gen_maps:
                y = tuple((x[k] + h) % r for k, h in zip(srcs, shifts))
                if y not in seen:
                    seen.add(y)
                    orbit.append(y)
                    frontier.append(y)
        # the stabilizer of the first element of the orbit, computed
        # as the stabilizer of the point 1 in the permutation
        # action of the generators on the orbit
        to_gap = {x: i for i, x in enumerate(orbit, 1)}
        perm_gens = [PermutationGroupElement([to_gap[tuple((x[k] + h) % r
                                                            for k, h in zip(srcs, shifts))]
                                              for x in orbit])
                     for srcs, shifts in gen_maps]
        OS = libgap.OrbitStabilizer(WG, 1, gens_WG, perm_gens)
        H = PermutationGroup(gap_group=OS["stabilizer"], domain=WG.domain())
        yield ell, H


def _signed_graph_orbits(n, r, connected=False):
    r"""
    Iterate over the orbits of the signed graphs with ``n`` vertices.

    The signed graphs are acted upon by the hyperoctahedral group
    `W(r, n) = C_r \wr S_n` on the domain `\{1, \ldots, rn\}`: the
    induced permutation of the blocks relabels the vertices, whereas a
    sign change in a block `i` multiplies the labels of the edges
    adjacent to vertex `i` with `\zeta`.

    INPUT:

    - ``n`` -- positive integer; the number of vertices

    - ``r`` -- positive integer; the order of the cyclic group `C_r`

    - ``connected`` -- boolean (default: ``False``); whether the
      underlying graphs should be connected

    OUTPUT: triples ``(G, ell, H)``, where ``G`` is a graph on the
    vertices `1, \ldots, n`, ``ell`` is a tuple of edge labels in
    ``range(r)``, indexed by the sorted edges of ``G``, and ``H`` is
    the stabilizer of the corresponding signed graph, a subgroup of
    `_wreath_group(r, n)`.

    ALGORITHM:

    We first iterate over the isomorphism classes of the underlying
    graphs.  Since the action of `W(r, n)` maps a signed graph to a
    signed graph with an isomorphic underlying graph, the orbits of
    signed graphs with underlying graph `G` are the orbits of the edge
    labelings of `G` under the subgroup `C_r^n \rtimes \operatorname{Aut}(G)`
    of `W(r, n)`.  These orbits are enumerated using breadth first
    search with the generators of this subgroup, and the stabilizers
    are computed with GAP.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _signed_graph_orbits
        sage: for G, ell, H in _signed_graph_orbits(2, 2):
        ....:     print(G.edges(labels=False), ell, H.cardinality())
        [] () 8
        [(1, 2)] (0,) 4
    """
    W = _wreath_group(r, n)
    # the generators of the cyclic groups C_r acting on the blocks
    rotations = [PermutationGroupElement(tuple(range(i * r + 1, i * r + r + 1)))
                 for i in range(n)]

    if connected:
        underlying = [G.canonical_label().relabel(range(1, n + 1), inplace=False)
                      for G in graphs.nauty_geng("%s -c" % n)]
    else:
        underlying = [G.canonical_label().relabel(range(1, n + 1), inplace=False)
                      for G in graphs(n)]

    for G in underlying:
        E = sorted((min(u, v), max(u, v))
                   for u, v in G.edge_iterator(labels=False))
        # the subgroup of the wreath product preserving the underlying
        # graph: the lifts of its automorphisms together with the
        # rotations of the blocks
        A = G.automorphism_group()
        gens = [_lift_permutation(sigma, r) for sigma in A.gens()]
        gens.extend(rotations)
        WG = W.subgroup(gens)

        for ell, H in _labeling_orbits(WG, E, r):
            yield G, ell, H


def _seidel_graph_orbits(n, r):
    r"""
    Iterate over the switching orbits of the graphs with ``n`` vertices.

    A graph is a labeling of the pairs of `\{1, \ldots, n\}` with
    labels in ``range(r)``.  The graphs are acted upon by the
    hyperoctahedral group `W(r, n) = C_r \wr S_n`: the induced
    permutation of the blocks relabels the vertices, whereas a sign
    change in a block `i` adds `1` to the labels of the pairs adjacent
    to vertex `i`.

    For `r = 2`, identifying the label `1` of a pair with the presence
    of an edge, this is the action of `W(2, n)` on the simple graphs
    with `n` vertices by relabelings and switchings.

    INPUT:

    - ``n`` -- positive integer; the number of vertices

    - ``r`` -- positive integer; the order of the cyclic group `C_r`

    OUTPUT: pairs ``(ell, H)``, where ``ell`` is a tuple of labels in
    ``range(r)``, indexed by the sorted pairs of `\{1, \ldots, n\}`, and
    ``H`` is the stabilizer of the corresponding graph, a subgroup of
    `_wreath_group(r, n)`.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import _seidel_graph_orbits
        sage: [(ell, H.cardinality()) for ell, H in _seidel_graph_orbits(3, 2)]
        [((0, 0, 0), 12), ((0, 0, 1), 12)]
    """
    W = _wreath_group(r, n)
    E = sorted(itertools.combinations(range(1, n + 1), 2))
    yield from _labeling_orbits(W, E, r)


class SignedGraphSpecies(LazyHyperoctahedralSpeciesElement, UniqueRepresentation,
                         metaclass=InheritComparisonClasscallMetaclass):
    r"""
    The species of signed graphs.

    A signed graph is a simple graph in which every edge carries a
    label in `Zmod(r)`, thought of as a power of a primitive `r`-th
    root of unity.  The `r`-species of signed graphs assigns to a free
    `C_r`-set the set of all signed graphs whose vertices are its
    `C_r`-orbits.  A relabeling permutes the vertices, whereas a sign
    change multiplies the labels of the adjacent edges.

    Since the stabilizer of a signed graph is in general not a wreath
    product, this is an example of an `r`-species which is not a type 1
    substitution of an ordinary species.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
        sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
        sage: S = L.SignedGraphs()
        sage: S[:3]
        [1,
         X°,
         E_2(X°) + X°(E_2)]

        The orders of the stabilizer groups of the molecular components
        of a homogeneous component do not depend on the order in which
        the components are computed, in contrast to the choice of their
        generators::

        sage: [sorted((M.permutation_group()[0].order(), c) for M, c in S[n])
        ....:  for n in [3, 4]]
        [[(4, 1), (8, 1), (12, 2), (48, 1)],
         [(4, 1), (4, 1), (4, 2), (8, 1), (8, 3), (12, 1), (16, 1), (16, 1), (24, 2), (32, 1), (32, 1), (48, 2), (384, 1)]]

    The isomorphism types of signed graphs are signed graphs with
    vertices `1, \ldots, n`::

        sage: sorted(str(G.edges()) for G in S.isotypes(3))
        ['[(1, 2, 0), (1, 3, 0), (2, 3, 0)]', '[(1, 2, 0), (1, 3, 0), (2, 3, 1)]',
         '[(1, 3, 0), (2, 3, 0)]', '[(2, 3, 0)]', '[]']

    Every signed graph decomposes uniquely into connected signed graphs::

        sage: E = LazyCombinatorialSpecies(QQ, "Z").Sets()
        sage: all(E(L.SignedGraphs(connected=True))[n] == S[n] for n in range(5))
        True

    There are `(r+1)^{\binom{n}{2}}` signed graphs with `n` vertices, so
    the generating series has a closed form::

        sage: S.generating_series().truncate(7)
        1 + 1/2*X + 3/8*X^2 + 9/16*X^3 + 243/128*X^4 + 19683/1280*X^5 + 1594323/5120*X^6
        sage: L.SignedGraphs(connected=True).generating_series().truncate(7)
        1/2*X + 1/4*X^2 + 5/12*X^3 + 13/8*X^4 + 1151/80*X^5 + 9103/30*X^6

    The number of isomorphism types of signed graphs with `n` vertices
    is the number of switching classes of signed graphs::

        sage: S.isotype_generating_series()[:6]
        [1, 1, 2, 5, 18, 100]

    For `r = 1` we recover the species of simple graphs::

        sage: L1 = LazyHyperoctahedralSpecies(QQ, 1)
        sage: Gs = L1.SignedGraphs()
        sage: G = LazyCombinatorialSpecies(QQ, "X").Graphs()
        sage: all(Gs.isotype_generating_series()[n] == G.isotype_generating_series()[n]
        ....:      for n in range(6))
        True

        sage: Gs[3]
        2*E_3(X) + 2*X*E_2(X)

    TESTS::

        sage: TestSuite(S).run(skip=['_test_category', '_test_pickling'])
        sage: TestSuite(L.SignedGraphs(connected=True)).run(skip=['_test_category', '_test_pickling'])
    """
    def __init__(self, parent, connected=False):
        r"""
        Initialize the species of signed graphs.

        INPUT:

        - ``parent`` -- a lazy species ring

        - ``connected`` -- boolean; whether the graphs should be
          connected

        TESTS::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L = LazyHyperoctahedralSpecies(QQ, 2)
            sage: S = L.SignedGraphs()
            sage: TestSuite(S).run(skip=['_test_category', '_test_pickling'])

            sage: S is L.SignedGraphs()
            True

            sage: S == L.SignedGraphs(connected=True)
            False

            sage: L.SignedGraphs(True) is L.SignedGraphs(connected=1)
            True
        """
        P = parent._laurent_poly_ring
        self._connected = connected = bool(connected)

        def coefficient(n):
            if not n:
                return P.one() if not connected else P.zero()
            return sum(P(H) for _, _, H
                       in _signed_graph_orbits(n, P._r, connected))

        S = parent(coefficient)
        super().__init__(parent, S._coeff_stream)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L = LazyHyperoctahedralSpecies(QQ, 2)
            sage: L.SignedGraphs()  # indirect doctest
            Signed graph species

            sage: L.SignedGraphs(connected=True)
            Connected signed graph species
        """
        if self._connected:
            return "Connected signed graph species"
        return "Signed graph species"

    def isotypes(self, labels):
        r"""
        Iterate over the isomorphism types of signed graphs with the
        given number of vertices.

        The isomorphism types are signed graphs with vertices
        `1, \ldots, n`, whose edge labels are in `Zmod(r)`.

        INPUT:

        - ``labels`` -- the number of vertices

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L = LazyHyperoctahedralSpecies(QQ, 2)
            sage: sorted(str(G.edges()) for G in L.SignedGraphs().isotypes(2))
            ['[(1, 2, 0)]', '[]']
            sage: sorted(str(G.edges()) for G in L.SignedGraphs(connected=True).isotypes(2))
            ['[(1, 2, 0)]']

            sage: list(L.SignedGraphs(connected=True).isotypes(0))
            []
            sage: list(L.SignedGraphs().isotypes(0))
            [Graph on 0 vertices]
        """
        if labels not in ZZ:
            raise NotImplementedError("isotypes with given labels are currently not supported")
        if not labels:
            if self._connected:
                return
            yield Graph([], immutable=True)
            return
        r = self.parent()._laurent_poly_ring._r
        for G, ell, _ in _signed_graph_orbits(labels, r, self._connected):
            E = sorted((min(u, v), max(u, v))
                       for u, v in G.edge_iterator(labels=False))
            result = Graph([(u, v, Zmod(r)(k)) for (u, v), k in zip(E, ell)])
            result.add_vertices(range(1, labels + 1))
            yield result.copy(immutable=True)

    def generating_series(self):
        r"""
        Return the generating series of the species of signed graphs.

        There are `(r+1)^{\binom{n}{2}}` signed graphs with `n` vertices,
        so the coefficient of `X^n` is
        `\frac{(r+1)^{\binom{n}{2}}}{n! r^n}`.

        The generating series of the species of connected signed graphs
        is its logarithm.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: L.SignedGraphs().generating_series().truncate(7)
            1 + 1/2*X + 3/8*X^2 + 9/16*X^3 + 243/128*X^4 + 19683/1280*X^5 + 1594323/5120*X^6

            sage: L.SignedGraphs(connected=True).generating_series().truncate(7)
            1/2*X + 1/4*X^2 + 5/12*X^3 + 13/8*X^4 + 1151/80*X^5 + 9103/30*X^6

            sage: LazyHyperoctahedralSpecies(QQ, 3).SignedGraphs().generating_series().truncate(5)
            1 + 1/3*X + 2/9*X^2 + 32/81*X^3 + 512/243*X^4

        TESTS::

            sage: L.SignedGraphs(connected=True).generating_series().exp().truncate(7) == L.SignedGraphs().generating_series().truncate(7)
            True
        """
        P = self.parent()
        r = P._laurent_poly_ring._r
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        s = L(lambda n: ZZ(r + 1) ** ZZ(n).binomial(2) / ZZ(n).factorial() / ZZ(r) ** n)
        if self._connected:
            return s.log()
        return s


class SeidelGraphSpecies(LazyHyperoctahedralSpeciesElement, UniqueRepresentation,
                         metaclass=InheritComparisonClasscallMetaclass):
    r"""
    The species of Seidel graphs.

    A Seidel graph is a complete graph in which every edge carries a
    label in `Zmod(r)`.  The `r`-species of Seidel graphs assigns to a
    free `C_r`-set the set of all graphs whose vertices are its
    `C_r`-orbits.  A relabeling permutes the vertices, whereas a sign
    change adds `1` to the labels of the adjacent edges.

    For `r = 2`, identifying the label `1` of a pair of vertices with
    the presence of an edge, this is the species of simple graphs
    together with the action of the hyperoctahedral group `W(2, n)` by
    relabelings and switchings: a switch at a vertex complements the
    edges adjacent to it.  The isomorphism types are the switching
    classes of graphs, i.e., the two-graphs.

    Since switching does not preserve connectedness, connected graphs do
    not form a subspecies.  Since the stabilizer of a graph is in
    general not a wreath product, this is an example of an `r`-species
    which is not a type 1 substitution of an ordinary species.

    EXAMPLES::

        sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
        sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
        sage: G = L.SeidelGraphs()
        sage: G[:3]
        [1, X°, X°(E_2)]

        The orders of the stabilizer groups of the molecular components
        of a homogeneous component do not depend on the order in which
        the components are computed, in contrast to the choice of their
        generators::

        sage: [sorted((M.permutation_group()[0].order(), c) for M, c in G[n])
        ....:  for n in [3, 4]]
        [[(12, 2)], [(8, 1), (48, 2)]]

    The isomorphism types are representatives of the switching classes
    of simple graphs, i.e., of the two-graphs.  They are complete graphs
    with vertices `1, \ldots, n` whose edges carry labels in `Zmod(r)`;
    for `r = 2` the label `1` of a pair of vertices indicates the
    presence of an edge::

        sage: sorted(str(g.edges()) for g in G.isotypes(3))
        ['[(1, 2, 0), (1, 3, 0), (2, 3, 0)]', '[(1, 2, 0), (1, 3, 0), (2, 3, 1)]']

    There are `r^{\binom{n}{2}}` graphs with `n` vertices, so the
    generating series has a closed form::

        sage: G.generating_series().truncate(7)
        1 + 1/2*X + 1/4*X^2 + 1/6*X^3 + 1/6*X^4 + 4/15*X^5 + 32/45*X^6

        sage: G.isotype_generating_series()[:7]
        [1, 1, 1, 2, 3, 7, 16]

    For `r = 1` this is the species of sets::

        sage: LazyHyperoctahedralSpecies(QQ, 1).SeidelGraphs()[:4]
        [1, X, E_2(X), E_3(X)]

    TESTS::

        sage: TestSuite(G).run(skip=['_test_category', '_test_pickling'])
    """
    def __init__(self, parent):
        r"""
        Initialize the species of Seidel graphs.

        INPUT:

        - ``parent`` -- a lazy species ring

        TESTS::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L = LazyHyperoctahedralSpecies(QQ, 2)
            sage: G = L.SeidelGraphs()
            sage: TestSuite(G).run(skip=['_test_category', '_test_pickling'])

            sage: G is L.SeidelGraphs()
            True
        """
        P = parent._laurent_poly_ring

        def coefficient(n):
            if not n:
                return P.one()
            return sum(P(H) for _, H in _seidel_graph_orbits(n, P._r))

        S = parent(coefficient)
        super().__init__(parent, S._coeff_stream)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: LazyHyperoctahedralSpecies(QQ, 2).SeidelGraphs()  # indirect doctest
            Seidel graph species
        """
        return "Seidel graph species"

    def isotypes(self, labels):
        r"""
        Iterate over the isomorphism types of Seidel graphs with the
        given number of vertices.

        The isomorphism types are complete graphs with vertices
        `1, \ldots, n`, whose edge labels are in `Zmod(r)`.  For `r = 2`,
        identifying the label `1` of a pair of vertices with the
        presence of an edge, these are representatives of the switching
        classes of simple graphs, i.e., of the two-graphs.

        INPUT:

        - ``labels`` -- the number of vertices

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L = LazyHyperoctahedralSpecies(QQ, 2)
            sage: sorted(str(G.edges()) for G in L.SeidelGraphs().isotypes(2))
            ['[(1, 2, 0)]']
            sage: sorted(str(G.edges()) for G in L.SeidelGraphs().isotypes(4))
            ['[(1, 2, 0), (1, 3, 0), (1, 4, 0), (2, 3, 0), (2, 4, 0), (3, 4, 0)]',
             '[(1, 2, 0), (1, 3, 0), (1, 4, 0), (2, 3, 0), (2, 4, 0), (3, 4, 1)]',
             '[(1, 2, 0), (1, 3, 0), (1, 4, 0), (2, 3, 1), (2, 4, 1), (3, 4, 1)]']

            sage: list(L.SeidelGraphs().isotypes(0))
            [Graph on 0 vertices]
        """
        if labels not in ZZ:
            raise NotImplementedError("isotypes with given labels are currently not supported")
        if not labels:
            yield Graph([], immutable=True)
            return
        r = self.parent()._laurent_poly_ring._r
        E = sorted(itertools.combinations(range(1, labels + 1), 2))
        for ell, _ in _seidel_graph_orbits(labels, r):
            result = Graph([(u, v, Zmod(r)(k)) for (u, v), k in zip(E, ell)])
            result.add_vertices(range(1, labels + 1))
            yield result.copy(immutable=True)

    def generating_series(self):
        r"""
        Return the generating series of the species of Seidel graphs.

        There are `r^{\binom{n}{2}}` graphs with `n` vertices, so the
        coefficient of `X^n` is `\frac{r^{\binom{n}{2}}}{n! r^n}`.

        EXAMPLES::

            sage: from sage.rings.lazy_species_hyperoctahedral import LazyHyperoctahedralSpecies
            sage: L.<X> = LazyHyperoctahedralSpecies(QQ, 2)
            sage: L.SeidelGraphs().generating_series().truncate(7)
            1 + 1/2*X + 1/4*X^2 + 1/6*X^3 + 1/6*X^4 + 4/15*X^5 + 32/45*X^6

            sage: LazyHyperoctahedralSpecies(QQ, 3).SeidelGraphs().generating_series().truncate(5)
            1 + 1/3*X + 1/6*X^2 + 1/6*X^3 + 3/8*X^4
        """
        P = self.parent()
        r = P._laurent_poly_ring._r
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        return L(lambda n: ZZ(r) ** ZZ(n).binomial(2) / ZZ(n).factorial() / ZZ(r) ** n)
