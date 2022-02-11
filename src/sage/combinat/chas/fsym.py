r"""
Poirier-Reutenauer Hopf algebra of standard tableaux

AUTHORS:

- Franco Saliola (2012): initial implementation
- Travis Scrimshaw (2018-04-11): added missing doctests and reorganization
"""

# ****************************************************************************
#       Copyright (C) 2012 Franco Saliola   <saliola at gmail.com>
#                     2018 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.bindable_class import BindableClass
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.realizations import Category_realization_of_parent
from sage.categories.hopf_algebras import HopfAlgebras

from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.tableau import Tableau, StandardTableaux
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.composition import Composition


class FSymBasis_abstract(CombinatorialFreeModule, BindableClass):
    r"""
    Abstract base class for graded bases of `FSym` and of `FSym^*`
    indexed by standard tableaux.

    This must define the following attributes:

    - ``_prefix`` -- the basis prefix
    """
    def __init__(self, alg, graded=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = algebras.FSym(QQ).G()
            sage: TestSuite(G).run()  # long time

        Checks for the antipode::

            sage: FSym = algebras.FSym(QQ)
            sage: G = FSym.G()
            sage: for b in G.basis(degree=3):
            ....:     print("%s : %s" % (b, b.antipode()))
            G[123] : -G[1|2|3]
            G[13|2] : -G[13|2]
            G[12|3] : -G[12|3]
            G[1|2|3] : -G[123]

            sage: F = FSym.dual().F()
            sage: for b in F.basis(degree=3):
            ....:     print("%s : %s" % (b, b.antipode()))
            F[123] : -F[1|2|3]
            F[13|2] : -F[13|2]
            F[12|3] : -F[12|3]
            F[1|2|3] : -F[123]
        """
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         StandardTableaux(),
                                         category=FSymBases(alg),
                                         bracket="", prefix=self._prefix)

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - elements of the algebra `FSym` over a base ring with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: G = algebras.FSym(GF(7)).G(); G
            Hopf algebra of standard tableaux over the Finite Field of size 7
             in the Fundamental basis

        Elements of `FSym` canonically coerce in::

            sage: x, y = G([[1]]), G([[1,3],[2]])
            sage: G.coerce(x + y) == x + y
            True

        Elements of `FSym` over `\ZZ` coerce in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: H = algebras.FSym(ZZ).G()
            sage: Hx, Hy = H([[1]]), H([[1,3],[2]])
            sage: z = G.coerce(Hx+Hy); z
            G[1] + G[13|2]
            sage: z.parent() is G
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the
        elements of `FSym` over `\GF{7}` do not coerce
        to the same algebra over `\ZZ`::

            sage: H.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Hopf algebra of standard tableaux
             over the Finite Field of size 7 in the Fundamental basis
             to Hopf algebra of standard tableaux over the Integer Ring in the Fundamental basis

        TESTS::

            sage: G = algebras.FSym(ZZ).G()
            sage: H = algebras.FSym(QQ).G()
            sage: G.has_coerce_map_from(H)
            False
            sage: H.has_coerce_map_from(G)
            True
            sage: G.has_coerce_map_from(QQ)
            False
            sage: H.has_coerce_map_from(QQ)
            True
            sage: G.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
            sage: G = algebras.FSym(QQ).G()
            sage: TG = G.dual_basis()
            sage: TG.coerce_map_from(G) is None
            True
            sage: G.coerce_map_from(TG) is None
            True
        """
        # free symmetric functions in the same variables
        # over any base that coerces in:
        if isinstance(R, FSymBasis_abstract):
            FSym = self.realization_of()
            if R.realization_of() == FSym:
                return True
            if (isinstance(R.realization_of(), FreeSymmetricFunctions) !=
                    isinstance(FSym, FreeSymmetricFunctions)):
                # If they are dual bases, then no coercion
                return False
            if not self.base_ring().has_coerce_map_from(R.base_ring()):
                return False
            if self._realization_name() == R._realization_name():
                # The same basis

                def coerce_base_ring(self, x):
                    return self._from_dict(x.monomial_coefficients())
                return coerce_base_ring
            # Otherwise lift that basis up and then coerce over
            target = getattr(FSym, R._realization_name())()
            return self._coerce_map_via([target], R)
        return super(FSymBasis_abstract, self)._coerce_map_from_(R)

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: G = algebras.FSym(QQ).G()
            sage: G.some_elements()
            [G[], G[1], G[12], G[1] + G[1|2], G[] + 1/2*G[1]]
        """
        u = self.one()
        o = self([[1]])
        s = self.base_ring().an_element()
        return [u, o, self([[1, 2]]), o + self([[1], [2]]), u + s * o]

    def _repr_term(self, phi):
        r"""
        The string representation of a basis element.

        EXAMPLES:

        We use a compact notation for standard tableaux::

            sage: FSym = algebras.FSym(QQ)
            sage: G = FSym.G()
            sage: G.zero()
            0
            sage: G._repr_term(StandardTableau([[1],[2],[3],[4]]))
            'G[1|2|3|4]'
            sage: G[[1,3,5],[2,4]]
            G[135|24]
        """
        return "{}[{}]".format(self._prefix,
                               "|".join("".join(map(str, block))
                                        for block in phi))


class FSymBases(Category_realization_of_parent):
    r"""
    The category of graded bases of `FSym` and `FSym^*` indexed
    by standard tableaux.
    """
    def super_categories(self):
        """
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.fsym import FSymBases
            sage: FSym = algebras.FSym(ZZ)
            sage: bases = FSymBases(FSym)
            sage: bases.super_categories()
            [Category of realizations of Hopf algebra of standard tableaux over the Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring
                 and Category of graded coalgebras over Integer Ring,
             Category of graded connected hopf algebras with basis over Integer Ring]
        """
        R = self.base().base_ring()
        return [self.base().Realizations(),
                HopfAlgebras(R).Graded().Realizations(),
                HopfAlgebras(R).Graded().WithBasis().Graded().Connected()]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `FSym`.

            EXAMPLES::

                sage: FSym = algebras.FSym(ZZ)
                sage: FSym.G()
                Hopf algebra of standard tableaux over the Integer Ring
                 in the Fundamental basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._realization_name())

        def __getitem__(self, key):
            r"""
            Override the ``__getitem__`` method to allow passing a standard
            tableau in a nonstandard form (e.g., as a tuple of rows instead
            of a list of rows; or as a single row for a single-rowed tableau).

            EXAMPLES:

            Construct the basis element indexed by a standard tableau by
            passing data that defines the standard tableau::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: G[[1,3],[2]]
                G[13|2]
                sage: G[(1,3),(2,)]
                G[13|2]
                sage: G[[1,3],[2]].leading_support() in StandardTableaux()
                True
                sage: G[1,2,3]
                G[123]
            """
            try:
                return self.monomial(self._indices(list(key)))
            except (TypeError, ValueError):
                return self.monomial(self._indices([key]))

        def basis(self, degree=None):
            r"""
            The basis elements (optionally: of the specified degree).

            OUTPUT: Family

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: TG = FSym.G()
                sage: TG.basis()
                Lazy family (Term map from Standard tableaux to Hopf algebra of standard tableaux
                 over the Rational Field in the Fundamental basis(i))_{i in Standard tableaux}
                sage: TG.basis().keys()
                Standard tableaux
                sage: TG.basis(degree=3).keys()
                Standard tableaux of size 3
                sage: TG.basis(degree=3).list()
                [G[123], G[13|2], G[12|3], G[1|2|3]]
            """
            from sage.sets.family import Family
            if degree is None:
                return Family(self._indices, self.monomial)
            else:
                return Family(StandardTableaux(degree), self.monomial)

        @cached_method
        def one_basis(self):
            r"""
            Return the basis index corresponding to `1`.

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: TG = FSym.G()
                sage: TG.one_basis()
                []
            """
            return self._indices([])

        def duality_pairing(self, x, y):
            r"""
            The canonical pairing between `FSym` and `FSym^*`.

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: F = G.dual_basis()
                sage: t1 = StandardTableau([[1,3,5],[2,4]])
                sage: t2 = StandardTableau([[1,3],[2,5],[4]])
                sage: G.duality_pairing(G[t1], F[t2])
                0
                sage: G.duality_pairing(G[t1], F[t1])
                1
                sage: G.duality_pairing(G[t2], F[t2])
                1
                sage: F.duality_pairing(F[t2], G[t2])
                1

                sage: z = G[[1,3,5],[2,4]]
                sage: all(F.duality_pairing(F[p1] * F[p2], z) == c
                ....:     for ((p1, p2), c) in z.coproduct())
                True

            TESTS:

            If ``x`` is zero, then the output still has the right
            type::

                sage: z = G.duality_pairing(G.zero(), F.zero()); z
                0
                sage: parent(z)
                Rational Field
            """
            y = self.dual_basis()(y)
            return self.base_ring().sum(coeff * y[t] for (t, coeff) in x)

        def duality_pairing_matrix(self, basis, degree):
            r"""
            The matrix of scalar products between elements of `FSym` and
            elements of `FSym^*`.

            INPUT:

            - ``basis`` -- a basis of the dual Hopf algebra
            - ``degree`` -- a non-negative integer

            OUTPUT:

            - the matrix of scalar products between the basis ``self`` and the
              basis ``basis`` in the dual Hopf algebra of degree ``degree``

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: G.duality_pairing_matrix(G.dual_basis(), 3)
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]
            """
            from sage.matrix.constructor import matrix
            keys = self.basis(degree=degree).keys()
            return matrix(self.base_ring(),
                          [[self.duality_pairing(self[s], basis[t])
                            for t in keys] for s in keys])

        def degree_on_basis(self, t):
            """
            Return the degree of a standard tableau in the algebra
            of free symmetric functions.

            This is the size of the tableau ``t``.

            EXAMPLES::

                sage: G = algebras.FSym(QQ).G()
                sage: t = StandardTableau([[1,3],[2]])
                sage: G.degree_on_basis(t)
                3
                sage: u = StandardTableau([[1,3,4,5],[2]])
                sage: G.degree_on_basis(u)
                5
            """
            return t.size()

    class ElementMethods:
        def duality_pairing(self, other):
            r"""
            Compute the pairing between ``self`` and an element ``other``
            of the dual.

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: F = G.dual_basis()
                sage: elt = G[[1,3],[2]] - 3*G[[1,2],[3]]
                sage: elt.duality_pairing(F[[1,3],[2]])
                1
                sage: elt.duality_pairing(F[[1,2],[3]])
                -3
                sage: elt.duality_pairing(F[[1,2]])
                0
            """
            return self.parent().duality_pairing(self, other)


class FreeSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The free symmetric functions.

    The *free symmetric functions* is a combinatorial Hopf algebra
    defined using tableaux and denoted `FSym`.

    Consider the Hopf algebra `FQSym`
    (:class:`~sage.combinat.fqsym.FreeQuasisymmetricFunctions`)
    over a commutative ring `R`, and its bases `(F_w)` and `(G_w)`
    (where `w`, in both cases, ranges over all permutations in all
    symmetric groups `S_0, S_1, S_2, \ldots`).
    For each word `w`, let `P(w)` be the P-tableau of `w` (that
    is, the first of the two tableaux obtained by applying the
    RSK algorithm to `w`; see :meth:`~sage.combinat.rsk.RSK`).
    If `t` is a standard tableau of size `n`, then we define
    `\mathcal{G}_t \in FQSym` to be the sum of the `F_w` with
    `w` ranging over all permutations of `\{1, 2, \ldots, n\}`
    satisfying `P(w) = t`. Equivalently, `\mathcal{G}_t` is the
    sum of the `G_w` with `w` ranging over all permutations of
    `\{1, 2, \ldots, n\}` satisfying `Q(w) = t` (where `Q(w)`
    denotes the Q-tableau of `w`).

    The `R`-linear span of the `\mathcal{G}_t` (for `t` ranging
    over all standard tableaux) is a Hopf subalgebra of `FQSym`,
    denoted by `FSym` and known as the *free symmetric functions*
    or the *Poirier-Reutenauer Hopf algebra of tableaux*. It has been
    introduced in [PoiReu95]_, where it was denoted by
    `(\ZZ T, \ast, \delta)`. (What we call `\mathcal{G}_t`
    has just been called `t` in [PoiReu95]_.)
    The family `(\mathcal{G}_t)` (with `t` ranging over all standard
    tableaux) is a basis of `FSym`, called the *Fundamental basis*.

    EXAMPLES:

    As explained above, `FSym` is constructed as a Hopf subalgebra of
    `FQSym`::

        sage: G = algebras.FSym(QQ).G()
        sage: F = algebras.FQSym(QQ).F()
        sage: G[[1,3],[2]]
        G[13|2]
        sage: G[[1,3],[2]].to_fqsym()
        G[2, 1, 3] + G[3, 1, 2]
        sage: F(G[[1,3],[2]])
        F[2, 1, 3] + F[2, 3, 1]

    This embedding is a Hopf algebra morphism::

        sage: all(F(G[t1] * G[t2]) == F(G[t1]) * F(G[t2])
        ....:     for t1 in StandardTableaux(2)
        ....:     for t2 in StandardTableaux(3))
        True

        sage: FF = F.tensor_square()
        sage: all(FF(G[t].coproduct()) == F(G[t]).coproduct()
        ....:     for t in StandardTableaux(4))
        True

    There is a Hopf algebra map from `FSym` onto the Hopf algebra
    of symmetric functions, which maps a tableau `t` to the Schur
    function indexed by the shape of `t`::

        sage: TG = algebras.FSym(QQ).G()
        sage: t = StandardTableau([[1,3],[2,4],[5]])
        sage: TG[t]
        G[13|24|5]
        sage: TG[t].to_symmetric_function()
        s[2, 2, 1]
    """
    def __init__(self, base_ring):
        r"""
        TESTS::

            sage: FSym = algebras.FSym(QQ)
            sage: TestSuite(FSym).run()
        """
        cat = HopfAlgebras(base_ring).Graded().Connected()
        Parent.__init__(self, base=base_ring, category=cat.WithRealizations())

    _shorthands = ['G']

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the Fundamental basis).

        EXAMPLES::

            sage: FSym = algebras.FSym(QQ)
            sage: FSym.a_realization()
            Hopf algebra of standard tableaux over the Rational Field
             in the Fundamental basis
        """
        return self.Fundamental()

    def dual(self):
        r"""
        Return the dual Hopf algebra of `FSym`.

        EXAMPLES::

            sage: algebras.FSym(QQ).dual()
            Dual Hopf algebra of standard tableaux over the Rational Field
        """
        return FreeSymmetricFunctions_Dual(self.base_ring())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.FSym(QQ)
            Hopf algebra of standard tableaux over the Rational Field
        """
        return "Hopf algebra of standard tableaux over the %s" % self.base_ring()

    class Fundamental(FSymBasis_abstract):
        r"""
        The Hopf algebra of tableaux on the Fundamental basis.

        EXAMPLES::

            sage: FSym = algebras.FSym(QQ)
            sage: TG = FSym.G()
            sage: TG
            Hopf algebra of standard tableaux over the Rational Field
             in the Fundamental basis

        Elements of the algebra look like::

            sage: TG.an_element()
            2*G[] + 2*G[1] + 3*G[12]

        TESTS::

            sage: FSym = algebras.FSym(QQ)
            sage: TG = FSym.G()
            sage: TestSuite(TG).run()
        """
        _prefix = "G"

        def _coerce_map_from_(self, R):
            r"""
            Return ``True`` if there is a coercion from ``R`` into ``self``
            and ``False`` otherwise.

            The things that coerce into ``self`` are

            - elements of the algebra `FSym` over a base ring
              with a coercion map into ``self.base_ring()``
            - non-commutative symmetric functions over a base ring with
              a coercion map into ``self.base_ring()``

            EXAMPLES:

            There exists a morphism from `NCSF` to `FSym`::

                sage: G = algebras.FSym(QQ).G()
                sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                sage: G(R[3,1,2,2,1])
                G[123|46|58|7|9] + G[123|46|58|79] + G[123|468|5|7|9]
                 + G[123|468|57|9] + G[123|468|579] + G[123|468|59|7]
                 + G[1236|478|5|9] + G[1236|478|59] + G[1236|48|5|7|9]
                 + G[1236|48|59|7] + G[12368|4|5|7|9] + G[12368|47|5|9]
                 + G[12368|47|59] + G[12368|479|5] + G[12368|49|5|7]
                 + G[1238|46|5|7|9] + G[1238|46|57|9] + G[1238|46|59|7]
                 + G[1238|469|5|7] + G[1238|469|57]
                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: G(S[2,1,2])
                G[12|35|4] + G[123|45] + G[12345] + G[1235|4]
                 + G[1245|3] + G[125|3|4] + G[125|34]
                sage: G(R(S[3,1,2,2])) == G(S[3,1,2,2])
                True

            This mapping is a Hopf algebra morphism::

                sage: all(G(R[a1] * R[a2]) == G(R[a1]) * G(R[a2])
                ....:     for a1 in Compositions(2)
                ....:     for a2 in Compositions(4))
                True

                sage: R2 = R.tensor_square()
                sage: phi = R2.module_morphism(
                ....:               lambda x: tensor([G(R[x[0]]), G(R[x[1]])]),
                ....:               codomain=G.tensor_square())
                sage: all(phi(R[p].coproduct()) == G(R[p]).coproduct()
                ....:     for p in Compositions(4))
                True

                sage: all(G(S[a1] * S[a2]) == G(S[a1]) * G(S[a2])
                ....:     for a1 in Compositions(2)
                ....:     for a2 in Compositions(4))
                True

                sage: S2 = S.tensor_square()
                sage: psi = S2.module_morphism(
                ....:               lambda x: tensor([G(S[x[0]]), G(S[x[1]])]),
                ....:               codomain=G.tensor_square())
                sage: all(psi(S[p].coproduct()) == G(S[p]).coproduct()
                ....:     for p in Compositions(4))
                True
            """
            if hasattr(R, "realization_of"):
                if not self.base_ring().has_coerce_map_from(R.base_ring()):
                    return False
                A = R.realization_of()
                # NSym to FSym
                from sage.combinat.ncsf_qsym.ncsf import NonCommutativeSymmetricFunctions
                if isinstance(A, NonCommutativeSymmetricFunctions):
                    ribbon = A.ribbon()
                    if R is ribbon:
                        ST = self._indices

                        def R_to_G_on_basis(alpha):
                            return self.sum_of_monomials(ST(t) for t in StandardTableaux(alpha.size())
                                                         if descent_composition(t) == alpha)
                        return ribbon.module_morphism(R_to_G_on_basis, codomain=self)
                    return self._coerce_map_via([ribbon], R)
            return super(FreeSymmetricFunctions.Fundamental, self)._coerce_map_from_(R)

        def dual_basis(self):
            r"""
            Return the dual basis to ``self``.

            EXAMPLES::

                sage: G = algebras.FSym(QQ).G()
                sage: G.dual_basis()
                Dual Hopf algebra of standard tableaux over the Rational Field
                 in the FundamentalDual basis
            """
            return self.realization_of().dual().F()

        @cached_method
        def product_on_basis(self, t1, t2):
            r"""
            Return the product of basis elements indexed by ``t1`` and ``t2``.

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: t1 = StandardTableau([[1,2], [3]])
                sage: t2 = StandardTableau([[1,2,3]])
                sage: G.product_on_basis(t1, t2)
                G[12456|3] + G[1256|3|4] + G[1256|34] + G[126|35|4]

                sage: t1 = StandardTableau([[1],[2]])
                sage: t2 = StandardTableau([[1,2]])
                sage: G.product_on_basis(t1, t2)
                G[134|2] + G[14|2|3]

                sage: t1 = StandardTableau([[1,2],[3]])
                sage: t2 = StandardTableau([[1],[2]])
                sage: G.product_on_basis(t1, t2)
                G[12|3|4|5] + G[12|34|5] + G[124|3|5] + G[124|35]
            """
            n = t1.size()
            m = n + t2.size()
            tableaux = []
            for t in StandardTableaux(m):
                if t.restrict(n) == t1 and standardize(t.anti_restrict(n).rectify()) == t2:
                    tableaux.append(t)
            return self.sum_of_monomials(tableaux)

        @cached_method
        def coproduct_on_basis(self, t):
            r"""
            Return the coproduct of the basis element indexed by ``t``.

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: t = StandardTableau([[1,2,5], [3,4]])
                sage: G.coproduct_on_basis(t)
                G[] # G[125|34] + G[1] # G[12|34] + G[1] # G[124|3]
                 + G[1|2] # G[13|2] + G[12] # G[12|3] + G[12] # G[123]
                 + G[12|34] # G[1] + G[123] # G[12] + G[125|34] # G[]
                 + G[13|2] # G[1|2] + G[13|2] # G[12] + G[134|2] # G[1]
            """
            # we use the duality to compute this
            n = t.size()
            L = []
            dual_basis = self.dual_basis()
            for i in range(n + 1):
                for t1 in StandardTableaux(i):
                    for t2 in StandardTableaux(n - i):
                        coeff = (dual_basis[t1] * dual_basis[t2])[t]
                        if coeff:
                            L.append(((t1, t2), coeff))
            TT = self.tensor_square()
            return TT.sum_of_terms(L)

        class Element(FSymBasis_abstract.Element):
            def to_fqsym(self):
                r"""
                Return the image of ``self`` under the natural inclusion
                map to `FQSym`.

                EXAMPLES::

                    sage: FSym = algebras.FSym(QQ)
                    sage: G = FSym.G()
                    sage: t = StandardTableau([[1,3],[2,4],[5]])
                    sage: G[t].to_fqsym()
                    G[2, 1, 5, 4, 3] + G[3, 1, 5, 4, 2] + G[3, 2, 5, 4, 1]
                     + G[4, 1, 5, 3, 2] + G[4, 2, 5, 3, 1]
                """
                from sage.combinat.fqsym import FreeQuasisymmetricFunctions
                R = self.parent().base_ring()
                G = FreeQuasisymmetricFunctions(R).G()
                return G(self)

            def to_symmetric_function(self):
                r"""
                Return the image of ``self`` under the natural projection
                map to `Sym`.

                The natural projection map `FSym \to Sym` sends each
                standard tableau `t` to the Schur function `s_\lambda`,
                where `\lambda` is the shape of `t`.
                This map is a surjective Hopf algebra homomorphism.

                EXAMPLES::

                    sage: FSym = algebras.FSym(QQ)
                    sage: G = FSym.G()
                    sage: t = StandardTableau([[1,3],[2,4],[5]])
                    sage: G[t].to_symmetric_function()
                    s[2, 2, 1]
                """
                s = SymmetricFunctions(self.parent().base_ring()).s()
                return s.sum_of_terms((t.shape(), coeff) for t, coeff in self)

    G = Fundamental


class FreeSymmetricFunctions_Dual(UniqueRepresentation, Parent):
    r"""
    The Hopf dual `FSym^*` of the free symmetric functions `FSym`.

    See :class:`FreeSymmetricFunctions` for the definition of the
    latter.

    Recall that the fundamental basis of `FSym` consists of the
    elements `\mathcal{G}_t` for `t` ranging over all standard
    tableaux. The dual basis of this is called the *dual
    fundamental basis* of `FSym^*`, and is denoted by
    `(\mathcal{G}_t^*)`.
    The Hopf dual `FSym^*` is isomorphic to the Hopf algebra
    `(\ZZ T, \ast', \delta')` from [PoiReu95]_; the
    isomorphism sends a basis element `\mathcal{G}_t^*` to `t`.

    EXAMPLES::

        sage: FSym = algebras.FSym(QQ)
        sage: TF = FSym.dual().F()
        sage: TF[1,2] * TF[[1],[2]]
        F[12|3|4] + F[123|4] + F[124|3] + F[13|2|4] + F[134|2] + F[14|2|3]
        sage: TF[[1,2],[3]].coproduct()
        F[] # F[12|3] + F[1] # F[1|2] + F[12] # F[1] + F[12|3] # F[]

    The Hopf algebra `FSym^*` is a Hopf quotient of `FQSym`;
    the canonical projection sends `F_w` (for a permutation `w`)
    to `\mathcal{G}_{Q(w)}^*`, where `Q(w)` is the Q-tableau of
    `w`. This projection is implemented as a coercion::

        sage: FQSym = algebras.FQSym(QQ)
        sage: F = FQSym.F()
        sage: TF(F[[1, 3, 2]])
        F[12|3]
        sage: TF(F[[5, 1, 4, 2, 3]])
        F[135|2|4]
    """
    def __init__(self, base_ring):
        r"""
        Initialize ``self``.

        TESTS::

            sage: FSymD = algebras.FSym(QQ).dual()
            sage: TestSuite(FSymD).run()
        """
        cat = HopfAlgebras(base_ring).Graded().Connected()
        Parent.__init__(self, base=base_ring, category=cat.WithRealizations())

    _shorthands = ['F']

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the Fundamental
        dual basis).

        EXAMPLES::

            sage: FSym = algebras.FSym(QQ).dual()
            sage: FSym.a_realization()
            Dual Hopf algebra of standard tableaux over the Rational Field
             in the FundamentalDual basis
        """
        return self.FundamentalDual()

    def dual(self):
        r"""
        Return the dual Hopf algebra of ``self``, which is `FSym`.

        EXAMPLES::

            sage: D = algebras.FSym(QQ).dual()
            sage: D.dual()
            Hopf algebra of standard tableaux over the Rational Field
        """
        return FreeSymmetricFunctions(self.base_ring())

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: algebras.FSym(QQ).dual()
            Dual Hopf algebra of standard tableaux over the Rational Field
        """
        return "Dual Hopf algebra of standard tableaux over the %s" % self.base_ring()

    class FundamentalDual(FSymBasis_abstract):
        r"""
        The dual to the Hopf algebra of tableaux,
        on the fundamental dual basis.

        EXAMPLES::

            sage: FSym = algebras.FSym(QQ)
            sage: TF = FSym.dual().F()
            sage: TF
            Dual Hopf algebra of standard tableaux over the Rational Field
             in the FundamentalDual basis

        Elements of the algebra look like::

            sage: TF.an_element()
            2*F[] + 2*F[1] + 3*F[12]

        TESTS::

            sage: FSym = algebras.FSym(QQ)
            sage: TF = FSym.dual().F()
            sage: TestSuite(TF).run()
        """
        _prefix = "F"

        def _coerce_map_from_(self, R):
            r"""
            Return ``True`` if there is a coercion from ``R`` into ``self``
            and ``False`` otherwise.

            The things that coerce into ``self`` are

            - elements of the algebra `FSym^*` over a base ring
              with a coercion map into ``self.base_ring()``
            - symmetric functions over a base ring with a coercion
              map into ``self.base_ring()``
            - elements of the algebra `FQSym` over a base ring with
              a coercion map into ``self.base_ring()``

            EXAMPLES:

            `FSym^*` is a quotient Hopf algebra of `FQSym`: the basis
            element `F_\sigma` indexed by a permutation `\sigma` is
            mapped to the tableau `Q(\sigma)`::

                sage: TF = algebras.FSym(QQ).dual().F()
                sage: SF = algebras.FQSym(QQ).F()
                sage: TF(SF([3,1,4,5,2]))
                F[134|25]
                sage: SG = algebras.FQSym(QQ).G()
                sage: TF(SG([3,1,4,5,2]))
                F[125|34]

            This mapping is a Hopf algebra morphism::

                sage: all(TF(SF[p1] * SF[p2]) == TF(SF[p1]) * TF(SF[p2])
                ....:     for p1 in Permutations(2)
                ....:     for p2 in Permutations(3))
                True

                sage: SSym2 = SF.tensor_square()
                sage: phi = SSym2.module_morphism(
                ....:               lambda x: tensor([TF(SF[x[0]]), TF(SF[x[1]])]),
                ....:               codomain=TF.tensor_square())
                sage: all(phi(SF[p].coproduct()) == TF(SF[p]).coproduct()
                ....:     for p in Permutations(4))
                True

            There is also an injective Hopf algebra morphism
            `Sym \to FSym^*` (adjoint to the projection `FSym \to Sym`
            implemented as
            :meth:`FreeSymmetricFunctions.Fundamental.Element.to_symmetric_function`)
            that sends each Schur function `s_\lambda` to the sum of
            all standard tableaux of shape `\lambda`::

                sage: Sym = SymmetricFunctions(QQ)
                sage: s = Sym.schur()
                sage: TF = algebras.FSym(QQ).dual().F()
                sage: TF(s[2,1])
                F[12|3] + F[13|2]
                sage: TF(s[2,2,1])
                F[12|34|5] + F[12|35|4] + F[13|24|5] + F[13|25|4] + F[14|25|3]
                sage: h = Sym.h()
                sage: TF(h[2,1])
                F[12|3] + F[123] + F[13|2]

            This mapping is a Hopf algebra morphism::

                sage: all(TF(s[p1] * s[p2]) == TF(s[p1]) * TF(s[p2])
                ....:     for p1 in Partitions(2)
                ....:     for p2 in Partitions(3))
                True

                sage: s2 = s.tensor_square()
                sage: phi = s2.module_morphism(
                ....:               lambda x: tensor([TF(s[x[0]]), TF(s[x[1]])]),
                ....:               codomain=TF.tensor_square())
                sage: all(phi(s[p].coproduct()) == TF(s[p]).coproduct()
                ....:     for p in Partitions(4))
                True
            """
            if hasattr(R, "realization_of"):
                if not self.base_ring().has_coerce_map_from(R.base_ring()):
                    return False
                A = R.realization_of()
                # FQSym to FSym^*
                from sage.combinat.fqsym import FreeQuasisymmetricFunctions
                if isinstance(A, FreeQuasisymmetricFunctions):
                    F = A.F()
                    if R is F:
                        def F_to_SF_on_basis(sigma):
                            return self.monomial(sigma.right_tableau())
                        return F.module_morphism(F_to_SF_on_basis, codomain=self)
                    return self._coerce_map_via([F], R)

                # Sym to FSym^*
                if isinstance(A, SymmetricFunctions):
                    s = A.s()
                    if R is s:
                        def s_to_F_on_basis(mu):
                            return self.sum_of_monomials(StandardTableaux(mu))
                        return s.module_morphism(s_to_F_on_basis, codomain=self)
                    return self._coerce_map_via([s], R)
            return super(FreeSymmetricFunctions_Dual.FundamentalDual, self)._coerce_map_from_(R)

        def dual_basis(self):
            r"""
            Return the dual basis to ``self``.

            EXAMPLES::

                sage: F = algebras.FSym(QQ).dual().F()
                sage: F.dual_basis()
                Hopf algebra of standard tableaux over the Rational Field
                 in the Fundamental basis
            """
            return self.realization_of().dual().G()

        @cached_method
        def product_on_basis(self, t1, t2):
            r"""
            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: TF = FSym.dual().F()
                sage: t1 = StandardTableau([[1,2]])
                sage: TF.product_on_basis(t1, t1)
                F[12|34] + F[123|4] + F[1234] + F[124|3] + F[13|24] + F[134|2]
                sage: t0 = StandardTableau([])
                sage: TF.product_on_basis(t1, t0) == TF[t1] == TF.product_on_basis(t0, t1)
                True
            """
            z = []
            n = t1.size()
            m = t2.size()
            npmp1 = n + m + 1
            ST = self._indices
            from itertools import combinations
            for I in combinations(range(1, npmp1), n):
                J = [j for j in range(1, npmp1) if (j not in I)]
                tt1 = [[I[x - 1] for x in row] for row in t1]
                tt2 = [tuple([J[x - 1] for x in row]) for row in t2]
                z.append(ST(Tableau(tt1).slide_multiply(tt2)))
            return self.sum_of_monomials(z)

        @cached_method
        def coproduct_on_basis(self, t):
            r"""
            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: TF = FSym.dual().F()
                sage: t = StandardTableau([[1,2,5], [3,4]])
                sage: TF.coproduct_on_basis(t)
                F[] # F[125|34] + F[1] # F[134|2] + F[12] # F[123]
                 + F[12|3] # F[12] + F[12|34] # F[1] + F[125|34] # F[]
            """
            terms = [(t.restrict(i), standardize(t.anti_restrict(i).rectify()))
                     for i in range(t.size() + 1)]
            return self.tensor_square().sum_of_monomials(terms)

        class Element(FSymBasis_abstract.Element):
            def to_quasisymmetric_function(self):
                r"""
                Return the image of ``self`` under the canonical projection
                `FSym^* \to QSym` to the ring of quasi-symmetric functions.

                This projection is the adjoint of the canonical injection
                `NSym \to FSym` (see
                :meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_fsym`).
                It sends each tableau `t` to the fundamental quasi-symmetric
                function `F_\alpha`, where `\alpha` is the descent composition
                of `t`.

                EXAMPLES::

                    sage: F = algebras.FSym(QQ).dual().F()
                    sage: F[[1,3,5],[2,4]].to_quasisymmetric_function()
                    F[1, 2, 2]
                """
                from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
                QF = QuasiSymmetricFunctions(self.base_ring()).Fundamental()
                return QF.sum_of_terms((descent_composition(t), coeff)
                                       for t, coeff in self)

    F = FundamentalDual


# some utility functions for tableaux

def standardize(t):
    r"""
    Return the standard tableau corresponding to a given
    semistandard tableau ``t`` with no repeated entries.

    .. NOTE::

        This is an optimized version of :meth:`Tableau.standardization`
        for computations in `FSym` by using the assumption of no
        repeated entries in ``t``.

    EXAMPLES::

        sage: from sage.combinat.chas.fsym import standardize
        sage: t = Tableau([[1,3,5,7],[2,4,8],[9]])
        sage: standardize(t)
        [[1, 3, 5, 6], [2, 4, 7], [8]]
        sage: t = Tableau([[3,8,9,15],[5,10,12],[133]])
        sage: standardize(t)
        [[1, 3, 4, 7], [2, 5, 6], [8]]

    TESTS:

    This returns an equal tableau if already standard::

        sage: t = Tableau([[1,3,4,5],[2,6,7],[8]])
        sage: standardize(t)
        [[1, 3, 4, 5], [2, 6, 7], [8]]
        sage: standardize(t) == t
        True
    """
    A = sorted(sum(t, ()))
    std = {j: i + 1 for i, j in enumerate(A)}
    ST = StandardTableaux()
    return ST([[std[i] for i in row] for row in t])


def ascent_set(t):
    """
    Return the ascent set of a standard tableau ``t``
    (encoded as a sorted list).

    The *ascent set* of a standard tableau `t` is defined as
    the set of all entries `i` of `t` such that the number `i+1`
    either appears to the right of `i` or appears in a row above
    `i` or does not appear in `t` at all.

    EXAMPLES::

        sage: from sage.combinat.chas.fsym import ascent_set
        sage: t = StandardTableau([[1,3,4,7],[2,5,6],[8]])
        sage: ascent_set(t)
        [2, 3, 5, 6, 8]
        sage: ascent_set(StandardTableau([]))
        []
        sage: ascent_set(StandardTableau([[1, 2, 3]]))
        [1, 2, 3]
        sage: ascent_set(StandardTableau([[1, 2, 4], [3]]))
        [1, 3, 4]
        sage: ascent_set([[1, 3, 5], [2, 4]])
        [2, 4, 5]
    """
    row_locations = {}
    for (i, row) in enumerate(t):
        for entry in row:
            row_locations[entry] = i
    n = len(row_locations)
    if not n:
        return []
    ascents = [n]
    for i in range(1, n):
        # ascent means i+1 appears to the right or above
        x = row_locations[i]
        u = row_locations[i + 1]
        if u <= x:
            ascents.append(i)
    return sorted(ascents)


def descent_set(t):
    """
    Return the descent set of a standard tableau ``t``
    (encoded as a sorted list).

    The *descent set* of a standard tableau `t` is defined as
    the set of all entries `i` of `t` such that the number `i+1`
    appears in a row below `i` in `t`.

    EXAMPLES::

        sage: from sage.combinat.chas.fsym import descent_set
        sage: t = StandardTableau([[1,3,4,7],[2,5,6],[8]])
        sage: descent_set(t)
        [1, 4, 7]
        sage: descent_set(StandardTableau([]))
        []
        sage: descent_set(StandardTableau([[1, 2, 3]]))
        []
        sage: descent_set(StandardTableau([[1, 2, 4], [3]]))
        [2]
        sage: descent_set([[1, 3, 5], [2, 4]])
        [1, 3]
    """
    ascents = set(ascent_set(t))
    n = sum(len(row) for row in t)
    return [i for i in range(1, n) if i not in ascents]


def descent_composition(t):
    """
    Return the descent composition of a standard tableau ``t``.

    This is the composition of the size of `t` whose partial
    sums are the elements of the descent set of ``t`` (see
    :meth:`descent_set`).

    EXAMPLES::

        sage: from sage.combinat.chas.fsym import descent_composition
        sage: t = StandardTableau([[1,3,4,7],[2,5,6],[8]])
        sage: descent_composition(t)
        [1, 3, 3, 1]
        sage: descent_composition([[1, 3, 5], [2, 4]])
        [1, 2, 2]
    """
    n = sum(len(row) for row in t)
    return Composition(from_subset=(descent_set(t), n))
