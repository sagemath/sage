# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Free Quasi-symmetric functions

AUTHORS:

- Frédéric Chapoton, Darij Grinberg (2017)
"""

# ****************************************************************************
#       Copyright (C) 2010-2015 Frédéric Chapoton <chapoton@unistra.fr>,
#       Copyright (C) 2017      Darij Grinberg <dgrinber at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations, Permutation
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
from sage.combinat.words.word import Word

lazy_import('sage.groups.perm_gps.permgroup_element', 'PermutationGroupElement')


class FQSymBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for bases of FQSym.

    This must define two attributes:

    - ``_prefix`` -- the basis prefix
    - ``_basis_name`` -- the name of the basis and must match one
      of the names that the basis can be constructed from FQSym
    """

    def __init__(self, alg):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: TestSuite(algebras.FQSym(QQ).F()).run()  # long time
        """
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         Permutations(),
                                         category=FQSymBases(alg),
                                         bracket='', prefix=self._prefix)

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - free quasi-symmetric functions over a base with
          a coercion map into ``self.base_ring()``
        - free symmetric functions over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: F = algebras.FQSym(GF(7)).F(); F
            Free Quasi-symmetric functions over Finite Field of size 7 in the F basis

        Elements of the free quasi-symmetric functions canonically coerce in::

            sage: x, y, z = F([1]), F([2,1]), F([1,3,2])
            sage: F.coerce(x+y) == x+y
            True

        The free quasi-symmetric functions over `\ZZ` coerces in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: G = algebras.FQSym(ZZ).F()
            sage: Gx, Gy = G([1]), G([2,1])
            sage: z = F.coerce(Gx+Gy); z
            F[1] + F[2, 1]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so free
        quasi-symmetric functions over `\GF{7}` does not coerce
        to the same algebra over `\ZZ`::

            sage: G.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free Quasi-symmetric functions
             over Finite Field of size 7 in the F basis to
             Free Quasi-symmetric functions over Integer Ring in the F basis

        Check that `FSym` bases coerce in::

            sage: FSym = algebras.FSym(ZZ)
            sage: TG = FSym.G()
            sage: t = StandardTableau([[1,3],[2,4],[5]])
            sage: F(TG[t])
            F[2, 1, 5, 4, 3] + F[2, 5, 1, 4, 3] + F[2, 5, 4, 1, 3]
             + F[5, 2, 1, 4, 3] + F[5, 2, 4, 1, 3]
            sage: algebras.FQSym(QQ)(TG[t])
            F[2, 1, 5, 4, 3] + F[2, 5, 1, 4, 3] + F[2, 5, 4, 1, 3]
             + F[5, 2, 1, 4, 3] + F[5, 2, 4, 1, 3]
            sage: G7 = algebras.FQSym(GF(7)).G()
            sage: G7(TG[[1,2],[3,4]])
            G[2, 4, 1, 3] + G[3, 4, 1, 2]

        TESTS::

            sage: F = algebras.FQSym(ZZ).F()
            sage: G = algebras.FQSym(QQ).F()
            sage: F.has_coerce_map_from(G)
            False
            sage: G.has_coerce_map_from(F)
            True
            sage: F.has_coerce_map_from(QQ)
            False
            sage: G.has_coerce_map_from(QQ)
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # free quasi-symmetric functions in the same variables
        # over any base that coerces in:
        if isinstance(R, FQSymBasis_abstract):
            if R.realization_of() == self.realization_of():
                return True
            if not self.base_ring().has_coerce_map_from(R.base_ring()):
                return False
            if self._basis_name == R._basis_name:  # The same basis
                def coerce_base_ring(self, x):
                    return self._from_dict(x.monomial_coefficients())
                return coerce_base_ring
            # Otherwise lift that basis up and then coerce over
            target = getattr(self.realization_of(), R._basis_name)()
            return self._coerce_map_via([target], R)

        # FSym coerces in:
        from sage.combinat.chas.fsym import FreeSymmetricFunctions
        if isinstance(R, FreeSymmetricFunctions.Fundamental):
            if not self.base_ring().has_coerce_map_from(R.base_ring()):
                return False
            G = self.realization_of().G()
            P = G._indices

            def G_to_G_on_basis(t):
                return G.sum_of_monomials(P(sigma) for sigma in Permutations(t.size())
                                          if sigma.right_tableau() == t)
            phi = R.module_morphism(G_to_G_on_basis, codomain=G)
            if self is G:
                return phi
            else:
                return self.coerce_map_from(G) * phi

        return super()._coerce_map_from_(R)

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = algebras.FQSym(QQ)
            sage: F = A.F()
            sage: F.an_element()
            F[1] + 2*F[1, 2] + 2*F[2, 1]
            sage: G = A.G()
            sage: G.an_element()
            G[1] + 2*G[1, 2] + 2*G[2, 1]
            sage: M = A.M()
            sage: M.an_element()
            M[1] + 2*M[1, 2] + 4*M[2, 1]
        """
        o = self.monomial(Permutation([1]))
        return o + 2 * o * o


class FreeQuasisymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The free quasi-symmetric functions.

    The Hopf algebra `FQSym` of free quasi-symmetric functions
    over a commutative ring `R` is the free `R`-module with basis
    indexed by all permutations (i.e., the indexing set is
    the disjoint union of all symmetric groups).
    Its product is determined by the shifted shuffles of two
    permutations, whereas its coproduct is given by splitting
    a permutation (regarded as a word) into two (at every
    possible point) and standardizing the two pieces.
    This Hopf algebra was introduced in [MR]_.
    See [GriRei18]_ (Chapter 8) for a treatment using modern
    notations.

    In more detail:
    For each `n \geq 0`, consider the symmetric group `S_n`.
    Let `S` be the disjoint union of the `S_n` over all
    `n \geq 0`.
    Then, `FQSym` is the free `R`-module with basis
    `(F_w)_{w \in S}`.
    This `R`-module is graded, with the `n`-th graded
    component being spanned by all `F_w` for `w \in S_n`.
    A multiplication is defined on `FQSym` as follows:
    For any two permutations `u \in S_k` and `v \in S_l`,
    we set

    .. MATH::

        F_u F_v = \sum F_w ,

    where the sum is over all shuffles of `u` with `v[k]`.
    Here, the permutations `u` and `v` are regarded as words
    (by writing them in one-line notation), and `v[k]` means
    the word obtained from `v` by increasing each letter by
    `k` (for example, `(1,4,2,3)[5] = (6,9,7,8)`); and the
    shuffles `w` are translated back into permutations.
    This defines an associative multiplication on `FQSym`;
    its unity is `F_e`, where `e` is the identity
    permutation in `S_0`.

    In Section 1.3 of [AguSot05]_, Aguiar and Sottile construct a
    different basis of `FQSym`. Their basis, called the
    *monomial basis* and denoted by `(\mathcal{M}_u)`,
    is also indexed by permutations. It is connected to the
    above F-basis by the relation

    .. MATH::

        F_u = \sum_v \mathcal{M}_v ,

    where the sum ranges over all permutations `v` such that each
    inversion of `u` is an inversion of `v`. (An *inversion* of a
    permutation `w` means a pair `(i, j)` of positions satisfying
    `i < j` and `w(i) > w(j)`.) The above relation yields a
    unitriangular change-of-basis matrix, and thus can be used to
    compute the `\mathcal{M}_u` by Mobius inversion.

    Another classical basis of `FQSym` is `(G_w)_{w \in S}`,
    where `G_w = F_{w^{-1}}`.
    This is just a relabeling of the basis `(F_w)_{w \in S}`,
    but is a more natural choice from some viewpoints.

    The algebra `FQSym` is often identified with ("realized as") a
    subring of the ring of all bounded-degree noncommutative power
    series in countably many indeterminates (i.e., elements in
    `R \langle \langle x_1, x_2, x_3, \ldots \rangle \rangle` of bounded
    degree). Namely, consider words over the alphabet `\{1, 2, 3, \ldots\}`;
    every noncommutative power series is an infinite `R`-linear
    combination of these words.
    Consider the `R`-linear map that sends each `G_u` to the sum of
    all words whose standardization (also known as "standard
    permutation"; see
    :meth:`~sage.combinat.words.finite_word.FiniteWord_class.standard_permutation`)
    is `u`. This map is an injective `R`-algebra homomorphism, and
    thus embeds `FQSym` into the latter ring.

    As an associative algebra, `FQSym` has the richer structure
    of a dendriform algebra. This means that the associative
    product ``*`` is decomposed as a sum of two binary operations

    .. MATH::

        x y = x \succ y + x \prec y

    that satisfy the axioms:

    .. MATH::

        (x \succ y) \prec z = x \succ (y \prec z),

    .. MATH::

        (x \prec y) \prec z = x \prec (y z),

    .. MATH::

        (x y) \succ z = x \succ (y \succ z).

    These two binary operations are defined similarly to the
    (associative) product above: We set

    .. MATH::

        F_u \prec F_v = \sum F_w ,

    where the sum is now over all shuffles of `u` with `v[k]`
    whose first letter is taken from `u` (rather than from
    `v[k]`). Similarly,

    .. MATH::

        F_u \succ F_v = \sum F_w ,

    where the sum is over all remaining shuffles of `u` with
    `v[k]`.

    .. TODO::

        Decide what `1 \prec 1` and `1 \succ 1` are.

    .. NOTE::

        The usual binary operator ``*`` is used for the
        associative product.

    EXAMPLES::

        sage: F = algebras.FQSym(ZZ).F()
        sage: x,y,z = F([1]), F([1,2]), F([1,3,2])
        sage: (x * y) * z
        F[1, 2, 3, 4, 6, 5] + ...

    The product of `FQSym` is associative::

        sage: x * (y * z) == (x * y) * z
        True

    The associative product decomposes into two parts::

        sage: x * y == F.prec(x, y) + F.succ(x, y)
        True

    The axioms of a dendriform algebra hold::

        sage: F.prec(F.succ(x, y), z) == F.succ(x, F.prec(y, z))
        True
        sage: F.prec(F.prec(x, y), z) == F.prec(x, y * z)
        True
        sage: F.succ(x * y, z) == F.succ(x, F.succ(y, z))
        True

    `FQSym` is also known as the Malvenuto-Reutenauer algebra::

        sage: algebras.MalvenutoReutenauer(ZZ)
        Free Quasi-symmetric functions over Integer Ring

    REFERENCES:

    - [MR]_
    - [LR1998]_
    - [GriRei18]_
    """

    def __init__(self, R):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.FQSym(QQ); A
            Free Quasi-symmetric functions over Rational Field
            sage: TestSuite(A).run()  # long time (3s)

            sage: F = algebras.FQSym(QQ)
            sage: TestSuite(F).run() # long time (3s)
        """
        category = HopfAlgebras(R).Graded().Connected()
        Parent.__init__(self, base=R, category=category.WithRealizations())

        # Bases
        F = self.F()
        G = self.G()

        F.module_morphism(G._F_to_G_on_basis,
                          codomain=G, category=category).register_as_coercion()
        G.module_morphism(G._G_to_F_on_basis,
                          codomain=F, category=category).register_as_coercion()

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: algebras.FQSym(QQ)  # indirect doctest
            Free Quasi-symmetric functions over Rational Field
        """
        s = "Free Quasi-symmetric functions over {}"
        return s.format(self.base_ring())

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the F-basis).

        EXAMPLES::

            sage: FQSym = algebras.FQSym(QQ)
            sage: FQSym.a_realization()
            Free Quasi-symmetric functions over Rational Field in the F basis
        """
        return self.F()

    _shorthands = tuple(['F', 'G', 'M'])

    class F(FQSymBasis_abstract):
        r"""
        The F-basis of `FQSym`.

        This is the basis `(F_w)`, with `w` ranging over all
        permutations. See the documentation of
        :class:`FreeQuasisymmetricFunctions` for details.

        EXAMPLES::

            sage: FQSym = algebras.FQSym(QQ)
            sage: FQSym.F()
            Free Quasi-symmetric functions over Rational Field in the F basis
        """
        _prefix = "F"
        _basis_name = "F"

        def _element_constructor_(self, x):
            r"""
            Convert ``x`` into ``self``.

            EXAMPLES::

                sage: R = algebras.FQSym(QQ).F()
                sage: x, y, z = R([1]), R([2,1]), R([3,2,1])
                sage: R(x)
                F[1]
                sage: R(x+4*y)
                F[1] + 4*F[2, 1]
                sage: R(1)
                F[]

                sage: D = algebras.FQSym(ZZ).F()
                sage: X, Y, Z = D([1]), D([2,1]), D([3,2,1])
                sage: R(X-Y).parent()
                Free Quasi-symmetric functions over Rational Field in the F basis

                sage: R([1, 3, 2])
                F[1, 3, 2]
                sage: R(Permutation([1, 3, 2]))
                F[1, 3, 2]
                sage: R(SymmetricGroup(4)(Permutation([1,3,4,2])))
                F[1, 3, 4, 2]
            """
            if isinstance(x, (list, tuple, PermutationGroupElement)):
                x = Permutation(x)
            try:
                P = x.parent()
                if isinstance(P, FreeQuasisymmetricFunctions.F):
                    if P is self:
                        return x
                    return self.element_class(self, x.monomial_coefficients())
            except AttributeError:
                pass
            return CombinatorialFreeModule._element_constructor_(self, x)

        def __getitem__(self, r):
            r"""
            The default implementation of ``__getitem__`` interprets
            the input as a tuple, which in case of permutations
            is interpreted as cycle notation, even though the input
            looks like a one-line notation.
            We override this method to amend this.

            EXAMPLES::

                sage: F = algebras.FQSym(QQ).F()
                sage: F[3, 2, 1]
                F[3, 2, 1]
                sage: F[1]
                F[1]
            """
            if isinstance(r, tuple):
                r = list(r)
            elif r == 1:
                r = [1]
            return super().__getitem__(r)

        def degree_on_basis(self, t):
            """
            Return the degree of a permutation in
            the algebra of free quasi-symmetric functions.

            This is the size of the permutation (i.e., the `n`
            for which the permutation belongs to `S_n`).

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: u = Permutation([2,1])
                sage: A.degree_on_basis(u)
                2
            """
            return len(t)

        def product_on_basis(self, x, y):
            r"""
            Return the `*` associative product of two permutations.

            This is the shifted shuffle of `x` and `y`.

            .. SEEALSO::

                :meth:`succ_product_on_basis`, :meth:`prec_product_on_basis`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = Permutation([1])
                sage: A.product_on_basis(x, x)
                F[1, 2] + F[2, 1]
            """
            return self.sum_of_monomials(u for u in x.shifted_shuffle(y))

        def succ_product_on_basis(self, x, y):
            r"""
            Return the `\succ` product of two permutations.

            This is the shifted shuffle of `x` and `y` with the additional
            condition that the first letter of the result comes from `y`.

            The usual symbol for this operation is `\succ`.

            .. SEEALSO::

                - :meth:`product_on_basis`, :meth:`prec_product_on_basis`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = Permutation([1,2])
                sage: A.succ_product_on_basis(x, x)
                F[3, 1, 2, 4] + F[3, 1, 4, 2] + F[3, 4, 1, 2]
                sage: y = Permutation([])
                sage: A.succ_product_on_basis(x, y) == 0
                True
                sage: A.succ_product_on_basis(y, x) == A(x)
                True

            TESTS::

                sage: u = A.one().support()[0] # this is F[]
                sage: A.succ_product_on_basis(x, u)
                0
                sage: A.succ_product_on_basis(u, x)
                F[1, 2]
                sage: A.succ_product_on_basis(u, u)
                Traceback (most recent call last):
                ...
                ValueError: products | < | and | > | are not defined
            """
            if not y:
                if not x:
                    raise ValueError("products | < | and | > | are not defined")
                else:
                    return self.zero()
            basis = self.basis()
            if not x:
                return basis[y]
            K = basis.keys()
            n = len(x)
            shy = Word([a + n for a in y])
            shy0 = shy[0]
            return self.sum_of_monomials(K([shy0] + list(u))
                                         for u in Word(x).shuffle(Word(shy[1:])))

        def prec_product_on_basis(self, x, y):
            r"""
            Return the `\prec` product of two permutations.

            This is the shifted shuffle of `x` and `y` with the additional
            condition that the first letter of the result comes from `x`.

            The usual symbol for this operation is `\prec`.

            .. SEEALSO::

                :meth:`product_on_basis`, :meth:`succ_product_on_basis`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = Permutation([1,2])
                sage: A.prec_product_on_basis(x, x)
                F[1, 2, 3, 4] + F[1, 3, 2, 4] + F[1, 3, 4, 2]
                sage: y = Permutation([])
                sage: A.prec_product_on_basis(x, y) == A(x)
                True
                sage: A.prec_product_on_basis(y, x) == 0
                True

            TESTS::

                sage: u = A.one().support()[0] # this is F[]
                sage: A.prec_product_on_basis(x, u)
                F[1, 2]
                sage: A.prec_product_on_basis(u, x)
                0
                sage: A.prec_product_on_basis(u, u)
                Traceback (most recent call last):
                ...
                ValueError: products | < | and | > | are not defined
            """
            if not x and not y:
                raise ValueError("products | < | and | > | are not defined")
            if not x:
                return self.zero()
            basis = self.basis()
            if not y:
                return basis[x]
            K = basis.keys()
            n = len(x)
            shy = Word([a + n for a in y])
            x0 = x[0]
            return self.sum_of_monomials(K([x0] + list(u))
                                         for u in Word(x[1:]).shuffle(shy))

        def coproduct_on_basis(self, x):
            r"""
            Return the coproduct of `F_{\sigma}` for `\sigma` a permutation
            (here, `\sigma` is ``x``).

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = A([1])
                sage: ascii_art(A.coproduct(A.one()))  # indirect doctest
                1 # 1

                sage: ascii_art(A.coproduct(x))  # indirect doctest
                1 # F    + F    # 1
                     [1]    [1]

                sage: A = algebras.FQSym(QQ).F()
                sage: x, y, z = A([1]), A([2,1]), A([3,2,1])
                sage: A.coproduct(z)
                F[] # F[3, 2, 1] + F[1] # F[2, 1] + F[2, 1] # F[1]
                + F[3, 2, 1] # F[]
            """
            if not x:
                return self.one().tensor(self.one())
            return sum(self(Word(x[:i]).standard_permutation()).tensor(
                self(Word(x[i:]).standard_permutation()))
                for i in range(len(x) + 1))

        class Element(FQSymBasis_abstract.Element):
            def to_symmetric_group_algebra(self, n=None):
                """
                Return the element of a symmetric group algebra
                corresponding to the element ``self`` of `FQSym`.

                INPUT:

                - ``n`` -- integer (default: the maximal degree of ``self``);
                  the rank of the target symmetric group algebra

                EXAMPLES::

                    sage: A = algebras.FQSym(QQ).F()
                    sage: x = A([1,3,2,4]) + 5/2 * A([1,2,4,3])
                    sage: x.to_symmetric_group_algebra()
                    5/2*[1, 2, 4, 3] + [1, 3, 2, 4]
                    sage: x.to_symmetric_group_algebra(n=7)
                    5/2*[1, 2, 4, 3, 5, 6, 7] + [1, 3, 2, 4, 5, 6, 7]
                    sage: a = A.zero().to_symmetric_group_algebra(); a
                    0
                    sage: parent(a)
                    Symmetric group algebra of order 0 over Rational Field

                    sage: y = A([1,3,2,4]) + 5/2 * A([2,1])
                    sage: y.to_symmetric_group_algebra()
                    [1, 3, 2, 4] + 5/2*[2, 1, 3, 4]
                    sage: y.to_symmetric_group_algebra(6)
                    [1, 3, 2, 4, 5, 6] + 5/2*[2, 1, 3, 4, 5, 6]
                """
                if not self:
                    if n is None:
                        n = 0
                    return SymmetricGroupAlgebra(self.base_ring(), n).zero()
                m = self.maximal_degree()
                if n is None:
                    n = m
                elif n < m:
                    raise ValueError("n must be at least the maximal degree")

                SGA = SymmetricGroupAlgebra(self.base_ring(), n)
                return SGA._from_dict({Permutations(n)(key): c for (key, c) in self})

    class G(FQSymBasis_abstract):
        r"""
        The G-basis of `FQSym`.

        This is the basis `(G_w)`, with `w` ranging over all
        permutations. See the documentation of
        :class:`FreeQuasisymmetricFunctions` for details.

        EXAMPLES::

            sage: FQSym = algebras.FQSym(QQ)
            sage: G = FQSym.G(); G
            Free Quasi-symmetric functions over Rational Field in the G basis

            sage: G([3, 1, 2]).coproduct()
            G[] # G[3, 1, 2] + G[1] # G[2, 1] + G[1, 2] # G[1]
             + G[3, 1, 2] # G[]

            sage: G([3, 1, 2]) * G([2, 1])
            G[3, 1, 2, 5, 4] + G[4, 1, 2, 5, 3] + G[4, 1, 3, 5, 2]
             + G[4, 2, 3, 5, 1] + G[5, 1, 2, 4, 3] + G[5, 1, 3, 4, 2]
             + G[5, 1, 4, 3, 2] + G[5, 2, 3, 4, 1] + G[5, 2, 4, 3, 1]
             + G[5, 3, 4, 2, 1]
        """
        _prefix = "G"
        _basis_name = "G"

        def _element_constructor_(self, x):
            r"""
            Convert ``x`` into ``self``.

            EXAMPLES::

                sage: R = algebras.FQSym(QQ).G()
                sage: x, y, z = R([1]), R([2,1]), R([3,2,1])
                sage: R(x)
                G[1]
                sage: R(x+4*y)
                G[1] + 4*G[2, 1]
                sage: R(1)
                G[]

                sage: D = algebras.FQSym(ZZ).G()
                sage: X, Y, Z = D([1]), D([2,1]), D([3,2,1])
                sage: R(X-Y).parent()
                Free Quasi-symmetric functions over Rational Field in the G basis

                sage: R([1, 3, 2])
                G[1, 3, 2]
                sage: R(Permutation([1, 3, 2]))
                G[1, 3, 2]
                sage: R(SymmetricGroup(4)(Permutation([1,3,4,2])))
                G[1, 3, 4, 2]

                sage: RF = algebras.FQSym(QQ).F()
                sage: R(RF([2, 3, 4, 1]))
                G[4, 1, 2, 3]
                sage: R(RF([3, 2, 4, 1]))
                G[4, 2, 1, 3]
                sage: DF = algebras.FQSym(ZZ).F()
                sage: D(DF([2, 3, 4, 1]))
                G[4, 1, 2, 3]
                sage: R(DF([2, 3, 4, 1]))
                G[4, 1, 2, 3]
                sage: RF(R[2, 3, 4, 1])
                F[4, 1, 2, 3]
            """
            if isinstance(x, (list, tuple, PermutationGroupElement)):
                x = Permutation(x)
            try:
                P = x.parent()
                if isinstance(P, FreeQuasisymmetricFunctions.G):
                    if P is self:
                        return x
                    return self.element_class(self, x.monomial_coefficients())
            except AttributeError:
                pass
            return CombinatorialFreeModule._element_constructor_(self, x)

        def __getitem__(self, r):
            r"""
            The default implementation of ``__getitem__`` interprets
            the input as a tuple, which in case of permutations
            is interpreted as cycle notation, even though the input
            looks like a one-line notation.
            We override this method to amend this.

            EXAMPLES::

                sage: G = algebras.FQSym(QQ).G()
                sage: G[3, 2, 1]
                G[3, 2, 1]
                sage: G[1]
                G[1]
            """
            if isinstance(r, tuple):
                r = list(r)
            elif r == 1:
                r = [1]
            return super().__getitem__(r)

        def _G_to_F_on_basis(self, w):
            r"""
            Return `G_w` in terms of the F basis.

            INPUT:

            - ``w`` -- a permutation

            OUTPUT: an element of the F basis

            TESTS::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: F = FQSym.F()
                sage: G = FQSym.G()
                sage: F(G[3, 2, 1] - 4 * G[4, 2, 1, 3])
                F[3, 2, 1] - 4*F[3, 2, 4, 1]
                sage: all(F(G._G_to_F_on_basis(w)) == G[w] for i in range(5)
                ....:     for w in Permutations(i))
                True
                sage: G[3, 2, 1] == F[3, 2, 1]
                True
                sage: G[4, 2, 1, 3] == F[3, 2, 4, 1]
                True
                sage: G[4, 2, 1, 3] == F[4, 2, 1, 3]
                False
            """
            F = self.realization_of().F()
            return F.basis()[w.inverse()]

        def _F_to_G_on_basis(self, w):
            r"""
            Return `F_w` in terms of the G basis.

            INPUT:

            - ``w`` -- a permutation

            OUTPUT: an element of the G basis

            TESTS::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: F = FQSym.F()
                sage: G = FQSym.G()
                sage: G(F[3, 2, 1] - 4 * F[4, 2, 1, 3])
                G[3, 2, 1] - 4*G[3, 2, 4, 1]
                sage: all(G(G._F_to_G_on_basis(w)) == F[w] for i in range(5)
                ....:     for w in Permutations(i))
                True
                sage: F[3, 2, 1] == G[3, 2, 1]
                True
                sage: F[4, 2, 1, 3] == G[3, 2, 4, 1]
                True
                sage: F[4, 2, 1, 3] == G[4, 2, 1, 3]
                False
            """
            return self.basis()[w.inverse()]

        def degree_on_basis(self, t):
            """
            Return the degree of a permutation in
            the algebra of free quasi-symmetric functions.

            This is the size of the permutation (i.e., the `n`
            for which the permutation belongs to `S_n`).

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).G()
                sage: u = Permutation([2,1])
                sage: A.degree_on_basis(u)
                2
            """
            return len(t)

    class M(FQSymBasis_abstract):
        r"""
        The M-basis of `FQSym`.

        This is the Monomial basis `(\mathcal{M}_w)`, with `w` ranging
        over all permutations. See the documentation of :class:`FQSym`
        for details.

        EXAMPLES::

            sage: FQSym = algebras.FQSym(QQ)
            sage: M = FQSym.M(); M
            Free Quasi-symmetric functions over Rational Field in the Monomial basis

            sage: M([3, 1, 2]).coproduct()
            M[] # M[3, 1, 2] + M[1] # M[1, 2] + M[3, 1, 2] # M[]
            sage: M([3, 2, 1]).coproduct()
            M[] # M[3, 2, 1] + M[1] # M[2, 1] + M[2, 1] # M[1]
             + M[3, 2, 1] # M[]

            sage: M([1, 2]) * M([1])
            M[1, 2, 3] + 2*M[1, 3, 2] + M[2, 3, 1] + M[3, 1, 2]
        """
        _prefix = "M"
        _basis_name = "Monomial"

        def __init__(self, alg):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: M = algebras.FQSym(QQ).M()
                sage: TestSuite(M).run(elements=M.some_elements()[:-1])  # long time
            """
            FQSymBasis_abstract.__init__(self, alg)

            F = self.realization_of().F()
            phi = F.module_morphism(self._F_to_M_on_basis, codomain=self,
                                    unitriangular='lower')
            phi.register_as_coercion()
            phi_i = self.module_morphism(self._M_to_F_on_basis, codomain=F,
                                         unitriangular='lower')
            phi_i.register_as_coercion()

        def _element_constructor_(self, x):
            r"""
            Convert ``x`` into ``self``.

            EXAMPLES::

                sage: R = algebras.FQSym(QQ).M()
                sage: x, y, z = R([1]), R([2,1]), R([3,2,1])
                sage: R(x)
                M[1]
                sage: R(x+4*y)
                M[1] + 4*M[2, 1]
                sage: R(1)
                M[]

                sage: D = algebras.FQSym(ZZ).M()
                sage: X, Y, Z = D([1]), D([2,1]), D([3,2,1])
                sage: R(X-Y).parent()
                Free Quasi-symmetric functions over Rational Field in the Monomial basis

                sage: R([1, 3, 2])
                M[1, 3, 2]
                sage: R(Permutation([1, 3, 2]))
                M[1, 3, 2]
                sage: R(SymmetricGroup(4)(Permutation([1,3,4,2])))
                M[1, 3, 4, 2]

                sage: RF = algebras.FQSym(QQ).F()
                sage: R(RF([2, 3, 4, 1]))
                M[2, 3, 4, 1] + M[2, 4, 3, 1] + M[3, 2, 4, 1] + M[3, 4, 2, 1]
                 + M[4, 2, 3, 1] + M[4, 3, 2, 1]
                sage: R(RF([3, 2, 4, 1]))
                M[3, 2, 4, 1] + M[4, 2, 3, 1] + M[4, 3, 2, 1]
                sage: DF = algebras.FQSym(ZZ).F()
                sage: D(DF([2, 3, 4, 1]))
                M[2, 3, 4, 1] + M[2, 4, 3, 1] + M[3, 2, 4, 1] + M[3, 4, 2, 1]
                 + M[4, 2, 3, 1] + M[4, 3, 2, 1]
                sage: R(DF([2, 3, 4, 1]))
                M[2, 3, 4, 1] + M[2, 4, 3, 1] + M[3, 2, 4, 1] + M[3, 4, 2, 1]
                 + M[4, 2, 3, 1] + M[4, 3, 2, 1]
                sage: RF(R[2, 3, 4, 1])
                F[2, 3, 4, 1] - F[2, 4, 3, 1] - F[3, 2, 4, 1] + F[4, 3, 2, 1]

                sage: RG = algebras.FQSym(QQ).G()
                sage: R(RG([4, 1, 2, 3]))
                M[2, 3, 4, 1] + M[2, 4, 3, 1] + M[3, 2, 4, 1] + M[3, 4, 2, 1]
                 + M[4, 2, 3, 1] + M[4, 3, 2, 1]
                sage: R(RG([4, 2, 1, 3]))
                M[3, 2, 4, 1] + M[4, 2, 3, 1] + M[4, 3, 2, 1]
                sage: DG = algebras.FQSym(ZZ).G()
                sage: D(DG([4, 1, 2, 3]))
                M[2, 3, 4, 1] + M[2, 4, 3, 1] + M[3, 2, 4, 1] + M[3, 4, 2, 1]
                 + M[4, 2, 3, 1] + M[4, 3, 2, 1]
                sage: R(DG([4, 1, 2, 3]))
                M[2, 3, 4, 1] + M[2, 4, 3, 1] + M[3, 2, 4, 1] + M[3, 4, 2, 1]
                 + M[4, 2, 3, 1] + M[4, 3, 2, 1]
                sage: RG(R[2, 3, 4, 1])
                G[4, 1, 2, 3] - G[4, 1, 3, 2] - G[4, 2, 1, 3] + G[4, 3, 2, 1]
            """
            if isinstance(x, (list, tuple, PermutationGroupElement)):
                x = Permutation(x)
            try:
                P = x.parent()
                if isinstance(P, FreeQuasisymmetricFunctions.M):
                    if P is self:
                        return x
                    return self.element_class(self, x.monomial_coefficients())
            except AttributeError:
                pass
            return CombinatorialFreeModule._element_constructor_(self, x)

        def __getitem__(self, r):
            r"""
            The default implementation of ``__getitem__`` interprets
            the input as a tuple, which in case of permutations
            is interpreted as cycle notation, even though the input
            looks like a one-line notation.
            We override this method to amend this.

            EXAMPLES::

                sage: M = algebras.FQSym(QQ).M()
                sage: M[3, 2, 1]
                M[3, 2, 1]
                sage: M[1]
                M[1]
            """
            if isinstance(r, tuple):
                r = list(r)
            elif r == 1:
                r = [1]
            return super().__getitem__(r)

        def _F_to_M_on_basis(self, w):
            r"""
            Return `F_w` in terms of the M basis.

            INPUT:

            - ``w`` -- a permutation

            OUTPUT: an element of the M basis

            TESTS::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: F = FQSym.F()
                sage: M = FQSym.M()
                sage: M(F[3, 2, 1] - 4 * F[4, 2, 1, 3])
                M[3, 2, 1] - 4*M[4, 2, 1, 3] - 4*M[4, 3, 1, 2] - 4*M[4, 3, 2, 1]
                sage: all(M(M._F_to_M_on_basis(w)) == F[w] for i in range(5)
                ....:     for w in Permutations(i))
                True
                sage: F[3, 2, 1] == M[3, 2, 1]
                True
                sage: F[4, 2, 1, 3] == M[3, 2, 4, 1]
                False
            """
            return self.sum_of_monomials(w.permutohedron_greater(side='left'))

        def _M_to_F_on_basis(self, w):
            r"""
            Return `\mathcal{M}_w` in terms of the F basis.

            INPUT:

            - ``w`` -- a permutation

            OUTPUT: an element of the F basis

            ALGORITHM:

            If `w` is any permutation in `S_n`, then

            .. MATH::

                \mathcal{M}_w = \sum_u (-1)^{j(w, u)} F_u,

            where the sum ranges over all permutations `u \in S_n`
            obtained as follows:

            * Let `v = w^{-1}`.

            * Subdivide the list `(v(1), v(2), \ldots, v(n))` into
              an arbitrary number of nonempty blocks (by putting
              dividers between adjacent entries) in such a way that
              each block is strictly increasing (i.e., each descent
              of `v` is followed by a divider, but not every
              divider must necessarily follow a descent).

            * Reverse the order of entries in each block.

            * Remove the dividers. The resulting list is the
              one-line notation `(x(1), x(2), \ldots, x(n))` of
              some permutation `x \in S_n`.

            * Set `u = x^{-1}`. Also, let `j(w, u)` be `n` minus
              the number of blocks in our subdivision.

            This formula is equivalent to the formula (1.13) in
            [AguSot05]_, since Corollary 3.2.8 in [BB2005]_ expresses
            the Mobius function of the weak order.

            TESTS::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: F = FQSym.F()
                sage: M = FQSym.M()
                sage: F(M[3, 2, 1] - 4 * F[4, 2, 1, 3])
                F[3, 2, 1] - 4*F[4, 2, 1, 3]
                sage: all(F(M._M_to_F_on_basis(w)) == M[w] for i in range(5)  # indirect doctest
                ....:     for w in Permutations(i))
                True
                sage: all(M(F(M[w])) == M[w] for i in range(5)
                ....:     for w in Permutations(i))
                True
                sage: M[3, 2, 1] == F[3, 2, 1]
                True
                sage: M[3, 2, 4, 1] == F[4, 2, 1, 3]
                False
                sage: F(M[[]]) == F[[]]
                True
            """
            F = self.realization_of().F()
            if len(w) <= 1:
                return F.monomial(w)

            w_i = w.inverse()
            w_i = w_i[:]
            n = len(w_i)
            des = tuple([0] + [g for g in range(1, n)
                               if w_i[g - 1] > w_i[g]] + [n])
            non_des = [g for g in range(1, n) if w_i[g - 1] < w_i[g]]
            # Now, des is a list of all descents of w_i and also 0 and n,
            # whereas non_des is a list of all non-descents of w_i.

            Perms = self.basis().keys()

            R = self.base_ring()
            one = R.one()
            mine = -one

            dc = {w: one}
            from itertools import combinations
            for k in range(len(non_des)):
                kk = k + len(des)
                for extra_des in combinations(non_des, k):
                    breakpoints = sorted(des + extra_des)
                    # so that kk == len(breakpoints)
                    p = sum([w_i[breakpoints[g]: breakpoints[g + 1]][::-1]
                             for g in range(kk - 1)],
                            [])
                    u = Perms(p).inverse()
                    dc[u] = one if n % 2 != kk % 2 else mine

            return F._from_dict(dc)

        def degree_on_basis(self, t):
            """
            Return the degree of a permutation in
            the algebra of free quasi-symmetric functions.

            This is the size of the permutation (i.e., the `n`
            for which the permutation belongs to `S_n`).

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).M()
                sage: u = Permutation([2,1])
                sage: A.degree_on_basis(u)
                2
            """
            return len(t)

        def coproduct_on_basis(self, x):
            r"""
            Return the coproduct of `\mathcal{M}_{\sigma}` for `\sigma`
            a permutation (here, `\sigma` is ``x``).

            This uses Theorem 3.1 in [AguSot05]_.

            EXAMPLES::

                sage: M = algebras.FQSym(QQ).M()
                sage: x = M([1])
                sage: ascii_art(M.coproduct(M.one()))  # indirect doctest
                1 # 1

                sage: ascii_art(M.coproduct(x))  # indirect doctest
                1 # M    + M    # 1
                     [1]    [1]

                sage: M.coproduct(M([2, 1, 3]))
                M[] # M[2, 1, 3] + M[2, 1, 3] # M[]
                sage: M.coproduct(M([2, 3, 1]))
                M[] # M[2, 3, 1] + M[1, 2] # M[1] + M[2, 3, 1] # M[]
                sage: M.coproduct(M([3, 2, 1]))
                M[] # M[3, 2, 1] + M[1] # M[2, 1] + M[2, 1] # M[1]
                + M[3, 2, 1] # M[]
                sage: M.coproduct(M([3, 4, 2, 1]))
                M[] # M[3, 4, 2, 1] + M[1, 2] # M[2, 1] + M[2, 3, 1] # M[1]
                 + M[3, 4, 2, 1] # M[]
                sage: M.coproduct(M([3, 4, 1, 2]))
                M[] # M[3, 4, 1, 2] + M[1, 2] # M[1, 2] + M[3, 4, 1, 2] # M[]
            """
            n = len(x)
            if not n:
                return self.one().tensor(self.one())
            return sum(self(Word(x[:i]).standard_permutation()).tensor(
                self(Word(x[i:]).standard_permutation()))
                for i in range(n + 1)
                if (i == 0 or i == n or min(x[:i]) > max(x[i:])))

        class Element(FQSymBasis_abstract.Element):
            def star_involution(self):
                r"""
                Return the image of the element ``self`` of `FQSym`
                under the star involution.

                See
                :meth:`FQSymBases.ElementMethods.star_involution`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`omega_involution`, :meth:`psi_involution`

                EXAMPLES::

                    sage: FQSym = algebras.FQSym(ZZ)
                    sage: M = FQSym.M()
                    sage: M[[2,3,1]].star_involution()
                    M[3, 1, 2]
                    sage: M[[]].star_involution()
                    M[]

                TESTS::

                    sage: F = FQSym.F()
                    sage: all(M(F[w]).star_involution() == M(F[w].star_involution())
                    ....:     for w in Permutations(4))
                    True
                """
                # See the FQSymBases.ElementMethods.star_involution doc
                # for the formula we're using here.
                M = self.parent()
                return M._from_dict({w.complement().reverse(): c for (w, c) in self},
                                    remove_zeros=False)


class FQSymBases(Category_realization_of_parent):
    r"""
    The category of graded bases of `FQSym` indexed by permutations.
    """

    def __init__(self, base):
        r"""
        Initialize the bases of an `FQSym`.

        INPUT:

        - ``base`` -- an instance of `FQSym`

        TESTS::

            sage: from sage.combinat.fqsym import FQSymBases
            sage: FQSym = algebras.FQSym(ZZ)
            sage: bases = FQSymBases(FQSym)
            sage: FQSym.F() in bases
            True
        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.fqsym import FQSymBases
            sage: FQSym = algebras.FQSym(ZZ)
            sage: FQSymBases(FQSym)
            Category of bases of Free Quasi-symmetric functions over Integer Ring
        """
        return "Category of bases of {}".format(self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.fqsym import FQSymBases
            sage: FQSym = algebras.FQSym(ZZ)
            sage: bases = FQSymBases(FQSym)
            sage: bases.super_categories()
            [Category of realizations of Free Quasi-symmetric functions over Integer Ring,
             Join of Category of realizations of Hopf algebras over Integer Ring
               and Category of graded algebras over Integer Ring
               and Category of graded coalgebras over Integer Ring,
             Category of graded connected Hopf algebras with basis over Integer Ring]
        """
        R = self.base().base_ring()
        return [self.base().Realizations(),
                HopfAlgebras(R).Graded().Realizations(),
                HopfAlgebras(R).Graded().WithBasis().Graded().Connected(),
                ]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `FQSym`.

            EXAMPLES::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: FQSym.F()
                Free Quasi-symmetric functions over Integer Ring in the F basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, p):
            """
            Return the basis element indexed by ``p``.

            INPUT:

            - ``p`` -- a permutation

            EXAMPLES::

                sage: R = algebras.FQSym(QQ).F()
                sage: R[[1, 3, 2]]
                F[1, 3, 2]
                sage: R[Permutation([1, 3, 2])]
                F[1, 3, 2]
                sage: R[SymmetricGroup(4)(Permutation([1,3,4,2]))]
                F[1, 3, 4, 2]
            """
            return self.monomial(Permutation(p))

        def basis(self, degree=None):
            r"""
            The basis elements (optionally: of the specified degree).

            OUTPUT: Family

            EXAMPLES::

                sage: FQSym = algebras.FQSym(QQ)
                sage: G = FQSym.G()
                sage: G.basis()
                Lazy family (Term map from Standard permutations to Free Quasi-symmetric functions over Rational Field in the G basis(i))_{i in Standard permutations}
                sage: G.basis().keys()
                Standard permutations
                sage: G.basis(degree=3).keys()
                Standard permutations of 3
                sage: G.basis(degree=3).list()
                [G[1, 2, 3], G[1, 3, 2], G[2, 1, 3], G[2, 3, 1], G[3, 1, 2], G[3, 2, 1]]
            """
            from sage.sets.family import Family
            if degree is None:
                return Family(self._indices, self.monomial)
            else:
                return Family(Permutations(degree), self.monomial)

        def is_field(self, proof=True):
            """
            Return whether this `FQSym` is a field.

            EXAMPLES::

                sage: F = algebras.FQSym(QQ).F()
                sage: F.is_field()
                False
            """
            return False

        def is_commutative(self):
            """
            Return whether this `FQSym` is commutative.

            EXAMPLES::

                sage: F = algebras.FQSym(ZZ).F()
                sage: F.is_commutative()
                False
            """
            return self.base_ring().is_zero()

        def some_elements(self):
            """
            Return some elements of the free quasi-symmetric functions.

            EXAMPLES::

                sage: A = algebras.FQSym(QQ)
                sage: F = A.F()
                sage: F.some_elements()
                [F[], F[1], F[1, 2] + F[2, 1], F[] + F[1, 2] + F[2, 1]]
                sage: G = A.G()
                sage: G.some_elements()
                [G[], G[1], G[1, 2] + G[2, 1], G[] + G[1, 2] + G[2, 1]]
                sage: M = A.M()
                sage: M.some_elements()
                [M[], M[1], M[1, 2] + 2*M[2, 1], M[] + M[1, 2] + 2*M[2, 1]]
            """
            u = self.one()
            o = self.monomial(Permutation([1]))
            x = o * o
            y = u + x
            return [u, o, x, y]

        @cached_method
        def one_basis(self):
            """
            Return the index of the unit.

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: A.one_basis()
                []
            """
            Perm = self.basis().keys()
            return Perm([])

        @lazy_attribute
        def succ(self):
            r"""
            Return the `\succ` product.

            On the F-basis of ``FQSym``, this product is determined by
            `F_x \succ F_y = \sum F_z`, where the sum ranges over all `z`
            in the shifted shuffle of `x` and `y` with the additional
            condition that the first letter of the result comes from `y`.

            The usual symbol for this operation is `\succ`.

            .. SEEALSO::

                :meth:`~sage.categories.magmas.Magmas.ParentMethods.product`,
                :meth:`prec`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = A([1])
                sage: A.succ(x, x)
                F[2, 1]
                sage: y = A([3,1,2])
                sage: A.succ(x, y)
                F[4, 1, 2, 3] + F[4, 2, 1, 3] + F[4, 2, 3, 1]
                sage: A.succ(y, x)
                F[4, 3, 1, 2]
            """
            try:
                suc = self.succ_product_on_basis
            except AttributeError:
                return self.succ_by_coercion
            return self._module_morphism(self._module_morphism(suc, position=0,
                                                               codomain=self),
                                         position=1)

        def succ_by_coercion(self, x, y):
            r"""
            Return `x \succ y`, computed using coercion to the F-basis.

            See :meth:`succ` for the definition of the objects involved.

            EXAMPLES::

                sage: G = algebras.FQSym(ZZ).G()
                sage: G.succ(G([1]), G([2, 3, 1])) # indirect doctest
                G[2, 3, 4, 1] + G[3, 2, 4, 1] + G[4, 2, 3, 1]
            """
            F = self.realization_of().a_realization()
            return self(F.succ(F(x), F(y)))

        @lazy_attribute
        def prec(self):
            r"""
            Return the `\prec` product.

            On the F-basis of ``FQSym``, this product is determined by
            `F_x \prec F_y = \sum F_z`, where the sum ranges over all `z`
            in the shifted shuffle of `x` and `y` with the additional
            condition that the first letter of the result comes from `x`.

            The usual symbol for this operation is `\prec`.

            .. SEEALSO::

                :meth:`~sage.categories.magmas.Magmas.ParentMethods.product`,
                :meth:`succ`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = A([2,1])
                sage: A.prec(x, x)
                F[2, 1, 4, 3] + F[2, 4, 1, 3] + F[2, 4, 3, 1]
                sage: y = A([2,1,3])
                sage: A.prec(x, y)
                F[2, 1, 4, 3, 5] + F[2, 4, 1, 3, 5] + F[2, 4, 3, 1, 5]
                 + F[2, 4, 3, 5, 1]
                sage: A.prec(y, x)
                F[2, 1, 3, 5, 4] + F[2, 1, 5, 3, 4] + F[2, 1, 5, 4, 3]
                 + F[2, 5, 1, 3, 4] + F[2, 5, 1, 4, 3] + F[2, 5, 4, 1, 3]
            """
            try:
                pre = self.prec_product_on_basis
            except AttributeError:
                return self.prec_by_coercion
            return self._module_morphism(self._module_morphism(pre, position=0,
                                                               codomain=self),
                                         position=1)

        def prec_by_coercion(self, x, y):
            r"""
            Return `x \prec y`, computed using coercion to the F-basis.

            See :meth:`prec` for the definition of the objects involved.

            EXAMPLES::

                sage: G = algebras.FQSym(ZZ).G()
                sage: a = G([1])
                sage: b = G([2, 3, 1])
                sage: G.prec(a, b) + G.succ(a, b) == a * b # indirect doctest
                True
            """
            F = self.realization_of().a_realization()
            return self(F.prec(F(x), F(y)))

        def from_symmetric_group_algebra(self, x):
            """
            Return the element of `FQSym` corresponding to the element
            `x` of a symmetric group algebra.

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: SGA4 = SymmetricGroupAlgebra(QQ, 4)
                sage: x = SGA4([1,3,2,4]) + 5/2 * SGA4([1,2,4,3])
                sage: A.from_symmetric_group_algebra(x)
                5/2*F[1, 2, 4, 3] + F[1, 3, 2, 4]
                sage: A.from_symmetric_group_algebra(SGA4.zero())
                0
            """
            return self._from_dict({Permutation(key): c for (key, c) in x})

    class ElementMethods:
        def omega_involution(self):
            r"""
            Return the image of the element ``self`` of `FQSym`
            under the omega involution.

            The `\omega` involution is defined as the
            linear map `FQSym \to FQSym` that sends each basis
            element `F_u` of the F-basis of `FQSym`
            to the basis element `F_{u \circ w_0}`, where `w_0` is
            the longest word (i.e., `w_0(i) = n + 1 - i`) in the
            symmetric group `S_n` that contains `u`. The `\omega`
            involution is a graded algebra automorphism and a
            coalgebra anti-automorphism of `FQSym`. Every
            permutation `u \in S_n` satisfies

            .. MATH::

                \omega(F_u) = F_{u \circ w_0}, \qquad
                \omega(G_u) = G_{w_0 \circ u},

            where standard notations for classical bases of `FQSym`
            are being used (that is, `F` for the F-basis, and
            `G` for the G-basis).
            In other words, writing permutations in one-line notation,
            we have

            .. MATH::

                \omega(F_{(u_1, u_2, \ldots, u_n)})
                = F_{(u_n, u_{n-1}, \ldots, u_1)}, \qquad
                \omega(G_{(u_1, u_2, \ldots, u_n)})
                = G_{(n+1-u_1, n+1-u_2, \ldots, n+1-u_n)}.

            If we also consider the `\omega` involution
            (:meth:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.omega_involution`)
            of the quasisymmetric functions (by slight abuse
            of notation), and if we let `\pi` be the canonical
            projection `FQSym \to QSym`, then
            `\pi \circ \omega = \omega \circ \pi`.

            Additionally, consider the `\psi` involution
            (:meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution`)
            of the noncommutative symmetric functions, and if we let
            `\iota` be the canonical inclusion `NSym \to FQSym`,
            then `\omega \circ \iota = \iota \circ \psi`.

            .. TODO::

                Duality?

            .. SEEALSO::

                :meth:`psi_involution`, :meth:`star_involution`

            EXAMPLES::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: F = FQSym.F()
                sage: F[[2,3,1]].omega_involution()
                F[1, 3, 2]
                sage: (3*F[[1]] - 4*F[[]] + 5*F[[1,2]]).omega_involution()
                -4*F[] + 3*F[1] + 5*F[2, 1]
                sage: G = FQSym.G()
                sage: G[[2,3,1]].omega_involution()
                G[2, 1, 3]
                sage: M = FQSym.M()
                sage: M[[2,3,1]].omega_involution()
                -M[1, 2, 3] - M[2, 1, 3] - M[3, 1, 2]

            The omega involution is an algebra homomorphism::

                sage: (F[1,2] * F[1]).omega_involution()
                F[2, 1, 3] + F[2, 3, 1] + F[3, 2, 1]
                sage: F[1,2].omega_involution() * F[1].omega_involution()
                F[2, 1, 3] + F[2, 3, 1] + F[3, 2, 1]

            The omega involution intertwines the antipode
            and the inverse of the antipode::

                sage: all( F(I).antipode().omega_involution().antipode()
                ....:      == F(I).omega_involution()
                ....:      for I in Permutations(4) )
                True

            Testing the `\pi \circ \omega = \omega \circ \pi` relation
            noticed above::

                sage: all( M[I].omega_involution().to_qsym()
                ....:      == M[I].to_qsym().omega_involution()
                ....:      for I in Permutations(4) )
                True

            Testing the `\omega \circ \iota = \iota \circ \psi` relation::

                sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                sage: S = NSym.S()
                sage: all( S[I].psi_involution().to_fqsym() == S[I].to_fqsym().omega_involution()
                ....:      for I in Compositions(4) )
                True

            .. TODO::

                Check further commutative squares.
            """
            # Convert to the F-basis, there apply the reversal
            # componentwise, then convert back.
            parent = self.parent()
            F = parent.realization_of().F()
            dct = {I.reverse(): coeff for (I, coeff) in F(self)}
            return parent(F._from_dict(dct, remove_zeros=False))

        def psi_involution(self):
            r"""
            Return the image of the element ``self`` of `FQSym`
            under the psi involution.

            The `\psi` involution is defined as the
            linear map `FQSym \to FQSym` that sends each basis
            element `F_u` of the F-basis of `FQSym`
            to the basis element `F_{w_0 \circ u}`, where `w_0` is
            the longest word (i.e., `w_0(i) = n + 1 - i`) in the
            symmetric group `S_n` that contains `u`. The `\psi`
            involution is a graded coalgebra automorphism and
            an algebra anti-automorphism of `FQSym`. Every
            permutation `u \in S_n` satisfies

            .. MATH::

                \psi(F_u) = F_{w_0 \circ u}, \qquad
                \psi(G_u) = G_{u \circ w_0},

            where standard notations for classical bases of `FQSym`
            are being used (that is, `F` for the F-basis, and
            `G` for the G-basis). In other words, writing
            permutations in one-line notation, we have

            .. MATH::

                \psi(F_{(u_1, u_2, \ldots, u_n)})
                = F_{(n+1-u_1, n+1-u_2, \ldots, n+1-u_n)}, \qquad
                \psi(G_{(u_1, u_2, \ldots, u_n)})
                = G_{(u_n, u_{n-1}, \ldots, u_1)}.

            If we also consider the `\psi` involution
            (:meth:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.psi_involution`)
            of the quasisymmetric functions (by slight abuse of
            notation), and if we let `\pi` be the canonical
            projection `FQSym \to QSym`, then
            `\pi \circ \psi = \psi \circ \pi`.

            Additionally, consider the `\omega` involution
            (:meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.omega_involution`)
            of the noncommutative symmetric functions, and if we let
            `\iota` be the canonical inclusion `NSym \to FQSym`,
            then `\psi \circ \iota = \iota \circ \omega`.

            .. TODO::

                Duality?

            .. SEEALSO::

                :meth:`omega_involution`, :meth:`star_involution`

            EXAMPLES::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: F = FQSym.F()
                sage: F[[2,3,1]].psi_involution()
                F[2, 1, 3]
                sage: (3*F[[1]] - 4*F[[]] + 5*F[[1,2]]).psi_involution()
                -4*F[] + 3*F[1] + 5*F[2, 1]
                sage: G = FQSym.G()
                sage: G[[2,3,1]].psi_involution()
                G[1, 3, 2]
                sage: M = FQSym.M()
                sage: M[[2,3,1]].psi_involution()
                -M[1, 2, 3] - M[1, 3, 2] - M[2, 3, 1]

            The `\psi` involution intertwines the antipode
            and the inverse of the antipode::

                sage: all( F(I).antipode().psi_involution().antipode()
                ....:      == F(I).psi_involution()
                ....:      for I in Permutations(4) )
                True

            Testing the `\pi \circ \psi = \psi \circ \pi` relation above::

                sage: all( M[I].psi_involution().to_qsym()
                ....:      == M[I].to_qsym().psi_involution()
                ....:      for I in Permutations(4) )
                True

            Testing the `\psi \circ \iota = \iota \circ \omega` relation::

                sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                sage: S = NSym.S()
                sage: all( S[I].omega_involution().to_fqsym() == S[I].to_fqsym().psi_involution()
                ....:      for I in Compositions(4) )
                True

            .. TODO::

                Check further commutative squares.
            """
            # Convert to the F-basis, there apply the complement
            # componentwise, then convert back.
            parent = self.parent()
            F = parent.realization_of().F()
            dct = {I.complement(): coeff for (I, coeff) in F(self)}
            return parent(F._from_dict(dct, remove_zeros=False))

        def star_involution(self):
            r"""
            Return the image of the element ``self`` of `FQSym`
            under the star involution.

            The star involution is defined as the
            linear map `FQSym \to FQSym` that sends each basis
            element `F_u` of the F-basis of `FQSym`
            to the basis element `F_{w_0 \circ u \circ w_0}`, where
            `w_0` is the longest word (i.e., `w_0(i) = n + 1 - i`)
            in the symmetric group `S_n` that contains `u`.
            The star involution is a graded Hopf algebra
            anti-automorphism of `FQSym`.
            It is denoted by `f \mapsto f^*`. Every permutation
            `u \in S_n` satisfies

            .. MATH::

                (F_u)^* = F_{w_0 \circ u \circ w_0}, \qquad
                (G_u)^* = G_{w_0 \circ u \circ w_0}, \qquad
                (\mathcal{M}_u)^* = \mathcal{M}_{w_0 \circ u \circ w_0},

            where standard notations for classical bases of `FQSym`
            are being used (that is, `F` for the F-basis,
            `G` for the G-basis, and `\mathcal{M}` for the Monomial
            basis). In other words, writing permutations in one-line
            notation, we have

            .. MATH::

                (F_{(u_1, u_2, \ldots, u_n)})^*
                = F_{(n+1-u_n, n+1-u_{n-1}, \ldots, n+1-u_1)}, \qquad
                (G_{(u_1, u_2, \ldots, u_n)})^*
                = G_{(n+1-u_n, n+1-u_{n-1}, \ldots, n+1-u_1)},

            and

            .. MATH::

                (\mathcal{M}_{(u_1, u_2, \ldots, u_n)})^*
                = \mathcal{M}_{(n+1-u_n, n+1-u_{n-1}, \ldots, n+1-u_1)}.

            Let us denote the star involution by `(\ast)` as well.

            If we also denote by `(\ast)` the star involution of
            of the quasisymmetric functions
            (:meth:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution`)
            and if we let `\pi : FQSym \to QSym` be the canonical
            projection then `\pi \circ (\ast) = (\ast) \circ \pi`.
            Similar for the noncommutative symmetric functions
            (:meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution`)
            with `\pi : NSym \to FQSym` being the canonical inclusion
            and the word quasisymmetric functions
            (:meth:`~sage.combinat.chas.wqsym.WordQuasiSymmetricFunctions.Bases.ElementMethods.star_involution`)
            with `\pi : FQSym \to WQSym` the canonical inclusion.

            .. TODO::

                Duality?

            .. SEEALSO::

                :meth:`omega_involution`, :meth:`psi_involution`

            EXAMPLES::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: F = FQSym.F()
                sage: F[[2,3,1]].star_involution()
                F[3, 1, 2]
                sage: (3*F[[1]] - 4*F[[]] + 5*F[[1,2]]).star_involution()
                -4*F[] + 3*F[1] + 5*F[1, 2]
                sage: G = FQSym.G()
                sage: G[[2,3,1]].star_involution()
                G[3, 1, 2]
                sage: M = FQSym.M()
                sage: M[[2,3,1]].star_involution()
                M[3, 1, 2]

            The star involution commutes with the antipode::

                sage: all( F(I).antipode().star_involution()
                ....:      == F(I).star_involution().antipode()
                ....:      for I in Permutations(4) )
                True

            Testing the `\pi \circ (\ast) = (\ast) \circ \pi` relation::

                sage: all( M[I].star_involution().to_qsym()
                ....:      == M[I].to_qsym().star_involution()
                ....:      for I in Permutations(4) )
                True

            Similar for `NSym`::

                sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                sage: S = NSym.S()
                sage: all( S[I].star_involution().to_fqsym() == S[I].to_fqsym().star_involution()
                ....:      for I in Compositions(4) )
                True

            Similar for `WQSym`::

                sage: WQSym = algebras.WQSym(ZZ)
                sage: all( F(I).to_wqsym().star_involution()
                ....:      == F(I).star_involution().to_wqsym()
                ....:      for I in Permutations(4) )
                True

            .. TODO::

                Check further commutative squares.
            """
            # Convert to the F-basis, there apply the reversal and
            # complement componentwise, then convert back.
            parent = self.parent()
            F = parent.realization_of().F()
            dct = {I.complement().reverse(): coeff for (I, coeff) in F(self)}
            return parent(F._from_dict(dct, remove_zeros=False))

        def to_symmetric_group_algebra(self, n=None):
            """
            Return the element of a symmetric group algebra
            corresponding to the element ``self`` of `FQSym`.

            INPUT:

            - ``n`` -- integer (default: the maximal degree of ``self``);
              the rank of the target symmetric group algebra

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).G()
                sage: x = A([1,3,2,4]) + 5/2 * A([2,3,4,1])
                sage: x.to_symmetric_group_algebra()
                [1, 3, 2, 4] + 5/2*[4, 1, 2, 3]
            """
            F = self.parent().realization_of().F()
            return F(self).to_symmetric_group_algebra(n=n)

        def to_wqsym(self):
            r"""
            Return the image of ``self`` under the canonical
            inclusion map `FQSym \to WQSym`.

            The canonical inclusion map `FQSym \to WQSym` is
            an injective homomorphism of Hopf algebras. It sends
            a basis element `G_w` of `FQSym` to the sum of
            basis elements `\mathbf{M}_u` of `WQSym`, where `u`
            ranges over all packed words whose standardization
            is `w`.

            .. SEEALSO::

                :class:`WordQuasiSymmetricFunctions` for a
                definition of `WQSym`.

            EXAMPLES::

                sage: G = algebras.FQSym(QQ).G()
                sage: x = G[1, 3, 2]
                sage: x.to_wqsym()
                M[{1}, {3}, {2}] + M[{1, 3}, {2}]
                sage: G[1, 2].to_wqsym()
                M[{1}, {2}] + M[{1, 2}]
                sage: F = algebras.FQSym(QQ).F()
                sage: F[3, 1, 2].to_wqsym()
                M[{3}, {1}, {2}] + M[{3}, {1, 2}]
                sage: G[2, 3, 1].to_wqsym()
                M[{3}, {1}, {2}] + M[{3}, {1, 2}]
            """
            parent = self.parent()
            FQSym = parent.realization_of()
            G = FQSym.G()
            from sage.combinat.chas.wqsym import WordQuasiSymmetricFunctions
            M = WordQuasiSymmetricFunctions(parent.base_ring()).M()
            OSP = M.basis().keys()
            from sage.combinat.words.finite_word import word_to_ordered_set_partition

            def to_wqsym_on_G_basis(w):
                # Return the image of `G_w` under the inclusion
                # map `FQSym \to WQSym`.
                dc = w.inverse().descents_composition()
                res = M.zero()
                for comp in dc.finer():
                    v = w.destandardize(comp)
                    res += M[OSP(word_to_ordered_set_partition(v))]
                return res
            return M.linear_combination((to_wqsym_on_G_basis(w), coeff)
                                        for w, coeff in G(self))

        def to_qsym(self):
            r"""
            Return the image of ``self`` under the canonical
            projection `FQSym \to QSym`.

            The canonical projection `FQSym \to QSym` is a
            surjective homomorphism of Hopf algebras. It sends a
            basis element `F_w` of `FQSym` to the basis element
            `F_{\operatorname{Comp} w}` of the fundamental basis
            of `QSym`, where `\operatorname{Comp} w` stands for
            the descent composition
            (:meth:`sage.combinat.permutation.Permutation.descents_composition`)
            of the permutation `w`.

            .. SEEALSO::

                :class:`QuasiSymmetricFunctions` for a
                definition of `QSym`.

            EXAMPLES::

                sage: G = algebras.FQSym(QQ).G()
                sage: x = G[1, 3, 2]
                sage: x.to_qsym()
                F[2, 1]
                sage: G[2, 3, 1].to_qsym()
                F[1, 2]
                sage: F = algebras.FQSym(QQ).F()
                sage: F[2, 3, 1].to_qsym()
                F[2, 1]
                sage: (F[2, 3, 1] + F[1, 3, 2] + F[1, 2, 3]).to_qsym()
                2*F[2, 1] + F[3]
                sage: F2 = algebras.FQSym(GF(2)).F()
                sage: F2[2, 3, 1].to_qsym()
                F[2, 1]
                sage: (F2[2, 3, 1] + F2[1, 3, 2] + F2[1, 2, 3]).to_qsym()
                F[3]
            """
            parent = self.parent()
            FQSym = parent.realization_of()
            F = FQSym.F()
            from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
            QF = QuasiSymmetricFunctions(parent.base_ring()).F()
            return QF.sum_of_terms((w.descents_composition(), coeff)
                                   for w, coeff in F(self))
