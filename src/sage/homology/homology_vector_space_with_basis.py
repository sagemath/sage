# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.graphs          (because all doctests use the catalogs simplicial_complexes, cubical_complexes)
"""
Homology and cohomology with a basis

This module provides homology and cohomology vector spaces suitable
for computing cup products and cohomology operations.

REFERENCES:

- [GDR2003]_
- [GDR1999]_

AUTHORS:

- John H. Palmieri, Travis Scrimshaw (2015-09)
"""
########################################################################
#       Copyright (C) 2015 John H. Palmieri <palmieri@math.washington.edu>
#                          Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
########################################################################

from sage.algebras.steenrod.steenrod_algebra import SteenrodAlgebra
from sage.misc.cachefunc import cached_method
from sage.categories.algebras import Algebras
from sage.categories.category import Category
from sage.categories.left_modules import LeftModules
from sage.categories.right_modules import RightModules
from sage.categories.modules import Modules
from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.sets.family import Family

try:
    from sage.topology.simplicial_complex import SimplicialComplex
    from sage.topology.simplicial_set import SimplicialSet_arbitrary
    from sage.topology.cubical_complex import CubicalComplex
except ImportError:
    SimplicialComplex = SimplicialSet_arbitrary = CubicalComplex = ()


class HomologyVectorSpaceWithBasis(CombinatorialFreeModule):
    r"""
    Homology (or cohomology) vector space.

    This provides enough structure to allow the computation of cup
    products and cohomology operations. See the class
    :class:`CohomologyRing` (which derives from this) for examples.

    It also requires field coefficients (hence the "VectorSpace" in
    the name of the class).

    .. NOTE::

        This is not intended to be created directly by the user, but
        instead via the methods
        :meth:`~sage.topology.cell_complex.GenericCellComplex.homology_with_basis` and
        :meth:`~sage.topology.cell_complex.GenericCellComplex.cohomology_ring`
        for the class of :class:`cell
        complexes<sage.topology.cell_complex.GenericCellComplex>`.

    INPUT:

    - ``base_ring`` -- must be a field
    - ``cell_complex`` -- the cell complex whose homology we are
      computing
    - ``cohomology`` -- boolean (default: ``False``); if ``True``, return
      the cohomology as a module
    - ``category`` -- (optional) a subcategory of modules with basis

    EXAMPLES:

    Homology classes are denoted by ``h_{d,i}`` where ``d`` is the
    degree of the homology class and ``i`` is their index in the list
    of basis elements in that degree. Cohomology classes are denoted
    ``h^{1,0}``::

        sage: RP2 = cubical_complexes.RealProjectivePlane()
        sage: RP2.homology_with_basis(GF(2))
        Homology module of Cubical complex with 21 vertices and 81 cubes
         over Finite Field of size 2
        sage: RP2.cohomology_ring(GF(2))
        Cohomology ring of Cubical complex with 21 vertices and 81 cubes
         over Finite Field of size 2
        sage: simplicial_complexes.Torus().homology_with_basis(QQ)
        Homology module of Minimal triangulation of the torus
         over Rational Field

    To access a basis element, use its degree and index (0 or 1 in the 1st
    cohomology group of a torus)::

        sage: H = simplicial_complexes.Torus().cohomology_ring(QQ)
        sage: H.basis(1)
        Finite family {(1, 0): h^{1,0}, (1, 1): h^{1,1}}
        sage: x = H.basis()[1,0]; x
        h^{1,0}
        sage: y = H.basis()[1,1]; y
        h^{1,1}
        sage: 2*x-3*y
        2*h^{1,0} - 3*h^{1,1}

    You can compute cup products of cohomology classes::

        sage: x.cup_product(y)
        -h^{2,0}
        sage: y.cup_product(x)
        h^{2,0}
        sage: x.cup_product(x)
        0

    This works with simplicial, cubical, and `\Delta`-complexes, and
    also simplicial sets::

        sage: Torus_c = cubical_complexes.Torus()
        sage: H = Torus_c.cohomology_ring(GF(2))
        sage: x,y = H.basis(1)
        sage: x.cup_product(x)
        0
        sage: x.cup_product(y)
        h^{2,0}
        sage: y.cup_product(y)
        0

        sage: Klein_d = delta_complexes.KleinBottle()
        sage: H = Klein_d.cohomology_ring(GF(2))
        sage: u,v = sorted(H.basis(1))
        sage: u.cup_product(u)
        h^{2,0}
        sage: u.cup_product(v)
        0
        sage: v.cup_product(v)
        h^{2,0}

    An isomorphism between the rings for the cubical model and the
    `\Delta`-complex model can be obtained by sending `x` to `u+v`,
    `y` to `v`. ::

        sage: # needs sage.groups
        sage: X = simplicial_sets.RealProjectiveSpace(6)
        sage: H_X = X.cohomology_ring(GF(2))
        sage: a = H_X.basis()[1,0]
        sage: a**6
        h^{6,0}
        sage: a**7
        0

    All products of positive-dimensional elements in a suspension
    should be zero::

        sage: # needs sage.groups
        sage: Y = X.suspension()
        sage: H_Y = Y.cohomology_ring(GF(2))
        sage: b = H_Y.basis()[2,0]
        sage: b**2
        0
        sage: B = sorted(H_Y.basis())[1:]
        sage: B
        [h^{2,0}, h^{3,0}, h^{4,0}, h^{5,0}, h^{6,0}, h^{7,0}]
        sage: import itertools
        sage: [a*b for (a,b) in itertools.combinations(B, 2)]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    The basis elements in the simplicial complex case have been chosen
    differently; apply the change of basis `x \mapsto a + b`, `y \mapsto
    b` to see the same product structure. ::

        sage: Klein_s = simplicial_complexes.KleinBottle()
        sage: H = Klein_s.cohomology_ring(GF(2))
        sage: a,b = H.basis(1)
        sage: a.cup_product(a)
        0
        sage: a.cup_product(b)
        h^{2,0}
        sage: (a+b).cup_product(a+b)
        h^{2,0}
        sage: b.cup_product(b)
        h^{2,0}
    """
    def __init__(self, base_ring, cell_complex, cohomology=False, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RP2 = simplicial_complexes.ProjectivePlane()
            sage: H = RP2.homology_with_basis(QQ)
            sage: TestSuite(H).run()
            sage: H = RP2.homology_with_basis(GF(2))
            sage: TestSuite(H).run()
            sage: H = RP2.cohomology_ring(GF(2))
            sage: TestSuite(H).run()
            sage: H = RP2.cohomology_ring(GF(5))
            sage: TestSuite(H).run()
            sage: H = simplicial_complexes.ComplexProjectivePlane().cohomology_ring()
            sage: TestSuite(H).run()
        """
        # phi is the associated chain contraction.
        # M is the homology chain complex.
        phi, M = cell_complex.algebraic_topological_model(base_ring)
        if cohomology:
            phi = phi.dual()
            # We only need the rank of M in each degree, and since
            # we're working over a field, we don't need to dualize M
            # if working with cohomology.
        category = Modules(base_ring).WithBasis().Graded().FiniteDimensional().or_subcategory(category)
        self._contraction = phi
        self._complex = cell_complex
        self._cohomology = cohomology
        self._graded_indices = {deg: range(M.free_module_rank(deg))
                                for deg in range(cell_complex.dimension()+1)}
        indices = [(deg, i) for deg in self._graded_indices
                   for i in self._graded_indices[deg]]
        CombinatorialFreeModule.__init__(self, base_ring, indices, category=category)

    def basis(self, d=None):
        """
        Return (the degree ``d`` homogeneous component of) the basis
        of this graded vector space.

        INPUT:

        - ``d`` -- (optional) the degree

        EXAMPLES::

            sage: RP2 = simplicial_complexes.ProjectivePlane()
            sage: H = RP2.homology_with_basis(QQ)
            sage: H.basis()
            Finite family {(0, 0): h_{0,0}}
            sage: H.basis(0)
            Finite family {(0, 0): h_{0,0}}
            sage: H.basis(1)
            Finite family {}
            sage: H.basis(2)
            Finite family {}
        """
        if d is None:
            return Family(self._indices, self.monomial)
        else:
            indices = [(d, i) for i in self._graded_indices.get(d, [])]
            return Family(indices, self.monomial)

    def degree_on_basis(self, i):
        r"""
        Return the degree of the basis element indexed by ``i``.

        EXAMPLES::

            sage: H = simplicial_complexes.Torus().homology_with_basis(GF(7))
            sage: H.degree_on_basis((2,0))
            2
        """
        return i[0]

    def contraction(self):
        r"""
        The chain contraction associated to this homology computation.

        That is, to work with chain representatives of homology
        classes, we need the chain complex `C` associated to the cell
        complex, the chain complex `H` of its homology (with trivial
        differential), chain maps `\pi: C \to H` and `\iota: H \to C`,
        and a chain contraction `\phi` giving a chain homotopy between
        `1_C` and `\iota \circ \pi`.

        OUTPUT: `\phi`

        See :class:`~sage.homology.chain_homotopy.ChainContraction` for information
        about chain contractions, and see
        :func:`~sage.homology.algebraic_topological_model.algebraic_topological_model`
        for the construction of this particular chain contraction `\phi`.

        EXAMPLES::

            sage: H = simplicial_complexes.Simplex(2).homology_with_basis(QQ)
            sage: H.contraction()
            Chain homotopy between:
              Chain complex endomorphism of Chain complex with at most 3 nonzero terms over Rational Field
              and Chain complex endomorphism of Chain complex with at most 3 nonzero terms over Rational Field

        From the chain contraction, one can also recover the maps `\pi`
        and `\iota`::

            sage: phi = H.contraction()
            sage: phi.pi()
            Chain complex morphism:
              From: Chain complex with at most 3 nonzero terms over Rational Field
              To: Chain complex with at most 1 nonzero terms over Rational Field
            sage: phi.iota()
            Chain complex morphism:
              From: Chain complex with at most 1 nonzero terms over Rational Field
              To: Chain complex with at most 3 nonzero terms over Rational Field
        """
        return self._contraction

    def complex(self):
        """
        The cell complex whose homology is being computed.

        EXAMPLES::

            sage: H = simplicial_complexes.Simplex(2).homology_with_basis(QQ)
            sage: H.complex()
            The 2-simplex
        """
        return self._complex

    def _repr_(self):
        """
        EXAMPLES::

            sage: simplicial_complexes.Torus().homology_with_basis(QQ)
            Homology module of Minimal triangulation of the torus
             over Rational Field
        """
        if self._cohomology:
            base = "Cohomology"
        else:
            base = "Homology"
        return base + " module of {} over {}".format(self._complex, self.base_ring())

    def _repr_term(self, i):
        """
        Return ``'h_{i[0],i[1]}'`` for homology, ``'h^{i[0],i[1]}'`` for
        cohomology, for the basis element indexed by ``i``.

        EXAMPLES::

            sage: H = simplicial_complexes.Torus().homology_with_basis(QQ)
            sage: H.basis()[1,0] # indirect doctest
            h_{1,0}
            sage: latex(H.basis()[1,1]) # indirect doctest
            h_{1,1}
            sage: co = simplicial_complexes.KleinBottle().cohomology_ring(GF(2))
            sage: co.basis()[1,0] # indirect doctest
            h^{1,0}
        """
        sym = '^' if self._cohomology else '_'
        return 'h{}{{{},{}}}'.format(sym, i[0], i[1])

    _latex_term = _repr_term

    @cached_method
    def _to_cycle_on_basis(self, i):
        r"""
        Return the (co)cycle representative of the basis element
        indexed by ``i``.

        .. SEEALSO::

            :meth:`HomologyVectorSpaceWithBasis.Element.to_cocycle`

        EXAMPLES::

            sage: S2 = simplicial_complexes.Sphere(2)
            sage: H = S2.homology_with_basis(QQ)
            sage: H._to_cycle_on_basis((2,0))
            -(0, 1, 2) + (0, 1, 3) - (0, 2, 3) + (1, 2, 3)

            sage: S2.cohomology_ring(QQ)._to_cycle_on_basis((2,0))
            \chi_(1, 2, 3)
            sage: S2.cohomology_ring(QQ)._to_cycle_on_basis((0,0))
            \chi_(0,) + \chi_(1,) + \chi_(2,) + \chi_(3,)

            sage: RP3 = simplicial_complexes.RealProjectiveSpace(3)
            sage: H = RP3.cohomology_ring(GF(2))
            sage: H._to_cycle_on_basis((0,0))
            \chi_(1,) + \chi_(2,) + \chi_(3,) + \chi_(4,) + \chi_(5,) + \chi_(6,)
             + \chi_(7,) + \chi_(8,) + \chi_(9,) + \chi_(10,) + \chi_(11,)
            sage: H._to_cycle_on_basis((1,0))
            \chi_(2, 4) + \chi_(2, 5) + \chi_(2, 8) + \chi_(2, 10)
            + \chi_(3, 4) + \chi_(3, 6) + \chi_(3, 8) + \chi_(3, 9)
            + \chi_(4, 5) + \chi_(4, 6) + \chi_(4, 11) + \chi_(5, 7)
            + \chi_(5, 9) + \chi_(6, 7) + \chi_(6, 10) + \chi_(7, 8)
            + \chi_(9, 10)
            sage: H._to_cycle_on_basis((2,0))
            \chi_(3, 5, 9) + \chi_(3, 6, 10) + \chi_(3, 9, 10)
            + \chi_(4, 5, 7) + \chi_(4, 5, 9) + \chi_(4, 6, 7) + \chi_(6, 7, 10)
            sage: H._to_cycle_on_basis((3,0))
            \chi_(5, 6, 7, 8)
        """
        vec = self.contraction().iota().in_degree(i[0]).column(i[1])
        chains = self.complex().n_chains(i[0], self.base_ring(),
                                         cochains=self._cohomology)
        return chains.from_vector(vec)

    def dual(self):
        r"""
        Return the dual space.

        If ``self`` is homology, return the cohomology ring. If
        ``self`` is cohomology, return the homology as a vector space.

        EXAMPLES::

            sage: T = simplicial_complexes.Torus()
            sage: hom = T.homology_with_basis(GF(2))
            sage: coh = T.cohomology_ring(GF(2))
            sage: hom.dual() is coh
            True
            sage: coh.dual() is hom
            True
        """
        if is_GF2(self.base_ring()):
            if self._cohomology:
                return HomologyVectorSpaceWithBasis_mod2(self.base_ring(),
                                                         self.complex())
            return CohomologyRing_mod2(self.base_ring(), self.complex())
        if self._cohomology:
            return HomologyVectorSpaceWithBasis(self.base_ring(),
                                                self.complex(),
                                                not self._cohomology)
        return CohomologyRing(self.base_ring(), self.complex())

    def _test_duality(self, **options):
        r"""
        Test if the ordered bases for homology and cohomology are compatible.
        Return nothing if the test succeeds.

        This checks whether each evaluation map `H^n \otimes H_n \to
        k` is represented by the identity matrix, in terms of the
        chosen bases.

        TESTS::

            sage: T = simplicial_complexes.Torus()
            sage: K = T.suspension()
            sage: K.set_immutable()
            sage: H = K.cohomology_ring(QQ)
            sage: H._test_duality()

            sage: simplicial_complexes.RandomComplex(8, 2, .2).homology_with_basis(GF(2))._test_duality()
            sage: simplicial_complexes.RandomComplex(8, 2, .4).homology_with_basis(GF(2))._test_duality()
            sage: simplicial_complexes.RandomComplex(8, 2, .6).homology_with_basis(GF(2))._test_duality()

            sage: simplicial_complexes.RandomComplex(12, 3, .5).homology_with_basis(GF(2))._test_duality() # long time
        """
        tester = self._tester(**options)
        dual = self.dual()
        dims = [a[0] for a in self._indices]
        for dim in range(max(max(dims), tester._max_runs) + 1):
            n = len(self.basis(dim))
            m = matrix(n, n, [a.eval(b) for a in self.basis(dim) for b in dual.basis(dim)])
            tester.assertEqual(m, 1, f"error in dimension {dim}")

    class Element(CombinatorialFreeModule.Element):
        def to_cycle(self):
            r"""
            (Co)cycle representative of this homogeneous (co)homology class.

            EXAMPLES::

                sage: S2 = simplicial_complexes.Sphere(2)
                sage: H = S2.homology_with_basis(QQ)
                sage: h20 = H.basis()[2,0]; h20
                h_{2,0}
                sage: h20.to_cycle()
                -(0, 1, 2) + (0, 1, 3) - (0, 2, 3) + (1, 2, 3)

            Chains are written as linear combinations of simplices
            `\sigma`. Cochains are written as linear combinations of
            characteristic functions `\chi_{\sigma}` for those
            simplices::

                sage: S2.cohomology_ring(QQ).basis()[2,0].to_cycle()
                \chi_(1, 2, 3)
                sage: S2.cohomology_ring(QQ).basis()[0,0].to_cycle()
                \chi_(0,) + \chi_(1,) + \chi_(2,) + \chi_(3,)
            """
            if not self.is_homogeneous():
                raise ValueError("only defined for homogeneous elements")
            return sum(c * self.parent()._to_cycle_on_basis(i) for i, c in self)

        def eval(self, other):
            r"""
            Evaluate ``self`` at ``other``.

            INPUT:

            - ``other`` -- an element of the dual space; if ``self``
              is an element of cohomology in dimension `n`, then
              ``other`` should be an element of homology in dimension
              `n`, and vice versa

            This just calls the :meth:`~sage.homology.chains.Cochains.Element.eval`
            method on the representing chains and cochains.

            EXAMPLES::

                sage: T = simplicial_complexes.Torus()
                sage: homology = T.homology_with_basis(QQ)
                sage: cohomology = T.cohomology_ring(QQ)
                sage: a1, a2 = homology.basis(1)
                sage: alpha1, alpha2 = cohomology.basis(1)
                sage: a1.to_cycle()
                (0, 3) - (0, 6) + (3, 6)
                sage: alpha1.to_cycle()
                -\chi_(1, 3) - \chi_(1, 4) - \chi_(2, 3) - \chi_(2, 4) - \chi_(2, 5) + \chi_(3, 6)
                sage: a1.eval(alpha1)
                1
                sage: alpha2.to_cycle()
                \chi_(1, 3) + \chi_(1, 4) + \chi_(1, 6) + \chi_(2, 4) - \chi_(4, 5) + \chi_(5, 6)
                sage: alpha2.eval(a1)
                0
                sage: (2 * alpha2).eval(a1 + a2)
                2
            """
            if not self or not other:
                return self.base_ring().zero()
            if self.parent()._cohomology:
                return self.to_cycle().eval(other.to_cycle())
            else:
                return other.to_cycle().eval(self.to_cycle())


class HomologyVectorSpaceWithBasis_mod2(HomologyVectorSpaceWithBasis):
    r"""
    Homology vector space mod 2.

    Based on :class:`HomologyVectorSpaceWithBasis`, with Steenrod
    operations included.

    .. NOTE::

        This is not intended to be created directly by the user, but
        instead via the method
        :meth:`~sage.topology.cell_complex.GenericCellComplex.homology_with_basis`
        for the class of :class:`cell
        complexes<sage.topology.cell_complex.GenericCellComplex>`.

    .. TODO::

        Implement Steenrod operations on (co)homology at odd primes,
        and thereby implement this class over `\GF{p}` for any `p`.

    INPUT:

    - ``base_ring`` -- must be the field ``GF(2)``
    - ``cell_complex`` -- the cell complex whose homology we are
      computing
    - ``category`` -- (optional) a subcategory of modules with basis

    This does not include the ``cohomology`` argument present for
    :class:`HomologyVectorSpaceWithBasis`: use
    :class:`CohomologyRing_mod2` for cohomology.

    EXAMPLES:

    Mod 2 cohomology operations are defined on both the left and the
    right::

        sage: # needs sage.groups
        sage: RP4 = simplicial_sets.RealProjectiveSpace(5)
        sage: H = RP4.homology_with_basis(GF(2))
        sage: x4 = H.basis()[4,0]
        sage: x4 * Sq(1)
        h_{3,0}
        sage: Sq(1) * x4
        h_{3,0}
        sage: Sq(2) * x4
        h_{2,0}
        sage: Sq(3) * x4
        h_{1,0}
        sage: Sq(0,1) * x4
        h_{1,0}
        sage: x4 * Sq(0,1)
        h_{1,0}
        sage: Sq(3) * x4
        h_{1,0}
        sage: x4 * Sq(3)
        0
    """
    def __init__(self, base_ring, cell_complex, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: H = simplicial_complexes.Torus().homology_with_basis(GF(2))
            sage: TestSuite(H).run()
            sage: H = simplicial_complexes.Sphere(3).homology_with_basis(GF(2))
            sage: TestSuite(H).run()
        """
        if not is_GF2(base_ring):
            raise ValueError
        category = Modules(base_ring).WithBasis().Graded().FiniteDimensional().or_subcategory(category)
        category = Category.join((category,
                                  LeftModules(SteenrodAlgebra(2)),
                                  RightModules(SteenrodAlgebra(2))))
        HomologyVectorSpaceWithBasis.__init__(self, base_ring,
                                              cell_complex,
                                              cohomology=False,
                                              category=category)

    class Element(HomologyVectorSpaceWithBasis.Element):

        def _acted_upon_(self, a, self_on_left):
            r"""
            Define multiplication of ``self`` by ``a``, an
            element of the Steenrod algebra.

            INPUT:

            - ``a`` -- an element of the mod 2 Steenrod algebra
            - ``self_on_left`` -- ``True`` if we are computing ``self * a``,
              otherwise ``a * self``

            Algorithm: use the action of the Steenrod algebra `A` on
            cohomology to construct the action on homology. That is,
            given a right action of `A` on `H^*`,

            .. MATH::

                \phi_L: H^* \otimes A \to H^*

            we define (a la Boardman [Boa1982]_, p. 190)

            .. MATH::

                S'' \phi_L: A \otimes H_* \to H_*

            using the formula

            .. MATH::

                \langle (S'' \phi) (f \otimes a), x \rangle
                = \langle f, \phi_L (a \otimes x) \rangle,

            for `f \in H_m`, `a \in A^n`, and `x \in
            H^{m-n}`. Somewhat more succinctly, we define the action `f
            \cdot a` by

            .. MATH::

                (f \cdot a) (x) = f (a \cdot x).

            So given `f` (a.k.a. ``self``) and `a`, we compute `f (a
            \cdot x)` for all basis elements `x` in `H^{m-n}`,
            yielding a vector indexed by those basis elements. Since
            our basis for homology is dual to the basis for
            cohomology, we can then use the homology basis to convert
            the vector to an element of `H_{m-n}`.

            This gives a right module structure. To get a left module
            structure, use the right module structure after applying
            the antipode to `a`.

            EXAMPLES::

                sage: # needs sage.groups
                sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
                sage: H = RP5.homology_with_basis(GF(2))
                sage: x5 = list(H.basis(5))[0]
                sage: Sq(1) * x5
                0
                sage: Sq(2) * x5
                h_{3,0}
                sage: x5 * Sq(2)
                h_{3,0}

            TESTS::

                sage: # needs sage.groups
                sage: RP4 = simplicial_sets.RealProjectiveSpace(5)
                sage: H = RP4.homology_with_basis(GF(2))
                sage: x4 = H.basis()[4,0]
                sage: (Sq(1) * Sq(2)) * x4 != 0
                True
                sage: (Sq(1) * Sq(2)) * x4 == Sq(1) * (Sq(2) * x4)
                True
                sage: x4 * (Sq(2) * Sq(1)) == (x4 * Sq(2)) * Sq(1)
                True
                sage: 1 * x4
                h_{4,0}
                sage: x4 * 0
                0
            """
            # Handle field elements first.
            ret = CombinatorialFreeModule.Element._acted_upon_(self, a, self_on_left)
            if ret is not None:  # did the scalar action
                return ret
            m = self.degree()
            n = a.degree()
            if m <= n:
                return self.parent().zero()

            if not self_on_left: # i.e., module element on left
                a = a.antipode()
            P = self.parent()
            B = list(P.basis(m-n))
            return P._from_dict({x.support()[0]: self.eval(a * x)
                                 for x in sorted(self.parent().dual().basis(m-n))})


class CohomologyRing(HomologyVectorSpaceWithBasis):
    """
    The cohomology ring.

    .. NOTE::

        This is not intended to be created directly by the user, but
        instead via the
        :meth:`cohomology ring<sage.topology.cell_complex.GenericCellComplex.cohomology_ring>`
        of a :class:`cell
        complex<sage.topology.cell_complex.GenericCellComplex>`.

    INPUT:

    - ``base_ring`` -- must be a field
    - ``cell_complex`` -- the cell complex whose homology we are
      computing
    - ``category`` -- (optional) a subcategory of modules with basis

    EXAMPLES::

        sage: CP2 = simplicial_complexes.ComplexProjectivePlane()
        sage: H = CP2.cohomology_ring(QQ)
        sage: H.basis(2)
        Finite family {(2, 0): h^{2,0}}
        sage: x = H.basis(2)[2,0]

    The product structure is the cup product::

        sage: x.cup_product(x)
        -h^{4,0}
        sage: x * x
        -h^{4,0}
    """
    def __init__(self, base_ring, cell_complex, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RP2 = simplicial_complexes.ProjectivePlane()
            sage: H = RP2.cohomology_ring(GF(5))
            sage: TestSuite(H).run()
            sage: T = simplicial_complexes.Torus()
            sage: H = T.cohomology_ring(QQ)
            sage: TestSuite(H).run()
        """
        if category is None:
            category = Algebras(base_ring).WithBasis().Graded().FiniteDimensional()
        HomologyVectorSpaceWithBasis.__init__(self, base_ring, cell_complex, True, category)

    def _repr_(self):
        """
        EXAMPLES::

            sage: simplicial_complexes.Torus().cohomology_ring(QQ)
            Cohomology ring of Minimal triangulation of the torus
             over Rational Field
        """
        return "Cohomology ring of {} over {}".format(self._complex, self.base_ring())

    @cached_method
    def one(self):
        """
        The multiplicative identity element.

        EXAMPLES::

            sage: H = simplicial_complexes.Torus().cohomology_ring(QQ)
            sage: H.one()
            h^{0,0}
            sage: all(H.one() * x == x == x * H.one() for x in H.basis())
            True
        """
        one = self.base_ring().one()
        d = {(0, i): one for i in self._graded_indices[0]}
        return self._from_dict(d, remove_zeros=False)

    @cached_method
    def product_on_basis(self, li, ri):
        r"""
        The cup product of the basis elements indexed by ``li`` and ``ri``
        in this cohomology ring.

        INPUT:

        - ``li``, ``ri`` -- index of a cohomology class

        .. SEEALSO::

            :meth:`CohomologyRing.Element.cup_product` -- the
            documentation for this method describes the algorithm.

        EXAMPLES::

            sage: RP3 = simplicial_complexes.RealProjectiveSpace(3)
            sage: H = RP3.cohomology_ring(GF(2))
            sage: c = H.basis()[1,0]
            sage: c.cup_product(c).cup_product(c) # indirect doctest
            h^{3,0}

            sage: T = simplicial_complexes.Torus()
            sage: x,y = T.cohomology_ring(QQ).basis(1)
            sage: x.cup_product(y)
            -h^{2,0}
            sage: x.cup_product(x)
            0

            sage: one = T.cohomology_ring(QQ).basis()[0,0]
            sage: x.cup_product(one)
            h^{1,0}
            sage: one.cup_product(y) == y
            True
            sage: one.cup_product(one)
            h^{0,0}
            sage: x.cup_product(y) + y.cup_product(x)
            0

        This also works with cubical complexes::

            sage: T = cubical_complexes.Torus()
            sage: x,y = T.cohomology_ring(QQ).basis(1)
            sage: x.cup_product(y)
            h^{2,0}
            sage: x.cup_product(x)
            0

        `\Delta`-complexes::

            sage: T_d = delta_complexes.Torus()
            sage: a,b = T_d.cohomology_ring(QQ).basis(1)
            sage: a.cup_product(b)
            h^{2,0}
            sage: b.cup_product(a)
            -h^{2,0}
            sage: RP2 = delta_complexes.RealProjectivePlane()
            sage: w = RP2.cohomology_ring(GF(2)).basis()[1,0]
            sage: w.cup_product(w)
            h^{2,0}

        and simplicial sets::

            sage: from sage.topology.simplicial_set_examples import RealProjectiveSpace
            sage: RP5 = RealProjectiveSpace(5)                                          # needs sage.groups
            sage: x = RP5.cohomology_ring(GF(2)).basis()[1,0]                           # needs sage.groups
            sage: x**4                                                                  # needs sage.groups
            h^{4,0}

        A non-connected example::

            sage: K = cubical_complexes.Torus().disjoint_union(cubical_complexes.Torus())
            sage: a,b,c,d = K.cohomology_ring(QQ).basis(1)
            sage: x,y = K.cohomology_ring(QQ).basis(0)
            sage: a.cup_product(x) == a
            True
            sage: a.cup_product(y)
            0
        """
        scomplex = self.complex()
        base_ring = self.base_ring()
        deg_left = li[0]
        deg_right = ri[0]
        deg_tot = deg_left + deg_right
        left_cycle = self._to_cycle_on_basis(li)
        right_cycle = self._to_cycle_on_basis(ri)
        n_chains_left = scomplex.n_chains(deg_left, base_ring)
        n_chains_right = scomplex.n_chains(deg_right, base_ring)

        result = {}
        H = scomplex.homology_with_basis(base_ring)
        for gamma_index in H._graded_indices.get(deg_tot, []):
            gamma_coeff = base_ring.zero()
            for cell, coeff in H._to_cycle_on_basis((deg_tot, gamma_index)):
                for (c, left_cell, right_cell) in scomplex.alexander_whitney(cell, deg_left):
                    if c:
                        left = n_chains_left(left_cell)
                        right = n_chains_right(right_cell)
                        gamma_coeff += c * coeff * left_cycle.eval(left) * right_cycle.eval(right)
            if gamma_coeff != base_ring.zero():
                result[(deg_tot, gamma_index)] = gamma_coeff
        return self._from_dict(result, remove_zeros=False)

    class Element(HomologyVectorSpaceWithBasis.Element):
        def cup_product(self, other):
            r"""
            Return the cup product of this element and ``other``.

            Algorithm: see González-Díaz and Réal [GDR2003]_, p. 88.
            Given two cohomology classes, lift them to cocycle
            representatives via the chain contraction for this
            complex, using
            :meth:`~HomologyVectorSpaceWithBasis.Element.to_cycle`. In
            the sum of their dimensions, look at all of the homology
            classes `\gamma`: lift each of those to a cycle
            representative, apply the Alexander-Whitney diagonal map
            to each cell in the cycle, evaluate the two cocycles on
            these factors, and multiply. The result is the value of
            the cup product cocycle on this homology class. After this
            has been done for all homology classes, since homology and
            cohomology are dual, one can tell which cohomology class
            corresponds to the cup product.

            .. SEEALSO::

                :meth:`CohomologyRing.product_on_basis`

            EXAMPLES::

                sage: RP3 = simplicial_complexes.RealProjectiveSpace(3)
                sage: H = RP3.cohomology_ring(GF(2))
                sage: c = H.basis()[1,0]
                sage: c.cup_product(c)
                h^{2,0}
                sage: c * c * c
                h^{3,0}

            We can also take powers::

                sage: RP2 = simplicial_complexes.RealProjectivePlane()
                sage: a = RP2.cohomology_ring(GF(2)).basis()[1,0]
                sage: a**0
                h^{0,0}
                sage: a**1
                h^{1,0}
                sage: a**2
                h^{2,0}
                sage: a**3
                0

            A non-connected example::

                sage: K = cubical_complexes.Torus().disjoint_union(cubical_complexes.Sphere(2))
                sage: a,b = K.cohomology_ring(QQ).basis(2)
                sage: a**0
                h^{0,0} + h^{0,1}
            """
            return self * other


class CohomologyRing_mod2(CohomologyRing):
    r"""
    The mod 2 cohomology ring.

    Based on :class:`CohomologyRing`, with Steenrod operations included.

    .. NOTE::

        This is not intended to be created directly by the user, but
        instead via the
        :meth:`cohomology ring<sage.topology.cell_complex.GenericCellComplex.cohomology_ring>`
        of a :class:`cell
        complex<sage.topology.cell_complex.GenericCellComplex>`.

    .. TODO::

        Implement Steenrod operations on (co)homology at odd primes,
        and thereby implement this class over `\GF{p}` for any `p`.

    INPUT:

    - ``base_ring`` -- must be the field ``GF(2)``
    - ``cell_complex`` -- the cell complex whose homology we are
      computing

    EXAMPLES:

    Mod 2 cohomology operations are defined on both the left and the
    right::

        sage: CP2 = simplicial_complexes.ComplexProjectivePlane()
        sage: Hmod2 = CP2.cohomology_ring(GF(2))
        sage: y = Hmod2.basis(2)[2,0]
        sage: y.Sq(2)
        h^{4,0}

        sage: # needs sage.groups
        sage: Y = simplicial_sets.RealProjectiveSpace(6).suspension()
        sage: H_Y = Y.cohomology_ring(GF(2))
        sage: b = H_Y.basis()[2,0]
        sage: b.Sq(1)
        h^{3,0}
        sage: b.Sq(2)
        0
        sage: c = H_Y.basis()[4,0]
        sage: c.Sq(1)
        h^{5,0}
        sage: c.Sq(2)
        h^{6,0}
        sage: c.Sq(3)
        h^{7,0}
        sage: c.Sq(4)
        0

    Cohomology can be viewed as a left module over the Steenrod
    algebra, and also as a right module::

        sage: # needs sage.groups
        sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
        sage: H = RP4.cohomology_ring(GF(2))
        sage: x = H.basis()[1,0]
        sage: Sq(0,1) * x
        h^{4,0}
        sage: Sq(3) * x
        0
        sage: x * Sq(3)
        h^{4,0}
    """
    def __init__(self, base_ring, cell_complex):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RP2 = simplicial_complexes.ProjectivePlane()
            sage: H = RP2.cohomology_ring(GF(2))
            sage: TestSuite(H).run()
        """
        if not is_GF2(base_ring):
            raise ValueError("the base ring must be GF(2)")
        category = Algebras(base_ring).WithBasis().Graded().FiniteDimensional()
        category = Category.join((category,
                                  LeftModules(SteenrodAlgebra(2)),
                                  RightModules(SteenrodAlgebra(2))))
        CohomologyRing.__init__(self, base_ring, cell_complex, category=category)

    class Element(CohomologyRing.Element):
        def Sq(self, i):
            r"""
            Return the result of applying `Sq^i` to this element.

            INPUT:

            - ``i`` -- nonnegative integer

            .. WARNING::

                The main implementation is only for simplicial
                complexes and simplicial sets; cubical complexes are
                converted to simplicial complexes first. Note that
                this converted complex may be large and so computations
                may be slow. There is no implementation for
                `\Delta`-complexes.

            This cohomology operation is only defined in
            characteristic 2. Odd primary Steenrod operations are not
            implemented.

            Algorithm: see González-Díaz and Réal [GDR1999]_,
            Corollary 3.2.

            EXAMPLES::

                sage: RP2 = simplicial_complexes.RealProjectiveSpace(2)
                sage: x = RP2.cohomology_ring(GF(2)).basis()[1,0]
                sage: x.Sq(1)
                h^{2,0}

                sage: K = RP2.suspension()
                sage: K.set_immutable()
                sage: y = K.cohomology_ring(GF(2)).basis()[2,0]
                sage: y.Sq(1)
                h^{3,0}

                sage: # long time
                sage: # needs sage.groups
                sage: RP4 = simplicial_complexes.RealProjectiveSpace(4)
                sage: H = RP4.cohomology_ring(GF(2))
                sage: x = H.basis()[1,0]
                sage: y = H.basis()[2,0]
                sage: z = H.basis()[3,0]
                sage: x.Sq(1) == y
                True
                sage: z.Sq(1)
                h^{4,0}

            This calculation is much faster with simplicial sets (on
            one machine, 20 seconds with a simplicial complex, 4 ms
            with a simplicial set). ::

                sage: RP4_ss = simplicial_sets.RealProjectiveSpace(4)                   # needs sage.groups
                sage: z_ss = RP4_ss.cohomology_ring(GF(2)).basis()[3,0]                 # needs sage.groups
                sage: z_ss.Sq(1)                                                        # needs sage.groups
                h^{4,0}

            TESTS::

                sage: RP_cubical = cubical_complexes.RealProjectivePlane()
                sage: x = RP_cubical.cohomology_ring(GF(2)).basis()[1,0]
                sage: x.Sq(1)
                h^{2,0}
                sage: T = delta_complexes.Torus()
                sage: x = T.cohomology_ring(GF(2)).basis()[1,0]
                sage: x.Sq(1)
                Traceback (most recent call last):
                ...
                NotImplementedError: Steenrod squares are not implemented for this type of cell complex
            """
            P = self.parent()
            scomplex = P.complex()
            if isinstance(scomplex, CubicalComplex):
                # Convert cubical complex to simplicial complex, and
                # convert self to basis element in the new complex's
                # cohomology ring.
                scomplex = SimplicialComplex(scomplex, is_mutable=False)
                P = scomplex.cohomology_ring(self.base_ring())
                self = P.sum_of_terms(self.monomial_coefficients().items())
            if not isinstance(scomplex, (SimplicialComplex, SimplicialSet_arbitrary)):
                print(scomplex, isinstance(scomplex, SimplicialComplex))
                raise NotImplementedError('Steenrod squares are not implemented for '
                                          'this type of cell complex')
            scomplex = P.complex()
            base_ring = P.base_ring()
            if not is_GF2(base_ring):
                # This should never happen: the class should only be
                # instantiated in characteristic 2.
                raise ValueError('Steenrod squares are only defined in characteristic 2')
            # We keep the same notation as in [GDR1999].
            # The trivial cases:
            if i == 0:
                # Sq^0 is the identity.
                return self

            # Construct each graded component of ``self``
            ret = P.zero()
            H = scomplex.homology_with_basis(base_ring)
            deg_comp = {}
            for index, coeff in self:
                d = deg_comp.get(index[0], {})
                d[index] = coeff
                deg_comp[index[0]] = d

            # Do the square on each graded component of ``self``.
            for j in deg_comp:
                # Make it into an actual element
                m = j + i
                if not P._graded_indices.get(m, []) or i > j:
                    continue
                elt = P._from_dict(deg_comp[j], remove_zeros=False)
                if i == j:
                    ret += elt.cup_product(elt)
                    continue

                n = j - i
                # Now assemble the indices over which the sums take place.
                # S(n) is defined to be floor((m+1)/2) + floor(n/2).
                S_n = (m+1) // 2 + n // 2
                if n == 0:
                    sums = [[S_n]]
                else:
                    sums = [[i_n] + l for i_n in range(S_n, m+1)
                            for l in sum_indices(n-1, i_n, S_n)]
                # At this point, 'sums' is a list of lists of the form
                # [i_n, i_{n-1}, ..., i_0]. (It is reversed from the
                # obvious order because this is closer to the order in
                # which the face maps will be applied.)  Now we sum over
                # these, according to the formula in [GDR1999], Corollary 3.2.
                result = {}
                cycle = elt.to_cycle()
                n_chains = scomplex.n_chains(j, base_ring)
                for gamma_index in H._graded_indices.get(m, []):
                    gamma_coeff = base_ring.zero()
                    for cell, coeff in H._to_cycle_on_basis((m, gamma_index)):
                        for indices in sums:
                            indices = list(indices)
                            left = cell
                            right = cell
                            # Since we are working with a simplicial complex, 'cell' is a simplex.
                            if not m % 2:
                                left_endpoint = m
                                while indices:
                                    right_endpoint = indices[0] - 1
                                    for k in range(left_endpoint, indices.pop(0), -1):
                                        left = scomplex.face(left, k)
                                    try:
                                        left_endpoint = indices[0] - 1
                                        for k in range(right_endpoint, indices.pop(0), -1):
                                            right = scomplex.face(right, k)
                                    except IndexError:
                                        pass
                                for k in range(right_endpoint, -1, -1):
                                    right = scomplex.face(right, k)
                            else:
                                right_endpoint = m
                                while indices:
                                    left_endpoint = indices[0] - 1
                                    try:
                                        for k in range(right_endpoint, indices.pop(0), -1):
                                            right = scomplex.face(right, k)
                                        right_endpoint = indices[0] - 1
                                    except IndexError:
                                        pass
                                    for k in range(left_endpoint, indices.pop(0), -1):
                                        left = scomplex.face(left, k)
                                for k in range(right_endpoint, -1, -1):
                                    right = scomplex.face(right, k)

                            if ((hasattr(left, 'is_nondegenerate')
                                 and left.is_nondegenerate()
                                 and right.is_nondegenerate())
                                    or not hasattr(left, 'is_nondegenerate')):
                                left = n_chains(left)
                                right = n_chains(right)
                                gamma_coeff += coeff * cycle.eval(left) * cycle.eval(right)
                    if gamma_coeff != base_ring.zero():
                        result[(m, gamma_index)] = gamma_coeff
                ret += P._from_dict(result, remove_zeros=False)
            return ret

        def _acted_upon_(self, a, self_on_left):
            r"""
            Define multiplication of ``self`` by ``a``, an
            element of the Steenrod algebra.

            INPUT:

            - ``a`` -- an element of the mod 2 Steenrod algebra
            - ``self_on_left`` -- ``True`` if we are computing ``self * a``,
              otherwise ``a * self``

            Algorithm: for left multiplication by ``a``, since we have
            :meth:`Sq` to compute multiplication by a single generator
            `Sq^i`, first convert ``a`` to the Serre-Cartan basis ---
            that is, sums of products of the elements `Sq^i` --- and
            then apply :meth:`Sq` repeatedly. Right multiplication by
            ``a`` is the same as left multiplication by the antipode
            applied to ``a``.

            EXAMPLES::

                sage: # needs sage.groups
                sage: SRP4 = simplicial_sets.RealProjectiveSpace(4).suspension()
                sage: H = SRP4.cohomology_ring(GF(2))
                sage: x = H.basis()[2,0]
                sage: Sq(0,1) * x
                h^{5,0}
                sage: Sq(3) * x
                0
                sage: x * Sq(3)
                h^{5,0}

            TESTS::

                sage: # needs sage.groups
                sage: (Sq(2) * Sq(1)) * x == Sq(2) * (Sq(1) * x)
                True
                sage: x * (Sq(1) * Sq(2)) == (x * Sq(1)) * Sq(2)
                True
                sage: x * 1
                h^{2,0}
                sage: 0 * x
                0
            """
            # Handle field elements first.
            ret = CombinatorialFreeModule.Element._acted_upon_(self, a, self_on_left)
            if ret is not None:  # did the scalar action
                return ret
            if self_on_left: # i.e., module element on left
                a = a.antipode()
            b = a.change_basis('adem')
            ans = self.parent().zero()
            mono_dict = b.monomial_coefficients()
            for seq in mono_dict:
                x = self
                for i in reversed(seq):
                    x = x.Sq(i)
                ans += mono_dict[seq] * x
            return ans

    def steenrod_module_map(self, deg_domain, deg_codomain, side='left'):
        r"""
        Return a component of the module structure map `A \otimes
        H \to H`, where `H` is this cohomology ring and `A` is the
        Steenrod algebra.

        INPUT:

        - ``deg_domain`` -- the degree of the domain in the cohomology
          ring

        - ``deg_codomain`` -- the degree of the codomain in the
          cohomology ring

        - ``side`` -- (default: ``'left'``) whether we are computing
          the action as a left module action or a right module

        We will write this with respect to the left action;
        for the right action, just switch all of the the tensors.
        Writing `m` for ``deg_domain`` and `n` for ``deg_codomain``, this
        returns `A^{n-m} \otimes H^{m} \to H^{n}`, one single
        component of the map making `H` into an `A`-module.

        .. WARNING::

           This is only implemented in characteristic two. The main
           implementation is only for simplicial complexes and simplicial
           sets; cubical complexes are converted to simplicial complexes
           first. Note that this converted complex may be large and so
           computations may be slow. There is no implementation for
           `\Delta`-complexes.

        ALGORITHM:

        Use the Milnor basis for the truncated Steenrod
        algebra `A`, and for cohomology, use the basis with which it
        is equipped. For each pair of basis elements `a` and `h`,
        compute the product `a \otimes h`, and use this to assemble a
        matrix defining the action map via multiplication on the
        appropriate side. That is, if ``side`` is ``'left'``, return a
        matrix suitable for multiplication on the left, etc.

        EXAMPLES::

            sage: # needs sage.groups
            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: H = RP4.cohomology_ring(GF(2))
            sage: H.steenrod_module_map(1, 2)
            [1]
            sage: H.steenrod_module_map(1, 3)
            [0]
            sage: H.steenrod_module_map(1, 4, 'left')
            [1 0]
            sage: H.steenrod_module_map(1, 4, 'right')
            [1]
            [1]

        Products of projective spaces::

            sage: RP3 = simplicial_sets.RealProjectiveSpace(3)
            sage: K = RP3.product(RP3)
            sage: H = K.cohomology_ring(GF(2))
            sage: H
            Cohomology ring of RP^3 x RP^3 over Finite Field of size 2

        There is one column for each element `a \otimes b`, where `a`
        is a basis element for the Steenrod algebra and `b` is a basis
        element for the cohomology algebra.  There is one row for each
        basis element of the cohomology algebra. Unfortunately, the
        chosen basis for this truncated polynomial algebra is not the
        monomial basis::

            sage: x1, x2 = H.basis(1)
            sage: x1 * x1
            h^{2,0} + h^{2,1}
            sage: x2 * x2
            h^{2,2}
            sage: x1 * x2
            h^{2,0}

            sage: H.steenrod_module_map(1, 2)
            [1 0]
            [1 0]
            [0 1]
            sage: H.steenrod_module_map(1, 3, 'left')
            [0 0]
            [0 0]
            [0 0]
            [0 0]
            sage: H.steenrod_module_map(1, 3, 'right')
            [0 0 0 0]
            [0 0 0 0]
            sage: H.steenrod_module_map(2, 3)
            [0 0 0]
            [1 1 0]
            [0 0 0]
            [0 0 0]
        """
        side = side.lower()
        if side not in ['right', 'left']:
            raise ValueError('side must be either "left" or "right"')
        base_ring = self.base_ring()
        char = base_ring.characteristic()
        # Use sorted(...) so that there are reproducible orders on the bases.
        A_basis = sorted(SteenrodAlgebra(char).basis(deg_codomain - deg_domain))
        H_basis_dom = sorted(self.basis(deg_domain))
        H_basis_cod = self.basis(deg_codomain).keys()
        H_basis_cod = {elt: idx for idx, elt in enumerate(sorted(H_basis_cod))}
        entries = []
        for a in A_basis:
            for h in H_basis_dom:
                vec = vector(base_ring, len(H_basis_cod))
                if side == 'left':
                    x = a * h
                else:
                    x = h * a
                monos = x.monomial_coefficients()
                for seq in monos:
                    vec[H_basis_cod[seq]] = monos[seq]
                entries.extend(vec)
        # We built the matrix column by column, so now we take the
        # transpose.
        if side == 'left':
            return matrix(base_ring, len(A_basis) * len(H_basis_dom),
                          len(H_basis_cod), entries).transpose()
        return matrix(base_ring, len(A_basis) * len(H_basis_dom),
                          len(H_basis_cod), entries)


def sum_indices(k, i_k_plus_one, S_k_plus_one):
    r"""
    This is a recursive function for computing the indices for the
    nested sums in González-Díaz and Réal [GDR1999]_, Corollary 3.2.

    In the paper, given indices `i_n`, `i_{n-1}`, ..., `i_{k+1}`,
    given `k`, and given `S(k+1)`, the number `S(k)` is defined to be

    .. MATH::

        S(k) = -S(k+1) + floor(k/2) + floor((k+1)/2) + i_{k+1},

    and `i_k` ranges from `S(k)` to `i_{k+1}-1`. There are two special
    cases: if `k=0`, then `i_0 = S(0)`. Also, the initial case of
    `S(k)` is `S(n)`, which is set in the method :meth:`Sq` before
    calling this function. For this function, given `k`, `i_{k+1}`,
    and `S(k+1)`, return a list consisting of the allowable possible
    indices `[i_k, i_{k-1}, ..., i_1, i_0]` given by the above
    formula.

    INPUT:

    - ``k`` -- nonnegative integer
    - ``i_k_plus_one`` -- the positive integer `i_{k+1}`
    - ``S_k_plus_one`` -- the integer `S(k+1)`

    EXAMPLES::

        sage: from sage.homology.homology_vector_space_with_basis import sum_indices
        sage: sum_indices(1, 3, 3)
        [[1, 0], [2, 1]]
        sage: sum_indices(0, 4, 2)
        [[2]]
    """
    S_k = -S_k_plus_one + k//2 + (k+1)//2 + i_k_plus_one
    if k == 0:
        return [[S_k]]
    return [[i_k] + l for i_k in range(S_k, i_k_plus_one)
            for l in sum_indices(k-1, i_k, S_k)]

def is_GF2(R):
    r"""
    Return ``True`` iff ``R`` is isomorphic to the field `\GF{2}`.

    EXAMPLES::

        sage: from sage.homology.homology_vector_space_with_basis import is_GF2
        sage: is_GF2(GF(2))
        True
        sage: is_GF2(GF(2, impl='ntl'))
        True
        sage: is_GF2(GF(3))
        False
    """
    return 2 == R.characteristic() == R.cardinality()
