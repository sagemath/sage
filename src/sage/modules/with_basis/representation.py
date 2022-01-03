"""
Representations Of A Semigroup

AUTHORS:

- Travis Scrimshaw (2015-11-21): Initial version
- Siddharth Singh  (2020-03-21): Signed Representation
"""

##############################################################################
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.misc.abstract_method import abstract_method
from sage.structure.element import Element
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.modules import Modules


class Representation_abstract(CombinatorialFreeModule):
    """
    Abstract base class for representations of semigroups.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``base_ring`` -- a commutative ring
    """
    def __init__(self, semigroup, base_ring, *args, **opts):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: T = G.trivial_representation()
            sage: TestSuite(T).run()
        """
        self._semigroup = semigroup
        self._semigroup_algebra = semigroup.algebra(base_ring)
        CombinatorialFreeModule.__init__(self, base_ring, *args, **opts)

    def semigroup(self):
        """
        Return the semigroup whose representation ``self`` is.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: R = Representation(G, M, on_basis)
            sage: R.semigroup()
            Symmetric group of order 4! as a permutation group
        """
        return self._semigroup

    def semigroup_algebra(self):
        """
        Return the semigroup algebra whose representation ``self`` is.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: R = Representation(G, M, on_basis)
            sage: R.semigroup_algebra()
            Symmetric group algebra of order 4 over Rational Field
        """
        return self._semigroup_algebra

    @abstract_method
    def side(self):
        """
        Return whether ``self`` is a left, right, or two-sided representation.

        OUTPUT:

        - the string ``"left"``, ``"right"``, or ``"twosided"``

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R.side()
            'left'
        """

    def invariant_module(self, S=None, **kwargs):
        r"""
        Return the submodule of ``self`` invariant under the action of ``S``.

        For a semigroup `S` acting on a module `M`, the invariant
        submodule is given by

        .. MATH::

            M^S = \{m \in M : s \cdot m = m \forall s \in S\}.

        INPUT:

        - ``S`` -- a finitely-generated semigroup (default: the semigroup
          this is a representation of)
        - ``action`` -- a function (default: :obj:`operator.mul`)
        - ``side`` -- ``'left'`` or ``'right'`` (default: :meth:`side()`);
          which side of ``self`` the elements of ``S`` acts

        .. NOTE::

            Two sided actions are considered as left actions for the
            invariant module.

        OUTPUT:

        - :class:`~sage.modules.with_basis.invariant.FiniteDimensionalInvariantModule`

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: M = S3.regular_representation()
            sage: I = M.invariant_module()
            sage: [I.lift(b) for b in I.basis()]
            [() + (2,3) + (1,2) + (1,2,3) + (1,3,2) + (1,3)]

        We build the `D_4`-invariant representation inside of the regular
        representation of `S_4`::

            sage: D4 = groups.permutation.Dihedral(4)
            sage: S4 = SymmetricGroup(4)
            sage: R = S4.regular_representation()
            sage: I = R.invariant_module(D4)
            sage: [I.lift(b) for b in I.basis()]
            [() + (2,4) + (1,2)(3,4) + (1,2,3,4) + (1,3) + (1,3)(2,4) + (1,4,3,2) + (1,4)(2,3),
             (3,4) + (2,3,4) + (1,2) + (1,2,4) + (1,3,2) + (1,3,2,4) + (1,4,3) + (1,4,2,3),
             (2,3) + (2,4,3) + (1,2,3) + (1,2,4,3) + (1,3,4,2) + (1,3,4) + (1,4,2) + (1,4)]
        """
        if S is None:
            S = self.semigroup()
        side = kwargs.pop('side', self.side())
        if side == "twosided":
            side = "left"

        return super().invariant_module(S, side=side, **kwargs)

    def twisted_invariant_module(self, chi, G=None, **kwargs):
        r"""
        Create the isotypic component of the action of ``G`` on
        ``self`` with irreducible character given by ``chi``.

        .. SEEALSO::

            - :class:`~sage.modules.with_basis.invariant.FiniteDimensionalTwistedInvariantModule`

        INPUT:

        - ``chi`` -- a list/tuple of character values or an instance
          of :class:`~sage.groups.class_function.ClassFunction_gap`
        - ``G`` -- a finitely-generated semigroup (default: the semigroup
          this is a representation of)

        This also accepts the group to be the first argument to be the group.

        OUTPUT:

        - :class:`~sage.modules.with_basis.invariant.FiniteDimensionalTwistedInvariantModule`

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: R = G.regular_representation(QQ)
            sage: T = R.twisted_invariant_module([2,0,-1])
            sage: T.basis()
            Finite family {0: B[0], 1: B[1], 2: B[2], 3: B[3]}
            sage: [T.lift(b) for b in T.basis()]
            [() - (1,2,3), -(1,2,3) + (1,3,2), (2,3) - (1,2), -(1,2) + (1,3)]

        We check the different inputs work

            sage: R.twisted_invariant_module([2,0,-1], G) is T
            True
            sage: R.twisted_invariant_module(G, [2,0,-1]) is T
            True
        """
        from sage.categories.groups import Groups
        if G is None:
            G = self.semigroup()
        elif chi in Groups():
            G, chi = chi, G
        side = kwargs.pop('side', self.side())
        if side == "twosided":
            side = "left"

        return super().twisted_invariant_module(G, chi, side=side, **kwargs)

class Representation(Representation_abstract):
    """
    Representation of a semigroup.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``module`` -- a module with a basis
    - ``on_basis`` -- function which takes as input ``g``, ``m``, where
      ``g`` is an element of the semigroup and ``m`` is an element of the
      indexing set for the basis, and returns the result of ``g`` acting
      on ``m``
    - ``side`` -- (default: ``"left"``) whether this is a
      ``"left"`` or ``"right"`` representation

    EXAMPLES:

    We construct the sign representation of a symmetric group::

        sage: G = SymmetricGroup(4)
        sage: M = CombinatorialFreeModule(QQ, ['v'])
        sage: from sage.modules.with_basis.representation import Representation
        sage: on_basis = lambda g,m: M.term(m, g.sign())
        sage: R = Representation(G, M, on_basis)
        sage: x = R.an_element(); x
        2*B['v']
        sage: c,s = G.gens()
        sage: c,s
        ((1,2,3,4), (1,2))
        sage: c * x
        -2*B['v']
        sage: s * x
        -2*B['v']
        sage: c * s * x
        2*B['v']
        sage: (c * s) * x
        2*B['v']

    This extends naturally to the corresponding group algebra::

        sage: A = G.algebra(QQ)
        sage: s,c = A.algebra_generators()
        sage: c,s
        ((1,2,3,4), (1,2))
        sage: c * x
        -2*B['v']
        sage: s * x
        -2*B['v']
        sage: c * s * x
        2*B['v']
        sage: (c * s) * x
        2*B['v']
        sage: (c + s) * x
        -4*B['v']

    REFERENCES:

    - :wikipedia:`Group_representation`
    """
    def __init__(self, semigroup, module, on_basis, side="left", **kwargs):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: R = Representation(G, M, on_basis)
            sage: R._test_representation()

            sage: G = CyclicPermutationGroup(3)
            sage: M = algebras.Exterior(QQ, 'x', 3)
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.prod([M.monomial((g(j+1)-1,)) for j in m]) #cyclically permute generators
            sage: from sage.categories.algebras import Algebras
            sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
            sage: r = R.an_element(); r
            1 + 2*x0 + x0*x1 + 3*x1
            sage: r*r
            1 + 4*x0 + 2*x0*x1 + 6*x1
            sage: x0, x1, x2 = M.gens()
            sage: s = R(x0*x1)
            sage: g = G.an_element()
            sage: g*s
            x1*x2
            sage: g*R(x1*x2)
            -x0*x2
            sage: g*r
            1 + 2*x1 + x1*x2 + 3*x2
            sage: g^2*r
            1 + 3*x0 - x0*x2 + 2*x2

            sage: G = SymmetricGroup(4)
            sage: A = SymmetricGroup(4).algebra(QQ)
            sage: from sage.categories.algebras import Algebras
            sage: from sage.modules.with_basis.representation import Representation
            sage: action = lambda g,x: A.monomial(g*x)
            sage: category = Algebras(QQ).WithBasis().FiniteDimensional()
            sage: R = Representation(G, A, action, 'left', category=category)
            sage: r = R.an_element(); r
            () + (2,3,4) + 2*(1,3)(2,4) + 3*(1,4)(2,3)
            sage: r^2
            14*() + 2*(2,3,4) + (2,4,3) + 12*(1,2)(3,4) + 3*(1,2,4) + 2*(1,3,2) + 4*(1,3)(2,4) + 5*(1,4,3) + 6*(1,4)(2,3)
            sage: g = G.an_element(); g
            (2,3,4)
            sage: g*r
            (2,3,4) + (2,4,3) + 2*(1,3,2) + 3*(1,4,3)
        """
        try:
            self.product_on_basis = module.product_on_basis
        except AttributeError:
            pass

        category = kwargs.pop('category', Modules(module.base_ring()).WithBasis())
        
        if side not in ["left", "right"]:
            raise ValueError('side must be "left" or "right"')
        
        self._left_repr = (side == "left")
        self._on_basis = on_basis
        self._module = module
        
        indices = module.basis().keys()
        
        if 'FiniteDimensional' in module.category().axioms():
            category = category.FiniteDimensional()
        
        Representation_abstract.__init__(self, semigroup, module.base_ring(), indices,
                                         category=category, **module.print_options())

    def _test_representation(self, **options):
        """
        Check (on some elements) that ``self`` is a representation of the
        given semigroup.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R._test_representation()

            sage: G = CoxeterGroup(['A',4,1], base_ring=ZZ)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.term(m, (-1)**g.length())
            sage: R = Representation(G, M, on_basis, side="right")
            sage: R._test_representation(max_runs=500)
        """
        from sage.misc.functional import sqrt
        tester = self._tester(**options)
        S = tester.some_elements()
        L = []
        max_len = int(sqrt(tester._max_runs)) + 1
        for i,x in enumerate(self._semigroup):
            L.append(x)
            if i >= max_len:
                break
        for x in L:
            for y in L:
                for elt in S:
                    if self._left_repr:
                        tester.assertEqual(x*(y*elt), (x*y)*elt)
                    else:
                        tester.assertEqual((elt*y)*x, elt*(y*x))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = Permutations(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: Representation(P, M, on_basis)
            Representation of Standard permutations of 4 indexed by {'v'}
             over Rational Field
        """
        return "Representation of {} indexed by {} over {}".format(
            self._semigroup, self.basis().keys(), self.base_ring())

    def _repr_term(self, b):
        """
        Return a string representation of a basis index ``b`` of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: R = SGA.regular_representation()
            sage: all(R._repr_term(b) == SGA._repr_term(b) for b in SGA.basis().keys())
            True
        """
        return self._module._repr_term(b)

    def _latex_term(self, b):
        """
        Return a LaTeX representation of a basis index ``b`` of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: R = SGA.regular_representation()
            sage: all(R._latex_term(b) == SGA._latex_term(b) for b in SGA.basis().keys())
            True
        """
        return self._module._latex_term(b)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: A = G.algebra(ZZ)
            sage: R = A.regular_representation()
            sage: x = A.an_element(); x
            () + (1,3) + 2*(1,3)(2,4) + 3*(1,4,3,2)
            sage: R(x)
            () + (1,3) + 2*(1,3)(2,4) + 3*(1,4,3,2)
        """
        if isinstance(x, Element) and x.parent() is self._module:
            return self._from_dict(x.monomial_coefficients(copy=False), remove_zeros=False)
        return super(Representation, self)._element_constructor_(x)

    def product_by_coercion(self, left, right):
        """
        Return the product of ``left`` and ``right`` by passing to ``self._module``
        and then building a new element of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.KleinFour()
            sage: E = algebras.Exterior(QQ,'e',4)
            sage: on_basis = lambda g,m: E.monomial(m) # the trivial representation
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G, E, on_basis)
            sage: r = R.an_element(); r
            1 + 2*e0 + 3*e1 + e1*e2
            sage: g = G.an_element();
            sage: g*r == r
            True
            sage: r*r
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
             'Representation of The Klein 4 group of order 4, as a permutation
             group indexed by Subsets of {0, 1, 2, 3} over Rational Field' and 
             'Representation of The Klein 4 group of order 4, as a permutation 
             group indexed by Subsets of {0, 1, 2, 3} over Rational Field'
            
            sage: from sage.categories.algebras import Algebras
            sage: category = Algebras(QQ).FiniteDimensional().WithBasis()
            sage: T = Representation(G, E, on_basis, category=category)
            sage: t = T.an_element(); t
            1 + 2*e0 + 3*e1 + e1*e2
            sage: g*t == t
            True
            sage: t*t
            1 + 4*e0 + 4*e0*e1*e2 + 6*e1 + 2*e1*e2

        """
        M = self._module

        # Multiply in self._module
        p = M._from_dict(left._monomial_coefficients, False, False) * M._from_dict(right._monomial_coefficients, False, False)

        # Convert from a term in self._module to a term in self
        return self._from_dict(p.monomial_coefficients(copy=False), False, False)

    def side(self):
        """
        Return whether ``self`` is a left or a right representation.

        OUTPUT:

        - the string ``"left"`` or ``"right"``

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R.side()
            'left'
            sage: S = G.regular_representation(side="right")
            sage: S.side()
            'right'
        """
        return "left" if self._left_repr else "right"


    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: G = groups.misc.WeylGroup(['B',2], prefix='s')
                sage: R = G.regular_representation()
                sage: s1,s2 = G.gens()
                sage: x = R.an_element(); x
                2*s2*s1*s2 + s1*s2 + 3*s2 + 1
                sage: 2 * x
                4*s2*s1*s2 + 2*s1*s2 + 6*s2 + 2
                sage: s1 * x
                2*s2*s1*s2*s1 + 3*s1*s2 + s1 + s2
                sage: s2 * x
                s2*s1*s2 + 2*s1*s2 + s2 + 3

                sage: G = groups.misc.WeylGroup(['B',2], prefix='s')
                sage: R = G.regular_representation(side="right")
                sage: s1,s2 = G.gens()
                sage: x = R.an_element(); x
                2*s2*s1*s2 + s1*s2 + 3*s2 + 1
                sage: x * s1
                2*s2*s1*s2*s1 + s1*s2*s1 + 3*s2*s1 + s1
                sage: x * s2
                2*s2*s1 + s1 + s2 + 3

                sage: G = groups.misc.WeylGroup(['B',2], prefix='s')
                sage: R = G.regular_representation()
                sage: R.base_ring()
                Integer Ring
                sage: A = G.algebra(ZZ)
                sage: s1,s2 = A.algebra_generators()
                sage: x = R.an_element(); x
                2*s2*s1*s2 + s1*s2 + 3*s2 + 1
                sage: s1 * x
                2*s2*s1*s2*s1 + 3*s1*s2 + s1 + s2
                sage: s2 * x
                s2*s1*s2 + 2*s1*s2 + s2 + 3
                sage: (2*s1 - s2) * x
                4*s2*s1*s2*s1 - s2*s1*s2 + 4*s1*s2 + 2*s1 + s2 - 3
                sage: (3*s1 + s2) * R.zero()
                0

                sage: A = G.algebra(QQ)
                sage: s1,s2 = A.algebra_generators()
                sage: a = 1/2 * s1
                sage: a * x
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *:
                 'Algebra of Weyl Group of type ['B', 2] ... over Rational Field'
                 and 'Left Regular Representation of Weyl Group of type ['B', 2] ... over Integer Ring'

            Check that things that coerce into the group (algebra) also have
            an action::

                sage: D4 = groups.permutation.Dihedral(4)
                sage: S4 = SymmetricGroup(4)
                sage: S4.has_coerce_map_from(D4)
                True
                sage: R = S4.regular_representation()
                sage: D4.an_element() * R.an_element()
                2*(2,4) + 3*(1,2,3,4) + (1,3) + (1,4,2,3)
            """
            if isinstance(scalar, Element):
                P = self.parent()
                sP = scalar.parent()
                if sP is P._semigroup:
                    if not self:
                        return self
                    if self_on_left == P._left_repr:
                        scalar = ~scalar
                    return P.linear_combination(((P._on_basis(scalar, m), c)
                                                 for m,c in self), not self_on_left)

                if sP is P._semigroup_algebra:
                    if not self:
                        return self
                    ret = P.zero()
                    for ms,cs in scalar:
                        if self_on_left == P._left_repr:
                            ms = ~ms
                        ret += P.linear_combination(((P._on_basis(ms, m), cs*c)
                                                    for m,c in self), not self_on_left)
                    return ret

                if P._semigroup.has_coerce_map_from(sP):
                    scalar = P._semigroup(scalar)
                    return self._acted_upon_(scalar, self_on_left)

                # Check for scalars first before general coercion to the semigroup algebra.
                # This will result in a faster action for the scalars.
                ret = CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)
                if ret is not None:
                    return ret

                if P._semigroup_algebra.has_coerce_map_from(sP):
                    scalar = P._semigroup_algebra(scalar)
                    return self._acted_upon_(scalar, self_on_left)

                return None

            return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)

class RegularRepresentation(Representation):
    r"""
    The regular representation of a semigroup.

    The left regular representation of a semigroup `S` over a commutative
    ring `R` is the semigroup ring `R[S]` equipped with the left
    `S`-action `x b_y = b_{xy}`, where `(b_z)_{z \in S}` is the natural
    basis of `R[S]` and `x,y \in S`.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``base_ring`` -- the base ring for the representation
    - ``side`` -- (default: ``"left"``) whether this is a
      ``"left"`` or ``"right"`` representation

    REFERENCES:

    - :wikipedia:`Regular_representation`
    """
    def __init__(self, semigroup, base_ring, side="left"):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: TestSuite(R).run()
        """
        if side == "left":
            on_basis = self._left_on_basis
        else:
            on_basis = self._right_on_basis
        module = semigroup.algebra(base_ring)
        Representation.__init__(self, semigroup, module, on_basis, side)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: G.regular_representation()
            Left Regular Representation of Dihedral group of order 8
             as a permutation group over Integer Ring
            sage: G.regular_representation(side="right")
            Right Regular Representation of Dihedral group of order 8
             as a permutation group over Integer Ring
        """
        if self._left_repr:
            base = "Left Regular Representation"
        else:
            base = "Right Regular Representation"
        return base + " of {} over {}".format(self._semigroup, self.base_ring())

    def _left_on_basis(self, g, m):
        """
        Return the left action of ``g`` on ``m``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R._test_representation()  # indirect doctest
        """
        return self.monomial(g*m)

    def _right_on_basis(self, g, m):
        """
        Return the right action of ``g`` on ``m``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation(side="right")
            sage: R._test_representation()  # indirect doctest
        """
        return self.monomial(m*g)

class TrivialRepresentation(Representation_abstract):
    """
    The trivial representation of a semigroup.

    The trivial representation of a semigroup `S` over a commutative ring
    `R` is the `1`-dimensional `R`-module on which every element of `S`
    acts by the identity.

    This is simultaneously a left and right representation.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``base_ring`` -- the base ring for the representation

    REFERENCES:

    - :wikipedia:`Trivial_representation`
    """
    def __init__(self, semigroup, base_ring):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.PGL(2, 3)
            sage: V = G.trivial_representation()
            sage: TestSuite(V).run()
        """
        cat = Modules(base_ring).WithBasis().FiniteDimensional()
        Representation_abstract.__init__(self, semigroup, base_ring, ['v'], category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: G.trivial_representation()
            Trivial representation of Dihedral group of order 8
             as a permutation group over Integer Ring
        """
        return "Trivial representation of {} over {}".format(self._semigroup,
                                                             self.base_ring())

    def side(self):
        """
        Return that ``self`` is a two-sided representation.

        OUTPUT:

        - the string ``"twosided"``

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.trivial_representation()
            sage: R.side()
            'twosided'
        """
        return "twosided"

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(QQ, 3)
                sage: V = SGA.trivial_representation()
                sage: x = V.an_element()
                sage: 2 * x
                4*B['v']
                sage: all(x * b == x for b in SGA.basis())
                True
                sage: all(b * x == x for b in SGA.basis())
                True
                sage: z = V.zero()
                sage: all(b * z == z for b in SGA.basis())
                True

                sage: H = groups.permutation.Dihedral(5)
                sage: G = SymmetricGroup(5)
                sage: G.has_coerce_map_from(H)
                True
                sage: R = G.trivial_representation(QQ)
                sage: H.an_element() * R.an_element()
                2*B['v']

                sage: AG = G.algebra(QQ)
                sage: AG.an_element() * R.an_element()
                14*B['v']

                sage: AH = H.algebra(ZZ)
                sage: AG.has_coerce_map_from(AH)
                True
                sage: AH.an_element() * R.an_element()
                14*B['v']
            """
            if isinstance(scalar, Element):
                P = self.parent()
                if P._semigroup.has_coerce_map_from(scalar.parent()):
                    return self
                if P._semigroup_algebra.has_coerce_map_from(scalar.parent()):
                    if not self:
                        return self
                    scalar = P._semigroup_algebra(scalar)
                    d = self.monomial_coefficients(copy=True)
                    d['v'] *= sum(scalar.coefficients())
                    return P._from_dict(d)
            return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)


class SignRepresentation_abstract(Representation_abstract):
    """
    Generic implementation of a sign representation.

    The sign representation of a semigroup `S` over a commutative ring
    `R` is the `1`-dimensional `R`-module on which every element of `S`
    acts by 1 if order of element is even (including 0) or -1 if order of element if odd.

    This is simultaneously a left and right representation.

    INPUT:

    - ``permgroup`` -- a permgroup
    - ``base_ring`` -- the base ring for the representation
    - ``sign_function`` -- a function which returns 1 or -1 depending on the elements sign

    REFERENCES:

    - :wikipedia:`Representation_theory_of_the_symmetric_group`
    """

    def __init__(self, group, base_ring, sign_function=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.PGL(2, 3)
            sage: V = G.sign_representation()
            sage: TestSuite(V).run()
        """
        self.sign_function = sign_function
        if sign_function is None:
            try:
                self.sign_function = self._default_sign
            except AttributeError:
                raise TypeError("a sign function must be given")

        cat = Modules(base_ring).WithBasis().FiniteDimensional()

        Representation_abstract.__init__(self, group, base_ring, ["v"], category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: G.sign_representation()
            Sign representation of Dihedral group of order 8
             as a permutation group over Integer Ring
        """
        return "Sign representation of {} over {}".format(
            self._semigroup, self.base_ring()
        )

    def side(self):
        """
        Return that ``self`` is a two-sided representation.

        OUTPUT:

        - the string ``"twosided"``

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.sign_representation()
            sage: R.side()
            'twosided'
        """
        return "twosided"

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: G = PermutationGroup(gens=[(1,2,3), (1,2)])
                sage: S = G.sign_representation()
                sage: x = S.an_element(); x
                2*B['v']
                sage: s,c = G.gens(); c
                (1,2,3)
                sage: s*x
                -2*B['v']
                sage: s*x*s
                2*B['v']
                sage: s*x*s*s*c
                -2*B['v']
                sage: A = G.algebra(ZZ)
                sage: s,c = A.algebra_generators()
                sage: c
                (1,2,3)
                sage: s
                (1,2)
                sage: c*x
                2*B['v']
                sage: c*c*x
                2*B['v']
                sage: c*x*s
                -2*B['v']
                sage: c*x*s*s
                2*B['v']
                sage: (c+s)*x
                0
                sage: (c-s)*x
                4*B['v']

                sage: H = groups.permutation.Dihedral(4)
                sage: G = SymmetricGroup(4)
                sage: G.has_coerce_map_from(H)
                True
                sage: R = G.sign_representation()
                sage: H.an_element() * R.an_element()
                -2*B['v']

                sage: AG = G.algebra(ZZ)
                sage: AH = H.algebra(ZZ)
                sage: AG.has_coerce_map_from(AH)
                True
                sage: AH.an_element() * R.an_element()
                -2*B['v']
            """
            if isinstance(scalar, Element):
                P = self.parent()
                if P._semigroup.has_coerce_map_from(scalar.parent()):
                    scalar = P._semigroup(scalar)
                    return self if P.sign_function(scalar) > 0 else -self

                # We need to check for scalars first
                ret = CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)
                if ret is not None:
                    return ret

                if P._semigroup_algebra.has_coerce_map_from(scalar.parent()):
                    if not self:
                        return self
                    sum_scalar_coeff = 0
                    scalar = P._semigroup_algebra(scalar)
                    for ms, cs in scalar:
                        sum_scalar_coeff += P.sign_function(ms) * cs
                    return sum_scalar_coeff * self

                return None

            return CombinatorialFreeModule.Element._acted_upon_(
                self, scalar, self_on_left
            )


class SignRepresentationPermgroup(SignRepresentation_abstract):
    """
    The sign representation for a permutation group.

    EXAMPLES::

        sage: G = groups.permutation.PGL(2, 3)
        sage: V = G.sign_representation()
        sage: TestSuite(V).run()
    """
    def _default_sign(self, elem):
        """
        Return the sign of the element
        
        INPUT:
      
        - ``elem`` -- the element of the group

        EXAMPLES::
        
            sage: G = groups.permutation.PGL(2, 3)
            sage: V = G.sign_representation()
            sage: elem = G.an_element()
            sage: elem
            (1,2,4,3)
            sage: V._default_sign(elem)
            -1
        """
        return elem.sign()


class SignRepresentationMatrixGroup(SignRepresentation_abstract):
    """
    The sign representation for a matrix group.

    EXAMPLES::

        sage: G = groups.permutation.PGL(2, 3)
        sage: V = G.sign_representation()
        sage: TestSuite(V).run()
    """
    def _default_sign(self, elem):
        """
        Return the sign of the element

        INPUT:

        - ``elem`` -- the element of the group

        EXAMPLES::

            sage: G = GL(2, QQ)
            sage: V = G.sign_representation()
            sage: m = G.an_element()
            sage: m
            [1 0]
            [0 1]
            sage: V._default_sign(m)
            1
        """
        return 1 if elem.matrix().det() > 0 else -1


class SignRepresentationCoxeterGroup(SignRepresentation_abstract):
    """
    The sign representation for a Coxeter group.

    EXAMPLES::

        sage: G = WeylGroup(["A", 1, 1])
        sage: V = G.sign_representation()
        sage: TestSuite(V).run()
    """
    def _default_sign(self, elem):
        """
        Return the sign of the element

        INPUT:

        - ``elem`` -- the element of the group

        EXAMPLES::

            sage: G = WeylGroup(["A", 1, 1])
            sage: elem = G.an_element()
            sage: V = G.sign_representation()
            sage: V._default_sign(elem)
            1
        """
        return -1 if elem.length() % 2 else 1
