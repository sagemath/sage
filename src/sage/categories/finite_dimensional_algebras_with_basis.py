# sage_setup: distribution = sagemath-categories
r"""
Finite dimensional algebras with basis

.. TODO::

    Quotients of polynomial rings.

    Quotients in general.

    Matrix rings.

REFERENCES:

- [CR1962]_
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011-2015 Nicolas M. Thiéry <nthiery at users.sf.net>
#                2011-2015 Franco Saliola <saliola@gmail.com>
#                2014-2015 Aladin Virmaux <aladin.virmaux at u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

import operator
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.algebras import Algebras
from sage.categories.associative_algebras import AssociativeAlgebras
from sage.categories.tensor import TensorProductsCategory


class FiniteDimensionalAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    r"""
    The category of finite dimensional algebras with a distinguished basis.

    EXAMPLES::

        sage: C = FiniteDimensionalAlgebrasWithBasis(QQ); C
        Category of finite dimensional algebras with basis over Rational Field
        sage: C.super_categories()
        [Category of algebras with basis over Rational Field,
         Category of finite dimensional magmatic algebras with basis over Rational Field]
        sage: C.example()                                                               # needs sage.modules
        An example of a finite dimensional algebra with basis:
        the path algebra of the Kronecker quiver
        (containing the arrows a:x->y and b:x->y) over Rational Field

    TESTS::

        sage: # needs sage.graphs sage.modules
        sage: TestSuite(C).run()
        sage: C is Algebras(QQ).FiniteDimensional().WithBasis()
        True
        sage: C is Algebras(QQ).WithBasis().FiniteDimensional()
        True
    """

    class ParentMethods:

        @cached_method
        def radical_basis(self):
            r"""
            Return a basis of the Jacobson radical of this algebra.

            .. NOTE::

               This implementation handles algebras over fields of
               characteristic zero (using Dixon's lemma) or fields of
               characteristic `p` in which we can compute `x^{1/p}`
               [FR1985]_, [Eb1989]_.

            OUTPUT: list of elements of ``self``

            .. SEEALSO:: :meth:`radical`, :class:`Algebras.Semisimple`

            EXAMPLES::

                sage: # needs sage.graphs sage.modules
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: A.radical_basis()
                (a, b)

            We construct the group algebra of the Klein Four-Group
            over the rationals::

                sage: A = KleinFourGroup().algebra(QQ)                                  # needs sage.groups sage.modules

            This algebra belongs to the category of finite dimensional
            algebras over the rationals::

                sage: A in Algebras(QQ).FiniteDimensional().WithBasis()                 # needs sage.groups sage.modules
                True

            Since the field has characteristic `0`, Maschke's Theorem
            tells us that the group algebra is semisimple. So its
            radical is the zero ideal::

                sage: A in Algebras(QQ).Semisimple()                                    # needs sage.groups sage.modules
                True
                sage: A.radical_basis()                                                 # needs sage.groups sage.modules
                ()

            Let's work instead over a field of characteristic `2`::

                sage: A = KleinFourGroup().algebra(GF(2))                               # needs sage.groups sage.modules
                sage: A in Algebras(GF(2)).Semisimple()                                 # needs sage.groups sage.modules
                False
                sage: A.radical_basis()                                                 # needs sage.groups sage.modules
                (() + (1,2)(3,4), (3,4) + (1,2)(3,4), (1,2) + (1,2)(3,4))

            We now implement the algebra `A = K[x] / (x^p-1)`, where `K`
            is a finite field of characteristic `p`, and check its
            radical; alas, we currently need to wrap `A` to make it a
            proper :class:`ModulesWithBasis`::

                sage: # needs sage.modules
                sage: class AnAlgebra(CombinatorialFreeModule):
                ....:     def __init__(self, F):
                ....:         R.<x> = PolynomialRing(F)
                ....:         I = R.ideal(x**F.characteristic()-F.one())
                ....:         self._xbar = R.quotient(I).gen()
                ....:         basis_keys = [self._xbar**i for i in range(F.characteristic())]
                ....:         CombinatorialFreeModule.__init__(self, F, basis_keys,
                ....:                 category=Algebras(F).FiniteDimensional().WithBasis())
                ....:     def one(self):
                ....:         return self.basis()[self.base_ring().one()]
                ....:     def product_on_basis(self, w1, w2):
                ....:         return self.from_vector(vector(w1*w2))
                sage: AnAlgebra(GF(3)).radical_basis()                                  # needs sage.libs.pari
                (B[1] + 2*B[xbar^2], B[xbar] + 2*B[xbar^2])
                sage: AnAlgebra(GF(16,'a')).radical_basis()                             # needs sage.rings.finite_rings
                (B[1] + B[xbar],)
                sage: AnAlgebra(GF(49,'a')).radical_basis()                             # needs sage.rings.finite_rings
                (B[1] + 6*B[xbar^6], B[xbar] + 6*B[xbar^6], B[xbar^2] + 6*B[xbar^6],
                 B[xbar^3] + 6*B[xbar^6], B[xbar^4] + 6*B[xbar^6], B[xbar^5] + 6*B[xbar^6])

            We compute the radical basis in a subalgebra using
            the inherited product::

                sage: # needs sage.modules
                sage: scoeffs = {('a','e'): {'a':1}, ('b','e'): {'a':1, 'b':1},
                ....:            ('c','d'): {'a':1}, ('c','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: MS = MatrixSpace(QQ, 5)
                sage: A = MS.subalgebra([bg.adjoint_matrix() for bg in L.lie_algebra_generators()])
                sage: A.radical_basis()
                (B[1], B[2], B[3], B[4], B[5])

            TESTS::

                sage: # needs sage.groups sage.modules
                sage: A = KleinFourGroup().algebra(GF(2))
                sage: A.radical_basis()
                (() + (1,2)(3,4), (3,4) + (1,2)(3,4), (1,2) + (1,2)(3,4))
                sage: A = KleinFourGroup().algebra(QQ, category=Monoids())
                sage: A.radical_basis.__module__
                'sage.categories.finite_dimensional_algebras_with_basis'
                sage: A.radical_basis()
                ()
            """
            F = self.base_ring()
            if not F.is_field():
                raise NotImplementedError("the base ring must be a field")
            p = F.characteristic()
            from sage.matrix.constructor import matrix
            from sage.modules.free_module_element import vector

            if p == 0:
                B = self.basis()
                product_on_basis = self.product_on_basis
                if product_on_basis is NotImplemented:
                    def product_on_basis(i, j):
                        return B[i] * B[j]

                keys = B.keys()
                cache = [{(i, j): c for i in keys for j, c in product_on_basis(y, i)}
                         for y in keys]
                mat = [[sum(x.get((j, i), 0) * c for (i,j), c in y.items())
                        for x in cache]
                       for y in cache]

                mat = matrix(self.base_ring(), mat)
                rad_basis = mat.kernel().basis()

            else:
                # TODO: some finite field elements in Sage have both an
                # ``nth_root`` method and a ``pth_root`` method (such as ``GF(9,'a')``),
                # some only have a ``nth_root`` element such as ``GF(2)``
                # I imagine that ``pth_root`` would be fastest, but it is not
                # always available....
                if hasattr(self.base_ring().one(), 'nth_root'):
                    def root_fcn(s, x):
                        return x.nth_root(s)

                else:
                    def root_fcn(s, x):
                        return x ** (1 / s)

                s = 1
                n = self.dimension()
                B = [b.on_left_matrix() for b in self.basis()]
                while s <= n:
                    # we use that p_{AB}(x) = p_{BA}(x) here
                    data = [[None]*(len(B)+1) for _ in B]
                    for i, b in enumerate(B):
                        for j, bb in enumerate(B[i:], start=i):
                            val = (-1)**s * (b*bb).charpoly()[n-s]
                            data[i][j] = data[j][i] = val
                        data[i][-1] = (-1)**s * b.charpoly()[n-s]
                    C = matrix(data).left_kernel().basis()
                    if 1 < s < F.order():
                        C = [vector(F, [root_fcn(s, ci) for ci in c]) for c in C]
                    B = [sum(ci * b for (ci, b) in zip(c, B)) for c in C]
                    s = p * s
                e = vector(self.one())
                rad_basis = [b * e for b in B]

            return tuple([self.from_vector(vec) for vec in rad_basis])

        @cached_method
        def radical(self):
            r"""
            Return the Jacobson radical of ``self``.

            This uses :meth:`radical_basis`, whose default
            implementation handles algebras over fields of
            characteristic zero or fields of characteristic `p` in
            which we can compute `x^{1/p}`.

            .. SEEALSO:: :meth:`radical_basis`, :meth:`semisimple_quotient`

            EXAMPLES::

                sage: # needs sage.graphs sage.modules
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: radical = A.radical(); radical
                Radical of An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field

            The radical is an ideal of `A`, and thus a finite
            dimensional non unital associative algebra::

                sage: # needs sage.graphs sage.modules
                sage: from sage.categories.associative_algebras import AssociativeAlgebras
                sage: radical in AssociativeAlgebras(QQ).WithBasis().FiniteDimensional()
                True
                sage: radical in Algebras(QQ)
                False

                sage: # needs sage.graphs sage.modules
                sage: radical.dimension()
                2
                sage: radical.basis()
                Finite family {0: B[0], 1: B[1]}
                sage: radical.ambient() is A
                True
                sage: [c.lift() for c in radical.basis()]
                [a, b]

            .. TODO::

                - Tell Sage that the radical is in fact an ideal;
                - Pickling by construction, as ``A.center()``;
                - Lazy evaluation of ``_repr_``.

            TESTS::

                sage: # needs sage.graphs sage.modules
                sage: TestSuite(radical).run()
            """
            category = AssociativeAlgebras(self.category().base_ring()).WithBasis().FiniteDimensional().Subobjects()
            radical = self.submodule(self.radical_basis(),
                                     category=category,
                                     already_echelonized=True)
            radical.rename("Radical of {}".format(self))
            return radical

        @cached_method
        def semisimple_quotient(self):
            """
            Return the semisimple quotient of ``self``.

            This is the quotient of ``self`` by its radical.

            .. SEEALSO:: :meth:`radical`

            EXAMPLES::

                sage: # needs sage.graphs sage.modules
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: a,b,x,y = sorted(A.basis())
                sage: S = A.semisimple_quotient(); S
                Semisimple quotient of An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: S in Algebras(QQ).Semisimple()
                True
                sage: S.basis()
                Finite family {'x': B['x'], 'y': B['y']}
                sage: xs,ys = sorted(S.basis())
                sage: (xs + ys) * xs
                B['x']

            Sanity check: the semisimple quotient of the `n`-th
            descent algebra of the symmetric group is of dimension the
            number of partitions of `n`::

                sage: [ DescentAlgebra(QQ,n).B().semisimple_quotient().dimension()      # needs sage.combinat sage.groups sage.modules
                ....:   for n in range(6) ]
                [1, 1, 2, 3, 5, 7]
                sage: [Partitions(n).cardinality() for n in range(10)]                  # needs sage.combinat
                [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]

            .. TODO::

               - Pickling by construction, as ``A.semisimple_quotient()``?
               - Lazy evaluation of ``_repr_``

            TESTS::

                sage: TestSuite(S).run()                                                # needs sage.graphs sage.modules
            """
            ring = self.base_ring()
            category = Algebras(ring).WithBasis().FiniteDimensional().Quotients().Semisimple()
            result = self.quotient_module(self.radical(), category=category)
            result.rename("Semisimple quotient of {}".format(self))
            return result

        @cached_method
        def center_basis(self):
            r"""
            Return a basis of the center of ``self``.

            OUTPUT: list of elements of ``self``

            .. SEEALSO:: :meth:`center`

            EXAMPLES::

                sage: # needs sage.graphs sage.modules
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: A.center_basis()
                (x + y,)
            """
            return self.annihilator_basis(self.algebra_generators(), self.bracket)

        @cached_method
        def center(self):
            r"""
            Return the center of ``self``.

            .. SEEALSO:: :meth:`center_basis`

            EXAMPLES::

                sage: # needs sage.graphs sage.modules
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: center = A.center(); center
                Center of An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: center in Algebras(QQ).WithBasis().FiniteDimensional().Commutative()
                True
                sage: center.dimension()
                1
                sage: center.basis()
                Finite family {0: B[0]}
                sage: center.ambient() is A
                True
                sage: [c.lift() for c in center.basis()]
                [x + y]

            The center of a semisimple algebra is semisimple::

                sage: A = DihedralGroup(6).algebra(QQ)                                  # needs sage.groups sage.modules
                sage: A.center() in Algebras(QQ).Semisimple()                           # needs sage.groups sage.modules
                True

            .. TODO::

                - Pickling by construction, as ``A.center()``?
                - Lazy evaluation of ``_repr_``

            TESTS::

                sage: TestSuite(center).run()                                           # needs sage.graphs sage.modules
            """
            category = Algebras(self.base_ring()).FiniteDimensional().Subobjects().Commutative().WithBasis()
            if self in Algebras.Semisimple:
                category = category.Semisimple()
            center = self.submodule(self.center_basis(),
                                    category=category,
                                    already_echelonized=True)
            center.rename("Center of {}".format(self))
            return center

        def subalgebra(self, gens, category=None, *args, **opts):
            r"""
            Return the subalgebra of ``self`` generated by ``gens``.

            EXAMPLES::

                sage: # needs sage.modules
                sage: scoeffs = {('a','e'): {'a':1}, ('b','e'): {'a':1, 'b':1},
                ....:            ('c','d'): {'a':1}, ('c','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: MS = MatrixSpace(QQ, 5)
                sage: A = MS.subalgebra([bg.adjoint_matrix() for bg in L.lie_algebra_generators()])
                sage: A.dimension()
                7

                sage: # needs sage.modules
                sage: L.<x,y,z> = LieAlgebra(GF(3), {('x','z'): {'x':1, 'y':1}, ('y','z'): {'y':1}})
                sage: MS = MatrixSpace(L.base_ring(), L.dimension())
                sage: gens = [b.adjoint_matrix() for b in L.basis()]
                sage: A = MS.subalgebra(gens)
                sage: A.dimension()
                5
            """
            # add the unit to make sure it is unital
            basis = []
            new_elts = [self(g) for g in gens] + [self.one()]
            while new_elts:
                basis = self.echelon_form(basis + new_elts)
                trailsupp = {b.trailing_support(): b for b in basis}
                sortsupp = sorted(trailsupp)
                new_elts = []
                # We (re)implement the reduction here
                for b in basis:
                    for bp in basis:
                        elt = b * bp
                        for s in sortsupp:
                            c = elt[s]
                            if c:
                                elt -= c / trailsupp[s].trailing_coefficient() * trailsupp[s]
                        if elt:
                            new_elts.append(elt)
            C = FiniteDimensionalAlgebrasWithBasis(self.category().base_ring())
            category = C.Subobjects().or_subcategory(category)
            return self.submodule(basis, check=False, already_echelonized=True,
                                  category=category)

        def ideal_submodule(self, gens, side='left', category=None, *args, **opts):
            r"""
            Return the ``side`` ideal of ``self`` generated by ``gens``
            as a submodule.

            .. TODO::

                This is not generally compatible with the implementation of
                the ideals. This method should be folded into the ``ideal``
                method after the corresponding classes are refactored to
                be compatible.

            EXAMPLES::

                sage: # needs sage.modules
                sage: scoeffs = {('a','e'): {'a':1}, ('b','e'): {'a':1, 'b':1},
                ....:            ('c','d'): {'a':1}, ('c','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: MS = MatrixSpace(QQ, 5)
                sage: I = MS.ideal_submodule([bg.adjoint_matrix() for bg in L.lie_algebra_generators()])
                sage: I.dimension()
                25
            """
            C = AssociativeAlgebras(self.category().base_ring()).WithBasis().FiniteDimensional()
            category = C.Subobjects().or_subcategory(category)
            if gens in self:
                gens = [self(gens)]
            if side == 'left':
                return self.submodule([b * self(g) for b in self.basis() for g in gens],
                                      category=category, *args, **opts)
            if side == 'right':
                return self.submodule([self(g) * b for b in self.basis() for g in gens],
                                      category=category, *args, **opts)
            if side == 'twosided':
                return self.submodule([b * self(g) * bp for b in self.basis()
                                       for bp in self.basis() for g in gens],
                                      category=category, *args, **opts)
            raise ValueError("side must be either 'left', 'right', or 'twosided'")

        def principal_ideal(self, a, side='left', *args, **opts):
            r"""
            Construct the ``side`` principal ideal generated by ``a``.

            INPUT:

            - ``a`` -- an element
            - ``side`` -- ``left`` (default) or ``right``
            - ``coerce`` -- ignored, for compatibility with categories

            EXAMPLES:

            In order to highlight the difference between left and
            right principal ideals, our first example deals with a non
            commutative algebra::

                sage: # needs sage.graphs sage.modules
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: x, y, a, b = A.basis()

            In this algebra, multiplication on the right by `x`
            annihilates all basis elements but `x`::

                sage: x*x, y*x, a*x, b*x                                                # needs sage.graphs sage.modules
                (x, 0, 0, 0)

            so the left ideal generated by `x` is one-dimensional::

                sage: Ax = A.principal_ideal(x, side='left'); Ax                        # needs sage.graphs sage.modules
                Free module generated by {0} over Rational Field
                sage: [B.lift() for B in Ax.basis()]                                    # needs sage.graphs sage.modules
                [x]

            Multiplication on the left by `x` annihilates
            only `x` and fixes the other basis elements::

                sage: x*x, x*y, x*a, x*b                                                # needs sage.graphs sage.modules
                (x, 0, a, b)

            so the right ideal generated by `x` is 3-dimensional::

                sage: xA = A.principal_ideal(x, side='right'); xA                       # needs sage.graphs sage.modules
                Free module generated by {0, 1, 2} over Rational Field
                sage: [B.lift() for B in xA.basis()]                                    # needs sage.graphs sage.modules
                [x, a, b]

            .. SEEALSO::

                - :meth:`peirce_summand`
            """
            opts.pop("coerce", None)
            return self.submodule([(a * b if side == 'right' else b * a)
                                   for b in self.basis()], *args, **opts)

        @cached_method
        def orthogonal_idempotents_central_mod_radical(self):
            r"""
            Return a family of orthogonal idempotents of ``self`` that project
            on the central orthogonal idempotents of the semisimple quotient.

            OUTPUT:

            - a list of orthogonal idempotents obtained by lifting the central
              orthogonal idempotents of the semisimple quotient.

            ALGORITHM:

            The orthogonal idempotents of `A` are obtained by lifting the
            central orthogonal idempotents of the semisimple quotient
            `\overline{A}`.

            Namely, let `(\overline{f_i})` be the central orthogonal
            idempotents of the semisimple quotient of `A`. We
            recursively construct orthogonal idempotents of `A` by the
            following procedure: assuming `(f_i)_{i < n}` is a set of
            already constructed orthogonal idempotent, we construct
            `f_k` by idempotent lifting of `(1-f) g (1-f)`, where `g`
            is any lift of `\overline{e_k}` and `f=\sum_{i<k} f_i`.

            See [CR1962]_ for correctness and termination proofs.

            .. SEEALSO::

                - :meth:`Algebras.SemiSimple.FiniteDimensional.WithBasis.ParentMethods.central_orthogonal_idempotents`
                - :meth:`idempotent_lift`

            EXAMPLES::

                sage: # needs sage.graphs sage.modules
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: A.orthogonal_idempotents_central_mod_radical()                    # needs sage.rings.number_field
                (x, y)

            ::

                sage: # needs sage.modules sage.rings.number_field
                sage: Z12 = Monoids().Finite().example(); Z12
                An example of a finite multiplicative monoid: the integers modulo 12
                sage: A = Z12.algebra(QQ)
                sage: idempotents = A.orthogonal_idempotents_central_mod_radical()
                sage: sorted(idempotents, key=str)
                [-B[0] + 1/2*B[4] + 1/2*B[8],
                 1/2*B[4] - 1/2*B[8],
                 1/2*B[9] + 1/2*B[3] - B[0],
                 1/2*B[9] - 1/2*B[3],
                 1/4*B[1] + 1/4*B[11] - 1/4*B[5] - 1/4*B[7],
                 1/4*B[1] - 1/2*B[9] + 1/4*B[5] - 1/4*B[7] + 1/2*B[3] - 1/4*B[11],
                 1/4*B[1] - 1/2*B[9] - 1/2*B[3] + 1/4*B[11] + 1/4*B[5] + 1/4*B[7] + B[0] - 1/2*B[4] - 1/2*B[8],
                 1/4*B[1] - 1/4*B[5] + 1/4*B[7] - 1/4*B[11] - 1/2*B[4] + 1/2*B[8],
                 B[0]]
                sage: sum(idempotents) == 1
                True
                sage: all(e*e == e for e in idempotents)
                True
                sage: all(e*f == 0 and f*e == 0
                ....:     for e in idempotents for f in idempotents if e != f)
                True

            This is best tested with::

                sage: A.is_identity_decomposition_into_orthogonal_idempotents(idempotents)          # needs sage.graphs sage.modules sage.rings.number_field
                True

            We construct orthogonal idempotents for the algebra of the
            `0`-Hecke monoid::

                sage: # needs sage.combinat sage.graphs sage.groups sage.modules
                sage: from sage.monoids.hecke_monoid import HeckeMonoid
                sage: A = HeckeMonoid(SymmetricGroup(4)).algebra(QQ)
                sage: idempotents = A.orthogonal_idempotents_central_mod_radical()
                sage: A.is_identity_decomposition_into_orthogonal_idempotents(idempotents)
                True
            """
            one = self.one()
            # Construction of the orthogonal idempotents
            idempotents = []
            f = self.zero()
            for g in self.semisimple_quotient().central_orthogonal_idempotents():
                fi = self.idempotent_lift((one - f) * g.lift() * (one - f))
                idempotents.append(fi)
                f = f + fi
            return tuple(idempotents)

        def idempotent_lift(self, x):
            r"""
            Lift an idempotent of the semisimple quotient into an idempotent of ``self``.

            Let `A` be this finite dimensional algebra and `\pi` be
            the projection `A \rightarrow \overline{A}` on its
            semisimple quotient. Let `\overline{x}` be an idempotent
            of `\overline A`, and `x` any lift thereof in `A`. This
            returns an idempotent `e` of `A` such that `\pi(e)=\pi(x)`
            and `e` is a polynomial in `x`.

            INPUT:

            - ``x`` -- an element of `A` that projects on an idempotent
              `\overline x` of the semisimple quotient of `A`.
              Alternatively one may give as input the idempotent
              `\overline{x}`, in which case some lift thereof will be
              taken for `x`.

            OUTPUT: the idempotent `e` of ``self``

            ALGORITHM:

            Iterate the formula `1 - (1 - x^2)^2` until having an
            idempotent.

            See [CR1962]_ for correctness and termination proofs.

            EXAMPLES::

                sage: # needs sage.graphs sage.modules
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example()
                sage: S = A.semisimple_quotient()
                sage: A.idempotent_lift(S.basis()['x'])
                x
                sage: A.idempotent_lift(A.basis()['y'])
                y

            .. TODO::

                Add some non trivial example
            """
            if not self.is_parent_of(x):
                x = x.lift()
            p = self.semisimple_quotient().retract(x)
            if p * p != p:
                raise ValueError("%s does not retract to an idempotent." % p)
            x_prev = None
            one = self.one()
            while x != x_prev:
                tmp = x
                x = (one - (one - x**2)**2)
                x_prev = tmp
            return x

        @cached_method
        def cartan_invariants_matrix(self):
            r"""
            Return the Cartan invariants matrix of the algebra.

            OUTPUT: a matrix of nonnegative integers

            Let `A` be this finite dimensional algebra and
            `(S_i)_{i\in I}` be representatives of the right simple
            modules of `A`. Note that their adjoints `S_i^*` are
            representatives of the left simple modules.

            Let `(P^L_i)_{i\in I}` and `(P^R_i)_{i\in I}` be
            respectively representatives of the corresponding
            indecomposable projective left and right modules of `A`.
            In particular, we assume that the indexing is consistent
            so that `S_i^*=\operatorname{top} P^L_i` and
            `S_i=\operatorname{top} P^R_i`.

            The *Cartan invariant matrix* `(C_{i,j})_{i,j\in I}` is a
            matrix of nonnegative integers that encodes much of the
            representation theory of `A`; namely:

            - `C_{i,j}` counts how many times `S_i^*\otimes S_j`
              appears as composition factor of `A` seen as a bimodule
              over itself;

            - `C_{i,j}=\dim Hom_A(P^R_j, P^R_i)`;

            - `C_{i,j}` counts how many times `S_j` appears as
              composition factor of `P^R_i`;

            - `C_{i,j}=\dim Hom_A(P^L_i, P^L_j)`;

            - `C_{i,j}` counts how many times `S_i^*` appears as
              composition factor of `P^L_j`.

            In the commutative case, the Cartan invariant matrix is
            diagonal. In the context of solving systems of
            multivariate polynomial equations of dimension zero, `A`
            is the quotient of the polynomial ring by the ideal
            generated by the equations, the simple modules correspond
            to the roots, and the numbers `C_{i,i}` give the
            multiplicities of those roots.

            .. NOTE::

                For simplicity, the current implementation assumes
                that the index set `I` is of the form
                `\{0,\dots,n-1\}`. Better indexations will be possible
                in the future.

            ALGORITHM:

            The Cartan invariant matrix of `A` is computed from the
            dimension of the summands of its Peirce decomposition.

            .. SEEALSO::

                - :meth:`peirce_decomposition`
                - :meth:`isotypic_projective_modules`

            EXAMPLES:

            For a semisimple algebra, in particular for group algebras
            in characteristic zero, the Cartan invariants matrix is
            the identity::

                sage: A3 = SymmetricGroup(3).algebra(QQ)                                # needs sage.combinat sage.groups sage.modules
                sage: A3.cartan_invariants_matrix()                                     # needs sage.combinat sage.groups sage.modules
                [1 0 0]
                [0 1 0]
                [0 0 1]

            For the path algebra of a quiver, the Cartan invariants
            matrix counts the number of paths between two vertices::

                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example()
                sage: A.cartan_invariants_matrix()                                      # needs sage.modules sage.rings.number_field
                [1 2]
                [0 1]

            In the commutative case, the Cartan invariant matrix is diagonal::

                sage: Z12 = Monoids().Finite().example(); Z12
                An example of a finite multiplicative monoid: the integers modulo 12
                sage: A = Z12.algebra(QQ)                                               # needs sage.modules
                sage: A.cartan_invariants_matrix()                                      # needs sage.modules sage.rings.number_field
                [1 0 0 0 0 0 0 0 0]
                [0 1 0 0 0 0 0 0 0]
                [0 0 2 0 0 0 0 0 0]
                [0 0 0 1 0 0 0 0 0]
                [0 0 0 0 2 0 0 0 0]
                [0 0 0 0 0 1 0 0 0]
                [0 0 0 0 0 0 1 0 0]
                [0 0 0 0 0 0 0 2 0]
                [0 0 0 0 0 0 0 0 1]

            With the algebra of the `0`-Hecke monoid::

                sage: # needs sage.combinat sage.groups sage.modules
                sage: from sage.monoids.hecke_monoid import HeckeMonoid
                sage: A = HeckeMonoid(SymmetricGroup(4)).algebra(QQ)
                sage: A.cartan_invariants_matrix()                                      # needs sage.rings.number_field
                [1 0 0 0 0 0 0 0]
                [0 2 1 0 1 1 0 0]
                [0 1 1 0 1 0 0 0]
                [0 0 0 1 0 1 1 0]
                [0 1 1 0 1 0 0 0]
                [0 1 0 1 0 2 1 0]
                [0 0 0 1 0 1 1 0]
                [0 0 0 0 0 0 0 1]
            """
            from sage.matrix.constructor import Matrix
            from sage.rings.integer_ring import ZZ
            A_quo = self.semisimple_quotient()
            idempotents_quo = A_quo.central_orthogonal_idempotents()
            # Dimension of simple modules
            dim_simples = [A_quo.principal_ideal(e).dimension().sqrt()
                          for e in idempotents_quo]
            # Orthogonal idempotents
            idempotents = self.orthogonal_idempotents_central_mod_radical()

            def C(i, j):
                summand = self.peirce_summand(idempotents[i], idempotents[j])
                return summand.dimension() / (dim_simples[i] * dim_simples[j])
            m = Matrix(ZZ, len(idempotents), C)
            m.set_immutable()
            return m

        def isotypic_projective_modules(self, side='left'):
            r"""
            Return the isotypic projective ``side`` ``self``-modules.

            Let `P_i` be representatives of the indecomposable
            projective ``side``-modules of this finite dimensional
            algebra `A`, and `S_i` be the associated simple modules.

            The regular ``side`` representation of `A` can be
            decomposed as a direct sum `A = \bigoplus_i Q_i` where
            each `Q_i` is an isotypic projective module; namely `Q_i`
            is the direct sum of `\dim S_i` copies of the
            indecomposable projective module `P_i`. This decomposition
            is not unique.

            The isotypic projective modules are constructed as
            `Q_i=e_iA`, where the `(e_i)_i` is the decomposition of
            the identity into orthogonal idempotents obtained by
            lifting the central orthogonal idempotents of the
            semisimple quotient of `A`.

            INPUT:

            - ``side`` -- ``'left'`` or ``'right'`` (default: ``'left'``)

            OUTPUT: list of subspaces of ``self``

            EXAMPLES::

                sage: # needs sage.graphs sage.modules sage.rings.number_field
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: Q = A.isotypic_projective_modules(side='left'); Q
                [Free module generated by {0} over Rational Field,
                 Free module generated by {0, 1, 2} over Rational Field]
                sage: [[x.lift() for x in Qi.basis()]
                ....:  for Qi in Q]
                [[x],
                 [y, a, b]]

            We check that the sum of the dimensions of the isotypic
            projective modules is the dimension of ``self``::

                sage: sum([Qi.dimension() for Qi in Q]) == A.dimension()                # needs sage.graphs sage.modules sage.rings.number_field
                True

            .. SEEALSO::

                - :meth:`orthogonal_idempotents_central_mod_radical`
                - :meth:`peirce_decomposition`
            """
            return [self.principal_ideal(e, side) for e in
                    self.orthogonal_idempotents_central_mod_radical()]

        @cached_method
        def peirce_summand(self, ei, ej):
            r"""
            Return the Peirce decomposition summand `e_i A e_j`.

            INPUT:

            - ``self`` -- an algebra `A`

            - ``ei``, ``ej`` -- two idempotents of `A`

            OUTPUT: `e_i A e_j`, as a subspace of `A`

            .. SEEALSO::

                - :meth:`peirce_decomposition`
                - :meth:`principal_ideal`

            EXAMPLES::

                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example()
                sage: idemp = A.orthogonal_idempotents_central_mod_radical()            # needs sage.rings.number_field
                sage: A.peirce_summand(idemp[0], idemp[1])                              # needs sage.rings.number_field
                Free module generated by {0, 1} over Rational Field
                sage: A.peirce_summand(idemp[1], idemp[0])                              # needs sage.rings.number_field
                Free module generated by {} over Rational Field

            We recover the `2\times2` block of `\QQ[S_4]`
            corresponding to the unique simple module of dimension `2`
            of the symmetric group `S_4`::

                sage: # needs sage.combinat sage.groups sage.rings.number_field
                sage: A4 = SymmetricGroup(4).algebra(QQ)
                sage: e = A4.central_orthogonal_idempotents()[2]
                sage: A4.peirce_summand(e, e)
                Free module generated by {0, 1, 2, 3} over Rational Field

            TESTS:

            We check each idempotent belong to its own Peirce summand
            (see :issue:`24687`)::

                sage: # needs sage.combinat sage.groups sage.rings.number_field
                sage: from sage.monoids.hecke_monoid import HeckeMonoid
                sage: M = HeckeMonoid(SymmetricGroup(4))
                sage: A = M.algebra(QQ)
                sage: Idms = A.orthogonal_idempotents_central_mod_radical()
                sage: all(A.peirce_summand(e, e).retract(e)
                ....:     in A.peirce_summand(e, e) for e in Idms)
                True
            """
            B = self.basis()
            phi = self.module_morphism(on_basis=lambda k: ei * B[k] * ej,
                                       codomain=self, triangular='lower')
            ideal = phi.matrix(side='right').image()

            return self.submodule([self.from_vector(v) for v in ideal.basis()],
                                  already_echelonized=True)

        def peirce_decomposition(self, idempotents=None, check=True):
            r"""
            Return a Peirce decomposition of ``self``.

            Let `(e_i)_i` be a collection of orthogonal idempotents of
            `A` with sum `1`. The *Peirce decomposition* of `A` is the
            decomposition of `A` into the direct sum of the subspaces
            `e_i A e_j`.

            With the default collection of orthogonal idempotents, one has

            .. MATH::

                \dim e_i A e_j = C_{i,j} \dim S_i \dim S_j

            where `(S_i)_i` are the simple modules of `A` and
            `(C_{i,j})_{i, j}` is the Cartan invariants matrix.

            INPUT:

            - ``idempotents`` -- list of orthogonal idempotents
              `(e_i)_{i=0,\ldots,n}` of the algebra that sum to `1`
              (default: the idempotents returned by
              :meth:`orthogonal_idempotents_central_mod_radical`)

            - ``check`` -- boolean (default: ``True``); whether to check that
              the idempotents are indeed orthogonal and idempotent and
              sum to `1`

            OUTPUT:

            A list of lists `l` such that ``l[i][j]`` is the subspace
            `e_i A e_j`.

            .. SEEALSO::

                - :meth:`orthogonal_idempotents_central_mod_radical`
                - :meth:`cartan_invariants_matrix`

            EXAMPLES::

                sage: # needs sage.graphs sage.groups sage.modules sage.rings.number_field
                sage: A = Algebras(QQ).FiniteDimensional().WithBasis().example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: A.orthogonal_idempotents_central_mod_radical()
                (x, y)
                sage: decomposition = A.peirce_decomposition(); decomposition
                [[Free module generated by {0} over Rational Field,
                  Free module generated by {0, 1} over Rational Field],
                 [Free module generated by {} over Rational Field,
                  Free module generated by {0} over Rational Field]]
                sage: [ [[x.lift() for x in decomposition[i][j].basis()]
                ....:    for j in range(2)]
                ....:   for i in range(2)]
                [[[x], [a, b]],
                 [[], [y]]]

            We recover that the group algebra of the symmetric group
            `S_4` is a block matrix algebra::

                sage: # needs sage.groups sage.modules sage.rings.number_field
                sage: A = SymmetricGroup(4).algebra(QQ)
                sage: decomposition = A.peirce_decomposition()  # long time
                sage: [[decomposition[i][j].dimension()         # long time (4s)
                ....:   for j in range(len(decomposition))]
                ....:  for i in range(len(decomposition))]
                [[9, 0, 0, 0, 0],
                 [0, 9, 0, 0, 0],
                 [0, 0, 4, 0, 0],
                 [0, 0, 0, 1, 0],
                 [0, 0, 0, 0, 1]]

            The dimension of each block is `d^2`, where `d` is the
            dimension of the corresponding simple module of `S_4`. The
            latter are given by::

                sage: [p.standard_tableaux().cardinality() for p in Partitions(4)]      # needs sage.combinat
                [1, 3, 2, 3, 1]
            """
            if idempotents is None:
                idempotents = self.orthogonal_idempotents_central_mod_radical()
            if check:
                if not self.is_identity_decomposition_into_orthogonal_idempotents(idempotents):
                    raise ValueError("Not a decomposition of the identity into orthogonal idempotents")
            return [[self.peirce_summand(ei, ej) for ej in idempotents]
                    for ei in idempotents]

        def is_identity_decomposition_into_orthogonal_idempotents(self, l):
            r"""
            Return whether ``l`` is a decomposition of the identity
            into orthogonal idempotents.

            INPUT:

            - ``l`` -- list or iterable of elements of ``self``

            EXAMPLES::

                sage: # needs sage.graphs sage.modules
                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: x,y,a,b = A.algebra_generators(); x,y,a,b
                (x, y, a, b)
                sage: A.is_identity_decomposition_into_orthogonal_idempotents([A.one()])
                True
                sage: A.is_identity_decomposition_into_orthogonal_idempotents([x, y])
                True
                sage: A.is_identity_decomposition_into_orthogonal_idempotents([x + a, y - a])
                True

            Here the idempotents do not sum up to `1`::

                sage: A.is_identity_decomposition_into_orthogonal_idempotents([x])                  # needs sage.graphs sage.modules
                False

            Here `1+x` and `-x` are neither idempotent nor orthogonal::

                sage: A.is_identity_decomposition_into_orthogonal_idempotents([1 + x, -x])          # needs sage.graphs sage.modules
                False

            With the algebra of the `0`-Hecke monoid::

                sage: # needs sage.combinat sage.groups sage.modules sage.rings.number_field
                sage: from sage.monoids.hecke_monoid import HeckeMonoid
                sage: A = HeckeMonoid(SymmetricGroup(4)).algebra(QQ)
                sage: idempotents = A.orthogonal_idempotents_central_mod_radical()
                sage: A.is_identity_decomposition_into_orthogonal_idempotents(idempotents)
                True

            Here are some more counterexamples:

            1. Some orthogonal elements summing to `1` but not being
               idempotent::

                sage: # needs sage.libs.pari sage.modules
                sage: class PQAlgebra(CombinatorialFreeModule):
                ....:     def __init__(self, F, p):
                ....:         # Construct the quotient algebra F[x] / p,
                ....:         # where p is a univariate polynomial.
                ....:         R = parent(p); x = R.gen()
                ....:         I = R.ideal(p)
                ....:         self._xbar = R.quotient(I).gen()
                ....:         basis_keys = [self._xbar**i for i in range(p.degree())]
                ....:         CombinatorialFreeModule.__init__(self, F, basis_keys,
                ....:                 category=Algebras(F).FiniteDimensional().WithBasis())
                ....:     def x(self):
                ....:         return self(self._xbar)
                ....:     def one(self):
                ....:         return self.basis()[self.base_ring().one()]
                ....:     def product_on_basis(self, w1, w2):
                ....:         return self.from_vector(vector(w1*w2))
                sage: R.<x> = PolynomialRing(QQ)
                sage: A = PQAlgebra(QQ, x**3 - x**2 + x + 1); y = A.x()
                sage: a, b = y, 1 - y
                sage: A.is_identity_decomposition_into_orthogonal_idempotents((a, b))
                False

               For comparison::

                sage: # needs sage.libs.pari sage.modules
                sage: A = PQAlgebra(QQ, x**2 - x); y = A.x()
                sage: a, b = y, 1-y
                sage: A.is_identity_decomposition_into_orthogonal_idempotents((a, b))
                True
                sage: A.is_identity_decomposition_into_orthogonal_idempotents((a, A.zero(), b))
                True
                sage: A = PQAlgebra(QQ, x**3 - x**2 + x - 1); y = A.x()
                sage: a = (y**2 + 1) / 2
                sage: b = 1 - a
                sage: A.is_identity_decomposition_into_orthogonal_idempotents((a, b))
                True

            2. Some idempotents summing to 1 but not orthogonal::

                sage: # needs sage.libs.pari sage.modules
                sage: R.<x> = PolynomialRing(GF(2))
                sage: A = PQAlgebra(GF(2), x)
                sage: a = A.one()
                sage: A.is_identity_decomposition_into_orthogonal_idempotents((a,))
                True
                sage: A.is_identity_decomposition_into_orthogonal_idempotents((a, a, a))
                False

            3. Some orthogonal idempotents not summing to the identity::

                sage: # needs sage.libs.pari sage.modules
                sage: A.is_identity_decomposition_into_orthogonal_idempotents((a,a))
                False
                sage: A.is_identity_decomposition_into_orthogonal_idempotents(())
                False
            """
            return (self.sum(l) == self.one()
                    and all(e*e == e for e in l)
                    and all(e*f == 0 and f*e == 0 for i, e in enumerate(l)
                                                  for f in l[:i]))

        @cached_method
        def is_commutative(self) -> bool:
            """
            Return whether ``self`` is a commutative algebra.

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: S4 = SymmetricGroupAlgebra(QQ, 4)
                sage: S4.is_commutative()
                False
                sage: S2 = SymmetricGroupAlgebra(QQ, 2)
                sage: S2.is_commutative()
                True
            """
            B = list(self.basis())
            try:  # See if 1 is a basis element, if so, remove it
                B.remove(self.one())
            except ValueError:
                pass
            return all(b*bp == bp*b for i,b in enumerate(B) for bp in B[i+1:])

    class ElementMethods:

        def to_matrix(self, base_ring=None, action=operator.mul, side='left'):
            """
            Return the matrix of the action of ``self`` on the algebra.

            INPUT:

            - ``base_ring`` -- the base ring for the matrix to be constructed
            - ``action`` -- a bivariate function (default: :func:`operator.mul`)
            - ``side`` -- ``'left'`` or ``'right'`` (default: ``'left'``)

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
                sage: a = QS3([2,1,3])
                sage: a.to_matrix(side='left')
                [0 0 1 0 0 0]
                [0 0 0 0 1 0]
                [1 0 0 0 0 0]
                [0 0 0 0 0 1]
                [0 1 0 0 0 0]
                [0 0 0 1 0 0]
                sage: a.to_matrix(side='right')
                [0 0 1 0 0 0]
                [0 0 0 1 0 0]
                [1 0 0 0 0 0]
                [0 1 0 0 0 0]
                [0 0 0 0 0 1]
                [0 0 0 0 1 0]
                sage: a.to_matrix(base_ring=RDF, side='left')
                [0.0 0.0 1.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 0.0 1.0 0.0]
                [1.0 0.0 0.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 0.0 0.0 1.0]
                [0.0 1.0 0.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 1.0 0.0 0.0]

            AUTHORS: Mike Hansen, ...
            """
            basis = self.parent().basis()
            action_left = action
            if side == 'right':
                action = lambda x: action_left(basis[x], self)
            else:
                action = lambda x: action_left(self, basis[x])
            endo = self.parent().module_morphism(on_basis=action, codomain=self.parent())
            return endo.matrix(base_ring=base_ring)

        _matrix_ = to_matrix  # For temporary backward compatibility
        on_left_matrix = to_matrix

        def __invert__(self):
            r"""
            Return the inverse of ``self`` if it exists, and
            otherwise raise an error.

            .. WARNING::

                This always returns the inverse or fails on elements
                that are not invertible when the base ring is a field.
                In other cases, it may fail to find an inverse even
                if one exists if we cannot solve a linear system of
                equations over (the fraction field of) the base ring.

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
                sage: P = Permutation
                sage: a = 3 * QS3(P([1,2,3])) + QS3(P([1,3,2])) + QS3(P([2,1,3]))
                sage: b = ~a; b
                9/20*[1, 2, 3] - 7/40*[1, 3, 2] - 7/40*[2, 1, 3]
                 + 3/40*[2, 3, 1] + 3/40*[3, 1, 2] - 1/20*[3, 2, 1]
                sage: a * b
                [1, 2, 3]
                sage: ~b == a
                True

                sage: # needs sage.groups sage.modules
                sage: a = 3 * QS3.one()
                sage: b = ~a
                sage: b * a == QS3.one()
                True
                sage: b == 1/3 * QS3.one()
                True
                sage: ~b == a
                True

                sage: R.<t> = QQ[]
                sage: RS3 = SymmetricGroupAlgebra(R, 3)                                 # needs sage.groups sage.modules
                sage: a = RS3(P([1,2,3])) - RS3(P([1,3,2])) + RS3(P([2,1,3])); ~a       # needs sage.groups sage.modules
                -1/2*[1, 3, 2] + 1/2*[2, 1, 3] + 1/2*[2, 3, 1] + 1/2*[3, 1, 2]

            Some examples on elements that do not have an inverse::

                sage: c = 2 * QS3(P([1,2,3])) + QS3(P([1,3,2])) + QS3(P([2,1,3]))       # needs sage.groups sage.modules
                sage: ~c                                                                # needs sage.groups sage.modules
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= 2*[1, 2, 3] + [1, 3, 2] + [2, 1, 3])

                sage: # needs sage.groups sage.modules
                sage: ZS3 = SymmetricGroupAlgebra(ZZ, 3)
                sage: aZ = 3 * ZS3(P([1,2,3])) + ZS3(P([1,3,2])) + ZS3(P([2,1,3]))
                sage: ~aZ
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= 3*[1, 2, 3] + [1, 3, 2] + [2, 1, 3])
                sage: x = 2 * ZS3.one()
                sage: ~x
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= 2*[1, 2, 3])

            TESTS:

            An algebra that does not define ``one_basis()``::

                sage: # needs sage.combinat sage.groups sage.modules
                sage: I = DescentAlgebra(QQ, 3).I()
                sage: a = 3 * I.one()
                sage: ~a == 1/3 * I.one()
                True
            """
            alg = self.parent()
            R = alg.base_ring()
            ob = None
            try:
                ob = alg.one_basis()
            except (AttributeError, TypeError, ValueError):
                pass
            if ob is not None:
                mc = self.monomial_coefficients(copy=False)
                if len(mc) == 1 and ob in mc:
                    try:
                        return alg.term(ob, R(~mc[ob]))
                    except (ValueError, TypeError):
                        raise ValueError("cannot invert self (= %s)" % self)

            e = alg.one().to_vector()
            A = self.to_matrix()
            try:
                inv = A.solve_right(e)
                inv.change_ring(R)
                return alg.from_vector(inv)
            except (ValueError, TypeError):
                raise ValueError("cannot invert self (= %s)" % self)

    class Cellular(CategoryWithAxiom_over_base_ring):
        r"""
        Cellular algebras.

        Let `R` be a commutative ring. A `R`-algebra `A` is a
        *cellular algebra* if it has a *cell datum*, which is
        a tuple `(\Lambda, i, M, C)`, where `\Lambda` is finite
        poset with order `\ge`, if `\mu \in \Lambda` then `T(\mu)`
        is a finite set and

        .. MATH::

            C \colon \coprod_{\mu\in\Lambda}T(\mu) \times T(\mu)
              \longrightarrow A; (\mu,s,t) \mapsto c^\mu_{st}
              \text{ is an injective map}

        such that the following holds:

        * The set `\{c^\mu_{st}\mid \mu\in\Lambda, s,t\in T(\mu)\}` is a
          basis of `A`.
        * If `a \in A` and `\mu\in\Lambda, s,t \in T(\mu)` then:

          .. MATH::

              a c^\mu_{st} = \sum_{u\in T(\mu)} r_a(s,u) c^\mu_{ut}
              \pmod{A^{>\mu}},

          where `A^{>\mu}` is spanned by

          .. MATH::

              \{ c^\nu_{ab} \mid \nu > \mu \text{ and } a,b \in T(\nu) \}.

          Moreover, the scalar `r_a(s,u)` depends only on `a`, `s` and
          `u` and, in particular, is independent of `t`.

        * The map `\iota \colon A \longrightarrow A; c^\mu_{st} \mapsto
          c^\mu_{ts}` is an algebra anti-isomorphism.

        A *cellular  basis* for `A` is any basis of the form
        `\{c^\mu_{st} \mid \mu \in \Lambda, s,t \in T(\mu)\}`.

        Note that in particular, the scalars `r_a(u, s)` in the second
        condition do not depend on `t`.

        REFERENCES:

        - [GrLe1996]_
        - [KX1998]_
        - [Mat1999]_
        - :wikipedia:`Cellular_algebra`
        - http://webusers.imj-prg.fr/~bernhard.keller/ictp2006/lecturenotes/xi.pdf
        """
        class ParentMethods:
            def _test_cellular(self, **options):
                """
                Check that ``self`` satisfies the properties of a
                cellular algebra.

                EXAMPLES::

                    sage: S = SymmetricGroupAlgebra(QQ, 3)                              # needs sage.combinat sage.groups sage.modules
                    sage: S._test_cellular()                                            # needs sage.combinat sage.groups sage.modules
                """
                tester = self._tester(**options)
                cell_basis = self.cellular_basis()
                B = cell_basis.basis()
                P = self.cell_poset()
                for mu in P:
                    C = self.cell_module_indices(mu)
                    for s in C:
                        t = C[0]
                        vals = []
                        basis_elt = B[(mu, s, t)]
                        for a in B:
                            elt = a * basis_elt
                            tester.assertTrue( all(P.lt(i[0], mu) or i[2] == t
                                                   for i in elt.support()) )
                            vals.append([elt[(mu, u, t)] for u in C])
                        for t in C[1:]:
                            basis_elt = B[(mu, s, t)]
                            for i,a in enumerate(B):
                                elt = a * basis_elt
                                tester.assertTrue( all(P.lt(i[0], mu) or i[2] == t
                                                       for i in elt.support()) )
                                tester.assertEqual(vals[i], [elt[(mu, u, t)]
                                                             for u in C])

            @abstract_method
            def cell_poset(self):
                """
                Return the cell poset of ``self``.

                EXAMPLES::

                    sage: S = SymmetricGroupAlgebra(QQ, 4)                              # needs sage.groups sage.modules
                    sage: S.cell_poset()                                                # needs sage.groups sage.modules
                    Finite poset containing 5 elements
                """

            @abstract_method
            def cell_module_indices(self, mu):
                r"""
                Return the indices of the cell module of ``self``
                indexed by ``mu`` .

                This is the finite set `M(\lambda)`.

                EXAMPLES::

                    sage: S = SymmetricGroupAlgebra(QQ, 3)                              # needs sage.groups sage.modules
                    sage: S.cell_module_indices([2,1])                                  # needs sage.groups sage.modules
                    Standard tableaux of shape [2, 1]
                """

            @abstract_method(optional=True)
            def _to_cellular_element(self, i):
                """
                Return the image of the basis index ``i`` in the
                cellular basis of ``self``.

                EXAMPLES::

                    sage: S = SymmetricGroupAlgebra(QQ, 3)                                          # needs sage.groups sage.modules
                    sage: S._to_cellular_element   # no implementation currently uses this          # needs sage.groups sage.modules
                    NotImplemented
                """

            @abstract_method(optional=True)
            def _from_cellular_index(self, x):
                """
                Return the image in ``self`` from the index of the
                cellular basis ``x``.

                EXAMPLES::

                    sage: # needs sage.combinat sage.groups sage.modules
                    sage: S = SymmetricGroupAlgebra(QQ, 3)
                    sage: mu = Partition([2,1])
                    sage: s = StandardTableau([[1,2],[3]])
                    sage: t = StandardTableau([[1,3],[2]])
                    sage: S._from_cellular_index((mu, s, t))
                    1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1]
                """

            def cellular_involution(self, x):
                """
                Return the cellular involution of ``x`` in ``self``.

                EXAMPLES::

                    sage: S = SymmetricGroupAlgebra(QQ, 3)                              # needs sage.groups sage.modules
                    sage: for b in S.basis(): b, S.cellular_involution(b)               # needs sage.groups sage.modules
                    ([1, 2, 3], [1, 2, 3])
                    ([1, 3, 2], 49/48*[1, 3, 2] + 7/48*[2, 3, 1]
                                - 7/48*[3, 1, 2] - 1/48*[3, 2, 1])
                    ([2, 1, 3], [2, 1, 3])
                    ([2, 3, 1], -7/48*[1, 3, 2] - 1/48*[2, 3, 1]
                                 + 49/48*[3, 1, 2] + 7/48*[3, 2, 1])
                    ([3, 1, 2], 7/48*[1, 3, 2] + 49/48*[2, 3, 1]
                                 - 1/48*[3, 1, 2] - 7/48*[3, 2, 1])
                    ([3, 2, 1], -1/48*[1, 3, 2] - 7/48*[2, 3, 1]
                                 + 7/48*[3, 1, 2] + 49/48*[3, 2, 1])
                """
                C = self.cellular_basis()
                if C is self:
                    M = x.monomial_coefficients(copy=False)
                    return self._from_dict({(i[0], i[2], i[1]): M[i] for i in M},
                                           remove_zeros=False)
                return self(C(x).cellular_involution())

            @cached_method
            def cells(self):
                """
                Return the cells of ``self``.

                EXAMPLES::

                    sage: S = SymmetricGroupAlgebra(QQ, 3)                              # needs sage.groups sage.modules
                    sage: dict(S.cells())                                               # needs sage.groups sage.modules
                    {[1, 1, 1]: Standard tableaux of shape [1, 1, 1],
                     [2, 1]: Standard tableaux of shape [2, 1],
                     [3]: Standard tableaux of shape [3]}
                """
                from sage.sets.family import Family
                return Family(self.cell_poset(), self.cell_module_indices)

            def cellular_basis(self):
                """
                Return the cellular basis of ``self``.

                EXAMPLES::

                    sage: S = SymmetricGroupAlgebra(QQ, 3)                              # needs sage.groups sage.modules
                    sage: S.cellular_basis()                                            # needs sage.groups sage.modules
                    Cellular basis of Symmetric group algebra of order 3
                     over Rational Field
                """
                from sage.algebras.cellular_basis import CellularBasis
                return CellularBasis(self)

            def cell_module(self, mu, **kwds):
                """
                Return the cell module indexed by ``mu``.

                EXAMPLES::

                    sage: S = SymmetricGroupAlgebra(QQ, 3)                              # needs sage.groups sage.modules
                    sage: S.cell_module(Partition([2,1]))                               # needs sage.combinat sage.groups sage.modules
                    Cell module indexed by [2, 1] of Cellular basis of
                     Symmetric group algebra of order 3 over Rational Field
                """
                from sage.modules.with_basis.cell_module import CellModule
                return CellModule(self.cellular_basis(), mu, **kwds)

            @cached_method
            def simple_module_parameterization(self):
                r"""
                Return a parameterization of the simple modules of ``self``.

                The set of simple modules are parameterized by
                `\lambda \in \Lambda` such that the cell module
                bilinear form `\Phi_{\lambda} \neq 0`.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: TL = TemperleyLiebAlgebra(5, 30, QQ)  # semisimple
                    sage: len(TL.radical_basis())
                    0
                    sage: TL.simple_module_parameterization()
                    (1, 3, 5)

                    sage: # needs sage.modules
                    sage: TL = TemperleyLiebAlgebra(5, 1, QQ)  # not semisimple
                    sage: len(TL.radical_basis())
                    24
                    sage: TL.simple_module_parameterization()
                    (1, 3, 5)

                    sage: # needs sage.modules
                    sage: TL = TemperleyLiebAlgebra(6, 30, QQ)  # semisimple
                    sage: all(TL.cell_module(la).dimension()
                    ....:     == TL.cell_module(la).simple_module().dimension()
                    ....:     for la in TL.simple_module_parameterization())
                    True
                    sage: TL.simple_module_parameterization()
                    (0, 2, 4, 6)

                    sage: # needs sage.modules
                    sage: TL = TemperleyLiebAlgebra(6, 0, QQ)  # not semisimple
                    sage: TL.simple_module_parameterization()
                    (2, 4, 6)
                """
                return tuple([mu for mu in self.cell_poset()
                              if self.cell_module(mu).nonzero_bilinear_form()])

        class ElementMethods:
            def cellular_involution(self):
                """
                Return the cellular involution on ``self``.

                EXAMPLES::

                    sage: # needs sage.groups sage.modules
                    sage: S = SymmetricGroupAlgebra(QQ, 4)
                    sage: elt = S([3,1,2,4])
                    sage: ci = elt.cellular_involution(); ci
                    7/48*[1, 3, 2, 4] + 49/48*[2, 3, 1, 4]
                     - 1/48*[3, 1, 2, 4] - 7/48*[3, 2, 1, 4]
                    sage: ci.cellular_involution()
                    [3, 1, 2, 4]
                """
                return self.parent().cellular_involution(self)

        class TensorProducts(TensorProductsCategory):
            """
            The category of cellular algebras constructed by tensor
            product of cellular algebras.
            """
            @cached_method
            def extra_super_categories(self):
                """
                Tensor products of cellular algebras are cellular.

                EXAMPLES::

                    sage: cat = Algebras(QQ).FiniteDimensional().WithBasis()
                    sage: cat.Cellular().TensorProducts().extra_super_categories()
                    [Category of finite dimensional cellular algebras with basis
                     over Rational Field]
                """
                return [self.base_category()]

            class ParentMethods:
                @cached_method
                def cell_poset(self):
                    """
                    Return the cell poset of ``self``.

                    EXAMPLES::

                        sage: # needs sage.groups sage.modules
                        sage: S2 = SymmetricGroupAlgebra(QQ, 2)
                        sage: S3 = SymmetricGroupAlgebra(QQ, 3)
                        sage: T = S2.tensor(S3)
                        sage: T.cell_poset()                                            # needs sage.combinat sage.graphs
                        Finite poset containing 6 elements
                    """
                    ret = self._sets[0].cell_poset()
                    for A in self._sets[1:]:
                        ret = ret.product(A.cell_poset())
                    return ret

                def cell_module_indices(self, mu):
                    r"""
                    Return the indices of the cell module of ``self``
                    indexed by ``mu`` .

                    This is the finite set `M(\lambda)`.

                    EXAMPLES::

                        sage: # needs sage.groups sage.modules
                        sage: S2 = SymmetricGroupAlgebra(QQ, 2)
                        sage: S3 = SymmetricGroupAlgebra(QQ, 3)
                        sage: T = S2.tensor(S3)
                        sage: T.cell_module_indices(([1,1], [2,1]))
                        The Cartesian product of (Standard tableaux of shape [1, 1],
                                                  Standard tableaux of shape [2, 1])
                    """
                    from sage.categories.cartesian_product import cartesian_product
                    return cartesian_product([self._sets[i].cell_module_indices(x)
                                              for i,x in enumerate(mu)])

                @lazy_attribute
                def cellular_involution(self):
                    """
                    Return the image of the cellular involution of the basis
                    element indexed by ``i``.

                    EXAMPLES::

                        sage: # needs sage.groups sage.modules
                        sage: S2 = SymmetricGroupAlgebra(QQ, 2)
                        sage: S3 = SymmetricGroupAlgebra(QQ, 3)
                        sage: T = S2.tensor(S3)
                        sage: for b in T.basis(): b, T.cellular_involution(b)
                        ([1, 2] # [1, 2, 3], [1, 2] # [1, 2, 3])
                        ([1, 2] # [1, 3, 2],
                         49/48*[1, 2] # [1, 3, 2] + 7/48*[1, 2] # [2, 3, 1]
                          - 7/48*[1, 2] # [3, 1, 2] - 1/48*[1, 2] # [3, 2, 1])
                        ([1, 2] # [2, 1, 3], [1, 2] # [2, 1, 3])
                        ([1, 2] # [2, 3, 1],
                         -7/48*[1, 2] # [1, 3, 2] - 1/48*[1, 2] # [2, 3, 1]
                          + 49/48*[1, 2] # [3, 1, 2] + 7/48*[1, 2] # [3, 2, 1])
                        ([1, 2] # [3, 1, 2],
                         7/48*[1, 2] # [1, 3, 2] + 49/48*[1, 2] # [2, 3, 1]
                          - 1/48*[1, 2] # [3, 1, 2] - 7/48*[1, 2] # [3, 2, 1])
                        ([1, 2] # [3, 2, 1],
                         -1/48*[1, 2] # [1, 3, 2] - 7/48*[1, 2] # [2, 3, 1]
                          + 7/48*[1, 2] # [3, 1, 2] + 49/48*[1, 2] # [3, 2, 1])
                        ([2, 1] # [1, 2, 3], [2, 1] # [1, 2, 3])
                        ([2, 1] # [1, 3, 2],
                         49/48*[2, 1] # [1, 3, 2] + 7/48*[2, 1] # [2, 3, 1]
                          - 7/48*[2, 1] # [3, 1, 2] - 1/48*[2, 1] # [3, 2, 1])
                        ([2, 1] # [2, 1, 3], [2, 1] # [2, 1, 3])
                        ([2, 1] # [2, 3, 1],
                         -7/48*[2, 1] # [1, 3, 2] - 1/48*[2, 1] # [2, 3, 1]
                          + 49/48*[2, 1] # [3, 1, 2] + 7/48*[2, 1] # [3, 2, 1])
                        ([2, 1] # [3, 1, 2],
                         7/48*[2, 1] # [1, 3, 2] + 49/48*[2, 1] # [2, 3, 1]
                          - 1/48*[2, 1] # [3, 1, 2] - 7/48*[2, 1] # [3, 2, 1])
                        ([2, 1] # [3, 2, 1],
                         -1/48*[2, 1] # [1, 3, 2] - 7/48*[2, 1] # [2, 3, 1]
                          + 7/48*[2, 1] # [3, 1, 2] + 49/48*[2, 1] # [3, 2, 1])
                    """
                    if self.cellular_basis() is self:
                        def func(x):
                            M = x.monomial_coefficients(copy=False)
                            return self._from_dict({(i[0], i[2], i[1]): M[i] for i in M},
                                                   remove_zeros=False)
                        return self.module_morphism(function=func, codomain=self)

                    def on_basis(i):
                        return self._tensor_of_elements([A.basis()[i[j]].cellular_involution()
                                                         for j, A in enumerate(self._sets)])
                    return self.module_morphism(on_basis, codomain=self)

                @cached_method
                def _to_cellular_element(self, i):
                    """
                    Return the image of the basis index ``i`` in the
                    cellular basis of ``self``.

                    EXAMPLES::

                        sage: # needs sage.groups sage.modules
                        sage: S2 = SymmetricGroupAlgebra(QQ, 2)
                        sage: S3 = SymmetricGroupAlgebra(QQ, 3)
                        sage: T = S2.tensor(S3)
                        sage: all(T(T._to_cellular_element(k)).leading_support() == k
                        ....:     for k in T.basis().keys())
                        True
                    """
                    C = [A.cellular_basis() for A in self._sets]
                    elts = [C[j](self._sets[j].basis()[ij]) for j, ij in enumerate(i)]
                    from sage.categories.tensor import tensor
                    T = tensor(C)
                    temp = T._tensor_of_elements(elts)
                    B = self.cellular_basis()
                    M = temp.monomial_coefficients(copy=False)

                    def convert_index(i):
                        mu = []
                        s = []
                        t = []
                        for a, b, c in i:
                            mu.append(a)
                            s.append(b)
                            t.append(c)
                        C = self.cell_module_indices(mu)
                        return (tuple(mu), C(s), C(t))
                    return B._from_dict({convert_index(i): M[i] for i in M},
                                        remove_zeros=False)

                @cached_method
                def _from_cellular_index(self, x):
                    """
                    Return the image in ``self`` from the index of the
                    cellular basis ``x``.

                    EXAMPLES::

                        sage: # needs sage.groups sage.modules
                        sage: S2 = SymmetricGroupAlgebra(QQ, 2)
                        sage: S3 = SymmetricGroupAlgebra(QQ, 3)
                        sage: T = S2.tensor(S3)
                        sage: C = T.cellular_basis()
                        sage: all(C(T._from_cellular_index(k)).leading_support() == k
                        ....:     for k in C.basis().keys())
                        True
                    """
                    elts = [A(A.cellular_basis().basis()[ (x[0][i], x[1][i], x[2][i]) ])
                            for i,A in enumerate(self._sets)]
                    return self._tensor_of_elements(elts)

    class SubcategoryMethods:
        @cached_method
        def Cellular(self):
            """
            Return the full subcategory of the cellular objects
            of ``self``.

            .. SEEALSO:: :wikipedia:`Cellular_algebra`

            EXAMPLES::

                sage: Algebras(QQ).FiniteDimensional().WithBasis().Cellular()
                Category of finite dimensional cellular algebras with basis
                 over Rational Field

            TESTS::

                sage: cat = Algebras(QQ).FiniteDimensional().WithBasis()
                sage: TestSuite(cat.Cellular()).run()
                sage: HopfAlgebras(QQ).FiniteDimensional().WithBasis().Cellular.__module__
                'sage.categories.finite_dimensional_algebras_with_basis'
            """
            return self._with_axiom('Cellular')
