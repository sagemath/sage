r"""
Modules With Basis

AUTHORS:

- Nicolas M. Thiery (2008-2014): initial revision, axiomatization
- Jason Bandlow and Florent Hivert (2010): Triangular Morphisms
- Christian Stump (2010): :trac:`9648` module_morphism's to a wider class
  of codomains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.misc.cachefunc import cached_method
from sage.misc.sage_itertools import max_cmp, min_cmp
from sage.categories.homsets import HomsetsCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.tensor import tensor, TensorProductsCategory
from sage.categories.dual import DualObjectsCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.modules import Modules
from sage.structure.element import Element, parent
from sage.misc.lazy_import import lazy_import
lazy_import('sage.modules.with_basis.morphism',
            ['ModuleMorphismByLinearity',
             'ModuleMorphismFromMatrix',
             'ModuleMorphismFromFunction',
             'DiagonalModuleMorphism',
             'TriangularModuleMorphismByLinearity',
             'TriangularModuleMorphismFromFunction'])

class ModulesWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of modules with a distinguished basis.

    The elements are represented by expanding them in the distinguished basis.
    The morphisms are not required to respect the distinguished basis.

    EXAMPLES::

        sage: ModulesWithBasis(ZZ)
        Category of modules with basis over Integer Ring
        sage: ModulesWithBasis(ZZ).super_categories()
        [Category of modules over Integer Ring]

    If the base ring is actually a field, this constructs instead the
    category of vector spaces with basis::

        sage: ModulesWithBasis(QQ)
        Category of vector spaces with basis over Rational Field

        sage: ModulesWithBasis(QQ).super_categories()
        [Category of modules with basis over Rational Field,
         Category of vector spaces over Rational Field]

    Let `X` and `Y` be two modules with basis. We can build `Hom(X,Y)`::

        sage: X = CombinatorialFreeModule(QQ, [1,2]); X.__custom_name = "X"
        sage: Y = CombinatorialFreeModule(QQ, [3,4]); Y.__custom_name = "Y"
        sage: H = Hom(X, Y); H
        Set of Morphisms from X to Y in Category of finite dimensional vector spaces with basis over Rational Field

    The simplest morphism is the zero map::

        sage: H.zero()         # todo: move this test into module once we have an example
        Generic morphism:
          From: X
          To:   Y

    which we can apply to elements of `X`::

        sage: x = X.monomial(1) + 3 * X.monomial(2)
        sage: H.zero()(x)
        0

    TESTS::

        sage: f = H.zero().on_basis()
        sage: f(1)
        0
        sage: f(2)
        0

    EXAMPLES:

    We now construct a more interesting morphism by extending a
    function by linearity::

        sage: phi = H(on_basis = lambda i: Y.monomial(i+2)); phi
        Generic morphism:
          From: X
          To:   Y
        sage: phi(x)
        B[3] + 3*B[4]

    We can retrieve the function acting on indices of the basis::

        sage: f = phi.on_basis()
        sage: f(1), f(2)
        (B[3], B[4])

    `Hom(X,Y)` has a natural module structure (except for the zero,
    the operations are not yet implemented though). However since the
    dimension is not necessarily finite, it is not a module with
    basis; but see :class:`FiniteDimensionalModulesWithBasis` and
    :class:`GradedModulesWithBasis`::

        sage: H in ModulesWithBasis(QQ), H in Modules(QQ)
        (False, True)

    Some more playing around with categories and higher order homsets::

        sage: H.category()
        Category of homsets of modules with basis over Rational Field
        sage: Hom(H, H).category()
        Category of endsets of homsets of modules with basis over Rational Field

    .. TODO:: ``End(X)`` is an algebra.

    TESTS::

        sage: TestSuite(ModulesWithBasis(ZZ)).run()
    """

    def _call_(self, x):
        """
        Construct a module with basis (resp. vector space) from the data in ``x``.

        EXAMPLES::

            sage: CZ = ModulesWithBasis(ZZ); CZ
            Category of modules with basis over Integer Ring
            sage: CQ = ModulesWithBasis(QQ); CQ
            Category of vector spaces with basis over Rational Field

        ``x`` is returned unchanged if it is already in this category::

            sage: CZ(CombinatorialFreeModule(ZZ, ('a','b','c')))
            Free module generated by {'a', 'b', 'c'} over Integer Ring
            sage: CZ(ZZ^3)
            Ambient free module of rank 3 over the principal ideal domain Integer Ring

        If needed (and possible) the base ring is changed appropriately::

            sage: CQ(ZZ^3)                       # indirect doctest
            Vector space of dimension 3 over Rational Field

        If ``x`` itself is not a module with basis, but there is a
        canonical one associated to it, the later is returned::

            sage: CQ(AbelianVariety(Gamma0(37))) # indirect doctest
            Vector space of dimension 4 over Rational Field
        """
        try:
            M = x.free_module()
            if M.base_ring() != self.base_ring():
                M = M.change_ring(self.base_ring())
        except (TypeError, AttributeError) as msg:
            raise TypeError("%s\nunable to coerce x (=%s) into %s"%(msg,x,self))
        return M

    def is_abelian(self):
        """
        Returns whether this category is abelian.

        This is the case if and only if the base ring is a field.

        EXAMPLES::

            sage: ModulesWithBasis(QQ).is_abelian()
            True
            sage: ModulesWithBasis(ZZ).is_abelian()
            False
        """
        return self.base_ring().is_field()

    FiniteDimensional = LazyImport('sage.categories.finite_dimensional_modules_with_basis', 'FiniteDimensionalModulesWithBasis')
    Graded = LazyImport('sage.categories.graded_modules_with_basis', 'GradedModulesWithBasis')

    class ParentMethods:
        @cached_method
        def basis(self):
            """
            Return the basis of ``self``.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
                sage: F.basis()
                Finite family {'a': B['a'], 'c': B['c'], 'b': B['b']}

            ::

                sage: QS3 = SymmetricGroupAlgebra(QQ,3)
                sage: list(QS3.basis())
                [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
            """
            from sage.combinat.family import Family
            return Family(self._indices, self.monomial)

        def module_morphism(self, on_basis=None, matrix=None, function=None,
                            diagonal=None, triangular=None, unitriangular=False,
                            **keywords):
            r"""
            Construct a module morphism from ``self`` to ``codomain``.

            Let ``self`` be a module `X` with a basis indexed by `I`.
            This constructs a morphism `f: X \to Y` by linearity from
            a map `I \to Y` which is to be its restriction to the
            basis `(x_i)_{i \in I}` of `X`. Some variants are possible
            too.

            INPUT:

            - ``self`` -- a parent `X` in ``ModulesWithBasis(R)`` with
              basis `x=(x_i)_{i\in I}`.

            Exactly one of the four following options must be
            specified in order to define the morphism:

            - ``on_basis`` -- a function `f` from `I` to `Y`
            - ``diagonal`` -- a function `d` from `I` to `R`
            - ``function`` -- a function `f` from `X` to `Y`
            - ``matrix``   -- a matrix of size `\dim X \times \dim Y` or `\dim Y \times \dim X`

            Further options include:

            - ``codomain`` -- the codomain `Y` of the morphism (default:
              ``f.codomain()`` if it's defined; otherwise it must be specified)

            - ``category`` -- a category or ``None`` (default: `None``)

            - ``zero`` -- the zero of the codomain (default: ``codomain.zero()``);
              can be used (with care) to define affine maps.
              Only meaningful with ``on_basis``.

            - ``position`` -- a non-negative integer specifying which
              positional argument in used as the input of the function `f`
              (default: 0); this is currently only used with ``on_basis``.

            - ``triangular`` --  (default: ``None``) ``"upper"`` or
              ``"lower"`` or ``None``:

              * ``"upper"`` - if the
                :meth:`~ModulesWithBasis.ElementMethods.leading_support`
                of the image of the basis vector `x_i` is `i`, or

              * ``"lower"`` - if the
                :meth:`~ModulesWithBasis.ElementMethods.trailing_support`
                of the image of the basis vector `x_i` is `i`.

            - ``unitriangular`` -- (default: ``False``) a boolean.
              Only meaningful for a triangular morphism.
              As a shorthand, one may use ``unitriangular="lower"``
              for ``triangular="lower", unitriangular=True``.

            - ``side`` -- "left" or "right" (default: "left")
              Only meaningful for a morphism built from a matrix.

            EXAMPLES:

            With the ``on_basis`` option, this returns a function `g`
            obtained by extending `f` by linearity on the
            ``position``-th positional argument. For example, for
            ``position == 1`` and a ternary function `f`, one has:

            .. MATH::

                g\left( a,\ \sum_i \lambda_i x_i,\ c \right)
                = \sum_i \lambda_i f(a, i, c).

            ::

                sage: X = CombinatorialFreeModule(QQ, [1,2,3]);   X.rename("X")
                sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4]); Y.rename("Y")
                sage: phi = X.module_morphism(lambda i: Y.monomial(i) + 2*Y.monomial(i+1), codomain = Y)
                sage: x = X.basis(); y = Y.basis()
                sage: phi(x[1] + x[3])
                B[1] + 2*B[2] + B[3] + 2*B[4]

                sage: phi
                Generic morphism:
                From: X
                To:   Y

            By default, the category is the first of
            ``Modules(R).WithBasis().FiniteDimensional()``,
            ``Modules(R).WithBasis()``, ``Modules(R)``, and
            ``CommutativeAdditiveMonoids()`` that contains both the
            domain and the codomain::

                sage: phi.category_for()
                Category of finite dimensional vector spaces with basis over Rational Field

            With the ``zero`` argument, one can define affine morphisms::

                sage: phi = X.module_morphism(lambda i: Y.monomial(i) + 2*Y.monomial(i+1),
                ....:                         codomain = Y, zero = 10*y[1])
                sage: phi(x[1] + x[3])
                11*B[1] + 2*B[2] + B[3] + 2*B[4]

            In this special case, the default category is ``Sets()``::

                sage: phi.category_for()
                Category of sets

            One can construct morphisms with the base ring as codomain::

                sage: X = CombinatorialFreeModule(ZZ,[1,-1])
                sage: phi = X.module_morphism( on_basis=lambda i: i, codomain=ZZ )
                sage: phi( 2 * X.monomial(1) + 3 * X.monomial(-1) )
                -1
                sage: phi.category_for()
                Category of commutative additive semigroups
                sage: phi.category_for() # todo: not implemented (ZZ is currently not in Modules(ZZ))
                Category of modules over Integer Ring

            Or more generaly any ring admitting a coercion map from the base ring::

                sage: phi = X.module_morphism(on_basis=lambda i: i, codomain=RR )
                sage: phi( 2 * X.monomial(1) + 3 * X.monomial(-1) )
                -1.00000000000000
                sage: phi.category_for()
                Category of commutative additive semigroups
                sage: phi.category_for() # todo: not implemented (RR is currently not in Modules(ZZ))
                Category of modules over Integer Ring

                sage: phi = X.module_morphism(on_basis=lambda i: i, codomain=Zmod(4) )
                sage: phi( 2 * X.monomial(1) + 3 * X.monomial(-1) )
                3

                sage: phi = Y.module_morphism(on_basis=lambda i: i, codomain=Zmod(4) )
                Traceback (most recent call last):
                ...
                ValueError: codomain(=Ring of integers modulo 4) should be a module over the base ring of the domain(=Y)

            On can also define module morphisms between free modules
            over different base rings; here we implement the natural
            map from `X = \RR^2` to `Y = \CC`::

                sage: X = CombinatorialFreeModule(RR,['x','y'])
                sage: Y = CombinatorialFreeModule(CC,['z'])
                sage: x = X.monomial('x')
                sage: y = X.monomial('y')
                sage: z = Y.monomial('z')
                sage: def on_basis( a ):
                ....:     if a == 'x':
                ....:         return CC(1) * z
                ....:     elif a == 'y':
                ....:         return CC(I) * z
                sage: phi = X.module_morphism( on_basis=on_basis, codomain=Y )
                sage: v = 3 * x + 2 * y; v
                3.00000000000000*B['x'] + 2.00000000000000*B['y']
                sage: phi(v)
                (3.00000000000000+2.00000000000000*I)*B['z']
                sage: phi.category_for()
                Category of commutative additive semigroups
                sage: phi.category_for() # todo: not implemented (CC is currently not in Modules(RR)!)
                Category of vector spaces over Real Field with 53 bits of precision

                sage: Y = CombinatorialFreeModule(CC['q'],['z'])
                sage: z = Y.monomial('z')
                sage: phi = X.module_morphism( on_basis=on_basis, codomain=Y )
                sage: phi(v)
                (3.00000000000000+2.00000000000000*I)*B['z']

            Of course, there should be a coercion between the
            respective base rings of the domain and the codomain for
            this to be meaningful::

                sage: Y = CombinatorialFreeModule(QQ,['z'])
                sage: phi = X.module_morphism( on_basis=on_basis, codomain=Y )
                Traceback (most recent call last):
                ...
                ValueError: codomain(=Free module generated by {'z'} over Rational Field) should be a module over the base ring of the domain(=Free module generated by {'x', 'y'} over Real Field with 53 bits of precision)

                sage: Y = CombinatorialFreeModule(RR['q'],['z'])
                sage: phi = Y.module_morphism( on_basis=on_basis, codomain=X )
                Traceback (most recent call last):
                ...
                ValueError: codomain(=Free module generated by {'x', 'y'} over Real Field with 53 bits of precision) should be a module over the base ring of the domain(=Free module generated by {'z'} over Univariate Polynomial Ring in q over Real Field with 53 bits of precision)


            With the ``diagonal=d`` argument, this constructs the
            module morphism `g` such that

            .. MATH::

                `g(x_i) = d(i) y_i`

            This assumes that the respective bases `x` and `y` of `X`
            and `Y` have the same index set `I`::

                sage: X = CombinatorialFreeModule(ZZ, [1,2,3]); X.rename("X")
                sage: phi = X.module_morphism(diagonal=factorial, codomain=X)
                sage: x = X.basis()
                sage: phi(x[1]), phi(x[2]), phi(x[3])
                (B[1], 2*B[2], 6*B[3])

            See also: :class:`sage.modules.with_basis.morphism.DiagonalModuleMorphism`.

            With the ``matrix=m`` argument, this constructs the module
            morphism whose matrix in the distinguished basis of `X`
            and `Y` is `m`::

                sage: X = CombinatorialFreeModule(ZZ, [1,2,3]); X.rename("X"); x = X.basis()
                sage: Y = CombinatorialFreeModule(ZZ, [3,4]); Y.rename("Y"); y = Y.basis()
                sage: m = matrix([[0,1,2],[3,5,0]])
                sage: phi = X.module_morphism(matrix=m, codomain=Y)
                sage: phi(x[1])
                3*B[4]
                sage: phi(x[2])
                B[3] + 5*B[4]


            See also: :class:`sage.modules.with_basis.morphism.ModuleMorphismFromMatrix`.

            With ``triangular="upper"``, the constructed module morphism is
            assumed to be upper triangular; that is its matrix in the
            distinguished basis of `X` and `Y` would be upper triangular with
            invertible elements on its diagonal. This is used to compute
            preimages and to invert the morphism::

                sage: I = range(1,200)
                sage: X = CombinatorialFreeModule(QQ, I); X.rename("X"); x = X.basis()
                sage: Y = CombinatorialFreeModule(QQ, I); Y.rename("Y"); y = Y.basis()
                sage: f = Y.sum_of_monomials * divisors
                sage: phi = X.module_morphism(f, triangular="upper", codomain = Y)
                sage: phi(x[2])
                B[1] + B[2]
                sage: phi(x[6])
                B[1] + B[2] + B[3] + B[6]
                sage: phi(x[30])
                B[1] + B[2] + B[3] + B[5] + B[6] + B[10] + B[15] + B[30]
                sage: phi.preimage(y[2])
                -B[1] + B[2]
                sage: phi.preimage(y[6])
                B[1] - B[2] - B[3] + B[6]
                sage: phi.preimage(y[30])
                -B[1] + B[2] + B[3] + B[5] - B[6] - B[10] - B[15] + B[30]
                sage: (phi^-1)(y[30])
                -B[1] + B[2] + B[3] + B[5] - B[6] - B[10] - B[15] + B[30]

            Since :trac:`8678`, one can also define a triangular
            morphism from a function::

                sage: X = CombinatorialFreeModule(QQ, [0,1,2,3,4]); x = X.basis()
                sage: from sage.modules.with_basis.morphism import TriangularModuleMorphismFromFunction
                sage: def f(x): return x + X.term(0, sum(x.coefficients()))
                sage: phi = X.module_morphism(function=f, codomain=X, triangular="upper")
                sage: phi(x[2] + 3*x[4])
                4*B[0] + B[2] + 3*B[4]
                sage: phi.preimage(_)
                B[2] + 3*B[4]

            For details and further optional arguments, see
            :class:`sage.modules.with_basis.morphism.TriangularModuleMorphism`.

            .. WARNING::

                As a temporary measure, until multivariate morphisms
                are implemented, the constructed morphism is in
                ``Hom(codomain, domain, category)``. This is only
                correct for unary functions.

            .. TODO::

               - Should codomain be ``self`` by default in the
                 diagonal, triangular, and matrix cases?

               - Support for diagonal morphisms between modules not
                 sharing the same index set

            TESTS::

                sage: X = CombinatorialFreeModule(ZZ, [1,2,3]); X.rename("X")
                sage: phi = X.module_morphism(codomain=X)
                Traceback (most recent call last):
                ...
                ValueError: module_morphism() takes exactly one option
                out of `matrix`, `on_basis`, `function`, `diagonal`

            ::

                sage: X = CombinatorialFreeModule(ZZ, [1,2,3]); X.rename("X")
                sage: phi = X.module_morphism(diagonal=factorial, matrix=matrix(), codomain=X)
                Traceback (most recent call last):
                ...
                ValueError: module_morphism() takes exactly one option
                out of `matrix`, `on_basis`, `function`, `diagonal`

            ::

                sage: X = CombinatorialFreeModule(ZZ, [1,2,3]); X.rename("X")
                sage: phi = X.module_morphism(matrix=factorial, codomain=X)
                Traceback (most recent call last):
                ...
                ValueError: matrix (=factorial) should be a matrix

            ::

                sage: X = CombinatorialFreeModule(ZZ, [1,2,3]); X.rename("X")
                sage: phi = X.module_morphism(diagonal=3, codomain=X)
                Traceback (most recent call last):
                ...
                ValueError: diagonal (=3) should be a function

            """
            if not len([x for x in [matrix, on_basis, function, diagonal] if x is not None]) == 1:
                raise ValueError("module_morphism() takes exactly one option out of `matrix`, `on_basis`, `function`, `diagonal`")
            if matrix is not None:
                return ModuleMorphismFromMatrix(domain=self, matrix=matrix, **keywords)
            if diagonal is not None:
                return DiagonalModuleMorphism(domain=self, diagonal=diagonal, **keywords)
            if unitriangular in ["upper", "lower"] and triangular is None:
                triangular = unitriangular
                unitriangular = True
            if triangular is not None:
                if on_basis is not None:
                    return TriangularModuleMorphismByLinearity(
                        domain=self, on_basis=on_basis,
                        triangular=triangular, unitriangular=unitriangular,
                        **keywords)
                else:
                    return TriangularModuleMorphismFromFunction(
                        domain=self, function=function,
                        triangular=triangular, unitriangular=unitriangular,
                        **keywords)
            if on_basis is not None:
                return ModuleMorphismByLinearity(
                    domain=self, on_basis=on_basis, **keywords)
            else:
                return ModuleMorphismFromFunction( # Or just SetMorphism?
                    domain=self, function=function, **keywords)

        _module_morphism = module_morphism

        def _repr_(self):
            """
            EXAMPLES::

                sage: class FooBar(CombinatorialFreeModule): pass
                sage: C = FooBar(QQ, (1,2,3)); C # indirect doctest
                Free module generated by {1, 2, 3} over Rational Field

                sage: C._name = "foobar"; C
                foobar over Rational Field

                sage: C.rename("barfoo"); C
                barfoo

                sage: class FooBar(Parent):
                ....:     def basis(self): return Family({1:"foo", 2:"bar"})
                ....:     def base_ring(self): return QQ
                sage: FooBar(category = ModulesWithBasis(QQ))
                Free module generated by [1, 2] over Rational Field
            """
            if hasattr(self, "_name"):
                name = self._name
            else:
                name = "Free module generated by {}".format(self.basis().keys())
            return name + " over {}".format(self.base_ring())

        def echelon_form(self, elements):
            r"""
            Return a basis in echelon form of the subspace spanned by
            a finite set of elements.

            INPUT:

            - ``elements`` -- a list or finite iterable of elements of ``self``.

            OUTPUT:

            A list of elements of ``self`` whose expressions as
            vectors form a matrix in echelon form. If ``base_ring`` is
            specified, then the calculation is achieved in this base
            ring.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x")
                sage: x = X.basis()
                sage: V = X.echelon_form([x[0]-x[1], x[0]-x[2],x[1]-x[2]]); V
                [x[0] - x[2], x[1] - x[2]]
                sage: matrix(map(vector, V))
                [ 1  0 -1]
                [ 0  1 -1]

            ::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3,4])
                sage: B = F.basis()
                sage: elements = [B[1]-17*B[2]+6*B[3], B[1]-17*B[2]+B[4]]
                sage: F.echelon_form(elements)
                [B[1] - 17*B[2] + B[4], 6*B[3] - B[4]]

            ::

                sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
                sage: a,b,c = F.basis()
                sage: F.echelon_form([8*a+b+10*c, -3*a+b-c, a-b-c])
                [B['a'] + B['c'], B['b'] + 2*B['c']]
            """
            from sage.matrix.constructor import matrix
            mat = matrix([g._vector_() for g in elements])
            mat.echelonize()
            return [self.from_vector(vec) for vec in mat if vec]

        def submodule(self, gens,
                      check=True, already_echelonized=False, category=None):
            r"""
            The submodule spanned by a finite set of elements.

            INPUT:

            - ``gens`` -- a list or family of elements of ``self``

            - ``check`` -- (default: ``True``) whether to verify that the
               elements of ``gens`` are in ``self``.

            - ``already_echelonized`` -- (default: ``False``) whether
               the elements of ``gens`` are already in (not necessarily
               reduced) echelon form.

            If ``already_echelonized`` is ``False``, then the
            generators are put in reduced echelon form using
            :meth:`echelonize`, and reindexed by `0,1,...`.

            .. WARNING::

                At this point, this method only works for finite
                dimensional submodules and if matrices can be
                echelonized over the base ring.

            The basis of the submodule uses the same index set as the
            generators, and the lifting map sends `y_i` to `gens[i]`.


            .. SEEALSO::

                 - :meth:`ModulesWithBasis.FiniteDimensional.ParentMethods.quotient_module`
                 - :class:`sage.modules.with_basis.subquotient.SubmoduleWithBasis`

            EXAMPLES:

            We construct a submodule of the free `\QQ`-module generated by
            `x_0, x_1, x_2`. The submodule is spanned by `y_0 = x_0 - x_1` and
            `y_1 - x_1 - x_2`, and its basis elements are indexed by `0` and `1`::

                sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x")
                sage: x = X.basis()
                sage: gens = [x[0] - x[1], x[1] - x[2]]; gens
                [x[0] - x[1], x[1] - x[2]]
                sage: Y = X.submodule(gens, already_echelonized=True)
                sage: Y.print_options(prefix='y'); Y
                Free module generated by {0, 1} over Rational Field
                sage: y = Y.basis()
                sage: y[1]
                y[1]
                sage: y[1].lift()
                x[1] - x[2]
                sage: Y.retract(x[0]-x[2])
                y[0] + y[1]
                sage: Y.retract(x[0])
                Traceback (most recent call last):
                ...
                ValueError: x[0] is not in the image

            By using a family to specify a basis of the submodule, we obtain a
            submodule whose index set coincides with the index set of the family::

                sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x")
                sage: x = X.basis()
                sage: gens = Family({1 : x[0] - x[1], 3: x[1] - x[2]}); gens
                Finite family {1: x[0] - x[1], 3: x[1] - x[2]}
                sage: Y = X.submodule(gens, already_echelonized=True)
                sage: Y.print_options(prefix='y'); Y
                Free module generated by {1, 3} over Rational Field
                sage: y = Y.basis()
                sage: y[1]
                y[1]
                sage: y[1].lift()
                x[0] - x[1]
                sage: y[3].lift()
                x[1] - x[2]
                sage: Y.retract(x[0]-x[2])
                y[1] + y[3]
                sage: Y.retract(x[0])
                Traceback (most recent call last):
                ...
                ValueError: x[0] is not in the image

            It is not necessary that the generators of the submodule form
            a basis (an explicit basis will be computed)::

                sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x")
                sage: x = X.basis()
                sage: gens = [x[0] - x[1], 2*x[1] - 2*x[2], x[0] - x[2]]; gens
                [x[0] - x[1], 2*x[1] - 2*x[2], x[0] - x[2]]
                sage: Y = X.submodule(gens, already_echelonized=False)
                sage: Y.print_options(prefix='y')
                sage: Y
                Free module generated by {0, 1} over Rational Field
                sage: [b.lift() for b in Y.basis()]
                [x[0] - x[2], x[1] - x[2]]

            We now implement by hand the center of the algebra of the
            symmetric group `S_3`::

                sage: S3 = SymmetricGroup(3)
                sage: S3A = S3.algebra(QQ)
                sage: basis = S3A.annihilator_basis(S3A.algebra_generators(), S3A.bracket)
                sage: basis
                [(), (2,3) + (1,2) + (1,3), (1,2,3) + (1,3,2)]
                sage: center = S3A.submodule(basis,
                ....:                        category=AlgebrasWithBasis(QQ).Subobjects(),
                ....:                        already_echelonized=True)
                sage: center
                Free module generated by {0, 1, 2} over Rational Field
                sage: center in Algebras
                True
                sage: center.print_options(prefix='c')
                sage: c = center.basis()
                sage: c[1].lift()
                (2,3) + (1,2) + (1,3)
                sage: c[0]^2
                c[0]
                sage: e = 1/6*(c[0]+c[1]+c[2])
                sage: e.is_idempotent()
                True

            Of course, this center is best constructed using::

                sage: center = S3A.center()

            TESTS::

                sage: TestSuite(Y).run()
                sage: TestSuite(center).run()
            """
            if not already_echelonized:
                gens = self.echelon_form(gens)
            from sage.modules.with_basis.subquotient import SubmoduleWithBasis
            return SubmoduleWithBasis(gens, ambient=self, category=category)

        def tensor(*parents):
            """
            Return the tensor product of the parents.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example(); A.rename("A")
                sage: A.tensor(A,A)
                A # A # A
                sage: A.rename(None)
            """
            return parents[0].__class__.Tensor(parents, category = tensor.category_from_parents(parents))


    class ElementMethods:
        # TODO: Define the appropriate element methods here (instead of in
        # subclasses).  These methods should be consistent with those on
        # polynomials.

#         def _neg_(self):
#             """
#             Default implementation of negation by trying to multiply by -1.
#             TODO: doctest
#             """
#             return self._lmul_(-self.parent().base_ring().one(), self)


        def support_of_term(self):
            """
            Return the support of ``self``, where ``self`` is a monomial
            (possibly with coefficient).

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1,2,3,4]); X.rename("X")
                sage: X.monomial(2).support_of_term()
                2
                sage: X.term(3, 2).support_of_term()
                3

            An exception is raised if ``self`` has more than one term::

                sage: (X.monomial(2) + X.monomial(3)).support_of_term()
                Traceback (most recent call last):
                ...
                ValueError: B[2] + B[3] is not a single term
            """
            if len(self) == 1:
                return self.support()[0]
            else:
                raise ValueError("%s is not a single term"%(self))

        def leading_support(self, cmp=None):
            r"""
            Return the maximal element of the support of ``self``. Note
            that this may not be the term which actually appears first when
            ``self`` is printed.

            If the default ordering of the basis elements is not what is
            desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, ``0`` if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + 4*X.monomial(3)
                sage: x.leading_support()
                3
                sage: def cmp(x,y): return y-x
                sage: x.leading_support(cmp=cmp)
                1

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_support()
                [3]
            """
            return max_cmp(self.support(), cmp)


        def leading_item(self, cmp=None):
            r"""
            Return the pair ``(k, c)`` where

            .. MATH::

                c \cdot (\mbox{the basis element indexed by } k)

            is the leading term of ``self``.

            Here 'leading term' means that the corresponding basis element is
            maximal.  Note that this may not be the term which actually appears
            first when ``self`` is printed.  If the default term ordering is not
            what is desired, a comparison function, ``cmp(x,y)``, can be
            provided.  This should return a negative value if ``x < y``, ``0``
            if ``x == y`` and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + 4*X.monomial(3)
                sage: x.leading_item()
                (3, 4)
                sage: def cmp(x,y): return y-x
                sage: x.leading_item(cmp=cmp)
                (1, 3)

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_item()
                ([3], -5)
            """
            k = self.leading_support(cmp=cmp)
            return k, self[k]

        def leading_monomial(self, cmp=None):
            r"""
            Return the leading monomial of ``self``.

            This is the monomial whose corresponding basis element is
            maximal. Note that this may not be the term which actually appears
            first when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, ``0`` if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.leading_monomial()
                B[3]
                sage: def cmp(x,y): return y-x
                sage: x.leading_monomial(cmp=cmp)
                B[1]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_monomial()
                s[3]
            """
            return self.parent().monomial( self.leading_support(cmp=cmp) )

        def leading_coefficient(self, cmp=None):
            r"""
            Returns the leading coefficient of ``self``.

            This is the coefficient of the term whose corresponding basis element is
            maximal. Note that this may not be the term which actually appears
            first when ``self`` is printed.  If the default term ordering is not
            what is desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X")
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.leading_coefficient()
                1
                sage: def cmp(x,y): return y-x
                sage: x.leading_coefficient(cmp=cmp)
                3

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_coefficient()
                -5
            """
            return self.leading_item(cmp=cmp)[1]

        def leading_term(self, cmp=None):
            r"""
            Return the leading term of ``self``.

            This is the term whose corresponding basis element is
            maximal. Note that this may not be the term which actually appears
            first when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.leading_term()
                B[3]
                sage: def cmp(x,y): return y-x
                sage: x.leading_term(cmp=cmp)
                3*B[1]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_term()
                -5*s[3]
            """
            return self.parent().term(*self.leading_item(cmp=cmp))

        def trailing_support(self, cmp=None):
            r"""
            Return the minimal element of the support of ``self``. Note
            that this may not be the term which actually appears last when
            ``self`` is printed.

            If the default ordering of the basis elements is not what is
            desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, `0` if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + 4*X.monomial(3)
                sage: x.trailing_support()
                1
                sage: def cmp(x,y): return y-x
                sage: x.trailing_support(cmp=cmp)
                3

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_support()
                [1]
            """
            return min_cmp(self.support(), cmp)

        def trailing_item(self, cmp=None):
            r"""
            Returns the pair ``(c, k)`` where ``c*self.parent().monomial(k)``
            is the trailing term of ``self``.

            This is the monomial whose corresponding basis element is
            minimal. Note that this may not be the term which actually appears
            last when ``self`` is printed.  If the default term ordering is not
            what is desired, a comparison function ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.trailing_item()
                (1, 3)
                sage: def cmp(x,y): return y-x
                sage: x.trailing_item(cmp=cmp)
                (3, 1)

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_item()
                ([1], 2)
            """
            k = self.trailing_support(cmp=cmp)
            return k, self[k]

        def trailing_monomial(self, cmp=None):
            r"""
            Return the trailing monomial of ``self``.

            This is the monomial whose corresponding basis element is
            minimal. Note that this may not be the term which actually appears
            last when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.trailing_monomial()
                B[1]
                sage: def cmp(x,y): return y-x
                sage: x.trailing_monomial(cmp=cmp)
                B[3]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_monomial()
                s[1]
            """
            return self.parent().monomial( self.trailing_support(cmp=cmp) )

        def trailing_coefficient(self, cmp=None):
            r"""
            Return the trailing coefficient of ``self``.

            This is the coefficient of the monomial whose corresponding basis element is
            minimal. Note that this may not be the term which actually appears
            last when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.trailing_coefficient()
                3
                sage: def cmp(x,y): return y-x
                sage: x.trailing_coefficient(cmp=cmp)
                1

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_coefficient()
                2
            """

            return self.trailing_item(cmp=cmp)[1]

        def trailing_term(self, cmp=None):
            r"""
            Return the trailing term of ``self``.

            This is the term whose corresponding basis element is
            minimal. Note that this may not be the term which actually appears
            last when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.trailing_term()
                3*B[1]
                sage: def cmp(x,y): return y-x
                sage: x.trailing_term(cmp=cmp)
                B[3]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_term()
                2*s[1]
            """
            return self.parent().term( *self.trailing_item(cmp=cmp) )

        def map_coefficients(self, f):
            """
            Mapping a function on coefficients.

            INPUT:

            - ``f`` -- an endofunction on the coefficient ring of the
              free module

            Return a new element of ``self.parent()`` obtained by applying the
            function ``f`` to all of the coefficients of ``self``.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
                sage: B = F.basis()
                sage: f = B['a'] - 3*B['c']
                sage: f.map_coefficients(lambda x: x+5)
                6*B['a'] + 2*B['c']

            Killed coefficients are handled properly::

                sage: f.map_coefficients(lambda x: 0)
                0
                sage: list(f.map_coefficients(lambda x: 0))
                []

            ::

                sage: s = SymmetricFunctions(QQ).schur()
                sage: a = s([2,1])+2*s([3,2])
                sage: a.map_coefficients(lambda x: x*2)
                2*s[2, 1] + 4*s[3, 2]
            """
            return self.parent().sum_of_terms( (m, f(c)) for m,c in self )

        def map_support(self, f):
            """
            Mapping a function on the support.

            INPUT:

            - ``f`` -- an endofunction on the indices of the free module

            Return a new element of ``self.parent()`` obtained by
            applying the function ``f`` to all of the objects indexing
            the basis elements.

            EXAMPLES::

                sage: B = CombinatorialFreeModule(ZZ, [-1, 0, 1])
                sage: x = B.an_element(); x
                2*B[-1] + 2*B[0] + 3*B[1]
                sage: x.map_support(lambda i: -i)
                3*B[-1] + 2*B[0] + 2*B[1]

            ``f`` needs not be injective::

                sage: x.map_support(lambda i: 1)
                7*B[1]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: a = s([2,1])+2*s([3,2])
                sage: a.map_support(lambda x: x.conjugate())
                s[2, 1] + 2*s[2, 2, 1]

            TESTS::

                sage: B.zero()      # This actually failed at some point!!! See #8890
                0

                sage: y = B.zero().map_support(lambda i: i/0); y
                0
                sage: y.parent() is B
                True
            """
            return self.parent().sum_of_terms( (f(m), c) for m,c in self )

        def map_support_skip_none(self, f):
            """
            Mapping a function on the support.

            INPUT:

            - ``f`` -- an endofunction on the indices of the free module

            Returns a new element of ``self.parent()`` obtained by
            applying the function `f` to all of the objects indexing
            the basis elements.

            EXAMPLES::

                sage: B = CombinatorialFreeModule(ZZ, [-1, 0, 1])
                sage: x = B.an_element(); x
                2*B[-1] + 2*B[0] + 3*B[1]
                sage: x.map_support_skip_none(lambda i: -i if i else None)
                3*B[-1] + 2*B[1]

            ``f`` needs not be injective::

                sage: x.map_support_skip_none(lambda i: 1 if i else None)
                5*B[1]

            TESTS::

                sage: y = x.map_support_skip_none(lambda i: None); y
                0
                sage: y.parent() is B
                True
            """
            return self.parent().sum_of_terms( (fm,c) for (fm,c) in ((f(m), c) for m,c in self) if fm is not None)

        def map_item(self, f):
            """
            Mapping a function on items.

            INPUT:

            - ``f`` -- a function mapping pairs ``(index, coeff)`` to
              other such pairs

            Return a new element of ``self.parent()`` obtained by
            applying the function `f` to all items ``(index, coeff)`` of
            ``self``.

            EXAMPLES::

                sage: B = CombinatorialFreeModule(ZZ, [-1, 0, 1])
                sage: x = B.an_element(); x
                2*B[-1] + 2*B[0] + 3*B[1]
                sage: x.map_item(lambda i, c: (-i, 2*c))
                6*B[-1] + 4*B[0] + 4*B[1]

            ``f`` needs not be injective::

                sage: x.map_item(lambda i, c: (1, 2*c))
                14*B[1]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = lambda m,c: (m.conjugate(), 2*c)
                sage: a = s([2,1]) + s([1,1,1])
                sage: a.map_item(f)
                2*s[2, 1] + 2*s[3]
            """
            return self.parent().sum_of_terms( f(m,c) for m,c in self )

        def tensor(*elements):
            """
            Return the tensor product of its arguments, as an element of
            the tensor product of the parents of those elements.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example()
                sage: (a,b,c) = A.algebra_generators()
                sage: a.tensor(b, c)
                B[word: a] # B[word: b] # B[word: c]

            FIXME: is this a policy that we want to enforce on all parents?
            """
            assert(all(isinstance(element, Element) for element in elements))
            parents = [parent(element) for element in elements]
            return tensor(parents)._tensor_of_elements(elements) # good name???

    class Homsets(HomsetsCategory):
        class ParentMethods:
            def __call_on_basis__(self, **options):
                """
                Construct a morphism in this homset from a function defined
                on the basis.

                INPUT:

                - ``on_basis`` -- a function from the indices of the
                  basis of the domain of ``self`` to the codomain of
                  ``self``

                This method simply delegates the work to
                :meth:`ModulesWithBasis.ParentMethods.module_morphism`. It
                is used by :meth:`Homset.__call__` to handle the
                ``on_basis`` argument, and will disapear as soon as
                the logic will be generalized.

                EXAMPLES::

                    sage: X = CombinatorialFreeModule(QQ, [1,2,3]);   X.rename("X")
                    sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4]); Y.rename("Y")
                    sage: H = Hom(X, Y)
                    sage: x = X.basis()

                    sage: phi = H(on_basis = lambda i: Y.monomial(i) + 2*Y.monomial(i+1)) # indirect doctest
                    sage: phi
                    Generic morphism:
                    From: X
                    To:   Y
                    sage: phi(x[1] + x[3])
                    B[1] + 2*B[2] + B[3] + 2*B[4]

                Diagonal functions can be constructed using the ``diagonal`` option::

                    sage: X = CombinatorialFreeModule(QQ, [1,2,3,4]); X.rename("X")
                    sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4], key="Y"); Y.rename("Y")
                    sage: H = Hom(X, Y)
                    sage: x = X.basis()
                    sage: phi = H(diagonal = lambda x: x^2)
                    sage: phi(x[1] + x[2] + x[3])
                    B[1] + 4*B[2] + 9*B[3]

                TESTS::

                As for usual homsets, the argument can be a Python function::

                    sage: phi = H(lambda x: Y.zero())
                    sage: phi
                    Generic morphism:
                      From: X
                      To:   Y
                    sage: phi(x[1] + x[3])
                    0

               We check that the homset category is properly set up::

                    sage: X = CombinatorialFreeModule(QQ, [1,2,3]);   X.rename("X")
                    sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4]); Y.rename("Y")
                    sage: H = Hom(X, Y)
                    sage: H.zero().category_for()
                    Category of finite dimensional vector spaces with basis over Rational Field
                """
                return self.domain().module_morphism(codomain = self.codomain(),
                                                     **options)

    class MorphismMethods:
        @cached_method
        def on_basis(self):
            """
            Return the action of this morphism on basis elements.

            OUTPUT:

            - a function from the indices of the basis of the domain to
              the codomain

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1,2,3]);   X.rename("X")
                sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4]); Y.rename("Y")
                sage: H = Hom(X, Y)
                sage: x = X.basis()

                sage: f = H(lambda x: Y.zero()).on_basis()
                sage: f(2)
                0

                sage: f = lambda i: Y.monomial(i) + 2*Y.monomial(i+1)
                sage: g = H(on_basis = f).on_basis()
                sage: g(2)
                B[2] + 2*B[3]
                sage: g == f
                True
            """
            return self._on_basis

        def _on_basis(self, i):
            """
            Return the image of ``self`` on the basis element indexed by ``i``.

            INPUT:

            - ``i`` -- the index of an element of the basis of the domain of ``self``

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1,2,3]); X.rename("X")
                sage: phi = End(X)(lambda x: 2*x)
                sage: phi._on_basis(3)
                2*B[3]
            """
            return self(self.domain().monomial(i))

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of modules with basis constructed by cartesian products
        of modules with basis.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: ModulesWithBasis(QQ).CartesianProducts().extra_super_categories()
                [Category of vector spaces with basis over Rational Field]
                sage: ModulesWithBasis(QQ).CartesianProducts().super_categories()
                [Category of Cartesian products of modules with basis over Rational Field,
                 Category of vector spaces with basis over Rational Field,
                 Category of Cartesian products of vector spaces over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:

            def _an_element_(self):
                """
                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: B = HopfAlgebrasWithBasis(QQ).example(); B
                    An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                    sage: A.an_element(), B.an_element()
                    (2*B[word: ] + 2*B[word: a] + 3*B[word: b], B[()] + 4*B[(1,2,3)] + 2*B[(1,3)])
                    sage: cartesian_product((A, B, A)).an_element()           # indirect doctest
                    2*B[(0, word: )] + 2*B[(0, word: a)] + 3*B[(0, word: b)]
                """
                from cartesian_product import cartesian_product
                return cartesian_product([module.an_element() for module in self.modules])

    class TensorProducts(TensorProductsCategory):
        """
        The category of modules with basis constructed by tensor product of
        modules with basis.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: ModulesWithBasis(QQ).TensorProducts().extra_super_categories()
                [Category of vector spaces with basis over Rational Field]
                sage: ModulesWithBasis(QQ).TensorProducts().super_categories()
                [Category of tensor products of modules with basis over Rational Field,
                 Category of vector spaces with basis over Rational Field,
                 Category of tensor products of vector spaces over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            """
            Implements operations on tensor products of modules with basis.
            """
            pass

        class ElementMethods:
            """
            Implements operations on elements of tensor products of modules
            with basis.
            """

            def apply_multilinear_morphism(self, f, codomain = None):
                r"""
                Return the result of applying the morphism induced by ``f``
                to ``self``.

                INPUT:

                - ``f`` -- a multilinear morphism from the component
                  modules of the parent tensor product to any module

                - ``codomain`` -- the codomain of ``f`` (optional)

                By the universal property of the tensor product, ``f``
                induces a linear morphism from `self.parent()` to the
                target module. Returns the result of applying that
                morphism to ``self``.

                The codomain is used for optimizations purposes
                only. If it's not provided, it's recovered by calling
                ``f`` on the zero input.

                EXAMPLES:

                We start with simple (admittedly not so interesting)
                examples, with two modules `A` and `B`::

                    sage: A = CombinatorialFreeModule(ZZ, [1,2], prefix="A"); A.rename("A")
                    sage: B = CombinatorialFreeModule(ZZ, [3,4], prefix="B"); B.rename("B")

                and `f` the bilinear morphism `(a,b) \mapsto b \otimes a`
                from `A \times B` to `B \otimes A`::

                    sage: def f(a,b):
                    ....:     return tensor([b,a])

                Now, calling applying `f` on `a \otimes b` returns the same
                as `f(a,b)`::

                    sage: a = A.monomial(1) + 2 * A.monomial(2); a
                    A[1] + 2*A[2]
                    sage: b = B.monomial(3) - 2 * B.monomial(4); b
                    B[3] - 2*B[4]
                    sage: f(a,b)
                    B[3] # A[1] + 2*B[3] # A[2] - 2*B[4] # A[1] - 4*B[4] # A[2]
                    sage: tensor([a,b]).apply_multilinear_morphism(f)
                    B[3] # A[1] + 2*B[3] # A[2] - 2*B[4] # A[1] - 4*B[4] # A[2]

                `f` may be a bilinear morphism to any module over the
                base ring of `A` and `B`. Here the codomain is `\ZZ`::

                    sage: def f(a,b):
                    ....:     return sum(a.coefficients(), 0) * sum(b.coefficients(), 0)
                    sage: f(a,b)
                    -3
                    sage: tensor([a,b]).apply_multilinear_morphism(f)
                    -3

                Mind the `0` in the sums above; otherwise `f` would
                not return `0` in `\ZZ`::

                    sage: def f(a,b):
                    ....:     return sum(a.coefficients()) * sum(b.coefficients())
                    sage: type(f(A.zero(), B.zero()))
                    <type 'int'>

                Which would be wrong and break this method::

                    sage: tensor([a,b]).apply_multilinear_morphism(f)
                    Traceback (most recent call last):
                    ...
                    AttributeError: 'int' object has no attribute 'parent'

                Here we consider an example where the codomain is a
                module with basis with a different base ring::

                    sage: C = CombinatorialFreeModule(QQ, [(1,3),(2,4)], prefix="C"); C.rename("C")
                    sage: def f(a,b):
                    ....:     return C.sum_of_terms( [((1,3), QQ(a[1]*b[3])), ((2,4), QQ(a[2]*b[4]))] )
                    sage: f(a,b)
                    C[(1, 3)] - 4*C[(2, 4)]
                    sage: tensor([a,b]).apply_multilinear_morphism(f)
                    C[(1, 3)] - 4*C[(2, 4)]

                 We conclude with a real life application, where we
                 check that the antipode of the Hopf algebra of
                 Symmetric functions on the Schur basis satisfies its
                 defining formula::

                    sage: Sym = SymmetricFunctions(QQ)
                    sage: s = Sym.schur()
                    sage: def f(a,b): return a*b.antipode()
                    sage: x = 4*s.an_element(); x
                    8*s[] + 8*s[1] + 12*s[2]
                    sage: x.coproduct().apply_multilinear_morphism(f)
                    8*s[]
                    sage: x.coproduct().apply_multilinear_morphism(f) == x.counit()
                    True

                We recover the constant term of `x`, as desired.

                .. TODO::

                    Extract a method to linearize a multilinear
                    morphism, and delegate the work there.
                """
                K = self.parent().base_ring()
                modules = self.parent()._sets
                if codomain is None:
                    try:
                        codomain = f.codomain()
                    except AttributeError:
                        codomain = f(*[module.zero() for module in modules]).parent()
                if codomain in ModulesWithBasis(K):
                    return codomain.linear_combination((f(*[module.monomial(t) for (module,t) in zip(modules, m)]), c)
                                                       for m,c in self)
                else:
                    return sum((c * f(*[module.monomial(t) for (module,t) in zip(modules, m)])
                                for m,c in self),
                               codomain.zero())

    class DualObjects(DualObjectsCategory):

        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: ModulesWithBasis(ZZ).DualObjects().extra_super_categories()
                [Category of modules over Integer Ring]
                sage: ModulesWithBasis(QQ).DualObjects().super_categories()
                [Category of duals of vector spaces over Rational Field, Category of duals of modules with basis over Rational Field]
            """
            return [Modules(self.base_category().base_ring())]
