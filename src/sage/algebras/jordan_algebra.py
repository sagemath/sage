# sage.doctest: needs sage.combinat sage.modules
r"""
Jordan Algebras

AUTHORS:

- Travis Scrimshaw (2014-04-02): initial version
- Travis Scrimshaw (2023-05-09): added the 27 dimensional exceptional
  Jordan algebra
"""

#*****************************************************************************
#  Copyright (C) 2014, 2023 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import AlgebraElement
from sage.structure.richcmp import richcmp
from sage.categories.magmatic_algebras import MagmaticAlgebras
from sage.misc.cachefunc import cached_method
from sage.structure.element import Matrix
from sage.modules.free_module import FreeModule
from sage.matrix.constructor import matrix
from sage.sets.family import Family


class JordanAlgebra(UniqueRepresentation, Parent):
    r"""
    A Jordan algebra.

    A *Jordan algebra* is a magmatic algebra (over a commutative ring
    `R`) whose multiplication satisfies the following axioms:

    - `xy = yx`, and
    - `(xy)(xx) = x(y(xx))` (the Jordan identity).

    See [Ja1971]_, [Ch2012]_, and [McC1978]_, for example.

    These axioms imply that a Jordan algebra is power-associative and the
    following generalization of Jordan's identity holds [Al1947]_:
    `(x^m y) x^n = x^m (y x^n)` for all `m, n \in \ZZ_{>0}`.

    Let `A` be an associative algebra over a ring `R` in which `2` is
    invertible. We construct a Jordan algebra `A^+` with ground set `A`
    by defining the multiplication as

    .. MATH::

        x \circ y = \frac{xy + yx}{2}.

    Often the multiplication is written as `x \circ y` to avoid confusion
    with the product in the associative algebra `A`. We note that if `A` is
    commutative then this reduces to the usual multiplication in `A`.

    Jordan algebras constructed in this fashion, or their subalgebras,
    are called *special*. All other Jordan algebras are called *exceptional*.

    Jordan algebras can also be constructed from a module `M` over `R` with
    a symmetric bilinear form `(\cdot, \cdot) : M \times M \to R`.
    We begin with the module `M^* = R \oplus M` and define multiplication
    in `M^*` by

    .. MATH::

        (\alpha + x) \circ (\beta + y) =
        \underbrace{\alpha \beta + (x,y)}_{\in R}
        + \underbrace{\beta x + \alpha y}_{\in M},

    where `\alpha, \beta \in R` and `x,y \in M`.

    INPUT:

    Can be either an associative algebra `A` or a symmetric bilinear
    form given as a matrix (possibly followed by, or preceded by, a base
    ring argument).

    EXAMPLES:

    We let the base algebra `A` be the free algebra on 3 generators::

        sage: F.<x,y,z> = FreeAlgebra(QQ)
        sage: J = JordanAlgebra(F); J
        Jordan algebra of Free Algebra on 3 generators (x, y, z) over Rational Field
        sage: a,b,c = map(J, F.gens())
        sage: a*b
        1/2*x*y + 1/2*y*x
        sage: b*a
        1/2*x*y + 1/2*y*x

    Jordan algebras are typically non-associative::

        sage: (a*b)*c
        1/4*x*y*z + 1/4*y*x*z + 1/4*z*x*y + 1/4*z*y*x
        sage: a*(b*c)
        1/4*x*y*z + 1/4*x*z*y + 1/4*y*z*x + 1/4*z*y*x

    We check the Jordan identity::

        sage: (a*b)*(a*a) == a*(b*(a*a))
        True
        sage: x = a + c
        sage: y = b - 2*a
        sage: (x*y)*(x*x) == x*(y*(x*x))
        True

    Next we construct a Jordan algebra from a symmetric bilinear form::

        sage: m = matrix([[-2,3],[3,4]])
        sage: J.<a,b,c> = JordanAlgebra(m); J
        Jordan algebra over Integer Ring given by the symmetric bilinear form:
        [-2  3]
        [ 3  4]
        sage: a
        1 + (0, 0)
        sage: b
        0 + (1, 0)
        sage: x = 3*a - 2*b + c; x
        3 + (-2, 1)

    We again show that Jordan algebras are usually non-associative::

        sage: (x*b)*b
        -6 + (7, 0)
        sage: x*(b*b)
        -6 + (4, -2)

    We verify the Jordan identity::

        sage: y = -a + 4*b - c
        sage: (x*y)*(x*x) == x*(y*(x*x))
        True

    The base ring, while normally inferred from the matrix, can also
    be explicitly specified::

        sage: J.<a,b,c> = JordanAlgebra(m, QQ); J
        Jordan algebra over Rational Field given by the symmetric bilinear form:
        [-2  3]
        [ 3  4]
        sage: J.<a,b,c> = JordanAlgebra(QQ, m); J # either order work
        Jordan algebra over Rational Field given by the symmetric bilinear form:
        [-2  3]
        [ 3  4]

    REFERENCES:

    - :wikipedia:`Jordan_algebra`
    - [Ja1971]_
    - [Ch2012]_
    - [McC1978]_
    - [Al1947]_
    """
    @staticmethod
    def __classcall_private__(self, arg0, arg1=None, names=None):
        """
        Choose the correct parent based upon input.

        TESTS:

        We check arguments with passing in an associative algebra::

            sage: cat = Algebras(QQ).WithBasis().FiniteDimensional()
            sage: C = CombinatorialFreeModule(QQ, ['x','y','z'], category=cat)
            sage: J1 = JordanAlgebra(C, names=['a','b','c'])
            sage: J2.<a,b,c> = JordanAlgebra(C)
            sage: J1 is J2
            True

        We check with passing in a symmetric bilinear form::

            sage: m = matrix([[0,1],[1,1]])
            sage: J1 = JordanAlgebra(m)
            sage: J2 = JordanAlgebra(QQ, m)
            sage: J3 = JordanAlgebra(m, QQ)
            sage: J1 is J2
            False
            sage: J2 is J3
            True
            sage: J4 = JordanAlgebra(ZZ, m)
            sage: J1 is J4
            True
            sage: m = matrix(QQ, [[0,1],[1,1]])
            sage: J1 = JordanAlgebra(m)
            sage: J1 is J2
            True
        """
        if names is not None:
            if isinstance(names, str):
                names = names.split(',')
            names = tuple(names)

        if arg1 is None:
            if not isinstance(arg0, Matrix):
                from sage.algebras.octonion_algebra import OctonionAlgebra
                if isinstance(arg0, OctonionAlgebra):
                    return ExceptionalJordanAlgebra(arg0)
                if arg0.base_ring().characteristic() == 2:
                    raise ValueError("the base ring cannot have characteristic 2")
                return SpecialJordanAlgebra(arg0, names)
            arg0, arg1 = arg0.base_ring(), arg0
        elif isinstance(arg0, Matrix):
            arg0, arg1 = arg1, arg0

        # arg0 is the base ring and arg1 is a matrix
        if not arg1.is_symmetric():
            raise ValueError("the bilinear form is not symmetric")

        arg1 = arg1.change_ring(arg0) # This makes a copy
        arg1.set_immutable()
        return JordanAlgebraSymmetricBilinear(arg0, arg1, names=names)

    def _test_jordan_relations(self, **options):
        r"""
        Test the Jordan algebra relations.

        The Jordan algebra relations are

        - `xy = yx`, and
        - `(xy)(xx) = x(y(xx))` (the Jordan identity).

        EXAMPLES::

            sage: O = OctonionAlgebra(GF(7), 1, 3, 4)
            sage: J = JordanAlgebra(O)
            sage: J._test_jordan_relations()
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        from sage.misc.misc import some_tuples
        for x, y in some_tuples(S, 2, tester._max_runs):
            tester.assertEqual(x * y, y * x)
            tester.assertEqual((x * y) * (x * x), x * (y * (x * x)))


class SpecialJordanAlgebra(JordanAlgebra):
    r"""
    A (special) Jordan algebra `A^+` from an associative algebra `A`.
    """
    def __init__(self, A, names=None):
        """
        Initialize ``self``.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: TestSuite(J).run()
            sage: J.category()
            Category of commutative unital algebras with basis over Rational Field
        """
        R = A.base_ring()
        C = MagmaticAlgebras(R)
        if A not in C.Associative():
            raise ValueError("A is not an associative algebra")

        self._A = A
        cat = C.Commutative()
        if A in C.Unital():
            cat = cat.Unital()
        if A in C.WithBasis():
            cat = cat.WithBasis()
        if A in C.FiniteDimensional():
            cat = cat.FiniteDimensional()

        Parent.__init__(self, base=R, names=names, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: JordanAlgebra(F)
            Jordan algebra of Free Algebra on 3 generators (x, y, z) over Rational Field
        """
        return "Jordan algebra of {}".format(self._A)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J(5)
            5
            sage: elt = J(x + 2*x*y); elt
            x + 2*x*y
            sage: elt.parent() is J
            True
        """
        return self.element_class(self, self._A(x))

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J.an_element()
            2 + 2*x + 3*y
        """
        return self.element_class(self, self._A.an_element())

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J.basis()
            Lazy family (Term map(i))_{i in Free monoid on 3 generators (x, y, z)}
        """
        B = self._A.basis()
        return Family(B.keys(), lambda x: self.element_class(self, B[x]), name="Term map")

    algebra_generators = basis

    # TODO: Keep this until we can better handle R.<...> shorthand
    def gens(self) -> tuple:
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: cat = Algebras(QQ).WithBasis().FiniteDimensional()
            sage: C = CombinatorialFreeModule(QQ, ['x','y','z'], category=cat)
            sage: J = JordanAlgebra(C)
            sage: J.gens()
            (B['x'], B['y'], B['z'])

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J.gens()
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite set
        """
        return tuple(self.algebra_generators())

    @cached_method
    def zero(self):
        """
        Return the element `0`.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J.zero()
            0
        """
        return self.element_class(self, self._A.zero())

    @cached_method
    def one(self):
        """
        Return the element `1` if it exists.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J.one()
            1
        """
        return self.element_class(self, self._A.one())

    class Element(AlgebraElement):
        """
        An element of a special Jordan algebra.
        """
        def __init__(self, parent, x):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: TestSuite(a + 2*b - c).run()
            """
            self._x = x
            AlgebraElement.__init__(self, parent)

        def _repr_(self):
            """
            Return a string representation of ``self``.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: a + 2*b - c
                x + 2*y - z
            """
            return repr(self._x)

        def _latex_(self):
            """
            Return a latex representation of ``self``.

            EXAMPLES::

                sage: F.<x0,x1,x2> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: latex(a + 2*b - c)
                x_{0} + 2 x_{1} - x_{2}
            """
            from sage.misc.latex import latex
            return latex(self._x)

        def __bool__(self) -> bool:
            """
            Return if ``self`` is nonzero.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: bool(a + 2*b - c)
                True
            """
            return bool(self._x)

        def __eq__(self, other):
            """
            Check equality.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: elt = a + 2*b - c
                sage: elt == elt
                True
                sage: elt == x
                False
                sage: elt == 2*b
                False
            """
            if not isinstance(other, SpecialJordanAlgebra.Element):
                return False
            if other.parent() != self.parent():
                return False
            return self._x == other._x

        def __ne__(self, other):
            """
            Check inequality.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: elt = a + 2*b - c
                sage: elt != elt
                False
                sage: elt != x
                True
                sage: elt != 2*b
                True
            """
            return not self == other

        def _add_(self, other):
            """
            Add ``self`` and ``other``.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: a + 2*b
                x + 2*y
            """
            return self.__class__(self.parent(), self._x + other._x)

        def _neg_(self):
            """
            Negate ``self``.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: -(a + 2*b)
                -x - 2*y
            """
            return self.__class__(self.parent(), -self._x)

        def _sub_(self, other):
            """
            Subtract ``other`` from ``self``.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: a - 2*b
                x - 2*y
            """
            return self.__class__(self.parent(), self._x - other._x)

        def _mul_(self, other):
            """
            Multiply ``self`` and ``other``.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: (a + 2*b) * (c - b)
                -1/2*x*y + 1/2*x*z - 1/2*y*x - 2*y^2 + y*z + 1/2*z*x + z*y

                sage: F.<x,y,z> = FreeAlgebra(GF(3))
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: (a + 2*b) * (c - b)
                x*y + 2*x*z + y*x + y^2 + y*z + 2*z*x + z*y
            """
            x = self._x
            y = other._x
            # This is safer than dividing by 2
            R = self.parent().base_ring()
            return self.__class__(self.parent(), (x*y + y*x) * ~R(2))

        def _lmul_(self, other):
            """
            Multiply ``self`` by the scalar ``other`` on the left.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: (a + b) * 2
                2*x + 2*y
            """
            return self.__class__(self.parent(), self._x * other)

        def _rmul_(self, other):
            """
            Multiply ``self`` and the scalar ``other`` by the right
            action.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: 2 * (a + b)
                2*x + 2*y
            """
            return self.__class__(self.parent(), other * self._x)

        def monomial_coefficients(self, copy=True):
            """
            Return a dictionary whose keys are indices of basis elements in
            the support of ``self`` and whose values are the corresponding
            coefficients.

            INPUT:

            - ``copy`` -- boolean (default: ``True``); if ``self`` is
              internally represented by a dictionary ``d``, then make a copy of
              ``d``; if ``False``, then this can cause undesired behavior by
              mutating ``d``

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: elt = a + 2*b - c
                sage: elt.monomial_coefficients()
                {x: 1, y: 2, z: -1}
            """
            return self._x.monomial_coefficients(copy)


class JordanAlgebraSymmetricBilinear(JordanAlgebra):
    r"""
    A Jordan algebra given by a symmetric bilinear form `m`.
    """
    def __init__(self, R, form, names=None):
        """
        Initialize ``self``.

        TESTS::

            sage: m = matrix([[-2,3],[3,4]])
            sage: J = JordanAlgebra(m)
            sage: TestSuite(J).run()
        """
        self._form = form
        self._M = FreeModule(R, form.ncols())
        cat = MagmaticAlgebras(R).Commutative().Unital().FiniteDimensional().WithBasis()
        Parent.__init__(self, base=R, names=names, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: m = matrix([[-2,3],[3,4]])
            sage: JordanAlgebra(m)
            Jordan algebra over Integer Ring given by the symmetric bilinear form:
            [-2  3]
            [ 3  4]
        """
        return "Jordan algebra over {} given by the symmetric bilinear" \
               " form:\n{}".format(self.base_ring(), self._form)

    def _element_constructor_(self, *args):
        """
        Construct an element of ``self`` from ``s``.

        Here ``s`` can be a pair of an element of `R` and an
        element of `M`, or an element of `R`, or an element of
        `M`, or an element of a(nother) Jordan algebra given
        by a symmetric bilinear form.

        EXAMPLES::

            sage: m = matrix([[0,1],[1,1]])
            sage: J = JordanAlgebra(m)
            sage: J(2)
            2 + (0, 0)
            sage: J((-4, (2, 5)))
            -4 + (2, 5)
            sage: J((-4, (ZZ^2)((2, 5))))
            -4 + (2, 5)
            sage: J(2, (-2, 3))
            2 + (-2, 3)
            sage: J(2, (ZZ^2)((-2, 3)))
            2 + (-2, 3)
            sage: J(-1, 1, 0)
            -1 + (1, 0)
            sage: J((ZZ^2)((1, 3)))
            0 + (1, 3)

            sage: m = matrix([[2]])
            sage: J = JordanAlgebra(m)
            sage: J(2)
            2 + (0)
            sage: J((-4, (2,)))
            -4 + (2)
            sage: J(2, (-2,))
            2 + (-2)
            sage: J(-1, 1)
            -1 + (1)
            sage: J((ZZ^1)((3,)))
            0 + (3)

            sage: m = Matrix(QQ, [])
            sage: J = JordanAlgebra(m)
            sage: J(2)
            2 + ()
            sage: J((-4, ()))
            -4 + ()
            sage: J(2, ())
            2 + ()
            sage: J(-1)
            -1 + ()
            sage: J((ZZ^0)(()))
            0 + ()
        """
        R = self.base_ring()
        if len(args) == 1:
            s = args[0]

            if isinstance(s, JordanAlgebraSymmetricBilinear.Element):
                if s.parent() is self:
                    return s
                return self.element_class(self, R(s._s), self._M(s._v))

            if isinstance(s, (list, tuple)):
                if len(s) != 2:
                    raise ValueError("must be length 2")
                return self.element_class(self, R(s[0]), self._M(s[1]))

            if s in self._M:
                return self.element_class(self, R.zero(), self._M(s))

            return self.element_class(self, R(s), self._M.zero())

        if len(args) == 2 and (isinstance(args[1], (list, tuple)) or args[1] in self._M):
            return self.element_class(self, R(args[0]), self._M(args[1]))

        if len(args) == self._form.ncols() + 1:
            return self.element_class(self, R(args[0]), self._M(args[1:]))

        raise ValueError("unable to construct an element from the given data")

    def _coerce_map_from_base_ring(self):
        """
        Return a coercion map from the base ring of ``self``.

        TESTS::

            sage: J = JordanAlgebra(Matrix([[0, 1], [1, 1]]))
            sage: J._coerce_map_from_base_ring()
            Conversion map:
              From: Integer Ring
              To:   Jordan algebra over Integer Ring given by the symmetric bilinear form:
            [0 1]
            [1 1]
            sage: J.coerce_map_from(ZZ)
            Coercion map:
              From: Integer Ring
              To:   Jordan algebra over Integer Ring given by the symmetric bilinear form:
            [0 1]
            [1 1]
        """
        # Return a DefaultConvertMap_unique; this can pass additional
        # arguments to _element_constructor_, unlike the map returned
        # by UnitalAlgebras.ParentMethods._coerce_map_from_base_ring.
        return self._generic_coerce_map(self.base_ring())

    @cached_method
    def basis(self):
        """
        Return a basis of ``self``.

        The basis returned begins with the unity of `R` and continues with
        the standard basis of `M`.

        EXAMPLES::

            sage: m = matrix([[0,1],[1,1]])
            sage: J = JordanAlgebra(m)
            sage: J.basis()
            Family (1 + (0, 0), 0 + (1, 0), 0 + (0, 1))
        """
        R = self.base_ring()
        ret = (self.element_class(self, R.one(), self._M.zero()),)
        ret += tuple(self.element_class(self, R.zero(), x)
                     for x in self._M.basis())
        return Family(ret)

    algebra_generators = basis

    def gens(self) -> tuple:
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: m = matrix([[0,1],[1,1]])
            sage: J = JordanAlgebra(m)
            sage: J.gens()
            (1 + (0, 0), 0 + (1, 0), 0 + (0, 1))
        """
        return tuple(self.algebra_generators())

    @cached_method
    def zero(self):
        """
        Return the element 0.

        EXAMPLES::

            sage: m = matrix([[0,1],[1,1]])
            sage: J = JordanAlgebra(m)
            sage: J.zero()
            0 + (0, 0)
        """
        return self.element_class(self, self.base_ring().zero(), self._M.zero())

    @cached_method
    def one(self):
        """
        Return the element 1 if it exists.

        EXAMPLES::

            sage: m = matrix([[0,1],[1,1]])
            sage: J = JordanAlgebra(m)
            sage: J.one()
            1 + (0, 0)
        """
        return self.element_class(self, self.base_ring().one(), self._M.zero())

    class Element(AlgebraElement):
        """
        An element of a Jordan algebra defined by a symmetric bilinear form.
        """
        def __init__(self, parent, s, v):
            """
            Initialize ``self``.

            TESTS::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: TestSuite(a + 2*b - c).run()
            """
            self._s = s
            self._v = v
            AlgebraElement.__init__(self, parent)

        def _repr_(self):
            """
            Return a string representation of ``self``.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: a + 2*b - c
                1 + (2, -1)
            """
            return "{} + {}".format(self._s, self._v)

        def _latex_(self):
            r"""
            Return a latex representation of ``self``.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: latex(a + 2*b - c)
                1 + \left(2,\,-1\right)
            """
            from sage.misc.latex import latex
            return "{} + {}".format(latex(self._s), latex(self._v))

        def __bool__(self) -> bool:
            """
            Return if ``self`` is nonzero.

            TESTS::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: bool(1)
                True
                sage: bool(b)
                True
                sage: bool(a + 2*b - c)
                True
            """
            return bool(self._s) or bool(self._v)

        def __eq__(self, other):
            """
            Check equality.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: x = 4*a - b + 3*c
                sage: x == J((4, (-1, 3)))
                True
                sage: a == x
                False

                sage: m = matrix([[-2,3],[3,4]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: 4*a - b + 3*c == x
                False
            """
            if not isinstance(other, JordanAlgebraSymmetricBilinear.Element):
                return False
            if other.parent() != self.parent():
                return False
            return self._s == other._s and self._v == other._v

        def __ne__(self, other):
            """
            Check inequality.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: x = 4*a - b + 3*c
                sage: x != J((4, (-1, 3)))
                False
                sage: a != x
                True

                sage: m = matrix([[-2,3],[3,4]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: 4*a - b + 3*c != x
                True
            """
            return not self == other

        def _add_(self, other):
            """
            Add ``self`` and ``other``.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: a + b
                1 + (1, 0)
                sage: b + c
                0 + (1, 1)
            """
            return self.__class__(self.parent(), self._s + other._s, self._v + other._v)

        def _neg_(self):
            """
            Negate ``self``.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: -(a + b - 2*c)
                -1 + (-1, 2)
            """
            return self.__class__(self.parent(), -self._s, -self._v)

        def _sub_(self, other):
            """
            Subtract ``other`` from ``self``.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: a - b
                1 + (-1, 0)
                sage: b - c
                0 + (1, -1)
            """
            return self.__class__(self.parent(), self._s - other._s, self._v - other._v)

        def _mul_(self, other):
            """
            Multiply ``self`` and ``other``.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: (4*a - b + 3*c)*(2*a + 2*b - c)
                12 + (6, 2)

                sage: m = matrix([[-2,3],[3,4]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: (4*a - b + 3*c)*(2*a + 2*b - c)
                21 + (6, 2)
            """
            P = self.parent()
            return self.__class__(P,
                                  self._s * other._s
                                   + (self._v * P._form * other._v.column())[0],
                                  other._s * self._v + self._s * other._v)

        def _lmul_(self, other):
            """
            Multiply ``self`` by the scalar ``other`` on the left.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: (a + b - c) * 2
                2 + (2, -2)
            """
            return self.__class__(self.parent(), self._s * other, self._v * other)

        def _rmul_(self, other):
            """
            Multiply ``self`` with the scalar ``other`` by the right
            action.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: 2 * (a + b - c)
                2 + (2, -2)
            """
            return self.__class__(self.parent(), other * self._s, other * self._v)

        def monomial_coefficients(self, copy=True):
            """
            Return a dictionary whose keys are indices of basis elements in
            the support of ``self`` and whose values are the corresponding
            coefficients.

            INPUT:

            - ``copy`` -- ignored

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: elt = a + 2*b - c
                sage: elt.monomial_coefficients()
                {0: 1, 1: 2, 2: -1}
            """
            d = {0: self._s}
            for i,c in enumerate(self._v):
                d[i+1] = c
            return d

        def trace(self):
            r"""
            Return the trace of ``self``.

            The trace of an element `\alpha + x \in M^*` is given by
            `t(\alpha + x) = 2 \alpha`.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: x = 4*a - b + 3*c
                sage: x.trace()
                8
            """
            return 2 * self._s

        def norm(self):
            r"""
            Return the norm of ``self``.

            The norm of an element `\alpha + x \in M^*` is given by
            `n(\alpha + x) = \alpha^2 - (x, x)`.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: x = 4*a - b + 3*c; x
                4 + (-1, 3)
                sage: x.norm()
                13
            """
            return self._s * self._s - (self._v * self.parent()._form
                                        * self._v.column())[0]

        def bar(self):
            r"""
            Return the result of the bar involution of ``self``.

            The bar involution `\bar{\cdot}` is the `R`-linear
            endomorphism of `M^*` defined by `\bar{1} = 1` and
            `\bar{x} = -x` for `x \in M`.

            EXAMPLES::

                sage: m = matrix([[0,1],[1,1]])
                sage: J.<a,b,c> = JordanAlgebra(m)
                sage: x = 4*a - b + 3*c
                sage: x.bar()
                4 + (1, -3)

            We check that it is an algebra morphism::

                sage: y = 2*a + 2*b - c
                sage: x.bar() * y.bar() == (x*y).bar()
                True
            """
            return self.__class__(self.parent(), self._s, -self._v)


class ExceptionalJordanAlgebra(JordanAlgebra):
    r"""
    The exceptional `27` dimensional Jordan algebra as self-adjoint
    `3 \times 3` matrix over an octonion algebra.

    Let `\mathbf{O}` be the :class:`OctonionAlgebra` over a commutative
    ring `R` of characteristic not equal to `2`. The *exceptional Jordan
    algebra* `\mathfrak{h}_3(\mathbf{O})` is a `27` dimensional free
    `R`-module spanned by the matrices

    .. MATH::

        \begin{bmatrix}
        \alpha & x & y \\
        x^* & \beta & z \\
        y^* & z^* & \gamma
        \end{bmatrix}

    for `\alpha, \beta, \gamma \in R` and `x, y, z \in \mathbf{O}`,
    with multiplication given by the usual symmetrizer operation
    `X \circ Y = \frac{1}{2}(XY + YX)`.

    These are also known as *Albert algebras* due to the work of
    Abraham Adrian Albert on these algebras over `\RR`.

    EXAMPLES:

    We construct an exceptional Jordan algebra over `\QQ` and perform
    some basic computations::

        sage: O = OctonionAlgebra(QQ)
        sage: J = JordanAlgebra(O)
        sage: gens = J.gens()
        sage: gens[1]
        [0 0 0]
        [0 1 0]
        [0 0 0]
        sage: gens[3]
        [0 1 0]
        [1 0 0]
        [0 0 0]
        sage: gens[1] * gens[3]
        [  0 1/2   0]
        [1/2   0   0]
        [  0   0   0]

    The Lie algebra of derivations of the exceptional Jordan algebra
    is isomorphic to the simple Lie algebra of type `F_4`. We verify
    that we the derivation module has the correct dimension::

        sage: len(J.derivations_basis())  # long time
        52
        sage: LieAlgebra(QQ, cartan_type='F4').dimension()
        52

    REFERENCES:

    - :wikipedia:`Albert_algebra`
    - :wikipedia:`Jordan_algebra#Examples`
    - :wikipedia:`Hurwitz's_theorem_(composition_algebras)#Applications_to_Jordan_algebras`
    - `<https://math.ucr.edu/home/baez/octonions/octonions.pdf>`_
    """
    def __init__(self, O):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: J = JordanAlgebra(O)
            sage: TestSuite(J).run()  # long time

            sage: O = OctonionAlgebra(QQ, 1, -2, 9)
            sage: J = JordanAlgebra(O)
            sage: TestSuite(J).run()  # long time

            sage: R.<x, y> = GF(11)[]
            sage: O = OctonionAlgebra(R, 1, x + y, 9)
            sage: J = JordanAlgebra(O)
            sage: TestSuite(J).run()  # long time

            sage: O = OctonionAlgebra(ZZ)
            sage: J = JordanAlgebra(O)
            Traceback (most recent call last):
            ...
            ValueError: 2 must be invertible
        """
        self._O = O
        R = O.base_ring()

        if not R(2).is_unit():
            raise ValueError("2 must be invertible")
        self._half = R(2).inverse_of_unit()

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        Onames = list(O.variable_names())
        Onames.extend(Onames[3] + Onames[i] for i in range(3))
        self._repr_poly_ring = PolynomialRing(R, Onames)

        cat = MagmaticAlgebras(R).Unital().FiniteDimensional().WithBasis()
        Parent.__init__(self, base=R, category=cat)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: JordanAlgebra(O)
            Exceptional Jordan algebra constructed from Octonion algebra
             over Rational Field
        """
        return "Exceptional Jordan algebra constructed from {}".format(self._O)

    def _element_constructor_(self, x):
        r"""
        Construct an element of ``self`` from ``s``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: J = JordanAlgebra(O)
            sage: J(2)
            [2 0 0]
            [0 2 0]
            [0 0 2]
            sage: J([2, 3, 0, -1, O.basis()[3], 2])
            [ 2 -1  k]
            [-1  3  2]
            [-k  2  0]
        """
        R = self.base_ring()
        try:
            x = R(x)
            zero = self._O.zero()
            return self.element_class(self, [x, x, x, zero, zero, zero])
        except (ValueError, TypeError):
            pass
        x = list(x)
        if len(x) != 6:
            raise ValueError("invalid data to construct an element")
        R = self.base_ring()
        for i in range(3):
            x[i] = R(x[i])
            x[3+i] = self._O(x[3+i])
        return self.element_class(self, x)

    def _test_multiplication_self_adjoint(self, **options):
        r"""
        Test that `(XY + YX) / 2` is self-adjoint.

        EXAMPLES::

            sage: O = OctonionAlgebra(GF(7), 1, 3, 4)
            sage: J = JordanAlgebra(O)
            sage: J._test_multiplication_self_adjoint()
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        data_pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        zerO = self._O.zero()
        from sage.misc.misc import some_tuples
        for x, y in some_tuples(S, 2, tester._max_runs):
            SD = x._data
            OD = y._data
            X = [[SD[0], SD[3], SD[4]],
                 [SD[3].conjugate(), SD[1], SD[5]],
                 [SD[4].conjugate(), SD[5].conjugate(), SD[2]]]
            Y = [[OD[0], OD[3], OD[4]],
                  [OD[3].conjugate(), OD[1], OD[5]],
                  [OD[4].conjugate(), OD[5].conjugate(), OD[2]]]
            for r, c in data_pairs:
                if r != c:
                    val = sum(X[r][i] * Y[i][c] + Y[r][i] * X[i][c] for i in range(3)) * self._half
                    val_opp = sum(X[c][i] * Y[i][r] + Y[c][i] * X[i][r] for i in range(3)) * self._half
                    tester.assertEqual(val, val_opp.conjugate())
                else:
                    val = sum(X[r][i] * Y[i][c] + Y[r][i] * X[i][c] for i in range(3)) * self._half
                    tester.assertEqual(val.imag_part(), zerO)

    @cached_method
    def basis(self):
        r"""
        Return a basis of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: J = JordanAlgebra(O)
            sage: B = J.basis()
            sage: B[::6]
            ([1 0 0]
             [0 0 0]
             [0 0 0],
             [ 0  k  0]
             [-k  0  0]
             [ 0  0  0],
             [ 0  0  i]
             [ 0  0  0]
             [-i  0  0],
             [  0   0  lk]
             [  0   0   0]
             [-lk   0   0],
             [  0   0   0]
             [  0   0  li]
             [  0 -li   0])
            sage: len(B)
            27
        """
        import itertools
        R = self.base_ring()
        OB = self._O.basis()
        base = [R.zero()] * 3 + [self._O.zero()] * 3
        ret = []
        for i in range(3):
            temp = list(base)
            temp[i] = R.one()
            ret.append(self.element_class(self, temp))
        for i in range(3):
            for b in OB:
                temp = list(base)
                temp[3+i] = b
                ret.append(self.element_class(self, temp))
        return Family(ret)

    algebra_generators = basis

    def gens(self) -> tuple:
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: J = JordanAlgebra(O)
            sage: G = J.gens()
            sage: G[0]
            [1 0 0]
            [0 0 0]
            [0 0 0]
            sage: G[5]
            [ 0  j  0]
            [-j  0  0]
            [ 0  0  0]
            sage: G[22]
            [ 0  0  0]
            [ 0  0  k]
            [ 0 -k  0]
        """
        return tuple(self.algebra_generators())

    @cached_method
    def zero(self):
        r"""
        Return the additive identity.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: J = JordanAlgebra(O)
            sage: J.zero()
            [0 0 0]
            [0 0 0]
            [0 0 0]
        """
        Rz = self.base_ring().zero()
        Oz = self._O.zero()
        return self.element_class(self, (Rz, Rz, Rz, Oz, Oz, Oz))

    @cached_method
    def one(self):
        r"""
        Return multiplicative identity.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: J = JordanAlgebra(O)
            sage: J.one()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: all(J.one() * b == b for b in J.basis())
            True
        """
        one = self.base_ring().one()
        zero = self._O.zero()
        return self.element_class(self, (one, one, one, zero, zero, zero))

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: J = JordanAlgebra(O)
            sage: J.some_elements()
            [[6/5   0   0]
             [  0 6/5   0]
             [  0   0 6/5],
             [1 0 0]
             [0 1 0]
             [0 0 1],
             [0 0 0]
             [0 0 0]
             [0 0 0],
             [0 0 0]
             [0 1 0]
             [0 0 0],
             [ 0  j  0]
             [-j  0  0]
             [ 0  0  0],
             [  0   0  lj]
             [  0   0   0]
             [-lj   0   0],
             [      0       0       0]
             [      0       1  1/2*lj]
             [      0 -1/2*lj       0],
             [        1         0  j + 2*li]
             [        0         1         0]
             [-j - 2*li         0         1],
             [      1  j + lk       l]
             [-j - lk       0  i + lj]
             [     -l -i - lj       0],
             [     1  3/2*l    2*k]
             [-3/2*l      0  5/2*j]
             [  -2*k -5/2*j      0]]

            sage: O = OctonionAlgebra(GF(3))
            sage: J = JordanAlgebra(O)
            sage: J.some_elements()
            [[-1  0  0]
             [ 0 -1  0]
             [ 0  0 -1],
             [1 0 0]
             [0 1 0]
             [0 0 1],
             [0 0 0]
             [0 0 0]
             [0 0 0],
             [0 0 0]
             [0 1 0]
             [0 0 0],
             [ 0  j  0]
             [-j  0  0]
             [ 0  0  0],
             [  0   0  lj]
             [  0   0   0]
             [-lj   0   0],
             [  0   0   0]
             [  0   1 -lj]
             [  0  lj   0],
             [      1       0  j - li]
             [      0       1       0]
             [-j + li       0       1],
             [      1  j + lk       l]
             [-j - lk       0  i + lj]
             [     -l -i - lj       0],
             [ 1  0 -k]
             [ 0  0  j]
             [ k -j  0]]
        """
        B = self.basis()
        S = [self.an_element(), self.one(), self.zero(),
             B[1], B[5], B[17], B[1] + self._half*B[25],
             self.one() + B[13] + 2*B[16]]
        S.append(sum(B[::5]))
        S.append(sum(self._half * ind * b for ind, b in enumerate(B[::7], start=2)))
        return S

    class Element(AlgebraElement):
        r"""
        An element of an exceptional Jordan algebra.
        """
        def __init__(self, parent, data):
            """
            Initialize ``self``.

            TESTS::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: elt = sum(J.basis())
                sage: TestSuite(elt).run()
            """
            self._data = tuple(data)
            AlgebraElement.__init__(self, parent)

        def _to_print_matrix(self):
            r"""
            Return ``self`` as a matrix for printing.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: elt = J([2, 3, 0, -1 + O.basis()[2], O.basis()[3], -O.basis()[5] + 5*O.basis()[7]])
                sage: elt._to_print_matrix()
                [         2      j - 1          k]
                [    -j - 1          3 -li + 5*lk]
                [        -k  li - 5*lk          0]
            """
            PR = self.parent()._repr_poly_ring
            gens = [PR.one()] + list(PR.gens())
            data = [PR(self._data[i]) for i in range(3)]
            data.extend(PR.sum(c * g for c, g in zip(self._data[3+i].vector(), gens))
                        for i in range(3))
            # add the conjugates
            for i in range(1, 8):
                gens[i] = -gens[i]
            data.extend(PR.sum(c * g for c, g in zip(self._data[3+i].vector(), gens))
                        for i in range(3))
            return matrix(PR, [[data[0], data[3], data[4]], [data[6], data[1], data[5]], [data[7], data[8], data[2]]])

        def _repr_(self):
            r"""
            Return a string representation of ``self``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: J.an_element()
                [6/5   0   0]
                [  0 6/5   0]
                [  0   0 6/5]
            """
            return repr(self._to_print_matrix())

        def _latex_(self):
            r"""
            Return a latex representation of ``self``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: latex(J.an_element())
                \left(\begin{array}{rrr}
                \frac{6}{5} & 0 & 0 \\
                0 & \frac{6}{5} & 0 \\
                0 & 0 & \frac{6}{5}
                \end{array}\right)
            """
            from sage.misc.latex import latex
            return latex(self._to_print_matrix())

        def _ascii_art_(self):
            r"""
            Return an ascii art representation of ``self``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: ascii_art(J.an_element())
                [6/5   0   0]
                [  0 6/5   0]
                [  0   0 6/5]
            """
            from sage.typeset.ascii_art import ascii_art
            return ascii_art(self._to_print_matrix())

        def _unicode_art_(self):
            r"""
            Return a unicode art representation of ``self``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: unicode_art(J.an_element())
                ⎛6/5   0   0⎞
                ⎜  0 6/5   0⎟
                ⎝  0   0 6/5⎠
            """
            from sage.typeset.unicode_art import unicode_art
            return unicode_art(self._to_print_matrix())

        def __bool__(self) -> bool:
            """
            Return if ``self`` is nonzero.

            TESTS::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: all(bool(b) for b in J.basis())
                True
                sage: bool(J.zero())
                False
            """
            return any(d for d in self._data)

        def _richcmp_(self, other, op):
            r"""
            Rich comparison of ``self`` with ``other`` by ``op``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: x = sum(J.basis()[::6])
                sage: y = sum(J.basis()[::5])
                sage: x == x
                True
                sage: x == y
                False
                sage: x < y
                True
                sage: x != J.zero()
                True
            """
            return richcmp(self._data, other._data, op)

        def _add_(self, other):
            """
            Add ``self`` and ``other``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: x = sum(J.basis()[::6])
                sage: y = sum(J.basis()[::5])
                sage: x + x
                [          2         2*k  2*i + 2*lk]
                [       -2*k           0        2*li]
                [-2*i - 2*lk       -2*li           0]
                sage: x + y
                [           2   j + k + lk   i + l + lk]
                [ -j - k - lk            0  i + li + lj]
                [ -i - l - lk -i - li - lj            0]
            """
            return self.__class__(self.parent(), [a + b for a, b in zip(self._data, other._data)])

        def _neg_(self):
            """
            Negate ``self``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: x = sum(J.basis()[::6])
                sage: -x
                [     -1      -k -i - lk]
                [      k       0     -li]
                [ i + lk      li       0]
            """
            return self.__class__(self.parent(), [-c for c in self._data])

        def _sub_(self, other):
            r"""
            Subtract ``other`` from ``self``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: x = sum(J.basis()[::6])
                sage: y = sum(J.basis()[::5])
                sage: x - x
                [0 0 0]
                [0 0 0]
                [0 0 0]
                sage: x - y
                [           0  -j + k - lk   i - l + lk]
                [  j - k + lk            0 -i + li - lj]
                [ -i + l - lk  i - li + lj            0]
            """
            return self.__class__(self.parent(), [a - b for a, b in zip(self._data, other._data)])

        def _mul_(self, other):
            """
            Multiply ``self`` and ``other``.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: x = sum(J.basis()[::7])
                sage: y = sum(J.basis()[::11])
                sage: x * y
                [                    1  -1/2*j + 1/2*l + 1/2  1/2*k + 1/2*lk + 1/2]
                [  1/2*j - 1/2*l + 1/2                     0                -1/2*l]
                [-1/2*k - 1/2*lk + 1/2                 1/2*l                     0]
            """
            P = self.parent()
            SD = self._data
            OD = other._data
            X = [[SD[0], SD[3], SD[4]],
                 [SD[3].conjugate(), SD[1], SD[5]],
                 [SD[4].conjugate(), SD[5].conjugate(), SD[2]]]
            Y = [[OD[0], OD[3], OD[4]],
                  [OD[3].conjugate(), OD[1], OD[5]],
                  [OD[4].conjugate(), OD[5].conjugate(), OD[2]]]
            # we do a simplified multiplication for the diagonal entries since
            # we have, e.g., \alpha * \alpha' + (x (x')^* + x' x^* + y (y')^* + y' y^*) / 2
            ret = [X[0][0] * Y[0][0] + (X[0][1] * Y[1][0]).real_part() + (X[0][2] * Y[2][0]).real_part(),
                   X[1][1] * Y[1][1] + (X[1][0] * Y[0][1]).real_part() + (X[1][2] * Y[2][1]).real_part(),
                   X[2][2] * Y[2][2] + (X[2][0] * Y[0][2]).real_part() + (X[2][1] * Y[1][2]).real_part()]
            ret += [sum(X[r][i] * Y[i][c] + Y[r][i] * X[i][c] for i in range(3)) * P._half
                    for r, c in [(0, 1), (0, 2), (1, 2)]]
            return self.__class__(P, ret)

        def _lmul_(self, other):
            r"""
            Multiply ``self`` by the scalar ``other`` on the left.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: elt = sum(2 * b for b in J.basis()[::6]); elt
                [          2         2*k  2*i + 2*lk]
                [       -2*k           0        2*li]
                [-2*i - 2*lk       -2*li           0]
                sage: elt * 2
                [          4         4*k  4*i + 4*lk]
                [       -4*k           0        4*li]
                [-4*i - 4*lk       -4*li           0]
            """
            return self.__class__(self.parent(), [c * other for c in self._data])

        def _rmul_(self, other):
            r"""
            Multiply ``self`` with the scalar ``other`` by the right
            action.

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: elt = sum(b * 2 for b in J.basis()[::6]); elt
                [          2         2*k  2*i + 2*lk]
                [       -2*k           0        2*li]
                [-2*i - 2*lk       -2*li           0]
                sage: (1/2) * elt
                [      1       k  i + lk]
                [     -k       0      li]
                [-i - lk     -li       0]
            """
            return self.__class__(self.parent(), [other * c for c in self._data])

        def monomial_coefficients(self, copy=True):
            r"""
            Return a dictionary whose keys are indices of basis elements in
            the support of ``self`` and whose values are the corresponding
            coefficients.

            INPUT:

            - ``copy`` -- ignored

            EXAMPLES::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: elt = sum(~QQ(ind) * b for ind, b in enumerate(J.basis()[::6], start=1)); elt
                [              1           1/2*k  1/3*i + 1/4*lk]
                [         -1/2*k               0          1/5*li]
                [-1/3*i - 1/4*lk         -1/5*li               0]
                sage: elt.monomial_coefficients()
                {0: 1, 6: 1/2, 12: 1/3, 18: 1/4, 24: 1/5}

            TESTS::

                sage: O = OctonionAlgebra(QQ)
                sage: J = JordanAlgebra(O)
                sage: all(b.monomial_coefficients() == {i: 1} for i,b in enumerate(J.basis()))
                True
            """
            ret = {}
            for i in range(3):
                if self._data[i]:
                    ret[i] = self._data[i]
                mc = self._data[3+i].monomial_coefficients()
                for k, coeff in mc.items():
                    ret[3+i*8+k] = coeff
            return ret
