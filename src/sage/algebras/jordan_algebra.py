r"""
Jordan Algebras

AUTHORS:

- Travis Scrimshaw (2014-04-02): initial version
"""

#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import AlgebraElement
from sage.categories.magmatic_algebras import MagmaticAlgebras
from sage.misc.cachefunc import cached_method
from sage.structure.element import is_Matrix
from sage.modules.free_module import FreeModule
from sage.sets.family import Family


class JordanAlgebra(Parent, UniqueRepresentation):
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
        + \underbrace{\beta x + \alpha y}_{\in M}

    where `\alpha, \beta \in R` and `x,y \in M`.

    INPUT:

    Can be either an associative algebra `A` or a symmetric bilinear
    form given as a matrix (possibly followed by, or preceded by, a base
    ring argument)

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
            if not is_Matrix(arg0):
                if arg0.base_ring().characteristic() == 2:
                    raise ValueError("the base ring cannot have characteristic 2")
                return SpecialJordanAlgebra(arg0, names)
            arg0, arg1 = arg0.base_ring(), arg0
        elif is_Matrix(arg0):
            arg0, arg1 = arg1, arg0

        # arg0 is the base ring and arg1 is a matrix
        if not arg1.is_symmetric():
            raise ValueError("the bilinear form is not symmetric")

        arg1 = arg1.change_ring(arg0) # This makes a copy
        arg1.set_immutable()
        return JordanAlgebraSymmetricBilinear(arg0, arg1, names=names)

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
    def gens(self):
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

        def __bool__(self):
            """
            Return if ``self`` is non-zero.

            EXAMPLES::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: J = JordanAlgebra(F)
                sage: a,b,c = map(J, F.gens())
                sage: bool(a + 2*b - c)
                True
            """
            return bool(self._x)

        __nonzero__ = __bool__

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

            - ``copy`` -- (default: ``True``) if ``self`` is internally
              represented by a dictionary ``d``, then make a copy of ``d``;
              if ``False``, then this can cause undesired behavior by
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

    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: m = matrix([[0,1],[1,1]])
            sage: J = JordanAlgebra(m)
            sage: J.basis()
            Family (1 + (0, 0), 0 + (1, 0), 0 + (0, 1))
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

        def __bool__(self):
            """
            Return if ``self`` is non-zero.

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

        __nonzero__ = __bool__

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

