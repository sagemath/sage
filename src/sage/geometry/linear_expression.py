"""
Linear Expressions

A linear expression is just a linear polynomial in some (fixed)
variables (allowing a nonzero constant term). This class only implements
linear expressions for others to use.

EXAMPLES::

    sage: from sage.geometry.linear_expression import LinearExpressionModule
    sage: L.<x,y,z> = LinearExpressionModule(QQ);  L
    Module of linear expressions in variables x, y, z over Rational Field
    sage: x + 2*y + 3*z + 4
    x + 2*y + 3*z + 4
    sage: L(4)
    0*x + 0*y + 0*z + 4

You can also pass coefficients and a constant term to construct linear
expressions::

    sage: L([1, 2, 3], 4)
    x + 2*y + 3*z + 4
    sage: L([(1, 2, 3), 4])
    x + 2*y + 3*z + 4
    sage: L([4, 1, 2, 3])   # note: constant is first in single-tuple notation
    x + 2*y + 3*z + 4

The linear expressions are a module over the base ring, so you can
add them and multiply them with scalars::

    sage: m = x + 2*y + 3*z + 4
    sage: 2*m
    2*x + 4*y + 6*z + 8
    sage: m+m
    2*x + 4*y + 6*z + 8
    sage: m-m
    0*x + 0*y + 0*z + 0
"""

from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp
from sage.structure.element import ModuleElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method


class LinearExpression(ModuleElement):
    """
    A linear expression.

    A linear expression is just a linear polynomial in some (fixed)
    variables.

    EXAMPLES::

        sage: from sage.geometry.linear_expression import LinearExpressionModule
        sage: L.<x,y,z> = LinearExpressionModule(QQ)
        sage: m = L([1, 2, 3], 4); m
        x + 2*y + 3*z + 4
        sage: m2 = L([(1, 2, 3), 4]); m2
        x + 2*y + 3*z + 4
        sage: m3 = L([4, 1, 2, 3]); m3   # note: constant is first in single-tuple notation
        x + 2*y + 3*z + 4
        sage: m == m2
        True
        sage: m2 == m3
        True
        sage: L.zero()
        0*x + 0*y + 0*z + 0
        sage: a = L([12, 2/3, -1], -2)
        sage: a - m
        11*x - 4/3*y - 4*z - 6
        sage: LZ.<x,y,z> = LinearExpressionModule(ZZ)
        sage: a - LZ([2, -1, 3], 1)
        10*x + 5/3*y - 4*z - 3
    """
    def __init__(self, parent, coefficients, constant, check=True):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: linear = L([1, 2, 3], 4)   # indirect doctest
            sage: linear.parent() is L
            True

            sage: TestSuite(linear).run()
        """
        super(LinearExpression, self).__init__(parent)
        self._coeffs = coefficients
        self._const = constant
        if check:
            if self._coeffs.parent() is not self.parent().ambient_module():
                raise ValueError("coefficients are not in the ambient module")
            if not self._coeffs.is_immutable():
                raise ValueError("coefficients are not immutable")
            if self._const.parent() is not self.parent().base_ring():
                raise ValueError("the constant is not in the base ring")

    def A(self):
        """
        Return the coefficient vector.

        OUTPUT:

        The coefficient vector of the linear expression.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: linear = L([1, 2, 3], 4);  linear
            x + 2*y + 3*z + 4
            sage: linear.A()
            (1, 2, 3)
            sage: linear.b()
            4
        """
        return self._coeffs

    def b(self):
        """
        Return the constant term.

        OUTPUT:

        The constant term of the linear expression.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: linear = L([1, 2, 3], 4);  linear
            x + 2*y + 3*z + 4
            sage: linear.A()
            (1, 2, 3)
            sage: linear.b()
            4
        """
        return self._const

    constant_term = b

    def coefficients(self):
        """
        Return all coefficients.

        OUTPUT:

        The constant (as first entry) and coefficients of the linear
        terms (as subsequent entries) in a list.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: linear = L([1, 2, 3], 4);  linear
            x + 2*y + 3*z + 4
            sage: linear.coefficients()
            [4, 1, 2, 3]
        """
        return [self._const] + list(self._coeffs)

    dense_coefficient_list = coefficients

    def monomial_coefficients(self, copy=True):
        """
        Return a dictionary whose keys are indices of basis elements in
        the support of ``self`` and whose values are the corresponding
        coefficients.

        INPUT:

        - ``copy`` -- ignored

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: linear = L([1, 2, 3], 4)
            sage: sorted(linear.monomial_coefficients().items(), key=lambda x: str(x[0]))
            [(0, 1), (1, 2), (2, 3), ('b', 4)]
        """
        zero = self.parent().base_ring().zero()
        d = {i: v for i, v in enumerate(self._coeffs) if v != zero}
        if self._const != zero:
            d['b'] = self._const
        return d

    def _repr_vector(self, variable='x'):
        """
        Return a string representation.

        INPUT:

        - ``variable`` -- string; the name of the variable vector

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: L([1, 2, 3], 4)._repr_vector()
            '(1, 2, 3) x + 4 = 0'
            sage: L([-1, -2, -3], -4)._repr_vector('u')
            '(-1, -2, -3) u - 4 = 0'
        """
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        constant = repr(self._const)
        if not atomic_repr:
            constant = '({0})'.format(constant)
        constant = '+ {0}'.format(constant).replace('+ -', '- ')
        return '{0} {1} {2} = 0'.format(repr(self._coeffs), variable, constant)

    def _repr_linear(self, include_zero=True, include_constant=True, multiplication='*'):
        """
        Return a representation as a linear polynomial.

        INPUT:

        - ``include_zero`` -- whether to include terms with zero
          coefficient

        - ``include_constant`` -- whether to include the constant
          term

        - ``multiplication`` -- string (optional, default: ``*``); the
          multiplication symbol to use

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: L([1, 2, 3], 4)._repr_linear()
            'x + 2*y + 3*z + 4'
            sage: L([-1, -2, -3], -4)._repr_linear()
            '-x - 2*y - 3*z - 4'
            sage: L([0, 0, 0], 1)._repr_linear()
            '0*x + 0*y + 0*z + 1'
            sage: L([0, 0, 0], 0)._repr_linear()
            '0*x + 0*y + 0*z + 0'

            sage: R.<u,v> = QQ[]
            sage: L.<x,y,z> = LinearExpressionModule(R)
            sage: L([-u+v+1, -3*u-2, 3], -4*u+v)._repr_linear()
            '(-u + v + 1)*x + (-3*u - 2)*y + 3*z - 4*u + v'

            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: L([1, 0, 3], 0)._repr_linear()
            'x + 0*y + 3*z + 0'
            sage: L([1, 0, 3], 0)._repr_linear(include_zero=False)
            'x + 3*z'
            sage: L([1, 0, 3], 1)._repr_linear(include_constant=False, multiplication='.')
            'x + 0.y + 3.z'
            sage: L([1, 0, 3], 1)._repr_linear(include_zero=False, include_constant=False)
            'x + 3*z'
            sage: L([0, 0, 0], 0)._repr_linear(include_zero=False)
            '0'
        """
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        names = [multiplication + n for n in self.parent()._names]
        terms = list(zip(self._coeffs, names))
        if include_constant:
            terms += [(self._const, '')]
        if not include_zero:
            terms = [t for t in terms if t[0] != 0]
        if len(terms) == 0:
            return '0'
        summands = []
        for coeff, name in terms:
            coeff = str(coeff)
            if not atomic_repr and name != '' and any(c in coeff for c in ['+', '-']):
                coeff = '({0})'.format(coeff)
            summands.append(coeff + name)
        s = ' ' + ' + '.join(summands)
        s = s.replace(' + -', ' - ')
        s = s.replace(' 1' + multiplication, ' ')
        s = s.replace(' -1' + multiplication, ' -')
        return s[1:]

    _repr_ = _repr_linear

    def _add_(self, other):
        """
        Add two linear expressions.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: a = L([1, 2, 3], 4)
            sage: b = L([-1, 3, -3], 0)
            sage: a + b
            0*x + 5*y + 0*z + 4
            sage: a - b
            2*x - y + 6*z + 4
        """
        const = self._const + other._const
        coeffs = self._coeffs + other._coeffs
        coeffs.set_immutable()
        return self.__class__(self.parent(), coeffs, const)

    def _lmul_(self, scalar):
        """
        Multiply a linear expression by a scalar.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: a = L([1, 2, 3], 4);  a
            x + 2*y + 3*z + 4
            sage: 2 * a
            2*x + 4*y + 6*z + 8
            sage: a * 2
            2*x + 4*y + 6*z + 8
            sage: -a
            -x - 2*y - 3*z - 4
            sage: RDF(1) * a
            1.0*x + 2.0*y + 3.0*z + 4.0

        TESTS::

            sage: a._lmul_(2)
            2*x + 4*y + 6*z + 8
        """
        const = scalar * self._const
        coeffs = scalar * self._coeffs
        coeffs.set_immutable()
        return self.__class__(self.parent(), coeffs, const)

    def _acted_upon_(self, scalar, self_on_left):
        """
        Action by scalars that do not live in the base ring.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: a = x + 2*y + 3*z + 4
            sage: a * RDF(3)
            3.0*x + 6.0*y + 9.0*z + 12.0
        """
        base_ring = scalar.base_ring()
        parent = self.parent().change_ring(base_ring)
        changed = parent(self)
        return changed._rmul_(scalar)

    def change_ring(self, base_ring):
        """
        Change the base ring of this linear expression.

        INPUT:

        - ``base_ring`` -- a ring; the new base ring

        OUTPUT:

        A new linear expression over the new base ring.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: a = x + 2*y + 3*z + 4;  a
            x + 2*y + 3*z + 4
            sage: a.change_ring(RDF)
            1.0*x + 2.0*y + 3.0*z + 4.0
        """
        P = self.parent()
        if P.base_ring() is base_ring:
            return self
        return P.change_ring(base_ring)(self)

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x> = LinearExpressionModule(QQ)
            sage: hash(L([0,1])) == hash((1,))
            True
        """
        return hash(self._coeffs) ^ hash(self._const)

    def _richcmp_(self, other, op):
        """
        Compare two linear expressions.

        INPUT:

        - ``other`` -- another linear expression (will be enforced by
          the coercion framework)

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x> = LinearExpressionModule(QQ)
            sage: x == L([0, 1])
            True
            sage: x == x + 1
            False

            sage: M.<x> = LinearExpressionModule(ZZ)
            sage: L.gen(0) == M.gen(0)   # because there is a conversion
            True
            sage: L.gen(0) == L(M.gen(0))   # this is the conversion
            True

            sage: x == 'test'
            False
        """
        return richcmp((self._coeffs, self._const),
                       (other._coeffs, other._const), op)

    def evaluate(self, point):
        """
        Evaluate the linear expression.

        INPUT:

        - ``point`` -- list/tuple/iterable of coordinates; the
          coordinates of a point

        OUTPUT:

        The linear expression `Ax + b` evaluated at the point `x`.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y> = LinearExpressionModule(QQ)
            sage: ex = 2*x + 3* y + 4
            sage: ex.evaluate([1,1])
            9
            sage: ex([1,1])    # syntactic sugar
            9
            sage: ex([pi, e])
            2*pi + 3*e + 4
        """
        try:
            point = self.parent().ambient_module()(point)
        except TypeError:
            from sage.matrix.constructor import vector
            point = vector(point)
        return self._coeffs * point + self._const

    __call__ = evaluate


class LinearExpressionModule(Parent, UniqueRepresentation):
    """
    The module of linear expressions.

    This is the module of linear polynomials which is the parent for
    linear expressions.

    EXAMPLES::

        sage: from sage.geometry.linear_expression import LinearExpressionModule
        sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))
        sage: L
        Module of linear expressions in variables x, y, z over Rational Field
        sage: L.an_element()
        x + 0*y + 0*z + 0
    """
    Element = LinearExpression

    def __init__(self, base_ring, names=tuple()):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))
            sage: type(L)
            <class 'sage.geometry.linear_expression.LinearExpressionModule_with_category'>
            sage: L.base_ring()
            Rational Field

            sage: TestSuite(L).run()

            sage: L = LinearExpressionModule(QQ)
            sage: TestSuite(L).run()
        """
        from sage.categories.modules import Modules
        super(LinearExpressionModule, self).__init__(base_ring, category=Modules(base_ring).WithBasis().FiniteDimensional())
        self._names = names

    @cached_method
    def basis(self):
        """
        Return a basis of ``self``.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))
            sage: list(L.basis())
            [x + 0*y + 0*z + 0,
             0*x + y + 0*z + 0,
             0*x + 0*y + z + 0,
             0*x + 0*y + 0*z + 1]
        """
        from sage.sets.family import Family
        gens = self.gens()
        d = {i: g for i, g in enumerate(gens)}
        d['b'] = self.element_class(self, self.ambient_module().zero(),
                                    self.base_ring().one())
        return Family(list(range(len(gens))) + ['b'], lambda i: d[i])

    @cached_method
    def ngens(self):
        """
        Return the number of linear variables.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))
            sage: L.ngens()
            3
        """
        return len(self._names)

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        OUTPUT:

        A tuple of linear expressions, one for each linear variable.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))
            sage: L.gens()
            (x + 0*y + 0*z + 0, 0*x + y + 0*z + 0, 0*x + 0*y + z + 0)
        """
        from sage.matrix.constructor import identity_matrix
        identity = identity_matrix(self.base_ring(), self.ngens())
        return tuple(self(e, 0) for e in identity.rows())

    def gen(self, i):
        """
        Return the `i`-th generator.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        A linear expression.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))
            sage: L.gen(0)
            x + 0*y + 0*z + 0
        """
        return self.gens()[i]

    def _element_constructor_(self, arg0, arg1=None):
        """
        The element constructor.

        This is part of the Sage parent/element framework.

        TESTS::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))

        Construct from coefficients and constant term::

            sage: L._element_constructor_([1, 2, 3], 4)
            x + 2*y + 3*z + 4
            sage: L._element_constructor_(vector(ZZ, [1, 2, 3]), 4)
            x + 2*y + 3*z + 4

        Construct constant linear expression term::

            sage: L._element_constructor_(4)
            0*x + 0*y + 0*z + 4

        Construct from list/tuple/iterable::

            sage: L._element_constructor_(vector([4, 1, 2, 3]))
            x + 2*y + 3*z + 4

        Construct from a pair ``(coefficients, constant)``::

            sage: L([(1, 2, 3), 4])
            x + 2*y + 3*z + 4

        Construct from linear expression::

            sage: M = LinearExpressionModule(ZZ, ('u', 'v', 'w'))
            sage: m = M([1, 2, 3], 4)
            sage: L._element_constructor_(m)
            x + 2*y + 3*z + 4
        """
        R = self.base_ring()
        if arg1 is None:
            if arg0 in R:
                const = arg0
                coeffs = self.ambient_module().zero()
            elif isinstance(arg0, LinearExpression):
                # Construct from linear expression
                const = arg0.b()
                coeffs = arg0.A()
            elif isinstance(arg0, (list, tuple)) and len(arg0) == 2 and isinstance(arg0[0], (list, tuple)):
                # Construct from pair
                coeffs = arg0[0]
                const = arg0[1]
            else:
                # Construct from list/tuple/iterable::
                try:
                    arg0 = arg0.dense_coefficient_list()
                except AttributeError:
                    arg0 = list(arg0)
                const = arg0[0]
                coeffs = arg0[1:]
        else:
            # arg1 is not None, construct from coefficients and constant term
            coeffs = list(arg0)
            const = arg1
        coeffs = self.ambient_module()(coeffs)
        coeffs.set_immutable()
        const = R(const)
        return self.element_class(self, coeffs, const)

    def random_element(self):
        """
        Return a random element.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x,y,z> = LinearExpressionModule(QQ)
            sage: L.random_element() in L
            True
        """
        A = self.ambient_module().random_element()
        b = self.base_ring().random_element()
        return self(A, b)

    @cached_method
    def ambient_module(self):
        """
        Return the ambient module.

        .. SEEALSO::

            :meth:`ambient_vector_space`

        OUTPUT:

        The domain of the linear expressions as a free module over the
        base ring.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))
            sage: L.ambient_module()
            Vector space of dimension 3 over Rational Field
            sage: M = LinearExpressionModule(ZZ, ('r', 's'))
            sage: M.ambient_module()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: M.ambient_vector_space()
            Vector space of dimension 2 over Rational Field
        """
        from sage.modules.free_module import FreeModule
        return FreeModule(self.base_ring(), self.ngens())

    @cached_method
    def ambient_vector_space(self):
        """
        Return the ambient vector space.

        .. SEEALSO::

            :meth:`ambient_module`

        OUTPUT:

        The vector space (over the fraction field of the base ring)
        where the linear expressions live.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L = LinearExpressionModule(QQ, ('x', 'y', 'z'))
            sage: L.ambient_vector_space()
            Vector space of dimension 3 over Rational Field
            sage: M = LinearExpressionModule(ZZ, ('r', 's'))
            sage: M.ambient_module()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: M.ambient_vector_space()
            Vector space of dimension 2 over Rational Field
        """
        from sage.modules.free_module import VectorSpace
        field = self.base_ring().fraction_field()
        return VectorSpace(field, self.ngens())

    def _coerce_map_from_(self, P):
        """
        Return whether there is a coercion.

        TESTS::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x> = LinearExpressionModule(QQ)
            sage: M.<y> = LinearExpressionModule(ZZ)
            sage: L.coerce_map_from(M)
            Coercion map:
              From: Module of linear expressions in variable y over Integer Ring
              To:   Module of linear expressions in variable x over Rational Field
            sage: M.coerce_map_from(L)

            sage: M.coerce_map_from(ZZ)
            Coercion map:
              From: Integer Ring
              To:   Module of linear expressions in variable y over Integer Ring
            sage: M.coerce_map_from(QQ)
        """
        if self.base().has_coerce_map_from(P):
            return True
        try:
            return self.ngens() == P.ngens() and \
                self.base().has_coerce_map_from(P.base())
        except AttributeError:
            pass
        return super(LinearExpressionModule, self)._coerce_map_from_(P)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: L.<x> = LinearExpressionModule(QQ);  L
            Module of linear expressions in variable x over Rational Field
       """
        return 'Module of linear expressions in variable{2} {0} over {1}'.format(
            ', '.join(self._names), self.base_ring(), 's' if self.ngens() > 1 else '')

    def change_ring(self, base_ring):
        """
        Return a new module with a changed base ring.

        INPUT:

        - ``base_ring`` -- a ring; the new base ring

        OUTPUT:

        A new linear expression over the new base ring.

        EXAMPLES::

            sage: from sage.geometry.linear_expression import LinearExpressionModule
            sage: M.<y> = LinearExpressionModule(ZZ)
            sage: L = M.change_ring(QQ);  L
            Module of linear expressions in variable y over Rational Field

        TESTS::

            sage: L.change_ring(QQ) is L
            True
        """
        return LinearExpressionModule(base_ring, self._names)
