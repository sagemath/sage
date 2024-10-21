# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.combinat
"""
Lie Algebra Elements

AUTHORS:

- Travis Scrimshaw (2013-05-04): Initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport Py_EQ, Py_NE, Py_GT, Py_GE

from sage.misc.repr import repr_lincomb
from sage.structure.element cimport have_same_parent, parent
from sage.structure.coerce cimport coercion_model
from sage.cpython.wrapperdescr cimport wrapperdescr_fastcall
from sage.structure.element_wrapper cimport ElementWrapper
from sage.structure.richcmp cimport richcmp, richcmp_not_equal
from sage.data_structures.blas_dict cimport axpy, add, negate, scal

# TODO: Do we want a dense version?
cdef class LieAlgebraElement(IndexedFreeModuleElement):
    """
    A Lie algebra element.
    """
    # Need to bypass the coercion model
    def __mul__(left, right):
        """
        If we are multiplying two nonzero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: y*x
            x*y - z

        Check that actions work::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: d = L.basis()
            sage: M = L.chargeless_representation(1/2, 3/4)
            sage: d[-5] * M.basis()[10]
            -47/4*v[5]

        TESTS::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: int(3) * x
            3*x
            sage: x * int(3)
            3*x
            sage: y * x.lift()
            x*y - z
            sage: y.lift() * x
            x*y - z
        """
        try:
            # Try the normal coercion first
            return wrapperdescr_fastcall(IndexedFreeModuleElement.__mul__,
                                         left, (right,), <object>NULL)
        except TypeError:
            pass

        # Lift up to the UEA and try multiplication there
        # We will eventually want to lift stuff up anyways,
        #   so just do it here.
        if isinstance(left, LieAlgebraElement):
            left = (<LieAlgebraElement> left).lift()
        if isinstance(right, LieAlgebraElement):
            right = (<LieAlgebraElement> right).lift()
        return left * right

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` in ``codomain`` under the
        map that sends the generators of the parent of ``self``
        to the elements of the tuple ``im_gens``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            doctest:warning...:
            FutureWarning: The Hall basis has not been fully proven correct, but currently no bugs are known
            See https://github.com/sagemath/sage/issues/16823 for details.
            sage: elt = Lyn.an_element()
            sage: elt._im_gens_(H, H.gens())
            x + y + z
            sage: x, y, z = Lyn.gens()
            sage: elt = x.bracket(y) + y.bracket(z) + z.bracket(x); elt
            [x, y] - [x, z] + [y, z]
            sage: elt._im_gens_(H, H.gens())
            [x, y] - [x, z] + [y, z]
            sage: xx, yy, zz = H.gens()
            sage: elt._im_gens_(H, [xx, yy, xx])
            0
            sage: elt._im_gens_(H, [xx, yy, yy+zz])
            -[x, z] + [y, z]
            sage: L2 = LieAlgebra(QQ, 'a,b')
            sage: Lyn2 = L2.Lyndon()
            sage: a, b = Lyn2.gens()
            sage: elt._im_gens_(Lyn2, [a, b, a.bracket(b)])
            -[a, [a, b]] + [a, b] - [[a, b], b]
        """
        s = codomain.zero()
        if not self:  # If we are 0
            return s
        names = self.parent().variable_names()
        if base_map is None:
            def base_map(x):
                return x

        return codomain.sum(base_map(c) * t._im_gens_(codomain, im_gens, names)
                            for t, c in self._monomial_coefficients.items())

    cpdef lift(self):
        """
        Lift ``self`` to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: x.lift().parent() == L.universal_enveloping_algebra()
            True

        TESTS::

            sage: L = lie_algebras.pwitt(GF(5), 5); L
            The 5-Witt Lie algebra over Finite Field of size 5
            sage: x = L.basis()[2]
            sage: y = L.basis()[3]
            sage: x.lift()
            b2
            sage: y.lift()
            b3
            sage: x * y
            b2*b3
            sage: y * x
            b2*b3 + b0

            sage: L = lie_algebras.regular_vector_fields(QQ)
            sage: L.an_element()
            d[-1] + d[0] - 3*d[1]
            sage: L.an_element().lift()
            PBW[-1] + PBW[0] - 3*PBW[1]
        """
        UEA = self._parent.universal_enveloping_algebra()
        s = UEA.zero()
        if not self:
            return s
        # Special hook for when the index set of the parent of ``self``
        #   does not match the generators index set of the UEA.
        if hasattr(self._parent, '_UEA_names_map'):
            names_map = self._parent._UEA_names_map
            gen_dict = UEA.gens_dict()
            for t, c in self._monomial_coefficients.items():
                s += c * gen_dict[names_map[t]]
        else:
            try:
                gen_dict = UEA.algebra_generators()
            except (TypeError, AttributeError):
                gen_dict = UEA.gens_dict()
            for t, c in self._monomial_coefficients.items():
                s += c * gen_dict[t]
        return s

cdef class LieAlgebraElementWrapper(ElementWrapper):
    """
    Wrap an element as a Lie algebra element.

    TESTS:

    We check comparisons::

        sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
        sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
        True

    The next doctests show similar behavior, although on elements of
    other classes::

        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 3)
        sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
        True

        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
        sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
        True

    Check inequality::

        sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
        sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
        False
        sage: L.zero() == 0
        True
        sage: L.zero() != 0
        False

    The next doctests show similar behavior, although on elements of
    other classes::

        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 3)
        sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
        False
        sage: L.an_element()
        X + Y + Z
        sage: L.an_element() == 0
        False
        sage: L.an_element() != 0
        True

        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
        sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
        False
        sage: L.zero() == 0
        True
        sage: L.zero() != 0
        False
        sage: L.zero() >= 0
        True
        sage: L.zero() < 0
        False

    We check the display of elements::

        sage: R = FreeAlgebra(QQ, 3, 'x')
        sage: L.<l0,l1,l2> = LieAlgebra(associative=R.gens())
        sage: elt = l0 + l1
        sage: elt
        x0 + x1
        sage: latex(elt)
        x_{0} + x_{1}

        sage: s = SymmetricFunctions(QQ).s()
        sage: L = LieAlgebra(associative=s)
        sage: P = Partition([4,2,2,1])
        sage: x = L.basis()[P]
        sage: ascii_art(x)
        s
         ****
         **
         **
         *
        sage: unicode_art(x)
        s
         ┌┬┬┬┐
         ├┼┼┴┘
         ├┼┤
         ├┼┘
         └┘
    """
    def __bool__(self):
        """
        Return if ``self`` is nonzero.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: bool(L.zero())
            False
            sage: bool(x + y)
            True
        """
        return bool(self.value)

    cpdef _add_(self, right):
        """
        Add ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x + y
            x + y
        """
        return type(self)(self._parent, self.value + right.value)

    cpdef _sub_(self, right):
        """
        Subtract ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x - y
            x - y
        """
        return type(self)(self._parent, self.value - right.value)

    # Need to bypass the coercion model
    def __mul__(left, right):
        """
        If we are multiplying two nonzero elements, automatically
        lift up to the universal enveloping algebra.

        .. TODO::

            Write more tests for this method.

        EXAMPLES::

            sage: S = SymmetricGroup(3).algebra(QQ)
            sage: L = LieAlgebra(associative=S)
            sage: x = L.gen(2); x
            (1,2,3)
            sage: y = L.gen(3); y
            (2,3)
            sage: u = x*3; u
            3*(1,2,3)
            sage: parent(u) == L
            True
            sage: u = x*(3/2); u
            3/2*(1,2,3)
            sage: parent(u) == L
            True
            sage: elt = x*y - y*x; elt
            b4 - b5
            sage: xp, yp = x.lift_associative(), y.lift_associative()
            sage: eltp = xp*yp - yp*xp; eltp
            -(1,2) + (1,3)
            sage: G = list(S.basis())
            sage: G[4] - G[5]
            -(1,2) + (1,3)

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L.<x,y> = LieAlgebra(associative=S.gens())
            sage: int(3) * x
            3*(1,2,3)
            sage: y * int(3)
            3*(1,2)
        """
        try:
            # Try the normal coercion first
            return wrapperdescr_fastcall(ElementWrapper.__mul__,
                                         left, (right,), <object>NULL)
        except TypeError:
            pass

        # Lift up to the UEA and try multiplication there
        # We will eventually want to lift stuff up anyways,
        #   so just do it here.
        if isinstance(left, LieAlgebraElementWrapper):
            left = (<LieAlgebraElementWrapper> left).lift()
        if isinstance(right, LieAlgebraElementWrapper):
            right = (<LieAlgebraElementWrapper> right).lift()
        return left * right

    def __truediv__(self, x):
        """
        Division by coefficients.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 3)
            sage: x = L.an_element(); x
            p1
            sage: x / 2
            1/2*p1
        """
        return self * (~x)

    cpdef _acted_upon_(self, scalar, bint self_on_left):
        """
        Return the action of a scalar on ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: 3*x
            3*x
            sage: parent(3*x) == parent(x)
            True
            sage: x / 2
            1/2*x
            sage: y * (1/2)
            1/2*y
            sage: y * 1/2
            1/2*y
            sage: 1/2 * y
            1/2*y
            sage: QQ(1/2) * y
            1/2*y
        """
        # This was copied and IDK if it still applies (TCS):
        # With the current design, the coercion model does not have
        # enough information to detect apriori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        scalar_parent = parent(scalar)
        if scalar_parent != self._parent.base_ring():
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if self._parent.base_ring().has_coerce_map_from(scalar_parent):
                scalar = self._parent.base_ring()(scalar)
            else:
                return None
        if self_on_left:
            return type(self)(self._parent, self.value * scalar)
        return type(self)(self._parent, scalar * self.value)

    def __neg__(self):
        """
        Return the negation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: -x
            -x
        """
        return type(self)(self._parent, -self.value)

    def __getitem__(self, i):
        """
        Redirect the ``__getitem__()`` to the wrapped element.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
            sage: m = L.gen(0)
            sage: m[0,0]
            0
            sage: m[0][1]
            1
        """
        return self.value.__getitem__(i)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = G.algebra(QQ)
            sage: L = LieAlgebra(associative=S)
            sage: x = L.an_element() + L.basis()[G.one()]
            sage: x
            2*() + (2,3) + (1,2) + (1,2,3) + (1,3,2) + (1,3)
            sage: sorted(x)
            [((), 2), ((2,3), 1), ((1,2), 1), ((1,2,3), 1),
             ((1,3,2), 1), ((1,3), 1)]
        """
        cdef dict d = self.value.monomial_coefficients(copy=False)
        yield from d.items()


# TODO: Also used for vectors, find a better name
cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    """
    Lie algebra element wrapper around a matrix.
    """
    def __init__(self, parent, value):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation='matrix')
            sage: z = L.z()
            sage: z.value.is_immutable()
            True
        """
        value.set_immutable()  # Make the matrix immutable for hashing
        LieAlgebraElementWrapper.__init__(self, parent, value)


cdef class LieSubalgebraElementWrapper(LieAlgebraElementWrapper):
    r"""
    Wrap an element of the ambient Lie algebra as an element.
    """
    def __init__(self, parent, value):
        """
        Initialize ``self``.

        TESTS::

            sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}})
            sage: S = L.subalgebra([X, Y])
            sage: TestSuite(S(X)).run()
        """
        LieAlgebraElementWrapper.__init__(self, parent, value)
        self._monomial_coefficients = None

    def __getitem__(self, i):
        r"""
        Return the coefficient of ``self`` indexed by ``i``.

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}})
            sage: S = L.subalgebra([X, Y])
            sage: el = S(2*Y + 9*Z)
            sage: el[1]
            2
            sage: el[2]
            9
        """
        if self._monomial_coefficients is None:
            # This sets _monomial_coefficients
            self.monomial_coefficients(copy=False)
        try:
            return self._monomial_coefficients[i]
        except KeyError:
            return self.parent().base_ring().zero()

    def _bracket_(self, x):
        """
        Return the Lie bracket ``[self, x]``.

        Assumes ``x`` and ``self`` have the same parent.

        INPUT:

        - ``x`` -- an element of the same Lie subalgebra as ``self``

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}})
            sage: S = L.subalgebra([X, Y])
            sage: S(X)._bracket_(S(Y))
            Z
        """
        x_lift = (<LieSubalgebraElementWrapper> x).value
        return type(self)(self._parent, self.value._bracket_(x_lift))

    def _vector_(self, sparse=False, order=None):
        r"""
        Return the vector in ``g.module()`` corresponding to the
        element ``self`` of ``g`` (where ``g`` is the parent of ``self``).

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: S = L.subalgebra([X, Y])
            sage: S.basis()
            Family (X, Y, 3*Z)
            sage: S(2*Y + 9*Z).to_vector()
            (0, 2, 9)
            sage: S2 = L.subalgebra([Y, Z])
            sage: S2.basis()
            Family (Y, Z)
            sage: S2(2*Y + 9*Z).to_vector()
            (0, 2, 9)

        TESTS::

            sage: L.<X,Y> = LieAlgebra(ZZ, abelian=True)
            sage: S = L.subalgebra(X)
            sage: S(X).to_vector() in S.module()
            True
            sage: S(X).to_vector().parent() is S.module()
            True
        """
        return self._parent.module()(self.value.to_vector(sparse=sparse))

    to_vector = _vector_

    cpdef dict monomial_coefficients(self, bint copy=True):
        r"""
        Return a dictionary whose keys are indices of basis elements
        in the support of ``self`` and whose values are the
        corresponding coefficients.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``self`` is internally
          represented by a dictionary ``d``, then make a copy of ``d``;
          if ``False``, then this can cause undesired behavior by
          mutating ``d``

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: S = L.subalgebra([X, Y])
            sage: S(2*Y + 9*Z).monomial_coefficients()
            {1: 2, 2: 3}
            sage: S2 = L.subalgebra([Y, Z])
            sage: S2(2*Y + 9*Z).monomial_coefficients()
            {0: 2, 1: 9}
        """
        cdef Py_ssize_t k
        if self._monomial_coefficients is None:
            sm = self.parent().module()
            v = sm.coordinate_vector(self.to_vector())
            self._monomial_coefficients = {k: v[k] for k in range(len(v))
                                           if v[k]}
        if copy:
            return dict(self._monomial_coefficients)
        return self._monomial_coefficients

    cpdef _add_(self, right):
        """
        Add ``self`` and ``rhs``.

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: S = L.subalgebra([X, Y])
            sage: a = S(2*Y + 12*Z)
            sage: b = S(X + 2*Y)
            sage: (a + b).monomial_coefficients()
            {0: 1, 1: 4, 2: 4}
            sage: a.monomial_coefficients()        # We set a._monomial_coefficients
            {1: 2, 2: 4}
            sage: b.monomial_coefficients()        # We set b._monomial_coefficients
            {0: 1, 1: 2}
            sage: (a + b).monomial_coefficients()  # This is now computed from a and b
            {0: 1, 1: 4, 2: 4}
        """
        cdef LieSubalgebraElementWrapper ret, other = <LieSubalgebraElementWrapper> right
        ret = type(self)(self._parent, self.value + other.value)
        if self._monomial_coefficients is not None and other._monomial_coefficients is not None:
            mc = add(self._monomial_coefficients, other._monomial_coefficients)
            ret._monomial_coefficients = mc
        return ret

    cpdef _sub_(self, right):
        """
        Subtract ``self`` and ``rhs``.

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: S = L.subalgebra([X, Y])
            sage: a = S(2*Y + 12*Z)
            sage: b = S(X + 2*Y)
            sage: (a - b).monomial_coefficients()
            {0: -1, 2: 4}
            sage: a.monomial_coefficients()        # We set a._monomial_coefficients
            {1: 2, 2: 4}
            sage: b.monomial_coefficients()        # We set b._monomial_coefficients
            {0: 1, 1: 2}
            sage: (a - b).monomial_coefficients()  # This is now computed from a and b
            {0: -1, 2: 4}
        """
        cdef LieSubalgebraElementWrapper ret, other = <LieSubalgebraElementWrapper> right
        ret = type(self)(self._parent, self.value - other.value)
        if self._monomial_coefficients is not None and other._monomial_coefficients is not None:
            mc = axpy(-1, other._monomial_coefficients, self._monomial_coefficients)
            ret._monomial_coefficients = mc
        return ret

    cpdef _acted_upon_(self, scalar, bint self_on_left):
        """
        Return the action of a scalar on ``self``.

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: S = L.subalgebra([X, Y])
            sage: a = S(2*Y + 12*Z)
            sage: (2*a).monomial_coefficients()
            {1: 4, 2: 8}
            sage: a.monomial_coefficients()      # We set a._monomial_coefficients
            {1: 2, 2: 4}
            sage: (2*a).monomial_coefficients()  # This is now computed from a
            {1: 4, 2: 8}
        """
        # This was copied and IDK if it still applies (TCS):
        # With the current design, the coercion model does not have
        # enough information to detect apriori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        scalar_parent = parent(scalar)
        if scalar_parent != self._parent.base_ring():
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if self._parent.base_ring().has_coerce_map_from(scalar_parent):
                scalar = self._parent.base_ring()(scalar)
            else:
                return None
        cdef LieSubalgebraElementWrapper ret
        if self_on_left:
            ret = type(self)(self._parent, self.value * scalar)
        else:
            ret = type(self)(self._parent, scalar * self.value)
        if self._monomial_coefficients is not None:
            ret._monomial_coefficients = scal(scalar, self._monomial_coefficients, self_on_left)
        return ret

cdef class StructureCoefficientsElement(LieAlgebraMatrixWrapper):
    """
    An element of a Lie algebra given by structure coefficients.
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: x - 3/2 * y
            x - 3/2*y
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self._parent._print_options['scalar_mult'],
                            repr_monomial=self._parent._repr_generator,
                            strip_one=True)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: elt = x - 3/2 * y
            sage: latex(elt)
            x - \frac{3}{2} y
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self._parent._print_options['scalar_mult'],
                            latex_scalar_mult=self._parent._print_options['latex_scalar_mult'],
                            repr_monomial=self._parent._latex_term,
                            is_latex=True, strip_one=True)

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: ascii_art(x - 3/2 * y)
            x - 3/2*y
        """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art(repr_lincomb(self._sorted_items_for_printing(),
                                      scalar_mult=ascii_art(self._parent._print_options['scalar_mult']),
                                      repr_monomial=ascii_art,
                                      strip_one=True))

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: unicode_art(x - 3/2 * y)
            x - 3/2·y
        """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art(repr_lincomb(self._sorted_items_for_printing(),
                                        scalar_mult='·',
                                        strip_one=True))

    cpdef bracket(self, right):
        """
        Return the Lie bracket ``[self, right]``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}})
            sage: x.bracket(y)
            z
            sage: y.bracket(x)
            -z
            sage: (x + y - z).bracket(x - y + z)
            -2*y - 2*z
        """
        if not have_same_parent(self, right):
            self, right = coercion_model.canonical_coercion(self, right)
        return self._bracket_(right)

    # We need this method because the LieAlgebra.bracket method (from the
    #   category) calls this, where we are guaranteed to have the same parent.
    cpdef _bracket_(self, right):
        """
        Return the Lie bracket ``[self, right]``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}})
            sage: x._bracket_(y)
            z
            sage: y._bracket_(x)
            -z
        """
        P = self._parent
        cdef dict s_coeff = P._s_coeff
        d = P.dimension()
        cdef list ret = [P.base_ring().zero()]*d
        cdef int i1, i2, i3
        cdef StructureCoefficientsElement rt = <StructureCoefficientsElement> right
        for i1 in range(d):
            c1 = self.value[i1]
            if not c1:
                continue
            for i2 in range(d):
                c2 = rt.value[i2]
                if not c2:
                    continue
                prod_c1_c2 = c1 * c2
                if (i1, i2) in s_coeff:
                    v = s_coeff[i1, i2]
                    for i3 in range(d):
                        ret[i3] += prod_c1_c2 * v[i3]
                elif (i2, i1) in s_coeff:
                    v = s_coeff[i2, i1]
                    for i3 in range(d):
                        ret[i3] -= prod_c1_c2 * v[i3]
        return type(self)(P, P._M(ret))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: elt = x - 3/2 * y
            sage: list(elt)
            [('x', 1), ('y', -3/2)]
        """
        zero = self._parent.base_ring().zero()
        I = self._parent._indices
        cdef int i
        for i, v in enumerate(self.value):
            if v != zero:
                yield (I[i], v)

    cpdef _vector_(self, bint sparse=False, order=None):
        """
        Return ``self`` as a vector.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: a = x + 3*y - z/5
            sage: vector(a)
            (1, 3, -1/5)
        """
        if sparse:
            return self.value.sparse_vector()
        return self.value

    cpdef to_vector(self, bint sparse=False, order=None):
        """
        Return ``self`` as a vector.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: a = x + 3*y - z/2
            sage: a.to_vector()
            (1, 3, -1/2)
        """
        return self._vector_(sparse=sparse)

    def lift(self):
        """
        Return the lift of ``self`` to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: elt = x - 3/2 * y
            sage: l = elt.lift(); l
            x - 3/2*y
            sage: l.parent()
            Noncommutative Multivariate Polynomial Ring in x, y
             over Rational Field, nc-relations: {y*x: x*y - x}
        """
        UEA = self._parent.universal_enveloping_algebra()
        gens = UEA.gens()
        return UEA.sum(c * gens[i] for i, c in self.value.items())

    cpdef dict monomial_coefficients(self, bint copy=True):
        """
        Return the monomial coefficients of ``self`` as a dictionary.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: a = 2*x - 3/2*y + z
            sage: a.monomial_coefficients()
            {'x': 2, 'y': -3/2, 'z': 1}
            sage: a = 2*x - 3/2*z
            sage: a.monomial_coefficients()
            {'x': 2, 'z': -3/2}
        """
        I = self._parent._indices
        return {I[i]: v for i, v in self.value.items()}

    def __getitem__(self, i):
        """
        Return the coefficient of the basis element indexed by ``i``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: elt = x - 3/2 * y
            sage: elt['y']
            -3/2
        """
        return self.value[self._parent._indices.index(i)]


cdef class UntwistedAffineLieAlgebraElement(Element):
    """
    An element of an untwisted affine Lie algebra.
    """
    def __init__(self, parent, dict t_dict, c_coeff, d_coeff):
        """
        Initialize ``self``.

        TESTS::

            sage: L = lie_algebras.Affine(QQ, ['A',2,1])
            sage: x = L.an_element()
            sage: TestSuite(x).run()
        """
        Element.__init__(self, parent)
        self._t_dict = t_dict
        self._c_coeff = c_coeff
        self._d_coeff = d_coeff
        self._hash = -1

    def __reduce__(self):
        """
        Used in pickling.

        TESTS::

            sage: L = lie_algebras.Affine(QQ, ['B',3,1])
            sage: x = L.an_element()
            sage: loads(dumps(x)) == x
            True
        """
        return (_build_untwisted_affine_element,
                (self._parent, self._t_dict, self._c_coeff, self._d_coeff))

    def _repr_generic(self, style, coeff, t_disp, mult, tensor_symb):
        """
        Return a representation of ``self`` based on ``style``.

        INPUT:

        - ``style`` -- a function for how to convert the objects
        - ``coeff`` -- a function for how to display the coefficients
        - ``t_disp`` -- a function for how to display the powers of `t`
        - ``mult`` -- the multiplication symbol; must be compatible
          with ``style``
        - ``tensor_symb`` -- the tensor symbol; must be compatible
          with ``style``

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['B', 3, 1])
            sage: elt = L.an_element()
            sage: elt._repr_generic(str, str, lambda t: "T^{}".format(t), '.', '(x)')
            '(E[alpha[3]] + E[alpha[2]] + E[alpha[1]] + h1 + h2 + h3 + E[-alpha[3]]
             + E[-alpha[2]] + E[-alpha[1]])(x)T^0 + (E[-alpha[1] - 2*alpha[2]
             - 2*alpha[3]])(x)T^1 + (E[alpha[1] + 2*alpha[2] + 2*alpha[3]])(x)T^-1 + c + d'
        """
        ret = style('')
        mult = style(mult)
        tensor_symb = style(tensor_symb)
        for t, g in self._t_dict.items():
            if ret:
                ret += style(' + ')
            if coeff == str:
                # We need to special case this because of the necessary added
                #   comma by Python
                ret += "({})".format(g) + tensor_symb + style(t_disp(t))
            else:
                ret += coeff((g,)) + tensor_symb + style(t_disp(t))
        if self._c_coeff != 0:
            if ret:
                ret += style(' + ')
            if self._c_coeff != 1:
                ret += coeff(self._c_coeff) + mult + style('c')
            else:
                ret += style('c')

        if self._d_coeff != 0:
            if ret:
                ret += style(' + ')
            if self._d_coeff != 1:
                ret += coeff(self._d_coeff) + mult + style('d')
            else:
                ret += style('d')

        if not ret:
            return style('0')
        return ret

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: list(L.lie_algebra_generators())
            [(E[alpha[1]])#t^0,
             (E[-alpha[1]])#t^0,
             (h1)#t^0,
             (E[-alpha[1]])#t^1,
             (E[alpha[1]])#t^-1,
             c,
             d]
            sage: L.an_element()
            (E[alpha[1]] + h1 + E[-alpha[1]])#t^0
             + (E[-alpha[1]])#t^1 + (E[alpha[1]])#t^-1
             + c + d
            sage: L.zero()
            0

            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: e1 + 2*f1 - h1 + e0 + 3*c - 2*d
            (E[alpha[1]] - h1 + 2*E[-alpha[1]])#t^0 + (E[-alpha[1]])#t^1
             + 3*c + -2*d
        """
        return self._repr_generic(str, str, lambda t: "t^{}".format(t), '*', '#')

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: [latex(g) for g in L.lie_algebra_generators()]
            [\left(E_{\alpha_{1}}\right) \otimes t^{0},
             \left(E_{-\alpha_{1}}\right) \otimes t^{0},
             \left(E_{\alpha^\vee_{1}}\right) \otimes t^{0},
             \left(E_{-\alpha_{1}}\right) \otimes t^{1},
             \left(E_{\alpha_{1}}\right) \otimes t^{-1},
             c,
             d]
            sage: latex(L.an_element())
            \left(E_{\alpha_{1}} + E_{\alpha^\vee_{1}} + E_{-\alpha_{1}}\right) \otimes t^{0}
             + \left(E_{-\alpha_{1}}\right) \otimes t^{1}
             + \left(E_{\alpha_{1}}\right) \otimes t^{-1}
             + c + d
            sage: latex(L.zero())
            0

            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: latex(e1 + 2*f1 - h1 + e0 + 3*c - 2*d)
            \left(E_{\alpha_{1}} - E_{\alpha^\vee_{1}} + 2 E_{-\alpha_{1}}\right) \otimes t^{0}
             + \left(E_{-\alpha_{1}}\right) \otimes t^{1} + 3 c + -2 d
        """
        from sage.misc.latex import latex
        return self._repr_generic(str, latex, lambda t: "t^{{{}}}".format(t), ' ', ' \\otimes ')

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: unicode_art([g for g in L.lie_algebra_generators()])
            [ ( alpha[1] )⊗t⁰, ( -alpha[1] )⊗t⁰, ( alphacheck[1] )⊗t⁰, ( -alpha[1] )⊗t¹,
            <BLANKLINE>
             ( alpha[1] )⊗t⁻¹, c, d ]
            sage: unicode_art(L.an_element())
            ( alpha[1] + alphacheck[1] + -alpha[1] )⊗t⁰ + ( -alpha[1] )⊗t¹ + ( alpha[1] )⊗
            <BLANKLINE>
            t⁻¹ + c + d
            sage: unicode_art(L.zero())
            0

            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: unicode_art(e1 + 2*f1 - h1 + e0 + 3*c - 2*d)
            ( alpha[1] - alphacheck[1] + 2·-alpha[1] )⊗t⁰ + ( -alpha[1] )⊗t¹ + 3⋅c + -2⋅d
        """
        from sage.typeset.unicode_art import unicode_art, unicode_superscript
        return self._repr_generic(unicode_art, unicode_art,
                                  lambda t: "t" + unicode_superscript(t),
                                  unicode_art('⋅'), unicode_art('⊗'))

    cpdef dict t_dict(self):
        r"""
        Return the ``dict``, whose keys are powers of `t` and values are
        elements of the classical Lie algebra, of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: x = L.an_element()
            sage: x.t_dict()
            {-1: E[alpha[1]],
             0: E[alpha[1]] + h1 + E[-alpha[1]],
             1: E[-alpha[1]]}
        """
        return self._t_dict.copy()

    cpdef c_coefficient(self):
        r"""
        Return the coefficient of `c` of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: x = L.an_element() - 3 * L.c()
            sage: x.c_coefficient()
            -2
        """
        return self._c_coeff

    cpdef d_coefficient(self):
        r"""
        Return the coefficient of `d` of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: x = L.an_element() + L.d()
            sage: x.d_coefficient()
            2
        """
        return self._d_coeff

    cpdef _richcmp_(self, other, int op):
        """
        Return the rich comparison of ``self`` with ``other``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['C',2,1])
            sage: x = L.an_element()
            sage: c = L.basis()['c']
            sage: d = L.basis()['d']
            sage: c == d
            False
            sage: x != c
            True
            sage: 2*c - d == c + c - d
            True
            sage: x - c != x - c
            False
            sage: x - c != x - d
            True
        """
        if op != Py_EQ and op != Py_NE:
            return NotImplemented
        cdef UntwistedAffineLieAlgebraElement rt = <UntwistedAffineLieAlgebraElement> other
        return richcmp((self._t_dict, self._c_coeff, self._d_coeff),
                       (rt._t_dict, rt._c_coeff, rt._d_coeff),
                       op)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: asl = lie_algebras.Affine(QQ, ['A',4,1])
            sage: x = asl.an_element()
            sage: hash(x) == hash(x)
            True
            sage: hash(asl.zero())
            0
        """
        if not self:
            self._hash = 0
        if self._hash == -1:
            self._hash = hash((tuple([self._t_dict[i] for i in sorted(self._t_dict)]),
                               self._c_coeff, self._d_coeff))
        return self._hash

    def __bool__(self):
        """
        Return ``self`` as a boolean.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['C',2,1])
            sage: x = L.an_element()
            sage: bool(x)
            True
            sage: bool(L.zero())
            False
        """
        return bool(self._t_dict) or bool(self._c_coeff) or bool(self._d_coeff)

    cpdef _add_(self, other):
        """
        Add ``self`` and ``other``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: e0.bracket(e1) + d + c + 3*d
            (-h1)#t^1 + c + 4*d
        """
        cdef UntwistedAffineLieAlgebraElement rt = <UntwistedAffineLieAlgebraElement> other
        return type(self)(self._parent, add(self._t_dict, rt._t_dict),
                          self._c_coeff + rt._c_coeff,
                          self._d_coeff + rt._d_coeff)

    cpdef _sub_(self, other):
        """
        Subtract ``self`` and ``other``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: d - e1 + c - 3*d
            (-E[alpha[1]])#t^0 + c + -2*d
            sage: 4*c - e0.bracket(f0)
            (h1)#t^0
            sage: 4*c - e0.bracket(f0) - h1
            0
            sage: 4*c - e0.bracket(f0) - h1 == L.zero()
            True
            sage: e1 - f1
            (E[alpha[1]] - E[-alpha[1]])#t^0
        """
        cdef UntwistedAffineLieAlgebraElement rt = <UntwistedAffineLieAlgebraElement> other
        return type(self)(self._parent, axpy(-1, rt._t_dict, self._t_dict),
                          self._c_coeff - rt._c_coeff,
                          self._d_coeff - rt._d_coeff)

    cpdef _neg_(self):
        """
        Negate ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: x = e0.bracket(e1) + d + e1 + c + 3*d
            sage: -x + e1
            (h1)#t^1 + -1*c + -4*d
        """
        return type(self)(self._parent, negate(self._t_dict),
                          -self._c_coeff, -self._d_coeff)

    cpdef _acted_upon_(self, scalar, bint self_on_left):
        """
        Return ``self`` acted upon by ``scalar``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: x = e1 + f0.bracket(f1) + 3*c - 2/5 * d
            sage: x
            (E[alpha[1]])#t^0 + (h1)#t^-1 + 3*c + -2/5*d
            sage: -2 * x
            (-2*E[alpha[1]])#t^0 + (-2*h1)#t^-1 + -6*c + 4/5*d
        """
        # This was copied and IDK if it still applies (TCS):
        # With the current design, the coercion model does not have
        # enough information to detect apriori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        scalar_parent = parent(scalar)
        if scalar_parent != self._parent.base_ring():
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if self._parent.base_ring().has_coerce_map_from(scalar_parent):
                scalar = self._parent.base_ring()(scalar)
            else:
                return None
        return type(self)(self._parent, scal(scalar, self._t_dict, self_on_left),
                          scalar * self._c_coeff,
                          scalar * self._d_coeff)

    cpdef monomial_coefficients(self, bint copy=True):
        """
        Return the monomial coefficients of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['C',2,1])
            sage: x = L.an_element()
            sage: sorted(x.monomial_coefficients(), key=str)
            [(-2*alpha[1] - alpha[2], 1),
             (-alpha[1], 0),
             (-alpha[2], 0),
             (2*alpha[1] + alpha[2], -1),
             (alpha[1], 0),
             (alpha[2], 0),
             (alphacheck[1], 0),
             (alphacheck[2], 0),
             'c',
             'd']
        """
        cdef dict d = {}
        for t, g in self._t_dict.items():
            for k, c in g.monomial_coefficients(copy=False).items():
                d[k, t] = c
        if self._c_coeff:
            d['c'] = self._c_coeff
        if self._d_coeff:
            d['d'] = self._d_coeff
        return d

    cpdef bracket(self, right):
        """
        Return the Lie bracket ``[self, right]``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',1,1])
            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: e0.bracket(f0)
            (-h1)#t^0 + 4*c
            sage: e1.bracket(0)
            0
            sage: e1.bracket(1)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
             'Affine Kac-Moody algebra of ['A', 1] in the Chevalley basis'
             and 'Integer Ring'
        """
        if not have_same_parent(self, right):
            self, right = coercion_model.canonical_coercion(self, right)
        return self._bracket_(right)

    cpdef _bracket_(self, y):
        """
        Return the Lie bracket ``[self, y]``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',1,1])
            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: al = RootSystem(['A',1]).root_lattice().simple_roots()
            sage: x = L.basis()[al[1], 5]
            sage: y = L.basis()[-al[1], -3]
            sage: z = L.basis()[-al[1], -5]
            sage: x._bracket_(y)
            (h1)#t^2
            sage: x._bracket_(z)
            (h1)#t^0 + 20*c
            sage: x._bracket_(e1)
            0
            sage: x._bracket_(f1)
            (h1)#t^5
            sage: x._bracket_(h1)
            (-2*E[alpha[1]])#t^5
            sage: x._bracket_(d)
            (-5*E[alpha[1]])#t^5
            sage: all(c._bracket_(g) == 0 for g in L.lie_algebra_generators())
            True
        """
        if not self or not y:
            return self._parent.zero()

        cdef dict d = {}
        cdef UntwistedAffineLieAlgebraElement rt = <UntwistedAffineLieAlgebraElement>(y)
        c = self._parent.base_ring().zero()
        for tl, gl in self._t_dict.items():
            # d contribution from the left
            if rt._d_coeff:
                if tl in d:
                    d[tl] -= rt._d_coeff * gl * tl
                else:
                    d[tl] = -rt._d_coeff * gl * tl
                if not d[tl]:
                    del d[tl]
            # main bracket of the central extension
            for tr, gr in rt._t_dict.items():
                b = gl.bracket(gr)
                if b:
                    if tl+tr in d:
                        d[tl+tr] += b
                    else:
                        d[tl+tr] = b
                    if not d[tl+tr]:
                        del d[tl+tr]
                if tl + tr == 0:
                    c += gl.killing_form(gr) * tl

        # d contribution from the right
        if self._d_coeff:
            for tr, gr in rt._t_dict.items():
                if tr in d:
                    d[tr] += self._d_coeff * gr * tr
                else:
                    d[tr] = self._d_coeff * gr * tr
                if not d[tr]:
                    del d[tr]

        return type(self)(self._parent, d, c,
                          self._parent.base_ring().zero())

    cpdef canonical_derivation(self):
        r"""
        Return the canonical derivation `d` applied to ``self``.

        The canonical derivation `d` is defined as

        .. MATH::

            d(a \otimes t^m + \alpha c) = a \otimes m t^m.

        Another formulation is by `d = t \frac{d}{dt}`.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['E',6,1])
            sage: al = RootSystem(['E',6]).root_lattice().simple_roots()
            sage: x = L.basis()[al[2]+al[3]+2*al[4]+al[5],5] + 4*L.c() + L.d()
            sage: x.canonical_derivation()
            (5*E[alpha[2] + alpha[3] + 2*alpha[4] + alpha[5]])#t^5
        """
        cdef dict d = {tl: tl * gl
                       for tl, gl in self._t_dict.items() if tl != 0}
        zero = self._parent.base_ring().zero()
        return type(self)(self._parent, d, zero, zero)


def _build_untwisted_affine_element(P, t_dict, c, d):
    """
    Used to unpickle an element.

    EXAMPLES::

        sage: L = lie_algebras.Affine(QQ, ['A',2,1])
        sage: from sage.algebras.lie_algebras.lie_algebra_element import _build_untwisted_affine_element
        sage: _build_untwisted_affine_element(L, {}, 0, 0) == L.zero()
        True
        sage: x = L.an_element()
        sage: loads(dumps(x)) == x  # indirect doctest
        True
    """
    return P.element_class(P, t_dict, c, d)


class FreeLieAlgebraElement(LieAlgebraElement):
    """
    An element of a free Lie algebra.
    """
    # TODO: Move to the category/lift morphism?
    # TODO: Don't override the LieAlgebraElement.lift or should we move
    #    LieAlgebraElement.lift because it is for a specific implementation?
    def lift(self):
        """
        Lift ``self`` to the universal enveloping algebra.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: x,y,z = Lyn.gens()
            sage: a = Lyn([z, [[x, y], x]]); a
            [x, [x, [y, z]]] + [x, [[x, z], y]] - [[x, y], [x, z]]
            sage: a.lift()
            x^2*y*z - 2*x*y*x*z + y*x^2*z - z*x^2*y + 2*z*x*y*x - z*y*x^2
        """
        UEA = self.parent().universal_enveloping_algebra()
        s = UEA.zero()
        if not self:
            return s
        gen_dict = UEA.gens_dict()
        for t, c in self._monomial_coefficients.items():
            s += c * t.lift(gen_dict)
        return s

    def list(self):
        """
        Return ``self`` as a list of pairs ``(m, c)`` where ``m`` is a
        basis key (i.e., a key of one of the basis elements)
        and ``c`` is its coefficient.

        This list is sorted from highest to lowest degree.

        EXAMPLES::

            sage: L.<x, y> = LieAlgebra(QQ)
            sage: elt = x + L.bracket(y, x)
            sage: elt.list()
            [([x, y], -1), (x, 1)]
        """
        def k(x):
            y = x[0]
            return (-y._grade, y) if isinstance(y, GradedLieBracket) else (-1, y)
        return sorted((<dict>self._monomial_coefficients).items(), key=k)

    def _bracket_(self, y):
        """
        Return the Lie bracket ``[self, y]``.

        EXAMPLES::

            sage: L.<x, y> = LieAlgebra(QQ)
            sage: L.bracket(x, y)
            [x, y]
            sage: L.bracket(x, x)
            0
            sage: L.bracket(x, L.bracket(x, y))
            [x, [x, y]]
            sage: L.bracket(y, L.bracket(y, x))
            [[x, y], y]
        """
        if not self or not y:
            return self.parent().zero()

        cdef dict d = {}
        zero = self.base_ring().zero()
        for ml, cl in self._monomial_coefficients.items():  # The left monomials
            for mr, cr in y._monomial_coefficients.items():  # The right monomials
                if ml == mr:
                    continue
                if ml < mr:  # Make sure ml < mr
                    a, b = ml, mr
                else:
                    a, b = mr, ml
                    cr = -cr
                for b_elt, coeff in self.parent()._rewrite_bracket(a, b).items():
                    d[b_elt] = d.get(b_elt, zero) + cl * cr * coeff
                    if d[b_elt] == zero:
                        del d[b_elt]

        if not d:
            return self.parent().zero()
        return type(self)(self.parent(), d)


#####################################################################
# Helper classes for free Lie algebras

cdef class LieObject(SageObject):
    """
    Abstract base class for :class:`LieGenerator` and :class:`LieBracket`.
    """
    cpdef tuple to_word(self):
        """
        Return the word ("flattening") of ``self``.

        If ``self`` is a tree of Lie brackets, this word is
        usually obtained by "forgetting the brackets".

        TESTS::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieObject
            sage: x = LieObject()
            sage: x.to_word()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


cdef class LieGenerator(LieObject):
    """
    A wrapper around an object so it can ducktype with and do
    comparison operations with :class:`LieBracket`.
    """
    def __init__(self, name, index):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: x = LieGenerator('x', 0)
            sage: TestSuite(x).run()
        """
        self._word = (name,)
        self._name = name
        self._index_word = (index,)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: x = LieGenerator('x', 0)
            sage: loads(dumps(x)) == x
            True
        """
        return (LieGenerator, (self._name, self._index_word[0]))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: LieGenerator('x', 0)
            x
        """
        return self._name

    _latex_ = _repr_

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: x = LieGenerator('x', 0)
            sage: hash(x) == hash('x')
            True
        """
        return hash(self._name)

    def __richcmp__(self, rhs, int op):
        """
        Compare equals.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: x == y
            False
            sage: x < y
            True
            sage: y < x
            False
            sage: z = LieGenerator('x', 0)
            sage: x == z
            True
            sage: z = LieBracket(x, y)
            sage: x < z
            True
        """
        if isinstance(rhs, LieBracket):
            if op == Py_NE:
                return True
            if op == Py_EQ:
                return False
            return NotImplemented
            # This causes the comparison to be done ``rhs``
            # (different subclasses of ``LieBracket`` may
            # compare differently with objects in ``LieGenerator``).
            # (Python automatically tries to check ``rhs > self``
            # when the comparison ``self < rhs`` returns a
            # NotImplemented error.)
        if isinstance(rhs, LieGenerator):
            return richcmp(self._index_word[0], <LieGenerator>(rhs)._index_word[0], op)
        return op == Py_NE

    def _im_gens_(self, codomain, im_gens, names):
        """
        Return the image of ``self`` in ``codomain`` under the
        map that sends the generators of the parent of ``self``
        to the elements of the tuple ``im_gens``.
        ``names`` should be the list of all names of the domain
        where ``self`` comes from.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: im = H(Lyn.lie_algebra_generators()['x']); im # indirect doctest
            x
            sage: im.parent() is H
            True
        """
        return im_gens[names.index(self._name)]

    cpdef tuple to_word(self):
        """
        Return the word ("flattening") of ``self``.

        If ``self`` is a tree of Lie brackets, this word is
        usually obtained by "forgetting the brackets".

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: x = LieGenerator('x', 0)
            sage: x.to_word()
            ('x',)
        """
        return self._word

    cpdef lift(self, dict UEA_gens_dict):
        """
        Lift ``self`` to the universal enveloping algebra.

        ``UEA_gens_dict`` should be the dictionary for the
        generators of the universal enveloping algebra.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: x,y,z = Lyn.gens()
            sage: x.lift()
            x
            sage: x.lift().parent()
            Free Algebra on 3 generators (x, y, z) over Rational Field
        """
        return UEA_gens_dict[self._name]

cdef class LieBracket(LieObject):
    """
    An abstract Lie bracket (formally, just a binary tree).
    """
    def __init__(self, LieObject l, LieObject r):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: z = LieBracket(x, y)
            sage: TestSuite(z).run()
        """
        self._left = l
        self._right = r
        self._word = ()
        self._index_word = self._left._index_word + self._right._index_word
        self._hash = -1

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: z = LieBracket(x, y)
            sage: loads(dumps(z)) == z
            True
        """
        return (LieBracket, (self._left, self._right))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: LieBracket(x, y)
            [x, y]
        """
        return "[{!s}, {!s}]".format(self._left, self._right)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: z = LieBracket(x, y)
            sage: latex(z)
            \left[ x , y \right]
        """
        from sage.misc.latex import latex
        return "\\left[" + latex(self._left) + "," + latex(self._right) + "\\right]"

    def __getitem__(self, i):
        r"""
        Return the `i`-th item of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: z = LieBracket(x, y)
            sage: z[0]
            x
            sage: z[1] is y
            True
            sage: z[2]
            Traceback (most recent call last):
            ...
            IndexError: must be either 0 or 1
        """
        if i == 0:
            return self._left
        if i == 1:
            return self._right
        raise IndexError("must be either 0 or 1")

    def __richcmp__(self, rhs, int op):
        """
        Check equality.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: z = LieGenerator('z', 2)
            sage: b = LieBracket(x, y)
            sage: c = LieBracket(y, x)
            sage: b == c
            False
            sage: b == x
            False
            sage: b < x
            False
            sage: b < z
            False
            sage: a = LieBracket(x, y)
            sage: a == b
            True
            sage: a == [x, y]
            True
            sage: c = LieBracket(x, z)
            sage: b < c
            True
        """
        cdef LieBracket right
        if isinstance(rhs, LieBracket):
            right = <LieBracket>(rhs)
            return richcmp([self._left, self._right], [right._left, right._right], op)
        if isinstance(rhs, LieGenerator):
            # Check this is right as in LieGenerator.__richcmp__
            return op == Py_NE or op == Py_GT or op == Py_GE
        if isinstance(rhs, list):
            return richcmp([self._left, self._right], rhs, op)
        return op == Py_NE

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: z = LieGenerator('z', 2)
            sage: b = LieBracket(x, y)
            sage: hash(b) == hash(b)
            True
        """
        if self._hash == -1:
            self._hash = hash((self._left, self._right))
        return self._hash

    def _im_gens_(self, codomain, im_gens, names):
        """
        Return the image of ``self`` in ``codomain`` under the
        map that sends the generators of the parent of ``self``
        to the elements of the tuple ``im_gens``.
        ``names`` should be the list of all names of the domain
        where ``self`` comes from.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: x,y,z = Lyn.gens()
            sage: im = H(Lyn([z, [[x, y], x]])); im # indirect doctest
            -[z, [x, [x, y]]]
            sage: im.parent() is H
            True
        """
        return codomain.bracket(self._left._im_gens_(codomain, im_gens, names),
                                self._right._im_gens_(codomain, im_gens, names))

    cpdef lift(self, dict UEA_gens_dict):
        """
        Lift ``self`` to the universal enveloping algebra.

        ``UEA_gens_dict`` should be the dictionary for the
        generators of the universal enveloping algebra.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: x,y,z = Lyn.gens()
            sage: a = Lyn([z, [[x, y], x]]); a
            [x, [x, [y, z]]] + [x, [[x, z], y]] - [[x, y], [x, z]]
            sage: a.lift() # indirect doctest
            x^2*y*z - 2*x*y*x*z + y*x^2*z - z*x^2*y + 2*z*x*y*x - z*y*x^2
        """
        if isinstance(self._left, LieBracket):
            l = self._left.lift(UEA_gens_dict)
        else:
            l = UEA_gens_dict[self._left._name]

        if isinstance(self._right, LieBracket):
            r = self._right.lift(UEA_gens_dict)
        else:
            r = UEA_gens_dict[self._right._name]

        return l*r - r*l

    cpdef tuple to_word(self):
        """
        Return the word ("flattening") of ``self``.

        If ``self`` is a tree of Lie brackets, this word is
        usually obtained by "forgetting the brackets".

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: b = LieBracket(x, y)
            sage: c = LieBracket(b, x)
            sage: c.to_word()
            ('x', 'y', 'x')
        """
        if not self._word:
            self._word = self._left.to_word() + self._right.to_word()
        return self._word

cdef class GradedLieBracket(LieBracket):
    """
    A Lie bracket (:class:`LieBracket`) for a graded Lie algebra.

    Unlike the vanilla Lie bracket class, this also stores a
    degree, and uses it as a first criterion when comparing
    graded Lie brackets.
    (Graded Lie brackets still compare greater than Lie
    generators.)
    """
    def __init__(self, LieObject l, LieObject r, grade):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, GradedLieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: b = GradedLieBracket(x, y, 2)
            sage: TestSuite(b).run()
        """
        self._grade = grade
        LieBracket.__init__(self, l, r)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, GradedLieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: b = GradedLieBracket(x, y, 2)
            sage: loads(dumps(b)) == b
            True
        """
        return (type(self), (self._left, self._right, self._grade))

    def __richcmp__(self, rhs, int op):
        """
        Check less than.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, GradedLieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: z = LieGenerator('z', 2)
            sage: b = GradedLieBracket(x, y, 2)
            sage: b < x
            False
            sage: b < z
            False
            sage: c = GradedLieBracket(x, z, 2)
            sage: b < c
            True
            sage: c = GradedLieBracket(x, z, 1)
            sage: b < c
            False
        """
        cdef GradedLieBracket right
        if isinstance(rhs, GradedLieBracket):
            right = <GradedLieBracket>(rhs)
            if self._grade != right._grade:
                return richcmp_not_equal(self._grade, right._grade, op)
            return richcmp([self._left, self._right], [right._left, right._right], op)
        if isinstance(rhs, LieGenerator):
            return op == Py_NE or op == Py_GT or op == Py_GE
        return op == Py_NE

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, GradedLieBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: z = LieGenerator('z', 2)
            sage: b = GradedLieBracket(x, y, 2)
            sage: hash(b) == hash(b)
            True
        """
        if self._hash == -1:
            self._hash = hash((self._grade, self._left, self._right))
        return self._hash

cdef class LyndonBracket(GradedLieBracket):
    """
    A Lie bracket (:class:`LieBracket`) tailored for the Lyndon
    basis.

    The order on these brackets is defined by `l < r`
    if `w(l) < w(r)`, where `w(l)` is the word corresponding
    to `l`.
    (This is also true if one or both of `l` and `r` is a
    :class:`LieGenerator`.)
    """
    def __richcmp__(self, rhs, op):
        """
        Compare ``self`` and ``rhs``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LyndonBracket
            sage: x,y,z = [LieGenerator(letter, ind) for letter,ind in zip(['x', 'y', 'z'],range(3))]
            sage: LyndonBracket(x, LyndonBracket(y, z, 2), 3) < LyndonBracket(LyndonBracket(y, z, 2), x, 3)
            True
        """
        if not isinstance(rhs, LieObject):
            return op == Py_NE
        return richcmp(self._index_word, <LieObject>(rhs)._index_word, op)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LyndonBracket
            sage: x = LieGenerator('x', 0)
            sage: y = LieGenerator('y', 1)
            sage: b = LyndonBracket(x, y, 2)
            sage: hash(b) == hash((0, 1))
            True
        """
        if self._hash == -1:
            self._hash = hash(self._index_word)
        return self._hash
