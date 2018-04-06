# -*- coding: utf-8 -*-
"""
Lie Algebra Elements

AUTHORS:

- Travis Scrimshaw (2013-05-04): Initial implementation
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from cpython.object cimport Py_EQ, Py_NE

from sage.misc.misc import repr_lincomb
from sage.combinat.free_module import CombinatorialFreeModule
from sage.structure.element cimport have_same_parent, coercion_model, parent
from sage.cpython.wrapperdescr cimport wrapperdescr_fastcall
from sage.structure.element_wrapper cimport ElementWrapper
from sage.structure.richcmp cimport richcmp
from sage.data_structures.blas_dict cimport axpy, add, negate, scal

# TODO: Do we want a dense version?
cdef class LieAlgebraElement(IndexedFreeModuleElement):
    """
    A Lie algebra element.
    """
    # Need to bypass the coercion model
    def __mul__(left, right):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: y*x
            x*y - z

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

    #def _im_gens_(self, codomain, im_gens):
    #    """
    #    Return the image of ``self`` in ``codomain`` under the map that sends
    #    the images of the generators of the parent of ``self`` to the
    #    tuple of elements of ``im_gens``.
    #
    #    EXAMPLES::
    #    """
    #    s = codomain.zero()
    #    if not self: # If we are 0
    #        return s
    #    names = self.parent().variable_names()
    #    return codomain.sum(c * t._im_gens_(codomain, im_gens, names)
    #                        for t, c in self._monomial_coefficients.iteritems())

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
            b2*b3 - b0

            sage: L = lie_algebras.regular_vector_fields(QQ)
            sage: L.an_element()
            d[-1] + d[0] - 3*d[1]
            sage: L.an_element().lift()
            PBW[-1] + PBW[0] - 3*PBW[1]
        """
        UEA = self._parent.universal_enveloping_algebra()
        try:
            gen_dict = UEA.algebra_generators()
        except (TypeError, AttributeError):
            gen_dict = UEA.gens_dict()
        s = UEA.zero()
        if not self:
            return s
        # Special hook for when the index set of the parent of ``self``
        #   does not match the generators index set of the UEA.
        if hasattr(self._parent, '_UEA_names_map'):
            names_map = self._parent._UEA_names_map
            for t, c in self._monomial_coefficients.iteritems():
                s += c * gen_dict[names_map[t]]
        else:
            for t, c in self._monomial_coefficients.iteritems():
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
    """

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x + y
            x + y
        """
        return repr(self.value)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x')
            sage: L.<x0,x1,x2> = LieAlgebra(associative=R.gens())
            sage: latex(x0 + x1)
            x_{0} + x_{1}
        """
        from sage.misc.latex import latex
        return latex(self.value)

    def _ascii_art_(self):
        """
        Return an ascii art representation of ``self``.

        EXAMPLES::

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
        """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art(self.value)

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: L = LieAlgebra(associative=s)
            sage: P = Partition([4,2,2,1])
            sage: x = L.basis()[P]
            sage: unicode_art(x)
            s
             ┌┬┬┬┐
             ├┼┼┴┘
             ├┼┤
             ├┼┘
             └┘
        """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art(self.value)

    def __nonzero__(self):
        """
        Return if ``self`` is non-zero.

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
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        .. TODO::

            Write more tests for this method.

        EXAMPLES::

            sage: S = SymmetricGroup(3).algebra(QQ)
            sage: L = LieAlgebra(associative=S)
            sage: x = L.gen(2); x
            (1,2,3)
            sage: y = L.gen(1); y
            (1,2)
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
            (2,3) - (1,3)
            sage: G = list(S.basis())
            sage: G[4] - G[5]
            (2,3) - (1,3)

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

    def __div__(self, x):
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
                scalar = self._parent.base_ring()( scalar )
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
            sage: list(x)
            [((2,3), 1), ((1,2), 1), ((1,3), 1),
             ((1,2,3), 1), ((1,3,2), 1), ((), 2)]
        """
        return self.value.monomial_coefficients(copy=False).iteritems()

# TODO: Also used for vectors, find a better name
cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    """
    Lie algebra element wrapper around a matrix.
    """
    def __init__(self, parent, value):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
            sage: z = L.z()
            sage: z.value.is_immutable()
            True
        """
        value.set_immutable() # Make the matrix immutable for hashing
        LieAlgebraElementWrapper.__init__(self, parent, value)

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
            x - \frac{3}{2}y
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self._parent._print_options['scalar_mult'],
                            latex_scalar_mult=self._parent._print_options['latex_scalar_mult'],
                            repr_monomial=self._parent._latex_term,
                            is_latex=True, strip_one=True)

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
        for i,v in enumerate(self.value):
            if v != zero:
                yield (I[i], v)

    cpdef to_vector(self):
        """
        Return ``self`` as a vector.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: a = x + 3*y - z/2
            sage: a.to_vector()
            (1, 3, -1/2)
        """
        return self.value

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
        return UEA.sum(c * gens[i] for i, c in self.value.iteritems())

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
        return {I[i]: v for i,v in self.value.iteritems()}

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
        ret = ' + '.join('({})#t^{}'.format(g, t)
                         for t,g in self._t_dict.iteritems())
        if self._c_coeff != 0:
            if ret:
                ret += ' + '
            if self._c_coeff != 1:
                ret += repr(self._c_coeff) + '*c'
            else:
                ret += 'c'

        if self._d_coeff != 0:
            if ret:
                ret += ' + '
            if self._d_coeff != 1:
                ret += repr(self._d_coeff) + '*d'
            else:
                ret += 'd'

        if not ret:
            return '0'
        return ret

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: [latex(g) for g in L.lie_algebra_generators()]
            [(E_{\alpha_{1}}) \otimes t^{0},
             (E_{-\alpha_{1}}) \otimes t^{0},
             (E_{\alpha^\vee_{1}}) \otimes t^{0},
             (E_{-\alpha_{1}}) \otimes t^{1},
             (E_{\alpha_{1}}) \otimes t^{-1},
             c,
             d]
            sage: latex(L.an_element())
            (E_{\alpha_{1}} + E_{\alpha^\vee_{1}} + E_{-\alpha_{1}}) \otimes t^{0}
             + (E_{-\alpha_{1}}) \otimes t^{1} + (E_{\alpha_{1}}) \otimes t^{-1}
             + c + d
            sage: latex(L.zero())
            0

            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: latex(e1 + 2*f1 - h1 + e0 + 3*c - 2*d)
            (E_{\alpha_{1}} - E_{\alpha^\vee_{1}} + 2E_{-\alpha_{1}}) \otimes t^{0}
             + (E_{-\alpha_{1}}) \otimes t^{1} + 3 c + -2 d
        """
        from sage.misc.latex import latex
        ret = ' + '.join('({}) \otimes t^{{{}}}'.format(latex(g), t)
                         for t,g in self._t_dict.iteritems())
        if self._c_coeff != 0:
            if ret:
                ret += ' + '
            if self._c_coeff != 1:
                ret += latex(self._c_coeff) + ' c'
            else:
                ret += 'c'

        if self._d_coeff != 0:
            if ret:
                ret += ' + '
            if self._d_coeff != 1:
                ret += latex(self._d_coeff) + ' d'
            else:
                ret += 'd'

        if not ret:
            return '0'
        return ret

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

    def __nonzero__(self):
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
            sage: e0.bracket(e1) + d + e1 + c + 3*d
            (E[alpha[1]])#t^0 + (-h1)#t^1 + c + 4*d
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
            sage: e0.bracket(e1) + d - e1 + c - 3*d
            (-E[alpha[1]])#t^0 + (-h1)#t^1 + c + -2*d
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
            sage: -x
            (-E[alpha[1]])#t^0 + (h1)#t^1 + -1*c + -4*d
        """
        return type(self)(self._parent, negate(self._t_dict),
                          -self._c_coeff, -self._d_coeff)

    cpdef _acted_upon_(self, x, bint self_on_left):
        """
        Return ``self`` acted upon by ``x``.

        EXAMPLES::

            sage: L = lie_algebras.Affine(QQ, ['A',1,1])
            sage: e1,f1,h1,e0,f0,c,d = list(L.lie_algebra_generators())
            sage: x = e1 + f0.bracket(f1) + 3*c - 2/5 * d
            sage: x
            (E[alpha[1]])#t^0 + (h1)#t^-1 + 3*c + -2/5*d
            sage: -2 * x
            (-2*E[alpha[1]])#t^0 + (-2*h1)#t^-1 + -6*c + 4/5*d
        """
        return type(self)(self._parent, scal(x, self._t_dict, self_on_left),
                          x * self._c_coeff,
                          x * self._d_coeff)

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
        for t,g in self._t_dict.iteritems():
            for k,c in g.monomial_coefficients(copy=False).iteritems():
                d[k,t] = c
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

        gd = self._parent._g.basis()
        cdef dict d = {}
        cdef UntwistedAffineLieAlgebraElement rt = <UntwistedAffineLieAlgebraElement>(y)
        c = self._parent.base_ring().zero()
        for tl,gl in self._t_dict.iteritems():
            # d contribution from the left
            if rt._d_coeff:
                if tl in d:
                    d[tl] -= rt._d_coeff * gl * tl
                else:
                    d[tl] = -rt._d_coeff * gl * tl
                if not d[tl]:
                    del d[tl]
            # main bracket of the central extension
            for tr,gr in rt._t_dict.iteritems():
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
            for tr,gr in rt._t_dict.iteritems():
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
        cdef dict d = {tl: tl * gl for tl,gl in self._t_dict.iteritems() if tl != 0}
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

