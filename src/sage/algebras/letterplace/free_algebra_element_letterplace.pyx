"""
Weighted homogeneous elements of free algebras, in letterplace implementation

AUTHOR:

- Simon King (2011-03-23): Github issue :issue:`7797`
"""

# ****************************************************************************
#       Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.perm_gps.permgroup_named import CyclicPermutationGroup
from sage.libs.singular.function import lib, singular_function
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from cpython.object cimport PyObject_RichCompare

# Define some singular functions
lib("freegb.lib")
poly_reduce = singular_function("NF")

#####################
# Free algebra elements
cdef class FreeAlgebraElement_letterplace(AlgebraElement):
    """
    Weighted homogeneous elements of a free associative unital algebra
    (letterplace implementation).

    EXAMPLES::

        sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
        sage: x+y
        x + y
        sage: x*y !=y*x
        True
        sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
        sage: (y^3).reduce(I)
        y*y*y
        sage: (y^3).normal_form(I)
        y*y*z - y*z*y + y*z*z

    Here is an example with nontrivial degree weights::

        sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
        sage: I = F*[x*y-y*x, x^2+2*y*z, (x*y)^2-z^2]*F
        sage: x.degree()
        2
        sage: y.degree()
        1
        sage: z.degree()
        3
        sage: (x*y)^3
        x*y*x*y*x*y
        sage: ((x*y)^3).normal_form(I)
        z*z*y*x
        sage: ((x*y)^3).degree()
        9
    """
    def __init__(self, A, x, check=True):
        """
        INPUT:

        - ``A`` - a free associative unital algebra in letterplace implementation
        - ``x`` -- a homogeneous polynomial that can be coerced into the currently
          used polynomial ring of `A`
        - ``check`` -- boolean (default: ``True``); do not attempt the
          above coercion (for internal use only)

        TESTS::

            sage: from sage.algebras.letterplace.free_algebra_element_letterplace import FreeAlgebraElement_letterplace
            sage: F.<x,y,z> = FreeAlgebra(GF(3), implementation='letterplace')
            sage: F.set_degbound(2)
            sage: P = F.current_ring()
            sage: F.set_degbound(4)
            sage: P == F.current_ring()
            False
            sage: p = FreeAlgebraElement_letterplace(F,P.1*P.3+2*P.0*P.4); p
            -x*y + y*x
            sage: loads(dumps(p)) == p
            True
        """
        cdef FreeAlgebra_letterplace P = A
        if check:
            if not x.is_homogeneous():
                raise ValueError("free algebras based on Letterplace can currently only work with weighted homogeneous elements")
            P.set_degbound(x.degree())
            x = P._current_ring(x)
        AlgebraElement.__init__(self, P)
        self._poly = x

    def __reduce__(self):
        """
        Pickling.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: loads(dumps(x*y*x)) == x*y*x   # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._poly)

    def __copy__(self):
        """
        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: copy(x*y*z+z*y*x) == x*y*z+z*y*x   # indirect doctest
            True
        """
        self._poly = (<FreeAlgebra_letterplace>self._parent)._current_ring(self._poly)
        return self.__class__(self._parent, self._poly, check=False)

    def __hash__(self):
        """
        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: set([x*y*z, z*y+x*z,x*y*z])  # indirect doctest
            {x*z + z*y, x*y*z}
        """
        return hash(self._poly)

    def __iter__(self):
        """
        Iterate over the pairs "tuple of exponents, coefficient".

        EXAMPLES::

            sage: F.<w,x,y,z> = FreeAlgebra(GF(3), implementation='letterplace')
            sage: p = x*y-z^2
            sage: sorted(p)   # indirect doctest
            [((0, 0, 0, 1, 0, 0, 0, 1), 2), ((0, 1, 0, 0, 0, 0, 1, 0), 1)]
        """
        cdef dict d = self._poly.monomial_coefficients()
        yield from d.items()

    def _repr_(self):
        """
        TESTS::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
            sage: -(a+b*(z+1)-c)^2     # indirect doctest
            -a*a + (4*z + 4)*a*b + a*c + (4*z + 4)*b*a + (2*z + 1)*b*b + (z + 1)*b*c + c*a + (z + 1)*c*b - c*c

        It is possible to change the names temporarily::

            sage: from sage.structure.parent_gens import localvars
            sage: with localvars(F, ['w', 'x','y']):
            ....:     print(a+b*(z+1)-c)
            w + (z + 1)*x - y
            sage: print(a+b*(z+1)-c)
            a + (z + 1)*b - c
        """
        cdef list L = []
        cdef FreeAlgebra_letterplace P = self._parent
        if P._base._repr_option('element_is_atomic'):
            for E, c in zip(self._poly.exponents(), self._poly.coefficients()):
                monstr = P.exponents_to_string(E)
                if monstr:
                    if c == 1:
                        if L:
                            L.extend(['+', monstr])
                        else:
                            L.append(monstr)
                    elif c == -1:
                        if L:
                            L.extend(['-', monstr])
                        else:
                            L.append('-' + monstr)
                    else:
                        if L:
                            if c >= 0:
                                L.extend(['+', repr(c) + '*' + monstr])
                            else:
                                L.extend(['-', repr(-c) + '*' + monstr])
                        else:
                            L.append(repr(c) + '*' + monstr)
                else:
                    if c >= 0:
                        if L:
                            L.extend(['+', repr(c)])
                        else:
                            L.append(repr(c))
                    else:
                        if L:
                            L.extend(['-', repr(-c)])
                        else:
                            L.append(repr(c))
        else:
            for E, c in zip(self._poly.exponents(), self._poly.coefficients()):
                monstr = P.exponents_to_string(E)
                if monstr:
                    if c == 1:
                        if L:
                            L.extend(['+', monstr])
                        else:
                            L.append(monstr)
                    elif c == -1:
                        if L:
                            L.extend(['-', monstr])
                        else:
                            L.append('-' + monstr)
                    else:
                        if L:
                            L.extend(['+', '(' + repr(c) + ')*' + monstr])
                        else:
                            L.append('(' + repr(c) + ')*' + monstr)
                else:
                    if L:
                        L.extend(['+', repr(c)])
                    else:
                        L.append(repr(c))
        if L:
            return ' '.join(L)
        return '0'

    def _latex_(self):
        r"""
        TESTS::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace', degrees=[1,2,3])
            sage: -(a*b*(z+1)-c)^2
            (2*z + 1)*a*b*a*b + (z + 1)*a*b*c + (z + 1)*c*a*b - c*c
            sage: latex(-(a*b*(z+1)-c)^2)     # indirect doctest
            \left(2 z + 1\right) a b a b + \left(z + 1\right) a b c + \left(z + 1\right) c a b - c c
        """
        cdef list L = []
        cdef FreeAlgebra_letterplace P = self._parent
        from sage.misc.latex import latex
        if P._base._repr_option('element_is_atomic'):
            for E, c in zip(self._poly.exponents(), self._poly.coefficients()):
                monstr = P.exponents_to_latex(E)
                if monstr:
                    if c == 1:
                        if L:
                            L.extend(['+', monstr])
                        else:
                            L.append(monstr)
                    elif c == -1:
                        if L:
                            L.extend(['-', monstr])
                        else:
                            L.append('-' + monstr)
                    else:
                        if L:
                            if c >= 0:
                                L.extend(['+', repr(latex(c)) + ' ' + monstr])
                            else:
                                L.extend(['-', repr(latex(-c)) + ' ' + monstr])
                        else:
                            L.append(repr(latex(c)) + ' ' + monstr)
                else:
                    if c >= 0:
                        if L:
                            L.extend(['+', repr(latex(c))])
                        else:
                            L.append(repr(latex(c)))
                    else:
                        if L:
                            L.extend(['-', repr(latex(-c))])
                        else:
                            L.append(repr(c))
        else:
            for E, c in zip(self._poly.exponents(), self._poly.coefficients()):
                monstr = P.exponents_to_latex(E)
                if monstr:
                    if c == 1:
                        if L:
                            L.extend(['+', monstr])
                        else:
                            L.append(monstr)
                    elif c == -1:
                        if L:
                            L.extend(['-', monstr])
                        else:
                            L.append('-' + monstr)
                    else:
                        if L:
                            L.extend(['+', '\\left(' + repr(latex(c)) + '\\right) ' + monstr])
                        else:
                            L.append('\\left(' + repr(latex(c)) + '\\right) ' + monstr)
                else:
                    if L:
                        L.extend(['+', repr(latex(c))])
                    else:
                        L.append(repr(latex(c)))
        if L:
            return ' '.join(L)
        return '0'

    def degree(self):
        """
        Return the degree of this element.

        .. NOTE::

            Generators may have a positive integral degree weight. All
            elements must be weighted homogeneous.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: ((x+y+z)^3).degree()
            3
            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: ((x*y+z)^3).degree()
            9
        """
        return self._poly.degree()

    def letterplace_polynomial(self):
        """
        Return the commutative polynomial that is used internally to represent this free algebra element.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: ((x+y-z)^2).letterplace_polynomial()
            x*x_1 + x*y_1 - x*z_1 + y*x_1 + y*y_1 - y*z_1 - z*x_1 - z*y_1 + z*z_1

        If degree weights are used, the letterplace polynomial is
        homogenized by slack variables::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: ((x*y+z)^2).letterplace_polynomial()
            x*x__1*y_2*x_3*x__4*y_5 + x*x__1*y_2*z_3*x__4*x__5 + z*x__1*x__2*x_3*x__4*y_5 + z*x__1*x__2*z_3*x__4*x__5
        """
        return self._poly

    def lm(self):
        """
        The leading monomial of this free algebra element.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: ((2*x+3*y-4*z)^2*(5*y+6*z)).lm()
            x*x*y
            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: ((2*x*y+z)^2).lm()
            x*y*x*y
        """
        return FreeAlgebraElement_letterplace(self._parent, self._poly.lm())

    def lt(self):
        """
        The leading term (monomial times coefficient) of this free algebra
        element.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: ((2*x+3*y-4*z)^2*(5*y+6*z)).lt()
            20*x*x*y
            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: ((2*x*y+z)^2).lt()
            4*x*y*x*y
        """
        return FreeAlgebraElement_letterplace(self._parent, self._poly.lt())

    def lc(self):
        """
        The leading coefficient of this free algebra element, as element
        of the base ring.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: ((2*x+3*y-4*z)^2*(5*y+6*z)).lc()
            20
            sage: ((2*x+3*y-4*z)^2*(5*y+6*z)).lc().parent() is F.base()
            True
            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: ((2*x*y+z)^2).lc()
            4
        """
        return self._poly.lc()

    def __bool__(self):
        """
        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: bool(x)      # indirect doctest
            True
            sage: bool(F.zero())
            False
        """
        return bool(self._poly)

    def lm_divides(self, FreeAlgebraElement_letterplace p):
        """
        Tell whether or not the leading monomial of ``self`` divides the
        leading monomial of another element.

        .. NOTE::

            A free algebra element `p` divides another one `q` if there are
            free algebra elements `s` and `t` such that `spt = q`.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: ((2*x*y+z)^2*z).lm()
            x*y*x*y*z
            sage: (y*x*y-y^4).lm()
            y*x*y
            sage: (y*x*y-y^4).lm_divides((2*x*y+z)^2*z)
            True
        """
        if self._parent is not p._parent:
            raise TypeError("the two arguments must be elements in the same free algebra")
        cdef FreeAlgebra_letterplace A = self._parent
        P = A._current_ring
        p_poly = p._poly = P(p._poly)
        s_poly = self._poly = P(self._poly)
        cdef int p_d = p_poly.degree()
        cdef int s_d = s_poly.degree()
        if s_d > p_d:
            return False
        cdef int i
        if P.monomial_divides(s_poly, p_poly):
            return True
        realngens = A._commutative_ring.ngens()
        CG = CyclicPermutationGroup(P.ngens())
        for i in range(p_d - s_d):
            s_poly = s_poly * CG[realngens]
            if P.monomial_divides(s_poly, p_poly):
                return True
        return False

    cpdef _richcmp_(self, other, int op):
        """
        Implement comparisons, using the Cython richcmp convention.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: p = ((2*x+3*y-4*z)^2*(5*y+6*z))
            sage: p - p.lt() < p     # indirect doctest
            True
        """
        left = (<FreeAlgebraElement_letterplace>self)._poly
        right = (<FreeAlgebraElement_letterplace>other)._poly
        return PyObject_RichCompare(left, right, op)

    ################################
    # Arithmetic
    cpdef _neg_(self):
        """
        TESTS::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
            sage: -((z+2)*a^2*b+3*c^3)  # indirect doctest
            (4*z + 3)*a*a*b + (2)*c*c*c
            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: -(3*x*y+2*z^2)
            -3*x*y - 2*z*z
        """
        return FreeAlgebraElement_letterplace(self._parent, -self._poly,
                                              check=False)

    cpdef _add_(self, other):
        """
        Addition, under the side condition that either one summand
        is zero, or both summands have the same degree.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: x+y    # indirect doctest
            x + y
            sage: x+1
            Traceback (most recent call last):
            ...
            ArithmeticError: can only add elements of the same weighted degree
            sage: x+0
            x
            sage: 0+x
            x
        """
        if not other:
            return self
        if not self:
            return other
        cdef FreeAlgebraElement_letterplace right = other
        if right._poly.degree() != self._poly.degree():
            raise ArithmeticError("can only add elements of the same weighted degree")
        # update the polynomials
        cdef FreeAlgebra_letterplace A = self._parent
        self._poly = A._current_ring(self._poly)
        right._poly = A._current_ring(right._poly)
        return FreeAlgebraElement_letterplace(self._parent,
                                              self._poly + right._poly,
                                              check=False)

    cpdef _sub_(self, other):
        """
        Difference, under the side condition that either one summand
        is zero or both have the same weighted degree.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: x*y-y*x     # indirect doctest
            x*y - y*x
            sage: x-1
            Traceback (most recent call last):
            ...
            ArithmeticError: can only subtract elements of the same degree
            sage: x-0
            x
            sage: 0-x
            -x

        Here is an example with non-trivial degree weights::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: x*y+z
            x*y + z
        """
        if not other:
            return self
        if not self:
            return -other
        cdef FreeAlgebraElement_letterplace right = other
        if right._poly.degree() != self._poly.degree():
            raise ArithmeticError("can only subtract elements of the same degree")
        # update the polynomials
        cdef FreeAlgebra_letterplace A = self._parent
        self._poly = A._current_ring(self._poly)
        right._poly = A._current_ring(right._poly)
        return FreeAlgebraElement_letterplace(self._parent,
                                              self._poly - right._poly,
                                              check=False)

    cpdef _lmul_(self, Element right):
        """
        Multiplication from the right with an element of the base ring.

        TESTS::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
            sage: (a+b)*(z+1)    # indirect doctest
            (z + 1)*a + (z + 1)*b
        """
        return FreeAlgebraElement_letterplace(self._parent,
                                              self._poly._lmul_(right),
                                              check=False)

    cpdef _rmul_(self, Element left):
        """
        Multiplication from the left with an element of the base ring.

        TESTS::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
            sage: (z+1)*(a+b)   # indirect doctest
            (z + 1)*a + (z + 1)*b
        """
        return FreeAlgebraElement_letterplace(self._parent,
                                              self._poly._rmul_(left),
                                              check=False)

    cpdef _mul_(self, other):
        """
        Product of two free algebra elements in letterplace implementation.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: (x*y+z)*z   # indirect doctest
            x*y*z + z*z
        """
        cdef FreeAlgebraElement_letterplace left = self
        cdef FreeAlgebraElement_letterplace right = other
        cdef FreeAlgebra_letterplace A = left._parent
        A.set_degbound(left._poly.degree() + right._poly.degree())
        # we must put the polynomials into the same ring
        left._poly = A._current_ring(left._poly)
        right._poly = A._current_ring(right._poly)
        realngens = A._commutative_ring.ngens()
        CG = CyclicPermutationGroup(A._current_ring.ngens())
        rshift = right._poly * CG[left._poly.degree() * realngens]
        return FreeAlgebraElement_letterplace(A, left._poly * rshift,
                                              check=False)

    def __pow__(FreeAlgebraElement_letterplace self, int n, k):
        """
        TESTS::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
            sage: (a+z*b)^3    # indirect doctest
            a*a*a + (z)*a*a*b + (z)*a*b*a + (z + 3)*a*b*b + (z)*b*a*a + (z + 3)*b*a*b + (z + 3)*b*b*a + (4*z + 3)*b*b*b
        """
        cdef FreeAlgebra_letterplace A = self._parent
        if n < 0:
            raise ValueError("negative exponents are not allowed")
        if n == 0:
            return FreeAlgebraElement_letterplace(A, A._current_ring(1),
                                                  check=False)
        if n == 1:
            return self
        A.set_degbound(self._poly.degree() * n)
        cdef MPolynomial_libsingular p, q
        self._poly = A._current_ring(self._poly)
        cdef int d = self._poly.degree()
        q = p = self._poly
        realngens = A._commutative_ring.ngens()
        cdef int i
        CG = CyclicPermutationGroup(A._current_ring.ngens())
        for i in range(1, n):
            q = q * CG[d * realngens]
            p *= q
        return FreeAlgebraElement_letterplace(A, p, check=False)

    # Groebner related stuff
    def reduce(self, G):
        """
        Reduce this element by a list of elements or by a
        twosided weighted homogeneous ideal.

        INPUT:

        Either a list or tuple of weighted homogeneous elements of the
        free algebra, or an ideal of the free algebra, or an ideal in
        the commutative polynomial ring that is currently used to
        implement the multiplication in the free algebra.

        OUTPUT: the twosided reduction of this element by the argument

        .. NOTE::

            This may not be the normal form of this element, unless
            the argument is a twosided Groebner basis up to the degree
            of this element.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: p = y^2*z*y^2+y*z*y*z*y

        We compute the letterplace version of the Groebner basis
        of `I` with degree bound 4::

            sage: G = F._reductor_(I.groebner_basis(4).gens(),4)
            sage: G.ring() is F.current_ring()
            True

        Since the element `p` is of degree 5, it is no surprise
        that its reductions with respect to the original generators
        of `I` (of degree 2), or with respect to `G` (Groebner basis
        with degree bound 4), or with respect to the Groebner basis
        with degree bound 5 (which yields its normal form) are
        pairwise different::

            sage: p.reduce(I)
            y*y*z*y*y + y*z*y*z*y
            sage: p.reduce(G)
            y*y*z*z*y + y*z*y*z*y - y*z*z*y*y + y*z*z*z*y
            sage: p.normal_form(I)
            y*y*z*z*z + y*z*y*z*z - y*z*z*y*z + y*z*z*z*z
            sage: p.reduce(I) != p.reduce(G) != p.normal_form(I) != p.reduce(I)
            True
        """
        cdef FreeAlgebra_letterplace P = self._parent
        if not isinstance(G, (list, tuple)):
            if G == P:
                return P.zero()
            if not (isinstance(G, MPolynomialIdeal) and G.ring() == P._current_ring):
                G = G.gens()
        C = P.current_ring()
        cdef int selfdeg = self._poly.degree()
        if isinstance(G, MPolynomialIdeal):
            gI = G
        else:
            gI = P._reductor_(G, selfdeg)  # C.ideal(g,coerce=False)
        from sage.libs.singular.option import LibSingularOptions
        libsingular_options = LibSingularOptions()
        bck = (libsingular_options['redTail'], libsingular_options['redSB'])
        libsingular_options['redTail'] = True
        libsingular_options['redSB'] = True
        poly = poly_reduce(C(self._poly), gI, ring=C,
                           attributes={gI: {"isSB": 1}})
        libsingular_options['redTail'] = bck[0]
        libsingular_options['redSB'] = bck[1]
        return FreeAlgebraElement_letterplace(P, poly, check=False)

    def normal_form(self, I):
        """
        Return the normal form of this element with respect to
        a twosided weighted homogeneous ideal.

        INPUT:

        - ``I`` -- a twosided homogeneous ideal of the parent `F` of
          this element, `x`

        OUTPUT: the normal form of `x` wrt. `I`

        .. NOTE::

            The normal form is computed by reduction with respect
            to a Groebnerbasis of `I` with degree bound `deg(x)`.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: (x^5).normal_form(I)
            -y*z*z*z*x - y*z*z*z*y - y*z*z*z*z

        We verify two basic properties of normal forms: The
        difference of an element and its normal form is contained
        in the ideal, and if two elements of the free algebra
        differ by an element of the ideal then they have the same
        normal form::

            sage: x^5 - (x^5).normal_form(I) in I
            True
            sage: (x^5+x*I.0*y*z-3*z^2*I.1*y).normal_form(I) == (x^5).normal_form(I)
            True

        Here is an example with non-trivial degree weights::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[1,2,3])
            sage: I = F*[x*y-y*x+z, y^2+2*x*z, (x*y)^2-z^2]*F
            sage: ((x*y)^3).normal_form(I)
            z*z*y*x - z*z*z
            sage: (x*y)^3-((x*y)^3).normal_form(I) in I
            True
            sage: ((x*y)^3+2*z*I.0*z+y*I.1*z-x*I.2*y).normal_form(I) == ((x*y)^3).normal_form(I)
            True
        """
        if self._parent != I.ring():
            raise ValueError("cannot compute normal form wrt an ideal that does not belong to %s" % self._parent)
        sdeg = self._poly.degree()
        return self.reduce(self._parent._reductor_(I.groebner_basis(degbound=sdeg).gens(), sdeg))
