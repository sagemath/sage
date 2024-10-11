r"""
Elements of multivariate Laurent polynomial rings
"""
# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer cimport Integer
from sage.structure.element cimport CommutativeAlgebraElement, Element, ModuleElement, RingElement
from sage.structure.element import coerce_binop, parent
from sage.structure.factorization import Factorization
from sage.misc.derivative import multi_derivative
from sage.rings.polynomial.polydict cimport monomial_exponent
from sage.matrix.matrix0 cimport Matrix
from sage.rings.infinity import Infinity, minus_infinity


cdef class LaurentPolynomial_mpair(LaurentPolynomial):
    """
    Multivariate Laurent polynomials.
    """

    def __init__(self, parent, x, mon=None, reduce=True):
        """
        Currently, one can only create LaurentPolynomials out of dictionaries
        and elements of the base ring.

        INPUT:

        - ``parent`` -- a SageMath parent

        - ``x`` -- an element or dictionary or anything the underlying
          polynomial ring accepts

        - ``mon`` -- (default: ``None``) a tuple specifying the shift
          in the exponents

        - ``reduce`` -- boolean (default: ``True``)

        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: f = L({(-1,-1):1}); f
            w^-1*z^-1
            sage: f = L({(1,1):1}); f
            w*z
            sage: f =  L({(-1,-1):1, (1,3):4}); f
            4*w*z^3 + w^-1*z^-1
            sage: L(1/2)
            1/2

        TESTS:

        Check that :issue:`19538` is fixed::

            sage: R = LaurentPolynomialRing(QQ,'x2,x0')
            sage: S = LaurentPolynomialRing(QQ,'x',3)
            sage: f = S.coerce_map_from(R)
            sage: f(R.gen(0) + R.gen(1)^2)
            x0^2 + x2
            sage: _.parent()
            Multivariate Laurent Polynomial Ring in x0, x1, x2 over Rational Field

        ::

            sage: from sage.rings.polynomial.laurent_polynomial_mpair import LaurentPolynomial_mpair
            sage: LaurentPolynomial_mpair(L, {(1,2): 1/42}, mon=(-3, -3))
            1/42*w^-2*z^-1

        :issue:`22398`::

            sage: LQ = LaurentPolynomialRing(QQ, 'x0, x1, x2, y0, y1, y2, y3, y4, y5')
            sage: LZ = LaurentPolynomialRing(ZZ, 'x0, x1, x2, y0, y1, y2, y3, y4, y5')
            sage: LQ.inject_variables()
            Defining x0, x1, x2, y0, y1, y2, y3, y4, y5
            sage: x2^-1*y0*y1*y2*y3*y4*y5 + x1^-1*x2^-1*y0*y1*y3*y4 + x0^-1 in LZ
            True
            sage: x2^-1*y0*y1*y2*y3*y4*y5 + x1^-1*x2^-1*y0*y1*y3*y4 + x0^-1*x1^-1*y0*y3 + x0^-1 in LZ
            True

        Check that input is not modified::

            sage: LQ.<x,y> = LaurentPolynomialRing(QQ)
            sage: D = {(-1, 1): 1}
            sage: k = tuple(D)[0]
            sage: v = D[k]
            sage: type(k), type(v)
            (<... 'tuple'>, <class 'sage.rings.integer.Integer'>)
            sage: LQ(D)
            x^-1*y
            sage: tuple(D)[0] is k
            True
            sage: D[k] is v
            True
        """
        if isinstance(x, PolyDict):
            x = x.dict()
        if mon is not None:
            if isinstance(mon, ETuple):
                self._mon = mon
            else:
                self._mon = ETuple(mon)
        else:
            if isinstance(x, dict):
                self._mon = ETuple({}, int(parent.ngens()))
                D = {}
                for k, x_k in x.iteritems():  # ETuple-ize keys, set _mon
                    if not isinstance(k, (tuple, ETuple)) or len(k) != parent.ngens():
                        self._mon = ETuple({}, int(parent.ngens()))
                        break
                    if isinstance(k, tuple):
                        k = ETuple(k)
                    D[k] = x_k
                    self._mon = self._mon.emin(k)  # point-wise min of _mon and k
                else:
                    x = D
                if not self._mon.is_constant():  # factor out _mon
                    x = {k.esub(self._mon): x_k for k, x_k in x.iteritems()}
            elif (isinstance(x, LaurentPolynomial_mpair) and
                  parent.variable_names() == x.parent().variable_names()):
                self._mon = ( < LaurentPolynomial_mpair > x)._mon
                x = ( < LaurentPolynomial_mpair > x)._poly
            else:  # since x should coerce into parent, _mon should be (0,...,0)
                self._mon = ETuple({}, int(parent.ngens()))
        self._poly = parent._R(x)
        CommutativeAlgebraElement.__init__(self, parent)

    def __reduce__(self):
        """
        TESTS::

            sage: R = LaurentPolynomialRing(QQ, 2, 'x')
            sage: R.<x1,x2> = LaurentPolynomialRing(QQ)
            sage: loads(dumps(x1)) == x1 # indirect doctest
            True
            sage: z = x1/x2
            sage: loads(dumps(z)) == z
            True
        """
        return self._parent, (self._poly, self._mon)

    def __hash__(self):
        r"""
        TESTS:

        Test that the hash is non-constant (see also :issue:`27914`)::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: len({hash(w^i*z^j) for i in [-2..2] for j in [-2..2]})
            25

        Check that :issue:`20490` is fixed::

            sage: R.<a,b> = LaurentPolynomialRing(ZZ)
            sage: p = a*~a
            sage: p._fraction_pair()
            (a, a)
            sage: p == R.one()
            True
            sage: hash(p)
            1

        Check that :issue:`23864` is fixed (compatibility with integers, rationals
        and polynomial rings)::

            sage: L = LaurentPolynomialRing(QQ, 'x0,x1,x2')
            sage: hash(L.zero())
            0
            sage: hash(L.one())
            1
            sage: hash(-L.one())
            -2
            sage: hash(L(1/2)) == hash(1/2)
            True

            sage: R = PolynomialRing(QQ, 'x0,x1,x2')
            sage: x0,x1,x2 = R.gens()
            sage: hash(x0) == hash(L(x0))
            True
            sage: hash(1 - 7*x0 + x1*x2) == hash(L(1 - 7*x0 + x1*x2))
            True

        Check that :issue:`27914` is fixed::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: Lw = LaurentPolynomialRing(QQ, 'w')
            sage: Lz = LaurentPolynomialRing(QQ, 'z')
            sage: all(hash(w^k) == hash(Lw(w^k))
            ....:     and hash(z^k) == hash(Lz(z^k)) for k in (-5..5))
            True
            sage: p = w^-1 + 2 + w
            sage: hash(p) == hash(Lw(p))
            True
        """
        # we reimplement the hash from multipolynomial to handle negative exponents
        # (see multi_polynomial.pyx)
        cdef long result = 0
        cdef long exponent
        cdef list var_name_hash = [hash(v) for v in self._parent.variable_names()]
        cdef int p
        cdef int n = len(var_name_hash)
        cdef long c_hash
        for m, c in self._poly.iterator_exp_coeff():
            c_hash = hash(c)
            if c_hash != 0:
                for p in range(n):
                    exponent = m[p] + self._mon[p]
                    if exponent > 0:
                        c_hash = (1000003 * c_hash) ^ var_name_hash[p]
                        c_hash = (1000003 * c_hash) ^ exponent
                    elif exponent < 0:
                        c_hash = (1000003 * c_hash) ^ var_name_hash[p]
                        c_hash = (700005 * c_hash) ^ exponent
                result += c_hash

        return result

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` under the morphism defined by
        ``im_gens`` in ``codomain``.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(ZZ)
            sage: M.<u,v> = LaurentPolynomialRing(ZZ)
            sage: phi = L.hom([u,v])
            sage: phi(x^2*~y -5*y**3)            # indirect doctest
            -5*v^3 + u^2*v^-1

        TESTS:

        check compatibility with  :issue:`26105`::

            sage: # needs sage.rings.finite_rings
            sage: F.<t> = GF(4)
            sage: LF.<a,b> = LaurentPolynomialRing(F)
            sage: rho = LF.hom([b,a], base_map=F.frobenius_endomorphism())
            sage: s = t*~a + b +~t*(b**-3)*a**2; rs = rho(s); rs
            a + (t + 1)*b^-1 + t*a^-3*b^2
            sage: s == rho(rs)
            True
        """
        p = self._poly
        m = self._mon
        if base_map is not None:
            p = p.map_coefficients(base_map)
        from sage.misc.misc_c import prod
        return codomain(p(im_gens) * prod(ig**m[im_gens.index(ig)] for ig in im_gens))

    cdef _normalize(self, i=None):
        r"""
        Remove the common monomials from ``self._poly`` and store
        them in ``self._mon``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x*y + 2*y*x^2 + y  # indirect doctest
            sage: f.factor() # Notice the y has been factored out.
            (y) * (2*x^2 + x + 1)

        Check that :issue:`23864` has been fixed::

            sage: hash(L.zero())
            0
        """
        if not self._poly:
            self._mon = ETuple({}, int(self._parent.ngens()))
            return

        # cdef dict D = <dict> self._poly._mpoly_dict_recursive(
        #                                <tuple> self._parent.variable_names(),
        #                                self._parent.base_ring()
        #                                )
        cdef dict D = <dict > self._poly.dict()

        cdef ETuple e
        if i is None:
            e = None
            for k in D:
                if e is None:
                    e = <ETuple > k
                else:
                    e = e.emin(k)
            if not e.is_constant():
                self._poly = <ModuleElement > (self._poly // self._poly._parent({e: 1}))
                self._mon = self._mon.eadd(e)
        else:
            e = None
            for k in D:
                if e is None or k[i] < e:
                    e = <ETuple > k[i]
            if e > 0:
                self._poly = <ModuleElement > (self._poly // self._poly._parent.gen(i))
                self._mon = self._mon.eadd_p(e, i)

    cdef _compute_polydict(self):
        """
        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1 +3
            sage: a.dict()  # indirect doctest
            {(0, 0): 3, (2, -1): 1}
        """
        # cdef dict D = <dict> self._poly._mpoly_dict_recursive(self._parent.variable_names(),
        #                                                      self._parent.base_ring())
        cdef dict D = <dict > self._poly.dict()
        cdef dict DD
        if self._mon.is_constant():
            self._prod = PolyDict(D)
            return
        DD = {}
        for k in D:
            DD[k.eadd(self._mon)] = D[k]
        self._prod = PolyDict(DD)

    def is_unit(self):
        """
        Return ``True`` if ``self`` is a unit.

        The ground ring is assumed to be an integral domain.

        This means that the Laurent polynomial is a monomial
        with unit coefficient.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: (x*y/2).is_unit()
            True
            sage: (x + y).is_unit()
            False
            sage: (L.zero()).is_unit()
            False
            sage: (L.one()).is_unit()
            True

            sage: L.<x,y> = LaurentPolynomialRing(ZZ)
            sage: (2*x*y).is_unit()
            False
        """
        coeffs = self.coefficients()
        if len(coeffs) != 1:
            return False
        return coeffs[0].is_unit()

    def _repr_(self):
        """
        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x^2 + x*y/2 + 2*y^-1
            sage: f._repr_()
            'x^2 + 1/2*x*y + 2*y^-1'
        """
        if self._prod is None:
            self._compute_polydict()
        try:
            key = self.parent().term_order().sortkey
        except AttributeError:
            key = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self._prod.poly_repr(self.parent().variable_names(),
                                    atomic_coefficients=atomic, sortkey=key)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1+3; a
            w^2*z^-1 + 3
            sage: latex(a)
            w^{2} z^{-1} + 3

        TESTS::

            sage: L.<lambda2, y2> = LaurentPolynomialRing(QQ)
            sage: latex(1/lambda2 + y2^(-3))
            \lambda_{2}^{-1} + y_{2}^{-3}
        """
        if self._prod is None:
            self._compute_polydict()
        try:
            key = self.parent().term_order().sortkey
        except AttributeError:
            key = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self._prod.latex(self.parent().latex_variable_names(),
                                atomic_coefficients=atomic, sortkey=key)

    cpdef long number_of_terms(self) except -1:
        """
        Return the number of nonzero coefficients of ``self``.

        Also called weight, hamming weight or sparsity.

        EXAMPLES::

            sage: R.<x, y> = LaurentPolynomialRing(ZZ)
            sage: f = x^3 - y
            sage: f.number_of_terms()
            2
            sage: R(0).number_of_terms()
            0
            sage: f = (x+1/y)^100
            sage: f.number_of_terms()
            101

        The method :meth:`hamming_weight` is an alias::

            sage: f.hamming_weight()
            101
        """
        return self._poly.number_of_terms()

    def __invert__(LaurentPolynomial_mpair self):
        """
        Return the inverse of ``self``.

        This treats monomials specially so they remain Laurent
        polynomials; the inverse of any other polynomial is an element
        of the rational function field.

        TESTS::

            sage: L.<x,y> = LaurentPolynomialRing(ZZ)
            sage: f = ~x
            sage: parent(f)
            Multivariate Laurent Polynomial Ring in x, y over Integer Ring
            sage: parent(f.coefficients()[0]) is parent(f).base_ring()
            True
            sage: g = ~(2*x)
            sage: parent(g)
            Multivariate Laurent Polynomial Ring in x, y over Rational Field
            sage: parent(g.coefficients()[0]) is parent(g).base_ring()
            True
        """
        cdef ETuple e
        if self._poly.is_term():
            (e, c), = self.dict().items()
            e = e.emul(-1)
            P = self._parent
            try:
                c = c.inverse_of_unit()
            except (AttributeError, ZeroDivisionError, ArithmeticError):
                c = ~c
                if c.parent() is not P.base_ring():
                    P = P.change_ring(c.parent())
            return P({e: c})
        return super().__invert__()

    def __pow__(LaurentPolynomial_mpair self, n, mod):
        """
        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x + y
            sage: f^2
            x^2 + 2*x*y + y^2
            sage: f^(-1)
            1/(x + y)

        TESTS:

        Check that :issue:`2952` is fixed::

            sage: R.<q> = QQ[]
            sage: L.<x,y,z> = LaurentPolynomialRing(R)
            sage: f = (x+y+z^-1)^2
            sage: f.substitute(z=1)
            x^2 + 2*x*y + y^2 + 2*x + 2*y + 1

        Check that using third argument raises an error::

            sage: L.<x,y,z> = LaurentPolynomialRing(R)
            sage: pow(x + y + z, 2, x)
            Traceback (most recent call last):
            ...
            NotImplementedError: pow() with a modulus is not implemented for this ring
        """
        cdef LaurentPolynomial_mpair ans
        if mod is not None:
            raise NotImplementedError(
                "pow() with a modulus is not implemented for this ring"
            )
        if n < 0:
            return ~(self ** -n)
        ans = self._new_c()
        ans._poly = self._poly ** n
        ans._mon = self._mon.emul(n)
        return ans

    def __getitem__(self, n):
        r"""
        Return the coefficient of `x^n = x_1^{n_1} \cdots x_k^{n_k}` where
        `n` is a tuple of length `k` and `k` is the number of variables.

        If the number of inputs is not equal to the number of variables, this
        raises a :exc:`TypeError`.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3 + x*z; f
            -x^6 + x*z - 7*x^-2*y^3 + 5*x^-2*y + x^-3*y^2
            sage: f[6,0,0]
            -1
            sage: f[-2,3,0]
            -7
            sage: f[-1,4,2]
            0
            sage: f[1,0,1]
            1
            sage: f[6]
            Traceback (most recent call last):
            ...
            TypeError: must have exactly 3 inputs
            sage: f[6,0]
            Traceback (most recent call last):
            ...
            TypeError: must have exactly 3 inputs
            sage: f[6,0,0,0]
            Traceback (most recent call last):
            ...
            TypeError: must have exactly 3 inputs
        """
        if isinstance(n, slice):
            raise TypeError("multivariate Laurent polynomials are not iterable")
        if not isinstance(n, tuple) or len(n) != self._parent.ngens():
            raise TypeError("must have exactly %s inputs" %
                            self.parent().ngens())
        cdef ETuple t = ETuple(n)
        if self._prod is None:
            self._compute_polydict()
        try:
            return self._prod[t]
        except KeyError:
            return self._parent.base_ring().zero()

    def __iter__(self):
        """
        Iterate through all terms by returning a list of the coefficient and
        the corresponding monomial.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: sorted(f) # indirect doctest
            [(-7, x^-2*y^3), (-1, x^6), (1, x^-3*y^2), (5, x^-2*y)]
        """
        P = self._parent
        one = P._R.one()
        if self._mon.is_constant():
            for exp, coeff in self._poly.iterator_exp_coeff():
                yield (coeff, P.element_class(P, one, exp))
        else:
            for exp, coeff in self._poly.iterator_exp_coeff():
                yield (coeff, P.element_class(P, one, exp.eadd(self._mon)))

    def _as_extended_polynomial(self):
        """
        This Laurent polynomial seen as a polynomial in twice as many variables,
        where half of the variables are the inverses of the other half.

        EXAMPLES::

            sage: L.<t1,t2> = LaurentPolynomialRing(QQ)
            sage: f = t1-t2^-2
            sage: f
            t1 - t2^-2
            sage: f._as_extended_polynomial()
            -t2inv^2 + t1
            sage: _.parent()
            Multivariate Polynomial Ring in t1, t1inv, t2, t2inv over Rational Field
        """
        dres = {}
        for (e, c) in self.dict().items():
            exps = []
            for t in e:
                if t > 0:
                    exps += [t, 0]
                else:
                    exps += [0, -t]
            dres[tuple(exps)] = c
        return self.parent()._extended_ring(dres)

    def iterator_exp_coeff(self):
        """
        Iterate over ``self`` as pairs of (ETuple, coefficient).

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: list(f.iterator_exp_coeff())
            [((6, 0), -1), ((-2, 3), -7), ((-2, 1), 5), ((-3, 2), 1)]
        """
        for exp, coeff in self._poly.iterator_exp_coeff():
            yield (exp.eadd(self._mon), coeff)

    def monomials(self):
        """
        Return the list of monomials in ``self``.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: sorted(f.monomials())
            [x^-3*y^2, x^-2*y, x^-2*y^3, x^6]
        """
        return [mon for coeff, mon in self]

    def monomial_coefficient(self, mon):
        """
        Return the coefficient in the base ring of the monomial ``mon`` in
        ``self``, where ``mon`` must have the same parent as ``self``.

        This function contrasts with the function :meth:`coefficient()`
        which returns the coefficient of a monomial viewing this
        polynomial in a polynomial ring over a base ring having fewer
        variables.

        INPUT:

        - ``mon`` -- a monomial

        .. SEEALSO::

            For coefficients in a base ring of fewer variables, see
            :meth:`coefficient()`.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: f.monomial_coefficient(x^-2*y^3)
            -7
            sage: f.monomial_coefficient(x^2)
            0

        TESTS::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = y^2 * x^-2
            sage: f.monomial_coefficient(x + y)
            Traceback (most recent call last):
            ...
            ValueError: not a monomial
        """
        if parent(mon) != self._parent:
            raise TypeError("input must have the same parent")
        cdef LaurentPolynomial_mpair m = <LaurentPolynomial_mpair > mon
        if m._prod is None:
            m._compute_polydict()
        if self._prod is None:
            self._compute_polydict()
        exp = monomial_exponent(m._prod)
        zero = self._parent.base_ring().zero()
        return self._prod.get(exp, zero)

    def constant_coefficient(self):
        """
        Return the constant coefficient of ``self``.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^2 + 5*x*y)*x^-3; f
            -x^6 - 7*x^-2*y^2 + 5*x^-2*y + x^-3*y^2
            sage: f.constant_coefficient()
            0
            sage: f = (x^3 + 2*x^-2*y+y^3)*y^-3; f
            x^3*y^-3 + 1 + 2*x^-2*y^-2
            sage: f.constant_coefficient()
            1
        """
        return self[(0,)*self._parent.ngens()]

    def coefficient(self, mon):
        r"""
        Return the coefficient of ``mon`` in ``self``, where ``mon`` must
        have the same parent as ``self``.

        The coefficient is defined as follows. If `f` is this polynomial, then
        the coefficient `c_m` is sum:

        .. MATH::

            c_m := \sum_T \frac{T}{m}

        where the sum is over terms `T` in `f` that are exactly divisible
        by `m`.

        A monomial `m(x,y)` 'exactly divides' `f(x,y)` if `m(x,y) | f(x,y)`
        and neither `x \cdot m(x,y)` nor `y \cdot m(x,y)` divides `f(x,y)`.

        INPUT:

        - ``mon`` -- a monomial

        OUTPUT: element of the parent of ``self``

        .. NOTE::

            To get the constant coefficient, call
            :meth:`constant_coefficient()`.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)

        The coefficient returned is an element of the parent of ``self``; in
        this case, ``P``. ::

            sage: f = 2 * x * y
            sage: c = f.coefficient(x*y); c
            2
            sage: c.parent()
            Multivariate Laurent Polynomial Ring in x, y over Rational Field

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^2 + 5*x*y)*x^-3; f
            -x^6 - 7*x^-2*y^2 + 5*x^-2*y + x^-3*y^2
            sage: f.coefficient(y)
            5*x^-2
            sage: f.coefficient(y^2)
            -7*x^-2 + x^-3
            sage: f.coefficient(x*y)
            0
            sage: f.coefficient(x^-2)
            -7*y^2 + 5*y
            sage: f.coefficient(x^-2*y^2)
            -7
            sage: f.coefficient(1)
            -x^6 - 7*x^-2*y^2 + 5*x^-2*y + x^-3*y^2
        """
        if mon.parent() is not self._parent:
            mon = self._parent(mon)
        cdef LaurentPolynomial_mpair m = <LaurentPolynomial_mpair > mon
        if self._prod is None:
            self._compute_polydict()
        if m._prod is None:
            m._compute_polydict()
        return self._parent(self._prod.coefficient(m.dict()))

    def coefficients(self):
        """
        Return the nonzero coefficients of ``self`` in a list.

        The returned list is decreasingly ordered by the term ordering
        of ``self.parent()``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ, order='degrevlex')
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.coefficients()
            [4, 3, 2, 1]
            sage: L.<x,y,z> = LaurentPolynomialRing(QQ,order='lex')
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.coefficients()
            [4, 1, 2, 3]
        """
        return self._poly.coefficients()

    def variables(self, sort=True):
        """
        Return a tuple of all variables occurring in ``self``.

        INPUT:

        - ``sort`` -- specifies whether the indices shall be sorted

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.variables()
            (z, y, x)
            sage: f.variables(sort=False) #random
            (y, z, x)
        """
        cdef dict d = self.dict()
        cdef tuple g = self._parent.gens()
        cdef Py_ssize_t nvars = len(g)
        cdef set vars = set()
        for k in d:
            vars.update(k.nonzero_positions())
            if len(vars) == nvars:
                break
        cdef list v = [g[i] for i in vars]
        if sort:
            v.sort()
        return tuple(v)

    cpdef dict dict(self):
        """
        Return ``self`` represented as a ``dict``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: sorted(f.dict().items())
            [((3, 1, 0), 3), ((4, 0, -2), 2), ((6, -7, 0), 1), ((7, 0, -1), 4)]
        """
        if self._prod is None:
            self._compute_polydict()
        return < dict > self._prod.dict()

    def _fraction_pair(self):
        """
        Return one representation of ``self`` as a pair
        ``(numerator, denominator)``.

        Here both the numerator and the denominator are polynomials.

        This is used for coercion into the fraction field.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f._fraction_pair()
            (4*x^7*y^7*z + 3*x^3*y^8*z^2 + 2*x^4*y^7 + x^6*z^2, y^7*z^2)
        """
        ring = self._parent._R
        numer = self._poly
        denom = ring.one()
        var = ring.gens()
        for i, j in enumerate(self._mon):
            if j > 0:
                numer *= var[i] ** j
            else:
                denom *= var[i] ** (-j)
        return (numer, denom)

    cpdef _add_(self, _right):
        """
        Return the Laurent polynomial ``self + right``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: g = y + z
            sage: f + g
            x + y + z + y^-1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        cdef LaurentPolynomial_mpair right = <LaurentPolynomial_mpair > _right
        ans._mon, a, b = self._mon.combine_to_positives(right._mon)
        if not a.is_constant():
            ans._poly = self._poly * self._poly._parent({a: 1})
        else:
            ans._poly = self._poly
        if not b.is_constant():
            ans._poly += right._poly * self._poly._parent({b: 1})
        else:
            ans._poly += right._poly
        return ans

    cpdef _sub_(self, _right):
        """
        Return the Laurent polynomial ``self - right``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: g = y + z + x
            sage: f - g
            -y - z + y^-1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        cdef LaurentPolynomial_mpair right = <LaurentPolynomial_mpair > _right
        cdef ETuple a, b
        ans._mon, a, b = self._mon.combine_to_positives(right._mon)
        if not a.is_constant():
            ans._poly = self._poly * self._poly._parent({a: 1})
        else:
            ans._poly = self._poly
        if not b.is_constant():
            ans._poly -= right._poly * self._poly._parent({b: 1})
        else:
            ans._poly -= right._poly
        return ans

    cpdef _div_(self, rhs):
        """
        Return the division of ``self`` by ``rhs``.

        If the denominator is not a unit,
        the result will be given in the fraction field.

        EXAMPLES::

            sage: R.<s,q,t> = LaurentPolynomialRing(QQ)
            sage: 1/s
            s^-1
            sage: 1/(s*q)
            s^-1*q^-1
            sage: 1/(s+q)
            1/(s + q)
            sage: (1/(s+q)).parent()
            Fraction Field of Multivariate Polynomial Ring in s, q, t over Rational Field
            sage: (1/(s*q)).parent()
            Multivariate Laurent Polynomial Ring in s, q, t over Rational Field
            sage: (s+q)/(q^2*t^(-2))
            s*q^-2*t^2 + q^-1*t^2
        """
        cdef LaurentPolynomial_mpair right = <LaurentPolynomial_mpair > rhs
        if right.is_zero():
            raise ZeroDivisionError
        if right._poly.is_term():
            return self * ~right
        else:
            return RingElement._div_(self, rhs)

    def is_monomial(self):
        """
        Return ``True`` if ``self`` is a monomial.

        EXAMPLES::

            sage: k.<y,z> = LaurentPolynomialRing(QQ)
            sage: z.is_monomial()
            True
            sage: k(1).is_monomial()
            True
            sage: (z+1).is_monomial()
            False
            sage: (z^-2909).is_monomial()
            True
            sage: (38*z^-2909).is_monomial()
            False
        """
        return self._poly.is_monomial()

    cpdef _neg_(self):
        """
        Return ``-self``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: -f
            -x - y^-1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        ans._mon = self._mon
        ans._poly = -self._poly
        return ans

    cpdef _lmul_(self, Element right):
        """
        Return ``self * right`` where ``right`` is in ``self``'s base ring.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: f*(1/2)
            1/2*x + 1/2*y^-1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        ans._mon = self._mon
        ans._poly = self._poly * right
        return ans

    cpdef _rmul_(self, Element left):
        """
        Return ``left * self`` where ``left`` is in ``self``'s base ring.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: (1/2)*f
            1/2*x + 1/2*y^-1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        ans._mon = self._mon
        ans._poly = left * self._poly
        return ans

    cpdef _mul_(self, right):
        """
        Return ``self * right``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: g = y + z
            sage: f*g
            x*y + x*z + 1 + y^-1*z
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        ans._mon = self._mon.eadd(( < LaurentPolynomial_mpair > right)._mon)
        ans._poly = self._poly * ( < LaurentPolynomial_mpair > right)._poly
        return ans

    cpdef _floordiv_(self, right):
        """
        Perform division with remainder and return the quotient.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x^3 + y^-3
            sage: g = y + x
            sage: f // g                                                                # needs sage.libs.singular
            x^5*y^-3 - x^4*y^-2 + x^3*y^-1

            sage: h = x + y^(-1)
            sage: f // h                                                                # needs sage.libs.singular
            x^2 - x*y^-1 + y^-2
            sage: h * (f // h) == f                                                     # needs sage.libs.singular
            True
            sage: f // 1
            x^3 + y^-3
            sage: 1 // f                                                                # needs sage.libs.singular
            0

        TESTS:

        Check that :issue:`19357` is fixed::

            sage: x // y
            x*y^-1

        Check that :issue:`21999` is fixed::

            sage: L.<a,b> = LaurentPolynomialRing(QQbar)                                # needs sage.rings.number_field
            sage: (a+a*b) // a                                                          # needs sage.libs.singular sage.rings.number_field
            b + 1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        cdef LaurentPolynomial_mpair rightl = <LaurentPolynomial_mpair > right
        self._normalize()
        rightl._normalize()
        ans._mon = self._mon.esub(rightl._mon)
        ans._poly = self._poly // rightl._poly
        return ans

    @coerce_binop
    def quo_rem(self, right):
        """
        Divide this Laurent polynomial by ``right`` and return a quotient and
        a remainder.

        INPUT:

        - ``right`` -- a Laurent polynomial

        OUTPUT: a pair of Laurent polynomials

        EXAMPLES::

            sage: R.<s, t> = LaurentPolynomialRing(QQ)
            sage: (s^2 - t^2).quo_rem(s - t)                                            # needs sage.libs.singular
            (s + t, 0)
            sage: (s^-2 - t^2).quo_rem(s - t)                                           # needs sage.libs.singular
            (s + t, -s^2 + s^-2)
            sage: (s^-2 - t^2).quo_rem(s^-1 - t)                                        # needs sage.libs.singular
            (t + s^-1, 0)

        TESTS:

        Verify that :issue:`31257` is fixed::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = LaurentPolynomialRing(QQ)
            sage: q, r = (1/x).quo_rem(y)
            sage: q, r
            (x^-1*y^-1, 0)
            sage: q*y + r == 1/x
            True
            sage: q, r = (x^-2 - y^2).quo_rem(x - y)
            sage: q*(x - y) + r == x^-2 - y^2
            True
        """
        # make copies of self and right so that the input can be normalized
        # without affecting the objects that were passed to the method
        cdef LaurentPolynomial_mpair selfl = self._new_c()
        selfl._poly = self._poly
        selfl._mon = self._mon
        cdef LaurentPolynomial_mpair rightl = self._new_c()
        rightl._poly = (< LaurentPolynomial_mpair > right)._poly
        rightl._mon = (< LaurentPolynomial_mpair > right)._mon

        selfl._normalize()
        rightl._normalize()
        q, r = selfl._poly.quo_rem(rightl._poly)
        ql = LaurentPolynomial_mpair(self._parent, q,
                                     mon=selfl._mon.esub(rightl._mon))
        rl = LaurentPolynomial_mpair(self._parent, r,
                                     mon=selfl._mon)
        ql._normalize()
        rl._normalize()
        return (ql, rl)

    cpdef _richcmp_(self, right, int op):
        """
        Compare two polynomials in a `LaurentPolynomialRing` based on the term
        order from the parent ring.  If the parent ring does not specify a term
        order then only comparison by equality is supported.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: g = y + z
            sage: f == f
            True
            sage: f == g
            False
            sage: f == 2
            False
        """
        if self._prod is None:
            self._compute_polydict()
        if (< LaurentPolynomial_mpair > right)._prod is None:
            (< LaurentPolynomial_mpair > right)._compute_polydict()

        try:
            sortkey = self._parent.term_order().sortkey
        except AttributeError:
            sortkey = None

        return self._prod.rich_compare(( < LaurentPolynomial_mpair > right)._prod,
                                       op, sortkey)

    def exponents(self):
        """
        Return a list of the exponents of ``self``.

        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1 + 3; a
            w^2*z^-1 + 3
            sage: e = a.exponents()
            sage: e.sort(); e
            [(0, 0), (2, -1)]
        """
        return [a.eadd(self._mon) for a in self._poly.exponents()]

    def degree(self, x=None):
        r"""
        Return the degree of ``self``.

        INPUT:

        - ``x`` -- (default: ``None``) a generator of the parent ring

        OUTPUT:

        If ``x`` is ``None``, return the total degree of ``self``.
        If ``x`` is a given generator of the parent ring,
        the output is the maximum degree of ``x`` in ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.degree()
            6
            sage: f.degree(x)
            7
            sage: f.degree(y)
            1
            sage: f.degree(z)
            0

        The zero polynomial is defined to have degree `-\infty`::

            sage: R.<x, y, z> = LaurentPolynomialRing(ZZ)
            sage: R.zero().degree()
            -Infinity
            sage: R.zero().degree(x)
            -Infinity
            sage: R.zero().degree(x) == R.zero().degree(y) == R.zero().degree(z)
            True

        TESTS::

            sage: R.<x, y, z> = LaurentPolynomialRing(ZZ)
            sage: f = x + y + z
            sage: f.degree(1)
            Traceback (most recent call last):
            ...
            TypeError: 1 is not a generator of parent
        """
        # The zero polynomial is defined to have degree -Infinity
        if self.is_zero():
            return minus_infinity

        if x is None:
            return self._poly.total_degree() + sum(self._mon)

        # Get the index of the generator or error
        cdef tuple g = <tuple > self._parent.gens()
        cdef Py_ssize_t i
        try:
            i = g.index(x)
        except ValueError:  # not in the tuple
            raise TypeError(f"{x} is not a generator of parent")
        return self._poly.degree(self._parent._R.gens()[i]) + self._mon[i]

    def valuation(self, x=None):
        r"""
        Return the valuation of ``self``.

        If ``x`` is ``None``, the returned valuation is the minimal total degree
        of the monomials occurring in ``self``. Geometrically, this is the order
        of vanishing of ``self`` at the generic point of the blow-up of the
        point `(0,0,\ldots,0)`.

        If ``x`` is not ``None``, then it must be a generator. In that case, the
        minimum degree of that generator occurring in ``self`` is returned.
        Geometrically, this is the order of vanishing of ``self`` at the generic
        point of the curve `x = 0`.

        INPUT:

        - ``x`` -- (optional) a generator; if given, return the valuation
          with respect to this generator

        EXAMPLES::

            sage: R.<x,y> = LaurentPolynomialRing(ZZ)
            sage: f = 2*x^2*y^-3 - 13*x^-1*y^-3 + 2*x^2*y^-5 - 2*x^-3*y^2
            sage: f.valuation()
            -4
            sage: f.valuation(x)
            -3
            sage: f.valuation(y)
            -5
            sage: R.zero().valuation()
            +Infinity

        TESTS:

        If supplied, ``x`` must be a generator::

            sage: R.<x,y> = LaurentPolynomialRing(ZZ)
            sage: f = 1 + x + x^2*y^-1
            sage: f.valuation(1)
            Traceback (most recent call last):
            ...
            TypeError: 1 is not a generator of parent
        """
        # Valuation of zero polynomial is defined to be +Infinity
        if self.is_zero():
            return Infinity

        # When x is None find the minimal valuation by finding the minimal
        # valuation of the sum of exponents
        if x is None:
            return Integer(min(sum(e) for e in self.exponents()))

        # Get the index of the generator or error
        cdef tuple g = <tuple > self._parent.gens()
        cdef Py_ssize_t i
        try:
            i = g.index(x)
        except ValueError:  # not in the tuple
            raise TypeError(f"{x} is not a generator of parent")

        # Find the minimal valuation of x by checking each term
        return Integer(min(e[i] for e in self.exponents()))

    def has_inverse_of(self, i):
        """
        INPUT:

        - ``i`` -- the index of a generator of ``self.parent()``

        OUTPUT:

        Return ``True`` if ``self`` contains a monomial including the inverse of
        ``self.parent().gen(i)``, ``False`` otherwise.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.has_inverse_of(0)
            False
            sage: f.has_inverse_of(1)
            True
            sage: f.has_inverse_of(2)
            True
        """
        if (not isinstance(i, (int, Integer))) or (i < 0) or (i >= self._parent.ngens()):
            raise TypeError("argument is not the index of a generator")
        if self._mon[i] < 0:
            self._normalize(i)
            if self._mon[i] < 0:
                return True
            return False
        return False

    def has_any_inverse(self):
        """
        Return ``True`` if ``self`` contains any monomials with a negative
        exponent, ``False`` otherwise.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.has_any_inverse()
            True
            sage: g = x^2 + y^2
            sage: g.has_any_inverse()
            False
        """
        for m in self._mon.nonzero_values(sort=False):
            if m < 0:
                return True
        return False

    def __call__(self, *x, **kwds):
        """
        Compute value of ``self`` at ``x``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + 2*y + 3*z
            sage: f(1,1,1)
            6
            sage: f = x^-1 + y + z
            sage: f(0,1,1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        TESTS::

            sage: f = x + 2*y + 3*z
            sage: f(2)
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match the number of generators in parent
            sage: f(2,0)
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match the number of generators in parent
            sage: f( (1,1,1) )
            6
        """
        if kwds:
            f = self.subs(**kwds)
            if x:  # More than 1 non-keyword argument
                return f(*x)
            else:
                return f

        cdef int l = len(x)

        if l == 1 and isinstance(x[0], (tuple, list)):
            x = x[0]
            l = len(x)

        if l != self._parent.ngens():
            raise TypeError("number of arguments does not match the number"
                            " of generators in parent")

        # Check to make sure that we aren't dividing by zero
        cdef Py_ssize_t m
        for m in range(l):
            if x[m] == 0:
                if self.has_inverse_of(m):
                    raise ZeroDivisionError

        ans = self._poly(*x)
        if ans:
            for m in self._mon.nonzero_positions():
                ans *= x[m]**self._mon[m]

        return ans

    def subs(self, in_dict=None, **kwds):
        """
        Substitute some variables in this Laurent polynomial.

        Variable/value pairs for the substitution may be given
        as a dictionary or via keyword-value pairs. If both are
        present, the latter take precedence.

        INPUT:

        - ``in_dict`` -- dictionary (optional)

        - ``**kwds`` -- keyword arguments

        OUTPUT: a Laurent polynomial

        EXAMPLES::

            sage: L.<x, y, z> = LaurentPolynomialRing(QQ)
            sage: f = x + 2*y + 3*z
            sage: f.subs(x=1)
            2*y + 3*z + 1
            sage: f.subs(y=1)
            x + 3*z + 2
            sage: f.subs(z=1)
            x + 2*y + 3
            sage: f.subs(x=1, y=1, z=1)
            6

            sage: f = x^-1
            sage: f.subs(x=2)
            1/2
            sage: f.subs({x: 2})
            1/2

            sage: f = x + 2*y + 3*z
            sage: f.subs({x: 1, y: 1, z: 1})
            6
            sage: f.substitute(x=1, y=1, z=1)
            6

        TESTS::

            sage: f = x + 2*y + 3*z
            sage: f(q=10)
            x + 2*y + 3*z

            sage: x.subs({x: 2}, x=1)
            1

            sage: f.subs({1: 2}, x=1)
            3*z + 5
        """
        cdef list variables = list(self._parent.gens())
        cdef Py_ssize_t i
        for i in range(len(variables)):
            if str(variables[i]) in kwds:
                variables[i] = kwds[str(variables[i])]
            elif in_dict:
                if variables[i] in in_dict:
                    variables[i] = in_dict[variables[i]]
                elif i in in_dict:
                    variables[i] = in_dict[i]
        return self(tuple(variables))

    def is_constant(self):
        r"""
        Return whether this Laurent polynomial is constant.

        EXAMPLES::

            sage: L.<a, b> = LaurentPolynomialRing(QQ)
            sage: L(0).is_constant()
            True
            sage: L(42).is_constant()
            True
            sage: a.is_constant()
            False
            sage: (1/b).is_constant()
            False
        """
        return (self._mon == ETuple({}, int(self._parent.ngens())) and
                self._poly.is_constant())

    def _symbolic_(self, R):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: R.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x^3 + y/x
            sage: g = f._symbolic_(SR); g
            (x^4 + y)/x
            sage: g(x=2, y=2)
            9
            sage: g = SR(f)
            sage: g(x=2, y=2)
            9
        """
        d = {repr(g): R.var(g) for g in self._parent.gens()}
        return self.subs(**d)

    def derivative(self, *args):
        r"""
        The formal derivative of this Laurent polynomial, with respect
        to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global :func:`derivative` function for more
        details.

        .. SEEALSO::

            :meth:`_derivative`

        EXAMPLES::

            sage: R = LaurentPolynomialRing(ZZ,'x, y')
            sage: x, y = R.gens()
            sage: t = x**4*y + x*y + y + x**(-1) + y**(-3)
            sage: t.derivative(x, x)
            12*x^2*y + 2*x^-3
            sage: t.derivative(y, 2)
            12*y^-5
        """
        return multi_derivative(self, args)

    # add .diff(), .differentiate() as aliases for .derivative()
    diff = differentiate = derivative

    def _derivative(self, var=None):
        """
        Compute formal derivative of this Laurent polynomial with
        respect to the given variable.

        If ``var`` is among the generators of this ring, the derivative
        is with respect to the generator. Otherwise, ``_derivative(var)`` is called
        recursively for each coefficient of this polynomial.

        .. SEEALSO:: :meth:`derivative`

        EXAMPLES::

            sage: R = LaurentPolynomialRing(ZZ,'x, y')
            sage: x, y = R.gens()
            sage: t = x**4*y+x*y+y+x**(-1)+y**(-3)
            sage: t._derivative(x)
            4*x^3*y + y - x^-2
            sage: t._derivative(y)
            x^4 + x + 1 - 3*y^-4

            sage: R = LaurentPolynomialRing(QQ['z'],'x')
            sage: z = R.base_ring().gen()
            sage: x = R.gen()
            sage: t = 33*z*x**4+x**(-1)
            sage: t._derivative(z)
            33*x^4
            sage: t._derivative(x)
            -x^-2 + 132*z*x^3
        """
        if var is None:
            raise ValueError("must specify which variable to differentiate "
                             "with respect to")
        P = self._parent
        cdef list gens = list(P.gens())

        # check if var is one of the generators
        try:
            index = gens.index(var)
        except ValueError:
            # call _derivative() recursively on coefficients
            return P({m: c._derivative(var)
                      for (m, c) in self.dict().iteritems()})

        # compute formal derivative with respect to generator
        cdef dict d = {}
        for m, c in self.dict().iteritems():
            if m[index] != 0:
                new_m = [u for u in m]
                new_m[index] += -1
                d[ETuple(new_m)] = m[index] * c
        return P(d)

    def is_univariate(self):
        """
        Return ``True`` if this is a univariate or constant Laurent polynomial,
        and ``False`` otherwise.

        EXAMPLES::

            sage: R.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = (x^3 + y^-3)*z
            sage: f.is_univariate()
            False
            sage: g = f(1, y, 4)
            sage: g.is_univariate()
            True
            sage: R(1).is_univariate()
            True
        """
        return len(self.variables()) < 2

    def univariate_polynomial(self, R=None):
        """
        Return a univariate polynomial associated to this
        multivariate polynomial.

        INPUT:

        - ``R`` -- (default: ``None``) a univariate Laurent polynomial ring

        If this polynomial is not in at most one variable, then a
        :exc:`ValueError` exception is raised.  The new polynomial is over
        the same base ring as the given :class:`LaurentPolynomial` and in the
        variable ``x`` if no ring ``R`` is provided.

        EXAMPLES::

            sage: R.<x, y> = LaurentPolynomialRing(ZZ)
            sage: f = 3*x^2 - 2*y^-1 + 7*x^2*y^2 + 5
            sage: f.univariate_polynomial()
            Traceback (most recent call last):
            ...
            TypeError: polynomial must involve at most one variable
            sage: g = f(10, y); g
            700*y^2 + 305 - 2*y^-1
            sage: h = g.univariate_polynomial(); h
            -2*y^-1 + 305 + 700*y^2
            sage: h.parent()
            Univariate Laurent Polynomial Ring in y over Integer Ring
            sage: g.univariate_polynomial(LaurentPolynomialRing(QQ,'z'))
            -2*z^-1 + 305 + 700*z^2

        Here's an example with a constant multivariate polynomial::

            sage: g = R(1)
            sage: h = g.univariate_polynomial(); h
            1
            sage: h.parent()
            Univariate Laurent Polynomial Ring in x over Integer Ring
        """
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        v = self.variables()
        if len(v) > 1:
            raise TypeError("polynomial must involve at most one variable")
        elif len(v) == 1:
            x = v[0]
            i = self._parent.gens().index(x)
            x = str(x)
        else:
            x = 'x'
            i = 0

        # construct ring if none
        if R is None:
            R = LaurentPolynomialRing(self.base_ring(), x)

        return R({m[i]: c for m, c in self.dict().iteritems()})

    def monomial_reduction(self):
        """
        Factor ``self`` into a polynomial and a monomial.

        OUTPUT:

        A tuple ``(p, v)`` where ``p`` is the underlying polynomial and ``v``
        is a monomial.

        EXAMPLES::

            sage: R.<x, y> = LaurentPolynomialRing(QQ)
            sage: f = y / x + x^2 / y + 3 * x^4 * y^-2
            sage: f.monomial_reduction()
            (3*x^5 + x^3*y + y^3, x^-1*y^-2)
            sage: f = y * x + x^2 / y + 3 * x^4 * y^-2
            sage: f.monomial_reduction()
             (3*x^3 + y^3 + x*y, x*y^-2)
            sage: x.monomial_reduction()
            (1, x)
            sage: (y^-1).monomial_reduction()
            (1, y^-1)
        """
        self._normalize()
        ring = self._parent
        g = ring.gens()
        mon = ring.prod(g[i] ** j for i, j in enumerate(self._mon))
        return (self._poly, mon)

    def factor(self):
        """
        Return a Laurent monomial (the unit part of the factorization) and a factored multi-polynomial.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.factor()
            (x^3*y^-7*z^-2) * (4*x^4*y^7*z + 3*y^8*z^2 + 2*x*y^7 + x^3*z^2)

        TESTS:

        Tests for :issue:`29173`::

            sage: L.<a, b> = LaurentPolynomialRing(ZZ, 'a, b')
            sage: (a*b + a + b + 1).factor()
            (b + 1) * (a + 1)
            sage: ((a^-1)*(a*b + a + b + 1)).factor()
            (a^-1) * (b + 1) * (a + 1)
            sage: L(-12).factor()
            -1 * 2^2 * 3
        """
        pf = self._poly.factor()

        if self._poly.degree() == 0:
            # Factorization is broken for polynomials, see
            # https://github.com/sagemath/sage/issues/20214
            return pf

        u = self.parent(pf.unit())

        cdef tuple g = <tuple > self._parent.gens()
        for i in self._mon.nonzero_positions():
            u *= g[i] ** self._mon[i]

        cdef list f = []
        cdef dict d
        for t in pf:
            d = <dict > (t[0].dict())
            if len(d) == 1:  # monomials are units
                u *= self.parent(d) ** t[1]
            else:
                f.append((self.parent(d), t[1]))

        return Factorization(f, unit=u)

    def is_square(self, root=False):
        r"""
        Test whether this Laurent polynomial is a square.

        INPUT:

        - ``root`` -- boolean (default: ``False``); if set to ``True``
          then return a pair ``(True, sqrt)`` with ``sqrt`` a square
          root of this Laurent polynomial when it exists or
          ``(False, None)``

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: p = 1 + x*y + z^-3
            sage: (p**2).is_square()
            True
            sage: (p**2).is_square(root=True)
            (True, x*y + 1 + z^-3)

            sage: x.is_square()
            False
            sage: x.is_square(root=True)
            (False, None)

            sage: (x**-4 * (1 + z)).is_square(root=False)
            False
            sage: (x**-4 * (1 + z)).is_square(root=True)
            (False, None)
        """
        self._normalize()
        if not self._mon.is_multiple_of(2):
            return (False, None) if root else False

        cdef LaurentPolynomial_mpair ans

        if not root:
            return self._poly.is_square(root=False)
        else:
            (pans, root) = self._poly.is_square(root=True)
            if not pans:
                return (False, None)

            mon = self._mon.escalar_div(2)
            ans = self._new_c()
            ans._mon = mon
            ans._poly = root
            return (True, ans)

    cpdef rescale_vars(self, dict d, h=None, new_ring=None):
        r"""
        Rescale variables in a Laurent polynomial.

        INPUT:

        - ``d`` -- a ``dict`` whose keys are the generator indices
          and values are the coefficients; so a pair ``(i, v)``
          means `x_i \mapsto v x_i`
        - ``h`` -- (optional) a map to be applied to coefficients
          done after rescaling
        - ``new_ring`` -- (optional) a new ring to map the result into

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: p = x^-2*y + x*y^-2
            sage: p.rescale_vars({0: 2, 1: 3})
            2/9*x*y^-2 + 3/4*x^-2*y
            sage: F = GF(2)
            sage: p.rescale_vars({0: 3, 1: 7}, new_ring=L.change_ring(F))
            x*y^-2 + x^-2*y

        Test for :issue:`30331`::

            sage: F.<z> = CyclotomicField(3)                                            # needs sage.rings.number_field
            sage: p.rescale_vars({0: 2, 1: z}, new_ring=L.change_ring(F))               # needs sage.rings.number_field
            2*z*x*y^-2 + 1/4*z*x^-2*y
        """
        cdef int i
        cdef dict df
        cdef ETuple v
        cdef LaurentPolynomial_mpair ans

        if self._prod is None:
            self._compute_polydict()

        df = dict(self._prod.__repn)  # This makes a copy for us to manipulate
        if new_ring is None:
            R = self._parent._base
        else:
            R = new_ring._base
        if h is None:
            for v in df:
                val = df[v]
                for i in d:
                    val *= d[i]**v[i]
                df[v] = val
        else:
            for v in df:
                val = df[v]
                for i in d:
                    val *= d[i]**v[i]
                df[v] = R(h(val))

        ans = <LaurentPolynomial_mpair > self._new_c()
        ans._prod = PolyDict(df)
        ans._mon = self._mon
        if new_ring is None:
            S = self._poly._parent
        else:
            S = self._poly._parent.change_ring(R)
        ans._poly = <MPolynomial > S({v.esub(ans._mon): df[v] for v in df})
        if new_ring is not None:
            return new_ring(ans)
        return ans

    cpdef toric_coordinate_change(self, M, h=None, new_ring=None):
        r"""
        Apply a matrix to the exponents in a Laurent polynomial.

        For efficiency, we implement this directly, rather than as a substitution.

        The optional argument ``h`` is a map to be applied to coefficients.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: p = 2*x^2 + y - x*y
            sage: p.toric_coordinate_change(Matrix([[1,-3], [1,1]]))
            2*x^2*y^2 - x^-2*y^2 + x^-3*y
            sage: F = GF(2)
            sage: p.toric_coordinate_change(Matrix([[1,-3], [1,1]]),
            ....:                           new_ring=L.change_ring(F))
            x^-2*y^2 + x^-3*y
        """
        cdef int n, i, j, x
        cdef dict d, dr
        cdef ETuple v
        cdef LaurentPolynomial_mpair ans
        cdef list mon, exp
        cdef Matrix mat = M

        n = self._parent.ngens()
        if mat.dimensions() != (n, n):
            raise ValueError("the matrix M must be a {k} x {k} matrix".format(k=n))

        if not self:
            if new_ring is None:
                return self._parent.zero()
            else:
                return new_ring.zero()

        if self._prod is None:
            self._compute_polydict()

        d = self._prod.__repn
        dr = {}
        mon = [0] * n
        for v in d:
            # Make a copy of mon as this might be faster than creating the data from scratch.
            # We will set every entry, so no need to clear the data.
            exp = list(mon)
            for j in range(n):
                x = 0
                for i in range(n):
                    if not mat.get_is_zero_unsafe(j, i):
                        x += (< int > v[i]) * int(mat.get_unsafe(j, i))
                if x < (< int > mon[j]):
                    mon[j] = x
                exp[j] = x
            dr[ETuple(exp)] = d[v]

        if h is not None:
            for v in dr:
                dr[v] = self._parent._base(h(dr[v]))

        ans = <LaurentPolynomial_mpair > self._new_c()
        ans._prod = PolyDict(dr)
        ans._mon = ETuple(mon)
        ans._poly = <MPolynomial > self._poly._parent({v.esub(ans._mon): dr[v] for v in dr})
        if new_ring is not None:
            return new_ring(ans)
        return ans

    cpdef toric_substitute(self, v, v1, a, h=None, new_ring=None):
        r"""
        Perform a single-variable substitution up to a toric coordinate change.

        The optional argument ``h`` is a map to be applied to coefficients.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: p = x + y
            sage: p.toric_substitute((2,3), (-1,1), 2)
            1/2*x^3*y^3 + 2*x^-2*y^-2
            sage: F = GF(5)
            sage: p.toric_substitute((2,3), (-1,1), 2, new_ring=L.change_ring(F))
            3*x^3*y^3 + 2*x^-2*y^-2

        TESTS:

        Tests for :issue:`30331`::

            sage: L.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: p = x + y
            sage: F.<z> = CyclotomicField(3)                                            # needs sage.rings.number_field
            sage: p.toric_substitute((2,3), (-1,1), z, new_ring=L.change_ring(F))       # needs sage.rings.number_field
            (-z - 1)*x^3*y^3 + z*x^-2*y^-2

            sage: P.<x> = LaurentPolynomialRing(QQ, 1)
            sage: u = x - 1
            sage: v = u.toric_substitute((-1,), (-1,), 1)
            sage: v.is_zero()
            True
        """
        cdef dict d, dr
        cdef ETuple ve, v1e, w, w1, mon
        cdef LaurentPolynomial_mpair ans
        cdef int t

        if self._prod is None:
            self._compute_polydict()

        d = self._prod.__repn
        dr = {}
        ve = ETuple(v)
        v1e = ETuple(v1)
        mon = self._mon
        if h is not None:
            d = dict(d)  # Make a copy so we can manipulate it
            for w in d:
                d[w] = h(d[w])
        for w in d:
            x = d[w]
            t = w.dotprod(v1e)
            w1 = w.eadd_scaled(ve, -t)
            if w1 in dr:
                dr[w1] += x * a**t
            else:
                dr[w1] = x * a**t
            mon = mon.emin(w1)
        for v in tuple(dr.keys()):
            if not dr[v]:
                del dr[v]

        if new_ring is None:
            S = self._poly._parent
        else:
            S = self._poly._parent.change_ring(new_ring._base)
        ans = <LaurentPolynomial_mpair > self._new_c()
        ans._prod = PolyDict(dr)
        ans._mon = mon
        ans._poly = <MPolynomial > S({v.esub(ans._mon): dr[v] for v in dr})
        if new_ring is not None:
            return new_ring(ans)
        return ans

    @coerce_binop
    def divides(self, other):
        """
        Check if ``self`` divides ``other``.

        EXAMPLES::

            sage: R.<x,y> = LaurentPolynomialRing(QQ)
            sage: f1 = x^-2*y^3 - 9 - 1/14*x^-1*y - 1/3*x^-1
            sage: h = 3*x^-1 - 3*x^-2*y - 1/2*x^-3*y^2 - x^-3*y + x^-3
            sage: f2 = f1 * h
            sage: f3 = f2 + x * y
            sage: f1.divides(f2)
            True
            sage: f1.divides(f3)
            False
            sage: f1.divides(3)
            False
        """
        p = self.monomial_reduction()[0]
        q = other.monomial_reduction()[0]
        return p.divides(q)
