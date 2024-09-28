r"""
Base class for elements of multivariate polynomial rings
"""

# ********************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ********************************************************************

from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.structure.coerce cimport coercion_model
from sage.misc.derivative import multi_derivative
from itertools import chain

from sage.misc.misc_c import prod


def is_MPolynomial(x):
    from sage.misc.superseded import deprecation
    deprecation(32709, "the function is_MPolynomial is deprecated; use isinstance(x, sage.rings.polynomial.multi_polynomial.MPolynomial) instead")

    return isinstance(x, MPolynomial)


from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.categories.map cimport Map
from sage.rings.rational_field import QQ

from sage.rings.polynomial.polydict cimport ETuple
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class MPolynomial(CommutativePolynomial):

    # -------------------------
    # Some standard conversions
    # -------------------------
    def _scalar_conversion(self, R):
        r"""
        TESTS::

            sage: # needs sage.rings.real_mpfr
            sage: ZZ(RR['x,y'](0))  # indirect doctest
            0
            sage: ZZ(RR['x,y'](0.5))
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
            sage: ZZ(RR['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Integer Ring

            sage: # needs sage.rings.real_mpfr
            sage: RR(RR['x,y'](0))  # indirect doctest
            0.000000000000000
            sage: RR(ZZ['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Real Field with 53 bits of precision

            sage: # needs sage.rings.real_mpfr
            sage: CC(RR['x,y'](0))  # indirect doctest
            0.000000000000000
            sage: CC(ZZ['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Complex Field with 53 bits of precision

            sage: # needs sage.rings.real_mpfr
            sage: RDF(RR['x,y'](0))
            0.0
            sage: RDF(ZZ['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Real Double Field

            sage: # needs sage.rings.real_mpfr
            sage: CDF(RR['x,y'](0))  # indirect doctest
            0.0
            sage: CDF(ZZ['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Complex Double Field

            sage: # needs sage.libs.flint sage.rings.real_mpfr
            sage: a = RR['x,y'](1)
            sage: RBF(a)
            1.000000000000000
            sage: RIF(a)
            1
            sage: CBF(a)
            1.000000000000000
            sage: CIF(a)
            1
            sage: CBF(RR['x,y'](1))  # indirect doctest
            1.000000000000000
            sage: CBF(ZZ['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Complex ball field with 53 bits of precision

            sage: x = polygen(QQ)
            sage: A.<u> = NumberField(x^3 - 2)                                          # needs sage.rings.number_field
            sage: A(A['x,y'](u))                                                        # needs sage.rings.number_field
            u
        """
        if self.degree() <= 0:
            return R(self.constant_coefficient())
        raise TypeError(f"unable to convert non-constant polynomial {self} to {R}")

    _real_double_ = _scalar_conversion
    _complex_double_ = _scalar_conversion
    _mpfr_ = _scalar_conversion
    _complex_mpfr_ = _scalar_conversion
    _real_mpfi_ = _scalar_conversion
    _complex_mpfi_ = _scalar_conversion
    _arb_ = _scalar_conversion
    _acb_ = _scalar_conversion
    _integer_ = _scalar_conversion
    _algebraic_ = _scalar_conversion
    _number_field_ = _scalar_conversion

    def __int__(self):
        r"""
        TESTS::

            sage: type(RR['x,y'])
            <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_domain_with_category'>
            sage: type(RR['x, y'](0))
            <class 'sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict'>

            sage: int(RR['x,y'](0))  # indirect doctest
            0
            sage: int(RR['x,y'](10))
            10
            sage: int(ZZ['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to <class 'int'>

            sage: ZZ(RR['x,y'](0))  # indirect doctest
            0
            sage: ZZ(RR['x,y'](0.5))                                                    # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
            sage: ZZ(RR['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Integer Ring
        """
        return self._scalar_conversion(int)

    def __float__(self):
        r"""
        TESTS::

            sage: float(RR['x,y'](0))  # indirect doctest
            0.0
            sage: float(ZZ['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to <class 'float'>
        """
        return self._scalar_conversion(float)

    def _rational_(self):
        r"""
        TESTS::

            sage: QQ(RR['x,y'](0.5))  # indirect doctest
            1/2
            sage: QQ(RR['x,y'].gen(0))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Rational Field
        """
        from sage.rings.rational_field import QQ
        return self._scalar_conversion(QQ)

    def _symbolic_(self, R):
        r"""
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: R.<x,y> = QQ[]
            sage: f = x^3 + y
            sage: g = f._symbolic_(SR); g
            x^3 + y
            sage: g(x=2, y=2)
            10
            sage: g = SR(f)
            sage: g(x=2, y=2)
            10
        """
        d = dict([(repr(g), R.var(g)) for g in self.parent().gens()])
        return self.subs(**d)

    def _polynomial_(self, R):
        var = R.variable_name()
        if var in self._parent.variable_names():
            return R(self.polynomial(self._parent(var)))
        return R([self])

    def coefficients(self):
        r"""
        Return the nonzero coefficients of this polynomial in a list.

        The returned list is decreasingly ordered by the term ordering
        of ``self.parent()``, i.e. the list of coefficients matches the list
        of monomials returned by
        :meth:`sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular.monomials`.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='degrevlex')
            sage: f = 23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]
            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: f = 23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

        Test the same stuff with base ring `\ZZ` -- different implementation::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 3, order='degrevlex')
            sage: f = 23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]
            sage: R.<x,y,z> = PolynomialRing(ZZ, 3, order='lex')
            sage: f = 23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

        AUTHOR:

        - Didier Deshommes
        """
        d = self.dict()
        return [d[i] for i in self.exponents()]

    def truncate(self, var, n):
        r"""
        Return a new multivariate polynomial obtained from ``self`` by
        deleting all terms that involve the given variable to a power
        at least ``n``.
        """
        cdef int ind
        R = self.parent()
        G = R.gens()
        Z = list(G)
        try:
            ind = Z.index(var)
        except ValueError:
            raise ValueError("var must be one of the generators of the parent polynomial ring.")
        d = self.dict()
        return R(dict([(k, c) for k, c in d.iteritems() if k[ind] < n]))

    def _fast_callable_(self, etb):
        r"""
        Given an :class:`ExpressionTreeBuilder`, return an :class:`Expression` representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y','z'])
            sage: K.<x,y,z> = QQ[]
            sage: v = -6/5*x*y*z + 2*y*z^2 - x
            sage: v._fast_callable_(etb)
            add(add(add(0, mul(-6/5, mul(mul(ipow(v_0, 1), ipow(v_1, 1)), ipow(v_2, 1)))), mul(2, mul(ipow(v_1, 1), ipow(v_2, 2)))), mul(-1, ipow(v_0, 1)))

        TESTS::

            sage: v = K(0)
            sage: vf = fast_callable(v)
            sage: type(v(0r, 0r, 0r))
            <class 'sage.rings.rational.Rational'>
            sage: type(vf(0r, 0r, 0r))
            <class 'sage.rings.rational.Rational'>
            sage: K.<x,y,z> = QQ[]
            sage: from sage.ext.fast_eval import fast_float
            sage: fast_float(K(0)).op_list()
            [('load_const', 0.0), 'return']
            sage: fast_float(K(17)).op_list()
            [('load_const', 0.0), ('load_const', 17.0), 'add', 'return']
            sage: fast_float(y).op_list()
            [('load_const', 0.0), ('load_const', 1.0), ('load_arg', 1), ('ipow', 1), 'mul', 'add', 'return']
        """
        my_vars = self.parent().variable_names()
        x = [etb.var(v) for v in my_vars]
        n = len(x)

        expr = etb.constant(self.base_ring().zero())
        for (m, c) in self.dict().iteritems():
            monom = prod([ x[i]**m[i] for i in range(n) if m[i] != 0],
                             etb.constant(c))
            expr = expr + monom
        return expr

    def derivative(self, *args):
        r"""
        The formal derivative of this polynomial, with respect to
        variables supplied in ``args``.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global function :func:`derivative` for more details.

        .. SEEALSO:: :meth:`._derivative`

        EXAMPLES:

        Polynomials implemented via Singular::

            sage: # needs sage.libs.singular
            sage: R.<x, y> = PolynomialRing(FiniteField(5))
            sage: f = x^3*y^5 + x^7*y
            sage: type(f)
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: f.derivative(x)
            2*x^6*y - 2*x^2*y^5
            sage: f.derivative(y)
            x^7

        Generic multivariate polynomials::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: S.<x, y> = PolynomialRing(R)
            sage: f = (t^2 + O(t^3))*x^2*y^3 + (37*t^4 + O(t^5))*x^3
            sage: type(f)
            <class 'sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict'>
            sage: f.derivative(x)       # with respect to x
            (2*t^2 + O(t^3))*x*y^3 + (111*t^4 + O(t^5))*x^2
            sage: f.derivative(y)       # with respect to y
            (3*t^2 + O(t^3))*x^2*y^2
            sage: f.derivative(t)       # with respect to t (recurses into base ring)
            (2*t + O(t^2))*x^2*y^3 + (148*t^3 + O(t^4))*x^3
            sage: f.derivative(x, y)    # with respect to x and then y
            (6*t^2 + O(t^3))*x*y^2
            sage: f.derivative(y, 3)    # with respect to y three times
            (6*t^2 + O(t^3))*x^2
            sage: f.derivative()        # can't figure out the variable
            Traceback (most recent call last):
            ...
            ValueError: must specify which variable to differentiate with respect to

        Polynomials over the symbolic ring (just for fun....)::

            sage: # needs sage.symbolic
            sage: x = var("x")
            sage: S.<u, v> = PolynomialRing(SR)
            sage: f = u*v*x
            sage: f.derivative(x) == u*v
            True
            sage: f.derivative(u) == v*x
            True
        """
        return multi_derivative(self, args)

    def polynomial(self, var):
        r"""
        Let ``var`` be one of the variables of the parent of ``self``.  This
        returns ``self`` viewed as a univariate polynomial in ``var`` over the
        polynomial ring generated by all the other variables of the parent.

        EXAMPLES::

            sage: R.<x,w,z> = QQ[]
            sage: f = x^3 + 3*w*x + w^5 + (17*w^3)*x + z^5
            sage: f.polynomial(x)
            x^3 + (17*w^3 + 3*w)*x + w^5 + z^5
            sage: parent(f.polynomial(x))
            Univariate Polynomial Ring in x
             over Multivariate Polynomial Ring in w, z over Rational Field

            sage: f.polynomial(w)
            w^5 + 17*x*w^3 + 3*x*w + z^5 + x^3
            sage: f.polynomial(z)
            z^5 + w^5 + 17*x*w^3 + x^3 + 3*x*w
            sage: R.<x,w,z,k> = ZZ[]
            sage: f = x^3 + 3*w*x + w^5 + (17*w^3)*x + z^5 +x*w*z*k + 5
            sage: f.polynomial(x)
            x^3 + (17*w^3 + w*z*k + 3*w)*x + w^5 + z^5 + 5
            sage: f.polynomial(w)
            w^5 + 17*x*w^3 + (x*z*k + 3*x)*w + z^5 + x^3 + 5
            sage: f.polynomial(z)
            z^5 + x*w*k*z + w^5 + 17*x*w^3 + x^3 + 3*x*w + 5
            sage: f.polynomial(k)
            x*w*z*k + w^5 + z^5 + 17*x*w^3 + x^3 + 3*x*w + 5
            sage: R.<x,y> = GF(5)[]
            sage: f = x^2 + x + y
            sage: f.polynomial(x)
            x^2 + x + y
            sage: f.polynomial(y)
            y + x^2 + x
        """
        cdef int ind
        R = self._parent
        cdef list Z = list(R.gens())
        cdef Py_ssize_t i
        cdef dict c, w
        cdef list v
        try:
            ind = Z.index(var)
        except ValueError:
            raise ValueError("var must be one of the generators of the parent polynomial ring")

        if len(Z) <= 1:
            return self.univariate_polynomial()

        del Z[ind]

        # Make polynomial ring over all variables except var.
        S = PolynomialRing(R.base_ring(), Z)
        ring = S[var]
        if not self:
            return ring(0)

        d = self.degree(var)
        B = ring.base_ring()
        w = {remove_from_tuple(e, ind): val
             for e, val in self.dict().iteritems() if not e[ind]}
        v = [B(w)]  # coefficients that don't involve var
        z = var
        for i in range(1,d+1):
            c = <dict> self.coefficient(z).dict()
            w = {remove_from_tuple(e, ind): val for e, val in c.iteritems()}
            v.append(B(w))
            z *= var
        return ring(v)

    cpdef dict _mpoly_dict_recursive(self, tuple vars=None, base_ring=None):
        r"""
        Return a ``dict`` of coefficient entries suitable for construction
        of a ``MPolynomial_polydict`` with the given variables.

        EXAMPLES::

            sage: R = Integers(10)['x,y,z']['t,s']
            sage: t,s = R.gens()
            sage: x,y,z = R.base_ring().gens()
            sage: (x+y+2*z*s+3*t)._mpoly_dict_recursive(('z','t','s'))
            {(0, 0, 0): x + y, (0, 1, 0): 3, (1, 0, 1): 2}

        TESTS::

            sage: # needs sage.rings.padics
            sage: R = Qp(7)['x,y,z,t,p']; S = ZZ['x,z,t']['p']
            sage: R(S.0)
            p
            sage: R = QQ['x,y,z,t,p']; S = ZZ['x']['y,z,t']['p']
            sage: z = S.base_ring().gen(1)
            sage: R(z)
            z
            sage: R = QQ['x,y,z,t,p']; S = ZZ['x']['y,z,t']['p']
            sage: z = S.base_ring().gen(1); p = S.0; x = S.base_ring().base_ring().gen()
            sage: R(z+p)
            z + p
            sage: R = Qp(7)['x,y,z,p']; S = ZZ['x']['y,z,t']['p']  # shouldn't work, but should throw a better error
            sage: R(S.0)
            p

        See :issue:`2601`::

            sage: R.<a,b,c> = PolynomialRing(QQ, 3)
            sage: a._mpoly_dict_recursive(('c', 'b', 'a'))
            {(0, 0, 1): 1}
            sage: testR.<a,b,c> = PolynomialRing(QQ,3)
            sage: id_ringA = ideal([a^2 - b, b^2 - c, c^2 - a])
            sage: id_ringB = ideal(id_ringA.gens()).change_ring(PolynomialRing(QQ,'c,b,a'))
        """
        if not self:
            return {}

        if vars is None:
            vars = self._parent.variable_names_recursive()
        cdef tuple my_vars = self._parent.variable_names()
        if vars == my_vars:
            return <dict> self.dict()
        elif my_vars[-1] not in vars:
            x = base_ring(self) if base_ring is not None else self
            const_ix = ETuple((0,)*len(vars))
            return { const_ix: x }
        elif not set(my_vars).issubset(set(vars)):
            # we need to split it up
            p = self.polynomial(self._parent.gen(len(my_vars)-1))
            if not isinstance(p, MPolynomial):
                # Not a multivariate polynomial, so it must be a univariate
                return (<Polynomial> p)._mpoly_dict_recursive(vars, base_ring)
            return (<MPolynomial> p)._mpoly_dict_recursive(vars, base_ring)

        cdef dict D = {}
        cdef list mapping = [vars.index(z) for z in my_vars]
        cdef list new_map
        cdef Py_ssize_t m = min(mapping)
        cdef tuple prev_vars = vars[:m]
        cdef list tmp
        cdef ETuple postfix
        cdef Py_ssize_t k
        cdef dict mpoly
        if prev_vars:
            new_map = list(mapping)
            for k in range(len(mapping)):
                new_map[k] -= m
            tmp = [0] * (len(vars) - m)
            try:
                for ix,a in self.dict().iteritems():
                    for k in range(len(my_vars)):
                        tmp[new_map[k]] = ix[k]
                    postfix = ETuple(tmp)
                    mpoly = <dict> a._mpoly_dict_recursive(prev_vars, base_ring)
                    for prefix,b in mpoly.iteritems():
                        D[prefix+postfix] = b
                return D

            except AttributeError:
                pass

        if base_ring is self.base_ring():
            base_ring = None

        tmp = [0] * len(vars)
        for ix,a in self.dict().iteritems():
            for k in range(len(my_vars)):
                tmp[mapping[k]] = ix[k]
            if base_ring is not None:
                a = base_ring(a)
            D[ETuple(tmp)] = a
        return D

    cdef long _hash_c(self) except -1:
        r"""
        This hash incorporates the variable name in an effort to respect the obvious inclusions
        into multi-variable polynomial rings.

        The tuple algorithm is borrowed from http://effbot.org/zone/python-hash.htm.

        EXAMPLES::

            sage: T.<y>=QQ[]
            sage: R.<x>=ZZ[]
            sage: S.<x,y>=ZZ[]
            sage: hash(S.0)==hash(R.0)  # respect inclusions into mpoly rings (with matching base rings)
            True
            sage: hash(S.1)==hash(T.0)  # respect inclusions into mpoly rings (with unmatched base rings)
            True
            sage: hash(S(12))==hash(12)  # respect inclusions of the integers into an mpoly ring
            True
            sage: # the point is to make for more flexible dictionary look ups
            sage: d={S.0:12}
            sage: d[R.0]
            12
            sage: # or, more to the point, make subs in fraction field elements work
            sage: f=x/y
            sage: f.subs({x:1})
            1/y

        TESTS:

        Verify that :issue:`16251` has been resolved, i.e., polynomials with
        unhashable coefficients are unhashable::

            sage: K.<a> = Qq(9)                                                         # needs sage.rings.padics
            sage: R.<t,s> = K[]                                                         # needs sage.rings.padics
            sage: hash(t)                                                               # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'sage.rings.padics.qadic_flint_CR.qAdicCappedRelativeElement'
        """
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        var_name_hash = [hash(v) for v in self._parent.variable_names()]
        cdef long c_hash
        for m,c in self.dict().iteritems():
            #  I'm assuming (incorrectly) that hashes of zero indicate that the element is 0.
            # This assumption is not true, but I think it is true enough for the purposes and it
            # it allows us to write fast code that omits terms with 0 coefficients.  This is
            # important if we want to maintain the '==' relationship with sparse polys.
            c_hash = hash(c)
            if c_hash != 0: # this is always going to be true, because we are sparse (correct?)
                # Hash (self[i], gen_a, exp_a, gen_b, exp_b, gen_c, exp_c, ...) as a tuple according to the algorithm.
                # I omit gen,exp pairs where the exponent is zero.
                result_mon = c_hash
                for p in m.nonzero_positions():
                    result_mon = (1000003 * result_mon) ^ var_name_hash[p]
                    result_mon = (1000003 * result_mon) ^ m[p]
                result += result_mon
        if result == -1:
            return -2
        return result

    # you may have to replicate this boilerplate code in derived classes if you override
    # __richcmp__.  The python documentation at  https://docs.python.org/api/type-structs.html
    # explains how __richcmp__, __hash__, and __cmp__ are tied together.
    def __hash__(self):
        return self._hash_c()

    def args(self):
        r"""
        Return the names of the arguments of ``self``, in the
        order they are accepted from call.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: x.args()
            (x, y)
        """
        return self._parent.gens()

    def homogenize(self, var='h'):
        r"""
        Return the homogenization of this polynomial.

        The polynomial itself is returned if it is homogeneous already.
        Otherwise, the monomials are multiplied with the smallest powers of
        ``var`` such that they all have the same total degree.

        INPUT:

        - ``var`` -- a variable in the polynomial ring (as a string, an element of
          the ring, or a zero-based index in the list of variables) or a name
          for a new variable (default: ``'h'``)

        OUTPUT:

        If ``var`` specifies a variable in the polynomial ring, then a
        homogeneous element in that ring is returned. Otherwise, a homogeneous
        element is returned in a polynomial ring with an extra last variable
        ``var``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = x^2 + y + 1 + 5*x*y^10
            sage: f.homogenize()
            5*x*y^10 + x^2*h^9 + y*h^10 + h^11

        The parameter ``var`` can be used to specify the name of the variable::

            sage: g = f.homogenize('z'); g
            5*x*y^10 + x^2*z^9 + y*z^10 + z^11
            sage: g.parent()
            Multivariate Polynomial Ring in x, y, z over Rational Field

        However, if the polynomial is homogeneous already, then that parameter
        is ignored and no extra variable is added to the polynomial ring::

            sage: f = x^2 + y^2
            sage: g = f.homogenize('z'); g
            x^2 + y^2
            sage: g.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        If you want the ring of the result to be independent of whether the
        polynomial is homogenized, you can use ``var`` to use an existing
        variable to homogenize::

            sage: R.<x,y,z> = QQ[]
            sage: f = x^2 + y^2
            sage: g = f.homogenize(z); g
            x^2 + y^2
            sage: g.parent()
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: f = x^2 - y
            sage: g = f.homogenize(z); g
            x^2 - y*z
            sage: g.parent()
            Multivariate Polynomial Ring in x, y, z over Rational Field

        The parameter ``var`` can also be given as a zero-based index in the
        list of variables::

            sage: g = f.homogenize(2); g
            x^2 - y*z

        If the variable specified by ``var`` is not present in the polynomial,
        then setting it to 1 yields the original polynomial::

            sage: g(x,y,1)
            x^2 - y

        If it is present already, this might not be the case::

            sage: g = f.homogenize(x); g
            x^2 - x*y
            sage: g(1,y,z)
            -y + 1

        In particular, this can be surprising in positive characteristic::

            sage: R.<x,y> = GF(2)[]
            sage: f = x + 1
            sage: f.homogenize(x)
            0

        TESTS::

            sage: R = PolynomialRing(QQ, 'x', 5)
            sage: p = R.random_element()
            sage: q1 = p.homogenize()
            sage: q2 = p.homogenize()
            sage: q1.parent() is q2.parent()
            True
        """
        P = self.parent()

        if self.is_homogeneous():
            return self

        if isinstance(var, str):
            V = list(P.variable_names())
            try:
                i = V.index(var)
                return self._homogenize(i)
            except ValueError:
                P = PolynomialRing(P.base_ring(), len(V)+1, V + [var], order=P.term_order())
                return P(self)._homogenize(len(V))

        elif isinstance(var, MPolynomial) and \
             ((<MPolynomial>var)._parent is P or (<MPolynomial>var)._parent == P):
            V = list(P.gens())
            try:
                i = V.index(var)
                return self._homogenize(i)
            except ValueError:
                P = P.change_ring(names=P.variable_names() + [str(var)])
                return P(self)._homogenize(len(V))

        elif isinstance(var, int) or isinstance(var, Integer):
            if 0 <= var < P.ngens():
                return self._homogenize(var)
            else:
                raise TypeError("Variable index %d must be < parent(self).ngens()." % var)
        else:
            raise TypeError("Parameter var must be either a variable, a string or an integer.")

    def is_homogeneous(self):
        r"""
        Return ``True`` if ``self`` is a homogeneous polynomial.

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial import MPolynomial
            sage: P.<x, y> = PolynomialRing(QQ, 2)
            sage: MPolynomial.is_homogeneous(x+y)
            True
            sage: MPolynomial.is_homogeneous(P(0))
            True
            sage: MPolynomial.is_homogeneous(x+y^2)
            False
            sage: MPolynomial.is_homogeneous(x^2 + y^2)
            True
            sage: MPolynomial.is_homogeneous(x^2 + y^2*x)
            False
            sage: MPolynomial.is_homogeneous(x^2*y + y^2*x)
            True

        .. NOTE::

            This is a generic implementation which is likely overridden by
            subclasses.
        """
        M = self.monomials()
        if M==[]:
            return True
        d = M.pop().degree()
        for m in M:
            if m.degree() != d:
                return False
        else:
            return True

    def homogeneous_components(self):
        r"""
        Return the homogeneous components of this polynomial.

        OUTPUT: a dictionary mapping degrees to homogeneous polynomials

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: (x^3 + 2*x*y^3 + 4*y^3 + y).homogeneous_components()
            {1: y, 3: x^3 + 4*y^3, 4: 2*x*y^3}
            sage: R.zero().homogeneous_components()
            {}

        In case of weighted term orders, the polynomials are homogeneous with
        respect to the weights::

             sage: S.<a,b,c> = PolynomialRing(ZZ, order=TermOrder('wdegrevlex', (1,2,3)))
             sage: (a^6 + b^3 + b*c + a^2*c + c + a + 1).homogeneous_components()
             {0: 1, 1: a, 3: c, 5: a^2*c + b*c, 6: a^6 + b^3}
        """
        cdef ETuple e
        from collections import defaultdict
        d = defaultdict(dict)
        if self._parent.term_order()._weights:
            for c, m in self:
                d[m.degree()][m.exponents()[0]] = c
        else:
            # Otherwise it is unweighted, so we use a faster implementation
            for e, c in self.iterator_exp_coeff():
                d[e.unweighted_degree()][e] = c
        return {k: self._parent(d[k]) for k in d}

    cpdef _mod_(self, other):
        r"""
        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = (x^2*y + 2*x - 3)
            sage: g = (x + 1)*f
            sage: g % f                                                                 # needs sage.libs.singular
            0

            sage: (g+1) % f                                                             # needs sage.libs.singular
            1

            sage: M = x*y
            sage: N = x^2*y^3
            sage: M.divides(N)                                                          # needs sage.libs.singular
            True
        """
        try:
            quo_rem = self.quo_rem
        except AttributeError:
            raise NotImplementedError
        else:
            q, r = quo_rem(other)
            return r

    def change_ring(self, R):
        r"""
        Return this polynomial with coefficients converted to ``R``.

        INPUT:

        - ``R`` -- a ring or morphism; if a morphism, the coefficients
          are mapped to the codomain of ``R``

        OUTPUT: a new polynomial with the base ring changed to ``R``

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = x^3 + 3/5*y + 1
            sage: f.change_ring(GF(7))
            x^3 + 2*y + 1
            sage: g = x^2 + 5*y
            sage: g.change_ring(GF(5))
            x^2

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = GF(9,'a')[]
            sage: (x+2*y).change_ring(GF(3))
            x - y

        ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(7^2)
            sage: R.<x,y> = F[]
            sage: f = x^2 + a^2*y^2 + a*x + a^3*y
            sage: g = f.change_ring(F.frobenius_endomorphism()); g
            x^2 + (-a - 2)*y^2 + (-a + 1)*x + (2*a + 2)*y
            sage: g.change_ring(F.frobenius_endomorphism()) == f
            True

        ::

            sage: # needs sage.rings.number_field
            sage: K.<z> = CyclotomicField(3)
            sage: R.<x,y> = K[]
            sage: f = x^2 + z*y
            sage: f.change_ring(K.embeddings(CC)[1])
            x^2 + (-0.500000000000000 - 0.866025403784438*I)*y

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = CyclotomicField(5)
            sage: R.<x,y> = K[]
            sage: f = x^2 + w*y
            sage: f.change_ring(K.embeddings(QQbar)[1])
            x^2 + (-0.8090169943749474? + 0.5877852522924731?*I)*y

        TESTS:

        Check that :issue:`25022` is fixed::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: K.<x,y> = ZZ[]
            sage: (x*y).change_ring(SR).monomials()
            [x*y]

        Check that :issue:`36832` is fixed::

            sage: F = GF(11)
            sage: phi = Hom(F,F).an_element()
            sage: R.<x,y> = F[]
            sage: x.change_ring(phi)
            x
        """
        if isinstance(R, Map):
            return self.map_coefficients(R)
        return self.parent().change_ring(R)(self.dict())

    def is_symmetric(self, group=None):
        r"""
        Return whether this polynomial is symmetric.

        INPUT:

        - ``group`` -- (default: symmetric group) if set, test whether the
          polynomial is invariant with respect to the given permutation group

        EXAMPLES::

            sage: # needs sage.groups
            sage: R.<x,y,z> = QQ[]
            sage: p = (x+y+z)**2 - 3 * (x+y)*(x+z)*(y+z)
            sage: p.is_symmetric()
            True
            sage: (x + y - z).is_symmetric()
            False
            sage: R.one().is_symmetric()
            True
            sage: p = (x-y)*(y-z)*(z-x)
            sage: p.is_symmetric()
            False
            sage: p.is_symmetric(AlternatingGroup(3))
            True

            sage: R.<x,y> = QQ[]
            sage: ((x + y)**2).is_symmetric()                                           # needs sage.groups
            True
            sage: R.one().is_symmetric()                                                # needs sage.groups
            True
            sage: (x + 2*y).is_symmetric()                                              # needs sage.groups
            False

        An example with a GAP permutation group (here the quaternions)::

            sage: R = PolynomialRing(QQ, 'x', 8)
            sage: x = R.gens()
            sage: p = sum(prod(x[i] for i in e)
            ....:         for e in [(0,1,2), (0,1,7), (0,2,7), (1,2,7),
            ....:                   (3,4,5), (3,4,6), (3,5,6), (4,5,6)])
            sage: p.is_symmetric(libgap.TransitiveGroup(8, 5))                          # needs sage.groups
            True
            sage: p = sum(prod(x[i] for i in e)
            ....:     for e in [(0,1,2), (0,1,7), (0,2,7), (1,2,7),
            ....:               (3,4,5), (3,4,6), (3,5,6)])
            sage: p.is_symmetric(libgap.TransitiveGroup(8, 5))                          # needs sage.groups
            False

        TESTS::

            sage: R = PolynomialRing(QQ, 'x', 3)
            sage: R.one().is_symmetric(3)                                               # needs sage.groups
            Traceback (most recent call last):
            ...
            ValueError: argument must be a permutation group

            sage: R.one().is_symmetric(SymmetricGroup(4))                               # needs sage.groups
            Traceback (most recent call last):
            ...
            ValueError: invalid data to initialize a permutation
        """
        n = self.parent().ngens()
        if n <= 1:
            return True

        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        S = SymmetricGroup(n)
        if group is None:
            gens = S.gens()
        else:
            try:
                # for Sage group
                gens = group.gens()
            except AttributeError:
                # for GAP group
                try:
                    gens = group.GeneratorsOfGroup()
                except AttributeError:
                    raise ValueError("argument must be a permutation group")
            gens = [S(g) for g in gens]

        cdef dict coeffs = self.dict()
        zero = self.base_ring().zero()
        return all(coeffs.get(g._act_on_etuple_on_position(e), zero) == coeff
                   for e, coeff in coeffs.items() for g in gens)

    def _gap_(self, gap):
        r"""
        Return a representation of ``self`` in the GAP interface.

        INPUT:

        - ``gap`` -- a GAP or libgap instance

        TESTS:

        Multivariate polynomial over integers::

            sage: # needs sage.libs.gap
            sage: R.<x,y,z> = ZZ[]
            sage: gap(-x*y + 3*z)   # indirect doctest
            -x*y+3*z
            sage: gap(R.zero())     # indirect doctest
            0
            sage: (x+y+z)._gap_(libgap)
            x+y+z

            sage: g = gap(x - y + 3*x*y*z)                                              # needs sage.libs.gap
            sage: R(g)                                                                  # needs sage.libs.gap
            3*x*y*z + x - y

            sage: g = libgap(5*x - y*z)                                                 # needs sage.libs.gap
            sage: R(g)                                                                  # needs sage.libs.gap
            -y*z + 5*x

        Multivariate polynomial over a cyclotomic field::

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: F.<zeta> = CyclotomicField(8)
            sage: P.<x,y> = F[]
            sage: p = zeta + zeta^2*x + zeta^3*y + (1+zeta)*x*y
            sage: gap(p)     # indirect doctest
            (1+E(8))*x*y+E(4)*x+E(8)^3*y+E(8)
            sage: libgap(p)  # indirect doctest
            (1+E(8))*x*y+E(4)*x+E(8)^3*y+E(8)

        Multivariate polynomial over a polynomial ring over a cyclotomic field::

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: S.<z> = F[]
            sage: P.<x,y> = S[]
            sage: p = zeta + zeta^2*x*z + zeta^3*y*z^2 + (1+zeta)*x*y*z
            sage: gap(p)     # indirect doctest
            ((1+E(8))*z)*x*y+E(4)*z*x+E(8)^3*z^2*y+E(8)
            sage: libgap(p)  # indirect doctest
            ((1+E(8))*z)*x*y+E(4)*z*x+E(8)^3*z^2*y+E(8)
        """
        R = gap(self.parent())
        variables = R.IndeterminatesOfPolynomialRing()
        return self(*variables)

    def _libgap_(self):
        r"""
        TESTS::

            sage: R.<x,y,z> = ZZ[]
            sage: libgap(-x*y + 3*z)   # indirect doctest                               # needs sage.libs.gap
            -x*y+3*z
            sage: libgap(R.zero())     # indirect doctest                               # needs sage.libs.gap
            0
        """
        from sage.libs.gap.libgap import libgap
        return self._gap_(libgap)

    def _magma_init_(self, magma):
        r"""
        Return a Magma string representation of ``self`` valid in the
        given magma session.

        EXAMPLES::

            sage: # optional - magma, needs sage.rings.finite_rings
            sage: k.<b> = GF(25); R.<x,y> = k[]
            sage: f = y*x^2*b + x*(b+1) + 1
            sage: magma = Magma()                       # so var names same below
            sage: magma(f)
            b*x^2*y + b^22*x + 1
            sage: f._magma_init_(magma)
            '_sage_[...]!((_sage_[...]!(_sage_[...]))*_sage_[...]^2*_sage_[...]+(_sage_[...]!(_sage_[...] + 1))*_sage_[...]+(_sage_[...]!(1))*1)'

        A more complicated nested example::

            sage: # optional - magma
            sage: R.<x,y> = QQ[]; S.<z,w> = R[]; f = (2/3)*x^3*z + w^2 + 5
            sage: f._magma_init_(magma)
            '_sage_[...]!((_sage_[...]!((1/1)*1))*_sage_[...]^2+(_sage_[...]!((2/3)*_sage_[...]^3))*_sage_[...]+(_sage_[...]!((5/1)*1))*1)'
            sage: magma(f)
            w^2 + 2/3*x^3*z + 5
        """
        R = magma(self.parent())
        g = R.gen_names()
        v = []
        for m, c in zip(self.monomials(), self.coefficients()):
            v.append('(%s)*%s'%( c._magma_init_(magma),
                                 m._repr_with_changed_varnames(g)))
        if len(v) == 0:
            s = '0'
        else:
            s = '+'.join(v)

        return '%s!(%s)'%(R.name(), s)

    def _giac_init_(self):
        r"""
        Return a Giac string representation of this polynomial.

        TESTS::

            sage: # needs sage.libs.giac
            sage: R.<x,y,z> = GF(101)['e,i'][]
            sage: f = R('e*i') * x + y^2
            sage: f._giac_init_()
            '((1)*1)*sageVARy^2+((1)*sageVARe*sageVARi)*sageVARx'
            sage: giac(f)
            sageVARy^2+sageVARe*sageVARi*sageVARx
            sage: giac(R.zero())
            0
        """
        g = ['sageVAR' + x for x in self.parent().variable_names()]
        s = '+'.join('(%s)*%s' % (c._giac_init_(),
                                  m._repr_with_changed_varnames(g))
                     for c, m in self)
        return s if s else '0'

    def gradient(self):
        r"""
        Return a list of partial derivatives of this polynomial,
        ordered by the variables of ``self.parent()``.

        EXAMPLES::

           sage: P.<x,y,z> = PolynomialRing(ZZ, 3)
           sage: f = x*y + 1
           sage: f.gradient()
           [y, x, 0]
        """
        return [ self.derivative(var) for var in self.parent().gens() ]

    def jacobian_ideal(self):
        r"""
        Return the Jacobian ideal of the polynomial ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = x^3 + y^3 + z^3
            sage: f.jacobian_ideal()
            Ideal (3*x^2, 3*y^2, 3*z^2) of
             Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.parent().ideal(self.gradient())

    def newton_polytope(self):
        r"""
        Return the Newton polytope of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = 1 + x*y + x^3 + y^3
            sage: P = f.newton_polytope(); P                                            # needs sage.geometry.polyhedron
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: P.is_simple()                                                         # needs sage.geometry.polyhedron
            True

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: R(0).newton_polytope()                                                # needs sage.geometry.polyhedron
            The empty polyhedron in ZZ^0
            sage: R(1).newton_polytope()                                                # needs sage.geometry.polyhedron
            A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: R(x^2+y^2).newton_polytope().integral_points()                        # needs sage.geometry.polyhedron
            ((0, 2), (1, 1), (2, 0))
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        e = self.exponents()
        P = Polyhedron(vertices = e, base_ring=ZZ)
        return P

    def __iter__(self):
        r"""
        Facilitates iterating over the monomials of ``self``,
        returning tuples of the form ``(coeff, mon)`` for each
        nonzero monomial.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3)
            sage: f = 3*x^3*y + 16*x + 7
            sage: [(c,m) for c,m in f]
            [(3, x^3*y), (16, x), (7, 1)]
            sage: f = P.random_element(12,14)
            sage: sum(c*m for c,m in f) == f
            True
        """
        for exp, coeff in self.iterator_exp_coeff():
            yield (coeff, self.monomial(exp))

    def iterator_exp_coeff(self, as_ETuples=True):
        r"""
        Iterate over ``self`` as pairs of ((E)Tuple, coefficient).

        INPUT:

        - ``as_ETuples`` -- boolean (default: ``True``); if ``True``, iterate over
          pairs whose first element is an :class:`ETuple`, otherwise as a tuples

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: f = a*c^3 + a^2*b + 2*b^4
            sage: list(f.iterator_exp_coeff())
            [((0, 4, 0), 2), ((1, 0, 3), 1), ((2, 1, 0), 1)]
            sage: list(f.iterator_exp_coeff(as_ETuples=False))
            [((0, 4, 0), 2), ((1, 0, 3), 1), ((2, 1, 0), 1)]

            sage: R.<a,b,c> = PolynomialRing(QQ, 3, order='lex')
            sage: f = a*c^3 + a^2*b + 2*b^4
            sage: list(f.iterator_exp_coeff())
            [((2, 1, 0), 1), ((1, 0, 3), 1), ((0, 4, 0), 2)]
        """
        for exp in self.exponents():
            yield (exp, self.monomial_coefficient(exp))

    def content(self):
        r"""
        Return the content of this polynomial.  Here, we define content as
        the gcd of the coefficients in the base ring.

        .. SEEALSO::

            :meth:`content_ideal`

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = 4*x+6*y
            sage: f.content()
            2
            sage: f.content().parent()
            Integer Ring

        TESTS:

        Since :issue:`10771`, the gcd in QQ restricts to the gcd in ZZ::

            sage: R.<x,y> = QQ[]
            sage: f = 4*x+6*y
            sage: f.content(); f.content().parent()
            2
            Rational Field
        """
        from sage.arith.misc import GCD as gcd
        return gcd(self.coefficients())

    def content_ideal(self):
        r"""
        Return the content ideal of this polynomial, defined as the ideal
        generated by its coefficients.

        .. SEEALSO::

            :meth:`content`

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = 2*x*y + 6*x - 4*y + 2
            sage: f.content_ideal()
            Principal ideal (2) of Integer Ring
            sage: S.<z,t> = R[]
            sage: g = x*z + y*t
            sage: g.content_ideal()
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return self.base_ring().ideal(self.coefficients())

    def is_generator(self):
        r"""
        Return ``True`` if this polynomial is a generator of its parent.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: x.is_generator()
            True
            sage: (x + y - y).is_generator()
            True
            sage: (x*y).is_generator()
            False
            sage: R.<x,y> = QQ[]
            sage: x.is_generator()
            True
            sage: (x + y - y).is_generator()
            True
            sage: (x*y).is_generator()
            False
        """
        return (self in self.parent().gens())

    def map_coefficients(self, f, new_base_ring=None):
        r"""
        Return the polynomial obtained by applying ``f`` to the nonzero
        coefficients of ``self``.

        If ``f`` is a :class:`sage.categories.map.Map`, then the resulting
        polynomial will be defined over the codomain of ``f``. Otherwise, the
        resulting polynomial will be over the same ring as ``self``. Set
        ``new_base_ring`` to override this behaviour.

        INPUT:

        - ``f`` -- a callable that will be applied to the coefficients of ``self``

        - ``new_base_ring`` -- (optional) if given, the resulting polynomial
          will be defined over this ring

        EXAMPLES::

            sage: k.<a> = GF(9); R.<x,y> = k[];  f = x*a + 2*x^3*y*a + a                # needs sage.rings.finite_rings
            sage: f.map_coefficients(lambda a: a + 1)                                   # needs sage.rings.finite_rings
            (-a + 1)*x^3*y + (a + 1)*x + (a + 1)

        Examples with different base ring::

            sage: # needs sage.rings.finite_rings
            sage: R.<r> = GF(9); S.<s> = GF(81)
            sage: h = Hom(R,S)[0]; h
            Ring morphism:
              From: Finite Field in r of size 3^2
              To:   Finite Field in s of size 3^4
              Defn: r |--> 2*s^3 + 2*s^2 + 1
            sage: T.<X,Y> = R[]
            sage: f = r*X + Y
            sage: g = f.map_coefficients(h); g
            (-s^3 - s^2 + 1)*X + Y
            sage: g.parent()
            Multivariate Polynomial Ring in X, Y over Finite Field in s of size 3^4
            sage: h = lambda x: x.trace()
            sage: g = f.map_coefficients(h); g
            X - Y
            sage: g.parent()
            Multivariate Polynomial Ring in X, Y over Finite Field in r of size 3^2
            sage: g = f.map_coefficients(h, new_base_ring=GF(3)); g
            X - Y
            sage: g.parent()
            Multivariate Polynomial Ring in X, Y over Finite Field of size 3
        """
        R = self.parent()
        if new_base_ring is not None:
            R = R.change_ring(new_base_ring)
        elif isinstance(f, Map):
            R = R.change_ring(f.codomain())
        return R(dict([(k,f(v)) for (k,v) in self.dict().items()]))

    def _norm_over_nonprime_finite_field(self):
        r"""
        Given a multivariate polynomial over a nonprime finite field
        `\GF{p^e}`, compute the norm of the polynomial down to `\GF{p}`.

        This is the product of the conjugates by the Frobenius action
        on coefficients, where Frobenius acts by `p`-th power.

        This is (currently) an internal function used in factoring over finite
        fields.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(9)
            sage: R.<x,y> = PolynomialRing(k)
            sage: f = (x-a) * (y-a)
            sage: f._norm_over_nonprime_finite_field()
            x^2*y^2 - x^2*y - x*y^2 - x^2 + x*y - y^2 + x + y + 1
        """
        P = self.parent()
        k = P.base_ring()
        if not k.is_field() and k.is_finite():
            raise TypeError("k must be a finite field")
        p = k.characteristic()
        e = k.degree()
        v = [self] + [self.map_coefficients(k.hom([k.gen()**(p**i)])) for i in range(1,e)]
        return prod(v).change_ring(k.prime_subfield())

    def sylvester_matrix(self, right, variable = None):
        r"""
        Given two nonzero polynomials ``self`` and ``right``, return the Sylvester
        matrix of the polynomials with respect to a given variable.

        Note that the Sylvester matrix is not defined if one of the polynomials
        is zero.

        INPUT:

        - ``self``, ``right`` -- multivariate polynomials
        - ``variable`` -- (optional) compute the Sylvester matrix with respect
          to this variable. If ``variable`` is not provided, the first variable
          of the polynomial ring is used.

        OUTPUT: the Sylvester matrix of ``self`` and ``right``

        EXAMPLES::

            sage: R.<x, y> = PolynomialRing(ZZ)
            sage: f = (y + 1)*x + 3*x**2
            sage: g = (y + 2)*x + 4*x**2
            sage: M = f.sylvester_matrix(g, x)                                          # needs sage.modules
            sage: M                                                                     # needs sage.modules
            [    3 y + 1     0     0]
            [    0     3 y + 1     0]
            [    4 y + 2     0     0]
            [    0     4 y + 2     0]

        If the polynomials share a non-constant common factor then the
        determinant of the Sylvester matrix will be zero::

            sage: M.determinant()                                                       # needs sage.modules
            0

            sage: f.sylvester_matrix(1 + g, x).determinant()                            # needs sage.modules
            y^2 - y + 7

        If both polynomials are of positive degree with respect to ``variable``, the
        determinant of the Sylvester matrix is the resultant::

            sage: f = R.random_element(4) or (x^2 * y^2)
            sage: g = R.random_element(4) or (x^2 * y^2)
            sage: f.sylvester_matrix(g, x).determinant() == f.resultant(g, x)           # needs sage.libs.singular sage.modules
            True

        TESTS:

        The variable is optional::

            sage: f = x + y
            sage: g = x + y
            sage: f.sylvester_matrix(g)                                                 # needs sage.modules
            [1 y]
            [1 y]

        Polynomials must be defined over compatible base rings::

            sage: # needs sage.modules
            sage: K.<x, y> = QQ[]
            sage: f = x + y
            sage: L.<x, y> = ZZ[]
            sage: g = x + y
            sage: R.<x, y> = GF(25, 'a')[]                                              # needs sage.rings.finite_rings
            sage: h = x + y
            sage: f.sylvester_matrix(g, 'x')
            [1 y]
            [1 y]
            sage: g.sylvester_matrix(h, 'x')                                            # needs sage.rings.finite_rings
            [1 y]
            [1 y]
            sage: f.sylvester_matrix(h, 'x')                                            # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Multivariate Polynomial Ring in x, y over Rational Field' and
            'Multivariate Polynomial Ring in x, y over Finite Field in a of size 5^2'
            sage: K.<x, y, z> = QQ[]
            sage: f = x + y
            sage: L.<x, z> = QQ[]
            sage: g = x + z
            sage: f.sylvester_matrix(g)
            [1 y]
            [1 z]

        Corner cases::

            sage: # needs sage.modules
            sage: K.<x, y> = QQ[]
            sage: f = x^2 + 1
            sage: g = K(0)
            sage: f.sylvester_matrix(g)
            Traceback (most recent call last):
            ...
            ValueError: The Sylvester matrix is not defined for zero polynomials
            sage: g.sylvester_matrix(f)
            Traceback (most recent call last):
            ...
            ValueError: The Sylvester matrix is not defined for zero polynomials
            sage: g.sylvester_matrix(g)
            Traceback (most recent call last):
            ...
            ValueError: The Sylvester matrix is not defined for zero polynomials
            sage: K(3).sylvester_matrix(x^2)
            [3 0]
            [0 3]
            sage: K(3).sylvester_matrix(K(4))
            []
        """

        # This code is almost exactly the same as that of
        # sylvester_matrix() in polynomial_element.pyx.

        from sage.matrix.constructor import matrix

        if self.parent() != right.parent():
            a, b = coercion_model.canonical_coercion(self,right)
            if variable:
                variable = a.parent()(variable)
            #We add the variable in case right is a multivariate polynomial
            return a.sylvester_matrix(b, variable)

        if not variable:
            variable = self.parent().gen()

        #coerce the variable to a polynomial
        if variable.parent() != self.parent():
            variable = self.parent()(variable)

        if self.is_zero() or right.is_zero():
            raise ValueError("The Sylvester matrix is not defined for zero polynomials")

        m = self.degree(variable)
        n = right.degree(variable)

        M = matrix(self.parent(), m + n, m + n)

        r = 0
        offset = 0
        for _ in range(n):
            for c in range(m, -1, -1):
                M[r, m - c + offset] = self.coefficient({variable:c})
            offset += 1
            r += 1

        offset = 0
        for _ in range(m):
            for c in range(n, -1, -1):
                M[r, n - c + offset] = right.coefficient({variable:c})
            offset += 1
            r += 1

        return M

    def discriminant(self,variable):
        r"""
        Return the discriminant of ``self`` with respect to the given variable.

        INPUT:

        - ``variable`` -- the variable with respect to which we compute
          the discriminant

        OUTPUT: an element of the base ring of the polynomial ring

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = 4*x*y^2 + 1/4*x*y*z + 3/2*x*z^2 - 1/2*z^2
            sage: f.discriminant(x)                                                     # needs sage.libs.singular
            1
            sage: f.discriminant(y)                                                     # needs sage.libs.singular
            -383/16*x^2*z^2 + 8*x*z^2
            sage: f.discriminant(z)                                                     # needs sage.libs.singular
            -383/16*x^2*y^2 + 8*x*y^2

        Note that, unlike the univariate case, the result lives in
        the same ring as the polynomial::

            sage: R.<x,y> = QQ[]
            sage: f = x^5*y + 3*x^2*y^2 - 2*x + y - 1
            sage: f.discriminant(y)                                                     # needs sage.libs.singular
            x^10 + 2*x^5 + 24*x^3 + 12*x^2 + 1
            sage: f.polynomial(y).discriminant()                                        # needs sage.libs.pari sage.modules
            x^10 + 2*x^5 + 24*x^3 + 12*x^2 + 1
            sage: f.discriminant(y).parent() == f.polynomial(y).discriminant().parent()             # needs sage.libs.singular sage.modules
            False

        TESTS:

        Test polynomials over QQbar (:issue:`25265`)::

            sage: # needs sage.rings.number_field
            sage: R.<x,y> = QQbar[]
            sage: f = x^5*y + 3*x^2*y^2 - 2*x + y - 1
            sage: f.discriminant(y)                                                     # needs sage.libs.singular
            x^10 + 2*x^5 + 24*x^3 + 12*x^2 + 1

        AUTHOR: Miguel Marco
        """
        if self.is_zero():
            return self.parent().zero()
        n = self.degree(variable)
        d = self.derivative(variable)
        k = d.degree(variable)

        r = n % 4
        u = -1 # (-1)**(n*(n-1)/2)
        if r == 0 or r == 1:
            u = 1
        an = self.coefficient(variable**n)**(n - k - 2)
        return self.parent()(u * self.resultant(d, variable) * an)

    def subresultants(self, other, variable=None):
        r"""
        Return the nonzero subresultant polynomials of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a polynomial

        OUTPUT: list of polynomials in the same ring as ``self``

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: p = (y^2 + 6)*(x - 1) - y*(x^2 + 1)
            sage: q = (x^2 + 6)*(y - 1) - x*(y^2 + 1)
            sage: p.subresultants(q, y)
            [2*x^6 - 22*x^5 + 102*x^4 - 274*x^3 + 488*x^2 - 552*x + 288,
             -x^3 - x^2*y + 6*x^2 + 5*x*y - 11*x - 6*y + 6]
            sage: p.subresultants(q, x)
            [2*y^6 - 22*y^5 + 102*y^4 - 274*y^3 + 488*y^2 - 552*y + 288,
             x*y^2 + y^3 - 5*x*y - 6*y^2 + 6*x + 11*y - 6]
        """
        R = self.parent()
        if variable is None:
            x = R.gen(0)
        else:
            x = variable
        p = self.polynomial(x)
        q = other.polynomial(x)
        return [R(f) for f in p.subresultants(q)]

    def macaulay_resultant(self, *args):
        r"""
        This is an implementation of the Macaulay resultant. It computes
        the resultant of universal polynomials as well as polynomials
        with constant coefficients. This is a project done in
        sage days 55. It's based on the implementation in Maple by
        Manfred Minimair, which in turn is based on the references [CLO], [Can], [Mac].
        It calculates the Macaulay resultant for a list of Polynomials,
        up to sign!

        AUTHORS:

        - Hao Chen, Solomon Vishkautsan (7-2014)

        INPUT:

        - ``args`` -- list of `n-1` homogeneous polynomials in `n` variables;
          works when ``args[0]`` is the list of polynomials, or ``args`` is
          itself the list of polynomials

        OUTPUT: the Macaulay resultant

        EXAMPLES:

        The number of polynomials has to match the number of variables::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: y.macaulay_resultant(x + z)                                           # needs sage.modules
            Traceback (most recent call last):
            ...
            TypeError: number of polynomials(= 2) must equal number of variables (= 3)

        The polynomials need to be all homogeneous::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: y.macaulay_resultant([x + z, z + x^3])                                # needs sage.modules
            Traceback (most recent call last):
            ...
            TypeError: resultant for non-homogeneous polynomials is not supported

        All polynomials must be in the same ring::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: S.<x,y> = PolynomialRing(QQ, 2)
            sage: y.macaulay_resultant(z + x, z)                                        # needs sage.modules
            Traceback (most recent call last):
            ...
            TypeError: not all inputs are polynomials in the calling ring

        The following example recreates Proposition 2.10 in Ch.3 of Using Algebraic Geometry::

            sage: K.<x,y> = PolynomialRing(ZZ, 2)
            sage: flist, R = K._macaulay_resultant_universal_polynomials([1,1,2])
            sage: flist[0].macaulay_resultant(flist[1:])                                # needs sage.modules
            u2^2*u4^2*u6 - 2*u1*u2*u4*u5*u6 + u1^2*u5^2*u6 - u2^2*u3*u4*u7 + u1*u2*u3*u5*u7
            + u0*u2*u4*u5*u7 - u0*u1*u5^2*u7 + u1*u2*u3*u4*u8 - u0*u2*u4^2*u8 - u1^2*u3*u5*u8
            + u0*u1*u4*u5*u8 + u2^2*u3^2*u9 - 2*u0*u2*u3*u5*u9 + u0^2*u5^2*u9
            - u1*u2*u3^2*u10 + u0*u2*u3*u4*u10 + u0*u1*u3*u5*u10 - u0^2*u4*u5*u10
            + u1^2*u3^2*u11 - 2*u0*u1*u3*u4*u11 + u0^2*u4^2*u11

        The following example degenerates into the determinant of a `3\times 3` matrix::

            sage: K.<x,y> = PolynomialRing(ZZ, 2)
            sage: flist, R = K._macaulay_resultant_universal_polynomials([1,1,1])
            sage: flist[0].macaulay_resultant(flist[1:])                                # needs sage.modules
            -u2*u4*u6 + u1*u5*u6 + u2*u3*u7 - u0*u5*u7 - u1*u3*u8 + u0*u4*u8

        The following example is by Patrick Ingram (:arxiv:`1310.4114`)::

            sage: U = PolynomialRing(ZZ,'y',2); y0,y1 = U.gens()
            sage: R = PolynomialRing(U,'x',3); x0,x1,x2 = R.gens()
            sage: f0 = y0*x2^2 - x0^2 + 2*x1*x2
            sage: f1 = y1*x2^2 - x1^2 + 2*x0*x2
            sage: f2 = x0*x1 - x2^2
            sage: f0.macaulay_resultant(f1, f2)                                         # needs sage.modules
            y0^2*y1^2 - 4*y0^3 - 4*y1^3 + 18*y0*y1 - 27

        a simple example with constant rational coefficients::

            sage: R.<x,y,z,w> = PolynomialRing(QQ, 4)
            sage: w.macaulay_resultant([z, y, x])                                       # needs sage.modules
            1

        an example where the resultant vanishes::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: (x + y).macaulay_resultant([y^2, x])                                  # needs sage.modules
            0

        an example of bad reduction at a prime ``p = 5``::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: y.macaulay_resultant([x^3 + 25*y^2*x, 5*z])                           # needs sage.libs.pari sage.modules
            125

        The input can given as an unpacked list of polynomials::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: y.macaulay_resultant(x^3 + 25*y^2*x, 5*z)                             # needs sage.libs.pari sage.modules
            125

        an example when the coefficients live in a finite field::

            sage: F = FiniteField(11)
            sage: R.<x,y,z,w> = PolynomialRing(F, 4)
            sage: z.macaulay_resultant([x^3, 5*y, w])                                   # needs sage.modules sage.rings.finite_rings
            4

        example when the denominator in the algorithm vanishes(in this case
        the resultant is the constant term of the quotient of
        char polynomials of numerator/denominator)::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: y.macaulay_resultant([x + z, z^2])                                    # needs sage.libs.pari sage.modules
            -1

        When there are only 2 polynomials, the Macaulay resultant degenerates to the traditional resultant::

            sage: R.<x> = PolynomialRing(QQ, 1)
            sage: f = x^2 + 1; g = x^5 + 1
            sage: fh = f.homogenize()
            sage: gh = g.homogenize()
            sage: RH = fh.parent()
            sage: f.resultant(g) == fh.macaulay_resultant(gh)                           # needs sage.modules
            True
        """
        if len(args) == 1 and isinstance(args[0],list):
            return self.parent().macaulay_resultant(self, *args[0])
        return self.parent().macaulay_resultant(self, *args)

    def denominator(self):
        r"""
        Return a denominator of ``self``.

        First, the lcm of the denominators of the entries of ``self``
        is computed and returned. If this computation fails, the
        unit of the parent of ``self`` is returned.

        Note that some subclasses may implement its own denominator
        function.

        .. warning::

           This is not the denominator of the rational function
           defined by ``self``, which would always be 1 since ``self`` is a
           polynomial.

        EXAMPLES:

        First we compute the denominator of a polynomial with
        integer coefficients, which is of course 1.

        ::

            sage: R.<x,y> = ZZ[]
            sage: f = x^3 + 17*y + x + y
            sage: f.denominator()
            1

        Next we compute the denominator of a polynomial over a number field.

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: R.<x,y> = NumberField(symbolic_expression(x^2+3),'a')['x,y']
            sage: f = (1/17)*x^19 + (1/6)*y - (2/3)*x + 1/3; f
            1/17*x^19 - 2/3*x + 1/6*y + 1/3
            sage: f.denominator()
            102

        Finally, we try to compute the denominator of a polynomial with
        coefficients in the real numbers, which is a ring whose elements do
        not have a denominator method.

        ::

            sage: # needs sage.rings.real_mpfr
            sage: R.<a,b,c> = RR[]
            sage: f = a + b + RR('0.3'); f
            a + b + 0.300000000000000
            sage: f.denominator()
            1.00000000000000

        Check that the denominator is an element over the base whenever the base
        has no denominator function. This closes :issue:`9063`::

            sage: R.<a,b,c> = GF(5)[]
            sage: x = R(0)
            sage: x.denominator()
            1
            sage: type(x.denominator())
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: type(a.denominator())
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: from sage.rings.polynomial.multi_polynomial_element import MPolynomial
            sage: isinstance(a / b, MPolynomial)
            False
            sage: isinstance(a.numerator() / a.denominator(), MPolynomial)
            True
        """
        if self.degree() == -1:
            return self.base_ring().one()
        x = self.coefficients()
        try:
            d = x[0].denominator()
            for y in x:
                d = d.lcm(y.denominator())
            return d
        except(AttributeError):
            return self.base_ring().one()

    def numerator(self):
        r"""
        Return a numerator of ``self``, computed as ``self * self.denominator()``.

        Note that some subclasses may implement its own numerator
        function.

        .. warning::

           This is not the numerator of the rational function
           defined by ``self``, which would always be ``self`` since ``self`` is a
           polynomial.

        EXAMPLES:

        First we compute the numerator of a polynomial with
        integer coefficients, which is of course ``self``.

        ::

            sage: R.<x, y> = ZZ[]
            sage: f = x^3 + 17*x + y + 1
            sage: f.numerator()
            x^3 + 17*x + y + 1
            sage: f == f.numerator()
            True

        Next we compute the numerator of a polynomial over a number field.

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: R.<x,y> = NumberField(symbolic_expression(x^2+3), 'a')['x,y']
            sage: f = (1/17)*y^19 - (2/3)*x + 1/3; f
            1/17*y^19 - 2/3*x + 1/3
            sage: f.numerator()
            3*y^19 - 34*x + 17
            sage: f == f.numerator()
            False

        We try to compute the numerator of a polynomial with coefficients in
        the finite field of 3 elements.

        ::

            sage: K.<x,y,z> = GF(3)['x, y, z']
            sage: f = 2*x*z + 2*z^2 + 2*y + 1; f
            -x*z - z^2 - y + 1
            sage: f.numerator()
            -x*z - z^2 - y + 1

        We check that the computation the numerator and denominator
        are valid.

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: K = NumberField(symbolic_expression('x^3+2'), 'a')['x']['s,t']
            sage: f = K.random_element()
            sage: f.numerator() / f.denominator() == f
            True
            sage: R = RR['x,y,z']
            sage: f = R.random_element()
            sage: f.numerator() / f.denominator() == f
            True
        """
        return self * self.denominator()

    def lift(self, I):
        r"""
        Given an ideal `I = (f_1,\dots,f_r)` that contains ``self``,
        find `s_1,\dots,s_r` such that ``self`` `= s_1 f_1 + ... + s_r f_r`.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: A.<x,y> = PolynomialRing(CC, 2, order='degrevlex')
            sage: I = A.ideal([x^10 + x^9*y^2, y^8 - x^2*y^7 ])
            sage: f = x*y^13 + y^12
            sage: M = f.lift(I); M                                                      # needs sage.libs.singular
            [y^7, x^7*y^2 + x^8 + x^5*y^3 + x^6*y + x^3*y^4 + x^4*y^2 + x*y^5 + x^2*y^3 + y^4]
            sage: sum(map(mul, zip(M, I.gens()))) == f                                  # needs sage.libs.singular
            True
        """
        raise NotImplementedError

    def inverse_mod(self, I):
        r"""
        Return an inverse of ``self`` modulo the polynomial ideal `I`,
        namely a multivariate polynomial `f` such that
        ``self * f - 1`` belongs to `I`.

        INPUT:

        - ``I`` -- an ideal of the polynomial ring in which ``self`` lives

        OUTPUT: a multivariate polynomial representing the inverse of ``f`` modulo `I`

        EXAMPLES::

           sage: R.<x1,x2> = QQ[]
           sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
           sage: f = x1 + 3*x2^2; g = f.inverse_mod(I); g                               # needs sage.libs.singular
           1/16*x1 + 3/16
           sage: (f*g).reduce(I)                                                        # needs sage.libs.singular
           1

        Test a non-invertible element::

           sage: R.<x1,x2> = QQ[]
           sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
           sage: f = x1 + x2
           sage: f.inverse_mod(I)                                                       # needs sage.libs.singular
           Traceback (most recent call last):
           ...
           ArithmeticError: element is non-invertible
        """
        P = self.parent()
        B  = I.gens()
        try:
            XY = P.one().lift((self,) + tuple(B))
            return P(XY[0])
        except ValueError:
            raise ArithmeticError("element is non-invertible")

    def weighted_degree(self, *weights):
        r"""
        Return the weighted degree of ``self``, which is the maximum weighted
        degree of all monomials in ``self``; the weighted degree of a monomial
        is the sum of all powers of the variables in the monomial, each power
        multiplied with its respective weight in ``weights``.

        This method is given for convenience. It is faster to use polynomial
        rings with weighted term orders and the standard ``degree`` function.

        INPUT:

        - ``weights`` -- either individual numbers, an iterable or a dictionary,
          specifying the weights of each variable. If it is a dictionary, it
          maps each variable of ``self`` to its weight. If it is a sequence of
          individual numbers or a tuple, the weights are specified in the order
          of the generators as given by ``self.parent().gens()``.

        EXAMPLES::

            sage: R.<x,y,z> = GF(7)[]
            sage: p = x^3 + y + x*z^2
            sage: p.weighted_degree({z:0, x:1, y:2})
            3
            sage: p.weighted_degree(1, 2, 0)
            3
            sage: p.weighted_degree((1, 4, 2))
            5
            sage: p.weighted_degree((1, 4, 1))
            4
            sage: p.weighted_degree(2**64, 2**50, 2**128)
            680564733841876926945195958937245974528
            sage: q = R.random_element(100, 20)
            sage: q.weighted_degree(1, 1, 1) == q.total_degree()
            True

        You may also work with negative weights

        ::

            sage: p.weighted_degree(-1, -2, -1)
            -2

        Note that only integer weights are allowed

        ::

            sage: p.weighted_degree(x, 1, 1)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to Integer Ring
            sage: p.weighted_degree(2/1, 1, 1)
            6

        The :meth:`weighted_degree` coincides with the :meth:`degree` of a weighted
        polynomial ring, but the latter is faster.

        ::

            sage: K = PolynomialRing(QQ, 'x,y', order=TermOrder('wdegrevlex', (2,3)))
            sage: p = K.random_element(10)
            sage: p.degree() == p.weighted_degree(2,3)
            True

        TESTS::

            sage: # needs sage.modules
            sage: R = PolynomialRing(QQ, 'a', 5)
            sage: f = R.random_element(terms=20)
            sage: w = random_vector(ZZ,5)
            sage: d1 = f.weighted_degree(w)
            sage: d2 = (f*1.0).weighted_degree(w)
            sage: d1 == d2
            True
        """
        if self.is_zero():
            #Corner case, note that the degree of zero is an Integer
            return Integer(-1)

        if len(weights) ==  1:
            # First unwrap it if it is given as one element argument
            weights = weights[0]

        if isinstance(weights, dict):
            weights = [weights[g] for g in self.parent().gens()]

        weights = [Integer(w) for w in weights]

        # Go through each monomial, calculating the weight
        cdef int n = self.parent().ngens()
        cdef int i, j
        cdef Integer deg
        cdef Integer l
        cdef tuple m
        A = self.exponents(as_ETuples=False)
        l = Integer(0)
        m = <tuple>(A[0])
        for i in range(n):
            l += weights[i]*m[i]
        deg = l
        for j in range(1,len(A)):
            l = Integer(0)
            m = <tuple>A[j]
            for i in range(n):
                l += weights[i]*m[i]
            if deg < l:
                deg = l
        return deg

    def gcd(self, other):
        r"""
        Return a greatest common divisor of this polynomial and ``other``.

        INPUT:

        - ``other`` -- a polynomial with the same parent as this polynomial

        EXAMPLES::

            sage: Q.<z> = Frac(QQ['z'])
            sage: R.<x,y> = Q[]
            sage: r = x*y - (2*z-1)/(z^2+z+1) * x + y/z
            sage: p = r * (x + z*y - 1/z^2)
            sage: q = r * (x*y*z + 1)
            sage: gcd(p, q)
            (z^3 + z^2 + z)*x*y + (-2*z^2 + z)*x + (z^2 + z + 1)*y

        Polynomials over polynomial rings are converted to a simpler polynomial
        ring with all variables to compute the gcd::

            sage: A.<z,t> = ZZ[]
            sage: B.<x,y> = A[]
            sage: r = x*y*z*t + 1
            sage: p = r * (x - y + z - t + 1)
            sage: q = r * (x*z - y*t)
            sage: gcd(p, q)                                                             # needs sage.libs.singular
            z*t*x*y + 1
            sage: _.parent()
            Multivariate Polynomial Ring in x, y over
             Multivariate Polynomial Ring in z, t over Integer Ring

        Some multivariate polynomial rings have no gcd implementation::

            sage: R.<x,y> = GaussianIntegers()[]                                        # needs sage.rings.number_field
            sage: x.gcd(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: GCD is not implemented for multivariate polynomials over
            Gaussian Integers generated by I in Number Field in I with defining polynomial x^2 + 1 with I = 1*I

        TESTS::

            sage: Pol = QQ['x']['x','y']
            sage: Pol.one().gcd(1)
            1

            sage: P = PolynomialRing(QQ, 'x', 0)
            sage: P.gens()
            ()
        """
        flatten = self._parent.flattening_morphism()
        tgt = flatten.codomain()
        if tgt is not self._parent and tgt._has_singular:
            g = flatten(self).gcd(flatten(other))
            return flatten.section()(g)

        try:
            self._parent._singular_().set_ring()
            g = self._singular_().gcd(other._singular_())
            return self._parent(g)
        except (TypeError, AttributeError):
            pass

        gens = self.parent().gens()
        if not gens:
            # no variables
            base = self.parent().base_ring()
            return base(self).gcd(base(other))

        x = gens[-1]
        uniself = self.polynomial(x)
        unibase = uniself.base_ring()
        try:
            doit = unibase._gcd_univariate_polynomial
        except AttributeError:
            raise NotImplementedError("GCD is not implemented for multivariate polynomials over {}".format(self._parent._mpoly_base_ring()))
        else:
            return self.parent()(doit(uniself, other.polynomial(x)))

    def nth_root(self, n):
        r"""
        Return a `n`-th root of this element.

        If there is no such root, a :exc:`ValueError` is raised.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: a = 32 * (x*y + 1)^5 * (x+y+z)^5
            sage: a.nth_root(5)
            2*x^2*y + 2*x*y^2 + 2*x*y*z + 2*x + 2*y + 2*z
            sage: b = x + 2*y + 3*z
            sage: b.nth_root(42)
            Traceback (most recent call last):
            ...
            ValueError: not a 42nd power

            sage: R.<x,y> = QQ[]
            sage: S.<z,t> = R[]
            sage: T.<u,v> = S[]
            sage: p = (1 + x*u + y + v) * (1 + z*t)
            sage: (p**3).nth_root(3)
            (x*z*t + x)*u + (z*t + 1)*v + (y + 1)*z*t + y + 1
            sage: (p**3).nth_root(3).parent() is p.parent()
            True
            sage: ((1+x+z+t)**2).nth_root(3)
            Traceback (most recent call last):
            ...
            ValueError: not a 3rd power
        """
        R = self.parent()
        phi = R.flattening_morphism()
        S = phi.codomain()
        p = phi(self)

        V = p.variables()
        if not V:
            # constant
            root = self.constant_coefficient().nth_root(n)
            return phi.section()(S(root))
        elif len(V) == 1:
            # univariate
            U = PolynomialRing(S.base_ring(), str(V[0]))
            pU = U(p)
        else:
            # specialize one variable
            # (in order to call the univariate case)
            U0 = PolynomialRing(S.base_ring(), [str(v) for v in V[:-1]])
            U = U0[str(V[-1])]
            pU = U(p)

        # recursive call
        root = pU.nth_root(n)
        return phi.section()(S(root))

    def is_square(self, root=False):
        r"""
        Test whether this polynomial is a square.

        INPUT:

        - ``root`` -- if set to ``True``, return a pair ``(True, root)``
          where ``root`` is a square root or ``(False, None)`` if
          it is not a square.

        EXAMPLES::

            sage: R.<a,b> = QQ[]
            sage: a.is_square()
            False
            sage: ((1+a*b^2)^2).is_square()
            True
            sage: ((1+a*b^2)^2).is_square(root=True)
            (True, a*b^2 + 1)
        """
        try:
            sqrt = self.nth_root(2)
        except ValueError:
            return (False,None) if root else False
        else:
            return (True,sqrt) if root else True

    def specialization(self, D=None, phi=None):
        r"""
        Specialization of this polynomial.

        Given a family of polynomials defined over a polynomial ring. A specialization
        is a particular member of that family. The specialization can be specified either
        by a dictionary or a :class:`SpecializationMorphism`.

        INPUT:

        - ``D`` -- dictionary (optional)

        - ``phi`` -- :class:`SpecializationMorphism` (optional)

        OUTPUT: a new polynomial

        EXAMPLES::

            sage: R.<c> = PolynomialRing(QQ)
            sage: S.<x,y> = PolynomialRing(R)
            sage: F = x^2 + c*y^2
            sage: F.specialization({c:2})
            x^2 + 2*y^2

        ::

            sage: S.<a,b> = PolynomialRing(QQ)
            sage: P.<x,y,z> = PolynomialRing(S)
            sage: RR.<c,d> = PolynomialRing(P)
            sage: f = a*x^2 + b*y^3 + c*y^2 - b*a*d + d^2 - a*c*b*z^2
            sage: f.specialization({a:2, z:4, d:2})
            (y^2 - 32*b)*c + b*y^3 + 2*x^2 - 4*b + 4

        Check that we preserve multi- versus uni-variate::

            sage: R.<l> = PolynomialRing(QQ, 1)
            sage: S.<k> = PolynomialRing(R)
            sage: K.<a, b, c> = PolynomialRing(S)
            sage: F = a*k^2 + b*l + c^2
            sage: F.specialization({b:56, c:5}).parent()
            Univariate Polynomial Ring in a over Univariate Polynomial Ring in k
            over Multivariate Polynomial Ring in l over Rational Field
        """
        if D is None:
            if phi is None:
                raise ValueError("either the dictionary or the specialization must be provided")
        else:
            from sage.rings.polynomial.flatten import SpecializationMorphism
            phi = SpecializationMorphism(self.parent(),D)
        return phi(self)

    def reduced_form(self, **kwds):
        r"""
        Return a reduced form of this polynomial.

        The algorithm is from Stoll and Cremona's "On the Reduction Theory of
        Binary Forms" [CS2003]_. This takes a two variable homogeneous polynomial and
        finds a reduced form. This is a `SL(2,\ZZ)`-equivalent binary form
        whose covariant in the upper half plane is in the fundamental domain.
        If the polynomial has multiple roots, they are removed and the algorithm
        is applied to the portion without multiple roots.

        This reduction should also minimize the sum of the squares of the coefficients,
        but this is not always the case.  By default the coefficient minimizing
        algorithm in [HS2018]_ is applied. The coefficients can be minimized
        either with respect to the sum of their squares or the maximum of their
        global heights.

        A portion of the algorithm uses Newton's method to find a solution to
        a system of equations. If Newton's method fails to converge to a point
        in the upper half plane, the function will use the less precise `z_0`
        covariant from the `Q_0` form as defined on page 7 of [CS2003]_.
        Additionally, if this polynomial has
        a root with multiplicity at least half the total degree of the polynomial,
        then we must also use the `z_0` covariant. See [CS2003]_ for details.

        Note that, if the covariant is within ``error_limit`` of the boundary
        but outside the fundamental domain, our function will erroneously move
        it to within the fundamental domain, hence our conjugation will be off
        by 1. If you don't want this to happen, decrease your ``error_limit``
        and increase your precision.

        Implemented by Rebecca Lauren Miller as part of GSOC 2016. Smallest
        coefficients added by Ben Hutz July 2018.

        INPUT: keyword arguments:

        - ``prec`` -- integer (default: 300); sets the precision

        - ``return_conjugation`` -- boolean (default: ``True``); whether to
          return element of `SL(2, \ZZ)`

        - ``error_limit`` -- sets the error tolerance (default: 0.000001)

        - ``smallest_coeffs`` -- boolean (default: ``True``); whether to find the
          model with smallest coefficients

        - ``norm_type`` -- either ``'norm'`` or ``'height'``; what type of norm
          to use for smallest coefficients

        - ``emb`` -- (optional) embedding of based field into ``CC``

        OUTPUT:

        - a polynomial (reduced binary form)

        - a matrix (element of `SL(2, \ZZ)`)

        .. TODO::

            When Newton's Method doesn't converge to a root in the upper half plane.
            Now we just return `z_0`. It would be better to modify and find the unique root
            in the upper half plane.

        EXAMPLES::

            sage: R.<x,h> = PolynomialRing(QQ)
            sage: f = 19*x^8 - 262*x^7*h + 1507*x^6*h^2 - 4784*x^5*h^3 + 9202*x^4*h^4\
             -10962*x^3*h^5 + 7844*x^2*h^6 - 3040*x*h^7 + 475*h^8
            sage: f.reduced_form(prec=200, smallest_coeffs=False)                       # needs sage.modules sage.rings.complex_interval_field
            (
            -x^8 - 2*x^7*h + 7*x^6*h^2 + 16*x^5*h^3 + 2*x^4*h^4 - 2*x^3*h^5 + 4*x^2*h^6 - 5*h^8,
            <BLANKLINE>
            [ 1 -2]
            [ 1 -1]
            )

        An example where the multiplicity is too high::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = x^3 + 378666*x^2*y - 12444444*x*y^2 + 1234567890*y^3
            sage: j = f * (x-545*y)^9
            sage: j.reduced_form(prec=200, smallest_coeffs=False)                       # needs sage.modules sage.rings.complex_interval_field
            Traceback (most recent call last):
            ...
            ValueError: cannot have a root with multiplicity >= 12/2

        An example where Newton's Method does not find the right root::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: F = x^6 + 3*x^5*y - 8*x^4*y^2 - 2*x^3*y^3 - 44*x^2*y^4 - 8*x*y^5
            sage: F.reduced_form(smallest_coeffs=False, prec=400)                       # needs sage.modules sage.rings.complex_interval_field
            Traceback (most recent call last):
            ...
            ArithmeticError: Newton's method converged to z not in the upper half plane

        An example with covariant on the boundary, therefore a non-unique form::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: F = 5*x^2*y - 5*x*y^2 - 30*y^3
            sage: F.reduced_form(smallest_coeffs=False)                                 # needs sage.modules sage.rings.complex_interval_field
            (
                                        [1 1]
            5*x^2*y + 5*x*y^2 - 30*y^3, [0 1]
            )

        An example where precision needs to be increased::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: F = (-16*x^7 - 114*x^6*y - 345*x^5*y^2 - 599*x^4*y^3
            ....:      - 666*x^3*y^4 - 481*x^2*y^5 - 207*x*y^6 - 40*y^7)
            sage: F.reduced_form(prec=50, smallest_coeffs=False)                        # needs sage.modules sage.rings.complex_interval_field
            Traceback (most recent call last):
            ...
            ValueError: accuracy of Newton's root not within tolerance(0.000012... > 1e-06),
            increase precision
            sage: F.reduced_form(prec=100, smallest_coeffs=False)                       # needs sage.modules sage.rings.complex_interval_field
            (
                                                                  [-1 -1]
            -x^5*y^2 - 24*x^3*y^4 - 3*x^2*y^5 - 2*x*y^6 + 16*y^7, [ 1  0]
            )

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: F = - 8*x^4 - 3933*x^3*y - 725085*x^2*y^2 - 59411592*x*y^3 - 1825511633*y^4
            sage: F.reduced_form(return_conjugation=False)                              # needs sage.modules sage.rings.complex_interval_field
            x^4 + 9*x^3*y - 3*x*y^3 - 8*y^4

        ::

            sage: R.<x,y> = QQ[]
            sage: F = -2*x^3 + 2*x^2*y + 3*x*y^2 + 127*y^3
            sage: F.reduced_form()                                                      # needs sage.modules sage.rings.complex_interval_field
            (
                                                   [1 4]
            -2*x^3 - 22*x^2*y - 77*x*y^2 + 43*y^3, [0 1]
            )

        ::

            sage: R.<x,y> = QQ[]
            sage: F = -2*x^3 + 2*x^2*y + 3*x*y^2 + 127*y^3
            sage: F.reduced_form(norm_type='height')                                    # needs sage.modules sage.rings.complex_interval_field
            (
                                                    [5 4]
            -58*x^3 - 47*x^2*y + 52*x*y^2 + 43*y^3, [1 1]
            )

        ::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: F = x^4 + x^3*y*z + y^2*z
            sage: F.reduced_form()                                                      # needs sage.modules sage.rings.complex_interval_field
            Traceback (most recent call last):
            ...
            ValueError: (=x^3*y*z + x^4 + y^2*z) must have two variables

        ::

            sage: R.<x,y> = PolynomialRing(ZZ)
            sage: F = - 8*x^6 - 3933*x^3*y - 725085*x^2*y^2 - 59411592*x*y^3 - 99*y^6
            sage: F.reduced_form(return_conjugation=False)                              # needs sage.modules sage.rings.complex_interval_field
            Traceback (most recent call last):
            ...
            ValueError: (=-8*x^6 - 99*y^6 - 3933*x^3*y - 725085*x^2*y^2 -
            59411592*x*y^3) must be homogeneous

        ::

            sage: R.<x,y> = PolynomialRing(RR)
            sage: F = (217.992172373276*x^3 + 96023.1505442490*x^2*y
            ....:      + 1.40987971253579e7*x*y^2 + 6.90016027113216e8*y^3)
            sage: F.reduced_form(smallest_coeffs=False)  # tol 1e-8                     # needs sage.modules sage.rings.complex_interval_field
            (
            -39.5673942565918*x^3 + 111.874026298523*x^2*y
             + 231.052762985229*x*y^2 - 138.380829811096*y^3,
            <BLANKLINE>
            [-147 -148]
            [   1    1]
            )

        ::

            sage: R.<x,y> = PolynomialRing(CC)                                          # needs sage.rings.real_mpfr
            sage: F = ((0.759099196558145 + 0.845425869641446*CC.0)*x^3                 # needs sage.rings.real_mpfr
            ....:      + (84.8317207268542 + 93.8840848648033*CC.0)*x^2*y
            ....:      + (3159.07040755858 + 3475.33037377779*CC.0)*x*y^2
            ....:      + (39202.5965389079 + 42882.5139724962*CC.0)*y^3)
            sage: F.reduced_form(smallest_coeffs=False)  # tol 1e-11                    # needs sage.modules sage.rings.complex_interval_field sage.rings.real_mpfr
            (
            (-0.759099196558145 - 0.845425869641446*I)*x^3
            + (-0.571709908900118 - 0.0418133346027929*I)*x^2*y
            + (0.856525964330103 - 0.0721403997649759*I)*x*y^2
            + (-0.965531044130330 + 0.754252314465703*I)*y^3,
            <BLANKLINE>
            [-1 37]
            [ 0 -1]
            )
        """
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector
        from sage.rings.complex_interval_field import ComplexIntervalField
        from sage.rings.real_mpfr import RealField

        if self.parent().ngens() != 2:
            raise ValueError("(=%s) must have two variables"%self)
        if not self.is_homogeneous():
            raise ValueError("(=%s) must be homogeneous"%self)

        prec = kwds.get('prec', 300)
        return_conjugation  =kwds.get('return_conjugation', True)
        error_limit = kwds.get('error_limit', 0.000001)
        emb = kwds.get('emb', None)

        # getting a numerical approximation of the roots of our polynomial
        CF = ComplexIntervalField(prec=prec)  # keeps trac of our precision error
        RF = RealField(prec=prec)
        R = self.parent()
        x,y = R.gens()

        # finding quadratic Q_0, gives us our covariant, z_0
        from sage.rings.polynomial.binary_form_reduce import covariant_z0
        try:
            z, th = covariant_z0(self, prec=prec, emb=emb, z0_cov=True)
        except ValueError:# multiple roots
            F = self.lc()*prod([p for p,e in self.factor()])
            z, th = covariant_z0(F, prec=prec, emb=emb, z0_cov=True)
        z = CF(z)
        # this moves z_0 to our fundamental domain using the three steps laid
        # out in the algorithm by [CS2003]
        # this is found in section 5 of their paper
        M = matrix(QQ, [[1,0], [0,1]]) # used to keep track of how our z is moved.
        zc = z.center()
        while zc.real() < RF(-0.5) or zc.real() >= RF(0.5) or (zc.real() <= RF(0) and zc.abs() < RF(1))\
         or (zc.real() > RF(0) and zc.abs() <= RF(1)):
            if (zc.real() < RF(-0.5)) or (zc.real() >= RF(0.5)):
                # moves z into fundamental domain by m
                m = zc.real().round() # finds amount to move z's real part by
                Qm = QQ(m)
                M = M * matrix(QQ, [[1,Qm], [0,1]]) # move
                z -= m  # M.inverse()*z is supposed to move z by m
            elif (zc.real() <= RF(0) and zc.abs() < RF(1)) or (zc.real() > RF(0) and zc.abs() <= RF(1)): # flips z
                z = -1/z
                M = M * matrix(QQ, [[0,-1], [1,0]])# multiply on left because we are taking inverse matrices
            zc = z.center()

        smallest_coeffs = kwds.get('smallest_coeffs', True)
        if smallest_coeffs:
            # since we are searching anyway, don't need the 'true' reduced covariant
            from sage.rings.polynomial.binary_form_reduce import smallest_poly
            norm_type = kwds.get('norm_type', 'norm')
            sm_F, sm_m = smallest_poly(self(tuple(M * vector([x,y]))), prec=prec, norm_type=norm_type, emb=emb)
            M = M*sm_m
        else:
            # solve the minimization problem for 'true' covariant
            z, th = covariant_z0(self(tuple(M * vector([x,y]))), prec=prec, emb=emb)
            z = CF(z)
            zc = z.center()
            # moves our z to fundamental domain as before
            while zc.real() < RF(-0.5) or zc.real() >= RF(0.5) or (zc.real() <= RF(0) and zc.abs() < RF(1))\
             or (zc.real() > RF(0) and zc.abs() <= RF(1)):
                if (zc.real() < RF(-0.5)) or (zc.real() >= RF(0.5)):
                    # moves z into fundamental domain by m
                    m = zc.real().round() # finds amount to move z's real part by
                    Qm = QQ(m)
                    M = M * matrix(QQ, [[1,Qm], [0,1]]) # move
                    z -= m  # M.inverse()*z is supposed to move z by m
                elif (zc.real() <= RF(0) and zc.abs() < RF(1)) or (zc.real() > RF(0) and zc.abs() <= RF(1)): # flips z
                    z = -1/z
                    M = M * matrix(QQ, [[0,-1], [1,0]])# multiply on left because we are taking inverse matrices
                zc = z.center()

        if return_conjugation:
            return (self(tuple(M * vector([x,y]))), M)
        return self(tuple(M * vector([x,y])))

    def is_unit(self):
        r"""
        Return ``True`` if ``self`` is a unit, that is, has a
        multiplicative inverse.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x,y> = QQbar[]
            sage: (x + y).is_unit()
            False
            sage: R(0).is_unit()
            False
            sage: R(-1).is_unit()
            True
            sage: R(-1 + x).is_unit()
            False
            sage: R(2).is_unit()
            True

        Check that :issue:`22454` is fixed::

            sage: _.<x,y> = Zmod(4)[]
            sage: (1 + 2*x).is_unit()
            True
            sage: (x*y).is_unit()
            False
            sage: _.<x,y> = Zmod(36)[]
            sage: (7+ 6*x + 12*y - 18*x*y).is_unit()
            True
        """
        # EXERCISE (Atiyah-McDonald, Ch 1): Let `A[x]` be a polynomial
        # ring in one variable. Then `f=\sum a_i x^i \in A[x]` is a unit\
        # if and only if `a_0` is a unit and `a_1,\ldots, a_n` are nilpotent.
        # (Also noted in Dummit and Foote, "Abstract Algebra", 1991,
        # Section 7.3 Exercise 33).
        # Also f is nilpotent if and only if all a_i are nilpotent.
        # This generalizes easily to the multivariate case, by considering
        # K[x,y,...] as K[x][y]...
        if not self.constant_coefficient().is_unit():
            return False
        cdef dict d = self.dict()
        cdef ETuple zero_key = ETuple({}, int(self.parent().ngens()))
        d.pop(zero_key, None)
        return all(d[k].is_nilpotent() for k in d)

    def is_nilpotent(self):
        r"""
        Return ``True`` if ``self`` is nilpotent, i.e., some power of ``self``
        is 0.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]                                                     # needs sage.rings.number_field
            sage: (x + y).is_nilpotent()                                                # needs sage.rings.number_field
            False
            sage: R(0).is_nilpotent()                                                   # needs sage.rings.number_field
            True
            sage: _.<x,y> = Zmod(4)[]
            sage: (2*x).is_nilpotent()
            True
            sage: (2 + y*x).is_nilpotent()
            False
            sage: _.<x,y> = Zmod(36)[]
            sage: (4 + 6*x).is_nilpotent()
            False
            sage: (6*x + 12*y + 18*x*y + 24*(x^2+y^2)).is_nilpotent()
            True
        """
        # EXERCISE (Atiyah-McDonald, Ch 1): Let `A[x]` be a polynomial
        # ring in one variable. Then `f=\sum a_i x^i \in A[x]` is
        # nilpotent if and only if `a_0,\ldots, a_n` are nilpotent.
        # (Also noted in Dummit and Foote, "Abstract Algebra", 1991,
        # Section 7.3 Exercise 33).
        # This generalizes easily to the multivariate case, by considering
        # K[x,y,...] as K[x][y]...
        d = self.dict()
        return all(c.is_nilpotent() for c in d.values())

    def _test_subs(self, tester=None, **options):
        r"""
        Run some tests using the ``subs`` method.

        TESTS::

            sage: R.<x,y> = QQbar[]                                                     # needs sage.rings.number_field
            sage: (x + y)._test_subs()                                                  # needs sage.rings.number_field
        """
        if tester is None:
            tester = self._tester(**options)

        gens = self.parent().gens()

        if gens:
            # substituting all variables (in a polynomial ring with variables) with 0
            d = {str(gen): self.base_ring().zero() for gen in gens}
            tester.assertEqual(self.subs(**d).parent(), self.parent().base_ring())

            # substituting all variables (in a polynomial ring with variables)
            # with elements of another ring
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            other = PolynomialRing(self.parent().base_ring(), 'other', len(gens))
            other_gens = other.gens()
            d = {str(gen): other_gen for gen, other_gen in zip(gens, other_gens)}
            tester.assertEqual(self.subs(**d).parent(), other)

        if len(gens) > 1:
            # substituting one variable (in a polynomial ring with variables) with 0
            d = {str(gens[0]): self.base_ring().zero()}
            tester.assertEqual(self.subs(**d).parent(), self.parent())

            # test error checking: partial substitution by elements
            # from another ring is not allowed
            d = {str(gens[0]): other_gens[0]}
            with tester.assertRaises((ValueError, TypeError)):
                self.subs(**d)

    def is_lorentzian(self, explain=False):
        r"""
        Return whether this is a Lorentzian polynomial.

        INPUT:

        - ``explain`` -- boolean (default: ``False``); if ``True``
          return a tuple whose first element is the boolean result of the test,
          and the second element is a string describing the reason the test failed,
          or ``None`` if the test succeeded.

        Lorentzian polynomials are a class of polynomials connected with the area
        of discrete convex analysis.  A polynomial `f` with positive real coefficients
        is Lorentzian if:

        - `f` is homogeneous;

        - the support of `f` is `M`-convex

        - `f` has degree less than `2`, or if its degree is at least two,
          the collection of sequential partial derivatives of `f` which are
          quadratic forms have Gram matrices with at most one positive eigenvalue.

        Note in particular that the zero polynomial is Lorentzian.  Examples of
        Lorentzian polynomials include homogeneous stable polynomials, volume
        polynomials of convex bodies and projective varieties, and Schur polynomials
        after renormalizing the coefficient of each monomial `x^\alpha` by `1/\alpha!`.

        EXAMPLES:

        Renormalized Schur polynomials are Lorentzian, but not in general if the
        renormalization is skipped::

            sage: P.<x,y> = QQ[]
            sage: p = (x^2 / 2) + x*y + (y^2 / 2)
            sage: p.is_lorentzian()
            True
            sage: p = x^2 + x*y + y^2
            sage: p.is_lorentzian()
            False

        Homogeneous linear forms and constant polynomials with positive
        coefficients are Lorentzian, as well as the zero polynomial::

            sage: p = x + 2*y
            sage: p.is_lorentzian()
            True
            sage: p = P(5)
            sage: p.is_lorentzian()
            True
            sage: P.zero().is_lorentzian()
            True

        Inhomogeneous polynomials and polynomials with negative coefficients
        are not Lorentzian::

            sage: p = x^2 + 2*x + y^2
            sage: p.is_lorentzian()
            False
            sage: p = 2*x^2 - y^2
            sage: p.is_lorentzian()
            False

        It is an error to check if a polynomial is Lorentzian if its base ring
        is not a subring of the real numbers, as the notion is not defined in
        this case::

            sage: # needs sage.rings.real_mpfr
            sage: Q.<z,w> = CC[]
            sage: q = z^2 + w^2
            sage: q.is_lorentzian()
            Traceback (most recent call last):
            ...
            NotImplementedError: is_lorentzian only implemented for real polynomials

        The method can give a reason for a polynomial failing to be Lorentzian::

            sage: p = x^2 + 2*x + y^2
            sage: p.is_lorentzian(explain=True)
            (False, 'inhomogeneous')

        REFERENCES:

        For full definitions and related discussion, see [BrHu2019]_ and
        [HMMS2019]_.  The second reference gives the characterization of
        Lorentzian polynomials applied in this implementation explicitly.
        """
        from sage.rings.imaginary_unit import I

        # function to handle return value when reason requested
        def result(val, explanation=None):
            return (val, explanation) if explain else val

        try:
            # this would better be handled by a category of RealFields()
            self.base_ring()(I)
        except (ValueError, TypeError):
            pass
        else:
            raise NotImplementedError("is_lorentzian only implemented for real polynomials")

        if self.is_zero():
            return result(True)

        if not self.is_homogeneous():
            return result(False, "inhomogeneous")

        if any(coeff < 0 for coeff in self.coefficients()):
            return result(False, "negative coefficient")

        # for degree <= 1, homogeneous with positive coefficients is sufficient
        if self.degree() <= 1:
            return result(True)

        # check support is M-convex
        if not _is_M_convex_(self.exponents()):
            return result(False, "not M-convex")

        # compute quadratic forms coming from a sequence of partial derivatives
        if self.degree() == 2:
            quadratic_derivs = set([self])
        else:
            from sage.combinat.integer_lists.invlex import IntegerListsLex

            gens = self.parent().gens()
            quadratic_derivs = set()
            multi_exponents = IntegerListsLex(self.degree() - 2, length=len(gens))
            for alpha in multi_exponents:
                # construct list [gen_1, exp_1, ..., gen_n, exp_n] for derivative function
                d_list = chain(*zip(gens, alpha))
                d = self.derivative(*d_list)
                quadratic_derivs.add(d)

        # check derivative quadratic forms have at most one positive eigenvalue
        for deriv in quadratic_derivs:
            from sage.quadratic_forms.quadratic_form import QuadraticForm
            G = QuadraticForm(deriv).Gram_matrix()
            spectrum = sorted(G.eigenvalues(), reverse=True)
            if len(spectrum) > 1 and spectrum[1] > 0:
                return result(False, "multiple positive eigenvalues")

        return result(True)


def _is_M_convex_(points):
    r"""
    Return whether ``points`` represents a set of integer lattice points
    which are M-convex.

    Utility function for method ``is_lorentzian``, which would more properly
    fit with code related to discrete convex geometry, generalized permutahedra,
    or polymatroids, which are not currently implemented in Sage.

    INPUT:

    - ``points`` -- iterable for a list of integer lattice points of the
      same dimension

    Examples of M-convex sets include the vertices of a matroid polytope, and the
    support sets of Schur polynomials.

    EXAMPLES:

    The following points represent the vertices of a matroid polytope (indicator
    vectors of the bases) of rank `2` on five elements::

        sage: from sage.rings.polynomial.multi_polynomial import _is_M_convex_
        sage: P = [[1,1,0,0], [1,0,1,0], [0,1,1,0], [0,1,0,1], [0,0,1,1]]
        sage: _is_M_convex_(P)
        True

    These points are the support of the Schur polynomial in three variables for
    the partition `(2,2)`::

        sage: P = [[2,2,0], [2,0,2], [0,2,2], [2,1,1], [1,2,1], [1,1,2]]
        sage: _is_M_convex_(P)
        True

    The following are not examples of `M`-convex sets of points::

        sage: P = [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
        sage: _is_M_convex_(P)
        False

        sage: P = [[0, 1, 2], [2, 1]]
        sage: _is_M_convex_(P)
        Traceback (most recent call last):
        ...
        ValueError: input points are not the same dimension

        sage: P = [[0, 0.5, 1], [1, 1.5, 2]]
        sage: _is_M_convex_(P)
        Traceback (most recent call last):
        ...
        ValueError: input points are not integer lattice points

    REFERENCES:

    See [BrHu2019]_ for a definition of M-convexity.
    """
    points_set = set(map(tuple, points))
    if not points_set:
        return True
    dim = len(next(iter(points_set)))
    if any(len(p) != dim for p in points_set):
        raise ValueError("input points are not the same dimension")
    if any(entry not in ZZ for p in points_set for entry in p):
        raise ValueError("input points are not integer lattice points")
    for p1 in points_set:
        list_p1 = list(p1)
        for p2 in points_set:
            if p2 == p1:
                continue
            delta = list(x2 - x1 for x1, x2 in zip(p1, p2))
            for i in range(dim):
                if p2[i] > p1[i]:
                    # modify list_p1 to represent point p1 + e_i - e_j for various i, j
                    list_p1[i] += 1  # add e_i
                    # check exchange condition is satisfied by some index j
                    for j in range(dim):
                        if p2[j] < p1[j]:
                            list_p1[j] -= 1  # subtract e_j
                            exch = tuple(list_p1)  # p1 + e_i - e_j
                            list_p1[j] += 1  # add e_j again
                            if tuple(exch) in points_set:
                                break
                    else:
                        return False
                    list_p1[i] -= 1  # subtract e_i
                    # list_p1 should now have same entries as p1 again
    return True


cdef remove_from_tuple(e, int ind):
    w = list(e)
    del w[ind]
    if len(w) == 1:
        return w[0]
    else:
        return tuple(w)


cdef class MPolynomial_libsingular(MPolynomial):
    r"""
    Abstract base class for :class:`~sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular`.

    This class is defined for the purpose of :func:`isinstance` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: from sage.rings.polynomial.multi_polynomial import MPolynomial_libsingular
        sage: R1.<x> = QQ[]
        sage: isinstance(x, MPolynomial_libsingular)
        False
        sage: R2.<y,z> = QQ[]
        sage: isinstance(y, MPolynomial_libsingular)                                    # needs sage.libs.singular
        True

    By design, there is a unique direct subclass::

        sage: len(sage.rings.polynomial.multi_polynomial.MPolynomial_libsingular.__subclasses__()) <= 1
        True
    """
