r"""
Base class for elements of multivariate polynomial rings
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.structure.element cimport coercion_model
from sage.misc.derivative import multi_derivative
from sage.rings.infinity import infinity
from sage.structure.element cimport Element

from sage.misc.all import prod

def is_MPolynomial(x):
    return isinstance(x, MPolynomial)

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.categories.map cimport Map
from sage.categories.morphism cimport Morphism
from sage.modules.free_module_element import vector
from sage.rings.rational_field import QQ
from sage.arith.misc import gcd
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.real_mpfr import RealField_class,RealField

cdef class MPolynomial(CommutativeRingElement):

    ####################
    # Some standard conversions
    ####################
    def __int__(self):
        """
        TESTS::

            sage: type(RR['x,y'])
            <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_domain_with_category'>
            sage: type(RR['x, y'](0))
            <class 'sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict'>

            sage: int(RR['x,y'](0)) # indirect doctest
            0
            sage: int(RR['x,y'](10))
            10
            sage: int(RR['x,y'].gen())
            Traceback (most recent call last):
            ...
            TypeError...
        """
        if self.degree() <= 0:
            return int(self.constant_coefficient())
        else:
            raise TypeError

    def __long__(self):
        """
        TESTS::

            sage: long(RR['x,y'](0)) # indirect doctest
            0L
        """
        if self.degree() <= 0:
            return long(self.constant_coefficient())
        else:
            raise TypeError

    def __float__(self):
        """
        TESTS::

            sage: float(RR['x,y'](0)) # indirect doctest
            0.0
        """
        if self.degree() <= 0:
            return float(self.constant_coefficient())
        else:
            raise TypeError

    def _mpfr_(self, R):
        """
        TESTS::

            sage: RR(RR['x,y'](0)) # indirect doctest
            0.000000000000000
        """
        if self.degree() <= 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _complex_mpfr_field_(self, R):
        """
        TESTS::

            sage: CC(RR['x,y'](0)) # indirect doctest
            0.000000000000000
        """
        if self.degree() <= 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _complex_double_(self, R):
        """
        TESTS::

            sage: CDF(RR['x,y'](0)) # indirect doctest
            0.0
        """
        if self.degree() <= 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _real_double_(self, R):
        """
        TESTS::

            sage: RR(RR['x,y'](0)) # indirect doctest
            0.000000000000000
        """
        if self.degree() <= 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _rational_(self):
        """
        TESTS::

            sage: QQ(RR['x,y'](0)) # indirect doctest
            0
            sage: QQ(RR['x,y'](0.5)) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError...
        """
        if self.degree() <= 0:
            from sage.rings.rational import Rational
            return Rational(repr(self))
        else:
            raise TypeError

    def _integer_(self, ZZ=None):
        """
        TESTS::

            sage: ZZ(RR['x,y'](0)) # indirect doctest
            0
            sage: ZZ(RR['x,y'](0.0))
            0
            sage: ZZ(RR['x,y'](0.5))
            Traceback (most recent call last):
            ...
            TypeError...
        """
        if self.degree() <= 0:
            from sage.rings.integer import Integer
            return Integer(repr(self))
        else:
            raise TypeError

    def _symbolic_(self, R):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = x^3 + y
            sage: g = f._symbolic_(SR); g
            x^3 + y
            sage: g(x=2,y=2)
            10

            sage: g = SR(f)
            sage: g(x=2,y=2)
            10
        """
        d = dict([(repr(g), R.var(g)) for g in self.parent().gens()])
        return self.subs(**d)

    def _polynomial_(self, R):
        var = R.variable_name()
        if var in self._parent.variable_names():
            return R(self.polynomial(self._parent(var)))
        else:
            return R([self])

    def coefficients(self):
        """
        Return the nonzero coefficients of this polynomial in a list.
        The returned list is decreasingly ordered by the term ordering
        of ``self.parent()``, i.e. the list of coefficients matches the list
        of monomials returned by
        :meth:`sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular.monomials`.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3,order='degrevlex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]
            sage: R.<x,y,z> = PolynomialRing(QQ,3,order='lex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

        Test the same stuff with base ring `\ZZ` -- different implementation::

            sage: R.<x,y,z> = PolynomialRing(ZZ,3,order='degrevlex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]
            sage: R.<x,y,z> = PolynomialRing(ZZ,3,order='lex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

        AUTHOR:

        - Didier Deshommes
        """
        degs = self.exponents()
        d = self.dict()
        return  [ d[i] for i in degs ]

    def truncate(self, var, n):
        """
        Returns a new multivariate polynomial obtained from self by
        deleting all terms that involve the given variable to a power
        at least n.
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

    def _fast_float_(self, *vars):
        """
        Returns a quickly-evaluating function on floats.

        EXAMPLES::

            sage: K.<x,y,z> = QQ[]
            sage: f = (x+2*y+3*z^2)^2 + 42
            sage: f(1, 10, 100)
            901260483
            sage: ff = f._fast_float_()
            sage: ff(0, 0, 1)
            51.0
            sage: ff(0, 1, 0)
            46.0
            sage: ff(1, 10, 100)
            901260483.0
            sage: ff_swapped = f._fast_float_('z', 'y', 'x')
            sage: ff_swapped(100, 10, 1)
            901260483.0
            sage: ff_extra = f._fast_float_('x', 'A', 'y', 'B', 'z', 'C')
            sage: ff_extra(1, 7, 10, 13, 100, 19)
            901260483.0

        Currently, we use a fairly unoptimized method that evaluates one
        monomial at a time, with no sharing of repeated computations and
        with useless additions of 0 and multiplications by 1::

            sage: list(ff)
            ['push 0.0', 'push 12.0', 'load 1', 'load 2', 'dup', 'mul', 'mul', 'mul', 'add', 'push 4.0', 'load 0', 'load 1', 'mul', 'mul', 'add', 'push 42.0', 'add', 'push 1.0', 'load 0', 'dup', 'mul', 'mul', 'add', 'push 9.0', 'load 2', 'dup', 'mul', 'dup', 'mul', 'mul', 'add', 'push 6.0', 'load 0', 'load 2', 'dup', 'mul', 'mul', 'mul', 'add', 'push 4.0', 'load 1', 'dup', 'mul', 'mul', 'add']

        TESTS::

            sage: from sage.ext.fast_eval import fast_float
            sage: list(fast_float(K(0), old=True))
            ['push 0.0']
            sage: list(fast_float(K(17), old=True))
            ['push 0.0', 'push 17.0', 'add']
            sage: list(fast_float(y, old=True))
            ['push 0.0', 'push 1.0', 'load 1', 'mul', 'add']
        """
        from sage.ext.fast_eval import fast_float_arg, fast_float_constant
        my_vars = self.parent().variable_names()
        vars = list(vars)
        if len(vars) == 0:
            indices = range(len(my_vars))
        else:
            indices = [vars.index(v) for v in my_vars]
        x = [fast_float_arg(i) for i in indices]

        n = len(x)
        expr = fast_float_constant(0)
        for (m,c) in self.dict().iteritems():
            monom = prod([ x[i]**m[i] for i in range(n) if m[i] != 0], fast_float_constant(c))
            expr = expr + monom
        return expr

    def _fast_callable_(self, etb):
        """
        Given an ExpressionTreeBuilder, return an Expression representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y','z'])
            sage: K.<x,y,z> = QQ[]
            sage: v = K.random_element(degree=3, terms=4); v
            -6/5*x*y*z + 2*y*z^2 - x
            sage: v._fast_callable_(etb)
            add(add(add(0, mul(-6/5, mul(mul(ipow(v_0, 1), ipow(v_1, 1)), ipow(v_2, 1)))), mul(2, mul(ipow(v_1, 1), ipow(v_2, 2)))), mul(-1, ipow(v_0, 1)))

        TESTS::

            sage: v = K(0)
            sage: vf = fast_callable(v)
            sage: type(v(0r, 0r, 0r))
            <type 'sage.rings.rational.Rational'>
            sage: type(vf(0r, 0r, 0r))
            <type 'sage.rings.rational.Rational'>
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

        expr = etb.constant(self.base_ring()(0))
        for (m, c) in self.dict().iteritems():
            monom = prod([ x[i]**m[i] for i in range(n) if m[i] != 0],
                             etb.constant(c))
            expr = expr + monom
        return expr

    def derivative(self, *args):
        r"""
        The formal derivative of this polynomial, with respect to
        variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more details.

        .. SEEALSO:: :meth:`._derivative`

        EXAMPLES:

        Polynomials implemented via Singular::

            sage: R.<x, y> = PolynomialRing(FiniteField(5))
            sage: f = x^3*y^5 + x^7*y
            sage: type(f)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
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
            sage: f.derivative(x)   # with respect to x
            (2*t^2 + O(t^3))*x*y^3 + (111*t^4 + O(t^5))*x^2
            sage: f.derivative(y)   # with respect to y
            (3*t^2 + O(t^3))*x^2*y^2
            sage: f.derivative(t)   # with respect to t (recurses into base ring)
            (2*t + O(t^2))*x^2*y^3 + (148*t^3 + O(t^4))*x^3
            sage: f.derivative(x, y) # with respect to x and then y
            (6*t^2 + O(t^3))*x*y^2
            sage: f.derivative(y, 3) # with respect to y three times
            (6*t^2 + O(t^3))*x^2
            sage: f.derivative()    # can't figure out the variable
            Traceback (most recent call last):
            ...
            ValueError: must specify which variable to differentiate with respect to

        Polynomials over the symbolic ring (just for fun....)::

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
        """
        Let var be one of the variables of the parent of self.  This
        returns self viewed as a univariate polynomial in var over the
        polynomial ring generated by all the other variables of the parent.

        EXAMPLES::

            sage: R.<x,w,z> = QQ[]
            sage: f = x^3 + 3*w*x + w^5 + (17*w^3)*x + z^5
            sage: f.polynomial(x)
            x^3 + (17*w^3 + 3*w)*x + w^5 + z^5
            sage: parent(f.polynomial(x))
            Univariate Polynomial Ring in x over Multivariate Polynomial Ring in w, z over Rational Field

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
            sage: R.<x,y>=GF(5)[]
            sage: f=x^2+x+y
            sage: f.polynomial(x)
            x^2 + x + y
            sage: f.polynomial(y)
            y + x^2 + x
        """
        cdef int ind
        R = self.parent()
        G = R.gens()
        Z = list(G)
        try:
            ind = Z.index(var)
        except ValueError:
            raise ValueError("var must be one of the generators of the parent polynomial ring.")

        if R.ngens() <= 1:
            return self.univariate_polynomial()

        other_vars = Z
        del other_vars[ind]

        # Make polynomial ring over all variables except var.
        S = R.base_ring()[tuple(other_vars)]
        ring = S[var]
        if not self:
            return ring(0)

        d = self.degree(var)
        B = ring.base_ring()
        w = dict([(remove_from_tuple(e, ind), val) for e, val in self.dict().iteritems() if not e[ind]])
        v = [B(w)]  # coefficients that don't involve var
        z = var
        for i in range(1,d+1):
            c = self.coefficient(z).dict()
            w = dict([(remove_from_tuple(e, ind), val) for e, val in c.iteritems()])
            v.append(B(w))
            z *= var
        return ring(v)

    def _mpoly_dict_recursive(self, vars=None, base_ring=None):
        """
        Return a dict of coefficient entries suitable for construction of a MPolynomial_polydict
        with the given variables.

        EXAMPLES::

            sage: R = Integers(10)['x,y,z']['t,s']
            sage: t,s = R.gens()
            sage: x,y,z = R.base_ring().gens()
            sage: (x+y+2*z*s+3*t)._mpoly_dict_recursive(['z','t','s'])
            {(0, 0, 0): x + y, (0, 1, 0): 3, (1, 0, 1): 2}

        TESTS::

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
            sage: R = Qp(7)['x,y,z,p']; S = ZZ['x']['y,z,t']['p'] # shouldn't work, but should throw a better error
            sage: R(S.0)
            p

        See :trac:`2601`::

            sage: R.<a,b,c> = PolynomialRing(QQ, 3)
            sage: a._mpoly_dict_recursive(['c', 'b', 'a'])
            {(0, 0, 1): 1}
            sage: testR.<a,b,c> = PolynomialRing(QQ,3)
            sage: id_ringA = ideal([a^2-b,b^2-c,c^2-a])
            sage: id_ringB = ideal(id_ringA.gens()).change_ring(PolynomialRing(QQ,'c,b,a'))
        """
        from polydict import ETuple
        if not self:
            return {}

        if vars is None:
            vars = self.parent().variable_names_recursive()
        vars = list(vars)
        my_vars = self.parent().variable_names()
        if vars == list(my_vars):
            return self.dict()
        elif not my_vars[-1] in vars:
            x = base_ring(self) if base_ring is not None else self
            const_ix = ETuple((0,)*len(vars))
            return { const_ix: x }
        elif not set(my_vars).issubset(set(vars)):
            # we need to split it up
            return self.polynomial(self.parent().gen(len(my_vars)-1))._mpoly_dict_recursive(vars, base_ring)
        else:
            D = {}
            m = min([vars.index(z) for z in my_vars])
            prev_vars = vars[:m]
            var_range = range(len(my_vars))
            if len(prev_vars) > 0:
                mapping = [vars.index(v) - len(prev_vars) for v in my_vars]
                tmp = [0] * (len(vars) - len(prev_vars))
                try:
                    for ix,a in self.dict().iteritems():
                        for k in var_range:
                            tmp[mapping[k]] = ix[k]
                        postfix = ETuple(tmp)
                        mpoly = a._mpoly_dict_recursive(prev_vars, base_ring)
                        for prefix,b in mpoly.iteritems():
                            D[prefix+postfix] = b
                    return D

                except AttributeError:
                    pass

            if base_ring is self.base_ring():
                base_ring = None

            mapping = [vars.index(v) for v in my_vars]
            tmp = [0] * len(vars)
            for ix,a in self.dict().iteritems():
                for k in var_range:
                    tmp[mapping[k]] = ix[k]
                if base_ring is not None:
                    a = base_ring(a)
                D[ETuple(tmp)] = a
            return D

    cdef long _hash_c(self) except -1:
        """
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

        Verify that :trac:`16251` has been resolved, i.e., polynomials with
        unhashable coefficients are unhashable::

            sage: K.<a> = Qq(9)
            sage: R.<t,s> = K[]
            sage: hash(t)
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
    # __richcmp__.  The python documentation at  http://docs.python.org/api/type-structs.html
    # explains how __richcmp__, __hash__, and __cmp__ are tied together.
    def __hash__(self):
        return self._hash_c()

    def args(self):
        r"""
        Returns the named of the arguments of self, in the
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

        if isinstance(var, basestring):
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
        Return ``True`` if self is a homogeneous polynomial.

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

    cpdef _mod_(self, other):
        """
        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = (x^2*y + 2*x - 3)
            sage: g = (x + 1)*f
            sage: g % f
            0

            sage: (g+1) % f
            1

            sage: M = x*y
            sage: N = x^2*y^3
            sage: M.divides(N)
            True
        """
        q,r = self.quo_rem(other)
        return r

    def change_ring(self, R):
        """
        Return a copy of this polynomial but with coefficients in ``R``,
        if at all possible.

        INPUT:

        - ``R`` -- a ring or morphism.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = x^3 + 3/5*y + 1
            sage: f.change_ring(GF(7))
            x^3 + 2*y + 1

        ::

            sage: R.<x,y> = GF(9,'a')[]
            sage: (x+2*y).change_ring(GF(3))
            x - y

        ::

            sage: K.<z> = CyclotomicField(3)
            sage: R.<x,y> = K[]
            sage: f = x^2 + z*y
            sage: f.change_ring(K.embeddings(CC)[1])
            x^2 + (-0.500000000000000 + 0.866025403784439*I)*y
        """
        if isinstance(R, Morphism):
        #if we're given a hom of the base ring extend to a poly hom
            if R.domain() == self.base_ring():
                R = self.parent().hom(R, self.parent().change_ring(R.codomain()))
            return R(self)
        else:
            return self.parent().change_ring(R)(self)

    def _gap_(self, gap):
        """
        Return a representation of ``self`` in the GAP interface

        INPUT:

        - ``gap`` -- a GAP or libgap instance

        TESTS:

        Multivariate polynomial over integers::

            sage: R.<x,y,z> = ZZ[]
            sage: gap(-x*y + 3*z)   # indirect doctest
            -x*y+3*z
            sage: gap(R.zero())     # indirect doctest
            0
            sage: (x+y+z)._gap_(libgap)
            x+y+z

            sage: g = gap(x - y + 3*x*y*z)
            sage: R(g)
            3*x*y*z + x - y

            sage: g = libgap(5*x - y*z)
            sage: R(g)
            -y*z + 5*x

        Multivariate polynomial over a cyclotomic field::

            sage: F.<zeta> = CyclotomicField(8)
            sage: P.<x,y> = F[]
            sage: p = zeta + zeta^2*x + zeta^3*y + (1+zeta)*x*y
            sage: gap(p)     # indirect doctest
            (1+E(8))*x*y+E(4)*x+E(8)^3*y+E(8)
            sage: libgap(p)  # indirect doctest
            (1+E(8))*x*y+E(4)*x+E(8)^3*y+E(8)

        Multivariate polynomial over a polynomial ring over a cyclotomic field::

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
            sage: libgap(-x*y + 3*z)   # indirect doctest
            -x*y+3*z
            sage: libgap(R.zero())     # indirect doctest
            0
        """
        from sage.libs.gap.libgap import libgap
        return self._gap_(libgap)

    def _magma_init_(self, magma):
        """
        Returns a Magma string representation of self valid in the
        given magma session.

        EXAMPLES::

            sage: k.<b> = GF(25); R.<x,y> = k[]
            sage: f = y*x^2*b + x*(b+1) + 1
            sage: magma = Magma()                       # so var names same below
            sage: magma(f)                              # optional - magma
            b*x^2*y + b^22*x + 1
            sage: f._magma_init_(magma)                 # optional - magma
            '_sage_[...]!((_sage_[...]!(_sage_[...]))*_sage_[...]^2*_sage_[...]+(_sage_[...]!(_sage_[...] + 1))*_sage_[...]+(_sage_[...]!(1))*1)'

        A more complicated nested example::

            sage: R.<x,y> = QQ[]; S.<z,w> = R[]; f = (2/3)*x^3*z + w^2 + 5
            sage: f._magma_init_(magma)               # optional - magma
            '_sage_[...]!((_sage_[...]!((1/1)*1))*_sage_[...]^2+(_sage_[...]!((2/3)*_sage_[...]^3))*_sage_[...]+(_sage_[...]!((5/1)*1))*1)'
            sage: magma(f)                            # optional - magma
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


    def gradient(self):
        r"""
        Return a list of partial derivatives of this polynomial,
        ordered by the variables of ``self.parent()``.

        EXAMPLES::

           sage: P.<x,y,z> = PolynomialRing(ZZ,3)
           sage: f = x*y + 1
           sage: f.gradient()
           [y, x, 0]
        """
        return [ self.derivative(var) for var in self.parent().gens() ]

    def jacobian_ideal(self):
        r"""
        Return the Jacobian ideal of the polynomial self.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = x^3 + y^3 + z^3
            sage: f.jacobian_ideal()
            Ideal (3*x^2, 3*y^2, 3*z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.parent().ideal(self.gradient())

    def newton_polytope(self):
        """
        Return the Newton polytope of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = 1 + x*y + x^3 + y^3
            sage: P = f.newton_polytope()
            sage: P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: P.is_simple()
            True

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: R(0).newton_polytope()
            The empty polyhedron in ZZ^0
            sage: R(1).newton_polytope()
            A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: R(x^2+y^2).newton_polytope().integral_points()
            ((0, 2), (1, 1), (2, 0))
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        e = self.exponents()
        P = Polyhedron(vertices = e, base_ring=ZZ)
        return P

    def __iter__(self):
        """
        Facilitates iterating over the monomials of self,
        returning tuples of the form ``(coeff, mon)`` for each
        non-zero monomial.

        .. NOTE::

            This function creates the entire list upfront because Cython
            doesn't (yet) support iterators.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3)
            sage: f = 3*x^3*y + 16*x + 7
            sage: [(c,m) for c,m in f]
            [(3, x^3*y), (16, x), (7, 1)]
            sage: f = P.random_element(12,14)
            sage: sum(c*m for c,m in f) == f
            True
        """
        L = zip(self.coefficients(), self.monomials())
        return iter(L)

    def content(self):
        """
        Returns the content of this polynomial.  Here, we define content as
        the gcd of the coefficients in the base ring.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = 4*x+6*y
            sage: f.content()
            2
            sage: f.content().parent()
            Integer Ring

        TESTS:

        Since :trac:`10771`, the gcd in QQ restricts to the gcd in ZZ::

            sage: R.<x,y> = QQ[]
            sage: f = 4*x+6*y
            sage: f.content(); f.content().parent()
            2
            Rational Field

        """
        from sage.arith.all import gcd
        return gcd(self.coefficients())

    def is_generator(self):
        r"""
        Returns ``True`` if this polynomial is a generator of its
        parent.

        EXAMPLES::

            sage: R.<x,y>=ZZ[]
            sage: x.is_generator()
            True
            sage: (x+y-y).is_generator()
            True
            sage: (x*y).is_generator()
            False
            sage: R.<x,y>=QQ[]
            sage: x.is_generator()
            True
            sage: (x+y-y).is_generator()
            True
            sage: (x*y).is_generator()
            False
        """
        return (self in self.parent().gens())

    def map_coefficients(self, f, new_base_ring=None):
        """
        Returns the polynomial obtained by applying ``f`` to the non-zero
        coefficients of self.

        If ``f`` is a :class:`sage.categories.map.Map`, then the resulting
        polynomial will be defined over the codomain of ``f``. Otherwise, the
        resulting polynomial will be over the same ring as self. Set
        ``new_base_ring`` to override this behaviour.

        INPUT:

        - ``f`` -- a callable that will be applied to the coefficients of self.

        - ``new_base_ring`` (optional) -- if given, the resulting polynomial
          will be defined over this ring.

        EXAMPLES::

            sage: k.<a> = GF(9); R.<x,y> = k[];  f = x*a + 2*x^3*y*a + a
            sage: f.map_coefficients(lambda a : a + 1)
            (-a + 1)*x^3*y + (a + 1)*x + (a + 1)

        Examples with different base ring::

            sage: R.<r> = GF(9); S.<s> = GF(81)
            sage: h = Hom(R,S)[0]; h
            Ring morphism:
              From: Finite Field in r of size 3^2
              To:   Finite Field in s of size 3^4
              Defn: r |--> 2*s^3 + 2*s^2 + 1
            sage: T.<X,Y> = R[]
            sage: f = r*X+Y
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
        """
        Given a multivariate polynomial over a nonprime finite field
        `\GF{p**e}`, compute the norm of the polynomial down to `\GF{p}`, which
        is the product of the conjugates by the Frobenius action on
        coefficients, where Frobenius acts by p-th power.

        This is (currently) an internal function used in factoring over finite
        fields.

        EXAMPLES::

            sage: k.<a> = GF(9)
            sage: R.<x,y> = PolynomialRing(k)
            sage: f = (x-a)*(y-a)
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
        """
        Given two nonzero polynomials self and right, returns the Sylvester
        matrix of the polynomials with respect to a given variable.

        Note that the Sylvester matrix is not defined if one of the polynomials
        is zero.

        INPUT:

        - self , right: multivariate polynomials
        - variable: optional, compute the Sylvester matrix with respect to this
          variable. If variable is not provided, the first variable of the
          polynomial ring is used.

        OUTPUT:

        - The Sylvester matrix of self and right.

        EXAMPLES::

            sage: R.<x, y> = PolynomialRing(ZZ)
            sage: f = (y + 1)*x + 3*x**2
            sage: g = (y + 2)*x + 4*x**2
            sage: M = f.sylvester_matrix(g, x)
            sage: M
            [    3 y + 1     0     0]
            [    0     3 y + 1     0]
            [    4 y + 2     0     0]
            [    0     4 y + 2     0]

        If the polynomials share a non-constant common factor then the
        determinant of the Sylvester matrix will be zero::

            sage: M.determinant()
            0

            sage: f.sylvester_matrix(1 + g, x).determinant()
            y^2 - y + 7

        If both polynomials are of positive degree with respect to variable, the
        determinant of the Sylvester matrix is the resultant::

            sage: f = R.random_element(4)
            sage: g = R.random_element(4)
            sage: f.sylvester_matrix(g, x).determinant() == f.resultant(g, x)
            True

        TEST:

        The variable is optional::

            sage: f = x + y
            sage: g = x + y
            sage: f.sylvester_matrix(g)
            [1 y]
            [1 y]

        Polynomials must be defined over compatible base rings::

            sage: K.<x, y> = QQ[]
            sage: f = x + y
            sage: L.<x, y> = ZZ[]
            sage: g = x + y
            sage: R.<x, y> = GF(25, 'a')[]
            sage: h = x + y
            sage: f.sylvester_matrix(g, 'x')
            [1 y]
            [1 y]
            sage: g.sylvester_matrix(h, 'x')
            [1 y]
            [1 y]
            sage: f.sylvester_matrix(h, 'x')
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Multivariate Polynomial Ring in x, y over Rational Field' and 'Multivariate Polynomial Ring in x, y over Finite Field in a of size 5^2'
            sage: K.<x, y, z> = QQ[]
            sage: f = x + y
            sage: L.<x, z> = QQ[]
            sage: g = x + z
            sage: f.sylvester_matrix(g)
            [1 y]
            [1 z]

        Corner cases::

            sage: K.<x ,y>=QQ[]
            sage: f = x^2+1
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

    def macaulay_resultant(self, *args):
        r"""
        This is an implementation of the Macaulay Resultant. It computes
        the resultant of universal polynomials as well as polynomials
        with constant coefficients. This is a project done in
        sage days 55. It's based on the implementation in Maple by
        Manfred Minimair, which in turn is based on the references [CLO], [Can], [Mac].
        It calculates the Macaulay resultant for a list of Polynomials,
        up to sign!

        AUTHORS:

        - Hao Chen, Solomon Vishkautsan (7-2014)

        INPUT:

        - ``args`` -- a list of `n-1` homogeneous polynomials in `n` variables.
                  works when ``args[0]`` is the list of polynomials,
                  or ``args`` is itself the list of polynomials

        OUTPUT:

        - the macaulay resultant

        EXAMPLES:

        The number of polynomials has to match the number of variables::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: y.macaulay_resultant(x+z)
            Traceback (most recent call last):
            ...
            TypeError: number of polynomials(= 2) must equal number of variables (= 3)

        The polynomials need to be all homogeneous::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: y.macaulay_resultant([x+z, z+x^3])
            Traceback (most recent call last):
            ...
            TypeError: resultant for non-homogeneous polynomials is not supported

        All polynomials must be in the same ring::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: S.<x,y> = PolynomialRing(QQ, 2)
            sage: y.macaulay_resultant(z+x,z)
            Traceback (most recent call last):
            ...
            TypeError: not all inputs are polynomials in the calling ring

        The following example recreates Proposition 2.10 in Ch.3 of Using Algebraic Geometry::

            sage: K.<x,y> = PolynomialRing(ZZ, 2)
            sage: flist,R = K._macaulay_resultant_universal_polynomials([1,1,2])
            sage: flist[0].macaulay_resultant(flist[1:])
            u2^2*u4^2*u6 - 2*u1*u2*u4*u5*u6 + u1^2*u5^2*u6 - u2^2*u3*u4*u7 + u1*u2*u3*u5*u7 + u0*u2*u4*u5*u7 - u0*u1*u5^2*u7 + u1*u2*u3*u4*u8 - u0*u2*u4^2*u8 - u1^2*u3*u5*u8 + u0*u1*u4*u5*u8 + u2^2*u3^2*u9 - 2*u0*u2*u3*u5*u9 + u0^2*u5^2*u9 - u1*u2*u3^2*u10 + u0*u2*u3*u4*u10 + u0*u1*u3*u5*u10 - u0^2*u4*u5*u10 + u1^2*u3^2*u11 - 2*u0*u1*u3*u4*u11 + u0^2*u4^2*u11

        The following example degenerates into the determinant of a `3*3` matrix::

            sage: K.<x,y> = PolynomialRing(ZZ, 2)
            sage: flist,R = K._macaulay_resultant_universal_polynomials([1,1,1])
            sage: flist[0].macaulay_resultant(flist[1:])
            -u2*u4*u6 + u1*u5*u6 + u2*u3*u7 - u0*u5*u7 - u1*u3*u8 + u0*u4*u8

        The following example is by Patrick Ingram(arxiv:1310.4114)::

            sage: U = PolynomialRing(ZZ,'y',2); y0,y1 = U.gens()
            sage: R = PolynomialRing(U,'x',3); x0,x1,x2 = R.gens()
            sage: f0 = y0*x2^2 - x0^2 + 2*x1*x2
            sage: f1 = y1*x2^2 - x1^2 + 2*x0*x2
            sage: f2 = x0*x1 - x2^2
            sage: f0.macaulay_resultant(f1,f2)
            y0^2*y1^2 - 4*y0^3 - 4*y1^3 + 18*y0*y1 - 27

        a simple example with constant rational coefficients::

            sage: R.<x,y,z,w> = PolynomialRing(QQ,4)
            sage: w.macaulay_resultant([z,y,x])
            1

        an example where the resultant vanishes::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: (x+y).macaulay_resultant([y^2,x])
            0

        an example of bad reduction at a prime ``p = 5``::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: y.macaulay_resultant([x^3+25*y^2*x,5*z])
            125

        The input can given as an unpacked list of polynomials::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: y.macaulay_resultant(x^3+25*y^2*x,5*z)
            125

        an example when the coefficients live in a finite field::

            sage: F = FiniteField(11)
            sage: R.<x,y,z,w> = PolynomialRing(F,4)
            sage: z.macaulay_resultant([x^3,5*y,w])
            4

        example when the denominator in the algorithm vanishes(in this case
        the resultant is the constant term of the quotient of
        char polynomials of numerator/denominator)::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: y.macaulay_resultant([x+z, z^2])
            -1

        when there are only 2 polynomials, macaulay resultant degenerates to the traditional resultant::

            sage: R.<x> = PolynomialRing(QQ,1)
            sage: f =  x^2+1; g = x^5+1
            sage: fh = f.homogenize()
            sage: gh = g.homogenize()
            sage: RH = fh.parent()
            sage: f.resultant(g) == fh.macaulay_resultant(gh)
            True

        """
        if len(args) == 1 and isinstance(args[0],list):
            return self.parent().macaulay_resultant(self, *args[0])
        return self.parent().macaulay_resultant(self, *args)

    def denominator(self):
        """
        Return a denominator of self.

        First, the lcm of the denominators of the entries of self
        is computed and returned. If this computation fails, the
        unit of the parent of self is returned.

        Note that some subclases may implement its own denominator
        function.

        .. warning::

           This is not the denominator of the rational function
           defined by self, which would always be 1 since self is a
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

            sage: R.<x,y> = NumberField(symbolic_expression(x^2+3)  ,'a')['x,y']
            sage: f = (1/17)*x^19 + (1/6)*y - (2/3)*x + 1/3; f
            1/17*x^19 - 2/3*x + 1/6*y + 1/3
            sage: f.denominator()
            102

        Finally, we try to compute the denominator of a polynomial with
        coefficients in the real numbers, which is a ring whose elements do
        not have a denominator method.

        ::

            sage: R.<a,b,c> = RR[]
            sage: f = a + b + RR('0.3'); f
            a + b + 0.300000000000000
            sage: f.denominator()
            1.00000000000000

        Check that the denominator is an element over the base whenever the base
        has no denominator function. This closes #9063

        ::

            sage: R.<a,b,c> = GF(5)[]
            sage: x = R(0)
            sage: x.denominator()
            1
            sage: type(x.denominator())
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: type(a.denominator())
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
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
        """
        Return a numerator of self computed as self * self.denominator()

        Note that some subclases may implement its own numerator
        function.

        .. warning::

           This is not the numerator of the rational function
           defined by self, which would always be self since self is a
           polynomial.

        EXAMPLES:

        First we compute the numerator of a polynomial with
        integer coefficients, which is of course self.

        ::

            sage: R.<x, y> = ZZ[]
            sage: f = x^3 + 17*x + y + 1
            sage: f.numerator()
            x^3 + 17*x + y + 1
            sage: f == f.numerator()
            True

        Next we compute the numerator of a polynomial over a number field.

        ::

            sage: R.<x,y> = NumberField(symbolic_expression(x^2+3)  ,'a')['x,y']
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
        are valid

        ::

            sage: K=NumberField(symbolic_expression('x^3+2'),'a')['x']['s,t']
            sage: f=K.random_element()
            sage: f.numerator() / f.denominator() == f
            True
            sage: R=RR['x,y,z']
            sage: f=R.random_element()
            sage: f.numerator() / f.denominator() == f
            True
        """
        return self * self.denominator()

    def lift(self, I):
        """
        given an ideal ``I = (f_1,...,f_r)`` and some ``g (== self)`` in ``I``,
        find ``s_1,...,s_r`` such that ``g = s_1 f_1 + ... + s_r f_r``.

        EXAMPLE::

            sage: A.<x,y> = PolynomialRing(CC,2,order='degrevlex')
            sage: I = A.ideal([x^10 + x^9*y^2, y^8 - x^2*y^7 ])
            sage: f = x*y^13 + y^12
            sage: M = f.lift(I)
            sage: M
            [y^7, x^7*y^2 + x^8 + x^5*y^3 + x^6*y + x^3*y^4 + x^4*y^2 + x*y^5 + x^2*y^3 + y^4]
            sage: sum( map( mul , zip( M, I.gens() ) ) ) == f
            True
        """
        raise NotImplementedError

    def inverse_mod(self, I):
        """
        Returns an inverse of self modulo the polynomial ideal `I`,
        namely a multivariate polynomial `f` such that
        ``self * f - 1`` belongs to `I`.

        INPUT:
         - ``I`` -- an ideal of the polynomial ring in which self lives

        OUTPUT:

         - a multivariate polynomial representing the inverse of ``f`` modulo ``I``

        EXAMPLES::

           sage: R.<x1,x2> = QQ[]
           sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
           sage: f = x1 + 3*x2^2; g = f.inverse_mod(I); g
           1/16*x1 + 3/16
           sage: (f*g).reduce(I)
           1

        Test a non-invertible element::

           sage: R.<x1,x2> = QQ[]
           sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
           sage: f = x1 + x2
           sage: f.inverse_mod(I)
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
        """
        Return the weighted degree of ``self``, which is the maximum weighted
        degree of all monomials in ``self``; the weighted degree of a monomial
        is the sum of all powers of the variables in the monomial, each power
        multiplied with its respective weight in ``weights``.

        This method is given for convenience. It is faster to use polynomial
        rings with weighted term orders and the standard ``degree`` function.

        INPUT:

        - ``weights`` - Either individual numbers, an iterable or a dictionary,
          specifying the weights of each variable. If it is a dictionary, it
          maps each variable of ``self`` to its weight. If it is a sequence of
          individual numbers or a tuple, the weights are specified in the order
          of the generators as given by ``self.parent().gens()``:

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
            sage: q = R.random_element(100, 20) #random
            sage: q.weighted_degree(1, 1, 1) == q.total_degree()
            True

        You may also work with negative weights

        ::

            sage: p.weighted_degree(-1, -2, -1)
            -2

        Note that only integer weights are allowed

        ::

            sage: p.weighted_degree(x,1,1)
            Traceback (most recent call last):
            ...
            TypeError
            sage: p.weighted_degree(2/1,1,1)
            6

        The ``weighted_degree`` coincides with the ``degree`` of a weighted
        polynomial ring, but the later is faster.

        ::

            sage: K = PolynomialRing(QQ, 'x,y', order=TermOrder('wdegrevlex', (2,3)))
            sage: p = K.random_element(10)
            sage: p.degree() == p.weighted_degree(2,3)
            True

        TESTS::

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
        """
        Return a greatest common divisor of this polynomial and ``other``.

        INPUT:

        - ``other`` -- a polynomial with the same parent as this polynomial

        EXAMPLES::

            sage: Q.<z> = Frac(QQ['z'])
            sage: R.<x,y> = Q[]
            sage: r = x*y - (2*z-1)/(z^2+z+1) * x + y/z
            sage: p = r * (x + z*y - 1/z^2)
            sage: q = r * (x*y*z + 1)
            sage: gcd(p,q)
            (z^3 + z^2 + z)*x*y + (-2*z^2 + z)*x + (z^2 + z + 1)*y

        Polynomials over polynomial rings are converted to a simpler polynomial
        ring with all variables to compute the gcd::

            sage: A.<z,t> = ZZ[]
            sage: B.<x,y> = A[]
            sage: r = x*y*z*t+1
            sage: p = r * (x - y + z - t + 1)
            sage: q = r * (x*z - y*t)
            sage: gcd(p,q)
            z*t*x*y + 1
            sage: _.parent()
            Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring in z, t over Integer Ring

        Some multivariate polynomial rings have no gcd implementation::

            sage: R.<x,y> =GaussianIntegers()[]
            sage: x.gcd(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: GCD is not implemented for multivariate polynomials over Gaussian Integers in Number Field in I with defining polynomial x^2 + 1
        """
        variables = self._parent.variable_names_recursive()
        if len(variables) > self._parent.ngens():
            base = self._parent._mpoly_base_ring()
            d1 = self._mpoly_dict_recursive()
            d2 = other._mpoly_dict_recursive()
            ring = PolynomialRing(base, variables)
            try:
                return self._parent(ring(d1).gcd(ring(d2)))
            except (AttributeError, NotImplementedError):
                pass

        try:
            self._parent._singular_().set_ring()
            g = self._singular_().gcd(other._singular_())
            return self._parent(g)
        except (TypeError, AttributeError):
            pass

        x = self._parent.gens()[-1]
        uniself = self.polynomial(x)
        unibase = uniself.base_ring()
        if hasattr(unibase, "_gcd_univariate_polynomial"):
            return self._parent(unibase._gcd_univariate_polynomial(uniself, other.polynomial(x)))
        else:
            raise NotImplementedError("GCD is not implemented for multivariate polynomials over {}".format(self._parent._mpoly_base_ring()))

    def nth_root(self, n):
        r"""
        Return a `n`-th root of this element.

        This method relies on factorization.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: a = 32 * (x*y + 1)^5 * (x+y+z)^5
            sage: a.nth_root(5)
            2*x^2*y + 2*x*y^2 + 2*x*y*z + 2*x + 2*y + 2*z
            sage: b = x + 2*y + 3*z
            sage: b.nth_root(42)
            Traceback (most recent call last):
            ...
            ValueError: (x + 2*y + 3*z)^(1/42) does not lie in
            Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        # note: this code is duplicated in
        # sage.rings.polynomial.polynomial_element.Polynomial.nth_root
        from sage.rings.integer_ring import ZZ

        n = ZZ.coerce(n)

        if n <= 0:
            raise ValueError("n (={}) must be positive".format(n))
        elif n.is_one() or self.is_zero():
            return self
        elif self.degree() % n:
            raise ValueError("({})^(1/{}) does not lie in {}".format(self, n, self.parent()))
        else:
            f = self.factor()
            u = self.base_ring()(f.unit())

            if u.is_one():
                ans = self.parent().one()
            else:
                # try to compute a n-th root of the unit in the
                # base ring. the `nth_root` method thus has to be
                # implemented in the base ring.
                try:
                    ans = self.parent(u.nth_root(n))
                except AttributeError:
                    raise NotImplementedError("nth root not implemented for {}".format(u.parent()))

            for (v, exp) in f:
                if exp % n:
                    raise ValueError("({})^(1/{}) does not lie in {}".format(self, n, self.parent()))
                ans *= v ** (exp // n)

            return ans

    def reduced_form(self, prec=300, return_conjugation=True, error_limit=0.000001):
        r"""
        Returns a reduced form of this polynomial.

        The algorithm is from Stoll and Cremona's "On the Reduction Theory of Binary Forms" [SC]_.
        This takes a two variable homogenous polynomial and finds a reduced form. This is a
        `SL(2,\ZZ)`-equivalent binary form whose covariant in the upper half plane is in the fundamental
        domain. This should also minimize the sum of the squares of the coefficients,
        but this is not always the case.

        A portion of the algorithm uses Newton's method to find a solution to a system of equations.
        If Newton's method fails to converge to a point in the upper half plane, the function
        will use the less precise `Q_0` covariant as defined in [SC]_. Additionally, if this polynomial has
        a root with multiplicity at lease half the total degree of the polynomial, then
        we must also use the `Q_0` covariant. See [SC]_ for details.

        Note that, if the covariant is within ``error_limit`` of the boundry but outside
        the fundamental domain, our function will erroneously move it to within the
        fundamental domain, hence our conjugation will be off by 1. If you don't want
        this to happen, decrease your ``error_limit`` and increase your precision.

        Implemented by Rebecca Lauren Miller as part of GSOC 2016.

        INPUT:

        - ``prec`` --  integer, sets the precision (default:300)

        - ``return_conjugation`` -- boolean. Returns element of `SL(2, \ZZ)` (default:True)

        - ``error_limit`` -- sets the error tolerance (default:0.000001)

        OUTPUT:

            - a polynomial (reduced binary form)

            - a matrix (element of `SL(2, \ZZ)`)

        TODO: When Newton's Method doesn't converge to a root in the upper half plane.
            Now we just return z0. It would be better to modify and find the unique root
            in the upper half plane.


        REFERENCES:

        .. [SC] Michael Stoll and John E. Cremona. On The Reduction Theory of Binary Forms.
            Journal für die reine und angewandte Mathematik, 565 (2003), 79-99.

        EXAMPLES::

            sage: R.<x,h> = PolynomialRing(QQ)
            sage: f = 19*x^8 - 262*x^7*h + 1507*x^6*h^2 - 4784*x^5*h^3 + 9202*x^4*h^4\
             -10962*x^3*h^5 + 7844*x^2*h^6 - 3040*x*h^7 + 475*h^8
            sage: f.reduced_form(prec=200)
            (
            -x^8 - 2*x^7*h + 7*x^6*h^2 + 16*x^5*h^3 + 2*x^4*h^4 - 2*x^3*h^5 + 4*x^2*h^6 - 5*h^8,
            <BLANKLINE>
            [ 1 -2]
            [ 1 -1]
            )

        An example were the multiplicity is too high::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = x^3 + 378666*x^2*y - 12444444*x*y^2 + 1234567890*y^3
            sage: j = f * (x-545*y)^9
            sage: j.reduced_form(prec=200)
            (
            x^12 + 374553*x^11*y - 1587470292*x^10*y^2 + 2960311881270*x^9*y^3 - 3189673382015880*x^8*y^4
            + 2180205736473134502*x^7*y^5 - 972679603186995463284*x^6*y^6 + 278555935048988817910176*x^5*y^7
            - 47339497613591564056277355*x^4*y^8 + 3719790227462793441137663545*x^3*y^9
            + 4017321423785434880978464176*x^2*y^10 + 1605293849731195593699202674738*x*y^11
            - 2738526775493743375819069013598582*y^12,
            <BLANKLINE>
            [ 1 66]
            [ 0  1]
            )

        An example where Newton's Method doesnt find the right root::

            sage: R.<x,h> = PolynomialRing(QQ)
            sage: f = 234*x^11*h + 104832*x^10*h^2 + 21346884*x^9*h^3 + 2608021728*x^8*h^4\
             + 212413000410*x^7*h^5 + 12109691106162*x^6*h^6 + 493106447396862*x^5*h^7\
             + 14341797993350646*x^4*h^8 + 291976289803277118*x^3*h^9 +3962625618555930456*x^2*h^10\
             + 32266526239647689652*x*h^11 + 119421058057217196228*h^12
            sage: f.reduced_form(prec=600) # long time
            (
            234*x^11*h - 702*x^10*h^2 + 234*x^9*h^3 - 1638*x^8*h^4 + 17550*x^7*h^5 - 35568*x^6*h^6
            - 42120*x^5*h^7 - 248508*x^4*h^8 + 35802*x^3*h^9 + 23868*x^2*h^10 - 936*x*h^11 - 468*h^12,
            <BLANKLINE>
            [  1 -41]
            [  0   1]
            )

        An example with covariant on the boundary, therefore a non-unique form also a_0 is 0::

            sage: R.<x,h> = PolynomialRing(QQ)
            sage: g = -1872*x^5*h - 1375452*x^4*h^2 - 404242956*x^3*h^3 - 59402802888*x^2*h^4\
             -4364544068352*x*h^5 - 128270946360960*h^6
            sage: g.reduced_form()
            (
            -1872*x^5*h + 468*x^4*h^2 + 2340*x^3*h^3 - 2340*x^2*h^4 - 468*x*h^5 + 1872*h^6,
            <BLANKLINE>
            [ -1 147]
            [  0  -1]
            )

        An example where precision needs to be increased::

            sage: R.<x,h> = PolynomialRing(QQ)
            sage: f = -1872*x^5*h - 1375452*x^4*h^2 - 404242956*x^3*h^3 - 59402802888*x^2*h^4\
             -4364544068352*x*h^5 - 128270946360960*h^6
            sage: f.reduced_form(prec=200)
            Traceback (most recent call last):
            ...
            ValueError: accuracy of Newton's root not within tolerance(1.5551623876686905873160660564410782587973928631765344695031 > 1e-06), increase precision
            sage: f.reduced_form(prec=400)
            (
            -1872*x^5*h + 468*x^4*h^2 + 2340*x^3*h^3 - 2340*x^2*h^4 - 468*x*h^5 + 1872*h^6,
            <BLANKLINE>
            [ -1 147]
            [  0  -1]
            )

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: F = - 8*x^4 - 3933*x^3*y - 725085*x^2*y^2 - 59411592*x*y^3 - 1825511633*y^4
            sage: F.reduced_form(return_conjugation=False)
            x^4 + 9*x^3*y - 3*x*y^3 - 8*y^4

        ::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: F = x^4 + x^3*y*z + y^2*z
            sage: F.reduced_form()
            Traceback (most recent call last):
            ...
            ValueError: (=x^3*y*z + x^4 + y^2*z) must have two variables

        ::

            sage: R.<x,y> = PolynomialRing(ZZ)
            sage: F = - 8*x^6 - 3933*x^3*y - 725085*x^2*y^2 - 59411592*x*y^3 - 99*y^6
            sage: F.reduced_form(return_conjugation=False)
            Traceback (most recent call last):
            ...
            ValueError: (=-8*x^6 - 99*y^6 - 3933*x^3*y - 725085*x^2*y^2 -
            59411592*x*y^3) must be homogenous

        ::

            sage: R.<x,y> = PolynomialRing(RR)
            sage: F = 217.992172373276*x^3 + 96023.1505442490*x^2*y + 1.40987971253579e7*x*y^2\
            + 6.90016027113216e8*y^3
            sage: F.reduced_form()
            (
            -39.5673942565918*x^3 + 111.874026298523*x^2*y + 231.052762985229*x*y^2 - 138.380829811096*y^3,
            <BLANKLINE>
            [-147 -148]
            [   1    1]
            )

        ::

            sage: R.<x,y> = PolynomialRing(CC)
            sage: F = (0.759099196558145 + 0.845425869641446*CC.0)*x^3 + (84.8317207268542 + 93.8840848648033*CC.0)*x^2*y\
            + (3159.07040755858 + 3475.33037377779*CC.0)*x*y^2 + (39202.5965389079 + 42882.5139724962*CC.0)*y^3
            sage: F.reduced_form()
            (
            (-0.759099196558145 - 0.845425869641446*I)*x^3 + (-0.571709908900118 - 0.0418133346027929*I)*x^2*y
            + (0.856525964330103 - 0.0721403997649759*I)*x*y^2 + (-0.965531044130330 + 0.754252314465703*I)*y^3,
            <BLANKLINE>
            [-1 37]
            [ 0 -1]
            )
        """
        from sage.matrix.constructor import matrix
        from sage.calculus.functions import jacobian

        if self.parent().ngens() != 2:
            raise ValueError("(=%s) must have two variables"%self)
        if self.is_homogeneous() != True:
            raise ValueError("(=%s) must be homogenous"%self)

        #getting a numerical approximation of the roots of our polynomial
        CF = ComplexIntervalField(prec=prec) # keeps trac of our precision error
        RF = RealField(prec=prec)
        R = self.parent()
        S = PolynomialRing(R.base_ring(),'z')
        phi = R.hom([S.gen(0), 1], S)# dehomogenization
        F = phi(self).quo_rem(gcd(phi(self), phi(self).derivative()))[0] # removes multiple roots
        from sage.rings.polynomial.complex_roots import complex_roots
        roots = [p for p,e in complex_roots(F, min_prec=prec)]
        #roots = F.roots(ring=CF, multiplicities=False)

        #finding quadratic Q_0, gives us our convariant, z_0
        dF = F.change_ring(CF).derivative()
        n = F.degree()
        R = PolynomialRing(CF,'x,y')
        x,y = R.gens()
        Q  = []
        # finds Stoll and Cremona's Q_0
        for j in range(len(roots)):
            k = (1/(dF(roots[j]).abs()**(2/(n-2)))) * ((x-(roots[j]*y)) * (x-(roots[j].conjugate()*y)))
            Q.append(k)
        # this is Q_o , always positive def as long as F HAS DISTINCT ROOTS
        q = sum([Q[i] for i in range(len(Q))])
        A = q.monomial_coefficient(x**2)
        B = q.monomial_coefficient(x*y)
        C = q.monomial_coefficient(y**2)
        # need positive root this will be our first z
        try:
            z = (-B + ((B**2)-(4*A*C)).sqrt())/(2*A)# this is z_o
        except ValueError:
            raise ValueError("not enough precision")
        if z.imag()<0:
            z = (-B - ((B**2)-(4*A*C)).sqrt())/(2*A)

        # this moves z to our fundamental domain using the three steps laid out in the algorithim by [SC]
        # this is found in section 5 of their paper
        M = matrix(QQ, [[1,0], [0,1]]) # used to keep track of how our z is moved.
        zc = z.center()
        while zc.real() < RF(-0.5) or zc.real() >= RF(0.5) or (zc.real() <= RF(0) and zc.abs() < RF(1))\
         or (zc.real() > RF(0) and zc.abs() <= RF(1)):
            if zc.real() < RF(-0.5): # moves z into fundamental domain by m
                m = zc.real().abs().round() # finds amount to move z's real part by
                Qm = QQ(m)
                M = M * matrix(QQ, [[1,-Qm], [0,1]]) # move
                z += m  # M.inverse()*z is supposed to move z by m
            elif zc.real()>=RF(0.5): #moves z into fundamental domain by m
                m = zc.real().round()
                Qm = QQ(m)
                M = M * matrix(QQ, [[1,Qm], [0,1]])  #move z
                z -= m
            elif (zc.real() <= RF(0) and zc.abs() < RF(1)) or (zc.real() > RF(0) and zc.abs() <= RF(1)): # flips z
                z = -1/z
                M = M * matrix(QQ, [[0,-1], [1,0]])# multiply on left because we are taking inverse matrices
            zc = z.center()
        z0 = z
        # creates and solves equations 4.4 in [SC], gives us a new z
        x,y = self.parent().gens()
        F = S(phi(self(tuple((M * vector([x, y])))))) # New self, S pushes it to polyomial ring
        #L1 = F.roots(ring=CF, multiplicities=True)
        L1 = complex_roots(F, min_prec=prec)
        L=[]
        newton = True
        err = z.diameter()
        # making sure multiplicity isn't too large using convergence conditions in paper
        for p,e in L1:
            if e > self.degree()/2:
                newton = False
                break
            for l in range(e):
                L.append(p)
        if newton:
            a = 0
            c = 0
            RCF = PolynomialRing(CF, 'u,t')
            u,t = RCF.gens()
            for j in range(len(L)):
                b = u**2/((t-(L[j])) * (t-(L[j].conjugate()))+ u**2)
                a += b
                d = (t-(L[j].real()))/((t-(L[j])) * (t-(L[j].conjugate())) + u**2)
                c += d
            #Newton's Method, to to find solutions. Error bound is while less than diameter of our z
            err = z.diameter()
            zz = z.diameter()
            n = F.degree()
            g1 = a.numerator() - n/2*a.denominator()
            g2 = c.numerator()
            G = vector([g1, g2])
            J = jacobian(G, [u,t])
            v0 = vector([z.imag(), z.real()]) #z0 as starting point
            #finds our correct z
            while err <= zz:
                NJ = J.subs({u:v0[0], t:v0[1]})
                NJinv = NJ.inverse()
                #inverse for CIF matrix seems to return fractions not CIF elements, fix them
                if NJinv.base_ring() != CF:
                    NJinv = matrix(CF,2,2,[CF(zw.numerator()/zw.denominator()) for zw in NJinv.list()])
                w = z
                v0 = v0 - NJinv*G.subs({u:v0[0], t:v0[1]})
                z = v0[1].constant_coefficient() + v0[0].constant_coefficient()*CF.gen(0)
                err = z.diameter() # precision
                zz = (w - z).abs() #difference in w and z
            else:
                if err > error_limit:
                    raise ValueError("accuracy of Newton's root not within tolerance(%s > %s), increase precision"%(err, error_limit))
            zc = z.center()
            # moves our z to fundamental domain as before
            if  zc.imag()> z.diameter():
                while zc.real() < RF(-0.5) or zc.real() >= RF(0.5) or (zc.real() <= RF(0) and zc.abs() < RF(1))\
                 or (zc.real() > RF(0) and zc.abs() <= RF(1)):
                    if zc.real() < RF(-0.5):
                        m = zc.real().abs().round()
                        Qm = QQ(m)
                        M = M*matrix(QQ, [[1,-Qm], [0,1]])
                        z += m  # M.inverse()*Z is supposed to move z by m
                    elif zc.real() >= RF(0.5): #else if
                        m = zc.real().round()
                        Qm = QQ(m)
                        M = M * matrix(QQ, [[1,Qm], [0,1]])
                        z -= m
                    elif (zc.real() <= RF(0) and zc.abs() < RF(1)) or (zc.real() > RF(0) and zc.abs() <= RF(1)):
                        z = -1/z
                        M = M * matrix(QQ, [[0,-1],[ 1,0]])
                    zc = z.center()
            else:
                #verbose("Warning: Newton's method converged to z not in the upper half plane.", level=1)
                z = z0
        else:
            # means multiplicity of roots is too high, need to use Q0 reduced form, z0
            z=z0
        err = z.diameter()
        x,y = self.parent().gens()
        if return_conjugation:
            return (self(tuple(M * vector([x,y]))), M)
        return self(tuple(M * vector([x,y])))

cdef remove_from_tuple(e, int ind):
    w = list(e)
    del w[ind]
    if len(w) == 1:
        return w[0]
    else:
        return tuple(w)
