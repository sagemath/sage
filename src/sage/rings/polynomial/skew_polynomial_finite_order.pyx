# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.combinat sage.rings.finite_rings
r"""
Univariate dense skew polynomials over a field with a finite order automorphism

AUTHOR::

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests
  and refactored classes and methods
"""

# ***************************************************************************
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************
from sage.structure.parent cimport Parent
from sage.structure.element cimport RingElement
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.skew_polynomial_element cimport SkewPolynomial_generic_dense


cdef class SkewPolynomial_finite_order_dense(SkewPolynomial_generic_dense):
    def __init__(self, parent, x=None, int check=1, int construct=0, **kwds):
        """
        This method constructs a generic dense skew polynomial over a field equipped
        with an automorphism of finite order.

        INPUT:

        - ``parent`` -- parent of ``self``

        - ``x`` -- list of coefficients from which ``self`` can be constructed

        - ``check`` -- flag variable to normalize the polynomial

        - ``construct`` -- boolean (default: ``False``)

        TESTS::

            sage: R.<t> = GF(5^3)
            sage: Frob = R.frobenius_endomorphism()
            sage: S.<x> = R['x', Frob]; S
            Ore Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5

        We create a skew polynomial from a list::

            sage: S([t,1])
            x + t

        from another skew polynomial::

            sage: S(x^2 + t)
            x^2 + t

        from a constant::

            sage: x = S(t^2 + 1); x
            t^2 + 1
            sage: x.parent() is S
            True
        """
        SkewPolynomial_generic_dense.__init__ (self, parent, x, check, construct, **kwds)
        self._norm = None
        self._charpoly = None
        self._optbound = None

    cdef _matphir_c(self):
        r"""
        Return the matrix of the multiplication by `X^r` on
        the quotient `K[X,\sigma] / K[X,\sigma]*self`.
        """
        from sage.matrix.constructor import matrix
        parent = self._parent
        cdef Py_ssize_t i, j, col, exp, n
        cdef Py_ssize_t d = self.degree()
        cdef Py_ssize_t r = parent._order
        cdef k = parent.base_ring()
        cdef RingElement zero = k(0)
        cdef RingElement one = k(1)
        cdef list line, phir = []
        if r < d:
            for i from 0 <= i < d-r:
                line = d * [zero]
                line[r+i] = one
                phir.append(line)
            col = d-r
            exp = d
        else:
            col = 0
            exp = r
        cdef SkewPolynomial_finite_order_dense powx = <SkewPolynomial_finite_order_dense>self._new_c([zero,one], parent)
        cdef SkewPolynomial_finite_order_dense v
        if (exp % 2 == 1):
            v = <SkewPolynomial_finite_order_dense>self._new_c([zero,one], parent)
            if self.degree() == 1:
                v %= self
        else:
            v = <SkewPolynomial_finite_order_dense>self._new_c([one], parent)
        exp = exp >> 1
        n = 1
        while exp != 0:
            powx = (powx.conjugate(n) * powx) % self
            n = n << 1
            if (exp % 2 == 1):
                v = v.conjugate(n)
                v = (v * powx) % self
            exp = exp >> 1
        line = v.list()
        line += (d - len(line)) * [zero]
        phir.append(line)
        for j from col+1 <= j < d:
            v <<= 1
            v = v.conjugate(1) % self
            line = v._coeffs[:] + (d - len(v._coeffs)) * [zero]
            phir.append(line)
        return matrix(k, phir)

    cdef _matmul_c(self):
        r"""
        Return the matrix of the multiplication by ``self`` on
        `K[X,\sigma]` considered as a free module over `K[X^r]`
        (here `r` is the order of `\sigma`).
        """
        from sage.matrix.constructor import matrix
        cdef Py_ssize_t i, j, deb, k, r = self.parent()._order
        cdef Py_ssize_t d = self.degree ()
        cdef Parent base_ring = <Parent?>self.parent().base_ring()
        cdef RingElement minusone = <RingElement?>base_ring(-1)
        cdef RingElement zero = <RingElement?>base_ring(0)
        cdef Polk = PolynomialRing (base_ring, 'xr')
        cdef list M = []
        cdef list l = self.list()
        for j in range(r):
            for i in range(r):
                if i < j:
                    pol = [zero]
                    deb = i-j+r
                else:
                    pol = []
                    deb = i-j
                for k from deb <= k <= d by r:
                    pol.append(l[k])
                M.append(Polk(pol))
            for i from 0 <= i <= d:
                l[i] = self._parent.twisting_morphism()(l[i])
        return matrix(Polk, r, r, M)

    def reduced_trace(self, var=None):
        r"""
        Return the reduced trace of this skew polynomial.

        INPUT:

        - ``var`` -- string or ``False`` or ``None`` (default: ``None``);
          the variable name. If ``False``, return the list of coefficients.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: a = x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: tr = a.reduced_trace(); tr
            3*z + 4

        The reduced trace lies in the center of `S`, which is a univariate
        polynomial ring in the variable `z = x^3` over `\GF{5}`::

            sage: tr.parent()
            Univariate Polynomial Ring in z over Finite Field of size 5
            sage: tr.parent() is S.center()
            True

        We can use explicit conversion to view ``tr`` as a skew polynomial::

            sage: S(tr)
            3*x^3 + 4

        By default, the name of the central variable is usually ``z`` (see
        :meth:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_order.center`
        for more details about this).
        However, the user can specify a different variable name if desired::

            sage: a.reduced_trace(var='u')
            3*u + 4

        When passing in ``var=False``, a tuple of coefficients (instead of
        an actual polynomial) is returned::

            sage: a.reduced_trace(var=False)
            (4, 3)

        TESTS:

        We check that the reduced trace is additive::

            sage: a = S.random_element(degree=5)
            sage: b = S.random_element(degree=7)
            sage: a.reduced_trace() + b.reduced_trace() == (a+b).reduced_trace()
            True

        .. SEEALSO::

            :meth:`reduced_norm`, :meth:`reduced_charpoly`
        """
        order = self.parent()._order
        twisting_morphism = self.parent().twisting_morphism()
        coeffs = []
        for i in range(0, self.degree()+1, order):
            tr = c = self._coeffs[i]
            for _ in range(order-1):
                tr = c + twisting_morphism(tr)
            coeffs.append(tr)
        if var is False:
            return tuple(coeffs)
        Z = self.parent().center(name=var)
        return Z(coeffs)

    def reduced_norm(self, var=None):
        r"""
        Return the reduced norm of this skew polynomial.

        INPUT:

        - ``var`` -- string or ``False`` or ``None`` (default: ``None``);
          the variable name. If ``False``, return the list of coefficients.

        .. NOTE::

            The result is cached.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: a = x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: N = a.reduced_norm(); N
            z^3 + 4*z^2 + 4

        The reduced norm lies in the center of `S`, which is a univariate
        polynomial ring in the variable `z = x^3` over `\GF{5}`::

            sage: N.parent()
            Univariate Polynomial Ring in z over Finite Field of size 5
            sage: N.parent() is S.center()
            True

        We can use explicit conversion to view ``N`` as a skew polynomial::

            sage: S(N)
            x^9 + 4*x^6 + 4

        By default, the name of the central variable is usually ``z`` (see
        :meth:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_order.center`
        for more details about this).
        However, the user can specify a different variable name if desired::

            sage: a.reduced_norm(var='u')
            u^3 + 4*u^2 + 4

        When passing in ``var=False``, a tuple of coefficients (instead of
        an actual polynomial) is returned::

            sage: a.reduced_norm(var=False)
            (4, 0, 4, 1)

        TESTS:

        We check that `N` is a multiple of `a`::

            sage: S(N).is_right_divisible_by(a)
            True
            sage: S(N).is_left_divisible_by(a)
            True

        We check that the reduced norm is a multiplicative map::

            sage: a = S.random_element(degree=5)
            sage: b = S.random_element(degree=7)
            sage: a.reduced_norm() * b.reduced_norm() == (a*b).reduced_norm()
            True

        We check that the reduced norm is correctly computed for a
        constant polynomial::

            sage: c = k.random_element()
            sage: S(c).reduced_norm() == c.norm()
            True

        ALGORITHM:

        If `r` (= the order of the twist map) is small compared
        to `d` (= the degree of this skew polynomial), the reduced
        norm is computed as the determinant of the multiplication
        by `P` (= this skew polynomial) acting on `K[X,\sigma]`
        (= the underlying skew ring) viewed as a free module of
        rank `r` over `K[X^r]`.

        Otherwise, the reduced norm is computed as the characteristic
        polynomial of the left multiplication by `X` on the quotient
        `K[X,\sigma] / K[X,\sigma] P` (which is a `K`-vector space
        of dimension `d`).

        .. SEEALSO::

            :meth:`reduced_trace`, :meth:`reduced_charpoly`
        """
        if self._norm is None:
            if self.is_zero():
                self._norm = 0
            else:
                parent = self._parent
                section = parent._embed_constants.section()
                exp = (parent.base_ring().cardinality() - 1) / (parent._constants.cardinality() - 1)
                order = self.parent()._order
                lc = section(self.leading_coefficient()**exp)
                if self.degree() == 0:
                    self._norm = (lc,)
                elif order < self.degree():
                    M = self._matmul_c()
                    self._norm = tuple(lc*section(x) for x in M.determinant().monic().list())
                else:
                    charpoly = self._matphir_c().characteristic_polynomial()
                    self._norm = tuple(lc*section(x) for x in charpoly.list())
        if var is False:
            return self._norm
        center = self.parent().center(name=var)
        return center(self._norm)

    def reduced_charpoly(self, var=None):
        r"""
        Return the reduced characteristic polynomial of this
        skew polynomial.

        INPUT:

        - ``var`` -- string, a pair of strings or ``None``
          (default: ``None``); the variable names used for the
          characteristic polynomial and the center

        .. NOTE::

            The result is cached.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<u> = k['u', Frob]
            sage: a = u^3 + (2*t^2 + 3)*u^2 + (4*t^2 + t + 4)*u + 2*t^2 + 2
            sage: chi = a.reduced_charpoly()
            sage: chi
            x^3 + (2*z + 1)*x^2 + (3*z^2 + 4*z)*x + 4*z^3 + z^2 + 1

        The reduced characteristic polynomial has coefficients in the center
        of `S`, which is itself a univariate polynomial ring in the variable
        `z = u^3` over `\GF{5}`. Hence it appears as a bivariate polynomial::

            sage: chi.parent()
            Univariate Polynomial Ring in x over Univariate Polynomial Ring in z over Finite Field of size 5

        The constant coefficient of the reduced characteristic polynomial is
        the reduced norm, up to a sign::

            sage: chi[0] == -a.reduced_norm()
            True

        Its coefficient of degree `\deg(a) - 1` is the opposite of the reduced
        trace::

            sage: chi[2] == -a.reduced_trace()
            True

        By default, the name of the variable of the reduced characteristic
        polynomial is ``x`` and the name of central variable is usually ``z``
        (see :meth:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_order.center`
        for more details about this).
        The user can speciify different names if desired::

            sage: a.reduced_charpoly(var='T')  # variable name for the caracteristic polynomial
            T^3 + (2*z + 1)*T^2 + (3*z^2 + 4*z)*T + 4*z^3 + z^2 + 1

            sage: a.reduced_charpoly(var=('T', 'c'))
            T^3 + (2*c + 1)*T^2 + (3*c^2 + 4*c)*T + 4*c^3 + c^2 + 1

        .. SEEALSO::

            :meth:`reduced_trace`, :meth:`reduced_norm`
        """
        if self._charpoly is None:
            parent = self._parent
            section = parent._embed_constants.section()
            M = self._matmul_c()
            chi = M.charpoly()
            self._charpoly = [tuple(c.list()) for c in chi.list()]
            if self._norm is not None:
                self._norm = self._charpoly[-1]
        varcenter = None
        if var is None:
            varcharpoly = 'x'
        elif isinstance(var, (tuple, list)) and len(var) == 2:
            (varcharpoly, varcenter) = var
        else:
            varcharpoly = var
        center = self.parent().center(name=varcenter)
        coeffs = [center(c) for c in self._charpoly]
        return PolynomialRing(center, name=varcharpoly)(coeffs)

    def is_central(self) -> bool:
        r"""
        Return ``True`` if this skew polynomial lies in the center.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]

            sage: x.is_central()
            False
            sage: (t*x^3).is_central()
            False
            sage: (x^6 + x^3).is_central()
            True
        """
        center = self.parent()._working_center
        try:
            center(self)
            return True
        except ValueError:
            return False

    def bound(self):
        r"""
        Return a bound of this skew polynomial (i.e. a multiple
        of this skew polynomial lying in the center).

        .. NOTE::

            Since `b` is central, it divides a skew polynomial
            on the left iff it divides it on the right

        ALGORITHM:

        #. Sage first checks whether self is itself in the
           center. It if is, it returns ``self``

        #. If an optimal bound was previously computed and
           cached, Sage returns it

        #. Otherwise, Sage returns the reduced norm of ``self``

        As a consequence, the output of this function may depend
        on previous computations (an example is given below).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: Z = S.center(); Z
            Univariate Polynomial Ring in z over Finite Field of size 5

            sage: a = x^2 + (4*t + 2)*x + 4*t^2 + 3
            sage: b = a.bound(); b
            z^2 + z + 4

        We observe that the bound is explicitly given as an element of the
        center (which is a univariate polynomial ring in the variable `z`).
        We can use conversion to send it in the skew polynomial ring::

            sage: S(b)
            x^6 + x^3 + 4

        We check that `b` is divisible by `a`::

            sage: S(b).is_right_divisible_by(a)
            True
            sage: S(b).is_left_divisible_by(a)
            True

        Actually, `b` is the reduced norm of `a`::

            sage: b == a.reduced_norm()
            True

        Now, we compute the optimal bound of `a` and see that
        it affects the behaviour of ``bound()``::

            sage: a.optimal_bound()
            z + 3
            sage: a.bound()
            z + 3

        TESTS:

        We check that when the input skew polynomial lies in
        the center, the output is the skew polynomial itself::

            sage: a = Z.random_element(degree=4)
            sage: b = S(a).bound()
            sage: a == b
            True
        """
        center = self.parent().center()
        try:
            return center(self)
        except ValueError:
            pass
        if self._optbound is not None:
            return center(self._optbound)
        return self.reduced_norm()

    def optimal_bound(self):
        r"""
        Return the optimal bound of this skew polynomial (i.e.
        the monic multiple of this skew polynomial of minimal
        degree lying in the center).

        .. NOTE::

            The result is cached.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: Z = S.center(); Z
            Univariate Polynomial Ring in z over Finite Field of size 5

            sage: a = x^2 + (4*t + 2)*x + 4*t^2 + 3
            sage: b = a.optimal_bound(); b
            z + 3

        We observe that the bound is explicitly given as an element of the
        center (which is a univariate polynomial ring in the variable `z`).
        We can use conversion to send it in the skew polynomial ring::

            sage: S(b)
            x^3 + 3

        We check that `b` is divisible by `a`::

            sage: S(b).is_right_divisible_by(a)
            True
            sage: S(b).is_left_divisible_by(a)
            True
        """
        center = self.parent().center()
        if self._optbound is None:
            try:
                self._optbound = center(self).monic()
            except ValueError:
                bound = self._matphir_c().minimal_polynomial()
                section = self._parent._embed_constants.section()
                self._optbound = [section(x) for x in bound.list()]
        return center(self._optbound)

     # TODO:
     # fast multiplication
     # reduced characteristic polynomial
