"""
Base class for finite field elements

AUTHORS::

- David Roe (2010-1-14) -- factored out of sage.structure.element

"""

from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.integer import Integer

def is_FiniteFieldElement(x):
    """
    Returns if x is a finite field element.

    EXAMPLE::

        sage: from sage.rings.finite_rings.element_base import is_FiniteFieldElement
        sage: is_FiniteFieldElement(1)
        False
        sage: is_FiniteFieldElement(IntegerRing())
        False
        sage: is_FiniteFieldElement(GF(5)(2))
        True
    """
    from sage.rings.finite_rings.finite_field_base import is_FiniteField
    return isinstance(x, Element) and is_FiniteField(x.parent())

cdef class FiniteRingElement(CommutativeRingElement):
    def _nth_root_common(self, n, all, algorithm, cunningham):
        """
        This function exists to reduce code duplication between finite field
        nth roots and integer_mod nth roots.

        The inputs are described there.

        TESTS::

            sage: a = Zmod(17)(13)
            sage: a._nth_root_common(4, True, "Johnston", False)
            [3, 5, 14, 12]
            sage: a._nth_root_common(4, True, "Johnston", cunningham = True) # optional - cunningham
            [3, 5, 14, 12]
        """
        K = self.parent()
        q = K.order()
        if self.is_one():
            gcd = n.gcd(q-1)
            if gcd == 1:
                if all: return [self]
                else: return self
            else:
                # the following may eventually be improved to not need a multiplicative generator.
                g = K.multiplicative_generator()
                q1overn = (q-1)//gcd
                nthroot = g**q1overn
                return [nthroot**a for a in range(gcd)] if all else nthroot
        n = n % (q-1)
        if n == 0:
            if all: return []
            else: raise ValueError("no nth root")
        gcd, alpha, beta = n.xgcd(q-1) # gcd = alpha*n + beta*(q-1), so 1/n = alpha/gcd (mod q-1)
        if gcd == 1:
            return [self**alpha] if all else self**alpha
        n = gcd
        q1overn = (q-1)//n
        if self**q1overn != 1:
            if all: return []
            else: raise ValueError("no nth root")
        self = self**alpha
        if cunningham:
            from sage.rings.factorint import factor_cunningham
            F = factor_cunningham(n)
        else:
            F = n.factor()
        from sage.groups.generic import discrete_log
        if algorithm is None or algorithm == 'Johnston':
            g = K.multiplicative_generator()
            for r, v in F:
                k, h = (q-1).val_unit(r)
                z = h * (-h).inverse_mod(r**v)
                x = (1 + z) // r**v
                if k == 1:
                    self = self**x
                else:
                    t = discrete_log(self**h, g**(r**v*h), r**(k-v), operation='*')
                    self = self**x * g**(-z*t)
            if all:
                nthroot = g**q1overn
                L = [self]
                for i in range(1,n):
                    self *= nthroot
                    L.append(self)
                return L
            else:
                return self
        else:
            raise ValueError("unknown algorithm")

cdef class FinitePolyExtElement(FiniteRingElement):
    """
    Elements represented as polynomials modulo a given ideal.

    TESTS::

        sage: k.<a> = GF(64)
        sage: TestSuite(a).run()
    """
    def _im_gens_(self, codomain, im_gens):
        """
        Used for applying homomorphisms of finite fields.

        EXAMPLES::

            sage: k.<a> = FiniteField(73^2, 'a')
            sage: K.<b> = FiniteField(73^4, 'b')
            sage: phi = k.hom([ b^(73*73+1) ]) # indirect doctest
            sage: phi(0)
            0
            sage: phi(a)
            7*b^3 + 13*b^2 + 65*b + 71

            sage: phi(a+3)
            7*b^3 + 13*b^2 + 65*b + 1
        """
        ## NOTE: see the note in sage/rings/number_field_element.pyx,
        ## in the comments for _im_gens_ there -- something analogous
        ## applies here.
        return codomain(self.polynomial()(im_gens[0]))

    def minpoly(self,var='x'):
        """
        Returns the minimal polynomial of this element
        (over the corresponding prime subfield).

        EXAMPLES::

            sage: k.<a> = FiniteField(19^2)
            sage: parent(a)
            Finite Field in a of size 19^2
            sage: b=a**20;p=b.charpoly("x");p
            x^2 + 15*x + 4
            sage: factor(p)
            (x + 17)^2
            sage: b.minpoly('x')
            x + 17
        """
        p=self.charpoly(var);
        for q in p.factor():
            if q[0](self)==0:
                return q[0]
        # This shouldn't be reached, but you never know!
        raise ArithmeticError("Could not find the minimal polynomial")

        ## We have two names for the same method
        ## for compatibility with sage.matrix
    def minimal_polynomial(self,var='x'):
        """
        Returns the minimal polynomial of this element
        (over the corresponding prime subfield).

        EXAMPLES::

            sage: k.<a> = FiniteField(3^4)
            sage: parent(a)
            Finite Field in a of size 3^4
            sage: b=a**20;p=charpoly(b,"y");p
            y^4 + 2*y^2 + 1
            sage: factor(p)
            (y^2 + 1)^2
            sage: b.minimal_polynomial('y')
            y^2 + 1
        """
        return self.minpoly(var)

    def _vector_(self, reverse=False):
        """
        Return a vector in self.parent().vector_space() matching
        self. The most significant bit is to the right.

        INPUT:

        - ``reverse`` -- reverse the order of the bits
          from little endian to big endian.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: e = a^2 + 1
            sage: v = vector(e)
            sage: v
            (1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: k(v)
            a^2 + 1

            sage: k.<a> = GF(3^16)
            sage: e = 2*a^2 + 1
            sage: v = vector(e)
            sage: v
            (1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: k(v)
            2*a^2 + 1

        You can also compute the vector in the other order::

            sage: e._vector_(reverse=True)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1)
        """
        #vector(foo) might pass in ZZ
        if isinstance(reverse, Parent):
            raise TypeError("Base field is fixed to prime subfield.")

        k = self.parent()

        v = self.polynomial().list()

        ret = [v[i] for i in range(len(v))]

        for i in range(k.degree() - len(ret)):
            ret.append(0)

        if reverse:
            ret = list(reversed(ret))
        return k.vector_space()(ret)

    def _matrix_(self, reverse=False):
        """
        Return the matrix of right multiplication by the element on
        the power basis `1, x, x^2, \ldots, x^{d-1}` for the field
        extension.  Thus the \emph{rows} of this matrix give the images
        of each of the `x^i`.

        INPUT:

        - ``reverse`` - if True act on vectors in reversed order

        EXAMPLE::

            sage: k.<a> = GF(2^4)
            sage: a._vector_(reverse=True), a._matrix_(reverse=True) * a._vector_(reverse=True)
            ((0, 0, 1, 0), (0, 1, 0, 0))
            sage: vector(a), matrix(a) * vector(a)
            ((0, 1, 0, 0), (0, 0, 1, 0))
        """
        import sage.matrix.matrix_space

        K = self.parent()
        a = K.gen()
        x = K(1)
        d = K.degree()

        columns = []

        if not reverse:
            l = xrange(d)
        else:
            l = reversed(range(d))

        for i in l:
            columns.append( (self * x)._vector_() )
            x *= a

        k = K.base_ring()
        M = sage.matrix.matrix_space.MatrixSpace(k, d)

        if reverse:
            return M(columns)
        else:
            return M(columns).transpose()
    def _latex_(self):
        r"""
        Return the latex representation of self, which is just the
        latex representation of the polynomial representation of self.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: b._latex_()
            'b'
            sage: (b^2+1)._latex_()
            'b + 4'
        """
        if self.parent().degree()>1:
            return self.polynomial()._latex_()
        else:
            return str(self)

    def _pari_(self, var=None):
        r"""
        Return PARI representation of this finite field element.

        INPUT:

        - ``var`` -- (default: ``None``) optional variable string

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: a._pari_()
            a
            sage: a._pari_('b')
            b
            sage: t = 3*a^2 + 2*a + 4
            sage: t_string = t._pari_init_('y')
            sage: t_string
            'Mod(Mod(3, 5)*y^2 + Mod(2, 5)*y + Mod(4, 5), Mod(1, 5)*y^3 + Mod(3, 5)*y + Mod(3, 5))'
            sage: type(t_string)
            <type 'str'>
            sage: t_element = t._pari_('b')
            sage: t_element
            3*b^2 + 2*b + 4
            sage: t_element.parent()
            Interface to the PARI C library
        """
        if var is None:
            var = self.parent().variable_name()
        from sage.libs.pari.all import pari
        ffgen = self._parent.modulus()._pari_with_name(var).ffgen()
        polypari = self.polynomial()._pari_with_name()
        # Add ffgen - ffgen to ensure that we really get an FFELT
        return polypari.subst("x", ffgen) + ffgen - ffgen

    def _pari_init_(self, var=None):
        r"""
        Return a string that defines this element when evaluated in PARI.

        INPUT:

        - ``var`` - default: ``None`` - a string for a new variable name to use.

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b._pari_init_()
            'Mod(Mod(1, 5)*b, Mod(1, 5)*b^2 + Mod(4, 5)*b + Mod(2, 5))'
            sage: (2*b+3)._pari_init_()
            'Mod(Mod(2, 5)*b + Mod(3, 5), Mod(1, 5)*b^2 + Mod(4, 5)*b + Mod(2, 5))'

        TESTS:

        The following tests against a bug fixed in :trac:`11530`::

            sage: F.<d> = GF(3^4)
            sage: F.modulus()
            x^4 + 2*x^3 + 2
            sage: d._pari_init_()
            'Mod(Mod(1, 3)*d, Mod(1, 3)*d^4 + Mod(2, 3)*d^3 + Mod(2, 3))'
            sage: (d^2+2*d+1)._pari_init_("p")
            'Mod(Mod(1, 3)*p^2 + Mod(2, 3)*p + Mod(1, 3), Mod(1, 3)*p^4 + Mod(2, 3)*p^3 + Mod(2, 3))'
            sage: d._pari_()
            d

            sage: K.<M> = GF(2^8)
            sage: K.modulus()
            x^8 + x^4 + x^3 + x^2 + 1
            sage: (M^3+1)._pari_init_()
            'Mod(Mod(1, 2)*M^3 + Mod(1, 2), Mod(1, 2)*M^8 + Mod(1, 2)*M^4 + Mod(1, 2)*M^3 + Mod(1, 2)*M^2 + Mod(1, 2))'
            sage: M._pari_init_(var='foo')
            'Mod(Mod(1, 2)*foo, Mod(1, 2)*foo^8 + Mod(1, 2)*foo^4 + Mod(1, 2)*foo^3 + Mod(1, 2)*foo^2 + Mod(1, 2))'
        """
        if var is None:
            var = self.parent().variable_name()
        g = self.parent().modulus()._pari_with_name(var)
        f = self.polynomial()._pari_with_name(var)
        return 'Mod({0}, {1})'.format(f, g)

    def charpoly(self, var='x', algorithm='matrix'):
        """
        Return the characteristic polynomial of self as a polynomial with given variable.

        INPUT:

        - ``var`` - string (default: 'x')

        - ``algorithm`` - string (default: 'matrix')

          - 'matrix' - return the charpoly computed from the matrix of
            left multiplication by self

          - 'pari' -- use pari's charpoly routine on polymods, which
            is not very good except in small cases

        The result is not cached.

        EXAMPLES::

            sage: k.<a> = GF(19^2)
            sage: parent(a)
            Finite Field in a of size 19^2
            sage: a.charpoly('X')
            X^2 + 18*X + 2
            sage: a^2 + 18*a + 2
            0
            sage: a.charpoly('X', algorithm='pari')
            X^2 + 18*X + 2
        """
        if algorithm == 'matrix':
            return self._matrix_().charpoly(var)
        elif algorithm == 'pari':
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self.parent().prime_subfield(), var)
            return R(self._pari_().charpoly('x').lift())
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)

    def norm(self):
        """
        Return the norm of self down to the prime subfield.

        This is the product of the Galois conjugates of self.

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b.norm()
            2
            sage: b.charpoly('t')
            t^2 + 4*t + 2

        Next we consider a cubic extension::

            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: a.norm()
            2
            sage: a.charpoly('t')
            t^3 + 3*t + 3
            sage: a * a^5 * (a^25)
            2
        """
        f = self.charpoly('x')
        n = f[0]
        if f.degree() % 2 != 0:
            return -n
        else:
            return n

    def trace(self):
        """
        Return the trace of this element, which is the sum of the
        Galois conjugates.

        EXAMPLES::

            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: a.trace()
            0
            sage: a.charpoly('t')
            t^3 + 3*t + 3
            sage: a + a^5 + a^25
            0
            sage: z = a^2 + a + 1
            sage: z.trace()
            2
            sage: z.charpoly('t')
            t^3 + 3*t^2 + 2*t + 2
            sage: z + z^5 + z^25
            2
        """
        return self.parent().prime_subfield()(self._pari_().trace().lift())

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of this field element.

        EXAMPLE::

            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: a.multiplicative_order()
            124
            sage: (a^8).multiplicative_order()
            31
            sage: S(0).multiplicative_order()
            Traceback (most recent call last):
            ...
            ArithmeticError: Multiplicative order of 0 not defined.
        """
        if self.is_zero():
            raise ArithmeticError("Multiplicative order of 0 not defined.")
        n = self._parent.order() - 1
        F = self._parent.factored_unit_order()[0]
        order = Integer(1)
        for p, e in F:
            # Determine the power of p that divides the order.
            a = self**(n//(p**e))
            while a != 1:
                order = order * p
                a = a**p

        return order

    def additive_order(self):
        """
        Return the additive order of this finite field element.

        EXAMPLES::

            sage: k.<a> = FiniteField(2^12, 'a')
            sage: b = a^3 + a + 1
            sage: b.additive_order()
            2
            sage: k(0).additive_order()
            1
        """
        if self.is_zero():
            from sage.rings.integer import Integer
            return Integer(1)
        return self.parent().characteristic()

    def is_square(self):
        """
        Returns ``True`` if and only if this element is a perfect square.

        EXAMPLES::

            sage: k.<a> = FiniteField(9, impl='givaro', modulus='primitive')
            sage: a.is_square()
            False
            sage: (a**2).is_square()
            True
            sage: k.<a> = FiniteField(4, impl='ntl', modulus='primitive')
            sage: (a**2).is_square()
            True
            sage: k.<a> = FiniteField(17^5, impl='pari_ffelt', modulus='primitive')
            sage: a.is_square()
            False
            sage: (a**2).is_square()
            True

        ::

            sage: k(0).is_square()
            True
        """
        K = self.parent()
        if K.characteristic() == 2:
            return True
        n = K.order() - 1
        a = self**(n // 2)
        return a == 1 or a == 0

    def square_root(self, extend=False, all=False):
        """
        The square root function.

        INPUT:


        -  ``extend`` -- bool (default: ``True``); if ``True``, return a
           square root in an extension ring, if necessary. Otherwise, raise a
           ValueError if the root is not in the base ring.

           .. WARNING::

               This option is not implemented!

        -  ``all`` - bool (default: ``False``); if ``True``, return all
           square roots of ``self``, instead of just one.

        .. WARNING::

           The ``'extend'`` option is not implemented (yet).

        EXAMPLES::

            sage: F = FiniteField(7^2, 'a')
            sage: F(2).square_root()
            4
            sage: F(3).square_root()
            2*a + 6
            sage: F(3).square_root()**2
            3
            sage: F(4).square_root()
            2
            sage: K = FiniteField(7^3, 'alpha', impl='pari_ffelt')
            sage: K(3).square_root()
            Traceback (most recent call last):
            ...
            ValueError: must be a perfect square.
        """
        try:
            return self.nth_root(2, extend=extend, all=all)
        except ValueError:
            raise ValueError("must be a perfect square.")

    def sqrt(self, extend=False, all = False):
        """
        See :meth:square_root().

        EXAMPLES::

            sage: k.<a> = GF(3^17)
            sage: (a^3 - a - 1).sqrt()
            a^16 + 2*a^15 + a^13 + 2*a^12 + a^10 + 2*a^9 + 2*a^8 + a^7 + a^6 + 2*a^5 + a^4 + 2*a^2 + 2*a + 2
        """
        return self.square_root(extend=extend, all=all)

    def nth_root(self, n, extend = False, all = False, algorithm=None, cunningham=False):
        r"""
        Returns an `n`\th root of ``self``.

        INPUT:

        - ``n`` - integer `\geq 1`

        - ``extend`` - bool (default: ``False``); if ``True``, return an `n`\th
          root in an extension ring, if necessary. Otherwise, raise a
          ValueError if the root is not in the base ring.  Warning:
          this option is not implemented!

        - ``all`` - bool (default: ``False``); if ``True``, return all `n`\th
          roots of ``self``, instead of just one.

        - ``algorithm`` - string (default: ``None``); 'Johnston' is the only
          currently supported option.  For IntegerMod elements, the problem
          is reduced to the prime modulus case using CRT and `p`-adic logs,
          and then this algorithm used.

        OUTPUT:

        If self has an `n`\th root, returns one (if ``all`` is ``False``) or a
        list of all of them (if ``all`` is ``True``).
        Otherwise, raises a ``ValueError`` (if ``extend`` is ``False``)
        or a ``NotImplementedError`` (if ``extend`` is ``True``).

        .. warning::

           The ``extend`` option is not implemented (yet).

        EXAMPLES::

            sage: K = GF(31)
            sage: a = K(22)
            sage: K(22).nth_root(7)
            13
            sage: K(25).nth_root(5)
            5
            sage: K(23).nth_root(3)
            29

            sage: K.<a> = GF(625)
            sage: (3*a^2+a+1).nth_root(13)**13
            3*a^2 + a + 1

            sage: k.<a> = GF(29^2)
            sage: b = a^2 + 5*a + 1
            sage: b.nth_root(11)
            3*a + 20
            sage: b.nth_root(5)
            Traceback (most recent call last):
            ...
            ValueError: no nth root
            sage: b.nth_root(5, all = True)
            []
            sage: b.nth_root(3, all = True)
            [14*a + 18, 10*a + 13, 5*a + 27]

            sage: k.<a> = GF(29^5)
            sage: b = a^2 + 5*a + 1
            sage: b.nth_root(5)
            19*a^4 + 2*a^3 + 2*a^2 + 16*a + 3
            sage: b.nth_root(7)
            Traceback (most recent call last):
            ...
            ValueError: no nth root
            sage: b.nth_root(4, all=True)
            []

        TESTS::

            sage: for p in [2,3,5,7,11]:  # long time, random because of PARI warnings
            ....:     for n in [2,5,10]:
            ....:         q = p^n
            ....:         K.<a> = GF(q)
            ....:         for r in (q-1).divisors():
            ....:             if r == 1: continue
            ....:             x = K.random_element()
            ....:             y = x^r
            ....:             assert y.nth_root(r)^r == y
            ....:             assert (y^41).nth_root(41*r)^(41*r) == y^41
            ....:             assert (y^307).nth_root(307*r)^(307*r) == y^307
            sage: k.<a> = GF(4)
            sage: a.nth_root(0,all=True)
            []
            sage: k(1).nth_root(0,all=True)
            [a, a + 1, 1]

        ALGORITHMS:

        - The default is currently an algorithm described in the following paper:

        Johnston, Anna M. A generalized qth root algorithm. Proceedings of the tenth annual ACM-SIAM symposium on Discrete algorithms. Baltimore, 1999: pp 929-930.

        AUTHOR:

        - David Roe (2010-02-13)
        """
        if self.is_zero():
            if n <= 0:
                if all: return []
                else: raise ValueError
            if all: return [self]
            else: return self
        if n < 0:
            self = ~self
            n = -n
        elif n == 0:
            if self == 1:
                if all: return [a for a in self.parent().list() if a != 0]
                else: return self
            else:
                if all: return []
                else: raise ValueError
        if extend:
            raise NotImplementedError
        from sage.rings.integer import Integer
        n = Integer(n)
        return self._nth_root_common(n, all, algorithm, cunningham)

    def pth_power(self, int k = 1):
        """
        Return the `(p^k)^{th}` power of self, where `p` is the
        characteristic of the field.

        INPUT:

        - ``k`` - integer (default: 1, must fit in C int type)

        Note that if `k` is negative, then this computes the appropriate root.

        EXAMPLES::

            sage: F.<a> = GF(29^2)
            sage: z = a^2 + 5*a + 1
            sage: z.pth_power()
            19*a + 20
            sage: z.pth_power(10)
            10*a + 28
            sage: z.pth_power(-10) == z
            True
            sage: F.<b> = GF(2^12)
            sage: y = b^3 + b + 1
            sage: y == (y.pth_power(-3))^(2^3)
            True
            sage: y.pth_power(2)
            b^7 + b^6 + b^5 + b^4 + b^3 + b
        """
        p = self.additive_order()
        n = self.parent().degree()
        return self**(p**(k % n))

    frobenius = pth_power

    def pth_root(self, int k = 1):
        """
        Return the `(p^k)^{th}` root of self, where `p` is the characteristic
        of the field.

        INPUT:

        - ``k`` - integer (default: 1, must fit in C int type)

        Note that if `k` is negative, then this computes the appropriate power.

        EXAMPLES::

            sage: F.<b> = GF(2^12)
            sage: y = b^3 + b + 1
            sage: y == (y.pth_root(3))^(2^3)
            True
            sage: y.pth_root(2)
            b^11 + b^10 + b^9 + b^7 + b^5 + b^4 + b^2 + b
        """
        return self.pth_power(-k)

