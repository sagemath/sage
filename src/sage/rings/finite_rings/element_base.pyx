# sage.doctest: needs sage.rings.finite_rings
"""
Base class for finite field elements

AUTHORS:

- David Roe (2010-01-14): factored out of sage.structure.element
- Sebastian Oehms (2018-07-19): added :meth:`conjugate` (see :issue:`26761`)
"""

# ****************************************************************************
#       Copyright (C) 2010 David Roe <roed@math.harvard.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer


def is_FiniteFieldElement(x):
    """
    Return ``True`` if ``x`` is a finite field element.

    This function is deprecated.

    EXAMPLES::

        sage: from sage.rings.finite_rings.element_base import is_FiniteFieldElement
        sage: is_FiniteFieldElement(1)
        doctest:...: DeprecationWarning: the function is_FiniteFieldElement is deprecated; use isinstance(x, sage.structure.element.FieldElement) and x.parent().is_finite() instead
        See https://github.com/sagemath/sage/issues/32664 for details.
        False
        sage: is_FiniteFieldElement(IntegerRing())
        False
        sage: is_FiniteFieldElement(GF(5)(2))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(32664, "the function is_FiniteFieldElement is deprecated; use isinstance(x, sage.structure.element.FieldElement) and x.parent().is_finite() instead")

    from sage.rings.finite_rings.finite_field_base import FiniteField
    return isinstance(x, Element) and isinstance(x.parent(), FiniteField)


cdef class FiniteRingElement(CommutativeRingElement):
    def _nth_root_common(self, n, all, algorithm, cunningham):
        """
        This function exists to reduce code duplication between finite field
        `n`-th roots and ``integer_mod`` `n`-th roots. It assumes that ``self``
        is a field element.

        The inputs are described there.

        TESTS::

            sage: a = Zmod(17)(13)
            sage: sorted(a._nth_root_common(4, True, "Johnston", False))
            [3, 5, 12, 14]
            sage: sorted(a._nth_root_common(4, True, "Johnston", cunningham=True))  # optional - cunningham_tables
            [3, 5, 12, 14]

        Test various prime powers::

            sage: p = 5^5*10000000100 + 1
            sage: a = GF(p)(3)**(5^7)
            sage: for e in range(20):
            ....:     r = a._nth_root_common(5^e, False, "Johnston", False)
            ....:     assert r**(5^e) == a

        Test very large modulus (assumed impossible to factor in reasonable time)::

            sage: p = 2^1024 + 643
            sage: a = GF(p, proof=False)(3)**(29*283*3539)
            sage: r = a._nth_root_common(29*283*3539*12345, False, "Johnston", False)
            sage: r**(29*283*3539*12345) == a
            True
        """
        K = self.parent()
        q = K.order()
        gcd = n.gcd(q-1)
        if self.is_one():
            if gcd == 1:
                return [self] if all else self
            nthroot = K.zeta(gcd)
            return [nthroot**a for a in range(gcd)] if all else nthroot
        if gcd == q-1:
            if all:
                return []
            raise ValueError("no nth root")
        gcd, alpha, _ = n.xgcd(q-1)  # gcd = alpha*n + beta*(q-1), so 1/n = alpha/gcd (mod q-1)
        if gcd == 1:
            return [self**alpha] if all else self**alpha

        n = gcd
        q1overn = (q-1)//n
        if self**q1overn != 1:
            if all:
                return []
            raise ValueError("no nth root")
        self = self**alpha
        if cunningham:
            from sage.rings.factorint import factor_cunningham
            F = factor_cunningham(n)
        else:
            F = n.factor()
        from sage.groups.generic import discrete_log
        if algorithm is None or algorithm == 'Johnston':
            # In the style of the Adleman-Manders-Miller algorithm,
            # we will use small order elements instead of a multiplicative
            # generator, which can be expensive to compute.
            for r, v in F:
                # 0 < v <= k
                k, h = (q-1).val_unit(r)
                hinv = (-h).inverse_mod(r**v)
                z = h * hinv
                x = (1 + z) // r**v
                if k == v:
                    self = self**x
                else:
                    # We need an element of order r^k (g^h in Johnston's article)
                    # self^x differs from the actual nth root by an element of
                    # order dividing r^(k-v)
                    gh = K.zeta(r**k)
                    t = discrete_log(self**h, gh**(r**v), r**(k-v), operation='*')
                    self = self**x * gh**(-hinv*t)
            if all:
                nthroot = K.zeta(n)
                L = [self]
                for i in range(1, n):
                    self *= nthroot
                    L.append(self)
                return L
            return self
        else:
            raise ValueError("unknown algorithm")

    def to_bytes(self, byteorder='big'):
        r"""
        Return an array of bytes representing an integer.

        Internally relies on the python ``int.to_bytes()`` method.
        Length of byte array is determined from the field's order.

        INPUT:

        - ``byteorder`` -- string (default: ``'big'``); determines the byte order of
          ``input_bytes``; can only be ``'big'`` or ``'little'``

        EXAMPLES::

            sage: F = GF(65537)
            sage: a = F(8726)
            sage: a.to_bytes()
            b'\x00"\x16'
            sage: a.to_bytes(byteorder='little')
            b'\x16"\x00'
        """
        length = (self.parent().order().nbits() + 7) // 8
        return int(self).to_bytes(length=length, byteorder=byteorder)

cdef class FinitePolyExtElement(FiniteRingElement):
    """
    Elements represented as polynomials modulo a given ideal.

    TESTS::

        sage: k.<a> = GF(64)
        sage: TestSuite(a).run()
    """
    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Used for applying homomorphisms of finite fields.

        EXAMPLES::

            sage: k.<a> = FiniteField(73^2)
            sage: K.<b> = FiniteField(73^4)
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
        f = self.polynomial()
        if base_map is not None:
            Cx = codomain['x']
            f = Cx([base_map(c) for c in f])
        return codomain(f(im_gens[0]))

    def minpoly(self, var='x', algorithm='pari'):
        """
        Return the minimal polynomial of this element
        (over the corresponding prime subfield).

        INPUT:

        - ``var`` -- string (default: ``'x'``)

        - ``algorithm`` -- string (default: ``'pari'``):

          - ``'pari'`` -- use pari's minpoly

          - ``'matrix'`` -- return the minpoly computed from the matrix of
            left multiplication by self

        EXAMPLES::

            sage: from sage.rings.finite_rings.element_base import FinitePolyExtElement
            sage: k.<a> = FiniteField(19^2)
            sage: parent(a)
            Finite Field in a of size 19^2
            sage: b=a**20
            sage: p=FinitePolyExtElement.minpoly(b,"x", algorithm='pari')
            sage: q=FinitePolyExtElement.minpoly(b,"x", algorithm='matrix')
            sage: q == p
            True
            sage: p
            x + 17
        """
        if self.polynomial().degree() == 0:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self.parent().prime_subfield(), var)
            return R.gen() - self.polynomial()[0]

        if algorithm == 'pari':
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self.parent().prime_subfield(), var)
            return R(self.__pari__().minpoly('x').lift())
        
        if algorithm == 'matrix':
            return self.matrix().minpoly(var)

        raise ValueError("unknown algorithm '%s'" % algorithm)

    # We have two names for the same method
    # for compatibility with sage.matrix
    def minimal_polynomial(self, var='x'):
        """
        Return the minimal polynomial of this element
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

    def __getitem__(self, n):
        r"""
        Return the `n`-th coefficient of this finite field element when
        written as a polynomial in the generator.

        EXAMPLES::

            sage: x = polygen(GF(19))
            sage: F.<i> = GF(19^2, modulus=x^2+1)
            sage: a = 5 + 7*i
            sage: a[0]
            5
            sage: a[1]
            7

        ::

            sage: b = F(11)
            sage: b[0]
            11
            sage: b[1]
            0

        TESTS::

            sage: # needs sage.modules
            sage: F,t = GF(random_prime(99)^randrange(2,99), 't').objgen()
            sage: a = F.random_element()
            sage: all(a[i] == a.polynomial()[i] for i in range(F.degree()))
            True
            sage: a == sum(a[i]*t^i for i in range(F.degree()))
            True
        """
        if n < 0 or n >= self.parent().degree():
            raise IndexError("index must lie between 0 and the degree minus 1")
        return self.polynomial()[n]

    def list(self):
        r"""
        Return the list of coefficients (in little-endian) of this
        finite field element when written as a polynomial in the
        generator.

        Equivalent to calling ``list()`` on this element.

        EXAMPLES::

            sage: x = polygen(GF(71))
            sage: F.<u> = GF(71^7, modulus=x^7 + x + 1)
            sage: a = 3 + u + 3*u^2 + 3*u^3 + 7*u^4
            sage: a.list()
            [3, 1, 3, 3, 7, 0, 0]
            sage: a.list() == list(a) == [a[i] for i in range(F.degree())]
            True

        The coefficients returned are those of a fully reduced
        representative of the finite field element::

            sage: b = u^777
            sage: b.list()
            [9, 69, 4, 27, 40, 10, 56]
            sage: (u.polynomial()^777).list()
            [0, 0, 0, 0, ..., 0, 1]

        TESTS::

            sage: # needs sage.modules
            sage: R.<x> = GF(17)[]
            sage: F.<t> = GF(17^60)
            sage: a = F.random_element()
            sage: a == R(a.list())(t)
            True
            sage: list(a) == a.list()
            True
        """
        return self.polynomial().padded_list(self.parent().degree())

    def __iter__(self):
        r"""
        Return an iterator over the coefficients of this finite field
        element, in the same order as :meth:`list`.

        EXAMPLES::

            sage: x = polygen(GF(19))
            sage: F.<i> = GF(19^2, modulus=x^2+1)
            sage: a = 5 + 7*i
            sage: it = iter(a)
            sage: next(it)
            5
            sage: next(it)
            7
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
            sage: list(a)   # implicit doctest
            [5, 7]
            sage: tuple(a)  # implicit doctest
            (5, 7)
            sage: b = F(11)
            sage: list(b)   # implicit doctest
            [11, 0]
            sage: tuple(b)  # implicit doctest
            (11, 0)
            sage: list(b.polynomial())
            [11]

        TESTS::

            sage: # needs sage.modules
            sage: F = GF(random_prime(333)^randrange(111,999),'t')
            sage: a = F.random_element()
            sage: list(a) == a.list()  # implicit doctest
            True

        ::

            sage: # needs sage.modules
            sage: F.<t> = GF(17^60)
            sage: a = F.random_element()
            sage: a == sum(c*t^i for i,c in enumerate(a))  # implicit doctest
            True

        ::

            sage: # needs sage.modules
            sage: F.<t> = GF((2^127 - 1)^10, 't')
            sage: a = F.random_element()
            sage: a == sum(c*t^i for i,c in enumerate(a))  # implicit doctest
            True
        """
        return iter(self.list())

    def _vector_(self, reverse=False):
        """
        Return a vector matching this element in the vector space attached
        to the parent.  The most significant coefficient is to the right.

        INPUT:

        - ``reverse`` -- reverse the order of the bits
          from little endian to big endian

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
        # vector(foo) might pass in ZZ
        if isinstance(reverse, Parent):
            raise TypeError("Base field is fixed to prime subfield.")

        k = self.parent()
        ret = self.polynomial().padded_list(k.degree())

        if reverse:
            ret.reverse()
        return k.vector_space(map=False)(ret)

    def matrix(self, reverse=False):
        r"""
        Return the matrix of left multiplication by the element on
        the power basis `1, x, x^2, \ldots, x^{d-1}` for the field
        extension.

        Thus the \emph{columns} of this matrix give the images
        of each of the `x^i`.

        INPUT:

        - ``reverse`` -- if ``True``, act on vectors in reversed order

        EXAMPLES::

            sage: # needs sage.modules
            sage: k.<a> = GF(2^4)
            sage: b = k.random_element()
            sage: vector(a*b) == a.matrix() * vector(b)
            True
            sage: (a*b)._vector_(reverse=True) == a.matrix(reverse=True) * b._vector_(reverse=True)
            True
        """
        K = self.parent()
        a = K.gen()
        x = K(1)
        d = K.degree()

        columns = []

        for i in range(d):
            columns.append( (self * x)._vector_(reverse=reverse) )
            x *= a

        if reverse:
            columns.reverse()

        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(K.base_ring(), d)

        return M(columns).transpose()

    def _latex_(self):
        r"""
        Return the latex representation of ``self``, which is just the
        latex representation of the polynomial representation of ``self``.

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
        return str(self)

    def __pari__(self, var=None):
        r"""
        Return PARI representation of this finite field element.

        INPUT:

        - ``var`` -- (default: ``None``) optional variable string

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: a.__pari__()
            a
            sage: a.__pari__('b')
            b
            sage: t = 3*a^2 + 2*a + 4
            sage: t_string = t._pari_init_('y')
            sage: t_string
            'Mod(Mod(3, 5)*y^2 + Mod(2, 5)*y + Mod(4, 5), Mod(1, 5)*y^3 + Mod(3, 5)*y + Mod(3, 5))'
            sage: type(t_string)
            <... 'str'>
            sage: t_element = t.__pari__('b')
            sage: t_element
            3*b^2 + 2*b + 4
            sage: type(t_element)
            <class 'cypari2.gen.Gen'>
        """
        if var is None:
            var = self.parent().variable_name()
        ffgen = self._parent.modulus()._pari_with_name(var).ffgen()
        polypari = self.polynomial()._pari_with_name()
        # Add ffgen - ffgen to ensure that we really get an FFELT
        return polypari.subst("x", ffgen) + ffgen - ffgen

    def _pari_init_(self, var=None):
        r"""
        Return a string that defines this element when evaluated in PARI.

        INPUT:

        - ``var`` -- (default: ``None``) a string for a new variable name to use

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b._pari_init_()
            'Mod(Mod(1, 5)*b, Mod(1, 5)*b^2 + Mod(4, 5)*b + Mod(2, 5))'
            sage: (2*b+3)._pari_init_()
            'Mod(Mod(2, 5)*b + Mod(3, 5), Mod(1, 5)*b^2 + Mod(4, 5)*b + Mod(2, 5))'

        TESTS:

        The following tests against a bug fixed in :issue:`11530`::

            sage: F.<d> = GF(3^4)
            sage: F.modulus()
            x^4 + 2*x^3 + 2
            sage: d._pari_init_()
            'Mod(Mod(1, 3)*d, Mod(1, 3)*d^4 + Mod(2, 3)*d^3 + Mod(2, 3))'
            sage: (d^2+2*d+1)._pari_init_("p")
            'Mod(Mod(1, 3)*p^2 + Mod(2, 3)*p + Mod(1, 3), Mod(1, 3)*p^4 + Mod(2, 3)*p^3 + Mod(2, 3))'
            sage: d.__pari__()
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

    def charpoly(self, var='x', algorithm='pari'):
        """
        Return the characteristic polynomial of ``self`` as a polynomial with given variable.

        INPUT:

        - ``var`` -- string (default: ``'x'``)

        - ``algorithm`` -- string (default: ``'pari'``):

          - ``'pari'`` -- use pari's charpoly

          - ``'matrix'`` -- return the charpoly computed from the matrix of
            left multiplication by ``self``

        The result is not cached.

        EXAMPLES::

            sage: from sage.rings.finite_rings.element_base import FinitePolyExtElement
            sage: k.<a> = FiniteField(19^2)
            sage: parent(a)
            Finite Field in a of size 19^2
            sage: b = a**20
            sage: p = FinitePolyExtElement.charpoly(b, "x", algorithm='pari')
            sage: q = FinitePolyExtElement.charpoly(b, "x", algorithm='matrix')         # needs sage.modules
            sage: q == p                                                                # needs sage.modules
            True
            sage: p
            x^2 + 15*x + 4
            sage: factor(p)
            (x + 17)^2
            sage: b.minpoly('x')
            x + 17
        """
        if algorithm == 'pari':
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self.parent().prime_subfield(), var)
            return R(self.__pari__().charpoly('x').lift())

        if algorithm == 'matrix':
            return self.matrix().charpoly(var)

        raise ValueError("unknown algorithm '%s'" % algorithm)

    def norm(self):
        """
        Return the norm of ``self`` down to the prime subfield.

        This is the product of the Galois conjugates of ``self``.

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
        return -n if f.degree() % 2 else n

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
        return self.parent().prime_subfield()(self.__pari__().trace().lift())

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of this field element.

        EXAMPLES::

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
            return Integer(1)
        return self.parent().characteristic()

    def is_square(self):
        """
        Return ``True`` if and only if this element is a perfect square.

        EXAMPLES::

            sage: k.<a> = FiniteField(9, impl='givaro', modulus='primitive')            # needs sage.libs.linbox
            sage: a.is_square()                                                         # needs sage.libs.linbox
            False
            sage: (a**2).is_square()                                                    # needs sage.libs.linbox
            True
            sage: k.<a> = FiniteField(4, impl='ntl', modulus='primitive')               # needs sage.libs.ntl
            sage: (a**2).is_square()                                                    # needs sage.libs.ntl
            True
            sage: k.<a> = FiniteField(17^5, impl='pari_ffelt', modulus='primitive')     # needs sage.libs.pari
            sage: a.is_square()                                                         # needs sage.libs.pari
            False
            sage: (a**2).is_square()                                                    # needs sage.libs.pari
            True

        ::

            sage: k(0).is_square()                                                      # needs sage.libs.linbox
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

        - ``extend`` -- boolean (default: ``True``); if ``True``, return a
          square root in an extension ring, if necessary. Otherwise, raise a
          :exc:`ValueError` if the root is not in the base ring.

           .. WARNING::

               This option is not implemented!

        - ``all`` -- boolean (default: ``False``); if ``True``, return all
          square roots of ``self``, instead of just one

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

    def sqrt(self, extend=False, all=False):
        """
        See :meth:`square_root`.

        EXAMPLES::

            sage: k.<a> = GF(3^17)
            sage: (a^3 - a - 1).sqrt()
            a^16 + 2*a^15 + a^13 + 2*a^12 + a^10 + 2*a^9 + 2*a^8 + a^7 + a^6 + 2*a^5 + a^4 + 2*a^2 + 2*a + 2
        """
        return self.square_root(extend=extend, all=all)

    def nth_root(self, n, extend=False, all=False, algorithm=None, cunningham=False):
        r"""
        Return an `n`-th root of ``self``.

        INPUT:

        - ``n`` -- integer `\geq 1`

        - ``extend`` -- boolean (default: ``False``); if ``True``, return an
          `n`-th root in an extension ring, if necessary. Otherwise, raise a
          :exc:`ValueError` if the root is not in the base ring.  Warning:
          this option is not implemented!

        - ``all`` -- boolean (default: ``False``); if ``True``, return all `n`-th
          roots of ``self``, instead of just one

        - ``algorithm`` -- string (default: ``None``); ``'Johnston'`` is the
          only currently supported option.  For IntegerMod elements, the problem
          is reduced to the prime modulus case using CRT and `p`-adic logs,
          and then this algorithm used.

        OUTPUT:

        If ``self`` has an `n`-th root, returns one (if ``all`` is ``False``) or a
        list of all of them (if ``all`` is ``True``).
        Otherwise, raises a :exc:`ValueError` (if ``extend`` is ``False``)
        or a :exc:`NotImplementedError` (if ``extend`` is ``True``).

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

        ALGORITHM:

        The default is currently an algorithm described in [Joh1999]_.

        AUTHOR:

        - David Roe (2010-02-13)
        """
        if self.is_zero():
            if n <= 0:
                if all:
                    return []
                raise ValueError
            return [self] if all else self
        if n < 0:
            self = ~self
            n = -n
        elif n == 0:
            if self == 1:
                if all:
                    return [a for a in self.parent().list() if a != 0]
                return self
            else:
                if all:
                    return []
                raise ValueError
        if extend:
            raise NotImplementedError
        n = Integer(n)
        return self._nth_root_common(n, all, algorithm, cunningham)

    def pth_power(self, int k=1):
        """
        Return the `(p^k)`-th power of self, where `p` is the
        characteristic of the field.

        INPUT:

        - ``k`` -- integer (default: 1, must fit in C int type)

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

    def pth_root(self, int k=1):
        """
        Return the `(p^k)`-th root of self, where `p` is the characteristic
        of the field.

        INPUT:

        - ``k`` -- integer (default: 1, must fit in C int type)

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

    def conjugate(self):
        """
        This methods returns the result of the Frobenius morphism
        in the case where the field is a quadratic extension, say
        `GF(q^2)`, where `q=p^k` is a prime power and `p` the
        characteristic of the field.

        OUTPUT:

        Instance of this class representing the image under
        the Frobenius morphisms.

        EXAMPLES::

            sage: F.<a> = GF(16)
            sage: b = a.conjugate(); b
            a + 1
            sage: a == b.conjugate()
            True

            sage: F.<a> = GF(27)
            sage: a.conjugate()
            Traceback (most recent call last):
            ...
            TypeError: cardinality of the field must be a square number

        TESTS:

        Check that :issue:`26761` is fixed::

            sage: # needs sage.libs.gap
            sage: G32 = GU(3,2)
            sage: g1, g2 = G32.gens()
            sage: m1 = g1.matrix()
            sage: m1.is_unitary()
            True
            sage: G32(m1) == g1
            True
        """
        k2 = self.parent().degree()
        if k2 % 2:
            raise TypeError("cardinality of the field must be a square number")
        k = k2 / 2

        return self.pth_power(k=k)

    def to_integer(self, reverse=False):
        r"""
        Return an integer representation of this finite field element
        obtained by lifting its representative polynomial to `\ZZ` and
        evaluating it at the characteristic `p`.

        If ``reverse`` is set to ``True`` (default: ``False``),
        the list of coefficients is reversed prior to evaluation.

        Inverse of :meth:`sage.rings.finite_rings.finite_field_base.FiniteField.from_integer`.

        EXAMPLES::

            sage: F.<t> = GF(7^5)
            sage: F(5).to_integer()
            5
            sage: t.to_integer()
            7
            sage: (t^2).to_integer()
            49
            sage: (t^2+1).to_integer()
            50
            sage: (t^2+t+1).to_integer()
            57

        ::

            sage: F.<t> = GF(2^8)
            sage: u = F.from_integer(0xd1)
            sage: bin(u.to_integer(False))
            '0b11010001'
            sage: bin(u.to_integer(True))
            '0b10001011'

        TESTS::

            sage: # needs sage.modules
            sage: p = random_prime(2^99)
            sage: k = randrange(2,10)
            sage: F.<t> = GF((p, k))
            sage: rev = bool(randrange(2))
            sage: u = F.random_element()
            sage: 0 <= u.to_integer(rev) < F.cardinality()
            True
            sage: F.from_integer(u.to_integer(rev), rev) == u
            True
            sage: n = randrange(F.cardinality())
            sage: F.from_integer(n, rev).to_integer(rev) == n
            True
        """
        if not reverse:
            try:
                return self._integer_representation()
            except AttributeError:
                pass
        p = self.parent().characteristic()
        f = self.polynomial().change_ring(ZZ)
        if reverse:
            f = f.reverse(self.parent().degree() - 1)
        return f(p)

    def to_bytes(self, byteorder='big'):
        r"""
        Return an array of bytes representing an integer.

        Internally relies on the python ``int.to_bytes()`` method.
        Length of byte array is determined from the field's order.

        INPUT:

        - ``byteorder`` -- string (default: ``'big'``); determines the byte order of
          the output; can only be ``'big'`` or ``'little'``

        EXAMPLES::

            sage: F.<z5> = GF(3^5)
            sage: a = z5^4 + 2*z5^3 + 1
            sage: a.to_bytes()
            b'\x88'

        ::

            sage: F.<z3> = GF(163^3)
            sage: a = 136*z3^2 + 10*z3 + 125
            sage: a.to_bytes()
            b'7)\xa3'
        """
        length = (self.parent().order().nbits() + 7) // 8
        return self.to_integer().to_bytes(length=length, byteorder=byteorder)

cdef class Cache_base(SageObject):
    cpdef FinitePolyExtElement fetch_int(self, number):
        r"""
        Given an integer less than `p^n` with base `2`
        representation `a_0 + a_1 \cdot 2 + \cdots + a_k 2^k`, this returns
        `a_0 + a_1 x + \cdots + a_k x^k`, where `x` is the
        generator of this finite field.

        EXAMPLES::

            sage: k.<a> = GF(2^48)
            sage: k._cache.fetch_int(2^33 + 2 + 1)                                      # needs sage.libs.ntl
            a^33 + a + 1
        """
        raise NotImplementedError("this must be implemented by subclasses")
