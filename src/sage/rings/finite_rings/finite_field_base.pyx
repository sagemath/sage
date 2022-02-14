# -*- coding: utf-8 -*-
"""
Base Classes for Finite Fields

TESTS::

    sage: K.<a> = NumberField(x^2 + 1)
    sage: F = K.factor(3)[0][0].residue_field()
    sage: loads(dumps(F)) == F
    True

AUTHORS:

- Adrien Brochard, David Roe, Jeroen Demeyer, Julian Rüth, Niles Johnson,
  Peter Bruin, Travis Scrimshaw, Xavier Caruso: initial version

"""
#*****************************************************************************
#       Copyright (C) 2009-2019 David Roe <roed.math@gmail.com>
#       Copyright (C) 2010      Niles Johnson <nilesj@gmail.com>
#       Copyright (C) 2011      Jeroen Demeyer <jdemeyer@cage.ugent.be>
#       Copyright (C) 2012      Adrien Brochard <aaa.brochard@gmail.com>
#       Copyright (C) 2012      Travis Scrimshaw <tscrim@ucdavis.edu>
#       Copyright (C) 2012      Xavier Caruso <xavier.caruso@normalesup.org>
#       Copyright (C) 2013      Peter Bruin <P.Bruin@warwick.ac.uk>
#       Copyright (C) 2014-2022 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
#*****************************************************************************

cimport cython
from cysignals.signals cimport sig_check
from cpython.array cimport array

from sage.categories.finite_fields import FiniteFields
from sage.structure.parent cimport Parent
from sage.misc.persist import register_unpickle_override
from sage.misc.cachefunc import cached_method
from sage.misc.prandom import randrange
from sage.rings.integer cimport Integer
import sage.rings.abc
from sage.misc.superseded import deprecation_cython as deprecation

# Copied from sage.misc.fast_methods, used in __hash__() below.
cdef int SIZEOF_VOID_P_SHIFT = 8*sizeof(void *) - 4

cdef class FiniteField(Field):
    """
    Abstract base class for finite fields.

    TESTS::

        sage: GF(997).is_finite()
        True
    """
    def __init__(self, base, names, normalize, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: K = GF(7); K
            Finite Field of size 7
            sage: loads(K.dumps()) == K
            True
            sage: GF(7^10, 'a')
            Finite Field in a of size 7^10
            sage: K = GF(7^10, 'a'); K
            Finite Field in a of size 7^10
            sage: loads(K.dumps()) == K
            True
        """
        if category is None:
            category = FiniteFields()
        Field.__init__(self, base, names, normalize, category)

    # The methods __hash__ and __richcmp__ below were copied from
    # sage.misc.fast_methods.WithEqualityById; we cannot inherit from
    # this since Cython does not support multiple inheritance.

    def __hash__(self):
        """
        The hash provided by this class coincides with that of ``<class 'object'>``.

        TESTS::

            sage: F.<a> = FiniteField(2^3)
            sage: hash(F) == hash(F)
            True
            sage: hash(F) == object.__hash__(F)
            True

        """
        # This is the default hash function in Python's object.c:
        cdef long x
        cdef size_t y = <size_t><void *>self
        y = (y >> 4) | (y << SIZEOF_VOID_P_SHIFT)
        x = <long>y
        if x==-1:
            x = -2
        return x

    def __richcmp__(self, other, int m):
        """
        Compare ``self`` with ``other``.

        Finite fields compare equal if and only if they are identical.
        In particular, they are not equal unless they have the same
        cardinality, modulus, variable name and implementation.

        EXAMPLES::

            sage: x = polygen(GF(3))
            sage: F = FiniteField(3^2, 'c', modulus=x^2+1)
            sage: F == F
            True
            sage: F == FiniteField(3^3, 'c')
            False
            sage: F == FiniteField(3^2, 'c', modulus=x^2+x+2)
            False
            sage: F == FiniteField(3^2, 'd')
            False
            sage: F == FiniteField(3^2, 'c', impl='pari_ffelt')
            False
        """
        if self is other:
            if m == 2: # ==
                return True
            elif m == 3: # !=
                return False
            else:
                # <= or >= or NotImplemented
                return m==1 or m==5 or NotImplemented
        else:
            if m == 2:
                return False
            elif m == 3:
                return True
            else:
                return NotImplemented

    def is_perfect(self):
        r"""
        Return whether this field is perfect, i.e., every element has a `p`-th
        root. Always returns ``True`` since finite fields are perfect.

        EXAMPLES::

            sage: GF(2).is_perfect()
            True

        """
        return True

    def _repr_(self):
        """
        String representation of this finite field.

        EXAMPLES::

            sage: k = GF(127)
            sage: k # indirect doctest
            Finite Field of size 127

            sage: k.<b> = GF(2^8)
            sage: k
            Finite Field in b of size 2^8

            sage: k.<c> = GF(2^20)
            sage: k
            Finite Field in c of size 2^20

            sage: k.<d> = GF(7^20)
            sage: k
            Finite Field in d of size 7^20
        """
        if self.absolute_degree()>1:
            return "Finite Field in %s of size %s^%s"%(self.variable_name(),self.characteristic(),self.absolute_degree())
        else:
            return "Finite Field of size %s"%(self.characteristic())

    def _latex_(self):
        r"""
        Returns a string denoting the name of the field in LaTeX.

        The :func:`~sage.misc.latex.latex` function calls the
        ``_latex_`` attribute when available.

        EXAMPLES:

        The ``latex`` command parses the string::

            sage: GF(81, 'a')._latex_()
            '\\Bold{F}_{3^{4}}'
            sage: latex(GF(81, 'a'))
            \Bold{F}_{3^{4}}
            sage: GF(3)._latex_()
            '\\Bold{F}_{3}'
            sage: latex(GF(3))
            \Bold{F}_{3}
        """
        if self.absolute_degree() > 1:
            e = "^{%s}"%self.absolute_degree()
        else:
            e = ""
        return "\\Bold{F}_{%s%s}"%(self.characteristic(), e)

    def _gap_init_(self):
        """
        Return string that initializes the GAP version of
        this finite field.

        EXAMPLES::

            sage: GF(9,'a')._gap_init_()
            'GF(9)'
        """
        return 'GF(%s)'%self.order()

    def _magma_init_(self, magma):
        """
        Return string representation of ``self`` that Magma can
        understand.

        EXAMPLES::

            sage: GF(97,'a')._magma_init_(magma)            # optional - magma
            'GF(97)'
            sage: GF(9,'a')._magma_init_(magma)             # optional - magma
            'SageCreateWithNames(ext<GF(3)|_sage_[...]![GF(3)!2,GF(3)!2,GF(3)!1]>,["a"])'
            sage: magma(GF(9,'a'))                          # optional - magma
            Finite field of size 3^2
            sage: magma(GF(9,'a')).1                        # optional - magma
            a
        """
        if self.absolute_degree() == 1:
            return 'GF(%s)'%self.order()
        B = self.base_ring()
        p = self.polynomial()
        s = "ext<%s|%s>"%(B._magma_init_(magma),p._magma_init_(magma))
        return magma._with_names(s, self.variable_names())

    def _macaulay2_init_(self, macaulay2=None):
        """
        Return the string representation of this finite field that Macaulay2
        can understand.

        Note that, in the case of a prime field, this returns ``ZZ/p`` instead
        of the Galois field ``GF p``, since computations in polynomial rings
        over ``ZZ/p`` are faster in Macaulay2 (as of 2019).

        EXAMPLES::

            sage: macaulay2(GF(97, 'a'))       # indirect doctest, optional - macaulay2
            ZZ
            --
            97
            sage: macaulay2(GF(49, 'a'))       # indirect doctest, optional - macaulay2
            GF 49

        TESTS:

        The variable name is preserved (:trac:`28566`)::

            sage: K = macaulay2(GF(49, 'b'))  # optional - macaulay2
            sage: K.gens()                    # optional - macaulay2
            {b}
            sage: K._sage_()                  # optional - macaulay2
            Finite Field in b of size 7^2
        """
        if self.is_prime_field():
            return "ZZ/%s" % self.order()
        return "GF(%s,Variable=>symbol %s)" % (self.order(),
                                               self.variable_name())

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __iter__(self):
        """
        Iterate over all elements of this finite field.

        EXAMPLES::

            sage: k.<a> = FiniteField(9, impl="pari")
            sage: list(k)
            [0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2]

        Partial iteration of a very large finite field::

            sage: p = next_prime(2^64)
            sage: k.<a> = FiniteField(p^2, impl="pari")
            sage: it = iter(k); it
            <generator object at ...>
            sage: [next(it) for i in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        TESTS:

        Check that the generic implementation works in all cases::

            sage: L = []
            sage: from sage.rings.finite_rings.finite_field_base import FiniteField
            sage: for impl in ("givaro", "pari", "ntl"):
            ....:     k = GF(8, impl=impl, names="z")
            ....:     print(list(FiniteField.__iter__(k)))
            [0, 1, z, z + 1, z^2, z^2 + 1, z^2 + z, z^2 + z + 1]
            [0, 1, z, z + 1, z^2, z^2 + 1, z^2 + z, z^2 + z + 1]
            [0, 1, z, z + 1, z^2, z^2 + 1, z^2 + z, z^2 + z + 1]
        """
        cdef Py_ssize_t n = self.absolute_degree()
        cdef unsigned long lim  # maximum value for coefficients
        try:
            lim = (<unsigned long>self.characteristic()) - 1
        except OverflowError:
            # If the characteristic is too large to represent in an
            # "unsigned long", it is reasonable to assume that we won't
            # be able to iterate over all elements. Typically n >= 2,
            # so that would mean >= 2^64 elements on a 32-bit system.
            # Just in case that we iterate to the end anyway, we raise
            # an exception at the end.
            lim = <unsigned long>(-1)

        elt = self.zero()
        one = self.one()
        x = self.absolute_gen()

        # coeffs[] is an array of coefficients of the current finite
        # field element (as unsigned longs)
        #
        # Stack represents an element
        #   sum_{i=0}^{n-1}  coeffs[i] x^i
        # as follows:
        #   stack[k] = sum_{i=k}^{n-1} coeffs[i] x^(i-k)
        #
        # This satisfies the recursion
        #   stack[k-1] = x * stack[k] + coeffs[k-1]
        #
        # Finally, elt is a shortcut for stack[k]
        #
        cdef list stack = [elt] * n
        cdef array coeffsarr = array('L', [0] * n)
        cdef unsigned long* coeffs = coeffsarr.data.as_ulongs

        yield elt  # zero
        cdef Py_ssize_t k = 0
        while k < n:
            # Find coefficients of next element
            coeff = coeffs[k]
            if coeff >= lim:
                # We cannot increase coeffs[k], so we wrap around to 0
                # and try to increase the next coefficient
                coeffs[k] = 0
                k += 1
                continue
            coeffs[k] = coeff + 1

            # Now compute and yield the finite field element
            sig_check()
            elt = stack[k] + one
            stack[k] = elt
            # Fix lower elements of stack, until k == 0
            while k > 0:
                elt *= x
                k -= 1
                stack[k] = elt

            yield elt

        if lim == <unsigned long>(-1):
            raise NotImplementedError("iterating over all elements of a large finite field is not supported")

    def fetch_int(self, n):
        r"""
        Return the element of ``self`` that equals `n` under the condition that
        :meth:`gen()` is set to the characteristic of the finite field ``self``.

        INPUT:

        - `n` -- integer. Must not be negative, and must be less than the
          cardinality of ``self``.

        EXAMPLES::

            sage: p = 4091
            sage: F = GF(p^4, 'a')
            sage: n = 100*p^3 + 37*p^2 + 12*p + 6
            sage: F.fetch_int(n)
            100*a^3 + 37*a^2 + 12*a + 6
            sage: F.fetch_int(n) in F
            True

        TESTS::

            sage: F = GF(19^5)
            sage: F.fetch_int(0)
            0
            sage: _.parent()
            Finite Field in ... of size 19^5
            sage: F.fetch_int(-5)
            Traceback (most recent call last):
            ...
            TypeError: n must be between 0 and self.order()
            sage: F.fetch_int(F.cardinality())
            Traceback (most recent call last):
            ...
            TypeError: n must be between 0 and self.order()
        """
        n = Integer(n)
        if (n < 0) or (n >= self.order()):
            raise TypeError("n must be between 0 and self.order()")
        if n == 0:
            return self.zero()
        cdef list digs = n.digits(base=self.characteristic())
        g = self.gen()
        return self.sum(self(d) * g**i for i, d in enumerate(digs) if d)

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        Return ``True`` if the map from self to codomain sending
        ``self.0`` to the unique element of ``im_gens`` is a valid field
        homomorphism. Otherwise, return ``False``.

        EXAMPLES::

        Between prime fields::

            sage: k0 = FiniteField(73, modulus='primitive')
            sage: k1 = FiniteField(73)
            sage: k0._is_valid_homomorphism_(k1, (k1(5),) )
            True
            sage: k1._is_valid_homomorphism_(k0, (k0(1),) )
            True

        Now for extension fields::

            sage: k.<a> = FiniteField(73^2)
            sage: K.<b> = FiniteField(73^3)
            sage: L.<c> = FiniteField(73^4)
            sage: k0._is_valid_homomorphism_(k, (k(5),) )
            True
            sage: k.hom([c]) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism

            sage: k.hom([c^(73*73+1)])
            Ring morphism:
            From: Finite Field in a of size 73^2
            To:   Finite Field in c of size 73^4
            Defn: a |--> 7*c^3 + 13*c^2 + 65*c + 71

            sage: k.hom([b])
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism
        """
        #if self.characteristic() != codomain.characteristic():
        #    raise ValueError("no map from %s to %s" % (self, codomain))
        # When the base is not just Fp, we want to ensure that there's a
        # coercion map from the base rather than just checking the characteristic
        if base_map is None and not codomain.has_coerce_map_from(self.base_ring()):
            return False
        if len(im_gens) != 1:
            raise ValueError("only one generator for finite fields")

        f = self.modulus()
        if base_map is not None:
            f = f.map_coefficients(base_map)
        return f(im_gens[0]).is_zero()

    def _Hom_(self, codomain, category=None):
        """
        Return the set of homomorphisms from ``self`` to ``codomain``
        in ``category``.

        This function is implicitly called by the ``Hom`` method or
        function.

        EXAMPLES::

            sage: K.<a> = GF(25); K
            Finite Field in a of size 5^2
            sage: K.Hom(K) # indirect doctest
            Automorphism group of Finite Field in a of size 5^2
        """
        from sage.rings.finite_rings.homset import FiniteFieldHomset
        if category.is_subcategory(FiniteFields()):
            return FiniteFieldHomset(self, codomain, category)
        return super(FiniteField, self)._Hom_(codomain, category)

    def _squarefree_decomposition_univariate_polynomial(self, f):
        """
        Return the square-free decomposition of this polynomial.  This is a
        partial factorization into square-free, coprime polynomials.

        This is a helper method for
        :meth:`sage.rings.polynomial.squarefree_decomposition`.

        INPUT:

        - ``f`` -- a univariate non-zero polynomial over this field

        ALGORITHM; [Coh1993]_, algorithm 3.4.2 which is basically the algorithm in
        [Yun1976]_ with special treatment for powers divisible by `p`.

        EXAMPLES::

            sage: K.<a> = GF(3^2)
            sage: R.<x> = K[]
            sage: f = x^243+2*x^81+x^9+1
            sage: f.squarefree_decomposition()
            (x^27 + 2*x^9 + x + 1)^9
            sage: f = x^243+a*x^27+1
            sage: f.squarefree_decomposition()
            (x^9 + (2*a + 1)*x + 1)^27

        TESTS::

            sage: for K in [GF(2^18,'a'), GF(3^2,'a'), GF(47^3,'a')]:
            ....:     R.<x> = K[]
            ....:     if K.characteristic() < 5: m = 4
            ....:     else: m = 1
            ....:     for _ in range(m):
            ....:         f = (R.random_element(4)^3*R.random_element(m)^(m+1))(x^6)
            ....:         F = f.squarefree_decomposition()
            ....:         assert F.prod() == f
            ....:         for i in range(len(F)):
            ....:             assert gcd(F[i][0], F[i][0].derivative()) == 1
            ....:             for j in range(len(F)):
            ....:                 if i == j: continue
            ....:                 assert gcd(F[i][0], F[j][0]) == 1
            ....:
        """
        from sage.structure.factorization import Factorization
        if f.degree() == 0:
            return Factorization([], unit=f[0])

        factors = []
        p = self.characteristic()
        unit = f.leading_coefficient()
        T0 = f.monic()
        e = 1
        if T0.degree() > 0:
            der = T0.derivative()
            while der.is_zero():
                T0 = T0.parent()([T0[p*i].pth_root() for i in range(T0.degree()//p + 1)])
                if T0 == 1:
                    raise RuntimeError
                der = T0.derivative()
                e = e*p
            T = T0.gcd(der)
            V = T0 // T
            k = 0
            while T0.degree() > 0:
                k += 1
                if p.divides(k):
                    T = T // V
                    k += 1
                W = V.gcd(T)
                if W.degree() < V.degree():
                    factors.append((V // W, e*k))
                    V = W
                    T = T // V
                    if V.degree() == 0:
                        if T.degree() == 0:
                            break
                        # T is of the form sum_{i=0}^n t_i X^{pi}
                        T0 = T0.parent()([T[p*i].pth_root() for i in range(T.degree()//p + 1)])
                        der = T0.derivative()
                        e = p*e
                        while der.is_zero():
                            T0 = T0.parent()([T0[p*i].pth_root() for i in range(T0.degree()//p + 1)])
                            der = T0.derivative()
                            e = p*e
                        T = T0.gcd(der)
                        V = T0 // T
                        k = 0
                else:
                    T = T//V

        return Factorization(factors, unit=unit, sort=False)

    def _any_embedding(self, codomain):
        r"""
        Return an embedding of this field into ``codomain``.

        EXAMPLES::

            sage: GF(2)._any_embedding(GF(2^2))
            Ring morphism:
              From: Finite Field of size 2
              To:   Finite Field in z2 of size 2^2
              Defn: 1 |--> 1
            sage: GF(2^2)._any_embedding(GF(2^4))
            Ring morphism:
              From: Finite Field in z2 of size 2^2
              To:   Finite Field in z4 of size 2^4
              Defn: z2 |--> z4^2 + z4
            sage: set_random_seed(0)
            sage: k = GF(3^2).extension(3, absolute=False); c = k.gen()
            sage: f = k._any_embedding(GF(3^12)); f
            Ring morphism:
              From: Finite Field in z6 of size 3^6 over its base
              To:   Finite Field in z12 of size 3^12
              Defn: z6 |--> z12^11 + 2*z12^9 + 2*z12^8 + 2*z12^7 + z12^5 + z12^3 + 2*z12^2 + 1
            sage: c.minpoly().change_ring(f)(f(c))
            0

        TESTS:

        Trivial extensions can be embedded into each other::

            sage: k = GF(2).extension(1, absolute=False)
            sage: k._any_embedding(k)
            sage: k._any_embedding(GF(2))
            sage: GF(2)._any_embedding(k)

        ::

            sage: k = GF(2).extension(1, absolute=False).extension(1, absolute=False)
            sage: k._any_embedding(k)
            sage: GF(2)._any_embedding(k)
            sage: k._any_embedding(GF(2))
        """
        if codomain.has_coerce_map_from(self):
            return codomain.coerce_map_from(self)

        base_hom = self.base_ring()._any_embedding(codomain)

        if self is self.prime_subfield():
            im_gens = []
        else:
            minpoly = self.gen().minpoly().change_ring(base_hom)
            im_gens = [minpoly.any_root()]

        return self.hom(im_gens, codomain=codomain, base_map=base_hom)

    def zeta_order(self):
        """
        Return the order of the distinguished root of unity in ``self``.

        EXAMPLES::

            sage: GF(9,'a').zeta_order()
            8
            sage: GF(9,'a').zeta()
            a
            sage: GF(9,'a').zeta().multiplicative_order()
            8
        """
        return self.order() - 1

    def zeta(self, n=None):
        """
        Return an element of multiplicative order ``n`` in this finite
        field. If there is no such element, raise ``ValueError``.

        .. WARNING::

            In general, this returns an arbitrary element of the correct
            order. There are no compatibility guarantees:
            ``F.zeta(9)^3`` may not be equal to ``F.zeta(3)``.

        EXAMPLES::

            sage: k = GF(7)
            sage: k.zeta()
            3
            sage: k.zeta().multiplicative_order()
            6
            sage: k.zeta(3)
            2
            sage: k.zeta(3).multiplicative_order()
            3
            sage: k = GF(49, 'a')
            sage: k.zeta().multiplicative_order()
            48
            sage: k.zeta(6)
            3
            sage: k.zeta(5)
            Traceback (most recent call last):
            ...
            ValueError: no 5th root of unity in Finite Field in a of size 7^2

        Even more examples::

            sage: GF(9,'a').zeta_order()
            8
            sage: GF(9,'a').zeta()
            a
            sage: GF(9,'a').zeta(4)
            a + 1
            sage: GF(9,'a').zeta()^2
            a + 1

        This works even in very large finite fields, provided that ``n``
        can be factored (see :trac:`25203`)::

            sage: k.<a> = GF(2^2000)
            sage: p = 8877945148742945001146041439025147034098690503591013177336356694416517527310181938001
            sage: z = k.zeta(p)
            sage: z
            a^1999 + a^1996 + a^1995 + a^1994 + ... + a^7 + a^5 + a^4 + 1
            sage: z ^ p
            1
        """
        if n is None:
            return self.multiplicative_generator()

        n = Integer(n)
        grouporder = self.order() - 1
        co_order = grouporder // n
        if co_order * n != grouporder:
            raise ValueError("no {}th root of unity in {}".format(n, self))

        # If the co_order is small or we know a multiplicative
        # generator, use a multiplicative generator
        mg = self.multiplicative_generator
        if mg.cache is not None or co_order <= 500000:
            return mg() ** co_order
        return self._element_of_factored_order(n.factor())

    @cached_method
    def multiplicative_generator(self):
        """
        Return a primitive element of this finite field, i.e. a
        generator of the multiplicative group.

        You can use :meth:`multiplicative_generator()` or
        :meth:`primitive_element()`, these mean the same thing.

        .. WARNING::

           This generator might change from one version of Sage to another.

        EXAMPLES::

            sage: k = GF(997)
            sage: k.multiplicative_generator()
            7
            sage: k.<a> = GF(11^3)
            sage: k.primitive_element()
            a
            sage: k.<b> = GF(19^32)
            sage: k.multiplicative_generator()
            b + 4

        TESTS:

        Check that large characteristics work (:trac:`11946`)::

            sage: p = 10^20 + 39
            sage: x = polygen(GF(p))
            sage: K.<a> = GF(p^2, modulus=x^2+1)
            sage: K.multiplicative_generator()
            a + 12
        """
        if self.absolute_degree() == 1:
            from sage.arith.all import primitive_root
            return self(primitive_root(self.order()))
        F, = self.factored_unit_order()
        return self._element_of_factored_order(F)

    primitive_element = multiplicative_generator

    def _element_of_factored_order(self, F):
        """
        Return an element of ``self`` of order ``n`` where ``n`` is
        given in factored form.

        INPUT:

        - ``F`` -- the factorization of the required order. The order
          must be a divisor of ``self.order() - 1`` but this is not
          checked.

        EXAMPLES::

            sage: k.<a> = GF(16, modulus=cyclotomic_polynomial(5))
            sage: k._element_of_factored_order(factor(15))
            a^2 + a + 1
            sage: k._element_of_factored_order(factor(5))
            a^3
            sage: k._element_of_factored_order(factor(3))
            a^3 + a^2 + 1
            sage: k._element_of_factored_order(factor(30))
            Traceback (most recent call last):
            ...
            AssertionError: no element found
        """
        n = Integer(1)
        cdef list primes = []
        for p, e in F:
            primes.append(p)
            n *= p ** e

        N = self.order() - 1
        c = N // n

        # We check whether (x + g)^c has the required order, where
        # x runs through the finite field.
        # This has the advantage that g is the first element we try,
        # so if that was a chosen to be a multiplicative generator,
        # we are done immediately. Second, the PARI finite field
        # iterator gives all the constant elements first, so we try
        # (g+(constant))^c before anything else.
        g = self.gen(0)
        for x in self:
            a = (g + x) ** c
            if not a:
                continue
            if all(a ** (n // p) != 1 for p in primes):
                return a
        raise AssertionError("no element found")

    def ngens(self):
        """
        The number of generators of the finite field.  Always 1.

        EXAMPLES::

            sage: k = FiniteField(3^4, 'b')
            sage: k.ngens()
            1
        """
        return 1

    def is_field(self, proof = True):
        """
        Return whether or not the finite field is a field, i.e.,
        always return ``True``.

        EXAMPLES::

            sage: k.<a> = FiniteField(3^4)
            sage: k.is_field()
            True
        """
        return True

    def order(self):
        """
        Return the order of this finite field.

        EXAMPLES::

            sage: GF(997).order()
            997
        """
        return self.characteristic()**self.absolute_degree()

    # cached because constructing the Factorization is slow;
    # see trac #11628.
    @cached_method
    def factored_order(self):
        """
        Returns the factored order of this field.  For compatibility with
        :mod:`~sage.rings.finite_rings.integer_mod_ring`.

        EXAMPLES::

            sage: GF(7^2,'a').factored_order()
            7^2
        """
        from sage.structure.factorization import Factorization
        return Factorization([(self.characteristic(), self.absolute_degree())])

    @cached_method
    def factored_unit_order(self):
        """
        Returns the factorization of ``self.order()-1``, as a 1-tuple.

        The format is for compatibility with
        :mod:`~sage.rings.finite_rings.integer_mod_ring`.

        EXAMPLES::

            sage: GF(7^2,'a').factored_unit_order()
            (2^4 * 3,)

        TESTS:

        Check that :trac:`31686` is fixed::

            sage: p = 1100585370631
            sage: F = GF(p^24, 'a')
            sage: F.factored_unit_order()
            (2^6 * 3^2 * 5 * 7 * 11 * 13 * 17 * 53 * 97 * 229 * 337 * 421
             * 3929 * 215417 * 249737 * 262519 * 397897 * 59825761 * 692192057
             * 12506651939 * 37553789761 * 46950147799 * 172462808473 * 434045140817
             * 81866093016401 * 617237859576697 * 659156729361017707
             * 268083135725348991493995910983015600019336657
             * 90433843562394341719266736354746485652016132372842876085423636587989263202299569913,)
        """
        from sage.structure.factorization import Factorization
        from sage.rings.polynomial.cyclotomic import cyclotomic_value as cv
        p, d = self.characteristic(), self.degree()
        F = Factorization(f for n in d.divisors() for f in cv(n, p).factor())
        return (F,)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        Same as :meth:`order`.

        EXAMPLES::

            sage: GF(997).cardinality()
            997
        """
        return self.order()

    __len__ = cardinality

    def is_prime_field(self):
        """
        Return ``True`` if ``self`` is a prime field, i.e., has degree 1.

        EXAMPLES::

            sage: GF(3^7, 'a').is_prime_field()
            False
            sage: GF(3, 'a').is_prime_field()
            True
        """
        return self.absolute_degree() == 1

    @cached_method
    def modulus(self):
        r"""
        Return the minimal polynomial of the generator of this field over its
        base field.

        The minimal polynomial of an element `a` in a field is the
        unique monic irreducible polynomial of smallest degree with
        coefficients in the base field that has `a` as a root. In
        finite absolute field extensions, `\GF{p^n}`, the base field is
        `\GF{p}`.

        OUTPUT:

        - a monic polynomial over the base ring in the variable `x`.

        EXAMPLES::

            sage: F.<a> = GF(7^2); F
            Finite Field in a of size 7^2
            sage: F.polynomial_ring()
            Univariate Polynomial Ring in a over Finite Field of size 7
            sage: f = F.modulus(); f
            x^2 + 6*x + 3
            sage: f(a)
            0

        Although `f` is irreducible over the base field, we can double-check
        whether or not `f` factors in `F` as follows. The command
        ``F['x'](f)`` coerces `f` as a polynomial with coefficients in `F`.
        (Instead of a polynomial with coefficients over the base field.)

        ::

            sage: f.factor()
            x^2 + 6*x + 3
            sage: F['x'](f).factor()
            (x + a + 6) * (x + 6*a)

        Here is an example with a degree 3 extension::

            sage: G.<b> = GF(7^3); G
            Finite Field in b of size 7^3
            sage: g = G.modulus(); g
            x^3 + 6*x^2 + 4
            sage: g.degree(); G.degree()
            3
            3

        For prime fields, this returns `x - 1` unless a custom modulus
        was given when constructing this field::

            sage: k = GF(199)
            sage: k.modulus()
            x + 198
            sage: var('x')
            x
            sage: k = GF(199, modulus=x+1)
            sage: k.modulus()
            x + 1

        The given modulus is always made monic::

            sage: k.<a> = GF(7^2, modulus=2*x^2-3, impl="pari_ffelt")
            sage: k.modulus()
            x^2 + 2

        TESTS:

        We test the various finite field implementations::

            sage: GF(2, impl="modn").modulus()
            x + 1
            sage: GF(2, impl="givaro").modulus()
            x + 1
            sage: GF(2, impl="ntl").modulus()
            x + 1
            sage: GF(2, impl="modn", modulus=x).modulus()
            x
            sage: GF(2, impl="givaro", modulus=x).modulus()
            x
            sage: GF(2, impl="ntl", modulus=x).modulus()
            x
            sage: GF(13^2, 'a', impl="givaro", modulus=x^2+2).modulus()
            x^2 + 2
            sage: GF(13^2, 'a', impl="pari_ffelt", modulus=x^2+2).modulus()
            x^2 + 2
            sage: set_random_seed(0); GF(4).extension(2, absolute=False).modulus()
            x^2 + z2*x + z2

        """
        # Normally, this is set by the constructor of the implementation
        try:
            return self._modulus
        except AttributeError:
            pass

        if self.base() is not self:
            raise NotImplementedError("finite field does not implement modulus yet")

        from sage.rings.all import PolynomialRing
        from .finite_field_constructor import GF
        R = PolynomialRing(GF(self.characteristic()), 'x')
        self._modulus = R((-1,1))  # Polynomial x - 1
        return self._modulus

    def polynomial(self, name=None):
        """
        Return the minimal polynomial of the generator of this field over
        its base field.

        INPUT:

        - ``name`` -- a variable name to use for the polynomial. By
          default, use the name given when constructing this field.

        OUTPUT:

        - a monic polynomial over the base field in the variable ``name``.

        .. SEEALSO::

            Except for the ``name`` argument, this is identical to the
            :meth:`modulus` method.

        EXAMPLES::

            sage: k.<a> = FiniteField(9)
            sage: k.polynomial('x')
            x^2 + 2*x + 2
            sage: k.polynomial()
            a^2 + 2*a + 2

            sage: F = FiniteField(9, 'a', impl='pari_ffelt')
            sage: F.polynomial()
            a^2 + 2*a + 2

            sage: F = FiniteField(7^20, 'a', impl='pari_ffelt')
            sage: f = F.polynomial(); f
            a^20 + a^12 + 6*a^11 + 2*a^10 + 5*a^9 + 2*a^8 + 3*a^7 + a^6 + 3*a^5 + 3*a^3 + a + 3
            sage: f(F.gen())
            0

            sage: k.<a> = GF(2^20, impl='ntl')
            sage: k.polynomial()
            a^20 + a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1
            sage: k.polynomial('FOO')
            FOO^20 + FOO^10 + FOO^9 + FOO^7 + FOO^6 + FOO^5 + FOO^4 + FOO + 1
            sage: a^20
            a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1
        """
        if name is None:
            name = self.variable_name()
        return self.modulus().change_variable_name(name)

    def unit_group_exponent(self):
        """
        The exponent of the unit group of the finite field.  For a
        finite field, this is always the order minus 1.

        EXAMPLES::

            sage: k = GF(2^10, 'a')
            sage: k.order()
            1024
            sage: k.unit_group_exponent()
            1023
        """
        return self.order() - 1


    def random_element(self, *args, **kwds):
        r"""
        A random element of the finite field.  Passes arguments to
        ``random_element()`` function of underlying vector space.

        EXAMPLES::

            sage: k = GF(19^4, 'a')
            sage: k.random_element().parent() is k
            True

        Passes extra positional or keyword arguments through::

            sage: k.random_element(prob=0)
            0

        """
        if self.absolute_degree() == 1:
            return self(randrange(self.order()))
        v = self.vector_space(map=False).random_element(*args, **kwds)
        return self(v)

    def some_elements(self):
        """
        Return a collection of elements of this finite field for use in unit
        testing.

        EXAMPLES::

            sage: k = GF(2^8,'a')
            sage: k.some_elements() # random output
            [a^4 + a^3 + 1, a^6 + a^4 + a^3, a^5 + a^4 + a, a^2 + a]
        """
        return [self.random_element() for i in range(4)]

    @cached_method
    def polynomial_ring(self, variable_name=None):
        """
        Return the polynomial ring over the base field in the same variable as
        this finite field.

        EXAMPLES::

            sage: k.<alpha> = FiniteField(3^4)
            sage: k.polynomial_ring()
            Univariate Polynomial Ring in alpha over Finite Field of size 3
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        if variable_name is None:
            variable_name = self.variable_name()
        return PolynomialRing(self.base_ring(), variable_name)

    cpdef _coerce_map_from_(self, R):
        r"""
        Canonical coercion to this finite field.

        TESTS::

            sage: k.<a> = GF(2^8)
            sage: a + 1
            a + 1
            sage: a + int(1)
            a + 1
            sage: a + GF(2)(1)
            a + 1

            sage: k.<a> = GF(3^8)
            sage: a + 1
            a + 1
            sage: a + int(1)
            a + 1
            sage: a + GF(3)(1)
            a + 1

            sage: k = GF(4, 'a')
            sage: k._coerce_(GF(2)(1))
            1
            sage: k._coerce_(k.0)
            a
            sage: k._coerce_(3)
            1
            sage: k._coerce_(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Finite Field in a of size 2^2

            sage: FiniteField(16)._coerce_(FiniteField(4).0)
            z4^2 + z4

            sage: FiniteField(8, 'a')._coerce_(FiniteField(4, 'a').0)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Finite Field in a of size 2^2 to Finite Field in a of size 2^3

            sage: FiniteField(8, 'a')._coerce_(FiniteField(7, 'a')(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Finite Field of size 7 to Finite Field in a of size 2^3

        There is no coercion from a `p`-adic ring to its residue field::

            sage: R.<a> = Zq(81); k = R.residue_field()
            sage: k.has_coerce_map_from(R)
            False
        """
        from sage.rings.integer_ring import ZZ
        if R is int or R is long or R is ZZ:
            return True
        if isinstance(R, sage.rings.abc.IntegerModRing) and self.characteristic().divides(R.characteristic()):
            return R.hom((self.one(),), check=False)
        if is_FiniteField(R):
            if R is self:
                return True
            from .residue_field import ResidueField_generic
            if isinstance(R, ResidueField_generic):
                return False
            if R.characteristic() == self.characteristic():
                if R.absolute_degree() == 1:
                    return R.hom((self.one(),), check=False)
                elif (R.absolute_degree().divides(self.absolute_degree())
                      and hasattr(self, '_prefix') and hasattr(R, '_prefix')):
                    return R.hom((self.gen() ** ((self.order() - 1)//(R.order() - 1)),))

    cpdef _convert_map_from_(self, R):
        """
        Conversion from p-adic fields.

        EXAMPLES::

            sage: K.<a> = Qq(49); k = K.residue_field()
            sage: k.convert_map_from(K)
            Reduction morphism:
              From: 7-adic Unramified Extension Field in a defined by x^2 + 6*x + 3
              To:   Finite Field in a0 of size 7^2

        Check that :trac:`8240 is resolved::

            sage: R.<a> = Zq(81); k = R.residue_field()
            sage: k.convert_map_from(R)
            Reduction morphism:
              From: 3-adic Unramified Extension Ring in a defined by x^4 + 2*x^3 + 2
              To:   Finite Field in a0 of size 3^4
        """
        from sage.rings.padics.padic_generic import pAdicGeneric, ResidueReductionMap
        if isinstance(R, pAdicGeneric) and R.residue_field() is self:
            return ResidueReductionMap._create_(R, self)

    def construction(self):
        """
        Return the construction of this finite field, as a ``ConstructionFunctor``
        and the base field.

        EXAMPLES::

            sage: v = GF(3^3).construction(); v
            (AlgebraicExtensionFunctor, Finite Field of size 3)
            sage: v[0].polys[0]
            3
            sage: v = GF(2^1000, 'a').construction(); v[0].polys[0]
            a^1000 + a^5 + a^4 + a^3 + 1

        The implementation is taken into account, by :trac:`15223`::

            sage: k = FiniteField(9, 'a', impl='pari_ffelt')
            sage: F, R = k.construction()
            sage: F(R) is k
            True

        """
        from sage.categories.pushout import AlgebraicExtensionFunctor
        try:
            kwds = {'impl': self._factory_data[2][4]}
        except (AttributeError, IndexError, TypeError):
            kwds = {}
        kwds["implementations"] = ["GF"]
        kwds["absolute"] = False
        if self.base() is self:
            # this is not of type FiniteField_prime_modn
            from sage.rings.integer import Integer
            return AlgebraicExtensionFunctor([self.polynomial()], [None], [None], **kwds), self.base_ring()
        elif hasattr(self, '_prefix'):
            return (AlgebraicExtensionFunctor([self.degree()], [self.variable_name()], [None],
                                              prefix=self._prefix, **kwds),
                    self.base_ring())
        else:
            return (AlgebraicExtensionFunctor([self.polynomial()], [self.variable_name()], [None], **kwds),
                    self.base_ring())

    def extension(self, modulus, name=None, names=None, map=False,
                  *, latex_name=None, latex_names=None, absolute=None, implementation=None, **kwds):
        """
        Return an extension of this finite field.

        INPUT:

        - ``modulus`` -- a polynomial with coefficients in ``self``,
          or an integer.

        - ``name`` or ``names`` -- string: the name of the generator
          in the new extension

        - ``latex_name`` or ``latex_names`` -- string: latex name of
          the generator in the new extension

        - ``map`` -- boolean (default: ``False``): if ``False``,
          return just the extension `E`; if ``True``, return a pair
          `(E, f)`, where `f` is an embedding of ``self`` into `E`.

        - ``absolute`` -- boolean (default: ``None``): if ``True``,
          return the absolute finite field of this extension.
          This is currently the default behavior when ``modulus`` is
          an integer, but we will be changing the default to always be relative
          after a deprecation period.

        - ``implementation`` -- string, either ``"GF"`` or ``"PQR"``.
          In the first case we use relative extensions of finite fields,
          while in the second we use the generic ``PolynomialQuotientRing``.
          The default is currently ``"PQR"``, but this is only for backward
          compatibility.  It will be changed to ``"GF"`` eventually in the
          case that the modulus is irreducible.

        - ``**kwds``: further keywords, passed to the finite field
          constructor.

        OUTPUT:

        An extension of the given modulus, or pseudo-Conway of the
        given degree if ``modulus`` is an integer.

        EXAMPLES::

            sage: k = GF(2)
            sage: R.<x> = k[]
            sage: k.extension(x^1000 + x^5 + x^4 + x^3 + 1, 'a', implementation="GF")
            Finite Field in a of size 2^1000
            sage: k = GF(3^4)
            sage: R.<x> = k[]
            sage: k.extension(3, absolute=True)
            Finite Field in z12 of size 3^12

        The base ring will have a coercion to the extension if it is constructed
        as a relative extension or if you don't specify variable names (so that it
        is part of the same algebraic closure)::

            sage: K = k.extension(2, 'a', absolute=False)
            sage: k.is_subring(K)
            True
            sage: K = k.extension(2, absolute=True)
            sage: k.is_subring(K)
            True
            sage: K = k.extension(2, 'a', absolute=True)
            sage: k.is_subring(K)
            False

        An example using the ``map`` argument::

            sage: F = GF(5)
            sage: E, f = F.extension(2, 'b', map=True, absolute=True)
            sage: E
            Finite Field in b of size 5^2
            sage: f
            Ring morphism:
              From: Finite Field of size 5
              To:   Finite Field in b of size 5^2
              Defn: 1 |--> 1
            sage: f.parent()
            Set of field embeddings from Finite Field of size 5 to Finite Field in b of size 5^2

        Extensions of non-prime finite fields by polynomials are now supported::

            sage: l.<b> = k.extension(x^5 + x^2 + x - 1, implementation="GF"); l
            Finite Field in b of size 3^20 over its base

        A relative extensions::

            sage: l.extension(2, absolute=False)
            Finite Field in z40 of size 3^40 over its base

        TESTS:

        We check that :trac:`18915` is fixed::

            sage: F = GF(2)
            sage: F.extension(int(3), 'a', absolute=True)
            Finite Field in a of size 2^3

            sage: F = GF(2 ** 4, 'a')
            sage: F.extension(int(3), 'aa', absolute=True)
            Finite Field in aa of size 2^12
        """
        from .finite_field_constructor import GF
        from sage.rings.all import PolynomialRing
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        from sage.rings.integer import Integer
        from sage.misc.superseded import deprecation

        if name is None and names is not None:
            name = names
        if names is None and name is not None:
            names = (name,)
        if latex_name is None and latex_names is not None:
            latex_name = latex_names

        if implementation not in [None, "GF", "PQR"]:
            raise ValueError("Unrecognized implementation %s" % implementation)

        if isinstance(modulus, (list, tuple)):
            modulus = PolynomialRing(self, 'x')(modulus)

        if isinstance(modulus, (int, Integer)):
            relative_degree = Integer(modulus)
            modulus = None
            is_field = True
            if implementation == "PQR":
                # We could find a defining polynomial and use it to create a PQR....
                raise NotImplementedError("PQR implementation incompatible with just specifying extension degree")
            if absolute is None:
                if self.base() is not self:
                    deprecation(28485, "extension() will produce an absolute extension since you did not specify the `absolute` keyword argument. In the future this will default to a relative extension. If you want to keep the current behaviour, use or extension(…, absolute=True), otherwise use extension(…, absolute=False)")
                # Once the above deprecation warning has been removed, add the following one.
                # (I know this is silly, but we can't think of a better way to do this)
                # from sage.misc.superseded import deprecation
                # if absolute:
                #     deprecation(28485, "the absolute keyword is deprecated. Use .extension(…).absolute_field() instead.")
                # else:
                #     deprecation(28485, "the absolute keyword is deprecated. There is no need to specify it explicitly to create a relative extension.")
                absolute = True
            implementation = "GF"
        elif is_Polynomial(modulus):
            relative_degree = modulus.degree()
            modulus = modulus.change_ring(self)
            is_field = modulus.is_irreducible()
            if implementation is None:
                implementation = "PQR"
                if is_field:
                    deprecation(28485, "extension() will produce a polynomial quotient ring since you did not specify the `implementation` keyword. In the future this will change to producing a relative extension of finite fields. If you want to keep the current behavior, use extension(…, implementation='PQR'), otherwise use extension(…, implementation='GF')")
            if implementation == "GF":
                if not is_field:
                    raise ValueError("GF implementation only available for field extensions")
                if absolute is None:
                    absolute = False
            else:
                if absolute:
                    raise NotImplementedError("absolute extensions not available for polynomial quotient rings")
        else:
            raise TypeError("unsupported modulus for finite field extension")

        absolute_degree = self.absolute_degree() * relative_degree

        if not absolute and absolute_degree == 1:
            # To get a trivial extension, we need to provide an explicit modulus
            modulus = PolynomialRing(self, 'x').gen()

        # Create the actual (relative) extension ring E
        if implementation == "GF":
            E = GF(self.characteristic() ** absolute_degree, name=name, modulus=modulus, base=self, **kwds)
        else:
            E = Field.extension(self, modulus, name=name, latex_name=latex_name, **kwds)
        if map:
            self_to_E = E.coerce_map_from(self)
            if self_to_E is None:
                assert absolute, "created a relative extension but there is no coercion embedding"
                # This code block can be dropped once absolute has been removed.
                self_absolute, absolute_to_self, self_to_absolute = self.absolute_field(map=True)
                absolute_to_E = self_absolute.hom([self_absolute.modulus().change_ring(E).any_root()], E)
                base_hom = absolute_to_E * self_to_absolute * self.coerce_map_from(self.base_ring())
                self_to_E = self.hom([self.modulus().change_ring(base_hom).any_root()], base_map=base_hom)

        if absolute:
            F = E.absolute_field(map=map, names=names)
            if map:
                F, F_to_E, E_to_F = F
                
                return F, E_to_F * self_to_E
            else:
                return F
        else:
            if map:
                return E, self_to_E
            else:
                return E

    def absolute_field(self, map=False, names=None):
        r"""
        Return an isomorphic finite field that is an absolute extension of the
        underlying prime field.

        EXAMPLES::

            TODO: Needs testing
            
        """
        if self.base() is not self.base_ring():
            raise NotImplementedError("this relative finite field extension does not implement absolute_field() yet")

        E = self.change(names=names) 
        if map:
            return (E, E.hom([self.gen()], check=False), self.hom([E.gen()], check=False))
        else:
            return E

    def change(self, **kwds):
        r"""
        Return a finite field with some characteristics changed.

        EXAMPLES::

            sage: GF(4).change(name='a')
            Finite Field in a of size 2^2

        All the parameters that
        :meth:`sage.rings.finite_rings.finite_field_constructor.FiniteFieldConstructor`
        takes are in principle supported::

            sage: GF(2).change(order=3)
            Finite Field of size 3

        However, some things might not work out immediately as expected::

            sage: GF(64).change(order=3)
            Traceback (most recent call last):
            ...
            ValueError: degree of the modulus does not match the relative degree of finite field

            sage: GF(64).change(order=3, modulus=None)
            Finite Field of size 3

        TESTS:

        Check that this syntactic sugar works as expected::

            sage: l.<b> = GF(4).change(); l
            Finite Field in b of size 2^2

        Arguments are checked::

            sage: GF(2).change(p=3)
            Traceback (most recent call last):
            ...
            ValueError: can only change arguments that GF understands when changing a finite field

        """
        from sage.all import GF
        return GF._change(self, **kwds)
        # TODO: Add test: switching trivial stuff back and forth works, say the name.

    @cached_method
    def _compatible_family(self):
        """
        Return a family of elements of this field that generate each subfield in a compatible way.

        OUTPUT:

        - A dictionary `D` so that if `n` is a positive integer dividing the degree of this field then
        ``D[n] = (a, f)`` where `a` generates the subfield of order `p^n` and `f` is the minimal polynomial of `a`.
        Moreover, if `a` and `b` are elements in this family of degree `m` and `n` respectively and `m` divides `n`
        then `a = b^{(p^n-1)/(p^m-1)}`.

        EXAMPLES::

            sage: k.<a> = GF(3^72)
            sage: F = k._compatible_family()
            sage: all(f(b) == 0 for (b, f) in F.values())
            True
            sage: all(f.degree() == n for (n, (b, f)) in F.items())
            True
            sage: D = 72.divisors()
            sage: for (m,n) in zip(D, D):
            ....:     if (n/m) in [2,3]:
            ....:         b, c = F[m][0], F[n][0]
            ....:         assert c^((3^n-1)//(3^m-1)) == b

        TESTS::

            sage: k = GF(9).extension(3, absolute=False)
            sage: F = k._compatible_family()
            sage: sorted(F)
            [1, 2, 3, 6]
        """
        p = self.characteristic()
        # We try to use the appropriate power of the generator,
        # as in the definition of Conway polynomials.
        # this can fail if the generator is not a primitive element,
        # but it can succeed sometimes even
        # if the generator is not primitive.
        g = self.gen()
        f = self.modulus()
        d = self.absolute_degree()
        D = list(reversed(d.divisors()[:-1]))
        P = d.support()
        def make_family(gen, poly):
            if poly.degree() != d:
                return False, {}
            fam = {d: (gen, poly)}
            for n in D:
                for l in P:
                    if l*n in fam:
                        a, _ = fam[l*n]
                        b = a**((p**(l*n) - 1)//(p**n - 1))
                        bpoly = b.minimal_polynomial()
                        if bpoly.degree() != n:
                            return False, fam
                        fam[n] = (b, bpoly)
                        break
            return True, fam
        while True:
            ok, fam = make_family(g, f)
            if ok:
                return fam
            g = self.random_element()
            f = g.minimal_polynomial()

    def subfield(self, degree, name=None, map=False):
        """
        Return the subfield of the field of given absolute degree.

        The inclusion maps between these subfields will always commute, but they are only added as coercion maps
        if the following condition holds for the generator `g` of the field, where `d` is the degree of this field
        over the prime field:

        The element `g^{(p^d - 1)/(p^n - 1)}` generates the subfield of degree `n` for all divisors `n` of `d`.

        The resulting finite field will be an absolute extension of the prime field, even if the input is a relative finite field.

        INPUT:

        - ``degree`` -- integer; degree of the subfield

        - ``name`` -- string (optional); name of the generator of the subfield.

        - ``map`` -- boolean (default ``False``); whether to also return the inclusion map

        EXAMPLES::

            sage: k = GF(2^21)
            sage: k.subfield(3)
            Finite Field in z3 of size 2^3
            sage: k0 = k.subfield(7, 'a'); k0
            Finite Field in a of size 2^7
            sage: k.coerce_map_from(k0)
            Ring morphism:
              From: Finite Field in a of size 2^7
              To:   Finite Field in z21 of size 2^21
              Defn: a |--> z21^20 + z21^19 + z21^17 + z21^15 + z21^14 + z21^6 + z21^4 + z21^3 + z21
            sage: k.subfield(8)
            Traceback (most recent call last):
            ...
            ValueError: no subfield of order 2^8

        TESTS:

        We check that :trac:`23801` is resolved::

            sage: k.<a> = GF(5^240)
            sage: l, inc = k.subfield(3, 'z', map=True); l
            Finite Field in z of size 5^3
            sage: inc
            Ring morphism:
              From: Finite Field in z of size 5^3
              To:   Finite Field in a of size 5^240
              Defn: z |--> ...

        There may be no coercion in the case that we can't ensure
        compatibility with larger fields::

            sage: R.<x> = GF(2)[]
            sage: k.<a> = GF(2^12, modulus=x^12 + x^3 + 1)
            sage: a.multiplicative_order()
            45
            sage: l = k.subfield(3)
            sage: k.has_coerce_map_from(l)
            False

        But there is still a compatibility among the generators chosen for the subfields::

            sage: ll, iinc = k.subfield(6, 'w', map=True)
            sage: t = iinc(ll.gen())^((2^6-1)/(2^3-1))
            sage: t.minimal_polynomial() == l.modulus()
            True

            sage: S = GF(37^16).subfields()
            sage: len(S) == len(16.divisors())
            True
            sage: all(f is not None for (l, f) in S)
            True

            sage: S = GF(2^93).subfields()
            sage: len(S) == len(93.divisors())
            True
            sage: all(f is not None for (l, f) in S)
            True

        We choose a default variable name::

            sage: GF(3^8, 'a').subfield(4)
            Finite Field in a4 of size 3^4

        Check that it works for relative finite fields::

            sage: set_random_seed(0)
            sage: k = GF(9).extension(3, absolute=False)
            sage: k.subfield(3, map=True)
            (Finite Field in z63 of size 3^3,
             Ring morphism:
               From: Finite Field in z3 of size 3^3
               To:   Finite Field in z6 of size 3^6 over its base
               Defn: z3 |--> (2*z2 + 1) - z6 + z2*z6^2)
        """
        from .finite_field_constructor import GF
        p = self.characteristic()
        n = self.absolute_degree()
        if not n % degree == 0:
            raise ValueError("no subfield of order {}^{}".format(p, degree))

        if degree == 1:
            K = self.prime_subfield()
            inc = self.coerce_map_from(K)
        elif degree == n:
            K = self
            inc = self.coerce_map_from(self)
        elif hasattr(self, '_prefix') and name is None:
            K = GF((p, degree), prefix=self._prefix)
            a = self.gen()**((p**n-1)//(p**degree - 1))
            inc = K.hom([a], codomain=self, check=False)
        else:
            if name is None:
                name = self.variable_name().rstrip('0123456789') + str(degree)
            fam = self._compatible_family()
            a, modulus = fam[degree]
            K = GF((p, degree), modulus=modulus, name=name)
            inc = K.hom([a], codomain=self, check=False)
            if fam[n][0] == self.gen():
                try: # to register a coercion map, embedding of K to self
                    self.register_coercion(inc)
                except AssertionError: # coercion already exists
                    pass
        if map:
            return K, inc
        else:
            return K

    def subfields(self, degree=0, name=None):
        """
        Return all subfields of the given ``degree``, or of all possible
        degrees if ``degree`` is `0`.

        The subfields are returned as absolute fields together with an
        embedding into this field.

        INPUT:

        - ``degree`` -- (default: `0`) an integer

        - ``name`` -- a string, a dictionary or ``None``:

          - If ``degree`` is nonzero, then ``name`` must be a string
            (or ``None``, if this is a pseudo-Conway extension),
            and will be the variable name of the returned field.
          - If ``degree`` is zero, the dictionary should have keys the divisors
            of the degree of this field, with the desired variable name for the
            field of that degree as an entry.
          - As a shortcut, you can provide a string and the degree of each
            subfield will be appended for the variable name of that subfield.
          - If ``None``, uses the prefix of this field.

        OUTPUT:

        A list of pairs ``(K, e)``, where ``K`` ranges over the subfields of
        this field and ``e`` gives an embedding of ``K`` into this field.

        EXAMPLES::

            sage: k = GF(2^21)
            sage: k.subfields()
            [(Finite Field of size 2,
              Ring morphism:
                  From: Finite Field of size 2
                  To:   Finite Field in z21 of size 2^21
                  Defn: 1 |--> 1),
             (Finite Field in z3 of size 2^3,
              Ring morphism:
                  From: Finite Field in z3 of size 2^3
                  To:   Finite Field in z21 of size 2^21
                  Defn: z3 |--> z21^20 + z21^19 + z21^17 + z21^15 + z21^11 + z21^9 + z21^8 + z21^6 + z21^2),
             (Finite Field in z7 of size 2^7,
              Ring morphism:
                  From: Finite Field in z7 of size 2^7
                  To:   Finite Field in z21 of size 2^21
                  Defn: z7 |--> z21^20 + z21^19 + z21^17 + z21^15 + z21^14 + z21^6 + z21^4 + z21^3 + z21),
             (Finite Field in z21 of size 2^21,
              Identity endomorphism of Finite Field in z21 of size 2^21)]
        """
        n = self.degree()

        if degree != 0:
            if not n % degree == 0:
                return []
            else:
                K, inc = self.subfield(degree, name=name, map=True)
                return [(K, inc)]

        divisors = n.divisors()

        if name is None:
            if hasattr(self, '_prefix'):
                name = self._prefix
            else:
                name = self.variable_name()
        if isinstance(name, str):
            name = {m: name + str(m) for m in divisors}
        elif not isinstance(name, dict):
            raise ValueError("name must be None, a string or a dictionary indexed by divisors of the degree")

        pairs = []
        for m in divisors:
            K, inc = self.subfield(m, name=name[m], map=True)
            pairs.append((K, inc))
        return pairs

    @cached_method
    def algebraic_closure(self, name='z', **kwds):
        """
        Return an algebraic closure of this field.

        INPUT:

        - ``name`` -- string (default: 'z'): prefix to use for
          variable names of subfields

        - ``implementation`` -- string (optional): specifies how to
          construct the algebraic closure.  The only value supported
          at the moment is ``'pseudo_conway'``.  For more details, see
          :mod:`~sage.rings.algebraic_closure_finite_field`.

        OUTPUT:

        An algebraic closure of ``self``.  Note that mathematically
        speaking, this is only unique up to *non-unique* isomorphism.
        To obtain canonically defined algebraic closures, one needs an
        algorithm that also provides a canonical isomorphism between
        any two algebraic closures constructed using the algorithm.

        This non-uniqueness problem can in principle be solved by
        using *Conway polynomials*; see for example
        :wikipedia:`Conway_polynomial_(finite_fields)`. These have
        the drawback that computing them takes a long time.  Therefore
        Sage implements a variant called *pseudo-Conway polynomials*,
        which are easier to compute but do not determine an algebraic
        closure up to unique isomorphism.

        The output of this method is cached, so that within the same
        Sage session, calling it multiple times will return the same
        algebraic closure (i.e. the same Sage object).  Despite this,
        the non-uniqueness of the current implementation means that
        coercion and pickling cannot work as one might expect.  See
        below for an example.

        EXAMPLES::

            sage: F = GF(5).algebraic_closure()
            sage: F
            Algebraic closure of Finite Field of size 5
            sage: F.gen(3)
            z3

        The default name is 'z' but you can change it through the option
        ``name``::

            sage: Ft = GF(5).algebraic_closure('t')
            sage: Ft.gen(3)
            t3

        Because Sage currently only implements algebraic closures
        using a non-unique definition (see above), it is currently
        impossible to implement pickling in such a way that a pickled
        and unpickled element compares equal to the original::

            sage: F = GF(7).algebraic_closure()
            sage: x = F.gen(2)
            sage: loads(dumps(x)) == x
            False

        .. NOTE::

            This is currently only implemented for prime fields.

        TESTS::

            sage: GF(5).algebraic_closure() is GF(5).algebraic_closure()
            True
        """
        from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
        return AlgebraicClosureFiniteField(self, name, **kwds)

    @cached_method
    def is_conway(self):
        """
        Return whether this field is defined by a Conway polynomial.

        EXAMPLES::

            sage: GF(5^3, 'a').is_conway()
            True
            sage: GF(5^3, 'a', modulus='adleman-lenstra').is_conway()
            False
            sage: GF(next_prime(2^16, 2), 'a').is_conway()
            False

        TESTS:

        Prime fields are defined by the Conway polynomial `x - 1`::

            sage: GF(2).is_conway()
            True

        """
        if self.base_ring() is not self.base():
            # Note that it is not entirely clear what could be implemented here:
            # * We could just use the definition of a Conway polynomial over
            #   self.base_ring() with some fixed ordering on that base. While
            #   that should provide us with a sound definition it is probably
            #   not a very useful notion.
            # * Somewhat more interesting would be the question whether the
            #   base ring is defined by a Conway polynomial and our modulus
            #   comes from the Conway polynomial of the absolute field.
            raise NotImplementedError("relative conway polynomials not supported yet")
        from .conway_polynomials import conway_polynomial, exists_conway_polynomial
        p = self.characteristic()
        n = self.degree()
        return (exists_conway_polynomial(p, n)
                and self.polynomial() == self.polynomial_ring()(conway_polynomial(p, n)))

    def relative_degree(self, base=None):
        r"""
        Return the degree of this field over ``base``.

        When base is not specified, this returns the degree over the base
        field from which this field was generated or the degree over the prime
        field if this field was not created as an extension.

        EXAMPLES::

            sage: GF(8).relative_degree()
            3
            sage: GF(4).extension(2, absolute=False).relative_degree()
            2
            sage: GF(4).extension(2, absolute=False).relative_degree(base=GF(2))
            4
            sage: GF(4).relative_degree(GF(4))
            1
            sage: GF(16).relative_degree(GF(4))
            2

        """
        if base is None:
            base = self.base_ring()
        if base is self.base_ring():
            return self.modulus().degree()
        if base is self.prime_subfield():
            return self.absolute_degree()
        if not base.is_subring(self):
            raise ValueError("base is not a subfield")
        return self.absolute_degree() // base.absolute_degree()

    degree = relative_degree

    def absolute_degree(self):
        r"""
        Return the absolute degree of this field over its prime field.

        EXAMPLES::

            sage: GF(2).absolute_degree()
            1
            sage: GF(8).absolute_degree()
            3
            sage: GF(4).extension(2, absolute=False).absolute_degree()
            4

        """
        from sage.rings.integer_ring import ZZ
        if self.base() is self:
            return ZZ(1)
        return self.relative_degree() * self.base_ring().absolute_degree()

    def frobenius_endomorphism(self, n=1):
        """
        INPUT:

        -  ``n`` -- an integer (default: 1)

        OUTPUT:

        The `n`-th power of the absolute arithmetic Frobenius
        endomorphism on this finite field.

        EXAMPLES::

            sage: k.<t> = GF(3^5)
            sage: Frob = k.frobenius_endomorphism(); Frob
            Frobenius endomorphism t |--> t^3 on Finite Field in t of size 3^5

            sage: a = k.random_element()
            sage: Frob(a) == a^3
            True

        We can specify a power::

            sage: k.frobenius_endomorphism(2)
            Frobenius endomorphism t |--> t^(3^2) on Finite Field in t of size 3^5

        The result is simplified if possible::

            sage: k.frobenius_endomorphism(6)
            Frobenius endomorphism t |--> t^3 on Finite Field in t of size 3^5
            sage: k.frobenius_endomorphism(5)
            Identity endomorphism of Finite Field in t of size 3^5

        Comparisons work::

            sage: k.frobenius_endomorphism(6) == Frob
            True
            sage: from sage.categories.morphism import IdentityMorphism
            sage: k.frobenius_endomorphism(5) == IdentityMorphism(k)
            True

        AUTHOR:

        - Xavier Caruso (2012-06-29)
        """
        if self.base_ring() is not self.base():
            # It should be straightforward to adapt
            # FrobeniusEndomorphism_finite_field to work in the general case,
            # we just have not done it yet.
            raise NotImplementedError("Frobenius endomorphism not implemented for relative extensions yet")
        from sage.rings.finite_rings.hom_finite_field import FrobeniusEndomorphism_finite_field
        return FrobeniusEndomorphism_finite_field(self, n)

    def galois_group(self):
        r"""
        Return the Galois group of this finite field, a cyclic group generated by Frobenius.

        EXAMPLES::

            sage: G = GF(3^6).galois_group()
            sage: G
            Galois group C6 of GF(3^6)
            sage: F = G.gen()
            sage: F^2
            Frob^2
            sage: F^6
            1
        """
        from sage.rings.finite_rings.galois_group import GaloisGroup_GF
        return GaloisGroup_GF(self)

    def dual_basis(self, basis=None, check=True):
        r"""
        Return the dual basis of ``basis``, or the dual basis of the power
        basis if no basis is supplied.

        If `e = \{e_0, e_1, ..., e_{n-1}\}` is a basis of
        `\GF{p^n}` as a vector space over `\GF{p}`, then the dual basis of `e`,
        `d = \{d_0, d_1, ..., d_{n-1}\}`, is the unique basis such that
        `\mathrm{Tr}(e_i d_j) = \delta_{i,j}, 0 \leq i,j \leq n-1`, where
        `\mathrm{Tr}` is the trace function.

        INPUT:

        - ``basis`` -- (default: ``None``): a basis of the finite field
          ``self``, `\GF{p^n}`, as a vector space over the base field
          `\GF{p}`. Uses the power basis `\{x^i : 0 \leq i \leq n-1\}` as
          input if no basis is supplied, where `x` is the generator of
          ``self``.

        - ``check`` -- (default: ``True``): verifies that ``basis`` is
          a valid basis of ``self``.

        ALGORITHM:

        The algorithm used to calculate the dual basis comes from pages
        110--111 of [McE1987]_.

        Let `e = \{e_0, e_1, ..., e_{n-1}\}` be a basis of `\GF{p^n}` as a
        vector space over `\GF{p}` and `d = \{d_0, d_1, ..., d_{n-1}\}` be the
        dual basis of `e`. Since `e` is a basis, we can rewrite any
        `d_c, 0 \leq c \leq n-1`, as
        `d_c = \beta_0 e_0 + \beta_1 e_1 + ... + \beta_{n-1} e_{n-1}`, for some
        `\beta_0, \beta_1, ..., \beta_{n-1} \in \GF{p}`. Using properties of
        the trace function, we can rewrite the `n` equations of the form
        `\mathrm{Tr}(e_i d_c) = \delta_{i,c}` and express the result as the
        matrix vector product:
        `A [\beta_0, \beta_1, ..., \beta_{n-1}] = i_c`, where the `i,j`-th
        element of `A` is `\mathrm{Tr(e_i e_j)}` and `i_c` is the `i`-th
        column of the `n \times n` identity matrix. Since `A` is an invertible
        matrix, `[\beta_0, \beta_1, ..., \beta_{n-1}] = A^{-1} i_c`, from
        which we can easily calculate `d_c`.

        EXAMPLES::

            sage: F.<a> = GF(2^4)
            sage: F.dual_basis(basis=None, check=False)
            [a^3 + 1, a^2, a, 1]

        We can test that the dual basis returned satisfies the defining
        property of a dual basis:
        `\mathrm{Tr}(e_i d_j) = \delta_{i,j}, 0 \leq i,j \leq n-1` ::

            sage: F.<a> = GF(7^4)
            sage: e = [4*a^3, 2*a^3 + a^2 + 3*a + 5,
            ....:      3*a^3 + 5*a^2 + 4*a + 2, 2*a^3 + 2*a^2 + 2]
            sage: d = F.dual_basis(e, check=True); d
            [3*a^3 + 4*a^2 + 6*a + 2, a^3 + 6*a + 5,
            3*a^3 + 6*a^2 + 2*a + 5, 5*a^2 + 4*a + 3]
            sage: vals = [[(x * y).trace() for x in e] for y in d]
            sage: matrix(vals) == matrix.identity(4)
            True

        We can test that if `d` is the dual basis of `e`, then `e` is the dual
        basis of `d`::

            sage: F.<a> = GF(7^8)
            sage: e = [a^0, a^1, a^2, a^3, a^4, a^5, a^6, a^7]
            sage: d = F.dual_basis(e, check=False); d
            [6*a^6 + 4*a^5 + 4*a^4 + a^3 + 6*a^2 + 3,
            6*a^7 + 4*a^6 + 4*a^5 + 2*a^4 + a^2,
            4*a^6 + 5*a^5 + 5*a^4 + 4*a^3 + 5*a^2 + a + 6,
            5*a^7 + a^6 + a^4 + 4*a^3 + 4*a^2 + 1,
            2*a^7 + 5*a^6 + a^5 + a^3 + 5*a^2 + 2*a + 4,
            a^7 + 2*a^6 + 5*a^5 + a^4 + 5*a^2 + 4*a + 4,
            a^7 + a^6 + 2*a^5 + 5*a^4 + a^3 + 4*a^2 + 4*a + 6,
            5*a^7 + a^6 + a^5 + 2*a^4 + 5*a^3 + 6*a]
            sage: F.dual_basis(d)
            [1, a, a^2, a^3, a^4, a^5, a^6, a^7]

        We cannot calculate the dual basis if ``basis`` is not a valid basis.
        ::

            sage: F.<a> = GF(2^3)
            sage: F.dual_basis([a], check=True)
            Traceback (most recent call last):
            ...
            ValueError: basis length should be 3, not 1

            sage: F.dual_basis([a^0, a, a^0 + a], check=True)
            Traceback (most recent call last):
            ...
            ValueError: value of 'basis' keyword is not a basis

        AUTHOR:

        - Thomas Gagne (2015-06-16)
        """
        if self.base_ring() is not self.base():
            raise NotImplementedError("Dual basis not implemented for relative extensions yet")

        from sage.matrix.constructor import matrix

        if basis is None:
            basis = [self.gen() ** i for i in range(self.degree())]
            check = False

        if check:
            if len(basis) != self.degree():
                msg = 'basis length should be {0}, not {1}'
                raise ValueError(msg.format(self.degree(), len(basis)))
            V = self.vector_space(map=False)
            vec_reps = [V(b) for b in basis]
            if matrix(vec_reps).is_singular():
                raise ValueError("value of 'basis' keyword is not a basis")

        entries = [(bi * bj).trace() for bi in basis for bj in basis]
        B = matrix(self.base_ring(), self.degree(), entries).inverse()
        return [sum(x * y for x, y in zip(col, basis))
                for col in B.columns()]

cdef class FiniteFieldAbsolute(FiniteField):
    r"""
    Base class for finite fields that are not relative, i.e., prime fields and
    absolute extensions thereof.

    EXAMPLES::

        sage: k = GF(4)

    TESTS::

        sage: from sage.rings.finite_rings.finite_field_base import FiniteFieldAbsolute
        sage: isinstance(k, FiniteFieldAbsolute)
        True
    """
    def absolute_gen(self):
        return self.gen()

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(GF(5), verify=True)
            # Verified
            GF(5)
            sage: sage_input(GF(32, 'a'), verify=True)
            # Verified
            R.<x> = GF(2)[]
            GF(2^5, 'a', x^5 + x^2 + 1)
            sage: K = GF(125, 'b')
            sage: sage_input((K, K), verify=True)
            # Verified
            R.<x> = GF(5)[]
            GF_5_3 = GF(5^3, 'b', x^3 + 3*x + 3)
            (GF_5_3, GF_5_3)
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: GF(81, 'a')._sage_input_(SageInputBuilder(), False)
            {call: {atomic:GF}({binop:** {atomic:3} {atomic:4}}, {atomic:'a'}, {binop:+ {binop:+ {binop:** {gen:x {constr_parent: {subscr: {call: {atomic:GF}({atomic:3})}[{atomic:'x'}]} with gens: ('x',)}} {atomic:4}} {binop:* {atomic:2} {binop:** {gen:x {constr_parent: {subscr: {call: {atomic:GF}({atomic:3})}[{atomic:'x'}]} with gens: ('x',)}} {atomic:3}}}} {atomic:2}})}

        """
        if self.base_ring() is self:
            v = sib.name('GF')(sib.int(self.characteristic()))
            name = 'GF_%d' % self.characteristic()
        elif self.base_ring() is self.base():
            v = sib.name('GF')(sib.int(self.characteristic()) ** sib.int(self.degree()),
                               self.variable_name(),
                               self.modulus())
            name = 'GF_%d_%d' % (self.characteristic(), self.degree())
        else:
            raise NotImplementedError("not implemented for this extension yet")
        sib.cache(self, v, name)
        return v

    @cached_method
    def free_module(self, base=None, basis=None, map=None, subfield=None):
        """
        Return the vector space over the subfield isomorphic to this
        finite field as a vector space, along with the isomorphisms.

        INPUT:

        - ``base`` -- a subfield of or a morphism into this finite field.
          If not given, the base ring is assumed. A subfield means
          a finite field with coercion to this finite field.

        - ``basis`` -- a basis of the finite field as a vector space
          over the subfield. If not given, one is chosen automatically.

        - ``map`` -- boolean (default: ``True``); if ``True``, isomorphisms
          from and to the vector space are also returned.

        The ``basis`` maps to the standard basis of the vector space
        by the isomorphisms.

        OUTPUT: if ``map`` is ``False``,

        - vector space over the subfield or the domain of the morphism,
          isomorphic to this finite field.

        and if ``map`` is ``True``, then also

        - an isomorphism from the vector space to the finite field.

        - the inverse isomorphism to the vector space from the finite field.

        EXAMPLES::

            sage: GF(27,'a').vector_space(map=False)
            Vector space of dimension 3 over Finite Field of size 3

            sage: F = GF(8)
            sage: E = GF(64)
            sage: V, from_V, to_V = E.vector_space(F, map=True)
            sage: V
            Vector space of dimension 2 over Finite Field in z3 of size 2^3
            sage: to_V(E.gen())
            (0, 1)
            sage: all(from_V(to_V(e)) == e for e in E)
            True
            sage: all(to_V(e1 + e2) == to_V(e1) + to_V(e2) for e1 in E for e2 in E)
            True
            sage: all(to_V(c * e) == c * to_V(e) for e in E for c in F)
            True

            sage: basis = [E.gen(), E.gen() + 1]
            sage: W, from_W, to_W = E.vector_space(F, basis, map=True)
            sage: all(from_W(to_W(e)) == e for e in E)
            True
            sage: all(to_W(c * e) == c * to_W(e) for e in E for c in F)
            True
            sage: all(to_W(e1 + e2) == to_W(e1) + to_W(e2) for e1 in E for e2 in E)  # long time
            True
            sage: to_W(basis[0]); to_W(basis[1])
            (1, 0)
            (0, 1)

            sage: F = GF(9, 't', modulus=(x^2+x-1))
            sage: E = GF(81)
            sage: h = Hom(F,E).an_element()
            sage: V, from_V, to_V = E.vector_space(h, map=True)
            sage: V
            Vector space of dimension 2 over Finite Field in t of size 3^2
            sage: V.base_ring() is F
            True
            sage: all(from_V(to_V(e)) == e for e in E)
            True
            sage: all(to_V(e1 + e2) == to_V(e1) + to_V(e2) for e1 in E for e2 in E)
            True
            sage: all(to_V(h(c) * e) == c * to_V(e) for e in E for c in F)
            True

        In a tower of finite fields::

            sage: k.<a> = GF(9)
            sage: l.<b> = k.extension(3, absolute=False)
            sage: m.<c> = l.extension(4, absolute=False)

            sage: m.vector_space(base=m, map=False)
            Vector space of dimension 1 over Finite Field in c of size 3^24 over its base
            sage: m.vector_space(base=l, map=False)
            Vector space of dimension 4 over Finite Field in b of size 3^6 over its base
            sage: m.vector_space(base=k, map=False)
            Vector space of dimension 12 over Finite Field in a of size 3^2
            sage: m.vector_space(base=m.prime_subfield(), map=False)
            

        """
        from sage.modules.all import VectorSpace
        from sage.categories.morphism import is_Morphism
        if subfield is not None:
            if base is not None:
                raise ValueError
            deprecation(28481, "The subfield keyword argument has been renamed to base")
            base = subfield
        if map is None:
            deprecation(28481, "The default value for map will be changing to True.  To keep the current behavior, explicitly pass map=False.")
            map = False

        if base is None:
            base = self.base_ring()
            s = self.relative_degree()
            if not hasattr(self, '_vector_space'):
                self._vector_space = VectorSpace(base, s)
            V = self._vector_space
            inclusion_map = None
        elif is_Morphism(base):
            inclusion_map = base
            base = inclusion_map.domain()
            s = self.absolute_degree() // base.absolute_degree()
            V = VectorSpace(base, s)
        elif base.is_subring(self):
            s = self.absolute_degree() // base.absolute_degree()
            V = VectorSpace(base, s)
            inclusion_map = None
        else:
            raise ValueError("{} is not a subfield".format(base))

        if map is False: # shortcut
            return V

        if inclusion_map is None:
            inclusion_map = self.coerce_map_from(base)

        from sage.modules.all import vector
        from sage.matrix.all import matrix
        from .maps_finite_field import (
            MorphismVectorSpaceToFiniteField, MorphismFiniteFieldToVectorSpace)

        E = self
        F = base

        alpha = E.gen()
        beta = F.gen()

        if basis is None:
            basis = [alpha**i for i in range(s)] # of E over F

        F_basis = [beta**i for i in range(F.degree())]

        # E_basis_alpha is the implicit basis of E over the prime subfield
        E_basis_beta = [inclusion_map(F_basis[i]) * basis[j]
                        for j in range(s)
                        for i in range(F.absolute_degree())]

        C = matrix(E.prime_subfield(), E.absolute_degree(), E.absolute_degree(),
                   [E_basis_beta[i]._vector_() for i in range(E.absolute_degree())])
        C.set_immutable()
        Cinv = C.inverse()
        Cinv.set_immutable()

        phi = MorphismVectorSpaceToFiniteField(V, self, C)
        psi = MorphismFiniteFieldToVectorSpace(self, V, Cinv)

        return V, phi, psi



def unpickle_FiniteField_ext(_type, order, variable_name, modulus, kwargs):
    r"""
    Used to unpickle extensions of finite fields. Now superseded (hence no
    doctest), but kept around for backward compatibility.
    """
    return _type(order, variable_name, modulus, **kwargs)


def unpickle_FiniteField_prm(_type, order, variable_name, kwargs):
    r"""
    Used to unpickle finite prime fields. Now superseded (hence no doctest),
    but kept around for backward compatibility.
    """
    return _type(order, variable_name, **kwargs)

register_unpickle_override(
    'sage.rings.ring', 'unpickle_FiniteField_prm', unpickle_FiniteField_prm)


def is_FiniteField(R):
    r"""
    Return whether the implementation of ``R`` has the interface provided by
    the standard finite field implementation.

    EXAMPLES::

        sage: from sage.rings.finite_rings.finite_field_base import is_FiniteField
        sage: is_FiniteField(GF(9,'a'))
        True
        sage: is_FiniteField(GF(next_prime(10^10)))
        True

    Note that the integers modulo n are not backed by the finite field type::

        sage: is_FiniteField(Integers(7))
        False
    """
    return isinstance(R, FiniteField)
