"""
Base Classes for Finite Fields

TESTS::

    sage: K.<a> = NumberField(x^2 + 1)
    sage: F = K.factor(3)[0][0].residue_field()
    sage: loads(dumps(F)) == F
    True

AUTHORS:

- Adrien Brochard, David Roe, Jeroen Demeyer, Julian Rueth, Niles Johnson,
  Peter Bruin, Travis Scrimshaw, Xavier Caruso: initial version

"""
#*****************************************************************************
#       Copyright (C) 2009 David Roe <roed@math.harvard.edu>
#       Copyright (C) 2010 Niles Johnson <nilesj@gmail.com>
#       Copyright (C) 2011 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#       Copyright (C) 2012 Adrien Brochard <aaa.brochard@gmail.com>
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
#       Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#       Copyright (C) 2013 Peter Bruin <P.Bruin@warwick.ac.uk>
#       Copyright (C) 2014 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.finite_fields import FiniteFields
from sage.structure.parent cimport Parent
from sage.misc.cachefunc import cached_method
from sage.misc.prandom import randrange

# Copied from sage.misc.fast_methods, used in __hash__() below.
cdef int SIZEOF_VOID_P_SHIFT = 8*sizeof(void *) - 4

cdef class FiniteFieldIterator:
    r"""
    An iterator over a finite field. This should only be used when the field
    is an extension of a smaller field which already has a separate iterator
    function.
    """
    cdef object iter
    cdef FiniteField parent

    def __init__(self,FiniteField parent):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: k = iter(FiniteField(9, 'a', impl='pari_ffelt')) # indirect doctest
            sage: isinstance(k, sage.rings.finite_rings.finite_field_base.FiniteFieldIterator)
            True
            sage: k = iter(FiniteField(16, 'a', impl='ntl')) # indirect doctest
            sage: isinstance(k, sage.rings.finite_rings.finite_field_base.FiniteFieldIterator)
            True
        """
        self.parent = parent
        self.iter = iter(self.parent.vector_space())

    def __next__(self):
        r"""
        Return the next element in the iterator.

        EXAMPLE::

            sage: k = iter(FiniteField(9, 'a', impl='pari_ffelt'))
            sage: next(k) # indirect doctest
            0
        """
        return self.parent(next(self.iter))

    def __iter__(self):
        """
        Return ``self`` since this is an iterator class.

        EXAMPLES::

            sage: K.<a> = GF(7^4)
            sage: K.list()[:7]
            [0, a, a^2, a^3, 2*a^2 + 3*a + 4, 2*a^3 + 3*a^2 + 4*a, 3*a^3 + a^2 + 6*a + 1]
            sage: K.<a> = GF(5^9)
            sage: for x in K:
            ....:     if x == a+3: break
            ....:     print(x)
            0
            1
            2
            3
            4
            a
            a + 1
            a + 2
        """
        return self


cdef class FiniteField(Field):
    """
    Abstract base class for finite fields.
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
        The hash provided by this class coincides with that of ``<type 'object'>``.

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

    def __repr__(self):
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
        if self.degree()>1:
            return "Finite Field in %s of size %s^%s"%(self.variable_name(),self.characteristic(),self.degree())
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
        if self.degree() > 1:
            e = "^{%s}"%self.degree()
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
        if self.degree() == 1:
            return 'GF(%s)'%self.order()
        B = self.base_ring()
        p = self.polynomial()
        s = "ext<%s|%s>"%(B._magma_init_(magma),p._magma_init_(magma))
        return magma._with_names(s, self.variable_names())

    def _macaulay2_init_(self):
        """
        Returns the string representation of ``self`` that Macaulay2 can
        understand.

        EXAMPLES::

            sage: GF(97,'a')._macaulay2_init_()
            'GF 97'

            sage: macaulay2(GF(97, 'a'))       # optional - macaulay2
            GF 97
            sage: macaulay2(GF(49, 'a'))       # optional - macaulay2
            GF 49
        """
        return "GF %s"%(self.order())

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
        if self.degree() == 1:
            v = sib.name('GF')(sib.int(self.characteristic()))
            name = 'GF_%d' % self.characteristic()
        else:
            v = sib.name('GF')(sib.int(self.characteristic()) ** sib.int(self.degree()),
                               self.variable_name(),
                               self.modulus())
            name = 'GF_%d_%d' % (self.characteristic(), self.degree())
        sib.cache(self, v, name)
        return v

    def __iter__(self):
        """
        Return an iterator over the elements of this finite field. This generic
        implementation uses the fairly simple :class:`FiniteFieldIterator`
        class; derived classes may implement their own more sophisticated
        replacements.

        EXAMPLES::

            sage: k = FiniteField(8, 'a', impl='pari_ffelt')
            sage: i = iter(k); i # indirect doctest
            <sage.rings.finite_rings.finite_field_base.FiniteFieldIterator object at ...>
            sage: next(i)
            0
            sage: list(k) # indirect doctest
            [0, 1, a, a + 1, a^2, a^2 + 1, a^2 + a, a^2 + a + 1]
        """
        return FiniteFieldIterator(self)

    def _is_valid_homomorphism_(self, codomain, im_gens):
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
        if self.characteristic() != codomain.characteristic():
            raise ValueError("no map from %s to %s" % (self, codomain))
        if len(im_gens) != 1:
            raise ValueError("only one generator for finite fields")

        return self.modulus()(im_gens[0]).is_zero()

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

        ALGORITHM; [Coh]_, algorithm 3.4.2 which is basically the algorithm in
        [Yun]_ with special treatment for powers divisible by `p`.

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

        REFERENCES:

        .. [Coh] \H. Cohen, A Course in Computational Algebraic Number
           Theory.  Springer-Verlag, 1993.

        .. [Yun] Yun, David YY. On square-free decomposition algorithms.
           In Proceedings of the third ACM symposium on Symbolic and algebraic
           computation, pp. 26-35. ACM, 1976.

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

    def gen(self):
        r"""
        Return a generator of this field (over its prime field). As this is an
        abstract base class, this just raises a ``NotImplementedError``.

        EXAMPLES::

            sage: K = GF(17)
            sage: sage.rings.finite_rings.finite_field_base.FiniteField.gen(K)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

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
        Returns an element of multiplicative order ``n`` in this
        finite field, if there is one.  Raises a ``ValueError`` if there
        is not.

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

        Even more examples::

            sage: GF(9,'a').zeta_order()
            8
            sage: GF(9,'a').zeta()
            a
            sage: GF(9,'a').zeta(4)
            a + 1
            sage: GF(9,'a').zeta()^2
            a + 1
        """
        z = self.multiplicative_generator()
        if n is None:
            return z
        else:
            import sage.rings.integer
            n = sage.rings.integer.Integer(n)
            m = z.multiplicative_order()
            if m % n != 0:
                raise ValueError("No %sth root of unity in self" % n)
            return z**(m // n)

    def multiplicative_generator(self):
        """
        Return a primitive element of this finite field, i.e. a generator
        of the multiplicative group.

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
        from sage.arith.all import primitive_root

        if self.__multiplicative_generator is not None:
            return self.__multiplicative_generator
        if self.degree() == 1:
            self.__multiplicative_generator = self(primitive_root(self.order()))
            return self.__multiplicative_generator
        n = self.order() - 1
        g = self.gen(0)
        # We check whether x+g is a multiplicative generator, where
        # x runs through the finite field.
        # This has the advantage that g is the first element we try,
        # so we always get g as generator if possible.  Second, the
        # PARI finite field iterator gives all the constant elements
        # first, so we try g+(constant) before anything else.
        for x in self:
            a = g+x
            if a != 0 and a.multiplicative_order() == n:
                self.__multiplicative_generator = a
                return a

    primitive_element = multiplicative_generator

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
        Returns whether or not the finite field is a field, i.e.,
        always returns ``True``.

        EXAMPLES::

            sage: k.<a> = FiniteField(3^4)
            sage: k.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return ``True`` since a finite field is finite.

        EXAMPLES::

            sage: GF(997).is_finite()
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
        return self.characteristic()**self.degree()

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
        return Factorization([(self.characteristic(), self.degree())])

    def factored_unit_order(self):
        """
        Returns the factorization of ``self.order()-1``, as a 1-element list.

        The format is for compatibility with
        :mod:`~sage.rings.finite_rings.integer_mod_ring`.

        EXAMPLES::

            sage: GF(7^2,'a').factored_unit_order()
            [2^4 * 3]
        """
        if self.__factored_unit_order is None:
            self.__factored_unit_order = [(self.order()-1).factor()]
        return self.__factored_unit_order

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
        return self.degree()==1

    def modulus(self):
        r"""
        Return the minimal polynomial of the generator of ``self`` over
        the prime finite field.

        The minimal polynomial of an element `a` in a field is the
        unique monic irreducible polynomial of smallest degree with
        coefficients in the base field that has `a` as a root. In
        finite field extensions, `\GF{p^n}`, the base field is `\GF{p}`.

        OUTPUT:

        - a monic polynomial over `\GF{p}` in the variable `x`.

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
        """
        # Normally, this is set by the constructor of the implementation
        try:
            return self._modulus
        except AttributeError:
            pass

        from sage.rings.all import PolynomialRing
        from finite_field_constructor import GF
        R = PolynomialRing(GF(self.characteristic()), 'x')
        self._modulus = R((-1,1))  # Polynomial x - 1
        return self._modulus

    def polynomial(self, name=None):
        """
        Return the minimal polynomial of the generator of ``self`` over
        the prime finite field.

        INPUT:

        - ``name`` -- a variable name to use for the polynomial. By
          default, use the name given when constructing this field.

        OUTPUT:

        - a monic polynomial over `\GF{p}` in the variable ``name``.

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
            sage: k.random_element()
            a^3 + 3*a^2 + 6*a + 9

        Passes extra positional or keyword arguments through::

            sage: k.random_element(prob=0)
            0

        """
        if self.degree() == 1:
            return self(randrange(self.order()))
        v = self.vector_space().random_element(*args, **kwds)
        return self(v)

    def some_elements(self):
        """
        Returns a collection of elements of this finite field for use in unit
        testing.

        EXAMPLES::

            sage: k = GF(2^8,'a')
            sage: k.some_elements() # random output
            [a^4 + a^3 + 1, a^6 + a^4 + a^3, a^5 + a^4 + a, a^2 + a]
        """
        return [self.random_element() for i in range(4)]

    def polynomial_ring(self, variable_name=None):
        """
        Returns the polynomial ring over the prime subfield in the
        same variable as this finite field.

        EXAMPLES::

            sage: k.<alpha> = FiniteField(3^4)
            sage: k.polynomial_ring()
            Univariate Polynomial Ring in alpha over Finite Field of size 3
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.finite_rings.finite_field_constructor import GF

        if variable_name is None and self.__polynomial_ring is not None:
            return self.__polynomial_ring
        else:
            if variable_name is None:
                self.__polynomial_ring = PolynomialRing(GF(self.characteristic()), self.variable_name())
                return self.__polynomial_ring
            else:
                return PolynomialRing(GF(self.characteristic()), variable_name)

    def vector_space(self):
        """
        Return the vector space over the prime subfield isomorphic
        to this finite field as a vector space.

        EXAMPLES::

            sage: GF(27,'a').vector_space()
            Vector space of dimension 3 over Finite Field of size 3
        """
        if self.__vector_space is not None:
            return self.__vector_space
        else:
            import sage.modules.all
            V = sage.modules.all.VectorSpace(self.prime_subfield(),self.degree())
            self.__vector_space = V
            return V

    cpdef _coerce_map_from_(self, R):
        r"""
        Canonical coercion to ``self``.

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
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.finite_rings.finite_field_base import is_FiniteField
        from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing
        if R is int or R is long or R is ZZ:
            return True
        if is_IntegerModRing(R) and self.characteristic().divides(R.characteristic()):
            return R.hom((self.one(),), check=False)
        if is_FiniteField(R):
            if R is self:
                return True
            from residue_field import ResidueField_generic
            if isinstance(R, ResidueField_generic):
                return False
            if R.characteristic() == self.characteristic():
                if R.degree() == 1:
                    return R.hom((self.one(),), check=False)
                elif (R.degree().divides(self.degree())
                      and hasattr(self, '_prefix') and hasattr(R, '_prefix')):
                    return R.hom((self.gen() ** ((self.order() - 1)//(R.order() - 1)),))

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
        """
        from sage.categories.pushout import AlgebraicExtensionFunctor
        if self.degree() == 1:
            # this is not of type FiniteField_prime_modn
            from sage.rings.integer import Integer
            return AlgebraicExtensionFunctor([self.polynomial()], [None], [None]), self.base_ring()
        elif hasattr(self, '_prefix'):
            return (AlgebraicExtensionFunctor([self.degree()], [self.variable_name()], [None],
                                              prefix=self._prefix),
                    self.base_ring())
        else:
            return (AlgebraicExtensionFunctor([self.polynomial()], [self.variable_name()], [None]),
                    self.base_ring())

    def extension(self, modulus, name=None, names=None, map=False, embedding=None, **kwds):
        """
        Return an extension of this finite field.

        INPUT:

        - ``modulus`` -- a polynomial with coefficients in ``self``,
          or an integer.

        - ``name`` -- string: the name of the generator in the new
          extension

        - ``map`` -- boolean (default: ``False``): if ``False``,
          return just the extension `E`; if ``True``, return a pair
          `(E, f)`, where `f` is an embedding of ``self`` into `E`.

        - ``embedding`` -- currently not used; for compatibility with
          other ``AlgebraicExtensionFunctor`` calls.

        - ``**kwds``: further keywords, passed to the finite field
          constructor.

        OUTPUT:

        An extension of the given modulus, or pseudo-Conway of the
        given degree if ``modulus`` is an integer.

        EXAMPLES::

            sage: k = GF(2)
            sage: R.<x> = k[]
            sage: k.extension(x^1000 + x^5 + x^4 + x^3 + 1, 'a')
            Finite Field in a of size 2^1000
            sage: k = GF(3^4)
            sage: R.<x> = k[]
            sage: k.extension(3)
            Finite Field in z12 of size 3^12

        An example using the ``map`` argument::

            sage: F = GF(5)
            sage: E, f = F.extension(2, 'b', map=True)
            sage: E
            Finite Field in b of size 5^2
            sage: f
            Ring morphism:
              From: Finite Field of size 5
              To:   Finite Field in b of size 5^2
              Defn: 1 |--> 1
            sage: f.parent()
            Set of field embeddings from Finite Field of size 5 to Finite Field in b of size 5^2

        Extensions of non-prime finite fields by polynomials are not yet
        supported: we fall back to generic code::

            sage: k.extension(x^5 + x^2 + x - 1)
            Univariate Quotient Polynomial Ring in x over Finite Field in z4 of size 3^4 with modulus x^5 + x^2 + x + 2

        TESTS:

        We check that :trac:`18915` is fixed::

            sage: F = GF(2)
            sage: F.extension(int(3), 'a')
            Finite Field in a of size 2^3

            sage: F = GF(2 ** 4, 'a')
            sage: F.extension(int(3), 'aa')
            Finite Field in aa of size 2^12
        """
        from finite_field_constructor import GF
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        from sage.rings.integer import Integer
        if name is None and names is not None:
            name = names
        if self.degree() == 1:
            if isinstance(modulus, (int, Integer)):
                E = GF(self.characteristic()**modulus, name=name, **kwds)
            elif isinstance(modulus, (list, tuple)):
                E = GF(self.characteristic()**(len(modulus) - 1), name=name, modulus=modulus, **kwds)
            elif is_Polynomial(modulus):
                if modulus.change_ring(self).is_irreducible():
                    E = GF(self.characteristic()**(modulus.degree()), name=name, modulus=modulus, **kwds)
                else:
                    E = Field.extension(self, modulus, name=name, embedding=embedding)
        elif isinstance(modulus, (int, Integer)):
            E = GF(self.order()**modulus, name=name, **kwds)
        else:
            E = Field.extension(self, modulus, name=name, embedding=embedding)
        if not map:
            return E
        # Use the canonical map if it exists.
        f = E.coerce_map_from(self)
        if f is None:
            from sage.categories.homset import Hom
            f = Hom(self, E).an_element()
        return (E, f)

    def subfields(self, degree=0, name=None):
        """
        Return all subfields of ``self`` of the given ``degree``,
        or all possible degrees if ``degree`` is `0`.

        The subfields are returned as absolute fields together with
        an embedding into ``self``.

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
        this field and ``e`` gives an embedding of ``K`` into ``self``.

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
              Ring morphism:
                  From: Finite Field in z21 of size 2^21
                  To:   Finite Field in z21 of size 2^21
                  Defn: z21 |--> z21)]
        """
        from sage.rings.integer import Integer
        from finite_field_constructor import GF
        p = self.characteristic()
        n = self.degree()
        if degree != 0:
            degree = Integer(degree)
            if not degree.divides(n):
                return []
            elif hasattr(self, '_prefix'):
                K = GF(p**degree, name=name, prefix=self._prefix)
                return [(K, self.coerce_map_from(K))]
            elif degree == 1:
                K = GF(p)
                return [(K, self.coerce_map_from(K))]
            else:
                gen = self.gen()**((self.order() - 1)//(p**degree - 1))
                K = GF(p**degree, modulus=gen.minimal_polynomial(), name=name)
                return [(K, K.hom((gen,)))]
        else:
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
            return [self.subfields(m, name=name[m])[0] for m in divisors]

    @cached_method
    def algebraic_closure(self, name='z', **kwds):
        """
        Return an algebraic closure of ``self``.

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
        using *Conway polynomials*; see for example [CP]_.  These have
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

        EXAMPLE::

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

        REFERENCE:

        .. [CP] Wikipedia entry on Conway polynomials,
           :wikipedia:`Conway_polynomial_(finite_fields)`

        TEST::

            sage: GF(5).algebraic_closure() is GF(5).algebraic_closure()
            True

        """
        from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
        return AlgebraicClosureFiniteField(self, name, **kwds)

    @cached_method
    def is_conway(self):
        """
        Return ``True`` if self is defined by a Conway polynomial.

        EXAMPLES:

            sage: GF(5^3, 'a').is_conway()
            True
            sage: GF(5^3, 'a', modulus='adleman-lenstra').is_conway()
            False
            sage: GF(next_prime(2^16, 2), 'a').is_conway()
            False
        """
        from conway_polynomials import conway_polynomial, exists_conway_polynomial
        p = self.characteristic()
        n = self.degree()
        return (exists_conway_polynomial(p, n)
                and self.polynomial() == self.polynomial_ring()(conway_polynomial(p, n)))

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
        from sage.rings.finite_rings.hom_finite_field import FrobeniusEndomorphism_finite_field
        return FrobeniusEndomorphism_finite_field(self, n)

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
        110--111 of [FFCSE1987]_.

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

        REFERENCES:

        .. [FFCSE1987] Robert J. McEliece. Finite Fields for Computer
           Scientists and Engineers. Kluwer Academic Publishers, 1987.

        AUTHOR:

        - Thomas Gagne (2015-06-16)
        """
        from sage.matrix.constructor import matrix

        if basis == None:
            basis = [self.gen()**i for i in range(self.degree())]
            check = False

        if check:
            if len(basis) != self.degree():
                msg = 'basis length should be {0}, not {1}'
                raise ValueError(msg.format(self.degree(), len(basis)))
            V = self.vector_space()
            vec_reps = [V(b) for b in basis]
            if matrix(vec_reps).is_singular():
                raise ValueError('value of \'basis\' keyword is not a basis')

        entries = [(basis[i] * basis[j]).trace() for i in range(self.degree())
                    for j in range(self.degree())]
        B = matrix(self.base_ring(), self.degree(), entries).inverse()
        db = [sum(map(lambda x: x[0] * x[1], zip(col, basis)))
              for col in B.columns()]
        return db

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


def is_FiniteField(x):
    """
    Return ``True`` if ``x`` is of type finite field, and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.rings.finite_rings.finite_field_base import is_FiniteField
        sage: is_FiniteField(GF(9,'a'))
        True
        sage: is_FiniteField(GF(next_prime(10^10)))
        True

    Note that the integers modulo n are not of type finite field,
    so this function returns ``False``::

        sage: is_FiniteField(Integers(7))
        False
    """
    return isinstance(x, FiniteField)
