# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.rings.number_field
"""
Base class of number fields

AUTHORS:

- William Stein (2007-09-04): initial version
"""


def is_NumberField(x):
    """
    Return ``True`` if ``x`` is of number field type.

    This function is deprecated.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_base import is_NumberField
        sage: x = polygen(ZZ)
        sage: is_NumberField(NumberField(x^2 + 1, 'a'))
        doctest:...: DeprecationWarning: the function is_NumberField is deprecated; use
        isinstance(x, sage.rings.number_field.number_field_base.NumberField) instead
        See https://github.com/sagemath/sage/issues/35283 for details.
        True
        sage: is_NumberField(QuadraticField(-97, 'theta'))
        True
        sage: is_NumberField(CyclotomicField(97))
        True

    Note that the rational numbers ``QQ`` are a number field.::

        sage: is_NumberField(QQ)
        True
        sage: is_NumberField(ZZ)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(35283,
                "the function is_NumberField is deprecated; use "
                "isinstance(x, sage.rings.number_field.number_field_base.NumberField) instead")

    return isinstance(x, NumberField)


from sage.rings.ring cimport Field


cdef class NumberField(Field):
    r"""
    Base class for all number fields.

    TESTS::

        sage: x = polygen(ZZ)
        sage: k = NumberField(x^2 + 1, 'i')
        sage: k == loads(dumps(k))
        True

        sage: z = polygen(QQ)
        sage: K.<theta, beta> = NumberField([z^3 - 3, z^2 + 1])
        sage: K.is_finite()
        False
        sage: K.order()
        +Infinity
    """
    # This token docstring is mostly there to prevent Sphinx from pasting in
    # the docstring of the __init__ method inherited from IntegralDomain, which
    # is rather confusing.
    def _pushout_(self, other):
        r"""
        If ``self`` and/or ``other`` are embedded, use this embedding to
        discover a common parent.

        Currently embeddings into ``AA`` and ``QQbar`` are supported.

        TESTS:

        Pushout is implemented for number field embedded in ``AA``::

            sage: x = polygen(ZZ)
            sage: K.<a> = NumberField(x^2 - 3, embedding=AA(3)**(1/2))
            sage: L.<b> = NumberField(x^2 - 2, embedding=AA(2)**(1/2))
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.explain(K,L,operator.add)
            Coercion on left operand via
                Generic morphism:
                  From: Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
                  To:   Algebraic Real Field
                  Defn: a -> 1.732050807568878?
            Coercion on right operand via
                Generic morphism:
                  From: Number Field in b with defining polynomial x^2 - 2 with b = 1.414213562373095?
                  To:   Algebraic Real Field
                  Defn: b -> 1.414213562373095?
            Arithmetic performed after coercions.
            Result lives in Algebraic Real Field
            Algebraic Real Field

        As a consequence, operations and comparisons work nicely::

            sage: a + b
            3.146264369941973?
            sage: a < b
            False
            sage: 3*a < 4*b
            True

        Using number field with other classes::

            sage: K.<cbrt2> = NumberField(x^3 - 2, embedding=AA(2)**(1/3))
            sage: (cbrt2 + a) * b
            4.231287179063857?
            sage: b + QQbar(-3).sqrt()
            1.414213562373095? + 1.732050807568878?*I

        Pushout is implemented for number field embedded in ``QQbar``::

            sage: Km2.<sqrtm2> = NumberField(x^2 + 2, embedding=QQbar(-2).sqrt())
            sage: b + sqrtm2
            1.414213562373095? + 1.414213562373095?*I
            sage: sqrtm2 + b
            1.414213562373095? + 1.414213562373095?*I
            sage: sqrtm2 + AA(3).sqrt()
            1.732050807568878? + 1.414213562373095?*I
        """
        # Use the embedding of ``self``, if it exists.
        if self._embedding:
            codomain_self = self._embedding.codomain()
        else:
            codomain_self = self

        # Use the embedding of ``other``, if it exists.
        if isinstance(other, NumberField):
            embedding = (<NumberField>other)._embedding
            if embedding:
                codomain_other = embedding.codomain()
            else:
                codomain_other = other
        else:
            codomain_other = other

        try:
            from sage.rings.qqbar import AA, QQbar
        except ImportError:
            pass
        else:
            if codomain_self is AA and codomain_other is AA:
                return AA
            if codomain_self in (AA, QQbar) and codomain_other in (AA, QQbar):
                return QQbar

    def ring_of_integers(self, *args, **kwds):
        r"""
        Synonym for :meth:`maximal_order`.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.ring_of_integers()
            Gaussian Integers generated by a in Number Field in a with defining polynomial x^2 + 1
        """
        return self.maximal_order(*args, **kwds)

    def OK(self, *args, **kwds):
        r"""
        Synonym for :meth:`maximal_order`.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: NumberField(x^3 - 2,'a').OK()
            Maximal Order generated by a in Number Field in a with defining polynomial x^3 - 2
        """
        return self.maximal_order(*args, **kwds)

    def maximal_order(self):
        r"""
        Return the maximal order, i.e., the ring of integers of this
        number field.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: NumberField(x^3 - 2,'b').maximal_order()
            Maximal Order generated by b in Number Field in b with defining polynomial x^3 - 2
        """
        raise NotImplementedError

    def is_absolute(self):
        r"""
        Return ``True`` if ``self`` is viewed as a single extension over `\QQ`.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: K.<a> = NumberField(x^3 + 2)
            sage: K.is_absolute()
            True
            sage: y = polygen(K)
            sage: L.<b> = NumberField(y^2 + 1)
            sage: L.is_absolute()
            False
            sage: QQ.is_absolute()
            True
        """
        raise NotImplementedError

    def signature(self):
        r"""
        Return `(r_1, r_2)`, where `r_1` and `r_2` are the number of real embeddings
        and pairs of complex embeddings of this field, respectively.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: NumberField(x^3 - 2, 'a').signature()
            (1, 1)
        """
        raise NotImplementedError

    def degree(self):
        """
        Return the degree of this number field.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: NumberField(x^3 + 9, 'a').degree()
            3
        """
        raise NotImplementedError

    def discriminant(self):
        """
        Return the discriminant of this number field.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: NumberField(x^3 + 9, 'a').discriminant()
            -243
        """
        raise NotImplementedError

    def minkowski_bound(self):
        r"""
        Return the Minkowski bound associated to this number field.

        This is a bound `B` so that every integral ideal is equivalent
        modulo principal fractional ideals to an integral ideal of
        norm at most `B`.

        .. SEEALSO::

            :meth:`~bach_bound`

        OUTPUT: symbolic expression or Rational

        EXAMPLES:

        The Minkowski bound for `\QQ[i]` tells us that the class
        number is 1::

            sage: # needs sage.symbolic
            sage: K = QQ[I]
            sage: B = K.minkowski_bound(); B
            4/pi
            sage: B.n()
            1.27323954473516

        We compute the Minkowski bound for `\QQ[\sqrt[3]{2}]`::

            sage: # needs sage.symbolic
            sage: K = QQ[2^(1/3)]
            sage: B = K.minkowski_bound(); B
            16/3*sqrt(3)/pi
            sage: B.n()
            2.94042077558289
            sage: int(B)
            2

        We compute the Minkowski bound for `\QQ[\sqrt{10}]`, which has class
        number 2::

            sage: # needs sage.symbolic
            sage: K = QQ[sqrt(10)]
            sage: B = K.minkowski_bound(); B
            sqrt(10)
            sage: int(B)
            3
            sage: K.class_number()
            2

        We compute the Minkowski bound for `\QQ[\sqrt{2}+\sqrt{3}]`::

            sage: # needs sage.symbolic
            sage: x = polygen(ZZ)
            sage: K.<y,z> = NumberField([x^2 - 2, x^2 - 3])
            sage: L.<w> = QQ[sqrt(2) + sqrt(3)]
            sage: B = K.minkowski_bound(); B
            9/2
            sage: int(B)
            4
            sage: B == L.minkowski_bound()
            True
            sage: K.class_number()
            1

        The bound of course also works for the rational numbers::

            sage: QQ.minkowski_bound()
            1
        """
        _, s = self.signature()
        n = self.absolute_degree()
        d = self.absolute_discriminant().abs().sqrt()
        from sage.symbolic.constants import pi
        if s > 0:
            return d * (4/pi)**s * n.factorial() / (n**n)
        else:
            return d * n.factorial() / (n**n)

    def bach_bound(self):
        r"""
        Return the Bach bound associated to this number field.

        Assuming the General Riemann Hypothesis, this is a bound `B` so
        that every integral ideal is equivalent modulo principal
        fractional ideals to an integral ideal of norm at most `B`.

        .. SEEALSO::

            :meth:`~minkowski_bound`

        OUTPUT: symbolic expression or the Integer 1

        EXAMPLES:

        We compute both the Minkowski and Bach bounds for a quadratic
        field, where the Minkowski bound is much better::

            sage: # needs sage.symbolic
            sage: K = QQ[sqrt(5)]
            sage: K.minkowski_bound()
            1/2*sqrt(5)
            sage: K.minkowski_bound().n()
            1.11803398874989
            sage: K.bach_bound()
            12*log(5)^2
            sage: K.bach_bound().n()
            31.0834847277628

        We compute both the Minkowski and Bach bounds for a bigger
        degree field, where the Bach bound is much better::

            sage: # needs sage.symbolic
            sage: K = CyclotomicField(37)
            sage: K.minkowski_bound().n()
            7.50857335698544e14
            sage: K.bach_bound().n()
            191669.304126267

        The bound of course also works for the rational numbers::

            sage: QQ.bach_bound()                                                       # needs sage.symbolic
            1
        """
        ans = 12 * abs(self.discriminant()).log()**2
        if ans == 0: # rational numbers
            from sage.rings.integer import Integer
            return Integer(1)
        return ans

    # Approximate embeddings for comparisons with respect to the order of RR or
    # CC

    def _init_embedding_approx(self):
        r"""
        Initialize the approximation of embeddings.

        This should be called only once.

        TESTS::

            sage: x = polygen(ZZ)
            sage: K.<a> = NumberField(x^3 - x^2 - x - 1, embedding=1)
            sage: K._get_embedding_approx(0)   # indirect doctest
            1.839286755214161?
        """

        if self._gen_approx is not None or self._embedding is None:
            return

        try:
            from sage.rings.qqbar import AA
        except ImportError:
            AA = None

        try:
            from sage.rings.real_lazy import RLF
        except ImportError:
            RLF = None

        codomain = self._embedding.codomain()
        if codomain is AA or codomain is RLF:
            self._gen_approx = []
            self._embedded_real = 1

    cpdef _get_embedding_approx(self, size_t i):
        r"""
        Return an interval approximation of the generator of this number field.

        OUTPUT: a real interval element with precision `53 \times 2^i`

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: p = x^5 - 3*x + 1
            sage: a_AA = AA.polynomial_root(p, RIF(0,1))
            sage: K.<a> = NumberField(p, embedding=a_AA)
            sage: K._get_embedding_approx(2)
            0.3347341419433526870750989624732833071257517550374285560578335863?
            sage: K._get_embedding_approx(1)
            0.33473414194335268707509896247329?
            sage: K._get_embedding_approx(1).str(style='brackets')
            '[0.334734141943352687075098962473280 .. 0.334734141943352687075098962473287]'

            sage: K._get_embedding_approx(2).prec()
            212
            sage: K._get_embedding_approx(1).prec()
            106
            sage: K._get_embedding_approx(0).prec()
            53

        If a real embedding is not specified, this method will result in an error::

            sage: N.<g> = NumberField(x^3 + 2)
            sage: N._get_embedding_approx(1)
            Traceback (most recent call last):
            ...
            ValueError: No embedding set. You need to specify a real embedding.


        .. SEEALSO::

            :class:` RealIntervalField_class <sage.rings.real_mpfi.RealIntervalField_class>`
        """
        if self._embedded_real and i < len(self._gen_approx):
            return self._gen_approx[i]

        cdef size_t j
        if self._embedded_real:
            j = len(self._gen_approx)
            from sage.rings.real_mpfi import RealIntervalField
            gen = self._embedding.gen_image()
            while j <= i:
                self._gen_approx.append(RealIntervalField(53 << j)(gen))
                j += 1
            return self._gen_approx[i]
        else:
            raise ValueError("No embedding set. You need to specify a real embedding.")

    def _matrix_charpoly(self, M, var):
        r"""
        Use PARI to compute the characteristic polynomial of ``self`` as a
        polynomial over the base ring.

        EXAMPLES::

            sage: x = QQ['x'].gen()
            sage: K.<a> = NumberField(x^2 - 2)
            sage: m = matrix(K, [[a-1, 2], [a, a+1]])
            sage: m.charpoly('Z')   # indirect doctest
            Z^2 - 2*a*Z - 2*a + 1
            sage: m.charpoly('a')(m) == 0   # indirect doctest
            True
            sage: m = matrix(K, [[0, a, 0], [-a, 0, 0], [0, 0, 0]])
            sage: m.charpoly('Z')   # indirect doctest
            Z^3 + 2*Z

         ::

            sage: L.<b> = K.extension(x^3 - a)
            sage: m = matrix(L, [[b+a, 1], [a, b^2-2]])
            sage: m.charpoly('Z')   # indirect doctest
            Z^2 + (-b^2 - b - a + 2)*Z + a*b^2 - 2*b - 2*a
            sage: m.charpoly('a')   # indirect doctest
            a^2 + (-b^2 - b - a + 2)*a + a*b^2 - 2*b - 2*a
            sage: m.charpoly('a')(m) == 0
            True

        ::

            sage: M.<c> = L.extension(x^2 - a*x + b)
            sage: m = matrix(M, [[a+b+c, 0, b], [0, c, 1], [a-1, b^2+1, 2]])
            sage: f = m.charpoly('Z'); f    # indirect doctest
            Z^3 + (-2*c - b - a - 2)*Z^2 + ((b + 2*a + 4)*c - b^2 + (-a + 2)*b + 2*a - 1)*Z + (b^2 + (a - 3)*b - 4*a + 1)*c + a*b^2 + 3*b + 2*a
            sage: f(m) == 0
            True
            sage: f.base_ring() is M
            True
        """
        paripoly = M.__pari__().charpoly()
        return self[var](paripoly)
