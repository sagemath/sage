r"""
Local components of modular forms

If `f` is a (new, cuspidal, normalised) modular eigenform, then one can
associate to `f` an *automorphic representation* `\pi_f` of the group
`\operatorname{GL}_2(\mathbf{A})` (where `\mathbf{A}` is the adele ring of
`\QQ`). This object factors as a restricted tensor product of components
`\pi_{f, v}` for each place of `\QQ`. These are infinite-dimensional
representations, but they are specified by a finite amount of data, and this
module provides functions which determine a description of the local factor
`\pi_{f, p}` at a finite prime `p`.

The functions in this module are based on the algorithms described in:

.. [LW11] David Loeffler and Jared Weinstein, *On the computation of local components of a newform*,
   Mathematics of Computation (to appear), 2011. `Online version
   <http://dx.doi.org/10.1090/S0025-5718-2011-02530-5>`_.

AUTHORS:

- David Loeffler
- Jared Weinstein
"""
from __future__ import absolute_import

import operator
from sage.structure.sage_object     import SageObject
from sage.rings.all                 import QQ, ZZ, Zmod, QQbar, PolynomialRing, polygen
from sage.modular.modform.element   import Newform
from sage.modular.dirichlet         import DirichletGroup
from sage.misc.cachefunc            import cached_method
from sage.misc.abstract_method      import abstract_method
from sage.structure.sequence        import Sequence

from .type_space                     import TypeSpace
from .smoothchar                     import SmoothCharacterGroupQp, SmoothCharacterGroupUnramifiedQuadratic, SmoothCharacterGroupRamifiedQuadratic

def LocalComponent(f, p, twist_factor=None):
    r"""
    Calculate the local component at the prime `p` of the automorphic
    representation attached to the newform `f`.

    INPUT:

    - ``f`` (:class:`~sage.modular.modform.element.Newform`) a newform of weight `k \ge 2`
    - ``p`` (integer) a prime
    - ``twist_factor`` (integer) an integer congruent to `k` modulo 2 (default: `k - 2`)

    .. note::

        The argument ``twist_factor`` determines the choice of normalisation: if it is
        set to `j \in \ZZ`, then the central character of `\pi_{f, \ell}` maps `\ell`
        to `\ell^j \varepsilon(\ell)` for almost all `\ell`, where `\varepsilon` is the
        Nebentypus character of `f`.

        In the analytic theory it is conventional to take `j = 0` (the "Langlands
        normalisation"), so the representation `\pi_f` is unitary; however, this is
        inconvenient for `k` odd, since in this case one needs to choose a square root of `p`
        and thus the map `f \to \pi_{f}` is not Galois-equivariant. Hence we use, by default, the
        "Hecke normalisation" given by `j = k - 2`. This is also the most natural normalisation
        from the perspective of modular symbols.

        We also adopt a slightly unusual definition of the principal series: we
        define `\pi(\chi_1, \chi_2)` to be the induction from the Borel subgroup of
        the character of the maximal torus `\begin{pmatrix} x & \\ & y
        \end{pmatrix} \mapsto \chi_1(a) \chi_2(b) |b|`, so its central character is
        `z \mapsto \chi_1(z) \chi_2(z) |z|`. Thus `\chi_1 \chi_2` is the
        restriction to `\QQ_p^\times` of the unique character of the id\'ele class
        group mapping `\ell` to `\ell^{k-1} \varepsilon(\ell)` for almost all `\ell`.
        This has the property that the *set* `\{\chi_1, \chi_2\}` also depends
        Galois-equivariantly on `f`.

    EXAMPLE::

        sage: Pi = LocalComponent(Newform('49a'), 7); Pi
        Smooth representation of GL_2(Q_7) with conductor 7^2
        sage: Pi.central_character()
        Character of Q_7*, of level 0, mapping 7 |--> 1
        sage: Pi.species()
        'Supercuspidal'
        sage: Pi.characters()
        [
        Character of unramified extension Q_7(s)* (s^2 + 6*s + 3 = 0), of level 1, mapping s |--> d, 7 |--> 1,
        Character of unramified extension Q_7(s)* (s^2 + 6*s + 3 = 0), of level 1, mapping s |--> -d, 7 |--> 1
        ]
    """
    p = ZZ(p)
    if not p.is_prime():
        raise ValueError( "p must be prime" )
    if not isinstance(f, Newform):
        raise TypeError( "f (=%s of type %s) should be a Newform object" % (f, type(f)) )

    r = f.level().valuation(p)
    if twist_factor is None:
        twist_factor = ZZ(f.weight() - 2)
    else:
        twist_factor = ZZ(twist_factor)
    if r == 0:
        return UnramifiedPrincipalSeries(f, p, twist_factor)
    c = ZZ(f.character().conductor()).valuation(p)
    if f[p] != 0:
        if c == r:
            return PrimitivePrincipalSeries(f, p, twist_factor)
        if c == 0 and r == 1:
            return PrimitiveSpecial(f, p, twist_factor)
    Xf = TypeSpace(f, p)
    if Xf.is_minimal():
        return PrimitiveSupercuspidal(f, p, twist_factor)
    else:
        raise NotImplementedError( "Form %s is not %s-primitive" % (f, p) )

class LocalComponentBase(SageObject):
    r"""
    Base class for local components of newforms. Not to be directly instantiated; use the :func:`~LocalComponent` constructor function.
    """

    def __init__(self, newform, prime, twist_factor):
        r"""
        Standard initialisation function.

        EXAMPLE::

            sage: LocalComponent(Newform('49a'), 7) # indirect doctest
            Smooth representation of GL_2(Q_7) with conductor 7^2
        """
        self._p = prime
        self._f = newform
        self._twist_factor = twist_factor

    @abstract_method
    def species(self):
        r"""
        The species of this local component, which is either 'Principal
        Series', 'Special' or 'Supercuspidal'.

        EXAMPLE::

            sage: from sage.modular.local_comp.local_comp import LocalComponentBase
            sage: LocalComponentBase(Newform('50a'), 3, 0).species()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method species at ...>
        """
        pass

    @abstract_method
    def check_tempered(self):
        r"""
        Check that this representation is quasi-tempered, i.e. `\pi \otimes
        |\det|^{j/2}` is tempered. It is well known that local components of
        modular forms are *always* tempered, so this serves as a useful check
        on our computations.

        EXAMPLE::

            sage: from sage.modular.local_comp.local_comp import LocalComponentBase
            sage: LocalComponentBase(Newform('50a'), 3, 0).check_tempered()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method check_tempered at ...>
        """
        pass

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: LocalComponent(Newform('50a'), 5)._repr_()
            'Smooth representation of GL_2(Q_5) with conductor 5^2'
        """
        return "Smooth representation of GL_2(Q_%s) with conductor %s^%s" % (self.prime(), self.prime(), self.conductor())

    def newform(self):
        r"""
        The newform of which this is a local component.

        EXAMPLE::

            sage: LocalComponent(Newform('50a'), 5).newform()
            q - q^2 + q^3 + q^4 + O(q^6)
        """
        return self._f

    def prime(self):
        r"""
        The prime at which this is a local component.

        EXAMPLE::

            sage: LocalComponent(Newform('50a'), 5).prime()
            5
        """
        return self._p

    def conductor(self):
        r"""
        The smallest `r` such that this representation has a nonzero vector fixed by the subgroup
        `\begin{pmatrix} * & * \\ 0 & 1\end{pmatrix} \pmod{p^r}`. This is equal to the power of `p` dividing the level of the corresponding newform.

        EXAMPLE::

            sage: LocalComponent(Newform('50a'), 5).conductor()
            2
        """
        return self.newform().level().valuation(self.prime())

    def coefficient_field(self):
        r"""
        The field `K` over which this representation is defined. This is the field generated by the Hecke eigenvalues of the corresponding newform (over whatever base ring the newform is created).

        EXAMPLE::

            sage: LocalComponent(Newforms(50)[0], 3).coefficient_field()
            Rational Field
            sage: LocalComponent(Newforms(Gamma1(10), 3, base_ring=QQbar)[0], 5).coefficient_field()
            Algebraic Field
            sage: LocalComponent(Newforms(DirichletGroup(5).0, 7,names='c')[0], 5).coefficient_field()
            Number Field in c0 with defining polynomial x^2 + (5*zeta4 + 5)*x - 88*zeta4 over its base field
        """
        return self.newform().hecke_eigenvalue_field()

    def twist_factor(self):
        r"""
        The unique `j` such that `\begin{pmatrix} p & 0 \\ 0 & p\end{pmatrix}`
        acts as multiplication by `p^j` times a root of unity.

        There are various conventions for this; see the documentation of the
        :func:`~LocalComponent` constructor function for more information.

        The twist factor should have the same parity as the weight of the form,
        since otherwise the map sending `f` to its local component won't be
        Galois equivariant.

        EXAMPLE::

            sage: LocalComponent(Newforms(50)[0], 3).twist_factor()
            0
            sage: LocalComponent(Newforms(50)[0], 3, twist_factor=173).twist_factor()
            173
        """
        return self._twist_factor

    def central_character(self):
        r"""
        Return the central character of this representation. This is the
        restriction to `\QQ_p^\times` of the unique smooth character `\omega`
        of `\mathbf{A}^\times / \QQ^\times` such that `\omega(\varpi_\ell) =
        \ell^j \varepsilon(\ell)` for all primes `\ell \nmid Np`, where
        `\varpi_\ell` is a uniformiser at `\ell`, `\varepsilon` is the
        Nebentypus character of the newform `f`, and `j` is the twist factor
        (see the documentation for :func:`~LocalComponent`).

        EXAMPLES::

            sage: LocalComponent(Newform('27a'), 3).central_character()
            Character of Q_3*, of level 0, mapping 3 |--> 1

            sage: LocalComponent(Newforms(Gamma1(5), 5, names='c')[0], 5).central_character()
            Character of Q_5*, of level 1, mapping 2 |--> c0 + 1, 5 |--> 125

            sage: LocalComponent(Newforms(DirichletGroup(24)([1, -1,-1]), 3, names='a')[0], 2).central_character()
            Character of Q_2*, of level 3, mapping 7 |--> 1, 5 |--> -1, 2 |--> -2
        """
        from sage.arith.all import crt
        chi = self.newform().character()
        f = self.prime() ** self.conductor()
        N = self.newform().level() // f
        G = DirichletGroup(f, self.coefficient_field())
        chip = G([chi(crt(ZZ(x), 1, f, N)) for x in G.unit_gens()]).primitive_character()
        a = crt(1, self.prime(), f, N)

        if chip.conductor() == 1:
            return SmoothCharacterGroupQp(self.prime(), self.coefficient_field()).character(0, [chi(a) * self.prime()**self.twist_factor()])
        else:
            return SmoothCharacterGroupQp(self.prime(), self.coefficient_field()).character(chip.conductor().valuation(self.prime()), list((~chip).values_on_gens()) + [chi(a) * self.prime()**self.twist_factor()])

    def __eq__(self, other):
        r"""
        Comparison function.

        EXAMPLE::

            sage: Pi = LocalComponent(Newform("50a"), 5)
            sage: Pi == LocalComponent(Newform("50a"), 3)
            False
            sage: Pi == LocalComponent(Newform("50b"), 5)
            False
            sage: Pi == QQ
            False
            sage: Pi == None
            False
            sage: Pi == loads(dumps(Pi))
            True
        """
        return (isinstance(other, LocalComponentBase)
                and self.prime() == other.prime()
                and self.newform() == other.newform()
                and self.twist_factor() == other.twist_factor())

    def __ne__(self, other):
        """
        Return True if ``self != other``.

        EXAMPLE::

            sage: Pi = LocalComponent(Newform("50a"), 5)
            sage: Pi != LocalComponent(Newform("50a"), 3)
            True
            sage: Pi != LocalComponent(Newform("50b"), 5)
            True
            sage: Pi != QQ
            True
            sage: Pi != None
            True
            sage: Pi != loads(dumps(Pi))
            False
        """
        return not (self == other)


class PrincipalSeries(LocalComponentBase):
    r"""
    A principal series representation. This is an abstract base class, not to
    be instantiated directly; see the subclasses
    :class:`~UnramifiedPrincipalSeries` and :class:`~PrimitivePrincipalSeries`.
    """

    def species(self):
        r"""
        The species of this local component, which is either 'Principal
        Series', 'Special' or 'Supercuspidal'.

        EXAMPLE::

            sage: LocalComponent(Newform('50a'), 3).species()
            'Principal Series'
        """
        return "Principal Series"

    def check_tempered(self):
        r"""
        Check that this representation is tempered (after twisting by
        `|\det|^{j/2}`), i.e. that `|\chi_1(p)| = |\chi_2(p)| = p^{(j + 1)/2}`.
        This follows from the Ramanujan--Petersson conjecture, as proved by
        Deligne.

        EXAMPLE::

            sage: LocalComponent(Newform('49a'), 3).check_tempered()
        """
        c1, c2 = self.characters()
        K = c1.base_ring()
        p = self.prime()
        w = QQbar(p)**((1 + self.twist_factor()) / 2)
        for sigma in K.embeddings(QQbar):
            assert sigma(c1(p)).abs() == sigma(c2(p)).abs() == w

    @abstract_method
    def characters(self):
        r"""
        Return the two characters `(\chi_1, \chi_2)` such this representation
        `\pi_{f, p}` is equal to the principal series `\pi(\chi_1, \chi_2)`.

        EXAMPLE::

            sage: from sage.modular.local_comp.local_comp import PrincipalSeries
            sage: PrincipalSeries(Newform('50a'), 3, 0).characters()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method characters at ...>
        """
        pass

class UnramifiedPrincipalSeries(PrincipalSeries):
    r"""
    An unramified principal series representation of `{\rm GL}_2(\QQ_p)`
    (corresponding to a form whose level is not divisible by `p`).

    EXAMPLE::

        sage: Pi = LocalComponent(Newform('50a'), 3)
        sage: Pi.conductor()
        0
        sage: type(Pi)
        <class 'sage.modular.local_comp.local_comp.UnramifiedPrincipalSeries'>
        sage: TestSuite(Pi).run()
    """

    def satake_polynomial(self):
        r"""
        Return the Satake polynomial of this representation, i.e.~the polynomial whose roots are `\chi_1(p), \chi_2(p)`
        where this representation is `\pi(\chi_1, \chi_2)`. Concretely, this is the polynomial

        .. math::

            X^2 - p^{(j - k + 2)/2} a_p(f) X + p^{j + 1} \varepsilon(p)`.

        An error will be raised if `j \ne k \bmod 2`.

        EXAMPLES::

            sage: LocalComponent(Newform('11a'), 17).satake_polynomial()
            X^2 + 2*X + 17
            sage: LocalComponent(Newform('11a'), 17, twist_factor = -2).satake_polynomial()
            X^2 + 2/17*X + 1/17
        """
        p = self.prime()
        return PolynomialRing(self.coefficient_field(), 'X')([
                self.central_character()(p)*p,
                -self.newform()[p] * p**((self.twist_factor() - self.newform().weight() + 2)/2),
                1
              ])

    def characters(self):
        r"""
        Return the two characters `(\chi_1, \chi_2)` such this representation
        `\pi_{f, p}` is equal to the principal series `\pi(\chi_1, \chi_2)`.
        These are the unramified characters mapping `p` to the roots of the Satake polynomial,
        so in most cases (but not always) they will be defined over an
        extension of the coefficient field of self.

        EXAMPLES::

            sage: LocalComponent(Newform('11a'), 17).characters()
            [
            Character of Q_17*, of level 0, mapping 17 |--> d,
            Character of Q_17*, of level 0, mapping 17 |--> -d - 2
            ]
            sage: LocalComponent(Newforms(Gamma1(5), 6, names='a')[1], 3).characters()
            [
            Character of Q_3*, of level 0, mapping 3 |--> -3/2*a1 + 12,
            Character of Q_3*, of level 0, mapping 3 |--> -3/2*a1 - 12
            ]
        """
        f = self.satake_polynomial()
        if not f.is_irreducible():
            # This can happen; see the second example above
            d = f.roots()[0][0]
        else:
            d = self.coefficient_field().extension(f, 'd').gen()
        G = SmoothCharacterGroupQp(self.prime(), d.parent())
        return Sequence([G.character(0, [d]), G.character(0, [self.newform()[self.prime()] - d])], cr=True, universe=G)

class PrimitivePrincipalSeries(PrincipalSeries):
    r"""
    A ramified principal series of the form `\pi(\chi_1, \chi_2)`
    where `\chi_1` is unramified but `\chi_2` is not.

    EXAMPLE::

        sage: Pi = LocalComponent(Newforms(Gamma1(13), 2, names='a')[0], 13)
        sage: type(Pi)
        <class 'sage.modular.local_comp.local_comp.PrimitivePrincipalSeries'>
        sage: TestSuite(Pi).run()
    """

    def characters(self):
        r"""
        Return the two characters `(\chi_1, \chi_2)` such that the local component `\pi_{f, p}` is the induction of the character `\chi_1 \times \chi_2` of the Borel subgroup.

        EXAMPLE::

            sage: LocalComponent(Newforms(Gamma1(13), 2, names='a')[0], 13).characters()
            [
            Character of Q_13*, of level 0, mapping 13 |--> 3*a0 + 2,
            Character of Q_13*, of level 1, mapping 2 |--> a0 + 2, 13 |--> -3*a0 - 7
            ]
        """
        G = SmoothCharacterGroupQp(self.prime(), self.coefficient_field())
        chi1 = G.character(0, [self.newform()[self.prime()]])
        chi2 = G.character(0, [self.prime()]) * self.central_character() / chi1
        return Sequence([chi1, chi2], cr=True, universe=G)

class PrimitiveSpecial(LocalComponentBase):
    r"""
    A primitive special representation: that is, the Steinberg representation
    twisted by an unramified character. All such representations have conductor
    1.

    EXAMPLES::

        sage: Pi = LocalComponent(Newform('37a'), 37)
        sage: Pi.species()
        'Special'
        sage: Pi.conductor()
        1
        sage: type(Pi)
        <class 'sage.modular.local_comp.local_comp.PrimitiveSpecial'>
        sage: TestSuite(Pi).run()
    """

    def species(self):
        r"""
        The species of this local component, which is either 'Principal
        Series', 'Special' or 'Supercuspidal'.

        EXAMPLE::

            sage: LocalComponent(Newform('37a'), 37).species()
            'Special'
        """
        return "Special"

    def characters(self):
        r"""
        Return the defining characters of this representation. In this case, it
        will return the unique unramified character `\chi` of `\QQ_p^\times`
        such that this representation is equal to `\mathrm{St} \otimes \chi`,
        where `\mathrm{St}` is the Steinberg representation (defined as the
        quotient of the parabolic induction of the trivial character by its
        trivial subrepresentation).

        EXAMPLES:

        Our first example is the newform corresponding to an elliptic curve of
        conductor `37`. This is the nontrivial quadratic twist of Steinberg,
        corresponding to the fact that the elliptic curve has non-split
        multiplicative reduction at 37::

            sage: LocalComponent(Newform('37a'), 37).characters()
            [Character of Q_37*, of level 0, mapping 37 |--> -1]

        We try an example in odd weight, where the central character isn't
        trivial::

            sage: Pi = LocalComponent(Newforms(DirichletGroup(21)([-1, 1]), 3, names='j')[0], 7); Pi.characters()
            [Character of Q_7*, of level 0, mapping 7 |--> -1/2*j0^2 - 7/2]
            sage: Pi.characters()[0] ^2 == Pi.central_character()
            True

        An example using a non-standard twist factor::

            sage: Pi = LocalComponent(Newforms(DirichletGroup(21)([-1, 1]), 3, names='j')[0], 7, twist_factor=3); Pi.characters()
            [Character of Q_7*, of level 0, mapping 7 |--> -7/2*j0^2 - 49/2]
            sage: Pi.characters()[0]^2 == Pi.central_character()
            True
        """

        return [SmoothCharacterGroupQp(self.prime(), self.coefficient_field()).character(0, [self.newform()[self.prime()] * self.prime() ** ((self.twist_factor() - self.newform().weight() + 2)/2)])]

    def check_tempered(self):
        r"""
        Check that this representation is tempered (after twisting by
        `|\det|^{j/2}` where `j` is the twist factor). Since local components
        of modular forms are always tempered, this is a useful check on our
        calculations.

        EXAMPLE::

            sage: Pi = LocalComponent(Newforms(DirichletGroup(21)([-1, 1]), 3, names='j')[0], 7)
            sage: Pi.check_tempered()
        """
        c1 = self.characters()[0]
        K = c1.base_ring()
        p = self.prime()
        w = QQbar(p)**(self.twist_factor() / ZZ(2))
        for sigma in K.embeddings(QQbar):
            assert sigma(c1(p)).abs() == w

class PrimitiveSupercuspidal(LocalComponentBase):
    r"""
    A primitive supercuspidal representation. Except for some excpetional cases
    when `p = 2` which we do not implement here, such representations are
    parametrized by smooth characters of tamely ramified quadratic extensions
    of `\QQ_p`.

    EXAMPLES::

        sage: f = Newform("50a")
        sage: Pi = LocalComponent(f, 5)
        sage: type(Pi)
        <class 'sage.modular.local_comp.local_comp.PrimitiveSupercuspidal'>
        sage: Pi.species()
        'Supercuspidal'
        sage: TestSuite(Pi).run()
    """

    def species(self):
        r"""
        The species of this local component, which is either 'Principal
        Series', 'Special' or 'Supercuspidal'.

        EXAMPLE::

            sage: LocalComponent(Newform('49a'), 7).species()
            'Supercuspidal'
        """
        return "Supercuspidal"

    @cached_method
    def type_space(self):
        r"""
        Return a :class:`~sage.modular.local_comp.type_space.TypeSpace` object
        describing the (homological) type space of this newform, which we know
        is dual to the type space of the local component.

        EXAMPLE::

            sage: LocalComponent(Newform('49a'), 7).type_space()
            6-dimensional type space at prime 7 of form q + q^2 - q^4 + O(q^6)
        """
        return TypeSpace(self.newform(), self.prime())

    def characters(self):
        r"""
        Return the two conjugate characters of `K^\times`, where `K` is some
        quadratic extension of `\QQ_p`, defining this representation. This is
        fully implemented only in the case where the power of `p` dividing the
        level of the form is even, in which case `K` is the unique unramified
        quadratic extension of `\QQ_p`.

        EXAMPLES:

        The first example from _[LW11]::

            sage: f = Newform('50a')
            sage: Pi = LocalComponent(f, 5)
            sage: chars = Pi.characters(); chars
            [
            Character of unramified extension Q_5(s)* (s^2 + 4*s + 2 = 0), of level 1, mapping s |--> d, 5 |--> 1,
            Character of unramified extension Q_5(s)* (s^2 + 4*s + 2 = 0), of level 1, mapping s |--> -d - 1, 5 |--> 1
            ]
            sage: chars[0].base_ring()
            Number Field in d with defining polynomial x^2 + x + 1

        These characters are interchanged by the Frobenius automorphism of `\mathbb{F}_{25}`::

            sage: chars[0] == chars[1]**5
            True

        A more complicated example (higher weight and nontrivial central character)::

            sage: f = Newforms(GammaH(25, [6]), 3, names='j')[0]; f
            q + j0*q^2 + 1/3*j0^3*q^3 - 1/3*j0^2*q^4 + O(q^6)
            sage: Pi = LocalComponent(f, 5)
            sage: Pi.characters()
            [
            Character of unramified extension Q_5(s)* (s^2 + 4*s + 2 = 0), of level 1, mapping s |--> d, 5 |--> 5,
            Character of unramified extension Q_5(s)* (s^2 + 4*s + 2 = 0), of level 1, mapping s |--> -d - 1/3*j0^3, 5 |--> 5
            ]
            sage: Pi.characters()[0].base_ring()
            Number Field in d with defining polynomial x^2 + 1/3*j0^3*x - 1/3*j0^2 over its base field

        .. warning::

            The above output isn't actually the same as in Example 2 of
            _[LW11], due to an error in the published paper (correction
            pending) -- the published paper has the inverses of the above
            characters.

        A higher level example::

            sage: f = Newform('81a', names='j'); f
            q + j0*q^2 + q^4 - j0*q^5 + O(q^6)
            sage: LocalComponent(f, 3).characters()  # long time (12s on sage.math, 2012)
            [
            Character of unramified extension Q_3(s)* (s^2 + 2*s + 2 = 0), of level 2, mapping -2*s |--> -2*d - j0, 4 |--> 1, 3*s + 1 |--> -j0*d - 2, 3 |--> 1,
            Character of unramified extension Q_3(s)* (s^2 + 2*s + 2 = 0), of level 2, mapping -2*s |--> 2*d + j0, 4 |--> 1, 3*s + 1 |--> j0*d + 1, 3 |--> 1
            ]

        In the ramified case, it's not fully implemented, and just returns a
        string indicating which ramified extension is being considered::

            sage: Pi = LocalComponent(Newform('27a'), 3)
            sage: Pi.characters()
            'Character of Q_3(sqrt(-3))'
            sage: Pi = LocalComponent(Newform('54a'), 3)
            sage: Pi.characters()
            'Character of Q_3(sqrt(3))'
        """
        T = self.type_space()
        if self.conductor() % 2 == 0:

            G = SmoothCharacterGroupUnramifiedQuadratic(self.prime(), self.coefficient_field())
            n = self.conductor() // 2
            g = G.quotient_gen(n)
            m = g.matrix().change_ring(ZZ).list()
            tr = (~T.rho(m)).trace()

            # The inverse is needed here because T is the *homological* type space,
            # which is dual to the cohomological one that defines the local component.

            X = polygen(self.coefficient_field())
            theta_poly = X**2 - (-1)**n*tr*X + self.central_character()(g.norm())
            if theta_poly.is_irreducible():
                F = self.coefficient_field().extension(theta_poly, "d")
                G = G.base_extend(F)
            chi1, chi2 = [G.extend_character(n, self.central_character(), x[0]) for x in theta_poly.roots(G.base_ring())]

            # Consistency checks
            assert chi1.restrict_to_Qp() == chi2.restrict_to_Qp() == self.central_character()
            assert chi1*chi2 == chi1.parent().compose_with_norm(self.central_character())

            return Sequence([chi1, chi2], check=False, cr=True)

        else:
            # The ramified case.

            p = self.prime()

            if p == 2:
                # The ramified 2-adic representations aren't classified by admissible pairs. Die.
                raise NotImplementedError( "Computation with ramified 2-adic representations not implemented" )

            if p % 4 == 3:
                a = ZZ(-1)
            else:
                a = ZZ(Zmod(self.prime()).quadratic_nonresidue())

            tr1 = (~T.rho([0,1,a*p, 0])).trace()
            tr2 = (~T.rho([0,1,p,0])).trace()

            if tr1 == tr2 == 0:
                # This *can* happen. E.g. if the central character satisfies
                # chi(-1) = -1, then we have theta(pi) + theta(-pi) = theta(pi)
                # * (1 + -1) = 0. In this case, one can presumably identify
                # the character and the extension by some more subtle argument
                # but I don't know of a good way to automate the process.
                raise NotImplementedError( "Can't identify ramified quadratic extension -- both traces zero" )
            elif tr1 == 0:
                return "Character of Q_%s(sqrt(%s))" % (p, p)

            elif tr2 == 0:
                return "Character of Q_%s(sqrt(%s))" % (p, a*p)

            else:
                # At least one of the traces is *always* 0, since the type
                # space has to be isomorphic to its twist by the (ramified
                # quadratic) character corresponding to the quadratic
                # extension.
                raise RuntimeError( "Can't get here!" )

    def check_tempered(self):
        r"""
        Check that this representation is tempered (after twisting by
        `|\det|^{j/2}` where `j` is the twist factor). Since local components
        of modular forms are always tempered, this is a useful check on our
        calculations.

        Since the computation of the characters attached to this representation
        is not implemented in the odd-conductor case, a NotImplementedError
        will be raised for such representations.

        EXAMPLE::

            sage: LocalComponent(Newform("50a"), 5).check_tempered()
            sage: LocalComponent(Newform("27a"), 3).check_tempered() # not tested
        """
        if self.conductor() % 2:
            raise NotImplementedError
        c1, c2 = self.characters()
        K = c1.base_ring()
        p = self.prime()
        w = QQbar(p)**(self.twist_factor() / ZZ(2))
        for sigma in K.embeddings(QQbar):
            assert c1(p).abs() == c2(p).abs() == w
