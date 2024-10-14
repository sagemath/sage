r"""
Localization

Localization is an important ring construction tool. Whenever you have
to extend a given integral domain such that it contains the inverses
of a finite set of elements but should allow non injective homomorphic
images this construction will be needed. See the example on
Ariki-Koike algebras below for such an application.

EXAMPLES::

    sage: # needs sage.modules
    sage: LZ = Localization(ZZ, (5,11))
    sage: m = matrix(LZ, [[5, 7], [0,11]])
    sage: m.parent()
    Full MatrixSpace of 2 by 2 dense matrices over Integer Ring localized at (5, 11)
    sage: ~m      # parent of inverse is different: see documentation of m.__invert__
    [  1/5 -7/55]
    [    0  1/11]
    sage: _.parent()
    Full MatrixSpace of 2 by 2 dense matrices over Rational Field
    sage: mi = matrix(LZ, ~m)
    sage: mi.parent()
    Full MatrixSpace of 2 by 2 dense matrices over Integer Ring localized at (5, 11)
    sage: mi == ~m
    True

The next example defines the most general ring containing the coefficients of the irreducible
representations of the Ariki-Koike algebra corresponding to the three colored permutations on
three elements::

    sage: R.<u0, u1, u2, q> = ZZ[]
    sage: u = [u0, u1, u2]
    sage: S = Set(u)
    sage: I = S.cartesian_product(S)
    sage: add_units = u + [q, q + 1] + [ui - uj for ui, uj in I if ui != uj]
    sage: add_units += [q*ui - uj for ui, uj in I if ui != uj]
    sage: L = R.localization(tuple(add_units)); L                                       # needs sage.libs.pari
    Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring localized at
     (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1, u2*q - u1, u2*q - u0,
      u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)

Define the representation matrices (of one of the three dimensional irreducible representations)::

    sage: # needs sage.libs.pari sage.modules
    sage: m1 = matrix(L, [[u1, 0, 0], [0, u0, 0], [0, 0, u0]])
    sage: m2 = matrix(L, [[(u0*q - u0)/(u0 - u1), (u0*q - u1)/(u0 - u1), 0],
    ....:                 [(-u1*q + u0)/(u0 - u1), (-u1*q + u1)/(u0 - u1), 0],
    ....:                 [0, 0, -1]])
    sage: m3 = matrix(L, [[-1, 0, 0],
    ....:                 [0, u0*(1 - q)/(u1*q - u0), q*(u1 - u0)/(u1*q - u0)],
    ....:                 [0, (u1*q^2 - u0)/(u1*q - u0), (u1*q^ 2 - u1*q)/(u1*q - u0)]])
    sage: m1.base_ring() == L
    True

Check relations of the Ariki-Koike algebra::

    sage: # needs sage.libs.pari sage.modules
    sage: m1*m2*m1*m2 == m2*m1*m2*m1
    True
    sage: m2*m3*m2 == m3*m2*m3
    True
    sage: m1*m3 == m3*m1
    True
    sage: m1**3 - (u0+u1+u2)*m1**2 + (u0*u1+u0*u2+u1*u2)*m1 - u0*u1*u2 == 0
    True
    sage: m2**2 - (q-1)*m2 - q == 0
    True
    sage: m3**2 - (q-1)*m3 - q == 0
    True
    sage: ~m1 in m1.parent()
    True
    sage: ~m2 in m2.parent()
    True
    sage: ~m3 in m3.parent()
    True

Obtain specializations in positive characteristic::

    sage: # needs sage.libs.pari sage.modules
    sage: Fp = GF(17)
    sage: f = L.hom((3,5,7,11), codomain=Fp); f
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring localized at
            (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1, u2*q - u1, u2*q - u0,
             u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Finite Field of size 17
      Defn: u0 |--> 3
            u1 |--> 5
            u2 |--> 7
            q |--> 11
    sage: mFp1 = matrix({k: f(v) for k, v in m1.dict().items()}); mFp1
    [5 0 0]
    [0 3 0]
    [0 0 3]
    sage: mFp1.base_ring()
    Finite Field of size 17
    sage: mFp2 = matrix({k: f(v) for k, v in m2.dict().items()}); mFp2
    [ 2  3  0]
    [ 9  8  0]
    [ 0  0 16]
    sage: mFp3 = matrix({k: f(v) for k, v in m3.dict().items()}); mFp3
    [16  0  0]
    [ 0  4  5]
    [ 0  7  6]

Obtain specializations in characteristic 0::

    sage: # needs sage.libs.pari
    sage: fQ = L.hom((3,5,7,11), codomain=QQ); fQ
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring
            localized at (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1,
            u2*q - u1, u2*q - u0, u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Rational Field
      Defn: u0 |--> 3
            u1 |--> 5
            u2 |--> 7
            q |--> 11

    sage: # needs sage.libs.pari sage.modules sage.rings.finite_rings
    sage: mQ1 = matrix({k: fQ(v) for k, v in m1.dict().items()}); mQ1
    [5 0 0]
    [0 3 0]
    [0 0 3]
    sage: mQ1.base_ring()
    Rational Field
    sage: mQ2 = matrix({k: fQ(v) for k, v in m2.dict().items()}); mQ2
    [-15 -14   0]
    [ 26  25   0]
    [  0   0  -1]
    sage: mQ3 = matrix({k: fQ(v) for k, v in m3.dict().items()}); mQ3
    [    -1      0      0]
    [     0 -15/26  11/26]
    [     0 301/26 275/26]

    sage: # needs sage.libs.pari sage.libs.singular
    sage: S.<x, y, z, t> = QQ[]
    sage: T = S.quo(x + y + z)
    sage: F = T.fraction_field()
    sage: fF = L.hom((x, y, z, t), codomain=F); fF
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring
            localized at (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1,
            u2*q - u1, u2*q - u0, u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Fraction Field of Quotient of Multivariate Polynomial Ring in x, y, z, t
            over Rational Field by the ideal (x + y + z)
      Defn: u0 |--> -ybar - zbar
            u1 |--> ybar
            u2 |--> zbar
            q |--> tbar
    sage: mF1 = matrix({k: fF(v) for k, v in m1.dict().items()}); mF1                   # needs sage.modules
    [        ybar            0            0]
    [           0 -ybar - zbar            0]
    [           0            0 -ybar - zbar]
    sage: mF1.base_ring() == F                                                          # needs sage.modules
    True

TESTS::

    sage: TestSuite(L).run()                                                            # needs sage.libs.pari sage.libs.singular sage.modules

AUTHORS:

- Sebastian Oehms 2019-12-09: initial version.
- Sebastian Oehms 2022-03-05: fix some corner cases and add :meth:`factor` (:issue:`33463`)
"""


# ****************************************************************************
#       Copyright (C) 2019 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.integral_domains import IntegralDomains
from sage.rings.ring import IntegralDomain
from sage.structure.element import IntegralDomainElement


def normalize_extra_units(base_ring, add_units, warning=True):
    """
    Function to normalize input data.

    The given list will be replaced by a list of the involved prime factors
    (if possible).

    INPUT:

    - ``base_ring`` -- an instance of :class:`IntegralDomain`
    - ``add_units`` -- list of elements from base ring
    - ``warning`` -- boolean (default: ``True``); to suppress a warning which
      is thrown if no normalization was possible

    OUTPUT: list of all prime factors of the elements of the given list

    EXAMPLES::

        sage: from sage.rings.localization import normalize_extra_units
        sage: normalize_extra_units(ZZ, [3, -15, 45, 9, 2, 50])
        [2, 3, 5]
        sage: P.<x,y,z> = ZZ[]
        sage: normalize_extra_units(P,                                                  # needs sage.libs.pari
        ....:                       [3*x, z*y**2, 2*z, 18*(x*y*z)**2, x*z, 6*x*z, 5])
        [2, 3, 5, z, y, x]
        sage: P.<x,y,z> = QQ[]
        sage: normalize_extra_units(P,                                                  # needs sage.libs.pari
        ....:                       [3*x, z*y**2, 2*z, 18*(x*y*z)**2, x*z, 6*x*z, 5])
        [z, y, x]

        sage: # needs sage.libs.singular
        sage: R.<x, y> = ZZ[]
        sage: Q.<a, b> = R.quo(x**2 - 5)
        sage: p = b**2 - 5
        sage: p == (b-a)*(b+a)
        True
        sage: normalize_extra_units(Q, [p])                                             # needs sage.libs.pari
        doctest:...: UserWarning: Localization may not be represented uniquely
        [b^2 - 5]
        sage: normalize_extra_units(Q, [p], warning=False)                              # needs sage.libs.pari
        [b^2 - 5]
    """
    # convert to base ring
    add_units = [base_ring(n) for n in add_units]

    # split down to prime factors if possible
    add_units_result = []
    for n in add_units:
        try:
            if n.is_unit():
                continue
            F = list(n.factor())
            add_units_result += [f[0] for f in F]
        except (NotImplementedError, AttributeError):
            # if :meth:`is_unit` or :meth:`factor` are not available we can't do any more.
            if warning:
                from warnings import warn
                warn('Localization may not be represented uniquely')
            add_units_result = add_units
            break

    return sorted(set(add_units_result))


class LocalizationElement(IntegralDomainElement):
    """
    Element class for localizations of integral domains.

    INPUT:

    - ``parent`` -- instance of :class:`Localization`
    - ``x`` -- instance of :class:`FractionFieldElement` whose parent is the
      fraction field of the parent's base ring

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: from sage.rings.localization import LocalizationElement
        sage: P.<x,y,z> = GF(5)[]
        sage: L = P.localization((x, y*z - x))
        sage: LocalizationElement(L, 4/(y*z-x)**2)
        (-1)/(y^2*z^2 - 2*x*y*z + x^2)
        sage: _.parent()
        Multivariate Polynomial Ring in x, y, z over Finite Field of size 5
         localized at (x, y*z - x)
    """

    def __init__(self, parent, x):
        """
        Python constructor for the element class for localizations of integral domains.

        EXAMPLES::

            sage: from sage.rings.localization import LocalizationElement
            sage: P.<x> = RR[]
            sage: L = Localization(P, x**2 + x + 1)                                     # needs sage.libs.pari
            sage: l = LocalizationElement(L, (x**2+1)/(x**2+x+1))                       # needs sage.libs.pari
            sage: l._value == (x**2+1)/(x**2+x+1)                                       # needs sage.libs.pari
            True
        """
        IntegralDomainElement.__init__(self, parent)
        self._value = x

    def _repr_(self):
        """
        How to print ``self``.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: from sage.rings.localization import LocalizationElement
            sage: P.<x> = CC[]
            sage: L = Localization(P, x**2 + x + 1)
            sage: l = LocalizationElement(L, (x**2+1)/(x**2+x+1))
            sage: l._repr_() == str(l)
            True

            sage: R.<X, Y> = ZZ[]
            sage: L.<x, y> = R.localization(X - Y)                                      # needs sage.libs.singular
            sage: x*y/(x-y)                                                             # needs sage.libs.singular
            x*y/(x - y)
        """
        s = "%s" % self._value
        L = self.parent()
        names = L._names
        if names:
            n = len(names)
            bnames = L.base_ring()._names
            if bnames != names and n == len(bnames):
                # replace separate names (see :issue:`33482`)
                for i in range(n):
                    s = s.replace(bnames[i], names[i])
        return s

    def _add_(left, right):
        """
        Compute addition with another instance of ``self`` (via `+` operator).

        EXAMPLES::

            sage: L = Localization(ZZ, (5,11))      # indirect doctest
            sage: L(1/5) + L(2/11)
            21/55
        """
        return left.parent()._fraction_to_element(left._value + right._value)

    def _sub_(left, right):
        """
        Compute subtraction with another instance of ``self`` (via `-` operator).

        EXAMPLES::

            sage: L = Localization(ZZ, (5,11))
            sage: L(2) - L(1/121)                  # indirect doctest
            241/121
        """
        return left.parent()._fraction_to_element(left._value - right._value)

    def _mul_(left, right):
        """
        Compute multiplication with another instance of ``self`` (via `*` operator).

        EXAMPLES::

            sage: L = Localization(ZZ, (5,11))
            sage: L(1/55) * L(2/25)               # indirect doctest
            2/1375
        """
        return left.parent()._fraction_to_element(left._value * right._value)

    def _div_(left, right):
        """
        Compute division with another instance of ``self`` (via `/` operator).

        EXAMPLES::

            sage: L = Localization(ZZ, (5,11))
            sage: L(1/5) / L(5/11)                # indirect doctest
            11/25
        """
        return left.parent()._fraction_to_element(left._value / right._value)

    def _rmul_(self, c):
        """
        Compute right multiplication with an instance of the base ring of ``self`` (via `*` operator).

        EXAMPLES::

            sage: L = Localization(ZZ, (5,11))
            sage: L(2/11) * 11                    # indirect doctest
            2
        """
        return self.parent()._fraction_to_element(c * self._value)

    def _lmul_(self, c):
        """
        Compute left multiplication with an instance of the base ring of ``self`` (via `*` operator).

        EXAMPLES::

            sage: L = Localization(ZZ, (5,11))
            sage: 7 * L(3/5)                      # indirect doctest
            21/5
        """
        return self.parent()._fraction_to_element(self._value * c)

    def factor(self, proof=None):
        r"""
        Return the factorization of this polynomial.

        INPUT:

        - ``proof`` -- (optional) if given it is passed to the
          corresponding method of the numerator of ``self``

        EXAMPLES::

            sage: P.<X, Y> = QQ['x, y']
            sage: L = P.localization(X - Y)
            sage: x, y = L.gens()                                                       # needs sage.libs.singular
            sage: p = (x^2 - y^2)/(x-y)^2                                               # needs sage.libs.singular
            sage: p.factor()                                                            # needs sage.libs.singular
            (1/(x - y)) * (x + y)
        """
        num = self._value.numerator()
        den = self._value.denominator()
        if proof is not None:
            F = num.factor(proof=proof)
        else:
            F = num.factor()
        P = self.parent()
        fac = [(P(f), e) for (f, e) in F]
        from sage.structure.factorization import Factorization
        return Factorization(fac, unit=~P(den)*F.unit())

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: L = Localization(R, x**2 + 1)                                         # needs sage.libs.pari
            sage: f = L.hom([5], codomain=Localization(ZZ, 26))   # indirect doctest    # needs sage.libs.pari
            sage: f(x/(x**2+1))                                                         # needs sage.libs.pari
            5/26
        """
        return self._value._im_gens_(codomain, im_gens, base_map=base_map)

    def numerator(self):
        """
        Return the numerator of ``self``.

        EXAMPLES::

            sage: L = ZZ.localization((3,5))
            sage: L(7/15).numerator()
            7
        """
        return self._value.numerator()

    def denominator(self):
        """
        Return the denominator of ``self``.

        EXAMPLES::

            sage: L = Localization(ZZ, (3,5))
            sage: L(7/15).denominator()
            15
        """
        return self._value.denominator()

    def is_unit(self):
        """
        Return ``True`` if ``self`` is a unit.

        EXAMPLES::

            sage: # needs sage.libs.pari sage.singular
            sage: P.<x,y,z> = QQ[]
            sage: L = P.localization((x, y*z))
            sage: L(y*z).is_unit()
            True
            sage: L(z).is_unit()
            True
            sage: L(x*y*z).is_unit()
            True
        """
        return self.parent()._cut_off_extra_units_from_base_ring_element(self._value.numerator()).is_unit()

    def inverse_of_unit(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: P.<x,y,z> = ZZ[]
            sage: L = Localization(P, x*y*z)
            sage: L(x*y*z).inverse_of_unit()                                            # needs sage.libs.singular
            1/(x*y*z)
            sage: L(z).inverse_of_unit()                                                # needs sage.libs.singular
            1/z
        """
        parent = self.parent()
        if not self.is_unit():
            raise ArithmeticError("element is not a unit")
        return parent.element_class(parent, ~(parent._fraction_field(self)))

    def _richcmp_(self, other, op):
        """
        EXAMPLES::

           sage: # needs sage.libs.singular
           sage: P.<x,y,z> = GF(7)[]
           sage: L = Localization(P, (x, y, z))
           sage: L(1/x) < L(3/(x*y*z)**3)
           False
           sage: ~L(y*z/x) == L(x/(y*z))
           True
        """
        sval = self._value
        oval = other._value
        return sval._richcmp_(oval, op)

    def __hash__(self):
        """
        Return the hash of the corresponding fraction field element.

        EXAMPLES::

            sage: L = ZZ.localization(5)
            sage: l5 = L(5); l7 = L(7)
            sage: {l5: ~l5, l7: 7}              # indirect doctest
            {5: 1/5, 7: 7}
        """
        return hash(self._value)

    def _rational_(self):
        r"""
        Convert ``self``  to a rational.

        This is only possible if its base ring is the ring of integers.

        OUTPUT: a rational

        TESTS::

            sage: L = ZZ.localization(5)
            sage: cp3 = cyclotomic_polynomial(3).change_ring(L)
            sage: cp3.splitting_field('t')      #   indirect doctest                    # needs sage.libs.pari sage.rings.number_field
            Number Field in t with defining polynomial x^2 + x + 1
        """
        from sage.rings.rational_field import QQ
        if not self._value.parent() == QQ:
            raise ValueError('{} is not a rational'.format(self))
        return self._value

    def _integer_(self, Z=None):
        r"""
        Convert ``self``  to an integer.

        This is only possible if its base ring is the ring of integers and
        the denominator of ``self`` is one.

        OUTPUT: integer

        TESTS::

            sage: L = ZZ.localization(5)
            sage: L(5) in ZZ                  # indirect doctest
            True
        """
        from sage.rings.rational_field import QQ
        if not self._value.parent() == QQ:
            raise ValueError('{} is not a rational'.format(self))
        return self._value._integer_(Z=Z)


class Localization(IntegralDomain, UniqueRepresentation):
    r"""
    The localization generalizes the construction of the field of fractions of
    an integral domain to an arbitrary ring. Given a (not necessarily
    commutative) ring `R` and a subset `S` of `R`, there exists a ring
    `R[S^{-1}]` together with the ring homomorphism `R \longrightarrow R[S^{-1}]`
    that "inverts" `S`; that is, the homomorphism maps elements in `S` to unit
    elements in `R[S^{-1}]` and, moreover, any ring homomorphism from `R` that
    "inverts" `S` uniquely factors through `R[S^{-1}]`.

    The ring `R[S^{-1}]` is called the *localization* of `R` with respect to
    `S`. For example, if `R` is a commutative ring and `f` an element in `R`,
    then the localization consists of elements of the form
    `r/f, r\in R, n \geq 0` (to be precise, `R[f^{-1}] = R[t]/(ft-1)`).

    The above text is taken from `Wikipedia`. The construction here used for
    this class relies on the construction of the field of fraction and is
    therefore restricted to integral domains.

    Accordingly, this class is inherited from :class:`IntegralDomain` and can
    only be used in that context. Furthermore, the base ring should support
    :meth:`sage.structure.element.CommutativeRingElement.divides` and the exact
    division operator ``//`` (:meth:`sage.structure.element.Element.__floordiv__`)
    in order to guarantee a successful application.

    INPUT:

    - ``base_ring`` -- an instance of :class:`Ring` allowing the construction
      of :meth:`fraction_field` (that is an integral domain)
    - ``extra_units`` -- tuple of elements of ``base_ring`` which should be
      turned into units
    - ``names`` -- passed to :class:`IntegralDomain`
    - ``normalize`` -- boolean (default: ``True``); passed to :class:`IntegralDomain`
    - ``category`` -- (default: ``None``) passed to :class:`IntegralDomain`
    - ``warning`` -- boolean (default: ``True``); to suppress a warning which
      is thrown if ``self`` cannot be represented uniquely

    REFERENCES:

    - :wikipedia:`Ring_(mathematics)#Localization`

    EXAMPLES::

        sage: L = Localization(ZZ, (3,5))
        sage: 1/45 in L
        True
        sage: 1/43 in L
        False

        sage: Localization(L, (7,11))
        Integer Ring localized at (3, 5, 7, 11)
        sage: _.is_subring(QQ)
        True

        sage: L(~7)
        Traceback (most recent call last):
        ...
        ValueError: factor 7 of denominator is not a unit

        sage: Localization(Zp(7), (3, 5))                                               # needs sage.rings.padics
        Traceback (most recent call last):
        ...
        ValueError: all given elements are invertible in
        7-adic Ring with capped relative precision 20

        sage: # needs sage.libs.pari
        sage: R.<x> = ZZ[]
        sage: L = R.localization(x**2 + 1)
        sage: s = (x+5)/(x**2+1)
        sage: s in L
        True
        sage: t = (x+5)/(x**2+2)
        sage: t in L
        False
        sage: L(t)
        Traceback (most recent call last):
        ...
        TypeError: fraction must have unit denominator
        sage: L(s) in R
        False
        sage: y = L(x)
        sage: g = L(s)
        sage: g.parent()
        Univariate Polynomial Ring in x over Integer Ring localized at (x^2 + 1,)
        sage: f = (y+5)/(y**2+1); f
        (x + 5)/(x^2 + 1)
        sage: f == g
        True
        sage: (y+5)/(y**2+2)
        Traceback (most recent call last):
        ...
        ValueError: factor x^2 + 2 of denominator is not a unit

        sage: Lau.<u, v> = LaurentPolynomialRing(ZZ)                                    # needs sage.modules
        sage: LauL = Lau.localization(u + 1)                                            # needs sage.modules
        sage: LauL(~u).parent()                                                         # needs sage.modules
        Multivariate Polynomial Ring in u, v over Integer Ring localized at (v, u, u + 1)

    More examples will be shown typing ``sage.rings.localization?``

    TESTS:

    Check that :issue:`33463` is fixed::

        sage: R = ZZ.localization(5)
        sage: R.localization(~5)
        Integer Ring localized at (5,)
    """

    Element = LocalizationElement

    def __init__(self, base_ring, extra_units, names=None, normalize=True, category=None, warning=True):
        """
        Python constructor of Localization.

        TESTS::

            sage: L = Localization(ZZ, (3, 5))
            sage: TestSuite(L).run()

            sage: R.<x> = ZZ[]
            sage: L = R.localization(x**2 + 1)                                          # needs sage.libs.pari
            sage: TestSuite(L).run()
        """
        if type(extra_units) is tuple:
            extra_units = list(extra_units)
        if not type(extra_units) is list:
            extra_units = [extra_units]

        from sage.rings.polynomial.laurent_polynomial_ring_base import LaurentPolynomialRing_generic

        if isinstance(base_ring, LaurentPolynomialRing_generic):
            extra_units += list(base_ring.gens())
            base_ring = base_ring.polynomial_ring()

        if isinstance(base_ring, Localization):
            # don't allow recursive constructions
            extra_units = [u for u in extra_units
                           if ~u not in base_ring._extra_units]  # :issue:`33463`
            extra_units += base_ring._extra_units
            base_ring = base_ring.base_ring()

        extra_units = normalize_extra_units(base_ring, extra_units, warning=warning)

        if not extra_units:
            raise ValueError('all given elements are invertible in %s' % (base_ring))

        if category is None:
            # since by construction the base ring must contain non units self must be infinite
            category = IntegralDomains().Infinite()

        IntegralDomain.__init__(self, base_ring, names=names, normalize=normalize, category=category)
        self._extra_units = tuple(extra_units)
        self._fraction_field = base_ring.fraction_field()
        self._populate_coercion_lists_()

    def _repr_(self):
        """
        How to print ``self``.

        EXAMPLES::

            sage: R.<a> = GF(3)[]
            sage: Localization(R, a**2 - 1)                                             # needs sage.libs.pari
            Univariate Polynomial Ring in a over Finite Field of size 3
             localized at (a + 1, a + 2)
        """
        return "%s localized at %s" % (self.base(), self._extra_units)

    def _element_constructor_(self, x):
        """
        Make sure x is a valid member of self, and return the constructed element.

        EXAMPLES::

            sage: L = Localization(ZZ, (5, 2))
            sage: L._element_constructor_(1/25)
            1/25
            sage: L._element_constructor_(1/20)
            1/20
            sage: L._element_constructor_(1/10)
            1/10
        """
        if isinstance(x, LocalizationElement):
            x = x._value
        else:
            x = self._fraction_field(x)
        return self._fraction_to_element(x)

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: L = Localization(R, x**2 + 1)                                         # needs sage.libs.pari
            sage: L.hom([5])   # indirect doctest                                       # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: images of some localized elements fail to be units

            sage: L.hom([5], codomain=Localization(ZZ, 26))   # indirect doctest        # needs sage.libs.pari
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
                    localized at (x^2 + 1,)
              To:   Integer Ring localized at (2, 13)
              Defn: x |--> 5

        TESTS::

            sage: # needs sage.libs.pari
            sage: phi = R.hom([5])
            sage: L._is_valid_homomorphism_(ZZ, [3], base_map=phi)
            Traceback (most recent call last):
            ...
            ValueError: given base_map is not compatible with im_gens
            sage: L._is_valid_homomorphism_(ZZ, [5], base_map=phi)
            Traceback (most recent call last):
            ...
            ValueError: images of some localized elements fail to be units
            sage: phi = R.hom([5], codomain=QQ)
            sage: L._is_valid_homomorphism_(ZZ, [5], base_map=phi)
            Traceback (most recent call last):
            ...
            ValueError: codomain of base_map must be Integer Ring
            sage: L._is_valid_homomorphism_(QQ, [5], base_map=phi)
            True
        """
        B = self.base_ring()
        if base_map is not None:
            if base_map.domain() is not B:
                raise ValueError('domain of base_map must be %s' % B)
            if base_map.codomain() is not codomain.base_ring():
                raise ValueError('codomain of base_map must be %s' % codomain.base_ring())
            bas_gens = B.gens()
            if im_gens and not all(base_map(g) == im_gens[bas_gens.index(g)] for g in bas_gens):
                raise ValueError('given base_map is not compatible with im_gens')
            im_gens = [base_map(g) for g in bas_gens]
            if not all(base_map(au).is_unit() for au in self._extra_units):
                raise ValueError('images of some localized elements fail to be units')
            return B._is_valid_homomorphism_(codomain, im_gens, base_map=None)
        else:
            if B._is_valid_homomorphism_(codomain, im_gens, base_map=base_map):
                phi = B.hom(im_gens, base_map=base_map)
                if not all(phi(au).is_unit() for au in self._extra_units):
                    raise ValueError('images of some localized elements fail to be units')
                return True
            return False

    def ngens(self):
        """
        Return the number of generators of ``self``
        according to the same method for the base ring.

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: Localization(R, (x**2 + 1, y - 1)).ngens()                            # needs sage.libs.pari
            2

            sage: Localization(ZZ, 2).ngens()
            1
        """
        return self.base_ring().ngens()

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self`` which is
        the ``i``-th generator of the base ring.

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: R.localization((x**2 + 1, y - 1)).gen(0)                              # needs sage.libs.pari
            x

            sage: ZZ.localization(2).gen(0)
            1
        """
        return self(self.base_ring().gen(i))

    def gens(self):
        """
        Return a tuple whose entries are the generators for this
        object, in order.

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: Localization(R, (x**2 + 1, y - 1)).gens()                             # needs sage.libs.pari
            (x, y)

            sage: Localization(ZZ, 2).gens()
            (1,)
        """
        return tuple(self(g) for g in self.base_ring().gens())

    def _cut_off_extra_units_from_base_ring_element(self, x):
        """
        Return a factor of x not divided by any additional unit of ``self``.

        INPUT:

        - ``x`` -- an element of the base ring of ``self``

        OUTPUT:

        A factor of ``x`` not divided by any additional unit of ``self`` as element
        of the base ring of ``self``.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: L = Localization(P, (x, y*z))                                         # needs sage.libs.pari
            sage: L._cut_off_extra_units_from_base_ring_element(x*y*z)                  # needs sage.libs.pari
            1
            sage: L._cut_off_extra_units_from_base_ring_element(x*z)                    # needs sage.libs.pari
            1

        TESTS:

        Check that :issue:`33463` is fixed::

            sage: L = ZZ.localization(5)
            sage: L(0).is_unit()
            False
        """
        if x.is_zero() or x.numerator().is_unit():
            # treat corner cases
            return x
        add_units = self._extra_units
        res = x
        for au in add_units:
            if au.divides(x):
                # recursion must terminate by reducing the number of factors
                res = self._cut_off_extra_units_from_base_ring_element(x // au)
                if res.is_unit():
                    return res
        return res

    def _fraction_to_element(self, x):
        """
        Check if the given element of the fraction field is contained in ``self``
        and construct it as an element of ``self`` in case the answer is true.

        INPUT:

        - ``x`` -- an element of the fraction field of the base ring

        OUTPUT: an instance of the element class of ``self`` representing `x`

        EXAMPLES::

            sage: # needs sage.libs.pari sage.libs.singular
            sage: P.<x,y,z> = QQ[]
            sage: d = x**2 + y**2 + z**2
            sage: L = Localization(P, d)
            sage: L._fraction_to_element((x+y+z)/d)
            (x + y + z)/(x^2 + y^2 + z^2)
            sage: _ in L
            True

        TESTS::

            sage: TestSuite(L).run()                                                    # needs sage.libs.pari sage.libs.singular
        """
        potential_non_unit_denom = self._cut_off_extra_units_from_base_ring_element(x.denominator())
        if potential_non_unit_denom.is_unit():
            return self.element_class(self, x)
        raise ValueError("factor %s of denominator is not a unit" % potential_non_unit_denom)

    def _coerce_map_from_(self, S):
        """
        The only things that coerce into this ring are:

        - the base ring
        - other localizations

        EXAMPLES::

            sage: # needs sage.libs.pari sage.libs.singular
            sage: P.<x,y,z> = QQ[]
            sage: L = Localization(P, y*z)
            sage: M = Localization(P, (x, y, z))
            sage: M._coerce_map_from_(L)
            True
            sage: L._coerce_map_from_(M)
            False
            sage: Q.<u, v, w> = ZZ[]
            sage: N = Localization(Q, v*w)
            sage: L._coerce_map_from_(N)
            True
            sage: N._coerce_map_from_(M)
            False
            sage: O = Localization(L, x**2 + 1)
            sage: O._coerce_map_from_(M)
            False
            sage: O._coerce_map_from_(L)
            True
        """
        if S is self.base_ring():
            return True
        elif self.base_ring().has_coerce_map_from(S):
            return True
        elif isinstance(S, Localization):
            return all(self(p).is_unit() for p in S._extra_units)

    def fraction_field(self):
        """
        Return the fraction field of ``self``.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: R.<a> = GF(5)[]
            sage: L = Localization(R, (a**2 - 3, a))
            sage: L.fraction_field()
            Fraction Field of Univariate Polynomial Ring in a over Finite Field of size 5
            sage: L.is_subring(_)
            True
        """
        return self._fraction_field

    def characteristic(self):
        """
        Return the characteristic of ``self``.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: R.<a> = GF(5)[]
            sage: L = R.localization((a**2 - 3, a))
            sage: L.characteristic()
            5
        """
        return self.base_ring().characteristic()

    def krull_dimension(self):
        """
        Return the Krull dimension of this localization.

        Since the current implementation just allows integral domains as base ring
        and localization at a finite set of elements the spectrum of ``self``
        is open in the irreducible spectrum of its base ring.
        Therefore, by density we may take the dimension from there.

        EXAMPLES::

            sage: R = ZZ.localization((2, 3))
            sage: R.krull_dimension()
            1
        """
        return self.base_ring().krull_dimension()

    def is_field(self, proof=True):
        """
        Return ``True`` if this ring is a field.

        INPUT:

        - ``proof`` -- boolean (default: ``True``); determines what to do in
          unknown cases

        ALGORITHM:

        If the parameter ``proof`` is set to ``True``, the returned value is
        correct but the method might throw an error.  Otherwise, if it is set
        to ``False``, the method returns ``True`` if it can establish that
        ``self`` is a field and ``False`` otherwise.

        EXAMPLES::

            sage: R = ZZ.localization((2, 3))
            sage: R.is_field()
            False
        """
        if proof:
            try:
                if self.krull_dimension() > 0:
                    return False
            except NotImplementedError:
                pass
        return super().is_field(proof=proof)
