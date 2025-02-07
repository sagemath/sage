# sage.doctest: needs sage.libs.pari sage.rings.padics
r"""
Overconvergent `p`-adic modular forms for small primes

This module implements computations of Hecke operators and `U_p`-eigenfunctions
on `p`-adic overconvergent modular forms of tame level 1, where `p` is one of
the primes `\{2, 3, 5, 7, 13\}`, using the algorithms described in [Loe2007]_.

- [Loe2007]_

AUTHORS:

- David Loeffler (August 2008): initial version
- David Loeffler (March 2009): extensively reworked
- Lloyd Kilford (May 2009): add
  :meth:`~sage.modular.overconvergent.genus0.OverconvergentModularFormsSpace.slopes`
  method
- David Loeffler (June 2009): miscellaneous bug fixes and usability improvements

The Theory
~~~~~~~~~~

Let `p` be one of the above primes, so `X_0(p)` has genus 0, and let

.. MATH::

    f_p = \sqrt[p-1]{\frac{\Delta(pz)}{\Delta(z)}}

(an `\eta`-product of level `p` -- see module :mod:`sage.modular.etaproducts`).
Then one can show that `f_p` gives an isomorphism `X_0(p) \to \mathbb{P}^1`.
Furthermore, if we work over `\CC_p`, the `r`-overconvergent locus on `X_0(p)`
(or of `X_0(1)`, via the canonical subgroup lifting), corresponds to the
`p`-adic disc

.. MATH::

    |f_p|_p \le p^{\frac{12r}{p-1}}.

(This is Theorem 1 of [Loe2007]_.)

Hence if we fix an element `c` with `|c| = p^{-\frac{12r}{p-1}}`, the space
`S_k^\dagger(1, r)` of overconvergent `p`-adic modular forms has an orthonormal
basis given by the functions `(cf)^n`.  So any element can be written in the
form `E_k \times \sum_{n \ge 0} a_n (cf)^n`, where `a_n \to 0` as `N \to
\infty`, and any such sequence `a_n` defines a unique overconvergent form.

One can now find the matrix of Hecke operators in this basis, either by
calculating `q`-expansions, or (for the special case of `U_p`) using a
recurrence formula due to Kolberg.

An Extended Example
~~~~~~~~~~~~~~~~~~~

We create a space of 3-adic modular forms::

    sage: M = OverconvergentModularForms(3, 8, 1/6, prec=60)

Creating an element directly as a linear combination of basis vectors.

.. link

::

    sage: f1 = M.3 + M.5; f1.q_expansion()
    27*q^3 + 1055916/1093*q^4 + 19913121/1093*q^5 + 268430112/1093*q^6 + ...
    sage: f1.coordinates(8)
    [0, 0, 0, 1, 0, 1, 0, 0]

We can coerce from elements of classical spaces of modular forms:

.. link

::

    sage: f2 = M(CuspForms(3, 8).0); f2
    3-adic overconvergent modular form of weight-character 8 with q-expansion
     q + 6*q^2 - 27*q^3 - 92*q^4 + 390*q^5 - 162*q^6 ...

We express this in a basis, and see that the coefficients go to zero very fast:

.. link

::

    sage: [x.valuation(3) for x in f2.coordinates(60)]
    [+Infinity, -1, 3, 6, 10, 13, 18, 20, 24, 27, 31, 34, 39, 41, 45, 48, 52, 55, 61,
     62, 66, 69, 73, 76, 81, 83, 87, 90, 94, 97, 102, 104, 108, 111, 115, 118, 124, 125,
     129, 132, 136, 139, 144, 146, 150, 153, 157, 160, 165, 167, 171, 174, 178, 181,
     188, 188, 192, 195, 199, 202]

This form has more level at `p`, and hence is less overconvergent:

.. link

::

    sage: f3 = M(CuspForms(9, 8).0); [x.valuation(3) for x in f3.coordinates(60)]
    [+Infinity, -1, -1, 0, -4, -4, -2, -3, 0, 0, -1, -1, 1, 0, 3, 3, 3, 3, 5, 3, 7, 7,
     6, 6, 8, 7, 10, 10, 8, 8, 10, 9, 12, 12, 12, 12, 14, 12, 17, 16, 15, 15, 17, 16,
     19, 19, 18, 18, 20, 19, 22, 22, 22, 22, 24, 21, 25, 26, 24, 24]

An error will be raised for forms which are not sufficiently overconvergent:

.. link

::

    sage: M(CuspForms(27, 8).0)
    Traceback (most recent call last):
    ...
    ValueError: Form is not overconvergent enough (form is only 1/12-overconvergent)

Let's compute some Hecke operators. Note that the coefficients of this matrix are `p`-adically tiny:

.. link

::

    sage: M.hecke_matrix(3, 4).change_ring(Qp(3, prec=1))
    [        1 + O(3)                0                0                0]
    [               0   2*3^3 + O(3^4)   2*3^3 + O(3^4)     3^2 + O(3^3)]
    [               0   2*3^7 + O(3^8)   2*3^8 + O(3^9)     3^6 + O(3^7)]
    [               0 2*3^10 + O(3^11) 2*3^10 + O(3^11)  2*3^9 + O(3^10)]

We compute the eigenfunctions of a 4x4 truncation:

.. link

::

    sage: efuncs = M.eigenfunctions(4)
    sage: for i in [1..3]:
    ....:     print(efuncs[i].q_expansion(prec=4).change_ring(Qp(3, prec=20)))
    (1 + O(3^20))*q
     + (2*3 + 3^15 + 3^16 + 3^17 + 2*3^19 + 2*3^20 + O(3^21))*q^2
     + (2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9
        + 2*3^10 + 2*3^11 + 2*3^12 + 2*3^13 + 2*3^14 + 2*3^15
        + 2*3^16 + 3^17 + 2*3^18 + 2*3^19 + 3^21 + 3^22 + O(3^23))*q^3
     + O(q^4)
    (1 + O(3^20))*q
     + (3 + 2*3^2 + 3^3 + 3^4 + 3^12 + 3^13 + 2*3^14
         + 3^15 + 2*3^17 + 3^18 + 3^19 + 3^20 + O(3^21))*q^2
     + (3^7 + 3^13 + 2*3^14 + 2*3^15 + 3^16 + 3^17 + 2*3^18
         + 3^20 + 2*3^21 + 2*3^22 + 2*3^23 + 2*3^25 + O(3^27))*q^3
     + O(q^4)
    (1 + O(3^20))*q
     + (2*3 + 3^3 + 2*3^4 + 3^6 + 2*3^8 + 3^9 + 3^10
         + 2*3^11 + 2*3^13 + 3^16 + 3^18 + 3^19 + 3^20 + O(3^21))*q^2
     + (3^9 + 2*3^12 + 3^15 + 3^17 + 3^18 + 3^19 + 3^20
         + 2*3^22 + 2*3^23 + 2*3^27 + 2*3^28 + O(3^29))*q^3
     + O(q^4)

The first eigenfunction is a classical cusp form of level 3:

.. link

::

    sage: (efuncs[1] - M(CuspForms(3, 8).0)).valuation()
    13

The second is an Eisenstein series!

.. link

::

    sage: (efuncs[2] - M(EisensteinForms(3, 8).1)).valuation()
    10

The third is a genuinely new thing (not a classical modular form at all); the
coefficients are almost certainly not algebraic over `\QQ`. Note that the slope
is 9, so Coleman's classicality criterion (forms of slope `< k-1` are
classical) does not apply.

.. link

::

    sage: a3 = efuncs[3].q_expansion()[3]; a3
    3^9 + 2*3^12 + 3^15 + 3^17 + 3^18 + 3^19 + 3^20 + 2*3^22 + 2*3^23 + 2*3^27
     + 2*3^28 + 3^32 + 3^33 + 2*3^34 + 3^38 + 2*3^39 + 3^40 + 2*3^41 + 3^44 + 3^45
     + 3^46 + 2*3^47 + 2*3^48 + 3^49 + 3^50 + 2*3^51 + 2*3^52 + 3^53 + 2*3^54 + 3^55
     + 3^56 + 3^57 + 2*3^58 + 2*3^59 + 3^60 + 2*3^61 + 2*3^63 + 2*3^64 + 3^65 + 2*3^67
     + 3^68 + 2*3^69 + 2*3^71 + 3^72 + 2*3^74 + 3^75 + 3^76 + 3^79 + 3^80 + 2*3^83
     + 2*3^84 + 3^85 + 2*3^87 + 3^88 + 2*3^89 + 2*3^90 + 2*3^91 + 3^92 + O(3^98)
    sage: efuncs[3].slope()
    9
"""
# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                     2008-9 David Loeffler <d.loeffler.01@cantab.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from typing import Iterator

import weakref

import sage.rings.abc

from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.special import diagonal_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.verbose import verbose
from sage.modular.arithgroup.all import Gamma0_class, Gamma1_class
from sage.modular.dirichlet import trivial_character
from sage.modular.etaproducts import EtaProduct
from sage.modular.modform.element import ModularFormElement
from sage.modular.modform.hecke_operator_on_qexp import hecke_operator_on_qexp
from sage.modules.free_module_element import vector
from sage.modules.module import Module
from sage.rings.big_oh import O
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.element import Vector, ModuleElement
from sage.structure.richcmp import richcmp
from .weightspace import WeightSpace_constructor as WeightSpace, WeightCharacter

lazy_import('sage.modular.modform.j_invariant', 'j_invariant_qexp')
lazy_import('sage.rings.padics.factory', 'Qp', as_='pAdicField')


__ocmfdict = {}


####################
# Factory function #
####################

def OverconvergentModularForms(prime, weight, radius, base_ring=QQ, prec=20, char=None):
    r"""
    Create a space of overconvergent `p`-adic modular forms of level
    `\Gamma_0(p)`, over the given base ring. The base ring need not be a
    `p`-adic ring (the spaces we compute with typically have bases over
    `\QQ`).

    INPUT:

    - ``prime`` -- a prime number `p`, which must be one of the primes `\{2, 3,
      5, 7, 13\}`, or the congruence subgroup `\Gamma_0(p)` where `p` is one of
      these primes

    - ``weight`` -- integer (which at present must be 0 or `\ge 2`), the
      weight

    - ``radius`` -- a rational number in the interval `\left( 0, \frac{p}{p+1}
      \right)`, the radius of overconvergence

    - ``base_ring`` -- (default: `\QQ`), a ring over which to compute; this
      need not be a `p`-adic ring

    - ``prec`` -- integer (default: 20); the number of `q`-expansion terms to
      compute

    - ``char`` -- a Dirichlet character modulo `p` or ``None`` (the default);
      here ``None`` is interpreted as the trivial character modulo `p`

    The character `\chi` and weight `k` must satisfy `(-1)^k = \chi(-1)`, and
    the base ring must contain an element `v` such that
    `{\rm ord}_p(v) = \frac{12 r}{p-1}` where `r` is the radius of
    overconvergence (and `{\rm ord}_p` is normalised so `{\rm ord}_p(p) = 1`).

    EXAMPLES::

        sage: OverconvergentModularForms(3, 0, 1/2)
        Space of 3-adic 1/2-overconvergent modular forms
         of weight-character 0 over Rational Field
        sage: OverconvergentModularForms(3, 16, 1/2)
        Space of 3-adic 1/2-overconvergent modular forms
         of weight-character 16 over Rational Field
        sage: OverconvergentModularForms(3, 3, 1/2, char=DirichletGroup(3,QQ).0)
        Space of 3-adic 1/2-overconvergent modular forms
         of weight-character (3, 3, [-1]) over Rational Field
    """
    if isinstance(prime, (Gamma0_class, Gamma1_class)):
        prime = prime.level()
    else:
        prime = ZZ(prime)
    if char is None:
        char = trivial_character(prime, base_ring=QQ)
    if int(prime) not in [2, 3, 5, 7, 13]:
        raise ValueError("p must be one of {2, 3, 5, 7, 13}")
    key = (prime, weight, radius, base_ring, prec, char)
    if key in __ocmfdict:
        w = __ocmfdict[key]
        M = w()
        if M is not None:
            return M
    M = OverconvergentModularFormsSpace(*key)
    __ocmfdict[key] = weakref.ref(M)
    return M


#########################
# Main class definition #
#########################

class OverconvergentModularFormsSpace(Module):
    r"""
    A space of overconvergent modular forms of level `\Gamma_0(p)`,
    where `p` is a prime such that `X_0(p)` has genus 0.

    Elements are represented as power series, with a formal power series `F`
    corresponding to the modular form `E_k^\ast \times F(g)` where `E_k^\ast`
    is the `p`-deprived Eisenstein series of weight-character `k`, and `g` is a
    uniformiser of `X_0(p)` normalised so that the `r`-overconvergent region
    `X_0(p)_{\ge r}` corresponds to `|g| \le 1`.

    TESTS::

        sage: x = polygen(ZZ, 'x')
        sage: K.<w> = Qp(13).extension(x^2 - 13); M = OverconvergentModularForms(13, 20, radius=1/2, base_ring=K)
        sage: M is loads(dumps(M))
        True
    """

    ###############
    # Init script #
    ###############

    def __init__(self, prime, weight, radius, base_ring, prec, char):
        r"""
        Create a space of overconvergent `p`-adic modular forms of level
        `\Gamma_0(p)`, over the given base ring.

        The base ring need not be a
        `p`-adic ring (the spaces we compute with typically have bases over
        `\QQ`).

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2); M
            Space of 3-adic 1/2-overconvergent modular forms of weight-character 0 over Rational Field
            sage: TestSuite(M).run()
        """
        self._p = prime

        if not (base_ring == QQ or isinstance(base_ring, sage.rings.abc.pAdicField)):
            raise TypeError("Base ring must be QQ or a p-adic field")

        if base_ring != QQ and base_ring.prime() != self._p:
            raise TypeError("Residue characteristic of base ring (=%s) must be %s" % (base_ring, self._p))

        if isinstance(weight, WeightCharacter):
            self._wtchar = weight
        else:
            self._wtchar = WeightSpace(prime, base_ring=char.base_ring())(weight, char, algebraic=True)

        if not self._wtchar.is_even():
            raise ValueError("Weight-character must be even")

        Module.__init__(self, base_ring)

        self._prec = prec

        self._qsr = PowerSeriesRing(base_ring, 'q', prec)  # q-series ring
        self._gsr = PowerSeriesRing(base_ring, 'g', prec)  # g-adic expansions, g = c*f

        self._cached_recurrence_matrix = None
        self._set_radius(radius)
        self._basis_cache = [self._wtchar.pAdicEisensteinSeries(self._qsr, self.prec())]
        self._uniformiser = self._qsr(EtaProduct(prime, {prime: 24/ZZ(prime-1), ZZ(1): -24/ZZ(prime-1)}).qexp(self.prec()))

        for i in range(1, self.prec()):
            self._basis_cache.append(self._basis_cache[-1] * self._uniformiser * self._const)

    #####################################
    # Methods called by the init script #
    #####################################

    def _set_radius(self, radius):
        r"""
        Set the radius of overconvergence to be `r`, where `r` is a rational
        number in the interval `0 < r < \frac{p}{p+1}`.

        This only makes sense if the base ring contains an element of
        normalised valuation `\frac{12r}{p-1}`. If this valuation is an
        integer, we use the appropriate power of `p`. Otherwise, we assume the
        base ring has a ``uniformiser`` method and take an appropriate power of
        the uniformiser, raising an error if no such element exists.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 2, 1/2) # indirect doctest
            sage: M._set_radius(1/3); M
            Space of 3-adic 1/3-overconvergent modular forms of weight-character 2
             over Rational Field

            sage: x = polygen(ZZ, 'x')
            sage: L.<w> = Qp(3).extension(x^5 - 3)
            sage: OverconvergentModularForms(3, 2, 1/30, base_ring=L).normalising_factor()  # indirect doctest
            w + O(w^101)

            sage: OverconvergentModularForms(3, 2, 1/40, base_ring=L)
            Traceback (most recent call last):
            ...
            ValueError: no element of base ring (=3-adic Eisenstein Extension ...)
            has normalised valuation 3/20
        """

        p = ZZ(self.prime())

        if (radius < 0 or radius > p/(p+1)):
            raise ValueError("radius (=%s) must be between 0 and p/(p+1)" % radius)
        d = 12/(p-1)*radius
        if d.is_integral():
            self._const = p ** ZZ(d)
            self._radius = radius
        else:
            try:
                pi = self.base_ring().uniformiser()
                e = d / pi.normalized_valuation()
            except AttributeError:  # base ring isn't a p-adic ring
                pi = p
                e = d
            if not e.is_integral():
                raise ValueError("no element of base ring (=%s) has normalised valuation %s" % (self.base_ring(), radius * 12 / (p-1)))
            self._radius = radius
            self._const = pi ** ZZ(e)

    ##############################################
    # Boring functions that access internal data #
    ##############################################

    def is_exact(self):
        r"""
        Return ``True`` if elements of this space are represented exactly.

        This would mean that there is no precision loss when doing arithmetic.
        As this is never true for overconvergent modular forms spaces,
        this method returns ``False``.

        EXAMPLES::

            sage: OverconvergentModularForms(13, 12, 0).is_exact()
            False
        """
        return False

    def change_ring(self, ring):
        r"""
        Return the space corresponding to ``self`` but over the given base ring.

        EXAMPLES::

            sage: M = OverconvergentModularForms(2, 0, 1/2)
            sage: M.change_ring(Qp(2))
            Space of 2-adic 1/2-overconvergent modular forms of weight-character 0
             over 2-adic Field with ...
        """
        return OverconvergentModularForms(self.prime(), self.weight(), self.radius(), ring, self.prec(), self.character())

    def base_extend(self, ring):
        r"""
        Return the base extension of ``self`` to the given base ring.

        There must be a canonical map to this ring from the current
        base ring, otherwise a :exc:`TypeError` will be raised.

        EXAMPLES::

            sage: M = OverconvergentModularForms(2, 0, 1/2, base_ring=Qp(2))
            sage: x = polygen(ZZ, 'x')
            sage: M.base_extend(Qp(2).extension(x^2 - 2, names='w'))
            Space of 2-adic 1/2-overconvergent modular forms of weight-character 0
             over 2-adic Eisenstein Extension ...
            sage: M.base_extend(QQ)
            Traceback (most recent call last):
            ...
            TypeError: Base extension of self (over '2-adic Field with capped
            relative precision 20') to ring 'Rational Field' not defined.
        """
        if ring.has_coerce_map_from(self.base_ring()):
            return self.change_ring(ring)
        else:
            raise TypeError("Base extension of self (over '%s') to ring '%s' not defined." % (self.base_ring(), ring))

    def _an_element_(self):
        r"""
        Return an element of this space (used by the coercion machinery).

        EXAMPLES::

            sage: OverconvergentModularForms(3, 2, 1/3, prec=4).an_element()  # indirect doctest
            3-adic overconvergent modular form of weight-character 2
             with q-expansion 9*q + 216*q^2 + 2430*q^3 + O(q^4)
        """
        return self.element_class(self, self._gsr.an_element())

    def character(self):
        r"""
        Return the character of ``self``.

        For overconvergent forms, the weight and the character are unified into
        the concept of a weight-character, so this returns exactly the same
        thing as :meth:`weight`.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 0, 1/2).character()
            0
            sage: type(OverconvergentModularForms(3, 0, 1/2).character())
            <class '...weightspace.AlgebraicWeight'>
            sage: OverconvergentModularForms(3, 3, 1/2, char=DirichletGroup(3,QQ).0).character()
            (3, 3, [-1])
        """
        return self._wtchar

    def weight(self):
        r"""
        Return the weight of ``self``.

        For overconvergent forms, the weight and the character are unified into
        the concept of a weight-character, so this returns exactly the same
        thing as :meth:`character`.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 0, 1/2).weight()
            0
            sage: type(OverconvergentModularForms(3, 0, 1/2).weight())
            <class '...weightspace.AlgebraicWeight'>
            sage: OverconvergentModularForms(3, 3, 1/2, char=DirichletGroup(3,QQ).0).weight()
            (3, 3, [-1])
        """
        return self._wtchar

    def normalising_factor(self):
        r"""
        Return the normalising factor of ``self``.

        The normalising factor `c` such that `g = c f` is a parameter for the
        `r`-overconvergent disc in `X_0(p)`, where `f` is the standard
        uniformiser.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<w> = Qp(7).extension(x^2 - 7)
            sage: OverconvergentModularForms(7, 0, 1/4, base_ring=L).normalising_factor()
            w + O(w^41)
        """
        return self._const

    def __eq__(self, other):
        r"""
        Check whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 12, 1/2) == ModularForms(3, 12)
            False
            sage: OverconvergentModularForms(3, 0, 1/2) == OverconvergentModularForms(3, 0, 1/3)
            False
            sage: OverconvergentModularForms(3, 0, 1/2) == OverconvergentModularForms(3, 0, 1/2, base_ring=Qp(3))
            False
            sage: OverconvergentModularForms(3, 0, 1/2) == OverconvergentModularForms(3, 0, 1/2)
            True
        """
        if not isinstance(other, OverconvergentModularFormsSpace):
            return False
        else:
            return self._params() == other._params()

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 12, 1/2) != ModularForms(3, 12)
            True
            sage: OverconvergentModularForms(3, 0, 1/2) != OverconvergentModularForms(3, 0, 1/3)
            True
            sage: OverconvergentModularForms(3, 0, 1/2) != OverconvergentModularForms(3, 0, 1/2, base_ring=Qp(3))
            True
            sage: OverconvergentModularForms(3, 0, 1/2) != OverconvergentModularForms(3, 0, 1/2)
            False
        """
        return not (self == other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: h1 = hash(OverconvergentModularForms(3, 12, 1/2))
            sage: h2 = hash(OverconvergentModularForms(3, 12, 1/2))
            sage: h3 = hash(OverconvergentModularForms(3, 0, 1/2))
            sage: h1 == h2 and h1 != h3
            True
        """
        return hash(self._params())

    def _params(self):
        r"""
        Return the parameters that define this module uniquely: prime, weight,
        character, radius of overconvergence and base ring. Mostly used for
        pickling.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<w> = Qp(7).extension(x^2 - 7)
            sage: OverconvergentModularForms(7, 0, 1/4, base_ring=L)._params()
            (7,
             0,
             1/4,
             7-adic Eisenstein Extension Field in w defined by x^2 - 7,
             20,
             Dirichlet character modulo 7 of conductor 1 mapping 3 |--> 1)
        """
        return (self.prime(), self.weight().k(), self.radius(), self.base_ring(), self.prec(), self.weight().chi())

    def __reduce__(self):
        r"""
        Return the function and arguments used to construct ``self``. Used for pickling.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<w> = Qp(7).extension(x^2 - 7)
            sage: OverconvergentModularForms(7, 0, 1/4, base_ring=L).__reduce__()
            (<function OverconvergentModularForms at ...>,
             (7,
              0,
              1/4,
              7-adic Eisenstein Extension Field in w defined by x^2 - 7,
              20,
              Dirichlet character modulo 7 of conductor 1 mapping 3 |--> 1))
        """
        return (OverconvergentModularForms, self._params())

    def gen(self, i):
        r"""
        Return the `i`-th module generator of ``self``.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 2, 1/2, prec=4)
            sage: M.gen(0)
            3-adic overconvergent modular form of weight-character 2
             with q-expansion 1 + 12*q + 36*q^2 + 12*q^3 + O(q^4)
            sage: M.gen(1)
            3-adic overconvergent modular form of weight-character 2
             with q-expansion 27*q + 648*q^2 + 7290*q^3 + O(q^4)
            sage: M.gen(30)
            3-adic overconvergent modular form of weight-character 2
             with q-expansion O(q^4)
        """

        return OverconvergentModularFormElement(self, gexp=self._gsr.gen()**i)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 0, 1/2)._repr_()
            'Space of 3-adic 1/2-overconvergent modular forms of weight-character 0 over Rational Field'
        """
        return "Space of %s-adic %s-overconvergent modular forms of weight-character %s over %s" % (self.prime(), self.radius(), self.weight(), self.base_ring())

    def prime(self):
        r"""
        Return the residue characteristic of ``self``.

        This is the prime `p` such that this is a `p`-adic space.

        EXAMPLES::

            sage: OverconvergentModularForms(5, 12, 1/3).prime()
            5
        """
        return self._p

    def radius(self):
        r"""
        The radius of overconvergence of this space.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 0, 1/3).radius()
            1/3
        """
        return self._radius

    def gens(self) -> Iterator:
        r"""
        Return a generator object that iterates over the (infinite) set of
        basis vectors of ``self``.

        EXAMPLES::

            sage: o = OverconvergentModularForms(3, 12, 1/2)
            sage: t = o.gens()
            sage: next(t)
            3-adic overconvergent modular form of weight-character 12 with q-expansion
            1 - 32760/61203943*q - 67125240/61203943*q^2 - ...
            sage: next(t)
            3-adic overconvergent modular form of weight-character 12 with q-expansion
            27*q + 19829193012/61203943*q^2 + 146902585770/61203943*q^3 + ...
        """
        i = 0
        while True:
            yield self.gen(i)
            i += 1

    def prec(self):
        r"""
        Return the series precision of ``self``.

        Note that this is different from the `p`-adic precision of the base ring.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 0, 1/2).prec()
            20
            sage: OverconvergentModularForms(3, 0, 1/2, prec=40).prec()
            40
        """
        return self._prec

    #####################################
    # Element construction and coercion #
    #####################################

    def _element_constructor_(self, input):
        r"""
        Create an element of this space. Allowable inputs are:

        - elements of compatible spaces of modular forms or overconvergent
          modular forms

        - arbitrary power series in `q`

        - lists of elements of the base ring (interpreted as vectors in the
          basis given by :meth:`gens`).

        Precision may be specified by padding lists at the end with zeros;
        inputs with a higher precision than the set precision of this space
        will be rounded.

        EXAMPLES:

        From a `q`-expansion::

            sage: M = OverconvergentModularForms(3, 0, 1/2, prec=5)
            sage: R.<q> = QQ[[]]
            sage: f = M(q + q^2 - q^3 + O(q^16)); f
            3-adic overconvergent modular form of weight-character 0 with q-expansion q + q^2 - q^3 + O(q^5)
            sage: M.coordinate_vector(f)
            (0, 1/27, -11/729, 173/19683, -3172/531441)

        From a list or a vector::

            sage: M([1,0,1])
            3-adic overconvergent modular form of weight-character 0 with q-expansion 1 + 729*q^2 + O(q^3)
            sage: M([1,0,1,0,0])
            3-adic overconvergent modular form of weight-character 0 with q-expansion 1 + 729*q^2 + 17496*q^3 + 236196*q^4 + O(q^5)
            sage: f = M([1,0,1,0,0]); v = M.coordinate_vector(f); v
            (1, 0, 1, 0, 0)
            sage: M(v) == f
            True

        From a classical modular form::

            sage: f = CuspForms(Gamma0(3), 12).0; f
            q - 176*q^4 + 2430*q^5 + O(q^6)
            sage: fdag = OverconvergentModularForms(3, 12, 1/3, prec=8)(f); fdag
            3-adic overconvergent modular form of weight-character 12 with q-expansion q - 176*q^4 + 2430*q^5 - 5832*q^6 - 19336*q^7 + O(q^8)
            sage: fdag.parent().coordinate_vector(f)*(1 + O(3^2))
            (0, 3^-2 + O(3^0), 2*3^-3 + 2*3^-2 + O(3^-1), 3^-4 + 3^-3 + O(3^-2), 2 + 3 + O(3^2), 2*3 + 3^2 + O(3^3), 2*3^4 + 2*3^5 + O(3^6), 3^5 + 3^6 + O(3^7))
            sage: OverconvergentModularForms(3, 6, 1/3)(f)
            Traceback (most recent call last):
            ...
            TypeError: Cannot create an element of 'Space of 3-adic ...' from element of incompatible space 'Cuspidal subspace ...'

        We test that zero elements are handled properly::

            sage: M(0)
            3-adic overconvergent modular form of weight-character 0 with q-expansion O(q^5)
            sage: M(O(q^3))
            3-adic overconvergent modular form of weight-character 0 with q-expansion O(q^3)

        We test coercion between spaces of different precision::

            sage: M10 = OverconvergentModularForms(3, 0, 1/2, prec=10)
            sage: f = M10.1
            sage: M(f)
            3-adic overconvergent modular form of weight-character 0 with q-expansion 27*q + 324*q^2 + 2430*q^3 + 13716*q^4 + O(q^5)
            sage: M10(M(f))
            3-adic overconvergent modular form of weight-character 0 with q-expansion 27*q + 324*q^2 + 2430*q^3 + 13716*q^4 + O(q^5)
        """
        if isinstance(input, int):
            input = ZZ(input)

        if isinstance(input, OverconvergentModularFormElement):
            return self._coerce_from_ocmf(input)

        elif isinstance(input, ModularFormElement):
            if ((input.level() == 1 or input.level().prime_factors() == [self.prime()])
                    and input.weight() == self.weight().k()
                    and input.character().primitive_character() == self.weight().chi().primitive_character()):
                p = ZZ(self.prime())
                nu = (input.level() == 1 and p/(p+1)) or (1 / (p + 1) * p**(2 - input.level().valuation(p)))
                if self.radius() > nu:
                    raise ValueError("Form is not overconvergent enough (form is only %s-overconvergent)" % nu)
                else:
                    return self(self._qsr(input.q_expansion(self.prec())))
            else:
                raise TypeError("Cannot create an element of '%s' from element of incompatible space '%s'" % (self, input.parent()))

        elif isinstance(input, (list, tuple, Vector)):
            v = list(input)
            n = len(v)
            return OverconvergentModularFormElement(self, gexp=self._gsr(v).add_bigoh(n), qexp=None)

        elif self._qsr.has_coerce_map_from(input.parent()):
            return OverconvergentModularFormElement(self, gexp=None, qexp=self._qsr(input))

        else:
            raise TypeError("Don't know how to create an overconvergent modular form from %s" % input)

    @cached_method
    def zero(self):
        """
        Return the zero of this space.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<w> = Qp(13).extension(x^2 - 13)
            sage: M = OverconvergentModularForms(13, 20, radius=1/2, base_ring=K)
            sage: K.zero()
            0
        """
        return self(0)

    def _coerce_from_ocmf(self, f):
        r"""
        Try to convert the overconvergent modular form `f` into an element of ``self``.

        An error will be raised if this is obviously nonsense.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2)
            sage: MM = M.base_extend(Qp(3))
            sage: R.<q> = Qp(3)[[]]; f = MM(q + O(q^2)); f
            3-adic overconvergent modular form of weight-character 0 with q-expansion (1 + O(3^20))*q + O(q^2)
            sage: M._coerce_from_ocmf(f)
            3-adic overconvergent modular form of weight-character 0 with q-expansion q + O(q^2)
            sage: f in M # indirect doctest
            True
        """
        prime, weight, radius, base_ring, prec, char = f.parent()._params()
        if (prime, weight, char) != (self.prime(), self.weight().k(), self.weight().chi()):
            raise TypeError("Cannot create an element of '%s' from element of incompatible space '%s'" % (self, input.parent()))
        return self(self._qsr(f.q_expansion()))

    def _coerce_map_from_(self, other):
        r"""
        Canonical coercion of ``x`` into ``self``.

        Here the possibilities for ``x`` are more restricted.

        TESTS::

            sage: M = OverconvergentModularForms(3, 0, 1/2)
            sage: MM = M.base_extend(Qp(3))
            sage: MM.has_coerce_map_from(M) # indirect doctest
            True
            sage: MM.coerce(M.1)
            3-adic overconvergent modular form of weight-character 0 with q-expansion (3^3 + O(3^23))*q + (3^4 + 3^5 + O(3^24))*q^2 ...
            sage: M.has_coerce_map_from(MM)
            False
            sage: M.coerce(1)
            3-adic overconvergent modular form of weight-character 0 with q-expansion 1 + O(q^20)
        """
        if (isinstance(other, OverconvergentModularFormsSpace) and
                self.base_ring().has_coerce_map_from(other.base_ring())):
            return True
        return self.base_ring().has_coerce_map_from(other)

    def coordinate_vector(self, x):
        r"""
        Write ``x`` as a vector with respect to the basis given by ``self.basis()``.

        Here ``x`` must be an element of this space or something that can be
        converted into one. If ``x`` has precision less than the default precision
        of ``self``, then the returned vector will be shorter.

        EXAMPLES::

            sage: M = OverconvergentModularForms(Gamma0(3), 0, 1/3, prec=4)
            sage: M.coordinate_vector(M.gen(2))
            (0, 0, 1, 0)
            sage: q = QQ[['q']].gen(); M.coordinate_vector(q - q^2 + O(q^4))
            (0, 1/9, -13/81, 74/243)
            sage: M.coordinate_vector(q - q^2 + O(q^3))
            (0, 1/9, -13/81)
        """
        if hasattr(x, 'base_ring') and x.base_ring() != self.base_ring():
            return self.base_extend(x.base_ring()).coordinate_vector(x)

        if x.parent() != self:
            x = self(x)

        return vector(self.base_ring(), x.gexp().padded_list(x.gexp().prec()))

    ##########################################################
    # Pointless routines required by parent class definition #
    ##########################################################

    def ngens(self):
        r"""
        The number of generators of ``self`` (as a module over its base ring), i.e. infinity.

        EXAMPLES::

            sage: M = OverconvergentModularForms(2, 4, 1/6)
            sage: M.ngens()
            +Infinity
        """
        return Infinity

    #####################################
    # Routines with some actual content #
    #####################################

    def hecke_operator(self, f, m):
        r"""
        Given an element `f` and an integer `m`, calculates the Hecke operator
        `T_m` acting on `f`.

        The input may be either a "bare" power series, or an
        :class:`OverconvergentModularFormElement` object; the return value will be of
        the same type.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2)
            sage: f = M.1
            sage: M.hecke_operator(f, 3)
            3-adic overconvergent modular form of weight-character 0 with q-expansion
            2430*q + 265356*q^2 + 10670373*q^3 + 249948828*q^4 + 4113612864*q^5
             + 52494114852*q^6 + O(q^7)
            sage: M.hecke_operator(f.q_expansion(), 3)
            2430*q + 265356*q^2 + 10670373*q^3 + 249948828*q^4 + 4113612864*q^5
             + 52494114852*q^6 + O(q^7)
        """

        # This should just be an instance of hecke_operator_on_qexp but that
        # won't accept arbitrary power series as input, although it's clearly
        # supposed to, which seems rather to defy the point but never mind...

        if f.parent() is self:
            return self(self.hecke_operator(f.q_expansion(), m))
        elif isinstance(f, OverconvergentModularFormElement):
            if f.parent() is self.base_extend(f.parent().base_ring()):
                return f.parent().hecke_operator(f, m)
            else:
                raise TypeError("Not an element of this space")
        else:
            return hecke_operator_on_qexp(f, m, self.weight().k(), eps=self.weight().chi())

    def _convert_to_basis(self, qexp):
        r"""
        Given a `q`-expansion, converts it to a vector in the basis of this
        space, to the maximum possible precision (which is the minimum of the
        `q`-adic precision of the `q`-expansion and the precision of ``self``).

        EXAMPLES::

            sage: M = OverconvergentModularForms(2, 0, 1/2)
            sage: R.<q> = QQ[[]]
            sage: M._convert_to_basis(q + q^2 + O(q^4))
            1/64*g - 23/4096*g^2 + 201/65536*g^3 + O(g^4)
        """
        n = min(qexp.prec(), self.prec())
        x = qexp
        g = self._gsr.gen()
        answer = self._gsr(0)
        for i in range(n):
            assert x.valuation() >= i
            answer += (x[i] / self._basis_cache[i][i]) * g**i
            x = x - self._basis_cache[i] * answer[i]
        return answer + O(g**n)

    def hecke_matrix(self, m, n, use_recurrence=False, exact_arith=False, side='left'):
        r"""
        Calculate the matrix of the `T_m` operator, truncated to `n \times n`.

        INPUT:

        - ``m`` -- integer; determines the operator `T_m`

        - ``n`` -- integer; truncate the matrix in the basis of this space
          to an `n \times n` matrix

        - ``use_recurrence`` -- boolean (default: ``False``); whether to use
          Kolberg style recurrences. If ``False``, use naive `q`-expansion
          arguments.

        - ``exact_arith`` -- boolean (default: ``True``); whether to do the
          computation to be done with rational arithmetic, even if the base ring
          is an inexact `p`-adic ring.

          This is useful as there can be precision loss issues (particularly
          with ``use_recurrence=False``).

        - ``side`` -- ``'left'`` (default) or ``'right'``; if ``'left'``, the
          operator acts on the left on column vectors

        EXAMPLES::

            sage: OverconvergentModularForms(2, 0, 1/2).hecke_matrix(2, 4)
            [    1     0     0     0]
            [    0    24    64     0]
            [    0    32  1152  4608]
            [    0     0  3072 61440]
            sage: o = OverconvergentModularForms(2, 12, 1/2, base_ring=pAdicField(2))
            sage: o.hecke_matrix(2, 3) * (1 + O(2^2))
            [        1 + O(2^2)                  0                  0]
            [                 0       2^3 + O(2^5)       2^6 + O(2^8)]
            [                 0       2^4 + O(2^6) 2^7 + 2^8 + O(2^9)]
            sage: o = OverconvergentModularForms(2, 12, 1/2, base_ring=pAdicField(2))
            sage: o.hecke_matrix(2, 3, exact_arith=True)
            [                             1                              0                              0]
            [                             0               33881928/1414477                             64]
            [                             0 -192898739923312/2000745183529             1626332544/1414477]

        Side switch::

            sage: OverconvergentModularForms(2, 0, 1/2).hecke_matrix(2, 4, side='right')
            [    1     0     0     0]
            [    0    24    32     0]
            [    0    64  1152  3072]
            [    0     0  4608 61440]
        """
        if exact_arith and not self.base_ring().is_exact():
            return self.change_ring(QQ).hecke_matrix(m, n, use_recurrence, side=side)

        M = MatrixSpace(self.base_ring(), n)
        mat = M(0)
        for j in range(min(n, self.prime())):
            l = self._convert_to_basis(self.hecke_operator(self._basis_cache[j], m))
            for i in range(n):
                try:
                    mat[i, j] = l[i]
                except IndexError:
                    if not self.weight().is_zero():
                        raise ValueError("n is too large for current precision")
                    else:
                        if i <= self.prime() * j:
                            raise ValueError("n is too large computing initial conds: can't work out u[%s, %s]" % (i, j))
                        else:
                            mat[i, j] = 0  # computations are exact for weight 0, and we know these terms are zero
        if use_recurrence:
            if m != self.prime():
                raise ValueError("Recurrence method not valid when m != p")
            for j in range(self.prime(), n):
                # can only apply recurrence if have i,j both >= p.
                if j >= self.prec():
                    for i in range(self.prime()):
                        if self.weight() != 0:
                            raise ValueError("n is too large for current precision")
                        else:
                            if j <= self.prime() * i:
                                raise ValueError("n is too large computing initial conds: can't work out u[%s,%s]" % (i, j))
                        mat[i, j] = 0

                else:
                    l = self._convert_to_basis(self.hecke_operator(self._basis_cache[j], m))
                    for i in range(self.prime()):
                        mat[i, j] = l[i]
                for i in range(self.prime(), n):
                    for u in range(self.prime()):
                        for v in range(self.prime()):
                            mat[i, j] = mat[i, j] + mat[i-u-1, j-v-1]*self.recurrence_matrix()[u, v]

        else:
            if n * self.prime() > self.prec():
                raise ValueError("n is too large")
            for j in range(self.prime(), n):
                l = self._convert_to_basis(self.hecke_operator(self._basis_cache[j], m))
                for i in range(n):
                    mat[i, j] = l[i]

        if side == 'right':
            return mat.transpose()

        return mat

    def slopes(self, n, use_recurrence=False):
        r"""
        Compute the slopes of the `U_p` operator acting on ``self``, using an `n\times n` matrix.

        EXAMPLES::

            sage: OverconvergentModularForms(5, 2, 1/3, base_ring=Qp(5), prec=100).slopes(5)
            [0, 2, 5, 6, 9]
            sage: OverconvergentModularForms(2, 1, 1/3, char=DirichletGroup(4,QQ).0).slopes(5)
            [0, 2, 4, 6, 8]
        """
        if self.base_ring() == QQ:
            slopelist = self.cps_u(n).truncate().newton_slopes(self.prime())
        elif isinstance(self.base_ring(), sage.rings.abc.pAdicField):
            slopelist = self.cps_u(n).truncate().newton_slopes()
        else:
            print("slopes are only defined for base field QQ or a p-adic field")
        return [-i for i in slopelist]

    def eigenfunctions(self, n, F=None, exact_arith=True):
        """
        Calculate approximations to eigenfunctions of ``self``.

        These are the eigenfunctions of ``self.hecke_matrix(p, n)``, which
        are approximations to the true eigenfunctions. Returns a list
        of :class:`OverconvergentModularFormElement` objects, in increasing
        order of slope.

        INPUT:

        - ``n`` -- integer; the size of the matrix to use

        - ``F`` -- either ``None`` or a field over which to calculate eigenvalues. If the
          field is ``None``, the current base ring is used. If the base ring is not
          a `p`-adic ring, an error will be raised.

        - ``exact_arith`` -- boolean (default: ``True``); if ``True``, use exact
          rational arithmetic to calculate the matrix of the `U` operator and its
          characteristic power series, even when the base ring is an inexact
          `p`-adic ring. This is typically slower, but more numerically stable.

        NOTE: Try using ``set_verbose(1, 'sage/modular/overconvergent')`` to
        get more feedback on what is going on in this algorithm. For even more
        feedback, use 2 instead of 1.

        EXAMPLES::

            sage: X = OverconvergentModularForms(2, 2, 1/6).eigenfunctions(8, Qp(2, 100))
            sage: X[1]
            2-adic overconvergent modular form of weight-character 2 with q-expansion
            (1 + O(2^74))*q
             + (2^4 + 2^5 + 2^9 + 2^10 + 2^12 + 2^13 + 2^15 + 2^17 + 2^19 + 2^20
                 + 2^21 + 2^23 + 2^28 + 2^30 + 2^31 + 2^32 + 2^34 + 2^36 + 2^37
                 + 2^39 + 2^40 + 2^43 + 2^44 + 2^45 + 2^47 + 2^48 + 2^52 + 2^53
                 + 2^54 + 2^55 + 2^56 + 2^58 + 2^59 + 2^60 + 2^61 + 2^67 + 2^68
                 + 2^70 + 2^71 + 2^72 + 2^74 + 2^76 + O(2^78))*q^2
             + (2^2 + 2^7 + 2^8 + 2^9 + 2^12 + 2^13 + 2^16 + 2^17 + 2^21 + 2^23
                 + 2^25 + 2^28 + 2^33 + 2^34 + 2^36 + 2^37 + 2^42 + 2^45 + 2^47
                 + 2^49 + 2^50 + 2^51 + 2^54 + 2^55 + 2^58 + 2^60 + 2^61 + 2^67
                 + 2^71 + 2^72 + O(2^76))*q^3
             + (2^8 + 2^11 + 2^14 + 2^19 + 2^21 + 2^22 + 2^24 + 2^25 + 2^26
                 + 2^27 + 2^28 + 2^29 + 2^32 + 2^33 + 2^35 + 2^36 + 2^44 + 2^45
                 + 2^46 + 2^47 + 2^49 + 2^50 + 2^53 + 2^54 + 2^55 + 2^56 + 2^57
                 + 2^60 + 2^63 + 2^66 + 2^67 + 2^69 + 2^74 + 2^76 + 2^79 + 2^80
                 + 2^81 + O(2^82))*q^4
             + (2 + 2^2 + 2^9 + 2^13 + 2^15 + 2^17 + 2^19 + 2^21 + 2^23 + 2^26
                 + 2^27 + 2^28 + 2^30 + 2^33 + 2^34 + 2^35 + 2^36 + 2^37 + 2^38
                 + 2^39 + 2^41 + 2^42 + 2^43 + 2^45 + 2^58 + 2^59 + 2^60 + 2^61
                 + 2^62 + 2^63 + 2^65 + 2^66 + 2^68 + 2^69 + 2^71 + 2^72 + O(2^75))*q^5
             + (2^6 + 2^7 + 2^15 + 2^16 + 2^21 + 2^24 + 2^25 + 2^28 + 2^29 + 2^33
                 + 2^34 + 2^37 + 2^44 + 2^45 + 2^48 + 2^50 + 2^51 + 2^54 + 2^55
                 + 2^57 + 2^58 + 2^59 + 2^60 + 2^64 + 2^69 + 2^71 + 2^73 + 2^75
                 + 2^78 + O(2^80))*q^6 + (2^3 + 2^8 + 2^9 + 2^10 + 2^11 + 2^12
                 + 2^14 + 2^15 + 2^17 + 2^19 + 2^20 + 2^21 + 2^23 + 2^25 + 2^26
                 + 2^34 + 2^37 + 2^38 + 2^39 + 2^40 + 2^41 + 2^45 + 2^47 + 2^49
                 + 2^51 + 2^53 + 2^54 + 2^55 + 2^57 + 2^58 + 2^59 + 2^60 + 2^61
                 + 2^66 + 2^69 + 2^70 + 2^71 + 2^74 + 2^76 + O(2^77))*q^7
             + O(q^8)
            sage: [x.slope() for x in X]
            [0, 4, 8, 14, 16, 18, 26, 30]
        """

        if F is None:
            F = self.base_ring()

        if F.is_exact():
            # raise TypeError("cannot calculate eigenfunctions over exact base fields")
            F = pAdicField(self.prime(), 100)  # noqa:F821

        m = self.hecke_matrix(self.prime(), n, use_recurrence=True, exact_arith=exact_arith)
        cp = m.charpoly()
        eigenvalues = cp.roots(F)
        eigenfunctions = []
        verbose("Expected %s eigenvalues, got %s" % (n, len(eigenvalues)))
        for (r, d) in eigenvalues:
            if d != 1:
                continue

            mr = m.__pari__() - r.__pari__()
            # Annoying thing: r isn't quite as precise as it claims to be
            # (bug reported to sage-support list)
            while F(mr.matdet()) != 0:
                verbose("p-adic solver returned wrong result in slope %s; refining" % r.valuation(), level=2)
                r = r - cp(r)/cp.derivative()(r)
                mr2 = m.__pari__() - r.__pari__()
                if mr2.matdet().valuation(self.prime()) > mr.matdet().valuation(self.prime()):
                    mr = mr2
                else:
                    mr = None
                    break

            if mr is None:
                verbose("Unable to calculate exact root in slope %s" % r.valuation())
                continue

            # now calculate the kernel using PARI

            v = mr.matker()

            if repr(v) == "[;]":
                verbose("PARI returned empty eigenspace in slope %s" % r.valuation())
                continue
                # Can't happen? Does PARI always return a
                # nonempty kernel for matrices that have det
                # indistinguishable from 0?

            if v.ncols() != 1:
                verbose("PARI returned non-simple eigenspace in slope %s" % r.valuation())
                continue

            gexp = self._gsr(0)
            for i in range(v.nrows()):
                gexp += self._gsr.gen()**i * F(v[i, 0])
            gexp = gexp + O(self._gsr.gen()**int(v.nrows()))

            if gexp[0] != 0:
                gexp = gexp/gexp[0]
            elif gexp[1] != 0:
                gexp = gexp/gexp[1]/self._const
            # This is slightly subtle. We want all eigenfunctions to have q-exps in Z_p.
            # Normalising the q-term to be 1 doesn't work for the Eisenstein series if
            # we're in the 0 component of weight-character space. But normalising the const term
            # to 1 works as *none of the small primes we deal with are irregular*! :-)
            else:
                raise ValueError("Constant and linear terms both zero!")
                # if this gets called something is very wrong.

            efunc = OverconvergentModularFormElement(self.base_extend(F), gexp=gexp)
            efunc._notify_eigen(r)
            assert efunc.is_integral()
            # This sometimes fails if n is too large -- last row of matrix fills
            # up with garbage. I don't know why. XXX FIX THIS XXX
            eigenfunctions.append((r.valuation(), efunc))

        eigenfunctions.sort()  # sort by slope
        return [f for _, f in eigenfunctions]

    def recurrence_matrix(self, use_smithline=True):
        r"""
        Return the recurrence matrix satisfied by the coefficients of `U`.

        This is a matrix `R =(r_{rs})_{r,s=1,\dots,p}` such that `u_{ij} =
        \sum_{r,s=1}^p r_{rs} u_{i-r, j-s}`.

        Uses an elegant construction which the author believes to be due
        to Smithline. See [Loe2007]_.

        EXAMPLES::

            sage: OverconvergentModularForms(2, 0, 0).recurrence_matrix()
            [  48    1]
            [4096    0]
            sage: OverconvergentModularForms(2, 0, 1/2).recurrence_matrix()
            [48 64]
            [64  0]
            sage: OverconvergentModularForms(3, 0, 0).recurrence_matrix()
            [   270     36      1]
            [ 26244    729      0]
            [531441      0      0]
            sage: OverconvergentModularForms(5, 0, 0).recurrence_matrix()
            [     1575      1300       315        30         1]
            [   162500     39375      3750       125         0]
            [  4921875    468750     15625         0         0]
            [ 58593750   1953125         0         0         0]
            [244140625         0         0         0         0]
            sage: OverconvergentModularForms(7, 0, 0).recurrence_matrix()
            [       4018        8624        5915        1904         322          28           1]
            [     422576      289835       93296       15778        1372          49           0]
            [   14201915     4571504      773122       67228        2401           0           0]
            [  224003696    37882978     3294172      117649           0           0           0]
            [ 1856265922   161414428     5764801           0           0           0           0]
            [ 7909306972   282475249           0           0           0           0           0]
            [13841287201           0           0           0           0           0           0]
            sage: OverconvergentModularForms(13, 0, 0).recurrence_matrix()
            [         15145         124852         354536 ...
        """

        if self._cached_recurrence_matrix is not None:
            return self._cached_recurrence_matrix

        MM = OverconvergentModularForms(self.prime(), 0, 0, base_ring=QQ)
        m = MM._discover_recurrence_matrix(use_smithline=True).base_extend(self.base_ring())

        r = diagonal_matrix([self._const**i for i in range(self.prime())])
        self._cached_recurrence_matrix = (r**(-1)) * m * r
        self._cached_recurrence_matrix.set_immutable()
        return self._cached_recurrence_matrix

    def _discover_recurrence_matrix(self, use_smithline=True):
        r"""
        Do the hard work of calculating the recurrence matrix, which is cached
        to avoid doing this every time.

        EXAMPLES::

            sage: o = OverconvergentModularForms(3,12,0)
            sage: o._discover_recurrence_matrix() == o.recurrence_matrix()
            True
        """

        (f_ring, f) = PolynomialRing(self.base_ring(), "f").objgen()

        if use_smithline:
            # Compute Smithline's polynomial H_p
            jq = self._qsr(j_invariant_qexp(1 + self.prime()).shift(1).power_series())  # noqa:F821

            # avoid dividing by q so as not to instantiate a Laurent series
            h = self._uniformiser.shift(-1) * jq
            fi = self._qsr(1)
            coeffs = []
            for i in range(self.prime() + 2):
                if not h.valuation() >= i:
                    raise ValueError("Something strange is happening here")

                coeffs.append(h[i] / fi[i])
                h = h - coeffs[-1] * fi
                fi = fi*self._uniformiser
            SmiH = f_ring(coeffs)
            assert SmiH.degree() == self.prime() + 1
            xyring = PolynomialRing(self.base_ring(), ["x", "y"], 2)
            x, y = xyring.gens()
            cc = self.prime() ** (-12/(self.prime() - 1))
            bigI = x*SmiH(y*cc) - y*cc*SmiH(x)
            smallI = xyring(bigI / (x - cc*y))
            r = matrix(ZZ, self.prime(), self.prime())
            for i in range(self.prime()):
                for j in range(self.prime()):
                    r[i, j] = -smallI[i+1, j+1]
            return r
        else:
            # compute from U(f^j) for small j via Newton's identities
            # to be implemented when I can remember Newton's identities!
            raise NotImplementedError

    def cps_u(self, n, use_recurrence=False):
        r"""
        Compute the characteristic power series of `U_p` acting on ``self``, using
        an `n\times n` matrix.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 16, 1/2, base_ring=Qp(3)).cps_u(4)
            1 + O(3^20)
             + (2 + 2*3 + 2*3^2 + 2*3^4 + 3^5 + 3^6 + 3^7
                 + 3^11 + 3^12 + 2*3^14 + 3^16 + 3^18 + O(3^19))*T
             + (2*3^3 + 3^5 + 3^6 + 3^7 + 2*3^8 + 2*3^9 + 2*3^10
                 + 3^11 + 3^12 + 2*3^13 + 2*3^16 + 2*3^18 + O(3^19))*T^2
             + (2*3^15 + 2*3^16 + 2*3^19 + 2*3^20 + 2*3^21 + O(3^22))*T^3
             + (3^17 + 2*3^18 + 3^19 + 3^20 + 3^22 + 2*3^23 + 2*3^25 + 3^26 + O(3^27))*T^4
            sage: OverconvergentModularForms(3, 16, 1/2, base_ring=Qp(3), prec=30).cps_u(10)
            1 + O(3^20)
             + (2 + 2*3 + 2*3^2 + 2*3^4 + 3^5 + 3^6 + 3^7 + 2*3^15 + O(3^16))*T
             + (2*3^3 + 3^5 + 3^6 + 3^7 + 2*3^8 + 2*3^9 + 2*3^10
                 + 2*3^11 + 2*3^12 + 2*3^13 + 3^14 + 3^15 + O(3^16))*T^2
             + (3^14 + 2*3^15 + 2*3^16 + 3^17 + 3^18 + O(3^19))*T^3
             + (3^17 + 2*3^18 + 3^19 + 3^20 + 3^21 + O(3^24))*T^4
             + (3^29 + 2*3^32 + O(3^33))*T^5
             + (2*3^44 + O(3^45))*T^6
             + (2*3^59 + O(3^60))*T^7
             + (2*3^78 + O(3^79))*T^8

        .. NOTE::

            Uses the Hessenberg form of the Hecke matrix to compute
            the characteristic polynomial.  Because of the use of
            relative precision here this tends to give better
            precision in the `p`-adic coefficients.
        """
        m = self.hecke_matrix(self.prime(), n, use_recurrence)
        A = PowerSeriesRing(self.base_ring(), 'T')

        # From a conversation with David Loeffler, apparently self.base_ring()
        # is either the field of rational numbers or some p-adic field.  In the
        # first case we want to use the linbox algorithm, and in the second
        # case the Hessenberg form algorithm.
        #
        if self.base_ring().is_exact():
            g = A(m.charpoly('T').reverse())
        else:
            g = A(m.charpoly('T', "hessenberg").reverse())

        return g


class OverconvergentModularFormElement(ModuleElement):
    r"""
    A class representing an element of a space of overconvergent modular forms.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<w> = Qp(5).extension(x^7 - 5)
        sage: s = OverconvergentModularForms(5, 6, 1/21, base_ring=K).0
        sage: s == loads(dumps(s))
        True
    """

    def __init__(self, parent, gexp=None, qexp=None):
        r"""
        Create an element of this space.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 2, 1/6,prec=5).an_element() # indirect doctest
            3-adic overconvergent modular form of weight-character 2
             with q-expansion 3*q + 72*q^2 + 810*q^3 + 6096*q^4 + O(q^5)
        """
        ModuleElement.__init__(self, parent)

        self._p = self.parent().prime()
        # self.weight = self.parent().weight
        if (gexp is None and qexp is None) or (gexp is not None and qexp is not None):
            raise ValueError("Must supply exactly one of a q-expansion and a g-expansion")
        if gexp is not None:
            self._gexp = gexp.add_bigoh(self.parent().prec())
            self._qexp = sum([self.parent()._basis_cache[i] * gexp[i] for i in range(min(gexp.prec(), self.parent().prec()))])
            self._qexp = self._qexp.add_bigoh(self._gexp.prec())
        else:  # qexp is not None
            self._qexp = qexp.add_bigoh(self.parent().prec())
            self._gexp = self.parent()._convert_to_basis(self._qexp)

        self._is_eigen = False
        self._eigenvalue = None
        self._slope = None

    def _add_(self, other):
        r"""
        Add ``self`` to ``other``.

        Here ``other`` has the same parent as ``self``.

        EXAMPLES::

            sage: M = OverconvergentModularForms(2, 12, 1/6)
            sage: f = M.0
            sage: f + f # indirect doctest
            2-adic overconvergent modular form of weight-character 12 with q-expansion 2 - 131040/1414477*q ...
        """
        return OverconvergentModularFormElement(self.parent(), gexp=self.gexp() + other.gexp())

    def _lmul_(self, x):
        r"""
        Left multiplication by other.

        EXAMPLES::

            sage: M = OverconvergentModularForms(2, 12, 1/6)
            sage: f = M.0
            sage: 2*f # indirect doctest
            2-adic overconvergent modular form of weight-character 12 with q-expansion 2 - 131040/1414477*q ...
        """
        return OverconvergentModularFormElement(self.parent(), gexp=x * self.gexp())

    def _rmul_(self, x):
        r"""
        Right multiplication by other.

        EXAMPLES::

            sage: M = OverconvergentModularForms(2, 12, 1/6)
            sage: f = M.0
            sage: f * 3 # indirect doctest
            2-adic overconvergent modular form of weight-character 12 with q-expansion 3 - 196560/1414477*q ...
        """
        return OverconvergentModularFormElement(self.parent(), gexp=x * self.gexp())

    def prec(self):
        r"""
        Return the series expansion precision of this overconvergent modular form.

        This is not the same as the `p`-adic precision of the coefficients.

        EXAMPLES::

            sage: OverconvergentModularForms(5, 6, 1/3, prec=15).gen(1).prec()
            15
        """
        return self.gexp().prec()

    def is_eigenform(self):
        r"""
        Return ``True`` if this is an eigenform.

        At present this returns ``False`` unless this element was explicitly
        flagged as an eigenform, using the method :meth:`_notify_eigen`.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2)
            sage: f = M.eigenfunctions(3)[1]
            sage: f.is_eigenform()
            True
            sage: M.gen(4).is_eigenform()
            False
        """
        return self._is_eigen

    def slope(self):
        r"""
        Return the slope of this eigenform.

        This is the valuation of its `U_p`-eigenvalue.

        Raises an error unless this element was explicitly
        flagged as an eigenform, using the method :meth:`_notify_eigen`.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2)
            sage: f = M.eigenfunctions(3)[1]
            sage: f.slope()
            2
            sage: M.gen(4).slope()
            Traceback (most recent call last):
            ...
            TypeError: slope only defined for eigenfunctions
        """
        if not self.is_eigenform():
            raise TypeError("slope only defined for eigenfunctions")
        return self._slope

    def eigenvalue(self):
        r"""
        Return the `U_p`-eigenvalue of this eigenform.

        This raises an error unless this element was explicitly flagged
        as an eigenform, using the method :meth:`_notify_eigen`.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2)
            sage: f = M.eigenfunctions(3)[1]
            sage: f.eigenvalue()
            3^2 + 3^4 + 2*3^6 + 3^7 + 3^8 + 2*3^9 + 2*3^10 + 3^12 + 3^16 + 2*3^17
             + 3^18 + 3^20 + 2*3^21 + 3^22 + 2*3^23 + 3^25 + 3^26 + 2*3^27 + 2*3^29
             + 3^30 + 3^31 + 3^32 + 3^33 + 3^34 + 3^36 + 3^40 + 2*3^41 + 3^43 + 3^44
             + 3^45 + 3^46 + 3^48 + 3^49 + 3^50 + 2*3^51 + 3^52 + 3^54 + 2*3^57
             + 2*3^59 + 3^60 + 3^61 + 2*3^63 + 2*3^66 + 2*3^67 + 3^69 + 2*3^72
             + 3^74 + 2*3^75 + 3^76 + 2*3^77 + 2*3^78 + 2*3^80 + 3^81 + 2*3^82
             + 3^84 + 2*3^85 + 2*3^86 + 3^87 + 3^88 + 2*3^89 + 2*3^91 + 3^93 + 3^94
             + 3^95 + 3^96 + 3^98 + 2*3^99 + O(3^100)
            sage: M.gen(4).eigenvalue()
            Traceback (most recent call last):
            ...
            TypeError: eigenvalue only defined for eigenfunctions
        """
        if not self.is_eigenform():
            raise TypeError("eigenvalue only defined for eigenfunctions")
        return self._eigenvalue

    def q_expansion(self, prec=None):
        r"""
        Return the `q`-expansion of ``self``, to as high precision as it is known.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 4, 1/2).gen(0).q_expansion()
            1 - 120/13*q - 1080/13*q^2 - 120/13*q^3 - 8760/13*q^4 - 15120/13*q^5
             - 1080/13*q^6 - 41280/13*q^7 - 5400*q^8 - 120/13*q^9 - 136080/13*q^10
             - 159840/13*q^11 - 8760/13*q^12 - 263760/13*q^13 - 371520/13*q^14
             - 15120/13*q^15 - 561720/13*q^16 - 45360*q^17 - 1080/13*q^18
             - 823200/13*q^19 + O(q^20)
        """
        if prec is None:
            return self._qexp
        elif prec > self.prec():
            raise ValueError
        else:
            return self._qexp.add_bigoh(prec)

    def gexp(self):
        r"""
        Return the formal power series in `g` corresponding to ``self``.

        If this overconvergent modular form is `E_k^\ast \times F(g)`
        where `g` is the appropriately normalised parameter of `X_0(p)`,
        the result is `F`.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2)
            sage: f = M.eigenfunctions(3)[1]
            sage: f.gexp()
            (3^-3 + O(3^95))*g
            + (3^-1 + 1 + 2*3 + 3^2 + 2*3^3 + 3^5 + 3^7 + 3^10 + 3^11 + 3^14 + 3^15
                + 3^16 + 2*3^19 + 3^21 + 3^22 + 2*3^23 + 2*3^24 + 3^26 + 2*3^27
                + 3^29 + 3^31 + 3^34 + 2*3^35 + 2*3^36 + 3^38 + 2*3^39 + 3^41 + 2*3^42
                + 2*3^43 + 2*3^44 + 2*3^46 + 2*3^47 + 3^48 + 2*3^49 + 2*3^50 + 3^51
                + 2*3^54 + 2*3^55 + 2*3^56 + 3^57 + 2*3^58 + 2*3^59 + 2*3^60 + 3^61
                + 3^62 + 3^63 + 3^64 + 2*3^65 + 3^67 + 3^68 + 2*3^69 + 3^70 + 2*3^71
                + 2*3^74 + 3^76 + 2*3^77 + 3^78 + 2*3^79 + 2*3^80 + 3^84 + 2*3^85
                + 2*3^86 + 3^88 + 2*3^89 + 3^91 + 3^92 + 2*3^94 + 3^95 + O(3^97))*g^2
            + O(g^3)
        """
        return self._gexp

    def coordinates(self, prec=None):
        r"""
        Return the coordinates of this modular form in terms of the basis of this space.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2, prec=15)
            sage: f = (M.0 + M.3); f.coordinates()
            [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: f.coordinates(6)
            [1, 0, 0, 1, 0, 0]
            sage: OverconvergentModularForms(3, 0, 1/6)(f).coordinates(6)
            [1, 0, 0, 729, 0, 0]
            sage: f.coordinates(100)
            Traceback (most recent call last):
            ...
            ValueError: Precision too large for space
        """
        if prec > self.prec():
            raise ValueError("Precision too large for space")
        if prec is None:
            prec = self.prec()
        return self._gexp.padded_list(prec)

    def prime(self):
        r"""
        If this is a `p`-adic modular form, return `p`.

        EXAMPLES::

            sage: OverconvergentModularForms(2, 0, 1/2).an_element().prime()
            2
        """
        return self._p

    def _notify_eigen(self, eigenvalue):
        """
        Flag this element as an eigenform.

        It then remembers some extra data.

        EXAMPLES::

            sage: OverconvergentModularForms(3, 16, 1/3).eigenfunctions(4)  # indirect doctest
            [...]
        """
        self._is_eigen = True
        self._eigenvalue = eigenvalue
        self._slope = eigenvalue.normalized_valuation()

    def is_integral(self):
        r"""
        Test whether this element has `q`-expansion coefficients that are `p`-adically integral.

        This should always be the case with eigenfunctions, but sometimes
        if `n` is very large this breaks down for unknown reasons!

        EXAMPLES::

            sage: M = OverconvergentModularForms(2, 0, 1/3)
            sage: q = QQ[['q']].gen()
            sage: M(q - 17*q^2 + O(q^3)).is_integral()
            True
            sage: M(q - q^2/2 + 6*q^7  + O(q^9)).is_integral()
            False
        """
        for co in self.q_expansion().list():
            if (co * (1 + O(self.prime()))).valuation() < 0:  # have to force it into ZZ_p
                return False
        return True

    def _repr_(self) -> str:
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: o = OverconvergentModularForms(3, 0, 1/2)
            sage: o([1, 0, 1, 3])._repr_()
            '3-adic overconvergent modular form of weight-character 0 with q-expansion 1 + 729*q^2 + 76545*q^3 + O(q^4)'
        """
        return "%s-adic overconvergent modular form of weight-character %s with q-expansion %s" % (self.prime(), self.weight(), self.q_expansion())

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: o = OverconvergentModularForms(3, 0, 1/2)
            sage: o([1, 1, 1, 0, 0, 0, 0]) == o([2, 1, 0])
            False
            sage: o([1, 1, 1, 0, 0, 0, 0]) == o([1,1])
            True
        """
        return richcmp(self.gexp(), other.gexp(), op)

    def r_ord(self, r):
        r"""
        The `p`-adic valuation of the norm of ``self`` on the `r`-overconvergent region.

        EXAMPLES::

            sage: o = OverconvergentModularForms(3, 0, 1/2)
            sage: t = o([1, 1, 1/3])
            sage: t.r_ord(1/2)
            1
            sage: t.r_ord(2/3)
            3
        """
        ord = -Infinity
        p = self.prime()
        s = self.parent().radius()

        F = self.parent().base_ring()
        if not isinstance(F, sage.rings.abc.pAdicField):
            F = pAdicField(p)  # noqa:F821

        for i in range(self.prec()):
            ord = max(ord, 12/ZZ(p - 1)*i*(r - s) - F(self.gexp()[i]).normalized_valuation())

        return ord

    def valuation(self):
        r"""
        Return the `p`-adic valuation of this form.

        This is the minimum of the `p`-adic valuations of its coordinates.

        EXAMPLES::

            sage: M = OverconvergentModularForms(3, 0, 1/2)
            sage: (M.7).valuation()
            0
            sage: (3^18 * (M.2)).valuation()
            18
        """
        if isinstance(self.parent().base_ring(), sage.rings.abc.pAdicField):
            v = lambda u: u.normalized_valuation()
        else:
            v = lambda u: u.valuation(self.parent().prime())
        return min([v(x) for x in self.gexp().list()])

    def governing_term(self, r):
        r"""
        The degree of the series term with largest norm on the `r`-overconvergent region.

        EXAMPLES::

            sage: o = OverconvergentModularForms(3, 0, 1/2)
            sage: f = o.eigenfunctions(10)[1]
            sage: f.governing_term(1/2)
            1
        """
        p = self.prime()
        F = self.parent().base_ring()
        if not isinstance(F, sage.rings.abc.pAdicField):
            F = pAdicField(p)  # noqa:F821
        s = self.parent().radius()
        p = self.prime()

        for i in range(self.gexp().prec()):
            if 12/ZZ(p - 1)*i*(r - s) - F(self.gexp()[i]).normalized_valuation() == self.r_ord(r):
                return i
        raise RuntimeError("Can't get here")

    def valuation_plot(self, rmax=None):
        r"""
        Draw a graph depicting the growth of the norm of this overconvergent
        modular form as it approaches the boundary of the overconvergent
        region.

        EXAMPLES::

            sage: o = OverconvergentModularForms(3, 0, 1/2)
            sage: f = o.eigenfunctions(4)[1]
            sage: f.valuation_plot()                                                    # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.plot import plot

        if rmax is None:
            rmax = ZZ(self.prime()) / ZZ(1 + self.prime())
        return plot(self.r_ord, (0, rmax))

    def weight(self):
        r"""
        Return the weight of this overconvergent modular form.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: R = Qp(13).extension(x^2 - 13, names='a')
            sage: M = OverconvergentModularForms(13, 10, 1/2, base_ring=R)
            sage: M.gen(0).weight()
            10
        """
        return self.parent().weight()

    def additive_order(self):
        r"""
        Return the additive order of this element.

        This implements a required method for all
        elements deriving from :class:`sage.modules.ModuleElement`.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: R = Qp(13).extension(x^2 - 13, names='a')
            sage: M = OverconvergentModularForms(13, 10, 1/2, base_ring=R)
            sage: M.gen(0).additive_order()
            +Infinity
            sage: M(0).additive_order()
            1
        """
        from sage.rings.infinity import Infinity
        if self.is_zero():
            return ZZ(1)
        else:
            return Infinity

    def base_extend(self, R):
        r"""
        Return a copy of ``self`` but with coefficients in the given ring.

        EXAMPLES::

            sage: M = OverconvergentModularForms(7, 10, 1/2, prec=5)
            sage: f = M.1
            sage: f.base_extend(Qp(7, 4))
            7-adic overconvergent modular form of weight-character 10 with q-expansion
             (7 + O(7^5))*q + (6*7 + 4*7^2 + 7^3 + 6*7^4 + O(7^5))*q^2
             + (5*7 + 5*7^2 + 7^4 + O(7^5))*q^3 + (7^2 + 4*7^3 + 3*7^4 + 2*7^5
             + O(7^6))*q^4 + O(q^5)
        """
        S = self.parent().base_extend(R)
        return S(self)

    def __pari__(self):
        r"""
        Return the Pari object corresponding to ``self``.

        This is just the `q`-expansion of ``self`` as a formal power series.

        EXAMPLES::

            sage: f = OverconvergentModularForms(3, 0, 1/2).1
            sage: pari(f) # indirect doctest
            27*q + 324*q^2 + 2430*q^3 + 13716*q^4 + 64557*q^5 + 265356*q^6 + 983556*q^7 + 3353076*q^8 + 10670373*q^9 + 32031288*q^10 + 91455804*q^11 + 249948828*q^12 + 657261999*q^13 + 1669898592*q^14 + 4113612864*q^15 + 9853898292*q^16 + 23010586596*q^17 + 52494114852*q^18 + 117209543940*q^19 + O(q^20)
            sage: pari(f.base_extend(Qp(3))) # indirect doctest
            (3^3 + O(3^23))*q + (3^4 + 3^5 + O(3^24))*q^2 + (3^5 + 3^7 + O(3^25))*q^3 + (3^3 + 3^4 + 2*3^5 + 2*3^8 + O(3^23))*q^4 + (2*3^4 + 3^5 + 3^6 + 2*3^7 + 3^10 + O(3^24))*q^5 + (3^6 + 3^7 + 3^8 + 3^9 + 3^10 + 3^11 + O(3^26))*q^6 + (2*3^3 + 3^4 + 2*3^6 + 2*3^7 + 2*3^8 + 3^9 + 3^10 + 2*3^11 + 3^12 + O(3^23))*q^7 + (2*3^4 + 3^5 + 3^8 + 2*3^9 + 2*3^10 + 2*3^13 + O(3^24))*q^8 + (3^7 + 2*3^9 + 2*3^12 + 2*3^14 + O(3^27))*q^9 + (2*3^5 + 3^8 + 3^9 + 2*3^10 + 2*3^13 + 2*3^15 + O(3^25))*q^10 + (3^4 + 2*3^5 + 2*3^6 + 3^8 + 2*3^9 + 3^12 + 3^14 + 2*3^16 + O(3^24))*q^11 + (3^5 + 3^6 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^12 + 3^14 + 2*3^15 + 2*3^16 + 3^17 + O(3^25))*q^12 + (2*3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 2*3^11 + 3^13 + 2*3^14 + 2*3^17 + 3^18 + O(3^23))*q^13 + (2*3^4 + 2*3^6 + 2*3^7 + 3^8 + 2*3^9 + 3^10 + 3^12 + 3^14 + 2*3^15 + 2*3^16 + 3^18 + 3^19 + O(3^24))*q^14 + (2*3^6 + 3^7 + 3^9 + 3^10 + 3^11 + 2*3^14 + 3^15 + 2*3^16 + 3^17 + 3^18 + 3^20 + O(3^26))*q^15 + (3^3 + 2*3^4 + 2*3^7 + 2*3^8 + 3^9 + 3^10 + 2*3^11 + 3^12 + 2*3^14 + 2*3^15 + 3^17 + 3^18 + 2*3^19 + 2*3^20 + O(3^23))*q^16 + (2*3^5 + 2*3^7 + 2*3^8 + 3^10 + 3^11 + 2*3^12 + 2*3^13 + 3^14 + 3^15 + 3^17 + 2*3^18 + 3^19 + 2*3^21 + O(3^25))*q^17 + (3^8 + 3^9 + 2*3^10 + 2*3^11 + 3^12 + 3^14 + 3^15 + 3^16 + 3^17 + 2*3^21 + 3^22 + O(3^28))*q^18 + (2*3^3 + 3^5 + 2*3^6 + 2*3^8 + 2*3^9 + 3^11 + 2*3^12 + 3^13 + 3^14 + 2*3^15 + 3^16 + 3^17 + 2*3^18 + 3^19 + 2*3^21 + O(3^23))*q^19 + O(q^20)
        """
        return self.q_expansion().__pari__()


OverconvergentModularFormsSpace.Element = OverconvergentModularFormElement
