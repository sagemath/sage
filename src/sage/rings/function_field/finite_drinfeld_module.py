r"""
This module provides classes for finite Drinfeld modules
(`FiniteDrinfeldModule`) and their module action on the algebraic
closure of `\Fq` (`FiniteDrinfeldModuleAction`).

AUTHORS:

- Antoine Leudière (2022-04): initial version

Let `\tau` be the `\Fq`-linear Frobenius endomorphism of `\Fqbar`
defined by `x \mapsto x^q`. Let `L\{\tau\}` be the ring of Ore
polynomials in `\tau` with coefficients in L. Fix an element `\omega` in
`L` (global parameter). A finite Drinfeld module is an `\Fq`-algebra
morphism `\phi: \Fq[X] \to L\{\tau\]` such that:
    - the constant coefficient of `\phi(X)` is `\omega`,
    - there exists at least one `a \in \Fq[X]` such that `\phi(a)` has a
      non zero `\tau`-degree.

As an `\Fq[X]`-algebra morphism, a finite Drinfeld module is only
determined by the image of `X`.

Crucially, the Drinfeld module `\phi` gives rise to the `\Fq[X]`-module
law on `\Fqbar` defined by `(a, x) = \phi(a)(x)`.
"""

#*****************************************************************************
#       Copyright (C) 2022 Antoine Leudière <antoine.leudiere@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.action import Action
from sage.categories.homset import Hom
from sage.matrix.constructor import Matrix
from sage.misc.latex import latex
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.rings.morphism import RingHomomorphism_im_gens
from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_finite_field
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class FiniteDrinfeldModule(RingHomomorphism_im_gens):
    r"""
    Class for finite Drinfeld modules.

    INPUT:

    - ``polring`` -- the base polynomial ring
    - ``gen`` -- the generator of the Drinfeld module, i.e. the image of `X` in
      the Ore polynomial ring
    - ``characteristic`` -- the Fq[X]-characteristic of the Drinfeld
      module, i.e. the minimal polynomial in `polring` of the constant term of
      the generator

    EXAMPLES:

    .. RUBRIC:: Instantiation

    We must first create the base objects::

        sage: Fq = GF(3^2)
        sage: z2 = Fq.gen()
        sage: FqX.<X> = Fq[]
        sage: p = X^3 + (z2 + 2)*X^2 + (6*z2 + 1)*X + 3*z2 + 5
        sage: L = Fq.extension(6)
        sage: frob = L.frobenius_endomorphism(2)
        sage: Ltau.<t> = OrePolynomialRing(L, frob)
        sage: omega = p.roots(L, multiplicities=False)[0]
        sage: phi_X = omega + t + t^2

    Notice that we have freedom on choosing the polynomial `p`, but not
    `omega`. It is generally more useful this way. Then we instantiate the
    Drinfeld module::

        sage: phi = FiniteDrinfeldModule(FqX, phi_X, p)
        sage: phi
        Finite Drinfeld module:
          Polring:        Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          Ore polring:    Ore Polynomial Ring in t over Finite Field in z12 of size 3^12 twisted by z12 |--> z12^(3^2)
          Generator:      t^2 + t + z12^11 + z12^10 + z12^9 + 2*z12^5 + 2*z12^4 + z12^3 + 2*z12
          Characteristic: X^3 + (z2 + 2)*X^2 + X + 2

    .. RUBRIC:: Getters

    With getters, we can retrieve many basic objects associated to a Drinfeld
    module.

    First, we can retrieve the polynomial ring, the Ore polynomial ring, and
    the generator. Note that the class inherits `RingHomomorphism_im_gens`, so
    that `domain`, `codomain` and `im_gens` are available::

        sage: phi.polring()
        Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
        sage: phi.domain() is phi.polring()
        True
        sage: phi.ore_polring()
        Ore Polynomial Ring in t over Finite Field in z12 of size 3^12 twisted by z12 |--> z12^(3^2)
        sage: phi.codomain() is phi.ore_polring()
        True
        sage: phi.gen()
        t^2 + t + z12^11 + z12^10 + z12^9 + 2*z12^5 + 2*z12^4 + z12^3 + 2*z12
        sage: phi.im_gens()[0] is phi.gen()
        True

    We can retrieve `omega`, the constant term of the generator, and ensure
    that it is a root of `p`, the characteristic::

        sage: phi.constant_term()
        z12^11 + z12^10 + z12^9 + 2*z12^5 + 2*z12^4 + z12^3 + 2*z12
        sage: phi.constant_term() == omega
        True
        sage: phi.characteristic()
        X^3 + (z2 + 2)*X^2 + X + 2
        sage: phi.characteristic() == p
        True
        sage: phi.characteristic()(phi.constant_term())
        0

    We can retrieve the rank and the height (note that the height is always one
    here)::

        sage: phi.rank()
        2
        sage: phi.height()
        1

    And finally we can retrieve some rank-two specifics::

        sage: phi.j()  # j-invariant
        1
        sage: phi.g()  # Standard notation
        1
        sage: phi.delta()  # Standard notation
        1
        sage: phi(X) == phi.constant_term() + phi.g()*t + phi.delta()*t^2
        True

    .. RUBRIC:: Evaluation of the Drinfeld module

    By definition, a Drinfeld module is a ring morphism from an polynomial ring
    to an Ore polynomial ring. We can compute the images under this morphism
    using the standard `phi(...)` notation::

        sage: phi(X)
        t^2 + t + z12^11 + z12^10 + z12^9 + 2*z12^5 + 2*z12^4 + z12^3 + 2*z12
        sage: phi(X) == phi.gen()
        True
        sage: phi(1)
        1
        sage: phi(z2) == z2
        True
        sage: phi(X^2 + 1)
        t^4 + 2*t^3 + (2*z12^11 + 2*z12^10 + z12^9 + z12^6 + z12^5 + 2*z12^4 + z12^3 + 2*z12^2 + z12 + 2)*t^2 + (2*z12^8 + z12^7 + 2*z12^6 + z12^5 + z12^4 + z12 + 1)*t + 2*z12^11 + 2*z12^10 + z12^8 + z12^7 + 2*z12^6 + 2*z12^5 + z12^4 + 2*z12

    .. RUBRIC:: The module law induced by a Drinfeld module

    The most important feature of Drinfeld modules is that they endow any
    subextension of `Fqbar` with an `Fq[X]`-module law. This action is
    represented by the class `FiniteDrinfeldModuleAction`, which inherits
    `Action`. For the sake of simplicity, `phi` will only act on the base field
    (`L`) of its Ore polynomial ring. If you want to act on a bigger field, you
    can define a new Drinfeld module using the method `change_ring`.

        sage: action = phi._get_action_()
        sage: action
        Action on Finite Field in z12 of size 3^12 induced by Finite Drinfeld module:
          Polring:        Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          Ore polring:    Ore Polynomial Ring in t over Finite Field in z12 of size 3^12 twisted by z12 |--> z12^(3^2)
          Generator:      t^2 + t + z12^11 + z12^10 + z12^9 + 2*z12^5 + 2*z12^4 + z12^3 + 2*z12
          Characteristic: X^3 + (z2 + 2)*X^2 + X + 2
    
    The calculation of the action is simple. Careful, the evaluation of Ore
    polynomial is, at the moment, experimental::

        sage: x = L.gen() + 1
        sage: g = X^3 + X + 5
        sage: action(g, x)
        ...
        z12^11 + z12^10 + 2*z12^9 + z12^7 + z12^6 + z12^4 + 2*z12^2 + z12 + 1

    To change ring, use::

        sage: M = L.extension(5)
        sage: psi = phi.change_ring(M)
        sage: psi
        Finite Drinfeld module:
          Polring:        Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          Ore polring:    Ore Polynomial Ring in t over Finite Field in z60 of size 3^60 twisted by z60 |--> z60^(3^2)
          Generator:      t^2 + t + 2*z60^59 + z60^57 + 2*z60^56 + 2*z60^54 + 2*z60^53 + 2*z60^52 + z60^51 + z60^50 + 2*z60^48 + z60^47 + z60^46 + 2*z60^45 + 2*z60^44 + 2*z60^42 + 2*z60^41 + z60^40 + z60^39 + z60^38 + z60^37 + 2*z60^34 + z60^33 + z60^31 + 2*z60^29 + z60^27 + z60^26 + z60^25 + 2*z60^24 + z60^22 + 2*z60^21 + z60^19 + 2*z60^17 + 2*z60^16 + 2*z60^15 + z60^14 + z60^12 + z60^11 + 2*z60^10 + z60^8 + z60^7 + 2*z60^6 + 2*z60^5 + 2*z60 + 1
          Characteristic: X^3 + (z2 + 2)*X^2 + X + 2

    .. RUBRIC:: Morphisms and isogenies

    Being given an Ore polynomial `m`, we can decide if `m` is a morphism or
    isogeny of Drinfeld module whose domain is `phi`::

        sage: m = phi(X)
        sage: phi.is_morphism(m)
        True
        sage: phi.is_isogeny(m)
        True
        sage: phi.is_endomorphism(m)
        True
        sage: phi.is_automorphism(m)
        False
        sage: m = 0
        sage: phi.is_endomorphism(m)
        True
        sage: phi.is_automorphism(m)
        False
        sage: phi.is_isogeny(m)
        False
        sage: m = t^6
        sage: phi.is_endomorphism(m)
        True
        sage: phi.is_automorphism(m)
        False

    We implemented the Vélu formula for Drinfeld modules, in the sense that
    given `m`, we can compute (if it exists) the unique Drinfeld module `psi`
    such that `m` is an isogeny from `phi` to `psi`::

        sage: m = phi(X^2 + 1)
        sage: phi.velu(m)
        Finite Drinfeld module:
          Polring:        Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          Ore polring:    Ore Polynomial Ring in t over Finite Field in z12 of size 3^12 twisted by z12 |--> z12^(3^2)
          Generator:      t^2 + t + z12^11 + z12^10 + z12^9 + 2*z12^5 + 2*z12^4 + z12^3 + 2*z12
          Characteristic: X^3 + (z2 + 2)*X^2 + X + 2
        sage: phi.velu(m) == phi
        True
        sage: z12 = L.gen()
        sage: m = (z12^7 + z12^6 + 2*z12^4)*t^3
        sage: phi.is_isogeny(m)
        True
        sage: phi.is_endomorphism(m)
        False
        sage: phi.velu(m)
        Finite Drinfeld module:
          Polring:        Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          Ore polring:    Ore Polynomial Ring in t over Finite Field in z12 of size 3^12 twisted by z12 |--> z12^(3^2)
          Generator:      (2*z12^10 + z12^9 + 2*z12^7 + 2*z12^6 + z12^3 + 2*z12^2 + 2)*t^2 + (2*z12^9 + z12^7 + 2*z12^5 + z12 + 1)*t + z12^11 + z12^10 + z12^9 + 2*z12^5 + 2*z12^4 + z12^3 + 2*z12
          Characteristic: X^3 + (z2 + 2)*X^2 + X + 2
        sage: phi.velu(m) == phi
        False

    .. RUBRIC:: Complex multiplication of Drinfeld modules

    There are various methods to manipulate the complex multiplication theory
    of Drinfeld modules. We can compute the Frobenius norm, Frobenius trace and
    characteristic polynomial::

        sage: phi.frobenius_norm()
        X^6 + (2*z2 + 1)*X^5 + (2*z2 + 1)*X^4 + (2*z2 + 2)*X^3 + z2*X^2 + X + 1
        sage: phi.frobenius_trace()
        2*X^3 + (2*z2 + 1)*X^2 + 2*X + z2 + 2
        sage: phi.characteristic_polynomial()
        T^2 + (2*X^3 + (2*z2 + 1)*X^2 + 2*X + z2 + 2)*T + X^6 + (2*z2 + 1)*X^5 + (2*z2 + 1)*X^4 + (2*z2 + 2)*X^3 + z2*X^2 + X + 1

    With those methods, it is easy to decide if a Drinfeld module is ordinary
    or supersingular::

        sage: phi.is_ordinary()
        True
        sage: phi.is_supersingular()
        False

    .. NOTE::

    The general definition of a Drinfeld module is out of the scope of this
    implementation.

    ::

    You can see all available methods of `RingHomomorphism_im_gens` with
    `dir(sage.rings.morphism.RingHomomorphism_im_gens)`. Same for `Action`.

    .. SEEALSO::
        :mod:`sage.categories.action.Action`
        :mod:`sage.rings.polynomial.ore_polynomial_element`
        :mod:`sage.rings.polynomial.ore_polynomial_ring`
    """
    
    def __init__(self, polring, gen, characteristic):
        # Verifications
        if not isinstance(polring, PolynomialRing_dense_finite_field):
            raise TypeError('First argument must be a polynomial ring')
        if not characteristic in polring:
            raise TypeError('The characteristic must be in the polynomial ' \
                    'ring')
        if not characteristic.is_prime():
            raise ValueError('The characteristic must be irreducible')
        if not isinstance(gen, OrePolynomial):
            raise TypeError('The generator must be an Ore polynomial')
        # We define those for convenience before continuing
        FqX = polring
        Ltau = gen.parent()
        Fq = FqX.base_ring()
        L = Ltau.base_ring()
        _check_base_fields(Fq, L)
        if not Ltau.twisting_derivation() is None:
            raise ValueError('The Ore polynomial ring should have no ' \
                    'derivation')
        if Ltau.twisting_morphism().power() != Fq.degree():
            raise ValueError('The twisting morphism of the Ore polynomial ' \
                    'ring must be the Frobenius endomorphism of the base ' \
                    'field of the polynomial ring')
        if characteristic(gen[0]) != 0:
            raise ValueError('The constant coefficient of the generator ' \
                    'must be a root of the characteristic')
        if gen.is_constant():
            raise ValueError('The generator must not be constant')
        # Work
        self.__characteristic = characteristic
        super().__init__(Hom(FqX, Ltau), gen)

    #################
    # Private utils #
    #################

    def _Fq(self):
        return self.polring().base_ring()

    def _L(self):
        return self.ore_polring().base_ring()

    ###########
    # Methods #
    ###########

    def change_ring(self, R):
        # VERIFICATIONS
        if not R.is_field() and R.is_finite():
            raise TypeError('Argument must be a finite field')
        if not self.ore_polring().base_ring().is_subring(R):
            raise ValueError('The new field must be a finite field ' \
                    'extension of the base field of the Ore polynomial ring.')
        _check_base_fields(self.polring().base_ring(), R)
        # ACTUAL WORK
        new_frobenius = R.frobenius_endomorphism(self.frobenius().power())
        new_ore_polring = OrePolynomialRing(R, new_frobenius,
                names=self.ore_polring().variable_names())
        return FiniteDrinfeldModule(self.polring(),
                new_ore_polring(self.gen()), self.characteristic())

    def characteristic_polynomial(self):
        FqXT = PolynomialRing(self.polring(), 'T')
        return FqXT([self.frobenius_norm(), self.frobenius_trace(), 1])

    def delta(self):
        if self.rank() != 2:
            raise ValueError('The rank must equal 2')
        return self.gen()[2]

    def frobenius_norm(self):
        # Notations from Schost-Musleh:
        n = self._L().over(self._Fq()).degree_over(self._Fq())
        d = self.characteristic().degree()
        m = n // d
        norm = self._L().over(self._Fq())(self.delta()).norm()
        return ((-1)**n) * (self.characteristic()**m) / norm

    def frobenius_trace(self):
        # Notations from Schost-Musleh:
        n = self._L().over(self._Fq()).degree_over(self._Fq())
        B = self.frobenius_norm()
        t = self.ore_polring().gen()
        return self.invert(t**n + (self(B) // t**n))

    def g(self):
        if self.rank() != 2:
            raise ValueError('The rank must equal 2')
        return self.gen()[1]

    def height(self):
        return Integer(1)

    def invert(self, image):
        """
        Given an Ore polynomial `image` of the form `phi(c)`, find c.
        """
        if not image in self.ore_polring():
            raise TypeError('The tested image should be in the Ore ' \
                    'polynomial ring')
        if image in self._L():  # Only works if `image` is in the image of self
            return self._Fq()(image)
        r = self.rank()
        X = self.polring().gen()
        k = image.degree() // r
        m_lines = [[0 for _ in range(k+1)] for _ in range(k+1)]
        for i in range(k+1):
            phi_X_i = self(X**i)
            for j in range(i+1):
                m_lines[j][i] = phi_X_i[r*j]
        m = Matrix(m_lines)
        v = vector([list(image)[r*j] for j in range(k+1)])
        pre_image = self.polring()(list((m**(-1)) * v))
        if self(pre_image) == image:
            return pre_image
        else:
            return None

    def is_supersingular(self):
        return self.characteristic().divides(self.frobenius_trace())

    def is_ordinary(self):
        return not self.is_supersingular()

    def is_morphism(self, candidate):
        return candidate == 0 or self.is_isogeny(candidate)

    def is_isogeny(self, candidate):
        if not candidate in self.ore_polring():
            raise TypeError('The candidate must be in the Ore polynomial ' \
                    'ring')
        if candidate == 0:
            return False
        elif candidate in self.ore_polring().base_ring():
            return True
        else:
            return self.characteristic().degree().divides(candidate.valuation()) \
                and candidate.right_divides(candidate * self.gen())

    def is_endomorphism(self, candidate):
        return candidate == 0 or self == self.velu(candidate)

    def is_automorphism(self, candidate):
        if not candidate in self.ore_polring():
            raise TypeError('The candidate must be in the Ore polynomial ' \
                    'ring')
        return candidate != 0 and candidate in self._Fq()

    def j(self):
        if self.rank() != 2:
            raise ValueError('The j-invariant is only defined for rank 2 ' \
                    'Drinfeld modules')
        g = self.gen()[1]
        delta = self.gen()[2]
        q = self.polring().base_ring().order()
        return (g**(q+1)) / delta

    def rank(self):
        return self.gen().degree()

    def velu(self, candidate):
        if not candidate in self.ore_polring():
            raise TypeError('The candidate must be in the Ore polynomial ' \
                    'ring')
        # There are two main ways to give the result. The first way is
        # to return the Drinfeld module generated by the right-quotient
        # of `candidate * self(X)` right-divided by `candidate`. The
        # second way is to recursively find the coefficients (see
        # arXiv:2203.06970, Eq. 1.1). For now, the former is
        # implemented, as it is very easy to write.
        if candidate == 0:
            return None
        if not self.characteristic().degree().divides(candidate.valuation()):
            return None
        q, r = (candidate * self.gen()).right_quo_rem(candidate)
        if r != 0:
            return None
        else:
            return FiniteDrinfeldModule(self.polring(), q, self.characteristic())

    ##########################
    # Special Sage functions #
    ##########################

    def _get_action_(self):
        return FiniteDrinfeldModuleAction(self)

    def _latex_(self):
        return f'\\text{{Finite{{ }}Drinfeld{{ }}module{{ }}defined{{ }}by{{ }}}}\n' \
                f'\\begin{{align}}\n' \
                f'  {latex(self.polring())}\n' \
                f'  &\\to {latex(self.ore_polring())} \\\\\n' \
                f'  {latex(self.polring().gen())}\n' \
                f'  &\\mapsto {latex(self.gen())}\n' \
                f'\\end{{align}}\n' \
                f'\\text{{with{{ }}characteristic{{ }}}} ' \
                f'{latex(self.characteristic())}'

    def _repr_(self):
        return f'Finite Drinfeld module:\n' \
                f'  Polring:        {self.polring()}\n' \
                f'  Ore polring:    {self.ore_polring()}\n' \
                f'  Generator:      {self.gen()}\n' \
                f'  Characteristic: {self.characteristic()}'

    ###########
    # Getters #
    ###########

    def characteristic(self):
        return self.__characteristic

    def constant_term(self):
        return self.gen()[0]

    def frobenius(self):
        return self.ore_polring().twisting_morphism()

    def gen(self):
        [gen] = self.im_gens()
        return gen

    def ore_polring(self):
        return self.codomain()

    def polring(self):
        return self.domain()


class FiniteDrinfeldModuleAction(Action):

    def __init__(self, finite_drinfeld_module):
        # Verifications
        if not isinstance(finite_drinfeld_module, FiniteDrinfeldModule):
            raise TypeError('First argument must be a FiniteDrinfeldModule')
        # Work
        self.__finite_drinfeld_module = finite_drinfeld_module
        super().__init__(finite_drinfeld_module.polring(),
                finite_drinfeld_module.ore_polring().base_ring())

    ###########
    # Methods #
    ###########

    def finite_drinfeld_module(self):
        return self.__finite_drinfeld_module

    ##########################
    # Special Sage functions #
    ##########################

    def _latex_(self):
        return f'\\text{{Action{{ }}on{{ }}}}' \
                f'{latex(self.extension())}\\text{{{{ }}' \
                f'induced{{ }}by{{ }}}}{self.finite_drinfeld_module()}'

    def _repr_(self):
        return f'Action on {self.domain()} induced by ' \
                f'{self.finite_drinfeld_module()}'

    def _act_(self, g, x):
        return self.finite_drinfeld_module()(g)(x)


def _check_base_fields(Fq, L):
    if not (L.is_field() and L.is_finite() and Fq.is_subring(L)):
        raise ValueError(f'The base field of the Ore polynomial ring must ' \
                'be a finite field extension of the base field of the ' \
                'polynomial ring')
