# sage_setup: distribution = sagemath-categories
r"""
Class to flatten polynomial rings over polynomial ring

For example ``QQ['a','b'],['x','y']`` flattens to ``QQ['a','b','x','y']``.

EXAMPLES::

    sage: R = QQ['x']['y']['s','t']['X']
    sage: from sage.rings.polynomial.flatten import FlatteningMorphism
    sage: phi = FlatteningMorphism(R); phi
    Flattening morphism:
      From: Univariate Polynomial Ring in X
            over Multivariate Polynomial Ring in s, t
            over Univariate Polynomial Ring in y
            over Univariate Polynomial Ring in x over Rational Field
      To:   Multivariate Polynomial Ring in x, y, s, t, X over Rational Field
    sage: phi('x*y*s + t*X').parent()
    Multivariate Polynomial Ring in x, y, s, t, X over Rational Field

Authors:

Vincent Delecroix, Ben Hutz (July 2016): initial implementation
"""

# ****************************************************************************
#                  Copyright (C) 2016
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools

from sage.categories.homset import Homset
from sage.categories.morphism import Morphism
from sage.misc.cachefunc import cached_method
from .polynomial_ring_constructor import PolynomialRing
from .polynomial_ring import PolynomialRing_general
from .multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.fraction_field import FractionField_generic
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.polynomial.polydict import ETuple


class FlatteningMorphism(Morphism):
    r"""
    EXAMPLES::

        sage: R = QQ['a','b']['x','y','z']['t1','t2']
        sage: from sage.rings.polynomial.flatten import FlatteningMorphism
        sage: f = FlatteningMorphism(R)
        sage: f.codomain()
        Multivariate Polynomial Ring in a, b, x, y, z, t1, t2 over Rational Field
        sage: p = R('(a+b)*x + (a^2-b)*t2*(z+y)')
        sage: p
        ((a^2 - b)*y + (a^2 - b)*z)*t2 + (a + b)*x
        sage: f(p)
        a^2*y*t2 + a^2*z*t2 - b*y*t2 - b*z*t2 + a*x + b*x
        sage: f(p).parent()
        Multivariate Polynomial Ring in a, b, x, y, z, t1, t2 over Rational Field

    Also works when univariate polynomial ring are involved::

        sage: R = QQ['x']['y']['s','t']['X']
        sage: from sage.rings.polynomial.flatten import FlatteningMorphism
        sage: f = FlatteningMorphism(R)
        sage: f.codomain()
        Multivariate Polynomial Ring in x, y, s, t, X over Rational Field
        sage: p = R('((x^2 + 1) + (x+2)*y + x*y^3)*(s+t) + x*y*X')
        sage: p
        x*y*X + (x*y^3 + (x + 2)*y + x^2 + 1)*s + (x*y^3 + (x + 2)*y + x^2 + 1)*t
        sage: f(p)
        x*y^3*s + x*y^3*t + x^2*s + x*y*s + x^2*t + x*y*t + x*y*X + 2*y*s + 2*y*t + s + t
        sage: f(p).parent()
        Multivariate Polynomial Ring in x, y, s, t, X over Rational Field
    """
    def __init__(self, domain):
        """
        The Python constructor.

        EXAMPLES::

            sage: R = ZZ['a', 'b', 'c']['x', 'y', 'z']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R)
            Flattening morphism:
              From: Multivariate Polynomial Ring in x, y, z over Multivariate Polynomial Ring in a, b, c over Integer Ring
              To:   Multivariate Polynomial Ring in a, b, c, x, y, z over Integer Ring

        ::

            sage: R = ZZ['a']['b']['c']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R)
            Flattening morphism:
              From: Univariate Polynomial Ring in c over Univariate Polynomial Ring in b over Univariate Polynomial Ring in a over Integer Ring
              To:   Multivariate Polynomial Ring in a, b, c over Integer Ring

        ::

            sage: R = ZZ['a']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R)
            Flattening morphism:
              From: Multivariate Polynomial Ring in a, b over Univariate Polynomial Ring in a over Integer Ring
              To:   Multivariate Polynomial Ring in a, a0, b over Integer Ring

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<v> = NumberField(x^3 - 2)
            sage: R = K['x','y']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f(R('v*a*x^2 + b^2 + 1/v*y'))
            v*x^2*a + b^2 + (1/2*v^2)*y

        ::

            sage: # needs sage.rings.number_field
            sage: R = QQbar['x','y']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f(R('QQbar(sqrt(2))*a*x^2 + b^2 + QQbar(I)*y'))                       # needs sage.symbolic
            1.414213562373095?*x^2*a + b^2 + I*y

        ::

            sage: # needs sage.rings.number_field
            sage: R.<z> = PolynomialRing(QQbar, 1)
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f.domain(), f.codomain()
            (Multivariate Polynomial Ring in z over Algebraic Field,
             Multivariate Polynomial Ring in z over Algebraic Field)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<z> = PolynomialRing(QQbar)
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f.domain(), f.codomain()
            (Univariate Polynomial Ring in z over Algebraic Field,
             Univariate Polynomial Ring in z over Algebraic Field)

        TESTS::

            sage: Pol = QQ['x']['x0']['x']
            sage: fl = FlatteningMorphism(Pol)
            sage: fl
            Flattening morphism:
              From: Univariate Polynomial Ring in x over Univariate Polynomial Ring in x0 over Univariate Polynomial Ring in x over Rational Field
              To:   Multivariate Polynomial Ring in x, x0, x1 over Rational Field
            sage: p = Pol([[[1,2],[3,4]],[[5,6],[7,8]]])
            sage: fl.section()(fl(p)) == p
            True
        """
        if not isinstance(domain, PolynomialRing_general) and not isinstance(domain, MPolynomialRing_base):
            raise ValueError("domain should be a polynomial ring")

        ring = domain
        variables = []
        intermediate_rings = []

        while isinstance(ring, PolynomialRing_general) or isinstance(ring, MPolynomialRing_base):
            intermediate_rings.append(ring)
            v = ring.variable_names()
            variables.extend(reversed(v))
            ring = ring.base_ring()
        self._intermediate_rings = intermediate_rings
        variables.reverse()
        for i, a in enumerate(variables):
            if a in variables[:i]:
                for index in itertools.count():
                    b = a + str(index)
                    if b not in variables:  # not just variables[:i]!
                        break
                variables[i] = b
        if isinstance(domain, MPolynomialRing_base):
            codomain = PolynomialRing(ring, variables, len(variables))
        else:
            codomain = PolynomialRing(ring, variables)

        hom = Homset(domain, codomain, base=ring, check=False)
        Morphism.__init__(self, hom)
        self._repr_type_str = 'Flattening'

    def _call_(self, p):
        r"""
        Evaluate a flattening morphism.

        EXAMPLES::

            sage: R = QQ['a','b','c']['x','y','z']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: h = FlatteningMorphism(R)('2*a*x + b*z'); h
            2*a*x + b*z
            sage: h.parent()
            Multivariate Polynomial Ring in a, b, c, x, y, z over Rational Field

        TESTS::

            sage: R = QQ['x']['y']['s','t']
            sage: p = R('s*x + y*t + x^2*s + 1 + t')
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f._call_(p)
            x^2*s + x*s + y*t + t + 1
        """
        # If we are just specializing a univariate polynomial, then
        # the flattening morphism is the identity
        if self.codomain().ngens() == 1:
            return p

        p = {(): p}

        for ring in self._intermediate_rings:
            new_p = {}
            if isinstance(ring, PolynomialRing_general):
                for mon, pp in p.items():
                    assert pp.parent() is ring
                    for i, j in pp.monomial_coefficients().items():
                        new_p[(i,) + (mon)] = j
            elif isinstance(ring, MPolynomialRing_base):
                for mon, pp in p.items():
                    assert pp.parent() is ring
                    for mmon, q in pp.monomial_coefficients().items():
                        new_p[tuple(mmon) + mon] = q
            else:
                raise RuntimeError
            p = new_p

        return self.codomain()(p, check=False)

    @cached_method
    def section(self):
        """
        Inverse of this flattening morphism.

        EXAMPLES::

            sage: R = QQ['a','b','c']['x','y','z']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: h = FlatteningMorphism(R)
            sage: h.section()
            Unflattening morphism:
              From: Multivariate Polynomial Ring in a, b, c, x, y, z over Rational Field
              To:   Multivariate Polynomial Ring in x, y, z
                    over Multivariate Polynomial Ring in a, b, c over Rational Field

        ::

            sage: R = ZZ['a']['b']['c']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R).section()
            Unflattening morphism:
              From: Multivariate Polynomial Ring in a, b, c over Integer Ring
              To:   Univariate Polynomial Ring in c over Univariate Polynomial Ring in b
                    over Univariate Polynomial Ring in a over Integer Ring
        """
        return UnflatteningMorphism(self.codomain(), self.domain())

    def inverse(self):
        """
        Return the inverse of this flattening morphism.

        This is the same as calling :meth:`section`.

        EXAMPLES::

            sage: f = QQ['x,y']['u,v'].flattening_morphism()
            sage: f.inverse()
            Unflattening morphism:
              From: Multivariate Polynomial Ring in x, y, u, v over Rational Field
              To:   Multivariate Polynomial Ring in u, v
                    over Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self.section()


class UnflatteningMorphism(Morphism):
    r"""
    Inverses for :class:`FlatteningMorphism`.

    EXAMPLES::

        sage: R = QQ['c','x','y','z']
        sage: S = QQ['c']['x','y','z']
        sage: from sage.rings.polynomial.flatten import UnflatteningMorphism
        sage: f = UnflatteningMorphism(R, S)
        sage: g = f(R('x^2 + c*y^2 - z^2'));g
        x^2 + c*y^2 - z^2
        sage: g.parent()
        Multivariate Polynomial Ring in x, y, z
         over Univariate Polynomial Ring in c over Rational Field

    ::

        sage: R = QQ['a','b', 'x','y']
        sage: S = QQ['a','b']['x','y']
        sage: from sage.rings.polynomial.flatten import UnflatteningMorphism
        sage: UnflatteningMorphism(R, S)
        Unflattening morphism:
          From: Multivariate Polynomial Ring in a, b, x, y over Rational Field
          To:   Multivariate Polynomial Ring in x, y
                over Multivariate Polynomial Ring in a, b over Rational Field
    """

    def __init__(self, domain, codomain):
        """
        The Python constructor.

        EXAMPLES::

            sage: R = QQ['x']['y']['s','t']['X']
            sage: p = R.random_element()
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: g = f.section()
            sage: g(f(p)) == p
            True

        ::

            sage: R = QQ['a','b','x','y']
            sage: S = ZZ['a','b']['x','z']
            sage: from sage.rings.polynomial.flatten import UnflatteningMorphism
            sage: UnflatteningMorphism(R, S)
            Traceback (most recent call last):
            ...
            ValueError: rings must have same base ring

        ::

            sage: R = QQ['a','b','x','y']
            sage: S = QQ['a','b']['x','z','w']
            sage: from sage.rings.polynomial.flatten import UnflatteningMorphism
            sage: UnflatteningMorphism(R, S)
            Traceback (most recent call last):
            ...
            ValueError: rings must have the same number of variables
        """
        if not isinstance(domain, MPolynomialRing_base):
            raise ValueError("domain should be a multivariate polynomial ring")
        if not isinstance(codomain, PolynomialRing_general) and not isinstance(codomain, MPolynomialRing_base):
            raise ValueError("codomain should be a polynomial ring")

        ring = codomain
        intermediate_rings = []

        while True:
            is_polynomial_ring = isinstance(ring, PolynomialRing_general)
            if not (is_polynomial_ring or isinstance(ring, MPolynomialRing_base)):
                break
            intermediate_rings.append((ring, is_polynomial_ring))
            ring = ring.base_ring()

        if domain.base_ring() != intermediate_rings[-1][0].base_ring():
            raise ValueError("rings must have same base ring")
        if domain.ngens() != sum([R.ngens() for R, _ in intermediate_rings]):
            raise ValueError("rings must have the same number of variables")

        self._intermediate_rings = intermediate_rings

        hom = Homset(domain, codomain, base=ring, check=False)
        Morphism.__init__(self, hom)
        self._repr_type_str = 'Unflattening'

    def _call_(self, p):
        """
        Evaluate an unflattening morphism.

        TESTS::

            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: rings = [ZZ['x']['y']['a,b,c']]
            sage: rings += [GF(4)['x','y']['a','b']]                                    # needs sage.rings.finite_rings
            sage: rings += [AA['x']['a','b']['y'], QQbar['a1','a2']['t']['X','Y']]      # needs sage.rings.number_field
            sage: for R in rings:                                                       # needs sage.modules
            ....:    f = FlatteningMorphism(R)
            ....:    g = f.section()
            ....:    for _ in range(10):
            ....:        p = R.random_element()
            ....:        assert p == g(f(p))
            ....:        z = R.zero()
            ....:        assert z == g(f(z))
        """
        index = [0]
        for R, _ in reversed(self._intermediate_rings):
            index.append(index[-1] + len(R.gens()))
        newpol = [{} for _ in self._intermediate_rings]
        expo = sorted(p.exponents(), key=lambda e: tuple(reversed(e)))
        for i in range(len(expo)):
            cur_exp = expo[i]
            for l in range(len(self._intermediate_rings)):
                R, univariate = self._intermediate_rings[-1 - l]
                idx = index[l + 1]
                sub_exp = (cur_exp[index[l]] if univariate
                           else cur_exp[index[l]:idx])
                if l == 0:
                    newpol[l][sub_exp] = p[cur_exp]
                else:
                    newpol[l][sub_exp] = newpol[l - 1]
                    newpol[l - 1] = {}
                if (i == len(expo) - 1 or expo[i + 1][idx:] != cur_exp[idx:]):
                    newpol[l] = R(newpol[l], check=False)
                else:
                    break
        return R(newpol[-1], check=False)


class SpecializationMorphism(Morphism):
    r"""
    Morphisms to specialize parameters in (stacked) polynomial rings.

    EXAMPLES::

        sage: R.<c> = PolynomialRing(QQ)
        sage: S.<x,y,z> = PolynomialRing(R)
        sage: D = dict({c:1})
        sage: from sage.rings.polynomial.flatten import SpecializationMorphism
        sage: f = SpecializationMorphism(S, D)
        sage: g = f(x^2 + c*y^2 - z^2); g
        x^2 + y^2 - z^2
        sage: g.parent()
        Multivariate Polynomial Ring in x, y, z over Rational Field

    ::

        sage: R.<c> = PolynomialRing(QQ)
        sage: S.<z> = PolynomialRing(R)
        sage: from sage.rings.polynomial.flatten import SpecializationMorphism
        sage: xi = SpecializationMorphism(S, {c:0}); xi
        Specialization morphism:
          From: Univariate Polynomial Ring in z
                over Univariate Polynomial Ring in c over Rational Field
          To:   Univariate Polynomial Ring in z over Rational Field
        sage: xi(z^2+c)
        z^2

    ::

        sage: R1.<u,v> = PolynomialRing(QQ)
        sage: R2.<a,b,c> = PolynomialRing(R1)
        sage: S.<x,y,z> = PolynomialRing(R2)
        sage: D = dict({a:1, b:2, x:0, u:1})
        sage: from sage.rings.polynomial.flatten import SpecializationMorphism
        sage: xi = SpecializationMorphism(S, D); xi
        Specialization morphism:
          From: Multivariate Polynomial Ring in x, y, z
                over Multivariate Polynomial Ring in a, b, c
                over Multivariate Polynomial Ring in u, v over Rational Field
          To:   Multivariate Polynomial Ring in y, z over Univariate Polynomial Ring in c
                over Univariate Polynomial Ring in v over Rational Field
        sage: xi(a*(x*z+y^2)*u+b*v*u*(x*z+y^2)*y^2*c+c*y^2*z^2)
        2*v*c*y^4 + c*y^2*z^2 + y^2
    """

    def __init__(self, domain, D):
        """
        The Python constructor.

        EXAMPLES::

            sage: S.<x,y> = PolynomialRing(QQ)
            sage: D = dict({x:1})
            sage: from sage.rings.polynomial.flatten import SpecializationMorphism
            sage: phi = SpecializationMorphism(S, D); phi
            Specialization morphism:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Univariate Polynomial Ring in y over Rational Field
            sage: phi(x^2 + y^2)
            y^2 + 1

        ::

            sage: R.<a,b,c> = PolynomialRing(ZZ)
            sage: S.<x,y,z> = PolynomialRing(R)
            sage: from sage.rings.polynomial.flatten import SpecializationMorphism
            sage: xi = SpecializationMorphism(S, {a:1/2})
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        The following was fixed in :issue:`23811`::

            sage: R.<c> = RR[]
            sage: P.<z> = AffineSpace(R, 1)
            sage: H = End(P)
            sage: f = H([z^2 + c])
            sage: f.specialization({c:1})                                               # needs sage.modules
            Scheme endomorphism of
             Affine Space of dimension 1 over Real Field with 53 bits of precision
              Defn: Defined on coordinates by sending (z) to
                    (z^2 + 1.00000000000000)
        """
        if not isinstance(domain, PolynomialRing_general) and not isinstance(domain, MPolynomialRing_base):
            raise TypeError("domain should be a polynomial ring")

        # use only the generators that are in the stack somewhere,
        # and ignore the rest
        all_gens = domain.gens_dict_recursive()
        new_D = {}
        for gen in D:
            if str(gen) in all_gens:
                new_D[gen] = D[gen]
        D = new_D

        # _sub_specialization is a specialization morphism (recursive)
        # which is applied to the base Fraction field, or None if it's
        # any other base ring

        self._sub_specialization = None

        # We use this composition where "flat" is a flattened
        # polynomial ring.
        #
        #            phi       D       psi
        #     domain  →  flat  →  flat  →  R
        #        │         │               │
        #        └─────────┴───────────────┘
        # _flattening_morph     _eval_morph
        #             = phi       = psi ∘ D

        phi = FlatteningMorphism(domain)
        flat = phi.codomain()
        base = flat.base_ring()

        # Change domain of D to "flat" and ensure that the values lie
        # in the base ring.
        D = {phi(k): base(D[k]) for k in D}

        # Construct unflattened codomain R
        new_vars = []
        R = domain
        while isinstance(R, PolynomialRing_general) or isinstance(R, MPolynomialRing_base) or isinstance(R, FractionField_generic):
            if isinstance(R, FractionField_generic):
                # We've hit base_ring, so set _sub_specialization and exit the loop
                field_over = R.base()
                applicable_vars = {key: val for key, val in D.items()
                                   if key not in flat.gens()}
                # If there are any variables in D to set in _sub_specialization
                if applicable_vars:
                    # Coerce the generators to be in the right ring
                    # This un-does changing the domain of D to be in the flat base ring
                    tmp = {}
                    for var, val in applicable_vars.items():
                        for gstr, gen in field_over.gens_dict_recursive().items():
                            if str(var) == gstr:
                                tmp[gen] = val
                                break
                        else:
                            # Should have been caught earlier
                            raise NameError("argument " + str(var) + " is not a generator anywhere in the polynomial tower")
                    applicable_vars = tmp
                    self._sub_specialization = FractionSpecializationMorphism(R, applicable_vars)
                break
            # We're still in the polynomials, so keep track of the tower
            old = R.gens()
            new = [t for t in old if t not in D]
            force_multivariate = ((len(old) == 1) and isinstance(R, MPolynomialRing_base))
            new_vars.append((new, force_multivariate, old))
            R = R.base_ring()

        if self._sub_specialization:
            # The sub_specialization range will be different
            # if it applied some variables from D
            R = self._sub_specialization.codomain().fraction_field()

        # Construct unflattening map psi (only defined on the variables
        # of "flat" which are not involved in D)
        psi = dict()
        # Reconstruct the proper domain of this morphism
        # based on the sub_specialization domains
        new_domain = R
        for new, force_multivariate, old in reversed(new_vars):
            if self._sub_specialization:
                if force_multivariate:
                    new_domain = PolynomialRing(new_domain, old, len(old))
                else:
                    new_domain = PolynomialRing(new_domain, old)
            if not new:
                continue
            var_names = [str(var) for var in new]
            if force_multivariate:
                R = PolynomialRing(R, var_names, len(var_names))
            else:
                R = PolynomialRing(R, var_names)
            # Map variables in "new" to R
            psi.update(zip([phi(w) for w in new], R.gens()))

        # Fix domain of eval_morph
        # (note: phi's domain is correct)
        if self._sub_specialization:
            phi_prime = FlatteningMorphism(new_domain)
            flat_old = flat
            flat = phi_prime.codomain()
            base_prime = flat.base_ring()
            D = {phi(k): base_prime(D[k]) for k in D}
        else:
            # The bottom of our tower has not changed
            def flat_old(x):
                return x

        # Compose D with psi
        vals = []
        for t in flat.gens():
            if t in D:
                vals.append(R.coerce(D[t]))
            else:
                # Make sure keys are in the old domain
                # or else they won't match exactly
                vals.append(psi[flat_old(t)])

        self._flattening_morph = phi
        self._eval_morph = flat.hom(vals, R)
        self._repr_type_str = 'Specialization'
        Morphism.__init__(self, domain, R)

    def _call_(self, p):
        """
        Evaluate a specialization morphism.

        EXAMPLES::

            sage: R.<a,b,c> = PolynomialRing(ZZ)
            sage: S.<x,y,z> = PolynomialRing(R)
            sage: D = dict({a:1, b:2, c:3})
            sage: from sage.rings.polynomial.flatten import SpecializationMorphism
            sage: xi = SpecializationMorphism(S, D)
            sage: xi(a*x + b*y + c*z)
            x + 2*y + 3*z
        """
        flat = self._flattening_morph(p)
        if self._sub_specialization is not None:
            # The base_ring should be a fraction field, so
            # apply _sub_specialization to each coefficient
            # in the flattened polynomial
            tmp = {}
            for exponent, coefficient in flat.monomial_coefficients().items():
                # Fix the type of exponent from (a,) to a
                #     (necessary for R(tmp) later)
                if isinstance(exponent, ETuple) and len(exponent) == 1:
                    exponent = exponent[0]
                # Coefficient should be a fraction
                tmp[exponent] = self._sub_specialization._call_(coefficient)
            # tmp's parent should be the same construction as flat
            # but over _sub_specialization's codomain
            ring_constructor = flat.parent().construction()[0]
            fraction_type = self._sub_specialization.codomain()
            R = ring_constructor(fraction_type)
            flat = R(tmp)
        return self._eval_morph(flat)


class FractionSpecializationMorphism(Morphism):
    """
    A specialization morphism for fraction fields over (stacked) polynomial rings
    """
    def __init__(self, domain, D):
        """
        Initialize the morphism with a domain and dictionary of specializations.

        EXAMPLES::

            sage: R.<a,c> = QQ[]
            sage: S.<x,y> = R[]
            sage: from sage.rings.polynomial.flatten import FractionSpecializationMorphism
            sage: phi = FractionSpecializationMorphism(Frac(S), {c:3})
            sage: phi
            Fraction Specialization morphism:
              From: Fraction Field of Multivariate Polynomial Ring in x, y
                    over Multivariate Polynomial Ring in a, c over Rational Field
              To:   Fraction Field of Multivariate Polynomial Ring in x, y
                    over Univariate Polynomial Ring in a over Rational Field
        """
        if not isinstance(domain, FractionField_generic):
            raise TypeError("domain must be a fraction field")
        self._specialization = SpecializationMorphism(domain.base(), D)
        self._repr_type_str = 'Fraction Specialization'
        Morphism.__init__(self, domain, self._specialization.codomain().fraction_field())

    def _call_(self, p):
        """
        Evaluate a fraction specialization morphism.

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: S.<x,y,z> = R[]
            sage: from sage.rings.polynomial.flatten import FractionSpecializationMorphism
            sage: phi = FractionSpecializationMorphism(Frac(S), {a:3, b:2, c:-2})
            sage: spec = phi((a*x + b*y) / (c*z))
            sage: spec
            (3*x + 2*y)/(-2*z)
            sage: spec.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        if not isinstance(p, FractionFieldElement):
            raise TypeError("p must be a fraction field element")
        numerator = self._specialization._call_(p.numerator())
        denominator = self._specialization._call_(p.denominator())
        return numerator / denominator
