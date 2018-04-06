# -*- coding: utf-8 -*-
r"""
Class to flatten polynomial rings over polynomial ring

For example ``QQ['a','b'],['x','y']`` flattens to ``QQ['a','b','x','y']``.


EXAMPLES::

    sage: R = QQ['x']['y']['s','t']['X']
    sage: from sage.rings.polynomial.flatten import FlatteningMorphism
    sage: phi = FlatteningMorphism(R); phi
    Flattening morphism:
      From: Univariate Polynomial Ring in X over Multivariate Polynomial Ring in s, t over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
      To:   Multivariate Polynomial Ring in x, y, s, t, X over Rational Field
    sage: phi('x*y*s + t*X').parent()
    Multivariate Polynomial Ring in x, y, s, t, X over Rational Field


Authors:

Vincent Delecroix, Ben Hutz (July 2016): initial implementation
"""

#*****************************************************************************
#                  Copyright (C) 2016
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function

from sage.misc.cachefunc import cached_method

from sage.categories.morphism import Morphism
import six

from .polynomial_ring_constructor import PolynomialRing
from .polynomial_ring import is_PolynomialRing
from .multi_polynomial_ring_generic import is_MPolynomialRing

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
        The Python constructor

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
            Traceback (most recent call last):
            ...
            ValueError: clash in variable names

        ::

            sage: K.<v> = NumberField(x^3 - 2)
            sage: R = K['x','y']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f(R('v*a*x^2 + b^2 + 1/v*y'))
            (v)*x^2*a + b^2 + (1/2*v^2)*y

        ::

            sage: R = QQbar['x','y']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f(R('QQbar(sqrt(2))*a*x^2 + b^2 + QQbar(I)*y'))
            1.414213562373095?*x^2*a + b^2 + I*y

        ::

            sage: R.<z> = PolynomialRing(QQbar,1)
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f.domain(), f.codomain()
            (Multivariate Polynomial Ring in z over Algebraic Field,
             Multivariate Polynomial Ring in z over Algebraic Field)

        ::

            sage: R.<z> = PolynomialRing(QQbar)
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f.domain(), f.codomain()
            (Univariate Polynomial Ring in z over Algebraic Field,
             Univariate Polynomial Ring in z over Algebraic Field)
        """
        if not is_PolynomialRing(domain) and not is_MPolynomialRing(domain):
            raise ValueError("domain should be a polynomial ring")

        ring = domain
        variables = []
        intermediate_rings = []

        while is_PolynomialRing(ring) or is_MPolynomialRing(ring):
            intermediate_rings.append(ring)
            v = ring.variable_names()
            if any(vv in variables for vv in v):
                raise ValueError("clash in variable names")
            variables.extend(reversed(v))
            ring = ring.base_ring()
        self._intermediate_rings = intermediate_rings
        variables.reverse()
        if is_MPolynomialRing(domain):
            codomain = PolynomialRing(ring, variables, len(variables))
        else:
            codomain = PolynomialRing(ring, variables)

        Morphism.__init__(self, domain, codomain)
        self._repr_type_str = 'Flattening'

    def _call_(self, p):
        r"""
        Evaluate an flattening morphism.

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
         #If we are just specializing a univariate polynomial, then
         #the flattening morphism is the identity
        if self.codomain().ngens()==1:
            return p

        p = {(): p}

        for ring in self._intermediate_rings:
            new_p = {}
            if is_PolynomialRing(ring):
                for mon,pp in six.iteritems(p):
                    assert pp.parent() == ring
                    for i,j in six.iteritems(pp.dict()):
                        new_p[(i,)+(mon)] = j
            elif is_MPolynomialRing(ring):
                for mon,pp in six.iteritems(p):
                    assert pp.parent() == ring
                    for mmon,q in six.iteritems(pp.dict()):
                        new_p[tuple(mmon)+mon] = q
            else:
                raise RuntimeError
            p = new_p

        return self.codomain()(p)

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
              To:   Multivariate Polynomial Ring in x, y, z over Multivariate Polynomial Ring in a, b, c over Rational Field

        ::

            sage: R = ZZ['a']['b']['c']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R).section()
            Unflattening morphism:
              From: Multivariate Polynomial Ring in a, b, c over Integer Ring
              To:   Univariate Polynomial Ring in c over Univariate Polynomial Ring in b over Univariate Polynomial Ring in a over Integer Ring
        """
        phi= UnflatteningMorphism(self.codomain(), self.domain())
        return phi

class UnflatteningMorphism(Morphism):
    r"""
    Inverses for :class:`FlatteningMorphism`

    EXAMPLES::

        sage: R = QQ['c','x','y','z']
        sage: S = QQ['c']['x','y','z']
        sage: from sage.rings.polynomial.flatten import UnflatteningMorphism
        sage: f = UnflatteningMorphism(R, S)
        sage: g = f(R('x^2 + c*y^2 - z^2'));g
        x^2 + c*y^2 - z^2
        sage: g.parent()
        Multivariate Polynomial Ring in x, y, z over Univariate Polynomial Ring in c over Rational Field

    ::

        sage: R = QQ['a','b', 'x','y']
        sage: S = QQ['a','b']['x','y']
        sage: from sage.rings.polynomial.flatten import UnflatteningMorphism
        sage: UnflatteningMorphism(R, S)
        Unflattening morphism:
          From: Multivariate Polynomial Ring in a, b, x, y over Rational Field
          To:   Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring in a, b over Rational Field
    """

    def __init__(self, domain, codomain):
        """
        The Python constructor

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
        if not is_MPolynomialRing(domain):
            raise ValueError("domain should be a multivariate polynomial ring")
        if not is_PolynomialRing(codomain) and not is_MPolynomialRing(codomain):
            raise ValueError("codomain should be a polynomial ring")

        ring = codomain
        intermediate_rings = []

        while is_PolynomialRing(ring) or is_MPolynomialRing(ring):
            intermediate_rings.append(ring)
            ring = ring.base_ring()

        if domain.base_ring() != intermediate_rings[-1].base_ring():
            raise ValueError("rings must have same base ring")
        if domain.ngens() != sum([R.ngens() for R in intermediate_rings]):
            raise ValueError("rings must have the same number of variables")

        self._intermediate_rings = intermediate_rings
        self._intermediate_rings.reverse()

        Morphism.__init__(self, domain, codomain)
        self._repr_type_str = 'Unflattening'

    def _call_(self, p):
        """
        Evaluate an unflattening morphism.

        TESTS::

            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: for R in [ZZ['x']['y']['a,b,c'], GF(4)['x','y']['a','b'],
            ....:           AA['x']['a','b']['y'], QQbar['a1','a2']['t']['X','Y']]:
            ....:    f = FlatteningMorphism(R)
            ....:    g = f.section()
            ....:    for _ in range(10):
            ....:        p = R.random_element()
            ....:        assert p == g(f(p))
        """
        num = [len(R.gens()) for R in self._intermediate_rings]
        f = self.codomain().zero()
        for mon,pp in six.iteritems(p.dict()):
            ind = 0
            g = pp
            for i in range(len(num)):
                m = mon[ind:ind+num[i]]
                ind += num[i]
                R = self._intermediate_rings[i]
                if is_PolynomialRing(R):
                    m = m[0]
                g = R({m: g})
            f += g

        return f


class SpecializationMorphism(Morphism):
    r"""
    Morphisms to specialize parameters in (stacked) polynomial rings

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
              From: Univariate Polynomial Ring in z over Univariate Polynomial Ring in c over Rational Field
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
          From: Multivariate Polynomial Ring in x, y, z over Multivariate Polynomial Ring in a, b, c over Multivariate Polynomial Ring in u, v over Rational Field
          To:   Multivariate Polynomial Ring in y, z over Univariate Polynomial Ring in c over Univariate Polynomial Ring in v over Rational Field
        sage: xi(a*(x*z+y^2)*u+b*v*u*(x*z+y^2)*y^2*c+c*y^2*z^2)
        2*v*c*y^4 + c*y^2*z^2 + y^2
    """

    def __init__(self, domain, D):
        """
        The Python constructor

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

        The following was fixed in :trac:`23811`::

            sage: R.<c>=RR[]
            sage: P.<z>=AffineSpace(R,1)
            sage: H=End(P)
            sage: f=H([z^2+c])
            sage: f.specialization({c:1})
            Scheme endomorphism of Affine Space of dimension 1 over Real Field with 53 bits of precision
              Defn: Defined on coordinates by sending (z) to
                    (z^2 + 1.00000000000000)

        """
        if not is_PolynomialRing(domain) and not is_MPolynomialRing(domain):
            raise TypeError("domain should be a polynomial ring")

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
        while is_PolynomialRing(R) or is_MPolynomialRing(R):
            old = R.gens()
            new = [t for t in old if t not in D]
            force_multivariate = ((len(old) == 1) and is_MPolynomialRing(R))
            new_vars.append((new, force_multivariate))
            R = R.base_ring()

        # Construct unflattening map psi (only defined on the variables
        # of "flat" which are not involved in D)
        psi = dict()
        for new, force_multivariate in reversed(new_vars):
            if not new:
                continue
            # Pass in the names of the variables
            var_names = [str(var) for var in new]
            if force_multivariate:
                R = PolynomialRing(R, var_names, len(var_names))
            else:
                R = PolynomialRing(R, var_names)
            # Map variables in "new" to R
            psi.update(zip([phi(w) for w in new], R.gens()))

        # Compose D with psi
        vals = []
        for t in flat.gens():
            if t in D:
                vals.append(R.coerce(D[t]))
            else:
                vals.append(psi[t])

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
        return self._eval_morph(self._flattening_morph(p))
