r"""
The set of morphism between ChowRings

Space of ChowRing homomorphisms.

Derives from the space of quotient ring homomorphisms implement in
`sage.ring.homset.RingHomset_quo_ring` in order to explicitly allow
"no generators", e.g. im_gens = [].


EXAMPLES::

    sage: A.<h> = ChowRing('h', 1, 'h^3')
    sage: B.<k> = ChowRing('k', 1, 'k^6')
    sage: phi = B.hom([2*h], A); phi
    Ring morphism:
      From: Quotient of Multivariate Polynomial Ring in k over Rational Field by the ideal (k^6)
      To:   Quotient of Multivariate Polynomial Ring in h over Rational Field by the ideal (h^3)
      Defn: k |--> 2*h
    sage: phi(2)
    2
    sage: phi(k)
    2*h

    sage: A = ChowRing()
    sage: f = A.hom([], A)
    sage: f(1)
    1

    sage: B.<h> = ChowRing('h', 1, 'h^3')
    sage: f = B.hom([A(0)], A)
    sage: f(h)
    0
    sage: f(3)
    3

TESTS::

    sage: A = ChowRing()
    sage: B.<h> = ChowRing('h', 1, 'h^3')
    sage: H = B.Hom(A)
    sage: H == loads(dumps(H))
    True

    sage: TestSuite(phi).run()

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2023)
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2023 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.homset import RingHomset_generic
from sage.rings import morphism


class ChowRingHomSet(RingHomset_generic):
    """
    sage: CR.<h> = ChowRing('h', 1, 'h^2')
    sage: CS.<k> = ChowRing('k', 1, 'k^4')
    sage: cf = CS.hom([3*h], CR)
    sage: cg = loads(cf.dumps())
    sage: cf == cg
    True
    """
    Element = morphism.RingHomomorphism_from_quotient

    def _element_constructor_(self, x, base_map=None, check=True):
        """
        Construct an element of ``self`` from ``x``.
        """
        if isinstance(x, morphism.RingHomomorphism_from_quotient):
            phi = x._phi()
        else:
            pi = self.domain().cover()
            if not x:
                # We explicitly allow "no generators", e.g. x = [], which is
                # not the case in the original RingHomset_quo_ring implementation.
                phi = pi.domain().hom([], self.codomain(), base_map=base_map, check=check)
            else:
                phi = pi.domain().hom(x, base_map=base_map, check=check)
        return self.element_class(self, phi)
