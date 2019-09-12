r"""
General extensions of p-adic rings and fields; the base ring may also be an
extension.

These are implemented as proxy parents, backed by an absolute extension.
"""
#*****************************************************************************
#       Copyright (C) 2019 David Roe <roed.math@gmail.com>
#                          Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from .extension_element import pAdicGenericExtensionElement
from .padic_extension_generic import pAdicExtensionGeneric

# NotImplementedError: AlgebraFromMorphism shouldn't inherit from UniqueRepresentation
class pAdicGeneralExtension(AlgebraFromMorphism, pAdicExtensionGeneric):
    def __init__(self, exact_modulus, poly, prec, print_mode, shift_seed, names, implementation='FLINT'):
        self._exact_modulus = exact_modulus
        self._shift_seed = shift_seed
        self._implementation = 'proxy'
        defining_morphism = None # NotImplementedError
        AlgebraFromMorphism.__init__(self, defining_morphism, False)
        pAdicGeneric.__init__(self, poly, prec, print_mode, names, pAdicGenericExtensionElement)
        # Fix the following
        self._gen = None

    @cached_method
    def f(self):
        r"""
        Return the residual degree of this ring over its base ring.

        EXAMPLES::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.f()
            2

        """
        v = self.base_ring().valuation()
        w = v.mac_lane_approximants(self.modulus(), assume_squarefree=True, require_final_EF=True)
        if len(w) != 1:
            raise ValueError("defining polynomial is not irreducible")
        return w[0].F()

    def e(self):
        r"""
        Return the ramification degree of this ring over its base ring.

        EXAMPLES::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.e()
            1

        """
        return self.relative_degree() // self.f()

    def teichmuller(self, x, prec=None):
        R = self._backend()
        x = R(x) if prec is None else R(x, prec)
        return self(R.teichmuller(x))

    def _prec_type(self):
        return self._backend()._prec_type()

    def is_field(self):
        return self._backend().is_field()

    def random_element(self, **kwds):
        return self(self._backend().random_element(**kwds))

    def residue_ring(self, n):
        raise NotImplementedError

    def residue_class_field(self):
        raise NotImplementedError

    def inertia_subring(self):
        raise NotImplementedError

    def gen(self, n=0):
        if n != 0:
            raise IndexError("only one generator")
        return self._gen

    def uniformizer(self):
        return self(self._backend().uniformizer())

    def uniformizer_pow(self, n):
        return self(self._backend().uniformizer_pow(n))

    def _uniformizer_print(self):
        return self._backend()._uniformizer_print()

    def _unram_print(self):
        return self._backend()._unram_print()

    def has_pth_root(self):
        return self._backend().has_pth_root()

    def has_root_of_unity(self, n):
        return self._backend().has_root_of_unity(self, n)

    
