# -*- coding: utf-8 -*-
"""
Elements of general extensions of p-adic rings and fields; the base ring may also be an extension.

These are implemented as proxy elements, backed by an absolute extension.
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

from copy import deepcopy
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.morphism import RingHomomorphism
from sage.rings.ring_extension_element import RingExtensionElement
from .padic_generic_element import pAdicGenericElement

class pAdicGeneralExtensionElement(RingExtensionElement, pAdicGenericElement):
    # We start with the interesting functions; need to port these to two step extensions
    def polynomial(self, var='x'):
        raise NotImplementedError

    def _poly_rep(self):
        return self.polynomial().change_ring(self.parent()._FP_base())

    # Now a bunch of trivial functions
    def frobenius(self, arithmetic=True):
        return self.__class__(self.parent(), self._element.frobenius(arithmetic=arithmetic))

    def __lshift__(self, shift):
        return self.__class__(self.parent(), self._element << shift)

    def __rshift__(self, shift):
        return self.__class__(self.parent(), self._element >> shift)

    def lift_to_precision(self, absprec=None):
        return self.__class__(self.parent(), self._element.lift_to_precision(absprec=absprec))

    def expansion(self, n=None, lift_mode='simple', start_val=None):
        """
        Write docstring clarifying that this is in terms of the absolute extension
        """
        if lift_mode == 'teichmuller':
            wrap = lambda x: self.__class__(self.parent(), x)
        else:
            wrap = lambda x: x
        E = ExpansionIterable(self, self._element.expansion(n, lift_mode, start_val), wrap, lift_mode)
        if n is None:
            return E
        elif self.parent().is_field():
            v = self.valuation()
            return wrap(E[n-v])
        else:
            return wrap(E[n])

    def _ext_p_list(self, pos):
        return self._element._ext_p_list(pos)

    def unit_part(self):
        return self.__class__(self.parent(), self._element.unit_part())

    def _is_base_elt(self, p):
        return self.parent().prime() == p and self.prime_pow.deg == 1

    def residue(self, absprec=1, field=None, check_prec=True):
        raise NotImplementedError

    def __copy__(self):
        return self.__class__(self.parent(), self._element)

    def __deepcopy__(self):
        return self.__class__(self.parent(), deepcopy(self._element))

    # This and _div_ could be moved to a FieldExtension class
    def __invert__(self):
        return self.__class__(self.parent().fraction_field(), ~self._element)

    def _div_(self, right):
        return self.__class__(self.parent().fraction_field(), self._element / right._element)

    def __pow__(self, right, dummy):
        K = self.parent()
        if isinstance(right, (Integer, Rational, int, long)) and right < 0:
            K = K.fraction_field()
        return self.__class__(K, self._element**right)

    def _quo_rem(self, right):
        q, r = self._element._quo_rem(right._element)
        K = self.parent()
        C = self.__class__
        return C(K, q), C(K, r)

    def add_bigoh(self, absprec):
        return self.__class__(self.parent(), self._element.add_bigoh(absprec))

    def _is_exact_zero(self):
        return self._element._is_exact_zero()

    def _is_inexact_zero(self):
        return self._element._is_inexact_zero()

    def _is_zero_rep(self):
        return self._element._is_zero_rep()

    def is_zero(self, absprec=None):
        return self._element.is_zero(absprec)

    def __nonzero__(self):
        return self._element.__nonzero__()

    def is_equal_to(self, right, absprec=None):
        right = self.parent().coerce(right)
        return self._element.is_equal_to(right._element, absprec=absprec)

    def precision_absolute(self):
        return self._element.precision_absolute()

    def precision_relative(self):
        return self._element.precision_relative()

    def unit_part(self):
        return self.__class__(self.parent(), self._element.unit_part())

    def valuation(self):
        return self._element.valuation()

    def val_unit(self):
        v, u = self._element.val_unit()
        return v, self.__class__(self.parent(), u)

    def _cache_key(self):
        return self._element._cache_key()

    def __hash__(self):
        raise TypeError("unhashable type: 'sage.rings.padics.extension_element.pAdicGeneralExtensionElement'")

# NotImplementedError: coercion/conversion maps
# NotImplementedError: printing

class ExpansionIterable(object):
    def __init__(self, elt, exp, wrap, mode):
        self._elt = elt
        self._exp = exp
        self._wrap = wrap
        self._mode = mode

    def __iter__(self):
        for x in self._exp:
            yield self._wrap(x)

    def __len__(self):
        return len(self._exp)

    def __getitem__(self, n):
        return self._wrap(self._exp[n])

    def __repr__(self):
        if self._mode == 'simple':
            modestr = ""
        elif self.mode == 'smallest':
            modestr = " (balanced)"
        else:
            modestr = " (teichmuller)"
        p = self._elt.prime_pow.prime
        return "%s-adic expansion of %s%s"%(p, self._elt, modestr)

class pAdicGeneralMorphism(RingHomomorphism):
    """
    A homomorphism from a relative extension L/K of p-adic rings or fields to an arbitrary ring A.

    One can specify such a homomorphism by giving either
    * a homomorphism from K into A together with the image in A of the generator of L/K.
    * a homomorphism from the backend representation to A.

    TODO: When A is non-exact, describe precision.

    INPUT:

        - ``parent`` -- the Homset containing this morphism
        - ``im_gen`` -- the image of the generator of L/K (a length one list is also accepted)
        - ``base_hom`` -- the homomorphism from the base ring (defaults to the coercion map from K to A)
        - ``backend_hom`` -- the homomorphism on the two-step extension representating L/K
        - ``check`` -- whether to check that ``im_gen`` defines a valid homomorphism
    """
    def __init__(self, parent, im_gen=None, base_hom=None, backend_hom=None, check=True):
        # We specify either im_gen or backend_hom but not both
        L = parent.domain()
        K = L.base_ring()
        A = parent.codomain()
        Ax = A['x']
        B = L._backend()
        cat = parent.homset_category()
        if backend_hom is None:
            if base_hom is None:
                base_hom = A.coerce_map_from(K)
                if base_hom is None:
                    raise ValueError("Must specify homomorphism on the base ring")
            elif base_hom.domain() is not K or base_hom.codomain() is not A:
                raise ValueError("Base homomorphism does not have correct domain/codomain")
            if im_gen is None:
                if self.degree() > 1:
                    raise ValueError("Must specify the image of the generator")
            elif check:
                # Should check be done here or in the hom constructions down below?
                if im_gen.parent() is not A:
                    raise ValueError
                pol = Ax([base_hom(c) for c in L.defining_polynomial()])
                if pol(im_gen) != 0:
                    raise ValueError("relations do not all (canonically) map to 0 under map determined by images of generators")
            def get_image(pol):
                if im_gen is None:
                    assert pol.degree() < 1
                    return pol[0]
                else:
                    return pol(im_gen)
            if B.absolute_e() == 1:
                # B is unramified
                if B.absolute_f() == 1:
                    # B is a trivial extension of Qp or Zp
                    # Every hom from Qp or Zp is a natural one
                    backend_hom = B.hom(A)
                else:
                    u = L(B.gen())
                    u_image = get_image(Ax([base_hom(c) for c in u.polynomial()]))
                    backend_hom = B.hom([u_image], check=False)
            else:
                if B.absolute_f() == 1:
                    # B is a totally ramified extension of Qp or Zp
                    pi = L(B.uniformizer())
                    pi_image = get_image(Ax([base_hom(c) for c in pi.polynomial()]))
                    backend_hom = B.hom([pi_image], check=False)
                else:
                    # Two step extension
                    U = B.base_ring()
                    u = L(B(U.gen()))
                    pi = L(B.uniformizer())
                    u_image = get_image(Ax([base_hom(c) for c in u.polynomial()]))
                    pi_image = get_image(Ax([base_hom(c) for c in pi.polynomial()]))
                    backend_hom = B.hom([pi_image], base_hom=U.hom([u_image]), check=False)
        else:
            if im_gen is not None or base_hom is not None:
                raise NotImplementedError("Cannot specify both im_gen/base_hom and backend_hom")
