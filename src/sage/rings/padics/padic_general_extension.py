r"""
General extensions of p-adic rings and fields; the base ring may also be an
extension.

These are implemented as proxy parents, backed by an absolute extension.

EXAMPLES:

A trivial extension::

    sage: L.<a> = Qp(2).extension(x)
    sage: L
    2-adic Field with capped relative precision 20
    sage: L is Qp(2)
    False
    sage: a == 0
    True

A trivial extension of a trivial extension::

    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b - a)
    sage: M
    2-adic Field with capped relative precision 20
    sage: b == a
    True

An unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: L
    Field in a with defining polynomial (1 + O(2^20))*x^2 + (2 + O(2^21))*x + 2^2 + O(2^22) over its base
    sage: a^2 + 2*a + 4 == 0
    True
    sage: L.f()
    2

"""
# ****************************************************************************
#       Copyright (C)      2019 David Roe <roed.math@gmail.com>
#                     2019-2022 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from .padic_general_extension_element import pAdicGeneralExtensionElement
from .padic_extension_generic import pAdicExtensionGeneric
from sage.rings.ring_extension import RingExtensionWithGen


class pAdicGeneralExtension(RingExtensionWithGen, pAdicExtensionGeneric):
    def __init__(self, exact_modulus, poly, prec, print_mode, shift_seed, names, implementation='FLINT', category=None):
        base = poly.base_ring()
        self._exact_modulus = exact_modulus
        self._shift_seed = shift_seed
        self._implementation = 'proxy'
        # TODO: To make things work for now, we use the base's prime pow.
        self.prime_pow = base.prime_pow
        self._prec_type = base._prec_type
        category = category or base.category()

        pAdicExtensionGeneric.__init__(self, exact_modulus, poly, prec, print_mode, names, pAdicGeneralExtensionElement, category=category)

        if not self._exact_modulus.is_monic():
            raise NotImplementedError(f"defining modulus must be monic but {exact_modulus} is not")

        # Construct the backend, a p-adic ring that is not a general extension.
        if self.f() == 1 and self.e() == 1:
            # This is a trivial extension. The best backend is base ring
            # (possibly rewritten as an absolute extension.)
            assert self._exact_modulus.degree() == 1

            (backend, backend_to_base, base_to_backend) = base.absolute_ring(map=True)
            defining_morphism = base_to_backend
            gen = defining_morphism(base(-self._exact_modulus[0]))
        else:
            # The underlying Zp or Qp
            backend_base = self.ground_ring_of_tower()

            # The unramified part of this extension.
            if self.absolute_f() == 1:
                backend_unramified = backend_base
            else:
                from sage.all import Zq, Qq
                backend_unramified = self.ground_ring_of_tower().change(q=self.prime()**self.absolute_f(), names=names[2])

            # The totally ramified part of this extension.
            if self.absolute_e() == 1:
                backend = backend_unramified
            else:
                raise NotImplementedError("cannot construct general ramified extensions yet")

            # TODO: This won't work in general. When it works it should be correct.
            defining_morphism = self.base_ring().hom(backend)

            # TODO: The poly.change_ring() might not have enough precision.
            # TODO: The any_root() might not have enough precision.
            gen = poly.change_ring(defining_morphism).any_root()

        if backend is not backend.absolute_ring():
            raise NotImplementedError("relative backends are not supported for general p-adic extensions yet")

        self._backend = backend

        RingExtensionWithGen.__init__(self, defining_morphism=defining_morphism, gen=gen, names=[self.variable_name()], category=category)

    modulus = pAdicExtensionGeneric.defining_polynomial

    @cached_method
    def f(self):
        r"""
        Return the residual degree of this ring over its base ring.

        EXAMPLES:

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x - 2)
            sage: L.f()
            1

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.f()
            2

        """
        v = self.base_ring().exact_valuation()
        w = v.mac_lane_approximants(self._exact_modulus, assume_squarefree=True, require_final_EF=True)
        if len(w) != 1:
            raise ValueError("defining polynomial is not irreducible")
        return w[0].F()

    def absolute_f(self):
        r"""
        Return the absolute residue degree of this ring, i.e., the degree of
        the residue field over its prime subfield.

        EXAMPLES:

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x - 2)
            sage: L.absolute_f()
            1

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.absolute_f()
            2

        """
        return self.f() * self.base_ring().absolute_f()

    def e(self):
        r"""
        Return the ramification degree of this ring over its base ring.

        EXAMPLES::

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.e()
            1

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.e()
            1

        """
        return self._exact_modulus.degree() // self.f()

    def absolute_e(self):
        r"""
        Return the total degree of ramification of this ring.

        EXAMPLES::

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.absolute_e()
            1

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.absolute_e()
            1

        """
        return self.e() * self.base_ring().absolute_e()

    def absolute_ring(self, map=False, **kwds):
        r"""
        Return an absolute extension of the absolute base isomorphic to this
        field.

        Note that this might not be a simple extension. It might be a p-adic
        base ring for a trivial extension or a two step extension, i.e., a
        totally ramified extension given by an Eisenstein polynomial over an
        unramified extension.

        EXAMPLES:

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.absolute_ring()
            2-adic Field with capped relative precision 20

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.absolute_ring()
            2-adic Unramified Extension Field in a_u defined by x^2 + x + 1
            sage: M, M_to_L, L_to_M = L.absolute_ring(map=True)
            sage: M_to_L(L_to_M(L.gen())) == L.gen()
            True

        """
        if map:
            return self._backend, self.convert_map_from(self._backend), self._backend.convert_map_from(self)
        else:
            return self._backend

    def teichmuller(self, x, prec=None):
        R = self._backend()
        x = R(x) if prec is None else R(x, prec)
        return self(R.teichmuller(x))

    def _prec_type(self):
        return self._backend()._prec_type()

    def is_field(self):
        r"""
        Return whether this ring is a field.

        EXAMPLES::

            sage: L.<a> = Zp(2).extension(x + 2)
            sage: L.is_field()
            False

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.is_field()
            True

        """
        return self._backend.is_field()

    def random_element(self, **kwds):
        return self(self._backend().random_element(**kwds))

    def residue_ring(self, n):
        raise NotImplementedError

    def residue_class_field(self):
        r"""
        Return the residue class field of this ring.

        EXAMPLES::

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.residue_class_field()
            Finite Field of size 2

        A trivial extension of a trivial extension::

            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b - a)
            sage: M.residue_field()

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.residue_class_field()

        """
        return self.base_ring().residue_class_field().extension(self.f(), absolute=False)

    def inertia_subring(self):
        raise NotImplementedError

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
