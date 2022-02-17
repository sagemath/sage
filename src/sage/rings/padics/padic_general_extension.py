r"""
General extensions of p-adic rings and fields; the base ring may also be an
extension.

These are implemented as proxy parents, backed by an absolute extension.

EXAMPLES:

A trivial extension::

    sage: L.<a> = Qp(2).extension(x)
    sage: L
    2-adic Trivial Extension Field in a defined by x
    sage: a == 0
    True

A trivial extension of a trivial extension::

    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b - a)
    sage: M
    2-adic Trivial Extension Field in b defined by b
    sage: b == a
    True

An unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: L
    2-adic Unramified Extension Field in a defined by x^2 + 2*x + 4
    sage: a^2 + 2*a + 4 == 0
    True
    sage: L.f()
    2

An unramified extension given by a non-monic defining polynomial (currently, not supported, see #33362)::

    sage: L.<a> = Qp(2).extension(4*x^2 + 2*x + 1)
    Traceback (most recent call last):
    ...
    ValueError: G must be integral
    sage: a^2 + 2*a + 4 == 0  # not tested
    True
    sage: L.f()  # not tested
    2

An unramified extension given by a non-integral defining polynomial (currently, not supported, see #33362)::

    sage: L.<a> = Qp(2).extension(x^2 + x/4 + 1/16)
    Traceback (most recent call last):
    ...
    ValueError: G must be integral
    sage: a^2 + 2*a + 4 == 0  # not tested
    True
    sage: L.f()  # not tested
    2

A trivial extension of an unramified extension::

    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b - 2)
    sage: M
    2-adic Trivial Extension Field in b defined by b - 2 over its base field
    sage: M.f()
    1
    sage: M.absolute_f()
    2

An unramified extension of a trivial extension::

    sage: L.<a> = Qp(2).extension(x + 1)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 - b - a)
    sage: M
    2-adic Unramified Extension Field in b defined by b^2 - b + 1
    sage: M.f()
    2

An unramified extension of an unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + x + 1)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 + b + a)
    sage: M
    2-adic Unramified Extension Field in b defined by b^2 + b + a over its base field
    sage: M.f()
    2
    sage: M.absolute_f()
    4

::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 + a*b + 4)
    sage: M
    2-adic Unramified Extension Field in b defined by b^2 + a*b + 4 over its base field
    sage: M.f()
    2
    sage: M.absolute_f()
    4

A totally ramified extension not given by an Eisenstein polynomial::

    sage: L.<a> = Qp(2).extension(x^2 + 8)
    sage: L.e()
    2

A trivial extension of a totally ramified extension::

    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b - a)

A totally ramified extension of a trivial extension::

    sage: L.<a> = Qp(2).extension(x - 2)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 - a)
    sage: M.absolute_e(), M.absolute_f()
    (2, 1)

A totally ramified extension of an unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 - 8)
    sage: M.absolute_e(), M.absolute_f()
    (2, 2)

An unramified extension of a totally ramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 8)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 + a*b + a^2)
    sage: M.absolute_e(), M.absolute_f()
    (2, 2)

A totally ramified extension of a totally ramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 8)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 + 2*a)
    sage: M.absolute_e(), M.absolute_f()
    (4, 1)

A mixed case::

    sage: L.<a> = Qp(2).extension(x^4 + 8*x^2 + 64)
    sage: L.absolute_e(), L.absolute_f()
    (2, 2)

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
from sage.rings.ring_extension_conversion import backend_parent
from sage.rings.padics.pow_computer import PowComputer_class


class pAdicGeneralExtension(pAdicExtensionGeneric):
    r"""
    A general extension of a p-adic ring such as a relative extension or an
    extension not given by an unramified polynomial or an Eisenstein
    polynomial.

    EXAMPLES:

        sage: L.<a> = Qp(2).extension(x)

    """
    def modulus(self, *args, **kwds):
        return self.defining_polynomial(*args, **kwds)

    @staticmethod
    def _hom_to_backend(base, backend):
        if base is base.ground_ring_of_tower():
            return base.hom(backend)
        base_map = pAdicGeneralExtension._hom_to_backend(base.base_ring(), backend)

        # TODO: The poly.change_ring() might not have enough precision.
        # TODO: The any_root() might not have enough precision.
        modulus = base.modulus().change_ring(base_map)
        return base.hom([modulus.any_root()], codomain=backend, base_map=base_map)

    def _create_backend(self):
        r"""
        Return a backend for this extension, i.e., a p-adic ring that is not a
        general extension itself.
        """
        if self.f() == 1 and self.e() == 1:
            # This is a trivial extension. The best backend is base ring
            # (possibly rewritten as an absolute extension.)
            assert self._exact_modulus.degree() == 1

            (backend, backend_to_base, base_to_backend) = self.base_ring().absolute_ring(map=True)
            defining_morphism = base_to_backend
            gen = defining_morphism(self.base_ring()(-self._exact_modulus[0]))
        else:
            # The underlying Zp or Qp
            backend_base = self.ground_ring_of_tower()

            # The unramified part of this extension.
            if self.absolute_f() == 1:
                backend_unramified = backend_base
            else:
                backend_unramified = self.ground_ring_of_tower().change(q=self.prime()**self.absolute_f(), names=self._printer.unram_name)

            # The totally ramified part of this extension.
            if self.absolute_e() == 1:
                backend = backend_unramified
            else:
                # We construct the charpoly of the uniformizer and factor it
                # over the unramified part. Currently, we do this completely
                # naively in the corresponding number field which is terribly
                # slow.
                charpoly = self.exact_field().absolute_field('x').valuation(self.prime()).uniformizer().charpoly()

                assert charpoly.degree() == self.absolute_e() * self.absolute_f()

                charpoly = charpoly.change_ring(backend_unramified.exact_field())
                charpoly = backend_unramified.exact_field().valuation(self.prime()).montes_factorization(charpoly)

                assert all(f.degree() == self.absolute_e() for f,e in charpoly), f"charpoly of uniformizer did not factor as an approximate {self.absolute_f()}th power: {charpoly}"

                minpoly = charpoly[0][0]

                backend = backend_unramified.extension(minpoly, names='pi')

            defining_morphism = pAdicGeneralExtension._hom_to_backend(self.base_ring(), backend)

            # TODO: The poly.change_ring() might not have enough precision.
            # TODO: The any_root() might not have enough precision.
            gen = self.defining_polynomial().change_ring(defining_morphism).any_root()

        if backend is not backend.absolute_ring():
            raise NotImplementedError("relative backends are not supported for general p-adic extensions yet")

        return defining_morphism, gen

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
        return self.exact_valuation().F()

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
        return self.exact_valuation().E()

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

    def absolute_ring(self, map=False):
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
        return backend_parent(self, map=map)

    def teichmuller(self, x, prec=None):
        R = self._backend
        x = R(x) if prec is None else R(x, prec)
        return self(R.teichmuller(x))

    def _prec_type(self):
        return self._backend._prec_type()

    def is_field(self):
        r"""
        Return whether this ring is a field.

        EXAMPLES::

            sage: L.<a> = Zp(2).extension(x + 3)
            sage: L.is_field()
            False

            sage: L.<a> = Qp(2).extension(x)
            sage: L.is_field()
            True

        """
        return self._backend.is_field()

    def random_element(self, **kwds):
        return self(self._backend.random_element(**kwds))

    def residue_ring(self, n):
        raise NotImplementedError

    @cached_method
    def residue_class_field(self):
        r"""
        Return the residue class field of this ring.

        EXAMPLES::

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.residue_class_field()
            Trivial extension of Finite Field of size 2

        A trivial extension of a trivial extension::

            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b - a)
            sage: M.residue_field()
            Trivial extension of Trivial extension of Finite Field of size 2

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.residue_class_field()
            Finite Field in z2 of size 2^2 over its base

        A trivial extension of an unramified extension::

            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b - 2)
            sage: M.residue_field()
            Trivial extension of Finite Field in z2 of size 2^2 over its base

        An unramified extension of an unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b^2 + a*b + 4)
            sage: m = M.residue_field()
            sage: m
            Finite Field in z4 of size 2^4 over its base
            sage: m.base_ring() is L.residue_field()
            True
            sage: m.modulus()
            x^2 + x + z2

        """
        return self.base_ring().residue_class_field().extension(self.f(), absolute=False, implementation="GF", backend=self._backend.residue_class_field())

    def inertia_subring(self):
        if self.absolute_e() == 1:
            return self
        if self.f() == 1:
            return self.base_ring().inertia_subring()
        raise NotImplementedError("cannot compute inertia subring of this tower yet")

    def uniformizer(self):
        backend, from_backend, _ = backend_parent(self, map=True)
        return from_backend(backend.uniformizer())

    def uniformizer_pow(self, n):
        backend, from_backend, _ = backend_parent(self, map=True)
        return from_backend(backend.uniformizer_pow(n))

    def _uniformizer_print(self):
        return self._backend._uniformizer_print()

    def gen_unram(self):
        backend, from_backend, _ = backend_parent(self, map=True)
        return from_backend(backend.gen_unram())

    def _unram_print(self):
        return self._backend._unram_print()

    def has_pth_root(self):
        return self._backend.has_pth_root()

    def has_root_of_unity(self, n):
        return self._backend.has_root_of_unity(self, n)

    def integer_ring(self):
        raise NotImplementedError

    def construction(self, forbid_frac_field=None):
        # Prefer AlgebraicExtensionFunctor for pushout since FractionField
        # functor often does not work because there is no integer_ring.
        if forbid_frac_field is None:
            forbid_frac_field = True

        # TODO: Change prec of AlgebraicExtensionFunctor
        construction = pAdicExtensionGeneric.construction(self, forbid_frac_field=forbid_frac_field)
        return construction


class pAdicGeneralExtension_ring(pAdicGeneralExtension):
    pass


class pAdicGeneralExtension_field(pAdicGeneralExtension, RingExtensionWithGen):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation='FLINT', category=None):
        r"""

        TESTS::

            sage: L.<a> = Qp(2).extension(x)
            sage: from sage.rings.padics.padic_general_extension import pAdicGeneralExtension
            sage: isinstance(L, pAdicGeneralExtension)
            True
            sage: TestSuite(L).run()

        ::

        """
        base = approx_modulus.base_ring()
        self._exact_modulus = exact_modulus
        self._shift_seed = shift_seed
        self._implementation = 'proxy'
        self._prec_type = base._prec_type
        self.prime_pow = PowComputer_general(base.prime(), cache_limit=0, prec_cap=prec, ram_prec_cap=prec, in_field=base.is_field(), poly=approx_modulus)
        category = category or base.category()

        pAdicGeneralExtension.__init__(self, exact_modulus, approx_modulus, prec, print_mode, names, pAdicGeneralExtensionElement, category=category)

        if prec != self.base_ring().precision_cap():
            raise NotImplementedError("cannot change precision in general extension yet")

        if not self._exact_modulus.is_monic():
            raise NotImplementedError(f"defining modulus must be monic but {exact_modulus} is not")

        if not self._exact_modulus.is_squarefree():
            # We only check squarefreeness here. Irreducibility is checked
            # automatically, when the extensions of the valuations on base to
            # the ring are constructed. (If there is more than one extension,
            # i.e., the polynomial is not irreducible, exact_valuation() is
            # going to complain.)
            raise ValueError("polynomial must be irreducible but %r is not"%(polynomial,))

        defining_morphism, gen = self._create_backend()

        self._backend = gen.parent()
        self._prec = prec * self.e()

        RingExtensionWithGen.__init__(self, defining_morphism=defining_morphism, gen=gen, names=[self.variable_name()], category=category, import_methods=False)

    gen = RingExtensionWithGen.gen
    gens = RingExtensionWithGen.gens
    degree = RingExtensionWithGen.degree
    # Use the implementation of __reduce__ from the factory and ignore RingExtensionWithGen's override.
    __reduce__ = pAdicGeneralExtension.__reduce__


class PowComputer_general(PowComputer_class):
    pass
