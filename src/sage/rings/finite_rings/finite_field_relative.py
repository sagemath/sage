# -*- coding: utf-8 -*-
r"""
Relative Extensions of Finite Fields

The classes in this module mainly model extensions of extensions of finite
fields. In fact they are used for extensions that in their implementation rely
on an isomorphic absolute field; this includes e.g. all trivial extensions.

EXAMPLES::

    sage: k.<a> = GF(4)
    sage: R.<x> = k[]
    sage: l.<b> = k.extension(x^2 + x + a, implementation="GF"); l
    Finite Field in b of size 2^4 over its base
    sage: l.degree()
    2
    sage: l.absolute_degree()
    4

AUTHORS:

- Julian Rüth: initial version

"""
#*****************************************************************************
#       Copyright (C) 2019 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

from sage.misc.cachefunc import cached_method
from sage.misc.randstate import seed

from .finite_field_base import FiniteField
from .element_relative import FiniteField_relativeElement
from sage.rings.ring_extension import RingExtensionWithGen
from sage.rings.ring_extension_conversion import backend_parent

# TODO: Make sure that we run TestSuite for all the constellations that I can
# imagine. TestSuite should test something for every method that I touched in
# finite_field_base.pyx.
class FiniteField_relative(FiniteField, RingExtensionWithGen):
    r"""
    A finite field extension which delegates all the computations to an
    absolute finite field, the ``_backend``. This is used for relative
    extensions of (non-prime) finite fields and for trivial extensions of
    finite fields.

    EXAMPLES::

        sage: k.<a> = GF(4)
        sage: R.<x> = k[]
        sage: l.<b> = k.extension(x^2 + x + a, implementation="GF"); l
        Finite Field in b of size 2^4 over its base

    Extensions can also be constructed by just specifying degrees::

        sage: k = GF(4)
        sage: l = k.extension(3, absolute=False); l
        Finite Field in z6 of size 2^6 over its base

    A trivial extension::

        sage: k = GF(2)
        sage: m = k.extension(1, absolute=False); m
        Trivial extension of Finite Field of size 2

    TESTS::

        sage: from sage.rings.finite_rings.finite_field_relative import FiniteField_relative
        sage: isinstance(l, FiniteField_relative)
        True
        sage: isinstance(m, FiniteField_relative)
        True

    """
    def __init__(self, base, modulus, names, category=None, **kwds):
        r"""
        TESTS:

        Run the test suite for all relative extensions mentioned in the doctests of this file::

            sage: k.<a> = GF(4)
            sage: R.<x> = k[]
            sage: l.<b> = k.extension(x^2 + x + a, implementation="GF")

        # TODO: actually do that

        """
        order = base.order() ** modulus.degree()

        from sage.all import GF, FiniteFields, Hom
        if kwds.get('backend') is not None:
            backend = kwds['backend']
        else:
            backend = GF(order, names=["b%s"%(base.absolute_degree() * modulus.degree(),)], **kwds)

        assert modulus.base_ring() is base
        self._modulus = modulus
        assert names and names[0]
        category = category or FiniteFields()

        FiniteField.__init__(self, base, names, normalize=False, category=category)

        with seed(0): # We want the isomorphism with the backend to be deterministic
            defining_embedding = self.base_ring()._any_embedding(backend)
            gen = modulus.map_coefficients(defining_embedding).roots(multiplicities=False)[0]
        from .hom_finite_field import FiniteFieldHomomorphism_generic
        RingExtensionWithGen.__init__(self, defining_morphism=defining_embedding, gen=gen, names=names, import_methods=False, category=category)

        self.register_conversion(self.free_module(map=True)[1])

    # We want to use RingExtensionWithGen's free_module
    #free_module = RingExtensionWithGen.free_module

    def __reduce__(self):
        r"""
        TESTS::

            sage: k = GF(4).extension(2, absolute=False)
            sage: loads(dumps(k)) is k
            True

        """
        return self._factory_data[0].reduce_data(self)

    def absolute_gen(self):
        return self(self._backend.absolute_gen())

    def absolute_field(self, map=False, names=None):
        r"""
        Return an absolute extension of the prime field isomorphic to this field.

        INPUT:

        - ``map`` -- boolean, default ``False``), whether to return maps to and from the absolute field
        - ``names`` -- string, the variable name for the absolute field

        OUTPUT:

        If ``map`` is ``False``, an absolute finite field isomorphic to this one.  Otherwise,

        - ``absolute`` -- the absolute field
        - ``from_absolute`` -- an isomorphism from the absolute field to this field
        - ``to_absolute`` -- the inverse isomorphism

        EXAMPLES::

            sage: k = GF(9).extension(3, absolute=False)
            sage: k.absolute_field()
            Finite Field in b6 of size 3^6

            sage: l, f, g = k.absolute_field(map=True)
            sage: a = k.random_element(); f(g(a)) == a
            True
            sage: b = l.random_element(); g(f(b)) == b
            True
            sage: f
            Ring morphism:
              From: Finite Field in b6 of size 3^6
              To:   Finite Field in z6 of size 3^6 over its base
              Defn: b6 |--> 1 + (2*z2 + 1)*z6 + (z2 + 1)*z6^2
            sage: g
            Canonical morphism:
              From: Finite Field in z6 of size 3^6 over its base
              To:   Finite Field in b6 of size 3^6
        """
        backend, from_backend, to_backend = backend_parent(self, map=True)
        absolute = backend.absolute_field(map=map, names=names)
        if map:
            (absolute, absolute_to_backend, backend_to_absolute) = absolute
            return (absolute,
                from_backend * absolute_to_backend,
                backend_to_absolute * to_backend)
        else:
            return absolute

    def _compatible_family(self):
        backend = self._backend
        fam = backend._compatible_family()
        f = self.convert_map_from(backend)
        return {d: (f(b), bpoly) for (d, (b, bpoly)) in fam.items()}

    def characteristic(self):
        r"""
        Return `p`, the characteristic of this field.

        EXAMPLES::

            sage: k = GF(9).extension(3, absolute=False)
            sage: k.characteristic()
            3

        """
        return self._backend.characteristic()

    def _repr_topring(self, **options):
        """
        Return a printable representation of this field.

        EXAMPLES::

            sage: k = GF(9).extension(3, absolute=False)
            sage: k # indirect doctest
            Finite Field in z6 of size 3^6 over its base
        """
        return "Finite Field in %s of size %s^%s"%(self.variable_name(), self.characteristic(), self.absolute_degree())

    _repr_ = RingExtensionWithGen._repr_

    Element = FiniteField_relativeElement
