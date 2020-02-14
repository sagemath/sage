# -*- coding: utf-8 -*-
r"""
Relative Extensions of Finite Fields

The classes in this module mainly model extensions of extensions of finite
fields. In fact they are used for extensions that in their implementation rely
on an isomorphic absolute field; this includes e.g. all trivial extensions.

EXAMPLES::

    sage: k.<a> = GF(4)
    sage: R.<x> = k[]
    sage: l.<b> = k.extension(x^2 + x + a, absolute=False); l
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

from .finite_field_base import FiniteField
from .element_relative import FiniteField_relativeElement
from sage.rings.ring_extension import RingExtensionWithGen

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
        sage: l.<b> = k.extension(x^2 + x + a, absolute=False); l
        Finite Field in b of size 2^4 over its base

    Extensions can also be constructed by just specifying degrees::

        sage: k = GF(4)
        sage: l = k.extension(3, absolute=False); l
        Finite Field in z6 of size 2^6 over its base

    A trivial extension::

        sage: k = GF(2)
        sage: m = k.extension(1, absolute=False); m
        Finite Field in z1 of size 2 over its base

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
            sage: l.<b> = k.extension(x^2 + x + a, absolute=False)
        
        # TODO: actually do that

        """
        order = base.order() ** modulus.degree()

        from sage.all import GF, FiniteFields, Hom
        backend = GF(order, names=["b%s"%(base.absolute_degree() * modulus.degree(),)], **kwds)

        assert modulus.base_ring() is base
        self._modulus = modulus
        assert names and names[0]

        FiniteField.__init__(self, base, names, normalize=False, category=category or FiniteFields())

        defining_embedding = self.base_ring()._any_embedding(backend)
        gen = modulus.map_coefficients(defining_embedding).any_root()
        RingExtensionWithGen.__init__(self, defining_morphism=defining_embedding, gen=gen, names=names)

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

    def absolute_field(self, map=False, **kwds):
        r"""
        Return an absolute extension of the prime field isomorphic to this field.

        EXAMPLES::

            sage: k = GF(9).extension(3, absolute=False)
            sage: k.absolute_field()
            Finite Field in z6 of size 3^6

            sage: k.absolute_field(map=True)

        """
        backend = self._backend
        absolute = backend.absolute_field(map=map, **kwds)
        if map:
            (absolute, absolute_to_backend, backend_to_absolute) = absolute
            return (absolute,
                self.convert_map_from(backend) * absolute_to_backend,
                backend_to_absolute * backend.convert_map_from(self))
        else:
            return absolute

    def characteristic(self):
        r"""
        Return `p`, the characteristic of this field.

        EXAMPLES::

            sage: k = GF(9).extension(3, absolute=False)
            sage: k.characteristic()
            3

        """
        return self._backend.characteristic()

    def _repr_(self):
        """
        Return a printable representation of this field.

        EXAMPLES::

            sage: k = GF(9).extension(3, absolute=False)
            sage: k # indirect doctest
        """
        return "Finite Field in %s of size %s^%s over its base"%(self.variable_name(), self.characteristic(), self.absolute_degree())

    Element = FiniteField_relativeElement
