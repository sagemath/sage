r"""
Extension of rings

See :func:`RingExtension` for documentation.

AUTHOR:

- Xavier Caruso (2016)
"""

#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************


from sage.rings.integer_ring import IntegerRing
from sage.rings.ring_extension import RingExtension_class


def RingExtension(ring, base=None, defining_morphism=None):
    r"""
    Create a ring extension

    INPUT:

    - ``ring`` -- a commutative ring

    - ``base`` (optional) -- a commutative ring lying below ``ring``

    - ``defining_morphism`` (optional) -- an homomorphism of rings
    from ``base`` to ``ring``

    OUTPUT:

    The ring ``ring`` viewed as an algebra over ``base`` through
    the homomorphism ``defining_morphism``.

    If ``base`` is not set the canonical base of ``ring`` (which
    is available via ``ring.base()``) is used.

    If ``defining_morphism`` is not set and ``base`` coerces to
    ``ring`` the coercion map is used for defining the extension.


    CREATION OF EXTENSIONS:

    We create an extension of finite fields::

        sage: K = GF(5^2); z2 = K.gen()
        sage: L = GF(5^4); z4 = L.gen()
        sage: E = RingExtension(L,K); E
        Finite Field in z4 of size 5^4 viewed as an algebra over its base

    The base of ``E`` is accessible as follows::

        sage: E.base_ring()
        Finite Field in z2 of size 5^2
        sage: E.base_ring() is K
        True

    Here is an example where the base is implicit::

        sage: R.<x> = PolynomialRing(K)
        sage: RingExtension(R)
        Univariate Polynomial Ring in x over Finite Field in z2 of size 5^2 viewed as an algebra over its base


    ELEMENTS IN EXTENSIONS:

    We can create elements in the extension `E` using standard methods::

        sage: E.zero()
        0
        sage: E.one()
        1
        sage: E.an_element()
        0
        sage: E.random_element()  # random
        4*z4^2 + 4*z4 + 2

    Conversion also works::

        sage: a = z4 + 1
        sage: a.parent()
        Finite Field in z4 of size 5^4
        sage: aE = E(a); aE
        z4 + 1
        sage: aE.parent()
        Finite Field in z4 of size 5^4 viewed as an algebra over its base
        sage: aE.parent() is E
        True

    Of course all ring operations are available in E::

        sage: bE = (aE + 1) * (aE + 2); bE
        z4^2 + 1

    And the result stays in the extension::

        sage: bE.parent()
        Finite Field in z4 of size 5^4 viewed as an algebra over its base
        sage: bE.parent() is E
        True


    COERCION:

    A coercion map going from the base to the extension is set::

        sage: E.has_coerce_map_from(K)
        True
        sage: E.coerce_map_from(K)
        Ring morphism:
          From: Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4 viewed as an algebra over its base

    Therefore we can safely add an element of the base with an element
    of the extension, obtaining this way an element in the extension::

        sage: c = z2 + 3
        sage: s = aE + c; s
        z4^3 + z4^2 + 2*z4 + 2
        sage: s.parent()
        Finite Field in z4 of size 5^4 viewed as an algebra over its base
        sage: s.parent() is E
        True

    When the defining morphism of the extension ``E = L/K`` is a coercion
    map (see below), a coercion from ``E`` to ``L`` (acting as the identity)
    is set::

        sage: L.has_coerce_map_from(E)
        True
        sage: L.coerce_map_from(E)
        Ring morphism:
          From: Finite Field in z4 of size 5^4 viewed as an algebra over its base
          To:   Finite Field in z4 of size 5^4


    Coercions between extensions are implemented as follows: an extension
    L1/K1 coerces to another extension L2/K2 when L1 coerces to L2 and K2
    coerces to K1 (be careful at the direction) and the defining morphisms
    agree. For example::

        sage: K1 = GF(3^4); L1 = GF(3^8);  E1 = RingExtension(L1,K1)
        sage: K2 = GF(3^2); L2 = GF(3^16); E2 = RingExtension(L2,K2)
        sage: E2.has_coerce_map_from(E1)
        True

        sage: x1 = E1.random_element()
        sage: x2 = E2.random_element()
        sage: s = x1 + x2
        sage: s.parent()
        Finite Field in z16 of size 3^16 viewed as an algebra over its base
        sage: s.parent() is E2
        True


    DEFINING MORPHISM:

    When creating the extension ``E = L/K``, we have not specified the defining
    morphism. In that case, the coercion map from the base to the ring (if it
    exists) is used by default::

        sage: E.defining_morphism()
        Ring morphism:
          From: Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4
          Defn: z2 |--> z4^3 + z4^2 + z4 + 3
        sage: E.defining_morphism() == L.coerce_map_from(K)
        True

    When there is no coercion map from the base to the ring, an error is
    raised::

        sage: L2 = GF(7^2)
        sage: RingExtension(L2,K)
        Traceback (most recent call last):
        ...
        ValueError: No coercion map from Finite Field in z2 of size 5^2 to Finite Field in z2 of size 7^2

    It is also possible to use a custom defining morphism.
    As an example, let us create the algebra L/K viewed through the
    Frobenius endomorphism.
    We first create the Frobenius map phi : K -> L::

        sage: FrobK = K.frobenius_endomorphism()
        sage: phi = FrobK.post_compose(L.coerce_map_from(K))
        sage: phi
        Composite map:
          From: Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4
          Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                then
                  Ring morphism:
                  From: Finite Field in z2 of size 5^2
                  To:   Finite Field in z4 of size 5^4
                  Defn: z2 |--> z4^3 + z4^2 + z4 + 3

    And then the twisted extension::

        sage: ExtFrob = RingExtension(L, K, phi); ExtFrob
        Finite Field in z4 of size 5^4 viewed as an algebra over its base

    Explicitely composing with coercion maps (on the left and on the right)
    is not mandatory: Sage does it automatically for you if necessary.
    Consequenty the extension ``ExtFrob`` can be created in a simpler way
    as follows::

        sage: ExtFrob2 = RingExtension(L, K, FrobK)
        sage: ExtFrob is ExtFrob2
        True

    We draw the attention of the user that non canonical defining morphisms
    should be used with extreme care!
    As an example, in the example below, Sage defines a coercion map from
    ``K`` to ``ExtFrob`` but no coercion map from ``ExtFrob`` to ``L`` (since
    otherwise the Frobenius morphism from ``K`` to ``L`` would be promoted as
    a coercion map as well).

        sage: ExtFrob.has_coerce_map_from(K)
        True
        sage: ExtFrob.coerce_map_from(K)
        Ring morphism:
          From: Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4 viewed as an algebra over its base

        sage: L.has_coerce_map_from(ExtFrob)
        False

        sage: aE = ExtFrob(a)
        sage: aE + c
        4*z4^3 + 4*z4^2 + 2
        sage: a + c^5
        4*z4^3 + 4*z4^2 + 2


    AUTHOR:

    - Xavier Caruso (2016)
    """
    if defining_morphism is None:
        coerce = True
        if base is None:
            base = ring.base_ring()
            if base is None:
                # ZZ is initial in the category of rings
                base = IntegerRing()
        if ring.has_coerce_map_from(base):
            defining_morphism = ring.coerce_map_from(base)
        else:
            raise ValueError("No coercion map from %s to %s" % (base,ring))
    else:
        domain = defining_morphism.domain()
        codomain = defining_morphism.codomain()
        coerce = (defining_morphism == codomain.coerce_map_from(domain))
        if base is None:
            base = domain
        elif base is domain:
            pass
        elif domain.has_coerce_map_from(base):
            coercion = domain.coerce_map_from(base)
            defining_morphism = defining_morphism.pre_compose(coercion)
        else:
            raise ValueError("No coercion map from %s to %s" % (base,domain))
        if ring is codomain:
            pass
        elif ring.has_coerce_map_from(codomain):
            coercion = ring.coerce_map_from(codomain)
            defining_morphism = defining_morphism.post_compose(coercion)
        else:
            raise ValueError("No coercion map from %s to %s" % (codomain,ring))
    return RingExtension_class(defining_morphism, coerce)


#def TowerExtensions(*rings):
#    tower = rings[-1]
#    for i in range(len(rings)-2, -1, -1):
#        tower = RingExtension(rings[i], tower)
#    return tower
