#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#                  2022 Julian Rüth <julian.rueth@fsfe.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
#****************************************************************************


from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.categories.pushout import construction_tower
from sage.categories.map cimport Map, FormalCompositeMap
from sage.categories.morphism import IdentityMorphism
from sage.rings.ring_extension_element cimport RingExtensionElement


# For parents
#############

cpdef backend_parent(R, map=False):
    r"""
    Return the backend parent of ``R``.

    INPUT:

    - ``R`` -- a parent

    EXAMPLES:

        sage: from sage.rings.ring_extension_conversion import backend_parent

        sage: K.<a> = GF(5^2).over()  # over GF(5)
        sage: backend_parent(K)
        Finite Field in z2 of size 5^2
        sage: backend_parent(K) is GF(5^2)
        True
    """
    from sage.rings.ring_extension import RingExtension_generic
    if isinstance(R, RingExtension_generic):
        if map:
            return R._backend, R._from_backend_morphism, R._to_backend_morphism
        else:
            return R._backend
    else:
        if map:
            return R, IdentityMorphism(R), IdentityMorphism(R)
        else:
            return R


cpdef from_backend_parent(R, E):
    r"""
    Try to reconstruct a ring extension (somehow related to ``E``)
    whose backend is ``R``.

    INPUT:

    - ``R`` -- a parent

    - ``E`` -- a ring extension

    EXAMPLES::

        sage: from sage.rings.ring_extension_conversion import from_backend_parent

        sage: K.<a> = GF(5^2).over()   # over GF(5)
        sage: L.<b> = GF(5^4).over(K)

        sage: from_backend_parent(GF(5^2), K)
        Field in a with defining polynomial x^2 + 4*x + 2 over its base
        sage: from_backend_parent(GF(5^2), K) is K
        True

    Bases are recognized::

        sage: from_backend_parent(GF(5^2), L)
        Field in a with defining polynomial x^2 + 4*x + 2 over its base
        sage: from_backend_parent(GF(5^2), L) is K
        True

    And also certain constructions::

        sage: S.<x> = GF(5^2)[]
        sage: T = from_backend_parent(S, L)
        sage: T
        Univariate Polynomial Ring in x over Field in a with defining polynomial x^2 + 4*x + 2 over its base
        sage: T.base_ring() is K
        True
    """
    tower = construction_tower(R)
    bases_tower = [ parent for (_, parent) in tower ]
    for base in E.bases():
        backend = backend_parent(base)
        try:
            s = bases_tower.index(backend)
            ans = base
            for i in range(s, 0, -1):
                functor = tower[i][0]
                ans = functor(ans)
            return ans
        except ValueError:
            pass
    return R


# For elements
##############

cpdef backend_element(x):
    r"""
    Return the backend element of ``x``.

    INPUT:

    - ``x`` -- an element

    EXAMPLES:

        sage: from sage.rings.ring_extension_conversion import backend_element

        sage: K.<a> = GF(5^2).over()   # over GF(5)
        sage: backend_element(a)
        z2
        sage: backend_element(a) in GF(5^2)
        True
    """
    if isinstance(x, RingExtensionElement):
        return (<RingExtensionElement>x)._backend
    else:
        return x

cpdef from_backend_element(x, E):
    r"""
    Try to reconstruct an element in a ring extension (somehow
    related to ``E``) whose backend is ``x``.

    INPUT:

    - ``x`` -- an element

    - ``E`` -- a ring extension

    EXAMPLES::

        sage: from sage.rings.ring_extension_conversion import from_backend_element

        sage: K.<a> = GF(5^2).over()   # over GF(5)
        sage: L.<b> = GF(5^4).over(K)
        sage: z2 = GF(5^2).gen()

        sage: from_backend_element(z2, K)
        a
        sage: from_backend_element(z2, K).parent() is K
        True

    Bases are recognized::

        sage: from_backend_element(z2, L)
        a
        sage: from_backend_element(z2, L).parent() is K
        True

    And also certain constructions::

        sage: S.<x> = GF(5^2)[]
        sage: u = from_backend_element(x + z2, L); u
        x + a
        sage: u.parent()
        Univariate Polynomial Ring in x over Field in a with defining polynomial x^2 + 4*x + 2 over its base
        sage: u.base_ring() is K
        True
    """
    parent = from_backend_parent(x.parent(), E)

    if parent is None:
        return x

    if x.parent() is parent:
        return x

    if parent is backend_parent(parent):
        return parent(x)

    _, from_parent_backend, _ = backend_parent(parent, map=True)
    return from_parent_backend(x)


# For morphisms
###############

cpdef from_backend_morphism(f, E):
    r"""
    Try to reconstruct a morphism between ring extensions
    (somehow related to ``E``) whose backend is ``f``.

    INPUT:

    - ``x`` -- a morphism

    - ``E`` -- a ring extension

    EXAMPLES::

        sage: from sage.rings.ring_extension_conversion import from_backend_morphism

        sage: K.<a> = GF(5^2).over()  # over GF(5)
        sage: L.<b> = GF(5^6).over(K)

        sage: Frob = GF(5^2).frobenius_endomorphism()
        sage: from_backend_morphism(Frob, K)
        Ring endomorphism of Field in a with defining polynomial x^2 + 4*x + 2 over its base
          Defn: a |--> 1 - a

    Bases are recognized::

        sage: from_backend_morphism(Frob, L)
        Ring endomorphism of Field in a with defining polynomial x^2 + 4*x + 2 over its base
          Defn: a |--> 1 - a
    """
    from sage.rings.ring_extension_morphism import RingExtensionHomomorphism
    cdef domain = from_backend_parent(f.domain(), E)
    cdef codomain = from_backend_parent(f.codomain(), E)
    return RingExtensionHomomorphism(domain.Hom(codomain), f)


# Generic
#########

cpdef to_backend(arg):
    r"""
    Return the backend of ``arg``.

    EXAMPLES::

        sage: from sage.rings.ring_extension_conversion import to_backend

    This function accepts parents::

        sage: K.<a> = GF(5^2).over()   # over GF(5)
        sage: to_backend(K)
        Finite Field in z2 of size 5^2

    elements::

        sage: to_backend(a)
        z2

    morphisms::

        sage: f = K.hom([a^5])
        sage: to_backend(f)
        Ring morphism:
          From: Finite Field in z2 of size 5^2
          To:   Field in a with defining polynomial x^2 + 4*x + 2 over its base
          Defn: z2 |--> 1 - a

    list/tuple of them::

        sage: to_backend(([K, a], f))
        ([Finite Field in z2 of size 5^2, z2],
         Ring morphism:
           From: Finite Field in z2 of size 5^2
           To:   Field in a with defining polynomial x^2 + 4*x + 2 over its base
           Defn: z2 |--> 1 - a)

    and dictionaries::

        sage: to_backend({a: K})
        {z2: Finite Field in z2 of size 5^2}

    .. SEEALSO::

        :meth:`to_backend_parent`, :meth:`to_backend_element`, :meth:`to_backend_morphism`
    """
    from sage.rings.ring_extension import RingExtension_generic
    from sage.rings.ring_extension_morphism import RingExtensionHomomorphism, RingExtensionHomomorphism_baseinclusion
    if isinstance(arg, list):
        return [ to_backend(x) for x in arg ]
    elif isinstance(arg, tuple):
        return tuple([ to_backend(x) for x in arg ])
    elif isinstance(arg, dict):
        return { to_backend(key): to_backend(value) for (key, value) in arg.items() }
    elif isinstance(arg, RingExtension_generic):
        return arg._backend
    elif isinstance(arg, RingExtensionHomomorphism):
        return arg._backend
    elif isinstance(arg, RingExtensionHomomorphism_baseinclusion):
        return arg.codomain()._backend_defining_morphism
    elif isinstance(arg, RingExtensionElement):
        return arg._backend
    return arg

cpdef from_backend(arg, E):
    r"""
    Try to reconstruct something (somehow related to ``E``)
    whose backend is ``arg``.

    INPUT:

    - ``arg`` -- any argument

    - ``E`` -- a ring extension

    EXAMPLES::

        sage: from sage.rings.ring_extension_conversion import from_backend

    This function accepts parents::

        sage: K.<a> = GF(5^2).over()   # over GF(5)
        sage: from_backend(GF(5^2), K)
        Field in a with defining polynomial x^2 + 4*x + 2 over its base

    elements::

        sage: z2 = GF(5^2).gen()
        sage: from_backend(z2, K)
        a

    morphisms::

        sage: f = GF(5^2).frobenius_endomorphism()
        sage: from_backend(f, K)
        Ring endomorphism of Field in a with defining polynomial x^2 + 4*x + 2 over its base
          Defn: a |--> 1 - a

    list/tuple of them::

        sage: from_backend(([K, a], f), K)
        ([Field in a with defining polynomial x^2 + 4*x + 2 over its base, a],
         Ring endomorphism of Field in a with defining polynomial x^2 + 4*x + 2 over its base
           Defn: a |--> 1 - a)

    and dictionaries::

        sage: from_backend({a: K}, K)
        {a: Field in a with defining polynomial x^2 + 4*x + 2 over its base}

    .. SEEALSO::

        :meth:`from_backend_parent`, :meth:`from_backend_element`, :meth:`from_backend_morphism`
    """
    ans = None
    if isinstance(arg, list):
        ans = [ from_backend(x, E) for x in arg ]
    elif isinstance(arg, tuple):
        ans = tuple([ from_backend(x, E) for x in arg ])
    elif isinstance(arg, dict):
        ans = { from_backend(key, E): from_backend(value, E) for (key, value) in arg.items() }
    elif isinstance(arg, Parent):
        ans = from_backend_parent(arg, E)
    elif isinstance(arg, Map):
        ans = from_backend_morphism(arg, E)
    elif isinstance(arg, Element):
        ans = from_backend_element(arg, E)

    if ans is None:
        return arg
    else:
        return ans
