r"""
Helper classes for structural embeddings and isomorphisms of number fields

AUTHORS:

- Julian Rueth (2014-04-03): initial version

Consider the following fields `L` and `M`::

    sage: L.<a> = QuadraticField(2)
    sage: M.<a> = L.absolute_field()

Both produce the same extension of `\QQ`. However, they should not be
identical because `M` carries additional information::

    sage: L.structure()
    (Ring Coercion endomorphism of Number Field in a with defining polynomial x^2 - 2,
     Ring Coercion endomorphism of Number Field in a with defining polynomial x^2 - 2)
    sage: M.structure()
    (Isomorphism given by variable name change map:
      From: Number Field in a with defining polynomial x^2 - 2
      To:   Number Field in a with defining polynomial x^2 - 2,
     Isomorphism given by variable name change map:
      From: Number Field in a with defining polynomial x^2 - 2
      To:   Number Field in a with defining polynomial x^2 - 2)

This used to cause trouble with caching and made (absolute) number fields not
unique when they should have been. The underlying technical problem is that the
morphisms returned by ``structure()`` can only be defined once the fields in
question have been created. Therefore, these morphisms cannot be part of a key
which uniquely identifies a number field.

The classes defined in this file encapsulate information about these structure
morphisms which can be passed to the factory creating number fields. This makes
it possible to distinguish number fields which only differ in terms of these
structure morphisms::

    sage: L is M
    False
    sage: N.<a> = L.absolute_field()
    sage: M is N
    True

"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2014 Julian Rueth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation

class NumberFieldStructure(object):
    r"""
    Abstract base class encapsulating information about a number fields
    relation to other number fields.

    TESTS::

        sage: from sage.rings.number_field.structure import NumberFieldStructure
        sage: NumberFieldStructure()
        <sage.rings.number_field.structure.NumberFieldStructure object at 0x...>

    """
    def create_structure(self, field):
        r"""
        Return a tuple encoding structural information about ``field``.

        OUTPUT:

        Typically, the output is a pair of morphisms. The first one from
        ``field`` to a field from which ``field`` has been constructed and the
        second one its inverse. In this case, these morphisms are used as
        conversion maps between the two fields.

        TESTS::

            sage: from sage.rings.number_field.structure import NumberFieldStructure
            sage: NumberFieldStructure().create_structure(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError

        The morphisms created by this method are used as conversion maps::

            sage: K.<i> = QuadraticField(-1)
            sage: L.<j> = K.change_names()
            sage: isinstance(L._structure, NumberFieldStructure)
            True
            sage: from_L, to_L = L.structure()
            sage: L._convert_map_from_(K) is to_L
            True
            sage: L(i)
            j
            sage: K(j)
            i

        """
        raise NotImplementedError

class NumberFieldStructureFromUniqueField(UniqueRepresentation, NumberFieldStructure):
    r"""
    Abstract base class encapsulating information about a number fields
    relation to another number field.

    INPUT:

    - ``other`` -- a number field

    TESTS:

        sage: from sage.rings.number_field.structure import NumberFieldStructureFromUniqueField
        sage: NumberFieldStructureFromUniqueField(QQ)
        <sage.rings.number_field.structure.NumberFieldStructureFromUniqueField object at 0x...>

    Instances are cached through
    :class:`sage.structure.unique_representation.UniqueRepresentation`::

        sage: NumberFieldStructureFromUniqueField(QQ) is NumberFieldStructureFromUniqueField(QQ)
        True

    However, the caching (and the `==` operator of the objects) is on identity
    and not on equality::

        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2+1)
        sage: L = K.change_names('j').change_names('i')
        sage: K == L
        True
        sage: K is L # K and L differ in "structure", one is the "name-change" of the other
        False

        sage: NumberFieldStructureFromUniqueField(L) is NumberFieldStructureFromUniqueField(L)
        True
        sage: NumberFieldStructureFromUniqueField(K) is NumberFieldStructureFromUniqueField(L)
        False

    This important because otherwise caching of number fields would be broken::

        sage: R.<x> = QQ[]
        sage: from sage.rings.number_field.structure import NameChange
        sage: KK.<j> = NumberField(x^2+1, structure=NameChange(K))
        sage: LL.<j> = NumberField(x^2+1, structure=NameChange(L))
        sage: KK is LL
        False

    """
    @staticmethod
    def __classcall__(cls, other, *args, **kwargs):
        r"""
        To make sure that caching only depends on the id of ``other`` and not
        on its implementation of ``==``, we add ``id(other)`` to the cache key.

        TESTS:

        If ``other`` is a tuple, its second entry is ignored, it is always overwritten by the preprocessing::

            sage: from sage.rings.number_field.structure import NumberFieldStructureFromUniqueField
            sage: NumberFieldStructureFromUniqueField((QQ, 1)) is NumberFieldStructureFromUniqueField((QQ, 2)) is NumberFieldStructureFromUniqueField(QQ)
            True

        """
        if isinstance(other, tuple):
            other = other[0]
        return super(NumberFieldStructureFromUniqueField,cls).__classcall__(cls, (other, id(other)), *args, **kwargs)

    def __init__(self, other):
        r"""
        Initialization.

        TESTS::

            sage: from sage.rings.number_field.structure import NumberFieldStructureFromUniqueField
            sage: type(NumberFieldStructureFromUniqueField(QQ))
            <class 'sage.rings.number_field.structure.NumberFieldStructureFromUniqueField'>

        """
        self.other = other[0]

class NameChange(NumberFieldStructureFromUniqueField):
    r"""
    Structure for a number field created by a change in variable name.

    INPUT:

    - ``other`` -- the number field from which this field has been created.

    TESTS::

        sage: from sage.rings.number_field.structure import NameChange
        sage: K.<i> = QuadraticField(-1)
        sage: NameChange(K)
        <sage.rings.number_field.structure.NameChange object at 0x...>

    """
    def create_structure(self, field):
        r"""
        Return a pair of isomorphisms which send the generator of ``field`` to
        the generator of ``other`` and vice versa.

        TESTS::

            sage: CyclotomicField(5).absolute_field('a').structure() # indirect doctest
            (Isomorphism given by variable name change map:
              From: Number Field in a with defining polynomial x^4 + x^3 + x^2 + x + 1
              To:   Cyclotomic Field of order 5 and degree 4,
             Isomorphism given by variable name change map:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Number Field in a with defining polynomial x^4 + x^3 + x^2 + x + 1)

        """
        from . import maps
        return maps.NameChangeMap(field, self.other), maps.NameChangeMap(self.other, field)

class AbsoluteFromRelative(NumberFieldStructureFromUniqueField):
    r"""
    Structure for an absolute number field created from a relative number
    field.

    INPUT:

    - ``other`` -- the number field from which this field has been created.

    TESTS::

        sage: from sage.rings.number_field.structure import AbsoluteFromRelative
        sage: K.<a> = QuadraticField(2)
        sage: R.<x> = K[]
        sage: L.<b> = K.extension(x^2 - 3)
        sage: AbsoluteFromRelative(L)
        <sage.rings.number_field.structure.AbsoluteFromRelative object at 0x...>

    """
    def create_structure(self, field):
        r"""
        Return a pair of isomorphisms which go from ``field`` to ``other`` and
        vice versa.

        TESTS::

            sage: K.<a> = QuadraticField(2)
            sage: R.<x> = K[]
            sage: L.<b> = K.extension(x^2 - 3)
            sage: M.<c> = L.absolute_field()
            sage: M.structure() # indirect doctest
            (Isomorphism map:
              From: Number Field in c with defining polynomial x^4 - 10*x^2 + 1
              To:   Number Field in b with defining polynomial x^2 - 3 over its base field, Isomorphism map:
              From: Number Field in b with defining polynomial x^2 - 3 over its base field
              To:   Number Field in c with defining polynomial x^4 - 10*x^2 + 1)

        """
        from . import maps
        return maps.MapAbsoluteToRelativeNumberField(field, self.other), maps.MapRelativeToAbsoluteNumberField(self.other, field)

class RelativeFromAbsolute(NumberFieldStructureFromUniqueField):
    r"""
    Structure for a relative number field created from an absolute number
    field.

    INPUT:

    - ``other`` -- the (absolute) number field from which this field has been
      created.

    - ``gen`` -- the generator of the intermediate field

    TESTS::

        sage: from sage.rings.number_field.structure import RelativeFromAbsolute
        sage: RelativeFromAbsolute(QQ, 1/2)
        <sage.rings.number_field.structure.RelativeFromAbsolute object at 0x...>

    """
    def __init__(self, other, gen):
        r"""
        Initialization.

        TESTS::

            sage: from sage.rings.number_field.structure import RelativeFromAbsolute
            sage: type(RelativeFromAbsolute(QQ, 1/2))
            <class 'sage.rings.number_field.structure.RelativeFromAbsolute'>

        """
        NumberFieldStructureFromUniqueField.__init__(self, other)
        self.gen = gen

    def create_structure(self, field):
        r"""
        Return a pair of isomorphisms which go from ``field`` to ``other`` and
        vice versa.

        INPUT:

        - ``field`` -- a relative number field

        TESTS::

            sage: K.<a> = QuadraticField(2)
            sage: M.<b,a_> = K.relativize(-a)
            sage: M.structure() # indirect doctest
            (Relative number field morphism:
             From: Number Field in b with defining polynomial x + a_ over its base field
             To:   Number Field in a with defining polynomial x^2 - 2
             Defn: -a_ |--> a
                   a_ |--> -a, Ring morphism:
             From: Number Field in a with defining polynomial x^2 - 2
             To:   Number Field in b with defining polynomial x + a_ over its base field
             Defn: a |--> -a_)

        """
        # other     field
        #    \       /
        #     \   Q(gen)
        #      \ /
        #       Q
        other = self.other
        gen = self.gen

        # the isomorphism from left to right is easy since the generators of
        # other and field are the same
        other_to_field = other.hom([field.gen(0)], field, check=True)
        assert other_to_field(gen) == field(field.base_field().gen())

        # to go from right to left, we first define a map from Q(gen) to other
        base_hom = field.base_field().hom([gen], other)
        # and extend it to a map from field
        field_to_other = field.Hom(other)([other.gen()], base_hom=base_hom, check=True)

        return field_to_other, other_to_field

class RelativeFromRelative(NumberFieldStructureFromUniqueField):
    r"""
    Structure for a relative number field created from another relative number
    field.

    INPUT:

    - ``other`` -- the relative number field used in the construction, see
      :meth:`create_structure`; there this field will be called ``field_``.

    TESTS::

        sage: from sage.rings.number_field.structure import RelativeFromRelative
        sage: K.<i> = QuadraticField(-1)
        sage: R.<x> = K[]
        sage: L.<a> = K.extension(x^2 - 2)
        sage: RelativeFromRelative(L)
        <sage.rings.number_field.structure.RelativeFromRelative object at 0x...>

    """
    def create_structure(self, field):
        r"""
        Return a pair of isomorphisms which go from ``field`` to the relative
        number field (called ``other`` below) from which ``field`` has been
        created and vice versa.

        The isomorphism is created via the relative number field ``field_``
        which is identical to ``field`` but is equipped with an isomorphism to
        an absolute field which was used in the construction of ``field``.

        INPUT:

        - ``field`` -- a relative number field

        TESTS::

            sage: K.<i> = QuadraticField(-1)
            sage: R.<x> = K[]
            sage: L.<a> = K.extension(x^2 - 2)
            sage: M.<b,a> = L.relativize(a)
            sage: M.structure() # indirect doctest
            (Relative number field morphism:
              From: Number Field in b with defining polynomial x^2 - 2*a*x + 3 over its base field
              To:   Number Field in a with defining polynomial x^2 - 2 over its base field
              Defn: b |--> a - i
                    a |--> a, Relative number field morphism:
              From: Number Field in a with defining polynomial x^2 - 2 over its base field
              To:   Number Field in b with defining polynomial x^2 - 2*a*x + 3 over its base field
              Defn: a |--> a
                    i |--> -b + a)

        """
        # other and field_ are relative number fields which are isomorphic via
        # an absolute number field abs:
        #
        #       f      g
        # other -> abs -> field_ = field
        field_ = self.other
        g_, g = field_.structure()
        abs = g.domain()
        f_, f = abs.structure()
        other = f.domain()

        # field and field_ are identical except that field_.structure() returns
        # the isomorphism with abs and field returns an isomorphism with other
        # (which we construct in this method).
        #
        # First, we construct the isomorphism from other to field by embedding
        # other.base_field() into field.
        gf = g*f
        base = other.base_field()
        base_to_field = base.Hom(field)([gf(base.gen())])
        other_to_field = other.Hom(field)([gf(other.gen())], base_hom=base_to_field)

        # And its inverse, essentially the same construction:
        f_g_ = f_*g_
        base = field.base_field()
        base_to_other = base.Hom(other)([f_g_(base.gen())])
        field_to_other = field.Hom(other)([f_g_(field.gen())], base_hom=base_to_other)

        return field_to_other, other_to_field
