r"""
Support for internal use of number fields in backends for polyhedral computations
"""

# ****************************************************************************
#  Copyright (C) 2016-2022 Matthias Köppe <mkoeppe at math.ucdavis.edu>
#                2016-2018 Travis Scrimshaw
#                2017      Jeroen Demeyer
#                2018-2020 Jean-Philippe Labbé
#                2019      Vincent Delecroix
#                2019-2021 Jonathan Kliem
#                2019-2021 Sophia Elia
#                2020      Frédéric Chapoton
#                2022      Yuan Zhou
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from .base import Polyhedron_base


def _number_field_elements_from_algebraics_list_of_lists_of_lists(listss, **kwds):
    r"""
    Like `number_field_elements_from_algebraics`, but for a list of lists of lists.

    EXAMPLES::

        sage: rt2 = AA(sqrt(2)); rt2                                                    # needs sage.rings.number_field sage.symbolic
        1.414213562373095?
        sage: rt3 = AA(sqrt(3)); rt3                                                    # needs sage.rings.number_field sage.symbolic
        1.732050807568878?
        sage: from sage.geometry.polyhedron.base_number_field import _number_field_elements_from_algebraics_list_of_lists_of_lists
        sage: K, results, hom = _number_field_elements_from_algebraics_list_of_lists_of_lists([[[rt2], [1]], [[rt3]], [[1], []]]); results                      # needs sage.rings.number_field sage.symbolic
        [[[-a^3 + 3*a], [1]], [[a^2 - 2]], [[1], []]]
    """
    from sage.rings.qqbar import number_field_elements_from_algebraics
    numbers = []
    for lists in listss:
        for list in lists:
            numbers.extend(list)
    K, K_numbers, hom = number_field_elements_from_algebraics(numbers, **kwds)
    g = iter(K_numbers)
    return K, [ [ [ next(g) for _ in list ] for list in lists ] for lists in listss ], hom


class Polyhedron_base_number_field(Polyhedron_base):

    def _compute_data_lists_and_internal_base_ring(self, data_lists, convert_QQ, convert_NF):
        r"""
        Compute data lists in Normaliz or ``number_field`` backend format and the internal base ring of the data.

        EXAMPLES::

            sage: # optional - pynormaliz, needs sage.rings.number_field
            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],
            ....:                base_ring=AA, backend='normaliz')
            sage: def convert_QQ(ieqs, eqs):
            ....:     return [[1000*x for x in ieq] for ieq in ieqs], \
            ....:            [[1000*x for x in eq] for eq in eqs]
            sage: def convert_NF(ieqs, eqs):
            ....:     return ieqs, eqs
            sage: p._compute_data_lists_and_internal_base_ring(
            ....:     [[[1]], [[1/2]]], convert_QQ, convert_NF)
            (([[1000]], [[500]]), Rational Field)
            sage: p._compute_data_lists_and_internal_base_ring(
            ....:     [[[AA(1)]], [[1/2]]], convert_QQ, convert_NF)
            (([[1000]], [[500]]), Rational Field)
            sage: p._compute_data_lists_and_internal_base_ring(
            ....:     [[[AA(sqrt(2))]], [[1/2]]], convert_QQ, convert_NF)
            ([[[a]], [[1/2]]],
             Number Field in a with defining polynomial y^2 - 2 with a = 1.414213562373095?)

        TESTS::

            sage: K.<a> = QuadraticField(-5)                                            # needs sage.rings.number_field
            sage: p = Polyhedron(base_ring=K,       # indirect doctest  # optional - pynormaliz, needs sage.rings.number_field
            ....:                backend='normaliz',
            ....:                vertices=[(a,1/2), (2,0), (4,5/6)])
            Traceback (most recent call last):
            ...
            ValueError: invalid base ring: Number Field in a ... is not real embedded

        Checks that :issue:`30248` is fixed::

            sage: q = Polyhedron(base_ring=AA,      # indirect doctest  # optional - pynormaliz, needs sage.rings.number_field
            ....:                backend='normaliz',
            ....:                rays=[(0, 0, 1), (0, 1, -1), (1, 0, -1)]); q
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 1 vertex and 3 rays
            sage: -q                                                    # optional - pynormaliz, needs sage.rings.number_field
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 1 vertex and 3 rays
        """
        from sage.categories.number_fields import NumberFields
        from sage.rings.real_double import RDF

        if self.base_ring() in (QQ, ZZ):
            internal_base_ring = QQ
            internal_data_lists = convert_QQ(*data_lists)
        else:
            # Allows to re-iterate if K is QQ below when data_lists contain
            # iterators:
            data_lists = [tuple(_) for _ in data_lists]
            internal_data_lists = convert_NF(*data_lists)
            if self.base_ring() in NumberFields():
                if not RDF.has_coerce_map_from(self.base_ring()):
                    raise ValueError("invalid base ring: {} is a number field that is not real embedded".format(self.base_ring()))
                internal_base_ring = self.base_ring()
            else:
                K, internal_data_lists, hom = _number_field_elements_from_algebraics_list_of_lists_of_lists(internal_data_lists, embedded=True)
                internal_base_ring = K
                if K is QQ:
                    # Compute it with Normaliz, not QNormaliz
                    internal_data_lists = convert_QQ(*[ [ [ QQ(x) for x in v ] for v in l]
                                                   for l in data_lists ])
        return internal_data_lists, internal_base_ring
