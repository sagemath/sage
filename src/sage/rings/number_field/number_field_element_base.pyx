r"""
Number field elements (abstract base class)
"""

# ****************************************************************************
#       Copyright (C) 2023 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cdef class NumberFieldElement_base(FieldElement):
    r"""
    Abstract base class for :class:`~sage.rings.number_field.number_field_element.NumberFieldElement`

    This class is defined for the purpose of :func:`isinstance` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^3 + x + 1)                                          # needs sage.rings.number_field
        sage: isinstance(a, sage.rings.number_field.number_field_element_base.NumberFieldElement_base)                  # needs sage.rings.number_field
        True

    By design, there is a unique direct subclass::

        sage: len(sage.rings.number_field.number_field_element_base.NumberFieldElement_base.__subclasses__()) <= 1
        True
    """

    pass
