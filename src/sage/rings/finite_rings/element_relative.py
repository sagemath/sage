# -*- coding: utf-8 -*-
r"""
Elements of Relative Extensions of Finite Fields

EXAMPLES::

    sage: set_random_seed(4)
    sage: k = GF(4).extension(2, absolute=False)
    sage: k.random_element()
    (z2 + 1)*z4

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
from sage.structure.element import CommutativeRingElement
from sage.rings.ring_extension_element import RingExtensionWithBasisElement

class FiniteField_relativeElement(RingExtensionWithBasisElement):
    r"""
    Element of a
    :class:`sage.rings.finite_rings.finite_field_relative.FiniteField_relative`.

    EXAMPLES::

        sage: k = GF(9).extension(2, absolute=False)
        sage: a = k.gen(); a
        z4
        sage: k.modulus()(a) == 0
        True

    """
    def __init__(self, parent, x):
        r"""
        TESTS::

            sage: from sage.rings.finite_rings.element_relative import FiniteField_relativeElement
            sage: k = GF(8).extension(2, absolute=False)
            sage: isinstance(k.gen(), FiniteField_relativeElement)
            True

        """
#        if x.parent() is not parent:
#            raise ValueError("parent must be %s but it is %s"%(parent, x.parent()))
        RingExtensionWithBasisElement.__init__(self, parent, x)
