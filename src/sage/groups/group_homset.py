"""
Set of homomorphisms between two groups

"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.all import HomsetWithBase, Groups
import sage.rings.integer_ring

def is_GroupHomset(H):
    return isinstance(H, GroupHomset_generic)

def GroupHomset(G, H):
    return RingHomset_generic(G, H)


class GroupHomset_generic(HomsetWithBase):
    """
    This class will not work since morphism.GroupHomomorphism_coercion
    is undefined and morphism.GroupHomomorphism_im_gens is undefined.
    """
    def __init__(self, G, H):
        HomsetWithBase.__init__(self, G, H, Groups(), sage.rings.integer_ringer.ZZ)

    def _repr_(self):
        return "Set of Homomorphisms from %s to %s"%(self.domain(), self.codomain())

    def __call__(self, im_gens, check=True):
        """
        EXAMPLES:

        """
        try:
            return morphism.GroupHomomorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError) as err:
            raise TypeError("images (=%s) do not define a valid homomorphism"%im_gens)

    def natural_map(self):
        return morphism.GroupHomomorphism_coercion(self)





