"""
Jacobian of a general hyperelliptic curve
"""

# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.schemes.hyperelliptic_curves_smooth_model import (
    jacobian_homset_inert,
    jacobian_homset_ramified,
    jacobian_homset_split,
    jacobian_morphism,
)
from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic


class HyperellipticJacobian_generic(Jacobian_generic):
    """
    TODO
    """
    def dimension(self):
        """
        Return the dimension of this Jacobian.
        """
        return Integer(self.curve().genus())

    # TODO why is check passed here and not used
    def point(self, *mumford, check=True):
        try:
            return self.point_homset()(*mumford, check=check)
        except AttributeError:
            raise ValueError("Arguments must determine a valid Mumford divisor.")

    def _point_homset(self, *args, **kwds):
        # TODO: make a constructor for this??
        H = self.curve()
        if H.is_ramified():
            return jacobian_homset_ramified.HyperellipticJacobianHomsetRamified(*args, **kwds)
        elif  H.is_split():
            return jacobian_homset_split.HyperellipticJacobianHomsetSplit(*args, **kwds)
        return jacobian_homset_inert.HyperellipticJacobianHomsetInert(*args, **kwds)

    def _point(self, *args, **kwds):
        H = self.curve()
        if H.is_ramified():
            return jacobian_morphism.MumfordDivisorClassFieldRamified(*args, **kwds)
        elif  H.is_split():
            return jacobian_morphism.MumfordDivisorClassFieldSplit(*args, **kwds)
        return jacobian_morphism.MumfordDivisorClassFieldInert(*args, **kwds)

    # Stupid functions
    def zero(self):
        return self.point_homset().zero()

    @cached_method
    def order(self):
        return self.point_homset().order()

    def random_element(self, fast=True):
        return self.point_homset().random_element(fast=fast)
