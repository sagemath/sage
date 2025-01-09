r"""
    This module implements the arithemtic for Jacobians of
    hyperelliptic curves. 

    AUTHORS:

    - David Kohel (2006): initial version
    - Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                2024 Sabrina Kunzweiler, Gareth Ma, Giacomo Pope
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
    r"""
    This is the base class for Jacobians of hyperelliptic curves.

    We represent elements of the Jacobian by tuples of the form 
    `(u, v : n)`, where
    - (u,v) is the Mumford representative of a divisor `P_1 + ... + P_r`,
    - n is a non-negative integer

    This tuple represents the equivalence class

    ..MATH::

        [P_1 + ... + P_r + n \cdot \infty_+ + m\cdot \infty_- - D_\infty],
        
    where  `m = g - \deg(u) - n`, and `\infty_+`, \infty_-` are the 
    points at infinity of the hyperelliptic curve,

    ..MATH::
        D_\infty =
        \lceil g/2 \rceil \infty_+ + \lfloor g/2 \rfloor \infty_-.
        
    Here, `\infty_- = \infty_+`, if the hyperelliptic curve is ramified.
        
    Such a representation exists and is unique, unless the genus `g` is odd 
    and the curve is inert. 
        
    If the hyperelliptic curve is ramified or inert, then `n` can be deduced 
    from `\deg(u)` and `g`. In these cases, `n` is omitted in the description.
    """

    def dimension(self):
        """
        Return the dimension of this Jacobian.

        EXAMPLES: 
            
            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^2, x^4+1); H
            Hyperelliptic Curve over Rational Field defined by y^2 + (x^4 + 1)*y = x^2
            sage: J = Jacobian(H)
            sage: J.dimension()
            3
        """
        return Integer(self.curve().genus())

    def point(self, *mumford, check=True, **kwargs):
        r"""
        Return a point on the Jacobian, given:

        1. No arguments or the integer `0`; return `0 \in J`;

        2. A point `P` on `J = Jac(C)`, return `P`;

        3. A point `P` on the curve `H` such that `J = Jac(H)`;
           return `[P - P_0]`, where `P_0` is the distinguished point of `H`.
           By default, `P_0 = \infty`;

        4. Two points `P, Q` on the curve `H` such that `J = Jac(H)`;
           return `[P - Q]`;

        5. Polynomials `(u, v)` such that `v^2 + hv - f \equiv 0 \pmod u`;
           return `[(u(x), y - v(x))]`.    
        
        .. SEEALSO:: 
        
            :mod:`sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic`.
        """
        try:
            return self.point_homset()(*mumford, check=check)
        except AttributeError:
            raise ValueError("Arguments must determine a valid Mumford divisor.")

    def _point_homset(self, *args, **kwds):
        # TODO: make a constructor for this??
        H = self.curve()
        if H.is_ramified():
            return jacobian_homset_ramified.HyperellipticJacobianHomsetRamified(
                *args, **kwds
            )
        elif H.is_split():
            return jacobian_homset_split.HyperellipticJacobianHomsetSplit(*args, **kwds)
        return jacobian_homset_inert.HyperellipticJacobianHomsetInert(*args, **kwds)

    def _point(self, *args, **kwds):
        H = self.curve()
        if H.is_ramified():
            return jacobian_morphism.MumfordDivisorClassFieldRamified(*args, **kwds)
        elif H.is_split():
            return jacobian_morphism.MumfordDivisorClassFieldSplit(*args, **kwds)
        return jacobian_morphism.MumfordDivisorClassFieldInert(*args, **kwds)

    @cached_method
    def order(self):
        """
        Compute the order of the Jacobian.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic.order`.
        """
        return self.point_homset().order()

    def count_points(self, *args, **kwds):
        """
        .. SEEALSO::
        
            :meth:`sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic.count_points`.
        """
        return self.point_homset().count_points(*args, **kwds)

    def lift_u(self, *args, **kwds):
        """
        Return one or all points with given `u`-coordinate.

        .. SEEALSO::
        
            :meth:`sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic.lift_u`.
        """
        return self.point_homset().lift_u(*args, **kwds)

    def random_element(self, *args, **kwds):
        """
        Return a random element of the Jacobian.

        .. SEEALSO::
        
            :meth:`sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic.random_element`.
        """
        return self.point_homset().random_element(*args, **kwds)

    def points(self, *args, **kwds):
        """
        Return all points on the Jacobian.

        .. SEEALSO::
        
            :meth:`sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic.points`.
        """

        return self.point_homset().points(*args, **kwds)

    def list(self):
        """
        Return all points on the Jacobian.

        .. SEEALSO::
        
            :meth:`sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic.points`.
        """

        return self.point_homset().points()

    def __iter__(self):
        yield from self.list()

    rational_points = points
