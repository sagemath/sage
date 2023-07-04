"""
Vectors over callable symbolic rings

AUTHOR:
    -- Jason Grout (2010)

EXAMPLES::

    sage: f(r, theta, z) = (r*cos(theta), r*sin(theta), z)
    sage: f.parent()
    Vector space of dimension 3 over Callable function ring with arguments (r, theta, z)
    sage: f
    (r, theta, z) |--> (r*cos(theta), r*sin(theta), z)
    sage: f[0]
    (r, theta, z) |--> r*cos(theta)
    sage: f+f
    (r, theta, z) |--> (2*r*cos(theta), 2*r*sin(theta), 2*z)
    sage: 3*f
    (r, theta, z) |--> (3*r*cos(theta), 3*r*sin(theta), 3*z)
    sage: f*f # dot product
    (r, theta, z) |--> r^2*cos(theta)^2 + r^2*sin(theta)^2 + z^2
    sage: f.diff()(0,1,2) # the matrix derivative
    [cos(1)      0      0]
    [sin(1)      0      0]
    [     0      0      1]


TESTS::

    sage: f(u,v,w) = (2*u+v,u-w,w^2+u)
    sage: loads(dumps(f)) == f
    True


"""

#*****************************************************************************
#       Copyright (C) 2010 Jason Grout <jason-sage@creativetrax.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import free_module_element
from sage.symbolic.ring import SR


class Vector_callable_symbolic_dense(free_module_element.FreeModuleElement_generic_dense):
    def _repr_(self):
        """
        Returns the string representation of the vector

        EXAMPLES::

            sage: f(u,v,w) = (2*u+v,u-w,w^2+u)
            sage: f
            (u, v, w) |--> (2*u + v, u - w, w^2 + u)
            sage: r(t) = (cos(t), sin(t))
            sage: r
            t |--> (cos(t), sin(t))
        """
        ring = self.coordinate_ring()
        args = ring.arguments()
        repr_x=self.change_ring(SR)._repr_()
        if len(args) == 1:
            return "%s |--> %s" % (args[0], repr_x)
        else:
            args = ", ".join(map(str, args))
            return "(%s) |--> %s" % (args, repr_x)

    def _latex_(self):
        r"""
        Return the latex representation of the vector.

        EXAMPLES::

            sage: f(u,v,w) = (2*u+v,u-w,w^2+u)
            sage: f
            (u, v, w) |--> (2*u + v, u - w, w^2 + u)
            sage: latex(f)
            \left( u, v, w \right) \ {\mapsto} \ \left(2 \, u + v,\,u - w,\,w^{2} + u\right)
            sage: r(t) = (cos(t), sin(t))
            sage: r
            t |--> (cos(t), sin(t))
            sage: latex(r)
            t \ {\mapsto}\ \left(\cos\left(t\right),\,\sin\left(t\right)\right)
        """
        from sage.misc.latex import latex
        ring = self.coordinate_ring()
        args = ring.arguments()
        args = [latex(arg) for arg in args]
        latex_x = self.change_ring(SR)._latex_()
        if len(args) == 1:
            return r"%s \ {\mapsto}\ %s" % (args[0], latex_x)
        else:
            vars = ", ".join(args)
            return r"\left( %s \right) \ {\mapsto} \ %s" % (vars, latex_x)
