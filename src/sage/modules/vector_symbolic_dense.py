"""
Dense vectors over the symbolic ring

Implements dense vectors over the symbolic ring.

AUTHORS:

- Robert Bradshaw (2011-05-25): Added more element-wise simplification methods

- Joris Vankerschaver (2011-05-15): Initial version

EXAMPLES::

    sage: x, y = var('x, y')
    sage: u = vector([sin(x)^2 + cos(x)^2, log(2*y) + log(3*y)]); u
    (cos(x)^2 + sin(x)^2, log(3*y) + log(2*y))
    sage: type(u)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>
    sage: u.simplify_full()
    (1, log(3*y) + log(2*y))

TESTS:

Check that the outcome of arithmetic with symbolic vectors is again
a symbolic vector (:issue:`11549`)::

    sage: v = vector(SR, [1, 2])
    sage: w = vector(SR, [sin(x), 0])
    sage: type(v)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>
    sage: type(w)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>
    sage: type(v + w)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>
    sage: type(-v)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>
    sage: type(5*w)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>

Test pickling/unpickling::

    sage: u = vector(SR, [sin(x^2)])
    sage: loads(dumps(u)) == u
    True
"""

# ****************************************************************************
#       Copyright (C) 2011 Joris Vankerschaver (jv@caltech.edu)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules import free_module_element
from sage.symbolic.expression import Expression


def apply_map(phi):
    """
    Return a function that applies ``phi`` to its argument.

    EXAMPLES::

        sage: from sage.modules.vector_symbolic_dense import apply_map
        sage: v = vector([1,2,3])
        sage: f = apply_map(lambda x: x+1)
        sage: f(v)
        (2, 3, 4)
    """
    def apply(self, *args, **kwds):
        """
        Generic function used to implement common symbolic operations
        elementwise as methods of a vector.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: v = vector([sin(x)^2 + cos(x)^2, log(x*y), sin(x/(x^2 + x)), factorial(x+1)/factorial(x)])
            sage: v.simplify_trig()
            (1, log(x*y), sin(1/(x + 1)), factorial(x + 1)/factorial(x))
            sage: v.canonicalize_radical()
            (cos(x)^2 + sin(x)^2, log(x) + log(y), sin(1/(x + 1)), factorial(x + 1)/factorial(x))
            sage: v.simplify_rational()
            (cos(x)^2 + sin(x)^2, log(x*y), sin(1/(x + 1)), factorial(x + 1)/factorial(x))
            sage: v.simplify_factorial()
            (cos(x)^2 + sin(x)^2, log(x*y), sin(x/(x^2 + x)), x + 1)
            sage: v.simplify_full()
            (1, log(x*y), sin(1/(x + 1)), x + 1)

            sage: v = vector([sin(2*x), sin(3*x)])
            sage: v.simplify_trig()
            (2*cos(x)*sin(x), (4*cos(x)^2 - 1)*sin(x))
            sage: v.simplify_trig(False)
            (sin(2*x), sin(3*x))
            sage: v.simplify_trig(expand=False)
            (sin(2*x), sin(3*x))
        """
        return self.apply_map(lambda x: phi(x, *args, **kwds))
    apply.__doc__ += "\nSee Expression." + phi.__name__ + "() for optional arguments."
    return apply


class Vector_symbolic_dense(free_module_element.FreeModuleElement_generic_dense):
    pass


# Add elementwise methods.
for method in ['simplify', 'simplify_factorial',
               'simplify_log', 'simplify_rational',
               'simplify_trig', 'simplify_full', 'trig_expand',
               'canonicalize_radical', 'trig_reduce']:
    setattr(Vector_symbolic_dense, method, apply_map(getattr(Expression, method)))
