r"""
Vector Bundle Fiber Elements

The class :class:`VectorBundleFiberElement` implements vectors in the fiber of
a vector bundle.

AUTHORS:

- Michael Jung (2019): initial version
"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement

class VectorBundleFiberElement(FiniteRankFreeModuleElement):
    r"""
    Vector in a fiber of a vector bundle at the given point.

    INPUT:

    - ``parent`` -- :class:`~sage.manifolds.vector_bundle_fiber.VectorBundleFiber`;
      the fiber to which the vector belongs
    - ``name`` -- (default: ``None``) string; symbol given to the vector
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the vector; if ``None``, ``name`` will be used

    EXAMPLES:

    A vector `v` in a fiber of a rank 2 vector bundle::

        sage: M = Manifold(2, 'M', structure='top')
        sage: X.<x,y> = M.chart()
        sage: p = M((1,-1), name='p')
        sage: E = M.vector_bundle(2, 'E')
        sage: e = E.local_frame('e')
        sage: Ep = E.fiber(p)
        sage: v = Ep((-2,1), name='v'); v
        Vector v in the fiber of E at Point p on the 2-dimensional topological
         manifold M
        sage: v.display()
        v = -2 e_0 + e_1
        sage: v.parent()
        Fiber of E at Point p on the 2-dimensional topological manifold M
        sage: v in Ep
        True

    .. SEEALSO::

        :class:`~sage.tensor.modules.free_module_element.FiniteRankFreeModuleElement`
        for more documentation.
    """
    def __init__(self, parent, name=None, latex_name=None):
        r"""
        Construct a vector in the given fiber of a given vector bundle.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: p = M((1,-1), name='p')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: Ep = E.fiber(p)
            sage: v = Ep.element_class(Ep, name='v') ; v
            Vector v in the fiber of E at Point p on the 2-dimensional
             topological manifold M
            sage: v[:] = 5, -3/2
            sage: TestSuite(v).run()
        """
        FiniteRankFreeModuleElement.__init__(self, parent, name=name,
                                             latex_name=latex_name)
        # Extra data (with respect to FiniteRankFreeModuleElement):
        self._point = parent._point
        self._vbundle = parent._vbundle

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: p = M((1,-1), name='p')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: Ep = E.fiber(p)
            sage: v = Ep([-3,2], name='v')
            sage: v._repr_()
            'Vector v in the fiber of E at Point p on the 2-dimensional
             topological manifold M'
            sage: repr(v)  # indirect doctest
            'Vector v in the fiber of E at Point p on the 2-dimensional
             topological manifold M'
        """
        desc = "Vector "
        if self._name:
            desc += str(self._name) + " "
        desc += "in the fiber of {} at {}".format(self._vbundle._name,
                                                  self._point)
        return desc
