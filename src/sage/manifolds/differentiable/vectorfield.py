r"""
Vector Fields

Given two differentiable manifolds `U` and `M` over the same topological field
`K` and a differentiable map

.. MATH::

    \Phi:\ U \longrightarrow  M,

we define a *vector field along* `U` *with values on* `M` to be a
differentiable map

.. MATH::

    v:\ U  \longrightarrow TM

(`TM` being the tangent bundle of `M`) such that

.. MATH::

    \forall p \in U,\ v(p) \in T_{\Phi(p)}M.

The standard case of vector fields *on* a differentiable manifold corresponds
to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
being an immersion and `\Phi` being a curve in `M` (`U` is then an open
interval of `\RR`).

Vector fields are implemented via two classes: :class:`VectorFieldParal` and
:class:`VectorField`, depending respectively whether the manifold `M`
is parallelizable or not, i.e. whether the bundle `TM` is trivial or not.


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Marco Mancini (2015): parallelization of vector field plots
- Travis Scrimshaw (2016): review tweaks
- Eric Gourgoulhon (2017): vector fields inherit from multivector fields
- Eric Gourgoulhon (2018): dot and cross products, operators norm and curl

REFERENCES:

- [KN1963]_
- [Lee2013]_
- [ONe1983]_
- [BG1988]_
"""

# *****************************************************************************
#       Copyright (C) 2015, 2017 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement
from sage.manifolds.differentiable.multivectorfield import (
                                       MultivectorField, MultivectorFieldParal)
from sage.misc.decorators import options


class VectorField(MultivectorField):
    r"""
    Vector field along a differentiable manifold.

    An instance of this class is a vector field along a differentiable
    manifold `U` with values on a differentiable manifold `M`, via a
    differentiable map `U \rightarrow M`. More precisely, given a
    differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    a *vector field along* `U` *with values on* `M` is a differentiable map

    .. MATH::

        v:\ U  \longrightarrow TM

    (`TM` being the tangent bundle of `M`) such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M.

    The standard case of vector fields *on* a differentiable manifold
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are
    `\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is then an
    open interval of `\RR`).

    .. NOTE::

        If `M` is parallelizable, then
        :class:`~sage.manifolds.differentiable.vectorfield.VectorFieldParal`
        *must* be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `M\supset\Phi(U)`
    - ``name`` -- (default: ``None``) name given to the vector field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the vector
      field; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A vector field on a non-parallelizable 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_tu.<t,u> = V.chart()
        sage: transf = c_xy.transition_map(c_tu, (x+y, x-y), intersection_name='W',
        ....:                              restrictions1= x>0, restrictions2= t+u>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_tu.frame()
        sage: c_tuW = c_tu.restrict(W) ; eVW = c_tuW.frame()
        sage: v = M.vector_field(name='v') ; v
        Vector field v on the 2-dimensional differentiable manifold M
        sage: v.parent()
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M

    The vector field is first defined on the domain `U` by means of its
    components with respect to the frame ``eU``::

        sage: v[eU,:] = [-y, 1+x]

    The components with respect to the frame ``eV`` are then deduced
    by continuation of the components with respect to the frame ``eVW``
    on the domain `W = U \cap V`, expressed in terms on the coordinates
    covering `V`::

        sage: v[eV,0] = v[eVW,0,c_tuW].expr()
        sage: v[eV,1] = v[eVW,1,c_tuW].expr()

    At this stage, the vector field is fully defined on the whole manifold::

        sage: v.display(eU)
        v = -y ∂/∂x + (x + 1) ∂/∂y
        sage: v.display(eV)
        v = (u + 1) ∂/∂t + (-t - 1) ∂/∂u

    The vector field acting on scalar fields::

        sage: f = M.scalar_field({c_xy: (x+y)^2, c_tu: t^2}, name='f')
        sage: s = v(f) ; s
        Scalar field v(f) on the 2-dimensional differentiable manifold M
        sage: s.display()
        v(f): M → ℝ
        on U: (x, y) ↦ 2*x^2 - 2*y^2 + 2*x + 2*y
        on V: (t, u) ↦ 2*t*u + 2*t

    Some checks::

        sage: v(f) == f.differential()(v)
        True
        sage: v(f) == f.lie_der(v)
        True

    The result is defined on the intersection of the vector field's
    domain and the scalar field's one::

        sage: s = v(f.restrict(U)) ; s
        Scalar field v(f) on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s == v(f).restrict(U)
        True
        sage: s = v(f.restrict(W)) ; s
        Scalar field v(f) on the Open subset W of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): W → ℝ
           (x, y) ↦ 2*x^2 - 2*y^2 + 2*x + 2*y
           (t, u) ↦ 2*t*u + 2*t
        sage: s = v.restrict(U)(f) ; s
        Scalar field v(f) on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): U → ℝ
           (x, y) ↦ 2*x^2 - 2*y^2 + 2*x + 2*y
        on W: (t, u) ↦ 2*t*u + 2*t
        sage: s = v.restrict(U)(f.restrict(V)) ; s
        Scalar field v(f) on the Open subset W of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): W → ℝ
           (x, y) ↦ 2*x^2 - 2*y^2 + 2*x + 2*y
           (t, u) ↦ 2*t*u + 2*t
    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        r"""
        Construct a vector field with values on a non-parallelizable manifold.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``VectorField``, to fit with the category framework::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: XM = M.vector_field_module()
            sage: a = XM.element_class(XM, name='a'); a
            Vector field a on the 2-dimensional differentiable manifold M
            sage: a[c_xy.frame(),:] = [x, y]
            sage: a[c_uv.frame(),:] = [-u, -v]
            sage: TestSuite(a).run(skip='_test_pickling')

        Construction with ``DifferentiableManifold.vector_field``::

            sage: a1 = M.vector_field(name='a'); a1
            Vector field a on the 2-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True

        .. TODO::

            Fix ``_test_pickling`` (in the superclass :class:`TensorField`).
        """
        MultivectorField.__init__(self, vector_field_module, 1, name=name,
                                  latex_name=latex_name)
        # Initialization of derived quantities:
        MultivectorField._init_derived(self)
        # Initialization of list of quantities depending on self:
        self._init_dependencies()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._repr_()
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: repr(v)  # indirect doctest
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: v  # indirect doctest
            Vector field v on the 2-dimensional differentiable manifold M
        """
        description = "Vector field "
        if self._name is not None:
            description += self._name + " "
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same module.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: u = v._new_instance(); u
            Vector field on the 2-dimensional differentiable manifold M
            sage: u.parent() is v.parent()
            True
        """
        return type(self)(self._vmodule)

    def _init_dependencies(self):
        r"""
        Initialize list of quantities that depend on ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._init_dependencies()
        """
        self._lie_der_along_self = {}

    def _del_dependencies(self):
        r"""
        Clear list of quantities that depend on ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._del_dependencies()
        """
        if self._lie_der_along_self != {}:
            for tens in self._lie_der_along_self.values():
                del tens._lie_derivatives[id(self)]
            self._lie_der_along_self.clear()

    def __call__(self, scalar):
        r"""
        Action on a scalar field (or on a 1-form).

        INPUT:

        - ``scalar`` -- scalar field `f`

        OUTPUT:

        - scalar field representing the derivative of `f` along the vector
          field, i.e. `v^i \frac{\partial f}{\partial x^i}`

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: a = M.vector_field({c_xy.frame(): [x, y],
            ....:                     c_uv.frame(): [-u, -v]}, name='a')
            sage: f = M.scalar_field({c_xy: atan(x^2+y^2),
            ....:                     c_uv: pi/2-atan(u^2+v^2)}, name='f')
            sage: s = a.__call__(f); s
            Scalar field a(f) on the 2-dimensional differentiable manifold M
            sage: s.display()
            a(f): M → ℝ
            on U: (x, y) ↦ 2*(x^2 + y^2)/(x^4 + 2*x^2*y^2 + y^4 + 1)
            on V: (u, v) ↦ 2*(u^2 + v^2)/(u^4 + 2*u^2*v^2 + v^4 + 1)
            sage: s == f.differential()(a)
            True
        """
        if scalar._tensor_type == (0,1):
            # This is actually the action of the vector field on a 1-form,
            # as a tensor field of type (1,0):
            return scalar(self)
        if scalar._tensor_type != (0,0):
            raise TypeError("the argument must be a scalar field")
        resu = scalar.differential()(self)
        if not resu.is_immutable():
            if self._name is not None and scalar._name is not None:
                name = f"{self._name}({scalar._name})"
            else:
                name = None
            if self._latex_name is not None and scalar._latex_name is not None:
                latex_name = fr"{self._latex_name}\left({scalar._latex_name}\right)"
            else:
                latex_name = None
            resu.set_name(name=name, latex_name=latex_name)
        return resu

    @options(max_range=8, scale=1, color='blue')
    def plot(self, chart=None, ambient_coords=None, mapping=None,
             chart_domain=None, fixed_coords=None, ranges=None,
             number_values=None, steps=None,
             parameters=None, label_axes=True, **extra_options):
        r"""
        Plot the vector field in a Cartesian graph based on the coordinates
        of some ambient chart.

        The vector field is drawn in terms of two (2D graphics) or three
        (3D graphics) coordinates of a given chart, called hereafter the
        *ambient chart*.
        The vector field's base points `p` (or their images `\Phi(p)` by some
        differentiable mapping `\Phi`) must lie in the ambient chart's domain.

        INPUT:

        - ``chart`` -- (default: ``None``) the ambient chart (see above); if
          ``None``, the default chart of the vector field's domain is used

        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2
          or 3 coordinates of the ambient chart in terms of which the plot
          is performed; if ``None``, all the coordinates of the ambient
          chart are considered

        - ``mapping`` -- :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          (default: ``None``); differentiable map `\Phi` providing the link
          between the vector field's domain and the ambient chart ``chart``;
          if ``None``, the identity map is assumed

        - ``chart_domain`` -- (default: ``None``) chart on the vector field's
          domain to define the points at which vector arrows are to be plotted;
          if ``None``, the default chart of the vector field's domain is used

        - ``fixed_coords`` -- (default: ``None``) dictionary with keys the
          coordinates of ``chart_domain`` that are kept fixed and with values
          the value of these coordinates; if ``None``, all the coordinates of
          ``chart_domain`` are used

        - ``ranges`` -- (default: ``None``) dictionary with keys the
          coordinates of ``chart_domain`` to be used and values tuples
          ``(x_min, x_max)`` specifying the coordinate range for the plot;
          if ``None``, the entire coordinate range declared during the
          construction of ``chart_domain`` is considered (with ``-Infinity``
          replaced by ``-max_range`` and ``+Infinity`` by ``max_range``)

        - ``number_values`` -- (default: ``None``) either an integer or a
          dictionary with keys the coordinates of ``chart_domain`` to be
          used and values the number of values of the coordinate for sampling
          the part of the vector field's domain involved in the plot ; if
          ``number_values`` is a single integer, it represents the number of
          values for all coordinates; if ``number_values`` is ``None``, it is
          set to 9 for a 2D plot and to 5 for a 3D plot

        - ``steps`` -- (default: ``None``) dictionary with keys the
          coordinates of ``chart_domain`` to be used and values the step
          between each constant value of the coordinate; if ``None``, the
          step is computed from the coordinate range (specified in ``ranges``)
          and ``number_values``; on the contrary, if the step is provided
          for some coordinate, the corresponding number of values is deduced
          from it and the coordinate range

        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the coordinate expression
          of the vector field (see example below)

        - ``label_axes`` -- boolean (default: ``True``); determining whether
          the labels of the coordinate axes of ``chart`` shall be added to
          the graph; can be set to ``False`` if the graph is 3D and must be
          superposed with another graph

        - ``color`` -- (default: ``'blue'``) color of the arrows representing
          the vectors

        - ``max_range`` -- (default: 8) numerical value substituted to
          ``+Infinity`` if the latter is the upper bound of the range of a
          coordinate for which the plot is performed over the entire coordinate
          range (i.e. for which no specific plot range has been set in
          ``ranges``); similarly ``-max_range`` is the numerical valued
          substituted for ``-Infinity``

        - ``scale`` -- (default: 1) value by which the lengths of the arrows
          representing the vectors is multiplied

        - ``**extra_options`` -- extra options for the arrow plot, like
          ``linestyle``, ``width`` or ``arrowsize`` (see
          :func:`~sage.plot.arrow.arrow2d` and
          :func:`~sage.plot.plot3d.shapes.arrow3d` for details)

        OUTPUT:

        - a graphic object, either an instance of
          :class:`~sage.plot.graphics.Graphics` for a 2D plot (i.e. based on
          2 coordinates of ``chart``) or an instance of
          :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot (i.e.
          based on 3 coordinates of ``chart``)

        EXAMPLES:

        Plot of a vector field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: v = M.vector_field(-y, x, name='v')
            sage: v.display()
            v = -y ∂/∂x + x ∂/∂y
            sage: v.plot()                                                              # needs sage.plot
            Graphics object consisting of 80 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(-y, x, name='v')
            g = v.plot()
            sphinx_plot(g)

        Plot with various options::

            sage: v.plot(scale=0.5, color='green', linestyle='--', width=1,             # needs sage.plot
            ....:        arrowsize=6)
            Graphics object consisting of 80 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(-y, x, name='v')
            g = v.plot(scale=0.5, color='green', linestyle='--', width=1, arrowsize=6)
            sphinx_plot(g)

        ::

            sage: v.plot(max_range=4, number_values=5, scale=0.5)                       # needs sage.plot
            Graphics object consisting of 24 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(-y, x, name='v')
            g = v.plot(max_range=4, number_values=5, scale=0.5)
            sphinx_plot(g)

        Plot using parallel computation::

            sage: Parallelism().set(nproc=2)
            sage: v.plot(scale=0.5,  number_values=10, linestyle='--', width=1,         # needs sage.plot
            ....:        arrowsize=6)
            Graphics object consisting of 100 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(-y, x, name='v')
            g = v.plot(scale=0.5,  number_values=10, linestyle='--', width=1, arrowsize=6)
            sphinx_plot(g)

        ::

            sage: Parallelism().set(nproc=1)  # switch off parallelization

        Plots along a line of fixed coordinate::

            sage: v.plot(fixed_coords={x: -2})                                          # needs sage.plot
            Graphics object consisting of 9 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(-y, x, name='v')
            g = v.plot(fixed_coords={x: -2})
            sphinx_plot(g)

        ::

            sage: v.plot(fixed_coords={y: 1})                                           # needs sage.plot
            Graphics object consisting of 9 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(-y, x, name='v')
            g = v.plot(fixed_coords={y: 1})
            sphinx_plot(g)

        Let us now consider a vector field on a 4-dimensional manifold::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: v = M.vector_field((t/8)^2, -t*y/4, t*x/4, t*z/4, name='v')
            sage: v.display()
            v = 1/64*t^2 ∂/∂t - 1/4*t*y ∂/∂x + 1/4*t*x ∂/∂y + 1/4*t*z ∂/∂z

        We cannot make a 4D plot directly::

            sage: v.plot()
            Traceback (most recent call last):
            ...
            ValueError: the number of ambient coordinates must be either 2 or 3, not 4

        Rather, we have to select some coordinates for the plot, via
        the argument ``ambient_coords``. For instance, for a 3D plot::

            sage: v.plot(ambient_coords=(x, y, z), fixed_coords={t: 1},         # long time, needs sage.plot
            ....:        number_values=4)
            Graphics3d Object

        .. PLOT::

            M = Manifold(4, 'M')
            X = M.chart('t x y z') ; t,x,y,z = X[:]
            v = M.vector_field((t/8)**2, -t*y/4, t*x/4, t*z/4, name='v')
            sphinx_plot(v.plot(ambient_coords=(x, y, z), fixed_coords={t: 1},
                               number_values=4))

        ::

            sage: v.plot(ambient_coords=(x, y, t), fixed_coords={z: 0},         # long time, needs sage.plot
            ....:        ranges={x: (-2,2), y: (-2,2), t: (-1, 4)},
            ....:        number_values=4)
            Graphics3d Object

        .. PLOT::

            M = Manifold(4, 'M')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            v = M.vector_field((t/8)**2, -t*y/4, t*x/4, t*z/4, name='v')
            sphinx_plot(v.plot(ambient_coords=(x, y, t), fixed_coords={z: 0},
                               ranges={x: (-2,2), y: (-2,2), t: (-1, 4)},
                               number_values=4))

        or, for a 2D plot::

            sage: v.plot(ambient_coords=(x, y), fixed_coords={t: 1, z: 0})      # long time, needs sage.plot
            Graphics object consisting of 80 graphics primitives

        .. PLOT::

            M = Manifold(4, 'M')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            v = M.vector_field((t/8)**2, -t*y/4, t*x/4, t*z/4, name='v')
            g = v.plot(ambient_coords=(x, y), fixed_coords={t: 1, z: 0})
            sphinx_plot(g)

        ::

            sage: v.plot(ambient_coords=(x, t), fixed_coords={y: 1, z: 0})      # long time, needs sage.plot
            Graphics object consisting of 72 graphics primitives

        .. PLOT::

            M = Manifold(4, 'M')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            v = M.vector_field((t/8)**2, -t*y/4, t*x/4, t*z/4, name='v')
            g = v.plot(ambient_coords=(x, t), fixed_coords={y: 1, z: 0})
            sphinx_plot(g)

        An example of plot via a differential mapping: plot of a vector field
        tangent to a 2-sphere viewed in `\RR^3`::

            sage: S2 = Manifold(2, 'S^2')
            sage: U = S2.open_subset('U') # the open set covered by spherical coord.
            sage: XS.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: R3 = Manifold(3, 'R^3')
            sage: X3.<x,y,z> = R3.chart()
            sage: F = S2.diff_map(R3, {(XS, X3): [sin(th)*cos(ph),
            ....:                       sin(th)*sin(ph), cos(th)]}, name='F')
            sage: F.display() # the standard embedding of S^2 into R^3
            F: S^2 → R^3
            on U: (th, ph) ↦ (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: v = XS.frame()[1] ; v  # the coordinate vector ∂/∂phi
            Vector field ∂/∂ph on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: graph_v = v.plot(chart=X3, mapping=F, label_axes=False)               # needs sage.plot
            sage: graph_S2 = XS.plot(chart=X3, mapping=F, number_values=9)              # needs sage.plot
            sage: graph_v + graph_S2                                                    # needs sage.plot
            Graphics3d Object

        .. PLOT::

            S2 = Manifold(2, 'S^2')
            U = S2.open_subset('U')
            XS = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            th, ph = XS[:]
            R3 = Manifold(3, 'R^3')
            X3 = R3.chart('x y z')
            F = S2.diff_map(R3, {(XS, X3): [sin(th)*cos(ph), sin(th)*sin(ph),
                                            cos(th)]}, name='F')
            v = XS.frame()[1]
            graph_v = v.plot(chart=X3, mapping=F, label_axes=False)
            graph_S2 = XS.plot(chart=X3, mapping=F, number_values=9)
            sphinx_plot(graph_v + graph_S2)

        Note that the default values of some arguments of the method ``plot``
        are stored in the dictionary ``plot.options``::

            sage: v.plot.options  # random (dictionary output)
            {'color': 'blue', 'max_range': 8, 'scale': 1}

        so that they can be adjusted by the user::

            sage: v.plot.options['color'] = 'red'

        From now on, all plots of vector fields will use red as the default
        color. To restore the original default options, it suffices to type::

            sage: v.plot.reset()
        """
        from sage.rings.infinity import Infinity
        from sage.misc.functional import numerical_approx
        from sage.misc.latex import latex
        from sage.plot.graphics import Graphics
        from sage.manifolds.chart import RealChart
        from sage.manifolds.utilities import set_axes_labels
        from sage.parallel.decorate import parallel
        from sage.parallel.parallelism import Parallelism

        #
        # 1/ Treatment of input parameters
        #    -----------------------------
        max_range = extra_options.pop("max_range")
        scale = extra_options.pop("scale")
        color = extra_options.pop("color")
        if chart is None:
            chart = self._domain.default_chart()
        elif not isinstance(chart, RealChart):
            raise TypeError("{} is not a chart on a real ".format(chart) +
                            "manifold")
        if chart_domain is None:
            chart_domain = self._domain.default_chart()
        elif not isinstance(chart_domain, RealChart):
            raise TypeError("{} is not a chart on a ".format(chart_domain) +
                            "real manifold")
        elif not chart_domain.domain().is_subset(self._domain):
            raise ValueError("the domain of {} is not ".format(chart_domain) +
                             "included in the domain of {}".format(self))
        coords_full = tuple(chart_domain[:]) # all coordinates of chart_domain
        if fixed_coords is None:
            coords = coords_full
        else:
            fixed_coord_list = fixed_coords.keys()
            coords = []
            for coord in coords_full:
                if coord not in fixed_coord_list:
                    coords.append(coord)
            coords = tuple(coords)
        if ambient_coords is None:
            ambient_coords = chart[:]
        elif not isinstance(ambient_coords, tuple):
            ambient_coords = tuple(ambient_coords)
        nca = len(ambient_coords)
        if nca != 2 and nca != 3:
            raise ValueError("the number of ambient coordinates must be " +
                             "either 2 or 3, not {}".format(nca))
        if ranges is None:
            ranges = {}
        ranges0 = {}
        for coord in coords:
            if coord in ranges:
                ranges0[coord] = (numerical_approx(ranges[coord][0]),
                                  numerical_approx(ranges[coord][1]))
            else:
                bounds = chart_domain._bounds[coords_full.index(coord)]
                xmin0 = bounds[0][0]
                xmax0 = bounds[1][0]
                if xmin0 == -Infinity:
                    xmin = numerical_approx(-max_range)
                elif bounds[0][1]:
                    xmin = numerical_approx(xmin0)
                else:
                    xmin = numerical_approx(xmin0 + 1.e-3)
                if xmax0 == Infinity:
                    xmax = numerical_approx(max_range)
                elif bounds[1][1]:
                    xmax = numerical_approx(xmax0)
                else:
                    xmax = numerical_approx(xmax0 - 1.e-3)
                ranges0[coord] = (xmin, xmax)
        ranges = ranges0
        if number_values is None:
            if nca == 2: # 2D plot
                number_values = 9
            else:   # 3D plot
                number_values = 5
        if not isinstance(number_values, dict):
            number_values0 = {}
            for coord in coords:
                number_values0[coord] = number_values
            number_values = number_values0
        if steps is None:
            steps = {}
        for coord in coords:
            if coord not in steps:
                steps[coord] = (ranges[coord][1] - ranges[coord][0]) / \
                               (number_values[coord]-1)
            else:
                number_values[coord] = 1 + int(
                           (ranges[coord][1] - ranges[coord][0]) / steps[coord])
        #
        # 2/ Plots
        #    -----
        dom = chart_domain.domain()
        vector = self.restrict(dom)
        if vector.parent().destination_map() is dom.identity_map():
            if mapping is not None:
                vector = mapping.pushforward(vector)
                mapping = None
        nc = len(coords_full)
        ncp = len(coords)
        xx = [0] * nc
        if fixed_coords is not None:
            if len(fixed_coords) != nc - ncp:
                raise ValueError("bad number of fixed coordinates")
            for fc, val in fixed_coords.items():
                xx[coords_full.index(fc)] = val
        ind_coord = []
        for coord in coords:
            ind_coord.append(coords_full.index(coord))

        resu = Graphics()
        ind = [0] * ncp
        ind_max = [0] * ncp
        ind_max[0] = number_values[coords[0]]
        xmin = [ranges[cd][0] for cd in coords]
        step_tab = [steps[cd] for cd in coords]

        nproc = Parallelism().get('tensor')
        if nproc != 1 and nca == 2:
            # parallel plot construct : Only for 2D plot (at  moment) !

            # creation of the list of parameters
            list_xx = []

            while ind != ind_max:
                for i in range(ncp):
                    xx[ind_coord[i]] = xmin[i] + ind[i] * step_tab[i]

                if chart_domain.valid_coordinates(*xx, tolerance=1e-13,
                                                  parameters=parameters):

                    # needed a xx*1 to copy the list by value
                    list_xx.append(xx*1)

                # Next index:
                ret = 1
                for pos in range(ncp-1,-1,-1):
                    imax = number_values[coords[pos]] - 1
                    if ind[pos] != imax:
                        ind[pos] += ret
                        ret = 0
                    elif ret == 1:
                        if pos == 0:
                            ind[pos] = imax + 1 # end point reached
                        else:
                            ind[pos] = 0
                            ret = 1

            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_step = max(1, int(len(list_xx)/nproc/2))
            local_list = lol(list_xx,ind_step)

            # definition of the list of input parameters
            listParalInput = [(vector, dom, ind_part,
                               chart_domain, chart,
                               ambient_coords, mapping,
                               scale, color, parameters,
                               extra_options)
                              for ind_part in local_list]

            # definition of the parallel function
            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def add_point_plot(vector, dom, xx_list, chart_domain, chart,
                               ambient_coords, mapping, scale, color,
                               parameters, extra_options):
                count = 0
                for xx in xx_list:
                    point = dom(xx, chart=chart_domain)
                    part = vector.at(point).plot(chart=chart,
                                                 ambient_coords=ambient_coords,
                                                 mapping=mapping,scale=scale,
                                                 color=color, print_label=False,
                                                 parameters=parameters,
                                                 **extra_options)
                    if count == 0:
                        local_resu = part
                    else:
                        local_resu += part
                    count += 1
                return local_resu

            # parallel execution and reconstruction of the plot
            for ii, val in add_point_plot(listParalInput):
                resu += val

        else:
            # sequential plot
            while ind != ind_max:
                for i in range(ncp):
                    xx[ind_coord[i]] = xmin[i] + ind[i]*step_tab[i]
                if chart_domain.valid_coordinates(*xx, tolerance=1e-13,
                                                  parameters=parameters):
                    point = dom(xx, chart=chart_domain)
                    resu += vector.at(point).plot(chart=chart,
                                                  ambient_coords=ambient_coords,
                                                  mapping=mapping, scale=scale,
                                                  color=color, print_label=False,
                                                  parameters=parameters,
                                                  **extra_options)
                # Next index:
                ret = 1
                for pos in range(ncp-1, -1, -1):
                    imax = number_values[coords[pos]] - 1
                    if ind[pos] != imax:
                        ind[pos] += ret
                        ret = 0
                    elif ret == 1:
                        if pos == 0:
                            ind[pos] = imax + 1 # end point reached
                        else:
                            ind[pos] = 0
                            ret = 1

        if label_axes:
            if nca == 2:  # 2D graphic
                # We update the dictionary _extra_kwds (options to be passed
                # to show()), instead of using the method
                # Graphics.axes_labels() since the latter is not robust w.r.t.
                # graph addition
                resu._extra_kwds['axes_labels'] = [r'$'+latex(ac)+r'$'
                                                   for ac in ambient_coords]
            else: # 3D graphic
                labels = [str(ac) for ac in ambient_coords]
                resu = set_axes_labels(resu, *labels)
        return resu

    def bracket(self, other):
        """
        Return the Lie bracket ``[self, other]``.

        INPUT:

        - ``other`` -- a :class:`VectorField`

        OUTPUT: the :class:`VectorField` ``[self, other]``

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: v = -X.frame()[0] + 2*X.frame()[1] - (x^2 - y)*X.frame()[2]
            sage: w = (z + y) * X.frame()[1] - X.frame()[2]
            sage: vw = v.bracket(w); vw
            Vector field on the 3-dimensional differentiable manifold M
            sage: vw.display()
            (-x^2 + y + 2) ∂/∂y + (-y - z) ∂/∂z

        Some checks::

            sage: vw == - w.bracket(v)
            True
            sage: f = M.scalar_field({X: x+y*z})
            sage: vw(f) == v(w(f)) - w(v(f))
            True
            sage: vw == w.lie_derivative(v)
            True
        """
        # Call of the Schouten-Nijenhuis bracket
        return MultivectorField.bracket(self, other)

    def curl(self, metric=None):
        r"""
        Return the curl of ``self`` with respect to a given metric, assuming
        that the domain of ``self`` is 3-dimensional.

        If ``self`` is a vector field `v` on a 3-dimensional differentiable
        orientable manifold `M`, the curl of `v` with respect to a metric `g`
        on `M` is the vector field defined by

        .. MATH::

            \mathrm{curl}\, v = (*(\mathrm{d} v^\flat))^\sharp

        where `v^\flat` is the 1-form associated to `v` by the metric `g` (see
        :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.down`),
        `*(\mathrm{d} v^\flat)` is the Hodge dual with respect to `g` of the
        2-form `\mathrm{d} v^\flat` (exterior derivative of `v^\flat`) (see
        :meth:`~sage.manifolds.differentiable.diff_form.DiffForm.hodge_dual`)
        and
        `(*(\mathrm{d} v^\flat))^\sharp` is corresponding vector field by
        `g`-duality (see
        :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.up`).

        An alternative expression of the curl is

        .. MATH::

            (\mathrm{curl}\, v)^i = \epsilon^{ijk} \nabla_j v_k

        where `\nabla` is the Levi-Civita connection of `g` (cf.
        :class:`~sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection`)
        and `\epsilon` the volume 3-form (Levi-Civita tensor) of `g` (cf.
        :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.volume_form`)

        .. NOTE::

            The method ``curl`` is meaningful only if ``self`` is a vector
            field on a 3-dimensional manifold.

        INPUT:

        - ``metric`` -- (default: ``None``) the pseudo-Riemannian metric `g`
          involved in the definition of the curl; if none is provided, the
          domain of ``self`` is supposed to be endowed with a default metric
          (i.e. is supposed to be pseudo-Riemannian manifold, see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the curl

        OUTPUT: instance of :class:`VectorField` representing the curl of ``self``

        EXAMPLES:

        Curl of a vector field in the Euclidean 3-space::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: v = M.vector_field(-y, x, 0, name='v')
            sage: v.display()
            v = -y e_x + x e_y
            sage: s = v.curl(); s
            Vector field curl(v) on the Euclidean space E^3
            sage: s.display()
            curl(v) = 2 e_z

        The function :func:`~sage.manifolds.operators.curl` from the
        :mod:`~sage.manifolds.operators` module can be used instead of the
        method :meth:`curl`::

            sage: from sage.manifolds.operators import curl
            sage: curl(v) == s
            True

        If one prefers the notation ``rot`` over ``curl``, it suffices to do::

            sage: from sage.manifolds.operators import curl as rot
            sage: rot(v) == s
            True

        The curl of a gradient vanishes identically::

            sage: f = M.scalar_field(function('F')(x,y,z))
            sage: gradf = f.gradient()
            sage: gradf.display()
            d(F)/dx e_x + d(F)/dy e_y + d(F)/dz e_z
            sage: s = curl(gradf); s
            Vector field on the Euclidean space E^3
            sage: s.display()
            0
        """
        if self._domain.dim() < 3:
            raise ValueError("the curl is not defined in dimension lower " +
                             "than 3")
        default_metric = metric is None
        if default_metric:
            metric = self._domain.metric()
        der = self.down(metric).exterior_derivative()  # 2-form d(v^\flat)
        resu = der.hodge_dual(metric).up(metric)
        if self._name is not None:
            if default_metric:
                resu._name = "curl({})".format(self._name)
                resu._latex_name = r"\mathrm{curl}\left(" + self._latex_name + \
                                   r"\right)"
            else:
                resu._name = "curl_{}({})".format(metric._name, self._name)
                resu._latex_name = r"\mathrm{curl}_{" + metric._latex_name + \
                                   r"}\left(" + self._latex_name + r"\right)"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    def dot_product(self, other, metric=None):
        r"""
        Return the scalar product of ``self`` with another vector field (with
        respect to a given metric).

        If ``self`` is the vector field `u` and other is the vector field `v`,
        the *scalar product of* `u` *by* `v` with respect to a given
        pseudo-Riemannian metric `g` is the scalar field `s` defined by

        .. MATH::

            s = u\cdot v = g(u,v) = g_{ij} u^i v^j

        INPUT:

        - ``other`` -- a vector field, defined on the same domain as ``self``
        - ``metric`` -- (default: ``None``) the pseudo-Riemannian metric `g`
          involved in the definition of the scalar product; if none is
          provided, the domain of ``self`` is supposed to be endowed with a
          default metric (i.e. is supposed to be pseudo-Riemannian manifold,
          see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the scalar product

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`
          representing the scalar product of ``self`` by ``other``.

        EXAMPLES:

        Scalar product in the Euclidean plane::

            sage: M.<x,y> = EuclideanSpace()
            sage: u = M.vector_field(x, y, name='u')
            sage: v = M.vector_field(y, x, name='v')
            sage: s = u.dot_product(v); s
            Scalar field u.v on the Euclidean plane E^2
            sage: s.display()
            u.v: E^2 → ℝ
               (x, y) ↦ 2*x*y

        A shortcut alias of ``dot_product`` is ``dot``::

            sage: u.dot(v) == s
            True

        A test of orthogonality::

            sage: v[:] = -y, x
            sage: u.dot_product(v) == 0
            True

        Scalar product with respect to a metric that is not the default one::

            sage: h = M.riemannian_metric('h')
            sage: h[1,1], h[2,2] = 1/(1+y^2), 1/(1+x^2)
            sage: s = u.dot_product(v, metric=h); s
            Scalar field h(u,v) on the Euclidean plane E^2
            sage: s.display()
            h(u,v): E^2 → ℝ
               (x, y) ↦ -(x^3*y - x*y^3)/((x^2 + 1)*y^2 + x^2 + 1)

        Scalar product of two vector fields along a curve (a lemniscate of
        Gerono)::

            sage: R.<t> = manifolds.RealLine()
            sage: C = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='C')
            sage: u = C.tangent_vector_field(name='u')
            sage: u.display()
            u = cos(t) e_x + (2*cos(t)^2 - 1) e_y
            sage: I = C.domain(); I
            Real interval (0, 2*pi)
            sage: v = I.vector_field(cos(t), -1, dest_map=C, name='v')
            sage: v.display()
            v = cos(t) e_x - e_y
            sage: s = u.dot_product(v); s
            Scalar field u.v on the Real interval (0, 2*pi)
            sage: s.display()
            u.v: (0, 2*pi) → ℝ
               t ↦ sin(t)^2

        Scalar product between a vector field along the curve and a vector
        field on the ambient Euclidean plane::

            sage: e_x = M.cartesian_frame()[1]
            sage: s = u.dot_product(e_x); s
            Scalar field u.e_x on the Real interval (0, 2*pi)
            sage: s.display()
            u.e_x: (0, 2*pi) → ℝ
               t ↦ cos(t)
        """
        default_metric = metric is None
        if default_metric:
            metric = self._ambient_domain.metric()
        dest_map = self.parent().destination_map()
        if dest_map != metric.parent().base_module().destination_map():
            metric = metric.along(dest_map)
        if dest_map != other.parent().destination_map():
            other = other.along(dest_map)
        resu = metric(self, other)
        # From the above operation the name of resu is "g(u,v')" where
        # g = metric._name, u = self._name, v = other._name
        # For a default metric, we change it to "u.v":
        if (default_metric and self._name is not None and
            other._name is not None):
            resu._name = "{}.{}".format(self._name, other._name)
            resu._latex_name = "{" + self._latex_name + r"}\cdot{" + \
                               other._latex_name + "}"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    dot = dot_product

    def norm(self, metric=None):
        r"""
        Return the norm of ``self`` (with respect to a given metric).

        The *norm* of a vector field `v` with respect to a given
        pseudo-Riemannian metric `g` is the scalar field `\|v\|` defined by

        .. MATH::

            \|v\| = \sqrt{g(v,v)}

        .. NOTE::

            If the metric `g` is not positive definite, it may be that `\|v\|`
            takes imaginary values.

        INPUT:

        - ``metric`` -- (default: ``None``) the pseudo-Riemannian metric `g`
          involved in the definition of the norm; if none is
          provided, the domain of ``self`` is supposed to be endowed with a
          default metric (i.e. is supposed to be pseudo-Riemannian manifold,
          see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the norm

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`
          representing the norm of ``self``.

        EXAMPLES:

        Norm in the Euclidean plane::

            sage: M.<x,y> = EuclideanSpace()
            sage: v = M.vector_field(-y, x, name='v')
            sage: s = v.norm(); s
            Scalar field |v| on the Euclidean plane E^2
            sage: s.display()
            |v|: E^2 → ℝ
               (x, y) ↦ sqrt(x^2 + y^2)

        The global function :func:`~sage.misc.functional.norm` can be used
        instead of the method ``norm()``::

            sage: norm(v) == s
            True

        Norm with respect to a metric that is not the default one::

            sage: h = M.riemannian_metric('h')
            sage: h[1,1], h[2,2] = 1/(1+y^2), 1/(1+x^2)
            sage: s = v.norm(metric=h); s
            Scalar field |v|_h on the Euclidean plane E^2
            sage: s.display()
            |v|_h: E^2 → ℝ
               (x, y) ↦ sqrt((2*x^2 + 1)*y^2 + x^2)/(sqrt(x^2 + 1)*sqrt(y^2 + 1))

        Norm of the tangent vector field to a curve (a lemniscate of Gerono)::

            sage: R.<t> = manifolds.RealLine()
            sage: C = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='C')
            sage: v = C.tangent_vector_field()
            sage: v.display()
            C' = cos(t) e_x + (2*cos(t)^2 - 1) e_y
            sage: s = v.norm(); s
            Scalar field |C'| on the Real interval (0, 2*pi)
            sage: s.display()
            |C'|: (0, 2*pi) → ℝ
               t ↦ sqrt(4*cos(t)^4 - 3*cos(t)^2 + 1)
        """
        default_metric = metric is None
        if default_metric:
            metric = self._ambient_domain.metric()
        dest_map = self.parent().destination_map()
        if dest_map != metric.parent().base_module().destination_map():
            metric = metric.along(dest_map)
        resu = metric(self, self).sqrt()
        if self._name is not None:
            if default_metric:
                resu._name = "|{}|".format(self._name)
                resu._latex_name = r"\left\|" + self._latex_name + \
                                   r"\right\|"
            else:
                resu._name = "|{}|_{}".format(self._name, metric._name)
                resu._latex_name = r"\left\|" + self._latex_name + \
                                   r"\right\| _{" + metric._latex_name + "}"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    def cross_product(self, other, metric=None):
        r"""
        Return the cross product of ``self`` with another vector field (with
        respect to a given metric),  assuming that the domain of ``self`` is
        3-dimensional.

        If ``self`` is a vector field `u` on a 3-dimensional differentiable
        orientable manifold `M` and ``other`` is a vector field `v` on `M`,
        the *cross product* (also called *vector product*) *of* `u` *by* `v`
        with respect to a pseudo-Riemannian metric `g` on `M` is the vector
        field `w = u\times v` defined by

        .. MATH::

            w^i = \epsilon^i_{\phantom{i} jk} u^j v^k
                = g^{il} \epsilon_{ljk} u^j v^k

        where `\epsilon` is the volume 3-form (Levi-Civita tensor) of `g` (cf.
        :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.volume_form`)

        .. NOTE::

            The method ``cross_product`` is meaningful only if for vector fields on a
            3-dimensional manifold.

        INPUT:

        - ``other`` -- a vector field, defined on the same domain as ``self``
        - ``metric`` -- (default: ``None``) the pseudo-Riemannian metric `g`
          involved in the definition of the cross product; if none is
          provided, the domain of ``self`` is supposed to be endowed with a
          default metric (i.e. is supposed to be pseudo-Riemannian manifold,
          see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the cross product

        OUTPUT:

        - instance of :class:`VectorField` representing the cross product of
          ``self`` by ``other``.

        EXAMPLES:

        Cross product in the Euclidean 3-space::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: u = M.vector_field(-y, x, 0, name='u')
            sage: v = M.vector_field(x, y, 0, name='v')
            sage: w = u.cross_product(v); w
            Vector field u x v on the Euclidean space E^3
            sage: w.display()
            u x v = (-x^2 - y^2) e_z

        A shortcut alias of ``cross_product`` is ``cross``::

            sage: u.cross(v) == w
            True

        The cross product of a vector field with itself is zero::

            sage: u.cross_product(u).display()
            u x u = 0

        Cross product with respect to a metric that is not the default one::

            sage: h = M.riemannian_metric('h')
            sage: h[1,1], h[2,2], h[3,3] = 1/(1+y^2), 1/(1+z^2), 1/(1+x^2)
            sage: w = u.cross_product(v, metric=h); w
            Vector field on the Euclidean space E^3
            sage: w.display()
            -(x^2 + y^2)*sqrt(x^2 + 1)/(sqrt(y^2 + 1)*sqrt(z^2 + 1)) e_z

        Cross product of two vector fields along a curve (arc of a helix)::

            sage: R.<t> = manifolds.RealLine()
            sage: C = M.curve((cos(t), sin(t), t), (t, 0, 2*pi), name='C')
            sage: u = C.tangent_vector_field()
            sage: u.display()
            C' = -sin(t) e_x + cos(t) e_y + e_z
            sage: I = C.domain(); I
            Real interval (0, 2*pi)
            sage: v = I.vector_field(-cos(t), sin(t), 0, dest_map=C)
            sage: v.display()
            -cos(t) e_x + sin(t) e_y
            sage: w = u.cross_product(v); w
            Vector field along the Real interval (0, 2*pi) with values on the
             Euclidean space E^3
            sage: w.parent().destination_map()
            Curve C in the Euclidean space E^3
            sage: w.display()
            -sin(t) e_x - cos(t) e_y + (2*cos(t)^2 - 1) e_z

        Cross product between a vector field along the curve and a vector field
        on the ambient Euclidean space::

            sage: e_x = M.cartesian_frame()[1]
            sage: w = u.cross_product(e_x); w
            Vector field C' x e_x along the Real interval (0, 2*pi) with values
             on the Euclidean space E^3
            sage: w.display()
            C' x e_x = e_y - cos(t) e_z
        """
        if self._ambient_domain.dim() != 3:
            raise ValueError("the cross product is not defined in dimension " +
                             "different from 3")
        default_metric = metric is None
        if default_metric:
            metric = self._ambient_domain.metric()
        dest_map = self.parent().destination_map()
        if dest_map == metric.parent().base_module().destination_map():
            eps = metric.volume_form(1)
        else:
            eps = metric.volume_form(1).along(dest_map)
        if dest_map != other.parent().destination_map():
            other = other.along(dest_map)
        resu = eps.contract(1, 2, self.wedge(other), 0, 1) / 2
        # The result is named "u x v" only for a default metric:
        if (default_metric and self._name is not None and
            other._name is not None):
            resu._name = "{} x {}".format(self._name, other._name)
            resu._latex_name = "{" + self._latex_name + r"}\times{" + \
                               other._latex_name + "}"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    cross = cross_product

#******************************************************************************

class VectorFieldParal(FiniteRankFreeModuleElement, MultivectorFieldParal,
                       VectorField):
    r"""
    Vector field along a differentiable manifold, with values on a
    parallelizable manifold.

    An instance of this class is a vector field along a differentiable
    manifold `U` with values on a parallelizable manifold `M`, via a
    differentiable map `\Phi: U \rightarrow M`. More precisely, given
    a differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    a *vector field along* `U` *with values on* `M` is a differentiable map

    .. MATH::

        v:\ U  \longrightarrow TM

    (`TM` being the tangent bundle of `M`) such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M.

    The standard case of vector fields *on* a differentiable manifold
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases
    are `\Phi` being an immersion and `\Phi` being a curve in `M` (`U`
    is then an open interval of `\RR`).

    .. NOTE::

        If `M` is not parallelizable, then
        :class:`~sage.manifolds.differentiable.vectorfield.VectorField`
        *must* be used instead.

    INPUT:

    - ``vector_field_module`` -- free module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `M\supset\Phi(U)`
    - ``name`` -- (default: ``None``) name given to the vector field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the vector
      field; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A vector field on a parallelizable 3-dimensional manifold::

        sage: M = Manifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart()
        sage: v = M.vector_field(name='V') ; v
        Vector field V on the 3-dimensional differentiable manifold M
        sage: latex(v)
        V

    Vector fields are considered as elements of a module over the ring
    (algebra) of scalar fields on `M`::

        sage: v.parent()
        Free module X(M) of vector fields on the 3-dimensional differentiable
         manifold M
        sage: v.parent().base_ring()
        Algebra of differentiable scalar fields on the 3-dimensional
         differentiable manifold M
        sage: v.parent() is M.vector_field_module()
        True

    A vector field is a tensor field of rank 1 and of type `(1,0)`::

        sage: v.tensor_rank()
        1
        sage: v.tensor_type()
        (1, 0)

    Components of a vector field with respect to a given frame::

        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: v[0], v[1], v[2] = (1+y, 4*x*z, 9)  # components on M's default frame (e)
        sage: v.comp()
        1-index components w.r.t. Vector frame (M, (e_0,e_1,e_2))

    The totality of the components are accessed via the operator ``[:]``::

        sage: v[:] = (1+y, 4*x*z, 9)
        sage: v[:]
        [y + 1, 4*x*z, 9]

    The components are also read on the expansion on the frame ``e``,
    as provided by the method
    :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.display`::

        sage: v.display()  # expansion in the default frame
        V = (y + 1) e_0 + 4*x*z e_1 + 9 e_2

    A subset of the components can be accessed by using slice notation::

        sage: v[1:] = (-2, -x*y)
        sage: v[:]
        [y + 1, -2, -x*y]
        sage: v[:2]
        [y + 1, -2]

    Components in another frame::

        sage: f = M.vector_frame('f')
        sage: for i in range(3):
        ....:     v.set_comp(f)[i] = (i+1)**3 * c_xyz[i]
        sage: v.comp(f)[2]
        27*z
        sage: v[f, 2]  # equivalent to above
        27*z
        sage: v.display(f)
        V = x f_0 + 8*y f_1 + 27*z f_2

    One can set the components at the vector definition::

        sage: v = M.vector_field(1+y, 4*x*z, 9, name='V')
        sage: v.display()
        V = (y + 1) e_0 + 4*x*z e_1 + 9 e_2

    If the components regard a vector frame different from the default one,
    the vector frame has to be specified via the argument ``frame``::

        sage: v = M.vector_field(x, 8*y, 27*z, frame=f, name='V')
        sage: v.display(f)
        V = x f_0 + 8*y f_1 + 27*z f_2

    For providing the components in various frames, one may use a dictionary::

        sage: v = M.vector_field({e: [1+y, -2, -x*y], f: [x, 8*y, 27*z]},
        ....:                    name='V')
        sage: v.display(e)
        V = (y + 1) e_0 - 2 e_1 - x*y e_2
        sage: v.display(f)
        V = x f_0 + 8*y f_1 + 27*z f_2

    It is also possible to construct a vector field from a vector of symbolic
    expressions (or any other iterable)::

        sage: v = M.vector_field(vector([1+y, 4*x*z, 9]), name='V')
        sage: v.display()
        V = (y + 1) e_0 + 4*x*z e_1 + 9 e_2

    The range of the indices depends on the convention set for the manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: v = M.vector_field(1+y, 4*x*z, 9, name='V')
        sage: v[0]
        Traceback (most recent call last):
        ...
        IndexError: index out of range: 0 not in [1, 3]
        sage: v[1]  # OK
        y + 1

    A vector field acts on scalar fields (derivation along the vector field)::

        sage: M = Manifold(2, 'M')
        sage: c_cart.<x,y> = M.chart()
        sage: f = M.scalar_field(x*y^2, name='f')
        sage: v = M.vector_field(-y, x, name='v')
        sage: v.display()
        v = -y ∂/∂x + x ∂/∂y
        sage: v(f)
        Scalar field v(f) on the 2-dimensional differentiable manifold M
        sage: v(f).expr()
        2*x^2*y - y^3
        sage: latex(v(f))
        v\left(f\right)

    Example of a vector field associated with a non-trivial map `\Phi`;
    a vector field along a curve in `M`::

        sage: R = Manifold(1, 'R')
        sage: T.<t> = R.chart()  # canonical chart on R
        sage: Phi = R.diff_map(M, [cos(t), sin(t)], name='Phi') ; Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold R
         to the 2-dimensional differentiable manifold M
        sage: Phi.display()
        Phi: R → M
           t ↦ (x, y) = (cos(t), sin(t))
        sage: w = R.vector_field(-sin(t), cos(t), dest_map=Phi, name='w') ; w
        Vector field w along the 1-dimensional differentiable manifold R with
         values on the 2-dimensional differentiable manifold M
        sage: w.parent()
        Free module X(R,Phi) of vector fields along the 1-dimensional
         differentiable manifold R mapped into the 2-dimensional differentiable
         manifold M
        sage: w.display()
        w = -sin(t) ∂/∂x + cos(t) ∂/∂y

    Value at a given point::

        sage: p = R((0,), name='p') ; p
        Point p on the 1-dimensional differentiable manifold R
        sage: w.at(p)
        Tangent vector w at Point Phi(p) on the 2-dimensional differentiable
         manifold M
        sage: w.at(p).display()
        w = ∂/∂y
        sage: w.at(p) == v.at(Phi(p))
        True
    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        r"""
        Construct a vector field with values on a parallelizable manifold.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``VectorFieldParal``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()  # the parent
            sage: v = XM.element_class(XM, name='v'); v
            Vector field v on the 2-dimensional differentiable manifold M
            sage: v[:] = (-y, x)
            sage: v.display()
            v = -y ∂/∂x + x ∂/∂y
            sage: TestSuite(v).run()

        Construction via ``DifferentiableManifold.vector_field``::

            sage: u = M.vector_field(1+x, 1-y, name='u'); u
            Vector field u on the 2-dimensional differentiable manifold M
            sage: type(u) == type(v)
            True
            sage: u.parent() is v.parent()
            True
            sage: TestSuite(u).run()
        """
        FiniteRankFreeModuleElement.__init__(self, vector_field_module,
                                             name=name, latex_name=latex_name)
        # MultivectorFieldParal attributes:
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # VectorField attributes:
        self._vmodule = vector_field_module
        # Initialization of derived quantities:
        MultivectorFieldParal._init_derived(self)
        VectorField._init_derived(self)
        # Initialization of list of quantities depending on self:
        self._init_dependencies()

    def _repr_(self) :
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: v._repr_()
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: repr(v)  # indirect doctest
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: v  # indirect doctest
            Vector field v on the 2-dimensional differentiable manifold M
        """
        return VectorField._repr_(self)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same module.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: u = v._new_instance(); u
            Vector field on the 2-dimensional differentiable manifold M
            sage: u.parent() is v.parent()
            True
        """
        return type(self)(self._fmodule)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities.

        INPUT:

        - ``del_restrictions`` -- boolean (default: ``True``); determines whether
          the restrictions of ``self`` to subdomains are deleted

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: v._del_derived()
        """
        MultivectorFieldParal._del_derived(self,
                                           del_restrictions=del_restrictions)
        VectorField._del_derived(self)
        self._del_dependencies()

    def __call__(self, scalar):
        r"""
        Action on a scalar field.

        INPUT:

        - ``scalar`` -- scalar field `f`

        OUTPUT:

        - scalar field representing the derivative of `f` along the vector
          field, i.e. `v^i \frac{\partial f}{\partial x^i}`

        EXAMPLES:

        Action of a vector field on a scalar field on a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_cart.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: v = M.vector_field(-y, x)
            sage: v(f)
            Scalar field on the 2-dimensional differentiable manifold M
            sage: v(f).display()
            M → ℝ
            (x, y) ↦ 2*x^2*y - y^3
        """
        # This method enforces VectorField.__call__
        # instead of FiniteRankFreeModuleElement.__call__, which would have
        # been inheritated otherwise
        return VectorField.__call__(self, scalar)
