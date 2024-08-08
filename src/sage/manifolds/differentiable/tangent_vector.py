r"""
Tangent Vectors

The class :class:`TangentVector` implements tangent vectors to a differentiable
manifold.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- Chap. 3 of [Lee2013]_
"""

# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement
from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
from .scalarfield import DiffScalarField
from sage.misc.decorators import options


class TangentVector(FiniteRankFreeModuleElement):
    r"""
    Tangent vector to a differentiable manifold at a given point.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.tangent_space.TangentSpace`;
      the tangent space to which the vector belongs
    - ``name`` -- (default: ``None``) string; symbol given to the vector
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the vector; if ``None``, ``name`` will be used

    EXAMPLES:

    A tangent vector `v` on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: p = M.point((2,3), name='p')
        sage: Tp = M.tangent_space(p)
        sage: v = Tp((-2,1), name='v') ; v
        Tangent vector v at Point p on the 2-dimensional differentiable
         manifold M
        sage: v.display()
        v = -2 ∂/∂x + ∂/∂y
        sage: v.parent()
        Tangent space at Point p on the 2-dimensional differentiable manifold M
        sage: v in Tp
        True

    Tangent vectors can also be constructed via the manifold method
    :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.tangent_vector`::

        sage: v = M.tangent_vector(p, (-2, 1), name='v'); v
        Tangent vector v at Point p on the 2-dimensional differentiable
         manifold M
        sage: v.display()
        v = -2 ∂/∂x + ∂/∂y

    or via the method
    :meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.at`
    of vector fields::

        sage: vf = M.vector_field(x - 4*y/3, (x-y)^2, name='v')
        sage: v = vf.at(p); v
        Tangent vector v at Point p on the 2-dimensional differentiable
         manifold M
        sage: v.display()
        v = -2 ∂/∂x + ∂/∂y

    By definition, a tangent vector at `p\in M` is a *derivation at* `p` on
    the space `C^\infty(M)` of smooth scalar fields on `M`. Indeed  let us
    consider a generic scalar field `f`::

        sage: f = M.scalar_field(function('F')(x,y), name='f')
        sage: f.display()
        f: M → ℝ
           (x, y) ↦ F(x, y)

    The tangent vector `v` maps `f` to the real number
    `v^i \left. \frac{\partial F}{\partial x^i} \right|_p`::

        sage: v(f)
        -2*D[0](F)(2, 3) + D[1](F)(2, 3)
        sage: vdf(x, y) = v[0]*diff(f.expr(), x) + v[1]*diff(f.expr(), y)
        sage: X(p)
        (2, 3)
        sage: bool( v(f) == vdf(*X(p)) )
        True

    and if `g` is a second scalar field on `M`::

        sage: g = M.scalar_field(function('G')(x,y), name='g')

    then the product `f g` is also a scalar field on `M`::

        sage: (f*g).display()
        f*g: M → ℝ
           (x, y) ↦ F(x, y)*G(x, y)

    and we have the derivation law `v(f g) = v(f) g(p) + f(p) v(g)`::

        sage: bool( v(f*g) == v(f)*g(p) + f(p)*v(g) )
        True

    .. SEEALSO::

        :class:`~sage.tensor.modules.free_module_element.FiniteRankFreeModuleElement`
        for more documentation.
    """
    def __init__(self, parent, name=None, latex_name=None):
        r"""
        Construct a tangent vector.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp.element_class(Tp, name='v') ; v
            Tangent vector v at Point p on the 2-dimensional differentiable
             manifold M
            sage: v[:] = 5, -3/2
            sage: TestSuite(v).run()
        """
        FiniteRankFreeModuleElement.__init__(self, parent, name=name,
                                             latex_name=latex_name)
        # Extra data (with respect to FiniteRankFreeModuleElement):
        self._point = parent._point

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp([-3,2], name='v')
            sage: v._repr_()
            'Tangent vector v at Point p on the 2-dimensional differentiable manifold M'
            sage: repr(v)  # indirect doctest
            'Tangent vector v at Point p on the 2-dimensional differentiable manifold M'
        """
        from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace
        if isinstance(self._point.parent(), EuclideanSpace):
            desc = "Vector"
        else:
            desc = "Tangent vector"
        if self._name:
            desc += " " + str(self._name)
        desc += " at " + str(self._point)
        return desc

    @options(scale=1)
    def plot(self, chart=None, ambient_coords=None, mapping=None,
             color='blue', print_label=True, label=None, label_color=None,
             fontsize=10, label_offset=0.1, parameters=None, **extra_options):
        r"""
        Plot the vector in a Cartesian graph based on the coordinates of some
        ambient chart.

        The vector is drawn in terms of two (2D graphics) or three (3D graphics)
        coordinates of a given chart, called hereafter the *ambient chart*.
        The vector's base point `p` (or its image `\Phi(p)` by some
        differentiable mapping `\Phi`) must lie in the ambient chart's domain.
        If `\Phi` is different from the identity mapping, the vector
        actually depicted is `\mathrm{d}\Phi_p(v)`, where `v` is the current
        vector (``self``) (see the example of a vector tangent to the
        2-sphere below, where `\Phi: S^2 \to \RR^3`).

        INPUT:

        - ``chart`` -- (default: ``None``) the ambient chart (see above); if
          ``None``, it is set to the default chart of the open set containing
          the point at which the vector (or the vector image via the
          differential `\mathrm{d}\Phi_p` of ``mapping``) is defined

        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2
          or 3 coordinates of the ambient chart in terms of which the plot
          is performed; if ``None``, all the coordinates of the ambient
          chart are considered

        - ``mapping`` -- (default: ``None``)
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`;
          differentiable mapping `\Phi` providing the link between the
          point `p` at which the vector is defined and the ambient chart
          ``chart``: the domain of ``chart`` must contain `\Phi(p)`;
          if ``None``, the identity mapping is assumed

        - ``scale`` -- (default: 1) value by which the length of the arrow
          representing the vector is multiplied

        - ``color`` -- (default: ``'blue'``) color of the arrow representing the
          vector

        - ``print_label`` -- boolean (default: ``True``); determines whether a
          label is printed next to the arrow representing the vector

        - ``label`` -- (string; default: ``None``) label printed next to the
          arrow representing the vector; if ``None``, the vector's symbol is
          used, if any

        - ``label_color`` -- (default: ``None``) color to print the label;
          if ``None``, the value of ``color`` is used

        - ``fontsize`` -- (default: 10) size of the font used to print the
          label

        - ``label_offset`` -- (default: 0.1) determines the separation between
          the vector arrow and the label

        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the coordinate expression
          of ``self`` (see example below)

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

        Vector tangent to a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((2, 1), name='v') ; v
            Tangent vector v at Point p on the 2-dimensional differentiable
             manifold M

        Plot of the vector alone (arrow + label)::

            sage: v.plot()                                                              # needs sage.plot
            Graphics object consisting of 2 graphics primitives

        Plot atop of the chart grid::

            sage: X.plot() + v.plot()                                                   # needs sage.plot
            Graphics object consisting of 20 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            p = M((2,2), name='p'); Tp = M.tangent_space(p)
            v = Tp((2, 1), name='v')
            g = X.plot() + v.plot()
            sphinx_plot(g)

        Plots with various options::

            sage: X.plot() + v.plot(color='green', scale=2, label='V')                  # needs sage.plot
            Graphics object consisting of 20 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            p = M((2,2), name='p'); Tp = M.tangent_space(p)
            v = Tp((2, 1), name='v')
            g = X.plot() + v.plot(color='green', scale=2, label='V')
            sphinx_plot(g)

        ::

            sage: X.plot() + v.plot(print_label=False)                                  # needs sage.plot
            Graphics object consisting of 19 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            p = M((2,2), name='p'); Tp = M.tangent_space(p)
            v = Tp((2, 1), name='v')
            g = X.plot() + v.plot(print_label=False)
            sphinx_plot(g)

        ::

            sage: X.plot() + v.plot(color='green', label_color='black',                 # needs sage.plot
            ....:                   fontsize=20, label_offset=0.2)
            Graphics object consisting of 20 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            p = M((2,2), name='p'); Tp = M.tangent_space(p)
            v = Tp((2, 1), name='v')
            g = X.plot() + v.plot(color='green', label_color='black', fontsize=20, label_offset=0.2)
            sphinx_plot(g)

        ::

            sage: X.plot() + v.plot(linestyle=':', width=4, arrowsize=8,                # needs sage.plot
            ....:                   fontsize=20)
            Graphics object consisting of 20 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            p = M((2,2), name='p'); Tp = M.tangent_space(p)
            v = Tp((2, 1), name='v')
            g = X.plot() + v.plot(linestyle=':', width=4, arrowsize=8, fontsize=20)
            sphinx_plot(g)

        Plot with specific values of some free parameters::

            sage: var('a b')
            (a, b)
            sage: v = Tp((1+a, -b^2), name='v') ; v.display()
            v = (a + 1) ∂/∂x - b^2 ∂/∂y
            sage: X.plot() + v.plot(parameters={a: -2, b: 3})                           # needs sage.plot
            Graphics object consisting of 20 graphics primitives

        Special case of the zero vector::

            sage: v = Tp.zero() ; v
            Tangent vector zero at Point p on the 2-dimensional differentiable
             manifold M
            sage: X.plot() + v.plot()                                                   # needs sage.plot
            Graphics object consisting of 19 graphics primitives

        Vector tangent to a 4-dimensional manifold::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: p = M((0,1,2,3), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((5,4,3,2), name='v') ; v
            Tangent vector v at Point p on the 4-dimensional differentiable
             manifold M

        We cannot make a 4D plot directly::

            sage: v.plot()                                                              # needs sage.plot
            Traceback (most recent call last):
            ...
            ValueError: the number of coordinates involved in the plot must
             be either 2 or 3, not 4

        Rather, we have to select some chart coordinates for the plot, via
        the argument ``ambient_coords``. For instance, for a 2-dimensional
        plot in terms of the coordinates `(x, y)`::

            sage: v.plot(ambient_coords=(x,y))                                          # needs sage.plot
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            M = Manifold(4, 'M')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            p = M((0,1,2,3), name='p'); Tp = M.tangent_space(p)
            v = Tp((5,4,3,2), name='v')
            g = X.plot(ambient_coords=(x,y)) + v.plot(ambient_coords=(x,y))
            sphinx_plot(g)

        This plot involves only the components `v^x` and `v^y` of `v`.
        Similarly, for a 3-dimensional plot in terms of the coordinates
        `(t, x, y)`::

            sage: g = v.plot(ambient_coords=(t,x,z))                                    # needs sage.plot
            sage: print(g)                                                              # needs sage.plot
            Graphics3d Object

        This plot involves only the components `v^t`,  `v^x` and `v^z` of `v`.
        A nice 3D view atop the coordinate grid is obtained via::

            sage: (X.plot(ambient_coords=(t,x,z))  # long time                          # needs sage.plot
            ....:  + v.plot(ambient_coords=(t,x,z),
            ....:           label_offset=0.5, width=6))
            Graphics3d Object

        .. PLOT::

            M = Manifold(4, 'M')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            p = M((0,1,2,3), name='p'); Tp = M.tangent_space(p)
            v = Tp((5,4,3,2), name='v')
            g = X.plot(ambient_coords=(t,x,z)) + v.plot(ambient_coords=(t,x,z),
                       label_offset=0.5, width=6)
            sphinx_plot(g)

        An example of plot via a differential mapping: plot of a vector tangent
        to a 2-sphere viewed in `\RR^3`::

            sage: S2 = Manifold(2, 'S^2')
            sage: U = S2.open_subset('U') # the open set covered by spherical coord.
            sage: XS.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: R3 = Manifold(3, 'R^3')
            sage: X3.<x,y,z> = R3.chart()
            sage: F = S2.diff_map(R3, {(XS, X3): [sin(th)*cos(ph),
            ....:                                 sin(th)*sin(ph),
            ....:                                 cos(th)]}, name='F')
            sage: F.display() # the standard embedding of S^2 into R^3
            F: S^2 → R^3
            on U: (th, ph) ↦ (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: p = U.point((pi/4, 7*pi/4), name='p')
            sage: v = XS.frame()[1].at(p) ; v  # the coordinate vector ∂/∂phi at p
            Tangent vector ∂/∂ph at Point p on the 2-dimensional differentiable
             manifold S^2
            sage: graph_v = v.plot(mapping=F)                                           # needs sage.plot
            sage: graph_S2 = XS.plot(chart=X3, mapping=F, number_values=9)  # long time, needs sage.plot
            sage: graph_v + graph_S2                                        # long time, needs sage.plot
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
            p = U.point((pi/4, 7*pi/4), name='p')
            v = XS.frame()[1].at(p)
            graph_v = v.plot(mapping=F)
            graph_S2 = XS.plot(chart=X3, mapping=F, number_values=9)
            sphinx_plot(graph_v + graph_S2)
        """
        from sage.plot.arrow import arrow2d
        from sage.plot.text import text
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.shapes import arrow3d
        from sage.plot.plot3d.shapes2 import text3d
        from sage.misc.functional import numerical_approx
        from sage.manifolds.differentiable.chart import DiffChart

        scale = extra_options.pop("scale")

        #
        # The "effective" vector to be plotted
        #
        if mapping is None:
            eff_vector = self
            base_point = self._point
        else:
            #!# check
            # For efficiency, the method FiniteRankFreeModuleMorphism._call_()
            # is called instead of FiniteRankFreeModuleMorphism.__call__()
            eff_vector = mapping.differential(self._point)._call_(self)
            base_point = mapping(self._point)
        #
        # The chart w.r.t. which the vector is plotted
        #
        if chart is None:
            chart = base_point.parent().default_chart()
        elif not isinstance(chart, DiffChart):
            raise TypeError("{} is not a chart".format(chart))
        #
        # Coordinates of the above chart w.r.t. which the vector is plotted
        #
        if ambient_coords is None:
            ambient_coords = chart[:]  # all chart coordinates are used
        n_pc = len(ambient_coords)
        if n_pc != 2 and n_pc != 3:
            raise ValueError("the number of coordinates involved in the " +
                             "plot must be either 2 or 3, not {}".format(n_pc))
        # indices coordinates involved in the plot:
        ind_pc = [chart[:].index(pc) for pc in ambient_coords]
        #
        # Components of the vector w.r.t. the chart frame
        #
        basis = chart.frame().at(base_point)
        vcomp = eff_vector.comp(basis=basis)[:]
        xp = base_point.coord(chart=chart)
        #
        # The arrow
        #
        resu = Graphics()
        if parameters is None:
            coord_tail = [numerical_approx(xp[i]) for i in ind_pc]
            coord_head = [numerical_approx(xp[i] + scale*vcomp[i])
                          for i in ind_pc]
        else:
            coord_tail = [numerical_approx(xp[i].substitute(parameters))
                          for i in ind_pc]
            coord_head = [numerical_approx(
                           (xp[i] + scale*vcomp[i]).substitute(parameters))
                          for i in ind_pc]
        if coord_head != coord_tail:
            if n_pc == 2:
                resu += arrow2d(tailpoint=coord_tail, headpoint=coord_head,
                                color=color, **extra_options)
            else:
                resu += arrow3d(coord_tail, coord_head, color=color,
                                **extra_options)
        #
        # The label
        #
        if print_label:
            if label is None:
                if n_pc == 2 and self._latex_name is not None:
                    label = r'$' + self._latex_name + r'$'
                if n_pc == 3 and self._name is not None:
                    label = self._name
            if label is not None:
                xlab = [xh + label_offset for xh in coord_head]
                if label_color is None:
                    label_color = color
                if n_pc == 2:
                    resu += text(label, xlab, fontsize=fontsize,
                                 color=label_color)
                else:
                    resu += text3d(label, xlab, fontsize=fontsize,
                                   color=label_color)
        return resu

    def __call__(self, f):
        r"""
        Action on a scalar field (as a derivation) or on a linear form.

        INPUT:

        - ``f`` -- either a scalar field on the manifold of which ``self`` is
          defined or a linear form in the same tangent space as ``self``

        OUTPUT:

        - scalar (element of the manifold base field)

        EXAMPLES:

        Let us consider a tangent vector on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: p = M((2, 3), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((-1, 2))
            sage: v.display()
            -∂/∂x + 2 ∂/∂y

        The action of `v` on a scalar field `f`::

            sage: f = M.scalar_field(x*y^2, name='f')
            sage: v(f)
            15

        Check of the formula `v(f) = v^i \frac{\partial f}{\partial x^i}|_p`::

            sage: vdf(x, y) = v[1]*diff(f.expr(), x) + v[2]*diff(f.expr(), y)
            sage: vdf
            (x, y) |--> 4*x*y - y^2
            sage: X(p)
            (2, 3)
            sage: bool( v(f) == vdf(*X(p)) )
            True

        Case of a generic scalar field::

            sage: f = M.scalar_field(function('F')(x,y), name='f')
            sage: v(f)
            -D[0](F)(2, 3) + 2*D[1](F)(2, 3)

        Action of a tangent vector on a linear form on the same tangent space::

            sage: omega = Tp.linear_form()
            sage: omega[:] = 4, 1
            sage: omega.display()
            4 dx + dy
            sage: v(omega)
            -2

        Checks of he formula `v(\omega) = v^i \omega_i`::

            sage: bool( v(omega) == v[1]*omega[1] + v[2]*omega[2] )
            True

        Another check::

            sage: bool( v(omega) == omega(v) )
            True
        """
        if isinstance(f, FreeModuleAltForm):
            # Case of self acting on a linear form
            if f.tensor_type() != (0, 1):
                raise TypeError("the argument of __call__ must be a linear form, "
                                "not {}".format(f))
            return f(self)
        if not isinstance(f, DiffScalarField):
            raise TypeError("the argument of __call__ must be either a linear "
                            "form or a scalar field, not {}".format(f))
        # Case of self acting on a scalar field
        return f.differential().at(self._point)(self)
