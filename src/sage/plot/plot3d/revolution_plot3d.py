"""
Surfaces of revolution

AUTHORS:

- Oscar Gerardo Lazo Arjona (2010): initial version.
"""

#*****************************************************************************
#       Copyright (C) 2010 Oscar Gerardo Lazo Arjona algebraicamente@gmail.com
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
from sage.misc.decorators import rename_keyword

from sage.plot.plot3d.parametric_plot3d import parametric_plot3d

@rename_keyword(alpha='opacity')
def revolution_plot3d(curve,trange,phirange=None,parallel_axis='z',axis=(0,0),print_vector=False,show_curve=False,**kwds):
    r"""
    Return a plot of a revolved curve.

    There are three ways to call this function:

    - ``revolution_plot3d(f,trange)`` where `f` is a function located in the `x z` plane.

    - ``revolution_plot3d((f_x,f_z),trange)`` where `(f_x,f_z)` is a parametric curve on the `x z` plane.

    - ``revolution_plot3d((f_x,f_y,f_z),trange)`` where `(f_x,f_y,f_z)` can be any parametric curve.

    INPUT:

    - ``curve`` - A curve to be revolved, specified as a function, a 2-tuple or a 3-tuple.

    - ``trange`` - A 3-tuple `(t,t_{\min},t_{\max})` where t is the independent variable of the curve.

    - ``phirange`` - A 2-tuple of the form `(\phi_{\min},\phi_{\max})`, (default `(0,\pi)`) that specifies the angle in which the curve is to be revolved.

    - ``parallel_axis`` - A string (Either 'x', 'y', or 'z') that specifies the coordinate axis parallel to the revolution axis.

    - ``axis`` - A 2-tuple that specifies the position of the revolution axis. If parallel is:

        - 'z' - then axis is the point in which the revolution axis intersects the  `x y` plane.

        - 'x' - then axis is the point in which the revolution axis intersects the  `y z` plane.

        - 'y' - then axis is the point in which the revolution axis intersects the `x z` plane.

    - ``print_vector`` - If True, the parametrization of the surface of revolution will be printed.

    - ``show_curve`` - If True, the curve will be displayed.


    EXAMPLES:

    Let's revolve a simple function around different axes::

        sage: u = var('u')
        sage: f = u^2
        sage: revolution_plot3d(f, (u,0,2), show_curve=True, opacity=0.7).show(aspect_ratio=(1,1,1))

    .. PLOT::

        u = var('u')
        f = u**2
        P = revolution_plot3d(f, (u,0,2), show_curve=True, opacity=0.7).plot()
        sphinx_plot(P)

    If we move slightly the axis, we get a goblet-like surface::

        sage: revolution_plot3d(f, (u,0,2), axis=(1,0.2), show_curve=True, opacity=0.5).show(aspect_ratio=(1,1,1))

    .. PLOT::

        u = var('u')
        f = u**2
        P = revolution_plot3d(f, (u,0,2), axis=(1,0.2), show_curve=True, opacity=0.5).plot()
        sphinx_plot(P)

    A common problem in calculus books, find the volume within the following revolution solid::

        sage: line = u
        sage: parabola = u^2
        sage: sur1 = revolution_plot3d(line, (u,0,1), opacity=0.5, rgbcolor=(1,0.5,0), show_curve=True, parallel_axis='x')
        sage: sur2 = revolution_plot3d(parabola, (u,0,1), opacity=0.5, rgbcolor=(0,1,0), show_curve=True, parallel_axis='x')
        sage: (sur1+sur2).show()

    .. PLOT::

        u = var('u')
        line = u
        parabola = u**2
        sur1 = revolution_plot3d(line, (u,0,1), opacity=0.5, rgbcolor=(1,0.5,0), show_curve=True, parallel_axis='x')
        sur2 = revolution_plot3d(parabola, (u,0,1), opacity=0.5, rgbcolor=(0,1,0), show_curve=True, parallel_axis='x')
        P = sur1 + sur2
        sphinx_plot(P)

    Now let's revolve a parametrically defined circle. We can play with the topology of the surface by changing the axis,
    an axis in `(0,0)` (as the previous one) will produce a sphere-like surface::

        sage: u = var('u')
        sage: circle = (cos(u), sin(u))
        sage: revolution_plot3d(circle, (u,0,2*pi), axis=(0,0), show_curve=True, opacity=0.5).show(aspect_ratio=(1,1,1))

    .. PLOT::

        u = var('u')
        circle = (cos(u), sin(u))
        P = revolution_plot3d(circle, (u,0,2*pi), axis=(0,0), show_curve=True, opacity=0.5)
        sphinx_plot(P)

    An axis on `(0,y)` will produce a cylinder-like surface::

        sage: revolution_plot3d(circle, (u,0,2*pi), axis=(0,2), show_curve=True, opacity=0.5).show(aspect_ratio=(1,1,1))

    .. PLOT::

        u = var('u')
        circle = (cos(u), sin(u))
        P = revolution_plot3d(circle, (u,0,2*pi), axis=(0,2), show_curve=True, opacity=0.5)
        sphinx_plot(P)

    And any other axis will produce a torus-like surface::

        sage: revolution_plot3d(circle, (u,0,2*pi), axis=(2,0), show_curve=True, opacity=0.5).show(aspect_ratio=(1,1,1))

    .. PLOT::

        u = var('u')
        circle = (cos(u), sin(u))
        P = revolution_plot3d(circle, (u,0,2*pi), axis=(2,0), show_curve=True, opacity=0.5)
        sphinx_plot(P)

    Now, we can get another goblet-like surface by revolving a curve in 3d::

        sage: u = var('u')
        sage: curve = (u, cos(4*u), u^2)
        sage: P = revolution_plot3d(curve, (u,0,2), show_curve=True, parallel_axis='z',axis=(1,.2), opacity=0.5)
        sage: P.show(aspect_ratio=(1,1,1))

    .. PLOT::

        u = var('u')
        curve = (u, cos(4*u), u**2)
        P = revolution_plot3d(curve, (u,0,2), show_curve=True, parallel_axis='z', axis=(1,.2), opacity=0.5)
        sphinx_plot(P)

    A curvy curve with only a quarter turn::

        sage: u = var('u')
        sage: curve = (sin(3*u), .8*cos(4*u), cos(u))
        sage: revolution_plot3d(curve, (u,0,pi), (0,pi/2), show_curve=True, parallel_axis='z', opacity=0.5).show(aspect_ratio=(1,1,1),frame=False)

    .. PLOT::

        u = var('u')
        curve = (sin(3*u), .8*cos(4*u), cos(u))
        P = revolution_plot3d(curve, (u,0,pi), (0,pi/2), show_curve=True, parallel_axis='z', opacity=0.5)
        sphinx_plot(P)

    One can also color the surface using a coloring function of two
    parameters and a colormap as follows. Note that the coloring
    function must take values in the interval [0,1]. ::

        sage: u, phi = var('u,phi')
        sage: def cf(u,phi): return sin(phi+u) ^ 2
        sage: curve = (1+u^2/4, 0, u)
        sage: revolution_plot3d(curve, (u,-2,2), (0,2*pi), parallel_axis='z', color=(cf, colormaps.PiYG)).show(aspect_ratio=(1,1,1))

    .. PLOT::

        u, phi = var('u,phi')
        def cf(u, phi): return sin(phi+u) ** 2
        curve = (1+u**2/4, 0, u)
        P = revolution_plot3d(curve, (u,-2,2), (0,2*pi), parallel_axis='z', color=(cf, colormaps.PiYG))
        sphinx_plot(P)

    The first parameter of the coloring function will be identified with the
    parameter of the curve, and the second with the angle parameter.

    .. WARNING::

        This kind of coloring using a colormap can be visualized using
        Jmol, Tachyon (option ``viewer='tachyon'``) and Canvas3D
        (option ``viewer='canvas3d'`` in the notebook).

    Another colored example, illustrating that one can use (colormap, color function) instead of (color function, colormap)::

        sage: u, phi = var('u,phi')
        sage: def cf(u, phi): return float(2 * u / pi) % 1
        sage: curve = (sin(u), 0, u)
        sage: revolution_plot3d(curve, (u,0,pi), (0,2*pi), parallel_axis
        ....: ='z', color=(colormaps.brg, cf)).show(aspect_ratio=1)

    .. PLOT::

        u, phi = var('u,phi')
        def cf(u, phi): return float(2 * u / pi) % 1
        curve = (sin(u), 0, u)
        P = revolution_plot3d(curve, (u,0,pi), (0,2*pi), parallel_axis='z', color=(colormaps.brg, cf))
        sphinx_plot(P)
    """
    from sage.symbolic.ring import SR
    from sage.symbolic.constants import pi
    from sage.misc.functional import sqrt
    from sage.functions.trig import sin
    from sage.functions.trig import cos
    from sage.functions.trig import atan2

    if parallel_axis not in ['x', 'y', 'z']:
        raise ValueError("parallel_axis must be either 'x', 'y', or 'z'.")

    vart = trange[0]

    if str(vart) == 'phi':
        phi = SR.var('fi')
    else:
        phi = SR.var('phi')

    if phirange is None:  # this if-else provides a phirange
        phirange = (phi, 0, 2 * pi)
    elif len(phirange) == 3:
        phi = phirange[0]
    else:
        phirange = (phi, phirange[0], phirange[1])
        
    if isinstance(curve, (tuple, list)):
        #this if-else provides a vector v to be plotted
        #if curve is a tuple or a list of length 2, it is interpreted as a parametric curve
        #in the x-z plane.
        #if it is of length 3 it is interpreted as a parametric curve in 3d space

        if len(curve) == 2:
            x = curve[0]
            y = 0
            z = curve[1]
        elif len(curve) == 3:
            x = curve[0]
            y = curve[1]
            z = curve[2]
    else:
        x = vart
        y = 0
        z = curve

    phase = 0
    if parallel_axis == 'z':
        x0 = axis[0]
        y0 = axis[1]
        # (0,0) must be handled separately for the phase value
        if x0 != 0 or y0 != 0:
            phase = atan2(y - y0, x - x0)
        R = sqrt((x-x0)**2 + (y-y0)**2)
        v = (R*cos(phi+phase)+x0, R*sin(phi+phase)+y0, z)
    elif parallel_axis == 'x':
        y0 = axis[0]
        z0 = axis[1]
        # (0,0) must be handled separately for the phase value
        if z0 != 0 or y0 != 0:
            phase = atan2(z - z0, y - y0)
        R = sqrt((y-y0)**2 + (z-z0)**2)
        v = (x, R*cos(phi+phase)+y0, R*sin(phi+phase)+z0)
    elif parallel_axis == 'y':
        x0 = axis[0]
        z0 = axis[1]
        # (0,0) must be handled separately for the phase value
        if z0 != 0 or x0 != 0:
            phase = atan2(z - z0, x - x0)
        R = sqrt((x-x0)**2 + (z-z0)**2)
        v = (R*cos(phi+phase)+x0, y, R*sin(phi+phase)+z0)

    if print_vector:
        print(v)

    if show_curve:
        curveplot = parametric_plot3d((x, y, z), trange, thickness=2,
                                      rgbcolor=(1, 0, 0))
        return parametric_plot3d(v, trange, phirange, **kwds) + curveplot

    return parametric_plot3d(v, trange, phirange, **kwds)
