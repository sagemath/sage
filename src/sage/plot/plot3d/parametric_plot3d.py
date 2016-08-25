# -*- coding: utf-8 -*-
"""
Parametric Plots
"""
from __future__ import absolute_import

from .parametric_surface import ParametricSurface
from .shapes2 import line3d
from sage.arith.srange import xsrange, srange
from sage.structure.element import is_Vector

from sage.ext.fast_eval import fast_float, fast_float_constant


def parametric_plot3d(f, urange, vrange=None, plot_points="automatic",
                      boundary_style=None, **kwds):
    r"""
    Return a parametric three-dimensional space curve or surface.

    There are four ways to call this function:

    - ``parametric_plot3d([f_x, f_y, f_z], (u_min,
      u_max))``:
      `f_x, f_y, f_z` are three functions and
      `u_{\min}` and `u_{\max}` are real numbers

    - ``parametric_plot3d([f_x, f_y, f_z], (u, u_min,
      u_max))``:
      `f_x, f_y, f_z` can be viewed as functions of
      `u`

    - ``parametric_plot3d([f_x, f_y, f_z], (u_min, u_max),
      (v_min, v_max))``:
      `f_x, f_y, f_z` are each functions of two variables

    - ``parametric_plot3d([f_x, f_y, f_z], (u, u_min,
      u_max), (v, v_min, v_max))``:
      `f_x, f_y, f_z` can be viewed as functions of
      `u` and `v`


    INPUT:

    - ``f`` - a 3-tuple of functions or expressions, or vector of size 3

    - ``urange`` - a 2-tuple (u_min, u_max) or a 3-tuple
      (u, u_min, u_max)

    - ``vrange`` - (optional - only used for surfaces) a
      2-tuple (v_min, v_max) or a 3-tuple (v, v_min, v_max)

    - ``plot_points`` - (default: "automatic", which is
      75 for curves and [40,40] for surfaces) initial number of sample
      points in each parameter; an integer for a curve, and a pair of
      integers for a surface.

    - ``boundary_style`` - (default: None, no boundary) a dict that describes
      how to draw the boundaries of regions by giving options that are passed
      to the line3d command.

    - ``mesh`` - bool (default: False) whether to display
      mesh grid lines

    - ``dots`` - bool (default: False) whether to display
      dots at mesh grid points

    .. note::

       #. By default for a curve any points where `f_x`,
          `f_y`, or `f_z` do not evaluate to a real number
          are skipped.

       #. Currently for a surface `f_x`, `f_y`, and
          `f_z` have to be defined everywhere. This will change.

       #. mesh and dots are not supported when using the Tachyon ray tracer
          renderer.


    EXAMPLES: We demonstrate each of the four ways to call this
    function.


    #. A space curve defined by three functions of 1 variable:

       ::

           sage: parametric_plot3d( (sin, cos, lambda u: u/10), (0, 20))
           Graphics3d Object

       Note above the lambda function, which creates a callable Python
       function that sends `u` to `u/10`.

    #. Next we draw the same plot as above, but using symbolic
       functions:

       ::

           sage: u = var('u')
           sage: parametric_plot3d( (sin(u), cos(u), u/10), (u, 0, 20))
           Graphics3d Object

    #. We draw a parametric surface using 3 Python functions (defined
       using lambda):

       ::

           sage: f = (lambda u,v: cos(u), lambda u,v: sin(u)+cos(v), lambda u,v: sin(v))
           sage: parametric_plot3d(f, (0, 2*pi), (-pi, pi))
           Graphics3d Object

    #. The same surface, but where the defining functions are
       symbolic:

       ::

           sage: u, v = var('u,v')
           sage: parametric_plot3d((cos(u), sin(u) + cos(v), sin(v)), (u, 0, 2*pi), (v, -pi, pi))
           Graphics3d Object

    The surface, but with a mesh::

           sage: u, v = var('u,v')
           sage: parametric_plot3d((cos(u), sin(u) + cos(v), sin(v)), (u, 0, 2*pi), (v, -pi, pi), mesh=True)
           Graphics3d Object

    We increase the number of plot points, and make the surface green
    and transparent::

        sage: parametric_plot3d((cos(u), sin(u) + cos(v), sin(v)), (u, 0, 2*pi), (v, -pi, pi),
        ....: color='green', opacity=0.1, plot_points=[30,30])
        Graphics3d Object

    One can also color the surface using a coloring function and a colormap::

        sage: u,v = var('u,v')
        sage: def cf(u,v): return sin(u+v/2)**2
        sage: P = parametric_plot3d((cos(u), sin(u) + cos(v), sin(v)),
        ....:   (u, 0, 2*pi), (v, -pi, pi), color=(cf, colormaps.PiYG), plot_points=[60,60])
        sage: P.show(viewer='tachyon')

    Another example, a colored Möbius band::

        sage: cm = colormaps.ocean
        sage: def c(x,y): return sin(x*y)**2
        sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
        sage: MoebiusStrip(5,1,plot_points=200, color=(c,cm))
        Graphics3d Object

    Yet another colored example::

        sage: from sage.plot.plot3d.parametric_surface import ParametricSurface
        sage: cm = colormaps.autumn
        sage: def c(x,y): return sin(x*y)**2
        sage: def g(x,y): return x, y+sin(y), x**2 + y**2
        sage: ParametricSurface(g, (srange(-10,10,0.1), srange(-5,5.0,0.1)),color=(c,cm))
        Graphics3d Object

    .. WARNING::

        This kind of coloring using a colormap can be visualized
        using Jmol, Tachyon (option ``viewer='tachyon'``) and
        Canvas3D (option ``viewer='canvas3d'`` in the
        notebook).

    We call the space curve function but with polynomials instead of
    symbolic variables.

    ::

        sage: R.<t> = RDF[]
        sage: parametric_plot3d( (t, t^2, t^3), (t, 0, 3) )
        Graphics3d Object

    Next we plot the same curve, but because we use (0, 3) instead of
    (t, 0, 3), each polynomial is viewed as a callable function of one
    variable::

        sage: parametric_plot3d( (t, t^2, t^3), (0, 3) )
        Graphics3d Object

    We do a plot but mix a symbolic input, and an integer::

        sage: t = var('t')
        sage: parametric_plot3d( (1, sin(t), cos(t)), (t, 0, 3) )
        Graphics3d Object

    We specify a boundary style to show us the values of the function at its
    extrema::

        sage: u, v = var('u,v')
        sage: parametric_plot3d((cos(u), sin(u) + cos(v), sin(v)), (u, 0, pi), (v, 0, pi),
        ....:              boundary_style={"color": "black", "thickness": 2})
        Graphics3d Object

    We can plot vectors::

        sage: x,y = var('x,y')
        sage: parametric_plot3d(vector([x-y,x*y,x*cos(y)]), (x,0,2), (y,0,2))
        Graphics3d Object
        sage: t = var('t')
        sage: p = vector([1,2,3])
        sage: q = vector([2,-1,2])
        sage: parametric_plot3d(p*t+q, (t, 0, 2))
        Graphics3d Object


    Any options you would normally use to specify the appearance of a curve are
    valid as entries in the ``boundary_style`` dict.

    MANY MORE EXAMPLES:

    We plot two interlinked tori::

        sage: u, v = var('u,v')
        sage: f1 = (4+(3+cos(v))*sin(u), 4+(3+cos(v))*cos(u), 4+sin(v))
        sage: f2 = (8+(3+cos(v))*cos(u), 3+sin(v), 4+(3+cos(v))*sin(u))
        sage: p1 = parametric_plot3d(f1, (u,0,2*pi), (v,0,2*pi), texture="red")
        sage: p2 = parametric_plot3d(f2, (u,0,2*pi), (v,0,2*pi), texture="blue")
        sage: p1 + p2
        Graphics3d Object

    A cylindrical Star of David::

        sage: u,v = var('u v')
        sage: f_x = cos(u)*cos(v)*(abs(cos(3*v/4))^500 + abs(sin(3*v/4))^500)^(-1/260)*(abs(cos(4*u/4))^200 + abs(sin(4*u/4))^200)^(-1/200)
        sage: f_y = cos(u)*sin(v)*(abs(cos(3*v/4))^500 + abs(sin(3*v/4))^500)^(-1/260)*(abs(cos(4*u/4))^200 + abs(sin(4*u/4))^200)^(-1/200)
        sage: f_z = sin(u)*(abs(cos(4*u/4))^200 + abs(sin(4*u/4))^200)^(-1/200)
        sage: parametric_plot3d([f_x, f_y, f_z], (u, -pi, pi), (v, 0, 2*pi))
        Graphics3d Object

    Double heart::

        sage: u, v = var('u,v')
        sage: f_x = ( abs(v) - abs(u) - abs(tanh((1/sqrt(2))*u)/(1/sqrt(2))) + abs(tanh((1/sqrt(2))*v)/(1/sqrt(2))) )*sin(v)
        sage: f_y = ( abs(v) - abs(u) - abs(tanh((1/sqrt(2))*u)/(1/sqrt(2))) - abs(tanh((1/sqrt(2))*v)/(1/sqrt(2))) )*cos(v)
        sage: f_z = sin(u)*(abs(cos(4*u/4))^1 + abs(sin(4*u/4))^1)^(-1/1)
        sage: parametric_plot3d([f_x, f_y, f_z], (u, 0, pi), (v, -pi, pi))
        Graphics3d Object

    Heart::

        sage: u, v = var('u,v')
        sage: f_x = cos(u)*(4*sqrt(1-v^2)*sin(abs(u))^abs(u))
        sage: f_y = sin(u) *(4*sqrt(1-v^2)*sin(abs(u))^abs(u))
        sage: f_z = v
        sage: parametric_plot3d([f_x, f_y, f_z], (u, -pi, pi), (v, -1, 1), frame=False, color="red")
        Graphics3d Object

    Green bowtie::

        sage: u, v = var('u,v')
        sage: f_x = sin(u) / (sqrt(2) + sin(v))
        sage: f_y = sin(u) / (sqrt(2) + cos(v))
        sage: f_z = cos(u) / (1 + sqrt(2))
        sage: parametric_plot3d([f_x, f_y, f_z], (u, -pi, pi), (v, -pi, pi), frame=False, color="green")
        Graphics3d Object

    Boy's surface http://en.wikipedia.org/wiki/Boy's_surface

    ::

        sage: u, v = var('u,v')
        sage: fx = 2/3* (cos(u)* cos(2*v) + sqrt(2)* sin(u)* cos(v))* cos(u) / (sqrt(2) - sin(2*u)* sin(3*v))
        sage: fy = 2/3* (cos(u)* sin(2*v) - sqrt(2)* sin(u)* sin(v))* cos(u) / (sqrt(2) - sin(2*u)* sin(3*v))
        sage: fz = sqrt(2)* cos(u)* cos(u) / (sqrt(2) - sin(2*u)* sin(3*v))
        sage: parametric_plot3d([fx, fy, fz], (u, -2*pi, 2*pi), (v, 0, pi), plot_points = [90,90], frame=False, color="orange") # long time -- about 30 seconds
        Graphics3d Object

    Maeder's Owl (pretty but can't find an internet reference)::

        sage: u, v = var('u,v')
        sage: fx = v *cos(u) - 0.5* v^2 * cos(2* u)
        sage: fy = -v *sin(u) - 0.5* v^2 * sin(2* u)
        sage: fz = 4 *v^1.5 * cos(3 *u / 2) / 3
        sage: parametric_plot3d([fx, fy, fz], (u, -2*pi, 2*pi), (v, 0, 1),plot_points = [90,90], frame=False, color="purple")
        Graphics3d Object

    Bracelet::

        sage: u, v = var('u,v')
        sage: fx = (2 + 0.2*sin(2*pi*u))*sin(pi*v)
        sage: fy = 0.2*cos(2*pi*u) *3*cos(2*pi*v)
        sage: fz = (2 + 0.2*sin(2*pi*u))*cos(pi*v)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, pi/2), (v, 0, 3*pi/4), frame=False, color="gray")
        Graphics3d Object

    Green goblet

    ::

        sage: u, v = var('u,v')
        sage: fx = cos(u)*cos(2*v)
        sage: fy = sin(u)*cos(2*v)
        sage: fz = sin(v)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, pi), frame=False, color="green")
        Graphics3d Object

    Funny folded surface - with square projection::

        sage: u, v = var('u,v')
        sage: fx = cos(u)*sin(2*v)
        sage: fy = sin(u)*cos(2*v)
        sage: fz = sin(v)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi), frame=False, color="green")
        Graphics3d Object

    Surface of revolution of figure 8::

        sage: u, v = var('u,v')
        sage: fx = cos(u)*sin(2*v)
        sage: fy = sin(u)*sin(2*v)
        sage: fz = sin(v)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi), frame=False, color="green")
        Graphics3d Object

    Yellow Whitney's umbrella
    http://en.wikipedia.org/wiki/Whitney_umbrella::

        sage: u, v = var('u,v')
        sage: fx = u*v
        sage: fy = u
        sage: fz = v^2
        sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -1, 1), frame=False, color="yellow")
        Graphics3d Object

    Cross cap http://en.wikipedia.org/wiki/Cross-cap::

        sage: u, v = var('u,v')
        sage: fx = (1+cos(v))*cos(u)
        sage: fy = (1+cos(v))*sin(u)
        sage: fz = -tanh((2/3)*(u-pi))*sin(v)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi), frame=False, color="red")
        Graphics3d Object

    Twisted torus::

        sage: u, v = var('u,v')
        sage: fx = (3+sin(v)+cos(u))*cos(2*v)
        sage: fy = (3+sin(v)+cos(u))*sin(2*v)
        sage: fz = sin(u)+2*cos(v)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi), frame=False, color="red")
        Graphics3d Object

    Four intersecting discs::

        sage: u, v = var('u,v')
        sage: fx = v *cos(u) -0.5*v^2*cos(2*u)
        sage: fy = -v*sin(u) -0.5*v^2*sin(2*u)
        sage: fz = 4* v^1.5 *cos(3* u / 2) / 3
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 4*pi), (v, 0,2*pi), frame=False, color="red", opacity=0.7)
        Graphics3d Object

    Steiner surface/Roman's surface (see
    http://en.wikipedia.org/wiki/Roman_surface and
    http://en.wikipedia.org/wiki/Steiner_surface)::

        sage: u, v = var('u,v')
        sage: fx = (sin(2 * u) * cos(v) * cos(v))
        sage: fy = (sin(u) * sin(2 * v))
        sage: fz = (cos(u) * sin(2 * v))
        sage: parametric_plot3d([fx, fy, fz], (u, -pi/2, pi/2), (v, -pi/2,pi/2), frame=False, color="red")
        Graphics3d Object

    Klein bottle? (see http://en.wikipedia.org/wiki/Klein_bottle)::

        sage: u, v = var('u,v')
        sage: fx = (3*(1+sin(v)) + 2*(1-cos(v)/2)*cos(u))*cos(v)
        sage: fy = (4+2*(1-cos(v)/2)*cos(u))*sin(v)
        sage: fz = -2*(1-cos(v)/2) * sin(u)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi), frame=False, color="green")
        Graphics3d Object

    A Figure 8 embedding of the Klein bottle (see
    http://en.wikipedia.org/wiki/Klein_bottle)::

        sage: u, v = var('u,v')
        sage: fx = (2 + cos(v/2)* sin(u) - sin(v/2)* sin(2 *u))* cos(v)
        sage: fy = (2 + cos(v/2)* sin(u) - sin(v/2)* sin(2 *u))* sin(v)
        sage: fz = sin(v/2)* sin(u) + cos(v/2) *sin(2* u)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi), frame=False, color="red")
        Graphics3d Object

    Enneper's surface (see
    http://en.wikipedia.org/wiki/Enneper_surface)::

        sage: u, v = var('u,v')
        sage: fx = u -u^3/3  + u*v^2
        sage: fy = v -v^3/3  + v*u^2
        sage: fz = u^2 - v^2
        sage: parametric_plot3d([fx, fy, fz], (u, -2, 2), (v, -2, 2), frame=False, color="red")
        Graphics3d Object

    Henneberg's surface (see
    http://xahlee.org/surface/gallery_m.html)

    ::

        sage: u, v = var('u,v')
        sage: fx = 2*sinh(u)*cos(v) -(2/3)*sinh(3*u)*cos(3*v)
        sage: fy = 2*sinh(u)*sin(v) +(2/3)*sinh(3*u)*sin(3*v)
        sage: fz = 2*cosh(2*u)*cos(2*v)
        sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -pi/2, pi/2), frame=False, color="red")
        Graphics3d Object

    Dini's spiral

    ::

        sage: u, v = var('u,v')
        sage: fx = cos(u)*sin(v)
        sage: fy = sin(u)*sin(v)
        sage: fz = (cos(v)+log(tan(v/2))) + 0.2*u
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 12.4), (v, 0.1, 2),frame=False, color="red")
        Graphics3d Object

    Catalan's surface (see
    http://xahlee.org/surface/catalan/catalan.html)::

        sage: u, v = var('u,v')
        sage: fx = u-sin(u)*cosh(v)
        sage: fy = 1-cos(u)*cosh(v)
        sage: fz = 4*sin(1/2*u)*sinh(v/2)
        sage: parametric_plot3d([fx, fy, fz], (u, -pi, 3*pi), (v, -2, 2), frame=False, color="red")
        Graphics3d Object

    A Conchoid::

        sage: u, v = var('u,v')
        sage: k = 1.2; k_2 = 1.2; a = 1.5
        sage: f = (k^u*(1+cos(v))*cos(u), k^u*(1+cos(v))*sin(u), k^u*sin(v)-a*k_2^u)
        sage: parametric_plot3d(f, (u,0,6*pi), (v,0,2*pi), plot_points=[40,40], texture=(0,0.5,0))
        Graphics3d Object

    A Möbius strip::

        sage: u,v = var("u,v")
        sage: parametric_plot3d([cos(u)*(1+v*cos(u/2)), sin(u)*(1+v*cos(u/2)), 0.2*v*sin(u/2)], (u,0, 4*pi+0.5), (v,0, 0.3),plot_points=[50,50])
        Graphics3d Object

    A Twisted Ribbon

    ::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([3*sin(u)*cos(v), 3*sin(u)*sin(v), cos(v)], (u,0, 2*pi), (v, 0, pi),plot_points=[50,50])
        Graphics3d Object

    An Ellipsoid::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([3*sin(u)*cos(v), 2*sin(u)*sin(v), cos(u)], (u,0, 2*pi), (v, 0, 2*pi),plot_points=[50,50], aspect_ratio=[1,1,1])
        Graphics3d Object

    A Cone::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([u*cos(v), u*sin(v), u], (u, -1, 1), (v, 0, 2*pi+0.5), plot_points=[50,50])
        Graphics3d Object

    A Paraboloid::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([u*cos(v), u*sin(v), u^2], (u, 0, 1), (v, 0, 2*pi+0.4), plot_points=[50,50])
        Graphics3d Object

    A Hyperboloid::

        sage: u, v = var('u,v')
        sage: plot3d(u^2-v^2, (u, -1, 1), (v, -1, 1), plot_points=[50,50])
        Graphics3d Object

    A weird looking surface - like a Möbius band but also an O::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([sin(u)*cos(u)*log(u^2)*sin(v), (u^2)^(1/6)*(cos(u)^2)^(1/4)*cos(v), sin(v)], (u, 0.001, 1), (v, -pi, pi+0.2), plot_points=[50,50])
        Graphics3d Object

    A heart, but not a cardioid (for my wife)::

        sage: u, v = var('u,v')
        sage: p1 = parametric_plot3d([sin(u)*cos(u)*log(u^2)*v*(1-v)/2, ((u^6)^(1/20)*(cos(u)^2)^(1/4)-1/2)*v*(1-v), v^(0.5)], (u, 0.001, 1), (v, 0, 1), plot_points=[70,70], color='red')
        sage: p2 = parametric_plot3d([-sin(u)*cos(u)*log(u^2)*v*(1-v)/2, ((u^6)^(1/20)*(cos(u)^2)^(1/4)-1/2)*v*(1-v), v^(0.5)], (u, 0.001, 1), (v, 0, 1), plot_points=[70,70], color='red')
        sage: show(p1+p2, frame=False)

    A Hyperhelicoidal::

        sage: u = var("u")
        sage: v = var("v")
        sage: fx = (sinh(v)*cos(3*u))/(1+cosh(u)*cosh(v))
        sage: fy = (sinh(v)*sin(3*u))/(1+cosh(u)*cosh(v))
        sage: fz = (cosh(v)*sinh(u))/(1+cosh(u)*cosh(v))
        sage: parametric_plot3d([fx, fy, fz], (u, -pi, pi), (v, -pi, pi), plot_points = [50,50], frame=False, color="red")
        Graphics3d Object

    A Helicoid (lines through a helix,
    http://en.wikipedia.org/wiki/Helix)::

        sage: u, v = var('u,v')
        sage: fx = sinh(v)*sin(u)
        sage: fy = -sinh(v)*cos(u)
        sage: fz = 3*u
        sage: parametric_plot3d([fx, fy, fz], (u, -pi, pi), (v, -pi, pi), plot_points = [50,50], frame=False, color="red")
        Graphics3d Object

    Kuen's surface
    (http://virtualmathmuseum.org/Surface/kuen/kuen.html)::

        sage: fx = (2*(cos(u) + u*sin(u))*sin(v))/(1+ u^2*sin(v)^2)
        sage: fy = (2*(sin(u) - u*cos(u))*sin(v))/(1+ u^2*sin(v)^2)
        sage: fz = log(tan(1/2 *v)) + (2*cos(v))/(1+ u^2*sin(v)^2)
        sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0.01, pi-0.01), plot_points = [50,50], frame=False, color="green")
        Graphics3d Object

    A 5-pointed star::

        sage: fx = cos(u)*cos(v)*(abs(cos(1*u/4))^0.5 + abs(sin(1*u/4))^0.5)^(-1/0.3)*(abs(cos(5*v/4))^1.7 + abs(sin(5*v/4))^1.7)^(-1/0.1)
        sage: fy = cos(u)*sin(v)*(abs(cos(1*u/4))^0.5 + abs(sin(1*u/4))^0.5)^(-1/0.3)*(abs(cos(5*v/4))^1.7 + abs(sin(5*v/4))^1.7)^(-1/0.1)
        sage: fz = sin(u)*(abs(cos(1*u/4))^0.5 + abs(sin(1*u/4))^0.5)^(-1/0.3)
        sage: parametric_plot3d([fx, fy, fz], (u, -pi/2, pi/2), (v, 0, 2*pi), plot_points = [50,50], frame=False, color="green")
        Graphics3d Object

    A cool self-intersecting surface (Eppener surface?)::

        sage: fx = u - u^3/3 + u*v^2
        sage: fy = v - v^3/3 + v*u^2
        sage: fz = u^2 - v^2
        sage: parametric_plot3d([fx, fy, fz], (u, -25, 25), (v, -25, 25), plot_points = [50,50], frame=False, color="green")
        Graphics3d Object

    The breather surface
    (http://en.wikipedia.org/wiki/Breather_surface)::

        sage: fx = (2*sqrt(0.84)*cosh(0.4*u)*(-(sqrt(0.84)*cos(v)*cos(sqrt(0.84)*v)) - sin(v)*sin(sqrt(0.84)*v)))/(0.4*((sqrt(0.84)*cosh(0.4*u))^2 + (0.4*sin(sqrt(0.84)*v))^2))
        sage: fy = (2*sqrt(0.84)*cosh(0.4*u)*(-(sqrt(0.84)*sin(v)*cos(sqrt(0.84)*v)) + cos(v)*sin(sqrt(0.84)*v)))/(0.4*((sqrt(0.84)*cosh(0.4*u))^2 + (0.4*sin(sqrt(0.84)*v))^2))
        sage: fz = -u + (2*0.84*cosh(0.4*u)*sinh(0.4*u))/(0.4*((sqrt(0.84)*cosh(0.4*u))^2 + (0.4*sin(sqrt(0.84)*v))^2))
        sage: parametric_plot3d([fx, fy, fz], (u, -13.2, 13.2), (v, -37.4, 37.4), plot_points = [90,90], frame=False, color="green")
        Graphics3d Object

    TESTS::

        sage: u, v = var('u,v')
        sage: plot3d(u^2-v^2, (u, -1, 1), (u, -1, 1))
        Traceback (most recent call last):
        ...
        ValueError: range variables should be distinct, but there are duplicates


    From :trac:`2858`::

        sage: parametric_plot3d((u,-u,v), (u,-10,10),(v,-10,10))
        Graphics3d Object
        sage: f(u)=u; g(v)=v^2; parametric_plot3d((g,f,f), (-10,10),(-10,10))
        Graphics3d Object

    From :trac:`5368`::

        sage: x, y = var('x,y')
        sage: plot3d(x*y^2 - sin(x), (x,-1,1), (y,-1,1))
        Graphics3d Object

    """
    # TODO:
    #   * Surface -- behavior of functions not defined everywhere -- see note above
    #   * Iterative refinement


    # color_function -- (default: "automatic") how to determine the color of curves and surfaces
    # color_function_scaling -- (default: True) whether to scale the input to color_function
    # exclusions -- (default: "automatic") u points or (u,v) conditions to exclude.
    #         (E.g., exclusions could be a function e = lambda u, v: False if u < v else True
    # exclusions_style -- (default: None) what to draw at excluded points
    # max_recursion -- (default: "automatic") maximum number of recursive subdivisions,
    #                   when ...
    # mesh -- (default: "automatic") how many mesh divisions in each direction to draw
    # mesh_functions -- (default: "automatic") how to determine the placement of mesh divisions
    # mesh_shading -- (default: None) how to shade regions between mesh divisions
    # plot_range -- (default: "automatic") range of values to include

    if is_Vector(f):
        f = tuple(f)

    if isinstance(f, (list,tuple)) and len(f) > 0 and isinstance(f[0], (list,tuple)):
        return sum([parametric_plot3d(v, urange, vrange, plot_points=plot_points, **kwds) for v in f])

    if not isinstance(f, (tuple, list)) or len(f) != 3:
        raise ValueError("f must be a list, tuple, or vector of length 3")

    if vrange is None:
        if plot_points == "automatic":
            plot_points = 75
        G = _parametric_plot3d_curve(f, urange, plot_points=plot_points, **kwds)
    else:
        if plot_points == "automatic":
            plot_points = [40,40]
        G = _parametric_plot3d_surface(f, urange, vrange, plot_points=plot_points, boundary_style=boundary_style, **kwds)
    G._set_extra_kwds(kwds)
    return G

def _parametric_plot3d_curve(f, urange, plot_points, **kwds):
    r"""
    Return a parametric three-dimensional space curve.
    This function is used internally by the
    :func:`parametric_plot3d` command.

    There are two ways this function is invoked by
    :func:`parametric_plot3d`.

    - ``parametric_plot3d([f_x, f_y, f_z], (u_min,
      u_max))``:
      `f_x, f_y, f_z` are three functions and
      `u_{\min}` and `u_{\max}` are real numbers

    - ``parametric_plot3d([f_x, f_y, f_z], (u, u_min,
      u_max))``:
      `f_x, f_y, f_z` can be viewed as functions of
      `u`

    INPUT:

    - ``f`` - a 3-tuple of functions or expressions, or vector of size 3

    - ``urange`` - a 2-tuple (u_min, u_max) or a 3-tuple
      (u, u_min, u_max)

    - ``plot_points`` - (default: "automatic", which is 75) initial
      number of sample points in each parameter; an integer.

    EXAMPLES:

    We demonstrate each of the two ways of calling this.  See
    :func:`parametric_plot3d` for many more examples.

    We do the first one with a lambda function, which creates a
    callable Python function that sends `u` to `u/10`::

        sage: parametric_plot3d( (sin, cos, lambda u: u/10), (0, 20)) # indirect doctest
        Graphics3d Object

    Now we do the same thing with symbolic expressions::

        sage: u = var('u')
        sage: parametric_plot3d( (sin(u), cos(u), u/10), (u, 0, 20))
        Graphics3d Object
    """
    from sage.plot.misc import setup_for_eval_on_grid
    g, ranges = setup_for_eval_on_grid(f, [urange], plot_points)
    f_x,f_y,f_z = g
    w = [(f_x(u), f_y(u), f_z(u)) for u in xsrange(*ranges[0], include_endpoint=True)]
    return line3d(w, **kwds)

def _parametric_plot3d_surface(f, urange, vrange, plot_points, boundary_style, **kwds):
    r"""
    Return a parametric three-dimensional space surface.
    This function is used internally by the
    :func:`parametric_plot3d` command.

    There are two ways this function is invoked by
    :func:`parametric_plot3d`.

    - ``parametric_plot3d([f_x, f_y, f_z], (u_min, u_max),
      (v_min, v_max))``:
      `f_x, f_y, f_z` are each functions of two variables

    - ``parametric_plot3d([f_x, f_y, f_z], (u, u_min,
      u_max), (v, v_min, v_max))``:
      `f_x, f_y, f_z` can be viewed as functions of
      `u` and `v`

    INPUT:

    - ``f`` - a 3-tuple of functions or expressions, or vector of size 3

    - ``urange`` - a 2-tuple (u_min, u_max) or a 3-tuple
      (u, u_min, u_max)

    - ``vrange`` - a 2-tuple (v_min, v_max) or a 3-tuple
      (v, v_min, v_max)

    - ``plot_points`` - (default: "automatic", which is [40,40]
      for surfaces) initial number of sample points in each parameter;
      a pair of integers.

    - ``boundary_style`` - (default: None, no boundary) a dict that describes
      how to draw the boundaries of regions by giving options that are passed
      to the line3d command.

    EXAMPLES:

    We demonstrate each of the two ways of calling this.  See
    :func:`parametric_plot3d` for many more examples.

    We do the first one with lambda functions::

        sage: f = (lambda u,v: cos(u), lambda u,v: sin(u)+cos(v), lambda u,v: sin(v))
        sage: parametric_plot3d(f, (0, 2*pi), (-pi, pi)) # indirect doctest
        Graphics3d Object

    Now we do the same thing with symbolic expressions::

        sage: u, v = var('u,v')
        sage: parametric_plot3d((cos(u), sin(u) + cos(v), sin(v)), (u, 0, 2*pi), (v, -pi, pi), mesh=True)
        Graphics3d Object
    """
    from sage.plot.misc import setup_for_eval_on_grid
    g, ranges = setup_for_eval_on_grid(f, [urange,vrange], plot_points)
    urange = srange(*ranges[0], include_endpoint=True)
    vrange = srange(*ranges[1], include_endpoint=True)
    G = ParametricSurface(g, (urange, vrange), **kwds)

    if boundary_style is not None:
        for u in (urange[0], urange[-1]):
            G += line3d([(g[0](u,v), g[1](u,v), g[2](u,v)) for v in vrange], **boundary_style)
        for v in (vrange[0], vrange[-1]):
            G += line3d([(g[0](u,v), g[1](u,v), g[2](u,v)) for u in urange], **boundary_style)
    return G
