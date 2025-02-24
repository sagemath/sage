# sage.doctest: needs sage.symbolic
"""
Parametric plots
"""

from .parametric_surface import ParametricSurface
from .shapes2 import line3d
from sage.arith.srange import xsrange, srange
from sage.structure.element import Vector
from sage.misc.decorators import rename_keyword


@rename_keyword(alpha='opacity')
def parametric_plot3d(f, urange, vrange=None, plot_points='automatic',
                      boundary_style=None, **kwds):
    r"""
    Return a parametric three-dimensional space curve or surface.

    There are four ways to call this function:

    - ``parametric_plot3d([f_x, f_y, f_z], (u_min, u_max))``:
      `f_x, f_y, f_z` are three functions and
      `u_{\min}` and `u_{\max}` are real numbers

    - ``parametric_plot3d([f_x, f_y, f_z], (u, u_min, u_max))``:
      `f_x, f_y, f_z` can be viewed as functions of
      `u`

    - ``parametric_plot3d([f_x, f_y, f_z], (u_min, u_max),
      (v_min, v_max))``:
      `f_x, f_y, f_z` are each functions of two variables

    - ``parametric_plot3d([f_x, f_y, f_z], (u, u_min, u_max), (v, v_min, v_max))``:
      `f_x, f_y, f_z` can be viewed as functions of
      `u` and `v`


    INPUT:

    - ``f`` -- a 3-tuple of functions or expressions, or vector of size 3

    - ``urange`` -- a 2-tuple (u_min, u_max) or a 3-tuple
      (u, u_min, u_max)

    - ``vrange`` -- (optional, only used for surfaces) a
      2-tuple (v_min, v_max) or a 3-tuple (v, v_min, v_max)

    - ``plot_points`` -- (default: "automatic", which is
      75 for curves and [40,40] for surfaces) initial number of sample
      points in each parameter; an integer for a curve, and a pair of
      integers for a surface.

    - ``boundary_style`` -- (default: None, no boundary) a dict that describes
      how to draw the boundaries of regions by giving options that are passed
      to the line3d command.

    - ``mesh`` -- boolean (default: ``False``); whether to display
      mesh grid lines

    - ``dots`` -- boolean (default: ``False``); whether to display
      dots at mesh grid points

    .. NOTE::

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

           sage: parametric_plot3d((sin, cos, lambda u: u/10), (0,20))
           Graphics3d Object


        .. PLOT::

            sphinx_plot(parametric_plot3d((sin, cos, lambda u: u/10), (0,20)))

       Note above the lambda function, which creates a callable Python
       function that sends `u` to `u/10`.

    #. Next we draw the same plot as above, but using symbolic
       functions:

        ::

           sage: u = var('u')
           sage: parametric_plot3d((sin(u), cos(u), u/10), (u,0,20))
           Graphics3d Object


        .. PLOT::

            u = var('u')
            sphinx_plot(parametric_plot3d((sin(u), cos(u), u/10), (u,0,20)))

    #. We draw a parametric surface using 3 Python functions (defined
       using lambda):

        ::

           sage: f = (lambda u,v: cos(u), lambda u,v: sin(u)+cos(v), lambda u,v: sin(v))
           sage: parametric_plot3d(f, (0,2*pi), (-pi,pi))
           Graphics3d Object


        .. PLOT::

            f = (lambda u,v: cos(u), lambda u,v: sin(u)+cos(v), lambda u,v: sin(v))
            sphinx_plot(parametric_plot3d(f, (0,2*pi), (-pi,pi)))

    #. The same surface, but where the defining functions are
       symbolic:

        ::

           sage: u, v = var('u,v')
           sage: parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)), (u,0,2*pi), (v,-pi,pi))
           Graphics3d Object


        .. PLOT::

            u, v = var('u,v')
            sphinx_plot(parametric_plot3d((cos(u), sin(u)+cos(v) ,sin(v)), (u,0,2*pi), (v,-pi,pi)))

    The surface, but with a mesh::

           sage: u, v = var('u,v')
           sage: parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)), (u,0,2*pi), (v,-pi,pi), mesh=True)
           Graphics3d Object


    .. PLOT::

        u, v = var('u,v')
        sphinx_plot(parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)), (u,0,2*pi), (v,-pi,pi), mesh=True))

    We increase the number of plot points, and make the surface green
    and transparent::

        sage: parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)), (u,0,2*pi), (v,-pi,pi),
        ....: color='green', opacity=0.1, plot_points=[30,30])
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        sphinx_plot(parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)), (u,0,2*pi), (v,-pi,pi),
                    color='green', opacity=0.1, plot_points=[30,30]))

    One can also color the surface using a coloring function and a
    colormap as follows. Note that the coloring function must take
    values in the interval [0,1]. ::

        sage: u,v = var('u,v')
        sage: def cf(u, v): return sin(u+v/2)**2
        sage: P = parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)),
        ....:   (u,0,2*pi), (v,-pi,pi), color=(cf,colormaps.PiYG), plot_points=[60,60])
        sage: P.show(viewer='tachyon')

    .. PLOT::

        u,v = var('u,v')
        def cf(u, v): return sin(u+v/2)**2
        P = parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)),
            (u,0,2*pi), (v,-pi,pi), color=(cf,colormaps.PiYG), plot_points=[60,60])
        sphinx_plot(P)

    Another example, a colored Möbius band::

        sage: cm = colormaps.ocean
        sage: def c(x, y): return sin(x*y)**2
        sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
        sage: MoebiusStrip(5, 1, plot_points=200, color=(c,cm))
        Graphics3d Object

    .. PLOT::

        cm = colormaps.ocean
        def c(x, y): return sin(x*y)**2
        from sage.plot.plot3d.parametric_surface import MoebiusStrip
        sphinx_plot(MoebiusStrip(5, 1, plot_points=200, color=(c,cm)))

    Yet another colored example::

        sage: from sage.plot.plot3d.parametric_surface import ParametricSurface
        sage: cm = colormaps.autumn
        sage: def c(x, y): return sin(x*y)**2
        sage: def g(x, y): return x, y+sin(y), x**2 + y**2
        sage: ParametricSurface(g, (srange(-10,10,0.1), srange(-5,5.0,0.1)), color=(c,cm))
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.parametric_surface import ParametricSurface
        cm = colormaps.autumn
        def c(x, y): return sin(x*y)**2
        def g(x, y): return x, y+sin(y), x**2 + y**2
        sphinx_plot(ParametricSurface(g, (srange(-10,10,0.1), srange(-5,5.0,0.1)), color=(c,cm)))

    We call the space curve function but with polynomials instead of
    symbolic variables.

    ::

        sage: R.<t> = RDF[]
        sage: parametric_plot3d((t, t^2, t^3), (t,0,3))
        Graphics3d Object

    .. PLOT::

        t = var('t')
        R = RDF['t']
        sphinx_plot(parametric_plot3d((t, t**2, t**3), (t,0,3)))

    Next we plot the same curve, but because we use (0, 3) instead of
    (t, 0, 3), each polynomial is viewed as a callable function of one
    variable::

        sage: parametric_plot3d((t, t^2, t^3), (0,3))
        Graphics3d Object

    .. PLOT::

        t = var('t')
        R = RDF['t']
        sphinx_plot(parametric_plot3d((t, t**2, t**3), (0,3)))

    We do a plot but mix a symbolic input, and an integer::

        sage: t = var('t')
        sage: parametric_plot3d((1, sin(t), cos(t)), (t,0,3))
        Graphics3d Object

    .. PLOT::

        t = var('t')
        sphinx_plot(parametric_plot3d((1, sin(t), cos(t)), (t,0,3)))

    We specify a boundary style to show us the values of the function at its
    extrema::

        sage: u, v = var('u,v')
        sage: parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)), (u,0,pi), (v,0,pi),
        ....:              boundary_style={"color": "black", "thickness": 2})
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        P = parametric_plot3d((cos(u), sin(u)+cos(v), sin(v)), (u,0,pi), (v,0,pi),
                            boundary_style={"color":"black", "thickness":2})
        sphinx_plot(P)

    We can plot vectors::

        sage: x,y = var('x,y')
        sage: parametric_plot3d(vector([x-y, x*y, x*cos(y)]), (x,0,2), (y,0,2))
        Graphics3d Object

    .. PLOT::

        x,y = var('x,y')
        sphinx_plot(parametric_plot3d(vector([x-y, x*y, x*cos(y)]), (x,0,2), (y,0,2)))

    ::

        sage: t = var('t')
        sage: p = vector([1,2,3])
        sage: q = vector([2,-1,2])
        sage: parametric_plot3d(p*t+q, (t,0,2))
        Graphics3d Object

    .. PLOT::

        t = var('t')
        p = vector([1,2,3])
        q = vector([2,-1,2])
        sphinx_plot(parametric_plot3d(p*t+q, (t,0,2)))

    Any options you would normally use to specify the appearance of a curve are
    valid as entries in the ``boundary_style`` dict.

    MANY MORE EXAMPLES:

    We plot two interlinked tori::

        sage: u, v = var('u,v')
        sage: f1 = (4+(3+cos(v))*sin(u), 4+(3+cos(v))*cos(u), 4+sin(v))
        sage: f2 = (8+(3+cos(v))*cos(u), 3+sin(v), 4+(3+cos(v))*sin(u))
        sage: p1 = parametric_plot3d(f1, (u,0,2*pi), (v,0,2*pi), texture='red')
        sage: p2 = parametric_plot3d(f2, (u,0,2*pi), (v,0,2*pi), texture='blue')
        sage: p1 + p2
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f1 = (4+(3+cos(v))*sin(u), 4+(3+cos(v))*cos(u), 4+sin(v))
        f2 = (8+(3+cos(v))*cos(u), 3+sin(v), 4+(3+cos(v))*sin(u))
        p1 = parametric_plot3d(f1, (u,0,2*pi), (v,0,2*pi), texture='red')
        p2 = parametric_plot3d(f2, (u,0,2*pi), (v,0,2*pi), texture='blue')
        sphinx_plot(p1 + p2)

    A cylindrical Star of David::

        sage: u,v = var('u v')
        sage: K = (abs(cos(u))^200+abs(sin(u))^200)^(-1.0/200)
        sage: f_x = cos(u) * cos(v) * (abs(cos(3*v/4))^500+abs(sin(3*v/4))^500)^(-1/260) * K
        sage: f_y = cos(u) * sin(v) * (abs(cos(3*v/4))^500+abs(sin(3*v/4))^500)^(-1/260) * K
        sage: f_z = sin(u) * K
        sage: parametric_plot3d([f_x, f_y, f_z], (u, -pi, pi), (v, 0, 2*pi))
        Graphics3d Object

    .. PLOT::

        u,v = var('u v')
        K = (abs(cos(u))**200+abs(sin(u))**200)**(-1.0/200)
        f_x = cos(u) * cos(v) * (abs(cos(0.75*v))**500+abs(sin(0.75*v))**500)**(-1.0/260) * K
        f_y = cos(u)*sin(v)*(abs(cos(0.75*v))**500+abs(sin(0.75*v))**500)**(-1.0/260) * K
        f_z = sin(u) * K
        P = parametric_plot3d([f_x, f_y, f_z], (u, -pi, pi), (v, 0, 2*pi))
        sphinx_plot(P)

    Double heart::

        sage: u, v = var('u,v')
        sage: G1 = abs(sqrt(2)*tanh((u/sqrt(2))))
        sage: G2 = abs(sqrt(2)*tanh((v/sqrt(2))))
        sage: f_x = (abs(v) - abs(u) - G1 + G2)*sin(v)
        sage: f_y = (abs(v) - abs(u) - G1 - G2)*cos(v)
        sage: f_z = sin(u)*(abs(cos(u)) + abs(sin(u)))^(-1)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,pi), (v,-pi,pi))
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        G1 = abs(sqrt(2)*tanh((u/sqrt(2))))
        G2 = abs(sqrt(2)*tanh((v/sqrt(2))))
        f_x = (abs(v) - abs(u) - G1 + G2)*sin(v)
        f_y = (abs(v) - abs(u) - G1 - G2)*cos(v)
        f_z = sin(u)*(abs(cos(u)) + abs(sin(u)))**(-1)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,pi), (v,-pi,pi)))

    Heart::

        sage: u, v = var('u,v')
        sage: f_x = cos(u)*(4*sqrt(1-v^2)*sin(abs(u))^abs(u))
        sage: f_y = sin(u)*(4*sqrt(1-v^2)*sin(abs(u))^abs(u))
        sage: f_z = v
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-1,1), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = cos(u)*(4*sqrt(1-v**2)*sin(abs(u))**abs(u))
        f_y = sin(u) *(4*sqrt(1-v**2)*sin(abs(u))**abs(u))
        f_z = v
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-1,1), frame=False, color='red'))

    A Trefoil knot (:wikipedia:`Trefoil_knot`)::

        sage: u, v = var('u,v')
        sage: f_x = (4*(1+0.25*sin(3*v))+cos(u))*cos(2*v)
        sage: f_y = (4*(1+0.25*sin(3*v))+cos(u))*sin(2*v)
        sage: f_z = sin(u)+2*cos(3*v)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-pi,pi), frame=False, color='blue')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = (4*(1+0.25*sin(3*v))+cos(u))*cos(2*v)
        f_y = (4*(1+0.25*sin(3*v))+cos(u))*sin(2*v)
        f_z = sin(u)+2*cos(3*v)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-pi,pi), frame=False, color='blue'))

    Green bowtie::

        sage: u, v = var('u,v')
        sage: f_x = sin(u) / (sqrt(2) + sin(v))
        sage: f_y = sin(u) / (sqrt(2) + cos(v))
        sage: f_z = cos(u) / (1 + sqrt(2))
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-pi,pi), frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = sin(u) / (sqrt(2) + sin(v))
        f_y = sin(u) / (sqrt(2) + cos(v))
        f_z = cos(u) / (1 + sqrt(2))
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-pi,pi), frame=False, color='green'))

    Boy's surface (:wikipedia:`Boy%27s_surface` and https://mathcurve.com/surfaces/boy/boy.shtml)::

        sage: u, v = var('u,v')
        sage: K = cos(u) / (sqrt(2) - cos(2*u)*sin(3*v))
        sage: f_x = K * (cos(u)*cos(2*v)+sqrt(2)*sin(u)*cos(v))
        sage: f_y = K * (cos(u)*sin(2*v)-sqrt(2)*sin(u)*sin(v))
        sage: f_z = 3 * K * cos(u)
        sage: parametric_plot3d([f_x, f_y, f_z],      # long time
        ....:                   (u,-2*pi,2*pi),
        ....:                   (v,0,pi),
        ....:                   plot_points=[90,90],
        ....:                   frame=False,
        ....:                   color='orange')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        K = cos(u) / (sqrt(2) - cos(2*u)*sin(3*v))
        f_x = K * (cos(u)*cos(2*v)+sqrt(2)*sin(u)*cos(v))
        f_y = K * (cos(u)*sin(2*v)-sqrt(2)*sin(u)*sin(v))
        f_z = 3 * K * cos(u)
        P = parametric_plot3d([f_x, f_y, f_z], (u,-2*pi,2*pi), (v,0,pi),
                              plot_points=[90,90], frame=False, color='orange') # long time -- about 30 seconds
        sphinx_plot(P)

    Maeder's Owl also known as Bour's minimal surface (:wikipedia:`Bour%27s_minimal_surface`)::

        sage: u, v = var('u,v')
        sage: f_x = v*cos(u) - 0.5*v^2*cos(2*u)
        sage: f_y = -v*sin(u) - 0.5*v^2*sin(2*u)
        sage: f_z = 4 * v^1.5 * cos(3*u/2) / 3
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-2*pi,2*pi), (v,0,1),
        ....:                    plot_points=[90,90], frame=False, color='purple')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = v*cos(u) - 0.5*v**2*cos(2*u)
        f_y = -v*sin(u) - 0.5*v**2*sin(2*u)
        f_z = 4 * v**1.5 * cos(3*u/2) / 3
        P = parametric_plot3d([f_x, f_y, f_z], (u,-2*pi,2*pi), (v,0,1),
                              plot_points=[90,90], frame=False, color='purple')
        sphinx_plot(P)

    Bracelet::

        sage: u, v = var('u,v')
        sage: f_x = (2 + 0.2*sin(2*pi*u))*sin(pi*v)
        sage: f_y = 0.2 * cos(2*pi*u) * 3 * cos(2*pi*v)
        sage: f_z = (2 + 0.2*sin(2*pi*u))*cos(pi*v)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,pi/2), (v,0,3*pi/4), frame=False, color='gray')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = (2 + 0.2*sin(2*pi*u))*sin(pi*v)
        f_y = 0.2 * cos(2*pi*u) * 3* cos(2*pi*v)
        f_z = (2 + 0.2*sin(2*pi*u))*cos(pi*v)
        P = parametric_plot3d([f_x, f_y, f_z], (u,0,pi/2), (v,0,3*pi/4), frame=False, color='gray')
        sphinx_plot(P)


    Green goblet::

        sage: u, v = var('u,v')
        sage: f_x = cos(u) * cos(2*v)
        sage: f_y = sin(u) * cos(2*v)
        sage: f_z = sin(v)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,pi), frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = cos(u) * cos(2*v)
        f_y = sin(u) * cos(2*v)
        f_z = sin(v)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,pi), frame=False, color='green'))


    Funny folded surface - with square projection::

        sage: u, v = var('u,v')
        sage: f_x = cos(u) * sin(2*v)
        sage: f_y = sin(u) * cos(2*v)
        sage: f_z = sin(v)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = cos(u) * sin(2*v)
        f_y = sin(u) * cos(2*v)
        f_z = sin(v)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='green'))

    Surface of revolution of figure 8::

        sage: u, v = var('u,v')
        sage: f_x = cos(u) * sin(2*v)
        sage: f_y = sin(u) * sin(2*v)
        sage: f_z = sin(v)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = cos(u) * sin(2*v)
        f_y = sin(u) * sin(2*v)
        f_z = sin(v)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='green'))

    Yellow Whitney's umbrella (:wikipedia:`Whitney_umbrella`)::

        sage: u, v = var('u,v')
        sage: f_x = u*v
        sage: f_y = u
        sage: f_z = v^2
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-1,1), (v,-1,1), frame=False, color='yellow')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = u*v
        f_y = u
        f_z = v**2
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-1,1), (v,-1,1), frame=False, color='yellow'))

    Cross cap (:wikipedia:`Cross-cap`)::

        sage: u, v = var('u,v')
        sage: f_x = (1+cos(v)) * cos(u)
        sage: f_y = (1+cos(v)) * sin(u)
        sage: f_z = -tanh((2/3)*(u-pi)) * sin(v)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = (1+cos(v)) * cos(u)
        f_y = (1+cos(v)) * sin(u)
        f_z = -tanh((2.0/3.0)*(u-pi)) * sin(v)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='red'))

    Twisted torus::

        sage: u, v = var('u,v')
        sage: f_x = (3+sin(v)+cos(u)) * cos(2*v)
        sage: f_y = (3+sin(v)+cos(u)) * sin(2*v)
        sage: f_z = sin(u) + 2*cos(v)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = (3+sin(v)+cos(u)) * cos(2*v)
        f_y = (3+sin(v)+cos(u)) * sin(2*v)
        f_z = sin(u) + 2*cos(v)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='red'))

    Four intersecting discs::

        sage: u, v = var('u,v')
        sage: f_x = v*cos(u) - 0.5*v^2*cos(2*u)
        sage: f_y = -v*sin(u) - 0.5*v^2*sin(2*u)
        sage: f_z = 4 * v^1.5 * cos(3*u/2) / 3
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,4*pi), (v,0,2*pi), frame=False, color='red', opacity=0.7)
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = v*cos(u) - 0.5*v**2*cos(2*u)
        f_y = -v*sin(u) - 0.5*v**2*sin(2*u)
        f_z = 4 * v**1.5 * cos(3.0*u/2.0) /3
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,4*pi), (v,0,2*pi), frame=False, color='red', opacity=0.7))

    Steiner surface/Roman's surface (see
    :wikipedia:`Roman_surface` and
    :wikipedia:`Steiner_surface`)::

        sage: u, v = var('u,v')
        sage: f_x = (sin(2*u) * cos(v) * cos(v))
        sage: f_y = (sin(u) * sin(2*v))
        sage: f_z = (cos(u) * sin(2*v))
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-pi/2,pi/2), (v,-pi/2,pi/2), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = (sin(2*u) * cos(v) * cos(v))
        f_y = (sin(u) * sin(2*v))
        f_z = (cos(u) * sin(2*v))
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-pi/2,pi/2), (v,-pi/2,pi/2), frame=False, color='red'))

    Klein bottle? (see :wikipedia:`Klein_bottle`)::

        sage: u, v = var('u,v')
        sage: f_x = (3*(1+sin(v)) + 2*(1-cos(v)/2)*cos(u)) * cos(v)
        sage: f_y = (4+2*(1-cos(v)/2)*cos(u)) * sin(v)
        sage: f_z = -2 * (1-cos(v)/2) * sin(u)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = (3*(1+sin(v)) + 2*(1-cos(v)/2)*cos(u)) * cos(v)
        f_y = (4+2*(1-cos(v)/2)*cos(u)) * sin(v)
        f_z = -2 * (1-cos(v)/2) * sin(u)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='green'))

    A Figure 8 embedding of the Klein bottle (see
    :wikipedia:`Klein_bottle`)::

        sage: u, v = var('u,v')
        sage: f_x = (2+cos(v/2)*sin(u)-sin(v/2)*sin(2*u)) * cos(v)
        sage: f_y = (2+cos(v/2)*sin(u)-sin(v/2)*sin(2*u)) * sin(v)
        sage: f_z = sin(v/2)*sin(u) + cos(v/2)*sin(2*u)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = (2+cos(0.5*v)*sin(u)-sin(0.5*v)*sin(2*u)) * cos(v)
        f_y = (2+cos(0.5*v)*sin(u)-sin(0.5*v)*sin(2*u)) * sin(v)
        f_z = sin(v*0.5)*sin(u) + cos(v*0.5)*sin(2*u)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0,2*pi), frame=False, color='red'))

    Enneper's surface (see
    :wikipedia:`Enneper_surface`)::

        sage: u, v = var('u,v')
        sage: f_x = u - u^3/3 + u*v^2
        sage: f_y = v - v^3/3 + v*u^2
        sage: f_z = u^2 - v^2
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-2,2), (v,-2,2), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = u - u**3/3 + u*v**2
        f_y = v - v**3/3 + v*u**2
        f_z = u**2 - v**2
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-2,2), (v,-2,2), frame=False, color='red'))

    Henneberg's surface
    (see http://xahlee.org/surface/gallery_m.html)::

        sage: u, v = var('u,v')
        sage: f_x = 2*sinh(u)*cos(v) - (2/3)*sinh(3*u)*cos(3*v)
        sage: f_y = 2*sinh(u)*sin(v) + (2/3)*sinh(3*u)*sin(3*v)
        sage: f_z = 2 * cosh(2*u) * cos(2*v)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-1,1), (v,-pi/2,pi/2), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = 2.0*sinh(u)*cos(v) - (2.0/3.0)*sinh(3*u)*cos(3*v)
        f_y = 2.0*sinh(u)*sin(v) + (2.0/3.0)*sinh(3*u)*sin(3*v)
        f_z = 2.0 * cosh(2*u) * cos(2*v)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-1,1), (v,-pi/2,pi/2), frame=False, color='red'))

    Dini's spiral::

        sage: u, v = var('u,v')
        sage: f_x = cos(u) * sin(v)
        sage: f_y = sin(u) * sin(v)
        sage: f_z = (cos(v)+log(tan(v/2))) + 0.2*u
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,12.4), (v,0.1,2), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = cos(u) * sin(v)
        f_y = sin(u) * sin(v)
        f_z = (cos(v)+log(tan(v*0.5))) + 0.2*u
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,12.4), (v,0.1,2), frame=False, color='red'))

    Catalan's surface (see
    http://xahlee.org/surface/catalan/catalan.html)::

        sage: u, v = var('u,v')
        sage: f_x = u - sin(u)*cosh(v)
        sage: f_y = 1 - cos(u)*cosh(v)
        sage: f_z = 4 * sin(1/2*u) * sinh(v/2)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-pi,3*pi), (v,-2,2), frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = u - sin(u)*cosh(v)
        f_y = 1.0 - cos(u)*cosh(v)
        f_z = 4.0 * sin(0.5*u) * sinh(0.5*v)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-pi,3*pi), (v,-2,2), frame=False, color='red'))

    A Conchoid::

        sage: u, v = var('u,v')
        sage: k = 1.2; k_2 = 1.2; a = 1.5
        sage: f = (k^u*(1+cos(v))*cos(u), k^u*(1+cos(v))*sin(u), k^u*sin(v)-a*k_2^u)
        sage: parametric_plot3d(f, (u,0,6*pi), (v,0,2*pi), plot_points=[40,40], texture=(0,0.5,0))
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        k = 1.2; k_2 = 1.2; a = 1.5
        f = (k**u*(1+cos(v))*cos(u), k**u*(1+cos(v))*sin(u), k**u*sin(v)-a*k_2**u)
        sphinx_plot(parametric_plot3d(f, (u,0,6*pi), (v,0,2*pi), plot_points=[40,40], texture=(0,0.5,0)))

    A Möbius strip::

        sage: u,v = var("u,v")
        sage: parametric_plot3d([cos(u)*(1+v*cos(u/2)), sin(u)*(1+v*cos(u/2)), 0.2*v*sin(u/2)],
        ....:                   (u,0, 4*pi+0.5), (v,0, 0.3), plot_points=[50,50])
        Graphics3d Object

    .. PLOT::

        u,v = var("u,v")
        sphinx_plot(parametric_plot3d([cos(u)*(1+v*cos(u*0.5)), sin(u)*(1+v*cos(u*0.5)), 0.2*v*sin(u*0.5)],
                                      (u,0,4*pi+0.5), (v,0,0.3), plot_points=[50,50]))

    A Twisted Ribbon::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([3*sin(u)*cos(v), 3*sin(u)*sin(v), cos(v)],
        ....:                   (u,0,2*pi), (v,0,pi), plot_points=[50,50])
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        sphinx_plot(parametric_plot3d([3*sin(u)*cos(v), 3*sin(u)*sin(v), cos(v)],
                                      (u,0,2*pi), (v,0,pi), plot_points=[50,50]))

    An Ellipsoid::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([3*sin(u)*cos(v), 2*sin(u)*sin(v), cos(u)],
        ....:                   (u,0, 2*pi), (v, 0, 2*pi), plot_points=[50,50], aspect_ratio=[1,1,1])
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        sphinx_plot(parametric_plot3d([3*sin(u)*cos(v), 2*sin(u)*sin(v), cos(u)],
                                      (u,0,2*pi), (v,0,2*pi), plot_points=[50,50], aspect_ratio=[1,1,1]))

    A Cone::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([u*cos(v), u*sin(v), u], (u,-1,1), (v,0,2*pi+0.5), plot_points=[50,50])
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        sphinx_plot(parametric_plot3d([u*cos(v), u*sin(v), u], (u,-1,1), (v,0,2*pi+0.5), plot_points=[50,50]))

    A Paraboloid::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([u*cos(v), u*sin(v), u^2], (u,0,1), (v,0,2*pi+0.4), plot_points=[50,50])
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        sphinx_plot(parametric_plot3d([u*cos(v), u*sin(v), u**2], (u,0,1), (v,0,2*pi+0.4), plot_points=[50,50]))

    A Hyperboloid::

        sage: u, v = var('u,v')
        sage: plot3d(u^2-v^2, (u,-1,1), (v,-1,1), plot_points=[50,50])
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        sphinx_plot(plot3d(u**2-v**2, (u,-1,1), (v,-1,1), plot_points=[50,50]))

    A weird looking surface - like a Möbius band but also an O::

        sage: u, v = var('u,v')
        sage: parametric_plot3d([sin(u)*cos(u)*log(u^2)*sin(v), (u^2)^(1/6)*(cos(u)^2)^(1/4)*cos(v), sin(v)],
        ....:                   (u,0.001,1), (v,-pi,pi+0.2), plot_points=[50,50])
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        sphinx_plot(parametric_plot3d([sin(u)*cos(u)*log(u**2)*sin(v), (u**2)**(1.0/6.0)*(cos(u)**2)**(0.25)*cos(v), sin(v)],
                                      (u,0.001,1),
                                      (v,-pi,pi+0.2),
                                      plot_points=[50,50]))

    A heart, but not a cardioid (for my wife)::

        sage: u, v = var('u,v')
        sage: p1 = parametric_plot3d([sin(u)*cos(u)*log(u^2)*v*(1-v)/2, ((u^6)^(1/20)*(cos(u)^2)^(1/4)-1/2)*v*(1-v), v^(0.5)],
        ....:                        (u,0.001,1), (v,0,1), plot_points=[70,70], color='red')
        sage: p2 = parametric_plot3d([-sin(u)*cos(u)*log(u^2)*v*(1-v)/2, ((u^6)^(1/20)*(cos(u)^2)^(1/4)-1/2)*v*(1-v), v^(0.5)],
        ....:                        (u, 0.001,1), (v,0,1), plot_points=[70,70], color='red')
        sage: show(p1+p2)

    .. PLOT::

        u, v = var('u,v')
        p1 = parametric_plot3d([sin(u)*cos(u)*log(u**2)*v*(1-v)*0.5, ((u**6)**(1/20.0)*(cos(u)**2)**(0.25)-0.5)*v*(1-v), v**(0.5)],
                               (u,0.001,1), (v,0,1), plot_points=[70,70], color='red')
        p2 = parametric_plot3d([-sin(u)*cos(u)*log(u**2)*v*(1-v)*0.5, ((u**6)**(1/20.0)*(cos(u)**2)**(0.25)-0.5)*v*(1-v), v**(0.5)],
                               (u,0.001,1), (v,0,1), plot_points=[70,70], color='red')
        sphinx_plot(p1+p2)

    A Hyperhelicoidal::

        sage: u = var("u")
        sage: v = var("v")
        sage: f_x = (sinh(v)*cos(3*u)) / (1+cosh(u)*cosh(v))
        sage: f_y = (sinh(v)*sin(3*u)) / (1+cosh(u)*cosh(v))
        sage: f_z = (cosh(v)*sinh(u)) / (1+cosh(u)*cosh(v))
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-pi,pi), plot_points=[50,50], frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u = var("u")
        v = var("v")
        f_x = (sinh(v)*cos(3*u)) / (1+cosh(u)*cosh(v))
        f_y = (sinh(v)*sin(3*u)) / (1+cosh(u)*cosh(v))
        f_z = (cosh(v)*sinh(u)) / (1+cosh(u)*cosh(v))
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-pi,pi), plot_points=[50,50], frame=False, color='red'))

    A Helicoid (lines through a helix,
    :wikipedia:`Helix`)::

        sage: u, v = var('u,v')
        sage: f_x = sinh(v) * sin(u)
        sage: f_y = -sinh(v) * cos(u)
        sage: f_z = 3 * u
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-pi,pi), plot_points=[50,50], frame=False, color='red')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = sinh(v) * sin(u)
        f_y = -sinh(v) * cos(u)
        f_z = 3 * u
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-pi,pi), (v,-pi,pi), plot_points=[50,50], frame=False, color='red'))

    Kuen's surface
    (http://virtualmathmuseum.org/Surface/kuen/kuen.html)::

        sage: f_x = (2*(cos(u) + u*sin(u))*sin(v))/(1+ u^2*sin(v)^2)
        sage: f_y = (2*(sin(u) - u*cos(u))*sin(v))/(1+ u^2*sin(v)^2)
        sage: f_z = log(tan(1/2 *v)) + (2*cos(v))/(1+ u^2*sin(v)^2)
        sage: parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0.01,pi-0.01), plot_points=[50,50], frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = (2.0*(cos(u)+u*sin(u))*sin(v)) / (1.0+u**2*sin(v)**2)
        f_y = (2.0*(sin(u)-u*cos(u))*sin(v)) / (1.0+u**2*sin(v)**2)
        f_z = log(tan(0.5 *v)) + (2*cos(v))/(1.0+u**2*sin(v)**2)
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,0,2*pi), (v,0.01,pi-0.01), plot_points=[50,50], frame=False, color='green'))

    A 5-pointed star::

        sage: G1 = (abs(cos(u/4))^0.5+abs(sin(u/4))^0.5)^(-1/0.3)
        sage: G2 = (abs(cos(5*v/4))^1.7+abs(sin(5*v/4))^1.7)^(-1/0.1)
        sage: f_x = cos(u) * cos(v) * G1 * G2
        sage: f_y = cos(u) * sin(v) * G1 * G2
        sage: f_z = sin(u) * G1
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-pi/2,pi/2), (v,0,2*pi), plot_points=[50,50], frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        G1 = (abs(cos(u/4))**0.5+abs(sin(u/4))**0.5)**(-1/0.3)
        G2 = (abs(cos(5*v/4))**1.7+abs(sin(5*v/4))**1.7)**(-1/0.1)
        f_x = cos(u) * cos(v) * G1 * G2
        f_y = cos(u) * sin(v) * G1 * G2
        f_z = sin(u) * G1
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-pi/2,pi/2), (v,0,2*pi), plot_points=[50,50], frame=False, color='green'))

    A cool self-intersecting surface (Eppener surface?)::

        sage: f_x = u - u^3/3 + u*v^2
        sage: f_y = v - v^3/3 + v*u^2
        sage: f_z = u^2 - v^2
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-25,25), (v,-25,25), plot_points=[50,50], frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        f_x = u - u**3/3 + u*v**2
        f_y = v - v**3/3 + v*u**2
        f_z = u**2 - v**2
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-25,25), (v,-25,25), plot_points=[50,50], frame=False, color='green'))

    The breather surface
    (:wikipedia:`Breather_surface`)::

        sage: K = sqrt(0.84)
        sage: G = (0.4*((K*cosh(0.4*u))^2 + (0.4*sin(K*v))^2))
        sage: f_x = (2*K*cosh(0.4*u)*(-(K*cos(v)*cos(K*v)) - sin(v)*sin(K*v)))/G
        sage: f_y = (2*K*cosh(0.4*u)*(-(K*sin(v)*cos(K*v)) + cos(v)*sin(K*v)))/G
        sage: f_z = -u + (2*0.84*cosh(0.4*u)*sinh(0.4*u))/G
        sage: parametric_plot3d([f_x, f_y, f_z], (u,-13.2,13.2), (v,-37.4,37.4), plot_points=[90,90], frame=False, color='green')
        Graphics3d Object

    .. PLOT::

        u, v = var('u,v')
        K = sqrt(0.84)
        G = (0.4*((K*cosh(0.4*u))**2 + (0.4*sin(K*v))**2))
        f_x = (2*K*cosh(0.4*u)*(-(K*cos(v)*cos(K*v)) - sin(v)*sin(K*v)))/G
        f_y = (2*K*cosh(0.4*u)*(-(K*sin(v)*cos(K*v)) + cos(v)*sin(K*v)))/G
        f_z = -u + (2*0.84*cosh(0.4*u)*sinh(0.4*u))/G
        sphinx_plot(parametric_plot3d([f_x, f_y, f_z], (u,-13.2,13.2), (v,-37.4,37.4), plot_points=[90,90], frame=False, color='green'))

    TESTS::

        sage: u, v = var('u,v')
        sage: plot3d(u^2-v^2, (u,-1,1), (u,-1,1))
        Traceback (most recent call last):
        ...
        ValueError: range variables should be distinct, but there are duplicates


    From :issue:`2858`::

        sage: parametric_plot3d((u,-u,v), (u,-10,10),(v,-10,10))
        Graphics3d Object
        sage: f(u)=u; g(v)=v^2; parametric_plot3d((g,f,f), (-10,10),(-10,10))
        Graphics3d Object

    From :issue:`5368`::

        sage: x, y = var('x,y')
        sage: plot3d(x*y^2 - sin(x), (x,-1,1), (y,-1,1))
        Graphics3d Object
    """
    # TODO:
    #   * Surface -- behavior of functions not defined everywhere -- see note above
    #   * Iterative refinement

    # color_function -- (default: ``'automatic'``) how to determine the color of curves and surfaces
    # color_function_scaling -- boolean (default: ``True``); whether to scale the input to color_function
    # exclusions -- (default: ``'automatic'``) u points or (u,v) conditions to exclude.
    #         (E.g., exclusions could be a function e = lambda u, v: False if u < v else True
    # exclusions_style -- (default: ``None``) what to draw at excluded points
    # max_recursion -- (default: ``'automatic'``) maximum number of recursive subdivisions,
    #                   when ...
    # mesh -- (default: ``'automatic'``) how many mesh divisions in each direction to draw
    # mesh_functions -- (default: ``'automatic'``) how to determine the placement of mesh divisions
    # mesh_shading -- (default: ``None``) how to shade regions between mesh divisions
    # plot_range -- (default: ``'automatic'``) range of values to include

    if isinstance(f, Vector):
        f = tuple(f)

    if isinstance(f, (list, tuple)) and len(f) > 0 and isinstance(f[0], (list, tuple)):
        return sum([parametric_plot3d(v, urange, vrange, plot_points=plot_points, **kwds) for v in f])

    if not isinstance(f, (tuple, list)) or len(f) != 3:
        raise ValueError("f must be a list, tuple, or vector of length 3")

    if vrange is None:
        if plot_points == "automatic":
            plot_points = 75
        G = _parametric_plot3d_curve(f, urange, plot_points=plot_points, **kwds)
    else:
        if plot_points == "automatic":
            plot_points = [40, 40]
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

    - ``f`` -- a 3-tuple of functions or expressions, or vector of size 3

    - ``urange`` -- a 2-tuple (u_min, u_max) or a 3-tuple
      (u, u_min, u_max)

    - ``plot_points`` -- integer (default: 75, "automatic"); initial
      number of sample points in each parameter

    EXAMPLES:

    We demonstrate each of the two ways of calling this.  See
    :func:`parametric_plot3d` for many more examples.

    We do the first one with a lambda function, which creates a
    callable Python function that sends `u` to `u/10`::

        sage: parametric_plot3d((sin, cos, lambda u: u/10), (0,20)) # indirect doctest
        Graphics3d Object

    Now we do the same thing with symbolic expressions::

        sage: u = var('u')
        sage: parametric_plot3d((sin(u), cos(u), u/10), (u,0,20))
        Graphics3d Object
    """
    from sage.plot.misc import setup_for_eval_on_grid
    g, ranges = setup_for_eval_on_grid(f, [urange], plot_points)
    f_x, f_y, f_z = g
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

    - ``f`` -- a 3-tuple of functions or expressions, or vector of size 3

    - ``urange`` -- a 2-tuple (u_min, u_max) or a 3-tuple
      (u, u_min, u_max)

    - ``vrange`` -- a 2-tuple (v_min, v_max) or a 3-tuple
      (v, v_min, v_max)

    - ``plot_points`` -- (default: "automatic", which is [40,40]
      for surfaces) initial number of sample points in each parameter;
      a pair of integers.

    - ``boundary_style`` -- (default: None, no boundary) a dict that describes
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
    g, ranges = setup_for_eval_on_grid(f, [urange, vrange], plot_points)
    urange = srange(*ranges[0], include_endpoint=True)
    vrange = srange(*ranges[1], include_endpoint=True)
    G = ParametricSurface(g, (urange, vrange), **kwds)

    if boundary_style is not None:
        for u in (urange[0], urange[-1]):
            G += line3d([(g[0](u,v), g[1](u,v), g[2](u,v)) for v in vrange], **boundary_style)
        for v in (vrange[0], vrange[-1]):
            G += line3d([(g[0](u,v), g[1](u,v), g[2](u,v)) for u in urange], **boundary_style)
    return G
