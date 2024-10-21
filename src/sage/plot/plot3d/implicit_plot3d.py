# sage_setup: distribution = sagemath-plot
# sage.doctest: needs sage.symbolic
"""
Implicit plots
"""
from .implicit_surface import ImplicitSurface


def implicit_plot3d(f, xrange, yrange, zrange, **kwds):
    r"""
    Plot an isosurface of a function.

    INPUT:

    - ``f`` -- function

    - ``xrange`` -- a 2-tuple (x_min, x_max) or a 3-tuple (x, x_min, x_max)

    - ``yrange`` -- a 2-tuple (y_min, y_max) or a 3-tuple (y, y_min, y_max)

    - ``zrange`` -- a 2-tuple (z_min, z_max) or a 3-tuple (z, z_min, z_max)

    - ``plot_points`` -- (default: "automatic", which is 40) the number of
      function evaluations in each direction. (The number of cubes in the
      marching cubes algorithm will be one less than this). Can be a triple of
      integers, to specify a different resolution in each of x,y,z.

    - ``contour`` -- (default: 0) plot the isosurface ``f(x,y,z) == contour``.
      Can be a list, in which case multiple contours are plotted

    - ``region`` -- (default: ``None``) if region is given, it must be a Python
      callable. Only segments of the surface where ``region(x,y,z)`` returns a
      number `>0` will be included in the plot. (Note that returning a Python
      boolean is acceptable, since ``True == 1`` and ``False == 0``).

    EXAMPLES::

        sage: var('x,y,z')
        (x, y, z)

    A simple sphere::

        sage: implicit_plot3d(x^2+y^2+z^2==4, (x,-3,3), (y,-3,3), (z,-3,3))
        Graphics3d Object

    .. PLOT::

        var('x,y,z')
        F = x**2 + y**2 + z**2
        P = implicit_plot3d(F==4, (x,-3,3), (y,-3,3), (z,-3,3))
        sphinx_plot(P)

    A nested set of spheres with a hole cut out::

        sage: F = x^2 + y^2 + z^2
        sage: P = implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), plot_points=60,
        ....:                     contour=[1,3,5], color='aquamarine',
        ....:                     region=lambda x,y,z: x<=0.2 or y>=0.2 or z<=0.2)
        sage: P.show(viewer='tachyon')  # long time

    .. PLOT::

        var('x,y,z')
        F = x**2 + y**2 + z**2
        P = implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), plot_points=60, contour=[1,3,5],
                            region=lambda x,y,z: x<=0.2 or y>=0.2 or z<=0.2, color='aquamarine')
        sphinx_plot(P)

    A very pretty example, attributed to Douglas Summers-Stay (`archived page
    <http://web.archive.org/web/20080529033738/http://iat.ubalt.edu/summers/math/platsol.htm>`_)::

        sage: T = RDF(golden_ratio)
        sage: F = 2 - (cos(x+T*y) + cos(x-T*y) + cos(y+T*z)
        ....:           + cos(y-T*z) + cos(z-T*x) + cos(z+T*x))
        sage: r = 4.77
        sage: implicit_plot3d(F, (x,-r,r), (y,-r,r), (z,-r,r),
        ....:                 plot_points=40, color='darkkhaki').show(viewer='tachyon')

    .. PLOT::

        var('x,y,z')
        T = RDF(golden_ratio)
        F = 2 - (cos(x+T*y) + cos(x-T*y) + cos(y+T*z) + cos(y-T*z) + cos(z-T*x) + cos(z+T*x))
        r = 4.77
        V = implicit_plot3d(F, (x,-r,r), (y,-r,r), (z,-r,r), plot_points=40, color='darkkhaki')
        sphinx_plot(V)

    As I write this (but probably not as you read it), it's almost Valentine's
    day, so let's try a heart (from http://mathworld.wolfram.com/HeartSurface.html)

    ::

        sage: F = (x^2+9/4*y^2+z^2-1)^3 - x^2*z^3 - 9/(80)*y^2*z^3
        sage: r = 1.5
        sage: V = implicit_plot3d(F, (x,-r,r), (y,-r,r), (z,-r,r),
        ....:                     plot_points=80, color='red', smooth=False)
        sage: V.show(viewer='tachyon')  # long time

    .. PLOT::

        var('x,y,z')
        F = (x**2+9.0/4.0*y**2+z**2-1)**3 - x**2*z**3 - 9.0/(80)*y**2*z**3
        r = 1.5
        V = implicit_plot3d(F, (x,-r,r), (y,-r,r), (z,-r,r), plot_points=80, color='red', smooth=False)
        sphinx_plot(V)

    The same examples also work with the default Jmol viewer; for example::

        sage: T = RDF(golden_ratio)
        sage: F = 2 - (cos(x + T*y) + cos(x - T*y) + cos(y + T*z)
        ....:           + cos(y - T*z) + cos(z - T*x) + cos(z + T*x))
        sage: r = 4.77
        sage: implicit_plot3d(F, (x,-r,r), (y,-r,r), (z,-r,r),
        ....:                 plot_points=40, color='deepskyblue').show()

    Here we use ``smooth=True`` with a Tachyon graph::

        sage: implicit_plot3d(x^2 + y^2 + z^2, (x,-2,2), (y,-2,2), (z,-2,2),
        ....:                 contour=4, color='deepskyblue', smooth=True)
        Graphics3d Object

    .. PLOT::

        var('x,y,z')
        F = x**2 + y**2 + z**2
        P = implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), contour=4, color='deepskyblue', smooth=True)
        sphinx_plot(P)

    We explicitly specify a gradient function (in conjunction with ``smooth=True``)
    and invert the normals::

        sage: gx = lambda x, y, z: -(2*x + y^2 + z^2)
        sage: gy = lambda x, y, z: -(x^2 + 2*y + z^2)
        sage: gz = lambda x, y, z: -(x^2 + y^2 + 2*z)
        sage: implicit_plot3d(x^2+y^2+z^2, (x,-2,2), (y,-2,2), (z,-2,2),
        ....:                 contour=4, plot_points=40, smooth=True,
        ....:                 gradient=(gx, gy, gz)).show(viewer='tachyon')

    .. PLOT::

        var('x,y,z')
        gx = lambda x, y, z: -(2*x + y**2 + z**2)
        gy = lambda x, y, z: -(x**2 + 2*y + z**2)
        gz = lambda x, y, z: -(x**2 + y**2 + 2*z)
        P = implicit_plot3d(x**2+y**2+z**2, (x,-2,2), (y,-2,2), (z,-2,2), contour=4,
                           plot_points=40, smooth=True, gradient=(gx, gy, gz))
        sphinx_plot(P)

    A graph of two metaballs interacting with each other::

        sage: def metaball(x0, y0, z0):
        ....:     return 1 / ((x-x0)^2+(y-y0)^2+(z-z0)^2)
        sage: implicit_plot3d(metaball(-0.6,0,0) + metaball(0.6,0,0),
        ....:                 (x,-2,2), (y,-2,2), (z,-2,2),
        ....:                 plot_points=60, contour=2, color='seagreen')
        Graphics3d Object

    .. PLOT::

        var('x,y,z')
        def metaball(x0, y0, z0): return 1 / ((x-x0)**2+(y-y0)**2+(z-z0)**2)
        P = implicit_plot3d(metaball(-0.6,0,0) + metaball(0.6,0,0), (x,-2,2), (y,-2,2), (z,-2,2), plot_points=60, contour=2, color='seagreen')
        sphinx_plot(P)

    One can also color the surface using a coloring function and a
    colormap as follows. Note that the coloring function must take
    values in the interval `[0,1]`. ::

        sage: t = (sin(3*z)**2).function(x,y,z)
        sage: cm = colormaps.gist_rainbow
        sage: G = implicit_plot3d(x^2 + y^2 + z^2, (x,-2,2), (y,-2,2), (z,-2, 2),
        ....:                     contour=4, color=(t,cm), plot_points=100)
        sage: G.show(viewer='tachyon')  # long time

    .. PLOT::

        var('x,y,z')
        t = (sin(3*z)**2).function(x,y,z)
        cm = colormaps.gist_rainbow
        G = implicit_plot3d(x**2 + y**2 + z**2, (x,-2,2), (y,-2,2), (z,-2, 2),
                            contour=4, color=(t,cm), plot_points=60)
        sphinx_plot(G)

    Here is another colored example::

        sage: x, y, z = var('x,y,z')
        sage: t = (x).function(x,y,z)
        sage: cm = colormaps.PiYG
        sage: G = implicit_plot3d(x^4 + y^2 + z^2, (x,-2,2), (y,-2,2),(z,-2,2),
        ....:                     contour=4, color=(t,cm), plot_points=40)
        sage: G
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        t = (x).function(x,y,z)
        cm = colormaps.PiYG
        G = implicit_plot3d(x**4 + y**2 + z**2, (x,-2,2),
                           (y,-2,2),(z,-2,2), contour=4, color=(t,cm), plot_points=40)
        sphinx_plot(G)

    MANY MORE EXAMPLES:

    A kind of saddle::

        sage: implicit_plot3d(x^3 + y^2 - z^2, (x,-2,2), (y,-2,2), (z,-2,2),
        ....:                 plot_points=60, contour=0, color='lightcoral')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(x**3 + y**2 - z**2, (x,-2,2), (y,-2,2), (z,-2,2), plot_points=60, contour=0, color='lightcoral')
        sphinx_plot(G)

    A smooth surface with six radial openings::

        sage: implicit_plot3d(-(cos(x) + cos(y) + cos(z)),
        ....:                 (x,-4,4), (y,-4,4), (z,-4,4), color='orchid')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(-(cos(x) + cos(y) + cos(z)), (x,-4,4), (y,-4,4), (z,-4,4), color='orchid')
        sphinx_plot(G)

    A cube composed of eight conjoined blobs::

        sage: F = x^2 + y^2 + z^2 + cos(4*x) + cos(4*y) + cos(4*z) - 0.2
        sage: implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='mediumspringgreen')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = x**2 + y**2 + z**2 + cos(4*x) + cos(4*y) + cos(4*z) - 0.2
        G = implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='mediumspringgreen')
        sphinx_plot(G)

    A variation of the blob cube featuring heterogeneously sized blobs::

        sage: F = x^2 + y^2 + z^2 + sin(4*x) + sin(4*y) + sin(4*z) - 1
        sage: implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='lavenderblush')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = x**2 + y**2 + z**2 + sin(4*x) + sin(4*y) + sin(4*z) - 1
        G = implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='lavenderblush')
        sphinx_plot(G)

    A Klein bottle::

        sage: G = x^2 + y^2 + z^2
        sage: F = (G+2*y-1)*((G-2*y-1)^2-8*z^2) + 16*x*z*(G-2*y-1)
        sage: implicit_plot3d(F, (x,-3,3), (y,-3.1,3.1), (z,-4,4), color='moccasin')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = x**2 + y**2 + z**2
        F = (G+2*y-1)*((G-2*y-1)**2-8*z**2)+16*x*z*(G-2*y-1)
        G = implicit_plot3d(F, (x,-3,3), (y,-3.1,3.1), (z,-4,4), color='moccasin')
        sphinx_plot(G)

    A lemniscate::

        sage: F = 4*x^2*(x^2+y^2+z^2+z) + y^2*(y^2+z^2-1)
        sage: implicit_plot3d(F, (x,-0.5,0.5), (y,-1,1), (z,-1,1), color='deeppink')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = 4*x**2*(x**2+y**2+z**2+z) + y**2*(y**2+z**2-1)
        G = implicit_plot3d(F, (x,-0.5,0.5), (y,-1,1), (z,-1,1), color='deeppink')
        sphinx_plot(G)

    Drope::

        sage: implicit_plot3d(z - 4*x*exp(-x^2-y^2), (x,-2,2), (y,-2,2), (z,-1.7,1.7),
        ....:                 color='darkcyan')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(z - 4*x*exp(-x**2-y**2), (x,-2,2), (y,-2,2), (z,-1.7,1.7), color='darkcyan')
        sphinx_plot(G)

    A cube with a circular aperture on each face::

        sage: F = ((1/2.3)^2 * (x^2 + y^2 + z^2))^(-6) + ((1/2)^8 * (x^8 + y^8 + z^8))^6 - 1
        sage: implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='palevioletred')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = ((1/2.3)**2 * (x**2 + y**2 + z**2))**(-6) + ((1/2)**8 * (x**8 + y**8 + z**8))**6 - 1
        G = implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='palevioletred')
        sphinx_plot(G)

    A simple hyperbolic surface::

        sage: implicit_plot3d(x^2 + y - z^2, (x,-1,1), (y,-1,1), (z,-1,1),
        ....:                 color='darkslategray')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(x**2 + y - z**2, (x,-1,1), (y,-1,1), (z,-1,1), color='darkslategray')
        sphinx_plot(G)

    A hyperboloid::

        sage: implicit_plot3d(x^2 + y^2 - z^2 -0.3, (x,-2,2), (y,-2,2), (z,-1.8,1.8),
        ....:                 color='honeydew')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(x**2 + y**2 - z**2 -0.3, (x,-2,2), (y,-2,2), (z,-1.8,1.8), color='honeydew')
        sphinx_plot(G)

    Dupin cyclide (:wikipedia:`Dupin_cyclide`) ::

        sage: x, y, z, a, b, c, d = var('x,y,z,a,b,c,d')
        sage: a = 3.5
        sage: b = 3
        sage: c = sqrt(a^2 - b^2)
        sage: d = 2
        sage: F = (x^2 + y^2 + z^2 + b^2 - d^2)^2 - 4*(a*x-c*d)^2 - 4*b^2*y^2
        sage: implicit_plot3d(F, (x,-6,6), (y,-6,6), (z,-6,6), color='seashell')
        Graphics3d Object

    .. PLOT::

        x, y, z , a, b, c, d = var('x,y,z,a,b,c,d')
        a = 3.5
        b = 3
        c = sqrt(a**2 - b**2)
        d = 2
        F = (x**2 + y**2 + z**2 + b**2 - d**2)**2 - 4*(a*x-c*d)**2 - 4*b**2*y**2
        G = implicit_plot3d(F, (x,-6,6), (y,-6,6), (z,-6,6), color='seashell')
        sphinx_plot(G)

    Sinus::

        sage: implicit_plot3d(sin(pi*((x)^2+(y)^2))/2 + z, (x,-1,1), (y,-1,1), (z,-1,1),
        ....:                 color='rosybrown')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(sin(pi*((x)**2+(y)**2))/2 + z, (x,-1,1), (y,-1,1), (z,-1,1), color='rosybrown')
        sphinx_plot(G)

    A torus::

        sage: implicit_plot3d((sqrt(x*x+y*y)-3)^2 + z*z - 1, (x,-4,4), (y,-4,4), (z,-1,1),
        ....:                 color='indigo')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d((sqrt(x*x+y*y)-3)**2 + z*z - 1, (x,-4,4), (y,-4,4), (z,-1,1), color='indigo')
        sphinx_plot(G)

    An octahedron::

        sage: implicit_plot3d(abs(x) + abs(y) + abs(z) - 1, (x,-1,1), (y,-1,1), (z,-1,1),
        ....:                 color='olive')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(abs(x) + abs(y) + abs(z) - 1, (x,-1,1), (y,-1,1), (z,-1,1), color='olive')
        sphinx_plot(G)

    A cube::

        sage: implicit_plot3d(x^100 + y^100 + z^100 - 1, (x,-2,2), (y,-2,2), (z,-2,2),
        ....:                 color='lightseagreen')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(x**100 + y**100 + z**100 - 1, (x,-2,2), (y,-2,2), (z,-2,2), color='lightseagreen')
        sphinx_plot(G)

    Toupie::

        sage: implicit_plot3d((sqrt(x*x+y*y)-3)^3 + z*z - 1, (x,-4,4), (y,-4,4), (z,-6,6),
        ....:                 color='mintcream')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d((sqrt(x*x+y*y)-3)**3 + z*z - 1, (x,-4,4), (y,-4,4), (z,-6,6), color='mintcream')
        sphinx_plot(G)

    A cube with rounded edges::

        sage: F = x^4 + y^4 + z^4 - (x^2 + y^2 + z^2)
        sage: implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='mediumvioletred')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = x**4 + y**4 + z**4 - (x**2 + y**2 + z**2)
        G = implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='mediumvioletred')
        sphinx_plot(G)

    Chmutov::

        sage: F = x^4 + y^4 + z^4 - (x^2 + y^2 + z^2 - 0.3)
        sage: implicit_plot3d(F, (x,-1.5,1.5), (y,-1.5,1.5), (z,-1.5,1.5),
        ....:                 color='lightskyblue')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = x**4 + y**4 + z**4 - (x**2 + y**2 + z**2 - 0.3)
        G = implicit_plot3d(F, (x,-1.5,1.5), (y,-1.5,1.5), (z,-1.5,1.5), color='lightskyblue')
        sphinx_plot(G)

    Further Chmutov::

        sage: F = 2*(x^2*(3-4*x^2)^2+y^2*(3-4*y^2)^2+z^2*(3-4*z^2)^2) - 3
        sage: implicit_plot3d(F, (x,-1.3,1.3), (y,-1.3,1.3), (z,-1.3,1.3),
        ....:                 color='darksalmon')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = 2*(x**2*(3-4*x**2)**2+y**2*(3-4*y**2)**2+z**2*(3-4*z**2)**2) - 3
        G = implicit_plot3d(F, (x,-1.3,1.3), (y,-1.3,1.3), (z,-1.3,1.3), color='darksalmon')
        sphinx_plot(G)

    Clebsch surface::

        sage: F_1 = 81 * (x^3+y^3+z^3)
        sage: F_2 = 189 * (x^2*(y+z)+y^2*(x+z)+z^2*(x+y))
        sage: F_3 = 54 * x * y * z
        sage: F_4 = 126 * (x*y+x*z+y*z)
        sage: F_5 = 9 * (x^2+y^2+z^2)
        sage: F_6 = 9 * (x+y+z)
        sage: F = F_1 - F_2 + F_3 + F_4 - F_5 + F_6 + 1
        sage: implicit_plot3d(F, (x,-1,1), (y,-1,1), (z,-1,1), color='yellowgreen')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F_1 = 81 * (x**3+y**3+z**3)
        F_2 = 189 * (x**2*(y+z)+y**2*(x+z)+z**2*(x+y))
        F_3 = 54 * x * y * z
        F_4 = 126 * (x*y+x*z+y*z)
        F_5 = 9 * (x**2+y**2+z**2)
        F_6 = 9 * (x+y+z)
        F = F_1 - F_2 + F_3 + F_4 - F_5 + F_6 + 1
        G = implicit_plot3d(F, (x,-1,1), (y,-1,1), (z,-1,1), color='yellowgreen')
        sphinx_plot(G)

    Looks like a water droplet::

        sage: implicit_plot3d(x^2 + y^2 - (1-z)*z^2, (x,-1.5,1.5), (y,-1.5,1.5), (z,-1,1),
        ....:                 color='bisque')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(x**2 +y**2 -(1-z)*z**2, (x,-1.5,1.5), (y,-1.5,1.5), (z,-1,1), color='bisque')
        sphinx_plot(G)

    Sphere in a cage::

        sage: F = (x^8+z^30+y^8-(x^4 + z^50 + y^4 -0.3)) * (x^2+y^2+z^2-0.5)
        sage: implicit_plot3d(F, (x,-1.2,1.2), (y,-1.3,1.3), (z,-1.5,1.5),
        ....:                 color='firebrick')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = (x**8+z**30+y**8-(x**4 + z**50 + y**4 -0.3)) * (x**2+y**2+z**2-0.5)
        G = implicit_plot3d(F, (x,-1.2,1.2), (y,-1.3,1.3), (z,-1.5,1.5), color='firebrick')
        sphinx_plot(G)

    Ortho circle::

        sage: F = (((x^2+y^2-1)^2+z^2) * ((y^2+z^2-1)^2+x^2) * ((z^2+x^2-1)^2+y^2)
        ....:      - 0.075^2 * (1+3*(x^2+y^2+z^2)))
        sage: implicit_plot3d(F, (x,-1.5,1.5), (y,-1.5,1.5), (z,-1.5,1.5),
        ....:                 color='lemonchiffon')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = ((x**2+y**2-1)**2+z**2) * ((y**2+z**2-1)**2+x**2) * ((z**2+x**2-1)**2+y**2)-0.075**2 * (1+3*(x**2+y**2+z**2))
        G = implicit_plot3d(F, (x,-1.5,1.5), (y,-1.5,1.5), (z,-1.5,1.5), color='lemonchiffon')
        sphinx_plot(G)

    Cube sphere::

        sage: F = 12 - ((1/2.3)^2 *(x^2 + y^2 + z^2))^-6 - ((1/2)^8 * (x^8 + y^8 + z^8))^6
        sage: implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='rosybrown')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        F = 12 - ((1/2.3)**2 *(x**2 + y**2 + z**2))**-6 - ( (1/2)**8 * (x**8 + y**8 + z**8) )**6
        G = implicit_plot3d(F, (x,-2,2), (y,-2,2), (z,-2,2), color='rosybrown')
        sphinx_plot(G)

    Two cylinders intersect to make a cross::

        sage: implicit_plot3d((x^2+y^2-1) * (x^2+z^2-1) - 1, (x,-3,3), (y,-3,3), (z,-3,3),
        ....:                 color='burlywood')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d((x**2+y**2-1) * (x**2+z**2-1) - 1, (x,-3,3), (y,-3,3), (z,-3,3), color='burlywood')
        sphinx_plot(G)

    Three cylinders intersect in a similar fashion::

        sage: implicit_plot3d((x^2+y^2-1) * (x^2+z^2-1) * (y^2+z^2-1) - 1,
        ....:                 (x,-3,3), (y,-3,3), (z,-3,3), color='aqua')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d((x**2+y**2-1) * (x**2+z**2-1) * (y**2+z**2-1)-1, (x,-3,3), (y,-3,3), (z,-3,3), color='aqua')
        sphinx_plot(G)

    A sphere-ish object with twelve holes, four on each XYZ plane::

        sage: implicit_plot3d(3*(cos(x)+cos(y)+cos(z)) + 4*cos(x)*cos(y)*cos(z),
        ....:                 (x,-3,3), (y,-3,3), (z,-3,3), color='orangered')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(3*(cos(x)+cos(y)+cos(z)) + 4*cos(x)*cos(y)*cos(z), (x,-3,3), (y,-3,3), (z,-3,3), color='orangered')
        sphinx_plot(G)

    A gyroid::

        sage: implicit_plot3d(cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x),
        ....:                 (x,-4,4), (y,-4,4), (z,-4,4), color='sandybrown')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d(cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x), (x,-4,4), (y,-4,4), (z,-4,4), color='sandybrown')
        sphinx_plot(G)

    Tetrahedra::

        sage: implicit_plot3d((x^2+y^2+z^2)^2 + 8*x*y*z - 10*(x^2+y^2+z^2) + 25,
        ....:                 (x,-4,4), (y,-4,4), (z,-4,4), color='plum')
        Graphics3d Object

    .. PLOT::

        x, y, z = var('x,y,z')
        G = implicit_plot3d((x**2+y**2+z**2)**2 + 8*x*y*z - 10*(x**2+y**2+z**2) + 25, (x,-4,4), (y,-4,4), (z,-4,4), color='plum')
        sphinx_plot(G)

    TESTS:

    Test a separate resolution in the X direction; this should look like a
    regular sphere::

        sage: implicit_plot3d(x^2 + y^2 + z^2, (x,-2,2), (y,-2,2), (z,-2,2), plot_points=(10,40,40), contour=4)
        Graphics3d Object

    Test using different plot ranges in the different directions; each
    of these should generate half of a sphere.  Note that we need to use
    the ``aspect_ratio`` keyword to make it look right with the unequal
    plot ranges::

        sage: implicit_plot3d(x^2 + y^2 + z^2, (x,0,2), (y,-2,2), (z,-2,2), contour=4, aspect_ratio=1)
        Graphics3d Object

        sage: implicit_plot3d(x^2 + y^2 + z^2, (x,-2,2), (y,0,2), (z,-2,2), contour=4, aspect_ratio=1)
        Graphics3d Object

        sage: implicit_plot3d(x^2 + y^2 + z^2, (x,-2,2), (y,-2,2), (z,0,2), contour=4, aspect_ratio=1)
        Graphics3d Object

    Extra keyword arguments will be passed to show()::

        sage: implicit_plot3d(x^2 + y^2 + z^2, (x,-2,2), (y,-2,2), (z,-2,2), contour=4, viewer='tachyon')
        Graphics3d Object

    An implicit plot that does not include any surface in the view volume
    produces an empty plot::

        sage: implicit_plot3d(x^2 + y^2 + z^2 - 5000, (x,-2,2), (y,-2,2), (z,-2,2), plot_points=6)
        Graphics3d Object

    Make sure that implicit_plot3d does not error if the function cannot
    be symbolically differentiated::

        sage: implicit_plot3d(max_symbolic(x, y^2) - z, (x,-2,2), (y,-2,2), (z,-2,2), plot_points=6)
        Graphics3d Object

    TESTS:

    Check for :issue:`10599`::

        sage: var('x,y,z')
        (x, y, z)
        sage: M = matrix(3,[1,-1,-1,-1,3,1,-1,1,3])
        sage: v = 1/M.eigenvalues()[1]
        sage: implicit_plot3d(x^2+y^2+z^2==v, [x,-3,3], [y,-3,3],[z,-3,3])
        Graphics3d Object
    """
    # These options, related to rendering with smooth shading, are irrelevant
    # since IndexFaceSet does not support surface normals:
    # smooth: (default: ``False``) Whether to use vertex normals to produce a
    #   smooth-looking surface. False is slightly faster.
    # gradient: (default: None) If smooth is True (the default), then
    #   Tachyon rendering needs vertex normals. In that case, if gradient is None
    #   (the default), then we try to differentiate the function to get the
    #   gradient. If that fails, then we use central differencing on the scalar
    #   field. But it's also possible to specify the gradient; this must be either
    #   a single python callable that takes (x,y,z) and returns a tuple (dx,dy,dz)
    #   or a tuple of three callables that each take (x,y,z) and return dx, dy, dz
    #   respectively.

    G = ImplicitSurface(f, xrange, yrange, zrange, **kwds)
    G._set_extra_kwds(kwds)
    return G
