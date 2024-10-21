# sage_setup: distribution = sagemath-symbolics
r"""
Common parametrized surfaces in 3D.

AUTHORS::

- Joris Vankerschaver (2012-06-16)
"""
#*****************************************************************************
#       Copyright (C) 2010  Joris Vankerschaver <joris.vankerschaver@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.symbolic.constants import pi
from sage.functions.log import log
from sage.functions.trig import sin, cos, tan
from sage.functions.hyperbolic import cosh, tanh
from sage.symbolic.ring import var
from sage.geometry.riemannian_manifolds.parametrized_surface3d import \
    ParametrizedSurface3D


class SurfaceGenerators:
    """
    A class consisting of generators for several common parametrized surfaces
    in 3D.
    """
    @staticmethod
    def Catenoid(c=1, name='Catenoid'):
        r"""
        Return a catenoid surface, with parametric representation.

        .. MATH::

            \begin{aligned}
              x(u, v) & = c \cosh(v/c) \cos(u); \\
              y(u, v) & = c \cosh(v/c) \sin(u); \\
              z(u, v) & = v.
            \end{aligned}

        INPUT:

        - ``c`` -- surface parameter

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Catenoid`.

        EXAMPLES::

            sage: cat = surfaces.Catenoid(); cat
            Parametrized surface ('Catenoid') with equation (cos(u)*cosh(v), cosh(v)*sin(u), v)
            sage: cat.plot()                                                            # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        catenoid_eq = [c*cosh(v/c)*cos(u), c*cosh(v/c)*sin(u), v]
        coords = ((u, 0, 2 * pi), (v, -1, 1))
        return ParametrizedSurface3D(catenoid_eq, coords, name)

    @staticmethod
    def Crosscap(r=1, name='Crosscap'):
        r"""
        Return a crosscap surface, with parametrization.

        .. MATH::

            \begin{aligned}
              x(u, v) & = r(1 + \cos(v)) \cos(u); \\
              y(u, v) & = r(1 + \cos(v)) \sin(u); \\
              z(u, v) & = - r\tanh(u - \pi) \sin(v).
            \end{aligned}

        INPUT:

        - ``r`` -- surface parameter

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Cross-cap`.

        EXAMPLES::

            sage: crosscap = surfaces.Crosscap(); crosscap
            Parametrized surface ('Crosscap') with equation ((cos(v) + 1)*cos(u), (cos(v) + 1)*sin(u), -sin(v)*tanh(-pi + u))
            sage: crosscap.plot()                                                       # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        crosscap_eq = [r*(1+cos(v))*cos(u), r*(1+cos(v))*sin(u),
                       -tanh(u-pi)*r*sin(v)]
        coords = ((u, 0, 2 * pi), (v, 0, 2 * pi))
        return ParametrizedSurface3D(crosscap_eq, coords, name)

    @staticmethod
    def Dini(a=1, b=1, name="Dini's surface"):
        r"""
        Return Dini's surface, with parametrization.

        .. MATH::

            \begin{aligned}
              x(u, v) & = a \cos(u)\sin(v); \\
              y(u, v) & = a \sin(u)\sin(v); \\
              z(u, v) & = u + \log(\tan(v/2)) + \cos(v).
            \end{aligned}

        INPUT:

        - ``a``, ``b`` -- surface parameters

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Dini%27s_surface`.

        EXAMPLES::

            sage: dini = surfaces.Dini(a=3, b=4); dini
            Parametrized surface ('Dini's surface') with equation (3*cos(u)*sin(v), 3*sin(u)*sin(v), 4*u + 3*cos(v) + 3*log(tan(1/2*v)))
            sage: dini.plot()                                                           # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        dini_eq = [a*cos(u)*sin(v), a*sin(u)*sin(v),
                   a*(cos(v) + log(tan(v/2))) + b*u]
        coords = ((u, 0, 2 * pi), (v, 0, 2 * pi))
        return ParametrizedSurface3D(dini_eq, coords, name)

    @staticmethod
    def Ellipsoid(center=(0, 0, 0), axes=(1, 1, 1), name='Ellipsoid'):
        r"""
        Return an ellipsoid centered at ``center`` whose semi-principal axes
        have lengths given by the components of ``axes``. The
        parametrization of the ellipsoid is given by

        .. MATH::

            \begin{aligned}
              x(u, v) & = x_0 + a \cos(u) \cos(v); \\
              y(u, v) & = y_0 + b \sin(u) \cos(v); \\
              z(u, v) & = z_0 + c \sin(v).
            \end{aligned}

        INPUT:

        - ``center`` -- 3-tuple; coordinates of the center of the ellipsoid

        - ``axes`` -- 3-tuple; lengths of the semi-principal axes

        - ``name`` -- string; name of the ellipsoid

        For more information, see :wikipedia:`Ellipsoid`.

        EXAMPLES::

            sage: ell = surfaces.Ellipsoid(axes=(1, 2, 3)); ell
            Parametrized surface ('Ellipsoid') with equation (cos(u)*cos(v), 2*cos(v)*sin(u), 3*sin(v))
            sage: ell.plot()                                                            # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        x, y, z = center
        a, b, c = axes
        ellipsoid_parametric_eq = [x + a*cos(u)*cos(v),
                                   y + b*sin(u)*cos(v),
                                   z + c*sin(v)]
        coords = ((u, 0, 2 * pi), (v, -pi / 2, pi / 2))
        return ParametrizedSurface3D(ellipsoid_parametric_eq, coords, name)

    @staticmethod
    def Enneper(name="Enneper's surface"):
        r"""
        Return Enneper's surface, with parametrization.

        .. MATH::

            \begin{aligned}
              x(u, v) & = u(1 - u^2/3 + v^2)/3; \\
              y(u, v) & = -v(1 - v^2/3 + u^2)/3; \\
              z(u, v) & = (u^2 - v^2)/3.
            \end{aligned}

        INPUT:

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Enneper_surface`.

        EXAMPLES::

            sage: enn = surfaces.Enneper(); enn
            Parametrized surface ('Enneper's surface') with equation (-1/9*(u^2 - 3*v^2 - 3)*u, -1/9*(3*u^2 - v^2 + 3)*v, 1/3*u^2 - 1/3*v^2)
            sage: enn.plot()                                                            # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        enneper_eq = [u*(1-u**2/3+v**2)/3, -v*(1-v**2/3+u**2)/3, (u**2-v**2)/3]
        coords = ((u, -3, 3), (v, -3, 3))
        return ParametrizedSurface3D(enneper_eq, coords, name)

    @staticmethod
    def Helicoid(h=1, name='Helicoid'):
        r"""
        Return a helicoid surface, with parametrization.

        .. MATH::

            \begin{aligned}
              x(\rho, \theta) & = \rho \cos(\theta); \\
              y(\rho, \theta) & = \rho \sin(\theta); \\
              z(\rho, \theta) & = h\theta/(2\pi).
            \end{aligned}

        INPUT:

        - ``h`` -- distance along the z-axis between two
          successive turns of the helicoid

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Helicoid`.

        EXAMPLES::

            sage: helicoid = surfaces.Helicoid(h=2); helicoid
            Parametrized surface ('Helicoid') with equation (rho*cos(theta), rho*sin(theta), theta/pi)
            sage: helicoid.plot()                                                       # needs sage.plot
            Graphics3d Object
        """
        rho, theta = var('rho, theta')
        helicoid_eq = [rho*cos(theta), rho*sin(theta), h*theta/(2*pi)]
        coords = ((rho, -2, 2), (theta, 0, 2 * pi))
        return ParametrizedSurface3D(helicoid_eq, coords, name)

    @staticmethod
    def Klein(r=1, name="Klein bottle"):
        r"""
        Return the Klein bottle, in the figure-8 parametrization given by

        .. MATH::

            \begin{aligned}
              x(u, v) & = (r + \cos(u/2)\cos(v) - \sin(u/2)\sin(2v)) \cos(u); \\
              y(u, v) & = (r + \cos(u/2)\cos(v) - \sin(u/2)\sin(2v)) \sin(u); \\
              z(u, v) & = \sin(u/2)\cos(v) + \cos(u/2)\sin(2v).
            \end{aligned}

        INPUT:

        - ``r`` -- radius of the "figure-8" circle

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Klein_bottle`.

        EXAMPLES::

            sage: klein = surfaces.Klein(); klein
            Parametrized surface ('Klein bottle') with equation (-(sin(1/2*u)*sin(2*v) - cos(1/2*u)*sin(v) - 1)*cos(u), -(sin(1/2*u)*sin(2*v) - cos(1/2*u)*sin(v) - 1)*sin(u), cos(1/2*u)*sin(2*v) + sin(1/2*u)*sin(v))
            sage: klein.plot()                                                          # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        x = (r + cos(u/2)*sin(v) - sin(u/2)*sin(2*v))*cos(u)
        y = (r + cos(u/2)*sin(v) - sin(u/2)*sin(2*v))*sin(u)
        z = sin(u/2)*sin(v) + cos(u/2)*sin(2*v)
        klein_eq = [x, y, z]
        coords = ((u, 0, 2 * pi), (v, 0, 2 * pi))

        return ParametrizedSurface3D(klein_eq, coords, name)

    @staticmethod
    def MonkeySaddle(name="Monkey saddle"):
        r"""
        Return a monkey saddle surface, with equation

        .. MATH::

            z = x^3 - 3xy^2.

        INPUT:

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Monkey_saddle`.

        EXAMPLES::

            sage: saddle = surfaces.MonkeySaddle(); saddle
            Parametrized surface ('Monkey saddle') with equation (u, v, u^3 - 3*u*v^2)
            sage: saddle.plot()                                                         # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        monkey_eq = [u, v, u**3 - 3*u*v**2]
        coords = ((u, -2, 2), (v, -2, 2))

        return ParametrizedSurface3D(monkey_eq, coords, name)

    @staticmethod
    def Paraboloid(a=1, b=1, c=1, elliptic=True, name=None):
        r"""
        Return a paraboloid with equation.

        .. MATH::

            \frac{z}{c} = \pm \frac{x^2}{a^2} + \frac{y^2}{b^2}

        When the plus sign is selected, the paraboloid is elliptic. Otherwise
        the surface is a hyperbolic paraboloid.

        INPUT:

        - ``a``, ``b``, ``c`` -- surface parameters

        - ``elliptic`` -- boolean (default: ``True``); whether to create an
          elliptic or hyperbolic paraboloid

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Paraboloid`.

        EXAMPLES::

            sage: epar = surfaces.Paraboloid(1, 3, 2); epar
            Parametrized surface ('Elliptic paraboloid') with equation (u, v, 2*u^2 + 2/9*v^2)
            sage: epar.plot()                                                           # needs sage.plot
            Graphics3d Object

            sage: hpar = surfaces.Paraboloid(2, 3, 1, elliptic=False); hpar
            Parametrized surface ('Hyperbolic paraboloid') with equation (u, v, -1/4*u^2 + 1/9*v^2)
            sage: hpar.plot()                                                           # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        x = u
        y = v
        if elliptic:
            z = c*(v**2/b**2 + u**2/a**2)
        else:
            z = c*(v**2/b**2 - u**2/a**2)
        paraboloid_eq = [x, y, z]
        coords = ((u, -3, 3), (v, -3, 3))

        if name is None:
            if elliptic:
                name = "Elliptic paraboloid"
            else:
                name = "Hyperbolic paraboloid"

        return ParametrizedSurface3D(paraboloid_eq, coords, name)

    @staticmethod
    def Sphere(center=(0, 0, 0), R=1, name='Sphere'):
        r"""
        Return a sphere of radius ``R`` centered at ``center``.

        INPUT:

        - ``center`` -- 3-tuple; center of the sphere

        - ``R`` -- radius of the sphere

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Sphere`.

        EXAMPLES::

            sage: sphere = surfaces.Sphere(center=(0, 1, -1), R=2); sphere
            Parametrized surface ('Sphere') with equation (2*cos(u)*cos(v), 2*cos(v)*sin(u) + 1, 2*sin(v) - 1)
            sage: sphere.plot()                                                         # needs sage.plot
            Graphics3d Object

        Note that the radius of the sphere can be negative. The surface thus
        obtained is equal to the sphere (or part thereof) with positive radius,
        whose coordinate functions have been multiplied by -1. Compare for
        instant the first octant of the unit sphere with positive radius::

            sage: octant1 = surfaces.Sphere(R=1); octant1
            Parametrized surface ('Sphere') with equation (cos(u)*cos(v), cos(v)*sin(u), sin(v))
            sage: octant1.plot((0, pi/2), (0, pi/2))                                    # needs sage.plot
            Graphics3d Object

        with the first octant of the unit sphere with negative radius::

            sage: octant2 = surfaces.Sphere(R=-1); octant2
            Parametrized surface ('Sphere') with equation (-cos(u)*cos(v), -cos(v)*sin(u), -sin(v))
            sage: octant2.plot((0, pi/2), (0, pi/2))                                    # needs sage.plot
            Graphics3d Object
        """
        return SurfaceGenerators.Ellipsoid(center, (R, R, R), name)

    @staticmethod
    def Torus(r=2, R=3, name='Torus'):
        r"""
        Return a torus obtained by revolving a circle of radius ``r`` around
        a coplanar axis ``R`` units away from the center of the circle. The
        parametrization used is

        .. MATH::

            \begin{aligned}
              x(u, v) & = (R + r \cos(v)) \cos(u); \\
              y(u, v) & = (R + r \cos(v)) \sin(u); \\
              z(u, v) & = r \sin(v).
            \end{aligned}

        INPUT:

        - ``r``, ``R`` -- minor and major radius of the torus

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Torus`.

        EXAMPLES::

            sage: torus = surfaces.Torus(); torus
            Parametrized surface ('Torus') with equation ((2*cos(v) + 3)*cos(u), (2*cos(v) + 3)*sin(u), 2*sin(v))
            sage: torus.plot()                                                          # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        torus_eq = [(R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v)]
        coords = ((u, 0, 2 * pi), (v, 0, 2 * pi))
        return ParametrizedSurface3D(torus_eq, coords, name)

    @staticmethod
    def WhitneyUmbrella(name="Whitney's umbrella"):
        r"""
        Return Whitney's umbrella, with parametric representation

        .. MATH::

            x(u, v) = uv, \quad y(u, v) = u, \quad z(u, v) = v^2.

        INPUT:

        - ``name`` -- string; name of the surface

        For more information, see :wikipedia:`Whitney_umbrella`.

        EXAMPLES::

            sage: whitney = surfaces.WhitneyUmbrella(); whitney
            Parametrized surface ('Whitney's umbrella') with equation (u*v, u, v^2)
            sage: whitney.plot()                                                        # needs sage.plot
            Graphics3d Object
        """
        u, v = var('u, v')
        whitney_eq = [u*v, u, v**2]
        coords = ((u, -1, 1), (v, -1, 1))
        return ParametrizedSurface3D(whitney_eq, coords, name)


# Easy access to the surface generators
surfaces = SurfaceGenerators()
