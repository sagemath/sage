"""
Differential Geometry of Parametrized Surfaces.

AUTHORS:
        - Mikhail Malakhaltsev (2010-09-25): initial version
        - Joris Vankerschaver  (2010-10-25): implementation, doctests

"""
#*****************************************************************************
#       Copyright (C) 2010  Mikhail Malakhaltsev <mikarm@gmail.com>
#       Copyright (C) 2010  Joris Vankerschaver <joris.vankerschaver@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from itertools import product

from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.calculus.functional import diff
from sage.functions.other import sqrt
from sage.misc.cachefunc import cached_method
from sage.symbolic.ring import SR
from sage.symbolic.constants import pi
from sage.symbolic.assumptions import assume

def _simplify_full_rad(f):
    """
    Helper function to conveniently call :meth:`simplify_full` and
    :meth:`simplify_radical` in succession.

    INPUT:

     - ``f`` - a symbolic expression.

    EXAMPLES::

        sage: from sage.geometry.riemannian_manifolds.parametrized_surface3d import _simplify_full_rad
        sage: _simplify_full_rad(sqrt(x^2)/x)
        1

    """
    return f.simplify_full().simplify_radical()


class ParametrizedSurface3D(SageObject):
    r"""
    Class representing a parametrized two-dimensional surface in
    Euclidian three-space.  Provides methods for calculating the main
    geometrical objects related to such a surface, such as the first
    and the second fundamental form, the total (Gaussian) and the mean
    curvature, the geodesic curves, parallel transport, etc.


    INPUT:

     - ``surface_equation`` -- a 3-tuple of functions specifying a parametric
       representation of the surface.

     - ``variables`` -- a 2-tuple of intrinsic coordinates `(u, v)` on the
       surface, with `u` and `v` symbolic variables, or a 2-tuple of triples
       $(u, u_{min}, u_{max})$,
       $(v, v_{min}, v_{max})$ when the parameter range
       for the coordinates is known.

     - ``name`` -- name of the surface (optional).


    .. note::

       Throughout the documentation, we use the Einstein summation
       convention: whenever an index appears twice, once as a
       subscript, and once as a superscript, summation over that index
       is implied.  For instance, `g_{ij} g^{jk}` stands for `\sum_j
       g_{ij}g^{jk}`.


    EXAMPLES:

    We give several examples of standard surfaces in differential
    geometry.  First, let's construct an elliptic paraboloid by
    explicitly specifying its parametric equation::

        sage: u, v = var('u,v', domain='real')
        sage: eparaboloid = ParametrizedSurface3D((u, v, u^2 + v^2), (u, v),'elliptic paraboloid'); eparaboloid
        Parametrized surface ('elliptic paraboloid') with equation (u, v, u^2 + v^2)

    When the ranges for the intrinsic coordinates are known, they can be
    specified explicitly. This is mainly useful for plotting. Here we
    construct half of an ellipsoid::

        sage: u1, u2 = var ('u1, u2', domain='real');
        sage: coords = ((u1, -pi/2, pi/2), (u2, 0, pi))
        sage: ellipsoid_eq = (cos(u1)*cos(u2), 2*sin(u1)*cos(u2), 3*sin(u2))
        sage: ellipsoid = ParametrizedSurface3D(ellipsoid_eq, coords, 'ellipsoid'); ellipsoid
        Parametrized surface ('ellipsoid') with equation (cos(u1)*cos(u2), 2*cos(u2)*sin(u1), 3*sin(u2))
        sage: ellipsoid.plot()

    Standard surfaces can be constructed using the ``surfaces`` generator::

        sage: klein = surfaces.Klein(); klein
        Parametrized surface ('Klein bottle') with equation (-(sin(1/2*u)*sin(2*v) - cos(1/2*u)*sin(v) - 1)*cos(u), -(sin(1/2*u)*sin(2*v) - cos(1/2*u)*sin(v) - 1)*sin(u), cos(1/2*u)*sin(2*v) + sin(1/2*u)*sin(v))

    Latex representation of the surfaces::

        sage: u, v = var('u, v', domain='real')
        sage: sphere = ParametrizedSurface3D((cos(u)*cos(v), sin(u)*cos(v), sin(v)), (u, v), 'sphere')
        sage: print latex(sphere)
        \left(\cos\left(u\right) \cos\left(v\right), \cos\left(v\right) \sin\left(u\right), \sin\left(v\right)\right)
        sage: print sphere._latex_()
        \left(\cos\left(u\right) \cos\left(v\right), \cos\left(v\right) \sin\left(u\right), \sin\left(v\right)\right)
        sage: print sphere
        Parametrized surface ('sphere') with equation (cos(u)*cos(v), cos(v)*sin(u), sin(v))

    To plot a parametric surface, use the :meth:`plot` member function::

        sage: enneper = surfaces.Enneper(); enneper
        Parametrized surface ('Enneper's surface') with equation (-1/9*(u^2 - 3*v^2 - 3)*u, -1/9*(3*u^2 - v^2 + 3)*v, 1/3*u^2 - 1/3*v^2)
        sage: enneper.plot(aspect_ratio='automatic')

    We construct an ellipsoid whose axes are given by symbolic variables `a`,
    `b` and `c`, and find the natural frame of tangent vectors,
    expressed in intrinsic coordinates. Note that the result is a
    dictionary of vector fields::

        sage: a, b, c = var('a, b, c', domain='real')
        sage: u1, u2 = var('u1, u2', domain='real')
        sage: ellipsoid_eq = (a*cos(u1)*cos(u2), b*sin(u1)*cos(u2), c*sin(u2))
        sage: ellipsoid = ParametrizedSurface3D(ellipsoid_eq, (u1, u2), 'Symbolic ellipsoid'); ellipsoid
        Parametrized surface ('Symbolic ellipsoid') with equation (a*cos(u1)*cos(u2), b*cos(u2)*sin(u1), c*sin(u2))

        sage: ellipsoid.natural_frame()
        {1: (-a*cos(u2)*sin(u1), b*cos(u1)*cos(u2), 0), 2: (-a*cos(u1)*sin(u2), -b*sin(u1)*sin(u2), c*cos(u2))}

    We find the normal vector field to the surface.  The normal vector
    field is the vector product of the vectors of the natural frame,
    and is given by::

        sage: ellipsoid.normal_vector()
        (b*c*cos(u1)*cos(u2)^2, a*c*cos(u2)^2*sin(u1), a*b*cos(u2)*sin(u2))

    By default, the normal vector field is not normalized.  To obtain
    the unit normal vector field of the elliptic paraboloid, we put::

        sage: u, v = var('u,v', domain='real')
        sage: eparaboloid = ParametrizedSurface3D([u,v,u^2+v^2],[u,v],'elliptic paraboloid')
        sage: eparaboloid.normal_vector(normalized=True)
        (-2*u/sqrt(4*u^2 + 4*v^2 + 1), -2*v/sqrt(4*u^2 + 4*v^2 + 1), 1/sqrt(4*u^2 + 4*v^2 + 1))

    Now let us compute the coefficients of the first fundamental form of the torus::

        sage: u, v = var('u, v', domain='real')
        sage: a, b = var('a, b', domain='real')
        sage: torus = ParametrizedSurface3D(((a + b*cos(u))*cos(v),(a + b*cos(u))*sin(v), b*sin(u)),[u,v],'torus')
        sage: torus.first_fundamental_form_coefficients()
        {(1, 2): 0, (1, 1): b^2, (2, 1): 0, (2, 2): b^2*cos(u)^2 + 2*a*b*cos(u) + a^2}

    The first fundamental form can be used to compute the length of a
    curve on the surface.  For example, let us find the length of the
    curve $u^1 = t$, $u^2 = t$, $t \in [0,2\pi]$, on the ellipsoid
    with axes $a=1$, $b=1.5$ and $c=1$. So we take the curve::

        sage: t = var('t', domain='real')
        sage: u1 = t
        sage: u2 = t

    Then find the tangent vector::

        sage: du1 = diff(u1,t)
        sage: du2 = diff(u2,t)
        sage: du = vector([du1, du2]); du
        (1, 1)

    Once we specify numerical values for the axes of the ellipsoid, we can
    determine the numerical value of the length integral::

        sage: L = sqrt(ellipsoid.first_fundamental_form(du, du).substitute(u1=u1,u2=u2))
        sage: print numerical_integral(L.substitute(a=2, b=1.5, c=1),0,1)[0]
        2.00127905972

    We find the area of the sphere of radius $R$::

        sage: R = var('R', domain='real');
        sage: u, v = var('u,v', domain='real');
        sage: assume(R>0)
        sage: assume(cos(v)>0)
        sage: sphere = ParametrizedSurface3D([R*cos(u)*cos(v),R*sin(u)*cos(v),R*sin(v)],[u,v],'sphere')
        sage: integral(integral(sphere.area_form(),u,0,2*pi),v,-pi/2,pi/2)
        4*pi*R^2

    We can find an orthonormal frame field $\{e_1, e_2\}$ of a surface
    and calculate its structure functions.  Let us first determine the
    orthonormal frame field for the elliptic paraboloid::

        sage: u, v = var('u,v', domain='real')
        sage: eparaboloid = ParametrizedSurface3D([u,v,u^2+v^2],[u,v],'elliptic paraboloid')
        sage: eparaboloid.orthonormal_frame()
        {1: (1/sqrt(4*u^2 + 1), 0, 2*u/sqrt(4*u^2 + 1)), 2: (-4*u*v/(sqrt(4*u^2 + 4*v^2 + 1)*sqrt(4*u^2 + 1)), sqrt(4*u^2 + 1)/sqrt(4*u^2 + 4*v^2 + 1), 2*v/(sqrt(4*u^2 + 4*v^2 + 1)*sqrt(4*u^2 + 1)))}

    We can express the orthogonal frame field both in exterior
    coordinates (i.e. expressed as vector field fields in the ambient
    space $\RR^3$, the default) or in intrinsic coordinates
    (with respect to the natural frame).  Here we use intrinsic
    coordinates::

        sage: eparaboloid.orthonormal_frame(coordinates='int')
        {1: (1/sqrt(4*u^2 + 1), 0), 2: (-4*u*v/(sqrt(4*u^2 + 4*v^2 + 1)*sqrt(4*u^2 + 1)), sqrt(4*u^2 + 1)/sqrt(4*u^2 + 4*v^2 + 1))}

    Using the orthonormal frame in interior coordinates, we can calculate
    the structure functions $c^k_{ij}$ of the surface, defined by
    $[e_i,e_j] =  c^k_{ij} e_k$, where $[e_i, e_j]$ represents the Lie
    bracket of two frame vector fields $e_i, e_j$.  For the
    elliptic paraboloid, we get::

        sage: EE = eparaboloid.orthonormal_frame(coordinates='int')
        sage: E1 = EE[1]; E2 = EE[2]
        sage: CC = eparaboloid.frame_structure_functions(E1,E2)
        sage: CC[1,2,1].simplify_full()
        4*sqrt(4*u^2 + 4*v^2 + 1)*v/((16*u^4 + 4*(4*u^2 + 1)*v^2 + 8*u^2 + 1)*sqrt(4*u^2 + 1))

    We compute the Gaussian and mean curvatures of the sphere::

        sage: sphere = surfaces.Sphere(); sphere
        Parametrized surface ('Sphere') with equation (cos(u)*cos(v), cos(v)*sin(u), sin(v))
        sage: K = sphere.gauss_curvature(); K # Not tested -- see trac 12737
        1
        sage: H = sphere.mean_curvature(); H # Not tested -- see trac 12737
        -1

    We can easily generate a color plot of the Gaussian curvature of a surface.
    Here we deal with the ellipsoid::

        sage: u1, u2 = var('u1,u2', domain='real');
        sage: u = [u1,u2]
        sage: ellipsoid_equation(u1,u2) = [2*cos(u1)*cos(u2),1.5*cos(u1)*sin(u2),sin(u1)]
        sage: ellipsoid = ParametrizedSurface3D(ellipsoid_equation(u1,u2), [u1, u2],'ellipsoid')
        sage: # set intervals for variables and the number of division points
        sage: u1min, u1max = -1.5, 1.5
        sage: u2min, u2max = 0, 6.28
        sage: u1num, u2num = 10, 20
        sage: # make the arguments array
        sage: from numpy import linspace
        sage: u1_array = linspace(u1min, u1max, u1num)
        sage: u2_array = linspace(u2min, u2max, u2num)
        sage: u_array = [ (uu1,uu2) for uu1 in u1_array for uu2 in u2_array]
        sage: # Find the gaussian curvature
        sage: K(u1,u2) = ellipsoid.gauss_curvature()
        sage: # Make array of K values
        sage: K_array = [K(uu[0],uu[1]) for uu in u_array]
        sage: # Find minimum and max of the gauss curvature
        sage: K_max = max(K_array)
        sage: K_min = min(K_array)
        sage: # Make the array of color coefficients
        sage: cc_array = [ (ccc - K_min)/(K_max - K_min) for ccc in K_array ]
        sage: points_array = [ellipsoid_equation(u_array[counter][0],u_array[counter][1]) for counter in range(0,len(u_array)) ]
        sage: curvature_ellipsoid_plot = sum( point([xx for xx in points_array[counter]],color=hue(cc_array[counter]/2))  for counter in range(0,len(u_array)) )
        sage: curvature_ellipsoid_plot.show(aspect_ratio=1)

    We can find the principal curvatures and principal directions of the
    elliptic paraboloid::

        sage: u, v = var('u, v', domain='real')
        sage: eparaboloid = ParametrizedSurface3D([u, v, u^2+v^2], [u, v], 'elliptic paraboloid')
        sage: pd = eparaboloid.principal_directions(); pd
        [(2*sqrt(4*u^2 + 4*v^2 + 1)/(16*u^4 + 16*v^4 + 8*(4*u^2 + 1)*v^2 + 8*u^2 + 1), [(1, v/u)], 1), (2/sqrt(4*u^2 + 4*v^2 + 1), [(1, -u/v)], 1)]

    We extract the principal curvatures::

        sage: k1 = pd[0][0].simplify_full()
        sage: k1
        2*sqrt(4*u^2 + 4*v^2 + 1)/(16*u^4 + 16*v^4 + 8*(4*u^2 + 1)*v^2 + 8*u^2 + 1)
        sage: k2 = pd[1][0].simplify_full()
        sage: k2
        2/sqrt(4*u^2 + 4*v^2 + 1)

    and check them by comparison with the Gaussian and mean curvature
    expressed in terms of the principal curvatures::

        sage: K = eparaboloid.gauss_curvature().simplify_full()
        sage: K
        4/(16*u^4 + 16*v^4 + 8*(4*u^2 + 1)*v^2 + 8*u^2 + 1)
        sage: H = eparaboloid.mean_curvature().simplify_full()
        sage: H
        2*(2*u^2 + 2*v^2 + 1)/(4*u^2 + 4*v^2 + 1)^(3/2)
        sage: (K - k1*k2).simplify_full()
        0
        sage: (2*H - k1 - k2).simplify_full()
        0

    We can find the intrinsic (local coordinates) of the principal directions::

        sage: pd[0][1]
        [(1, v/u)]
        sage: pd[1][1]
        [(1, -u/v)]

    The ParametrizedSurface3D class also contains functionality to
    compute the coefficients of the second fundamental form, the shape
    operator, the rotation on the surface at a given angle, the
    connection coefficients.  One can also calculate numerically the
    geodesics and the parallel translation along a curve.

    Here we compute a number of geodesics on the sphere emanating
    from the point ``(1, 0, 0)``, in various directions. The geodesics
    intersect again in the antipodal point ``(-1, 0, 0)``, indicating
    that these points are conjugate::

        sage: S = surfaces.Sphere()
        sage: g1 = [c[-1] for c in S.geodesics_numerical((0,0),(1,0),(0,2*pi,100))]
        sage: g2 = [c[-1] for c in S.geodesics_numerical((0,0),(cos(pi/3),sin(pi/3)),(0,2*pi,100))]
        sage: g3 = [c[-1] for c in S.geodesics_numerical((0,0),(cos(2*pi/3),sin(2*pi/3)),(0,2*pi,100))]
        sage: (S.plot(opacity=0.3) + line3d(g1,color='red') + line3d(g2,color='red') + line3d(g3,color='red')).show()

    """


    def __init__(self, equation, variables, name=None):
        r"""
        See ``ParametrizedSurface3D`` for full documentation.

        .. note::

            The orientation of the surface is determined by the
            parametrization, that is, the natural frame with positive
            orientation is given by `\partial_1 \vec r`, `\partial_2 \vec
            r`.


        EXAMPLES::

            sage: u, v = var('u,v', domain='real')
            sage: eq = (3*u + 3*u*v^2 - u^3, 3*v + 3*u^2*v - v^3, 3*(u^2-v^2))
            sage: enneper = ParametrizedSurface3D(eq, (u, v),'Enneper Surface'); enneper
            Parametrized surface ('Enneper Surface') with equation (-u^3 + 3*u*v^2 + 3*u, 3*u^2*v - v^3 + 3*v, 3*u^2 - 3*v^2)

        """
        self.equation = tuple(equation)

        if len(variables[0]) > 0:
            self.variables_range = (variables[0][1:3], variables[1][1:3])
            self.variables_list  = (variables[0][0], variables[1][0])
        else:
            self.variables_range = None
            self.variables_list = variables

        self.variables = {1:self.variables_list[0],2:self.variables_list[1]}
        self.name = name


    def _latex_(self):
        r"""
        Returns the LaTeX representation of this parametrized surface.

        EXAMPLES::

            sage: u, v = var('u, v')
            sage: sphere = ParametrizedSurface3D((cos(u)*cos(v), sin(u)*cos(v), sin(v)), (u, v),'sphere')
            sage: latex(sphere)
            \left(\cos\left(u\right) \cos\left(v\right), \cos\left(v\right) \sin\left(u\right), \sin\left(v\right)\right)
            sage: sphere._latex_()
            \left(\cos\left(u\right) \cos\left(v\right), \cos\left(v\right) \sin\left(u\right), \sin\left(v\right)\right)

        """
        from sage.misc.latex import latex
        return latex(self.equation)


    def _repr_(self):
        r"""
        Returns the string representation of this parametrized surface.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: eq = (3*u + 3*u*v^2 - u^3, 3*v + 3*u^2*v - v^3, 3*(u^2-v^2))
            sage: enneper = ParametrizedSurface3D(eq,[u,v],'enneper_surface')
            sage: print enneper
            Parametrized surface ('enneper_surface') with equation (-u^3 + 3*u*v^2 + 3*u, 3*u^2*v - v^3 + 3*v, 3*u^2 - 3*v^2)
            sage: enneper._repr_()
            "Parametrized surface ('enneper_surface') with equation (-u^3 + 3*u*v^2 + 3*u, 3*u^2*v - v^3 + 3*v, 3*u^2 - 3*v^2)"

        """
        name = 'Parametrized surface'
        if self.name is not None:
            name += " ('%s')" % self.name
        s ='%(designation)s with equation %(eq)s' % \
            {'designation': name, 'eq': str(self.equation)}
        return s


    def point(self, coords):
        r"""
        Returns a point on the surface given its intrinsic coordinates.

        INPUT:

         - ``coords`` - 2-tuple specifying the intrinsic coordinates ``(u, v)`` of the point.

        OUTPUT:

         - 3-vector specifying the coordinates in `\RR^3` of the point.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: torus = ParametrizedSurface3D(((2 + cos(u))*cos(v),(2 + cos(u))*sin(v), sin(u)),[u,v],'torus')
            sage: torus.point((0, pi/2))
            (0, 3, 0)
            sage: torus.point((pi/2, pi))
            (-2, 0, 1)
            sage: torus.point((pi, pi/2))
            (0, 1, 0)

        """

        d = dict(zip(self.variables_list, coords))
        return vector([f.subs(d) for f in self.equation])


    def tangent_vector(self, coords, components):
        r"""
        Returns the components of a tangent vector given the intrinsic
        coordinates of the base point and the components of the vector
        in the intrinsic frame.

        INPUT:

         - ``coords`` - 2-tuple specifying the intrinsic coordinates ``(u, v)`` of the point.

         - ``components`` - 2-tuple specifying the components of the tangent vector in the intrinsic coordinate frame.

        OUTPUT:

         - 3-vector specifying the components in `\RR^3` of the vector.

        EXAMPLES:

        We compute two tangent vectors to Enneper's surface along the
        coordinate lines and check that their cross product gives the
        normal vector::

            sage: u, v = var('u,v', domain='real')
            sage: eq = (3*u + 3*u*v^2 - u^3, 3*v + 3*u^2*v - v^3, 3*(u^2-v^2))
            sage: e = ParametrizedSurface3D(eq, (u, v),'Enneper Surface')

            sage: w1 = e.tangent_vector((1, 2), (1, 0)); w1
            (12, 12, 6)
            sage: w2 = e.tangent_vector((1, 2), (0, 1)); w2
            (12, -6, -12)
            sage: w1.cross_product(w2)
            (-108, 216, -216)

            sage: n = e.normal_vector().subs({u: 1, v: 2}); n
            (-108, 216, -216)
            sage: n == w1.cross_product(w2)
            True

        """

        components = vector(components)
        d = dict(zip(self.variables_list, coords))
        jacobian = matrix([[f.diff(u).subs(d) for u in self.variables_list]
                           for f in self.equation])
        return jacobian * components


    def plot(self, urange=None, vrange=None, **kwds):
        r"""
        Enable easy plotting directly from the surface class.

        The optional keywords ``urange`` and ``vrange`` specify the range for
        the surface parameters `u` and `v`. If either of these parameters
        is ``None``, the method checks whether a parameter range was
        specified when the surface was created. If not, the default of
        $(0, 2 \pi)$ is used.

        INPUT:

         - ``urange`` - 2-tuple specifying the parameter range for `u`.
         - ``vrange`` - 2-tuple specifying the parameter range for `v`.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: eq = (3*u + 3*u*v^2 - u^3, 3*v + 3*u^2*v - v^3, 3*(u^2-v^2))
            sage: enneper = ParametrizedSurface3D(eq, (u, v), 'Enneper Surface')
            sage: enneper.plot((-5, 5), (-5, 5))

        """

        from sage.plot.plot3d.parametric_plot3d import parametric_plot3d

        if self.variables_range is None:
            if urange is None:
                urange = (0, 2*pi)
            if vrange is None:
                vrange = (0, 2*pi)
        else:
            if urange is None:
                urange = self.variables_range[0]
            if vrange is None:
                vrange = self.variables_range[1]

        urange3 = (self.variables[1],) + tuple(urange)
        vrange3 = (self.variables[2],) + tuple(vrange)
        P = parametric_plot3d(self.equation, urange3, vrange3, **kwds)

        return P


    @cached_method
    def natural_frame(self):
        """
        Returns the natural tangent frame on the parametrized surface.
        The vectors of this frame are tangent to the coordinate lines
        on the surface.

        OUTPUT:

        - The natural frame as a dictionary.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: eparaboloid = ParametrizedSurface3D((u, v, u^2+v^2), (u, v), 'elliptic paraboloid')
            sage: eparaboloid.natural_frame()
            {1: (1, 0, 2*u), 2: (0, 1, 2*v)}
        """

        dr1 = \
            vector([_simplify_full_rad( diff(f,self.variables[1]) )
                    for f in self.equation])
        dr2 = \
            vector([_simplify_full_rad( diff(f,self.variables[2]) )
                    for f in self.equation])

        return {1:dr1, 2:dr2}


    @cached_method
    def normal_vector(self, normalized=False):
        """
        Returns the normal vector field of the parametrized surface.

        INPUT:

          - ``normalized`` - default ``False`` - specifies whether the normal vector should be normalized.

        OUTPUT:

         - Normal vector field.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: eparaboloid = ParametrizedSurface3D((u, v, u^2 + v^2), (u, v), 'elliptic paraboloid')
            sage: eparaboloid.normal_vector(normalized=False)
            (-2*u, -2*v, 1)
            sage: eparaboloid.normal_vector(normalized=True)
            (-2*u/sqrt(4*u^2 + 4*v^2 + 1), -2*v/sqrt(4*u^2 + 4*v^2 + 1), 1/sqrt(4*u^2 + 4*v^2 + 1))

        """

        dr = self.natural_frame()
        normal = dr[1].cross_product(dr[2])

        if normalized:
            normal /= normal.norm()
        return _simplify_full_rad(normal)


    @cached_method
    def _compute_first_fundamental_form_coefficient(self, index):
        """
        Helper function to compute coefficients of the first fundamental form.

        Do not call this method directly; instead use
        ``first_fundamental_form_coefficient``.
        This method is cached, and expects its argument to be a list.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: eparaboloid = ParametrizedSurface3D((u, v, u^2+v^2), (u, v))
            sage: eparaboloid._compute_first_fundamental_form_coefficient((1,2))
            4*u*v

        """
        dr = self.natural_frame()
        return _simplify_full_rad(dr[index[0]]*dr[index[1]])


    def first_fundamental_form_coefficient(self, index):
        r"""
        Compute a single component $g_{ij}$ of the first fundamental form.  If
        the parametric representation of the surface is given by the vector
        function $\vec r(u^i)$, where $u^i$, $i = 1, 2$ are curvilinear
        coordinates, then $g_{ij} = \frac{\partial \vec r}{\partial u^i} \cdot \frac{\partial \vec r}{\partial u^j}$.

        INPUT:

         - ``index`` - tuple ``(i, j)`` specifying the index of the component $g_{ij}$.

        OUTPUT:

         - Component $g_{ij}$ of the first fundamental form

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: eparaboloid = ParametrizedSurface3D((u, v, u^2+v^2), (u, v))
            sage: eparaboloid.first_fundamental_form_coefficient((1,2))
            4*u*v

        When the index is invalid, an error is raised::

            sage: u, v = var('u, v', domain='real')
            sage: eparaboloid = ParametrizedSurface3D((u, v, u^2+v^2), (u, v))
            sage: eparaboloid.first_fundamental_form_coefficient((1,5))
            Traceback (most recent call last):
            ...
            ValueError: Index (1, 5) out of bounds.

        """
        index = tuple(sorted(index))
        if len(index) == 2 and all(i == 1 or i == 2 for i in index):
            return self._compute_first_fundamental_form_coefficient(index)
        else:
            raise ValueError("Index %s out of bounds." % str(index))

    def first_fundamental_form_coefficients(self):
        r"""
        Returns the coefficients of the first fundamental form as a dictionary.
        The keys are tuples $(i, j)$, where $i$ and $j$ range over $1, 2$,
        while the values are the corresponding coefficients $g_{ij}$.

        OUTPUT:

         - Dictionary of first fundamental form coefficients.

        EXAMPLES::

            sage: u, v = var('u,v', domain='real')
            sage: sphere = ParametrizedSurface3D((cos(u)*cos(v), sin(u)*cos(v), sin(v)), (u, v), 'sphere')
            sage: sphere.first_fundamental_form_coefficients()
            {(1, 2): 0, (1, 1): cos(v)^2, (2, 1): 0, (2, 2): 1}

        """
        coefficients = {}
        for index in product((1, 2), repeat=2):
            sorted_index = list(sorted(index))
            coefficients[index] = \
                self._compute_first_fundamental_form_coefficient(index)
        return coefficients


    def first_fundamental_form(self, vector1, vector2):
        r"""
        Evaluate the first fundamental form on two vectors expressed with
        respect to the natural coordinate frame on the surface. In other words,
        if the vectors are $v = (v^1, v^2)$ and $w = (w^1, w^2)$, calculate
        $g_{11} v^1 w^1 + g_{12}(v^1 w^2 + v^2 w^1) + g_{22} v^2 w^2$, with
        $g_{ij}$ the coefficients of the first fundamental form.

        INPUT:

         - ``vector1``, ``vector2`` - vectors on the surface.

        OUTPUT:

         - First fundamental form evaluated on the input vectors.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: v1, v2, w1, w2 = var('v1, v2, w1, w2', domain='real')
            sage: sphere = ParametrizedSurface3D((cos(u)*cos(v), sin(u)*cos(v), sin(v)), (u, v),'sphere')
            sage: sphere.first_fundamental_form(vector([v1,v2]),vector([w1,w2]))
            v1*w1*cos(v)^2 + v2*w2

            sage: vv = vector([1,2])
            sage: sphere.first_fundamental_form(vv,vv)
            cos(v)^2 + 4

            sage: sphere.first_fundamental_form([1,1],[2,1])
            2*cos(v)^2 + 1
        """
        gamma = self.first_fundamental_form_coefficients()
        return sum(gamma[(i,j)] * vector1[i - 1] * vector2[j - 1]
                   for i, j in product((1, 2), repeat=2))


    def area_form_squared(self):
        """
        Returns the square of the coefficient of the area form on the surface.
        In terms of the coefficients $g_{ij}$ (where $i, j = 1, 2$) of the
        first fundamental form, this invariant is given by
        $A^2 = g_{11}g_{22} - g_{12}^2$.

        See also :meth:`.area_form`.

        OUTPUT:

         - Square of the area form

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.area_form_squared()
            cos(v)^2

        """
        gamma = self.first_fundamental_form_coefficients()
        sq = gamma[(1,1)] * gamma[(2,2)] - gamma[(1,2)]**2
        return _simplify_full_rad(sq)


    def area_form(self):
        r"""
        Returns the coefficient of the area form on the surface.  In terms of
        the coefficients $g_{ij}$ (where $i, j = 1, 2$) of the first
        fundamental form, the coefficient of the area form is given by
        $A = \sqrt{g_{11}g_{22} - g_{12}^2}$.

        See also :meth:`.area_form_squared`.

        OUTPUT:

         - Coefficient of the area form

        EXAMPLES::

            sage: u, v = var('u,v', domain='real')
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.area_form()
            cos(v)

        """
        f = abs(sqrt(self.area_form_squared()))
        return _simplify_full_rad(f)


    def first_fundamental_form_inverse_coefficients(self):
        r"""
        Returns the coefficients $g^{ij}$ of the inverse of the fundamental
        form, as a dictionary.  The inverse coefficients are defined by
        $g^{ij} g_{jk} = \delta^i_k$  with $\delta^i_k$ the Kronecker
        delta.

        OUTPUT:

         - Dictionary of the inverse coefficients.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.first_fundamental_form_inverse_coefficients()
            {(1, 2): 0, (1, 1): cos(v)^(-2), (2, 1): 0, (2, 2): 1}

        """

        g = self.first_fundamental_form_coefficients()
        D = g[(1,1)] * g[(2,2)] - g[(1,2)]**2

        gi11 = _simplify_full_rad(g[(2,2)]/D)
        gi12 = _simplify_full_rad(-g[(1,2)]/D)
        gi21 = gi12
        gi22 = _simplify_full_rad(g[(1,1)]/D)

        return {(1,1): gi11, (1,2): gi12, (2,1): gi21, (2,2): gi22}


    def first_fundamental_form_inverse_coefficient(self, index):
        r"""
        Returns a specific component $g^{ij}$ of the inverse of the fundamental
        form.

        INPUT:

         - ``index`` - tuple ``(i, j)`` specifying the index of the component $g^{ij}$.

        OUTPUT:

         - Component of the inverse of the fundamental form.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.first_fundamental_form_inverse_coefficient((1, 2))
            0
            sage: sphere.first_fundamental_form_inverse_coefficient((1, 1))
            cos(v)^(-2)

        """

        index = tuple(sorted(index))
        if len(index) == 2 and all(i == 1 or i == 2 for i in index):
            return self.first_fundamental_form_inverse_coefficients()[index]
        else:
            raise ValueError("Index %s out of bounds." % str(index))



    @cached_method
    def rotation(self,theta):
        r"""
        Gives the matrix of the rotation operator over a given angle $\theta$
        with respect to the natural frame.

        INPUT:

         - ``theta`` - rotation angle

        OUTPUT:

         - Rotation matrix with respect to the natural frame.

        ALGORITHM:

        The operator of rotation over $\pi/2$ is $J^i_j = g^{ik}\omega_{jk}$,
        where $\omega$ is the area form.  The operator of rotation over an
        angle $\theta$ is $\cos(\theta) I + sin(\theta) J$.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: assume(cos(v)>0)
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')

        We first compute the matrix of rotation over $\pi/3$::

            sage: rotation = sphere.rotation(pi/3); rotation
            [                1/2 -1/2*sqrt(3)/cos(v)]
            [ 1/2*sqrt(3)*cos(v)                 1/2]

        We verify that three succesive rotations over $\pi/3$ yield minus the identity::

            sage: rotation^3
            [-1  0]
            [ 0 -1]

        """

        from sage.functions.trig import sin, cos

        gi = self.first_fundamental_form_inverse_coefficients()
        w12 = self.area_form()
        R11 = (cos(theta) + sin(theta)*gi[1,2]*w12).simplify_full()
        R12 = (- sin(theta)*gi[1,1]*w12).simplify_full()
        R21 = (sin(theta)*gi[2,2]*w12).simplify_full()
        R22 = (cos(theta) - sin(theta)*gi[2,1]*w12).simplify_full()
        return matrix([[R11,R12],[R21,R22]])


    @cached_method
    def orthonormal_frame(self, coordinates='ext'):
        r"""
        Returns the orthonormal frame field on the surface, expressed either
        in exterior coordinates (i.e. expressed as vector fields in the
        ambient space $\mathbb{R}^3$, the default) or interior coordinates
        (with respect to the natural frame)

        INPUT:

         - ``coordinates`` - either ``ext`` (default) or ``int``.

        OUTPUT:

         - Orthogonal frame field as a dictionary.

        ALGORITHM:

        We normalize the first vector $\vec e_1$ of the natural frame and then
        get the second frame vector as $\vec e_2 = [\vec n, \vec e_1]$, where
        $\vec n$ is the unit normal to the surface.

        EXAMPLES::

            sage: u, v = var('u,v', domain='real')
            sage: assume(cos(v)>0)
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v), sin(u)*cos(v), sin(v)], [u, v],'sphere')
            sage: frame = sphere.orthonormal_frame(); frame
            {1: (-sin(u), cos(u), 0), 2: (-cos(u)*sin(v), -sin(u)*sin(v), cos(v))}
            sage: (frame[1]*frame[1]).simplify_full()
            1
            sage: (frame[1]*frame[2]).simplify_full()
            0
            sage: frame[1] == sphere.orthonormal_frame_vector(1)
            True

        We compute the orthonormal frame with respect to the natural frame on
        the surface::

            sage: frame_int = sphere.orthonormal_frame(coordinates='int'); frame_int
            {1: (1/cos(v), 0), 2: (0, 1)}
            sage: sphere.first_fundamental_form(frame_int[1], frame_int[1])
            1
            sage: sphere.first_fundamental_form(frame_int[1], frame_int[2])
            0
            sage: sphere.first_fundamental_form(frame_int[2], frame_int[2])
            1

        """


        from sage.symbolic.constants import pi

        if coordinates not in ['ext', 'int']:
            raise ValueError("Coordinate system must be exterior ('ext') "
                             "or interior ('int').")

        c  = self.first_fundamental_form_coefficient([1,1])
        if coordinates == 'ext':
            f1 = self.natural_frame()[1]

            E1 = _simplify_full_rad(f1/sqrt(c))
            E2 = _simplify_full_rad(
                self.normal_vector(normalized=True).cross_product(E1))
        else:
            E1 =  vector([_simplify_full_rad(1/sqrt(c)), 0])
            E2 = (self.rotation(pi/2)*E1).simplify_full()
        return  {1:E1, 2:E2}


    def orthonormal_frame_vector(self, index, coordinates='ext'):
        r"""
        Returns a specific basis vector field of the orthonormal frame field on
        the surface, expressed in exterior or interior coordinates.  See
        :meth:`orthogonal_frame` for more details.

        INPUT:

         - ``index`` - index of the basis vector;
         - ``coordinates`` - either ``ext`` (default) or ``int``.

        OUTPUT:

         - Orthonormal frame vector field.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: assume(cos(v)>0)
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: V1 = sphere.orthonormal_frame_vector(1); V1
            (-sin(u), cos(u), 0)
            sage: V2 = sphere.orthonormal_frame_vector(2); V2
            (-cos(u)*sin(v), -sin(u)*sin(v), cos(v))
            sage: (V1*V1).simplify_full()
            1
            sage: (V1*V2).simplify_full()
            0

            sage: n = sphere.normal_vector(normalized=True)
            sage: (V1.cross_product(V2) - n).simplify_full()
            (0, 0, 0)
        """

        return self.orthonormal_frame(coordinates)[index]


    def lie_bracket(self, v, w):
        r"""
        Returns the Lie bracket of two vector fields that are tangent
        to the surface. The vector fields should be given in intrinsic
        coordinates, i.e. with respect to the natural frame.

        INPUT:

         - ``v`` and ``w`` - vector fields on the surface, expressed
           as pairs of functions or as vectors of length 2.

        OUTPUT:

         - The Lie bracket $[v, w]$.


        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: assume(cos(v)>0)
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.lie_bracket([u,v],[-v,u])
            (0, 0)

            sage: EE_int = sphere.orthonormal_frame(coordinates='int')
            sage: sphere.lie_bracket(EE_int[1],EE_int[2])
            (sin(v)/cos(v)^2, 0)
        """
        v = vector(SR, v)
        w = vector(SR, w)

        variables = self.variables_list
        Dv = matrix([[_simplify_full_rad(diff(component, u))
                      for u in variables] for component in v])
        Dw = matrix([[_simplify_full_rad(diff(component, u))
                      for u in variables] for component in w])
        return vector(Dv*w - Dw*v).simplify_full()


    def frame_structure_functions(self, e1, e2):
        r"""
        Returns the structure functions $c^k_{ij}$ for a frame field
        $e_1, e_2$, i.e. a pair of vector fields on the surface which are
        linearly independent at each point.  The structure functions are
        defined using the Lie bracket by $[e_i,e_j] = c^k_{ij}e_k$.

        INPUT:

         - ``e1``, ``e2`` - vector fields in intrinsic coordinates on
           the surface, expressed as pairs of functions, or as vectors of
           length 2.

        OUTPUT:

         - Dictionary of structure functions, where the key ``(i, j, k)`` refers to
           the structure function $c_{i,j}^k$.


        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: assume(cos(v) > 0)
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v), sin(u)*cos(v), sin(v)], [u, v], 'sphere')
            sage: sphere.frame_structure_functions([u, v], [-v, u])
            {(1, 2, 1): 0, (2, 1, 2): 0, (2, 2, 2): 0, (1, 2, 2): 0, (1, 1, 1): 0, (2, 1, 1): 0, (2, 2, 1): 0, (1, 1, 2): 0}

        We construct the structure functions of the orthonormal frame on the
        surface::

            sage: EE_int = sphere.orthonormal_frame(coordinates='int')
            sage: CC = sphere.frame_structure_functions(EE_int[1],EE_int[2]); CC
            {(1, 2, 1): sin(v)/cos(v), (2, 1, 2): 0, (2, 2, 2): 0, (1, 2, 2): 0, (1, 1, 1): 0, (2, 1, 1): -sin(v)/cos(v), (2, 2, 1): 0, (1, 1, 2): 0}
            sage: sphere.lie_bracket(EE_int[1],EE_int[2]) - CC[(1,2,1)]*EE_int[1] - CC[(1,2,2)]*EE_int[2]
            (0, 0)
            """
        e1 = vector(SR, e1)
        e2 = vector(SR, e2)

        lie_bracket = self.lie_bracket(e1, e2).simplify_full()
        transformation = matrix(SR, [e1, e2]).transpose()

        w = (transformation.inverse()*lie_bracket).simplify_full()

        return {(1,1,1): 0, (1,1,2): 0, (1,2,1): w[0], (1,2,2): w[1],
                (2,1,1): -w[0], (2,1,2): -w[1], (2,2,1): 0, (2,2,2): 0}


    @cached_method
    def _compute_second_order_frame_element(self, index):
        """
        Compute an element of the second order frame of the surface.  See
        :meth:`second_order_natural_frame` for more details.

        This method expects its arguments in tuple form for caching.
        As it does no input checking, it should not be called directly.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: paraboloid = ParametrizedSurface3D([u, v, u^2 + v^2], [u,v], 'paraboloid')
            sage: paraboloid._compute_second_order_frame_element((1, 2))
            (0, 0, 0)
            sage: paraboloid._compute_second_order_frame_element((2, 2))
            (0, 0, 2)

        """
        variables = [self.variables[i] for i in index]
        ddr_element = vector([_simplify_full_rad(diff(f, variables))
                              for f in self.equation])

        return ddr_element


    def second_order_natural_frame(self):
        r"""
        Returns the second-order frame of the surface, i.e. computes the
        second-order derivatives (with respect to the parameters on the
        surface) of the parametric expression $\vec r = \vec r(u^1,u^2)$
        of the surface.

        OUTPUT:

         - Dictionary where the keys are 2-tuples ``(i, j)`` and the values are the corresponding derivatives $r_{ij}$.

        EXAMPLES:

        We compute the second-order natural frame of the sphere::

            sage: u, v = var('u, v', domain='real')
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.second_order_natural_frame()
            {(1, 2): (sin(u)*sin(v), -cos(u)*sin(v), 0), (1, 1): (-cos(u)*cos(v), -cos(v)*sin(u), 0), (2, 1): (sin(u)*sin(v), -cos(u)*sin(v), 0), (2, 2): (-cos(u)*cos(v), -cos(v)*sin(u), -sin(v))}

        """

        vectors = {}
        for index in product((1, 2), repeat=2):
            sorted_index = tuple(sorted(index))
            vectors[index] = \
                self._compute_second_order_frame_element(sorted_index)
        return vectors


    def second_order_natural_frame_element(self, index):
        r"""
        Returns a vector in the second-order frame of the surface, i.e.
        computes the second-order derivatives of the parametric expression
        $\vec{r}$ of the surface with respect to the parameters listed in the
        argument.

        INPUT:

         - ``index`` - a 2-tuple ``(i, j)`` specifying the element of the second-order frame.

        OUTPUT:

         - The second-order derivative $r_{ij}$.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.second_order_natural_frame_element((1, 2))
            (sin(u)*sin(v), -cos(u)*sin(v), 0)

        """

        index = tuple(sorted(index))
        if len(index) == 2 and all(i == 1 or i == 2 for i in index):
            return self._compute_second_order_frame_element(index)
        else:
            raise ValueError("Index %s out of bounds." % str(index))

    @cached_method
    def _compute_second_fundamental_form_coefficient(self, index):
        """
        Compute a coefficient of the second fundamental form of the surface.
        See ``second_fundamental_form_coefficient`` for more details.

        This method expects its arguments in tuple form for caching.  As it
        does no input checking, it should not be called directly.

        EXAMPLES::

            sage: u, v = var('u,v', domain='real')
            sage: paraboloid = ParametrizedSurface3D([u, v, u^2+v^2], [u, v], 'paraboloid')
            sage: paraboloid._compute_second_fundamental_form_coefficient((1,1))
            2/sqrt(4*u^2 + 4*v^2 + 1)

        """
        N = self.normal_vector(normalized=True)
        v = self.second_order_natural_frame_element(index)
        return _simplify_full_rad(v*N)


    def second_fundamental_form_coefficient(self, index):
        r"""
        Returns the coefficient $h_{ij}$ of the second fundamental form
        corresponding to the index $(i, j)$.  If the equation of the surface
        is $\vec{r}(u^1, u^2)$, then $h_{ij} = \vec{r}_{u^i u^j} \cdot \vec{n}$,
        where $\vec{n}$ is the unit normal.

        INPUT:

         - ``index`` - a 2-tuple ``(i, j)``

        OUTPUT:

         - Component $h_{ij}$ of the second fundamental form.

        EXAMPLES::

            sage: u, v = var('u,v', domain='real')
            sage: assume(cos(v)>0)
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.second_fundamental_form_coefficient((1, 1))
            -cos(v)^2
            sage: sphere.second_fundamental_form_coefficient((2, 1))
            0

        """
        index = tuple(index)
        if len(index) == 2 and all(i == 1 or i == 2 for i in index):
            return self._compute_second_fundamental_form_coefficient(index)
        else:
            raise ValueError("Index %s out of bounds." % str(index))


    def second_fundamental_form_coefficients(self):
        """
        Returns the coefficients $h_{ij}$ of the second fundamental form as
        a dictionary, where the keys are the indices $(i, j)$ and the values
        are the corresponding components $h_{ij}$.

        When only one component is needed, consider instead the function
        :meth:`second_fundamental_form_coefficient`.

        OUTPUT:

        Dictionary of second fundamental form coefficients.

        EXAMPLES::

            sage: u, v = var('u, v', domain='real')
            sage: assume(cos(v)>0)
            sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
            sage: sphere.second_fundamental_form_coefficients()
            {(1, 2): 0, (1, 1): -cos(v)^2, (2, 1): 0, (2, 2): -1}

        """

        coefficients = {}
        for index in product((1, 2), repeat=2):
            coefficients[index] = \
                self._compute_second_fundamental_form_coefficient(index)
        return coefficients


    def second_fundamental_form(self,vector1,vector2):
        r"""
        Evaluates the second fundamental form on two vectors on the surface.
        If the vectors are given by $v=(v^1,v^2)$ and $w=(w^1,w^2)$, the
        result of this function is $h_{11} v^1 w^1 + h_{12}(v^1 w^2 + v^2 w^1) + h_{22} v^2 w^2$.

        INPUT:

         - ``vector1``, ``vector2`` - 2-tuples representing the input vectors.

        OUTPUT:

         - Value of the second fundamental form evaluated on the given vectors.

        EXAMPLES:

        We evaluate the second fundamental form on two symbolic vectors::

           sage: u, v = var('u, v', domain='real')
           sage: v1, v2, w1, w2 = var('v1, v2, w1, w2', domain='real')
           sage: assume(cos(v) > 0)
           sage: sphere = ParametrizedSurface3D([cos(u)*cos(v),sin(u)*cos(v),sin(v)],[u,v],'sphere')
           sage: sphere.second_fundamental_form(vector([v1, v2]), vector([w1, w2]))
           -v1*w1*cos(v)^2 - v2*w2

        We evaluate the second fundamental form on vectors with numerical
        components::

           sage: vect = vector([1,2])
           sage: sphere.second_fundamental_form(vect, vect)
           -cos(v)^2 - 4
           sage: sphere.second_fundamental_form([1,1], [2,1])
           -2*cos(v)^2 - 1

        """
        hh = self.second_fundamental_form_coefficients()
        return sum(hh[(i, j)] * vector1[i - 1] * vector2[j - 1]
                   for (i, j) in product((1, 2), repeat=2))


    def gauss_curvature(self):
        r"""
        Finds the gaussian curvature of the surface, given by
        $K = \frac{h_{11}h_{22} - h_{12}^2}{g_{11}g_{22} - g_{12}^2}$,
        where $g_{ij}$ and $h_{ij}$ are the coefficients of the first
        and second fundamental form, respectively.

        OUTPUT:

         - Gaussian curvature of the surface.

        EXAMPLES::

           sage: R = var('R')
           sage: assume(R>0)
           sage: u, v = var('u,v', domain='real')
           sage: assume(cos(v)>0)
           sage: sphere = ParametrizedSurface3D([R*cos(u)*cos(v),R*sin(u)*cos(v),R*sin(v)],[u,v],'sphere')
           sage: sphere.gauss_curvature()
           R^(-2)

        """
        hh = self.second_fundamental_form_coefficients()
        return _simplify_full_rad(
            (hh[(1,1)] * hh[(2,2)] - hh[(1,2)]**2)/self.area_form_squared())


    def mean_curvature(self):
        r"""
        Finds the mean curvature of the surface, given by
        $H = \frac{1}{2}\frac{g_{22}h_{11} - 2g_{12}h_{12} + g_{11}h_{22}}{g_{11}g_{22} - g_{12}^2}$,
        where $g_{ij}$ and $h_{ij}$ are the components of the first and second
        fundamental forms, respectively.

        OUTPUT:

         - Mean curvature of the surface

        EXAMPLES::

           sage: R = var('R')
           sage: assume(R>0)
           sage: u, v = var('u,v', domain='real')
           sage: assume(cos(v)>0)
           sage: sphere = ParametrizedSurface3D([R*cos(u)*cos(v),R*sin(u)*cos(v),R*sin(v)],[u,v],'sphere')
           sage: sphere.mean_curvature()
           -1/R

        """
        gg = self.first_fundamental_form_coefficients()
        hh = self.second_fundamental_form_coefficients()
        denom = 2*self.area_form_squared()
        numer = (gg[(2,2)]*hh[(1,1)] - 2*gg[(1,2)]*hh[(1,2)] +
                 gg[(1,1)]*hh[(2,2)]).simplify_full()
        return _simplify_full_rad(numer/denom)


    @cached_method
    def shape_operator_coefficients(self):
        r"""
        Returns the components of the shape operator of the surface as a
        dictionary. See ``shape_operator`` for more information.

        OUTPUT:

         - Dictionary where the keys are two-tuples ``(i, j)``, with values the
           corresponding component of the shape operator.

        EXAMPLES::

           sage: R = var('R')
           sage: u, v = var('u,v', domain='real')
           sage: assume(cos(v)>0)
           sage: sphere = ParametrizedSurface3D([R*cos(u)*cos(v),R*sin(u)*cos(v),R*sin(v)],[u,v],'sphere')
           sage: sphere.shape_operator_coefficients()
           {(1, 2): 0, (1, 1): -1/R, (2, 1): 0, (2, 2): -1/R}

        """

        gi = self.first_fundamental_form_inverse_coefficients()
        hh = self.second_fundamental_form_coefficients()

        sh_op11 = _simplify_full_rad(gi[(1,1)]*hh[(1,1)] + gi[(1,2)]*hh[(1,2)])
        sh_op12 = _simplify_full_rad(gi[(1,1)]*hh[(2,1)] + gi[(1,2)]*hh[(2,2)])
        sh_op21 = _simplify_full_rad(gi[(2,1)]*hh[(1,1)] + gi[(2,2)]*hh[(1,2)])
        sh_op22 = _simplify_full_rad(gi[(2,1)]*hh[(2,1)] + gi[(2,2)]*hh[(2,2)])

        return {(1,1): sh_op11, (1,2): sh_op12, (2,1): sh_op21, (2,2): sh_op22}


    def shape_operator(self):
        r"""
        Returns the shape operator of the surface as a matrix.  The shape
        operator is defined as the derivative of the Gauss map, and is
        computed here in terms of the first and second fundamental form by
        means of the Weingarten equations.

        OUTPUT:

         - Matrix of the shape operator

        EXAMPLES::

           sage: R = var('R')
           sage: assume(R>0)
           sage: u, v = var('u,v', domain='real')
           sage: assume(cos(v)>0)
           sage: sphere = ParametrizedSurface3D([R*cos(u)*cos(v),R*sin(u)*cos(v),R*sin(v)],[u,v],'sphere')
           sage: S = sphere.shape_operator(); S
           [-1/R    0]
           [   0 -1/R]

        The eigenvalues of the shape operator are the principal curvatures of
        the surface::

            sage: u, v = var('u,v', domain='real')
            sage: paraboloid = ParametrizedSurface3D([u, v, u^2+v^2], [u, v], 'paraboloid')
            sage: S = paraboloid.shape_operator(); S
            [2*(4*v^2 + 1)/(4*u^2 + 4*v^2 + 1)^(3/2)        -8*u*v/(4*u^2 + 4*v^2 + 1)^(3/2)]
            [       -8*u*v/(4*u^2 + 4*v^2 + 1)^(3/2) 2*(4*u^2 + 1)/(4*u^2 + 4*v^2 + 1)^(3/2)]
            sage: S.eigenvalues()
            [2*sqrt(4*u^2 + 4*v^2 + 1)/(16*u^4 + 16*v^4 + 8*(4*u^2 + 1)*v^2 + 8*u^2 + 1), 2/sqrt(4*u^2 + 4*v^2 + 1)]

        """

        shop = self.shape_operator_coefficients()
        shop_matrix=matrix([[shop[(1,1)],shop[(1,2)]],
                            [shop[(2,1)],shop[(2,2)]]])
        return shop_matrix


    def principal_directions(self):
        r"""
        Finds the principal curvatures and principal directions of the surface

        OUTPUT:

        For each principal curvature, returns a list of the form
        $(\rho, V, n)$, where $\rho$ is the principal curvature,
        $V$ is the corresponding principal direction, and $n$ is
        the multiplicity.

        EXAMPLES::

           sage: u, v = var('u, v', domain='real')
           sage: R, r = var('R,r', domain='real')
           sage: assume(R>r,r>0)
           sage: torus = ParametrizedSurface3D([(R+r*cos(v))*cos(u),(R+r*cos(v))*sin(u),r*sin(v)],[u,v],'torus')
           sage: torus.principal_directions()
           [(-cos(v)/(r*cos(v) + R), [(1, 0)], 1), (-1/r, [(0, 1)], 1)]

        """
        return self.shape_operator().eigenvectors_left()


    @cached_method
    def connection_coefficients(self):
        r"""
        Computes the connection coefficients or Christoffel symbols
        $\Gamma^k_{ij}$ of the surface. If the coefficients of the first
        fundamental form are given by $g_{ij}$ (where $i, j = 1, 2$), then
        $\Gamma^k_{ij} = \frac{1}{2} g^{kl} \left( \frac{\partial g_{li}}{\partial x^j}
        - \frac{\partial g_{ij}}{\partial x^l}
        + \frac{\partial g_{lj}}{\partial x^i} \right)$.
        Here, $(g^{kl})$ is the inverse of the matrix $(g_{ij})$, with
        $i, j = 1, 2$.

        OUTPUT:

        Dictionary of connection coefficients, where the keys are 3-tuples
        $(i,j,k)$ and the values are the corresponding coefficients
        $\Gamma^k_{ij}$.

        EXAMPLES::

           sage: r = var('r')
           sage: assume(r > 0)
           sage: u, v = var('u,v', domain='real')
           sage: assume(cos(v)>0)
           sage: sphere = ParametrizedSurface3D([r*cos(u)*cos(v),r*sin(u)*cos(v),r*sin(v)],[u,v],'sphere')
           sage: sphere.connection_coefficients()
           {(1, 2, 1): -sin(v)/cos(v), (2, 2, 2): 0, (1, 2, 2): 0, (2, 1, 1): -sin(v)/cos(v), (1, 1, 2): cos(v)*sin(v), (2, 2, 1): 0, (2, 1, 2): 0, (1, 1, 1): 0}

        """
        x = self.variables
        gg = self.first_fundamental_form_coefficients()
        gi = self.first_fundamental_form_inverse_coefficients()

        dg = {}
        for i,j,k in product((1, 2), repeat=3):
            dg[(i,j,k)] = _simplify_full_rad(gg[(j,k)].differentiate(x[i]))

        structfun={}
        for i,j,k in product((1, 2), repeat=3):
            structfun[(i,j,k)] = sum(gi[(k,s)]*(dg[(i,j,s)] + dg[(j,i,s)]
                                                -dg[(s,i,j)])/2
                                     for s in (1,2))
            structfun[(i,j,k)] = _simplify_full_rad(structfun[(i,j,k)])
        return structfun


    @cached_method
    def _create_geodesic_ode_system(self):
        r"""
        Helper method to create a fast floating-point version of the
        geodesic equations, used by :meth:`geodesics_numerical`.

        EXAMPLES::
           sage: p, q = var('p,q', domain='real')
           sage: sphere = ParametrizedSurface3D([cos(q)*cos(p),sin(q)*cos(p),sin(p)],[p,q],'sphere')
           sage: ode = sphere._create_geodesic_ode_system()
           sage: ode.function(0.0, (1.0, 0.0, 1.0, 1.0))
           [1.00000000000000, 1.00000000000000, -0.4546487134128409, 3.114815449309804]

        """
        from sage.ext.fast_eval import fast_float
        from sage.calculus.var import var
        from sage.gsl.ode import ode_solver

        u1 = self.variables[1]
        u2 = self.variables[2]
        v1, v2 = var('v1, v2', domain='real')

        C = self.connection_coefficients()

        dv1 = - C[(1,1,1)]*v1**2 - 2*C[(1,2,1)]*v1*v2 - C[(2,2,1)]*v2**2
        dv2 = - C[(1,1,2)]*v1**2 - 2*C[(1,2,2)]*v1*v2 - C[(2,2,2)]*v2**2
        fun1 = fast_float(dv1, str(u1), str(u2), str(v1), str(v2))
        fun2 = fast_float(dv2, str(u1), str(u2), str(v1), str(v2))

        geodesic_ode = ode_solver()
        geodesic_ode.function = \
                              lambda t, (u1, u2, v1, v2) : \
                              [v1, v2, fun1(u1, u2, v1, v2), fun2(u1, u2, v1, v2)]
        return geodesic_ode


    def geodesics_numerical(self, p0, v0, tinterval):
        r"""
        Numerical integration of the geodesic equations.  Explicitly, the
        geodesic equations are given by
        $\frac{d^2 u^i}{dt^2} + \Gamma^i_{jk} \frac{d u^j}{dt} \frac{d u^k}{dt} = 0$.

        Solving these equations gives the coordinates $(u^1, u^2)$ of
        the geodesic on the surface.  The coordinates in space can
        then be found by substituting $(u^1, u^2)$ into the vector
        $\vec{r}(u^1, u^2)$ representing the surface.

        ALGORITHM:

        The geodesic equations are integrated forward in time using
        the ode solvers from ``sage.gsl.ode``.  See the member
        function ``_create_geodesic_ode_system`` for more details.

        INPUT:

         - ``p0`` - 2-tuple with coordinates of the initial point.

         - ``v0`` - 2-tuple with components of the initial tangent vector to the geodesic.

         - ``tinterval`` - List ``[a, b, M]``, where ``(a,b)`` is the domain of the geodesic and ``M`` is the number of subdivision points used when returning the solution.

        OUTPUT:

        List of lists ``[t, [u1(t), u2(t)], [v1(t), v2(t)], [x1(t), x2(t), x3(t)]]``, where

         - ``t`` is a subdivision point;

         - ``[u1(t), u2(t)]`` are the intrinsic coordinates of the geodesic point;

         - ``[v1(t), v2(t)]`` are the intrinsic coordinates of the tangent vector to the geodesic;

         - ``[x1(t), x2(t), x3(t)]`` are the coordinates of the geodesic point in the three-dimensional space.

        EXAMPLES::

           sage: p, q = var('p,q', domain='real')
           sage: assume(cos(q)>0)
           sage: sphere = ParametrizedSurface3D([cos(q)*cos(p),sin(q)*cos(p),sin(p)],[p,q],'sphere')
           sage: geodesic = sphere.geodesics_numerical([0.0,0.0],[1.0,1.0],[0,2*pi,5])
           sage: times, points, tangent_vectors, ext_points = zip(*geodesic)

           sage: round4 = lambda vec: [N(x, digits=4) for x in vec]   # helper function to round to 4 digits
           sage: round4(times)
           [0.0000, 1.257, 2.513, 3.770, 5.027, 6.283]
           sage: [round4(p) for p in points]
           [[0.0000, 0.0000], [0.7644, 1.859], [-0.2876, 3.442], [-0.6137, 5.502], [0.5464, 6.937], [0.3714, 9.025]]
           sage: [round4(p) for p in ext_points]
           [[1.000, 0.0000, 0.0000], [-0.2049, 0.6921, 0.6921], [-0.9160, -0.2836, -0.2836], [0.5803, -0.5759, -0.5759], [0.6782, 0.5196, 0.5196], [-0.8582, 0.3629, 0.3629]]
        """

        u1 = self.variables[1]
        u2 = self.variables[2]

        solver = self._create_geodesic_ode_system()

        t_interval, n = tinterval[0:2], tinterval[2]
        solver.y_0 = [p0[0], p0[1], v0[0], v0[1]]
        solver.ode_solve(t_span=t_interval, num_points=n)

        parsed_solution = \
          [[vec[0], vec[1][0:2], vec[1][2:], self.point(vec[1])]
           for vec in solver.solution]

        return parsed_solution


    @cached_method
    def _create_pt_ode_system(self, curve, t):
        """
        Helper method to create a fast floating-point version of the parallel
        transport equations, used by ``parallel_translation_numerical``.

        INPUT:

         - ``curve`` - curve in intrinsic coordinates along which to do parallel transport.
         - ``t`` - curve parameter

        EXAMPLES::
           sage: p, q = var('p,q', domain='real')
           sage: sphere = ParametrizedSurface3D([cos(q)*cos(p),sin(q)*cos(p),sin(p)],[p,q],'sphere')
           sage: s = var('s')
           sage: ode = sphere._create_pt_ode_system((s, s), s)
           sage: ode.function(0.0, (1.0, 1.0))
           [-0.0, 0.0]

        """

        from sage.ext.fast_eval import fast_float
        from sage.calculus.var import var
        from sage.gsl.ode import ode_solver

        u1 = self.variables[1]
        u2 = self.variables[2]
        v1, v2 = var('v1, v2', domain='real')

        du1 = diff(curve[0], t)
        du2 = diff(curve[1], t)

        C = self.connection_coefficients()
        for coef in C:
            C[coef] = C[coef].subs({u1: curve[0], u2: curve[1]})

        dv1 = - C[(1,1,1)]*v1*du1 - C[(1,2,1)]*(du1*v2 + du2*v1) - \
            C[(2,2,1)]*du2*v2
        dv2 = - C[(1,1,2)]*v1*du1 - C[(1,2,2)]*(du1*v2 + du2*v1) - \
            C[(2,2,2)]*du2*v2
        fun1 = fast_float(dv1, str(t), str(v1), str(v2))
        fun2 = fast_float(dv2, str(t), str(v1), str(v2))

        pt_ode = ode_solver()
        pt_ode.function = lambda t, (v1, v2): [fun1(t, v1, v2), fun2(t, v1, v2)]
        return pt_ode


    def parallel_translation_numerical(self,curve,t,v0,tinterval):
        r"""
        Numerically solves the equations for parallel translation of a vector
        along a curve on the surface.  Explicitly, the equations for parallel
        translation are given by
        $\frac{d u^i}{dt} + u^j \frac{d c^k}{dt} \Gamma^i_{jk} = 0$,
        where $\Gamma^i_{jk}$ are the connection coefficients of the surface,
        the vector to be transported has components $u^j$ and the curve along
        which to transport has components $c^k$.

        ALGORITHM:

        The parallel transport equations are integrated forward in time using
        the ode solvers from ``sage.gsl.ode``. See :meth:`_create_pt_ode_system`
        for more details.

        INPUT:

         - ``curve`` - 2-tuple of functions which determine the curve with respect to
           the local coordinate system;

         - ``t`` - symbolic variable denoting the curve parameter;

         - ``v0`` - 2-tuple representing the initial vector;

         - ``tinterval`` - list ``[a, b, N]``, where ``(a, b)`` is the domain of the curve
           and ``N`` is the number of subdivision points.

        OUTPUT:

        The list consisting of lists ``[t, [v1(t), v2(t)]]``, where

         - ``t`` is a subdivision point;

         - ``[v1(t), v2(t)]`` is the list of coordinates of the vector parallel translated
           along the curve.

        EXAMPLES::

           sage: p, q = var('p,q', domain='real')
           sage: v = [p,q]
           sage: assume(cos(q)>0)
           sage: sphere = ParametrizedSurface3D([cos(q)*cos(p),sin(q)*cos(p),sin(p)],v,'sphere')
           sage: s = var('s')
           sage: vector_field = sphere.parallel_translation_numerical([s,s],s,[1.0,1.0],[0.0, pi/4, 5])
           sage: times, components = zip(*vector_field)

           sage: round4 = lambda vec: [N(x, digits=4) for x in vec]   # helper function to round to 4 digits
           sage: round4(times)
           [0.0000, 0.1571, 0.3142, 0.4712, 0.6283, 0.7854]
           sage: [round4(v) for v in components]
           [[1.000, 1.000], [0.9876, 1.025], [0.9499, 1.102], [0.8853, 1.238], [0.7920, 1.448], [0.6687, 1.762]]

        """

        u1 = self.variables[1]
        u2 = self.variables[2]

        solver = self._create_pt_ode_system(tuple(curve), t)

        t_interval, n = tinterval[0:2], tinterval[2]
        solver.y_0 = v0
        solver.ode_solve(t_span=t_interval, num_points=n)

        return solver.solution
