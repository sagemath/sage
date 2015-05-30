r"""
Plotting Functions

EXAMPLES::

    sage: def f(x,y):
    ...       return math.sin(y*y+x*x)/math.sqrt(x*x+y*y+.0001)
    ...
    sage: P = plot3d(f,(-3,3),(-3,3), adaptive=True, color=rainbow(60, 'rgbtuple'), max_bend=.1, max_depth=15)
    sage: P.show()

::

    sage: def f(x,y):
    ...       return math.exp(x/5)*math.sin(y)
    ...
    sage: P = plot3d(f,(-5,5),(-5,5), adaptive=True, color=['red','yellow'])
    sage: from sage.plot.plot3d.plot3d import axes
    sage: S = P + axes(6, color='black')
    sage: S.show()

We plot "cape man"::

    sage: S = sphere(size=.5, color='yellow')

::

    sage: from sage.plot.plot3d.shapes import Cone
    sage: S += Cone(.5, .5, color='red').translate(0,0,.3)

::

    sage: S += sphere((.45,-.1,.15), size=.1, color='white') + sphere((.51,-.1,.17), size=.05, color='black')
    sage: S += sphere((.45, .1,.15),size=.1, color='white') + sphere((.51, .1,.17), size=.05, color='black')
    sage: S += sphere((.5,0,-.2),size=.1, color='yellow')
    sage: def f(x,y): return math.exp(x/5)*math.cos(y)
    sage: P = plot3d(f,(-5,5),(-5,5), adaptive=True, color=['red','yellow'], max_depth=10)
    sage: cape_man = P.scale(.2) + S.translate(1,0,0)
    sage: cape_man.show(aspect_ratio=[1,1,1])

Or, we plot a very simple function indeed::

    sage: plot3d(pi, (-1,1), (-1,1))
    Graphics3d Object

AUTHORS:

- Tom Boothby: adaptive refinement triangles

- Josh Kantor: adaptive refinement triangles

- Robert Bradshaw (2007-08): initial version of this file

- William Stein (2007-12, 2008-01): improving 3d plotting

- Oscar Lazo, William Cauchois, Jason Grout (2009-2010): Adding coordinate transformations
"""


#TODO:
#    -- smooth triangles

#*****************************************************************************
#      Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
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


from tri_plot import TrianglePlot
from index_face_set import IndexFaceSet
from shapes import arrow3d
from base import Graphics3dGroup
from sage.plot.colors import rainbow
from texture import Texture

from sage.ext.fast_eval import fast_float_arg

from sage.functions.trig import cos, sin

class _Coordinates(object):
    """
    This abstract class encapsulates a new coordinate system for plotting.
    Sub-classes must implement the :meth:`transform` method which, given
    symbolic variables to use, generates a 3-tuple of functions in terms of
    those variables that can be used to find the cartesian (X, Y, and Z)
    coordinates for any point in this space.
    """
    def __init__(self, dep_var, indep_vars):
        """
        INPUT:

         - ``dep_var`` - The dependent variable (the function value will be
           substituted for this).

         - ``indep_vars`` - A list of independent variables (the parameters will be
           substituted for these).

        TESTS:

        Because the base :class:`_Coordinates` class automatically checks the
        initializing variables with the transform method, :class:`_Coordinates`
        cannot be instantiated by itself.  We test a subclass.

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates as arb
            sage: x,y,z=var('x,y,z')
            sage: arb((x+z,y*z,z), z, (x,y))
            Arbitrary Coordinates coordinate transform (z in terms of x, y)
        """
        import inspect
        all_vars=inspect.getargspec(self.transform).args[1:]
        if set(all_vars) != set(indep_vars + [dep_var]):
            raise ValueError('variables were specified incorrectly for this coordinate system; incorrect variables were %s'%list(set(all_vars).symmetric_difference(set(indep_vars+[dep_var]))))
        self.dep_var = dep_var
        self.indep_vars = indep_vars

    @property
    def _name(self):
        """
        A default name for a coordinate system.  Override this in a
        subclass to set a different name.

        TESTS::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates as arb
            sage: x,y,z=var('x,y,z')
            sage: c=arb((x+z,y*z,z), z, (x,y))
            sage: c._name
            'Arbitrary Coordinates'
        """
        return self.__class__.__name__

    def transform(self, **kwds):
        """
        Return the transformation for this coordinate system in terms of the
        specified variables (which should be keywords).

        TESTS::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates as arb
            sage: x,y,z=var('x,y,z')
            sage: c=arb((x+z,y*z,z), z, (x,y))
            sage: c.transform(x=1,y=2,z=3)
            (4, 6, 3)
        """
        raise NotImplementedError

    def to_cartesian(self, func, params=None):
        """
        Returns a 3-tuple of functions, parameterized over ``params``, that
        represents the cartesian coordinates of the value of ``func``.

        INPUT:

         - ``func`` - A function in this coordinate space. Corresponds to the
           independent variable.

         - ``params`` - The parameters of func. Corresponds to the dependent
           variables.

        EXAMPLE::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates
            sage: x, y, z = var('x y z')
            sage: T = _ArbitraryCoordinates((x + y, x - y, z), z,[x,y])
            sage: f(x, y) = 2*x+y
            sage: T.to_cartesian(f, [x, y])
            (x + y, x - y, 2*x + y)
            sage: [h(1,2) for h in T.to_cartesian(lambda x,y: 2*x+y)]
            [3.0, -1.0, 4.0]

        We try to return a function having the same variable names as
        the function passed in::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates
            sage: x, y, z = var('x y z')
            sage: T = _ArbitraryCoordinates((x + y, x - y, z), z,[x,y])
            sage: f(a, b) = 2*a+b
            sage: T.to_cartesian(f, [a, b])
            (a + b, a - b, 2*a + b)
            sage: t1,t2,t3=T.to_cartesian(lambda a,b: 2*a+b)
            sage: import inspect
            sage: inspect.getargspec(t1)
            ArgSpec(args=['a', 'b'], varargs=None, keywords=None, defaults=None)
            sage: inspect.getargspec(t2)
            ArgSpec(args=['a', 'b'], varargs=None, keywords=None, defaults=None)
            sage: inspect.getargspec(t3)
            ArgSpec(args=['a', 'b'], varargs=None, keywords=None, defaults=None)
            sage: def g(a,b): return 2*a+b
            sage: t1,t2,t3=T.to_cartesian(g)
            sage: inspect.getargspec(t1)
            ArgSpec(args=['a', 'b'], varargs=None, keywords=None, defaults=None)
            sage: t1,t2,t3=T.to_cartesian(2*a+b)
            sage: inspect.getargspec(t1)
            ArgSpec(args=['a', 'b'], varargs=None, keywords=None, defaults=None)

        If we cannot guess the right parameter names, then the
        parameters are named `u` and `v`::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates
            sage: x, y, z = var('x y z')
            sage: T = _ArbitraryCoordinates((x + y, x - y, z), z,[x,y])
            sage: t1,t2,t3=T.to_cartesian(operator.add)
            sage: inspect.getargspec(t1)
            ArgSpec(args=['u', 'v'], varargs=None, keywords=None, defaults=None)
            sage: [h(1,2) for h in T.to_cartesian(operator.mul)]
            [3.0, -1.0, 2.0]
            sage: [h(u=1,v=2) for h in T.to_cartesian(operator.mul)]
            [3.0, -1.0, 2.0]

        The output of the function `func` is coerced to a float when
        it is evaluated if the function is something like a lambda or
        python callable. This takes care of situations like f returning a
        singleton numpy array, for example.

            sage: from numpy import array
            sage: v_phi=array([ 0.,  1.57079637,  3.14159274, 4.71238911,  6.28318548])
            sage: v_theta=array([ 0.,  0.78539819,  1.57079637,  2.35619456,  3.14159274])
            sage: m_r=array([[ 0.16763356,  0.25683223,  0.16649297,  0.10594339, 0.55282422],
            ... [ 0.16763356,  0.19993708,  0.31403568,  0.47359696, 0.55282422],
            ... [ 0.16763356,  0.25683223,  0.16649297,  0.10594339, 0.55282422],
            ... [ 0.16763356,  0.19993708,  0.31403568,  0.47359696, 0.55282422],
            ... [ 0.16763356,  0.25683223,  0.16649297,  0.10594339, 0.55282422]])
            sage: import scipy.interpolate
            sage: f=scipy.interpolate.RectBivariateSpline(v_phi,v_theta,m_r)
            sage: spherical_plot3d(f,(0,2*pi),(0,pi))
            Graphics3d Object

        """
        from sage.symbolic.expression import is_Expression
        from sage.rings.real_mpfr import is_RealNumber
        from sage.rings.integer import is_Integer
        if params is not None and (is_Expression(func) or is_RealNumber(func) or is_Integer(func)):
            return self.transform(**{
                self.dep_var: func,
                self.indep_vars[0]: params[0],
                self.indep_vars[1]: params[1]
            })
        else:
            # func might be a lambda or a Python callable; this makes it slightly
            # more complex.
            import sage.symbolic.ring
            dep_var_dummy = sage.symbolic.ring.var(self.dep_var)
            indep_var_dummies = sage.symbolic.ring.var(','.join(self.indep_vars))
            transformation = self.transform(**{
                self.dep_var: dep_var_dummy,
                self.indep_vars[0]: indep_var_dummies[0],
                self.indep_vars[1]: indep_var_dummies[1]
            })
            if params is None:
                if callable(func):
                    params = _find_arguments_for_callable(func)
                    if params is None:
                        params=['u','v']
                else:
                    raise ValueError("function is not callable")
            def subs_func(t):
                # We use eval so that the lambda function has the same
                # variable names as the original function
                ll="""lambda {x},{y}: t.subs({{
                    dep_var_dummy: float(func({x}, {y})),
                    indep_var_dummies[0]: float({x}),
                    indep_var_dummies[1]: float({y})
                }})""".format(x=params[0], y=params[1])
                return eval(ll,dict(t=t, func=func, dep_var_dummy=dep_var_dummy,
                                    indep_var_dummies=indep_var_dummies))
            return map(subs_func, transformation)

    def __repr__(self):
        """
        Print out a coordinate system

        ::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates as arb
            sage: x,y,z=var('x,y,z')
            sage: c=arb((x+z,y*z,z), z, (x,y))
            sage: c
            Arbitrary Coordinates coordinate transform (z in terms of x, y)
            sage: c.__dict__['_name'] = 'My Special Coordinates'
            sage: c
            My Special Coordinates coordinate transform (z in terms of x, y)
        """
        return '%s coordinate transform (%s in terms of %s)' % \
          (self._name, self.dep_var, ', '.join(self.indep_vars))


import inspect

def _find_arguments_for_callable(func):
    """
    Find the names of arguments (that do not have default values) for
    a callable function, taking care of several special cases in Sage.
    If the parameters cannot be found, then return None.

    EXAMPLES::

        sage: from sage.plot.plot3d.plot3d import _find_arguments_for_callable
        sage: _find_arguments_for_callable(lambda x,y: x+y)
        ['x', 'y']
        sage: def f(a,b,c): return a+b+c
        sage: _find_arguments_for_callable(f)
        ['a', 'b', 'c']
        sage: _find_arguments_for_callable(lambda x,y,z=2: x+y+z)
        ['x', 'y']
        sage: def f(a,b,c,d=2,e=1): return a+b+c+d+e
        sage: _find_arguments_for_callable(f)
        ['a', 'b', 'c']
        sage: g(w,r,t)=w+r+t
        sage: _find_arguments_for_callable(g)
        ['w', 'r', 't']
        sage: a,b = var('a,b')
        sage: _find_arguments_for_callable(a+b)
        ['a', 'b']
        sage: _find_arguments_for_callable(operator.add)
    """
    if inspect.isfunction(func):
        f_args=inspect.getargspec(func)
        if f_args.defaults is None:
            params=f_args.args
        else:
            params=f_args.args[:-len(f_args.defaults)]
    else:
        try:
            f_args=inspect.getargspec(func.__call__)
            if f_args.defaults is None:
                params=f_args.args
            else:
                params=f_args.args[:-len(f_args.defaults)]
        except TypeError:
            # func.__call__ may be a built-in (or Cython) function
            if hasattr(func, 'arguments'):
                params=[repr(s) for s in func.arguments()]
            else:
                params=None
    return params


class _ArbitraryCoordinates(_Coordinates):
    """
    An arbitrary coordinate system.
    """
    _name = "Arbitrary Coordinates"

    def __init__(self, custom_trans, dep_var, indep_vars):
        """
        Initialize an arbitrary coordinate system.

        INPUT:

         - ``custom_trans`` - A 3-tuple of transformation
           functions.

         - ``dep_var`` - The dependent (function) variable.

         - ``indep_vars`` - a list of the two other independent
           variables.

        EXAMPLES::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates
            sage: x, y, z = var('x y z')
            sage: T = _ArbitraryCoordinates((x + y, x - y, z), z,[x,y])
            sage: f(x, y) = 2*x + y
            sage: T.to_cartesian(f, [x, y])
            (x + y, x - y, 2*x + y)
            sage: [h(1,2) for h in T.to_cartesian(lambda x,y: 2*x+y)]
            [3.0, -1.0, 4.0]
        """
        self.dep_var = str(dep_var)
        self.indep_vars = [str(i) for i in indep_vars]
        self.custom_trans = tuple(custom_trans)

    def transform(self, **kwds):
        """
        EXAMPLE::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates
            sage: x, y, z = var('x y z')
            sage: T = _ArbitraryCoordinates((x + y, x - y, z), x,[y,z])

            sage: T.transform(x=z,y=1)
            (z + 1, z - 1, z)
        """
        return tuple(t.subs(**kwds) for t in self.custom_trans)

class Spherical(_Coordinates):
    """
    A spherical coordinate system for use with ``plot3d(transformation=...)``
    where the position of a point is specified by three numbers:

    - the *radial distance* (``radius``) from the origin

    - the *azimuth angle* (``azimuth``) from the positive `x`-axis

    - the *inclination angle* (``inclination``) from the positive `z`-axis

    These three variables must be specified in the constructor.

    EXAMPLES:

    Construct a spherical transformation for a function for the radius
    in terms of the azimuth and inclination::

        sage: T = Spherical('radius', ['azimuth', 'inclination'])

    If we construct some concrete variables, we can get a
    transformation in terms of those variables::

        sage: r, phi, theta = var('r phi theta')
        sage: T.transform(radius=r, azimuth=theta, inclination=phi)
        (r*cos(theta)*sin(phi), r*sin(phi)*sin(theta), r*cos(phi))

    We can plot with this transform.  Remember that the dependent
    variable is the radius, and the independent variables are the
    azimuth and the inclination (in that order)::

        sage: plot3d(phi * theta, (theta, 0, pi), (phi, 0, 1), transformation=T)
        Graphics3d Object

    We next graph the function where the inclination angle is constant::

        sage: S=Spherical('inclination', ['radius', 'azimuth'])
        sage: r,theta=var('r,theta')
        sage: plot3d(3, (r,0,3), (theta, 0, 2*pi), transformation=S)
        Graphics3d Object

    See also :func:`spherical_plot3d` for more examples of plotting in spherical
    coordinates.
    """

    def transform(self, radius=None, azimuth=None, inclination=None):
        """
        A spherical coordinates transform.

        EXAMPLE::

            sage: T = Spherical('radius', ['azimuth', 'inclination'])
            sage: T.transform(radius=var('r'), azimuth=var('theta'), inclination=var('phi'))
            (r*cos(theta)*sin(phi), r*sin(phi)*sin(theta), r*cos(phi))
        """
        return (radius * sin(inclination) * cos(azimuth),
                radius * sin(inclination) * sin(azimuth),
                radius * cos(inclination))

class SphericalElevation(_Coordinates):
    """
    A spherical coordinate system for use with ``plot3d(transformation=...)``
    where the position of a point is specified by three numbers:

    - the *radial distance* (``radius``) from the origin

    - the *azimuth angle* (``azimuth``) from the positive `x`-axis

    - the *elevation angle* (``elevation``) from the `xy`-plane toward the
      positive `z`-axis

    These three variables must be specified in the constructor.

    EXAMPLES:

    Construct a spherical transformation for the radius
    in terms of the azimuth and elevation. Then, get a
    transformation in terms of those variables::

        sage: T = SphericalElevation('radius', ['azimuth', 'elevation'])
        sage: r, theta, phi = var('r theta phi')
        sage: T.transform(radius=r, azimuth=theta, elevation=phi)
        (r*cos(phi)*cos(theta), r*cos(phi)*sin(theta), r*sin(phi))

    We can plot with this transform.  Remember that the dependent
    variable is the radius, and the independent variables are the
    azimuth and the elevation (in that order)::

        sage: plot3d(phi * theta, (theta, 0, pi), (phi, 0, 1), transformation=T)
        Graphics3d Object

    We next graph the function where the elevation angle is constant. This
    should be compared to the similar example for the ``Spherical`` coordinate
    system::

        sage: SE=SphericalElevation('elevation', ['radius', 'azimuth'])
        sage: r,theta=var('r,theta')
        sage: plot3d(3, (r,0,3), (theta, 0, 2*pi), transformation=SE)
        Graphics3d Object

    Plot a sin curve wrapped around the equator::

        sage: P1=plot3d( (pi/12)*sin(8*theta), (r,0.99,1), (theta, 0, 2*pi), transformation=SE, plot_points=(10,200))
        sage: P2=sphere(center=(0,0,0), size=1, color='red', opacity=0.3)
        sage: P1+P2
        Graphics3d Object

    Now we graph several constant elevation functions alongside several constant
    inclination functions. This example illustrates the difference between the
    ``Spherical`` coordinate system and the ``SphericalElevation`` coordinate
    system::

        sage: r, phi, theta = var('r phi theta')
        sage: SE = SphericalElevation('elevation', ['radius', 'azimuth'])
        sage: angles = [pi/18, pi/12, pi/6]
        sage: P1 = [plot3d( a, (r,0,3), (theta, 0, 2*pi), transformation=SE, opacity=0.85, color='blue') for a in angles]

        sage: S = Spherical('inclination', ['radius', 'azimuth'])
        sage: P2 = [plot3d( a, (r,0,3), (theta, 0, 2*pi), transformation=S, opacity=0.85, color='red') for a in angles]
        sage: show(sum(P1+P2), aspect_ratio=1)

    See also :func:`spherical_plot3d` for more examples of plotting in spherical
    coordinates.
    """

    def transform(self, radius=None, azimuth=None, elevation=None):
        """
        A spherical elevation coordinates transform.

        EXAMPLE::

            sage: T = SphericalElevation('radius', ['azimuth', 'elevation'])
            sage: T.transform(radius=var('r'), azimuth=var('theta'), elevation=var('phi'))
            (r*cos(phi)*cos(theta), r*cos(phi)*sin(theta), r*sin(phi))
        """
        return (radius * cos(elevation) * cos(azimuth),
                radius * cos(elevation) * sin(azimuth),
                radius * sin(elevation))

class Cylindrical(_Coordinates):
    """
    A cylindrical coordinate system for use with ``plot3d(transformation=...)``
    where the position of a point is specified by three numbers:

    - the *radial distance* (``radius``) from the `z`-axis

    - the *azimuth angle* (``azimuth``) from the positive `x`-axis

    - the *height* or *altitude* (``height``) above the `xy`-plane

    These three variables must be specified in the constructor.

    EXAMPLES:

    Construct a cylindrical transformation for a function for ``height`` in terms of
    ``radius`` and ``azimuth``::

        sage: T = Cylindrical('height', ['radius', 'azimuth'])

    If we construct some concrete variables, we can get a transformation::

        sage: r, theta, z = var('r theta z')
        sage: T.transform(radius=r, azimuth=theta, height=z)
        (r*cos(theta), r*sin(theta), z)

    We can plot with this transform.  Remember that the dependent
    variable is the height, and the independent variables are the
    radius and the azimuth (in that order)::

        sage: plot3d(9-r^2, (r, 0, 3), (theta, 0, pi), transformation=T)
        Graphics3d Object

    We next graph the function where the radius is constant::

        sage: S=Cylindrical('radius', ['azimuth', 'height'])
        sage: theta,z=var('theta, z')
        sage: plot3d(3, (theta,0,2*pi), (z, -2, 2), transformation=S)
        Graphics3d Object

    See also :func:`cylindrical_plot3d` for more examples of plotting in cylindrical
    coordinates.
    """

    def transform(self, radius=None, azimuth=None, height=None):
        """
        A cylindrical coordinates transform.

        EXAMPLE::

            sage: T = Cylindrical('height', ['azimuth', 'radius'])
            sage: T.transform(radius=var('r'), azimuth=var('theta'), height=var('z'))
            (r*cos(theta), r*sin(theta), z)
        """
        return (radius * cos(azimuth),
                radius * sin(azimuth),
                height)

class TrivialTriangleFactory:
    """
    Class emulating behavior of :class:`~sage.plot.plot3d.tri_plot.TriangleFactory`
    but simply returning a list of vertices for both regular and
    smooth triangles.
    """
    def triangle(self, a, b, c, color = None):
        """
        Function emulating behavior of
        :meth:`~sage.plot.plot3d.tri_plot.TriangleFactory.triangle`
        but simply returning a list of vertices.

        INPUT:

        - ``a``, ``b``, ``c`` : triples (x,y,z) representing corners
          on a triangle in 3-space
        - ``color``: ignored

        OUTPUT:

        - the list ``[a,b,c]``

        TESTS::

            sage: from sage.plot.plot3d.plot3d import TrivialTriangleFactory
            sage: factory = TrivialTriangleFactory()
            sage: tri = factory.triangle([0,0,0],[0,0,1],[1,1,0])
            sage: tri
            [[0, 0, 0], [0, 0, 1], [1, 1, 0]]
        """
        return [a,b,c]
    def smooth_triangle(self, a, b, c, da, db, dc, color = None):
        """
        Function emulating behavior of
        :meth:`~sage.plot.plot3d.tri_plot.TriangleFactory.smooth_triangle`
        but simply returning a list of vertices.

        INPUT:

        - ``a``, ``b``, ``c`` : triples (x,y,z) representing corners
          on a triangle in 3-space
        - ``da``, ``db``, ``dc`` : ignored
        - ``color`` : ignored

        OUTPUT:

        - the list ``[a,b,c]``

        TESTS::

            sage: from sage.plot.plot3d.plot3d import TrivialTriangleFactory
            sage: factory = TrivialTriangleFactory()
            sage: sm_tri = factory.smooth_triangle([0,0,0],[0,0,1],[1,1,0],[0,0,1],[0,2,0],[1,0,0])
            sage: sm_tri
            [[0, 0, 0], [0, 0, 1], [1, 1, 0]]
        """
        return [a,b,c]

import parametric_plot3d
def plot3d(f, urange, vrange, adaptive=False, transformation=None, **kwds):
    """
    INPUT:


    -  ``f`` - a symbolic expression or function of 2
       variables

    -  ``urange`` - a 2-tuple (u_min, u_max) or a 3-tuple
       (u, u_min, u_max)

    -  ``vrange`` - a 2-tuple (v_min, v_max) or a 3-tuple
       (v, v_min, v_max)

    -  ``adaptive`` - (default: False) whether to use
       adaptive refinement to draw the plot (slower, but may look better).
       This option does NOT work in conjuction with a transformation
       (see below).

    -  ``mesh`` - bool (default: False) whether to display
       mesh grid lines

    -  ``dots`` - bool (default: False) whether to display
       dots at mesh grid points

    -  ``plot_points`` - (default: "automatic") initial number of sample
       points in each direction; an integer or a pair of integers


    - ``transformation`` - (default: None) a transformation to
      apply. May be a 3 or 4-tuple (x_func, y_func, z_func,
      independent_vars) where the first 3 items indicate a
      transformation to cartesian coordinates (from your coordinate
      system) in terms of u, v, and the function variable fvar (for
      which the value of f will be substituted). If a 3-tuple is
      specified, the independent variables are chosen from the range
      variables.  If a 4-tuple is specified, the 4th element is a list
      of independent variables.  ``transformation`` may also be a
      predefined coordinate system transformation like Spherical or
      Cylindrical.

    .. note::

       ``mesh`` and ``dots`` are not supported when using the Tachyon
       raytracer renderer.

    EXAMPLES: We plot a 3d function defined as a Python function::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2))
        Graphics3d Object

    We plot the same 3d function but using adaptive refinement::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2), adaptive=True)
        Graphics3d Object

    Adaptive refinement but with more points::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2), adaptive=True, initial_depth=5)
        Graphics3d Object

    We plot some 3d symbolic functions::

        sage: var('x,y')
        (x, y)
        sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
        Graphics3d Object
        sage: plot3d(sin(x*y), (x, -pi, pi), (y, -pi, pi))
        Graphics3d Object

    We give a plot with extra sample points::

        sage: var('x,y')
        (x, y)
        sage: plot3d(sin(x^2+y^2),(x,-5,5),(y,-5,5), plot_points=200)
        Graphics3d Object
        sage: plot3d(sin(x^2+y^2),(x,-5,5),(y,-5,5), plot_points=[10,100])
        Graphics3d Object

    A 3d plot with a mesh::

        sage: var('x,y')
        (x, y)
        sage: plot3d(sin(x-y)*y*cos(x),(x,-3,3),(y,-3,3), mesh=True)
        Graphics3d Object

    Two wobby translucent planes::

        sage: x,y = var('x,y')
        sage: P = plot3d(x+y+sin(x*y), (x,-10,10),(y,-10,10), opacity=0.87, color='blue')
        sage: Q = plot3d(x-2*y-cos(x*y),(x,-10,10),(y,-10,10),opacity=0.3,color='red')
        sage: P + Q
        Graphics3d Object

    We draw two parametric surfaces and a transparent plane::

        sage: L = plot3d(lambda x,y: 0, (-5,5), (-5,5), color="lightblue", opacity=0.8)
        sage: P = plot3d(lambda x,y: 4 - x^3 - y^2, (-2,2), (-2,2), color='green')
        sage: Q = plot3d(lambda x,y: x^3 + y^2 - 4, (-2,2), (-2,2), color='orange')
        sage: L + P + Q
        Graphics3d Object

    We draw the "Sinus" function (water ripple-like surface)::

        sage: x, y = var('x y')
        sage: plot3d(sin(pi*(x^2+y^2))/2,(x,-1,1),(y,-1,1))
        Graphics3d Object

    Hill and valley (flat surface with a bump and a dent)::

        sage: x, y = var('x y')
        sage: plot3d( 4*x*exp(-x^2-y^2), (x,-2,2), (y,-2,2))
        Graphics3d Object

    An example of a transformation::

        sage: r, phi, z = var('r phi z')
        sage: trans=(r*cos(phi),r*sin(phi),z)
        sage: plot3d(cos(r),(r,0,17*pi/2),(phi,0,2*pi),transformation=trans,opacity=0.87).show(aspect_ratio=(1,1,2),frame=False)

    An example of a transformation with symbolic vector::

        sage: cylindrical(r,theta,z)=[r*cos(theta),r*sin(theta),z]
        sage: plot3d(3,(theta,0,pi/2),(z,0,pi/2),transformation=cylindrical)
        Graphics3d Object

    Many more examples of transformations::

        sage: u, v, w = var('u v w')
        sage: rectangular=(u,v,w)
        sage: spherical=(w*cos(u)*sin(v),w*sin(u)*sin(v),w*cos(v))
        sage: cylindric_radial=(w*cos(u),w*sin(u),v)
        sage: cylindric_axial=(v*cos(u),v*sin(u),w)
        sage: parabolic_cylindrical=(w*v,(v^2-w^2)/2,u)

    Plot a constant function of each of these to get an idea of what it does::

        sage: A = plot3d(2,(u,-pi,pi),(v,0,pi),transformation=rectangular,plot_points=[100,100])
        sage: B = plot3d(2,(u,-pi,pi),(v,0,pi),transformation=spherical,plot_points=[100,100])
        sage: C = plot3d(2,(u,-pi,pi),(v,0,pi),transformation=cylindric_radial,plot_points=[100,100])
        sage: D = plot3d(2,(u,-pi,pi),(v,0,pi),transformation=cylindric_axial,plot_points=[100,100])
        sage: E = plot3d(2,(u,-pi,pi),(v,-pi,pi),transformation=parabolic_cylindrical,plot_points=[100,100])
        sage: @interact
        ... def _(which_plot=[A,B,C,D,E]):
        ...       show(which_plot)
        <html>...

    Now plot a function::

        sage: g=3+sin(4*u)/2+cos(4*v)/2
        sage: F = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=rectangular,plot_points=[100,100])
        sage: G = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=spherical,plot_points=[100,100])
        sage: H = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=cylindric_radial,plot_points=[100,100])
        sage: I = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=cylindric_axial,plot_points=[100,100])
        sage: J = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=parabolic_cylindrical,plot_points=[100,100])
        sage: @interact
        ... def _(which_plot=[F, G, H, I, J]):
        ...       show(which_plot)
        <html>...

    TESTS:

    Make sure the transformation plots work::

        sage: show(A + B + C + D + E)
        sage: show(F + G + H + I + J)

    Listing the same plot variable twice gives an error::

        sage: x, y = var('x y')
        sage: plot3d( 4*x*exp(-x^2-y^2), (x,-2,2), (x,-2,2))
        Traceback (most recent call last):
        ...
        ValueError: range variables should be distinct, but there are duplicates
    """
    if transformation is not None:
        params=None
        from sage.symbolic.callable import is_CallableSymbolicExpression
        # First, determine the parameters for f (from the first item of urange
        # and vrange, preferably).
        if len(urange) == 3 and len(vrange) == 3:
            params = (urange[0], vrange[0])
        elif is_CallableSymbolicExpression(f):
            params = f.variables()

        from sage.modules.vector_callable_symbolic_dense import Vector_callable_symbolic_dense
        if isinstance(transformation, (tuple, list,Vector_callable_symbolic_dense)):
            if len(transformation)==3:
                if params is None:
                    raise ValueError("must specify independent variable names in the ranges when using generic transformation")
                indep_vars = params
            elif len(transformation)==4:
                indep_vars = transformation[3]
                transformation = transformation[0:3]
            else:
                raise ValueError("unknown transformation type")
            # find out which variable is the function variable by
            # eliminating the parameter variables.
            all_vars = set(sum([list(s.variables()) for s in transformation],[]))
            dep_var=all_vars - set(indep_vars)
            if len(dep_var)==1:
                dep_var = dep_var.pop()
                transformation = _ArbitraryCoordinates(transformation, dep_var, indep_vars)
            else:
                raise ValueError("unable to determine the function variable in the transform")

        if isinstance(transformation, _Coordinates):
            R = transformation.to_cartesian(f, params)
            return parametric_plot3d.parametric_plot3d(R, urange, vrange, **kwds)
        else:
            raise ValueError('unknown transformation type')
    elif adaptive:
        P = plot3d_adaptive(f, urange, vrange, **kwds)
    else:
        u=fast_float_arg(0)
        v=fast_float_arg(1)
        P=parametric_plot3d.parametric_plot3d((u,v,f), urange, vrange, **kwds)
    P.frame_aspect_ratio([1.0,1.0,0.5])
    return P

def plot3d_adaptive(f, x_range, y_range, color="automatic",
                    grad_f=None,
                    max_bend=.5, max_depth=5, initial_depth=4, num_colors=128, **kwds):
    r"""
    Adaptive 3d plotting of a function of two variables.

    This is used internally by the plot3d command when the option
    ``adaptive=True`` is given.

    INPUT:


    -  ``f`` - a symbolic function or a Python function of
       3 variables.

    -  ``x_range`` - x range of values: 2-tuple (xmin,
       xmax) or 3-tuple (x,xmin,xmax)

    -  ``y_range`` - y range of values: 2-tuple (ymin,
       ymax) or 3-tuple (y,ymin,ymax)

    -  ``grad_f`` - gradient of f as a Python function

    -  ``color`` - "automatic" - a rainbow of num_colors
       colors

    -  ``num_colors`` - (default: 128) number of colors to
       use with default color

    -  ``max_bend`` - (default: 0.5)

    -  ``max_depth`` - (default: 5)

    -  ``initial_depth`` - (default: 4)

    -  ``**kwds`` - standard graphics parameters


    EXAMPLES:

    We plot `\sin(xy)`::

        sage: from sage.plot.plot3d.plot3d import plot3d_adaptive
        sage: x,y=var('x,y'); plot3d_adaptive(sin(x*y), (x,-pi,pi), (y,-pi,pi), initial_depth=5)
        Graphics3d Object
    """
    if initial_depth >= max_depth:
        max_depth = initial_depth

    from sage.plot.misc import setup_for_eval_on_grid
    g, ranges = setup_for_eval_on_grid(f, [x_range,y_range], plot_points=2)
    xmin,xmax = ranges[0][:2]
    ymin,ymax = ranges[1][:2]

    opacity = kwds.get('opacity',1)

    if color == "automatic":
        texture = rainbow(num_colors, 'rgbtuple')
    else:
        if isinstance(color, list):
            texture = color
        else:
            kwds['color'] = color
            texture = Texture(kwds)

    factory = TrivialTriangleFactory()
    plot = TrianglePlot(factory, g, (xmin, xmax), (ymin, ymax), g = grad_f,
                        min_depth=initial_depth, max_depth=max_depth,
                        max_bend=max_bend, num_colors = None)

    P = IndexFaceSet(plot._objects)
    if isinstance(texture, (list, tuple)):
        if len(texture) == 2:
            # do a grid coloring
            xticks = (xmax - xmin)/2**initial_depth
            yticks = (ymax - ymin)/2**initial_depth
            parts = P.partition(lambda x,y,z: (int((x-xmin)/xticks) + int((y-ymin)/yticks)) % 2)
        else:
            # do a topo coloring
            bounds = P.bounding_box()
            min_z = bounds[0][2]
            max_z = bounds[1][2]
            if max_z == min_z:
                span = 0
            else:
                span = (len(texture)-1) / (max_z - min_z)    # max to avoid dividing by 0
            parts = P.partition(lambda x,y,z: int((z-min_z)*span))
        all = []
        for k, G in parts.iteritems():
            G.set_texture(texture[k], opacity=opacity)
            all.append(G)
        P = Graphics3dGroup(all)
    else:
        P.set_texture(texture)

    P.frame_aspect_ratio([1.0,1.0,0.5])
    P._set_extra_kwds(kwds)
    return P

def spherical_plot3d(f, urange, vrange, **kwds):
    """
    Plots a function in spherical coordinates.  This function is
    equivalent to::

        sage: r,u,v=var('r,u,v')
        sage: f=u*v; urange=(u,0,pi); vrange=(v,0,pi)
        sage: T = (r*cos(u)*sin(v), r*sin(u)*sin(v), r*cos(v), [u,v])
        sage: plot3d(f, urange, vrange, transformation=T)
        Graphics3d Object

    or equivalently::

        sage: T = Spherical('radius', ['azimuth', 'inclination'])
        sage: f=lambda u,v: u*v; urange=(u,0,pi); vrange=(v,0,pi)
        sage: plot3d(f, urange, vrange, transformation=T)
        Graphics3d Object

    INPUT:

    - ``f`` - a symbolic expression or function of two variables.

    - ``urange`` - a 3-tuple (u, u_min, u_max), the domain of the azimuth variable.

    - ``vrange`` - a 3-tuple (v, v_min, v_max), the domain of the inclination variable.

    EXAMPLES:

    A sphere of radius 2::

        sage: x,y=var('x,y')
        sage: spherical_plot3d(2,(x,0,2*pi),(y,0,pi))
        Graphics3d Object

    The real and imaginary parts of a spherical harmonic with `l=2` and `m=1`::

        sage: phi, theta = var('phi, theta')
        sage: Y = spherical_harmonic(2, 1, theta, phi)
        sage: rea = spherical_plot3d(abs(real(Y)), (phi,0,2*pi), (theta,0,pi), color='blue', opacity=0.6)
        sage: ima = spherical_plot3d(abs(imag(Y)), (phi,0,2*pi), (theta,0,pi), color='red', opacity=0.6)
        sage: (rea + ima).show(aspect_ratio=1)  # long time (4s on sage.math, 2011)

    A drop of water::

        sage: x,y=var('x,y')
        sage: spherical_plot3d(e^-y,(x,0,2*pi),(y,0,pi),opacity=0.5).show(frame=False)

    An object similar to a heart::

        sage: x,y=var('x,y')
        sage: spherical_plot3d((2+cos(2*x))*(y+1),(x,0,2*pi),(y,0,pi),rgbcolor=(1,.1,.1))
        Graphics3d Object

    Some random figures:

    ::

        sage: x,y=var('x,y')
        sage: spherical_plot3d(1+sin(5*x)/5,(x,0,2*pi),(y,0,pi),rgbcolor=(1,0.5,0),plot_points=(80,80),opacity=0.7)
        Graphics3d Object

    ::

        sage: x,y=var('x,y')
        sage: spherical_plot3d(1+2*cos(2*y),(x,0,3*pi/2),(y,0,pi)).show(aspect_ratio=(1,1,1))
    """
    return plot3d(f, urange, vrange, transformation=Spherical('radius', ['azimuth', 'inclination']), **kwds)

def cylindrical_plot3d(f, urange, vrange, **kwds):
    """
    Plots a function in cylindrical coordinates.  This function is
    equivalent to::

        sage: r,u,v=var('r,u,v')
        sage: f=u*v; urange=(u,0,pi); vrange=(v,0,pi)
        sage: T = (r*cos(u), r*sin(u), v, [u,v])
        sage: plot3d(f, urange, vrange, transformation=T)
        Graphics3d Object

    or equivalently::

        sage: T = Cylindrical('radius', ['azimuth', 'height'])
        sage: f=lambda u,v: u*v; urange=(u,0,pi); vrange=(v,0,pi)
        sage: plot3d(f, urange, vrange, transformation=T)
        Graphics3d Object


    INPUT:

    - ``f`` - a symbolic expression or function of two variables,
      representing the radius from the `z`-axis.

    - ``urange`` - a 3-tuple (u, u_min, u_max), the domain of the
      azimuth variable.

    - ``vrange`` - a 3-tuple (v, v_min, v_max), the domain of the
      elevation (`z`) variable.

    EXAMPLES:

    A portion of a cylinder of radius 2::

        sage: theta,z=var('theta,z')
        sage: cylindrical_plot3d(2,(theta,0,3*pi/2),(z,-2,2))
        Graphics3d Object

    Some random figures:

    ::

        sage: cylindrical_plot3d(cosh(z),(theta,0,2*pi),(z,-2,2))
        Graphics3d Object

    ::

        sage: cylindrical_plot3d(e^(-z^2)*(cos(4*theta)+2)+1,(theta,0,2*pi),(z,-2,2),plot_points=[80,80]).show(aspect_ratio=(1,1,1))
    """
    return plot3d(f, urange, vrange, transformation=Cylindrical('radius', ['azimuth', 'height']), **kwds)

def axes(scale=1, radius=None, **kwds):
    """
    Creates basic axes in three dimensions.  Each axis is a three
    dimensional arrow object.

    INPUT:

    - ``scale`` - (default: 1) The length of the axes (all three
      will be the same).
    - ``radius`` - (default: .01) The radius of the axes as arrows.

    EXAMPLES::

        sage: from sage.plot.plot3d.plot3d import axes
        sage: S = axes(6, color='black'); S
        Graphics3d Object

    ::

        sage: T = axes(2, .5); T
        Graphics3d Object
    """
    if radius is None:
        radius = scale/100.0
    return Graphics3dGroup([arrow3d((0,0,0),(scale,0,0), radius, **kwds),
                            arrow3d((0,0,0),(0,scale,0), radius, **kwds),
                            arrow3d((0,0,0),(0,0,scale), radius, **kwds)])
