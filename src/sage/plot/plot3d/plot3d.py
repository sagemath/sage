# sage.doctest: needs sage.symbolic
r"""
Plotting functions

EXAMPLES::

    sage: x, y = var('x y')
    sage: W = plot3d(sin(pi*((x)^2 + (y)^2))/2, (x, -1, 1), (y, -1, 1),
    ....:            frame=False, color='purple', opacity=0.8)
    sage: S = sphere((0, 0, 0), size=0.3, color='red', aspect_ratio=[1,1,1])
    sage: show(W + S, figsize=8)

.. PLOT::

    x, y = var('x y')
    W = plot3d(sin(pi*((x)**2 + (y)**2))/2, (x, -1, 1), (y, -1, 1), frame=False, color='purple', opacity=0.8)
    S = sphere((0, 0, 0), size=0.3, color='red', aspect_ratio=[1, 1, 1])
    sphinx_plot(W + S)

::

    sage: def f(x, y):
    ....:     return math.sin(y^2 + x^2)/math.sqrt(x^2 + y^2 + 0.0001)
    sage: P = plot3d(f, (-3, 3),(-3, 3), adaptive=True,
    ....:            color=rainbow(60, 'rgbtuple'), max_bend=.1, max_depth=15)
    sage: P.show()

.. PLOT::

    def f(x, y): return math.sin(y**2 + x**2)/math.sqrt(x**2 + y**2 + 0.0001)
    P = plot3d(f, (-3, 3), (-3, 3), adaptive=True, color=rainbow(60, 'rgbtuple'), max_bend=.1, max_depth=15)
    sphinx_plot(P)

::

    sage: def f(x, y):
    ....:     return math.exp(x/5)*math.sin(y)
    ....:
    sage: P = plot3d(f, (-5, 5), (-5, 5), adaptive=True, color=['red', 'yellow'])
    sage: from sage.plot.plot3d.plot3d import axes
    sage: S = P + axes(6, color='black')
    sage: S.show()

.. PLOT::

    def f(x, y): return math.exp(x/5)*math.sin(y)
    P = plot3d(f, (-5, 5), (-5, 5), adaptive=True, color=['red', 'yellow'])
    from sage.plot.plot3d.plot3d import axes
    S = P + axes(6, color='black')
    sphinx_plot(S)

Here is an example using a colormap and a color function ``c``::

    sage: x, y = var('x y')
    sage: cm = colormaps.hsv
    sage: def c(x, y): return float((x + y + x*y)/15) % 1
    sage: plot3d(x*x + y*y, (x, -4, 4), (y, -4, 4), color=(c, cm))
    Graphics3d Object

.. PLOT::

    x, y = var('x y')
    cm = colormaps.hsv
    def c(x, y): return float((x + y + x*y)/15) % 1
    sphinx_plot(plot3d(x*x + y*y,(x, -4, 4), (y, -4, 4), color=(c, cm)))

Beware that the color function must take values between 0 and 1.

We plot "cape man"::

    sage: S = sphere(size=.5, color='yellow')

::

    sage: from sage.plot.plot3d.shapes import Cone
    sage: S += Cone(.5, .5, color='red').translate(0,0,.3)

::

    sage: S += sphere((.45, -.1, .15), size=.1, color='white')
    sage: S += sphere((.51,-.1,.17), size=.05, color='black')
    sage: S += sphere((.45, .1, .15), size=.1, color='white')
    sage: S += sphere((.51, .1,.17), size=.05, color='black')
    sage: S += sphere((.5, 0, -.2), size=.1, color='yellow')
    sage: def f(x, y): return math.exp(x/5)*math.cos(y)
    sage: P = plot3d(f, (-5, 5), (-5, 5), adaptive=True,
    ....:            color=['red','yellow'], max_depth=10)
    sage: cape_man = P.scale(.2) + S.translate(1, 0, 0)
    sage: cape_man.show(aspect_ratio=[1, 1, 1])

.. PLOT::

    S = sphere(size=.5, color='yellow')
    from sage.plot.plot3d.shapes import Cone
    S += Cone(.5, .5, color='red').translate(0,0,.3)
    S += sphere((.45, -.1, .15), size=.1, color='white') + sphere((.51,-.1,.17), size=.05, color='black')
    S += sphere((.45, .1, .15), size=.1, color='white') + sphere((.51, .1,.17), size=.05, color='black')
    S += sphere((.5, 0, -.2), size=.1, color='yellow')
    def f(x, y): return math.exp(x/5)*math.cos(y)
    P = plot3d(f, (-5, 5), (-5, 5), adaptive=True, color=['red','yellow'], max_depth=10)
    cape_man = P.scale(.2) + S.translate(1, 0, 0)
    cape_man.aspect_ratio([1, 1, 1])
    sphinx_plot(cape_man)

Or, we plot a very simple function indeed::

    sage: plot3d(pi, (-1,1), (-1,1))
    Graphics3d Object

.. PLOT::

    sphinx_plot(plot3d(pi, (-1,1), (-1,1)))

Transparent with fractional opacity value::

    sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2), opacity=8/10)
    Graphics3d Object

.. TODO::

    Add support for smooth triangles.

AUTHORS:

- Tom Boothby: adaptive refinement triangles

- Josh Kantor: adaptive refinement triangles

- Robert Bradshaw (2007-08): initial version of this file

- William Stein (2007-12, 2008-01): improving 3d plotting

- Oscar Lazo, William Cauchois, Jason Grout (2009-2010): Adding coordinate transformations
"""
# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.misc.sageinspect import sage_getargspec, is_function_or_cython_function
from sage.plot.colors import rainbow
from sage.plot.plot3d.base import Graphics3dGroup
from sage.plot.plot3d.index_face_set import IndexFaceSet
from sage.plot.plot3d.shapes import arrow3d
from sage.plot.plot3d.texture import Texture
from sage.plot.plot3d.tri_plot import TrianglePlot
from . import parametric_plot3d
lazy_import("sage.functions.trig", ["cos", "sin"])


class _Coordinates:
    """
    This abstract class encapsulates a new coordinate system for plotting.
    Sub-classes must implement the :meth:`transform` method which, given
    symbolic variables to use, generates a 3-tuple of functions in terms of
    those variables that can be used to find the Cartesian (X, Y, and Z)
    coordinates for any point in this space.

    INPUT:

    - ``dep_var`` -- the dependent variable (the function value will be

    - ``indep_vars`` -- list of independent variables (the parameters will be
      substituted for these)
    """
    def __init__(self, dep_var, indep_vars):
        """
        Initialize.

        TESTS:

        Because the base :class:`_Coordinates` class automatically checks the
        initializing variables with the transform method, :class:`_Coordinates`
        cannot be instantiated by itself.  We test a subclass::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates as arb
            sage: x,y,z=var('x,y,z')
            sage: arb((x+z,y*z,z), z, (x,y))
            Arbitrary Coordinates coordinate transform (z in terms of x, y)
        """
        all_vars = sage_getargspec(self.transform).args[1:]
        A = set(all_vars)
        B = set(indep_vars + [dep_var])
        if A != B:
            raise ValueError('variables were specified incorrectly for this '
                             'coordinate system; incorrect variables '
                             'were %s' % list(A.symmetric_difference(B)))
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
        Return a 3-tuple of functions, parameterized over ``params``, that
        represents the Cartesian coordinates of the value of ``func``.

        INPUT:

        - ``func`` -- function in this coordinate space; corresponds
          to the independent variable

        - ``params`` -- the parameters of ``func``; corresponds to
          the dependent variables

        EXAMPLES::

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
            sage: from sage.misc.sageinspect import sage_getargspec
            sage: sage_getargspec(t1)
            FullArgSpec(args=['a', 'b'], varargs=None, varkw=None, defaults=None,
                        kwonlyargs=[], kwonlydefaults=None, annotations={})
            sage: sage_getargspec(t2)
            FullArgSpec(args=['a', 'b'], varargs=None, varkw=None, defaults=None,
                        kwonlyargs=[], kwonlydefaults=None, annotations={})
            sage: sage_getargspec(t3)
            FullArgSpec(args=['a', 'b'], varargs=None, varkw=None, defaults=None,
                        kwonlyargs=[], kwonlydefaults=None, annotations={})

            sage: def g(a, b): return 2*a+b
            sage: t1,t2,t3=T.to_cartesian(g)
            sage: sage_getargspec(t1)
            FullArgSpec(args=['a', 'b'], varargs=None, varkw=None, defaults=None,
                        kwonlyargs=[], kwonlydefaults=None, annotations={})
            sage: t1,t2,t3=T.to_cartesian(2*a+b)
            sage: sage_getargspec(t1)
            FullArgSpec(args=['a', 'b'], varargs=None, varkw=None, defaults=None,
                        kwonlyargs=[], kwonlydefaults=None, annotations={})

        If we cannot guess the right parameter names, then the
        parameters are named `u` and `v`::

            sage: from sage.plot.plot3d.plot3d import _ArbitraryCoordinates
            sage: from sage.misc.sageinspect import sage_getargspec
            sage: x, y, z = var('x y z')
            sage: T = _ArbitraryCoordinates((x + y, x - y, z), z,[x,y])
            sage: t1,t2,t3=T.to_cartesian(operator.add)
            sage: sage_getargspec(t1)
            FullArgSpec(args=['u', 'v'], varargs=None, varkw=None, defaults=None,
                        kwonlyargs=[], kwonlydefaults=None, annotations={})
            sage: [h(1,2) for h in T.to_cartesian(operator.mul)]
            [3.0, -1.0, 2.0]
            sage: [h(u=1,v=2) for h in T.to_cartesian(operator.mul)]
            [3.0, -1.0, 2.0]

        The output of the function ``func`` is coerced to a float when
        it is evaluated if the function is something like a lambda or
        python callable. This takes care of situations like f returning a
        singleton numpy array, for example.

            sage: from numpy import array
            sage: v_phi=array([ 0.,  1.57079637,  3.14159274, 4.71238911,  6.28318548])
            sage: v_theta=array([ 0.,  0.78539819,  1.57079637,  2.35619456,  3.14159274])
            sage: m_r=array([[ 0.16763356,  0.25683223,  0.16649297,  0.10594339, 0.55282422],
            ....: [ 0.16763356,  0.19993708,  0.31403568,  0.47359696, 0.55282422],
            ....: [ 0.16763356,  0.25683223,  0.16649297,  0.10594339, 0.55282422],
            ....: [ 0.16763356,  0.19993708,  0.31403568,  0.47359696, 0.55282422],
            ....: [ 0.16763356,  0.25683223,  0.16649297,  0.10594339, 0.55282422]])
            sage: import scipy.interpolate
            sage: f=scipy.interpolate.RectBivariateSpline(v_phi,v_theta,m_r).ev
            sage: spherical_plot3d(f,(0,2*pi),(0,pi))
            Graphics3d Object
        """
        from sage.structure.element import Expression
        from sage.rings.real_mpfr import RealNumber
        from sage.rings.integer import Integer
        if params is not None and isinstance(func, (Expression,
                                                    RealNumber,
                                                    Integer)):
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
                    if not params:
                        params = ['u', 'v']
                else:
                    raise ValueError("function is not callable")

            def subs_func(t):
                # We use eval so that the lambda function has the same
                # variable names as the original function
                ll = f"""lambda {params[0]},{params[1]}: t.subs({{
                    dep_var_dummy: float(func({params[0]}, {params[1]})),
                    indep_var_dummies[0]: float({params[0]}),
                    indep_var_dummies[1]: float({params[1]})
                }})"""
                return eval(ll, {'t': t, 'func': func,
                                 'dep_var_dummy': dep_var_dummy,
                                 'indep_var_dummies': indep_var_dummies})
            return [subs_func(m) for m in transformation]

    def __repr__(self):
        """
        Print out a coordinate system.

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
        return '{} coordinate transform ({} in terms of {})'.format(self._name, self.dep_var, ', '.join(self.indep_vars))


def _find_arguments_for_callable(func):
    """
    Find the names of arguments (that do not have default values) for
    a callable function, taking care of several special cases in Sage.

    If the parameters cannot be found, then return ``[]``.

    EXAMPLES::

        sage: from sage.plot.plot3d.plot3d import _find_arguments_for_callable
        sage: _find_arguments_for_callable(lambda x,y: x+y)
        ['x', 'y']
        sage: def f(a, b, c): return a+b+c
        sage: _find_arguments_for_callable(f)
        ['a', 'b', 'c']
        sage: _find_arguments_for_callable(lambda x,y,z=2: x+y+z)
        ['x', 'y']
        sage: def f(a, b, c, d=2, e=1): return a+b+c+d+e
        sage: _find_arguments_for_callable(f)
        ['a', 'b', 'c']
        sage: g(w,r,t)=w+r+t
        sage: _find_arguments_for_callable(g)
        ['w', 'r', 't']
        sage: a,b = var('a,b')
        sage: _find_arguments_for_callable(a+b)
        ['a', 'b']
        sage: _find_arguments_for_callable(operator.add)
        []
    """
    if is_function_or_cython_function(func):
        pass
    elif hasattr(func, 'arguments'):
        # Might be a symbolic function with arguments
        return [repr(s) for s in func.arguments()]
    else:
        func = func.__call__

    f_args = sage_getargspec(func)
    if f_args.defaults is None:
        params = f_args.args
    else:
        params = f_args.args[:-len(f_args.defaults)]

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

        - ``custom_trans`` -- a 3-tuple of transformation
          functions

        - ``dep_var`` -- the dependent (function) variable

        - ``indep_vars`` -- list of the two other independent
          variables

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
        EXAMPLES::

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

    .. PLOT::

        r, phi, theta = var('r phi theta')
        T = Spherical('radius', ['azimuth', 'inclination'])
        sphinx_plot(plot3d(phi * theta, (theta, 0, pi), (phi, 0, 1), transformation=T))

    We next graph the function where the inclination angle is constant::

        sage: S = Spherical('inclination', ['radius', 'azimuth'])
        sage: r, theta = var('r,theta')
        sage: plot3d(3, (r,0,3), (theta, 0, 2*pi), transformation=S)
        Graphics3d Object

    .. PLOT::

        S = Spherical('inclination', ['radius', 'azimuth'])
        r, theta = var('r,theta')
        sphinx_plot(plot3d(r-r+3, (r,0,3), (theta, 0, 2*pi), transformation=S))

    See also :func:`spherical_plot3d` for more examples of plotting in spherical
    coordinates.
    """

    def transform(self, radius=None, azimuth=None, inclination=None):
        """
        A spherical coordinates transform.

        EXAMPLES::

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

    .. PLOT::

        T = SphericalElevation('radius', ['azimuth', 'elevation'])
        r, theta, phi = var('r theta phi')
        sphinx_plot(plot3d(phi * theta, (theta, 0, pi), (phi, 0, 1), transformation=T))

    We next graph the function where the elevation angle is constant. This
    should be compared to the similar example for the ``Spherical`` coordinate
    system::

        sage: SE = SphericalElevation('elevation', ['radius', 'azimuth'])
        sage: r, theta = var('r,theta')
        sage: plot3d(3, (r, 0, 3), (theta, 0, 2*pi), transformation=SE)
        Graphics3d Object

    .. PLOT::

        SE = SphericalElevation('elevation', ['radius', 'azimuth'])
        r, theta = var('r,theta')
        sphinx_plot(plot3d(3 + r - r, (r,0,3), (theta, 0, 2*pi), transformation=SE))

    Plot a sin curve wrapped around the equator::

        sage: P1 = plot3d((pi/12)*sin(8*theta), (r,0.99,1), (theta, 0, 2*pi),
        ....:             transformation=SE, plot_points=(10,200))
        sage: P2 = sphere(center=(0,0,0), size=1, color='red', opacity=0.3)
        sage: P1 + P2
        Graphics3d Object

    .. PLOT::

        r, theta = var('r,theta')
        SE = SphericalElevation('elevation', ['radius', 'azimuth'])
        P1 = plot3d( (pi/12)*sin(8*theta), (r, 0.99, 1), (theta, 0, 2*pi), transformation=SE, plot_points=(10,200))
        P2 = sphere(center=(0, 0, 0), size=1, color='red', opacity=0.3)
        sphinx_plot(P1 + P2)

    Now we graph several constant elevation functions alongside several constant
    inclination functions. This example illustrates the difference between the
    ``Spherical`` coordinate system and the ``SphericalElevation`` coordinate
    system::

        sage: r, phi, theta = var('r phi theta')
        sage: SE = SphericalElevation('elevation', ['radius', 'azimuth'])
        sage: angles = [pi/18, pi/12, pi/6]
        sage: P1 = [plot3d(a, (r,0,3), (theta, 0, 2*pi), transformation=SE,
        ....:              opacity=0.85, color='blue')
        ....:       for a in angles]

        sage: S = Spherical('inclination', ['radius', 'azimuth'])
        sage: P2 = [plot3d(a, (r,0,3), (theta, 0, 2*pi), transformation=S,
        ....:              opacity=0.85, color='red')
        ....:       for a in angles]
        sage: show(sum(P1+P2), aspect_ratio=1)

    .. PLOT::

        r, phi, theta = var('r phi theta')
        SE = SphericalElevation('elevation', ['radius', 'azimuth'])
        S = Spherical('inclination', ['radius', 'azimuth'])
        angles = [pi/18, pi/12, pi/6]
        P1=Graphics()
        P2=Graphics()
        for a in angles:
            P1 += plot3d( a, (r,0,3), (theta, 0, 2*pi), transformation=SE, opacity=0.85, color='blue')
            P2 += plot3d( a, (r,0,3), (theta, 0, 2*pi), transformation=S, opacity=0.85, color='red')
        sphinx_plot(P1+P2)

    See also :func:`spherical_plot3d` for more examples of plotting in spherical
    coordinates.
    """

    def transform(self, radius=None, azimuth=None, elevation=None):
        """
        A spherical elevation coordinates transform.

        EXAMPLES::

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

    .. PLOT::

        T = Cylindrical('height', ['radius', 'azimuth'])
        r, theta, z = var('r theta z')
        sphinx_plot(plot3d(9 - r**2, (r, 0, 3), (theta, 0, pi), transformation=T))

    We next graph the function where the radius is constant::

        sage: S = Cylindrical('radius', ['azimuth', 'height'])
        sage: theta, z = var('theta, z')
        sage: plot3d(3, (theta, 0, 2*pi), (z, -2, 2), transformation=S)
        Graphics3d Object

    .. PLOT::

        S = Cylindrical('radius', ['azimuth', 'height'])
        theta, z = var('theta, z')
        sphinx_plot(plot3d(3 + z - z, (theta, 0, 2*pi), (z, -2, 2), transformation=S))

    See also :func:`cylindrical_plot3d` for more examples of plotting in cylindrical
    coordinates.
    """

    def transform(self, radius=None, azimuth=None, height=None):
        """
        A cylindrical coordinates transform.

        EXAMPLES::

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
    def triangle(self, a, b, c, color=None):
        """
        Function emulating behavior of
        :meth:`~sage.plot.plot3d.tri_plot.TriangleFactory.triangle`
        but simply returning a list of vertices.

        INPUT:

        - ``a``, ``b``, ``c`` -- triples (x,y,z) representing corners
          on a triangle in 3-space

        - ``color`` -- ignored

        OUTPUT: the list ``[a,b,c]``

        TESTS::

            sage: from sage.plot.plot3d.plot3d import TrivialTriangleFactory
            sage: factory = TrivialTriangleFactory()
            sage: tri = factory.triangle([0,0,0],[0,0,1],[1,1,0])
            sage: tri
            [[0, 0, 0], [0, 0, 1], [1, 1, 0]]
        """
        return [a, b, c]

    def smooth_triangle(self, a, b, c, da, db, dc, color=None):
        """
        Function emulating behavior of
        :meth:`~sage.plot.plot3d.tri_plot.TriangleFactory.smooth_triangle`
        but simply returning a list of vertices.

        INPUT:

        - ``a``, ``b``, ``c`` -- triples (x,y,z) representing corners
          on a triangle in 3-space

        - ``da``, ``db``, ``dc`` -- ignored

        - ``color`` -- ignored

        OUTPUT: the list ``[a,b,c]``

        TESTS::

            sage: from sage.plot.plot3d.plot3d import TrivialTriangleFactory
            sage: factory = TrivialTriangleFactory()
            sage: sm_tri = factory.smooth_triangle([0,0,0],[0,0,1],[1,1,0],[0,0,1],[0,2,0],[1,0,0])
            sage: sm_tri
            [[0, 0, 0], [0, 0, 1], [1, 1, 0]]
        """
        return [a, b, c]


def plot3d(f, urange, vrange, adaptive=False, transformation=None, **kwds):
    """
    Plots a function in 3d.

    INPUT:

    - ``f`` -- a symbolic expression or function of 2
      variables

    - ``urange`` -- a 2-tuple (u_min, u_max) or a 3-tuple
      (u, u_min, u_max)

    - ``vrange`` -- a 2-tuple (v_min, v_max) or a 3-tuple
      (v, v_min, v_max)

    - ``adaptive`` -- boolean (default: ``False``); whether to use
      adaptive refinement to draw the plot (slower, but may look better).
      This option does NOT work in conjunction with a transformation
      (see below).

    - ``mesh`` -- boolean (default: ``False``); whether to display
      mesh grid lines

    - ``dots`` -- boolean (default: ``False``); whether to display
      dots at mesh grid points

    - ``plot_points`` -- (default: ``'automatic'``) initial number of sample
      points in each direction; an integer or a pair of integers

    - ``transformation`` -- (default: ``None``) a transformation to
      apply. May be a 3 or 4-tuple (x_func, y_func, z_func,
      independent_vars) where the first 3 items indicate a
      transformation to Cartesian coordinates (from your coordinate
      system) in terms of u, v, and the function variable fvar (for
      which the value of f will be substituted). If a 3-tuple is
      specified, the independent variables are chosen from the range
      variables.  If a 4-tuple is specified, the 4th element is a list
      of independent variables.  ``transformation`` may also be a
      predefined coordinate system transformation like Spherical or
      Cylindrical.

    .. NOTE::

       ``mesh`` and ``dots`` are not supported when using the Tachyon
       raytracer renderer.

    EXAMPLES: We plot a 3d function defined as a Python function::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2))
        Graphics3d Object

    .. PLOT::

        sphinx_plot(plot3d(lambda x, y: x**2 + y**2, (-2,2), (-2,2)))

    We plot the same 3d function but using adaptive refinement::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2), adaptive=True)
        Graphics3d Object

    .. PLOT::

        sphinx_plot(plot3d(lambda x, y: x**2 + y**2, (-2,2), (-2,2), adaptive=True))

    Adaptive refinement but with more points::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2), adaptive=True, initial_depth=5)
        Graphics3d Object

    .. PLOT::

        sphinx_plot(plot3d(lambda x, y: x**2 + y**2, (-2,2), (-2,2), adaptive=True, initial_depth=5))

    We plot some 3d symbolic functions::

        sage: var('x,y')
        (x, y)
        sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
        Graphics3d Object

    .. PLOT::

        var('x y')
        sphinx_plot(plot3d(x**2 + y**2, (x,-2,2), (y,-2,2)))

    ::

        sage: plot3d(sin(x*y), (x, -pi, pi), (y, -pi, pi))
        Graphics3d Object

    .. PLOT::

        var('x y')
        sphinx_plot(plot3d(sin(x*y), (x, -pi, pi), (y, -pi, pi)))

    We give a plot with extra sample points::

        sage: var('x,y')
        (x, y)
        sage: plot3d(sin(x^2 + y^2), (x,-5,5), (y,-5,5), plot_points=200)
        Graphics3d Object

    .. PLOT::

        var('x y')
        sphinx_plot(plot3d(sin(x**2 + y**2),(x,-5,5),(y,-5,5), plot_points=200))

    ::

        sage: plot3d(sin(x^2 + y^2), (x, -5, 5), (y, -5, 5), plot_points=[10, 100])
        Graphics3d Object

    .. PLOT::

        var('x y')
        sphinx_plot(plot3d(sin(x**2 + y**2), (x, -5, 5), (y, -5, 5), plot_points=[10,100]))

    A 3d plot with a mesh::

        sage: var('x,y')
        (x, y)
        sage: plot3d(sin(x - y)*y*cos(x), (x, -3, 3), (y, -3, 3), mesh=True)
        Graphics3d Object

    .. PLOT::

        var('x y')
        sphinx_plot(plot3d(sin(x - y)*y*cos(x), (x, -3, 3), (y, -3, 3), mesh=True))

    The same with thicker mesh lines (not supported in all viewers)::

        sage: var('x,y')
        (x, y)
        sage: plot3d(sin(x - y)*y*cos(x), (x, -3, 3), (y, -3, 3), mesh=True,
        ....:        thickness=2, viewer='threejs')
        Graphics3d Object

    Two wobby translucent planes::

        sage: x,y = var('x,y')
        sage: P = plot3d(x + y + sin(x*y), (x, -10, 10), (y, -10, 10),
        ....:            opacity=0.87, color='blue')
        sage: Q = plot3d(x - 2*y - cos(x*y),(x, -10, 10), (y, -10, 10),
        ....:            opacity=0.3, color='red')
        sage: P + Q
        Graphics3d Object

    .. PLOT::

        x,y = var('x,y')
        P = plot3d(x + y + sin(x*y), (x, -10, 10), (y, -10, 10), opacity=0.87, color='blue')
        Q = plot3d(x - 2*y - cos(x*y),(x, -10, 10), (y, -10, 10), opacity=0.3, color='red')
        sphinx_plot(P + Q)

    We draw two parametric surfaces and a transparent plane::

        sage: L = plot3d(lambda x,y: 0, (-5,5), (-5,5), color='lightblue', opacity=0.8)
        sage: P = plot3d(lambda x,y: 4 - x^3 - y^2, (-2,2), (-2,2), color='green')
        sage: Q = plot3d(lambda x,y: x^3 + y^2 - 4, (-2,2), (-2,2), color='orange')
        sage: L + P + Q
        Graphics3d Object

    .. PLOT::

        L = plot3d(lambda x,y: 0, (-5,5), (-5,5), color='lightblue', opacity=0.8)
        P = plot3d(lambda x,y: 4 - x**3 - y**2, (-2,2), (-2,2), color='green')
        Q = plot3d(lambda x,y: x**3 + y**2 - 4, (-2,2), (-2,2), color='orange')
        sphinx_plot(L + P + Q)

    We draw the "Sinus" function (water ripple-like surface)::

        sage: x, y = var('x y')
        sage: plot3d(sin(pi*(x^2 + y^2))/2, (x, -1, 1), (y, -1, 1))
        Graphics3d Object

    .. PLOT::

        x, y = var('x y')
        sphinx_plot(plot3d(sin(pi*(x**2 + y**2))/2, (x, -1, 1), (y, -1, 1)))

    Hill and valley (flat surface with a bump and a dent)::

        sage: x, y = var('x y')
        sage: plot3d(4*x*exp(-x^2 - y^2), (x, -2, 2), (y, -2, 2))
        Graphics3d Object

    .. PLOT::

        x, y = var('x y')
        sphinx_plot(plot3d( 4*x*exp(-x**2 - y**2), (x, -2, 2), (y, -2, 2)))

    An example of a transformation::

        sage: r, phi, z = var('r phi z')
        sage: trans = (r*cos(phi), r*sin(phi), z)
        sage: plot3d(cos(r), (r, 0, 17*pi/2), (phi, 0, 2*pi), transformation=trans, opacity=0.87).show(aspect_ratio=(1,1,2), frame=False)

    .. PLOT::

        r, phi, z = var('r phi z')
        trans = (r*cos(phi), r*sin(phi), z)
        P = plot3d(cos(r), (r, 0, 17*pi/2), (phi, 0, 2*pi), transformation=trans, opacity=0.87)
        P.aspect_ratio([1, 1, 2])
        sphinx_plot(P)

    An example of a transformation with symbolic vector::

        sage: cylindrical(r, theta, z) = [r*cos(theta), r*sin(theta), z]
        sage: plot3d(3, (theta, 0, pi/2), (z, 0, pi/2), transformation=cylindrical)
        Graphics3d Object

    .. PLOT::

        r, theta, z = var('r theta z')
        cylindrical = (r*cos(theta), r*sin(theta), z)
        P = plot3d(z-z+3, (theta, 0, pi/2), (z, 0, pi/2), transformation=cylindrical)
        sphinx_plot(P)

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
        ....: def _(which_plot=[A,B,C,D,E]):
        ....:     show(which_plot)
        ...Interactive function <function _ at ...> with 1 widget
          which_plot: Dropdown(description='which_plot', options=(Graphics3d Object, Graphics3d Object, Graphics3d Object, Graphics3d Object, Graphics3d Object), value=Graphics3d Object)

    Now plot a function::

        sage: g=3+sin(4*u)/2+cos(4*v)/2
        sage: F = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=rectangular,plot_points=[100,100])
        sage: G = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=spherical,plot_points=[100,100])
        sage: H = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=cylindric_radial,plot_points=[100,100])
        sage: I = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=cylindric_axial,plot_points=[100,100])
        sage: J = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=parabolic_cylindrical,plot_points=[100,100])
        sage: @interact
        ....: def _(which_plot=[F, G, H, I, J]):
        ....:     show(which_plot)
        ...Interactive function <function _ at ...> with 1 widget
          which_plot: Dropdown(description='which_plot', options=(Graphics3d Object, Graphics3d Object, Graphics3d Object, Graphics3d Object, Graphics3d Object), value=Graphics3d Object)

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

    Verify that :issue:`7423` is fixed::

        sage: f(x,y)=ln(x)
        sage: P=plot3d(f,(x,0,1),(y,0,1))
        sage: P
        Graphics3d Object
    """
    if transformation is not None:
        params = None
        from sage.structure.element import Expression
        # First, determine the parameters for f (from the first item of urange
        # and vrange, preferably).
        if len(urange) == 3 and len(vrange) == 3:
            params = (urange[0], vrange[0])
        elif isinstance(f, Expression) and f.is_callable():
            params = f.variables()

        from sage.modules.vector_callable_symbolic_dense import Vector_callable_symbolic_dense
        if isinstance(transformation, (tuple, list, Vector_callable_symbolic_dense)):
            if len(transformation) == 3:
                if params is None:
                    raise ValueError("must specify independent variable names in the ranges when using generic transformation")
                indep_vars = params
            elif len(transformation) == 4:
                indep_vars = transformation[3]
                transformation = transformation[0:3]
            else:
                raise ValueError("unknown transformation type")
            # find out which variable is the function variable by
            # eliminating the parameter variables.
            all_vars = set(sum([list(s.variables()) for s in transformation], []))
            dep_var = all_vars - set(indep_vars)
            if len(dep_var) == 1:
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
        arg1 = lambda u, v: u
        arg2 = lambda u, v: v
        P = parametric_plot3d.parametric_plot3d((arg1, arg2, f),
                                                urange,
                                                vrange,
                                                **kwds)
    P.frame_aspect_ratio([1.0, 1.0, 0.5])
    return P


def plot3d_adaptive(f, x_range, y_range, color='automatic',
                    grad_f=None,
                    max_bend=.5, max_depth=5, initial_depth=4, num_colors=128, **kwds):
    r"""
    Adaptive 3d plotting of a function of two variables.

    This is used internally by the plot3d command when the option
    ``adaptive=True`` is given.

    INPUT:

    - ``f`` -- a symbolic function or a Python function of 3 variables

    - ``x_range`` -- x range of values: 2-tuple (xmin,
      xmax) or 3-tuple (x,xmin,xmax)

    - ``y_range`` -- y range of values: 2-tuple (ymin, ymax) or 3-tuple
      (y,ymin,ymax)

    - ``grad_f`` -- gradient of f as a Python function

    - ``color`` -- "automatic"; a rainbow of num_colors colors

    - ``num_colors`` -- (default: 128) number of colors to use with default
      color

    - ``max_bend`` -- (default: 0.5)

    - ``max_depth`` -- (default: 5)

    - ``initial_depth`` -- (default: 4)

    - ``**kwds`` -- standard graphics parameters

    EXAMPLES:

    We plot `\sin(xy)`::

        sage: from sage.plot.plot3d.plot3d import plot3d_adaptive
        sage: x, y = var('x,y')
        sage: plot3d_adaptive(sin(x*y), (x, -pi, pi), (y, -pi, pi), initial_depth=5)
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.plot3d import plot3d_adaptive
        x, y = var('x,y')
        sphinx_plot(plot3d_adaptive(sin(x*y), (x,-pi,pi), (y,-pi,pi), initial_depth=5))
    """
    max_depth = max(max_depth, initial_depth)

    from sage.plot.misc import setup_for_eval_on_grid
    g, ranges = setup_for_eval_on_grid(f, [x_range, y_range], plot_points=2)
    xmin, xmax = ranges[0][:2]
    ymin, ymax = ranges[1][:2]

    opacity = float(kwds.get('opacity', 1))

    if color == "automatic":
        texture = rainbow(num_colors, 'rgbtuple')
    else:
        if isinstance(color, list):
            texture = color
        else:
            kwds['color'] = color
            texture = Texture(kwds)

    factory = TrivialTriangleFactory()
    plot = TrianglePlot(factory, g, (xmin, xmax), (ymin, ymax), g=grad_f,
                        min_depth=initial_depth, max_depth=max_depth,
                        max_bend=max_bend, num_colors=None)

    P = IndexFaceSet(plot._objects)
    if isinstance(texture, (list, tuple)):
        if len(texture) == 2:
            # do a grid coloring
            xticks = (xmax - xmin) / 2**initial_depth
            yticks = (ymax - ymin) / 2**initial_depth
            parts = P.partition(lambda x, y, z: (int((x - xmin) / xticks) + int((y - ymin) / yticks)) % 2)
        else:
            # do a topo coloring
            bounds = P.bounding_box()
            min_z = bounds[0][2]
            max_z = bounds[1][2]
            if max_z == min_z:
                span = 0
            else:
                span = (len(texture) - 1) / (max_z - min_z)  # max to avoid dividing by 0
            parts = P.partition(lambda x, y, z: int((z - min_z) * span))
        all = []
        for k, G in parts.items():
            G.set_texture(texture[k], opacity=opacity)
            all.append(G)
        P = Graphics3dGroup(all)
    else:
        P.set_texture(texture)

    P.frame_aspect_ratio([1.0, 1.0, 0.5])
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

    - ``f`` -- a symbolic expression or function of two variables

    - ``urange`` -- a 3-tuple (u, u_min, u_max), the domain of the azimuth variable

    - ``vrange`` -- a 3-tuple (v, v_min, v_max), the domain of the inclination variable

    EXAMPLES:

    A sphere of radius 2::

        sage: x,y = var('x,y')
        sage: spherical_plot3d(2, (x, 0, 2*pi), (y, 0, pi))
        Graphics3d Object

    .. PLOT::

        x, y = var('x,y')
        sphinx_plot(spherical_plot3d(x-x+2, (x, 0, 2*pi), (y, 0, pi)))

    The real and imaginary parts of a spherical harmonic with `l=2` and `m=1`::

        sage: phi, theta = var('phi, theta')
        sage: Y = spherical_harmonic(2, 1, theta, phi)
        sage: rea = spherical_plot3d(abs(real(Y)), (phi, 0, 2*pi), (theta, 0, pi), color='blue', opacity=0.6)
        sage: ima = spherical_plot3d(abs(imag(Y)), (phi, 0, 2*pi), (theta, 0, pi), color='red', opacity=0.6)
        sage: (rea + ima).show(aspect_ratio=1)  # long time (4s on sage.math, 2011)

    .. PLOT::

        phi, theta = var('phi, theta')
        Y = spherical_harmonic(2, 1, theta, phi)
        rea = spherical_plot3d(abs(real(Y)), (phi, 0, 2*pi), (theta, 0, pi), color='blue', opacity=0.6)
        ima = spherical_plot3d(abs(imag(Y)), (phi, 0, 2*pi), (theta, 0, pi), color='red', opacity=0.6)
        sphinx_plot(rea + ima)

    A drop of water::

        sage: x,y = var('x,y')
        sage: spherical_plot3d(e^-y, (x, 0, 2*pi), (y, 0, pi), opacity=0.5).show(frame=False)

    .. PLOT::

        x, y = var('x,y')
        sphinx_plot(spherical_plot3d(e**-y, (x, 0, 2*pi), (y, 0, pi), opacity=0.5))

    An object similar to a heart::

        sage: x,y = var('x,y')
        sage: spherical_plot3d((2 + cos(2*x))*(y + 1), (x, 0, 2*pi), (y, 0, pi), rgbcolor=(1, .1, .1))
        Graphics3d Object

    .. PLOT::

        x, y = var('x,y')
        sphinx_plot(spherical_plot3d((2 + cos(2*x))*(y + 1), (x, 0, 2*pi), (y, 0, pi), rgbcolor=(1, .1, .1)))

    Some random figures::

        sage: x,y = var('x,y')
        sage: spherical_plot3d(1 + sin(5*x)/5, (x, 0, 2*pi), (y, 0, pi), rgbcolor=(1, 0.5, 0), plot_points=(80, 80), opacity=0.7)
        Graphics3d Object

    .. PLOT::

        x,y = var('x,y')
        sphinx_plot(spherical_plot3d(1 + sin(5*x)/5, (x, 0, 2*pi), (y, 0, pi), rgbcolor=(1, 0.5, 0), plot_points=(80, 80), opacity=0.7))

    ::

        sage: x, y = var('x,y')
        sage: spherical_plot3d(1 + 2*cos(2*y), (x, 0, 3*pi/2), (y, 0, pi)).show(aspect_ratio=(1, 1, 1))

    .. PLOT::

        x, y = var('x,y')
        sphinx_plot(spherical_plot3d(1 + 2*cos(2*y), (x, 0, 3*pi/2), (y, 0, pi)))
    """
    return plot3d(f, urange, vrange, transformation=Spherical('radius', ['azimuth', 'inclination']), **kwds)


def cylindrical_plot3d(f, urange, vrange, **kwds):
    """
    Plots a function in cylindrical coordinates.  This function is
    equivalent to::

        sage: r, u, v = var('r,u,v')
        sage: f = u*v; urange = (u, 0, pi); vrange = (v, 0, pi)
        sage: T = (r*cos(u), r*sin(u), v, [u, v])
        sage: plot3d(f, urange, vrange, transformation=T)
        Graphics3d Object

    .. PLOT::

        r, u, v = var('r,u,v')
        f = u*v; urange = (u, 0, pi); vrange = (v, 0, pi)
        T = (r*cos(u), r*sin(u), v, [u,v])
        sphinx_plot(plot3d(f, urange, vrange, transformation=T))

    or equivalently::

        sage: T = Cylindrical('radius', ['azimuth', 'height'])
        sage: f=lambda u,v: u*v; urange=(u,0,pi); vrange=(v,0,pi)
        sage: plot3d(f, urange, vrange, transformation=T)
        Graphics3d Object

    INPUT:

    - ``f`` -- a symbolic expression or function of two variables,
      representing the radius from the `z`-axis

    - ``urange`` -- a 3-tuple (u, u_min, u_max), the domain of the
      azimuth variable

    - ``vrange`` -- a 3-tuple (v, v_min, v_max), the domain of the
      elevation (`z`) variable

    EXAMPLES:

    A portion of a cylinder of radius 2::

        sage: theta, z = var('theta,z')
        sage: cylindrical_plot3d(2, (theta, 0, 3*pi/2), (z, -2, 2))
        Graphics3d Object

    .. PLOT::

        theta, z = var('theta,z')
        sphinx_plot(cylindrical_plot3d(z-z+2, (theta, 0, 3*pi/2), (z, -2, 2)))

    Some random figures:

    ::

        sage: cylindrical_plot3d(cosh(z), (theta, 0, 2*pi), (z, -2, 2))
        Graphics3d Object

    .. PLOT::

        theta, z = var('theta,z')
        sphinx_plot(cylindrical_plot3d(cosh(z), (theta, 0, 2*pi), (z, -2, 2)))

    ::

        sage: cylindrical_plot3d(e^(-z^2)*(cos(4*theta) + 2) + 1, (theta, 0, 2*pi), (z, -2, 2), plot_points=[80, 80]).show(aspect_ratio=(1, 1, 1))

    .. PLOT::

        theta, z = var('theta,z')
        P = cylindrical_plot3d(e**(-z**2)*(cos(4*theta) + 2) + 1, (theta, 0, 2*pi), (z, -2, 2), plot_points=[80, 80])
        P.aspect_ratio([1, 1, 1])
        sphinx_plot(P)
    """
    return plot3d(f, urange, vrange, transformation=Cylindrical('radius', ['azimuth', 'height']), **kwds)


def axes(scale=1, radius=None, **kwds):
    """
    Create basic axes in three dimensions.  Each axis is a three
    dimensional arrow object.

    INPUT:

    - ``scale`` -- (default: 1) the length of the axes (all three will be the
      same)

    - ``radius`` -- (default: .01) the radius of the axes as arrows

    EXAMPLES::

        sage: from sage.plot.plot3d.plot3d import axes
        sage: S = axes(6, color='black'); S
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.plot3d import axes
        S = axes(6, color='black')
        sphinx_plot(S)

    ::

        sage: T = axes(2, .5); T
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.plot3d import axes
        T = axes(2, .5)
        sphinx_plot(T)
    """
    if radius is None:
        radius = scale / 100.0
    return Graphics3dGroup([arrow3d((0, 0, 0), (scale, 0, 0), radius, **kwds),
                            arrow3d((0, 0, 0), (0, scale, 0), radius, **kwds),
                            arrow3d((0, 0, 0), (0, 0, scale), radius, **kwds)])
