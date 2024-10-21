# sage_setup: distribution = sagemath-plot
# sage.doctest: needs sage.symbolic
"""
Plotting fields
"""
# ****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
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
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options
from sage.arith.srange import xsrange

# Below is the base class that is used to make 'field plots'.
# Its implementation is motivated by 'PlotField'.
# Currently it is used to make the functions 'plot_vector_field'
# and 'plot_slope_field'.
# TODO: use this to make these functions:
# 'plot_gradient_field' and 'plot_hamiltonian_field'


class PlotField(GraphicPrimitive):
    """
    Primitive class that initializes the
    PlotField graphics type
    """
    def __init__(self, xpos_array, ypos_array, xvec_array, yvec_array, options):
        """
        Create the graphics primitive PlotField.  This sets options
        and the array to be plotted as attributes.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: R=plot_slope_field(x + y, (x,0,1), (y,0,1), plot_points=2)
            sage: r=R[0]
            sage: r.options()['headaxislength']
            0
            sage: r.xpos_array
            [0.0, 0.0, 1.0, 1.0]
            sage: r.yvec_array
            masked_array(data=[0.0, 0.70710678118..., 0.70710678118...,
                               0.89442719...],
                         mask=[False, False, False, False],
                   fill_value=1e+20)

        TESTS:

        We test dumping and loading a plot::

            sage: x,y = var('x,y')
            sage: P = plot_vector_field((sin(x),cos(y)), (x,-3,3), (y,-3,3))
            sage: Q = loads(dumps(P))
        """
        self.xpos_array = xpos_array
        self.ypos_array = ypos_array
        self.xvec_array = xvec_array
        self.yvec_array = yvec_array
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Return a dictionary with the bounding box data.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: d = plot_vector_field((.01*x,x+y), (x,10,20), (y,10,20))[0].get_minmax_data()
            sage: d['xmin']
            10.0
            sage: d['ymin']
            10.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xpos_array, self.ypos_array, dict=True)

    def _allowed_options(self):
        """
        Return a dictionary with allowed options for PlotField.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: P=plot_vector_field((sin(x),cos(y)), (x,-3,3), (y,-3,3))
            sage: d=P[0]._allowed_options()
            sage: d['pivot']
            'Where the arrow should be placed in relation to the point (tail, middle, tip)'
        """
        return {'plot_points': 'How many points to use for plotting precision',
                'pivot': 'Where the arrow should be placed in relation to the point (tail, middle, tip)',
                'headwidth': 'Head width as multiple of shaft width, default is 3',
                'headlength': 'head length as multiple of shaft width, default is 5',
                'headaxislength': 'head length at shaft intersection, default is 4.5',
                'zorder': 'The layer level in which to draw',
                'color': 'The color of the arrows'}

    def _repr_(self):
        """
        String representation of PlotField graphics primitive.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: P=plot_vector_field((sin(x),cos(y)), (x,-3,3), (y,-3,3))
            sage: P[0]
            PlotField defined by a 20 x 20 vector grid

        TESTS:

        We check that :issue:`15052` is fixed
        (note that in general :issue:`15002` should be fixed)::

            sage: x,y=var('x,y')
            sage: P=plot_vector_field((sin(x), cos(y)), (x,-3,3), (y,-3,3), wrong_option='nonsense')
            sage: P[0].options()['plot_points']
            verbose 0 (...: primitive.py, options) WARNING: Ignoring option 'wrong_option'=nonsense
            verbose 0 (...: primitive.py, options)
            The allowed options for PlotField defined by a 20 x 20 vector grid are:
                color          The color of the arrows
                headaxislength head length at shaft intersection, default is 4.5
                headlength     head length as multiple of shaft width, default is 5
                headwidth      Head width as multiple of shaft width, default is 3
                pivot          Where the arrow should be placed in relation to the point (tail, middle, tip)
                plot_points    How many points to use for plotting precision
                zorder         The layer level in which to draw
            <BLANKLINE>
            20
        """
        return "PlotField defined by a {} x {} vector grid".format(
               self._options['plot_points'], self._options['plot_points'])

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: x,y = var('x,y')
            sage: P=plot_vector_field((sin(x),cos(y)), (x,-3,3), (y,-3,3))
        """
        options = self.options()
        quiver_options = options.copy()
        quiver_options.pop('plot_points')
        subplot.quiver(self.xpos_array, self.ypos_array,
                       self.xvec_array, self.yvec_array,
                       angles='xy', **quiver_options)


@options(plot_points=20, frame=True)
def plot_vector_field(f_g, xrange, yrange, **options):
    r"""
    ``plot_vector_field`` takes two functions of two variables xvar and yvar
    (for instance, if the variables are `x` and `y`, take `(f(x,y), g(x,y))`)
    and plots vector arrows of the function over the specified ranges, with
    xrange being of xvar between xmin and xmax, and yrange similarly
    (see below).

    ``plot_vector_field((f,g), (xvar,xmin,xmax), (yvar,ymin,ymax))``

    EXAMPLES:

    Plot some vector fields involving sin and cos::

        sage: x,y = var('x y')
        sage: plot_vector_field((sin(x),cos(y)), (x,-3,3), (y,-3,3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y = var('x y')
        g = plot_vector_field((sin(x),cos(y)), (x,-3,3), (y,-3,3))
        sphinx_plot(g)

    ::

        sage: plot_vector_field((y,(cos(x)-2) * sin(x)), (x,-pi,pi), (y,-pi,pi))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y = var('x y')
        g = plot_vector_field((y,(cos(x)-2) * sin(x)), (x,-pi,pi), (y,-pi,pi))
        sphinx_plot(g)

    Plot a gradient field::

        sage: u, v = var('u v')
        sage: f = exp(-(u^2 + v^2))
        sage: plot_vector_field(f.gradient(), (u,-2,2), (v,-2,2), color='blue')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        u, v = var('u v')
        f = exp(-(u**2 + v**2))
        g = plot_vector_field(f.gradient(), (u,-2,2), (v,-2,2), color='blue')
        sphinx_plot(g)

    Plot two orthogonal vector fields::

        sage: x,y = var('x,y')
        sage: a = plot_vector_field((x,y), (x,-3,3), (y,-3,3), color='blue')
        sage: b = plot_vector_field((y,-x), (x,-3,3), (y,-3,3), color='red')
        sage: show(a + b)

    .. PLOT::

        x,y = var('x,y')
        a = plot_vector_field((x,y), (x,-3,3), (y,-3,3), color='blue')
        b = plot_vector_field((y,-x), (x,-3,3), (y,-3,3), color='red')
        sphinx_plot(a + b)

    We ignore function values that are infinite or NaN::

        sage: x,y = var('x,y')
        sage: plot_vector_field((-x/sqrt(x^2+y^2),-y/sqrt(x^2+y^2)), (x,-10,10), (y,-10,10))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = plot_vector_field((-x/sqrt(x**2+y**2),-y/sqrt(x**2+y**2)), (x,-10,10), (y,-10,10))
        sphinx_plot(g)

    ::

        sage: x,y = var('x,y')
        sage: plot_vector_field((-x/sqrt(x+y),-y/sqrt(x+y)), (x,-10, 10), (y,-10,10))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = plot_vector_field((-x/sqrt(x+y),-y/sqrt(x+y)), (x,-10,10), (y,-10,10))
        sphinx_plot(g)

    Extra options will get passed on to show(), as long as they are valid::

        sage: plot_vector_field((x,y), (x,-2,2), (y,-2,2), xmax=10)
        Graphics object consisting of 1 graphics primitive
        sage: plot_vector_field((x,y), (x,-2,2), (y,-2,2)).show(xmax=10) # These are equivalent

    .. PLOT::

        x,y = var('x,y')
        g = plot_vector_field((x,y), (x,-2,2), (y,-2,2), xmax=10)
        sphinx_plot(g)
    """
    f, g = f_g
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    z, ranges = setup_for_eval_on_grid([f, g], [xrange, yrange],
                                       options['plot_points'])
    f, g = z

    xpos_array, ypos_array, xvec_array, yvec_array = [], [], [], []
    for x in xsrange(*ranges[0], include_endpoint=True):
        for y in xsrange(*ranges[1], include_endpoint=True):
            xpos_array.append(x)
            ypos_array.append(y)
            xvec_array.append(f(x, y))
            yvec_array.append(g(x, y))

    import numpy
    xvec_array = numpy.ma.masked_invalid(numpy.array(xvec_array, dtype=float))
    yvec_array = numpy.ma.masked_invalid(numpy.array(yvec_array, dtype=float))
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(PlotField(xpos_array, ypos_array,
                              xvec_array, yvec_array, options))
    return g


def plot_slope_field(f, xrange, yrange, **kwds):
    r"""
    ``plot_slope_field`` takes a function of two variables xvar and yvar
    (for instance, if the variables are `x` and `y`, take `f(x,y)`), and at
    representative points `(x_i,y_i)` between xmin, xmax, and ymin, ymax
    respectively, plots a line with slope `f(x_i,y_i)` (see below).

    ``plot_slope_field(f, (xvar,xmin,xmax), (yvar,ymin,ymax))``

    EXAMPLES:

    A logistic function modeling population growth::

        sage: x,y = var('x y')
        sage: capacity = 3 # thousand
        sage: growth_rate = 0.7 # population increases by 70% per unit of time
        sage: plot_slope_field(growth_rate * (1-y/capacity) * y, (x,0,5), (y,0,capacity*2))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x y')
        capacity = 3 # thousand
        growth_rate = 0.7 # population increases by 70% per unit of time
        g = plot_slope_field(growth_rate * (1-y/capacity) * y, (x,0,5), (y,0,capacity*2))
        sphinx_plot(g)

    Plot a slope field involving sin and cos::

        sage: x,y = var('x y')
        sage: plot_slope_field(sin(x+y) + cos(x+y), (x,-3,3), (y,-3,3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x y')
        g = plot_slope_field(sin(x+y)+cos(x+y), (x,-3,3), (y,-3,3))
        sphinx_plot(g)

    Plot a slope field using a lambda function::

        sage: plot_slope_field(lambda x,y: x + y, (-2,2), (-2,2))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x y')
        g = plot_slope_field(lambda x,y: x + y, (-2,2), (-2,2))
        sphinx_plot(g)

    TESTS:

    Verify that we're not getting warnings due to use of headless quivers
    (:issue:`11208`)::

        sage: x,y = var('x y')
        sage: import numpy # bump warnings up to errors for testing purposes
        sage: old_err = numpy.seterr('raise')
        sage: plot_slope_field(sin(x+y) + cos(x+y), (x,-3,3), (y,-3,3))
        Graphics object consisting of 1 graphics primitive
        sage: dummy_err = numpy.seterr(**old_err)
    """
    slope_options = {'headaxislength': 0,
                     'headlength': 1e-9,
                     'pivot': 'middle'}
    slope_options.update(kwds)

    from sage.misc.functional import sqrt
    from sage.misc.sageinspect import is_function_or_cython_function
    if is_function_or_cython_function(f):
        def norm_inverse(x, y):
            return 1 / sqrt(f(x, y)**2 + 1)

        def f_normalized(x, y):
            return f(x, y) * norm_inverse(x, y)
    else:
        norm_inverse = 1 / sqrt(f**2 + 1)
        f_normalized = f * norm_inverse
    return plot_vector_field((norm_inverse, f_normalized), xrange, yrange, **slope_options)
