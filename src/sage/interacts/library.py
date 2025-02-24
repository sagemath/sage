r"""
Sage Interacts

Sage interacts are applications of the `@interact decorator <../../sagenb/notebook/interact.html>`_.
They are conveniently accessible in the Sage Notebook via ``interacts.[TAB].[TAB]()``.
The first ``[TAB]`` lists categories and the second ``[TAB]`` reveals the interact examples.

EXAMPLES:

Invoked in the notebook, the following command will produce the fully formatted
interactive mathlet.  In the command line, it will simply return the underlying
HTML and Sage code which creates the mathlet::

    sage: interacts.calculus.taylor_polynomial()
    ...Interactive function <function taylor_polynomial at ...> with 3 widgets
      title: HTMLText(value='<h2>Taylor polynomial</h2>')
      f: EvalText(value='e^(-x)*sin(x)', description='$f(x)=$', layout=Layout(max_width='81em'))
      order: SelectionSlider(description='order', options=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), value=1)

AUTHORS:

- William Stein
- Harald Schilly, Robert Marik (2011-01-16): added many examples (#9623)
  partially based on work by Lauri Ruotsalainen
"""

# *****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Harald Schilly <harald.schilly@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from typing import Any, Callable

from sage.arith.misc import factor
from sage.arith.srange import srange
from sage.calculus.all import symbolic_expression
from sage.calculus.functional import derivative
from sage.calculus.integration import numerical_integral as integral_numerical
from sage.ext.fast_callable import fast_callable
from sage.functions.log import exp
from sage.misc.functional import sqrt
from sage.functions.trig import (acos, cos, sin, tan)
from sage.misc.decorators import sage_wraps
from sage.misc.functional import N
from sage.misc.latex import latex
from sage.misc.sage_eval import sage_eval
from sage.misc.table import table
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.circle", "circle")
lazy_import("sage.plot.complex_plot", "complex_plot")
lazy_import("sage.plot.disk", "disk")
lazy_import("sage.plot.graphics", "Graphics")
lazy_import("sage.plot.line", ["line", "line2d"])
lazy_import("sage.plot.matrix_plot", "matrix_plot")
lazy_import("sage.plot.plot", ["graphics_array", "parametric_plot", "plot"])
lazy_import("sage.plot.point", ["point", "points"])
lazy_import("sage.plot.polygon", "polygon2d")
lazy_import("sage.plot.text", "text")
from sage.repl.rich_output.pretty_print import (pretty_print, show)
from sage.rings.complex_double import CDF
from sage.rings.integer import Integer
from sage.symbolic.constants import pi
from sage.symbolic.relation import solve
from sage.symbolic.ring import SR

x = SR.var('x')

# It is important that this file is lazily imported for this to work
from sage.repl.user_globals import get_global

assert get_global  # to suppress pyflakes warning

from sage.repl.ipython_kernel.all_jupyter import (interact, checkbox,
    input_box, input_grid, range_slider, selector, slider, text_control)


def library_interact(
    decorator_target: Callable[..., Any] = None, **widgets: Callable[..., Any]
):
    r"""
    This is a decorator for using interacts in the Sage library.

    This is essentially a wrapper around the ``interact`` function.

    INPUT:

    - ``**widgets`` -- keyword arguments that are passed to the
      ``interact`` function to create the widgets. Each value must
      be a callable that returns a widget.

    EXAMPLES::

        sage: from sage.interacts.library import library_interact
        sage: @library_interact(n=lambda: slider(-5, 15, None, 5))
        ....: def f(n):
        ....:     print(n)
        sage: f()  # an interact appears if using the notebook, else code
        ...Interactive function <function f at ...> with 1 widget
        n: TransformIntSlider(value=5, description='n', max=15, min=-5)

    TESTS:

    Backwards compatibility::

        sage: from sage.interacts.library import library_interact
        sage: @library_interact
        ....: def f(n=slider(-5, 15, None, 5)):
        ....:     print(n)
        doctest:warning
        ...
        DeprecationWarning: Use decorator factory @library_interact(widgets) instead of @library_interact without any arguments.
        See https://github.com/sagemath/sage/issues/33382 for details.
        sage: f()  # an interact appears if using the notebook, else code
        ...Interactive function <function f at ...> with 1 widget
        n: TransformIntSlider(value=5, description='n', max=15, min=-5)

    .. NOTE::

        We use a callable to construct the widget so that the widget is only
        initialized when the function is called and not upon loading the module.
    """

    def decorator(f: Callable[..., Any]):
        @sage_wraps(f)
        def library_wrapper():
            widgets_constructed = {
                key: constructor() for key, constructor in widgets.items()
            }
            # This will display the interact, no need to return anything
            interact(f, **widgets_constructed)

        # Store the widget constructors for later use in tests
        library_wrapper._widgets = widgets
        return library_wrapper

    if decorator_target is None:
        return decorator
    else:
        from sage.misc.superseded import deprecation

        deprecation(
            33382,
            "Use decorator factory @library_interact(widgets) instead of @library_interact without any arguments.",
        )
        return decorator(decorator_target)


def html(obj):
    r"""
    Shorthand to pretty print HTML.

    EXAMPLES::

        sage: from sage.interacts.library import html
        sage: html("<h1>Hello world</h1>")
        <h1>Hello world</h1>
    """
    from sage.misc.html import html
    pretty_print(html(obj))


@library_interact(n=lambda: slider(range(10)), m=lambda: slider(range(10)))
def demo(n: int, m: int):
    r"""
    This is a demo interact that sums two numbers.

    INPUT:

    - ``n`` -- integer
    - ``m`` -- integer

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.demo()
        ...Interactive function <function demo at ...> with 2 widgets
          n: SelectionSlider(description='n', options=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), value=0)
          m: SelectionSlider(description='m', options=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), value=0)
    """
    print(n + m)


@library_interact(
    title=lambda: text_control("<h2>Taylor polynomial</h2>"),
    f=lambda: input_box(sin(x) * exp(-x), label="$f(x)=$"),
    order=lambda: slider(range(1, 13)),
)
def taylor_polynomial(title, f, order: int):
    r"""
    Illustrate the Taylor polynomial approximation
    of various orders around `x=0`.

    INPUT:

    - ``f`` -- function expression
    - ``order`` -- integer

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.taylor_polynomial()
        ...Interactive function <function taylor_polynomial at ...> with 3 widgets
          title: HTMLText(value='<h2>Taylor polynomial</h2>')
          f: EvalText(value='e^(-x)*sin(x)', description='$f(x)=$', layout=Layout(max_width='81em'))
          order: SelectionSlider(description='order', options=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), value=1)
    """
    x0 = 0
    p = plot(f, (x, -1, 5), thickness=2)
    dot = point((x0, f(x=x0)), pointsize=80, rgbcolor=(1, 0, 0))
    ft = f.taylor(x, x0, order)
    pt = plot(ft, (-1, 5), color='green', thickness=2)
    html(r'$f(x)\;=\;%s$' % latex(f))
    html(r'$\hat{f}(x;%s)\;=\;%s+\mathcal{O}(x^{%s})$' % (x0, latex(ft),
                                                          order + 1))
    show(dot + p + pt, ymin=-.5, ymax=1)


@library_interact(
    title=lambda: text_control("<h2>Definite integral</h2>"),
    f=lambda: input_box(default='3*x', label="$f(x)=$"),
    g=lambda: input_box(default='x^2', label="$g(x)=$"),
    interval=lambda: range_slider(-10, 10, default=(0, 3), label='Interval'),
    x_range=lambda: range_slider(-10, 10, default=(0, 3), label="plot range (x)"),
    selection=lambda: selector(
        ["f", "g", "f and g", "f - g"], default="f and g", label="Select"
    ),
)
def definite_integral(title, f, g, interval, x_range, selection):
    r"""
    This is a demo interact for plotting the definite integral of a function
    based on work by Lauri Ruotsalainen, 2010.

    INPUT:

    - ``function`` -- input box, function in x
    - ``interval`` -- interval for the definite integral
    - ``x_range`` -- range slider for plotting range
    - ``selection`` -- selector on how to visualize the integrals

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.definite_integral()
        ...Interactive function <function definite_integral at ...> with 6 widgets
          title: HTMLText(value='<h2>Definite integral</h2>')
          f: EvalText(value='3*x', description='$f(x)=$', layout=Layout(max_width='81em'))
          g: EvalText(value='x^2', description='$g(x)=$', layout=Layout(max_width='81em'))
          interval: IntRangeSlider(value=(0, 3), description='Interval', max=10, min=-10)
          x_range: IntRangeSlider(value=(0, 3), description='plot range (x)', max=10, min=-10)
          selection: Dropdown(description='Select', index=2,
                              options=('f', 'g', 'f and g', 'f - g'), value='f and g')
    """
    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    g = symbolic_expression(g).function(x)
    f_plot = Graphics()
    g_plot = Graphics()
    h_plot = Graphics()
    text = ""

    # Plot function f.
    if selection != "g":
        f_plot = plot(f(x), x, x_range, color='blue', thickness=1.5)

    # Color and calculate the area between f and the horizontal axis.
    if selection == "f" or selection == "f and g":
        f_plot += plot(f(x), x, interval, color='blue', fill=True, fillcolor='blue', fillalpha=0.15)
        text += r"$\int_{%.2f}^{%.2f}(\color{Blue}{f(x)})\,\mathrm{d}x=\int_{%.2f}^{%.2f}(%s)\,\mathrm{d}x=%.2f$" % (
            interval[0], interval[1],
            interval[0], interval[1],
            latex(f(x)),
            f(x).nintegrate(x, interval[0], interval[1])[0]
        )

    if selection == "f and g":
        text += r"<br/>"

    # Plot function g. Also color and calculate the area between g and the horizontal axis.
    if selection == "g" or selection == "f and g":
        g_plot = plot(g(x), x, x_range, color='green', thickness=1.5)
        g_plot += plot(g(x), x, interval, color='green', fill=True, fillcolor='yellow', fillalpha=0.5)
        text += r"$\int_{%.2f}^{%.2f}(\color{Green}{g(x)})\,\mathrm{d}x=\int_{%.2f}^{%.2f}(%s)\,\mathrm{d}x=%.2f$" % (
            interval[0], interval[1],
            interval[0], interval[1],
            latex(g(x)),
            g(x).nintegrate(x, interval[0], interval[1])[0]
        )

    # Plot function f-g. Also color and calculate the area between f-g and the horizontal axis.
    if selection == "f - g":
        g_plot = plot(g(x), x, x_range, color='green', thickness=1.5)
        g_plot += plot(g(x), x, interval, color='green', fill=f(x), fillcolor='red', fillalpha=0.15)
        h_plot = plot(f(x)-g(x), x, interval, color='red', thickness=1.5, fill=True, fillcolor='red', fillalpha=0.15)
        text = r"$\int_{%.2f}^{%.2f}(\color{Red}{f(x)-g(x)})\,\mathrm{d}x=\int_{%.2f}^{%.2f}(%s)\,\mathrm{d}x=%.2f$" % (
            interval[0], interval[1],
            interval[0], interval[1],
            latex(f(x)-g(x)),
            (f(x)-g(x)).nintegrate(x, interval[0], interval[1])[0]
        )

    show(f_plot + g_plot + h_plot, gridlines=True)
    html(text)


@library_interact(
    title=lambda: text_control("<h2>Derivative grapher</h2>"),
    function=lambda: input_box(default='x^5-3*x^3+1', label='Function:'),
    x_range=lambda: range_slider(-15, 15, 0.1, default=(-2, 2), label="Range (x)"),
    y_range=lambda: range_slider(-15, 15, 0.1, default=(-8, 6), label="Range (y)"),
)
def function_derivative(title, function, x_range, y_range):
    r"""
    This is a demo interact for plotting derivatives of a function based on work by
    Lauri Ruotsalainen, 2010.

    INPUT:

        - ``function`` -- input box, function in x
        - ``x_range`` -- range slider for plotting range
        - ``y_range`` -- range slider for plotting range

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.function_derivative()
        ...Interactive function <function function_derivative at ...> with 4 widgets
          title: HTMLText(value='<h2>Derivative grapher</h2>')
          function: EvalText(value='x^5-3*x^3+1', description='Function:', layout=Layout(max_width='81em'))
          x_range: FloatRangeSlider(value=(-2.0, 2.0), description='Range (x)', max=15.0, min=-15.0)
          y_range: FloatRangeSlider(value=(-8.0, 6.0), description='Range (y)', max=15.0, min=-15.0)
    """
    x = SR.var('x')
    f = symbolic_expression(function).function(x)
    df = derivative(f, x)
    ddf = derivative(df, x)
    plots = plot(f(x), x_range, thickness=1.5) + plot(df(x), x_range, color='green') + plot(ddf(x), x_range, color='red')
    if y_range == (0,0):
        show(plots, xmin=x_range[0], xmax=x_range[1])
    else:
        show(plots, xmin=x_range[0], xmax=x_range[1], ymin=y_range[0], ymax=y_range[1])

    html(r"<center>$\color{Blue}{f(x) = %s}$</center>" % latex(f(x)))
    html(r"<center>$\color{Green}{f'(x) = %s}$</center>" % latex(df(x)))
    html(r"<center>$\color{Red}{f''(x) = %s}$</center>" % latex(ddf(x)))


@library_interact(
    title=lambda: text_control("<h2>Difference quotient</h2>"),
    f=lambda: input_box(default='sin(x)', label="f(x)"),
    interval=lambda: range_slider(0, 10, 0.1, default=(0.0, 10.0), label='Range'),
    a=lambda: slider(0, 10, None, 5.5, label='$a$'),
    x0=lambda: slider(0, 10, None, 2.5, label="$x_0$ (start point)"),
)
def difference_quotient(title, f, interval, a, x0):
    r"""
    This is a demo interact for difference quotient based on work by
    Lauri Ruotsalainen, 2010.

    INPUT:

    - ``f`` -- input box, function in `x`
    - ``interval`` -- range slider for plotting
    - ``a`` -- slider for `a`
    - ``x0`` -- slider for starting point `x_0`

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.difference_quotient()
        ...Interactive function <function difference_quotient at ...> with 5 widgets
          title: HTMLText(value='<h2>Difference quotient</h2>')
          f: EvalText(value='sin(x)', description='f(x)', layout=Layout(max_width='81em'))
          interval: FloatRangeSlider(value=(0.0, 10.0), description='Range', max=10.0)
          a: IntSlider(value=5, description='$a$', max=10)
          x0: IntSlider(value=2, description='$x_0$ (start point)', max=10)
    """
    html('<h2>Difference Quotient</h2>')
    html('<div style="white-space: normal;">\
         <a href="https://en.wikipedia.org/wiki/Difference_quotient" target="_blank">\
         Wikipedia article about difference quotient</a></div>'
         )

    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    fmax = f.find_local_maximum(interval[0], interval[1])[0]
    fmin = f.find_local_minimum(interval[0], interval[1])[0]
    f_height = fmax - fmin
    measure_y = fmin - 0.1*f_height

    measure_0 = line2d([(x0, measure_y), (a, measure_y)], rgbcolor='black')
    measure_1 = line2d([(x0, measure_y + 0.02*f_height), (x0, measure_y-0.02*f_height)], rgbcolor='black')
    measure_2 = line2d([(a, measure_y + 0.02*f_height), (a, measure_y-0.02*f_height)], rgbcolor='black')
    text_x0 = text("x0", (x0, measure_y - 0.05*f_height), rgbcolor='black')
    text_a = text("a", (a, measure_y - 0.05*f_height), rgbcolor='black')
    measure = measure_0 + measure_1 + measure_2 + text_x0 + text_a

    tanf = symbolic_expression((f(x0)-f(a))*(x-a)/(x0-a)+f(a)).function(x)

    fplot = plot(f(x), x, interval[0], interval[1])
    tanplot = plot(tanf(x), x, interval[0], interval[1], rgbcolor='#FF0000')
    points = point([(x0, f(x0)), (a, f(a))], pointsize=20, rgbcolor='#005500')
    dashline = line2d([(x0, f(x0)), (x0, f(a)), (a, f(a))], rgbcolor='#005500', linestyle='--')
    html('<h2>Difference Quotient</h2>')
    show(fplot + tanplot + points + dashline + measure, xmin=interval[0], xmax=interval[1], ymin=fmin-0.2*f_height, ymax=fmax)
    html(r"<br>$\text{Line's equation:}$")
    html(r"$y = %s$<br>" % tanf(x))
    html(r"$\text{Slope:}$")
    html(r"$k = \frac{f(x_0)-f(a)}{x_0-a} = %s$<br>" % (N(derivative(tanf(x), x), digits=5)))


@library_interact(
    A=lambda: slider(-7, 7, 1, 1),
    B=lambda: slider(-7, 7, 1, 1),
    C=lambda: slider(-7, 7, 1, -2),
)
def quadratic_equation(A, B, C):
    r"""
    This is a demo interact for solving quadratic equations based on work by
    Lauri Ruotsalainen, 2010.

    INPUT:

    - ``A`` -- integer slider
    - ``B`` -- integer slider
    - ``C`` -- integer slider

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.quadratic_equation()
        ...Interactive function <function quadratic_equation at ...> with 3 widgets
          A: IntSlider(value=1, description='A', max=7, min=-7)
          B: IntSlider(value=1, description='B', max=7, min=-7)
          C: IntSlider(value=-2, description='C', max=7, min=-7)
    """
    x = SR.var('x')
    f = symbolic_expression(A*x**2 + B*x + C).function(x)
    html('<h2>The Solutions of the Quadratic Equation</h2>')
    html("$%s = 0$" % f(x))

    show(plot(f(x), x, (-10, 10), ymin=-10, ymax=10), aspect_ratio=1, figsize=4)

    d = B**2 - 4*A*C

    if d < 0:
        color = "Red"
        sol = r"\text{solution} \in \mathbb{C}"
    elif d == 0:
        color = "Blue"
        sol = -B/(2*A)
    else:
        color = "Green"
        a = (-B+sqrt(B**2-4*A*C))/(2*A)
        b = (-B-sqrt(B**2-4*A*C))/(2*A)
        sol = r"\begin{cases}%s\\%s\end{cases}" % (latex(a), latex(b))

    if B < 0:
        dis1 = "(%s)^2-4*%s*%s" % (B, A, C)
    else:
        dis1 = "%s^2-4*%s*%s" % (B, A, C)
    dis2 = r"\color{%s}{%s}" % (color, d)

    html("$Ax^2 + Bx + C = 0$")
    calc = r"$x = \frac{-B\pm\sqrt{B^2-4AC}}{2A} = " + \
           r"\frac{-%s\pm\sqrt{%s}}{2*%s} = " + \
           r"\frac{-%s\pm\sqrt{%s}}{%s} = %s$"
    html(calc % (B, dis1, A, B, dis2, (2*A), sol))

@library_interact(
    a0=lambda: slider(0, 360, 1, 30, label='A'),
    a1=lambda: slider(0, 360, 1, 180, label='B'),
    a2=lambda: slider(0, 360, 1, 300, label='C'),
)
def trigonometric_properties_triangle(a0, a1, a2):
    r"""
    This is an interact for demonstrating trigonometric properties
    in a triangle based on work by Lauri Ruotsalainen, 2010.

    INPUT:

    - ``a0`` -- angle
    - ``a1`` -- angle
    - ``a2`` -- angle

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.geometry.trigonometric_properties_triangle()
        ...Interactive function <function trigonometric_properties_triangle at ...> with 3 widgets
          a0: IntSlider(value=30, description='A', max=360)
          a1: IntSlider(value=180, description='B', max=360)
          a2: IntSlider(value=300, description='C', max=360)
    """
    import math

    # Returns the distance between points (x1,y1) and (x2,y2)
    def distance(x1_y1, x2_y2):
        (x1, y1) = x1_y1
        (x2, y2) = x2_y2
        return sqrt((x2-x1)**2 + (y2-y1)**2)

    # Returns an angle (in radians) when sides a and b
    # are adjacent and the side c is opposite to the angle
    def angle(a, b, c):
        a,b,c = map(float,[a,b,c])
        return acos(0.5 * (b**2 + c**2 - a**2) / (b * c))

    # Returns the area of a triangle when an angle alpha
    # and adjacent sides a and b are known
    def area(alpha, a, b):
        return 0.5 * a * b * sin(alpha)

    xy = [0]*3
    html('<h2>Trigonometric Properties of a Triangle</h2>')
    # Coordinates of the angles
    a = [math.radians(float(x)) for x in [a0, a1, a2]]
    for i in range(3):
        xy[i] = (cos(a[i]), sin(a[i]))

    # Side lengths (bc, ca, ab) corresponding to triangle vertices (a, b, c)
    al = [distance(xy[1], xy[2]), distance(xy[2], xy[0]), distance(xy[0], xy[1])]

    # The angles (a, b, c) in radians
    ak = [angle(al[0], al[1], al[2]), angle(al[1], al[2], al[0]), angle(al[2], al[0], al[1])]

    # The area of the triangle
    A = area(ak[0], al[1], al[2])

    unit_circle = circle((0, 0), 1, aspect_ratio=1)

    # Triangle
    triangle = line([xy[0], xy[1], xy[2], xy[0]], rgbcolor='black')
    triangle_points = point(xy, pointsize=30)

    # Labels of the angles drawn in a distance from points
    a_label = text("A", (xy[0][0]*1.07, xy[0][1]*1.07))
    b_label = text("B", (xy[1][0]*1.07, xy[1][1]*1.07))
    c_label = text("C", (xy[2][0]*1.07, xy[2][1]*1.07))
    labels = a_label + b_label + c_label

    show(unit_circle + triangle + triangle_points + labels, figsize=[5, 5], xmin=-1, xmax=1, ymin=-1, ymax=1)
    html(r"$\angle A = {%.3f}^{\circ},$ $\angle B = {%.3f}^{\circ},$ $\angle C = {%.3f}^{\circ}$"
         % (math.degrees(ak[0]), math.degrees(ak[1]), math.degrees(ak[2])))
    html(r"$AB = %.6f$, $BC = %.6f$, $CA = %.6f$" % (al[2], al[0], al[1]))
    html(r"Area of triangle $ABC = %.6f$" % A)


@library_interact(
    function=lambda: selector([(0, sin(x)), (1, cos(x)), (2, tan(x))]),
    x=lambda: slider(0, 2 * pi, 0.005 * pi, 0),
)
def unit_circle(function, x):
    r"""
    This is an interact for Sin, Cos and Tan in the Unit Circle
    based on work by Lauri Ruotsalainen, 2010.

    INPUT:

    - ``function`` -- select Sin, Cos or Tan
    - ``x`` -- slider to select angle in unit circle

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.geometry.unit_circle()
        ...Interactive function <function unit_circle at ...> with 2 widgets
          function: Dropdown(description='function', options=(('sin(x)', 0), ('cos(x)', 1), ('tan(x)', 2)), value=0)
          x: TransformFloatSlider(value=0.0, description='x', max=6.283185307179586, step=0.015707963267948967)
    """
    xy = (cos(x), sin(x))
    t = SR.var('t')
    html('<div style="white-space: normal;">Lines of the same color have\
         the same length</div>')

    # Unit Circle
    C = circle((0, 0), 1, figsize=[5, 5], aspect_ratio=1)
    C_line = line([(0, 0), (xy[0], xy[1])], rgbcolor='black')
    C_point = point((xy[0], xy[1]), pointsize=40, rgbcolor='green')
    C_inner = parametric_plot((cos(t), sin(t)), (t, 0, x + 0.001), color='green', thickness=3)
    C_outer = parametric_plot((0.1 * cos(t), 0.1 * sin(t)), (t, 0, x + 0.001), color='black')
    C_graph = C + C_line + C_point + C_inner + C_outer

    # Graphics related to the graph of the function
    G_line = line([(0, 0), (x, 0)], rgbcolor='green', thickness=3)
    G_point = point((x, 0), pointsize=30, rgbcolor='green')
    G_graph = G_line + G_point

    # Sine
    if function == 0:
        Gf = plot(sin(t), t, 0, 2*pi, axes_labels=("x", "sin(x)"))
        Gf_point = point((x, sin(x)), pointsize=30, rgbcolor='red')
        Gf_line = line([(x, 0),(x, sin(x))], rgbcolor='red')
        Cf_point = point((0, xy[1]), pointsize=40, rgbcolor='red')
        Cf_line1 = line([(0, 0), (0, xy[1])], rgbcolor='red', thickness=3)
        Cf_line2 = line([(0, xy[1]), (xy[0], xy[1])], rgbcolor='purple', linestyle='--')
    # Cosine
    elif function == 1:
        Gf = plot(cos(t), t, 0, 2*pi, axes_labels=("x", "cos(x)"))
        Gf_point = point((x, cos(x)), pointsize=30, rgbcolor='red')
        Gf_line = line([(x, 0), (x, cos(x))], rgbcolor='red')
        Cf_point = point((xy[0], 0), pointsize=40, rgbcolor='red')
        Cf_line1 = line([(0, 0), (xy[0], 0)], rgbcolor='red', thickness=3)
        Cf_line2 = line([(xy[0], 0), (xy[0], xy[1])], rgbcolor='purple', linestyle='--')
    # Tangent
    else:
        Gf = plot(tan(t), t, 0, 2*pi, ymin=-8, ymax=8, axes_labels=("x", "tan(x)"))
        Gf_point = point((x, tan(x)), pointsize=30, rgbcolor='red')
        Gf_line = line([(x, 0), (x, tan(x))], rgbcolor='red')
        Cf_point = point((1, tan(x)), pointsize=40, rgbcolor='red')
        Cf_line1 = line([(1, 0), (1, tan(x))], rgbcolor='red', thickness=3)
        Cf_line2 = line([(xy[0], xy[1]), (1, tan(x))], rgbcolor='purple', linestyle='--')

    C_graph += Cf_point + Cf_line1 + Cf_line2
    G_graph += Gf + Gf_point + Gf_line

    show(graphics_array([C_graph, G_graph]))


@library_interact(
    title=lambda: text_control("<h2>Special points in triangle</h2>"),
    a0=lambda: slider(0, 360, 1, 30, label='A'),
    a1=lambda: slider(0, 360, 1, 180, label='B'),
    a2=lambda: slider(0, 360, 1, 300, label='C'),
    show_median=lambda: checkbox(False, label='Medians'),
    show_pb=lambda: checkbox(False, label="Perpendicular Bisectors"),
    show_alt=lambda: checkbox(False, label='Altitudes'),
    show_ab=lambda: checkbox(False, label="Angle Bisectors"),
    show_incircle=lambda: checkbox(False, label='Incircle'),
    show_euler=lambda: checkbox(False, label="Euler's Line"),
)
def special_points(
    title,
    a0,
    a1,
    a2,
    show_median,
    show_pb,
    show_alt,
    show_ab,
    show_incircle,
    show_euler,
):
    r"""
    This interact demo shows special points in a triangle
    based on work by Lauri Ruotsalainen, 2010.

    INPUT:

    - ``a0`` -- angle
    - ``a1`` -- angle
    - ``a2`` -- angle
    - ``show_median`` -- checkbox
    - ``show_pb`` -- checkbox to show perpendicular bisectors
    - ``show_alt`` -- checkbox to show altitudes
    - ``show_ab`` -- checkbox to show angle bisectors
    - ``show_incircle`` -- checkbox to show incircle
    - ``show_euler`` -- checkbox to show euler's line

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.geometry.special_points()
        ...Interactive function <function special_points at ...> with 10 widgets
          title: HTMLText(value='<h2>Special points in triangle</h2>')
          a0: IntSlider(value=30, description='A', max=360)
          a1: IntSlider(value=180, description='B', max=360)
          a2: IntSlider(value=300, description='C', max=360)
          show_median: Checkbox(value=False, description='Medians')
          show_pb: Checkbox(value=False, description='Perpendicular Bisectors')
          show_alt: Checkbox(value=False, description='Altitudes')
          show_ab: Checkbox(value=False, description='Angle Bisectors')
          show_incircle: Checkbox(value=False, description='Incircle')
          show_euler: Checkbox(value=False, description="Euler's Line")
    """
    import math
    # Return the intersection point of the bisector of the angle <(A[a],A[c],A[b]) and the unit circle. Angles given in radians.

    def half(A, a, b, c):
        if (A[a] < A[b] and (A[c] < A[a] or A[c] > A[b])) or (A[a] > A[b] and (A[c] > A[a] or A[c] < A[b])):
            p = A[a] + 0.5 * (A[b] - A[a])
        else:
            p = A[b] + 0.5 * (2*pi - (A[b]-A[a]))
        return (math.cos(p), math.sin(p))

    # Returns the distance between points (x1,y1) and (x2,y2)
    def distance(x1_y1, x2_y2):
        (x1, y1) = x1_y1
        (x2, y2) = x2_y2
        return math.sqrt((x2-x1)**2 + (y2-y1)**2)

    # Returns the line (graph) going through points (x1,y1) and (x2,y2)
    def line_to_points(x1_y1, x2_y2, **plot_kwargs):
        (x1, y1) = x1_y1
        (x2, y2) = x2_y2
        return plot((y2-y1) / (x2-x1) * (x-x1) + y1, (x,-3,3), **plot_kwargs)

    # Coordinates of the angles
    a = [math.radians(float(x)) for x in [a0, a1, a2]]
    xy = [(math.cos(a[i]), math.sin(a[i])) for i in range(3)]

    # Labels of the angles drawn in a distance from points
    a_label = text("A", (xy[0][0]*1.07, xy[0][1]*1.07))
    b_label = text("B", (xy[1][0]*1.07, xy[1][1]*1.07))
    c_label = text("C", (xy[2][0]*1.07, xy[2][1]*1.07))
    labels = a_label + b_label + c_label

    C = circle((0, 0), 1, aspect_ratio=1)

    # Triangle
    triangle = line([xy[0], xy[1], xy[2], xy[0]], rgbcolor='black')
    triangle_points = point(xy, pointsize=30)

    # Side lengths (bc, ca, ab) corresponding to triangle vertices (a, b, c)
    ad = [distance(xy[1], xy[2]), distance(xy[2], xy[0]), distance(xy[0], xy[1])]

    # Midpoints of edges (bc, ca, ab)
    a_middle = [
        (0.5 * (xy[1][0] + xy[2][0]), 0.5 * (xy[1][1] + xy[2][1])),
        (0.5 * (xy[2][0] + xy[0][0]), 0.5 * (xy[2][1] + xy[0][1])),
        (0.5 * (xy[0][0] + xy[1][0]), 0.5 * (xy[0][1] + xy[1][1]))
    ]

    # Incircle
    perimeter = float(ad[0] + ad[1] + ad[2])
    incircle_center = (
        (ad[0]*xy[0][0] + ad[1]*xy[1][0] + ad[2]*xy[2][0]) / perimeter,
        (ad[0]*xy[0][1] + ad[1]*xy[1][1] + ad[2]*xy[2][1]) / perimeter
    )

    if show_incircle:
        s = 0.5 * perimeter
        incircle_r = math.sqrt((s - ad[0]) * (s - ad[1]) * (s - ad[2]) / s)
        incircle_graph = circle(incircle_center, incircle_r) + point(incircle_center)
    else:
        incircle_graph = Graphics()

    # Angle Bisectors
    if show_ab:
        a_ab = line([xy[0], half(a, 1, 2, 0)], rgbcolor='blue', alpha=0.6)
        b_ab = line([xy[1], half(a, 2, 0, 1)], rgbcolor='blue', alpha=0.6)
        c_ab = line([xy[2], half(a, 0, 1, 2)], rgbcolor='blue', alpha=0.6)
        ab_point = point(incircle_center, rgbcolor='blue', pointsize=28)
        ab_graph = a_ab + b_ab + c_ab + ab_point
    else:
        ab_graph = Graphics()

    # Medians
    if show_median:
        a_median = line([xy[0], a_middle[0]], rgbcolor='green', alpha=0.6)
        b_median = line([xy[1], a_middle[1]], rgbcolor='green', alpha=0.6)
        c_median = line([xy[2], a_middle[2]], rgbcolor='green', alpha=0.6)
        median_point = point(
            (
                (xy[0][0]+xy[1][0]+xy[2][0])/3.0,
                (xy[0][1]+xy[1][1]+xy[2][1])/3.0
            ), rgbcolor='green', pointsize=28)
        median_graph = a_median + b_median + c_median + median_point
    else:
        median_graph = Graphics()

    # Perpendicular Bisectors
    if show_pb:
        a_pb = line_to_points(a_middle[0], half(a, 1, 2, 0), rgbcolor='red', alpha=0.6)
        b_pb = line_to_points(a_middle[1], half(a, 2, 0, 1), rgbcolor='red', alpha=0.6)
        c_pb = line_to_points(a_middle[2], half(a, 0, 1, 2), rgbcolor='red', alpha=0.6)
        pb_point = point((0, 0), rgbcolor='red', pointsize=28)
        pb_graph = a_pb + b_pb + c_pb + pb_point
    else:
        pb_graph = Graphics()

    # Altitudes
    if show_alt:
        xA, xB, xC = xy[0][0], xy[1][0], xy[2][0]
        yA, yB, yC = xy[0][1], xy[1][1], xy[2][1]
        a_alt = plot(((xC-xB)*x+(xB-xC)*xA)/(yB-yC)+yA, (x,-3,3), rgbcolor='brown', alpha=0.6)
        b_alt = plot(((xA-xC)*x+(xC-xA)*xB)/(yC-yA)+yB, (x,-3,3), rgbcolor='brown', alpha=0.6)
        c_alt = plot(((xB-xA)*x+(xA-xB)*xC)/(yA-yB)+yC, (x,-3,3), rgbcolor='brown', alpha=0.6)
        alt_lx = (xA*xB*(yA-yB)+xB*xC*(yB-yC)+xC*xA*(yC-yA)-(yA-yB)*(yB-yC)*(yC-yA))/(xC*yB-xB*yC+xA*yC-xC*yA+xB*yA-xA*yB)
        alt_ly = (yA*yB*(xA-xB)+yB*yC*(xB-xC)+yC*yA*(xC-xA)-(xA-xB)*(xB-xC)*(xC-xA))/(yC*xB-yB*xC+yA*xC-yC*xA+yB*xA-yA*xB)
        alt_intersection = point((alt_lx, alt_ly), rgbcolor='brown', pointsize=28)
        alt_graph = a_alt + b_alt + c_alt + alt_intersection
    else:
        alt_graph = Graphics()

    # Euler's Line
    if show_euler:
        euler_graph = line_to_points(
            (0, 0),
            (
                (xy[0][0]+xy[1][0]+xy[2][0])/3.0,
                (xy[0][1]+xy[1][1]+xy[2][1])/3.0
            ),
            rgbcolor='purple',
            thickness=2,
            alpha=0.7
        )
    else:
        euler_graph = Graphics()

    show(
        C + triangle + triangle_points + labels + ab_graph + median_graph +
        pb_graph + alt_graph + incircle_graph + euler_graph,
        figsize=[5,5], xmin=-1, xmax=1, ymin=-1, ymax=1
    )


@library_interact(
    n=lambda: slider(2, 10000, 100, default=1000, label="Number of Tosses"),
    interval=lambda: range_slider(
        0.0, 1.0, default=(0.45, 0.55), label="Plotting range (y)"
    ),
)
def coin(n, interval):
    r"""
    This interact demo simulates repeated tosses of a coin,
    based on work by Lauri Ruotsalainen, 2010.

    The points give the cumulative percentage of tosses which
    are heads in a given run of the simulation, so that the
    point `(x,y)` gives the percentage of the first `x` tosses
    that were heads; this proportion should approach .5, of
    course, if we are simulating a fair coin.

    INPUT:

    - ``n`` -- number of tosses
    - ``interval`` -- plot range along vertical axis

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.statistics.coin()
        ...Interactive function <function coin at ...> with 2 widgets
          n: IntSlider(value=1000, description='Number of Tosses', max=10000, min=2, step=100)
          interval: FloatRangeSlider(value=(0.45, 0.55), description='Plotting range (y)', max=1.0)
    """
    from random import random
    c = []
    k = 0.0
    for i in range(1, n + 1):
        k += random()
        c.append((i, k/i))
    show(point(c[1:], gridlines=[None, [0.5]], pointsize=1), ymin=interval[0], ymax=interval[1])


@library_interact(
    title=lambda: text_control("<h2>Bisection method</h2>"),
    f=lambda: input_box("x^2-2", label="f(x)"),
    interval=lambda: range_slider(-5, 5, default=(0, 4), label='range'),
    d=lambda: slider(1, 8, 1, 3, label="$10^{-d}$ precision"),
    maxn=lambda: slider(0, 50, 1, 10, label="max iterations"),
)
def bisection_method(title, f, interval, d, maxn):
    r"""
    Interact explaining the bisection method, based on similar interact
    explaining secant method and Wiliam Stein's example from wiki.

    INPUT:

    - ``f`` -- function
    - ``interval`` -- range slider for the search interval
    - ``d`` -- slider for the precision (`10^{-d}`)
    - ``maxn`` -- max number of iterations

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.secant_method()
        ...Interactive function <function secant_method at ...> with 5 widgets
          title: HTMLText(value='<h2>Secant method for numerical root finding</h2>')
          f: EvalText(value='x^2-2', description='f(x)', layout=Layout(max_width='81em'))
          interval: IntRangeSlider(value=(0, 4), description='range', max=5, min=-5)
          d: IntSlider(value=3, description='10^-d precision', max=16, min=1)
          maxn: IntSlider(value=10, description='max iterations', max=15)
    """
    def _bisection_method(f, a, b, maxn, eps):
        intervals = [(a,b)]
        round = 1
        two = float(2)
        while True:
            c = (b+a)/two
            if abs(f(c)) < h or round >= maxn:
                break
            fa = f(a)
            fb = f(b)
            fc = f(c)
            if abs(fc) < eps:
                return c, intervals
            if fa*fc < 0:
                a, b = a, c
            elif fc*fb < 0:
                a, b = c, b
            else:
                raise ValueError("f must have a sign change in the interval (%s,%s)" % (a,b))
            intervals.append((a,b))
            round += 1
        return c, intervals

    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    a, b = interval
    h = 10**(-d)
    try:
        c, intervals = _bisection_method(f, float(a), float(b), maxn, h)
    except ValueError:
        print("f must have opposite sign at the endpoints of the interval")
        show(plot(f, a, b, color='red'), xmin=a, xmax=b)
    else:
        html(r"$\text{Precision }h = 10^{-d}=10^{-%s}=%.5f$" % (d, float(h)))
        html(r"${c = }%s$" % latex(c))
        html(r"${f(c) = }%s" % latex(f(c)))
        html(r"$%s \text{ iterations}" % len(intervals))
        P = plot(f, a, b, color='red')
        k = (P.ymax() - P.ymin()) / (1.5*len(intervals))
        L = sum(line([(c,k*i), (d,k*i)]) for i, (c,d) in enumerate(intervals) )
        L += sum(line([(c,k*i-k/4), (c,k*i+k/4)]) for i, (c,d) in enumerate(intervals) )
        L += sum(line([(d,k*i-k/4), (d,k*i+k/4)]) for i, (c,d) in enumerate(intervals) )
        show(P + L, xmin=a, xmax=b)


@library_interact(
    title=lambda: text_control("<h2>Secant method for numerical root finding</h2>"),
    f=lambda: input_box("x^2-2", label="f(x)"),
    interval=lambda: range_slider(-5, 5, default=(0, 4), label='range'),
    d=lambda: slider(1, 16, 1, 3, label="10^-d precision"),
    maxn=lambda: slider(0, 15, 1, 10, label="max iterations"),
)
def secant_method(title, f, interval, d, maxn):
    r"""
    Interact explaining the secant method, based on work by
    Lauri Ruotsalainen, 2010.
    Originally this is based on work by William Stein.

    INPUT:

    - ``f`` -- function
    - ``interval`` -- range slider for the search interval
    - ``d`` -- slider for the precision (10^-d)
    - ``maxn`` -- max number of iterations

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.secant_method()
        ...Interactive function <function secant_method at ...> with 5 widgets
          title: HTMLText(value='<h2>Secant method for numerical root finding</h2>')
          f: EvalText(value='x^2-2', description='f(x)', layout=Layout(max_width='81em'))
          interval: IntRangeSlider(value=(0, 4), description='range', max=5, min=-5)
          d: IntSlider(value=3, description='10^-d precision', max=16, min=1)
          maxn: IntSlider(value=10, description='max iterations', max=15)
    """
    def _secant_method(f, a, b, maxn, h):
        intervals = [(a,b)]
        round = 1
        while True:
            c = b-(b-a)*f(b)/(f(b)-f(a))
            if abs(f(c)) < h or round >= maxn:
                break
            a, b = b, c
            intervals.append((a,b))
            round += 1
        return c, intervals

    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    a, b = interval
    h = 10**(-d)
    if float(f(a)*f(b)) > 0:
        print("f must have opposite sign at the endpoints of the interval")
        show(plot(f, a, b, color='red'), xmin=a, xmax=b)
    else:
        c, intervals = _secant_method(f, float(a), float(b), maxn, h)
        html(r"$\text{Precision }h = 10^{-d}=10^{-%s}=%.5f$" % (d, float(h)))
        html(r"${c = }%s$" % latex(c))
        html(r"${f(c) = }%s" % latex(f(c)))
        html(r"$%s \text{ iterations}" % len(intervals))
        P = plot(f, a, b, color='red')
        k = (P.ymax() - P.ymin()) / (1.5*len(intervals))
        L = sum(line([(c,k*i), (d,k*i)]) for i, (c,d) in enumerate(intervals) )
        L += sum(line([(c,k*i-k/4), (c,k*i+k/4)]) for i, (c,d) in enumerate(intervals) )
        L += sum(line([(d,k*i-k/4), (d,k*i+k/4)]) for i, (c,d) in enumerate(intervals) )
        S = sum(line([(c,f(c)), (d,f(d)), (d-(d-c)*f(d)/(f(d)-f(c)), 0)], color='green') for (c, d) in intervals)
        show(P + L + S, xmin=a, xmax=b)


@library_interact(
    title=lambda: text_control("<h2>Newton method</h2>"),
    f=lambda: input_box("x^2 - 2"),
    c=lambda: slider(-10, 10, default=6, label="Start ($x$)"),
    d=lambda: slider(1, 16, 1, 3, label="$10^{-d}$ precision"),
    maxn=lambda: slider(0, 15, 1, 10, label="max iterations"),
    interval=lambda: range_slider(-10, 10, default=(0, 6), label='Interval'),
    list_steps=lambda: checkbox(default=False, label="List steps"),
)
def newton_method(title, f, c, d, maxn, interval, list_steps):
    r"""
    Interact explaining the Newton method, based on work by
    Lauri Ruotsalainen, 2010.
    Originally this is based on work by William Stein.

    INPUT:

    - ``f`` -- function
    - ``c`` -- starting position (`x`)
    - ``d`` -- slider for the precision (`10^{-d}`)
    - ``maxn`` -- max number of iterations
    - ``interval`` -- range slider for the search interval
    - ``list_steps`` -- checkbox, if ``True`` shows the steps numerically

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.newton_method()
        ...Interactive function <function newton_method at ...> with 7 widgets
          title: HTMLText(value='<h2>Newton method</h2>')
          f: EvalText(value='x^2 - 2', description='f', layout=Layout(max_width='81em'))
          c: IntSlider(value=6, description='Start ($x$)', max=10, min=-10)
          d: IntSlider(value=3, description='$10^{-d}$ precision', max=16, min=1)
          maxn: IntSlider(value=10, description='max iterations', max=15)
          interval: IntRangeSlider(value=(0, 6), description='Interval', max=10, min=-10)
          list_steps: Checkbox(value=False, description='List steps')
    """
    def _newton_method(f, c, maxn, h):
        midpoints = [c]
        round = 1
        while True:
            c = c-f(c)/f.derivative(x)(x=c)
            midpoints.append(c)
            if f(c-h)*f(c+h) < 0 or round == maxn:
                break
            round += 1
        return c, midpoints

    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    a, b = interval
    h = 10**(-d)
    c, midpoints = _newton_method(f, float(c), maxn, 0.5 * h)
    html(r"$\text{Precision } 2h = %s$" % latex(float(h)))
    html(r"${c = }%s$" % c)
    html(r"${f(c) = }%s" % latex(f(c)))
    html(r"$%s \text{ iterations}" % len(midpoints))
    if list_steps:
        s = [["$n$", "$x_n$", "$f(x_n)$", r"$f(x_n-h)\,f(x_n+h)$"]]
        for i, c in enumerate(midpoints):
            s.append([i+1, c, f(c), (c-h)*f(c+h)])
        pretty_print(table(s, header_row=True))
    else:
        P = plot(f, x, interval, color='blue')
        L = sum(line([(c, 0), (c, f(c))], color='green') for c in midpoints[:-1])
        for i in range(len(midpoints) - 1):
            L += line([(midpoints[i], f(midpoints[i])), (midpoints[i+1], 0)], color='red')
        show(P + L, xmin=interval[0], xmax=interval[1], ymin=P.ymin(), ymax=P.ymax())


@library_interact(
    title=lambda: text_control("<h2>Trapezoid integration</h2>"),
    f=lambda: input_box(default="x^2-5*x + 10", label="$f(x)=$"),
    n=lambda: slider(1, 100, 1, 5, label="# divisions"),
    interval_input=lambda: selector(
        ["from slider", "from keyboard"],
        label="Integration interval",
        buttons=True,
    ),
    interval_s=lambda: range_slider(-10, 10, default=(0, 8), label="slider: "),
    interval_g=lambda: input_grid(1, 2, default=[[0, 8]], label="keyboard: "),
    output_form=lambda: selector(
        ["traditional", "table", "none"], label="Computations form", buttons=True
    ),
)
def trapezoid_integration(
    title, f, n, interval_input, interval_s, interval_g, output_form
):
    r"""
    Interact explaining the trapezoid method for definite integrals.

    Based on work by
    Lauri Ruotsalainen, 2010 (based on the application "Numerical integrals with various rules"
    by Marshall Hampton and Nick Alexander)

    INPUT:

    - ``f`` -- function of variable x to integrate
    - ``n`` -- number of divisions
    - ``interval_input`` -- switches the input for interval between slider and keyboard
    - ``interval_s`` -- slider for interval to integrate
    - ``interval_g`` -- input grid for interval to integrate
    - ``output_form`` -- the computation is formatted in a traditional form, in a table or missing

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.trapezoid_integration()
        ...Interactive function <function trapezoid_integration at ...> with 7 widgets
          title: HTMLText(value='<h2>Trapezoid integration</h2>')
          f: EvalText(value='x^2-5*x + 10', description='$f(x)=$', layout=Layout(max_width='81em'))
          n: IntSlider(value=5, description='# divisions', min=1)
          interval_input: ToggleButtons(description='Integration interval', options=('from slider', 'from keyboard'), value='from slider')
          interval_s: IntRangeSlider(value=(0, 8), description='slider: ', max=10, min=-10)
          interval_g: Grid(value=[[0, 8]], children=(Label(value='keyboard: '), VBox(children=(EvalText(value='0', layout=Layout(max_width='5em')),)), VBox(children=(EvalText(value='8', layout=Layout(max_width='5em')),))))
          output_form: ToggleButtons(description='Computations form', options=('traditional', 'table', 'none'), value='traditional')
    """
    xs = []
    ys = []
    if interval_input == 'from slider':
        interval = interval_s
    else:
        interval = interval_g[0]
    h = float(interval[1]-interval[0])/n
    x = SR.var('x')
    f = symbolic_expression(f).function(x)

    trapezoids = Graphics()

    for i in range(n):
        xi = interval[0] + i*h
        yi = f(xi)
        trapezoids += line([[xi, 0], [xi, yi], [xi + h, f(xi + h)],[xi + h, 0],[xi, 0]], rgbcolor=(1, 0, 0))
        xs.append(xi)
        ys.append(yi)
    xs.append(xi + h)
    ys.append(f(xi + h))

    html(r'Function $f(x)=%s$' % latex(f(x)))
    show(plot(f, interval[0], interval[1]) + trapezoids,
         xmin=interval[0], xmax=interval[1])

    numeric_value = integral_numerical(f, interval[0], interval[1])[0]
    approx = h * (ys[0]/2 + sum([ys[i] for i in range(1,n)]) + ys[n]/2)

    html(r'Integral value to seven decimal places is: $\displaystyle\int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x} = %.6f$' % (
            interval[0], interval[1], N(numeric_value, digits=7))
         )

    if output_form == 'traditional':
        sum_formula_html = r"\frac {d}{2} \cdot \left[f(x_0) + %s + f(x_{%s})\right]" % (
            ' + '.join([ "2 f(x_{%s})" % i for i in range(1,n)]),
            n
            )
        sum_placement_html = r"\frac{%.2f}{2} \cdot \left[f(%.2f) + %s + f(%.2f)\right]" % (
            h,
            N(xs[0], digits=5),
            ' + '.join([ "2 f(%.2f)" % N(i, digits=5) for i in xs[1:-1]]),
            N(xs[n], digits=5)
            )
        sum_values_html = r"\frac{%.2f}{2} \cdot \left[%.2f + %s + %.2f\right]" % (
            h,
            N(ys[0], digits=5),
            ' + '.join([ r"2\cdot %.2f" % N(i, digits=5) for i in ys[1:-1]]),
            N(ys[n], digits=5)
            )

        html(r'''
            <div class="math">
            \begin{align*}
            \int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x}
                & \approx %s \\
                & = %s \\
                & = %s \\
                & = %s
            \end{align*}
            </div>
        ''' % (
                interval[0], interval[1],
                sum_formula_html, sum_placement_html, sum_values_html,
                N(approx, digits=7)
        ))
    elif output_form == 'table':
        s = [['$i$', '$x_i$', '$f(x_i)$', '$m$', r'$m\cdot f(x_i)$']]
        for i in range(0,n+1):
            if i == 0 or i == n:
                j = 1
            else:
                j = 2
            s.append([i, xs[i], ys[i],j,N(j*ys[i])])
        pretty_print(table(s, header_row=True))


@library_interact(
    title=lambda: text_control("<h2>Simpson integration</h2>"),
    f=lambda: input_box(default='x*sin(x)+x+1', label="$f(x)=$"),
    n=lambda: slider(2, 100, 2, 6, label="# divisions"),
    interval_input=lambda: selector(
        ["from slider", "from keyboard"],
        label="Integration interval",
        buttons=True,
    ),
    interval_s=lambda: range_slider(-10, 10, default=(0, 10), label="slider: "),
    interval_g=lambda: input_grid(1, 2, default=[[0, 10]], label="keyboard: "),
    output_form=lambda: selector(
        ["traditional", "table", "none"], label="Computations form", buttons=True
    ),
)
def simpson_integration(
    title,
    f,
    n,
    interval_input,
    interval_s,
    interval_g,
    output_form,
):
    r"""
    Interact explaining the simpson method for definite integrals.

    Based on work by
    Lauri Ruotsalainen, 2010 (based on the application "Numerical integrals with various rules"
    by Marshall Hampton and Nick Alexander)

    INPUT:

    - ``f`` -- function of variable x to integrate
    - ``n`` -- number of divisions (mult. of 2)
    - ``interval_input`` -- switches the input for interval between slider and keyboard
    - ``interval_s`` -- slider for interval to integrate
    - ``interval_g`` -- input grid for interval to integrate
    - ``output_form`` -- the computation is formatted in a traditional form, in a table or missing

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.simpson_integration()
        ...Interactive function <function simpson_integration at ...> with 7 widgets
          title: HTMLText(value='<h2>Simpson integration</h2>')
          f: EvalText(value='x*sin(x)+x+1', description='$f(x)=$', layout=Layout(max_width='81em'))
          n: IntSlider(value=6, description='# divisions', min=2, step=2)
          interval_input: ToggleButtons(description='Integration interval', options=('from slider', 'from keyboard'), value='from slider')
          interval_s: IntRangeSlider(value=(0, 10), description='slider: ', max=10, min=-10)
          interval_g: Grid(value=[[0, 10]], children=(Label(value='keyboard: '), VBox(children=(EvalText(value='0', layout=Layout(max_width='5em')),)), VBox(children=(EvalText(value='10', layout=Layout(max_width='5em')),))))
          output_form: ToggleButtons(description='Computations form', options=('traditional', 'table', 'none'), value='traditional')
    """
    x = SR.var("x")
    f = symbolic_expression(f).function(x)
    if interval_input == 'from slider':
        interval = interval_s
    else:
        interval = interval_g[0]

    def parabola(a, b, c):
        from sage.symbolic.relation import solve
        A, B, C = SR.var("A, B, C")
        K = solve([A*a[0]**2+B*a[0]+C == a[1], A*b[0]**2+B*b[0]+C == b[1], A*c[0]**2+B*c[0]+C == c[1]], [A, B, C], solution_dict=True)[0]
        f = K[A]*x**2+K[B]*x+K[C]
        return f
    xs = []
    ys = []
    dx = float(interval[1]-interval[0])/n

    for i in range(n+1):
        xs.append(interval[0] + i*dx)
        ys.append(f(x=xs[-1]))

    parabolas = Graphics()
    lines = Graphics()

    for i in range(0, n-1, 2):
        p = parabola((xs[i],ys[i]),(xs[i+1],ys[i+1]),(xs[i+2],ys[i+2]))
        parabolas += plot(p(x=x), (x, xs[i], xs[i+2]), color='red')
        lines += line([(xs[i],ys[i]), (xs[i],0), (xs[i+2],0)],color='red')
        lines += line([(xs[i+1],ys[i+1]), (xs[i+1],0)], linestyle='-.', color='red')

    lines += line([(xs[-1],ys[-1]), (xs[-1],0)], color='red')

    html(r'Function $f(x)=%s$' % latex(f(x)))

    show(plot(f(x), x, interval[0], interval[1]) + parabolas + lines,
         xmin=interval[0], xmax=interval[1])

    numeric_value = integral_numerical(f,interval[0],interval[1])[0]
    approx = dx/3 * (ys[0] + sum([4*ys[i] for i in range(1,n,2)]) + sum([2*ys[i] for i in range(2,n,2)]) + ys[n])

    html(r'Integral value to seven decimal places is: $\displaystyle\int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x} = %.6f$' %
         (interval[0],interval[1],
          N(numeric_value,digits=7)))

    if output_form == 'traditional':
        sum_formula_html = r"\frac{d}{3} \cdot \left[ f(x_0) + %s + f(x_{%s})\right]" % (
            ' + '.join(r"%s \cdot f(x_{%s})" % (i % 2 * (-2) + 4, i + 1)
                       for i in range(0,n-1)),
            n
            )

        sum_placement_html = r"\frac{%.2f}{3} \cdot \left[ f(%.2f) +  %s + f(%.2f)\right]" % (
            dx,
            N(xs[0],digits=5),
            ' + '.join(r"%s \cdot f(%.2f)" % (i % 2 * (-2) + 4, N(xk, digits=5))
                       for i, xk in enumerate(xs[1:-1])),
            N(xs[n],digits=5)
            )

        sum_values_html = r"\frac{%.2f}{3} \cdot \left[ %s %s %s\right]" % (
            dx,
            "%.2f + " % N(ys[0],digits=5),
            ' + '.join(r"%s \cdot %.2f" % (i % 2 * (-2) + 4, N(yk, digits=5))
                       for i, yk in enumerate(ys[1:-1])),
            " + %.2f" % N(ys[n],digits=5)
            )

        html(r'''
        <div class="math">
        \begin{align*}
        \int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x}
            & \approx %s \\
            & = %s \\
            & = %s \\
            & = %.6f
        \end{align*}
        </div>
        ''' % (
                interval[0], interval[1],
                sum_formula_html, sum_placement_html, sum_values_html,
                N(approx,digits=7)
                ))
    elif output_form == 'table':
        s = [['$i$', '$x_i$', '$f(x_i)$', '$m$', r'$m\cdot f(x_i)$']]
        for i in range(n + 1):
            if i == 0 or i == n:
                j = 1
            else:
                j = (i + 1) % 2 * (-2) + 4
            s.append([i, xs[i], ys[i], j, N(j*ys[i])])
        s.append(['', '', '', r'$\sum$', '$%s$' % latex(3/dx*approx)])
        pretty_print(table(s, header_row=True))
        html(r'$\int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x}\approx\frac {%.2f}{3}\cdot %s=%s$' %
             (interval[0], interval[1],dx,latex(3/dx*approx),latex(approx)))


@library_interact(
    title=lambda: text_control("<h2>Riemann integral with random sampling</h2>"),
    f=lambda: input_box("x^2+1", label='$f(x)=$', width=40),
    n=lambda: slider(1, 30, 1, 5, label="# divisions"),
    hr1=lambda: text_control("<hr>"),
    interval_input=lambda: selector(
        ["from slider", "from keyboard"],
        label="Integration interval",
        buttons=True,
    ),
    interval_s=lambda: range_slider(-5, 10, default=(0, 2), label="slider: "),
    interval_g=lambda: input_grid(1, 2, default=[[0, 2]], label="keyboard: "),
    hr2=lambda: text_control("<hr>"),
    list_table=lambda: checkbox(default=False, label="List table"),
    auto_update=lambda: False,
)
def riemann_sum(
    title,
    f,
    n,
    hr1,
    interval_input,
    interval_s,
    interval_g,
    hr2,
    list_table,
    auto_update=False,
):
    r"""
    Interact explaining the definition of Riemann integral.

    INPUT:

    - ``f`` -- function of variable x to integrate
    - ``n`` -- number of divisions
      - ``interval_input`` -- switches the input for interval between slider and keyboard
      - ``interval_s`` -- slider for interval to integrate
      - ``interval_g`` -- input grid for interval to integrate
    - ``list_table`` -- print table with values of the function

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.riemann_sum()
        Manual interactive function <function riemann_sum at ...> with 9 widgets
          title: HTMLText(value='<h2>Riemann integral with random sampling</h2>')
          f: EvalText(value='x^2+1',... description='$f(x)=$', layout=Layout(max_width='41em'))
          n: IntSlider(value=5, description='# divisions', max=30, min=1)
          hr1: HTMLText(value='<hr>')
          interval_input: ToggleButtons(description='Integration interval', options=('from slider', 'from keyboard'), value='from slider')
          interval_s: IntRangeSlider(value=(0, 2), description='slider: ', max=10, min=-5)
          interval_g: Grid(value=[[0, 2]], children=(Label(value='keyboard: '), VBox(children=(EvalText(value='0', layout=Layout(max_width='5em')),)), VBox(children=(EvalText(value='2', layout=Layout(max_width='5em')),))))
          hr2: HTMLText(value='<hr>')
          list_table: Checkbox(value=False, description='List table')

    AUTHORS:

    - Robert Marik (2010-08)
    """
    x = SR.var('x')
    from random import random
    if interval_input == 'from slider':
        a = interval_s[0]
        b = interval_s[1]
    else:
        a = interval_g[0][0]
        b = interval_g[0][1]
    func = symbolic_expression(f).function(x)
    division = [a]+[a+random()*(b-a) for i in range(n-1)]+[b]
    division = sorted(division)
    xs = [division[i]+random()*(division[i+1]-division[i]) for i in range(n)]
    ys = [func(x_val) for x_val in xs]
    rects = Graphics()
    for i in range(n):
        body = [[division[i], 0], [division[i], ys[i]],
                [division[i+1], ys[i]], [division[i+1], 0]]
        if ys[i].n() > 0:
            color_rect = 'green'
        else:
            color_rect = 'red'
        rects = rects + polygon2d(body, rgbcolor=color_rect, alpha=0.1)\
         + point((xs[i], ys[i]), rgbcolor=(1, 0, 0))\
         + line(body, rgbcolor='black', zorder=-1)
    html('<small>Adjust your data and click Update button. Click repeatedly for another random values.</small>')

    show(plot(func(x),(x,a,b),zorder=5) + rects)
    delka_intervalu = [division[i+1]-division[i] for i in range(n)]
    if list_table:
        pretty_print(table([
            ["$i$", "$[x_{i-1},x_i]$", r"$\eta_i$", r"$f(\eta_i)$", "$x_{i}-x_{i-1}$"]
        ] + [
            [i+1,[division[i],division[i+1]],xs[i],ys[i],delka_intervalu[i]] for i in range(n)
        ], header_row=True))

    html(r'Riemann sum: $\displaystyle\sum_{i=1}^{%s} f(\eta_i)(x_i-x_{i-1})=%s$ ' %
         (latex(n),latex(sum([ys[i]*delka_intervalu[i] for i in range(n)]))))
    html(r'Exact value of the integral $\displaystyle\int_{%s}^{%s}%s\,\mathrm{d}x=%s$' %
         (latex(a),latex(b),latex(func(x)),latex(integral_numerical(func(x),a,b)[0])))


x = SR.var("x")


@library_interact(
    f=lambda: sin(x),
    g=lambda: cos(x),
    xrange=lambda: range_slider(-3, 3, default=(0, 1), label='x-range'),
    yrange=lambda: "auto",
    a=lambda: 1,
    action=lambda: selector(
        [
            "f",
            "df/dx",
            "int f",
            "num f",
            "den f",
            "1/f",
            "finv",
            "f+a",
            "f-a",
            "f*a",
            "f/a",
            "f^a",
            "f(x+a)",
            "f(x*a)",
            "f+g",
            "f-g",
            "f*g",
            "f/g",
            "f(g)",
        ],
        width=15,
        nrows=5,
        label="h = ",
    ),
    do_plot=lambda: ("Draw Plots", True),
)
def function_tool(f, g, xrange, yrange, a, action, do_plot):
    r"""
    `Function Plotting Tool <http://wiki.sagemath.org/interact/calculus#Functiontool>`_
    (by William Stein (?))

    INPUT:

    - ``f`` -- function f(x)
    - ``g`` -- function g(x)
    - ``xrange`` -- range for plotting (x)
    - ``yrange`` -- range for plotting ('auto' is default, otherwise a tuple)
    - ``a`` -- factor ``a``
    - ``action`` -- select given operation on or combination of functions
    - ``do_plot`` -- if ``True``, a plot is drawn

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.function_tool()  # long time
        ...Interactive function <function function_tool at ...> with 7 widgets
          f: EvalText(value='sin(x)', description='f')
          g: EvalText(value='cos(x)', description='g')
          xrange: IntRangeSlider(value=(0, 1), description='x-range', max=3, min=-3)
          yrange: Text(value='auto', description='yrange')
          a: IntSlider(value=1, description='a', max=3, min=-1)
          action: ToggleButtons(description='h = ', options=('f', 'df/dx', 'int f', 'num f', 'den f', '1/f', 'finv', 'f+a', 'f-a', 'f*a', 'f/a', 'f^a', 'f(x+a)', 'f(x*a)', 'f+g', 'f-g', 'f*g', 'f/g', 'f(g)'), value='f')
          do_plot: Checkbox(value=True, description='Draw Plots')
    """
    x = SR.var('x')
    try:
        f = SR(f)
        g = SR(g)
        a = SR(a)
    except TypeError as msg:
        print(msg[-200:])
        print("Unable to make sense of f,g, or a as symbolic expressions in single variable x.")
        return
    if not (isinstance(xrange, tuple) and len(xrange) == 2):
        xrange = (0, 1)
    h = 0
    lbl = ''
    if action == 'f':
        h = f
        lbl = 'f'
    elif action == 'df/dx':
        h = f.derivative(x)
        lbl = r'\frac{\mathrm{d}f}{\mathrm{d}x}'
    elif action == 'int f':
        h = f.integrate(x)
        lbl = r'\int f \,\mathrm{d}x'
    elif action == 'num f':
        h = f.numerator()
        lbl = r'\text{numer(f)}'
    elif action == 'den f':
        h = f.denominator()
        lbl = r'\text{denom(f)}'
    elif action == '1/f':
        h = 1/f
        lbl = r'\frac{1}{f}'
    elif action == 'finv':
        h = solve(f == SR.var('y'), x)[0].rhs()
        lbl = 'f^{-1}(y)'
    elif action == 'f+a':
        h = f+a
        lbl = 'f + a'
    elif action == 'f-a':
        h = f-a
        lbl = 'f - a'
    elif action == 'f*a':
        h = f*a
        lbl = r'f \times a'
    elif action == 'f/a':
        h = f/a
        lbl = r'\frac{f}{a}'
    elif action == 'f^a':
        h = f**a
        lbl = 'f^a'
    elif action == 'f^a':
        h = f**a
        lbl = 'f^a'
    elif action == 'f(x+a)':
        h = f.subs(x=x+a)
        lbl = 'f(x+a)'
    elif action == 'f(x*a)':
        h = f.subs(x=x*a)
        lbl = 'f(xa)'
    elif action == 'f+g':
        h = f+g
        lbl = 'f + g'
    elif action == 'f-g':
        h = f-g
        lbl = 'f - g'
    elif action == 'f*g':
        h = f*g
        lbl = r'f \times g'
    elif action == 'f/g':
        h = f/g
        lbl = r'\frac{f}{g}'
    elif action == 'f(g)':
        h = f(g)
        lbl = 'f(g)'
    html('<center><font color="red">$f = %s$</font></center>' % latex(f))
    html('<center><font color="green">$g = %s$</font></center>' % latex(g))
    html('<center><font color="blue"><b>$h = %s = %s$</b></font></center>' % (lbl, latex(h)))
    if do_plot:
        P = plot(f, xrange, color='red', thickness=2) +  \
            plot(g, xrange, color='green', thickness=2) + \
            plot(h, xrange, color='blue', thickness=2)
        if yrange == 'auto':
            show(P, xmin=xrange[0], xmax=xrange[1])
        else:
            yrange = sage_eval(yrange)
            show(P, xmin=xrange[0], xmax=xrange[1], ymin=yrange[0], ymax=yrange[1])


@library_interact(
    expo=lambda: slider(-10, 10, 0.1, 2),
    c_real=lambda: slider(-2, 2, 0.01, 0.5, label="real part const."),
    c_imag=lambda: slider(-2, 2, 0.01, 0.5, label="imag part const."),
    iterations=lambda: slider(1, 100, 1, 20, label="# iterations"),
    zoom_x=lambda: range_slider(-2, 2, 0.01, (-1.5, 1.5), label="Zoom X"),
    zoom_y=lambda: range_slider(-2, 2, 0.01, (-1.5, 1.5), label="Zoom Y"),
    plot_points=lambda: slider(20, 400, 20, default=150, label="plot points"),
    dpi=lambda: slider(20, 200, 10, default=80, label='dpi'),
)
def julia(expo, c_real, c_imag, iterations, zoom_x, zoom_y, plot_points, dpi):
    r"""
    Julia Fractal, based on
    `Julia by Harald Schilly <http://wiki.sagemath.org/interact/fractal#Julia>`_.

    INPUT:

    - ``exponent`` -- exponent ``e`` in `z^e+c`
    - ``c_real`` -- real part of the constant ``c``
    - ``c_imag`` -- imaginary part of the constant ``c``
    - ``iterations`` -- number of iterations
    - ``zoom_x`` -- range slider for zoom in x direction
    - ``zoom_y`` -- range slider for zoom in y direction
    - ``plot_points`` -- number of points to plot
    - ``dpi`` -- dots-per-inch parameter for the plot

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.fractals.julia()  # long time
        ...Interactive function <function julia at ...> with 8 widgets
          expo: FloatSlider(value=2.0, description='expo', max=10.0, min=-10.0)
          c_real: FloatSlider(value=0.5, description='real part const.', max=2.0, min=-2.0, step=0.01)
          c_imag: FloatSlider(value=0.5, description='imag part const.', max=2.0, min=-2.0, step=0.01)
          iterations: IntSlider(value=20, description='# iterations', min=1)
          zoom_x: FloatRangeSlider(value=(-1.5, 1.5), description='Zoom X', max=2.0, min=-2.0, step=0.01)
          zoom_y: FloatRangeSlider(value=(-1.5, 1.5), description='Zoom Y', max=2.0, min=-2.0, step=0.01)
          plot_points: IntSlider(value=150, description='plot points', max=400, min=20, step=20)
          dpi: IntSlider(value=80, description='dpi', max=200, min=20, step=10)
    """
    z = SR.var('z')
    I = CDF.gen()
    f = symbolic_expression(z**expo + c_real + c_imag*I).function(z)
    ff_j = fast_callable(f, vars=[z], domain=CDF)

    from sage.interacts.library_cython import julia

    html('<h2>Julia Fractal</h2>')
    html(r'Recursive Formula: $z \leftarrow z^{%.2f} + (%.2f+%.2f*\mathbb{I})$' % (expo, c_real, c_imag))
    complex_plot(lambda z: julia(ff_j, z, iterations), zoom_x, zoom_y, plot_points=plot_points, dpi=dpi).show(frame=True, aspect_ratio=1)


@library_interact(
    expo=lambda: slider(-10, 10, 0.1, 2),
    iterations=lambda: slider(1, 100, 1, 20, label="# iterations"),
    zoom_x=lambda: range_slider(-2, 2, 0.01, (-2, 1), label="Zoom X"),
    zoom_y=lambda: range_slider(-2, 2, 0.01, (-1.5, 1.5), label="Zoom Y"),
    plot_points=lambda: slider(20, 400, 20, default=150, label="plot points"),
    dpi=lambda: slider(20, 200, 10, default=80, label='dpi'),
)
def mandelbrot(expo, iterations, zoom_x, zoom_y, plot_points, dpi):
    r"""
    Mandelbrot Fractal, based on
    `Mandelbrot by Harald Schilly <http://wiki.sagemath.org/interact/fractal#Mandelbrot>`_.

    INPUT:

    - ``exponent`` -- exponent ``e`` in `z^e+c`
    - ``iterations`` -- number of iterations
    - ``zoom_x`` -- range slider for zoom in x direction
    - ``zoom_y`` -- range slider for zoom in y direction
    - ``plot_points`` -- number of points to plot
    - ``dpi`` -- dots-per-inch parameter for the plot

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.fractals.mandelbrot()  # long time
        ...Interactive function <function mandelbrot at ...> with 6 widgets
          expo: FloatSlider(value=2.0, description='expo', max=10.0, min=-10.0)
          iterations: IntSlider(value=20, description='# iterations', min=1)
          zoom_x: FloatRangeSlider(value=(-2.0, 1.0), description='Zoom X', max=2.0, min=-2.0, step=0.01)
          zoom_y: FloatRangeSlider(value=(-1.5, 1.5), description='Zoom Y', max=2.0, min=-2.0, step=0.01)
          plot_points: IntSlider(value=150, description='plot points', max=400, min=20, step=20)
          dpi: IntSlider(value=80, description='dpi', max=200, min=20, step=10)
    """
    x, z, c = SR.var('x, z, c')
    f = symbolic_expression(z**expo + c).function(z, c)
    ff_m = fast_callable(f, vars=[z,c], domain=CDF)

    from sage.interacts.library_cython import mandel

    html('<h2>Mandelbrot Fractal</h2>')
    html(r'Recursive Formula: $z \leftarrow z^{%.2f} + c$ for $c \in \mathbb{C}$' % expo)
    complex_plot(lambda z: mandel(ff_m, z, iterations), zoom_x, zoom_y, plot_points=plot_points, dpi=dpi).show(frame=True, aspect_ratio=1)


@library_interact(
    N=lambda: slider(1, 500, 1, label="Number of iterations", default=100),
    rule_number=lambda: slider(0, 255, 1, default=110, label="Rule number"),
    size=lambda: slider(1, 11, step_size=1, default=6, label="size of graphic"),
)
def cellular_automaton(N, rule_number, size):
    r"""
    Yields a matrix showing the evolution of a
    `Wolfram's cellular automaton <http://mathworld.wolfram.com/CellularAutomaton.html>`_.

    `Based on work by Pablo Angulo <http://wiki.sagemath.org/interact/misc#CellularAutomata>`_.

    INPUT:

    - ``N`` -- iterations
    - ``rule_number`` -- rule number (0 to 255)
    - ``size`` -- size of the shown picture

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.fractals.cellular_automaton()  # long time
        ...Interactive function <function cellular_automaton at ...> with 3 widgets
          N: IntSlider(value=100, description='Number of iterations', max=500, min=1)
          rule_number: IntSlider(value=110, description='Rule number', max=255)
          size: IntSlider(value=6, description='size of graphic', max=11, min=1)
    """
    from sage.rings.integer import Integer
    if not 0 <= rule_number <= 255:
        raise ValueError('Invalid rule number')
    binary_digits = Integer(rule_number).digits(base=2)
    rule = binary_digits + [0]*(8-len(binary_digits))

    html('<h2>Cellular Automaton</h2>' +
         '<div style="white-space: normal;">"A cellular automaton is a collection of "colored" cells \
         on a grid of specified shape that evolves through a number of \
         discrete time steps according to a set of rules based on the \
         states of neighboring cells." &mdash; \
         <a target="_blank" href="http://mathworld.wolfram.com/CellularAutomaton.html">Mathworld,\
         Cellular Automaton</a></div>\
         <div>Rule %s expands to %s</div>' % (rule_number, ''.join(map(str,rule)))
        )

    from sage.interacts.library_cython import cellular
    M = cellular(rule, N)
    plot_M = matrix_plot(M, cmap='binary')
    plot_M.show(figsize=[size,size])


@library_interact(
    interval=lambda: range_slider(1, 4000, 10, default=(1, 1000), label='range'),
    show_factors=lambda: True,
    highlight_primes=lambda: True,
    show_curves=lambda: True,
    n=lambda: slider(1, 200, 1, default=89, label="number $n$"),
    dpi=lambda: slider(10, 300, 10, default=100, label='dpi'),
)
def polar_prime_spiral(interval, show_factors, highlight_primes, show_curves, n, dpi):
    r"""
    Polar Prime Spiral interact, based on work by David Runde.

    For more information about the factors in the spiral,
    `visit John Williamson's website <http://www.dcs.gla.ac.uk/~jhw/spirals/index.html>`_.

    INPUT:

    - ``interval`` -- range slider to specify start and end
    - ``show_factors`` -- if ``True``, show factors
    - ``highlight_primes`` -- if ``True``, prime numbers are highlighted
    - ``show_curves`` -- if ``True``, curves are plotted
    - ``n`` -- number `n`
    - ``dpi`` -- dots per inch resolution for plotting

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: sage.interacts.algebra.polar_prime_spiral()  # long time
        ...Interactive function <function polar_prime_spiral at ...> with 6 widgets
          interval: IntRangeSlider(value=(1, 1000), description='range', max=4000, min=1, step=10)
          show_factors: Checkbox(value=True, description='show_factors')
          highlight_primes: Checkbox(value=True, description='highlight_primes')
          show_curves: Checkbox(value=True, description='show_curves')
          n: IntSlider(value=89, description='number $n$', max=200, min=1)
          dpi: IntSlider(value=100, description='dpi', max=300, min=10, step=10)
    """
    html('<h2>Polar Prime Spiral</h2> \
          <div style="white-space: normal;">\
          For more information about the factors in the spiral, visit \
          <a href="http://www.dcs.gla.ac.uk/~jhw/spirals/index.html" target="_blank">\
          Number Spirals by John Williamson</a>.</div>'
        )

    start, end = interval
    from sage.ext.fast_eval import fast_float
    from math import ceil
    from sage.plot.colors import hue

    if start < 1 or end <= start:
        print("invalid start or end value")
        return
    if n > end:
        print("WARNING: n is greater than end value")
        return
    if n < start:
        print("n < start value")
        return
    nn = SR.var('nn')
    f1 = fast_float(sqrt(nn)*cos(2*pi*sqrt(nn)), 'nn')
    f2 = fast_float(sqrt(nn)*sin(2*pi*sqrt(nn)), 'nn')
    f = lambda x: (f1(x), f2(x))

    list = []
    list2 = []
    if not show_factors:
        for i in srange(start, end, include_endpoint=True):
            if Integer(i).is_pseudoprime():
                list.append(f(i-start+1))  # primes list
            else:
                list2.append(f(i-start+1))  # composites list
        P = points(list)
        R = points(list2, alpha=.1)  # faded composites
    else:
        for i in srange(start, end, include_endpoint=True):
            # Resize each of the dots depending of the number of factors of each number
            list.append(disk((f(i-start+1)),0.05*pow(2,len(factor(i))-1), (0,2*pi)))
            if Integer(i).is_pseudoprime() and highlight_primes:
                list2.append(f(i-start+1))
        P = Graphics()
        for g in list:
            P += g
        p_size = 5  # the orange dot size of the prime markers
        if not highlight_primes:
            list2 = [(f(n-start+1))]
        R = points(list2, hue=.1, pointsize=p_size)

    if n > 0:
        html('$n = %s$' % factor(n))

        p = 1
        # The X which marks the given n
        W1 = disk((f(n-start+1)), p, (pi/6, 2*pi/6), alpha=.1)
        W2 = disk((f(n-start+1)), p, (4*pi/6, 5*pi/6), alpha=.1)
        W3 = disk((f(n-start+1)), p, (7*pi/6, 8*pi/6), alpha=.1)
        W4 = disk((f(n-start+1)), p, (10*pi/6, 11*pi/6), alpha=.1)
        Q = W1 + W2 + W3 + W4

        n -= start - 1  # offset n for different start values to ensure accurate plotting
        if show_curves:
            begin_curve = 0
            t = SR.var('t')
            a = 1.0
            b = 0.0
            S = int(sqrt(n))
            if n <= S * (S + 1):
                c = n - S**2
            else:
                c = n - (S + 1)**2
            c2 = n - S * (S + 1)
            html('Pink Curve:  $n^2 + %s$' % c)
            html('Green Curve: $n^2 + n + %s$' % c2)
            m = SR.var('m')
            g = symbolic_expression(a*m**2+b*m+c).function(m)
            r = symbolic_expression(sqrt(g(m))).function(m)
            theta = symbolic_expression(r(m) - m*sqrt(a)).function(m)
            S1 = parametric_plot(((r(t))*cos(2*pi*(theta(t))),(r(t))*sin(2*pi*(theta(t)))),
                                 (begin_curve, ceil(sqrt(end-start))),
                                 color=hue(0.8), thickness=.3)  # pink line

            b = 1
            c = c2
            g = symbolic_expression(a*m**2+b*m+c).function(m)
            r = symbolic_expression(sqrt(g(m))).function(m)
            theta = symbolic_expression(r(m) - m*sqrt(a)).function(m)
            S2 = parametric_plot(((r(t))*cos(2*pi*(theta(t))),(r(t))*sin(2*pi*(theta(t)))),
                                 (begin_curve, ceil(sqrt(end-start))),
                                 color=hue(0.6), thickness=.3)  # green line

            show(R+P+S1+S2+Q, aspect_ratio=1, axes=False, dpi=dpi)
        else:
            show(R+P+Q, aspect_ratio=1, axes=False, dpi=dpi)
    else:
        show(R+P, aspect_ratio=1, axes=False, dpi=dpi)
