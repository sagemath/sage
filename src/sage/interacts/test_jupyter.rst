.. -*- coding: utf-8 -*-

.. linkall

Test the Sage interact library in Jupyter

We need to setup a proper test environment for widgets::

    sage: from ipywidgets.widgets.tests.utils import setup_test_comm
    sage: setup_test_comm()

Make sure that we use the Jupyter interacts::

    sage: from sage.repl.ipython_kernel.all_jupyter import *

Function to test library interacts. This will "display" the widgets
textually and it will actually run the interactive function. Graphics
which are produced by the interact are not displayed, but text is. ::

    sage: def test(lib_interact):
    ....:     f = lib_interact.f      # Underlying non-wrapped function
    ....:     widgets = lib_interact._widgets
    ....:     widgets_constructed = {key: constructor() for key, constructor in widgets.items()}
    ....:     w = interact(f, **widgets_constructed).widget  # "interactive" widget
    ....:     kwargs = {widget._kwarg: widget.get_interact_value()
    ....:             for widget in w.kwargs_widgets}
    ....:     return f(**kwargs)

This is just to test that failures in the interact are actually seen::

    sage: from sage.interacts.library import library_interact
    sage: @library_interact()
    ....: def failure():
    ....:     raise Exception("gotcha")
    sage: test(failure)
    Traceback (most recent call last):
    ...
    Exception: gotcha

Test all interacts from the Sage interact library::

    sage: test(interacts.algebra.polar_prime_spiral)  # long time
    ...Interactive function <function polar_prime_spiral at ...> with 6 widgets
      interval: IntRangeSlider(value=(1, 1000), description='range', max=4000, min=1, step=10)
      show_factors: Checkbox(value=True, description='show_factors')
      highlight_primes: Checkbox(value=True, description='highlight_primes')
      show_curves: Checkbox(value=True, description='show_curves')
      n: IntSlider(value=89, description='number $n$', max=200, min=1)
      dpi: IntSlider(value=100, description='dpi', max=300, min=10, step=10)
    <h2>Polar Prime Spiral</h2>           <div style="white-space: normal;">          For more information about the factors in the spiral, visit           <a href="http://www.dcs.gla.ac.uk/~jhw/spirals/index.html" target="_blank">          Number Spirals by John Williamson</a>.</div>
    \(n = 89\)
    Pink Curve: \(n^2 + 8\)
    Green Curve: \(n^2 + n + -1\)

    sage: test(interacts.calculus.taylor_polynomial)
    ...Interactive function <function taylor_polynomial at ...> with 3 widgets
      title: HTMLText(value='<h2>Taylor polynomial</h2>')
      f: EvalText(value='e^(-x)*sin(x)', description='$f(x)=$', layout=Layout(max_width='81em'))
      order: SelectionSlider(description='order', options=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), value=1)
    \(f(x)\;=\;e^{\left(-x\right)} \sin\left(x\right)\)
    \(\hat{f}(x;0)\;=\;x+\mathcal{O}(x^{2})\)

    sage: test(interacts.calculus.definite_integral)
    ...Interactive function <function definite_integral at ...> with 6 widgets
      title: HTMLText(value='<h2>Definite integral</h2>')
      f: EvalText(value='3*x', description='$f(x)=$', layout=Layout(max_width='81em'))
      g: EvalText(value='x^2', description='$g(x)=$', layout=Layout(max_width='81em'))
      interval: IntRangeSlider(value=(0, 3), description='Interval', max=10, min=-10)
      x_range: IntRangeSlider(value=(0, 3), description='plot range (x)', max=10, min=-10)
      selection: Dropdown(description='Select', index=2, options=('f', 'g', 'f and g', 'f - g'), value='f and g')
    \(\int_{0.00}^{3.00}(\color{Blue}{f(x)})\,\mathrm{d}x=\int_{0.00}^{3.00}(3 \, x)\,\mathrm{d}x=13.50\)<br/>\(\int_{0.00}^{3.00}(\color{Green}{g(x)})\,\mathrm{d}x=\int_{0.00}^{3.00}(x^{2})\,\mathrm{d}x=9.00\)

    sage: test(interacts.calculus.function_derivative)
    ...Interactive function <function function_derivative at ...> with 4 widgets
      title: HTMLText(value='<h2>Derivative grapher</h2>')
      function: EvalText(value='x^5-3*x^3+1', description='Function:', layout=Layout(max_width='81em'))
      x_range: FloatRangeSlider(value=(-2.0, 2.0), description='Range (x)', max=15.0, min=-15.0)
      y_range: FloatRangeSlider(value=(-8.0, 6.0), description='Range (y)', max=15.0, min=-15.0)
    <center>\(\color{Blue}{f(x) = x^{5} - 3 \, x^{3} + 1}\)</center>
    <center>\(\color{Green}{f'(x) = 5 \, x^{4} - 9 \, x^{2}}\)</center>
    <center>\(\color{Red}{f''(x) = 20 \, x^{3} - 18 \, x}\)</center>

    sage: test(interacts.calculus.difference_quotient)
    ...Interactive function <function difference_quotient at ...> with 5 widgets
      title: HTMLText(value='<h2>Difference quotient</h2>')
      f: EvalText(value='sin(x)', description='f(x)', layout=Layout(max_width='81em'))
      interval: FloatRangeSlider(value=(0.0, 10.0), description='Range', max=10.0)
      a: IntSlider(value=5, description='$a$', max=10)
      x0: IntSlider(value=2, description='$x_0$ (start point)', max=10)
    <h2>Difference Quotient</h2>
    <div style="white-space: normal;">         <a href="https://en.wikipedia.org/wiki/Difference_quotient" target="_blank">         Wikipedia article about difference quotient</a></div>
    <h2>Difference Quotient</h2>
    <br>\(\text{Line's equation:}\)
    \(y = 1/3*(x - 5)*(sin(5) - sin(2)) + sin(5)\)<br>
    \(\text{Slope:}\)
    \(k = \frac{f(x_0)-f(a)}{x_0-a} = -0.62274\)<br>

    sage: test(interacts.calculus.quadratic_equation)
    ...Interactive function <function quadratic_equation at ...> with 3 widgets
      A: IntSlider(value=1, description='A', max=7, min=-7)
      B: IntSlider(value=1, description='B', max=7, min=-7)
      C: IntSlider(value=-2, description='C', max=7, min=-7)
    <h2>The Solutions of the Quadratic Equation</h2>
    \(x^2 + x - 2 = 0\)
    \(Ax^2 + Bx + C = 0\)
    \(x = \frac{-B\pm\sqrt{B^2-4AC}}{2A} = \frac{-1\pm\sqrt{1^2-4*1*-2}}{2*1} = \frac{-1\pm\sqrt{\color{Green}{9}}}{2} = \begin{cases}1\\-2\end{cases}\)

    sage: test(interacts.calculus.secant_method)
    ...Interactive function <function secant_method at ...> with 5 widgets
      title: HTMLText(value='<h2>Secant method for numerical root finding</h2>')
      f: EvalText(value='x^2-2', description='f(x)', layout=Layout(max_width='81em'))
      interval: IntRangeSlider(value=(0, 4), description='range', max=5, min=-5)
      d: IntSlider(value=3, description='10^-d precision', max=16, min=1)
      maxn: IntSlider(value=10, description='max iterations', max=15)
    \(\text{Precision }h = 10^{-d}=10^{-3}=0.00100\)
    \({c = }1.4144038097709382\)
    \({f(c) = }0.0005381370945443109\)
    \(6 \text{ iterations}\)

    sage: test(interacts.calculus.newton_method)
    ...Interactive function <function newton_method at ...> with 7 widgets
      title: HTMLText(value='<h2>Newton method</h2>')
      f: EvalText(value='x^2 - 2', description='f', layout=Layout(max_width='81em'))
      c: IntSlider(value=6, description='Start ($x$)', max=10, min=-10)
      d: IntSlider(value=3, description='$10^{-d}$ precision', max=16, min=1)
      maxn: IntSlider(value=10, description='max iterations', max=15)
      interval: IntRangeSlider(value=(0, 6), description='Interval', max=10, min=-10)
      list_steps: Checkbox(value=False, description='List steps')
    \(\text{Precision } 2h = 0.001\)
    \({c = }1.4142141576301823\)
    \({f(c) = }1.6836416460996873 \times 10^{-06}\)
    \(6 \text{ iterations}\)

    sage: test(interacts.calculus.trapezoid_integration)
    ...Interactive function <function trapezoid_integration at ...> with 7 widgets
      title: HTMLText(value='<h2>Trapezoid integration</h2>')
      f: EvalText(value='x^2-5*x + 10', description='$f(x)=$', layout=Layout(max_width='81em'))
      n: IntSlider(value=5, description='# divisions', min=1)
      interval_input: ToggleButtons(description='Integration interval', options=('from slider', 'from keyboard'), value='from slider')
      interval_s: IntRangeSlider(value=(0, 8), description='slider: ', max=10, min=-10)
      interval_g: Grid(value=[[0, 8]], children=(Label(value='keyboard: '), VBox(children=(EvalText(value='0', layout=Layout(max_width='5em')),)), VBox(children=(EvalText(value='8', layout=Layout(max_width='5em')),))))
      output_form: ToggleButtons(description='Computations form', options=('traditional', 'table', 'none'), value='traditional')
    Function \(f(x)=x^{2} - 5 \, x + 10\)
    Integral value to seven decimal places is: \(\displaystyle\int_{0.00}^{8.00} {f(x) \, \mathrm{d}x} = 90.666667\)
    <BLANKLINE>
                <div class="math">
                \begin{align*}
                \int_{0.00}^{8.00} {f(x) \, \mathrm{d}x}
                    & \approx \frac {d}{2} \cdot \left[f(x_0) + 2 f(x_{1}) + 2 f(x_{2}) + 2 f(x_{3}) + 2 f(x_{4}) + f(x_{5})\right] \\
                    & = \frac{1.60}{2} \cdot \left[f(0.00) + 2 f(1.60) + 2 f(3.20) + 2 f(4.80) + 2 f(6.40) + f(8.00)\right] \\
                    & = \frac{1.60}{2} \cdot \left[10.00 + 2\cdot 4.56 + 2\cdot 4.24 + 2\cdot 9.04 + 2\cdot 18.96 + 34.00\right] \\
                    & = 94.08000
                \end{align*}
                </div>
    <BLANKLINE>

    sage: test(interacts.calculus.simpson_integration)
    ...Interactive function <function simpson_integration at ...> with 7 widgets
      title: HTMLText(value='<h2>Simpson integration</h2>')
      f: EvalText(value='x*sin(x)+x+1', description='$f(x)=$', layout=Layout(max_width='81em'))
      n: IntSlider(value=6, description='# divisions', min=2, step=2)
      interval_input: ToggleButtons(description='Integration interval', options=('from slider', 'from keyboard'), value='from slider')
      interval_s: IntRangeSlider(value=(0, 10), description='slider: ', max=10, min=-10)
      interval_g: Grid(value=[[0, 10]], children=(Label(value='keyboard: '), VBox(children=(EvalText(value='0', layout=Layout(max_width='5em')),)), VBox(children=(EvalText(value='10', layout=Layout(max_width='5em')),))))
      output_form: ToggleButtons(description='Computations form', options=('traditional', 'table', 'none'), value='traditional')
    Function \(f(x)=x \sin\left(x\right) + x + 1\)
    Integral value to seven decimal places is: \(\displaystyle\int_{0.00}^{10.00} {f(x) \, \mathrm{d}x} = 67.846694\)
    <BLANKLINE>
            <div class="math">
            \begin{align*}
            \int_{0.00}^{10.00} {f(x) \, \mathrm{d}x}
                & \approx \frac{d}{3} \cdot \left[ f(x_0) + 4 \cdot f(x_{1}) + 2 \cdot f(x_{2}) + 4 \cdot f(x_{3}) + 2 \cdot f(x_{4}) + 4 \cdot f(x_{5}) + f(x_{6})\right] \\
                & = \frac{1.67}{3} \cdot \left[ f(0.00) +  4 \cdot f(1.67) + 2 \cdot f(3.33) + 4 \cdot f(5.00) + 2 \cdot f(6.67) + 4 \cdot f(8.33) + f(10.00)\right] \\
                & = \frac{1.67}{3} \cdot \left[ 1.00 +  4 \cdot 4.33 + 2 \cdot 3.70 + 4 \cdot 1.21 + 2 \cdot 10.16 + 4 \cdot 16.73  + 5.56\right] \\
                & = 68.506699
            \end{align*}
            </div>
    <BLANKLINE>

    sage: test(interacts.calculus.bisection_method)
    ...Interactive function <function bisection_method at ...> with 5 widgets
      title: HTMLText(value='<h2>Bisection method</h2>')
      f: EvalText(value='x^2-2', description='f(x)', layout=Layout(max_width='81em'))
      interval: IntRangeSlider(value=(0, 4), description='range', max=5, min=-5)
      d: IntSlider(value=3, description='$10^{-d}$ precision', max=8, min=1)
      maxn: IntSlider(value=10, description='max iterations', max=50)
    \(\text{Precision }h = 10^{-d}=10^{-3}=0.00100\)
    \({c = }1.4140625\)
    \({f(c) = }-0.00042724609375\)
    \(9 \text{ iterations}\)

    sage: test(interacts.calculus.riemann_sum)
    ...Manual interactive function <function riemann_sum at ...> with 9 widgets
      title: HTMLText(value='<h2>Riemann integral with random sampling</h2>')
      f: EvalText(value='x^2+1',... description='$f(x)=$', layout=Layout(max_width='41em'))
      n: IntSlider(value=5, description='# divisions', max=30, min=1)
      hr1: HTMLText(value='<hr>')
      interval_input: ToggleButtons(description='Integration interval', options=('from slider', 'from keyboard'), value='from slider')
      interval_s: IntRangeSlider(value=(0, 2), description='slider: ', max=10, min=-5)
      interval_g: Grid(value=[[0, 2]], children=(Label(value='keyboard: '), VBox(children=(EvalText(value='0', layout=Layout(max_width='5em')),)), VBox(children=(EvalText(value='2', layout=Layout(max_width='5em')),))))
      hr2: HTMLText(value='<hr>')
      list_table: Checkbox(value=False, description='List table')
    <small>Adjust your data and click Update button. Click repeatedly for another random values.</small>
    Riemann sum: \(\displaystyle\sum_{i=1}^{5} f(\eta_i)(x_i-x_{i-1})=...\)
    Exact value of the integral \(\displaystyle\int_{0}^{2}x^{2} +
    1\,\mathrm{d}x=4.666666666666668\)

    sage: test(interacts.calculus.function_tool)
    ...Interactive function <function function_tool at ...> with 7 widgets
      f: EvalText(value='sin(x)', description='f')
      g: EvalText(value='cos(x)', description='g')
      xrange: IntRangeSlider(value=(0, 1), description='x-range', max=3, min=-3)
      yrange: Text(value='auto', description='yrange')
      a: IntSlider(value=1, description='a', max=3, min=-1)
      action: ToggleButtons(description='h = ', options=('f', 'df/dx', 'int f', 'num f', 'den f', '1/f', 'finv', 'f+a', 'f-a', 'f*a', 'f/a', 'f^a', 'f(x+a)', 'f(x*a)', 'f+g', 'f-g', 'f*g', 'f/g', 'f(g)'), value='f')
      do_plot: Checkbox(value=True, description='Draw Plots')
    <center><font color="red">\(f = \sin\left(x\right)\)</font></center>
    <center><font color="green">\(g = \cos\left(x\right)\)</font></center>
    <center><font color="blue"><b>\(h = f = \sin\left(x\right)\)</b></font></center>

    sage: test(interacts.fractals.mandelbrot)
    ...Interactive function <function mandelbrot at ...> with 6 widgets
      expo: FloatSlider(value=2.0, description='expo', max=10.0, min=-10.0)
      iterations: IntSlider(value=20, description='# iterations', min=1)
      zoom_x: FloatRangeSlider(value=(-2.0, 1.0), description='Zoom X', max=2.0, min=-2.0, step=0.01)
      zoom_y: FloatRangeSlider(value=(-1.5, 1.5), description='Zoom Y', max=2.0, min=-2.0, step=0.01)
      plot_points: IntSlider(value=150, description='plot points', max=400, min=20, step=20)
      dpi: IntSlider(value=80, description='dpi', max=200, min=20, step=10)
    <h2>Mandelbrot Fractal</h2>
    Recursive Formula: \(z \leftarrow z^{2.00} + c\) for \(c \in \mathbb{C}\)

    sage: test(interacts.fractals.julia)
    ...Interactive function <function julia at ...> with 8 widgets
      expo: FloatSlider(value=2.0, description='expo', max=10.0, min=-10.0)
      c_real: FloatSlider(value=0.5, description='real part const.', max=2.0, min=-2.0, step=0.01)
      c_imag: FloatSlider(value=0.5, description='imag part const.', max=2.0, min=-2.0, step=0.01)
      iterations: IntSlider(value=20, description='# iterations', min=1)
      zoom_x: FloatRangeSlider(value=(-1.5, 1.5), description='Zoom X', max=2.0, min=-2.0, step=0.01)
      zoom_y: FloatRangeSlider(value=(-1.5, 1.5), description='Zoom Y', max=2.0, min=-2.0, step=0.01)
      plot_points: IntSlider(value=150, description='plot points', max=400, min=20, step=20)
      dpi: IntSlider(value=80, description='dpi', max=200, min=20, step=10)
    <h2>Julia Fractal</h2>
    Recursive Formula: \(z \leftarrow z^{2.00} + (0.50+0.50*\mathbb{I})\)

    sage: test(interacts.fractals.cellular_automaton)
    ...Interactive function <function cellular_automaton at ...> with 3 widgets
      N: IntSlider(value=100, description='Number of iterations', max=500, min=1)
      rule_number: IntSlider(value=110, description='Rule number', max=255)
      size: IntSlider(value=6, description='size of graphic', max=11, min=1)
    <h2>Cellular Automaton</h2><div style="white-space: normal;">"A cellular automaton is a collection of "colored" cells          on a grid of specified shape that evolves through a number of          discrete time steps according to a set of rules based on the          states of neighboring cells." &mdash;          <a target="_blank" href="http://mathworld.wolfram.com/CellularAutomaton.html">Mathworld,         Cellular Automaton</a></div>         <div>Rule 110 expands to 01110110</div>

    sage: test(interacts.geometry.unit_circle)
    ...Interactive function <function unit_circle at ...> with 2 widgets
      function: Dropdown(description='function', options=(('sin(x)', 0), ('cos(x)', 1), ('tan(x)', 2)), value=0)
      x: TransformFloatSlider(value=0.0, description='x', max=6.283185307179586, step=0.015707963267948967)
    <div style="white-space: normal;">Lines of the same color have         the same length</div>

    sage: test(interacts.geometry.trigonometric_properties_triangle)
    ...Interactive function <function trigonometric_properties_triangle at ...> with 3 widgets
      a0: IntSlider(value=30, description='A', max=360)
      a1: IntSlider(value=180, description='B', max=360)
      a2: IntSlider(value=300, description='C', max=360)
    <h2>Trigonometric Properties of a Triangle</h2>
    \(\angle A = {60.000}^{\circ},\) \(\angle B = {45.000}^{\circ},\) \(\angle C = {75.000}^{\circ}\)
    \(AB = 1.931852\), \(BC = 1.732051\), \(CA = 1.414214\)
    Area of triangle \(ABC = 1.183013\)

    sage: test(interacts.geometry.special_points)
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

    sage: test(interacts.statistics.coin)
    ...Interactive function <function coin at ...> with 2 widgets
      n: IntSlider(value=1000, description='Number of Tosses', max=10000, min=2, step=100)
      interval: FloatRangeSlider(value=(0.45, 0.55), description='Plotting range (y)', max=1.0)

Test matrix control (see :issue:`27735`)::

    sage: @library_interact(A=lambda: matrix(QQ, 2, 2, range(4)))
    ....: def matrix_test(A):
    ....:     print(A)
    ....:     print(parent(A))
    sage: test(matrix_test)
    ...Interactive function <function matrix_test at ...> with 1 widget
      A: Grid(value=[[0, 1], [2, 3]], children=(Label(value='A'), VBox(children=(EvalText(value='0', layout=Layout(max_width='5em')), EvalText(value='2', layout=Layout(max_width='5em')))), VBox(children=(EvalText(value='1', layout=Layout(max_width='5em')), EvalText(value='3', layout=Layout(max_width='5em'))))))
    [0 1]
    [2 3]
    Full MatrixSpace of 2 by 2 dense matrices over Rational Field
