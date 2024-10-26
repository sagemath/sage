r"""
Mandelbrot and Julia sets

Plots the Mandelbrot and Julia sets for general polynomial maps in the complex
plane.

The Mandelbrot set is the set of complex numbers `c` for which the map
`f_c(z)` does not diverge when iterated from `z = 0`. This set of complex
numbers can be visualized by plotting each value for `c` in the complex plane.
The Mandelbrot set is often an example of a fractal when plotted in the complex
plane. For general one parameter families of polynomials, the mandelbrot set
is the parameter values for which the orbits of all critical points remains
bounded.

The Julia set for a given parameter `c` is the set of complex numbers for which
the function `f_c(z)` is bounded under iteration.

AUTHORS:

- Ben Barros
"""

#*****************************************************************************
#       Copyright (C) 2017 BEN BARROS <bbarros@slu.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.dynamics.complex_dynamics.mandel_julia_helper import (fast_mandelbrot_plot,
                                                                fast_external_ray,
                                                                convert_to_pixels,
                                                                get_line,
                                                                fast_julia_plot,
                                                                general_julia,
                                                                polynomial_mandelbrot,
                                                                julia_helper)

from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.colors", "Color")
from sage.repl.image import Image
from sage.functions.log import logb
from sage.functions.other import floor
from sage.rings.rational_field import QQ
from sage.rings.cc import CC
from sage.rings.complex_double import CDF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.misc.prandom import randint
from sage.calculus.var import var
from sage.rings.fraction_field import FractionField_generic
from sage.categories.function_fields import FunctionFields

lazy_import('sage.dynamics.arithmetic_dynamics.generic_ds', 'DynamicalSystem')

EPS = 0.00001

def mandelbrot_plot(f=None, **kwds):
    r"""
    Plot of the Mandelbrot set for a one parameter family of polynomial maps.

    The family `f_c(z)` must have parent ``R`` of the
    form ``R.<z,c> = CC[]``.

    REFERENCE:

    [Dev2005]_

    INPUT:

    - ``f`` -- map (default: ``z^2 + c``); polynomial family used to
      plot the Mandelbrot set

    - ``parameter`` -- variable (default: ``c``); parameter variable
      used to plot the Mandelbrot set

    - ``x_center`` -- double (default: ``-1.0``); Real part of center
      point

    - ``y_center`` -- double (default: ``0.0``); Imaginary part of
      center point

    - ``image_width`` -- double (default: ``4.0``); width of image
      in the complex plane

    - ``max_iteration`` -- long (default: ``500``); maximum number of
      iterations the map ``f_c(z)``

    - ``pixel_count`` -- long (default: ``500``); side length of
      image in number of pixels

    - ``base_color`` -- RGB color (default: ``[40, 40, 40]``); color
      used to determine the coloring of set

    - ``level_sep`` -- long (default: 1); number of iterations
      between each color level

    - ``number_of_colors`` -- long (default: 30); number of colors
      used to plot image

    - ``interact`` -- boolean (default: ``False``); controls whether
      plot will have interactive functionality

    OUTPUT: 24-bit RGB image of the Mandelbrot set in the complex plane

    EXAMPLES:

    ::

        sage: mandelbrot_plot()
        500x500px 24-bit RGB image

    ::

        sage: mandelbrot_plot(pixel_count=1000)
        1000x1000px 24-bit RGB image

    ::

        sage: mandelbrot_plot(x_center=-1.11, y_center=0.2283, image_width=1/128, # long time
        ....: max_iteration=2000, number_of_colors=500, base_color=[40, 100, 100])
        500x500px 24-bit RGB image

    To display an interactive plot of the Mandelbrot in the Notebook, set
    ``interact`` to ``True``. (This is only implemented for ``z^2 + c``)::

        sage: mandelbrot_plot(interact=True)
        ...interactive(children=(FloatSlider(value=0.0, description='Real center', max=1.0, min=-1.0, step=1e-05),
        FloatSlider(value=0.0, description='Imag center', max=1.0, min=-1.0, step=1e-05),
        FloatSlider(value=4.0, description='Width', max=4.0, min=1e-05, step=1e-05),
        IntSlider(value=500, description='Iterations', max=1000),
        IntSlider(value=500, description='Pixels', max=1000, min=10),
        IntSlider(value=1, description='Color sep', max=20, min=1),
        IntSlider(value=30, description='# Colors', min=1),
        ColorPicker(value='#ff6347', description='Base color'), Output()),
        _dom_classes=('widget-interact',))

    ::

        sage: mandelbrot_plot(interact=True, x_center=-0.75, y_center=0.25,
        ....: image_width=1/2, number_of_colors=75)
        ...interactive(children=(FloatSlider(value=-0.75, description='Real center', max=1.0, min=-1.0, step=1e-05),
        FloatSlider(value=0.25, description='Imag center', max=1.0, min=-1.0, step=1e-05),
        FloatSlider(value=0.5, description='Width', max=4.0, min=1e-05, step=1e-05),
        IntSlider(value=500, description='Iterations', max=1000),
        IntSlider(value=500, description='Pixels', max=1000, min=10),
        IntSlider(value=1, description='Color sep', max=20, min=1),
        IntSlider(value=75, description='# Colors', min=1),
        ColorPicker(value='#ff6347', description='Base color'), Output()),
        _dom_classes=('widget-interact',))

    Polynomial maps can be defined over a multivariate polynomial ring or a
    univariate polynomial ring tower::

        sage: R.<z,c> = CC[]
        sage: f = z^2 + c
        sage: mandelbrot_plot(f)
        500x500px 24-bit RGB image

    ::

        sage: B.<c> = CC[]
        sage: R.<z> = B[]
        sage: f = z^5 + c
        sage: mandelbrot_plot(f)  # long time
        500x500px 24-bit RGB image

    When the polynomial is defined over a multivariate polynomial ring it is
    necessary to specify the parameter variable (default parameter is ``c``)::

        sage: R.<a,b> = CC[]
        sage: f = a^2 + b^3
        sage: mandelbrot_plot(f, parameter=b)  # long time
        500x500px 24-bit RGB image

    Interact functionality is not implemented for general polynomial maps::

        sage: R.<z,c> = CC[]
        sage: f = z^3 + c
        sage: mandelbrot_plot(f, interact=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: interact only implemented for z^2 + c
    """
    parameter = kwds.pop("parameter", None)
    x_center = kwds.pop("x_center", 0.0)
    y_center = kwds.pop("y_center", 0.0)
    image_width = kwds.pop("image_width", 4.0)
    max_iteration = kwds.pop("max_iteration", None)
    pixel_count = kwds.pop("pixel_count", 500)
    level_sep = kwds.pop("level_sep", 1)
    number_of_colors = kwds.pop("number_of_colors", 30)
    interacts = kwds.pop("interact", False)
    base_color = kwds.pop("base_color", Color('tomato'))
    # Check if user specified maximum number of iterations
    given_iterations = True
    if max_iteration is None:
        # Set default to 500 for z^2 + c map
        max_iteration = 500
        given_iterations = False

    from ipywidgets.widgets import FloatSlider, IntSlider, ColorPicker, interact
    widgets = dict(
                   x_center=FloatSlider(min=-1.0, max=1.0, step=EPS,
                                          value=x_center, description="Real center"),
                   y_center=FloatSlider(min=-1.0, max=1.0, step=EPS,
                                          value=y_center, description="Imag center"),
                   image_width=FloatSlider(min=EPS, max=4.0, step=EPS,
                                             value=image_width, description='Width'),
                   max_iteration=IntSlider(min=0, max=1000,
                                             value=max_iteration, description='Iterations'),
                   pixel_count=IntSlider(min=10, max=1000,
                                           value=pixel_count, description='Pixels'),
                   level_sep=IntSlider(min=1, max=20,
                                         value=level_sep, description="Color sep"),
                   color_num=IntSlider(min=1, max=100,
                                         value=number_of_colors, description="# Colors"),
                   base_color=ColorPicker(value=Color(base_color).html_color(),
                                            description="Base color"),
                   )

    if f is None:
        # Quadratic map f = z^2 + c

        if interacts:
            return interact(**widgets).widget(fast_mandelbrot_plot)

        else:
            return fast_mandelbrot_plot(x_center, y_center, image_width,
             max_iteration, pixel_count, level_sep, number_of_colors,
             base_color)

    else:
        if parameter is None:
            c = var('c')
            parameter = c

        P = f.parent()

        if P.base_ring() is CC or P.base_ring() is CDF:
            if isinstance(P, FractionField_generic):
                raise NotImplementedError("coefficients must be polynomials in the parameter")
            gen_list = list(P.gens())
            parameter = gen_list.pop(gen_list.index(parameter))
            variable = gen_list.pop()

        elif P.base_ring().base_ring() is CC or P.base_ring().base_ring() is CDF:
            if isinstance(P.base_ring(), FractionField_generic):
                raise NotImplementedError("coefficients must be polynomials in the parameter")
            phi = P.flattening_morphism()
            f = phi(f)
            gen_list = list(f.parent().gens())
            parameter = gen_list.pop(gen_list.index(parameter))
            variable = gen_list.pop()

        elif P.base_ring() in FunctionFields():
            raise NotImplementedError("coefficients must be polynomials in the parameter")

        else:
            raise ValueError("base ring must be a complex field")

        if f == variable**2 + parameter:
            # Quadratic map f = z^2 + c
            if interacts:
                return interact(**widgets).widget(fast_mandelbrot_plot)

            else:
                return fast_mandelbrot_plot(x_center, y_center, image_width,
                 max_iteration, pixel_count, level_sep, number_of_colors,
                 base_color)
        else:
            if interacts:
                raise NotImplementedError("interact only implemented for z^2 + c")
            else:
                # Set default of max_iteration to 50 for general polynomial maps
                # This prevents the function from being very slow by default
                if not given_iterations:
                    max_iteration = 50

                # Mandelbrot of General Polynomial Map
                return polynomial_mandelbrot(f, parameter, x_center, y_center,
                 image_width, max_iteration, pixel_count, level_sep,
                 number_of_colors, base_color)


def external_ray(theta, **kwds):
    r"""
    Draws the external ray(s) of a given angle (or list of angles)
    by connecting a finite number of points that were approximated using
    Newton's method. The algorithm used is described in a paper by
    Tomoki Kawahira.

    REFERENCE:

    [Kaw2009]_

    INPUT:

    - ``theta`` -- double or list of doubles, angles between 0 and 1 inclusive

    kwds:

    - ``image`` -- 24-bit RGB image (default: ``None``); user specified
      image of Mandelbrot set

    - ``D`` -- long (default: ``25``); depth of the approximation.
      As ``D`` increases, the external ray gets closer to the boundary of the
      Mandelbrot set. If the ray doesn't reach the boundary of the Mandelbrot
      set, increase ``D``.

    - ``S`` -- long (default: ``10``); sharpness of the approximation.
      Adjusts the number of points used to approximate the external ray (number
      of points is equal to ``S*D``). If ray looks jagged, increase ``S``.

    - ``R`` -- long (default: ``100``); radial parameter. If ``R`` is
      large, the external ray reaches sufficiently close to infinity. If ``R`` is
      too small, Newton's method may not converge to the correct ray.

    - ``prec`` -- long (default: ``300``); specifies the bits of
      precision used by the Complex Field when using Newton's method to compute
      points on the external ray

    - ``ray_color`` -- RGB color (default: ``[255, 255, 255]``) color
      of the external ray(s)

    OUTPUT: 24-bit RGB image of external ray(s) on the Mandelbrot set

    EXAMPLES::

        sage: external_ray(1/3)
        500x500px 24-bit RGB image

    ::

        sage: external_ray(0.6, ray_color=[255, 0, 0])
        500x500px 24-bit RGB image

    ::

        sage: external_ray([0, 0.2, 0.4, 0.7])
        500x500px 24-bit RGB image

    ::

        sage: external_ray([i/5 for i in range(1,5)])
        500x500px 24-bit RGB image

    WARNING:

    If you are passing in an image, make sure you specify
    which parameters to use when drawing the external ray.
    For example, the following is incorrect::

        sage: M = mandelbrot_plot(x_center=0)  # not tested
        sage: external_ray(5/7, image=M)       # not tested
        500x500px 24-bit RGB image

    To get the correct external ray, we adjust our parameters::

        sage: M = mandelbrot_plot(x_center=0)
        sage: external_ray(5/7, x_center=0, image=M)
        500x500px 24-bit RGB image

    .. TODO::

        The ``copy()`` function for bitmap images needs to be implemented
        in Sage.
    """
    x_0 = kwds.get("x_center", -1)
    y_0 = kwds.get("y_center", 0)
    plot_width = kwds.get("image_width", 4)
    pixel_width = kwds.get("pixel_count", 500)
    depth = kwds.get("D", 25)
    sharpness = kwds.get("S", 10)
    radial_parameter = kwds.get("R", 100)
    precision = kwds.get("prec", 300)
    precision = max(precision, -logb(pixel_width * 0.001, 2).round() + 10)
    ray_color = kwds.get("ray_color", [255] * 3)
    image = kwds.get("image", None)
    if image is None:
        image = mandelbrot_plot(x_center=x_0, **kwds)

    # Make a copy of the bitmap image.
    old_pixel = image.pixels()
    M = Image('RGB', (pixel_width, pixel_width))
    pixel = M.pixels()
    for i in range(pixel_width):
        for j in range(pixel_width):
            pixel[i, j] = old_pixel[i, j]

    # Make sure that theta is a list so loop below works
    if not isinstance(theta, list):
        theta = [theta]

    # Check if theta is in the interval [0,1]
    for angle in theta:
        if angle < 0 or angle > 1:
            raise ValueError("values for theta must be in "
                             "the closed interval [0,1].")

    # Loop through each value for theta in list and plot the external ray.
    for angle in theta:
        E = fast_external_ray(angle, D=depth, S=sharpness, R=radial_parameter,
                              prec=precision, image_width=plot_width,
                              pixel_count=pixel_width)

        # Convert points to pixel coordinates.
        pixel_list = convert_to_pixels(E, x_0, y_0, plot_width, pixel_width)

        # Find the pixels between points in pixel_list.
        extra_points = []
        for i in range(len(pixel_list) - 1):
            if min(pixel_list[i + 1]) >= 0 and max(pixel_list[i + 1]) < pixel_width:
                for j in get_line(pixel_list[i], pixel_list[i + 1]):
                    extra_points.append(j)

        # Add these points to pixel_list to fill in gaps in the ray.
        pixel_list += extra_points

        # Remove duplicates from list.
        pixel_list = list(set(pixel_list))

        # Check if point is in window and if it is, plot it on the image to
        # create an external ray.
        for k in pixel_list:
            if max(k) < pixel_width and min(k) >= 0:
                pixel[int(k[0]), int(k[1])] = tuple(ray_color)
    return M


def kneading_sequence(theta):
    r"""
    Determines the kneading sequence for an angle theta in RR/ZZ which
    is periodic under doubling. We use the definition for the kneading
    sequence given in Definition 3.2 of [LS1994]_.

    INPUT:

    - ``theta`` -- a rational number with odd denominator

    OUTPUT: string representing the kneading sequence of theta in RR/ZZ

    REFERENCES:

    [LS1994]_

    EXAMPLES::

        sage: kneading_sequence(0)
        '*'

    ::

        sage: kneading_sequence(1/3)
        '1*'

    Since 1/3 and 7/3 are the same in RR/ZZ, they have the same kneading sequence::

        sage: kneading_sequence(7/3)
        '1*'

    We can also use (finite) decimal inputs, as long as the denominator in reduced form is odd::

        sage: kneading_sequence(1.2)
        '110*'

    Since rationals with even denominator are not periodic under doubling, we have not implemented kneading sequences for such rationals::

        sage: kneading_sequence(1/4)
        Traceback (most recent call last):
        ...
        ValueError: input must be a rational number with odd denominator
    """

    if theta not in QQ:
        raise TypeError('input must be a rational number with odd denominator')
    elif QQ(theta).valuation(2) < 0:
        raise ValueError('input must be a rational number with odd denominator')
    else:
        theta = QQ(theta)
        theta = theta - floor(theta)
        KS = []
        not_done = True
        left = theta/2
        right = (theta + 1)/2
        y = theta
        while not_done:
            if ((y < left) or (y > right)):
                KS.append('0')
            elif ((y > left) and (y < right)):
                KS.append('1')
            else:
                not_done = False
            y = 2*y - floor(2*y)
        KS_str = ''.join(KS) + '*'
    return KS_str

def julia_plot(f=None, **kwds):
    r"""
    Plots the Julia set of a given polynomial ``f``. Users can specify whether
    they would like to display the Mandelbrot side by side with the Julia set
    with the ``mandelbrot`` argument. If ``f`` is not specified, this method
    defaults to `f(z) = z^2-1`.

    The Julia set of a polynomial ``f`` is the set of complex numbers `z` for
    which the function `f(z)` is bounded under iteration. The Julia set can
    be visualized by plotting each point in the set in the complex plane.
    Julia sets are examples of fractals when plotted in the complex plane.

    ALGORITHM:

    Let `R_c = \bigl(1 + \sqrt{1 + 4|c|}\bigr)/2` if the polynomial is of the
    form `f(z) = z^2 + c`; otherwise, let `R_c = 2`.
    For every `p \in \mathbb{C}`, if `|f^{k}(p)| > R_c` for some `k \geq 0`,
    then `f^{n}(p) \to \infty`.  Let `N` be the maximum number of iterations.
    Compute the first `N` points on the orbit of `p` under `f`. If for
    any `k < N`, `|f^{k}(p)| > R_c`, we stop the iteration and assign a color
    to the point `p` based on how quickly `p` escaped to infinity under
    iteration of `f`. If `|f^{i}(p)| \leq R_c` for all `i \leq N`, we assume
    `p` is in the Julia set and assign the point `p` the color black.

    INPUT:

    - ``f`` -- input polynomial (default: ``z^2 - 1``)

    - ``period`` -- list (default: ``None``); returns the Julia set
      for a random `c` value with the given (formal) cycle structure

    - ``mandelbrot`` -- boolean (default: ``True``); when set to
      ``True``, an image of the Mandelbrot set is appended to the right of the
      Julia set

    - ``point_color`` -- RGB color (default: ``'tomato'``);
      color of the point `c` in the Mandelbrot set (any valid input for Color)

    - ``x_center`` -- double (default: ``-1.0``); real part
      of center point

    - ``y_center`` -- double (default: ``0.0``); imaginary part
      of center point

    - ``image_width`` -- double (default: ``4.0``); width of image
      in the complex plane

    - ``max_iteration`` -- long (default: ``500``); maximum number
      of iterations the map `f(z)`

    - ``pixel_count`` -- long (default: ``500``); side length of
      image in number of pixels

    - ``base_color`` -- hex color (default: ``'steelblue'``); color
      used to determine the coloring of set (any valid input for Color)

    - ``level_sep`` -- long (default: 1); number of iterations
      between each color level

    - ``number_of_colors`` -- long (default: 30); number of colors
      used to plot image

    - ``interact`` -- boolean (default: ``False``); controls whether
      plot will have interactive functionality

    OUTPUT: 24-bit RGB image of the Julia set in the complex plane

    .. TODO::

        Implement the side-by-side Mandelbrot-Julia plots for general one-parameter families
        of polynomials.

    EXAMPLES:

    The default ``f`` is `z^2 - 1`::

        sage: julia_plot()
        1001x500px 24-bit RGB image

    To display only the Julia set, set ``mandelbrot`` to ``False``::

        sage: julia_plot(mandelbrot=False)
        500x500px 24-bit RGB image

    ::

        sage: R.<z> = CC[]
        sage: f = z^3 - z + 1
        sage: julia_plot(f)  # long time
        500x500px 24-bit RGB image

    To display an interactive plot of the Julia set in the Notebook,
    set ``interact`` to ``True``. (This is only implemented for polynomials of
    the form ``f = z^2 + c``)::

        sage: julia_plot(interact=True)
        ...interactive(children=(FloatSlider(value=-1.0, description='Real c'...

        ::

        sage: R.<z> = CC[]
        sage: f = z^2 + 1/2
        sage: julia_plot(f,interact=True)
        ...interactive(children=(FloatSlider(value=0.5, description='Real c'...

    To return the Julia set of a random `c` value with (formal) cycle structure
    `(2,3)`, set ``period = [2,3]``::

        sage: julia_plot(period=[2,3])
        1001x500px 24-bit RGB image

    To return all of the Julia sets of `c` values with (formal) cycle structure
    `(2,3)`::

        sage: period = [2,3] # not tested
        ....: R.<c> = QQ[]
        ....: P.<x,y> = ProjectiveSpace(R,1)
        ....: f = DynamicalSystem([x^2+c*y^2, y^2])
        ....: L = f.dynatomic_polynomial(period).subs({x:0,y:1}).roots(ring=CC)
        ....: c_values = [k[0] for k in L]
        ....: for c in c_values:
        ....:     julia_plot(c)

    Polynomial maps can be defined over a polynomial ring or a fraction field,
    so long as ``f`` is polynomial::

        sage: R.<z> = CC[]
        sage: f = z^2 - 1
        sage: julia_plot(f)
        1001x500px 24-bit RGB image

    ::

        sage: R.<z> = CC[]
        sage: K = R.fraction_field(); z = K.gen()
        sage: f = z^2-1
        sage: julia_plot(f)
        1001x500px 24-bit RGB image

    Interact functionality is not implemented if the polynomial is not of the
    form `f = z^2 + c`::

        sage: R.<z> = CC[]
        sage: f = z^3 + 1
        sage: julia_plot(f, interact=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: The interactive plot is only implemented for ...
    """

    # extract keyword arguments
    period = kwds.pop("period", None)
    mandelbrot = kwds.pop("mandelbrot", True)
    point_color = kwds.pop("point_color", 'tomato')
    x_center = kwds.pop("x_center", 0.0)
    y_center = kwds.pop("y_center", 0.0)
    image_width = kwds.pop("image_width", 4.0)
    max_iteration = kwds.pop("max_iteration", 500)
    pixel_count = kwds.pop("pixel_count", 500)
    base_color = kwds.pop("base_color", 'steelblue')
    level_sep = kwds.pop("level_sep", 1)
    number_of_colors = kwds.pop("number_of_colors", 30)
    interacts = kwds.pop("interact", False)

    f_is_default_after_all = None

    if period: # pick a random c with the specified period
        R = PolynomialRing(CC, 'c')
        c = R.gen()
        x, y = ProjectiveSpace(R, 1, 'x,y').gens()
        F = DynamicalSystem([x**2 + c * y**2, y**2])
        L = F.dynatomic_polynomial(period).subs({x: 0, y: 1}).roots(ring=CC)
        c = L[randint(0, len(L) - 1)][0]

    base_color = Color(base_color)
    point_color = Color(point_color)

    EPS = 0.00001

    if f is not None and period is None: # f user-specified and no period given

        # try to coerce f to live in a polynomial ring
        S = PolynomialRing(CC, names='z')
        z = S.gen()
        try:
            f_poly = S(f)
        except TypeError:
            R = f.parent()
            if not (R.is_integral_domain() and
                   (CC.is_subring(R) or CDF.is_subring(R))):
                raise ValueError('given `f` must be a complex polynomial')
            raise NotImplementedError(
                'Julia sets not implemented for rational functions')

        if (f_poly - z*z) in CC: # f is specified and of the form z^2 + c.
            f_is_default_after_all = True
            c = f_poly - z*z
        else: # f is specified and not of the form z^2 + c
            if interacts:
                raise NotImplementedError(
                        "The interactive plot is only implemented for "
                        "polynomials of the form f = z^2 + c."
                            )
            else:
                return general_julia(f_poly, x_center, y_center,
                                    image_width, max_iteration,
                                    pixel_count, level_sep,
                                    number_of_colors, base_color)

    # otherwise we can use fast_julia_plot for z^2 + c
    if f_is_default_after_all or f is None or period is not None:

        # specify default c = -1 value if f and period were not specified
        if not f_is_default_after_all and period is None:
            c = -1

        c = CC(c)
        c_real = c.real()
        c_imag = c.imag()

        if interacts: # set widgets
            from ipywidgets.widgets import FloatSlider, IntSlider, \
                                           ColorPicker, interact
            widgets = dict(
                c_real=FloatSlider(min=-2.0, max=2.0, step=EPS,
                                  value=c_real, description="Real c"),
                c_imag=FloatSlider(min=-2.0, max=2.0, step=EPS,
                                  value=c_imag, description="Imag c"),
                x_center=FloatSlider(min=-1.0, max=1.0, step=EPS,
                                  value=x_center, description="Real center"),
                y_center=FloatSlider(min=-1.0, max=1.0, step=EPS,
                                  value=y_center, description="Imag center"),
                image_width=FloatSlider(min=EPS, max=4.0, step=EPS,
                                  value=image_width, description='Width'),
                max_iteration=IntSlider(min=0, max=1000,
                                  value=max_iteration, description='Iterations'),
                pixel_count=IntSlider(min=10, max=1000,
                                  value=pixel_count, description='Pixels'),
                level_sep=IntSlider(min=1, max=20,
                                  value=level_sep, description="Color sep"),
                color_num=IntSlider(min=1, max=100,
                                  value=number_of_colors, description="# Colors"),
                base_color=ColorPicker(value=base_color.html_color(),
                                         description="Base color"),
            )
            if mandelbrot:
                widgets["point_color"] = ColorPicker(value=point_color.html_color(),
                                                     description="Point color")
                return interact(**widgets).widget(julia_helper)
            else:
                return interact(**widgets).widget(fast_julia_plot)
        elif mandelbrot:  # non-interactive with mandelbrot
            return julia_helper(c_real, c_imag, x_center, y_center,
                                image_width, max_iteration, pixel_count,
                                level_sep, number_of_colors, base_color,
                                point_color)
        else:  # non-interactive without mandelbrot
            return fast_julia_plot(c_real, c_imag, x_center, y_center,
                                   image_width, max_iteration, pixel_count,
                                   level_sep, number_of_colors, base_color)
