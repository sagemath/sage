# cython: binding=True
# sage.doctest: needs sage.plot
r"""
Mandelbrot and Julia sets (Cython helper)

This is the helper file providing functionality for mandel_julia.py.

AUTHORS:

- Ben Barros
"""
# ****************************************************************************
#       Copyright (C) 2017 BEN BARROS <bbarros@slu.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.plot.colors import Color
from sage.repl.image import Image
from copy import copy
from cysignals.signals cimport sig_check
from sage.rings.complex_mpfr import ComplexField
from sage.functions.log import exp
from sage.symbolic.relation import solve
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.cc import CC
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.ext.fast_callable import fast_callable
from sage.calculus.all import symbolic_expression
from sage.symbolic.ring import SR
from sage.calculus.var import var
from sage.rings.fraction_field import FractionField_generic
from sage.categories.function_fields import FunctionFields
from cypari2.handle_error import PariError
from math import sqrt


def _color_to_RGB(color):
    """
    Convert a color to an RGB triple with values in the interval [0,255].

    EXAMPLES::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import _color_to_RGB
        sage: _color_to_RGB("aquamarine")
        (127, 255, 212)
        sage: _color_to_RGB(Color([0, 1/2, 1]))
        (0, 127, 255)
        sage: _color_to_RGB([0, 100, 200])
        (0, 100, 200)
    """
    if not isinstance(color, (list, tuple)):
        color = [int(255.0 * k) for k in Color(color)]
    return tuple(color)


cpdef fast_mandelbrot_plot(double x_center, double y_center,
                           double image_width, long max_iteration,
                           long pixel_count,
                           long level_sep, long color_num, base_color):
    r"""
    Plot the Mandelbrot set in the complex plane for the map `Q_c(z) = z^2 + c`.

    ALGORITHM:

    Let each pixel in the image be a point `c \in \mathbb{C}` and define the
    map `Q_c(z) = z^2 + c`. If `|Q_{c}^{k}(c)| > 2` for some `k \geq 0`, it
    follows that `Q_{c}^{n}(c) \to \infty`. Let `N` be the maximum number of
    iterations. Compute the first `N` points on the orbit of `0` under `Q_c`.
    If for any `k < N`, `|Q_{c}^{k}(0)| > 2`, we stop the iteration and assign
    a color to the point `c` based on how quickly `0` escaped to infinity under
    iteration of `Q_c`. If `|Q_{c}^{i}(0)| \leq 2` for all `i \leq N`, we assume
    `c` is in the Mandelbrot set and assign the point `c` the color black.

    INPUT:

    - ``x_center`` -- double; real part of the center point in the complex plane

    - ``y_center`` -- double; imaginary part of the center point in the complex
      plane

    - ``image_width`` -- double; width of the image in the complex plane

    - ``max_iteration`` -- long; maximum number of iterations the map `Q_c(z)`
      considered

    - ``pixel_count`` -- long; side length of image in number of pixels

    - ``level_sep`` -- long; number of iterations between each color level

    - ``color_num`` -- long; number of colors used to plot image

    - ``base_color`` -- list; RGB color used to determine the coloring of set

    OUTPUT: 24-bit RGB image of the Mandelbrot set in the complex plane

    EXAMPLES:

    Plot the Mandelbrot set with the center point `-1 + 0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandelbrot_plot
        sage: fast_mandelbrot_plot(-1, 0, 4, 500, 600, 1, 20, [40, 40, 40])
        600x600px 24-bit RGB image

    We can focus on smaller parts of the set by adjusting image_width::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandelbrot_plot
        sage: fast_mandelbrot_plot(-1.11, 0.2283, 1/128, 2000, 500, 1, 500, [40, 100, 100])
        500x500px 24-bit RGB image
    """

    cdef:
        M, pixel, color_list
        long i, j, col, row, level, iteration
        double x_corner, y_corner, step_size, x_coor, y_coor, new_x, new_y

    # Make sure image_width is positive
    image_width = abs(image_width)

    # Initialize an image to the color black and access the pixels
    M = Image("RGB", (pixel_count, pixel_count), 'black')
    pixel = M.pixels()

    # Take the given base color and create a list of evenly spaced
    # colors between the given base color and white. The number of
    # colors in the list depends on the variable color_num.
    base_color = _color_to_RGB(base_color)
    color_list = []
    for i in range(color_num):
        sig_check()
        color = [base_color[j] + i * (255 - base_color[j]) // color_num
                 for j in range(3)]
        color_list.append(tuple(color))

    # First, we determine the complex coordinates of the point in the top left
    # corner of the image. Then, we loop through each pixel in the image and
    # assign it complex coordinates relative to the image's top left corner.
    x_corner = x_center - image_width/2
    y_corner = y_center + image_width/2
    step_size = image_width / pixel_count
    for col in range(pixel_count):
        x_coor = x_corner + col*step_size
        for row in range(pixel_count):
            sig_check()
            y_coor = y_corner - row*step_size

            # We compute the orbit of 0 under the map Q(z) = z^2 + c
            # until we either reach the maximum number of iterations
            # or find a point in the orbit with modulus greater than 2
            new_x, new_y = 0.0, 0.0
            iteration = 0
            while (new_x**2 + new_y**2 <= 4.0 and iteration < max_iteration):
                sig_check()
                new_x, new_y = new_x**2 - new_y**2 + x_coor, \
                    2*new_x*new_y + y_coor
                iteration += 1

            # If the point escapes to infinity, assign the point a color
            # based on how fast it escapes. The more iterations it takes for
            # a point to escape to infinity, the lighter its color will be.
            # Otherwise, assume the point is in the Mandelbrot set and leave
            # it black.
            if iteration != max_iteration:
                # Assign each point a level based on its number of iterations.
                level = iteration // level_sep
                # Assign the pixel a color based on it's level. If we run out
                # of colors, assign it the last color in the list.
                if level < color_num:
                    pixel[col, row] = color_list[level]
                else:
                    pixel[col, row] = color_list[-1]
    return M


cpdef fast_external_ray(double theta, long D=30, long S=10, long R=100,
                        long pixel_count=500, double image_width=4,
                        long prec=300):
    r"""
    Return a list of points that approximate the external ray for a given angle.

    INPUT:

    - ``theta`` -- double; angle between 0 and 1 inclusive

    - ``D`` -- long (default: ``25``); depth of the approximation.
     As ``D`` increases, the external ray gets closer to the boundary of the
     Mandelbrot set.

    - ``S`` -- long (default: ``10``); sharpness of the approximation.
     Adjusts the number of points used to approximate the external ray (number
     of points is equal to ``S*D``).

    - ``R`` -- long (default: ``100``); radial parameter. If ``R`` is
     sufficiently large, the external ray reaches enough close to infinity.

    - ``pixel_count`` -- long (default: ``500``); side length of image
     in number of pixels

    - ``image_width`` -- double (default: ``4``); width of the image
     in the complex plane

    - ``prec`` -- long (default: ``300``); specifies the bits of
     precision used by the Complex Field when using Newton's method to compute
     points on the external ray

    OUTPUT: list of tuples of Real Interval Field Elements

    EXAMPLES::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_external_ray
        sage: fast_external_ray(0,S=1,D=1)
        [(100.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
          0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000),
         (9.51254777713729174697578576623132297117784691109499464854806785133621315075854778426714908,
          0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000)]


    ::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_external_ray
        sage: fast_external_ray(1/3,S=1,D=1)
        [(-49.9999999999999786837179271969944238662719726562500000000000000000000000000000000000000000,
          86.6025403784438765342201804742217063903808593750000000000000000000000000000000000000000000),
         (-5.50628047023173006234970878097113901879832542655926629309001652388544528575532346900138516,
          8.64947510053972513843999918917106032664030380426885745306040284140385975750462108180377187)]

    ::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_external_ray
        sage: fast_external_ray(0.75234,S=1,D=1)
        [(1.47021239172637052661229972727596759796142578125000000000000000000000000000000000000000000,
          -99.9891917935294287644865107722580432891845703125000000000000000000000000000000000000000000),
         (-0.352790406744857508500937144524776555433184352559852962308757189778284058275081335121601384,
          -9.98646630765023514178761177926164047797465369576787921409326037870837930920646860774032363)]
    """

    cdef:
        CF = ComplexField(prec)
        PI = CF.pi()
        I = CF.gen()
        r_m, t_m, temp_c, C_k, D_k, old_c, x, y, dist
        int k, j, t
        double difference, m
        double error = pixel_count * 0.0001

        double pixel_width = image_width / pixel_count

        # initialize list with c_0
        c_list = [CF(R*exp(2*PI*I*theta))]

    # Loop through each subinterval and approximate point on external ray.
    for k in range(1, D+1):
        for j in range(1, S+1):
            m = (k-1)*S + j
            r_m = CF(R**(2**(-m/S)))
            t_m = CF(r_m**(2**k) * exp(2*PI*I*theta * 2**k))
            temp_c = c_list[-1]
            difference = error

            # Repeat Newton's method until points are close together.
            while error <= difference:
                sig_check()
                old_c = temp_c
                # Recursive formula for iterates of q(z) = z^2 + c
                C_k, D_k = CF(old_c), CF(1)
                for t in range(k):
                    C_k, D_k = C_k**2 + old_c, CF(2)*D_k*C_k + CF(1)
                temp_c = old_c - (C_k - t_m) / D_k   # Newton map
                difference = abs(old_c) - abs(temp_c)

            dist = (2*C_k.abs()*(C_k.abs()).log()) / D_k.abs()
            if dist < pixel_width:
                break
            c_list.append(CF(temp_c))
        if dist < pixel_width:
            break

    # Convert Complex Field elements into tuples.
    for k in range(len(c_list)):
        x, y = c_list[k].real(), c_list[k].imag()
        c_list[k] = (x, y)

    return c_list

cpdef convert_to_pixels(point_list, double x_0, double y_0, double width,
                        long number_of_pixels):
    r"""
    Convert cartesian coordinates to pixels within a specified window.

    INPUT:

    - ``point_list`` -- list of tuples; points in cartesian coordinates

    - ``x_0`` -- double; x-coordinate of the center of the image

    - ``y_0`` -- double; y-coordinate of the center of the image

    - ``width`` -- double; width of visible window in cartesian coordinates

    - ``number_of_pixels`` -- long; width of image in pixels

    OUTPUT: list of tuples of integers representing pixels

    EXAMPLES::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import convert_to_pixels
        sage: convert_to_pixels([(-1,3),(0,-4),(5,0)], 0, 0, 12, 100)
        [(42, 25), (50, 83), (92, 50)]
    """
    cdef:
        k, pixel_list, x_corner, y_corner, step_size
        long x_pixel, y_pixel
    pixel_list = []

    # Compute top left corner of window and step size
    x_corner = x_0 - width/2
    y_corner = y_0 + width/2
    step_size = number_of_pixels / width

    # Convert each point in list to pixel coordinates
    for k in point_list:
        sig_check()
        x_pixel = round((k[0] - x_corner) * step_size)
        y_pixel = round((y_corner - k[1]) * step_size)
        pixel_list.append((x_pixel, y_pixel))
    return pixel_list

cpdef get_line(start, end):
    r"""
    Produce a list of pixel coordinates approximating a line from a starting
    point to an ending point using the Bresenham's Line Algorithm.

    REFERENCE:

    [Br2016]_

    INPUT:

    - ``start`` -- tuple; starting point of line

    - ``end`` -- tuple; ending point of line

    OUTPUT: list of tuples of integers approximating the line between two pixels

    EXAMPLES::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import get_line
        sage: get_line((0, 0), (3, 4))
        [(0, 0), (1, 1), (1, 2), (2, 3), (3, 4)]

    ::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import get_line
        sage: get_line((3, 4), (0, 0))
        [(3, 4), (2, 3), (1, 2), (1, 1), (0, 0)]
    """
    # Setup initial conditions
    cdef:
        long x1, x2, y1, y2, dx, dy, error, ystep, y
        is_steep, swapped, points
    x1, y1 = start
    x2, y2 = end
    dx, dy = x2 - x1, y2 - y1

    # Determine how steep the line is
    is_steep = abs(dy) > abs(dx)

    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2

    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True

    # Recalculate differentials
    dx, dy = x2 - x1, y2 - y1

    # Calculate error
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1

    # Iterate over bounding box generating points between start and end
    y = y1
    points = []
    for x in range(x1, x2 + 1):
        sig_check()
        coord = (y, x) if is_steep else (x, y)
        points.append(coord)
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx

    # Reverse the list if the coordinates were swapped
    if swapped:
        points.reverse()
    return points

# Commented out temporarily for safekeeping, but probably should be deleted
# def fast_julia_plot(double c_real, double c_imag,
#                     double x_center, double y_center, double image_width,
#                     long max_iteration, long pixel_count, long level_sep,
#                     long color_num, base_color):

cpdef fast_julia_plot(double c_real, double c_imag,
                      double x_center=0, double y_center=0,
                      double image_width=4,
                      long max_iteration=500, long pixel_count=500,
                      long level_sep=2,
                      long color_num=40, base_color=[50, 50, 50]):
    r"""
    Plot the Julia set for a given `c` value in the complex plane for the map `Q_c(z) = z^2 + c`.

    INPUT:

    - ``c_real`` -- double; Real part of `c` value that determines Julia set

    - ``c_imag`` -- double; Imaginary part of `c` value that determines Julia
      set

    - ``x_center`` -- double (default: ``0.0``); real part of center
      point

    - ``y_center`` -- double (default: ``0.0``); imaginary part of
      center point

    - ``image_width`` -- double (default: ``4.0``); width of image
      in the complex plane

    - ``max_iteration`` -- long (default: ``500``); maximum number of
      iterations the map ``Q_c(z)``

    - ``pixel_count`` -- long (default: ``500``); side length of
      image in number of pixels

    - ``level_sep`` -- long (default: ``2``); number of iterations
      between each color level

    - ``color_num`` -- long (default: ``40``); number of colors used
      to plot image

    - ``base_color`` -- RGB color (default: ``[50, 50, 50]``); color
      used to determine the coloring of set

    OUTPUT: 24-bit RGB image of the Julia set in the complex plane

    EXAMPLES:

    Plot the Julia set for `c=-1+0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_julia_plot
        sage: fast_julia_plot(-1, 0, 0, 0, 4, 500, 200, 1, 20, [40, 40, 40])
        200x200px 24-bit RGB image
    """

    cdef:
        color_list
        long i, j, col, row, level, iteration
        double x_corner, y_corner, step_size, x_coor, y_coor, new_x, new_y, escape_radius_squared

    # Make sure image_width is positive
    image_width = abs(image_width)

    # Calculate the escape radius. (Actually, we store the square of this radius.) If any
    # iterate is farther from the origin than this distance, then the orbit goes to infinity.
    escape_radius_squared = ((1.0 + sqrt(1.0 + 4.0*sqrt(c_real**2 + c_imag**2))))**2/4.0

    # Initialize an image to the color black and access the pixels
    J = Image("RGB", (pixel_count, pixel_count), 'black')
    Jp = J.pixels()

    # Take the given base color and create a list of evenly spaced
    # colors between the given base color and white. The number of
    # colors in the list depends on the variable color_num.
    base_color = _color_to_RGB(base_color)
    color_list = []
    for i in range(color_num):
        sig_check()
        color = [base_color[j] + i * (255 - base_color[j]) // color_num
                 for j in range(3)]
        color_list.append(tuple(color))

    # First, we determine the complex coordinates of the point in the top left
    # corner of the image. Then, we loop through each pixel in the image and
    # assign it complex coordinates relative to the image's top left corner.
    x_corner = x_center - image_width/2
    y_corner = y_center + image_width/2
    step_size = image_width / pixel_count
    for col in range(pixel_count):
        x_coor = x_corner + col*step_size
        for row in range(pixel_count):
            sig_check()
            y_coor = y_corner - row*step_size

            # We compute the orbit of each pixel under the map Q(z) = z^2 + c
            # until we either reach the maximum number of iterations
            # or find a point in the orbit with modulus greater than
            # the escape radius.
            new_x, new_y = x_coor, y_coor
            iteration = 0
            while (new_x**2 + new_y**2 <= escape_radius_squared and iteration < max_iteration):
                sig_check()
                new_x, new_y = new_x**2 - new_y**2 + c_real, \
                    2*new_x*new_y + c_imag
                iteration += 1

            # If the point escapes to infinity, assign the point a color
            # based on how fast it escapes. The more iterations it takes for
            # a point to escape to infinity, the lighter its color will be.
            # Otherwise, assume the point is in the Julia set and leave
            # it black.
            if iteration != max_iteration:
                # Assign each point a level based on its number of iterations.
                level = iteration // level_sep
                # Assign the pixel a color based on it's level. If we run out
                # of colors, assign it the last color in the list.
                if level < color_num:
                    Jp[col, row] = color_list[level]
                else:
                    Jp[col, row] = color_list[-1]

    return J

cpdef julia_helper(double c_real, double c_imag, double x_center=0,
                   double y_center=0, double image_width=4,
                   long max_iteration=500,
                   long pixel_count=500, long level_sep=2, long color_num=40,
                   base_color=[50, 50, 50], point_color=[255, 0, 0]):
    r"""
    Helper function that returns the image of a Julia set for a given
    `c` value side by side with the Mandelbrot set with a point denoting
    the `c` value.

    INPUT:

    - ``c_real`` -- double; Real part of `c` value that determines Julia set

    - ``c_imag`` -- double; Imaginary part of `c` value that determines Julia
      set

    - ``x_center`` -- double (default: ``0.0``); Real part of center
      point

    - ``y_center`` -- double (default: ``0.0``); Imaginary part of
      center point

    - ``image_width`` -- double (default: ``4.0``); width of image in
      the complex plane

    - ``max_iteration`` -- long (default: ``500``); maximum number of
      iterations the map ``Q_c(z)``

    - ``pixel_count`` -- long (default: ``500``); side length of
      image in number of pixels

    - ``level_sep`` -- long (default: ``2``); number of iterations
      between each color level

    - ``color_num`` -- long (default: ``40``); number of colors used
      to plot image

    - ``base_color`` -- RGB color (default: ``[50, 50, 50]``); color
      used to determine the coloring of set

    - ``point_color`` -- RGB color (default: ``[255, 0, 0]``); color
      of the point `c` in the Mandelbrot set

    OUTPUT: 24-bit RGB image of the Julia and Mandelbrot sets in the complex plane

    EXAMPLES:

    Plot the Julia set for `c=-1+0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import julia_helper
        sage: julia_helper(-1, 0, 0, 0, 4, 500, 200, 1, 20, [40, 40, 40], [255, 0, 0])
        401x200px 24-bit RGB image
    """
    cdef int i, j

    # Initialize the Julia set
    J = fast_julia_plot(c_real, c_imag, x_center, y_center, image_width,
                        max_iteration, pixel_count, level_sep,
                        color_num, base_color)
    Jp = J.pixels()

    # Initialize the image with Julia set on left side
    # Add white border between images
    G = Image("RGB", (2*pixel_count+1, pixel_count), 'white')
    Gp = G.pixels()
    for i in range(pixel_count):
        for j in range(pixel_count):
            Gp[i, j] = Jp[i, j]

    # Plot the Mandelbrot set on the right side
    M = fast_mandelbrot_plot(-1, 0, 4, 500, pixel_count, 1, 30, base_color)
    Mp = M.pixels()
    for i in range(pixel_count+1, 2*pixel_count):
        for j in range(pixel_count):
            Gp[i, j] = Mp[int(i-pixel_count), j]

    point_color = _color_to_RGB(point_color)

    # Add a cross representing c-value to the Mandelbrot set.
    CP = convert_to_pixels([(c_real, c_imag)], -1, 0, 4, pixel_count)
    for i in range(-3, 4):
        # Loop through x and y coordinates and check if they are in image
        if min(CP[0][0]+i, CP[0][1]) >= 0 and \
           max(CP[0][0]+i, CP[0][1]) < pixel_count:
            Gp[CP[0][0]+i+pixel_count+1, CP[0][1]] = tuple(point_color)
        if min(CP[0][0], CP[0][1]+i) >= 0 and \
           max(CP[0][0], CP[0][1]+i) < pixel_count:
            Gp[CP[0][0]+pixel_count+1, CP[0][1]+i] = tuple(point_color)

    return G

cpdef polynomial_mandelbrot(f, parameter=None, double x_center=0,
                            double y_center=0, image_width=4,
                            int max_iteration=50, int pixel_count=500,
                            int level_sep=1, int color_num=30,
                            base_color=Color('red')):
    r"""
    Plot the Mandelbrot set in the complex plane for a family of polynomial maps.

    INPUT:

    - ``f`` -- a one parameter family of polynomial maps defined over the
      multivariate polynomial ring in z, c over the Complex field

    - ``parameter`` -- designates which variable is used as the parameter
      If no parameter is provided, ``c`` will be used as the parameter

    - ``x_center`` -- double, real part of the center point in the complex plane

    - ``y_center`` -- double, imaginary part of the center point in the complex
      plane

    - ``image_width`` -- double, width of the image in the complex plane

    - ``max_iteration`` -- long, maximum number of iterations the map `f(z)`
      considered

    - ``pixel_count`` -- long, side length of image in number of pixels

    - ``level_sep`` -- long, number of iterations between each color level

    - ``color_num`` -- long, number of colors used to plot image

    - ``base_color`` -- list; RGB color used to determine the coloring of set

    OUTPUT: 24-bit RGB image of a Mandelbrot set in the complex plane

    EXAMPLES::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import polynomial_mandelbrot
        sage: R.<z,c> = CC[]
        sage: f = z^5 + c
        sage: polynomial_mandelbrot(f, pixel_count=100)
        100x100px 24-bit RGB image

    ::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import polynomial_mandelbrot
        sage: B.<c> = CC[]
        sage: R.<z> = B[]
        sage: f = z^4 - z + c
        sage: polynomial_mandelbrot(f, pixel_count=100)
        100x100px 24-bit RGB image

    ::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import polynomial_mandelbrot
        sage: B.<c> = CC[]
        sage: R.<z> = B[]
        sage: f = z^2*(z-c) + c
        sage: polynomial_mandelbrot(f, pixel_count=100)
        100x100px 24-bit RGB image
    """

    cdef:
        z, c, x, y, w, cr, ci, J, M, S, S2, R, P, phi, t, im, re, a_n, df
        pt, c_pts, cp_real, cp_imag, pixel, color_list, critical_pts
        f_real, f_imag, f_temp, escape_time, cf_list, cf
        int i, j, d, k, col, row, iteration
        double C, L, Rad, x_corner, y_corner, x_coor, y_coor, \
            step_size, new_x, new_y
        I = CDF.gen()
        constant_c = True

    if parameter is None:
        c = var('c')
        parameter = c

    P = f.parent()

    if P.base_ring() is CC:
        if isinstance(P, FractionField_generic):
            raise NotImplementedError("coefficients must be polynomials in the parameter")
        gen_list = list(P.gens())
        parameter = gen_list.pop(gen_list.index(parameter))
        variable = gen_list.pop()

    elif P.base_ring().base_ring() is CC:
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

    # Make sure image_width is positive
    image_width = abs(image_width)

    # Initialize an image to the color black and access the pixels
    M = Image("RGB", (pixel_count, pixel_count), 'black')
    pixel = M.pixels()

    # Take the given base color and create a list of evenly spaced
    # colors between the given base color and white. The number of
    # colors in the list depends on the variable color_num.
    if isinstance(base_color, Color):
        # Convert Color to RGB list
        base_color = [int(k*255) for k in base_color]
    color_list = []
    for i in range(color_num):
        sig_check()
        color_list.append(copy(base_color))
        for j in range(3):
            color_list[i][j] += i * (255 - color_list[i][j]) // color_num
        color_list[i] = tuple(color_list[i])

    # Split function into real and imaginary parts
    R = PolynomialRing(CC, [variable, parameter])
    if len(R.gens()) > 2:
        raise NotImplementedError("base ring must have only 2 variables")
    z, c = R.gens()
    f = R(str(f))
    S = PolynomialRing(f.base_ring(), 'x,y,J,cr,ci')
    x, y, J, cr, ci = S.gens()
    S2 = S.quotient_ring(J**2+1)
    phi = R.hom([x+y*J, cr+ci*J], S2)
    t = phi(f).lift()
    im = t.coefficient(J)
    re = t - im*J

    f_real = fast_callable(re, vars=[x, y, cr, ci, J], domain=RDF)
    f_imag = fast_callable(im, vars=[x, y, cr, ci, J], domain=RDF)

    # Compute critical points
    try:
        df = f.derivative(z).univariate_polynomial()
        critical_pts = df.roots(multiplicities=False)
        constant_c = True
    except (PariError, TypeError):
        constant_c = False

    # If c is in the constant term of the polynomial, then the critical points
    # will be independent of c.
    if constant_c:
        c_pts = [[pt.real(), pt.imag()] for pt in critical_pts]

        # Calculate escape condition
        d = f.degree(z)
        cf_list = f.coefficients()
        a_n = cf_list.pop(0).abs()
        C = 0
        for cf in cf_list:
            C += cf.abs()
        L = 1.00000000000001
        if d >= 2:
            Rad = max(1, 2*C / a_n, (2*L / a_n**(1/(d-1))))
        else:
            Rad = max(1, 2*C / a_n)

        # First, we determine the complex coordinates of the point in the top
        # left corner of the image. Then, we loop through each pixel in the
        # image and assign it complex coordinates relative to the image's top
        # left corner.
        x_corner = x_center - image_width/2
        y_corner = y_center + image_width/2
        step_size = image_width*1.0 / pixel_count
        for col in range(pixel_count):
            x_coor = x_corner + col*step_size
            for row in range(pixel_count):
                sig_check()
                y_coor = y_corner - row*step_size

                # Initialize escape_time to be maximum number of iterations.
                escape_time = max_iteration

                # Loop though each critical point. If all of the critical
                # points are bounded for a particular c value, then c is in
                # the mandelbrot set.
                for pt in c_pts:
                    cp_real = pt[0]
                    cp_imag = pt[1]

                    # Iterate the map f using the critical point as the initial
                    # value.
                    new_x, new_y = cp_real, cp_imag
                    iteration = 0
                    while new_x**2 + new_y**2 <= Rad**2 and iteration < escape_time:
                        sig_check()
                        new_x, new_y = f_real(new_x, new_y, x_coor, y_coor, 1), \
                            f_imag(new_x, new_y, x_coor, y_coor, 1)
                        iteration += 1

                    # For each point, we take the minimum number of iterations
                    # over all the critical points and use this minimum to
                    # color the point.
                    if iteration < escape_time:
                        escape_time = iteration

                # If the point escapes to infinity, assign the point a color
                # based on how fast it escapes. The more iterations it takes for
                # a point to escape to infinity, the lighter its color will be.
                # Otherwise, assume the point is in the Mandelbrot set and leave
                # it black.
                if iteration != max_iteration:
                    # Assign each point a level based on its number of iterations.
                    level = iteration // level_sep
                    # Assign the pixel a color based on it's level. If we run out
                    # of colors, assign it the last color in the list.
                    if level < color_num:
                        pixel[col, row] = color_list[level]
                    else:
                        pixel[col, row] = color_list[-1]

    # If the critical points of f depend on c, we must compute the different
    # critical points for each c.
    else:
        # Solve for critical points symbollically.
        with SR.temp_var() as w:
            df = f.derivative(z).polynomial(z).subs({z: w})
            critical_pts = solve(symbolic_expression(df) == 0, w)
        c_pts = []
        for pt in critical_pts:
            c_pts.append(fast_callable(pt.rhs(), vars=[c], domain=CDF))

        # Calculate degree of f
        d = f.degree(z)

        # First, we determine the complex coordinates of the point in the top
        # left corner of the image. Then, we loop through each pixel in the
        # image and assign it complex coordinates relative to the image's top
        # left corner.
        x_corner = x_center - image_width/2
        y_corner = y_center + image_width/2
        step_size = image_width*1.0 / pixel_count
        for col in range(pixel_count):
            x_coor = x_corner + col*step_size
            for row in range(pixel_count):

                sig_check()
                y_coor = y_corner - row*step_size

                # Initialize escape time to be maximum number of iterations
                escape_time = max_iteration

                # Calculate escape condition for each c value
                f_temp = f.subs({c: x_coor+y_coor*I})
                cf_list = f_temp.coefficients()
                a_n = cf_list.pop(0).abs()
                C = 0
                for cf in cf_list:
                    C += cf.abs()
                L = 1.00000000000001
                if d >= 2:
                    Rad = max(1, 2*C / a_n, (2*L / a_n**(1/(d-1))))
                else:
                    Rad = max(1, 2*C / a_n)

                for f_cp in c_pts:

                    # Compute real and imaginary critical point
                    cp_real = f_cp(x_coor + y_coor*I).real()
                    cp_imag = f_cp(x_coor + y_coor*I).imag()

                    # Iterate the map f using the critical point as the initial
                    # value.
                    new_x, new_y = cp_real, cp_imag
                    iteration = 0
                    while new_x**2 + new_y**2 <= Rad**2 and iteration < escape_time:
                        sig_check()
                        new_x, new_y = f_real(new_x, new_y, x_coor, y_coor, 1), \
                            f_imag(new_x, new_y, x_coor, y_coor, 1)
                        iteration += 1

                    # For each point, we take the minimum number of iterations
                    # over all the critical points and use this minimum to
                    # color the point.
                    if iteration < escape_time:
                        escape_time = iteration

                # If the point escapes to infinity, assign the point a color
                # based on how fast it escapes. The more iterations it takes for
                # a point to escape to infinity, the lighter its color will be.
                # Otherwise, assume the point is in the Mandelbrot set and leave
                # it black.
                if iteration != max_iteration:
                    # Assign each point a level based on its number of iterations.
                    level = iteration // level_sep
                    # Assign the pixel a color based on it's level. If we run out
                    # of colors, assign it the last color in the list.
                    if level < color_num:
                        pixel[col, row] = color_list[level]
                    else:
                        pixel[col, row] = color_list[-1]
    return M

cpdef general_julia(f, double x_center=0, double y_center=0, image_width=4,
                    int max_iteration=50, int pixel_count=500,
                    int level_sep=1, int color_num=30,
                    base_color=[50, 50, 50]):
    r"""
    Plot Julia sets for general polynomials.

    EXAMPLES::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import general_julia
        sage: from sage.plot.colors import Color
        sage: R.<z> = CC[]
        sage: f = z^3 - z + 1
        sage: general_julia(f)
        500x500px 24-bit RGB image

    ::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import general_julia
        sage: from sage.plot.colors import Color
        sage: R.<z> = CC[]
        sage: f = z^5 - 1
        sage: general_julia(f)
        500x500px 24-bit RGB image
    """

    cdef:
        M, pixel, color_list, z, a_n
        int i, j, d, k, col, row, iteration
        double C, L, Rad, x_corner, y_corner, x_coor, y_coor, step_size
        I = CDF.gen()

    # Make sure image_width is positive
    image_width = abs(image_width)

    # Initialize an image to the color black and access the pixels
    M = Image("RGB", (pixel_count, pixel_count), 'black')
    pixel = M.pixels()

    # Take the given base color and create a list of evenly spaced
    # colors between the given base color and white. The number of
    # colors in the list depends on the variable color_num.
    if isinstance(base_color, Color):
        # Convert Color to RGB list
        base_color = [int(k*255) for k in base_color]
    color_list = []
    for i in range(color_num):
        sig_check()
        color_list.append(copy(base_color))
        for j in range(3):
            color_list[i][j] += i * (255 - color_list[i][j]) // color_num
        color_list[i] = tuple(color_list[i])

    z = f.variables()[0]
    f_fast = fast_callable(f, vars=[z], domain=CDF)

    # Calculate escape condition for each c value
    d = f.degree(z)
    cf_list = f.coefficients(sparse=False)
    a_n = cf_list.pop(-1).abs()
    C = 0
    for cf in cf_list:
        C += cf.abs()
    L = 1.00000000000001
    if d >= 2:
        Rad = max(1, 2*C / a_n, (2*L / a_n**(1/(d-1))))
    else:
        Rad = max(1, 2*C / a_n)

    # First, we determine the complex coordinates of the point in the top left
    # corner of the image. Then, we loop through each pixel in the image and
    # assign it complex coordinates relative to the image's top left corner.
    x_corner = x_center - image_width/2
    y_corner = y_center + image_width/2
    step_size = image_width*1.0 / pixel_count
    for col in range(pixel_count):
        x_coor = x_corner + col*step_size
        for row in range(pixel_count):
            sig_check()
            y_coor = y_corner - row*step_size

            # We compute the orbit of c under the map f
            # until we either reach the maximum number of iterations
            # or find a point in the orbit with modulus greater than
            # some the escape condition (Rad)

            new_z = x_coor + y_coor*I
            iteration = 0
            while new_z.abs() <= Rad**2 and iteration < max_iteration:
                sig_check()
                new_z = f_fast(new_z)
                iteration += 1

            # If the point escapes to infinity, assign the point a color
            # based on how fast it escapes. The more iterations it takes for
            # a point to escape to infinity, the lighter its color will be.
            # Otherwise, assume the point is in the Mandelbrot set and leave
            # it black.
            if iteration != max_iteration:
                # Assign each point a level based on its number of iterations.
                level = iteration // level_sep
                # Assign the pixel a color based on it's level. If we run out
                # of colors, assign it the last color in the list.
                if level < color_num:
                    pixel[col, row] = color_list[level]
                else:
                    pixel[col, row] = color_list[-1]
    return M
