# cython: binding=True
r"""
Mandelbrot and Julia sets (Cython helper)

This is the helper file providing functionality for mandel_julia.py.

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

from __future__ import absolute_import, division
from sage.plot.colors import Color
from sage.repl.image import Image
from copy import copy
from cysignals.signals cimport sig_check
from sage.rings.complex_field import ComplexField
from sage.functions.log import exp, log
from sage.symbolic.constants import pi


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


def fast_mandelbrot_plot(double x_center, double y_center, double image_width,
                         long max_iteration, long pixel_count, long level_sep,
                         long color_num, base_color):
    r"""
    Plots the Mandelbrot set in the complex plane for the map `Q_c(z) = z^2 + c`.

    INPUT:

    - ``x_center`` -- double, real part of the center point in the complex plane.

    - ``y_center`` -- double, imaginary part of the center point in the complex plane.

    - ``image_width`` -- double, width of the image in the complex plane.

    - ``max_iteration`` -- long, maximum number of iterations the map `Q_c(z)` considered.

    - ``pixel_count`` -- long, side length of image in number of pixels.

    - ``level_sep`` -- long, number of iterations between each color level.

    - ``color_num`` -- long, number of colors used to plot image.

    - ``base_color`` -- list, RGB color used to determine the coloring of set.

    OUTPUT:

    24-bit RGB image of the Mandelbrot set in the complex plane

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
    cdef long i, j, col, row, level, color_value, iteration
    cdef double k, x_corner, y_corner, step_size, x_coor, y_coor, new_x, new_y

    # Make sure image_width is positive
    image_width = abs(image_width)

    # Initialize an image to the color black and access the pixels
    M = Image("RGB", (pixel_count,pixel_count), 'black')
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
                    pixel[col,row] = color_list[level]
                else:
                    pixel[col,row] = color_list[-1]
    return M

cpdef fast_external_ray(double theta, long D=30, long S=10, long R=100,
 long pixel_count=500, double image_width=4, long prec=300):
    r"""
    Returns a list of points that approximate the external ray for a given angle.

    INPUT:

    - ``theta`` -- double, angle between 0 and 1 inclusive.

    - ``D`` -- long (optional - default: ``25``) depth of the approximation. As ``D`` increases, the external ray gets closer to the boundary of the Mandelbrot set.

    - ``S`` -- long (optional - default: ``10``) sharpness of the approximation. Adjusts the number of points used to approximate the external ray (number of points is equal to ``S*D``).

    - ``R`` -- long (optional - default: ``100``) radial parameter. If ``R`` is sufficiently large, the external ray reaches enough close to infinity.

    - ``pixel_count`` -- long (optional - default: ``500``) side length of image in number of pixels.

    - ``image_width`` -- double (optional - default: ``4``) width of the image in the complex plane.

    - ``prec`` -- long (optional - default: ``300``) specifies the bits of precision used by the Complex Field when using Newton's method to compute points on the external ray.

    OUTPUT:

    List of tuples of Real Interval Field Elements

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
        c_0, r_m, t_m, temp_c, C_k, D_k, old_c, x, y, dist
        int k, j, t
        double difference, m
        double error = pixel_count * 0.0001

        double pixel_width = image_width / pixel_count

        # initialize list with c_0
        c_list = [CF(R*exp(2*PI*I*theta))]

    # Loop through each subinterval and approximate point on external ray.
    for k in range(1,D+1):
        for j in range(1,S+1):
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
        x,y = c_list[k].real(), c_list[k].imag()
        c_list[k] = (x, y)

    return c_list

cpdef convert_to_pixels(point_list, double x_0, double y_0, double width,
 long number_of_pixels):
    r"""
    Converts cartesian coordinates to pixels within a specified window.

    INPUT:

    - ``point_list`` -- list of tuples, points in cartesian coordinates.

    - ``x_0`` -- double, x-coordinate of the center of the image.

    - ``y_0`` -- double, y-coordinate of the center of the image.

    - ``width`` -- double, width of visible window in caresian coordinates.

    - ``number_of_pixels`` -- long, width of image in pixels.

    OUTPUT:

    List of tuples of integers representing pixels.

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
    Produces a list of pixel coordinates approximating a line from a starting
    point to an ending point using the Bresenham's Line Algorithm.

    REFERENCE:

    [Br2016]_

    INPUT:

    - ``start`` -- tuple, starting point of line.

    - ``end`` -- tuple, ending point of line.

    OUTPUT:

    List of tuples of integers approximating the line between two pixels.

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


def fast_julia_plot(double c_real, double c_imag,
                    double x_center, double y_center, double image_width,
                    long max_iteration, long pixel_count, long level_sep,
                    long color_num, base_color):
    r"""
    Plots the Julia set for a given `c` value in the complex plane for the map `Q_c(z) = z^2 + c`.

    INPUT:

    - ``c_real`` -- double, Real part of `c` value that determines Julia set.

    - ``c_imag`` -- double, Imaginary part of `c` value that determines Julia set.

    - ``x_center`` -- double, Real part of center point.

    - ``y_center`` -- double, Imaginary part of center point.

    - ``image_width`` -- double, width of image in the complex plane.

    - ``max_iteration`` -- long, maximum number of iterations the map ``Q_c(z)``.

    - ``pixel_count`` -- long, side length of image in number of pixels.

    - ``level_sep`` -- long, number of iterations between each color level.

    - ``color_num`` -- long, number of colors used to plot image.

    - ``base_color`` -- RGB color, color used to determine the coloring of set.

    OUTPUT:

    24-bit RGB image of the Julia set in the complex plane.

    EXAMPLES:

    Plot the Julia set for `c=-1+0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_julia_plot
        sage: fast_julia_plot(-1, 0, 0, 0, 4, 500, 200, 1, 20, [40, 40, 40])
        200x200px 24-bit RGB image
    """
    cdef long i, j, col, row, level, color_value, iteration
    cdef double k, x_corner, y_corner, step_size, x_coor, y_coor, new_x, new_y

    # Make sure image_width is positive
    image_width = abs(image_width)

    # Initialize an image to the color black and access the pixels
    J = Image("RGB", (pixel_count,pixel_count), 'black')
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
            # or find a point in the orbit with modulus greater than 2
            new_x, new_y = x_coor, y_coor
            iteration = 0
            while (new_x**2 + new_y**2 <= 4.0 and iteration < max_iteration):
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
                    Jp[col,row] = color_list[level]
                else:
                    Jp[col,row] = color_list[-1]

    return J


def julia_helper(double c_real, double c_imag,
                 double x_center, double y_center, double image_width,
                 long max_iteration, long pixel_count, long level_sep,
                 long color_num, base_color, point_color):
    r"""
    Helper function that returns the image of a Julia set for a given
    `c` value side by side with the Mandelbrot set with a point denoting
    the `c` value.

    INPUT:

    - ``c_real`` -- double, Real part of `c` value that determines Julia set.

    - ``c_imag`` -- double, Imaginary part of `c` value that determines Julia set.

    - ``x_center`` -- double, Real part of center point.

    - ``y_center`` -- double, Imaginary part of center point.

    - ``image_width`` -- double, width of image in the complex plane.

    - ``max_iteration`` -- long, maximum number of iterations the map ``Q_c(z)``.

    - ``pixel_count`` -- long, side length of image in number of pixels.

    - ``level_sep`` -- long, number of iterations between each color level.

    - ``color_num`` -- long, number of colors used to plot image.

    - ``base_color`` -- RGB, color used to determine the coloring of set.

    - ``point_color`` -- RGB color, color of the point `c` in the Mandelbrot set.

    OUTPUT:

    24-bit RGB image of the Julia and Mandelbrot sets in the complex plane.

    EXAMPLES:

    Plot the Julia set for `c=-1+0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import julia_helper
        sage: julia_helper(-1, 0, 0, 0, 4, 500, 200, 1, 20, [40, 40, 40], [255, 0, 0])
        401x200px 24-bit RGB image
    """
    cdef int i, j

    # Initialize the Julia set
    J = fast_julia_plot(c_real, c_imag, x_center, y_center, image_width,
     max_iteration, pixel_count, level_sep, color_num, base_color)
    Jp = J.pixels()

    # Initialize the image with Julia set on left side
    # Add white border between images
    G = Image("RGB", (2*pixel_count+1,pixel_count), 'white')
    Gp = G.pixels()
    for i in range(pixel_count):
        for j in range(pixel_count):
            Gp[i,j] = Jp[i,j]

    # Plot the Mandelbrot set on the right side
    M = fast_mandelbrot_plot(-1, 0, 4, 500, pixel_count, 1, 30, base_color)
    Mp = M.pixels()
    for i in range(pixel_count+1,2*pixel_count):
        for j in range(pixel_count):
            Gp[i,j] = Mp[int(i-pixel_count),j]

    point_color = _color_to_RGB(point_color)

    # Add a cross representing c-value to the Mandelbrot set.
    CP = convert_to_pixels([(c_real, c_imag)], -1, 0, 4, pixel_count)
    for i in range(-3,4):
        # Loop through x and y coordinates and check if they are in image
        if min(CP[0][0]+i, CP[0][1]) >= 0 and \
         max(CP[0][0]+i, CP[0][1]) < pixel_count:
            Gp[CP[0][0]+i+pixel_count+1, CP[0][1]] = tuple(point_color)
        if min(CP[0][0], CP[0][1]+i) >= 0 and \
         max(CP[0][0], CP[0][1]+i) < pixel_count:
            Gp[CP[0][0]+pixel_count+1, CP[0][1]+i] = tuple(point_color)

    return G
