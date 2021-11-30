"""
List Plots
"""

from sage.structure.element import is_Matrix
from sage.matrix.constructor import matrix
from sage.rings.real_double import RDF
from sage.misc.superseded import deprecation


def list_plot3d(v, interpolation_type='default', point_list=None, **kwds):
    r"""
    A 3-dimensional plot of a surface defined by the list `v`
    of points in 3-dimensional space.

    INPUT:

    - ``v`` - something that defines a set of points in 3 space:

      - a matrix

      - a list of 3-tuples

      - a list of lists (all of the same length) - this is treated the same as
        a matrix.

    OPTIONAL KEYWORDS:

    - ``interpolation_type`` - 'linear', 'clough' (CloughTocher2D), 'spline'

      'linear' will perform linear interpolation

      The option 'clough' will interpolate by using a piecewise cubic
      interpolating Bezier polynomial on each triangle, using a
      Clough-Tocher scheme.  The interpolant is guaranteed to be
      continuously differentiable.  The gradients of the interpolant
      are chosen so that the curvature of the interpolating surface is
      approximatively minimized.

      The option 'spline' interpolates using a bivariate B-spline.

      When v is a matrix the default is to use linear interpolation, when
      v is a list of points the default is 'clough'.

    - ``degree`` - an integer between 1 and 5, controls the degree of spline
      used for spline interpolation. For data that is highly oscillatory
      use higher values

    - ``point_list`` - If point_list=True is passed, then if the array
      is a list of lists of length three, it will be treated as an
      array of points rather than a 3xn array.

    - ``num_points`` - Number of points to sample interpolating
      function in each direction, when ``interpolation_type`` is not
      ``default``. By default for an `n\times n` array this is `n`.

    - ``**kwds`` - all other arguments are passed to the surface function

    OUTPUT: a 3d plot

    EXAMPLES:

    We plot a matrix that illustrates summation modulo `n`::

        sage: n = 5
        sage: list_plot3d(matrix(RDF, n, [(i+j)%n for i in [1..n] for j in [1..n]]))
        Graphics3d Object

    We plot a matrix of values of ``sin``::

        sage: pi = float(pi)
        sage: m = matrix(RDF, 6, [sin(i^2 + j^2) for i in [0,pi/5,..,pi] for j in [0,pi/5,..,pi]])
        sage: list_plot3d(m, color='yellow', frame_aspect_ratio=[1, 1, 1/3])
        Graphics3d Object

    Though it does not change the shape of the graph, increasing
    num_points can increase the clarity of the graph::

        sage: list_plot3d(m, color='yellow', frame_aspect_ratio=[1, 1, 1/3], num_points=40)
        Graphics3d Object

    We can change the interpolation type::

        sage: import warnings
        sage: warnings.simplefilter('ignore', UserWarning)
        sage: list_plot3d(m, color='yellow', interpolation_type='clough', frame_aspect_ratio=[1, 1, 1/3])
        Graphics3d Object

    We can make this look better by increasing the number of samples::

        sage: list_plot3d(m, color='yellow', interpolation_type='clough', frame_aspect_ratio=[1, 1, 1/3], num_points=40)
        Graphics3d Object

    Let us try a spline::

        sage: list_plot3d(m, color='yellow', interpolation_type='spline', frame_aspect_ratio=[1, 1, 1/3])
        Graphics3d Object

    That spline does not capture the oscillation very well; let's try a
    higher degree spline::

        sage: list_plot3d(m, color='yellow', interpolation_type='spline', degree=5, frame_aspect_ratio=[1, 1, 1/3])
        Graphics3d Object

    We plot a list of lists::

        sage: show(list_plot3d([[1, 1, 1, 1], [1, 2, 1, 2], [1, 1, 3, 1], [1, 2, 1, 4]]))

    We plot a list of points.  As a first example we can extract the
    (x,y,z) coordinates from the above example and make a list plot
    out of it. By default we do linear interpolation::

        sage: l = []
        sage: for i in range(6):
        ....:     for j in range(6):
        ....:         l.append((float(i*pi/5), float(j*pi/5), m[i, j]))
        sage: list_plot3d(l, color='red')
        Graphics3d Object

    Note that the points do not have to be regularly sampled. For example::

        sage: l = []
        sage: for i in range(-5, 5):
        ....:     for j in range(-5, 5):
        ....:         l.append((normalvariate(0, 1), normalvariate(0, 1), normalvariate(0, 1)))
        sage: L = list_plot3d(l, interpolation_type='clough', color='orange', num_points=100)
        sage: L
        Graphics3d Object

    Check that no NaNs are produced (see :trac:`13135`)::

        sage: any(math.isnan(c) for v in L.vertices() for c in v)
        False

    TESTS:

    We plot 0, 1, and 2 points::

        sage: list_plot3d([])
        Graphics3d Object

        sage: list_plot3d([(2, 3, 4)])
        Graphics3d Object

        sage: list_plot3d([(0, 0, 1), (2, 3, 4)])
        Graphics3d Object

    However, if two points are given with the same x,y coordinates but
    different z coordinates, an exception will be raised::

        sage: pts = [(-4/5, -2/5, -2/5), (-4/5, -2/5, 2/5), (-4/5, 2/5, -2/5), (-4/5, 2/5, 2/5), (-2/5, -4/5, -2/5), (-2/5, -4/5, 2/5), (-2/5, -2/5, -4/5), (-2/5, -2/5, 4/5), (-2/5, 2/5, -4/5), (-2/5, 2/5, 4/5), (-2/5, 4/5, -2/5), (-2/5, 4/5, 2/5), (2/5, -4/5, -2/5), (2/5, -4/5, 2/5), (2/5, -2/5, -4/5), (2/5, -2/5, 4/5), (2/5, 2/5, -4/5), (2/5, 2/5, 4/5), (2/5, 4/5, -2/5), (2/5, 4/5, 2/5), (4/5, -2/5, -2/5), (4/5, -2/5, 2/5), (4/5, 2/5, -2/5), (4/5, 2/5, 2/5)]
        sage: show(list_plot3d(pts, interpolation_type='clough'))
        Traceback (most recent call last):
        ...
        ValueError: points with same x,y coordinates and different
        z coordinates were given. Interpolation cannot handle this.

    Additionally we need at least 3 points to do the interpolation::

        sage: mat = matrix(RDF, 1, 2, [3.2, 1.550])
        sage: show(list_plot3d(mat, interpolation_type='clough'))
        Traceback (most recent call last):
        ...
        ValueError: we need at least 3 points to perform the interpolation

    TESTS::

        sage: P = list_plot3d([(0, 0, 1), (2, 3, 4)], texture='tomato')
        doctest:warning...:
        DeprecationWarning: please use 'color' instead of 'texture'
        See https://trac.sagemath.org/27084 for details.
    """
    import numpy
    if 'texture' in kwds:
        deprecation(27084, "please use 'color' instead of 'texture'")
        txtr = kwds.pop('texture')
        if txtr == "automatic":
            txtr = "lightblue"
        kwds['color'] = txtr
    if is_Matrix(v):
        if (interpolation_type == 'default' or
            interpolation_type == 'linear' and 'num_points' not in kwds):
            return list_plot3d_matrix(v, **kwds)
        else:
            data = [(i, j, v[i, j])
                    for i in range(v.nrows())
                    for j in range(v.ncols())]
            return list_plot3d_tuples(data, interpolation_type, **kwds)

    if isinstance(v, numpy.ndarray):
        return list_plot3d(matrix(v), interpolation_type, **kwds)

    if isinstance(v, list):
        if not v:
            # return empty 3d graphic
            from .base import Graphics3d
            return Graphics3d()
        elif len(v) == 1:
            # return a point
            from .shapes2 import point3d
            return point3d(v[0], **kwds)
        elif len(v) == 2:
            # return a line
            from .shapes2 import line3d
            return line3d(v, **kwds)
        elif isinstance(v[0], tuple) or point_list and len(v[0]) == 3:
            return list_plot3d_tuples(v, interpolation_type, **kwds)
        else:
            return list_plot3d_array_of_arrays(v, interpolation_type, **kwds)
    raise TypeError("v must be a matrix or list")


def list_plot3d_matrix(m, **kwds):
    """
    A 3-dimensional plot of a surface defined by a matrix ``M``
    defining points in 3-dimensional space.

    See :func:`list_plot3d` for full details.

    INPUT:

    - ``M`` -- a matrix

    OPTIONAL KEYWORDS:

    - ``**kwds`` - all other arguments are passed to the surface function

    OUTPUT: a 3d plot

    EXAMPLES:

    We plot a matrix that illustrates summation modulo `n`::

        sage: n = 5
        sage: list_plot3d(matrix(RDF, n, [(i+j)%n for i in [1..n] for j in [1..n]])) # indirect doctest
        Graphics3d Object

    The interpolation type for matrices is 'linear'; for other types
    use other :func:`list_plot3d` input types.

    We plot a matrix of values of `sin`::

        sage: pi = float(pi)
        sage: m = matrix(RDF, 6, [sin(i^2 + j^2) for i in [0,pi/5,..,pi] for j in [0,pi/5,..,pi]])
        sage: list_plot3d(m, color='yellow', frame_aspect_ratio=[1, 1, 1/3]) # indirect doctest
        Graphics3d Object
        sage: list_plot3d(m, color='yellow', interpolation_type='linear') # indirect doctest
        Graphics3d Object

    Here is a colored example, using a colormap and a coloring function
    which must take values in (0, 1)::

        sage: cm = colormaps.rainbow
        sage: n = 20
        sage: cf = lambda x, y: ((2*(x-y)/n)**2) % 1
        sage: list_plot3d(matrix(RDF, n, [cos(pi*(i+j)/n) for i in [1..n]
        ....:   for j in [1..n]]), color=(cf,cm))
        Graphics3d Object

    .. PLOT::

        cm = colormaps.rainbow
        cf = lambda x, y: ((2*(x-y)/20)**2) % 1
        expl = list_plot3d(matrix(RDF,20,20,[cos(pi*(i+j)/20) for i in range(1,21) for j in range(1,21)]),color=(cf,cm))
        sphinx_plot(expl)
    """
    from .parametric_surface import ParametricSurface

    def f(i, j):
        return (i, j, float(m[int(i), int(j)]))
    G = ParametricSurface(f, (list(range(m.nrows())), list(range(m.ncols()))),
                          **kwds)
    G._set_extra_kwds(kwds)
    return G


def list_plot3d_array_of_arrays(v, interpolation_type, **kwds):
    """
    A 3-dimensional plot of a surface defined by a list of lists ``v``
    defining points in 3-dimensional space.

    This is done by making the list of lists into a matrix and passing
    back to :func:`list_plot3d`.  See :func:`list_plot3d` for full
    details.

    INPUT:

    - ``v`` - a list of lists, all the same length
    - ``interpolation_type`` - (default: 'linear')

    OPTIONAL KEYWORDS:

    - ``**kwds`` - all other arguments are passed to the surface function

    OUTPUT: a 3d plot

    EXAMPLES:

    The resulting matrix does not have to be square::

        sage: show(list_plot3d([[1, 1, 1, 1], [1, 2, 1, 2], [1, 1, 3, 1]])) # indirect doctest

    The normal route is for the list of lists to be turned into a matrix
    and use :func:`list_plot3d_matrix`::

        sage: show(list_plot3d([[1, 1, 1, 1], [1, 2, 1, 2], [1, 1, 3, 1], [1, 2, 1, 4]]))

    With certain extra keywords (see :func:`list_plot3d_matrix`), this function
    will end up using :func:`list_plot3d_tuples`::

        sage: show(list_plot3d([[1, 1, 1, 1], [1, 2, 1, 2], [1, 1, 3, 1], [1, 2, 1, 4]], interpolation_type='spline'))
    """
    m = matrix(RDF, len(v), len(v[0]), v)
    G = list_plot3d(m, interpolation_type, **kwds)
    G._set_extra_kwds(kwds)
    return G


def list_plot3d_tuples(v, interpolation_type, **kwds):
    r"""
    A 3-dimensional plot of a surface defined by the list `v`
    of points in 3-dimensional space.

    INPUT:

    - ``v`` - something that defines a set of points in 3
      space, for example:

      - a matrix

        This will be if using an interpolation type other than
        'linear', or if using ``num_points`` with 'linear'; otherwise
        see :func:`list_plot3d_matrix`.

      - a list of 3-tuples

      - a list of lists (all of the same length, under same conditions
        as a matrix)

    OPTIONAL KEYWORDS:

    - ``interpolation_type`` - 'linear', 'clough' (CloughTocher2D), 'spline'

      'linear' will perform linear interpolation

      The option 'clough' will interpolate by using a piecewise cubic
      interpolating Bezier polynomial on each triangle, using a
      Clough-Tocher scheme.  The interpolant is guaranteed to be
      continuously differentiable.

      The option 'spline' interpolates using a bivariate B-spline.

      When v is a matrix the default is to use linear interpolation, when
      v is a list of points the default is 'clough'.

    - ``degree`` - an integer between 1 and 5, controls the degree of spline
      used for spline interpolation. For data that is highly oscillatory
      use higher values

    - ``point_list`` - If point_list=True is passed, then if the array
      is a list of lists of length three, it will be treated as an
      array of points rather than a `3\times n` array.

    - ``num_points`` - Number of points to sample interpolating
      function in each direction.  By default for an `n\times n`
      array this is `n`.

    - ``**kwds`` - all other arguments are passed to the
      surface function

    OUTPUT: a 3d plot

    EXAMPLES:

    All of these use this function; see :func:`list_plot3d` for other
    list plots::

        sage: pi = float(pi)
        sage: m = matrix(RDF, 6, [sin(i^2 + j^2) for i in [0,pi/5,..,pi] for j in [0,pi/5,..,pi]])
        sage: list_plot3d(m, color='yellow', interpolation_type='linear', num_points=5) # indirect doctest
        Graphics3d Object

    ::

        sage: list_plot3d(m, color='yellow', interpolation_type='spline', frame_aspect_ratio=[1, 1, 1/3])
        Graphics3d Object

    ::

        sage: show(list_plot3d([[1, 1, 1], [1, 2, 1], [0, 1, 3], [1, 0, 4]], point_list=True))

    ::

        sage: list_plot3d([(1, 2, 3), (0, 1, 3), (2, 1, 4), (1, 0, -2)], color='yellow', num_points=50)  # long time
        Graphics3d Object
    """
    from matplotlib import tri
    import numpy
    from random import random
    from scipy import interpolate
    from .plot3d import plot3d

    if len(v) < 3:
        raise ValueError("we need at least 3 points to perform the "
                         "interpolation")

    x = [float(p[0]) for p in v]
    y = [float(p[1]) for p in v]
    z = [float(p[2]) for p in v]

    # If the (x,y)-coordinates lie in a one-dimensional subspace, the
    # matplotlib Delaunay code segfaults.  Therefore, we compute the
    # correlation of the x- and y-coordinates and add small random
    # noise to avoid the problem if needed.
    corr_matrix = numpy.corrcoef(x, y)
    if not (-0.9 <= corr_matrix[0, 1] <= 0.9):
        ep = float(.000001)
        x = [float(p[0]) + random() * ep for p in v]
        y = [float(p[1]) + random() * ep for p in v]

    # If the list of data points has two points with the exact same
    # (x,y)-coordinate but different z-coordinates, then we sometimes
    # get segfaults.  The following block checks for this and raises
    # an exception if this is the case.
    # We also remove duplicate points (which matplotlib can't handle).
    # Alternatively, the code in the if block above which adds random
    # error could be applied to perturb the points.
    drop_list = []
    nb_points = len(x)
    for i in range(nb_points):
        for j in range(i + 1, nb_points):
            if x[i] == x[j] and y[i] == y[j]:
                if z[i] != z[j]:
                    raise ValueError("points with same x,y coordinates"
                                     " and different z coordinates were"
                                     " given. Interpolation cannot handle this.")
                elif z[i] == z[j]:
                    drop_list.append(j)
    x = [x[i] for i in range(nb_points) if i not in drop_list]
    y = [y[i] for i in range(nb_points) if i not in drop_list]
    z = [z[i] for i in range(nb_points) if i not in drop_list]

    xmin = float(min(x))
    xmax = float(max(x))
    ymin = float(min(y))
    ymax = float(max(y))

    num_points = kwds.get('num_points', int(4 * numpy.sqrt(len(x))))
    # arbitrary choice - assuming more or less a n x n grid of points
    # x should have n^2 entries. We sample 4 times that many points.

    if interpolation_type == 'linear':
        T = tri.Triangulation(x, y)
        f = tri.LinearTriInterpolator(T, z)
        j = complex(0, 1)
        from .parametric_surface import ParametricSurface

        def g(x, y):
            z = f(x, y)
            return (x, y, z)
        G = ParametricSurface(g, (list(numpy.r_[xmin:xmax:num_points * j]),
                                  list(numpy.r_[ymin:ymax:num_points * j])),
                              **kwds)
        G._set_extra_kwds(kwds)
        return G

    if interpolation_type == 'clough' or interpolation_type == 'default':

        points = [[x[i], y[i]] for i in range(len(x))]
        j = complex(0, 1)
        f = interpolate.CloughTocher2DInterpolator(points, z)
        from .parametric_surface import ParametricSurface

        def g(x, y):
            z = f([x, y])
            return (x, y, z)
        G = ParametricSurface(g, (list(numpy.r_[xmin:xmax:num_points * j]),
                                  list(numpy.r_[ymin:ymax:num_points * j])),
                              **kwds)
        G._set_extra_kwds(kwds)
        return G

    if interpolation_type == 'spline':
        kx = kwds['kx'] if 'kx' in kwds else 3
        ky = kwds['ky'] if 'ky' in kwds else 3
        if 'degree' in kwds:
            kx = kwds['degree']
            ky = kwds['degree']
        s = kwds.get('smoothing', len(x) - numpy.sqrt(2 * len(x)))
        s = interpolate.bisplrep(x, y, z, [int(1)] * len(x), xmin, xmax,
                                 ymin, ymax, kx=kx, ky=ky, s=s)

        def f(x, y):
            return interpolate.bisplev(x, y, s)
        return plot3d(f, (xmin, xmax), (ymin, ymax),
                      plot_points=[num_points, num_points], **kwds)
