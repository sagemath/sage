"""
Transformations
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.math cimport sin, cos, sqrt

include "point_c.pxi"

from sage.rings.real_double import RDF

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

pi = RDF.pi()


cdef class Transformation:
    def __init__(self, scale=(1,1,1),
                       rot=None,
                       trans=[0,0,0],
                       m=None):

        if scale is None:
            scale = (1,1,1)
        if trans is None:
            trans = [0,0,0]

        # TODO: determine for sure if x3d does scale or rotation first
        if m is not None:
            self.matrix = m

        else:
            m = matrix(RDF, 3, 3,
                      [scale[0], 0, 0, 0, scale[1], 0, 0, 0, scale[2]])

            if rot is not None:
                # rotate about v by theta
                vx, vy, vz, theta = rot
                m *= rotate_arbitrary((vx, vy, vz), theta)

            self.matrix = m.augment(matrix(RDF, 3, 1, list(trans))) \
                           .stack(matrix(RDF, 1, 4, [0,0,0,1]))

        # this raw data is used for optimized transformations
        m_data = self.matrix.list()
        cdef int i
        for i in range(12):
            self._matrix_data[i] = m_data[i]

    def get_matrix(self):
        return self.matrix.__copy__()

    def is_skew(self, eps=1e-5):
        dx, dy, dz = self.matrix[0:3, 0:3].columns()
        return abs(dx.dot_product(dy)) + abs(dx.dot_product(dz)) + abs(dy.dot_product(dz)) > eps

    def is_uniform(self, eps=1e-5):
        cols = self.matrix[0:3, 0:3].columns()
        lens = [col.dot_product(col) for col in cols]
        return abs(lens[0] - lens[1]) + abs(lens[0] - lens[2]) < eps

    def is_uniform_on(self, basis, eps=1e-5):
        basis = [vector(RDF, self.transform_vector(b)) for b in basis]
        a = basis.pop()
        len_a = a.dot_product(a)
        return max([abs(len_a - b.dot_product(b)) for b in basis]) < eps

    cpdef transform_point(self, x):
        cdef point_c res, P
        P.x, P.y, P.z = x
        point_c_transform(&res, self._matrix_data, P)
        return res.x, res.y, res.z

    cpdef transform_vector(self, v):
        cdef point_c res, P
        P.x, P.y, P.z = v
        point_c_stretch(&res, self._matrix_data, P)
        return res.x, res.y, res.z

    cpdef transform_bounding_box(self, box):
        cdef point_c lower, upper, res, temp
        cdef point_c bounds[2]
        bounds[0].x, bounds[0].y, bounds[0].z = box[0]
        bounds[1].x, bounds[1].y, bounds[1].z = box[1]
        point_c_transform(&lower, self._matrix_data, bounds[0])
        point_c_transform(&upper, self._matrix_data, bounds[0])
        cdef int i
        for i in range(1, 8):
            temp.x = bounds[ i & 1      ].x
            temp.y = bounds[(i & 2) >> 1].y
            temp.z = bounds[(i & 4) >> 2].z
            point_c_transform(&res, self._matrix_data, temp)
            point_c_lower_bound(&lower, lower, res)
            point_c_upper_bound(&upper, upper, res)
        return (lower.x, lower.y, lower.z), (upper.x, upper.y, upper.z)

    cdef void transform_point_c(self, point_c* res, point_c P) noexcept:
        point_c_transform(res, self._matrix_data, P)

    cdef void transform_vector_c(self, point_c* res, point_c P) noexcept:
        point_c_stretch(res, self._matrix_data, P)

    def __mul__(Transformation self, Transformation other):
        return Transformation(m = self.matrix * other.matrix)

    def __invert__(Transformation self):
        return Transformation(m=~self.matrix)

    def __call__(self, p):
        return self.transform_point(p)

    def max_scale(self):
        if self._svd is None:
            self._svd = self.matrix[0:3, 0:3].SVD()
        return self._svd[1][0,0]

    def avg_scale(self):
        if self._svd is None:
            self._svd = self.matrix[0:3, 0:3].SVD()
        return (self._svd[1][0,0] * self._svd[1][1,1] * self._svd[1][2,2]) ** (1/3.0)


def rotate_arbitrary(v, double theta):
    """
    Return a matrix that rotates the coordinate space about
    the axis v by the angle ``theta``.

    INPUT:

    - ``theta`` -- real number, the angle

    EXAMPLES::

        sage: from sage.plot.plot3d.transform import rotate_arbitrary

    Try rotating about the axes::

        sage: rotate_arbitrary((1,0,0), 1)
        [                1.0                 0.0                 0.0]
        [                0.0  0.5403023058681398  0.8414709848078965]
        [                0.0 -0.8414709848078965  0.5403023058681398]
        sage: rotate_arbitrary((0,1,0), 1)
        [ 0.5403023058681398                 0.0 -0.8414709848078965]
        [                0.0                 1.0                 0.0]
        [ 0.8414709848078965                 0.0  0.5403023058681398]
        sage: rotate_arbitrary((0,0,1), 1)
        [ 0.5403023058681398  0.8414709848078965                 0.0]
        [-0.8414709848078965  0.5403023058681398                 0.0]
        [                0.0                 0.0                 1.0]

    These next two should be the same (up to floating-point errors)::

        sage: rotate_arbitrary((1,1,1), 1)  # rel tol 1e-15
        [  0.6935348705787598   0.6390560643047186 -0.33259093488347846]
        [-0.33259093488347846   0.6935348705787598   0.6390560643047186]
        [  0.6390560643047186  -0.3325909348834784   0.6935348705787598]
        sage: rotate_arbitrary((1,1,1), -1)^(-1)  # rel tol 1e-15
        [  0.6935348705787598   0.6390560643047186 -0.33259093488347846]
        [-0.33259093488347846   0.6935348705787598   0.6390560643047186]
        [  0.6390560643047185 -0.33259093488347835   0.6935348705787598]

    Make sure it does the right thing...::

        sage: rotate_arbitrary((1,2,3), -1).det()
        1.0000000000000002
        sage: rotate_arbitrary((1,1,1), 2*pi/3) * vector(RDF, (1,2,3))  # rel tol 2e-15             # needs sage.symbolic
        (1.9999999999999996, 2.9999999999999996, 0.9999999999999999)
        sage: rotate_arbitrary((1,2,3), 5) * vector(RDF, (1,2,3))  # rel tol 2e-15
        (1.0000000000000002, 2.0, 3.000000000000001)
        sage: rotate_arbitrary((1,1,1), pi/7)^7  # rel tol 2e-15                        # needs sage.symbolic
        [-0.33333333333333337   0.6666666666666671   0.6666666666666665]
        [  0.6666666666666665 -0.33333333333333337   0.6666666666666671]
        [  0.6666666666666671   0.6666666666666667 -0.33333333333333326]

    AUTHORS:

       - Robert Bradshaw

    ALGORITHM:

        There is a formula. Where did it come from? Lets take
        a quick jaunt into Sage's calculus package...

        Setup some variables::

            sage: vx,vy,vz,theta = var('x y z theta')                                   # needs sage.symbolic

        Symbolic rotation matrices about X and Y axis::

            sage: def rotX(theta): return matrix(SR, 3, 3, [1, 0, 0,  0, cos(theta), -sin(theta), 0, sin(theta), cos(theta)])
            sage: def rotZ(theta): return matrix(SR, 3, 3, [cos(theta), -sin(theta), 0,  sin(theta), cos(theta), 0, 0, 0, 1])

        Normalizing `y` so that `|v|=1`. Perhaps there is a better
        way to tell Maxima that `x^2+y^2+z^2=1` which would make for
        a much cleaner calculation::

            sage: vy = sqrt(1-vx^2-vz^2)                                                # needs sage.symbolic

        Now we rotate about the `x`-axis so `v` is in the `xy`-plane::

            sage: t = arctan(vy/vz)+pi/2                                                # needs sage.symbolic
            sage: m = rotX(t)                                                           # needs sage.symbolic
            sage: new_y = vy*cos(t) - vz*sin(t)                                         # needs sage.symbolic

        And rotate about the `z` axis so `v` lies on the `x` axis::

            sage: s = arctan(vx/new_y) + pi/2                                           # needs sage.symbolic
            sage: m = rotZ(s) * m                                                       # needs sage.symbolic

        Rotating about `v` in our old system is the same as rotating
        about the `x`-axis in the new::

            sage: m = rotX(theta) * m                                                   # needs sage.symbolic

        Do some simplifying here to avoid blow-up::

            sage: m = m.simplify_rational()                                             # needs sage.symbolic

        Now go back to the original coordinate system::

            sage: m = rotZ(-s) * m                                                      # needs sage.symbolic
            sage: m = rotX(-t) * m                                                      # needs sage.symbolic

        And simplify every single entry (which is more effective that simplify
        the whole matrix like above)::

            sage: m.parent()([x.simplify_full() for x in m._list()])  # random  # long time, needs sage.symbolic
            [                                       -(cos(theta) - 1)*x^2 + cos(theta)              -(cos(theta) - 1)*sqrt(-x^2 - z^2 + 1)*x + sin(theta)*abs(z)      -((cos(theta) - 1)*x*z^2 + sqrt(-x^2 - z^2 + 1)*sin(theta)*abs(z))/z]
            [             -(cos(theta) - 1)*sqrt(-x^2 - z^2 + 1)*x - sin(theta)*abs(z)                           (cos(theta) - 1)*x^2 + (cos(theta) - 1)*z^2 + 1 -((cos(theta) - 1)*sqrt(-x^2 - z^2 + 1)*z*abs(z) - x*z*sin(theta))/abs(z)]
            [     -((cos(theta) - 1)*x*z^2 - sqrt(-x^2 - z^2 + 1)*sin(theta)*abs(z))/z -((cos(theta) - 1)*sqrt(-x^2 - z^2 + 1)*z*abs(z) + x*z*sin(theta))/abs(z)                                        -(cos(theta) - 1)*z^2 + cos(theta)]

        Re-expressing some entries in terms of y and resolving the absolute
        values introduced by eliminating y, we get the desired result.
    """
    cdef double x,y,z, len_v
    x,y,z = v
    len_v = sqrt(x*x+y*y+z*z)
    # normalize for an easier formula
    x /= len_v
    y /= len_v
    z /= len_v
    cdef double cos_t = cos(theta), sin_t = sin(theta)

    entries = [
        (1 - cos_t)*x*x + cos_t,
        sin_t*z - (cos_t - 1)*x*y,
       -sin_t*y + (1 - cos_t)*x*z,

       -sin_t*z + (1 - cos_t)*x*y,
        (1 - cos_t)*y*y + cos_t,
        sin_t*x - (cos_t - 1)*z*y,

        sin_t*y - (cos_t - 1)*x*z,
       -(cos_t - 1)*z*y - sin_t*x,
        (1 - cos_t)*z*z + cos_t        ]

    return matrix(RDF, 3, 3, entries)
