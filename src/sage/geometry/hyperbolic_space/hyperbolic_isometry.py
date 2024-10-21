# sage_setup: distribution = sagemath-symbolics
# sage.doctest: needs scipy
r"""
Hyperbolic Isometries

This module implements the abstract base class for isometries of
hyperbolic space of arbitrary dimension.  It also contains the
implementations for specific models of hyperbolic geometry.

The isometry groups of all implemented models are either matrix Lie
groups or are doubly covered by matrix Lie groups.  As such, the
isometry constructor takes a matrix as input.  However, since the
isometries themselves may not be matrices, quantities like the trace
and determinant are not directly accessible from this class.

AUTHORS:

- Greg Laun (2013): initial version

EXAMPLES:

We can construct isometries in the upper half plane model, abbreviated
UHP for convenience::

    sage: UHP = HyperbolicPlane().UHP()
    sage: UHP.get_isometry(matrix(2,[1,2,3,4]))
    Isometry in UHP
    [1 2]
    [3 4]
    sage: A = UHP.get_isometry(matrix(2,[0,1,1,0]))
    sage: A.inverse()
    Isometry in UHP
    [0 1]
    [1 0]
"""

# **********************************************************************
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# **********************************************************************

from copy import copy
from sage.categories.homset import Hom
from sage.categories.morphism import Morphism
from sage.misc.lazy_attribute import lazy_attribute
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.infinity import infinity
from sage.misc.latex import latex
from sage.rings.real_double import RDF
from sage.functions.other import imag
from sage.misc.functional import sqrt
from sage.functions.all import arccosh, sign

from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON
from sage.geometry.hyperbolic_space.hyperbolic_geodesic import HyperbolicGeodesic


class HyperbolicIsometry(Morphism):
    r"""
    Abstract base class for hyperbolic isometries.  This class should
    never be instantiated.

    INPUT:

    - ``A`` -- a matrix representing a hyperbolic isometry in the
      appropriate model

    EXAMPLES::

        sage: HyperbolicPlane().HM().get_isometry(identity_matrix(3))
        Isometry in HM
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """

    #####################
    # "Private" Methods #
    #####################

    def __init__(self, model, A, check=True):
        r"""
        See :class:`HyperbolicIsometry` for full documentation.

        EXAMPLES::

            sage: A = HyperbolicPlane().UHP().get_isometry(matrix(2, [0,1,-1,0]))
            sage: TestSuite(A).run(skip='_test_category')
        """
        if check:
            model.isometry_test(A)
        self._matrix = copy(A) # Make a copy of the potentially mutable matrix
        self._matrix.set_immutable() # Make it immutable
        Morphism.__init__(self, Hom(model, model))

    @lazy_attribute
    def _cached_isometry(self):
        r"""
        The representation of the current isometry used for
        calculations.  For example, if the current model uses the
        upper half plane, then ``_cached_isometry`` will
        hold the `SL(2,\RR)` representation of ``self.matrix()``.

        EXAMPLES::

            sage: A = HyperbolicPlane().HM().get_isometry(identity_matrix(3))
            sage: A._cached_isometry
            Isometry in UHP
            [1 0]
            [0 1]
        """
        R = self.domain().realization_of().a_realization()
        return self.to_model(R)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_isometry(identity_matrix(2))
            Isometry in UHP
            [1 0]
            [0 1]
        """
        return self._repr_type() + " in {0}\n{1}".format(self.domain().short_name(), self._matrix)

    def _repr_type(self):
        r"""
        Return the type of morphism.

        EXAMPLES::

            sage: A = HyperbolicPlane().UHP().get_isometry(identity_matrix(2))
            sage: A._repr_type()
            'Isometry'
        """
        return "Isometry"

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: A = HyperbolicPlane().UHP().get_isometry(identity_matrix(2))
            sage: latex(A)
            \pm \left(\begin{array}{rr}
            1 & 0 \\
            0 & 1
            \end{array}\right)

            sage: B = HyperbolicPlane().HM().get_isometry(identity_matrix(3))
            sage: latex(B)
            \left(\begin{array}{rrr}
            1 & 0 & 0 \\
            0 & 1 & 0 \\
            0 & 0 & 1
            \end{array}\right)
        """
        if self.domain().is_isometry_group_projective():
            return r"\pm " + latex(self._matrix)
        else:
            return latex(self._matrix)

    def __eq__(self, other):
        r"""
        Return ``True`` if the isometries are the same and ``False`` otherwise.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = UHP.get_isometry(identity_matrix(2))
            sage: B = UHP.get_isometry(-identity_matrix(2))
            sage: A == B
            True

            sage: HM = HyperbolicPlane().HM()
            sage: A = HM.random_isometry()
            sage: A == A
            True
        """
        if not isinstance(other, HyperbolicIsometry):
            return False
        test_matrix = bool((self.matrix() - other.matrix()).norm() < EPSILON)
        if self.domain().is_isometry_group_projective():
            A,B = self.matrix(), other.matrix() # Rename for simplicity
            m = self.matrix().ncols()
            A = A / sqrt(A.det(), m) # Normalized to have determinant 1
            B = B / sqrt(B.det(), m)
            test_matrix = ((A - B).norm() < EPSILON
                           or (A + B).norm() < EPSILON)
        return self.domain() is other.domain() and test_matrix

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = UHP.get_isometry(identity_matrix(2))
            sage: B = UHP.get_isometry(-identity_matrix(2))
            sage: hash(A) == hash(B)
            True

            sage: HM = HyperbolicPlane().HM()
            sage: A = HM.random_isometry()
            sage: hash(A) == hash(A)
            True
        """
        if self.domain().is_isometry_group_projective():
            # Special care must be taken for projective groups
            m = matrix(self._matrix.nrows(),
                       [abs(x) for x in self._matrix.list()])
            m.set_immutable()
        else:
            m = self._matrix
        return hash((self.domain(), self.codomain(), m))

    def __pow__(self, n):
        r"""
        EXAMPLES::

            sage: A = HyperbolicPlane().UHP().get_isometry(matrix(2,[3,1,2,1]))
            sage: A**3
            Isometry in UHP
            [41 15]
            [30 11]
        """
        return self.__class__(self.domain(), self._matrix**n)

    def __mul__(self, other):
        r"""
        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = UHP.get_isometry(Matrix(2,[5,2,1,2]))
            sage: B = UHP.get_isometry(Matrix(2,[3,1,1,2]))
            sage: B * A
            Isometry in UHP
            [16  8]
            [ 7  6]
            sage: A = UHP.get_isometry(Matrix(2,[5,2,1,2]))
            sage: p = UHP.get_point(2 + I)
            sage: A * p
            Point in UHP 8/17*I + 53/17

            sage: g = UHP.get_geodesic(2 + I, 4 + I)
            sage: A * g
            Geodesic in UHP from 8/17*I + 53/17 to 8/37*I + 137/37

            sage: A = diagonal_matrix([1, -1, 1])
            sage: A = HyperbolicPlane().HM().get_isometry(A)
            sage: A.preserves_orientation()
            False
            sage: p = HyperbolicPlane().HM().get_point((0, 1, sqrt(2)))
            sage: A * p
            Point in HM (0, -1, sqrt(2))
        """
        if isinstance(other, HyperbolicIsometry):
            other = other.to_model(self.codomain())
            return self.__class__(self.codomain(), self._matrix*other._matrix)
        from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPoint
        if isinstance(other, HyperbolicPoint):
            return self(other)
        if isinstance(other, HyperbolicGeodesic):
            return self.codomain().get_geodesic(self(other.start()), self(other.end()))

        raise NotImplementedError("multiplication is not defined between a "
                                  "hyperbolic isometry and {0}".format(other))

    def _call_(self, p):
        r"""
        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = UHP.get_isometry(Matrix(2,[5,2,1,2]))
            sage: p = UHP.get_point(2 + I)
            sage: A(p)
            Point in UHP 8/17*I + 53/17

            sage: A = diagonal_matrix([1, -1, 1])
            sage: A = HyperbolicPlane().HM().get_isometry(A)
            sage: A.preserves_orientation()
            False
            sage: p = HyperbolicPlane().HM().get_point((0, 1, sqrt(2)))
            sage: A(p)
            Point in HM (0, -1, sqrt(2))

            sage: I2 = UHP.get_isometry(identity_matrix(2))
            sage: p = UHP.random_point()
            sage: bool(UHP.dist(I2(p), p) < 10**-9)
            True
        """
        return self.codomain().get_point(self._matrix * vector(p._coordinates))

    #######################
    # Setters and Getters #
    #######################

    def matrix(self):
        r"""
        Return the matrix of the isometry.

        .. NOTE::

            We do not allow the ``matrix`` constructor to work as these may
            be elements of a projective group (ex. `PSL(n, \RR)`), so these
            isometries aren't true matrices.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_isometry(-identity_matrix(2)).matrix()
            [-1  0]
            [ 0 -1]
        """
        return self._matrix

    def __invert__(self):
        r"""
        Return the inverse of the isometry ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = UHP.get_isometry(matrix(2,[4,1,3,2]))
            sage: B = A.inverse()   # indirect doctest
            sage: A*B == UHP.get_isometry(identity_matrix(2))
            True
        """
        return self.__class__(self.domain(), self.matrix().__invert__())

    def is_identity(self):
        """
        Return ``True`` if ``self`` is the identity isometry.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_isometry(matrix(2,[4,1,3,2])).is_identity()
            False
            sage: UHP.get_isometry(identity_matrix(2)).is_identity()
            True
        """
        return self._matrix.is_one()

    def model(self):
        r"""
        Return the model to which ``self`` belongs.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_isometry(identity_matrix(2)).model()
            Hyperbolic plane in the Upper Half Plane Model

            sage: HyperbolicPlane().PD().get_isometry(identity_matrix(2)).model()
            Hyperbolic plane in the Poincare Disk Model

            sage: HyperbolicPlane().KM().get_isometry(identity_matrix(3)).model()
            Hyperbolic plane in the Klein Disk Model

            sage: HyperbolicPlane().HM().get_isometry(identity_matrix(3)).model()
            Hyperbolic plane in the Hyperboloid Model
        """
        return self.domain()

    def to_model(self, other):
        r"""
        Convert the current object to image in another model.

        INPUT:

        - ``other`` -- (a string representing) the image model

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: UHP = H.UHP()
            sage: PD = H.PD()
            sage: KM = H.KM()
            sage: HM = H.HM()

            sage: A = UHP.get_isometry(identity_matrix(2))
            sage: A.to_model(HM)
            Isometry in HM
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: A.to_model('HM')
            Isometry in HM
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: A = PD.get_isometry(matrix([[I, 0], [0, -I]]))
            sage: A.to_model(UHP)
            Isometry in UHP
            [ 0  1]
            [-1  0]
            sage: A.to_model(HM)
            Isometry in HM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]
            sage: A.to_model(KM)
            Isometry in KM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: A = HM.get_isometry(diagonal_matrix([-1, -1, 1]))
            sage: A.to_model('UHP')
            Isometry in UHP
            [ 0 -1]
            [ 1  0]
            sage: A.to_model('PD')
            Isometry in PD
            [-I  0]
            [ 0  I]
            sage: A.to_model('KM')
            Isometry in KM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]
        """
        if isinstance(other, str):
            other = getattr(self.domain().realization_of(), other)()
        if other is self.domain():
            return self
        phi = other.coerce_map_from(self.domain())
        return phi.convert_isometry(self)

    ###################
    # Boolean Methods #
    ###################

    def preserves_orientation(self):
        r"""
        Return ``True`` if ``self`` is orientation-preserving and ``False``
        otherwise.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = UHP.get_isometry(identity_matrix(2))
            sage: A.preserves_orientation()
            True
            sage: B = UHP.get_isometry(matrix(2,[0,1,1,0]))
            sage: B.preserves_orientation()
            False
        """
        return self._cached_isometry.preserves_orientation()

    def classification(self):
        r"""
        Classify the hyperbolic isometry as elliptic, parabolic,
        hyperbolic or a reflection.

        A hyperbolic isometry fixes two points on the boundary of
        hyperbolic space, a parabolic isometry fixes one point on the
        boundary of hyperbolic space, and an elliptic isometry fixes no
        points.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: H = UHP.get_isometry(matrix(2,[2,0,0,1/2]))
            sage: H.classification()
            'hyperbolic'

            sage: P = UHP.get_isometry(matrix(2,[1,1,0,1]))
            sage: P.classification()
            'parabolic'

            sage: E = UHP.get_isometry(matrix(2,[-1,0,0,1]))
            sage: E.classification()
            'reflection'
        """
        return self._cached_isometry.classification()

    def translation_length(self):
        r"""
        For hyperbolic elements, return the translation length;
        otherwise, raise a :exc:`ValueError`.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: H = UHP.get_isometry(matrix(2,[2,0,0,1/2]))
            sage: H.translation_length()
            2*arccosh(5/4)

        ::

            sage: f_1 = UHP.get_point(-1)
            sage: f_2 = UHP.get_point(1)
            sage: H = UHP.isometry_from_fixed_points(f_1, f_2)
            sage: p = UHP.get_point(exp(i*7*pi/8))
            sage: bool((p.dist(H*p) - H.translation_length()) < 10**-9)
            True
        """
        return self._cached_isometry.translation_length()

    def axis(self):
        r"""
        For a hyperbolic isometry, return the axis of the
        transformation; otherwise raise a :exc:`ValueError`.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: H = UHP.get_isometry(matrix(2,[2,0,0,1/2]))
            sage: H.axis()
            Geodesic in UHP from 0 to +Infinity

        It is an error to call this function on an isometry that is
        not hyperbolic::

            sage: P = UHP.get_isometry(matrix(2,[1,4,0,1]))
            sage: P.axis()
            Traceback (most recent call last):
            ...
            ValueError: the isometry is not hyperbolic: axis is undefined
        """
        if self.classification() not in ['hyperbolic',
                                         'orientation-reversing hyperbolic']:
            raise ValueError("the isometry is not hyperbolic: axis is undefined")
        return self.fixed_point_set()

    def fixed_point_set(self):
        r"""
        Return a list containing the fixed point set of
        orientation-preserving isometries.

        OUTPUT: list of hyperbolic points or a hyperbolic geodesic

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: H = KM.get_isometry(matrix([[5/3,0,4/3], [0,1,0], [4/3,0,5/3]]))
            sage: g = H.fixed_point_set(); g
            Geodesic in KM from (1, 0) to (-1, 0)
            sage: H(g.start()) == g.start()
            True
            sage: H(g.end()) == g.end()
            True
            sage: A = KM.get_isometry(matrix([[1,0,0], [0,-1,0], [0,0,1]]))
            sage: A.preserves_orientation()
            False
            sage: A.fixed_point_set()
            Geodesic in KM from (1, 0) to (-1, 0)

        ::

            sage: B = KM.get_isometry(identity_matrix(3))
            sage: B.fixed_point_set()
            Traceback (most recent call last):
            ...
            ValueError: the identity transformation fixes the entire hyperbolic plane
        """
        M = self.domain()
        pts = self._cached_isometry.fixed_point_set()
        if isinstance(pts, HyperbolicGeodesic):
            return pts.to_model(M)
        return [M(k) for k in pts]

    def fixed_geodesic(self):
        r"""
        If ``self`` is a reflection in a geodesic, return that geodesic.

        EXAMPLES::

            sage: A = HyperbolicPlane().PD().get_isometry(matrix([[0, 1], [1, 0]]))
            sage: A.fixed_geodesic()
            Geodesic in PD from -1 to 1
        """
        fps = self._cached_isometry.fixed_point_set()
        if not isinstance(fps, HyperbolicGeodesic):
            raise ValueError("isometries of type {0}".format(self.classification())
                             + " do not fix geodesics")
        return fps.to_model(self.domain())

    def repelling_fixed_point(self):
        r"""
        For a hyperbolic isometry, return the attracting fixed point;
        otherwise raise a :exc:`ValueError`.

        OUTPUT: a hyperbolic point

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = UHP.get_isometry(Matrix(2,[4,0,0,1/4]))
            sage: A.repelling_fixed_point()
            Boundary point in UHP 0
        """
        fp = self._cached_isometry.repelling_fixed_point()
        return self.domain().get_point(fp)

    def attracting_fixed_point(self):
        r"""
        For a hyperbolic isometry, return the attracting fixed point;
        otherwise raise a :exc:`ValueError`.

        OUTPUT: a hyperbolic point

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = UHP.get_isometry(Matrix(2,[4,0,0,1/4]))
            sage: A.attracting_fixed_point()
            Boundary point in UHP +Infinity
        """
        fp = self._cached_isometry.attracting_fixed_point()
        return self.domain().get_point(fp)

class HyperbolicIsometryUHP(HyperbolicIsometry):
    r"""
    Create a hyperbolic isometry in the UHP model.

    INPUT:

    - a matrix in `GL(2, \RR)`

    EXAMPLES::

        sage: HyperbolicPlane().UHP().get_isometry(identity_matrix(2))
        Isometry in UHP
        [1 0]
        [0 1]
    """
    def _call_(self, p): #UHP
        r"""
        Return image of ``p`` under the action of ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: I2 = UHP.get_isometry(identity_matrix(2))
            sage: p = UHP.random_point()
            sage: bool(UHP.dist(I2(p), p) < 10**-9)
            True
        """
        coords = p.coordinates()
        # We apply complex conjugation to the point for negative determinants
        # We check the coordinate is not equal to infinity. If we use !=, then
        #   it cannot determine it is not infinity, so it also returns False.
        if not (coords == infinity) and bool(self._matrix.det() < 0):
            coords = coords.conjugate()
        return self.codomain().get_point(moebius_transform(self._matrix, coords))

    def preserves_orientation(self): #UHP
        r"""
        Return ``True`` if ``self`` is orientation-preserving and ``False``
        otherwise.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = identity_matrix(2)
            sage: UHP.get_isometry(A).preserves_orientation()
            True
            sage: B = matrix(2,[0,1,1,0])
            sage: UHP.get_isometry(B).preserves_orientation()
            False
        """
        return bool(self._matrix.det() > 0)

    def classification(self): #UHP
        r"""
        Classify the hyperbolic isometry as elliptic, parabolic, or
        hyperbolic.

        A hyperbolic isometry fixes two points on the boundary of
        hyperbolic space, a parabolic isometry fixes one point on the
        boundary of hyperbolic space, and an elliptic isometry fixes
        no points.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_isometry(identity_matrix(2)).classification()
            'identity'

            sage: UHP.get_isometry(4*identity_matrix(2)).classification()
            'identity'

            sage: UHP.get_isometry(matrix(2,[2,0,0,1/2])).classification()
            'hyperbolic'

            sage: UHP.get_isometry(matrix(2, [0, 3, -1/3, 6])).classification()
            'hyperbolic'

            sage: UHP.get_isometry(matrix(2,[1,1,0,1])).classification()
            'parabolic'

            sage: UHP.get_isometry(matrix(2,[-1,0,0,1])).classification()
            'reflection'
        """
        A = self._matrix.n()
        A = A / (abs(A.det()).sqrt())
        tau = abs(A.trace())
        a = A.list()
        if A.det() > 0:
            tf = bool((a[0] - 1)**2 + a[1]**2 + a[2]**2 + (a[3] - 1)**2 < EPSILON)
            if tf:
                return 'identity'
            if tau - 2 < -EPSILON:
                return 'elliptic'
            if tau - 2 > -EPSILON and tau - 2 < EPSILON:
                return 'parabolic'
            if tau - 2 > EPSILON:
                return 'hyperbolic'
            raise ValueError("something went wrong with classification:" +
                             " trace is {}".format(A.trace()))
        # Otherwise The isometry reverses orientation
        if tau < EPSILON:
            return 'reflection'
        return 'orientation-reversing hyperbolic'

    def translation_length(self): #UHP
        r"""
        For hyperbolic elements, return the translation length;
        otherwise, raise a :exc:`ValueError`.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_isometry(matrix(2,[2,0,0,1/2])).translation_length()
            2*arccosh(5/4)

        ::

            sage: H = UHP.isometry_from_fixed_points(-1,1)
            sage: p = UHP.get_point(exp(i*7*pi/8))
            sage: Hp = H(p)
            sage: bool((UHP.dist(p, Hp) - H.translation_length()) < 10**-9)
            True
        """
        d = sqrt(self._matrix.det()**2)
        tau = sqrt((self._matrix / sqrt(d)).trace()**2)
        if self.classification() in ['hyperbolic', 'orientation-reversing hyperbolic']:
            return 2 * arccosh(tau / 2)
        raise TypeError("translation length is only defined for hyperbolic transformations")

    def fixed_point_set(self):  # UHP
        r"""
        Return a list or geodesic containing the fixed point set of
        orientation-preserving isometries.

        OUTPUT: list of hyperbolic points or a hyperbolic geodesic

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: H = UHP.get_isometry(matrix(2, [-2/3,-1/3,-1/3,-2/3]))
            sage: g = H.fixed_point_set(); g
            Geodesic in UHP from -1 to 1
            sage: H(g.start()) == g.start()
            True
            sage: H(g.end()) == g.end()
            True
            sage: A = UHP.get_isometry(matrix(2,[0,1,1,0]))
            sage: A.preserves_orientation()
            False
            sage: A.fixed_point_set()
            Geodesic in UHP from 1 to -1

        ::

            sage: B = UHP.get_isometry(identity_matrix(2))
            sage: B.fixed_point_set()
            Traceback (most recent call last):
            ...
            ValueError: the identity transformation fixes the entire hyperbolic plane
        """
        d = sqrt(self._matrix.det() ** 2)
        M = self._matrix / sqrt(d)
        tau = M.trace() ** 2
        M_cls = self.classification()
        if M_cls == 'identity':
            raise ValueError("the identity transformation fixes the entire "
                             "hyperbolic plane")

        pt = self.domain().get_point
        if M_cls == 'parabolic':
            if abs(M[1, 0]) < EPSILON:
                return [pt(infinity)]
            else:
                # boundary point
                return [pt((M[0,0] - M[1,1]) / (2*M[1,0]))]
        elif M_cls == 'elliptic':
            d = sqrt(tau - 4)
            return [pt((M[0,0] - M[1,1] + sign(M[1,0])*d) / (2*M[1,0]))]
        elif M_cls == 'hyperbolic':
            if M[1,0] != 0: #if the isometry doesn't fix infinity
                d = sqrt(tau - 4)
                p_1 = (M[0,0] - M[1,1]+d) / (2*M[1,0])
                p_2 = (M[0,0] - M[1,1]-d) / (2*M[1,0])
                return self.domain().get_geodesic(pt(p_1), pt(p_2))
            #else, it fixes infinity.
            p_1 = M[0,1] / (M[1,1] - M[0,0])
            p_2 = infinity
            return self.domain().get_geodesic(pt(p_1), pt(p_2))

        try:
            p, q = (M.eigenvectors_right()[k][1][0] for k in range(2))
        except IndexError:
            M = M.change_ring(RDF)
            p, q = (M.eigenvectors_right()[k][1][0] for k in range(2))

        pts = []
        if p[1] == 0:
            pts.append(infinity)
        else:
            p = p[0] / p[1]
            if imag(p) >= 0:
                pts.append(p)
        if q[1] == 0:
            pts.append(infinity)
        else:
            q = q[0] / q[1]
            if imag(q) >= 0:
                pts.append(q)
        pts = [pt(k) for k in pts]
        if len(pts) == 2:
            return self.domain().get_geodesic(*pts)
        return pts

    def repelling_fixed_point(self):  # UHP
        r"""
        Return the repelling fixed point.

        Otherwise, this raises a :exc:`ValueError`.

        OUTPUT: a hyperbolic point

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = matrix(2,[4,0,0,1/4])
            sage: UHP.get_isometry(A).repelling_fixed_point()
            Boundary point in UHP 0
        """
        if self.classification() not in ['hyperbolic',
                                         'orientation-reversing hyperbolic']:
            raise ValueError("repelling fixed point is defined only" +
                             "for hyperbolic isometries")
        v = self._matrix.eigenmatrix_right()[1].column(1)
        if v[1] == 0:
            return self.domain().get_point(infinity)
        return self.domain().get_point(v[0] / v[1])

    def attracting_fixed_point(self):  # UHP
        r"""
        Return the attracting fixed point.

        Otherwise, this raises a :exc:`ValueError`.

        OUTPUT: a hyperbolic point

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = matrix(2,[4,0,0,1/4])
            sage: UHP.get_isometry(A).attracting_fixed_point()
            Boundary point in UHP +Infinity
        """
        if self.classification() not in \
                ['hyperbolic', 'orientation-reversing hyperbolic']:
            raise ValueError("Attracting fixed point is defined only" +
                             "for hyperbolic isometries.")
        v = self._matrix.eigenmatrix_right()[1].column(0)
        if v[1] == 0:
            return self.domain().get_point(infinity)
        return self.domain().get_point(v[0] / v[1])

class HyperbolicIsometryPD(HyperbolicIsometry):
    r"""
    Create a hyperbolic isometry in the PD model.

    INPUT:

    - a matrix in `PU(1,1)`

    EXAMPLES::

        sage: HyperbolicPlane().PD().get_isometry(identity_matrix(2))
        Isometry in PD
        [1 0]
        [0 1]
    """
    def _call_(self, p): #PD
        r"""
        Return image of ``p`` under the action of ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: I2 = PD.get_isometry(identity_matrix(2))
            sage: q = PD.random_point()
            sage: bool(PD.dist(I2(q), q) < 10**-9)
            True
        """
        coords = p.coordinates()
        # We apply complex conjugation to the point for negative determinants
        if bool(self._matrix.det() < 0):
            coords = coords.conjugate()
        _image = moebius_transform(self._matrix, coords)
        return self.codomain().get_point(_image)

    def __mul__(self, other): #PD
        r"""
        Return image of ``p`` under the action of ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: X = PD.get_isometry(matrix([[3/4, -I/4], [-I/4, -3/4]]))
            sage: X*X
            Isometry in PD
            [   5/8  3/8*I]
            [-3/8*I    5/8]
        """
        if isinstance(other, HyperbolicIsometry):
            M = self._cached_isometry * other._cached_isometry
            return M.to_model('PD')
        return super().__mul__(other)

    def __pow__(self, n): #PD
        r"""
        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: X = PD.get_isometry(matrix([[3/4, -I/4], [-I/4, -3/4]]))
            sage: X^2
            Isometry in PD
            [   5/8  3/8*I]
            [-3/8*I    5/8]
        """
        return (self._cached_isometry**n).to_model('PD')

    def preserves_orientation(self): #PD
        """
        Return ``True`` if ``self`` preserves orientation and ``False``
        otherwise.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: PD.get_isometry(matrix([[-I, 0], [0, I]])).preserves_orientation()
            True
            sage: PD.get_isometry(matrix([[0, I], [I, 0]])).preserves_orientation()
            False
        """
        return bool(self._matrix.det() > 0) and HyperbolicIsometryPD._orientation_preserving(self._matrix)

    @staticmethod
    def _orientation_preserving(A): #PD
        r"""
        For a matrix ``A`` of a PD isometry, determine if it preserves
        orientation.

        This test is more involved than just checking the sign of
        the determinant.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometryPD
            sage: orient = HyperbolicIsometryPD._orientation_preserving
            sage: orient(matrix([[-I, 0], [0, I]]))
            True
            sage: orient(matrix([[0, I], [I, 0]]))
            False
        """
        return bool(A[1][0] == A[0][1].conjugate() and A[1][1] == A[0][0].conjugate()
                    and abs(A[0][0]) - abs(A[0][1]) != 0)

class HyperbolicIsometryKM(HyperbolicIsometry):
    r"""
    Create a hyperbolic isometry in the KM model.

    INPUT:

    - a matrix in `SO(2,1)`

    EXAMPLES::

        sage: HyperbolicPlane().KM().get_isometry(identity_matrix(3))
        Isometry in KM
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    def _call_(self, p): #KM
        r"""
        Return image of ``p`` under the action of ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: I3 = KM.get_isometry(identity_matrix(3))
            sage: v = KM.random_point()
            sage: bool(KM.dist(I3(v), v) < 10**-9)
            True
        """
        v = self._matrix * vector(list(p.coordinates()) + [1])
        if v[2] == 0:
            return self.codomain().get_point(infinity)
        return self.codomain().get_point(v[0:2] / v[2])

#####################################################################
## Helper functions


def moebius_transform(A, z):
    r"""
    Given a matrix ``A`` in `GL(2, \CC)` and a point ``z`` in the complex
    plane return the Möbius transformation action of ``A`` on ``z``.

    INPUT:

    - ``A`` -- a `2 \times 2` invertible matrix over the complex numbers
    - ``z`` -- a complex number or infinity

    OUTPUT: a complex number or infinity

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_model import moebius_transform
        sage: moebius_transform(matrix(2,[1,2,3,4]),2 + I)
        -2/109*I + 43/109
        sage: y = var('y')
        sage: moebius_transform(matrix(2,[1,0,0,1]),x + I*y)
        x + I*y

    The matrix must be square and `2 \times 2`::

        sage: moebius_transform(matrix([[3,1,2],[1,2,5]]),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring

        sage: moebius_transform(identity_matrix(3),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring

    The matrix can be symbolic or can be a matrix over the real
    or complex numbers, but must be provably invertible::

        sage: a,b,c,d = var('a,b,c,d')
        sage: moebius_transform(matrix(2,[a,b,c,d]),I)
        (I*a + b)/(I*c + d)
        sage: moebius_transform(matrix(2,[1,b,c,b*c+1]),I)
        (b + I)/(b*c + I*c + 1)
        sage: moebius_transform(matrix(2,[0,0,0,0]),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring
    """
    if A.ncols() == 2 == A.nrows() and A.det() != 0:
        a, b, c, d = A.list()
        if z == infinity:
            if c == 0:
                return infinity
            return a / c
        if c * z + d == 0:
            return infinity
        return (a * z + b) / (c * z + d)
    raise TypeError("A must be an invertible 2x2 matrix over the"
                    " complex numbers or a symbolic ring")
