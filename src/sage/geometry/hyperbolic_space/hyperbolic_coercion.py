"""
Coercion Maps Between Hyperbolic Plane Models

This module implements the coercion maps between different hyperbolic
plane models.

AUTHORS:

- Travis Scrimshaw (2014): initial version
"""

#***********************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.categories.morphism import Morphism
from sage.symbolic.constants import I
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.functions.other import real, imag
from sage.misc.functional import sqrt
from sage.misc.lazy_import import lazy_import
lazy_import('sage.misc.call', 'attrcall')


class HyperbolicModelCoercion(Morphism):
    """
    Abstract base class for morphisms between the hyperbolic models.
    """
    def _repr_type(self):
        """
        Return the type of morphism.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi._repr_type()
            'Coercion Isometry'
        """
        return "Coercion Isometry"

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: HM = HyperbolicPlane().HM()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi(PD.get_point(0.5+0.5*I))
            Point in UHP 2.00000000000000 + 1.00000000000000*I
            sage: psi = HM.coerce_map_from(UHP)
            sage: psi(UHP.get_point(I))
            Point in HM (0, 0, 1)

        It is an error to try to convert a boundary point to a model
        that doesn't support boundary points::

            sage: psi(UHP.get_point(infinity))
            Traceback (most recent call last):
            ...
            NotImplementedError: boundary points are not implemented for the Hyperboloid Model

        It is an error to try to convert a boundary point to a model
        that doesn't support boundary points::

            sage: psi(UHP(infinity))
            Traceback (most recent call last):
            ...
            NotImplementedError: boundary points are not implemented for the Hyperboloid Model
        """
        C = self.codomain()
        if not C.is_bounded() and self.domain().is_bounded() and x.is_boundary():
            msg = "boundary points are not implemented for the {}"
            raise NotImplementedError(msg.format(C.name()))

        y = self.image_coordinates(x.coordinates())
        if self.domain().is_bounded():
            bdry = x.is_boundary()
        else:
            bdry = C.boundary_point_in_model(y)

        return C.element_class(C, y, bdry, check=False, **x.graphics_options())

    def convert_geodesic(self, x):
        """
        Convert the geodesic ``x`` of the domain into a geodesic of
        the codomain.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi.convert_geodesic(PD.get_geodesic(0.5+0.5*I, -I))
            Geodesic in UHP from 2.00000000000000 + 1.00000000000000*I to 0
        """
        return self.codomain().get_geodesic(self(x.start()), self(x.end()),
                                            **x.graphics_options())

    def convert_isometry(self, x):
        """
        Convert the hyperbolic isometry ``x`` of the domain into an
        isometry of the codomain.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: HM = HyperbolicPlane().HM()
            sage: phi = HM.coerce_map_from(UHP)
            sage: I2 = UHP.get_isometry(identity_matrix(2))
            sage: phi.convert_isometry(I2)
            Isometry in HM
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        C = self.codomain()
        return C._Isometry(C, self.image_isometry_matrix(x._matrix), check=False)

    def __invert__(self):
        """
        Return the inverse coercion of ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = UHP.coerce_map_from(PD)
            sage: ~phi
            Coercion Isometry morphism:
              From: Hyperbolic plane in the Upper Half Plane Model
              To:   Hyperbolic plane in the Poincare Disk Model
        """
        return self.domain().coerce_map_from(self.codomain())

############
# From UHP #
############


class CoercionUHPtoPD(HyperbolicModelCoercion):
    """
    Coercion from the UHP to PD model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = PD.coerce_map_from(UHP)
            sage: phi.image_coordinates(I)
            0
        """
        if x == infinity:
            return I
        return (x - I) / (Integer(1) - I*x)

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = PD.coerce_map_from(UHP)
            sage: phi.image_isometry_matrix(identity_matrix(2))
            [1 0]
            [0 1]
        """
        if x.det() < 0:
            # x = I * x
            return matrix([[1,-I],[-I,1]]) * x * matrix([[1,I],[I,1]]).conjugate()/Integer(2)
        return matrix([[1,-I],[-I,1]]) * x * matrix([[1,I],[I,1]])/Integer(2)


class CoercionUHPtoKM(HyperbolicModelCoercion):
    """
    Coercion from the UHP to KM model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: KM = HyperbolicPlane().KM()
            sage: phi = KM.coerce_map_from(UHP)
            sage: phi.image_coordinates(3 + I)
            (6/11, 9/11)
        """
        if x == infinity:
            return (0, 1)
        return ((2*real(x))/(real(x)**2 + imag(x)**2 + 1),
                (real(x)**2 + imag(x)**2 - 1)/(real(x)**2 + imag(x)**2 + 1))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: KM = HyperbolicPlane().KM()
            sage: phi = KM.coerce_map_from(UHP)
            sage: phi.image_isometry_matrix(identity_matrix(2))
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return SL2R_to_SO21(x)


class CoercionUHPtoHM(HyperbolicModelCoercion):
    """
    Coercion from the UHP to HM model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: HM = HyperbolicPlane().HM()
            sage: phi = HM.coerce_map_from(UHP)
            sage: phi.image_coordinates(3 + I)
            (3, 9/2, 11/2)
        """
        return vector((real(x)/imag(x),
                      (real(x)**2 + imag(x)**2 - 1)/(2*imag(x)),
                      (real(x)**2 + imag(x)**2 + 1)/(2*imag(x))))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: HM = HyperbolicPlane().HM()
            sage: phi = HM.coerce_map_from(UHP)
            sage: phi.image_isometry_matrix(identity_matrix(2))
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return SL2R_to_SO21(x)

###########
# From PD #
###########


class CoercionPDtoUHP(HyperbolicModelCoercion):
    """
    Coercion from the PD to UHP model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: UHP = HyperbolicPlane().UHP()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi.image_coordinates(0.5+0.5*I)
            2.00000000000000 + 1.00000000000000*I
            sage: phi.image_coordinates(0)
            I
            sage: phi.image_coordinates(I)
            +Infinity
            sage: phi.image_coordinates(-I)
            0
        """
        if x == I:
            return infinity
        return (x + I)/(Integer(1) + I*x)

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES:

        We check that orientation-reversing isometries behave as they
        should::

            sage: PD = HyperbolicPlane().PD()
            sage: UHP = HyperbolicPlane().UHP()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi.image_isometry_matrix(matrix([[0,I],[I,0]]))
            [-1  0]
            [ 0 -1]
        """
        from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometryPD
        if not HyperbolicIsometryPD._orientation_preserving(x):
            return matrix([[1,I],[I,1]]) * x * matrix([[1,-I],[-I,1]]).conjugate() / Integer(2)
        return matrix([[1,I],[I,1]]) * x * matrix([[1,-I],[-I,1]]) / Integer(2)


class CoercionPDtoKM(HyperbolicModelCoercion):
    """
    Coercion from the PD to KM model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: KM = HyperbolicPlane().KM()
            sage: phi = KM.coerce_map_from(PD)
            sage: phi.image_coordinates(0.5+0.5*I)
            (0.666666666666667, 0.666666666666667)
        """
        return (2*real(x)/(Integer(1) + real(x)**2 + imag(x)**2),
                2*imag(x)/(Integer(1) + real(x)**2 + imag(x)**2))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: KM = HyperbolicPlane().KM()
            sage: phi = KM.coerce_map_from(PD)
            sage: phi.image_isometry_matrix(matrix([[0,I],[I,0]]))
            [-1  0  0]
            [ 0  1  0]
            [ 0  0 -1]
        """
        return SL2R_to_SO21(matrix(2, [1, I, I, 1]) * x *
                            matrix(2, [1, -I, -I, 1]) / Integer(2))


class CoercionPDtoHM(HyperbolicModelCoercion):
    """
    Coercion from the PD to HM model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: HM = HyperbolicPlane().HM()
            sage: phi = HM.coerce_map_from(PD)
            sage: phi.image_coordinates(0.5+0.5*I)
            (2.00000000000000, 2.00000000000000, 3.00000000000000)
        """
        return vector((2*real(x)/(1 - real(x)**2 - imag(x)**2),
                       2*imag(x)/(1 - real(x)**2 - imag(x)**2),
                       (real(x)**2 + imag(x)**2 + 1) /
                       (1 - real(x)**2 - imag(x)**2)))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: HM = HyperbolicPlane().HM()
            sage: phi = HM.coerce_map_from(PD)
            sage: phi.image_isometry_matrix(matrix([[0,I],[I,0]]))
            [-1  0  0]
            [ 0  1  0]
            [ 0  0 -1]
        """
        return SL2R_to_SO21(matrix(2, [1, I, I, 1]) * x *
                            matrix(2, [1, -I, -I, 1]) / Integer(2))

###########
# From KM #
###########


class CoercionKMtoUHP(HyperbolicModelCoercion):
    """
    Coercion from the KM to UHP model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: UHP = HyperbolicPlane().UHP()
            sage: phi = UHP.coerce_map_from(KM)
            sage: phi.image_coordinates((0, 0))
            I
            sage: phi.image_coordinates((0, 1))
            +Infinity
        """
        if tuple(x) == (0, 1):
            return infinity
        return (-x[0]/(x[1] - 1)
                + I*(-(sqrt(-x[0]**2 - x[1]**2 + 1) - x[0]**2 - x[1]**2 + 1)
                     / ((x[1] - 1)*sqrt(-x[0]**2 - x[1]**2 + 1) + x[1] - 1)))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: UHP = HyperbolicPlane().UHP()
            sage: phi = UHP.coerce_map_from(KM)
            sage: m = matrix([[5/3,0,4/3], [0,1,0], [4/3,0,5/3]])
            sage: phi.image_isometry_matrix(m)
            [2*sqrt(1/3)   sqrt(1/3)]
            [  sqrt(1/3) 2*sqrt(1/3)]
        """
        return SO21_to_SL2R(x)


class CoercionKMtoPD(HyperbolicModelCoercion):
    """
    Coercion from the KM to PD model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = PD.coerce_map_from(KM)
            sage: phi.image_coordinates((0, 0))
            0
        """
        return (x[0]/(1 + (1 - x[0]**2 - x[1]**2).sqrt())
                + I*x[1]/(1 + (1 - x[0]**2 - x[1]**2).sqrt()))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = PD.coerce_map_from(KM)
            sage: m = matrix([[5/3,0,4/3], [0,1,0], [4/3,0,5/3]])
            sage: phi.image_isometry_matrix(m)
            [2*sqrt(1/3)   sqrt(1/3)]
            [  sqrt(1/3) 2*sqrt(1/3)]
        """
        return (matrix(2,[1,-I,-I,1]) * SO21_to_SL2R(x) *
                matrix(2,[1,I,I,1])/Integer(2))


class CoercionKMtoHM(HyperbolicModelCoercion):
    """
    Coercion from the KM to HM model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: HM = HyperbolicPlane().HM()
            sage: phi = HM.coerce_map_from(KM)
            sage: phi.image_coordinates((0, 0))
            (0, 0, 1)
        """
        return (vector((2*x[0], 2*x[1], 1 + x[0]**2 + x[1]**2))
                / (1 - x[0]**2 - x[1]**2))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: HM = HyperbolicPlane().HM()
            sage: phi = HM.coerce_map_from(KM)
            sage: m = matrix([[5/3,0,4/3], [0,1,0], [4/3,0,5/3]])
            sage: phi.image_isometry_matrix(m)
            [5/3   0 4/3]
            [  0   1   0]
            [4/3   0 5/3]
        """
        return x

###########
# From HM #
###########


class CoercionHMtoUHP(HyperbolicModelCoercion):
    """
    Coercion from the HM to UHP model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: UHP = HyperbolicPlane().UHP()
            sage: phi = UHP.coerce_map_from(HM)
            sage: phi.image_coordinates( vector((0,0,1)) )
            I
        """
        return -((x[0]*x[2] + x[0]) + I*(x[2] + 1)) / ((x[1] - 1)*x[2]
                                        - x[0]**2 - x[1]**2 + x[1] - 1)

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: UHP = HyperbolicPlane().UHP()
            sage: phi = UHP.coerce_map_from(HM)
            sage: phi.image_isometry_matrix(identity_matrix(3))
            [1 0]
            [0 1]
        """
        return SO21_to_SL2R(x)


class CoercionHMtoPD(HyperbolicModelCoercion):
    """
    Coercion from the HM to PD model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = PD.coerce_map_from(HM)
            sage: phi.image_coordinates( vector((0,0,1)) )
            0
        """
        return x[0]/(1 + x[2]) + I*(x[1]/(1 + x[2]))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = PD.coerce_map_from(HM)
            sage: phi.image_isometry_matrix(identity_matrix(3))
            [1 0]
            [0 1]
        """
        return (matrix(2,[1,-I,-I,1]) * SO21_to_SL2R(x) *
                matrix(2,[1,I,I,1])/Integer(2))


class CoercionHMtoKM(HyperbolicModelCoercion):
    """
    Coercion from the HM to KM model.
    """
    def image_coordinates(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: KM = HyperbolicPlane().KM()
            sage: phi = KM.coerce_map_from(HM)
            sage: phi.image_coordinates( vector((0,0,1)) )
            (0, 0)
        """
        return (x[0]/(1 + x[2]), x[1]/(1 + x[2]))

    def image_isometry_matrix(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: KM = HyperbolicPlane().KM()
            sage: phi = KM.coerce_map_from(HM)
            sage: phi.image_isometry_matrix(identity_matrix(3))
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return x

#####################################################################
## Helper functions


def SL2R_to_SO21(A):
    r"""
    Given a matrix in `SL(2, \RR)` return its irreducible representation in
    `O(2,1)`.

    Note that this is not the only homomorphism, but it is the only one
    that works in the context of the implemented 2D hyperbolic geometry
    models.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_coercion import SL2R_to_SO21
        sage: A = SL2R_to_SO21(identity_matrix(2))
        sage: J = matrix([[1,0,0],[0,1,0],[0,0,-1]]) #Lorentzian Gram matrix
        sage: norm(A.transpose()*J*A - J) < 10**-4                                      # needs scipy
        True
    """
    a, b, c, d = (A/A.det().sqrt()).list()

    # Kill ~0 imaginary parts
    components = [
        a*d + b*c, a*c - b*d, a*c + b*d, a*b - c*d,
        Integer(1)/Integer(2)*a**2 - Integer(1)/Integer(2)*b**2 -
                Integer(1)/Integer(2)*c**2 + Integer(1)/Integer(2)*d**2,
        Integer(1)/Integer(2)*a**2 + Integer(1)/Integer(2)*b**2 -
                Integer(1)/Integer(2)*c**2 - Integer(1)/Integer(2)*d**2,
        a*b + c*d, Integer(1)/Integer(2)*a**2 -
                Integer(1)/Integer(2)*b**2 + Integer(1)/Integer(2)*c**2 -
        Integer(1)/Integer(2)*d**2, Integer(1)/Integer(2)*a**2 +
                Integer(1)/Integer(2)*b**2 + Integer(1)/Integer(2)*c**2 +
        Integer(1)/Integer(2)*d**2
    ]
    B = matrix(3, [real(comp) for comp in components])

    #B = B.apply_map(attrcall('real'))
    if A.det() > 0:
        return B
    else:
        # Orientation-reversing isometries swap the nappes of
        #  the lightcone.  This fixes that issue.
        return -B


def SO21_to_SL2R(M):
    r"""
    A homomorphism from `SO(2, 1)` to `SL(2, \RR)`.

    Note that this is not the only homomorphism, but it is the only one
    that works in the context of the implemented 2D hyperbolic geometry
    models.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_coercion import SO21_to_SL2R
        sage: (SO21_to_SL2R(identity_matrix(3)) - identity_matrix(2)).norm() < 10**-4   # needs scipy
        True
    """
    ####################################################################
    # SL(2,R) is the double cover of SO (2,1)^+, so we need to choose  #
    # a lift.  I have formulas for the absolute values of each entry   #
    # a,b ,c,d of the lift matrix(2,[a,b,c,d]), but we need to choose  #
    # one entry to be positive.  I choose d for no particular reason,  #
    # unless d = 0, then we choose c > 0.  The basic strategy for this #
    # function is to find the linear map induced by the SO(2,1)        #
    # element on the Lie algebra sl(2, R).  This corresponds to the    #
    # Adjoint action by a matrix A or -A in SL(2,R).  To find which    #
    # matrix let X,Y,Z be a basis for sl(2,R) and look at the images   #
    # of X,Y,Z as well as the second and third standard basis vectors  #
    # for 2x2 matrices (these are traceless, so are in the Lie         #
    # algebra).  These corresponds to AXA^-1 etc and give formulas     #
    # for the entries of A.                                            #
    ####################################################################
    (m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9) = M.list()
    d = sqrt(Integer(1)/Integer(2)*m_5 - Integer(1)/Integer(2)*m_6 -
             Integer(1)/Integer(2)*m_8 + Integer(1)/Integer(2)*m_9)
    if M.det() > 0:  # EPSILON?
        det_sign = 1
    elif M.det() < 0:  # EPSILON?
        det_sign = -1
    if d > 0:  # EPSILON?
        c = (-Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/d
        b = (-Integer(1)/Integer(2)*m_2 + Integer(1)/Integer(2)*m_3)/d
        ad = det_sign*1 + b*c  # ad - bc = pm 1
        a = ad/d
    else:  # d is 0, so we make c > 0
        c = sqrt(-Integer(1)/Integer(2)*m_5 - Integer(1)/Integer(2)*m_6 +
                 Integer(1)/Integer(2)*m_8 + Integer(1)/Integer(2)*m_9)
        d = (-Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/c
            # d = 0, so ad - bc = -bc = pm 1.
        b = - (det_sign*1)/c
        a = (Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/b
    A = matrix(2, [a, b, c, d])
    return A
