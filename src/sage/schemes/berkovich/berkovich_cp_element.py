# sage.doctest: needs sage.rings.padics
r"""
Elements of Berkovich space.

:class:`Berkovich_Element` is an abstract parent class for elements of any Berkovich space.

:class:`Berkovich_Element_Cp_Affine` and :class:`Berkovich_Element_Cp_Projective`
implement elements of Berkovich space over `\CC_p` and `P^1(\CC_p)`. Elements are
determined by specific data and fall into one of the four following types:

- Type I points are represented by a center.

- Type II points are represented by a center and a rational power of `p`.

- Type III points are represented by a center and a nonnegative real radius.

- Type IV points are represented by a finite list of centers and a finite list of
  nonnegative radii.

For an exposition of Berkovich space over `\CC_p`, see Chapter 6 of [Ben2019]_. For a more
involved exposition, see Chapter 1 and 2 of [BR2010]_.

AUTHORS:

- Alexander Galarraga (2020-06-22): initial implementation
"""

# *****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import sage.rings.abc

from sage.categories.function_fields import FunctionFields
from sage.misc.lazy_import import lazy_import
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealNumber, RR
from sage.schemes.projective.projective_point import SchemeMorphism_point_projective_field
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.structure.element import Element, Expression

lazy_import('sage.rings.padics.padic_generic_element', 'pAdicGenericElement')
lazy_import('sage.rings.padics.padic_base_generic', 'pAdicBaseGeneric')


class Berkovich_Element(Element):
    """
    The parent class for any element of a Berkovich space.
    """
    pass


class Berkovich_Element_Cp(Berkovich_Element):
    r"""
    The abstract parent class for any element of Berkovich space over `\CC_p`.

    This class should never be instantiated, instead use :class:`Berkovich_Element_Cp_Affine`
    or :class:`Berkovich_Element_Cp_Projective`.

    EXAMPLES::

        sage: B = Berkovich_Cp_Affine(3)
        sage: B(2)
        Type I point centered at 2 + O(3^20)

    ::

        sage: B(0, 1)
        Type II point centered at 0 of radius 3^0
    """

    def __init__(self, parent, center, radius=None, power=None, prec=20, space_type=None, error_check=True):
        """
        Initialization function.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(5)
            sage: B(4)
            Type I point centered at 4 + O(5^20)
        """
        from sage.rings.polynomial.polynomial_element import Polynomial
        from sage.rings.fraction_field_element import FractionFieldElement_1poly_field
        self._type = None

        # if radius is a list or a tuple, this is a type 4 point
        if isinstance(radius, (list, tuple)):
            if error_check:
                if not (isinstance(center, list) or isinstance(center, tuple)):
                    raise TypeError("center was passed a list but radius was not a list")
                if len(radius) != len(center):
                    raise ValueError("the same number of centers and radii "
                                     "must be specified to create "
                                     "a type IV point")
            self._center_lst = list(center)
            self._radius_lst = list(radius)
            self._prec = len(self._radius_lst)
            self._center_func = None
            self._radius_func = None
            self._type = 4
            self._radius = None
            self._center = None
            if not error_check:
                return

        elif isinstance(center, Element) and isinstance(radius, Element):
            from sage.rings.polynomial.multi_polynomial import MPolynomial
            if isinstance(center, MPolynomial):
                try:
                    center = center.univariate_polynomial()
                except AttributeError:
                    raise TypeError('center was %s, a multivariable polynomial' % center)

            # check if the radius and the center are functions
            center_func_check = center.parent() in FunctionFields() or \
                isinstance(center, (Polynomial, FractionFieldElement_1poly_field, Expression))
            radius_func_check = radius.parent() in FunctionFields() or \
                isinstance(radius, (Polynomial, FractionFieldElement_1poly_field, Expression))

            if center_func_check:
                # check that both center and radii are supported univariate function
                center_expr_check = False
                radius_expr_check = False
                if error_check:
                    if isinstance(center, Expression):
                        if len(center.variables()) != 1:
                            raise ValueError("an expression with %s " % (len(center.variables())) +
                                             "variables cannot define the centers approximating a type IV point")
                        else:
                            # we do this since .subs is currently buggy for polynomials but not expressions
                            center_expr_check = True
                    if not radius_func_check:
                        raise TypeError("center was passed a function but radius was not a function")
                    if isinstance(radius, Expression):
                        if len(radius.variables()) != 1:
                            raise ValueError("an expression with %s " % (len(radius.variables())) +
                                             "variables cannot define the radii approximating a type IV point")
                        else:
                            radius_expr_check = True
                else:
                    if isinstance(center, Expression):
                        center_expr_check = True
                    if isinstance(radius, Expression):
                        radius_expr_check = True
                self._type = 4
                self._prec = prec
                center_lst = []
                radius_lst = []
                self._center_func = center
                self._radius_func = radius
                if center_expr_check:
                    x = self._center_func.variables()[0]
                if radius_expr_check:
                    y = self._radius_func.variables()[0]
                for i in range(1, self._prec + 1):
                    if center_expr_check:
                        # we use .subs for expressions to avoid deprecation
                        center_lst.append(self._center_func.subs({x: i}))
                    else:
                        # .subs for polynomials is currently buggy
                        center_lst.append(self._center_func(i))
                    if radius_expr_check:
                        radius_lst.append(self._radius_func.subs({y: i}))
                    else:
                        radius_lst.append(self._radius_func(i))
                self._center_lst = center_lst
                self._radius_lst = radius_lst
                self._radius = None
                self._center = None
                if not error_check:
                    return

        if self._type == 4 and error_check:
            if space_type == "projective":
                for i in range(len(self._center_lst)):
                    center = self._center_lst[i]
                    radius = self._radius_lst[i]
                    # make sure the center is a point of projective space and not the point at infinity
                    if not isinstance(center, SchemeMorphism_point_projective_field):
                        try:
                            center = (self._base_space)(center)
                        except (TypeError, ValueError):
                            raise TypeError('could not convert %s to %s' % (center, self._base_space))
                    if self._base_type == 'padic field':
                        if not isinstance(center.scheme().base_ring(), sage.rings.abc.pAdicField):
                            if not isinstance(center.scheme().base_ring(), pAdicBaseGeneric):
                                try:
                                    center = (self._base_space)(center)
                                except (TypeError, ValueError):
                                    raise ValueError("could not convert %s to %s" % (center, self._base_space))
                            else:
                                # center is padic, not but an element of a scheme over a padic field.
                                # we convert to scheme over a padic field
                                center = ProjectiveSpace(center.scheme().base_ring().fraction_field(), 1)(center)
                        if center.scheme().base_ring().prime() != self._p:
                            raise ValueError("center must be an element of " +
                                             "%s not %s" % self._base_space, center.scheme())
                    else:
                        if center not in self._base_space:
                            try:
                                center = (self._base_space)(center)
                            except (TypeError, ValueError):
                                raise ValueError('could not convert %s to %s' % (center, self._base_space))
                    if center.scheme().ambient_space() != center.scheme():
                        raise ValueError("the center of a point of Berkovich space over " +
                                         "P^1(Cp(%s)) must be a point of Cp not %s" % (self._p, center.scheme()))
                    if center == (center.scheme())((1, 0)):
                        raise ValueError("the center of a disk approximating a type IV point of Berkovich " +
                                         "space cannot be centered at %s" % ((center.scheme())((1, 0))))
                    # since we are over a field, we can normalize coordinates. all code assumes normalized coordinates
                    center.normalize_coordinates()
                    # make sure the radius coerces into the reals
                    if not isinstance(radius, RealNumber):
                        if isinstance(radius, Expression):
                            radius = RR(radius)
                        elif RR.has_coerce_map_from(radius.parent()):
                            radius = RR(radius)
                        else:
                            raise TypeError("the radius of a disk approximating a type IV point" +
                                            "must coerce into the real numbers, %s does not coerce" % (radius))
                    if i != 0:
                        # check containment for the sequence of disks
                        previous_center = self._center_lst[i - 1]
                        previous_radius = self._radius_lst[i - 1]
                        dist = self._custom_abs(center[0] - previous_center[0])
                        if previous_radius < radius or dist > previous_radius:
                            raise ValueError("sequence of disks does not define a type IV point as " +
                                             "containment is not proper")
                    self._center_lst[i] = center
                    self._radius_lst[i] = radius
                return
            elif space_type == "affine":
                for i in range(len(self._center_lst)):
                    center = self._center_lst[i]
                    radius = self._radius_lst[i]
                    if self._base_type == 'padic field':
                        # make sure the center is in Cp
                        if not isinstance(center, pAdicGenericElement):
                            try:
                                center = (self._base_space)(center)
                            except (TypeError, ValueError):
                                raise TypeError("could not convert %s to %s" % (center, self._base_space))
                        elif not isinstance(center.parent(), sage.rings.abc.pAdicField):
                            # center is padic, not but an element of a padic field. we convert to padic field
                            center = (center.parent().fraction_field())(center)
                        if (center.parent()).prime() != self._p:
                            raise ValueError("center in %s, should be in %s") % (center.parent(), self._base_space)
                    else:
                        # make sure the center is in the appropriate number field
                        if center.parent() == self._base_space:
                            try:
                                center = (self._base_space)(center)
                            except (TypeError, ValueError):
                                raise ValueError('could not convert %s to %s' % (center, self._base_space))
                    # make sure the radius coerces into the reals
                    if not isinstance(radius, RealNumber):
                        if isinstance(radius, Expression):
                            radius = RR(radius)
                        elif RR.has_coerce_map_from(radius.parent()):
                            radius = RR(radius)
                            self._radius_lst[i] = radius
                        else:
                            raise ValueError("the radius of a disk approximating a type IV point must " +
                                             "coerce into the real numbers, %s does not coerce" % (radius))
                    if i != 0:
                        # check containment for the sequence of disks
                        previous_center = self._center_lst[i - 1]
                        previous_radius = self._radius_lst[i - 1]
                        dist = self._custom_abs(center - previous_center)
                        if previous_radius < radius or dist > previous_radius:
                            raise ValueError("sequence of disks does not define a type IV point as " +
                                             "containment is not proper")
                    self._center_lst[i] = center
                    self._radius_lst[i] = radius
                return
            else:
                raise ValueError("bad value %s passed to space_type. Do not initialize  " % (space_type) +
                                 "Berkovich_Element_Cp directly")

        # the point must now be type 1, 2, or 3, so we check that the center is of the appropriate type
        if error_check:
            if space_type == "projective":
                if not isinstance(center, SchemeMorphism_point_projective_field):
                    try:
                        center = (self._base_space)(center)
                    except (ValueError, TypeError):
                        raise TypeError("could not convert %s to %s" % (center, self._base_space))
                if self._base_type == 'padic field':
                    if not isinstance(center.scheme().base_ring(), sage.rings.abc.pAdicField):
                        if not isinstance(center.scheme().base_ring(), pAdicBaseGeneric):
                            try:
                                center = (self._base_space)(center)
                            except (TypeError, ValueError):
                                raise ValueError("could not convert %s to %s" % (center, self._base_space))
                        else:
                            # center is padic, not but an element of a scheme over a padic field.
                            # we convert to scheme over a padic field
                            field_scheme = ProjectiveSpace(center.scheme().base_ring().fraction_field(), 1)
                            try:
                                center = field_scheme(center)
                            except (TypeError, ValueError):
                                raise ValueError('could not convert %s to %s' % center, field_scheme)
                    if center.scheme().base_ring().prime() != self._p:
                        raise ValueError("center must be an element of " +
                                         "%s not %s" % self._base_space, center.scheme())
                else:
                    if center not in self._base_space:
                        try:
                            center = (self._base_space)(center)
                        except (TypeError, ValueError):
                            raise ValueError('could not convert %s to %s' % (center, self._base_space))
                if center.scheme().ambient_space() is not center.scheme():
                    raise ValueError("the center of a point of projective Berkovich space cannot be " +
                                     "a point of %s" % (center.scheme()))
                # since we are over a field, we normalize coordinates
                center.normalize_coordinates()
            elif space_type == 'affine':
                if self._base_type == 'padic field':
                    # make sure the center is in Cp
                    if not isinstance(center, pAdicGenericElement):
                        try:
                            center = (self._base_space)(center)
                        except (TypeError, ValueError):
                            raise TypeError("could not convert %s to %s" % (center, self._base_space))
                    elif not isinstance(center.parent(), sage.rings.abc.pAdicField):
                        # center is padic, not but an element of a padic field. we convert to padic field
                        center = (center.parent().fraction_field())(center)
                    if (center.parent()).prime() != self._p:
                        raise ValueError("center in %s, should be in %s") % (center.parent(), self._base_space)
                else:
                    # make sure the center is in the appropriate number field
                    if not (center.parent() == self._base_space):
                        try:
                            center = (self._base_space)(center)
                        except (TypeError, ValueError):
                            raise ValueError('could not convert %s to %s' % (center, self._base_space))
            else:
                raise ValueError("bad value %s passed to space_type. Do not initialize  " % (space_type) +
                                 "Berkovich_Element_Cp directly")

        self._center = center

        # since this point is not type IV, these are None
        self._center_func = None
        self._center_lst = None
        self._radius_lst = None
        self._radius_func = None

        if (radius is None and power is None) or radius == 0:
            self._type = 1
            self._radius = 0
            self._power = None
            return
        # In order to simplify our representation, type II and III points cannot be centered at infinity
        if space_type == "projective":
            # TODO use involution map to allow for infinity to be passed in as center
            if center[1] == 0:
                raise ValueError('type II and III points can not be centered at infinity')
        if power is not None:
            if error_check:
                try:
                    power = QQ(power)
                except TypeError:
                    raise TypeError("power must convert to rationals")
                if radius is not None:
                    if radius != RR(self._p**power):
                        raise ValueError("conflicting inputs for power and radius")
            self._power = power
            self._radius = RR(self._p**power)
            self._type = 2
            return
        if radius is not None:
            if isinstance(radius, Expression):
                try:
                    power = QQ(radius.log(self._p).expand_log())
                except TypeError:
                    pass
                try:
                    radius = RR(radius)
                    self._radius = radius
                except TypeError:
                    if len(radius.variables()) == 1:
                        raise ValueError('radius univariate function but center is constant. ' +
                                         'this does not define a type IV point')
                    raise TypeError("symbolic radius must be a real number")
            if (not isinstance(radius, RealNumber)) and power is None:
                if RR.has_coerce_map_from(radius.parent()):
                    self._radius = RR(radius)
                else:
                    raise TypeError("radius must coerce into real numbers")
            else:
                self._radius = radius
            if power is not None:
                self._power = power
                self._type = 2
                return
            power = RR(radius.log(self._p))
            if power.is_integer():
                self._power = QQ(power)
                self._type = 2
            else:
                self._type = 3
                self._power = power
            return

        raise ValueError('unknown error constructing point of Berkovich space over Cp')

    def _custom_abs(self, x):
        """
        Return the absolute value of ``x`` with respect to the norm on ``Cp``.

        Used to simplify code, as ``x`` may be a point of a number field
        or a `p`-adic field.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(QQ, 3)
            sage: Q1 = B(9)
            sage: Q1._custom_abs(Q1.center())
            1/9

        ::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: Q1 = B(9)
            sage: Q1._custom_abs(Q1.center())
            1/9
        """
        if self._base_type == 'padic field':
            return x.abs()
        if x.valuation(self._ideal) == Infinity:
            return 0
        if self._ideal in QQ:
            return self.prime()**(-x.valuation(self._ideal))
        return self.prime()**(-x.valuation(self._ideal) / self._ideal.absolute_ramification_index())

    def center_function(self):
        """
        Return the function defining the centers of disks in the approximation.

        Not defined unless this point is a type IV point created by using
        a univariate function to compute centers.

        OUTPUT: a univariate function

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(5)
            sage: L.<t> = PolynomialRing(Qp(5))
            sage: T = FractionField(L)
            sage: f = T(1/t)
            sage: R.<x> = RR[]
            sage: Y = FractionField(R)
            sage: g = (40*pi)/x                                                         # needs sage.symbolic
            sage: Q1 = B(f, g)                                                          # needs sage.symbolic
            sage: Q1.center_function()                                                  # needs sage.symbolic
            (1 + O(5^20))/((1 + O(5^20))*t)
        """
        if self.type_of_point() != 4:
            raise ValueError('center_function not defined for points which are not type IV')
        if self._center_func is None:
            raise ValueError('this type IV point does not have a center function')
        return self._center_func

    def radius_function(self):
        """
        Return the function defining the radii of disks in the approximation.

        Not defined unless this point is a type IV point created by using
        a univariate function to compute radii.

        OUTPUT: a univariate function

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(5)
            sage: L.<t> = PolynomialRing(Qp(5))
            sage: T = FractionField(L)
            sage: f = T(1/t)
            sage: R.<x> = RR[]
            sage: Y = FractionField(R)
            sage: g = (40*pi)/x                                                         # needs sage.symbolic
            sage: Q1 = B(f, g)                                                          # needs sage.symbolic
            sage: Q1.radius_function()                                                  # needs sage.symbolic
            40.0000000000000*pi/x
        """
        if self.type_of_point() != 4:
            raise ValueError('center_function not defined for points which are not type IV')
        if self._radius_func is None:
            raise ValueError('this type IV point does not have a radius function')
        return self._radius_func

    def precision(self):
        """
        Return the precision of a type IV point.

        This integer is the number of disks used in the approximation of the type IV point.
        Not defined for type I, II, or III points.

        OUTPUT: integer

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: d = B([2, 2, 2], [1.761, 1.123, 1.112])
            sage: d.precision()
            3

        TESTS::

            sage: d.precision == d.prec
            True
        """
        if self._type in [1, 2, 3]:
            raise AttributeError("type I, II, and III points do not have a precision")
        return self._prec

    prec = precision

    def ideal(self):
        r"""
        The ideal which defines an embedding of the ``base_ring`` into `\CC_p`.

        If this Berkovich space is backed by a `p`-adic field, then an embedding is
        already specified, and this returns ``None``.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(QQ, 3)
            sage: B(0).ideal()
            3

            ::

            sage: B = Berkovich_Cp_Projective(3)
            sage: B(0).ideal()
        """
        return self.parent().ideal()

    def power(self):
        r"""
        The power of ``p`` such that `p^\text{power} = \text{radius}`.

        For type II points, always in `\QQ`. For type III points,
        a real number. Not defined for type I or IV points.

        OUTPUT:

        - A rational for type II points.
        - A real number for type III points.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(1, 9)
            sage: Q1.power()
            2

        ::

            sage: Q2 = B(1, 4)
            sage: Q2.power()
            1.26185950714291
        """
        if self._type in [1, 4]:
            raise AttributeError("type I and IV points do not have a power")
        return self._power

    def radius(self):
        r"""
        Radius of the corresponding disk (or sequence of disks) in `\CC_p`.

        OUTPUT:

        - A nonnegative real number for type I, II, or III points.
        - A list of nonnegative real numbers for type IV points.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(1, 2/5)
            sage: Q1.radius()
            0.400000000000000

        ::

            sage: d = B([2, 2, 2], [1.761, 1.123, 1.112])
            sage: d.radius()
            [1.76100000000000, 1.12300000000000, 1.11200000000000]
        """
        if self._type == 4:
            return self._radius_lst
        return self._radius

    def diameter(self, basepoint=Infinity):
        r"""
        Generalized diameter function on Berkovich space.

        If the basepoint is infinity, the diameter is equal to
        the limit of the radii of the corresponding disks in `\CC_p`.

        If the basepoint is not infinity, the diameter
        is the Hsia kernel of this point with itself at
        basepoint ``basepoint``.

        INPUT:

        - ``basepoint`` -- (default: ``Infinity``) a point of the
          same Berkovich space as this point

        OUTPUT: a real number

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(3)
            sage: Q1.diameter()
            0

        ::

            sage: Q2 = B(1/2, 9)
            sage: Q2.diameter()
            9.00000000000000

        The diameter of a type IV point is the limit of the radii::

            sage: R.<x> = PolynomialRing(Qp(3))
            sage: f = R(2)
            sage: S.<y> = PolynomialRing(RR)
            sage: S = FractionField(S)
            sage: g = (y+1)/y
            sage: B(f,g).diameter()
            1.0

        ::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(1/81, 1)
            sage: Q2 = B(1/3)
            sage: Q1.diameter(Q2)
            0.00137174211248285

        ::

            sage: Q2.diameter(Q2)
            +infinity
        """
        if basepoint == Infinity:
            if self._type == 4:
                if self._radius_func is None:
                    return self._radius_lst[-1]
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                R = PolynomialRing(QQ, names='x')
                x = R.gens()[0]
                if isinstance(self._radius_func, Expression):
                    radius_func_variable = self._radius_func.variables()[0]
                    radius_expr = self._radius_func.subs({radius_func_variable: x})
                else:
                    radius_expr = self._radius_func(x)
                    from sage.symbolic.ring import SymbolicRing as SR
                    radius_expr = SR(RR)(radius_expr)
                return radius_expr.limit(x='oo')
            return self._radius
        if not isinstance(basepoint, Berkovich_Element_Cp):
            raise TypeError('basepoint must be a point of Berkovich space, not %s' % basepoint)
        if basepoint.parent() != self.parent():
            raise ValueError('basepoint must be a point of the same Berkovich space')
        return self.Hsia_kernel(self, basepoint)

    def path_distance_metric(self, other):
        r"""
        Return the path distance metric distance between this point and ``other``.

        Also referred to as the hyperbolic metric, or the big metric.

        On the set of type II, III and IV points, the path distance metric
        is a metric. Following Baker and Rumely, we extend
        the path distance metric to type I points `x`, `y` by `\rho(x,x) = 0` and `\rho(x,y) =
        \infty`. See [BR2010]_.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point

        OUTPUT: a finite or infinite real number

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(1/4, 4)
            sage: Q2 = B(1/4, 6)
            sage: Q1.path_distance_metric(Q2)
            0.369070246428542

        ::

            sage: Q3 = B(1)
            sage: Q3.path_distance_metric(Q1)
            +infinity

        ::

            sage: Q3.path_distance_metric(Q3)
            0
        """
        if not isinstance(other, type(self)):
            raise TypeError('other must be a point of Berkovich space. other was %s' % other)
        if self.parent() != other.parent():
            raise ValueError("other must be a point of the same Berkovich space")
        if self.type_of_point() == 1 or other.type_of_point() == 1:
            if self == other:
                return 0
            else:
                return RR(Infinity)
        return 2 * self.join(other).diameter().log(self.prime()) \
            - self.diameter().log(self.prime()) \
            - other.diameter().log(other.prime())

    big_metric = path_distance_metric

    hyperbolic_metric = path_distance_metric

    def Hsia_kernel(self, other, basepoint):
        """
        The Hsia kernel of this point and ``other``,
        with basepoint ``basepoint``.

        The Hsia kernel with arbitrary basepoint
        is a generalization of the Hsia kernel at infinity.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point
        - ``basepoint`` -- a point of the same Berkovich space as this point

        OUTPUT: a finite or infinite real number

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(2, 9)
            sage: Q2 = B(1/27, 1/27)
            sage: Q3 = B(1, 1/3)
            sage: Q1.Hsia_kernel(Q2, Q3)
            0.111111111111111

        ::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(2, 9)
            sage: Q2 = B(1/2)
            sage: Q3 = B(1/2)
            sage: Q1.Hsia_kernel(Q2, Q3)
            +infinity
        """
        if not isinstance(other, type(self)):
            raise TypeError('other must be a point of Berkovich space. other was %s' % other)
        if self.parent() != other.parent():
            raise ValueError("other must be a point of the same Berkovich space")
        if not isinstance(basepoint, type(self)):
            raise TypeError('basepoint must be a point of Berkovich space. basepoint was %s' % basepoint)
        if basepoint.parent() != self.parent():
            raise ValueError("basepoint must be a point of the same Berkovich space")
        if basepoint.type_of_point() == 1:
            if self == basepoint or other == basepoint:
                return RR(Infinity)
        return self.spherical_kernel(other) / \
            (self.spherical_kernel(basepoint) * other.spherical_kernel(basepoint))

    def small_metric(self, other):
        r"""
        Return the small metric distance between this point and ``other``.

        The small metric is an extension of twice
        the spherical distance on `P^1(\CC_p)`.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point

        OUTPUT: a real number

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(1/4, 4)
            sage: Q2 = B(1/4, 6)
            sage: Q1.small_metric(Q2)
            0.0833333333333333

        ::

            sage: B = Berkovich_Cp_Projective(QQ, 5)
            sage: Q1 = B(0, 1)
            sage: Q2 = B(99)
            sage: Q1.small_metric(Q2)
            1.00000000000000

        ::

            sage: Q3 = B(1/4, 4)
            sage: Q3.small_metric(Q2)
            1.75000000000000

        ::

            sage: Q2.small_metric(Q3)
            1.75000000000000
        """
        if not isinstance(other, Berkovich_Element_Cp):
            raise TypeError('other must be a point of affine Berkovich space. other was %s' % other)
        if self.parent() != other.parent():
            raise ValueError('other must be a point of the same Berkovich space')
        gauss = self.parent()(RR(0), RR(1))
        g_greater_than_s = gauss.gt(self)
        g_greater_than_o = gauss.gt(other)
        if g_greater_than_s and g_greater_than_o:
            return 2 * self.join(other, gauss).diameter() - self.diameter() - other.diameter()
        if not g_greater_than_s:
            new_self = self.involution_map()
        else:
            new_self = self
        if not g_greater_than_o:
            new_other = other.involution_map()
        else:
            new_other = other
        return 2 * new_self.join(new_other, gauss).diameter() \
            - new_self.diameter() - new_other.diameter()

    def potential_kernel(self, other, basepoint):
        """
        The potential kernel of this point with ``other``,
        with basepoint ``basepoint``.

        The potential kernel is the hyperbolic distance
        between ``basepoint`` and the join of this point
        with ``other`` relative to ``basepoint``.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point
        - ``basepoint`` -- a point of the same Berkovich space as this point

        OUTPUT: a finite or infinite real number

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(27, 1)
            sage: Q2 = B(1/3, 2)
            sage: Q3 = B(1/9, 1/2)
            sage: Q3.potential_kernel(Q1, Q2)
            0.369070246428543

        ::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(27, 1)
            sage: Q2 = B(1/3, 2)
            sage: Q3 = B(1/9, 1/2)
            sage: Q3.potential_kernel(Q1, Q2)
            0.369070246428543
        """
        if not isinstance(other, type(self)):
            raise TypeError('other must be a point of a Berkovich space, not %s' % other)
        if other.parent() != self.parent():
            raise ValueError('other must be a point of the same Berkovich space')
        if not isinstance(basepoint, type(self)):
            raise TypeError('basepoint must be a point of Berkovich line, not %s' % basepoint)
        if basepoint.parent() != self.parent():
            raise ValueError('basepoint must be a point of the same Berkovich space')
        return basepoint.path_distance_metric(self.join(other, basepoint))

    def spherical_kernel(self, other):
        r"""
        The spherical kernel of this point with ``other``.

        The spherical kernel is one possible extension of the spherical
        distance on `P^1(\CC_p)` to the projective Berkovich line.
        See [BR2010]_ for details.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point

        OUTPUT: a real number

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(2, 2)
            sage: Q2 = B(1/9, 1)
            sage: Q1.spherical_kernel(Q2)
            0.500000000000000

        ::

            sage: Q3 = B(2)
            sage: Q3.spherical_kernel(Q3)
            0
        """
        if not isinstance(other, type(self)):
            raise TypeError('other must be a point of Berkovich space, not %s' % other)
        if other.parent() != self.parent():
            raise ValueError('other must be a point of the same Berkovich space')
        gauss_point = self.parent()(ZZ(0), ZZ(1))
        w = self.join(other, gauss_point)
        dist = gauss_point.path_distance_metric(w)
        if dist == Infinity:
            return 0
        return self.prime()**(-dist)

    def Hsia_kernel_infinity(self, other):
        r"""
        Return the Hsia kernel at infinity of this point with ``other``.

        The Hsia kernel at infinity is the natural extension of the
        absolute value on `\CC_p` to Berkovich space.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point

        OUTPUT: a real number

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: Q1 = B(1/4, 4)
            sage: Q2 = B(1/4, 6)
            sage: Q1.Hsia_kernel_infinity(Q2)
            6.00000000000000

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^3 + 20)
            sage: ideal = A.ideal(-1/2*a^2 + a - 3)
            sage: B = Berkovich_Cp_Projective(A, ideal)
            sage: Q1 = B(4)
            sage: Q2 = B(0, 1.5)
            sage: Q1.Hsia_kernel_infinity(Q2)
            1.50000000000000
        """
        return self.join(other).diameter()

    def center(self):
        r"""
        Return the center of the corresponding disk (or sequence of disks)
        in `\CC_p`.

        OUTPUT: an element of the ``base`` of the parent Berkovich space

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: B(3, 1).center()
            3 + O(3^21)

        ::

            sage: C = Berkovich_Cp_Projective(3)
            sage: C(3, 1).center()
            (3 + O(3^21) : 1 + O(3^20))

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^3 + 20)
            sage: ideal = A.ideal(-1/2*a^2 + a - 3)
            sage: B = Berkovich_Cp_Projective(A, ideal)
            sage: B(a^2 + 4).center()
            (a^2 + 4 : 1)
        """
        if self._type == 4:
            return self._center_lst
        return self._center

    def type_of_point(self):
        r"""
        Return the type of this point of Berkovich space over `\CC_p`.

        OUTPUT: integer between 1 and 4 inclusive

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: B(1).type_of_point()
            1

        ::

            sage: B(0, 1).type_of_point()
            2
        """
        return ZZ(self._type)

    def prime(self):
        """
        The residue characteristic of the parent.

        OUTPUT: a prime integer

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: B(1).prime()
            3
        """
        return ZZ(self._p)

    def __ne__(self, other):
        """
        Non-equality operator.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(3, 3**(1/2))
            sage: Q2 = B(3, RR(3**(1/2)))
            sage: Q1 != Q2
            False
        """
        return not (self == other)

    def _repr_(self):
        """
        String representation of this point.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: B(2, 1)
            Type II point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 3^0
        """
        if self._type == 1:
            return "Type I point centered at " + format(self._center)
        elif self._type == 2:
            return "Type II point centered at " \
                + format(self._center) \
                + " of radius %s^%s" % (self._p, self._power)
        elif self._type == 3:
            return "Type III point centered at " \
                + format(self._center) + " of radius " \
                + format(self._radius)
        else:
            if self._center_func is not None and self._radius_func is not None:
                return "Type IV point of precision %s " % self._prec + \
                    "with centers given by %s and radii given by %s"\
                    % (self._center_func, self._radius_func)
            else:
                return "Type IV point of precision %s, approximated " % self._prec + \
                    "by disks centered at %s ... with radii %s ..." \
                    % (self._center_lst[:min(self._prec, 2)], self._radius_lst[:min(self._prec, 2)])

    def _latex_(self):
        r"""
        LaTeX representation of this point.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: latex(B(2, 1))
            \text{type 2 Point of } \text{Projective Berkovich line over }
            \Bold{C}_{3} \text{equivalent to the disk centered at
            (2 + O(3^20) : 1 + O(3^20)) of radius 1.00000000000000 in } \Bold{C}_3
        """
        from sage.misc.latex import latex
        if self._type == 1:
            text = r"the point %s of } \Bold{C}_%s" % (self._center, self._p)
        elif self._type in [2, 3]:
            text = r"the disk centered at %s of radius %s in } \Bold{C}_%s" \
                % (self._center, self._radius, self._p)
        else:
            text = "the sequence of disks with centers %s } " % self._center_lst[:2] + \
                r"\ldots \text{ and radii %s } \ldots" % self._radius_lst[:2]
        return r"\text{type %s Point of }" % (self._type) \
            + latex(self.parent()) + r"\text{equivalent to " + text


class Berkovich_Element_Cp_Affine(Berkovich_Element_Cp):
    r"""
    Element class of the Berkovich affine line over `\CC_p`.

    Elements are categorized into four types, represented by specific data:

    - Type I points are represented by a center in the ``base`` of the parent Berkovich space,
      which is `\QQ_p`, a finite extension of `\QQ_p`, or a number field.

    - Type II points are represented by a center in the ``base`` of the parent Berkovich space,
      and a rational power of `p`.

    - Type III points are represented by a center in the ``base`` of the parent Berkovich space,
      and a radius, a real number in `[0,\infty)`.

    - Type IV points are represented by a finite list of centers in the ``base`` of the parent
      Berkovich space and a finite list of radii in `[0,\infty)`. Type IV points can be created
      from univariate functions, allowing for arbitrary precision.

    INPUT:

    - ``center`` -- for type I, II, and III points, the center of the
      corresponding disk in `\CC_p`. If the parent Berkovich space was created using a number field
      `K`, then ``center`` must be an element of `K`. Otherwise, ``center`` must be an element of a
      `p`-adic field. For type IV points, can be a list of centers used to approximate the point or a
      univariate function that computes the centers (computation starts at 1).

    - ``radius`` -- (optional) For type I, II, and III points, the radius of the
      corresponding disk in ``Cp``. Must coerce into the real numbers. For type IV points,
      can be a list of radii used to approximate the point or a univariate function that
      computes the radii (computation starts at 1).

    - ``power`` -- (optional) Rational number. Used for constructing type II points; specifies
      the power of ``p`` such that `p^\text{power}` = radius

    - ``prec`` -- (default: 20) the number of disks to be used to approximate a type IV point

    - ``error_check`` -- boolean (default: ``True``); if error checking should be run on input. If
      input is correctly formatted, can be set to ``False`` for better performance.
      WARNING: with error check set to ``False``, any error in the input will lead to
      incorrect results.

    EXAMPLES:

    Type I points can be created by specifying the corresponding point of ``Cp``::

        sage: B = Berkovich_Cp_Affine(Qp(3))
        sage: B(4)
        Type I point centered at 1 + 3 + O(3^20)

    The center of a point can be an element of a finite extension of ``Qp``::

        sage: A.<t> = Qq(27)
        sage: B(1 + t)
        Type I point centered at (t + 1) + O(3^20)

    Type II and III points can be created by specifying a center and a radius::

        sage: B(2, 3**(1/2))                                                            # needs sage.symbolic
        Type II point centered at 2 + O(3^20) of radius 3^1/2

    ::

        sage: B(2, 1.6)
        Type III point centered at 2 + O(3^20) of radius 1.60000000000000

    Some type II points may be mistaken for type III points::

        sage: B(3, 3**0.5)                      # not tested
        Type III point centered at 3 + O(3^21) of radius 1.73205080756888

    To avoid these errors, specify the power instead of the radius::

        sage: B(3, power=RR(1/100000))
        Type II point centered at 3 + O(3^21) of radius 3^1/100000

    Type IV points can be constructed in a number of ways, the first being
    from a list of centers and radii used to approximate the point::

        sage: B([Qp(3)(2), Qp(3)(2), Qp(3)(2)], [1.761, 1.123, 1.112])
        Type IV point of precision 3, approximated by disks centered at
        [2 + O(3^20), 2 + O(3^20)] ... with radii [1.76100000000000, 1.12300000000000] ...

    Type IV points can be constructed from univariate functions, with arbitrary precision::

        sage: A.<t> = Qq(27)
        sage: R.<x> = PolynomialRing(A)
        sage: f = (1 + t)^2*x
        sage: S.<y> = PolynomialRing(RR)
        sage: S = FractionField(S)
        sage: g = (y + 1)/y
        sage: d = B(f, g, prec=100); d
        Type IV point of precision 100 with centers given by
        ((t^2 + 2*t + 1) + O(3^20))*x and radii given by (y + 1.00000000000000)/y

    For increased performance, ``error_check`` can be set to ``False``. WARNING: with error check set
    to ``False``, any error in the input will lead to incorrect results::

        sage: B(f, g, prec=100, error_check=False)
        Type IV point of precision 100 with centers given by
        ((t^2 + 2*t + 1) + O(3^20))*x and radii given by (y + 1.00000000000000)/y

    When creating a Berkovich space backed by a number field, points can be created similarly::

        sage: # needs sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: A.<a> = NumberField(x^3 + 20)
        sage: ideal = A.prime_above(3)
        sage: B = Berkovich_Cp_Projective(A, ideal)
        sage: Q1 = B(a); Q1
        Type I point centered at (a : 1)

    ::

        sage: B(a + 1, 3)                                                               # needs sage.rings.number_field
        Type II point centered at (a + 1 : 1) of radius 3^1

    TESTS::

        sage: A = Berkovich_Cp_Affine(3)
        sage: Q1 = A(3, 1); Q1
        Type II point centered at 3 + O(3^21) of radius 3^0
        sage: Q2 = A(2.5, 1); Q2
        Type II point centered at 1 + 2*3 + 3^2 + 3^3 + 3^4 + 3^5 + 3^6 + 3^7 +
        3^8 + 3^9 + 3^10 + 3^11 + 3^12 + 3^13 + 3^14 + 3^15 + 3^16 + 3^17 +
        3^18 + 3^19 + O(3^20) of radius 3^0
        sage: Q5 = A(3, 0); Q5
        Type I point centered at 3 + O(3^21)
        sage: A(Zp(3)(2), 2).center().parent() == A(Qp(3)(2), 2).center().parent()
        True
        sage: Q1 == Q2
        True
        sage: Q1 == Q5
        False
        sage: Q3 = A(Qp(3)(3), power=0, error_check=False); Q3
        Type II point centered at 3 + O(3^21) of radius 3^0
        sage: Q4 = A(3, 3**0); Q4
        Type II point centered at 3 + O(3^21) of radius 3^0
        sage: Q5 = A(3, power=1/2); Q5
        Type II point centered at 3 + O(3^21) of radius 3^1/2
        sage: Q6 = A(3, RR(3**(1/2))); Q6                                               # needs sage.symbolic
        Type III point centered at 3 + O(3^21) of radius 1.73205080756888
        sage: Q5 == Q6                                                                  # needs sage.symbolic
        True

        sage: k = Qp(5)
        sage: R.<x> = k[]
        sage: l.<w> = k.extension(x^2 - 5)
        sage: B = Berkovich_Cp_Affine(5)
        sage: B(w, power=1)
        Type II point centered at w + O(w^41) of radius 5^1

        sage: TestSuite(Q5).run()
    """

    def __init__(self, parent, center, radius=None, power=None, prec=20, error_check=True):
        """
        Initialization function.

        EXAMPLES::

            sage: A = Berkovich_Cp_Affine(17)
            sage: A(5, 1)
            Type II point centered at 5 + O(17^20) of radius 17^0
        """
        # we call Berkovich_Element_Cp constructor which is shared with projective Berkovich space
        # unless we are passed a point of projective Berkovich space
        Element.__init__(self, parent)
        self._p = parent.prime()
        self._base_space = parent.base()
        self._base_type = parent._base_type
        self._ideal = parent._ideal

        # if this is a point of projective Berkovich space, we raise an error
        if isinstance(center, Berkovich_Element_Cp_Projective):
            raise TypeError('use as_affine_point to convert to affine Berkovich space')

        Berkovich_Element_Cp.__init__(self, parent=parent, center=center, radius=radius, power=power,
                                      prec=prec, space_type='affine', error_check=error_check)

    def as_projective_point(self):
        r"""
        Return the corresponding point of projective Berkovich space.

        We identify affine Berkovich space with the subset `P^1_{\text{Berk}}(C_p) - \{(1 : 0)\}`.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(5)
            sage: B(5).as_projective_point()
            Type I point centered at (5 + O(5^21) : 1 + O(5^20))

        ::

            sage: B(0, 1).as_projective_point()
            Type II point centered at (0 : 1 + O(5^20)) of radius 5^0

        ::

            sage: L.<t> = PolynomialRing(Qp(5))
            sage: T = FractionField(L)
            sage: f = T(1/t)
            sage: R.<x> = RR[]
            sage: Y = FractionField(R)
            sage: g = (40*pi)/x                                                         # needs sage.symbolic
            sage: Q2 = B(f, g)                                                          # needs sage.symbolic
            sage: Q2.as_projective_point()                                              # needs sage.symbolic
            Type IV point of precision 20 with centers given by (1 + O(5^20))/((1 + O(5^20))*t)
            and radii given by 40.0000000000000*pi/x
        """
        from sage.schemes.berkovich.berkovich_space import Berkovich_Cp_Projective
        new_space = Berkovich_Cp_Projective(self.parent().base_ring(), self.parent().ideal())
        if self.type_of_point() == 1:
            return new_space(self.center())
        elif self.type_of_point() == 2:
            return new_space(self.center(), power=self.power())
        elif self.type_of_point() == 3:
            return new_space(self.center(), self.radius())
        if self._center_func is None:
            center = self.center()
        else:
            center = self.center_function()
        if self._radius_func is None:
            radius = self.radius()
        else:
            radius = self.radius_function()
        return new_space(center, radius, prec=self.prec())

    def __eq__(self, other):
        """
        Equality operator.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(1, RR(3**(1/2)))                                               # needs sage.symbolic
            sage: Q2 = B(1, 3**(1/2))                                                   # needs sage.symbolic
            sage: Q1 == Q2                                                              # needs sage.symbolic
            True

        ::

            sage: Q3 = B(1)
            sage: Q4 = B(4)
            sage: Q3 == Q4
            False

        ::

            sage: Q5 = B(1, 4)
            sage: Q1 == Q5                                                              # needs sage.symbolic
            False

        ::

            sage: Q1 == Q3                                                              # needs sage.symbolic
            False
        """
        if other is self:
            return True
        if not isinstance(other, Berkovich_Element_Cp_Affine):
            return False
        if other.parent() != self.parent():
            return False
        stype = self.type_of_point()
        otype = other.type_of_point()
        if stype == otype and stype == 1:
            return self.center() == other.center()
        elif stype == otype and stype == 4:
            raise NotImplementedError("Equality for type IV points not yet implemented")
        elif stype in [2, 3] and otype in [2, 3]:
            if self.radius() != other.radius():
                return False
            center_dist = self._custom_abs(self.center() - other.center())
            return center_dist <= self.radius()
        else:
            return False

    def __hash__(self):
        """
        Return the hash of this point.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(1, RR(3**(1/2)))                                               # needs sage.symbolic
            sage: Q2 = B(1, 3**(1/2))                                                   # needs sage.symbolic
            sage: hash(Q1) == hash(Q2)                                                  # needs sage.symbolic
            True

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^3 + 20)
            sage: ideal = A.ideal(-1/2*a^2 + a - 3)
            sage: B = Berkovich_Cp_Projective(A, ideal)
            sage: Q1 = B(a^2 + 1, 2)
            sage: Q2 = B(0, 2)
            sage: hash(Q1) == hash(Q2)
            True
        """
        if self.type_of_point() == 1:
            return hash(self.center())
        elif self.type_of_point() == 4:
            raise NotImplementedError('hash not defined for type IV points')
        return hash(self.radius())

    def lt(self, other):
        r"""
        Return ``True`` if this point is strictly less than ``other`` in the standard partial order.

        Roughly, the partial order corresponds to containment of
        the corresponding disks in ``Cp``.

        For example, let x and y be points of type II or III.
        If x has center `c_1` and radius `r_1` and y has center
        `c_2` and radius `r_2`, `x < y` if and only if `D(c_1,r_1)`
        is a subset of `D(c_2,r_2)` in `\CC_p`.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point

        OUTPUT:

        - ``True`` -- if this point is less than ``other`` in the standard partial order
        - ``False`` -- otherwise

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(5, 0.5)
            sage: Q2 = B(5, 1)
            sage: Q1.lt(Q2)
            True

        ::

            sage: Q3 = B(1)
            sage: Q1.lt(Q3)
            False

        TESTS::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(5)
            sage: Q1.lt(Q1)
            False

        ::

            sage: Q2 = B([4, 1/3], [5, 1])
            sage: Q1.lt(Q2)
            False

        ::

            sage: Q4 = B(0, 1)
            sage: Q1.lt(Q4)
            True

        ::

            sage: Q2.lt(Q4)
            False
        """
        if not isinstance(other, Berkovich_Element_Cp_Affine):
            raise TypeError('other must be a point of a projective Berkovich space, but was %s' % other)
        if self.parent() != other.parent():
            raise ValueError('other must be a point of the same projective Berkovich space')

        if self == other:
            return False
        if other.type_of_point() in [1, 4]:
            return False

        if self.type_of_point() == 4:
            center = self.center()[-1]
            dist = self._custom_abs(other.center() - center)
            return dist <= other.radius() and self.radius()[-1] <= other.radius()
        else:
            dist = self._custom_abs(self.center() - other.center())
            return dist <= other.radius() and self.radius() <= other.radius()

    def gt(self, other):
        r"""
        Return ``True`` if this point is strictly greater than ``other`` in the standard partial order.

        Roughly, the partial order corresponds to containment of
        the corresponding disks in `\CC_p`.

        For example, let x and y be points of type II or III.
        If x has center `c_1` and radius `r_1` and y has center
        `c_2` and radius `r_2`, `x < y` if and only if `D(c_1,r_1)`
        is a subset of `D(c_2,r_2)` in `\CC_p`.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point

        OUTPUT:

        - ``True`` -- if this point is greater than ``other`` in the standard partial order
        - ``False`` -- otherwise

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(QQ, 3)
            sage: Q1 = B(5, 3)
            sage: Q2 = B(5, 1)
            sage: Q1.gt(Q2)
            True

        ::

            sage: Q3 = B(1/27)
            sage: Q1.gt(Q3)
            False

        TESTS::

            sage: B = Berkovich_Cp_Affine(QQ, 3)
            sage: Q1 = B(5)
            sage: Q1.gt(Q1)
            False

        ::

            sage: Q2 = B(0, 1)
            sage: Q1.gt(Q2)
            False

        ::

            sage: Q3 = B([0, 3], [5, 1])
            sage: Q2.gt(Q3)
            True
        """
        if not isinstance(other, Berkovich_Element_Cp_Affine):
            raise TypeError('other must be a point of a projective Berkovich space, but was %s' % other)
        if self.parent() != other.parent():
            raise ValueError('other must be a point of the same projective Berkovich space')

        if self == other:
            return False
        if self.type_of_point() in [1, 4]:
            return False

        if other.type_of_point() == 4:
            center = other.center()[-1]
            dist = self._custom_abs(self.center() - center)
            return dist <= self.radius() and other.radius()[-1] <= self.radius()
        else:
            dist = self._custom_abs(self.center() - other.center())
            return dist <= self.radius() and other.radius() <= self.radius()

    def join(self, other, basepoint=Infinity):
        """
        Compute the join of this point and ``other`` with respect to ``basepoint``.

        The join is first point that lies on the intersection
        of the path from this point to ``basepoint`` and the path from ``other`` to
        ``basepoint``.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point
        - ``basepoint`` -- (default: ``Infinity``) a point of the same
          Berkovich space as this point or ``Infinity``

        OUTPUT: a point of the same Berkovich space

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(2, 1)
            sage: Q2 = B(2, 2)
            sage: Q1.join(Q2)
            Type III point centered at 2 + O(3^20) of radius 2.00000000000000

        ::

            sage: Q3 = B(5)
            sage: Q3.join(Q1)
            Type II point centered at 2 + 3 + O(3^20) of radius 3^0

        ::

            sage: Q3.join(Q1, basepoint=Q2)
            Type II point centered at 2 + O(3^20) of radius 3^0

        TESTS::

            sage: Q4 = B(1/3**8 + 2, 1)
            sage: Q2.join(Q4, basepoint=Q1)
            Type III point centered at 2 + O(3^20) of radius 2.00000000000000

        ::

            sage: Q5 = B(2, 1/9)
            sage: Q6 = B(1, 1/27)
            sage: Q4.join(Q5, basepoint=Q6)
            Type II point centered at 1 + O(3^20) of radius 3^0

        ::

            sage: Q7 = B(1/27, 1/27)
            sage: Q1.join(Q7, Q2)
            Type III point centered at 2 + O(3^20) of radius 2.00000000000000
        """
        # we error check and then pass to projective space to do the join
        if not isinstance(other, Berkovich_Element_Cp_Affine):
            raise TypeError('other must be a point of affine Berkovich space. other was %s' % other)
        if self.parent() != other.parent():
            raise ValueError('other must be a point of the same affine Berkovich space')
        if self.type_of_point() == 4 or other.type_of_point() == 4:
            raise NotImplementedError("join with type IV points not implemented")

        proj_self = self.as_projective_point()
        proj_other = other.as_projective_point()

        if basepoint == Infinity:
            return proj_self.join(proj_other).as_affine_point()

        if not isinstance(basepoint, Berkovich_Element_Cp_Affine):
            raise TypeError('basepoint must a point of affine Berkovich space. basepoint was %s' % basepoint)
        if basepoint.parent() != self.parent():
            raise ValueError("basepoint must be a point of the same affine Berkovich space")
        if basepoint.type_of_point() == 4:
            raise NotImplementedError("join not implemented for type IV basepoint")
        proj_basepoint = basepoint.as_projective_point()
        return proj_self.join(proj_other, proj_basepoint).as_affine_point()

    def involution_map(self):
        r"""
        Return the image of this point under the involution map.

        The involution map is the extension of the map ``z |-> 1/z``
        on `\CC_p` to Berkovich space.

        For affine Berkovich space, not defined for the type I
        point centered at 0.

        If zero is contained in every disk approximating a type IV point,
        then the image under the involution map is not defined. To avoid
        this error, increase precision.

        OUTPUT: a point of the same Berkovich space

        EXAMPLES:

        The involution map is 1/z on type I points::

            sage: B = Berkovich_Cp_Affine(3)
            sage: Q1 = B(1/2)
            sage: Q1.involution_map()
            Type I point centered at 2 + O(3^20)

        ::

            sage: Q2 = B(0, 1/3)
            sage: Q2.involution_map()
            Type II point centered at 0 of radius 3^1

        ::

            sage: Q3 = B(1/3, 1/3)
            sage: Q3.involution_map()
            Type II point centered at 3 + O(3^21) of radius 3^-3

        TESTS::

            sage: B = Berkovich_Cp_Affine(3)
            sage: B(0).involution_map()
            Traceback (most recent call last):
            ...
            ValueError: involution map not defined on affine type I point centered at 0

        ::

            sage: B(1/81, 1.5).involution_map()
            Type III point centered at 3^4 + O(3^24) of radius 0.000228623685413809

        ::

            sage: B([1, 2], [3, 1]).involution_map()
            Traceback (most recent call last):
            ...
            ValueError: precision of type IV is not high enough to define image

        ::

            sage: B([1/81, 10/81], [10, 9]).involution_map()
            Type IV point of precision 2, approximated by disks centered at [3^4 + O(3^24),
            3^4 + 2*3^6 + 2*3^7 + 2*3^10 + 2*3^11 + 2*3^14 + 2*3^15 + 2*3^18 + 2*3^19 + 2*3^22
            + 2*3^23 + O(3^24)] ... with radii [0.00152415790275873, 0.00137174211248285] ...
        """
        if self.type_of_point() == 1:
            if self.center() == 0:
                raise ValueError("involution map not defined on affine type I point centered at 0")
            return self.parent()(1 / self.center())

        zero = self.parent()(ZZ(0))
        radius = self.radius()

        if self.type_of_point() in [2, 3]:
            zero_contained_in_self = self.gt(zero)
            if zero_contained_in_self:
                if self.type_of_point() == 2:
                    power = self.power()
                    return self.parent()(ZZ(0), power=-power)
                return self.parent()(ZZ(0), RR(1 / radius))
            return self.parent()(1 / self.center(), RR(radius / (self._custom_abs(self.center())**2)))

        new_center_lst = []
        new_radius_lst = []
        for i in range(len(self.center())):
            berk_point = self.parent()(self.center()[i], self.radius()[i])
            zero_check = berk_point.gt(zero)
            if zero_check:
                continue
            else:
                new_center = 1 / self.center()[i]
                new_radius = self.radius()[i] / (self._custom_abs(self.center()[i])**2)
                new_center_lst.append(new_center)
                new_radius_lst.append(new_radius)
        if not new_center_lst:
            raise ValueError('precision of type IV is not high enough to define image')
        return self.parent()(new_center_lst, new_radius_lst, error_check=False)

    def contained_in_interval(self, start, end):
        """
        Check if this point is an element of the interval [``start``, ``end``].

        INPUT:

        - ``start`` -- a point of the same Berkovich space as this point
        - ``end`` -- a point of the same Berkovich space as this point

        OUTPUT:

        - ``True`` if this point is an element of [``start``, ``end``].
        - ``False`` otherwise.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective((3))
            sage: Q1 = B(2, 1)
            sage: Q2 = B(2, 4)
            sage: Q3 = B(1/3)
            sage: Q2.contained_in_interval(Q1, Q3.join(Q1))
            False

        ::

            sage: Q4 = B(1/81, 1)
            sage: Q2.contained_in_interval(Q1, Q4.join(Q1))
            True
        """
        if not isinstance(start, Berkovich_Element_Cp_Affine):
            raise TypeError("start must be a point of affine Berkovich space. start was %s" % start)
        if start.parent() != self.parent():
            raise ValueError("start must be a point of the same Berkovich space as this point")
        if not isinstance(end, Berkovich_Element_Cp_Affine):
            raise TypeError("end must be a point of affine Berkovich space. end was %s" % end)
        if end.parent() != self.parent():
            raise ValueError("end must be a point of the same Berkovich space as this point")

        proj_self = self.as_projective_point()
        proj_start = start.as_projective_point()
        proj_end = end.as_projective_point()
        return proj_self.contained_in_interval(proj_start, proj_end)


class Berkovich_Element_Cp_Projective(Berkovich_Element_Cp):
    r"""
    Element class of the Berkovich projective line over `\CC_p`.

    Elements are categorized into four types, represented by specific data:

    - Type I points are represented by a center in the ``base`` of the parent Berkovich space,
      which is projective space of dimension 1 over either `\QQ_p`, a finite extension of `\QQ_p`,
      or a number field.

    - Type II points are represented by a center in the ``base`` of the parent Berkovich space,
      and a rational power of `p`.

    - Type III points are represented by a center in the ``base`` of the parent Berkovich space,
      and by a radius, a real number, in `[0,\infty)`.

    - Type IV points are represented by a finite list of centers in the ``base`` of the parent
      Berkovich space and a finite list of radii in `[0,\infty)`.

    The projective Berkovich line is viewed as the one-point compactification of
    the affine Berkovich line. The projective Berkovich line therefore contains
    every point of the affine Berkovich line, along with a type I point centered
    at infinity.

    INPUT:

    - ``center`` -- for type I, II, and III points, the center of the
      corresponding disk in `P^1(\CC_p)`. If the parent Berkovich space was created using a number field
      `K`, then ``center`` can be an element of `P^1(K)`. Otherwise, ``center``
      must be an element of a projective space of dimension 1 over a `p`-adic field.
      For type IV points, can be a list of centers used to approximate the point or a
      univariate function that computes the centers (computation starts at 1).

    - ``radius`` -- (optional) For type I, II, and III points, the radius of the
      corresponding disk in `\CC_p`. Must coerce into the real numbers. For type IV points,
      can be a list of radii used to approximate the point or a univariate function that
      computes the radii (computation starts at 1).

    - ``power`` -- (optional) Rational number. Used for constructing type II points; specifies
      the power of ``p`` such that `p^\text{power}` = radius

    - ``prec`` -- (default: 20) the number of disks to be used to approximate a type IV point

    - ``error_check`` -- boolean (default: ``True``); if error checking should be run on input. If
      input is correctly formatted, can be set to ``False`` for better performance.
      WARNING: with error check set to ``False``, any error in the input will lead to
      incorrect results.

    EXAMPLES:

    Type I points can be created by specifying the corresponding point of `P^1(\CC_p)`::

        sage: S = ProjectiveSpace(Qp(5), 1)
        sage: P = Berkovich_Cp_Projective(S); P
        Projective Berkovich line over Cp(5) of precision 20

    ::

        sage: a = S(0, 1)
        sage: Q1 = P(a); Q1
        Type I point centered at (0 : 1 + O(5^20))

    ::

        sage: Q2 = P((1,0)); Q2
        Type I point centered at (1 + O(5^20) : 0)

    Type II and III points can be created by specifying a center and a radius::

        sage: Q3 = P((0,5), 5**(3/2)); Q3                                               # needs sage.symbolic
        Type II point centered at (0 : 1 + O(5^20)) of radius 5^3/2

    ::

        sage: Q4 = P(0, 3**(3/2)); Q4                                                   # needs sage.symbolic
        Type III point centered at (0 : 1 + O(5^20)) of radius 5.19615242270663

    Type IV points can be created from lists of centers and radii::

        sage: b = S((3,2))  # create centers
        sage: c = S((4,3))
        sage: d = S((2,3))
        sage: L = [b, c, d]
        sage: R = [1.761, 1.123, 1.112]
        sage: Q5 = P(L, R); Q5
        Type IV point of precision 3, approximated by disks centered at
        [(4 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + 2*5^10 +
         2*5^11 + 2*5^12 + 2*5^13 + 2*5^14 + 2*5^15 + 2*5^16 + 2*5^17 + 2*5^18 + 2*5^19 + O(5^20) :
         1 + O(5^20)), (3 + 3*5 + 5^2 + 3*5^3 + 5^4 + 3*5^5 + 5^6 + 3*5^7 + 5^8 + 3*5^9 +
         5^10 + 3*5^11 + 5^12 + 3*5^13 + 5^14 + 3*5^15 + 5^16 + 3*5^17 + 5^18 + 3*5^19 + O(5^20) :
         1 + O(5^20))] ... with radii [1.76100000000000, 1.12300000000000] ...

    Type IV points can also be created from univariate functions. Since the centers of
    the sequence of disks can not be the point at infinity in `P^1(\CC_p)`, only functions
    into `\CC_p` are supported::

        sage: L.<t> = PolynomialRing(Qp(5))
        sage: T = FractionField(L)
        sage: f = T(1/t)
        sage: R.<x> = RR[]
        sage: Y = FractionField(R)
        sage: g = (40*pi)/x                                                             # needs sage.symbolic
        sage: Q6 = P(f, g); Q6                                                          # needs sage.symbolic
        Type IV point of precision 20 with centers given by (1 + O(5^20))/((1 + O(5^20))*t)
         and radii given by 40.0000000000000*pi/x

    TESTS::

        sage: P((1,0), 3)
        Traceback (most recent call last):
        ...
        ValueError: type II and III points can not be centered at infinity

        sage: B = Berkovich_Cp_Projective(3)
        sage: Q1 = B(3)
        sage: TestSuite(Q1).run()
    """

    def __init__(self, parent, center, radius=None, power=None, prec=20, error_check=True):
        """
        Initialization function.

        EXAMPLES::

            sage: S = ProjectiveSpace(Qp(7), 1)
            sage: P = Berkovich_Cp_Projective(S)
            sage: P(0,1)
            Type II point centered at (0 : 1 + O(7^20)) of radius 7^0
        """
        # if we are given a point of Affine Berkovich Space, we do the conversion
        # otherwise we call the Berkovich_Element_Cp constructor with space_type="projective"
        Element.__init__(self, parent)
        self._p = parent.prime()
        self._base_space = parent.base()
        self._base_type = parent._base_type
        self._ideal = parent._ideal

        # conversion from Affine points is handled in this constructor
        if isinstance(center, Berkovich_Element_Cp_Affine):
            raise TypeError('use as_projective_point to convert to projective Berkovich space')

        Berkovich_Element_Cp.__init__(self, parent=parent, center=center, radius=radius, power=power,
                                      prec=prec, space_type='projective', error_check=error_check)

    def as_affine_point(self):
        """
        Return the corresponding affine point after dehomogenizing at infinity.

        OUTPUT: a point of affine Berkovich space

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(5)
            sage: B(5).as_affine_point()
            Type I point centered at 5 + O(5^21)

        ::

            sage: Q = B(0, 1).as_affine_point(); Q
            Type II point centered at 0 of radius 5^0
            sage: Q.parent()
            Affine Berkovich line over Cp(5) of precision 20

        ::

            sage: L.<t> = PolynomialRing(Qp(5))
            sage: T = FractionField(L)
            sage: f = T(1/t)
            sage: R.<x> = RR[]
            sage: Y = FractionField(R)
            sage: g = (40*pi)/x                                                         # needs sage.symbolic
            sage: Q2 = B(f, g)                                                          # needs sage.symbolic
            sage: Q2.as_affine_point()                                                  # needs sage.symbolic
            Type IV point of precision 20 with centers given by (1 + O(5^20))/((1 + O(5^20))*t)
            and radii given by 40.0000000000000*pi/x
        """
        if self.center()[1] == 0:
            raise ValueError('cannot convert infinity to affine Berkovich space')
        from sage.schemes.berkovich.berkovich_space import Berkovich_Cp_Affine
        new_space = Berkovich_Cp_Affine(self.parent().base_ring(), self.parent().ideal())
        if self.type_of_point() in [1, 2, 3]:
            center = self.center()[0]
            if self.type_of_point() == 1:
                return new_space(center)
            elif self.type_of_point() == 2:
                return new_space(center, power=self.power())
            elif self.type_of_point() == 3:
                return new_space(center, self.radius())
        if self._center_func is None:
            center = [i[0] for i in self.center()]
        else:
            center = self.center_function()
        if self._radius_func is None:
            radius = self.radius()
        else:
            radius = self.radius_function()
        return new_space(center, radius, prec=self.prec())

    def __eq__(self, other):
        """
        Equality operator.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B([2, 2], RR(3**(1/2)))                                          # needs sage.symbolic
            sage: Q2 = B([1, 1], 3**(1/2))                                              # needs sage.symbolic
            sage: Q1 == Q2                                                              # needs sage.symbolic
            True

        ::

            sage: Q3 = B(1)
            sage: Q4 = B(4)
            sage: Q3 == Q4
            False

        ::

            sage: Q5 = B(1, 4)
            sage: Q1 == Q5                                                              # needs sage.symbolic
            False

        ::

            sage: Q1 == Q3                                                              # needs sage.symbolic
            False
        """
        if other is self:
            return True
        if not isinstance(other, Berkovich_Element_Cp_Projective):
            return False
        if other.parent() != self.parent():
            return False
        stype = self.type_of_point()
        otype = other.type_of_point()
        if stype == otype and stype == 1:
            return self.center() == other.center()
        elif stype == otype and stype == 4:
            raise NotImplementedError("equality for type IV points not implemented")
        elif stype in [2, 3] and otype in [2, 3]:
            if self.radius() != other.radius():
                return False
            scent = self.center()[0]
            ocent = other.center()[0]
            center_dist = self._custom_abs(scent - ocent)
            return center_dist <= self.radius()
        else:
            return False

    def __hash__(self):
        """
        Return the hash of this point.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: P = ProjectiveSpace(B.base_ring(), 1)
            sage: Q1 = B(P.point([2, 2], False), RR(3**(1/2)))                          # needs sage.symbolic
            sage: Q2 = B([1, 1], 3**(1/2))                                              # needs sage.symbolic
            sage: hash(Q1) == hash(Q2)                                                  # needs sage.symbolic
            True

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^3 + 20)
            sage: ideal = A.ideal(-1/2*a^2 + a - 3)
            sage: B = Berkovich_Cp_Projective(A, ideal)
            sage: Q1 = B(a^2 + 1, 2)
            sage: Q2 = B(0, 2)
            sage: hash(Q1) == hash(Q2)
            True
        """
        if self.type_of_point() == 1:
            return hash(self.center())
        elif self.type_of_point() == 4:
            raise ValueError('hash not defined for type IV points')
        return hash(self.radius())

    def lt(self, other):
        r"""
        Return ``True`` if this point is strictly less than ``other`` in the standard partial order.

        Roughly, the partial order corresponds to containment of
        the corresponding disks in `\CC_p`.

        For example, let x and y be points of type II or III.
        If x has center `c_1` and radius `r_1` and y has center
        `c_2` and radius `r_2`, `x < y` if and only if `D(c_1,r_1)`
        is a subset of `D(c_2,r_2)` in `\CC_p`.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point

        OUTPUT:

        - ``True`` -- if this point is less than ``other`` in the standard partial order
        - ``False`` -- otherwise

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(5, 0.5)
            sage: Q2 = B(5, 1)
            sage: Q1.lt(Q2)
            True

        ::

            sage: Q3 = B(1)
            sage: Q1.lt(Q3)
            False

        TESTS::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(5)
            sage: Q1.lt(Q1)
            False

        ::

            sage: Q2 = B([4, 1/3], [5, 1])
            sage: Q1.lt(Q2)
            False

        ::

            sage: Q3 = B((1,0))
            sage: Q4 = B(0, 1)
            sage: Q3.lt(Q4)
            False

        ::

            sage: Q4.lt(Q3)
            True

        ::

            sage: Q1.lt(Q4)
            True

        ::

            sage: Q2.lt(Q4)
            False
        """
        if not isinstance(other, Berkovich_Element_Cp_Projective):
            raise TypeError('other must be a point of a projective Berkovich space, but was %s' % other)
        if self.parent() != other.parent():
            raise ValueError('other must be a point of the same projective Berkovich space')

        if self == other:
            return False

        # infinity is maximal with respect to the standard partial order
        infinity = self.parent()((1, 0))
        if self == infinity:
            return False
        if other == infinity:
            return True

        if other.type_of_point() in [1, 4]:
            return False
        if self.type_of_point() == 4:
            center = self.center()[-1]
            dist = self._custom_abs(other.center()[0] - center[0])
            return dist <= other.radius() and self.radius()[-1] <= other.radius()
        else:
            dist = self._custom_abs(self.center()[0] - other.center()[0])
            return dist <= other.radius() and self.radius() <= other.radius()

    def gt(self, other):
        r"""
        Return ``True`` if this point is strictly greater than ``other`` in the standard partial order.

        Roughly, the partial order corresponds to containment of
        the corresponding disks in `\CC_p`.

        For example, let x and y be points of type II or III.
        If x has center `c_1` and radius `r_1` and y has center
        `c_2` and radius `r_2`, `x < y` if and only if `D(c_1, r_1)`
        is a subset of `D(c_2, r_2)` in `\CC_p`.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point

        OUTPUT:

        - ``True`` -- if this point is greater than ``other`` in the standard partial order
        - ``False`` -- otherwise

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(QQ, 3)
            sage: Q1 = B(5, 3)
            sage: Q2 = B(5, 1)
            sage: Q1.gt(Q2)
            True

        ::

            sage: Q3 = B(1/27)
            sage: Q1.gt(Q3)
            False

        TESTS::

            sage: B = Berkovich_Cp_Projective(QQ, 3)
            sage: Q1 = B(5)
            sage: Q1.gt(Q1)
            False

        ::

            sage: Q2 = B(0, 1)
            sage: Q1.gt(Q2)
            False

        ::

            sage: Q3 = B([0, 3], [5, 1])
            sage: Q2.gt(Q3)
            True

        ::

            sage: Q4 = B((1,0))
            sage: Q4.gt(Q2)
            True

        ::

            sage: Q1.gt(Q4)
            False
        """
        if not isinstance(other, Berkovich_Element_Cp_Projective):
            raise TypeError('other must be a point of a projective Berkovich space, but was %s' % other)
        if self.parent() != other.parent():
            raise ValueError('other must be a point of the same projective Berkovich space')

        if self == other:
            return False
        # infinity is maximal with respect to the standard partial order
        infinity = self.parent()((1, 0))
        if self == infinity:
            return True
        if other == infinity:
            return False

        if self.type_of_point() in [1, 4]:
            return False
        if other.type_of_point() == 4:
            center = other.center()[-1]
            dist = self._custom_abs(self.center()[0] - center[0])
            return dist <= self.radius() and other.radius()[-1] <= self.radius()
        else:
            dist = self._custom_abs(self.center()[0] - other.center()[0])
            return dist <= self.radius() and other.radius() <= self.radius()

    def join(self, other, basepoint=Infinity):
        """
        Compute the join of this point and ``other``, with respect to ``basepoint``.

        The join is first point that lies on the intersection
        of the path from this point to ``basepoint`` and the path from ``other`` to
        ``basepoint``.

        INPUT:

        - ``other`` -- a point of the same Berkovich space as this point
        - ``basepoint`` -- (default: ``Infinity``) a point of the same
          Berkovich space as this point, or ``Infinity``

        OUTPUT: a point of the same Berkovich space

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(2, 1)
            sage: Q2 = B(2, 2)
            sage: Q1.join(Q2)
            Type III point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 2.00000000000000

        ::

            sage: Q3 = B(5)
            sage: Q3.join(Q1)
            Type II point centered at (2 + 3 + O(3^20) : 1 + O(3^20)) of radius 3^0

        ::

            sage: Q3.join(Q1, basepoint=Q2)
            Type II point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 3^0

        TESTS::

            sage: Q4 = B(1/3**8 + 2, 1)
            sage: Q2.join(Q4, basepoint=Q1)
            Type III point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 2.00000000000000
            sage: Q5 = B(2, 1/9)
            sage: Q6 = B(1, 1/27)
            sage: Q4.join(Q5, basepoint=Q6)
            Type II point centered at (1 + O(3^20) : 1 + O(3^20)) of radius 3^0
            sage: Q7 = B(1/27, 1/27)
            sage: Q1.join(Q7, Q2)
            Type III point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 2.00000000000000
            sage: Q1.join(Q2, Q7)
            Type III point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 2.00000000000000
            sage: Q8 = B(0, power=1/3)
            sage: Q9 = B(0, power=1/2)
            sage: Q8.join(Q9)
            Type II point centered at (0 : 1 + O(3^20)) of radius 3^1/2

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^3 + 20)
            sage: ideal = A.prime_above(3)
            sage: C = Berkovich_Cp_Projective(A, ideal)
            sage: Q10 = C(a, 1/9)
            sage: Q10.join(Q9)
            Traceback (most recent call last):
            ...
            ValueError: other must be a point of the same projective Berkovich line
            sage: Q11 = C(0, 1/3)
            sage: Q11.join(Q10)
            Type II point centered at (0 : 1) of radius 3^0
        """
        if not isinstance(other, Berkovich_Element_Cp_Projective):
            raise TypeError('other must be a point of a projective Berkovich line, instead was %s' % other)
        if other.parent() != self.parent():
            raise ValueError('other must be a point of the same projective Berkovich line')

        # if either self or other is type IV, we use the last disk in the approximation
        if self.type_of_point() == 4:
            new_center = self.center()[-1]
            new_radius = self.radius()[-1]
            return self.parent()(new_center, new_radius).join(other)
        if other.type_of_point() == 4:
            new_center = other.center()[-1]
            new_radius = other.radius()[-1]
            return self.join(self.parent()(new_center, new_radius))

        # we deal with the point at infinity as a special case
        infty = self.parent()((1, 0))

        if basepoint == Infinity or basepoint == infty:
            if self == infty or other == infty:
                return infty
            dist = self._custom_abs(self.center()[0] - other.center()[0])
            maximum = max(dist, self.radius(), other.radius())
            # optimize for when self or other are type II
            if maximum == self.radius() and self.type_of_point() == 2:
                return self.parent()(self.center(), power=self.power())
            if maximum == other.radius() and other.type_of_point() == 2:
                return self.parent()(self.center(), power=other.power())
            return self.parent()(self.center(), maximum)

        if not isinstance(basepoint, Berkovich_Element_Cp_Projective):
            raise TypeError('basepoint must be a point of a projective Berkovich line, instead was %s' % basepoint)
        if basepoint.parent() != self.parent():
            raise ValueError("basepoint must be a point of the same Berkovich projective line")

        # if the basepoint is type IV, we use the last disk in the approximation
        if basepoint.type_of_point() == 4:
            new_center = other.center()[-1]
            new_radius = other.radius()[-1]
            return self.join(other, self.parent()(new_center, new_radius))

        if self == infty:
            return other.join(basepoint)
        if other == infty:
            return self.join(basepoint)

        b_ge_s = basepoint.gt(self) or basepoint == self
        b_lt_s = basepoint.lt(self)
        b_ge_o = basepoint.gt(other) or basepoint == other
        b_lt_o = basepoint.lt(other)
        s_ge_o = self.gt(other) or self == other
        s_lt_o = self.lt(other)

        # we deal with all the cases where self and other are not comparable first
        if not (s_lt_o or s_ge_o):
            if not (b_ge_o or b_lt_o):
                if not (b_ge_s or b_lt_s):
                    # case where none of the points are comparable
                    dist_b_s = self._custom_abs(self.center()[0] - basepoint.center()[0])
                    dist_b_o = self._custom_abs(other.center()[0] - basepoint.center()[0])
                    return self.parent()(basepoint.center(),
                                         min(max(dist_b_o, other.radius(), basepoint.radius()),
                                             max(dist_b_s, self.radius(), basepoint.radius())))

                # case where self and basepoint are comparable
                else:
                    if b_ge_s:
                        return basepoint
                    else:
                        return self

            # case where other and basepoint are comparable
            else:
                if b_ge_o:
                    return basepoint
                else:
                    return other

        # now the cases where self > other
        elif s_ge_o:
            if not (b_ge_s or b_lt_s):
                return self
            if b_ge_s:
                return self
            if b_ge_o:
                return basepoint
            if b_lt_o:
                return other

        # join is symmetric, so we flip self and other so that self > other
        else:
            return other.join(self, basepoint)

    def involution_map(self):
        r"""
        Return the image of this point under the involution map.

        The involution map is the extension of the map ``z |-> 1/z``
        on `P^1(\CC_p)` to Berkovich space.

        If zero is contained in every disk approximating a type IV point,
        then the image under the involution map is not defined. To avoid
        this error, increase precision.

        OUTPUT: a point of the same Berkovich space

        EXAMPLES:

        The involution map is 1/z on type I points::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(1/2)
            sage: Q1.involution_map()
            Type I point centered at (2 + O(3^20) : 1 + O(3^20))

        ::

            sage: Q2 = B(0, 1/3)
            sage: Q2.involution_map()
            Type II point centered at (0 : 1 + O(3^20)) of radius 3^1

        ::

            sage: Q3 = B(1/3, 1/3)
            sage: Q3.involution_map()
            Type II point centered at (3 + O(3^21) : 1 + O(3^20)) of radius 3^-3

        TESTS::

            sage: B = Berkovich_Cp_Projective(3)
            sage: B((1,0)).involution_map()
            Type I point centered at (0 : 1 + O(3^20))

        ::

            sage: B(0).involution_map()
            Type I point centered at (1 + O(3^20) : 0)

        ::

            sage: B(1/81, 1.5).involution_map()
            Type III point centered at (3^4 + O(3^24) : 1 + O(3^20)) of radius 0.000228623685413809

        ::

            sage: B([1, 2], [3, 1]).involution_map()
            Traceback (most recent call last):
            ...
            ValueError: precision of type IV is not high enough to define image

        ::

            sage: B([1/81, 10/81], [10, 9]).involution_map()
            Type IV point of precision 2, approximated by disks centered at
            [(3^4 + O(3^24) : 1 + O(3^20)), (3^4 + 2*3^6 + 2*3^7 + 2*3^10 + 2*3^11 +
            2*3^14 + 2*3^15 + 2*3^18 + 2*3^19 + 2*3^22 + 2*3^23 + O(3^24) : 1 + O(3^20))]
            ... with radii [0.00152415790275873, 0.00137174211248285] ...
        """
        infty = self.parent()((1, 0))
        zero = self.parent()(0)

        if self.type_of_point() == 1:
            if self == infty:
                return zero
            if self == zero:
                return infty
            return self.parent()(1 / self.center()[0])

        if self.type_of_point() in [2, 3]:
            zero_contained_in_self = self.gt(zero)
            if zero_contained_in_self:
                if self.type_of_point() == 2:
                    power = self.power()
                    return self.parent()(ZZ(0), power=-power)
                return self.parent()(ZZ(0), 1 / self.radius())
            return self.parent()(1 / self.center()[0], self.radius() / (self._custom_abs(self.center()[0])**2))

        new_center_lst = []
        new_radius_lst = []
        for i in range(len(self.center())):
            berk_point = self.parent()(self.center()[i], self.radius()[i])
            zero_check = berk_point.gt(zero)
            if zero_check:
                continue
            else:
                new_center = 1 / self.center()[i][0]
                new_radius = self.radius()[i] / (self._custom_abs(self.center()[i][0])**2)
                new_center_lst.append(new_center)
                new_radius_lst.append(new_radius)
        if not new_center_lst:
            raise ValueError('precision of type IV is not high enough to define image')
        return self.parent()(new_center_lst, new_radius_lst)

    def contained_in_interval(self, start, end):
        """
        Check if this point is an element of the interval [``start``, ``end``].

        INPUT:

        - ``start`` -- a point of the same Berkovich space as this point
        - ``end`` -- a point of the same Berkovich space as this point

        OUTPUT:

        - ``True`` if this point is an element of [``start``, ``end``].
        - ``False`` otherwise.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(2, 1)
            sage: Q2 = B(2, 4)
            sage: Q3 = B(1/3)
            sage: Q2.contained_in_interval(Q1, Q3.join(Q1))
            False

        ::

            sage: Q4 = B(1/81, 1)
            sage: Q2.contained_in_interval(Q1, Q4.join(Q1))
            True

        TESTS::

            sage: B = Berkovich_Cp_Projective(3)
            sage: infty = B((1, 0))
            sage: zero = B(0)
            sage: gauss = B(0, 1)
            sage: infty.contained_in_interval(zero, gauss)
            False

        ::

            sage: Q1 = B(1, 3)
            sage: infty.contained_in_interval(gauss, Q1)
            False

        ::

            sage: zero.contained_in_interval(infty, gauss)
            False

        ::

            sage: gauss.contained_in_interval(zero, infty)
            True

        ::

            sage: Q2 = B(81, 1/3)
            sage: gauss.contained_in_interval(infty, Q2)
            True
        """
        if not isinstance(start, Berkovich_Element_Cp_Projective):
            raise TypeError("start must be a point of Berkovich space")
        if start.parent() != self.parent():
            raise ValueError("start must be a point of the same Berkovich space as this point")
        if not isinstance(end, Berkovich_Element_Cp_Projective):
            raise TypeError("start must be a point of Berkovich space")
        if end.parent() != self.parent():
            raise ValueError("start must be a point of the same Berkovich space as this point")

        # we treat infinity as a special case
        infty = self.parent()((1, 0))
        zero = self.parent()(ZZ(0))
        if self == infty:
            if start == zero or end == zero:
                return end == infty or start == infty
            return (self.involution_map()).contained_in_interval(start.involution_map(),
                                                                 end.involution_map())
        if start == infty or end == infty:
            if self == zero:
                return end == zero or start == zero
            if start == zero or end == zero:
                gauss = self.parent()(ZZ(0), ZZ(1))
                return self.contained_in_interval(start, gauss) or self.contained_in_interval(gauss, end)
            return self.involution_map().contained_in_interval(start.involution_map(),
                                                               end.involution_map())
        join = start.join(end)
        j_ge_s = join.gt(self) or join == self
        s_ge_start = self.gt(start) or self == start
        s_ge_end = self.gt(end) or self == end
        return j_ge_s and (s_ge_end or s_ge_start)
