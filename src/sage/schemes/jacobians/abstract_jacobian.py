"""
Jacobians of curves

This module defines the base class of Jacobians as an abstract scheme.

AUTHORS:

- William Stein (2005)
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.schemes import Jacobians
from sage.categories.fields import Fields
from sage.categories.schemes import AbelianVarieties
_Fields = Fields()
from sage.schemes.generic.scheme import Scheme, is_Scheme
from sage.structure.richcmp import richcmp_method, richcmp

_Fields = Fields()


def is_Jacobian(J):
    """
    Return ``True`` if `J` is of type ``Jacobian_generic``.

    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian, is_Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: J = Jacobian(C)
        sage: is_Jacobian(J)
        ...
        DeprecationWarning: Use Jacobian_generic directly
        See https://github.com/sagemath/sage/issues/35467 for details.
        True

    ::

        sage: E = EllipticCurve('37a1')
        sage: is_Jacobian(E)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(35467, "Use Jacobian_generic directly")
    return isinstance(J, Jacobian_generic)


def Jacobian(C):
    """
    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: Jacobian(C)
        Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
    """
    try:
        return C.jacobian()
    except (AttributeError, TypeError):
        return Jacobian_generic(C)


@richcmp_method
class Jacobian_generic(Scheme):
    """
    Base class for Jacobians of projective curves.

    The input must be a projective curve over a field.

    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: J = Jacobian(C); J
        Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
    """
    def __init__(self, C, category=None):
        """
        Initialize.

        TESTS::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve(x^3 + y^3 + z^3)
            sage: J = Jacobian_generic(C); J
            Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
            sage: type(J)
            <class 'sage.schemes.jacobians.abstract_jacobian.Jacobian_generic_with_category'>

        Note: this is an abstract parent, so we skip element tests::

            sage: TestSuite(J).run(skip=["_test_an_element",\
                                         "_test_elements",\
                                         "_test_elements_eq_reflexive",\
                                         "_test_elements_eq_symmetric",\
                                         "_test_elements_eq_transitive",\
                                         "_test_elements_neq",\
                                         "_test_some_elements"])

        ::

            sage: Jacobian_generic(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Argument (=Integer Ring) must be a scheme.
            sage: Jacobian_generic(P2)
            Traceback (most recent call last):
            ...
            ValueError: C (=Projective Space of dimension 2 over Rational Field)
            must have dimension 1.

        ::

            sage: P2.<x, y, z> = ProjectiveSpace(Zmod(6), 2)
            sage: C = Curve(x + y + z, P2)
            sage: Jacobian_generic(C)
            Traceback (most recent call last):
            ...
            TypeError: C (=Projective Plane Curve over Ring of integers modulo 6
            defined by x + y + z) must be defined over a field.
        """
        if not isinstance(C, Scheme):
            raise TypeError("Argument (=%s) must be a scheme." % C)
        if C.base_ring() not in _Fields:
            raise TypeError("C (=%s) must be defined over a field." % C)
        if C.dimension() != 1:
            raise ValueError("C (=%s) must have dimension 1." % C)
        self.__curve = C
        # I am not sure how to deal with non-smooth curves...
        # TODO: Isolate the smooth case into a Jacobian_AbelianVariety or something
        category = Jacobians(C.base_scheme()).or_subcategory(category) if C.is_smooth() else None
        Scheme.__init__(self, C.base_scheme(), category=category)

    def __richcmp__(self, J, op):
        """
        Compare the Jacobian ``self`` to `J`.  If `J` is a Jacobian, then
        ``self`` and `J` are equal if and only if their curves are equal.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J1 = Jacobian(Curve(x^3 + y^3 + z^3))
            sage: J1 == J1
            True
            sage: J1 == P2
            False
            sage: J1 != P2
            True
            sage: J2 = Jacobian(Curve(x + y + z))
            sage: J1 == J2
            False
            sage: J1 != J2
            True
        """
        if not isinstance(J, Jacobian_generic):
            return NotImplemented
        return richcmp(self.curve(), J.curve(), op)

    def _repr_(self):
        """
        Return a string representation of this Jacobian.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J = Jacobian(Curve(x^3 + y^3 + z^3)); J
            Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
            sage: J._repr_()
            'Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3'
        """
        return "Jacobian of %s" % self.__curve

    def _point(self):
        """
        Return the Hom-set from some affine scheme to ``self``.

        OUTPUT:

        This method always raises a :exc:`NotImplementedError`; it is
        only abstract.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J = Jacobian(Curve(x^3 + y^3 + z^3))
            sage: J._point()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def curve(self):
        """
        Return the curve of which ``self`` is the Jacobian.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J = Jacobian(Curve(x^3 + y^3 + z^3))
            sage: J.curve()
            Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
        """
        return self.__curve

    base_curve = curve

    def change_ring(self, R):
        r"""
        Return the Jacobian over the ring `R`.

        INPUT:

        - ``R`` -- a field; the new base ring

        OUTPUT: the Jacobian over the ring `R`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3 - 10*x + 9)
            sage: Jac = H.jacobian();   Jac
            Jacobian of Hyperelliptic Curve over Rational Field
             defined by y^2 = x^3 - 10*x + 9
            sage: Jac.change_ring(RDF)
            Jacobian of Hyperelliptic Curve over Real Double Field
             defined by y^2 = x^3 - 10.0*x + 9.0
        """
        return self.curve().change_ring(R).jacobian()

    def base_extend(self, R):
        r"""
        Return the natural extension of ``self`` over `R`.

        INPUT:

        - ``R`` -- a field; the new base field

        OUTPUT: the Jacobian over the ring `R`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3 - 10*x + 9)
            sage: Jac = H.jacobian();   Jac
            Jacobian of Hyperelliptic Curve over Rational Field
             defined by y^2 = x^3 - 10*x + 9

            sage: # needs sage.rings.number_field
            sage: F.<a> = QQ.extension(x^2 + 1)
            sage: Jac.base_extend(F)
            Jacobian of Hyperelliptic Curve over Number Field in a with defining
             polynomial x^2 + 1 defined by y^2 = x^3 - 10*x + 9
        """
        if R not in _Fields:
            raise ValueError('Not a field: ' + str(R))
        if self.base_ring() is R:
            return self
        if not R.has_coerce_map_from(self.base_ring()):
            raise ValueError('no natural map from the base ring (=%s) to R (=%s)!'
                             % (self.base_ring(), R))
        return self.change_ring(R)
