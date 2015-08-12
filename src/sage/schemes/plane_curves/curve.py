"""
Generic plane curves
"""

from sage.misc.all import latex


from sage.schemes.generic.algebraic_scheme import (
    AlgebraicScheme_subscheme, AlgebraicScheme_subscheme_projective)

from sage.schemes.generic.divisor_group import DivisorGroup

from sage.schemes.generic.divisor import Divisor_curve

class Curve_generic(AlgebraicScheme_subscheme):
    r"""
    EXAMPLES::

        sage: A.<x,y,z> = AffineSpace(QQ,3)
        sage: C = Curve([x-y,z-2])
        sage: loads(C.dumps()) == C
        True
    """

    def _repr_(self):
        """
        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ,3)
            sage: C = Curve([x-y,z-2])
            sage: C
            Affine Space Curve over Rational Field defined by x - y, z - 2

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: C = Curve(x-y)
            sage: C
            Projective Curve over Rational Field defined by x - y
        """
        return "%s Curve over %s defined by %s"%(
            self._repr_type(), self.base_ring(), ', '.join([str(x) for x in self.defining_polynomials()]))

    def _repr_type(self):
        return "Generic"

    def _latex_(self):
        """
        EXAMPLES::

            sage: x,y,z = PolynomialRing(QQ, 3, names='x,y,z').gens()
            sage: C = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: latex(C)
            - x^{3} + y^{2} z - 17 x z^{2} + y z^{2}

            sage: A2 = AffineSpace(2, QQ, names=['x','y'])
            sage: x, y = A2.coordinate_ring().gens()
            sage: C = Curve(y^2 - x^3 - 17*x + y)
            sage: latex(C)
            - x^{3} + y^{2} - 17 x + y
        """
        return latex(self.defining_polynomial())

    def defining_polynomial(self):
        """
        Return the defining polynomial of the curve.

        EXAMPLES::

            sage: x,y,z = PolynomialRing(QQ, 3, names='x,y,z').gens()
            sage: C = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: C.defining_polynomial()
            -x^3 + y^2*z - 17*x*z^2 + y*z^2
        """
        return self.defining_polynomials()[0]

    def divisor_group(self, base_ring=None):
        r"""
        Return the divisor group of the curve.

        INPUT:

        - ``base_ring`` -- the base ring of the divisor
          group. Usually, this is `\ZZ` (default) or `\QQ`.

        OUTPUT:

        The divisor group of the curve.

        EXAMPLES::

            sage: x,y,z = PolynomialRing(QQ, 3, names='x,y,z').gens()
            sage: C  = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: Cp = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: C.divisor_group() is Cp.divisor_group()
            True
        """
        return DivisorGroup(self, base_ring)

    def divisor(self, v, base_ring=None, check=True, reduce=True):
        r"""
        Return the divisor specified by ``v``.

        .. WARNING::

            The coefficients of the divisor must be in the base ring
            and the terms must be reduced. If you set ``check=False``
            and/or ``reduce=False`` it is your responsibility to pass
            a valid object ``v``.
        """
        return Divisor_curve(v, check=check, reduce=reduce, parent=self.divisor_group(base_ring))

    def genus(self):
        """
        The geometric genus of the curve.
        """
        return self.geometric_genus()

    def geometric_genus(self):
        r"""
        Return the geometric genus of the curve.

        This is by definition the
        genus of the normalization of the projective closure of the
        curve over the algebraic closure of the base field; the base
        field must be a prime field.

        .. NOTE::

            This calls Singular's genus command.

        EXAMPLE:

        Examples of projective curves. ::

            sage: P2 = ProjectiveSpace(2, GF(5), names=['x','y','z'])
            sage: x, y, z = P2.coordinate_ring().gens()
            sage: C = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: C.geometric_genus()
                  1
            sage: C = Curve(y^2*z - x^3)
            sage: C.geometric_genus()
                  0
            sage: C = Curve(x^10 + y^7*z^3 + z^10)
            sage: C.geometric_genus()
                  3

        Examples of affine curves. ::

            sage: x, y = PolynomialRing(GF(5), 2, 'xy').gens()
            sage: C = Curve(y^2 - x^3 - 17*x + y)
            sage: C.geometric_genus()
            1
            sage: C = Curve(y^2 - x^3)
            sage: C.geometric_genus()
            0
            sage: C = Curve(x^10 + y^7 + 1)
            sage: C.geometric_genus()
            3

        """
        try:
            return self.__genus
        except AttributeError:
            self.__genus = self.defining_ideal().genus()
            return self.__genus

    def union(self, other):
        """
        Return the union of ``self`` and ``other``.

        EXAMPLES::

            sage: x,y,z = PolynomialRing(QQ, 3, names='x,y,z').gens()
            sage: C1 = Curve(z - x)
            sage: C2 = Curve(y - x)
            sage: C1.union(C2).defining_polynomial()
            x^2 - x*y - x*z + y*z
        """
        from constructor import Curve
        return Curve(AlgebraicScheme_subscheme.union(self, other))

    __add__ = union

class Curve_generic_projective(Curve_generic, AlgebraicScheme_subscheme_projective):
    pass
