# sage_setup: distribution = sagemath-modules
r"""
Derivations of function fields: rational
"""

# ****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2011-2017 Julian Rüth <julian.rueth@gmail.com>
#                     2017      Alyson Deines
#                     2017-2019 Kwankyu Lee
#                     2018-2019 Travis Scrimshaw
#                     2019      Brent Baccala
#                     2022      Xavier Caruso
#                     2022      Frédéric Chapoton
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .derivations import FunctionFieldDerivation


class FunctionFieldDerivation_rational(FunctionFieldDerivation):
    """
    Derivations on rational function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: K.derivation()
        d/dx
    """
    def __init__(self, parent, u=None):
        """
        Initialize a derivation.

        INPUT:

        - ``parent`` -- the parent of this derivation

        - ``u`` -- a parameter describing the derivation

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: TestSuite(d).run()

        The parameter ``u`` can be the name of the variable::

            sage: K.derivation(x)
            d/dx

        or a list of length one whose unique element is the image
        of the generator of the underlying function field::

            sage: K.derivation([x^2])
            x^2*d/dx
        """
        FunctionFieldDerivation.__init__(self, parent)
        if u is None or u == parent.domain().gen():
            self._u = parent.codomain().one()
        elif u == 0 or isinstance(u, (list, tuple)):
            if u == 0 or len(u) == 0:
                self._u = parent.codomain().zero()
            elif len(u) == 1:
                self._u = parent.codomain()(u[0])
            else:
                raise ValueError("the length does not match")
        else:
            raise ValueError("you must pass in either a name of a variable or a list of coefficients")

    def _call_(self, x):
        """
        Compute the derivation of ``x``.

        INPUT:

        - ``x`` -- element of the rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d(x)  # indirect doctest
            1
            sage: d(x^3)
            3*x^2
            sage: d(1/x)
            -1/x^2
        """
        f = x.numerator()
        g = x.denominator()
        numerator = f.derivative() * g - f * g.derivative()
        if numerator.is_zero():
            return self.codomain().zero()
        else:
            v = numerator / g**2
            defining_morphism = self.parent()._defining_morphism
            if defining_morphism is not None:
                v = defining_morphism(v)
            return self._u * v

    def _add_(self, other):
        """
        Return the sum of this derivation and ``other``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d + d
            2*d/dx
        """
        return type(self)(self.parent(), [self._u + other._u])

    def _lmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d
            d/dx
            sage: x * d
            x*d/dx
        """
        return type(self)(self.parent(), [factor * self._u])
