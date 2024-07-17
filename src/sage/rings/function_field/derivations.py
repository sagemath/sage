r"""
Derivations of function fields

For global function fields, which have positive characteristics, the higher
derivation is available::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                          # needs sage.rings.function_field
    sage: h = L.higher_derivation()                                                     # needs sage.rings.function_field
    sage: h(y^2, 2)                                                                     # needs sage.rings.function_field
    ((x^7 + 1)/x^2)*y^2 + x^3*y

AUTHORS:

- William Stein (2010): initial version

- Julian Rüth (2011-09-14, 2014-06-23, 2017-08-21): refactored class hierarchy; added
  derivation classes; morphisms to/from fraction fields

- Kwankyu Lee (2017-04-30): added higher derivations and completions
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

from sage.rings.derivation import RingDerivationWithoutTwist


class FunctionFieldDerivation(RingDerivationWithoutTwist):
    r"""
    Base class for derivations on function fields.

    A derivation on `R` is a map `R \to R` with
    `D(\alpha+\beta)=D(\alpha)+D(\beta)` and `D(\alpha\beta)=\beta
    D(\alpha)+\alpha D(\beta)` for all `\alpha,\beta\in R`.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: d = K.derivation()
        sage: d
        d/dx
    """
    def __init__(self, parent):
        r"""
        Initialize a derivation.

        INPUT:

        - ``parent`` -- the differential module in which this
          derivation lives

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: TestSuite(d).run()
        """
        RingDerivationWithoutTwist.__init__(self, parent)
        self.__field = parent.domain()

    def is_injective(self) -> bool:
        r"""
        Return ``False`` since a derivation is never injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d.is_injective()
            False
        """
        return False

    def _rmul_(self, factor):
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
        return self._lmul_(factor)
