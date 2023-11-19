"""
Univariate differential polynomial rings

This module provides the
:class:`~sage.rings.polynomial.differential_polynomial_ring.DifferentialPolynomialRing`.
In the class hierarchy in Sage, the locution *Differential Polynomial* is used
for a Ore polynomial without twisting morphism.

.. SEEALSO::

    :class:`~sage.rings.polynomial.ore_polynomial_ring.OrePolynomialRing`

AUTHOR:

- Xavier Caruso, Raphaël Pagès (2023-06): Initial version
"""

# ***************************************************************************
#    Copyright (C) 2023 Xavier Caruso <xavier.caruso@normalesup.org>
#                  2023 Raphaël Pagès <raphael.pages@math.u-bordeaux.fr>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

class DifferentialPolynomialRing(OrePolynomialRing):
    r"""
    Class for generic differential operators.

    TESTS::

        sage: R.<t> = ZZ[]
        sage: S.<D> = OrePolynomialRing(R, t*R.derivation())
        sage: TestSuite(S).run()
    """

    def __init__(self, base_ring, morphism, derivation, name, sparse, category=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``derivation`` -- a derivation of the base ring

        - ``name`` -- string or list of strings representing the name of
          the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``category`` -- a category

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: derivation = t*R.derivation()
            sage: S.<D> = OrePolynomialRing(R,derivation)
            sage: S.category()
            Category of algebras over Univariate Polynomial Ring in t over Integer Ring
            sage: D*t
            t*D + t

        We test that S has the good type::

            sage: type(S)
            <class 'sage.rings.polynomial.differential_polynomial_ring.DifferentialPolynomialRing_with_category'>
        """
        if morphism is not None:
            raise NotImplementedError
        if self.Element is None:
            import sage.rings.polynomial.differential_polynomial_element
            self.Element = sage.rings.polynomial.differential_polynomial_element.DifferentialPolynomial_generic_dense
        OrePolynomialRing.__init__(self, base_ring, None, derivation, name, sparse, category)
