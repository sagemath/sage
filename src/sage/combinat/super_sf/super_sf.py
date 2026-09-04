r"""
Supersymmetric functions, with their realizations

AUTHORS:

- Shriya M
"""
# ****************************************************************************
#       Copyright (C) 2026 Shriya M <25shriya at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.rational_field import QQ
from sage.categories.fields import Fields
from sage.categories.rings import Rings
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.combinat.super_sf.powersum import SupersymFunctionAlgebra_powersum
from sage.combinat.super_sf.hom_el import SupersymFunctionAlgebra_hom_el
from sage.combinat.super_sf.schur import SupersymFunctionAlgebra_schur

class SuperSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The abstract class of commutative supersymmetric functions.

    A *supersymmetric function* is a symmetric function `f(\mathbb{x}|\mathbb{y})`
    defined over the variable sets `\mathbb{x}=(x_1, \ldots, x_m)` and `\mathbb{y}=(y_1, \ldots, y_n)`
    such that substituting `x_m = t` and `y_n = -t` results in `f(\mathbb{x'}|\mathbb{y'})`, where
    `\mathbb{x}=(x_1, \ldots, x_{m-1})` and `\mathbb{y}=(y_1, \ldots, y_{n-1})`.

    INPUT:

    - ``R`` -- commutative ring

    REFERENCES:

    - [BHS25]_

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: R = SuperSymmetricFunctions(QQ)
        sage: TestSuite(R).run(skip="_test_fraction_field")
    """
    def __init__(self, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: R = SuperSymmetricFunctions(QQ)
            sage: TestSuite(R).run(skip="_test_fraction_field")
        """
        if R not in Rings().Commutative():
            raise ValueError("the base ring must be a commutative ring")
        cat = GradedHopfAlgebras(R).Commutative().Cocommutative()
        if R in PrincipalIdealDomains():
            cat &= UniqueFactorizationDomains()
        Parent.__init__(self, base=R, category=cat.WithRealizations())

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the powersum basis).

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: Sym = SuperSymmetricFunctions(QQ)
            sage: Sym.a_realization()
            Supersymmetric functions over Rational Field in the powersum basis
        """
        return self.powersum()

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: R = SuperSymmetricFunctions(QQ)
            sage: R._repr_()
            'Supersymmetric functions over Rational Field'
        """
        return "Supersymmetric functions over %s" % self._base

    def powersum(self):
        r"""
        The powersum basis of supersymmetric functions.

        .. SEEALSO::

            :mod:`sage.combinat.super_sf.powersum`

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: s.powersum()
            Supersymmetric functions over Rational Field in the powersum basis
        """
        return SupersymFunctionAlgebra_powersum(self)

    p = powersum

    def homogeneous(self):
        r"""
        The homogeneous basis of supersymmetric functions.

        .. SEEALSO::

            :mod:`sage.combinat.super_sf.hom_el`

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: s.homogeneous()
            Supersymmetric functions over Rational Field in the homogeneous basis
        """
        return SupersymFunctionAlgebra_hom_el(self)

    h = homogeneous

    def elementary(self):
        r"""
        The elementary basis of supersymmetric functions.

        .. SEEALSO::

            :mod:`sage.combinat.super_sf.hom_el`

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: s.elementary()
            Supersymmetric functions over Rational Field in the elementary basis
        """
        return SupersymFunctionAlgebra_hom_el(self, basis_name='elementary')

    e = elementary

    def schur(self):
        r"""
        The Schur basis of supersymmetric functions.

        .. SEEALSO::

            :mod:`sage.combinat.super_sf.schur`

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: s.schur()
            Supersymmetric functions over Rational Field in the Schur basis
        """
        return SupersymFunctionAlgebra_schur(self)

    s = schur
