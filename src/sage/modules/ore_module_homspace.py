r"""
Space of morphisms between Ore modules

AUTHOR:

- Xavier Caruso (2024-10)
"""

# ***************************************************************************
#    Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.homset import HomsetWithBase
from sage.matrix.matrix_space import MatrixSpace

from sage.modules.ore_module import OreModule
from sage.modules.ore_module_morphism import OreModuleMorphism

class OreModule_homspace(UniqueRepresentation, HomsetWithBase):
    r"""
    Class for hom spaces between Ore modules.
    """
    Element = OreModuleMorphism

    def __init__(self, domain, codomain, category=None):
        r"""
        Initialize this homspace.

        INPUT:

        - ``domain`` -- a Ore module

        - ``codomain`` -- a Ore module

        - ``category`` (default: ``None``) -- the category in which
          the morphisms are

        TESTS::

            sage: K.<z> = GF(7^2)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X + z)
            sage: N = S.quotient_module(X^2 + z)
            sage: Hom(M, N)  # indirect doctest
            Set of Morphisms
            from Ore module of rank 1 over Finite Field in z of size 7^2 twisted by z |--> z^7
            to Ore module of rank 2 over Finite Field in z of size 7^2 twisted by z |--> z^7
            in Category of enumerated finite dimensional Ore modules with basis over Finite Field in z of size 7^2 twisted by z |--> z^7
            sage: End(M)     # indirect doctest
            Set of Morphisms
            from Ore module of rank 1 over Finite Field in z of size 7^2 twisted by z |--> z^7
            to Ore module of rank 1 over Finite Field in z of size 7^2 twisted by z |--> z^7
            in Category of enumerated finite dimensional Ore modules with basis over Finite Field in z of size 7^2 twisted by z |--> z^7

        ::

            sage: V = M.module()
            sage: Hom(M, V)
            Traceback (most recent call last):
            ...
            ValueError: codomain must be a Ore module
        """
        if not isinstance(domain, OreModule):
            raise ValueError("domain must be a Ore module")
        if not isinstance(codomain, OreModule):
            raise ValueError("codomain must be a Ore module")
        if domain.ore_ring(action=False) is not codomain.ore_ring(action=False):
            raise ValueError("domain and codomain must be defined over the same ring with same twisting maps")
        super().__init__(domain, codomain, category)
        base = domain.base_ring()
        self._matrix_space = MatrixSpace(base, domain.dimension(), codomain.dimension())

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return a morphism in this homspace.

        TESTS::

            sage: K.<z> = GF(7^2)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X + z)
            sage: H = End(M)
            sage: H(2)
            Ore module endomorphism of Ore module of rank 1 over Finite Field in z of size 7^2 twisted by z |--> z^7
        """
        return self.element_class(self, *args, **kwds)

    def matrix_space(self):
        r"""
        Return the matrix space used to represent the
        morphisms in this homspace.

        EXAMPLES::

            sage: K.<z> = GF(7^2)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^3 + z*X + 1)
            sage: End(M).matrix_space()
            Full MatrixSpace of 3 by 3 dense matrices over Finite Field in z of size 7^2

        ::

            sage: N = S.quotient_module(X^2 + z)
            sage: Hom(M, N).matrix_space()
            Full MatrixSpace of 3 by 2 dense matrices over Finite Field in z of size 7^2
        """
        return self._matrix_space
