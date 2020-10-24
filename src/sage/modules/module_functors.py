"""
Module Functors

AUTHORS:

- Travis Scrimshaw (2017-10): Initial implementation of
  :class:`QuotientModuleFunctor`
"""

#*****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
 
##############################################################
# Construction functor for quotient modules
##############################################################

from sage.categories.pushout import ConstructionFunctor
from sage.categories.modules import Modules

class QuotientModuleFunctor(ConstructionFunctor):
    r"""
    Construct the quotient of a module by a submodule.

    INPUT:

    - ``relations`` -- a module

    .. NOTE::

        This construction functor keeps track of the basis of defining
        ``relations``. It can only be applied to free modules into which
        this basis coerces.

    EXAMPLES::

        sage: A = (1/2)*ZZ^2
        sage: B = 2*ZZ^2
        sage: Q = A / B
        sage: F = Q.construction()[0]
        sage: F
        QuotientModuleFunctor
        sage: F(A) == Q
        True
        
    The modules are constructed from the cover not the ambient module::
    
        sage: F(B.ambient_module()) == Q
        False

    We can construct quotients from different modules::

        sage: F((1/2)*ZZ^2)
        Finitely generated module V/W over Integer Ring with invariants (4, 4)
        sage: F(ZZ^2)
        Finitely generated module V/W over Integer Ring with invariants (2, 2)
        sage: F(2*ZZ^2)
        Finitely generated module V/W over Integer Ring with invariants ()

    This functor is used for constructing pushouts::

        sage: A = ZZ^3
        sage: x,y,z = A.basis()
        sage: A1 = A.submodule([x])
        sage: A2 = A.submodule([y, 2*x])
        sage: B1 = A.submodule([])
        sage: B2 = A.submodule([2*x])
        sage: Q1 = A1 / B1
        sage: Q2 = A2 / B2
        sage: q3 = Q1.an_element() + Q2.an_element()
    """
    rank = 5 # ranking of functor, not rank of module

    def __init__(self, relations):
        """
        Initialization of ``self``.

        TESTS::

            sage: from sage.modules.module_functors import QuotientModuleFunctor
            sage: B = (2/3)*ZZ^2
            sage: F = QuotientModuleFunctor(B)
            sage: TestSuite(F).run()
        """
        R = relations.category().base_ring()
        ConstructionFunctor.__init__(self, Modules(R), Modules(R))
        self._relations = relations

    def relations(self):
        """
        Return the defining relations of ``self``.

        EXAMPLES::

            sage: A = (ZZ**2) / span([[4,0],[0,3]], ZZ)
            sage: A.construction()[0].relations()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [4 0]
            [0 3]
        """
        return self._relations

    def _apply_functor(self, ambient):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: A = ZZ^3
            sage: B = 2 * A
            sage: C = 4 * A
            sage: D = B / C
            sage: F = D.construction()[0]
            sage: D == F(D.construction()[1])
            True
        """
        return ambient.quotient(self._relations)

    def __eq__(self, other):
        """
        The quotient functor ``self`` is equal to ``other`` if
        it is a :class:`QuotientModuleFunctor` and the relations
        subspace are equal.

        EXAMPLES::

            sage: F1 = ((ZZ^3) / (4*ZZ^3)).construction()[0]
            sage: F2 = ((2*ZZ^3) / (4*ZZ^3)).construction()[0]
            sage: F1 == F2
            True
            sage: F3 = ((ZZ^3) / (8*ZZ^3)).construction()[0]
            sage: F1 == F3
            False
        """
        if not isinstance(other, QuotientModuleFunctor):
            return False
        return self._relations == other._relations

    def __ne__(self, other):
        r"""
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: F1 = ((ZZ^3) / (4*ZZ^3)).construction()[0]
            sage: F2 = ((2*ZZ^3) / (4*ZZ^3)).construction()[0]
            sage: F1 != F2
            False
            sage: F3 = ((ZZ^3) / (8*ZZ^3)).construction()[0]
            sage: F1 != F3
            True
        """
        return not (self == other)

    def merge(self, other):
        r"""
        Merge the construction functors ``self`` and ``other``.

        EXAMPLES::

            sage: A = ZZ^3
            sage: x,y,z = A.basis()
            sage: A1 = A.submodule([x])
            sage: A2 = A.submodule([y, 2*x])
            sage: B1 = A.submodule([])
            sage: B2 = A.submodule([2*x])
            sage: Q1 = A1 / B1
            sage: Q2 = A2 / B2
            sage: F1 = Q1.construction()[0]
            sage: F2 = Q2.construction()[0]
            sage: F3 = F1.merge(F2)
            sage: q3 = Q1.an_element() + Q2.an_element()
            sage: q3.parent() == F3(A1 + A2)
            True

            sage: G = A1.construction()[0]; G
            SubspaceFunctor
            sage: F1.merge(G)
            sage: F2.merge(G)
        """
        if isinstance(other, QuotientModuleFunctor):
            return QuotientModuleFunctor(self._relations + other._relations)

