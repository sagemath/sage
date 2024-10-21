# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.pari
"""
Submodules of spaces of modular forms

EXAMPLES::

    sage: M = ModularForms(Gamma1(13),2); M
    Modular Forms space of dimension 13 for
     Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: M.eisenstein_subspace()
    Eisenstein subspace of dimension 11 of Modular Forms space of dimension 13 for
     Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: M == loads(dumps(M))
    True
    sage: M.cuspidal_subspace()
    Cuspidal subspace of dimension 2 of Modular Forms space of dimension 13 for
     Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from .space import ModularFormsSpace

import sage.modular.hecke.submodule


class ModularFormsSubmodule(ModularFormsSpace,
                            sage.modular.hecke.submodule.HeckeSubmodule):
    """
    A submodule of an ambient space of modular forms.
    """
    def __init__(self, ambient_module, submodule, dual=None, check=False):
        """
        INPUT:

        - ``ambient_module`` -- ModularFormsSpace
        - ``submodule`` -- a submodule of the ambient space
        - ``dual_module`` -- (default: ``None``) ignored
        - ``check`` -- boolean (default: ``False``); whether to check that the
          submodule is Hecke equivariant

        EXAMPLES::

          sage: M = ModularForms(Gamma1(13),2); M
          Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
          sage: M.eisenstein_subspace()
          Eisenstein subspace of dimension 11 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
        """
        A = ambient_module
        sage.modular.hecke.submodule.HeckeSubmodule.__init__(self, A, submodule, check=check)
        ModularFormsSpace.__init__(self, A.group(), A.weight(),
                                         A.character(), A.base_ring())

    def _repr_(self):
        """
        EXAMPLES::

          sage: ModularForms(Gamma1(13),2).eisenstein_subspace()._repr_()
          'Eisenstein subspace of dimension 11 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field'
        """
        return "Modular Forms subspace of dimension %s of %s" % (self.dimension(), self.ambient_module())

    def _compute_coefficients(self, element, X):
        """
        Compute all coefficients of the modular form element in ``self`` for
        indices in X.

        TODO: Implement this function.

        EXAMPLES::

            sage: M = ModularForms(6,4).cuspidal_subspace()
            sage: M._compute_coefficients( M.basis()[0], range(1,100) )
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _compute_q_expansion_basis(self, prec):
        """
        Compute ``q_expansions`` to precision ``prec`` for each element in
        ``self.basis()``.

        EXAMPLES::

            sage: M = ModularForms(Gamma1(13),2); M
            Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
            sage: S = M.eisenstein_subspace(); S
            Eisenstein subspace of dimension 11 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
            sage: S._compute_q_expansion_basis(5)
            [1 + O(q^5),
             q + O(q^5),
             q^2 + O(q^5),
             q^3 + O(q^5),
             q^4 + O(q^5),
             O(q^5),
             O(q^5),
             O(q^5),
             O(q^5),
             O(q^5),
             O(q^5)]
        """
        A = self.ambient_module()
        return [A._q_expansion(element=f.element(), prec=prec)
                for f in self.basis()]


# TODO
class ModularFormsSubmoduleWithBasis(ModularFormsSubmodule):
    pass
