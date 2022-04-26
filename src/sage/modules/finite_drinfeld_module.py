r"""
<Short one-line summary that ends with no period>

<Paragraph description>

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Antoine Leudière (2022-04-26): initial version

"""

#*****************************************************************************
#       Copyright (C) 2022 Antoine Leudière <antoine.leudiere@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.action import Action
from sage.categories.homset import Hom
from sage.rings.morphism import RingHomomorphism_im_gens
from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_finite_field
from sage.misc.latex import latex


class FiniteDrinfeldModule(RingHomomorphism_im_gens):
    
    def __init__(self, polring, gen):
        # VERIFICATIONS
        # Check `polring` is an Fq[X]:
        # See docstrings of `PolynomialRing_dense_finite_field` and
        # `is_PolynomialRing`.
        isinstance(polring, PolynomialRing_dense_finite_field)
        # Check `gen` is an Ore polynomial:
        if not isinstance(gen, OrePolynomial):
            raise TypeError('The generator must be an Ore polynomial')
        # Now we can define those for convenience:
        FqX = polring
        Ltau = gen.parent()
        Fq = FqX.base_ring()
        L = Ltau.base_ring()
        # Check the Ore polynomial ring is an L{tau} with L a finite
        # field extension of Fq:
        _check_base_fields(Fq, L)
        if not Ltau.twisting_derivation() is None:
            raise ValueError('The Ore polynomial ring should have no ' \
                    'derivation')
        # Check the frobenius is x -> x^q:
        if Ltau.twisting_morphism().power() != Fq.degree():
            raise ValueError('The twisting morphism of the Ore polynomial ' \
                    'ring must be the Frobenius endomorphism of the base ' \
                    'field of the polynomial ring')
        # The generator is not constant:
        if gen.is_constant():
            raise ValueError('The generator must not be constant')
        # ACTUAL WORK
        super().__init__(Hom(FqX, Ltau), gen)

    ###########
    # Methods #
    ###########

    def change_ring(self, R):
        # VERIFICATIONS
        if not R.is_field() and R.is_finite():
            raise TypeError('Argument must be a finite field')
        if not self.ore_polring().base_ring().is_subring(R):
            raise ValueError('The new field must be a finite field ' \
                    'extension of the base field of the Ore polynomial ring.')
        _check_base_fields(self.polring().base_ring(), R)
        # ACTUAL WORK
        new_frobenius = R.frobenius_endomorphism(self.frobenius().power())
        new_ore_polring = OrePolynomialRing(R, new_frobenius,
                names=self.ore_polring().variable_names())
        return FiniteDrinfeldModule(self.polring(), new_ore_polring(self.gen()))

    def rank(self):
        return self.gen().degree()

    ##########################
    # Special Sage functions #
    ##########################

    def _get_action_(self, extension):
        return FiniteDrinfeldModuleAction(self, extension)

    def _latex_(self):
        return f'\\text{{Finite{{ }}{latex(self.polring())}-Drinfeld{{ }}' \
                f'module{{ }}defined{{ }}over{{ }}}}' \
                f'{latex(self.ore_polring().base_ring())}\\text{{{{ }}' \
                f'by{{ }}}}\n' \
                f'\\begin{{align}}\n' \
                f'  {latex(self.polring())}\n' \
                f'  &\\to {latex(self.ore_polring())} \\\\\n' \
                f'  {latex(self.polring().gen())}\n' \
                f'  &\\mapsto {latex(self.gen())}\n' \
                f'\\end{{align}}'

    def _repr_(self):
        return f'Finite Drinfeld module from {self.polring()} over ' \
                f'{self.ore_polring().base_ring()} defined by {self.gen()}.'

    ###########
    # Getters #
    ###########

    def frobenius(self):
        return self.ore_polring().twisting_morphism()

    def gen(self):
        [gen] = self.im_gens()
        return gen

    def ore_polring(self):
        return self.codomain()

    def polring(self):
        return self.domain()

class FiniteDrinfeldModuleAction(Action):

    def __init__(self, finite_drinfeld_module, extension):
        # VERIFICATIONS
        if not isinstance(finite_drinfeld_module, FiniteDrinfeldModule):
            raise TypeError('First argument must be a FiniteDrinfeldModule')
        if not (extension.is_field() and extension.is_finite() and
                finite_drinfeld_module.polring().base_ring().is_subring(extension)):
            raise ValueError('The extension must be a finite field ' \
                    'extension of the base field of the Ore polynomial ring')
        # WORK
        self.__finite_drinfeld_module = finite_drinfeld_module
        super().__init__(finite_drinfeld_module.polring(), extension)

    ###########
    # Methods #
    ###########

    def extension(self):
        return self.codomain()

    def finite_drinfeld_module(self):
        return self.__finite_drinfeld_module

    ##########################
    # Special Sage functions #
    ##########################

    def _latex_(self):
        phi = self.finite_drinfeld_module()
        return f'\\text{{Drinfeld{{ }}module{{ }}action{{ }}' \
                f'on{{ }}}}{latex(self.extension())}\\text{{{{ }}' \
                f'induced{{ }}by{{ }}}}{latex(phi)}'

    def _repr_(self):
        return f'Action on {self.domain()} induced by the ' \
                f'{self.finite_drinfeld_module()}'

    def _act_(self, g, x):
        return self.finite_drinfeld_module().change_ring(self.extension())(g)(x)


def _check_base_fields(Fq, L):
    if not (L.is_field() and L.is_finite() and Fq.is_subring(L)):
        raise ValueError(f'The base field of the Ore polynomial ring must ' \
                'be a finite field extension of the base field of the ' \
                'polynomial ring')
