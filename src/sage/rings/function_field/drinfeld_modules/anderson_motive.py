# *****************************************************************************
#        Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

import operator

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.latex import latex
from sage.misc.functional import log

from sage.categories.homset import Homset
from sage.categories.anderson_motives import AndersonMotives
from sage.structure.factorization import Factorization

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField_1poly_field

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix, block_diagonal_matrix

from sage.modules.ore_module import OreModule, OreSubmodule, OreQuotientModule
from sage.modules.ore_module_element import OreModuleElement
from sage.modules.ore_module import OreAction
from sage.modules.ore_module import normalize_names

from sage.rings.function_field.drinfeld_modules.anderson_motive_morphism import DrinfeldToAnderson, AndersonToDrinfeld



class AndersonMotiveElement(OreModuleElement):
    def image(self, integral=None):
        if integral is None:
            integral = self.parent().is_effective()
        return super().image(integral=integral)


class AndersonMotive_general(OreModule):
    Element = AndersonMotiveElement

    @staticmethod
    def __classcall_private__(self, category, tau, twist=0, names=None, normalize=True):
        K = category.base()
        AK = category.base_combined()

        # We normalize the inputs
        twist = ZZ(twist)
        tau = tau.change_ring(AK)
        if normalize:
            divisor = category.divisor()
            exponent = Infinity
            for entry in tau.list():
                if not entry:
                    continue
                e = 0
                while entry.degree() > 0 and e < exponent:
                    entry, R = entry.quo_rem(divisor)
                    if R:
                        break
                    e += 1
                exponent = e
                if exponent == 0:
                    break
            if exponent is not Infinity and exponent > 0:
                denom = divisor ** exponent
                tau = tau.parent()([entry // denom for entry in tau.list()])
                twist -= exponent

        names = normalize_names(names, tau.nrows())
        denominator = Factorization([(category.divisor(), twist)])
        ore = category._ore_polring

        #if (isinstance(K, FractionField_1poly_field)
        #    and category.constant_coefficient() == K.gen()):
        #    from sage.rings.function_field.drinfeld_modules.anderson_motive_rational import AndersonMotive_rational
        #    cls = AndersonMotive_rational
        #else:
        cls = AndersonMotive_general

        return cls.__classcall__(cls, tau, ore, denominator, names, category)

    def __init__(self, mat, ore, denominator, names, category) -> None:
        OreModule.__init__(self, mat, ore, denominator, names, category)
        self._initialize_attributes()

    def _initialize_attributes(self):
        category = self._category
        self._A = A = category.function_ring()
        self._t_name = A.variable_name()
        self._Fq = Fq = A.base_ring()
        self._q = Fq.cardinality()
        self._deg = ZZ(log(self._q, Fq.characteristic()))
        self._K = self._base = K = category.base()
        self._theta = category.constant_coefficient()
        self._AK = base = category.base_combined()
        self._t = base.gen()
        self._tau = self._pseudohom.matrix()
        if self._denominator:
            self._twist = self._denominator[0][1]
        else:
            self._twist = 0
        self._submodule_class = AndersonSubMotive
        self._quotientModule_class = AndersonQuotientMotive

    @lazy_attribute
    def _dettau(self):
        det = self._tau.det()
        return det.leading_coefficient(), det.degree()

    def _repr_(self):
        s = "Anderson motive "
        if self._names is None:
            s += "of rank %s " % self.rank()
        else:
            s += "<" + ", ".join(self._names) + "> "
        s += "over %s" % self._AK
        return s

    def _latex_(self):
        if self._names is None:
            s = "\\texttt{Anderson motive of rank } %s" % self.rank()
            s += "\\texttt{ over } %s" % latex(self._AK)
        else:
            s = "\\left<" + ", ".join(self._latex_names) + "\\right>"
            s += "_{%s}" % latex(self._AK)
        return s

    def twist(self, n, names):
        return AndersonMotive_general(self._category, self._tau, self._twist + ZZ(n),
                                      names, normalize=False)

    def _Hom_(self, codomain, category):
        from sage.rings.function_field.drinfeld_modules.anderson_motive_morphism import AndersonMotive_homspace
        return AndersonMotive_homspace(self, codomain)

    def hodge_pink_weights(self):
        S = self._tau.smith_form(transformation=False)
        return [-self._twist + S[i,i].degree() for i in range(self.rank())]

    def is_effective(self):
        return self._twist <= 0

    def ore_variable(self):
        return self._category._ore_polring.gen()

    def ore_polring(self, names=None, action=True):
        if names is None:
            names = self._category._ore_variable_name
        S = self._ore_category.ore_ring(names)
        if action:
            self._unset_coercions_used()
            self.register_action(OreAction(S, self, True, operator.mul))
        return S


class AndersonMotive_drinfeld(AndersonMotive_general):
    def __init__(self, phi, names):
        category = AndersonMotives(phi.category())
        AK = category.base_combined()
        r = phi.rank()
        tau = matrix(AK, r)
        P = phi.gen()
        tau[r-1, 0] = (AK.gen() - P[0]) / P[r]
        for i in range(1, r):
            tau[i-1, i] = 1
            tau[r-1, i] = -P[i]/P[r]
        AndersonMotive_general.__init__(self, tau, category._ore_polring, None, names, category)
        Ktau = phi.ore_polring()
        self.register_coercion(DrinfeldToAnderson(Homset(Ktau, self), phi))
        try:
            Ktau.register_conversion(AndersonToDrinfeld(Homset(self, Ktau), phi))
        except AssertionError:
            pass
        self._drinfeld_module = phi

    def drinfeld_module(self):
        return self._drinfeld_module


class AndersonSubMotive(AndersonMotive_general, OreSubmodule):
    def __init__(self, ambient, submodule, names):
        OreSubmodule.__init__(self, ambient, submodule, names)
        self._initialize_attributes()


class AndersonQuotientMotive(AndersonMotive_general, OreQuotientModule):
    def __init__(self, cover, submodule, names):
        OreQuotientModule.__init__(self, cover, submodule, names)
        self._initialize_attributes()
