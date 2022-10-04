r"""
Drinfeld module morphisms

This module provides the class
:class:`sage.rings.function_fields.drinfeld_module.morphism.DrinfeldModuleMorphism`.

AUTHORS:
- Antoine Leudière (2022-04)
"""

# *****************************************************************************
#        Copyright (C) 2022 Antoine Leudière <antoine.leudiere@inria.fr>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.latex import latex
from sage.categories.morphism import Morphism
from sage.structure.unique_representation import UniqueRepresentation


class DrinfeldModuleMorphism(Morphism, UniqueRepresentation,
                             metaclass=InheritComparisonClasscallMetaclass):
    r"""
    This class represents a Drinfeld module morphism.

    Let `\phi, \psi` be two Drinfeld modules with base `\gamma: \mathbb{F}_q[X]
    \to K`. A *morphism of Drinfeld modules `\phi \to \psi`* is an Ore
    polynomial `f \in K\{\tau\}` such that `f \phi_a = \psi_a f` for every `a
    \in \mathbb{F}_q[X]`. In our case, this is equivalent to `f \phi_X = \psi_X
    f`. An *isogeny* is a non-zero morphism.

    To create a morphism object, the user should never explicitly
    instantiate :class:`DrinfeldModuleMorphism`, but rather call the parent
    homset with the defining Ore polynomial::

        sage: Fq = GF(25)
        sage: FqX.<X> = Fq[]
        sage: K.<z12> = Fq.extension(6)
        sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
        sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
        sage: t = phi.ore_polring().gen()
        sage: ore_pol = t + 2*z12^11 + 4*z12^9 + 2*z12^8 + 2*z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + 4*z12^2 + 4*z12 + 4
        sage: psi = phi.velu(ore_pol)
        sage: morphism = Hom(phi, psi)(ore_pol)
        sage: morphism
        Drinfeld Module morphism:
              From (gen): z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
              To (gen):   (z12^11 + 3*z12^10 + z12^9 + z12^7 + z12^5 + 4*z12^4 + 4*z12^3 + z12^2 + 1)*t^2 + (2*z12^11 + 4*z12^10 + 2*z12^8 + z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + z12^2 + z12 + 4)*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
              Defn:       t + 2*z12^11 + 4*z12^9 + 2*z12^8 + 2*z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + 4*z12^2 + 4*z12 + 4

    The input Ore polynomial must indeed define a morphism::

        sage: morphism = Hom(phi, psi)(1)
        Traceback (most recent call last):
        ...
        ValueError: Ore polynomial does not define a morphism

    We can get basic data on the morphism::

        sage: morphism.domain()
        Drinfeld module defined by X |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 5^2
          To:   Finite Field in z12 of size 5^12
          Defn: X |--> 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
        sage: morphism.domain() is phi
        True

        sage: morphism.codomain()
        Drinfeld module defined by X |--> (z12^11 + 3*z12^10 + z12^9 + z12^7 + z12^5 + 4*z12^4 + 4*z12^3 + z12^2 + 1)*t^2 + (2*z12^11 + 4*z12^10 + 2*z12^8 + z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + z12^2 + z12 + 4)*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 5^2
          To:   Finite Field in z12 of size 5^12
          Defn: X |--> 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
        sage: morphism.codomain() is psi
        True

        sage: morphism.ore_polynomial()
        t + 2*z12^11 + 4*z12^9 + 2*z12^8 + 2*z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + 4*z12^2 + 4*z12 + 4
        sage: morphism.ore_polynomial() is ore_pol
        True

    We can check various properties::

        sage: morphism.is_zero()
        False
        sage: morphism.is_isogeny()
        True
        sage: morphism.is_endomorphism()
        False
        sage: morphism.is_isomorphism()
        False

    TESTS::

        sage: morphism.parent() == Hom(phi, psi)
        True
        sage: phi.frobenius_endomorphism().parent() == End(phi)
        True
        sage: End(phi)(0).parent() == End(phi)
        True

    .. NOTE::

        For the sake of completeness, we explain how the user can
        directly instantiate the class, even though this should never be
        explicitly done::

            sage: from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
            sage: DrinfeldModuleMorphism(Hom(phi, psi), ore_pol)
            Drinfeld Module morphism:
              From (gen): z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
              To (gen):   (z12^11 + 3*z12^10 + z12^9 + z12^7 + z12^5 + 4*z12^4 + 4*z12^3 + z12^2 + 1)*t^2 + (2*z12^11 + 4*z12^10 + 2*z12^8 + z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + z12^2 + z12 + 4)*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
              Defn:       t + 2*z12^11 + 4*z12^9 + 2*z12^8 + 2*z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + 4*z12^2 + 4*z12 + 4
            sage: DrinfeldModuleMorphism(Hom(phi, psi), ore_pol) is morphism
            True
    """

    @staticmethod
    def __classcall_private__(cls, parent, x):
        """
        Create the morphism.

        INPUT:

        - ``cls`` -- DrinfeldModuleMorphism
        - ``parent`` -- the Drinfeld module homset
        - ``x`` -- the Ore polynomial defining the morphism or a
          DrinfeldModuleMorphism

        OUTPUT: the morphism object

        TESTS::

            sage: Fq = GF(2)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(FqX, [z6, 1, 1])
            sage: psi = DrinfeldModule(FqX, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism is Hom(phi, psi)(morphism)
            True

            sage: from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
            sage: morphism = DrinfeldModuleMorphism(Sets(), t + 1)
            Traceback (most recent call last):
            ...
            TypeError: parent should be a DrinfeldModuleHomset
        """
        from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
        if not isinstance(parent, DrinfeldModuleHomset):
            raise TypeError('parent should be a DrinfeldModuleHomset')
        domain = parent.domain()
        codomain = parent.codomain()
        # NOTE: it used to be x.parent() is parent, but this was False.
        # DrinfeldModuleHomset inherits Homset, which does NOT inherit
        # UniqueRepresentation
        if x.parent() == parent:  # x is a DrinfeldModuleMorphism
            ore_pol = x.ore_polynomial()
        else:  # x is an Ore polynomial
            ore_pol = domain.ore_polring()(x)
        if ore_pol * domain.gen() != codomain.gen() * ore_pol:
            raise ValueError('Ore polynomial does not define a morphism')
        return cls.__classcall__(cls, parent, ore_pol)

    def __init__(self, parent, ore_pol):
        r"""
        Initialize ``self``.

        INPUT:

        - ``parent`` -- the Drinfeld module homset
        - ``ore_pol`` -- the Ore polynomial that defines the morphism

        TESTS::

            sage: Fq = GF(2)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(FqX, [z6, 1, 1])
            sage: psi = DrinfeldModule(FqX, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism._domain is phi
            True
            sage: morphism._codomain is psi
            True
            sage: morphism._ore_polynomial == t + z6^5 + z6^2 + 1
            True
        """

        super().__init__(parent)
        self._domain = parent.domain()
        self._codomain = parent.codomain()
        self._ore_polynomial = ore_pol

    def _latex_(self):
        r"""
        Return a LaTeX representation of the morphism.

        OUTPUT: a string

        EXAMPLES::

            sage: Fq = GF(2)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(FqX, [z6, 1, 1])
            sage: psi = DrinfeldModule(FqX, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: latex(morphism)
            \begin{array}{l}
            \text{Drinfeld{ }module{ }morphism:}\\
            \text{{ }{ }From (gen):{ }}t^{2} + t + z_{6}\\
            \text{{ }{ }To (gen):{ }}{ }{ }t^{2} + \left(z_{6}^{4} + z_{6}^{2} + 1\right) t + z_{6}\\
            \text{{ }{ }Defn:{ }}t + z_{6}^{5} + z_{6}^{2} + 1
            \end{array}
        """
        return f'\\begin{{array}}{{l}}\n' \
               f'\\text{{Drinfeld{{ }}module{{ }}morphism:}}\\\\\n' \
               f'\\text{{{{ }}{{ }}From (gen):{{ }}}}'\
               f'{latex(self._domain.gen())}\\\\\n' \
               f'\\text{{{{ }}{{ }}To (gen):{{ }}}}{{ }}{{ }}' \
               f'{latex(self._codomain.gen())}\\\\\n' \
               f'\\text{{{{ }}{{ }}Defn:{{ }}}}' \
               f'{latex(self._ore_polynomial)}\n' \
               f'\\end{{array}}'

    def _repr_(self):
        r"""
        Return a string representation of the morphism.

        OUTPUT: a string

        EXAMPLES::

            sage: Fq = GF(2)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(FqX, [z6, 1, 1])
            sage: psi = DrinfeldModule(FqX, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism
            Drinfeld Module morphism:
              From (gen): t^2 + t + z6
              To (gen):   t^2 + (z6^4 + z6^2 + 1)*t + z6
              Defn:       t + z6^5 + z6^2 + 1
        """
        return f'Drinfeld Module morphism:\n' \
               f'  From (gen): {self._domain.gen()}\n'  \
               f'  To (gen):   {self._codomain.gen()}\n' \
               f'  Defn:       {self._ore_polynomial}'

    def is_zero(self):
        r"""
        Return ``True`` whether the morphism is the zero morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(FqX, [z6, 1, 1])
            sage: psi = DrinfeldModule(FqX, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_zero()
            False

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_zero()
            True
        """
        return self._ore_polynomial.is_zero()

    def is_isogeny(self):
        r"""
        Return ``True`` whether the morphism is an isogeny.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(FqX, [z6, 1, 1])
            sage: psi = DrinfeldModule(FqX, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_isogeny()
            True

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_isogeny()
            False

            sage: identity_morphism = End(phi)(1)
            sage: identity_morphism.is_isogeny()
            True

            sage: frobenius_endomorphism = phi.frobenius_endomorphism()
            sage: frobenius_endomorphism.is_isogeny()
            True
        """
        return not self.is_zero()

    def is_isomorphism(self):
        r"""
        Return ``True`` whether the morphism is an isomorphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(FqX, [z6, 1, 1])
            sage: psi = DrinfeldModule(FqX, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_isomorphism()
            False

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_isomorphism()
            False

            sage: identity_morphism = End(phi)(1)
            sage: identity_morphism.is_isomorphism()
            True

            sage: frobenius_endomorphism = phi.frobenius_endomorphism()
            sage: frobenius_endomorphism.is_isomorphism()
            False
        """
        return self.is_isogeny() and self._ore_polynomial.degree() == 0

    def ore_polynomial(self):
        r"""
        Return the Ore polynomial that defines the morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(FqX, [z6, 1, 1])
            sage: psi = DrinfeldModule(FqX, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: ore_pol = morphism.ore_polynomial()
            sage: ore_pol
            t + z6^5 + z6^2 + 1

            sage: ore_pol * phi(X) == psi(X) * ore_pol
            True
        """
        return self._ore_polynomial
