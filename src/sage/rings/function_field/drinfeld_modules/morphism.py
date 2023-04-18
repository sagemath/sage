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
    This class represents Drinfeld `\mathbb{F}_q[T]`-module morphisms.

    Let `\phi` and `\psi` be two Drinfeld `\mathbb{F}_q[T]`-modules over
    a field `K`. A *morphism of Drinfeld modules* `\phi \to \psi` is an
    Ore polynomial `f \in K\{\tau\}` such that `f \phi_a = \psi_a f` for
    every `a \in \mathbb{F}_q[T]`. In our case, this is equivalent to `f
    \phi_T = \psi_T f`. An *isogeny* is a nonzero morphism.

    To create a morphism object, the user should never explicitly
    instantiate :class:`DrinfeldModuleMorphism`, but rather call the
    parent homset with the defining Ore polynomial::

        sage: Fq = GF(4)
        sage: A.<T> = Fq[]
        sage: K.<z> = Fq.extension(3)
        sage: phi = DrinfeldModule(A, [z, z^2 + z, z^2 + z])
        sage: t = phi.ore_polring().gen()
        sage: ore_pol = t + z^5 + z^3 + z + 1
        sage: psi = phi.velu(ore_pol)
        sage: morphism = Hom(phi, psi)(ore_pol)
        sage: morphism
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> (z^2 + z)*t^2 + (z^2 + z)*t + z
          To:   Drinfeld module defined by T |--> (z^5 + z^2 + z + 1)*t^2 + (z^4 + z + 1)*t + z
          Defn: t + z^5 + z^3 + z + 1


    The given Ore polynomial must indeed define a morphism::

        sage: morphism = Hom(phi, psi)(1)
        Traceback (most recent call last):
        ...
        ValueError: Ore polynomial does not define a morphism

    One can get basic data on the morphism::

        sage: morphism.domain()
        Drinfeld module defined by T |--> (z^2 + z)*t^2 + (z^2 + z)*t + z
        sage: morphism.domain() is phi
        True

        sage: morphism.codomain()
        Drinfeld module defined by T |--> (z^5 + z^2 + z + 1)*t^2 + (z^4 + z + 1)*t + z
        sage: morphism.codomain() is psi
        True

    ::

        sage: morphism.ore_polynomial()
        t + z^5 + z^3 + z + 1
        sage: morphism.ore_polynomial() is ore_pol
        True

    One can check various properties::

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

    For the sake of completeness, we explain how the user can directly
    instantiate the class, even though this should never be explicitly
    done::

        sage: from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
        sage: DrinfeldModuleMorphism(Hom(phi, psi), ore_pol)
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> (z^2 + z)*t^2 + (z^2 + z)*t + z
          To:   Drinfeld module defined by T |--> (z^5 + z^2 + z + 1)*t^2 + (z^4 + z + 1)*t + z
          Defn: t + z^5 + z^3 + z + 1
        sage: DrinfeldModuleMorphism(Hom(phi, psi), ore_pol) is morphism
        True
    """

    @staticmethod
    def __classcall_private__(cls, parent, x):
        """
        Create the morphism.

        INPUT:

        - ``cls`` -- the class ``DrinfeldModuleMorphism``

        - ``parent`` -- the Drinfeld module homset

        - ``x`` -- the Ore polynomial defining the morphism or a
          DrinfeldModuleMorphism

        TESTS::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism is Hom(phi, psi)(morphism)
            True

        ::

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
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
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

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: latex(morphism)
            t + z_{6}^{5} + z_{6}^{2} + 1
        """
        return f'{latex(self._ore_polynomial)}'

    def _repr_(self):
        r"""
        Return a string representation of the morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> t^2 + t + z6
              To:   Drinfeld module defined by T |--> t^2 + (z6^4 + z6^2 + 1)*t + z6
              Defn: t + z6^5 + z6^2 + 1
        """
        if self.is_identity():
            return f'Identity morphism of {self._domain}'
        elif self.is_endomorphism():
            return f'Endomorphism of {self._domain}\n' \
                   f'  Defn: {self._ore_polynomial}'
        else:
            return f'Drinfeld Module morphism:\n' \
                   f'  From: {self._domain}\n'  \
                   f'  To:   {self._codomain}\n' \
                   f'  Defn: {self._ore_polynomial}'

    def __hash__(self):
        r"""
        Return a hash of ``self``.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: hash(morphism)  # random
            -4214883752078138009
        """
        return hash((self._domain, self._codomain, self._ore_polynomial))

    def is_zero(self):
        r"""
        Return ``True`` whether the morphism is the zero morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_zero()
            False

        ::

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_zero()
            True
        """
        return self._ore_polynomial.is_zero()

    def is_identity(self):
        r"""
        Return ``True`` whether the morphism is the identity morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: morphism = End(phi)(1)
            sage: morphism.is_identity()
            True

        ::

            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_identity()
            False
        """
        return self._ore_polynomial == 1

    def is_isogeny(self):
        r"""
        Return ``True`` whether the morphism is an isogeny.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_isogeny()
            True

        ::

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_isogeny()
            False

        ::

            sage: identity_morphism = End(phi)(1)
            sage: identity_morphism.is_isogeny()
            True

        ::

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
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_isomorphism()
            False

        ::

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_isomorphism()
            False

        ::

            sage: identity_morphism = End(phi)(1)
            sage: identity_morphism.is_isomorphism()
            True

        ::

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
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: ore_pol = morphism.ore_polynomial()
            sage: ore_pol
            t + z6^5 + z6^2 + 1

        ::

            sage: ore_pol * phi(T) == psi(T) * ore_pol
            True
        """
        return self._ore_polynomial
