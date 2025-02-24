# sage.doctest: needs sage.rings.finite_rings
r"""
Set of morphisms between two Drinfeld modules

This module provides the class
:class:`sage.rings.function_field.drinfeld_module.homset.DrinfeldModuleHomset`.

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

import operator

from sage.categories.drinfeld_modules import DrinfeldModules
from sage.categories.homset import Homset
from sage.categories.action import Action
from sage.misc.latex import latex
from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
from sage.structure.parent import Parent


class DrinfeldModuleMorphismAction(Action):
    r"""
    Action of the function ring on the homset of a Drinfeld module.

    EXAMPLES::

        sage: Fq = GF(5)
        sage: A.<T> = Fq[]
        sage: K.<z> = Fq.extension(3)
        sage: phi = DrinfeldModule(A, [z, 1, z])
        sage: psi = DrinfeldModule(A, [z, z^2 + 4*z + 3, 2*z^2 + 4*z + 4])
        sage: H = Hom(phi, psi)
        sage: t = phi.ore_variable()
        sage: f = H(t + 2)

    Left action::

        sage: (T + 1) * f
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> z*t^2 + t + z
          To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^2 + (z^2 + 4*z + 3)*t + z
          Defn: (2*z^2 + 4*z + 4)*t^3 + (2*z + 1)*t^2 + (2*z^2 + 4*z + 2)*t + 2*z + 2

    Right action currently does not work (it is a known bug, due to an
    incompatibility between multiplication of morphisms and the coercion
    system)::

        sage: f * (T + 1)
        Traceback (most recent call last):
        ...
        TypeError: right (=T + 1) must be a map to multiply it by Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> z*t^2 + t + z
          To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^2 + (z^2 + 4*z + 3)*t + z
          Defn: t + 2
    """
    def __init__(self, A, H, is_left, op):
        r"""
        Initialize this action.

        INPUT:

        - ``A`` -- the function ring of the underlying Drinfeld module

        - ``H`` -- a homset between Drinfeld modules

        - ``is_left`` -- boolean

        - ``op`` -- an operator

        TESTS::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)

            sage: from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleMorphismAction
            sage: left_action = DrinfeldModuleMorphismAction(A, H, True, operator.mul)
            sage: TestSuite(left_action).run(skip='_test_pickling')

            sage: right_action = DrinfeldModuleMorphismAction(A, H, False, operator.mul)
            sage: TestSuite(right_action).run(skip='_test_pickling')
        """
        Action.__init__(self, A, H, is_left, op)
        if is_left:
            self._phi = H.codomain()
        else:
            self._phi = H.domain()

    def _act_(self, a, f):
        r"""
        Return the action of `a` on `f`.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: t = phi.ore_variable()
            sage: f = phi.hom(t + 1)
            sage: T*f  # indirect doctest
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> z*t^3 + t^2 + z
              To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^3 + (3*z^2 + 2*z + 2)*t^2 + (2*z^2 + 3*z + 4)*t + z
              Defn: (2*z^2 + 4*z + 4)*t^4 + (z + 1)*t^3 + t^2 + (2*z^2 + 4*z + 4)*t + z
        """
        u = f.ore_polynomial()
        if self._is_left:
            u = self._phi(a) * u
        else:
            u = u * self._phi(a)
        return f.parent()(u)


class DrinfeldModuleHomset(Homset):
    r"""
    This class implements the set of morphisms between two Drinfeld
    `\mathbb{F}_q[T]`-modules.

    INPUT:

    - ``X`` -- the domain

    - ``Y`` -- the codomain

    EXAMPLES::

        sage: Fq = GF(27)
        sage: A.<T> = Fq[]
        sage: K.<z6> = Fq.extension(2)
        sage: phi = DrinfeldModule(A, [z6, z6, 2])
        sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
        sage: H = Hom(phi, psi)
        sage: H
        Set of Drinfeld module morphisms
         from (gen) 2*t^2 + z6*t + z6
           to (gen) 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6

    ::

        sage: from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
        sage: isinstance(H, DrinfeldModuleHomset)
        True

    There is a simpler syntax for endomorphisms sets::

        sage: E = End(phi)
        sage: E
        Set of Drinfeld module morphisms from (gen) 2*t^2 + z6*t + z6 to (gen) 2*t^2 + z6*t + z6
        sage: E is Hom(phi, phi)
        True

    The domain and codomain must have the same Drinfeld modules
    category::

        sage: rho = DrinfeldModule(A, [Frac(A)(T), 1])
        sage: Hom(phi, rho)
        Traceback (most recent call last):
        ...
        ValueError: Drinfeld modules must be in the same category

    ::

        sage: sigma = DrinfeldModule(A, [1, z6, 2])
        sage: Hom(phi, sigma)
        Traceback (most recent call last):
        ...
        ValueError: Drinfeld modules must be in the same category

    One can create morphism objects by calling the homset::

        sage: identity_morphism = E(1)
        sage: identity_morphism
        Identity morphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6

    ::

        sage: t = phi.ore_polring().gen()
        sage: frobenius_endomorphism = E(t^6)
        sage: frobenius_endomorphism
        Endomorphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
          Defn: t^6

    ::

        sage: isogeny = H(t + 1)
        sage: isogeny
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
          To:   Drinfeld module defined by T |--> 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
          Defn: t + 1

    And one can test if an Ore polynomial defines a morphism using the
    ``in`` syntax::

        sage: 1 in H
        False
        sage: t^6 in H
        False
        sage: t + 1 in H
        True
        sage: 1 in E
        True
        sage: t^6 in E
        True
        sage: t + 1 in E
        False

    This also works if the candidate is a morphism object::

        sage: isogeny in H
        True
        sage: E(0) in E
        True
        sage: identity_morphism in H
        False
        sage: frobenius_endomorphism in H
        False
    """
    Element = DrinfeldModuleMorphism

    def __init__(self, X, Y, category=None, check=True):
        """
        Initialize ``self``.

        INPUT:

        - ``X`` -- the domain of the homset

        - ``Y`` -- the codomain of the homset

        - ``category`` -- (default: ``None``) the Drinfeld modules category of
          the domain and codomain

        - ``check`` -- boolean (default: ``True``); check the validity of the
          category

        TESTS::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: H.domain() is phi
            True
            sage: H.codomain() is psi
            True
        """
        if category is None:
            category = X.category()
        if check:
            if X.category() != Y.category() \
                    or not isinstance(X.category(), DrinfeldModules):
                raise ValueError('Drinfeld modules must be in the same category')
            if category != X.category():
                raise ValueError('category should be DrinfeldModules')
        base = category.base()
        super().__init__(X, Y, category=category, base=base, check=check)
        A = X.function_ring()
        self.register_action(DrinfeldModuleMorphismAction(A, self, True, operator.mul))
        # ARGH: the next line does not work
        # because Map.__mul__ does not call the coercion system
        self.register_action(DrinfeldModuleMorphismAction(A, self, False, operator.mul))
        if X is Y:
            self.register_coercion(A)

    def _latex_(self):
        r"""
        Return a LaTeX representation of the homset.

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: latex(H)
            \text{Set{ }of{ }Drinfeld{ }module{ }morphisms{ }from{ }(gen){ }}2 t^{2} + z_{6} t + z_{6}\text{{ }to{ }(gen){ }}2 t^{2} + \left(2 z_{6}^{5} + 2 z_{6}^{4} + 2 z_{6} + 1\right) t + z_{6}
        """
        return f'\\text{{Set{{ }}of{{ }}Drinfeld{{ }}module{{ }}morphisms' \
               f'{{ }}from{{ }}(gen){{ }}}}{latex(self.domain().gen())}' \
               f'\\text{{{{ }}to{{ }}(gen){{ }}}}'\
               f'{latex(self.codomain().gen())}'

    def _repr_(self):
        r"""
        Return a string representation of the homset.

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: H
            Set of Drinfeld module morphisms from (gen) 2*t^2 + z6*t + z6 to (gen) 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
        """
        return f'Set of Drinfeld module morphisms from (gen) '\
               f'{self.domain().gen()} to (gen) {self.codomain().gen()}'

    def __contains__(self, x):
        r"""
        Return ``True`` if the input defines a morphism in the homset.

        INPUT:

        - ``x`` -- an Ore polynomial or a Drinfeld module morphism

        EXAMPLES:

        In the next examples, the input is an Ore polynomial::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: E = End(phi)
            sage: t = phi.ore_polring().gen()
            sage: 1 in H
            False
            sage: t^6 in H
            False
            sage: t + 1 in H
            True
            sage: 1 in E
            True
            sage: t^6 in E
            True
            sage: t + 1 in E
            False

        Whereas the input is now a Drinfeld module morphism::

            sage: isogeny = H(t + 1)
            sage: isogeny in H
            True
            sage: E(0) in E
            True
            sage: identity_morphism = E(1)
            sage: identity_morphism in H
            False
            sage: frobenius_endomorphism = phi.frobenius_endomorphism()
            sage: frobenius_endomorphism in H
            False
        """
        try:
            x = self(x)
            return True
        except (AttributeError, ValueError, TypeError):
            return False

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return the Drinfeld module morphism defined by the given Ore
        polynomial.

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: E = End(phi)
            sage: t = phi.ore_polring().gen()
            sage: identity_morphism = E(1)
            sage: identity_morphism
            Identity morphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6

        ::

            sage: scalar_multiplication = E(T)
            sage: scalar_multiplication
            Endomorphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
              Defn: 2*t^2 + z6*t + z6

        ::

            sage: frobenius_endomorphism = E(t^6)
            sage: frobenius_endomorphism
            Endomorphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
              Defn: t^6

        ::

            sage: isogeny = H(t + 1)
            sage: isogeny
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
              To:   Drinfeld module defined by T |--> 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
              Defn: t + 1
        """
        # NOTE: This used to be self.element_class(self, ...), but this
        # would call __init__ instead of __classcall_private__. This
        # seems to work, but I don't know what I'm doing.
        return DrinfeldModuleMorphism(self, *args, **kwds)
