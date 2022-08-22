r"""
The left-module action induced by a Drinfeld module

This module provides the class
:class:`sage.rings.function_field.drinfeld_module.action.DrinfeldModuleAction`.

AUTHORS:

- Antoine Leudière (2022-04)
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
from sage.misc.latex import latex
from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule


class DrinfeldModuleAction(Action):
    r"""
    This class represents the left `\Fq[X]`-module action induced by a
    Drinfeld `\Fq[X]`-module defined over an `\Fq[X]`-field `K`.

    Let `L/K` be a field extension, let `x \in L`, let `P \in \Fq[X]`;
    the action is defined as `(P, a) \mapsto \phi_P(a)`, where
    `\phi_P(a)`. In this implementation, `L` is `K`.

    The action is instanciated as follows. Note that the user should
    never explicitely instanciate the class `DrinfeldModuleAction`::

    INPUT: a Drinfeld module

    EXAMPLES:

        sage: Fq.<z2> = GF(11)
        sage: FqX.<X> = Fq[]
        sage: K.<z> = Fq.extension(2)
        sage: phi = DrinfeldModule(FqX, [z, 0, 0, 1])
        sage: action = phi.action()
        sage: action
        Action on Finite Field in z of size 11^2 induced by Drinfeld module defined by X |--> t^3 + z over Finite Field in z of size 11^2

    The action on elements is computed as follows::

        sage: P = X + 1
        sage: a = z
        sage: action(P, a)
        ...
        4*z + 2
        sage: action(0, K.random_element())
        0
        sage: action(FqX.random_element(), 0)
        0

    To act on a field larger than `K`, one can change the ring of the
    Drinfeld module, then create the action::

        sage: extended_action = phi.change_ring(K.extension(2)).action()
        sage: extended_action
        Action on Finite Field in z4 of size 11^4 induced by Drinfeld module defined by X |--> t + 10*z4^3 + 4*z4^2 + 5*z4 + 5 over Finite Field in z4 of size 11^4

    Finally, given a Drinfeld module action, it is easy to recover the
    corresponding Drinfeld module::

        sage: action.drinfeld_module() is phi
        True
    """

    def __init__(self, drinfeld_module):
        if not isinstance(drinfeld_module, DrinfeldModule):
            raise TypeError('input must be a DrinfeldModule')
        self._drinfeld_module = drinfeld_module
        super().__init__(drinfeld_module.function_ring(),
                drinfeld_module.base_ring())

    def _act_(self, pol, x):
        r"""
        Return ``pol * x``, where ``*`` is the action.

        INPUT:

        - ``pol`` -- a polynomial in the function ring of the Drinfeld
          module
        - ``x`` -- an element in the base ring of the Drinfeld module

        OUTPUT: an element in the base ring of the Drinfeld module.

        EXAMPLES:

            sage: Fq.<z2> = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(2)
            sage: phi = DrinfeldModule(FqX, [z, 0, 0, 1])
            sage: action = phi.action()
            sage: P = X + 1
            sage: a = z
            sage: action(P, a)
            4*z + 2
            sage: action(0, K.random_element())
            0
            sage: action(FqX.random_element(), 0)
            0
        """
        if pol not in self._drinfeld_module.function_ring():
            raise TypeError('first input must be in the function ring')
        if x not in self._drinfeld_module.base_ring():
            raise TypeError('second input must be in the base ring')
        return self._drinfeld_module(pol)(x)

    def _latex_(self):
        r"""
        Return a LaTeX representation of the action.

        OUTPUT: a string

        EXAMPLES:

            sage: Fq.<z2> = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(2)
            sage: phi = DrinfeldModule(FqX, [z, 0, 0, 1])
            sage: action = phi.action()
            sage: latex(action)
            \text{Action{ }on{ }}\Bold{F}_{11^{2}}\text{{ }induced{ }by{ }}Drinfeld module defined by X |--> t^3 + z over Finite Field in z of size 11^2
        """
        return f'\\text{{Action{{ }}on{{ }}}}' \
                f'{latex(self._drinfeld_module.base_ring())}\\text{{{{ }}' \
                f'induced{{ }}by{{ }}}}{self._drinfeld_module}'

    def _repr_(self):
        r"""
        Return a string representation of the action.

        OUTPUT: a string

        EXAMPLES:

            sage: Fq.<z2> = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(2)
            sage: phi = DrinfeldModule(FqX, [z, 0, 0, 1])
            sage: action = phi.action()
            sage: action
            Action on Finite Field in z of size 11^2 induced by Drinfeld module defined by X |--> t^3 + z over Finite Field in z of size 11^2
        """
        return f'Action on {self._drinfeld_module.base_ring()} induced by ' \
                f'{self._drinfeld_module}'

    def drinfeld_module(self):
        r"""
        Return the Drinfeld module associated to the action.

        OUTPUT: a Drinfeld module

            sage: Fq.<z2> = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(2)
            sage: phi = DrinfeldModule(FqX, [z, 0, 0, 1])
            sage: action = phi.action()
            sage: action.drinfeld_module() is phi
            True
        """
        return self._drinfeld_module
