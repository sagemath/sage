# -*- coding: utf-8 -*-
r"""
Modular Forms with Character

EXAMPLES::

    sage: eps = DirichletGroup(13).0
    sage: M = ModularForms(eps^2, 2); M
    Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2

    sage: S = M.cuspidal_submodule(); S
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2
    sage: S.modular_symbols()
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2

We create a spaces associated to Dirichlet characters of modulus 225::

    sage: e = DirichletGroup(225).0
    sage: e.order()
    6
    sage: e.base_ring()
    Cyclotomic Field of order 60 and degree 16
    sage: M = ModularForms(e,3)

Notice that the base ring is "minimized"::

    sage: M
    Modular Forms space of dimension 66, character [zeta6, 1] and weight 3
    over Cyclotomic Field of order 6 and degree 2

If we don't want the base ring to change, we can explicitly specify it::

    sage: ModularForms(e, 3, e.base_ring())
    Modular Forms space of dimension 66, character [zeta6, 1] and weight 3
    over Cyclotomic Field of order 60 and degree 16

Next we create a space associated to a Dirichlet character of order 20::

    sage: e = DirichletGroup(225).1
    sage: e.order()
    20
    sage: e.base_ring()
    Cyclotomic Field of order 60 and degree 16
    sage: M = ModularForms(e,17); M
    Modular Forms space of dimension 484, character [1, zeta20] and
    weight 17 over Cyclotomic Field of order 20 and degree 8

We compute the Eisenstein subspace, which is fast even though
the dimension of the space is large (since an explicit basis
of `q`-expansions has not been computed yet).

::

    sage: M.eisenstein_submodule()
    Eisenstein subspace of dimension 8 of Modular Forms space of
    dimension 484, character [1, zeta20] and weight 17 over Cyclotomic Field of order 20 and degree 8

    sage: M.cuspidal_submodule()
    Cuspidal subspace of dimension 476 of Modular Forms space of dimension 484, character [1, zeta20] and weight 17 over Cyclotomic Field of order 20 and degree 8

TESTS::

    sage: m = ModularForms(DirichletGroup(20).1,5)
    sage: m == loads(dumps(m))
    True
    sage: type(m)
    <class 'sage.modular.modform.ambient_eps.ModularFormsAmbient_eps_with_category'>
"""
from __future__ import absolute_import

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import sage.rings.all as rings

import sage.modular.arithgroup.all as arithgroup
import sage.modular.dirichlet as dirichlet
import sage.modular.modsym.modsym as modsym

from .ambient import ModularFormsAmbient

from . import ambient_R
from . import cuspidal_submodule
from . import eisenstein_submodule

class ModularFormsAmbient_eps(ModularFormsAmbient):
    """
    A space of modular forms with character.
    """
    def __init__(self, character, weight=2, base_ring=None):
        """
        Create an ambient modular forms space with character.

        .. note::

           The base ring must be of characteristic 0.  The ambient_R
           Python module is used for computing in characteristic p,
           which we view as the reduction of characteristic 0.

        INPUT:

        - ``weight`` - int

        - ``character`` - dirichlet.DirichletCharacter

        - ``base_ring`` - base field

        EXAMPLES::

            sage: m = ModularForms(DirichletGroup(11).0,3); m
            Modular Forms space of dimension 3, character [zeta10] and weight 3 over Cyclotomic Field of order 10 and degree 4
            sage: type(m)
            <class 'sage.modular.modform.ambient_eps.ModularFormsAmbient_eps_with_category'>
        """
        if not dirichlet.is_DirichletCharacter(character):
            raise TypeError("character (=%s) must be a Dirichlet character"%character)
        if base_ring is None: base_ring=character.base_ring()
        if character.base_ring() != base_ring:
            character = character.change_ring(base_ring)
        if base_ring.characteristic() != 0:
            raise ValueError("the base ring must have characteristic 0.")
        group = arithgroup.Gamma1(character.modulus())
        base_ring = character.base_ring()
        ModularFormsAmbient.__init__(self, group, weight, base_ring, character)

    def _repr_(self):
        """
        String representation of this space with character.

        EXAMPLES::

            sage: m = ModularForms(DirichletGroup(8).1,2)
            sage: m._repr_()
            'Modular Forms space of dimension 2, character [1, -1] and weight 2 over Rational Field'

        You can rename the space with the rename command.

        ::

            sage: m.rename('Modforms of level 8')
            sage: m
            Modforms of level 8
        """
        try:
            d = self.dimension()
        except NotImplementedError:
            d = "(unknown)"
        return "Modular Forms space of dimension %s, character %s and weight %s over %s"%(
            d, self.character()._repr_short_(), self.weight(), self.base_ring())

    def cuspidal_submodule(self):
        """
        Return the cuspidal submodule of this ambient space of modular forms.

        EXAMPLES::

            sage: eps = DirichletGroup(4).0
            sage: M = ModularForms(eps, 5); M
            Modular Forms space of dimension 3, character [-1] and weight 5 over Rational Field
            sage: M.cuspidal_submodule()
            Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3, character [-1] and weight 5 over Rational Field
        """
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_eps(self)
        return self.__cuspidal_submodule

    def change_ring(self, base_ring):
        """
        Return space with same defining parameters as this ambient
        space of modular symbols, but defined over a different base
        ring.

        EXAMPLES::

            sage: m = ModularForms(DirichletGroup(13).0^2,2); m
            Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2
            sage: m.change_ring(CyclotomicField(12))
            Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 12 and degree 4

        It must be possible to change the ring of the underlying Dirichlet character::

            sage: m.change_ring(QQ)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce zeta6 to a rational
        """
        if self.base_ring() == base_ring:
            return self
        return ambient_R.ModularFormsAmbient_R(self, base_ring = base_ring)

    def modular_symbols(self, sign=0):
        """
        Return corresponding space of modular symbols with given sign.

        EXAMPLES::

            sage: eps = DirichletGroup(13).0
            sage: M = ModularForms(eps^2, 2)
            sage: M.modular_symbols()
            Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
            sage: M.modular_symbols(1)
            Modular Symbols space of dimension 3 and level 13, weight 2, character [zeta6], sign 1, over Cyclotomic Field of order 6 and degree 2
            sage: M.modular_symbols(-1)
            Modular Symbols space of dimension 1 and level 13, weight 2, character [zeta6], sign -1, over Cyclotomic Field of order 6 and degree 2
            sage: M.modular_symbols(2)
            Traceback (most recent call last):
            ...
            ValueError: sign must be -1, 0, or 1
        """
        sign = rings.Integer(sign)
        try:
            return self.__modsym[sign]
        except AttributeError:
            self.__modsym = {}
        except KeyError:
            pass
        self.__modsym[sign] = modsym.ModularSymbols(\
            self.character(),
            weight = self.weight(),
            sign = sign,
            base_ring = self.base_ring())
        return self.__modsym[sign]

    def eisenstein_submodule(self):
        """
        Return the submodule of this ambient module with character that is
        spanned by Eisenstein series.  This is the Hecke stable complement
        of the cuspidal submodule.

        EXAMPLES::

            sage: m = ModularForms(DirichletGroup(13).0^2,2); m
            Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2
            sage: m.eisenstein_submodule()
            Eisenstein subspace of dimension 2 of Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2
        """
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = eisenstein_submodule.EisensteinSubmodule_eps(self)
        return self.__eisenstein_submodule

    def hecke_module_of_level(self, N):
        r"""
        Return the Hecke module of level N corresponding to self, which is the
        domain or codomain of a degeneracy map from self. Here N must be either
        a divisor or a multiple of the level of self, and a multiple of the
        conductor of the character of self.

        EXAMPLES::

            sage: M = ModularForms(DirichletGroup(15).0, 3); M.character().conductor()
            3
            sage: M.hecke_module_of_level(3)
            Modular Forms space of dimension 2, character [-1] and weight 3 over Rational Field
            sage: M.hecke_module_of_level(5)
            Traceback (most recent call last):
            ...
            ValueError: conductor(=3) must divide M(=5)
            sage: M.hecke_module_of_level(30)
            Modular Forms space of dimension 16, character [-1, 1] and weight 3 over Rational Field
        """
        from . import constructor
        if N % self.level() == 0:
            return constructor.ModularForms(self.character().extend(N), self.weight(), self.base_ring(), prec=self.prec())
        elif self.level() % N == 0:
            return constructor.ModularForms(self.character().restrict(N), self.weight(), self.base_ring(), prec=self.prec())
        else:
            raise ValueError("N (=%s) must be a divisor or a multiple of the level of self (=%s)" % (N, self.level()))

    ####################################################################
    # Computations of Dimensions
    ####################################################################
    def _dim_cuspidal(self):
        """
        Return the dimension of the cuspidal subspace, computed using a dimension formula.

        EXAMPLES::

            sage: m = ModularForms(DirichletGroup(389,CyclotomicField(4)).0,3); m._dim_cuspidal()
            64
        """
        try:
            return self.__the_dim_cuspidal
        except AttributeError:
            self.__the_dim_cuspidal = self.group().dimension_cusp_forms(
                self.weight(), eps=self.character())
        return self.__the_dim_cuspidal

    def _dim_eisenstein(self):
        """
        Return the dimension of the Eisenstein subspace of this space, computed
        using a dimension formula.

        EXAMPLES::

            sage: m = ModularForms(DirichletGroup(13).0,7); m
            Modular Forms space of dimension 8, character [zeta12] and weight 7 over Cyclotomic Field of order 12 and degree 4
            sage: m._dim_eisenstein()
            2
            sage: m._dim_cuspidal()
            6
        """
        try:
            return self.__the_dim_eisenstein
        except AttributeError:
            self.__the_dim_eisenstein = self.group().dimension_eis(self.weight(), eps=self.character())
        return self.__the_dim_eisenstein

    def _dim_new_cuspidal(self):
        """
        Return the dimension of the new cuspidal subspace, computed
        using a dimension formula.

        EXAMPLES::

            sage: m = ModularForms(DirichletGroup(33).0,7); m
            Modular Forms space of dimension 26, character [-1, 1] and weight 7 over Rational Field
            sage: m._dim_new_cuspidal()
            20
            sage: m._dim_cuspidal()
            22
        """
        try:
            return self.__the_dim_new_cuspidal
        except AttributeError:
            self.__the_dim_new_cuspidal = self.group().dimension_new_cusp_forms(self.weight(), eps=self.character())
        return self.__the_dim_new_cuspidal

    def _dim_new_eisenstein(self):
        """
        Return the dimension of the new Eisenstein subspace, computed
        by enumerating all Eisenstein series of the appropriate level.

        EXAMPLES::

            sage: m = ModularForms(DirichletGroup(36).0,5); m
            Modular Forms space of dimension 28, character [-1, 1] and weight 5 over Rational Field
            sage: m._dim_new_eisenstein()
            2
            sage: m._dim_eisenstein()
            8
        """
        try:
            return self.__the_dim_new_eisenstein
        except AttributeError:
            if self.character().is_trivial() and self.weight() == 2:
                if rings.is_prime(self.level()):
                    d = 1
                else:
                    d = 0
            else:
                E = self.eisenstein_series()
                d = len([g for g in E if g.new_level() == self.level()])
            self.__the_dim_new_eisenstein = d
        return self.__the_dim_new_eisenstein


