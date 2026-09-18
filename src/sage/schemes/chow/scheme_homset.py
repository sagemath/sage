# -*- coding: utf-8 -*-
r"""
The set of morphisms between ChowSchemes

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""
# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import HomsetWithBase
from sage.structure.factory import UniqueFactory
from sage.schemes.chow.scheme import PointChowScheme
from sage.schemes.chow.morphism import ChowSchemeMorphism
from sage.schemes.chow.scheme import is_chowScheme
from sage.categories.map import Map
from sage.categories.all import Rings
_Rings = Rings() # type: ignore


def is_ChowSchemeHomset(H):
    r"""
    Test whether ``H`` is a chowscheme SHom-set.
    """
    return isinstance(H, ChowSchemeHomset_generic)


class ChowSchemeHomsetFactory(UniqueFactory):
    r"""
    Factory for SHom-sets of schemes.
    """

    def create_key_and_extra_args(self, X, Y, category=None,
                                  base=None, check=True):
        r"""
        Create a key that uniquely determines the SHom-set.

        INPUT:

        - ``X`` -- a chowscheme. The domain of the morphisms.

        - ``Y`` -- a chowscheme. The codomain of the morphisms.

        - ``category`` -- a category for the SHom-sets (default: chowschemes
            over given base).

        - ``base`` -- a chowscheme. The base chowscheme of domain and codomain.

        - ``check`` -- boolean (default: ``True``).
        """
        if base is None:
            base = PointChowScheme
        if not is_chowScheme(base):
            raise ValueError("Expect a ChowScheme as base")
        if not category:
            from .schemes import ChowSchemes
            category = ChowSchemes(base)
        key = tuple([id(X), id(Y), category])
        extra = {'X': X, 'Y': Y, 'base_chowscheme': base, 'check': check}
        return key, extra

    def create_object(self, version, key, **extra_args):
        r"""
        Create a :class:`SchemeHomset_generic`.

        INPUT:

        - ``version`` -- object version. Currently not used.

        - ``key`` -- a key created by :meth:`create_key_and_extra_args`.

        - ``extra_args`` -- a dictionary of extra keyword arguments.
        """
        category = key[2]
        X = extra_args.pop('X')
        Y = extra_args.pop('Y')
        base_chowscheme = extra_args.pop('base_chowscheme')
        try:
            return X._homset(X, Y, category=category, base=base_chowscheme,
                             **extra_args)
        except AttributeError:
            return ChowSchemeHomset_generic(X, Y, category=category,
                                            base=base_chowscheme, **extra_args)


ChowSchemeHomset = ChowSchemeHomsetFactory('chow.scheme_homset.ChowSchemeHomset')


class ChowSchemeHomset_generic(HomsetWithBase):
    r"""
    The base class for SHom-sets of ChowSchemes.

    INPUT:

    - ``X`` -- a chowscheme. The domain of the SHom-set.

    - ``Y`` -- a chowscheme. The codomain of the SHom-set.

    - ``category`` -- a category (optional). The category of the
        SHom-set.

    - ``check`` -- boolean (optional, default=``True``). Whether to
        check the defining data for consistency.
    """
    Element = ChowSchemeMorphism

    def __call__(self, *args, **kwds):
        r"""
        Make SHom-sets callable.

        See the ``_call_()`` method of the derived class. All
        arguments are handed through.
        """
        from sage.structure.parent import Set_generic
        return Set_generic.__call__(self, *args, **kwds)

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string.
        """
        s = 'Set of morphisms'
        s += '\n  From: %s' % self.domain()
        s += '\n  To:   %s' % self.codomain()
        return s

    def natural_map(self):
        r"""
        Return a natural map in the SHom space.

        OUTPUT:

        A :class:`ChowSchemeMorphism` if there is a natural map from
        domain to codomain. Otherwise, a ``NotImplementedError`` is
        raised.
        """
        X = self.domain()
        Y = self.codomain()
        if Y == X.base_chowscheme():
            return X.base_morphism()
        raise NotImplementedError

    def _element_constructor_(self, x, check=True):
        r"""
        Construct a chowscheme morphism.

        INPUT:

        - `x` -- a ring morphism, or a list or a tuple that define a
          ring morphism.

        - ``check`` -- boolean (default: ``True``) passed onto
          functions called by this one to be more careful about input
          argument type checking.
        """
        if isinstance(x, (list, tuple)):
            return self.domain()._morphism(self, x)

        if isinstance(x, str):
            return self.domain()._morphism(self, x)

        if isinstance(x, Map) and x.category_for().is_subcategory(_Rings):
            return ChowSchemeMorphism(self, x)

        raise TypeError("x must be a ring homomorphism, list or tuple")

    def zero_element(self):
        r"""
        Backward compatibility alias for self.zero()

        """
        n = self.codomain().ngens() if self.codomain() else 0
        return self([0] * n)
