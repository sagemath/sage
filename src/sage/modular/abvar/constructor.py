"""
Constructors for certain modular abelian varieties

AUTHORS:

- William Stein (2007-03)
"""
###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

import weakref

from sage.rings.integer import Integer

from sage.modular.arithgroup.all import CongruenceSubgroupBase, Gamma0
from sage.modular.modsym.space import ModularSymbolsSpace
from .abvar_newform import ModularAbelianVariety_newform
import sage.modular.modform.element
from sage.modular.abvar import abvar

_cache = {}


def _get(key):
    """
    Return the cached abelian variety with given key. This is used
    internally by the abelian varieties constructor.

    INPUT:

    - ``key`` -- hashable

    EXAMPLES::

        sage: sage.modular.abvar.constructor._saved('a', J0(37))
        Abelian variety J0(37) of dimension 2
        sage: sage.modular.abvar.constructor._get('a')
        Abelian variety J0(37) of dimension 2
        sage: sage.modular.abvar.constructor._get('b')
        Traceback (most recent call last):
        ...
        ValueError: element not in cache
    """
    if key in _cache:
        z = _cache[key]()
        if z is not None:
            return z
    raise ValueError("element not in cache")


def _saved(key, J):
    """
    Return the cached abelian variety with given key. This is used
    internally by the abelian varieties constructor.

    INPUT:

    - ``key`` -- hashable

    - ``J`` -- modular abelian variety

    OUTPUT: ``J`` -- returns the modabvar, to make code that uses
    this simpler

    EXAMPLES::

        sage: sage.modular.abvar.constructor._saved('37', J0(37))
        Abelian variety J0(37) of dimension 2
    """
    _cache[key] = weakref.ref(J)
    return J


def J0(N):
    """
    Return the Jacobian `J_0(N)` of the modular curve
    `X_0(N)`.

    EXAMPLES::

        sage: J0(389)
        Abelian variety J0(389) of dimension 32

    The result is cached::

        sage: J0(33) is J0(33)
        True
    """
    key = 'J0(%s)' % N
    try:
        return _get(key)
    except ValueError:
        from sage.modular.arithgroup.all import Gamma0
        J = Gamma0(N).modular_abelian_variety()
        return _saved(key, J)


def J1(N):
    """
    Return the Jacobian `J_1(N)` of the modular curve
    `X_1(N)`.

    EXAMPLES::

        sage: J1(389)
        Abelian variety J1(389) of dimension 6112
    """
    key = 'J1(%s)' % N
    try:
        return _get(key)
    except ValueError:
        from sage.modular.arithgroup.all import Gamma1
        return _saved(key, Gamma1(N).modular_abelian_variety())


def JH(N, H):
    """
    Return the Jacobian `J_H(N)` of the modular curve
    `X_H(N)`.

    EXAMPLES::

        sage: JH(389,[16])
        Abelian variety JH(389,[16]) of dimension 64
    """
    key = 'JH(%s,%s)' % (N,H)
    try:
        return _get(key)
    except ValueError:
        from sage.modular.arithgroup.all import GammaH
        return _saved(key, GammaH(N, H).modular_abelian_variety())


def AbelianVariety(X):
    """
    Create the abelian variety corresponding to the given defining
    data.

    INPUT:

    - ``X`` -- integer, string, newform, modsym space,
      congruence subgroup or tuple of congruence subgroups

    OUTPUT: a modular abelian variety

    EXAMPLES::

        sage: AbelianVariety(Gamma0(37))
        Abelian variety J0(37) of dimension 2
        sage: AbelianVariety('37a')
        Newform abelian subvariety 37a of dimension 1 of J0(37)
        sage: AbelianVariety(Newform('37a'))
        Newform abelian subvariety 37a of dimension 1 of J0(37)
        sage: AbelianVariety(ModularSymbols(37).cuspidal_submodule())
        Abelian variety J0(37) of dimension 2
        sage: AbelianVariety((Gamma0(37), Gamma0(11)))
        Abelian variety J0(37) x J0(11) of dimension 3
        sage: AbelianVariety(37)
        Abelian variety J0(37) of dimension 2
        sage: AbelianVariety([1,2,3])
        Traceback (most recent call last):
        ...
        TypeError: X must be an integer, string, newform, modsym space, congruence subgroup or tuple of congruence subgroups
    """
    if isinstance(X, (int, Integer)):
        X = Gamma0(X)
    if isinstance(X, CongruenceSubgroupBase):
        X = X.modular_symbols().cuspidal_submodule()
    elif isinstance(X, str):
        from sage.modular.modform.constructor import Newform
        f = Newform(X, names='a')
        return ModularAbelianVariety_newform(f, internal_name=True)
    elif isinstance(X, sage.modular.modform.element.Newform):
        return ModularAbelianVariety_newform(X)

    if isinstance(X, ModularSymbolsSpace):
        return abvar.ModularAbelianVariety_modsym(X)

    if isinstance(X, (tuple,list)) and all(isinstance(G, CongruenceSubgroupBase) for G in X):
        return abvar.ModularAbelianVariety(X)

    raise TypeError("X must be an integer, string, newform, modsym space, congruence subgroup or tuple of congruence subgroups")
