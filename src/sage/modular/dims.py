# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.pari
r"""
Dimensions of spaces of modular forms

AUTHORS:

- William Stein

- Jordi Quer

ACKNOWLEDGEMENT: The dimension formulas and implementations in this
module grew out of a program that Bruce Kaskel wrote (around 1996)
in PARI, which Kevin Buzzard subsequently extended. I (William
Stein) then implemented it in C++ for Hecke. I also implemented it
in Magma. Also, the functions for dimensions of spaces with
nontrivial character are based on a paper (that has no proofs) by
Cohen and Oesterlé [CO1977]_. The formulas for `\Gamma_H(N)` were found
and implemented by Jordi Quer.

The formulas here are more complete than in Hecke or Magma.

Currently the input to each function below is an integer and either a Dirichlet
character `\varepsilon` or a finite index subgroup of `\SL_2(\ZZ)`.
If the input is a Dirichlet character `\varepsilon`, the dimensions are for
subspaces of `M_k(\Gamma_1(N), \varepsilon)`, where `N` is the modulus of
`\varepsilon`.

These functions mostly call the methods ``dimension_cusp_forms``,
``dimension_modular_forms`` and so on of the corresponding congruence subgroup
classes.

REFERENCES:

.. [CO1977] \H. Cohen, J. Oesterlé, *Dimensions des espaces de formes
   modulaires*, p. 69-78 in Modular functions in one variable VI.
   Lecture Notes in Math. 627, Springer-Verlag, NewYork, 1977.
"""
# ****************************************************************************
#       Copyright (C) 2004-2008 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import factor, is_prime, valuation
from sage.misc.misc_c import prod
from sage.modular.arithgroup.all import (Gamma0, Gamma1, ArithmeticSubgroup,
                                         GammaH_class)
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer import Integer
from sage.rings.rational_field import frac

from . import dirichlet

##########################################################################
# Helper functions for calculating dimensions of spaces of modular forms
##########################################################################


def eisen(p):
    """
    Return the Eisenstein number `n` which is the numerator of `(p-1)/12`.

    INPUT:

    - ``p`` -- a prime

    OUTPUT: integer

    EXAMPLES::

        sage: [(p, sage.modular.dims.eisen(p)) for p in prime_range(24)]
        [(2, 1), (3, 1), (5, 1), (7, 1), (11, 5), (13, 1), (17, 4),
         (19, 3), (23, 11)]
    """
    if not is_prime(p):
        raise ValueError("p must be prime")
    return frac(p - 1, 12).numerator()

##########################################################################
# Formula of Cohen-Oesterlé for dim S_k(Gamma_1(N),eps).  REF:
# Springer Lecture notes in math, volume 627, pages 69--78.  The
# functions CO_delta and CO_nu, which were first written by Kevin
# Buzzard, are used only by the function CohenOesterle.
##########################################################################


def CO_delta(r, p, N, eps):
    r"""
    This is used as an intermediate value in computations related to
    the paper of Cohen-Oesterlé.

    INPUT:

    - ``r`` -- positive integer

    - ``p`` -- a prime

    - ``N`` -- positive integer

    - ``eps`` -- character

    OUTPUT: element of the base ring of the character

    EXAMPLES::

        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CO_delta(1,5,7,eps^3)
        2
    """
    K = eps.base_ring()
    if p % 4 == 3:
        return K.zero()
    if p == 2:
        if r == 1:
            return K.one()
        return K.zero()
    # interesting case: p=1(mod 4).
    # omega is a primitive 4th root of unity mod p.
    omega = (IntegerModRing(p).unit_gens()[0])**((p - 1) // 4)
    # this n is within a p-power root of a "local" 4th root of 1 modulo p.
    n = Mod(int(omega.crt(Mod(1, N // (p**r)))), N)
    n = n**(p**(r - 1))   # this is correct now
    t = eps(n)
    if t == K.one():
        return K(2)
    if t == -K.one():
        return K(-2)
    return K.zero()


def CO_nu(r, p, N, eps):
    r"""
    This is used as an intermediate value in computations related to
    the paper of Cohen-Oesterlé.

    INPUT:

    - ``r`` -- positive integer

    - ``p`` -- a prime

    - ``N`` -- positive integer

    - ``eps`` -- character

    OUTPUT: element of the base ring of the character

    EXAMPLES::

        sage: G.<eps> = DirichletGroup(7)
        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CO_nu(1,7,7,eps)
        -1
    """
    K = eps.base_ring()
    if p % 3 == 2:
        return K.zero()
    if p == 3:
        if r == 1:
            return K.one()
        return K.zero()
    # interesting case: p=1(mod 3)
    # omega is a cube root of 1 mod p.
    omega = (IntegerModRing(p).unit_gens()[0])**((p - 1) // 3)
    n = Mod(omega.crt(Mod(1, N // (p**r))), N)
    # within a p-power root of a "local" cube root of 1 mod p.
    n = n**(p**(r - 1))  # this is right now
    t = eps(n)
    if t == K.one():
        return K(2)
    return -K.one()


def CohenOesterle(eps, k):
    r"""
    Compute the Cohen-Oesterlé function associate to eps, `k`.

    This is a summand in the formula for the dimension of the space of
    cusp forms of weight `2` with character `\varepsilon`.

    INPUT:

    - ``eps`` -- Dirichlet character

    - ``k`` -- integer

    OUTPUT: element of the base ring of eps

    EXAMPLES::

        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CohenOesterle(eps, 2)
        -2/3
        sage: sage.modular.dims.CohenOesterle(eps, 4)
        -1
    """
    N = eps.modulus()
    facN = factor(N)
    f = eps.conductor()
    gamma_k = 0
    if k % 4 == 2:
        gamma_k = frac(-1, 4)
    elif k % 4 == 0:
        gamma_k = frac(1, 4)
    mu_k = 0
    if k % 3 == 2:
        mu_k = frac(-1, 3)
    elif k % 3 == 0:
        mu_k = frac(1, 3)

    def _lambda(r, s, p):
        """
        Used internally by the CohenOesterle function.

        INPUT:

        - ``r``, ``s``, ``p`` -- integers

        OUTPUT: integer

        EXAMPLES: (indirect doctest)

        ::

            sage: # needs sage.rings.number_field
            sage: K = CyclotomicField(3)
            sage: eps = DirichletGroup(7*43, K).0^2
            sage: sage.modular.dims.CohenOesterle(eps, 2)
            -4/3
        """
        if 2 * s <= r:
            if r % 2 == 0:
                return p**(r // 2) + p**((r // 2) - 1)
            return 2 * p**((r - 1) // 2)
        return 2 * p**(r - s)

    K = eps.base_ring()
    return K(frac(-1, 2) *
             prod(_lambda(r, valuation(f, p), p) for p, r in facN) +
             gamma_k * K.prod(CO_delta(r, p, N, eps) for p, r in facN) +
             mu_k * K.prod(CO_nu(r, p, N, eps) for p, r in facN))


####################################################################
# Functions exported to the global namespace.
# These have very flexible inputs.
####################################################################

def dimension_new_cusp_forms(X, k=2, p=0):
    """
    Return the dimension of the new (or `p`-new) subspace of
    cusp forms for the character or group `X`.

    INPUT:

    - ``X`` -- integer, congruence subgroup or Dirichlet
      character

    - ``k`` -- weight (integer)

    - ``p`` -- 0 or a prime

    EXAMPLES::

        sage: from sage.modular.dims import dimension_new_cusp_forms
        sage: dimension_new_cusp_forms(100,2)
        1

        sage: dimension_new_cusp_forms(Gamma0(100),2)
        1
        sage: dimension_new_cusp_forms(Gamma0(100),4)
        5

        sage: dimension_new_cusp_forms(Gamma1(100),2)
        141
        sage: dimension_new_cusp_forms(Gamma1(100),4)
        463

        sage: dimension_new_cusp_forms(DirichletGroup(100).1^2,2)
        2
        sage: dimension_new_cusp_forms(DirichletGroup(100).1^2,4)
        8

        sage: sum(dimension_new_cusp_forms(e,3) for e in DirichletGroup(30))
        12
        sage: dimension_new_cusp_forms(Gamma1(30),3)
        12

    Check that :issue:`12640` is fixed::

        sage: dimension_new_cusp_forms(DirichletGroup(1)(1), 12)
        1
        sage: dimension_new_cusp_forms(DirichletGroup(2)(1), 24)
        1
    """
    if isinstance(X, GammaH_class):
        return X.dimension_new_cusp_forms(k, p=p)
    elif isinstance(X, dirichlet.DirichletCharacter):
        N = X.modulus()
        if N <= 2:
            return Gamma0(N).dimension_new_cusp_forms(k, p=p)
        else:
            # Gamma1(N) for N<=2 just returns Gamma0(N), which has no
            # eps parameter. See trac #12640.
            return Gamma1(N).dimension_new_cusp_forms(k, eps=X, p=p)
    elif isinstance(X, (int, Integer)):
        return Gamma0(X).dimension_new_cusp_forms(k, p=p)
    raise TypeError(f"X (={X}) must be an integer, a Dirichlet character or a congruence subgroup of type Gamma0, Gamma1 or GammaH")


def dimension_cusp_forms(X, k=2):
    r"""
    The dimension of the space of cusp forms for the given congruence
    subgroup or Dirichlet character.

    INPUT:

    - ``X`` -- congruence subgroup or Dirichlet character
      or integer

    - ``k`` -- weight (integer)

    EXAMPLES::

        sage: from sage.modular.dims import dimension_cusp_forms
        sage: dimension_cusp_forms(5,4)
        1

        sage: dimension_cusp_forms(Gamma0(11),2)
        1
        sage: dimension_cusp_forms(Gamma1(13),2)
        2

        sage: dimension_cusp_forms(DirichletGroup(13).0^2,2)
        1
        sage: dimension_cusp_forms(DirichletGroup(13).0,3)
        1

        sage: dimension_cusp_forms(Gamma0(11),2)
        1
        sage: dimension_cusp_forms(Gamma0(11),0)
        0
        sage: dimension_cusp_forms(Gamma0(1),12)
        1
        sage: dimension_cusp_forms(Gamma0(1),2)
        0
        sage: dimension_cusp_forms(Gamma0(1),4)
        0

        sage: dimension_cusp_forms(Gamma0(389),2)
        32
        sage: dimension_cusp_forms(Gamma0(389),4)
        97
        sage: dimension_cusp_forms(Gamma0(2005),2)
        199
        sage: dimension_cusp_forms(Gamma0(11),1)
        0

        sage: dimension_cusp_forms(Gamma1(11),2)
        1
        sage: dimension_cusp_forms(Gamma1(1),12)
        1
        sage: dimension_cusp_forms(Gamma1(1),2)
        0
        sage: dimension_cusp_forms(Gamma1(1),4)
        0

        sage: dimension_cusp_forms(Gamma1(389),2)
        6112
        sage: dimension_cusp_forms(Gamma1(389),4)
        18721
        sage: dimension_cusp_forms(Gamma1(2005),2)
        159201

        sage: dimension_cusp_forms(Gamma1(11),1)
        0

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_cusp_forms(e,2)
        0
        sage: dimension_cusp_forms(e^2,2)
        1

    Check that :issue:`12640` is fixed::

        sage: dimension_cusp_forms(DirichletGroup(1)(1), 12)
        1
        sage: dimension_cusp_forms(DirichletGroup(2)(1), 24)
        5
    """
    if isinstance(X, dirichlet.DirichletCharacter):
        N = X.modulus()
        if N <= 2:
            return Gamma0(N).dimension_cusp_forms(k)
        else:
            return Gamma1(N).dimension_cusp_forms(k, X)
    elif isinstance(X, ArithmeticSubgroup):
        return X.dimension_cusp_forms(k)
    elif isinstance(X, (int, Integer)):
        return Gamma0(X).dimension_cusp_forms(k)
    raise TypeError("argument 1 must be a Dirichlet character, an integer "
                    "or a finite index subgroup of SL2Z")


def dimension_eis(X, k=2):
    """
    The dimension of the space of Eisenstein series for the given
    congruence subgroup.

    INPUT:

    - ``X`` -- congruence subgroup or Dirichlet character
      or integer

    - ``k`` -- integer; weight

    EXAMPLES::

        sage: from sage.modular.dims import dimension_eis
        sage: dimension_eis(5,4)
        2

        sage: dimension_eis(Gamma0(11),2)
        1
        sage: dimension_eis(Gamma1(13),2)
        11
        sage: dimension_eis(Gamma1(2006),2)
        3711

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_eis(e,2)
        0
        sage: dimension_eis(e^2,2)
        2

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_eis(e,2)
        0
        sage: dimension_eis(e^2,2)
        2
        sage: dimension_eis(e,13)
        2

        sage: G = DirichletGroup(20)
        sage: dimension_eis(G.0,3)
        4
        sage: dimension_eis(G.1,3)
        6
        sage: dimension_eis(G.1^2,2)
        6

        sage: G = DirichletGroup(200)
        sage: e = prod(G.gens(), G(1))
        sage: e.conductor()
        200
        sage: dimension_eis(e,2)
        4

        sage: from sage.modular.dims import dimension_modular_forms
        sage: dimension_modular_forms(Gamma1(4), 11)
        6
    """
    if isinstance(X, ArithmeticSubgroup):
        return X.dimension_eis(k)
    elif isinstance(X, dirichlet.DirichletCharacter):
        return Gamma1(X.modulus()).dimension_eis(k, X)
    elif isinstance(X, (int, Integer)):
        return Gamma0(X).dimension_eis(k)
    raise TypeError(f"argument in dimension_eis must be an integer, a Dirichlet character, or a finite index subgroup of SL2Z (got {X})")


def dimension_modular_forms(X, k=2):
    r"""
    The dimension of the space of cusp forms for the given congruence
    subgroup (either `\Gamma_0(N)`, `\Gamma_1(N)`, or
    `\Gamma_H(N)`) or Dirichlet character.

    INPUT:

    - ``X`` -- congruence subgroup or Dirichlet character

    - ``k`` -- integer; weight

    EXAMPLES::

        sage: from sage.modular.dims import dimension_modular_forms
        sage: dimension_modular_forms(Gamma0(11),2)
        2
        sage: dimension_modular_forms(Gamma0(11),0)
        1
        sage: dimension_modular_forms(Gamma1(13),2)
        13
        sage: dimension_modular_forms(GammaH(11, [10]), 2)
        10
        sage: dimension_modular_forms(GammaH(11, [10]))
        10
        sage: dimension_modular_forms(GammaH(11, [10]), 4)
        20
        sage: e = DirichletGroup(20).1
        sage: dimension_modular_forms(e,3)
        9
        sage: from sage.modular.dims import dimension_cusp_forms
        sage: dimension_cusp_forms(e,3)
        3
        sage: from sage.modular.dims import dimension_eis
        sage: dimension_eis(e,3)
        6
        sage: dimension_modular_forms(11,2)
        2
    """
    if isinstance(X, (int, Integer)):
        return Gamma0(X).dimension_modular_forms(k)
    elif isinstance(X, ArithmeticSubgroup):
        return X.dimension_modular_forms(k)
    elif isinstance(X, dirichlet.DirichletCharacter):
        return Gamma1(X.modulus()).dimension_modular_forms(k, eps=X)
    else:
        raise TypeError("argument 1 must be an integer, a Dirichlet character "
                        "or an arithmetic subgroup")


def sturm_bound(level, weight=2):
    r"""
    Return the Sturm bound for modular forms with given level and weight.

    For more details, see the documentation for the ``sturm_bound`` method
    of :class:`sage.modular.arithgroup.CongruenceSubgroup` objects.

    INPUT:

    - ``level`` -- integer (interpreted as a level for `\Gamma0`) or a
      congruence subgroup

    - ``weight`` -- integer `\geq 2` (default: 2)

    EXAMPLES::

        sage: from sage.modular.dims import sturm_bound
        sage: sturm_bound(11,2)
        2
        sage: sturm_bound(389,2)
        65
        sage: sturm_bound(1,12)
        1
        sage: sturm_bound(100,2)
        30
        sage: sturm_bound(1,36)
        3
        sage: sturm_bound(11)
        2
    """
    if isinstance(level, ArithmeticSubgroup):
        if level.is_congruence():
            return level.sturm_bound(weight)
        raise ValueError("no Sturm bound defined for noncongruence subgroups")
    if isinstance(level, (int, Integer)):
        return Gamma0(level).sturm_bound(weight)
