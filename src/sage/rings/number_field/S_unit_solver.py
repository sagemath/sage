# sage.doctest: needs sage.rings.number_field sage.rings.padics
r"""
Solver for the `S`-unit equation `x + y = 1`

Inspired by works of Tzanakis--de Weger, Baker--Wustholz and Smart, we use the
LLL methods to implement an algorithm that returns all `S`-unit solutions to
the equation `x + y = 1`.

EXAMPLES::

    sage: from sage.rings.number_field.S_unit_solver import solve_S_unit_equation, eq_up_to_order
    sage: x = polygen(ZZ, 'x')
    sage: K.<xi> = NumberField(x^2 + x + 1)
    sage: S = K.primes_above(3)
    sage: expected = [((0, 1), (4, 0), xi + 2, -xi - 1),
    ....:             ((1, -1), (0, -1), 1/3*xi + 2/3, -1/3*xi + 1/3),
    ....:             ((1, 0), (5, 0), xi + 1, -xi),
    ....:             ((2, 0), (5, 1), xi, -xi + 1)]
    sage: sols = solve_S_unit_equation(K, S, 200)
    sage: eq_up_to_order(sols, expected)
    True

.. TODO::

    - Use Cython to improve timings on the sieve

REFERENCES:

- [MR2016]_
- [Sma1995]_
- [Sma1998]_
- [Yu2007]_
- [AKMRVW]_

AUTHORS:

- Alejandra Alvarado, Angelos Koutsianas, Beth Malmskog, Christopher Rasmussen,
  David Roe, Christelle Vincent, Mckenzie West (2018-04-25 to 2018-11-09):
  original version
"""


# ****************************************************************************
#       Copyright (C) 2020 Alejandra Alvarado <aalvarado2 at eiu.edu>
#                          Angelos Koutsianas <koutsis.jr at gmail.com>
#                          Beth Malmskog <beth.malmskog at gmail.com>
#                          Christopher Rasmussen <crasmussen at wesleyan.edu>
#                          Christelle Vincent <christelle.vincent at uvm.edu>
#                          Mckenzie West <westmr at uwec.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.rings.infinity import Infinity
from sage.symbolic.ring import SR
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RealField, RR
from sage.rings.complex_mpfr import ComplexField
from sage.functions.log import exp
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import is_real_place, refine_embedding
from sage.rings.number_field.unit_group import UnitGroup
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.finite_rings.integer_mod import mod
from sage.rings.padics.factory import Qp
from sage.combinat.combination import Combinations
from sage.misc.misc_c import prod
from sage.arith.functions import lcm
from sage.arith.misc import gcd, CRT, factorial
from sage.matrix.constructor import matrix, identity_matrix, vector, block_matrix, zero_matrix
from sage.modules.free_module_element import zero_vector
from itertools import combinations_with_replacement
from copy import copy
import itertools


def column_Log(SUK, iota, U, prec=106):
    r"""
    Return the log vector of ``iota``; i.e., the logs of all the valuations.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``iota`` -- an element of ``K``
    - ``U`` -- list of places (finite or infinite) of ``K``
    - ``prec`` -- the precision of the real field (default: 106)

    OUTPUT: the log vector as a list of real numbers

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import column_Log
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: S = tuple(K.primes_above(3))
        sage: SUK = UnitGroup(K, S=S)
        sage: phi_complex = K.places()[1]
        sage: v_fin = S[0]
        sage: U = [phi_complex, v_fin]
        sage: column_Log(SUK, xi^2, U)  # abs tol 1e-29
        [1.464816384890812968648768625966, -2.197224577336219382790490473845]

    REFERENCES:

    - [Sma1995]_ p. 823
    """
    R = RealField(prec)

    return [ R(SUK.number_field().abs_val(v, iota, prec)).log() for v in U]


def c3_func(SUK, prec=106):
    r"""
    Return the constant `c_3` from [AKMRVW]_.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``prec`` -- the precision of the real field (default: 106)

    OUTPUT: the constant `c_3`, as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c3_func
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))

        sage: c3_func(SUK)  # abs tol 1e-29
        0.4257859134798034746197327286726

    .. NOTE::

        The numerator should be as close to 1 as possible, especially as the rank of the `S`-units grows large

    REFERENCES:

    - [AKMRVW]_ :arxiv:`1903.00977`
    """
    R = RealField(prec)

    all_places = list(SUK.primes()) + SUK.number_field().places(prec)
    Possible_U = Combinations(all_places, SUK.rank())
    c1 = R(1)  # guarantees final c1 >= 1
    for U in Possible_U:
        # first, build the matrix C_{i,U}
        columns_of_C = [column_Log(SUK, unit, U, prec)
                        for unit in SUK.fundamental_units()]
        C = matrix(SUK.rank(), SUK.rank(), columns_of_C)
        # Is it invertible?
        if abs(C.determinant()) > 10**(-10):
            poss_c1 = C.inverse().apply_map(abs).norm(Infinity)
            c1 = R(max(poss_c1, c1))
    return R(0.9999999) / (c1 * SUK.rank())


def c4_func(SUK, v, A, prec=106):
    r"""
    Return the constant `c_4` from Smart's TCDF paper, [Sma1995]_.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of ``SUK.number_field().places(prec)``)
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` -- the precision of the real field (default: 106)

    OUTPUT: the constant `c_4`, as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c4_func
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c4_func(SUK, phi_real, A)
        1.000000000000000000000000000000

        sage: c4_func(SUK, phi_complex, A)
        1.000000000000000000000000000000

        sage: c4_func(SUK, v_fin, A)
        1.000000000000000000000000000000

    REFERENCES:

    - [Sma1995]_ p. 824
    """
    return max(SUK.number_field().abs_val(v, alpha, prec) for alpha in A)


def beta_k(betas_and_ns):
    r"""
    Return a pair `[\beta_k,|beta_k|_v]`, where `\beta_k` has the smallest
    nonzero valuation in absolute value of the list ``betas_and_ns``.

    INPUT:

    - ``betas_and_ns`` -- list of pairs ``[beta,val_v(beta)]`` outputted from
      the function where ``beta`` is an element of ``SUK.fundamental_units()``

    OUTPUT:

    The pair ``[beta_k,v(beta_k)]``, where ``beta_k`` is an element of ``K``
    and ``val_v(beta_k)`` is a integer.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import beta_k
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: betas = [[beta, beta.valuation(v_fin)] for beta in SUK.fundamental_units()]
        sage: beta_k(betas)
        [xi, 1]

    REFERENCES:

    - [Sma1995]_ pp. 824-825
    """
    for pair in betas_and_ns:
        if abs( pair[1] ) != 0:
            good_pair = pair
            break
    for pair in betas_and_ns:
        if ( abs(pair[1]) != 0 and abs(pair[1]) < abs(good_pair[1]) ):
            good_pair = pair
    return good_pair


def mus(SUK, v):
    r"""
    Return a list `[\mu]`, for `\mu` defined in [AKMRVW]_.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a finite place of ``K``

    OUTPUT: list ``[mus]`` where each ``mu`` is an element of ``K``

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import mus
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: mus(SUK, v_fin)
        [xi^2 - 2]

    REFERENCES:

    - [AKMRVW]_
    """
    betas = SUK.fundamental_units()
    beta_and_ns = [[beta,beta.valuation(v)] for beta in betas]
    if all(pair[1] == 0 for pair in beta_and_ns):
        return betas
    else:
        good_pair = beta_k(beta_and_ns)
        temp = [(beta[0]**good_pair[1])*(good_pair[0]**(-beta[1])) for beta in beta_and_ns]
        temp.remove(1)
        return temp


def possible_mu0s(SUK, v):
    r"""
    Return a list `[\mu_0]` of all possible `\mu_0` values defined in [AKMRVW]_.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a finite place of ``K``

    OUTPUT: list ``[mu0s]`` where each ``mu0`` is an element of ``K``

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import possible_mu0s
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: S = tuple(K.primes_above(3))
        sage: SUK = UnitGroup(K, S=S)
        sage: v_fin = S[0]

        sage: possible_mu0s(SUK, v_fin)
        [-1, 1]

    .. NOTE::

        `n_0` is the valuation of the coefficient `\alpha_d` of the `S`-unit equation such that `|\alpha_d \tau_d|_v = 1`
        We have set `n_0 = 0` here since the coefficients are roots of unity
        `\alpha_0` is not defined in the paper, we set it to be 1

    REFERENCES:

    - [AKMRVW]_
    - [Sma1995]_ pp. 824-825, but we modify the definition of ``sigma`` (``sigma_tilde``) to make it easier to code
    """
    beta_and_ns = [[beta,beta.valuation(v)] for beta in SUK.fundamental_units()]
    betak, nk = beta_k(beta_and_ns)
    ns = [beta[1] for beta in beta_and_ns if beta[0] != betak]
    betas = [beta[0] for beta in beta_and_ns if beta[0] != betak]
    mu0s = []
    for rs in combinations_with_replacement(range(abs(nk)), len(betas)):
        # n_0 = valuation_v of one of the coefficients of the equation = 0 for x + y = 1 p. 824
        n_rs = zip(ns, rs)
        sigma_tilde = -(sum([n_r[0]*n_r[1] for n_r in n_rs]))
        if sigma_tilde % nk == 0:
            beta_rs = zip(betas, rs)
            temp_prod = prod([beta_r[0]**beta_r[1] for beta_r in beta_rs]) * betak**(sigma_tilde/nk)
            for alpha0 in SUK.roots_of_unity():
                if alpha0*temp_prod not in mu0s:
                    mu0s.append(alpha0*temp_prod)
    return mu0s


def Yu_a1_kappa1_c1(p, dK, ep):
    r"""
    Compute the constants a(1), kappa1, and c(1) of [Yu2007]_.

    INPUT:

    - ``p`` -- a rational prime number
    - ``dK`` -- the absolute degree of some number field `K`
    - ``ep`` -- the absolute ramification index of some prime ``frak_p`` of `K` lying above `p`

    OUTPUT:

    The constants a(1), kappa1, and c(1).

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import Yu_a1_kappa1_c1
        sage: Yu_a1_kappa1_c1(5, 10, 3)
        (16, 20, 319)

    REFERENCES:

    - [Yu2007]_
    """

    # For readability, we compute a(1) and kappa1 first.

    if p == 2:
        a1 = 32
        kappa1 = 40
    elif p == 3:
        a1 = 16
        kappa1 = 20
    else:
        if ep >= 2:
            a1 = 16
            kappa1 = 20
        else:
            a1 = 8*(p-1)/(p-2)
            kappa1 = 10

    # Next we compute c(1), which has more cases to consider.

    if p == 2:
        c1 = 160
    elif p == 3:
        if dK == 1:
            c1 = 537
        else:
            c1 = 759
    elif p == 5:
        if ep == 1:
            c1 = 1473
        else:
            c1 = 319
    elif p % 4 == 1:
        if ep == 1:
            c1 = 1473
        else:
            c1 = 1502
    else:
        # p > 5 and p % 4 == 3
        if ep == 1:
            if dK == 1:
                c1 = 1288
            else:
                c1 = 1282
        else:
            c1 = 2190

    return a1, kappa1, c1


def Yu_condition_115(K, v):
    r"""
    Return ``True`` or ``False``, as the number field ``K`` and the finite place ``v`` satisfy condition (1.15) of [Yu2007]_.

    INPUT:

    - ``K`` -- a number field
    - ``v`` -- a finite place of ``K``

    OUTPUT:

    ``True`` if (1.15) is satisfied, otherwise ``False``.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import Yu_condition_115
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 + 5)
        sage: v2 = K.primes_above(2)[0]
        sage: v11 = K.primes_above(11)[0]
        sage: Yu_condition_115(K, v2)
        False
        sage: Yu_condition_115(K, v11)
        True

    REFERENCES:

    - [Yu2007]_ p. 188
    """

    p = v.smallest_integer()
    f = v.residue_class_degree()
    w = K.number_of_roots_of_unity()

    # Determine q.

    if p == 2:
        q = 3
    else:
        q = 2

    # Check the condition.

    if q == 2:
        if p**f % 4 == 1:
            return True
        if w % 4 == 0:
            return True
    else:
        if w % 3 == 0:
            return True

    return False


def Yu_modified_height(mu, n, v, prec=106):
    r"""
    Return the value of h(n)(mu) as appearing in [Yu2007]_ equation (1.21).

    INPUT:

    - ``mu`` -- an element of a field K
    - ``n`` -- number of mu_j to be considered in Yu's Theorem
    - ``v`` -- a place of K
    - ``prec`` -- the precision of the real field

    OUTPUT:

    The value `h_p(mu)`.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 + 5)
        sage: v11 = K.primes_above(11)[0]
        sage: from sage.rings.number_field.S_unit_solver import Yu_modified_height
        sage: Yu_modified_height(a, 3, v11)
        0.8047189562170501873003796666131

    If mu is a root of unity, the output is not zero. ::
        sage: Yu_modified_height(-1, 3, v11)
        0.03425564675426243634374205111379

    REFERENCES:

    - [Yu2007]_ p. 192
    """

    R = RealField(prec)

    K = v.number_field()
    dK = K.degree()
    p = v.smallest_integer()
    ep = v.ramification_index()
    fp = v.residue_class_degree()

    a1, kappa1, c1 = Yu_a1_kappa1_c1(p, dK, ep)

    h0 = mu.global_height(prec)
    h1 = R( fp * R(p).log() / (kappa1 * (n + 4) * dK) )

    if h0 > h1:
        return h0
    else:
        return h1


def Omega_prime(dK, v, mu_list, prec=106):
    r"""
    Return the constant `\Omega'` appearing in [AKMRVW]_.

    INPUT:

    - ``dK`` -- the degree of a number field `K`
    - ``v`` -- a finite place of `K`
    - ``mu_list`` -- list of nonzero elements of `K`. It is assumed that the
      sublist ``mu_list[1:]`` is multiplicatively independent
    - ``prec`` -- the precision of the real field

    OUTPUT: the constant `\Omega'`

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import mus, Omega_prime
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(6)))
        sage: v = K.primes_above(3)[0]
        sage: mu_list = [-1] + mus(SUK, v)
        sage: dK = K.degree()
        sage: Omega_prime(dK, v, mu_list)
        0.000487349679922696

    REFERENCES:

    - [AKMRVW]_ :arxiv:`1903.00977`
    """

    R = RealField(prec)

    omega_prime = R(1)
    for mu in mu_list[1:]:
        omega_prime *= mu.global_height()

    n = len(mu_list)

    omega_prime *= Yu_modified_height(mu_list[0], n, v, prec)

    return omega_prime


def Yu_C1_star(n, v, prec=106):
    r"""
    Return the constant `C_1^*` appearing in [Yu2007]_ (1.23).

    INPUT:

    - ``n`` -- the number of generators of a multiplicative subgroup of a field `K`
    - ``v`` -- a finite place of `K` (a fractional ideal)
    - ``prec`` -- the precision of the real field

    OUTPUT: the constant `C_1^*` as a real number

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 + 5)
        sage: v11 = K.primes_above(11)[0]
        sage: from sage.rings.number_field.S_unit_solver import Yu_C1_star
        sage: Yu_C1_star(1,v11)
        2.154667761574516556114215527020e6

    REFERENCES:

    - [Yu2007]_ p.189,193
    """

    R = RealField(prec)

    K = v.number_field()
    dK = K.absolute_degree()

    p = v.smallest_integer()
    ep = v.ramification_index()
    fp = v.residue_class_degree()

    if p == 2:
        q = 3
    else:
        q = 2

    w = K.number_of_roots_of_unity()
    u = ZZ(w).valuation(q)

    a_paren_1, kappa1, c_paren_1 = Yu_a1_kappa1_c1(p, dK, ep)

    C1 = R(1)
    C1 *= c_paren_1
    C1 *= a_paren_1**n
    C1 *= (n**n * (n+1)**(n+1))/factorial(n)
    C1 *= p**fp/(q**u)
    C1 *= ( dK / (fp * R(p).log()) )**(n+2)
    C1 *= R(max( dK, exp(1) )).log()
    C1 *= max( R(exp(4)*(n+1)*dK).log(), ep, fp * R(p).log() )

    C1_star = R((n+1) * C1)

    return C1_star


def Yu_bound(SUK, v, prec=106):
    r"""
    Return `c_8` such that `c_8 \geq exp(2)/\log(2)` and `ord_p (\Theta - 1) < c_8 \log B`,
    where `\Theta = \prod_{j=1}^n \alpha_j^{b_j}` and `B \geq \max_j |b_j|` and `B \geq 3`.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a finite place of `K` (a fractional ideal)
    - ``prec`` -- the precision of the real field

    OUTPUT: the constant `c_8` as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import Yu_bound
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 + 11)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(6)))
        sage: v = K.primes_above(3)[0]
        sage: Yu_bound(SUK, v)
        9.03984381033128e9

    REFERENCES:

    - [Sma1995]_ p. 825
    - [Yu2007]_ p. 189--193 esp. Theorem 1
    - [AKMRVW]_ :arxiv:`1903.00977`
    """

    # We are using Theorem 1 of "p-adic logarithmic forms and group varieties III" by Kunrui Yu.

    # We require the assumption of (1.18): B \geq max {|b_1|,...,|b_n|,3}

    # To be sure that the Lemma of Petho-de Weger is applicable in a later function, we always return a value >= exp(2)/log(2).

    R = RealField(prec)
    p = v.smallest_integer()

    K = SUK.number_field()
    dK = K.absolute_degree()

    mu_free_gens = mus(SUK, v)
    poss_mu0 = possible_mu0s(SUK, v)
    n = 1 + len(mu_free_gens)

    if Yu_condition_115(K,v):
        largest_Omega_prime = R(0)
        for mu0 in poss_mu0:
            current_Omega_prime = Omega_prime(dK, v, [mu0] + mu_free_gens[:], prec)
            largest_Omega_prime = max( current_Omega_prime, largest_Omega_prime )
        C1star = Yu_C1_star(n, v, prec)
        return max( exp(R(2))/R(2).log(), largest_Omega_prime * C1star)
    else:
        # K and v don't satisfy the theorem hypotheses, and we must move to a quadratic extension L.
        # For justification of this next bound, see [AKMRVW].
        x = SR.var('x')
        if p == 2:
            L_over_K = K.extension(x**2 + x + 1, 'xi0')
        else:
            L_over_K = K.extension(x**2 + 1, 'xi0')

        # pick any prime vL over v
        vL_0 = L_over_K.primes_above(v)[0]
        e_vL_v = vL_0.relative_ramification_index()

        # build absolute versions of L and vL

        L = L_over_K.absolute_field('xi_L')
        vL_gens = tuple( [L(z) for z in vL_0.gens()] )
        vL = L.fractional_ideal( vL_gens )

        dL = L.degree()

        largest_Omega_prime = R(0)
        for mu0 in poss_mu0:
            current_Omega_prime = Omega_prime(dL, vL, [mu0] + mu_free_gens[:], prec)
            largest_Omega_prime = max( current_Omega_prime, largest_Omega_prime )
        C1star = Yu_C1_star(n, vL, prec)
        return max(exp(R(2))/R(2).log(), e_vL_v * largest_Omega_prime * C1star)


def K0_func(SUK, A, prec=106):
    r"""
    Return the constant `K_0` from [AKMRVW]_.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``A`` -- the set of the products of the coefficients of the `S`-unit equation with each root of unity of `K`
    - ``prec`` -- the precision of the real field (default: 106)

    OUTPUT: the constant `K_0`, a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import K0_func
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 + 11)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(6)))
        sage: v = K.primes_above(3)[0]
        sage: K0_func(SUK, K.roots_of_unity())
        8.84763586062272e12

    REFERENCES:

    - [Sma1995]_ p. 824
    - [AKMRVW]_ :arxiv:`1903.00977`
    """
    R = RealField(prec)

    K0 = R(1)

    c3 = c3_func(SUK, prec)

    for v_l in SUK.primes():
        e_l = v_l.residue_class_degree()
        Norm_v_l = v_l.absolute_norm()

        c5_l = c3/(e_l * R(Norm_v_l).log())

        c8_l = Yu_bound(SUK, v_l, prec)

        K0_l = (2 * c8_l)/(e_l * c5_l) * R(c8_l / (e_l * c5_l)).log()

        K0 = max(K0, K0_l)

    return K0


def c11_func(SUK, v, A, prec=106):
    r"""
    Return the constant `c_{11}` from Smart's TCDF paper, [Sma1995]_.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a place of `K`, finite (a fractional ideal) or infinite (element of ``SUK.number_field().places(prec)``)
    - ``A`` -- the set of the product of the coefficients of the `S`-unit equation with each root of unity of `K`
    - ``prec`` -- the precision of the real field (default: 106)

    OUTPUT: the constant `c_{11}`, a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c11_func
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c11_func(SUK, phi_real, A)  # abs tol 1e-29
        3.255848343572896153455615423662

        sage: c11_func(SUK, phi_complex, A)  # abs tol 1e-29
        6.511696687145792306911230847323

    REFERENCES:

    - [Sma1995]_ p. 825
    """
    R = RealField(prec)
    if is_real_place(v):
        return R(4*c4_func(SUK, v, A, prec)).log() / c3_func(SUK, prec)
    else:
        return 2*R(4*(c4_func(SUK, v, A, prec)).sqrt()).log() / c3_func(SUK, prec)


def c13_func(SUK, v, prec=106):
    r"""
    Return the constant `c_{13}` from Smart's TCDF paper, [Sma1995]_.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- an infinite place of ``K`` (element of ``SUK.number_field().places(prec)``)
    - ``prec`` -- the precision of the real field (default: 106)

    OUTPUT: the constant `c_{13}`, as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c13_func
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]

        sage: c13_func(SUK, phi_real)  # abs tol 1e-29
        0.4257859134798034746197327286726

        sage: c13_func(SUK, phi_complex)  # abs tol 1e-29
        0.2128929567399017373098663643363

    It is an error to input a finite place. ::

        sage: phi_finite = K.primes_above(3)[0]
        sage: c13_func(SUK, phi_finite)
        Traceback (most recent call last):
        ...
        TypeError: Place must be infinite


    REFERENCES:

    - [Sma1995]_ p. 825
    """
    try:
        v.codomain()
    except AttributeError:
        raise TypeError('Place must be infinite')
    if is_real_place(v):
        return c3_func(SUK, prec)
    else:
        return c3_func(SUK, prec)/2


def K1_func(SUK, v, A, prec=106):
    r"""
    Return the constant `K_1` from Smart's TCDF paper, [Sma1995]_.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- an infinite place of `K` (element of ``SUK.number_field().places(prec)``)
    - ``A`` -- list of all products of each potential `a`, `b` in the `S`-unit
      equation `ax + by + 1 = 0` with each root of unity of `K`
    - ``prec`` -- the precision of the real field (default: 106)

    OUTPUT: the constant `K_1`, a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import K1_func
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: K1_func(SUK, phi_real, A)
        4.483038368145048508970350163578e16

        sage: K1_func(SUK, phi_complex, A)
        2.073346189067285101984136298965e17

    REFERENCES:

    - [Sma1995]_ p. 825
    """
    R = RealField(prec)

    # [Sma1995]_ p. 825
    if is_real_place(v):
        c11 = R(4*c4_func(SUK, v, A, prec)).log() / c3_func(SUK, prec)
    else:
        c11 = 2*( R(4*(c4_func(SUK,v, A, prec)).sqrt()).log() ) / c3_func(SUK, prec)

    # [Sma1995]_ p. 825
    if is_real_place(v):
        c12 = R(2 * c4_func(SUK, v, A, prec))
    else:
        c12 = R(2 * c4_func(SUK, v, A, prec).sqrt())

    # [Sma1998]_ p. 225, Theorem A.1
    d = SUK.number_field().degree()
    t = SUK.rank()
    Baker_C = R(18 * factorial(t+2) * (t+1)**(t+2) * (32*d)**(t+3) * R(2*(t+1) * d).log())

    def hprime(SUK, alpha, v):
        # [Sma1998]_ p. 225
        return R(max(alpha.global_height(), 1/SUK.number_field().degree(), abs( v(alpha).log() ) / SUK.number_field().degree()))

    # [Sma1995]_ p. 825 and [Sma1998]_ p. 225, Theorem A.1
    c14 = Baker_C * prod([hprime(SUK, alpha, v) for alpha in SUK.gens_values()])

    # [Sma1995]_ p. 825
    c13 = c13_func(SUK,v,prec)
    w = len(SUK.roots_of_unity())
    c15 = (2/c13)*(c12.log()+c14*(((t+1)*w*c14/c13).log()))

    return max([c11, c15])


def minimal_vector(A, y, prec=106):
    r"""
    INPUT:

    - ``A`` -- a square `n` by `n` non-singular integer matrix whose rows generate a lattice `\mathcal L`
    - ``y`` -- a row (1 by `n`) vector with integer coordinates
    - ``prec`` -- precision of real field (default: 106)

    OUTPUT: a lower bound for the square of

    .. MATH::

        \ell (\mathcal L,\vec y) =
        \begin{cases}
        \displaystyle\min_{\vec x\in\mathcal L}\Vert\vec x-\vec y\Vert &, \vec y\not\in\mathcal L. \\
        \displaystyle\min_{0\neq\vec x\in\mathcal L}\Vert\vec x\Vert &,\vec y\in\mathcal L.
        \end{cases}`

    ALGORITHM:

    The algorithm is based on V.9 and V.10 of [Sma1998]_

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import minimal_vector
        sage: B = matrix(ZZ, 2, [1,1,1,0])
        sage: y = vector(ZZ, [2,1])
        sage: minimal_vector(B, y)
        1/2

    ::

        sage: B = random_matrix(ZZ, 3)
        sage: while not B.determinant():
        ....:     B = random_matrix(ZZ, 3)
        sage: B # random
        [-2 -1 -1]
        [ 1  1 -2]
        [ 6  1 -1]
        sage: y = vector([1, 2, 100])
        sage: minimal_vector(B, y)  # random
        15/28
    """
    if A.is_singular():
        raise ValueError('The matrix A is singular')

    R = RealField(prec)

    n = len(y)
    c1 = 2**(n-1)
    ALLL = A.LLL()
    ALLLinv = ALLL.inverse()
    ybrace = [ abs(R(a-a.round())) for a in y * ALLLinv if (a-a.round()) != 0]

    v = ALLL.rows()[0]
    if len(ybrace) == 0:
        return v.dot_product(v) / c1
    else:
        sigma = ybrace[len(ybrace)-1]
        return v.dot_product(v) * sigma / c1


def reduction_step_complex_case(place, B0, list_of_gens, torsion_gen, c13):
    r"""
    INPUT:

    - ``place`` -- (ring morphism) an infinite place of a number field `K`
    - ``B0`` -- the initial bound
    - ``list_of_gens`` -- set of generators of the free part of the group
    - ``torsion_gen`` -- an element of the torsion part of the group
    - ``c13`` -- a positive real number

    OUTPUT: a tuple consisting of:

    1. a new upper bound, an integer
    2. a boolean value, ``True`` if we have to increase precision, otherwise ``False``

    .. NOTE::

        The constant ``c13``  in Section 5, [AKMRVW]_
        This function does handle both real and non-real infinite places.

    REFERENCES:

    See [Sma1998]_, [AKMRVW]_.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import reduction_step_complex_case
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField([x^3 - 2])
        sage: SK = sum([K.primes_above(p) for p in [2,3,5]],[])
        sage: G = [g for g in K.S_unit_group(S=SK).gens_values()
        ....:        if g.multiplicative_order() == Infinity]
        sage: p1 = K.places(prec=100)[1]
        sage: reduction_step_complex_case(p1, 10^5, G, -1, 2)
        (18, False)
    """
    prec = place.codomain().precision()
    R = RealField(prec)
    CF = ComplexField(prec)
    n = len(list_of_gens)
    w = torsion_gen.multiplicative_order()

    real_part_log_gens = [ R(CF(place(g).log()).real_part()) for g in list_of_gens]
    imag_part_log_gens = [ R(CF(place(g).log()).imag_part()) for g in list_of_gens]
    real_part_log_gens += [R(0)]
    imag_part_log_gens += [2*R.pi()/w]

    abs_log_parts = [abs(part) for part in real_part_log_gens]+[abs(part) for part in imag_part_log_gens]
    max_part_log = max(abs_log_parts)

    npi = []
    # we collect the list of indices of log(g) which are not pure imaginary
    # if this list is empty, we have to take a special case
    for i in range(len(real_part_log_gens)):
        lg = real_part_log_gens[i]
        if abs(lg) > 2**(-place.codomain().precision()):
            npi.append(i)
    # someday make this a separate function
    if not npi:
        # this is the pure imaginary case.
        # we have only imaginary numbers

        C = ZZ(1)
        S = n * B0**2
        T = (n+w+n*w)*B0 / 2
        finish = False
        while not finish:
            A = identity_matrix(ZZ, n+1)
            A[n] = vector([(g * C).round() for g in imag_part_log_gens])

            if A.is_singular():
                C = ZZ(2*C)
            else:
                # We have to work with rows because of the .LLL() function

                A = A.transpose()

                # Note that l is the lower bound on the square of the magnitude of the shortest nonzero vector in the lattice generated by A
                l = minimal_vector(A, zero_vector(ZZ,n+1))
                # Checking hypotheses of Lemma 5.3 in our paper:

                if l <= T**2+S:
                    C = ZZ(2*C)
                    # Need to check precision: must be at least two more than the number of digits in largest entry in A to ensure that we get true rounding--
                    if prec < R(C*max_part_log).log()/R(2).log()+3:
                        return 0, True
                else:
                    # Need to check precision: must be at least two more than the number of digits in largest entry in A to ensure that we get true rounding--
                    if prec < R(C*max_part_log).log()/R(2).log()+3:
                        return 0, True
                    else:
                        Bnew = ((R(C * 2).log() - ((l**2-S).sqrt()-T)).log() / c13).round()
                        finish = True
                        return max(4,w,Bnew), False
    elif is_real_place(place):
        # this is the case when we are working with a real embedding, we get savings here
        C = R(1)
        S = (n-1) * B0**2
        w = place.domain().number_of_roots_of_unity()
        T = (n*B0+1)/R(2)
        finish = False

        while not finish:

            A = copy(identity_matrix(ZZ, n+1))
            # We redefine the imaginary parts in case any generator was negative
            new_imag_part_log_gens = [0 for i in imag_part_log_gens[:-1]]+[imag_part_log_gens[-1]]
            A[n-1] = vector([(g*C).round() for g in real_part_log_gens])
            A[n] = vector([(g*C).round() for g in new_imag_part_log_gens])

            if A.is_singular():
                C *= 2
            else:
                # We apply Lemma 5.3 from [AKMRVW]
                A = A.transpose()
                l = minimal_vector(A, zero_vector(ZZ,n+1))
                # Note that l is the a lower bound on the square of the magnitude of the shortest nonzero vector in the lattice generated by A
                # Checking hypothesis of lemma 5.3 in [AKMRVW]
                if l <= T**2 + S:
                    C *= 2
                    # Need to check precision: must be at least two more than the number of digits in largest entry in A to ensure that we get true rounding--
                    if prec < R(C*max_part_log).log()/R(2).log()+3:
                        return 0, True
                else:
                    # Need to check precision: must be at least two more than the number of digits in largest entry in A to ensure that we get true rounding--
                    if prec < R(C*max_part_log).log()/R(2).log()+3:
                        return 0, True
                    else:
                        Bnew = ((R(C * 2).log() - ((l-S).sqrt()-T).log()) / c13).round()
                        finish = True
                        return max(4,w,Bnew), False

    else:

        # the case when the real part is not 0 for all log(a_i), see Lemma 5.2 in [AKMRVW]
        C = R(1)
        S = (n-1) * B0**2
        w = place.domain().number_of_roots_of_unity()
        T = (n+w+n*w)*B0/R(2).sqrt()
        finish = False

        # we reorder the generators to that the real part of the last non-torsion generator is not 0:
        if n-1 not in npi:
            new_last_gen_index = npi[0]
            old_last_gen_real = real_part_log_gens[n-1]
            old_last_gen_imag = imag_part_log_gens[n-1]
            real_part_log_gens[n-1] = real_part_log_gens[new_last_gen_index]
            imag_part_log_gens[n-1] = imag_part_log_gens[new_last_gen_index]
            real_part_log_gens[new_last_gen_index] = old_last_gen_real
            imag_part_log_gens[new_last_gen_index] = old_last_gen_imag

        while not finish:

            A = copy(identity_matrix(ZZ, n+1))
            A[n-1] = vector([(g*C).round() for g in real_part_log_gens])
            A[n] = vector([(g*C).round() for g in imag_part_log_gens])

            if A.is_singular():
                C *= 2
            else:
                # We apply Lemma 5.2 from [AKMRVW]
                A = A.transpose()
                l = minimal_vector(A, zero_vector(ZZ,n+1))
                # Note that l is the a lower bound on the square of the magnitude of the shortest nonzero vector in the lattice generated by A
                # Checking hypothesis of lemma 5.2 in [AKMRVW]
                if l <= T**2 + S:
                    C *= 2
                    # Need to check precision: must be at least two more than the number of digits in largest entry in A to ensure that we get true rounding--
                    if prec < R(C*max_part_log).log()/R(2).log()+3:
                        return 0, True
                else:
                    # Need to check precision: must be at least two more than the number of digits in largest entry in A to ensure that we get true rounding--
                    if prec < R(C*max_part_log).log()/R(2).log()+3:
                        return 0, True
                    else:
                        Bnew = ((R(C * 2).log() - ((l-S).sqrt()-T).log()) / c13).round()
                        finish = True
                        return max(4,w,Bnew), False


def cx_LLL_bound(SUK, A, prec=106):
    r"""
    Return the maximum of all of the `K_1`'s as they are LLL-optimized for each infinite place `v`.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``A`` -- list of all products of each potential `a`, `b` in the `S`-unit
      equation `ax + by + 1 = 0` with each root of unity of `K`
    - ``prec`` -- precision of real field (default: 106)

    OUTPUT: a bound for the exponents at the infinite place, as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import cx_LLL_bound
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: A = K.roots_of_unity()

        sage: cx_LLL_bound(SUK, A) # long time
        35
    """
    cx_LLL = 0
    # initialize a bound, a bad guess, as we iterate over the places of the number field, we will replace its value with the largest complex LLL bound we've found across the places
    for v in SUK.number_field().places(prec=prec):
        prec_v = prec
        c13_LLL = c13_func(SUK, v, prec_v)
        cx_bound = K1_func(SUK, v, A, prec_v)
        # cx_bound is the LLL bound according to this place, it will be replaced as LLL gives us smaller bounds
        new_bound, inc_prec = reduction_step_complex_case(v, cx_bound, SUK.fundamental_units(), SUK.zeta(), c13_LLL)
        while inc_prec:
            v = refine_embedding(v)
            c13_LLL = c13_func(SUK, v, prec_v)
            cx_bound = K1_func(SUK, v, A, prec_v)
            new_bound, inc_prec = reduction_step_complex_case(v, cx_bound, SUK.fundamental_units(), SUK.zeta(), c13_LLL)
        counter = 0
        while abs(cx_bound - new_bound) > .5*cx_bound and counter < 15:
            # We fear a loop that is not convergent, this is the purpose of the counter
            # Repeat complex LLL until we get essentially no change from it
            cx_bound = min(cx_bound, new_bound)
            new_bound, inc_prec = reduction_step_complex_case(v, cx_bound, SUK.fundamental_units(), SUK.zeta(), c13_LLL)
            while inc_prec:
                v = refine_embedding(v)
                c13_LLL = c13_func(SUK, v, prec_v)
                new_bound, inc_prec = reduction_step_complex_case(v, cx_bound, SUK.fundamental_units(), SUK.zeta(), c13_LLL)
            counter += 1

        cx_bound = min(cx_bound, new_bound)
        # for this place the complex LLL bound is cx_bound
        cx_LLL = max(cx_bound, cx_LLL)
        # compare this value with the complex LLL bounds we have found for the previous places, if it is bigger, replace that bound
    return cx_LLL


def log_p(a, prime, prec):
    r"""
    INPUT:

    - ``a`` -- an element of a number field `K`
    - ``prime`` -- a prime ideal of the number field `K`
    - ``prec`` -- positive integer

    OUTPUT:

    An element of `K` which is congruent to the ``prime``-adic logarithm of `a`
    with respect to ``prime`` modulo ``p^prec``, where `p` is the rational
    prime below ``prime``.

    .. NOTE::

        Here we take into account the other primes in `K` above `p` in order to
        get coefficients with small values.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import log_p
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 + 14)
        sage: p1 = K.primes_above(3)[0]
        sage: p1
        Fractional ideal (3, a + 1)
        sage: log_p(a+2, p1, 20)
        8255385638/3*a + 15567609440/3

    ::

        sage: K.<a> = NumberField(x^4 + 14)
        sage: p1 = K.primes_above(5)[0]
        sage: p1
        Fractional ideal (5, a + 1)
        sage: log_p(1/(a^2-4), p1, 30)
        -42392683853751591352946/25*a^3 - 113099841599709611260219/25*a^2 -
        8496494127064033599196/5*a - 18774052619501226990432/25
    """
    if a == 0:
        raise ValueError('a is the zero element')

    if a.valuation(prime) != 0:
        raise ValueError('The valuation of a with respect to prime is not zero')

    K = prime.ring()
    p = prime.smallest_integer()

    # In order to get an approximation with small coefficients we have to take into account the other primes above p
    # with negative valuation.  For example, say prime2 is another (principal ideal) prime above p, and a=(unit)(prime2)^(-k) for some unit and k
    # a positive integer, and let tilde(a):=a(prime2)^k.  Then log_p(a)=log_p(tilde(a))-k(log_p(prime2)), where the series representations
    # of these two logs will have smaller coefficients.

    primes = [(-(a.valuation(pr)),pr) for pr in K.primes_above(p) if a.valuation(pr) < 0]
    local_terms = []

    for (val, pr) in primes:
        # for its pair in primes we find an element in K such that it is divisible only by pr and not by any other ideal above p. Then we take this element in the correct exponent

        if pr.is_principal():
            local_terms.append(pr.gens_reduced()[0]**val)
        else:
            local_terms.append(pr.gens()[1]**val)

    return log_p_series_part(a*prod(local_terms), prime, prec) - sum([log_p_series_part(b, prime, prec) for b in local_terms])


def log_p_series_part(a, prime, prec):
    r"""
    INPUT:

    - ``a`` -- an element of a number field `K`
    - ``prime`` -- a prime ideal of the number field `K`
    - ``prec`` -- positive integer

    OUTPUT:

    The ``prime``-adic logarithm of `a` and accuracy ``p^prec``, where `p` is
    the rational prime below ``prime``.

    ALGORITHM:

    The algorithm is based on the algorithm on page 30 of [Sma1998]_.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import log_p_series_part
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 5)
        sage: p1 = K.primes_above(3)[0]
        sage: p1
        Fractional ideal (3)
        sage: log_p_series_part(a^2 - a + 1, p1, 30)
        120042736778562*a + 263389019530092

    ::

        sage: K.<a> = NumberField(x^4 + 14)
        sage: p1 = K.primes_above(5)[0]
        sage: p1
        Fractional ideal (5, a + 1)
        sage: log_p_series_part(1/(a^2-4), p1, 30)
        5628940883264585369224688048459896543498793204839654215019548600621221950915106576555819252366183605504671859902129729380543157757424169844382836287443485157589362653561119898762509175000557196963413830027960725069496503331353532893643983455103456070939403472988282153160667807627271637196608813155377280943180966078/1846595723557147156151786152499366687569722744011302407020455809280594038056223852568951718462474153951672335866715654153523843955513167531739386582686114545823305161128297234887329119860255600972561534713008376312342295724191173957260256352612807316114669486939448006523889489471912384033203125*a^2 + 2351432413692022254066438266577100183514828004415905040437326602004946930635942233146528817325416948515797296867947688356616798913401046136899081536181084767344346480810627200495531180794326634382675252631839139904967037478184840941275812058242995052383261849064340050686841429735092777331963400618255005895650200107/1846595723557147156151786152499366687569722744011302407020455809280594038056223852568951718462474153951672335866715654153523843955513167531739386582686114545823305161128297234887329119860255600972561534713008376312342295724191173957260256352612807316114669486939448006523889489471912384033203125
    """
    if a.valuation(prime) != 0:
        raise ValueError('The valuation of a with respect to prime is not zero')
    K = prime.ring()
    g = K.gen()
    p = prime.smallest_integer()
    f = prime.residue_class_degree()
    e = prime.absolute_ramification_index()
    q = p**f - 1
    R = RealField(prec)

    divisor = q.divisors()
    order = min(d for d in divisor if (a**d - 1).valuation(prime) > 0)
    gamma = a**order
    t = 0
    while (gamma-1).valuation(prime) <= e:
        t += 1
        gamma = gamma**p
    prec += t
    # since later we divide by p^t, we must increase the precision by t at this point.
    m = (gamma-1).valuation(prime) / e
    n = Integer(1)
    step = 10 ** (R(prec).log()/R(10).log()).floor()
    while n < (R(n).log()/R(p).log() + prec)/m:
        n += step
    # could use smaller stepsize to get actual smallest integer n, however this seems to run faster.
    w = (R(prec).log()/R(p).log()).floor()
    gamma = sum([ZZ(gi % (p**(prec+w))) * g**i
                 if gi.valuation(p) >= 0 else
                 ZZ((gi * p**(-gi.valuation(p))) % (p**(prec+w-gi.valuation(p)))) * p**(gi.valuation(p)) * g**i
                 for i,gi in enumerate(gamma) if gi != 0])

    beta = 0
    delta = 1 - gamma
    for i in range(1, n+1):
        beta -= delta / i
        delta *= (1 - gamma)
        delta = sum([ZZ(di % (p**(prec+w))) * g**b
                     if di.valuation(p) >= 0 else
                     ZZ((di * p**(-di.valuation(p))) % (p**(prec + w - di.valuation(p)))) * p**(di.valuation(p)) * g**b
                     for b,di in enumerate(delta) if di != 0])
    beta = beta / (order * p**t)

    # we try to make the coefficients small

    logp = 0
    for i,b in enumerate(beta.list()):
        val = b.valuation(p)
        if val < 0:
            t = b * p**(-val)
            t = ZZ(mod(t, p**(prec-val)))
            t = t * p**val
        else:
            t = ZZ(mod(b, p**prec))
        logp = logp + t * g**i

    return logp


def defining_polynomial_for_Kp(prime, prec=106):
    r"""
    INPUT:

    - ``prime`` -- a prime ideal of a number field `K`
    - ``prec`` -- a positive natural number (default: 106)

    OUTPUT: a polynomial with integer coefficients that is equivalent ``mod p^prec`` to a defining polynomial for the completion of `K` associated to the specified prime

    .. NOTE::

        `K` has to be an absolute extension

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import defining_polynomial_for_Kp
        sage: K.<a> = QuadraticField(2)
        sage: p2 = K.prime_above(7); p2
        Fractional ideal (-2*a + 1)
        sage: defining_polynomial_for_Kp(p2, 10)
        x + 266983762

    ::

        sage: K.<a> = QuadraticField(-6)
        sage: p2 = K.prime_above(2); p2
        Fractional ideal (2, a)
        sage: defining_polynomial_for_Kp(p2, 100)
        x^2 + 6
        sage: p5 = K.prime_above(5); p5
        Fractional ideal (5, a + 2)
        sage: defining_polynomial_for_Kp(p5, 100)
        x + 3408332191958133385114942613351834100964285496304040728906961917542037
    """
    K = prime.ring()
    if not K.is_absolute():
        raise ValueError('The number field is not an absolute extension')

    theta = K.gen()
    f = K.defining_polynomial()
    p = prime.smallest_integer()
    e = prime.absolute_ramification_index()

    N = prec
    while True:
        RQp = Qp(p, prec=N, type='capped-rel', print_mode='series')

        # We factor f in Integers(p**(precision)) using the factorization in Qp

        factors = f.change_ring(RQp).factor()

        # We are going to find which factor of f is related to the prime ideal 'prime'

        L = [g.change_ring(ZZ) for g, _ in factors]
        A = [g for g in L if (g(theta)).valuation(prime) >= e*N/2]

        # We narrow down the list unitl only one value remains

        if len(A) == 1:
            return A[0].change_ring(Integers(p**prec)).change_ring(ZZ)
        else:
            N += 1


def embedding_to_Kp(a, prime, prec):
    r"""
    INPUT:

    - ``a`` -- an element of a number field `K`
    - ``prime`` -- a prime ideal of `K`
    - ``prec`` -- a positive natural number

    OUTPUT:

    An element of `K` that is equivalent to `a` modulo ``p^(prec)`` and the generator of `K` appears with exponent less than `e \cdot f`, where `p` is the rational prime below ``prime`` and `e`, `f` are the ramification index and residue degree, respectively.

    .. NOTE::

        `K` has to be an absolute number field

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import embedding_to_Kp
        sage: K.<a> = QuadraticField(17)
        sage: p = K.prime_above(13); p
        Fractional ideal (-a + 2)
        sage: embedding_to_Kp(a-3, p, 15)
        -20542890112375827

    ::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^4 - 2)
        sage: p = K.prime_above(7); p
        Fractional ideal (-a^2 + a - 1)
        sage: embedding_to_Kp(a^3 - 3, p, 15)
        -1261985118949117459462968282807202378
    """
    K = prime.ring()
    if not K.is_absolute():
        raise ValueError('K has to be an absolute extension')

    g = defining_polynomial_for_Kp(prime, prec).change_ring(QQ)
    gen = K.gen()
    f = K(a).lift()

    return K( sum([b*gen**j for j,b in enumerate(f.mod(g))]) )


def p_adic_LLL_bound_one_prime(prime, B0, M, M_logp, m0, c3, prec=106):
    r"""
    INPUT:

    - ``prime`` -- a prime ideal of a number field `K`
    - ``B0`` -- the initial bound
    - ``M`` -- list of elements of `K`, the `\mu_i`'s from Lemma IX.3 of [Sma1998]_
    - ``M_logp`` -- the `p`-adic logarithm of elements in `M`
    - ``m0`` -- an element of `K`, this is `\mu_0` from Lemma IX.3 of [Sma1998]_
    - ``c3`` -- a positive real constant
    - ``prec`` -- the precision of the calculations (default: 106), i.e.,
      values are known to ``O(p^prec)``

    OUTPUT: a pair consisting of:

    1. a new upper bound, an integer
    2. a boolean value, ``True`` if we have to increase precision, otherwise ``False``

    .. NOTE::

        The constant `c_5` is the constant `c_5` at the page 89 of [Sma1998]_ which is equal to the constant `c_{10}` at the page 139 of [Sma1995]_.
        In this function, the `c_i` constants are in line with [Sma1998]_, but generally differ from the constants in [Sma1995]_ and other parts of this code.

    EXAMPLES:

    This example indicates a case where we must increase precision::

        sage: from sage.rings.number_field.S_unit_solver import p_adic_LLL_bound_one_prime
        sage: prec = 50
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^3 - 3)
        sage: S = tuple(K.primes_above(3))
        sage: SUK = UnitGroup(K, S=S)
        sage: v = S[0]
        sage: A = SUK.roots_of_unity()
        sage: K0_old = 9.4755766731093e17
        sage: Mus = [a^2 - 2]
        sage: Log_p_Mus = [185056824593551109742400*a^2
        ....:               + 1389583284398773572269676*a + 717897987691852588770249]
        sage: mu0 = K(-1)
        sage: c3_value = 0.42578591347980
        sage: m0_Kv_new, increase_prec = p_adic_LLL_bound_one_prime(v, K0_old, Mus, Log_p_Mus,
        ....:                                                       mu0, c3_value, prec)
        sage: m0_Kv_new
        0
        sage: increase_prec
        True

    And now we increase the precision to make it all work::

        sage: prec = 106
        sage: K0_old = 9.475576673109275443280257946930e17
        sage: Log_p_Mus = [1029563604390986737334686387890424583658678662701816*a^2
        ....:               + 661450700156368458475507052066889190195530948403866*a]
        sage: c3_value = 0.4257859134798034746197327286726
        sage: m0_Kv_new, increase_prec = p_adic_LLL_bound_one_prime(v, K0_old, Mus, Log_p_Mus,
        ....:                                                       mu0, c3_value, prec)
        sage: m0_Kv_new
        476
        sage: increase_prec
        False
    """
    if any(g.valuation(prime) != 0 for g in M+[m0]):
        raise ValueError('There is an element with nonzero valuation')

    K = prime.ring()
    w = K.number_of_roots_of_unity()
    p = prime.smallest_integer()
    f = prime.residue_class_degree()
    e = prime.absolute_ramification_index()
    R = RealField(prec)
    c5 = c3 / (f*e*R(p).log())
    theta = K.gen()

    # if M is empty then it is easy to give an upper bound
    if len(M) == 0:
        if m0 != 1:
            return max(4,w, R(max(R(p).log()*f*(m0-1).valuation(prime)/c3, 0)).floor()), False
        else:
            return 0, False
    # we evaluate the p-adic logarithms of m0 and we embed it in the completion of K with respect to prime

    m0_logp = log_p(m0, prime, prec)
    m0_logp = embedding_to_Kp(m0_logp, prime, prec)
    n = len(M_logp)
    # Below we implement paragraph VI.4.2 of [Sma1998], pages 89-93

    # we evaluate the order of discriminant of theta

    Theta = [theta**i for i in range(K.absolute_degree())]
    ordp_Disc = (K.disc(Theta)).valuation(p)
    # We evaluate Lambda

    c8 = min(min(a.valuation(p) for a in g) for g in M_logp)
    lam = p**c8

    # we apply lemma VI.5 of [Sma1998] page 90
    # c6 is 0 here because we seek to solve the equation x+y=1, so our set A
    # is contained in the roots of unity of K

    # In one very extreme case (p = 2 and all other constants as small as possible),
    # low_bound = 1/c5 is not quite enough to give strict inequality. So we add 1 to be safe.

    low_bound = (1/c5).round() + 1
    for a in m0_logp:
        if a != 0 and c8 > a.valuation(p):
            B1 = (c8 + ordp_Disc/2) / c5
            if B1 > low_bound:
                return max(4,w,RR(B1).floor()), False
            else:
                return max(4,w,low_bound), False

    c8 = min([a.valuation(p) for a in m0_logp] + [c8])
    B = [g/lam for g in M_logp]
    b0 = m0_logp / lam
    c9 = c8 + ordp_Disc/2

    # We evaluate 'u' and we construct the matrix A

    m = e * f
    u = 1
    while True:
        if prec <= u + c8:
            return 0, True

        # We construct the matrix A as a block matrix

        A11 = identity_matrix(ZZ, n)
        A12 = zero_matrix(ZZ, n, m)
        A21 = zero_matrix(ZZ, n, m)
        A22 = p**u * identity_matrix(ZZ, m)
        for i,b in enumerate(B):
            A21[i] = vector([mod(b[j],p**u) for j in range(m)])
        A = block_matrix( [[A11,A12], [A21.transpose(),A22]] )

        y = zero_vector(ZZ, n+m)
        for i in range(m):
            y[i+n] = -mod(b0[i], p**u)
        # This refers to c10 from Smart
        c10squared = minimal_vector(A.transpose(), y)
        if c10squared > n * B0**2:
            B2 = (u+c9) / c5
            if B2 > low_bound:
                return max(4,w,R(B2).floor()),False
            else:
                return max(4,w,low_bound),False
        else:
            u += 1


def p_adic_LLL_bound(SUK, A, prec=106):
    r"""
    Return the maximum of all of the `K_0`'s as they are LLL-optimized for each
    finite place `v`.

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``A`` -- list of all products of each potential `a`, `b` in the `S`-unit
      equation `ax + by + 1 = 0` with each root of unity of `K`
    - ``prec`` -- precision for `p`-adic LLL calculations (default: 106)

    OUTPUT:

    A bound for the max of exponents in the case that extremal place is finite (see [Sma1995]_) as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import p_adic_LLL_bound
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: A = SUK.roots_of_unity()
        sage: prec = 100
        sage: p_adic_LLL_bound(SUK, A, prec)
        89
    """
    S = SUK.primes()
    K0_old = K0_func(SUK, A, prec)
    LLL_K0_by_finite_place = []
    for i,v in enumerate(S):
        # Kv_old = K0_by_finite_place[0]
        Mus0 = possible_mu0s(SUK, v)
        Mus = mus(SUK, v)
        Log_p_Mus = [log_p(a, v, prec) for a in Mus]
        local_prec = prec
        val = 0
        for m0 in Mus0:
            m0_Kv_old = K0_old
            m0_Kv_new, increase_precision = p_adic_LLL_bound_one_prime(v, m0_Kv_old, Mus, Log_p_Mus, m0, c3_func(SUK, local_prec), local_prec)
            while increase_precision:
                local_prec *= 2
                Log_p_Mus = [log_p(a, v, local_prec) for a in Mus]
                Log_p_Mus = [embedding_to_Kp(a, v, local_prec) for a in Log_p_Mus]
                m0_Kv_new, increase_precision = p_adic_LLL_bound_one_prime(v, m0_Kv_old, Mus, Log_p_Mus, m0, c3_func(SUK, local_prec), local_prec)

            while m0_Kv_new < m0_Kv_old:
                m0_Kv_old = m0_Kv_new
                m0_Kv_new, increase_precision = p_adic_LLL_bound_one_prime(v, m0_Kv_old, Mus, Log_p_Mus, m0, c3_func(SUK,local_prec), local_prec)
                while increase_precision:
                    local_prec *= 2
                    Log_p_Mus = [log_p(a, v, local_prec) for a in Mus]
                    Log_p_Mus = [embedding_to_Kp(a, v, local_prec) for a in Log_p_Mus]
                    m0_Kv_new,increase_precision = p_adic_LLL_bound_one_prime(v, m0_Kv_old, Mus, Log_p_Mus, m0, c3_func(SUK, local_prec), local_prec)

            val = max(m0_Kv_old, val)

        LLL_K0_by_finite_place.append(val)
    return max(LLL_K0_by_finite_place)


def split_primes_large_lcm(SUK, bound):
    r"""
    Return a list `L` of rational primes `q` which split completely in `K` and
    which have desirable properties (see NOTE).

    INPUT:

    - ``SUK`` -- the `S`-unit group of an absolute number field `K`
    - ``bound`` -- positive integer

    OUTPUT: list `L` of rational primes `q`, with the following properties:

    - each prime `q` in `L` splits completely in `K`
    - if `Q` is a prime in `S` and `q` is the rational
      prime below `Q`, then `q` is **not** in `L`
    - the value `\lcm(\{ q-1 : q \in L \})` is greater than or equal to ``2*bound + 1``.

    .. NOTE::

        - A series of compatible exponent vectors for the primes in `L` will
          lift to **at most** one integer exponent vector whose entries
          `a_i` satisfy `|a_i|` is less than or equal to ``bound``.

        - The ordering of this set is not very intelligent for the purposes
          of the later sieving processes.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import split_primes_large_lcm
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3*x + 1)
        sage: S = K.primes_above(3)
        sage: SUK = UnitGroup(K, S=tuple(S))
        sage: split_primes_large_lcm(SUK, 200)
        [17, 19, 37, 53]

    With a tiny bound, Sage may ask you to increase the bound.

    ::

        sage: from sage.rings.number_field.S_unit_solver import split_primes_large_lcm
        sage: K.<xi> = NumberField(x^2 + 163)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(23)))
        sage: split_primes_large_lcm(SUK, 8)
        Traceback (most recent call last):
        ...
        ValueError: Not enough split primes found. Increase bound.
    """

    K = SUK.number_field()
    # we recover the rational primes below S:
    S0 = set(prime_ideal.smallest_integer() for prime_ideal in SUK.primes())
    split_prime_list = K.completely_split_primes(4*bound + 4)
    lcm_list = []
    L = 1
    while L < 2*bound + 1:
        if split_prime_list == []:
            # Need More Primes!
            raise ValueError('Not enough split primes found. Increase bound.')
        q = split_prime_list.pop(0)
        # only use q if it is *not* below a prime in S -- that is,
        # only if q does *not* appear in S0.
        if q not in S0:
            L = lcm(L, q-1)
            lcm_list.append(q)
    return lcm_list


def sieve_ordering(SUK, q):
    r"""
    Return ordered data for running sieve on the primes in ``SUK`` over the rational prime `q`.

    INPUT:

    - ``SUK`` -- the `S`-unit group of a number field `K`
    - ``q`` -- a rational prime number which splits completely in `K`

    OUTPUT: list of tuples;
    ``[ideals_over_q, residue_fields, rho_images, product_rho_orders]``, where:

    1. ``ideals_over_q`` is a list of the `d = [K:\QQ]` ideals in `K` over `q`
    2. ``residue_fields[i]`` is the residue field of ``ideals_over_q[i]``
    3. ``rho_images[i]`` is a list of the reductions of the generators in of the `S`-unit group, modulo ``ideals_over_q[i]``
    4. ``product_rho_orders[i]`` is the product of the multiplicative orders of the elements in ``rho_images[i]``

    .. NOTE::

        - The list ``ideals_over_q`` is sorted so that the product of orders is smallest for ``ideals_over_q[0]``, as this will make the later sieving steps more efficient.
        - The primes of `S` must not lie over `q`.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import sieve_ordering
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3*x + 1)
        sage: SUK = K.S_unit_group(S=3)
        sage: sieve_data = list(sieve_ordering(SUK, 19))
        sage: sieve_data[0]
        (Fractional ideal (-2*xi^2 + 3),
         Fractional ideal (-xi + 3),
         Fractional ideal (2*xi + 1))

        sage: sieve_data[1]
        (Residue field of Fractional ideal (-2*xi^2 + 3),
         Residue field of Fractional ideal (-xi + 3),
         Residue field of Fractional ideal (2*xi + 1))

        sage: sieve_data[2]
        ([18, 12, 16, 8], [18, 16, 10, 4], [18, 10, 12, 10])

        sage: sieve_data[3]
        (648, 2916, 3888)
    """

    K = SUK.number_field()
    rho = SUK.gens_values()
    d = K.absolute_degree()
    primes_over_q = K.primes_above(q)
    # q must split completely.
    if len(primes_over_q) != d:
        raise ValueError('The prime q is not completely split.')

    for P in SUK.primes():
        if P in primes_over_q:
            raise ValueError('There is a prime in S over q.')

    q_data = []
    for Qi in primes_over_q:
        resfield = Qi.residue_field()
        rho_mod_Qi = [resfield(rho_j) for rho_j in rho]
        orderprod = prod(rho_ij.multiplicative_order() for rho_ij in rho_mod_Qi)
        q_data.append([Qi, resfield, rho_mod_Qi, orderprod])
    q_data.sort(key=lambda X: [X[3],X[0],X[1],X[2]])
    # zip() will change the list of n list of length m to m tuples of length n
    return zip(*q_data)


def clean_rfv_dict(rfv_dictionary):
    r"""
    Given a residue field vector dictionary, remove some impossible keys and entries.

    INPUT:

    - ``rfv_dictionary`` -- dictionary whose keys are exponent vectors and
      whose values are residue field vectors

    OUTPUT: none, but it removes some keys from the input dictionary

    .. NOTE::

        - The keys of a residue field vector dictionary are exponent vectors modulo `q-1` for some prime `q`.
        - The values are residue field vectors. It is known that a residue field vector
          which comes from a solution to the `S`-unit equation cannot have 1 in any entry.

    EXAMPLES:

    In this example, we use a truncated list generated when solving the `S`-unit equation in the case that `K` is defined by the
    polynomial `x^2+x+1` and `S` consists of the primes above 3::

        sage: from sage.rings.number_field.S_unit_solver import clean_rfv_dict
        sage: rfv_dict = {(1, 3): [3, 2], (3, 0): [6, 6], (5, 4): [3, 6], (2, 1): [4, 6],
        ....:             (5, 1): [3, 1], (2, 5): [1, 5], (0, 3): [1, 6]}
        sage: len(rfv_dict)
        7
        sage: clean_rfv_dict(rfv_dict)
        sage: len(rfv_dict)
        4
        sage: rfv_dict
        {(1, 3): [3, 2], (2, 1): [4, 6], (3, 0): [6, 6], (5, 4): [3, 6]}
    """

    for a, val in list(rfv_dictionary.items()):
        if 1 in val:
            rfv_dictionary.pop(a)


def construct_rfv_to_ev(rfv_dictionary, q, d, verbose=False):
    r"""
    Return a reverse lookup dictionary, to find the exponent vectors associated
    to a given residue field vector.

    INPUT:

    - ``rfv_dictionary`` -- dictionary whose keys are exponent vectors and
      whose values are the associated residue field vectors
    - ``q`` -- a prime (assumed to split completely in the relevant number field)
    - ``d`` -- the number of primes in `K` above the rational prime `q`
    - ``verbose`` -- boolean (default: ``False``); flag to indicate more
      detailed output is desired

    OUTPUT:

    A dictionary ``P`` whose keys are residue field vectors and whose values are lists of all exponent vectors
    which correspond to the given residue field vector.

    .. NOTE::

        - For example, if ``rfv_dictionary[e0] = r0``, then ``P[r0]`` is a list which contains ``e0``.
        - During construction, some residue field vectors can be eliminated as coming from
          solutions to the `S`-unit equation. Such vectors are dropped from the keys of the dictionary ``P``.

    EXAMPLES:

    In this example, we use a truncated list generated when solving the `S`-unit equation in the case that `K` is defined by the
    polynomial `x^2+x+1` and `S` consists of the primes above 3::

        sage: from sage.rings.number_field.S_unit_solver import construct_rfv_to_ev
        sage: rfv_dict = {(1, 3): [3, 2], (3, 0): [6, 6], (5, 4): [3, 6], (2, 1): [4, 6],
        ....:             (4, 0): [4, 2], (1, 2): [5, 6]}
        sage: construct_rfv_to_ev(rfv_dict,7,2,False)
        {(3, 2): [(1, 3)], (4, 2): [(4, 0)], (4, 6): [(2, 1)], (5, 6): [(1, 2)]}
    """

    # The keys in P are just the possible first entries of a residue field vector.
    # The values (all empty lists now) will be added in the next step.

    P = {(v,) : [] for v in range(2, q)}

    # Step 1. Populate the empty lists in P[(v,)].
    # Loop through the keys in rfv_dictionary. For each, look at the output rf_vector.
    # Find the key in P which matches the first entry of the rf_vector.
    # Dump the **rest** of the rf_vector into a pair [exp_vec, rf_vec[1:]],
    # and append this pair into the dictionary P at the key (rf_vec[0], ).

    # Now, P[(v,)] = [ [a_0, e_0], [a_1, e_1], ...]
    #
    # The relationship between v, a_i, and e_i is as follows:
    #
    # a_i is an exponent vector, whose associated residue field vector is the
    # concatenation of v with e_i.

    for exponent_vector in rfv_dictionary:
        residue_field_vector = rfv_dictionary[exponent_vector]
        rf_vector_start = (residue_field_vector[0], )
        rf_vector_end = residue_field_vector[1:]
        P[rf_vector_start].append([exponent_vector, rf_vector_end])

    if verbose:
        print("Populated P. Currently it has ", len(P), "keys.")

    # Step 2: We build a new dictionary, P_new, from P.
    #
    # This is a step that will be repeated, once for each of the d primes over q.
    #
    # P is a dictionary whose keys are tuples of length m, representing the beginning of known residue field vectors.
    #
    # For any such beginning `s`,
    #
    # P[s] = [ [a_0, e_0], [a_1, e_1], ...]
    #
    # where for any exponent vector a_i, the associated residue field vector is the concatenation s + e_i.
    #
    # The dictionary P_new is constructed from the dictionary P. The new keys will be tuples of length m + 1.
    #
    # During the construction, we look for impossible entries for S-unit solutions, and drop them from the dictionary as needed.

    for j in range(d-1):
        if verbose:
            print("Constructing ", j, " th place of the residue field vectors, out of ", d-1, " total.")
        P_new = {}
        garbage = {}

        # we loop over each key of P.
        for rf_vector_start in P:

            # each key of P provides q-2 possible keys for P_new, which we introduce and assign an empty list.
            for w in range(2, q):
                new_rf_vector_start = tuple(list(rf_vector_start) + [w])
                P_new[new_rf_vector_start] = []

            # we populate P_new[ new_rf_vector_start ] using P[rf_vector_start]
            for exponent_vector, rf_vector_end in P[rf_vector_start]:
                new_rf_vector_end = rf_vector_end[1:]
                w = rf_vector_end[0]
                new_rf_vector_start = tuple(list(rf_vector_start) + [w])
                P_new[new_rf_vector_start].append([exponent_vector, new_rf_vector_end])

        if verbose:
            print("P_new is populated with ", len(P_new), " keys.")

        # we now loop over the keys of P_new, looking for incompatible entries.

        for rf_vector_start in P_new:
            # the final entry of rf_vector_start or rf_vector_complement_start must be < (q+3)/2.
            # No loss to insist that it is rf_vector_start.
            if rf_vector_start[-1] < (q+3)/2:
                # we find the complement to rf_vector_start:
                rf_vector_complement_start = tuple([ q+1-j for j in rf_vector_start])
                if P_new[ rf_vector_start ] == [] or P_new[rf_vector_complement_start] == []:
                    # these can't be solutions. Mark them for deletion.
                    garbage[rf_vector_start] = True
                    garbage[rf_vector_complement_start] = True

        # garbage removal
        for rf_vector_start in garbage:
            P_new.pop(rf_vector_start, 0)

        if verbose:
            print("After removing incompatible entries, P_new is down to ", len(P_new), " keys.")

        # Time to move on to the next dictionary.
        P = P_new.copy()

    # Now, we just clean up P.
    for residue_field_vector in P:
        # at this instant, P[ residue_field_vector ] is a list of pairs: [ [a0,e0], ... ]
        # We only care about the exponent vectors a0,...
        P[residue_field_vector] = [a[0] for a in P[residue_field_vector]]

    if verbose:
        print("Returning dictionary P with ", len(P), " keys.")
    return P.copy()


def construct_comp_exp_vec(rfv_to_ev_dict, q):
    r"""
    Construct a dictionary associating complement vectors to residue field vectors.

    INPUT:

    - ``rfv_to_ev_dict`` -- dictionary whose keys are residue field vectors and
      whose values are lists of exponent vectors with the associated residue
      field vector
    - ``q`` -- the characteristic of the residue field

    OUTPUT: a dictionary whose typical key is an exponent vector ``a``, and
    whose associated value is a list of complementary exponent vectors to ``a``

    EXAMPLES:

    In this example, we use the list generated when solving the `S`-unit equation in the case that `K` is defined by the
    polynomial `x^2+x+1` and `S` consists of the primes above 3

    ::

        sage: from sage.rings.number_field.S_unit_solver import construct_comp_exp_vec
        sage: rfv_to_ev_dict = {(6, 6): [(3, 0)], (5, 6): [(1, 2)], (5, 4): [(5, 3)],
        ....:                   (6, 2): [(5, 5)], (2, 5): [(0, 1)], (5, 5): [(3, 4)],
        ....:                   (4, 4): [(0, 2)], (6, 3): [(1, 4)], (3, 6): [(5, 4)],
        ....:                   (2, 2): [(0, 4)], (3, 5): [(1, 0)], (6, 4): [(1, 1)],
        ....:                   (3, 2): [(1, 3)], (2, 6): [(4, 5)], (4, 5): [(4, 3)],
        ....:                   (2, 3): [(2, 3)], (4, 2): [(4, 0)], (6, 5): [(5, 2)],
        ....:                   (3, 3): [(3, 2)], (5, 3): [(5, 0)], (4, 6): [(2, 1)],
        ....:                   (3, 4): [(3, 5)], (4, 3): [(0, 5)], (5, 2): [(3, 1)],
        ....:                   (2, 4): [(2, 0)]}
        sage: construct_comp_exp_vec(rfv_to_ev_dict, 7)
        {(0, 1): [(1, 4)],
         (0, 2): [(0, 2)],
         (0, 4): [(3, 0)],
         (0, 5): [(4, 3)],
         (1, 0): [(5, 0)],
         (1, 1): [(2, 0)],
         (1, 2): [(1, 3)],
         (1, 3): [(1, 2)],
         (1, 4): [(0, 1)],
         (2, 0): [(1, 1)],
         (2, 1): [(4, 0)],
         (2, 3): [(5, 2)],
         (3, 0): [(0, 4)],
         (3, 1): [(5, 4)],
         (3, 2): [(3, 4)],
         (3, 4): [(3, 2)],
         (3, 5): [(5, 3)],
         (4, 0): [(2, 1)],
         (4, 3): [(0, 5)],
         (4, 5): [(5, 5)],
         (5, 0): [(1, 0)],
         (5, 2): [(2, 3)],
         (5, 3): [(3, 5)],
         (5, 4): [(3, 1)],
         (5, 5): [(4, 5)]}
    """

    comp_exp_vec_dict = {}
    for residue_field_vector in rfv_to_ev_dict:
        rf_vector_complement = tuple([q + 1 - j for j in residue_field_vector])
        exponent_vector_list = rfv_to_ev_dict[ residue_field_vector ][:]
        exponent_vector_complement_list = rfv_to_ev_dict[rf_vector_complement][:]
        for exponent_vector in exponent_vector_list:
            comp_exp_vec_dict[exponent_vector] = exponent_vector_complement_list
    return comp_exp_vec_dict


def drop_vector(ev, p, q, complement_ev_dict):
    r"""
    Determine if the exponent vector, ``ev``, may be removed from the complement dictionary during construction.
    This will occur if ``ev`` is not compatible with an exponent vector mod `q-1`.

    INPUT:

    - ``ev`` -- an exponent vector modulo `p-1`
    - ``p`` -- the prime such that ``ev`` is an exponent vector modulo `p-1`
    - ``q`` -- a prime, distinct from `p`, that is a key in the ``complement_ev_dict``
    - ``complement_ev_dict`` -- dictionary of dictionaries, whose keys are primes
      ``complement_ev_dict[q]`` is a dictionary whose keys are exponent vectors modulo `q-1`
      and whose values are lists of complementary exponent vectors modulo `q-1`

    OUTPUT: ``True`` if ``ev`` may be dropped from the complement exponent
    vector dictionary, and ``False`` if not

    .. NOTE::

        - If ``ev`` is not compatible with any of the vectors modulo `q-1`, then it can no longer correspond to a solution
          of the `S`-unit equation. It returns ``True`` to indicate that it should be removed.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import drop_vector
        sage: drop_vector((1, 2, 5), 7, 11, {11: {(1, 1, 3): [(1, 1, 3), (2, 3, 4)]}})
        True

    ::

        sage: P = {3: {(1, 0, 0): [(1, 0, 0), (0, 1, 0)],
        ....:          (0, 1, 0): [(1, 0, 0), (0, 1, 0)]},
        ....:      7: {(0, 3, 4): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        ....:          (1, 2, 4): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        ....:          (0, 1, 2): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        ....:          (0, 5, 4): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        ....:          (1, 4, 2): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        ....:          (1, 0, 4): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        ....:          (0, 3, 2): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        ....:          (1, 0, 0): [(0, 5, 4), (0, 3, 2), (0, 1, 0)],
        ....:          (1, 2, 0): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        ....:          (0, 1, 0): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        ....:          (0, 5, 0): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        ....:          (1, 2, 2): [(0, 5, 4), (0, 3, 2), (0, 1, 0)],
        ....:          (1, 4, 0): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        ....:          (1, 0, 2): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        ....:          (1, 4, 4): [(0, 5, 4), (0, 3, 2), (0, 1, 0)]}}
        sage: drop_vector((0, 1, 0), 3, 7, P)
        False
    """
    # returns True if it is OK to drop exp_vec given the current comp_exp_vec dictionary associated to some q.
    # returns False otherwise
    # loop over the possible compatible vectors in the other modulus
    g = gcd(p-1, q-1)
    for compatible_exp_vec in compatible_vectors(ev, p-1, q-1, g):
        # do they appear in the other dictionary?
        if compatible_exp_vec in complement_ev_dict[q]:
            # OK, but the complements need to be compatible, too!
            ev_complement_list = complement_ev_dict[p][ev]
            for ev_comp in ev_complement_list:
                for compatible_cv in compatible_vectors(ev_comp, p-1, q-1, g):
                    if compatible_cv in complement_ev_dict[q][compatible_exp_vec]:
                        return False
    return True


def construct_complement_dictionaries(split_primes_list, SUK, verbose=False):
    r"""
    Construct the complement exponent vector dictionaries.

    INPUT:

    - ``split_primes_list`` -- list of rational primes which split completely
      in the number field `K`
    - ``SUK`` -- the `S`-unit group for a number field `K`
    - ``verbose`` -- boolean to provide additional feedback (default: ``False``)

    OUTPUT:

    A dictionary of dictionaries. The keys coincide with the primes in
    ``split_primes_list``. For each `q`, ``comp_exp_vec[q]`` is a dictionary
    whose keys are exponent vectors modulo `q-1`, and whose values are lists of
    exponent vectors modulo `q-1`.

    If `w` is an exponent vector in ``comp_exp_vec[q][v]``, then the residue field vectors modulo `q` for
    `v` and `w` sum to ``[1,1,...,1]``

    .. NOTE::

        - The data of ``comp_exp_vec`` will later be lifted to `\ZZ` to look for true `S`-Unit equation solutions.
        - During construction, the various dictionaries are compared to each other several times to
          eliminate as many mod `q` solutions as possible.
        - The authors acknowledge a helpful discussion with Norman Danner which helped formulate this code.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import construct_complement_dictionaries
        sage: x = polygen(ZZ, 'x')
        sage: f = x^2 + 5
        sage: H = 10
        sage: K.<xi> = NumberField(f)
        sage: SUK = K.S_unit_group(S=K.primes_above(H))
        sage: split_primes_list = [3, 7]
        sage: actual = construct_complement_dictionaries(split_primes_list, SUK)
        sage: expected = {3: {(0, 1, 0): [(1, 0, 0), (0, 1, 0)],
        ....:                 (1, 0, 0): [(1, 0, 0), (0, 1, 0)]},
        ....:             7: {(0, 1, 0): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        ....:                 (0, 1, 2): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        ....:                 (0, 3, 2): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        ....:                 (0, 3, 4): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        ....:                 (0, 5, 0): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        ....:                 (0, 5, 4): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        ....:                 (1, 0, 0): [(0, 5, 4), (0, 3, 2), (0, 1, 0)],
        ....:                 (1, 0, 2): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        ....:                 (1, 0, 4): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        ....:                 (1, 2, 0): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        ....:                 (1, 2, 2): [(0, 5, 4), (0, 3, 2), (0, 1, 0)],
        ....:                 (1, 2, 4): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        ....:                 (1, 4, 0): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        ....:                 (1, 4, 2): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        ....:                 (1, 4, 4): [(0, 5, 4), (0, 3, 2), (0, 1, 0)]}}
        sage: all(set(actual[p][vec]) == set(expected[p][vec])
        ....:     for p in [3, 7] for vec in expected[p])
        True
    """
    # We initialize some dictionaries.

    rho = SUK.gens_values()
    rho_length = len(rho)
    rho_images_dict = {}
    rho_orders_dict = {}

    K = SUK.number_field()
    for q in split_primes_list:
        ideals_over_q, residue_fields, rho_images, product_rho_orders = sieve_ordering(SUK, q)
        rho_images_dict[q] = rho_images
        rho_orders_dict[q] = product_rho_orders

    nK = K.absolute_degree()
    w0 = rho[0].multiplicative_order()

    # We build a dictionary of dictionaries.
    # rfv_to_ev[q] is the 'mod q' residue field vector to exponent vector dictionary.

    rfv_to_ev = {}

    # We build a second dictionary of dictionaries.
    # comp_exp_vec[q] is the dictionary mod q which assigns to each exponent vector
    # a list of 'complementary' exponent vectors.

    comp_exp_vec = {}

    q0 = split_primes_list[0]

    if verbose:
        print("Using the following primes: ", split_primes_list)
    for q in split_primes_list:
        rho_images = rho_images_dict[q]
        if verbose:
            print("q = ", q)

        def epsilon_q(a, i):
            # a is an exponent vector
            # i is an index for one of the primes over q
            # returns the value of rho_j^a_j inside the
            # residue field of Qi. (Necessarily isomorphic to F_q.)
            # rho_images[i][j] == rho[j] modulo Q[i]
            eps_value = rho_images[i][0]**a[0]
            for j in range(1, rho_length):
                eps_value *= rho_images[i][j]**a[j]
            return eps_value

        if verbose:
            print("The evaluation function epsilon has been defined using rho_images = ", rho_images)
        # Now, we run through the vectors in the iterator, but only keep the ones
        # which are compatible with the previously constructed dictionaries. That is,
        # in order to keep an exp_vec mod q, there must exist a compatible exp_vec mod p
        # in the keys of the rfv_to_ev[p] dictionary for each completely split prime
        # p appearing prior to q in split_primes_list.

        if q == q0:
            # for the first prime, there is no filtering possible, and we just build the exponent vector
            # iterator.

            # This should consist of all vectors (a0,...,a_{t-1}), where
            # a0 is in the range 0 .. w_0 - 1 and
            # aj is in the range 0 .. q - 2   (for j > 0)
            ranges = [range(w0)] + [range(q-1) for _ in range(rho_length-1)]
            ev_iterator = itertools.product(*ranges)

            # With the iterator built, we construct the exponent vector to residue field dictionary.

            ev_to_rfv_dict = {ev : [epsilon_q(ev, i) for i in range(nK)] for ev in ev_iterator}

            if verbose:
                print("The residue field dictionary currently has ", len(ev_to_rfv_dict), " exponent vector keys.")
        else:
            ev_to_rfv_dict = {}
            # We use compatibility requirements to keep the size of the dictionary down.
            # Later on, we'll compare all dictionaries pairwise. But for now, we just
            # check against the first.

            # That is, rather than loop over every possible exponent vector mod q-1,
            # we only consider those evs which are compatible with the mod q0 - 1 vectors.

            # Loop over exponent vectors modulo q0 - 1
            g = gcd(q0-1, q-1)
            for exp_vec_mod_q0 in comp_exp_vec[q0]:
                # Loop only over exponent vectors modulo q-1 which are compatible with exp_vec_mod_q0
                for exp_vec in compatible_vectors(exp_vec_mod_q0, q0-1, q-1, g):
                    # fill the dictionary with the residue field vectors using the evaluation function.
                    ev_to_rfv_dict[exp_vec] = [epsilon_q(exp_vec, i) for i in range(nK)]

        if verbose:
            print("The residue field dictionary currently has ", len(ev_to_rfv_dict), " exponent vector keys.")
        # At this point, we now have a dictionary ev_to_rfv_dict, which attaches
        # to each exponent vector a 'residue field vector,' which is a tuple of the
        # nK values epsilon_q(a,0),...,epsilon_q(a,nK-1).

        clean_rfv_dict( ev_to_rfv_dict )

        if verbose:
            print("clean_rfv_dict executed.")
            print("The residue field dictionary currently has ", len(ev_to_rfv_dict), " exponent vector keys.")
        # We essentially construct an inverse dictionary: one whose keys are residue field vectors,
        # and whose values are the exponent vectors that yield each key

        rfv_to_ev[q] = construct_rfv_to_ev(ev_to_rfv_dict, q, nK, verbose=verbose)

        if verbose:
            print("construct_rfv_to_ev executed.")
            print("The rfv_to_ev dictionary currently has ", len(rfv_to_ev[q]), "rfv keys.")

        comp_exp_vec[q] = construct_comp_exp_vec(rfv_to_ev[q], q)

        if verbose:
            print("construct_comp_exp_vec executed.")
            print("Size of comp_exp_vec[q]: ", len(comp_exp_vec[q]))

        # Now that we have a new dictionary, we compare all the dictionaries pairwise,
        # looking for opportunities to remove 'impossible' solutions.

        for p in comp_exp_vec.keys():
            if p == q:
                continue
            if verbose:
                print("Comparing dictionaries for p = ", p, "and q = ", q, ".")

            old_size_p = len(comp_exp_vec[p])

            if verbose:
                print("Size of comp_exp_vec[p] is: ", old_size_p, ".")
                cv_size = ((q-1)/gcd(p-1, q-1)) ** (rho_length - 1)
                print("Length of compatible_vectors: ", cv_size, ".")
                print("Product: ", old_size_p*cv_size)

            for exp_vec in list(comp_exp_vec[p]):
                if drop_vector(exp_vec, p, q, comp_exp_vec):
                    comp_exp_vec[p].pop(exp_vec)

            if verbose:
                print("Shrunk dictionary p from ", old_size_p, " to ", len(comp_exp_vec[p]))

            # Now, repeat, but swap p and q.

            old_size_q = len(comp_exp_vec[q])

            if verbose:
                print("Size of comp_exp_vec[q] is: ", old_size_q, ".")
                cv_size = ((p-1)/gcd(p-1, q-1)) ** (rho_length - 1)
                print("Length of compatible_vectors: ", cv_size, ".")
                print("Product: ", old_size_q * cv_size)

            for exp_vec in list(comp_exp_vec[q]):
                if drop_vector(exp_vec, q, p, comp_exp_vec):
                    comp_exp_vec[q].pop(exp_vec)

            if verbose:
                print("Shrunk dictionary q from ", old_size_q, " to ", len(comp_exp_vec[q]))

    return comp_exp_vec


def compatible_vectors_check(a0, a1, g, l):
    r"""
    Given exponent vectors with respect to two moduli, determine if they are compatible.

    INPUT:

    - ``a0`` -- an exponent vector modulo ``m0``
    - ``a1`` -- an exponent vector modulo ``m1`` (must have the same length as ``a0``)
    - ``g`` -- the gcd of ``m0`` and ``m1``
    - ``l`` -- the length of ``a0`` and of ``a1``

    OUTPUT: ``True`` if there is an integer exponent vector a satisfying

    .. MATH::

        \begin{aligned}
        a[0] &== a0[0] == a1[0]\\
        a[1:] &== a0[1:] \mod m_0\\
        a[1:] &== a1[1:] \mod m_1
        \end{aligned}

    and ``False`` otherwise.

    .. NOTE::

        - Exponent vectors must agree exactly in the first coordinate.
        - If exponent vectors are different lengths, an error is raised.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import compatible_vectors_check
        sage: a0 = (3, 1, 8, 11)
        sage: a1 = (3, 5, 6, 13)
        sage: a2 = (5, 5, 6, 13)
        sage: compatible_vectors_check(a0, a1, gcd(12, 22), 4r)
        True
        sage: compatible_vectors_check(a0, a2, gcd(12, 22), 4r)
        False
    """
    # exponent vectors must agree exactly in the 0th coordinate.
    return a0[0] == a1[0] and all((x0 - x1) % g == 0 for x0,x1 in zip(itertools.islice(a0, 1, l), itertools.islice(a1, 1, l)))


def compatible_vectors(a, m0, m1, g):
    r"""
    Given an exponent vector ``a`` modulo ``m0``, return an iterator over the exponent vectors for the modulus ``m1``, such that a lift to the lcm modulus exists.

    INPUT:

    - ``a`` -- an exponent vector for the modulus ``m0``
    - ``m0`` -- positive integer (specifying the modulus for ``a``)
    - ``m1`` -- positive integer (specifying the alternate modulus)
    - ``g`` -- the gcd of ``m0`` and ``m1``

    OUTPUT: list of exponent vectors modulo ``m1`` which are compatible with ``a``

    .. NOTE::

        Exponent vectors must agree exactly in the 0th position in order to be
        compatible.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import compatible_vectors
        sage: a = (3, 1, 8, 1)
        sage: list(compatible_vectors(a, 18, 12, gcd(18,12)))
        [(3, 1, 2, 1),
         (3, 1, 2, 7),
         (3, 1, 8, 1),
         (3, 1, 8, 7),
         (3, 7, 2, 1),
         (3, 7, 2, 7),
         (3, 7, 8, 1),
         (3, 7, 8, 7)]

    The order of the moduli matters. ::

        sage: len(list(compatible_vectors(a, 18, 12, gcd(18,12))))
        8
        sage: len(list(compatible_vectors(a, 12, 18, gcd(18,12))))
        27
    """
    # recall that the 0th entry must be an exact match.
    ranges = [[a[0]]] + [range(a[i] % g, (a[i] % g) + m1, g) for i in range(1, len(a))]
    return itertools.product(*ranges)


def compatible_systems(split_prime_list, complement_exp_vec_dict):
    r"""
    Given dictionaries of complement exponent vectors for various primes that
    split in `K`, compute all possible compatible systems.

    INPUT:

    - ``split_prime_list`` -- list of rational primes that split completely in `K`
    - ``complement_exp_vec_dict`` -- dictionary of dictionaries; the keys are
      primes from ``split_prime_list``

    OUTPUT: list of compatible systems of exponent vectors

    .. NOTE::

        - For any `q` in ``split_prime_list``, ``complement_exp_vec_dict[q]`` is a dictionary whose keys are exponent vectors modulo `q-1`
          and whose values are lists of exponent vectors modulo `q-1` which are complementary to the key.

        - An item in ``system_list`` has the form ``[ [v0, w0], [v1, w1], ..., [vk, wk] ]``, where::

          - ``qj = split_prime_list[j]``
          - ``vj`` and ``wj`` are complementary exponent vectors modulo ``qj - 1``
          - the pairs are all simultaneously compatible.

        - Let ``H = lcm( qj - 1 : qj in split_primes_list )``. Then for any compatible system, there is at most one pair of integer
          exponent vectors ``[v, w]`` such that::

          - every entry of ``v`` and ``w`` is bounded in absolute value by ``H``
          - for any ``qj``, ``v`` and ``vj`` agree modulo ``(qj - 1)``
          - for any ``qj``, ``w`` and ``wj`` agree modulo ``(qj - 1)``

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import compatible_systems
        sage: split_primes_list = [3, 7]
        sage: checking_dict = {3: {(0, 1, 0): [(1, 0, 0)]}, 7: {(0, 1, 0): [(1, 0, 0)]}}
        sage: compatible_systems(split_primes_list, checking_dict)
        [[[(0, 1, 0), (1, 0, 0)], [(0, 1, 0), (1, 0, 0)]]]
    """
    S0 = split_prime_list

    system_list = []
    if len(S0) == 1:
        q = S0[0]
        for exponent_vector in complement_exp_vec_dict[q]:
            for complementary_vector in complement_exp_vec_dict[q][exponent_vector]:
                pair = [[exponent_vector, complementary_vector]]
                system_list.append(pair)
    elif len(S0) > 1:
        S1 = S0[:-1]
        old_systems = compatible_systems(S1, complement_exp_vec_dict)
        q = S0[-1]
        gcds = [gcd(q-1, qj-1) for qj in S1]
        for exp_vec in complement_exp_vec_dict[q]:
            l = len(exp_vec)
            for comp_vec in complement_exp_vec_dict[q][exp_vec]:
                for old_system in old_systems:
                    if all((compatible_vectors_check(exp_vec, exp_vec_qj, g, l) and
                            compatible_vectors_check(comp_vec, comp_vec_qj, g, l))
                           for g, (exp_vec_qj, comp_vec_qj) in zip(gcds, old_system)):
                        # build the new system and append it to the list.
                        new_system = old_system + [[exp_vec, comp_vec]]
                        system_list.append(new_system)
    return system_list


def compatible_system_lift(compatible_system, split_primes_list):
    r"""
    Given a compatible system of exponent vectors and complementary exponent
    vectors, return a lift to the integers.

    INPUT:

    - ``compatible_system`` -- list of pairs ``[ [v0, w0], [v1, w1], .., [vk, wk] ]``
      where [vi, wi] is a pair of complementary exponent vectors modulo ``qi - 1``,
      and all pairs are compatible
    - ``split_primes_list`` -- list of primes ``[ q0, q1, .., qk ]``

    OUTPUT: a pair of vectors ``[v, w]`` satisfying:

    1. ``v[0] == vi[0]`` for all ``i``
    2. ``w[0] == wi[0]`` for all ``i``
    3. ``v[j] == vi[j]`` modulo ``qi - 1`` for all ``i`` and all ``j > 0``
    4. ``w[j] == wi[j]`` modulo ``qi - 1`` for all ``i`` and all `j > 0``
    5. every entry of ``v`` and ``w`` is bounded by ``L/2`` in absolute value, where ``L`` is the least common multiple of ``{qi - 1 : qi in split_primes_list }``

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import compatible_system_lift
        sage: split_primes_list = [3, 7]
        sage: comp_sys = [[(0, 1, 0), (0, 1, 0)], [(0, 3, 4), (0, 1, 2)]]
        sage: compatible_system_lift(comp_sys, split_primes_list)
        [(0, 3, -2), (0, 1, 2)]
    """

    if len(split_primes_list) != len(compatible_system):
        raise ValueError("The number of primes does not match the length of the given exponent vectors.")

    # the first entries are already determined.
    exponent_vector_lift = [ZZ(compatible_system[0][0][0])]
    complement_vector_lift = [ZZ(compatible_system[0][1][0])]

    # fill in exponent_vector_lift
    moduli_list = [q-1 for q in split_primes_list]
    L = lcm(moduli_list)

    t = len(compatible_system[0][0])
    for i in range(1,t):
        exp_coord_residues = [pair[0][i] for pair in compatible_system]
        comp_coord_residues = [pair[1][i] for pair in compatible_system]

        ev_lift_coordinate = CRT(exp_coord_residues, moduli_list)
        cv_lift_coordinate = CRT(comp_coord_residues, moduli_list)

        # these values lie in the range [0, L-1], so we must shift them if they are bigger than L/2.

        if ev_lift_coordinate > L/2:
            ev_lift_coordinate -= L
        if cv_lift_coordinate > L/2:
            cv_lift_coordinate -= L

        exponent_vector_lift.append(ev_lift_coordinate)
        complement_vector_lift.append(cv_lift_coordinate)

    return [tuple(exponent_vector_lift), tuple(complement_vector_lift)]


def solutions_from_systems(SUK, bound, cs_list, split_primes_list):
    r"""
    Lift compatible systems to the integers and return the `S`-unit equation
    solutions that the lifts yield.

    INPUT:

    - ``SUK`` -- the group of `S`-units where we search for solutions
    - ``bound`` -- a bound for the entries of all entries of all lifts
    - ``cs_list`` -- list of compatible systems of exponent vectors modulo
      `q-1` for various primes `q`
    - ``split_primes_list`` -- list of primes giving the moduli of the exponent
      vectors in ``cs_list``

    OUTPUT:

    A list of solutions to the S-unit equation. Each solution is a list:

    1. an exponent vector over the integers, ``ev``
    2. an exponent vector over the integers, ``cv``
    3. the S-unit corresponding to ``ev``, ``iota_exp``
    4. the S-unit corresponding to ``cv``, ``iota_comp``

    .. NOTE::

        - Every entry of ``ev`` is less than or equal to bound in absolute value
        - every entry of ``cv`` is less than or equal to bound in absolute value
        - ``iota_exp + iota_comp == 1``

    EXAMPLES:

    Given a single compatible system, a solution can be found. ::

        sage: from sage.rings.number_field.S_unit_solver import solutions_from_systems
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^2 - 15)
        sage: SUK = K.S_unit_group(S=K.primes_above(2))
        sage: split_primes_list = [7, 17]
        sage: a_compatible_system = [[[(0, 0, 5), (0, 0, 5)], [(0, 0, 15), (0, 0, 15)]]]
        sage: solutions_from_systems(SUK, 20, a_compatible_system, split_primes_list)
        [((0, 0, -1), (0, 0, -1), 1/2, 1/2)]
    """
    solutions = []

    for system in cs_list:
        ev, cv = compatible_system_lift(system, split_primes_list)
        if all(abs(x) <= bound for x in ev[1:] + cv[1:]):
            # the entries are all below the bound, so there is nothing left to do
            # except construct the elements and see if they are solutions to
            # the S-unit equation
            iota_exp = SUK.exp( ev )
            iota_comp = SUK.exp( cv )
            if iota_exp + iota_comp == 1:
                sol = ( ev, cv, iota_exp, iota_comp )
                solutions.append( sol )

    return solutions


def clean_sfs(sfs_list):
    r"""
    Given a list of `S`-unit equation solutions, remove trivial redundancies.

    INPUT:

    - ``sfs_list`` -- list of solutions to the `S`-unit equation

    OUTPUT: list of solutions to the `S`-unit equation

    .. NOTE::

        The function looks for cases where `x + y = 1` and `y + x = 1` appear
        as separate solutions, and removes one.

    EXAMPLES:

    The function is not dependent on the number field and removes redundancies in any list. ::

        sage: from sage.rings.number_field.S_unit_solver import clean_sfs
        sage: sols = [((1, 0, 0), (0, 0, 1), -1, 2), ((0, 0, 1), (1, 0, 0), 2, -1)]
        sage: clean_sfs( sols )
        [((1, 0, 0), (0, 0, 1), -1, 2)]
    """
    # given the output from solutions_from_systems,
    # look for trivial redundancies: swapping exp_vec, comp_vec, particularly.
    new_sfs = []
    for entry in sfs_list:
        swapped_entry = (entry[1], entry[0], entry[3], entry[2])
        if entry not in new_sfs and swapped_entry not in new_sfs:
            new_sfs.append(entry)
    return new_sfs


def sieve_below_bound(K, S, bound=10, bump=10, split_primes_list=[], verbose=False):
    r"""
    Return all solutions to the `S`-unit equation `x + y = 1` over `K` with
    exponents below the given bound.

    INPUT:

    - ``K`` -- a number field (an absolute extension of the rationals)
    - ``S`` -- list of finite primes of `K`
    - ``bound`` -- positive integer upper bound for exponents, solutions with
      exponents having absolute value below this bound will be found (default: 10)
    - ``bump`` -- positive integer by which the minimum LCM will be increased
      if not enough split primes are found in sieving step (default: 10)
    - ``split_primes_list`` -- list of rational primes that split completely in
      the extension `K/\QQ`, used for sieving. For complete list of solutions
      should have lcm of `\{(p_i-1)\} for primes `p_i` greater than bound
      (default: ``[]``).
    - ``verbose`` -- an optional parameter allowing the user to print
      information during the sieving process (default: ``False``)

    OUTPUT:

    A list of tuples `[(A_1, B_1, x_1, y_1), (A_2, B_2, x_2, y_2), \dots (A_n, B_n, x_n, y_n)]` such that:

    1. The first two entries are tuples `A_i = (a_0, a_1, \dots, a_t)` and `B_i = (b_0, b_1, \dots, b_t)` of exponents.
    2. The last two entries are `S`-units `x_i` and `y_i` in `K` with `x_i + y_i = 1`.
    3. If the default generators for the `S`-units of `K` are `(\rho_0, \rho_1, \dots, \rho_t)`,
       then these satisfy `x_i = \prod(\rho_i)^{(a_i)}` and `y_i = \prod(\rho_i)^{(b_i)}`.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import sieve_below_bound, eq_up_to_order
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^2 + x + 1)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: S = SUK.primes()
        sage: sols = sieve_below_bound(K, S, 10)
        sage: expected = [((1, -1), (0, -1), 1/3*xi + 2/3, -1/3*xi + 1/3),
        ....:             ((0, 1), (4, 0), xi + 2, -xi - 1),
        ....:             ((2, 0), (5, 1), xi, -xi + 1),
        ....:             ((1, 0), (5, 0), xi + 1, -xi)]
        sage: eq_up_to_order(sols, expected)
        True
    """
    SUK = UnitGroup(K, S=tuple(S))
    initial_bound = bound

    while not split_primes_list:
        try:
            split_primes_list = split_primes_large_lcm(SUK, initial_bound)
        except ValueError:
            initial_bound += bump
            print("Couldn't find enough split primes. Bumping to ", initial_bound)

    if not K.is_absolute():
        raise ValueError("K must be an absolute extension.")

    complement_exp_vec_dict = construct_complement_dictionaries(split_primes_list, SUK, verbose=verbose)

    cs_list = compatible_systems(split_primes_list, complement_exp_vec_dict)

    sfs_list = solutions_from_systems(SUK, bound, cs_list, split_primes_list)

    S_unit_solutions = clean_sfs(sfs_list)

    return S_unit_solutions


def solve_S_unit_equation(K, S, prec=106, include_exponents=True, include_bound=False, proof=None, verbose=False):
    r"""
    Return all solutions to the `S`-unit equation `x + y = 1` over `K`.

    INPUT:

    - ``K`` -- a number field (an absolute extension of the rationals)
    - ``S`` -- list of finite primes of `K`
    - ``prec`` -- precision used for computations in real, complex, and `p`-adic
      fields (default: 106)
    - ``include_exponents`` -- whether to include the exponent vectors in the
      returned value (default: ``True``)
    - ``include_bound`` -- whether to return the final computed bound
      (default: ``False``)
    - ``verbose`` -- whether to print information during the sieving step
      (default: ``False``)

    OUTPUT:

    A list of tuples `[(A_1, B_1, x_1, y_1), (A_2, B_2, x_2, y_2), \dots (A_n, B_n, x_n, y_n)]` such that:

    1. The first two entries are tuples `A_i = (a_0, a_1, \dots, a_t)` and `B_i = (b_0, b_1, \dots, b_t)` of exponents.  These will be omitted if ``include_exponents`` is ``False``.
    2. The last two entries are `S`-units `x_i` and `y_i` in `K` with `x_i + y_i = 1`.
    3. If the default generators for the `S`-units of `K` are `(\rho_0, \rho_1, \dots, \rho_t)``, then these satisfy `x_i = \prod(\rho_i)^{(a_i)}` and `y_i = \prod(\rho_i)^{(b_i)}`.

    If ``include_bound``, will return a pair ``(sols, bound)`` where ``sols`` is as above and ``bound`` is the bound used for the entries in the exponent vectors.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import solve_S_unit_equation, eq_up_to_order
        sage: x = polygen(ZZ, 'x')
        sage: K.<xi> = NumberField(x^2 + x + 1)
        sage: S = K.primes_above(3)
        sage: sols = solve_S_unit_equation(K, S, 200)
        sage: expected = [((0, 1), (4, 0), xi + 2, -xi - 1),
        ....:             ((1, -1), (0, -1), 1/3*xi + 2/3, -1/3*xi + 1/3),
        ....:             ((1, 0), (5, 0), xi + 1, -xi),
        ....:             ((2, 0), (5, 1), xi, -xi + 1)]
        sage: eq_up_to_order(sols, expected)
        True

    In order to see the bound as well, use the optional parameter ``include_bound``::

        sage: solutions, bound = solve_S_unit_equation(K, S, 100, include_bound=True)
        sage: bound
        7

    You can omit the exponent vectors::

        sage: sols = solve_S_unit_equation(K, S, 200, include_exponents=False)
        sage: expected = [(xi + 2, -xi - 1), (1/3*xi + 2/3, -1/3*xi + 1/3),
        ....:             (-xi, xi + 1), (-xi + 1, xi)]
        sage: set(frozenset(a) for a in sols) == set(frozenset(b) for b in expected)
        True

    It is an error to use values in S that are not primes in K::

        sage: solve_S_unit_equation(K, [3], 200)
        Traceback (most recent call last):
        ...
        ValueError: S must consist only of prime ideals,
        or a single element from which a prime ideal can be constructed.

    We check the case that the rank is 0::

        sage: K.<xi> = NumberField(x^2 + x + 1)
        sage: solve_S_unit_equation(K, [])
        [((1,), (5,), xi + 1, -xi)]
    """

    # Checks to make sure inputs are legal
    # K must be an absolute extension:
    if not K.is_absolute():
        raise ValueError("K must be an absolute extension.")
    # S must be a finite set of primes
    try:
        SUK = UnitGroup(K, proof=proof, S=tuple(S))
    except Exception:
        raise ValueError("S must consist only of prime ideals, or a single element from which a prime ideal can be constructed.")

    # Gather the roots of unity of the number field
    A = K.roots_of_unity()
    w = K.number_of_roots_of_unity()
    if SUK.rank() == 0:
        # Since the rank is 0, K is imaginary quadratic and S is empty
        # Only possibilities are combinations of roots of unity
        # and this can only occur when there are 6 roots of unity, when
        # (1+sqrt(-3))/2 + (1-sqrt(-3))/2 = 1 is the unique solution.
        if len(A) == 6:
            S_unit_solutions = [((ZZ(1),), (ZZ(5),), A[0], A[-2])]
        else:
            S_unit_solutions = []
    else:
        # First find a bound using the LLL reduction method
        # That bound must exceed both 4 and w. (See [AKMRVW].)
        all_LLL_bounds = [4, w]
        all_LLL_bounds += [cx_LLL_bound(SUK, A, prec)]
        if S:
            # only need p-adic bound when S nonempty
            all_LLL_bounds.append(p_adic_LLL_bound(SUK, A, prec))

        # Take the largest of all of the bounds we found
        final_LLL_bound = max(all_LLL_bounds)

        if verbose:
            print("The LLL bound is: ", final_LLL_bound)

        # Use the sieve to more easily find all bounds
        S_unit_solutions = sieve_below_bound(K, list(S), final_LLL_bound, verbose=verbose)

    if not include_exponents:
        S_unit_solutions = [sol[2:] for sol in S_unit_solutions]
    if include_bound:
        return S_unit_solutions, final_LLL_bound
    else:
        return S_unit_solutions


def eq_up_to_order(A, B):
    """
    If ``A`` and ``B`` are lists of four-tuples ``[a0,a1,a2,a3]`` and ``[b0,b1,b2,b3]``,
    check that there is some reordering so that either ``ai=bi`` for all ``i`` or
    ``a0==b1``, ``a1==b0``, ``a2==b3``, ``a3==b2``.

    The entries must be hashable.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import eq_up_to_order
        sage: L = [(1,2,3,4), (5,6,7,8)]
        sage: L1 = [L[1], L[0]]
        sage: L2 = [(2,1,4,3), (6,5,8,7)]
        sage: eq_up_to_order(L, L1)
        True
        sage: eq_up_to_order(L, L2)
        True
        sage: eq_up_to_order(L, [(1,2,4,3), (5,6,8,7)])
        False
    """
    # does not look very optimal
    Adup = set(A + [(a[1],a[0],a[3],a[2]) for a in A])
    Bdup = set(B + [(b[1],b[0],b[3],b[2]) for b in B])
    return Adup == Bdup
