r"""
Automorphism groups of dynamical systems of the projective line

AUTHORS:

- Xander Faber, Michelle Manes, Bianca Viray: algorithm and original code
  "Computing Conjugating Sets and Automorphism Groups of Rational Functions" by
  Xander Faber, Michelle Manes, and Bianca Viray [FMV]_.

- Joao de Faria, Ben Hutz, Bianca Thompson (11-2013): adaptation for inclusion in Sage

- Alexander Galarraga (7-2021): Added helper functions for conjugating set
"""
# ****************************************************************************
#       Copyright (C) 2012
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy, deepcopy
from itertools import permutations, combinations, product

from sage.arith.functions import lcm
from sage.arith.misc import CRT, divisors, gcd, is_square
from sage.combinat.permutation import Arrangements
from sage.combinat.subset import Subsets
from sage.matrix.constructor import matrix
from sage.misc.functional import sqrt
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.parallel.use_fork import p_iter_fork
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.sets.primes import Primes
from sage.sets.set import Set
from sage.structure.element import Matrix

lazy_import('sage.rings.number_field.number_field', 'NumberField')


def automorphism_group_QQ_fixedpoints(rational_function, return_functions=False, iso_type=False):
    r"""
    Compute the automorphism group for ``rational_function`` via the method of
    fixed points.

    ALGORITHM:

    See Algorithm 3 in Faber-Manes-Viray [FMV]_.

    INPUT:

    - ``rational_function`` -- Rational Function defined over `\ZZ` or `\QQ`

    - ``return_functions`` -- boolean value; ``True`` will return elements in
      the automorphism group as linear fractional transformations. ``False``
      will return elements as `PGL_2` matrices.

    - ``iso_type`` -- boolean; ``True`` will cause the classification of the
      finite automorphism group to also be returned

    OUTPUT: list of automorphisms that make up the automorphism group
    of ``rational_function``

    EXAMPLES::

        sage: F.<z> = PolynomialRing(QQ)
        sage: rational_function = (z^2 - 2*z - 2)/(-2*z^2 - 2*z + 1)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_QQ_fixedpoints
        sage: automorphism_group_QQ_fixedpoints(rational_function, True)
        [z, 1/z, -z - 1, -z/(z + 1), (-z - 1)/z, -1/(z + 1)]

    ::

        sage: F.<z> = PolynomialRing(QQ)
        sage: rational_function = (z^2 + 2*z)/(-2*z - 1)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_QQ_fixedpoints
        sage: automorphism_group_QQ_fixedpoints(rational_function)
        [
          [1 0]  [-1 -1]  [-2  0]  [0 2]  [-1 -1]  [ 0 -1]
          [0 1], [ 0  1], [ 2  2], [2 0], [ 1  0], [ 1  1]
        ]

    ::

        sage: F.<z> = PolynomialRing(QQ)
        sage: rational_function = (z^2 - 4*z - 3)/(-3*z^2 - 2*z + 2)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_QQ_fixedpoints
        sage: automorphism_group_QQ_fixedpoints(rational_function, True, True)
        ([z, (-z - 1)/z, -1/(z + 1)], 'Cyclic of order 3')
    """

    if rational_function.parent().is_field():
        K = rational_function.parent()
        R = K.ring()
    else:
        R = rational_function.parent()
        K = R.fraction_field()

    F = R.base_ring()

    if F != QQ and F != ZZ:
        raise TypeError("coefficient ring is not the rational numbers or the integers")

    z = R.gen(0)
    phi = R.fraction_field()(rational_function)

    f = phi.numerator()
    g = phi.denominator()

    #scale f,g so both have integer coefficients
    N = lcm(f.denominator(),g.denominator())
    f = f*N
    g = g*N
    N = gcd(gcd(f.coefficients()), gcd(g.coefficients()))
    f = f/N
    g = g/N

    d = max(f.degree(), g.degree())

    h = f - g*z

    if return_functions:
        elements = [z]
    else:
        elements = [matrix(F, 2, [1,0,0,1])]

    rational_roots = h.roots(multiplicities=False)

    min_poly = 1

    #check if infinity is a fixed point
    if g.degree() < d: #then infinity is a fixed point
        #find elements in W of the form (infinity, y)
        #where W is the set of F-rational points (x,y) such that
        #x is fixed by phi and phi(y)=x
        for T in g.roots(multiplicities=False):
            alpha = T
            zeta = -1
            s = (zeta*z + alpha*(1 - zeta))
            if s(phi(z)) == phi(s(z)):
                if return_functions:
                    elements.append(s)
                else:
                    elements.append(matrix(F, 2, [zeta, alpha*(1-zeta), 0, 1]))

    for S in h.roots():
        min_poly = min_poly*(z - S[0])**(S[1])

        if g.degree() < d: #then infinity is a fixed point so (infinity, S[0])
            alpha = S[0]  # is in Z_(1,1)**2
            zeta = -1
            s = (zeta*z + alpha*(1 - zeta))
            if s(phi(z)) == phi(s(z)):
                if return_functions:
                    elements.append(s)
                else:
                    elements.append(matrix(F, 2, [zeta, alpha*(1-zeta), 0, 1]))

        #now compute points in W
        preimage = f - g*S[0]
        if preimage.degree() < d: #infinity is in W
            zeta = -1
            alpha = S[0]
            s = (zeta*z + alpha*(1 - zeta))
            if s(phi(z)) == phi(s(z)):
                if return_functions:
                    elements.append(s)
                else:
                    elements.append(matrix(F, 2, [zeta, alpha*(1-zeta), 0, 1]))
        for T in preimage.roots(multiplicities=False):
            if T != S[0]:
                zeta = -1
                alpha = S[0]
                beta = T
                s = ( (alpha - zeta*beta)*z - (alpha*beta)*(1 - zeta))/((1 - zeta)*z + (alpha*zeta - beta))
                if s(phi(z)) == phi(s(z)):
                    if return_functions:
                        elements.append(s)
                    else:
                        elements.append(matrix(F, 2,
                            [(alpha - zeta*beta), - (alpha*beta)*(1 - zeta),
                             (1 - zeta), (alpha*zeta - beta)]))

    #first look at rational fixed points
    #Subsets is ok since we just needed unordered pairs
    for S in Subsets(rational_roots, 2):
        zeta = -1
        alpha = S[0]
        beta = S[1]
        s = ( (alpha - zeta*beta)*z - (alpha*beta)*(1 - zeta))/((1 - zeta)*z + (alpha*zeta - beta))
        if s(phi(z)) == phi(s(z)):
            if return_functions:
                elements.append(s)
            else:
                elements.append(matrix(F, 2,
                    [(alpha - zeta*beta), - (alpha*beta)*(1 - zeta),
                     (1 - zeta), (alpha*zeta - beta)]))

    # now consider 2-periodic points
    psi = phi(phi(z))
    f2 = psi.numerator()
    g2 = psi.denominator()
    period2_points = [x for x in (f2 - z*g2).roots(multiplicities=False)
                      if x not in rational_roots]
    for S in Subsets(period2_points, 2):
        zeta = -1
        alpha = S[0]
        beta = S[1]
        s = ( (alpha - zeta*beta)*z - (alpha*beta)*(1 - zeta))/((1 - zeta)*z + (alpha*zeta - beta))
        if s(phi(z)) == phi(s(z)):
            if return_functions:
                elements.append(s)
            else:
                elements.append(matrix(F, 2,
                    [(alpha - zeta*beta), - (alpha*beta)*(1 - zeta),
                     (1 - zeta), (alpha*zeta - beta)]))
    if g2.degree() < f2.degree() and g.degree() == d: #infinity has period 2
        for alpha in period2_points:
            zeta = -1
            s = (zeta*z + alpha*(1 - zeta))
            if s(phi(z)) == phi(s(z)):
                if return_functions:
                    elements.append(s)
                else:
                    elements.append(matrix(F, 2, [zeta, alpha*(1-zeta), 0, 1]))
    factors = (f2 - z*g2).factor()
    L1 = NumberField(z**2 + 1,'i')
    i = L1.gen(0)
    L2 = NumberField(z**2 + 3,'isqrt3')
    isqrt3 = L2.gen(0)
    for psi in factors:
        if psi[0].degree() == 2:
            a = psi[0][2]
            b = psi[0][1]
            c = psi[0][0]
            disc = b**2 - 4*a*c
            s = (-b*z - 2*c)/(2*a*z + b)
            if s(phi(z)) == phi(s(z)):
                if return_functions:
                    elements.append(K(s))
                else:
                    elements.append(matrix(F, 2, [-b,-2*c, 2*a, b]))
            if is_square(-disc): #psi[0] generates Q(i)
                alpha = psi[0].change_ring(L1).roots()[0][0]
                beta = alpha.trace() - alpha
                for zeta in [i, -i]:
                    a = (alpha - zeta*beta)/(1 - zeta)
                    d = (alpha*zeta - beta)/(1 - zeta)
                    if a in F and d in F:
                        a = F(a)
                        d = F(d)
                        b = F(-alpha*beta)
                        s = (a * z + b) / (z + d)
                        if s(phi(z)) == phi(s(z)):
                            if return_functions:
                                elements.append(K(s))
                            else:
                                elements.append(matrix(F, 2, [a,b, 1, d]))
            elif is_square(-3*disc): #psi[0] generates Q(zeta_3)
                alpha = psi[0].change_ring(L2).roots()[0][0]
                beta = alpha.trace() - alpha
                for zeta in [F(1)/F(2)*(1 + isqrt3), F(1)/F(2)*(1 - isqrt3),F(1)/F(2)*(-1 + isqrt3), F(1)/F(2)*(-1 - isqrt3)]:
                    a = (alpha - zeta*beta)/(1 - zeta)
                    d = (alpha*zeta - beta)/(1 - zeta)
                    if a in F and d in F:
                        a = F(a)
                        d = F(d)
                        b = F(-alpha*beta)
                        s = (a * z + b) / (z + d)
                        if s(phi(z)) == phi(s(z)):
                            if return_functions:
                                elements.append(K(s))
                            else:
                                elements.append(matrix(F, 2, [a,b, 1, d]))

    if iso_type:
        return elements, which_group(elements)
    return elements


def height_bound(polynomial):
    r"""
    Compute the maximum height of the coefficients of an automorphism.

    This bounds sets the termination criteria for the Chinese Remainder Theorem step.

    Let `f` be a square-free polynomial with coefficients in `K`
    Let `F` be an automorphism of `\mathbb{P}^1_{Frac(R)}` that permutes the roots of `f`
    This function returns a bound on the height of `F`,
    when viewed as an element of `\mathbb{P}^3`

    In [FMV]_ it is proven that `ht(F) <= 6^{[K:Q]}*M`, where `M` is the Mahler measure of `f`
    M is bounded above by `H(f)`, so we return the floor of `6*H(f)`
    (since `ht(F)` is an integer)

    INPUT:

    - ``polynomial`` -- a univariate polynomial

    OUTPUT: a positive integer

    EXAMPLES::

        sage: R.<z> = PolynomialRing(QQ)
        sage: f = z^3 + 2*z + 6
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import height_bound
        sage: height_bound(f)
        413526
    """
    # first check that polynomial is over QQ or ZZ
    K = polynomial.parent()

    if K.is_field():
        R = K.ring()
    else:
        R = K
    F = R.base_ring()

    if F != QQ and F != ZZ:
        raise TypeError("coefficient ring is not the rational numbers or the integers")

    # scale polynomial so that it has integer coefficients with gcd 1
    # this ensures that H(f) = H_infinity(f)
    f = R(polynomial)
    f = f*f.denominator()
    f = f/(gcd(f.coefficients()))

    # compute the infinite height
    L2norm_sq = sum([a**2 for a in f.coefficients()])

    return (6*(L2norm_sq)**3)

def PGL_repn(rational_function):
    r"""
    Take a linear fraction transformation and represent it as a 2x2 matrix.

    INPUT:

    - ``rational_function`` -- a linear fraction transformation

    OUTPUT: a 2x2 matrix representing ``rational_function``

    EXAMPLES::

        sage: R.<z> = PolynomialRing(QQ)
        sage: f = (2*z-1)/(3-z)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import PGL_repn
        sage: PGL_repn(f)
        [-2  1]
        [ 1 -3]
    """
    if isinstance(rational_function, Matrix):
        return rational_function
    K = rational_function.parent()
    F = K.base_ring()
    if not K.is_field():
        return matrix(F, 2, [rational_function[1], rational_function[0], 0, 1])
    else:
        f = rational_function.numerator()
        g = rational_function.denominator()
        return matrix(F, 2, [f[1], f[0], g[1], g[0]])

def PGL_order(A):
    r"""
    Find the multiplicative order of a linear fractional transformation that
    has a finite order as an element of `PGL_2(R)`.

    ``A`` can be represented either as a rational function or a 2x2 matrix

    INPUT:

    - ``A`` -- a linear fractional transformation

    OUTPUT: a positive integer

    EXAMPLES::

        sage: M = matrix([[0,2], [2,0]])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import PGL_order
        sage: PGL_order(M)
        2

    ::

        sage: R.<x> = PolynomialRing(QQ)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import PGL_order
        sage: PGL_order(-1/x)
        2
    """

    n = 1
    AA = PGL_repn(A)
    B = copy(AA)
    while B[0][0] != B[1][1] or B[0][1] != 0 or B[1][0] != 0:
        n = n + 1
        B = AA*B

    return n

def CRT_helper(automorphisms, moduli):
    r"""
    Lift the given list of automorphisms to `Zmod(M)`.

    Given a list of automorphisms over various `Zmod(p^k)` find a list
    of automorphisms over `Zmod(M)` where `M=\prod p^k` that surjects
    onto every tuple of automorphisms from the various `Zmod(p^k)`.

    INPUT:

    - ``automorphisms`` -- list of lists of automorphisms over various `Zmod(p^k)`

    - ``moduli`` -- list of the various `p^k`

    OUTPUT: list of automorphisms over `Zmod(M)`

    EXAMPLES::

        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import CRT_helper
        sage: CRT_helper([[matrix([[4,0], [0,1]]), matrix([[0,1], [1,0]])]], [5])
        ([
        [4 0]  [0 1]
        [0 1], [1 0]
        ], 5)
    """
    if len(automorphisms) > 2:
        temp, modulus = CRT_helper(
            [automorphisms[i] for i in range(len(automorphisms)) if i != 0],
            [moduli[i] for i in range(len(moduli)) if i != 0])
    elif len(automorphisms) == 2:
        temp = automorphisms[1]
        modulus = moduli[1]
    else:
        return automorphisms[0], moduli[0]

    autos = []
    for B in temp:
        for C in automorphisms[0]:
            A = matrix(Integers(modulus*moduli[0]), 2,
                [CRT(B[0][0].lift(), C[0][0].lift(), modulus, moduli[0]),
                 CRT(B[0][1].lift(), C[0][1].lift(), modulus, moduli[0]),
                 CRT(B[1][0].lift(), C[1][0].lift(), modulus, moduli[0]),
                 CRT(B[1][1].lift(), C[1][1].lift(), modulus, moduli[0])])
            autos.append(A)

    return autos, modulus*moduli[0]

def CRT_automorphisms(automorphisms, order_elts, degree, moduli):
    r"""
    Compute a maximal list of automorphisms over `Zmod(M)`.

    Given a list of automorphisms over various `Zmod(p^k)`, a list of the
    elements orders, an integer degree, and a list of the `p^k` values compute
    a maximal list of automorphisms over `Zmod(M)`, such that for every `j` in ``len(moduli)``,
    each element reduces mod ``moduli[j]`` to one of the elements in ``automorphisms[j]`` that
    has order = ``degree``

    INPUT:

    - ``automorphisms`` -- list of lists of automorphisms over various `Zmod(p^k)`

    - ``order_elts`` -- list of lists of the orders of the elements of ``automorphisms``

    - ``degree`` -- positive integer

    - ``moduli`` -- list of prime powers, i.e., `p^k`

    OUTPUT: list containing a list of automorphisms over `Zmod(M)` and the
    product of the moduli

    EXAMPLES::

        sage: aut = [[matrix([[1,0], [0,1]]), matrix([[0,1], [1,0]])]]
        sage: ords = [[1,2]]
        sage: degree = 2
        sage: mods = [5]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import CRT_automorphisms
        sage: CRT_automorphisms(aut,ords,degree,mods)
        ([
        [0 1]
        [1 0]
        ], 5)
    """
    # restrict to automorphisms of degree `degree`
    degree_d_autos = []
    for j in range(len(automorphisms)):
        L = automorphisms[j]
        degree_d_autos.append(
            [L[i] for i in range(len(L)) if order_elts[j][i] == degree])

    # get list of CRT'ed automorphisms
    return CRT_helper(degree_d_autos, moduli)

def valid_automorphisms(automorphisms_CRT, rational_function, ht_bound, M,
                        return_functions=False):
    r"""
    Check if automorphism mod `p^k` lifts to automorphism over `\ZZ`.

    Checks whether an element that is an automorphism of ``rational_function`` modulo `p^k` for various
    `p` s and `k` s can be lifted to an automorphism over `\ZZ`. It uses the fact that every
    automorphism has height at most ``ht_bound``

    INPUT:

    - ``automorphisms`` -- list of lists of automorphisms over various `Zmod(p^k)`

    - ``rational_function`` -- a one variable rational function

    - ``ht_bound`` -- positive integer

    - ``M`` -- positive integer, a product of prime powers

    - ``return_functions`` -- boolean (default: ``False``)

    OUTPUT: list of automorphisms over `\ZZ`

    EXAMPLES::

        sage: R.<z> = PolynomialRing(QQ)
        sage: F = z^2
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import valid_automorphisms
        sage: valid_automorphisms([matrix(GF(5), [[0,1],[1,0]])], F, 48, 5, True)
        [1/z]
    """
    z = rational_function.parent().gen(0)
    valid_auto = []

    for A in automorphisms_CRT:
        init_lift = [x.lift() for x in A.list()]  # lift coefficients of A
        # multiply lift by appropriate scalar matrices and adjust (mod M)
        # to find an element of minimal height. These will have
        # coefficients in [-M/2, M/2)
        for scalar in M.coprime_integers(M):
            new_lift = [scalar*x - (scalar*x/M).round()*M
                        for x in init_lift]
            g = gcd(new_lift)
            new_lift = [x // g for x in new_lift]
            if all(abs(x) <= ht_bound for x in new_lift):
                a, b, c, d = new_lift
                f = (a*z + b) / (c*z + d)
                if rational_function(f(z)) == f(rational_function(z)):
                    if return_functions:
                        valid_auto.append(f)
                    else:
                        valid_auto.append(matrix(ZZ,2,2,new_lift))
                    break

    return valid_auto


def remove_redundant_automorphisms(automorphisms, order_elts, moduli, integral_autos):
    r"""
    If an element of `Aut_{F_p}` has been lifted to `\QQ`
    remove that element from `Aut_{F_p}`.

    We don't want to attempt to lift that element again unnecessarily.

    INPUT:

    - ``automorphisms`` -- list of lists of automorphisms

    - ``order_elts`` -- list of lists of the orders of the elements of ``automorphisms``

    - ``moduli`` -- list of prime powers

    - ``integral_autos`` -- list of known automorphisms

    OUTPUT: list of automorphisms

    EXAMPLES::

        sage: auts = [[matrix([[1,0],[0,1]]), matrix([[6,0],[0,1]]), matrix([[0,1],[1,0]]),
        ....:          matrix([[6,1],[1,1]]), matrix([[1,1],[1,6]]), matrix([[0,6],[1,0]]),
        ....:          matrix([[1,6],[1,1]]), matrix([[6,6],[1,6]])]]
        sage: ord_elts = [[1, 2, 2, 2, 2, 2, 4, 4]]
        sage: mods = [7]
        sage: R.<x> = PolynomialRing(QQ)
        sage: int_auts = [-1/x]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import remove_redundant_automorphisms
        sage: remove_redundant_automorphisms(auts, ord_elts, mods, int_auts)
        [[
        [1 0]  [6 0]  [0 1]  [6 1]  [1 1]  [1 6]  [6 6]
        [0 1], [0 1], [1 0], [1 1], [1 6], [1 1], [1 6]
        ]]
    """
    to_del = []

    for i in range(len(automorphisms)):
        p = moduli[i]
        to_del_temp = []
        for psi in integral_autos:
            #The return_functions boolean determines if the automorphisms
            #are matrices or linear fractional transformations
            if isinstance(psi, Matrix):
                ppsi = psi.change_ring(GF(p))
                B = [ppsi[0,0], ppsi[0,1], ppsi[1,0], psi[1,1]]
            else:
                ff = psi.numerator().change_ring(GF(p))
                gg = psi.denominator().change_ring(GF(p))
                B = [ff[1],ff[0],gg[1],gg[0]]
            for j in range(len(automorphisms[i])):
                A = automorphisms[i][j]
                M = matrix(GF(p), [B, [A[0][0], A[0][1], A[1][0], A[1][1]]])
                if M.rank() == 1:
                    to_del_temp.append(j)
                    break
        to_del.append(to_del_temp)

    for i in range(len(to_del)):
        to_del[i].sort()
        to_del[i].reverse()
        for j in to_del[i]:
            del automorphisms[i][j]
            del order_elts[i][j]

    return automorphisms


def automorphism_group_QQ_CRT(rational_function, prime_lower_bound=4, return_functions=True, iso_type=False):
    r"""
    Determines the complete group of rational automorphisms (under the conjugation action
    of `PGL(2,\QQ)`) for a rational function of one variable.

    See [FMV]_ for details.

    INPUT:

    - ``rational_function`` -- a rational function of a univariate polynomial ring over `\QQ`

    - ``prime_lower_bound`` -- (default: 4) a positive integer; a lower bound for the primes to use for
      the Chinese Remainder Theorem step

    - ``return_functions`` -- boolean (default: ``True``); ``True`` returns
      linear fractional transformations ``False`` returns elements of `PGL(2,\QQ)`

    - ``iso_type`` -- boolean (default: ``False``); ``True`` returns the
      isomorphism type of the automorphism group

    OUTPUT: a complete list of automorphisms of ``rational_function``

    EXAMPLES::

        sage: R.<z> = PolynomialRing(QQ)
        sage: f = (3*z^2 - 1)/(z^3 - 3*z)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_QQ_CRT
        sage: sorted(automorphism_group_QQ_CRT(f, 4, True))
        [-1/z,
         1/z,
         (-z - 1)/(z - 1),
         (-z + 1)/(z + 1),
         (z - 1)/(z + 1),
         (z + 1)/(z - 1),
         -z,
         z]

    ::

        sage: R.<z> = PolynomialRing(QQ)
        sage: f = (3*z^2 - 1)/(z^3 - 3*z)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_QQ_CRT
        sage: sorted(automorphism_group_QQ_CRT(f, 4, False))
        [
        [-1 -1]  [-1  0]  [-1  1]  [ 0 -1]  [0 1]  [ 1 -1]  [1 0]  [ 1  1]
        [ 1 -1], [ 0  1], [ 1  1], [ 1  0], [1 0], [ 1  1], [0 1], [ 1 -1]
        ]
    """
    if rational_function.parent().is_field():
        K = rational_function.parent()
        R = K.ring()
    else:
        R = rational_function.parent()
        K = R.fraction_field()

    F = R.base_ring()

    if F != QQ and F != ZZ:
        raise TypeError("coefficient ring is not the rational numbers or the integers")

    z = R.gen(0)
    phi = K(rational_function)

    f = phi.numerator()
    g = phi.denominator()

    #scale f,g so both have integer coefficients
    N = lcm(f.denominator(),g.denominator())
    f = f*N
    g = g*N
    N = gcd(gcd(f.coefficients()), gcd(g.coefficients()))
    f = f/N
    g = g/N

    d = max(f.degree(), g.degree())

    if d == 1:
        raise ValueError("rational function has degree 1")

    #badprimes is an integer divisible by every prime p such that either
    #    1) phi has bad reduction at p or
    #    2) the reduction map fails to be injective
    badprimes = (gcd(f[d],g[d])*f.resultant(g)*6)
    #6 is because over Q, Aut(phi) has order dividing 12
    #when generalizing to a number field K, 6 should be replaced with
    # 2*gcd(2*[K:Q] + 1, d^3 - d)

    #Determining the set that is used to obtain the height bound
    h = R(prod(x[0] for x in (R(f - g*z)).factor()))# take minimal polynomial of fixed points
    if h.degree() == 2: #if there are only 2 finite fixed points, take preimage of fixed points
        h = h[2]*f**2 + h[1]*f*g + h[0]*g**2
    elif h.degree() == 1: #if there is just 1 finite fixed point, take preimages under phi^2
        psi = phi(phi(z))
        f2 = psi.numerator()
        g2 = psi.denominator()
        N = lcm(f2.denominator(),g2.denominator())
        f2 = f2*N
        g2 = g2*N
        N = gcd(gcd(f2.coefficients()), gcd(g2.coefficients()))
        f2 = f2/N
        g2 = g2/N
        h = h[1]*f2 + h[0]*g2

    MaxH = height_bound(h)
    congruence = 1
    primes = Primes()
    p = primes.next(ZZ(prime_lower_bound))
    primepowers = []
    automorphisms = []
    orderaut = []
    orderelts = []

    if return_functions:
        elements = [z]
    else:
        elements = [matrix(ZZ, 2, [1,0,0,1])]

    badorders = [1, 12]# order 12 not possible over Q, even though 4 and 6 are

    #over QQ, elts of PGL_2 of finite order can only have order dividing 6 or 4,
    # and the finite subgroups can only be cyclic or dihedral (Beauville) so
    # the only possible groups are C_n, D_2n for n|6 or n|4
    # all of these groups have order dividing 24
    while (congruence < (2*MaxH**2)) and len(elements) < gcd(orderaut + [24]):
        if badprimes % p != 0:  #prime of good reduction
            # compute automorphisms mod p
            phi_p = f.change_ring(GF(p))/g.change_ring(GF(p))
            sorted_automorphisms = automorphism_group_FF(phi_p)
            sorted_automorphisms.sort(key=PGL_order)
            orders = [PGL_order(A) for A in sorted_automorphisms]

            automorphisms.append(sorted_automorphisms)
            orderaut.append(len(automorphisms[-1]))
            orderelts.append(orders)
            primepowers.append(p)

            # check if we already found 8 or 12 automorphisms
            # and the gcd of orders over Fp and 24 is 24
            # or if the gcd is equal to the number of automorphisms we have
            if (len(elements) == gcd(orderaut + [24])) or \
                (gcd(orderaut + [24]) == 24 and
                (len(elements) == 12 or len(elements) == 8)):
                    if iso_type:
                        return elements, which_group(elements)
                    return elements
            else:
                N = gcd(orderaut + [12])  # all orders of elements divide N
                for order in divisors(N):
                    if order in badorders:
                        continue
                    # range over all orders
                    # that are possible over QQ such that we haven't already
                    # found all elements of that order

                    # First count number of elements of particular order
                    numeltsoffixedorder = []
                    for L in orderelts:
                        numeltsoffixedorder.append(L.count(order))
                    numelts = min(numeltsoffixedorder)
                    # Have some elts of fixed order mod p for each p
                    if numelts != 0:
                        # CRT order d elements together and check if
                        # they are an automorphism
                        autos, M = CRT_automorphisms(automorphisms,
                                orderelts, order, primepowers)
                        temp = valid_automorphisms(autos, phi, MaxH, M,
                                            return_functions)
                        elements.extend(temp)

                        if (len(elements) == gcd(orderaut + [24])):
                            #found enough automorphisms
                                if iso_type:
                                    return elements, which_group(elements)
                                return elements
                        elif numelts <= (len(temp)):
                            badorders.append(order)
                            # found all elements of order 'order;
                        elif len(temp) != 0:
                            # found some elements of order 'order'
                            # if an element of Aut_{F_p} has been lifted to QQ
                            # remove that element from Aut_{F_p} so we don't
                            # attempt to lift that element again unnecessarily
                            automorphisms = remove_redundant_automorphisms(automorphisms,
                                orderelts, primepowers, temp)
                            if order == 4: #have some elements of order 4
                                # so possible aut group is Z/4 or D_4
                                badorders.extend([3, 6])
                            elif order == 3 or order == 6:#have some elements of
                                # order 3 or 6 so possible aut groups are Z/3,
                                # D_3, Z/6, or D_6
                                badorders.append(4)
                    else: #no elements of order d in some F_v
                        for m in divisors(N):
                            if m % order == 0:
                                badorders.append(m)
                                #no elements of that order or any order that
                                # is a multiple of it
                if all(order in badorders for order in divisors(N)):
                    #found all elements of every possible order
                        if iso_type:
                            return (elements, which_group(elements))
                        return elements
            congruence = congruence * p

        p = primes.next(p)

    if iso_type:
        return elements, which_group(elements)
    return elements


def automorphism_group_FF(rational_function, absolute=False, iso_type=False, return_functions=False):
    r"""
    This function computes automorphism groups over finite fields.

    ALGORITHM:

    See Algorithm 4 in Faber-Manes-Viray [FMV]_.

    INPUT:

    - ``rational_function`` -- a rational function defined over the fraction field
        of a polynomial ring in one variable with finite field coefficients

    - ``absolute`` -- boolean (default: ``False``); ``True`` returns the
      absolute automorphism group and a field of definition

    - ``iso_type`` -- boolean (default: ``False``); ``True`` returns the
      isomorphism type of the automorphism group

    - ``return_functions`` -- boolean (default: ``False``); ``True`` returns
      linear fractional transformations ``False`` returns elements of `PGL(2)`

    OUTPUT: list of automorphisms of ``rational_function``

    EXAMPLES::

        sage: R.<x> = PolynomialRing(GF(5^2, 't'))
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_FF
        sage: automorphism_group_FF((x^2+x+1)/(x+1))
        [
        [1 0]  [4 3]
        [0 1], [0 1]
        ]

    ::

        sage: R.<x> = PolynomialRing(GF(2^5, 't'))
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_FF
        sage: automorphism_group_FF(x^(5), True, False, True)
        [Univariate Polynomial Ring in w over Finite Field in b of size 2^5, [w, 1/w]]

    ::

        sage: R.<x> = PolynomialRing(GF(2^5, 't'))
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_FF
        sage: automorphism_group_FF(x^(5), False, False, True)
        [x, 1/x]
    """

    if not absolute:
        G = automorphism_group_FF_alg3(rational_function)
    else:
        G = automorphism_group_FF_alg2(rational_function)

    if not return_functions:
        if absolute:
            R = G[1][0].parent()
            if R.is_field():
                R = R.ring()
            G[1] = [matrix(R.base_ring(),[[R(g.numerator())[1],R(g.numerator())[0]],[R(g.denominator())[1],R(g.denominator())[0]]]) for g in G[1]]
        else:
            R = G[0].parent()
            if R.is_field():
                R = R.ring()
            G = [matrix(R.base_ring(),[[R(g.numerator())[1],R(g.numerator())[0]],[R(g.denominator())[1],R(g.denominator())[0]]]) for g in G]

    if not iso_type:
        return G
    elif not absolute:
        return G, which_group(G)
    else:
        return G, which_group(G[1])

def field_descent(sigma, y):
    r"""
    Function for descending an element in a field `E` to a subfield `F`.

    Here `F`, `E` must be finite fields or number fields. This function determines
    the unique image of subfield which is ``y`` by the embedding ``sigma`` if it exists.
    Otherwise returns ``None``.
    This functionality is necessary because Sage does not keep track of subfields.

    INPUT:

    - ``sigma`` -- an embedding sigma: `F` -> `E` of fields

    - ``y`` --an element of the field `E`

    OUTPUT: the unique element of the subfield if it exists, otherwise ``None``

    EXAMPLES::

        sage: R = GF(11^2,'b')
        sage: RR = GF(11)
        sage: s = RR.Hom(R)[0]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import field_descent
        sage: field_descent(s, R(1))
        1
    """
    F = sigma.domain()
    a = F.gen()

    p = F.characteristic()
    r = F.degree()
    if p != 0 and y**(p**r) != y:
        return

    K = F.prime_subfield()
    R = PolynomialRing(K,'X')
    f = R(sigma(a).polynomial().coefficients(sparse=False))
    g = R(y.polynomial().coefficients(sparse=False))

    x = F(0)
    quotient, remainder = g.quo_rem(f)
    if not remainder.is_constant():
        return
    else:
        x = x + F(remainder)

    steps = 1
    while not quotient.is_constant():
        quotient, remainder = quotient.quo_rem(f)
        if not remainder.is_constant():
            return
        else:
            x = x + F(remainder)*a**(steps)
            steps += 1

    return x + F(quotient)*a**(steps)

def rational_function_coefficient_descent(rational_function, sigma, poly_ring):
    r"""
    Function for descending the coefficients of a rational function from field `E`
    to a subfield `F`.

    Here `F`, `E` must be finite fields or number fields.
    It determines the unique rational function in fraction field of
    ``poly_ring`` which is the image of ``rational_function`` by ``sigma``,
    if it exists, and otherwise returns ``None``.

    INPUT:

    - ``rational_function``--a rational function with coefficients in a field `E`

    - ``sigma`` -- a field embedding sigma: `F` -> `E`

    - ``poly_ring`` -- a polynomial ring `R` with coefficients in `F`

    OUTPUT: a rational function with coefficients in the fraction field of ``poly_ring``
    if it exists, and otherwise ``None``

    EXAMPLES::

        sage: T.<z> = PolynomialRing(GF(11^2,'b'))
        sage: S.<y> = PolynomialRing(GF(11))
        sage: s = S.base_ring().hom(T.base_ring())
        sage: f = (3*z^3 - z^2)/(z-1)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import rational_function_coefficient_descent
        sage: rational_function_coefficient_descent(f,s,S)
        (3*y^3 + 10*y^2)/(y + 10)
    """

    if rational_function.parent().is_field():
        S = rational_function.parent().ring()
    else:
        S = rational_function.parent()

    if rational_function == S(0):
        return poly_ring(0)

    num = S(rational_function.numerator())
    denom = S(rational_function.denominator())
    f = num.coefficients()
    fe = num.exponents()
    g = denom.coefficients()
    ge = denom.exponents()
    #force the cancellation of common coefficient factors by scaling by f[-1]
    ff = [ field_descent(sigma, x/f[-1]) for x in f]
    gg = [ field_descent(sigma, x/f[-1]) for x in g]
    if None in ff or None in gg:
        return

    z = poly_ring.gen(0)
    numer = sum(poly_ring(ff[i]) * z**fe[i] for i in range(len(ff)))
    denom = sum(poly_ring(gg[i]) * z**ge[i] for i in range(len(gg)))
    return numer / denom


def rational_function_coerce(rational_function, sigma, S_polys):
    r"""
    Function for coercing a rational function defined over a ring `R` to have
    coefficients in a second ring ``S_polys``.

    The fraction field of polynomial ring ``S_polys`` will contain the new rational function.

    INPUT:

    - ``rational_function`` -- rational function with coefficients in `R`

    - ``sigma`` -- a ring homomorphism sigma: `R` -> ``S_polys``

    - ``S_polys`` -- a polynomial ring

    OUTPUT: a rational function with coefficients in ``S_polys``

    EXAMPLES::

        sage: R.<y> = PolynomialRing(QQ)
        sage: S.<z> = PolynomialRing(ZZ)
        sage: s = S.hom([z],R)
        sage: f = (3*z^2 + 1)/(z^3-1)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import rational_function_coerce
        sage: rational_function_coerce(f,s,R)
        (3*y^2 + 1)/(y^3 - 1)
    """
    if rational_function.parent().is_field():
        R = rational_function.parent().ring()
    else:
        R = rational_function.parent()

    f = R(rational_function.numerator()).coefficients(sparse=False)
    g = R(rational_function.denominator()).coefficients(sparse=False)

    if g == [R(1)]:
        return S_polys([sigma(a) for a in f]) # allows for coercion of polynomials
    else:
        return S_polys([sigma(a) for a in f]) / S_polys([sigma(b) for b in g])

def rational_function_reduce(rational_function):
    r"""
    Force Sage to divide out common factors in numerator and denominator
    of rational function.

    INPUT:

    - ``rational_function`` -- rational function `= F/G` in univariate polynomial ring

    OUTPUT: rational function -- `(F/\gcd(F,G)) / (G/\gcd(F,G))`

    EXAMPLES::

        sage: R.<z> = PolynomialRing(GF(7))
        sage: f = ((z-1)*(z^2+z+1))/((z-1)*(z^3+1))
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import rational_function_reduce
        sage: rational_function_reduce(f)
        (z^2 + z + 1)/(z^3 + 1)
    """
    phi = rational_function
    F = phi.numerator()
    G = phi.denominator()
    comm_factor = gcd(F,G)
    return (F.quo_rem(comm_factor)[0]) / (G.quo_rem(comm_factor)[0])

def three_stable_points(rational_function, invariant_list):
    r"""
    Implementation of Algorithm 1 for automorphism groups from
    Faber-Manes-Viray [FMV]_.

    INPUT:

    - ``rational_function`` -- rational function `\phi` defined over finite
      field `E`

    - ``invariant_list`` -- list of at least `3` points of `\mathbb{P}^1(E)` that
      is stable under `Aut_{\phi}(E)`

    OUTPUT: list of automorphisms

    EXAMPLES::

        sage: R.<z> = PolynomialRing(GF(5^2,'t'))
        sage: f = z^3
        sage: L = [[0,1],[4,1],[1,1],[1,0]]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import three_stable_points
        sage: three_stable_points(f,L)
        [z, 4*z, 1/z, 4/z]
    """
    # define ground field and ambient function field
    if rational_function.parent().is_field():
        K = rational_function.parent()
        R = K.ring()
    else:
        R = rational_function.parent()
        K = R.fraction_field()

    z = R.gen(0)
    phi = K(rational_function)

    T = invariant_list

    automorphisms = []
    for t in permutations(range(len(T)),3):
        a = (T[0][0]*T[1][1]*T[2][1]*T[t[0]][0]*T[t[1]][0]*T[t[2]][1] -
            T[0][0]*T[1][1]*T[2][1]*T[t[0]][0]*T[t[1]][1]*T[t[2]][0] -
            T[0][1]*T[1][0]*T[2][1]*T[t[0]][0]*T[t[1]][0]*T[t[2]][1] +
            T[0][1]*T[1][0]*T[2][1]*T[t[0]][1]*T[t[1]][0]*T[t[2]][0] +
            T[0][1]*T[1][1]*T[2][0]*T[t[0]][0]*T[t[1]][1]*T[t[2]][0] -
            T[0][1]*T[1][1]*T[2][0]*T[t[0]][1]*T[t[1]][0]*T[t[2]][0])

        b = (T[0][0]*T[1][0]*T[2][1]*T[t[0]][0]*T[t[1]][1]*T[t[2]][0] -
            T[0][0]*T[1][0]*T[2][1]*T[t[0]][1]*T[t[1]][0]*T[t[2]][0] -
            T[0][0]*T[1][1]*T[2][0]*T[t[0]][0]*T[t[1]][0] * T[t[2]][1] +
            T[0][0]*T[1][1]*T[2][0]*T[t[0]][1]*T[t[1]][0]*T[t[2]][0] +
            T[0][1]*T[1][0]*T[2][0]*T[t[0]][0]*T[t[1]][0]*T[t[2]][1] -
            T[0][1]*T[1][0]*T[2][0]*T[t[0]][0]*T[t[1]][1]*T[t[2]][0])

        c = (T[0][0]*T[1][1]*T[2][1]*T[t[0]][1]*T[t[1]][0] * T[t[2]][1] -
            T[0][0]*T[1][1]*T[2][1]*T[t[0]][1]*T[t[1]][1]*T[t[2]][0] -
            T[0][1]*T[1][0]*T[2][1]*T[t[0]][0]*T[t[1]][1]*T[t[2]][1] +
            T[0][1]*T[1][0]*T[2][1]*T[t[0]][1]*T[t[1]][1]*T[t[2]][0] +
            T[0][1]*T[1][1]*T[2][0]*T[t[0]][0]*T[t[1]][1]*T[t[2]][1] -
            T[0][1]*T[1][1]*T[2][0]*T[t[0]][1]*T[t[1]][0]*T[t[2]][1])

        d = (T[0][0]*T[1][0]*T[2][1]*T[t[0]][0]*T[t[1]][1]*T[t[2]][1] -
            T[0][0]*T[1][0]*T[2][1]*T[t[0]][1]*T[t[1]][0] * T[t[2]][1] -
            T[0][0]*T[1][1]*T[2][0]*T[t[0]][0]*T[t[1]][1]*T[t[2]][1] +
            T[0][0]*T[1][1]*T[2][0]*T[t[0]][1]*T[t[1]][1]*T[t[2]][0] +
            T[0][1]*T[1][0]*T[2][0]*T[t[0]][1]*T[t[1]][0] * T[t[2]][1] -
            T[0][1]*T[1][0]*T[2][0]*T[t[0]][1]*T[t[1]][1]*T[t[2]][0])

        if a*d - b*c != 0:
            s = K(a*z + b) / K(c*z + d)
            if s(phi(z)) == phi(s(z)) and s not in automorphisms:
                automorphisms.append(s)
    return automorphisms

def automorphism_group_FF_alg2(rational_function):
    r"""
    Implementation of algorithm for determining the absolute automorphism
    group over a finite field, given an invariant set, see [FMV]_.

    INPUT:

    - ``rational_function`` -- a rational function defined over a finite field

    OUTPUT: absolute automorphism group of ``rational_function`` and a ring of definition

    EXAMPLES::

        sage: R.<z> = PolynomialRing(GF(7^2,'t'))
        sage: f = (3*z^3 - z^2)/(z-1)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_FF_alg2
        sage: automorphism_group_FF_alg2(f)
        [Univariate Polynomial Ring in w over Finite Field in b of size 7^2, [w, 5/w]]

    ::

        sage: R.<z> = PolynomialRing(GF(5^3,'t'))
        sage: f = (3456*z^(4))
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_FF_alg2
        sage: automorphism_group_FF_alg2(f)
        [Univariate Polynomial Ring in w over Finite Field in b of size 5^6,
         [w,
          (3*b^5 + 4*b^4 + 3*b^2 + 2*b + 1)*w,
          (2*b^5 + b^4 + 2*b^2 + 3*b + 3)*w,
          1/w,
          (3*b^5 + 4*b^4 + 3*b^2 + 2*b + 1)/w,
          (2*b^5 + b^4 + 2*b^2 + 3*b + 3)/w]]
    """
    # define ground field and ambient function field
    if rational_function.parent().is_field():
        K = rational_function.parent()
        R = K.ring()
    else:
        R = rational_function.parent()
        K = R.fraction_field()

    F = R.base_ring()
    if not F.is_finite() or not F.is_field():
        raise TypeError("coefficient ring is not a finite field")
    p = F.characteristic()
    z = R.gen(0)
    phi = K(rational_function)
    f = phi.numerator()
    g = phi.denominator()
    D = max(f.degree(), g.degree())

    # Build an invariant set for phi
    fix = f(z) - z*g(z)
    factor_list = fix.factor()
    minimal_fix_poly = R(prod(x[0] for x in factor_list))
    n = sum(x[0].degree() for x in factor_list) + bool(fix.degree() < D+1)

    if n >= 3:
        T_poly = minimal_fix_poly
        infinity_check = bool(fix.degree() < D+1)
    elif n == 2:
        # Infinity is a fixed point
        if bool(fix.degree() < D+1):
            y = fix.roots(multiplicities=False)[0]
            preimage = g*(f(z) - y*g(z))
            infinity_check = 1
        # Infinity is not a fixed point
        else:
            C = minimal_fix_poly.coefficients(sparse=False)
            preimage = C[2]*f(z)**2 + C[1]*f(z)*g(z) + C[0]*g(z)**2
            infinity_check = bool(preimage.degree() < 2*D)

        T_poly = R(prod(x[0] for x in preimage.factor()))

    else: #case n=1
        # Infinity is the fixed point
        if bool(fix.degree() < D+1):
            minimal_preimage = R(prod(x[0] for x in g.factor()))
            if minimal_preimage.degree() + 1 >= 3:
                T_poly = minimal_preimage
                infinity_check = 1
            else:
                T_poly = R(prod(x[0] for x in phi(phi(z)).denominator().factor() ) )
                infinity_check = 1

        # Infinity is not a fixed point
        else:
            y = fix.roots(multiplicities=False)[0]
            preimage = R(f(z) - y*g(z))
            minimal_preimage = R(prod(x[0] for x in preimage.factor()))
            if minimal_preimage.degree() + bool(preimage.degree() < D) >= 3:
                T_poly = minimal_preimage
                infinity_check = bool(preimage.degree() < D)
            else:
                preimage2 = R(phi(phi(z)).numerator() - y*phi(phi(z)).denominator())
                T_poly = R(prod(x[0] for x in preimage2.factor() ) )
                infinity_check = bool(preimage2.degree() < D**2)

    # Define a field of definition for the absolute automorphism group
    r = lcm([x[0].degree() for x in T_poly.factor()])*F.degree()
    E = GF(p**r,'b')
    sigma = F.Hom(E)[0]
    S = PolynomialRing(E, 'w')
    E_poly = rational_function_coerce(T_poly, sigma, S)

    T = [ [alpha, E(1)] for alpha in E_poly.roots(ring=E, multiplicities=False)]
    if infinity_check == 1:
        T.append([E(1),E(0)])

    # Coerce phi into the larger ring and call Algorithm 1
    Phi = rational_function_coerce(phi, sigma, S)
    return [S, three_stable_points(Phi, T)]


def order_p_automorphisms(rational_function, pre_image):
    r"""
    Determine the order-p automorphisms given the input data.

    This is algorithm 4 in Faber-Manes-Viray [FMV]_.

    INPUT:

    - ``rational_function`` -- rational function defined over finite field `F`

    - ``pre_image`` -- set of triples `[x, L, f]`, where `x` is an `F`-rational
        fixed point of ``rational_function``, `L` is the list of `F`-rational
        pre-images of `x` (excluding `x`), and `f` is the polynomial defining
        the full set of pre-images of `x` (again excluding `x` itself)

    OUTPUT: set of automorphisms of order `p` defined over `F`

    EXAMPLES::

        sage: R.<x> = PolynomialRing(GF(11))
        sage: f = x^11
        sage: L = [[[0, 1], [], 1], [[10, 1], [], 1], [[9, 1], [], 1],
        ....: [[8, 1], [],1], [[7, 1], [], 1], [[6, 1], [], 1], [[5, 1], [], 1],
        ....: [[4, 1], [], 1],[[3, 1], [], 1], [[2, 1], [], 1], [[1, 1], [], 1],
        ....: [[1, 0], [], 1]]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import order_p_automorphisms
        sage: order_p_automorphisms(f,L)
        [x/(x + 1), 6*x/(x + 6), 3*x/(x + 3), 7*x/(x + 7), 9*x/(x + 9), 10*x/(x
        + 10), 5*x/(x + 5), 8*x/(x + 8), 4*x/(x + 4), 2*x/(x + 2), 10/(x + 2),
        (5*x + 10)/(x + 7), (2*x + 10)/(x + 4), (6*x + 10)/(x + 8), (8*x +
        10)/(x + 10), (9*x + 10)/x, (4*x + 10)/(x + 6), (7*x + 10)/(x + 9), (3*x
        + 10)/(x + 5), (x + 10)/(x + 3), (10*x + 7)/(x + 3), (4*x + 7)/(x + 8),
        (x + 7)/(x + 5), (5*x + 7)/(x + 9), (7*x + 7)/x, (8*x + 7)/(x + 1), (3*x
        + 7)/(x + 7), (6*x + 7)/(x + 10), (2*x + 7)/(x + 6), 7/(x + 4), (9*x +
        2)/(x + 4), (3*x + 2)/(x + 9), 2/(x + 6), (4*x + 2)/(x + 10), (6*x +
        2)/(x + 1), (7*x + 2)/(x + 2), (2*x + 2)/(x + 8), (5*x + 2)/x, (x +
        2)/(x + 7), (10*x + 2)/(x + 5), (8*x + 6)/(x + 5), (2*x + 6)/(x + 10),
        (10*x + 6)/(x + 7), (3*x + 6)/x, (5*x + 6)/(x + 2), (6*x + 6)/(x + 3),
        (x + 6)/(x + 9), (4*x + 6)/(x + 1), 6/(x + 8), (9*x + 6)/(x + 6), (7*x +
        8)/(x + 6), (x + 8)/x, (9*x + 8)/(x + 8), (2*x + 8)/(x + 1), (4*x +
        8)/(x + 3), (5*x + 8)/(x + 4), 8/(x + 10), (3*x + 8)/(x + 2), (10*x +
        8)/(x + 9), (8*x + 8)/(x + 7), (6*x + 8)/(x + 7), 8/(x + 1), (8*x +
        8)/(x + 9), (x + 8)/(x + 2), (3*x + 8)/(x + 4), (4*x + 8)/(x + 5), (10*x
        + 8)/x, (2*x + 8)/(x + 3), (9*x + 8)/(x + 10), (7*x + 8)/(x + 8), (5*x +
        6)/(x + 8), (10*x + 6)/(x + 2), (7*x + 6)/(x + 10), 6/(x + 3), (2*x +
        6)/(x + 5), (3*x + 6)/(x + 6), (9*x + 6)/(x + 1), (x + 6)/(x + 4), (8*x
        + 6)/x, (6*x + 6)/(x + 9), (4*x + 2)/(x + 9), (9*x + 2)/(x + 3), (6*x +
        2)/x, (10*x + 2)/(x + 4), (x + 2)/(x + 6), (2*x + 2)/(x + 7), (8*x +
        2)/(x + 2), 2/(x + 5), (7*x + 2)/(x + 1), (5*x + 2)/(x + 10), (3*x +
        7)/(x + 10), (8*x + 7)/(x + 4), (5*x + 7)/(x + 1), (9*x + 7)/(x + 5),
        7/(x + 7), (x + 7)/(x + 8), (7*x + 7)/(x + 3), (10*x + 7)/(x + 6), (6*x
        + 7)/(x + 2), (4*x + 7)/x, (2*x + 10)/x, (7*x + 10)/(x + 5), (4*x +
        10)/(x + 2), (8*x + 10)/(x + 6), (10*x + 10)/(x + 8), 10/(x + 9), (6*x +
        10)/(x + 4), (9*x + 10)/(x + 7), (5*x + 10)/(x + 3), (3*x + 10)/(x + 1),
        x + 1, x + 2, x + 4, x + 8, x + 5, x + 10, x + 9, x + 7, x + 3, x + 6]
    """
    # define ground field and ambient function field
    if rational_function.parent().is_field():
        K = rational_function.parent()
        R = K.ring()
    else:
        R = rational_function.parent()
        K = R.fraction_field()

    z = R.gen(0)
    phi = K(rational_function)
    F = R.base_ring()
    q = F.cardinality()
    p = F.characteristic()
    r = (q-1) / (p-1) # index of F_p^\times inside F^\times

    # Compute the threshold r2 for determining which algorithm to use
    if len(pre_image) > 1:
        r2 = len(pre_image)
        case = 'fix'
    elif len(pre_image[0][1]) > 0:
        r2 = len(pre_image[0][1])
        case = 'F-pre_images'
    else:
        factor_list = pre_image[0][2].factor()
        r2 = sum(x[0].degree() for x in factor_list)
        # Note that infinity is F-rational, so covered by preceding case
        case = 'all pre_images'

    automorphisms_p = []

    if r2 >= r or r2 == 0:
        # Note that r2 == 0 corresponds to having a unique F-rational fixed point
        # that is totally ramified

        for guy in pre_image:
            pt = guy[0]
            zeta = F.multiplicative_generator()
            alpha = zeta**r

            if pt == [F(1),F(0)]:
                for j in range(r):
                    s = z + zeta**j
                    if s(phi(z)) == phi(s(z)):
                        for i in range(p-1):
                            automorphisms_p.append(z+alpha**i*zeta**j)

            else:
                u = F(1) / (z - pt[0])
                u_inv = pt[0] + F(1)/z
                for j in range(r):
                    s = u_inv( u(z) + zeta**j )
                    if s(phi(z)) == phi(s(z)):
                        for i in range(p-1):
                            automorphisms_p.append(u_inv( u(z) + alpha**i*zeta**j) )

    elif r2 < r:

        if case == 'fix':
            T = [x[0] for x in pre_image]
        elif case == 'F-pre_images':
            T = [x for x in pre_image[0][1]]
        else:
            T = []

        # loop over all F-rational pre-images
        for guy in pre_image:
            pt = guy[0]
            # treat case of multiple F-rational fixed points or
            #     1 F-rational fixed point with F-rational pre-images
            if T != []:
                M = [t for t in T if t != pt]
                m = len(M)
                if pt == [F(1),F(0)]:
                    for i in range(1, m):
                        s = z + M[i][0] - M[0][0]
                        if s(phi(z)) == phi(s(z)):
                            automorphisms_p.append(s)
                else:
                    u = F(1) / (z - pt[0])
                    u_inv = pt[0] + F(1)/z
                    for i in range(1,m):
                        if M[0] == [F(1),F(0)]:
                            uy1 = 0
                        else:
                            uy1 = u(M[0][0])
                        if M[i] == [F(1),F(0)]:
                            uy2 = 0
                        else:
                            uy2 = u(M[i][0])
                        s = u_inv( u(z) + uy2 - uy1 )
                        if s(phi(z)) == phi(s(z)):
                            automorphisms_p.append(s)
            elif not T:
                # create the extension field generated by pre-images of the unique fixed point
                T_poly = pre_image[0][2]
                e = lcm([x[0].degree() for x in T_poly.factor()])*F.degree()
                E = GF(p**e, 'b')
                sigma = F.Hom(E)[0]
                S = PolynomialRing(E, 'w')
                w = S.gen(0)
                E_poly = rational_function_coerce(T_poly, sigma, S)
                # List of roots permuted by elements of order p
                # Since infinity is F-rational, it won't appear in this list
                T = [ [alpha, E(1)] for alpha in E_poly.roots(ring=E, multiplicities=False)]

                # coerce the rational function and fixed point into E
                Phi = rational_function_coerce(phi, sigma, S)
                Pt = [sigma(pt[0]), sigma(pt[1])]

                m = len(T)
                if Pt == [E(1),E(0)]:
                    for i in range(1, m):
                        s = w + T[i][0] - T[0][0]
                        if s(Phi(w)) == Phi(s(w)):
                            automorphisms_p.append(rational_function_coefficient_descent(s, sigma, R))
                else:
                    u = E(1) / (w - Pt[0])
                    u_inv = Pt[0] + E(1)/w
                    for i in range(1,m):
                        uy1 = u(T[0][0])
                        uy2 = u(T[i][0])
                        s = u_inv( u(w) + uy2 - uy1 )
                        if s(Phi(w)) == Phi(s(w)):
                            s = rational_function_reduce(s)
                            automorphisms_p.append(rational_function_coefficient_descent(s,sigma,R))

    return automorphisms_p

def automorphisms_fixing_pair(rational_function, pair, quad):
    r"""
    Compute the set of automorphisms with order prime to the characteristic
    that fix the pair, excluding the identity.

    INPUT:

    - ``rational_function`` -- rational function defined over finite field `E`

    - ``pair`` -- a pair of points of `\mathbb{P}^1(E)`

    - ``quad`` -- boolean; an indicator if this is a quadratic pair of points

    OUTPUT: set of automorphisms with order prime to characteristic defined over `E` that fix
    the pair, excluding the identity

    EXAMPLES::

        sage: R.<z> = PolynomialRing(GF(7^2, 't'))
        sage: f = (z^2 + 5*z + 5)/(5*z^2 + 5*z + 1)
        sage: L = [[4, 1], [2, 1]]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphisms_fixing_pair
        sage: sorted(automorphisms_fixing_pair(f, L, False))
        [6/(z + 1), (6*z + 6)/z]
    """
    # define ground field and ambient function field
    if rational_function.parent().is_field():
        K = rational_function.parent()
        R = K.ring()
    else:
        R = rational_function.parent()
        K = R.fraction_field()

    z = R.gen(0)
    phi = K(rational_function)
    E = R.base_ring()
    f = phi.numerator()
    g = phi.denominator()
    D = max(f.degree(), g.degree())

    #assumes the second coordinate of the point is 1
    if pair[0] == [1,0]:
        u = K(z - pair[1][0])
        u_inv = K(z + pair[1][0])
    elif pair[1] == [1,0]:
        u = K(E(1) / (z - pair[0][0]))
        u_inv = K( (pair[0][0]*z + 1) / z )
    else:
        u = K( (z - pair[1][0]) / (z - pair[0][0]) )
        u_inv = K( (pair[0][0]*z - pair[1][0] ) / (z - 1) )

    automorphisms_prime_to_p = []
    # Quadratic automorphisms have order dividing q+1 and D, D-1, or D+1
    if quad:
        #need sqrt to get the cardinality of the base field and not the
        #degree 2 extension
        q = sqrt(E.cardinality())
        zeta = (E.multiplicative_generator())**(q-1)
        for j in [-1,0,1]:
            g = gcd(q+1, D + j)
            xi = zeta**( (q+1) / g )
            for i in range(1,g):
                s = u_inv(xi**i*u(z))
                if s(phi(z)) == phi(s(z)):
                    automorphisms_prime_to_p.append(rational_function_reduce(s))

    # rational automorphisms have order dividing q-1 and D, D-1, or D+1
    else:
        q = E.cardinality()
        zeta = E.multiplicative_generator()
        for j in [-1,0,1]:
            g = gcd(q-1, D + j)
            xi = zeta**( (q-1) / g )
            for i in range(1,g):
                s = u_inv(xi**i*u(z))
                if s(phi(z)) == phi(s(z)):
                    automorphisms_prime_to_p.append(rational_function_reduce(s))

    return list(set(automorphisms_prime_to_p))

def automorphism_group_FF_alg3(rational_function):
    r"""
    Implementation of Algorithm 3 in the paper by Faber/Manes/Viray [FMV]_
    for computing the automorphism group over a finite field.

    INPUT:

    - ``rational_function`` -- a rational function defined over a finite field `F`

    OUTPUT: list of `F`-rational automorphisms of ``rational_function``

    EXAMPLES::

        sage: R.<z> = PolynomialRing(GF(5^3,'t'))
        sage: f = 3456*z^4
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_FF_alg3
        sage: automorphism_group_FF_alg3(f)
        [z, 1/z]
    """
    # define ground field and ambient function field
    if rational_function.parent().is_field():
        K = rational_function.parent()
        R = K.ring()
    else:
        R = rational_function.parent()
        K = R.fraction_field()

    F = R.base_ring()
    if not F.is_finite() or not F.is_field():
        raise TypeError("coefficient ring is not a finite field")
    p = F.characteristic()
    q = F.cardinality()
    z = R.gen(0)
    phi = K(rational_function)
    f = phi.numerator()
    g = phi.denominator()
    D = max(f.degree(), g.degree())

    # For use in the quadratic extension parts of the algorithm
    E = GF(p**(2 * F.degree()), 'b')
    sigma = F.Hom(E)[0]
    S = PolynomialRing(E, 'w')
    Phi = rational_function_coerce(phi, sigma, S)

    # Compute the set of distinct F-rational and F-quadratic
    # factors of the fixed point polynomial
    fix = R(f(z) - z*g(z))
    linear_fix = gcd(fix, z**q - z)
    quad_temp = fix.quo_rem(linear_fix)[0]
    residual = gcd(quad_temp, z**q - z)
    while residual.degree() > 0:
        quad_temp = quad_temp.quo_rem(residual)[0]
        residual = gcd(quad_temp, z**q - z)
    quadratic_fix = gcd(quad_temp, z**(q**2) - z).factor()

    # Compute the set of distinct F-rational fixed points
    linear_fix_pts = [[ x, F(1)] for x in linear_fix.roots(multiplicities=False)]
    if bool(fix.degree() < D+1):
        linear_fix_pts.append( [F(1),F(0)] )
    n1 = len(linear_fix_pts)

    # Coerce quadratic factors into a quadratic extension
    quad_fix_factors = [ rational_function_coerce(poly[0], sigma, S) for poly in quadratic_fix]
    n2 = 2*len(quad_fix_factors)

    # Collect pre-image data as a list L with entries in the form
    # [fixed point y, F-rational pre-images z != y, polynomial defining the pre-images]
    # Note that we remove the fixed point from its pre-image set and its polynomial
    pre_images = []
    for y in linear_fix_pts:
        if y == [F(1),F(0)]:
            Fpre = [ [x,F(1)] for x in g.roots(multiplicities=False) ]
            pre_images.append([y, Fpre, g])
        else:
            Fpre = [ [x,F(1)] for x in (f - y[0]*g).roots(multiplicities=False) if x != y[0]]
            if y[0] == 0 and f.degree() < g.degree():
                Fpre.append([F(1), F(0)]) # infinity is a pre-image of 0
            elif f.degree() == g.degree() and f.leading_coefficient() == y[0]*g.leading_coefficient():
                Fpre.append([F(1), F(0)]) # infinity is a pre-image of y[0]
            # remove y[0] as a root of pre-image polynomial
            h = (f - y[0]*g).quo_rem(z-y[0])[0]
            h_common = gcd(h, z-y[0])
            while h_common.degree() > 0:
                h = h.quo_rem(z-y[0])[0]
                h_common = gcd(h,z-y[0])
            pre_images.append([y, Fpre, h])

    # Initialize the set of automorphisms to contain the identity
    automorphisms = [R(z)]
    automorphisms_quad = []

    # order p elements
    # An F-rational fixed point has orbit length 1 or p under the action of an element of
    # order p. An F-quadratic fixed point has orbit length p. The set of F-rational
    # pre-images of fixed points decomposes as a union of orbits of length p.
    if n1 % p == 1 and n2 % p == 0 and sum(len(x[1]) for x in pre_images) % p == 0:
        # Compute total number of distinct fixed points as a final check for order p auts
        factor_list = fix.factor()
        n = sum(x[0].degree() for x in factor_list) + bool(fix.degree() < D+1)
        if n % p == 1:
            automorphisms = automorphisms + order_p_automorphisms(phi, pre_images)

    # nontrivial elements with order prime to p #
    # case of 2 F-rational fixed points
    for pt_pair in combinations(linear_fix_pts, 2):
        x = pt_pair[0]
        y = pt_pair[1]
        automorphisms = automorphisms + automorphisms_fixing_pair(phi, [x,y], False)

    # case of 1 F-rational fixed point and an F-rational pre-image
    for y in pre_images:
        for x in y[1]:
            automorphisms = automorphisms + automorphisms_fixing_pair(phi, [x,y[0]], False)

    # case of a pair of quadratic fixed points
    for h in quad_fix_factors:
        quad_fix_pts = [ [x,E(1)] for x in h.roots(multiplicities=False)]
        automorphisms_quad = automorphisms_quad + automorphisms_fixing_pair(Phi, quad_fix_pts, True)

    phi_2 = phi(phi(z))
    f_2 = phi_2.numerator()
    g_2 = phi_2.denominator()

    period_2 = (f_2(z) - z*g_2(z)).quo_rem(fix)[0]
    factor_list_2 = period_2.factor()
    linear_period_2_pts = [[ x, F(1)] for x in period_2.roots(multiplicities=False)]
    if bool(period_2.degree() < D**2-D):
        linear_period_2_pts.append( [F(1),F(0)] )
    quad_period_2_factors = [rational_function_coerce(poly[0], sigma, S) for poly in factor_list_2 if poly[0].degree() == 2]
    # n2 = n1 + 2*len(quad_fix_factors)

    # case of a pair of F-rational period 2 points
    linear_period_2_pairs = []
    while linear_period_2_pts:
        x = linear_period_2_pts.pop(-1)
        if x[1] == 1 and g(x[0]) != 0:
            y = [phi(x[0]), F(1)]
        elif x[1] == 1 or f.degree() > g.degree():
            y = [F(1), F(0)]
        elif f.degree() == g.degree():
            y = [f.leading_coefficient() / g.leading_coefficient(), F(1)]
        else:
            y = [F(0), F(1)]

        if x != y:
            linear_period_2_pts.remove(y)
            linear_period_2_pairs.append([x,y])

    for pt_pair in linear_period_2_pairs:
        automorphisms = automorphisms + automorphisms_fixing_pair(phi, pt_pair, False)

    # case of a pair of quadratic period 2 points
    for h in quad_period_2_factors:
        pt_pair = [ [x,E(1)] for x in h.roots(multiplicities=False)]
        if Phi(pt_pair[0][0]) == pt_pair[1][0]:
            automorphisms_quad = automorphisms_quad + automorphisms_fixing_pair(Phi, pt_pair, True)

    # Descend coefficients of the quadratic guys back to the base field
    for s in automorphisms_quad:
        automorphisms.append(rational_function_coefficient_descent(s, sigma, R))

    return automorphisms


def which_group(list_of_elements):
    r"""
    Given a finite subgroup of `PGL_2` determine its isomorphism class.

    This function makes heavy use of the classification of finite subgroups of `PGL(2,K)`.

    INPUT:

    - ``list_of_elements`` -- a finite list of elements of `PGL(2,K)`
      that we know a priori form a group

    OUTPUT: string; the isomorphism type of the group

    EXAMPLES::

        sage: R.<x> = PolynomialRing(GF(7,'t'))
        sage: G = [x, 6*x/(x + 1), 6*x + 6, 1/x, (6*x + 6)/x, 6/(x + 1)]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import which_group
        sage: which_group(G)
        'Dihedral of order 6'
    """
    if isinstance(list_of_elements[-1], Matrix):
        R = PolynomialRing(list_of_elements[-1].base_ring(),'z')
        z = R.gen(0)
        G = [(t[0,0]*z+t[0,1])/(t[1,0]*z+t[1,1]) for t in list_of_elements]
    else:
        G = list_of_elements

    n = ZZ(len(G))

    # invalid input
    if n == 0:
        raise ValueError("group must have at least one element")

    # define ground field and ambient function field
    rational_function = G[-1]

    if rational_function.parent().is_field():
        K = rational_function.parent()
        R = K.ring()
    else:
        R = rational_function.parent()
        K = R.fraction_field()

    z = R.gen(0)
    p = K.characteristic()

    # factor n = mp^e; set e = 0 and m = n if p = 0 (Sage sets 0^0 = 1)
    if p > 0:
        m = n.prime_to_m_part(p)
        e = ZZ(n/m).exact_log(p)
    else:
        m = n
        e = 0

    # Determine if G is cyclic or dihedral.
    # This determines the maximal cyclic subgroup and the maximal cyclic
    # p-regular subgroup. Algorithm terminates if the order of this subgroup agrees with
    # the order of the group.
    max_reg_cyclic = [1, z, [z]]    # initialize order of cyclic p-regular subgroup and generator
    discard = []    # list of elements already considered

    for g in G:
        if g not in discard:
            H = [g]
            for i in range(n-1):
                h = g(H[-1])
                H.append(h)
            H = list(set(H))
            if len(H) == n:
                return 'Cyclic of order {0}'.format(n)
            if len(H) > max_reg_cyclic[0] and gcd(len(H), p) != p:
                max_reg_cyclic = [len(H), g, H]
            discard = list(set(discard + H)) # adjoin all new elements to discard

    n_reg = max_reg_cyclic[0]
    # Test for dihedral subgroup. A subgroup of index 2 is always normal, so the
    # presence of a cyclic subgroup H of index 2 indicates the group is either
    # H x Z/2Z or dihedral. The former occurs only if H has order 1 or 2, both of
    # which are dihedral.
    if 2*n_reg == n:
        for g in G:
            if g not in max_reg_cyclic[2]:
                return 'Dihedral of order {0}'.format(n)
    # Check the p-irregular cases. There is overlap in these cases when p^e = 2,
    # which is dihedral and so already dealt with above. By the classification theorem,
    # these are either p-semi-elementary, PGL(2,q), PSL(2,q), or A_5 when p=3. The latter
    # case is already covered by the remaining sporadic cases below.
    if e > 0:
        if n_reg == m: # p-semi-elementary
            return '{0}-semi-elementary of order {1}'.format(p, n)
        if n_reg == m / (p**e - 1) and m == p**(2*e) - 1:    # PGL(2)
            return 'PGL(2,{0})'.format(p**e)
        if n_reg == m / (p**e - 1) and m == (1/2)*(p**(2*e) - 1):    # PSL(2)
            return 'PSL(2,{0})'.format(p**e)

    # Treat sporadic cases
    if n == 12:
        return ['A_4']
    elif n == 24:
        return ['S_4']
    else:
        return ['A_5']


def conjugating_set_initializer(f, g):
    r"""
    Return a conjugation invariant set together with information
    to reduce the combinatorics of checking all possible conjugations.

    This function constructs the invariant pair (``source``, ``possible_targets``)
    necessary for the conjugating set algorithm described in [FMV2014]_.
    Let `f` and `g` be dynamical systems on `\mathbb{P}^n`.
    An invariant pair is a pair of two sets `U`, `V` such that
    `|U| = |V|` and for all `\phi \in PGL` such that `f^\phi = g`,
    `\phi(u) \in V` for all `u \in U`. Invariant pairs can be used
    to determine all conjugations from `f` to `g`. For details
    in the `\mathbb{P}^1` case, see [FMV2014]_.

    Additionally, this function keeps track of multipliers to reduce the combinatorics.
    This information is then passed to ``conjugating_set_helper`` or
    ``is_conjugate_helper``, which check all possible conjugations determined
    by the invariant pair.

    Do not call this function directly, instead use ``f.conjugating_set(g)``.

    INPUT:

    - ``f`` -- a rational function of degree at least 2, and the same
      degree as ``g``

    - ``g`` -- a nonconstant rational function of the same
      degree as ``f``

    OUTPUT:

    A tuple of the form (``source``, ``possible_targets``).

    - ``source`` -- a conjugation invariant set of `n+2` points of the domain of `f`,
      of which no `n+1` are linearly dependent. Used to specify a possible conjugation
      from `f` to `g`.

    - ``possible_targets`` -- list of tuples of the form (``points``, ``repeated``). ``points``
      is a list of ``points`` which are possible targets for point(s) in ``source``. ``repeated``
      specifies how many points in ``source`` have points in ``points`` as their possible target.

    EXAMPLES:

    We check that ``source`` has no `n+1` linearly dependent points, and that
    ``possible_targets`` tracks multiplier information::

        sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
        sage: f = DynamicalSystem([
        ....:         8*x^7 - 35*x^4*y^3 - 35*x^4*z^3 - 7*x*y^6 - 140*x*y^3*z^3 - 7*x*z^6,
        ....:         -7*x^6*y - 35*x^3*y^4 - 140*x^3*y*z^3 + 8*y^7 - 35*y^4*z^3 - 7*y*z^6,
        ....:         -7*x^6*z - 140*x^3*y^3*z - 35*x^3*z^4 - 7*y^6*z - 35*y^3*z^4 + 8*z^7])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import conjugating_set_initializer
        sage: source, possible_targets = conjugating_set_initializer(f, f)
        sage: P.is_linearly_independent(source, 3)
        True
        sage: f.multiplier(possible_targets[0][0][0], 1) == f.multiplier(source[0], 1)
        True
    """
    n = f.domain().dimension_relative()

    L = Set(f.periodic_points(1))
    K = Set(g.periodic_points(1))
    P = f.codomain().ambient_space()
    if len(L) != len(K):  # checks maps have the same number of fixed points
        return []

    # we store fixed points an multipliers in dictionaries
    # to avoid recalculating them
    mult_to_point_L = {}
    mult_to_point_K = {}
    point_to_mult_L = {}
    point_to_mult_K = {}

    # as we will calculate preimages, we differentiate points by their 'level'
    # which is how many preimages of a fixed point they are, i.e. a fixed point
    # has level 0, a preimage of a fixed point level 1, etc.
    level = 0

    # initializing the dictionaries
    for i in range(len(L)):
        mult_L = f.multiplier(L[i], 1).charpoly()
        mult_K = g.multiplier(K[i], 1).charpoly()
        tup_L = (mult_L, level)
        tup_K = (mult_K, level)
        if tup_L not in mult_to_point_L:
            mult_to_point_L[tup_L] = [L[i]]
        else:
            mult_to_point_L[tup_L] += [L[i]]
        if tup_K not in mult_to_point_K:
            mult_to_point_K[tup_K] = [K[i]]
        else:
            mult_to_point_K[tup_K] += [K[i]]
        point_to_mult_L[L[i]] = (mult_L, level)
        point_to_mult_K[K[i]] = (mult_K, level)

    # we keep a dictionary which tracks how often a (multiplier, level) pair
    # is repeated. As points can only be sent to points with the same (multiplier, level)
    # pair, the less times a (multiplier, level) pair are repeated the better the
    # combinatorics
    repeated_mult_L = {}
    for mult_L in mult_to_point_L:
        repeated = len(mult_to_point_L[mult_L])
        if mult_L not in mult_to_point_K:
            return []
        elif len(mult_to_point_K[mult_L]) != repeated:
            return []
        if repeated not in repeated_mult_L:
            repeated_mult_L[repeated] = [mult_to_point_L[mult_L]]
        else:
            repeated_mult_L[repeated] += [mult_to_point_L[mult_L]]
    more = True

    # the n+2 points to be used to specify PGL conjugations
    source = []

    # a list of tuples of the form ((multiplier, level), repeat) where the
    # (multiplier, level) pair specifies the possible targets of a point in source and repeat
    # specifies how many points in source have that (multiplier, level) pair
    corresponding = []

    # we now greedily look for a set of n+2 points, of which no n+1 are linearly dependent,
    # and we make sure to add the points with the best combinatorics first.
    # this check sometimes fails, i.e. sometimes there is a subset with the
    # desired property which is not found. however, this check is very fast and if it
    # does find a subset, then the subset will most likely minimize the combinatorics
    # of checking conjugations
    tup = greedy_independence_check(P, repeated_mult_L, point_to_mult_L)
    if tup is not None:
        more = False
        source, corresponding = tup

    else:
        # loop_repeated_mult stores the points to find preimages of
        loop_repeated_mult = deepcopy(repeated_mult_L)
        # next_repeated_mult stores the points to find preimages of on the next loop
        next_repeated_mult = {}
        found_no_more = True

    # if we don't find enough points, we go to preimages
    while more:
        level += 1
        # we calculate preimages, starting with preimages with the best
        # expected combinatorics
        for r in sorted(loop_repeated_mult.keys()):
            for point_lst_L in loop_repeated_mult[r]:
                old_tup_L = point_to_mult_L[point_lst_L[0]]
                point_lst_K = mult_to_point_K[old_tup_L]
                mult_L = old_tup_L[0]
                Tl = []
                Tk = []
                # first we calculate preimages
                for pnt in point_lst_L:
                    for preimage in f.rational_preimages(pnt):
                        if preimage != pnt:
                            Tl.append(preimage)
                for pnt in point_lst_K:
                    for preimage in g.rational_preimages(pnt):
                        if preimage != pnt:
                            Tk.append(preimage)
                if len(Tl) != len(Tk):
                    return []
                if Tl:
                    found_no_more = False
                    new_tup_L = (mult_L, level)
                    new_tup_K = (mult_L, level)
                    # we update dictionaries with the new preimages
                    mult_to_point_L[new_tup_L] = Tl
                    mult_to_point_K[new_tup_K] = Tk
                    for i in range(len(Tl)):
                        point_to_mult_L[Tl[i]] = new_tup_L
                        point_to_mult_K[Tk[i]] = new_tup_K
                    repeated = len(Tl)
                    if repeated not in repeated_mult_L:
                        repeated_mult_L[repeated] = [Tl]
                    else:
                        repeated_mult_L[repeated] += [Tl]
                    if repeated not in next_repeated_mult:
                        next_repeated_mult[repeated] = [Tl]
                    else:
                        next_repeated_mult[repeated] += [Tl]
                    # we again do a greedy check for a subset of n+2 points, of which no n+1
                    # are linearly dependent
                    tup = greedy_independence_check(P, repeated_mult_L, point_to_mult_L)
                    if tup is not None:
                        more = False
                        source, corresponding = tup
                if not more:
                    break
            if not more:
                break

        # if no more preimages can be found, we must check all subsets
        # of size n+2 to see if there is a subset in which no n+1 points
        # are linearly dependent
        if found_no_more:
            # we construct a list of all the possible sources points
            all_points = []
            # we order the list by how many repeated multipliers each point has
            # in an attempt to reduce the combinatorics of checking conjugations
            for r in sorted(repeated_mult_L.keys()):
                for point_lst in repeated_mult_L[r]:
                    all_points += point_lst
            # this loop is quite long, so we break after finding the
            # first subset with the desired property. There is,
            # however, no guarantee that the subset we found minimizes
            # the combinatorics when checking conjugations
            for subset in Subsets(range(len(all_points)), n+2):
                source = []
                for i in subset:
                    source.append(all_points[i])
                if P.is_linearly_independent(source, n+1):
                    more = False
                    corresponding = []
                    mult_only = []
                    for i in subset:
                        mult = point_to_mult_L[all_points[i]]
                        if mult in mult_only:
                            corresponding[mult_only.index(mult)][1] += 1
                        else:
                            corresponding.append([mult, 1])
                            mult_only.append(mult)
                    break
            # if we iterated over all subsets of size n+2, and did not find one
            # in which all subsets of size n+1 are linearly independent,
            # then we fail as we cannot specify conjugations
            if more:
                raise ValueError('no more rational preimages; try extending the base field and trying again')

        # if we need to add more preimages, we update loop dictionaries
        if more:
            loop_repeated_mult = deepcopy(next_repeated_mult)
            next_repeated_mult = {}
            found_no_more = True

    # we build a list of iterators in order to loop over the product of those iterators
    possible_targets = []
    for tup in corresponding:
        possible_targets.append([mult_to_point_K[tup[0]], tup[1]])
    return source, possible_targets


def greedy_independence_check(P, repeated_mult, point_to_mult):
    r"""
    Return an invariant pair together with information
    to reduce the combinatorics of checking all possible conjugations.

    Let `f` and `g` be dynamical systems on `\mathbb{P}^n`.
    An invariant pair is a pair of two sets `U`, `V` such that
    `|U| = |V|` and for all `\phi \in PGL` such that `f^\phi = g`,
    `\phi(u) \in V` for all `u \in U`. Invariant pairs can be used
    to determine all conjugations from `f` to `g`. For details
    in the `\mathbb{P}^1` case, see [FMV2014]_.

    This function may sometimes fail to find the invariant pair
    set even though one exists. It is useful, however, as it is fast
    and returns a set which usually minimizes the combinatorics of
    checking all conjugations.

    INPUT:

    - ``P`` -- a projective space

    - ``repeated_mult`` -- dictionary of integers to lists of points of
      the projective space ``P``. The list of points should be conjugation
      invariant. The keys are considered as weights, and this function attempts
      to minimize the total weight

    - ``point_to_mult`` -- dictionary of points of ``P`` to tuples of the form
      (multiplier, level), where multiplier is the characteristic polynomial
      of the multiplier of the point, and level is the number of preimages
      taken to find the point

    OUTPUT:

    If no set of `n+2` points of which all subsets of size `n+1` are linearly
    independent can be found, then ``None`` is returned.

    Otherwise, a tuple of the form (``source``, ``corresponding``) is returned.

    - ``source`` -- the set `U` of the conjugation invariant pair. A set of `n+2` points
      of the domain of `f`, of which no `n+1` are linearly dependent

    - ``corresponding`` -- list of tuples of the form ((multiplier, level), repeat) where the
      (multiplier, level) pair is the multiplier of a point in ``source`` and repeat
      specifies how many points in source have that (multiplier, level) pair. This
      information specifies the set `V` of the invariant pair.

    EXAMPLES::

        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import greedy_independence_check
        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: repeated_mult = {2: [[P((0, 1)), P((1, 0))]], 1: [[P((1, 1))]]}
        sage: point_to_mult = {P((0, 1)): (x, 0), P((1, 0)): (x, 0), P((1, 1)): (x - 2, 0)}
        sage: greedy_independence_check(P, repeated_mult, point_to_mult)
        ([(1 : 1), (0 : 1), (1 : 0)], [[(x - 2, 0), 1], [(x, 0), 2]])
    """
    n = P.dimension_relative()
    source = []
    corresponding = []
    for r in sorted(repeated_mult.keys()):
        for point_lst in repeated_mult[r]:
            for point in point_lst:
                if len(source) == n+1:
                    independent = P.is_linearly_independent(source + [point], n+1)
                else:
                    independent = P.is_linearly_independent(source + [point])
                if independent:
                    source.append(point)
                    mult = point_to_mult[point]
                    # if another point with this multiplier and level pair is in S
                    # then the multiplier level pair will be the last element of corresponding
                    if len(corresponding) != 0:
                        if corresponding[-1][0] == mult:
                            corresponding[-1][1] += 1
                        else:
                            corresponding.append([mult, 1])
                    else:
                        corresponding.append([mult, 1])
                if len(source) == n+2:
                    return source, corresponding

def conjugating_set_helper(f, g, num_cpus, source, possible_targets):
    r"""
    Return the set of elements in PGL over the base ring
    that conjugates ``f`` to ``g``.

    This function takes as input the invariant pair
    and multiplier data from ``conjugating_set_initializer``.

    Do not call this function directly, instead use ``f.conjugate_set(g)``.

    INPUT:

    - ``f`` -- a rational function of degree at least 2, and the same
      degree as ``g``

    - ``g`` -- a rational function of the same degree as ``f``

    - ``num_cpus`` -- the number of threads to run in parallel

    - ``source`` -- list of `n+2` conjugation invariant points, of which
      no `n+1` are linearly dependent

    - ``possible_targets`` -- list of tuples of the form (``points``, ``repeated``). ``points``
      is a list of ``points`` which are possible targets for point(s) in ``source``. ``repeated``
      specifies how many points in ``source`` have points in ``points`` as their possible target.

    OUTPUT: list of elements of PGL which conjugate ``f`` to ``g``

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: f = DynamicalSystem([x^2, y^2])
        sage: source = [P((1, 1)), P((0, 1)), P((1, 0))]
        sage: possible_targets = [[[P((1, 1))], 1], [[P((0, 1)), P((1, 0))], 2]]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import conjugating_set_helper
        sage: sorted(conjugating_set_helper(f, f, 2, source, possible_targets))
        [
        [0 1]  [1 0]
        [1 0], [0 1]
        ]
    """
    Conj = []
    P = f.domain().ambient_space()
    n = f.domain().dimension_relative()

    subset_iterators = []

    for lst in possible_targets:
        subset_iterators.append(Subsets(range(len(lst[0])), lst[1]))

    # helper function for parallelization
    # given a list of tuples which specify indices of possible target
    # points in possible_targets, check all arrangements of those
    # possible target points and if any of them define a conjugation
    # which sends f to g, return those conjugations as a list
    def find_conjugations_subset(tuples):
        conj = []
        for tup in tuples:
            target_set = []
            for i in range(len(tup)):
                for j in tup[i]:
                    target_set.append(possible_targets[i][0][j])

            # if there is a subset of n+1 points which is linearly dependent,
            # we don't need to check any of these arrangements
            if P.is_linearly_independent(target_set, n+1):
                subset_arrangements = []
                for subset in tup:
                    subset_arrangements.append(Arrangements(subset, len(subset)))
                for tup in product(*subset_arrangements):
                    current_target = []
                    for i in range(len(tup)):
                        for j in tup[i]:
                            current_target.append(possible_targets[i][0][j])
                    phi = P.point_transformation_matrix(current_target, source)
                    if f.conjugate(phi) == g:
                        conj.append(phi)
        return conj

    # helper function for parallelization
    # given a list of tuples which specify indices of possible target
    # points in possible_targets, check all possible target points and
    # if any of them define a conjugation which sends f to g, return
    # those conjugations as a list
    def find_conjugations_arrangement(tuples):
        conj = []
        for tup in tuples:
            current_target = []
            for i in range(len(tup)):
                for j in tup[i]:
                    current_target.append(possible_targets[i][0][j])
            phi = P.point_transformation_matrix(current_target, source)
            if f.conjugate(phi) == g:
                conj.append(phi)
        return conj

    if num_cpus > 1:
        all_subsets = list(product(*subset_iterators))
        parallel_data = []

        # if there are enough subsets, we can divide up the work based on subsets
        # and check linear independence in parallel
        if len(all_subsets) > num_cpus:
            for i in range(num_cpus):
                start = (len(all_subsets) * i) // num_cpus
                end = (len(all_subsets) * (i+1)) // num_cpus
                tuples = all_subsets[start:end]
                parallel_data.append(([tuples], {}))

            X = p_iter_fork(num_cpus)
            for ret in X(find_conjugations_subset, parallel_data):
                if ret[1]:
                    Conj += ret[1]
        # otherwise, we need to first check linear independence of the subsets
        # and then build a big list of all the arrangements to split among
        # the threads
        else:
            good_targets = []
            for tup in product(*subset_iterators):
                target_set = []
                for i in range(len(tup)):
                    for j in tup[i]:
                        target_set.append(possible_targets[i][0][j])
                if P.is_linearly_independent(target_set, n+1):
                    good_targets.append(tup)
            all_arrangements = []
            for tup in good_targets:
                subset_arrangements = []
                for subset in tup:
                    subset_arrangements.append(Arrangements(subset, len(subset)))
                all_arrangements += list(product(*subset_arrangements))
            parallel_data = []
            for i in range(num_cpus):
                start = (len(all_arrangements) * i) // num_cpus
                end = (len(all_arrangements) * (i+1)) // num_cpus
                tuples = all_arrangements[start:end]
                parallel_data.append(([tuples], {}))
            X = p_iter_fork(num_cpus)
            for ret in X(find_conjugations_arrangement, parallel_data):
                if ret[1]:
                    Conj += ret[1]
    else:
        Conj = find_conjugations_subset(product(*subset_iterators))
    return Conj

def is_conjugate_helper(f, g, num_cpus, source, possible_targets):
    r"""
    Return if ``f`` is conjugate to ``g``.

    This function takes as input the invariant pair
    and multiplier data from ``conjugating_set_initializer``.

    Do not call this function directly, instead use ``f.is_conjugate(g)``.

    INPUT:

    - ``f`` -- a rational function of degree at least 2, and the same
      degree as ``g``

    - ``g`` -- a rational function of the same degree as ``f``

    - ``num_cpus`` -- the number of threads to run in parallel

    - ``source`` -- list of `n+2` conjugation invariant points, of which
      no `n+1` are linearly dependent

    - ``possible_targets`` -- list of tuples of the form (``points``, ``repeated``). ``points``
      is a list of ``points`` which are possible targets for point(s) in ``source``. ``repeated``
      specifies how many points in ``source`` have points in ``points`` as their possible target.

    OUTPUT: ``True`` if ``f`` is conjugate to ``g``, ``False`` otherwise

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: f = DynamicalSystem([x^2, y^2])
        sage: source = [P((1, 1)), P((0, 1)), P((1, 0))]
        sage: possible_targets = [[[P((1, 1))], 1], [[P((0, 1)), P((1, 0))], 2]]
        sage: from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import is_conjugate_helper
        sage: is_conjugate_helper(f, f, 2, source, possible_targets)
        True
    """
    is_conj = False
    P = f.domain().ambient_space()
    n = f.domain().dimension_relative()

    subset_iterators = []

    for lst in possible_targets:
        subset_iterators.append(Subsets(range(len(lst[0])), lst[1]))

    # helper function for parallelization
    # given a list of tuples which specify indices of possible target
    # points in possible_targets, check all arrangements of those
    # possible target points and if any of them define a conjugation
    # which sends f to g, return True
    def find_conjugations_subset(tuples):
        for tup in tuples:
            target_set = []
            for i in range(len(tup)):
                for j in tup[i]:
                    target_set.append(possible_targets[i][0][j])

            # if there is a subset of n+1 points which is linearly dependent,
            # we don't need to check any of these arrangements
            if P.is_linearly_independent(target_set, n+1):
                subset_arrangements = []
                for subset in tup:
                    subset_arrangements.append(Arrangements(subset, len(subset)))
                for tup in product(*subset_arrangements):
                    current_target = []
                    for i in range(len(tup)):
                        for j in tup[i]:
                            current_target.append(possible_targets[i][0][j])
                    phi = P.point_transformation_matrix(current_target, source)
                    if f.conjugate(phi) == g:
                        return True
        return False

    # helper function for parallelization
    # given a list of tuples which specify indices of possible target points
    # in possible_targets, check all possible target points
    # and if any of them define a conjugation which sends f to g, return True
    def find_conjugations_arrangement(tuples):
        for tup in tuples:
            current_target = []
            for i in range(len(tup)):
                for j in tup[i]:
                    current_target.append(possible_targets[i][0][j])
            phi = P.point_transformation_matrix(current_target, source)
            if f.conjugate(phi) == g:
                return True
        return False

    if num_cpus > 1:
        all_subsets = list(product(*subset_iterators))
        parallel_data = []

        # if there are enough subsets, we can divide up the work based on subsets
        # and check linear independence in parallel
        if len(all_subsets) > num_cpus:
            for i in range(num_cpus):
                start = (len(all_subsets) * i) // num_cpus
                end = (len(all_subsets) * (i+1)) // num_cpus
                tuples = all_subsets[start:end]
                parallel_data.append(([tuples], {}))

            X = p_iter_fork(num_cpus)
            for ret in X(find_conjugations_subset, parallel_data):
                if ret[1]:
                    is_conj = True
                    break
        # otherwise, we need to first check linear independence of the subsets
        # and then build a big list of all the arrangements to split among
        # the threads
        else:
            good_targets = []
            for tup in product(*subset_iterators):
                target_set = []
                for i in range(len(tup)):
                    for j in tup[i]:
                        target_set.append(possible_targets[i][0][j])
                if P.is_linearly_independent(target_set, n+1):
                    good_targets.append(tup)
            all_arrangements = []
            for tup in good_targets:
                subset_arrangements = []
                for subset in tup:
                    subset_arrangements.append(Arrangements(subset, len(subset)))
                all_arrangements += list(product(*subset_arrangements))
            parallel_data = []
            for i in range(num_cpus):
                start = (len(all_arrangements) * i) // num_cpus
                end = (len(all_arrangements) * (i+1)) // num_cpus
                tuples = all_arrangements[start:end]
                parallel_data.append(([tuples], {}))
            X = p_iter_fork(num_cpus)
            for ret in X(find_conjugations_arrangement, parallel_data):
                if ret[1]:
                    is_conj = True
                    break
    else:
        is_conj = find_conjugations_subset(product(*subset_iterators))
    return is_conj
