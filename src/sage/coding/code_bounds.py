# sage_setup: distribution = sagemath-modules
r"""
Bounds for parameters of codes

This module provided some upper and lower bounds for the parameters
of codes.

AUTHORS:

- David Joyner (2006-07): initial implementation.

- William Stein (2006-07): minor editing of docs and code (fixed bug
  in elias_bound_asymp)

- David Joyner (2006-07): fixed dimension_upper_bound to return an
  integer, added example to elias_bound_asymp.

- " (2009-05): removed all calls to Guava but left it as an option.

- Dima Pasechnik (2012-10): added LP bounds.

Let `F` be a finite set of size `q`.
A subset `C` of `V=F^n` is called a code of length `n`.
Often one considers the case where `F` is a finite field,
denoted by `\GF{q}`.  Then `V` is an `F`-vector space.  A subspace
of `V` (with the standard basis) is called a linear code of length `n`. If its
dimension is denoted `k` then we typically store a basis of `C` as a `k\times
n` matrix (the rows are the basis vectors). If `F=\GF{2}` then `C` is called a
binary code. If `F` has `q` elements then `C` is called a `q`-ary code. The
elements of a code `C` are called codewords. The information rate of `C` is


.. MATH::

     R={\frac{\log_q\vert C\vert}{n}},


where `\vert C\vert` denotes the number of elements of `C`. If `{\bf
v}=(v_1,v_2,...,v_n)`, `{\bf w}=(w_1,w_2,...,w_n)` are elements of `V=F^n` then
we define


.. MATH::

     d({\bf v},{\bf w}) =\vert\{i\ \vert\ 1\leq i\leq n,\ v_i\not= w_i\}\vert


to be the Hamming distance between `{\bf v}` and `{\bf w}`. The function
`d:V\times V\rightarrow \Bold{N}` is called the Hamming metric. The weight of
an element (in the Hamming metric) is `d({\bf v},{\bf 0})`,
where `0` is a distinguished element of `F`;
in particular it is `0` of the field if `F` is a field.
The minimum distance of
a linear code is the smallest nonzero weight of a codeword in `C`.  The
relatively minimum distance is denoted


.. MATH::

     \delta = d/n.

A linear code with length `n`, dimension `k`, and minimum distance `d` is
called an `[n,k,d]_q`-code and `n,k,d` are called its parameters. A (not
necessarily linear) code `C` with length `n`, size `M=|C|`, and minimum
distance `d` is called an `(n,M,d)_q`-code (using parentheses instead of square
brackets). Of course, `k=\log_q(M)` for linear codes.

What is the "best" code of a given length?
Let `A_q(n,d)` denote the largest `M` such that there exists a
`(n,M,d)` code in `F^n`. Let `B_q(n,d)` (also denoted `A^{lin}_q(n,d)`) denote
the largest `k` such that there exists a `[n,k,d]` code in `F^n`. (Of course,
`A_q(n,d)\geq B_q(n,d)`.) Determining `A_q(n,d)` and `B_q(n,d)` is one of the
main problems in the theory of error-correcting codes. For more details see
[HP2003]_ and [Lin1999]_.

These quantities related to solving a generalization of the
childhood game of "20 questions".

GAME: Player 1 secretly chooses a number from `1` to
`M` (`M` is large but fixed). Player 2 asks a
series of "yes/no questions" in an attempt to determine that
number. Player 1 may lie at most `e` times
(`e\geq 0` is fixed). What is the minimum number of "yes/no
questions" Player 2 must ask to (always) be able to correctly
determine the number Player 1 chose?

If feedback is not allowed (the only situation considered here),
call this minimum number `g(M,e)`.

Lemma: For fixed `e` and `M`, `g(M,e)` is
the smallest `n` such that `A_2(n,2e+1)\geq M`.

Thus, solving the solving a generalization of the game of "20
questions" is equivalent to determining `A_2(n,d)`! Using
Sage, you can determine the best known estimates for this number in
2 ways:

1. Indirectly, using ``best_known_linear_code_www(n, k, F)``,
   which connects to the website http://www.codetables.de by Markus Grassl;

2. ``codesize_upper_bound(n,d,q)``, ``dimension_upper_bound(n,d,q)``,
   and ``best_known_linear_code(n, k, F)``.

The output of :func:`best_known_linear_code`,
:func:`best_known_linear_code_www`, or :func:`dimension_upper_bound` would
give only special solutions to the GAME because the bounds are applicable
to only linear codes. The output of :func:`codesize_upper_bound` would give
the best possible solution, that may belong to a linear or nonlinear code.

This module implements:

- ``codesize_upper_bound(n,d,q)``, for the best known (as of May, 2006) upper
  bound `A(n,d)` for the size of a code of length `n`, minimum distance `d`
  over a field of size `q`.

- ``dimension_upper_bound(n,d,q)``, an upper bound `B(n,d)=B_q(n,d)` for the
  dimension of a linear code of length `n`, minimum distance `d` over a field
  of size `q`.

- ``gilbert_lower_bound(n,q,d)``, a lower bound for number of elements in the
  largest code of min distance `d` in `\GF{q}^n`.

- ``gv_info_rate(n,delta,q)``, `log_q(GLB)/n`, where GLB is the Gilbert lower
  bound and `\delta = d/n`.

- ``gv_bound_asymp(delta,q)``, asymptotic analog of Gilbert lower bound.

- ``plotkin_upper_bound(n,q,d)``

- ``plotkin_bound_asymp(delta,q)``, asymptotic analog of Plotkin bound.

- ``griesmer_upper_bound(n,q,d)``

- ``elias_upper_bound(n,q,d)``

- ``elias_bound_asymp(delta,q)``, asymptotic analog of Elias bound.

- ``hamming_upper_bound(n,q,d)``

- ``hamming_bound_asymp(delta,q)``, asymptotic analog of Hamming bound.

- ``singleton_upper_bound(n,q,d)``

- ``singleton_bound_asymp(delta,q)``, asymptotic analog of Singleton bound.

- ``mrrw1_bound_asymp(delta,q)``, "first" asymptotic
  McEliese-Rumsey-Rodemich-Welsh bound for the information rate.

-  Delsarte (a.k.a. Linear Programming (LP)) upper bounds.

PROBLEM: In this module we shall typically either (a) seek bounds on `k`, given
`n`, `d`, `q`, (b) seek bounds on `R`, `\delta`, `q` (assuming `n` is
"infinity").

.. TODO::

    - Johnson bounds for binary codes.

    - mrrw2_bound_asymp(delta,q), "second" asymptotic
      McEliese-Rumsey-Rodemich-Welsh bound for the information rate.
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner <wdj@usna.edu>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import binomial, is_prime_power
from sage.features.gap import GapPackage
from sage.misc.functional import sqrt, log
from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF

from .delsarte_bounds import (delsarte_bound_hamming_space,
                              delsarte_bound_additive_hamming_space)

lazy_import('sage.libs.gap.libgap', 'libgap')


def _check_n_q_d(n, q, d, field_based=True):
    r"""
    Check that the length `n`, alphabet size `q` and minimum distance `d` type
    check and make sense for a code over a field.

    More precisely, this checks that the parameters are positive
    integers, that `q` is a prime power for codes over a field, or,
    more generally, that `q` is of size at least 2, and that `n >= d`.
    This raises a :exc:`ValueError` otherwise.

    TESTS::

        sage: # needs sage.libs.pari
        sage: from sage.coding.code_bounds import _check_n_q_d
        sage: _check_n_q_d(20, 16, 5)
        True
        sage: _check_n_q_d(20, 16, 6, field_based=False)
        True
        sage: _check_n_q_d(20, 21, 16)
        Traceback (most recent call last):
        ...
        ValueError: The alphabet size does not make sense for a code over a field
        sage: _check_n_q_d(20, -21, 16)
        Traceback (most recent call last):
        ...
        ValueError: The alphabet size must be an integer >1
        sage: _check_n_q_d(20, 2, 26)
        Traceback (most recent call last):
        ...
        ValueError: The length or minimum distance does not make sense
    """
    if q not in ZZ or q < 2:
        raise ValueError("The alphabet size must be an integer >1")
    if field_based and not is_prime_power(q):
        raise ValueError("The alphabet size does not make sense for a code over a field")
    if not (0 < d <= n and n in ZZ and d in ZZ):
        raise ValueError("The length or minimum distance does not make sense")
    return True


def codesize_upper_bound(n, d, q, algorithm=None):
    r"""
    Return an upper bound on the number of codewords in a (possibly non-linear)
    code.

    This function computes the minimum value of the upper bounds of Singleton,
    Hamming, Plotkin, and Elias.

    If ``algorithm="gap"``, then this returns the best known upper
    bound `A(n,d)=A_q(n,d)` for the size of a code of length `n`,
    minimum distance `d` over a field of size `q`. The function first
    checks for trivial cases (like `d=1` or `n=d`), and if the value
    is in the built-in table. Then it calculates the minimum value
    of the upper bound using the algorithms of Singleton, Hamming,
    Johnson, Plotkin and Elias. If the code is binary,
    `A(n, 2\ell-1) = A(n+1,2\ell)`, so the function
    takes the minimum of the values obtained from all algorithms for the
    parameters `(n, 2\ell-1)` and `(n+1, 2\ell)`. This
    wraps GUAVA's (i.e. GAP's package Guava) ``UpperBound(n, d, q)``.

    If ``algorithm="LP"``, then this returns the Delsarte (a.k.a. Linear
    Programming) upper bound.

    EXAMPLES::

        sage: codes.bounds.codesize_upper_bound(10, 3, 2)
        93
        sage: codes.bounds.codesize_upper_bound(24, 8, 2, algorithm='LP')               # needs sage.numerical.mip
        4096
        sage: codes.bounds.codesize_upper_bound(10, 3, 2, algorithm='gap')              # optional - gap_package_guava
        85
        sage: codes.bounds.codesize_upper_bound(11, 3, 4, algorithm=None)               # needs sage.symbolic
        123361
        sage: codes.bounds.codesize_upper_bound(11, 3, 4, algorithm='gap')              # optional - gap_package_guava
        123361
        sage: codes.bounds.codesize_upper_bound(11, 3, 4, algorithm='LP')               # needs sage.numerical.mip
        109226

    TESTS:

    Make sure :issue:`22961` is fixed::

        sage: codes.bounds.codesize_upper_bound(19, 10, 2)
        20
        sage: codes.bounds.codesize_upper_bound(19, 10, 2, algorithm='gap')            # optional - gap_package_guava
        20

    Meaningless parameters are rejected::

        sage: codes.bounds.codesize_upper_bound(10, -20, 6)
        Traceback (most recent call last):
        ...
        ValueError: The length or minimum distance does not make sense
    """
    _check_n_q_d(n, q, d, field_based=False)
    if algorithm == "gap":
        GapPackage("guava", spkg='gap_packages').require()
        libgap.load_package('guava')
        return int(libgap.UpperBound(n, d, q))
    if algorithm == "LP":
        return int(delsarte_bound_hamming_space(n, d, q))
    else:
        eub = elias_upper_bound(n, q, d)
        hub = hamming_upper_bound(n, q, d)
        pub = plotkin_upper_bound(n, q, d)
        sub = singleton_upper_bound(n, q, d)
        return min([eub, hub, pub, sub])


def dimension_upper_bound(n, d, q, algorithm=None):
    r"""
    Return an upper bound for the dimension of a linear code.

    Return an upper bound `B(n,d) = B_q(n,d)` for the
    dimension of a linear code of length `n`, minimum distance `d` over a
    field of size `q`.

    Parameter ``algorithm`` has the same meaning as in
    :func:`codesize_upper_bound`

    EXAMPLES::

        sage: codes.bounds.dimension_upper_bound(10,3,2)                                # needs sage.libs.pari sage.symbolic
        6
        sage: codes.bounds.dimension_upper_bound(30,15,4)                               # needs sage.libs.pari sage.symbolic
        13
        sage: codes.bounds.dimension_upper_bound(30,15,4,algorithm='LP')                # needs sage.libs.pari sage.numerical.mip
        12

    TESTS:

    Meaningless code parameters are rejected::

        sage: codes.bounds.dimension_upper_bound(13,3,6)                                # needs sage.libs.pari
        Traceback (most recent call last):
        ...
        ValueError: The alphabet size does not make sense for a code over a field
    """
    _check_n_q_d(n, q, d)
    q = ZZ(q)
    if algorithm == "LP":
        return delsarte_bound_additive_hamming_space(n, d, q)
    # algorithm == None or algorithm == "gap":
    return int(ZZ(codesize_upper_bound(n, d, q, algorithm=algorithm)).log(q))


def volume_hamming(n, q, r):
    r"""
    Return the number of elements in a Hamming ball.

    Return the number of elements in a Hamming ball of radius `r` in
    `\GF{q}^n`.

    EXAMPLES::

        sage: codes.bounds.volume_hamming(10,2,3)
        176
    """
    return sum([binomial(n, i) * (q-1)**i
                for i in range(r+1)])


def gilbert_lower_bound(n, q, d):
    r"""
    Return the Gilbert-Varshamov lower bound.

    Return the Gilbert-Varshamov lower bound for number of elements in a largest code of
    minimum distance d in `\GF{q}^n`. See :wikipedia:`Gilbert-Varshamov_bound`

    EXAMPLES::

        sage: codes.bounds.gilbert_lower_bound(10,2,3)
        128/7
    """
    _check_n_q_d(n, q, d, field_based=False)
    ans = q**n/volume_hamming(n,q,d-1)
    return ans


def plotkin_upper_bound(n,q,d, algorithm=None):
    r"""
    Return the Plotkin upper bound.

    Return the Plotkin upper bound for the number of elements in a largest
    code of minimum distance `d` in `\GF{q}^n`.
    More precisely this is a generalization of Plotkin's result for `q=2`
    to bigger `q` due to Berlekamp.

    The ``algorithm="gap"`` option wraps Guava's ``UpperBoundPlotkin``.

    EXAMPLES::

        sage: codes.bounds.plotkin_upper_bound(10,2,3)
        192
        sage: codes.bounds.plotkin_upper_bound(10,2,3,algorithm='gap')  # optional - gap_package_guava
        192
    """
    _check_n_q_d(n, q, d, field_based=False)
    if algorithm == "gap":
        GapPackage("guava", spkg='gap_packages').require()
        libgap.load_package("guava")
        return QQ(libgap.UpperBoundPlotkin(n, d, q))
    else:
        t = 1 - 1/q
        if (q == 2) and (n == 2*d) and (d % 2 == 0):
            return 4*d
        elif (q == 2) and (n == 2*d + 1) and (d % 2 == 1):
            return 4*d + 4
        elif d > t*n:
            return int(d/( d - t*n))
        elif d < t*n + 1:
            fact = (d-1) / t
            from sage.rings.real_mpfr import RR
            if RR(fact) == RR(int(fact)):
                fact = int(fact) + 1
            return int(d/( d - t * fact)) * q**(n - fact)


def griesmer_upper_bound(n,q,d,algorithm=None):
    r"""
    Return the Griesmer upper bound.

    Return the Griesmer upper bound for the number of elements in a
    largest linear code of minimum distance `d` in `\GF{q}^n`, cf. [HP2003]_.
    If the method is "gap", it wraps GAP's ``UpperBoundGriesmer``.

    The bound states:

    .. MATH::

        `n\geq \sum_{i=0}^{k-1} \lceil d/q^i \rceil.`

    EXAMPLES:

    The bound is reached for the ternary Golay codes::

        sage: codes.bounds.griesmer_upper_bound(12,3,6)                                 # needs sage.libs.pari
        729
        sage: codes.bounds.griesmer_upper_bound(11,3,5)                                 # needs sage.libs.pari
        729

    ::

        sage: codes.bounds.griesmer_upper_bound(10,2,3)                                 # needs sage.libs.pari
        128
        sage: codes.bounds.griesmer_upper_bound(10,2,3,algorithm='gap')                 # optional - gap_package_guava, needs sage.libs.pari
        128

    TESTS::

        sage: codes.bounds.griesmer_upper_bound(11,3,6)                                 # needs sage.libs.pari
        243
        sage: codes.bounds.griesmer_upper_bound(11,3,6)                                 # needs sage.libs.pari
        243
    """
    _check_n_q_d(n, q, d)
    if algorithm == "gap":
        GapPackage("guava", spkg='gap_packages').require()
        libgap.load_package("guava")
        return QQ(libgap.UpperBoundGriesmer(n, d, q))
    else:
        # To compute the bound, we keep summing up the terms on the RHS
        # until we start violating the inequality.
        from sage.arith.misc import integer_ceil as ceil
        den = 1
        s = 0
        k = 0
        while s <= n:
            s += ceil(d/den)
            den *= q
            k = k + 1
        return q**(k-1)


def elias_upper_bound(n,q,d,algorithm=None):
    r"""
    Return the Elias upper bound.

    Return the Elias upper bound for number of elements in the largest
    code of minimum distance `d` in `\GF{q}^n`, cf. [HP2003]_.
    If ``algorithm="gap"``, it wraps GAP's ``UpperBoundElias``.

    EXAMPLES::

        sage: codes.bounds.elias_upper_bound(10,2,3)
        232
        sage: codes.bounds.elias_upper_bound(10,2,3,algorithm='gap')  # optional - gap_package_guava
        232
    """
    _check_n_q_d(n, q, d, field_based=False)
    r = 1-1/q
    if algorithm == "gap":
        GapPackage("guava", spkg='gap_packages').require()
        libgap.load_package("guava")
        return QQ(libgap.UpperBoundElias(n, d, q))
    else:
        def ff(n, d, w, q):
            return r*n*d*q**n/((w**2-2*r*n*w+r*n*d)*volume_hamming(n,q,w))

    I = (i for i in range(1, int(r*n) + 1) if i**2 - 2*r*n*i + r*n*d > 0)
    bnd = min([ff(n, d, w, q) for w in I])
    return int(bnd)


def hamming_upper_bound(n,q,d):
    r"""
    Return the Hamming upper bound.

    Return the Hamming upper bound for number of elements in the
    largest code of length `n` and minimum distance `d` over alphabet
    of size `q`.

    The Hamming bound (also known as the sphere packing bound) returns
    an upper bound on the size of a code of length `n`, minimum distance
    `d`, over an alphabet of size `q`. The Hamming bound is obtained by
    dividing the contents of the entire Hamming space
    `q^n` by the contents of a ball with radius
    `floor((d-1)/2)`. As all these balls are disjoint, they can never
    contain more than the whole vector space.


    .. MATH::

         M \leq \frac{q^n}{V(n,e)},


    where `M` is the maximum number of codewords and `V(n,e)` is
    equal to the contents of a ball of radius `e`. This bound is useful
    for small values of `d`. Codes for which equality holds are called
    perfect. See e.g. [HP2003]_.

    EXAMPLES::

        sage: codes.bounds.hamming_upper_bound(10,2,3)
        93
    """
    _check_n_q_d(n, q, d, field_based=False)
    return int((q**n)/(volume_hamming(n, q, int((d-1)/2))))


def singleton_upper_bound(n, q, d):
    r"""
    Return the Singleton upper bound.

    Return the Singleton upper bound for number of elements in a
    largest code of minimum distance `d` in `\GF{q}^n`.

    This bound is based on the shortening of codes. By shortening an
    `(n, M, d)` code `d-1` times, an `(n-d+1,M,1)` code
    results, with `M \leq q^n-d+1`. Thus


    .. MATH::

         M \leq q^{n-d+1}.


    Codes that meet this bound are called maximum distance separable
    (MDS).

    EXAMPLES::

        sage: codes.bounds.singleton_upper_bound(10,2,3)
        256
    """
    _check_n_q_d(n, q, d, field_based=False)
    return q**(n - d + 1)


def gv_info_rate(n, delta, q):
    r"""
    The Gilbert-Varshamov lower bound for information rate.

    The Gilbert-Varshamov lower bound for information rate of a `q`-ary code of
    length `n` and minimum distance `n\delta`.

    EXAMPLES::

        sage: RDF(codes.bounds.gv_info_rate(100,1/4,3))  # abs tol 1e-15                # needs sage.libs.pari sage.symbolic
        0.36704992608261894
    """
    q = ZZ(q)
    return log(gilbert_lower_bound(n,q,int(n*delta)),q)/n


def entropy(x, q=2):
    """
    Compute the entropy at `x` on the `q`-ary symmetric channel.

    INPUT:

    - ``x`` -- real number in the interval `[0, 1]`

    - ``q`` -- (default: 2) integer greater than 1; this is the base of the
      logarithm

    EXAMPLES::

        sage: codes.bounds.entropy(0, 2)
        0
        sage: codes.bounds.entropy(1/5,4).factor()                                      # needs sage.symbolic
        1/10*(log(3) - 4*log(4/5) - log(1/5))/log(2)
        sage: codes.bounds.entropy(1, 3)                                                # needs sage.symbolic
        log(2)/log(3)

    Check that values not within the limits are properly handled::

        sage: codes.bounds.entropy(1.1, 2)
        Traceback (most recent call last):
        ...
        ValueError: The entropy function is defined only for x in the interval [0, 1]
        sage: codes.bounds.entropy(1, 1)
        Traceback (most recent call last):
        ...
        ValueError: The value q must be an integer greater than 1
    """
    if x < 0 or x > 1:
        raise ValueError("The entropy function is defined only for x in the"
                " interval [0, 1]")
    q = ZZ(q)   # This will error out if q is not an integer
    if q < 2:   # Here we check that q is actually at least 2
        raise ValueError("The value q must be an integer greater than 1")
    if x == 0:
        return 0
    if x == 1:
        return log(q-1,q)
    H = x*log(q-1,q)-x*log(x,q)-(1-x)*log(1-x,q)
    return H


def entropy_inverse(x, q=2):
    """
    Find the inverse of the `q`-ary entropy function at the point ``x``.

    INPUT:

    - ``x`` -- real number in the interval `[0, 1]`

    - ``q`` -- (default: 2) integer greater than 1; this is the base of the
      logarithm

    OUTPUT:

    Real number in the interval `[0, 1-1/q]`. The function has multiple
    values if we include the entire interval `[0, 1]`; hence only the
    values in the above interval is returned.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: from sage.coding.code_bounds import entropy_inverse
        sage: entropy_inverse(0.1)                                                      # needs scipy
        0.012986862055...
        sage: entropy_inverse(1)
        1/2
        sage: entropy_inverse(0, 3)
        0
        sage: entropy_inverse(1, 3)
        2/3
    """
    # No nice way to compute the inverse. We resort to root finding.
    if x < 0 or x > 1:
        raise ValueError("The inverse entropy function is defined only for "
                         "x in the interval [0, 1]")
    q = ZZ(q)   # This will error out if q is not an integer
    if q < 2:   # Here we check that q is actually at least 2
        raise ValueError("The value q must be an integer greater than 1")

    eps = 4.5e-16 # find_root has about this as the default xtol
    ymax = 1 - 1/q
    if x <= eps:
        return 0
    if x >= 1-eps:
        return ymax

    # find_root will error out if the root can not be found
    from sage.numerical.optimize import find_root
    f = lambda y: entropy(y, q) - x
    return find_root(f, 0, ymax)


def gv_bound_asymp(delta, q):
    """
    The asymptotic Gilbert-Varshamov bound for the information rate, R.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: RDF(codes.bounds.gv_bound_asymp(1/4,2))                                   # needs sage.libs.pari
        0.18872187554086...
        sage: f = lambda x: codes.bounds.gv_bound_asymp(x,2)
        sage: plot(f,0,1)                                                               # needs sage.libs.pari sage.plot
        Graphics object consisting of 1 graphics primitive
    """
    return 1 - entropy(delta, q)


def hamming_bound_asymp(delta, q):
    """
    The asymptotic Hamming bound for the information rate.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: RDF(codes.bounds.hamming_bound_asymp(1/4,2))                              # needs sage.libs.pari
        0.456435556800...
        sage: f = lambda x: codes.bounds.hamming_bound_asymp(x,2)
        sage: plot(f,0,1)                                                               # needs sage.libs.pari sage.plot
        Graphics object consisting of 1 graphics primitive
    """
    return 1 - entropy(delta / 2, q)


def singleton_bound_asymp(delta, q):
    """
    The asymptotic Singleton bound for the information rate.

    EXAMPLES::

        sage: codes.bounds.singleton_bound_asymp(1/4,2)
        3/4
        sage: f = lambda x: codes.bounds.singleton_bound_asymp(x,2)
        sage: plot(f,0,1)                                                               # needs sage.plot
        Graphics object consisting of 1 graphics primitive
    """
    return 1 - delta


def plotkin_bound_asymp(delta, q):
    r"""
    The asymptotic Plotkin bound for the information rate.

    This only makes sense when `0 < \delta < 1-1/q`.

    EXAMPLES::

        sage: codes.bounds.plotkin_bound_asymp(1/4,2)
        1/2
    """
    r = 1 - 1 / q
    return 1 - delta / r


def elias_bound_asymp(delta, q):
    r"""
    The asymptotic Elias bound for the information rate.

    This only makes sense when `0 < \delta < 1-1/q`.

    EXAMPLES::

        sage: codes.bounds.elias_bound_asymp(1/4,2)                                     # needs sage.symbolic
        0.39912396330...
    """
    r = 1 - 1 / q
    return RDF((1-entropy(r-sqrt(r*(r-delta)), q)))


def mrrw1_bound_asymp(delta, q):
    r"""
    The first asymptotic McEliese-Rumsey-Rodemich-Welsh bound.

    This only makes sense when `0 < \delta < 1-1/q`.

    EXAMPLES::

        sage: codes.bounds.mrrw1_bound_asymp(1/4,2)   # abs tol 4e-16                   # needs sage.symbolic
        0.3545789026652697
    """
    return RDF(entropy((q-1-delta*(q-2)-2*sqrt((q-1)*delta*(1-delta)))/q,q))
