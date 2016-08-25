r"""
Bounds for Parameters of Codes

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

Let `F` be a finite field (we denote the finite field with
`q` elements by `\GF{q}`).
A subset `C` of `V=F^n` is called a code of
length `n`. A subspace of `V` (with the standard
basis) is called a linear code of length `n`. If its
dimension is denoted `k` then we typically store a basis of
`C` as a `k\times n` matrix (the rows are the basis
vectors). If `F=\GF{2}` then `C` is
called a binary code. If `F` has `q` elements
then `C` is called a `q`-ary code. The elements
of a code `C` are called codewords. The information rate
of `C` is


.. math::

     R={\frac{\log_q\vert C\vert}{n}},


where `\vert C\vert` denotes the number of elements of
`C`. If `{\bf v}=(v_1,v_2,...,v_n)`,
`{\bf w}=(w_1,w_2,...,w_n)` are vectors in
`V=F^n` then we define


.. math::

     d({\bf v},{\bf w}) =\vert\{i\ \vert\ 1\leq i\leq n,\ v_i\not= w_i\}\vert


to be the Hamming distance between `{\bf v}` and
`{\bf w}`. The function
`d:V\times V\rightarrow \Bold{N}` is called the Hamming
metric. The weight of a vector (in the Hamming metric) is
`d({\bf v},{\bf 0})`. The minimum distance of a linear
code is the smallest non-zero weight of a codeword in `C`.
The relatively minimum distance is denoted


.. math::

     \delta = d/n.

A linear code with length
`n`, dimension `k`, and minimum distance
`d` is called an `[n,k,d]_q`-code and
`n,k,d` are called its parameters. A (not necessarily
linear) code `C` with length `n`, size
`M=|C|`, and minimum distance `d` is called an
`(n,M,d)_q`-code (using parentheses instead of square
brackets). Of course, `k=\log_q(M)` for linear codes.

What is the "best" code of a given length? Let `F` be a
finite field with `q` elements. Let `A_q(n,d)`
denote the largest `M` such that there exists a
`(n,M,d)` code in `F^n`. Let `B_q(n,d)`
(also denoted `A^{lin}_q(n,d)`) denote the largest
`k` such that there exists a `[n,k,d]` code in
`F^n`. (Of course, `A_q(n,d)\geq B_q(n,d)`.)
Determining `A_q(n,d)` and `B_q(n,d)` is one of
the main problems in the theory of error-correcting codes.

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

1. Indirectly, using best_known_linear_code_www(n, k, F),
    which connects to the website http://www.codetables.de by Markus Grassl;

2. codesize_upper_bound(n,d,q), dimension_upper_bound(n,d,q),
    and best_known_linear_code(n, k, F).

The output of :func:`best_known_linear_code`,
:func:`best_known_linear_code_www`, or :func:`dimension_upper_bound` would
give only special solutions to the GAME because the bounds are applicable
to only linear codes. The output of :func:`codesize_upper_bound` would give
the best possible solution, that may belong to a linear or nonlinear code.

This module implements:

-  codesize_upper_bound(n,d,q), for the best known (as of May,
   2006) upper bound A(n,d) for the size of a code of length n,
   minimum distance d over a field of size q.

-  dimension_upper_bound(n,d,q), an upper bound
   `B(n,d)=B_q(n,d)` for the dimension of a linear code of
   length n, minimum distance d over a field of size q.

-  gilbert_lower_bound(n,q,d), a lower bound for number of
   elements in the largest code of min distance d in
   `\GF{q}^n`.

-  gv_info_rate(n,delta,q), `log_q(GLB)/n`, where GLB is
   the Gilbert lower bound and delta = d/n.

-  gv_bound_asymp(delta,q), asymptotic analog of Gilbert lower
   bound.

-  plotkin_upper_bound(n,q,d)

-  plotkin_bound_asymp(delta,q), asymptotic analog of Plotkin
   bound.

-  griesmer_upper_bound(n,q,d)

-  elias_upper_bound(n,q,d)

-  elias_bound_asymp(delta,q), asymptotic analog of Elias bound.

-  hamming_upper_bound(n,q,d)

-  hamming_bound_asymp(delta,q), asymptotic analog of Hamming
   bound.

-  singleton_upper_bound(n,q,d)

-  singleton_bound_asymp(delta,q), asymptotic analog of Singleton
   bound.

-  mrrw1_bound_asymp(delta,q), "first" asymptotic
   McEliese-Rumsey-Rodemich-Welsh bound for the information rate.

-  Delsarte (a.k.a. Linear Programming (LP)) upper bounds.

PROBLEM: In this module we shall typically either (a) seek bounds
on k, given n, d, q, (b) seek bounds on R, delta, q (assuming n is
"infinity").

TODO:

- Johnson bounds for binary codes.

- mrrw2_bound_asymp(delta,q), "second" asymptotic
  McEliese-Rumsey-Rodemich-Welsh bound for the information rate.

REFERENCES:

- C. Huffman, V. Pless, Fundamentals of error-correcting codes,
  Cambridge Univ. Press, 2003.
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2006 David Joyner <wdj@usna.edu>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.interfaces.all import gap
from sage.rings.all import QQ, RR, ZZ, RDF
from sage.arith.all import factorial
from sage.functions.all import log, sqrt
from sage.misc.decorators import rename_keyword
from .delsarte_bounds import delsarte_bound_hamming_space, \
                delsarte_bound_additive_hamming_space

@rename_keyword(deprecation=6094, method="algorithm")
def codesize_upper_bound(n,d,q,algorithm=None):
    r"""
    This computes the minimum value of the upper bound using the
    methods of Singleton, Hamming, Plotkin, and Elias.

    If algorithm="gap" then this returns the best known upper
    bound `A(n,d)=A_q(n,d)` for the size of a code of length n,
    minimum distance d over a field of size q. The function first
    checks for trivial cases (like d=1 or n=d), and if the value
    is in the built-in table. Then it calculates the minimum value
    of the upper bound using the algorithms of Singleton, Hamming,
    Johnson, Plotkin and Elias. If the code is binary,
    `A(n, 2\ell-1) = A(n+1,2\ell)`, so the function
    takes the minimum of the values obtained from all algorithms for the
    parameters `(n, 2\ell-1)` and `(n+1, 2\ell)`. This
    wraps GUAVA's (i.e. GAP's package Guava) UpperBound( n, d, q ).

    If algorithm="LP" then this returns the Delsarte (a.k.a. Linear
    Programming) upper bound.

    EXAMPLES::

        sage: codes.bounds.codesize_upper_bound(10,3,2)
        93
        sage: codes.bounds.codesize_upper_bound(24,8,2,algorithm="LP")
        4096
        sage: codes.bounds.codesize_upper_bound(10,3,2,algorithm="gap")  # optional - gap_packages (Guava package)
        85
        sage: codes.bounds.codesize_upper_bound(11,3,4,algorithm=None)
        123361
        sage: codes.bounds.codesize_upper_bound(11,3,4,algorithm="gap")  # optional - gap_packages (Guava package)
        123361
        sage: codes.bounds.codesize_upper_bound(11,3,4,algorithm="LP")
        109226

    """
    if algorithm=="gap":
        gap.load_package('guava')
        return int(gap.eval("UpperBound(%s,%s,%s)"%( n, d, q )))
    if algorithm=="LP":
        return int(delsarte_bound_hamming_space(n,d,q))
    else:
        eub = elias_upper_bound(n,q,d)
        gub = griesmer_upper_bound(n,q,d)
        hub = hamming_upper_bound(n,q,d)
        pub = plotkin_upper_bound(n,q,d)
        sub = singleton_upper_bound(n,q,d)
        return min([eub,gub,hub,pub,sub])

@rename_keyword(deprecation=6094, method="algorithm")
def dimension_upper_bound(n,d,q,algorithm=None):
    r"""
    Returns an upper bound `B(n,d) = B_q(n,d)` for the
    dimension of a linear code of length n, minimum distance d over a
    field of size q.
    Parameter "algorithm" has the same meaning as in :func:`codesize_upper_bound`

    EXAMPLES::

        sage: codes.bounds.dimension_upper_bound(10,3,2)
        6
        sage: codes.bounds.dimension_upper_bound(30,15,4)
        13
        sage: codes.bounds.dimension_upper_bound(30,15,4,algorithm="LP")
        12

    """
    q = ZZ(q)
    if algorithm=="LP":
        return delsarte_bound_additive_hamming_space(n,d,q)

    else:       # algorithm==None or algorithm=="gap":
        return int(log(codesize_upper_bound(n,d,q,algorithm=algorithm),q))


def volume_hamming(n,q,r):
    r"""
    Returns number of elements in a Hamming ball of radius r in `\GF{q}^n`.
    Agrees with Guava's SphereContent(n,r,GF(q)).

    EXAMPLES::

        sage: codes.bounds.volume_hamming(10,2,3)
        176
    """
    ans=sum([factorial(n)/(factorial(i)*factorial(n-i))*(q-1)**i for i in range(r+1)])
    return ans

def gilbert_lower_bound(n,q,d):
    r"""
    Returns lower bound for number of elements in the largest code of
    minimum distance d in `\GF{q}^n`.

    EXAMPLES::

        sage: codes.bounds.gilbert_lower_bound(10,2,3)
        128/7
    """
    ans=q**n/volume_hamming(n,q,d-1)
    return ans

@rename_keyword(deprecation=6094, method="algorithm")
def plotkin_upper_bound(n,q,d, algorithm=None):
    r"""
    Returns Plotkin upper bound for number of elements in the largest
    code of minimum distance d in `\GF{q}^n`.

    The algorithm="gap" option wraps Guava's UpperBoundPlotkin.

    EXAMPLES::

        sage: codes.bounds.plotkin_upper_bound(10,2,3)
        192
        sage: codes.bounds.plotkin_upper_bound(10,2,3,algorithm="gap")  # optional - gap_packages (Guava package)
        192
    """
    if algorithm=="gap":
        ans=gap.eval("UpperBoundPlotkin(%s,%s,%s)"%(n,d,q))
        #print "calling Guava ..."
        return QQ(ans)
    else:
        t = 1 - 1/q
        if (q==2) and (n == 2*d) and (d%2 == 0):
            return 4*d
        elif (q==2) and (n == 2*d + 1) and (d%2 == 1):
            return 4*d + 4
        elif d > t*n:
            return int(d/( d - t*n))
        elif d < t*n + 1:
            fact = (d-1) / t
            if RR(fact)==RR(int(fact)):
                fact = int(fact) + 1
            return int(d/( d - t * fact)) * q**(n - fact)

@rename_keyword(deprecation=6094, method="algorithm")
def griesmer_upper_bound(n,q,d,algorithm=None):
    r"""
    Returns the Griesmer upper bound for number of elements in the
    largest code of minimum distance d in `\GF{q}^n`.
    Wraps GAP's UpperBoundGriesmer.

    EXAMPLES::

        sage: codes.bounds.griesmer_upper_bound(10,2,3)
        128
        sage: codes.bounds.griesmer_upper_bound(10,2,3,algorithm="gap")  # optional - gap_packages (Guava package)
        128
    """
    if algorithm=="gap":
        #print "calling Guava ..."
        ans=gap.eval("UpperBoundGriesmer(%s,%s,%s)"%(n,d,q))
        return QQ(ans)
    else:
        den = 1
        s = 0
        k = 0
        add = 0
        while s <= n:
            if not(add == 1):
                if d%den==0:
                    add = int(d/den)
                else:
                    add = int(d/den)+1
            s = s + add
            den = den * q
            k = k + 1
        return q**(k-1)


@rename_keyword(deprecation=6094, method="algorithm")
def elias_upper_bound(n,q,d,algorithm=None):
    r"""
    Returns the Elias upper bound for number of elements in the largest
    code of minimum distance d in `\GF{q}^n`. Wraps
    GAP's UpperBoundElias.

    EXAMPLES::

        sage: codes.bounds.elias_upper_bound(10,2,3)
        232
        sage: codes.bounds.elias_upper_bound(10,2,3,algorithm="gap")  # optional - gap_packages (Guava package)
        232

    """
    r = 1-1/q
    if algorithm=="gap":
        #print "calling Guava ..."
        ans=gap.eval("UpperBoundElias(%s,%s,%s)"%(n,d,q))
        return QQ(ans)
    else:
        def ff(n,d,w,q):
            return r*n*d*q**n/((w**2-2*r*n*w+r*n*d)*volume_hamming(n,q,w));
    def get_list(n,d,q):
        I = []
        for i in range(1,int(r*n)+1):
            if i**2-2*r*n*i+r*n*d>0:
                I.append(i)
        return I
    I = get_list(n,d,q)
    bnd = min([ff(n,d,w,q) for w in I])
    return int(bnd)

def hamming_upper_bound(n,q,d):
    r"""
    Returns the Hamming upper bound for number of elements in the
    largest code of minimum distance d in `\GF{q}^n`.
    Wraps GAP's UpperBoundHamming.

    The Hamming bound (also known as the sphere packing bound) returns
    an upper bound on the size of a code of length n, minimum distance
    d, over a field of size q. The Hamming bound is obtained by
    dividing the contents of the entire space
    `\GF{q}^n` by the contents of a ball with radius
    floor((d-1)/2). As all these balls are disjoint, they can never
    contain more than the whole vector space.


    .. math::

         M \leq {q^n \over V(n,e)},



    where M is the maximum number of codewords and `V(n,e)` is
    equal to the contents of a ball of radius e. This bound is useful
    for small values of d. Codes for which equality holds are called
    perfect.

    EXAMPLES::

        sage: codes.bounds.hamming_upper_bound(10,2,3)
        93
    """
    return int((q**n)/(volume_hamming(n, q, int((d-1)/2))))

def singleton_upper_bound(n,q,d):
    r"""
    Returns the Singleton upper bound for number of elements in the
    largest code of minimum distance d in `\GF{q}^n`.
    Wraps GAP's UpperBoundSingleton.

    This bound is based on the shortening of codes. By shortening an
    `(n, M, d)` code d-1 times, an `(n-d+1,M,1)` code
    results, with `M \leq q^n-d+1`. Thus


    .. math::

         M \leq q^{n-d+1}.



    Codes that meet this bound are called maximum distance separable
    (MDS).

    EXAMPLES::

        sage: codes.bounds.singleton_upper_bound(10,2,3)
        256
    """
    return q**(n - d + 1)

def gv_info_rate(n,delta,q):
    """
    GV lower bound for information rate of a q-ary code of length n
    minimum distance delta\*n

    EXAMPLES::

        sage: RDF(codes.bounds.gv_info_rate(100,1/4,3))  # abs tol 1e-15
        0.36704992608261894
    """
    q = ZZ(q)
    ans=log(gilbert_lower_bound(n,q,int(n*delta)),q)/n
    return ans

def entropy(x, q=2):
    """
    Computes the entropy at `x` on the `q`-ary symmetric channel.

    INPUT:

    - ``x`` - real number in the interval `[0, 1]`.

    - ``q`` - (default: 2) integer greater than 1. This is the base of the
      logarithm.

    EXAMPLES::

        sage: codes.bounds.entropy(0, 2)
        0
        sage: codes.bounds.entropy(1/5,4)
        1/5*log(3)/log(4) - 4/5*log(4/5)/log(4) - 1/5*log(1/5)/log(4)
        sage: codes.bounds.entropy(1, 3)
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
    Find the inverse of the ``q``-ary entropy function at the point ``x``.

    INPUT:

    - ``x`` -- real number in the interval `[0, 1]`.

    - ``q`` - (default: 2) integer greater than 1. This is the base of the
      logarithm.

    OUTPUT:

    Real number in the interval `[0, 1-1/q]`. The function has multiple
    values if we include the entire interval `[0, 1]`; hence only the
    values in the above interval is returned.

    EXAMPLES::

        sage: from sage.coding.code_bounds import entropy_inverse
        sage: entropy_inverse(0.1)
        0.012986862055848683
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

    eps  = 4.5e-16 # find_root has about this as the default xtol
    ymax = 1 - 1/q
    if x <= eps:
        return 0
    if x >= 1-eps:
        return ymax

    # find_root will error out if the root can not be found
    from sage.numerical.optimize import find_root
    f = lambda y: entropy(y, q) - x
    return find_root(f, 0, ymax)

def gv_bound_asymp(delta,q):
    """
    Computes the asymptotic GV bound for the information rate, R.

    EXAMPLES::

        sage: RDF(codes.bounds.gv_bound_asymp(1/4,2))
        0.18872187554086...
        sage: f = lambda x: codes.bounds.gv_bound_asymp(x,2)
        sage: plot(f,0,1)
        Graphics object consisting of 1 graphics primitive
    """
    return (1-entropy(delta,q))


def hamming_bound_asymp(delta,q):
    """
    Computes the asymptotic Hamming bound for the information rate.

    EXAMPLES::

        sage: RDF(codes.bounds.hamming_bound_asymp(1/4,2))
        0.456435556800...
        sage: f = lambda x: codes.bounds.hamming_bound_asymp(x,2)
        sage: plot(f,0,1)
        Graphics object consisting of 1 graphics primitive
    """
    return (1-entropy(delta/2,q))

def singleton_bound_asymp(delta,q):
    """
    Computes the asymptotic Singleton bound for the information rate.

    EXAMPLES::

        sage: codes.bounds.singleton_bound_asymp(1/4,2)
        3/4
        sage: f = lambda x: codes.bounds.singleton_bound_asymp(x,2)
        sage: plot(f,0,1)
        Graphics object consisting of 1 graphics primitive
    """
    return (1-delta)

def plotkin_bound_asymp(delta,q):
    """
    Computes the asymptotic Plotkin bound for the information rate,
    provided `0 < \delta < 1-1/q`.

    EXAMPLES::

        sage: codes.bounds.plotkin_bound_asymp(1/4,2)
        1/2
    """
    r = 1-1/q
    return (1-delta/r)

def elias_bound_asymp(delta,q):
    """
    Computes the asymptotic Elias bound for the information rate,
    provided `0 < \delta < 1-1/q`.

    EXAMPLES::

        sage: codes.bounds.elias_bound_asymp(1/4,2)
        0.39912396330...
    """
    r = 1-1/q
    return RDF((1-entropy(r-sqrt(r*(r-delta)), q)))

def mrrw1_bound_asymp(delta,q):
    """
    Computes the first asymptotic McEliese-Rumsey-Rodemich-Welsh bound
    for the information rate, provided `0 < \delta < 1-1/q`.

    EXAMPLES::

        sage: codes.bounds.mrrw1_bound_asymp(1/4,2)   # abs tol 4e-16
        0.3545789026652697
    """
    return RDF(entropy((q-1-delta*(q-2)-2*sqrt((q-1)*delta*(1-delta)))/q,q))
