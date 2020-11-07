# distutils: language = c
# distutils: libraries = fes
# sage_setup: distribution = sage-fes

"""
Binding for the FES library

Finding solutions of systems of boolean equations by exhaustive
search, via the fes library. This is usually (much) faster than
computing a Groebner basis, except in special cases where the latter
is particularly easy.

The FES library is presently only able to deal with polynomials in 64
variables. Performing a full exhaustive search over 64 variables will
take a **long** time. The number of variables can be artificially
reduced to 64 by specializing some of them.

Note that the FES library **requires** at least of the equations to be
non-linear.

AUTHORS:

- Charles Bouillaguet (2012-12-20) : initial version

EXAMPLES:

Random Degree-2 System::

    sage: from sage.libs.fes import exhaustive_search                            # optional - FES
    sage: n = 16                                                                 # optional - FES
    sage: R = BooleanPolynomialRing(n, 'x')                                      # optional - FES
    sage: solution = dict(zip(R.gens(), VectorSpace(GF(2), n).random_element())) # optional - FES

    sage: F = [ R.random_element() for i in range(n+10)  ]                       # optional - FES
    sage: G = [ f + f.subs(solution) for f in F]                                 # optional - FES
    sage: sols = exhaustive_search(G)                                            # optional - FES
    sage: solution in sols                                                       # optional - FES
    True
    sage: len(sols)                                                              # optional - FES
    1

Cyclic benchmark::

    sage: from sage.rings.ideal import Cyclic                 # optional - FES
    sage: from sage.libs.fes import exhaustive_search         # optional - FES
    sage: n = 10                                              # optional - FES
    sage: R = BooleanPolynomialRing(n, 'x')                   # optional - FES
    sage: sols = exhaustive_search( Cyclic(R) )               # optional - FES
    sage: len(sols)                                           # optional - FES
    1
    sage: set(sols[0]) == set(R.gens())                       # optional - FES
    True

REFERENCES:

- [BCCCNSY2010]_
"""

#*****************************************************************************
#       Copyright (C) 2012 Charles Bouillaguet <charles.bouillaguet@lifl.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.stdint cimport uint64_t
from cysignals.memory cimport sig_calloc, sig_free
from cysignals.signals cimport sig_on, sig_off

cdef extern from "fes_interface.h":
    ctypedef int (*solution_callback_t)(void *, uint64_t)

    void exhaustive_search_wrapper(int n, int n_eqs, int degree, int ***coeffs, solution_callback_t callback, void* callback_state, int verbose)


from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

from sage.structure.parent cimport Parent
from sage.structure.sequence import Sequence
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.polynomial.pbori.pbori import BooleanPolynomial, BooleanPolynomialRing
from sage.arith.all import binomial
from sage.combinat.subset import Subsets

from sage.matrix.all import *
from sage.modules.all import *



class InternalState:
    verbose = False
    sols = []
    max_sols = 0


cdef int report_solution(void *_state, uint64_t i):
#   This is the callback function which is invoked by the fes library
#   each time a solution is found. ``i`` describes the solution (least
#   significant bit gives the value of the first variable, most
#   significant bit gives the value of the last variable). ``state`` is
#   the pointer passed to the fes library initially.
    cdef object state = <object> _state
    state.sols.append(i)
    if state.verbose:
        print("fes: solution {0} / {1} found : {2:x}".format(len(state.sols), state.max_sols, i))
    if (state.max_sols > 0 and state.max_sols >= len(state.sols)):
        return 1  #stop the library
    return 0 # keep going


def exhaustive_search(eqs,  max_sols=Infinity, verbose=False):
    r"""
    Invokes the fes library to solve a system of boolean equation by
    exhaustive search.

    INPUT:

    - ``eqs`` -- list of boolean equations

    - ``max_sols`` -- stop after this many solutions are found

    - ``verbose`` -- whether the library should display information about
      its work

    NOTE:

        Using this function requires the optional package FES to be installed.

    EXAMPLES:

    A very simple example::

        sage: from sage.libs.fes import exhaustive_search               # optional - FES
        sage: R.<x,y,z> = BooleanPolynomialRing()                       # optional - FES
        sage: exhaustive_search( [x*y + z + y + 1, x+y+z, x*y + y*z] )  # optional - FES
        [{y: 0, z: 1, x: 1}]

    Another very simple example::

        sage: R.<x,y,z,t,w> = BooleanPolynomialRing()                   # optional - FES
        sage: f = x*y*z + z + y + 1                                     # optional - FES
        sage: sorted( exhaustive_search( [f, x+y+z, x*y + y*z] ) )      # optional - FES
        [{w: 0, t: 0, y: 0, z: 1, x: 1}, {w: 0, t: 1, y: 0, z: 1, x: 1}, {w: 1, t: 0, y: 0, z: 1, x: 1}, {w: 1, t: 1, y: 0, z: 1, x: 1}]

    An example with several solutions::

        sage: R.<x,y,z, t> = BooleanPolynomialRing()                    # optional - FES
        sage: eqs = [x*z + y*t + 1,   x*y + y*z + z*t + t*x , x+y+z+1]  # optional - FES
        sage: sorted( exhaustive_search( eqs ) )                        # optional - FES
        [{t: 0, y: 1, z: 1, x: 1}, {t: 1, y: 1, z: 0, x: 0}]

    We voluntarily limit the number of solutions returned by the library::

        sage: eqs =  [x*z + y*t + 1,   x*y + y*z + z*t + t*x , x+y+z+1] # optional - FES
        sage: exhaustive_search(eqs, max_sols=1 )                       # optional - FES
        [{t: 0, y: 1, z: 1, x: 1}]

    .. SEEALSO::

        This function should return the same solutions as
        :meth:`~sage.rings.polynomial.pbori.BooleanPolynomialIdeal.variety`
        on a :class:`~sage.rings.polynomial.pbori.BooleanPolynomialIdeal`.

    .. NOTE::

        The order in which the solutions are returned is implementation
        dependent ; it may vary accross machines, and may vary accross
        calls to the function.

    """
    if eqs == []:
        raise ValueError("No equations, no solutions")

    eqs = Sequence(eqs)
    R = eqs.ring()
    if R.base_ring() != GF(2):
        raise ValueError("FES only deals with equations over GF(2)")

    degree = max( [f.degree() for f in eqs] )
    if degree <= 1:
        raise ValueError("the FES library requires equations to be non-linear")
    n = R.ngens()


    # ------- initialize a data-structure to communicate the equations to the library
    cdef int ***coeffs = <int ***> sig_calloc(len(eqs), sizeof(int **))
    for e,f in enumerate(eqs):
        coeffs[e] = <int **> sig_calloc(degree+1, sizeof(int *))
        for d in range(degree+1):
            coeffs[e][d] = <int *> sig_calloc(binomial(n,d), sizeof(int))

        for m in f:  # we enumerate the monomials of f
            d = m.degree()
            temp_n = n
            int_idx = 0
            prev_var = -1
            for idx in sorted(m.iterindex()):  # iterate over all the variables in the monomial
                j = idx - prev_var - 1
                for k in range(j):
                    int_idx += binomial(temp_n-k-1, d-1)
                prev_var = idx
                d -= 1
                temp_n -= 1 + j
            coeffs[e][ m.degree() ][int_idx] = 1

    internal_state = InternalState()
    internal_state.verbose = verbose
    internal_state.sols = []
    if max_sols == Infinity:
        internal_state.max_sols = 0
    else:
        internal_state.max_sols = max_sols

    # ------- runs the library
    sig_on()
    exhaustive_search_wrapper(n, len(eqs), degree, coeffs, report_solution, <void *> internal_state, verbose)
    sig_off()

    # ------- frees memory occupied by the dense representation of the equations
    if coeffs != NULL:
        for e in range(len(eqs)):
            if coeffs[e] != NULL:
                for d in range(degree+1):
                    if coeffs[e][d] != NULL:
                        sig_free(coeffs[e][d])
                sig_free(coeffs[e])
        sig_free(coeffs)

    # ------ convert (packed) solutions to suitable format
    dict_sols = []
    for i in internal_state.sols:
        sol = dict([])
        for x in R.gens():
            sol[x] = GF(2)(i & 1)
            i = i >> 1
        dict_sols.append( sol )

    return dict_sols


def find_coordinate_change(As, max_tries=64):
    """
    Tries to find a linear change of coordinates such that certain
    coefficients of the quadratic forms As become zero. This in turn
    makes the exhaustive search faster

    EXAMPLES:

    Testing that the function works::

        sage: from sage.libs.fes import find_coordinate_change    # optional - FES
        sage: n = 40                                              # optional - FES
        sage: K = GF(2)                                           # optional - FES
        sage: R = BooleanPolynomialRing(n, 'x')                   # optional - FES
        sage: foo = [ MatrixSpace(GF(2), n, n).random_element() ] # optional - FES
        sage: f = [ M + M.T for M in foo ]                        # optional - FES
        sage: S = find_coordinate_change( f )                     # optional - FES
        sage: g = [ S.T*M*S for M in f ]                          # optional - FES
        sage: vector ([ M[0,1] for M in g[:32]]).is_zero()        # optional - FES
        True
    """
    n = As[0].nrows()
    m = min(32, len(As))
    system = matrix(GF(2), m, n-2)
    S = identity_matrix(GF(2), n)
    Bs = As

    for foo in range(max_tries):
        for i in range(m):
            system[i] = Bs[i][0,2:]
            RHS = vector(GF(2), m, [ Bs[i][0,1] for i in range(m) ])
        try:
            solution = system.solve_right( RHS )
            S_prime = identity_matrix(GF(2), n)
            S_prime[1,2:] = solution
            return S*S_prime.T
        except ValueError:
            S.randomize()
            while not S.is_invertible():
                S.randomize()
            Bs = [ S.T*M*S for M in As ]
            print("trying again...")
    raise ValueError("Could not find suitable coordinate change")




def prepare_polynomials(f):
    """
    Finds a linear combination of the equations that is faster to solve by FES

    INPUT:

    - ``f`` -- list of boolean equations

    EXAMPLES:

    We check that the function does what it is supposed to do::

        sage: from sage.libs.fes import prepare_polynomials                    # optional - FES
        sage: n = 35                                                           # optional - FES
        sage: K = GF(2)                                                        # optional - FES
        sage: R = BooleanPolynomialRing(n, 'x')                                # optional - FES
        sage: f = [ sum( [K.random_element() * R.gen(i) * R.gen(j) for i in range(n) for j in range(i)] ) \
               + sum( [K.random_element() * R.gen(i)  for i in range(n) ] ) \
               + K.random_element() for l in range(n) ]                        # optional - FES
        sage: g = prepare_polynomials(f)                                       # optional - FES
        sage: [x.lm() for x in g]                                          # optional - FES, random
        0
    """
    if f == []:
        return []
    excess = len(f) - 32

    if excess <= 0:
        return f

    # now, try to cancel the `excess` first monomials
    s = Sequence(f)

    # switch to degrevlex if not already there
    R = s.ring()
    if R.term_order() != 'degrevlex':
        R2 = BooleanPolynomialRing( R.ngens(), R.variable_names(), order='degrevlex')
        s = Sequence( [R2(g) for g in s] )

    monomials_in_s = list( s.monomials() )
    monomials_in_s.sort(reverse=True)

    m = matrix(R.base_ring(), [ [ g.monomial_coefficient(m) for m in monomials_in_s[:excess] ] for g in s ])
    # now find the linear combinations of the equations that kills the first `excess` monomials in all but `excess` equations
    # todo, this is very likely suboptimal, but m.echelonize() does not return the transformation...
    P, L, U = m.LU()
    S = (P*L).I
    result = Sequence( S * vector(s) )
    result.reverse()
    return result
