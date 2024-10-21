# sage_setup: distribution = sagemath-flint
# sage.doctest: needs sage.geometry.polyhedron sage.libs.linbox sage.modules sage.rings.number_field
r"""
Enumeration of totally real fields: relative extensions

This module contains functions to enumerate primitive extensions `L / K`, where
`K` is a given totally real number field, with given degree and small root
discriminant. This is a relative analogue of the problem described in
:mod:`sage.rings.number_field.totallyreal`, and we use a similar approach
based on a relative version of Hunter's theorem.

In this first simple example, we compute the totally real quadratic
fields of `F = \QQ(\sqrt{2})` of discriminant `\le 2000`.

::

    sage: ZZx.<x> = ZZ[]
    sage: F.<t> = NumberField(x^2 - 2)
    sage: enumerate_totallyreal_fields_rel(F, 2, 2000)
    [[1600, x^4 - 6*x^2 + 4, xF^2 + xF - 1]]

There is indeed only one such extension, given by `F(\sqrt{5})`.

Next, we list all totally real quadratic extensions of `\QQ(\sqrt 5)`
with root discriminant `\le 10`.

::

    sage: F.<t> = NumberField(x^2 - 5)
    sage: ls = enumerate_totallyreal_fields_rel(F, 2, 10^4)
    sage: ls # random (the second factor is platform-dependent)
    [[725, x^4 - x^3 - 3*x^2 + x + 1, xF^2 + (-1/2*t - 7/2)*xF + 1],
     [1125, x^4 - x^3 - 4*x^2 + 4*x + 1, xF^2 + (-1/2*t - 7/2)*xF + 1/2*t + 3/2],
     [1600, x^4 - 6*x^2 + 4, xF^2 - 2],
     [2000, x^4 - 5*x^2 + 5, xF^2 - 1/2*t - 5/2],
     [2225, x^4 - x^3 - 5*x^2 + 2*x + 4, xF^2 + (-1/2*t + 1/2)*xF - 3/2*t - 7/2],
     [2525, x^4 - 2*x^3 - 4*x^2 + 5*x + 5, xF^2 + (-1/2*t - 1/2)*xF - 1/2*t - 5/2],
     [3600, x^4 - 2*x^3 - 7*x^2 + 8*x + 1, xF^2 - 3],
     [4225, x^4 - 9*x^2 + 4, xF^2 + (-1/2*t - 1/2)*xF - 3/2*t - 9/2],
     [4400, x^4 - 7*x^2 + 11, xF^2 - 1/2*t - 7/2],
     [4525, x^4 - x^3 - 7*x^2 + 3*x + 9, xF^2 + (-1/2*t - 1/2)*xF - 3],
     [5125, x^4 - 2*x^3 - 6*x^2 + 7*x + 11, xF^2 + (-1/2*t - 1/2)*xF - t - 4],
     [5225, x^4 - x^3 - 8*x^2 + x + 11, xF^2 + (-1/2*t - 1/2)*xF - 1/2*t - 7/2],
     [5725, x^4 - x^3 - 8*x^2 + 6*x + 11, xF^2 + (-1/2*t + 1/2)*xF - 1/2*t - 7/2],
     [6125, x^4 - x^3 - 9*x^2 + 9*x + 11, xF^2 + (-1/2*t + 1/2)*xF - t - 4],
     [7225, x^4 - 11*x^2 + 9, xF^2 + (-1)*xF - 4],
     [7600, x^4 - 9*x^2 + 19, xF^2 - 1/2*t - 9/2],
     [7625, x^4 - x^3 - 9*x^2 + 4*x + 16, xF^2 + (-1/2*t - 1/2)*xF - 4],
     [8000, x^4 - 10*x^2 + 20, xF^2 - t - 5],
     [8525, x^4 - 2*x^3 - 8*x^2 + 9*x + 19, xF^2 + (-1)*xF - 1/2*t - 9/2],
     [8725, x^4 - x^3 - 10*x^2 + 2*x + 19, xF^2 + (-1/2*t - 1/2)*xF - 1/2*t - 9/2],
     [9225, x^4 - x^3 - 10*x^2 + 7*x + 19, xF^2 + (-1/2*t + 1/2)*xF - 1/2*t - 9/2]]
    sage: [ f[0] for f in ls ]
    [725, 1125, 1600, 2000, 2225, 2525, 3600, 4225, 4400, 4525, 5125, 5225, 5725, 6125, 7225, 7600, 7625, 8000, 8525, 8725, 9225]

    sage: [NumberField(ZZx(x[1]), 't').is_galois() for x in ls]                         # needs sage.groups
    [False, True, True, True, False, False, True, True, False, False, False, False, False, True, True, False, False, True, False, False, False]

Eight out of 21 such fields are Galois (with Galois group `C_4`
or `C_2 \times C_2`); the others have Galois closure of degree 8
(with Galois group `D_8`).

Finally, we compute the cubic extensions of `\QQ(\zeta_7)^+` with
discriminant `\le 17 \times 10^9`.

::

    sage: F.<t> = NumberField(ZZx([1,-4,3,1]))
    sage: F.disc()
    49
    sage: enumerate_totallyreal_fields_rel(F, 3, 17*10^9)  # not tested, too long time (258s on sage.math, 2013)
    [[16240385609L, x^9 - x^8 - 9*x^7 + 4*x^6 + 26*x^5 - 2*x^4 - 25*x^3 - x^2 + 7*x + 1, xF^3 + (-t^2 - 4*t + 1)*xF^2 + (t^2 + 3*t - 5)*xF + 3*t^2 + 11*t - 5]]    # 32-bit
    [[16240385609, x^9 - x^8 - 9*x^7 + 4*x^6 + 26*x^5 - 2*x^4 - 25*x^3 - x^2 + 7*x + 1, xF^3 + (-t^2 - 4*t + 1)*xF^2 + (t^2 + 3*t - 5)*xF + 3*t^2 + 11*t - 5]]     # 64-bit

TESTS:

Check that :issue:`27646` is fixed::

    sage: L = enumerate_totallyreal_fields_all(6,435000) # long time

AUTHORS:

- John Voight (2007-11-03): initial version
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein and John Voight
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import binomial
from sage.arith.misc import GCD as gcd
from sage.arith.misc import divisors
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.number_field.totallyreal_data import ZZx, lagrange_degree_3, int_has_small_square_divisor, hermite_constant
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.totallyreal import weed_fields, odlyzko_bound_totallyreal, enumerate_totallyreal_fields_prim
from sage.libs.pari.all import pari
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

import math
import sys


def integral_elements_in_box(K, C):
    r"""
    Return all integral elements of the totally real field `K` whose
    embeddings lie *numerically* within the bounds specified by the
    list ``C``.  The output is architecture dependent, and one may want
    to expand the bounds that define ``C`` by some epsilon.

    INPUT:

    - ``K`` -- a totally real number field
    - ``C`` -- list ``[[lower, upper], ...]`` of lower and upper bounds,
      for each embedding

    EXAMPLES::

        sage: x = polygen(QQ)
        sage: K.<alpha> = NumberField(x^2 - 2)
        sage: eps = 10e-6
        sage: C = [[0-eps, 5+eps], [0-eps, 10+eps]]
        sage: ls = sage.rings.number_field.totallyreal_rel.integral_elements_in_box(K, C)
        sage: sorted(a.trace() for a in ls)
        [0, 2, 4, 4, 4, 6, 6, 6, 6, 8, 8, 8, 10, 10, 10, 10, 12, 12, 14]
        sage: len(ls)
        19

        sage: v = sage.rings.number_field.totallyreal_rel.integral_elements_in_box(K, C)
        sage: sorted(v)
        [0, -alpha + 2, 1, -alpha + 3, 2, 3, alpha + 2, 4, alpha + 3, 5, alpha + 4,
         2*alpha + 3, alpha + 5, 2*alpha + 4, alpha + 6, 2*alpha + 5, 2*alpha + 6,
         3*alpha + 5, 2*alpha + 7]

    A cubic field::

        sage: x = polygen(QQ)
        sage: K.<a> = NumberField(x^3 - 16*x +16)
        sage: eps = 10e-6
        sage: C = [[0-eps,5+eps]]*3
        sage: v = sage.rings.number_field.totallyreal_rel.integral_elements_in_box(K, C)

    Note that the output is platform dependent (sometimes a 5 is listed
    below, and sometimes it isn't)::

        sage: sorted(v)
        [-1/2*a + 2, 1/4*a^2 + 1/2*a, 0, 1, 2, 3, 4,...-1/4*a^2 - 1/2*a + 5,
         1/2*a + 3, -1/4*a^2 + 5]
    """
    d = K.degree()
    Foo = K.real_embeddings()
    B = K.reduced_basis()

    import numpy
    import numpy.linalg
    L = numpy.array([[v(b) for b in B] for v in Foo])
    Linv = numpy.linalg.inv(L)
    Vi = [[C[0][0]], [C[0][1]]]
    for i in range(1, d):
        Vi = sum([[v + [C[i][0]], v + [C[i][1]]] for v in Vi], [])
    V = numpy.matrix(Linv) * numpy.matrix(Vi).transpose()
    j = 0
    while j < 2**d:
        for i in range(d):
            if V[i, j] < V[i, j+1]:
                V[i, j] = math.floor(V[i, j])
                V[i, j+1] = math.ceil(V[i, j+1])
            else:
                V[i, j] = math.ceil(V[i, j])
                V[i, j+1] = math.floor(V[i, j+1])
        j += 2
    W0 = (Linv*numpy.array([Vi[0]] * d)).transpose()
    W = (Linv*numpy.array([Vi[2**i] for i in range(d)])).transpose()
    for j in range(d):
        for i in range(d):
            if W[i,j] < W0[i,j]:
                W[i,j] = math.floor(W[i,j])
                W0[i,j] = math.ceil(W0[i,j])
            else:
                W[i,j] = math.ceil(W[i,j])
                W0[i,j] = math.floor(W0[i,j])
    M = [[int(V[i,j]) for i in range(V.shape[0])] for j in range(V.shape[1])]
    M += [[int(W0[i,j]) for j in range(W0.shape[0])] for i in range(W0.shape[0])]
    M += [[int(W[i,j]) for j in range(W.shape[1])] for i in range(W.shape[0])]

    from sage.matrix.constructor import matrix
    M = (matrix(IntegerRing(),len(M),len(M[0]), M).transpose()).columns()

    i = 0
    while i < len(M):
        j = i+1
        while j < len(M):
            if M[i] == M[j]:
                M.pop(j)
            else:
                j += 1
        i += 1

    from sage.geometry.lattice_polytope import LatticePolytope
    P = LatticePolytope(M)

    try:
        pts = P.points()
    except ValueError:
        return []

    S = []
    for p in pts:
        theta = sum(a * b for a, b in zip(p.list(), B))
        if all((C[i][0] <= Foo[i](theta) <= C[i][1])
               for i in range(d)):
            S.append(theta)

    return S


#********************************************************************************
# Main class
#********************************************************************************

eps_global = 10**(-6)


class tr_data_rel:
    r"""
    This class encodes the data used in the enumeration of totally real
    fields for relative extensions.

    We do not give a complete description here.  For more information,
    see the attached functions; all of these are used internally by the
    functions in totallyreal_rel.py, so see that file for examples and
    further documentation.
    """

    def __init__(self, F, m, B, a=None):
        r"""
        Initialization routine (constructor).

        INPUT:

        - ``F`` -- number field; the base field
        - ``m`` -- integer; the relative degree
        - ``B`` -- integer; the discriminant bound
        - ``a`` -- list (default: ``[]``); the coefficient list to begin with,
          corresponding to ``a[len(a)]*x^n + ... + a[0]x^(n-len(a))``

        OUTPUT:

        the data initialized to begin enumeration of totally real fields
        with base field `F`, degree `n`, discriminant bounded by `B`, and starting
        with coefficients `a`.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: F.<t> = NumberField(x^2 - 2)
            sage: T = sage.rings.number_field.totallyreal_rel.tr_data_rel(F, 2, 2000)
        """
        if a is None:  # don't make the stupid noob mistake of putting a=[]
            a = []     # in the function signature above.

        # Initialize constants.
        self.m = m
        d = F.degree()
        self.d = d
        self.n = n = m*d
        self.B = B
        self.gamma = hermite_constant(self.n-self.d)

        self.F = F
        self.Z_F = F.maximal_order()
        self.Foo = F.real_embeddings()
        self.dF = abs(F.disc())
        self.Fx = PolynomialRing(F, 'xF')

        self.beta = [[]]*m
        self.gnk = [[]]*m

        self.trace_elts = []

        Z_Fbasis = self.Z_F.basis()

        # Initialize variables.
        if not a:
            # No starting input, all polynomials will be found; initialize to zero.
            self.a = [0]*m + [1]
            self.amaxvals = [[]]*m
            anm1s = [[i] for i in range(0,m//2+1)]
            for i in range(1,self.d):
                for j in range(len(anm1s)):
                    anm1s[j] = [anm1s[j] + [i] for i in range(m)]
                anm1s = sum(anm1s, [])
            anm1s = [sum([Z_Fbasis[i] * aa[i] for i in range(self.d)]) for aa in anm1s]
            # Minimize trace in class.
            import numpy
            for i in range(len(anm1s)):
                Q = [[v(m*x) for v in self.Foo] + [0] for x in Z_Fbasis] + [[v(anm1s[i]) for v in self.Foo] + [10**6]]
                pari_string = '[' + ';'.join(','.join("%s" % ii for ii in row) for row in zip(*Q)) + ']'
                adj = pari(pari_string).qflll()[self.d]
                anm1s[i] += sum([m*Z_Fbasis[ii]*int(adj[ii])//int(adj[self.d]) for ii in range(self.d)])

            self.amaxvals[m-1] = anm1s
            self.a[m-1] = self.amaxvals[m-1].pop()
            self.k = m-2

            bl = math.ceil(1.7719*self.n)
            br = max([1./m*(am1**2).trace() +
                      self.gamma*(1./(m**d)*self.B/self.dF)**(1./(self.n-d))
                      for am1 in anm1s])
            br = math.floor(br)
            T2s = self.F._positive_integral_elements_with_trace([bl, br])
            self.trace_elts.append([bl, br, T2s])

        elif len(a) <= m + 1:
            # First few coefficients have been specified.
            # The value of k is the largest index of the coefficients of a which is
            # currently unknown; e.g., if k == -1, then we can iterate
            # over polynomials, and if k == n-1, then we have finished iterating.
            if a[len(a)-1] != 1:
                raise ValueError("a[len(a)-1](=%s) must be 1 so polynomial is monic" % a[len(a)-1])

            raise NotImplementedError("These have not been checked.")

            k = m-len(a)
            self.k = k
            a = [0]*(k+1) + a
            self.amaxvals = [[]]*m
            for i in range(n+1):
                self.a[i] = a[i]

            # Bounds come from an application of Lagrange multipliers in degrees 2,3.
            self.b_lower = [-1./m*(v(self.a[m-1]) +
                                   (m-1.)*math.sqrt(v(self.a[m-1])**2 - 2.*(1+1./(m-1))*v(self.a[m-2]))) for v in self.Foo]
            self.b_upper = [-1./m*(v(self.a[m-1]) -
                                   (m-1.)*math.sqrt(v(self.a[m-1])**2 - 2.*(1+1./(m-1))*v(self.a[m-2]))) for v in self.Foo]
            if k < m - 2:
                bminmax = [lagrange_degree_3(n,v(self.a[m-1]),v(self.a[m-2]),v(self.a[m-3])) for v in self.Foo]
                self.b_lower = bminmax[0]
                self.b_upper = bminmax[1]

            # Annoying, but must reverse coefficients for numpy.
            gnk = [binomial(j,k+2)*a[j] for j in range(k+2,n+1)]
            self.beta[k+1] = [[self.b_lower] + numpy.roots([v(gnk[i]) for i in range(len(gnk))].reverse()).tolist().sort() + [self.b_upper] for v in self.Foo]

            # Now to really initialize gnk.
            self.gnk[k+1] = [[0] + [binomial(j,k+1)*v(a[j])
                                    for j in range(k+2,m+1)] for v in self.Foo]
        else:
            # Bad input!
            raise ValueError("a has length %s > m+1" % len(a))

    def incr(self, f_out, verbose=False, haltk=0):
        r"""
        'Increment' the totally real data to the next
        value which satisfies the bounds essentially given by Rolle's
        theorem, and return the next polynomial in the sequence
        ``f_out``.

        The default or usual case just increments the constant
        coefficient; then inductively, if this is outside of the
        bounds we increment the next higher coefficient, and so on.

        If there are no more coefficients to be had, returns the zero
        polynomial.

        INPUT:

        - ``f_out`` -- integer sequence; to be written with the coefficients of
          the next polynomial
        - ``verbose`` -- boolean or nonnegative integer (default: ``False``);
          print verbosely computational details. It prints extra information if
          ``verbose`` is set to ``2`` or more.
        - ``haltk`` -- integer; the level at which to halt the inductive
          coefficient bounds

        OUTPUT: the successor polynomial as a coefficient list
        """
        import numpy

        m = self.m
        k = self.k
        d = self.d

        # If k == -1, we have a full polynomial, so we add 1 to the constant coefficient.
        if k == -1:
            if len(self.amaxvals[0]) > 0 and self.amaxvals[0]:
                self.a[0] = self.amaxvals[0].pop()
                for i in range(0,m):
                    f_out[i] = self.a[i]
                return
            else:
                if verbose:
                    print("  finished")

                # Already reached maximum, so "carry the 1" to find the next value of k.
                k += 1
                while k < m and len(self.amaxvals[k]) == 0:
                    k += 1
                if k < m:
                    self.a[k] = self.amaxvals[k].pop()
                    k -= 1

        # If we are working through an initialization routine, treat that.
        elif haltk and k == haltk-1:
            if len(self.maxvals[k]) == 0:
                k += 1
                while k <= m-1 and len(self.amaxvals[k]) == 0:
                    k += 1
                if k < m:
                    self.a[k] = self.amaxvals[k].pop()
                    k -= 1

        # If in the previous step we finished all possible values of
        # the lastmost coefficient, so we must compute bounds on the next coefficient.
        # Recall k == n-1 implies iteration is complete.
        while k < m-1:
            # maxoutflag flags a required abort along the way
            maxoutflag = False

            # Recall k == -1 means all coefficients are good to go.
            while k >= 0 and (not haltk or k >= haltk):
                if verbose:
                    print(k, ":", end="")
                    for i in range(self.m + 1):
                        print(self.a[i], end="")
                    print("")

                if k == m - 2:
                    # We only know the value of a[n-1], the trace.
                    bl = max(math.ceil(1.7719*self.n), ((self.a[m-1]**2).trace()*1./m))
                    br = 1./m*(self.a[m-1]**2).trace() + \
                         self.gamma*(1./(m**d)*self.B/self.dF)**(1./(self.n-d))
                    br = math.floor(br)

                    # Check for trivially empty.
                    if bl > br:
                        if verbose:
                            print(" ", br, ">", bl)
                        maxoutflag = 1
                        break

                    if verbose >= 2:
                        print("  bl, br:", bl, br)

                    # Enumerate all elements of Z_F with T_2 <= br
                    T2s = []
                    trace_elts_found = False
                    for tre in self.trace_elts:
                        if tre[0] <= bl and tre[1] >= br:
                            trace_elts_found = True
                            if verbose >= 2:
                                print("  found copy!")
                            T2s.extend(theta for theta in tre[2]
                                       if bl <= theta.trace() <= br)
                            break
                    if not trace_elts_found:
                        T2s = self.F._positive_integral_elements_with_trace([bl,br])
                        self.trace_elts.append([bl,br,T2s])

                    # Now ensure that T2 satisfies the correct parity condition
                    am2s = []
                    for t2 in T2s:
                        am2 = (self.a[m-1]**2-t2)/2
                        if am2.is_integral():
                            ispositive = True
                            for v in self.Foo:
                                ispositive = ispositive and v((m-1)*self.a[m-1]**2-2*m*am2) > 0
                            if ispositive:
                                am2s.append(am2)

                    if verbose >= 2:
                        print("  am2s:", am2s)

                    # If none survive, break!
                    if len(am2s) == 0:
                        if verbose:
                            print("  finished")
                        maxoutflag = 1
                        break

                    self.amaxvals[m-2] = am2s
                    self.a[m-2] = self.amaxvals[m-2].pop()

                    # Initialize the second derivative.
                    self.b_lower = [-1./m*(v(self.a[m-1]) +
                                           (m-1.)*math.sqrt(v(self.a[m-1])**2 - 2.*(1+1./(m-1))*v(self.a[m-2]))) for v in self.Foo]
                    self.b_upper = [-1./m*(v(self.a[m-1]) -
                                           (m-1.)*math.sqrt(v(self.a[m-1])**2 - 2.*(1+1./(m-1))*v(self.a[m-2]))) for v in self.Foo]
                    self.beta[k] = [[self.b_lower[i], -self.Foo[i](self.a[m-1])/m, self.b_upper[i]] for i in range(d)]
                    self.gnk[k] = [0, (m-1)*self.a[m-1], m*(m-1)/2]

                    if verbose >= 2:
                        print("  betak:", self.beta[k])
                else:
                    # Compute the roots of the derivative.
                    self.gnk[k+1][0] = self.a[k+1]
                    gnk = self.gnk[k+1]
                    self.beta[k] = [numpy.roots([v(gnk[len(gnk)-1-i]) for i in range(len(gnk))]).tolist() for v in self.Foo]

                    try:
                        for i in range(d):
                            self.beta[k][i].sort()
                    except TypeError:
                        if verbose:
                            print("  betak:", self.beta[k])
                        maxoutflag = True
                        break

                    # Check for double roots
                    for i in range(len(self.beta[k][0])-1):
                        if abs(self.beta[k][0][i] - self.beta[k][0][i+1]) < 2*eps_global:
                            # This happens reasonably infrequently, so calling
                            # the Python routine should be sufficiently fast...
                            f = self.Fx(self.gnk[k+1])
                            df = self.Fx(self.gnk[k+2])
                            if gcd(f,df) != 1:
                                if verbose:
                                    print("  gnk has multiple factor!")
                                maxoutflag = True
                                break
                    if maxoutflag:
                        break

                    if k == m-3:
                        self.b_lower = [-1./m*(v(self.a[m-1]) +
                                               (m-1.)*math.sqrt(v(self.a[m-1])**2 - 2.*(1+1./(m-1))*v(self.a[m-2]))) for v in self.Foo]
                        self.b_upper = [-1./m*(v(self.a[m-1]) -
                                               (m-1.)*math.sqrt(v(self.a[m-1])**2 - 2.*(1+1./(m-1))*v(self.a[m-2]))) for v in self.Foo]
                    elif k == m-4:
                        # New bounds from Lagrange multiplier in degree 3.
                        bminmax = [lagrange_degree_3(m,v(self.a[m-1]),v(self.a[m-2]),v(self.a[m-3])) for v in self.Foo]
                        self.b_lower = [bminmax[i][0] for i in range(len(bminmax))]
                        self.b_upper = [bminmax[i][1] for i in range(len(bminmax))]

                    self.beta[k] = [[self.b_lower[i]] + self.beta[k][i] + [self.b_upper[i]] for i in range(len(self.beta[k]))]

                    if verbose >= 2:
                        print("  betak:", self.beta[k])

                    # Compute next g_(m-(k+1)), k times the formal integral of g_(m-k).
                    self.gnk[k] = [self.F.primitive_element()*0] + [self.gnk[k+1][i-1]*(k+1)/i for i in range(1,m-k+1)]
                    gnk = self.gnk[k]
                    gnks = [[v(gnk[len(gnk)-1-i])
                             for i in range(len(gnk))]
                            for v in self.Foo]
                    gnkm1 = self.gnk[k+1]
                    gnkm1s = [[v(gnkm1[len(gnkm1)-1-i])
                               for i in range(len(gnkm1))]
                              for v in self.Foo]
                    mk = m - (k + 1)

                    if verbose >= 2:
                        print("  gnk:", self.gnk[k])
                        print("  gnks:", gnks)

                    # Compute upper and lower bounds which guarantee one retains
                    # a polynomial with all real roots.
                    betak = self.beta[k]
                    akmin = [-numpy.polyval(gnks[j], betak[j][mk+1]) -
                             abs(numpy.polyval(gnkm1s[j], betak[j][mk+1]))*eps_global for j in range(self.d)]
                    for i in range(1,(mk+1)//2+1):
                        # Use the fact that f(z) <= f(x)+|f'(x)|eps if |x-z| < eps
                        # for sufficiently small eps, f(z) = 0, and f''(z) < 0.
                        akmin = [max(akmin[j],
                                     -numpy.polyval(gnks[j], betak[j][mk+1-2*i]) -
                                     abs(numpy.polyval(gnkm1s[j], betak[j][mk+1-2*i])*eps_global)) for j in range(self.d)]

                    akmax = [-numpy.polyval(gnks[j], betak[j][mk]) +
                             abs(numpy.polyval(gnkm1s[j], betak[j][mk]))*eps_global for j in range(self.d)]
                    for i in range(1, mk//2+1):
                        akmax = [min(akmax[j],
                                     -numpy.polyval(gnks[j], betak[j][mk-2*i]) +
                                     abs(numpy.polyval(gnkm1s[j], betak[j][mk-2*i])*eps_global)) for j in range(self.d)]

                    if verbose >= 2:
                        print("  akmin:", akmin)
                        print("  akmax:", akmax)

                    for i in range(self.d):
                        if akmin[i] > akmax[i]:
                            if verbose:
                                print(" ", akmin[i], ">", akmax[i])
                            maxoutflag = 1
                            break
                    if maxoutflag:
                        break

                    self.amaxvals[k] = integral_elements_in_box(self.F, [[akmin[i],akmax[i]] for i in range(d)])
                    if k == 0:
                        a0s = [0, -sum([self.a[i] for i in range(1,m+1)]),
                               -sum([self.a[i]*(-1)**i for i in range(1,m+1)]),
                               -sum([self.a[i]*2**i for i in range(1,m+1)]),
                               -sum([self.a[i]*(-2)**i for i in range(1,m+1)])]
                        for a0 in a0s:
                            try:
                                self.amaxvals[0].remove(a0)
                            except Exception:
                                pass

                    if verbose:
                        print("  amaxvals[k]:", self.amaxvals[k])
                    if len(self.amaxvals[k]) == 0:
                        if verbose:
                            print("  finished")
                        maxoutflag = True
                        break
                    self.a[k] = self.amaxvals[k].pop()

                self.k -= 1
                k -= 1

            if not maxoutflag:
                self.k = k
                for i in range(m):
                    f_out[i] = self.a[i]
                return
            else:
                k += 1
                while k < m and len(self.amaxvals[k]) == 0:
                    k += 1
                if k < m:
                    self.a[k] = self.amaxvals[k].pop()
                    k -= 1

        # k == n-1, so iteration is complete; return the zero polynomial (of degree n+1).
        self.k = k
        f_out[m] = 0
        return


# ****************************************************************************
#                                 Main routine
# ****************************************************************************

def enumerate_totallyreal_fields_rel(F, m, B, a=[], verbose=0,
                                     return_seqs=False,
                                     return_pari_objects=True):
    r"""
    This function enumerates (primitive) totally real field extensions of
    degree `m>1` of the totally real field F with discriminant `d \leq B`;
    optionally one can specify the first few coefficients, where the sequence ``a``
    corresponds to a polynomial by

    ::

        a[d]*x^n + ... + a[0]*x^(n-d)

    if ``length(a) = d+1``, so in particular always ``a[d] = 1``.

    .. NOTE::

        This is guaranteed to give all primitive such fields, and
        seems in practice to give many imprimitive ones.

    INPUT:

    - ``F`` -- number field; the base field
    - ``m`` -- integer; the degree
    - ``B`` -- integer; the discriminant bound
    - ``a`` -- list (default: ``[]``); the coefficient list to begin with
    - ``verbose`` -- boolean or nonnegative integer or string (default: 0);
      give a verbose description of the computations being performed. If
      ``verbose`` is set to ``2`` or more then it outputs some extra
      information. If ``verbose`` is a string then it outputs to a file
      specified by ``verbose``.
    - ``return_seqs`` -- boolean (default: ``False``); if ``True``, then return
      the polynomials as sequences (for easier exporting to a file). This
      also returns a list of four numbers, as explained in the OUTPUT
      section below.
    - ``return_pari_objects`` -- boolean (default: ``True``); if both
      ``return_seqs`` and ``return_pari_objects`` are ``False`` then it returns
      the elements as Sage objects; otherwise it returns PARI objects.

    OUTPUT:

    - the list of fields with entries ``[d,fabs,f]``, where ``d`` is the
      discriminant, ``fabs`` is an absolute defining polynomial, and ``f``
      is a defining polynomial relative to `F`, sorted by discriminant.

    - if ``return_seqs`` is ``True``, then the first field of the list is
      a list containing the count of four items as explained below

      - the first entry gives the number of polynomials tested
      - the second entry gives the number of polynomials with its
        discriminant having a large enough square divisor
      - the third entry is the number of irreducible polynomials
      - the fourth entry is the number of irreducible polynomials with
        discriminant at most `B`

    EXAMPLES::

        sage: ZZx.<x> = ZZ[]
        sage: F.<t> = NumberField(x^2 - 2)
        sage: enumerate_totallyreal_fields_rel(F, 1, 2000)
        [[1, [-2, 0, 1], xF - 1]]
        sage: enumerate_totallyreal_fields_rel(F, 2, 2000)
        [[1600, x^4 - 6*x^2 + 4, xF^2 + xF - 1]]
        sage: enumerate_totallyreal_fields_rel(F, 2, 2000, return_seqs=True)
        [[9, 6, 5, 0], [[1600, [4, 0, -6, 0, 1], [-1, 1, 1]]]]

    TESTS:

    Each of the outputs must be elements of Sage if ``return_pari_objects``
    is set to ``False``::

        sage: type(enumerate_totallyreal_fields_rel(F, 2, 2000)[0][1])
        <class 'cypari2.gen.Gen'>
        sage: enumerate_totallyreal_fields_rel(F, 2, 2000, return_pari_objects=False)[0][0].parent()
        Integer Ring
        sage: enumerate_totallyreal_fields_rel(F, 2, 2000, return_pari_objects=False)[0][1].parent()
        Univariate Polynomial Ring in x over Rational Field
        sage: enumerate_totallyreal_fields_rel(F, 2, 2000, return_pari_objects=False)[0][2].parent()
        Univariate Polynomial Ring in xF over Number Field in t with defining polynomial x^2 - 2
        sage: enumerate_totallyreal_fields_rel(F, 2, 2000, return_seqs=True)[1][0][1][0].parent()
        Rational Field

    AUTHORS:

    - John Voight (2007-11-01)
    """

    if not isinstance(m, Integer):
        try:
            m = Integer(m)
        except TypeError:
            raise TypeError("cannot coerce m (= %s) to an integer" % m)
    if (m < 1):
        raise ValueError("m must be at least 1.")

    n = F.degree()*m

    # Initialize
    S = {}        # dictionary of the form {(d, fabs): f, ...}
    dB_odlyzko = odlyzko_bound_totallyreal(n)
    dB = math.ceil(40000*dB_odlyzko**n)
    counts = [0,0,0,0]

    # Trivial case
    if m == 1:
        g = pari(F.defining_polynomial()).polrecip().Vec()
        if return_seqs:
            return [[0,0,0,0], [1, [-1, 1], g]]
        elif return_pari_objects:
            return [[1, g, pari('xF-1')]]
        else:
            Px = PolynomialRing(QQ, 'xF')
            return [[ZZ(1), [QQ(_) for _ in g], Px.gen()-1]]

    if verbose:
        saveout = sys.stdout
        if isinstance(verbose, str):
            fsock = open(verbose, 'w')
            sys.stdout = fsock
        # Else, print to screen
    f_out = [0]*m + [1]
    T = tr_data_rel(F, m, B, a)
    if verbose == 2:
        T.incr(f_out,verbose)
    else:
        T.incr(f_out)

    Fx = PolynomialRing(F, 'xF')

    nfF = pari(str(F.defining_polynomial()).replace('x',
                                                    str(F.primitive_element())))
    parit = pari(str(F.primitive_element()))

    while f_out[m] != 0:
        counts[0] += 1
        if verbose:
            print("==>", f_out, end="")

        f_str = ''
        for i in range(len(f_out)):
            f_str += '(' + str(f_out[i]) + ')*x^' + str(i)
            if i < len(f_out)-1:
                f_str += '+'
        nf = pari(f_str)
        if nf.poldegree('t') == 0:
            nf = nf.subst('x', 'x-t')
        nf = nf.polresultant(nfF, parit)  # either monic or -monic
        if nf[n] == -1:
            nf *= -1
        d = nf.poldisc()
        #counts[0] += 1
        if d > 0 and nf.polsturm() == n:
            da = int_has_small_square_divisor(Integer(d))
            if d > dB or d <= B*da:
                counts[1] += 1
                if nf.polisirreducible():
                    counts[2] += 1
                    [zk,d] = nf.nfbasis_d()

                    if d <= B:
                        if verbose:
                            print("has discriminant", d, end="")

                        # Find a minimal lattice element
                        counts[3] += 1
                        ng = pari([nf,zk]).polredabs()

                        # Check if K is contained in the list.
                        if (d, ng) in S:
                            if verbose:
                                print("but is not new")
                        else:
                            if verbose:
                                print("and is new!")
                            S[(d, ng)] = Fx(f_out)
                    else:
                        if verbose:
                            print("has discriminant", abs(d), "> B")
                else:
                    if verbose:
                        print("is not absolutely irreducible")
            else:
                if verbose:
                    print("has discriminant", abs(d), "with no large enough square divisor")
        else:
            if verbose:
                if d == 0:
                    print("is not squarefree")
                else:
                    print("is not totally real")
        if verbose == 2:
            T.incr(f_out,verbose=verbose)
        else:
            T.incr(f_out)

    # In the application of Smyth's theorem above, we exclude finitely
    # many possibilities which we must now throw back in.
    if m == 2:
        if Fx([-1,1,1]).is_irreducible():
            K = F.extension(Fx([-1,1,1]), 'tK')
            Kabs = K.absolute_field('tKabs')
            Kabs_pari = pari(Kabs.defining_polynomial())
            d = K.absolute_discriminant()
            if abs(d) <= B:
                ng = Kabs_pari.polredabs()
                S[(d, ng)] = Fx([-1,1,1])
        elif F.degree() == 2:
            for ff in [[1,-7,13,-7,1],[1,-8,14,-7,1]]:
                f = Fx(ff).factor()[0][0]
                K = F.extension(f, 'tK')
                Kabs = K.absolute_field('tKabs')
                Kabs_pari = pari(Kabs.defining_polynomial())
                d = K.absolute_discriminant()
                if abs(d) <= B:
                    ng = Kabs_pari.polredabs()
                    S[(d, ng)] = f
    elif m == 3:
        if Fx([-1, 6, -5, 1]).is_irreducible():
            K = F.extension(Fx([-1, 6, -5, 1]), 'tK')
            Kabs = K.absolute_field('tKabs')
            Kabs_pari = pari(Kabs.defining_polynomial())
            d = K.absolute_discriminant()
            if abs(d) <= B:
                ng = Kabs_pari.polredabs()
                S[(d, ng)] = Fx([-1, 6, -5, 1])

    # Convert S to a sorted list of triples [d, fabs, f], taking care
    # not to use the comparison operators on PARI polynomials.
    S = [[s[0], s[1], t] for s, t in S.items()]
    S.sort(key=lambda x: (x[0], [QQ(cf) for cf in x[1].polrecip().Vec()]))

    # Now check for isomorphic fields
    weed_fields(S)

    # Output.
    if verbose:
        print("=" * 80)
        print("Polynomials tested: {}".format(counts[0]))
        print("Polynomials with discriminant with large enough square"
              " divisor: {}".format(counts[1]))
        print("Irreducible polynomials: {}".format(counts[2]))
        print("Polynomials with nfdisc <= B: {}".format(counts[3]))
        for i in range(len(S)):
            print(S[i])
        if isinstance(verbose, str):
            fsock.close()
        sys.stdout = saveout

    # Make sure to return elements that belong to Sage
    if return_seqs:
        return [[ZZ(x) for x in counts],
                [[s[0], [QQ(x) for x in s[1].polrecip().Vec()],
                  s[2].coefficients(sparse=False)]
                 for s in S]]
    elif return_pari_objects:
        return S
    else:
        Px = PolynomialRing(QQ, 'x')
        return [[s[0], Px([QQ(_) for _ in s[1].list()]), s[2]] for s in S]


def enumerate_totallyreal_fields_all(n, B, verbose=0, return_seqs=False,
                                     return_pari_objects=True):
    r"""
    Enumerate *all* totally real fields of degree ``n`` with discriminant
    at most ``B``, primitive or otherwise.

    INPUT:

    - ``n`` -- integer; the degree
    - ``B`` -- integer; the discriminant bound
    - ``verbose`` -- boolean or nonnegative integer or string (default: 0);
      give a verbose description of the computations being performed. If
      ``verbose`` is set to ``2`` or more, it outputs some extra information.
      If ``verbose`` is a string, it outputs to a file specified by ``verbose``.
    - ``return_seqs`` -- boolean (default: ``False``); if ``True``, then return
      the polynomials as sequences (for easier exporting to a file). This
      also returns a list of four numbers, as explained in the OUTPUT
      section below.
    - ``return_pari_objects`` -- boolean (default: ``True``); if both
      ``return_seqs`` and ``return_pari_objects`` are ``False`` then it
      returns the elements as Sage objects; otherwise it returns PARI objects.

    EXAMPLES::

        sage: enumerate_totallyreal_fields_all(4, 2000)
        [[725, x^4 - x^3 - 3*x^2 + x + 1],
        [1125, x^4 - x^3 - 4*x^2 + 4*x + 1],
        [1600, x^4 - 6*x^2 + 4],
        [1957, x^4 - 4*x^2 - x + 1],
        [2000, x^4 - 5*x^2 + 5]]
        sage: enumerate_totallyreal_fields_all(1, 10)
        [[1, x - 1]]

    TESTS:

    Each of the outputs must be elements of Sage if ``return_pari_objects``
    is set to ``False``::

        sage: enumerate_totallyreal_fields_all(2, 10)
        [[5, x^2 - x - 1], [8, x^2 - 2]]
        sage: type(enumerate_totallyreal_fields_all(2, 10)[0][1])
        <class 'cypari2.gen.Gen'>
        sage: enumerate_totallyreal_fields_all(2, 10, return_pari_objects=False)[0][1].parent()
        Univariate Polynomial Ring in x over Rational Field

    In practice most of these will be found by
    :func:`~sage.rings.number_field.totallyreal.enumerate_totallyreal_fields_prim`,
    which is guaranteed to return all primitive fields but often returns
    many non-primitive ones as well. For instance, only one of the five
    fields in the example above is primitive, but
    :func:`~sage.rings.number_field.totallyreal.enumerate_totallyreal_fields_prim`
    finds four out of the five (the exception being `x^4 - 6x^2 + 4`).

    The following was fixed in :issue:`13101`::

        sage: enumerate_totallyreal_fields_all(8, 10^6)  # long time (about 2 s)
        []
    """
    S = []
    counts = [0, 0, 0, 0]
    div_n = divisors(n)
    if len(div_n) > 4:
        raise ValueError("Only implemented for n = p*q with p,q prime")
    for d in div_n:
        if 1 < d < n:
            Sds = enumerate_totallyreal_fields_prim(d, int(math.floor((1.*B)**(1.*d/n))), verbose=verbose)
            for i in range(len(Sds)):
                if verbose:
                    print("=" * 80)
                    print("Taking F =", Sds[i][1])
                F = NumberField(ZZx(Sds[i][1]), 't')
                T = enumerate_totallyreal_fields_rel(F, n/d, B, verbose=verbose, return_seqs=return_seqs)
                if return_seqs:
                    for k in range(4):
                        counts[k] += T[0][k]
                    S += [[t[0],pari(t[1]).Polrev()] for t in T[1]]
                else:
                    S += [[t[0],t[1]] for t in T]
                for E in enumerate_totallyreal_fields_prim(n/d, int(math.floor((1.*B)**(1./d)/(1.*Sds[i][0])**(n*1./d**2)))):
                    for EF in F.composite_fields(NumberField(ZZx(E[1]), 'u')):
                        if EF.degree() == n and EF.disc() <= B:
                            S.append([EF.disc(), pari(EF.absolute_polynomial())])
    S += enumerate_totallyreal_fields_prim(n, B, verbose=verbose)
    S.sort(key=lambda x: (x[0], [QQ(cf) for cf in x[1].polrecip().Vec()]))
    weed_fields(S)

    # Output.
    if verbose:
        saveout = sys.stdout
        if isinstance(verbose, str):
            fsock = open(verbose, 'w')
            sys.stdout = fsock
        # Else, print to screen
        print("=" * 80)
        print("Polynomials tested: {}".format(counts[0]))
        print("Polynomials with discriminant with large enough square"
              " divisor: {}".format(counts[1]))
        print("Irreducible polynomials: {}".format(counts[2]))
        print("Polynomials with nfdisc <= B: {}".format(counts[3]))
        for i in range(len(S)):
            print(S[i])
        if isinstance(verbose, str):
            fsock.close()
        sys.stdout = saveout

    # Make sure to return elements that belong to Sage
    if return_seqs:
        return [[ZZ(_) for _ in counts],
                [[ZZ(s[0]), [QQ(_) for _ in s[1].polrecip().Vec()]] for s in S]]
    elif return_pari_objects:
        return S
    else:
        Px = PolynomialRing(QQ, 'x')
        return [[ZZ(s[0]), Px([QQ(_) for _ in s[1].list()])]
                for s in S]
