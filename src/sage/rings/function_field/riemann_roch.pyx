# distutils: language = c++
r"""
Utility functions for Riemann-Roch computations with divisors represented as pairs of ideals.

AUTHORS:

- Kwankyu Lee (2017-2022): initial versions
- Vincent Macri (2025-10-29): Cythonization and refactoring
"""
# ****************************************************************************
#       Copyright (C) 2017-2022 Kwankyu Lee <ekwankyu@gmail.com>
#                     2025      Vincent Macri <vincent.macri@ucalgary.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.list cimport list as linked_list

from sage.arith.functions import lcm  # Should we use LCM_list or make lcm cpdef?
from sage.matrix.constructor import matrix
from sage.matrix.matrix cimport Matrix

# Used to avoid a bunch of Python overhead.
# Should be unnecessary if https://github.com/cython/cython/issues/3679 is fixed
ctypedef pair[int, int] ipair

cpdef _divisor_to_inverted_ideals(divisor):
    """
    Return inverted ideal representation of divisor.

    TESTS::

        sage: from sage.rings.function_field import riemann_roch
        sage: k = GF(17)
        sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: pl = C([2,8,1]).place()
        sage: dS, ds = riemann_roch._divisor_to_inverted_ideals(2*pl)
        sage: (~dS).divisor() + (~ds).divisor() == 2*pl
        True
    """
    F = divisor.parent()._field
    O = F.maximal_order()
    Oinf = F.maximal_order_infinite()

    # The ideal I is the inverse of the product of prime ideals attached
    # to the finite places in the divisor while the ideal J corresponds
    # to the infinite places in the divisor.
    I = O.unit_ideal()
    J = Oinf.unit_ideal()

    for p, m in divisor._data.items():
        if p.is_infinite_place():
            J *= pow(p.prime_ideal(), -m)
        else:
            I *= pow(p.prime_ideal(), -m)

    return I, J

def _riemann_roch_divisor(divisor):
    r"""
    Wrapper around :func:`_riemann_roch_ideals` that first
    converts the inputted divisor to a pair of inverted ideals.
    """
    # The function _riemann_roch_ideals is basically to compute
    # the intersection of the ideals I and J.
    I, J = _divisor_to_inverted_ideals(divisor)
    return _riemann_roch_ideals(divisor.parent()._field, I, J)

cpdef _riemann_roch_ideals(F, I, J):
    r"""
    Wrapper around :func:`_riemann_roch_basis` that takes a function field
    and a pair of (inverted) ideals.
    """
    _, _, to = F.free_module(map=True)
    gens = list(I.gens_over_base())
    C = matrix([to(v) for v in gens])
    B = matrix([to(b) for b in J.gens_over_base()])
    return _riemann_roch_basis(F.degree(), C, B.inverse(), gens, False)


def _short_circuit_riemann_roch_ideals(F, I, J):
    r"""
    Wrapper around :func:`_riemann_roch_basis` that takes a function field
    and a pair of (inverted) ideals. Performs a short-circuited Riemann-Roch
    computation that returns after finding a single basis element.
    """
    _, _, to = F.free_module(map=True)
    gens = list(I.gens_over_base())
    C = matrix([to(v) for v in gens])
    B = matrix([to(b) for b in J.gens_over_base()])
    return _riemann_roch_basis(F.degree(), C, B.inverse(), gens, True)


def _short_circuit_riemann_roch_matrices(int n, Matrix C, Matrix B_inv, gens):
    r"""
    Wrapper around :func:`_riemann_roch_basis` where matrix conversion is
    already done. Used by Unique Hess implementation which can do some
    caching with the computation of ``B_inv``.

    INPUT:

    - ``n`` -- degree of the function field
    - ``C`` -- matrix representation of finite ideal
    - ``B_inv`` -- matrix inverse of matrix representation of infinite ideal
    - ``gens`` -- list of generators for the finite ideal
    """
    return _riemann_roch_basis(n, C, B_inv, list(gens), True)


cdef _riemann_roch_basis(int n, Matrix C, Matrix B_inv, list gens, bint short_circuit):
    """
    Return a basis of the Riemann-Roch space of the divisor.

    Modifies gens


    This implements Hess' algorithm 6.1 in [Hes2002]_
    """
    # The basic idea of this algorithm is to compute the intersection of the ideals I and J.

    # Step 1: construct matrix M of rational functions in x such that
    # M * B == C where B = [b1,b1,...,bn], C =[v1,v2,...,vn]
    M = C * B_inv

    # Step 1.5: get the denominator d of M and set mat = d * M
    den = lcm([e.denominator() for e in M.list()])
    R = den.parent()  # polynomial ring
    one = R.one()
    mat = matrix(R, n, [e.numerator() for e in (den * M).list()])

    # Step 2: transform mat to a weak Popov form, together with gens

    # initialise pivot_row and conflicts list
    cdef vector[vector[ipair]] pivot_row = vector[vector[ipair]](n)
    cdef vector[int] conflicts = vector[int](0)
    conflicts.reserve(n)

    cdef int i, j, bestp, best, c, d, ideg, jdeg
    for i in range(n):
        bestp = -1
        best = -1
        for c in range(n):
            d = mat[i, c].degree()
            if d >= best:
                bestp = c
                best = d

        if short_circuit and best <= den.degree():
            return gens[i]

        if best >= 0:
            pivot_row[bestp].push_back(ipair(i, best))

            if pivot_row[bestp].size() > 1:
                conflicts.push_back(bestp)

    # while there is a conflict, do a simple transformation
    while not conflicts.empty():
        c = conflicts.back()
        conflicts.pop_back()
        row = pivot_row[c]

        tmp = pivot_row[c].back()
        i, ideg = tmp.first, tmp.second
        pivot_row[c].pop_back()

        tmp = pivot_row[c].back()
        j, jdeg = tmp.first, tmp.second
        pivot_row[c].pop_back()

        if jdeg > ideg:
            i, j = j, i
            ideg, jdeg = jdeg, ideg

        coeff = - mat[i, c].lc() / mat[j, c].lc()
        s = coeff * one.shift(ideg - jdeg)

        mat.add_multiple_of_row(i, j, s)
        gens[i] += s * gens[j]

        pivot_row[c].push_back(ipair(j, jdeg))

        bestp = -1
        best = -1
        for c in range(n):
            d = mat[i, c].degree()
            if d >= best:
                bestp = c
                best = d

        if short_circuit and best <= den.degree():
            return gens[i]

        if best >= 0:
            pivot_row[bestp].push_back(ipair(i, best))
            if pivot_row[bestp].size() > 1:
                conflicts.push_back(bestp)

    if short_circuit:
        return None

    # Step 3: build a Riemann-Roch basis from the data in mat and gens.
    # Note that the values mat[i,j].degree() - den.degree() are known as
    # invariants of M.
    basis = []
    for j in range(n):
        i, ideg = pivot_row[j][0].first, pivot_row[j][0].second
        gi = gens[i]
        basis.extend(one.shift(k) * gi
                     for k in range(den.degree() - ideg + 1))
    # Done!
    return basis
