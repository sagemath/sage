"""
Arithmetic Engine for Polynomials as Tuples
"""
# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

###########
#   API   #
###########

cpdef inline tuple poly_to_tup(MPolynomial_libsingular poly):
    r"""
    Convert a polynomial object into the internal representation as tuple of
    ``(ETuple exp, NumberFieldElement coeff)`` pairs.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: R.<x, y> = PolynomialRing(QQ)
        sage: poly_to_tup(x**2 + 1)
        (((2, 0), 1), ((0, 0), 1))
        sage: poly_to_tup(x**2*y**4 - 4/5*x*y**2 + 1/3 * y)
        (((2, 4), 1), ((1, 2), -4/5), ((0, 1), 1/3))
    """
    return tuple(poly.dict().items())

cpdef inline MPolynomial_libsingular _tup_to_poly(tuple eq_tup, MPolynomialRing_libsingular parent):
    r"""
    Return a polynomial object from its tuple of pairs representation.

    Inverse of :meth:`poly_to_tup`:

    - ``poly_to_tup(tup_to_poly(eq_tup, ring)) == eq_tup`` and
    - ``tup_to_poly(poly_to_tup(eq), eq.parent()) == eq``.

    .. NOTE::

        Assumes ``all(parent.ngens() == len(exp_tup) for exp_tup, c in eq_tup)``.
        This method is meant for internal use.

    .. WARNING::

      Unsafe for client use, since it avoids implicit casting and
      it may lead to segmentation faults.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.poly_tup_engine import _tup_to_poly
        sage: K = CyclotomicField(20)
        sage: R.<x, y> = PolynomialRing(K)
        sage: poly_tup = (((2, 0), K.one()), ((0, 0), K.one()))
        sage: _tup_to_poly(poly_tup, parent=R)
        x^2 + 1
        sage: poly = x**2*y**4 - 4/5*x*y**2 + 1/3 * y
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: _tup_to_poly(poly_to_tup(poly), parent=R) == poly
        True

    TESTS::

        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup, _tup_to_poly
        sage: R.<x, y, z> = PolynomialRing(CyclotomicField(20))
        sage: r = R.random_element()
        sage: _tup_to_poly(poly_to_tup(r), parent=R) == r
        True
        sage: R = PolynomialRing(QQ, 'fx', 100)
        sage: r = R.random_element()
        sage: _tup_to_poly(poly_to_tup(r), parent=R) == r
        True
    """
    return parent._element_constructor_(dict(eq_tup), check=False)

cdef inline tuple _flatten_coeffs(tuple eq_tup):
    r"""
    Flatten cyclotomic coefficients to a representation as a tuple of rational
    coefficients.

    This is used to avoid pickling cyclotomic coefficient objects, which fails
    with new PARI settings introduced in :issue:`30537`.
    """
    cdef list flat = []
    cdef NumberFieldElement_absolute cyc_coeff
    for exp, cyc_coeff in eq_tup:
        flat.append((exp, tuple(cyc_coeff._coefficients())))
    return tuple(flat)

cpdef tuple _unflatten_coeffs(field, tuple eq_tup):
    r"""
    Restore cyclotomic coefficient object from its tuple of rational
    coefficients representation.

    Used to circumvent pickling issue introduced by PARI settings
    in :issue:`30537`.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.poly_tup_engine import _unflatten_coeffs
        sage: fm = FusionRing("A2", 2).get_fmatrix()
        sage: p = fm._poly_ring.random_element()
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: flat_poly_tup = list()
        sage: for exp, cyc_coeff in poly_to_tup(p):
        ....:     flat_poly_tup.append((exp, tuple(cyc_coeff._coefficients())))
        sage: flat_poly_tup = tuple(flat_poly_tup)
        sage: _unflatten_coeffs(fm.field(), flat_poly_tup) == poly_to_tup(p)
        True
    """
    cdef list unflat = []
    for exp, coeff_tup in eq_tup:
        unflat.append((exp, field(list(coeff_tup))))
    return tuple(unflat)

#################################
#   Useful private predicates   #
#################################

cdef inline int has_appropriate_linear_term(tuple eq_tup) noexcept:
    r"""
    Determine whether the given tuple of pairs (of length 2) contains
    an *appropriate* linear term.

    In this context, a linear term is said to be *appropriate* if
    it is in the largest variable in the given polynomial (w.r.t.
    the degrevlex ordering), the monomial in which the linear term
    appears is univariate, and the linear term is not a common factor in
    the polynomial.

    OUTPUT:

    If the given polynomial contains an appropriate linear term, this method
    returns the index of the monomial in which the term appears.

    Otherwise, the method returns `-1`.
    """
    max_var = degrees(eq_tup).nonzero_positions()[0]
    cdef ETuple m
    cdef int i
    for i in range(2):
        m = eq_tup[i][0]
        if m._nonzero == 1 and m._data[1] == 1 and m._data[0] == max_var and eq_tup[(i+1) % 2][0][max_var] == 0:
            return i
    return -1

######################
#   "Change rings"   #
######################

cpdef inline tup_to_univ_poly(tuple eq_tup, univ_poly_ring):
    r"""
    Given a tuple of pairs representing a univariate polynomial and a univariate
    polynomial ring, return a univariate polynomial object.

    Each pair in the tuple is assumed to be of the form ``(ETuple, coeff)``,
    where ``coeff`` is an element of ``univ_poly_ring.base_ring()``.

    EXAMPLES::

      sage: from sage.algebras.fusion_rings.poly_tup_engine import tup_to_univ_poly
      sage: from sage.rings.polynomial.polydict import ETuple
      sage: K = CyclotomicField(56)
      sage: poly_tup = ((ETuple([0, 3, 0]), K(2)), (ETuple([0, 1, 0]), K(-1)), (ETuple([0, 0, 0]), K(-2/3)))
      sage: R = K['b']
      sage: tup_to_univ_poly(poly_tup, R)
      2*b^3 - b - 2/3

    TESTS::

      sage: poly_tup = ((ETuple([0, 0, 0]), K(-1/5)), )
      sage: tup_to_univ_poly(poly_tup, R)
      -1/5
    """
    cdef ETuple exp
    cdef NumberFieldElement_absolute c
    return univ_poly_ring({exp._data[1] if exp._nonzero else 0: c for exp, c in eq_tup})

cpdef inline tuple resize(tuple eq_tup, dict idx_map, int nvars):
    r"""
    Return a tuple representing a polynomial in a ring with
    ``len(sorted_vars)`` generators.

    This method is used for creating polynomial objects with the
    "right number" of variables for computing Groebner bases of the
    partitioned equations graph and for adding constraints ensuring certain
    F-symbols are nonzero.

    EXAMPLES::

      sage: from sage.algebras.fusion_rings.poly_tup_engine import resize
      sage: from sage.rings.polynomial.polydict import ETuple
      sage: K = CyclotomicField(56)
      sage: poly_tup = ((ETuple([0, 3, 0, 2]), K(2)), (ETuple([0, 1, 0, 1]), K(-1)), (ETuple([0, 0, 0, 0]), K(-2/3)))
      sage: idx_map = {1: 0, 3: 1}
      sage: resize(poly_tup, idx_map, 2)
      (((3, 2), 2), ((1, 1), -1), ((0, 0), -2/3))

      sage: R = PolynomialRing(K, 'fx', 20)
      sage: R.inject_variables()
      Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13, fx14, fx15, fx16, fx17, fx18, fx19
      sage: sparse_poly = R(fx0**2 * fx17 + fx3)
      sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup, _tup_to_poly
      sage: S.<x, y, z> = PolynomialRing(K)
      sage: _tup_to_poly(resize(poly_to_tup(sparse_poly), {0:0, 3:1, 17:2}, 3), parent=S)
      x^2*z + y
    """
    cdef ETuple exp, new_e
    cdef NumberFieldElement_absolute c
    cdef list resized = list()
    for exp, c in eq_tup:
        new_e = ETuple({idx_map[pos]: d for pos, d in exp.sparse_iter()}, nvars)
        resized.append((new_e, c))
    return tuple(resized)

###########################
#   Convenience methods   #
###########################

cdef inline ETuple degrees(tuple poly_tup):
    r"""
    Return the maximal degree of each variable in the polynomial.
    """
    # Deal with the empty tuple, representing the zero polynomial
    if not poly_tup:
        return ETuple()
    cdef ETuple max_degs
    cdef int i
    max_degs = <ETuple> (<tuple> poly_tup[0])[0]
    for i in range(1, len(poly_tup)):
        max_degs = max_degs.emax(<ETuple> (<tuple> poly_tup[i])[0])
    return max_degs

cpdef list get_variables_degrees(list eqns, int nvars):
    r"""
    Find maximum degrees for each variable in equations.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.poly_tup_engine import get_variables_degrees
        sage: R.<x, y, z> = PolynomialRing(QQ)
        sage: polys = [x**2 + 1, x*y*z**2 - 4*x*y, x*z**3 - 4/3*y + 1]
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: get_variables_degrees([poly_to_tup(p) for p in polys], 3)
        [2, 1, 3]
    """
    if not eqns:
        return [0]*nvars
    cdef ETuple max_deg
    cdef int i
    max_deg = degrees(eqns[0])
    for i in range(1, len(eqns)):
        max_deg = max_deg.emax(degrees(<tuple>(eqns[i])))
    cdef list dense = [0] * len(max_deg)
    for i in range(max_deg._nonzero):
        dense[max_deg._data[2*i]] = max_deg._data[2*i+1]
    return dense

cpdef list variables(tuple eq_tup):
    """
    Return indices of all variables appearing in ``eq_tup``.

    EXAMPLES::

      sage: from sage.algebras.fusion_rings.poly_tup_engine import variables
      sage: from sage.rings.polynomial.polydict import ETuple
      sage: poly_tup = ((ETuple([0, 3, 0]), 2), (ETuple([0, 1, 0]), -1), (ETuple([0, 0, 0]), -2/3))
      sage: variables(poly_tup)
      [1]
      sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
      sage: R.<x, y, z> = PolynomialRing(QQ)
      sage: variables(poly_to_tup(x*2*y + y**3 - 4/3*x))
      [0, 1]
      sage: variables(poly_to_tup(R(1/4)))
      []
    """
    return degrees(eq_tup).nonzero_positions()

cpdef constant_coeff(tuple eq_tup, field):
    r"""
    Return the constant coefficient of the polynomial represented by
    given tuple.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.poly_tup_engine import constant_coeff
        sage: from sage.rings.polynomial.polydict import ETuple
        sage: poly_tup = ((ETuple([0, 3, 0]), 2), (ETuple([0, 1, 0]), -1), (ETuple([0, 0, 0]), -2/3))
        sage: constant_coeff(poly_tup, QQ)
        -2/3
        sage: R.<x, y, z> = PolynomialRing(QQ)
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: constant_coeff(poly_to_tup(x**5 + x*y*z - 9), QQ)
        -9
    """
    cdef ETuple exp
    for exp, coeff in eq_tup:
        if exp.is_constant():
            return coeff
    return field.zero()

cpdef tuple apply_coeff_map(tuple eq_tup, coeff_map):
    """
    Apply ``coeff_map`` to coefficients.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.poly_tup_engine import apply_coeff_map
        sage: sq = lambda x : x**2
        sage: R.<x, y, z> = PolynomialRing(ZZ)
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup, _tup_to_poly
        sage: _tup_to_poly(apply_coeff_map(poly_to_tup(x + 2*y + 3*z), sq), parent=R)
        x + 4*y + 9*z
    """
    cdef ETuple exp
    cdef list new_tup = []
    for exp, coeff in eq_tup:
        new_tup.append((exp, coeff_map(coeff)))
    return tuple(new_tup)

# cpdef inline bint tup_fixes_sq(tuple eq_tup):
cdef inline bint tup_fixes_sq(tuple eq_tup) noexcept:
    r"""
    Determine if given equation fixes the square of a variable.

    An equation fixes the sq of a variable if it is of the form `a*x^2 + c`
    for *nonzero* constants `a`, `c`.
    """
    if len(eq_tup) != 2:
        return False
    # To access _attributes, we must cdef ETuple
    cdef ETuple lm = eq_tup[0][0]
    if lm._nonzero != 1 or lm._data[1] != 2:
        return False
    cdef ETuple tm = eq_tup[1][0]
    if tm._nonzero != 0:
        return False
    return True

######################
#   Simplification   #
######################

cdef dict subs_squares(dict eq_dict, KSHandler known_sq):
    r"""
    Substitute for known squares into a given polynomial.

    INPUT:

    - ``eq_dict`` -- dictionary of ``(ETuple, coeff)`` pairs representing
      a polynomial
    - ``known_sq`` -- dictionary of ``(int i, NumberFieldElement a)`` pairs
      such that `x_i^2 - a = 0`

    OUTPUT:

    A dictionary of ``(ETuple, coeff)`` pairs.
    """
    cdef dict subbed, new_e
    cdef ETuple exp
    cdef int idx, power
    subbed = dict()
    for exp, coeff in eq_dict.items():
        new_e = dict()
        for idx, power in exp.sparse_iter():
            if known_sq.contains(idx):
                coeff *= pow(known_sq.get(idx), power // 2)
                # New power is 1 if power is odd
                if power & True:
                    new_e[idx] = 1
            else:
                new_e[idx] = power
        exp = ETuple(new_e, len(exp))
        # If exponent tuple is already present in dictionary, coefficients are added
        if exp in subbed:
            subbed[exp] += coeff
        else:
            subbed[exp] = coeff
    return subbed

cdef dict remove_gcf(dict eq_dict, ETuple nonz):
    r"""
    Return a dictionary of ``(ETuple, coeff)`` pairs describing the
    polynomial ``eq / GCF(eq)``.

    The input ``nonz`` is an ``ETuple`` indicating the positions of
    variables known to be nonzero. The entries of ``nonz`` are assumed to
    be some relatively large number, like 100.
    """
    # Find common variables, filtered according to known nonzeros
    cdef ETuple common_powers, exp
    cdef NumberFieldElement_absolute c
    common_powers = nonz
    for exp, c in eq_dict.items():
        common_powers = common_powers.emin(exp)
    cdef dict ret = {}
    for exp, c in eq_dict.items():
        ret[exp.esub(common_powers)] = c
    return ret

cdef tuple to_monic(dict eq_dict, one):
    """
    Return tuple of pairs ``(ETuple, coeff)`` describing the monic polynomial
    associated to ``eq_dict``.

    Here, the leading coefficient is chosen according to the degree reverse
    lexicographic ordering (default for multivariate polynomial rings).
    """
    if not eq_dict:
        return ()
    cdef list ord_monoms = sorted(eq_dict, key=monom_sortkey)
    cdef ETuple lm = ord_monoms[-1]
    cdef NumberFieldElement_absolute lc = eq_dict[lm]
    if not lc:
        return ()
    cdef list ret = [(lm, one)]
    inv_lc = lc.inverse_of_unit()
    cdef int i, n
    n = len(ord_monoms)
    for i in range(n-1):
        ret.append((ord_monoms[n-2-i], inv_lc * eq_dict[ord_monoms[n-2-i]]))
    return tuple(ret)

cdef tuple reduce_poly_dict(dict eq_dict, ETuple nonz, KSHandler known_sq, NumberFieldElement_absolute one):
    """
    Return a tuple describing a monic polynomial with no known nonzero
    gcf and no known squares.
    """
    if not eq_dict:
        return ()
    cdef dict sq_rmvd = subs_squares(eq_dict, known_sq)
    cdef dict gcf_rmvd = remove_gcf(sq_rmvd, nonz)
    return to_monic(gcf_rmvd, one)

####################
#   Substitution   #
####################

cpdef dict compute_known_powers(max_degs, dict val_dict, one):
    """
    Pre-compute powers of known values for efficiency when preparing to
    substitute into a list of polynomials.

    INPUT:

    - ``max_deg`` -- an ``ETuple`` indicating the maximal degree of
      each variable
    - ``val_dict`` -- dictionary of ``(var_idx, poly_tup)`` key-value pairs
    - ``poly_tup`` -- tuple of ``(ETuple, coeff)`` pairs representing a
      multivariate polynomial

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.poly_tup_engine import compute_known_powers
        sage: R.<x, y, z> = PolynomialRing(QQ)
        sage: polys = [x**3 + 1, x**2*y + z**3, y**2 - 3*y]
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: known_val = { 0 : poly_to_tup(R(-1)), 2 : poly_to_tup(y**2) }
        sage: from sage.algebras.fusion_rings.poly_tup_engine import get_variables_degrees
        sage: max_deg = get_variables_degrees([poly_to_tup(p) for p in polys], 3)
        sage: compute_known_powers(max_deg, known_val, R.base_ring().one())
        {0: [(((0, 0, 0), 1),),
             (((0, 0, 0), -1),),
             (((0, 0, 0), 1),),
             (((0, 0, 0), -1),)],
        2: [(((0, 0, 0), 1),),
            (((0, 2, 0), 1),),
            (((0, 4, 0), 1),),
            (((0, 6, 0), 1),)]}
    """
    assert (max_degs and max(max_degs) <= 100) or True, "cannot substitute for degree larger than 100"
    cdef ETuple max_deg = ETuple(list(max_degs))
    max_deg = max_deg.emin(ETuple({idx: 100 for idx in val_dict}, len(max_deg)))
    cdef dict known_powers
    # Get polynomial unit as tuple to initialize list elements
    cdef tuple one_tup = ((max_deg._new(), one), )
    cdef int d, power, var_idx
    known_powers = {var_idx: [one_tup]*(d+1) for var_idx, d in max_deg.sparse_iter()}
    for var_idx, d in max_deg.sparse_iter():
        for power in range(d):
            known_powers[var_idx][power+1] = tup_mul(known_powers[var_idx][power], val_dict[var_idx])
    return known_powers

cdef dict subs(tuple poly_tup, dict known_powers, one):
    """
    Substitute given variables into a polynomial tuple.
    """
    cdef dict subbed = {}
    cdef ETuple exp, m, shifted_exp
    cdef int var_idx, power
    cdef tuple temp
    for exp, coeff in poly_tup:
        # Get polynomial unit as tuple
        temp = ((exp._new(), one), )
        for var_idx, power in exp.sparse_iter():
            if var_idx in known_powers:
                exp = exp.eadd_p(-power, var_idx)
                temp = tup_mul(temp, known_powers[var_idx][power])
        for m, c in temp:
            shifted_exp = exp.eadd(m)
            if shifted_exp in subbed:
                subbed[shifted_exp] += coeff * c
            else:
                subbed[shifted_exp] = coeff * c
    return subbed

cdef tuple tup_mul(tuple p1, tuple p2):
    r"""
    Multiplication of two polynomial tuples using schoolbook multiplication.
    """
    cdef dict prod = {}
    cdef ETuple xi, yj, shifted_exp
    for xi, ai in p1:
        for yj, bj in p2:
            shifted_exp = xi.eadd(yj)
            if shifted_exp in prod:
                prod[shifted_exp] += ai * bj
            else:
                prod[shifted_exp] = ai * bj
    return tuple(prod.items())

###############
#   Sorting   #
###############

cdef tuple monom_sortkey(ETuple exp):
    r"""
    Produce a sortkey for a monomial exponent with respect to degree
    reversed lexicographic ordering.
    """
    cdef int deg = exp.unweighted_degree()
    # for i in range(exp._nonzero):
    # exp._data[2*i+1] = -exp._data[2*i+1]
    cdef ETuple rev = exp.reversed().emul(-1)
    return (deg, rev)

cpdef tuple poly_tup_sortkey(tuple eq_tup):
    r"""
    Return the sortkey of a polynomial represented as a tuple of
    ``(ETuple, coeff)`` pairs with respect to the degree
    lexicographical term order.

    Using this key to sort polynomial tuples results in comparing polynomials
    term by term (we assume the tuple representation is sorted so that the
    leading term with respect to the degree reverse lexicographical order
    comes first). For each term, we first compare degrees, then the monomials
    themselves. Different polynomials can have the same sortkey.

    EXAMPLES::

     sage: F = CyclotomicField(20)
     sage: zeta20 = F.gen()
     sage: R.<x, y, z> = PolynomialRing(F)
     sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_tup_sortkey, poly_to_tup
     sage: p = (zeta20 + 1)*x^2 + (zeta20^3 + 6)*x*z + (zeta20^2 + 7*zeta20)*z^2 + (2/3*zeta20 + 1/4)*x + y
     sage: p1 = poly_to_tup(p); p1
     (((2, 0, 0), zeta20 + 1),
      ((1, 0, 1), zeta20^3 + 6),
      ((0, 0, 2), zeta20^2 + 7*zeta20),
      ((1, 0, 0), 2/3*zeta20 + 1/4),
      ((0, 1, 0), 1))
        sage: poly_tup_sortkey(p1)
        (2, 0, 2, 2, 0, 1, -2, 1, 2, -2, 2, 1, 0, 1, 1, -1, 1)
    """
    cdef ETuple exp
    cdef int i
    cdef list key = []
    for exp, c in eq_tup:
        # Compare by term degree
        key.append(exp.unweighted_degree())
        # Next compare by term w.r.t. lex order
        for i in range(exp._nonzero):
            # key.append(exp._length-1-exp._data[2*(nnz-i-1)])
            # key.append(-exp._data[2*(nnz-i-1)+1])
            key.append(-exp._data[2*i])
            key.append(exp._data[2*i+1])
    return tuple(key)
