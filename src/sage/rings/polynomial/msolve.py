r"""
Solution of polynomial systems using msolve

`msolve <https://msolve.lip6.fr/>`_ is a multivariate polynomial system solver
based on Gröbner bases.

This module provide implementations of some operations on polynomial ideals
based on msolve. :

It additionally provides a function for computing sample points
per connected components of semi-algebraic sets defined by a single inequality
or inequation. Note that it does not guarantee uniqueness; in particular, there
can be multiple points in the output belonging to the same connected component.

Note that the :ref:`optional package msolve <spkg_msolve>` must be installed.

.. SEEALSO::

    - :mod:`sage.features.msolve`
    - :mod:`sage.rings.polynomial.multi_polynomial_ideal`

AUTHORS:
- Marc Mezzarobba (2022) -- initial version
- Edern Gillot (2026) -- sample points per connected components of semi-algebraic sets
"""

# ****************************************************************************
#       Copyright (C) 2022 Marc Mezzarobba
#                     2026 Edern Gillot
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import tempfile
import subprocess
import itertools

import sage.structure.proof.proof

from sage.features.msolve import msolve
from sage.misc.converting_dict import KeyConvertingDict
from sage.misc.sage_eval import sage_eval
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.rational_field import QQ
from sage.rings.real_arb import RealBallField
from sage.rings.real_double import RealDoubleField_class
from sage.rings.real_mpfr import RealField_class
from sage.rings.real_mpfi import RealIntervalField_class, RealIntervalField
from sage.structure.sequence import Sequence
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.functions.generalized import sign
from sage.functions.other import floor, ceil


def _run_msolve(ideal, options):
    r"""
    Internal utility function
    """

    base = ideal.base_ring()
    if not (base is QQ or isinstance(base, FiniteField) and
            base.is_prime_field() and base.characteristic() < 2**31):
        raise NotImplementedError(f"unsupported base field: {base}")

    # Run msolve

    drlpolring = ideal.ring().change_ring(order='degrevlex')
    polys = ideal.change_ring(drlpolring).gens()
    msolve_in = tempfile.NamedTemporaryFile(mode='w',
                                            encoding='ascii', delete=False)
    command = [msolve().absolute_filename(), "-f", msolve_in.name] + options
    try:
        print(",".join(drlpolring.variable_names()), file=msolve_in)
        print(base.characteristic(), file=msolve_in)
        print(*(pol._repr_().replace(" ", "") for pol in polys),
                sep=',\n', file=msolve_in)
        msolve_in.close()
        msolve_out = subprocess.run(command, capture_output=True, text=True)
    finally:
        os.unlink(msolve_in.name)
    msolve_out.check_returncode()

    return msolve_out.stdout


def groebner_basis_degrevlex(ideal, proof=True):
    r"""
    Compute a degrevlex Gröbner basis using msolve

    EXAMPLES::

        sage: from sage.rings.polynomial.msolve import groebner_basis_degrevlex

        sage: R.<a,b,c> = PolynomialRing(GF(101), 3, order='lex')
        sage: I = sage.rings.ideal.Katsura(R,3)
        sage: gb = groebner_basis_degrevlex(I); gb # optional - msolve
        [a + 2*b + 2*c - 1, b*c - 19*c^2 + 10*b + 40*c,
        b^2 - 41*c^2 + 20*b - 20*c, c^3 + 28*c^2 - 37*b + 13*c]
        sage: gb.universe() is R # optional - msolve
        False
        sage: gb.universe().term_order() # optional - msolve
        Degree reverse lexicographic term order
        sage: ideal(gb).transformed_basis(other_ring=R) # optional - msolve
        [c^4 + 38*c^3 - 6*c^2 - 6*c, 30*c^3 + 32*c^2 + b - 14*c,
        a + 2*b + 2*c - 1]

    Gröbner bases over the rationals require `proof=False`::

        sage: R.<x, y> = PolynomialRing(QQ, 2)
        sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])
        sage: I.groebner_basis(algorithm='msolve') # optional - msolve
        Traceback (most recent call last):
        ...
        ValueError: msolve relies on heuristics; please use proof=False
        sage: I.groebner_basis(algorithm='msolve', proof=False) # optional - msolve
        [x*y - 1, x^2 + y^2 - 4*x - 2*y + 4, y^3 - 2*y^2 + x + 4*y - 4]

    TESTS::

        sage: R.<foo, bar> = PolynomialRing(GF(536870909), 2)
        sage: I = Ideal([ foo^2 - 1, bar^2 - 1 ])
        sage: I.groebner_basis(algorithm='msolve') # optional - msolve
        [bar^2 - 1, foo^2 - 1]
    """

    if ideal.base_ring() is QQ and sage.structure.proof.proof.get_flag(proof, "polynomial"):
        raise ValueError("msolve relies on heuristics; please use proof=False")

    drlpolring = ideal.ring().change_ring(order='degrevlex')
    msolve_out = _run_msolve(ideal, ["-g", "2"])
    gbasis = sage_eval(msolve_out[:-2], locals=drlpolring.gens_dict())
    return Sequence(gbasis)


def variety(ideal, ring, *, proof=True):
    r"""
    Compute the variety of a zero-dimensional ideal using msolve.

    Part of the initial implementation was loosely based on the example
    interfaces available as part of msolve, with the authors' permission.

    EXAMPLES::

        sage: from sage.rings.polynomial.msolve import variety
        sage: p = 536870909
        sage: R.<x, y> = PolynomialRing(GF(p), 2, order='lex')
        sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])
        sage: sorted(variety(I, GF(p^2), proof=False), key=lambda d: str(sorted(d.items()))) # optional - msolve
        [{x: 1, y: 1},
         {x: 254228855*z2 + 114981228, y: 232449571*z2 + 402714189},
         {x: 267525699, y: 473946006},
         {x: 282642054*z2 + 154363985, y: 304421338*z2 + 197081624}]

    TESTS::

        sage: p = 536870909
        sage: R.<x, y> = PolynomialRing(GF(p), 2, order='lex')
        sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])

        sage: sorted(I.variety(algorithm='msolve', proof=False), key=lambda d: str(sorted(d.items()))) # optional - msolve
        [{x: 1, y: 1}, {x: 267525699, y: 473946006}]

        sage: K.<a> = GF(p^2)
        sage: sorted(I.variety(K, algorithm='msolve', proof=False), key=lambda d: str(sorted(d.items()))) # optional - msolve
        [{x: 1, y: 1},
         {x: 118750849*a + 194048031, y: 510295713*a + 18174854},
         {x: 267525699, y: 473946006},
         {x: 418120060*a + 75297182, y: 26575196*a + 44750050}]

        sage: R.<x, y> = PolynomialRing(GF(2147483659), 2, order='lex')
        sage: ideal([x, y]).variety(algorithm='msolve', proof=False)
        Traceback (most recent call last):
        ...
        NotImplementedError: unsupported base field: Finite Field of size 2147483659

        sage: R.<x, y> = PolynomialRing(QQ, 2, order='lex')
        sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])

        sage: I.variety(algorithm='msolve', proof=False) # optional - msolve
        [{x: 1, y: 1}]
        sage: I.variety(RealField(100), algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.7692923542386314152404094643, y: 0.36110308052864737763464656216},
         {x: 1.0000000000000000000000000000, y: 1.0000000000000000000000000000}]
        sage: I.variety(RealIntervalField(100), algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.76929235423863141524040946434?, y: 0.361103080528647377634646562159?},
         {x: 1, y: 1}]
        sage: I.variety(RBF, algorithm='msolve', proof=False) # optional - msolve
        [{x: [2.76929235423863 +/- 2.08e-15], y: [0.361103080528647 +/- 4.53e-16]},
         {x: 1.000000000000000, y: 1.000000000000000}]
        sage: I.variety(RDF, algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.7692923542386314, y: 0.36110308052864737}, {x: 1.0, y: 1.0}]
        sage: I.variety(AA, algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.769292354238632?, y: 0.3611030805286474?},
         {x: 1.000000000000000?, y: 1.000000000000000?}]
        sage: I.variety(QQbar, algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.769292354238632?, y: 0.3611030805286474?},
         {x: 1, y: 1},
         {x: 0.11535382288068429? + 0.5897428050222055?*I, y: 0.3194484597356763? - 1.633170240915238?*I},
         {x: 0.11535382288068429? - 0.5897428050222055?*I, y: 0.3194484597356763? + 1.633170240915238?*I}]
        sage: I.variety(ComplexField(100))
        [{y: 1.0000000000000000000000000000, x: 1.0000000000000000000000000000},
         {y: 0.36110308052864737763464656216, x: 2.7692923542386314152404094643},
         {y: 0.31944845973567631118267671892 - 1.6331702409152376561188467320*I, x: 0.11535382288068429237979526783 + 0.58974280502220550164728074602*I},
         {y: 0.31944845973567631118267671892 + 1.6331702409152376561188467320*I, x: 0.11535382288068429237979526783 - 0.58974280502220550164728074602*I}]

        sage: Ideal(x^2 + y^2 - 1, x - y).variety(RBF, algorithm='msolve', proof=False) # optional - msolve
        [{x: [-0.707106781186547 +/- 6.29e-16], y: [-0.707106781186547 +/- 6.29e-16]},
         {x: [0.707106781186547 +/- 6.29e-16], y: [0.707106781186547 +/- 6.29e-16]}]
        sage: sorted(Ideal(x^2 - 1, y^2 - 1).variety(QQ, algorithm='msolve', proof=False), key=lambda d: str(sorted(d.items()))) # optional - msolve
        [{x: -1, y: -1}, {x: -1, y: 1}, {x: 1, y: -1}, {x: 1, y: 1}]
        sage: Ideal(x^2-1, y^2-2).variety(CC, algorithm='msolve', proof=False) # optional - msolve
        [{x: 1.00000000000000, y: 1.41421356237310},
         {x: -1.00000000000000, y: 1.41421356237309},
         {x: 1.00000000000000, y: -1.41421356237309},
         {x: -1.00000000000000, y: -1.41421356237310}]

        sage: Ideal([x, y, x + y]).variety(algorithm='msolve', proof=False) # optional - msolve
        [{x: 0, y: 0}]

        sage: Ideal([x, y, x + y - 1]).variety(algorithm='msolve', proof=False) # optional - msolve
        []
        sage: Ideal([x, y, x + y - 1]).variety(RR, algorithm='msolve', proof=False) # optional - msolve
        []

        sage: Ideal([x*y - 1]).variety(QQbar, algorithm='msolve', proof=False) # optional - msolve
        Traceback (most recent call last):
        ...
        ValueError: positive-dimensional ideal

        sage: R.<x, y> = PolynomialRing(RR, 2, order='lex')
        sage: Ideal(x, y).variety(algorithm='msolve', proof=False)
        Traceback (most recent call last):
        ...
        NotImplementedError: unsupported base field: Real Field with 53 bits of precision

        sage: R.<x, y> = PolynomialRing(QQ, 2, order='lex')
        sage: Ideal(x, y).variety(ZZ, algorithm='msolve', proof=False)
        Traceback (most recent call last):
        ...
        ValueError: no coercion from base field Rational Field to output ring Integer Ring
    """

    proof = sage.structure.proof.proof.get_flag(proof, "polynomial")
    if proof:
        raise ValueError("msolve relies on heuristics; please use proof=False")

    base = ideal.base_ring()
    if ring is None:
        ring = base
    if not ring.has_coerce_map_from(base):
        raise ValueError(
            f"no coercion from base field {base} to output ring {ring}")

    if isinstance(ring, (RealIntervalField_class, RealBallField,
                         RealField_class, RealDoubleField_class)):
        parameterization = False
        options = ["-p", str(ring.precision())]
    else:
        parameterization = True
        options = ["-P", "1"]

    msolve_out = _run_msolve(ideal, options)

    # Interpret output

    try:
        data = sage_eval(msolve_out[:-2])
    except SyntaxError:
        raise NotImplementedError(f"unsupported msolve output format: {data}")

    dim = data[0]
    if dim == -1:
        return []
    elif dim > 0:
        raise ValueError("positive-dimensional ideal")
    else:
        assert dim.is_zero()

    out_ring = ideal.ring().change_ring(ring)

    if parameterization:

        def to_poly(p, d=1, *, upol=PolynomialRing(base, 't')):
            assert len(p[1]) == p[0] + 1 or p == [-1, [0]]
            return upol(p[1])/d

        try:
            char, nvars, deg, vars, _, [one, [elim, den, param]] = data[1]
        except (IndexError, ValueError):
            raise NotImplementedError(
                f"unsupported msolve output format: {data}")
        assert char == ideal.base_ring().characteristic()
        assert one.is_one()
        assert len(vars) == nvars
        ringvars = out_ring.variable_names()
        assert sorted(vars[:len(ringvars)]) == sorted(ringvars)
        vars = [out_ring(name) for name in vars[:len(ringvars)]]
        elim = to_poly(elim)
        # Criterion suggested by Mohab Safey El Din to avoid cases where there
        # is no rational parameterization or where the one returned by msolve
        # has a significant probability of being incorrect.
        if deg >= char > 0 or 0 < char <= 2**17 and deg != elim.degree():
            raise NotImplementedError(f"characteristic {char} too small")
        den = to_poly(den)
        # As of msolve 0.4.4, param is of the form [pol, denom] in char 0, but
        # [pol] in char p > 0. My understanding is that both cases will
        # eventually use the same format, so let's not be too picky.
        param = [to_poly(*f) for f in param]
        elim_roots = elim.roots(ring, multiplicities=False)
        variety = []
        for rt in elim_roots:
            den_of_rt = den(rt)
            point = [-p(rt) / den_of_rt for p in param]
            if len(param) != len(vars):
                point.append(rt)
            assert len(point) == len(vars)
            variety.append(point)

    else:

        if len(data[1]) < 2 or len(data[1]) != data[1][0] + 1:
            raise NotImplementedError(
                f"unsupported msolve output format: {data}")
        if isinstance(ring, (RealIntervalField_class, RealBallField)):
            to_out_ring = ring
        else:
            assert isinstance(ring, (RealField_class, RealDoubleField_class))
            myRIF = RealIntervalField(ring.precision())
            to_out_ring = lambda iv: ring.coerce(myRIF(iv).center())
        vars = out_ring.gens()
        variety = [[to_out_ring(iv) for iv in point]
                   for l in data[1][1:]
                   for point in l]

    return [KeyConvertingDict(out_ring, zip(vars, point)) for point in variety]

def _format_output_msolve_grobner(ms_output):
    r"""
    Internal utility function:

    Converts a msolve grobner basis format string into a list
    """

    sols = []
    is_sol_reached = False
    for l2 in ms_output.splitlines():
        if l2 == '#Leading ideal data':
            is_sol_reached = True
        if is_sol_reached and l2 != '' and l2[0] == '-':
            is_sol_reached = False
        if not is_sol_reached:
            print(l2)
            continue
        l = ''
        for c in l2:
            if c not in ['[', ']', '\n', ',', ':']:
                l += c
        if l == '' or l[0] == "#":
            continue
        sols.append(l)
    return sols

def _format_output_msolve_intervals(ms_output):
    r"""
    Internal utility function:

    Converts a msolve isolation intervals format string into a list
    """

    sols = ""
    is_sol_reached = False
    for l2 in ms_output.splitlines():
        if l2.startswith('[0,') or l2.startswith('[1,') or l2.startswith('[-1'):
            is_sol_reached = True
        if is_sol_reached and l2 != '' and l2[0] == '-':
            is_sol_reached = False
        if not is_sol_reached:
            print(l2)
            continue
        sols += l2
    sols = sols.replace("\n", "").replace(":", "")
    return sage_eval(sols)

def _is_smooth(poly, threads, msolve_verbose):
    r"""
    Function that checks for the smoothness of V(poly) where poly has
    rational coefficients.

    INPUT:

        - ``poly`` -- polynomial with rational coefficients

        - ``threads`` -- integer; number of threads to be used by msolve in
        computation

        - ``msolve_verbose`` -- 0, 1 or 2; msolve parameter for explicit
        description of computations, ranging from 0 (no desc) to 2 (full desc)

    OUTPUT:

        - boolean; answer to the question "Is V(poly) smooth?"

    EXAMPLES::

        sage: from sage.rings.polynomial.msolve import _is_smooth
        sage: R.<x,y> = QQ[]
        sage: f = x^2 + y^2
        sage: _is_smooth(f,1,0) # optional - msolve
        False

    ::

        sage: from sage.rings.polynomial.msolve import _is_smooth
        sage: R.<x,y,z> = QQ[]
        sage: f = x^2 + y^2 + z^2 - 1
        sage: _is_smooth(f,1,0) # optional - msolve
        True
    """

    R = poly.parent()
    variables = list(R.gens())
    input_list = [poly.derivative(variables[i]) for i in range(len(variables))]
    input_list.append(poly)
    ms_out = sage.rings.polynomial.msolve._run_msolve(R.ideal(input_list), ['-g','1','-v',f"{msolve_verbose}",'-t',f"{threads}"])
    gb = _format_output_msolve_grobner(ms_out)
    if gb == ['1']:
        return True
    return False

def _grp_random_matrix(n, changevar=True):
    r"""
    Internal Function

    Function that generates a random (n x n) change of variables matrix, such
    that all partial inverses B_k exist, for 0 <= k <= n.
    Under paper notation: such that it satisfies hypothesis (A2).

    INPUT:

        - ``n`` -- integer; number of variables

        - ``changevar`` -- boolean (default ``True``); uses A = Identity if
        set to ``False``

    OUTPUT:

        - list of matrices; format [A, B_0 (= A^{-1}), B_1, ..., B_n]
    """

    while True:
        try:
            if changevar:
                A = matrix([[ZZ.random_element(1,100,"uniform") for j in range(n)] for i in range(n)])
            else:
                A = matrix.identity(n)
            list_of_matrices = [A] + [A[list(range(k,n)), list(range(k,n))].inverse() for k in range(n)]
            return list_of_matrices
        except ZeroDivisionError:
            pass

def _derivative_order(poly):
    r"""
    Internal Function

    Function that re-labels the variables of the input polynomial such that its
    partial derivatives have increasing degree.

    INPUT:

        - ``poly`` -- polynomial with rational coefficients

    OUTPUT:

        - polynomial; re-labeled polynomial with partial derivatives of
        increasing degree
    """

    from sage.combinat.words.word import Word
    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    from sage.symbolic.ring import SR
    R = poly.parent()
    variables = list(R.gens())
    n = len(poly.variables())
    der_deg = {variables[i] : poly.derivative(variables[i]).degree() for i in range(n)}
    new_var = sorted(der_deg, key=lambda k: der_deg[k])
    perm = Word(variables).standard_permutation() / Word(new_var).standard_permutation()
    inv_perm = Word(new_var).standard_permutation() / Word(variables).standard_permutation()
    sigma = (SymmetricGroup(range(n)))([perm[i]-1 for i in range(n)])
    inv_sigma = (SymmetricGroup(range(n)))([inv_perm[i]-1 for i in range(n)])
    return(inv_sigma, R(SR(poly(*sigma(R.gens())))))

def _remove_absent_variable(P, x):
    r"""
    Internal Function

    Function that efficiently changes the parent ring of a polynomial to one
    with less variables.

    INPUT:

        - ``P`` -- multivariate polynomial with rational coefficients;
        not tested for univariate polynomials.

        - ``x`` -- tuple or list; list of variables such that
        each variable appears in P.parent() and not in P itself.

    OUTPUT:

        - polynomial; P, but such that P.parent() does not contain x anymore.
    """

    from sage.rings.polynomial.polydict import PolyDict, ETuple

    # if x is just a variable, or is a list, make it a tuple
    if isinstance(x, list):
        x = tuple(x)
    if not isinstance(x, tuple):
        x = (x,)

    # parent ring and ring with variables removed
    pring = P.parent()
    ring = pring.remove_var(*x)

    # get positions of variables to keep
    var_idx = []
    y = pring.gens()
    for i in range(pring.ngens()):
        yi = y[i]
        if yi not in x:
            var_idx.append(i)

    # retrieve dictionary {monomial : coeff}
    Pdic = P.monomial_coefficients()

    # remove variables
    Qdic = {ETuple([et[i] for i in var_idx]) : coeff for et, coeff in Pdic.items()}
    Q = ring(Qdic)

    return Q

def _critical_points(f, threads, msolve_verbose, precision, k, n, list_of_matrices, sigma, variables, der_list):
    r"""
    Internal Function

    Function computing the critical points of projection of V(f) on X_k-axis.

    INPUT:

        - ``f`` -- polynomial with rational coefficients

        - ``threads`` -- integer; number of threads to be used by msolve in
        computation

        - ``msolve_verbose`` -- 0, 1 or 2; msolve parameter for explicit
        description of computations, ranging from 0 (no desc) to 2 (full desc)

        - ``precision`` -- integer; number of bits of precision used by msolve
        for real root approximation

        - ``k`` -- integer; integer satisfying 0 <= k < n

        - ``n`` -- integer; total number of variables

        - ``list_of_matrices`` -- list of matrices; change of variables
        matrix A and its B_k's, formatted as in the output of _grp_random_matrix

        - ``sigma`` -- list of integers; list of n-1 integers

        - ``variables`` -- list; list of all variables

        - ``der_list`` -- list of polynomials; list of partial derivatives of f
        with respect to every element in `variables`

    OUTPUT:

        - list of k-1 polynomials; substitution expressions for the
        first k-1 variables

        - Isolation intervals with rational endpoints for the last (n-k)
        coordinates of A applied to each critical point of the projection of
        V^A on the X_k-axis, where the first coordinates have been
        instantiated to sigma[0], ..., sigma[k-1]. This follows the msolve
        output format for the -P0 flag.
    """

    # Computing the actual values to substitue into x_0, ..., x_{k-1}
    if k == 0:
        substitution = []
    else:
        inv_left = list_of_matrices[1][list(range(k)), list(range(k))].inverse()
        right = list_of_matrices[1][list(range(k)), list(range(k,n))]
        sigma_k = matrix(sigma[:k]).transpose()
        vars_k = matrix(variables[k:]).transpose()
        substitution = (inv_left * (sigma_k - (right * vars_k))).coefficients()

    # Computing the system equivalent to f^A, df^A/dX_{k+1}, ..., df^A/dX_n.
    if k == n-1:
        der_system = []
    else:
        left = matrix(der_list[k+1:])
        top_right = (list_of_matrices[0])[list(range(k+1)), list(range(k+1,n))]
        inv_right = list_of_matrices[k+2]
        der_system = (left + (matrix(der_list[:k+1]) * top_right * inv_right)).coefficients()
    input_system = [f] + der_system

    # Subsituting the former in the latter, with the right parent ring.
    input_system = [p.subs({variables[i] : substitution[i] for i in range(k)}) for p in input_system]
    if k == n-1:
        Rk = PolynomialRing(QQ, variables[k:], n-k)
        input_system = [Rk(p) for p in input_system]
    else:
        input_system = [_remove_absent_variable(p,variables[:k]) for p in input_system]

    # Calling msolve to solve the above system.
    ms_out = sage.rings.polynomial.msolve._run_msolve((input_system[0].parent()).ideal(input_system), ['-P', '0', '-v', f"{msolve_verbose}", '-t', f"{threads}", '-p', f"{precision}"])
    sol = _format_output_msolve_intervals(ms_out)

    return substitution, sol

def _rough_eval(point,poly):
    r"""
    Internal Function

    Function computing the isolation interval that a polynomial takes on a box
    approximating a point. Although MPFI technically does it already, it is not
    precise enough for us in most examples

    INPUT:

        - ``point`` -- list of lists; list of isolation intervals,
        in msolve approximation format

        - ``poly`` -- polynomial with rational coefficients

    OUTPUT:

        - two rational numbers; lower and upper bounds for the value that `poly`
        can take on the approximation box of `point`.
    """

    varss = poly.parent().gens()
    if len(poly.variables()) == 0:
        return QQ(poly), QQ(poly)

    rg = range(len(varss))
    sign_list = [sign(coord[0]) for coord in point]
    modified_point = []
    for i in range(len(point)):
        if sign_list[i] == -1:
            modified_point.append([point[i][1], point[i][0]])
        else:
            modified_point.append([point[i][0], point[i][1]])
    min_out_poly = 0
    max_out_poly = 0
    if len(varss) == 1:
        poly = PolynomialRing(QQ, 1, varss)(poly)
        varss = PolynomialRing(QQ, 1, varss).gens()
    for coeff,monom in poly:
        signn = coeff*monom.subs({varss[i] : sign_list[i] for i in rg})
        if sign(signn) == -1:
            min_out_poly += coeff*monom.subs({varss[i] : modified_point[i][1] for i in rg})
            max_out_poly += coeff*monom.subs({varss[i] : modified_point[i][0] for i in rg})
        else:
            min_out_poly += coeff*monom.subs({varss[i] : modified_point[i][0] for i in rg})
            max_out_poly += coeff*monom.subs({varss[i] : modified_point[i][1] for i in rg})
    return min_out_poly, max_out_poly

def _matrix_box(n, point, matrix, substitution):
        r"""
        Internal Function

        Function computing the isolation box of a real point after
        transformation by a matrix.

        INPUT:

            - ``n`` -- integer; size of the matrix, number of coordinates

            - ``point`` -- list of lists; list of isolation intervals,
            in msolve approximation format

            - ``matrix`` -- matrix; (n x n) matrix

            - ``substitution`` -- list of k-1 polynomials; substitution
            expressions for the first k-1 variables

        OUTPUT:

            - list; list of isolation intervals for the last n-k+1 coordinates
            of matrix*point, in msolve approximation format.
        """

        first_coords_point = [_rough_eval(point, item) for item in substitution]
        extended_point = first_coords_point + point
        vertices = [list(item) for item in itertools.product(*extended_point)]
        changed_vertices = [matrix*vector(vertex) for vertex in vertices]
        new_box = [[min(item[i] for item in changed_vertices),max(item[i] for item in changed_vertices)] for i in range(n)]
        for coord in new_box[len(substitution):]:
            if sign(coord[0]) != sign(coord[1]):
                raise ValueError("Not precise enough to ensure coordinate sign after transformation. Consider increasing the precision.")
        return new_box[len(substitution):]

def _do_boxes_intersect(point_list):
    r"""
    Internal Function

    Checks whether any point approximation box in a list intersects another

    INPUT:

        - ``point_list`` -- list of lists of lists; list of lists of
        isolation intervals, in msolve approximation format

    OUTPUT:

       - boolean; ``True`` if any box intersects another, ``False`` otherwise.
    """

    if point_list == []:
        return False
    n = len(point_list[0])
    all_vertices = [[list(item) for item in itertools.product(*point)] for point in point_list]
    for i in range(len(point_list)):
        vertices = all_vertices[i]
        for vertex in vertices:
            for j in range(i+1, len(point_list)):
                counter = 0
                for k in range(n):
                    if point_list[j][k][0] <= vertex[k] <= point_list[j][k][1]:
                        counter += 1
                if counter == n:
                    return True
    return False

def _transverse_intersection(poly, vars, point, low_prec, threads, msolve_verbose, precision):
    r"""
    Internal Function

    Function computing points 'to the left' and 'to the right' of the critical
    point, by means of the transverse line and real root approximation.

    INPUT:

        - ``poly`` -- polynomial with rational coefficients

        - ``vars`` -- list; list of all variables of the parent ring of the
        polynomial

        - ``point`` -- list of lists; list of isolation intervals,
        in msolve approximation format

        - ``low_prec`` list of lists; ``point``, but with coordinates at a
        lower precision

        - ``threads`` -- integer; number of threads to be used by msolve in
        computation

        - ``msolve_verbose`` -- 0, 1 or 2; msolve parameter for explicit
        description of computations, ranging from 0 (no desc) to 2 (full desc)

        - ``precision`` -- integer; number of bits of precision used by msolve
        for real root approximation

    OUTPUT:

        - Two points (in msolve format), each being to the 'left' and the
        'right' of the critical point on the transverse line, and sufficiently
        close to be in the right connected component.
    """

    UnivarRing = PolynomialRing(QQ, "ttttt")
    MultivarRing = PolynomialRing(QQ, 1, UnivarRing.variable_name())
    approx = [RealIntervalField(prec=5*precision)(item).simplest_rational(False, False) for item in point]
    list_to_sub = [UnivarRing("ttttt") + approx[0]]
    if len(approx) != 1:
        list_to_sub += approx[1:]
    # Changing base ring to a multivariate one because _run_msolve does not work
    # on univariate parents.
    transverse_poly = MultivarRing(poly.subs({vars[i] : list_to_sub[i] for i in range(len(vars))}))

    ms_out = sage.rings.polynomial.msolve._run_msolve(MultivarRing.ideal(transverse_poly), ['-P','0','-v',f"{msolve_verbose}",'-t',f"{threads}",'-p',f"{precision}"])
    inter = _format_output_msolve_intervals(ms_out)

    inter_lambda_values = inter[1][1]

    if inter_lambda_values == []:
        raise ValueError("Coordinates not precise enough to compute a good intersection line. Consider increasing the precision")

    allowed_lambda_interval = [low_prec[0] - approx[0], low_prec[1] - approx[0]]
    endpoints = list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(inter_lambda_values))))
    sorted_endpoints = sorted(endpoints, key=lambda x: (abs(x), x))

    if (allowed_lambda_interval[0] <= sorted_endpoints[0] <= allowed_lambda_interval[1]) \
            and (allowed_lambda_interval[0] <= sorted_endpoints[1] <= allowed_lambda_interval[1]):
        if len(sorted_endpoints) == 2:
            lambd = ceil(abs(sorted_endpoints[1]))+1
        else:
            lambd = RealIntervalField(prec=5*precision)(abs(sorted_endpoints[1]), abs(sorted_endpoints[2])).simplest_rational(True,True)
    else:
        raise ValueError("Coordinates not precise enough to compute a good intersection line. Consider increasing the precision")

    if len(sorted_endpoints) != 2 and allowed_lambda_interval[0] <= sorted_endpoints[2] <= allowed_lambda_interval[1]:
        raise ValueError("Isolation box not precise enough to gurantee a single intersection point of the transverse line inside it. Consider increasing the precision")

    right_pt = [approx[0] + lambd]
    left_pt = [approx[0] - lambd]
    if len(approx) != 1:
        right_pt += approx[1:]
        left_pt += approx[1:]
    return left_pt, right_pt

def _smooth_points_per_component(poly, threads, msolve_verbose, precision, inequation, isempty, changevar=True):
    r"""
    Internal Function

    Main function computing the points per connected components.

    INPUT:

        - ``poly`` -- polynomial with rational coefficients

        - ``threads`` -- integer; number of threads to be used by msolve in
        computation

        - ``msolve_verbose`` -- 0, 1 or 2; msolve parameter for explicit
        description of computations, ranging from 0 (no desc) to 2 (full desc)

        - ``precision`` -- integer; number of bits of precision used by msolve
        for real root approximation

        - ``inequation`` -- boolean (default ``True``); computes points per
        connected components of {x in R^n : f(x) =/= 0} if set to ``True``, and
        of {x in R^n : f(x) > 0} if set to ``False``

        - ``isempty`` -- boolean (default ``False``); if set to ``True``,
        computation stops as soon as a point in the set is computed

        - ``changevar`` -- boolean (default ``True``); uses A = Identity and
        sigma = [1,...,1] if set to ``False``

    OUTPUT:

        - list of lists; each sublist corresponds to the coordinates of a
        rational point
    """

    # Setting up correct parent rings and variables
    R = poly.parent()
    variables = list(R.gens())
    n = len(variables)

    # Re-naming variables to have the increasing partial derivatives degree
    inv_permutation, poly = _derivative_order(poly)
    f = R(poly)

    # Generating the change of variable matrix A, with its partial inverses
    list_of_matrices = _grp_random_matrix(n,changevar)

    # Generating the specialisation point sigma
    if changevar:
        sigma = [ZZ.random_element(1,100,"uniform") for i in range(n-1)]
    else:
        sigma = [1 for i in range(n-1)]

    # Pre-computing partial derivatives of f to avoid unnecessary computations
    der_list = [f.derivative(variables[i]) for i in range(n)]

    # Initialising final solutions list and number of solutions
    Sols = []

    # Computing f^A
    fA = f.subs({variables[i] : (list_of_matrices[0]*vector(variables))[i] for i in range(n)})

    # Main for loop
    for k in range(n):
        # print(f"k = {k}\n")

        # Obtaining approximations to critical points
        substitution, crit = _critical_points(f, threads, msolve_verbose,2*precision,k, n, list_of_matrices, sigma, variables, der_list)

        Solsk = []
        #print(crit)

        # In case we have infinitely many of them
        if crit[0] > 0:
            return _smooth_points_per_component(poly, threads, msolve_verbose, precision, inequation, isempty, changevar)

        # In case we have finitely many of them, and at least one
        if crit[0] != -1 and len(crit) < 3 and crit[1][1] != []:

            # Substituting the first variables & fixing parent ring issues
            f_sub = f.subs({variables[i] : substitution[i] for i in range(k)})
            fA_sub = fA.subs({variables[i] : sigma[i] for i in range(k)})

            if k == n-1:
                f_sub = f_sub.univariate_polynomial()
                fA_sub = fA_sub.univariate_polynomial()
                substitution = [item.univariate_polynomial() for item in substitution]
            else:
                f_sub = _remove_absent_variable(f_sub, variables[:k])
                fA_sub = _remove_absent_variable(fA_sub, variables[:k])
                substitution = [_remove_absent_variable(item, variables[:k]) for item in substitution]
            variabless = f_sub.parent().gens()

            # Looping over each computed point to obtain A^-1 * point
            if changevar:
                A_inv_list = []

                for point in crit[1][1]:

                    # Checking whether the sign of each coordinate is known
                    for coord in point:
                        if sign(coord[0]) != sign(coord[1]):
                            raise ValueError("Msolve not precise enough to determine coordinate sign. Consider increasing the precision.")

                    # Computing approximation box of original (no A) critical point
                    A_inv_point = _matrix_box(n,point,list_of_matrices[1],substitution)

                    # Checking whether the sign of each coordinate is still known
                    for coord in A_inv_point:
                        if sign(coord[0]) != sign(coord[1]):
                            raise ValueError("Msolve not precise enough to determine coordinate sign. Consider increasing the precision.")

                    A_inv_list.append(A_inv_point)

                     # Checking whether the new approximation boxes intersect
                    if _do_boxes_intersect(A_inv_list):
                        raise ValueError("Msolve not precise enough to isolate critical points. Consider increasing the precision.")

            else:
                A_inv_list = crit[1][1]

            dfAdxk = fA_sub.derivative(variabless[0])
            for point in A_inv_list:
                # Checking whether we do not have an exact point. If we actually
                # do, we can skip verification steps.
                if [item[0] for item in point] != [item[1] for item in point]:

                    # Checking if df^A/dx_k = 0 in the msolve approximation box
                    dfA_interval = _rough_eval(point, dfAdxk)

                    # Checking whether 0 is in the interval
                    if (dfA_interval[0] == 0 or dfA_interval[1] == 0 or sign(dfA_interval[0]) != sign(dfA_interval[1])):
                        raise ValueError("Msolve not precise enough to guarantee non-zero derivative. Consider increasing the precision.")

                    # Computing rougher approximation of that point in the x_k coordinate only
                    coord_low_prec = [floor(2**(precision)*point[0][0])/2**(precision), ceil(2**(precision)*point[0][1])/2**(precision)]

                    # Checking whether f vanishes on the low and high x_k
                    low_f = fA_sub.subs({variabless[0] : coord_low_prec[0]})
                    high_f = fA_sub.subs({variabless[0] : coord_low_prec[1]})
                    if k == n-1:
                        low_inter = [QQ(low_f), QQ(low_f)]
                        high_inter = [QQ(high_f), QQ(high_f)]
                    elif k == n-2:
                        low_f = low_f.univariate_polynomial()
                        high_f = high_f.univariate_polynomial()
                        low_inter = _rough_eval(point[1:], low_f)
                        high_inter = _rough_eval(point[1:], high_f)
                    else:
                        low_f = _remove_absent_variable(low_f, variabless[0])
                        high_f = _remove_absent_variable(high_f, variabless[0])
                        low_inter = _rough_eval(point[1:], low_f)
                        high_inter = _rough_eval(point[1:], high_f)

                    # Checking whether 0 is in any interval
                    if (low_inter[0] == 0 or low_inter[1] == 0 or sign(low_inter[0]) != sign(low_inter[1])) \
                    or (high_inter[0] == 0 or high_inter[1] == 0 or sign(high_inter[0]) != sign(high_inter[1])):
                        raise ValueError("Approximation is not sufficiently precise. Consider increasing the precision.")
                else:
                    coord_low_prec = [point[0][0], point[0][1]]
                # Computing the points lying on the transverse line in the
                # corresponding connected components
                left, right = _transverse_intersection(fA_sub, variabless, point, coord_low_prec, threads, msolve_verbose, precision)

                left = list_of_matrices[0]*vector(sigma[:k] + left)
                right = list_of_matrices[0]*vector(sigma[:k] + right)

                # Sanity check for the whole procedure
                if left == 0 or right == 0 or sign(f.subs({variables[i]: left[i] for i in range(n)})) == sign(f.subs({variables[i]: right[i] for i in range(n)})):
                    raise ValueError("Sanity check failed, something went wrong.")

                # Reverting to original coordinates
                left = inv_permutation(list(left))
                right = inv_permutation(list(right))

                if inequation or f.subs({variables[i] : left[i] for i in range(n)}) > 0:
                    if isempty:
                        return [left]
                    Solsk.append(left)
                if inequation or f.subs({variables[i] : right[i] for i in range(n)}) > 0:
                    if isempty:
                        return [right]
                    Solsk.append(right)

        Sols += [Solsk]

    final_point = list(list_of_matrices[0]*vector(sigma + [0]))
    if inequation or f.subs({variables[i] : final_point[i] for i in range(n)}) > 0:
        if isempty:
            return [final_point]
        Sols += [[final_point]]
    Sols = [x for xs in Sols for x in xs]
    return Sols

def points_per_components_single_inequality(poly, threads, msolve_verbose, precision=128, inequation=True, isempty=False, changevar=True, proof=False):
    r"""
    Function computing points per connected component of a semi-algebraic set
    defined by a single polynomial inequation.

    INPUT:

        - ``poly`` -- polynomial with rational coefficients

        - ``threads`` -- integer; number of threads to be used by msolve in
        computation

        - ``msolve_verbose`` -- 0, 1 or 2; msolve parameter for explicit
        description of computations, ranging from 0 (no desc) to 2 (full desc)

        - ``precision`` -- integer; number of bits of precision used by msolve
        for real root approximation

        - ``inequation`` -- boolean (default ``True``); computes points per
        connected components of {x in R^n : f(x) =/= 0} if set to ``True``, and
        of {x in R^n : f(x) > 0} if set to ``False``

        - ``isempty`` -- boolean (default ``False``); if set to ``True``,
        computation stops as soon as a point in the set is computed

        - ``changevar`` -- boolean (default ``True``); uses A = Identity and
        sigma = [1,...,1] if set to ``False``

        - ``proof`` -- boolean (default ``False``); inner sagemath variable to
        recall that the algorithm is based on msolve, which relies on
        heuristics. Returns an appropriate error message if set to ``True``

    OUTPUT:

        - list of lists; each sublist corresponds to the coordinates of a
        rational point

    EXAMPLES::

        sage: from sage.rings.polynomial.msolve import points_per_components_single_inequality
        sage: R.<x,y> = QQ[]
        sage: f = 2*x^2 - 7*x*y + 5*y^2 - 3*x + y - 2
        sage: points_per_components_single_inequality(f, 1, 0, 32, True, False, False) # optional - msolve
        [[-81327242789970644619/18230669955817908734, -151863759432751750973/52864310410056725586], [-63096572834152735885/18230669955817908734, -151863759432751750973/52864310410056725586], [-31255250995381558863/18941441967027853418, -2031568568823642965792129/2244655580804164026831606], [-12313809028353705445/18941441967027853418, -2031568568823642965792129/2244655580804164026831606], [1, -90794324783669909401/65802721238427014701], [1, 40811117693184120001/65802721238427014701], [1, 9172089397301669399/15819514147941225301], [1, 40811117693184120001/15819514147941225301], [1, 0]]

    ::

        sage: from sage.rings.polynomial.msolve import points_per_components_single_inequality
        sage: R.<x,y> = QQ[]
        sage: f = x^2 + y^2
        sage: points_per_components_single_inequality(f, 1, 0, 128, True, False, False) # optional - msolve
        Traceback (most recent call last):
        ...
        ValueError: Input polynomial does not define a smooth hypersurface, this case is not yet implemented.
    """

    if sage.structure.proof.proof.get_flag(proof, "polynomial"):
        raise ValueError("msolve relies on heuristics; please use proof=False.")

    if _is_smooth(poly, threads, msolve_verbose):
        return _smooth_points_per_component(poly, threads, msolve_verbose, precision, inequation, isempty, changevar)
    else:
        raise ValueError("Input polynomial does not define a smooth hypersurface, this case is not yet implemented.")
