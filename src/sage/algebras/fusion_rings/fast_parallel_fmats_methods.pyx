"""
Fast F-Matrix Methods
"""
# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython
from sage.algebras.fusion_rings.poly_tup_engine cimport (
    compute_known_powers,
    get_variables_degrees, variables,
    poly_to_tup, _tup_to_poly,
    subs, subs_squares, reduce_poly_dict, resize,
    _flatten_coeffs, _unflatten_coeffs,
    has_appropriate_linear_term,
    resize
)
from sage.algebras.fusion_rings.shm_managers cimport KSHandler, FvarsHandler
from sage.rings.number_field.number_field_element cimport NumberFieldElement_absolute
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular, MPolynomialRing_libsingular
from sage.rings.polynomial.polydict cimport ETuple

from ctypes import cast, py_object
from itertools import product
from sage.rings.ideal import Ideal
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

##########################
#   Fast class methods   #
##########################

cpdef _solve_for_linear_terms(factory, list eqns=None):
    r"""
    Solve for a linear term occurring in a two-term equation, and for
    variables appearing in univariate single-term equations.

    EXAMPLES::

        sage: f = FusionRing("D3", 1).get_fmatrix(inject_variables=True, new=True)
        creating variables fx1..fx27
        Defining fx0, ..., fx26
        sage: f._reset_solver_state()
        sage: f.ideal_basis = [fx0**3, fx0 + fx3**4, fx2**2 - fx3, fx2 - fx3**2, fx4 - fx2]
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: f.ideal_basis = [poly_to_tup(p) for p in f.ideal_basis]
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_tup_sortkey
        sage: f.ideal_basis.sort(key=poly_tup_sortkey)
        sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
        sage: n = f._poly_ring.ngens()
        sage: f._fvars = FvarsHandler(n, f._field, f._idx_to_sextuple, init_data=f._fvars)
        sage: from sage.algebras.fusion_rings.fast_parallel_fmats_methods import _solve_for_linear_terms
        sage: _solve_for_linear_terms(f)
        True
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[0]])
        0
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[2]])
        fx4
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[3]])
        fx3
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[4]])
        fx4
    """
    if eqns is None:
        eqns = factory.ideal_basis

    cdef bint linear_terms_exist = False
    cdef ETuple exp, rhs_exp
    cdef int max_var
    cdef tuple eq_tup
    cdef FvarsHandler fvars = factory._fvars
    cdef NumberFieldElement_absolute coeff, other
    cdef tuple rhs_coeff
    for eq_tup in eqns:
        if len(eq_tup) == 1:
            vars = variables(eq_tup)
            if len(vars) == 1 and not factory._solved[vars[0]]:
                fvars[factory._idx_to_sextuple[vars[0]]] = tuple()
                factory._solved[vars[0]] = True
                linear_terms_exist = True

                # TESTS:
                # s = factory._idx_to_sextuple[vars[0]]
                # factory.test_fvars[s] = tuple()
                # assert factory.test_fvars[s] == fvars[s], "OG value {}, Shared: {}".format(fvars[s], factory.test_fvars[s])
        if len(eq_tup) == 2:
            idx = has_appropriate_linear_term(eq_tup)
            if idx < 0:
                continue
            # The chosen term is guaranteed to be univariate in the largest variable
            exp = eq_tup[idx][0]
            max_var = exp._data[0]
            if not factory._solved[max_var]:
                rhs_exp = eq_tup[(idx+1) % 2][0]
                coeff = factory._field(list(eq_tup[(idx+1) % 2][1]))
                other = factory._field(list(eq_tup[idx][1]))
                rhs_coeff = tuple((-coeff / other)._coefficients())
                fvars[factory._idx_to_sextuple[max_var]] = ((rhs_exp, rhs_coeff), )
                factory._solved[max_var] = True
                linear_terms_exist = True

                # TESTS:
                # s = factory._idx_to_sextuple[max_var]
                # factory.test_fvars[s] = ((rhs_exp, rhs_coeff), )
                # assert _unflatten_coeffs(factory._field, factory.test_fvars[s]) == fvars[s], "OG value {}, Shared: {}".format(factory.test_fvars[s], fvars[s])
    return linear_terms_exist

cpdef _backward_subs(factory, bint flatten=True):
    r"""
    Perform backward substitution on ``self.ideal_basis``, traversing
    variables in reverse lexicographical order.

    EXAMPLES::

        sage: f = FusionRing("D3", 1).get_fmatrix(inject_variables=True, new=True)
        creating variables fx1..fx27
        Defining fx0, ..., fx26
        sage: f._reset_solver_state()
        sage: f.ideal_basis = [fx0**3, fx0 + fx3**4, fx2**2 - fx3, fx2 - fx3**2, fx4 - fx2]
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: f.ideal_basis = [poly_to_tup(p) for p in f.ideal_basis]
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_tup_sortkey
        sage: f.ideal_basis.sort(key=poly_tup_sortkey)
        sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
        sage: n = f._poly_ring.ngens()
        sage: f._fvars = FvarsHandler(n, f._field, f._idx_to_sextuple, init_data=f._fvars)
        sage: from sage.algebras.fusion_rings.fast_parallel_fmats_methods import _solve_for_linear_terms
        sage: _solve_for_linear_terms(f)
        True
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[0]])
        0
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[2]])
        fx4
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[3]])
        fx3
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[4]])
        fx4
        sage: f.ideal_basis.append(poly_to_tup(fx4-1))
        sage: f.ideal_basis.sort(key=poly_tup_sortkey)
        sage: _solve_for_linear_terms(f)
        True
        sage: from sage.algebras.fusion_rings.fast_parallel_fmats_methods import _backward_subs
        sage: _backward_subs(f)
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[0]])
        0
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[2]])
        1
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[3]])
        fx3
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[4]])
        1
    """
    one = factory._field.one()
    fvars = factory._fvars
    solved = factory._solved
    cdef KSHandler _ks = factory._ks
    cdef dict idx_to_sextuple = factory._idx_to_sextuple
    cdef int nvars = len(idx_to_sextuple)
    for i in range(nvars-1, -1, -1):
        sextuple = idx_to_sextuple[i]
        rhs = fvars[sextuple]
        d = {var_idx: fvars[idx_to_sextuple[var_idx]]
             for var_idx in variables(rhs) if solved[var_idx]}
        if d:
            kp = compute_known_powers(get_variables_degrees([rhs], nvars), d, one)
            res = tuple(subs_squares(subs(rhs, kp, one), _ks).items())
            if flatten:
                res = _flatten_coeffs(res)
            fvars[sextuple] = res


cdef _fmat(fvars, _Nk_ij, id_anyon, a, b, c, d, x, y):
    """
    Cython version of fmat class method. Using cdef for fastest dispatch
    """
    if _Nk_ij(a, b, x) == 0 or _Nk_ij(x, c, d) == 0 or _Nk_ij(b, c, y) == 0 or _Nk_ij(a, y, d) == 0:
        return 0
    # Some known F-symbols
    if a == id_anyon:
        return int(x == b and y == d)
    if b == id_anyon:
        return int(x == a and y == c)
    if c == id_anyon:
        return int(x == d and y == b)
    return fvars[a, b, c, d, x, y]

######################################
#   Fast fusion coefficients cache   #
######################################

# from sage.misc.cachefunc import cached_function
# cdef dict _Nk_ij = dict()

# cpdef _Nk_ij(factory, proc):
# cdef int coeff
# for a, b, c in product(factory._FR.basis(), repeat=3):
# try:
# coeff = (a*b).monomial_coefficients(copy=False)[c.weight()]
# except:
# coeff = 0
# _Nk_ij[a, b, c] = coeff

# cpdef int _Nk_ij(a, b, c):
# try:
# return (a*b).monomial_coefficients(copy=False)[c.weight()]
# except KeyError:
# return 0
#
# _Nk_ij = cached_function(_Nk_ij, name='_Nk_ij')

###############
#   Mappers   #
###############

cdef req_cy(tuple basis, r_matrix, dict fvars, Nk_ij, id_anyon, tuple sextuple):
    """
    Given an FMatrix factory and a sextuple, return a hexagon equation
    as a polynomial object.
    """
    a, b, c, d, e, g = sextuple
    # To add typing we need to ensure all fmats.fmat are of the same type?
    # Return fmats._poly_ring.zero() and fmats._poly_ring.one() instead of 0 and 1?
    lhs = r_matrix(a, c, e, base_coercion=False) * _fmat(fvars, Nk_ij, id_anyon, a, c, b, d, e, g) * r_matrix(b, c, g, base_coercion=False)
    rhs = 0
    for f in basis:
        rhs += _fmat(fvars, Nk_ij, id_anyon, c, a, b, d, e, f) * r_matrix(f, c, d, base_coercion=False) * _fmat(fvars, Nk_ij, id_anyon, a, b, c, d, f, g)
    return lhs-rhs


@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef get_reduced_hexagons(factory, tuple mp_params):
    """
    Set up and reduce the hexagon equations corresponding to this worker.
    """
    # Set up multiprocessing parameters
    cdef list worker_results = list()
    cdef int child_id, n_proc
    cdef unsigned long i
    child_id, n_proc, output = mp_params
    cdef tuple sextuple, red

    # Pre-compute common parameters for speed
    cdef tuple basis = tuple(factory._FR.basis())
    cdef dict fvars
    cdef bint must_zip_up = False
    if not output:
        fvars = {s: factory._poly_ring.gen(i) for i, s in factory._idx_to_sextuple.items()}
    else:
        # Handle both cyclotomic and orthogonal solution method
        for k, v in factory._fvars.items():
            must_zip_up = isinstance(v, tuple)
            break
        if must_zip_up:
            fvars = {k: _tup_to_poly(v, parent=factory._poly_ring) for k, v in factory._fvars.items()}
        else:
            fvars = factory._fvars
    r_matrix = factory._FR.r_matrix
    _Nk_ij = factory._FR.Nk_ij
    id_anyon = factory._FR.one()
    _field = factory._field
    cdef NumberFieldElement_absolute one = _field.one()
    cdef ETuple _nnz = factory._nnz
    _ks = factory._ks

    # Computation loop
    it = product(basis, repeat=6)
    for i in range(len(basis)**6):
        sextuple = next(it)
        if i % n_proc == child_id:
            he = req_cy(basis, r_matrix, fvars, _Nk_ij, id_anyon, sextuple)
            if he:
                red = reduce_poly_dict(he.monomial_coefficients(), _nnz, _ks, one)

                # Avoid pickling cyclotomic coefficients
                red = _flatten_coeffs(red)

                worker_results.append(red)

    return collect_eqns(worker_results)

cdef MPolynomial_libsingular feq_cy(tuple basis, fvars, Nk_ij, id_anyon, zero, tuple nonuple, bint prune=False):
    r"""
    Given an FMatrix factory and a nonuple, return a pentagon equation
    as a polynomial object.
    """
    a, b, c, d, e, f, g, k, l = nonuple
    lhs = _fmat(fvars, Nk_ij, id_anyon, f, c, d, e, g, l) * _fmat(fvars, Nk_ij, id_anyon, a, b, l, e, f, k)
    if lhs == 0 and prune:
        # it is believed that if lhs=0, the equation carries no new information
        return zero
    rhs = zero
    for h in basis:
        rhs += _fmat(fvars, Nk_ij, id_anyon, a, b, c, g, f, h)*_fmat(fvars, Nk_ij, id_anyon, a, h, d, e, g, k)*_fmat(fvars, Nk_ij, id_anyon, b, c, d, k, h, l)
    return lhs - rhs


@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef get_reduced_pentagons(factory, tuple mp_params):
    r"""
    Set up and reduce the pentagon equations corresponding to this worker.
    """
    # Set up multiprocessing parameters
    cdef list worker_results = list()
    cdef int child_id, n_proc
    child_id, n_proc, _ = mp_params
    cdef unsigned long i
    cdef tuple nonuple, red
    cdef MPolynomial_libsingular pe

    # Pre-compute common parameters for speed
    cdef tuple basis = tuple(factory._FR.basis())
    # Handle both cyclotomic and orthogonal solution method
    cdef bint must_zip_up
    for k, v in factory._fvars.items():
        must_zip_up = isinstance(v, tuple)
        break
    cdef dict fvars
    if must_zip_up:
        fvars = {k: _tup_to_poly(v, parent=factory._poly_ring) for k, v in factory._fvars.items()}
    else:
        fvars = factory._fvars
    _Nk_ij = factory._FR.Nk_ij
    id_anyon = factory._FR.one()
    _field = factory._field
    cdef NumberFieldElement_absolute one = _field.one()
    cdef MPolynomial_libsingular zero = factory._poly_ring.zero()
    cdef KSHandler _ks = factory._ks
    factory._nnz = factory._get_known_nonz()
    cdef ETuple _nnz = factory._nnz

    # Computation loop
    it = product(basis, repeat=9)
    for i in range(len(basis)**9):
        nonuple = next(it)
        if i % n_proc == child_id:
            pe = feq_cy(basis, fvars, _Nk_ij, id_anyon, zero, nonuple, prune=True)
            if pe:
                red = reduce_poly_dict(pe.monomial_coefficients(), _nnz, _ks, one)

                # Avoid pickling cyclotomic coefficients
                red = _flatten_coeffs(red)

                worker_results.append(red)
    return collect_eqns(worker_results)

cdef list update_reduce(factory, list eqns):
    r"""
    Substitute known values, known squares, and reduce.
    """
    cdef list res = list()
    cdef tuple eq_tup, red, unflat
    cdef dict eq_dict

    # Pre-compute common parameters for speed
    _field = factory._field
    one = _field.one()
    cdef KSHandler _ks = factory._ks
    # Update reduction params
    factory._nnz = factory._get_known_nonz()
    factory._kp = compute_known_powers(factory._var_degs, factory._get_known_vals(), factory._field.one())
    cdef dict _kp = factory._kp
    cdef ETuple _nnz = factory._nnz

    for i in range(len(eqns)):
        eq_tup = eqns[i]
        # Construct cyclotomic field elts from list repn
        unflat = _unflatten_coeffs(_field, eq_tup)

        eq_dict = subs(unflat, _kp, one)
        red = reduce_poly_dict(eq_dict, _nnz, _ks, one)

        # Avoid pickling cyclotomic coefficients
        red = _flatten_coeffs(red)

        res.append(red)
    return collect_eqns(res)

cdef list compute_gb(factory, tuple args):
    r"""
    Compute the reduced Groebner basis for given equations iterable.
    """
    cdef list res = list()
    cdef list eqns, sorted_vars
    eqns, term_order = args
    # Define smaller poly ring in component vars
    sorted_vars = []
    cdef tuple eq_tup
    cdef int fx
    for eq_tup in eqns:
        for fx in variables(eq_tup):
            sorted_vars.append(fx)
    sorted_vars = sorted(set(sorted_vars))
    cdef MPolynomialRing_libsingular R = PolynomialRing(factory._FR.field(), len(sorted_vars), 'a', order=term_order)

    # Zip tuples into R and compute Groebner basis
    cdef idx_map = {old: new for new, old in enumerate(sorted_vars)}
    nvars = len(sorted_vars)
    F = factory.field()
    cdef list polys = list()
    for eq_tup in eqns:
        eq_tup = _unflatten_coeffs(F, eq_tup)
        polys.append(_tup_to_poly(resize(eq_tup, idx_map, nvars), parent=R))
    gb = Ideal(sorted(polys)).groebner_basis(algorithm='libsingular:slimgb')

    # Change back to fmats poly ring and append to temp_eqns
    cdef dict inv_idx_map = {v: k for k, v in idx_map.items()}
    cdef tuple t
    nvars = factory._poly_ring.ngens()
    for p in gb:
        t = resize(poly_to_tup(p), inv_idx_map, nvars)

        # Avoid pickling cyclotomic coefficients
        t = _flatten_coeffs(t)

        res.append(t)
    return collect_eqns(res)

################
#   Reducers   #
################

cdef inline list collect_eqns(list eqns):
    r"""
    Helper function for returning processed results back to parent process.

    Trivial reducer: simply collects objects with the same key in the worker.
    This method is only useful when called after :meth:`executor`, whose
    function argument appends output to the ``worker_results`` list.
    """
    # Discard the zero polynomial
    reduced = set(eqns) - set([tuple()])
    return list(reduced)

##############################
#   Parallel code executor   #
##############################

# Hard-coded module __dict__-style attribute with visible cdef methods
cdef dict mappers = {
    "get_reduced_hexagons": get_reduced_hexagons,
    "get_reduced_pentagons": get_reduced_pentagons,
    "update_reduce": update_reduce,
    "compute_gb": compute_gb,
    "pent_verify": pent_verify
    }

cpdef executor(tuple params):
    r"""
    Execute a function defined in this module
    (``sage.algebras.fusion_rings.fast_parallel_fmats_methods``) in a worker
    process, and supply the factory parameter by constructing a reference
    to the ``FMatrix`` object in the worker's memory address space from
    its ``id``.

    INPUT:

    - ``params`` -- tuple ``((fn_name, fmats_id), fn_args)`` where
      ``fn_name`` is the name of the function to be executed, ``fmats_id``
      is the ``id`` of the :class:`FMatrix` object, and ``fn_args`` is a
      tuple containing all arguments to be passed to the function ``fn_name``.

    .. NOTE::

        When the parent process is forked, each worker gets a copy of
        every  global variable. The virtual memory address of object `X` in
        the parent process equals the *virtual* memory address of the copy
        of object `X` in each worker, so we may construct references to
        forked copies of `X` using an ``id`` obtained in the parent process.

    TESTS::

        sage: from sage.algebras.fusion_rings.fast_parallel_fmats_methods import executor
        sage: fmats = FusionRing("A1", 3).get_fmatrix()
        sage: fmats._reset_solver_state()
        sage: params = (('get_reduced_hexagons', id(fmats)), (0, 1, True))
        sage: len(executor(params)) == 63
        True
        sage: fmats = FusionRing("E6", 1).get_fmatrix()
        sage: fmats._reset_solver_state()
        sage: params = (('get_reduced_hexagons', id(fmats)), (0, 1, False))
        sage: len(executor(params)) == 6
        True
    """
    (fn_name, fmats_id), args = params
    # Construct a reference to global FMatrix object in this worker's memory
    fmats_obj = cast(fmats_id, py_object).value
    # Bind module method to FMatrix object in worker process, and call the method
    return mappers[fn_name](fmats_obj, args)

####################
#   Verification   #
####################

cdef feq_verif(factory, worker_results, fvars, Nk_ij, id_anyon, tuple nonuple, float tol=5e-8):
    r"""
    Check the pentagon equation corresponding to the given nonuple.
    """
    a, b, c, d, e, f, g, k, l = nonuple
    cdef float diff, lhs, rhs

    lhs = _fmat(fvars, Nk_ij, id_anyon, f, c, d, e, g, l)*_fmat(fvars, Nk_ij, id_anyon, a, b, l, e, f, k)
    rhs = 0.0
    for h in factory._FR.basis():
        rhs += _fmat(fvars, Nk_ij, id_anyon, a, b, c, g, f, h)*_fmat(fvars, Nk_ij, id_anyon, a, h, d, e, g, k)*_fmat(fvars, Nk_ij, id_anyon, b, c, d, k, h, l)
    diff = lhs - rhs
    if diff > tol or diff < -tol:
        worker_results.append(diff)


@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef pent_verify(factory, tuple mp_params):
    r"""
    Generate all the pentagon equations assigned to this process,
    and reduce them.
    """
    child_id, n_proc, verbose = mp_params
    cdef tuple nonuple
    cdef unsigned long long i
    cdef list worker_results = list()

    # Pre-compute common parameters for speed
    Nk_ij = factory._FR.Nk_ij
    cdef dict fvars = factory._fvars
    id_anyon = factory._FR.one()
    for i, nonuple in enumerate(product(factory._FR.basis(), repeat=9)):
        if i % n_proc == child_id:
            feq_verif(factory, worker_results, fvars, Nk_ij, id_anyon, nonuple)
        if i % 50000000 == 0 and i and verbose:
            print("{:5d}m equations checked... {} potential misses so far...".format(i // 1000000, len(worker_results)))
