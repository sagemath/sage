"""
Fast Fusion Ring Methods for Computing Braid Group Representations
"""
# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from ctypes import cast, py_object
cimport cython
from sage.algebras.fusion_rings.fast_parallel_fmats_methods cimport _fmat

from sage.rings.qqbar import QQbar

###############
#   Mappers   #
###############

cdef mid_sig_ij(fusion_ring, row, col, a, b):
    r"""
    Compute the (xi, yi), (xj, yj) entry of generator braiding the middle two
    strands in the tree b -> xi # yi -> (a # a) # (a # a), which results in
    a sum over j of trees b -> xj # yj -> (a # a) # (a # a)

    .. WARNING::

        This method assumes F-matrices are orthogonal.
    """
    # Pre-compute common parameters for efficiency
    _fvars = fusion_ring.get_fmatrix()._fvars
    _Nk_ij = fusion_ring.Nk_ij
    one = fusion_ring.one()

    xi, yi = row
    xj, yj = col
    entry = 0
    cdef list basis = list(fusion_ring.basis())
    for c in basis:
        for d in basis:
            # #Warning: We assume F-matrices are orthogonal!!! (using transpose for inverse)
            f1 = _fmat(_fvars, _Nk_ij, one, a, a, yi, b, xi, c)
            f2 = _fmat(_fvars, _Nk_ij, one, a, a, a, c, d, yi)
            f3 = _fmat(_fvars, _Nk_ij, one, a, a, a, c, d, yj)
            f4 = _fmat(_fvars, _Nk_ij, one, a, a, yj, b, xj, c)
            r = fusion_ring.r_matrix(a, a, d)
            entry += f1 * f2 * r * f3 * f4
    return entry

cdef odd_one_out_ij(fusion_ring, xi, xj, a, b):
    r"""
    Compute the `xi`, `xj` entry of the braid generator on the two right-most
    strands, corresponding to the tree b -> (xi # a) -> (a # a) # a, which
    results in a sum over j of trees b -> xj -> (a # a) # (a # a)

    .. WARNING::

        This method assumes F-matrices are orthogonal.
    """
    # Pre-compute common parameters for efficiency
    _fvars = fusion_ring.get_fmatrix()._fvars
    _Nk_ij = fusion_ring.Nk_ij
    one = fusion_ring.one()

    entry = 0
    for c in fusion_ring.basis():
        # #Warning: We assume F-matrices are orthogonal!!! (using transpose for inverse)
        f1 = _fmat(_fvars, _Nk_ij, one, a, a, a, b, xi, c)
        f2 = _fmat(_fvars, _Nk_ij, one, a, a, a, b, xj, c)
        r = fusion_ring.r_matrix(a, a, c)
        entry += f1 * r * f2
    return entry

# Cache methods (manually for cdef methods)
cdef odd_one_out_ij_cache = dict()
cdef mid_sig_ij_cache = dict()

cdef cached_mid_sig_ij(fusion_ring, row, col, a, b):
    r"""
    Cached version of :meth:`mid_sig_ij`.
    """
    if (row, col, a, b) in mid_sig_ij_cache:
        return mid_sig_ij_cache[row, col, a, b]
    entry = mid_sig_ij(fusion_ring, row, col, a, b)
    mid_sig_ij_cache[row, col, a, b] = entry
    return entry

cdef cached_odd_one_out_ij(fusion_ring, xi, xj, a, b):
    r"""
    Cached version of :meth:`odd_one_out_ij`.
    """
    if (xi, xj, a, b) in odd_one_out_ij_cache:
        return odd_one_out_ij_cache[xi, xj, a, b]
    entry = odd_one_out_ij(fusion_ring, xi, xj, a, b)
    odd_one_out_ij_cache[xi, xj, a, b] = entry
    return entry


@cython.nonecheck(False)
@cython.cdivision(True)
cdef sig_2k(fusion_ring, tuple args):
    r"""
    Compute entries of the `2k`-th braid generator
    """
    # Pre-compute common parameters for efficiency
    _fvars = fusion_ring.get_fmatrix()._fvars
    _Nk_ij = fusion_ring.Nk_ij
    one = fusion_ring.one()

    cdef int child_id, n_proc
    child_id, n_proc, fn_args = args
    k, a, b, n_strands = fn_args
    cdef int ctr = -1
    cdef list worker_results = list()
    # Get computational basis
    cdef list comp_basis = fusion_ring.get_computational_basis(a, b, n_strands)
    cdef dict basis_dict = {elt: i for i, elt in enumerate(comp_basis)}
    cdef int dim = len(comp_basis)
    cdef set coords = set()
    cdef int i
    # Avoid pickling cyclotomic field element objects
    cdef bint must_flatten_coeff = fusion_ring.fvars_field() != QQbar
    cdef list basis = list(fusion_ring.basis())
    for i in range(dim):
        for f in basis:
            for e in basis:
                for q in basis:
                    # Distribute work amongst processes
                    ctr += 1
                    if ctr % n_proc != child_id:
                        continue

                    # Compute appropriate possible nonzero row index
                    nnz_pos = list(comp_basis[i])
                    nnz_pos[k-1] = f
                    nnz_pos[k] = e
                    # Handle the special case k = 1
                    if k > 1:
                        nnz_pos[n_strands//2+k-2] = q
                    nnz_pos = tuple(nnz_pos)

                    # Skip repeated entries when k = 1
                    if nnz_pos in comp_basis and (basis_dict[nnz_pos], i) not in coords:
                        m, l = comp_basis[i][:n_strands//2], comp_basis[i][n_strands//2:]
                        # A few special cases
                        top_left = m[0]
                        if k >= 3:
                            top_left = l[k-3]
                        root = b
                        if k - 1 < len(l):
                            root = l[k-1]

                        # Handle the special case k = 1
                        if k == 1:
                            entry = cached_mid_sig_ij(fusion_ring, m[:2], (f, e), a, root)

                            # Avoid pickling cyclotomic field element objects
                            if must_flatten_coeff:
                                entry = entry.list()

                            worker_results.append(((basis_dict[nnz_pos], i), entry))
                            coords.add((basis_dict[nnz_pos], i))
                            continue

                        entry = 0
                        for p in fusion_ring.basis():
                            f1 = _fmat(_fvars, _Nk_ij, one, top_left, m[k-1], m[k], root, l[k-2], p)
                            f2 = _fmat(_fvars, _Nk_ij, one, top_left, f, e, root, q, p)
                            entry += f1 * cached_mid_sig_ij(fusion_ring, (m[k-1], m[k]), (f, e), a, p) * f2

                        # Avoid pickling cyclotomic field element objects
                        if must_flatten_coeff:
                            entry = entry.list()

                        worker_results.append(((basis_dict[nnz_pos], i), entry))
    return worker_results


@cython.nonecheck(False)
@cython.cdivision(True)
cdef odd_one_out(fusion_ring, tuple args):
    r"""
    Compute entries of the rightmost braid generator, in case we have an
    odd number of strands.
    """
    # Pre-compute common parameters for efficiency
    _fvars = fusion_ring.get_fmatrix()._fvars
    _Nk_ij = fusion_ring.Nk_ij
    one = fusion_ring.one()

    cdef list worker_results = []
    cdef list nnz_pos_temp
    cdef tuple nnz_pos
    cdef int child_id, n_proc, i
    child_id, n_proc, fn_args = args
    a, b, n_strands = fn_args
    cdef int ctr = -1
    # Get computational basis
    cdef list comp_basis = fusion_ring.get_computational_basis(a, b, n_strands)
    cdef dict basis_dict = {elt: i for i, elt in enumerate(comp_basis)}
    cdef int dim = len(comp_basis)

    # Avoid pickling cyclotomic field element objects
    cdef bint must_flatten_coeff = fusion_ring.fvars_field() != QQbar

    cdef list basis = list(fusion_ring.basis())
    for i in range(dim):
        for f in basis:
            for q in basis:
                # Distribute work amongst processes
                ctr += 1
                if ctr % n_proc != child_id:
                    continue

                # Compute appropriate possible nonzero row index
                nnz_pos_temp = list(comp_basis[i])
                nnz_pos_temp[n_strands//2-1] = f
                # Handle small special case
                if n_strands > 3:
                    nnz_pos_temp[-1] = q
                nnz_pos = tuple(nnz_pos_temp)

                if nnz_pos in comp_basis:
                    m, l = comp_basis[i][:n_strands//2], comp_basis[i][n_strands//2:]

                    # Handle a couple of small special cases
                    if n_strands == 3:
                        entry = cached_odd_one_out_ij(fusion_ring, m[-1], f, a, b)

                        # Avoid pickling cyclotomic field element objects
                        if must_flatten_coeff:
                            entry = entry.list()

                        worker_results.append(((basis_dict[nnz_pos], i), entry))
                        continue
                    top_left = m[0]
                    if n_strands > 5:
                        top_left = l[-2]
                    root = b

                    # Compute relevant entry
                    entry = 0
                    for p in basis:
                        f1 = _fmat(_fvars, _Nk_ij, one, top_left, m[-1], a, root, l[-1], p)
                        f2 = _fmat(_fvars, _Nk_ij, one, top_left, f, a, root, q, p)
                        entry += f1 * cached_odd_one_out_ij(fusion_ring, m[-1], f, a, p) * f2

                    # Avoid pickling cyclotomic field element objects
                    if must_flatten_coeff:
                        entry = entry.list()

                    worker_results.append(((basis_dict[nnz_pos], i), entry))
    return worker_results

##############################
#   Parallel code executor   #
##############################

# Hard-coded module __dict__-style attribute with visible cdef methods
cdef dict mappers = {
    "sig_2k": sig_2k,
    "odd_one_out": odd_one_out
}

cpdef executor(tuple params):
    r"""
    Execute a function registered in this module's ``mappers``
    in a worker process, and supply the ``FusionRing`` parameter by
    constructing a reference to the FMatrix object in the worker's memory
    address space from its ``id``.

    .. NOTE::

        When the parent process is forked, each worker gets a copy of
        every  global variable. The virtual memory address of object `X` in
        the parent process equals the *virtual* memory address of the copy of
        object `X` in each worker, so we may construct references to forked
        copies of `X`.

    TESTS::

        sage: from sage.algebras.fusion_rings.fast_parallel_fusion_ring_braid_repn import executor
        sage: FR = FusionRing("A1", 4)
        sage: FR.fusion_labels(['idd', 'one', 'two', 'three', 'four'], inject_variables=True)
        sage: FR.get_fmatrix().find_orthogonal_solution(verbose=False)    # long time
        sage: params = (('sig_2k', id(FR)), (0, 1, (1, one, one, 5)))    # long time
        sage: len(executor(params)) == 13                         # long time
        True
        sage: from sage.algebras.fusion_rings.fast_parallel_fusion_ring_braid_repn import executor
        sage: FR = FusionRing("A1", 2)
        sage: FR.fusion_labels("a", inject_variables=True)
        sage: FR.get_fmatrix().find_orthogonal_solution(verbose=False)
        sage: params = (('odd_one_out', id(FR)), (0, 1, (a2, a2, 5)))
        sage: len(executor(params)) == 1
        True
    """
    (fn_name, fr_id), args = params
    # Construct a reference to global FMatrix object in this worker's memory
    fusion_ring_obj = cast(fr_id, py_object).value
    # Bind module method to FMatrix object in worker process, and call the method
    return mappers[fn_name](fusion_ring_obj, args)

######################################
#   Pickling circumvention helpers   #
######################################

cpdef _unflatten_entries(fusion_ring, list entries):
    r"""
    Restore cyclotomic coefficient object from its tuple of rational
    coefficients representation.

    Used to circumvent pickling issue introduced by PARI settings
    in :issue:`30537`.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.fast_parallel_fusion_ring_braid_repn import _unflatten_entries
        sage: fr = FusionRing("B2", 2)
        sage: F = fr.field()
        sage: coeff = [F.random_element() for i in range(2)]
        sage: entries = [((0, 0), coeff[0].list()), ((0, 1), coeff[1].list())]
        sage: _unflatten_entries(fr, entries)
        sage: all(cyc_elt_obj == c for (coord, cyc_elt_obj), c in zip(entries, coeff))
        True
    """
    F = fusion_ring.fvars_field()
    if F != QQbar:
        for i, (coord, entry) in enumerate(entries):
            entries[i] = (coord, F(entry))
