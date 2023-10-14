# distutils: libraries = flint
# distutils: depends = flint/fmpz_lll.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpz_lll_context_init_default(fmpz_lll_t fl)
    # Sets ``fl->delta``, ``fl->eta``, ``fl->rt`` and ``fl->gt`` to
    # their default values, 0.99, 0.51, `Z\_BASIS` and `APPROX` respectively.

    void fmpz_lll_context_init(fmpz_lll_t fl, double delta, double eta, rep_type rt, gram_type gt)
    # Sets ``fl->delta``, ``fl->eta``, ``fl->rt`` and ``fl->gt`` to
    # ``delta``, ``eta``, ``rt`` and ``gt`` (given as input)
    # respectively. ``delta`` and ``eta`` are the L^2 parameters.
    # ``delta`` and ``eta`` must lie in the intervals `(0.25, 1)` and
    # `(0.5, \sqrt{\mathtt{delta}})` respectively. The representation type is input
    # using ``rt`` and can have the values `Z\_BASIS` for a lattice basis and
    # `GRAM` for a Gram matrix. The Gram type to be used during computation can
    # be specified using ``gt`` which can assume the values `APPROX` and
    # `EXACT`. Note that ``gt`` has meaning only when ``rt`` is `Z\_BASIS`.

    void fmpz_lll_randtest(fmpz_lll_t fl, flint_rand_t state)
    # Sets ``fl->delta`` and ``fl->eta`` to random values in the interval
    # `(0.25, 1)` and `(0.5, \sqrt{\mathtt{delta}})` respectively. ``fl->rt`` is
    # set to `GRAM` or `Z\_BASIS` and ``fl->gt`` is set to `APPROX` or `EXACT`
    # in a pseudo random way.

    double fmpz_lll_heuristic_dot(const double * vec1, const double * vec2, slong len2, const fmpz_mat_t B, slong k, slong j, slong exp_adj)
    # Computes the dot product of two vectors of doubles ``vec1`` and
    # ``vec2``, which are respectively ``double`` approximations (up to
    # scaling by a power of 2) to rows ``k`` and ``j`` in the exact integer
    # matrix ``B``. If massive cancellation is detected an exact computation
    # is made.
    # The exact computation is scaled by `2^{-\mathtt{exp_adj}}`, where
    # ``exp_adj = r2 + r1`` where `r2` is the exponent for row ``j`` and
    # `r1` is the exponent for row ``k`` (i.e. row ``j`` is notionally
    # thought of as being multiplied by `2^{r2}`, etc.).
    # The final dot product computed by this function is then notionally the
    # return value times `2^{\mathtt{exp_adj}}`.

    int fmpz_lll_check_babai(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s, d_mat_t appB, int *expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
    # Performs floating point size reductions of the ``kappa``-th row of
    # ``B`` by all of the previous rows, uses d_mats ``mu`` and ``r``
    # for storing the GSO data. ``U`` is used to capture the unimodular
    # transformations if it is not `NULL`. The ``double`` array ``s`` will
    # contain the size of the ``kappa``-th row if it were moved into position
    # `i`. The d_mat ``appB`` is an approximation of ``B`` with each row
    # receiving an exponent stored in ``expo`` which gets populated only when
    # needed. The d_mat ``A->appSP`` is an approximation of the Gram matrix
    # whose entries are scalar products of the rows of ``B`` and is used when
    # ``fl->gt`` == `APPROX`. When ``fl->gt`` == `EXACT` the fmpz_mat
    # ``A->exactSP`` (the exact Gram matrix) is used. The index ``a`` is
    # the smallest row index which will be reduced from the ``kappa``-th row.
    # Index ``zeros`` is the number of zero rows in the matrix.
    # ``kappamax`` is the highest index which has been size-reduced so far,
    # and ``n`` is the number of columns you want to consider. ``fl`` is an
    # LLL (L^2) context object. The output is the value -1 if the process fails
    # (usually due to insufficient precision) or 0 if everything was successful.
    # These descriptions will be true for the future Babai procedures as well.

    int fmpz_lll_check_babai_heuristic_d(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s, d_mat_t appB, int *expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
    # Same as :func:`fmpz_lll_check_babai` but using the heuristic inner product
    # rather than a purely floating point inner product. The heuristic will
    # compute at full precision when there is cancellation.

    int fmpz_lll_check_babai_heuristic(int kappa, fmpz_mat_t B, fmpz_mat_t U, mpf_mat_t mu, mpf_mat_t r, mpf *s, mpf_mat_t appB, fmpz_gram_t A, int a, int zeros, int kappamax, int n, mpf_t tmp, mpf_t rtmp, flint_bitcnt_t prec, const fmpz_lll_t fl)
    # This function is like the ``mpf`` version of
    # :func:`fmpz_lll_check_babai_heuristic_d`. However, it also inherits some
    # temporary ``mpf_t`` variables ``tmp`` and ``rtmp``.

    int fmpz_lll_advance_check_babai(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s, d_mat_t appB, int *expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
    # This is a Babai procedure which is used when size reducing a vector beyond
    # an index which LLL has reached. ``cur_kappa`` is the index behind which
    # we can assume ``B`` is LLL reduced, while ``kappa`` is the vector to
    # be reduced. This procedure only size reduces the ``kappa``-th row by
    # vectors up to ``cur_kappa``, **not** ``kappa - 1``.

    int fmpz_lll_advance_check_babai_heuristic_d(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s, d_mat_t appB, int *expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
    # Same as :func:`fmpz_lll_advance_check_babai` but using the heuristic inner
    # product rather than a purely floating point inner product. The heuristic
    # will compute at full precision when there is cancellation.

    int fmpz_lll_shift(const fmpz_mat_t B)
    # Computes the largest number of non-zero entries after the diagonal in
    # ``B``.

    int fmpz_lll_d(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    # This is a mildly greedy version of floating point LLL using doubles only.
    # It tries the fast version of the Babai algorithm
    # (:func:`fmpz_lll_check_babai`). If that fails, then it switches to the
    # heuristic version (:func:`fmpz_lll_check_babai_heuristic_d`) for only one
    # loop and switches right back to the fast version. It reduces ``B`` in
    # place. ``U`` is the matrix used to capture the unimodular
    # transformations if it is not `NULL`. An exception is raised if `U` != `NULL`
    # and ``U->r`` != `d`, where `d` is the lattice dimension. ``fl`` is the
    # context object containing information containing the LLL parameters \delta
    # and \eta. The function can perform reduction on both the lattice basis as
    # well as its Gram matrix. The type of lattice representation can be
    # specified via the parameter ``fl->rt``. The type of Gram matrix to be
    # used in computation (approximate or exact) can also be specified through
    # the variable ``fl->gt`` (applies only if ``fl->rt`` == `Z\_BASIS`).

    int fmpz_lll_d_heuristic(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    # This LLL reduces ``B`` in place using doubles only. It is similar to
    # :func:`fmpz_lll_d` but only uses the heuristic inner products which
    # attempt to detect cancellations.

    int fmpz_lll_mpf2(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_lll_t fl)
    # This is LLL using ``mpf`` with the given precision, ``prec`` for the
    # underlying GSO. It reduces ``B`` in place like the other LLL functions.
    # The `mpf2` in the function name refers to the way the ``mpf_t``'s are
    # initialised.

    int fmpz_lll_mpf(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    # A wrapper of :func:`fmpz_lll_mpf2`. This currently begins with
    # `prec == D\_BITS`, then for the first 20 loops, increases the precision one
    # limb at a time. After 20 loops, it doubles the precision each time. There
    # is a proof that this will eventually work. The return value of this
    # function is 0 if the LLL is successful or -1 if the precision maxes out
    # before ``B`` is LLL-reduced.

    int fmpz_lll_wrapper(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    # A wrapper of the above procedures. It begins with the greediest version
    # (:func:`fmpz_lll_d`), then adapts to the version using heuristic inner
    # products only (:func:`fmpz_lll_d_heuristic`) if ``fl->rt`` == `Z\_BASIS` and
    # ``fl->gt`` == `APPROX`, and finally to the mpf version (:func:`fmpz_lll_mpf`)
    # if needed.
    # ``U`` is the matrix used to capture the unimodular
    # transformations if it is not `NULL`. An exception is raised if `U` != `NULL`
    # and ``U->r`` != `d`, where `d` is the lattice dimension. ``fl`` is the
    # context object containing information containing the LLL parameters \delta
    # and \eta. The function can perform reduction on both the lattice basis as
    # well as its Gram matrix. The type of lattice representation can be
    # specified via the parameter ``fl->rt``. The type of Gram matrix to be
    # used in computation (approximate or exact) can also be specified through
    # the variable ``fl->gt`` (applies only if ``fl->rt`` == `Z\_BASIS`).

    int fmpz_lll_d_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    # Same as :func:`fmpz_lll_d` but with a removal bound, ``gs_B``. The
    # return value is the new dimension of ``B`` if removals are desired.

    int fmpz_lll_d_heuristic_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    # Same as :func:`fmpz_lll_d_heuristic` but with a removal bound,
    # ``gs_B``. The return value is the new dimension of ``B`` if removals
    # are desired.

    int fmpz_lll_mpf2_with_removal(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_t gs_B, const fmpz_lll_t fl)
    # Same as :func:`fmpz_lll_mpf2` but with a removal bound, ``gs_B``. The
    # return value is the new dimension of ``B`` if removals are desired.

    int fmpz_lll_mpf_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    # A wrapper of :func:`fmpz_lll_mpf2_with_removal`. This currently begins
    # with `prec == D\_BITS`, then for the first 20 loops, increases the precision
    # one limb at a time. After 20 loops, it doubles the precision each time.
    # There is a proof that this will eventually work. The return value of this
    # function is the new dimension of ``B`` if removals are desired or -1 if
    # the precision maxes out before ``B`` is LLL-reduced.

    int fmpz_lll_wrapper_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    # A wrapper of the procedures implementing the base case LLL with the
    # addition of the removal boundary. It begins with the greediest version
    # (:func:`fmpz_lll_d_with_removal`), then adapts to the version using
    # heuristic inner products only (:func:`fmpz_lll_d_heuristic_with_removal`)
    # if ``fl->rt`` == `Z\_BASIS` and ``fl->gt`` == `APPROX`, and finally to the mpf
    # version (:func:`fmpz_lll_mpf_with_removal`) if needed.

    int fmpz_lll_d_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    # This is floating point LLL specialized to knapsack-type lattices. It
    # performs early size reductions occasionally which makes things faster in
    # the knapsack case. Otherwise, it is similar to
    # ``fmpz_lll_d_with_removal``.

    int fmpz_lll_wrapper_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    # A wrapper of the procedures implementing the LLL specialized to
    # knapsack-type lattices. It begins with the greediest version and the engine
    # of this version, (:func:`fmpz_lll_d_with_removal_knapsack`), then adapts
    # to the version using heuristic inner products only
    # (:func:`fmpz_lll_d_heuristic_with_removal`) if ``fl->rt`` == `Z\_BASIS` and
    # ``fl->gt`` == `APPROX`, and finally to the mpf version
    # (:func:`fmpz_lll_mpf_with_removal`) if needed.

    int fmpz_lll_with_removal_ulll(fmpz_mat_t FM, fmpz_mat_t UM, slong new_size, const fmpz_t gs_B, const fmpz_lll_t fl)
    # ULLL is a new style of LLL which adjoins an identity matrix to the
    # input lattice ``FM``, then scales the lattice down to ``new_size``
    # bits and reduces this augmented lattice. This tends to be more stable
    # numerically than traditional LLL which means higher dimensions can be
    # attacked using doubles. In each iteration a new identity matrix is adjoined
    # to the truncated lattice. ``UM`` is used to capture the unimodular
    # transformations, while ``gs_B`` and ``fl`` have the same role as in
    # the previous routines. The function is optimised for factoring polynomials.

    bint fmpz_lll_is_reduced_d(const fmpz_mat_t B, const fmpz_lll_t fl)
    bint fmpz_lll_is_reduced_mpfr(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec)
    bint fmpz_lll_is_reduced_d_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd)
    bint fmpz_lll_is_reduced_mpfr_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd, flint_bitcnt_t prec)
    # A non-zero return indicates the matrix is definitely reduced, that is, that
    # * :func:`fmpz_mat_is_reduced` or :func:`fmpz_mat_is_reduced_gram` (for the first two)
    # * :func:`fmpz_mat_is_reduced_with_removal` or :func:`fmpz_mat_is_reduced_gram_with_removal` (for the last two)
    # return non-zero. A zero return value is inconclusive.
    # The `_d` variants are performed in machine precision, while the `_mpfr` uses a precision of `prec` bits.

    bint fmpz_lll_is_reduced(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec)
    bint fmpz_lll_is_reduced_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd, flint_bitcnt_t prec)
    # The return from these functions is always conclusive: the functions
    # * :func:`fmpz_mat_is_reduced` or :func:`fmpz_mat_is_reduced_gram`
    # * :func:`fmpz_mat_is_reduced_with_removal` or :func:`fmpz_mat_is_reduced_gram_with_removal`
    # are optimzied by calling the above heuristics first and returning right away if they give a conclusive answer.

    void fmpz_lll_storjohann_ulll(fmpz_mat_t FM, slong new_size, const fmpz_lll_t fl)
    # Performs ULLL using :func:`fmpz_mat_lll_storjohann` as the LLL function.

    void fmpz_lll(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    # Reduces ``B`` in place according to the parameters specified by the
    # LLL context object ``fl``.
    # This is the main LLL function which should be called by the user. It
    # currently calls the ULLL algorithm (without removals). The ULLL function
    # in turn calls a LLL wrapper which tries to choose an optimal LLL algorithm,
    # starting with a version using just doubles (ULLL tries to maximise usage
    # of this), then a heuristic LLL followed by a full precision floating point
    # LLL if required.
    # ``U`` is the matrix used to capture the unimodular
    # transformations if it is not `NULL`. An exception is raised if `U` != `NULL`
    # and ``U->r`` != `d`, where `d` is the lattice dimension. ``fl`` is the
    # context object containing information containing the LLL parameters \delta
    # and \eta. The function can perform reduction on both the lattice basis as
    # well as its Gram matrix. The type of lattice representation can be
    # specified via the parameter ``fl->rt``. The type of Gram matrix to be
    # used in computation (approximate or exact) can also be specified through
    # the variable ``fl->gt`` (applies only if ``fl->rt`` == `Z\_BASIS`).

    int fmpz_lll_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    # Reduces ``B`` in place according to the parameters specified by the
    # LLL context object ``fl`` and removes vectors whose squared Gram-Schmidt
    # length is greater than the bound ``gs_B``. The return value is the new
    # dimension of ``B`` to be considered for further computation.
    # This is the main LLL with removals function which should be called by
    # the user. Like ``fmpz_lll`` it calls ULLL, but it also sets the
    # Gram-Schmidt bound to that supplied and does removals.
