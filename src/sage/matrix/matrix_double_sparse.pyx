from .matrix2 cimport Matrix
from .matrix_generic_sparse cimport Matrix_generic_sparse

cdef class Matrix_double_sparse(Matrix_generic_sparse):
    r"""
    Class for sparse RDF/CDF matrices.

    EXAMPLES::

        sage: A = matrix.random(RDF, ZZ.random_element(5), sparse=True)
        sage: A.__class__
        <class 'sage.matrix.matrix_double_sparse.Matrix_double_sparse'>
        sage: A = matrix.random(CDF, ZZ.random_element(5), sparse=True)
        sage: A.__class__
        <class 'sage.matrix.matrix_double_sparse.Matrix_double_sparse'>

    """
    def cholesky(self):
        r"""
        Returns the Cholesky decomposition of a Hermitian matrix.

        INPUT:

        A positive-definite matrix over ``RDF`` or ``CDF``.

        OUTPUT:

        For a matrix `A` the routine returns a lower triangular
        matrix `L` (over the same base ring as `A`) such that,

        .. MATH::

            A = LL^\ast

        where `L^\ast` is the conjugate-transpose. If the matrix is
        not positive-definite (for example, if it is not Hermitian)
        then a ``ValueError`` results.

        ALGORITHM:

        If cvxopt is available, we use its sparse `cholmod` routines
        to compute the factorization quickly. Otherwise, we fall back
        to the naive implementation used for dense matrices.

        EXAMPLES::

            sage: A = matrix(RDF, [[10, 0, 3, 0],
            ....:                  [ 0, 5, 0,-2],
            ....:                  [ 3, 0, 5, 0],
            ....:                  [ 0,-2, 0, 2]],
            ....:                 sparse=True)
            sage: L = A.cholesky()
            sage: L.is_triangular()
            True
            sage: (A - L*L.T).norm(1) < 1e-10
            True

        ::

            sage: A = matrix(CDF, [[        2,   4 + 2*I,   6 - 4*I],
            ....:                  [ -2*I + 4,        11, 10 - 12*I],
            ....:                  [  4*I + 6, 10 + 12*I,        37]])
            sage: L = A.cholesky()
            sage: L.is_triangular()
            True
            sage: (A - L*L.H).norm(1) < 1e-10
            True

        TESTS:

        Test the properties of a Cholesky factorization using "random"
        symmetric/Hermitian positive-definite matrices. We also
        compare with the factorizations obtained over ``RR`` or
        ``CC``, which (when cvxopt is available) use a different
        implementation. This ensures that both implementations return
        comparable answers::

            sage: n = ZZ.random_element(1,5)
            sage: A = matrix.random(RDF, n, sparse=True)
            sage: I = matrix.identity(RDF, n, sparse=True)
            sage: A = A*A.transpose() + I
            sage: L = A.cholesky()
            sage: (A - L*L.T).norm(1) < 1e-10
            True
            sage: B = A.change_ring(RR)
            sage: (B.cholesky() - L).norm(1) < 1e-10
            True

        ::

            sage: n = ZZ.random_element(1,5)
            sage: A = matrix.random(CDF, n, sparse=True)
            sage: I = matrix.identity(CDF, n, sparse=True)
            sage: A = A*A.conjugate_transpose() + I
            sage: L = A.cholesky()
            sage: (A - L*L.H).norm(1) < 1e-10
            True
            sage: B = A.change_ring(CC)
            sage: (B.cholesky() - L).norm(1) < 1e-10
            True
        """
        cdef Matrix L # output matrix

        L = self.fetch('cholesky')
        if L is not None:
            return L

        # The superclass method does this, so we should too...
        if not self.is_hermitian():
            raise ValueError("matrix is not Hermitian")

        try:
            from cvxopt import cholmod, spmatrix, matrix as cvxopt_matrix
        except ModuleNotFoundError:
            # Sage can be built with --disable-cvxopt, so we have to
            # handle the case where cvxopt is not present. The
            # superclass method is slow, but no longer raises an
            # error, so let's try that.
            L = super().cholesky()

        cdef list idx_pairs = self.nonzero_positions(copy=False)
        cdef list row_idxs = [r for (r, c) in idx_pairs]
        cdef list col_idxs = [c for (r, c) in idx_pairs]
        value_type = float
        type_code = 'd'
        from sage.rings.complex_double import CDF
        if self.base_ring() is CDF:
            value_type = complex
            type_code = 'z'
        cdef list values = [value_type(self[idx]) for idx in idx_pairs]

        cvx_self = spmatrix(values,
                            row_idxs,
                            col_idxs,
                            size=self.dimensions(),
                            tc=type_code)

        # Insist that our ordering (= permutation) is used, and then
        # choose the identity permutation. In SageMath, cholesky()
        # does not allow for permutations. The 'postorder' option
        # I only found documented in the cholmod module's docstring.
        cholmod.options['nmethods'] = 1
        cholmod.options['postorder'] = False
        id_order = cvxopt_matrix(range(self.nrows()))

        # WARNING: the getfactor() interface used below is undocumented:
        #
        #   https://groups.google.com/g/cvxopt/c/xQ-lR9ESijg/discussion
        #
        # And while upstream suggests that we should use a simplicial
        # factorization, the cvxopt docs themselves say that we should
        # use the default (=2) value of 'supernodal' for the PAP^T =
        # LL^T factorization that we want.
        #
        cvx_symbolic = cholmod.symbolic(cvx_self, p=id_order)

        # This overwrites cvx_symbolic with the numeric factorization before
        # passing it to getfactor().
        cholmod.numeric(cvx_self, cvx_symbolic)
        cvx_L = cholmod.getfactor(cvx_symbolic)

        # The (I[k],J[k]) entry of cvx_L has value V[k]. But beware that V
        # contains only the non-zero entries of the matrix; as a result, the
        # dict below contains keys only for those nonzero entries.
        L = self.matrix_space()({
          (cvx_L.I[k], cvx_L.J[k]): cvx_L.V[k]
          for k in range(len(cvx_L.V))
        })

        self.cache('cholesky', L)
        return L
