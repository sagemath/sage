r"""
Sparse Matrices with CMR
"""

from libc.stdint cimport SIZE_MAX

from cysignals.signals cimport sig_on, sig_off

from sage.rings.integer cimport Integer

from .args cimport MatrixArgs_init


cdef class Matrix_cmr_sparse(Matrix_sparse):
    pass


cdef class Matrix_cmr_chr_sparse(Matrix_cmr_sparse):
    r"""
    EXAMPLES::

        sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
        sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
        ....:                           [[1, 2, 3], [4, 0, 6]]); M
        [1 2 3]
        [4 0 6]
        sage: M.dict()
        {(0, 0): 1, (0, 1): 2, (0, 2): 3, (1, 0): 4, (1, 2): 6}

    Matrices of this class are always immutable::

        sage: M.is_immutable()
        True
        sage: copy(M) is M
        True
        sage: deepcopy(M) is M
        True

    This matrix class can only store very small integers::

        sage: Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 1, 3, sparse=True), [-128, 0, 127])
        [-128    0  127]
        sage: Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 1, 3, sparse=True), [126, 127, 128])
        Traceback (most recent call last):
        ...
        OverflowError: value too large to convert to char

    Arithmetic does not preserve the implementation class (even if the numbers would fit)::

        sage: M2 = M + M; M2
        [ 2  4  6]
        [ 8  0 12]
        sage: type(M2)
        <class 'sage.matrix.matrix_integer_sparse.Matrix_integer_sparse'>
        sage: M * 100
        [100 200 300]
        [400   0 600]
    """
    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        TESTS::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 2, 3], [4, 0, 6]]); M
            [1 2 3]
            [4 0 6]
            sage: TestSuite(M).run()
        """
        cdef dict d
        ma = MatrixArgs_init(parent, entries)
        d = ma.dict(coerce)
        self._init_from_dict(d, ma.nrows, ma.ncols)

    cdef _init_from_dict(self, dict d, int nrows, int ncols):
        if cmr == NULL:
            CMRcreateEnvironment(&cmr)
        CMR_CALL(CMRchrmatCreate(cmr, &self._mat, nrows, ncols, len(d)))
        for row in range(nrows):
            self._mat.rowSlice[row] = 0
        for (row, col), coeff in d.items():
            if coeff:
                self._mat.rowSlice[row + 1] += 1
        for row in range(nrows):
            self._mat.rowSlice[row + 1] += self._mat.rowSlice[row]
        for (row, col), coeff in d.items():
            if coeff:
                index = self._mat.rowSlice[row]
                self._mat.entryColumns[index] = col
                self._mat.entryValues[index] = coeff
                self._mat.rowSlice[row] = index + 1
        for row in reversed(range(nrows)):
            self._mat.rowSlice[row + 1] = self._mat.rowSlice[row]
        self._mat.rowSlice[0] = 0
        CMR_CALL(CMRchrmatSortNonzeros(cmr, self._mat))
        self.set_immutable()

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef size_t index
        CMR_CALL(CMRchrmatFindEntry(self._mat, i, j, &index))
        if index == SIZE_MAX:
            return Integer(0)
        return Integer(self._mat.entryValues[index])

    def _dict(self):
        """
        Return the underlying dictionary of ``self``.
        """
        cdef dict d
        d = {}
        for row in range(self._mat.numRows):
            for index in range(self._mat.rowSlice[row], self._mat.rowSlice[row + 1]):
                d[(row, self._mat.entryColumns[index])] = Integer(self._mat.entryValues[index])
        return d

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __dealloc__(self):
        CMRchrmatFree(cmr, &self._mat)

    def _test_change_ring(self, **options):
        return

    def _pickle(self):
        version = 0
        return self._dict(), version

    def _unpickle(self, data, int version):
        if version != 0:
            raise RuntimeError("unknown matrix version (=%s)"%version)
        self._init_from_dict(data, self._nrows, self._ncols)

    # CMR-specific methods. Other classes that want to provide these methods should create
    # a copy of themselves as an instance of this class and delegate to it.

    def is_unimodular(self):
        r"""
        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 0, 0], [0, 1, 0]]); M
            [1 0 0]
            [0 1 0]
            sage: M.is_unimodular()
            True
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 1, 0], [-1, 1, 1]]); M
            [ 1  1  0]
            [-1  1  1]
            sage: M.is_unimodular()
            False
        """
        cdef bint result

        sig_on()
        try:
            CMR_CALL(CMRtestUnimodularity(cmr, self._mat, &result))
        finally:
            sig_off()

        return result

    def is_strongly_unimodular(self):
        cdef bint result

        sig_on()
        try:
            CMR_CALL(CMRtestStrongUnimodularity(cmr, self._mat, &result))
        finally:
            sig_off()

        return result

    def modulus(self):
        cdef bint result
        cdef size_t k

        sig_on()
        try:
            CMR_CALL(CMRtestKmodularity(cmr, self._mat, &result, &k))
        finally:
            sig_off()

        if result:
            return Integer(k)
        else:
            return None

    def strong_modulus(self):
        cdef bint result
        cdef size_t k

        sig_on()
        try:
            CMR_CALL(CMRtestStrongKmodularity(cmr, self._mat, &result, &k))
        finally:
            sig_off()

        if result:
            return Integer(k)
        else:
            return None

    def is_k_modular(self, k):
        return self.modulus() <= k

    def is_strongly_k_modular(self, k):
        return self.strong_modulus() <= k

    def is_graphic(self, *, time_limit=60.0, certificate=False):
        r"""
        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, -1]]); M
            [ 1  0]
            [-1  1]
            [ 0 -1]
            sage: M.is_graphic()  # ?? should it not check 0/1-ness?
            True
            sage: result, certificate = M.is_graphic(certificate=True)
            sage: graph, forest_edges, coforest_edges = certificate
            sage: graph.vertices(sort=True)  # what is the meaning of the numbers?
            [1, 2, 7, 12]
            sage: graph.edges(sort=True, labels=False)
            [(1, 2), (1, 7), (1, 12), (2, 7), (7, 12)]
            sage: forest_edges    # indexed by rows of M
            ((1, 2), (7, 1), (12, 7))
            sage: coforest_edges  # indexed by cols of M
            ((2, 7), (1, 12))

        TESTS::

            sage: M.is_graphic(time_limit=0.0)
            Traceback (most recent call last):
            ...
            RuntimeError: Time limit exceeded
        """
        cdef bint result
        cdef CMR_GRAPH *graph = NULL
        cdef CMR_GRAPH_EDGE* forest_edges = NULL
        cdef CMR_GRAPH_EDGE* coforest_edges = NULL
        cdef CMR_SUBMAT* submatrix = NULL
        cdef CMR_GRAPHIC_STATISTICS stats

        sig_on()
        try:
            if certificate:
                CMR_CALL(CMRtestGraphicMatrix(cmr, self._mat, &result, &graph, &forest_edges,
                                              &coforest_edges, &submatrix, &stats, time_limit))
            else:
                CMR_CALL(CMRtestGraphicMatrix(cmr, self._mat, &result, NULL, NULL,
                                              NULL, NULL, &stats, time_limit))
        finally:
            sig_off()

        if not certificate:
            return result

        # Until we have a proper CMR Graph backend, we just create a Sage graph with whatever backend
        from sage.graphs.graph import Graph

        def vertices():
            i = CMRgraphNodesFirst(graph)
            while CMRgraphNodesValid(graph, i):
                yield i
                i = CMRgraphNodesNext(graph, i)

        def edge(e):
            return Integer(CMRgraphEdgeU(graph, e)), Integer(CMRgraphEdgeV(graph, e))

        def edges():
            i = CMRgraphEdgesFirst(graph)
            while CMRgraphEdgesValid(graph, i):
                e = CMRgraphEdgesEdge(graph, i)
                yield edge(e)
                i = CMRgraphEdgesNext(graph, i)

        if result:
            sage_graph = Graph([list(vertices()), list(edges())])
            sage_forest_edges = tuple(edge(forest_edges[row])
                                      for row in range(self.nrows()))
            sage_coforest_edges = tuple(edge(coforest_edges[column])
                                        for column in range(self.ncols()))
            return True, (sage_graph, sage_forest_edges, sage_coforest_edges)

        return False, None  # submatrix TBD

    def is_cographic(self, *, time_limit=60.0, certificate=False):
         raise NotImplementedError

    def is_network_matrix(self, *, time_limit=60.0, certificate=False):
        raise NotImplementedError

    def is_dual_network_matrix(self, *, time_limit=60.0, certificate=False):
        raise NotImplementedError

    def _is_binary_linear_matroid_regular(self, *, time_limit=60.0, certificate=False,
                                          use_direct_graphicness_test=True,
                                          series_parallel_ok=True,
                                          check_graphic_minors_planar=False,
                                          complete_tree='if_regular',
                                          construct_matrices=False,
                                          construct_transposes=False,
                                          construct_graphs=False):
        r"""
        This is an internal method because it should really be exposed
        as a method of :class:`Matroid`.

        INPUT:

        - ``certificate``: One of ``False``, ``True``, ``'if_regular'``, ``'if_not_regular'``

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [1, 1], [0, 1]]); M
            [1 0]
            [1 1]
            [0 1]
            sage: M._is_binary_linear_matroid_regular()
            True
        """
        cdef bint result
        cdef CMR_REGULAR_PARAMETERS params
        cdef CMR_REGULAR_STATISTICS stats
        cdef CMR_DEC *dec = NULL
        cdef CMR_MINOR *minor = NULL

        cdef CMR_DEC **pdec = &dec
        cdef CMR_MINOR **pminor = &minor

        cdef dict kwds = dict(use_direct_graphicness_test=use_direct_graphicness_test,
                              series_parallel_ok=series_parallel_ok,
                              check_graphic_minors_planar=check_graphic_minors_planar,
                              complete_tree=complete_tree,
                              construct_matrices=construct_matrices,
                              construct_transposes=construct_transposes,
                              construct_graphs=construct_graphs)

        _set_cmr_regular_parameters(&params, kwds)
        sig_on()
        try:
            CMR_CALL(CMRtestBinaryRegular(cmr, self._mat, &result, pdec, pminor,
                                          &params, &stats, time_limit))
        finally:
            sig_off()

        if not certificate:
            return result

        raise NotImplementedError

    def is_totally_unimodular(self, *, time_limit=60.0, certificate=False,
                              use_direct_graphicness_test=True,
                              series_parallel_ok=True,
                              check_graphic_minors_planar=False,
                              complete_tree='if_regular',
                              construct_matrices=False,
                              construct_transposes=False,
                              construct_graphs=False):
        r"""
        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: M.is_totally_unimodular()
            True

        """
        cdef bint result
        cdef CMR_TU_PARAMETERS params
        cdef CMR_TU_STATISTICS stats
        cdef CMR_DEC *dec = NULL
        cdef CMR_SUBMAT *submat = NULL

        cdef CMR_DEC **pdec = &dec
        cdef CMR_SUBMAT **psubmat = &submat

        cdef dict kwds = dict(use_direct_graphicness_test=use_direct_graphicness_test,
                              series_parallel_ok=series_parallel_ok,
                              check_graphic_minors_planar=check_graphic_minors_planar,
                              complete_tree=complete_tree,
                              construct_matrices=construct_matrices,
                              construct_transposes=construct_transposes,
                              construct_graphs=construct_graphs)

        _set_cmr_regular_parameters(&params.regular, kwds)
        sig_on()
        try:
            CMR_CALL(CMRtestTotalUnimodularity(cmr, self._mat, &result, pdec, psubmat,
                                               &params, &stats, time_limit))
        finally:
            sig_off()

        if not certificate:
            return result

        raise NotImplementedError


cdef _cmr_dec_construct(param):
    if not param:
        return CMR_DEC_CONSTRUCT_NONE
    if param == 'leaves':
        return CMR_DEC_CONSTRUCT_LEAVES
    return CMR_DEC_CONSTRUCT_ALL


cdef _set_cmr_regular_parameters(CMR_REGULAR_PARAMETERS *params, dict kwds):
    CMR_CALL(CMRparamsRegularInit(params))
    params.directGraphicness = kwds['use_direct_graphicness_test']
    params.seriesParallel = kwds['series_parallel_ok']
    params.planarityCheck = kwds['check_graphic_minors_planar']
    params.completeTree = kwds['complete_tree'] is True
    params.matrices = _cmr_dec_construct(kwds['construct_matrices'])
    params.transposes = _cmr_dec_construct(kwds['construct_transposes'])
    params.graphs = _cmr_dec_construct(kwds['construct_graphs'])
