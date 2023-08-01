# sage_setup: distribution = sagemath-cmr
r"""
Sparse Matrices with CMR
"""

from libc.stdint cimport SIZE_MAX

from cysignals.signals cimport sig_on, sig_off

from sage.libs.cmr.cmr cimport *
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.structure.element cimport Matrix

from .args cimport MatrixArgs_init
from .constructor import matrix
from .matrix_space import MatrixSpace
from .seymour_decomposition cimport create_DecompositionNode


cdef class Matrix_cmr_sparse(Matrix_sparse):
    r"""
    Base class for sparse matrices implemented in CMR
    """
    pass


cdef class Matrix_cmr_chr_sparse(Matrix_cmr_sparse):
    r"""
    Sparse matrices with 8-bit integer entries implemented in CMR

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
    def __init__(self, parent, entries=None, copy=None, bint coerce=True, immutable=True):
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
        self._init_from_dict(d, ma.nrows, ma.ncols, immutable=True)

    cdef _init_from_dict(self, dict d, int nrows, int ncols, bint immutable=True):
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
        if immutable:
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
        if self._root is None or self._root is self:
            # We own it, so we have to free it.
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

    @staticmethod
    def _from_data(data, immutable=True):
        if not isinstance(data, Matrix):
            data = matrix(ZZ, data, sparse=True)
        if not isinstance(data, Matrix_cmr_chr_sparse):
            data = Matrix_cmr_chr_sparse(data.parent(), data, immutable=immutable)
        return data

    @staticmethod
    cdef _from_cmr(CMR_CHRMAT *mat, bint immutable=False):
        r"""
        INPUT:

        - ``mat`` -- a ``CMR_CHRMAT``; after this call, it is owned by the created Python object

        OUTPUT: A :class:`Matrix_cmr_chr_sparse`

        """
        cdef Matrix_cmr_chr_sparse result
        ms = MatrixSpace(ZZ, mat.numRows, mat.numColumns, sparse=True)
        result = Matrix_cmr_chr_sparse.__new__(Matrix_cmr_chr_sparse, ms, immutable=immutable)
        result._mat = mat
        result._root = None
        return result

    @staticmethod
    def _network_matrix_from_digraph(digraph, forest_arcs=None, arcs=None, vertices=None):
        r"""
        Return the network matrix of ``digraph``, pivoted according to ``forest_arcs``.

        Its rows are indexed parallel to ``forest_arcs``.
        It is in "short tableau" form, i.e., the columns are indexed parallel
        to the elements of ``arcs`` that are not in ``forest_arcs``.

        .. NOTE::

            In [Sch1986]_, the columns are indexed by all arcs of the digraph,
            giving a "long tableau" form of the network matrix.

        INPUT:

        - ``digraph`` -- a :class:`DiGraph`

        - ``forest_arcs`` -- a sequence of arcs of the ``digraph`` or ``None`` (the default:
          use the labels of the ``arcs`` as a boolean value)

        - ``arcs`` -- a sequence of arcs of the digraph or ``None`` (the default:
          all arcs going out from the ``vertices``)

        - ``vertices`` -- a sequence of vertices of the digraph or ``None`` (the default:
          all vertices of the ``digraph``)

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse

        Defining the forest by arc labels::

            sage: D = DiGraph([[0, 1, 2, 3],
            ....:              [(0, 1, True), (0, 2, True), (1, 2), (1, 3, True), (2, 3)]])
            sage: M = Matrix_cmr_chr_sparse._network_matrix_from_digraph(D); M
            [ 1 -1]
            [-1  1]
            [ 1  0]

        Defining the forest by a separate list of forest arcs::

            sage: D = DiGraph([[0, 1, 2, 3],
            ....:              [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]])
            sage: T = [(0, 1), (0, 2), (1, 3)]
            sage: M = Matrix_cmr_chr_sparse._network_matrix_from_digraph(D, T); M
            [ 1 -1]
            [-1  1]
            [ 1  0]

        Prescribing an order for the arcs (columns)::

            sage: D = DiGraph([[0, 1, 2, 3],
            ....:              [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]])
            sage: T = [(0, 1), (0, 2), (1, 3)]
            sage: A = [(2, 3), (0, 1), (0, 2), (1, 2), (1, 3)]
            sage: M = Matrix_cmr_chr_sparse._network_matrix_from_digraph(D, T, arcs=A); M
            [ 1 -1]
            [-1  1]
            [ 1  0]

        TESTS::

            sage: D = DiGraph([[0, 1, 2], [(0, 1), (1, 2), (0, 2)]])
            sage: not_a_forest = [(0, 1), (1, 2), (0, 2)]
            sage: M = Matrix_cmr_chr_sparse._network_matrix_from_digraph(D, not_a_forest); M
            Traceback (most recent call last):
            ...
            ValueError: not a spanning forest
        """
        cdef CMR_GRAPH *cmr_digraph = NULL
        cdef dict vertex_to_cmr_node = {}
        cdef dict arc_to_cmr_edge = {}
        cdef CMR_GRAPH_EDGE cmr_edge
        cdef CMR_GRAPH_NODE cmr_node

        if cmr == NULL:
            CMRcreateEnvironment(&cmr)

        CMR_CALL(CMRgraphCreateEmpty(cmr, &cmr_digraph, digraph.num_verts(), digraph.num_edges()))

        if vertices is None:
            iter_vertices = digraph.vertex_iterator()
        else:
            iter_vertices = vertices
        for u in iter_vertices:
            CMR_CALL(CMRgraphAddNode(cmr, cmr_digraph, &cmr_node))
            vertex_to_cmr_node[u] = cmr_node

        vertices = vertex_to_cmr_node.keys()
        if arcs is None:
            arcs = digraph.edge_iterator(labels=False, vertices=vertices, ignore_direction=False)

        for a in arcs:
            u, v = a
            CMR_CALL(CMRgraphAddEdge(cmr, cmr_digraph, vertex_to_cmr_node[u],
                                     vertex_to_cmr_node[v], &cmr_edge))
            arc_to_cmr_edge[(u, v)] = cmr_edge

        cdef CMR_GRAPH_EDGE *cmr_forest_arcs = NULL
        cdef bool *cmr_arcs_reversed = NULL
        cdef CMR_CHRMAT *cmr_matrix = NULL
        cdef bool is_correct_forest
        cdef size_t num_forest_arcs

        cdef size_t mem_arcs
        if forest_arcs is not None:
            mem_arcs = len(forest_arcs)
        else:
            mem_arcs = len(vertices) - 1

        CMR_CALL(_CMRallocBlockArray(cmr, <void **> &cmr_forest_arcs, mem_arcs, sizeof(CMR_GRAPH_EDGE)))
        try:
            if forest_arcs is None:
                num_forest_arcs = 0
                for u, v, label in digraph.edge_iterator(labels=True, vertices=vertices,
                                                         ignore_direction=False):
                    if label:
                        if num_forest_arcs >= mem_arcs:
                            raise ValueError('not a spanning forest')
                        cmr_forest_arcs[num_forest_arcs] = arc_to_cmr_edge[(u, v)]
                        num_forest_arcs += 1
            else:
                num_forest_arcs = len(forest_arcs)
                for i, (u, v) in enumerate(forest_arcs):
                    cmr_forest_arcs[i] = arc_to_cmr_edge[(u, v)]

            CMR_CALL(_CMRallocBlockArray(cmr, <void **> &cmr_arcs_reversed, len(arc_to_cmr_edge), sizeof(bool)))
            try:
                for j in range(len(arc_to_cmr_edge)):
                    cmr_arcs_reversed[j] = <bool> False
                CMR_CALL(CMRcomputeNetworkMatrix(cmr, cmr_digraph, &cmr_matrix, NULL,
                                                 cmr_arcs_reversed, num_forest_arcs, cmr_forest_arcs,
                                                 0, NULL, &is_correct_forest))
                if not is_correct_forest:
                    raise ValueError('not a spanning forest')
            finally:
                CMR_CALL(_CMRfreeBlockArray(cmr, <void **> &cmr_arcs_reversed))
        finally:
            CMR_CALL(_CMRfreeBlockArray(cmr, <void **> &cmr_forest_arcs))

        return Matrix_cmr_chr_sparse._from_cmr(cmr_matrix)

    @staticmethod
    def one_sum(*summands):
        r"""
        Return a block-diagonal matrix constructed from the given matrices (summands).

        INPUT:

        - ``summands`` -- integer matrices or data from which integer matrices can be
          constructed

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]], [[1, 1], [-1, 0]])
            sage: unicode_art(M)
            ⎛ 1  0│ 0  0⎞
            ⎜-1  1│ 0  0⎟
            ⎜─────┼─────⎟
            ⎜ 0  0│ 1  1⎟
            ⎝ 0  0│-1  0⎠
            sage: M = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]], [[1]], [[-1]])
            sage: unicode_art(M)
            ⎛ 1  0│ 0│ 0⎞
            ⎜-1  1│ 0│ 0⎟
            ⎜─────┼──┼──⎟
            ⎜ 0  0│ 1│ 0⎟
            ⎜─────┼──┼──⎟
            ⎝ 0  0│ 0│-1⎠

        TESTS::

            sage: M = Matrix_cmr_chr_sparse.one_sum(); M
            []
            sage: M.parent()
            Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring

            sage: M = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]]); unicode_art(M)
            ⎛ 1  0⎞
            ⎝-1  1⎠
            sage: M.parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        """
        cdef Matrix_cmr_chr_sparse sum, summand
        cdef CMR_CHRMAT *sum_mat
        summands = iter(summands)
        try:
            first = next(summands)
        except StopIteration:
            return Matrix_cmr_chr_sparse._from_data({})
        sum = Matrix_cmr_chr_sparse._from_data(first, immutable=False)
        sum_mat = sum._mat
        row_subdivision = []
        column_subdivision = []
        for s in summands:
            row_subdivision.append(sum_mat.numRows)
            column_subdivision.append(sum_mat.numColumns)
            summand = Matrix_cmr_chr_sparse._from_data(s)
            CMR_CALL(CMRoneSum(cmr, sum_mat, summand._mat, &sum_mat))
        if sum_mat != sum._mat:
            sum = Matrix_cmr_chr_sparse._from_cmr(sum_mat, immutable=False)
        if row_subdivision or column_subdivision:
            sum.subdivide(row_subdivision, column_subdivision)
        sum.set_immutable()
        return sum

    def two_sum(first_mat, second_mat, column, row):
        r"""
        Return the 2-sum matrix constructed from the given matrices ``first_mat`` and ``second_mat``, with column index of the first matrix ``column`` and row index of the second matrix ``row``.
        Suppose that ``column`` indicates the last column and ``row`` indicates the first row, i.e, the first matrix is `M_1=\begin{bmatrix} A & a\end{bmatrix}` and the second matrix is `M_2=\begin{bmatrix} b^T \\ B\end{bmatrix}`. Then the two sum `M_1 \oplus_2 M_2 =\begin{bmatrix}A & ab^T\\ 0 & B\end{bmatrix}`.

        INPUT:

        - ``first_mat`` -- the first integer matrix
        - ``second_mat`` -- the second integer matrix
        - ``column`` -- the column index of the first integer matrix
        - ``row`` -- the row index of the first integer matrix

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M1 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                            [[1, 2, 3], [4, 5, 6]]); M1
            [1 2 3]
            [4 5 6]
            sage: M2 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                            [[7, 8, 9], [-1, -2, -3]]); M2
            [ 7  8  9]
            [-1 -2 -3]
            sage: Matrix_cmr_chr_sparse.two_sum(M1,M2,2,0)
            [ 1  2|21 24 27]
            [ 4  5|42 48 54]
            [-----+--------]
            [ 0  0|-1 -2 -3]
            sage: M1.two_sum(M2,1,1)
            [  1   3| -2  -4  -6]
            [  4   6| -5 -10 -15]
            [-------+-----------]
            [  0   0|  7   8   9]
        """
        cdef Matrix_cmr_chr_sparse sum, first, second
        cdef CMR_CHRMAT *sum_mat
        first = Matrix_cmr_chr_sparse._from_data(first_mat, immutable=False)
        second = Matrix_cmr_chr_sparse._from_data(second_mat, immutable=False)
        if column < 0 or column >= first._mat.numColumns:
            raise ValueError("First marker should be a column index of the first matrix")
        if row < 0 or row >= second._mat.numRows:
            raise ValueError("Second marker should be a row index of the second matrix")
        row_subdivision = []
        column_subdivision = []
        row_subdivision.append(first._mat.numRows)
        column_subdivision.append(first._mat.numColumns - 1)
        CMR_CALL(CMRtwoSum(cmr, first._mat, second._mat, CMRcolumnToElement(column), CMRrowToElement(row), &sum_mat))
        sum = Matrix_cmr_chr_sparse._from_cmr(sum_mat, immutable=False)
        if row_subdivision or column_subdivision:
            sum.subdivide(row_subdivision, column_subdivision)
        sum.set_immutable()
        return sum

    def three_sum(first_mat, second_mat, first_col_index1, first_col_index2, second_col_index1, second_col_index2):
        r"""
        Return the 3-sum matrix constructed from the given matrices ``first_mat`` and ``second_mat``, with 'first_col_index1'
        and 'first_col_index2' being the indices of the column vectors of the matrix, which are identical except for one row
        having a 0 in one column and the other a non-zero entry in that row. The method assumes the nonzero entry is one. The same assumptions
        are made for 'second_mat' and its input index variables.
        
        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M1 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
            ....:                          [[1, 0, -1, 0, 1], [1, 1, 0, -1, 1], [0, 0, 1, 1, 1], 
            ....:                           [1, 1, -1, 0, 0], [-1, -1, 0, 0,1]]); M1
            [ 1  0 -1  0  1]
            [ 1  1  0 -1  1]
            [ 0  0  1  1  1]
            [ 1  1 -1  0  0]
            [-1 -1  0  0  1]
            sage: M2 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
            ....:                      [[1, 1, 1, 1, 1], [1, 1, 1, 0, 0], [1, 0, 1, 1, 0],
            ....:                       [0, 0, 0, 1, 1], [1, 1, 0, 0, 1]]); M2
            [1 1 1 1 1]
            [1 1 1 0 0]
            [1 0 1 1 0]
            [0 0 0 1 1]
            [1 1 0 0 1]
            sage: M3 = Matrix_cmr_chr_sparse.three_sum(M1, M2, 0, 1, 0, 1); M3
            [ 0 -1  1  1  1  0]
            [ 1  1  1  0  0  0]
            [-1  0  0  1  1  0]
            [ 0  0  1 -1 -1  0]
            [-1  0  1  1  1  1]
            [-1  0  1  1  0  0]
            [-1  0  1  0  1  1]
            [ 1  0 -1  0  0  1]
        """
        fc = len(first_mat.columns())
        sc = len(second_mat.columns())
        fr = len(first_mat.rows())
        sr = len(second_mat.rows())
        if any([fc < 3, sc < 3, fr < 2, sr < 2]):
            raise ValueError('Some matrix is not large enough to perform a 3-sum')
        if any([first_col_index1 >= fc, first_col_index2 >= fc, second_col_index1 >= sc, second_col_index2 >= sc]):
            raise ValueError('Some column indicated exceeds its matrix size')
        first_col1 = first_mat.columns()[first_col_index1]
        first_col2 = first_mat.columns()[first_col_index2]
        second_col1 = second_mat.columns()[second_col_index1]
        second_col2 = second_mat.columns()[second_col_index2]
        fir_nrows = range(fr)
        sec_nrows = range(sr)
        valid1 = False
        valid2 = False
        for i in fir_nrows:
            if (first_col1[i] == 1 and first_col2[i] == 0) or (first_col1[i] == 0 and first_col2[i] == 1):
                subcol1 = tuple(first_col1[k] for k in fir_nrows if k != i)
                subcol2 = tuple(first_col2[k] for k in fir_nrows if k != i)
                if subcol1 == subcol2:
                    valid1 = True
                    first_row_index = i 
                    break
        for i in sec_nrows:
            if (second_col1[i] == 1 and second_col2[i] == 0) or (second_col1[i] == 0 and second_col2[i] == 1):
                subcol1 = tuple(second_col1[k] for k in sec_nrows if k != i)
                subcol2 = tuple(second_col2[k] for k in sec_nrows if k != i)
                if subcol1 == subcol2:
                    valid2 = True
                    second_row_index = i 
                    break
        if not (valid1 and valid2):
            raise ValueError('indicated columns of Matrices are not of appropriate form for 3-sum')
        first_subcol = first_mat.delete_rows([first_row_index]).columns()[first_col_index1]
        second_subcol = first_mat.delete_rows([second_row_index]).columns()[second_col_index1]
        first_submat = first_mat.delete_columns([first_col_index1, first_col_index2])
        second_submat = second_mat.delete_columns([second_col_index1, second_col_index2])
        first_row = first_submat.rows()[first_row_index]
        second_row = second_submat.rows()[second_row_index]
        first_submat = first_submat.delete_rows([first_row_index])
        second_submat = second_submat.delete_rows([second_row_index])
        first_subrows = first_submat.rows()
        second_subrows = second_submat.rows()
        upper_right_rows = first_subcol.tensor_product(second_row).rows()
        lower_left_rows = second_subcol.tensor_product(first_row).rows()
        n1 = len(first_submat.rows())
        n2 = len(second_submat.rows())
        row_list = []
        for i in range(n1):
            r = list(first_subrows[i])
            u = list(upper_right_rows[i])
            r.extend(u)
            row_list.append(r)
        for i in range(n2):
            r = list(lower_left_rows[i])
            u = list(second_subrows[i])
            r.extend(u)
            row_list.append(r)
        return Matrix_cmr_chr_sparse._from_data(row_list, immutable = False)  

    def delete_rows(self, indices):
        rows = self.rows()
        row_list = []
        n = len(rows)
        for i in indices:
            if i >= n:
                raise ValueError('Found index greater than matrix size')
            rows.pop(i)
        for r in rows:
            x = []
            for i in range(len(r)):
                x.append(r[i])
            row_list.append(x)
        return Matrix_cmr_chr_sparse._from_data(row_list, immutable = False)

    def delete_columns(self, indices):
        rows = self.rows()
        n = len(rows)
        row_list = []
        for i in indices:
            if i >= n:
                raise ValueError('Found index greater than matrix size')
        for r in rows:
            x = []
            for k in range(len(r)):
                if not (k in indices):
                    x.append(r[k])
            row_list.append(x)
        return Matrix_cmr_chr_sparse._from_data(row_list, immutable = False)

    def is_unimodular(self):
        r"""
        Return whether ``self`` is a unimodular matrix.

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
        cdef bool result

        sig_on()
        try:
            CMR_CALL(CMRtestUnimodularity(cmr, self._mat, &result))
        finally:
            sig_off()

        return <bint> result

    def is_strongly_unimodular(self):
        r"""
        Return whether ``self`` is a strongly unimodular matrix.

        A matrix is strongly unimodular if ``self`` and ``self.transpose()`` are both unimodular.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.is_unimodular()
            True
            sage: M.is_strongly_unimodular()
            False
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 0, 0], [0, 1, 0]]); M
            [1 0 0]
            [0 1 0]
            sage: M.is_strongly_unimodular()
            True
        """
        cdef bool result

        sig_on()
        try:
            CMR_CALL(CMRtestStrongUnimodularity(cmr, self._mat, &result))
        finally:
            sig_off()

        return <bint> result

    def modulus(self):
        r"""
        Return the integer `k` such that ``self`` is `k`-modular.

        A matrix `M` of rank `r` is `k`-modular if the following two conditions are satisfied:
        - for some column basis `B` of `M`, the greatest common divisor of the determinants of all `r`-by-`r` submatrices of `B` is `k`;
        - the matrix `X` such that `M=BX` is totally unimodular.

        OUTPUT:

        - ``k``: ``self`` is  `k`-modular
        - ``None``: ``self`` is not `k`-modular for any `k`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.modulus()
            1
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 1, 1], [0, 1, 3]]); M
            [1 1 1]
            [0 1 3]
            sage: M.modulus()
        """
        cdef bool result
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
        r"""
        Return the integer `k` such that ``self`` is strongly `k`-modular.

        A matrix `M` of rank-`r` is `k`-modular if the following two conditions are satisfied:
        - for some column basis `B` of `M`, the greatest common divisor of the determinants of all `r`-by-`r` submatrices of `B` is `k`;
        - the matrix `X` such that `M=BX` is totally unimodular.

        OUTPUT:

        - ``k``: ``self`` is  `k`-modular
        - ``None``: ``self`` is not `k`-modular for any `k`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.strong_modulus()
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 0, 0], [0, 1, 0]]); M
            [1 0 0]
            [0 1 0]
            sage: M.strong_modulus()
            1
        """
        cdef bool result
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
        r"""
        Return whether ``self`` is `k`-modular.

        A matrix `M` of rank-`r` is `k`-modular if the following two conditions are satisfied:
        - for some column basis `B` of `M`, the greatest common divisor of the determinants of all `r`-by-`r` submatrices of `B` is `k`;
        - the matrix `X` such that `M=BX` is totally unimodular.
        If the matrix has full row rank, it is `k`-modular if all the full rank minor of the matrix has determinant `0,\pm k`.
        The matrix is also called strictly `k`-modular.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.is_k_modular(1)
            True
            sage: M.is_k_modular(2)
            False
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 1, 1], [0, 1, 3]]); M
            [1 1 1]
            [0 1 3]
            sage: M.is_k_modular(1)
            False
        """
        result = self.modulus()
        if not result:
            return False
        else:
            return result == k

    def is_strongly_k_modular(self, k):
        r"""
        Return whether ``self`` is strongly `k`-modular.

        A matrix is strongly `k`-modular if ``self`` and ``self.transpose()`` are both `k`-modular.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.is_strongly_k_modular(1)
            False
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 0, 0], [0, 1, 0]]); M
            [1 0 0]
            [0 1 0]
            sage: M.is_strongly_k_modular(1)
            True
        """
        result = self.strong_modulus()
        if not result:
            return False
        else:
            return result == k

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
        cdef bool result
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
            return <bint> result

        if <bint> result:
            sage_graph = _sage_graph(graph)
            sage_forest_edges = tuple(_sage_edge(graph, forest_edges[row])
                                      for row in range(self.nrows()))
            sage_coforest_edges = tuple(_sage_edge(graph, coforest_edges[column])
                                        for column in range(self.ncols()))
            return True, (sage_graph, sage_forest_edges, sage_coforest_edges)

        return False, NotImplemented  # submatrix TBD

    def is_cographic(self, *, time_limit=60.0, certificate=False):
        r"""
        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 4, 9, sparse=True), [[1, 0, 0, 0, 1, -1, 1, 0, 0],
            ....:                                  [0, 1, 0, 0, 0, 1, -1, 1, 0], [0, 0, 1, 0, 0, 0, 1, -1, 1], 
            ....:                                  [0, 0, 0, 1, 1, 0, 0, 1, -1]]); M
            [ 1  0  0  0  1 -1  1  0  0]
            [ 0  1  0  0  0  1 -1  1  0]
            [ 0  0  1  0  0  0  1 -1  1]
            [ 0  0  0  1  1  0  0  1 -1]
            sage: M.is_cographic()
            True
        """
        cdef bool result
        cdef CMR_GRAPH *graph = NULL
        cdef CMR_GRAPH_EDGE* forest_edges = NULL
        cdef CMR_GRAPH_EDGE* coforest_edges = NULL
        cdef CMR_SUBMAT* submatrix = NULL
        cdef CMR_GRAPHIC_STATISTICS stats

        sig_on()
        try:
            if certificate:
                CMR_CALL(CMRtestCographicMatrix(cmr, self._mat, &result, &graph, &forest_edges,
                                              &coforest_edges, &submatrix, &stats, time_limit))
            else:
                CMR_CALL(CMRtestCographicMatrix(cmr, self._mat, &result, NULL, NULL,
                                              NULL, NULL, &stats, time_limit))
        finally:
            sig_off()

        if not certificate:
            return <bint> result

        if <bint> result:
            sage_graph = _sage_graph(graph)
            sage_forest_edges = tuple(_sage_edge(graph, forest_edges[row])
                                      for row in range(self.nrows()))
            sage_coforest_edges = tuple(_sage_edge(graph, coforest_edges[column])
                                        for column in range(self.ncols()))
            return True, (sage_graph, sage_forest_edges, sage_coforest_edges)

        return False, NotImplemented  # submatrix TBD


    def is_network_matrix(self, *, time_limit=60.0, certificate=False):
        r"""
        EXAMPLES:

        This is test ``Basic`` in CMR's ``test_network.cpp``::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 6, 7, sparse=True),
            ....:                           [[-1,  0,  0,  0,  1, -1,  0],
            ....:                            [ 1,  0,  0,  1, -1,  1,  0],
            ....:                            [ 0, -1,  0, -1,  1, -1,  0],
            ....:                            [ 0,  1,  0,  0,  0,  0,  1],
            ....:                            [ 0,  0,  1, -1,  1,  0,  1],
            ....:                            [ 0,  0, -1,  1, -1,  0,  0]])
            sage: M.is_network_matrix()
            True
            sage: result, certificate = M.is_network_matrix(certificate=True)
            sage: result, certificate
            (True,
             (Digraph on 7 vertices,
              ((9, 8), (3, 8), (3, 4), (5, 4), (4, 6), (0, 6)),
              ((3, 9), (5, 3), (4, 0), (0, 8), (9, 0), (4, 9), (5, 6))))
            sage: digraph, forest_arcs, coforest_arcs = certificate
            sage: list(digraph.edges(sort=True))
            [(0, 6, None), (0, 8, None),
             (3, 4, None), (3, 8, None), (3, 9, None),
             (4, 0, None), (4, 6, None), (4, 9, None),
             (5, 3, None), (5, 4, None), (5, 6, None),
             (9, 0, None), (9, 8, None)]
            sage: digraph.plot(edge_colors={'red': forest_arcs})                        # needs sage.plot
            Graphics object consisting of 21 graphics primitives
        """
        cdef bool result
        cdef CMR_GRAPH *digraph = NULL
        cdef CMR_GRAPH_EDGE* forest_arcs = NULL
        cdef CMR_GRAPH_EDGE* coforest_arcs = NULL
        cdef bool* arcs_reversed = NULL
        cdef CMR_SUBMAT* submatrix = NULL
        cdef CMR_NETWORK_STATISTICS stats

        sig_on()
        try:
            if certificate:
                CMR_CALL(CMRtestNetworkMatrix(cmr, self._mat, &result, &digraph, &forest_arcs,
                                              &coforest_arcs, &arcs_reversed, &submatrix, &stats,
                                              time_limit))
            else:
                CMR_CALL(CMRtestNetworkMatrix(cmr, self._mat, &result, NULL, NULL,
                                              NULL, NULL, NULL, &stats, time_limit))
        finally:
            sig_off()

        if not certificate:
            return <bint> result

        if <bint> result:
            sage_digraph = _sage_digraph(digraph, arcs_reversed)
            sage_forest_arcs = tuple(_sage_arc(digraph, forest_arcs[row], arcs_reversed[forest_arcs[row]])
                                      for row in range(self.nrows()))
            sage_coforest_arcs = tuple(_sage_arc(digraph, coforest_arcs[column], arcs_reversed[coforest_arcs[column]])
                                        for column in range(self.ncols()))
            return True, (sage_digraph, sage_forest_arcs, sage_coforest_arcs)

        return False, NotImplemented  # submatrix TBD

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
        Return whether the linear matroid of ``self`` over `\GF{2}` is regular.

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

            sage: MF = matroids.named_matroids.Fano(); MF
            Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
            sage: MFR = MF.representation().change_ring(ZZ); MFR
            [1 0 0 0 1 1 1]
            [0 1 0 1 0 1 1]
            [0 0 1 1 1 0 1]
            sage: MFR2 = block_diagonal_matrix(MFR, MFR, sparse=True); MFR2
            [1 0 0 0 1 1 1|0 0 0 0 0 0 0]
            [0 1 0 1 0 1 1|0 0 0 0 0 0 0]
            [0 0 1 1 1 0 1|0 0 0 0 0 0 0]
            [-------------+-------------]
            [0 0 0 0 0 0 0|1 0 0 0 1 1 1]
            [0 0 0 0 0 0 0|0 1 0 1 0 1 1]
            [0 0 0 0 0 0 0|0 0 1 1 1 0 1]
            sage: MS2 = MFR2.parent(); MS2
            Full MatrixSpace of 6 by 14 sparse matrices over Integer Ring
            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: MFR2cmr = Matrix_cmr_chr_sparse(MS2, MFR2)
            sage: result, certificate = MFR2cmr._is_binary_linear_matroid_regular(
            ....:                           certificate=True)
            sage: result, certificate
            (False, (OneSumNode (6×14) with 2 children, NotImplemented))
            sage: certificate[0].summands()[0].parent_rows_and_columns()
            ((0, 1, 2), (0, 4, 5, 6, 2, 3, 1))
            sage: certificate[0].summands()[1].parent_rows_and_columns()
            ((3, 4, 5), (7, 11, 12, 13, 9, 10, 8))
            sage: unicode_art(certificate[0])  # random (whether the left or the right branch has been followed)
            ╭OneSumNode (6×14) with 2 children╮
            │                                 │
            SeriesParallelReductionNode (3×7) UnknownNode (3×7)
            │
            ThreeConnectedIrregularNode (3×4)
            sage: result, certificate = MFR2cmr._is_binary_linear_matroid_regular(
            ....:                           certificate=True, complete_tree=True)
            sage: result, certificate
            (False, (OneSumNode (6×14) with 2 children, NotImplemented))
            sage: unicode_art(certificate[0])
            ╭OneSumNode (6×14) with 2 children╮
            │                                 │
            SeriesParallelReductionNode (3×7) SeriesParallelReductionNode (3×7)
            │                                 │
            ThreeConnectedIrregularNode (3×4) ThreeConnectedIrregularNode (3×4)

        TESTS:

        This is test ``NestedMinorPivotsTwoSeparation`` in CMR's ``test_regular.cpp``::

            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 11, 11, sparse=True),
            ....:                           [[1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1],
            ....:                            [0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0],
            ....:                            [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
            ....:                            [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
            ....:                            [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:                            [0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:                            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])
            sage: result, certificate = M._is_binary_linear_matroid_regular(
            ....:                           certificate=True, complete_tree=True)
            sage: result, certificate
            (True, GraphicNode (11×11))
            sage: unicode_art(certificate)
            GraphicNode (11×11)
            sage: result, certificate = M._is_binary_linear_matroid_regular(
            ....:                           certificate=True, complete_tree=True,
            ....:                           use_direct_graphicness_test=False)
            sage: result, certificate
            (True, TwoSumNode (11×11) with 2 children)
            sage: unicode_art(certificate)
            ╭──────────TwoSumNode (11×11) with 2 children
            │                 │
            GraphicNode (7×8) SeriesParallelReductionNode (5×4)
                              │
                              GraphicNode (4×4)
        """
        cdef bool result
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
            return <bint> result

        if <bint> result:
            return True, create_DecompositionNode(dec)

        return False, (create_DecompositionNode(dec),
                       NotImplemented)

    def is_totally_unimodular(self, *, time_limit=60.0, certificate=False,
                              use_direct_graphicness_test=True,
                              series_parallel_ok=True,
                              check_graphic_minors_planar=False,
                              complete_tree='if_regular',
                              construct_matrices=False,
                              construct_transposes=False,
                              construct_graphs=False):
        r"""
        Return whether ``self`` is a totally unimodular matrix.

        REFERENCES:

        - [Sch1986]_

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: M.is_totally_unimodular()
            True
            sage: M.is_totally_unimodular(certificate=True)
            (True, GraphicNode (3×2))

            sage: MF = matroids.named_matroids.Fano(); MF
            Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
            sage: MFR = MF.representation().change_ring(ZZ); MFR
            [1 0 0 0 1 1 1]
            [0 1 0 1 0 1 1]
            [0 0 1 1 1 0 1]
            sage: MFR2 = block_diagonal_matrix(MFR, MFR, sparse=True); MFR2
            [1 0 0 0 1 1 1|0 0 0 0 0 0 0]
            [0 1 0 1 0 1 1|0 0 0 0 0 0 0]
            [0 0 1 1 1 0 1|0 0 0 0 0 0 0]
            [-------------+-------------]
            [0 0 0 0 0 0 0|1 0 0 0 1 1 1]
            [0 0 0 0 0 0 0|0 1 0 1 0 1 1]
            [0 0 0 0 0 0 0|0 0 1 1 1 0 1]
            sage: MS2 = MFR2.parent(); MS2
            Full MatrixSpace of 6 by 14 sparse matrices over Integer Ring
            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: MFR2cmr = Matrix_cmr_chr_sparse(MS2, MFR2)
            sage: MFR2cmr.is_totally_unimodular(certificate=True, construct_matrices=True)
            (False, (None, ((0, 1, 2), (3, 4, 5))))
            sage: result, certificate = MFR2cmr.is_totally_unimodular(certificate=True,
            ....:                                                     complete_tree=True,
            ....:                                                     construct_matrices=True)
            sage: result, certificate
            (False, (None, ((0, 1, 2), (3, 4, 5))))
            sage: submatrix = MFR2.matrix_from_rows_and_columns(*certificate[1]); submatrix
            [0 1 1]
            [1 0 1]
            [1 1 0]
            sage: submatrix.determinant()
            2
        """
        cdef bool result
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
            return <bint> result

        if <bint> result:
            return True, create_DecompositionNode(dec)

        if submat == NULL:
            submat_tuple = None
        else:
            submat_tuple = (tuple(submat.rows[i] for i in range(submat.numRows)),
                            tuple(submat.columns[i] for i in range(submat.numColumns)))

        return False, (create_DecompositionNode(dec),
                       submat_tuple)

    def is_complement_totally_unimodular(self, *, time_limit=60.0, certificate=False,
                                         use_direct_graphicness_test=True,
                                         series_parallel_ok=True,
                                         check_graphic_minors_planar=False,
                                         complete_tree='if_regular',
                                         construct_matrices=False,
                                         construct_transposes=False,
                                         construct_graphs=False):
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


cdef _sage_edge(CMR_GRAPH *graph, CMR_GRAPH_EDGE e):
    return Integer(CMRgraphEdgeU(graph, e)), Integer(CMRgraphEdgeV(graph, e))


cdef _sage_graph(CMR_GRAPH *graph):
    # Until we have a proper CMR Graph backend, we just create a Sage graph with whatever backend
    from sage.graphs.graph import Graph

    def vertices():
        i = CMRgraphNodesFirst(graph)
        while CMRgraphNodesValid(graph, i):
            yield i
            i = CMRgraphNodesNext(graph, i)

    def edges():
        i = CMRgraphEdgesFirst(graph)
        while CMRgraphEdgesValid(graph, i):
            e = CMRgraphEdgesEdge(graph, i)
            yield _sage_edge(graph, e)
            i = CMRgraphEdgesNext(graph, i)

    return Graph([list(vertices()), list(edges())])


cdef _sage_arc(CMR_GRAPH *graph, CMR_GRAPH_EDGE e, bint reversed):
    if reversed:
        return Integer(CMRgraphEdgeV(graph, e)), Integer(CMRgraphEdgeU(graph, e))
    return Integer(CMRgraphEdgeU(graph, e)), Integer(CMRgraphEdgeV(graph, e))


cdef _sage_digraph(CMR_GRAPH *graph, bool *arcs_reversed):
    from sage.graphs.digraph import DiGraph

    def vertices():
        i = CMRgraphNodesFirst(graph)
        while CMRgraphNodesValid(graph, i):
            yield i
            i = CMRgraphNodesNext(graph, i)

    def arcs():
        i = CMRgraphEdgesFirst(graph)
        while CMRgraphEdgesValid(graph, i):
            e = CMRgraphEdgesEdge(graph, i)
            yield _sage_arc(graph, e, arcs_reversed[e])
            i = CMRgraphEdgesNext(graph, i)

    return DiGraph([list(vertices()), list(arcs())])
