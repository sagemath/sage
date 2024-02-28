# sage_setup: distribution = sagemath-cmr
r"""
Sparse Matrices with CMR
"""

# ****************************************************************************
#       Copyright (C) 2023 Matthias Koeppe
#                     2023 Luze Xu
#                     2023 Javier Santillan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

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

    EXAMPLES::

        sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_sparse, Matrix_cmr_chr_sparse
        sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
        ....:                           [[1, 2, 3], [4, 0, 6]])
        sage: isinstance(M, Matrix_cmr_sparse)
        True
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

    def matrix_from_rows_and_columns(self, rows, columns):
        """
        Return the matrix constructed from ``self`` from the given rows and
        columns.

        OUTPUT: A :class:`Matrix_cmr_chr_sparse`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = Matrix_cmr_chr_sparse(M, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows_and_columns([1], [0,2])
            [3 5]
            sage: A.matrix_from_rows_and_columns([1,2], [1,2])
            [4 5]
            [7 0]

        Note that row and column indices can be reordered::

            sage: A.matrix_from_rows_and_columns([2,1], [2,1])
            [0 7]
            [5 4]

        But the column indices can not be repeated::

            sage: A.matrix_from_rows_and_columns([1,1,1],[2,0])
            [5 3]
            [5 3]
            [5 3]
            sage: A.matrix_from_rows_and_columns([1,1,1],[2,0,0])
            Traceback (most recent call last):
            ...
            ValueError: The column indices can not be repeated
        """
        if not isinstance(rows, (list, tuple)):
            rows = list(rows)

        if not isinstance(columns, (list, tuple)):
            columns = list(columns)

        if len(list(set(columns))) != len(columns):
            raise ValueError("The column indices can not be repeated")

        if cmr == NULL:
            CMRcreateEnvironment(&cmr)

        cdef CMR_SUBMAT *submatrix = NULL
        cdef CMR_CHRMAT *cmr_submatrix = NULL

        try:
            CMR_CALL(CMRsubmatCreate(cmr, len(rows), len(columns), &submatrix))

            for i in range(submatrix.numRows):
                submatrix.rows[i] = rows[i]

            for j in range(submatrix.numColumns):
                submatrix.columns[j] = columns[j]

            CMR_CALL(CMRchrmatZoomSubmat(cmr, self._mat, submatrix, &cmr_submatrix))
        finally:
            CMR_CALL(CMRsubmatFree(cmr, &submatrix))

        return Matrix_cmr_chr_sparse._from_cmr(cmr_submatrix)

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
            CMR_CALL(CMRchrmatFree(cmr, &self._mat))

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
                CMR_CALL(CMRnetworkComputeMatrix(cmr, cmr_digraph, &cmr_matrix, NULL,
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

        The terminology "1-sum" is used in the context of Seymour's decomposition
        of totally unimodular matrices and regular matroids, see [Sch1986]_.

        .. SEEALSO:: :meth:`two_sum`, :meth:`three_sum`

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

    def two_sum(first_mat, second_mat, first_index, second_index, nonzero_block="top_right"):
        r"""
        Return the 2-sum matrix constructed from the given matrices ``first_mat`` and
        ``second_mat``, with column index of the first matrix ``column`` and row index
        of the second matrix ``row``.
        Suppose that ``column`` indicates the last column and ``row`` indicates the
        first row, i.e, the first matrix is `M_1=\begin{bmatrix} A & a\end{bmatrix}`
        and the second matrix is `M_2=\begin{bmatrix} b^T \\ B\end{bmatrix}`. Then
        the two sum `M_1 \oplus_2 M_2 =\begin{bmatrix}A & ab^T\\ 0 & B\end{bmatrix}`.

        INPUT:

        - ``first_mat`` -- the first integer matrix
        - ``second_mat`` -- the second integer matrix
        - ``first_index`` -- the column/row index of the first integer matrix
        - ``second_index`` -- the row/column index of the second integer matrix
        - ``nonzero_block`` -- either ``"top_right"`` (default) or ``"bottom_left"``;
          whether the nonzero block in the 2-sum matrix locates in the top right or bottom left.
          If ``nonzero_block="top_right"``,
          ``first_index`` is the column index of the first integer matrix,
          ``second_index`` is the row index of the second integer matrix.
          The outer product of the corresponding column and row
          form the nonzero top right block of the 2-sum matrix.
          If ``nonzero_block="bottom_left"``,
          ``first_index`` is the row index of the first integer matrix,
          ``second_index`` is the column index of the second integer matrix.
          The outer product of the corresponding row and column
          form the nonzero bottom left block of the 2-sum matrix.

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
            sage: Matrix_cmr_chr_sparse.two_sum(M1, M2, 2, 0)
            [ 1  2|21 24 27]
            [ 4  5|42 48 54]
            [-----+--------]
            [ 0  0|-1 -2 -3]
            sage: M1.two_sum(M2, 1, 1)
            [  1   3| -2  -4  -6]
            [  4   6| -5 -10 -15]
            [-------+-----------]
            [  0   0|  7   8   9]
            sage: M1.two_sum(M2, 1, 1, nonzero_block="bottom_right")
            Traceback (most recent call last):
            ...
            ValueError: ('Unknown two sum mode', 'bottom_right')
            sage: M1.two_sum(M2, 1, 1, nonzero_block="bottom_left")
            [  1   2   3|  0   0]
            [-----------+-------]
            [ 32  40  48|  7   9]
            [ -8 -10 -12| -1  -3]
        """
        cdef Matrix_cmr_chr_sparse sum, first, second
        cdef CMR_CHRMAT *sum_mat
        first = Matrix_cmr_chr_sparse._from_data(first_mat, immutable=False)
        second = Matrix_cmr_chr_sparse._from_data(second_mat, immutable=False)

        if nonzero_block not in ["top_right", "bottom_left"]:
            raise ValueError("Unknown two sum mode", nonzero_block)

        if nonzero_block == "top_right":
            column = first_index
            row = second_index
            if column < 0 or column >= first._mat.numColumns:
                raise ValueError("First marker should be a column index of the first matrix")
            if row < 0 or row >= second._mat.numRows:
                raise ValueError("Second marker should be a row index of the second matrix")
            row_subdivision = []
            column_subdivision = []
            row_subdivision.append(first._mat.numRows)
            column_subdivision.append(first._mat.numColumns - 1)
            first_marker = CMRcolumnToElement(column)
            second_marker = CMRrowToElement(row)
        else:
            row = first_index
            column = second_index
            if row < 0 or row >= first._mat.numRows:
                raise ValueError("First marker should be a Row index of the first matrix")
            if column < 0 or column >= second._mat.numColumns:
                raise ValueError("Second marker should be a column index of the second matrix")
            row_subdivision = []
            column_subdivision = []
            row_subdivision.append(first._mat.numRows - 1)
            column_subdivision.append(first._mat.numColumns)
            first_marker = CMRrowToElement(row)
            second_marker = CMRcolumnToElement(column)

        cdef int8_t characteristic = 0
        CMR_CALL(CMRtwoSum(cmr, first._mat, second._mat, first_marker, second_marker, characteristic, &sum_mat))
        sum = Matrix_cmr_chr_sparse._from_cmr(sum_mat, immutable=False)
        if row_subdivision or column_subdivision:
            sum.subdivide(row_subdivision, column_subdivision)
        sum.set_immutable()
        return sum

    def three_sum(first_mat, second_mat, first_col_index1, first_col_index2, second_col_index1, second_col_index2):
        r"""
        Return the 3-sum matrix constructed from the given matrices ``first_mat`` and
        ``second_mat``, with ``first_col_index1`` and ``first_col_index2`` being the
        indices of the column vectors of the matrix, which are identical except for
        one row having a 0 in one column and the other a non-zero entry in that row.
        The method assumes the nonzero entry is one. The same goes are made for ``second_mat``, ``second_col_index1``, and  ``second_col_index2``.

        The operation is defined in [Sch1986]_, Ch. 19.4.:=
        [first_mat  first_col  first_col]  ___|___   [second_mat  second_col  second_col]
        [first_row      0          1    ]     |   3  [second_row       0           1    ]
        -----  [         first_mat       first_col x second_row]
        -----  [second_col x first_row          second_col     ]

        INPUT:

        - ``first_mat`` -- integer matrix having two collumns which are identical in
          every entry except for one row in which one is 0 and the other is 1
        - ``second_mat`` -- integer matrix having two collumns which are identical in
          every entry except for one row in which one is 0 and the other is 1
        - ``first_col_index1`` -- index of a column in ``first_mat`` identical to some
          other column in every entry except for one row in which one is 0 and the other is 1
        - ``first_col_index2`` -- index of the other column which is identical to
          first_mat[first_col_index1] in every entry except for one row in which
          one is 0 and the other is 1
        - ``second_col_index1`` -- index of a column in ``second_mat`` identical to some
          other column in every entry except for one row in which one is 0 and the other is 1
        - ``first_col_index2`` -- index of the other column which is identical to
          second_mat[second_col_index1] in every entry except for one row in which
          one is 0 and the other is 1

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
            ....:                           [[1, 1, 1, 1, 1], [1, 1, 1, 0, 0],
            ....:                            [1, 0, 1, 1, 0], [0, 0, 0, 1, 1],
            ....:                            [1, 1, 0, 0, 1]]); M
            [1 1 1 1 1]
            [1 1 1 0 0]
            [1 0 1 1 0]
            [0 0 0 1 1]
            [1 1 0 0 1]
            sage: M3 = Matrix_cmr_chr_sparse.three_sum(M, M, 0, 1, 0, 1); M3
            [1 1 1 1 1 0]
            [1 0 0 1 1 0]
            [0 1 1 0 0 0]
            [0 0 1 1 1 0]
            [1 1 0 1 1 1]
            [1 1 0 1 0 0]
            [0 0 0 0 1 1]
            [1 1 0 0 0 1]
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
                    if i == fr:
                        first_row_index = i - 1
                    else:
                        first_row_index = i
                    break
        for i in sec_nrows:
            if (second_col1[i] == 1 and second_col2[i] == 0) or (second_col1[i] == 0 and second_col2[i] == 1):
                subcol1 = tuple(second_col1[k] for k in sec_nrows if k != i)
                subcol2 = tuple(second_col2[k] for k in sec_nrows if k != i)
                if subcol1 == subcol2:
                    valid2 = True
                    if i == sr:
                        second_row_index = i - 1
                    else:
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
        return Matrix_cmr_chr_sparse._from_data(row_list, immutable=False)

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
        return Matrix_cmr_chr_sparse._from_data(row_list, immutable=False)

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
        return Matrix_cmr_chr_sparse._from_data(row_list, immutable=False)

    def binary_pivot(self, row, column):
        cdef Matrix_cmr_chr_sparse result
        cdef size_t pivot_row = row
        cdef size_t pivot_column = column
        cdef CMR_CHRMAT *result_mat

        CMR_CALL(CMRchrmatBinaryPivot(cmr, self._mat, pivot_row, pivot_column, &result_mat))
        result = Matrix_cmr_chr_sparse._from_cmr(result_mat)
        return result

    def binary_pivots(self, rows, columns):
        cdef Matrix_cmr_chr_sparse result
        cdef size_t* pivot_rows
        cdef size_t* pivot_columns
        cdef CMR_CHRMAT *result_mat

        npivots = len(rows)
        if len(columns) != npivots:
            raise ValueError("The pivot rows and columns must have the same length")

        for i in range(npivots):
            pivot_rows[i] = rows[i]
            pivot_columns[i] = columns[i]

        CMR_CALL(CMRchrmatBinaryPivots(cmr, self._mat, npivots, pivot_rows, pivot_columns, &result_mat))
        result = Matrix_cmr_chr_sparse._from_cmr(result_mat)
        return result

    def ternary_pivot(self, row, column):
        cdef Matrix_cmr_chr_sparse result
        cdef size_t pivot_row = row
        cdef size_t pivot_column = column
        cdef CMR_CHRMAT *result_mat

        CMR_CALL(CMRchrmatTernaryPivot(cmr, self._mat, pivot_row, pivot_column, &result_mat))
        result = Matrix_cmr_chr_sparse._from_cmr(result_mat)
        return result

    def ternary_pivots(self, rows, columns):
        cdef Matrix_cmr_chr_sparse result
        cdef size_t* pivot_rows
        cdef size_t* pivot_columns
        cdef CMR_CHRMAT *result_mat

        cdef size_t npivots = len(rows)
        if len(columns) != npivots:
            raise ValueError("The pivot rows and columns must have the same length")

        for i in range(npivots):
            pivot_rows[i] = rows[i]
            pivot_columns[i] = columns[i]

        CMR_CALL(CMRchrmatTernaryPivots(cmr, self._mat, npivots, pivot_rows, pivot_columns, &result_mat))
        result = Matrix_cmr_chr_sparse._from_cmr(result_mat)
        return result

    def is_unimodular(self, time_limit=60.0):
        r"""
        Return whether ``self`` is a unimodular matrix.

        A nonsingular square matrix `A` is called unimodular if it is integral
        and has determinant `\pm1`, i.e., an element of
        `\mathop{\operatorname{GL}}_n(\ZZ)` [Sch1986]_, Ch. 4.3.

        A rectangular matrix `A` of full row rank is called unimodular if it
        is integral and every basis `B` of `A` has determinant `\pm1`.
        [Sch1986]_, Ch. 19.1.

        More generally, a matrix `A` of rank `r` is called unimodular if it is
        integral and for every submatrix `B` formed by `r` linearly independent columns,
        the greatest common divisor of the determinants of all `r`-by-`r`
        submatrices of `B` is `1`. [Sch1986]_, Ch. 21.4.

        .. SEEALSO:: :meth:`is_k_modular`, :meth:`is_strongly_unimodular`, :meth:`is_totally_unimodular`

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
        cdef CMR_INTMAT *int_mat = NULL
        cdef bool result

        sig_on()
        try:
            CMR_CALL(CMRchrmatToInt(cmr, self._mat, &int_mat))
            CMR_CALL(CMRunimodularTest(cmr, int_mat, &result, NULL, NULL, time_limit))
        finally:
            CMR_CALL(CMRintmatFree(cmr, &int_mat))
            sig_off()

        return <bint> result

    def is_strongly_unimodular(self, time_limit=60.0):
        r"""
        Return whether ``self`` is a strongly unimodular matrix.

        A matrix is strongly unimodular if ``self`` and ``self.transpose()`` are both unimodular.

        .. SEEALSO:: meth:`is_unimodular`, :meth:`is_strongly_k_modular`

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
        cdef CMR_INTMAT *int_mat = NULL
        cdef bool result

        sig_on()
        try:
            CMR_CALL(CMRchrmatToInt(cmr, self._mat, &int_mat))
            CMR_CALL(CMRunimodularTestStrong(cmr, int_mat, &result, NULL, NULL, time_limit))
        finally:
            CMR_CALL(CMRintmatFree(cmr, &int_mat))
            sig_off()

        return <bint> result

    def equimodulus(self, time_limit=60.0):
        r"""
        Return the integer `k` such that ``self`` is
        equimodular with determinant gcd `k`.

        A matrix `M` of rank `r` is equimodular with determinant gcd `k`
        if the following two conditions are satisfied:
        - for some column basis `B` of `M`, the greatest common divisor of
          the determinants of all `r`-by-`r` submatrices of `B` is `k`;
        - the matrix `X` such that `M=BX` is totally unimodular.

        OUTPUT:

        - ``k``: ``self`` is equimodular with determinant gcd `k`
        - ``None``: ``self`` is not equimodular for any `k`

        .. SEEALSO:: :meth:`is_k_equimodular`, :meth:`strong_equimodulus`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.equimodulus()
            1
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 1, 1], [0, 1, 3]]); M
            [1 1 1]
            [0 1 3]
            sage: M.equimodulus()
        """
        cdef CMR_INTMAT *int_mat = NULL
        cdef bool result
        cdef int64_t k = 0

        sig_on()
        try:
            CMR_CALL(CMRchrmatToInt(cmr, self._mat, &int_mat))
            CMR_CALL(CMRequimodularTest(cmr, int_mat, &result, &k, NULL, NULL, time_limit))
        finally:
            CMR_CALL(CMRintmatFree(cmr, &int_mat))
            sig_off()

        if result:
            return k
        else:
            return None

    def strong_equimodulus(self, time_limit=60.0):
        r"""
        Return the integer `k` such that ``self`` is
        strongly equimodular with determinant gcd `k`.

        Return whether ``self`` is strongly `k`-equimodular.

        A matrix is strongly equimodular if ``self`` and ``self.transpose()``
        are both equimodular, which implies that they are equimodular for
        the same determinant gcd `k`.
        A matrix `M` of rank-`r` is `k`-modular if the following two conditions
        are satisfied:
        - for some column basis `B` of `M`, the greatest common divisor of the
          determinants of all `r`-by-`r` submatrices of `B` is `k`;
        - the matrix `X` such that `M=BX` is totally unimodular.

        OUTPUT:

        - ``k``: ``self`` is  `k`-equimodular
        - ``None``: ``self`` is not `k`-equimodular for any `k`

        .. SEEALSO:: :meth:`is_strongly_k_equimodular`, :meth:`equimodulus`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.strong_equimodulus()
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 0, 0], [0, 1, 0]]); M
            [1 0 0]
            [0 1 0]
            sage: M.strong_equimodulus()
            1
        """
        cdef CMR_INTMAT *int_mat = NULL
        cdef bool result
        cdef int64_t k = 0

        sig_on()
        try:
            CMR_CALL(CMRchrmatToInt(cmr, self._mat, &int_mat))
            CMR_CALL(CMRequimodularTestStrong(cmr, int_mat, &result, &k, NULL, NULL, time_limit))
        finally:
            CMR_CALL(CMRintmatFree(cmr, &int_mat))
            sig_off()

        if result:
            return k
        else:
            return None

    def is_k_equimodular(self, k, time_limit=60.0):
        r"""
        Return whether ``self`` is equimodular with determinant gcd `k`.

        A matrix `M` of rank-`r` is `k`-equimodular if the following two
        conditions are satisfied:
        - for some column basis `B` of `M`, the greatest common divisor of
          the determinants of all `r`-by-`r` submatrices of `B` is `k`;
        - the matrix `X` such that `M=BX` is totally unimodular.

        If the matrix has full row rank, it is `k`-equimodular if
        every full rank minor of the matrix has determinant `0,\pm k`.

        .. NOTE::

            In parts of the literature, a matrix with the above properties
            is called *strictly* `k`-modular.

        .. SEEALSO:: :meth:`is_unimodular`, :meth:`is_strongly_k_equimodular`,
                     :meth:`equimodulus`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.is_k_equimodular(1)
            True
            sage: M.is_k_equimodular(2)
            False
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 1, 1], [0, 1, 3]]); M
            [1 1 1]
            [0 1 3]
            sage: M.is_k_equimodular(1)
            False
        """
        cdef CMR_INTMAT *int_mat = NULL
        cdef bool result
        cdef int64_t gcd_det = k

        sig_on()
        try:
            CMR_CALL(CMRchrmatToInt(cmr, self._mat, &int_mat))
            CMR_CALL(CMRequimodularTest(cmr, int_mat, &result, &gcd_det, NULL, NULL, time_limit))
        finally:
            CMR_CALL(CMRintmatFree(cmr, &int_mat))
            sig_off()

        return True if result else False

    def is_strongly_k_equimodular(self, k, time_limit=60.0):
        r"""
        Return whether ``self`` is strongly `k`-equimodular.

        A matrix is strongly `k`-equimodular if ``self`` and ``self.transpose()``
        are both `k`-equimodular.

        .. SEEALSO:: :meth:`is_k_equimodular`, :meth:`is_strongly_unimodular`,
                     :meth:`strong_equimodulus`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 3, sparse=True),
            ....:                           [[1, 0, 1], [0, 1, 1], [1, 2, 3]]); M
            [1 0 1]
            [0 1 1]
            [1 2 3]
            sage: M.is_strongly_k_equimodular(1)
            False
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 0, 0], [0, 1, 0]]); M
            [1 0 0]
            [0 1 0]
            sage: M.is_strongly_k_equimodular(1)
            True
        """
        cdef CMR_INTMAT *int_mat = NULL
        cdef bool result
        cdef int64_t gcd_det = k

        sig_on()
        try:
            CMR_CALL(CMRchrmatToInt(cmr, self._mat, &int_mat))
            CMR_CALL(CMRequimodularTestStrong(cmr, int_mat, &result, &gcd_det, NULL, NULL, time_limit))
        finally:
            CMR_CALL(CMRintmatFree(cmr, &int_mat))
            sig_off()

        return True if result else False

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
                CMR_CALL(CMRgraphicTestMatrix(cmr, self._mat, &result, &graph, &forest_edges,
                                              &coforest_edges, &submatrix, &stats, time_limit))
            else:
                CMR_CALL(CMRgraphicTestMatrix(cmr, self._mat, &result, NULL, NULL,
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
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 4, 9, sparse=True),
            ....:                           [[1, 0, 0, 0, 1, -1, 1, 0, 0],
            ....:                            [0, 1, 0, 0, 0, 1, -1, 1, 0],
            ....:                            [0, 0, 1, 0, 0, 0, 1, -1, 1],
            ....:                            [0, 0, 0, 1, 1, 0, 0, 1, -1]]); M
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
                CMR_CALL(CMRgraphicTestTranspose(cmr, self._mat, &result, &graph, &forest_edges,
                                              &coforest_edges, &submatrix, &stats, time_limit))
            else:
                CMR_CALL(CMRgraphicTestTranspose(cmr, self._mat, &result, NULL, NULL,
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
        cdef bool support_result
        cdef CMR_GRAPH *digraph = NULL
        cdef CMR_GRAPH_EDGE* forest_arcs = NULL
        cdef CMR_GRAPH_EDGE* coforest_arcs = NULL
        cdef bool* arcs_reversed = NULL
        cdef CMR_SUBMAT* submatrix = NULL
        cdef CMR_NETWORK_STATISTICS stats

        sig_on()
        try:
            if certificate:
                CMR_CALL(CMRnetworkTestMatrix(cmr, self._mat, &result, &support_result, &digraph, &forest_arcs,
                                              &coforest_arcs, &arcs_reversed, &submatrix, &stats,
                                              time_limit))
            else:
                CMR_CALL(CMRnetworkTestMatrix(cmr, self._mat, &result, &support_result, NULL, NULL,
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
                                          three_sum_pivot_children=False,
                                          three_sum_strategy=None,
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

            sage: MF = matroids.catalog.Fano(); MF
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
        cdef CMR_REGULAR_PARAMS params
        cdef CMR_REGULAR_STATS stats
        cdef CMR_MATROID_DEC *dec = NULL
        cdef CMR_MINOR *minor = NULL

        cdef CMR_MATROID_DEC **pdec = &dec
        cdef CMR_MINOR **pminor = &minor

        cdef dict kwds = dict(use_direct_graphicness_test=use_direct_graphicness_test,
                              series_parallel_ok=series_parallel_ok,
                              check_graphic_minors_planar=check_graphic_minors_planar,
                              complete_tree=complete_tree,
                              three_sum_pivot_children=three_sum_pivot_children,
                              three_sum_strategy=three_sum_strategy,
                              construct_graphs=construct_graphs)

        _set_cmr_regular_parameters(&params, kwds)
        sig_on()
        try:
            CMR_CALL(CMRregularTest(cmr, self._mat, &result, pdec, pminor,
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
                              three_sum_pivot_children=False,
                              three_sum_strategy=None,
                              construct_graphs=False):
        r"""
        Return whether ``self`` is a totally unimodular matrix.

        A matrix is totally unimodular if every subdeterminant is `0`, `1`, or `-1`.

        REFERENCES:

        - [Sch1986]_, Chapter 19

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

            sage: MF = matroids.catalog.Fano(); MF
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
            sage: MFR2cmr.is_totally_unimodular(certificate=True)
            (False, (None, ((0, 1, 2), (3, 4, 5))))
            sage: result, certificate = MFR2cmr.is_totally_unimodular(certificate=True,
            ....:                                                     complete_tree=True)
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
        cdef CMR_TU_PARAMS params
        cdef CMR_TU_STATS stats
        cdef CMR_MATROID_DEC *dec = NULL
        cdef CMR_SUBMAT *submat = NULL

        cdef CMR_MATROID_DEC **pdec = &dec
        cdef CMR_SUBMAT **psubmat = &submat

        if three_sum_pivot_children:
            raise NotImplementedError
        cdef dict kwds = dict(use_direct_graphicness_test=use_direct_graphicness_test,
                              series_parallel_ok=series_parallel_ok,
                              check_graphic_minors_planar=check_graphic_minors_planar,
                              complete_tree=complete_tree,
                              three_sum_pivot_children=three_sum_pivot_children,
                              three_sum_strategy=three_sum_strategy,
                              construct_graphs=construct_graphs)

        params.algorithm = CMR_TU_ALGORITHM_DECOMPOSITION
        _set_cmr_regular_parameters(&params.regular, kwds)
        sig_on()
        try:
            CMR_CALL(CMRtuTest(cmr, self._mat, &result, pdec, psubmat,
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


cdef _set_cmr_regular_parameters(CMR_REGULAR_PARAMS *params, dict kwds):
    CMR_CALL(CMRregularParamsInit(params))
    params.directGraphicness = kwds['use_direct_graphicness_test']
    params.seriesParallel = kwds['series_parallel_ok']
    params.planarityCheck = kwds['check_graphic_minors_planar']
    params.completeTree = kwds['complete_tree'] is True
    params.threeSumPivotChildren = kwds['three_sum_pivot_children']
    if kwds['three_sum_strategy'] is not None:
        params.threeSumStrategy = kwds['three_sum_strategy']
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
