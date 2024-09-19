# sage_setup: distribution = sagemath-cmr
r"""
Seymour's decomposition of totally unimodular matrices and regular matroids
"""

# ****************************************************************************
#       Copyright (C) 2023      Javier Santillan
#                     2023-2024 Matthias Koeppe
#                     2023-2024 Luze Xu
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
from sage.misc.cachefunc import cached_method
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.structure.sage_object cimport SageObject

from .constructor import Matrix
from .matrix_cmr_sparse cimport Matrix_cmr_chr_sparse, _sage_edges, _sage_graph, _set_cmr_seymour_parameters
from .matrix_space import MatrixSpace


cdef class DecompositionNode(SageObject):
    r"""
    Base class for nodes in Seymour's decomposition
    """

    def __cinit__(self, *args, **kwds):
        r"""
        Initialize the internal decomposition, a ``CMR_SEYMOUR_NODE``.
        """
        self._dec = NULL

    def __init__(self, matrix=None, row_keys=None, column_keys=None, base_ring=None):
        r"""
        Create a node in Seymour's decomposition.

        INPUT:

        - ``matrix`` -- the internal matrix representing the node.
        Convert to a :class:`Matrix_cmr_chr_sparse`.

        - ``row_keys`` -- a finite or enumerated family of arbitrary objects
        that index the rows of the matrix

        - ``column_keys`` -- a finite or enumerated family of arbitrary objects
        that index the columns of the matrix

        - ``base_ring`` -- the base ring of ``matrix`` representing the node.
        For Seymour decomposition node, the base ring is `\GF{2}` or `\GF{3}` or ``ZZ``.
        If the base ring is `\GF{2}`, the node and the matrix deal with
        the underlying binary linear matroid;
        if the base ring is `\GF(3)` or ``ZZ``, the node deals with
        the matrix decomposition.

        A :class:DecompositionNode is usually created with an internal decomposition,
        ``self._dec``, see :meth:create_DecompositionNode
        Such decomposition comes from the certificate of the
        totally unimodularity test, see
        :meth:`matrix_cmr_sparse.Matrix_cmr_chr_sparse.is_totally_unimodular`

        Another usage is to create a :class:UnknownNode from a matrix.
        A root dummy decomposition is created before completing
        the decomposition, see :meth:_set_root_dec, :meth:complete_decomposition

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, 1], [0, 1]]); node
            UnknownNode (2×2)
        """
        if matrix is None:
            self._matrix = None
        elif isinstance(matrix, Matrix_cmr_chr_sparse):
            self._matrix = matrix
        else:
            try:
                self._matrix = matrix._matrix_cmr()
            except (AttributeError, ImportError, TypeError):
                if base_ring is not None:
                    matrix = Matrix(matrix, ring=base_ring)
                else:
                    matrix = Matrix(matrix)
                self._matrix = Matrix_cmr_chr_sparse(matrix.parent(), matrix)
            else:
                if row_keys is None:
                    row_keys = matrix.codomain().basis().keys()
                if column_keys is None:
                    column_keys = matrix.domain().basis().keys()
        if row_keys is not None:
            self._set_row_keys(row_keys)
        if column_keys is not None:
            self._set_column_keys(column_keys)
        if base_ring is None:
            if self._matrix is not None:
                base_ring = self._matrix.parent().base_ring()
        self._base_ring = base_ring

    cdef _set_dec(self, CMR_SEYMOUR_NODE *dec):
        r"""
        Set the decomposition ``self._dec`` to ``dec``.
        If the value was previously set, then free it first.
        """
        if self._dec != NULL:
            # We own it, so we have to free it.
            CMR_CALL(CMRseymourRelease(cmr, &self._dec))
        if dec != NULL:
            CMR_CALL(CMRseymourCapture(cmr, dec))
        self._dec = dec

    cdef _set_root_dec(self):
        r"""
        Set the decomposition by creating a root ``CMR_SEYMOUR_NODE``
        based on the internal matrix representation ``self._matrix``.
        """
        cdef CMR_SEYMOUR_NODE *root
        cdef Matrix_cmr_chr_sparse matrix
        try:
            matrix = self.matrix()
        except:
            raise ValueError('no Matrix_cmr_chr_sparse matrix')
        base_ring = self.base_ring()
        if base_ring.characteristic() not in [0, 2, 3] :
            raise ValueError(f'only defined over binary or ternary, got {base_ring}')
        isTernary = base_ring.characteristic() != 2
        cdef CMR_CHRMAT *mat = matrix._mat

        sig_on()
        try:
            CMR_CALL(CMRseymourCreate(cmr, &root, isTernary, mat))
        finally:
            sig_off()
        self._set_dec(root)

    cdef _set_row_keys(self, row_keys):
        """
        Set the row keys with consistency checking: if the
        value was previously set, it must remain the same.
        """
        if row_keys is not None:
            row_keys = tuple(row_keys)
        if self._row_keys is not None and self._row_keys != row_keys:
            raise ValueError(f"inconsistent row keys: should be {self._row_keys} "
                             f"but got {row_keys}")
        if row_keys is not None and self._dec != NULL and self.nrows() != len(row_keys):
            raise ValueError(f"inconsistent row keys: should be of cardinality {self.nrows()} "
                             f"but got {row_keys}")
        self._row_keys = row_keys

    cdef _set_column_keys(self, column_keys):
        """
        Set the column keys with consistency checking: if the
        value was previously set, it must remain the same.
        """
        if column_keys is not None:
            column_keys = tuple(column_keys)
        if self._column_keys is not None and self._column_keys != column_keys:
            raise ValueError(f"inconsistent column keys: should be {self._column_keys} "
                             f"but got {column_keys}")
        if column_keys is not None and self._dec != NULL and self.ncols() != len(column_keys):
            raise ValueError(f"inconsistent column keys: should be of cardinality {self.ncols()} "
                             f"but got {column_keys}")
        self._column_keys = column_keys

    def __dealloc__(self):
        """
        Frees all the memory allocated for this node.
        """
        self._set_dec(NULL)

    def __hash__(self):
        """
        Return a hash of this node. It is the hash of the decomposition.
        """
        return <int>self._dec

    def nrows(self):
        r"""
        Return the number of rows of the internal matrix representing this node.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.nrows()
            3
            sage: certificate.ncols()
            2

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, 1], [0, 1]]); node
            UnknownNode (2×2)
            sage: node.nrows()
            2
            sage: node.ncols()
            2
        """
        if self._row_keys is not None:
            return len(self._row_keys)
        if self._dec != NULL:
            return CMRseymourNumRows(self._dec)
        if self._matrix is not None:
            return self._matrix.nrows()
        raise RuntimeError('nrows undefined')

    def ncols(self):
        r"""
        Return the number of columns of the internal matrix representing this node.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.nrows()
            3
            sage: certificate.ncols()
            2

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, 1], [0, 1]]); node
            UnknownNode (2×2)
            sage: node.nrows()
            2
            sage: node.ncols()
            2
        """
        if self._column_keys is not None:
            return len(self._column_keys)
        if self._dec != NULL:
            return CMRseymourNumColumns(self._dec)
        if self._matrix is not None:
            return self._matrix.ncols()
        raise RuntimeError('ncols undefined')

    def dimensions(self):
        r"""
        Return the number of rows and columns of this node.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.dimensions()
            (3, 2)

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, 1], [0, 1]]); node
            UnknownNode (2×2)
            sage: node.dimensions()
            (2, 2)
        """
        return self.nrows(), self.ncols()

    def base_ring(self):
        r"""
        Return the base ring of the matrix representing the node.
        """
        return self._base_ring

    def matrix(self):
        r"""
        Return a :class:`Matrix`.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.matrix()
            [ 1  0]
            [-1  1]
            [ 0  1]

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, 1], [0, 1]]); node
            UnknownNode (2×2)
            sage: node.matrix()
            [1 1]
            [0 1]
        """
        if self._matrix is not None:
            return self._matrix
        if self._dec is NULL:
            if isinstance(self, SumNode):
                return self.block_matrix_form()
            raise ValueError('Matrix and decomposition are both missing')
        cdef Matrix_cmr_chr_sparse result
        cdef CMR_CHRMAT *mat = CMRseymourGetMatrix(self._dec)
        if mat == NULL:
            return None
        ms = MatrixSpace(self.base_ring(), mat.numRows, mat.numColumns, sparse=True)
        result = Matrix_cmr_chr_sparse.__new__(Matrix_cmr_chr_sparse, ms)
        result._mat = mat
        result._root = self  # Matrix is owned by us
        self._matrix = result
        return result

    def row_keys(self):
        r"""
        Return the row keys of this node.

        OUTPUT: a tuple or ``None``

        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode(matrix(ZZ, [[1, 0, 1], [0, 1, 1]]),
            ....:                    row_keys='ab',
            ....:                    column_keys=range(3)); node
            UnknownNode (2×3)
            sage: node.row_keys()
            ('a', 'b')
        """
        return self._row_keys

    def column_keys(self):
        r"""
        Return the column keys of this node.

        OUTPUT: a tuple or ``None``

        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode(matrix(ZZ, [[1, 0, 1], [0, 1, 1]]),
            ....:                    row_keys='ab',
            ....:                    column_keys=range(3)); node
            UnknownNode (2×3)
            sage: node.column_keys()
            (0, 1, 2)
        """
        return self._column_keys

    def set_default_keys(self):
        r"""
        Set default row and column keys.

        .. SEEALSO:: :class:ElementKey

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.row_keys() is None
            True
            sage: certificate.set_default_keys()
            sage: certificate.row_keys()
            (r0, r1, r2)
        """
        row_keys = self.row_keys()
        column_keys = self.column_keys()
        if row_keys is None or column_keys is None:
            row_keys = tuple(ElementKey(f"r{i}") for i in range(self.nrows()))
            column_keys  = tuple(ElementKey(f"c{i}") for i in range(self.ncols()))
        elif not isinstance(row_keys[0], ElementKey):
            row_keys = tuple(ElementKey(key) for key in row_keys)
            column_keys  = tuple(ElementKey(key) for key in column_keys)
        self._row_keys = row_keys
        self._column_keys = column_keys

    @cached_method
    def morphism(self):
        r"""
        Create the matrix in the distinguished bases of the domain and codomain
        to build the module morphism.

        See also: :class:`sage.modules.with_basis.morphism.ModuleMorphismFromMatrix`.

        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode(matrix(ZZ, [[1, 0, 1], [0, 1, 1]]),
            ....:                    row_keys='ab',
            ....:                    column_keys=range(3)); node
            UnknownNode (2×3)
            sage: node.matrix()
            [1 0 1]
            [0 1 1]
            sage: node.morphism()._unicode_art_matrix()
            0 1 2
            a⎛1 0 1⎞
            b⎝0 1 1⎠
        """
        return Matrix(self.matrix(),
                      row_keys=self.row_keys(),
                      column_keys=self.column_keys())

    def as_ordered_tree(self):
        r"""
        Return the decomposition tree rooted at ``self``.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = matrix([[1, 0], [-1, 1], [0, 1]], sparse=True)
            sage: M2 = block_diagonal_matrix([M, M], sparse=True)
            sage: M2cmr = Matrix_cmr_chr_sparse(M2.parent(), M2); M2cmr
            [ 1  0  0  0]
            [-1  1  0  0]
            [ 0  1  0  0]
            [ 0  0  1  0]
            [ 0  0 -1  1]
            [ 0  0  0  1]
            sage: result, certificate = M2cmr.is_totally_unimodular(certificate=True)
            sage: T = certificate.as_ordered_tree(); T
            OneSumNode (6×4) with 2 children[GraphicNode (3×2)[], GraphicNode (3×2)[]]
            sage: unicode_art(T)
            ╭───────────OneSumNode (6×4) with 2 children
            │                 │
            GraphicNode (3×2) GraphicNode (3×2)
        """
        from sage.combinat.ordered_tree import LabelledOrderedTree
        return LabelledOrderedTree([child.as_ordered_tree() for child in self.child_nodes()],
                                   label=self)

    def plot(self, **kwds):
        r"""
        Plot the decomposition tree rooted at ``self``.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = matrix([[1, 0], [-1, 1], [0, 1]], sparse=True)
            sage: M2MT = block_diagonal_matrix([M, M, M.T], sparse=True)
            sage: M2MTcmr = Matrix_cmr_chr_sparse(M2MT.parent(), M2MT)
            sage: result, certificate = M2MTcmr.is_totally_unimodular(certificate=True)
            sage: T = certificate.as_ordered_tree()
            sage: T.plot()                                                              # needs sage.plot
            Graphics object consisting of 8 graphics primitives
        """
        return self.as_ordered_tree().plot(**kwds)

    def is_ternary(self):
        r"""
        Returns true iff the decomposition is over `\mathbb{F}_3`.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.is_ternary()
            True

            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [1, 1], [0, 1]]); M
            [1  0]
            [1  1]
            [0  1]
            sage: result, certificate = M._is_binary_linear_matroid_regular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.is_ternary()
            False
        """
        return <bint> CMRseymourIsTernary(self._dec)

    def nchildren(self):
        r"""
        Return the number of children of the node.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.nchildren()
            0
        """
        if self._child_nodes is not None:
            return len(self._child_nodes)
        if self._dec == NULL:
            return 0
        return CMRseymourNumChildren(self._dec)

    cdef _CMRelement_to_key(self, CMR_ELEMENT element):
        r"""
        Transform a ``CMRelement`` (row or column index implemented in cmr)
        to a row key or a column key.
        """
        if not CMRelementIsValid(element):
            raise ValueError('CMRelement index not valid. Extra row or column is detected.')
        if self.row_keys() is None or self.column_keys() is None:
            raise ValueError('row_keys and column_keys are required')
        if CMRelementIsRow(element):
            return self.row_keys()[CMRelementToRowIndex(element)]
        else:
            return self.column_keys()[CMRelementToColumnIndex(element)]

    def _create_child_node(self, index):
        r"""
        Return the child node of ``self`` corresponding to the ``index``,
        and the corresponding row and column keys in the parent node.

        OUTPUT: a tuple of (child node, child row keys, child column keys)
        """
        row_keys = self.row_keys()
        column_keys = self.column_keys()
        cdef CMR_SEYMOUR_NODE *child_dec = CMRseymourChild(self._dec, index)
        cdef CMR_ELEMENT *parent_rows = CMRseymourChildRowsToParent(self._dec, index)
        cdef CMR_ELEMENT *parent_columns = CMRseymourChildColumnsToParent(self._dec, index)
        child_nrows = CMRseymourNumRows(child_dec)
        child_ncols = CMRseymourNumColumns(child_dec)

        if parent_rows == NULL or all(parent_rows[i] == 0 for i in range(child_nrows)):
            raise ValueError(f"Child {index} does not have parents rows")
        parent_rows_tuple = tuple(parent_rows[i] for i in range(child_nrows))

        if parent_columns == NULL or all(parent_columns[i] == 0 for i in range(child_ncols)):
            raise ValueError(f"Child {index} does not have parents columns")
        parent_columns_tuple = tuple(parent_columns[i] for i in range(child_ncols))

        if row_keys is not None and column_keys is not None:
            child_row_keys = tuple(self._CMRelement_to_key(element)
                                   for element in parent_rows_tuple)
            child_column_keys = tuple(self._CMRelement_to_key(element)
                                      for element in parent_columns_tuple)
            child = create_DecompositionNode(child_dec, matrix=None,
                                             row_keys=child_row_keys,
                                             column_keys=child_column_keys,
                                             base_ring=self.base_ring())
        else:
            child_row_keys = tuple(CMRelementToRowIndex(element)
                                      for element in parent_rows_tuple)
            child_column_keys = tuple(CMRelementToColumnIndex(element)
                                         for element in parent_columns_tuple)
            child = create_DecompositionNode(child_dec, matrix=None,
                                             row_keys=child_row_keys,
                                             column_keys=child_column_keys,
                                             base_ring=self.base_ring())
        return child, child_row_keys, child_column_keys

    def _children(self):
        r"""
        Return a tuple of the tuples of children and their row and column keys.
        The underlying implementation of :meth:`child_nodes`
        and :meth:`child_indices`.
        """
        if self._child_nodes is not None:
            return self._child_nodes
        children_tuple = tuple(self._create_child_node(index)
                               for index in range(self.nchildren()))
        self._child_nodes = children_tuple
        return self._child_nodes

    def child_nodes(self):
        r"""
        Return a tuple of the children.

        The children are sorted by the inherited ordering from cmr, which
        is their appreance in the parent.

        In the case of :class:`SumNode`, this is the same as :meth:`~SumNode.summands`.

        For graphic or leaf nodes, it returns the empty tuple.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]],
            ....:                                   [[1, 1], [-1, 0]],
            ....:                                   [[1, 0], [0,1]]); M
            [ 1  0| 0  0| 0  0]
            [-1  1| 0  0| 0  0]
            [-----+-----+-----]
            [ 0  0| 1  1| 0  0]
            [ 0  0|-1  0| 0  0]
            [-----+-----+-----]
            [ 0  0| 0  0| 1  0]
            [ 0  0| 0  0| 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True); certificate
            OneSumNode (6×6) with 4 children
            sage: certificate.child_nodes()
            (GraphicNode (2×2), GraphicNode (2×2), GraphicNode (1×1), GraphicNode (1×1))

            sage: M2 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 2, sparse=True),
            ....:                            [[1, 1], [-1, 0]]); M2
            [ 1  1]
            [-1  0]
            sage: result, certificate = M2.is_totally_unimodular(certificate=True); certificate
            GraphicNode (2×2)
            sage: certificate.child_nodes()
            ()
        """
        return tuple(child[0] for child in self._children())

    def child_indices(self):
        r"""
        Return a tuple of the tuples of the row and column keys of children
        in the parent node.

        OUTPUT: a tuple of (row keys, column keys)

        If the number of children is 1, then output (row keys, column keys).

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: certificate.child_indices()
            ()

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = matrix([[1, 0], [-1, 1], [0, 1]], sparse=True)
            sage: M2 = block_diagonal_matrix([M, M], sparse=True)
            sage: M2cmr = Matrix_cmr_chr_sparse(M2.parent(), M2); M2cmr
            [ 1  0  0  0]
            [-1  1  0  0]
            [ 0  1  0  0]
            [ 0  0  1  0]
            [ 0  0 -1  1]
            [ 0  0  0  1]
            sage: result, certificate = M2cmr.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, OneSumNode (6×4) with 2 children)
            sage: C = certificate.summands(); C
            (GraphicNode (3×2), GraphicNode (3×2))
            sage: certificate.child_indices()[0]
            ((0, 1, 2), (0, 1))
            sage: certificate.child_indices()[1]
            ((3, 4, 5), (2, 3))

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 12, sparse=True),
            ....: [[1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            ....: [0, 0, 0, 1, -1, 0, 0, 0, 1 , 1, 1, 1],
            ....: [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
            ....: [ 1,  0,  1,  0,  0,  0,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  1,  1,  0,  0,  0,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0],
            ....: [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1]])
            sage: result, certificate = R12.is_totally_unimodular(certificate=True,
            ....:                           row_keys=['r1', 'r2', 'r3', 'r4', 'r5',
            ....:                                     'r6', 'r7', 'r8', 'r9'],
            ....:                           column_keys=['a','b','c','d','e','f',
            ....:                                        'g','h','i','j','k','l'])
            sage: C = certificate.child_nodes()[0]; C
            ThreeSumNode (9×12) with 2 children
            sage: certificate.child_indices()
            ((r1, i, r3, r4, r5, r6, r7, r8, r9), (a, b, c, d, e, f, g, h, r2, j, k, l))

            sage: M2 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 2, sparse=True),
            ....:                            [[1, 1], [-1, 0]]); M2
            [ 1  1]
            [-1  0]
            sage: result, certificate = M2.is_totally_unimodular(certificate=True); certificate
            GraphicNode (2×2)
            sage: certificate.child_indices()
            ()
        """
        if self.nchildren() == 1:
            child = self._children()[0]
            return (child[1], child[2])
        return tuple((child[1], child[2]) for child in self._children())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 6, 6, sparse=True),
            ....: [[1,0,1,1,0,0],[0,1,1,1,0,0],[1,0,1,0,1,1],
            ....: [0,-1,0,-1,1,1],[1,0,1,0,1,0],[0,-1,0,-1,0,1]])
            sage: result, certificate = R12.is_totally_unimodular(certificate=True,
            ....:                           three_sum_strategy="Wide_Wide",
            ....:                           row_keys=range(6),
            ....:                           column_keys='abcdef')
            sage: print(certificate)
            PivotsNode (6×6)
        """
        nrows, ncols = self.dimensions()
        return f'{self.__class__.__name__} ({nrows}×{ncols})'

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of the decomposition tree rooted at ``self``.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 6, 6, sparse=True),
            ....: [[1,0,1,1,0,0],[0,1,1,1,0,0],[1,0,1,0,1,1],
            ....: [0,-1,0,-1,1,1],[1,0,1,0,1,0],[0,-1,0,-1,0,1]])
            sage: result, certificate = R12.is_totally_unimodular(certificate=True,
            ....:                           three_sum_strategy="Wide_Wide",
            ....:                           row_keys=range(6),
            ....:                           column_keys='abcdef')
            sage: unicode_art(certificate)
                    PivotsNode (6×6)
                    │
            ╭─────────────ThreeSumNode (6×6) with 2 children
            │                   │
            CographicNode (4×5) GraphicNode (4×5)
        """
        return self.as_ordered_tree()._unicode_art_()

    def _ascii_art_(self):
        r"""
        Return a ascii art representation of the decomposition tree rooted at ``self``.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 6, 6, sparse=True),
            ....: [[1,0,1,1,0,0],[0,1,1,1,0,0],[1,0,1,0,1,1],
            ....: [0,-1,0,-1,1,1],[1,0,1,0,1,0],[0,-1,0,-1,0,1]])
            sage: result, certificate = R12.is_totally_unimodular(certificate=True,
            ....:                           three_sum_strategy="Wide_Wide",
            ....:                           row_keys=range(6),
            ....:                           column_keys='abcdef')
            sage: ascii_art(certificate)
                               PivotsNode (6×6)
                               |
                      _________ThreeSumNode (6×6) with 2 children
                     /                   /
                    CographicNode (4×5) GraphicNode (4×5)
        """
        return self.as_ordered_tree()._ascii_art_()

    def one_sum(*summands, **kwds):
        r"""
        Return a :class:OneSumNode constructed from the given nodes (summands).

        INPUT:

        - ``summands`` -- decomposition nodes :class:DecompositionNode

        - ``summand_ids`` -- a tuple or list of ids for summands

        - ``row_keys`` -- a finite or enumerated family of arbitrary objects
          that index the rows of the result.
          Must be consistent with the row keys of the summands

        - ``column_keys`` -- a finite or enumerated family of arbitrary objects
          that index the columns of the result.
          Must be consistent with the column keys of the summands

        Note that ``row_keys``, ``column_keys`` of ``summands`` are disjoint

        OUTPUT: A :class:`OneSumNode`

        The terminology "1-sum" is used in the context of Seymour's decomposition
        of totally unimodular matrices and regular matroids, see [Sch1986]_.

        .. SEEALSO:: :meth:`two_sum`, :meth:`three_sum`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: from sage.matrix.seymour_decomposition import DecompositionNode
            sage: M2 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]],
            ....:                                    [[1, 1], [-1, 0]])
            sage: result, certificate = M2.is_totally_unimodular(certificate=True,
            ....:                                                row_keys=range(4),
            ....:                                                column_keys='abcd')
            sage: certificate
            OneSumNode (4×4) with 2 children
            sage: certificate.summand_matrices()
            (
            [ 1  0]  [ 1  1]
            [-1  1], [-1  0]
            )
            sage: certificate.child_indices()
            (((0, 1), ('a', 'b')), ((2, 3), ('c', 'd')))
            sage: node = DecompositionNode.one_sum(*certificate.child_nodes())
            sage: node.summand_matrices()
            (
            [ 1  0]  [ 1  1]
            [-1  1], [-1  0]
            )
            sage: node.child_indices()
            (((0, 1), ('a', 'b')), ((2, 3), ('c', 'd')))

            sage: M3 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]],
            ....:                                    [[1]], [[-1]])
            sage: result, certificate = M3.is_totally_unimodular(certificate=True,
            ....:                                                row_keys=range(4),
            ....:                                                column_keys='abcd')
            sage: certificate
            OneSumNode (4×4) with 3 children
            sage: certificate.summand_matrices()
            (
            [ 1  0]
            [-1  1], [1], [-1]
            )
            sage: certificate.child_indices()
            (((0, 1), ('a', 'b')), ((2,), ('c',)), ((3,), ('d',)))
            sage: node = DecompositionNode.one_sum(*certificate.child_nodes())
            sage: node.summand_matrices()
            (
            [ 1  0]
            [-1  1], [1], [-1]
            )
            sage: node.child_indices()
            (((0, 1), ('a', 'b')), ((2,), ('c',)), ((3,), ('d',)))

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode(matrix(ZZ, [[1, 0, 1], [0, 1, 1]]),
            ....:                    row_keys='ab',
            ....:                    column_keys=range(3)); node
            UnknownNode (2×3)
            sage: node = DecompositionNode.one_sum(certificate, node, summand_ids=range(2))
            sage: node.summand_matrices()
            (
            [ 1  0| 0| 0]
            [-1  1| 0| 0]
            [-----+--+--]
            [ 0  0| 1| 0]
            [-----+--+--]  [1 0 1]
            [ 0  0| 0|-1], [0 1 1]
            )
            sage: node.child_indices()
            ((((0, 0), (0, 1), (0, 2), (0, 3)), ((0, 'a'), (0, 'b'), (0, 'c'), (0, 'd'))),
            (((1, 'a'), (1, 'b')), ((1, 0), (1, 1), (1, 2))))

        ``row_keys``, ``column_keys`` of ``summands`` are disjoint:

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: from sage.matrix.seymour_decomposition import DecompositionNode
            sage: M2 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]],
            ....:                                    [[1, 1], [-1, 0]])
            sage: result, certificate = M2.is_totally_unimodular(certificate=True,
            ....:                                                row_keys='aefg',
            ....:                                                column_keys='abcd')
            sage: node = DecompositionNode.one_sum(*certificate.child_nodes())
            Traceback (most recent call last):
            ...
            ValueError: keys must be disjoint, got summand_row_keys=('a', 'e'), summand_column_keys=('a', 'b')

            sage: result, certificate = M2.is_totally_unimodular(certificate=True,
            ....:                                                row_keys=range(4),
            ....:                                                column_keys='abcd')
            sage: node = DecompositionNode.one_sum(*certificate.child_nodes(),
            ....:                                  row_keys=range(4),
            ....:                                  column_keys='abce')
            Traceback (most recent call last):
            ...
            ValueError: inconsistent column_keys, got column_keys=('a', 'b', 'c', 'e'), should be a permutation of ['a', 'b', 'c', 'd']
        """
        summand_ids = kwds.pop('summand_ids', None)
        row_keys = kwds.pop('row_keys', None)
        column_keys = kwds.pop('column_keys', None)
        if kwds:
            raise ValueError(f'unknown keywords: {sorted(kwds)}')

        result = OneSumNode()
        summands = tuple(summands)
        if summand_ids is not None:
            summand_ids = tuple(summand_ids)
        else:
            summand_ids = tuple(None for summand in summands)
        # TODO: Make summands DecompositionNodes if not already
        # Check row_keys, column_keys of summands are disjoint. Otherwise error
        summands_row_keys = []
        summands_column_keys = []
        row_key_list = []
        column_key_list = []
        key_set = set()
        for summand, id in zip(summands, summand_ids):
            summand_row_keys = summand.row_keys()
            summand_column_keys = summand.column_keys()
            if id is not None:
                summand_row_keys = tuple((id, key) for key in summand_row_keys)
                summand_column_keys = tuple((id, key) for key in summand_column_keys)

            old_num_keys = len(key_set)
            row_key_list.extend(summand_row_keys)
            column_key_list.extend(summand_column_keys)
            key_set.update(summand_row_keys)
            key_set.update(summand_column_keys)
            if old_num_keys + len(summand_row_keys) + len(summand_column_keys) != len(key_set):
                raise ValueError(f'keys must be disjoint, '
                                 f'got {summand_row_keys=}, {summand_column_keys=}')
            summands_row_keys.append(summand_row_keys)
            summands_column_keys.append(summand_column_keys)

        if row_keys is not None:
            row_keys = tuple(row_keys)
            if set(row_keys) != set(row_key_list) or len(row_keys) != len(row_key_list):
                raise ValueError(f'inconsistent row_keys, '
                                 f'got {row_keys=}, should be a permutation of {row_key_list}')
        else:
            row_keys = tuple(row_key_list)
        if column_keys is not None:
            column_keys = tuple(column_keys)
            if set(column_keys) != set(column_key_list) or len(column_keys) != len(column_key_list):
                raise ValueError(f'inconsistent column_keys, '
                                 f'got {column_keys=}, should be a permutation of {column_key_list}')
        else:
            column_keys = tuple(column_key_list)

        result._child_nodes = tuple(zip(summands, summands_row_keys, summands_column_keys))
        result._row_keys = row_keys
        result._column_keys = column_keys
        return result

    def _regularity(self):
        r"""
        Return whether the decomposition node is regular (binary) or TU (ternary).
        If it is not determined, raise ValueError.
        """
        cdef int8_t regularity
        if self._dec != NULL:
            regularity = CMRseymourRegularity(self._dec)
            if regularity:
                return regularity > 0
        raise ValueError('It is not determined whether the decomposition node is regular/TU')

    def _graphicness(self):
        r"""
        Return whether the decomposition node is graphic (binary) or network (ternary).
        If it is not determined, raise ValueError.
        """
        cdef int8_t graphicness
        if self._dec != NULL:
            graphicness = CMRseymourGraphicness(self._dec)
            if graphicness:
                return graphicness > 0
        raise ValueError('It is not determined whether the decomposition node is graphic/network')

    def _cographicness(self):
        r"""
        Return whether the decomposition node is cographic (binary) or conetwork (ternary).
        If it is not determined, raise ValueError.
        """
        cdef int8_t cographicness
        if self._dec != NULL:
            cographicness = CMRseymourCographicness(self._dec)
            if cographicness:
                return cographicness > 0
        raise ValueError('It is not determined whether the decomposition node is cographic/conetwork')

    def _is_binary_linear_matroid_graphic(self, *, decomposition=False, **kwds):
        r"""
        Return whether the linear matroid of ``self`` over `\GF{2}` is graphic.
        If there is some entry not in `\{0, 1\}`, return ``False``.

        This method is based on Seymour's decomposition.
        The decomposition will stop once nongraphicness is detected.
        For direct graphicness check,

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.
            Matrix_cmr_chr_sparse._is_binary_linear_matroid_graphic`
            :meth:`UnknownNode._is_binary_linear_matroid_graphic`
            :meth:`_binary_linear_matroid_complete_decomposition`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: from sage.matrix.seymour_decomposition import DecompositionNode
            sage: M2 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [1, 1], [0, 1]],
            ....:                                    [[1, 1, 0], [0, 1, 1]])
            sage: result, certificate = M2.is_totally_unimodular(certificate=True,
            ....:                                                row_keys=range(5),
            ....:                                                column_keys='abcde')
            sage: node = DecompositionNode.one_sum(*certificate.child_nodes())
            sage: node._graphicness()
            Traceback (most recent call last):
            ...
            ValueError: It is not determined whether the decomposition node is graphic/network
            sage: result, decomposition = node._is_binary_linear_matroid_graphic(decomposition=True)
            sage: result
            True
            sage: unicode_art(decomposition)
            ╭──────────OneSumNode (5×5) with 2 children
            │                │
            PlanarNode (3×2) PlanarNode (2×3)
            sage: decomposition._graphicness()
            True
            sage: node._is_binary_linear_matroid_graphic(certificate=True)
            Traceback (most recent call last):
            ...
            NotImplementedError

            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(GF(2), 9, 9, sparse=True),
            ....:                           [[1, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 1, 1, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 1, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 1, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 0, 1, 0],
            ....:                            [0, 0, 0, 0, 0, 1, 0, 1, 1],
            ....:                            [0, 0, 0, 0, 0, 0, 1, 1, 1]])
            sage: result, certificate = M._is_binary_linear_matroid_regular(certificate=True,
            ....:                         row_keys=[f'r{i}' for i in range(9)],
            ....:                         column_keys=[f'c{i}' for i in range(9)])
            sage: node = DecompositionNode.one_sum(*certificate.child_nodes())
            sage: result, decomposition = node._is_binary_linear_matroid_graphic(decomposition=True)
            sage: result
            False
            sage: unicode_art(decomposition)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) UnknownNode (4×5)
            sage: decomposition.child_nodes()[1]._graphicness()
            False
        """
        certificate = kwds.get('certificate', False)
        try:
            result = self._graphicness()
            if not decomposition and not certificate:
                return result
            result = [result]
            if decomposition:
                result.append(self)
            if certificate:
                raise NotImplementedError
            return result
        except ValueError:
            # compute it... wait for CMR functions
            full_dec = self._binary_linear_matroid_complete_decomposition(
                                        stop_when_nongraphic=True,
                                        check_graphic_minors_planar=True)
            return full_dec._is_binary_linear_matroid_graphic(
                                       decomposition=decomposition,
                                       certificate=certificate)

    def _is_binary_linear_matroid_cographic(self, *, decomposition=False, **kwds):
        r"""
        Return whether the linear matroid of ``self`` over `\GF{2}` is cographic.
        If there is some entry not in `\{0, 1\}`, return ``False``.

        This method is based on Seymour's decomposition.
        The decomposition will stop once noncographicness is detected.
        For direct cographicness check,

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.
            Matrix_cmr_chr_sparse._is_binary_linear_matroid_cographic`
            :meth:`UnknownNode._is_binary_linear_matroid_cographic`
            :meth:`_binary_linear_matroid_complete_decomposition`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: from sage.matrix.seymour_decomposition import DecompositionNode
            sage: M2 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [1, 1], [0, 1]],
            ....:                                    [[1, 1, 0], [0, 1, 1]])
            sage: result, certificate = M2.is_totally_unimodular(certificate=True,
            ....:                                                row_keys=range(5),
            ....:                                                column_keys='abcde')
            sage: node = DecompositionNode.one_sum(*certificate.child_nodes())
            sage: node._graphicness()
            Traceback (most recent call last):
            ...
            ValueError: It is not determined whether the decomposition node is graphic/network
            sage: result, decomposition = node._is_binary_linear_matroid_cographic(decomposition=True)
            sage: result
            True
            sage: unicode_art(decomposition)
            ╭──────────OneSumNode (5×5) with 2 children
            │                │
            PlanarNode (3×2) PlanarNode (2×3)
            sage: decomposition._cographicness()
            True
            sage: node._is_binary_linear_matroid_cographic(certificate=True)
            Traceback (most recent call last):
            ...
            NotImplementedError

            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(GF(2), 9, 9, sparse=True),
            ....:                           [[1, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 1, 1, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 1, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 1, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 0, 1, 0],
            ....:                            [0, 0, 0, 0, 0, 1, 0, 1, 1],
            ....:                            [0, 0, 0, 0, 0, 0, 1, 1, 1]])
            sage: result, certificate = M._is_binary_linear_matroid_regular(certificate=True,
            ....:                         row_keys=[f'r{i}' for i in range(9)],
            ....:                         column_keys=[f'c{i}' for i in range(9)])
            sage: node = DecompositionNode.one_sum(*certificate.child_nodes())
            sage: result, decomposition = node._is_binary_linear_matroid_cographic(decomposition=True)
            sage: result
            False
            sage: unicode_art(decomposition)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) UnknownNode (4×5)
        """
        certificate = kwds.get('certificate', False)
        try:
            result = self._cographicness()
            if not decomposition and not certificate:
                return result
            result = [result]
            if decomposition:
                result.append(self)
            if certificate:
                raise NotImplementedError
            return result
        except ValueError:
            # compute it... wait for CMR functions
            full_dec = self._binary_linear_matroid_complete_decomposition(
                                        stop_when_noncographic=True,
                                        check_graphic_minors_planar=True)
            return full_dec._is_binary_linear_matroid_cographic(
                                       decomposition=decomposition,
                                       certificate=certificate)

    def _is_binary_linear_matroid_regular(self, *, decomposition=False, **kwds):
        r"""
        Return whether the linear matroid of ``self`` over `\GF{2}` is regular.
        If there is some entry not in `\{0, 1\}`, return ``False``.

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.
            Matrix_cmr_chr_sparse._is_binary_linear_matroid_regular`
            :meth:`_binary_linear_matroid_complete_decomposition`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(GF(2), 9, 9, sparse=True),
            ....:                           [[1, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 1, 1, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 1, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 1, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 0, 1, 0],
            ....:                            [0, 0, 0, 0, 0, 1, 0, 1, 1],
            ....:                            [0, 0, 0, 0, 0, 0, 1, 1, 1]])
            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode(M)
            sage: node._regularity()
            Traceback (most recent call last):
            ...
            ValueError: It is not determined whether the decomposition node is regular/TU
            sage: node._is_binary_linear_matroid_regular()
            True
            sage: C0 = node._binary_linear_matroid_complete_decomposition()
            sage: C0
            OneSumNode (9×9) with 2 children
            sage: result, decomposition = C0._is_binary_linear_matroid_graphic(decomposition=True)
            sage: result
            False
            sage: unicode_art(decomposition)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) CographicNode (4×5)
            sage: decomposition.child_nodes()[1]._graphicness()
            False
            sage: result, decomposition = C0._is_binary_linear_matroid_cographic(decomposition=True)
            sage: result
            False
            sage: unicode_art(decomposition)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) UnknownNode (4×5)
            sage: decomposition.child_nodes()[1]._graphicness()
            Traceback (most recent call last):
            ...
            ValueError: It is not determined whether the decomposition node is graphic/network
            sage: result, decomposition = C0._is_binary_linear_matroid_regular(decomposition=True)
            sage: result
            True
            sage: unicode_art(decomposition)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) CographicNode (4×5)
        """
        certificate = kwds.get('certificate', False)
        try:
            result = self._regularity()
            if not decomposition and not certificate:
                return result
            result = [result]
            if decomposition:
                result.append(self)
            if certificate:
                raise NotImplementedError
            return result
        except ValueError:
            # compute it... wait for CMR functions
            full_dec = self._binary_linear_matroid_complete_decomposition(
                                        stop_when_irregular=True,
                                        check_graphic_minors_planar=True)
            return full_dec._is_binary_linear_matroid_regular(
                                       decomposition=decomposition,
                                       certificate=certificate)

    def _binary_linear_matroid_complete_decomposition(self, *,
                                        time_limit=60.0,
                                        use_direct_graphicness_test=True,
                                        prefer_graphicness=True,
                                        series_parallel_ok=True,
                                        check_graphic_minors_planar=False,
                                        stop_when_irregular=False,
                                        stop_when_nongraphic=False,
                                        stop_when_noncographic=False,
                                        stop_when_nongraphic_and_noncographic=False,
                                        three_sum_pivot_children=False,
                                        three_sum_strategy=None,
                                        construct_leaf_graphs=False,
                                        construct_all_graphs=False):
        r"""
        Complete the Seymour's decomposition of ``self`` over `\GF{2}`.

        INPUT:

        - ``stop_when_irregular`` -- boolean
        Whether to stop decomposing once irregularity is determined.

        - ``stop_when_nongraphic`` -- boolean
        Whether to stop decomposing once non-graphicness is determined.

        - ``stop_when_noncographic`` -- boolean
        Whether to stop decomposing once non-cographicness is determined.

        - ``stop_when_nongraphic_and_noncographic`` -- boolean
        Whether to stop decomposing once non-graphicness and non-cographicness
        is determined.

          For a description of other parameters, see
          :meth:`sage.matrix.matrix_cmr_sparse._set_cmr_seymour_parameters`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(GF(2), 9, 9, sparse=True),
            ....:                           [[1, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 1, 1, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 1, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 1, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 0, 1, 0],
            ....:                            [0, 0, 0, 0, 0, 1, 0, 1, 1],
            ....:                            [0, 0, 0, 0, 0, 0, 1, 1, 1]])
            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode(M); node
            UnknownNode (9×9)
            sage: C0 = node._binary_linear_matroid_complete_decomposition()
            sage: C0
            OneSumNode (9×9) with 2 children
            sage: unicode_art(C0)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) CographicNode (4×5)
            sage: C1, C2 = C0.child_nodes()
            sage: C11 = C1._binary_linear_matroid_complete_decomposition(); C11
            GraphicNode (5×4)
            sage: unicode_art(C11)
            GraphicNode (5×4)
            sage: C1.matrix()
            [1 1 0 0]
            [1 1 1 0]
            [0 1 1 1]
            [1 0 0 1]
            [0 0 1 1]
            sage: C11.matrix()
            [1 1 0 0]
            [1 1 1 0]
            [0 1 1 1]
            [1 0 0 1]
            [0 0 1 1]
            sage: C22 = C2._binary_linear_matroid_complete_decomposition(); C22
            CographicNode (4×5)
            sage: unicode_art(C22)
            CographicNode (4×5)
            sage: C2.matrix()
            [1 1 1 0 0]
            [0 0 1 1 1]
            [0 1 0 1 1]
            [1 1 0 1 0]
            sage: C22.matrix()
            [1 1 1 0 0]
            [0 0 1 1 1]
            [0 1 0 1 1]
            [1 1 0 1 0]

        This is test ``TreeFlagsStopNoncographic`` in CMR's ``test_regular.cpp``.
        The default settings will not do the planarity check::

            sage: unicode_art(node)
            UnknownNode (9×9)
            sage: certificate1 = node._binary_linear_matroid_complete_decomposition(
            ....:                                           stop_when_noncographic=True,
            ....:                                           check_graphic_minors_planar=True)
            sage: unicode_art(certificate1)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) UnknownNode (4×5)
            sage: certificate2 = node._binary_linear_matroid_complete_decomposition(
            ....:                                           stop_when_noncographic=True)
            sage: unicode_art(certificate2)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) CographicNode (4×5)

            sage: certificate1 = node._binary_linear_matroid_complete_decomposition(
            ....:                                           stop_when_nongraphic=True,
            ....:                                           check_graphic_minors_planar=True)
            sage: unicode_art(certificate1)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) UnknownNode (4×5)
            sage: C1, C2 = certificate1.child_nodes()
            sage: C1._cographicness()
            False
            sage: C2._graphicness()
            False
            sage: certificate2 = node._binary_linear_matroid_complete_decomposition(
            ....:                                           stop_when_nongraphic=True)
            sage: unicode_art(certificate2)
            ╭───────────OneSumNode (9×9) with 2 children
            │                 │
            GraphicNode (5×4) UnknownNode (4×5)
            sage: C1, C2 = certificate2.child_nodes()
            sage: C1._cographicness()
            Traceback (most recent call last):
            ...
            ValueError: It is not determined whether the decomposition node is cographic/conetwork
            sage: C2._graphicness()
            False
        """
        cdef CMR_REGULAR_PARAMS params
        cdef CMR_REGULAR_STATS stats
        cdef CMR_SEYMOUR_NODE *clone = NULL

        cdef CMR_SEYMOUR_NODE **pclone = &clone

        if self._dec == NULL:
            base_ring = self.base_ring()
            if base_ring is None:
                from sage.rings.finite_rings.finite_field_constructor import GF
                self._base_ring = GF(2)
            elif base_ring.characteristic() != 2:
                raise ValueError(f'only defined over binary, got {base_ring}')
            self._set_root_dec()

        cdef dict kwds = dict(use_direct_graphicness_test=use_direct_graphicness_test,
                              prefer_graphicness=prefer_graphicness,
                              series_parallel_ok=series_parallel_ok,
                              check_graphic_minors_planar=check_graphic_minors_planar,
                              stop_when_irregular=stop_when_irregular,
                              stop_when_nongraphic=stop_when_nongraphic,
                              stop_when_noncographic=stop_when_noncographic,
                              stop_when_nongraphic_and_noncographic=stop_when_nongraphic_and_noncographic,
                              three_sum_pivot_children=three_sum_pivot_children,
                              three_sum_strategy=three_sum_strategy,
                              construct_leaf_graphs=construct_leaf_graphs,
                              construct_all_graphs=construct_all_graphs)
        _set_cmr_seymour_parameters(&params.seymour, kwds)

        sig_on()
        try:
            CMR_CALL(CMRseymourCloneUnknown(cmr, self._dec, pclone))
            CMR_CALL(CMRregularCompleteDecomposition(cmr, clone, &params, &stats, time_limit))
        finally:
            sig_off()
        node = create_DecompositionNode(clone, matrix=self.matrix(),
                                        row_keys=self.row_keys(),
                                        column_keys=self.column_keys(),
                                        base_ring=self.base_ring())
        return node

    def is_network_matrix(self, *, decomposition=False, **kwds):
        r"""

        """
        certificate = kwds.get('certificate', False)
        try:
            result = self._graphicness()
            if not decomposition and not certificate:
                return result
            result = [result]
            if decomposition:
                result.append(self)
            if certificate:
                raise NotImplementedError
            return result
        except ValueError:
            # compute it... wait for CMR functions
            raise NotImplementedError("Network Not Determined")

    def is_conetwork_matrix(self, *, decomposition=False, **kwds):
        r"""

        """
        certificate = kwds.get('certificate', False)
        try:
            result = self._cographicness()
            if not decomposition and not certificate:
                return result
            result = [result]
            if decomposition:
                result.append(self)
            if certificate:
                raise NotImplementedError
            return result
        except ValueError:
            # compute it... wait for CMR functions
            raise NotImplementedError("Conetwork Not Determined")

    def is_totally_unimodular(self, *, decomposition=False, **kwds):
        r"""

        """
        certificate = kwds.get('certificate', False)
        try:
            result = self._regularity()
            if not decomposition and not certificate:
                return result
            result = [result]
            if decomposition:
                result.append(self)
            if certificate:
                raise NotImplementedError
            return result
        except ValueError:
            # compute it... wait for CMR functions
            raise NotImplementedError("TU Not Determined")

    def nminors(self):
        r"""
        Return the number of minors of the node.
        """
        if self._minors is not None:
            return len(self._minors)
        if self._dec == NULL:
            return 0
        return CMRseymourNumMinors(self._dec)

    def _create_minor(self, index):
        r"""
        A minor of a represented matroid is another one obtained
        by deleting or contracting elements. In the reduced matrix representation,
        an element associated with a row (resp. column) can be contracted (resp. deleted)
        by removing the corresponding matrix row (resp. column).
        In order to delete a row element or contract a column element,
        one must first pivot such that the element type (row/column) is changed.

        A minor of a matroid represented by matrix `M` is represented
        by means of a ``CMR_MINOR`` object.
        It consists of a (potentially empty) array of pivots and a ``CMR_SUBMAT`` object
        indicating a submatrix `M'` of the matrix obtained from `M` after applying the pivots.
        Moreover, the type field indicates a certain structure of `M'`:

        - A determinant ``CMR_MINOR_TYPE_DETERMINANT``
          indicates that `M'` is a submatrix of `M` with `|\det(M')| \geq 2`.
          In particular, no pivots are applied.
        - A Fano ``CMR_MINOR_TYPE_FANO``
          indicates that `M'` represents the Fano matroid `F_7`.
        - A Fano-dual ``CMR_MINOR_TYPE_FANO_DUAL``
          indicates that `M'` represents the dual `F_7^\star` of the Fano matroid.
        - A K5 ``CMR_MINOR_TYPE_K5``
          indicates that `M'` represents the graphic matroid `M(K_5)` of the complete graph `K_5`.
        - A K5-dual ``CMR_MINOR_TYPE_K5_DUAL``
          indicates that `M'` represents the dual matroid `M(K_5)^\star` of the complete graph `K_5`.
        - A K33 ``CMR_MINOR_TYPE_K33``
          indicates that `M'` represents the graphic matroid `M(K_{3,3})`
          of the complete bipartite graph `K_{3,3}`.
        - A K33-dual ``CMR_MINOR_TYPE_K33_DUAL``
          indicates that `M'` represents the dual matroid `M(K_{3,3})^\star`
          of the complete bipartite graph `K_{3,3}`.
        """
        cdef CMR_MINOR * minor = CMRseymourMinor(self._dec, index)
        cdef CMR_MINOR_TYPE typ = CMRminorType(minor)
        cdef size_t npivots = CMRminorNumPivots(minor)
        cdef size_t *pivot_rows = CMRminorPivotRows(minor)
        cdef size_t *pivot_columns = CMRminorPivotColumns(minor)
        cdef CMR_SUBMAT *submat = CMRminorSubmatrix(minor)
        # TODO: consider the pivots information and the corresponding keys?
        pivots_tuple = tuple((pivot_rows[i], pivot_columns[i]) for i in range(npivots))
        submat_tuple = (tuple(submat.rows[i] for i in range(submat.numRows)),
                            tuple(submat.columns[i] for i in range(submat.numColumns)))
        import sage.matroids.matroids_catalog as matroids
        from sage.graphs.graph_generators import graphs
        from sage.matroids.matroid import Matroid

        if typ == CMR_MINOR_TYPE_FANO:
            return matroids.catalog.Fano()
        if typ == CMR_MINOR_TYPE_FANO_DUAL:
            return matroids.catalog.Fano().dual()
        if typ == CMR_MINOR_TYPE_K5:
            return matroids.CompleteGraphic(5)
        if typ == CMR_MINOR_TYPE_K5_DUAL:
            return matroids.CompleteGraphic(5).dual()
        if typ == CMR_MINOR_TYPE_K33:
            E = 'abcdefghi'
            G = graphs.CompleteBipartiteGraph(3, 3)
            return Matroid(groundset=E, graph=G, regular=True)
        if typ == CMR_MINOR_TYPE_K33_DUAL:
            return matroids.catalog.K33dual()
        if typ == CMR_MINOR_TYPE_DETERMINANT:
            return '|det| = 2 submatrix', submat_tuple
        if typ == CMR_MINOR_TYPE_ENTRY:
            return 'bad entry'
        if typ == CMR_MINOR_TYPE_CUSTOM:
            return 'custom'

    def minors(self):
        r"""
        Return a tuple of minors of ``self``.

        Currently, it only handles determinant 2 submatrix.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 9, sparse=True),
            ....:                           [[1, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 1, 1, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 1, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 1, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 0, 1, 0],
            ....:                            [0, 0, 0, 0, 0, 1, 0, 1, 1],
            ....:                            [0, 0, 0, 0, 0, 0, 1, 1, 1]])
            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: C0 = UnknownNode(M).complete_decomposition()
            sage: C1, C2 = C0.child_nodes()
            sage: C1.minors()
            (('|det| = 2 submatrix', ((4, 3, 1), (2, 3, 0))),)
            sage: C2.minors()
            (('|det| = 2 submatrix', ((2, 1, 0), (4, 2, 1))),)
        """
        if self._minors is not None:
            return self._minors
        minors_tuple = tuple(self._create_minor(index)
                               for index in range(self.nminors()))
        self._minors = minors_tuple
        return self._minors

    def complete_decomposition(self, *, time_limit=60.0,
                               use_direct_graphicness_test=True,
                               prefer_graphicness=True,
                               series_parallel_ok=True,
                               check_graphic_minors_planar=False,
                               stop_when_nonTU=False,
                               stop_when_nonnetwork=False,
                               stop_when_nonconetwork=False,
                               stop_when_nonnetwork_and_nonconetwork=False,
                               three_sum_pivot_children=False,
                               three_sum_strategy=None,
                               construct_leaf_graphs=False,
                               construct_all_graphs=False):
        r"""
        Complete the Seymour's decomposition of ``self`` over `\GF{3}` or `\ZZ`.

        INPUT:

        - ``stop_when_nonTU`` -- boolean
        Whether to stop decomposing once being non-TU is determined.

        - ``stop_when_nonnetwork`` -- boolean
        Whether to stop decomposing once being non-network is determined.

        - ``stop_when_nonconetwork`` -- boolean
        Whether to stop decomposing once being non-conetwork is determined.

        - ``stop_when_nonnetwork_and_nonconetwork`` -- boolean
        Whether to stop decomposing once not being network
        and not being conetwork is determined.

          For a description of other parameters, see
          :meth:`sage.matrix.matrix_cmr_sparse._set_cmr_seymour_parameters`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 9, sparse=True),
            ....:                           [[1, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 1, 1, 0, 0, 0, 0, 0, 0],
            ....:                            [1, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 1, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 1, 1, 0, 0, 0, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 1, 0, 0],
            ....:                            [0, 0, 0, 0, 1, 1, 0, 1, 0],
            ....:                            [0, 0, 0, 0, 0, 1, 0, 1, 1],
            ....:                            [0, 0, 0, 0, 0, 0, 1, 1, 1]])
            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode(M); node
            UnknownNode (9×9)
            sage: C0 = node.complete_decomposition()
            sage: unicode_art(C0)
            ╭OneSumNode (9×9) with 2 children─╮
            │                                 │
            ThreeConnectedIrregularNode (5×4) ThreeConnectedIrregularNode (4×5)
            sage: C1, C2 = C0.child_nodes()
            sage: C11 = C1.complete_decomposition(); C11
            ThreeConnectedIrregularNode (5×4)
            sage: unicode_art(C11)
            ThreeConnectedIrregularNode (5×4)
            sage: C1.matrix()
            [1 1 0 0]
            [1 1 1 0]
            [0 1 1 1]
            [1 0 0 1]
            [0 0 1 1]
            sage: C11.matrix()
            [1 1 0 0]
            [1 1 1 0]
            [0 1 1 1]
            [1 0 0 1]
            [0 0 1 1]
            sage: C22 = C2.complete_decomposition(); C22
            ThreeConnectedIrregularNode (4×5)
            sage: unicode_art(C22)
            ThreeConnectedIrregularNode (4×5)
            sage: C2.matrix()
            [1 1 1 0 0]
            [0 0 1 1 1]
            [0 1 0 1 1]
            [1 1 0 1 0]
            sage: C22.matrix()
            [1 1 1 0 0]
            [0 0 1 1 1]
            [0 1 0 1 1]
            [1 1 0 1 0]
        """
        cdef CMR_TU_PARAMS params
        cdef CMR_TU_STATS stats
        cdef CMR_SEYMOUR_NODE *clone = NULL

        cdef CMR_SEYMOUR_NODE **pclone = &clone

        if self._dec == NULL:
            self._set_root_dec()

        cdef dict kwds = dict(use_direct_graphicness_test=use_direct_graphicness_test,
                              prefer_graphicness=prefer_graphicness,
                              series_parallel_ok=series_parallel_ok,
                              check_graphic_minors_planar=check_graphic_minors_planar,
                              stop_when_irregular=stop_when_nonTU,
                              stop_when_nongraphic=stop_when_nonnetwork,
                              stop_when_noncographic=stop_when_nonconetwork,
                              stop_when_nongraphic_and_noncographic=stop_when_nonnetwork_and_nonconetwork,
                              three_sum_pivot_children=three_sum_pivot_children,
                              three_sum_strategy=three_sum_strategy,
                              construct_leaf_graphs=construct_leaf_graphs,
                              construct_all_graphs=construct_all_graphs)
        params.algorithm = CMR_TU_ALGORITHM_DECOMPOSITION
        params.ternary = True
        params.camionFirst = False
        _set_cmr_seymour_parameters(&params.seymour, kwds)

        sig_on()
        try:
            CMR_CALL(CMRseymourCloneUnknown(cmr, self._dec, pclone))
            CMR_CALL(CMRtuCompleteDecomposition(cmr, clone, &params, &stats, time_limit))
        finally:
            sig_off()
        node = create_DecompositionNode(clone, matrix=self.matrix(),
                                        row_keys=self.row_keys(),
                                        column_keys=self.column_keys(),
                                        base_ring=self.base_ring())
        return node


cdef class ThreeConnectedIrregularNode(DecompositionNode):

    pass


cdef class UnknownNode(DecompositionNode):
    r"""
    EXAMPLES::

        sage: from sage.matrix.seymour_decomposition import UnknownNode
        sage: node = UnknownNode([[1, 1], [0, 1]]); node
        UnknownNode (2×2)
        sage: node.matrix()
        [1 1]
        [0 1]
        sage: node = UnknownNode(matrix(ZZ, [[1, 0, 1], [0, 1, 1]])); node
        UnknownNode (2×3)
        sage: node.matrix()
        [1 0 1]
        [0 1 1]
        sage: node = UnknownNode(matrix(ZZ, [[1, 0, 1], [0, 1, 1]]),
        ....:                    row_keys='ab',
        ....:                    column_keys=range(3)); node
        UnknownNode (2×3)
        sage: node.matrix()
        [1 0 1]
        [0 1 1]
        sage: node.morphism()._unicode_art_matrix()
          0 1 2
        a⎛1 0 1⎞
        b⎝0 1 1⎠

    From a module morphism::

        sage: phi = matrix(ZZ, [[1, 0, 1], [0, 1, 1]],
        ....:              row_keys='ab', column_keys=range(3)); phi
        Generic morphism:
        From: Free module generated by {0, 1, 2} over Integer Ring
        To:   Free module generated by {'a', 'b'} over Integer Ring
        sage: node = UnknownNode(phi); node
        UnknownNode (2×3)
        sage: node.matrix()
        [1 0 1]
        [0 1 1]
    """

    def _is_binary_linear_matroid_graphic(self, *, decomposition=False, certificate=False, **kwds):
        r"""
        Return whether the linear matroid of ``self`` over `\GF{2}` is graphic.
        If there is some entry not in `\{0, 1\}`, return ``False``.

        This is an internal method because it should really be exposed
        as a method of :class:`Matroid`.

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.
            Matrix_cmr_chr_sparse._is_binary_linear_matroid_graphic`

        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, 0], [1, 1], [0, 1]]); node
            UnknownNode (3×2)
            sage: node.matrix()
            [1 0]
            [1 1]
            [0 1]
            sage: node._is_binary_linear_matroid_graphic()
            True
            sage: result, certificate = node._is_binary_linear_matroid_graphic(certificate=True)
            sage: graph, forest_edges, coforest_edges = certificate
            sage: graph.vertices(sort=True)  # the numbers have no meaning
            [1, 2, 7, 12]
            sage: graph.edges(sort=True, labels=False)
            [(1, 2), (1, 7), (1, 12), (2, 7), (7, 12)]
            sage: forest_edges    # indexed by rows of M
            ((1, 2), (7, 1), (12, 7))
            sage: coforest_edges  # indexed by cols of M
            ((2, 7), (1, 12))
        """
        matrix = self.matrix()
        if not decomposition and not certificate:
            return matrix._is_binary_linear_matroid_graphic(**kwds)
        result, cert = matrix._is_binary_linear_matroid_graphic(certificate=True,
                                         row_keys=self.row_keys(),
                                         column_keys=self.column_keys(), **kwds)
        result = [result]
        if decomposition:
            graph, forest_edges, coforest_edges = cert
            node = GraphicNode(matrix, graph, forest_edges, coforest_edges)
            result.append(node)
        if certificate:
            result.append(cert)
        return result

    def _is_binary_linear_matroid_cographic(self, *, decomposition=False, certificate=False, **kwds):
        r"""
        Return whether the linear matroid of ``self`` over `\GF{2}` is cographic.
        If there is some entry not in `\{0, 1\}`, return ``False``.

        This is an internal method because it should really be exposed
        as a method of :class:`Matroid`.

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.
            Matrix_cmr_chr_sparse._is_binary_linear_matroid_cographic`

        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, 1, 0], [0, 1, 1]]); node
            UnknownNode (2×3)
            sage: node.matrix()
            [1 1 0]
            [0 1 1]
            sage: node._is_binary_linear_matroid_cographic()
            True
            sage: result, certificate = node._is_binary_linear_matroid_cographic(certificate=True)
            sage: graph, forest_edges, coforest_edges = certificate
            sage: graph.vertices(sort=True)  # the numbers have no meaning
            [1, 2, 7, 12]
            sage: graph.edges(sort=True, labels=False)
            [(1, 2), (1, 7), (1, 12), (2, 7), (7, 12)]
            sage: forest_edges    # indexed by rows of M
            ((1, 2), (7, 1))
            sage: coforest_edges  # indexed by cols of M
            ((2, 7), (1, 12), (1, 2))
        """
        matrix = self.matrix()
        if not decomposition and not certificate:
            return matrix._is_binary_linear_matroid_cographic(**kwds)
        result, cert = matrix._is_binary_linear_matroid_cographic(certificate=True,
                                         row_keys=self.row_keys(),
                                         column_keys=self.column_keys(), **kwds)
        result = [result]
        if decomposition:
            graph, forest_edges, coforest_edges = cert
            node = CographicNode(matrix, graph, forest_edges, coforest_edges)
            result.append(node)
        if certificate:
            result.append(cert)
        return result

    def is_network_matrix(self, *, decomposition=False, certificate=False, **kwds):
        r"""
        Return whether the matrix ``self`` over `\GF{3}` or `QQ` is a network matrix.
        If there is some entry not in `\{-1, 0, 1\}`, return ``False``.

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.
            Matrix_cmr_chr_sparse.is_network_matrix`

        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, 0], [-1, 1], [0, -1]]); node
            UnknownNode (3×2)
            sage: node.matrix()
            [ 1  0]
            [-1  1]
            [ 0 -1]
            sage: node.is_network_matrix()
            True
            sage: result, certificate = node.is_network_matrix(certificate=True)
            sage: graph, forest_edges, coforest_edges = certificate
            sage: graph
            Digraph on 4 vertices
            sage: graph.vertices(sort=True)  # the numbers have no meaning
            [1, 2, 7, 12]
            sage: graph.edges(sort=True, labels=False)
            [(2, 1), (2, 7), (7, 1), (7, 12), (12, 1)]
            sage: forest_edges    # indexed by rows of M
            ((2, 1), (7, 1), (7, 12))
            sage: coforest_edges  # indexed by cols of M
            ((2, 7), (12, 1))
        """
        matrix = self.matrix()
        if not decomposition and not certificate:
            return matrix.is_network_matrix(**kwds)
        result, cert = matrix.is_network_matrix(certificate=True,
                                         row_keys=self.row_keys(),
                                         column_keys=self.column_keys(), **kwds)
        result = [result]
        if decomposition:
            graph, forest_edges, coforest_edges = cert
            node = GraphicNode(matrix, graph, forest_edges, coforest_edges)
            result.append(node)
        if certificate:
            result.append(cert)
        return result

    def is_conetwork_matrix(self, *, decomposition=False, certificate=False, **kwds):
        r"""
        Return whether the matrix ``self`` over `\GF{3}` or `QQ` is a conetwork matrix.
        If there is some entry not in `\{-1, 0, 1\}`, return ``False``.

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.
            Matrix_cmr_chr_sparse.is_conetwork_matrix`

        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode([[1, -1, 0], [0, 1, -1]]); node
            UnknownNode (2×3)
            sage: node.matrix()
            [ 1 -1  0]
            [ 0  1 -1]
            sage: node.is_conetwork_matrix()
            True
            sage: result, certificate = node.is_conetwork_matrix(certificate=True)
            sage: graph, forest_edges, coforest_edges = certificate
            sage: graph
            Digraph on 4 vertices
            sage: graph.vertices(sort=True)  # the numbers have no meaning
            [1, 2, 7, 12]
            sage: graph.edges(sort=True, labels=False)
            [(2, 1), (2, 7), (7, 1), (7, 12), (12, 1)]
            sage: forest_edges    # indexed by rows of M
            ((2, 1), (7, 1))
            sage: coforest_edges  # indexed by cols of M
            ((2, 7), (12, 1), (2, 1))
        """
        matrix = self.matrix()
        if not decomposition and not certificate:
            return matrix.is_conetwork_matrix(**kwds)
        result, cert = matrix.is_conetwork_matrix(certificate=True,
                                         row_keys=self.row_keys(),
                                         column_keys=self.column_keys(), **kwds)
        result = [result]
        if decomposition:
            graph, forest_edges, coforest_edges = cert
            node = CographicNode(matrix, graph, forest_edges, coforest_edges)
            result.append(node)
        if certificate:
            result.append(cert)
        return result


cdef class SumNode(DecompositionNode):
    r"""
    Base class for 1-sum, 2-sum, and 3-sum nodes in Seymour's decomposition
    """

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M2 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]],
            ....:                                    [[1, 1], [-1, 0]])
            sage: result, certificate = M2.is_totally_unimodular(certificate=True,
            ....:                                                row_keys=range(4),
            ....:                                                column_keys='abcd')
            sage: print(certificate)
            OneSumNode (4×4) with 2 children
        """
        result = super()._repr_()
        result += f' with {self.nchildren()} children'
        return result

    def permuted_block_matrix(self):
        r"Return (Prow, BlockMatrix, Pcolumn) so that self.matrix() == Prow * BlockMatrix * Pcolumn ????"
        raise NotImplementedError

    summands = DecompositionNode.child_nodes

    def summand_matrices(self):
        r"""
        Return a tuple of matrices representing the child nodes.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M2 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]],
            ....:                                    [[1, 1], [-1, 0]])
            sage: result, certificate = M2.is_totally_unimodular(certificate=True,
            ....:                                                row_keys=range(4),
            ....:                                                column_keys='abcd')
            sage: certificate.summand_matrices()
            (
            [ 1  0]  [ 1  1]
            [-1  1], [-1  0]
            )
        """
        return tuple(s.matrix() for s in self.summands())


cdef class OneSumNode(SumNode):

    def block_matrix_form(self):
        r"""
        Return the block matrix representing the one sum node.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]], [[1, 1], [-1, 0]])
            sage: result, certificate = M.is_totally_unimodular(certificate=True); certificate
            OneSumNode (4×4) with 2 children
            sage: certificate.summand_matrices()
            (
            [ 1  0]  [ 1  1]
            [-1  1], [-1  0]
            )
            sage: certificate.block_matrix_form()
            [ 1  0| 0  0]
            [-1  1| 0  0]
            [-----+-----]
            [ 0  0| 1  1]
            [ 0  0|-1  0]

            sage: M3 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]],
            ....:                                    [[1, 1], [-1, 0]],
            ....:                                    [[1, 0], [0, 1]]); M3
            [ 1  0| 0  0| 0  0]
            [-1  1| 0  0| 0  0]
            [-----+-----+-----]
            [ 0  0| 1  1| 0  0]
            [ 0  0|-1  0| 0  0]
            [-----+-----+-----]
            [ 0  0| 0  0| 1  0]
            [ 0  0| 0  0| 0  1]
            sage: result, certificate = M3.is_totally_unimodular(certificate=True); certificate
            OneSumNode (6×6) with 4 children
            sage: certificate.summand_matrices()
            (
            [ 1  0]  [ 1  1]
            [-1  1], [-1  0], [1], [1]
            )
            sage: certificate.block_matrix_form()
            [ 1  0| 0  0| 0| 0]
            [-1  1| 0  0| 0| 0]
            [-----+-----+--+--]
            [ 0  0| 1  1| 0| 0]
            [ 0  0|-1  0| 0| 0]
            [-----+-----+--+--]
            [ 0  0| 0  0| 1| 0]
            [-----+-----+--+--]
            [ 0  0| 0  0| 0| 1]
        """
        return Matrix_cmr_chr_sparse.one_sum(*self.summand_matrices())

    @staticmethod
    def check(result_matrix, summand_matrices, summand_parent_rows_and_columns):
        r"""
        Check that ``result_matrix`` is a 1-sum of ``summand_matrices``.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: from sage.matrix.seymour_decomposition import OneSumNode

            sage: M2 = Matrix_cmr_chr_sparse.one_sum([[1, 0], [-1, 1]],
            ....:                                    [[1, 1], [-1, 0]])
            sage: result, certificate = M2.is_totally_unimodular(certificate=True); certificate
            OneSumNode (4×4) with 2 children
            sage: OneSumNode.check(M2,
            ....:                  certificate.summand_matrices(),
            ....:                  certificate.child_indices())

        Symbolic identities::

            sage: from sage.matrix.seymour_decomposition import OneSumNode
            sage: R.<x,y> = QQ[]
            sage: A = matrix([[x, 0], [-x, 1]])
            sage: B = matrix([[x, y], [-x, 0]])
            sage: A1B = block_diagonal_matrix([A, B])
            sage: OneSumNode.check(A1B, [A, B], [([0, 1], [0, 1]),
            ....:                                ([2, 3], [2, 3])])

        Using program analysis::

            sage: # optional - cutgeneratingfunctionology
            sage: R.<x,y,z> = ParametricRealField({x: 1}, {y: -1}, {z: 0})  # true example
            sage: A = matrix([[x, 0], [-x, 1]])
            sage: B = matrix([[x, y], [-x, 0]])
            sage: A1B = matrix([[z, 0, 0, 0], [-x, z, 0, 0], [], []])
            sage: OneSumNode.check(A1B, [A, B], [([0, 1], [0, 1]),
            ....:                                ([2, 3], [2, 3])])
            sage: # side-effect: R stores polynomial identities
        """
        # TODO: Check that summand_parent_rows_and_columns form partitions of rows and columns
        for matrix, rows_and_columns in zip(summand_matrices, summand_parent_rows_and_columns):
            assert result_matrix.matrix_from_rows_and_columns(*rows_and_columns) == matrix
        # TODO: Check zero blocks


cdef class TwoSumNode(SumNode):

    def block_matrix_form(self):
        r"""
        Return the block matrix representing the two sum node.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M2 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
            ....:                            [[1, 1, 1, 1, 1], [1, 1, 1, 0, 0],
            ....:                             [1, 0, 1, 1, 0], [1, 0, 0, 1, 1],
            ....:                             [1, 1, 0, 0, 1]]); M2
            [1 1 1 1 1]
            [1 1 1 0 0]
            [1 0 1 1 0]
            [1 0 0 1 1]
            [1 1 0 0 1]
            sage: M3 = Matrix_cmr_chr_sparse.two_sum(M2, M2, 0, 1); M3
            [1 1 1 1|1 1 1 0 0]
            [1 1 0 0|1 1 1 0 0]
            [0 1 1 0|1 1 1 0 0]
            [0 0 1 1|1 1 1 0 0]
            [1 0 0 1|1 1 1 0 0]
            [-------+---------]
            [0 0 0 0|1 1 1 1 1]
            [0 0 0 0|1 0 1 1 0]
            [0 0 0 0|1 0 0 1 1]
            [0 0 0 0|1 1 0 0 1]
            sage: result, certificate = M3.is_totally_unimodular(certificate=True); certificate
            TwoSumNode (9×9) with 2 children

            sage: K33 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 4, sparse=True),
            ....:                            [[1, 1, 0, 0], [1, 1, 1, 0],
            ....:                             [1, 0, 0,-1], [0, 1, 1, 1],
            ....:                             [0, 0, 1, 1]]); K33
            [ 1  1  0  0]
            [ 1  1  1  0]
            [ 1  0  0 -1]
            [ 0  1  1  1]
            [ 0  0  1  1]
            sage: K33_dual = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 4, 5, sparse=True),
            ....:                            [[1, 1, 1, 0, 0], [1, 1, 0, 1, 0],
            ....:                             [0, 1, 0, 1, 1], [0, 0,-1, 1, 1]]); K33_dual
            [ 1  1  1  0  0]
            [ 1  1  0  1  0]
            [ 0  1  0  1  1]
            [ 0  0 -1  1  1]
            sage: M = Matrix_cmr_chr_sparse.two_sum(K33, K33_dual, 0, 0,
            ....:                                   nonzero_block="bottom_left"); M
            [ 1  1  1  0| 0  0  0  0]
            [ 1  0  0 -1| 0  0  0  0]
            [ 0  1  1  1| 0  0  0  0]
            [ 0  0  1  1| 0  0  0  0]
            [-----------+-----------]
            [ 1  1  0  0| 1  1  0  0]
            [ 1  1  0  0| 1  0  1  0]
            [ 0  0  0  0| 1  0  1  1]
            [ 0  0  0  0| 0 -1  1  1]
            sage: result1, certificate1 = M.is_totally_unimodular(certificate=True); certificate1
            TwoSumNode (8×8) with 2 children
            sage: certificate1.summand_matrices()
            (
            [ 1  1  1  0]
            [ 1  0  0 -1]  [ 1  1  1  0  0]
            [ 0  1  1  1]  [ 1  1  0  1  0]
            [ 0  0  1  1]  [ 0  1  0  1  1]
            [ 1  1  0  0], [ 0  0 -1  1  1]
            )
            sage: certificate1.block_matrix_form()
            [ 1  1  1  0| 0  0  0  0]
            [ 1  0  0 -1| 0  0  0  0]
            [ 0  1  1  1| 0  0  0  0]
            [ 0  0  1  1| 0  0  0  0]
            [-----------+-----------]
            [ 1  1  0  0| 1  1  0  0]
            [ 1  1  0  0| 1  0  1  0]
            [ 0  0  0  0| 1  0  1  1]
            [ 0  0  0  0| 0 -1  1  1]
            sage: certificate1.child_indices()
            (((0, 1, 2, 3, 4), (0, 1, 2, 3)), ((4, 5, 6, 7), (0, 4, 5, 6, 7)))
            sage: M_perm = M.matrix_from_rows_and_columns([4, 6, 5, 7, 0, 1, 2, 3], range(M.ncols()))
            sage: M_perm
            [ 1  1  0  0  1  1  0  0]
            [ 0  0  0  0  1  0  1  1]
            [ 1  1  0  0  1  0  1  0]
            [ 0  0  0  0  0 -1  1  1]
            [ 1  1  1  0  0  0  0  0]
            [ 1  0  0 -1  0  0  0  0]
            [ 0  1  1  1  0  0  0  0]
            [ 0  0  1  1  0  0  0  0]
            sage: result2, certificate2 = M_perm.is_totally_unimodular(certificate=True)
            sage: certificate2.summand_matrices()
            (
            [ 1  1  1  0]
            [ 1  0  0 -1]  [ 1  1  1  0  0]
            [ 0  1  1  1]  [ 0  1  0  1  1]
            [ 0  0  1  1]  [ 1  1  0  1  0]
            [ 1  1  0  0], [ 0  0 -1  1  1]
            )
            sage: certificate2.block_matrix_form()
            [ 1  1  1  0| 0  0  0  0]
            [ 1  0  0 -1| 0  0  0  0]
            [ 0  1  1  1| 0  0  0  0]
            [ 0  0  1  1| 0  0  0  0]
            [-----------+-----------]
            [ 1  1  0  0| 1  1  0  0]
            [ 0  0  0  0| 1  0  1  1]
            [ 1  1  0  0| 1  0  1  0]
            [ 0  0  0  0| 0 -1  1  1]
            sage: certificate2.child_indices()
            (((4, 5, 6, 7, 0), (0, 1, 2, 3)), ((0, 1, 2, 3), (0, 4, 5, 6, 7)))
        """
        M1, M2 = self.summand_matrices()
        return Matrix_cmr_chr_sparse.two_sum(M1, M2, M1.nrows() - 1, 0, "bottom_left")

cdef class ThreeSumNode(SumNode):

    def _children(self):
        r"""
        Return a tuple of the tuples of the two children
        and their row and column keys.

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.Matrix_cmr_chr_sparse.three_sum`

        TESTS:

        This is test ``WideWideR12`` and ``MixedMixedR12`` in CMR's ``test_tu.cpp``::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 6, 6, sparse=True),
            ....: [[1,0,1,1,0,0],[0,1,1,1,0,0],[1,0,1,0,1,1],
            ....: [0,-1,0,-1,1,1],[1,0,1,0,1,0],[0,-1,0,-1,0,1]])
            sage: result, certificate = R12.is_totally_unimodular(certificate=True,
            ....:                           three_sum_strategy="Wide_Wide",
            ....:                           row_keys=range(6),
            ....:                           column_keys='abcdef')
            sage: certificate.child_indices()
            ((0, 1, 2, 3, a, 5), (4, b, c, d, e, f))
            sage: C = certificate.child_nodes()[0]
            sage: C1, C2 = C.child_nodes()
            sage: C1.matrix()
            [ 0  0  1 -1 -1]
            [ 1  1  1  0  0]
            [ 0  1  0  1  1]
            [-1  0 -1  0  1]
            sage: C2.matrix()
            [ 1  0  1 -1  0]
            [ 0  0  1  0  1]
            [-1 -1  0  1  1]
            [-1 -1  0  0  1]
            sage: C.child_indices()
            (((0, 1, a, 3), (b, c, d, e, +3+e)), ((0, 2, 3, 5), (+0+d, d, 4, e, f)))
            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: node = UnknownNode(R12,
            ....:                    row_keys=range(6),
            ....:                    column_keys='abcdef'); node
            UnknownNode (6×6)
            sage: C0 = node.complete_decomposition(
            ....:                            three_sum_strategy="Wide_Wide",
            ....:                            )
            sage: C0
            PivotsNode (6×6)
            sage: unicode_art(C0)
                    PivotsNode (6×6)
                    │
            ╭─────────────ThreeSumNode (6×6) with 2 children
            │                   │
            CographicNode (4×5) GraphicNode (4×5)
            sage: unicode_art(node)
            UnknownNode (6×6)

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12_large = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 12, sparse=True),
            ....: [[1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            ....: [0, 0, 0, 1, -1, 0, 0, 0, 1 , 1, 1, 1],
            ....: [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
            ....: [ 1,  0,  1,  0,  0,  0,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  1,  1,  0,  0,  0,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0],
            ....: [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1]])
            sage: result, certificate = R12_large.is_totally_unimodular(certificate=True,
            ....:                                 three_sum_strategy="Wide_Wide",
            ....:                                 row_keys=range(9),
            ....:                                 column_keys='abcdefghijkl')
            sage: C = certificate.child_nodes()[0]; C
            ThreeSumNode (9×12) with 2 children
            sage: C1, C2 = C.child_nodes()
            sage: C1.matrix()
            [ 0  0  1  1  1  1  1]
            [ 1  1  0  0  0 -1 -1]
            [ 1  0 -1  0 -1 -1 -1]
            [ 0  1  1  0  1  0  0]
            [ 0  0  0 -1 -1  0 -1]
            sage: C2.matrix()
            [ 1  0  0  0  0  1 -1  0 -1]
            [ 0  0  1 -1  0 -1  1  0  1]
            [-1 -1  1  0  1 -1  1  0  1]
            [-1 -1  0  1  1  0  0  0  0]
            [-1 -1  0  0  0  0  1  1  1]
            [-1 -1  0  0  0  0  1  1  0]
            sage: C.row_keys()
            (0, i, 2, 3, 4, 5, 6, 7, 8)
            sage: C.column_keys()
            (a, b, c, d, e, f, g, h, 1, j, k, l)
            sage: C.child_indices()[0]
            ((i, 2, 7, 8, 3), (g, h, j, k, l, d, -3+d))
            sage: C.child_indices()[1]
            ((i, 0, 3, 4, 5, 6), (+i+k, k, a, b, c, d, e, f, 1))

            sage: result, certificate = R12.is_totally_unimodular(certificate=True,
            ....:                           three_sum_strategy="Mixed_Mixed")
            sage: C1, C2 = certificate.child_nodes()
            sage: C1.matrix()
            [ 1  0  1  1  0]
            [ 0  1  1  1  0]
            [ 1  0  1  0  1]
            [ 0 -1  0 -1  1]
            sage: C2.matrix()
            [ 1  1  0  0]
            [ 1  0  1  1]
            [ 0 -1  1  1]
            [ 1  0  1  0]
            [ 0 -1  0  1]
            sage: certificate.child_indices()[0]
            ((r0, r1, r2, r3), (c0, c1, c2, c3, +r2+r3))
            sage: certificate.child_indices()[1]
            ((+c0+c3, r2, r3, r4, r5), (c0, c3, c4, c5))

            sage: result, certificate = R12.is_totally_unimodular(certificate=True,
            ....:                           three_sum_strategy="Mixed_Mixed",
            ....:                           row_keys=range(6),
            ....:                           column_keys='abcdef')
            sage: C1, C2 = certificate.child_nodes()
            sage: C1.matrix()
            [ 1  0  1  1  0]
            [ 0  1  1  1  0]
            [ 1  0  1  0  1]
            [ 0 -1  0 -1  1]
            sage: C2.matrix()
            [ 1  1  0  0]
            [ 1  0  1  1]
            [ 0 -1  1  1]
            [ 1  0  1  0]
            [ 0 -1  0  1]
            sage: certificate.child_indices()[0]
            ((0, 1, 2, 3), (a, b, c, d, +2+3))
            sage: certificate.child_indices()[1]
            ((+a+d, 2, 3, 4, 5), (a, d, e, f))
        """
        if self._child_nodes is not None:
            return self._child_nodes

        if self.nchildren() != 2:
            raise ValueError(f"ThreeSumNode has exactly two children not {self.nchildren()}!")

        self.set_default_keys()

        cdef CMR_SEYMOUR_NODE *child1_dec = CMRseymourChild(self._dec, 0)
        cdef CMR_ELEMENT *parent_rows1 = CMRseymourChildRowsToParent(self._dec, 0)
        cdef CMR_ELEMENT *parent_columns1 = CMRseymourChildColumnsToParent(self._dec, 0)
        cdef CMR_CHRMAT *mat1 = CMRseymourGetMatrix(child1_dec)

        cdef CMR_SEYMOUR_NODE *child2_dec = CMRseymourChild(self._dec, 1)
        cdef CMR_ELEMENT *parent_rows2 = CMRseymourChildRowsToParent(self._dec, 1)
        cdef CMR_ELEMENT *parent_columns2 = CMRseymourChildColumnsToParent(self._dec, 1)
        cdef CMR_CHRMAT *mat2 = CMRseymourGetMatrix(child2_dec)

        cdef size_t index1, index2

        child1_nrows = CMRseymourNumRows(child1_dec)
        child1_ncols = CMRseymourNumColumns(child1_dec)

        if self.is_concentrated_rank(): # Mixed_Mixed
            child1_row_keys = tuple(self._CMRelement_to_key(parent_rows1[i])
                                    for i in range(child1_nrows))
            child1_column_keys = tuple(self._CMRelement_to_key(parent_columns1[i])
                                    for i in range(child1_ncols - 1))

            row1_index = child1_nrows - 2
            CMR_CALL(CMRchrmatFindEntry(mat1, row1_index, child1_ncols-1, &index1))
            if index1 == SIZE_MAX:
                eps1 = Integer(0)
            else:
                eps1 = Integer(mat1.entryValues[index1])
            if eps1 != 1:
                raise ValueError(f"First child in the Mixed_Mixed Three Sum "
                                 f"has 1 in the entry  "
                                 f"row {row1_index} and column {child1_ncols-1} "
                                 f"but got {eps1}")

            row2_index = child1_nrows - 1
            CMR_CALL(CMRchrmatFindEntry(mat1, row2_index, child1_ncols-1, &index2))
            if index2 == SIZE_MAX:
                eps2 = Integer(0)
            else:
                eps2 = Integer(mat1.entryValues[index2])
            if eps2 != 1 and eps2 != -1:
                raise ValueError(f"First child in the Mixed_Mixed Three Sum "
                                 f"has 1 or -1 in the entry  "
                                 f"row {row2_index} and column {child1_ncols-1} "
                                 f"but got {eps2}")

            extra_key = ElementKey((eps1, child1_row_keys[row1_index],
                                    eps2, child1_row_keys[row2_index]),
                                    composition=True)
            child1_column_keys += (extra_key,)
        else: # Wide_Wide
            child1_row_keys = tuple(self._CMRelement_to_key(parent_rows1[i])
                                    for i in range(child1_nrows))
            child1_column_keys = tuple(self._CMRelement_to_key(parent_columns1[i])
                                    for i in range(child1_ncols - 1))

            row_index = child1_nrows - 1
            column_index = child1_ncols - 2
            CMR_CALL(CMRchrmatFindEntry(mat1, row_index, child1_ncols-1, &index1))
            if index1 == SIZE_MAX:
                eps1 = Integer(0)
            else:
                eps1 = Integer(mat1.entryValues[index1])
            if eps1 != 1 and eps1 != -1:
                raise ValueError(f"First child in the Wide_Wide Three Sum "
                                 f"has 1 or -1 in the entry  "
                                 f"row {row_index} and column {child1_ncols-1} "
                                 f"but got {eps1}")

            extra_key = ElementKey((1, child1_column_keys[column_index],
                                    eps1, child1_row_keys[row_index]),
                                    composition=True)
            child1_column_keys += (extra_key,)

        child1 = create_DecompositionNode(child1_dec, matrix=None,
                                          row_keys=child1_row_keys,
                                          column_keys=child1_column_keys,
                                          base_ring=self.base_ring())

        child2_nrows = CMRseymourNumRows(child2_dec)
        child2_ncols = CMRseymourNumColumns(child2_dec)

        if self.is_concentrated_rank(): # Mixed_Mixed
            child2_row_keys = tuple(self._CMRelement_to_key(parent_rows2[i])
                                    for i in range(1, child2_nrows))
            child2_column_keys = tuple(self._CMRelement_to_key(parent_columns2[i])
                                       for i in range(child2_ncols))

            CMR_CALL(CMRchrmatFindEntry(mat2, 0, 0, &index1))
            if index1 == SIZE_MAX:
                eps1 = Integer(0)
            else:
                eps1 = Integer(mat1.entryValues[index1])
            if eps1 != 1 and eps1 != -1:
                raise ValueError(f"Second child in the Mixed_Mixed Three Sum "
                                 f"has 1 or -1 in the entry  "
                                 f"row {0} and column {0} "
                                 f"but got {eps1}")

            CMR_CALL(CMRchrmatFindEntry(mat2, 0, 1, &index2))
            if index2 == SIZE_MAX:
                eps2 = Integer(0)
            else:
                eps2 = Integer(mat1.entryValues[index2])
            if eps2 != 1:
                raise ValueError(f"Second child in the Mixed_Mixed Three Sum "
                                 f"has 1 in the entry  "
                                 f"row {0} and column {1} "
                                 f"but got {eps2}")

            extra_key = ElementKey((eps1, child2_column_keys[0],
                                    eps2, child2_column_keys[1]),
                                    composition=True)
            child2_row_keys = (extra_key,) + child2_row_keys
        else: # Wide_Wide
            child2_row_keys = tuple(self._CMRelement_to_key(parent_rows2[i])
                                    for i in range(child2_nrows))

            CMR_CALL(CMRchrmatFindEntry(mat2, 0, 0, &index1))
            if index1 == SIZE_MAX:
                eps1 = Integer(0)
            else:
                eps1 = Integer(mat1.entryValues[index1])

            if eps1 != 1 and eps1 != -1:
                raise ValueError(f"Second child in the Wide_Wide Three Sum "
                                 f"has 1 or -1 in the entry  "
                                 f"row {0} and column {0} "
                                 f"but got {eps1}")

            child2_column_keys = tuple(self._CMRelement_to_key(parent_columns2[i])
                                       for i in range(1, child2_ncols))
            extra_key = ElementKey((1, child2_column_keys[0],
                                    eps1, child2_row_keys[0]),
                                    composition=True)
            child2_column_keys = (extra_key,) + child2_column_keys

        child2 = create_DecompositionNode(child2_dec, matrix=None,
                                          row_keys=child2_row_keys,
                                          column_keys=child2_column_keys,
                                          base_ring=self.base_ring())

        self._child_nodes = ((child1, child1_row_keys, child1_column_keys),
                             (child2, child2_row_keys, child2_column_keys))
        return self._child_nodes

    def is_distributed_ranks(self):
        r"""
        Check whether the three sum node ``self`` is formed with
        ``three_sum_strategy="distributed_ranks"`` or ``"Wide_Wide"``.

        The matrix representing the first child is
        `M_1=\begin{bmatrix} A & a_2 & a_2\\ a_1^T & 0 & \epsilon_2\end{bmatrix}`
        and the matrix representing the second child is
        `M_2=\begin{bmatrix} \epsilon_1 & 0 & b_2^T\\ b_1 & b_1 & B\end{bmatrix}`,
        where `\epsilon_1`, `\epsilon_2` are `1` or `-1`.
        And the matrix representing ``self`` is a permutation of
        `M_1 \oplus_3 M_2 = \begin{bmatrix} A & a_2 b_2^T \\ b_1 a_1^T & B\end{bmatrix}`.

        ``distributed_ranks`` is named after the two rank 1 off-diagonal blocks.
        ``Wide_Wide`` is named after the structure of the two children.

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.Matrix_cmr_chr_sparse.three_sum_wide_wide`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12_large = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 12, sparse=True),
            ....: [[1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            ....: [0, 0, 0, 1, -1, 0, 0, 0, 1 , 1, 1, 1],
            ....: [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
            ....: [ 1,  0,  1,  0,  0,  0,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  1,  1,  0,  0,  0,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0],
            ....: [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1]])
            sage: result, certificate = R12_large.is_totally_unimodular(certificate=True)
            sage: C = certificate.child_nodes()[0]; C
            ThreeSumNode (9×12) with 2 children
            sage: C.is_distributed_ranks()
            True
            sage: C.is_concentrated_rank()
            False
        """
        return <bint> CMRseymourThreeSumDistributedRanks(self._dec)

    def is_concentrated_rank(self):
        r"""
        Check whether the three sum node ``self`` is formed with
        ``three_sum_strategy="concentrated_rank"`` or ``"Mixed_Mixed"``.

        The matrix representing the first child is
        `M_1=\begin{bmatrix} A & 0 \\ a_1^T & 1\\ a_2^T & \epsilon_2\end{bmatrix}`
        and the matrix representing the second child is
        `M_2=\begin{bmatrix} \epsilon_1 & 1 & 0\\ b_1 & b_2 & B\end{bmatrix}`,
        where `\epsilon_1`, `\epsilon_2` are `1` or `-1`.
        And the matrix representing ``self`` is a permutation of
        `M_1 \oplus_3 M_2 = \begin{bmatrix} A & 0 \\ C & B\end{bmatrix}`,
        where `\begin{bmatrix}a_1^T \\ a_2^T\end{bmatrix}=\begin{bmatrix} C_1 & \bar{C}\end{bmatrix}`,
        `\begin{bmatrix}b_1 & b_2\end{bmatrix}=\begin{bmatrix} \bar{C}\\ C_2\end{bmatrix}`,
        `C_12 = C_2 \bar{C}^{-1} C_1`,
        `C=\begin{bmatrix} C_1 & \bar{C} \\ C_{12} & C_2\end{bmatrix}`, i.e.,
        `C=\begin{bmatrix}b_1 & b_2\end{bmatrix}\bar{C}^{-1}\begin{bmatrix}a_1^T\\ a_2^T\end{bmatrix}`

        ``concentrated_rank`` is named after the rank 2 off-diagonal block.
        ``Mixed_Mixed`` is named after the structure of the two children.

        .. SEEALSO::

            :meth:`sage.matrix.matrix_cmr_sparse.Matrix_cmr_chr_sparse.three_sum_mixed_mixed`

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12_large = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 12, sparse=True),
            ....: [[1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            ....: [0, 0, 0, 1, -1, 0, 0, 0, 1 , 1, 1, 1],
            ....: [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
            ....: [ 1,  0,  1,  0,  0,  0,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  1,  1,  0,  0,  0,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0],
            ....: [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1]])
            sage: result, certificate = R12_large.is_totally_unimodular(certificate=True,
            ....:                                 three_sum_strategy="Mixed_Mixed")
            sage: C = certificate; C
            ThreeSumNode (9×12) with 2 children
            sage: C.is_distributed_ranks()
            False
            sage: C.is_concentrated_rank()
            True
        """
        return <bint> CMRseymourThreeSumConcentratedRank(self._dec)

    def block_matrix_form(self):
        r"""
        Return the block matrix constructed from the three sum of children.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 6, 6, sparse=True),
            ....: [[1,0,1,1,0,0],[0,1,1,1,0,0],[1,0,1,0,1,1],
            ....: [0,-1,0,-1,1,1],[1,0,1,0,1,0],[0,-1,0,-1,0,1]])
            sage: R12
            [ 1  0  1  1  0  0]
            [ 0  1  1  1  0  0]
            [ 1  0  1  0  1  1]
            [ 0 -1  0 -1  1  1]
            [ 1  0  1  0  1  0]
            [ 0 -1  0 -1  0  1]
            sage: result, certificate = R12.is_totally_unimodular(certificate=True)
            sage: C = certificate.child_nodes()[0]; C
            ThreeSumNode (6×6) with 2 children
            sage: C.matrix()
            [ 1  0  0  1 -1  0]
            [ 0  1  1  1  0  0]
            [ 1  0  0  0  0  1]
            [ 0 -1  0 -1  1  1]
            [-1  0  1  0  1  0]
            [ 0 -1  0 -1  0  1]
            sage: C.summand_matrices()
            (
            [ 0  0  1 -1 -1]  [ 1  0  1 -1  0]
            [ 1  1  1  0  0]  [ 0  0  1  0  1]
            [ 0  1  0  1  1]  [-1 -1  0  1  1]
            [-1  0 -1  0  1], [-1 -1  0  0  1]
            )
            sage: C.block_matrix_form()
            [ 0  0  1 -1  1  0]
            [ 1  1  1  0  0  0]
            [ 0  1  0  1 -1  0]
            [ 0  0  0  1  0  1]
            [ 1  0  1  0  1  1]
            [ 1  0  1  0  0  1]

            sage: result, certificate = R12.is_totally_unimodular(certificate=True,
            ....:                           three_sum_strategy="Mixed_Mixed")
            sage: C = certificate; C
            ThreeSumNode (6×6) with 2 children
            sage: C.matrix()
            [ 1  0  1  1  0  0]
            [ 0  1  1  1  0  0]
            [ 1  0  1  0  1  1]
            [ 0 -1  0 -1  1  1]
            [ 1  0  1  0  1  0]
            [ 0 -1  0 -1  0  1]
            sage: C.summand_matrices()
            (
            [ 1  1  0  0]
            [ 1  0  1  1  0]  [ 1  0  1  1]
            [ 0  1  1  1  0]  [ 0 -1  1  1]
            [ 1  0  1  0  1]  [ 1  0  1  0]
            [ 0 -1  0 -1  1], [ 0 -1  0  1]
            )
            sage: C.child_indices()
            (((r0, r1, r2, r3), (c0, c1, c2, c3, +r2+r3)),
            ((+c0+c3, r2, r3, r4, r5), (c0, c3, c4, c5)))
            sage: C.block_matrix_form()
            [ 1  0  1  1  0  0]
            [ 0  1  1  1  0  0]
            [ 1  0  1  0  1  1]
            [ 0 -1  0 -1  1  1]
            [ 1  0  1  0  1  0]
            [ 0 -1  0 -1  0  1]
        """
        M1, M2 = self.summand_matrices()
        if self.is_distributed_ranks():
            three_sum_strategy = 'distributed_ranks'
        else:
            three_sum_strategy = 'concentrated_rank'
        return Matrix_cmr_chr_sparse.three_sum(M1, M2,
                                               three_sum_strategy=three_sum_strategy,
                                               sign_verify=True)


cdef class BaseGraphicNode(DecompositionNode):

    def __init__(self, matrix=None,
                 graph=None, forest_edges=None, coforest_edges=None,
                 row_keys=None, column_keys=None, base_ring=None):
        super().__init__(matrix=matrix, row_keys=row_keys, column_keys=column_keys,
                         base_ring=base_ring)
        self._graph = graph
        self._forest_edges = forest_edges
        self._coforest_edges = coforest_edges

    def graph(self):
        r"""
        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: G = certificate.graph(); G
            Graph on 4 vertices
            sage: G.vertices(sort=True)
            [1, 2, 7, 12]
            sage: G.edges(sort=True)
            [(1, 2, None), (1, 7, None), (1, 12, None), (2, 7, None), (7, 12, None)]
        """
        if self._graph is not None:
            return self._graph
        self._graph = _sage_graph(CMRseymourGraph(self._dec))
        return self._graph

    def forest_edges(self):
        r"""
        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, GraphicNode (3×2))
            sage: certificate.forest_edges()
            ((1, 2), (7, 1), (12, 7))

        Starting with a morphism::

            sage: from sage.matrix.seymour_decomposition import UnknownNode
            sage: phi = matrix(ZZ, [[1, 0], [1, 1], [0, 1]],
            ....:              row_keys=['a', 'b', 'c'], column_keys=['v', 'w'])
            sage: phi; phi._unicode_art_matrix()
            Generic morphism:
            From: Free module generated by {'v', 'w'} over Integer Ring
            To:   Free module generated by {'a', 'b', 'c'} over Integer Ring
              v w
            a⎛1 0⎞
            b⎜1 1⎟
            c⎝0 1⎠
            sage: phi_node = UnknownNode(phi)
            sage: is_graphic, rephined_node = phi_node._is_binary_linear_matroid_graphic(decomposition=True)
            sage: is_graphic, rephined_node
            (True, GraphicNode (3×2))
            sage: rephined_node.forest_edges()
            {'a': (1, 2), 'b': (7, 1), 'c': (12, 7)}
            sage: phi_node  # still in the dark about graphicness
            UnknownNode (3×2)
        """
        if self._forest_edges is not None:
            return self._forest_edges
        cdef CMR_GRAPH *graph = CMRseymourGraph(self._dec)
        cdef size_t num_edges = CMRseymourGraphSizeForest(self._dec)
        cdef CMR_GRAPH_EDGE *edges = CMRseymourGraphForest(self._dec)
        self._forest_edges = _sage_edges(graph, edges, num_edges, self.row_keys())
        return self._forest_edges

    def coforest_edges(self):
        if self._coforest_edges is not None:
            return self._coforest_edges
        cdef CMR_GRAPH *graph = CMRseymourGraph(self._dec)
        cdef size_t num_edges = CMRseymourGraphSizeCoforest(self._dec)
        cdef CMR_GRAPH_EDGE *edges = CMRseymourGraphCoforest(self._dec)
        self._coforest_edges = _sage_edges(graph, edges, num_edges, self.column_keys())
        return self._coforest_edges


cdef class GraphicNode(BaseGraphicNode):

    pass


cdef class CographicNode(BaseGraphicNode):
    @cached_method
    def graph(self):
        r"""
        Actually the cograph of matrix, in the case where it is not graphic.
        """
        return _sage_graph(CMRseymourCograph(self._dec))


cdef class PlanarNode(BaseGraphicNode):
    @cached_method
    def cograph(self):
        r"""
        Return the cograph of matrix.
        """
        return _sage_graph(CMRseymourCograph(self._dec))


cdef class SeriesParallelReductionNode(DecompositionNode):

    def core(self):
        r"""
        Return the core of ``self``.

        A :class:SeriesParallelReductionNode indicates that `M`
        arises from a smaller matrix `M'` (called the core)
        by successively adding zero rows/columns,
        unit rows/columns or duplicates of existing rows/columns
        (potentially scaled with `-1`).

        Note that such series-parallel reductions preserve total unimodularity
        and binary regularity.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 6, sparse=True),
            ....:                           [[1, 1, 1, 1, 1, 0], [1, 1, 1, 0, 0, 0],
            ....:                            [1, 0, 1, 1, 0, 1] ,[1, 0, 0, 1, 1, 0],
            ....:                            [1, 1, 0, 0, 1, 0]]); M
            [1 1 1 1 1 0]
            [1 1 1 0 0 0]
            [1 0 1 1 0 1]
            [1 0 0 1 1 0]
            [1 1 0 0 1 0]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: result, certificate
            (True, SeriesParallelReductionNode (5×6))
            sage: certificate.core()
            [1 1 1 1 1]
            [1 1 1 0 0]
            [1 0 1 1 0]
            [1 0 0 1 1]
            [1 1 0 0 1]
        """
        return self.child_nodes()[0].matrix()


cdef class R10Node(DecompositionNode):
    r"""
    Special R10 Node.
    Only two possible 5 by 5 matrices up to row and column permutations and negations.

    EXAMPLES::

        sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
        sage: R10 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
        ....:                            [[1, 0, 0, 1, 1],
        ....:                             [1, 1, 0, 0, 1],
        ....:                             [0, 1, 1, 0, 1],
        ....:                             [0, 0, 1, 1, 1],
        ....:                             [1, 1, 1, 1, 1]]); R10
        [1 0 0 1 1]
        [1 1 0 0 1]
        [0 1 1 0 1]
        [0 0 1 1 1]
        [1 1 1 1 1]
        sage: result, certificate = R10.is_totally_unimodular(certificate=True)
        sage: result
        True
        sage: R10.is_network_matrix()
        False
        sage: R10.is_conetwork_matrix()
        False
        sage: result, certificate = R10._is_binary_linear_matroid_regular(certificate=True)
        sage: result
        True
        sage: certificate._is_binary_linear_matroid_graphic()
        False
        sage: certificate._is_binary_linear_matroid_cographic()
        False

        sage: R10 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
        ....:                            [[ 1,-1, 0, 0,-1],
        ....:                             [-1, 1,-1, 0, 0],
        ....:                             [ 0,-1, 1,-1, 0],
        ....:                             [ 0, 0,-1, 1,-1],
        ....:                             [-1, 0, 0,-1, 1]]); R10
        [ 1 -1  0  0 -1]
        [-1  1 -1  0  0]
        [ 0 -1  1 -1  0]
        [ 0  0 -1  1 -1]
        [-1  0  0 -1  1]
        sage: result, certificate = R10.is_totally_unimodular(certificate=True)
        sage: result
        True
        sage: R10.is_network_matrix()
        False
        sage: R10.is_conetwork_matrix()
        False

        sage: R10 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
        ....:                            [[1, 1, 0, 0, 1],
        ....:                             [1, 1,-1, 0, 0],
        ....:                             [0, 1,-1,-1, 0],
        ....:                             [0, 0, 1, 1, 1],
        ....:                             [1, 0, 0, 1, 1]]); R10
        [ 1  1  0  0  1]
        [ 1  1 -1  0  0]
        [ 0  1 -1 -1  0]
        [ 0  0  1  1  1]
        [ 1  0  0  1  1]
        sage: result, certificate = R10.is_totally_unimodular(certificate=True)
        sage: result
        True
        sage: R10.is_network_matrix()
        False
        sage: R10.is_conetwork_matrix()
        False

        sage: R10 = Matrix_cmr_chr_sparse(MatrixSpace(GF(2), 5, 5, sparse=True),
        ....:                            [[1, 1, 0, 0, 1],
        ....:                             [1, 1,-1, 0, 0],
        ....:                             [0, 1,-1,-1, 0],
        ....:                             [0, 0, 1, 1, 1],
        ....:                             [1, 0, 0, 1, 1]]); R10
        [1 1 0 0 1]
        [1 1 1 0 0]
        [0 1 1 1 0]
        [0 0 1 1 1]
        [1 0 0 1 1]
        sage: result, certificate = R10._is_binary_linear_matroid_regular(certificate=True)
        sage: result
        True
        sage: certificate
        R10Node (5×5) a reduced matrix representation of R10 matroid
        sage: certificate._is_binary_linear_matroid_graphic()
        False
        sage: certificate._is_binary_linear_matroid_cographic()
        False
    """

    @cached_method
    def _matroid(self):
        r"""
        Return the R10 matroid represented by ``self``.
        ``self.matrix()`` is the reduced matrix representation of the matroid.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R10 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
            ....:                            [[1, 1, 1, 1, 1],
            ....:                             [1, 1, 1, 0, 0],
            ....:                             [1, 0, 1, 1, 0],
            ....:                             [1, 0, 0, 1, 1],
            ....:                             [1, 1, 0, 0, 1]]); R10
            [1 1 1 1 1]
            [1 1 1 0 0]
            [1 0 1 1 0]
            [1 0 0 1 1]
            [1 1 0 0 1]
            sage: result, certificate = R10.is_totally_unimodular(certificate=True,
            ....:                                                 row_keys=range(5),
            ....:                                                 column_keys='abcde')
            sage: certificate._matroid()
            R10: Regular matroid of rank 5 on 10 elements with 162 bases
            sage: certificate._isomorphism()
            {'a': 0,
            'b': 1,
            'c': 2,
            'd': 'a',
            'e': 'b',
            'f': 4,
            'g': 'c',
            'h': 'e',
            'i': 3,
            'j': 'd'}
            sage: result, certificate = R10.is_totally_unimodular(certificate=True)
            sage: certificate._matroid()
            R10: Regular matroid of rank 5 on 10 elements with 162 bases
            sage: certificate._isomorphism()
            {'a': 0,
            'b': 1,
            'c': 2,
            'd': 5,
            'e': 6,
            'f': 4,
            'g': 7,
            'h': 9,
            'i': 3,
            'j': 8}

            sage: R10 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 5, 5, sparse=True),
            ....:                            [[ 1,-1, 0, 0,-1],
            ....:                             [-1, 1,-1, 0, 0],
            ....:                             [ 0,-1, 1,-1, 0],
            ....:                             [ 0, 0,-1, 1,-1],
            ....:                             [-1, 0, 0,-1, 1]]); R10
            [ 1 -1  0  0 -1]
            [-1  1 -1  0  0]
            [ 0 -1  1 -1  0]
            [ 0  0 -1  1 -1]
            [-1  0  0 -1  1]
            sage: result, certificate = R10.is_totally_unimodular(certificate=True,
            ....:                                                 row_keys=range(5),
            ....:                                                 column_keys='abcde')
            sage: certificate._matroid()
            R10: Regular matroid of rank 5 on 10 elements with 162 bases
            sage: certificate._isomorphism() # The isomorphism is not unique
            {'a': 0,
            'b': 1,
            'c': 2,
            'd': 'c',
            'e': 'a',
            'f': 4,
            'g': 'b',
            'h': 3,
            'i': 'e',
            'j': 'd'}
        """
        from sage.matroids.constructor import Matroid
        if self.row_keys() is None or self.column_keys() is None:
            M = Matroid(reduced_matrix=self.matrix(), regular=True)
        else:
            M = Matroid(reduced_morphism=self.morphism(), regular=True)
        M.rename("R10: " + repr(M))
        return M

    def _isomorphism(self):
        r"""
        Return one isomorphism between the R10 matroid
        and ``self._matroid()``.
        """
        import sage.matroids.matroids_catalog as matroids
        M0 = matroids.catalog.R10()
        result, isomorphism = M0.is_isomorphic(self._matroid(), certificate=True)
        if result is False:
            raise ValueError("This is not a R10 node")
        else:
            return isomorphism

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        result = super()._repr_()
        result += f' a reduced matrix representation of R10 matroid'
        return result


cdef class PivotsNode(DecompositionNode):

    def npivots(self):
        r"""
        Return the number of pivots in ``self``.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 12, sparse=True),
            ....: [[1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            ....: [0, 0, 0, 1, -1, 0, 0, 0, 1 , 1, 1, 1],
            ....: [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
            ....: [ 1,  0,  1,  0,  0,  0,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  1,  1,  0,  0,  0,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0],
            ....: [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1]])
            sage: result, certificate = R12.is_totally_unimodular(certificate=True)
            sage: certificate
            PivotsNode (9×12)
            sage: certificate.npivots()
            1
        """
        return CMRseymourNumPivots(self._dec)

    @cached_method
    def pivot_rows_and_columns(self):
        r"""
        Return a tuple of the row and column indices of all pivot entries.

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 12, sparse=True),
            ....: [[1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            ....: [0, 0, 0, 1, -1, 0, 0, 0, 1 , 1, 1, 1],
            ....: [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
            ....: [ 1,  0,  1,  0,  0,  0,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  1,  1,  0,  0,  0,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0],
            ....: [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1]])
            sage: result, certificate = R12.is_totally_unimodular(certificate=True)
            sage: certificate
            PivotsNode (9×12)
            sage: certificate.pivot_rows_and_columns()
            ((1, 8),)
        """
        cdef size_t *pivot_rows = CMRseymourPivotRows(self._dec)
        cdef size_t *pivot_columns = CMRseymourPivotColumns(self._dec)

        return tuple((pivot_rows[i], pivot_columns[i]) for i in range(self.npivots()))

    def _children(self):
        r"""
        Return a tuple of the tuples of children and their row and column keys.
        The underlying implementation of :meth:`child_nodes`
        and :meth:`child_indices`.

        If row and column keys are not given, set the default keys.
        See alse :meth:set_default_keys

        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: R12 = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 9, 12, sparse=True),
            ....: [[1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            ....: [0, 0, 0, 1, -1, 0, 0, 0, 1 , 1, 1, 1],
            ....: [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
            ....: [ 1,  0,  1,  0,  0,  0,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  1,  1,  0,  0,  0,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0],
            ....: [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1],
            ....: [ 0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0],
            ....: [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1]])
            sage: result, certificate = R12.is_totally_unimodular(certificate=True)
            sage: certificate
            PivotsNode (9×12)
            sage: certificate.row_keys()
            sage: certificate.child_nodes()
            (ThreeSumNode (9×12) with 2 children,)
            sage: certificate.row_keys()
            (r0, r1, r2, r3, r4, r5, r6, r7, r8)
        """
        if self._child_nodes is not None:
            return self._child_nodes
        self.set_default_keys()
        children_tuple = tuple(self._create_child_node(index)
                               for index in range(self.nchildren()))
        self._child_nodes = children_tuple
        return self._child_nodes


cdef class SymbolicNode(DecompositionNode):

    def __init__(self, symbol, *, row_keys=None, column_keys=None, base_ring=None):
        r"""
        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import SymbolicNode
            sage: X = SymbolicNode('X', row_keys='abc', column_keys=range(6)); X
            SymbolicNode X (3×6)
            sage: XX = X.one_sum(X)
            Traceback (most recent call last):
            ...
            ValueError: keys must be disjoint...
            sage: XX = X.one_sum(X, summand_ids=(0, 1)); XX
            OneSumNode (6×12) with 2 children
            sage: XX.row_keys()
            ((0, 'a'), (0, 'b'), (0, 'c'), (1, 'a'), (1, 'b'), (1, 'c'))
            sage: T = XX.as_ordered_tree(); T
            OneSumNode (6×12) with 2 children[SymbolicNode X (3×6)[],
                                              SymbolicNode X (3×6)[]]
            sage: unicode_art(T)
            ╭────────────────OneSumNode (6×12) with 2 children
            │                    │
            SymbolicNode X (3×6) SymbolicNode X (3×6)
            sage: Y = SymbolicNode('Y', row_keys='de', column_keys='fg'); Y
            SymbolicNode Y (2×2)
            sage: XY = X.one_sum(Y); XY
            OneSumNode (5×8) with 2 children
        """
        super().__init__(row_keys=row_keys, column_keys=column_keys, base_ring=base_ring)
        self._symbol = symbol

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matrix.seymour_decomposition import SymbolicNode
            sage: X = SymbolicNode('X', row_keys='abc', column_keys=range(6))
            sage: print(X)
            SymbolicNode X (3×6)
        """
        nrows, ncols = self.dimensions()
        symbol = self.symbol()
        return f'{self.__class__.__name__} {symbol} ({nrows}×{ncols})'

    def matrix(self):
        raise ValueError('symbolic nodes are not backed by CMR matrices')

    def symbol(self):
        return self._symbol


cdef class ElementKey:

    cdef frozenset _key
    cdef bint _composition

    def __init__(self, keys, composition=False):
        r"""
        Create the element key for a row or column index
        of :class:`DecompositionNode`.

        INPUT:

        - ``keys`` -- a row/column key or a tuple
          (`\pm 1`, row/column key, `\pm 1`, row/column key).

        - ``composition`` -- ``True`` or ``False`` (default).
          whether the key is a composition key or not.
          If ``False``, ``self._key`` is a frozenset with a row/column key.
          If ``True``, ``self._key`` is a frozenset with two tuples,
          where each tuple has a sign value and a row/column key.
          For example, ``frozenset((1,'a'), (-1,'7'))``.
        """
        if composition:
            sign1, key1, sign2, key2 = keys
            self._key = frozenset([(sign1, key1), (sign2, key2)])
            self._composition = True
        else:
            self._key = frozenset((keys,))
            self._composition = False

    @property
    def key(self):
        return self._key

    def __hash__(self):
        """
        Return a hash of this element key
        """
        return hash(self._key)

    def __eq__(self, other):
        if isinstance(other, ElementKey):
            return self._key == other._key
        return False

    def __repr__(self):
        """
        Return a string representation of ``self``.
        The composition key is sorted by the string of keys.
        """
        if self._composition:
            sorted_key = sorted(self._key, key=lambda x: (str(x[1]), x[0]))
            return "".join(['+'+str(a[1]) if a[0] == 1 else '-'+str(a[1]) for a in sorted_key])
        else:
            return "".join([str(a) for a in self._key])


cdef _class(CMR_SEYMOUR_NODE *dec):
    r"""
    Return the class of the decomposition `CMR_SEYMOUR_NODE`.

    INPUT:

    - ``dec`` -- a ``CMR_SEYMOUR_NODE``
    """
    cdef CMR_SEYMOUR_NODE_TYPE typ = CMRseymourType(dec)

    if typ == CMR_SEYMOUR_NODE_TYPE_ONE_SUM:
        return OneSumNode
    if typ == CMR_SEYMOUR_NODE_TYPE_TWO_SUM:
        return TwoSumNode
    if typ == CMR_SEYMOUR_NODE_TYPE_THREE_SUM:
        return ThreeSumNode
    if typ == CMR_SEYMOUR_NODE_TYPE_GRAPH:
        return GraphicNode
    if typ == CMR_SEYMOUR_NODE_TYPE_COGRAPH:
        return CographicNode
    if typ == CMR_SEYMOUR_NODE_TYPE_PLANAR:
        return PlanarNode
    if typ == CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL:
        return SeriesParallelReductionNode
    if typ == CMR_SEYMOUR_NODE_TYPE_PIVOTS:
        return PivotsNode
    if typ == CMR_SEYMOUR_NODE_TYPE_IRREGULAR:
        return ThreeConnectedIrregularNode
    if typ == CMR_SEYMOUR_NODE_TYPE_UNKNOWN:
        return UnknownNode
    if typ == CMR_SEYMOUR_NODE_TYPE_R10:
        return R10Node
    raise NotImplementedError


cdef create_DecompositionNode(CMR_SEYMOUR_NODE *dec, matrix=None, row_keys=None, column_keys=None, base_ring=None):
    r"""
    Create an instance of a subclass of :class:`DecompositionNode`.

    INPUT:

    - ``dec`` -- a ``CMR_SEYMOUR_NODE``
    """
    if dec == NULL:
        return None
    cdef DecompositionNode result = <DecompositionNode> _class(dec)(
        matrix, row_keys=row_keys, column_keys=column_keys, base_ring=base_ring)
    result._set_dec(dec)
    return result
