r"""
Seymour's decomposition of totally unimodular matrices and regular matroids
"""

from sage.libs.cmr.cmr cimport *
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.structure.sage_object cimport SageObject

from .matrix_cmr_sparse cimport Matrix_cmr_chr_sparse, _sage_graph
from .matrix_space import MatrixSpace


cdef class DecompositionNode(SageObject):

    cdef _set_dec(self, CMR_DEC *dec, root):
        if self._root is None:
            CMR_CALL(CMRdecFree(cmr, &self._dec))
        self._dec = dec
        self._root = root

    def __dealloc__(self):
        self._set_dec(NULL, None)

    @cached_method
    def matrix(self):
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
            (True, GraphicNode with 0 children)
            sage: certificate.matrix() is None
            True

            sage: result, certificate = M.is_totally_unimodular(certificate=True,
            ....:                                               construct_matrices=True)
            sage: result, certificate
            (True, GraphicNode with 0 children)
            sage: certificate.matrix()
            [ 1  0]
            [-1  1]
            [ 0  1]
        """
        cdef Matrix_cmr_chr_sparse result
        cdef CMR_CHRMAT *mat = CMRdecGetMatrix(self._dec)
        if mat == NULL:
            return None
        ms = MatrixSpace(ZZ, mat.numRows, mat.numColumns, sparse=True)
        result = Matrix_cmr_chr_sparse.__new__(Matrix_cmr_chr_sparse, ms)
        result._mat = mat
        result._root = self._root
        return result

    def parent_rows_and_columns(self):
        r"""
        EXAMPLES::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 3, 2, sparse=True),
            ....:                           [[1, 0], [-1, 1], [0, 1]]); M
            [ 1  0]
            [-1  1]
            [ 0  1]
            sage: result, certificate = M.is_totally_unimodular(certificate=True)
            sage: certificate.parent_rows_and_columns()
            (None, None)
        """

        cdef size_t *parent_rows = CMRdecRowsParent(self._dec)
        cdef size_t *parent_columns = CMRdecColumnsParent(self._dec)
        if parent_rows == NULL:
            parent_rows_tuple = None
        else:
            parent_rows_tuple = tuple(parent_rows[i] for i in range(CMRdecNumRows(self._dec)))
        if parent_columns == NULL:
            parent_columns_tuple = None
        else:
            parent_columns_tuple = tuple(parent_columns[i] for i in range(CMRdecNumColumns(self._dec)))

        return parent_rows_tuple, parent_columns_tuple

    def plot(self):
        raise NotImplementedError

    @cached_method
    def _children(self):
        return tuple(create_DecompositionNode(CMRdecChild(self._dec, index),
                                              self._root or self)
                     for index in range(CMRdecNumChildren(self._dec)))

    def _repr_(self):
        return f'{self.__class__.__name__} with {len(self._children())} children'


cdef class SumNode(DecompositionNode):

    def permuted_block_matrix(self):
        r"Return (Prow, BlockMatrix, Pcolumn) so that self.matrix() == Prow * BlockMatrix * Pcolumn ????"
        raise NotImplementedError

    summands = DecompositionNode._children


cdef class OneSumNode(SumNode):

    pass


cdef class TwoSumNode(SumNode):

    pass


cdef class ThreeSumNode(SumNode):

    pass


cdef class BaseGraphicNode(DecompositionNode):

    @cached_method
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
            (True, GraphicNode with 0 children)
            sage: G = certificate.graph(); G
            Graph on 4 vertices
            sage: G.vertices(sort=True)
            [1, 2, 7, 12]
            sage: G.edges(sort=True)
            [(1, 2, None), (1, 7, None), (1, 12, None), (2, 7, None), (7, 12, None)]
        """
        return _sage_graph(CMRdecGraph(self._dec))


cdef class GraphicNode(BaseGraphicNode):

    pass


cdef class CographicNode(BaseGraphicNode):

    pass


cdef class PlanarNode(BaseGraphicNode):

    pass


cdef class SpecialLeafNode(DecompositionNode):


    pass


cdef class K33Node(LeafNode):

    def graph(self):
        raise NotImplementedError

    pass


cdef _class(CMR_DEC *dec):
    k = CMRdecIsSum(dec, NULL, NULL)
    if k == 1:
        return OneSumNode
    if k == 2:
        return TwoSumNode
    if k == 3:
        return ThreeSumNode
    if CMRdecIsGraphicLeaf(dec):
        if CMRdecIsCographicLeaf(dec):
            return PlanarNode
        return GraphicNode
    if CMRdecIsCographicLeaf(dec):
        return CographicNode
    # More TBD
    return DecompositionNode


cdef create_DecompositionNode(CMR_DEC *dec, root=None):
    r"""
    Create an instance of a subclass of :class:`DecompositionNode`.

    INPUT:

    - ``dec`` -- a ``CMR_DEC``
    - ``root`` -- a :class:`DecompositionNode` or ``None``.
      If ``None``, ``dec`` will be owned by the returned instance.
      If non-``None``, ``dec`` is owned by that instance.
    """
    if dec == NULL:
        return None
    cdef DecompositionNode result = <DecompositionNode> _class(dec)()
    result._set_dec(dec, root)
    return result
