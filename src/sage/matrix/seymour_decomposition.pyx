r"""
Seymour's decomposition of totally unimodular matrices and regular matroids
"""

from sage.structure.sage_object cimport SageObject


cdef class DecompositionNode(SageObject):

    def matrix(self):
        raise NotImplementedError

    def plot(self):
        raise NotImplementedError

    def _repr_(self):
        # print it
        raise NotImplementedError


cdef class SumNode(DecompositionNode):

    def permuted_block_matrix(self):
        r"Return (Prow, BlockMatrix, Pcolumn) so that self.matrix() == Prow * BlockMatrix * Pcolumn ????"
        raise NotImplementedError

    def summands(self):
        raise NotImplementedError


cdef class OneSumNode(SumNode):

    pass


cdef class BaseGraphicNode(DecompositionNode):

    pass


cdef class GraphicNode(BaseGraphicNode):

    def graph(self):
        raise NotImplementedError


cdef class CographicNode(BaseGraphicNode):

    pass


cdef class PlanarNode(BaseGraphicNode):

    pass


cdef class LeafNode(DecompositionNode):


    pass


cdef class K33Node(LeafNode):

    def graph(self):
        raise NotImplementedError

    pass
