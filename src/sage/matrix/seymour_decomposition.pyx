r"""
Seymour's decomposition of totally unimodular matrices and regular matroids
"""

from sage.libs.cmr.cmr cimport *
from sage.misc.cachefunc import cached_method
from sage.structure.sage_object cimport SageObject


cdef class DecompositionNode(SageObject):

    cdef _set_dec(self, CMR_DEC *dec, root):
        if self._root is None:
            CMR_CALL(CMRdecFree(cmr, &self._dec))
        self._dec = dec
        self._root = root

    def __dealloc__(self):
        self._set_dec(NULL, None)

    def matrix(self):
        raise NotImplementedError

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


cdef _class(CMR_DEC *dec):
    k = CMRdecIsSum(dec, NULL, NULL)
    if k == 1:
        return OneSumNode
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
    cdef DecompositionNode result = <DecompositionNode> _class(dec)()
    result._set_dec(dec, root)
    return result
