"""
Overrides to unpickle old matrix groups
"""

from sage.structure.sage_object import register_unpickle_override

from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
from sage.groups.matrix_gps.group_element_gap import MatrixGroupElement_gap
from sage.groups.matrix_gps.linear import GL, LinearMatrixGroup_generic


class LegacyMatrixGroup(FinitelyGeneratedMatrixGroup_gap):

    def __setstate__(self, state):
        """
        Restore from old pickle.

        EXAMPLES::

            sage: from sage.groups.matrix_gps.pickling_overrides import LegacyMatrixGroup
            sage: state = dict()
            sage: state['_MatrixGroup_gap__n'] = 2
            sage: state['_MatrixGroup_gap__R'] = GF(3)
            sage: state['_gensG'] = [ matrix(GF(3), [[1,2],[0,1]]) ]
            sage: M = LegacyMatrixGroup.__new__(LegacyMatrixGroup)
            sage: M.__setstate__(state)
            sage: M
            Matrix group over Finite Field of size 3 with 1 generators (
            [1 2]
            [0 1]
            )
        """
        matrix_gens = state['_gensG']
        ring = state['_MatrixGroup_gap__R']
        degree = state['_MatrixGroup_gap__n']
        from sage.libs.gap.libgap import libgap
        libgap_group = libgap.Group(libgap(matrix_gens))
        self.__init__(degree, ring, libgap_group)


register_unpickle_override(
    'sage.groups.matrix_gps.matrix_group', 'MatrixGroup_gens_finite_field',
    LegacyMatrixGroup)


class LegacyMatrixGroupElement(MatrixGroupElement_gap):

    def __setstate__(self, state):
        """
        Restore from old pickle.

        EXAMPLES::

            sage: from sage.groups.matrix_gps.pickling_overrides import LegacyMatrixGroup, LegacyMatrixGroupElement
            sage: state = dict()
            sage: state['_MatrixGroup_gap__n'] = 2
            sage: state['_MatrixGroup_gap__R'] = GF(3)
            sage: state['_gensG'] = [ matrix(GF(3), [[1,2],[0,1]]) ]
            sage: M = LegacyMatrixGroup.__new__(LegacyMatrixGroup)
            sage: M.__setstate__(state)
            sage: M
            Matrix group over Finite Field of size 3 with 1 generators (
            [1 2]
            [0 1]
            )
            sage: state = [ M, {'_MatrixGroupElement__mat':matrix(GF(3), [[1,2],[0,1]])} ]
            sage: m = LegacyMatrixGroupElement.__new__(LegacyMatrixGroupElement)
            sage: m.__setstate__(state)
            sage: m
            [1 2]
            [0 1]
        """
        parent = state[0]
        m = state[1]['_MatrixGroupElement__mat']
        m = parent.matrix_space()(m)
        self.__init__(parent, m, check=False)


register_unpickle_override(
    'sage.groups.matrix_gps.matrix_group_element', 'MatrixGroupElement',
    LegacyMatrixGroupElement)


class LegacyGeneralLinearGroup(LinearMatrixGroup_generic):

    def __setstate__(self, state):
        """
        Restore from old pickle.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: from sage.groups.matrix_gps.pickling_overrides import LegacyGeneralLinearGroup
            sage: state = dict()
            sage: state['_MatrixGroup_gap__n'] = 2
            sage: state['_MatrixGroup_gap__R'] = ZZ
            sage: M = Group.__new__(LegacyGeneralLinearGroup)
            sage: M.__setstate__(state)
            sage: M
            General Linear Group of degree 2 over Integer Ring
        """
        ring = state['_MatrixGroup_gap__R']
        n = state['_MatrixGroup_gap__n']
        G = GL(n, ring)
        self.__init__(G.degree(), G.base_ring(), G._special, G._name_string, G._latex_string)


register_unpickle_override(
    'sage.groups.matrix_gps.general_linear', 'GeneralLinearGroup_finite_field',
    LegacyGeneralLinearGroup)
