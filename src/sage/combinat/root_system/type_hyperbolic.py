"""
Hyperbolic Coxeter types.
"""
# ****************************************************************************
#       Copyright (C) 2025 Samy Mekkati <samy.mekkati.1@ens.etsmtl.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.root_system.coxeter_type import CoxeterType
from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix
from sage.combinat.root_system.hyperbolic_coxeter_matrices import (
    hyperbolic_coxeter_matrices,
    mcmullen_notation)


class CoxeterType_Hyperbolic(CoxeterType):
    r"""
    Hyperbolic Coxeter type
    """
    def __init__(self, data):
        """
        EXAMPLES::

            sage: C = CoxeterType(["Hyperbolic", (142, 1, 3)])
            sage: C
            Coxeter type with Humphrey's datum
            (Page : 142, Column : 1, Row : 3)
            sage: C = CoxeterType(["Dh", 8])
            sage: C
            Coxeter type of ['Dh', 8] with Humphrey's datum
            (Page : 144, Column : 1, Row : 2)

        TESTS::

            sage: TestSuite(C).run()
        """
        if data[0] == "Hyperbolic":
            self._position = tuple(data[1])

            if self._position not in hyperbolic_coxeter_matrices:
                raise ValueError(f"position {self._position} is not a valid")

            if self._position in mcmullen_notation:
                self._prefix = mcmullen_notation[self._position][0]
                self._index = mcmullen_notation[self._position][1]
        else:
            self._prefix = data[0]
            self._index = data[1]

            if (self._prefix, self._index) in mcmullen_notation:
                self._position = mcmullen_notation[(self._prefix, self._index)]

        super().__init__()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        This method returns a string that describes the Coxeter type,
        including its Humphrey's datum and McMullen's notation if defined.

        EXAMPLES::

            sage: C = CoxeterType(["Hyperbolic", (142, 1, 1)])
            sage: C
            Coxeter type with Humphrey's datum
            (Page : 142, Column : 1, Row : 1)
            sage: C2 = CoxeterType(["Ah", 6])
            sage: C2
            Coxeter type of ['Ah', 6] with Humphrey's datum
            (Page : 143, Column : 2, Row : 3)
        """
        a, b, c = self._position

        if hasattr(self, "_prefix") and hasattr(self, "_index"):
            return (
                f"Coxeter type of ['{self._prefix}', {self._index}] "
                f"with Humphrey's datum (Page : {a}, Column : {b}, Row : {c})"
            )
        else:
            return (
                f"Coxeter type with Humphrey's datum "
                f"(Page : {a}, Column : {b}, Row : {c})"
            )

    def rank(self):
        """
        Return the rank of ``self``.

        This is the number of nodes of the associated Coxeter graph.

        EXAMPLES::

            sage: CoxeterType(['Hyperbolic', (144,1,10)]).rank()
            10
            sage: CoxeterType(['Hyperbolic', (142,2,2)]).rank()
            4
            sage: CoxeterType(['L', 633]).rank()
            4
        """
        return len(hyperbolic_coxeter_matrices[self._position])

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: CoxeterType(["Hyperbolic", (142, 2, 1)]).coxeter_matrix()
            [1 4 2 2]
            [4 1 4 2]
            [2 4 1 3]
            [2 2 3 1]
            sage: CoxeterType(["Dh", 9]).coxeter_matrix()
            [1 2 3 2 2 2 2 2 2]
            [2 1 3 2 2 2 2 2 2]
            [3 3 1 3 2 2 2 2 2]
            [2 2 3 1 3 2 2 2 2]
            [2 2 2 3 1 2 3 2 2]
            [2 2 2 2 2 1 3 2 2]
            [2 2 2 2 3 3 1 3 2]
            [2 2 2 2 2 2 3 1 3]
            [2 2 2 2 2 2 2 3 1]
        """
        return CoxeterMatrix(hyperbolic_coxeter_matrices[self._position])

    def humphreys_reference(self):
        """
        Return a string with the reference to Humphreys'
        Reflection groups and Coxeter groups.

        The reference is given by the page, column and row of
        the table in the book.

        EXAMPLES::

            sage: C = CoxeterType(["Hyperbolic", (142, 1, 3)])
            sage: C.humphreys_reference()
            'Page : 142, Column : 1, Row : 3'
            sage: CoxeterType(["Dh", 8]).humphreys_reference()
            'Page : 144, Column : 1, Row : 2'
        """
        return 'Page : {}, Column : {}, Row : {}'.format(*self._position)

    def coxeter_graph(self):
        """
        Return the Coxeter graph associated to ``self``.

        EXAMPLES::

            sage: C = CoxeterType(["Hyperbolic",(141, 2, 1)])
            sage: C.coxeter_graph()
            Graph on 4 vertices
            sage: C2 = CoxeterType(["Hyperbolic",(144, 1, 1)])
            sage: C2.coxeter_graph()
            Graph on 8 vertices
            sage: C3 = CoxeterType(["Dh", 9])
            sage: C3.coxeter_graph()
            Graph on 9 vertices
        """
        return self.coxeter_matrix().coxeter_graph()

    def is_hyperbolic(self):
        """
        Return ``True`` since ``self`` is a hyperbolic Coxeter type.

        EXAMPLES::
            sage: CoxeterType(["Hyperbolic", (142, 1, 3)]).is_hyperbolic()
            True
            sage: CoxeterType(["Bh", 5]).is_hyperbolic()
            True
        """
        return True

    def index_set(self):
        """
        Return the index set for ``self``.

        This is the list of the nodes of the associated Coxeter graph.

        EXAMPLES::

            sage: CoxeterType(["Hyperbolic",(144, 1, 1)]).index_set()
            (1, 2, 3, 4, 5, 6, 7, 8)
            sage: CoxeterType(["Hyperbolic",(143, 1, 6)]).index_set()
            (1, 2, 3, 4, 5, 6)
            sage: CoxeterType(["Q", 4]).index_set()
            (1, 2, 3, 4)
        """
        return self.coxeter_matrix().index_set()

    def is_affine(self):
        """
        Return ``False`` because hyperbolic Coxeter graphs are never affine.

        EXAMPLES::

            sage: CoxeterType(["Hyperbolic",(142, 3, 4)]).is_affine()
            False
            sage: CoxeterType(["Eh", 6]).is_affine()
            False
        """
        return False

    def is_finite(self):
        """
        Return ``False`` because hyperbolic Coxeter graphs are never finite.

        EXAMPLES::

            sage: CoxeterType(["Hyperbolic",(142, 3, 4)]).is_finite()
            False
            sage: CoxeterType(["Dh", 8]).is_finite()
            False
        """
        return False

    def is_crystallographic(self):
        """
        Return whether ``self`` is crystallographic.

        EXAMPLES::

            sage: CoxeterType(["Hyperbolic",(142, 3, 4)]).is_crystallographic()
            True
            sage: CoxeterType(["Hyperbolic",(141, 1, 1)]).is_crystallographic()
            False
        """
        return self.coxeter_matrix().is_crystallographic()

    def __eq__(self, other):
        """
        Return whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: C1 = CoxeterType(["Hyperbolic", (141, 1, 4)])
            sage: C2 = CoxeterType(["K", 53])
            sage: C1 == C2
            True
            sage: C3 = CoxeterType(["Hyperbolic", (141, 1, 1)])
            sage: C1 == C3
            False
        """
        if isinstance(other, CoxeterType_Hyperbolic):
            if self.coxeter_matrix() == other.coxeter_matrix():
                return True
        return False