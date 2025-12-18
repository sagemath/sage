r"""
Vector partitions

AUTHORS:

- Amritanshu Prasad (2013): initial version
- Shriya M (2022): added new parameters such as ``distinct``, ``parts`` and
  ``is_repeatable``
"""
# ****************************************************************************
#       Copyright (C) 2013 Amritanshu Prasad <amri@imsc.res.in>
#                     2022 Shriya M <25shriya@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.combinat import CombinatorialElement
from sage.combinat.partition import Partition


def find_min(vect):
    """
    Return a string of ``0``'s with one ``1`` at the location where the list
    ``vect`` has its last entry which is not equal to ``0``.

    INPUT:

    - ``vec`` -- list of integers

    OUTPUT:

    A list of the same length with ``0``'s everywhere, except for a ``1``
    at the last position where ``vec`` has an entry not equal to ``0``.

    EXAMPLES::

        sage: from sage.combinat.vector_partition import find_min
        sage: find_min([2, 1])
        [0, 1]
        sage: find_min([2, 1, 0])
        [0, 1, 0]
    """
    i = len(vect)
    while vect[i-1] == 0 and i > 0:
        i = i-1
    min = [0]*len(vect)
    if i > 0:
        min[i-1] = 1
    return min


def IntegerVectorsIterator(vect, min=None):
    """
    Return an iterator over the list of integer vectors which are componentwise
    less than or equal to ``vect``, and lexicographically greater than or equal
    to ``min``.

    INPUT:

    - ``vect`` -- list of nonnegative integers
    - ``min`` -- list of nonnegative integers dominated elementwise by ``vect``

    OUTPUT:

    A list in lexicographic order of all integer vectors (as lists) which are
    dominated elementwise by ``vect`` and are greater than or equal to ``min`` in
    lexicographic order.

    EXAMPLES::

        sage: from sage.combinat.vector_partition import IntegerVectorsIterator
        sage: list(IntegerVectorsIterator([1, 1]))
        [[0, 0], [0, 1], [1, 0], [1, 1]]

        sage: list(IntegerVectorsIterator([1, 1], min = [1, 0]))
        [[1, 0], [1, 1]]
    """
    vect = list(vect)
    if not vect:
        yield []
    else:
        if min is None:
            min = [0] * len(vect)
        if vect < min:
            return
        else:
            for vec in IntegerVectorsIterator(vect[1:], min=min[1:]):
                yield [min[0]] + vec
            for j in range(min[0] + 1, vect[0] + 1):
                for vec in IntegerVectorsIterator(vect[1:]):
                    yield [j] + vec


class VectorPartition(CombinatorialElement):
    r"""
    A vector partition is a multiset of integer vectors.
    """
    @staticmethod
    def __classcall_private__(cls, vecpar):
        """
        Create a vector partition.

        EXAMPLES::

            sage: VectorPartition([[3, 2, 1], [2, 2, 1]])
            [[2, 2, 1], [3, 2, 1]]

        The parent class is the class of vector partitions of the sum of the
        vectors in ``vecpar``::

            sage: V = VectorPartition([[3, 2, 1], [2, 2, 1]])
            sage: V.parent()._vec
            (5, 4, 2)
        """
        vec = [sum([vec[i] for vec in vecpar]) for i in range(len(vecpar[0]))]
        P = VectorPartitions(vec)
        return P(vecpar)

    def __init__(self, parent, vecpar):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: elt =  VectorPartition([[3, 2, 1], [2, 2, 1]])
            sage: TestSuite(elt).run()
        """
        CombinatorialElement.__init__(self, parent, sorted(vecpar))

    def sum(self):
        """
        Return the sum vector as a list.

        EXAMPLES::

            sage: V = VectorPartition([[3, 2, 1], [2, 2, 1]])
            sage: V.sum()
            [5, 4, 2]
        """
        return list(self.parent()._vec)

    def partition_at_vertex(self, i):
        """
        Return the partition obtained by sorting the ``i``-th elements of
        the vectors in the vector partition.

        EXAMPLES::

            sage: V = VectorPartition([[1, 2, 1], [2, 4, 1]])
            sage: V.partition_at_vertex(1)
            [4, 2]
        """
        return Partition(sorted([vec[i] for vec in self._list], reverse=True))


class VectorPartitions(UniqueRepresentation, Parent):
    r"""
    Class of all vector partitions of ``vec`` with all parts greater than
    or equal to ``min`` in lexicographic order, with parts from ``parts``.

    A vector partition of ``vec`` is a list of vectors with nonnegative
    integer entries whose sum is ``vec``.

    INPUT:

    - ``vec`` -- integer vector
    - ``min`` -- integer vector dominated elementwise by ``vec``
    - ``parts`` -- finite list of possible parts
    - ``distinct`` -- boolean, set to ``True`` if only vector partitions with
      distinct parts are enumerated
    - ``is_repeatable`` -- boolean function on ``parts`` which gives ``True``
      in parts that can be repeated

    EXAMPLES:

    If ``min`` is not specified, then the class of all vector partitions of
    ``vec`` is created::

        sage: VP = VectorPartitions([2, 2])
        sage: for vecpar in VP:
        ....:     print(vecpar)
        [[0, 1], [0, 1], [1, 0], [1, 0]]
        [[0, 1], [0, 1], [2, 0]]
        [[0, 1], [1, 0], [1, 1]]
        [[0, 1], [2, 1]]
        [[0, 2], [1, 0], [1, 0]]
        [[0, 2], [2, 0]]
        [[1, 0], [1, 2]]
        [[1, 1], [1, 1]]
        [[2, 2]]

    If ``distinct`` is set to be True, then distinct part partitions are created::

        sage: VP = VectorPartitions([2,2], distinct = True)
        sage: list(VP)
        [[[0, 1], [1, 0], [1, 1]],
         [[0, 1], [2, 1]],
         [[0, 2], [2, 0]],
         [[1, 0], [1, 2]],
         [[2, 2]]]

    If ``min`` is specified, then the class consists of only those vector
    partitions whose parts are all greater than or equal to ``min`` in
    lexicographic order::

        sage: VP = VectorPartitions([2, 2], min = [1, 0])
        sage: for vecpar in VP:
        ....:     print(vecpar)
        [[1, 0], [1, 2]]
        [[1, 1], [1, 1]]
        [[2, 2]]
        sage: VP = VectorPartitions([2, 2], min = [1, 0], distinct = True)
        sage: for vecpar in VP:
        ....:     print(vecpar)
        [[1, 0], [1, 2]]
        [[2, 2]]

    If ``parts`` is specified, then the class consists only of those vector partitions
    whose parts are from ``parts``::

        sage: Vec_Par = VectorPartitions([2,2], parts=[[0,1],[1,0],[1,1]])
        sage: list(Vec_Par)
        [[[0, 1], [0, 1], [1, 0], [1, 0]], [[0, 1], [1, 0], [1, 1]], [[1, 1], [1, 1]]]

    If ``is_repeatable`` is specified, then the parts which satisfy the boolean function
    ``is_repeatable`` are allowed to be repeated::


        sage: Vector_Partitions = VectorPartitions([2,2], parts=[[0,1],[1,0],[1,1]], is_repeatable=lambda vec: sum(vec)%2!=0)
        sage: list(Vector_Partitions)
        [[[0, 1], [0, 1], [1, 0], [1, 0]], [[0, 1], [1, 0], [1, 1]]]
    """
    @staticmethod
    def __classcall_private__(cls, vec, min=None, parts=None, distinct=False, is_repeatable=None):
        r"""
        Create the class of vector partitions of ``vec`` where all parts
        are greater than or equal to the vector ``min``.

        EXAMPLES::

            sage: VP1 = VectorPartitions([2, 1])
            sage: VP2 = VectorPartitions((2, 1), min = [0, 1])
            sage: VP1 is VP2
            True
        """
        if min is None:
            min = find_min(vec)  # tuple([0 for v in vec[:-1]]+[1])
        if parts is None:
            parts = list(IntegerVectorsIterator(vec, min=min))
        if [0]*len(vec) in parts:
            parts.remove([0]*len(vec))
        if min in parts:
            min_index = parts.index(min)
            parts = parts[min_index:]
        parts = list(parts)
        for part_index in range(len(parts)):
            parts[part_index] = tuple(parts[part_index])
        return super().__classcall__(cls, tuple(vec), tuple(min), tuple(parts), distinct, is_repeatable)

    def __init__(self, vec, min=None, parts=None, distinct=False, is_repeatable=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: VP = VectorPartitions([2, 2])
            sage: TestSuite(VP).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._vec = vec
        self._min = min
        self._parts = parts
        self._distinct = distinct
        self._is_repeatable = is_repeatable

    def _element_constructor_(self, vecpar):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: VP = VectorPartitions([2, 2])
            sage: elt = VP([[1, 0], [1, 2]]); elt
            [[1, 0], [1, 2]]
            sage: elt.parent() is VP
            True
        """
        return self.element_class(self, vecpar)

    Element = VectorPartition

    def __iter__(self):
        r"""
        Iterator for vector partitions.

        EXAMPLES::

            sage: VP = VectorPartitions([2, 2])
            sage: VP.cardinality()
            9
        """
        if all(coord == 0 for coord in self._vec):
            yield self.element_class(self, []) # the zero vector has only the empty partition
        else:
            for part in self._parts: # choose the first part
                if tuple(part) == self._vec:
                    yield self.element_class(self, [list(part)])
                elif any(part[i] > self._vec[i] for i in range(len(self._vec))):
                    pass
                else:# recursively find all possibilities for the rest of the vector partition
                    new_vec = tuple(self._vec[i]-part[i] for i in range(len(self._vec)))
                    i = self._parts.index(part)
                    if self._is_repeatable is None:
                        if self._distinct:
                            new_parts = self._parts[i+1:]
                        else:
                            new_parts = self._parts[i:]
                    else:
                        if self._is_repeatable(part):
                            new_parts = self._parts[i:]
                        else:
                            new_parts = self._parts[i+1:]
                    for vecpar in VectorPartitions(new_vec, min=self._min, parts=new_parts, distinct=self._distinct, is_repeatable=self._is_repeatable):
                        yield self.element_class(self, [list(part)] + list(vecpar))
