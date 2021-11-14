r"""
Quiver mutation types

AUTHORS:

- Gregg Musiker (2012, initial version)
- Christian Stump (2012, initial version)
- Hugh Thomas (2012, initial version)
"""
# ****************************************************************************
#       Copyright (C) 2011 Gregg Musiker <gmusiker@gmail.com>
#                          Christian Stump <christian.stump@gmail.com>
#                          Hugh Thomas <hugh@math.unb.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************
from __future__ import annotations
import os
import pickle

from sage.structure.sage_object import SageObject
from copy import copy
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.graphs.all import Graph, DiGraph
from sage.arith.all import binomial, euler_phi
from sage.misc.misc_c import prod
from sage.matrix.constructor import matrix


class QuiverMutationTypeFactory(SageObject):
    def __call__(self, *args):
        """
        For a detailed description, see :meth:`QuiverMutationType`.

        EXAMPLES::

            sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import QuiverMutationTypeFactory
            sage: QuiverMutationTypeFactory()
            QuiverMutationType
        """
        # get data as arguments or as list/tuple
        if len(args) == 1:
            data = args[0]
        else:
            data = args

        # data is a QuiverMutationType
        if isinstance(data, QuiverMutationType_Irreducible):
            return data
        elif isinstance(data, QuiverMutationType_Reducible):
            return data

        # check that data is a tuple or list
        if isinstance(data, tuple) and data:
            pass
        elif isinstance(data, list) and data:
            data = tuple(data)
        else:
            _mutation_type_error(data)

        # check for reducible types
        if all(type(data_component) in [list, tuple, QuiverMutationType_Irreducible] for data_component in data):
            if len(data) == 1:
                return QuiverMutationType(data[0])
            else:
                data = tuple(QuiverMutationType(comp) for comp in data)
                return QuiverMutationType_Reducible(*data)

        # check for irreducible types
        if len(data) == 2:
            data = (data[0], data[1], None)
        elif len(data) == 3:
            pass
        else:
            _mutation_type_error(data)

        if isinstance(data[2], list):
            data = (data[0], data[1], tuple(data[2]))
        if isinstance(data[1], list):
            data = (data[0], tuple(data[1]), data[2])

        # mutation type casting
        if data == ('D', 2, None):
            return QuiverMutationType(('A', 1, None), ('A', 1, None))
        elif data == ('D', 3, None):
            data = ('A', 3, None)
        elif data == ('C', 2, None):
            data = ('B', 2, None)
        elif data == ('E', 9, None):
            data = ('E', 8, 1)
        elif data[0] == 'A' and data[2] == 1 and isinstance(data[1], tuple) and len(data[1]) == 2 and min(data[1]) == 0:
            if max(data[1]) == 0:
                pass
            elif max(data[1]) == 1:
                data = ('A', 1, None)
            elif max(data[1]) == 2:
                return QuiverMutationType(('A', 1, None), ('A', 1, None))
            elif max(data[1]) == 3:
                data = ('A', 3, None)
            else:
                data = ('D', max(data[1]), None)
        elif data[0] == 'GR' and data[2] is None and isinstance(data[1], tuple) and len(data[1]) == 2 and data[1][1] > data[1][0]:
            if min(data[1]) > max(data[1]) / 2 and max(data[1]) != min(data[1]) + 1:
                data = (data[0], (max(data[1]) - min(data[1]), max(data[1])), data[2])
            if min(data[1]) == 2 and max(data[1]) > 3:
                data = ('A', max(data[1]) - 3, None)
            elif data[1] == (3, 6):
                data = ('D', 4, None)
            elif data[1] == (3, 7):
                data = ('E', 6, None)
            elif data[1] == (3, 8):
                data = ('E', 8, None)
            elif data[1] == (3, 9):
                data = ('E', 8, [1, 1])
            elif data[1] == (4, 8):
                data = ('E', 7, [1, 1])
        elif data == ('TR', 1, None):
            data = ('A', 1, None)
        elif data == ('TR', 2, None):
            data = ('A', 3, None)
        elif data == ('TR', 3, None):
            data = ('D', 6, None)
        elif data == ('TR', 4, None):
            data = ('E', 8, (1, 1))
        # mutation type casting from Kac conventions
        elif data == ('A', 1, 1):
            data = ('A', (1, 1), 1)
        elif data[0] == 'B' and data[2] == 1:
            if data[1] == 2:
                data = ('CC', 2, 1)
            elif data[1] > 2:
                data = ('BD', data[1], 1)
        elif data[0] == 'B' and data[2] == -1:
            if data[1] == 2:
                data = ('BB', 2, 1)
            elif data[1] > 2:
                data = ('CD', data[1], 1)
        elif data[0] == 'C' and data[1] > 1 and data[2] == 1:
            data = ('CC', data[1], 1)
        elif data[0] == 'C' and data[1] > 1 and data[2] == -1:
            data = ('BB', data[1], 1)
        elif data == ('A', 2, 2):
            data = ('BC', 1, 1)
        elif data[0] == 'A' and data[1] in ZZ and data[1] > 1 and data[1] % 2 == 0 and data[2] == 2:
            data = ('BC', data[1] // 2, 1)
        elif data[0] == 'A' and data[1] in ZZ and data[1] > 3 and data[1] % 2 and data[2] == 2:
            data = ('CD', (data[1] + 1) // 2, 1)
        # We think of ('A',3, 2) as ('D',3, 2)
        elif data == ('A', 3, 2):
            data = ('BB', 2, 1)
        elif data[0] == 'D' and data[1] in ZZ and data[1] > 2 and data[2] == 2:
            data = ('BB', data[1] - 1, 1)
        elif data == ('E', 6, 2):
            data = ('F', 4, -1)
        elif data == ('D', 4, 3):
            data = ('G', 2, -1)
        elif data == ('F', 4, (2, 1)):
            data = ('F', 4, (1, 2))
        elif data == ('G', 2, (3, 1)):
            data = ('G', 2, (1, 3))
        elif data[0] == 'T' and data[2] is None:
            data = (data[0], tuple(sorted(data[1])), data[2])
            r, p, q = data[1]
            if r == 1:
                data = ('A', p + q - 1, None)
            elif r == p == 2:
                data = ('D', q + 2, None)
            elif r == 2 and p == 3:
                if q in (3, 4, 5):
                    data = ('E', q + 3, None)
                elif q == 6:
                    data = ('E', 8, 1)
                else:
                    data = ('E', q + 3, None)
            elif r == 2 and p == q == 4:
                data = ('E', 7, 1)
            elif r == p == q == 3:
                data = ('E', 6, 1)
        elif data[0] == 'R2' and data[2] is None and all(data[1][i] in ZZ and data[1][i] > 0 for i in [0, 1]):
            data = (data[0], tuple(sorted(data[1])), data[2])
            if data[1] == (1, 1):
                data = ('A', 2, None)
            elif data[1] == (1, 2):
                data = ('B', 2, None)
            elif data[1] == (1, 3):
                data = ('G', 2, None)
            elif data[1] == (1, 4):
                data = ('BC', 1, 1)
            elif data[1] == (2, 2):
                data = ('A', (1, 1), 1)

        # setting the parameters and returning the mutation type
        letter, rank, twist = data
        if not isinstance(letter, str):
            _mutation_type_error(data)
        if isinstance(rank, list):
            rank = tuple(rank)
        if isinstance(twist, list):
            twist = tuple(twist)
        return QuiverMutationType_Irreducible(letter, rank, twist)

    def _repr_(self) -> str:
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: QuiverMutationType  # indirect doctest
            QuiverMutationType
        """
        return "QuiverMutationType"

    def samples(self, finite=None, affine=None, elliptic=None,
                mutation_finite=None):
        """
        Return a sample of the available quiver mutations types.

        INPUT:

        - ``finite``

        - ``affine``

        - ``elliptic``

        - ``mutation_finite``

        All four input keywords default values are ``None``. If
        set to ``True`` or ``False``, only these samples are returned.

        EXAMPLES::

            sage: QuiverMutationType.samples()
            [['A', 1], ['A', 5], ['B', 2], ['B', 5], ['C', 3],
             ['C', 5], [ ['A', 1], ['A', 1] ], ['D', 5], ['E', 6],
             ['E', 7], ['E', 8], ['F', 4], ['G', 2],
             ['A', [1, 1], 1], ['A', [4, 5], 1], ['D', 4, 1],
             ['BB', 5, 1], ['E', 6, [1, 1]], ['E', 7, [1, 1]],
             ['R2', [1, 5]], ['R2', [3, 5]], ['E', 10], ['BE', 5],
             ['GR', [3, 10]], ['T', [3, 3, 4]]]

            sage: QuiverMutationType.samples(finite=True)
            [['A', 1], ['A', 5], ['B', 2], ['B', 5], ['C', 3],
             ['C', 5], [ ['A', 1], ['A', 1] ], ['D', 5], ['E', 6],
             ['E', 7], ['E', 8], ['F', 4], ['G', 2]]

            sage: QuiverMutationType.samples(affine=True)
            [['A', [1, 1], 1], ['A', [4, 5], 1], ['D', 4, 1], ['BB', 5, 1]]

            sage: QuiverMutationType.samples(elliptic=True)
            [['E', 6, [1, 1]], ['E', 7, [1, 1]]]

            sage: QuiverMutationType.samples(mutation_finite=False)
            [['R2', [1, 5]], ['R2', [3, 5]], ['E', 10], ['BE', 5],
             ['GR', [3, 10]], ['T', [3, 3, 4]]]
        """
        result = self._samples()
        if finite is not None:
            result = [t for t in result if t.is_finite() == finite]
        if affine is not None:
            result = [t for t in result if t.is_affine() == affine]
        if elliptic is not None:
            result = [t for t in result if t.is_elliptic() == elliptic]
        if mutation_finite is not None:
            result = [t for t in result
                      if t.is_mutation_finite() == mutation_finite]
        return result

    @cached_method
    def _samples(self):
        """
        Return a list of sample of available Cartan types.

        EXAMPLES::

            sage: X = QuiverMutationType._samples()
        """
        finite_types = \
            [QuiverMutationType(t) for t in [['A', 1], ['A', 5], ['B', 2], ['B', 5],
                                             ['C', 3], ['C', 5], ['D', 2], ['D', 5],
                                             ["E", 6], ["E", 7], ["E", 8], ["F", 4],
                                             ["G", 2]]]
        affine_types = \
            [QuiverMutationType(t) for t in [['A', [1, 1], 1], ['A', [4, 5], 1], ['D', 4, 1], ['BB', 5, 1]]]
        elliptic_types = \
            [QuiverMutationType(t) for t in [['E', 6, [1, 1]], ['E', 7, [1, 1]]]]
        mutation_finite_types = \
            [QuiverMutationType(t) for t in [['R2', (1, 5)], ['R2', (3, 5)]]]
        mutation_infinite_types = \
            [QuiverMutationType(t) for t in [['E', 10], ['BE', 5], ['GR', (3, 10)], ['T', (3, 3, 4)]]]

        return finite_types + affine_types + elliptic_types + mutation_finite_types + mutation_infinite_types


QuiverMutationType = QuiverMutationTypeFactory()


QuiverMutationType.__doc__ = \
    r"""

*Quiver mutation types* can be seen as a slight generalization of
 *generalized Cartan types*.

Background on generalized Cartan types can be found at

        :wikipedia:`Generalized_Cartan_matrix`

For the compendium on the cluster algebra and quiver package in Sage see [MS2011]_

A `B`-matrix is a skew-symmetrizable `( n \times n )`-matrix `M`.
I.e., there exists an invertible diagonal matrix `D` such that `DM` is
skew-symmetric.  `M` can be encoded as a *quiver* by having a directed
edge from vertex `i` to vertex `j` with label `(a,b)` if `a = M_{i,j}
> 0` and `b = M_{j,i} < 0`.  We consider quivers up to *mutation
equivalence*.

To a quiver mutation type we can associate a *generalized Cartan type*
by sending `M` to the generalized Cartan matrix `C(M)` obtained by
replacing all positive entries by their negatives and adding `2`'s on
the main diagonal.

``QuiverMutationType`` constructs a quiver mutation type object. For
more detail on the possible different types, please see the
compendium.

INPUT:

The input consists either of a quiver mutation type, or of a
``letter`` (a string), a ``rank`` (one integer or a list/tuple of
integers), and an optional ``twist`` (an integer or a list of
integers).  There are several different naming conventions for quiver
mutation types.

- Finite type -- ``letter`` is a Dynkin type (A-G), and ``rank`` is
  the rank.

- Affine type -- there is more than one convention for naming affine
  types.

    * Kac's notation: ``letter`` is a Dynkin type, ``rank`` is the
      rank of the associated finite Dynkin diagram, and ``twist`` is the
      twist, which could be 1, 2, or 3.  In the special case of affine
      type A, there is more than one quiver mutation type associated to
      the Cartan type.  In this case only, ``rank`` is a pair of integers
      (i,j), giving the number of edges pointing clockwise and the number
      of edges pointing counter-clockwise.  The total number of vertices
      is given by i+j in this case.

    * Naive notation: ``letter`` is one of 'BB', 'BC', 'BD', 'CC',
      'CD'.  The name specifies the two ends of the diagram, which are
      joined by a path.  The total number of vertices is given by
      ``rank +1`` (to match the indexing people expect because these
      are affine types).  In general, ``rank`` must be large enough
      for the picture to make sense, but we accept ``letter`` is
      ``BC`` and ``rank=1``.

    * Macdonald notation: for the dual of an untwisted affine type
      (such as ['C', 6, 1]), we accept a twist of -1 (i.e.,
      ['C',6,-1]).

- Elliptic type -- ``letter`` is a Dynkin type, ``rank`` is the rank
  of the finite Dynkin diagram, and ``twist`` is a tuple of two
  integers.  We follow Saito's notation.

- Other shapes:

    * Rank 2: ``letter`` is 'R2', and ``rank`` is a pair of integers
      specifying the label on the unique edge.

    * Triangle: ``letter`` is ``TR``, and ``rank`` is the number of
      vertices along a side.

    * T: This defines a quiver shaped like a T.  ``letter`` is 'T',
      and the ``rank`` is a triple, whose entries specify the number
      of vertices along each path from the branch point (counting the
      branch point).

    * Grassmannian: This defines the cluster algebra (without
      coefficients) corresponding to the cluster algebra with
      coefficients which is the coordinate ring of a Grassmannian.
      ``letter`` is 'GR'.  ``rank`` is a pair of integers (`k`, `n`)
      with 'k' < 'n' specifying the Grassmannian of `k`-planes in
      `n`-space.  This defines a quiver given by a (k-1) x (n-k-1)
      grid where each square is cyclically oriented.

    * Exceptional mutation finite quivers: The two exceptional
      mutation finite quivers, found by Derksen-Owen, have ``letter``
      as 'X' and ``rank`` 6 or 7, equal to the number of vertices.

    * AE, BE, CE, DE: Quivers are built of one end which looks like
      type (affine A), B, C, or D, and the other end which looks like
      type E (i.e., it consists of two antennae, one of length one,
      and one of length two).  ``letter`` is 'AE', 'BE', 'CE', or
      'DE', and ``rank`` is the total number of vertices.  Note that
      'AE' is of a slightly different form and requires ``rank`` to be
      a pair of integers (i,j) just as in the case of affine type A.
      See Exercise 4.3 in Kac's book Infinite Dimensional Lie Algebras
      for more details.

    * Infinite type E: It is also possible to obtain infinite-type E
      quivers by specifying ``letter`` as 'E' and ``rank`` as the
      number of vertices.

REFERENCES:

- A good reference for finite and affine Dynkin diagrams, including
  Kac's notation, is the :wikipedia:`Dynkin_diagram`.

- A good reference for the skew-symmetrizable elliptic diagrams is
  "Cluster algebras of finite mutation type via unfolding" by
  A. Felikson, M. Shapiro, and P. Tumarkin, [FST2012]_.

EXAMPLES:

Finite types::

    sage: QuiverMutationType('A', 1)
    ['A', 1]
    sage: QuiverMutationType('A',5)
    ['A', 5]

    sage: QuiverMutationType('B', 2)
    ['B', 2]
    sage: QuiverMutationType('B',5)
    ['B', 5]

    sage: QuiverMutationType('C', 2)
    ['B', 2]
    sage: QuiverMutationType('C',5)
    ['C', 5]

    sage: QuiverMutationType('D', 2)
    [ ['A', 1], ['A', 1] ]
    sage: QuiverMutationType('D',3)
    ['A', 3]
    sage: QuiverMutationType('D',4)
    ['D', 4]

    sage: QuiverMutationType('E',6)
    ['E', 6]

    sage: QuiverMutationType('G', 2)
    ['G', 2]

    sage: QuiverMutationType('A',(1,0), 1)
    ['A', 1]

    sage: QuiverMutationType('A',(2,0), 1)
    [ ['A', 1], ['A', 1] ]

    sage: QuiverMutationType('A',(7,0), 1)
    ['D', 7]

Affine types::

    sage: QuiverMutationType('A',(1, 1), 1)
    ['A', [1, 1], 1]
    sage: QuiverMutationType('A',(2,4), 1)
    ['A', [2, 4], 1]

    sage: QuiverMutationType('BB', 2, 1)
    ['BB', 2, 1]
    sage: QuiverMutationType('BB',4, 1)
    ['BB', 4, 1]

    sage: QuiverMutationType('CC', 2, 1)
    ['CC', 2, 1]
    sage: QuiverMutationType('CC',4, 1)
    ['CC', 4, 1]

    sage: QuiverMutationType('BC', 1, 1)
    ['BC', 1, 1]
    sage: QuiverMutationType('BC',5, 1)
    ['BC', 5, 1]

    sage: QuiverMutationType('BD',3, 1)
    ['BD', 3, 1]
    sage: QuiverMutationType('BD',5, 1)
    ['BD', 5, 1]

    sage: QuiverMutationType('CD',3, 1)
    ['CD', 3, 1]
    sage: QuiverMutationType('CD',5, 1)
    ['CD', 5, 1]

    sage: QuiverMutationType('D',4, 1)
    ['D', 4, 1]
    sage: QuiverMutationType('D',6, 1)
    ['D', 6, 1]

    sage: QuiverMutationType('E',6, 1)
    ['E', 6, 1]
    sage: QuiverMutationType('E',7, 1)
    ['E', 7, 1]
    sage: QuiverMutationType('E',8, 1)
    ['E', 8, 1]

    sage: QuiverMutationType('F',4, 1)
    ['F', 4, 1]
    sage: QuiverMutationType('F',4,-1)
    ['F', 4, -1]

    sage: QuiverMutationType('G', 2, 1)
    ['G', 2, 1]
    sage: QuiverMutationType('G', 2,-1)
    ['G', 2, -1]
    sage: QuiverMutationType('A',3, 2) == QuiverMutationType('D',3, 2)
    True

Affine types using Kac's Notation::

    sage: QuiverMutationType('A', 1, 1)
    ['A', [1, 1], 1]
    sage: QuiverMutationType('B',5, 1)
    ['BD', 5, 1]
    sage: QuiverMutationType('C',5, 1)
    ['CC', 5, 1]
    sage: QuiverMutationType('A', 2, 2)
    ['BC', 1, 1]
    sage: QuiverMutationType('A',7, 2)
    ['CD', 4, 1]
    sage: QuiverMutationType('A',8, 2)
    ['BC', 4, 1]
    sage: QuiverMutationType('D',6, 2)
    ['BB', 5, 1]
    sage: QuiverMutationType('E',6, 2)
    ['F', 4, -1]
    sage: QuiverMutationType('D',4,3)
    ['G', 2, -1]

Elliptic types::

    sage: QuiverMutationType('E',6,[1, 1])
    ['E', 6, [1, 1]]
    sage: QuiverMutationType('F',4,[2, 1])
    ['F', 4, [1, 2]]
    sage: QuiverMutationType('G', 2,[3,3])
    ['G', 2, [3, 3]]

Mutation finite types:

Rank 2 cases::

    sage: QuiverMutationType('R2',(1, 1))
    ['A', 2]
    sage: QuiverMutationType('R2',(1, 2))
    ['B', 2]
    sage: QuiverMutationType('R2',(1,3))
    ['G', 2]
    sage: QuiverMutationType('R2',(1,4))
    ['BC', 1, 1]
    sage: QuiverMutationType('R2',(1,5))
    ['R2', [1, 5]]
    sage: QuiverMutationType('R2',(2, 2))
    ['A', [1, 1], 1]
    sage: QuiverMutationType('R2',(3,5))
    ['R2', [3, 5]]

Exceptional Derksen-Owen quivers::

    sage: QuiverMutationType('X',6)
    ['X', 6]


(Mainly) mutation infinite types:

Infinite type E::

    sage: QuiverMutationType('E',9)
    ['E', 8, 1]
    sage: QuiverMutationType('E', 10)
    ['E', 10]
    sage: QuiverMutationType('E', 12)
    ['E', 12]

    sage: QuiverMutationType('AE',(2,3))
    ['AE', [2, 3]]
    sage: QuiverMutationType('BE',5)
    ['BE', 5]
    sage: QuiverMutationType('CE',5)
    ['CE', 5]
    sage: QuiverMutationType('DE',6)
    ['DE', 6]

Grassmannian types::

    sage: QuiverMutationType('GR',(2,4))
    ['A', 1]
    sage: QuiverMutationType('GR',(2,6))
    ['A', 3]
    sage: QuiverMutationType('GR',(3,6))
    ['D', 4]
    sage: QuiverMutationType('GR',(3,7))
    ['E', 6]
    sage: QuiverMutationType('GR',(3,8))
    ['E', 8]
    sage: QuiverMutationType('GR',(3, 10))
    ['GR', [3, 10]]

Triangular types::

    sage: QuiverMutationType('TR', 2)
    ['A', 3]
    sage: QuiverMutationType('TR',3)
    ['D', 6]
    sage: QuiverMutationType('TR',4)
    ['E', 8, [1, 1]]
    sage: QuiverMutationType('TR',5)
    ['TR', 5]

T types::

    sage: QuiverMutationType('T',(1, 1, 1))
    ['A', 1]
    sage: QuiverMutationType('T',(1, 1,4))
    ['A', 4]
    sage: QuiverMutationType('T',(1,4,4))
    ['A', 7]
    sage: QuiverMutationType('T',(2, 2, 2))
    ['D', 4]
    sage: QuiverMutationType('T',(2, 2,4))
    ['D', 6]
    sage: QuiverMutationType('T',(2,3,3))
    ['E', 6]
    sage: QuiverMutationType('T',(2,3,4))
    ['E', 7]
    sage: QuiverMutationType('T',(2,3,5))
    ['E', 8]
    sage: QuiverMutationType('T',(2,3,6))
    ['E', 8, 1]
    sage: QuiverMutationType('T',(2,3,7))
    ['E', 10]
    sage: QuiverMutationType('T',(3,3,3))
    ['E', 6, 1]
    sage: QuiverMutationType('T',(3,3,4))
    ['T', [3, 3, 4]]

Reducible types::

    sage: QuiverMutationType(['A',3],['B',4])
    [ ['A', 3], ['B', 4] ]
"""


class QuiverMutationType_abstract(UniqueRepresentation, SageObject):
    """
    EXAMPLES::

        sage: mut_type1 = QuiverMutationType('A',5)
        sage: mut_type2 = QuiverMutationType('A',5)
        sage: mut_type3 = QuiverMutationType('A',6)
        sage: mut_type1 == mut_type2
        True
        sage: mut_type1 == mut_type3
        False
    """

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: QuiverMutationType(['A', 2]) # indirect doctest
            ['A', 2]
        """
        return self._description

    def plot(self, circular=False, directed=True):
        """
        Return the plot of the underlying graph or digraph of ``self``.

        INPUT:

        - ``circular`` -- (default:``False``) if ``True``, the
          circular plot is chosen, otherwise >>spring<< is used.

        - ``directed`` -- (default: ``True``) if ``True``, the
          directed version is shown, otherwise the undirected.

        EXAMPLES::

            sage: QMT = QuiverMutationType(['A',5])
            sage: pl = QMT.plot()
            sage: pl = QMT.plot(circular=True)
        """
        return self.standard_quiver().plot(circular=circular, directed=directed)

    def show(self, circular=False, directed=True):
        """
        Show the plot of the underlying digraph of ``self``.

        INPUT:

        - ``circular`` -- (default:``False``) if ``True``, the
          circular plot is chosen, otherwise >>spring<< is used.

        - ``directed`` -- (default: ``True``) if ``True``, the
          directed version is shown, otherwise the undirected.

        TESTS::

            sage: QMT = QuiverMutationType(['A',5])
            sage: QMT.show() # long time
        """
        self.plot(circular=circular, directed=directed).show()

    def letter(self):
        """
        Return the classification letter of ``self``.

        EXAMPLES::

            sage: mut_type = QuiverMutationType( ['A',5] ); mut_type
            ['A', 5]
            sage: mut_type.letter()
            'A'

            sage: mut_type = QuiverMutationType( ['BC',5, 1] ); mut_type
            ['BC', 5, 1]
            sage: mut_type.letter()
            'BC'

            sage: mut_type = QuiverMutationType(['A',3],['B',3]); mut_type
            [ ['A', 3], ['B', 3] ]
            sage: mut_type.letter()
            'A x B'

            sage: mut_type = QuiverMutationType(['A',3],['B',3],['X',6]); mut_type
            [ ['A', 3], ['B', 3], ['X', 6] ]
            sage: mut_type.letter()
            'A x B x X'
        """
        return self._letter

    def rank(self):
        """
        Return the rank in the standard quiver of ``self``.

        The rank is the number of vertices.

        EXAMPLES::

            sage: mut_type = QuiverMutationType( ['A',5] ); mut_type
            ['A', 5]
            sage: mut_type.rank()
            5

            sage: mut_type = QuiverMutationType( ['A',[4,5], 1] ); mut_type
            ['A', [4, 5], 1]
            sage: mut_type.rank()
            9

            sage: mut_type = QuiverMutationType( ['BC',5, 1] ); mut_type
            ['BC', 5, 1]
            sage: mut_type.rank()
            6

            sage: mut_type = QuiverMutationType(['A',3],['B',3]); mut_type
            [ ['A', 3], ['B', 3] ]
            sage: mut_type.rank()
            6

            sage: mut_type = QuiverMutationType(['A',3],['B',3],['X',6]); mut_type
            [ ['A', 3], ['B', 3], ['X', 6] ]
            sage: mut_type.rank()
            12
        """
        return self._rank

    @cached_method
    def b_matrix(self):
        """
        Return the B-matrix of the standard quiver of ``self``.

        The conventions for B-matrices agree with Fomin-Zelevinsky (up
        to a reordering of the simple roots).

        EXAMPLES::

            sage: mut_type = QuiverMutationType( ['A',5] ); mut_type
            ['A', 5]
            sage: mut_type.b_matrix()
            [ 0  1  0  0  0]
            [-1  0 -1  0  0]
            [ 0  1  0  1  0]
            [ 0  0 -1  0 -1]
            [ 0  0  0  1  0]

            sage: mut_type = QuiverMutationType(['A',3],['B',3]); mut_type
            [ ['A', 3], ['B', 3] ]
            sage: mut_type.b_matrix()
            [ 0  1  0  0  0  0]
            [-1  0 -1  0  0  0]
            [ 0  1  0  0  0  0]
            [ 0  0  0  0  1  0]
            [ 0  0  0 -1  0 -1]
            [ 0  0  0  0  2  0]
        """
        return _edge_list_to_matrix(self._digraph.edges(), list(range(self._rank)), [])

    @cached_method
    def standard_quiver(self):
        """
        Return the standard quiver of ``self``.

        EXAMPLES::

            sage: mut_type = QuiverMutationType( ['A',5] ); mut_type
            ['A', 5]
            sage: mut_type.standard_quiver()
            Quiver on 5 vertices of type ['A', 5]

            sage: mut_type = QuiverMutationType( ['A',[5,3], 1] ); mut_type
            ['A', [3, 5], 1]
            sage: mut_type.standard_quiver()
            Quiver on 8 vertices of type ['A', [3, 5], 1]

            sage: mut_type = QuiverMutationType(['A',3],['B',3]); mut_type
            [ ['A', 3], ['B', 3] ]
            sage: mut_type.standard_quiver()
            Quiver on 6 vertices of type [ ['A', 3], ['B', 3] ]

            sage: mut_type = QuiverMutationType(['A',3],['B',3],['X',6]); mut_type
            [ ['A', 3], ['B', 3], ['X', 6] ]
            sage: mut_type.standard_quiver()
            Quiver on 12 vertices of type [ ['A', 3], ['B', 3], ['X', 6] ]
        """
        from .quiver import ClusterQuiver
        Q = ClusterQuiver(self._digraph)
        Q._mutation_type = self
        return Q

    @cached_method
    def cartan_matrix(self):
        """
        Return the Cartan matrix of ``self``.

        Note that (up to a reordering of the simple roots) the convention for
        the definition of Cartan matrix, used here and elsewhere in Sage,
        agrees with the conventions of Kac, Fulton-Harris, and
        Fomin-Zelevinsky, but disagrees with the convention of Bourbaki.
        The `(i,j)` entry is `2(\\alpha_i,\\alpha_j)/(\\alpha_i,\\alpha_i)`.

        EXAMPLES::

            sage: mut_type = QuiverMutationType(['A',5]); mut_type
            ['A', 5]
            sage: mut_type.cartan_matrix()
            [ 2 -1  0  0  0]
            [-1  2 -1  0  0]
            [ 0 -1  2 -1  0]
            [ 0  0 -1  2 -1]
            [ 0  0  0 -1  2]

            sage: mut_type = QuiverMutationType(['A',3],['B',3]); mut_type
            [ ['A', 3], ['B', 3] ]
            sage: mut_type.cartan_matrix()
            [ 2 -1  0  0  0  0]
            [-1  2 -1  0  0  0]
            [ 0 -1  2  0  0  0]
            [ 0  0  0  2 -1  0]
            [ 0  0  0 -1  2 -1]
            [ 0  0  0  0 -2  2]
        """
        # as soon as CartanMatrix is implemented we should use it here:
        # from sage.combinat.root_system.cartan_matrix import CartanMatrix
        cmat = copy(self.b_matrix())
        for i, j in cmat.nonzero_positions():
            a = cmat[i, j]
            if a > 0:
                cmat[i, j] = -a
        for i in range(self._rank):
            cmat[i, i] = 2
        # return CartanMatrix(cmat)
        return cmat

    def is_irreducible(self):
        """
        Return ``True`` if ``self`` is irreducible.

        EXAMPLES::

            sage: mt = QuiverMutationType(['A', 2])
            sage: mt.is_irreducible()
            True
        """
        return self._info['irreducible']

    def is_mutation_finite(self):
        """
        Return ``True`` if ``self`` is of finite mutation type.

        This means that its mutation class has only finitely many
        different B-matrices.

        EXAMPLES::

            sage: mt = QuiverMutationType(['D',5, 1])
            sage: mt.is_mutation_finite()
            True
        """
        return self._info['mutation_finite']

    def is_simply_laced(self):
        """
        Return ``True`` if ``self`` is simply laced.

        This means that the only arrows that appear in the quiver of
        ``self`` are single unlabelled arrows.

        EXAMPLES::

            sage: mt = QuiverMutationType(['A', 2])
            sage: mt.is_simply_laced()
            True

            sage: mt = QuiverMutationType(['B', 2])
            sage: mt.is_simply_laced()
            False

            sage: mt = QuiverMutationType(['A',(1, 1), 1])
            sage: mt.is_simply_laced()
            False
        """
        return self._info['simply_laced']

    def is_skew_symmetric(self):
        """
        Return ``True`` if the B-matrix of ``self`` is skew-symmetric.

        EXAMPLES::

            sage: mt = QuiverMutationType(['A', 2])
            sage: mt.is_skew_symmetric()
            True

            sage: mt = QuiverMutationType(['B', 2])
            sage: mt.is_skew_symmetric()
            False

            sage: mt = QuiverMutationType(['A',(1, 1), 1])
            sage: mt.is_skew_symmetric()
            True
        """
        return self._info['skew_symmetric']

    def is_finite(self):
        """
        Return ``True`` if ``self`` is of finite type.

        This means that the cluster algebra associated to ``self`` has
        only a finite number of cluster variables.

        EXAMPLES::

            sage: mt = QuiverMutationType(['A', 2])
            sage: mt.is_finite()
            True

            sage: mt = QuiverMutationType(['A',[4, 2], 1])
            sage: mt.is_finite()
            False
        """
        return self._info['finite']

    def is_affine(self):
        """
        Return ``True`` if ``self`` is of affine type.

        EXAMPLES::

            sage: mt = QuiverMutationType(['A', 2])
            sage: mt.is_affine()
            False

            sage: mt = QuiverMutationType(['A',[4, 2], 1])
            sage: mt.is_affine()
            True
        """
        if self.is_irreducible():
            return self._info['affine']
        else:
            return False

    def is_elliptic(self):
        """
        Return ``True`` if ``self`` is of elliptic type.

        EXAMPLES::

            sage: mt = QuiverMutationType(['A', 2])
            sage: mt.is_elliptic()
            False

            sage: mt = QuiverMutationType(['E',6,[1, 1]])
            sage: mt.is_elliptic()
            True
        """
        if self.is_irreducible():
            return self._info['elliptic']
        else:
            return False

    def properties(self):
        """
        Print a scheme of all properties of ``self``.

        Most properties have natural definitions for either irreducible or
        reducible types.  ``affine`` and ``elliptic`` are only defined for
        irreducible types.

        EXAMPLES::

            sage: mut_type = QuiverMutationType(['A',3]); mut_type
            ['A', 3]
            sage: mut_type.properties()
            ['A', 3] has rank 3 and the following properties:
                - irreducible:       True
                - mutation finite:   True
                - simply-laced:      True
                - skew-symmetric:    True
                - finite:            True
                - affine:            False
                - elliptic:          False

            sage: mut_type = QuiverMutationType(['B',3]); mut_type
            ['B', 3]
            sage: mut_type.properties()
            ['B', 3] has rank 3 and the following properties:
                - irreducible:       True
                - mutation finite:   True
                - simply-laced:      False
                - skew-symmetric:    False
                - finite:            True
                - affine:            False
                - elliptic:          False

            sage: mut_type = QuiverMutationType(['B',3, 1]); mut_type
            ['BD', 3, 1]
            sage: mut_type.properties()
            ['BD', 3, 1] has rank 4 and the following properties:
                - irreducible:       True
                - mutation finite:   True
                - simply-laced:      False
                - skew-symmetric:    False
                - finite:            False
                - affine:            True
                - elliptic:          False

            sage: mut_type = QuiverMutationType(['E',6,[1, 1]]); mut_type
            ['E', 6, [1, 1]]
            sage: mut_type.properties()
            ['E', 6, [1, 1]] has rank 8 and the following properties:
                - irreducible:       True
                - mutation finite:   True
                - simply-laced:      False
                - skew-symmetric:    True
                - finite:            False
                - affine:            False
                - elliptic:          True

            sage: mut_type = QuiverMutationType(['A',3],['B',3]); mut_type
            [ ['A', 3], ['B', 3] ]
            sage: mut_type.properties()
            [ ['A', 3], ['B', 3] ] has rank 6 and the following properties:
                - irreducible:       False
                - mutation finite:   True
                - simply-laced:      False
                - skew-symmetric:    False
                - finite:            True

            sage: mut_type = QuiverMutationType('GR',[4,9]); mut_type
            ['GR', [4, 9]]
            sage: mut_type.properties()
            ['GR', [4, 9]] has rank 12 and the following properties:
                - irreducible:       True
                - mutation finite:   False
                - simply-laced:      True
                - skew-symmetric:    True
                - finite:            False
                - affine:            False
                - elliptic:          False
        """
        txt = '{} has rank {} and the following properties:'
        print(txt.format(self, self.rank()))
        s = "\t- {} {}"
        print(s.format('irreducible:      ', self.is_irreducible()))
        print(s.format('mutation finite:  ', self.is_mutation_finite()))
        print(s.format('simply-laced:     ', self.is_simply_laced()))
        print(s.format('skew-symmetric:   ', self.is_skew_symmetric()))
        print(s.format('finite:           ', self.is_finite()))
        if self.is_irreducible():
            print(s.format('affine:           ', self.is_affine()))
            print(s.format('elliptic:         ', self.is_elliptic()))


class QuiverMutationType_Irreducible(QuiverMutationType_abstract):
    """
    The mutation type for a cluster algebra or a quiver. Should not be
    called directly, but through QuiverMutationType.
    """

    def __init__(self, letter, rank, twist=None):
        """
        Should not be called directly but through QuiverMutationType.

        INPUT:

        - ``letter`` -- the letter of the mutation type
        - ``rank`` -- the rank of the mutation type
        - ``twist`` -- the twist of the mutation type

        EXAMPLES::

            sage: QuiverMutationType('A',5)
            ['A', 5]

            sage: QuiverMutationType('A',[4,5], 1)
            ['A', [4, 5], 1]

            sage: QuiverMutationType('BB',5, 1)
            ['BB', 5, 1]

            sage: QuiverMutationType('X',6)
            ['X', 6]
        """
        # _rank and _bi_rank are initialized
        self._rank = None
        self._bi_rank = None

        # _graph and _digraph are initialized
        self._graph = Graph()
        self._digraph = DiGraph()

        # _info is initialized
        self._info = {}
        self._info['irreducible'] = True
        self._info['mutation_finite'] = False
        self._info['simply_laced'] = False
        self._info['skew_symmetric'] = False
        self._info['finite'] = False
        self._info['affine'] = False
        self._info['elliptic'] = False
        self._info['irreducible_components'] = False

        if isinstance(rank, tuple):
            rank = list(rank)
        if isinstance(twist, tuple):
            twist = list(twist)

        # _letter/twist is the input letter/twist
        self._letter = letter
        self._twist = twist

        data = [letter, rank, twist]

        # type A (finite and affine)
        if letter == 'A':
            if twist is None and rank in ZZ and rank > 0:
                self._rank = rank
                self._info['mutation_finite'] = True
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                self._info['finite'] = True
            elif twist == 1 and isinstance(rank, list) and len(rank) == 2 and all(ri in ZZ and ri >= 0 for ri in rank) and rank != [0, 0]:
                if isinstance(rank, tuple):
                    rank = list(rank)
                    data[1] = rank
                rank = sorted(rank)
                self._bi_rank = rank
                self._rank = sum(self._bi_rank)
                self._info['mutation_finite'] = True
                if self._rank > 2:
                    self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                if rank[0] > 0:
                    self._info['affine'] = True
                elif rank[0] == 0:
                    self._info['finite'] = True
            else:
                _mutation_type_error(data)
            # types ['A', 1] and ['A',[0, 1], 1] need to be treated on
            # itself (as there is no edge)
            if twist is None and self._rank == 1 or twist == 1 and self._rank == 1:
                self._graph.add_vertex(0)
            # type ['A',[1, 1], 1] needs to be treated on itself as well
            # (as there is a double edge)
            elif twist == 1 and self._bi_rank[0] == 1 and self._bi_rank[1] == 1:
                self._graph.add_edge(0, 1, 2)
            else:
                for i in range(self._rank - 1):
                    self._graph.add_edge(i, i + 1, 1)
                if twist == 1:
                    self._digraph.add_edge(self._rank - 1, 0, 1)
                    for i in range(self._rank - 1):
                        if i < (2 * self._bi_rank[0]) and i % 2 == 0:
                            self._digraph.add_edge(i + 1, i, 1)
                        else:
                            self._digraph.add_edge(i, i + 1, 1)

        # type B (finite)
        elif letter == 'B':
            if twist is None and rank in ZZ and rank > 1:
                self._rank = rank
                self._info['mutation_finite'] = True
                self._info['finite'] = True
            else:
                _mutation_type_error(data)
            for i in range(rank - 2):
                self._graph.add_edge(i, i+1, 1)
            if (rank % 2 == 0):
                self._graph.add_edge(rank-2, rank-1, (1, -2))
            else:
                self._graph.add_edge(rank-2, rank-1, (2, -1))

        # type C (finite)
        elif letter == 'C':
            if twist is None and rank in ZZ and rank > 1:
                self._rank = rank
                self._info['mutation_finite'] = True
                self._info['finite'] = True
            else:
                _mutation_type_error(data)
            for i in range(rank - 2):
                self._graph.add_edge(i, i+1, 1)
            if (rank % 2 == 0):
                self._graph.add_edge(rank-2, rank-1, (2, -1))
            else:
                self._graph.add_edge(rank-2, rank-1, (1, -2))

        # type BB (affine)
        elif letter == 'BB':
            if twist == 1 and rank in ZZ and rank > 1:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
            else:
                _mutation_type_error(data)
            for i in range(rank - 2):
                self._graph.add_edge(i, i+1, 1)
            if rank % 2 == 0:
                self._graph.add_edge(rank-2, rank-1, (1, -2))
            else:
                self._graph.add_edge(rank-2, rank-1, (2, -1))
            self._graph.add_edge(rank, 0, (1, -2))

        # type CC (affine)
        elif letter == 'CC':
            if twist == 1 and rank in ZZ and rank > 1:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
            else:
                _mutation_type_error(data)
            for i in range(rank - 2):
                self._graph.add_edge(i, i+1, 1)
            if rank % 2 == 0:
                self._graph.add_edge(rank-2, rank-1, (2, -1))
            else:
                self._graph.add_edge(rank-2, rank-1, (1, -2))
            self._graph.add_edge(rank, 0, (2, -1))

        # type BC (affine)
        elif letter == 'BC':
            if twist == 1 and rank in ZZ and rank >= 1:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
            else:
                _mutation_type_error(data)
            if rank == 1:
                self._graph.add_edge(0, 1, (1, -4))
            else:
                for i in range(rank - 2):
                    self._graph.add_edge(i, i+1, 1)
                if (rank % 2 == 0):
                    self._graph.add_edge(rank-2, rank-1, (2, -1))
                else:
                    self._graph.add_edge(rank-2, rank-1, (1, -2))
                if twist == 1:
                    self._graph.add_edge(rank, 0, (1, -2))

        # type BD (affine)
        elif letter == 'BD':
            if twist == 1 and rank in ZZ and rank > 2:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
            else:
                _mutation_type_error(data)
            for i in range(rank - 2):
                self._graph.add_edge(i, i+1, 1)
            if (rank % 2 == 0):
                self._graph.add_edge(rank-2, rank-1, (1, -2))
            else:
                self._graph.add_edge(rank-2, rank-1, (2, -1))
            if twist == 1:
                self._graph.add_edge(rank, 1, 1)

        # type CD (affine)
        elif letter == 'CD':
            if twist == 1 and rank in ZZ and rank > 2:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
            else:
                _mutation_type_error(data)
            for i in range(rank - 2):
                self._graph.add_edge(i, i+1, 1)
            if (rank % 2 == 0):
                self._graph.add_edge(rank-2, rank-1, (2, -1))
            else:
                self._graph.add_edge(rank-2, rank-1, (1, -2))
            if twist == 1:
                self._graph.add_edge(rank, 1, 1)

        # type D (finite and affine)
        elif letter == 'D':
            if rank in ZZ and rank > 3 and twist is None:
                self._rank = rank
                self._info['mutation_finite'] = True
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                self._info['finite'] = True
            elif twist == 1 and rank in ZZ and rank > 3:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                self._info['affine'] = True
            else:
                _mutation_type_error(data)
            for i in range(rank - 2):
                self._graph.add_edge(i, i+1, 1)

            self._graph.add_edge(rank-3, rank-1, 1)
            if twist is not None:
                self._graph.add_edge(rank, 1, 1)

        # type E (finite, affine and elliptic)
        elif letter == 'E':
            if rank in [6, 7, 8] and twist is None:
                self._rank = rank
                self._info['mutation_finite'] = True
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                self._info['finite'] = True
                if rank == 6:
                    self._graph.add_edges([(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)])
                elif rank == 7:
                    self._graph.add_edges([(0, 1), (1, 2), (2, 3),
                                           (3, 4), (4, 5), (2, 6)])
                elif rank == 8:
                    self._graph.add_edges([(0, 1), (1, 2), (2, 3),
                                           (3, 4), (4, 5), (5, 6), (2, 7)])
            elif rank in [6, 7, 8] and twist == 1:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                self._info['affine'] = True
                if rank == 6:
                    self._graph.add_edges([(0, 1), (1, 2), (2, 3), (3, 4), (2, 5), (5, 6)])
                elif rank == 7:
                    self._graph.add_edges([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (3, 7)])
                elif rank == 8:
                    self._graph.add_edges([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (2, 8)])
            elif rank in [6, 7, 8] and twist == [1, 1]:
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['skew_symmetric'] = True
                self._info['elliptic'] = True
                if rank == 6:
                    self._digraph.add_edges([(0, 1, 1), (1, 2, 1), (3, 2, 1), (3, 4, 1), (5, 6, 1), (6, 7, 1), (5, 1, 1), (2, 5, 2), (5, 3, 1), (6, 2, 1)])
                elif rank == 7:
                    self._digraph.add_edges([(1, 0, 1), (1, 2, 1), (2, 3, 1), (4, 3, 1), (4, 5, 1),
                                             (6, 5, 1), (7, 8, 1), (3, 7, 2), (7, 2, 1), (7, 4, 1), (8, 3, 1)])
                elif rank == 8:
                    self._digraph.add_edges([(0, 1, 1), (1, 9, 1), (3, 9, 1), (3, 4, 1), (2, 8, 1), (2, 1, 1),
                                             (9, 2, 2), (2, 3, 1), (8, 9, 1), (5, 4, 1), (5, 6, 1), (7, 6, 1)])
            # type E (mutation infinite)
            elif rank > 9 and twist is None:
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                self._rank = rank
                for i in range(rank-2):
                    self._graph.add_edge(i, i+1, 1)
                self._graph.add_edge(2, rank-1)
            else:
                _mutation_type_error(data)

        # type AE (mutation infinite)
        elif letter == 'AE':
            if isinstance(rank, list) and len(rank) == 2 and all(rank[i] in ZZ and rank[i] > 0 for i in [0, 1]) and twist is None:
                if isinstance(rank, tuple):
                    rank = list(rank)
                    data[1] = rank
                rank = sorted(rank)
                self._bi_rank = rank
                self._rank = sum(self._bi_rank) + 1
                if self._rank > 3:
                    self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                if self._bi_rank == [1, 1]:
                    self._graph.add_edges([(0, 1, 2), (1, 2, None)])
                else:
                    self._digraph.add_edge(self._rank - 2, 0)
                    for i in range(self._rank-2):
                        if i < (2 * self._bi_rank[0]) and i % 2 == 0:
                            self._digraph.add_edge(i + 1, i)
                        else:
                            self._digraph.add_edge(i, i + 1)
                    self._digraph.add_edge(self._rank-2, self._rank-1)
            else:
                _mutation_type_error(data)

        # type BE (mutation infinite)
        elif letter == 'BE':
            if rank > 4 and twist is None:
                self._rank = rank
                for i in range(rank-3):
                    self._graph.add_edge(i, i+1)
                self._graph.add_edge(2, rank-1)
                if rank % 2 == 0:
                    self._graph.add_edge(rank-3, rank-2, (2, -1))
                else:
                    self._graph.add_edge(rank-3, rank-2, (1, -2))
            else:
                _mutation_type_error(data)

        # type CE (mutation infinite)
        elif letter == 'CE':
            if rank > 4 and twist is None:
                self._rank = rank
                for i in range(rank-3):
                    self._graph.add_edge(i, i+1)
                self._graph.add_edge(2, rank-1)
                if rank % 2 == 0:
                    self._graph.add_edge(rank-3, rank-2, (1, -2))
                else:
                    self._graph.add_edge(rank-3, rank-2, (2, -1))
            else:
                _mutation_type_error(data)

        # type DE (mutation infinite)
        elif letter == 'DE':
            if rank > 5 and twist is None:
                self._rank = rank
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                for i in range(rank-3):
                    self._graph.add_edge(i, i+1)
                self._graph.add_edge(2, rank-2)
                self._graph.add_edge(rank-4, rank-1)
            else:
                _mutation_type_error(data)

        # type F (finite, affine, and elliptic)
        elif letter == 'F':
            if rank == 4 and twist is None:
                self._rank = rank
                self._info['mutation_finite'] = True
                self._info['finite'] = True
                self._graph.add_edges([(0, 1, None), (1, 2, (2, -1)), (2, 3, None)])
            elif rank == 4 and twist == 1:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
                self._graph.add_edges([(0, 1, None), (1, 2, None),
                                       (2, 3, (1, -2)), (3, 4, None)])
            elif rank == 4 and twist == -1:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
                self._graph.add_edges([(0, 1, None), (1, 2, None),
                                       (2, 3, (2, -1)), (3, 4, None)])
            elif rank == 4 and (twist == [1, 2]):
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['elliptic'] = True
                self._digraph.add_edges([(0, 1, None), (1, 2, None),
                                         (2, 3, (2, -1)), (4, 2, (1, -2)),
                                         (3, 4, 2), (4, 5, None), (5, 3, None)])
            elif rank == 4 and (twist == [2, 1]):
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['elliptic'] = True
                self._digraph.add_edges([(0, 1, None), (1, 2, None),
                                         (2, 3, (1, -2)), (4, 2, (2, -1)),
                                         (3, 4, 2), (4, 5, None), (5, 3, None)])
            elif rank == 4 and twist == [2, 2]:
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['elliptic'] = True
                self._digraph.add_edges([(0, 1, None), (1, 2, None),
                                         (3, 1, None), (2, 3, 2),
                                         (4, 2, (2, -1)), (3, 4, (1, -2)),
                                         (5, 4, None)])
            elif rank == 4 and twist == [1, 1]:
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['elliptic'] = True
                self._digraph.add_edges([(0, 1, None), (1, 2, None),
                                         (3, 1, None), (2, 3, 2), (4, 2, (1, -2)),
                                         (3, 4, (2, -1)), (5, 4, None)])
            else:
                _mutation_type_error(data)

        # type G (finite, affine, and elliptic)
        elif letter == 'G':
            if rank == 2 and twist is None:
                self._rank = rank
                self._info['mutation_finite'] = True
                self._info['finite'] = True
                self._graph.add_edges([(0, 1, (1, -3))])
            elif rank == 2 and twist == -1:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
                self._graph.add_edges([(0, 1, None), (1, 2, (1, -3))])
            elif rank == 2 and twist == 1:
                self._rank = rank + 1
                self._info['mutation_finite'] = True
                self._info['affine'] = True
                self._graph.add_edges([(0, 1, None), (1, 2, (3, -1))])
            elif rank == 2 and (twist == [1, 3]):
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['elliptic'] = True
                self._digraph.add_edges([(0, 1, None), (1, 2, (3, -1)),
                                         (3, 1, (1, -3)), (2, 3, 2)])
            elif rank == 2 and (twist == [3, 1]):
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['elliptic'] = True
                self._digraph.add_edges([(0, 1, None), (1, 2, (1, -3)),
                                         (3, 1, (3, -1)), (2, 3, 2)])
            elif rank == 2 and twist == [3, 3]:
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['elliptic'] = True
                self._digraph.add_edges([(1, 0, None), (0, 2, 2), (3, 0, (3, -1)),
                                         (2, 1, None), (2, 3, (1, -3))])
            elif rank == 2 and twist == [1, 1]:
                self._rank = rank + 2
                self._info['mutation_finite'] = True
                self._info['elliptic'] = True
                self._digraph.add_edges([(1, 0, None), (0, 2, 2), (3, 0, (1, -3)),
                                         (2, 1, None), (2, 3, (3, -1))])
            else:
                _mutation_type_error(data)

        # type GR (mutation infinite)
        elif letter == 'GR':
            if twist is None and isinstance(rank, list) and len(rank) == 2 and all(rank[i] in ZZ and rank[i] > 0 for i in [0, 1]) and rank[1] - 1 > rank[0] > 1:
                gr_rank = (rank[0]-1, rank[1]-rank[0]-1)
                self._rank = prod(gr_rank)
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                a, b = gr_rank
                for i in range(a):
                    for j in range(b):
                        if i < a-1:
                            if (i+j) % 2 == 0:
                                self._digraph.add_edge(i*b+j, (i+1)*b+j)
                            else:
                                self._digraph.add_edge((i+1)*b+j, i*b+j)
                        if j < b-1:
                            if (i+j) % 2 == 0:
                                self._digraph.add_edge(i*b+j+1, i*b+j)
                            else:
                                self._digraph.add_edge(i*b+j, i*b+j+1)
            else:
                _mutation_type_error(data)

        # type R2 (rank 2 finite mutation types)
        elif letter == 'R2':
            if twist is None and isinstance(rank, list) and len(rank) == 2 and all(rank[i] in ZZ and rank[i] > 0 for i in [0, 1]):
                rank = sorted(rank)
                b, c = rank
                self._rank = 2
                if b == c:
                    self._info['skew_symmetric'] = True
                self._graph.add_edge(0, 1, (b, -c))
            else:
                _mutation_type_error(data)

        # type T
        elif letter == 'T':
            if twist is None and isinstance(rank, list) and len(rank) == 3 and all(rank[i] in ZZ and rank[i] > 0 for i in [0, 1, 2]):
                if isinstance(rank, tuple):
                    rank = list(rank)
                    data[1] = rank
                rank = sorted(rank)
                self._rank = sum(rank) - 2
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                r, p, q = rank
                for i in range(q-1):
                    if i == 0:
                        self._graph.add_edge(0, 1)
                        self._graph.add_edge(0, r)
                        self._graph.add_edge(0, r+p-1)
                    else:
                        if i < r-1:
                            self._graph.add_edge(i, i+1)
                        if i < p-1:
                            self._graph.add_edge(i+r-1, i+r)
                        self._graph.add_edge(i+r+p-2, i+r+p-1)
            else:
                _mutation_type_error(data)

        # type TR (mutation infinite if rank > 2)
        elif letter == 'TR':
            # type ['TR', 1] needs to be treated on itself (as there is no edge)
            if twist is None and rank == 1:
                self._graph.add_vertex(0)
            elif twist is None and rank > 1:
                self._rank = rank*(rank+1)//2
                self._info['simply_laced'] = True
                self._info['skew_symmetric'] = True
                level = 0
                while level < rank:
                    nr = rank*level-sum(range(level))
                    for i in range(nr, nr+rank-level-1):
                        self._digraph.add_edge(i, i+1)
                        self._digraph.add_edge(i+rank-level, i)
                        self._digraph.add_edge(i+1, i+rank-level)
                    level += 1
            else:
                _mutation_type_error(data)

        # type X
        elif letter == 'X':
            if rank in [6, 7] and twist is None:
                self._rank = rank
                self._info['mutation_finite'] = True
                self._info['skew_symmetric'] = True
                self._digraph.add_edges([(0, 1, 2), (1, 2, None), (2, 0, None),
                                         (2, 3, None), (3, 4, 2), (4, 2, None),
                                         (2, 5, None)])
                if rank == 7:
                    self._digraph.add_edges([(5, 6, 2), (6, 2, None)])
            else:
                _mutation_type_error(data)

        # otherwise, an error is raised
        else:
            _mutation_type_error(data)

        # in the bipartite case, the digraph is constructed from the graph
        if not self._digraph:
            if self._graph.is_bipartite():
                self._digraph = _bipartite_graph_to_digraph(self._graph)
            else:
                raise ValueError('The QuiverMutationType does not have '
                                 'a Coxeter diagram.')

        # in the other cases, the graph is constructed from the digraph
        if not self._graph:
            self._graph = self._digraph.to_undirected()

        # _description is as for CartanType
        if twist:
            self._description = str([letter, rank, twist])
        else:
            self._description = str([letter, rank])

    def irreducible_components(self):
        """
        Return a list of all irreducible components of ``self``.

        EXAMPLES::

            sage: mut_type = QuiverMutationType('A',3); mut_type
            ['A', 3]
            sage: mut_type.irreducible_components()
            (['A', 3],)
        """
        return tuple([self])

    @cached_method
    def class_size(self):
        r"""
        If it is known, the size of the mutation class of all quivers
        which are mutation equivalent to the standard quiver of
        ``self`` (up to isomorphism) is returned.

        Otherwise, ``NotImplemented`` is returned.

        Formula for finite type A is taken from Torkildsen - Counting
        cluster-tilted algebras of type `A_n`.
        Formulas for affine type A and finite type D are taken from Bastian,
        Prellberg, Rubey, Stump - Counting the number of elements in the
        mutation classes of `\widetilde A_n` quivers.
        Formulas for finite and affine types B and C are
        proven but not yet published.
        Conjectural formulas for several other non-simply-laced affine types
        are implemented.
        Exceptional Types (finite, affine, and elliptic) E, F, G, and X are
        hardcoded.

        EXAMPLES::

            sage: mut_type = QuiverMutationType( ['A',5] ); mut_type
            ['A', 5]
            sage: mut_type.class_size()
            19

            sage: mut_type = QuiverMutationType( ['A',[10,3], 1] ); mut_type
            ['A', [3, 10], 1]
            sage: mut_type.class_size()
            142120

            sage: mut_type = QuiverMutationType( ['B',6] ); mut_type
            ['B', 6]
            sage: mut_type.class_size()
            132

            sage: mut_type = QuiverMutationType( ['BD',6, 1] ); mut_type
            ['BD', 6, 1]
            sage: mut_type.class_size()
            Warning: This method uses a formula which has not been proved correct.
            504

        Check that :trac:`14048` is fixed::

            sage: mut_type = QuiverMutationType( ['F',4,(2, 1)] )
            sage: mut_type.class_size()
            90
        """
        if not self.is_mutation_finite():
            return infinity

        # type A (finite and affine)
        if self._letter == 'A':
            # the formula is taken from Torkildsen - Counting
            # cluster-tilted algebras of type A
            if self.is_finite():
                n = self._rank
                a = binomial(2*(n+1), n+1) // (n+2)
                if n % 2 == 1:
                    a += binomial(n+1, (n+1)//2)
                if n % 3 == 0:
                    a += 2 * binomial(2*n//3, n//3)
                return a // (n+3)
            # the formula is taken from Bastian, Prellberg, Rubey, Stump
            elif self.is_affine():
                i, j = self._bi_rank
                i = ZZ(i)
                j = ZZ(j)
                n = i+j
                f = euler_phi
                if i == j:
                    return (binomial(2 * i, i) +
                            sum(f(k) * binomial(2 * i // k, i // k)**2
                                for k in i.divisors()
                                if k in j.divisors()) // n) // 4
                else:
                    return sum(f(k) * binomial(2 * i // k, i // k) *
                               binomial(2 * j // k, j // k)
                               for k in i.divisors()
                               if k in j.divisors()) // (2 * n)

        # types B and C (finite and affine)
        elif self._letter in ['B', 'C']:
            # this formula is proven but nowhere published correctness
            # is clear enough that I don't think a warning is needed
            if self.is_finite():
                n = self._rank
                return binomial(2 * n, n) // (n + 1)

        elif self._letter in ['BB', 'CC']:
            # these two formulas are not yet proven
            print("Warning: This method uses a formula "
                  "which has not been proved correct.")
            if self.is_affine():
                if self._twist == 1:
                    n = self._rank - 1
                    N = binomial(2 * n - 1, n - 1)
                    if n % 2:
                        return N
                    else:
                        return N + binomial(n - 1, n // 2 - 1)

        # type BC (affine)
        elif self._letter == 'BC':
            # this formula is not yet proven
            print("Warning: This method uses a formula "
                  "which has not been proved correct.")
            if self.is_affine():
                if self._twist == 1:
                    n = self._rank - 1
                    return binomial(2 * n, n)

        # types BD and CD (affine)
        elif self._letter in ['BD', 'CD']:
            # this formula is not yet proven
            print("Warning: This method uses a formula "
                  "which has not been proved correct.")
            if self.is_affine():
                if self._twist == 1:
                    n = self._rank - 2
                    return 2 * binomial(2 * n, n)

        # type D (finite and affine)
        elif self._letter == 'D':
            # the formula is taken from Bastian, Prellberg, Rubey, Stump
            if self.is_finite():
                if self._rank == 4:
                    return 6
                else:
                    f = euler_phi
                    n = ZZ(self._rank)
                    return sum(f(n // k) * binomial(2 * k, k)
                               for k in n.divisors()) // (2 * n)
            # this formula is not yet proven
            elif self.is_affine():
                n = self._rank - 3
                if n == 2:
                    return 9
                else:
                    print("Warning: This method uses a formula "
                          "which has not been proved correct.")
                    if n % 2:
                        return 2 * binomial(2 * n, n)
                    else:
                        return 2 * binomial(2 * n, n) + binomial(n, n // 2)

        # the exceptional types are hard-coded
        # type E (finite, affine and elliptic)
        elif self._letter == 'E':
            if self.is_finite():
                if self._rank == 6:
                    return 67
                elif self._rank == 7:
                    return 416
                elif self._rank == 8:
                    return 1574
            elif self.is_affine():
                if self._rank == 7:
                    return 132
                elif self._rank == 8:
                    return 1080
                elif self._rank == 9:
                    return 7560
            elif self.is_elliptic():
                if self._rank == 8:
                    return 49
                elif self._rank == 9:
                    return 506
                elif self._rank == 10:
                    return 5739

        # type F
        elif self._letter == 'F':
            if self.is_finite():
                return 15
            elif self.is_affine():
                return 60
            elif self.is_elliptic():
                if self._twist == [1, 2]:
                    return 90
                if self._twist == [1, 1] or self._twist == [2, 2]:
                    return 35

        # type G
        elif self._letter == 'G':
            if self.is_finite():
                return 2
            elif self.is_affine():
                return 6
            elif self.is_elliptic():
                if self._twist == [1, 3]:
                    return 7
                if self._twist == [1, 1] or self._twist == [3, 3]:
                    return 2

        # type X
        elif self._letter == 'X':
            if self._rank == 6:
                return 5
            elif self._rank == 7:
                return 2

        # otherwise the size is returned to be unknown
        else:
            print("Size unknown")
            return NotImplemented

    def dual(self):
        """
        Return the QuiverMutationType which is dual to ``self``.

        EXAMPLES::

            sage: mut_type = QuiverMutationType('A',5); mut_type
            ['A', 5]
            sage: mut_type.dual()
            ['A', 5]

            sage: mut_type = QuiverMutationType('B',5); mut_type
            ['B', 5]
            sage: mut_type.dual()
            ['C', 5]
            sage: mut_type.dual().dual()
            ['B', 5]
            sage: mut_type.dual().dual() == mut_type
            True
        """
        letter = self.letter()
        # the self-dual cases
        if letter != 'BC' and letter[0] in ['B', 'C']:
            if letter == 'BB':
                letter = 'CC'
            elif letter == 'CC':
                letter = 'BB'
            elif letter[0] == 'B':
                letter = 'C' + letter[1:]
            elif letter[0] == 'C':
                letter = 'B' + letter[1:]
            rank = self._rank
            if self.is_affine():
                rank -= 1
            twist = self._twist
            return QuiverMutationType(letter, rank, twist)
        # the cases F and G have non-trivial duality in some cases
        elif letter in ['F', 'G']:
            if self.is_finite():
                return self
            elif self.is_affine():
                rank = self._rank - 1
                twist = - self._twist
            elif self.is_elliptic():
                twist = self._twist
                rank = self._rank - 2
                if letter == 'F':
                    if self._twist == [2, 2]:
                        twist == [1, 1]
                    if self._twist == [1, 1]:
                        twist == [2, 2]
                if letter == 'G':
                    if self._twist == [3, 3]:
                        twist = [1, 1]
                    elif self._twist == [1, 1]:
                        twist = [3, 3]
            else:
                rank = self._rank
            return QuiverMutationType(letter, rank, twist)
        else:
            return self


class QuiverMutationType_Reducible(QuiverMutationType_abstract):
    """
    The mutation type for a cluster algebra or a quiver. Should not be
    called directly, but through QuiverMutationType.  Inherits from
    QuiverMutationType_abstract.
    """

    def __init__(self, *args):
        """
        Should not be called directly, but through QuiverMutationType.

        INPUT:

        - ``data`` -- a list each of whose entries is a
          QuiverMutationType_Irreducible

        EXAMPLES::

            sage: QuiverMutationType(['A',4],['B',6])
            [ ['A', 4], ['B', 6] ]
        """
        data = args
        if len(data) < 2 or not all(isinstance(comp, QuiverMutationType_Irreducible) for comp in data):
            return _mutation_type_error(data)

        # _info is initialized
        self._info = {}
        self._info['irreducible'] = False
        self._info['mutation_finite'] = all(comp.is_mutation_finite()
                                            for comp in data)
        self._info['simply_laced'] = all(comp.is_simply_laced()
                                         for comp in data)
        self._info['skew_symmetric'] = all(comp.is_skew_symmetric()
                                           for comp in data)
        self._info['finite'] = all(comp.is_finite() for comp in data)
        self._info['irreducible_components'] = copy(data)

        #  letter and rank are initialized
        self._letter = ''
        self._rank = 0

        # graph and digraph are initialized
        self._graph = Graph()
        self._digraph = DiGraph()

        for comp in data:
            if self._letter:
                self._letter += ' x '
            self._letter += comp._letter
            self._rank += comp._rank
            self._graph = self._graph.disjoint_union(comp._graph,
                                                     labels='integers')
            self._digraph = self._digraph.disjoint_union(comp._digraph,
                                                         labels='integers')
        self._graph.name('')
        self._digraph.name('')

        # _description is as for CartanType
        self._description = "[ "
        comps = self.irreducible_components()
        for i in range(len(comps)):
            if i > 0:
                self._description += ", "
            self._description += comps[i]._description
        self._description += " ]"

    def irreducible_components(self):
        """
        Return a list of all irreducible components of ``self``.

        EXAMPLES::

            sage: mut_type = QuiverMutationType('A',3); mut_type
            ['A', 3]
            sage: mut_type.irreducible_components()
            (['A', 3],)

            sage: mut_type = QuiverMutationType(['A',3],['B',3]); mut_type
            [ ['A', 3], ['B', 3] ]
            sage: mut_type.irreducible_components()
            (['A', 3], ['B', 3])

            sage: mut_type = QuiverMutationType(['A',3],['B',3],['X',6])
            sage: mut_type
            [ ['A', 3], ['B', 3], ['X', 6] ]
            sage: mut_type.irreducible_components()
            (['A', 3], ['B', 3], ['X', 6])
        """
        return self._info['irreducible_components']

    @cached_method
    def class_size(self):
        """
        If it is known, the size of the mutation class of all quivers
        which are mutation equivalent to the standard quiver of
        ``self`` (up to isomorphism) is returned.

        Otherwise, ``NotImplemented`` is returned.

        EXAMPLES::

            sage: mut_type = QuiverMutationType(['A',3],['B',3]); mut_type
            [ ['A', 3], ['B', 3] ]
            sage: mut_type.class_size()
            20

            sage: mut_type = QuiverMutationType(['A',3],['B',3],['X',6])
            sage: mut_type
            [ ['A', 3], ['B', 3], ['X', 6] ]
            sage: mut_type.class_size()
            100
        """
        if not self.is_mutation_finite():
            return infinity

        components = []
        multiplicities = []
        for x in self.irreducible_components():
            if not components.count(x):
                components.append(x)
                multiplicities.append(1)
            else:
                y = components.index(x)
                multiplicities[y] += 1

        sizes = [x.class_size() for x in components]
        if NotImplemented in sizes:
            print("Size unknown")
            return NotImplemented
        else:
            return prod(binomial(sizes[i] + multiplicities[i] - 1,
                                 multiplicities[i])
                        for i in range(len(sizes)))

    def dual(self):
        """
        Return the QuiverMutationType which is dual to ``self``.

        EXAMPLES::

            sage: mut_type = QuiverMutationType(['A',5],['B',6],['C',5],['D',4]); mut_type
            [ ['A', 5], ['B', 6], ['C', 5], ['D', 4] ]
            sage: mut_type.dual()
            [ ['A', 5], ['C', 6], ['B', 5], ['D', 4] ]
        """
        comps = self.irreducible_components()
        return QuiverMutationType([comp.dual() for comp in comps])


def _construct_classical_mutation_classes(n) -> dict[tuple, list | set]:
    r"""
    Return a dict with keys being tuples representing regular
    QuiverMutationTypes of the given rank, and with values being lists
    or sets containing all mutation equivalent quivers as dig6 data.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _construct_classical_mutation_classes
        sage: rank_2_classes = _construct_classical_mutation_classes(2) # long time
        sage: for mut_class in sorted(rank_2_classes.keys(),key=str): # long time
        ....:   print("{} {}".format(mut_class, rank_2_classes[mut_class]))
        ('A', (1, 1), 1) [('AO', (((0, 1), (2, -2)),))]
        ('A', 2) [('AO', ())]
        ('B', 2) [('AO', (((0, 1), (1, -2)),)), ('AO', (((0, 1), (2, -1)),))]
        ('BC', 1, 1) [('AO', (((0, 1), (1, -4)),)),
        ('AO', (((0, 1), (4, -1)),))]
    """
    from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
    data: dict[tuple, set | list] = {}

    # finite A
    data[('A', n)] = ClusterQuiver(['A', n]).mutation_class(data_type='dig6')
    # affine A
    for j in range(1, n//2+1):
        data[('A', (n-j, j), 1)] = ClusterQuiver(['A', [n-j, j], 1]).mutation_class(data_type='dig6')
    # finite B
    if n > 1:
        data[('B', n)] = ClusterQuiver(['B', n]).mutation_class(data_type='dig6')
    # affine B
    if n > 2:
        data[('BB', n-1, 1)] = ClusterQuiver(['BB', n-1, 1]).mutation_class(data_type='dig6')
    # finite C
    if n > 2:
        data[('C', n)] = ClusterQuiver(['C', n]).mutation_class(data_type='dig6')
    # affine C
    if n > 1:
        data[('BC', n-1, 1)] = ClusterQuiver(['BC', n-1, 1]).mutation_class(data_type='dig6')
    # affine CC
    if n > 2:
        data[('CC', n-1, 1)] = ClusterQuiver(['CC', n-1, 1]).mutation_class(data_type='dig6')
    # affine BD
    if n > 3:
        data[('BD', n-1, 1)] = ClusterQuiver(['BD', n-1, 1]).mutation_class(data_type='dig6')
    # affine CD
    if n > 3:
        data[('CD', n-1, 1)] = ClusterQuiver(['CD', n-1, 1]).mutation_class(data_type='dig6')
    # finite D
    if n > 3:
        data[('D', n)] = ClusterQuiver(['D', n]).mutation_class(data_type='dig6')
    # affine D
    if n > 4:
        data[('D', n-1, 1)] = ClusterQuiver(['D', n-1, 1]).mutation_class(data_type='dig6')

    return data


def _construct_exceptional_mutation_classes(n) -> dict[tuple, list | set]:
    r"""
    Return a dict with keys being tuples representing exceptional
    QuiverMutationTypes of the given rank, and with values being lists
    or sets containing all mutation equivalent quivers as dig6 data.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _construct_exceptional_mutation_classes
        sage: rank_3_exceptional = _construct_exceptional_mutation_classes(3) # long time
        sage: for mut_class in sorted(rank_3_exceptional.keys(), key=str): # long time
        ....:   print("{} {}".format(mut_class, rank_3_exceptional[mut_class]))
        ('G', 2, -1) [('BH?', (((1, 2), (1, -3)),)),
        ('BGO', (((2, 1), (3, -1)),)), ('BW?', (((0, 1), (3, -1)),)),
        ('BP?', (((0, 1), (1, -3)),)),
        ('BP_', (((0, 1), (1, -3)), ((2, 0), (3, -1)))),
        ('BP_', (((0, 1), (3, -1)), ((1, 2), (1, -3)), ((2, 0), (2, -2))))]
        ('G', 2, 1) [('BH?', (((1, 2), (3, -1)),)),
        ('BGO', (((2, 1), (1, -3)),)), ('BW?', (((0, 1), (1, -3)),)),
        ('BP?', (((0, 1), (3, -1)),)),
        ('BKO', (((1, 0), (3, -1)), ((2, 1), (1, -3)))),
        ('BP_', (((0, 1), (2, -2)), ((1, 2), (1, -3)), ((2, 0), (3, -1))))]
    """
    from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
    data: dict[tuple, list | set] = {}
    # finite E
    if n in [6, 7, 8]:
        data[('E', n)] = ClusterQuiver(['E', n]).mutation_class(data_type='dig6')
    # affine E
    if n in [7, 8, 9]:
        data[('E', n - 1, 1)] = ClusterQuiver(['E', n - 1, 1]).mutation_class(data_type='dig6')
    # elliptic E
    if n in [8, 9, 10]:
        data[('E', n - 2, (1, 1))] = ClusterQuiver(['E', n - 2, [1, 1]]).mutation_class(data_type='dig6')
    # finite F
    if n == 4:
        data[('F', 4)] = ClusterQuiver(['F', 4]).mutation_class(data_type='dig6')
    # affine F
    if n == 5:
        data[('F', 4, 1)] = ClusterQuiver(['F', 4, 1]).mutation_class(data_type='dig6')
        data[('F', 4, -1)] = ClusterQuiver(['F', 4, -1]).mutation_class(data_type='dig6')
    # finite G
    if n == 2:
        data[('G', 2)] = ClusterQuiver(['G', 2]).mutation_class(data_type='dig6')
    # affine G
    if n == 3:
        data[('G', 2, 1)] = ClusterQuiver(['G', 2, 1]).mutation_class(data_type='dig6')
        data[('G', 2, -1)] = ClusterQuiver(['G', 2, -1]).mutation_class(data_type='dig6')
    # elliptic G
    if n == 4:
        data[('G', 2, (1, 3))] = ClusterQuiver(['G', 2, (1, 3)]).mutation_class(data_type='dig6')
        data[('G', 2, (1, 1))] = ClusterQuiver(['G', 2, (1, 1)]).mutation_class(data_type='dig6')
        data[('G', 2, (3, 3))] = ClusterQuiver(['G', 2, (3, 3)]).mutation_class(data_type='dig6')
    # X
    if n in [6, 7]:
        data[('X', n)] = ClusterQuiver(['X', n]).mutation_class(data_type='dig6')
    # elliptic F
    if n == 6:
        data[('F', 4, (1, 2))] = ClusterQuiver(['F', 4, (1, 2)]).mutation_class(data_type='dig6')
        data[('F', 4, (1, 1))] = ClusterQuiver(['F', 4, (1, 1)]).mutation_class(data_type='dig6')
        data[('F', 4, (2, 2))] = ClusterQuiver(['F', 4, (2, 2)]).mutation_class(data_type='dig6')

    return data


def _save_data_dig6(n, types='ClassicalExceptional', verbose=False):
    """
    Save all exceptional mutation classes as dig6 data into the file ``exc_classes_n.dig6`` in the folder ``DOT_SAGE``.

    TESTS::

        sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import save_quiver_data
        sage: save_quiver_data(2) # indirect doctest
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('A', 1)]
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('A', (1, 1), 1), ('A', 2), ('B', 2), ('BC', 1, 1), ('G', 2)]

        sage: save_quiver_data(2,up_to=False) # indirect doctest
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('A', (1, 1), 1), ('A', 2), ('B', 2), ('BC', 1, 1), ('G', 2)]

        sage: save_quiver_data(2,up_to=False, types='Classical') # indirect doctest
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('A', (1, 1), 1), ('A', 2), ('B', 2), ('BC', 1, 1)]

        sage: save_quiver_data(2,up_to=False, types='Exceptional') # indirect doctest
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('G', 2)]

        sage: save_quiver_data(2,up_to=False, verbose=False) # indirect doctest
    """
    data = {}
    possible_types = ['Classical', 'ClassicalExceptional', 'Exceptional']
    if types not in possible_types:
        raise ValueError('The third input must be either ClassicalExceptional'
                         ' (default), Classical, or Exceptional.')

    if types in possible_types[:2]:
        data.update(_construct_classical_mutation_classes(n))
    if types in possible_types[1:]:
        data.update(_construct_exceptional_mutation_classes(n))

    from sage.env import DOT_SAGE
    from sage.misc.misc import sage_makedirs
    types_path = os.path.join(DOT_SAGE, 'cluster_algebra_quiver')
    types_file = os.path.join(types_path, 'mutation_classes_%s.dig6' % n)
    sage_makedirs(types_path)
    from sage.misc.temporary_file import atomic_write
    with atomic_write(types_file, binary=True) as f:
        pickle.dump(data, f)
    if verbose:
        keys = sorted(data, key=str)
        print("\nThe following types are saved to file", types_file, "and will now be used to determine quiver mutation types:")
        print(keys)


def save_quiver_data(n, up_to=True, types='ClassicalExceptional', verbose=True):
    r"""
    Save mutation classes of certain quivers of ranks up to and equal
    to ``n`` or equal to ``n`` to
    ``DOT_SAGE/cluster_algebra_quiver/mutation_classes_n.dig6``.

    This data will then be used to determine quiver mutation types.

    INPUT:

    - ``n``: the rank (or the upper limit on the rank) of the mutation
      classes that are being saved.

    - ``up_to`` -- (default:``True``) if ``True``, saves data for
      ranks smaller than or equal to ``n``. If ``False``, saves data
      for rank exactly ``n``.

    - ``types`` -- (default:'ClassicalExceptional') if all, saves data
      for both exceptional mutation-finite quivers and for classical
      quiver. The input 'Exceptional' or 'Classical' is also allowed
      to save only part of this data.

    TESTS::

        sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import save_quiver_data
        sage: save_quiver_data(2)
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('A', 1)]
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('A', (1, 1), 1), ('A', 2), ('B', 2), ('BC', 1, 1), ('G', 2)]

        sage: save_quiver_data(2,up_to=False)
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('A', (1, 1), 1), ('A', 2), ('B', 2), ('BC', 1, 1), ('G', 2)]

        sage: save_quiver_data(2,up_to=False, types='Classical')
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('A', (1, 1), 1), ('A', 2), ('B', 2), ('BC', 1, 1)]

        sage: save_quiver_data(2,up_to=False, types='Exceptional')
        <BLANKLINE>
        The following types are saved to file ... and will now be used to determine quiver mutation types:
        [('G', 2)]

        sage: save_quiver_data(2,up_to=False, verbose=False)
    """
    from sage.combinat.cluster_algebra_quiver.mutation_type import load_data
    if up_to is True:
        ranks = range(1, n + 1)
    elif up_to is False:
        ranks = [n]
    for i in ranks:
        _save_data_dig6(i, types=types, verbose=verbose)
    # we finally clear the load_data
    load_data.clear_cache()


def _bipartite_graph_to_digraph(g) -> DiGraph:
    """
    Return a digraph obtained from a bipartite graph ``g`` by choosing one
    set of the bipartition to be the set of sinks and the other to be the
    set of sources.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _bipartite_graph_to_digraph
        sage: G = Graph([(1, 2)])
        sage: _bipartite_graph_to_digraph(G)
        Digraph on 2 vertices
    """
    sources = g.bipartite_sets()[0]
    dg = DiGraph()
    dg.add_vertices(g)
    for a, b, c in g.edge_iterator():
        if a in sources:
            dg.add_edge(a, b, c)
        else:
            dg.add_edge(b, a, c)
    return dg


def _is_mutation_type(data) -> bool:
    """
    Return ``True`` if ``data`` is a QuiverMutationType.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _is_mutation_type
        sage: _is_mutation_type ( [ 'A', 2 ] )
        True
        sage: _is_mutation_type ( [ 'P', 1 ] )
        False
    """
    try:
        QuiverMutationType(data)
        return True
    except Exception:
        return False


def _mutation_type_error(data):
    r"""
    Output an error message because data which is not a valid quiver mutation
    type has been passed to QuiverMutationType.

    EXAMPLES::

        sage: QuiverMutationType( 'Christian', 'Stump' ) # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: ['Christian', 'Stump'] is not a valid quiver mutation type
                Finite types have the form [ '?', n ] for type ? and rank n
                Affine type A has the form [ 'A', [ i, j ], 1 ] for rank i+j
                Affine type ? has the form [ '?', k, \pm 1 ] for rank k+1
                Elliptic type ? has the form [ '?', k, [i, j] ] (1 <= i,j <= 3) for rank k+2
                For correct syntax in other types, please consult the documentation.
    """
    if data[2] is None:
        del data[2]
    return_str = str(data) + ' is not a valid quiver mutation type'
    return_str += '\n            Finite types have the form [ \'?\', n ] for type ? and rank n'
    return_str += '\n            Affine type A has the form [ \'A\', [ i, j ], 1 ] for rank i+j'
    return_str += '\n            Affine type ? has the form [ \'?\', k, \\pm 1 ] for rank k+1'
    return_str += '\n            Elliptic type ? has the form [ \'?\', k, [i, j] ] (1 <= i,j <= 3) for rank k+2'
    return_str += '\n            For correct syntax in other types, please consult the documentation.'

    raise ValueError(return_str)


def _edge_list_to_matrix(edges, nlist, mlist) -> matrix:
    r"""
    Return the matrix obtained from the edge list of a quiver.

    INPUT:

    - ``edges`` -- the list of edges
    - ``nlist`` -- the list of mutable vertices of the quiver
    - ``mlist`` -- the list of frozen vertices of the quiver

    OUTPUT:

    An `(n+m) \times n` matrix corresponding to the edge-list.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _edge_list_to_matrix
        sage: G = QuiverMutationType(['A', 2])._digraph
        sage: _edge_list_to_matrix(G.edges(), [0, 1], [])
        [ 0  1]
        [-1  0]

        sage: G2 = DiGraph([('a', 'b', 1)])
        sage: _edge_list_to_matrix(G2.edges(), ['a', 'b'], [])
        [ 0  1]
        [-1  0]

        sage: G3 = DiGraph([('a', 'b', 1), ('b', 'c', 1)])
        sage: _edge_list_to_matrix(G3.edges(), ['a', 'b'], ['c'])
        [ 0  1]
        [-1  0]
        [ 0 -1]
    """
    n = len(nlist)
    nmlist = nlist + mlist
    nm = len(nmlist)
    M = matrix(ZZ, nm, n, sparse=True)
    for v1, v2, label in edges:
        if label is None:
            a, b = 1, -1
        elif label in ZZ:
            a, b = label, -label
        else:
            a, b = label
        if v1 in nlist:
            M[nmlist.index(v2), nmlist.index(v1)] = b
        if v2 in nlist:
            M[nmlist.index(v1), nmlist.index(v2)] = a
    return M
