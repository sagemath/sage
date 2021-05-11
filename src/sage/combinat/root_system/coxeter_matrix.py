"""
Coxeter Matrices
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#                     2015 Jean-Philippe Labbe <labbe at math.huji.ac.il>
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

from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.matrix.matrix_generic_dense import Matrix_generic_dense
from sage.graphs.graph import Graph
from sage.rings.all import ZZ, QQ, RR
from sage.rings.infinity import infinity
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.coxeter_type import CoxeterType


class CoxeterMatrix(CoxeterType, metaclass=ClasscallMetaclass):
    r"""
    A Coxeter matrix.

    A Coxeter matrix `M = (m_{ij})_{i,j \in I}` is a matrix encoding
    a Coxeter system `(W, S)`, where the relations are given by
    `(s_i s_j)^{m_{ij}}`. Thus `M` is symmetric and has entries
    in `\{1, 2, 3, \ldots, \infty\}` with `m_{ij} = 1` if and only
    if `i = j`.

    We represent `m_{ij} = \infty` by any number `m_{ij} \leq -1`. In
    particular, we can construct a bilinear form `B = (b_{ij})_{i,j \in I}`
    from `M` by

    .. MATH::

        b_{ij} = \begin{cases}
        m_{ij} & m_{ij} < 0\ (\text{i.e., } m_{ij} = \infty), \\
        -\cos\left( \frac{\pi}{m_{ij}} \right) & \text{otherwise}.
        \end{cases}

    EXAMPLES::

        sage: CoxeterMatrix(['A', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
        sage: CoxeterMatrix(['B', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: CoxeterMatrix(['C', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: CoxeterMatrix(['D', 4])
        [1 3 2 2]
        [3 1 3 3]
        [2 3 1 2]
        [2 3 2 1]

        sage: CoxeterMatrix(['E', 6])
        [1 2 3 2 2 2]
        [2 1 2 3 2 2]
        [3 2 1 3 2 2]
        [2 3 3 1 3 2]
        [2 2 2 3 1 3]
        [2 2 2 2 3 1]

        sage: CoxeterMatrix(['F', 4])
        [1 3 2 2]
        [3 1 4 2]
        [2 4 1 3]
        [2 2 3 1]

        sage: CoxeterMatrix(['G', 2])
        [1 6]
        [6 1]

    By default, entries representing `\infty` are given by `-1`
    in the Coxeter matrix::

        sage: G = Graph([(0,1,None), (1,2,4), (0,2,oo)])
        sage: CoxeterMatrix(G)
        [ 1  3 -1]
        [ 3  1  4]
        [-1  4  1]

    It is possible to give a number `\leq -1` to represent an infinite label::

        sage: CoxeterMatrix([[1,-1],[-1,1]])
        [ 1 -1]
        [-1  1]
        sage: CoxeterMatrix([[1,-3/2],[-3/2,1]])
        [   1 -3/2]
        [-3/2    1]
    """
    @staticmethod
    def __classcall_private__(cls, data=None, index_set=None, coxeter_type=None,
                              cartan_type=None, coxeter_type_check=True):
        r"""
        A Coxeter matrix can we created via a graph, a Coxeter type, or
        a matrix.

        .. NOTE::

            To disable the Coxeter type check, use the optional argument
            ``coxeter_type_check = False``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',1,1],['a','b'])
            sage: C2 = CoxeterMatrix([[1, -1], [-1, 1]])
            sage: C3 = CoxeterMatrix(matrix([[1, -1], [-1, 1]]), [0, 1])
            sage: C == C2 and C == C3
            True

        Check with `\infty` because of the hack of using `-1` to represent
        `\infty` in the Coxeter matrix::

            sage: G = Graph([(0, 1, 3), (1, 2, oo)])
            sage: W1 = CoxeterMatrix([[1, 3, 2], [3, 1, -1], [2, -1, 1]])
            sage: W2 = CoxeterMatrix(G)
            sage: W1 == W2
            True
            sage: CoxeterMatrix(W1.coxeter_graph()) == W1
            True

        The base ring of the matrix depends on the entries given::

            sage: CoxeterMatrix([[1,-1],[-1,1]])._matrix.base_ring()
            Integer Ring
            sage: CoxeterMatrix([[1,-3/2],[-3/2,1]])._matrix.base_ring()
            Rational Field
            sage: CoxeterMatrix([[1,-1.5],[-1.5,1]])._matrix.base_ring()
            Real Field with 53 bits of precision
        """
        if not data:
            if coxeter_type:
                data = CoxeterType(coxeter_type)
            elif cartan_type:
                data = CoxeterType(CartanType(cartan_type))

        # Special cases with no arguments passed
        if not data:
            data = []
            n = 0
            index_set = tuple()
            coxeter_type = None
            base_ring = ZZ
            mat = typecall(cls, MatrixSpace(base_ring, n, sparse=False), data, coxeter_type, index_set)
            mat._subdivisions = None

            return mat

        if isinstance(data, CoxeterMatrix):  # Initiate from itself
            return data

        # Initiate from a graph:
        # TODO: Check if a CoxeterDiagram once implemented
        if isinstance(data, Graph):
            return cls._from_graph(data, coxeter_type_check)

        # Get the Coxeter type
        coxeter_type = None
        from sage.combinat.root_system.cartan_type import CartanType_abstract
        if isinstance(data, CartanType_abstract):
            coxeter_type = data.coxeter_type()
        else:
            try:
                coxeter_type = CoxeterType(data)
            except (TypeError, ValueError, NotImplementedError):
                pass

        # Initiate from a Coxeter type
        if coxeter_type:
            return cls._from_coxetertype(coxeter_type)

        # TODO:: remove when oo is possible in matrices.
        n = len(data[0])
        data = [x if x != infinity else -1 for r in data for x in r]
        data = matrix(n, n, data)
        # until here

        # Get the index set
        if index_set:
            index_set = tuple(index_set)
        else:
            index_set = tuple(range(1, n + 1))
        if len(set(index_set)) != n:
                raise ValueError("the given index set is not valid")

        return cls._from_matrix(data, coxeter_type, index_set, coxeter_type_check)

    def __init__(self, parent, data, coxeter_type, index_set):
        """
        Initialize ``self``.

        TESTS::

            sage: C = CoxeterMatrix(['A', 2, 1])
            sage: TestSuite(C).run(skip=["_test_category", "_test_change_ring"])
        """
        self._matrix = Matrix_generic_dense(parent, data, False, True)
        self._matrix.set_immutable()

        if self._matrix.base_ring() not in [ZZ, QQ]:
            self._is_cyclotomic = False
        else:
            self._is_cyclotomic = True
        self._coxeter_type = coxeter_type

        if self._coxeter_type is not None:
            if self._coxeter_type.is_finite():
                self._is_finite = True
                self._is_affine = False
            elif self._coxeter_type.is_affine():
                self._is_finite = False
                self._is_affine = True
            else:
                self._is_finite = False
                self._is_affine = False
        else:
            self._is_finite = False
            self._is_affine = False

        self._index_set = index_set
        self._rank = self._matrix.nrows()

        self._dict = {(self._index_set[i], self._index_set[j]): self._matrix[i, j]
                      for i in range(self._rank) for j in range(self._rank)}

        for i, key in enumerate(self._index_set):
            self._dict[key] = {key2: self._matrix[i, j]
                               for j, key2 in enumerate(self._index_set)}

    @classmethod
    def _from_matrix(cls, data, coxeter_type, index_set, coxeter_type_check):
        """
        Initiate the Coxeter matrix from a matrix.

        TESTS::

            sage: CM = CoxeterMatrix([[1,2],[2,1]]); CM
            [1 2]
            [2 1]
            sage: CM = CoxeterMatrix([[1,-1],[-1,1]]); CM
            [ 1 -1]
            [-1  1]
            sage: CM = CoxeterMatrix([[1,-1.5],[-1.5,1]]); CM
            [ 1.00000000000000 -1.50000000000000]
            [-1.50000000000000  1.00000000000000]
            sage: CM = CoxeterMatrix([[1,-3/2],[-3/2,1]]); CM
            [   1 -3/2]
            [-3/2    1]
            sage: CM = CoxeterMatrix([[1,-3/2,5],[-3/2,1,-1],[5,-1,1]]); CM
            [   1 -3/2    5]
            [-3/2    1   -1]
            [   5   -1    1]
            sage: CM = CoxeterMatrix([[1,-3/2,5],[-3/2,1,oo],[5,oo,1]]); CM
            [   1 -3/2    5]
            [-3/2    1   -1]
            [   5   -1    1]
        """
        # Check that the data is valid
        check_coxeter_matrix(data)

        M = matrix(data)
        n = M.ncols()

        base_ring = M.base_ring()

        if not coxeter_type:
            if n == 1:
                coxeter_type = CoxeterType(['A', 1])
            elif coxeter_type_check:
                coxeter_type = recognize_coxeter_type_from_matrix(M, index_set)
            else:
                coxeter_type = None

        raw_data = M.list()

        mat = typecall(cls, MatrixSpace(base_ring, n, sparse=False), raw_data,
                       coxeter_type, index_set)
        mat._subdivisions = M._subdivisions

        return mat

    @classmethod
    def _from_graph(cls, graph, coxeter_type_check):
        """
        Initiate the Coxeter matrix from a graph.

        TESTS::

            sage: CoxeterMatrix(CoxeterMatrix(['A',4,1]).coxeter_graph())
            [1 3 2 2 3]
            [3 1 3 2 2]
            [2 3 1 3 2]
            [2 2 3 1 3]
            [3 2 2 3 1]
            sage: CoxeterMatrix(CoxeterMatrix(['B',4,1]).coxeter_graph())
            [1 2 3 2 2]
            [2 1 3 2 2]
            [3 3 1 3 2]
            [2 2 3 1 4]
            [2 2 2 4 1]
            sage: CoxeterMatrix(CoxeterMatrix(['F',4]).coxeter_graph())
            [1 3 2 2]
            [3 1 4 2]
            [2 4 1 3]
            [2 2 3 1]

            sage: G = Graph()
            sage: G.add_edge([0,1,oo])
            sage: CoxeterMatrix(G)
            [ 1 -1]
            [-1  1]
            sage: H = Graph()
            sage: H.add_edge([0,1,-1.5])
            sage: CoxeterMatrix(H)
            [ 1.00000000000000 -1.50000000000000]
            [-1.50000000000000  1.00000000000000]
        """
        verts = sorted(graph.vertices())
        index_set = tuple(verts)
        n = len(index_set)

        # Setup the basis matrix as all 2 except 1 on the diagonal
        data = []
        for i in range(n):
            data += [[]]
            for j in range(n):
                if i == j:
                    data[-1] += [ZZ.one()]
                else:
                    data[-1] += [2]

        for e in graph.edges():
            label = e[2]
            if label is None:
                label = 3
            elif label == infinity:
                label = -1
            elif label not in ZZ and label > -1:
                raise ValueError("invalid Coxeter graph label")
            elif label == 0 or label == 1:
                raise ValueError("invalid Coxeter graph label")
            i = verts.index(e[0])
            j = verts.index(e[1])
            data[j][i] = data[i][j] = label

        return cls._from_matrix(data, None, index_set, coxeter_type_check)

    @classmethod
    def _from_coxetertype(cls, coxeter_type):
        """
        Initiate the Coxeter matrix from a Coxeter type.

        TESTS::

            sage: CoxeterMatrix(['A',4]).coxeter_type()
            Coxeter type of ['A', 4]
            sage: CoxeterMatrix(['A',4,1]).coxeter_type()
            Coxeter type of ['A', 4, 1]
            sage: CoxeterMatrix(['D',4,1]).coxeter_type()
            Coxeter type of ['D', 4, 1]
        """
        index_set = coxeter_type.index_set()
        n = len(index_set)
        reverse = {index_set[i]: i for i in range(n)}
        data = [[1 if i == j else 2 for j in range(n)] for i in range(n)]
        for (i, j, l) in coxeter_type.coxeter_graph().edge_iterator():
            if l == infinity:
                l = -1
            data[reverse[i]][reverse[j]] = l
            data[reverse[j]][reverse[i]] = l

        return cls._from_matrix(data, coxeter_type, index_set, False)

    @classmethod
    def samples(self, finite=None, affine=None, crystallographic=None, higher_rank=None):
        """
        Return a sample of the available Coxeter types.

        INPUT:

        - ``finite`` -- (default: ``None``) a boolean or ``None``

        - ``affine`` -- (default: ``None``) a boolean or ``None``

        - ``crystallographic`` -- (default: ``None``) a boolean or ``None``

        - ``higher_rank`` -- (default: ``None``) a boolean or ``None``

        The sample contains all the exceptional finite and affine
        Coxeter types, as well as typical representatives of the
        infinite families.

        Here the ``higher_rank`` term denotes non-finite, non-affine,
        Coxeter groups (including hyperbolic types).

        .. TODO:: Implement the hyperbolic and compact hyperbolic in the samples.

        EXAMPLES::

            sage: [CM.coxeter_type() for CM in CoxeterMatrix.samples()]
            [
            Coxeter type of ['A', 1], Coxeter type of ['A', 5],
            <BLANKLINE>
            Coxeter type of ['B', 5], Coxeter type of ['D', 4],
            <BLANKLINE>
            Coxeter type of ['D', 5], Coxeter type of ['E', 6],
            <BLANKLINE>
            Coxeter type of ['E', 7], Coxeter type of ['E', 8],
            <BLANKLINE>
            Coxeter type of ['F', 4], Coxeter type of ['H', 3],
            <BLANKLINE>
            Coxeter type of ['H', 4], Coxeter type of ['I', 10],
            <BLANKLINE>
            Coxeter type of ['A', 2, 1], Coxeter type of ['B', 5, 1],
            <BLANKLINE>
            Coxeter type of ['C', 5, 1], Coxeter type of ['D', 5, 1],
            <BLANKLINE>
            Coxeter type of ['E', 6, 1], Coxeter type of ['E', 7, 1],
            <BLANKLINE>
            Coxeter type of ['E', 8, 1], Coxeter type of ['F', 4, 1],
            <BLANKLINE>
                                                                      [ 1 -1 -1]
                                                                      [-1  1 -1]
            Coxeter type of ['G', 2, 1], Coxeter type of ['A', 1, 1], [-1 -1  1],
            <BLANKLINE>
                     [ 1 -2  3  2]
            [1 2 3]  [-2  1  2  3]
            [2 1 7]  [ 3  2  1 -8]
            [3 7 1], [ 2  3 -8  1]
            ]

        The finite, affine and crystallographic options allow
        respectively for restricting to (non) finite, (non) affine,
        and (non) crystallographic Cartan types::

            sage: [CM.coxeter_type() for CM in CoxeterMatrix.samples(finite=True)]
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
             Coxeter type of ['B', 5], Coxeter type of ['D', 4],
             Coxeter type of ['D', 5], Coxeter type of ['E', 6],
             Coxeter type of ['E', 7], Coxeter type of ['E', 8],
             Coxeter type of ['F', 4], Coxeter type of ['H', 3],
             Coxeter type of ['H', 4], Coxeter type of ['I', 10]]

            sage: [CM.coxeter_type() for CM in CoxeterMatrix.samples(affine=True)]
            [Coxeter type of ['A', 2, 1], Coxeter type of ['B', 5, 1],
             Coxeter type of ['C', 5, 1], Coxeter type of ['D', 5, 1],
             Coxeter type of ['E', 6, 1], Coxeter type of ['E', 7, 1],
             Coxeter type of ['E', 8, 1], Coxeter type of ['F', 4, 1],
             Coxeter type of ['G', 2, 1], Coxeter type of ['A', 1, 1]]

            sage: [CM.coxeter_type() for CM in CoxeterMatrix.samples(crystallographic=True)]
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
             Coxeter type of ['B', 5], Coxeter type of ['D', 4],
             Coxeter type of ['D', 5], Coxeter type of ['E', 6],
             Coxeter type of ['E', 7], Coxeter type of ['E', 8],
             Coxeter type of ['F', 4], Coxeter type of ['A', 2, 1],
             Coxeter type of ['B', 5, 1], Coxeter type of ['C', 5, 1],
             Coxeter type of ['D', 5, 1], Coxeter type of ['E', 6, 1],
             Coxeter type of ['E', 7, 1], Coxeter type of ['E', 8, 1],
             Coxeter type of ['F', 4, 1], Coxeter type of ['G', 2, 1]]

            sage: CoxeterMatrix.samples(crystallographic=False)
            [
                     [1 3 2 2]
            [1 3 2]  [3 1 3 2]                    [ 1 -1 -1]  [1 2 3]
            [3 1 5]  [2 3 1 5]  [ 1 10]  [ 1 -1]  [-1  1 -1]  [2 1 7]
            [2 5 1], [2 2 5 1], [10  1], [-1  1], [-1 -1  1], [3 7 1],
            <BLANKLINE>
            [ 1 -2  3  2]
            [-2  1  2  3]
            [ 3  2  1 -8]
            [ 2  3 -8  1]
            ]

        .. TODO:: add some reducible Coxeter types (suggestions?)

        TESTS::

            sage: for ct in CoxeterMatrix.samples(): TestSuite(ct).run()
        """
        result = self._samples()
        if crystallographic is not None:
            result = [t for t in result if t.is_crystallographic() == crystallographic]
        if finite is not None:
            result = [t for t in result if t.is_finite() == finite]
        if affine is not None:
            result = [t for t in result if t.is_affine() == affine]
        if higher_rank is not None:
            result = [t for t in result if not t.is_affine() and not t.is_finite()]
        return result

    @cached_method
    def _samples(self):
        """
        Return a sample of all implemented Coxeter types.

        .. NOTE:: This is intended to be used through :meth:`samples`.

        EXAMPLES::

            sage: [CM.coxeter_type() for CM in CoxeterMatrix._samples()]
            [
            Coxeter type of ['A', 1], Coxeter type of ['A', 5],
            <BLANKLINE>
            Coxeter type of ['B', 5], Coxeter type of ['D', 4],
            <BLANKLINE>
            Coxeter type of ['D', 5], Coxeter type of ['E', 6],
            <BLANKLINE>
            Coxeter type of ['E', 7], Coxeter type of ['E', 8],
            <BLANKLINE>
            Coxeter type of ['F', 4], Coxeter type of ['H', 3],
            <BLANKLINE>
            Coxeter type of ['H', 4], Coxeter type of ['I', 10],
            <BLANKLINE>
            Coxeter type of ['A', 2, 1], Coxeter type of ['B', 5, 1],
            <BLANKLINE>
            Coxeter type of ['C', 5, 1], Coxeter type of ['D', 5, 1],
            <BLANKLINE>
            Coxeter type of ['E', 6, 1], Coxeter type of ['E', 7, 1],
            <BLANKLINE>
            Coxeter type of ['E', 8, 1], Coxeter type of ['F', 4, 1],
            <BLANKLINE>
                                                                      [ 1 -1 -1]
                                                                      [-1  1 -1]
            Coxeter type of ['G', 2, 1], Coxeter type of ['A', 1, 1], [-1 -1  1],
            <BLANKLINE>
                     [ 1 -2  3  2]
            [1 2 3]  [-2  1  2  3]
            [2 1 7]  [ 3  2  1 -8]
            [3 7 1], [ 2  3 -8  1]
            ]
        """
        finite = [CoxeterMatrix(t) for t in [['A', 1], ['A', 5], ['B', 5],
                                             ['D', 4], ['D', 5], ['E', 6], ['E', 7],
                                             ['E', 8], ['F', 4], ['H', 3], ['H', 4],
                                             ['I', 10]]]

        affine = [CoxeterMatrix(t) for t in [['A', 2, 1], ['B', 5, 1],
                                             ['C', 5, 1], ['D', 5, 1], ['E', 6, 1],
                                             ['E', 7, 1], ['E', 8, 1], ['F', 4, 1],
                                             ['G', 2, 1], ['A', 1, 1]]]

        higher_matrices = [[[1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
                           [[1, 2, 3], [2, 1, 7], [3, 7, 1]],
                           [[1, -2, 3, 2], [-2, 1, 2, 3], [3, 2, 1, -8], [2, 3, -8, 1]]]

        higher = [CoxeterMatrix(m) for m in higher_matrices]

        return finite + affine + higher

    def relabel(self, relabelling):
        """
        Return a relabelled copy of this Coxeter matrix.

        INPUT:

        - ``relabelling`` -- a function (or dictionary)

        OUTPUT:

        an isomorphic Coxeter type obtained by relabelling the nodes of
        the Coxeter graph. Namely, the node with label ``i`` is
        relabelled ``f(i)`` (or, by ``f[i]`` if ``f`` is a dictionary).

        EXAMPLES::

            sage: CoxeterMatrix(['F',4]).relabel({ 1:2, 2:3, 3:4, 4:1})
            [1 4 2 3]
            [4 1 3 2]
            [2 3 1 2]
            [3 2 2 1]
            sage: CoxeterMatrix(['F',4]).relabel(lambda x: x+1 if x<4 else 1)
            [1 4 2 3]
            [4 1 3 2]
            [2 3 1 2]
            [3 2 2 1]
        """
        if isinstance(relabelling, dict):
            data = [[self[relabelling[i]][relabelling[j]]
                     for j in self.index_set()] for i in self.index_set()]
        else:
            data = [[self[relabelling(i)][relabelling(j)]
                     for j in self.index_set()] for i in self.index_set()]

        return CoxeterMatrix(data)

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: C = CoxeterMatrix(['A',4])
            sage: M = loads(dumps(C))
            sage: M._index_set
            (1, 2, 3, 4)
        """
        if self._coxeter_type:
            return (CoxeterMatrix, (self._coxeter_type,))
        return (CoxeterMatrix, (self._matrix, self._index_set))

    def _repr_(self):
        """
        String representation of the Coxeter matrix.

        EXAMPLES::

            sage: CM = CoxeterMatrix(['A',3]); CM
            [1 3 2]
            [3 1 3]
            [2 3 1]
            sage: CM = CoxeterMatrix([[1,-3/2],[-3/2,1]]); CM
            [   1 -3/2]
            [-3/2    1]
        """
        return repr(self._matrix)

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: CM = CoxeterMatrix(['A',3])
            sage: CM._repr_option('ascii_art')
            True
        """
        if key == 'ascii_art' or key == 'element_ascii_art':
            return self._matrix.nrows() > 1
        return super(CoxeterMatrix, self)._repr_option(key)

    def _latex_(self):
        r"""
        Latex representation of the Coxeter matrix.

        EXAMPLES::

            sage: CM = CoxeterMatrix(['A',3])
            sage: latex(CM)
            \left(\begin{array}{rrr}
            1 & 3 & 2 \\
            3 & 1 & 3 \\
            2 & 3 & 1
            \end{array}\right)
        """
        return self._matrix._latex_()

    def __iter__(self):
        """
        Return an iterator for the rows of the Coxeter matrix.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,8],[8,1]])
            sage: next(CM.__iter__())
            (1, 8)
        """
        return iter(self._matrix)

    def __getitem__(self, key):
        """
        Return a dictionary of labels adjacent to a node or
        the label of an edge in the Coxeter graph.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,-2],[-2,1]])
            sage: CM = CoxeterMatrix([[1,-2],[-2,1]], ['a','b'])
            sage: CM['a']
            {'a': 1, 'b': -2}
            sage: CM['b']
            {'a': -2, 'b': 1}
            sage: CM['a','b']
            -2
            sage: CM['a','a']
            1
        """
        return self._dict[key]

    def __hash__(self):
        r"""
        Return hash of the Coxeter matrix.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1, -2], [-2, 1]], ['a', 'b'])
            sage: CM.__hash__()
            -337812865737895661  # 64-bit
            153276691            # 32-bit
            sage: CM = CoxeterMatrix([[1, -3], [-3, 1]], ['1', '2'])
            sage: CM.__hash__()
            -506719298606843492  # 64-bit
            -1917568612          # 32-bit
        """
        return hash(self._matrix)

    def __eq__(self, other):
        r"""
        Return if ``self`` and ``other`` are equal.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,-2],[-2,1]],['a','b'])
            sage: CM2 = CoxeterMatrix([[1,-2],[-2,1]],['1','2'])
            sage: CM == CM2
            True
            sage: CM == matrix(CM)
            False
            sage: CM3 = CoxeterMatrix([[1,-3],[-3,1]],['1','2'])
            sage: CM == CM3
            False
        """
        return isinstance(other, CoxeterMatrix) and self._matrix == other._matrix

    def __ne__(self, other):
        """
        Return if ``self`` and ``other`` are not equal.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,-2],[-2,1]],['a','b'])
            sage: CM2 = CoxeterMatrix([[1,-2],[-2,1]],['1','2'])
            sage: CM != CM2
            False
            sage: matrix(CM) != CM
            True
            sage: CM3 = CoxeterMatrix([[1,-3],[-3,1]],['1','2'])
            sage: CM != CM3
            True
        """
        return not (self == other)

    def _matrix_(self, R=None):
        """
        Return ``self`` as a matrix over the ring ``R``.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,-3],[-3,1]])
            sage: matrix(CM)
            [ 1 -3]
            [-3  1]
            sage: matrix(RR, CM)
            [ 1.00000000000000 -3.00000000000000]
            [-3.00000000000000  1.00000000000000]
        """
        if R is not None:
            return self._matrix.change_ring(R)
        else:
            return self._matrix

    ##########################################################################
    # Coxeter type methods

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',1,1])
            sage: C.index_set()
            (0, 1)
            sage: C = CoxeterMatrix(['E',6])
            sage: C.index_set()
            (1, 2, 3, 4, 5, 6)
        """
        return self._index_set

    def coxeter_type(self):
        """
        Return the Coxeter type of ``self`` or ``self`` if unknown.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',4,1])
            sage: C.coxeter_type()
            Coxeter type of ['A', 4, 1]

        If the Coxeter type is unknown::

            sage: C = CoxeterMatrix([[1,3,4], [3,1,-1], [4,-1,1]])
            sage: C.coxeter_type()
            [ 1  3  4]
            [ 3  1 -1]
            [ 4 -1  1]
        """
        if self._coxeter_type is None:
            return self
        return self._coxeter_type

    def rank(self):
        r"""
        Return the rank of ``self``.

        EXAMPLES::

            sage: CoxeterMatrix(['C',3]).rank()
            3
            sage: CoxeterMatrix(["A2","B2","F4"]).rank()
            8
        """
        return self._rank

    def coxeter_matrix(self):
        r"""
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: CoxeterMatrix(['C',3]).coxeter_matrix()
            [1 3 2]
            [3 1 4]
            [2 4 1]
        """
        return self

    def bilinear_form(self, R=None):
        r"""
        Return the bilinear form of ``self``.

        EXAMPLES::

            sage: CoxeterType(['A', 2, 1]).bilinear_form()
            [   1 -1/2 -1/2]
            [-1/2    1 -1/2]
            [-1/2 -1/2    1]
            sage: CoxeterType(['H', 3]).bilinear_form()
            [                      1                    -1/2                       0]
            [                   -1/2                       1 1/2*E(5)^2 + 1/2*E(5)^3]
            [                      0 1/2*E(5)^2 + 1/2*E(5)^3                       1]
            sage: C = CoxeterMatrix([[1,-1,-1],[-1,1,-1],[-1,-1,1]])
            sage: C.bilinear_form()
            [ 1 -1 -1]
            [-1  1 -1]
            [-1 -1  1]
        """
        return CoxeterType.bilinear_form(self, R=R)

    @cached_method
    def coxeter_graph(self):
        """
        Return the Coxeter graph of ``self``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',3])
            sage: C.coxeter_graph()
            Graph on 3 vertices

            sage: C = CoxeterMatrix([['A',3],['A',1]])
            sage: C.coxeter_graph()
            Graph on 4 vertices
        """
        n = self.rank()
        I = self.index_set()

        def val(x):
            return infinity if x == -1 else x
        G = Graph([(I[i], I[j], val((self._matrix)[i, j]))
                   for i in range(n) for j in range(i)
                   if self._matrix[i, j] not in [1, 2]],
                  format='list_of_edges')
        G.add_vertices(I)
        return G.copy(immutable=True)

    def is_simply_laced(self):
        """
        Return if ``self`` is simply-laced.

        A Coxeter matrix is simply-laced if all non-diagonal entries are
        either 2 or 3.

        EXAMPLES::

            sage: cm = CoxeterMatrix([[1,3,3,3], [3,1,3,3], [3,3,1,3], [3,3,3,1]])
            sage: cm.is_simply_laced()
            True
        """
        # We include 1 in this list to account for the diagonal
        L = [1, 2, 3]
        return all(x in L for row in self for x in row)

    def is_crystallographic(self):
        """
        Return whether ``self`` is crystallographic.

        A Coxeter matrix is crystallographic if all non-diagonal entries
        are either 2, 3, 4, or 6.

        EXAMPLES::

            sage: CoxeterMatrix(['F',4]).is_crystallographic()
            True
            sage: CoxeterMatrix(['H',3]).is_crystallographic()
            False
        """
        # We include 1 in this list to account for the diagonal
        L = [1, 2, 3, 4, 6]
        return all(x in L for row in self for x in row)

    def is_irreducible(self):
        """
        Return whether ``self`` is irreducible.

        A Coxeter matrix is irreducible if the Coxeter graph is connected.

        EXAMPLES::

            sage: CoxeterMatrix([['F',4],['A',1]]).is_irreducible()
            False
            sage: CoxeterMatrix(['H',3]).is_irreducible()
            True
        """
        return self.coxeter_graph().is_connected()

    def is_finite(self):
        """
        Return if ``self`` is a finite type or ``False`` if unknown.

        EXAMPLES::

            sage: M = CoxeterMatrix(['C',4])
            sage: M.is_finite()
            True
            sage: M = CoxeterMatrix(['D',4,1])
            sage: M.is_finite()
            False
            sage: M = CoxeterMatrix([[1, -1], [-1, 1]])
            sage: M.is_finite()
            False
        """
        return self._is_finite

    def is_affine(self):
        """
        Return if ``self`` is an affine type or ``False`` if unknown.

        EXAMPLES::

            sage: M = CoxeterMatrix(['C',4])
            sage: M.is_affine()
            False
            sage: M = CoxeterMatrix(['D',4,1])
            sage: M.is_affine()
            True
            sage: M = CoxeterMatrix([[1, 3],[3,1]])
            sage: M.is_affine()
            False
            sage: M = CoxeterMatrix([[1, -1, 7], [-1, 1, 3], [7, 3, 1]])
            sage: M.is_affine()
            False
        """
        return self._is_affine


#####################################################################
# Type check functions


def recognize_coxeter_type_from_matrix(coxeter_matrix, index_set):
    """
    Return the Coxeter type of ``coxeter_matrix`` if known,
    otherwise return ``None``.

    EXAMPLES:

    Some infinite ones::

        sage: C = CoxeterMatrix([[1,3,2],[3,1,-1],[2,-1,1]])
        sage: C.is_finite()  # indirect doctest
        False
        sage: C = CoxeterMatrix([[1,-1,-1],[-1,1,-1],[-1,-1,1]])
        sage: C.is_finite()  # indirect doctest
        False

    Some finite ones::

        sage: m = matrix(CoxeterMatrix(['D', 4]))
        sage: CoxeterMatrix(m).is_finite()  # indirect doctest
        True
        sage: m = matrix(CoxeterMatrix(['H', 4]))
        sage: CoxeterMatrix(m).is_finite()  # indirect doctest
        True

        sage: CoxeterMatrix(CoxeterType(['A',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 10]
        sage: CoxeterMatrix(CoxeterType(['B',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10]
        sage: CoxeterMatrix(CoxeterType(['C',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10]
        sage: CoxeterMatrix(CoxeterType(['D',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['D', 10]
        sage: CoxeterMatrix(CoxeterType(['E',6]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 6]
        sage: CoxeterMatrix(CoxeterType(['E',7]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 7]
        sage: CoxeterMatrix(CoxeterType(['E',8]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 8]
        sage: CoxeterMatrix(CoxeterType(['F',4]).coxeter_graph()).coxeter_type()
        Coxeter type of ['F', 4]
        sage: CoxeterMatrix(CoxeterType(['G',2]).coxeter_graph()).coxeter_type()
        Coxeter type of ['G', 2]
        sage: CoxeterMatrix(CoxeterType(['H',3]).coxeter_graph()).coxeter_type()
        Coxeter type of ['H', 3]
        sage: CoxeterMatrix(CoxeterType(['H',4]).coxeter_graph()).coxeter_type()
        Coxeter type of ['H', 4]
        sage: CoxeterMatrix(CoxeterType(['I',100]).coxeter_graph()).coxeter_type()
        Coxeter type of ['I', 100]

    Some affine graphs::

        sage: CoxeterMatrix(CoxeterType(['A',1,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 1, 1]
        sage: CoxeterMatrix(CoxeterType(['A',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['B',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['C',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['C', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['D',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['D', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['E',6,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 6, 1]
        sage: CoxeterMatrix(CoxeterType(['E',7,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 7, 1]
        sage: CoxeterMatrix(CoxeterType(['E',8,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 8, 1]
        sage: CoxeterMatrix(CoxeterType(['F',4,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['F', 4, 1]
        sage: CoxeterMatrix(CoxeterType(['G',2,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['G', 2, 1]

    TESTS:

    Check that we detect relabellings::

        sage: M = CoxeterMatrix([[1,2,3],[2,1,6],[3,6,1]], index_set=['a', 'b', 'c'])
        sage: M.coxeter_type()
        Coxeter type of ['G', 2, 1] relabelled by {0: 'a', 1: 'b', 2: 'c'}

        sage: from sage.combinat.root_system.coxeter_matrix import recognize_coxeter_type_from_matrix
        sage: for C in CoxeterMatrix.samples():
        ....:     relabelling_perm = Permutations(C.index_set()).random_element()
        ....:     relabelling_dict = {C.index_set()[i]: relabelling_perm[i] for i in range(C.rank())}
        ....:     relabeled_matrix = C.relabel(relabelling_dict)._matrix
        ....:     recognized_type = recognize_coxeter_type_from_matrix(relabeled_matrix, relabelling_perm)
        ....:     if C.is_finite() or C.is_affine():
        ....:         assert recognized_type == C.coxeter_type()

    We check the rank 2 cases (:trac:`20419`)::

        sage: for i in range(2, 10):
        ....:     M = matrix([[1,i],[i,1]])
        ....:     CoxeterMatrix(M).coxeter_type()
        Coxeter type of A1xA1 relabelled by {1: 2}
        Coxeter type of ['A', 2]
        Coxeter type of ['B', 2]
        Coxeter type of ['I', 5]
        Coxeter type of ['G', 2]
        Coxeter type of ['I', 7]
        Coxeter type of ['I', 8]
        Coxeter type of ['I', 9]
        sage: CoxeterMatrix(matrix([[1,-1],[-1,1]]), index_set=[0,1]).coxeter_type()
        Coxeter type of ['A', 1, 1]

    Check that this works for reducible types with relabellings
    (:trac:`24892`)::

        sage: CM = CoxeterMatrix([[1,2,5],[2,1,2],[5,2,1]]); CM
        [1 2 5]
        [2 1 2]
        [5 2 1]
        sage: CM.coxeter_type()
        Coxeter type of I5 relabelled by {1: 1, 2: 3}xA1 relabelled by {1: 2}
    """
    # First, we build the Coxeter graph of the group without the edge labels
    n = ZZ(coxeter_matrix.nrows())
    G = Graph([index_set,
               [(index_set[i], index_set[j], coxeter_matrix[i, j])
                for i in range(n) for j in range(i, n)
                if coxeter_matrix[i, j] not in [1, 2]]],
              format='vertices_and_edges')

    types = []
    for S in G.connected_components_subgraphs():
        r = S.num_verts()
        # Handle the special cases first
        if r == 1:
            types.append(CoxeterType(['A', 1]).relabel({1: S.vertices()[0]}))
            continue
        if r == 2:  # Type B2, G2, or I_2(p)
            e = S.edge_labels()[0]
            if e == 3:  # Can't be 2 because it is connected
                ct = CoxeterType(['A', 2])
            elif e == 4:
                ct = CoxeterType(['B', 2])
            elif e == 6:
                ct = CoxeterType(['G', 2])
            elif e > 0 and e < float('inf'):  # Remaining non-affine types
                ct = CoxeterType(['I', e])
            else:  # Otherwise it is infinite dihedral group Z_2 \ast Z_2
                ct = CoxeterType(['A', 1, 1])
            if not ct.is_affine():
                types.append(ct.relabel({1: S.vertices()[0],
                                         2: S.vertices()[1]}))
            else:
                types.append(ct.relabel({0: S.vertices()[0],
                                         1: S.vertices()[1]}))
            continue

        test = [['A', r], ['B', r], ['A', r - 1, 1]]
        if r >= 3:
            if r == 3:
                test += [['G', 2, 1], ['H', 3]]
            test.append(['C', r - 1, 1])
        if r >= 4:
            if r == 4:
                test += [['F', 4], ['H', 4]]
            test += [['D', r], ['B', r - 1, 1]]
        if r >= 5:
            if r == 5:
                test.append(['F', 4, 1])
            test.append(['D', r - 1, 1])
        if r == 6:
            test.append(['E', 6])
        elif r == 7:
            test += [['E', 7], ['E', 6, 1]]
        elif r == 8:
            test += [['E', 8], ['E', 7, 1]]
        elif r == 9:
            test.append(['E', 8, 1])

        found = False
        for ct in test:
            ct = CoxeterType(ct)
            T = ct.coxeter_graph()
            iso, match = T.is_isomorphic(S, certificate=True, edge_labels=True)
            if iso:
                types.append(ct.relabel(match))
                found = True
                break
        if not found:
            return None

    return CoxeterType(types)

#####################################################################
# Other functions


def check_coxeter_matrix(m):
    """
    Check if ``m`` represents a generalized Coxeter matrix and raise
    and error if not.

    EXAMPLES::

        sage: from sage.combinat.root_system.coxeter_matrix import check_coxeter_matrix
        sage: m = matrix([[1,3,2],[3,1,-1],[2,-1,1]])
        sage: check_coxeter_matrix(m)

        sage: m = matrix([[1,3],[3,1],[2,-1]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: not a square matrix

        sage: m = matrix([[1,3,2],[3,1,-1],[2,-1,2]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: the matrix diagonal is not all 1

        sage: m = matrix([[1,3,3],[3,1,-1],[2,-1,1]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: the matrix is not symmetric

        sage: m = matrix([[1,3,1/2],[3,1,-1],[1/2,-1,1]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: invalid Coxeter label 1/2

        sage: m = matrix([[1,3,1],[3,1,-1],[1,-1,1]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: invalid Coxeter label 1
    """
    mat = matrix(m)
    if not mat.is_square():
        raise ValueError("not a square matrix")
    for i, row in enumerate(m):
        if mat[i, i] != 1:
            raise ValueError("the matrix diagonal is not all 1")
        for j, val in enumerate(row[i + 1:]):
            if val != m[j + i + 1][i]:
                raise ValueError("the matrix is not symmetric")
            if val not in ZZ:
                if val > -1 and val in RR and val != infinity:
                    raise ValueError("invalid Coxeter label {}".format(val))
            else:
                if val == 1 or val == 0:
                    raise ValueError("invalid Coxeter label {}".format(val))


def coxeter_matrix_as_function(t):
    """
    Return the Coxeter matrix, as a function.

    INPUT:

    - ``t`` -- a Cartan type

    EXAMPLES::

        sage: from sage.combinat.root_system.coxeter_matrix import coxeter_matrix_as_function
        sage: f = coxeter_matrix_as_function(['A',4])
        sage: matrix([[f(i,j) for j in range(1,5)] for i in range(1,5)])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
    """
    m = CartanType(t).coxeter_matrix()
    return lambda i, j: m[i, j]
