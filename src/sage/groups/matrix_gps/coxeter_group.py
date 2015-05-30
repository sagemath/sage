"""
Coxeter Groups As Matrix Groups

This implements a general Coxeter group as a matrix group by using the
reflection representation.

AUTHORS:

- Travis Scrimshaw (2013-08-28): Initial version
"""

##############################################################################
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.coxeter_groups import CoxeterGroups

from sage.combinat.root_system.cartan_type import CartanType, CartanType_abstract
from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_generic
from sage.groups.matrix_gps.group_element import MatrixGroupElement_generic
from sage.graphs.graph import Graph
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField


class CoxeterMatrixGroup(FinitelyGeneratedMatrixGroup_generic, UniqueRepresentation):
    r"""
    A Coxeter group represented as a matrix group.

    Let `(W, S)` be a Coxeter system. We construct a vector space `V`
    over `\RR` with a basis of `\{ \alpha_s \}_{s \in S}` and inner product

    .. MATH::

        B(\alpha_s, \alpha_t) = -\cos\left( \frac{\pi}{m_{st}} \right)

    where we have `B(\alpha_s, \alpha_t) = -1` if `m_{st} = \infty`. Next we
    define a representation `\sigma_s : V \to V` by

    .. MATH::

        \sigma_s \lambda = \lambda - 2 B(\alpha_s, \lambda) \alpha_s.

    This representation is faithful so we can represent the Coxeter group `W`
    by the set of matrices `\sigma_s` acting on `V`.

    INPUT:

    - ``data`` -- a Coxeter matrix or graph or a Cartan type
    - ``base_ring`` -- (default: the universal cyclotomic field) the base
      ring which contains all values `\cos(\pi/m_{ij})` where `(m_{ij})_{ij}`
      is the Coxeter matrix
    - ``index_set`` -- (optional) an indexing set for the generators

    For more on creating Coxeter groups, see
    :meth:`~sage.combinat.root_system.coxeter_group.CoxeterGroup`.

    .. TODO::

        Currently the label `\infty` is implemented as `-1` in the Coxeter
        matrix.

    EXAMPLES:

    We can create Coxeter groups from Coxeter matrices::

        sage: W = CoxeterGroup([[1, 6, 3], [6, 1, 10], [3, 10, 1]])
        sage: W
        Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [ 1  6  3]
        [ 6  1 10]
        [ 3 10  1]
        sage: W.gens()
        (
        [                 -1 -E(12)^7 + E(12)^11                   1]
        [                  0                   1                   0]
        [                  0                   0                   1],
        <BLANKLINE>
        [                  1                   0                   0]
        [-E(12)^7 + E(12)^11                  -1     E(20) - E(20)^9]
        [                  0                   0                   1],
        <BLANKLINE>
        [              1               0               0]
        [              0               1               0]
        [              1 E(20) - E(20)^9              -1]
        )
        sage: m = matrix([[1,3,3,3], [3,1,3,2], [3,3,1,2], [3,2,2,1]])
        sage: W = CoxeterGroup(m)
        sage: W.gens()
        (
        [-1  1  1  1]  [ 1  0  0  0]  [ 1  0  0  0]  [ 1  0  0  0]
        [ 0  1  0  0]  [ 1 -1  1  0]  [ 0  1  0  0]  [ 0  1  0  0]
        [ 0  0  1  0]  [ 0  0  1  0]  [ 1  1 -1  0]  [ 0  0  1  0]
        [ 0  0  0  1], [ 0  0  0  1], [ 0  0  0  1], [ 1  0  0 -1]
        )
        sage: a,b,c,d = W.gens()
        sage: (a*b*c)^3
        [ 5  1 -5  7]
        [ 5  0 -4  5]
        [ 4  1 -4  4]
        [ 0  0  0  1]
        sage: (a*b)^3
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        sage: b*d == d*b
        True
        sage: a*c*a == c*a*c
        True

    We can create the matrix representation over different base rings and with
    different index sets. Note that the base ring must contain all
    `2*\cos(\pi/m_{ij})` where `(m_{ij})_{ij}` is the Coxeter matrix::

        sage: W = CoxeterGroup(m, base_ring=RR, index_set=['a','b','c','d'])
        sage: W.base_ring()
        Real Field with 53 bits of precision
        sage: W.index_set()
        ('a', 'b', 'c', 'd')

        sage: CoxeterGroup(m, base_ring=ZZ)
        Coxeter group over Integer Ring with Coxeter matrix:
        [1 3 3 3]
        [3 1 3 2]
        [3 3 1 2]
        [3 2 2 1]
        sage: CoxeterGroup([[1,4],[4,1]], base_ring=QQ)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert sqrt(2) to a rational

    Using the well-known conversion between Coxeter matrices and Coxeter
    graphs, we can input a Coxeter graph. Following the standard convention,
    edges with no label (i.e. labelled by ``None``) are treated as 3::

        sage: G = Graph([(0,3,None), (1,3,15), (2,3,7), (0,1,3)])
        sage: W = CoxeterGroup(G); W
        Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [ 1  3  2  3]
        [ 3  1  2 15]
        [ 2  2  1  7]
        [ 3 15  7  1]
        sage: G2 = W.coxeter_graph()
        sage: CoxeterGroup(G2) is W
        True

    Because there currently is no class for `\ZZ \cup \{ \infty \}`, labels
    of `\infty` are given by `-1` in the Coxeter matrix::

        sage: G = Graph([(0,1,None), (1,2,4), (0,2,oo)])
        sage: W = CoxeterGroup(G)
        sage: W.coxeter_matrix()
        [ 1  3 -1]
        [ 3  1  4]
        [-1  4  1]

    We can also create Coxeter groups from Cartan types using the
    ``implementation`` keyword::

        sage: W = CoxeterGroup(['D',5], implementation="reflection")
        sage: W
        Finite Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [1 3 2 2 2]
        [3 1 3 2 2]
        [2 3 1 3 3]
        [2 2 3 1 2]
        [2 2 3 2 1]
        sage: W = CoxeterGroup(['H',3], implementation="reflection")
        sage: W
        Finite Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]
    """
    @staticmethod
    def __classcall_private__(cls, data, base_ring=None, index_set=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: W1 = CoxeterGroup(['A',2], implementation="reflection", base_ring=UniversalCyclotomicField())
            sage: W2 = CoxeterGroup([[1,3],[3,1]], index_set=(1,2))
            sage: W1 is W2
            True
            sage: G1 = Graph([(1,2)])
            sage: W3 = CoxeterGroup(G1)
            sage: W1 is W3
            True
            sage: G2 = Graph([(1,2,3)])
            sage: W4 = CoxeterGroup(G2)
            sage: W1 is W4
            True

        Check with `\infty` because of the hack of using `-1` to represent
        `\infty` in the Coxeter matrix::

            sage: G = Graph([(0, 1, 3), (1, 2, oo)])
            sage: W1 = CoxeterGroup(matrix([[1, 3, 2], [3,1,-1], [2,-1,1]]))
            sage: W2 = CoxeterGroup(G)
            sage: W1 is W2
            True
            sage: CoxeterGroup(W1.coxeter_graph()) is W1
            True
        """
        if isinstance(data, CartanType_abstract):
            if index_set is None:
                index_set = data.index_set()
            data = data.coxeter_matrix()
        elif isinstance(data, Graph):
            G = data
            n = G.num_verts()

            # Setup the basis matrix as all 2 except 1 on the diagonal
            data = matrix(ZZ, [[2]*n]*n)
            for i in range(n):
                data[i, i] = ZZ.one()

            verts = G.vertices()
            for e in G.edges():
                m = e[2]
                if m is None:
                    m = 3
                elif m == infinity or m == -1:  # FIXME: Hack because there is no ZZ\cup\{\infty\}
                    m = -1
                elif m <= 1:
                    raise ValueError("invalid Coxeter graph label")
                i = verts.index(e[0])
                j = verts.index(e[1])
                data[j, i] = data[i, j] = m

            if index_set is None:
                index_set = G.vertices()
        else:
            try:
                data = matrix(data)
            except (ValueError, TypeError):
                data = CartanType(data).coxeter_matrix()
            if not data.is_symmetric():
                raise ValueError("the Coxeter matrix is not symmetric")
            if any(d != 1 for d in data.diagonal()):
                raise ValueError("the Coxeter matrix diagonal is not all 1")
            if any(val <= 1 and val != -1 for i, row in enumerate(data.rows())
                   for val in row[i+1:]):
                raise ValueError("invalid Coxeter label")

            if index_set is None:
                index_set = range(data.nrows())

        if base_ring is None:
            base_ring = UniversalCyclotomicField()
        data.set_immutable()
        return super(CoxeterMatrixGroup, cls).__classcall__(cls,
                     data, base_ring, tuple(index_set))

    def __init__(self, coxeter_matrix, base_ring, index_set):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup([[1,3,2],[3,1,3],[2,3,1]])
            sage: TestSuite(W).run() # long time
            sage: W = CoxeterGroup([[1,3,2],[3,1,4],[2,4,1]], base_ring=QQbar)
            sage: TestSuite(W).run() # long time
            sage: W = CoxeterGroup([[1,3,2],[3,1,6],[2,6,1]])
            sage: TestSuite(W).run(max_runs=30) # long time
            sage: W = CoxeterGroup([[1,3,2],[3,1,-1],[2,-1,1]])
            sage: TestSuite(W).run(max_runs=30) # long time

        We check that :trac:`16630` is fixed::

            sage: CoxeterGroup(['D',4], base_ring=QQ).category()
            Category of finite coxeter groups
            sage: CoxeterGroup(['H',4], base_ring=QQbar).category()
            Category of finite coxeter groups
            sage: F = CoxeterGroups().Finite()
            sage: all(CoxeterGroup([letter,i]) in F
            ....:     for i in range(2,5) for letter in ['A','B','D'])
            True
            sage: all(CoxeterGroup(['E',i]) in F for i in range(6,9))
            True
            sage: CoxeterGroup(['F',4]).category()
            Category of finite coxeter groups
            sage: CoxeterGroup(['G',2]).category()
            Category of finite coxeter groups
            sage: all(CoxeterGroup(['H',i]) in F for i in range(3,5))
            True
            sage: all(CoxeterGroup(['I',i]) in F for i in range(2,5))
            True
        """
        self._matrix = coxeter_matrix
        self._index_set = index_set
        n = ZZ(coxeter_matrix.nrows())
        # Compute the matrix with entries `2 \cos( \pi / m_{ij} )`.
        MS = MatrixSpace(base_ring, n, sparse=True)
        MC = MS._get_matrix_class()
        # FIXME: Hack because there is no ZZ \cup \{ \infty \}: -1 represents \infty
        if base_ring is UniversalCyclotomicField():
            val = lambda x: base_ring.gen(2*x) + ~base_ring.gen(2*x) if x != -1 else base_ring(2)
        else:
            from sage.functions.trig import cos
            from sage.symbolic.constants import pi
            val = lambda x: base_ring(2*cos(pi / x)) if x != -1 else base_ring(2)
        gens = [MS.one() + MC(MS, entries={(i, j): val(coxeter_matrix[i, j])
                                           for j in range(n)},
                              coerce=True, copy=True)
                for i in range(n)]
        # Compute the matrix with entries `- \cos( \pi / m_{ij} )`.
        # This describes the bilinear form corresponding to this
        # Coxeter system, and might lead us out of our base ring.
        base_field = base_ring.fraction_field()
        MS2 = MatrixSpace(base_field, n, sparse=True)
        MC2 = MS2._get_matrix_class()
        self._bilinear = MC2(MS2, entries={(i, j): val(coxeter_matrix[i, j]) / base_field(-2)
                                           for i in range(n) for j in range(n)
                                           if coxeter_matrix[i, j] != 2},
                             coerce=True, copy=True)
        self._bilinear.set_immutable()
        category = CoxeterGroups()
        # Now we shall see if the group is finite, and, if so, refine
        # the category to ``category.Finite()``. Otherwise the group is
        # infinite and we refine the category to ``category.Infinite()``.
        is_finite = self._finite_recognition()
        if is_finite:
            category = category.Finite()
        else:
            category = category.Infinite()
        FinitelyGeneratedMatrixGroup_generic.__init__(self, n, base_ring,
                                                      gens, category=category)

    def _finite_recognition(self):
        """
        Return ``True`` if and only if the type is finite.

        This is an auxiliary function used during the initialisation.

        EXAMPLES:

        Some infinite ones::

            sage: F = CoxeterGroups().Finite()
            sage: W = CoxeterGroup([[1,3,2],[3,1,-1],[2,-1,1]])
            sage: W in F  # indirect doctest
            False
            sage: W = CoxeterGroup([[1,-1,-1],[-1,1,-1],[-1,-1,1]])
            sage: W in F  # indirect doctest
            False

        Some finite ones::

            sage: CoxeterGroup(['D',4], base_ring=QQ) in F  # indirect doctest
            True
            sage: CoxeterGroup(['H',4]) in F  # indirect doctest
            True
        """
        # First, we build the Coxeter graph of the group, without the
        # edge labels.
        coxeter_matrix = self._matrix
        n = ZZ(coxeter_matrix.nrows())
        G = Graph([(i, j) for i in range(n) for j in range(i)
                   if coxeter_matrix[i, j] not in [1, 2]])
        # Coxeter graphs of finite Coxeter groups are forests
        if not(G.is_forest()):
            return False
        comps = G.connected_components()
        finite = True
        # The group is finite if and only if for every connected
        # component ``comp`` of its Coxeter graph, the submatrix of
        # the Coxeter matrix corresponding to ``comp`` is one of the
        # type-A,B,D,E,F,H,I matrices (up to permutation). So we
        # shall check this condition on every ``comp``.
        for comp in comps:
            l = len(comp)
            if l == 1:
                # Any `1 \times 1` Coxeter matrix gives a finite group.
                continue  # A1
            elif l == 2:
                # A dihedral group is finite iff there is no `\infty`
                # in its Coxeter matrix.
                c0, c1 = comp
                if coxeter_matrix[c0, c1] > 0:
                    continue  # I2
                return False
            elif l == 3:
                # The `3`-node case. The finite groups to check for
                # here are `A_3`, `B_3` and `H_3`.
                c0, c1, c2 = comp
                s = sorted([coxeter_matrix[c0, c1],
                            coxeter_matrix[c0, c2],
                            coxeter_matrix[c1, c2]])
                if s[1] == 3 and s[2] in [3, 4, 5]:
                    continue  # A3, B3, H3
                return False
            elif l == 4:
                # The `4`-node case. The finite groups to check for
                # here are `A_4`, `B_4`, `D_4`, `F_4` and `H_4`.
                c0, c1, c2, c3 = comp
                u = [coxeter_matrix[c0, c1],
                     coxeter_matrix[c0, c2],
                     coxeter_matrix[c0, c3],
                     coxeter_matrix[c1, c2],
                     coxeter_matrix[c1, c3],
                     coxeter_matrix[c2, c3]]
                s = sorted(u)
                # ``s`` is the list of all off-diagonal entries of
                # the ``comp``-submatrix of the Coxeter matrix,
                # sorted in increasing order.
                if s[3:5] == [3, 3]:
                    if s[5] == 3:
                        continue  # A4, D4
                    if s[5] in [4, 5]:
                        u0 = u[0] + u[1] + u[2]
                        u1 = u[0] + u[3] + u[4]
                        u2 = u[1] + u[3] + u[5]
                        u3 = u[2] + u[4] + u[5]
                        ss = sorted([u0, u1, u2, u3])
                        if ss in [[7, 7, 9, 9], [7, 8, 8, 9],
                                  [7, 8, 9, 10]]:
                            continue  # F4, B4, H4
                return False
            else:
                # The case of `l \geq 5` nodes. The finite
                # groups to check for here are `A_l`, `B_l`, `D_l`,
                # and `E_l` (for `l = 6, 7, 8`).

                # Checking that the Coxeter matrix of the subgroup
                # corresponding to the vertices ``comp`` has all its
                # off-diagonal entries equal to 2, 3 or at most once 4
                found_a_4 = False
                for j in range(l):
                    for i in range(j):
                        coxeter_entry = coxeter_matrix[comp[i], comp[j]]
                        if coxeter_entry in [2, 3]:
                            continue
                        if coxeter_entry == 4 and not found_a_4:
                            found_a_4 = True
                            continue
                        return False

                G0 = G.subgraph(comp)
                if found_a_4:
                    # The case when a `4` has been found in the
                    # Coxeter matrix. This needs only to be checked
                    # against `B_l`. We use the observation that
                    # the group is `B_l` if and only if the Coxeter
                    # graph is an `l`-path (i.e., has diameter
                    # `l - 1`) and the `4` corresponds to one of
                    # its two outermost edges.
                    diameter = G0.diameter()
                    if diameter != l - 1:
                        return False

                    ecc = sorted(((u, v) for (v, u) in G0.eccentricity(with_labels=True).items()))
                    left_end = ecc[-1][1]
                    right_end = ecc[-2][1]
                    left_almost_end = G0.neigbors(left_end)[0]
                    right_almost_end = G0.neigbors(right_end)[0]
                    if (coxeter_matrix[left_end, left_almost_end] == 4
                        or coxeter_matrix[right_end, right_almost_end] == 4):
                        continue  # Bl
                    return False

                # Now, all off-diagonal entries of the Coxeter matrix
                # are 2's and 3's. We need to check our group against
                # `A_l`, `D_l` and `E_l`. Knowing that the Coxeter
                # graph is a tree, we can use its vertex
                # eccentricities to check this.
                ecc = sorted(G0.eccentricity())
                if ecc[-1] == l - 1:
                    continue  # Al
                if ecc[-3] == l - 2:
                    continue  # Dl
                if l <= 8 and ecc[-2] == l - 2 and ecc[-5] == l - 3:
                    continue  # El
                return False

        return True

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CoxeterGroup([[1,3,2],[3,1,3],[2,3,1]])
            Finite Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
            [1 3 2]
            [3 1 3]
            [2 3 1]
        """
        rep = "Finite " if self.is_finite() else ""
        rep += "Coxeter group over {} with Coxeter matrix:\n{}".format(self.base_ring(), self._matrix)
        return rep

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup([[1,3],[3,1]])
            sage: W.index_set()
            (0, 1)
            sage: W = CoxeterGroup([[1,3],[3,1]], index_set=['x', 'y'])
            sage: W.index_set()
            ('x', 'y')
            sage: W = CoxeterGroup(['H',3])
            sage: W.index_set()
            (1, 2, 3)
        """
        return self._index_set

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup([[1,3],[3,1]])
            sage: W.coxeter_matrix()
            [1 3]
            [3 1]
            sage: W = CoxeterGroup(['H',3])
            sage: W.coxeter_matrix()
            [1 3 2]
            [3 1 5]
            [2 5 1]
        """
        return self._matrix

    def coxeter_graph(self):
        """
        Return the Coxeter graph of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['H',3], implementation="reflection")
            sage: G = W.coxeter_graph(); G
            Graph on 3 vertices
            sage: G.edges()
            [(1, 2, None), (2, 3, 5)]
            sage: CoxeterGroup(G) is W
            True
            sage: G = Graph([(0, 1, 3), (1, 2, oo)])
            sage: W = CoxeterGroup(G)
            sage: W.coxeter_graph() == G
            True
            sage: CoxeterGroup(W.coxeter_graph()) is W
            True
        """
        G = Graph()
        G.add_vertices(self.index_set())
        for i, row in enumerate(self._matrix.rows()):
            for j, val in enumerate(row[i+1:]):
                if val == 3:
                    G.add_edge(self._index_set[i], self._index_set[i+1+j])
                elif val > 3:
                    G.add_edge(self._index_set[i], self._index_set[i+1+j], val)
                elif val == -1:  # FIXME: Hack because there is no ZZ\cup\{\infty\}
                    G.add_edge(self._index_set[i], self._index_set[i+1+j], infinity)
        return G

    def bilinear_form(self):
        r"""
        Return the bilinear form associated to ``self``.

        Given a Coxeter group `G` with Coxeter matrix `M = (m_{ij})_{ij}`,
        the associated bilinear form `A = (a_{ij})_{ij}` is given by

        .. MATH::

            a_{ij} = -\cos\left( \frac{\pi}{m_{ij}} \right).

        If `A` is positive definite, then `G` is of finite type (and so
        the associated Coxeter group is a finite group). If `A` is
        positive semidefinite, then `G` is affine type.

        EXAMPLES::

            sage: W = CoxeterGroup(['D',4])
            sage: W.bilinear_form()
            [   1 -1/2    0    0]
            [-1/2    1 -1/2 -1/2]
            [   0 -1/2    1    0]
            [   0 -1/2    0    1]
        """
        return self._bilinear

    def is_finite(self):
        """
        Return ``True`` if this group is finite.

        EXAMPLES::

            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2],[3,1,l],[2,l,1]]).is_finite()]
            ....:
            [2, 3, 4, 5]
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2,2],[3,1,l,2],[2,l,1,3],[2,2,3,1]]).is_finite()]
            ....:
            [2, 3, 4]
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2,2,2], [3,1,3,3,2], [2,3,1,2,2],
            ....:                [2,3,2,1,l], [2,2,2,l,1]]).is_finite()]
            ....:
            [2, 3]
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2,2,2], [3,1,2,3,3], [2,2,1,l,2],
            ....:                [2,3,l,1,2], [2,3,2,2,1]]).is_finite()]
            ....:
            [2, 3]
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2,2,2,2], [3,1,l,2,2,2], [2,l,1,3,l,2],
            ....:                [2,2,3,1,2,2], [2,2,l,2,1,3], [2,2,2,2,3,1]]).is_finite()]
            ....:
            [2, 3]
        """
        # Finite Coxeter groups are marked as finite in
        # their ``__init__`` method, so we can just check
        # the category of ``self``.
        return "Finite" in self.category().axioms()

    def order(self):
        """
        Return the order of ``self``.

        If the Coxeter group is finite, this uses an iterator.

        EXAMPLES::

            sage: W = CoxeterGroup([[1,3],[3,1]])
            sage: W.order()
            6
            sage: W = CoxeterGroup([[1,-1],[-1,1]])
            sage: W.order()
            +Infinity
        """
        if self.is_finite():
            return self._cardinality_from_iterator()
        return infinity

    def canonical_representation(self):
        r"""
        Return the canonical faithful representation of ``self``, which
        is ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup([[1,3],[3,1]])
            sage: W.canonical_representation() is W
            True
        """
        return self

    def simple_reflection(self, i):
        """
        Return the simple reflection `s_i`.

        INPUT:

        - ``i`` -- an element from the index set

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation="reflection")
            sage: W.simple_reflection(1)
            [-1  1  0]
            [ 0  1  0]
            [ 0  0  1]
            sage: W.simple_reflection(2)
            [ 1  0  0]
            [ 1 -1  1]
            [ 0  0  1]
            sage: W.simple_reflection(3)
            [ 1  0  0]
            [ 0  1  0]
            [ 0  1 -1]
        """
        if not i in self._index_set:
            raise ValueError("%s is not in the index set %s" % (i, self.index_set()))
        return self.gen(self._index_set.index(i))

    class Element(MatrixGroupElement_generic):
        """
        A Coxeter group element.
        """
        def has_right_descent(self, i):
            r"""
            Return whether ``i`` is a right descent of ``self``.

            A Coxeter system `(W, S)` has a root system defined as
            `\{ w(\alpha_s) \}_{w \in W}` and we define the positive
            (resp. negative) roots `\alpha = \sum_{s \in S} c_s \alpha_s`
            by all `c_s \geq 0` (resp. `c_s \leq 0`). In particular, we note
            that if `\ell(w s) > \ell(w)` then `w(\alpha_s) > 0` and if
            `\ell(ws) < \ell(w)` then `w(\alpha_s) < 0`.
            Thus `i \in I` is a right descent if `w(\alpha_{s_i}) < 0`
            or equivalently if the matrix representing `w` has all entries
            of the `i`-th column being non-positive.

            INPUT:

            - ``i`` -- an element in the index set

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], implementation="reflection")
                sage: a,b,c = W.gens()
                sage: elt = b*a*c
                sage: map(lambda i: elt.has_right_descent(i), [1, 2, 3])
                [True, False, True]
            """
            i = self.parent()._index_set.index(i)
            col = self.matrix().column(i)
            return all(x <= 0 for x in col)

        def canonical_matrix(self):
            r"""
            Return the matrix of ``self`` in the canonical faithful
            representation, which is ``self`` as a matrix.

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], implementation="reflection")
                sage: a,b,c = W.gens()
                sage: elt = a*b*c
                sage: elt.canonical_matrix()
                [ 0  0 -1]
                [ 1  0 -1]
                [ 0  1 -1]
            """
            return self.matrix()
