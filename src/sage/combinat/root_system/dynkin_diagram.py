# sage.doctest: needs sage.combinat sage.graphs
"""
Dynkin diagrams

AUTHORS:

- Travis Scrimshaw (2012-04-22): Nicolas M. Thiery moved Cartan matrix creation
  to here and I cached results for speed.

- Travis Scrimshaw (2013-06-11): Changed inputs of Dynkin diagrams to handle
  other Dynkin diagrams and graphs. Implemented remaining Cartan type methods.

- Christian Stump, Travis Scrimshaw (2013-04-11): Added Cartan matrix as
  possible input for Dynkin diagrams.
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
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
from sage.misc.lazy_import import lazy_import
from sage.structure.element import Matrix
from sage.graphs.digraph import DiGraph
from sage.combinat.root_system.cartan_type import CartanType, CartanType_abstract

lazy_import('sage.combinat.root_system.cartan_matrix', 'CartanMatrix')


def DynkinDiagram(*args, **kwds):
    r"""
    Return the Dynkin diagram corresponding to the input.

    INPUT:

    The input can be one of the following:

    - empty to obtain an empty Dynkin diagram
    - a Cartan type
    - a Cartan matrix
    - a Cartan matrix and an indexing set

    One can also input an indexing set by passing a tuple using the optional
    argument ``index_set``.

    The edge multiplicities are encoded as edge labels. For the corresponding
    Cartan matrices, this uses the convention in Hong and Kang, Kac,
    Fulton and Harris, and crystals. This is the **opposite** convention
    in Bourbaki and Wikipedia's Dynkin diagram (:wikipedia:`Dynkin_diagram`).
    That is for `i \neq j`::

        i <--k-- j <==> a_ij = -k
                   <==> -scalar(coroot[i], root[j]) = k
                   <==> multiple arrows point from the longer root
                        to the shorter one

    For example, in type `C_2`, we have::

        sage: C2 = DynkinDiagram(['C',2]); C2
        O=<=O
        1   2
        C2
        sage: C2.cartan_matrix()
        [ 2 -2]
        [-1  2]

    However Bourbaki would have the Cartan matrix as:

    .. MATH::

        \begin{bmatrix}
        2 & -1 \\
        -2 & 2
        \end{bmatrix}.

    EXAMPLES::

        sage: DynkinDiagram(['A', 4])
        O---O---O---O
        1   2   3   4
        A4

        sage: DynkinDiagram(['A',1],['A',1])
        O
        1
        O
        2
        A1xA1

        sage: R = RootSystem("A2xB2xF4")
        sage: DynkinDiagram(R)
        O---O
        1   2
        O=>=O
        3   4
        O---O=>=O---O
        5   6   7   8
        A2xB2xF4

        sage: R = RootSystem("A2xB2xF4")
        sage: CM = R.cartan_matrix(); CM
        [ 2 -1| 0  0| 0  0  0  0]
        [-1  2| 0  0| 0  0  0  0]
        [-----+-----+-----------]
        [ 0  0| 2 -1| 0  0  0  0]
        [ 0  0|-2  2| 0  0  0  0]
        [-----+-----+-----------]
        [ 0  0| 0  0| 2 -1  0  0]
        [ 0  0| 0  0|-1  2 -1  0]
        [ 0  0| 0  0| 0 -2  2 -1]
        [ 0  0| 0  0| 0  0 -1  2]
        sage: DD = DynkinDiagram(CM); DD
        O---O
        1   2
        O=>=O
        3   4
        O---O=>=O---O
        5   6   7   8
        A2xB2xF4
        sage: DD.cartan_matrix()
        [ 2 -1  0  0  0  0  0  0]
        [-1  2  0  0  0  0  0  0]
        [ 0  0  2 -1  0  0  0  0]
        [ 0  0 -2  2  0  0  0  0]
        [ 0  0  0  0  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -2  2 -1]
        [ 0  0  0  0  0  0 -1  2]

    We can also create Dynkin diagrams from arbitrary Cartan matrices::

        sage: C = CartanMatrix([[2, -3], [-4, 2]])
        sage: DynkinDiagram(C)
        Dynkin diagram of rank 2
        sage: C.index_set()
        (0, 1)
        sage: CI = CartanMatrix([[2, -3], [-4, 2]], [3, 5])
        sage: DI = DynkinDiagram(CI)
        sage: DI.index_set()
        (3, 5)
        sage: CII = CartanMatrix([[2, -3], [-4, 2]])
        sage: DII = DynkinDiagram(CII, ('y', 'x'))
        sage: DII.index_set()
        ('x', 'y')

    .. SEEALSO::

        :func:`CartanType` for a general discussion on Cartan
        types and in particular node labeling conventions.

    TESTS:

    Check that :issue:`15277` is fixed by not having edges from 0s::

        sage: CM = CartanMatrix([[2,-1,0,0],[-3,2,-2,-2],[0,-1,2,-1],[0,-1,-1,2]])
        sage: CM
        [ 2 -1  0  0]
        [-3  2 -2 -2]
        [ 0 -1  2 -1]
        [ 0 -1 -1  2]
        sage: CM.dynkin_diagram().edges(sort=True)
        [(0, 1, 3),
         (1, 0, 1),
         (1, 2, 1),
         (1, 3, 1),
         (2, 1, 2),
         (2, 3, 1),
         (3, 1, 2),
         (3, 2, 1)]
    """
    if len(args) == 0:
        return DynkinDiagram_class()
    mat = args[0]
    if isinstance(mat, Matrix):
        mat = CartanMatrix(*args)
    if isinstance(mat, CartanMatrix):
        if mat.cartan_type() is not mat:
            try:
                return mat.cartan_type().dynkin_diagram()
            except AttributeError:
                ct = CartanType(*args)
                raise ValueError("Dynkin diagram data not yet hardcoded for type %s" % ct)
        if len(args) > 1:
            index_set = tuple(args[1])
        elif "index_set" in kwds:
            index_set = tuple(kwds["index_set"])
        else:
            index_set = mat.index_set()
        D = DynkinDiagram_class(index_set=index_set)
        for (i, j) in mat.nonzero_positions():
            if i != j:
                D.add_edge(index_set[i], index_set[j], -mat[j, i])
        return D
    ct = CartanType(*args)
    try:
        return ct.dynkin_diagram()
    except AttributeError:
        raise ValueError("Dynkin diagram data not yet hardcoded for type %s" % ct)


class DynkinDiagram_class(DiGraph, CartanType_abstract):
    """
    A Dynkin diagram.

    .. SEEALSO::

        :func:`DynkinDiagram()`

    INPUT:

    - ``t`` -- a Cartan type, Cartan matrix, or ``None``

    EXAMPLES::

        sage: DynkinDiagram(['A', 3])
        O---O---O
        1   2   3
        A3
        sage: C = CartanMatrix([[2, -3], [-4, 2]])
        sage: DynkinDiagram(C)
        Dynkin diagram of rank 2
        sage: C.dynkin_diagram().cartan_matrix() == C
        True

    TESTS:

    Check that the correct type is returned when copied::

        sage: d = DynkinDiagram(['A', 3])
        sage: type(copy(d))
        <class 'sage.combinat.root_system.dynkin_diagram.DynkinDiagram_class'>

    We check that :issue:`14655` is fixed::

        sage: cd = copy(d)
        sage: cd.add_vertex(4)
        sage: d.vertices(sort=True) != cd.vertices(sort=True)
        True

    Implementation note: if a Cartan type is given, then the nodes
    are initialized from the index set of this Cartan type.
    """

    def __init__(self, t=None, index_set=None, odd_isotropic_roots=[],
                 **options):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: d = DynkinDiagram(["A", 3])
            sage: TestSuite(d).run()
        """
        if isinstance(t, DiGraph):
            if isinstance(t, DynkinDiagram_class):
                self._cartan_type = t._cartan_type
                self._odd_isotropic_roots = tuple(odd_isotropic_roots)
            else:
                self._cartan_type = None
                self._odd_isotropic_roots = ()
            DiGraph.__init__(self, data=t, **options)
            return

        DiGraph.__init__(self, **options)
        self._cartan_type = t
        self._odd_isotropic_roots = tuple(odd_isotropic_roots)
        if index_set is not None:
            self.add_vertices(index_set)
        elif t is not None:
            self.add_vertices(t.index_set())

    def _repr_(self, compact=False):
        """
        EXAMPLES::

            sage: DynkinDiagram(['G',2])     # indirect doctest
              3
            O=<=O
            1   2
            G2
        """
        ct = self.cartan_type()
        result = ct.ascii_art() + "\n" if hasattr(ct, "ascii_art") else ""

        if ct is None or isinstance(ct, CartanMatrix):
            return result+"Dynkin diagram of rank %s" % self.rank()
        else:
            return result+"%s" % ct._repr_(compact=True)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method.

        Override rich output because :meth:`_repr_` outputs ascii
        art. The proper fix will be in :issue:`18328`.

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: E8 = WeylCharacterRing('E8')
            sage: E8.dynkin_diagram()._rich_repr_(dm)
            OutputAsciiArt container
        """
        OutputAsciiArt = display_manager.types.OutputAsciiArt
        OutputPlainText = display_manager.types.OutputPlainText
        if OutputAsciiArt in display_manager.supported_output():
            return OutputAsciiArt(self._repr_())
        else:
            return OutputPlainText(self._repr_())

    def _latex_(self, scale=0.5):
        r"""
        Return a latex representation of this Dynkin diagram.

        EXAMPLES::

            sage: latex(DynkinDiagram(['A',3,1]))
            \begin{tikzpicture}[scale=0.5]
            \draw (-1,0) node[anchor=east] {$A_{3}^{(1)}$};
            \draw (0 cm,0) -- (4 cm,0);
            \draw (0 cm,0) -- (2.0 cm, 1.2 cm);
            \draw (2.0 cm, 1.2 cm) -- (4 cm, 0);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (2.0 cm, 1.2 cm) circle (.25cm) node[anchor=south east]{$0$};
            \end{tikzpicture}
        """
        if self.cartan_type() is None:
            return "Dynkin diagram of rank {}".format(self.rank())

        from sage.graphs.graph_latex import setup_latex_preamble
        setup_latex_preamble()

        ret = "\\begin{{tikzpicture}}[scale={}]\n".format(scale)
        ret += "\\draw (-1,0) node[anchor=east] {{${}$}};\n".format(self.cartan_type()._latex_())
        ret += self.cartan_type()._latex_dynkin_diagram()
        ret += "\\end{tikzpicture}"
        return ret

    def _matrix_(self):
        """
        Return a regular matrix from ``self``.

        EXAMPLES::

            sage: M = DynkinDiagram(['C',3])._matrix_(); M
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
            sage: type(M)
            <class 'sage.combinat.root_system.cartan_matrix.CartanMatrix'>
        """
        return self.cartan_matrix()._matrix_()

    def add_edge(self, i, j, label=1):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.dynkin_diagram import DynkinDiagram_class
            sage: d = DynkinDiagram_class(CartanType(['A',3]))
            sage: sorted(d.edges(sort=True))
            []
            sage: d.add_edge(2, 3)
            sage: sorted(d.edges(sort=True))
            [(2, 3, 1), (3, 2, 1)]
        """
        DiGraph.add_edge(self, i, j, label)
        if not self.has_edge(j,i):
            self.add_edge(j,i,1)

    def __hash__(self):
        """
        EXAMPLES::

            sage: d = CartanType(['A',3]).dynkin_diagram()
            sage: dv = d.vertices(sort=True)
            sage: hash(d) == hash((d.cartan_type(), tuple(dv), tuple(d.edge_iterator(dv))))
            True
        """
        # Should assert for immutability!

        #return hash(self.cartan_type(), self.vertices(sort=True), tuple(self.edges(sort=True)))
        # FIXME: self.edges() currently tests at some point whether
        # self is a vertex of itself which causes an infinite
        # recursion loop. Current workaround: call self.edge_iterator directly
        verts = self.vertices(sort=True)
        return hash((self.cartan_type(), tuple(verts), tuple(self.edge_iterator(verts))))

    @staticmethod
    def an_instance():
        """
        Return an example of Dynkin diagram.

        EXAMPLES::

            sage: from sage.combinat.root_system.dynkin_diagram import DynkinDiagram_class
            sage: g = DynkinDiagram_class.an_instance()
            sage: g
            Dynkin diagram of rank 3
            sage: g.cartan_matrix()
            [ 2 -1 -1]
            [-2  2 -1]
            [-1 -1  2]
        """
        # hyperbolic Dynkin diagram of Exercise 4.9 p. 57 of Kac Infinite Dimensional Lie Algebras.
        g = DynkinDiagram()
        g.add_vertices([1,2,3])
        g.add_edge(1,2,2)
        g.add_edge(1,3)
        g.add_edge(2,3)
        return g

    ##########################################################################
    # Cartan type methods

    @cached_method
    def index_set(self):
        """
        EXAMPLES::

            sage: DynkinDiagram(['C',3]).index_set()
            (1, 2, 3)
            sage: DynkinDiagram("A2","B2","F4").index_set()
            (1, 2, 3, 4, 5, 6, 7, 8)
        """
        return tuple(self.vertices(sort=True))

    def cartan_type(self):
        """
        EXAMPLES::

            sage: DynkinDiagram("A2","B2","F4").cartan_type()
            A2xB2xF4
        """
        return self._cartan_type

    def rank(self):
        r"""
        Return the index set for this Dynkin diagram.

        EXAMPLES::

            sage: DynkinDiagram(['C',3]).rank()
            3
            sage: DynkinDiagram("A2","B2","F4").rank()
            8
        """
        return self.num_verts()

    def dynkin_diagram(self):
        """
        EXAMPLES::

            sage: DynkinDiagram(['C',3]).dynkin_diagram()
            O---O=<=O
            1   2   3
            C3
        """
        return self

    @cached_method
    def cartan_matrix(self):
        r"""
        Return the Cartan matrix for this Dynkin diagram.

        EXAMPLES::

            sage: DynkinDiagram(['C',3]).cartan_matrix()
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
        """
        return CartanMatrix(self)

    def dual(self):
        r"""
        Return the dual Dynkin diagram, obtained by reversing all edges.

        EXAMPLES::

            sage: D = DynkinDiagram(['C',3])
            sage: D.edges(sort=True)
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]
            sage: D.dual()
            O---O=>=O
            1   2   3
            B3
            sage: D.dual().edges(sort=True)
            [(1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]
            sage: D.dual() == DynkinDiagram(['B',3])
            True

        TESTS::

            sage: D = DynkinDiagram(['A',0]); D
            A0
            sage: D.edges(sort=True)
            []
            sage: D.dual()
            A0
            sage: D.dual().edges(sort=True)
            []
            sage: D = DynkinDiagram(['A',1])
            sage: D.edges(sort=True)
            []
            sage: D.dual()
            O
            1
            A1
            sage: D.dual().edges(sort=True)
            []
        """
        result = DynkinDiagram_class(None, odd_isotropic_roots=self._odd_isotropic_roots)
        result.add_vertices(self.vertices(sort=True))
        for source, target, label in self.edges(sort=False):
            result.add_edge(target, source, label)
        result._cartan_type = self._cartan_type.dual() if self._cartan_type is not None else None
        return result

    def relabel(self, *args, **kwds):
        """
        Return the relabelled Dynkin diagram of ``self``.

        INPUT: see :meth:`~sage.graphs.generic_graph.GenericGraph.relabel`

        There is one difference: the default value for ``inplace`` is
        ``False`` instead of ``True``.

        EXAMPLES::

            sage: D = DynkinDiagram(['C',3])
            sage: D.relabel({1:0, 2:4, 3:1})
            O---O=<=O
            0   4   1
            C3 relabelled by {1: 0, 2: 4, 3: 1}
            sage: D
            O---O=<=O
            1   2   3
            C3

            sage: _ = D.relabel({1:0, 2:4, 3:1}, inplace=True)
            sage: D
            O---O=<=O
            0   4   1
            C3 relabelled by {1: 0, 2: 4, 3: 1}

            sage: D = DynkinDiagram(['A', [1,2]])
            sage: Dp = D.relabel({-1:4, 0:-3, 1:3, 2:2})
            sage: Dp
            O---X---O---O
            4   -3  3   2
            A1|2 relabelled by {-1: 4, 0: -3, 1: 3, 2: 2}
            sage: Dp.odd_isotropic_roots()
            (-3,)

            sage: D = DynkinDiagram(['D', 5])
            sage: G, perm = D.relabel(range(5), return_map=True)
            sage: G
                    O 4
                    |
                    |
            O---O---O---O
            0   1   2   3
            D5 relabelled by {1: 0, 2: 1, 3: 2, 4: 3, 5: 4}
            sage: perm
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4}

            sage: perm = D.relabel(range(5), return_map=True, inplace=True)
            sage: D
                    O 4
                    |
                    |
            O---O---O---O
            0   1   2   3
            D5 relabelled by {1: 0, 2: 1, 3: 2, 4: 3, 5: 4}
            sage: perm
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4}
        """
        return_map = kwds.pop("return_map", False)
        inplace = kwds.pop("inplace", False)
        if inplace:
            G = self
        else:
            # We need to copy self because we want to return the
            # permutation and that works when relabelling in place.
            G = self.copy()

        perm = DiGraph.relabel(G, *args, inplace=True, return_map=True, **kwds)
        new_odds = [perm[i] for i in self._odd_isotropic_roots]
        G._odd_isotropic_roots = tuple(new_odds)
        if self._cartan_type is not None:
            G._cartan_type = self._cartan_type.relabel(perm.__getitem__)
        if return_map:
            if inplace:
                return perm
            else:
                return G, perm
        else:
            return G

    def subtype(self, index_set):
        """
        Return a subtype of ``self`` given by ``index_set``.

        A subtype can be considered the Dynkin diagram induced from
        the Dynkin diagram of ``self`` by ``index_set``.

        EXAMPLES::

            sage: D = DynkinDiagram(['A',6,2]); D
            O=<=O---O=<=O
            0   1   2   3
            BC3~
            sage: D.subtype([1,2,3])
            Dynkin diagram of rank 3
        """
        return self.cartan_matrix().subtype(index_set).dynkin_diagram()

    def is_finite(self):
        """
        Check if ``self`` corresponds to a finite root system.

        EXAMPLES::

            sage: CartanType(['F',4]).dynkin_diagram().is_finite()
            True
            sage: D = DynkinDiagram(CartanMatrix([[2, -4], [-3, 2]]))
            sage: D.is_finite()
            False
        """
        if self._cartan_type is not None:
            return self._cartan_type.is_finite()
        return self.cartan_matrix().is_finite()

    def is_affine(self):
        """
        Check if ``self`` corresponds to an affine root system.

        EXAMPLES::

            sage: CartanType(['F',4]).dynkin_diagram().is_affine()
            False
            sage: D = DynkinDiagram(CartanMatrix([[2, -4], [-3, 2]]))
            sage: D.is_affine()
            False
        """
        if self._cartan_type is not None:
            return self._cartan_type.is_affine()
        return self.cartan_matrix().is_affine()

    def is_irreducible(self):
        """
        Check if ``self`` corresponds to an irreducible root system.

        EXAMPLES::

            sage: CartanType(['F',4]).dynkin_diagram().is_irreducible()
            True
            sage: CM = CartanMatrix([[2,-6],[-4,2]])
            sage: CM.dynkin_diagram().is_irreducible()
            True
            sage: CartanType("A2xB3").dynkin_diagram().is_irreducible()
            False
            sage: CM = CartanMatrix([[2,-6,0],[-4,2,0],[0,0,2]])
            sage: CM.dynkin_diagram().is_irreducible()
            False
        """
        if self._cartan_type is not None:
            return self._cartan_type.is_irreducible()
        return self.connected_components_number() == 1

    def is_crystallographic(self):
        """
        Implement :meth:`CartanType_abstract.is_crystallographic`.

        A Dynkin diagram always corresponds to a crystallographic root system.

        EXAMPLES::

            sage: CartanType(['F',4]).dynkin_diagram().is_crystallographic()
            True

        TESTS::

            sage: CartanType(['G',2]).dynkin_diagram().is_crystallographic()
            True
        """
        return True

    def symmetrizer(self):
        """
        Return the symmetrizer of the corresponding Cartan matrix.

        EXAMPLES::

            sage: d = DynkinDiagram()
            sage: d.add_edge(1,2,3)
            sage: d.add_edge(2,3)
            sage: d.add_edge(3,4,3)
            sage: d.symmetrizer()
            Finite family {1: 9, 2: 3, 3: 3, 4: 1}

        TESTS:

        We check that :issue:`15740` is fixed::

            sage: d = DynkinDiagram()
            sage: d.add_edge(1,2,3)
            sage: d.add_edge(2,3)
            sage: d.add_edge(3,4,3)
            sage: L = d.root_system().root_lattice()
            sage: al = L.simple_roots()
            sage: al[1].associated_coroot()
            alphacheck[1]
            sage: al[1].reflection(al[2])
            alpha[1] + 3*alpha[2]
        """
        return self.cartan_matrix().symmetrizer()

    def odd_isotropic_roots(self):
        """
        Return the odd isotropic roots of ``self``.

        EXAMPLES::

            sage: g = DynkinDiagram(['A',4])
            sage: g.odd_isotropic_roots()
            ()
            sage: g = DynkinDiagram(['A',[4,3]])
            sage: g.odd_isotropic_roots()
            (0,)
        """
        return self._odd_isotropic_roots

    def __getitem__(self, i):
        r"""
        With a tuple (i,j) as argument, returns the scalar product
        `\langle \alpha^\vee_i, \alpha_j\rangle`.

        Otherwise, behaves as the usual ``DiGraph.__getitem__``

        EXAMPLES:

        We use the `C_4` Dynkin diagram as a Cartan matrix::

            sage: g = DynkinDiagram(['C',4])
            sage: matrix([[g[i,j] for j in range(1,5)] for i in range(1,5)])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -2]
            [ 0  0 -1  2]

        The neighbors of a node can still be obtained in the usual way::

            sage: [g[i] for i in range(1,5)]
            [[2], [1, 3], [2, 4], [3]]
        """
        if not isinstance(i, tuple):
            return DiGraph.__getitem__(self,i)
        [i,j] = i
        if i == j:
            if i in self._odd_isotropic_roots:
                return 0
            return 2
        elif self.has_edge(j, i):
            return -self.edge_label(j, i)
        else:
            return 0

    def column(self, j):
        """
        Return the `j`-th column `(a_{i,j})_i` of the
        Cartan matrix corresponding to this Dynkin diagram, as a container
        (or iterator) of tuples `(i, a_{i,j})`.

        EXAMPLES::

            sage: g = DynkinDiagram(["B",4])
            sage: [ (i,a) for (i,a) in g.column(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        val = 2 if j not in self._odd_isotropic_roots else 0
        return [(j,val)] + [(i,-m) for (j1, i, m) in self.outgoing_edges(j)]

    def row(self, i):
        """
        Return the `i`-th row `(a_{i,j})_j` of the
        Cartan matrix corresponding to this Dynkin diagram, as a container
        (or iterator) of tuples `(j, a_{i,j})`.

        EXAMPLES::

            sage: g = DynkinDiagram(["C",4])
            sage: [ (i,a) for (i,a) in g.row(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        val = 2 if i not in self._odd_isotropic_roots else 0
        return [(i,val)] + [(j,-m) for (j, i1, m) in self.incoming_edges(i)]

    @cached_method
    def coxeter_diagram(self):
        r"""
        Construct the Coxeter diagram of ``self``.

        .. SEEALSO:: :meth:`CartanType_abstract.coxeter_diagram`

        EXAMPLES::

            sage: cm = CartanMatrix([[2,-5,0],[-2,2,-1],[0,-1,2]])
            sage: D = cm.dynkin_diagram()
            sage: G = D.coxeter_diagram(); G
            Graph on 3 vertices
            sage: G.edges(sort=True)
            [(0, 1, +Infinity), (1, 2, 3)]

            sage: ct = CartanType([['A',2,2], ['B',3]])
            sage: ct.coxeter_diagram()
            Graph on 5 vertices
            sage: ct.dynkin_diagram().coxeter_diagram() == ct.coxeter_diagram()
            True
        """
        from sage.rings.infinity import infinity
        scalarproducts_to_order = {0: 2, 1: 3, 2: 4, 3: 6}
        from sage.graphs.graph import Graph
        coxeter_diagram = Graph(multiedges=False)
        I = self.index_set()
        coxeter_diagram.add_vertices(I)
        for i in I:
            for j in self.neighbors_out(i):
                # avoid adding the edge twice
                if not coxeter_diagram.has_edge(i,j):
                    val = scalarproducts_to_order.get(self[i,j]*self[j,i], infinity)
                    coxeter_diagram.add_edge(i,j, val)
        return coxeter_diagram.copy(immutable=True)


def precheck(t, letter=None, length=None, affine=None, n_ge=None, n=None):
    """
    EXAMPLES::

        sage: from sage.combinat.root_system.dynkin_diagram import precheck
        sage: ct = CartanType(['A',4])
        sage: precheck(ct, letter='C')
        Traceback (most recent call last):
        ...
        ValueError: t[0] must be = 'C'
        sage: precheck(ct, affine=1)
        Traceback (most recent call last):
        ...
        ValueError: t[2] must be = 1
        sage: precheck(ct, length=3)
        Traceback (most recent call last):
        ...
        ValueError: len(t) must be = 3
        sage: precheck(ct, n=3)
        Traceback (most recent call last):
        ...
        ValueError: t[1] must be = 3
        sage: precheck(ct, n_ge=5)
        Traceback (most recent call last):
        ...
        ValueError: t[1] must be >= 5
    """
    if letter is not None:
        if t[0] != letter:
            raise ValueError("t[0] must be = '%s'" % letter)

    if length is not None:
        if len(t) != length:
            raise ValueError("len(t) must be = %s" % length)

    if affine is not None:
        try:
            if t[2] != affine:
                raise ValueError("t[2] must be = %s" % affine)
        except IndexError:
            raise ValueError("t[2] must be = %s" % affine)

    if n_ge is not None:
        if t[1] < n_ge:
            raise ValueError("t[1] must be >= %s" % n_ge)

    if n is not None:
        if t[1] != n:
            raise ValueError("t[1] must be = %s" % n)
