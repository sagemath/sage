"""
Root system data for reducible Cartan types
"""
# ****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.combinat.root_system.cartan_type import CartanType_abstract, CartanType_simple, CartanType_finite, CartanType_simply_laced, CartanType_crystallographic
from sage.matrix.constructor import block_diagonal_matrix
from sage.sets.family import Family
from . import ambient_space
import sage.combinat.root_system as root_system
from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp_method, richcmp, rich_to_bool


@richcmp_method
class CartanType(SageObject, CartanType_abstract):
    r"""
    A class for reducible Cartan types.

    Reducible root systems are ones that can be factored as direct
    products. Strictly speaking type `D_2` (corresponding to
    orthogonal groups of degree 4) is reducible since it is
    isomorphic to `A_1\times A_1`. However type `D_2` is not built
    using this class for our purposes.

    INPUT:

    - ``types`` -- list of simple Cartan types

    EXAMPLES::

        sage: t1, t2 = [CartanType(x) for x in (['A',1], ['B',2])]
        sage: CartanType([t1, t2])
        A1xB2
        sage: t = CartanType("A2xB2")

    A reducible Cartan type is finite (resp. crystallographic,
    simply laced) if all its components are::

        sage: t.is_finite()
        True
        sage: t.is_crystallographic()
        True
        sage: t.is_simply_laced()
        False

    This is implemented by inserting the appropriate abstract
    super classes (see :meth:`~sage.combinat.root_system.cartan_type.CartanType_abstract._add_abstract_superclass`)::

        sage: t.__class__.mro()
        [<class 'sage.combinat.root_system.type_reducible.CartanType_with_superclass'>, <class 'sage.combinat.root_system.type_reducible.CartanType'>, <class 'sage.structure.sage_object.SageObject'>, <class 'sage.combinat.root_system.cartan_type.CartanType_finite'>, <class 'sage.combinat.root_system.cartan_type.CartanType_crystallographic'>, <class 'sage.combinat.root_system.cartan_type.CartanType_abstract'>, <class 'object'>]

    The index set of the reducible Cartan type is obtained by
    relabelling successively the nodes of the Dynkin diagrams of
    the components by 1,2,...::

        sage: t = CartanType(["A",4], ["BC",5,2], ["C",3])
        sage: t.index_set()
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

        sage: t.dynkin_diagram()                                                        # needs sage.graphs
        O---O---O---O
        1   2   3   4
        O=<=O---O---O---O=<=O
        5   6   7   8   9   10
        O---O=<=O
        11   12   13
        A4xBC5~xC3
    """

    def __init__(self, types):
        """
        Initialize ``self``.

        TESTS:

        Internally, this relabelling is stored as a dictionary::

            sage: t = CartanType(["A",4], ["BC",5,2], ["C",3])
            sage: sorted(t._index_relabelling.items())
            [((0, 1), 1), ((0, 2), 2), ((0, 3), 3), ((0, 4), 4),
             ((1, 0), 5), ((1, 1), 6), ((1, 2), 7), ((1, 3), 8), ((1, 4), 9), ((1, 5), 10),
             ((2, 1), 11), ((2, 2), 12), ((2, 3), 13)]

        Similarly, the attribute ``_shifts`` specifies by how much the
        indices of the bases of the ambient spaces of the components
        are shifted in the ambient space of this Cartan type::

            sage: t = CartanType("A2xB2")
            sage: t._shifts
            [0, 3, 5]
            sage: A = t.root_system().ambient_space(); A
            Ambient space of the Root system of type A2xB2
            sage: A.ambient_spaces()
            [Ambient space of the Root system of type ['A', 2], Ambient space of the Root system of type ['B', 2]]
            sage: x = A.ambient_spaces()[0]([2,1,0]); x
            (2, 1, 0)
            sage: A.inject_weights(0,x)
            (2, 1, 0, 0, 0)
            sage: x = A.ambient_spaces()[1]([1,0]); x
            (1, 0)
            sage: A.inject_weights(1,x)
            (0, 0, 0, 1, 0)

        More tests::

            sage: TestSuite(t).run()
        """
        self._types = types
        self.affine = False
        indices = (None,) + tuple( (i, j)
                                   for i in range(len(types))
                                   for j in types[i].index_set() )
        self._indices = indices
        self._index_relabelling = dict((indices[i], i) for i in range(1, len(indices)))

        self._spaces = [t.root_system().ambient_space() for t in types]
        if all(l is not None for l in self._spaces):
            self._shifts = [sum(l.dimension() for l in self._spaces[:k])
                            for k in range(len(types)+1)]

        self.tools = root_system.type_reducible
        # a direct product of finite Cartan types is again finite;
        # idem for simply laced and crystallographic.
        super_classes = tuple( cls
                               for cls in (CartanType_finite, CartanType_simply_laced, CartanType_crystallographic)
                               if all(isinstance(t, cls) for t in types) )
        self._add_abstract_superclass(super_classes)

    def _repr_(self, compact=True):  # We should make a consistent choice here
        """
        EXAMPLES::

           sage: CartanType("A2","B2")    # indirect doctest
           A2xB2

           sage: CartanType("A2",CartanType("F4~").dual())
           A2xF4~*
        """
        return "x".join(t._repr_(compact=True) for t in self._types)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType("A4","B2","D8"))
            A_{4} \times B_{2} \times D_{8}
        """
        return " \\times ".join(x._latex_() for x in self.component_types())

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: ct0 = CartanType(['A',1],['B',2])
            sage: ct1 = CartanType(['A',2],['B',3])
            sage: hash(ct0) != hash(ct1)
            True
        """
        return hash(repr(self._types))

    def __richcmp__(self, other, op):
        """
        Rich comparison.

        EXAMPLES::

            sage: ct1 = CartanType(['A',1],['B',2])
            sage: ct2 = CartanType(['B',2],['A',1])
            sage: ct3 = CartanType(['A',4])
            sage: ct1 == ct1
            True
            sage: ct1 == ct2
            False
            sage: ct1 == ct3
            False

        TESTS:

        Check that :issue:`20418` is fixed::

            sage: ct = CartanType(["A2", "B2"])
            sage: ct == (1, 2, 1)
            False
        """
        if isinstance(other, CartanType_simple):
            return rich_to_bool(op, 1)
        if not isinstance(other, CartanType):
            return NotImplemented
        return richcmp(self._types, other._types, op)

    def component_types(self):
        """
        A list of Cartan types making up the reducible type.

        EXAMPLES::

            sage: CartanType(['A',2],['B',2]).component_types()
            [['A', 2], ['B', 2]]
        """
        return self._types

    def type(self):
        """
        Return ``"reducible"`` since the type is reducible.

        EXAMPLES::

            sage: CartanType(['A',2],['B',2]).type()
            'reducible'
        """
        return "reducible"

    def rank(self):
        """
        Return the rank of ``self``.

        EXAMPLES::

            sage: CartanType("A2","A1").rank()
            3
        """
        return sum(t.rank() for t in self._types)

    @cached_method
    def index_set(self):
        r"""
        Implement :meth:`CartanType_abstract.index_set`.

        For the moment, the index set is always of the form `\{1, \ldots, n\}`.

        EXAMPLES::

            sage: CartanType("A2","A1").index_set()
            (1, 2, 3)
        """
        return tuple(range(1, self.rank()+1))

    def cartan_matrix(self, subdivide=True):
        """
        Return the Cartan matrix associated with ``self``. By default
        the Cartan matrix is a subdivided block matrix showing the
        reducibility but the subdivision can be suppressed with
        the option ``subdivide = False``.

        EXAMPLES::

            sage: ct = CartanType("A2","B2")
            sage: ct.cartan_matrix()                                                    # needs sage.graphs
            [ 2 -1| 0  0]
            [-1  2| 0  0]
            [-----+-----]
            [ 0  0| 2 -1]
            [ 0  0|-2  2]
            sage: ct.cartan_matrix(subdivide=False)                                     # needs sage.graphs
            [ 2 -1  0  0]
            [-1  2  0  0]
            [ 0  0  2 -1]
            [ 0  0 -2  2]
            sage: ct.index_set() == ct.cartan_matrix().index_set()                      # needs sage.graphs
            True
        """
        from sage.combinat.root_system.cartan_matrix import CartanMatrix
        return CartanMatrix(block_diagonal_matrix([t.cartan_matrix() for t in self._types], subdivide=subdivide),
                            cartan_type=self, index_set=self.index_set())

    def dynkin_diagram(self):
        """
        Return a Dynkin diagram for type reducible.

        EXAMPLES::

            sage: dd = CartanType("A2xB2xF4").dynkin_diagram(); dd                      # needs sage.graphs
            O---O
            1   2
            O=>=O
            3   4
            O---O=>=O---O
            5   6   7   8
            A2xB2xF4
            sage: dd.edges(sort=True)                                                   # needs sage.graphs
            [(1, 2, 1), (2, 1, 1), (3, 4, 2), (4, 3, 1), (5, 6, 1),
             (6, 5, 1), (6, 7, 2), (7, 6, 1), (7, 8, 1), (8, 7, 1)]

            sage: CartanType("F4xA2").dynkin_diagram()                                  # needs sage.graphs
            O---O=>=O---O
            1   2   3   4
            O---O
            5   6
            F4xA2
        """
        from .dynkin_diagram import DynkinDiagram_class
        relabelling = self._index_relabelling
        g = DynkinDiagram_class(self)
        for i in range(len(self._types)):
            for [e1, e2, l] in self._types[i].dynkin_diagram().edges(sort=True):
                g.add_edge(relabelling[i,e1], relabelling[i,e2], label=l)
        return g

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        .. NOTE::

            The arguments ``label`` and ``dual`` is ignored.

        EXAMPLES::

            sage: print(CartanType("A2","B2")._latex_dynkin_diagram())
            {
            \draw (0 cm,0) -- (2 cm,0);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \pgftransformyshift{-3 cm}
            \draw (0 cm,0) -- (0 cm,0);
            \draw (0 cm, 0.1 cm) -- +(2 cm,0);
            \draw (0 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(1.2, 0)}, rotate=0] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            }
        """
        if label is None:
            label = lambda x: x
        types = self.component_types()
        relabelling = self._index_relabelling
        ret = "{\n"
        ret += "\\pgftransformyshift{-3 cm}\n".join(types[i]._latex_dynkin_diagram(
                    lambda x: label(relabelling[i,x]), node, node_dist=node_dist)
                    for i in range(len(types)))
        ret += "}"
        return ret

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of this reducible Cartan type.

        EXAMPLES::

            sage: print(CartanType("F4xA2").ascii_art(label = lambda x: x+2))
            O---O=>=O---O
            3   4   5   6
            O---O
            7   8

            sage: print(CartanType(["BC",5,2], ["A",4]).ascii_art())
            O=<=O---O---O---O=<=O
            1   2   3   4   5   6
            O---O---O---O
            7   8   9   10

            sage: print(CartanType(["A",4], ["BC",5,2], ["C",3]).ascii_art())
            O---O---O---O
            1   2   3   4
            O=<=O---O---O---O=<=O
            5   6   7   8   9   10
            O---O=<=O
            11   12   13
        """
        if label is None:
            label = lambda i: i
        types = self.component_types()
        relabelling = self._index_relabelling
        return "\n".join(types[i].ascii_art(lambda x: label(relabelling[i,x]), node)
                         for i in range(len(types)))

    @cached_method
    def is_finite(self):
        """
        EXAMPLES::

            sage: ct1 = CartanType(['A',2],['B',2])
            sage: ct1.is_finite()
            True
            sage: ct2 = CartanType(['A',2],['B',2,1])
            sage: ct2.is_finite()
            False

        TESTS::

            sage: isinstance(ct1, sage.combinat.root_system.cartan_type.CartanType_finite)
            True
            sage: isinstance(ct2, sage.combinat.root_system.cartan_type.CartanType_finite)
            False
        """
        return all(t.is_finite() for t in self.component_types())

    def is_irreducible(self):
        """
        Report that this Cartan type is not irreducible.

        EXAMPLES::

            sage: ct = CartanType(['A',2],['B',2])
            sage: ct.is_irreducible()
            False
        """
        return False

    def dual(self):
        """
        EXAMPLES::

            sage: CartanType("A2xB2").dual()
            A2xC2
        """
        return CartanType([t.dual() for t in self._types])

    def is_affine(self):
        """
        Report that this reducible Cartan type is not affine.

        EXAMPLES::

            sage: CartanType(['A',2],['B',2]).is_affine()
            False
        """
        return False

    @cached_method
    def coxeter_diagram(self):
        """
        Return the Coxeter diagram for ``self``.

        EXAMPLES::

            sage: cd = CartanType("A2xB2xF4").coxeter_diagram(); cd                     # needs sage.graphs
            Graph on 8 vertices
            sage: cd.edges(sort=True)                                                   # needs sage.graphs
            [(1, 2, 3), (3, 4, 4), (5, 6, 3), (6, 7, 4), (7, 8, 3)]

            sage: CartanType("F4xA2").coxeter_diagram().edges(sort=True)                # needs sage.graphs
            [(1, 2, 3), (2, 3, 4), (3, 4, 3), (5, 6, 3)]

            sage: cd = CartanType("A1xH3").coxeter_diagram(); cd                        # needs sage.graphs
            Graph on 4 vertices
            sage: cd.edges(sort=True)                                                   # needs sage.graphs
            [(2, 3, 3), (3, 4, 5)]
        """
        from sage.graphs.graph import Graph
        relabelling = self._index_relabelling
        g = Graph(multiedges=False)
        g.add_vertices(self.index_set())
        for i,t in enumerate(self._types):
            for [e1, e2, l] in t.coxeter_diagram().edges(sort=True):
                g.add_edge(relabelling[i,e1], relabelling[i,e2], label=l)
        return g


class AmbientSpace(ambient_space.AmbientSpace):
    """
    EXAMPLES::

        sage: RootSystem("A2xB2").ambient_space()
        Ambient space of the Root system of type A2xB2
    """

    def cartan_type(self):
        """
        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().cartan_type()
            A2xB2
        """
        return self.root_system.cartan_type()

    def component_types(self):
        """
        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().component_types()
            [['A', 2], ['B', 2]]
        """
        return self.root_system.cartan_type().component_types()

    def dimension(self):
        """
        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().dimension()
            5
        """
        return sum(v.dimension() for v in self.ambient_spaces())

    def ambient_spaces(self):
        """
        Return a list of the irreducible Cartan types of which the
        given reducible Cartan type is a product.

        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().ambient_spaces()
            [Ambient space of the Root system of type ['A', 2],
             Ambient space of the Root system of type ['B', 2]]
        """
        return [t.root_system().ambient_space() for t in self.component_types()]

    def inject_weights(self, i, v):
        r"""
        Produces the corresponding element of the lattice.

        INPUT:

        - ``i`` -- integer in ``range(self.components)``

        - ``v`` -- a vector in the `i`-th component weight lattice

        EXAMPLES::

            sage: V = RootSystem("A2xB2").ambient_space()
            sage: [V.inject_weights(i,V.ambient_spaces()[i].fundamental_weights()[1]) for i in range(2)]
            [(1, 0, 0, 0, 0), (0, 0, 0, 1, 0)]
            sage: [V.inject_weights(i,V.ambient_spaces()[i].fundamental_weights()[2]) for i in range(2)]
            [(1, 1, 0, 0, 0), (0, 0, 0, 1/2, 1/2)]
        """
        shift = self.root_system.cartan_type()._shifts[i]
        return self._from_dict({shift + k: c for k, c in v})

    @cached_method
    def simple_root(self, i):
        """
        EXAMPLES::

            sage: A = RootSystem("A1xB2").ambient_space()
            sage: A.simple_root(2)
            (0, 0, 1, -1)
            sage: A.simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 1)}
        """
        if i not in self.index_set():
            raise ValueError("{} is not in the index set".format(i))
        (i, j) = self.cartan_type()._indices[i]
        return self.inject_weights(i, self.ambient_spaces()[i].simple_root(j))

    @cached_method
    def simple_coroot(self, i):
        """
        EXAMPLES::

            sage: A = RootSystem("A1xB2").ambient_space()
            sage: A.simple_coroot(2)
            (0, 0, 1, -1)
            sage: A.simple_coroots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 2)}
        """
        if i not in self.index_set():
            raise ValueError("{} is not in the index set".format(i))
        (i, j) = self.cartan_type()._indices[i]
        return self.inject_weights(i, self.ambient_spaces()[i].simple_coroot(j))

    def positive_roots(self):
        """
        EXAMPLES::

            sage: RootSystem("A1xA2").ambient_space().positive_roots()
            [(1, -1, 0, 0, 0), (0, 0, 1, -1, 0), (0, 0, 1, 0, -1), (0, 0, 0, 1, -1)]
        """
        res = []
        for i, ambient_sp in enumerate(self.ambient_spaces()):
            res.extend(self.inject_weights(i, v) for v in ambient_sp.positive_roots())
        return res

    def negative_roots(self):
        """
        EXAMPLES::

            sage: RootSystem("A1xA2").ambient_space().negative_roots()
            [(-1, 1, 0, 0, 0), (0, 0, -1, 1, 0), (0, 0, -1, 0, 1), (0, 0, 0, -1, 1)]
        """
        ret = []
        for i, ambient_sp in enumerate(self.ambient_spaces()):
            ret.extend(self.inject_weights(i, v) for v in ambient_sp.negative_roots())
        return ret

    def fundamental_weights(self):
        """
        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0, 0, 0), 2: (1, 1, 0, 0, 0), 3: (0, 0, 0, 1, 0), 4: (0, 0, 0, 1/2, 1/2)}
        """
        fw = []
        for i, ambient_sp in enumerate(self.ambient_spaces()):
            fw.extend(self.inject_weights(i, v) for v in ambient_sp.fundamental_weights())
        return Family(dict([i,fw[i-1]] for i in range(1,len(fw)+1)))


CartanType.AmbientSpace = AmbientSpace
