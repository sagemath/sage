# sage.doctest: needs sage.combinat sage.modules
r"""
Fast Rank Two Crystals
"""
# ****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
#                          Ben Brubaker   <brubaker at math.mit.edu>
#                          Daniel Bump    <bump at match.stanford.edu>
#                          Justin Walker  <justin at mac.com>
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
# ***************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.classical_crystals import ClassicalCrystals
from sage.structure.element import Element, parent
from sage.combinat.root_system.cartan_type import CartanType
from sage.structure.richcmp import richcmp


class FastCrystal(UniqueRepresentation, Parent):
    r"""
    An alternative implementation of rank 2 crystals. The root
    operators are implemented in memory by table lookup. This means
    that in comparison with the
    :class:`~sage.combinat.crystals.tensor_product.CrystalsOfTableaux`, these
    crystals are slow to instantiate but faster for computation. Implemented
    for types `A_2`, `B_2`, and `C_2`.

    INPUT:

    - ``cartan_type`` -- the Cartan type and must be either type `A_2`, `B_2`, or `C_2`

    - ``shape`` -- a shape is of the form ``[l1,l2]`` where ``l1`` and ``l2``
      are either integers or (in type `B_2`) half integers such that
      ``l1 - l2`` is integral. It is assumed that ``l1 >= l2 >= 0``. If
      ``l1`` and ``l2`` are integers, this will produce a crystal
      isomorphic to the one obtained by
      ``crystals.Tableaux(type, shape=[l1,l2])``. Furthermore
      ``crystals.FastRankTwo(['B', 2], l1+1/2, l2+1/2)`` produces a crystal
      isomorphic to the following crystal ``T``::

          sage: C = crystals.Tableaux(['B',2], shape=[l1,l2])               # not tested
          sage: D = crystals.Spins(['B',2])                                 # not tested
          sage: T = crystals.TensorProduct(C, D, C.list()[0], D.list()[0])  # not tested

    - ``format`` -- (default: ``'string'``) the default representation of
      elements is in term of theBerenstein-Zelevinsky-Littelmann (BZL)
      strings ``[a1, a2, ...]`` described under metapost in
      :mod:`~sage.categories.crystals`. Alternative representations may be
      obtained by the options ``'dual_string'`` or ``'simple'``.
      In the ``'simple'`` format, the element is represented by and integer,
      and in the ``'dual_string'`` format, it is represented by the
      BZL string, but the underlying decomposition of the long Weyl group
      element into simple reflections is changed.

    TESTS::

        sage: C = crystals.FastRankTwo(['A',2],shape=[4,1])
        sage: C.cardinality()
        24
        sage: C.cartan_type()
        ['A', 2]
        sage: TestSuite(C).run()
        sage: C = crystals.FastRankTwo(['B',2],shape=[4,1])
        sage: C.cardinality()
        154
        sage: TestSuite(C).run()
        sage: C = crystals.FastRankTwo(['B',2],shape=[3/2,1/2])
        sage: C.cardinality()
        16
        sage: TestSuite(C).run()
        sage: C = crystals.FastRankTwo(['C',2],shape=[2,1])
        sage: C.cardinality()
        16
        sage: C = crystals.FastRankTwo(['C',2],shape=[3,1])
        sage: C.cardinality()
        35
        sage: TestSuite(C).run()

        sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
        sage: C.list()
        [[0, 0, 0],
         [1, 0, 0],
         [0, 1, 1],
         [0, 2, 1],
         [1, 2, 1],
         [0, 1, 0],
         [1, 1, 0],
         [2, 1, 0]]
    """
    @staticmethod
    def __classcall__(cls, cartan_type, shape, format='string'):
        """
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: C1 = crystals.FastRankTwo(['A',2],            shape=(4,1))
            sage: C2 = crystals.FastRankTwo(CartanType(['A',2]),shape=[4,1])
            sage: C1 is C2
            True
        """
        cartan_type = CartanType(cartan_type)
        shape = tuple(shape)
        if len(shape) > 2:
            raise ValueError("The shape must have length <=2")
        shape = shape + (0,)*(2-len(shape))
        return super().__classcall__(cls, cartan_type, shape, format)

    def __init__(self, ct, shape, format):
        """
        EXAMPLES::

            sage: C = crystals.FastRankTwo(['A',2],shape=[4,1]); C
            The fast crystal for A2 with shape [4,1]
            sage: TestSuite(C).run()
        """
        Parent.__init__(self, category=ClassicalCrystals())
#        super().__init__(category = FiniteEnumeratedSets())
        self._cartan_type = ct
        if ct[1] != 2:
            raise NotImplementedError

        l1 = shape[0]
        l2 = shape[1]

        # For safety, delpat and gampat should be immutable

        self.delpat = []
        self.gampat = []

        if ct[0] == 'A':
            self._type_a_init(l1, l2)
        elif ct[0] == 'B' or ct[0] == 'C':
            self._type_bc_init(l1, l2)
        else:
            raise NotImplementedError

        self.format = format
        self.size = len(self.delpat)
        self._rootoperators = []
        self.shape = shape

        for i in range(self.size):
            target = list(self.delpat[i])

            target[0] = target[0]-1
            e1 = None if target not in self.delpat else self.delpat.index(target)
            target[0] = target[0]+1+1
            f1 = None if target not in self.delpat else self.delpat.index(target)

            target = list(self.gampat[i])
            target[0] = target[0]-1
            e2 = None if target not in self.gampat else self.gampat.index(target)
            target[0] = target[0]+1+1
            f2 = None if target not in self.gampat else self.gampat.index(target)

            self._rootoperators.append([e1,f1,e2,f2])

        if int(2*l1) % 2 == 0:
            l1_str = "%d" % l1
            l2_str = "%d" % l2
        else:
            assert self._cartan_type[0] == 'B' and int(2*l2) % 2 == 1
            l1_str = "%d/2" % int(2*l1)
            l2_str = "%d/2" % int(2*l2)
        self.rename("The fast crystal for %s2 with shape [%s,%s]" % (ct[0],l1_str,l2_str))
        self.module_generators = [self(0)]
        # self._digraph = ClassicalCrystal.digraph(self)
        self._digraph = super().digraph()
        self._digraph_closure = self.digraph().transitive_closure()

    def _type_a_init(self, l1, l2):
        """
        EXAMPLES::

            sage: C = crystals.FastRankTwo(['A',2],shape=[1,1])
            sage: C.delpat # indirect doctest
            [[0, 0, 0], [0, 1, 0], [1, 1, 0]]
            sage: C.gampat
            [[0, 0, 0], [1, 0, 0], [0, 1, 1]]
        """
        for b in range(l2,-1,-1):
            for a in range(l1,l2-1,-1):
                for c in range(a,b-1,-1):
                    a3 = l1-a
                    a2 = l1+l2-a-b
                    a1 = a-c
                    b1 = max(a3,a2-a1)
                    b2 = a1+a3
                    b3 = min(a2-a3,a1)
                    self.delpat.append([a1,a2,a3])
                    self.gampat.append([b1,b2,b3])

    def _type_bc_init(self, l1, l2):
        """
        EXAMPLES::

            sage: C = crystals.FastRankTwo(['B',2],shape=[1])
            sage: len(C.delpat) # indirect doctest
            5
            sage: len(C.gampat)
            5
            sage: C = crystals.FastRankTwo(['C',2],shape=[1])
            sage: len(C.delpat)
            4
            sage: len(C.gampat)
            4
        """
        if self._cartan_type[0] == 'B':
            [m1, m2] = [l1+l2, l1-l2]
        else:
            [m1, m2] = [l1, l2]
        for b in range(m2,-1,-1):
            for a in range(m1,m2-1,-1):
                for c in range(b,a+1):
                    for d in range(c,-1,-1):
                        a1 = c-d
                        a2 = m1+m2+c-a-2*b
                        a3 = m1+m2-a-b
                        a4 = m1-a
                        b1 = max(a4,2*a3-a2,a2-2*a1)
                        b2 = max(a3, a1+a4, a1+2*a3-a2)
                        b3 = min(a2, 2*a2-2*a3+a4, 2*a1+a4)
                        b4 = min(a1, a2-a3, a3-a4)
                        if self._cartan_type[0] == 'B':
                            self.delpat.append([a1,a2,a3,a4])
                            self.gampat.append([b1,b2,b3,b4])
                        else:
                            self.gampat.append([a1,a2,a3,a4])
                            self.delpat.append([b1,b2,b3,b4])

    def __call__(self, value):
        """
        EXAMPLES::

            sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
            sage: C(0)
            [0, 0, 0]
            sage: C(1)
            [1, 0, 0]
            sage: x = C(0)
            sage: C(x) is x
            True
        """
        if parent(value) is self:
            return value
        return self.element_class(self, value, self.format)

    def digraph(self):
        """
        Return the digraph associated to ``self``.

        EXAMPLES::

            sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
            sage: C.digraph()
            Digraph on 8 vertices
        """
        return self._digraph

    def cmp_elements(self, x, y):
        r"""
        Return ``True`` if and only if there is a path from `x` to `y` in the
        crystal graph.

        Because the crystal graph is classical, it is a directed acyclic
        graph which can be interpreted as a poset. This function implements
        the comparison function of this poset.

        EXAMPLES::

            sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
            sage: x = C(0)
            sage: y = C(1)
            sage: C.cmp_elements(x,y)
            -1
            sage: C.cmp_elements(y,x)
            1
            sage: C.cmp_elements(x,x)
            0
        """
        assert x.parent() == self and y.parent() == self
        if self._digraph_closure.has_edge(x,y):
            return -1
        elif self._digraph_closure.has_edge(y,x):
            return 1
        else:
            return 0

    class Element(Element):
        def __init__(self, parent, value, format):
            """
            EXAMPLES::

                sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
                sage: c = C(0); c
                [0, 0, 0]
                sage: C[0].parent()
                The fast crystal for A2 with shape [2,1]
                sage: TestSuite(c).run()
            """
            Element.__init__(self, parent)
            self.value = value
            self.format = format

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: [v.weight() for v in crystals.FastRankTwo(['A',2], shape=[2,1])]
                [(2, 1, 0), (1, 2, 0), (1, 1, 1), (1, 0, 2), (0, 1, 2), (2, 0, 1), (1, 1, 1), (0, 2, 1)]
                sage: [v.weight() for v in crystals.FastRankTwo(['B',2], shape=[1,0])]
                [(1, 0), (0, 1), (0, 0), (0, -1), (-1, 0)]
                sage: [v.weight() for v in crystals.FastRankTwo(['B',2], shape=[1/2,1/2])]
                [(1/2, 1/2), (1/2, -1/2), (-1/2, 1/2), (-1/2, -1/2)]
                sage: [v.weight() for v in crystals.FastRankTwo(['C',2], shape=[1,0])]
                [(1, 0), (0, 1), (0, -1), (-1, 0)]
                sage: [v.weight() for v in crystals.FastRankTwo(['C',2], shape=[1,1])]
                [(1, 1), (1, -1), (0, 0), (-1, 1), (-1, -1)]
            """
            delpat = self.parent().delpat[self.value]
            if self.parent()._cartan_type[0] == 'A':
                delpat = delpat + [0,]
            [alpha1, alpha2] = self.parent().weight_lattice_realization().simple_roots()
            hwv = sum(self.parent().shape[i]*self.parent().weight_lattice_realization().monomial(i) for i in range(2))
            return hwv - (delpat[0]+delpat[2])*alpha1 - (delpat[1]+delpat[3])*alpha2

        def _repr_(self):
            """
            EXAMPLES::

                sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
                sage: C[0]._repr_()
                '[0, 0, 0]'
            """
            if self.format == "string":
                return repr(self.parent().delpat[self.value])
            elif self.format == "dual_string":
                return repr(self.parent().gampat[self.value])
            elif self.format == "simple":
                return repr(self.value)
            else:
                raise NotImplementedError

        def __hash__(self):
            r"""
            TESTS::

                sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
                sage: hash(C(0))
                0
            """
            return hash(self.value)

        def _richcmp_(self, other, op):
            """
            EXAMPLES::

                sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
                sage: D = crystals.FastRankTwo(['B',2],shape=[2,1])
                sage: C(0) == C(0)
                True
                sage: C(1) == C(0)
                False
                sage: C(0) == D(0)
                False

                sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
                sage: D = crystals.FastRankTwo(['B',2],shape=[2,1])
                sage: C(0) != C(0)
                False
                sage: C(1) != C(0)
                True
                sage: C(0) != D(0)
                True

                sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
                sage: C(1) < C(2)
                True
                sage: C(2) < C(1)
                False
                sage: C(2) > C(1)
                True
                sage: C(1) <= C(1)
                True
            """
            return richcmp(self.value, other.value, op)

        def e(self, i):
            """
            Return the action of `e_i` on ``self``.

            EXAMPLES::

                sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
                sage: C(1).e(1)
                [0, 0, 0]
                sage: C(0).e(1) is None
                True
            """
            assert i in self.index_set()
            if i == 1:
                r = self.parent()._rootoperators[self.value][0]
            else:
                r = self.parent()._rootoperators[self.value][2]
            return self.parent()(r) if r is not None else None

        def f(self, i):
            """
            Return the action of `f_i` on ``self``.

            EXAMPLES::

                sage: C = crystals.FastRankTwo(['A',2],shape=[2,1])
                sage: C(6).f(1)
                [1, 2, 1]
                sage: C(7).f(1) is None
                True
            """
            assert i in self.index_set()
            if i == 1:
                r = self.parent()._rootoperators[self.value][1]
            else:
                r = self.parent()._rootoperators[self.value][3]
            return self.parent()(r) if r is not None else None
