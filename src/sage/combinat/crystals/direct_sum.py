# sage.doctest: needs sage.combinat sage.modules
"""
Direct sum of crystals
"""
# ****************************************************************************
#       Copyright (C) 2010 Anne Schilling <anne at math.ucdavis.edu>
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
#                  http://www.gnu.org/licenses/
# ***************************************************************************

from sage.categories.category import Category
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.element import get_coercion_model


class DirectSumOfCrystals(DisjointUnionEnumeratedSets):
    r"""
    Direct sum of crystals.

    Given a list of crystals `B_0, \ldots, B_k` of the same Cartan type,
    one can form the direct sum `B_0 \oplus \cdots \oplus B_k`.

    INPUT:

    - ``crystals`` -- list of crystals of the same Cartan type
    - ``keepkey`` -- boolean

    The option ``keepkey`` is by default set to ``False``, assuming
    that the crystals are all distinct. In this case the elements of
    the direct sum are just represented by the elements in the
    crystals `B_i`.  If the crystals are not all distinct, one should
    set the ``keepkey`` option to ``True``.  In this case, the
    elements of the direct sum are represented as tuples `(i, b)`
    where `b \in B_i`.

    EXAMPLES::

        sage: C = crystals.Letters(['A',2])
        sage: C1 = crystals.Tableaux(['A',2],shape=[1,1])
        sage: B = crystals.DirectSum([C,C1])
        sage: B.list()
        [1, 2, 3, [[1], [2]], [[1], [3]], [[2], [3]]]
        sage: [b.f(1) for b in B]
        [2, None, None, None, [[2], [3]], None]
        sage: B.module_generators
        (1, [[1], [2]])

    ::

        sage: B = crystals.DirectSum([C,C], keepkey=True)
        sage: B.list()
        [(0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3)]
        sage: B.module_generators
        ((0, 1), (1, 1))
        sage: b = B( tuple([0,C(1)]) )
        sage: b
        (0, 1)
        sage: b.weight()
        (1, 0, 0)

    The following is required, because
    :class:`~sage.combinat.crystals.direct_sum.DirectSumOfCrystals`
    takes the same arguments as :class:`DisjointUnionEnumeratedSets`
    (which see for details).

    TESTS::

        sage: C = crystals.Letters(['A',2])
        sage: B = crystals.DirectSum([C,C], keepkey=True)
        sage: B
        Direct sum of the crystals Family (The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2])

        sage: TestSuite(C).run()
    """
    @staticmethod
    def __classcall_private__(cls, crystals, facade=True, keepkey=False, category=None):
        """
        Normalization of arguments; see :class:`UniqueRepresentation`.

        TESTS:

        We check that direct sum of crystals have unique representation::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: C = crystals.Letters(['A',2])
            sage: D1 = crystals.DirectSum([B, C])
            sage: D2 = crystals.DirectSum((B, C))
            sage: D1 is D2
            True
            sage: D3 = crystals.DirectSum([B, C, C])
            sage: D4 = crystals.DirectSum([D1, C])
            sage: D3 is D4
            True
        """
        if not isinstance(facade, bool) or not isinstance(keepkey, bool):
            raise TypeError
        # Normalize the facade-keepkey by giving keepkey dominance
        facade = not keepkey

        # We expand out direct sums of crystals
        ret = []
        for x in Family(crystals):
            if isinstance(x, DirectSumOfCrystals):
                ret += list(x.crystals)
            else:
                ret.append(x)
        category = Category.meet([Category.join(c.categories()) for c in ret])
        return super().__classcall__(cls,
            Family(ret), facade=facade, keepkey=keepkey, category=category)

    def __init__(self, crystals, facade, keepkey, category, **options):
        """
        TESTS::

            sage: C = crystals.Letters(['A',2])
            sage: B = crystals.DirectSum([C,C], keepkey=True)
            sage: B
            Direct sum of the crystals Family (The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2])
            sage: B.cartan_type()
            ['A', 2]

            sage: from sage.combinat.crystals.direct_sum import DirectSumOfCrystals
            sage: isinstance(B, DirectSumOfCrystals)
            True
        """
        DisjointUnionEnumeratedSets.__init__(self, crystals, keepkey=keepkey,
                                             facade=facade, category=category)
        self.rename("Direct sum of the crystals {}".format(crystals))
        self._keepkey = keepkey
        self.crystals = crystals
        if len(crystals) == 0:
            raise ValueError("the direct sum is empty")
        else:
            assert all(crystal.cartan_type() == crystals[0].cartan_type() for crystal in crystals)
            self._cartan_type = crystals[0].cartan_type()
        if keepkey:
            self.module_generators = tuple([self((i, b))
                                            for i, B in enumerate(crystals)
                                            for b in B.module_generators])
        else:
            self.module_generators = sum((tuple(B.module_generators) for B in crystals), ())

    def weight_lattice_realization(self):
        r"""
        Return the weight lattice realization used to express weights.

        The weight lattice realization is the common parent which all
        weight lattice realizations of the crystals of ``self`` coerce
        into.

        EXAMPLES::

            sage: LaZ = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: LaQ = RootSystem(['A',2,1]).weight_space(extended=True).fundamental_weights()
            sage: B = crystals.LSPaths(LaQ[1])
            sage: B.weight_lattice_realization()
            Extended weight space over the Rational Field of the Root system of type ['A', 2, 1]
            sage: C = crystals.AlcovePaths(LaZ[1])
            sage: C.weight_lattice_realization()
            Extended weight lattice of the Root system of type ['A', 2, 1]
            sage: D = crystals.DirectSum([B,C])
            sage: D.weight_lattice_realization()
            Extended weight space over the Rational Field of the Root system of type ['A', 2, 1]
        """
        cm = get_coercion_model()
        return cm.common_parent(*[crystal.weight_lattice_realization()
                                  for crystal in self.crystals])

    class Element(ElementWrapper):
        r"""
        A class for elements of direct sums of crystals.
        """

        def e(self, i):
            r"""
            Return the action of `e_i` on ``self``.

            EXAMPLES::

                sage: C = crystals.Letters(['A',2])
                sage: B = crystals.DirectSum([C,C], keepkey=True)
                sage: [[b, b.e(2)] for b in B]
                [[(0, 1), None], [(0, 2), None], [(0, 3), (0, 2)], [(1, 1), None], [(1, 2), None], [(1, 3), (1, 2)]]
            """
            v = self.value
            vn = v[1].e(i)
            if vn is None:
                return None
            else:
                return self.parent()(tuple([v[0],vn]))

        def f(self, i):
            r"""
            Return the action of `f_i` on ``self``.

            EXAMPLES::

                sage: C = crystals.Letters(['A',2])
                sage: B = crystals.DirectSum([C,C], keepkey=True)
                sage: [[b,b.f(1)] for b in B]
                [[(0, 1), (0, 2)], [(0, 2), None], [(0, 3), None], [(1, 1), (1, 2)], [(1, 2), None], [(1, 3), None]]
            """
            v = self.value
            vn = v[1].f(i)
            if vn is None:
                return None
            else:
                return self.parent()(tuple([v[0],vn]))

        def weight(self):
            r"""
            Return the weight of ``self``.

            EXAMPLES::

                sage: C = crystals.Letters(['A',2])
                sage: B = crystals.DirectSum([C,C], keepkey=True)
                sage: b = B( tuple([0,C(2)]) )
                sage: b
                (0, 2)
                sage: b.weight()
                (0, 1, 0)
            """
            return self.value[1].weight()

        def phi(self, i):
            r"""
            EXAMPLES::

                sage: C = crystals.Letters(['A',2])
                sage: B = crystals.DirectSum([C,C], keepkey=True)
                sage: b = B( tuple([0,C(2)]) )
                sage: b.phi(2)
                1
            """
            return self.value[1].phi(i)

        def epsilon(self, i):
            r"""
            EXAMPLES::

                sage: C = crystals.Letters(['A',2])
                sage: B = crystals.DirectSum([C,C], keepkey=True)
                sage: b = B( tuple([0,C(2)]) )
                sage: b.epsilon(2)
                0
            """
            return self.value[1].epsilon(i)
