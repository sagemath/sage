r"""
Example of a crystal
"""
#*****************************************************************************
#  Copyright (C) 2010 Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.graphs.all import DiGraph
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType

class HighestWeightCrystalOfTypeA(UniqueRepresentation, Parent):
    r"""
    An example of a crystal: the highest weight crystal of type `A_n`
    of highest weight `\omega_1`.

    The purpose of this class is to provide a minimal template for
    implementing crystals. See
    :class:`~sage.combinat.crystals.letters.CrystalOfLetters` for a
    full featured and optimized implementation.

    EXAMPLES::

        sage: C = Crystals().example()
        sage: C
        Highest weight crystal of type A_3 of highest weight omega_1
        sage: C.category()
        Category of classical crystals

    The elements of this crystal are in the set `\{1,\ldots,n+1\}`::

        sage: C.list()
        [1, 2, 3,  4]
        sage: C.module_generators[0]
        1

    The crystal operators themselves correspond to the elementary
    transpositions::

        sage: b = C.module_generators[0]
        sage: b.f(1)
        2
        sage: b.f(1).e(1) == b
        True

    Only the following basic operations are implemented:

    - :meth:`~sage.categories.crystals.Crystals.cartan_type` or an attribute _cartan_type
    - an attribute module_generators
    - :meth:`.Element.e`
    - :meth:`.Element.f`

    All the other usual crystal operations are inherited from the
    categories; for example::

        sage: C.cardinality()
        4

    TESTS::

        sage: C = Crystals().example()
        sage: TestSuite(C).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          running ._test_stembridge_local_axioms() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_stembridge_local_axioms() . . . pass
    """

    def __init__(self, n = 3):
        """
        EXAMPLES::

            sage: C = sage.categories.examples.crystals.HighestWeightCrystalOfTypeA(n=4)
            sage: C == Crystals().example(n=4)
            True
        """
        Parent.__init__(self, category = ClassicalCrystals())
        self.n = n
        self._cartan_type = CartanType(['A',n])
        self.module_generators = [ self(1) ]

    def _repr_(self):
        """
        EXAMPLES::

            sage: Crystals().example()
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        return "Highest weight crystal of type A_%s of highest weight omega_1"%(self.n)

    # temporary workaround while an_element is overridden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

    class Element(ElementWrapper):

        def e(self, i):
            r"""
            Returns the action of `e_i` on ``self``.

            EXAMPLES::

                sage: C = Crystals().example(4)
                sage: [[c,i,c.e(i)] for i in C.index_set() for c in C if c.e(i) is not None]
                [[2, 1, 1], [3, 2, 2], [4, 3, 3], [5, 4, 4]]
            """
            assert i in self.index_set()
            if self.value == i+1:
                return self.parent()(self.value-1)
            else:
                return None

        def f(self, i):
            r"""
            Returns the action of `f_i` on ``self``.

            EXAMPLES::

                sage: C = Crystals().example(4)
                sage: [[c,i,c.f(i)] for i in C.index_set() for c in C if c.f(i) is not None]
                [[1, 1, 2], [2, 2, 3], [3, 3, 4], [4, 4, 5]]
            """
            assert i in self.index_set()
            if self.value == i:
                return self.parent()(self.value+1)
            else:
                return None


class NaiveCrystal(UniqueRepresentation, Parent):
    r"""
    This is an example of a "crystal" which does not come from any kind of
    representation, designed primarily to test the Stembridge local rules with.
    The crystal has vertices labeled 0 through 5, with 0 the highest weight.

    The code here could also possibly be generalized to create a class that
    automatically builds a crystal from an edge-colored digraph, if someone
    feels adventurous.

    Currently, only the methods :meth:`highest_weight_vector`, :meth:`e`, and :meth:`f` are
    guaranteed to work.

    EXAMPLES::

        sage: C = Crystals().example(choice='naive')
        sage: C.highest_weight_vector()
        0
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: C = sage.categories.examples.crystals.NaiveCrystal()
            sage: C == Crystals().example(choice='naive')
            True
        """
        Parent.__init__(self, category = ClassicalCrystals())
        self.n = 2
        self._cartan_type = CartanType(['A',2])
        self.G = DiGraph(5)
        self.G.add_edges([ [0,1,1], [1,2,1], [2,3,1], [3,5,1],  [0,4,2], [4,5,2] ])
        self.module_generators = [ self(0) ]

    def __repr__(self):
        """
        EXAMPLES::

            sage: Crystals().example(choice='naive')
            A broken crystal, defined by digraph, of dimension five.
        """
        return "A broken crystal, defined by digraph, of dimension five."

    class Element(ElementWrapper):
        def e(self, i):
            r"""
            Returns the action of `e_i` on ``self``.

            EXAMPLES::

                sage: C = Crystals().example(choice='naive')
                sage: [[c,i,c.e(i)] for i in C.index_set() for c in [C(j) for j in [0..5]] if c.e(i) is not None]
                [[1, 1, 0], [2, 1, 1], [3, 1, 2], [5, 1, 3], [4, 2, 0], [5, 2, 4]]
            """
            assert i in self.index_set()
            for edge in self.parent().G.edges():
               if edge[1]==int(str(self)) and edge[2]==i:
                   return self.parent()(edge[0])
            return None

        def f(self, i):
            r"""
            Returns the action of `f_i` on ``self``.

            EXAMPLES::

                sage: C = Crystals().example(choice='naive')
                sage: [[c,i,c.f(i)] for i in C.index_set() for c in [C(j) for j in [0..5]] if c.f(i) is not None]
                [[0, 1, 1], [1, 1, 2], [2, 1, 3], [3, 1, 5], [0, 2, 4], [4, 2, 5]]
            """
            assert i in self.index_set()
            for edge in self.parent().G.edges_incident(int(str(self))):
                if edge[2] == i:
                    return self.parent()(edge[1])
            return None


