# sage.doctest: needs sage.combinat sage.modules
r"""
Highest weight crystals
"""

#*****************************************************************************
#       Copyright (C) 2009   Anne Schilling <anne at math.ucdavis.edu>
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
#****************************************************************************

from sage.categories.classical_crystals import ClassicalCrystals
from sage.structure.parent import Parent
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.crystals.tensor_product import TensorProductOfCrystals, \
    TensorProductOfRegularCrystalsElement

from sage.combinat.crystals.tensor_product import CrystalOfTableaux
from sage.combinat.crystals.alcove_path import CrystalOfAlcovePaths
from sage.combinat.crystals.littelmann_path import CrystalOfLSPaths
from sage.combinat.crystals.generalized_young_walls import CrystalOfGeneralizedYoungWalls
from sage.combinat.crystals.monomial_crystals import CrystalOfNakajimaMonomials
from sage.combinat.rigged_configurations.rc_crystal import CrystalOfRiggedConfigurations
from sage.rings.integer_ring import ZZ


def HighestWeightCrystal(dominant_weight, model=None):
    r"""
    Return the highest weight crystal of highest weight ``dominant_weight``
    of the given ``model``.

    INPUT:

    - ``dominant_weight`` -- a dominant weight
    - ``model`` -- (optional) if not specified, then we have the following
      default models:

      * types `A_n, B_n, C_n, D_n, G_2` - :class:`tableaux
        <sage.combinat.crystals.tensor_product.CrystalOfTableaux>`
      * types `E_{6,7}` - :class:`type E finite dimensional crystal
        <FiniteDimensionalHighestWeightCrystal_TypeE>`
      * all other types - :class:`LS paths
        <sage.combinat.crystals.littelmann_path.CrystalOfLSPaths>`

      otherwise can be one of the following:

      * ``'Tableaux'`` - :class:`KN tableaux
        <sage.combinat.crystals.tensor_product.CrystalOfTableaux>`
      * ``'TypeE'`` - :class:`type E finite dimensional crystal
        <FiniteDimensionalHighestWeightCrystal_TypeE>`
      * ``'NakajimaMonomials'`` - :class:`Nakajima monomials
        <sage.combinat.crystals.monomial_crystals.CrystalOfNakajimaMonomials>`
      * ``'LSPaths'`` - :class:`LS paths
        <sage.combinat.crystals.littelmann_path.CrystalOfLSPaths>`
      * ``'AlcovePaths'`` - :class:`alcove paths
        <sage.combinat.crystals.alcove_path.CrystalOfAlcovePaths>`
      * ``'GeneralizedYoungWalls'`` - :class:`generalized Young walls
        <sage.combinat.crystals.generalized_young_walls.CrystalOfGeneralizedYoungWalls>`
      * ``'RiggedConfigurations'`` - :class:`rigged configurations
        <sage.combinat.rigged_configurations.rc_crystal.CrystalOfRiggedConfigurations>`

    EXAMPLES::

        sage: La = RootSystem(['A',2]).weight_lattice().fundamental_weights()
        sage: wt = La[1] + La[2]
        sage: crystals.HighestWeight(wt)
        The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]

        sage: La = RootSystem(['C',2]).weight_lattice().fundamental_weights()
        sage: wt = 5*La[1] + La[2]
        sage: crystals.HighestWeight(wt)
        The crystal of tableaux of type ['C', 2] and shape(s) [[6, 1]]

        sage: La = RootSystem(['B',2]).weight_lattice().fundamental_weights()
        sage: wt = La[1] + La[2]
        sage: crystals.HighestWeight(wt)
        The crystal of tableaux of type ['B', 2] and shape(s) [[3/2, 1/2]]

    Some type `E` examples::

        sage: C = CartanType(['E',6])
        sage: La = C.root_system().weight_lattice().fundamental_weights()
        sage: T = crystals.HighestWeight(La[1])
        sage: T.cardinality()
        27
        sage: T = crystals.HighestWeight(La[6])
        sage: T.cardinality()
        27
        sage: T = crystals.HighestWeight(La[2])
        sage: T.cardinality()
        78
        sage: T = crystals.HighestWeight(La[4])
        sage: T.cardinality()
        2925
        sage: T = crystals.HighestWeight(La[3])
        sage: T.cardinality()
        351
        sage: T = crystals.HighestWeight(La[5])
        sage: T.cardinality()
        351

        sage: C = CartanType(['E',7])
        sage: La = C.root_system().weight_lattice().fundamental_weights()
        sage: T = crystals.HighestWeight(La[1])
        sage: T.cardinality()
        133
        sage: T = crystals.HighestWeight(La[2])
        sage: T.cardinality()
        912
        sage: T = crystals.HighestWeight(La[3])
        sage: T.cardinality()
        8645
        sage: T = crystals.HighestWeight(La[4])
        sage: T.cardinality()
        365750
        sage: T = crystals.HighestWeight(La[5])
        sage: T.cardinality()
        27664
        sage: T = crystals.HighestWeight(La[6])
        sage: T.cardinality()
        1539
        sage: T = crystals.HighestWeight(La[7])
        sage: T.cardinality()
        56

    An example with an affine type::

        sage: C = CartanType(['C',2,1])
        sage: La = C.root_system().weight_lattice().fundamental_weights()
        sage: T = crystals.HighestWeight(La[1])
        sage: sorted(T.subcrystal(max_depth=3), key=str)
        [(-Lambda[0] + 3*Lambda[1] - Lambda[2] - delta,),
         (-Lambda[0] + Lambda[1] + Lambda[2] - delta,),
         (-Lambda[1] + 2*Lambda[2] - delta,),
         (2*Lambda[0] - Lambda[1],),
         (Lambda[0] + Lambda[1] - Lambda[2],),
         (Lambda[0] - Lambda[1] + Lambda[2],),
         (Lambda[1],)]

    Using the various models::

        sage: La = RootSystem(['F',4]).weight_lattice().fundamental_weights()
        sage: wt = La[1] + La[4]
        sage: crystals.HighestWeight(wt)
        The crystal of LS paths of type ['F', 4] and weight Lambda[1] + Lambda[4]
        sage: crystals.HighestWeight(wt, model='NakajimaMonomials')
        Highest weight crystal of modified Nakajima monomials of
         Cartan type ['F', 4] and highest weight Lambda[1] + Lambda[4]
        sage: crystals.HighestWeight(wt, model='AlcovePaths')
        Highest weight crystal of alcove paths of type ['F', 4] and weight Lambda[1] + Lambda[4]
        sage: crystals.HighestWeight(wt, model='RiggedConfigurations')
        Crystal of rigged configurations of type ['F', 4] and weight Lambda[1] + Lambda[4]
        sage: La = RootSystem(['A',3,1]).weight_lattice().fundamental_weights()
        sage: wt = La[0] + La[2]
        sage: crystals.HighestWeight(wt, model='GeneralizedYoungWalls')
        Highest weight crystal of generalized Young walls of
         Cartan type ['A', 3, 1] and highest weight Lambda[0] + Lambda[2]

    TESTS:

    Check that the correct crystal is constructed for the fundamental weights::

        sage: for ct in CartanType.samples(finite=True, crystallographic=True):  # long time
        ....:     L = ct.root_system().weight_lattice()
        ....:     La = L.fundamental_weights()
        ....:     for model in ['Tableaux', 'NakajimaMonomials', 'AlcovePaths', 'RiggedConfigurations']:
        ....:         if model == 'Tableaux' and ct.type() in ["E", "F"]:
        ....:             continue
        ....:         for wt in La:
        ....:             C = crystals.HighestWeight(wt, model=model)
        ....:             assert L.weyl_dimension(wt) == C.cardinality(), "wrong cardinality in %s, weight %s" % (ct, wt)
        ....:             assert C.highest_weight_vector().weight() == wt, "wrong weight in %s, weight %s" % (ct, wt)

    Same thing for weights constructed from the simple roots::

        sage: for ct in CartanType.samples(finite=True, crystallographic=True):
        ....:     L = ct.root_system().root_space()
        ....:     La = L.fundamental_weights_from_simple_roots()
        ....:     for model in ['Tableaux', 'NakajimaMonomials', 'AlcovePaths', 'RiggedConfigurations']:
        ....:         if model == 'Tableaux' and ct.type() in ["E", "F"]:
        ....:             continue
        ....:         for wt in La:
        ....:             C1 = crystals.HighestWeight(wt.to_ambient().to_weight_space(ZZ), model=model)
        ....:             C2 = crystals.HighestWeight(wt, model=model)
        ....:             assert C1 == C2
    """
    cartan_type = dominant_weight.parent().cartan_type()
    if model is None:
        if cartan_type.is_finite():
            if cartan_type.type() == 'E':
                model = 'TypeE'
            elif cartan_type.type() in ['A','B','C','D','G']:
                model = 'Tableaux'
            else:
                model = 'LSPaths'
        else:
            model = 'LSPaths'

    # Make sure dominant_weight in the weight space
    if cartan_type.is_finite():
        dominant_weight = dominant_weight.to_ambient().to_weight_space(ZZ)

    if model == 'Tableaux':
        # we rely on the specific choice of positive roots here
        # except in type G_2, the fundamental weights are realized by
        # vectors with weakly decreasing nonnegative integer (or in
        # type B_n and D_n, half-integer) entries
        sh = dominant_weight.to_ambient().to_vector()
        if cartan_type.type() == "G":
            sh = (-sh)[2:0:-1]
        return CrystalOfTableaux(cartan_type, shape=sh)

    if model == 'TypeE':
        if not cartan_type.is_finite() or cartan_type.type() != 'E':
            raise ValueError("only for finite type E")
        if cartan_type.rank() == 6:
            return FiniteDimensionalHighestWeightCrystal_TypeE6(dominant_weight)
        elif cartan_type.rank() == 7:
            return FiniteDimensionalHighestWeightCrystal_TypeE7(dominant_weight)
        raise NotImplementedError

    if model == 'NakajimaMonomials':
        return CrystalOfNakajimaMonomials(cartan_type, dominant_weight)

    if model == 'LSPaths':
        # Make sure it's in the (extended) weight space
        if cartan_type.is_affine():
            P = dominant_weight.parent().root_system.weight_space(extended=True)
        else:
            P = dominant_weight.parent().root_system.weight_space()
        wt = P.sum_of_terms((i, c) for i,c in dominant_weight)
        return CrystalOfLSPaths(wt)

    if model == 'AlcovePaths':
        return CrystalOfAlcovePaths(dominant_weight, highest_weight_crystal=True)

    if model == 'GeneralizedYoungWalls':
        if not cartan_type.is_affine():
            raise ValueError("only for affine types")
        if cartan_type.type() != 'A':
            raise NotImplementedError("only for affine type A")
        # Make sure it's in the weight lattice
        P = dominant_weight.parent().root_system.weight_lattice(extended=True)
        wt = P.sum_of_terms((i, c) for i,c in dominant_weight)
        return CrystalOfGeneralizedYoungWalls(cartan_type.rank()-1, wt)

    if model == 'RiggedConfigurations':
        return CrystalOfRiggedConfigurations(cartan_type, dominant_weight)

    raise ValueError("invalid model")


class FiniteDimensionalHighestWeightCrystal_TypeE(TensorProductOfCrystals):
    """
    Commonalities for all finite dimensional type `E` highest weight crystals.

    Subclasses should setup an attribute column_crystal in their
    ``__init__`` method before calling the ``__init__`` method of this class.
    """

    def __init__(self, dominant_weight):
        """
        EXAMPLES::

            sage: C = CartanType(['E',6])
            sage: La = C.root_system().weight_lattice().fundamental_weights()
            sage: T = crystals.HighestWeight(2*La[2])
            sage: T.cartan_type()
            ['E', 6]
            sage: T.module_generators
            [[[(2, -1), (1,)], [(2, -1), (1,)]]]
            sage: T.cardinality()
            2430
            sage: T = crystals.HighestWeight(La[2])
            sage: T.cardinality()
            78
        """
        self._cartan_type = dominant_weight.parent().cartan_type()
        self._highest_weight = dominant_weight
        assert dominant_weight.is_dominant()
        self.rename()
        Parent.__init__(self, category=ClassicalCrystals())
        self.module_generators = [self.module_generator()]

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = CartanType(['E',6])
            sage: La = C.root_system().weight_lattice().fundamental_weights()
            sage: crystals.HighestWeight(2*La[2])
            Finite dimensional highest weight crystal of type ['E', 6] and highest weight 2*Lambda[2]
        """
        return "Finite dimensional highest weight crystal of type {} and highest weight {}".format(
                self._cartan_type, self._highest_weight)

    Element = TensorProductOfRegularCrystalsElement

    def module_generator(self):
        """
        Yield the module generator (or highest weight element) of the classical
        crystal of given dominant weight in ``self``.

        EXAMPLES::

            sage: C=CartanType(['E',6])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: T = crystals.HighestWeight(La[2])
            sage: T.module_generator()
            [[(2, -1), (1,)]]
            sage: T = crystals.HighestWeight(0*La[2])
            sage: T.module_generator()
            []

            sage: C=CartanType(['E',7])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: T = crystals.HighestWeight(La[1])
            sage: T.module_generator()
            [[(-7, 1), (7,)]]
        """
        dominant_weight = self._highest_weight
        tensor = sum(( [self.column_crystal[i]]*dominant_weight.coefficient(i) for i in dominant_weight.support()), [])
        return self._element_constructor_(*[B.module_generators[0] for B in tensor])


class FiniteDimensionalHighestWeightCrystal_TypeE6(FiniteDimensionalHighestWeightCrystal_TypeE):
    r"""
    Class of finite dimensional highest weight crystals of type `E_6`.

    EXAMPLES::

        sage: C=CartanType(['E',6])
        sage: La=C.root_system().weight_lattice().fundamental_weights()
        sage: T = crystals.HighestWeight(La[2]); T
        Finite dimensional highest weight crystal of type ['E', 6] and highest weight Lambda[2]
        sage: B1 = T.column_crystal[1]; B1
        The crystal of letters for type ['E', 6]
        sage: B6 = T.column_crystal[6]; B6
        The crystal of letters for type ['E', 6] (dual)
        sage: t = T(B6([-1]),B1([-1,3])); t
        [(-1,), (-1, 3)]
        sage: [t.epsilon(i) for i in T.index_set()]
        [2, 0, 0, 0, 0, 0]
        sage: [t.phi(i) for i in T.index_set()]
        [0, 0, 1, 0, 0, 0]
        sage: TestSuite(t).run()
    """

    def __init__(self, dominant_weight):
        """
        EXAMPLES::

            sage: C=CartanType(['E',6])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: p2=2*La[2]
            sage: p1=La[2]
            sage: p0=0*La[2]
            sage: T = crystals.HighestWeight(0*La[2])
            sage: T.cardinality()
            1
            sage: T = crystals.HighestWeight(La[2])
            sage: T.cardinality()
            78
            sage: T = crystals.HighestWeight(2*La[2])
            sage: T.cardinality()
            2430
        """
        B1 = CrystalOfLetters(['E',6])
        B6 = CrystalOfLetters(['E',6], dual=True)
        self.column_crystal = {1 : B1, 6 : B6,
                               4 : TensorProductOfCrystals(B1,B1,B1,generators=[[B1([-3,4]),B1([-1,3]),B1([1])]]),
                               3 : TensorProductOfCrystals(B1,B1,generators=[[B1([-1,3]),B1([1])]]),
                               5 : TensorProductOfCrystals(B6,B6,generators=[[B6([5,-6]),B6([6])]]),
                               2 : TensorProductOfCrystals(B6,B1,generators=[[B6([2,-1]),B1([1])]])}
        FiniteDimensionalHighestWeightCrystal_TypeE.__init__(self, dominant_weight)


class FiniteDimensionalHighestWeightCrystal_TypeE7(FiniteDimensionalHighestWeightCrystal_TypeE):
    r"""
    Class of finite dimensional highest weight crystals of type `E_7`.

    EXAMPLES::

        sage: C=CartanType(['E',7])
        sage: La=C.root_system().weight_lattice().fundamental_weights()
        sage: T = crystals.HighestWeight(La[1])
        sage: T.cardinality()
        133
        sage: B7 = T.column_crystal[7]; B7
        The crystal of letters for type ['E', 7]
        sage: t = T(B7([-5, 6]), B7([-2, 3])); t
        [(-5, 6), (-2, 3)]
        sage: [t.epsilon(i) for i in T.index_set()]
        [0, 1, 0, 0, 1, 0, 0]
        sage: [t.phi(i) for i in T.index_set()]
        [0, 0, 1, 0, 0, 1, 0]
        sage: TestSuite(t).run()
    """

    def __init__(self, dominant_weight):
        """
        EXAMPLES::

            sage: C=CartanType(['E',7])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: T = crystals.HighestWeight(0*La[1])
            sage: T.cardinality()
            1
            sage: T = crystals.HighestWeight(La[1])
            sage: T.cardinality()
            133
            sage: T = crystals.HighestWeight(2*La[1])
            sage: T.cardinality()
            7371
        """
        B = CrystalOfLetters(['E',7])
        self.column_crystal = {7 : B,
                               1 : TensorProductOfCrystals(B,B,generators=[[B([-7,1]),B([7])]]),
                               2 : TensorProductOfCrystals(B,B,B,generators=[[B([-1,2]),B([-7,1]),B([7])]]),
                               3 : TensorProductOfCrystals(B,B,B,B,generators=[[B([-2,3]),B([-1,2]),B([-7,1]),B([7])]]),
                               4 : TensorProductOfCrystals(B,B,B,B,generators=[[B([-5,4]),B([-6,5]),B([-7,6]),B([7])]]),
                               5 : TensorProductOfCrystals(B,B,B,generators=[[B([-6,5]),B([-7,6]),B([7])]]),
                               6 : TensorProductOfCrystals(B,B,generators=[[B([-7,6]),B([7])]])}
        FiniteDimensionalHighestWeightCrystal_TypeE.__init__(self, dominant_weight)
