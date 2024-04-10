# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Crystal of Rigged Configurations

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

We only consider the highest weight crystal structure, not the
Kirillov-Reshetikhin structure, and we extend this to symmetrizable types.
"""

# ****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
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
from itertools import repeat
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
from sage.combinat.rigged_configurations.rigged_configuration_element import (
     RiggedConfigurationElement, RCHighestWeightElement, RCHWNonSimplyLacedElement)
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition

# Note on implementation, this class is used for simply-laced types only


class CrystalOfRiggedConfigurations(UniqueRepresentation, Parent):
    r"""
    A highest weight crystal of rigged configurations.

    The crystal structure for finite simply-laced types is given
    in [CrysStructSchilling06]_. These were then shown to be the crystal
    operators in all finite types in [SS2015]_, all simply-laced and
    a large class of foldings of simply-laced types in [SS2015II]_,
    and all symmetrizable types (uniformly) in [SS2017]_.

    INPUT:

    - ``cartan_type`` -- (optional) a Cartan type or a Cartan type
      given as a folding

    - ``wt`` -- the highest weight vector in the weight lattice

    EXAMPLES:

    For simplicity, we display the rigged configurations horizontally::

        sage: RiggedConfigurations.options.display='horizontal'

    We start with a simply-laced finite type::

        sage: La = RootSystem(['A', 2]).weight_lattice().fundamental_weights()
        sage: RC = crystals.RiggedConfigurations(La[1] + La[2])
        sage: mg = RC.highest_weight_vector()
        sage: mg.f_string([1,2])
        0[ ]0   0[ ]-1
        sage: mg.f_string([1,2,2])
        0[ ]0   -2[ ][ ]-2
        sage: mg.f_string([1,2,2,2])
        sage: mg.f_string([2,1,1,2])
        -1[ ][ ]-1   -1[ ][ ]-1
        sage: RC.cardinality()
        8
        sage: T = crystals.Tableaux(['A', 2], shape=[2,1])
        sage: RC.digraph().is_isomorphic(T.digraph(), edge_labels=True)
        True

    We construct a non-simply-laced affine type::

        sage: La = RootSystem(['C', 3]).weight_lattice().fundamental_weights()
        sage: RC = crystals.RiggedConfigurations(La[2])
        sage: mg = RC.highest_weight_vector()
        sage: mg.f_string([2,3])
        (/)   1[ ]1   -1[ ]-1
        sage: T = crystals.Tableaux(['C', 3], shape=[1,1])
        sage: RC.digraph().is_isomorphic(T.digraph(), edge_labels=True)
        True

    We can construct rigged configurations using a diagram folding of
    a simply-laced type. This yields an equivalent but distinct crystal::

        sage: vct = CartanType(['C', 3]).as_folding()
        sage: RC = crystals.RiggedConfigurations(vct, La[2])
        sage: mg = RC.highest_weight_vector()
        sage: mg.f_string([2,3])
        (/)   0[ ]0   -1[ ]-1
        sage: T = crystals.Tableaux(['C', 3], shape=[1,1])
        sage: RC.digraph().is_isomorphic(T.digraph(), edge_labels=True)
        True

    We reset the global options::

        sage: RiggedConfigurations.options._reset()

    REFERENCES:

    - [SS2015]_
    - [SS2015II]_
    - [SS2017]_
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, wt=None, WLR=None):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: La = RootSystem(['A', 2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1])
            sage: RC2 = crystals.RiggedConfigurations(['A', 2], La[1])
            sage: RC3 = crystals.RiggedConfigurations(['A', 2], La[1], La[1].parent())
            sage: RC is RC2 and RC2 is RC3
            True

            sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
            sage: LaE = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1])
            sage: RCE = crystals.RiggedConfigurations(LaE[1])
            sage: RC is RCE
            False
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded

        if wt is None:
            wt = cartan_type
            cartan_type = wt.parent().cartan_type()
        else:
            if not isinstance(cartan_type, CartanTypeFolded):
                cartan_type = CartanType(cartan_type)

        if WLR is None:
            WLR = wt.parent()
        else:
            wt = WLR(wt)

        if isinstance(cartan_type, CartanTypeFolded):
            return CrystalOfNonSimplyLacedRC(cartan_type, wt, WLR)

        return super().__classcall__(cls, wt, WLR=WLR)

    def __init__(self, wt, WLR):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: La = RootSystem(['A', 2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1] + La[2])
            sage: TestSuite(RC).run()

            sage: La = RootSystem(['A', 2, 1]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[0])
            sage: TestSuite(RC).run() # long time
        """
        self._cartan_type = WLR.cartan_type()
        self._wt = wt
        self._rc_index = self._cartan_type.index_set()
        self._rc_index_inverse = {i: ii for ii,i in enumerate(self._rc_index)}
        # We store the Cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = self._cartan_type.cartan_matrix()
        if self._cartan_type.is_finite():
            category = ClassicalCrystals()
        else:
            category = (RegularCrystals(), HighestWeightCrystals(), InfiniteEnumeratedSets())
        Parent.__init__(self, category=category)
        n = self._cartan_type.rank()  # == len(self._cartan_type.index_set())
        self.module_generators = (self.element_class(self, partition_list=[[] for _ in repeat(None, n)]),)

    options = RiggedConfigurations.options

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['A', 3]).weight_lattice().fundamental_weights()
            sage: crystals.RiggedConfigurations(La[1])
            Crystal of rigged configurations of type ['A', 3] and weight Lambda[1]
        """
        return "Crystal of rigged configurations of type {0} and weight {1}".format(
                self._cartan_type, self._wt)

    def _element_constructor_(self, *lst, **options):
        """
        Construct a ``RiggedConfigurationElement``.

        Typically the user should not call this method since it does not check
        if it is an actual configuration in the crystal. Instead the user
        should use the iterator.

        EXAMPLES::

            sage: La = RootSystem(['A', 2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1] + La[2])
            sage: RC(partition_list=[[1],[1]], rigging_list=[[0],[-1]])
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]-1
            <BLANKLINE>
            sage: RC(partition_list=[[1],[2]])
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -2[ ][ ]-2
            <BLANKLINE>

        TESTS:

        Check that :issue:`17054` is fixed::

            sage: La = RootSystem(['A', 2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(4*La[1] + 4*La[2])
            sage: B = crystals.infinity.RiggedConfigurations(['A',2])
            sage: x = B.an_element().f_string([2,2,1,1,2,1,2,1])
            sage: ascii_art(x)
            -4[ ][ ][ ][ ]-4  -4[ ][ ][ ][ ]0
            sage: ascii_art(RC(x.nu()))
            0[ ][ ][ ][ ]-4  0[ ][ ][ ][ ]0
            sage: x == B.an_element().f_string([2,2,1,1,2,1,2,1])
            True
        """
        if isinstance(lst[0], (list, tuple)):
            lst = lst[0]

        if isinstance(lst[0], RiggedPartition):
            lst = [p._clone() for p in lst] # Make a deep copy
        elif isinstance(lst[0], RiggedConfigurationElement):
            lst = [p._clone() for p in lst[0]] # Make a deep copy

        return self.element_class(self, list(lst), **options)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number `p_i^{(a)}(\nu)` in ``self``.

        This assumes that `\gamma_a = 1` for all `a` and
        `(\alpha_a | \alpha_b ) = A_{ab}`.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row length

        TESTS::

            sage: La = RootSystem(['A', 2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1] + La[2])
            sage: elt = RC(partition_list=[[1],[2]])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 2)
            -2
        """
        vac_num = self._wt[self.index_set()[a]]

        for b,nu in enumerate(partitions):
            val = self._cartan_matrix[a,b]
            if val:
                if i == float('inf'):
                    vac_num -= val * sum(nu)
                else:
                    vac_num -= val * nu.get_num_cells_to_column(i)

        return vac_num

    def weight_lattice_realization(self):
        """
        Return the weight lattice realization used to express the weights
        of elements in ``self``.

        EXAMPLES::

            sage: La = RootSystem(['A', 2, 1]).weight_lattice(extended=True).fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[0])
            sage: RC.weight_lattice_realization()
            Extended weight lattice of the Root system of type ['A', 2, 1]
        """
        return self._wt.parent()

    Element = RCHighestWeightElement


class CrystalOfNonSimplyLacedRC(CrystalOfRiggedConfigurations):
    """
    Highest weight crystal of rigged configurations in non-simply-laced type.
    """

    def __init__(self, vct, wt, WLR):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: La = RootSystem(['C', 3]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1])
            sage: TestSuite(RC).run()
        """
        self._folded_ct = vct
        CrystalOfRiggedConfigurations.__init__(self, wt, WLR)

    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: La = RootSystem(['C', 2, 1]).weight_lattice().fundamental_weights()
            sage: vct = CartanType(['C', 2, 1]).as_folding()
            sage: RC = crystals.RiggedConfigurations(vct, La[0])
            sage: RC
            Crystal of rigged configurations of type ['C', 2, 1] and weight Lambda[0]
            sage: RC.virtual
            Crystal of rigged configurations of type ['A', 3, 1] and weight 2*Lambda[0]
        """
        P = self._folded_ct._folding.root_system().weight_lattice()
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        vwt = P.sum_of_terms((b, gamma[a]*c) for a,c in self._wt for b in sigma[a])
        return CrystalOfRiggedConfigurations(vwt)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number `p_i^{(a)}(\nu)` in ``self``.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row length

        TESTS::

            sage: La = RootSystem(['C', 3]).weight_lattice().fundamental_weights()
            sage: vct = CartanType(['C', 3]).as_folding()
            sage: RC = crystals.RiggedConfigurations(vct, La[2])
            sage: elt = RC(partition_list=[[], [1], [1]])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 1)
            0
            sage: RC._calc_vacancy_number(elt.nu(), 2, 1)
            -1
        """
        I = self.index_set()
        ia = I[a]
        vac_num = self._wt[ia]

        if i == float('inf'):
            return vac_num - sum(self._cartan_matrix[a,b] * sum(nu)
                                 for b,nu in enumerate(partitions))

        gamma = self._folded_ct.scaling_factors()
        g = gamma[ia]
        for b, nu in enumerate(partitions):
            ib = I[b]
            q = nu.get_num_cells_to_column(g*i, gamma[ib])
            vac_num -= self._cartan_matrix[a,b] * q / gamma[ib]

        return vac_num

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- a rigged configuration element

        EXAMPLES::

            sage: La = RootSystem(['C', 3]).weight_lattice().fundamental_weights()
            sage: vct = CartanType(['C', 3]).as_folding()
            sage: RC = crystals.RiggedConfigurations(vct, La[2])
            sage: elt = RC(partition_list=[[], [1], [1]]); elt
            <BLANKLINE>
            (/)
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: RC.to_virtual(elt)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -2[ ][ ]-2
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        gamma = [int(f) for f in self._folded_ct.scaling_factors()]
        sigma = self._folded_ct._orbit
        n = self._folded_ct._folding.rank()
        vindex = self._folded_ct._folding.index_set()
        partitions = [None] * n
        riggings = [None] * n
        for a, rp in enumerate(rc):
            for i in sigma[a]:
                k = vindex.index(i)
                partitions[k] = [row_len*gamma[a] for row_len in rp._list]
                riggings[k] = [rig_val*gamma[a] for rig_val in rp.rigging]
        return self.virtual.element_class(self.virtual, partition_list=partitions,
                                          rigging_list=riggings)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configuration of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- a virtual rigged configuration

        EXAMPLES::

            sage: La = RootSystem(['C', 3]).weight_lattice().fundamental_weights()
            sage: vct = CartanType(['C', 3]).as_folding()
            sage: RC = crystals.RiggedConfigurations(vct, La[2])
            sage: elt = RC(partition_list=[[0], [1], [1]])
            sage: elt == RC.from_virtual(RC.to_virtual(elt))
            True
        """
        gamma = list(self._folded_ct.scaling_factors())  # map(int, self._folded_ct.scaling_factors())
        sigma = self._folded_ct._orbit
        n = self._cartan_type.rank()
        partitions = [None] * n
        riggings = [None] * n
        vindex = self._folded_ct._folding.index_set()
        for a in range(n):
            index = vindex.index(sigma[a][0])
            partitions[a] = [row_len // gamma[a] for row_len in vrc[index]._list]
            riggings[a] = [rig_val / gamma[a] for rig_val in vrc[index].rigging]
        return self.element_class(self, partition_list=partitions, rigging_list=riggings)

    Element = RCHWNonSimplyLacedElement
