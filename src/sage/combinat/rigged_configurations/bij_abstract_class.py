r"""
Abstract classes for the rigged configuration bijections

This file contains two sets of classes, one for the bijection from KR tableaux to
rigged configurations and the other for the reverse bijection. We do this for
two reasons, one is because we can store a state in the bijection locally, so
we do not have to constantly pass it around between functions. The other is because
it makes the code easier to read in the \*_element.py files.

These classes are not meant to be used by the user and are only supposed to be
used internally to perform the bijections between
:class:`~sage.combinat.rigged_configurations.tensor_product_kr_tableaux.TensorProductOfKirillovReshetikhinTableaux`
and :class:`RiggedConfigurations`.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2011, 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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
#*****************************************************************************

from copy import deepcopy
from sage.misc.abstract_method import abstract_method

class KRTToRCBijectionAbstract:
    """
    Root abstract class for the bijection from KR tableaux to rigged configurations.

    This class holds the state of the bijection and generates the next state.
    This class should never be created directly.
    """

    def __init__(self, tp_krt):
        """
        Initialize the bijection by obtaining the important information from
        the KR tableaux.

        INPUT:

        - ``parent`` -- The parent of tensor product of KR tableaux

        EXAMPLES::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
            sage: bijection = KRTToRCBijectionTypeA(KRT(pathlist=[[3,1]]))
            sage: TestSuite(bijection).run()
        """
        self.tp_krt = tp_krt
        self.n = tp_krt.parent().cartan_type().classical().rank()
        self.ret_rig_con = tp_krt.parent().rigged_configurations()(partition_list=[[]] * self.n)
        # We allow this to be mutable to make the bijection easier to program.
        # Upon completing the bijection, this will be set to immutable.
        # Do not call this, the object could be in a mutable state and ultimately
        #   be placed in an unstable state.
        # The user will (and should) never know about this temporary mutable state.
        self.ret_rig_con._set_mutable()
        self.cur_dims = []
        self.cur_path = []
        # self.L = {}

    def __eq__(self, rhs):
        r"""
        Check equality.

        This is only here for pickling check. This is a temporary placeholder
        class, and as such, should never be compared.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
            sage: bijection = KRTToRCBijectionTypeA(KRT(pathlist=[[5,3]]))
            sage: bijection2 = KRTToRCBijectionTypeA(KRT(pathlist=[[5,3]]))
            sage: bijection == bijection2
            True
        """
        return isinstance(rhs, KRTToRCBijectionAbstract)

    def run(self, verbose=False):
        """
        Run the bijection from a tensor product of KR tableaux to a rigged
        configuration.

        INPUT:

        - ``tp_krt`` -- A tensor product of KR tableaux

        - ``verbose`` -- (Default: ``False``) Display each step in the
          bijection

        EXAMPLES::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
            sage: KRTToRCBijectionTypeA(KRT(pathlist=[[5,2]])).run()
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
        """
        if verbose:
            from sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element \
              import TensorProductOfKirillovReshetikhinTableauxElement

        for cur_crystal in reversed(self.tp_krt):
            target = cur_crystal.parent()._r
            # Iterate through the columns
            for col_number, cur_column in enumerate(reversed(cur_crystal.to_array(False))):
                self.cur_path.insert(0, []) # Prepend an empty list

                self.cur_dims.insert(0, [0, 1])

                for letter in reversed(cur_column):
                    self.cur_dims[0][0] = self._next_index(self.cur_dims[0][0], target)
                    val = letter.value # Convert from a CrystalOfLetter to an Integer

                    if verbose:
                        print("====================")
                        print(repr(TensorProductOfKirillovReshetikhinTableauxElement(self.tp_krt.parent(), self.cur_path)))
                        print("--------------------")
                        print(repr(self.ret_rig_con))
                        print("--------------------\n")

                    # Build the next state
                    self.cur_path[0].insert(0, [letter]) # Prepend the value
                    self.next_state(val)

                # If we've split off a column, we need to merge the current column
                #   to the current crystal tableau
                if col_number > 0:
                    if verbose:
                        print("====================")
                        print(repr(TensorProductOfKirillovReshetikhinTableauxElement(self.tp_krt.parent(), self.cur_path)))
                        print("--------------------")
                        print(repr(self.ret_rig_con))
                        print("--------------------\n")
                        print("Applying column merge")

                    for i, letter_singleton in enumerate(self.cur_path[0]):
                        self.cur_path[1][i].insert(0, letter_singleton[0])
                    self.cur_dims[1][1] += 1
                    self.cur_path.pop(0)
                    self.cur_dims.pop(0)

                    # And perform the inverse column splitting map on the RC
                    for a in range(self.n):
                        self._update_vacancy_nums(a)

        self.ret_rig_con.set_immutable() # Return it to immutable
        return self.ret_rig_con

    @abstract_method
    def next_state(self, val):
        r"""
        Build the next state in the bijection.

        INPUT:

        - ``val`` -- The value we are adding

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
            sage: bijection = KRTToRCBijectionTypeA(KRT(pathlist=[[5,3]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [3])
            sage: bijection.next_state(3)
            sage: bijection.ret_rig_con
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """

    def _update_vacancy_nums(self, a):
        r"""
        Update the vacancy numbers of a rigged partition.

        Helper function to (batch) update the vacancy numbers of the rigged
        partition at position `a` in the rigged configuration stored by this
        bijection.

        INPUT:

        - ``a`` -- The index of the partition to update

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
            sage: bijection = KRTToRCBijectionAbstract(KRT(pathlist=[[3,2]]))  
            sage: bijection._update_vacancy_nums(2)
        """
        # Check to make sure we have a valid index (currently removed)
        # If the current tableau is empty, there is nothing to do
        if not self.ret_rig_con[a]: # Check to see if we have vacancy numbers
            return

        # Setup the first block
        block_len = self.ret_rig_con[a][0]
        nu = self.ret_rig_con.nu()
        vac_num = self.ret_rig_con.parent()._calc_vacancy_number(nu, a, nu[a][0],
                                                                 dims=self.cur_dims)

        for i, row_len in enumerate(self.ret_rig_con[a]):
            # If we've gone to a different sized block, then update the
            #   values which change when moving to a new block size
            if block_len != row_len:
                vac_num = self.ret_rig_con.parent()._calc_vacancy_number(nu, a, row_len,
                                                                         dims=self.cur_dims)
                block_len = row_len
            self.ret_rig_con[a].vacancy_numbers[i] = vac_num

    def _update_partition_values(self, a):
        r"""
        Update the partition values of a rigged partition.

        Helper function to update the partition values of a given rigged
        partition row. This will go through all of our partition values and set
        them to our vacancy number if the corresponding row has been changed
        (indicated by being set to ``None``).

        INPUT:

        - ``a`` -- The index of the partition to update

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
            sage: bijection = KRTToRCBijectionAbstract(KRT(pathlist=[[5,2]]))
            sage: bijection._update_partition_values(2)
        """
        rigged_partition = self.ret_rig_con[a]
        for index, value in enumerate(rigged_partition.rigging):
            if value is None:
                rigged_partition.rigging[index] = rigged_partition.vacancy_numbers[index]
                if index > 0 and rigged_partition[index - 1] == rigged_partition[index] \
                  and rigged_partition.rigging[index - 1] < rigged_partition.rigging[index]:
                    # If we need to reorder
                    pos = 0
                    width = rigged_partition[index]
                    val = rigged_partition.rigging[index]
                    for i in reversed(range(index-1)):
                        if rigged_partition[i] > width or rigged_partition.rigging[i] >= val:
                            pos = i + 1
                            break

                    rigged_partition.rigging.pop(index)
                    rigged_partition.rigging.insert(pos, val)

    def _next_index(self, r, target):
        """
        Return the next index after ``r`` when performing a step
        in the bijection going towards ``target``.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
            sage: bijection = KRTToRCBijectionAbstract(KRT(pathlist=[[5,2]]))
            sage: bijection._next_index(1, 2)
            2
        """
        return r + 1

class RCToKRTBijectionAbstract:
    """
    Root abstract class for the bijection from rigged configurations to
    tensor product of Kirillov-Reshetikhin tableaux.

    This class holds the state of the bijection and generates the next state.
    This class should never be created directly.
    """

    def __init__(self, RC_element):
        """
        Initialize the bijection helper.

        INPUT:

        - ``RC_element`` -- The rigged configuration

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract
            sage: bijection = RCToKRTBijectionAbstract(RC(partition_list=[[1],[1],[1],[1]]))
            sage: TestSuite(bijection).run()
        """
        # Make a mutable clone of the rigged configuration for the bijection
        # This will be deleted when the bijection is completed
        self.rigged_con = RC_element.__copy__()
        self.n = RC_element.parent().cartan_type().classical().rank()
        self.KRT = RC_element.parent().tensor_product_of_kirillov_reshetikhin_tableaux()

        # Make a (deep) copy of the dimensions for the bijection
        self.cur_dims = [list(x[:]) for x in self.rigged_con.parent().dims]

        # Note that this implementation of the bijection is destructive to cur_partitions,
        #   therefore we will make a (deep) copy of the partitions.
        # TODO: Convert from cur_partitions to rigged_con
        self.cur_partitions = deepcopy(list(self.rigged_con)[:])

        # This is a dummy edge to start the process
        cp = RC_element.__copy__()
        cp.set_immutable()
        self._graph = [ [[], (cp, 0)] ]

    def __eq__(self, rhs):
        r"""
        Check equality.

        This is only here for pickling check. This is a temporary placeholder
        class, and as such, should never be compared.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
            sage: bijection = RCToKRTBijectionTypeA(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection2 = RCToKRTBijectionTypeA(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection == bijection2
            True
        """
        return isinstance(rhs, RCToKRTBijectionAbstract)

    def run(self, verbose=False, build_graph=False):
        """
        Run the bijection from rigged configurations to tensor product of KR
        tableaux.

        INPUT:

        - ``verbose`` -- (default: ``False``) display each step in the
          bijection
        - ``build_graph`` -- (default: ``False``) build the graph of each
          step of the bijection

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: x = RC(partition_list=[[1],[1],[1],[1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
            sage: RCToKRTBijectionTypeA(x).run()
            [[2], [5]]
            sage: bij = RCToKRTBijectionTypeA(x)
            sage: bij.run(build_graph=True)
            [[2], [5]]
            sage: bij._graph
            Digraph on 3 vertices
        """
        from sage.combinat.crystals.letters import CrystalOfLetters
        letters = CrystalOfLetters(self.rigged_con.parent()._cartan_type.classical())

        # This is technically bad, but because the first thing we do is append
        #   an empty list to ret_crystal_path, we correct this. We do it this
        #   way so that we do not have to remove an empty list after the
        #   bijection has been performed.
        ret_crystal_path = []

        for dim in self.rigged_con.parent().dims:
            ret_crystal_path.append([])

            # Iterate over each column
            for dummy_var in range(dim[1]):
                # Split off a new column if necessary
                if self.cur_dims[0][1] > 1:
                    if verbose:
                        print("====================")
                        print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
                        print("--------------------")
                        print(ret_crystal_path)
                        print("--------------------\n")
                        print("Applying column split")

                    self.cur_dims[0][1] -= 1
                    self.cur_dims.insert(0, [dim[0], 1])

                    # Perform the corresponding splitting map on rigged configurations
                    # All it does is update the vacancy numbers on the RC side
                    for a in range(self.n):
                        self._update_vacancy_numbers(a)

                    if build_graph:
                        y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                        self._graph.append([self._graph[-1][1], (y, len(self._graph)), 'ls'])

                while self.cur_dims[0][0]: # > 0:
                    if verbose:
                        print("====================")
                        print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
                        print("--------------------")
                        print(ret_crystal_path)
                        print("--------------------\n")

                    ht = self.cur_dims[0][0]
                    self.cur_dims[0][0] = self._next_index(ht)
                    b = self.next_state(ht)

                    # Make sure we have a crystal letter
                    ret_crystal_path[-1].append(letters(b)) # Append the rank

                    if build_graph:
                        y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                        self._graph.append([self._graph[-1][1], (y, len(self._graph)), letters(b)])

                self.cur_dims.pop(0) # Pop off the leading column

        if build_graph:
            self._graph.pop(0) # Remove the dummy at the start
            from sage.graphs.digraph import DiGraph
            from sage.graphs.dot2tex_utils import have_dot2tex
            self._graph = DiGraph(self._graph, format="list_of_edges")
            if have_dot2tex():
                self._graph.set_latex_options(format="dot2tex", edge_labels=True)
        return self.KRT(pathlist=ret_crystal_path)

    @abstract_method
    def next_state(self, height):
        """
        Build the next state in the bijection.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
            sage: bijection = RCToKRTBijectionTypeA(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection.next_state(1)
            5
            sage: bijection.cur_partitions
            [(/)
             , (/)
             , (/)
             , (/)
             ]
        """

    def _update_vacancy_numbers(self, a):
        r"""
        Update the vacancy numbers during the bijection.

        INPUT:

        - ``a`` -- The index of the partition to update

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract
            sage: bijection = RCToKRTBijectionAbstract(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection._update_vacancy_numbers(2)
        """

        # Nothing to do if there the rigged partition is empty
        if not self.cur_partitions[a]:
            return

        partition = self.cur_partitions[a]

        # Setup the first block
        block_len = partition[0]
        vac_num = self.rigged_con.parent()._calc_vacancy_number(self.cur_partitions,
                                                                a, partition[0],
                                                                dims=self.cur_dims)

        for i, row_len in enumerate(self.cur_partitions[a]):
            # If we've gone to a different sized block, then update the
            #   values which change when moving to a new block size
            if block_len != row_len:
                vac_num = self.rigged_con.parent()._calc_vacancy_number(self.cur_partitions,
                                                                        a, row_len,
                                                                        dims=self.cur_dims)
                block_len = row_len

            partition.vacancy_numbers[i] = vac_num

    def _find_singular_string(self, partition, last_size):
        r"""
        Return the index of the singular string or ``None`` if not found.

        Helper method to find a singular string at least as long as
        ``last_size``.

        INPUT:

        - ``partition`` -- The partition to look in

        - ``last_size`` -- The last size found

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract
            sage: bijection = RCToKRTBijectionAbstract(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection._find_singular_string(bijection.cur_partitions[2], 2)
            sage: bijection._find_singular_string(bijection.cur_partitions[2], 0)
            0
        """
        for i in reversed(range(len(partition))):
            if (partition[i] >= last_size
                    and partition.vacancy_numbers[i] == partition.rigging[i]):
                return i

    def _next_index(self, r):
        """
        Return the next index after ``r`` when performing a step
        in the bijection.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract
            sage: bijection = RCToKRTBijectionAbstract(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection._next_index(2)
            1
        """
        return r - 1

