r"""
Bijection classes for type `D_{n+1}^{(2)}`

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `D_{n+1}^{(2)}`.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 2], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import KRTToRCBijectionTypeDTwisted
    sage: bijection = KRTToRCBijectionTypeDTwisted(KRT(pathlist=[[-1,2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['D', 4, 2], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import RCToKRTBijectionTypeDTwisted
    sage: bijection = RCToKRTBijectionTypeDTwisted(RC())
    sage: TestSuite(bijection).run()
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

from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_A2_even import KRTToRCBijectionTypeA2Even
from sage.combinat.rigged_configurations.bij_type_A2_even import RCToKRTBijectionTypeA2Even
from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD

class KRTToRCBijectionTypeDTwisted(KRTToRCBijectionTypeD, KRTToRCBijectionTypeA2Even):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `D_{n+1}^{(2)}`.

    This inherits from type `C_n^{(1)}` and `D_n^{(1)}` because we use the
    same methods in some places.
    """
    def run(self, verbose=False):
        """
        Run the bijection from a tensor product of KR tableaux to a rigged
        configuration for type `D_{n+1}^{(2)}`.

        INPUT:

        - ``tp_krt`` -- A tensor product of KR tableaux

        - ``verbose`` -- (Default: ``False``) Display each step in the
          bijection

        EXAMPLES::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 2], [[3,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import KRTToRCBijectionTypeDTwisted
            sage: KRTToRCBijectionTypeDTwisted(KRT(pathlist=[[-1,3,2]])).run()
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
        """
        if verbose:
            from sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element \
              import TensorProductOfKirillovReshetikhinTableauxElement

        for cur_crystal in reversed(self.tp_krt):
            r = cur_crystal.parent().r()
            # Iterate through the columns
            for col_number, cur_column in enumerate(reversed(cur_crystal.to_array(False))):
                self.cur_path.insert(0, []) # Prepend an empty list

                # Check to see if we are a spinor column
                if r == self.n:
                    if verbose:
                        print("====================")
                        print(repr(TensorProductOfKirillovReshetikhinTableauxElement(self.tp_krt.parent(), self.cur_path)))
                        print("--------------------")
                        print(repr(self.ret_rig_con))
                        print("--------------------\n")
                        print("Applying doubling map")
                    self.doubling_map()

                self.cur_dims.insert(0, [0, 1])

                for letter in reversed(cur_column):
                    self.cur_dims[0][0] += 1
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

                # Check to see if we are a spinor column
                if r == self.n:
                    if verbose:
                        print("====================")
                        print(repr(TensorProductOfKirillovReshetikhinTableauxElement(self.tp_krt.parent(), self.cur_path)))
                        print("--------------------")
                        print(repr(self.ret_rig_con))
                        print("--------------------\n")
                        print("Applying halving map")
                    self.halving_map()

                # If we've split off a column, we need to merge the current column
                #   to the current crystal tableau
                if col_number > 0:
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

    def next_state(self, val):
        r"""
        Build the next state for type `D_{n+1}^{(2)}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 2], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import KRTToRCBijectionTypeDTwisted
            sage: bijection = KRTToRCBijectionTypeDTwisted(KRT(pathlist=[[-1,2]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [2])
            sage: bijection.next_state(2)
        """
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        if val == 'E':
            KRTToRCBijectionTypeA2Even.next_state(self, val)
            return
        elif val > 0:
            # If it is a regular value, we follow the A_n rules
            KRTToRCBijectionTypeA.next_state(self, val)
            if tableau_height >= n - 1:
                self._correct_vacancy_nums()
            return

        pos_val = -val

        if pos_val == 0:
            if len(self.ret_rig_con[pos_val - 1]) > 0:
                max_width = self.ret_rig_con[n-1][0]
            else:
                max_width = 1
            max_width = self.ret_rig_con[n-1].insert_cell(max_width)
            width_n = max_width + 1

            # Follow regular A_n rules
            for a in reversed(range(tableau_height, n-1)):
                max_width = self.ret_rig_con[a].insert_cell(max_width)
                self._update_vacancy_nums(a + 1)
                self._update_partition_values(a + 1)
            self._update_vacancy_nums(tableau_height)
            if tableau_height >= n - 1:
                self._correct_vacancy_nums()
            self._update_partition_values(tableau_height)
            if tableau_height > 0:
                self._update_vacancy_nums(tableau_height-1)
                self._update_partition_values(tableau_height-1)

            # Make the new string at n quasi-singular
            p = self.ret_rig_con[n-1]
            num_rows = len(p)
            for i in range(num_rows):
                if p._list[i] == width_n:
                    j = i+1
                    while j < num_rows and p._list[j] == width_n \
                      and p.vacancy_numbers[j] == p.rigging[j]:
                        j += 1
                    p.rigging[j-1] -= 1
                    break
            return

        case_S = [None] * n
        pos_val = -val

        # Always add a cell to the first singular value in the first
        #   tableau we are updating.
        if len(self.ret_rig_con[pos_val - 1]) > 0:
            max_width = self.ret_rig_con[pos_val - 1][0]
        else:
            max_width = 1

        # Add cells similar to type A_n but we move to the right until we
        #   reach the value of n-1
        for a in range(pos_val - 1, n-1):
            max_width = self.ret_rig_con[a].insert_cell(max_width)
            case_S[a] = max_width

        # Special case for n
        # If we find a quasi-singular string first, then we are in case (Q, S)
        #   otherwise we will find a singular string and insert 2 cells
        partition = self.ret_rig_con[n-1]
        num_rows = len(partition)
        case_QS = False
        for i in range(num_rows + 1):
            if i == num_rows:
                max_width = 0
                if case_QS:
                    partition._list.append(1)
                    partition.vacancy_numbers.append(None)
                    # Go through our partition until we find a length of greater than 1
                    j = len(partition._list) - 1
                    while j >= 0 and partition._list[j] == 1:
                        j -= 1
                    partition.rigging.insert(j + 1, None)
                    width_n = 1
                else:
                    # Go through our partition until we find a length of greater than 2
                    j = len(partition._list) - 1
                    while j >= 0 and partition._list[j] <= 2:
                        j -= 1
                    partition._list.insert(j+1, 2)
                    partition.vacancy_numbers.insert(j+1, None)
                    partition.rigging.insert(j+1, None)
                break
            elif partition._list[i] <= max_width:
                if partition.vacancy_numbers[i] == partition.rigging[i]:
                    max_width = partition._list[i]
                    if case_QS:
                        partition._list[i] += 1
                        width_n = partition._list[i]
                        partition.rigging[i] = None
                    else:
                        j = i - 1
                        while j >= 0 and partition._list[j] <= max_width + 2:
                            partition.rigging[j+1] = partition.rigging[j] # Shuffle it along
                            j -= 1
                        partition._list.pop(i)
                        partition._list.insert(j+1, max_width + 2)
                        partition.rigging[j+1] = None
                    break
                elif partition.vacancy_numbers[i] - 1 == partition.rigging[i] and not case_QS:
                    case_QS = True
                    partition._list[i] += 1
                    partition.rigging[i] = None
                    # No need to set max_width here since we will find a singular string

        # Now go back following the regular C_n (ish) rules
        for a in reversed(range(tableau_height, n-1)):
            if case_S[a] == max_width:
                self._insert_cell_case_S(self.ret_rig_con[a])
            else:
                max_width = self.ret_rig_con[a].insert_cell(max_width)
            self._update_vacancy_nums(a + 1)
            self._update_partition_values(a + 1)

        # Update the final rigged partitions
        self._update_vacancy_nums(tableau_height)
        if tableau_height >= n - 1:
            self._correct_vacancy_nums()
        self._update_partition_values(tableau_height)

        if pos_val <= tableau_height:
            for a in range(pos_val-1, tableau_height):
                self._update_vacancy_nums(a)
                self._update_partition_values(a)
            if pos_val > 1:
                self._update_vacancy_nums(pos_val - 2)
                self._update_partition_values(pos_val - 2)
        elif tableau_height > 0:
            self._update_vacancy_nums(tableau_height - 1)
            self._update_partition_values(tableau_height - 1)

        if case_QS:
            # Make the new string quasi-singular
            num_rows = len(partition)
            for i in range(num_rows):
                if partition._list[i] == width_n:
                    j = i+1
                    while j < num_rows and partition._list[j] == width_n \
                      and partition.vacancy_numbers[j] == partition.rigging[j]:
                        j += 1
                    partition.rigging[j-1] -= 1
                    break

class RCToKRTBijectionTypeDTwisted(RCToKRTBijectionTypeD, RCToKRTBijectionTypeA2Even):
    r"""
    Specific implementation of the bijection from rigged configurations to
    tensor products of KR tableaux for type `D_{n+1}^{(2)}`.
    """
    def run(self, verbose=False, build_graph=False):
        """
        Run the bijection from rigged configurations to tensor product of KR
        tableaux for type `D_{n+1}^{(2)}`.

        INPUT:

        - ``verbose`` -- (default: ``False``) display each step in the
          bijection
        - ``build_graph`` -- (default: ``False``) build the graph of each
          step of the bijection

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 2], [[3, 1]])
            sage: x = RC(partition_list=[[],[1],[1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import RCToKRTBijectionTypeDTwisted
            sage: RCToKRTBijectionTypeDTwisted(x).run()
            [[1], [3], [-2]]
            sage: bij = RCToKRTBijectionTypeDTwisted(x)
            sage: bij.run(build_graph=True)
            [[1], [3], [-2]]
            sage: bij._graph
            Digraph on 6 vertices
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
                    self.cur_dims[0][1] -= 1
                    self.cur_dims.insert(0, [dim[0], 1])

                    # Perform the corresponding splitting map on rigged configurations
                    # All it does is update the vacancy numbers on the RC side
                    for a in range(self.n):
                        self._update_vacancy_numbers(a)

                    if build_graph:
                        y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                        self._graph.append([self._graph[-1][1], (y, len(self._graph)), 'ls'])

                # Check to see if we are a spinor
                if dim[0] == self.n:
                    if verbose:
                        print("====================")
                        print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
                        print("--------------------")
                        print(ret_crystal_path)
                        print("--------------------\n")
                        print("Applying doubling map")
                    self.doubling_map()

                    if build_graph:
                        y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                        self._graph.append([self._graph[-1][1], (y, len(self._graph)), '2x'])

                while self.cur_dims[0][0] > 0:
                    if verbose:
                        print("====================")
                        print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
                        print("--------------------")
                        print(ret_crystal_path)
                        print("--------------------\n")

                    self.cur_dims[0][0] -= 1 # This takes care of the indexing
                    b = self.next_state(self.cur_dims[0][0])

                    # Make sure we have a crystal letter
                    ret_crystal_path[-1].append(letters(b)) # Append the rank

                    if build_graph:
                        y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                        self._graph.append([self._graph[-1][1], (y, len(self._graph)), letters(b)])

                self.cur_dims.pop(0) # Pop off the leading column

                # Check to see if we were a spinor
                if dim[0] == self.n:
                    if verbose:
                        print("====================")
                        print(repr(self.rigged_con.parent()(*self.cur_partitions)))
                        print("--------------------")
                        print(ret_crystal_path)
                        print("--------------------\n")
                        print("Applying halving map")
                    self.halving_map()

                    if build_graph:
                        y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                        self._graph.append([self._graph[-1][1], (y, len(self._graph)), '1/2x'])

        if build_graph:
            self._graph.pop(0) # Remove the dummy at the start
            from sage.graphs.digraph import DiGraph
            from sage.graphs.dot2tex_utils import have_dot2tex
            self._graph = DiGraph(self._graph)
            if have_dot2tex():
                self._graph.set_latex_options(format="dot2tex", edge_labels=True)

        return self.KRT(pathlist=ret_crystal_path)

    def next_state(self, height):
        r"""
        Build the next state for type `D_{n+1}^{(2)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['D', 4, 2], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import RCToKRTBijectionTypeDTwisted
            sage: bijection = RCToKRTBijectionTypeDTwisted(RC(partition_list=[[2],[2,2],[2,2]]))
            sage: bijection.next_state(0)
            -1
        """
        n = self.n
        ell = [None] * (2*n)
        case_S = [False] * n
        case_Q = False
        b = None

        # Calculate the rank and ell values

        last_size = 0
        for a in range(height, n-1):
            ell[a] = self._find_singular_string(self.cur_partitions[a], last_size)

            if ell[a] is None:
                b = a + 1
                break
            else:
                last_size = self.cur_partitions[a][ell[a]]

        if b is None:
            partition = self.cur_partitions[n-1]
            # Modified version of _find_singular_string()
            for i in reversed(range(len(partition))):
                if partition[i] >= last_size:
                    if partition.vacancy_numbers[i] == partition.rigging[i]:
                        if partition[i] == 1:
                            b = 'E'
                        else:
                            last_size = partition[i]
                        case_S[n-1] = True
                        ell[2*n-1] = i
                        break
                    elif partition.vacancy_numbers[i] - 1 == partition.rigging[i] and not case_Q:
                        case_Q = True
                        # Check if it is singular as well
                        block_size = partition[i]
                        for j in reversed(range(i)):
                            if partition[j] != block_size:
                                break
                            elif partition.vacancy_numbers[j] == partition.rigging[j]:
                                case_Q = False
                                break
                        if case_Q:
                            last_size = partition[i] + 1
                            ell[n-1] = i

            if ell[2*n-1] is None:
                if not case_Q:
                    b = n
                else:
                    b = 0

        if b is None:
            # Now go back
            for a in reversed(range(n-1)):
                if a >= height and self.cur_partitions[a][ell[a]] == last_size:
                    ell[n+a] = ell[a]
                    case_S[a] = True
                else:
                    ell[n+a] = self._find_singular_string(self.cur_partitions[a], last_size)

                    if ell[n + a] is None:
                        b = -(a + 2)
                        break
                    else:
                        last_size = self.cur_partitions[a][ell[n + a]]

        if b is None:
            b = -1

        # Determine the new rigged configuration by removing boxes from the
        #   selected string and then making the new string singular
        if case_S[0]:
            row_num = None
            row_num_bar = self.cur_partitions[0].remove_cell(ell[n], 2)
        else:
            row_num = self.cur_partitions[0].remove_cell(ell[0])
            row_num_bar = self.cur_partitions[0].remove_cell(ell[n])
        for a in range(1, n-1):
            if case_S[a]:
                row_num_next = None
                row_num_bar_next = self.cur_partitions[a].remove_cell(ell[n+a], 2)
            else:
                row_num_next = self.cur_partitions[a].remove_cell(ell[a])
                row_num_bar_next = self.cur_partitions[a].remove_cell(ell[n+a])

            self._update_vacancy_numbers(a - 1)
            if row_num is not None:
                self.cur_partitions[a-1].rigging[row_num] = self.cur_partitions[a-1].vacancy_numbers[row_num]
            if row_num_bar is not None:
                self.cur_partitions[a-1].rigging[row_num_bar] = self.cur_partitions[a-1].vacancy_numbers[row_num_bar]
            row_num = row_num_next
            row_num_bar = row_num_bar_next

        if case_Q:
            row_num_next = self.cur_partitions[n-1].remove_cell(ell[n-1])
            if case_S[n-1]:
                row_num_bar_next = self.cur_partitions[n-1].remove_cell(ell[2*n-1])
            else:
                row_num_bar_next = None
        elif case_S[n-1]:
            row_num_next = None
            row_num_bar_next = self.cur_partitions[n-1].remove_cell(ell[2*n-1], 2)
        else:
            row_num_next = None
            row_num_bar_next = None

        self._update_vacancy_numbers(n - 2)
        if row_num is not None:
            self.cur_partitions[n-2].rigging[row_num] = self.cur_partitions[n-2].vacancy_numbers[row_num]
        if row_num_bar is not None:
            self.cur_partitions[n-2].rigging[row_num_bar] = self.cur_partitions[n-2].vacancy_numbers[row_num_bar]

        self._update_vacancy_numbers(n - 1)
        if height == n:
            self._correct_vacancy_nums()
        if row_num_next is not None:
            self.cur_partitions[n-1].rigging[row_num_next] = self.cur_partitions[n-1].vacancy_numbers[row_num_next]
        if row_num_bar_next is not None:
            if case_Q:
                vac_num = self.cur_partitions[n-1].vacancy_numbers[row_num_bar_next]
                self.cur_partitions[n-1].rigging[row_num_bar_next] = vac_num
                block_len = self.cur_partitions[n-1][row_num_bar_next]
                j = row_num_bar_next + 1
                length = len(self.cur_partitions[n-1])
                # Find the place for the quasisingular rigging
                while j < length and self.cur_partitions[n-1][j] == block_len \
                  and self.cur_partitions[n-1].rigging[j] == vac_num:
                    j += 1
                self.cur_partitions[n-1].rigging[j-1] = vac_num - 1
            else:
                self.cur_partitions[n-1].rigging[row_num_bar_next] = self.cur_partitions[n-1].vacancy_numbers[row_num_bar_next]

        return(b)

