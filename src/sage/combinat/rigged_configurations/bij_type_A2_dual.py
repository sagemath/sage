# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Bijection classes for type `A_{2n}^{(2)\dagger}`

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `A_{2n}^{(2)\dagger}`.

AUTHORS:

- Travis Scrimshaw (2012-12-21): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(CartanType(['A', 4, 2]).dual(), [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_A2_dual import KRTToRCBijectionTypeA2Dual
    sage: bijection = KRTToRCBijectionTypeA2Dual(KRT(pathlist=[[2,1]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(CartanType(['A', 4, 2]).dual(), [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_A2_dual import RCToKRTBijectionTypeA2Dual
    sage: bijection = RCToKRTBijectionTypeA2Dual(RC(partition_list=[[],[]]))
    sage: TestSuite(bijection).run()
"""

# ****************************************************************************
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from sage.combinat.rigged_configurations.bij_type_C import KRTToRCBijectionTypeC
from sage.combinat.rigged_configurations.bij_type_C import RCToKRTBijectionTypeC
from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA

from sage.rings.rational_field import QQ


class KRTToRCBijectionTypeA2Dual(KRTToRCBijectionTypeC):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `A_{2n}^{(2)\dagger}`.

    This inherits from type `C_n^{(1)}` because we use the same methods in
    some places.
    """

    def next_state(self, val):
        r"""
        Build the next state for type `A_{2n}^{(2)\dagger}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(CartanType(['A', 4, 2]).dual(), [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A2_dual import KRTToRCBijectionTypeA2Dual
            sage: bijection = KRTToRCBijectionTypeA2Dual(KRT(pathlist=[[-1,2]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [2])
            sage: bijection.next_state(2)
        """
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        if val > 0:
            # If it is a regular value, we follow the A_n rules
            KRTToRCBijectionTypeA.next_state(self, val)
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
            self._update_partition_values(tableau_height)
            if tableau_height > 0:
                self._update_vacancy_nums(tableau_height-1)
                self._update_partition_values(tableau_height-1)

            # Make the new string at n quasi-singular
            p = self.ret_rig_con[n-1]
            for i in range(len(p)):
                if p._list[i] == width_n:
                    p.rigging[i] = p.rigging[i] - QQ(1)/QQ(2)
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
                            partition.rigging[j+1] = partition.rigging[j]  # Shuffle it along
                            j -= 1
                        partition._list.pop(i)
                        partition._list.insert(j+1, max_width + 2)
                        partition.rigging[j+1] = None
                    break
                elif partition.vacancy_numbers[i] - QQ(1)/QQ(2) == partition.rigging[i] and not case_QS:
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
        if tableau_height < n:
            self._update_vacancy_nums(tableau_height)
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
                    partition.rigging[i] = partition.rigging[i] - QQ(1)/QQ(2)
                    break


class RCToKRTBijectionTypeA2Dual(RCToKRTBijectionTypeC):
    r"""
    Specific implementation of the bijection from rigged configurations to
    tensor products of KR tableaux for type `A_{2n}^{(2)\dagger}`.
    """

    def next_state(self, height):
        r"""
        Build the next state for type `A_{2n}^{(2)\dagger}`.

        TESTS::

            sage: RC = RiggedConfigurations(CartanType(['A', 4, 2]).dual(), [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A2_dual import RCToKRTBijectionTypeA2Dual
            sage: bijection = RCToKRTBijectionTypeA2Dual(RC(partition_list=[[2],[2,2]]))
            sage: bijection.next_state(2)
            -1
        """
        height -= 1  # indexing
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
            # Special case for n
            for i in reversed(range(len(partition))):
                if partition[i] >= last_size:
                    if partition.vacancy_numbers[i] == partition.rigging[i]:
                        last_size = partition[i]
                        case_S[n-1] = True
                        ell[2*n-1] = i
                        break
                    elif partition.vacancy_numbers[i] - QQ(1)/QQ(2) == partition.rigging[i] and not case_Q:
                        case_Q = True
                        # This will never be singular
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
        if n > 1:
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

        if n > 1:
            self._update_vacancy_numbers(n - 2)
            if row_num is not None:
                self.cur_partitions[n-2].rigging[row_num] = self.cur_partitions[n-2].vacancy_numbers[row_num]
            if row_num_bar is not None:
                self.cur_partitions[n-2].rigging[row_num_bar] = self.cur_partitions[n-2].vacancy_numbers[row_num_bar]

        self._update_vacancy_numbers(n - 1)
        if row_num_next is not None:
            self.cur_partitions[n-1].rigging[row_num_next] = self.cur_partitions[n-1].vacancy_numbers[row_num_next]
        if row_num_bar_next is not None:
            if case_Q:
                # This will always be the largest value
                self.cur_partitions[n-1].rigging[row_num_bar_next] = self.cur_partitions[n-1].vacancy_numbers[row_num_bar_next] - QQ(1)/QQ(2)
            else:
                self.cur_partitions[n-1].rigging[row_num_bar_next] = self.cur_partitions[n-1].vacancy_numbers[row_num_bar_next]

        return b
