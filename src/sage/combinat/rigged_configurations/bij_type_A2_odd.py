# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Bijection classes for type `A_{2n-1}^{(2)}`.

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `A_{2n-1}^{(2)}`.

AUTHORS:

- Travis Scrimshaw (2012-12-21): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 5, 2], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_A2_odd import KRTToRCBijectionTypeA2Odd
    sage: bijection = KRTToRCBijectionTypeA2Odd(KRT(pathlist=[[-1,2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['A', 5, 2], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_A2_odd import RCToKRTBijectionTypeA2Odd
    sage: bijection = RCToKRTBijectionTypeA2Odd(RC(partition_list=[[],[],[]]))
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

from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA


class KRTToRCBijectionTypeA2Odd(KRTToRCBijectionTypeA):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `A_{2n-1}^{(2)}`.

    This inherits from type `A_n^{(1)}` because we use the same methods in
    some places.
    """

    def next_state(self, val):
        r"""
        Build the next state for type `A_{2n-1}^{(2)}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 5, 2], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A2_odd import KRTToRCBijectionTypeA2Odd
            sage: bijection = KRTToRCBijectionTypeA2Odd(KRT(pathlist=[[-2,3]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [3])
            sage: bijection.next_state(3)
        """
        # Note that we must subtract 1 from n to match the indices.
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        # If it is a regular value, we follow the A_n rules
        if val > 0:
            KRTToRCBijectionTypeA.next_state(self, val)
            return

        pos_val = -val

        # Always add a cell to the first singular value in the first
        #   tableau we are updating.
        if len(self.ret_rig_con[pos_val - 1]) > 0:
            max_width = self.ret_rig_con[pos_val - 1][0]
        else:
            max_width = 1

        # Add cells similar to type A_n but we move to the right until we
        #   reach the value of n
        for a in range(pos_val - 1, n):
            max_width = self.ret_rig_con[a].insert_cell(max_width)

        # Now go back following the regular A_n rules
        for a in reversed(range(tableau_height, n - 1)):
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


class RCToKRTBijectionTypeA2Odd(RCToKRTBijectionTypeA):
    r"""
    Specific implementation of the bijection from rigged configurations to
    tensor products of KR tableaux for type `A_{2n-1}^{(2)}`.
    """

    def next_state(self, height):
        r"""
        Build the next state for type `A_{2n-1}^{(2)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 5, 2], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A2_odd import RCToKRTBijectionTypeA2Odd
            sage: bijection = RCToKRTBijectionTypeA2Odd(RC(partition_list=[[1],[2,1],[2]]))
            sage: bijection.next_state(1)
            -2
        """
        height -= 1  # indexing
        n = self.n
        ell = [None] * (2*n)
        b = None

        # Calculate the rank and ell values

        last_size = 0
        for a in range(height, n):
            ell[a] = self._find_singular_string(self.cur_partitions[a], last_size)

            if ell[a] is None:
                b = a + 1
                break
            else:
                last_size = self.cur_partitions[a][ell[a]]

        if b is None:
            # Now go back
            for a in reversed(range(n - 1)):
                # Modified form of _find_singular_string()
                end = ell[a]
                if a < height:
                    end = len(self.cur_partitions[a])
                for i in reversed(range(end)):
                    if self.cur_partitions[a][i] >= last_size and \
                      self.cur_partitions[a].vacancy_numbers[i] == self.cur_partitions[a].rigging[i]:
                        ell[n + a] = i
                        break

                if ell[n + a] is None:
                    b = -(a + 2)
                    break
                else:
                    last_size = self.cur_partitions[a][ell[n + a]]

        if b is None:
            b = -1

        # Determine the new rigged configuration by removing a box from the selected
        #   string and then making the new string singular
        ret_row = self.cur_partitions[0].remove_cell(ell[0])
        ret_row_bar = self.cur_partitions[0].remove_cell(ell[n])
        for a in range(1, n - 1):
            ret_row_next = self.cur_partitions[a].remove_cell(ell[a])
            ret_row_bar_next = self.cur_partitions[a].remove_cell(ell[n + a])

            self._update_vacancy_numbers(a - 1)
            if ret_row is not None:
                self.cur_partitions[a-1].rigging[ret_row] = self.cur_partitions[a-1].vacancy_numbers[ret_row]
            if ret_row_bar is not None:
                self.cur_partitions[a-1].rigging[ret_row_bar] = self.cur_partitions[a-1].vacancy_numbers[ret_row_bar]

            ret_row = ret_row_next
            ret_row_bar = ret_row_bar_next

        ret_row_next = self.cur_partitions[n-1].remove_cell(ell[n-1])

        self._update_vacancy_numbers(n - 2)
        if ret_row is not None:
            self.cur_partitions[n-2].rigging[ret_row] = self.cur_partitions[n-2].vacancy_numbers[ret_row]
        if ret_row_bar is not None:
            self.cur_partitions[n-2].rigging[ret_row_bar] = self.cur_partitions[n-2].vacancy_numbers[ret_row_bar]

        self._update_vacancy_numbers(n - 1)
        if ret_row_next is not None:
            self.cur_partitions[n-1].rigging[ret_row_next] = self.cur_partitions[n-1].vacancy_numbers[ret_row_next]

        return b
