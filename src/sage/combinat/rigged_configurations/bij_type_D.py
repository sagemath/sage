# sage.doctest: needs sage.combinat sage.modules
r"""
Bijection classes for type `D_n^{(1)}`

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `D_n^{(1)}`.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
    sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[3, 2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
    sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[],[],[],[]]))
    sage: TestSuite(bijection).run()
"""

# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA


class KRTToRCBijectionTypeD(KRTToRCBijectionTypeA):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `D_n^{(1)}`.

    This inherits from type `A_n^{(1)}` because we use the same methods in
    some places.
    """

    def run(self, verbose=False):
        """
        Run the bijection from a tensor product of KR tableaux to a rigged
        configuration for type `D_n^{(1)}`.

        INPUT:

        - ``tp_krt`` -- a tensor product of KR tableaux

        - ``verbose`` -- (default: ``False``) display each step in the
          bijection

        EXAMPLES::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
            sage: KRTToRCBijectionTypeD(KRT(pathlist=[[-3,2]])).run()
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            2[ ]2
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
        """
        if verbose:
            from sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element \
              import TensorProductOfKirillovReshetikhinTableauxElement

        for cur_crystal in reversed(self.tp_krt):
            r = cur_crystal.parent().r()
            # Iterate through the columns
            for col_number, cur_column in enumerate(reversed(cur_crystal.to_array(False))):
                self.cur_path.insert(0, [])  # Prepend an empty list

                # Check to see if we are a spinor column
                if r >= self.n-1:
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
                    # This check is needed for the n-1 spin column
                    if self.cur_dims[0][0] < r:
                        self.cur_dims[0][0] += 1
                    val = letter.value  # Convert from a CrystalOfLetter to an Integer

                    if verbose:
                        print("====================")
                        print(repr(TensorProductOfKirillovReshetikhinTableauxElement(self.tp_krt.parent(), self.cur_path)))
                        print("--------------------")
                        print(repr(self.ret_rig_con))
                        print("--------------------\n")

                    # Build the next state
                    self.cur_path[0].insert(0, [letter])  # Prepend the value
                    self.next_state(val)

                # Check to see if we are a spinor column
                if r >= self.n-1:
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

        self.ret_rig_con.set_immutable()  # Return it to immutable
        return self.ret_rig_con

    def next_state(self, val):
        r"""
        Build the next state for type `D_n^{(1)}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
            sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[5,3]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [3])
            sage: bijection.next_state(3)
        """
        # Note that type D_n only contains the (absolute) values between 1 and n
        #   (unlike type A_n which is 1 to n+1). Thus we only need to subtract 1
        #   to match the indices.
        # Also note that we must subtract 1 from n to match the indices as well.
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        # If it is a regular value, we follow the A_n rules
        if val > 0:
            KRTToRCBijectionTypeA.next_state(self, val)

            # If we are inserting n-1 (which will only update nu[a] for a <=
            #   n-1), we need to update the vacancy number of nu[n] as well
            if val == n - 1:
                self._update_vacancy_nums(val)
            if tableau_height >= n - 2:
                self._correct_vacancy_nums()
            return

        pos_val = -val

        if pos_val == n:
            # Special case for `\overline{n}` and adding to make height `n`
            #    where we only update the vacancy numbers
            # This only occurs with `r = n - 1`
            if self.cur_dims[0][0] == n - 1 and tableau_height == n - 1:
                self._update_vacancy_nums(n-2)
                self._update_vacancy_nums(n-1)
                self._correct_vacancy_nums()
                return

            if len(self.ret_rig_con[n - 1]) > 0:
                max_width = self.ret_rig_con[n - 1][0] + 1
            else:
                max_width = 1

            # Update the last one and skip a rigged partition
            max_width = self.ret_rig_con[n - 1].insert_cell(max_width)

            for a in reversed(range(tableau_height, n - 2)):
                max_width = self.ret_rig_con[a].insert_cell(max_width)
                self._update_vacancy_nums(a + 1)
                self._update_partition_values(a + 1)

            # Update all other effected remaining values

            self._update_vacancy_nums(n - 1)
            if tableau_height >= n - 2:
                self._correct_vacancy_nums()
            self._update_partition_values(n - 1)

            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

            if tableau_height > 0:
                self._update_vacancy_nums(tableau_height - 1)
                self._update_partition_values(tableau_height - 1)

            self._update_vacancy_nums(n - 2)
            return

        # Always add a cell to the first singular value in the first
        #   tableau we are updating.
        if len(self.ret_rig_con[pos_val - 1]) > 0:
            max_width = self.ret_rig_con[pos_val - 1][0] + 1
        else:
            max_width = 1
        # Special case for `\overline{n-1}` to take the larger of the last two
        if pos_val == n - 1 and len(self.ret_rig_con[n - 1]) > 0 and \
          self.ret_rig_con[n - 1][0] + 1 > max_width:
            max_width = self.ret_rig_con[n - 1][0] + 1

        # Add cells similar to type A_n but we move to the right until we reach
        #   the value of n-2
        for a in range(pos_val - 1, n - 2):
            max_width = self.ret_rig_con[a].insert_cell(max_width)

        # Handle the special behavior near values of n
        if tableau_height <= n - 2:
            max_width2 = self.ret_rig_con[n - 2].insert_cell(max_width)
            max_width = self.ret_rig_con[n - 1].insert_cell(max_width)
            max_width = min(max_width2, max_width)
        elif pos_val <= self.cur_dims[0][0]:
            # Special case when the height will become n
            max_width = self.ret_rig_con[self.cur_dims[0][0] - 1].insert_cell(max_width)

        # Go back following the regular A_n rules
        if tableau_height <= n - 3:
            max_width = self.ret_rig_con[n - 3].insert_cell(max_width)

        self._update_vacancy_nums(n - 2)
        self._update_vacancy_nums(n - 1)
        if tableau_height >= n - 2:
            self._correct_vacancy_nums()
        self._update_partition_values(n - 2)
        self._update_partition_values(n - 1)

        for a in reversed(range(tableau_height, n - 3)):
            max_width = self.ret_rig_con[a].insert_cell(max_width)
            self._update_vacancy_nums(a + 1)
            self._update_partition_values(a + 1)

        # Update the final rigged partitions
        if tableau_height < n - 2:
            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

            if pos_val <= tableau_height:
                for a in range(pos_val - 1, tableau_height):
                    self._update_vacancy_nums(a)
                    self._update_partition_values(a)
                if pos_val > 1:
                    self._update_vacancy_nums(pos_val-2)
                    self._update_partition_values(pos_val-2)
            elif 0 < tableau_height:
                self._update_vacancy_nums(tableau_height - 1)
                self._update_partition_values(tableau_height - 1)
        elif pos_val <= n - 1:
            for a in range(pos_val - 1, n - 2):
                self._update_vacancy_nums(a)
                self._update_partition_values(a)
            if pos_val > 1:
                self._update_vacancy_nums(pos_val-2)
                self._update_partition_values(pos_val-2)

    def _correct_vacancy_nums(self):
        r"""
        Correct the vacancy numbers with special considerations for spinor
        columns.

        This should only be called when we are going to have a (left-most)
        spinor column of height `n-1` or `n` in type `D^{(1)_n`.

        This is a correction for the spinor column where we consider the
        weight `\overline{\Lambda_k}` where `k = n-1,n` during the spinor
        bijection. This adds 1 to each of the respective vacancy numbers
        to account for this.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
            sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[-1,4,3,2]]))
            sage: bijection.doubling_map()
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [2])
            sage: bijection.next_state(2) # indirect doctest
            sage: bijection.ret_rig_con
            <BLANKLINE>
            -2[ ]-2
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        pos = self.n - 2
        if self.cur_dims[0][0] == len(self.cur_path[0]):
            # The current r value is never greater than the height of the current column
            # Thus if we do not enter into this if block, then r < height and
            #   we adjust the (n-1)-th partition.
            pos += 1

        for i in range(len(self.ret_rig_con[pos]._list)):
            self.ret_rig_con[pos].vacancy_numbers[i] += 1

    def doubling_map(self):
        r"""
        Perform the doubling map of the rigged configuration at the current
        state of the bijection.

        This is the map `B(\Lambda) \hookrightarrow B(2 \Lambda)` which
        doubles each of the rigged partitions and updates the vacancy numbers
        accordingly.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[4,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
            sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[-1,4,3,2]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [2])
            sage: bijection.next_state(2)
            sage: bijection.ret_rig_con
            <BLANKLINE>
            -2[ ]-2
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: bijection.cur_dims
            [[0, 1]]
            sage: bijection.doubling_map()
            sage: bijection.ret_rig_con
            <BLANKLINE>
            -4[ ][ ]-4
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: bijection.cur_dims
            [[0, 2]]
        """
        for i in range(len(self.cur_dims)):
            self.cur_dims[i][1] *= 2
        for i in range(len(self.ret_rig_con)):
            for j in range(len(self.ret_rig_con[i])):
                self.ret_rig_con[i]._list[j] *= 2
                self.ret_rig_con[i].rigging[j] *= 2
                self.ret_rig_con[i].vacancy_numbers[j] *= 2

    def halving_map(self):
        r"""
        Perform the halving map of the rigged configuration at the current
        state of the bijection.

        This is the inverse map to `B(\Lambda) \hookrightarrow B(2 \Lambda)`
        which halves each of the rigged partitions and updates the vacancy
        numbers accordingly.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[4,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
            sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[-1,4,3,2]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [2])
            sage: bijection.next_state(2)
            sage: test = bijection.ret_rig_con
            sage: bijection.doubling_map()
            sage: bijection.halving_map()
            sage: test == bijection.ret_rig_con
            True
        """
        # Skip the first column since it is a spinor
        for i in range(1, len(self.cur_dims)):
            self.cur_dims[i][1] //= 2
        for i in range(len(self.ret_rig_con)):
            for j in range(len(self.ret_rig_con[i])):
                self.ret_rig_con[i]._list[j] //= 2
                self.ret_rig_con[i].rigging[j] //= 2
                self.ret_rig_con[i].vacancy_numbers[j] //= 2


class RCToKRTBijectionTypeD(RCToKRTBijectionTypeA):
    r"""
    Specific implementation of the bijection from rigged configurations to tensor products of KR tableaux for type `D_n^{(1)}`.
    """

    def run(self, verbose=False, build_graph=False):
        """
        Run the bijection from rigged configurations to tensor product of KR
        tableaux for type `D_n^{(1)}`.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); display each step in the
          bijection
        - ``build_graph`` -- boolean (default: ``False``); build the graph of each
          step of the bijection

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
            sage: x = RC(partition_list=[[1],[1],[1],[1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
            sage: RCToKRTBijectionTypeD(x).run()
            [[2], [-3]]
            sage: bij = RCToKRTBijectionTypeD(x)
            sage: bij.run(build_graph=True)
            [[2], [-3]]
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
                if dim[0] >= self.n - 1:
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

                    if dim[0] == self.n - 1:
                        if verbose:
                            print("====================")
                            print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
                            print("--------------------")
                            print(ret_crystal_path)
                            print("--------------------\n")
                        b = self.next_state(self.n)
                        if b == self.n:
                            b = -self.n
                        ret_crystal_path[-1].append(letters(b))  # Append the rank

                        if build_graph:
                            y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                            self._graph.append([self._graph[-1][1], (y, len(self._graph)), letters(b)])

                while self.cur_dims[0][0] > 0:
                    if verbose:
                        print("====================")
                        print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
                        print("--------------------")
                        print(ret_crystal_path)
                        print("--------------------\n")

                    self.cur_dims[0][0] -= 1  # This takes care of the indexing
                    b = self.next_state(self.cur_dims[0][0])

                    # Corrections for spinor
                    if dim[0] == self.n and b == -self.n \
                      and self.cur_dims[0][0] == self.n - 1:
                        b = -(self.n-1)

                    # Make sure we have a crystal letter
                    ret_crystal_path[-1].append(letters(b))  # Append the rank

                    if build_graph:
                        y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                        self._graph.append([self._graph[-1][1], (y, len(self._graph)), letters(b)])

                self.cur_dims.pop(0)  # Pop off the leading column

                # Check to see if we were a spinor
                if dim[0] >= self.n-1:
                    if verbose:
                        print("====================")
                        print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
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
            self._graph = DiGraph(self._graph, format='list_of_edges')
            if have_dot2tex():
                self._graph.set_latex_options(format='dot2tex', edge_labels=True)

        return self.KRT(pathlist=ret_crystal_path)

    def next_state(self, height):
        r"""
        Build the next state for type `D_n^{(1)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
            sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[],[1,1],[1],[1]]))
            sage: bijection.next_state(0)
            1
        """
        n = self.n
        ell = [None] * (2 * n - 2) # No `\bar{\ell}^{n-1}` and `\bar{\ell}^n`
        b = None

        # Calculate the rank and ell values

        last_size = 0
        for a in range(height, n - 2):
            ell[a] = self._find_singular_string(self.cur_partitions[a], last_size)

            if ell[a] is None:
                b = a + 1
                break
            else:
                last_size = self.cur_partitions[a][ell[a]]

        if height == n:
            # Special case from height `n` spinor with `r = n-1`
            ell[n - 2] = self._find_singular_string(self.cur_partitions[n - 2], last_size)
            if ell[n - 2] is not None:
                last_size = self.cur_partitions[n - 2][ell[n - 2]]
            else:
                b = -n
        elif height == n - 1:
            # Special case for height `n-1` spinor
            ell[n - 1] = self._find_singular_string(self.cur_partitions[n - 1], last_size)
            if ell[n - 1] is not None:
                last_size = self.cur_partitions[n - 1][ell[n - 1]]
            else:
                b = n
        elif b is None:
            # Do the special cases when we've reached n - 2
            ell[n - 2] = self._find_singular_string(self.cur_partitions[n - 2], last_size)
            ell[n - 1] = self._find_singular_string(self.cur_partitions[n - 1], last_size)

            if ell[n - 2] is not None:
                temp_size = self.cur_partitions[n - 2][ell[n - 2]]
                if ell[n - 1] is not None:
                    last_size = self.cur_partitions[n - 1][ell[n - 1]]
                    last_size = max(temp_size, last_size)
                else:
                    b = n
            else:
                if ell[n - 1] is not None:
                    b = -n
                else:
                    b = n - 1

        if b is None:
            # Now go back
            for a in reversed(range(n - 2)):
                # Modified form of _find_singular_string
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
        for a in range(1, n - 2):
            ret_row_next = self.cur_partitions[a].remove_cell(ell[a])
            ret_row_bar_next = self.cur_partitions[a].remove_cell(ell[n + a])

            self._update_vacancy_numbers(a - 1)
            if ret_row is not None:
                self.cur_partitions[a - 1].rigging[ret_row] = \
                  self.cur_partitions[a - 1].vacancy_numbers[ret_row]
            if ret_row_bar is not None:
                self.cur_partitions[a - 1].rigging[ret_row_bar] = \
                  self.cur_partitions[a - 1].vacancy_numbers[ret_row_bar]

            ret_row = ret_row_next
            ret_row_bar = ret_row_bar_next

        # Special behavior for all a > n-3
        ret_row_next = self.cur_partitions[n - 2].remove_cell(ell[n - 2])
        ret_row_bar_next = self.cur_partitions[n - 1].remove_cell(ell[n - 1])

        self._update_vacancy_numbers(n - 3)
        if ret_row is not None:
            self.cur_partitions[n - 3].rigging[ret_row] = \
              self.cur_partitions[n - 3].vacancy_numbers[ret_row]
        if ret_row_bar is not None:
            self.cur_partitions[n - 3].rigging[ret_row_bar] = \
              self.cur_partitions[n - 3].vacancy_numbers[ret_row_bar]

        self._update_vacancy_numbers(n - 2)
        if ret_row_next is not None:
            self.cur_partitions[n - 2].rigging[ret_row_next] = \
              self.cur_partitions[n - 2].vacancy_numbers[ret_row_next]

        self._update_vacancy_numbers(n - 1)
        if height >= n - 1:
            self._correct_vacancy_nums()
        if ret_row_bar_next is not None:
            self.cur_partitions[n - 1].rigging[ret_row_bar_next] = \
              self.cur_partitions[n - 1].vacancy_numbers[ret_row_bar_next]

        return b

    def doubling_map(self):
        r"""
        Perform the doubling map of the rigged configuration at the current
        state of the bijection.

        This is the map `B(\Lambda) \hookrightarrow B(2 \Lambda)` which
        doubles each of the rigged partitions and updates the vacancy numbers
        accordingly.

        TESTS::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[4, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
            sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[],[],[],[1]]))
            sage: bijection.cur_partitions
            [(/)
             , (/)
             , (/)
             , -1[ ]-1
             ]
            sage: bijection.doubling_map()
            sage: bijection.cur_partitions
            [(/)
             , (/)
             , (/)
             , -2[ ][ ]-2
             ]
        """
        # Skip the first column since it is a spinor
        for i in range(1, len(self.cur_dims)):
            self.cur_dims[i][1] *= 2
        for partition in self.cur_partitions:
            for j in range(len(partition)):
                partition._list[j] *= 2
                partition.rigging[j] *= 2
                partition.vacancy_numbers[j] *= 2

    def halving_map(self):
        r"""
        Perform the halving map of the rigged configuration at the current
        state of the bijection.

        This is the inverse map to `B(\Lambda) \hookrightarrow B(2 \Lambda)`
        which halves each of the rigged partitions and updates the vacancy
        numbers accordingly.

        TESTS::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[4, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
            sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[],[],[],[1]]))
            sage: test = bijection.cur_partitions
            sage: bijection.doubling_map()
            sage: bijection.halving_map()
            sage: test == bijection.cur_partitions
            True
        """
        for i in range(len(self.cur_dims)):
            self.cur_dims[i][1] //= 2
        for partition in self.cur_partitions:
            for j in range(len(partition)):
                partition._list[j] //= 2
                partition.rigging[j] //= 2
                partition.vacancy_numbers[j] //= 2

    def _correct_vacancy_nums(self):
        r"""
        Correct the vacancy numbers with special considerations for spinor
        columns.

        This should only be called when we are going to have a (left-most)
        spinor column of height `n-1` or `n`.

        This is a correction for the spinor column where we consider the
        weight `\overline{\Lambda_k}` where `k = n-1,n` during the spinor
        bijection. This adds 1 to each of the respective vacancy numbers
        to account for this.

        TESTS::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[4, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
            sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[],[],[],[1]]))
            sage: bijection.doubling_map()
            sage: bijection.next_state(4) # indirect doctest
            -4
        """
        n = self.n
        for i in range(len(self.cur_partitions[n-1]._list)):
            self.cur_partitions[n-1].vacancy_numbers[i] += 1
