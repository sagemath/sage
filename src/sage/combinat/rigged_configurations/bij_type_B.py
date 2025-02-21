# sage.doctest: needs sage.combinat sage.modules
r"""
Bijection classes for type `B_n^{(1)}`

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `B_n^{(1)}`.

AUTHORS:

- Travis Scrimshaw (2012-12-21): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['B', 3, 1], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_B import KRTToRCBijectionTypeB
    sage: bijection = KRTToRCBijectionTypeB(KRT(pathlist=[[-1,2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['B', 3, 1], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_B import RCToKRTBijectionTypeB
    sage: bijection = RCToKRTBijectionTypeB(RC(partition_list=[[],[],[]]))
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
from sage.combinat.rigged_configurations.bij_type_C import KRTToRCBijectionTypeC
from sage.combinat.rigged_configurations.bij_type_C import RCToKRTBijectionTypeC


class KRTToRCBijectionTypeB(KRTToRCBijectionTypeC):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `B_n^{(1)}`.
    """

    def run(self, verbose=False):
        """
        Run the bijection from a tensor product of KR tableaux to a rigged
        configuration.

        INPUT:

        - ``tp_krt`` -- a tensor product of KR tableaux

        - ``verbose`` -- (default: ``False``) display each step in the
          bijection

        EXAMPLES::

            sage: from sage.combinat.rigged_configurations.bij_type_B import KRTToRCBijectionTypeB
            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['B', 3, 1], [[2, 1]])
            sage: KRTToRCBijectionTypeB(KRT(pathlist=[[0,3]])).run()
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            -1[ ]-1
            <BLANKLINE>
            0[]0
            <BLANKLINE>
            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['B', 3, 1], [[3, 1]])
            sage: KRTToRCBijectionTypeB(KRT(pathlist=[[-2,3,1]])).run()
            <BLANKLINE>
            (/)
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[]0
            <BLANKLINE>

        TESTS:

        Check that :issue:`19384` is fixed::

            sage: RC = RiggedConfigurations(['B',3,1], [[3,1],[3,1]])
            sage: RC._test_bijection()
            sage: RC = RiggedConfigurations(['B',3,1], [[1,1],[3,1],[1,1]])
            sage: RC._test_bijection()
        """
        if verbose:
            from sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element \
                import TensorProductOfKirillovReshetikhinTableauxElement

        for cur_crystal in reversed(self.tp_krt):
            r = cur_crystal.parent().r()

            # Check if it is a spinor
            if r == self.n:
                # Perform the spinor bijection by converting to type A_{2n-1}^{(2)}
                #   doing the bijection there and pulling back
                from sage.combinat.rigged_configurations.bij_type_A2_odd import KRTToRCBijectionTypeA2Odd
                from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
                from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition

                if verbose:
                    print("====================")
                    if len(self.cur_path) == 0:
                        print(repr([]))  # Special case for displaying when the rightmost factor is a spinor
                    else:
                        print(repr(TensorProductOfKirillovReshetikhinTableauxElement(self.tp_krt.parent(), self.cur_path)))
                    print("--------------------")
                    print(repr(self.ret_rig_con))
                    print("--------------------\n")
                    print("Applying doubling map")

                # Convert to a type A_{2n-1}^{(2)} RC
                dims = self.cur_dims[:]
                dims.insert(0, [r, cur_crystal.parent().s()])
                KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 2*self.n-1, 2], dims)
                # Convert the n-th partition into a regular rigged partition
                self.ret_rig_con[-1] = RiggedPartition(self.ret_rig_con[-1]._list,
                                                       self.ret_rig_con[-1].rigging,
                                                       self.ret_rig_con[-1].vacancy_numbers)
                # Placeholder element
                elt = KRT(*[C.module_generators[0] for C in KRT.crystals])
                bij = KRTToRCBijectionTypeA2Odd(elt)
                bij.ret_rig_con = KRT.rigged_configurations()(*self.ret_rig_con, use_vacancy_numbers=True)
                bij.cur_path = self.cur_path
                bij.cur_dims = self.cur_dims
                for i in range(len(self.cur_dims)):
                    if bij.cur_dims[i][0] != self.n:
                        bij.cur_dims[i][1] *= 2
                for i in range(self.n-1):
                    for j in range(len(bij.ret_rig_con[i])):
                        bij.ret_rig_con[i]._list[j] *= 2
                        bij.ret_rig_con[i].rigging[j] *= 2
                        bij.ret_rig_con[i].vacancy_numbers[j] *= 2

                # Perform the type A_{2n-1}^{(2)} bijection
                r = cur_crystal.parent().r()
                # Iterate through the columns
                for col_number, cur_column in enumerate(reversed(cur_crystal.to_array(False))):
                    bij.cur_path.insert(0, [])  # Prepend an empty list
                    bij.cur_dims.insert(0, [0, 1])

                    # Note that we do not need to worry about iterating over columns
                    #   (see previous note about the data structure).
                    for letter in reversed(cur_column):
                        bij.cur_dims[0][0] += 1
                        val = letter.value  # Convert from a CrystalOfLetter to an Integer

                        if verbose:
                            print("====================")
                            print(repr(TensorProductOfKirillovReshetikhinTableauxElement(self.tp_krt.parent(), bij.cur_path)))
                            print("--------------------")
                            print(repr(bij.ret_rig_con))
                            print("--------------------\n")

                        # Build the next state
                        bij.cur_path[0].insert(0, [letter])  # Prepend the value
                        bij.next_state(val)

                    # If we've split off a column, we need to merge the current column
                    #   to the current crystal tableau
                    if col_number > 0:
                        for i, letter_singleton in enumerate(self.cur_path[0]):
                            bij.cur_path[1][i].insert(0, letter_singleton[0])
                        bij.cur_dims[1][1] += 1
                        bij.cur_path.pop(0)
                        bij.cur_dims.pop(0)

                        # And perform the inverse column splitting map on the RC
                        for a in range(self.n):
                            bij._update_vacancy_nums(a)

                if verbose:
                    print("====================")
                    print(repr(TensorProductOfKirillovReshetikhinTableauxElement(self.tp_krt.parent(), bij.cur_path)))
                    print("--------------------")
                    print(repr(bij.ret_rig_con))
                    print("--------------------\n")
                    print("Applying halving map")

                # Convert back to a type B_n^{(1)}
                for i in range(len(self.cur_dims)):
                    if bij.cur_dims[i][0] != self.n:
                        bij.cur_dims[i][1] //= 2
                for i in range(self.n-1):
                    for j in range(len(bij.ret_rig_con[i])):
                        bij.ret_rig_con[i]._list[j] //= 2
                        bij.ret_rig_con[i].rigging[j] //= 2
                        bij.ret_rig_con[i].vacancy_numbers[j] //= 2
                self.ret_rig_con = self.tp_krt.parent().rigged_configurations()(*bij.ret_rig_con, use_vacancy_numbers=True)
                # Make it mutable so we don't have to keep making copies, at the
                #   end of the bijection, we will make it immutable again
                self.ret_rig_con._set_mutable()
            else:
                # Perform the regular type B_n^{(1)} bijection
                # Iterate through the columns
                for col_number, cur_column in enumerate(reversed(cur_crystal.to_array(False))):
                    self.cur_path.insert(0, [])  # Prepend an empty list
                    self.cur_dims.insert(0, [0, 1])

                    # Note that we do not need to worry about iterating over columns
                    #   (see previous note about the data structure).
                    for letter in reversed(cur_column):
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
        self.ret_rig_con.set_immutable()  # Return it to immutable
        return self.ret_rig_con

    def next_state(self, val):
        r"""
        Build the next state for type `B_n^{(1)}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['B', 3, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_B import KRTToRCBijectionTypeB
            sage: bijection = KRTToRCBijectionTypeB(KRT(pathlist=[[-1,2]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [3])
            sage: bijection.next_state(3)
        """
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        # If it is a regular value, we follow the A_n rules
        if val > 0:
            KRTToRCBijectionTypeA.next_state(self, val)
            return

        pos_val = -val

        # Special case for 0
        if pos_val == 0:
            if len(self.ret_rig_con[pos_val - 1]) > 0:
                max_width = self.ret_rig_con[n-1][0]
            else:
                max_width = 1
            max_width = self.ret_rig_con[n-1].insert_cell(max_width)
            width_n = max_width + 1
            max_width = max_width // 2

            # Check to see if we need to make the new string quasi-singular
            if tableau_height != n-1:
                max_width = self.ret_rig_con[n-2].insert_cell(max_width)
            else:
                max_width = -1
            self._update_vacancy_nums(n - 1)
            self._update_partition_values(n - 1)

            # Check if we need to make the new string at n quasi-singular
            p = self.ret_rig_con[n-1]
            num_rows = len(p)
            # Note that max width is 1 less than the corresponding string length
            if max_width*2 + 1 != width_n:
                for i in range(num_rows):
                    if p._list[i] == width_n:
                        j = i+1
                        while j < num_rows and p._list[j] == width_n \
                          and p.vacancy_numbers[j] == p.rigging[j]:
                            j += 1
                        p.rigging[j-1] -= 1
                        break

            # Follow regular A_n rules
            for a in reversed(range(tableau_height, n-2)):
                max_width = self.ret_rig_con[a].insert_cell(max_width)
                self._update_vacancy_nums(a + 1)
                self._update_partition_values(a + 1)
            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)
            if tableau_height > 0:
                self._update_vacancy_nums(tableau_height-1)
                self._update_partition_values(tableau_height-1)
            return

        # Always add a cell to the first singular value in the first
        #   tableau we are updating.
        if len(self.ret_rig_con[pos_val - 1]) > 0:
            max_width = self.ret_rig_con[pos_val - 1][0] + 1
        else:
            max_width = 0

        # Add cells similar to type A_n but we move to the right until n
        for a in range(pos_val - 1, n - 1):
            max_width = self.ret_rig_con[a].insert_cell(max_width)

        # Handle the special behavior at n
        if pos_val != n:
            max_width = max_width * 2

        # Find either the quasi-singular string and the next largest singular,
        #   the largest singular string if smaller than the max width
        #   or the two largest singular strings
        singular_max_width = False
        case_QS = False
        # Note, case_QS and singular_max_width will never both be True
        p = self.ret_rig_con[n-1]
        num_rows = len(p)
        width_n = 0
        for i in range(num_rows + 1):
            if i == num_rows:
                if case_QS:
                    # If we are in case (QS), we will be adding a box
                    p._list.append(1)
                    p.vacancy_numbers.append(None)
                    p.rigging.append(None)
                    width_n = 1
                    max_width = 0
                elif not singular_max_width:
                    # If we have not found a (quasi)singular string, we must add 2 boxes
                    # Go through our partition until we find a length of greater than 2
                    j = len(p._list) - 1
                    while j >= 0 and p._list[j] <= 2:
                        j -= 1
                    p._list.insert(j+1, 2)
                    p.vacancy_numbers.insert(j+1, None)
                    p.rigging.insert(j+1, None)
                    max_width = 0
                break
            elif p.vacancy_numbers[i] == p.rigging[i]:
                if p._list[i] < max_width:
                    if singular_max_width:
                        width_n = p._list[i]
                        break

                    max_width = p._list[i]
                    if case_QS:
                        p._list[i] += 1
                        p.rigging[i] = None
                        width_n = max_width + 1
                    else:
                        # Add 2 boxes
                        j = i - 1
                        while j >= 0 and p._list[j] <= max_width + 2:
                            p.rigging[j+1] = p.rigging[j]  # Shuffle it along
                            j -= 1
                        p._list.pop(i)
                        p._list.insert(j+1, max_width + 2)
                        p.rigging[j+1] = None
                    break

                if p._list[i] == max_width and not singular_max_width:
                    p._list[i] += 1  # We always at least add a box to the first singular value
                    p.rigging[i] = None
                    if case_QS:
                        width_n = p._list[i]
                        break
                    singular_max_width = True
                elif p._list[i] == max_width + 1 and not case_QS:
                    # If we can't add 2 boxes, we must be in case (QS)
                    p._list[i] += 1
                    p.rigging[i] = None
                    width_n = max_width
                    case_QS = True
            elif p.vacancy_numbers[i] - 1 == p.rigging[i] and not case_QS and not singular_max_width and p._list[i] <= max_width:
                case_QS = True
                max_width = p._list[i]
                p._list[i] += 1
                p.rigging[i] = None

        if singular_max_width:
            # There are 2 possibilities, case (S) and case (QS), we might need
            #   to attempt both
            # Make a *deep* copy of the element
            cp = self.ret_rig_con.__copy__()
            for i, rp in enumerate(cp):
                cp[i] = rp._clone()
            # We attempt case (S) first
            self._insert_cell_case_S(p)

        max_width = max_width // 2

        # We need to do the next partition in order to determine the step at n
        if tableau_height != n-1:
            max_width = self.ret_rig_con[n-2].insert_cell(max_width)
        else:
            max_width = -1

        self._update_vacancy_nums(n - 1)
        self._update_partition_values(n - 1)

        # If we need to make the smaller added string quasisingular
        # Note that max width is 1 less than the corresponding string length
        if case_QS and max_width*2 + 1 != width_n:
            for i in range(num_rows):
                if p._list[i] == width_n:
                    j = i+1
                    while j < num_rows and p._list[j] == width_n \
                      and p.vacancy_numbers[j] == p.rigging[j]:
                        j += 1
                    p.rigging[j-1] -= 1
                    break

        # Continue back following the regular A_n rules
        for a in reversed(range(tableau_height, n - 2)):
            max_width = self.ret_rig_con[a].insert_cell(max_width)
            self._update_vacancy_nums(a + 1)
            self._update_partition_values(a + 1)

        # Update the final rigged partitions
        if tableau_height < n:
            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

        assert pos_val > 0
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

        if singular_max_width:
            try:
                self.ret_rig_con.check()
            except Exception:
                self.other_outcome(cp, pos_val, width_n)

    def other_outcome(self, rc, pos_val, width_n):
        r"""
        Do the other case `(QS)` possibility.

        This arises from the ambiguity when we found a singular string at the
        max width in `\nu^{(n)}`. We had first attempted case `(S)`, and if
        that resulted in an invalid rigged configuration, we now
        finish the bijection using case `(QS)`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['B',3,1], [[2,1],[1,2]])
            sage: rc = RC(partition_list=[[2,1], [2,1,1], [5,1]])
            sage: t = rc.to_tensor_product_of_kirillov_reshetikhin_tableaux()
            sage: t.to_rigged_configuration() == rc # indirect doctest
            True
        """
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1
        self.ret_rig_con = rc

        # We need to do the next partition in order to determine the step at n
        max_width = self.ret_rig_con[n-2].insert_cell(width_n // 2)

        # We now attempt case (QS)
        case_QS = False
        p = self.ret_rig_con[n-1]
        num_rows = len(p)
        for i in range(len(p._list)):
            if p._list[i] == width_n:
                p._list[i] += 1
                p.rigging[i] = None
                case_QS = True
                break
        if not case_QS:  # we have not added a box yet
            p._list.append(1)
            p.rigging.append(None)
            p.vacancy_numbers.append(None)
            case_QS = True
        width_n += 1

        self._update_vacancy_nums(n - 1)
        self._update_partition_values(n - 1)

        # If we need to make the smaller added string quasisingular
        # Note that max width is 1 less than the corresponding string length
        if case_QS and max_width*2 + 1 != width_n:
            for i in range(num_rows):
                if p._list[i] == width_n:
                    j = i+1
                    while j < num_rows and p._list[j] == width_n \
                      and p.vacancy_numbers[j] == p.rigging[j]:
                        j += 1
                    p.rigging[j-1] -= 1
                    break

        # Continue back following the regular A_n rules
        for a in reversed(range(tableau_height, n - 2)):
            max_width = self.ret_rig_con[a].insert_cell(max_width)
            self._update_vacancy_nums(a + 1)
            self._update_partition_values(a + 1)

        # Update the final rigged partitions
        if tableau_height < n:
            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

        assert pos_val > 0
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


class RCToKRTBijectionTypeB(RCToKRTBijectionTypeC):
    r"""
    Specific implementation of the bijection from rigged configurations to
    tensor products of KR tableaux for type `B_n^{(1)}`.
    """

    def run(self, verbose=False, build_graph=False):
        """
        Run the bijection from rigged configurations to tensor product of KR
        tableaux for type `B_n^{(1)}`.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); display each step in the
          bijection
        - ``build_graph`` -- boolean (default: ``False``); build the graph of each
          step of the bijection

        EXAMPLES::

            sage: RC = RiggedConfigurations(['B', 3, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_B import RCToKRTBijectionTypeB
            sage: RCToKRTBijectionTypeB(RC(partition_list=[[1],[1,1],[1]])).run()
            [[3], [0]]

            sage: RC = RiggedConfigurations(['B', 3, 1], [[3, 1]])
            sage: x = RC(partition_list=[[],[1],[1]])
            sage: RCToKRTBijectionTypeB(x).run()
            [[1], [3], [-2]]
            sage: bij = RCToKRTBijectionTypeB(x)
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

            # Check to see if we are a spinor
            if dim[0] == self.n:
                # Perform the spinor bijection by converting to type A_{2n-1}^{(2)}
                #   doing the bijection there and pulling back

                from sage.combinat.rigged_configurations.bij_type_A2_odd import RCToKRTBijectionTypeA2Odd
                from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
                from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition, RiggedPartitionTypeB

                # Convert to a type A_{2n-1}^{(2)} RC
                RC = RiggedConfigurations(['A', 2*self.n-1, 2], self.cur_dims)
                if verbose:
                    print("====================")
                    print(repr(RC(*self.cur_partitions, use_vacancy_numbers=True)))
                    print("--------------------")
                    print(ret_crystal_path)
                    print("--------------------\n")
                    print("Applying doubling map\n")
                # Convert the n-th partition into a regular rigged partition
                self.cur_partitions[-1] = RiggedPartition(self.cur_partitions[-1]._list,
                                                          self.cur_partitions[-1].rigging,
                                                          self.cur_partitions[-1].vacancy_numbers)

                bij = RCToKRTBijectionTypeA2Odd(RC(*self.cur_partitions, use_vacancy_numbers=True))
                for i in range(len(self.cur_dims)):
                    if bij.cur_dims[i][0] != self.n:
                        bij.cur_dims[i][1] *= 2
                for i in range(self.n-1):
                    for j in range(len(bij.cur_partitions[i])):
                        bij.cur_partitions[i]._list[j] *= 2
                        bij.cur_partitions[i].rigging[j] *= 2
                        bij.cur_partitions[i].vacancy_numbers[j] *= 2

                if build_graph:
                    y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                    self._graph.append([self._graph[-1][1], (y, len(self._graph)), '2x'])

                # Perform the type A_{2n-1}^{(2)} bijection

                # Iterate over each column
                for dummy_var in range(dim[1]):
                    # Split off a new column if necessary
                    if bij.cur_dims[0][1] > 1:
                        bij.cur_dims[0][1] -= 1
                        bij.cur_dims.insert(0, [dim[0], 1])

                        # Perform the corresponding splitting map on rigged configurations
                        # All it does is update the vacancy numbers on the RC side
                        for a in range(self.n):
                            bij._update_vacancy_numbers(a)

                        if build_graph:
                            y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                            self._graph.append([self._graph[-1][1], (y, len(self._graph)), 'ls'])

                    while bij.cur_dims[0][0]:  # > 0:
                        if verbose:
                            print("====================")
                            print(repr(RC(*bij.cur_partitions, use_vacancy_numbers=True)))
                            print("--------------------")
                            print(ret_crystal_path)
                            print("--------------------\n")

                        ht = bij.cur_dims[0][0]
                        bij.cur_dims[0][0] = bij._next_index(ht)
                        b = bij.next_state(ht)
                        # Make sure we have a crystal letter
                        ret_crystal_path[-1].append(letters(b))  # Append the rank

                        if build_graph:
                            y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                            self._graph.append([self._graph[-1][1], (y, len(self._graph)), letters(b)])

                    bij.cur_dims.pop(0)  # Pop off the leading column

                self.cur_dims.pop(0)  # Pop off the spin rectangle

                self.cur_partitions = bij.cur_partitions
                # Convert the n-th partition back into the special type B one
                self.cur_partitions[-1] = RiggedPartitionTypeB(self.cur_partitions[-1])

                # Convert back to a type B_n^{(1)}
                if verbose:
                    print("====================")
                    print(repr(self.rigged_con.parent()(*bij.cur_partitions, use_vacancy_numbers=True)))
                    print("--------------------")
                    print(ret_crystal_path)
                    print("--------------------\n")
                    print("Applying halving map\n")

                for i in range(self.n-1):
                    for j in range(len(self.cur_partitions[i])):
                        self.cur_partitions[i]._list[j] //= 2
                        self.cur_partitions[i].rigging[j] //= 2
                        self.cur_partitions[i].vacancy_numbers[j] //= 2

                if build_graph:
                    y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                    self._graph.append([self._graph[-1][1], (y, len(self._graph)), '1/2x'])
            else:
                # Perform the regular type B_n^{(1)} bijection

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
                            self._graph.append([self._graph[-1][1], (y, len(self._graph)), '2x'])

                    while self.cur_dims[0][0]:  # > 0:
                        if verbose:
                            print("====================")
                            print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
                            print("--------------------")
                            print(ret_crystal_path)
                            print("--------------------\n")

                        self.cur_dims[0][0] -= 1  # This takes care of the indexing
                        b = self.next_state(self.cur_dims[0][0])

                        # Make sure we have a crystal letter
                        ret_crystal_path[-1].append(letters(b))  # Append the rank

                        if build_graph:
                            y = self.rigged_con.parent()(*[x._clone() for x in self.cur_partitions], use_vacancy_numbers=True)
                            self._graph.append([self._graph[-1][1], (y, len(self._graph)), letters(b)])

                    self.cur_dims.pop(0)  # Pop off the leading column

        if build_graph:
            self._graph.pop(0)  # Remove the dummy at the start
            from sage.graphs.digraph import DiGraph
            from sage.graphs.dot2tex_utils import have_dot2tex
            self._graph = DiGraph(self._graph)
            if have_dot2tex():
                self._graph.set_latex_options(format='dot2tex', edge_labels=True)

        return self.KRT(pathlist=ret_crystal_path)

    def next_state(self, height):
        r"""
        Build the next state for type `B_n^{(1)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['B', 3, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_B import RCToKRTBijectionTypeB
            sage: bijection = RCToKRTBijectionTypeB(RC(partition_list=[[1],[1,1],[1]]))
            sage: bijection.next_state(0)
            0
        """
        n = self.n
        ell = [None] * (2*n)
        case_S = False
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

        # Special case for n
        if b is None:
            last_size = 2 * last_size - 1
            partition = self.cur_partitions[n-1]
            # Modified version of _find_singular_string()
            for i in reversed(range(len(partition))):
                if partition[i] == last_size \
                  and partition.vacancy_numbers[i] == partition.rigging[i]:
                    case_Q = True
                    ell[n-1] = i
                elif partition[i] > last_size:
                    if not case_Q and partition.vacancy_numbers[i] - 1 == partition.rigging[i]:
                        case_Q = True
                        # Check if the block is singular as well
                        block_size = partition[i]
                        for j in reversed(range(i)):
                            if partition[j] != block_size:
                                break
                            elif partition.vacancy_numbers[j] == partition.rigging[j]:
                                case_Q = False
                                ell[2*n-1] = j
                                last_size = partition[j]
                                case_S = True
                                break
                        if not case_Q:  # We found a singular string above the quasi-singular one
                            break
                        ell[n-1] = i
                        last_size = partition[i]
                        # Now check for case QS
                    elif partition.vacancy_numbers[i] == partition.rigging[i]:
                        ell[2*n-1] = i
                        last_size = partition[i]
                        case_S = True
                        break

            if ell[2*n-1] is None:
                if not case_Q:
                    b = n
                else:
                    b = 0

        if b is None:
            # Now go back
            last_size = (last_size + 1) // 2
            for a in reversed(range(n - 1)):
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

        # Determine the new rigged configuration by removing boxes from the
        #   selected string and then making the new string singular

        # Determine if we need to make the n-th string quasisingular
        make_quasisingular = case_Q and case_S and \
                (ell[2*n-2] is None
                 or self.cur_partitions[n-1][ell[2*n-1]]
                    < 2*self.cur_partitions[n-2][ell[2*n-2]])

        row_num = self.cur_partitions[0].remove_cell(ell[0])
        row_num_bar = self.cur_partitions[0].remove_cell(ell[n])
        for a in range(1, n-1):
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
            if case_S:
                row_num_next = self.cur_partitions[n-1].remove_cell(ell[n-1])
                row_num_bar_next = self.cur_partitions[n-1].remove_cell(ell[2*n-1])
            else:
                row_num_next = self.cur_partitions[n-1].remove_cell(ell[n-1])
                row_num_bar_next = None
        elif case_S:
            row_num_next = self.cur_partitions[n-1].remove_cell(ell[2*n-1], 2)
            row_num_bar_next = None
        else:
            row_num_next = None
            row_num_bar_next = None

        self._update_vacancy_numbers(n - 2)
        if row_num is not None:
            self.cur_partitions[n-2].rigging[row_num] = self.cur_partitions[n-2].vacancy_numbers[row_num]
        if row_num_bar is not None:
            self.cur_partitions[n-2].rigging[row_num_bar] = self.cur_partitions[n-2].vacancy_numbers[row_num_bar]

        self._update_vacancy_numbers(n - 1)
        if row_num_next is not None:
            self.cur_partitions[n-1].rigging[row_num_next] = self.cur_partitions[n-1].vacancy_numbers[row_num_next]
        if row_num_bar_next is not None:  # If we enter here, it means case (Q, S) holds
            vac_num = self.cur_partitions[n-1].vacancy_numbers[row_num_bar_next]
            self.cur_partitions[n-1].rigging[row_num_bar_next] = vac_num
            if make_quasisingular:
                block_len = self.cur_partitions[n-1][row_num_bar_next]
                j = row_num_bar_next + 1
                length = len(self.cur_partitions[n-1])
                # Find the place for the quasisingular rigging
                while j < length and self.cur_partitions[n-1][j] == block_len \
                  and self.cur_partitions[n-1].rigging[j] == vac_num:
                    j += 1
                self.cur_partitions[n-1].rigging[j-1] = vac_num - 1

        return b
