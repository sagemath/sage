"""
Python interface to partition backtrack functions

EXAMPLES::

    sage: import sage.groups.perm_gps.partn_ref.refinement_python

This module provides Python frontends to the Cython-based partition backtrack
functions. This allows one to write the three input functions
(all_children_are_equivalent, refine_and_return_invariant, and compare_structures)
in pure Python, and still use the Cython algorithms. Experimentation with
specific partition backtrack implementations no longer requires compilation, as
the input functions can be dynamically changed at runtime.

.. NOTE::

    This is not intended for production quality implementations of partition
    refinement, but instead for experimentation, learning, and use of the
    Python debugger.
"""

#*****************************************************************************
#       Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.memory cimport sig_malloc, sig_free

from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label cimport (
    get_aut_gp_and_can_lab, aut_gp_and_can_lab,
    deallocate_agcl_output)
from sage.groups.perm_gps.partn_ref.double_coset cimport double_coset
from sage.rings.integer cimport Integer


cdef class PythonPartitionStack:
    """
    Instances of this class wrap a (Cython) PartitionStack.
    """

    def __init__(self, int n):
        """
        Initialize a PartitionStack.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7) # implicit doctest
        """
        self.c_ps = PS_new(n, 1)

    def __dealloc__(self):
        """
        Deallocate the PartitionStack.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: del(P) # implicit doctest
        """
        PS_dealloc(self.c_ps)

    def __repr__(self):
        """
        Return a string representing the stack.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P # implicit doctest
            PythonPartitionStack of degree 7 and depth 0.
        """
        return "PythonPartitionStack of degree %d and depth %d."%(self.c_ps.degree, self.c_ps.depth)

    def display(self):
        """
        Print a representation of the stack.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.depth(1)
            1
            sage: P.set_level(2,1)
            sage: P.display()
            (0 1 2 3 4 5 6)
            (0 1 2|3 4 5 6)
        """
        PS_print(self.c_ps)

    def is_discrete(self):
        """
        Return whether the deepest partition consists only of singleton cells.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.is_discrete()
            False
            sage: [P.set_level(i,0) for i in range(7)]
            [None, None, None, None, None, None, None]
            sage: P.is_discrete()
            True
        """
        return PS_is_discrete(self.c_ps)

    def num_cells(self):
        """
        Return the number of cells in the deepest partition.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.num_cells()
            1
        """
        return PS_num_cells(self.c_ps)

    def move_min_to_front(self, int start, int end):
        """
        Make sure that the first element of the segment of entries i with
        start <= i <= end is minimal.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.set_entry(1,0)
            sage: P.set_entry(0,1)
            sage: P.display()
            (1 0 2 3 4 5 6)
            sage: P.move_min_to_front(0,1)
            sage: P.display()
            (0 1 2 3 4 5 6)
        """
        PS_move_min_to_front(self.c_ps, start, end)

    def __copy__(self):
        """

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: Q = copy(P)
            sage: P.display()
            (0 1 2 3 4 5 6)
            sage: Q.display()
            (0 1 2 3 4 5 6)
        """
        cdef PythonPartitionStack cpy
        cpy = PythonPartitionStack(self.c_ps.degree)
        PS_copy_from_to(self.c_ps, cpy.c_ps)
        return cpy

    def clear(self):
        """
        Set the current partition to the first shallower one, i.e. forget about
        boundaries between cells that are new to the current level.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.depth(1)
            1
            sage: P.set_level(2,1)
            sage: P.display()
            (0 1 2 3 4 5 6)
            (0 1 2|3 4 5 6)
            sage: P.clear()
            sage: P.display()
            (0 1 2 3 4 5 6)
            (0 1 2 3 4 5 6)
        """
        PS_clear(self.c_ps)

    def entries(self):
        """
        Return the entries array as a Python list of ints.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.entries()
            [0, 1, 2, 3, 4, 5, 6]
            sage: P.levels()
            [7, 7, 7, 7, 7, 7, -1]
        """
        cdef int i
        return [self.c_ps.entries[i] for i from 0 <= i < self.c_ps.degree]

    def set_entry(self, int i, int entry):
        """
        Set the `i`-th entry of the entries array to entry.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.set_entry(1,0)
            sage: P.set_entry(0,1)
            sage: P.display()
            (1 0 2 3 4 5 6)
        """
        self.c_ps.entries[i] = entry

    def get_entry(self, int i):
        """
        Get the `i`-th entry of the entries array.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.get_entry(0)
            0
        """
        return self.c_ps.entries[i]

    def levels(self):
        """
        Return the levels array as a Python list of ints.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.entries()
            [0, 1, 2, 3, 4, 5, 6]
            sage: P.levels()
            [7, 7, 7, 7, 7, 7, -1]
        """
        return [self.c_ps.levels[i] for i from 0 <= i < self.c_ps.degree]

    def set_level(self, int i, int level):
        """
        Set the `i`-th entry of the levels array to entry.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.depth(1)
            1
            sage: P.set_level(2,1)
            sage: P.display()
            (0 1 2 3 4 5 6)
            (0 1 2|3 4 5 6)
        """
        self.c_ps.levels[i] = level

    def get_level(self, int i):
        """
        Get the `i`-th entry of the levels array.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.get_level(0)
            7
        """
        return self.c_ps.levels[i]

    def depth(self, new=None):
        """
        Return the depth of the deepest partition in the stack, setting it to
        new if new is not None.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.depth()
            0
        """
        if new is not None:
            self.c_ps.depth = new
        return self.c_ps.depth

    def degree(self, new=None):
        """
        Return the degree of the partition stack, setting it to
        new if new is not None.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.degree()
            7
        """
        if new is not None:
            self.c_ps.degree = new
        return self.c_ps.degree

    def partition(self, int k):
        """
        Return the partition at level k, as a Python list of lists.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonPartitionStack
            sage: P = PythonPartitionStack(7)
            sage: P.depth(1)
            1
            sage: P.set_level(2,1)
            sage: P.partition(0)
            [[0, 1, 2, 3, 4, 5, 6]]
            sage: P.partition(1)
            [[0, 1, 2], [3, 4, 5, 6]]
        """
        cdef int i
        cdef list partition = [], cell = []
        for i from 0 <= i < self.c_ps.degree:
            cell.append(self.c_ps.entries[i])
            if self.c_ps.levels[i] <= k:
                partition.append(cell)
                if i < self.c_ps.degree:
                    cell = []
        return partition


class PythonObjectWrapper:
    """
    Instances of this class wrap a Python object and the refinement functions.
    """
    def __init__(self, obj, acae_fn, rari_fn, cs_fn, int degree):
        """
        Initialize a PythonObjectWrapper.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonObjectWrapper
            sage: def acae(a, b):
            ....:  return 0
            sage: def rari(a, b, c):
            ....:  return 0
            sage: def cs(a, b, c, d, e):
            ....:  return 0
            sage: from sage.groups.perm_gps.partn_ref.refinement_python import PythonObjectWrapper
            sage: P = PythonObjectWrapper(None, acae, rari, cs, 7) # implicit doctest
            sage: P.obj
            sage: P.degree
            7
            sage: P.acae_fn
            <function acae at ...>
            sage: P.rari_fn
            <function rari at ...>
            sage: P.cs_fn
            <function cs at ...>
        """
        self.degree = degree
        self.obj = obj
        self.acae_fn = acae_fn
        self.rari_fn = rari_fn
        self.cs_fn = cs_fn


cdef bint all_children_are_equivalent_python(PartitionStack *PS, void *S) noexcept:
    """
    Python conversion of all_children_are_equivalent function.
    """
    cdef PythonPartitionStack Py_PS = PythonPartitionStack(PS.degree)
    cdef object S_obj = <object> S
    PS_copy_from_to(PS, Py_PS.c_ps)
    return S_obj.acae_fn(Py_PS, S_obj.obj)

cdef int refine_and_return_invariant_python(PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len) noexcept:
    """
    Python conversion of refine_and_return_invariant function.
    """
    cdef PythonPartitionStack Py_PS = PythonPartitionStack(PS.degree)
    cdef object S_obj = <object> S
    PS_copy_from_to(PS, Py_PS.c_ps)
    cdef int i
    cdef list ctrb_py = [cells_to_refine_by[i] for i from 0 <= i < ctrb_len]
    return S_obj.rari_fn(Py_PS, S_obj.obj, ctrb_py)

cdef int compare_structures_python(int *gamma_1, int *gamma_2, void *S1, void *S2, int degree) noexcept:
    """
    Python conversion of compare_structures function.
    """
    cdef int i
    cdef object S1_obj = <object> S1, S2_obj = <object> S2
    cdef list gamma_1_py = [gamma_1[i] for i from 0 <= i < degree]
    cdef list gamma_2_py = [gamma_2[i] for i from 0 <= i < degree]
    return S1_obj.cs_fn(gamma_1_py, gamma_2_py, S1_obj.obj, S2_obj.obj, degree)


def aut_gp_and_can_lab_python(S, partition, n,
    all_children_are_equivalent,
    refine_and_return_invariant,
    compare_structures,
    canonical_label, base, order):
    """
    Call the automorphism group and canonical label function.

    INPUT:

    - ``S`` -- the object to examine
    - ``partition`` -- an ordered partition, as a list of lists
    - ``n`` -- the degree of the automorphism group to be computed

    ::

        all_children_are_equivalent -- Python function of "signature":
            bool all_children_are_equivalent(PythonPartitionStack, object)
        refine_and_return_invariant -- Python function of "signature":
            int refine_and_return_invariant(PythonPartitionStack, object, list)
        compare_structures -- Python function of "signature":
            int compare_structures(list, list, object, object)
        (see automorphism_group_canonical_label.pyx for more documentation)

    ::

        canonical_label -- boolean; whether to search for a canonical label
        base -- boolean; whether to return a base for the automorphism group
        order -- boolean; whether to return the order of the automorphism group

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_python import aut_gp_and_can_lab_python
        sage: def acae(a, b):
        ....:  return 0
        sage: def rari(a, b, c):
        ....:  return 0
        sage: def cs(a, b, c, d, e):
        ....:  return 0
        sage: aut_gp_and_can_lab_python(None, [[0,1,2,3],[4,5]], 6, acae, rari, cs, True, True, True)
        ([[0, 1, 3, 2, 4, 5],
          [0, 2, 1, 3, 4, 5],
          [1, 0, 2, 3, 4, 5],
          [0, 1, 2, 3, 5, 4]],
         [0, 1, 2, 3, 4, 5],
         [4, 0, 1, 2],
         48)
        sage: factorial(4)*factorial(2)
        48
    """
    obj_wrapper = PythonObjectWrapper(S, all_children_are_equivalent, refine_and_return_invariant, compare_structures, n)
    cdef aut_gp_and_can_lab *output
    cdef int i, j
    cdef Integer I

    cdef PartitionStack *part = PS_from_list(partition)
    if part is NULL:
        raise MemoryError

    output = get_aut_gp_and_can_lab(<void *> obj_wrapper, part, n,
        &all_children_are_equivalent_python,
        &refine_and_return_invariant_python,
        &compare_structures_python,
        canonical_label, NULL, NULL, NULL)

    list_of_gens = []
    for i in range(output.num_gens):
        list_of_gens.append([output.generators[j+i*n] for j in range(n)])
    return_tuple = [list_of_gens]
    if canonical_label:
        return_tuple.append([output.relabeling[i] for i in range(n)])
    if base:
        return_tuple.append([output.group.base_orbits[i][0] for i from 0 <= i < output.group.base_size])
    if order:
        I = Integer()
        SC_order(output.group, 0, I.value)
        return_tuple.append(I)
    PS_dealloc(part)
    deallocate_agcl_output(output)
    if len(return_tuple) == 1:
        return return_tuple[0]
    else:
        return tuple(return_tuple)


def double_coset_python(S1, S2, partition1, ordering2, n,
    all_children_are_equivalent,
    refine_and_return_invariant,
    compare_structures):
    """
    Call the double coset function.

    INPUT:

    - ``S1``, ``S2`` -- the objects to examine
    - ``partition1`` -- an ordered partition, as a list of lists
    - ``ordering2`` -- represents a partition of the points of ``S2``, as a relabeling of ``partition1``
    - ``n`` -- the degree

    ::

        all_children_are_equivalent -- Python function of "signature":
            bool all_children_are_equivalent(PythonPartitionStack, object)
        refine_and_return_invariant -- Python function of "signature":
            int refine_and_return_invariant(PythonPartitionStack, object, list)
        compare_structures -- Python function of "signature":
            int compare_structures(list, list, object, object)
        (see double_coset.pyx for more documentation)

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_python import double_coset_python
        sage: def acae(a, b):
        ....:     return 0
        sage: def rari(a, b, c):
        ....:     return 0
        sage: def cs(a, b, c, d, e):
        ....:     return 0
        sage: double_coset_python(None, None, [[0,1,2,3],[4,5]], [2,3,1,5,0,4], 6, acae, rari, cs)
        [1, 2, 3, 5, 0, 4]

        sage: def compare_lists(p1, p2, l1, l2, deg):
        ....:     for i in range(len(l1)):
        ....:         a1 = l1[p1[i]]
        ....:         a2 = l2[p2[i]]
        ....:         if a1 < a2: return -1
        ....:         if a1 > a2: return 1
        ....:     return 0

        sage: double_coset_python([0,0,1], [1,0,0], [[0,1,2]], [0,1,2], 3, acae, rari, compare_lists)
        [1, 2, 0]
    """
    obj_wrapper1 = PythonObjectWrapper(S1, all_children_are_equivalent, refine_and_return_invariant, compare_structures, n)
    obj_wrapper2 = PythonObjectWrapper(S2, all_children_are_equivalent, refine_and_return_invariant, compare_structures, n)

    cdef PartitionStack *part = PS_from_list(partition1)
    cdef int *ordering = <int *> sig_malloc(n * sizeof(int))
    cdef int *output = <int *> sig_malloc(n * sizeof(int))
    if part is NULL or ordering is NULL or output is NULL:
        PS_dealloc(part)
        sig_free(ordering)
        sig_free(output)
        raise MemoryError
    for i from 0 <= i < n:
        ordering[i] = ordering2[i]

    cdef bint isomorphic = double_coset(<void *> obj_wrapper1, <void *> obj_wrapper2,
        part, ordering, n,
        &all_children_are_equivalent_python,
        &refine_and_return_invariant_python,
        &compare_structures_python, NULL, NULL, output)

    PS_dealloc(part)
    sig_free(ordering)
    if isomorphic:
        output_py = [output[i] for i from 0 <= i < n]
    else:
        output_py = False
    sig_free(output)
    return output_py
