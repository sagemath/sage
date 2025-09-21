r"""
Canonical augmentation

This module implements a general algorithm for generating isomorphism classes
of objects. The class of objects in question must be some kind of structure
which can be built up out of smaller objects by a process of augmentation,
and for which an automorphism is a permutation in `S_n` for some `n`. This
process consists of starting with a finite number of "seed objects" and
building up to more complicated objects by a sequence of "augmentations." It
should be noted that the word "canonical" in the term canonical augmentation
is used loosely. Given an object `X`, one must define a canonical parent
`M(X)`, which is essentially an arbitrary choice.

The class of objects in question must satisfy the assumptions made in the
module ``automorphism_group_canonical_label``, in particular the three
custom functions mentioned there must be implemented:

A. ``refine_and_return_invariant``:

    Signature:

    ``int refine_and_return_invariant(PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len)``

B. ``compare_structures``:

    Signature:

    ``int compare_structures(int *gamma_1, int *gamma_2, void *S1, void *S2, int degree)``

C. ``all_children_are_equivalent``:

    Signature:

    ``bint all_children_are_equivalent(PartitionStack *PS, void *S)``


In the following functions there is frequently a mem_err input. This is a
pointer to an integer which must be set to a nonzero value in case of an
allocation failure. Other functions have an int return value which serves the
same purpose. The idea is that if a memory error occurs, the canonical
generator should still be able to iterate over the objects already generated
before it terminates.

More details about these functions can be found in that module. In addition,
several other functions must be implemented, which will make use of the
following::

    ctypedef struct iterator:
        void *data
        void *(*next)(void *data, int *degree, int *mem_err)

The following functions must be implemented for each specific type of object to
be generated. Each function following which takes a ``mem_err`` variable as
input should make use of this variable.

D. ``generate_children``:

    Signature:

    ``int generate_children(void *S, aut_gp_and_can_lab *group, iterator *it)``

    This function receives a pointer to an iterator ``it``. The iterator
    has two fields: ``data`` and ``next``. The function ``generate_children``
    should set these two fields, returning 1 to indicate a memory error, or 0
    for no error.

    The function that ``next`` points to takes ``data`` as an argument, and
    should return a (``void *``) pointer to the next object to be iterated. It
    also takes a pointer to an int, and must update that int to reflect the
    degree of each generated object. The objects to be iterated over should
    satisfy the property that if `\gamma` is an automorphism of the parent
    object `S`, then for any two child objects `C_1, C_2` given by the iterator,
    it is not the case that `\gamma(C_1) = C_2`, where in the latter `\gamma` is
    appropriately extended if necessary to operate on `C_1` and `C_2`. It is
    essential for this iterator to handle its own ``data``. If the ``next``
    function is called and no suitable object is yielded, a NULL pointer
    indicates a termination of the iteration. At this point, the data pointed to
    by the ``data`` variable should be cleared by the ``next`` function, because
    the iterator struct itself will be deallocated.

    The ``next`` function must check ``mem_err[0]`` before proceeding. If it is
    nonzero then the function should deallocate the iterator right away and
    return NULL to end the iteration. This ensures that the canonical
    augmentation software will finish iterating over the objects found before
    finishing, and the ``mem_err`` attribute of the ``canonical_generator_data``
    will reflect this.

    The objects which the iterator generates can be thought of as augmentations,
    which the following function must turn into objects.

E. ``apply_augmentation``:

    Signature:

    ``void *apply_augmentation(void *parent, void *aug, void *child, int *degree, bint *mem_err)``

    This function takes the ``parent``, applies the augmentation ``aug`` and
    returns a pointer to the corresponding child object (freeing aug if
    necessary). Should also update degree[0] to be the degree of the new child.

F. ``free_object``:

    Signature:

    ``void free_object(void *child)``

    This function is a simple deallocation function for children which are not
    canonically generated, and therefore rejected in the canonical augmentation
    process. They should deallocate the contents of ``child``.

G. ``free_iter_data``:

    Signature:

    ``void free_iter_data(void *data)``

    This function deallocates the data part of the iterator which is set up by
    ``generate_children``.

H. ``free_aug``:

    Signature:

    ``void free_aug(void *aug)``

    This function frees an augmentation as generated by the iterator returned
    by ``generate_children``.

I. ``canonical_parent``:

    Signature:

    ``void *canonical_parent(void *child, void *parent, int *permutation, int *degree, bint *mem_err)``

    Apply the ``permutation`` to the ``child``, determine an arbitrary but fixed
    parent, apply the inverse of ``permutation`` to that parent, and return the
    resulting object. Must also set the integer ``degree`` points to the
    degree of the returned object.

.. NOTE::

    It is a good idea to try to implement an augmentation scheme where the
    degree of objects on each level of the augmentation tree is constant. The
    iteration will be more efficient in this case, as the relevant work spaces
    will never need to be reallocated. Otherwise, one should at least strive to
    iterate over augmentations in such a way that all children of the same degree
    are given in the same segment of iteration.

EXAMPLES::

    sage: import sage.groups.perm_gps.partn_ref.canonical_augmentation

REFERENCE:

- [1] McKay, Brendan D. Isomorph-free exhaustive generation. J Algorithms,
  Vol. 26 (1998), pp. 306-324.
"""

#*****************************************************************************
#       Copyright (C) 2010 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.memory cimport sig_malloc, sig_free

from sage.groups.perm_gps.partn_ref.data_structures cimport*


cdef void *canonical_generator_next(void *can_gen_data, int *degree, bint *mem_err) noexcept:
    r"""
    This function is part of the iterator struct which will iterate over
    objects. Return value of ``NULL`` indicates termination.
    """
    # degree ignored!
    cdef canonical_generator_data *cgd = <canonical_generator_data *> can_gen_data
    cdef iterator *cur_iter
    cdef void *next_candidate
    cdef void *parent_cand
    cdef void *aug
    cdef int next_cand_deg, parent_cand_deg
    cdef PartitionStack *part
    cdef bint augmentation_is_canonical

    if cgd.level == 0:
        if cgd.mem_err:
            mem_err[0] = 1
        if cgd.dealloc:
            deallocate_cgd(cgd)
        return NULL

    while cgd.level < cgd.max_level:
        cur_iter = cgd.iterator_stack + cgd.level - 1
        # getting next candidate at level cgd.level
        aug = cur_iter.next( cur_iter.data, &cgd.degree_stack[cgd.level], &cgd.mem_err )
        if aug is NULL:
            # cur_iter has run out: backtrack by one, return the parent
            cgd.level -= 1
            return cgd.object_stack[cgd.level]
        else:
            next_candidate = cgd.apply_augmentation(cgd.object_stack[cgd.level-1],
                aug, cgd.object_stack[cgd.level], &cgd.degree_stack[cgd.level],
                &cgd.mem_err)
            cgd.object_stack[cgd.level] = next_candidate
            if cgd.mem_err:
                continue
            next_cand_deg = cgd.degree_stack[cgd.level]
            if cgd.agcl_work_spaces[cgd.level] is NULL:
                # allocate a work space if it hasn't been allocated already
                cgd.agcl_work_spaces[cgd.level] = allocate_agcl_work_space(next_cand_deg)
                cgd.ps_stack[cgd.level] = PS_new(next_cand_deg, 0)
            elif next_cand_deg != cgd.agcl_work_spaces[cgd.level].degree:
                # in case the degree has changed, we need to reallocate the work
                # space: this is why one should do augmentations in order of degree
                deallocate_agcl_work_space(cgd.agcl_work_spaces[cgd.level])
                deallocate_agcl_output(cgd.aut_gp_stack[cgd.level])
                cgd.aut_gp_stack[cgd.level] = NULL
                PS_dealloc(cgd.ps_stack[cgd.level])
                cgd.agcl_work_spaces[cgd.level] = allocate_agcl_work_space(next_cand_deg)
                cgd.ps_stack[cgd.level] = PS_new(next_cand_deg, 0)
            if cgd.agcl_work_spaces[cgd.level] is NULL or cgd.ps_stack[cgd.level] is NULL:
                cgd.mem_err = 1
                continue
            # see if next_candidate is canonically augmented
            part = cgd.ps_stack[cgd.level]
            PS_unit_partition(part)
            try:
                cgd.aut_gp_stack[cgd.level] = get_aut_gp_and_can_lab(next_candidate,
                    part, next_cand_deg, cgd.all_children_are_equivalent,
                    cgd.refine_and_return_invariant, cgd.compare_structures,
                    1, cgd.group, cgd.agcl_work_spaces[cgd.level], cgd.aut_gp_stack[cgd.level])
            except MemoryError:
                cgd.mem_err = 1
                continue
            parent_cand = cgd.canonical_parent(next_candidate, cgd.parent_stack[cgd.level-1],
                cgd.aut_gp_stack[cgd.level].relabeling, &parent_cand_deg, &cgd.mem_err)
            if cgd.mem_err:
                continue
            if parent_cand_deg != cgd.degree_stack[cgd.level-1]:
                augmentation_is_canonical = 0
            else:
                if cgd.dc_work_spaces[cgd.level-1] is NULL:
                    # allocate a work space if it hasn't been allocated already
                    cgd.dc_work_spaces[cgd.level-1] = allocate_dc_work_space(next_cand_deg)
                elif next_cand_deg != cgd.dc_work_spaces[cgd.level-1].degree:
                    # in case the degree has changed, we need to reallocate the work
                    # space: this is why one should do augmentations in order of degree
                    deallocate_dc_work_space(cgd.dc_work_spaces[cgd.level-1])
                    cgd.dc_work_spaces[cgd.level-1] = allocate_dc_work_space(next_cand_deg)
                if cgd.dc_work_spaces[cgd.level-1] is NULL:
                    cgd.mem_err = 1
                    continue
                part = cgd.ps_stack[cgd.level]
                PS_unit_partition(part)
                try:
                    augmentation_is_canonical = double_coset(cgd.object_stack[cgd.level-1],
                        parent_cand, part, NULL, next_cand_deg,
                        cgd.all_children_are_equivalent,
                        cgd.refine_and_return_invariant,
                        cgd.compare_structures,
                        cgd.aut_gp_stack[cgd.level].group, cgd.dc_work_spaces[cgd.level-1], NULL)
                except MemoryError:
                    cgd.mem_err = 1
                    continue
            if augmentation_is_canonical:
                # the object is canonically augmented, so we add it to the chain
                if cgd.level + 1 != cgd.max_level:
                    cgd.mem_err |= cgd.generate_children(next_candidate,
                        cgd.aut_gp_stack[cgd.level], cgd.iterator_stack+cgd.level)
                    if cgd.mem_err:
                        continue
                cgd.level += 1

    if cgd.level == cgd.max_level:
        # we're at the end of the chain, so just give the object
        cgd.level -= 1
        return cgd.object_stack[cgd.level]

cdef canonical_generator_data *allocate_cgd(int max_depth, int degree) noexcept:
    r"""
    Allocate the data part of the canonical generation iterator struct.
    """
    cdef canonical_generator_data *cgd = <canonical_generator_data *> sig_malloc(sizeof(canonical_generator_data))
    cdef PartitionStack *part
    if cgd is NULL:
        sig_free(cgd)
        return NULL
    cgd.object_stack = <void **> sig_malloc(max_depth * sizeof(void *))
    cgd.degree_stack = <int *> sig_malloc(max_depth * sizeof(int))
    cgd.iterator_stack = <iterator *> sig_malloc(max_depth * sizeof(iterator))
    cgd.aut_gp_stack = <aut_gp_and_can_lab **> sig_malloc(max_depth * sizeof(aut_gp_and_can_lab *))
    cgd.agcl_work_spaces = <agcl_work_space **> sig_malloc(max_depth * sizeof(agcl_work_space *))
    cgd.dc_work_spaces = <dc_work_space **> sig_malloc(max_depth * sizeof(dc_work_space *))
    cgd.ps_stack = <PartitionStack **> sig_malloc(max_depth * sizeof(PartitionStack *))
    cgd.aug_stack = <void **> sig_malloc(max_depth * sizeof(void *))
    cgd.parent_stack = <void **> sig_malloc(max_depth * sizeof(void *))
    part = PS_new(degree, 1)
    cdef agcl_work_space *agclws    = allocate_agcl_work_space(degree)
    cdef aut_gp_and_can_lab *output = allocate_agcl_output(degree)
    if cgd.object_stack     is NULL or cgd.degree_stack   is NULL or \
       cgd.iterator_stack   is NULL or cgd.aut_gp_stack   is NULL or \
       cgd.agcl_work_spaces is NULL or cgd.dc_work_spaces is NULL or \
       cgd.ps_stack         is NULL or cgd.aug_stack      is NULL or \
       cgd.parent_stack     is NULL or agclws is NULL or output is NULL:
        sig_free(cgd.object_stack)
        sig_free(cgd.degree_stack)
        sig_free(cgd.iterator_stack)
        sig_free(cgd.aut_gp_stack)
        sig_free(cgd.agcl_work_spaces)
        sig_free(cgd.dc_work_spaces)
        sig_free(cgd.ps_stack)
        sig_free(cgd.aug_stack)
        sig_free(cgd.parent_stack)
        sig_free(cgd)
        PS_dealloc(part)
        deallocate_agcl_work_space(agclws)
        deallocate_agcl_output(output)
        return NULL

    cdef int i
    cgd.allocd_levels = max_depth
    for i from 0 <= i < max_depth:
        cgd.agcl_work_spaces[i]    = NULL
        cgd.dc_work_spaces[i]      = NULL
        cgd.aut_gp_stack[i]        = NULL
        cgd.ps_stack[i]            = NULL
        cgd.aug_stack[i]           = NULL
        cgd.parent_stack[i]        = NULL
        cgd.object_stack[i]        = NULL
        cgd.iterator_stack[i].data = NULL
    cgd.agcl_work_spaces[0] = agclws
    cgd.aut_gp_stack[0]     = output
    cgd.ps_stack[0]         = part

    cgd.degree_stack[0] = degree
    return cgd


cdef void deallocate_cgd(canonical_generator_data *cgd) noexcept:
    r"""
    Deallocate the data part of the canonical generation iterator struct.
    """
    if cgd is NULL:
        return
    cdef int i
    for i from 0 <= i < cgd.allocd_levels:
        if cgd.agcl_work_spaces[i] is not NULL:
            deallocate_agcl_work_space(cgd.agcl_work_spaces[i])
        if cgd.ps_stack[i] is not NULL:
            PS_dealloc(cgd.ps_stack[i])
        if cgd.dc_work_spaces[i] is not NULL:
            deallocate_dc_work_space(cgd.dc_work_spaces[i])
        if cgd.aut_gp_stack[i] is not NULL:
            deallocate_agcl_output(cgd.aut_gp_stack[i])
        if cgd.object_stack[i] is not NULL:
            cgd.free_object(cgd.object_stack[i])
        if cgd.parent_stack[i] is not NULL:
            cgd.free_object(cgd.parent_stack[i])
        if cgd.aug_stack[i] is not NULL:
            cgd.free_aug(cgd.aug_stack[i])
        if cgd.iterator_stack[i].data is not NULL:
            cgd.free_iter_data(cgd.iterator_stack[i].data)
    sig_free(cgd.object_stack)
    sig_free(cgd.degree_stack)
    sig_free(cgd.iterator_stack)
    sig_free(cgd.aut_gp_stack)
    sig_free(cgd.agcl_work_spaces)
    sig_free(cgd.dc_work_spaces)
    sig_free(cgd.ps_stack)
    sig_free(cgd.aug_stack)
    sig_free(cgd.parent_stack)
    sig_free(cgd)


cdef iterator *setup_canonical_generator(int degree,
    bint (*all_children_are_equivalent)(PartitionStack *PS, void *S) noexcept,
    int (*refine_and_return_invariant)(PartitionStack *PS, void *S,
                                       int *cells_to_refine_by, int ctrb_len) noexcept,
    int (*compare_structures)(int *gamma_1, int *gamma_2, void *S1, void *S2,
                              int degree) noexcept,
    int (*generate_children)(void *, aut_gp_and_can_lab *, iterator *) noexcept,
    void *(*apply_augmentation)(void *, void *, void *, int *, bint *) noexcept,
    void (*free_object)(void *) noexcept,
    void (*free_iter_data)(void *) noexcept,
    void (*free_aug)(void *) noexcept,
    void *(*canonical_parent)(void *child, void *parent, int *permutation, int *degree, bint *mem_err) noexcept,
    int max_depth, bint reduce_children, iterator *cangen_prealloc) except NULL:
    """
    Canonical generation of isomorphism classes of objects.

    INPUT:

    - ``S`` -- pointer to the seed object

    - ``degree`` -- the degree of S

    - ``all_children_are_equivalent`` -- pointer to a function
        INPUT:
        - ``PS`` -- pointer to a partition stack
        - ``S`` -- pointer to the structure
        OUTPUT:
        bint; returns ``True`` if it can be determined that all refinements below
        the current one will result in an equivalent discrete partition
    - ``refine_and_return_invariant`` -- pointer to a function
        INPUT:
        - ``PS`` -- pointer to a partition stack
        - ``S`` -- pointer to the structure
        - ``alpha`` -- an array consisting of numbers, which indicate the starting
          positions of the cells to refine against (will likely be modified)
        OUTPUT:
        integer; returns an invariant under application of arbitrary permutations
    - ``compare_structures`` -- pointer to a function
        INPUT:
        - ``gamma_1``, ``gamma_2`` -- (list) permutations of the points of S1 and S2
        - ``S1``, ``S2`` -- pointers to the structures
        - ``degree`` -- degree of gamma_1 and 2
        OUTPUT:
        integer; 0 if gamma_1(S1) = gamma_2(S2), otherwise -1 or 1 (see docs for cmp),
        such that the set of all structures is well-ordered
    - ``generate_children`` -- pointer to a function
        INPUT:
        - ``S`` -- pointer to the structure
        - ``group`` -- pointer to an automorphism group (canonical relabeling is not guaranteed)
        - ``it`` -- preallocated iterator struct
        OUTPUT:
        iterator *; pointer to an iterator over inequivalent augmentations of S
    - ``apply_augmentation`` -- pointer to a function
        INPUT:
        - ``parent`` -- object to augment
        - ``aug`` -- the augmentation
        - ``child`` -- space to put the augmented object
        - ``degree`` -- pointer to an int, function should store the degree of the augmented object here
        - ``mem_err`` -- pointer where memory error can be reported
        OUTPUT: pointer to child
    - ``free_object`` -- pointer to a function
        INPUT:
        - ``child`` -- object to be freed
    - ``free_iter_data`` -- pointer to a function
        INPUT:
        - ``data`` -- data part of an iterator struct
    - ``free_aug`` -- pointer to a function
        INPUT:
        - ``aug`` -- augmentation to be freed
    - ``canonical_parent`` -- pointer to a function
        INPUT:
        - ``child`` -- pointer to the structure
        - ``parent`` -- space to store the canonical parent
        - ``permutation`` -- array representing a relabeling of the child
        - ``degree`` -- pointer to store the degree of the parent
        - ``mem_err`` -- pointer for indicating memory errors
        OUTPUT: pointer to the parent

    - ``max_depth`` -- maximum depth of augmentations to be made from the seed object S

    OUTPUT: a pointer to an iterator of objects
    """
    if max_depth <= 1:
        raise ValueError("maximum depth (%d) must be at least two" % max_depth)
    if reduce_children:
        raise NotImplementedError

    # Allocate memory for the arrays and check for failures:
    cdef iterator *canonical_generator
    cdef canonical_generator_data *cgd
    if cangen_prealloc is NULL:
        canonical_generator = <iterator *> sig_malloc(sizeof(iterator))
        cgd = allocate_cgd(max_depth, degree)
        if canonical_generator is NULL or cgd is NULL:
            sig_free(canonical_generator)
            deallocate_cgd(cgd)
            raise MemoryError
        cgd.dealloc = 1
    else:
        canonical_generator = cangen_prealloc
        cgd = <canonical_generator_data *> canonical_generator.data
        if cgd.degree_stack[0] != degree or cgd.allocd_levels < max_depth:
            deallocate_cgd(cgd)
            cgd = allocate_cgd(max_depth, degree)
            if cgd is NULL:
                raise MemoryError
        cgd.dealloc = 0
    canonical_generator.data = <void *> cgd
    canonical_generator.next = &canonical_generator_next

    cgd.max_level = max_depth
    cgd.reduce_children = reduce_children
    cgd.mem_err = 0

    cgd.all_children_are_equivalent = all_children_are_equivalent
    cgd.refine_and_return_invariant = refine_and_return_invariant
    cgd.compare_structures          = compare_structures

    cgd.generate_children  = generate_children
    cgd.apply_augmentation = apply_augmentation
    cgd.free_object        = free_object
    cgd.free_iter_data     = free_iter_data
    cgd.free_aug           = free_aug
    cgd.canonical_parent   = canonical_parent

    return canonical_generator

cdef iterator *start_canonical_generator(StabilizerChain *group, void *obj, int degree, iterator *canonical_generator) except NULL:
    r"""
    Given the containing group ``group`` and the seed object ``obj`` of degree
    ``degree``, initiate the canonical generator stored (and already allocated)
    at ``canonical_generator``.
    """
    cdef canonical_generator_data *cgd = <canonical_generator_data *> canonical_generator.data
    if obj is NULL:
        obj = cgd.object_stack[0]
    else:
        cgd.object_stack[0] = obj
    cgd.level = 1
    cgd.group = group
    PS_unit_partition(cgd.ps_stack[0])
    try:
        cgd.aut_gp_stack[0] = get_aut_gp_and_can_lab(obj, cgd.ps_stack[0], degree,
            cgd.all_children_are_equivalent,
            cgd.refine_and_return_invariant,
            cgd.compare_structures,
            0, group, cgd.agcl_work_spaces[0], cgd.aut_gp_stack[0])
    except MemoryError:
        cgd.mem_err = 1
    else:
        cgd.mem_err |= cgd.generate_children(obj, cgd.aut_gp_stack[0], cgd.iterator_stack)
    if cgd.mem_err:
        raise MemoryError

    return canonical_generator
