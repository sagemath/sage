#*****************************************************************************
#      Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/cdefs.pxi'
include 'sage/ext/stdsage.pxi'
from sage.misc.bitset cimport *
from libc.stdlib cimport rand

cdef extern from "flint/ulong_extras.h":
    int n_is_prime(unsigned long n)

cdef enum:
    # The following is for the automorphism group computation, says what the
    # length of fixed point and minimal cell representative arrays should be,
    # which are used to prune the search tree under known symmetries.
    # Increasing this may make certain automorphism_group_and_canonical_label
    # computations faster for incredibly large groups.
    len_of_fp_and_mcr = 100

cdef struct OrbitPartition:
    # Disjoint set data structure for representing the orbits of the generators
    # so far found. Also keeps track of the minimum elements of the cells and
    # their sizes.
    int degree
    int num_cells
    int *parent
    int *rank
    int *mcr # minimum cell representatives - only valid at the root of a cell
    int *size # also only valid at the root of a cell

cdef struct PartitionStack:
    # Representation of a node of the search tree. A sequence of partitions of
    # length depth + 1, each of which is finer than the last. Partition k is
    # represented as PS.entries in order, broken up immediately after each
    # entry of levels which is at most k.
    int *entries
    int *levels
    int depth
    int degree

cdef struct StabilizerChain:
    # A representation of a permutation group acting on 0, 1, ..., degree-1.
    int degree
    int base_size

    int *orbit_sizes
    int *num_gens      # dimension of generator cube on each level
    int *array_size    # size of space to hold generators on each level (number of permutations)

    int **base_orbits #
    int **parents     # three n*n squares, orbits and tree structures
    int **labels      #

    int **generators   # generators for each level,
    int **gen_inverses # and their inverses

    bitset_s gen_used
    bitset_s gen_is_id
    int *perm_scratch
    OrbitPartition *OP_scratch
