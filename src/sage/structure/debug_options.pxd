# sage_setup: distribution = sagemath-objects
cdef class DebugOptions_class:
    cdef public bint unique_parent_warnings
    cdef public bint refine_category_hash_check
    cdef public bint test_category_graph

cdef DebugOptions_class debug
