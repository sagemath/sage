# sage_setup: distribution = sagemath-combinat
from sage.structure.sage_object cimport SageObject

cdef class RiggedPartition(SageObject):
    cdef public list _list
    cdef public list vacancy_numbers
    cdef public list rigging
    cdef long _hash

    cpdef get_num_cells_to_column(self, int end_column, t=*)
    cpdef insert_cell(self, int max_width)
    cpdef remove_cell(self, row, int num_cells=*)

cdef class RiggedPartitionTypeB(RiggedPartition):
    pass

