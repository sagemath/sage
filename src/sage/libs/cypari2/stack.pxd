from .types cimport GEN, pari_sp
from .gen cimport gen

cdef void clear_stack()
cdef GEN deepcopy_to_python_heap(GEN x, pari_sp* address)
cdef gen new_gen(GEN x)
cdef gen new_gen_noclear(GEN x)
