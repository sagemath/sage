cimport cython

@cython.final
cdef class MemoryAllocator:
    cdef size_t n
    cdef size_t size
    cdef void ** pointers
    cdef void * static_pointers[16]  # If n <= 16, store pointers here

    cdef void * malloc(self, size_t size) except? NULL
    cdef void * calloc(self, size_t nmemb, size_t size) except? NULL
    cdef void * allocarray(self, size_t nmemb, size_t size) except? NULL
    cdef int resize(self, size_t new_size) except -1
    cdef inline int enlarge_if_needed(self) except -1
