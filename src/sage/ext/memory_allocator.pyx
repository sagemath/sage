from cysignals.memory cimport *

cdef extern from *:
    int unlikely(int) nogil  # Defined by Cython


cdef class MemoryAllocator:
    r"""
    An object for memory allocation, whose resources are freed upon
    ``__dealloc__``.

    EXAMPLES::

        sage: cython(
        ....: '''
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: cdef MemoryAllocator mem = MemoryAllocator()
        ....: for n in range(100):
        ....:     mem.malloc(n)
        ....:     mem.calloc(n, n)
        ....:     mem.allocarray(n, n)
        ....: ''')
    """
    def __cinit__(self):
        """
        EXAMPLES::

            sage: cython(
            ....: '''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: cdef MemoryAllocator mem = MemoryAllocator.__new__(MemoryAllocator)
            ....: mem.malloc(10000)
            ....: print(mem.n)
            ....: print(mem.size)
            ....: ''')
            1
            16
        """
        self.n = 0
        self.size = 16
        self.pointers = self.static_pointers

    cdef int resize(self, size_t new_size) except -1:
        r"""
        Resize the list of pointers to contain ``new_size`` elements.

        It is required that ``new_size`` is at least ``self.n``, but
        this condition is not checked.
        """
        cdef size_t i
        if self.pointers == self.static_pointers:
            # Case 1: allocate pointers for the first time
            self.pointers = <void **>check_allocarray(new_size, sizeof(void*))
            for i in range(self.n):
                self.pointers[i] = self.static_pointers[i]
        else:
            # Case 2: resize pointers
            self.pointers = <void **>check_reallocarray(self.pointers, new_size, sizeof(void*))
        self.size = new_size

    cdef inline int enlarge_if_needed(self) except -1:
        r"""
        Enlarge the list of pointers if needed such that there is at
        least one free entry.
        """
        if unlikely(self.n >= self.size):
            return self.resize(self.size * 2)

    cdef void * malloc(self, size_t size) except? NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        self.enlarge_if_needed()
        cdef void * val = check_malloc(size)
        self.pointers[self.n] = val
        self.n += 1
        return val

    cdef void * calloc(self, size_t nmemb, size_t size) except? NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        self.enlarge_if_needed()
        cdef void * val = check_calloc(nmemb, size)
        self.pointers[self.n] = val
        self.n += 1
        return val

    cdef void * allocarray(self, size_t nmemb, size_t size) except? NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        self.enlarge_if_needed()
        cdef void * val = check_allocarray(nmemb, size)
        self.pointers[self.n] = val
        self.n += 1
        return val

    def __dealloc__(self):
        r"""
        Free the allocated resources

        EXAMPLES::

            sage: from sage.ext.memory_allocator import MemoryAllocator
            sage: _ = MemoryAllocator()
        """
        cdef size_t i
        for i in range(self.n):
            sig_free(self.pointers[i])
        if self.pointers != self.static_pointers:
            sig_free(self.pointers)
