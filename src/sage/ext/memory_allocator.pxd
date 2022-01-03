cimport cython
from libc.stdint cimport uintptr_t


cdef extern from *:
    int unlikely(int) nogil  # defined by Cython


cdef inline void* align(void* ptr, size_t alignment):
    """
    Round up ``ptr`` to the nearest multiple of ``alignment``, which
    must be a power of 2
    """
    cdef uintptr_t x = <uintptr_t>ptr
    x = (x + alignment - 1) & ~(alignment - 1)
    return <void*>x


@cython.final
cdef class MemoryAllocator:
    cdef size_t n
    cdef size_t size
    cdef void** pointers
    cdef void* static_pointers[16]  # If n <= 16, store pointers here

    cdef void* malloc(self, size_t size) except? NULL
    cdef void* calloc(self, size_t nmemb, size_t size) except? NULL
    cdef void* allocarray(self, size_t nmemb, size_t size) except? NULL
    cdef void* realloc(self, void* ptr, size_t size) except? NULL
    cdef void* reallocarray(self, void* ptr, size_t nmemb,
                            size_t size) except? NULL

    cdef int resize(self, size_t new_size) except -1
    cdef void** find_pointer(self, void* ptr) except NULL

    cdef inline int enlarge_if_needed(self) except -1:
        r"""
        Enlarge the list of pointers if needed such that there is at
        least one free entry.
        """
        if unlikely(self.n >= self.size):
            return self.resize(self.size * 2)

    cdef inline void* aligned_malloc(self, size_t alignment,
                                     size_t size) except? NULL:
        r"""
        Returns new aligned pointer. Stores it to be automatically freed later.

        Alignment must be a power of two.

        .. NOTE::

            If you want to allocate multiple (small) aligned arrays with the
            same alignment and really want to be memory efficient, you can
            allocate one large aligned array instead.

        TESTS::

            sage: cython('''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: cdef MemoryAllocator mem = MemoryAllocator()
            ....: cdef void* ptr
            ....: for i in range(12):
            ....:     ptr = mem.aligned_malloc(2**i, 4048)
            ....:     assert <size_t> ptr == (<size_t> ptr) & ~(2**i-1)
            ....: ''')
            doctest:...: DeprecationWarning: this class is deprecated; use the class from the python package `memory_allocator`
            See https://trac.sagemath.org/31591 for details.
        """
        cdef size_t extra = alignment - 1
        return align(self.malloc(size + extra), alignment)

    cdef inline void* aligned_calloc(self, size_t alignment, size_t nmemb,
                                     size_t size) except? NULL:
        r"""
        Returns new aligned pointer. Stores it to be automatically freed later.

        Alignment must be a power of two.

        .. NOTE::

            If you want to allocate multiple (small) aligned arrays with the
            same alignment and really want to be memory efficient, you can
            allocate one large aligned array instead.

        TESTS::

            sage: cython('''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: def foo():
            ....:     cdef MemoryAllocator mem = MemoryAllocator()
            ....:     cdef void* ptr
            ....:     for i in range(12):
            ....:         ptr = mem.aligned_calloc(2**i, i, 2**i)
            ....:         assert <size_t> ptr == (<size_t> ptr) & ~(2**i-1)
            ....: ''')
            sage: foo()
            doctest:...: DeprecationWarning: this class is deprecated; use the class from the python package `memory_allocator`
            See https://trac.sagemath.org/31591 for details.
        """
        # Find extra such that (nmemb + extra) * size >= nmemb * size + alignment - 1
        # ⇔ extra * size >= alignment - 1
        # ⇔ extra >= ceil( (alignment - 1) / size)
        # ⇔ extra >= (alignment - 1 + size - 1) // size
        cdef size_t extra = (alignment + size - 2) // size
        return align(self.calloc(nmemb + extra, size), alignment)

    cdef inline void* aligned_allocarray(self, size_t alignment, size_t nmemb,
                                         size_t size) except? NULL:
        r"""
        Returns new aligned pointer. Stores it to be automatically freed later.

        Alignment must be a power of two.

        .. NOTE::

            If you want to allocate multiple (small) aligned arrays with the
            same alignment and really want to be memory efficient, you can
            allocate one large aligned array instead.

        TESTS::

            sage: cython('''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: def foo():
            ....:     cdef MemoryAllocator mem = MemoryAllocator()
            ....:     cdef void* ptr
            ....:     for i in range(12):
            ....:         ptr = mem.aligned_allocarray(2**i, i, 2**i)
            ....:         assert <size_t> ptr == (<size_t> ptr) & ~(2**i-1)
            ....: ''')
            sage: foo()  # random  # might raise deprecation warning
            sage: foo()
        """
        # Find extra such that (nmemb + extra) * size >= nmemb * size + alignment - 1
        # ⇔ extra * size >= alignment - 1
        # ⇔ extra >= ceil( (alignment - 1) / size)
        # ⇔ extra >= (alignment - 1 + size - 1) // size
        cdef size_t extra = (alignment + size - 2) // size
        return align(self.allocarray(nmemb + extra, size), alignment)
