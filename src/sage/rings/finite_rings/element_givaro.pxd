from sage.structure.element cimport Element, RingElement, ModuleElement
from sage.rings.finite_rings.element_base cimport FinitePolyExtElement

from sage.structure.parent  cimport Parent
from sage.structure.sage_object cimport SageObject

cdef extern from "givaro/givconfig.h":
    pass

cdef extern from "givaro/givrandom.h":
    ctypedef struct GivRandom "Givaro::GivRandom":
        pass

    GivRandom GivRandomSeeded  "Givaro::GivRandom"(unsigned long seed)

cdef extern from "givaro/givgfq.h":
    ctypedef struct intvec "std::vector<unsigned int>":
        void (* push_back)(int elem)

    ctypedef struct constintvec "const std::vector<unsigned int>"

    intvec intvec_factory "std::vector<unsigned int>"(int len)

cdef extern from "givaro/givgfq.h":

    ctypedef struct GivaroGfq "Givaro::GFqDom<int>":
        #attributes
        unsigned int one
        unsigned int zero

        # methods
        int (* mul)(int r, int a, int b)
        int (* add)(int r, int a, int b)
        int (* sub)(int r, int a, int b)
        int (* div)(int r, int a, int b)
        int (* inv)(int r, int x)
        int (* neg)(int r, int x)
        int (* mulin)(int a, int b)
        unsigned int (* characteristic)()
        unsigned int (* cardinality)()
        int (* exponent)()
        int (* random)(GivRandom gen, int res)
        int (* initi "init")(int res, int e)
        int (* initd "init")(int res, double e)
        int (* indeterminate)()
        int (* convert)(int r, int p)
        int (* read)(int r, int p)
        int (* axpyin)(int r, int a, int x)
        int (* axpy)(int r, int a, int b, int c)
        int (* axmy)(int r, int a, int b, int c)
        int (* maxpy)(int r, int a, int b, int c)
        bint (* isZero)(int e)
        bint (* isOne)(int e)
        bint (* isunit)(int e)

    GivaroGfq *gfq_factorypk "new Givaro::GFqDom<int>" (unsigned int p, unsigned int k)
    GivaroGfq *gfq_factorypkp "new Givaro::GFqDom<int>" (unsigned int p, unsigned int k, intvec poly)
    GivaroGfq *gfq_factorycopy "new Givaro::GFqDom<int>"(GivaroGfq orig)
    GivaroGfq  gfq_deref "*"(GivaroGfq *orig)
    void delete "delete "(void *o)
    int gfq_element_factory "Givaro::GFqDom<int>::Element"()

cdef class FiniteField_givaroElement(FinitePolyExtElement) #forward declaration

cdef class Cache_givaro(SageObject):
    cdef GivaroGfq *objectptr # C++ object
    cdef public object _array
    cdef FiniteField_givaroElement _zero_element
    cdef FiniteField_givaroElement _one_element
    cdef public int repr
    cdef bint _has_array
    cdef bint _is_conway
    cdef Parent parent
    cdef gen_array(self)
    cpdef int exponent(self)
    cpdef int order_c(self)
    cpdef int characteristic(self)
    cpdef FiniteField_givaroElement gen(self)
    cpdef FiniteField_givaroElement element_from_data(self, e)
    cdef FiniteField_givaroElement _new_c(self, int value)
    cpdef int int_to_log(self, int i) except -1
    cpdef int log_to_int(self, int i) except -1

cdef class FiniteField_givaro_iterator:
    cdef int iterator
    cdef Cache_givaro _cache

cdef class FiniteField_givaroElement(FinitePolyExtElement):
    cdef int element
    cdef Cache_givaro _cache
    cdef object _multiplicative_order
    cdef FiniteField_givaroElement _new_c(self, int value)


cdef inline FiniteField_givaroElement make_FiniteField_givaroElement(Cache_givaro cache, int x)
