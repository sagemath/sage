# sage_setup: distribution = sagemath-linbox
# distutils: extra_compile_args = GIVARO_CFLAGS
# distutils: include_dirs = GIVARO_INCDIR
# distutils: libraries = GIVARO_LIBRARIES FFLASFFPACK_LIBRARIES
# distutils: library_dirs = GIVARO_LIBDIR
# distutils: language = c++

from libc.stdint cimport uint32_t, uint64_t

from sage.libs.gmp.types cimport (mpz_t, mpz_srcptr, mpz_ptr,
                                  mpq_t, mpq_srcptr, mpq_ptr)

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream& write(const char*, int) except +

cdef extern from "gmp++/gmp++.h" namespace "Givaro":
    cdef cppclass Integer:
        Integer()
        Integer(int32_t)
        Integer(int64_t)
        Integer(uint32_t)
        Integer(uint64_t)
        Integer(Integer&)

        mpz_ptr get_mpz()
        mpz_srcptr get_mpz_const()

cdef extern from "givaro/givcategory.h" namespace "Givaro":
    cdef cppclass Sporadic:
        pass
    cdef cppclass Dense:
        pass
    cdef cppclass Sparse:
        pass

cdef extern from "givaro/zring.h":
    ## template<class _Element> class ZRing
    cdef cppclass ZRing "Givaro::ZRing<Givaro::Integer>":
        ctypedef Integer Element
        Element zero
        Element one
        Element mone

cdef extern from "givaro/modular-integral.h":
    cdef cppclass Modular_uint64 "Givaro::Modular<uint64_t>":
        ctypedef uint64_t Element
        Modular_uint64(int modulus)

        Element init(Element res, int v)
        Element inv(Element x, Element y)
        Element neg(Element x, Element y)
        Element mul(Element r, Element x, Element y)
        Element mulin(Element x, Element y)
        Element addin(Element x, Element y)
        Element invin(Element y)
        Element negin(Element y)
        int characteristic(int c)
        bint isZero(Element x)

        ostream& write(ostream&)

cdef extern from "givaro/modular-floating.h":
    cdef cppclass Modular_double "Givaro::Modular<double>":
        ctypedef double Element
        Modular_double(int modulus)

        Element init(Element res, int v)
        Element init(Element res, double v)
        Element inv(Element x, Element y)
        Element neg(Element x, Element y)
        Element mul(Element r, Element x, Element y)
        Element mulin(Element x, Element y)
        Element addin(Element x, Element y)
        Element invin(Element y)
        Element negin(Element y)
        int characteristic(int c)
        bint isZero(Element x)

        ostream& write(ostream&)

    cdef cppclass Modular_float "Givaro::Modular<float>":
        ctypedef float Element
        Modular_float(int modulus)

        Element init(Element res, int v)
        Element init(Element res, float v)
        Element inv(Element x, Element y)
        Element neg(Element x, Element y)
        Element mul(Element r, Element x, Element y)
        Element mulin(Element x, Element y)
        Element addin(Element x, Element y)
        Element invin(Element y)
        Element negin(Element y)
        int characteristic(int c)
        bint isZero(Element x)

        ostream& write(ostream&)

cdef extern from "givaro/givpoly1.h" namespace "Givaro":
    ## template < typename T, typename A=std::allocator<T> >
    ## class givvector : public __GIV_STANDARD_VECTOR<T,A>
    cdef cppclass givvector [T,ALLOCATOR=*]:
        T& operator[](size_t i)
        size_t size()

    ## template<class Domain, class StorageTag=Dense> class Poly1Dom
    cdef cppclass Poly1Dom[Domain,StorageClass=*]:
        Poly1Dom(Domain&)
