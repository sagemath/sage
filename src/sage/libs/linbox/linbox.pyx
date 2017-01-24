# distutils: extra_compile_args = LINBOX_CFLAGS
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: language = c++

r"""
Linbox interface
"""

## NOTE: The sig_on()/sig_off() stuff can't go in here -- it has to be in the
## code that calls these functions.  Otherwise strangely objects get left
## in an incorrect state.

from sage.libs.gmp.mpz cimport *
from sage.rings.integer cimport Integer
from sage.misc.misc import verbose, get_verbose

##########################################################################
## Sparse matrices modulo p.
##########################################################################

cdef extern from "linbox/linbox-sage.h":
    ctypedef struct vector_uint "std::vector<unsigned int>":
        void (*push_back)(unsigned int)
        int (*get "operator[]") (size_t i)
        int (*size)()

    int linbox_modn_sparse_matrix_rank(unsigned long modulus, size_t nrows, size_t ncols, void* rows, int reorder) #, int **pivots)
    vector_uint linbox_modn_sparse_matrix_solve(unsigned long modulus, size_t numrows, size_t numcols, void *a,  void *b, int method)

cdef class Linbox_modn_sparse:
    cdef set(self, int modulus, size_t nrows, size_t ncols, c_vector_modint *rows):
        self.rows = rows
        self.nrows = nrows
        self.ncols = ncols
        self.modulus = modulus

    cdef object rank(self, int gauss):
        #cdef int *_pivots
        cdef int r
        r = linbox_modn_sparse_matrix_rank(self.modulus, self.nrows, self.ncols, self.rows, gauss)

        #pivots = [_pivots[i] for i in range(r)]
        #free(_pivots)
        return r#, pivots

    cdef void solve(self, c_vector_modint **x, c_vector_modint *b, int method):
        """
        """
        cdef vector_uint X
        X = linbox_modn_sparse_matrix_solve(self.modulus, self.nrows, self.ncols, self.rows, b, method)

        for i from 0 <= i < X.size():
            set_entry(x[0], i, X.get(i))


##########################################################################
## Sparse matrices over ZZ
##########################################################################


##########################################################################
## Dense matrices over ZZ
##########################################################################

cdef extern from "linbox/linbox-sage.h":
    void linbox_integer_dense_minpoly(mpz_t* &minpoly, size_t &degree, size_t n, mpz_t** matrix)

    void linbox_integer_dense_charpoly(mpz_t* &charpoly, size_t &degree, size_t n, mpz_t** matrix)

    void linbox_integer_dense_delete_array(mpz_t* f)

    int linbox_integer_dense_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B, size_t A_nr, size_t A_nc, size_t B_nc)

    unsigned long linbox_integer_dense_rank(mpz_t** matrix, size_t nrows,
                                            size_t ncols)

    void linbox_integer_dense_det(mpz_t ans, mpz_t** matrix,
                             size_t nrows, size_t ncols)

    void linbox_integer_dense_smithform(mpz_t **v,
                                        mpz_t **matrix,
                                        size_t nrows, size_t ncols)

cdef class Linbox_integer_dense:
    cdef set(self, mpz_t** matrix, size_t nrows, size_t ncols):
        self.nrows = nrows
        self.ncols = ncols
        self.matrix = matrix

    def minpoly(self):
        """
        OUTPUT:

        coefficients of minpoly as a Python list
        """
        cdef mpz_t* poly
        cdef size_t degree
        verbose("using linbox poly comp")
        linbox_integer_dense_minpoly(poly, degree, self.nrows, self.matrix)
        verbose("computed poly -- now converting back to Sage")

        v = []
        cdef Integer k
        cdef size_t n
        for n from 0 <= n <= degree:
            k = Integer()
            mpz_set(k.value, poly[n])
            mpz_clear(poly[n])
            v.append(k)
        linbox_integer_dense_delete_array(poly)
        return v

    def charpoly(self):
        """
        OUTPUT:

        coefficients of charpoly or minpoly as a Python list
        """
        cdef mpz_t* poly
        cdef size_t degree
        verbose("using linbox poly comp")
        linbox_integer_dense_charpoly(poly, degree, self.nrows, self.matrix)
        verbose("computed poly -- now converting back to Sage")

        v = []
        cdef Integer k
        cdef size_t n
        for n from 0 <= n <= degree:
            k = Integer()
            mpz_set(k.value, poly[n])
            mpz_clear(poly[n])
            v.append(k)
        linbox_integer_dense_delete_array(poly)
        return v

    cdef matrix_matrix_multiply(self,
                                mpz_t **ans,
                                mpz_t **B,
                                size_t B_nr, size_t B_nc):
        cdef int e
        e = linbox_integer_dense_matrix_matrix_multiply(ans, self.matrix,  B, self.nrows, self.ncols, B_nc)
        if e:
            raise RuntimeError("error doing matrix matrix multiply over ZZ using linbox")


    cdef unsigned long rank(self) except -1:
        return linbox_integer_dense_rank(self.matrix, self.nrows, self.ncols)

    def det(self):
        """
        OUTPUT:

        determinant as a sage Integer
        """
        cdef Integer z
        z = Integer()
        linbox_integer_dense_det(z.value, self.matrix, self.nrows, self.ncols)
        return z

    def smithform(self):
        raise NotImplementedError
        #cdef mpz_t* v
        #linbox_integer_dense_smithform(&v, self.matrix, self.nrows, self.ncols)
        #z = []
        #cdef Integer k
        #cdef size_t n
        #for n from 0 <= n < self.ncols:
        #    k = Integer()
        #    mpz_set(k.value, v[n])
        #    mpz_clear(v[n])
        #    z.append(k)
        #linbox_integer_dense_delete_array(v)
        #return z

