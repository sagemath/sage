from sage.rings.integer cimport Integer
from sage.libs.gmp.mpz cimport *

cdef extern from "ratpoints.h":
    long RATPOINTS_MAX_DEGREE
    long RATPOINTS_ARRAY_SIZE
    long RATPOINTS_DEFAULT_SP1
    long RATPOINTS_DEFAULT_SP2
    long RATPOINTS_DEFAULT_NUM_PRIMES
    long RATPOINTS_DEFAULT_MAX_FORBIDDEN
    long RATPOINTS_DEFAULT_STURM
    long RATPOINTS_NON_SQUAREFREE
    long RATPOINTS_BAD_ARGS

    # for args flags:
    long RATPOINTS_NO_CHECK # when set, do not check whether the surviving
                            # x-coordinates give rise to rational points
    long RATPOINTS_NO_Y # when set, only list x coordinates instead of actual points
    long RATPOINTS_NO_REVERSE # when set, do not modify the mpz_t array
    long RATPOINTS_NO_JACOBI # when set, prevent use of Jacobi symbol test
    long RATPOINTS_VERBOSE # when set, print some output on what ratpoints is doing
    # define RATPOINTS_FLAGS_INPUT_MASK \
    # (RATPOINTS_NO_CHECK | RATPOINTS_NO_Y | RATPOINTS_NO_REVERSE | \
    #  RATPOINTS_NO_JACOBI | RATPOINTS_VERBOSE)


    ctypedef struct ratpoints_interval:
        double low
        double up
    ctypedef struct ratpoints_args:
        mpz_t *cof
        long degree
        long height
        ratpoints_interval *domain
        long num_inter
        long b_low
        long b_high
        long sp1
        long sp2
        long array_size
        long sturm
        long num_primes
        long max_forbidden
        unsigned int flags
        # from here: private data
        # mpz_t *work
        # void *se_buffer
        # void *se_next
        # void *ba_buffer
        # void *ba_next
        # int *int_buffer
        # int *int_next
        # void *sieve_list
    long find_points(ratpoints_args*, int proc(long, long, mpz_t, void*, int*), void*)
    void find_points_init(ratpoints_args*)
    long find_points_work(ratpoints_args*, int proc(long, long, mpz_t, void*, int*), void*)
    void find_points_clear(ratpoints_args*)

ctypedef struct point_list:
    long *xes
    mpz_t *ys
    long *zs
    long array_size
    long num_points
    long max_num_points

ctypedef struct info_struct_exists_only:
    int verbose

cdef int ratpoints_mpz_exists_only(mpz_t *, long, int, bint) except -1




