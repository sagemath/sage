# distutils: libraries = primecount
# distutils: language = c++

from libc.stdint cimport int64_t
from libcpp.string cimport string as cppstring

cdef extern from "primecount.hpp" namespace "primecount":
    int64_t pi(int64_t x)

    cppstring pi(const cppstring& x)

    int64_t nth_prime(int64_t n)

    int64_t phi(int64_t x, int64_t a)

    void set_num_threads(int num_threads)
    int get_num_threads()

    cppstring get_max_x()
    cppstring get_max_x(double alpha)

    cppstring primecount_version()
