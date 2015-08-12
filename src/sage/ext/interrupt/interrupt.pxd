# distutils: depends = INTERRUPT_DEPENDS
#
# NOTE: these functions are actually defined in "macros.h".
# However, we intentionally do not mention that file here, because
# every .pyx file using interrupt functions *must* also call
# import_sage__ext__interrupt__interrupt() which is done automatically
# if interrupt.pxi is included.
#
# A failure to include interrupt.pxi will therefore lead to a compiler
# error (since sig_on() will not be defined) or an ImportError due to a
# missing symbol but hopefully not an obscure segmentation fault.
#
cdef extern from *:
    int sig_on() nogil except 0
    int sig_str(char*) nogil except 0
    int sig_check() nogil except 0
    void sig_off() nogil
    void sig_retry() nogil  # Does not return
    void sig_error() nogil  # Does not return
    void sig_block() nogil
    void sig_unblock() nogil

    # Macros behaving exactly like sig_on, sig_str and sig_check but
    # which are *not* declared "except 0".  This is useful if some
    # low-level Cython code wants to do its own exception handling.
    int sig_on_no_except "sig_on"() nogil
    int sig_str_no_except "sig_str"(char*) nogil
    int sig_check_no_except "sig_check"() nogil

# This function does nothing, but it is declared cdef except *, so it
# can be used to make Cython check whether there is a pending exception
# (PyErr_Occurred() is non-NULL). To Cython, it will look like
# cython_check_exception() actually raised the exception.
cdef inline void cython_check_exception() nogil except *:
    pass


# Private stuff below, do not use directly
cdef extern from "interrupt/struct_signals.h":
    ctypedef struct sage_signals_t:
        int sig_on_count
        const char* s

cdef api:
    sage_signals_t _signals "_signals"
    void print_backtrace() nogil
    void _sig_on_interrupt_received() nogil
    void _sig_on_recover() nogil
    void _sig_off_warning(const char*, int) nogil
