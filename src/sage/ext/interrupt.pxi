#
# See c_lib/include/interrupt.h
#
cdef extern from 'interrupt.h':
    int sig_on() nogil except 0
    int sig_str(char*) nogil except 0
    int sig_check() nogil except 0
    int sig_on_no_except() nogil
    int sig_str_no_except(char*) nogil
    int sig_check_no_except() nogil
    void sig_off() nogil
    void sig_retry() nogil  # Does not return
    void sig_error() nogil  # Does not return
    void sig_block() nogil
    void sig_unblock() nogil
    void set_sage_signal_handler_message(char* s) nogil
    void cython_check_exception() nogil except *

    ctypedef struct sage_signals_t:
        int sig_on_count
        int interrupt_received
        int inside_signal_handler
        int block_sigint
        char* s
        int (*raise_exception)(int sig, const char* msg) except 0

    sage_signals_t _signals
