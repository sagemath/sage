# Declare system calls related to signal handling

cdef extern from "<stdlib.h>":
    void abort() nogil

cdef extern from "<signal.h>":
    ctypedef void *sigset_t
    # Renaming of this struct is necessary because Cython folds the
    # "struct" namespace into the normal namespace.
    struct Sigaction "sigaction":
        void (*sa_handler)(int)
        sigset_t sa_mask
        int sa_flags

    int SIG_BLOCK, SIG_UNBLOCK, SIG_SETMASK
    int SIGHUP, SIGINT, SIGQUIT, SIGILL, SIGABRT, SIGFPE, SIGKILL
    int SIGSEGV, SIGPIPE, SIGALRM, SIGTERM, SIGBUS
    int signal_raise "raise"(int signum)
    int sigprocmask(int how, sigset_t *set, sigset_t *oldset)
    int sigemptyset(sigset_t *set)
    int sigaddset(sigset_t *set, int signum)
    int sigaction(int signum, Sigaction *act, Sigaction *oldact)

cdef extern from "<sys/time.h>":
    struct timespec:
        long tv_sec
        long tv_nsec

cdef extern from "<sys/select.h>":
    ctypedef void *fd_set
    void FD_CLR(int fd, fd_set *set)
    bint FD_ISSET(int fd, fd_set *set)
    void FD_SET(int fd, fd_set *set)
    void FD_ZERO(fd_set *set)
    int FD_SETSIZE

    int pselect(int nfds, fd_set *readfds, fd_set *writefds,
                fd_set *exceptfds, timespec *timeout,
                sigset_t *sigmask)
