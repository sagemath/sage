/*
Interrupt and signal handling for Sage.

For documentation about how to use these, see the Developer's Guide.

This code distinguishes between two kinds of signals:

(1) interrupt-like signals: SIGINT, SIGALRM, SIGHUP.  The word
"interrupt" refers to any of these signals.  These need not be handled
immediately, we might handle them at a suitable later time, outside of
sig_block() and with the Python GIL acquired.  SIGINT raises a
KeyboardInterrupt (as usual in Python), SIGALRM raises AlarmInterrupt
(a custom exception inheriting from KeyboardInterrupt), while SIGHUP
raises SystemExit, causing Python to exit.  The latter signal also
redirects stdin from /dev/null, to cause interactive sessions to exit.

(2) critical signals: SIGQUIT, SIGILL, SIGABRT, SIGFPE, SIGBUS, SIGSEGV.
These are critical because they cannot be ignored.  If they happen
outside of sig_on(), we can only exit Sage with the dreaded
"unhandled SIG..." message.  Inside of sig_on(), they can be handled
and raise various exceptions (see sage/ext/c_lib.pyx).  SIGQUIT will
never be handled and always causes Sage to exit.


AUTHORS:

- William Stein, Martin Albrecht (2006): initial version

- Jeroen Demeyer (2010-10-03): almost complete rewrite (#9678)

- Jeroen Demeyer (2013-01-11): handle SIGHUP also (#13908)

- Jeroen Demeyer (2013-01-28): handle SIGQUIT also (#14029)

- Jeroen Demeyer (2013-05-13): handle SIGALRM also (#13311)

*/

/*****************************************************************************
 *       Copyright (C) 2006 William Stein <wstein@gmail.com>
 *                     2006 Martin Albrecht <malb@informatik.uni-bremen.de>
 *                     2010-2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

/* Whether or not to compile debug routines for the interrupt handling
 * code (0: disable, 1: enable).  Enabling will make the code slower.
 * The debug level itself needs to be set in c_lib/src/interrupt.c */
#define ENABLE_DEBUG_INTERRUPT 0


#ifndef C_LIB_INCLUDE_INTERRUPT_H
#define C_LIB_INCLUDE_INTERRUPT_H
#include <Python.h>
#include <setjmp.h>
#include <signal.h>

#ifdef __cplusplus
extern "C" {
#endif


/* Declare likely() and unlikely() as in Cython */
#ifdef __GNUC__
/* Test for GCC > 2.95 */
#if __GNUC__ > 2 || (__GNUC__ == 2 && (__GNUC_MINOR__ > 95))
#define HAVE_BUILTIN_EXPECT 1
#endif
#endif

#if HAVE_BUILTIN_EXPECT
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define likely(x)   (x)
#define unlikely(x) (x)
#endif

/* Interrupt debug level */
#if ENABLE_DEBUG_INTERRUPT
extern int sage_interrupt_debug_level;
#endif


/* Print a C backtrace if supported by libc */
void print_backtrace(void);

/* Print a message s and kill ourselves with signal sig */
void sigdie(int sig, const char* s);


/*
 * The signal handlers for Sage, one for interrupt-like signals
 * (SIGINT, SIGHUP) and one for critical signals like SIGSEGV.
 *
 * Inside sig_on() (i.e. when _signals.sig_on_count is positive), these
 * handlers raise an exception and jump back to sig_on().
 * Outside of sig_on(), sage_interrupt_handler() sets Python's
 * interrupt flag using PyErr_SetInterrupt(); sage_signal_handler()
 * terminates Sage.
 */
void sage_interrupt_handler(int sig);
void sage_signal_handler(int sig);

/*
 * Setup the signal handlers. It is safe to call this more than once.
 *
 * We do not handle SIGALRM since there is code to deal with
 * alarms in sage/misc/misc.py
 */
void setup_sage_signal_handler(void);


/**********************************************************************
 * SAGE_SIGNALS_T STRUCTURE                                           *
 **********************************************************************/

/* All the state of the signal handler is in this struct. */
struct sage_signals_t
{
    /* Reference counter for sig_on().
     * If this is strictly positive, we are inside a sig_on(). */
    volatile sig_atomic_t sig_on_count;

    /* If this is nonzero, it is a signal number of a non-critical
     * signal (e.g. SIGINT) which happened during a time when it could
     * not be handled.  This may be set when an interrupt occurs either
     * outside of sig_on() or inside sig_block().  To avoid race
     * conditions, this value may only be changed when all
     * interrupt-like signals are masked. */
    volatile sig_atomic_t interrupt_received;

    /* Are we currently handling a signal inside sage_signal_handler()?
     * This is set to 1 on entry in sage_signal_handler (not in
     * sage_interrupt_handler) and 0 in _sig_on_postjmp.  This is
     * needed to check for signals raised within the signal handler. */
    volatile sig_atomic_t inside_signal_handler;

    /* Non-zero if we currently are in a function such as malloc()
     * which blocks interrupts, zero normally.
     * See sig_block(), sig_unblock(). */
    volatile sig_atomic_t block_sigint;

    /* A jump buffer holding where to siglongjmp() after a signal has
     * been received. This is set by sig_on(). */
    sigjmp_buf env;

    /* External Cython function which actually raises the appropriate
     * exception depending on the signal number. Must be set
     * immediately after calling setup_sage_signal_handler(). */
    int (*raise_exception)(int sig, const char* msg);

    /* An optional string may be passed to the signal handler which
     * will be used as the text for the exception. This can be set
     * using sig_str() instead of sig_on() or it can be changed by
     * set_sage_signal_handler_message() declared below.
     */
    const char* s;
};

/*
 * The actual object (there is a unique copy of this throughout Sage).
 */
extern struct sage_signals_t _signals;


/**********************************************************************
 * IMPLEMENTATION OF SIG_ON/SIG_OFF                                   *
 **********************************************************************/

/*
 * Implementation of sig_on().  Applications should not use this
 * directly, use sig_on() or sig_str() instead.
 *
 * _sig_on_(message) is a macro which pretends to be a function.
 * Since this is declared as "cdef except 0", Cython will know that an
 * exception occured if the value of _sig_on_() is 0 (false).
 *
 * INPUT:
 *
 *  - message -- a string to be displayed as error message when the code
 *    between sig_on() and sig_off() fails and raises an exception.
 *
 * OUTPUT: zero if an exception occured, non-zero otherwise.
 *
 * The function sigsetjmp() in the _sig_on_() macro can return:
 *  - zero: this happens in the actual sig_on() call. sigsetjmp() sets
 *    up the address for the Sage signal handler to jump to.  The
 *    program continues normally.
 *  - a signal number (e.g. 2 for SIGINT), assumed to be strictly
 *    positive: the Sage signal handler handled a signal.  Since
 *    _sig_on_() will return 0 in this case, the Exception (raised by
 *    sage_signal_handler) will be detected by Cython.
 *  - a negative number: this is assumed to come from sig_retry().  In
 *    this case, the program continues as if nothing happened between
 *    sig_on() and sig_retry().
 *
 * We cannot simply put sigsetjmp() in a function, because when that
 * function returns, we would lose the stack frame to siglongjmp() to.
 * That's why we need this hackish macro.  We use the fact that || is
 * a short-circuiting operator (the second argument is only evaluated
 * if the first returns 0).
 */
#define _sig_on_(message) ( unlikely(_sig_on_prejmp(message, __FILE__, __LINE__)) || _sig_on_postjmp(sigsetjmp(_signals.env,0)) )

/* This will be called during _sig_on_postjmp() when an interrupt was
 * received *before* the call to sig_on(). */
void _sig_on_interrupt_received(void);

/*
 * Set message, return 0 if we need to sigsetjmp(), return 1 otherwise.
 */
static inline int _sig_on_prejmp(const char* message, const char* file, int line)
{
    _signals.s = message;
#if ENABLE_DEBUG_INTERRUPT
    if (sage_interrupt_debug_level >= 4)
    {
        fprintf(stderr, "sig_on (count = %i) at %s:%i\n", _signals.sig_on_count+1, file, line);
        fflush(stderr);
    }
#endif
    if (_signals.sig_on_count > 0)
    {
        _signals.sig_on_count++;
        return 1;
    }

    /* At this point, _signals.sig_on_count == 0 */
    return 0;
}


/* Cleanup after siglongjmp() (reset signal mask to the default, set
 * sig_on_count to zero) */
void _sig_on_recover(void);

/*
 * Process the return value of sigsetjmp().
 * Return 0 if there was an exception, 1 otherwise.
 */
static inline int _sig_on_postjmp(int jmpret)
{
    if (unlikely(jmpret > 0))
    {
        /* An exception occured */
        _sig_on_recover();
        return 0;
    }

    /* When we are here, it's either the original sig_on() call or we
     * got here after sig_retry(). */
    _signals.sig_on_count = 1;

    /* Check whether we received an interrupt before this point.
     * _signals.interrupt_received can only be set by the interrupt
     * handler if _signals.sig_on_count is zero.  Because of that and
     * because _signals.sig_on_count and _signals.interrupt_received are
     * volatile, we can safely evaluate _signals.interrupt_received here
     * without race conditions. */
    if (unlikely(_signals.interrupt_received))
    {
        _sig_on_interrupt_received();
        return 0;
    }

    return 1;
}

/* Give a warning that sig_off() was called without sig_on() */
void _sig_off_warning(const char* file, int line);

/*
 * Implementation of sig_off().  Applications should not use this
 * directly, use sig_off() instead.
 */
static inline void _sig_off_(const char* file, int line)
{
#if ENABLE_DEBUG_INTERRUPT
    if (sage_interrupt_debug_level >= 4)
    {
        fprintf(stderr, "sig_off (count = %i) at %s:%i\n", _signals.sig_on_count, file, line);
        fflush(stderr);
    }
#endif
    if (unlikely(_signals.sig_on_count <= 0))
    {
        _sig_off_warning(file, line);
    }
    else
    {
        --_signals.sig_on_count;
    }
}


/**********************************************************************
 * USER MACROS/FUNCTIONS                                              *
 **********************************************************************/

/* The actual macros which should be used in a program. */
#define sig_on()           _sig_on_(NULL)
#define sig_str(message)   _sig_on_(message)
#define sig_off()          _sig_off_(__FILE__, __LINE__)

/* sig_check() should be functionally equivalent to sig_on(); sig_off();
 * but much faster.  Essentially, it checks whether we missed any
 * interrupts.
 *
 * OUTPUT: zero if an interrupt occured, non-zero otherwise.
 */
static inline int sig_check()
{
    if (unlikely(_signals.interrupt_received) && _signals.sig_on_count == 0)
    {
        _sig_on_interrupt_received();
        return 0;
    }

    return 1;
}

/* Macros behaving exactly like sig_on, sig_str and sig_check
 * but which are *not* declared cdef except 0.  This is useful if some
 * Cython code wants to do its own exception handling. */
#define sig_on_no_except()           sig_on()
#define sig_str_no_except(message)   sig_str(message)
#define sig_check_no_except()        sig_check()


/*
 * Temporarily block interrupts from happening inside sig_on().  This
 * is meant to wrap malloc() for example.  sig_unblock() checks whether
 * an interrupt happened in the mean time.  If yes, the interrupt is
 * re-raised.
 *
 * NOTES:
 * - This only works inside sig_on()/sig_off().  Outside of sig_on(),
 *   interrupts behave as usual.  This is because we can't propagate
 *   Python exceptions from low-level C code.
 * - Other signals still go through, because we can't really ignore
 *   SIGSEGV for example.
 * - For efficiency reasons, currently these may NOT be nested.
 *   Nesting could be implemented like src/headers/pariinl.h in PARI.
 */
static inline void sig_block()
{
#if ENABLE_DEBUG_INTERRUPT
    if (_signals.block_sigint != 0)
    {
        fprintf(stderr, "\n*** WARNING *** sig_block() with sig_on_count = %i, block_sigint = %i\n", _signals.sig_on_count, _signals.block_sigint);
        print_backtrace();
    }
#endif
    _signals.block_sigint = 1;
}

static inline void sig_unblock()
{
#if ENABLE_DEBUG_INTERRUPT
    if (_signals.block_sigint != 1)
    {
        fprintf(stderr, "\n*** WARNING *** sig_unblock() with sig_on_count = %i, block_sigint = %i\n", _signals.sig_on_count, _signals.block_sigint);
        print_backtrace();
    }
#endif
    _signals.block_sigint = 0;

    if (unlikely(_signals.interrupt_received) && _signals.sig_on_count > 0)
        kill(getpid(), _signals.interrupt_received);  /* Re-raise the signal */
}


/*
 * Call this before raising an exception to set _signals.s. The string
 * s will be used as the text for the exception.  Note that s is not
 * copied, we just store the pointer.
 */
void set_sage_signal_handler_message(const char* s);


/*
 * Retry a failed computation starting from sig_on().  This is useful
 * for PARI: if PARI complains that it doesn't have enough memory, we
 * allocate a larger stack and retry the computation.
 */
static inline void sig_retry()
{
    /* If we're outside of sig_on(), we can't jump, so we can only bail
     * out */
    if (unlikely(_signals.sig_on_count <= 0))
    {
        fprintf(stderr, "sig_retry() without sig_on()\n");
        abort();
    }
    siglongjmp(_signals.env, -1);
}

/* Used in error callbacks from C code (in particular NTL and PARI).
 * This should be used after an exception has been raised to jump back
 * to sig_on() where the exception will be seen. */
static inline void sig_error()
{
    if (unlikely(_signals.sig_on_count <= 0))
    {
        fprintf(stderr, "sig_error() without sig_on()\n");
    }
    abort();
}


/*
 * This function does nothing, but it is declared cdef except *, so it
 * can be used to make Cython check whether there is a pending exception
 * (PyErr_Occurred() is non-NULL).
 * To Cython, it will look like cython_check_exception() actually
 * raised the exception.
 */
static inline void cython_check_exception() {return;}


#ifdef __cplusplus
}  /* extern "C" */
#endif
#endif /* C_LIB_INCLUDE_INTERRUPT_H */
