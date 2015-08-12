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
and raise various exceptions (see sage/ext/interrupt.pyx).  SIGQUIT
will never be handled and always causes Sage to exit.


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
 *                     2010-2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/


#ifndef INTERRUPT_MACROS_H
#define INTERRUPT_MACROS_H

#include <setjmp.h>
#include <signal.h>
#include "interrupt/struct_signals.h"

#ifdef __cplusplus
extern "C" {
#endif


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

/*
 * Set message, return 0 if we need to sigsetjmp(), return 1 otherwise.
 */
static inline int _sig_on_prejmp(const char* message, const char* file, int line)
{
    _signals.s = message;
#if ENABLE_DEBUG_INTERRUPT
    if (_signals.debug_level >= 4)
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


/*
 * Implementation of sig_off().  Applications should not use this
 * directly, use sig_off() instead.
 */
static inline void _sig_off_(const char* file, int line)
{
#if ENABLE_DEBUG_INTERRUPT
    if (_signals.debug_level >= 4)
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
static inline int sig_check(void)
{
    if (unlikely(_signals.interrupt_received) && _signals.sig_on_count == 0)
    {
        _sig_on_interrupt_received();
        return 0;
    }

    return 1;
}


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
static inline void sig_block(void)
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

static inline void sig_unblock(void)
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
 * Retry a failed computation starting from sig_on().  This is useful
 * for PARI: if PARI complains that it doesn't have enough memory, we
 * allocate a larger stack and retry the computation.
 */
static inline void sig_retry(void)
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
static inline void sig_error(void)
{
    if (unlikely(_signals.sig_on_count <= 0))
    {
        fprintf(stderr, "sig_error() without sig_on()\n");
    }
    abort();
}


#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif  /* ifndef INTERRUPT_MACROS_H */
