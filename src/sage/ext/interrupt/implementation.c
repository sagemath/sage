/*
Interrupt and signal handling for Sage

AUTHORS:

- William Stein, Martin Albrecht (2006): initial version

- Jeroen Demeyer (2010-10-03): almost complete rewrite (#9678)

- Jeroen Demeyer (2013-01-11): handle SIGHUP also (#13908)

- Jeroen Demeyer (2013-01-28): handle SIGQUIT also (#14029)

- Jeroen Demeyer (2013-05-13): handle SIGALRM also (#13311)

- Jeroen Demeyer (2015-03-21): move to Cython (#18027)

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

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
/* glibc has a backtrace() command since version 2.1 */
#ifdef __GLIBC__
#if (__GLIBC__ > 2) || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 1)
#define HAVE_BACKTRACE 1
#include <execinfo.h>
#endif
#endif
#ifdef __linux__
#include <sys/prctl.h>
#endif
#include <pari/pari.h>
#include "interrupt/struct_signals.h"
#include "interrupt/interrupt.h"


/* Interrupt debug level.  This only works if ENABLE_DEBUG_INTERRUPT
 * has been set to "1" in debug.h */
#if ENABLE_DEBUG_INTERRUPT
static int default_debug_level = 2;
static struct timeval sigtime;  /* Time of signal */
#endif

/* The _signals object (there is a unique copy of this throughout Sage) */
static sage_signals_t _signals;

/* The default signal mask during normal operation,
 * initialized by setup_sage_signal_handler(). */
static sigset_t default_sigmask;

/* default_sigmask with SIGHUP, SIGINT, SIGALRM added. */
static sigset_t sigmask_with_sigint;


static void do_raise_exception(int sig);
static void sigdie(int sig, const char* s);
static void print_backtrace(void);


/* Do whatever is needed to reset the CPU to a sane state after
 * handling a signal.  In particular on x86 CPUs, we need to clear
 * the FPU (this is needed after MMX instructions have been used or
 * if an interrupt occurs during an FPU computation).
 * Linux and OS X 10.6 do this as part of their signals implementation,
 * but Solaris doesn't.  Since this code is called only when handling a
 * signal (which should be very rare), it's better to play safe and
 * always execute this instead of special-casing based on the operating
 * system.
 * See http://trac.sagemath.org/sage_trac/ticket/12873
 */
static inline void reset_CPU()
{
#if defined(__i386__) || defined(__x86_64__)
    /* Clear FPU tag word */
    asm("emms");
#endif
}


/* Handler for SIGHUP, SIGINT, SIGALRM
 *
 * Inside sig_on() (i.e. when _signals.sig_on_count is positive), this
 * raises an exception and jumps back to sig_on().
 * Outside of sig_on(), we set Python's interrupt flag using
 * PyErr_SetInterrupt() */
static void sage_interrupt_handler(int sig)
{
#if ENABLE_DEBUG_INTERRUPT
    if (_signals.debug_level >= 1) {
        fprintf(stderr, "\n*** SIG %i *** %s sig_on\n", sig, (_signals.sig_on_count > 0) ? "inside" : "outside");
        if (_signals.debug_level >= 3) print_backtrace();
        fflush(stderr);
        /* Store time of this signal, unless there is already a
         * pending signal. */
        if (!_signals.interrupt_received) gettimeofday(&sigtime, NULL);
    }
#endif

    if (_signals.sig_on_count > 0)
    {
        if (!_signals.block_sigint && !PARI_SIGINT_block)
        {
            /* Raise an exception so Python can see it */
            do_raise_exception(sig);

            /* Jump back to sig_on() (the first one if there is a stack) */
            reset_CPU();
            siglongjmp(_signals.env, sig);
        }
    }
    else
    {
        /* Set the Python interrupt indicator, which will cause the
         * Python-level interrupt handler in sage/ext/interrupt.pyx to
         * be called. */
        PyErr_SetInterrupt();
    }

    /* If we are here, we cannot handle the interrupt immediately, so
     * we store the signal number for later use.  But make sure we
     * don't overwrite a SIGHUP or SIGTERM which we already received. */
    if (_signals.interrupt_received != SIGHUP && _signals.interrupt_received != SIGTERM)
    {
        _signals.interrupt_received = sig;
        PARI_SIGINT_pending = sig;
    }
}

/* Handler for SIGQUIT, SIGILL, SIGABRT, SIGFPE, SIGBUS, SIGSEGV
 *
 * Inside sig_on() (i.e. when _signals.sig_on_count is positive), this
 * raises an exception and jumps back to sig_on().
 * Outside of sig_on(), we terminate Sage. */
static void sage_signal_handler(int sig)
{
    sig_atomic_t inside = _signals.inside_signal_handler;
    _signals.inside_signal_handler = 1;

    if (inside == 0 && _signals.sig_on_count > 0 && sig != SIGQUIT)
    {
        /* We are inside sig_on(), so we can handle the signal! */
#if ENABLE_DEBUG_INTERRUPT
        if (_signals.debug_level >= 1) {
            fprintf(stderr, "\n*** SIG %i *** inside sig_on\n", sig);
            if (_signals.debug_level >= 3) print_backtrace();
            fflush(stderr);
            gettimeofday(&sigtime, NULL);
        }
#endif

        /* Raise an exception so Python can see it */
        do_raise_exception(sig);

        /* Jump back to sig_on() (the first one if there is a stack) */
        reset_CPU();
        siglongjmp(_signals.env, sig);
    }
    else
    {
        /* We are outside sig_on() and have no choice but to terminate Sage */

        /* Reset all signals to their default behaviour and unblock
         * them in case something goes wrong as of now. */
        signal(SIGHUP, SIG_DFL);
        signal(SIGINT, SIG_DFL);
        signal(SIGQUIT, SIG_DFL);
        signal(SIGILL, SIG_DFL);
        signal(SIGABRT, SIG_DFL);
        signal(SIGFPE, SIG_DFL);
        signal(SIGBUS, SIG_DFL);
        signal(SIGSEGV, SIG_DFL);
        signal(SIGALRM, SIG_DFL);
        signal(SIGTERM, SIG_DFL);
        sigprocmask(SIG_SETMASK, &default_sigmask, NULL);

        if (inside) sigdie(sig, "An error occured during signal handling.");

        /* Quit Sage with an appropriate message. */
        switch(sig)
        {
            case SIGQUIT:
                sigdie(sig, NULL);
                break;  /* This will not be reached */
            case SIGILL:
                sigdie(sig, "Unhandled SIGILL: An illegal instruction occurred in Sage.");
                break;  /* This will not be reached */
            case SIGABRT:
                sigdie(sig, "Unhandled SIGABRT: An abort() occurred in Sage.");
                break;  /* This will not be reached */
            case SIGFPE:
                sigdie(sig, "Unhandled SIGFPE: An unhandled floating point exception occurred in Sage.");
                break;  /* This will not be reached */
            case SIGBUS:
                sigdie(sig, "Unhandled SIGBUS: A bus error occurred in Sage.");
                break;  /* This will not be reached */
            case SIGSEGV:
                sigdie(sig, "Unhandled SIGSEGV: A segmentation fault occurred in Sage.");
                break;  /* This will not be reached */
        };
        sigdie(sig, "Unknown signal received.\n");
    }
}


/* This calls sig_raise_exception() to actually raise the exception. */
static void do_raise_exception(int sig)
{
#if ENABLE_DEBUG_INTERRUPT
    struct timeval raisetime;
    if (_signals.debug_level >= 2) {
        gettimeofday(&raisetime, NULL);
        long delta_ms = (raisetime.tv_sec - sigtime.tv_sec)*1000L + ((long)raisetime.tv_usec - (long)sigtime.tv_usec)/1000;
        fprintf(stderr, "do_raise_exception(sig=%i)\nPyErr_Occurred() = %p\nRaising Python exception %li ms after signal...\n",
            sig, PyErr_Occurred(), delta_ms);
        fflush(stderr);
    }
#endif

    /* Call Cython function to raise exception */
    sig_raise_exception(sig, _signals.s);
}


/* This will be called during _sig_on_postjmp() when an interrupt was
 * received *before* the call to sig_on(). */
static void _sig_on_interrupt_received()
{
    /* Momentarily block signals to avoid race conditions */
    sigset_t oldset;
    sigprocmask(SIG_BLOCK, &sigmask_with_sigint, &oldset);

    do_raise_exception(_signals.interrupt_received);
    _signals.sig_on_count = 0;
    _signals.interrupt_received = 0;
    PARI_SIGINT_pending = 0;

    sigprocmask(SIG_SETMASK, &oldset, NULL);
}

/* Cleanup after siglongjmp() (reset signal mask to the default, set
 * sig_on_count to zero) */
static void _sig_on_recover()
{
    _signals.block_sigint = 0;
    PARI_SIGINT_block = 0;
    _signals.sig_on_count = 0;
    _signals.interrupt_received = 0;
    PARI_SIGINT_pending = 0;

    /* Reset signal mask */
    sigprocmask(SIG_SETMASK, &default_sigmask, NULL);
    _signals.inside_signal_handler = 0;
}

/* Give a warning that sig_off() was called without sig_on() */
static void _sig_off_warning(const char* file, int line)
{
    char buf[320];
    snprintf(buf, sizeof(buf), "sig_off() without sig_on() at %s:%i", file, line);

    /* Raise a warning with the Python GIL acquired */
    PyGILState_STATE gilstate_save = PyGILState_Ensure();
    PyErr_WarnEx(PyExc_RuntimeWarning, buf, 2);
    PyGILState_Release(gilstate_save);

    print_backtrace();
}


static void setup_sage_signal_handler()
{
    /* Reset the _signals structure */
    memset(&_signals, 0, sizeof(_signals));

    /* Save the default signal mask */
    sigprocmask(SIG_BLOCK, NULL, &default_sigmask);

    /* Save the signal mask with non-critical signals blocked */
    sigprocmask(SIG_BLOCK, NULL, &sigmask_with_sigint);
    sigaddset(&sigmask_with_sigint, SIGHUP);
    sigaddset(&sigmask_with_sigint, SIGINT);
    sigaddset(&sigmask_with_sigint, SIGALRM);

    /* Install signal handlers */
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    /* Block non-critical signals during the signal handlers */
    sigemptyset(&sa.sa_mask);
    sigaddset(&sa.sa_mask, SIGHUP);
    sigaddset(&sa.sa_mask, SIGINT);
    sigaddset(&sa.sa_mask, SIGALRM);

    sa.sa_handler = sage_interrupt_handler;
    if (sigaction(SIGHUP, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGINT, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGALRM, &sa, NULL)) {perror("sigaction"); exit(1);}
    sa.sa_handler = sage_signal_handler;
    /* Allow signals during signal handling, we have code to deal with
     * this case. */
    sa.sa_flags |= SA_NODEFER;
    if (sigaction(SIGQUIT, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGILL, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGABRT, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGFPE, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGBUS, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGSEGV, &sa, NULL)) {perror("sigaction"); exit(1);}

#if ENABLE_DEBUG_INTERRUPT
    _signals.debug_level = default_debug_level;
    if (_signals.debug_level >= 1)
        fprintf(stderr, "Finished setting up interrupts\n");
#endif
}


static void print_sep()
{
    fprintf(stderr,
        "------------------------------------------------------------------------\n");
    fflush(stderr);
}

/* Print a backtrace if supported by libc */
static void print_backtrace()
{
    void* backtracebuffer[1024];
    fflush(stderr);
#ifdef HAVE_BACKTRACE
    int btsize = backtrace(backtracebuffer, 1024);
    backtrace_symbols_fd(backtracebuffer, btsize, 2);
    print_sep();
#endif
}

/* Print a backtrace using gdb */
static void print_enhanced_backtrace()
{
    /* Bypass Linux Yama restrictions on ptrace() to allow debugging */
    /* See https://www.kernel.org/doc/Documentation/security/Yama.txt */
#ifdef PR_SET_PTRACER
    prctl(PR_SET_PTRACER, PR_SET_PTRACER_ANY, 0, 0, 0);
#endif

    /* Flush all buffers before forking */
    fflush(stdout);
    fflush(stderr);

    pid_t parent_pid = getpid();
    pid_t pid = fork();

    if (pid < 0)
    {
        /* Failed to fork: no problem, just ignore */
        perror("fork");
        return;
    }

    if (pid == 0) { /* child */
        /* We deliberately put these variables on the stack to avoid
         * malloc() calls, the heap might be messed up! */
        char path[1024];
        char pid_str[32];
        char* argv[5];

        snprintf(path, sizeof(path), "%s/bin/sage-CSI", getenv("SAGE_LOCAL"));
        snprintf(pid_str, sizeof(pid_str), "%i", parent_pid);

        argv[0] = "sage-CSI";
        argv[1] = "--no-color";
        argv[2] = "--pid";
        argv[3] = pid_str;
        argv[4] = NULL;
        execv(path, argv);
        perror("Failed to execute sage-CSI");
        exit(2);
    }
    /* Wait for sage-CSI to finish */
    waitpid(pid, NULL, 0);

    print_sep();
}


/* Print a message s and kill ourselves with signal sig */
static void sigdie(int sig, const char* s)
{
    print_sep();
    print_backtrace();

#if ENABLE_DEBUG_INTERRUPT
    /* Interrupt debugging is enabled, don't do enhanced backtraces as
     * the user is probably using other debugging tools and we don't
     * want to interfere with that. */
#else
#ifndef __APPLE__
    /* See http://trac.sagemath.org/13889 for how Apple screwed this up */
    print_enhanced_backtrace();
#endif
#endif

    if (s) {
        fprintf(stderr,
            "%s\n"
            "This probably occurred because a *compiled* component of Sage has a bug\n"
            "in it and is not properly wrapped with sig_on(), sig_off().\n"
            "Sage will now terminate.\n", s);
        print_sep();
    }

    /* Suicide with signal ``sig`` */
    kill(getpid(), sig);

    /* We should be dead! */
    exit(128 + sig);
}


/* Finally include the macros and inline functions for use in
 * interrupt.pyx. These require some of the above functions, therefore
 * this include must come at the end of this file. */
#include "interrupt/macros.h"
