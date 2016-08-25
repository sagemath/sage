"""
PARI C-library interface

AUTHORS:

- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta

- William Stein (2006-03-06): added newtonpoly

- Justin Walker: contributed some of the function definitions

- Gonzalo Tornaria: improvements to conversions; much better error
  handling.

- Robert Bradshaw, Jeroen Demeyer, William Stein (2010-08-15):
  Upgrade to PARI 2.4.3 (:trac:`9343`)

- Jeroen Demeyer (2011-11-12): rewrite various conversion routines
  (:trac:`11611`, :trac:`11854`, :trac:`11952`)

- Peter Bruin (2013-11-17): split off this file from gen.pyx
  (:trac:`15185`)

- Jeroen Demeyer (2014-02-09): upgrade to PARI 2.7 (:trac:`15767`)

- Jeroen Demeyer (2014-09-19): upgrade to PARI 2.8 (:trac:`16997`)

- Jeroen Demeyer (2015-03-17): automatically generate methods from
  ``pari.desc`` (:trac:`17631` and :trac:`17860`)

EXAMPLES::

    sage: pari('5! + 10/x')
    (120*x + 10)/x
    sage: pari('intnum(x=0,13,sin(x)+sin(x^2) + x)')
    83.8179442684285  # 32-bit
    84.1818153922297  # 64-bit
    sage: f = pari('x^3-1')
    sage: v = f.factor(); v
    [x - 1, 1; x^2 + x + 1, 1]
    sage: v[0]   # indexing is 0-based unlike in GP.
    [x - 1, x^2 + x + 1]~
    sage: v[1]
    [1, 1]~

Arithmetic obeys the usual coercion rules::

    sage: type(pari(1) + 1)
    <type 'sage.libs.pari.gen.gen'>
    sage: type(1 + pari(1))
    <type 'sage.libs.pari.gen.gen'>

GUIDE TO REAL PRECISION AND THE PARI LIBRARY

The default real precision in communicating with the PARI library
is the same as the default Sage real precision, which is 53 bits.
Inexact Pari objects are therefore printed by default to 15 decimal
digits (even if they are actually more precise).

Default precision example (53 bits, 15 significant decimals)::

    sage: a = pari(1.23); a
    1.23000000000000
    sage: a.sin()
    0.942488801931698

Example with custom precision of 200 bits (60 significant
decimals)::

    sage: R = RealField(200)
    sage: a = pari(R(1.23)); a   # only 15 significant digits printed
    1.23000000000000
    sage: R(a)         # but the number is known to precision of 200 bits
    1.2300000000000000000000000000000000000000000000000000000000
    sage: a.sin()      # only 15 significant digits printed
    0.942488801931698
    sage: R(a.sin())   # but the number is known to precision of 200 bits
    0.94248880193169751002382356538924454146128740562765030213504

It is possible to change the number of printed decimals::

    sage: R = RealField(200)    # 200 bits of precision in computations
    sage: old_prec = pari.set_real_precision(60)  # 60 decimals printed
    sage: a = pari(R(1.23)); a
    1.23000000000000000000000000000000000000000000000000000000000
    sage: a.sin()
    0.942488801931697510023823565389244541461287405627650302135038
    sage: pari.set_real_precision(old_prec)  # restore the default printing behavior
    60

Unless otherwise indicated in the docstring, most Pari functions
that return inexact objects use the precision of their arguments to
decide the precision of the computation. However, if some of these
arguments happen to be exact numbers (integers, rationals, etc.),
an optional parameter indicates the precision (in bits) to which
these arguments should be converted before the computation. If this
precision parameter is missing, the default precision of 53 bits is
used. The following first converts 2 into a real with 53-bit
precision::

    sage: R = RealField()
    sage: R(pari(2).sin())
    0.909297426825682

We can ask for a better precision using the optional parameter::

    sage: R = RealField(150)
    sage: R(pari(2).sin(precision=150))
    0.90929742682568169539601986591174484270225497

Warning regarding conversions Sage - Pari - Sage: Some care must be
taken when juggling inexact types back and forth between Sage and
Pari. In theory, calling p=pari(s) creates a Pari object p with the
same precision as s; in practice, the Pari library's precision is
word-based, so it will go up to the next word. For example, a
default 53-bit Sage real s will be bumped up to 64 bits by adding
bogus 11 bits. The function p.python() returns a Sage object with
exactly the same precision as the Pari object p. So
pari(s).python() is definitely not equal to s, since it has 64 bits
of precision, including the bogus 11 bits. The correct way of
avoiding this is to coerce pari(s).python() back into a domain with
the right precision. This has to be done by the user (or by Sage
functions that use Pari library functions in gen.pyx). For
instance, if we want to use the Pari library to compute sqrt(pi)
with a precision of 100 bits::

    sage: R = RealField(100)
    sage: s = R(pi); s
    3.1415926535897932384626433833
    sage: p = pari(s).sqrt()
    sage: x = p.python(); x  # wow, more digits than I expected!
    1.7724538509055160272981674833410973484
    sage: x.prec()           # has precision 'improved' from 100 to 128?
    128
    sage: x == RealField(128)(pi).sqrt()  # sadly, no!
    False
    sage: R(x)               # x should be brought back to precision 100
    1.7724538509055160272981674833
    sage: R(x) == s.sqrt()
    True

Elliptic curves and precision: If you are working with elliptic
curves, you should set the precision for each method::

    sage: e = pari([0,0,0,-82,0]).ellinit()
    sage: eta1 = e.elleta(precision=100)[0]
    sage: eta1.sage()
    3.6054636014326520859158205642077267748
    sage: eta1 = e.elleta(precision=180)[0]
    sage: eta1.sage()
    3.60546360143265208591582056420772677481026899659802474544

Number fields and precision: TODO

TESTS:

Check that output from PARI's print command is actually seen by
Sage (:trac:`9636`)::

    sage: pari('print("test")')
    test

Check that ``default()`` works properly::

    sage: pari.default("debug")
    0
    sage: pari.default("debug", 3)
    sage: pari(2^67+1).factor()
    IFAC: cracking composite
            49191317529892137643
    IFAC: factor 6713103182899
            is prime
    IFAC: factor 7327657
            is prime
    IFAC: prime 7327657
            appears with exponent = 1
    IFAC: prime 6713103182899
            appears with exponent = 1
    IFAC: found 2 large prime (power) factors.
    [3, 1; 7327657, 1; 6713103182899, 1]
    sage: pari.default("debug", 0)
    sage: pari(2^67+1).factor()
    [3, 1; 7327657, 1; 6713103182899, 1]
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from .paridecl cimport *
from .paripriv cimport *
include "cysignals/signals.pxi"
cdef extern from *:
    int sig_on_count "cysigs.sig_on_count"

import sys

cimport libc.stdlib
from libc.stdio cimport *
cimport cython

include "cysignals/memory.pxi"
from sage.ext.memory import init_memory_functions
from sage.structure.parent cimport Parent
from sage.libs.gmp.all cimport *
from sage.libs.flint.fmpz cimport fmpz_get_mpz, COEFF_IS_MPZ, COEFF_TO_PTR
from sage.libs.flint.fmpz_mat cimport *

from sage.libs.pari.gen cimport gen, objtogen
from sage.libs.pari.handle_error cimport _pari_init_error_handling
from sage.misc.superseded import deprecation, deprecated_function_alias

# real precision in decimal digits: see documentation for
# get_real_precision() and set_real_precision().  This variable is used
# in gp to set the precision of input quantities (e.g. sqrt(2)), and for
# determining the number of digits to be printed.  It is *not* used as
# a "default precision" for internal computations, which always use
# the actual precision of arguments together (where relevant) with a
# "prec" parameter.  In ALL cases (for real computations) the prec
# parameter is a WORD precision and NOT decimal precision.  Pari reals
# with word precision w have bit precision (of the mantissa) equal to
# 32*(w-2) or 64*(w-2).
#
# Hence the only relevance of this parameter in Sage is (1) for the
# output format of components of objects of type
# 'sage.libs.pari.gen.gen'; (2) for setting the precision of pari
# variables created from strings (e.g. via sage: pari('1.2')).
#
# WARNING: Many pari library functions take a last parameter "prec"
# which should be a words precision.  In many cases this is redundant
# and is simply ignored.  In our wrapping of these functions we use
# the variable prec here for convenience only.
cdef long prec

#################################################################
# conversions between various real precision models
#################################################################

def prec_bits_to_dec(unsigned long prec_in_bits):
    r"""
    Convert from precision expressed in bits to precision expressed in
    decimal.

    EXAMPLES::

        sage: from sage.libs.pari.pari_instance import prec_bits_to_dec
        sage: prec_bits_to_dec(53)
        15
        sage: [(32*n, prec_bits_to_dec(32*n)) for n in range(1, 9)]
        [(32, 9),
        (64, 19),
        (96, 28),
        (128, 38),
        (160, 48),
        (192, 57),
        (224, 67),
        (256, 77)]
    """
    cdef double log_2 = 0.301029995663981
    return int(prec_in_bits*log_2)

def prec_dec_to_bits(unsigned long prec_in_dec):
    r"""
    Convert from precision expressed in decimal to precision expressed
    in bits.

    EXAMPLES::

        sage: from sage.libs.pari.pari_instance import prec_dec_to_bits
        sage: prec_dec_to_bits(15)
        49
        sage: [(n, prec_dec_to_bits(n)) for n in range(10, 100, 10)]
        [(10, 33),
        (20, 66),
        (30, 99),
        (40, 132),
        (50, 166),
        (60, 199),
        (70, 232),
        (80, 265),
        (90, 298)]
    """
    cdef double log_10 = 3.32192809488736
    return int(prec_in_dec*log_10)

cpdef long prec_bits_to_words(unsigned long prec_in_bits):
    r"""
    Convert from precision expressed in bits to pari real precision
    expressed in words. Note: this rounds up to the nearest word,
    adjusts for the two codewords of a pari real, and is
    architecture-dependent.

    EXAMPLES::

        sage: from sage.libs.pari.pari_instance import prec_bits_to_words
        sage: prec_bits_to_words(70)
        5   # 32-bit
        4   # 64-bit

    ::

        sage: [(32*n, prec_bits_to_words(32*n)) for n in range(1, 9)]
        [(32, 3), (64, 4), (96, 5), (128, 6), (160, 7), (192, 8), (224, 9), (256, 10)] # 32-bit
        [(32, 3), (64, 3), (96, 4), (128, 4), (160, 5), (192, 5), (224, 6), (256, 6)] # 64-bit
    """
    if not prec_in_bits:
        return prec
    cdef unsigned long wordsize = BITS_IN_LONG

    # This equals ceil(prec_in_bits/wordsize) + 2
    return (prec_in_bits - 1)//wordsize + 3

cpdef long prec_words_to_bits(long prec_in_words):
    r"""
    Convert from pari real precision expressed in words to precision
    expressed in bits. Note: this adjusts for the two codewords of a
    pari real, and is architecture-dependent.

    EXAMPLES::

        sage: from sage.libs.pari.pari_instance import prec_words_to_bits
        sage: prec_words_to_bits(10)
        256   # 32-bit
        512   # 64-bit
        sage: [(n, prec_words_to_bits(n)) for n in range(3, 10)]
        [(3, 32), (4, 64), (5, 96), (6, 128), (7, 160), (8, 192), (9, 224)]  # 32-bit
        [(3, 64), (4, 128), (5, 192), (6, 256), (7, 320), (8, 384), (9, 448)] # 64-bit
    """
    # see user's guide to the pari library, page 10
    return (prec_in_words - 2) * BITS_IN_LONG

cpdef long default_bitprec():
    r"""
    Return the default precision in bits.

    EXAMPLES::

        sage: from sage.libs.pari.pari_instance import default_bitprec
        sage: default_bitprec()
        64
    """
    return (prec - 2) * BITS_IN_LONG

def prec_dec_to_words(long prec_in_dec):
    r"""
    Convert from precision expressed in decimal to precision expressed
    in words. Note: this rounds up to the nearest word, adjusts for the
    two codewords of a pari real, and is architecture-dependent.

    EXAMPLES::

        sage: from sage.libs.pari.pari_instance import prec_dec_to_words
        sage: prec_dec_to_words(38)
        6   # 32-bit
        4   # 64-bit
        sage: [(n, prec_dec_to_words(n)) for n in range(10, 90, 10)]
        [(10, 4), (20, 5), (30, 6), (40, 7), (50, 8), (60, 9), (70, 10), (80, 11)] # 32-bit
        [(10, 3), (20, 4), (30, 4), (40, 5), (50, 5), (60, 6), (70, 6), (80, 7)] # 64-bit
    """
    return prec_bits_to_words(prec_dec_to_bits(prec_in_dec))

def prec_words_to_dec(long prec_in_words):
    r"""
    Convert from precision expressed in words to precision expressed in
    decimal. Note: this adjusts for the two codewords of a pari real,
    and is architecture-dependent.

    EXAMPLES::

        sage: from sage.libs.pari.pari_instance import prec_words_to_dec
        sage: prec_words_to_dec(5)
        28   # 32-bit
        57   # 64-bit
        sage: [(n, prec_words_to_dec(n)) for n in range(3, 10)]
        [(3, 9), (4, 19), (5, 28), (6, 38), (7, 48), (8, 57), (9, 67)] # 32-bit
        [(3, 19), (4, 38), (5, 57), (6, 77), (7, 96), (8, 115), (9, 134)] # 64-bit
    """
    return prec_bits_to_dec(prec_words_to_bits(prec_in_words))


# The unique running Pari instance.
cdef PariInstance pari_instance, P
pari_instance = PariInstance()
P = pari_instance   # shorthand notation

# Also a copy of PARI accessible from external pure python code.
pari = pari_instance


# Callbacks from PARI to print stuff using sys.stdout.write() instead
# of C library functions like puts().
cdef PariOUT sage_pariOut

cdef void sage_putchar(char c):
    cdef char s[2]
    s[0] = c
    s[1] = 0
    sys.stdout.write(s)
    # Let PARI think the last character was a newline,
    # so it doesn't print one when an error occurs.
    pari_set_last_newline(1)

cdef void sage_puts(char* s):
    sys.stdout.write(s)
    pari_set_last_newline(1)

cdef void sage_flush():
    sys.stdout.flush()


include 'auto_instance.pxi'

@cython.final
cdef class PariInstance(PariInstance_auto):
    def __init__(self, long size=1000000, unsigned long maxprime=500000):
        """
        Initialize the PARI system.

        INPUT:


        -  ``size`` -- long, the number of bytes for the initial
           PARI stack (see note below)

        -  ``maxprime`` -- unsigned long, upper limit on a
           precomputed prime number table (default: 500000)


        .. note::

           In Sage, the PARI stack is different than in GP or the
           PARI C library. In Sage, instead of the PARI stack
           holding the results of all computations, it *only* holds
           the results of an individual computation. Each time a new
           Python/PARI object is computed, it it copied to its own
           space in the Python heap, and the memory it occupied on the
           PARI stack is freed. Thus it is not necessary to make the
           stack very large. Also, unlike in PARI, if the stack does
           overflow, in most cases the PARI stack is automatically
           increased and the relevant step of the computation rerun.

           This design obviously involves some performance penalties
           over the way PARI works, but it scales much better and is
           far more robust for large projects.

        .. note::

           If you do not want prime numbers, put ``maxprime=2``, but be
           careful because many PARI functions require this table. If
           you get the error message "not enough precomputed primes",
           increase this parameter.
        """
        if avma:
            return  # pari already initialized.

        # PARI has a "real" stack size (parisize) and a "virtual" stack
        # size (parisizemax). The idea is that the real stack will be
        # used if possible, but the stack might be increased up to
        # the complete virtual stack. Therefore, it is not a problem to
        # set the virtual stack size to a large value. There are two
        # constraints for the virtual stack size:
        # 1) on 32-bit systems, even virtual memory can be a scarce
        #    resource since it is limited by 4GB (of which the kernel
        #    needs a significant part)
        # 2) the system should actually be able to handle a stack size
        #    as large as the complete virtual stack.
        # As a simple heuristic, we set the virtual stack to 1/4 of the
        # virtual memory.

        from sage.misc.memory_info import MemoryInfo
        mem = MemoryInfo()

        pari_init_opts(size, maxprime, INIT_DFTm)
        paristack_setsize(size, mem.virtual_memory_limit() // 4)

        # Disable PARI's stack overflow checking which is incompatible
        # with multi-threading.
        pari_stackcheck_init(NULL)

        _pari_init_error_handling()

        # pari_init_opts() overrides MPIR's memory allocation functions,
        # so we need to reset them.
        init_memory_functions()

        # Set printing functions
        global pariOut, pariErr

        pariOut = &sage_pariOut
        pariOut.putch = sage_putchar
        pariOut.puts = sage_puts
        pariOut.flush = sage_flush

        # Display only 15 digits
        self._real_precision = 15
        sd_format("g.15", d_SILENT)

        # Init global prec variable (PARI's precision is always a
        # multiple of the machine word size)
        global prec
        prec = prec_bits_to_words(64)

        # Disable pretty-printing
        GP_DATA.fmt.prettyp = 0

        # This causes PARI/GP to use output independent of the terminal
        # (which is want we want for the PARI library interface).
        GP_DATA.flags = gpd_TEST

        # Ensure that Galois groups are represented in a sane way,
        # see the polgalois section of the PARI users manual.
        global new_galois_format
        new_galois_format = 1

        # By default, factor() should prove primality of returned
        # factors. This not only influences the factor() function, but
        # also many functions indirectly using factoring.
        global factor_proven
        factor_proven = 1

        # Initialize some constants
        sig_on()
        self.PARI_ZERO = self.new_gen_noclear(gen_0)
        self.PARI_ONE = self.new_gen_noclear(gen_1)
        self.PARI_TWO = self.new_gen_noclear(gen_2)
        sig_off()

    def debugstack(self):
        r"""
        Print the internal PARI variables ``top`` (top of stack), ``avma``
        (available memory address, think of this as the stack pointer),
        ``bot`` (bottom of stack).

        EXAMPLE::

            sage: pari.debugstack()  # random
            top =  0x60b2c60
            avma = 0x5875c38
            bot =  0x57295e0
            size = 1000000
        """
        # We deliberately use low-level functions to minimize the
        # chances that something goes wrong here (for example, if we
        # are out of memory).
        printf("top =  %p\navma = %p\nbot =  %p\nsize = %lu\n", pari_mainstack.top, avma, pari_mainstack.bot, <unsigned long>pari_mainstack.rsize)
        fflush(stdout)

    def __dealloc__(self):
        """
        Deallocation of the Pari instance.

        NOTE:

        Usually this deallocation happens only when Sage quits.
        We do not provide a direct test, since usually there
        is only one Pari instance, and when artificially creating
        another instance, C-data are shared.

        The fact that Sage does not crash when quitting is an
        indirect doctest. See the discussion at :trac:`13741`.

        """
        if avma:
            pari_close()

    def __repr__(self):
        return "Interface to the PARI C library"

    def __hash__(self):
        return 907629390   # hash('pari')

    cpdef _coerce_map_from_(self, x):
        """
        Return ``True`` if ``x`` admits a coercion map into the
        PARI interface.

        This currently always returns ``True``.

        EXAMPLES::

            sage: pari._coerce_map_from_(ZZ)
            True
            sage: pari.coerce_map_from(ZZ)
            Call morphism:
              From: Integer Ring
              To:   Interface to the PARI C library
        """
        return True

    def __richcmp__(left, right, int op):
        """
        EXAMPLES::

            sage: pari == pari
            True
            sage: pari == gp
            False
            sage: pari == 5
            False
        """
        return (<Parent>left)._richcmp(right, op)

    def set_debug_level(self, level):
        """
        Set the debug PARI C library variable.
        """
        self.default('debug', int(level))

    def get_debug_level(self):
        """
        Set the debug PARI C library variable.
        """
        return int(self.default('debug'))

    def set_real_precision(self, long n):
        """
        Sets the PARI default real precision in decimal digits.

        This is used both for creation of new objects from strings and for
        printing. It is the number of digits *IN DECIMAL* in which real
        numbers are printed. It also determines the precision of objects
        created by parsing strings (e.g. pari('1.2')), which is *not* the
        normal way of creating new pari objects in Sage. It has *no*
        effect on the precision of computations within the pari library.

        Returns the previous PARI real precision.

        EXAMPLES::

            sage: pari.set_real_precision(60)
            15
            sage: pari('1.2')
            1.20000000000000000000000000000000000000000000000000000000000
            sage: pari.set_real_precision(15)
            60
        """
        prev = self._real_precision
        cdef bytes strn = str(n)
        sig_on()
        sd_realprecision(strn, d_SILENT)
        sig_off()
        self._real_precision = n
        return prev

    def get_real_precision(self):
        """
        Returns the current PARI default real precision.

        This is used both for creation of new objects from strings and for
        printing. It is the number of digits *IN DECIMAL* in which real
        numbers are printed. It also determines the precision of objects
        created by parsing strings (e.g. pari('1.2')), which is *not* the
        normal way of creating new pari objects in Sage. It has *no*
        effect on the precision of computations within the pari library.

        EXAMPLES::

            sage: pari.get_real_precision()
            15
        """
        return self._real_precision

    def set_series_precision(self, long n):
        global precdl
        precdl = n

    def get_series_precision(self):
        return precdl

    cdef inline void clear_stack(self):
        """
        Call ``sig_off()``. If we are leaving the outermost
        ``sig_on() ... sig_off()`` block, then clear the PARI stack.
        """
        global avma
        if sig_on_count <= 1:
            avma = pari_mainstack.top
        sig_off()

    cdef inline gen new_gen(self, GEN x):
        """
        Create a new gen wrapping `x`, then call ``clear_stack()``.
        Except if `x` is ``gnil``, then we return ``None`` instead.
        """
        cdef gen g
        if x is gnil:
            g = None
        else:
            g = self.new_gen_noclear(x)
        self.clear_stack()
        return g

    cdef inline gen new_gen_noclear(self, GEN x):
        """
        Create a new gen, but don't free any memory on the stack and don't
        call sig_off().
        """
        cdef pari_sp address
        cdef gen y = gen.__new__(gen)
        y.g = self.deepcopy_to_python_heap(x, &address)
        y.b = address
        y._parent = self
        # y.refers_to (a dict which is None now) is initialised as needed
        return y

    cdef gen new_gen_from_mpz_t(self, mpz_t value):
        """
        Create a new gen from a given MPIR-integer ``value``.

        EXAMPLES::

            sage: pari(42)       # indirect doctest
            42

        TESTS:

        Check that the hash of an integer does not depend on existing
        garbage on the stack (:trac:`11611`)::

            sage: foo = pari(2^(32*1024));  # Create large integer to put PARI stack in known state
            sage: a5 = pari(5);
            sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1));  # Dirty PARI stack
            sage: b5 = pari(5);
            sage: a5.__hash__() == b5.__hash__()
            True
        """
        sig_on()
        return self.new_gen(self._new_GEN_from_mpz_t(value))

    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value):
        r"""
        Create a new PARI ``t_INT`` from a ``mpz_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        cdef unsigned long limbs = mpz_size(value)

        cdef GEN z = cgeti(limbs + 2)
        # Set sign and "effective length"
        z[1] = evalsigne(mpz_sgn(value)) + evallgefint(limbs + 2)
        mpz_export(int_LSW(z), NULL, -1, sizeof(long), 0, 0, value)

        return z

    cdef inline GEN _new_GEN_from_fmpz_t(self, fmpz_t value):
        r"""
        Create a new PARI ``t_INT`` from a ``fmpz_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        if COEFF_IS_MPZ(value[0]):
            return self._new_GEN_from_mpz_t(COEFF_TO_PTR(value[0]))
        else:
            return stoi(value[0])

    cdef gen new_gen_from_int(self, int value):
        sig_on()
        return self.new_gen(stoi(value))

    cdef gen new_gen_from_mpq_t(self, mpq_t value):
        """
        Create a new gen from a given MPIR-rational ``value``.

        EXAMPLES::

            sage: pari(-2/3)
            -2/3
            sage: pari(QQ(42))
            42
            sage: pari(QQ(42)).type()
            't_INT'
            sage: pari(QQ(1/42)).type()
            't_FRAC'

        TESTS:

        Check that the hash of a rational does not depend on existing
        garbage on the stack (:trac:`11854`)::

            sage: foo = pari(2^(32*1024));  # Create large integer to put PARI stack in known state
            sage: a5 = pari(5/7);
            sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1));  # Dirty PARI stack
            sage: b5 = pari(5/7);
            sage: a5.__hash__() == b5.__hash__()
            True
        """
        sig_on()
        return self.new_gen(self._new_GEN_from_mpq_t(value))

    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value):
        r"""
        Create a new PARI ``t_INT`` or ``t_FRAC`` from a ``mpq_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        cdef GEN num = self._new_GEN_from_mpz_t(mpq_numref(value))
        if mpz_cmpabs_ui(mpq_denref(value), 1) == 0:
            # Denominator is 1, return the numerator (an integer)
            return num
        cdef GEN denom = self._new_GEN_from_mpz_t(mpq_denref(value))
        return mkfrac(num, denom)

    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum):
        """
        Note that degree + 1 = length, so that recognizing 0 is easier.

        varnum = 0 is the general choice (creates a variable in x).
        """
        cdef GEN z
        cdef int i

        sig_on()
        z = cgetg(length + 2, t_POL)
        z[1] = evalvarn(varnum)
        if length != 0:
            setsigne(z,1)
            for i from 0 <= i < length:
                set_gel(z,i+2, stoi(vals[i]))
        else:
            ## polynomial is zero
            setsigne(z,0)

        return self.new_gen(z)

    cdef gen new_gen_from_padic(self, long ordp, long relprec,
                                mpz_t prime, mpz_t p_pow, mpz_t unit):
        cdef GEN z
        sig_on()
        z = cgetg(5, t_PADIC)
        z[1] = evalprecp(relprec) + evalvalp(ordp)
        set_gel(z, 2, self._new_GEN_from_mpz_t(prime))
        set_gel(z, 3, self._new_GEN_from_mpz_t(p_pow))
        set_gel(z, 4, self._new_GEN_from_mpz_t(unit))
        return self.new_gen(z)

    def double_to_gen(self, x):
        cdef double dx
        dx = float(x)
        return self.double_to_gen_c(dx)

    cdef gen double_to_gen_c(self, double x):
        """
        Create a new gen with the value of the double x, using Pari's
        dbltor.

        EXAMPLES::

            sage: pari.double_to_gen(1)
            1.00000000000000
            sage: pari.double_to_gen(1e30)
            1.00000000000000 E30
            sage: pari.double_to_gen(0)
            0.E-15
            sage: pari.double_to_gen(-sqrt(RDF(2)))
            -1.41421356237310
        """
        # Pari has an odd concept where it attempts to track the accuracy
        # of floating-point 0; a floating-point zero might be 0.0e-20
        # (meaning roughly that it might represent any number in the
        # range -1e-20 <= x <= 1e20).

        # Pari's dbltor converts a floating-point 0 into the Pari real
        # 0.0e-307; Pari treats this as an extremely precise 0.  This
        # can cause problems; for instance, the Pari incgam() function can
        # be very slow if the first argument is very precise.

        # So we translate 0 into a floating-point 0 with 53 bits
        # of precision (that's the number of mantissa bits in an IEEE
        # double).

        sig_on()
        if x == 0:
            return self.new_gen(real_0_bit(-53))
        else:
            return self.new_gen(dbltor(x))

    cdef GEN double_to_GEN(self, double x):
        if x == 0:
            return real_0_bit(-53)
        else:
            return dbltor(x)

    def complex(self, re, im):
        """
        Create a new complex number, initialized from re and im.
        """
        cdef gen t0 = self(re)
        cdef gen t1 = self(im)
        sig_on()
        return self.new_gen(mkcomplex(t0.g, t1.g))

    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address):
        cdef size_t s = <size_t> gsizebyte(x)
        cdef pari_sp tmp_bot = <pari_sp> sig_malloc(s)
        cdef pari_sp tmp_top = tmp_bot + s
        address[0] = tmp_bot
        return gcopy_avma(x, &tmp_top)

    cdef gen new_ref(self, GEN g, gen parent):
        """
        Create a new gen pointing to the given GEN, which is allocated as a
        part of parent.g.

        .. note::

           As a rule, there should never be more than one sage gen
           pointing to a given Pari GEN. So that means there is only
           one case where this function should be used: when a
           complicated Pari GEN is allocated with a single gen
           pointing to it, and one needs a gen pointing to one of its
           components.

           For example, doing x = pari("[1,2]") allocates a gen pointing to
           the list [1,2], but x[0] has no gen wrapping it, so new_ref
           should be used there. Then parent would be x in this
           case. See __getitem__ for an example of usage.

        EXAMPLES::

            sage: pari("[[1,2],3]")[0][1] ## indirect doctest
            2
        """
        cdef gen p = gen.__new__(gen)
        p.g = g
        p.b = 0
        p._parent = self
        p.refers_to = {-1: parent}
        return p

    def __call__(self, s):
        """
        Create the PARI object obtained by evaluating s using PARI.

        EXAMPLES::

            sage: pari([2,3,5])
            [2, 3, 5]
            sage: pari(Matrix(2,2,range(4)))
            [0, 1; 2, 3]
            sage: pari(x^2-3)
            x^2 - 3

        ::

            sage: a = pari(1); a, a.type()
            (1, 't_INT')
            sage: a = pari(1/2); a, a.type()
            (1/2, 't_FRAC')
            sage: a = pari(1/2); a, a.type()
            (1/2, 't_FRAC')

        See :func:`pari` for more examples.
        """
        return objtogen(s)

    cdef GEN _new_GEN_from_fmpz_mat_t(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc):
        r"""
        Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
        from a ``mpz_t**``.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        cdef GEN x
        cdef GEN A = zeromatcopy(nr, nc)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_fmpz_t(fmpz_mat_entry(B,i,j))
                set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
        return A

    cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc):
        r"""
        Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
        from a ``mpz_t**`` and rotate the matrix 90 degrees
        counterclockwise.  So the resulting matrix will have ``nc`` rows
        and ``nr`` columns.  This is useful for computing the Hermite
        Normal Form because Sage and PARI use different definitions.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        cdef GEN x
        cdef GEN A = zeromatcopy(nc, nr)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_fmpz_t(fmpz_mat_entry(B,i,nc-j-1))
                set_gcoeff(A, j+1, i+1, x)  # A[j+1, i+1] = x (using 1-based indexing)
        return A

    cdef gen integer_matrix(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf):
        """
        EXAMPLES::

            sage: matrix(ZZ,2,[1..6])._pari_()   # indirect doctest
            [1, 2, 3; 4, 5, 6]
        """
        sig_on()
        cdef GEN g
        if permute_for_hnf:
            g = self._new_GEN_from_fmpz_mat_t_rotate90(B, nr, nc)
        else:
            g = self._new_GEN_from_fmpz_mat_t(B, nr, nc)
        return self.new_gen(g)

    cdef GEN _new_GEN_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc):
        cdef GEN x
        # Allocate zero matrix
        cdef GEN A = zeromatcopy(nr, nc)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_mpq_t(B[i][j])
                set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
        return A

    cdef gen rational_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc):
        """
        EXAMPLES::

            sage: matrix(QQ,2,[1..6])._pari_()   # indirect doctest
            [1, 2, 3; 4, 5, 6]
        """
        sig_on()
        cdef GEN g = self._new_GEN_from_mpq_t_matrix(B, nr, nc)
        return self.new_gen(g)

    cdef _coerce_c_impl(self, x):
        """
        Implicit canonical coercion into a PARI object.
        """
        try:
            return self(x)
        except (TypeError, AttributeError):
            raise TypeError("no canonical coercion of %s into PARI" % x)

    cdef _an_element_c_impl(self):  # override this in Cython
        return self.PARI_ZERO

    cpdef gen zero(self):
        """
        EXAMPLES::

            sage: pari.zero()
            0
        """
        return self.PARI_ZERO

    cpdef gen one(self):
        """
        EXAMPLES::

            sage: pari.one()
            1
        """
        return self.PARI_ONE

    def new_with_bits_prec(self, s, long precision):
        r"""
        pari.new_with_bits_prec(self, s, precision) creates s as a PARI
        gen with (at most) precision *bits* of precision.
        """
        cdef unsigned long old_prec
        old_prec = GP_DATA.fmt.sigd
        precision = prec_bits_to_dec(precision)
        if not precision:
            precision = old_prec
        self.set_real_precision(precision)
        x = self(s)
        self.set_real_precision(old_prec)
        return x

    cdef long get_var(self, v) except -2:
        """
        Convert ``v`` into a PARI variable number.

        If ``v`` is a PARI object, return the variable number of
        ``variable(v)``. If ``v`` is ``None`` or ``-1``, return -1.
        Otherwise, treat ``v`` as a string and return the number of
        the variable named ``v``.

        OUTPUT: a PARI variable number (varn) or -1 if there is no
        variable number.

        .. WARNING::

            You can easily create variables with garbage names using
            this function. This can actually be used as a feature, if
            you want variable names which cannot be confused with
            ordinary user variables.

        EXAMPLES:

        We test this function using ``Pol()`` which calls this function::

            sage: pari("[1,0]").Pol()
            x
            sage: pari("[2,0]").Pol('x')
            2*x
            sage: pari("[Pi,0]").Pol('!@#$%^&')
            3.14159265358979*!@#$%^&

        We can use ``varhigher()`` and ``varlower()`` to create
        temporary variables without a name. The ``"xx"`` below is just a
        string to display the variable, it doesn't create a variable
        ``"xx"``::

            sage: xx = pari.varhigher("xx")
            sage: pari("[x,0]").Pol(xx)
            x*xx

        Indeed, this is not the same as::

            sage: pari("[x,0]").Pol("xx")
            Traceback (most recent call last):
            ...
            PariError: incorrect priority in gtopoly: variable x <= xx

        TESTS:

        The following example caused Sage to crash before
        :trac:`20630`::

            sage: R.<theta> = QQ[]
            sage: K.<a> = NumberField(theta^2 + 1)
            sage: K.galois_group(type='pari')
            Galois group PARI group [2, -1, 1, "S2"] of degree 2 of the Number Field in a with defining polynomial theta^2 + 1

        """
        if v is None:
            return -1
        cdef long varno
        if isinstance(v, gen):
            sig_on()
            varno = gvar((<gen>v).g)
            sig_off()
            if varno < 0:
                return -1
            else:
                return varno
        if v == -1:
            return -1
        cdef bytes s = bytes(v)
        sig_on()
        varno = fetch_user_var(s)
        sig_off()
        return varno

    ############################################################
    # Initialization
    ############################################################

    def stacksize(self):
        r"""
        Return the current size of the PARI stack, which is `10^6`
        by default.  However, the stack size is automatically doubled
        when needed up to some maximum.

        .. SEEALSO::

            - :meth:`stacksizemax` to get the maximum stack size
            - :meth:`allocatemem` to change the current or maximum
              stack size

        EXAMPLES::

            sage: pari.stacksize()
            1000000
            sage: pari.allocatemem(2^18, silent=True)
            sage: pari.stacksize()
            262144
        """
        return pari_mainstack.size

    def stacksizemax(self):
        r"""
        Return the maximum size of the PARI stack, which is determined
        at startup in terms of available memory. Usually, the PARI
        stack size is (much) smaller than this maximum but the stack
        will be increased up to this maximum if needed.

        .. SEEALSO::

            - :meth:`stacksize` to get the current stack size
            - :meth:`allocatemem` to change the current or maximum
              stack size

        EXAMPLES::

            sage: pari.allocatemem(2^18, 2^26, silent=True)
            sage: pari.stacksizemax()
            67108864
        """
        return pari_mainstack.vsize

    def allocatemem(self, size_t s=0, size_t sizemax=0, *, silent=False):
        r"""
        Change the PARI stack space to the given size ``s`` (or double
        the current size if ``s`` is `0`) and change the maximum stack
        size to ``sizemax``.

        PARI tries to use only its current stack (the size which is set
        by ``s``), but it will increase its stack if needed up to the
        maximum size which is set by ``sizemax``.

        The PARI stack is never automatically shrunk.  You can use the
        command ``pari.allocatemem(10^6)`` to reset the size to `10^6`,
        which is the default size at startup.  Note that the results of
        computations using Sage's PARI interface are copied to the
        Python heap, so they take up no space in the PARI stack.
        The PARI stack is cleared after every computation.

        It does no real harm to set this to a small value as the PARI
        stack will be automatically doubled when we run out of memory.

        INPUT:

        - ``s`` - an integer (default: 0).  A non-zero argument is the
          size in bytes of the new PARI stack.  If `s` is zero, double
          the current stack size.

        - ``sizemax`` - an integer (default: 0).  A non-zero argument
          is the maximum size in bytes of the PARI stack.  If
          ``sizemax`` is 0, the maximum of the current maximum and
          ``s`` is taken.

        EXAMPLES::

            sage: pari.allocatemem(10^7)
            PARI stack size set to 10000000 bytes, maximum size set to 67108864
            sage: pari.allocatemem()  # Double the current size
            PARI stack size set to 20000000 bytes, maximum size set to 67108864
            sage: pari.stacksize()
            20000000
            sage: pari.allocatemem(10^6)
            PARI stack size set to 1000000 bytes, maximum size set to 67108864

        The following computation will automatically increase the PARI
        stack size::

            sage: a = pari('2^100000000')

        ``a`` is now a Python variable on the Python heap and does not
        take up any space on the PARI stack.  The PARI stack is still
        large because of the computation of ``a``::

            sage: pari.stacksize()
            16000000

        Setting a small maximum size makes this fail::

            sage: pari.allocatemem(10^6, 2^22)
            PARI stack size set to 1000000 bytes, maximum size set to 4194304
            sage: a = pari('2^100000000')
            Traceback (most recent call last):
            ...
            PariError: _^s: the PARI stack overflows (current size: 1000000; maximum size: 4194304)
            You can use pari.allocatemem() to change the stack size and try again

        TESTS:

        Do the same without using the string interface and starting
        from a very small stack size::

            sage: pari.allocatemem(1, 2^26)
            PARI stack size set to 1024 bytes, maximum size set to 67108864
            sage: a = pari(2)^100000000
            sage: pari.stacksize()
            16777216

        We do not allow ``sizemax`` less than ``s``::

            sage: pari.allocatemem(10^7, 10^6)
            Traceback (most recent call last):
            ...
            ValueError: the maximum size (10000000) should be at least the stack size (1000000)
        """
        if s == 0:
            s = pari_mainstack.size * 2
            if s < pari_mainstack.size:
                raise OverflowError("cannot double stack size")
        elif s < 1024:
            s = 1024  # arbitrary minimum size
        if sizemax == 0:
            # For the default sizemax, use the maximum of current
            # sizemax and the given size s.
            if pari_mainstack.vsize > s:
                sizemax = pari_mainstack.vsize
            else:
                sizemax = s
        elif sizemax < s:
            raise ValueError("the maximum size ({}) should be at least the stack size ({})".format(s, sizemax))
        sig_on()
        paristack_setsize(s, sizemax)
        sig_off()
        if not silent:
            print("PARI stack size set to {} bytes, maximum size set to {}".
                format(self.stacksize(), self.stacksizemax()))

    def pari_version(self):
        return str(PARIVERSION)

    def init_primes(self, unsigned long M):
        """
        Recompute the primes table including at least all primes up to M
        (but possibly more).

        EXAMPLES::

            sage: pari.init_primes(200000)

        We make sure that ticket :trac:`11741` has been fixed::

            sage: pari.init_primes(2^30)
            Traceback (most recent call last):
            ...
            ValueError: Cannot compute primes beyond 436273290
        """
        # Hardcoded bound in PARI sources
        if M > 436273290:
            raise ValueError("Cannot compute primes beyond 436273290")

        if M <= maxprime():
            return
        sig_on()
        initprimetable(M)
        sig_off()

    def primes(self, n=None, end=None):
        """
        Return a pari vector containing the first `n` primes, the primes
        in the interval `[n, end]`, or the primes up to `end`.

        INPUT:

        Either

        - ``n`` -- integer

        or

        - ``n`` -- list or tuple `[a, b]` defining an interval of primes

        or

        - ``n, end`` -- start and end point of an interval of primes

        or

        - ``end`` -- end point for the list of primes

        OUTPUT: a PARI list of prime numbers

        EXAMPLES::

            sage: pari.primes(3)
            [2, 3, 5]
            sage: pari.primes(10)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
            sage: pari.primes(20)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
            sage: len(pari.primes(1000))
            1000
            sage: pari.primes(11,29)
            [11, 13, 17, 19, 23, 29]
            sage: pari.primes((11,29))
            [11, 13, 17, 19, 23, 29]
            sage: pari.primes(end=29)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
            sage: pari.primes(10^30, 10^30 + 100)
            [1000000000000000000000000000057, 1000000000000000000000000000099]

        TESTS::

            sage: pari.primes(0)
            []
            sage: pari.primes(-1)
            []
            sage: pari.primes(end=1)
            []
            sage: pari.primes(end=-1)
            []
            sage: pari.primes(3,2)
            []
        """
        cdef gen t0, t1
        if end is None:
            t0 = objtogen(n)
            sig_on()
            return self.new_gen(primes0(t0.g))
        elif n is None:
            t0 = self.PARI_TWO  # First prime
        else:
            t0 = objtogen(n)
        t1 = objtogen(end)
        sig_on()
        return self.new_gen(primes_interval(t0.g, t1.g))

    def primes_up_to_n(self, n):
        deprecation(20216, "pari.primes_up_to_n(n) is deprecated, use pari.primes(end=n) instead")
        return self.primes(end=n)

    prime_list = deprecated_function_alias(20216, primes)

    nth_prime = deprecated_function_alias(20216, PariInstance_auto.prime)

    euler = PariInstance_auto.Euler
    pi = PariInstance_auto.Pi

    def polchebyshev(self, long n, v=None):
        """
        Chebyshev polynomial of the first kind of degree `n`,
        in the variable `v`.

        EXAMPLES::

            sage: pari.polchebyshev(7)
            64*x^7 - 112*x^5 + 56*x^3 - 7*x
            sage: pari.polchebyshev(7, 'z')
            64*z^7 - 112*z^5 + 56*z^3 - 7*z
            sage: pari.polchebyshev(0)
            1
        """
        sig_on()
        return self.new_gen(polchebyshev1(n, self.get_var(v)))

    # Deprecated by upstream PARI: do not remove this deprecated alias
    # as long as it exists in PARI.
    poltchebi = deprecated_function_alias(18203, polchebyshev)

    def factorial(self, long n):
        """
        Return the factorial of the integer n as a PARI gen.

        EXAMPLES::

            sage: pari.factorial(0)
            1
            sage: pari.factorial(1)
            1
            sage: pari.factorial(5)
            120
            sage: pari.factorial(25)
            15511210043330985984000000
        """
        sig_on()
        return self.new_gen(mpfact(n))

    def polsubcyclo(self, long n, long d, v=None):
        """
        polsubcyclo(n, d, v=x): return the pari list of polynomial(s)
        defining the sub-abelian extensions of degree `d` of the
        cyclotomic field `\QQ(\zeta_n)`, where `d`
        divides `\phi(n)`.

        EXAMPLES::

            sage: pari.polsubcyclo(8, 4)
            [x^4 + 1]
            sage: pari.polsubcyclo(8, 2, 'z')
            [z^2 + 2, z^2 - 2, z^2 + 1]
            sage: pari.polsubcyclo(8, 1)
            [x - 1]
            sage: pari.polsubcyclo(8, 3)
            []
        """
        cdef gen plist
        sig_on()
        plist = self.new_gen(polsubcyclo(n, d, self.get_var(v)))
        if typ(plist.g) != t_VEC:
            return pari.vector(1, [plist])
        else:
            return plist

    polcyclo_eval = deprecated_function_alias(20217, PariInstance_auto.polcyclo)

    def setrand(self, seed):
        """
        Sets PARI's current random number seed.

        INPUT:

        - ``seed`` -- either a strictly positive integer or a GEN of
          type ``t_VECSMALL`` as output by ``getrand()``

        This should not be called directly; instead, use Sage's global
        random number seed handling in ``sage.misc.randstate``
        and call ``current_randstate().set_seed_pari()``.

        EXAMPLES::

            sage: pari.setrand(50)
            sage: a = pari.getrand()
            sage: pari.setrand(a)
            sage: a == pari.getrand()
            True

        TESTS:

        Check that invalid inputs are handled properly (:trac:`11825`)::

            sage: pari.setrand("foobar")
            Traceback (most recent call last):
            ...
            PariError: incorrect type in setrand (t_POL)
        """
        cdef gen t0 = self(seed)
        sig_on()
        setrand(t0.g)
        sig_off()

    def vector(self, long n, entries=None):
        """
        vector(long n, entries=None): Create and return the length n PARI
        vector with given list of entries.

        EXAMPLES::

            sage: pari.vector(5, [1, 2, 5, 4, 3])
            [1, 2, 5, 4, 3]
            sage: pari.vector(2, [x, 1])
            [x, 1]
            sage: pari.vector(2, [x, 1, 5])
            Traceback (most recent call last):
            ...
            IndexError: length of entries (=3) must equal n (=2)
        """
        cdef gen v = self._empty_vector(n)
        if entries is not None:
            if len(entries) != n:
                raise IndexError("length of entries (=%s) must equal n (=%s)"%\
                      (len(entries), n))
            for i, x in enumerate(entries):
                v[i] = x
        return v

    cdef gen _empty_vector(self, long n):
        cdef gen v
        sig_on()
        v = self.new_gen(zerovec(n))
        return v

    def matrix(self, long m, long n, entries=None):
        """
        matrix(long m, long n, entries=None): Create and return the m x n
        PARI matrix with given list of entries.
        """
        cdef long i, j, k
        cdef gen A
        cdef gen x

        sig_on()
        A = self.new_gen(zeromatcopy(m,n))
        if entries is not None:
            if len(entries) != m*n:
                raise IndexError("len of entries (=%s) must be %s*%s=%s"%(len(entries),m,n,m*n))
            k = 0
            A.refers_to = {}
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    x = pari(entries[k])
                    A.refers_to[(i,j)] = x
                    (<GEN>(A.g)[j+1])[i+1] = <long>(x.g)
                    k = k + 1
        return A

    def genus2red(self, P, P0=None):
        """
        Let `P` be a polynomial with integer coefficients.
        Determines the reduction of the (proper, smooth) genus 2
        curve `C/\QQ`, defined by the hyperelliptic equation `y^2 = P`.
        The special syntax ``genus2red([P,Q])`` is also allowed, where
        the polynomials `P` and `Q` have integer coefficients, to
        represent the model `y^2 + Q(x)y = P(x)`.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: pari.genus2red([-5*x^5, x^3 - 2*x^2 - 2*x + 1])
            [1416875, [2, -1; 5, 4; 2267, 1], x^6 - 240*x^4 - 2550*x^3 - 11400*x^2 - 24100*x - 19855, [[2, [2, [Mod(1, 2)]], []], [5, [1, []], ["[V] page 156", [3]]], [2267, [2, [Mod(432, 2267)]], ["[I{1-0-0}] page 170", []]]]]
        """
        cdef gen t0 = objtogen(P)
        sig_on()
        return self.new_gen(genus2red(t0.g, NULL))

    def List(self, x=None):
        """
        Create an empty list or convert `x` to a list.

        EXAMPLES::

            sage: pari.List(range(5))
            List([0, 1, 2, 3, 4])
            sage: L = pari.List()
            sage: L
            List([])
            sage: L.listput(42, 1)
            42
            sage: L
            List([42])
            sage: L.listinsert(24, 1)
            24
            sage: L
            List([24, 42])
        """
        if x is None:
            sig_on()
            return self.new_gen(listcreate())
        cdef gen t0 = objtogen(x)
        sig_on()
        return self.new_gen(gtolist(t0.g))


cdef inline void INT_to_mpz(mpz_ptr value, GEN g):
    """
    Store a PARI ``t_INT`` as an ``mpz_t``.
    """
    if typ(g) != t_INT:
        pari_err(e_TYPE, <char*>"conversion to mpz", g)

    cdef long size = lgefint(g) - 2
    mpz_import(value, size, -1, sizeof(long), 0, 0, int_LSW(g))

    if signe(g) < 0:
        mpz_neg(value, value)

cdef void INTFRAC_to_mpq(mpq_ptr value, GEN g):
    """
    Store a PARI ``t_INT`` or ``t_FRAC`` as an ``mpq_t``.
    """
    if typ(g) == t_FRAC:
        INT_to_mpz(mpq_numref(value), gel(g, 1))
        INT_to_mpz(mpq_denref(value), gel(g, 2))
    elif typ(g) == t_INT:
        INT_to_mpz(mpq_numref(value), g)
        mpz_set_ui(mpq_denref(value), 1)
    else:
        pari_err(e_TYPE, <char*>"conversion to mpq", g)
