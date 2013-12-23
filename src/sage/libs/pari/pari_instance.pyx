"""
PARI C-library interface

AUTHORS:

- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta

- William Stein (2006-03-06): added newtonpoly

- Justin Walker: contributed some of the function definitions

- Gonzalo Tornaria: improvements to conversions; much better error
  handling.

- Robert Bradshaw, Jeroen Demeyer, William Stein (2010-08-15):
  Upgrade to PARI 2.4.3 (#9343)

- Jeroen Demeyer (2011-11-12): rewrite various conversion routines
  (#11611, #11854, #11952)

- Peter Bruin (2013-11-17): split off this file from gen.pyx (#15185)


EXAMPLES::

    sage: pari('5! + 10/x')
    (120*x + 10)/x
    sage: pari('intnum(x=0,13,sin(x)+sin(x^2) + x)')
    85.1885681951527
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

The default real precision in communicating with the Pari library
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
curves and want to compute with a precision other than the default
53 bits, you should use the precision parameter of ellinit()::

    sage: R = RealField(150)
    sage: e = pari([0,0,0,-82,0]).ellinit(precision=150)
    sage: eta1 = e.elleta()[0]
    sage: R(eta1)
    3.6054636014326520859158205642077267748102690

Number fields and precision: TODO

TESTS:

Check that output from PARI's print command is actually seen by
Sage (ticket #9636)::

    sage: pari('print("test")')
    test

"""

include 'pari_err.pxi'
include 'sage/ext/stdsage.pxi'
include 'sage/ext/python.pxi'
include 'sage/ext/interrupt.pxi'

import sys

cimport libc.stdlib
cimport cython

from sage.structure.parent cimport Parent

from sage.libs.pari.gen cimport gen, objtogen
from sage.libs.pari.handle_error cimport pari_error_string, \
        _pari_init_error_handling, _pari_check_warning, \
        _pari_handle_exception, _pari_err_recover

# so Galois groups are represented in a sane way
# See the polgalois section of the PARI users manual.
new_galois_format = 1

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

def prec_words_to_bits(long prec_in_words):
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

# PariInstance.__init__ must not create gen objects because their parent is not constructed yet
pari_catch_sig_on()
pari_instance.PARI_ZERO = pari_instance.new_gen_noclear(gen_0)
pari_instance.PARI_ONE  = pari_instance.new_gen_noclear(gen_1)
pari_instance.PARI_TWO  = pari_instance.new_gen_noclear(gen_2)
pari_catch_sig_off()

# Also a copy of PARI accessible from external pure python code.
pari = pari_instance


cdef unsigned long num_primes

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

cdef PariOUT sage_pariErr

cdef void sage_pariErr_putchar(char c):
    cdef char s[2]
    s[0] = c
    s[1] = 0
    global pari_error_string
    pari_error_string += str(s)
    pari_set_last_newline(1)

cdef void sage_pariErr_puts(char *s):
    global pari_error_string
    pari_error_string += str(s)
    pari_set_last_newline(1)

cdef void sage_pariErr_flush():
    pass


@cython.final
cdef class PariInstance(sage.structure.parent_base.ParentWithBase):
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
        if bot:
            return  # pari already initialized.

        global num_primes, top, bot

        # The size here doesn't really matter, because we will allocate
        # our own stack anyway. We ask PARI not to set up signal and
        # error handlers.
        pari_init_opts(10000, maxprime, INIT_DFTm)

        _pari_init_error_handling()

        num_primes = maxprime

        # Free the PARI stack and allocate our own (using Cython)
        pari_free(<void*>bot); bot = 0
        init_stack(size)

        # Set printing functions
        global pariOut, pariErr

        pariOut = &sage_pariOut
        pariOut.putch = sage_putchar
        pariOut.puts = sage_puts
        pariOut.flush = sage_flush

        pariErr = &sage_pariErr
        pariErr.putch = sage_pariErr_putchar
        pariErr.puts = sage_pariErr_puts
        pariErr.flush = sage_pariErr_flush

        # Display only 15 digits
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
        global avma, top, bot
        printf("top =  %p\navma = %p\nbot =  %p\nsize = %lu\n", top, avma, bot, <unsigned long>(top) - <unsigned long>(bot))
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
        if bot:
            sage_free(<void*>bot)
        global top, bot, avma
        top = 0
        bot = 0
        avma = 0
        pari_close()

    def __repr__(self):
        return "Interface to the PARI C library"

    def __hash__(self):
        return 907629390   # hash('pari')

    cdef has_coerce_map_from_c_impl(self, x):
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

    def default(self, variable, value=None):
        if not value is None:
            return self('default(%s, %s)'%(variable, value))
        return self('default(%s)'%variable)

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
        prev = self.get_real_precision()
        cdef bytes strn = str(n)
        pari_catch_sig_on()
        sd_realprecision(strn, d_SILENT)
        pari_catch_sig_off()
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
        pari_catch_sig_on()
        cdef long prev = itos(sd_realprecision(NULL, d_RETURN))
        pari_catch_sig_off()
        return prev

    def set_series_precision(self, long n):
        global precdl
        precdl = n

    def get_series_precision(self):
        return precdl

    cdef inline void clear_stack(self):
        """
        Call ``pari_catch_sig_off()``, and clear the entire PARI stack
        if we are leaving the outermost ``pari_catch_sig_on() ...
        pari_catch_sig_off()`` block.

        """
        global top, avma
        if _signals.sig_on_count <= 1:
            avma = top
        pari_catch_sig_off()

    cdef inline gen new_gen(self, GEN x):
        """
        Create a new gen wrapping `x`, then call ``clear_stack()``.
        """
        cdef gen g = _new_gen(x)
        self.clear_stack()
        return g

    cdef inline gen new_gen_noclear(self, GEN x):
        """
        Create a new gen, but don't free any memory on the stack and don't
        call pari_catch_sig_off().
        """
        z = _new_gen(x)
        return z

    cdef gen new_gen_from_mpz_t(self, mpz_t value):
        """
        Create a new gen from a given MPIR-integer ``value``.

        EXAMPLES::

            sage: pari(42)       # indirect doctest
            42

        TESTS:

        Check that the hash of an integer does not depend on existing
        garbage on the stack (#11611)::

            sage: foo = pari(2^(32*1024));  # Create large integer to put PARI stack in known state
            sage: a5 = pari(5);
            sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1));  # Dirty PARI stack
            sage: b5 = pari(5);
            sage: a5.__hash__() == b5.__hash__()
            True
        """
        pari_catch_sig_on()
        return self.new_gen(self._new_GEN_from_mpz_t(value))

    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value):
        r"""
        Create a new PARI ``t_INT`` from a ``mpz_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``pari_catch_sig_on()`` before and ``pari_catch_sig_off()`` after.
        """
        cdef unsigned long limbs = mpz_size(value)

        cdef GEN z = cgeti(limbs + 2)
        # Set sign and "effective length"
        z[1] = evalsigne(mpz_sgn(value)) + evallgefint(limbs + 2)
        mpz_export(int_LSW(z), NULL, -1, sizeof(long), 0, 0, value)

        return z

    cdef gen new_gen_from_int(self, int value):
        pari_catch_sig_on()
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
        garbage on the stack (#11854)::

            sage: foo = pari(2^(32*1024));  # Create large integer to put PARI stack in known state
            sage: a5 = pari(5/7);
            sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1));  # Dirty PARI stack
            sage: b5 = pari(5/7);
            sage: a5.__hash__() == b5.__hash__()
            True
        """
        pari_catch_sig_on()
        return self.new_gen(self._new_GEN_from_mpq_t(value))

    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value):
        r"""
        Create a new PARI ``t_INT`` or ``t_FRAC`` from a ``mpq_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``pari_catch_sig_on()`` before and ``pari_catch_sig_off()`` after.
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

        pari_catch_sig_on()
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
        pari_catch_sig_on()
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

        pari_catch_sig_on()
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
        pari_catch_sig_on()
        return self.new_gen(mkcomplex(t0.g, t1.g))

    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address):
        return deepcopy_to_python_heap(x, address)

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
        cdef gen p = PY_NEW(gen)

        p.b = 0
        p._parent = self
        p._refers_to = {-1:parent}
        p.g = g

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

    cdef GEN _new_GEN_from_mpz_t_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc):
        r"""
        Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
        from a ``mpz_t**``.

        For internal use only; this directly uses the PARI stack.
        One should call ``pari_catch_sig_on()`` before and ``pari_catch_sig_off()`` after.
        """
        cdef GEN x
        cdef GEN A = zeromatcopy(nr, nc)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_mpz_t(B[i][j])
                set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
        return A

    cdef GEN _new_GEN_from_mpz_t_matrix_rotate90(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc):
        r"""
        Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
        from a ``mpz_t**`` and rotate the matrix 90 degrees
        counterclockwise.  So the resulting matrix will have ``nc`` rows
        and ``nr`` columns.  This is useful for computing the Hermite
        Normal Form because Sage and PARI use different definitions.

        For internal use only; this directly uses the PARI stack.
        One should call ``pari_catch_sig_on()`` before and ``pari_catch_sig_off()`` after.
        """
        cdef GEN x
        cdef GEN A = zeromatcopy(nc, nr)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_mpz_t(B[i][nc-j-1])
                set_gcoeff(A, j+1, i+1, x)  # A[j+1, i+1] = x (using 1-based indexing)
        return A

    cdef gen integer_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf):
        """
        EXAMPLES::

            sage: matrix(ZZ,2,[1..6])._pari_()   # indirect doctest
            [1, 2, 3; 4, 5, 6]
        """
        pari_catch_sig_on()
        cdef GEN g
        if permute_for_hnf:
            g = self._new_GEN_from_mpz_t_matrix_rotate90(B, nr, nc)
        else:
            g = self._new_GEN_from_mpz_t_matrix(B, nr, nc)
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
        pari_catch_sig_on()
        cdef GEN g = self._new_GEN_from_mpq_t_matrix(B, nr, nc)
        return self.new_gen(g)

    cdef _coerce_c_impl(self, x):
        """
        Implicit canonical coercion into a PARI object.
        """
        try:
            return self(x)
        except (TypeError, AttributeError):
            raise TypeError("no canonical coercion of %s into PARI"%x)

    cdef _an_element_c_impl(self):  # override this in Cython
        return self.PARI_ZERO

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



    cdef long get_var(self, v):
        """
        Converts a Python string into a PARI variable reference number. Or
        if v = -1, returns -1.
        """
        if v != -1:
            s = str(v)
            return fetch_user_var(s)
        return -1

    ############################################################
    # Initialization
    ############################################################

    def stacksize(self):
        r"""
        Returns the current size of the PARI stack, which is `10^6`
        by default.  However, the stack size is automatically doubled
        when needed.  It can also be set directly using
        ``pari.allocatemem()``.

        EXAMPLES::

            sage: pari.stacksize()
            1000000

        """
        global top, bot
        return (<size_t>(top) - <size_t>(bot))

    def _allocate_huge_mem(self):
        r"""
        This function tries to allocate an impossibly large amount of
        PARI stack in order to test ``init_stack()`` failing.

        TESTS::

            sage: pari.allocatemem(10^6, silent=True)
            sage: pari._allocate_huge_mem()
            Traceback (most recent call last):
            ...
            MemoryError: Unable to allocate ... (instead, allocated 1000000 bytes)

        Test that the PARI stack is sane after this failure::

            sage: a = pari('2^10000000')
            sage: pari.allocatemem(10^6, silent=True)
        """
        # Since size_t is unsigned, this should wrap over to a very
        # large positive number.
        init_stack(<size_t>(-4096))

    def _setup_failed_retry(self):
        r"""
        This function pretends that the PARI stack is larger than it
        actually is such that allocatemem(0) will allocate much more
        than double the current stack.  That allocation will then fail.
        This function is meant to be used only in this doctest.

        TESTS::

            sage: pari.allocatemem(4096, silent=True)
            sage: pari._setup_failed_retry()
            sage: try:  # Any result is fine as long as it's not a segfault
            ....:     a = pari('2^1000000')
            ....: except MemoryError:
            ....:     pass
            sage: pari.allocatemem(10^6, silent=True)
        """
        global top
        # Pretend top is at a very high address
        top = <pari_sp>(<size_t>(-16))

    def allocatemem(self, s=0, silent=False):
        r"""
        Change the PARI stack space to the given size (or double the
        current size if ``s`` is `0`).

        If `s = 0` and insufficient memory is avaible to double, the
        PARI stack will be enlarged by a smaller amount.  In any case,
        a ``MemoryError`` will be raised if the requested memory cannot
        be allocated.

        The PARI stack is never automatically shrunk.  You can use the
        command ``pari.allocatemem(10^6)`` to reset the size to `10^6`,
        which is the default size at startup.  Note that the results of
        computations using Sage's PARI interface are copied to the
        Python heap, so they take up no space in the PARI stack.
        The PARI stack is cleared after every computation.

        It does no real harm to set this to a small value as the PARI
        stack will be automatically doubled when we run out of memory.
        However, it could make some calculations slower (since they have
        to be redone from the start after doubling the stack until the
        stack is big enough).

        INPUT:

        - ``s`` - an integer (default: 0).  A non-zero argument should
          be the size in bytes of the new PARI stack.  If `s` is zero,
          try to double the current stack size.

        EXAMPLES::

            sage: pari.allocatemem(10^7)
            PARI stack size set to 10000000 bytes
            sage: pari.allocatemem()  # Double the current size
            PARI stack size set to 20000000 bytes
            sage: pari.stacksize()
            20000000
            sage: pari.allocatemem(10^6)
            PARI stack size set to 1000000 bytes

        The following computation will automatically increase the PARI
        stack size::

            sage: a = pari('2^100000000')

        ``a`` is now a Python variable on the Python heap and does not
        take up any space on the PARI stack.  The PARI stack is still
        large because of the computation of ``a``::

            sage: pari.stacksize()
            16000000
            sage: pari.allocatemem(10^6)
            PARI stack size set to 1000000 bytes
            sage: pari.stacksize()
            1000000

        TESTS:

        Do the same without using the string interface and starting
        from a very small stack size::

            sage: pari.allocatemem(1)
            PARI stack size set to 1024 bytes
            sage: a = pari(2)^100000000
            sage: pari.stacksize()
            16777216
        """
        s = long(s)
        if s < 0:
            raise ValueError("Stack size must be nonnegative")
        init_stack(s)
        if not silent:
            print "PARI stack size set to", self.stacksize(), "bytes"

    def pari_version(self):
        return str(PARIVERSION)

    def init_primes(self, _M):
        """
        Recompute the primes table including at least all primes up to M
        (but possibly more).

        EXAMPLES::

            sage: pari.init_primes(200000)

        We make sure that ticket #11741 has been fixed, and double check to
        make sure that diffptr has not been freed::

            sage: pari.init_primes(2^62)
            Traceback (most recent call last):
            ...
            PariError: not enough memory                  # 64-bit
            OverflowError: long int too large to convert  # 32-bit
            sage: pari.init_primes(200000)
        """
        cdef unsigned long M
        cdef char *tmpptr
        M = _M
        global diffptr, num_primes
        if M <= num_primes:
            return
        pari_catch_sig_on()
        tmpptr = initprimes(M)
        pari_catch_sig_off()
        pari_free(<void*> diffptr)
        num_primes = M
        diffptr = tmpptr


    ##############################################
    ## Support for GP Scripts
    ##############################################

    def read(self, bytes filename):
        r"""
        Read a script from the named filename into the interpreter.  The
        functions defined in the script are then available for use from
        Sage/PARI.  The result of the last command in ``filename`` is
        returned.

        EXAMPLES:

        Create a gp file::

            sage: import tempfile
            sage: gpfile = tempfile.NamedTemporaryFile(mode="w")
            sage: gpfile.file.write("mysquare(n) = {\n")
            sage: gpfile.file.write("    n^2;\n")
            sage: gpfile.file.write("}\n")
            sage: gpfile.file.write("polcyclo(5)\n")
            sage: gpfile.file.flush()

        Read it in Sage, we get the result of the last line::

            sage: pari.read(gpfile.name)
            x^4 + x^3 + x^2 + x + 1

        Call the function defined in the gp file::

            sage: pari('mysquare(12)')
            144
        """
        pari_catch_sig_on()
        return self.new_gen(gp_read_file(filename))


    ##############################################

    def _primelimit(self):
        """
        Return the number of primes already computed
        in this Pari instance.

        EXAMPLES:
            sage: pari._primelimit()
            500000
            sage: pari.init_primes(600000)
            sage: pari._primelimit()
            600000
        """
        global num_primes
        from sage.rings.all import ZZ
        return ZZ(num_primes)

    def prime_list(self, long n):
        """
        prime_list(n): returns list of the first n primes

        To extend the table of primes use pari.init_primes(M).

        INPUT:


        -  ``n`` - C long


        OUTPUT:


        -  ``gen`` - PARI list of first n primes


        EXAMPLES::

            sage: pari.prime_list(0)
            []
            sage: pari.prime_list(-1)
            []
            sage: pari.prime_list(3)
            [2, 3, 5]
            sage: pari.prime_list(10)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
            sage: pari.prime_list(20)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
            sage: len(pari.prime_list(1000))
            1000
        """
        if n >= 2:
            self.nth_prime(n)
        pari_catch_sig_on()
        return self.new_gen(primes(n))

    def primes_up_to_n(self, long n):
        """
        Return the primes <= n as a pari list.

        EXAMPLES::

            sage: pari.primes_up_to_n(1)
            []
            sage: pari.primes_up_to_n(20)
            [2, 3, 5, 7, 11, 13, 17, 19]
        """
        if n <= 1:
            return pari([])
        self.init_primes(n+1)
        return self.prime_list(pari(n).primepi())

##         cdef long k
##         k = (n+10)/math.log(n)
##         p = 2
##         while p <= n:
##             p = self.nth_prime(k)
##             k = 2
##         v = self.prime_list(k)
##         return v[:pari(n).primepi()]

    def __nth_prime(self, long n):
        """
        nth_prime(n): returns the n-th prime, where n is a C-int
        """
        global num_primes

        if n <= 0:
            raise ValueError, "nth prime meaningless for non-positive n (=%s)"%n
        cdef GEN g
        pari_catch_sig_on()
        g = prime(n)
        return self.new_gen(g)


    def nth_prime(self, long n):
        from sage.libs.pari.all import PariError
        try:
            return self.__nth_prime(n)
        except PariError:
            self.init_primes(max(2*num_primes,20*n))
            return self.nth_prime(n)

    def euler(self, unsigned long precision=0):
        """
        Return Euler's constant to the requested real precision (in bits).

        EXAMPLES::

            sage: pari.euler()
            0.577215664901533
            sage: pari.euler(precision=100).python()
            0.577215664901532860606512090082...
        """
        pari_catch_sig_on()
        return self.new_gen(mpeuler(prec_bits_to_words(precision)))

    def pi(self, unsigned long precision=0):
        """
        Return the value of the constant pi to the requested real precision
        (in bits).

        EXAMPLES::

            sage: pari.pi()
            3.14159265358979
            sage: pari.pi(precision=100).python()
            3.1415926535897932384626433832...
        """
        pari_catch_sig_on()
        return self.new_gen(mppi(prec_bits_to_words(precision)))

    def pollegendre(self, long n, v=-1):
        """
        pollegendre(n, v=x): Legendre polynomial of degree n (n C-integer),
        in variable v.

        EXAMPLES::

            sage: pari.pollegendre(7)
            429/16*x^7 - 693/16*x^5 + 315/16*x^3 - 35/16*x
            sage: pari.pollegendre(7, 'z')
            429/16*z^7 - 693/16*z^5 + 315/16*z^3 - 35/16*z
            sage: pari.pollegendre(0)
            1
        """
        pari_catch_sig_on()
        return self.new_gen(pollegendre(n, self.get_var(v)))

    def poltchebi(self, long n, v=-1):
        """
        poltchebi(n, v=x): Chebyshev polynomial of the first kind of degree
        n, in variable v.

        EXAMPLES::

            sage: pari.poltchebi(7)
            64*x^7 - 112*x^5 + 56*x^3 - 7*x
            sage: pari.poltchebi(7, 'z')
            64*z^7 - 112*z^5 + 56*z^3 - 7*z
            sage: pari.poltchebi(0)
            1
        """
        pari_catch_sig_on()
        return self.new_gen(polchebyshev1(n, self.get_var(v)))

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
        pari_catch_sig_on()
        return self.new_gen(mpfact(n))

    def polcyclo(self, long n, v=-1):
        """
        polcyclo(n, v=x): cyclotomic polynomial of degree n, in variable
        v.

        EXAMPLES::

            sage: pari.polcyclo(8)
            x^4 + 1
            sage: pari.polcyclo(7, 'z')
            z^6 + z^5 + z^4 + z^3 + z^2 + z + 1
            sage: pari.polcyclo(1)
            x - 1
        """
        pari_catch_sig_on()
        return self.new_gen(polcyclo(n, self.get_var(v)))

    def polcyclo_eval(self, long n, v):
        """
        polcyclo_eval(n, v): value of the nth cyclotomic polynomial at value v.

        EXAMPLES::

            sage: pari.polcyclo_eval(8, 2)
            17
            sage: cyclotomic_polynomial(8)(2)
            17
        """
        cdef gen t0 = self(v)
        pari_catch_sig_on()
        return self.new_gen(polcyclo_eval(n, t0.g))

    def polsubcyclo(self, long n, long d, v=-1):
        """
        polsubcyclo(n, d, v=x): return the pari list of polynomial(s)
        defining the sub-abelian extensions of degree `d` of the
        cyclotomic field `\QQ(\zeta_n)`, where `d`
        divides `\phi(n)`.

        EXAMPLES::

            sage: pari.polsubcyclo(8, 4)
            [x^4 + 1]
            sage: pari.polsubcyclo(8, 2, 'z')
            [z^2 - 2, z^2 + 1, z^2 + 2]
            sage: pari.polsubcyclo(8, 1)
            [x - 1]
            sage: pari.polsubcyclo(8, 3)
            []
        """
        cdef gen plist
        pari_catch_sig_on()
        plist = self.new_gen(polsubcyclo(n, d, self.get_var(v)))
        if typ(plist.g) != t_VEC:
            return pari.vector(1, [plist])
        else:
            return plist
        #return self.new_gen(polsubcyclo(n, d, self.get_var(v)))

    def polzagier(self, long n, long m):
        pari_catch_sig_on()
        return self.new_gen(polzag(n, m))

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
            sage: a = pari.getrand(); a
            Vecsmall([...])
            sage: pari.setrand(a)
            sage: a == pari.getrand()
            True

        TESTS:

        Check that invalid inputs are handled properly (#11825)::

            sage: pari.setrand(0)
            Traceback (most recent call last):
            ...
            PariError: incorrect type in setrand
            sage: pari.setrand("foobar")
            Traceback (most recent call last):
            ...
            PariError: incorrect type in setrand
        """
        cdef gen t0 = self(seed)
        pari_catch_sig_on()
        setrand(t0.g)
        pari_catch_sig_off()

    def getrand(self):
        """
        Returns PARI's current random number seed.

        OUTPUT:

        GEN of type t_VECSMALL

        EXAMPLES::

            sage: pari.setrand(50)
            sage: a = pari.getrand(); a
            Vecsmall([...])
            sage: pari.setrand(a)
            sage: a == pari.getrand()
            True
        """
        pari_catch_sig_on()
        return self.new_gen(getrand())

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
                raise IndexError, "length of entries (=%s) must equal n (=%s)"%\
                      (len(entries), n)
            for i, x in enumerate(entries):
                v[i] = x
        return v

    cdef gen _empty_vector(self, long n):
        cdef gen v
        pari_catch_sig_on()
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

        pari_catch_sig_on()
        A = self.new_gen(zeromatcopy(m,n))
        if entries is not None:
            if len(entries) != m*n:
                raise IndexError, "len of entries (=%s) must be %s*%s=%s"%(len(entries),m,n,m*n)
            k = 0
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    x = pari(entries[k])
                    A._refers_to[(i,j)] = x
                    (<GEN>(A.g)[j+1])[i+1] = <long>(x.g)
                    k = k + 1
        return A


cdef int init_stack(size_t requested_size) except -1:
    r"""
    Low-level Cython function to allocate the PARI stack.  This
    function should not be called directly, use ``pari.allocatemem()``
    instead.
    """
    global top, bot, avma

    cdef size_t old_size = <size_t>(top) - <size_t>(bot)

    cdef size_t new_size
    cdef size_t max_size = <size_t>(-1)
    if (requested_size == 0):
        if old_size < max_size/2:
            # Double the stack
            new_size = 2*old_size
        elif old_size < 4*(max_size/5):
            # We cannot possibly double since we already use over half
            # the addressable memory: take the average of current and
            # maximum size
            new_size = max_size/2 + old_size/2
        else:
            # We already use 80% of the addressable memory => give up
            raise MemoryError("Unable to enlarge PARI stack (instead, kept the stack at %s bytes)"%(old_size))
    else:
        new_size = requested_size

    # Align size to 16 bytes and take 1024 bytes as a minimum
    new_size = (new_size/16)*16
    if (new_size < 1024):
        new_size = 1024

    # Disable interrupts
    sig_on()
    sig_block()

    # If this is non-zero, the size we failed to allocate
    cdef size_t failed_size = 0

    try:
        # Free the current stack
        if bot:
            libc.stdlib.free(<void*>bot)

        # Allocate memory for new stack.
        bot = <pari_sp> libc.stdlib.malloc(new_size)

        # If doubling failed, instead add 25% to the current stack size.
        # We already checked that we use less than 80% of the maximum value
        # for s, so this will not overflow.
        if (bot == 0) and (requested_size == 0):
            new_size = (old_size/64)*80
            bot = <pari_sp> libc.stdlib.malloc(new_size)

        if not bot:
            failed_size = new_size
            # We lost our PARI stack and are not able to allocate the
            # requested size.  If we just raise an exception now, we end up
            # *without* a PARI stack which is not good.  We will raise an
            # exception later, after allocating *some* PARI stack.
            new_size = old_size
            while new_size >= 1024:  # hope this never fails!
                bot = <pari_sp> libc.stdlib.malloc(new_size)
                if bot: break
                new_size = (new_size/32)*16

        if not bot:
            top = 0
            avma = 0
            raise SystemError("Unable to allocate PARI stack, all subsequent PARI computations will crash")

        top = bot + new_size
        avma = top

        if failed_size:
            raise MemoryError("Unable to allocate %s bytes for the PARI stack (instead, allocated %s bytes)"%(failed_size, new_size))

        return 0
    finally:
        sig_unblock()
        sig_off()


cdef GEN deepcopy_to_python_heap(GEN x, pari_sp* address):
    cdef size_t s = <size_t> gsizebyte(x)
    cdef pari_sp tmp_bot, tmp_top

    tmp_bot = <pari_sp> sage_malloc(s)
    tmp_top = tmp_bot + s
    address[0] = tmp_bot
    return gcopy_avma(x, &tmp_top)

cdef gen _new_gen(GEN x):
    cdef GEN h
    cdef pari_sp address
    cdef gen y
    h = deepcopy_to_python_heap(x, &address)
    y = PY_NEW(gen)
    y.init(h, address)
    return y
