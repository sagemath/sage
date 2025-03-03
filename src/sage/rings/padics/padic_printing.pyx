# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"""
`p`-adic Printing

This file contains code for printing `p`-adic elements.

It has been moved here to prevent code duplication and make finding
the relevant code easier.

AUTHORS:

- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2008-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.list cimport *
from sage.libs.gmp.mpz cimport *
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool
from sage.rings.padics.misc import trim_zeros
from sage.misc.latex import latex_variable_name

import sys

from sage.rings.integer cimport Integer

cdef enum print_modes:
    terse
    series
    val_unit
    digits
    bars


def pAdicPrinter(ring, options={}):
    """
    Create a :class:`pAdicPrinter`.

    INPUT:

    - ``ring`` -- a `p`-adic ring or field

    - ``options`` -- dictionary, with keys in ``'mode'``, ``'pos'``,
      ``'ram_name'``, ``'unram_name'``, ``'var_name'``, ``'max_ram_terms'``,
      ``'max_unram_terms'``, ``'max_terse_terms'``, ``'sep'``, ``'alphabet'``; see
      :class:`pAdicPrinter_class` for the meanings of these keywords.

    EXAMPLES::

        sage: from sage.rings.padics.padic_printing import pAdicPrinter
        sage: R = Zp(5)
        sage: pAdicPrinter(R, {'sep': '&'})
        series printer for 5-adic Ring with capped relative precision 20
    """
    for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet', 'show_prec']:
        if option not in options:
            options[option] = None
    return pAdicPrinter_class(ring, **options)


class pAdicPrinterDefaults(SageObject):
    """
    This class stores global defaults for `p`-adic printing.
    """
    def __init__(self, mode='series', pos=True, max_ram_terms=-1, max_unram_terms=-1, max_terse_terms=-1, sep="|", alphabet=None):
        r"""
        Instances of this class store global defaults used in
        determining printing options during the creation of `p`-adic
        rings and fields.  One instance stored in padic_printing
        stores the globally relevant default values.

        See pAdicPrinter_class for details on the meanings of these
        inputs.

        TESTS::

            sage: from sage.rings.padics.padic_printing import pAdicPrinterDefaults
            sage: D = pAdicPrinterDefaults(sep='&'); D.sep()
            '&'
        """
        self._mode = mode
        self._pos = bool(pos)
        if not -1 <= max_ram_terms <= sys.maxsize:
            raise ValueError("max_ram_terms must be positive and fit in a long")
        self._max_ram_terms = int(max_ram_terms)
        if not -1 <= max_unram_terms <= sys.maxsize:
            raise ValueError("max_unram_terms must be positive and fit in a long")
        self._max_unram_terms = int(max_unram_terms)
        if not -1 <= max_terse_terms <= sys.maxsize:
            raise ValueError("max_terse_terms must be positive and fit in a long")
        self._max_terse_terms = int(max_terse_terms)
        self._sep = sep
        if alphabet is None:
            self._alphabet = ('0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z')
        else:
            self._alphabet = alphabet

    def mode(self, mode=None):
        r"""
        Set the default printing mode.

        ``mode=None`` returns the current value.

        The allowed values for mode are: ``'val-unit'``, ``'series'``,
        ``'terse'``, ``'digits'`` and ``'bars'``.

        EXAMPLES::

            sage: padic_printing.mode('terse')
            sage: padic_printing.mode()
            'terse'
            sage: Qp(7)(100)
            100 + O(7^20)
            sage: padic_printing.mode('series')
            sage: Qp(11)(100)
            1 + 9*11 + O(11^20)
            sage: padic_printing.mode('val-unit')
            sage: Qp(13)(130)
            13 * 10 + O(13^21)
            sage: padic_printing.mode('digits')
            sage: repr(Qp(17)(100))
            '...5F'
            sage: repr(Qp(17)(1000))
            '...37E'
            sage: padic_printing.mode('bars')
            sage: repr(Qp(19)(1000))
            '...2|14|12'

            sage: padic_printing.mode('series')
        """
        if mode is None:
            return self._mode
        else:
            if mode in ['val-unit','series','terse','digits','bars']:
                self._mode = mode
            else:
                raise ValueError("invalid printing mode")

    def allow_negatives(self, neg=None):
        r"""
        Controls whether or not to display a balanced representation.

        ``neg=None`` returns the current value.

        EXAMPLES::

            sage: padic_printing.allow_negatives(True)
            sage: padic_printing.allow_negatives()
            True
            sage: Qp(29)(-1)
            -1 + O(29^20)
            sage: Qp(29)(-1000)
            -14 - 5*29 - 29^2 + O(29^20)
            sage: padic_printing.allow_negatives(False)
        """
        if neg is None:
            return not self._pos
        else:
            self._pos = not neg

    def max_series_terms(self, max=None):
        r"""
        Controls the maximum number of terms shown when printing in
        ``'series'``, ``'digits'`` or ``'bars'`` mode.

        ``max=None`` returns the current value.

        ``max=-1`` encodes 'no limit.'

        EXAMPLES::

            sage: padic_printing.max_series_terms(2)
            sage: padic_printing.max_series_terms()
            2
            sage: Qp(31)(1000)
            8 + 31 + ... + O(31^20)
            sage: padic_printing.max_series_terms(-1)
            sage: Qp(37)(100000)
            26 + 37 + 36*37^2 + 37^3 + O(37^20)
        """
        if max is None:
            return self._max_ram_terms
        else:
            self._max_ram_terms = int(max)

    def max_unram_terms(self, max=None):
        r"""
        For rings with non-prime residue fields, controls how many
        terms appear in the coefficient of each ``pi^n`` when printing in
        ``'series'`` or ``'bar'`` modes.

        ``max=None`` returns the current value.

        ``max=-1`` encodes 'no limit.'

        EXAMPLES::

            sage: padic_printing.max_unram_terms(2)
            sage: padic_printing.max_unram_terms()
            2
            sage: Zq(5^6, 5, names='a')([1,2,3,-1])^17                                  # needs sage.libs.ntl
            (3*a^4 + ... + 3) + (a^5 + ... + a)*5 + (3*a^3 + ... + 2)*5^2 + (3*a^5 + ... + 2)*5^3 + (4*a^5 + ... + 4)*5^4 + O(5^5)

            sage: padic_printing.max_unram_terms(-1)
        """
        if max is None:
            return self._max_unram_terms
        else:
            self._max_unram_terms = int(max)

    def max_poly_terms(self, max=None):
        r"""
        Controls the number of terms appearing when printing
        polynomial representations in ``'terse'`` or ``'val-unit'`` modes.

        ``max=None`` returns the current value.

        ``max=-1`` encodes 'no limit.'

        EXAMPLES::

            sage: padic_printing.max_poly_terms(3)
            sage: padic_printing.max_poly_terms()
            3
            sage: padic_printing.mode('terse')
            sage: Zq(7^5, 5, names='a')([2,3,4])^8                                      # needs sage.libs.ntl
            2570 + 15808*a + 9018*a^2 + ... + O(7^5)

            sage: padic_printing.max_poly_terms(-1)
            sage: padic_printing.mode('series')
        """
        if max is None:
            return self._max_terse_terms
        else:
            self._max_terse_terms = int(max)

    def sep(self, sep=None):
        r"""
        Controls the separator used in ``'bars'`` mode.

        ``sep=None`` returns the current value.

        EXAMPLES::

            sage: padic_printing.sep('][')
            sage: padic_printing.sep()
            ']['
            sage: padic_printing.mode('bars')
            sage: repr(Qp(61)(-1))
            '...60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60'

            sage: padic_printing.sep('|')
            sage: padic_printing.mode('series')
        """
        if sep is None:
            return self._sep
        else:
            self._sep = str(sep)

    def alphabet(self, alphabet=None):
        r"""
        Controls the alphabet used to translate `p`-adic digits into
        strings (so that no separator need be used in ``'digits'`` mode).

        ``alphabet`` should be passed in as a list or tuple.

        ``alphabet=None`` returns the current value.

        EXAMPLES::

            sage: padic_printing.alphabet("abc")
            sage: padic_printing.mode('digits')
            sage: repr(Qp(3)(1234))
            '...bcaacab'

            sage: padic_printing.mode('series')
            sage: padic_printing.alphabet(('0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'))
        """
        if alphabet is None:
            return self._alphabet
        else:
            self._alphabet = tuple(alphabet)


_printer_defaults = pAdicPrinterDefaults()


cdef class pAdicPrinter_class(SageObject):
    """
    This class stores the printing options for a specific `p`-adic ring
    or field, and uses these to compute the representations of
    elements.
    """
    def __init__(self, ring, mode, pos, ram_name, unram_name, var_name, max_ram_terms, max_unram_terms, max_terse_terms, sep, alphabet, show_prec):
        """
        Initialize a :class:`pAdicPrinter`.

        INPUT:

        - ``ring`` -- the ring or field to which this :class:`pAdicPrinter` is
          attached

        - ``mode`` -- the allowed values for mode are: ``'val-unit'``,
          ``'series'``, ``'terse'``, ``'digits'`` and ``'bars'``:

          - ``'val-unit'`` -- elements are displayed as a power of the
            uniformizer times a unit, which is displayed in terse mode

          - ``'series'`` -- elements are displayed as power series in the
            uniformizer

          - ``'terse'`` -- for base rings and fields, elements are just
            displayed as an integer lift.  For extensions rings and fields,
            elements are displayed as a polynomial in the generator of the
            extension.

          - ``'digits'`` -- used only for small primes and totally ramified
            extensions (or trivial extensions), elements are displayed as just
            a string of `p`-adic digits, encoded using the 'alphabet' parameter

          - ``'bars'`` -- like ``'digits'``, but uses a separator in order to
            print a more canonical representation for each digit. This change
            allows the use of this printing mode for unramified extensions and
            extensions with larger primes.

        - ``pos`` -- if ``True`` then integers in the range [0,... p-1] will be
          used; if ``False`` integers in the range [(1-p)/2,..., p/2] will be used

        - ``ram_name`` -- the string used to represent the uniformizer

        - ``unram_name`` -- the string used to represent the trivial lift of a
          generator of the residue field over the prime field

        - ``var_name`` -- the string used to represent the user-specified
          generator of this extension ring or field

        - ``max_ram_terms`` -- controls the maximum number of terms shown when
          printing in ``'series'``, ``'digits'`` or ``'bars'`` mode

        - ``max_unram_terms`` -- for rings with non-prime residue fields,
          controls how many terms appear in the coefficient of each pi^n when
          printing in ``'series'`` or ``'bar'`` modes

        - ``max_terse_terms`` -- controls the number of terms appearing when
          printing polynomial representations in ``'terse'`` or ``'val-unit'``
          modes

        - ``sep`` -- controls the separator used in ``'bars'`` mode

        - ``alphabet`` -- controls the alphabet used to translate `p`-adic digits
          into strings (so that no separator need be used in ``'digits'`` mode)

        - ``show_prec`` -- Specify how the precision is printed; it can be
          ``'none'``, ``'bigoh'`` or ``'dots'`` (the latter being not available
          for all modes)

        TESTS::

            sage: R = Qp(7, print_mode='bars', print_sep='&')  # indirect doctest

            sage: R = Zp(5, print_mode='digits', print_max_terms=10)
            Traceback (most recent call last):
            ...
            ValueError: max_ram_terms must be unset when show_prec is 'dots'
        """
        global _printer_defaults
        self.ring = ring
        self.prime_pow = ring.prime_pow
        from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
        self.base = isinstance(ring, pAdicBaseGeneric)
        if alphabet is None:
            self.alphabet = _printer_defaults._alphabet
        else:
            self.alphabet = alphabet
        # note that self.pos is reset to True if mode == 'digits'
        if pos is None:
            self.pos = _printer_defaults._pos
        else:
            self.pos = pos
        if mode is None:
            mode = _printer_defaults._mode
        if mode == 'val-unit':
            self.mode = val_unit
        elif mode == 'series':
            self.mode = series
        elif mode == 'terse':
            self.mode = terse
        elif mode == 'digits':
            if len(self.alphabet) < self.prime_pow.prime or (not self.base and ring.absolute_f() != 1):
                raise ValueError("digits printing mode only usable for totally ramified extensions with p at most the length of the alphabet (default 62).  Try using print_mode = 'bars' instead.")
            else:
                self.mode = digits
                self.pos = True
        elif mode == 'bars':
            self.mode = bars
        else:
            raise ValueError("printing mode must be one of 'val-unit', 'series', 'terse', 'digits' or 'bars'")
        if ram_name is None:
            self.ram_name = ring._uniformizer_print()
        else:
            self.ram_name = ram_name
        if self.ram_name is None:
            self.latex_ram_name = None
        else:
            self.latex_ram_name = latex_variable_name(self.ram_name)
        if unram_name is None:
            self.unram_name = ring._unram_print()
        else:
            self.unram_name = unram_name
        if self.unram_name is None:
            self.latex_unram_name = None
        else:
            self.latex_unram_name = latex_variable_name(self.unram_name)
        if var_name is None:
            self.var_name = ring.variable_name()
        else:
            self.var_name = var_name
        if self.var_name is None:
            self.latex_var_name = None
        else:
            self.latex_var_name = latex_variable_name(self.var_name)
        if sep is None:
            self.sep = _printer_defaults._sep
        else:
            self.sep = sep
        if max_ram_terms is not None:
            self.max_ram_terms = max_ram_terms
            if self.max_ram_terms < -1:
                raise ValueError("max_ram_terms must be positive and fit in a long")
        else:
            self.max_ram_terms = _printer_defaults._max_ram_terms
        if max_unram_terms is not None:
            self.max_unram_terms = max_unram_terms
            if self.max_unram_terms < -1:
                raise ValueError("max_unram_terms must be positive and fit in a long")
        else:
            self.max_unram_terms = _printer_defaults._max_unram_terms
        if max_terse_terms is not None:
            self.max_terse_terms = max_terse_terms
            if self.max_terse_terms < -1:
                raise ValueError("max_terse_terms must be positive and fit in a long")
        else:
            self.max_terse_terms = _printer_defaults._max_terse_terms
        from sage.rings.padics.factory import _canonicalize_show_prec
        self.show_prec = _canonicalize_show_prec(self.ring._prec_type(), mode, show_prec)

        # Incompatibilities
        if self.show_prec == 'dots' and max_ram_terms >= 0:
            raise ValueError("max_ram_terms must be unset when show_prec is 'dots'")

    def __reduce__(self):
        """
        Pickling.

        TESTS::

            sage: R = Zp(5, print_mode='bars', print_sep='&'); P = loads(dumps(R._printer))
            sage: R._printer == P
            True
            sage: P._sep()
            '&'
        """
        return pAdicPrinter, (self.ring,
                              {'mode': self._print_mode(),
                               'pos': self.pos,
                               'ram_name': self.ram_name,
                               'unram_name': self.unram_name,
                               'var_name': self.var_name,
                               'max_ram_terms': self.max_ram_terms,
                               'max_unram_terms': self.max_unram_terms,
                               'max_terse_terms': self.max_terse_terms,
                               'sep':self.sep,
                               'alphabet': self.alphabet,
                               'show_prec': self.show_prec})

    def __richcmp__(self, other, op):
        """
        Comparison.

        TESTS::

            sage: R = Zp(5)
            sage: S = Zp(5,print_mode='bars')
            sage: R._printer == S._printer
            False
        """
        if not isinstance(other, pAdicPrinter_class):
            return NotImplemented
        return self.richcmp_modes(other, op)

    def richcmp_modes(pAdicPrinter_class self,
                      pAdicPrinter_class other, int op):
        """
        Return a comparison of the printing modes of ``self`` and ``other``.

        Return 0 if and only if all relevant modes are equal
        (``max_unram_terms`` is irrelevant if the ring is totally ramified
        over the base, for example). This does not check if the rings are
        equal (to prevent infinite recursion in the comparison
        functions of `p`-adic rings), but it does check if the primes
        are the same (since the prime affects whether ``pos`` is
        relevant).

        EXAMPLES::

            sage: R = Qp(7, print_mode='digits', print_pos=True)
            sage: S = Qp(7, print_mode='digits', print_pos=False)
            sage: R._printer == S._printer
            True
            sage: R = Qp(7)
            sage: S = Qp(7, print_mode='val-unit')
            sage: R == S
            False
            sage: R._printer < S._printer
            True
        """
        lx = self.mode
        rx = other.mode
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        p = self.ring.prime()
        q = other.ring.prime()
        if p != q:
            return richcmp_not_equal(p, q, op)

        if p != 2 and (self.mode == terse or self.mode == series or self.mode == val_unit or self.mode == bars):
            lx = self.pos
            rx = other.pos
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

        if self.mode != digits:
            lx = self.ram_name
            rx = other.ram_name
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

        if self.mode == bars:
            lx = self.sep
            rx = other.sep
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

        if self.mode == digits:
            lx = self.alphabet[:p]
            rx = other.alphabet[:q]
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

        if self.mode == series or self.mode == digits or self.mode == bars:
            lx = self.max_ram_terms
            rx = other.max_ram_terms
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

        f = max(self.ring.absolute_f(), other.ring.absolute_f())

        if f > 1 and (self.mode == series or self.mode == bars):
            lx = self.unram_name
            rx = other.unram_name
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

            lx = self.max_unram_terms
            rx = other.max_unram_terms
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

        f = max(self.ring.relative_degree(), other.ring.relative_degree())

        if f > 1 and self.mode == terse:
            lx = self.var_name
            rx = other.var_name
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)
            lx = self.max_terse_terms
            rx = other.max_terse_terms
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)
            lx = self.show_prec
            rx = other.show_prec
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

        return rich_to_bool(op, 0)

    def _repr_(self):
        """
        Representation of this printer.

        EXAMPLES::

            sage: Zp(5)._printer  # indirect doctest
            series printer for 5-adic Ring with capped relative precision 20
        """
        return "%s printer for %s"%(self._print_mode(), self.ring)

    def __enter__(self):
        """
        Used for context printing.

        EXAMPLES::

            sage: from sage.rings.padics.padic_printing import pAdicPrinter
            sage: R = Zp(5,5); a = R(-1); a
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
            sage: with pAdicPrinter(R, {'pos': False}): a
            -1 + O(5^5)
        """
        self.old = self.ring._printer
        self.ring._printer = self

    def dict(self):
        """
        Return a dictionary storing all of ``self``'s printing options.

        EXAMPLES::

            sage: D = Zp(5)._printer.dict(); D['sep']
            '|'
        """
        return {'mode': self._print_mode(), 'pos': self.pos, 'ram_name': self.ram_name, 'unram_name': self.unram_name, 'var_name': self.var_name, 'max_ram_terms': self.max_ram_terms, 'max_unram_terms': self.max_unram_terms, 'max_terse_terms': self.max_terse_terms, 'sep': self.sep, 'alphabet': self.alphabet, 'show_prec': self.show_prec}

    def __exit__(self, type, value, traceback):
        """
        Used for context printing.

        EXAMPLES::

            sage: from sage.rings.padics.padic_printing import pAdicPrinter
            sage: R = Zp(5,5); a = R(-1); a
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
            sage: with pAdicPrinter(R, {'pos': False}): a
            -1 + O(5^5)
            sage: a
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
        """
        self.ring._printer = self.old

    def _pos(self):
        """
        Access ``self.pos``.

        EXAMPLES::

            sage: R = Zp(5); R._printer._pos()
            True
        """
        return self.pos

    def _sep(self):
        """
        Access ``self.sep``.

        EXAMPLES::

            sage: R = Zp(5); R._printer._sep()
            '|'
        """
        return self.sep

    def _alphabet(self):
        """
        Access ``self.pos``.

        EXAMPLES::

            sage: R = Zp(17, print_mode='digits'); R._printer._alphabet()
            ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G')
        """
        return self.alphabet

    def _max_ram_terms(self):
        """
        Access ``self.max_ram_terms``.

        EXAMPLES::

            sage: R = Zp(5); R._printer._max_ram_terms()
            -1
        """
        return self.max_ram_terms

    def _max_unram_terms(self):
        """
        Access ``self.max_unram_terms``.

        EXAMPLES::

            sage: R = Zp(5); R._printer._max_unram_terms()
            -1
        """
        return self.max_unram_terms

    def _max_terse_terms(self):
        """
        Access ``self.max_terse_terms``.

        EXAMPLES::

            sage: R = Zp(5); R._printer._max_terse_terms()
            -1
        """
        return self.max_terse_terms

    def _show_prec(self):
        """
        Access ``self.show_prec``.

        EXAMPLES::

            sage: R = ZpFP(5); R._printer._show_prec()
            'none'
        """
        return self.show_prec

    def _ring(self):
        """
        Access ``self.ring``.

        EXAMPLES::

            sage: R = Zp(5,5); R._printer._ring()
            5-adic Ring with capped relative precision 5
        """
        return self.ring

    def _uniformizer_name(self):
        """
        Access ``self.ram_name``.

        EXAMPLES::

            sage: R = Zp(5,5); R._printer._uniformizer_name()
            '5'
        """
        return self.ram_name

    def _ram_name(self):
        """
        Access ``self.ram_name``.

        EXAMPLES::

            sage: R = Zp(5,5); R._printer._ram_name()
            '5'
        """
        return self.ram_name

    def _print_mode(self):
        """
        Access ``self.mode``.

        EXAMPLES::

            sage: R = Zp(5); R._printer._print_mode()
            'series'
        """
        if self.mode == val_unit:
            return 'val-unit'
        elif self.mode == series:
            return 'series'
        elif self.mode == terse:
            return 'terse'
        elif self.mode == digits:
            return 'digits'
        elif self.mode == bars:
            return 'bars'

    def _base_p_list(self, value, pos):
        """
        Return a list of integers forming the base p expansion of value.

        If pos is True, these integers will be in the range
        [0,..., p-1]; if pos is False, they will be in the range
        [(1-p)/2,..., p/2].

        The first entry will be the coefficient of p^0, etc.

        EXAMPLES::

            sage: P = Zp(17)._printer
            sage: P._base_p_list(1298734,True)
            [2, 15, 5, 9, 15]
            sage: P._base_p_list(1298734,False)
            [2, -2, 6, -8, -1, 1]
            sage: P._base_p_list(Zp(17)(1298734),True)
            [2, 15, 5, 9, 15]
            sage: P._base_p_list(Zp(17)(1298734),False)
            [2, -2, 6, -8, -1, 1]
        """
        return self.base_p_list(value, pos)

    cdef base_p_list(self, value, bint pos):
        """
        Return a list of integers forming the base p expansion of value.

        If pos is True, these integers will be in the range
        [0,... p-1]; if po is False, they will be in the range
        [(1-p)/2,..., p/2].

        The first entry will be the coefficient of p^0, etc.

        EXAMPLES::

            sage: P = Zp(17)._printer
            sage: P._base_p_list(1298734,True)  # indirect doctest
            [2, 15, 5, 9, 15]
            sage: P._base_p_list(1298734,False)
            [2, -2, 6, -8, -1, 1]
        """
        if isinstance(value, Integer):
            from sage.rings.padics.padic_capped_relative_element import base_p_list
            return base_p_list(value, pos, self.prime_pow)
        elif pos:
            return trim_zeros(list(value.unit_part().expansion()))
        else:
            return trim_zeros(list(value.unit_part().expansion(lift_mode='smallest')))

    def repr_gen(self, elt, do_latex, pos=None, mode=None, ram_name=None):
        """
        The entry point for printing an element.

        INPUT:

        - ``elt`` -- a `p`-adic element of the appropriate ring to print

        - ``do_latex`` -- whether to return a latex representation or
          a normal one

        EXAMPLES::

            sage: R = Zp(5,5); P = R._printer; a = R(-5); a
            4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6)
            sage: P.repr_gen(a, False, pos=False)
            '-5 + O(5^6)'
            sage: P.repr_gen(a, False, ram_name='p')
            '4*p + 4*p^2 + 4*p^3 + 4*p^4 + 4*p^5 + O(p^6)'
        """
        cdef int _mode
        cdef bint _pos
        if mode is None:
            _mode = self.mode
        elif mode == 'val-unit':
            _mode = val_unit
        elif mode == 'series':
            _mode = series
        elif mode == 'terse':
            _mode = terse
        elif mode == 'digits':
            _mode = digits
        elif mode == 'bars':
            _mode = bars
        else:
            raise ValueError("printing mode must be one of 'val-unit', 'series', 'terse', 'bars', or 'digits'")
        if pos is None:
            _pos = self.pos
        else:
            _pos = pos
        if ram_name is None:
            if do_latex:
                pprint = self.latex_ram_name
            else:
                pprint = self.ram_name
        else:
            pprint = str(ram_name)
            if do_latex:
                pprint = latex_variable_name(pprint)
        return self._repr_gen(elt, do_latex, _pos, _mode, pprint)

    cdef _repr_gen(self, pAdicGenericElement elt, bint do_latex, bint pos, int mode, ram_name):
        r"""
        Print a string representation of the element.  See ``__init__`` for
        more details on print modes.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-rel','val-unit'); a = R(364); a  # indirect doctest
            7 * 52 + O(7^5)
            sage: print(a.str('terse'))
            364 + O(7^5)
            sage: print(a.str('series'))
            3*7 + 7^3 + O(7^5)
            sage: K = Qp(7,4,'capped-rel','val-unit'); a = K(364); a
            7 * 52 + O(7^5)
            sage: print(a.str('series'))
            3*7 + 7^3 + O(7^5)
            sage: padic_printing.sep('')
            sage: K = Qp(7, print_mode='digits')
            sage: repr(K(1/2))
            '...33333333333333333334'
            sage: repr(K(1/42))
            '...5555555555555555555.6'
            sage: padic_printing.sep('|')
            sage: repr(Qp(97, print_mode='bars')(1/13))
            '...29|82|7|44|74|59|67|14|89|52|22|37|29|82|7|44|74|59|67|15'

        TESTS:

        Check that :issue:`24843` is resolved::

            sage: R = Zp(2, print_mode='digits', show_prec=True)
            sage: repr(R(0,10))
            '...0000000000'
        """
        s = ""
        if self.show_prec == "dots":
            unknown_digit = "?"
        else:
            unknown_digit = "0"
        if elt._is_exact_zero():
            return "0"
        elif elt._is_inexact_zero():
            prec = elt.precision_absolute()
            if mode == digits and self.show_prec == "dots":
                if prec > 0:
                    s = (self.alphabet[0] * prec)
                else:
                    s = unknown_digit + "." + (unknown_digit * (-prec)) + self.alphabet[0]
            elif mode == bars and self.show_prec == "dots":
                if self.base or self._ring().absolute_f() == 1:
                    zero = '0'
                else:
                    zero = '[]'
                if prec > 0:
                    L = [zero] * prec
                    s = self.sep.join(L)
                else:
                    L = ['.'] + ([unknown_digit] * (-prec)) + [zero]
                    s = self.sep + self.sep.join(L)
        elif mode == val_unit:
            if do_latex:
                if elt.valuation() == 0:
                    s = self._repr_spec(elt, do_latex, pos, terse, 0, ram_name)
                elif elt.valuation() == 1:
                    s = "%s \\cdot %s"%(ram_name, self._repr_spec(elt.unit_part(), do_latex, pos, terse, 1, ram_name))
                else:
                    s = "%s^{%s} \\cdot %s"%(ram_name, elt.valuation(), self._repr_spec(elt.unit_part(), do_latex, pos, terse, 1, ram_name))
            else:
                if elt.valuation() == 0:
                    s = self._repr_spec(elt, do_latex, pos, terse, 0, ram_name)
                elif elt.valuation() == 1:
                    s = "%s * %s"%(ram_name, self._repr_spec(elt.unit_part(), do_latex, pos, terse, 1, ram_name))
                else:
                    s = "%s^%s * %s"%(ram_name, elt.valuation(), self._repr_spec(elt.unit_part(), do_latex, pos, terse, 1, ram_name))
        elif mode == digits:
            n = elt.valuation()
            if self.base:
                L = self.base_p_list(elt, True)
            else:
                L = elt._ext_p_list(True)
            if self.max_ram_terms == -1:
                lenL = elt.precision_relative()
            else:
                lenL = min(elt.precision_relative(), max(self.max_ram_terms, -n))
            if self.show_prec == "dots":
                if len(L) < lenL:
                    L += [0] * (lenL - len(L))
            if len(L) > lenL:
                L = L[:lenL]
            L.reverse()
            # The following step should work since mode is only allowed to be digits in the case of totally ramified extensions
            # with primes smaller than the length of the alphabet
            L = [self.alphabet[a] for a in L]
            if n > 0:
                L += [self.alphabet[0]]*n
            elif n < 0:
                L = [unknown_digit]*(1 - n - lenL) + L
                L = L[:n] + ['.'] + L[n:]
            s = "".join(L)
        elif mode == bars:
            n = elt.valuation()
            if self.base:
                L = self.base_p_list(elt, self.pos)
            else:
                L = elt._ext_p_list(self.pos)
            if self.max_ram_terms == -1:
                lenL = elt.precision_relative()
            else:
                lenL = min(elt.precision_relative(), max(self.max_ram_terms, -n))
            if self.show_prec == "dots":
                if len(L) < lenL:
                    if self.base or self._ring().absolute_f() == 1:
                        L += [0]*(lenL - len(L))
                    else:
                        L += [[]]*(lenL - len(L))
            if len(L) > lenL:
                L = L[:lenL]
            L.reverse()
            if self.base or self._ring().absolute_f() == 1 or self.max_unram_terms == -1:
                L = [str(a) for a in L]
            else:
                if self.max_unram_terms == 0:
                    L = ['[...]' if len(a) > 0 else '[]' for a in L]
                elif self.max_unram_terms == 1:
                    L = ["[..., %s]"%(a[-1]) if len(a) > 1 else str(a) for a in L]
                else:
                    L = ["[%s,..., " % (a[0]) + ", ".join(str(b) for b in a[1-self.max_unram_terms:]) + "]" if len(a) > 2 else str(a) for a in L]
            if n > 0:
                if self.base or self._ring().absolute_f() == 1:
                    L += ['0']*n
                else:
                    L += ['[]']*n
            elif n < 0:
                L = ['0']*(min(-n, elt.precision_relative()) - len(L)) + L
                L = [unknown_digit]*(-n - len(L)) + L
                L = L[:n] + ['.'] + L[n:]
            if L[0] == '.':
                s = self.sep + self.sep.join(L)
            else:
                s = self.sep.join(L)
        else: # mode == terse or series
            s = self._repr_spec(elt, do_latex, pos, mode, 0, ram_name)

        if self.show_prec == "dots":
            s = "..." + s
        if self.show_prec == "bigoh":
            if s == "":
                s = "O(%s" % ram_name
            else:
                s += " + O(%s" % ram_name
            if elt.precision_absolute() == 1:
                s += ")"
            else:
                if do_latex:
                    s += "^{%s})"%(elt.precision_absolute())
                else:
                    s += "^%s)"%(elt.precision_absolute())
        if s == "": s = "0"
        return s

    cdef _repr_spec(self, pAdicGenericElement elt, bint do_latex, bint pos, int mode, bint paren, ram_name):
        """
        A function used by repr_gen for terse and series printing.

        Should not be called if elt is an exact or inexact zero.
        """
        cdef Integer lift_z, pprec
        cdef int ZZ_pEX
        cdef Py_ssize_t i
        cdef long val
        cdef bint integral
        var_name = self.latex_var_name if do_latex else self.var_name
        if self.base:
            if mode == terse:
                v = elt.valuation()
                if v >= 0:
                    lift_z = <Integer> elt.lift()
                    pprec = self.prime_pow.pow_Integer(mpz_get_ui((<Integer>elt.precision_absolute()).value))
                else:
                    lift_z = <Integer> elt.unit_part().lift()
                    pprec = self.prime_pow.pow_Integer(mpz_get_ui((<Integer>elt.precision_relative()).value))
                mpz_mod(lift_z.value, lift_z.value, pprec.value)
                if not pos:
                    if lift_z > pprec / 2:
                        mpz_sub(lift_z.value, lift_z.value, pprec.value)
                    else:
                        paren = 0
                else:
                    paren = 0
                # if v>=0, _terse_frac only uses the first input.
                # if v<0, _terse_frac doesn't use the first input at all, and expects the unit part in the third input.
                s = self._terse_frac(lift_z, v, lift_z, ram_name, do_latex)
                if paren and not do_latex:
                    return "(%s)"%(s)
                else:
                    return s
            else: # mode == series
                slist = self.base_p_list(elt, pos)
                slist, ellipsis = self._truncate_list(slist, self.max_ram_terms, 0)
                s = ""
                exp = elt.valuation()
                for a in slist:
                    if a != 0:
                        if a < 0:
                            if len(s) == 0:
                                s = "-"
                            else:
                                s += " - "
                            a = -a
                        elif len(s) != 0:
                            s += " + "
                        s += self._co_dot_var(a, ram_name, exp, do_latex)
                    exp += 1
                if ellipsis:
                    s += self._plus_ellipsis(do_latex)
        else: # not self.base
            if mode == terse:
                if elt.parent()._implementation == 'FLINT':
                    poly, k = elt._flint_rep_abs()
                    L = [repr(a) for a in poly.coefficients(sparse=False)]
                    ZZ_pEX = 1
                elif elt.parent()._implementation == 'Polynomial':
                    poly = elt._poly_rep()
                    if do_latex:
                        L = [a._latex_() for a in poly.coefficients(sparse=False)]
                    else:
                        L = [repr(a) for a in poly.coefficients(sparse=False)]
                    L, ellipsis = self._truncate_list(L, self.max_terse_terms, '0')
                    s = ''
                    for i, a in enumerate(L):
                        if a == '0':
                            continue
                        if (i > 0) and ('+' in a or '-' in a[1:]):
                            a = '(' + a + ')'
                        if a[0] == '-':
                            if s:
                                s += ' - '
                            else:
                                s = '-'
                            if a[1:] == '1':
                                s += self._var(var_name, i, do_latex)
                            else:
                                s += a[1:] + self._dot_var(var_name, i, do_latex)
                        else:
                            if s:
                                s += ' + '
                            if a == '1':
                                s += self._var(var_name, i, do_latex)
                            else:
                                s += a + self._dot_var(var_name, i, do_latex)
                    if ellipsis:
                        s += self._plus_ellipsis(do_latex)
                    if paren and ' ' in s:
                        s = '(' + s + ')'
                    return s
                else:
                    if elt.parent().is_capped_relative():
                        poly, k = elt._ntl_rep_abs()
                        s = repr(poly)
                    else:
                        s = repr(elt._ntl_rep())
                        k = 0
                    L = s.split("] [") # this splits a ZZ_pEX into the ZZ_pX components
                    ZZ_pEX = L[0].count("[") # will equal 2 if elt was a ZZ_pEX element, 1 if it was a ZZ_pX element
                    L[0] = L[0].replace("[","")
                    L[-1] = L[-1].replace("]","")
                    L = L[0].split()
                if ZZ_pEX == 2:
                    L = [a.split() for a in L]
                    L = [[("" if b == "0" else b) for b in a] for a in L]
                    L, ellipsis = self._truncate_list(L, self.max_ram_terms, "")
                    raise NotImplementedError
                else:
                    L = [("" if b == "0" else b) for b in L]
                    L, ellipsis = self._truncate_list(L, self.max_terse_terms, "")
                    s = ""
                    pn = self.prime_pow.pow_Integer(mpz_get_ui((<Integer>(elt.precision_absolute()-k)).value))
                    if elt.parent().is_capped_relative():
                        pk = self.prime_pow.pow_Integer(mpz_get_ui((<Integer>-k).value))
                        if k >= 0:
                            integral = True
                        else:
                            integral = False
                    else:
                        pk = Integer(1)
                        integral = True
                    for i from 0 <= i < len(L):
                        if L[i] != "":
                            a = Integer(L[i])
                            if not pos and 2*a > pn:
                                if integral:
                                    arep = pn - a
                                else:
                                    a = (pn - a) / pk
                                    v, u = a.val_unit(self.prime_pow.prime)
                                    arep = self._terse_frac(a, v, u, ram_name, do_latex)
                                if s == "":
                                    s = "-%s"%(arep)
                                    s += self._dot_var(var_name, i, do_latex)
                                elif a == 1:
                                    s += " - "
                                    s += self._var(var_name, i, do_latex)
                                else:
                                    s += " - %s" % (arep)
                                    s += self._dot_var(var_name, i, do_latex)
                            elif a == pk:
                                if s:
                                    s += " + "
                                s += self._var(var_name, i, do_latex)
                            else:
                                a = a / pk
                                v, u = a.val_unit(self.prime_pow.prime)
                                arep = self._terse_frac(a, v, u, ram_name, do_latex)
                                if s == "":
                                    s = "%s" % arep
                                else:
                                    s += " + %s" % arep
                                s += self._dot_var(var_name, i, do_latex)
                    if ellipsis:
                        s += self._plus_ellipsis(do_latex)
            else: # series
                s = ""
                L = elt._ext_p_list(pos)
                val = elt.valuation_c()
                # since elt was not supposed to be zero, this should give a non-empty list.
                if len(L) == 0:
                    raise RuntimeError("repr_spec called on zero")
                R = elt.parent()
                f = R.absolute_f()
                if f > 1: # unramified part to the extension
                    if self.unram_name is None:
                        raise RuntimeError("need to have specified a name for the unramified variable")
                    L, ellipsis = self._truncate_list(L, self.max_ram_terms, [])
                    for i from 0 <= i < len(L):
                        unram_name = self.latex_unram_name if do_latex else self.unram_name
                        term = self._print_unram_term(L[i], do_latex, unram_name, self.max_unram_terms, 0, 0)
                        if len(term) > 0:
                            exp = i + val
                            if (not do_latex and term.find(" ") != -1) or (do_latex and (term.find(" + ") != -1 or term.find(" - ") != -1)):
                                if len(s) > 0:
                                    s += " + (" + term + ")"
                                else:
                                    s = "(" + term + ")"
                                s += self._dot_var(ram_name, exp, do_latex)
                            else: # since len(term) > 0, len(L[i]) == 1.
                                if term == "1":
                                    if len(s) > 0:
                                        s += " + "
                                    s += self._var(ram_name, exp, do_latex)
                                    continue
                                elif term == "-1":
                                    if len(s) > 0:
                                        s += " - "
                                    else:
                                        s = "-"
                                    s += self._var(ram_name, exp, do_latex)
                                    continue
                                elif len(s) > 0:
                                    if term[0] == "-":
                                        s += " - " + term[1:]
                                    else:
                                        s += " + " + term
                                    s += self._dot_var(ram_name, exp, do_latex)
                                else:
                                    s = term
                                    s += self._dot_var(ram_name, exp, do_latex)
                    if ellipsis:
                        s += self._plus_ellipsis(do_latex)
                else: # f = 1, so no unramified printing required
                    L, ellipsis = self._truncate_list(L, self.max_ram_terms, 0)
                    s = self._print_list_as_poly(L, do_latex, ram_name, val, 1)
                    if ellipsis:
                        s += self._plus_ellipsis(do_latex)
        if paren and s.find(" ") != -1:
            s = "(" + s + ")"
        return s

    cdef _var(self, x, exp, do_latex):
        """
        Return a representation of 'x^exp', latexed if necessary.
        """
        if exp == 0:
            return "1"
        if exp == 1:
            return str(x)
        if do_latex:
            return "%s^{%s}"%(x, exp)
        else:
            return "%s^%s"%(x, exp)

    cdef _dot_var(self, x, exp, do_latex):
        """
        Return a representation of '*x^exp', latexed if necessary.
        """
        if exp == 0:
            return ""
        if exp == 1:
            if do_latex:
                return " \\cdot %s"%(x)
            else:
                return "*%s"%(x)
        if do_latex:
            return " \\cdot %s^{%s}"%(x, exp)
        else:
            return "*%s^%s"%(x, exp)

    cdef _co_dot_var(self, co, x, exp, do_latex):
        """
        Return a representation of 'co*x^exp', latexed if necessary.

        co should be greater than 0
        """
        if exp == 0:
            return "%s"%co
        if exp == 1:
            if co == 1:
                return "%s"%x
            if do_latex:
                return "%s \\cdot %s"%(co, x)
            else:
                return "%s*%s"%(co, x)
        if co == 1:
            if do_latex:
                return "%s^{%s}"%(x, exp)
            else:
                return "%s^%s"%(x, exp)
        if do_latex:
            return "%s \\cdot %s^{%s}"%(co, x, exp)
        else:
            return "%s*%s^%s"%(co, x, exp)

    cdef _plus_ellipsis(self, bint do_latex):
        """
        Return a representation of '+ ...', latexed if necessary.
        """
        if do_latex:
            return " + \\cdots"
        else:
            return " + ..."

    cdef _ellipsis(self, bint do_latex):
        """
        Return a representation of '...', latexed if necessary.
        """
        if do_latex:
            return "\\cdots"
        else:
            return "..."

    cdef _truncate_list(self, L, max_terms, zero):
        """
        Take a list L of coefficients and returns a list with at most max_terms nonzero terms.

        INPUT:

        - ``L`` -- list

        - ``max_terms`` -- nonnegative integer (or -1, in which case
          no truncation occurs)

        - ``zero`` -- what should be considered zero, usually 0 or []

        OUTPUT:

        - truncated list; later terms are removed
        - boolean; whether any truncation occurred
        """
        if max_terms == -1 or len(L) == 0:
            return list(L), False
        cdef Py_ssize_t i, nonzero_index
        cdef Py_ssize_t count = 0
        ans = []
        for i, c in enumerate(L):
            if c != zero:
                count += 1
                if count > max_terms:
                    return ans[:nonzero_index+1], True
                nonzero_index = i
            ans.append(c)
        return ans, False

    cdef _print_unram_term(self, L, bint do_latex, polyname, long max_unram_terms, long expshift, bint increasing):
        """
        Return a string representation of L when considered as a polynomial,
        truncating to at most ``max_unram_terms`` nonzero terms.

        INPUT:

        - ``L`` -- a list of coefficients

        - ``do_latex`` -- whether to print latex-style

        - ``polyname`` -- the name for the variable

        - ``max_unram_terms`` -- a maximum number of terms before
              truncation occurs and an ellipsis is added.  -1
              indicates no truncation should happen.

        - ``expshift`` -- a shift for all the exponents of the variable

        - ``increasing`` -- whether to order the exponents in increasing fashion
        """
        s = ""
        cdef Py_ssize_t j, newj
        cdef long exp, count = 0
        if increasing:
            for j from 0 <= j < len(L):
                exp = j + expshift
                if L[j] != 0:
                    if max_unram_terms == 0:
                        return "(" + self._ellipsis(do_latex) + ")"
                    elif max_unram_terms == 1:
                        s = self._co_dot_var(L[j], polyname, exp, do_latex)
                        s += self._plus_ellipsis(do_latex)
                        return s
                    else:
                        count += 1
                        if count == max_unram_terms:  # this will never trigger if max_unram_terms == -1
                            newj = len(L) - 1
                            while L[newj] == 0:
                                newj -= 1
                            if newj != j:
                                exp = newj + expshift
                                s += self._plus_ellipsis(do_latex)
                            s = self._print_term_of_poly(s, L[newj], do_latex, polyname, exp)
                            return s
                        else:
                            s = self._print_term_of_poly(s, L[j], do_latex, polyname, exp)
        else:
            for j from len(L) - 1 >= j >= 0:
                exp = j + expshift
                if L[j] != 0:
                    if max_unram_terms == 0:
                        return "(" + self._ellipsis(do_latex) + ")"
                    elif max_unram_terms == 1:
                        s = self._co_dot_var(L[j], polyname, exp, do_latex)
                        s += self._plus_ellipsis(do_latex)
                        return s
                    else:
                        count += 1
                        if count == max_unram_terms:  # this will never trigger if max_unram_terms == -1
                            newj = 0
                            while L[newj] == 0:
                                newj += 1
                            if newj != j:
                                exp = newj + expshift
                                s += self._plus_ellipsis(do_latex)
                            s = self._print_term_of_poly(s, L[newj], do_latex, polyname, exp)
                            return s
                        else:
                            s = self._print_term_of_poly(s, L[j], do_latex, polyname, exp)
        return s

    cdef _terse_frac(self, a, v, u, ram_name, bint do_latex):
        """
        Return a representation of 'a=u/ram_name^v', latexed if necessary.
        """
        if do_latex:
            if v >= 0:
                arep = a._latex_()
            elif v == -1:
                arep = "\\frac{%s}{%s}"%(u, ram_name)
            else:
                arep = "\\frac{%s}{%s^{%s}}"%(u, ram_name, -v)
        else:
            if v >= 0:
                arep = str(a)
            elif v == -1:
                arep = "%s/%s"%(u, ram_name)
            else:
                arep = "%s/%s^%s"%(u, ram_name, -v)
        return arep

    cdef _print_list_as_poly(self, L, bint do_latex, polyname, long expshift, bint increasing):
        """
        Print a list ``L`` as a polynomial.

        INPUT:

        - ``L`` -- a list of coefficients

        - ``do_latex`` -- whether to print latex-style

        - ``polyname`` -- the name for the variable

        - ``expshift`` -- a shift for all the exponents of the variable

        - ``increasing`` -- whether to order the exponents in increasing fashion
        """
        s = ""
        cdef Py_ssize_t j
        cdef long exp
        if increasing:
            for j from 0 <= j < len(L):
                exp = j + expshift
                s = self._print_term_of_poly(s, L[j], do_latex, polyname, exp)
        else:
            for j from len(L) - 1 >= j >= 0:
                exp = j + expshift
                s = self._print_term_of_poly(s, L[j], do_latex, polyname, exp)
        return s

    cdef _print_term_of_poly(self, s, coeff, bint do_latex, polyname, long exp):
        """
        Appends +coeff*polyname^exp to s, latexed if necessary.
        """
        if coeff < 0:
            if len(s) > 0:
                s += " - "
            else:
                s = "-"
            coeff = -coeff
            s += self._co_dot_var(coeff, polyname, exp, do_latex)
        elif coeff > 0:
            if len(s) > 0:
                s += " + "
            s += self._co_dot_var(coeff, polyname, exp, do_latex)
        return s
