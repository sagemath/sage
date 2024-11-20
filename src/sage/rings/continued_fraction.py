r"""
Continued fractions

A continued fraction is a representation of a real number in terms of a sequence
of integers denoted `[a_0; a_1, a_2, \ldots]`. The well known decimal expansion
is another way of representing a real number by a sequence of integers. The
value of a continued fraction is defined recursively as:

.. MATH::

    [a_0; a_1, a_2, \ldots] = a_0 + \frac{1}{[a_1; a_2, \ldots]} = a_0 +
    \frac{\displaystyle 1}
    {\displaystyle a_1 + \frac{\displaystyle 1}
    {\displaystyle a_2 + \frac{\displaystyle 1}
    {\ldots}}}

In this expansion, all coefficients `a_n` are integers and only the value `a_0`
may be non positive. Note that `a_0` is nothing else but the floor (this remark
provides a way to build the continued fraction expansion from a given real
number). As examples

.. MATH::

    \frac{45}{38} = 1 + \frac{\displaystyle 1}
    {\displaystyle 5 + \frac{\displaystyle 1}
    {\displaystyle 2 + \frac{\displaystyle 1}
    {\displaystyle 3}}}

.. MATH::

    \pi = 3 + \frac{\displaystyle 1}
    {\displaystyle 7 + \frac{\displaystyle 1}
    {\displaystyle 15 + \frac{\displaystyle 1}
    {\displaystyle 1 + \frac{\displaystyle 1}
    {\displaystyle 292 + \frac{\displaystyle 1}
    {\ldots}}}}}

It is quite remarkable that

- any real number admits a unique continued fraction expansion
- finite expansions correspond to rationals
- ultimately periodic expansions correspond to quadratic numbers (ie numbers of
  the form `a + b \sqrt{D}` with `a` and `b` rationals and `D` square free
  positive integer)
- two real numbers `x` and `y` have the same tail (up to a shift) in their
  continued fraction expansion if and only if there are integers `a,b,c,d` with
  `|ad - bc| = 1` and such that `y = (ax + b) / (cx + d)`.

Moreover, the rational numbers obtained by truncation of the expansion of a real
number gives its so-called best approximations. For more informations on
continued fractions, you may have a look at :wikipedia:`Continued_fraction`.

EXAMPLES:

If you want to create the continued fraction of some real number you may either
use its method continued_fraction (if it exists) or call
:func:`continued_fraction`::

    sage: (13/27).continued_fraction()
    [0; 2, 13]
    sage: 0 + 1/(2 + 1/13)
    13/27

    sage: continued_fraction(22/45)
    [0; 2, 22]
    sage: 0 + 1/(2 + 1/22)
    22/45

    sage: continued_fraction(pi)                                                        # needs sage.symbolic
    [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
    sage: continued_fraction_list(pi, nterms=5)                                         # needs sage.symbolic
    [3, 7, 15, 1, 292]

    sage: x = polygen(ZZ, 'x')
    sage: K.<cbrt5> = NumberField(x^3 - 5, embedding=1.709)                             # needs sage.rings.number_field
    sage: continued_fraction(cbrt5)                                                     # needs sage.rings.number_field
    [1; 1, 2, 2, 4, 3, 3, 1, 5, 1, 1, 4, 10, 17, 1, 14, 1, 1, 3052, 1, ...]

It is also possible to create a continued fraction from a list of partial
quotients::

    sage: continued_fraction([-3,1,2,3,4,1,2])
    [-3; 1, 2, 3, 4, 1, 2]

Even infinite::

    sage: w = words.ThueMorseWord([1,2]); w                                             # needs sage.combinat
    word: 1221211221121221211212211221211221121221...
    sage: continued_fraction(w)                                                         # needs sage.combinat
    [1; 2, 2, 1, 2, 1, 1, 2, 2, 1...]

To go back and forth between the value (as a real number) and the partial
quotients (seen as a finite or infinite list) you can use the methods
``quotients`` and ``value``::

    sage: cf = (13/27).continued_fraction()
    sage: cf.quotients()
    [0, 2, 13]
    sage: cf.value()
    13/27

    sage: cf = continued_fraction(pi)                                                   # needs sage.symbolic
    sage: cf.quotients()                                                                # needs sage.symbolic
    lazy list [3, 7, 15, ...]
    sage: cf.value()                                                                    # needs sage.symbolic
    pi

    sage: # needs sage.combinat
    sage: w = words.FibonacciWord([1,2])
    sage: cf = continued_fraction(w)
    sage: cf.quotients()
    word: 1211212112112121121211211212112112121121...
    sage: v = cf.value(); v
    1.387954587967143?
    sage: v.n(digits=100)
    1.387954587967142336919313859873185477878152452498532271894917289826418577622648932169885237034242967
    sage: v.continued_fraction()
    [1; 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2...]

Recall that quadratic numbers correspond to ultimately periodic continued
fractions. For them special methods give access to preperiod and period::

    sage: # needs sage.rings.number_field
    sage: K.<sqrt2> = QuadraticField(2)
    sage: cf = continued_fraction(sqrt2); cf
    [1; (2)*]
    sage: cf.value()
    sqrt2
    sage: cf.preperiod()
    (1,)
    sage: cf.period()
    (2,)

    sage: cf = (3*sqrt2 + 1/2).continued_fraction(); cf                                 # needs sage.rings.number_field
    [4; (1, 2, 1, 7)*]

    sage: cf = continued_fraction([(1,2,3),(1,4)]); cf
    [1; 2, 3, (1, 4)*]
    sage: cf.value()                                                                    # needs sage.rings.number_field
    -2/23*sqrt2 + 36/23

On the following we can remark how the tail may change even in the same
quadratic field::

    sage: for i in range(20): print(continued_fraction(i*sqrt2))                        # needs sage.rings.number_field
    [0]
    [1; (2)*]
    [2; (1, 4)*]
    [4; (4, 8)*]
    [5; (1, 1, 1, 10)*]
    [7; (14)*]
    ...
    [24; (24, 48)*]
    [25; (2, 5, 6, 5, 2, 50)*]
    [26; (1, 6, 1, 2, 3, 2, 26, 2, 3, 2, 1, 6, 1, 52)*]

Nevertheless, the tail is preserved under invertible integer homographies::

    sage: # needs sage.modular sage.rings.number_field
    sage: apply_homography =  lambda m,z: (m[0,0]*z + m[0,1]) / (m[1,0]*z + m[1,1])
    sage: m1 = SL2Z([60,13,83,18])
    sage: m2 = SL2Z([27,80,28,83])
    sage: a = sqrt2/3
    sage: a.continued_fraction()
    [0; 2, (8, 4)*]
    sage: b = apply_homography(m1, a)
    sage: b.continued_fraction()
    [0; 1, 2, 1, 1, 1, 1, 6, (8, 4)*]
    sage: c = apply_homography(m2, a)
    sage: c.continued_fraction()
    [0; 1, 26, 1, 2, 2, (8, 4)*]
    sage: d = apply_homography(m1**2*m2**3, a)
    sage: d.continued_fraction()
    [0; 1, 2, 1, 1, 1, 1, 5, 2, 1, 1, 1, 1, 5, 26, 1, 2, 1, 26, 1, 2, 1, 26, 1, 2, 2, (8, 4)*]

.. TODO::

    - Improve numerical approximation (the method
      :meth:`~ContinuedFraction_base._mpfr_` is quite slow compared to the
      same method for an element of a number field)

    - Make a class for generalized continued fractions of the form `a_0 +
      b_0/(a_1 + b_1/(...))` (the standard continued fractions are when all
      `b_n= 1` while the Hirzebruch-Jung continued fractions are the one for
      which `b_n = -1` for all `n`). See
      :wikipedia:`Generalized_continued_fraction`.

    - look at the function ContinuedFractionApproximationOfRoot in GAP

AUTHORS:

- Vincent Delecroix (2014): cleaning, refactorisation, documentation from the
  old implementation in ``contfrac`` (:issue:`14567`).
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2014-2020 Vincent Delecroix <20100.delecroix@gmail.com>
#       Copyright (C) 2020 Frédéric Chapoton <chapoton@math.unistra.fr>
#       Copyright (C) 2016 Moritz Firsching <moritz@math.fu-berlin.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import numbers

import sage.rings.abc

from sage.misc.lazy_import import lazy_import
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.richcmp import rich_to_bool, richcmp_method
from sage.structure.sage_object import SageObject

lazy_import('sage.combinat.words.abstract_word', 'Word_class')
lazy_import('sage.combinat.words.finite_word', 'FiniteWord_class')
lazy_import('sage.combinat.words.infinite_word', 'InfiniteWord_class')
lazy_import('sage.combinat.words.word', 'Word')


ZZ_0 = Integer(0)
ZZ_1 = Integer(1)
ZZ_m1 = Integer(-1)
ZZ_2 = Integer(2)


def last_two_convergents(x):
    """
    Given the list ``x`` that consists of numbers, return the two last
    convergents `p_{n-1}, q_{n-1}, p_n, q_n`.

    This function is principally used to compute the value of a ultimately
    periodic continued fraction.

    OUTPUT: a 4-tuple of Sage integers

    EXAMPLES::

        sage: from sage.rings.continued_fraction import last_two_convergents
        sage: last_two_convergents([])
        (0, 1, 1, 0)
        sage: last_two_convergents([0])
        (1, 0, 0, 1)
        sage: last_two_convergents([-1,1,3,2])
        (-1, 4, -2, 9)

    TESTS::

        sage: all(type(x) is Integer for x in last_two_convergents([]))
        True
    """
    p0, p1 = ZZ_0, ZZ_1
    q0, q1 = ZZ_1, ZZ_0
    for a in x:
        p0, p1 = p1, a * p1 + p0
        q0, q1 = q1, a * q1 + q0
    return p0, q0, p1, q1


def rat_interval_cf_list(r1, r2):
    r"""
    Return the common prefix of the rationals ``r1`` and ``r2`` seen as
    continued fractions.

    OUTPUT: list of Sage integers

    EXAMPLES::

        sage: from sage.rings.continued_fraction import rat_interval_cf_list
        sage: rat_interval_cf_list(257/113, 5224/2297)
        [2, 3, 1, 1, 1, 4]
        sage: for prec in range(10,54):                                                 # needs sage.rings.real_interval_field
        ....:     R = RealIntervalField(prec)
        ....:     for _ in range(100):
        ....:         x = R.random_element() * R.random_element() + R.random_element() / 100
        ....:         l = x.lower().exact_rational()
        ....:         u = x.upper().exact_rational()
        ....:         if l.floor() != u.floor():
        ....:             continue
        ....:         cf = rat_interval_cf_list(l,u)
        ....:         a = continued_fraction(cf).value()
        ....:         b = continued_fraction(cf+[1]).value()
        ....:         if a > b:
        ....:             a,b = b,a
        ....:         assert a <= l
        ....:         assert b >= u
    """
    l = []
    c1 = r1.floor()
    c2 = r2.floor()
    while c1 == c2:
        l.append(c1)
        r1 -= c1
        if not r1:
            break
        r2 -= c2
        if not r2:
            break

        r1 = ~r1
        r2 = ~r2

        c1 = r1.floor()
        c2 = r2.floor()
    return l


@richcmp_method
class ContinuedFraction_base(SageObject):
    r"""
    Base class for (standard) continued fractions.

    If you want to implement your own continued fraction, simply derived from
    this class and implement the following methods:

    - ``def quotient(self, n)``: return the ``n``-th quotient of ``self`` as a
      Sage integer

    - ``def length(self)``: the number of partial quotients of ``self`` as a
      Sage integer or ``Infinity``.

    and optionally:

    - ``def value(self)``: return the value of ``self`` (an exact real number)

    This base class will provide:

    - computation of convergents in :meth:`convergent`, :meth:`numerator` and
      :meth:`denominator`

    - comparison with other continued fractions (see :meth:`__richcmp__`)

    - elementary arithmetic function :meth:`floor`, :meth:`ceil`, :meth:`sign`

    - accurate numerical approximations :meth:`_mpfr_`

    All other methods, in particular the ones involving binary operations like
    sum or product, rely on the optional method :meth:`value` (and not on
    convergents) and may fail at execution if it is not implemented.
    """
    def __init__(self):
        r"""
        INPUT:

        - ``parent`` -- the parent of ``self``

        TESTS::

            sage: TestSuite(continued_fraction(3)).run()
        """
        self._pn = [ZZ_0, ZZ_1]
        self._qn = [ZZ_1, ZZ_0]

    def str(self, nterms=10, unicode=False, join=True):
        r"""
        Return a string representing this continued fraction.

        INPUT:

        - ``nterms`` -- the maximum number of terms to use

        - ``unicode`` -- (default: ``False``) whether to use unicode character

        - ``join`` -- (default: ``True``) if ``False`` instead of returning a
          string return a list of string, each of them representing a line

        EXAMPLES::

            sage: print(continued_fraction(pi).str())                                   # needs sage.symbolic
                                         1
            3 + ----------------------------------------------------
                                            1
                 7 + -----------------------------------------------
                                               1
                      15 + -----------------------------------------
                                                 1
                            1 + ------------------------------------
                                                     1
                                 292 + -----------------------------
                                                       1
                                        1 + ------------------------
                                                          1
                                             1 + -------------------
                                                            1
                                                  1 + --------------
                                                               1
                                                       2 + ---------
                                                            1 + ...
            sage: print(continued_fraction(pi).str(nterms=1))                           # needs sage.symbolic
            3 + ...
            sage: print(continued_fraction(pi).str(nterms=2))                           # needs sage.symbolic
                    1
            3 + ---------
                 7 + ...

            sage: print(continued_fraction(243/354).str())
                       1
            -----------------------
                         1
             1 + ------------------
                            1
                  2 + -------------
                              1
                       5 + --------
                                 1
                            3 + ---
                                 2
            sage: continued_fraction(243/354).str(join=False)
            ['           1           ',
             '-----------------------',
             '             1         ',
             ' 1 + ------------------',
             '                1      ',
             '      2 + -------------',
             '                  1    ',
             '           5 + --------',
             '                     1 ',
             '                3 + ---',
             '                     2 ']

            sage: print(continued_fraction(243/354).str(unicode=True))
                       1
            ───────────────────────
                         1
             1 + ──────────────────
                            1
                  2 + ─────────────
                              1
                       5 + ────────
                                 1
                            3 + ───
                                 2
        """
        nterms = int(nterms)
        if nterms < 0:
            raise ValueError("nterms must be positive")

        if unicode:
            import unicodedata
            frac = unicodedata.lookup('BOX DRAWINGS LIGHT HORIZONTAL')
        else:
            frac = '-'

        # get the partial quotients
        cf = [self.quotient(i) for i in range(nterms)]
        continued = self.quotient(nterms) is not Infinity
        while cf[-1] is Infinity:
            cf.pop()

        # write the lines starting from the end
        a = cf.pop()
        lines = [' {} + ... '.format(a) if continued else ' {} '.format(a)]
        w = len(lines[0])
        for a in reversed(cf):
            s = ' {} + '.format(a) if a else ' '
            w1 = len(s)
            s += frac * w
            lines.append(s)
            lines.append(' ' * w1 + '{:^{width}}'.format(1, width=w))
            w += w1

        # change the order
        lines.reverse()

        # remove extra whitespaces and equalize the width
        if len(lines) == 1:
            lines[0] = lines[0].strip()
        else:
            lines[0] = lines[0][1:]
            lines[1] = lines[1][1:]
            w = len(lines[0])
            for i, l in enumerate(lines):
                lines[i] = ' ' * (w - len(l)) + l

        return '\n'.join(lines) if join else lines

    def _ascii_art_(self):
        r"""
        EXAMPLES::

            sage: ascii_art(continued_fraction(43/30))
                      1
            1 + -------------
                        1
                 2 + --------
                           1
                      3 + ---
                           4
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self.str(unicode=False, join=False))

    def _unicode_art_(self):
        r"""
        EXAMPLES::

            sage: unicode_art(continued_fraction(43/30))
                      1
            1 + ─────────────
                        1
                 2 + ────────
                           1
                      3 + ───
                           4
        """
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(self.str(unicode=True, join=False))

    def _latex_(self, nterms=10):
        r"""
        EXAMPLES::

            sage: cf_pi = continued_fraction(pi)                                        # needs sage.symbolic
            sage: latex(cf_pi)                                                          # needs sage.symbolic
            3
            + \frac{\displaystyle 1}{\displaystyle 7
            + \frac{\displaystyle 1}{\displaystyle 15
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle 292
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle 2
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle \dots}}}}}}}}}}
            sage: print(cf_pi._latex_(nterms=3))                                        # needs sage.symbolic
            3
            + \frac{\displaystyle 1}{\displaystyle 7
            + \frac{\displaystyle 1}{\displaystyle 15
            + \frac{\displaystyle 1}{\displaystyle \dots}}}
        """
        nterms = int(nterms)
        if nterms <= 0:
            raise ValueError("nterms must be positive")
        v = [self.quotient(i) for i in range(nterms + 1)]
        s = str(v[0]) + '\n'
        for i in range(1, min(len(v) - 1, nterms)):
            s += '+ \\frac{\\displaystyle 1}{\\displaystyle %s\n' % v[i]
        s += '+ \\frac{\\displaystyle 1}{\\displaystyle \\dots'
        s += '}' * (len(v) - 1)
        return s

    def __abs__(self):
        """
        Return absolute value of ``self``.

        EXAMPLES::

            sage: a = continued_fraction(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: abs(a)
            [0; 22, 1, 7, 2]
            sage: QQ(abs(a))
            17/389
        """
        if self.quotient(0) >= 0:
            return self
        return -self

    def __richcmp__(self, other, op):
        """
        Rich comparison.

        EXAMPLES::

            sage: a = continued_fraction(-17/389)
            sage: b = continued_fraction(1/389)
            sage: c = continued_fraction([(),(1,)])     # the golden ratio
            sage: d = continued_fraction([(-1,),(1,)])
            sage: d < a and a < b and b < c
            True
            sage: d >= a
            False
            sage: d == d
            True

            sage: a == 'nothing'
            False
        """
        if not isinstance(other, ContinuedFraction_base):
            return NotImplemented

        i = 0
        while True:
            a = self.quotient(i)
            b = other.quotient(i)
            if a > b:
                return rich_to_bool(op, -1 if i % 2 else 1)
            if b > a:
                return rich_to_bool(op, 1 if i % 2 else -1)
            if a == ZZ_0 and b == ZZ_0 and i:  # rational case
                return rich_to_bool(op, 0)
            i += 1

    def _mpfr_(self, R):
        r"""
        Return a correctly-rounded numerical approximation of ``self``
        in the real mpfr ring ``R``.

        EXAMPLES::

            sage: continued_fraction(1/2).n()                                           # needs sage.rings.real_mpfr
            0.500000000000000
            sage: continued_fraction([0,4]).n()                                         # needs sage.rings.real_mpfr
            0.250000000000000
            sage: continued_fraction([12,1,3,4,2,2,3,1,2]).n(digits=4)                  # needs sage.rings.real_mpfr
            12.76

            sage: continued_fraction(12/7).n(digits=13) == (12/7).n(digits=13)          # needs sage.rings.real_mpfr
            True
            sage: continued_fraction(-14/333).n(digits=21) == (-14/333).n(digits=21)    # needs sage.rings.real_mpfr
            True

            sage: # needs sage.symbolic
            sage: a = (106*pi - 333) / (355 - 113*pi) - 292
            sage: cf = continued_fraction(a); cf
            [0; 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84, 2, 1, ...]
            sage: cf.n(digits=3)
            0.635
            sage: cf.n(digits=4)
            0.6346
            sage: cf.n(digits=5)
            0.63459
            sage: cf.n(digits=6)
            0.634591
            sage: cf.n(digits=7)
            0.6345910
            sage: cf.n(digits=8)
            0.63459101

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2, 'a', embedding=1.25)
            sage: b = 504/253*a^2 + 635/253*a + 661/253
            sage: cf = continued_fraction(b); cf
            [8; 1, 14, 1, 10, 2, 1, 4, 12, 2, 3, 2, 1, 3, 4, 1, 1, 2, 14, 3, ...]
            sage: cf.n(digits=3)
            8.94
            sage: cf.n(digits=6)
            8.93715
            sage: cf.n(digits=7)
            8.937154
            sage: cf.n(digits=8)
            8.9371541
            sage: cf.n(digits=9)
            8.93715414
            sage: cf.n(digits=10)
            8.937154138
            sage: cf.n(digits=11)
            8.9371541378

        TESTS:

        Check that the rounding works as expected (at least in the
        rational case)::

            sage: fields = []
            sage: for prec in [17, 24, 53, 128, 256]:                                   # needs sage.rings.real_mpfr
            ....:     for rnd in ['RNDN', 'RNDD', 'RNDU', 'RNDZ', 'RNDA']:
            ....:         fields.append(RealField(prec=prec, rnd=rnd))
            sage: for n in range(3000):  # long time, not tested, known bug (see :issue:`29957`)
            ....:     a = QQ.random_element(num_bound=2^(n%100))
            ....:     if a.denominator() % 8 == 0:  # not precise enough  # :issue:`29957`
            ....:         continue
            ....:     cf = continued_fraction(a)
            ....:     for R in fields:
            ....:         try:
            ....:             assert R(cf) == R(a)
            ....:         except ZeroDivisionError:  # :issue:`29957`
            ....:             pass
        """
        # 1. integer case
        if self.quotient(1) is Infinity:
            return R(self.quotient(0))

        rnd = R.rounding_mode()

        # 2. negative numbers: reduce to the positive case
        if self.quotient(0) < 0:
            sgn = -1
            self = -self
            # Adjust rounding for change in sign
            if rnd == 'RNDD':
                rnd = 'RNDA'
            elif rnd == 'RNDU':
                rnd = 'RNDZ'
        else:
            sgn = 1

        # 3. positive non integer
        if self.quotient(0) == 0:  # 0 <= self < 1
            N = R.prec() + self.quotient(1).nbits() - 1
            if self.quotient(2) is Infinity and self.quotient(1) % (1 << (self.quotient(1).nbits() - 1)) == 0:
                # if self is of the form [0; 2^N] then we need the following
                N -= 1
        else:  # self > 1
            N = R.prec() - self.quotient(0).nbits()

        # even/odd convergents are respectively below/above
        k = 0
        p_even = self.numerator(2 * k)
        p_odd = self.numerator(2 * k + 1)
        q_even = self.denominator(2 * k)
        q_odd = self.denominator(2 * k + 1)
        m_even = (p_even << N) // q_even      # floor((2^N p_even) / q_even)
        m_odd = (p_odd << N + q_odd - 1) // q_odd  # ceil((2^N p_odd) / q_odd)
        while m_odd - m_even > 1:
            k += 1
            p_even = self.numerator(2 * k)
            p_odd = self.numerator(2 * k + 1)
            q_even = self.denominator(2 * k)
            q_odd = self.denominator(2 * k + 1)
            m_even = (p_even << N) // q_even
            m_odd = ((p_odd << N) + q_odd - 1) // q_odd

        assert m_odd.nbits() == R.prec() or m_even.nbits() == R.prec()

        if m_even == m_odd:  # no need to worry (we have an exact number)
            return R(sgn * m_even) >> N

        # check ordering
        # m_even/2^N <= p_even/q_even <= self <= p_odd/q_odd <= m_odd/2^N
        assert m_odd == m_even + 1
        assert m_even / (ZZ_1 << N) <= p_even / q_even
        assert p_even / q_even <= p_odd / q_odd
        assert p_odd / q_odd <= m_odd / (ZZ_1 << N)

        if rnd == 'RNDN':  # round to the nearest
            # in order to find the nearest approximation we possibly need to
            # augment our precision on convergents.
            while True:
                assert not (p_odd << (N+1) <= (2*m_odd-1) * q_odd) or not (p_even << (N+1) >= (2*m_even+1) * q_even)
                if p_odd << (N+1) <= (2*m_odd-1) * q_odd:
                    return R(sgn * m_even) >> N
                if p_even << (N+1) >= (2*m_even+1) * q_even:
                    return R(sgn * m_odd) >> N
                k += 1
                p_even = self.numerator(2*k)
                p_odd = self.numerator(2*k+1)
                q_even = self.denominator(2*k)
                q_odd = self.denominator(2*k+1)
        elif rnd == 'RNDU' or rnd == 'RNDA':  # round up
            return R(sgn * m_odd) >> N
        elif rnd == 'RNDD' or rnd == 'RNDZ':  # round down
            return R(sgn * m_even) >> N
        else:
            raise ValueError("%s unknown rounding mode" % rnd)

    def __float__(self):
        """
        EXAMPLES::

            sage: a = continued_fraction(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: float(a)                                                              # needs sage.rings.real_mpfr
            -0.043701799485861184
            sage: float(-17/389)
            -0.043701799485861184
        """
        from sage.rings.real_mpfr import RR
        return float(self._mpfr_(RR))

    def numerator(self, n):
        """
        Return the numerator of the `n`-th partial convergent of ``self``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: c = continued_fraction(pi); c
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
            sage: c.numerator(0)
            3
            sage: c.numerator(12)
            80143857
            sage: c.numerator(152)
            3943771611212266962743738812600748213157266596588744951727393497446921245353005283
        """
        n = Integer(n)

        p = self._pn
        q = self._qn

        if n < -2:
            raise ValueError("n must be at least -2")

        for k in range(len(p), n + 3):
            x = self.quotient(k-2)
            if x is Infinity and k != 2:
                return p[-1]
            p.append(x*p[k-1] + p[k-2])
            q.append(x*q[k-1] + q[k-2])

        return p[n+2]

    p = numerator

    def denominator(self, n):
        """
        Return the denominator of the ``n``-th partial convergent of ``self``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: c = continued_fraction(pi); c
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
            sage: c.denominator(0)
            1
            sage: c.denominator(12)
            25510582
            sage: c.denominator(152)
            1255341492699841451528811722575401081588363886480089431843026103930863337221076748
        """
        self.numerator(n)   # ! silent computation of qn
        if len(self._qn) < n+3:
            return self._qn[-1]
        return self._qn[n+2]

    q = denominator

    def convergent(self, n):
        """
        Return the ``n``-th partial convergent to ``self``.

        EXAMPLES::

            sage: a = continued_fraction(pi); a                                         # needs sage.symbolic
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
            sage: a.convergent(3)                                                       # needs sage.symbolic
            355/113
            sage: a.convergent(15)                                                      # needs sage.symbolic
            411557987/131002976
        """
        return self.numerator(n) / self.denominator(n)

    def convergents(self):
        """
        Return the list of partial convergents of ``self``.

        If ``self`` is an infinite continued fraction, then the object returned
        is a :class:`~sage.misc.lazy_list.lazy_list_generic` which
        behave like an infinite list.

        EXAMPLES::

            sage: a = continued_fraction(23/157); a
            [0; 6, 1, 4, 1, 3]
            sage: a.convergents()
            [0, 1/6, 1/7, 5/34, 6/41, 23/157]

        .. TODO::

            Add an example with infinite list.
        """
        if self.length() == Infinity:
            from itertools import count

            from sage.misc.lazy_list import lazy_list
            return lazy_list(self.numerator(n) / self.denominator(n)
                             for n in count())
        return [self.numerator(n) / self.denominator(n)
                for n in range(len(self))]

    def quotients(self):
        r"""
        Return the list of partial quotients of ``self``.

        If ``self`` is an infinite continued fraction, then the object returned
        is a :class:`~sage.misc.lazy_list.lazy_list_generic` which behaves
        like an infinite list.

        EXAMPLES::

            sage: a = continued_fraction(23/157); a
            [0; 6, 1, 4, 1, 3]
            sage: a.quotients()
            [0, 6, 1, 4, 1, 3]

        .. TODO::

            Add an example with infinite list.
        """
        if self.length() == Infinity:
            from itertools import count

            from sage.misc.lazy_list import lazy_list
            return lazy_list(self.quotient(n) for n in count())
        return [self.quotient(n) for n in range(len(self))]

    def __getitem__(self, n):
        r"""
        Return the ``n``-th partial quotient of ``self`` or a continued fraction
        associated to a sublist of the partial quotients of ``self``.

        TESTS::

            sage: cf1 = continued_fraction(pi); cf1                                     # needs sage.symbolic
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
            sage: cf2 = continued_fraction(QuadraticField(2).gen()); cf2                # needs sage.rings.number_field
            [1; (2)*]
            sage: cf3 = continued_fraction(4/17); cf3
            [0; 4, 4]
            sage: cf1[3:17]                                                             # needs sage.symbolic
            [1; 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2]
            sage: cf2[:10]                                                              # needs sage.rings.number_field
            [1; 2, 2, 2, 2, 2, 2, 2, 2, 2]
            sage: cf3[1:16]
            [4; 4]

        Be careful that the truncation of an infinite continued fraction might
        be shorter by one::

            sage: len(continued_fraction(golden_ratio)[:8])                             # needs sage.symbolic
            7
        """
        if isinstance(n, slice):
            quots = self.quotients()[n]
            if n.stop is not None:
                quots = list(quots)
            return continued_fraction(quots)

        try:
            n = n.__index__()
        except (AttributeError, ValueError):
            raise ValueError("n (=%s) should be an integer" % n)
        if n < 0:
            raise ValueError("n (=%s) should be positive" % n)
        q = self.quotient(n)
        if q is Infinity:
            raise IndexError("index out of range")

        return q

    def __iter__(self):
        r"""
        Iterate over the partial quotient of ``self``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: cf = continued_fraction(pi)
            sage: i = iter(cf)
            sage: [next(i) for _ in range(10)]
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1]
            sage: [next(i) for _ in range(10)]
            [3, 1, 14, 2, 1, 1, 2, 2, 2, 2]
            sage: [next(i) for _ in range(10)]
            [1, 84, 2, 1, 1, 15, 3, 13, 1, 4]
        """
        yield self.quotient(0)
        i = 1
        while True:
            q = self.quotient(i)
            if q is Infinity:
                break
            yield q
            i += 1

    def __int__(self):
        """
        EXAMPLES::

            sage: a = continued_fraction(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: int(a)
            -1
        """
        return int(self.quotient(0))

    def sign(self):
        r"""
        Return the sign of ``self`` as an Integer.

        The sign is defined to be ``0`` if ``self`` is ``0``, ``1`` if ``self``
        is positive and ``-1`` if ``self`` is negative.

        EXAMPLES::

            sage: continued_fraction(tan(pi/7)).sign()                                  # needs sage.symbolic
            1
            sage: continued_fraction(-34/2115).sign()
            -1
            sage: continued_fraction([0]).sign()
            0
        """
        if self.quotient(0) == 0:
            if self.quotient(1) is Infinity:
                return ZZ_0
            return ZZ_1
        return self.quotient(0).sign()

    def floor(self):
        r"""
        Return the floor of ``self``.

        EXAMPLES::

            sage: cf = continued_fraction([2,1,2,3])
            sage: cf.floor()
            2
        """
        return self.quotient(0)

    def ceil(self):
        r"""
        Return the ceil of ``self``.

        EXAMPLES::

            sage: cf = continued_fraction([2,1,3,4])
            sage: cf.ceil()
            3
        """
        if self.length() == 1:
            return self.quotient(0)
        return self.quotient(0)+1

    def __bool__(self):
        """
        Return ``False`` if ``self`` is zero.

        EXAMPLES::

            sage: continued_fraction(0).is_zero()    # indirect doctest
            True
            sage: continued_fraction(1).is_zero()    # indirect doctest
            False
            sage: continued_fraction([(),(1,)]).is_zero()     # indirect doctest
            False
            sage: continued_fraction([(0,),(1,2)]).is_zero()  # indirect doctest
            False
        """
        return bool(self.quotient(0)) or self.quotient(1) is not Infinity

    def is_zero(self):
        r"""
        Test whether ``self`` is zero.

        EXAMPLES::

            sage: continued_fraction(0).is_zero()
            True
            sage: continued_fraction((0,1)).is_zero()
            False
            sage: continued_fraction(-1/2).is_zero()
            False
            sage: continued_fraction(pi).is_zero()                                      # needs sage.symbolic
            False
        """
        return self.quotient(0) == ZZ_0 and self.quotient(1) is Infinity

    def is_one(self):
        r"""
        Test whether ``self`` is one.

        EXAMPLES::

            sage: continued_fraction(1).is_one()
            True
            sage: continued_fraction(5/4).is_one()
            False
            sage: continued_fraction(0).is_one()
            False
            sage: continued_fraction(pi).is_one()                                       # needs sage.symbolic
            False
        """
        return self.quotient(0) == ZZ_1 and self.quotient(1) is Infinity

    def is_minus_one(self):
        r"""
        Test whether ``self`` is minus one.

        EXAMPLES::

            sage: continued_fraction(-1).is_minus_one()
            True
            sage: continued_fraction(1).is_minus_one()
            False
            sage: continued_fraction(0).is_minus_one()
            False
            sage: continued_fraction(-2).is_minus_one()
            False
            sage: continued_fraction([-1,1]).is_minus_one()
            False
        """
        return self.quotient(0) == ZZ_m1 and self.quotient(1) is Infinity

    def additive_order(self):
        """
        Return the additive order of this continued fraction,
        which we defined to be the additive order of its value.

        EXAMPLES::

            sage: continued_fraction(-1).additive_order()
            +Infinity
            sage: continued_fraction(0).additive_order()
            1
        """
        return Infinity if self else ZZ_1

    def multiplicative_order(self):
        """
        Return the multiplicative order of this continued fraction,
        which we defined to be the multiplicative order of its value.

        EXAMPLES::

            sage: continued_fraction(-1).multiplicative_order()
            2
            sage: continued_fraction(1).multiplicative_order()
            1
            sage: continued_fraction(pi).multiplicative_order()                         # needs sage.symbolic
            +Infinity
        """
        if self.is_zero():
            return Infinity
        if self.is_one():
            return ZZ_1
        if self.is_minus_one():
            return ZZ_2
        return Infinity

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        """
        Return a numerical approximation of this continued fraction with
        ``prec`` bits (or decimal ``digits``) of precision.

        INPUT:

        - ``prec`` -- precision in bits

        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)

        - ``algorithm`` -- ignored for continued fractions

        If neither ``prec`` nor ``digits`` is given, the default
        precision is 53 bits (roughly 16 digits).

        EXAMPLES::

            sage: w = words.FibonacciWord([1,3])                                        # needs sage.combinat
            sage: cf = continued_fraction(w); cf                                        # needs sage.combinat
            [1; 3, 1, 1, 3, 1, 3, 1, 1, 3, 1, 1, 3, 1, 3, 1, 1, 3, 1, 3...]
            sage: cf.numerical_approx(prec=53)                                          # needs sage.combinat
            1.28102513329557

        The method `n` is a shortcut to this one::

            sage: cf.n(digits=25)                                                       # needs sage.combinat
            1.281025133295569815552930
            sage: cf.n(digits=33)                                                       # needs sage.combinat
            1.28102513329556981555293038097590
        """
        from sage.arith.numerical_approx import digits_to_bits, numerical_approx_generic
        if prec is None:
            prec = digits_to_bits(digits)
        return numerical_approx_generic(self, prec)

    n = numerical_approx

    def apply_homography(self, a, b, c, d, forward_value=False):
        """
        Return the continued fraction of `(ax + b)/(cx + d)`.

        This is computed using Gosper's algorithm, see
        :mod:`~sage.rings.continued_fraction_gosper`.

        INPUT:

        - ``a``, ``b``, ``c``, ``d`` -- integers

        - ``forward_value`` -- boolean (default: ``False``); whether the
          returned continued fraction is given the symbolic value of
          `(a x + b)/(cx + d)` and not only the list of partial quotients
          obtained from Gosper's algorithm

        EXAMPLES::

            sage: (5 * 13/6 - 2) / (3 * 13/6 - 4)
            53/15
            sage: continued_fraction(13/6).apply_homography(5, -2,  3, -4).value()
            53/15

        We demonstrate now the effect of the optional argument ``forward_value``::

            sage: cf = continued_fraction(pi)                                           # needs sage.symbolic
            sage: h1 = cf.apply_homography(35, -27, 12, -5); h1                         # needs sage.symbolic
            [2; 1, 1, 6, 3, 1, 2, 1, 5, 3, 1, 1, 1, 1, 9, 12, 1, 1, 1, 3...]
            sage: h1.value()                                                            # needs sage.symbolic
            2.536941776086946?

            sage: h2 = cf.apply_homography(35, -27, 12, -5, forward_value=True); h2     # needs sage.symbolic
            [2; 1, 1, 6, 3, 1, 2, 1, 5, 3, 1, 1, 1, 1, 9, 12, 1, 1, 1, 3...]
            sage: h2.value()                                                            # needs sage.symbolic
            (35*pi - 27)/(12*pi - 5)

        TESTS::

            sage: CF = [continued_fraction(x) for x in [sqrt(2), AA(3).sqrt(),          # needs sage.rings.number_field sage.symbolic
            ....:       AA(3)**(1/3), QuadraticField(37).gen(), pi, 113/27,
            ....:       [3,1,2,2], words.FibonacciWord([1,3])]]
            sage: for _ in range(100):  # not tested, known bug (see :issue:`32086`)
            ....:     cf = choice(CF)
            ....:     forward_value = choice([True, False])
            ....:     a = ZZ.random_element(-30, 30)
            ....:     b = ZZ.random_element(-30, 30)
            ....:     c = ZZ.random_element(-30, 30)
            ....:     d = ZZ.random_element(-30, 30)
            ....:     if not c and not d:
            ....:         continue
            ....:     cf_gosper = cf.apply_homography(a, b, c, d, forward_value)
            ....:     x = cf.value()
            ....:     cf_hom = continued_fraction((a*x + b) / (c*x + d))
            ....:     assert cf_gosper[:30] == cf_hom[:30]

            sage: continued_fraction(13/25).apply_homography(0, 1, 25, -13)
            Traceback (most recent call last):
            ...
            ValueError: continued fraction can not represent infinity

            sage: continued_fraction(pi).apply_homography(0, 1, 0, 0)                   # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero

        REFERENCES:

        - [Gos1972]_
        - [Knu1998]_ Exercise 4.5.3.15
        - [LS1998]_
        """
        from .continued_fraction_gosper import gosper_iterator

        if not all(isinstance(x, numbers.Integral) for x in (a, b, c, d)):
            raise TypeError("a, b, c and d must be integral")

        a = ZZ(a)
        b = ZZ(b)
        c = ZZ(c)
        d = ZZ(d)
        if not c and not d:
            raise ZeroDivisionError("division by zero")

        _i = iter(gosper_iterator(a, b, c, d, self))
        value = None
        if isinstance(self, ContinuedFraction_periodic):
            l = list(_i)
            preperiod_length = _i.output_preperiod_length
            preperiod = l[:preperiod_length]
            period = l[preperiod_length:]
            return continued_fraction((preperiod, period), value)
        else:
            if forward_value:
                if isinstance(self, ContinuedFraction_real):
                    # self._x0 is an exact number, ie with infinite precision
                    x = self._x0
                    value = (a * x + b) / (c * x + d)
                elif isinstance(self, ContinuedFraction_infinite) and self._value is not None:
                    # if present, self._value is an exact number, ie with infinite precision
                    x = self._value
                    value = (a * x + b) / (c * x + d)

            from sage.misc.lazy_list import lazy_list
            return continued_fraction(lazy_list(_i), value)

    def __neg__(self):
        """
        Return the additive inverse of ``self``.

        EXAMPLES::

            sage: -continued_fraction(e)                                                # needs sage.symbolic
            [-3; 3, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 14...]
            sage: -continued_fraction(sqrt(7))                                          # needs sage.symbolic
            [-3; 2, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4...]
        """
        return self.apply_homography(-1, 0, 0, 1)

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: ~continued_fraction(e)                                                # needs sage.symbolic
            [0; 2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1...]
            sage: ~continued_fraction(sqrt(7))                                          # needs sage.symbolic
            [0; 2, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1...]
        """
        return self.apply_homography(0, 1, 1, 0)


class ContinuedFraction_periodic(ContinuedFraction_base):
    r"""
    Continued fraction associated with rational or quadratic number.

    A rational number has a finite continued fraction expansion (or ultimately
    0). The one of a quadratic number, ie a number of the form `a + b \sqrt{D}`
    with `a` and `b` rational, is ultimately periodic.

    .. NOTE::

        This class stores a tuple ``_x1`` for the preperiod and a tuple ``_x2``
        for the period. In the purely periodic case ``_x1`` is empty while in
        the rational case ``_x2`` is the tuple ``(0,)``.
    """
    def __init__(self, x1, x2=None, check=True):
        r"""
        INPUT:

        - ``x1`` -- tuple of integers

        - ``x2`` -- tuple of integers

        TESTS::

            sage: cf = continued_fraction((1,1,2,3,4)) # indirect doctest
            sage: loads(dumps(cf)) == cf
            True
            sage: cf = continued_fraction((1,5),(3,2)) # indirect doctest
            sage: loads(dumps(cf)) == cf
            True
        """
        ContinuedFraction_base.__init__(self)
        self._x1 = tuple(x1)
        if not x2:
            self._x2 = (Infinity,)
        else:
            self._x2 = tuple(x2)

    def period(self):
        r"""
        Return the periodic part of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K.<sqrt3> = QuadraticField(3)
            sage: cf = continued_fraction(sqrt3); cf
            [1; (1, 2)*]
            sage: cf.period()
            (1, 2)
            sage: for k in xsrange(2,40):
            ....:     if not k.is_square():
            ....:         s = QuadraticField(k).gen()
            ....:         cf = continued_fraction(s)
            ....:         print('%2d %d %s' % (k, len(cf.period()), cf))
             2 1 [1; (2)*]
             3 2 [1; (1, 2)*]
             5 1 [2; (4)*]
             6 2 [2; (2, 4)*]
             7 4 [2; (1, 1, 1, 4)*]
             8 2 [2; (1, 4)*]
            10 1 [3; (6)*]
            11 2 [3; (3, 6)*]
            12 2 [3; (2, 6)*]
            13 5 [3; (1, 1, 1, 1, 6)*]
            14 4 [3; (1, 2, 1, 6)*]
            ...
            35 2 [5; (1, 10)*]
            37 1 [6; (12)*]
            38 2 [6; (6, 12)*]
            39 2 [6; (4, 12)*]
        """
        if self._x2[0] is Infinity:
            return ()
        return self._x2

    def preperiod(self):
        r"""
        Return the preperiodic part of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K.<sqrt3> = QuadraticField(3)
            sage: cf = continued_fraction(sqrt3); cf
            [1; (1, 2)*]
            sage: cf.preperiod()
            (1,)
            sage: cf = continued_fraction(sqrt3/7); cf
            [0; 4, (24, 8)*]
            sage: cf.preperiod()
            (0, 4)
        """
        return self._x1

    def quotient(self, n):
        r"""
        Return the ``n``-th partial quotient of ``self``.

        EXAMPLES::

            sage: cf = continued_fraction([(12,5),(1,3)])
            sage: [cf.quotient(i) for i in range(10)]
            [12, 5, 1, 3, 1, 3, 1, 3, 1, 3]
        """
        n = int(n)
        if n < 0:
            raise ValueError("n (=%d) should be positive" % n)
        if n < len(self._x1):
            return self._x1[n]
        return self._x2[(n-len(self._x1)) % len(self._x2)]

    def length(self):
        r"""
        Return the number of partial quotients of ``self``.

        EXAMPLES::

            sage: continued_fraction(2/5).length()
            3
            sage: cf = continued_fraction([(0,1),(2,)]); cf
            [0; 1, (2)*]
            sage: cf.length()
            +Infinity
        """
        if len(self._x2) > 1 or self._x2[0] is not Infinity:
            return Infinity
        return Integer(len(self._x1))

    def preperiod_length(self):
        r"""
        Return the number of partial quotients of the preperiodic part of ``self``.

        EXAMPLES::

            sage: continued_fraction(2/5).preperiod_length()
            3
            sage: cf = continued_fraction([(0,1),(2,)]); cf
            [0; 1, (2)*]
            sage: cf.preperiod_length()
            2
        """
        return Integer(len(self._x1))

    def period_length(self):
        r"""
        Return the number of partial quotients of the preperiodic part of ``self``.

        EXAMPLES::

            sage: continued_fraction(2/5).period_length()
            1
            sage: cf = continued_fraction([(0,1),(2,)]); cf
            [0; 1, (2)*]
            sage: cf.period_length()
            1
        """
        return Integer(len(self._x2))

    def __richcmp__(self, other, op):
        r"""
        Rich comparison.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: a = continued_fraction([(0,),(1,2,3,1,2,3,1)]); a.n()
            0.694249167819459
            sage: b = continued_fraction([(0,),(1,2,3)]); b.n()
            0.694254176766073
            sage: c = continued_fraction([(0,1),(2,3)]); c.n()
            0.696140478029631
            sage: d = continued_fraction([(0,1,2),(3,)]); d.n()
            0.697224362268005
            sage: a < b and a < c and a < d
            True
            sage: b < c and b < d and c < d
            True
            sage: b == c
            False
            sage: c > c
            False
            sage: b >= d
            False
        """
        if isinstance(other, ContinuedFraction_periodic):
            n = max(len(self._x1) + 2 * len(self._x2),
                    len(other._x1) + 2 * len(other._x2))
            for i in range(n):
                a = self.quotient(i)
                b = other.quotient(i)
                if a > b:
                    return rich_to_bool(op, -1 if i % 2 else 1)
                if b > a:
                    return rich_to_bool(op, 1 if i % 2 else -1)
            return rich_to_bool(op, 0)

        return ContinuedFraction_base.__richcmp__(self, other, op)

    def value(self):
        r"""
        Return the value of ``self`` as a quadratic number (with square free
        discriminant).

        EXAMPLES:

        Some purely periodic examples::

            sage: cf = continued_fraction([(),(2,)]); cf
            [(2)*]
            sage: v = cf.value(); v                                                     # needs sage.rings.number_field
            sqrt2 + 1
            sage: v.continued_fraction()                                                # needs sage.rings.number_field
            [(2)*]

            sage: cf = continued_fraction([(),(1,2)]); cf
            [(1, 2)*]
            sage: v = cf.value(); v                                                     # needs sage.rings.number_field
            1/2*sqrt3 + 1/2
            sage: v.continued_fraction()                                                # needs sage.rings.number_field
            [(1, 2)*]

        The number ``sqrt3`` that appear above is actually internal to the
        continued fraction. In order to be access it from the console::

            sage: cf.value().parent().inject_variables()                                # needs sage.rings.number_field
            Defining sqrt3
            sage: sqrt3                                                                 # needs sage.rings.number_field
            sqrt3
            sage: ((sqrt3+1)/2).continued_fraction()                                    # needs sage.rings.number_field
            [(1, 2)*]

        Some ultimately periodic but non periodic examples::

            sage: cf = continued_fraction([(1,),(2,)]); cf
            [1; (2)*]
            sage: v = cf.value(); v                                                     # needs sage.rings.number_field
            sqrt2
            sage: v.continued_fraction()                                                # needs sage.rings.number_field
            [1; (2)*]

            sage: cf = continued_fraction([(1,3),(1,2)]); cf
            [1; 3, (1, 2)*]
            sage: v = cf.value(); v                                                     # needs sage.rings.number_field
            -sqrt3 + 3
            sage: v.continued_fraction()                                                # needs sage.rings.number_field
            [1; 3, (1, 2)*]

            sage: cf = continued_fraction([(-5,18), (1,3,1,5)])
            sage: cf.value().continued_fraction() == cf                                 # needs sage.rings.number_field
            True
            sage: cf = continued_fraction([(-1,),(1,)])
            sage: cf.value().continued_fraction() == cf                                 # needs sage.rings.number_field
            True

        TESTS::

            sage: a1 = ((0,1),(2,3))
            sage: a2 = ((-12,1,1),(2,3,2,4))
            sage: a3 = ((1,),(1,2))
            sage: a4 = ((-2,2),(1,124,13))
            sage: a5 = ((0,),(1,))
            sage: for a in a1,a2,a3,a4,a5:                                              # needs sage.rings.number_field
            ....:     cf = continued_fraction(a)
            ....:     assert cf.value().continued_fraction() == cf
        """
        if self._x1 and self._x1[0] < 0:
            return -(-self).value()

        if self._x2[0] is Infinity:
            return self._rational_()

        # determine the equation for the purely periodic cont. frac. determined
        # by self._x2
        p0, q0, p1, q1 = last_two_convergents(self._x2)

        # now x is one of the root of the equation
        #   q1 x^2 + (q0 - p1) x - p0 = 0
        from sage.misc.functional import squarefree_part
        from sage.rings.number_field.number_field import QuadraticField
        D = (q0-p1)**2 + 4*q1*p0
        DD = squarefree_part(D)
        Q = QuadraticField(DD, 'sqrt%d' % DD)
        x = ((p1 - q0) + (D/DD).sqrt() * Q.gen()) / (2*q1)

        # we add the preperiod
        p0, q0, p1, q1 = last_two_convergents(self._x1)
        return (p1*x + p0) / (q1*x + q0)

    def _repr_(self):
        r"""
        TESTS::

            sage: a = continued_fraction(pi.n()); a                                     # needs sage.symbolic
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
            sage: a.rename('continued fraction of pi')                                  # needs sage.symbolic
            sage: a                                                                     # needs sage.symbolic
            continued fraction of pi

            sage: continued_fraction([(0,1),(2,)])
            [0; 1, (2)*]
            sage: continued_fraction([(),(1,3)])
            [(1, 3)*]
        """
        if self._x2[0] is Infinity:  # rational case
            if len(self._x1) == 1:
                return '[%d]' % self._x1[0]
            return '[%d; ' % self._x1[0] + ', '.join(str(a) for a in self._x1[1:]) + ']'

        period = '(' + ', '.join(str(a) for a in self._x2) + ')*'
        if not self._x1:  # purely periodic case
            return '[' + period + ']'

        if len(self._x1) == 1:
            return '[%d; ' % self._x1[0] + period + ']'
        return '[%d; ' % self._x1[0] + ', '.join(str(a) for a in self._x1[1:]) + ', ' + period + ']'

    def __len__(self):
        """
        Return the number of terms in this continued fraction.

        EXAMPLES::

            sage: len(continued_fraction([1,2,3,4,5]) )
            5

            sage: len(continued_fraction(([],[1])))
            Traceback (most recent call last):
            ...
            ValueError: the length is infinite
        """
        if self._x2[0] is Infinity:  # rational case
            return len(self._x1)
        raise ValueError("the length is infinite")

    def _rational_(self):
        """
        EXAMPLES::

            sage: a = continued_fraction(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: a._rational_()
            -17/389
            sage: QQ(a)
            -17/389

            sage: a = continued_fraction([2,3,4,5,2])
            sage: QQ(a)
            344/149

            sage: a = continued_fraction(([2,3],[1]))
            sage: QQ(a)
            Traceback (most recent call last):
            ...
            ValueError: this is not a rational!
        """
        if self._x2[0] is not Infinity:
            raise ValueError("this is not a rational!")
        n = len(self)
        return self.numerator(n-1) / self.denominator(n-1)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: a = continued_fraction(-17/389)
            sage: latex(a)
            -1
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle 21
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle 7
            + \frac{\displaystyle 1}{\displaystyle 2
            }}}}}
            sage: K.<a> = QuadraticField(1234)                                          # needs sage.rings.number_field
            sage: cf = continued_fraction(a); cf                                        # needs sage.rings.number_field
            [35; (7, 1, 3, 1, 4, 4, 2, 9, 1, 1, 2, 3, 1, 1, 34, 1, 1, 3, 2, 1, 1, 9, 2, 4, 4, 1, 3, 1, 7, 70)*]
            sage: latex(cf)                                                             # needs sage.rings.number_field
            35
            + \frac{\displaystyle 1}{\displaystyle 7
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle 3
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle 4
            + \frac{\displaystyle 1}{\displaystyle 4
            + \frac{\displaystyle 1}{\displaystyle 2
            + \frac{\displaystyle 1}{\displaystyle 9
            + \frac{\displaystyle 1}{\displaystyle 1
            + \frac{\displaystyle 1}{\displaystyle \dots}}}}}}}}}}
        """
        if self._x2[0] is not Infinity:
            return super()._latex_()
        v = self._x1
        if len(v) == 0:
            return '0'
        s = str(v[0]) + '\n'
        for i in range(1, len(v)):
            s += '+ \\frac{\\displaystyle 1}{\\displaystyle %s\n' % v[i]
        s += '}'*(len(v)-1)
        return s

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: a = continued_fraction(13/25)
            sage: ~a == continued_fraction(25/13)
            True
            sage: a.value() * (~a).value()
            1

            sage: a = continued_fraction(-17/253)
            sage: ~a == continued_fraction(-253/17)
            True
            sage: a.value() * (~a).value()
            1

            sage: # needs sage.rings.number_field
            sage: K.<sqrt5> = QuadraticField(5)
            sage: a1 = (sqrt5+1)/2
            sage: c1 = a1.continued_fraction(); c1
            [(1)*]
            sage: ~c1
            [0; (1)*]
            sage: c1.value() * (~c1).value()
            1

            sage: c2 = (sqrt5/3 + 1/7).continued_fraction(); c2                         # needs sage.rings.number_field
            [0; 1, (7, 1, 17, ..., 1, 2)*]
            sage: c2.value() * (~c2).value()                                            # needs sage.rings.number_field
            1
        """
        if not self:
            raise ZeroDivisionError("rational division by zero")
        if self._x1:
            if self._x1[0] < 0:
                return -~-self
            if self._x1[0] == 0:
                return self.__class__(self._x1[1:], self._x2)
        return self.__class__((0,) + self._x1, self._x2)

    def __neg__(self):
        """
        Return the additive inverse of ``self``.

        TESTS::

            sage: quots1 = [(0,),(1,),(2,),(0,1),(1,1),(2,1),(1,2),
            ....:           (0,1,1),(1,1,1),(1,1,1,1),(1,2)]
            sage: for q in quots1:
            ....:     cf = continued_fraction(q)
            ....:     ncf = -cf
            ....:     nncf = -ncf
            ....:     assert cf == nncf
            ....:     assert ncf.value() == -cf.value()
            ....:     assert cf.length() < 2 or cf.quotients()[-1] != 1
            ....:     assert ncf.length() < 2 or ncf.quotients()[-1] != 1
            ....:     assert nncf.length() < 2 or nncf.quotients()[-1] != 1

            sage: quots2 = [((),(1,)), ((), (1,2)), ((0,),(1,)),
            ....:           ((),(2,1)), ((3,),(2,1))]
            sage: for q in quots2:                                                      # needs sage.rings.number_field
            ....:     cf = continued_fraction(q)
            ....:     ncf = -cf
            ....:     nncf = -ncf
            ....:     assert cf == nncf
            ....:     assert ncf.value() == -cf.value()
        """
        x1 = self._x1
        x2 = self._x2

        if x2[0] is Infinity:
            if len(x1) == 1:
                xx1 = (-x1[0],)

            elif x1[1] == 1:
                xx1 = (-x1[0] - 1, x1[2] + 1) + x1[3:]

            elif len(x1) == 2 and x1[1] == 2:
                xx1 = (-x1[0] - 1, ZZ_2)

            else:
                xx1 = (-x1[0] - 1, ZZ_1, x1[1] - 1) + x1[2:]

            return self.__class__(xx1, x2)

        # to make the quadratic case work, we need 3 elements in x1
        if len(x1) < 3:
            x1 += x2
            if len(x1) < 3:
                x1 += x2
                if len(x1) < 3:
                    x1 += x2

        if x1[1] == 1:
            xx1 = (-x1[0] - 1, x1[2] + 1) + x1[3:]
        else:
            xx1 = (-x1[0] - 1, ZZ_1, x1[1] - 1) + x1[2:]
        xx1, xx2 = check_and_reduce_pair(xx1, x2)
        return self.__class__(xx1, xx2)


class ContinuedFraction_real(ContinuedFraction_base):
    r"""
    Continued fraction of a real (exact) number.

    This class simply wraps a real number into an attribute (that can be
    accessed through the method :meth:`value`). The number is assumed to be
    irrational.

    EXAMPLES::

        sage: cf = continued_fraction(pi); cf                                           # needs sage.symbolic
        [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
        sage: cf.value()                                                                # needs sage.symbolic
        pi

        sage: cf = continued_fraction(e); cf                                            # needs sage.symbolic
        [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, ...]
        sage: cf.value()                                                                # needs sage.symbolic
        e
    """
    def __init__(self, x):
        r"""
        INPUT:

        - ``x`` -- the real number from which we want the continued fraction

        TESTS::

            sage: TestSuite(continued_fraction(pi)).run()                               # needs sage.symbolic
        """
        ContinuedFraction_base.__init__(self)
        self._x0 = x

        from .real_mpfi import RealIntervalField
        self._xa = RealIntervalField(53)(self._x0)
        # an approximation of the last element of the orbit under the
        # Gauss map

        self._quotients = []

    def length(self):
        r"""
        Return infinity.

        EXAMPLES::

            sage: continued_fraction(pi).length()                                       # needs sage.symbolic
            +Infinity
        """
        return Infinity

    def __len__(self):
        r"""
        TESTS::

            sage: len(continued_fraction(pi))                                           # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: the length is infinite!
        """
        raise ValueError("the length is infinite!")

    def __richcmp__(self, other, op):
        r"""
        Rich comparison.

        EXAMPLES::

            sage: continued_fraction(pi) > continued_fraction(e)                        # needs sage.symbolic
            True
            sage: continued_fraction(pi) > continued_fraction(e+4)                      # needs sage.symbolic
            False
        """
        try:
            if self.value() == other.value():
                return rich_to_bool(op, 0)
            if self.value() - other.value() > 0:
                return rich_to_bool(op, 1)
            return rich_to_bool(op, -1)
        except Exception:
            return ContinuedFraction_base.__richcmp__(self, other, op)

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: continued_fraction(pi)  # indirect doctest                            # needs sage.symbolic
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
        """
        return '[%d; ' % self.quotient(0) + ', '.join(str(self.quotient(i)) for i in range(1, 20)) + ", ...]"

    def quotient(self, n):
        r"""
        Return the ``n``-th quotient of ``self``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: cf = continued_fraction(pi)
            sage: cf.quotient(27)
            13
            sage: cf.quotient(2552)
            152
            sage: cf.quotient(10000)            # long time
            5

        The algorithm is not efficient with element of the symbolic ring and,
        if possible, one can always prefer number fields elements. The reason is
        that, given a symbolic element ``x``, there is no automatic way to
        evaluate in ``RIF`` an expression of the form ``(a*x+b)/(c*x+d)`` where
        both the numerator and the denominator are extremely small::

            sage: # needs sage.symbolic
            sage: a1 = pi
            sage: c1 = continued_fraction(a1)
            sage: p0 = c1.numerator(12); q0 = c1.denominator(12)
            sage: p1 = c1.numerator(13); q1 = c1.denominator(13)
            sage: num = (q0*a1 - p0); num.n()
            1.49011611938477e-8
            sage: den = (q1*a1 - p1); den.n()
            -2.98023223876953e-8
            sage: a1 = -num/den
            sage: RIF(a1)
            [-infinity .. +infinity]

        The same computation with an element of a number field instead of
        ``pi`` gives a very satisfactory answer::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a2> = NumberField(x^3 - 2, embedding=1.25)
            sage: c2 = continued_fraction(a2)
            sage: p0 = c2.numerator(111); q0 = c2.denominator(111)
            sage: p1 = c2.numerator(112); q1 = c2.denominator(112)
            sage: num = (q0*a2 - p0); num.n()
            -4.56719261665907e46
            sage: den = (q1*a2 - p1); den.n()
            -3.65375409332726e47
            sage: a2 = -num/den
            sage: b2 = RIF(a2); b2
            1.002685823312715?
            sage: b2.absolute_diameter()
            8.88178419700125e-16

        The consequence is that the precision needed with ``c1`` grows when we
        compute larger and larger partial quotients::

            sage: # needs sage.symbolic
            sage: c1.quotient(100)
            2
            sage: c1._xa.parent()
            Real Interval Field with 353 bits of precision
            sage: c1.quotient(200)
            3
            sage: c1._xa.parent()
            Real Interval Field with 753 bits of precision
            sage: c1.quotient(300)
            5
            sage: c1._xa.parent()
            Real Interval Field with 1053 bits of precision

            sage: # needs sage.rings.number_field
            sage: c2.quotient(200)
            6
            sage: c2._xa.parent()
            Real Interval Field with 53 bits of precision
            sage: c2.quotient(500)
            1
            sage: c2._xa.parent()
            Real Interval Field with 53 bits of precision
            sage: c2.quotient(1000)
            1
            sage: c2._xa.parent()
            Real Interval Field with 53 bits of precision
        """
        x = self._xa

        if len(self._quotients) > 1 and n >= len(self._quotients) and self._quotients[-1] == 0:
            return ZZ_0

        for k in range(len(self._quotients), n+1):
            if x.lower().is_infinity() or x.upper().is_infinity() or x.lower().floor() != x.upper().floor():

                def orbit(z):
                    return -(self.denominator(k-2)*z-self.numerator(k-2))/(self.denominator(k-1)*z-self.numerator(k-1))
                x = x.parent()(orbit(self._x0))

                # It may happen that the above line fails to give an
                # approximation with the expected number of digits (see the
                # examples). In that case, we augment the precision.
                while x.lower().is_infinity() or x.upper().is_infinity() or x.lower().floor() != x.upper().floor():
                    from .real_mpfi import RealIntervalField
                    self._prec = x.parent().prec() + 100
                    x = RealIntervalField(self._prec)(orbit(self._x0))

            self._quotients.append(x.unique_floor())
            x = (x - x.unique_floor())
            if not x:
                self._quotients.append(ZZ_0)
                return ZZ_0
            x = ~x

        self._xa = x
        return self._quotients[n]

    def value(self):
        r"""
        Return the value of ``self`` (the number from which it was built).

        EXAMPLES::

            sage: cf = continued_fraction(e)                                            # needs sage.symbolic
            sage: cf.value()                                                            # needs sage.symbolic
            e
        """
        return self._x0


class ContinuedFraction_infinite(ContinuedFraction_base):
    r"""
    A continued fraction defined by an infinite sequence of partial quotients.

    EXAMPLES::

        sage: t = continued_fraction(words.ThueMorseWord([1,2])); t                     # needs sage.combinat
        [1; 2, 2, 1, 2, 1, 1, 2, 2, 1...]
        sage: t.n(digits=100)                                                           # needs sage.combinat
        1.422388736882785488341547116024565825306879108991711829311892452916456747272565883312455412962072042

    We check that comparisons work well::

        sage: t > continued_fraction(1) and t < continued_fraction(3/2)                 # needs sage.combinat
        True
        sage: t < continued_fraction(1) or t > continued_fraction(2)                    # needs sage.combinat
        False

    Can also be called with a ``value`` option::

        sage: def f(n):
        ....:     if n % 3 == 2: return 2*(n+1)//3
        ....:     return 1
        sage: w = Word(f, alphabet=NN); w                                               # needs sage.combinat
        word: 1,1,2,1,1,4,1,1,6,1,1,8,1,1,10,1,1,12,1,1,14,1,1,16,1,1,18,1,1,20,1,1,22,1,1,24,1,1,26,1,...
        sage: cf = continued_fraction(w, value=e-1); cf                                 # needs sage.combinat sage.symbolic
        [1; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1...]

    In that case a small check is done on the input::

        sage: cf = continued_fraction(w, value=pi)                                      # needs sage.combinat sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: value evaluates to 3.141592653589794? while the continued
        fraction evaluates to 1.718281828459046? in Real Interval Field
        with 53 bits of precision.
    """
    def __init__(self, w, value=None, check=True):
        r"""
        INPUT:

        - ``parent`` -- a parent

        - ``w`` -- an infinite list

        - ``value`` -- an optional known value

        - ``check`` -- whether the constructor checks the input (default: ``True``)

        TESTS::

            sage: w = words.FibonacciWord(['a','b'])                                    # needs sage.combinat
            sage: continued_fraction(w)                                                 # needs sage.combinat
            Traceback (most recent call last):
            ...
            ValueError: the sequence must consist of integers

            sage: from itertools import count
            sage: w = Word(count(), length='infinite')                                  # needs sage.combinat
            sage: continued_fraction(w)                                                 # needs sage.combinat
            [0; 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19...]

            sage: w = Word(count(), length='unknown')                                   # needs sage.combinat
            sage: continued_fraction(w)                                                 # needs sage.combinat
            Traceback (most recent call last):
            ...
            ValueError: word with unknown length cannot be converted to
            continued fractions

            sage: continued_fraction(words.FibonacciWord([0,1]))                        # needs sage.combinat
            Traceback (most recent call last):
            ...
            ValueError: only the first partial quotient can be null

            sage: w = words.ThueMorseWord([int(1), int(2)])                             # needs sage.combinat
            sage: t = continued_fraction(w)                                             # needs sage.combinat
            sage: type(t.quotient(1))                                                   # needs sage.combinat
            <class 'sage.rings.integer.Integer'>
        """
        ContinuedFraction_base.__init__(self)
        self._w = w

        if check:
            for i in range(10):
                k = w[i]
                if not isinstance(k, Integer):
                    try:
                        k = Integer(w[i])
                    except (TypeError, ValueError):
                        raise ValueError("the sequence must consist of"
                                         " integers")
                    self.quotient = self._Integer_quotient

                if not k and i:
                    raise ValueError("only the first partial quotient can"
                                     " be null")

        if check and value is not None:
            from sage.rings.real_mpfi import RealIntervalField
            R = RealIntervalField(53)
            x = R(value)
            y = R(self)
            if x.lower() > y.lower() or x.upper() < y.upper():
                raise ValueError("value evaluates to %s while the continued fraction evaluates to %s in %s." % (x, y, R))

        self._value = value

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: w = words.FibonacciWord([3,13])                                       # needs sage.combinat
            sage: cf = continued_fraction(w); cf                                        # needs sage.combinat
            [3; 13, 3, 3, 13, 3, 13, 3, 3, 13, 3, 3, 13, 3, 13, 3, 3, 13, 3, 13...]
        """
        return "[" + str(self._w[0]) + "; " + ", ".join(map(str, self._w[1:20])) + "...]"

    def length(self):
        r"""
        Return infinity.

        EXAMPLES::

            sage: w = words.FibonacciWord([3,13])                                       # needs sage.combinat
            sage: cf = continued_fraction(w)                                            # needs sage.combinat
            sage: cf.length()                                                           # needs sage.combinat
            +Infinity
        """
        return Infinity

    def quotient(self, n):
        r"""
        Return the ``n``-th partial quotient of ``self``.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: # needs sage.combinat
            sage: w = words.FibonacciWord([1,3])
            sage: cf = continued_fraction(w)
            sage: cf.quotient(0)
            1
            sage: cf.quotient(1)
            3
            sage: cf.quotient(2)
            1
        """
        return self._w[n]

    def quotients(self):
        r"""
        Return the infinite list from which this continued fraction was built.

        EXAMPLES::

            sage: w = words.FibonacciWord([1,5])                                        # needs sage.combinat
            sage: cf = continued_fraction(w)                                            # needs sage.combinat
            sage: cf.quotients()                                                        # needs sage.combinat
            word: 1511515115115151151511511515115115151151...
        """
        return self._w

    def _Integer_quotient(self, n):
        r"""
        Return the ``n``-th partial quotient of ``self``.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: # needs sage.combinat
            sage: w = words.ThueMorseWord([int(1), int(2)])
            sage: t = continued_fraction(w)
            sage: t.quotient(0)
            1
            sage: t.quotient(1)
            2
            sage: type(t.quotient(1))      # indirect doctest
            <class 'sage.rings.integer.Integer'>
        """
        return Integer(self._w[n])

    def value(self):
        r"""
        Return the value of ``self``.

        If this value was provided on initialization, just return this value
        otherwise return an element of the real lazy field.

        EXAMPLES::

            sage: def f(n):
            ....:     if n % 3 == 2: return 2*(n+1)//3
            ....:     return 1
            sage: w = Word(f, alphabet=NN); w                                           # needs sage.combinat
            word: 1,1,2,1,1,4,1,1,6,1,1,8,1,1,10,1,1,12,1,1,14,1,1,16,1,1,18,1,1,20,1,1,22,1,1,24,1,1,26,1,...
            sage: cf = continued_fraction(w, value=e-1); cf                             # needs sage.combinat sage.symbolic
            [1; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1...]
            sage: cf.value()                                                            # needs sage.combinat sage.symbolic
            e - 1

            sage: w = words.FibonacciWord([2,5])                                        # needs sage.combinat
            sage: cf = continued_fraction(w); cf                                        # needs sage.combinat
            [2; 5, 2, 2, 5, 2, 5, 2, 2, 5, 2, 2, 5, 2, 5, 2, 2, 5, 2, 5...]
            sage: cf.value()                                                            # needs sage.combinat
            2.184951302409338?
        """
        if self._value is not None:
            return self._value
        else:
            from sage.rings.real_lazy import RLF
            if self._w[0] < 0:
                return -RLF(-self)
            return RLF(self)

    def __neg__(self):
        """
        Return the opposite of ``self``.

        EXAMPLES::

            sage: -continued_fraction(words.FibonacciWord([2,5]))                       # needs sage.combinat
            [-3; 1, 4, 2, 2, 5, 2, 5, 2, 2, 5, 2, 2, 5, 2, 5, 2, 2, 5, 2...]

            sage: from sage.misc.lazy_list import lazy_list
            sage: l = lazy_list(lambda n: (n**2)%17)
            sage: -continued_fraction(l)                                                # needs sage.combinat
            [-1; 5, 9, 16, 8, 2, 15, 13, 13, 15, 2, 8, 16, 9, 4, 1, 0, 1, 4, 9...]
        """
        from sage.combinat.words.word import Word
        _w = self._w
        if _w[1] == 1:
            _w = Word((-_w[0] - 1, _w[2] + 1)).concatenate(Word(_w[3:]))
        else:
            _w = Word((-_w[0] - 1, ZZ_1, _w[1] - 1)).concatenate(Word(_w[2:]))
        return self.__class__(_w)


def check_and_reduce_pair(x1, x2=None):
    r"""
    There are often two ways to represent a given continued fraction. This
    function makes it canonical.

    In the very special case of the number `0` we return the pair
    ``((0,),(0,))``.

    TESTS::

        sage: from sage.rings.continued_fraction import check_and_reduce_pair
        sage: check_and_reduce_pair([])
        Traceback (most recent call last):
        ...
        ValueError: continued fraction can not represent infinity
        sage: check_and_reduce_pair([-1,1])
        ((0,), (+Infinity,))
        sage: check_and_reduce_pair([1,1,1])
        ((1, 2), (+Infinity,))
        sage: check_and_reduce_pair([1,3],[2,3])
        ((1,), (3, 2))
        sage: check_and_reduce_pair([1,2,3],[2,3,2,3,2,3])
        ((1,), (2, 3))
        sage: check_and_reduce_pair([1,2],[])
        ((1, 2), (+Infinity,))
    """
    y1 = [Integer(x) for x in x1]

    if x2 is None or not x2 or x2[0] is Infinity:
        y2 = [Infinity]
        if not y1:
            raise ValueError("continued fraction can not represent infinity")
        elif len(y1) > 1 and y1[-1] == 1:
            y1.pop()
            y1[-1] += 1

    else:
        y2 = [Integer(x) for x in x2]
        if any(b <= ZZ_0 for b in y2):
            raise ValueError("the elements of the period cannot be negative")

    # add possibly some element of x1 into the period
    while y1 and y1[-1] == y2[-1]:
        y1.pop()
        y2.insert(0, y2.pop())

    # some special cases to treat
    if len(y2) == 1 and y2[0] == 0:
        if not y1:
            y1 = [ZZ_0]
        elif len(y1) > 1 and y1[-1] == 1:
            y1.pop()
            y1[-1] += 1

    # check that y2 is not a pure power (in a very naive way!!)
    n2 = len(y2)
    for i in range(1, (n2 + 2) // 2):
        if n2 % i == 0 and y2[:-i] == y2[i:]:
            y2 = y2[:i]
            break

    # check that at the end y1 has no zeros in it
    for i in range(1, len(y1)):
        if y1[i] <= 0:
            raise ValueError("all quotient except the first must be positive")

    return tuple(y1), tuple(y2)


def continued_fraction_list(x, type='std', partial_convergents=False,
                            bits=None, nterms=None):
    r"""
    Return the (finite) continued fraction of ``x`` as a list.

    The continued fraction expansion of ``x`` are the coefficients `a_i` in

    .. MATH::

        x = a_0 + 1/(a_1 + 1/(...))

    with `a_0` integer and `a_1`, `...` positive integers. The Hirzebruch-Jung
    continued fraction is the one for which the `+` signs are replaced with `-`
    signs

    .. MATH::

        x = a_0 - 1/(a_1 - 1/(...))

    .. SEEALSO::

        :func:`continued_fraction`

    INPUT:

    - ``x`` -- exact rational or floating-point number; the number to
      compute the continued fraction of

    - ``type`` -- either ``'std'`` (default) for standard continued fractions or
      ``'hj'`` for Hirzebruch-Jung ones

    - ``partial_convergents`` -- boolean; whether to return the
      partial convergents

    - ``bits`` -- an optional integer that specify a precision for the real
      interval field that is used internally

    - ``nterms`` -- integer; the upper bound on the number of terms in
      the continued fraction expansion to return

    OUTPUT:

    A lits of integers, the coefficients in the continued fraction expansion of
    ``x``. If ``partial_convergents`` is set to ``True``, then return a pair
    containing the coefficient list and the partial convergents list is
    returned.

     EXAMPLES::

        sage: continued_fraction_list(45/19)
        [2, 2, 1, 2, 2]
        sage: 2 + 1/(2 + 1/(1 + 1/(2 + 1/2)))
        45/19

        sage: continued_fraction_list(45/19, type='hj')
        [3, 2, 3, 2, 3]
        sage: 3 - 1/(2 - 1/(3 - 1/(2 - 1/3)))
        45/19

    Specifying ``bits`` or ``nterms`` modify the length of the output::

        sage: # needs sage.symbolic
        sage: continued_fraction_list(e, bits=20)
        [2, 1, 2, 1, 1, 4, 2]
        sage: continued_fraction_list(sqrt(2) + sqrt(3), bits=30)
        [3, 6, 1, 5, 7, 2]
        sage: continued_fraction_list(pi, bits=53)
        [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
        sage: continued_fraction_list(log(3/2), nterms=15)
        [0, 2, 2, 6, 1, 11, 2, 1, 2, 2, 1, 4, 3, 1, 1]
        sage: continued_fraction_list(tan(sqrt(pi)), nterms=20)
        [-5, 9, 4, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1, 1, 1, 2, 4, 3, 1, 63]

    When the continued fraction is infinite (ie ``x`` is an irrational number)
    and the parameters ``bits`` and ``nterms`` are not specified then a warning
    is raised::

        sage: continued_fraction_list(sqrt(2))                                          # needs sage.symbolic
        doctest:...: UserWarning: the continued fraction of sqrt(2) seems infinite,
        return only the first 20 terms
        [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: continued_fraction_list(sqrt(4/19))                                       # needs sage.symbolic
        doctest:...: UserWarning: the continued fraction of 2*sqrt(1/19) seems infinite,
        return only the first 20 terms
        [0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 16]

    An examples with the list of partial convergents::

        sage: continued_fraction_list(RR(pi), partial_convergents=True)                 # needs sage.symbolic
        ([3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3],
         [(3, 1),
          (22, 7),
          (333, 106),
          (355, 113),
          (103993, 33102),
          (104348, 33215),
          (208341, 66317),
          (312689, 99532),
          (833719, 265381),
          (1146408, 364913),
          (4272943, 1360120),
          (5419351, 1725033),
          (80143857, 25510582),
          (245850922, 78256779)])

    TESTS::

        sage: continued_fraction_list(1 + 10^-10, nterms=3)
        [1, 10000000000]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, nterms=3)                    # needs sage.symbolic
        [1, 100000000000000000000, 2688]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, nterms=5)                    # needs sage.symbolic
        [1, 100000000000000000000, 2688, 8, 1]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, nterms=5)                    # needs sage.symbolic
        [1, 100000000000000000000, 2688, 8, 1]

    Fixed :issue:`18901`::

        sage: a = 1.575709393346379
        sage: type(a)                                                                   # needs sage.rings.real_mpfr
        <class 'sage.rings.real_mpfr.RealLiteral'>
        sage: continued_fraction_list(a)
        [1, 1, 1, 2, 1, 4, 18, 1, 5, 2, 25037802, 7, 1, 3, 1, 28, 1, 8, 2]

    Check that this works for arb elements (:issue:`20069`)::

        sage: continued_fraction(RBF(e))                                                # needs sage.symbolic
        [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12]
    """
    from .rational_field import QQ

    try:
        return x.continued_fraction_list(type=type)
    except AttributeError:
        pass

    if bits is not None:
        from .real_mpfi import RealIntervalField
        x = RealIntervalField(bits)(x)

    if type == "hj":
        l = QQ(x).continued_fraction_list("hj")
        # The C-code in sage.rings.rational is much more faster than the pure
        # Python below
        # v = []
        # while True:
        #    div, mod = divmod(x.numerator(), x.denominator())
        #    if mod == 0:
        #        v.append(div)
        #        break
        #    v.append(div+1)
        #    if nterms is not None and len(v) >= nterms:
        #        break
        #    x = 1/(div+1-x)
        # return v
        if nterms is None:
            return l
        return l[:nterms]

    if type != "std":
        raise ValueError("type must be either \"std\" or \"hj\"")

    cf = None

    try:
        from sage.rings.real_mpfr import RealLiteral
    except ImportError:
        RealLiteral = ()
    if isinstance(x, RealLiteral):
        from sage.rings.real_mpfi import RealIntervalField
        x = RealIntervalField(x.prec())(x)
    if isinstance(x.parent(), (sage.rings.abc.RealIntervalField, sage.rings.abc.RealBallField)):
        cf = continued_fraction(rat_interval_cf_list(
            x.lower().exact_rational(),
            x.upper().exact_rational()))

    if cf is None:
        try:
            cf = continued_fraction(x)
        except ValueError:
            pass

    if cf is None:
        raise ValueError("does not know how to compute the continued fraction of %s" % x)

    if nterms:
        limit = min(cf.length(), nterms)
    elif cf.length() != Infinity:
        limit = cf.length()
    else:
        import warnings
        warnings.warn("the continued fraction of %s seems infinite, return only the first 20 terms" % x)
        limit = 20
    if partial_convergents:
        return ([cf.quotient(i) for i in range(limit)],
                [(cf.numerator(i), cf.denominator(i)) for i in range(limit)])
    return [cf.quotient(i) for i in range(limit)]


def continued_fraction(x, value=None):
    r"""
    Return the continued fraction of `x`.

    INPUT:

    - ``x`` -- a number or a list of partial quotients (for finite
      development) or two list of partial quotients (preperiod and period
      for ultimately periodic development)

    EXAMPLES:

    A finite continued fraction may be initialized by a number or by its list of
    partial quotients::

        sage: continued_fraction(12/571)
        [0; 47, 1, 1, 2, 2]
        sage: continued_fraction([3,2,1,4])
        [3; 2, 1, 4]

    It can be called with elements defined from symbolic values, in which case
    the partial quotients are evaluated in a lazy way::

        sage: c = continued_fraction(golden_ratio); c                                   # needs sage.symbolic
        [1; 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...]
        sage: c.convergent(12)                                                          # needs sage.symbolic
        377/233
        sage: fibonacci(14)/fibonacci(13)                                               # needs sage.libs.pari
        377/233

        sage: # needs sage.symbolic
        sage: continued_fraction(pi)
        [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
        sage: c = continued_fraction(pi); c
        [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
        sage: a = c.convergent(3); a
        355/113
        sage: a.n()
        3.14159292035398
        sage: pi.n()
        3.14159265358979

        sage: # needs sage.symbolic
        sage: continued_fraction(sqrt(2))
        [1; 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, ...]
        sage: continued_fraction(tan(1))
        [1; 1, 1, 3, 1, 5, 1, 7, 1, 9, 1, 11, 1, 13, 1, 15, 1, 17, 1, 19, ...]
        sage: continued_fraction(tanh(1))
        [0; 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, ...]
        sage: continued_fraction(e)
        [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, ...]

    If you want to play with quadratic numbers (such as ``golden_ratio`` and
    ``sqrt(2)`` above), it is much more convenient to use number fields as
    follows since preperiods and periods are computed::

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: K.<sqrt5> = NumberField(x^2 - 5, embedding=2.23)
        sage: my_golden_ratio = (1 + sqrt5)/2
        sage: cf = continued_fraction((1+sqrt5)/2); cf
        [(1)*]
        sage: cf.convergent(12)
        377/233
        sage: cf.period()
        (1,)
        sage: cf = continued_fraction(2/3+sqrt5/5); cf
        [1; 8, (1, 3, 1, 1, 3, 9)*]
        sage: cf.preperiod()
        (1, 8)
        sage: cf.period()
        (1, 3, 1, 1, 3, 9)

        sage: # needs sage.rings.number_field
        sage: L.<sqrt2> = NumberField(x^2 - 2, embedding=1.41)
        sage: cf = continued_fraction(sqrt2); cf
        [1; (2)*]
        sage: cf.period()
        (2,)
        sage: cf = continued_fraction(sqrt2/3); cf
        [0; 2, (8, 4)*]
        sage: cf.period()
        (8, 4)

    It is also possible to go the other way around, build a ultimately periodic
    continued fraction from its preperiod and its period and get its value
    back::

        sage: cf = continued_fraction([(1,1), (2,8)]); cf
        [1; 1, (2, 8)*]
        sage: cf.value()                                                                # needs sage.rings.number_field
        2/11*sqrt5 + 14/11

    It is possible to deal with higher degree number fields but in that case the
    continued fraction expansion is known to be aperiodic::

        sage: K.<a> = NumberField(x^3 - 2, embedding=1.25)                              # needs sage.rings.number_field
        sage: cf = continued_fraction(a); cf                                            # needs sage.rings.number_field
        [1; 3, 1, 5, 1, 1, 4, 1, 1, 8, 1, 14, 1, 10, 2, 1, 4, 12, 2, 3, ...]

    Note that initial rounding can result in incorrect trailing partial
    quotients::

        sage: continued_fraction(RealField(39)(e))                                      # needs sage.symbolic
        [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 2]

    Note the value returned for floating point number is the continued fraction
    associated to the rational number you obtain with a conversion::

        sage: for _ in range(10):
        ....:     x = RR.random_element()
        ....:     cff = continued_fraction(x)
        ....:     cfe = QQ(x).continued_fraction()
        ....:     assert cff == cfe, "%s %s %s"%(x,cff,cfe)

    TESTS:

    Fixed :issue:`18901`. For RealLiteral, continued_fraction calls
    continued_fraction_list::

        sage: continued_fraction(1.575709393346379)
        [1; 1, 1, 2, 1, 4, 18, 1, 5, 2, 25037802, 7, 1, 3, 1, 28, 1, 8, 2]

    Constants in symbolic subrings work like constants in ``SR``::

        sage: SCR = SR.subring(no_variables=True)                                       # needs sage.symbolic
        sage: continued_fraction(SCR(pi))                                               # needs sage.symbolic
        [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
    """

    if isinstance(x, ContinuedFraction_base):
        return x

    try:
        return x.continued_fraction()
    except AttributeError:
        pass

    # input for finite or ultimately periodic partial quotient expansion
    if isinstance(x, FiniteWord_class):
        x = list(x)

    if isinstance(x, (list, tuple)):
        if len(x) == 2 and isinstance(x[0], (list, tuple)) and isinstance(x[1], (list, tuple)):
            x1 = tuple(Integer(a) for a in x[0])
            x2 = tuple(Integer(a) for a in x[1])
            x1, x2 = check_and_reduce_pair(x1, x2)
        else:
            x1, x2 = check_and_reduce_pair(x)
        return ContinuedFraction_periodic(x1, x2)

    # input for infinite partial quotient expansion
    from sage.misc.lazy_list import lazy_list_generic
    if isinstance(x, (lazy_list_generic, InfiniteWord_class)):
        return ContinuedFraction_infinite(x, value)

    if isinstance(x, Word_class):
        raise ValueError("word with unknown length cannot be converted "
                         "to continued fractions")

    # input for numbers
    # TODO: the approach used below might be not what the user expects as we
    # have currently in sage (version 6.8)
    #
    #     sage: RR.random_element() in QQ
    #     True
    #
    # But, be careful with real literals
    #
    #     sage: a = 1.575709393346379
    #     sage: a in QQ
    #     False
    from .rational_field import QQ
    if x in QQ:
        return QQ(x).continued_fraction()

    is_real = False
    try:
        is_real = x.is_real()
    except AttributeError:
        pass

    if is_real is False:
        from .real_mpfi import RealIntervalField
        # we cannot rely on the answer of .is_real() for elements of the
        # symbolic ring. The thing below is a dirty temporary hack.
        RIF = RealIntervalField(53)
        try:
            RIF(x)
            is_real = True
        except (AttributeError, ValueError):
            pass

    if is_real is False:
        raise ValueError("the number %s does not seem to be a real number" % x)

    if x.parent().is_exact():
        return ContinuedFraction_real(x)

    # We treat the Symbolic Ring and its subrings separately.  They hold all constants and
    # are not exact.
    if isinstance(x.parent(), sage.rings.abc.SymbolicRing):
        return ContinuedFraction_real(x)

    return continued_fraction(continued_fraction_list(x))


def convergents(x):
    r"""
    Return the (partial) convergents of the number ``x``.

    EXAMPLES::

        sage: from sage.rings.continued_fraction import convergents
        sage: convergents(143/255)
        [0, 1, 1/2, 4/7, 5/9, 9/16, 14/25, 23/41, 60/107, 143/255]
    """
    return continued_fraction(x).convergents()
