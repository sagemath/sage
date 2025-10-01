# sage_setup: distribution = sagemath-repl
"""
Check tolerance when parsing docstrings
"""

# ****************************************************************************
#       Copyright (C) 2012-2018 David Roe <roed.math@gmail.com>
#                     2012      Robert Bradshaw <robertwb@gmail.com>
#                     2012      William Stein <wstein@gmail.com>
#                     2013      R. Andrew Ohana
#                     2013      Volker Braun
#                     2013-2018 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2016-2021 Frédéric Chapoton
#                     2017-2018 Erik M. Bray
#                     2020      Marc Mezzarobba
#                     2020-2023 Matthias Koeppe
#                     2022      John H. Palmieri
#                     2022      Sébastien Labbé
#                     2023      Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import re
from sage.doctest.rif_tol import RIFtol, add_tolerance
from sage.doctest.marked_output import MarkedOutput


# Regex pattern for float without the (optional) leading sign
float_without_sign = r'((\d*\.?\d+)|(\d+\.?))([eE][+-]?\d+)?'


# Regular expression for floats
float_regex = re.compile(r'\s*([+-]?\s*' + float_without_sign + r')')


class ToleranceExceededError(BaseException):
    pass


def check_tolerance_real_domain(want: MarkedOutput, got: str) -> tuple[str, str]:
    """
    Compare want and got over real domain with tolerance

    INPUT:

    - ``want`` -- a string, what you want
    - ``got`` -- a string, what you got

    OUTPUT:

    The strings to compare, but with matching float numbers replaced by asterisk.

    EXAMPLES::

        sage: from sage.doctest.check_tolerance import check_tolerance_real_domain
        sage: from sage.doctest.marked_output import MarkedOutput
        sage: check_tolerance_real_domain(
        ....:     MarkedOutput('foo:0.2').update(abs_tol=0.3),
        ....:     'bar:0.4')
        ['foo:*', 'bar:*']
        sage: check_tolerance_real_domain(
        ....:     MarkedOutput('foo:0.2').update(abs_tol=0.3),
        ....:     'bar:0.6')
        Traceback (most recent call last):
        ...
        sage.doctest.check_tolerance.ToleranceExceededError
    """
    # First check that the number of occurrences of floats appearing match
    want_str = [g[0] for g in float_regex.findall(want)]
    got_str = [g[0] for g in float_regex.findall(got)]
    if len(want_str) != len(got_str):
        raise ToleranceExceededError()

    # Then check the numbers
    want_values = [RIFtol(g) for g in want_str]
    want_intervals = [add_tolerance(v, want) for v in want_values]
    got_values = [RIFtol(g) for g in got_str]
    # The doctest is not successful if one of the "want" and "got"
    # intervals have an empty intersection
    if not all(a.overlaps(b) for a, b in zip(want_intervals, got_values)):
        raise ToleranceExceededError()

    # Then check the part of the doctests without the numbers
    # Continue the check process with floats replaced by stars
    want = float_regex.sub('*', want)
    got = float_regex.sub('*', got)
    return [want, got]


# match 1.0 or 1.0 + I or 1.0 + 2.0*I
real_plus_optional_imag = ''.join([
    r'\s*(?P<real>[+-]?\s*',
    float_without_sign,
    r')(\s*(?P<real_imag_coeff>[+-]\s*',
    float_without_sign,
    r')\*I|\s*(?P<real_imag_unit>[+-])\s*I)?',
])


# match - 2.0*I
only_imag = ''.join([
    r'\s*(?P<only_imag>[+-]?\s*',
    float_without_sign,
    r')\*I',
])


# match I or -I (no digits), require a non-word part before and after for specificity
imaginary_unit = r'(?P<unit_imag_pre>^|\W)(?P<unit_imag>[+-]?)I(?P<unit_imag_post>$|\W)'


complex_regex = re.compile(''.join([
    '(',
    only_imag,
    '|',
    imaginary_unit,
    '|',
    real_plus_optional_imag,
    ')',
]))


def complex_match_to_real_and_imag(m: re.Match) -> tuple[str, str]:
    """
    Extract real and imaginary part from match

    INPUT:

    - ``m`` -- match from ``complex_regex``

    OUTPUT:

    Pair of real and complex parts (as string)

    EXAMPLES::

        sage: from sage.doctest.check_tolerance import complex_match_to_real_and_imag, complex_regex
        sage: complex_match_to_real_and_imag(complex_regex.match('1.0'))
        ('1.0', '0')
        sage: complex_match_to_real_and_imag(complex_regex.match('-1.0 - I'))
        ('-1.0', '-1')
        sage: complex_match_to_real_and_imag(complex_regex.match('1.0 - 3.0*I'))
        ('1.0', '- 3.0')
        sage: complex_match_to_real_and_imag(complex_regex.match('1.0*I'))
        ('0', '1.0')
        sage: complex_match_to_real_and_imag(complex_regex.match('- 2.0*I'))
        ('0', '- 2.0')
        sage: complex_match_to_real_and_imag(complex_regex.match('-I'))
        ('0', '-1')
        sage: for match in complex_regex.finditer('[1, -1, I, -1, -I]'):
        ....:     print(complex_match_to_real_and_imag(match))
        ('1', '0')
        ('-1', '0')
        ('0', '1')
        ('-1', '0')
        ('0', '-1')
        sage: for match in complex_regex.finditer('[1, -1.3, -1.5 + 0.1*I, 0.5 - 0.1*I, -1.5*I]'):
        ....:     print(complex_match_to_real_and_imag(match))
        ('1', '0')
        ('-1.3', '0')
        ('-1.5', '+ 0.1')
        ('0.5', '- 0.1')
        ('0', '-1.5')
    """
    real = m.group('real')
    if real is not None:
        real_imag_coeff = m.group('real_imag_coeff')
        real_imag_unit = m.group('real_imag_unit')
        if real_imag_coeff is not None:
            return (real, real_imag_coeff)
        elif real_imag_unit is not None:
            return (real, real_imag_unit + '1')
        else:
            return (real, '0')
    only_imag = m.group('only_imag')
    if only_imag is not None:
        return ('0', only_imag)
    unit_imag = m.group('unit_imag')
    if unit_imag is not None:
        return ('0', unit_imag + '1')
    assert False, 'unreachable'


def complex_star_repl(m: re.Match):
    """
    Replace the complex number in the match with '*'
    """
    if m.group('unit_imag') is not None:
        # preserve the matched non-word part
        return ''.join([
            (m.group('unit_imag_pre') or '').strip(),
            '*',
            (m.group('unit_imag_post') or '').strip(),
        ])
    else:
        return '*'


def check_tolerance_complex_domain(want: MarkedOutput, got: str) -> tuple[str, str]:
    """
    Compare want and got over complex domain with tolerance

    INPUT:

    - ``want`` -- a string, what you want
    - ``got`` -- a string, what you got

    OUTPUT:

    The strings to compare, but with matching complex numbers replaced by asterisk.

    EXAMPLES::

        sage: from sage.doctest.check_tolerance import check_tolerance_complex_domain
        sage: from sage.doctest.marked_output import MarkedOutput
        sage: check_tolerance_complex_domain(
        ....:     MarkedOutput('foo:[0.2 + 0.1*I]').update(abs_tol=0.3),
        ....:     'bar:[0.4]')
        ['foo:[*]', 'bar:[*]']
        sage: check_tolerance_complex_domain(
        ....:     MarkedOutput('foo:-0.5 - 0.1*I').update(abs_tol=2),
        ....:     'bar:1')
        ['foo:*', 'bar:*']
        sage: check_tolerance_complex_domain(
        ....:     MarkedOutput('foo:[1.0*I]').update(abs_tol=0.3),
        ....:     'bar:[I]')
        ['foo:[*]', 'bar:[*]']
        sage: check_tolerance_complex_domain(MarkedOutput('foo:0.2 + 0.1*I').update(abs_tol=0.3), 'bar:0.6')
        Traceback (most recent call last):
        ...
        sage.doctest.check_tolerance.ToleranceExceededError
    """
    want_str = []
    for match in complex_regex.finditer(want):
        want_str.extend(complex_match_to_real_and_imag(match))
    got_str = []
    for match in complex_regex.finditer(got):
        got_str.extend(complex_match_to_real_and_imag(match))
    if len(want_str) != len(got_str):
        raise ToleranceExceededError()

    # Then check the numbers
    want_values = [RIFtol(g) for g in want_str]
    want_intervals = [add_tolerance(v, want) for v in want_values]
    got_values = [RIFtol(g) for g in got_str]
    # The doctest is not successful if one of the "want" and "got"
    # intervals have an empty intersection
    if not all(a.overlaps(b) for a, b in zip(want_intervals, got_values)):
        raise ToleranceExceededError()

    # Then check the part of the doctests without the numbers
    # Continue the check process with floats replaced by stars
    want = complex_regex.sub(complex_star_repl, want)
    got = complex_regex.sub(complex_star_repl, got)
    return [want, got]
