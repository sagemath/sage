# sage_setup: distribution = sagemath-repl
"""
Helper for attaching metadata to doctest output.
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

from typing import TypedDict

from typing_extensions import NotRequired, Unpack


class MarkedOutputType(TypedDict):
    random: NotRequired[bool]
    rel_tol: NotRequired[float]
    abs_tol: NotRequired[float]
    tol: NotRequired[float]
    bitness_32: NotRequired[str]
    bitness_64: NotRequired[str]


class MarkedOutput(str):
    """
    A subclass of string with context for whether another string
    matches it.

    EXAMPLES::

        sage: from sage.doctest.marked_output import MarkedOutput
        sage: s = MarkedOutput("abc")
        sage: s.rel_tol
        0
        sage: s.update(rel_tol = .05)
        'abc'
        sage: s.rel_tol
        0.0500000000000000

        sage: MarkedOutput("56 µs")
        '56 \xb5s'
    """

    random = False
    rel_tol = 0
    abs_tol = 0
    tol = 0
    bitness_32 = ""
    bitness_64 = ""

    def update(self, **kwargs: Unpack[MarkedOutputType]):
        """
        EXAMPLES::

            sage: from sage.doctest.marked_output import MarkedOutput
            sage: s = MarkedOutput("0.0007401")
            sage: s.update(abs_tol = .0000001)
            '0.0007401'
            sage: s.rel_tol
            0
            sage: s.abs_tol
            1.00000000000000e-7
        """
        self.__dict__.update(kwargs)
        return self

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: from sage.doctest.marked_output import MarkedOutput
            sage: s = MarkedOutput("0.0007401")
            sage: s.update(abs_tol = .0000001)
            '0.0007401'
            sage: t = loads(dumps(s)) # indirect doctest
            sage: t == s
            True
            sage: t.abs_tol
            1.00000000000000e-7
        """
        return make_marked_output, (str(self), self.__dict__)

    def __eq__(self, value: object) -> bool:
        if isinstance(value, MarkedOutput):
            return str(self) == str(value) and self.__dict__ == value.__dict__
        return False


def make_marked_output(s, D):
    """
    Auxiliary function for pickling.

    EXAMPLES::

        sage: from sage.doctest.marked_output import make_marked_output
        sage: s = make_marked_output("0.0007401", {'abs_tol':.0000001})
        sage: s
        '0.0007401'
        sage: s.abs_tol
        1.00000000000000e-7
    """
    ans = MarkedOutput(s)
    ans.__dict__.update(D)
    return ans
