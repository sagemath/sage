# -*- encoding: utf-8 -*-
"""
Parsing docstrings

This module contains functions and classes that parse docstrings.

AUTHORS:

- David Roe (2012-03-27) -- initial version, based on Robert Bradshaw's code.

- Jeroen Demeyer(2014-08-28) -- much improved handling of tolerances
  using interval arithmetic (:trac:`16889`).
"""

# ****************************************************************************
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import re
import doctest
from collections import defaultdict
from sage.repl.preparse import preparse, strip_string_literals
from functools import reduce


from .external import available_software

float_regex = re.compile(r'\s*([+-]?\s*((\d*\.?\d+)|(\d+\.?))([eE][+-]?\d+)?)')
optional_regex = re.compile(r'(arb216|arb218|py2|py3|long time|not implemented|not tested|known bug)|([^ a-z]\s*optional\s*[:-]*((\s|\w|[.])*))')
# Version 4.65 of glpk prints the warning "Long-step dual simplex will
# be used" frequently. When Sage uses a system installation of glpk
# which has not been patched, we need to ignore that message.
# See :trac:`29317`.
glpk_simplex_warning_regex = re.compile(r'(Long-step dual simplex will be used)')
# :trac:`31204` -- suppress warning about ld and OS version for dylib files.
ld_warning_regex = re.compile(r'.*dylib.*was built for newer macOS version.*than being linked.*')
find_sage_prompt = re.compile(r"^(\s*)sage: ", re.M)
find_sage_continuation = re.compile(r"^(\s*)\.\.\.\.:", re.M)
find_python_continuation = re.compile(r"^(\s*)\.\.\.([^\.])", re.M)
python_prompt = re.compile(r"^(\s*)>>>", re.M)
# The following are used to allow ... at the beginning of output
ellipsis_tag = "<TEMP_ELLIPSIS_TAG>"
continuation_tag = "<TEMP_CONTINUATION_TAG>"
random_marker = re.compile('.*random', re.I)
tolerance_pattern = re.compile(r'\b((?:abs(?:olute)?)|(?:rel(?:ative)?))? *?tol(?:erance)?\b( +[0-9.e+-]+)?')
backslash_replacer = re.compile(r"""(\s*)sage:(.*)\\\ *
\ *(((\.){4}:)|((\.){3}))?\ *""")

_RIFtol = None


def RIFtol(*args):
    """
    Create an element of the real interval field used for doctest tolerances.

    It allows large numbers like 1e1000, it parses strings with spaces
    like ``RIF(" - 1 ")`` out of the box and it carries a lot of
    precision. The latter is useful for testing libraries using
    arbitrary precision but not guaranteed rounding such as PARI. We use
    1044 bits of precision, which should be good to deal with tolerances
    on numbers computed with 1024 bits of precision.

    The interval approach also means that we do not need to worry about
    rounding errors and it is also very natural to see a number with
    tolerance as an interval.

    EXAMPLES::

        sage: from sage.doctest.parsing import RIFtol
        sage: RIFtol(-1, 1)
        0.?
        sage: RIFtol(" - 1 ")
        -1
        sage: RIFtol("1e1000")
        1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000?e1000
    """
    global _RIFtol
    if _RIFtol is None:
        try:
            # We need to import from sage.all to avoid circular imports.
            from sage.all import RealIntervalField
        except ImportError:
            from warnings import warn
            warn("RealIntervalField not available, ignoring all tolerance specifications in doctests")

            def fake_RIFtol(*args):
                return 0
            _RIFtol = fake_RIFtol
        else:
            _RIFtol = RealIntervalField(1044)
    return _RIFtol(*args)

# This is the correct pattern to match ISO/IEC 6429 ANSI escape sequences:
#
ansi_escape_sequence = re.compile(r'(\x1b[@-Z\\-~]|\x1b\[.*?[@-~]|\x9b.*?[@-~])')


_long_repr_re = re.compile(r'([+-]?[0-9]+)[lL]')
def normalize_long_repr(s):
    """
    Simple conversion from Python 2 representation of ``long`` ints (that
    is, integers with the ``L``) suffix, to the Python 3 representation
    (same number, without the suffix, since Python 3 doesn't have a
    distinct ``long`` type).

    Note: This just uses a simple regular expression that can't distinguish
    representations of long objects from strings containing a long repr.

    EXAMPLES::

        sage: from sage.doctest.parsing import normalize_long_repr
        sage: normalize_long_repr('10L')
        '10'
        sage: normalize_long_repr('[10L, -10L, +10L, "ALL"]')
        '[10, -10, +10, "ALL"]'
    """
    return _long_repr_re.sub(lambda m: m.group(1), s)


# Collection of fixups applied in the SageOutputChecker.  Each element in this
# this list a pair of functions applied to the actual test output ('g' for
# "got") and the expected test output ('w' for "wanted").  The first function
# should be a simple fast test on the expected and/or actual output to
# determine if a fixup should be applied.  The second function is the actual
# fixup, which is applied if the test function passes.  In most fixups only one
# of the expected or received outputs are normalized, depending on the
# application.
# For example, on Python 3 we strip all u prefixes from unicode strings in the
# expected output, because we never expect to see those on Python 3.
_repr_fixups = [
    (lambda g, w: 'L' in w or 'l' in w,
     lambda g, w: (g, normalize_long_repr(w)))
]


def parse_optional_tags(string):
    """
    Return a set consisting of the optional tags from the following
    set that occur in a comment on the first line of the input string.

    - 'long time'
    - 'not implemented'
    - 'not tested'
    - 'known bug'
    - 'py2'
    - 'py3'
    - 'arb216'
    - 'arb218'
    - 'optional: PKG_NAME' -- the set will just contain 'PKG_NAME'

    EXAMPLES::

        sage: from sage.doctest.parsing import parse_optional_tags
        sage: parse_optional_tags("sage: magma('2 + 2')# optional: magma")
        {'magma'}
        sage: parse_optional_tags("sage: #optional -- mypkg")
        {'mypkg'}
        sage: parse_optional_tags("sage: print(1)  # parentheses are optional here")
        set()
        sage: parse_optional_tags("sage: print(1)  # optional")
        {''}
        sage: sorted(list(parse_optional_tags("sage: #optional -- foo bar, baz")))
        ['bar', 'foo']
        sage: parse_optional_tags("sage: #optional -- foo.bar, baz")
        {'foo.bar'}
        sage: sorted(list(parse_optional_tags("    sage: factor(10^(10^10) + 1) # LoNg TiME, NoT TeSTED; OptioNAL -- P4cka9e")))
        ['long time', 'not tested', 'p4cka9e']
        sage: parse_optional_tags("    sage: raise RuntimeError # known bug")
        {'bug'}
        sage: sorted(list(parse_optional_tags("    sage: determine_meaning_of_life() # long time, not implemented")))
        ['long time', 'not implemented']

    We don't parse inside strings::

        sage: parse_optional_tags("    sage: print('  # long time')")
        set()
        sage: parse_optional_tags("    sage: print('  # long time')  # not tested")
        {'not tested'}

    UTF-8 works::

         sage: parse_optional_tags("'ěščřžýáíéďĎ'")
         set()
    """
    safe, literals, state = strip_string_literals(string)
    first_line = safe.split('\n', 1)[0]
    if '#' not in first_line:
        return set()
    comment = first_line[first_line.find('#')+1:]
    comment = comment[comment.index('(')+1 : comment.rindex(')')]
    # strip_string_literals replaces comments
    comment = "#" + (literals[comment]).lower()

    tags = []
    for m in optional_regex.finditer(comment):
        cmd = m.group(1)
        if cmd == 'known bug':
            tags.append('bug') # so that such tests will be run by sage -t ... -only-optional=bug
        elif cmd:
            tags.append(cmd)
        else:
            tags.extend(m.group(3).split() or [""])
    return set(tags)


def parse_tolerance(source, want):
    """
    Return a version of ``want`` marked up with the tolerance tags
    specified in ``source``.

    INPUT:

    - ``source`` -- a string, the source of a doctest
    - ``want`` -- a string, the desired output of the doctest

    OUTPUT:

    - ``want`` if there are no tolerance tags specified; a
      :class:`MarkedOutput` version otherwise.

    EXAMPLES::

        sage: from sage.doctest.parsing import parse_tolerance
        sage: marked = parse_tolerance("sage: s.update(abs_tol = .0000001)", "")
        sage: type(marked)
        <class 'str'>
        sage: marked = parse_tolerance("sage: s.update(tol = 0.1); s.rel_tol # abs tol     0.01 ", "")
        sage: marked.tol
        0
        sage: marked.rel_tol
        0
        sage: marked.abs_tol
        0.010000000000000000000...?
    """
    safe, literals, state = strip_string_literals(source)
    first_line = safe.split('\n', 1)[0]
    if '#' not in first_line:
        return want
    comment = first_line[first_line.find('#')+1:]
    comment = comment[comment.index('(')+1 : comment.rindex(')')]
    # strip_string_literals replaces comments
    comment = literals[comment]
    if random_marker.search(comment):
        want = MarkedOutput(want).update(random=True)
    else:
        m = tolerance_pattern.search(comment)
        if m:
            rel_or_abs, epsilon = m.groups()
            if epsilon is None:
                epsilon = RIFtol("1e-15")
            else:
                epsilon = RIFtol(epsilon)
            if rel_or_abs is None:
                want = MarkedOutput(want).update(tol=epsilon)
            elif rel_or_abs.startswith('rel'):
                want = MarkedOutput(want).update(rel_tol=epsilon)
            elif rel_or_abs.startswith('abs'):
                want = MarkedOutput(want).update(abs_tol=epsilon)
            else:
                raise RuntimeError
    return want


def pre_hash(s):
    """
    Prepends a string with its length.

    EXAMPLES::

        sage: from sage.doctest.parsing import pre_hash
        sage: pre_hash("abc")
        '3:abc'
    """
    return "%s:%s" % (len(s), s)


def get_source(example):
    """
    Return the source with the leading 'sage: ' stripped off.

    EXAMPLES::

        sage: from sage.doctest.parsing import get_source
        sage: from sage.doctest.sources import DictAsObject
        sage: example = DictAsObject({})
        sage: example.sage_source = "2 + 2"
        sage: example.source = "sage: 2 + 2"
        sage: get_source(example)
        '2 + 2'
        sage: example = DictAsObject({})
        sage: example.source = "3 + 3"
        sage: get_source(example)
        '3 + 3'
    """
    return getattr(example, 'sage_source', example.source)

def reduce_hex(fingerprints):
    """
    Return a symmetric function of the arguments as hex strings.

    The arguments should be 32 character strings consisting of hex
    digits: 0-9 and a-f.

    EXAMPLES::

        sage: from sage.doctest.parsing import reduce_hex
        sage: reduce_hex(["abc", "12399aedf"])
        '0000000000000000000000012399a463'
        sage: reduce_hex(["12399aedf","abc"])
        '0000000000000000000000012399a463'
    """
    from operator import xor
    res = reduce(xor, (int(x, 16) for x in fingerprints), 0)
    if res < 0:
        res += 1 << 128
    return "%032x" % res


class MarkedOutput(str):
    """
    A subclass of string with context for whether another string
    matches it.

    EXAMPLES::

        sage: from sage.doctest.parsing import MarkedOutput
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
    def update(self, **kwds):
        """
        EXAMPLES::

            sage: from sage.doctest.parsing import MarkedOutput
            sage: s = MarkedOutput("0.0007401")
            sage: s.update(abs_tol = .0000001)
            '0.0007401'
            sage: s.rel_tol
            0
            sage: s.abs_tol
            1.00000000000000e-7
        """
        self.__dict__.update(kwds)
        return self

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: from sage.doctest.parsing import MarkedOutput
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


def make_marked_output(s, D):
    """
    Auxiliary function for pickling.

    EXAMPLES::

        sage: from sage.doctest.parsing import make_marked_output
        sage: s = make_marked_output("0.0007401", {'abs_tol':.0000001})
        sage: s
        '0.0007401'
        sage: s.abs_tol
        1.00000000000000e-7
    """
    ans = MarkedOutput(s)
    ans.__dict__.update(D)
    return ans


class OriginalSource(object):
    r"""
    Context swapping out the pre-parsed source with the original for
    better reporting.

    EXAMPLES::

        sage: from sage.doctest.sources import FileDocTestSource
        sage: from sage.doctest.control import DocTestDefaults
        sage: from sage.env import SAGE_SRC
        sage: import os
        sage: filename = os.path.join(SAGE_SRC,'sage','doctest','forker.py')
        sage: FDS = FileDocTestSource(filename,DocTestDefaults())
        sage: doctests, extras = FDS.create_doctests(globals())
        sage: ex = doctests[0].examples[0]
        sage: ex.sage_source
        'doctest_var = 42; doctest_var^2\n'
        sage: ex.source
        'doctest_var = Integer(42); doctest_var**Integer(2)\n'
        sage: from sage.doctest.parsing import OriginalSource
        sage: with OriginalSource(ex):
        ....:     ex.source
        'doctest_var = 42; doctest_var^2\n'
    """
    def __init__(self, example):
        """
        Swaps out the source for the sage_source of a doctest example.

        INPUT:

        - ``example`` -- a :class:`doctest.Example` instance

        EXAMPLES::

            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: from sage.doctest.parsing import OriginalSource
            sage: OriginalSource(ex)
            <sage.doctest.parsing.OriginalSource object at ...>
        """
        self.example = example

    def __enter__(self):
        r"""
        EXAMPLES::

            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: from sage.doctest.parsing import OriginalSource
            sage: with OriginalSource(ex): # indirect doctest
            ....:     ex.source
            'doctest_var = 42; doctest_var^2\n'
        """
        if hasattr(self.example, 'sage_source'):
            self.old_source, self.example.source = self.example.source, self.example.sage_source

    def __exit__(self, *args):
        r"""
        EXAMPLES::

            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: from sage.doctest.parsing import OriginalSource
            sage: with OriginalSource(ex): # indirect doctest
            ....:     ex.source
            'doctest_var = 42; doctest_var^2\n'
            sage: ex.source # indirect doctest
            'doctest_var = Integer(42); doctest_var**Integer(2)\n'
        """
        if hasattr(self.example, 'sage_source'):
            self.example.source = self.old_source


class SageDocTestParser(doctest.DocTestParser):
    """
    A version of the standard doctest parser which handles Sage's
    custom options and tolerances in floating point arithmetic.
    """
    def __init__(self, optional_tags=(), long=False):
        r"""
        INPUT:

        - ``optional_tags`` -- a list or tuple of strings.
        - ``long`` -- boolean, whether to run doctests marked as taking a
          long time.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageDocTestParser
            sage: DTP = SageDocTestParser(('sage','magma','guava'))
            sage: ex = DTP.parse("sage: 2 + 2\n")[1]
            sage: ex.sage_source
            '2 + 2\n'
            sage: ex = DTP.parse("sage: R.<x> = ZZ[]")[1]
            sage: ex.source
            "R = ZZ['x']; (x,) = R._first_ngens(1)\n"

        TESTS::

            sage: TestSuite(DTP).run()
        """
        self.long = long
        self.optionals = defaultdict(int)  # record skipped optional tests
        if optional_tags is True:  # run all optional tests
            self.optional_tags = True
            self.optional_only = False
        else:
            self.optional_tags = set(optional_tags)
            if 'sage' in self.optional_tags:
                self.optional_only = False
                self.optional_tags.remove('sage')
            else:
                self.optional_only = True

    def __eq__(self, other):
        """
        Comparison.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageDocTestParser
            sage: DTP = SageDocTestParser(('sage','magma','guava'), True)
            sage: DTP2 = SageDocTestParser(('sage','magma','guava'), False)
            sage: DTP == DTP2
            False
        """
        if not isinstance(other, SageDocTestParser):
            return False
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Test for non-equality.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageDocTestParser
            sage: DTP = SageDocTestParser(('sage','magma','guava'), True)
            sage: DTP2 = SageDocTestParser(('sage','magma','guava'), False)
            sage: DTP != DTP2
            True
        """
        return not (self == other)

    def parse(self, string, *args):
        r"""
        A Sage specialization of :class:`doctest.DocTestParser`.

        INPUT:

        - ``string`` -- the string to parse.
        - ``name`` -- optional string giving the name identifying string,
          to be used in error messages.

        OUTPUT:

        - A list consisting of strings and :class:`doctest.Example`
          instances.  There will be at least one string between
          successive examples (exactly one unless or long or optional
          tests are removed), and it will begin and end with a string.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageDocTestParser
            sage: DTP = SageDocTestParser(('sage','magma','guava'))
            sage: example = 'Explanatory text::\n\n    sage: E = magma("EllipticCurve([1, 1, 1, -10, -10])") # optional: magma\n\nLater text'
            sage: parsed = DTP.parse(example)
            sage: parsed[0]
            'Explanatory text::\n\n'
            sage: parsed[1].sage_source
            'E = magma("EllipticCurve([1, 1, 1, -10, -10])") # optional: magma\n'
            sage: parsed[2]
            '\nLater text'

        If the doctest parser is not created to accept a given
        optional argument, the corresponding examples will just be
        removed::

            sage: DTP2 = SageDocTestParser(('sage',))
            sage: parsed2 = DTP2.parse(example)
            sage: parsed2
            ['Explanatory text::\n\n', '\nLater text']

        You can mark doctests as having a particular tolerance::

            sage: example2 = 'sage: gamma(1.6) # tol 2.0e-11\n0.893515349287690'
            sage: ex = DTP.parse(example2)[1]
            sage: ex.sage_source
            'gamma(1.6) # tol 2.0e-11\n'
            sage: ex.want
            '0.893515349287690\n'
            sage: type(ex.want)
            <class 'sage.doctest.parsing.MarkedOutput'>
            sage: ex.want.tol
            2.000000000000000000...?e-11

        You can use continuation lines::

            sage: s = "sage: for i in range(4):\n....:     print(i)\n....:\n"
            sage: ex = DTP2.parse(s)[1]
            sage: ex.source
            'for i in range(Integer(4)):\n    print(i)\n'

        Sage currently accepts backslashes as indicating that the end
        of the current line should be joined to the next line.  This
        feature allows for breaking large integers over multiple lines
        but is not standard for Python doctesting.  It's not
        guaranteed to persist, but works in Sage 5.5::

            sage: n = 1234\
            ....:     5678
            sage: print(n)
            12345678
            sage: type(n)
            <class 'sage.rings.integer.Integer'>

        It also works without the line continuation::

            sage: m = 8765\
            4321
            sage: print(m)
            87654321

        Test that :trac:`26575` is resolved::

            sage: example3 = 'sage: Zp(5,4,print_mode="digits")(5)\n...00010'
            sage: parsed3 = DTP.parse(example3)
            sage: dte = parsed3[1]
            sage: dte.sage_source
            'Zp(5,4,print_mode="digits")(5)\n'
            sage: dte.want
            '...00010\n'
        """
        # Hack for non-standard backslash line escapes accepted by the current
        # doctest system.
        m = backslash_replacer.search(string)
        while m is not None:
            next_prompt = find_sage_prompt.search(string,m.end())
            g = m.groups()
            if next_prompt:
                future = string[m.end():next_prompt.start()] + '\n' + string[next_prompt.start():]
            else:
                future = string[m.end():]
            string = string[:m.start()] + g[0] + "sage:" + g[1] + future
            m = backslash_replacer.search(string,m.start())

        replace_ellipsis = not python_prompt.search(string)
        if replace_ellipsis:
            # There are no >>> prompts, so we can allow ... to begin the output
            # We do so by replacing ellipses with a special tag, then putting them back after parsing
            string = find_python_continuation.sub(r"\1" + ellipsis_tag + r"\2", string)
        string = find_sage_prompt.sub(r"\1>>> sage: ", string)
        string = find_sage_continuation.sub(r"\1...", string)
        res = doctest.DocTestParser.parse(self, string, *args)
        filtered = []
        for item in res:
            if isinstance(item, doctest.Example):
                optional_tags = parse_optional_tags(item.source)
                if optional_tags:
                    for tag in optional_tags:
                        self.optionals[tag] += 1
                    if (('not implemented' in optional_tags) or
                            ('not tested' in optional_tags)):
                        continue

                    if 'long time' in optional_tags:
                        if self.long:
                            optional_tags.remove('long time')
                        else:
                            continue

                    if self.optional_tags is not True:
                        extra = optional_tags - self.optional_tags # set difference
                        if extra:
                            if not available_software.issuperset(extra):
                                continue
                elif self.optional_only:
                    self.optionals['sage'] += 1
                    continue
                if replace_ellipsis:
                    item.want = item.want.replace(ellipsis_tag, "...")
                    if item.exc_msg is not None:
                        item.exc_msg = item.exc_msg.replace(ellipsis_tag, "...")
                item.want = parse_tolerance(item.source, item.want)
                if item.source.startswith("sage: "):
                    item.sage_source = item.source[6:]
                    if item.sage_source.lstrip().startswith('#'):
                        continue
                    item.source = preparse(item.sage_source)
            filtered.append(item)
        return filtered

class SageOutputChecker(doctest.OutputChecker):
    r"""
    A modification of the doctest OutputChecker that can check
    relative and absolute tolerance of answers.

    EXAMPLES::

        sage: from sage.doctest.parsing import SageOutputChecker, MarkedOutput, SageDocTestParser
        sage: import doctest
        sage: optflag = doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS
        sage: DTP = SageDocTestParser(('sage','magma','guava'))
        sage: OC = SageOutputChecker()
        sage: example2 = 'sage: gamma(1.6) # tol 2.0e-11\n0.893515349287690'
        sage: ex = DTP.parse(example2)[1]
        sage: ex.sage_source
        'gamma(1.6) # tol 2.0e-11\n'
        sage: ex.want
        '0.893515349287690\n'
        sage: type(ex.want)
        <class 'sage.doctest.parsing.MarkedOutput'>
        sage: ex.want.tol
        2.000000000000000000...?e-11
        sage: OC.check_output(ex.want, '0.893515349287690', optflag)
        True
        sage: OC.check_output(ex.want, '0.8935153492877', optflag)
        True
        sage: OC.check_output(ex.want, '0', optflag)
        False
        sage: OC.check_output(ex.want, 'x + 0.8935153492877', optflag)
        False
    """
    def human_readable_escape_sequences(self, string):
        r"""
        Make ANSI escape sequences human readable.

        EXAMPLES::

            sage: print('This is \x1b[1mbold\x1b[0m text')
            This is <CSI-1m>bold<CSI-0m> text

        TESTS::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: OC = SageOutputChecker()
            sage: teststr = '-'.join([
            ....:     'bold\x1b[1m',
            ....:     'red\x1b[31m',
            ....:     'oscmd\x1ba'])
            sage: OC.human_readable_escape_sequences(teststr)
            'bold<CSI-1m>-red<CSI-31m>-oscmd<ESC-a>'
        """
        def human_readable(match):
            ansi_escape = match.group(1)
            assert len(ansi_escape) >= 2
            if len(ansi_escape) == 2:
                return '<ESC-'+ansi_escape[1]+'>'
            else:
                return '<CSI-'+ansi_escape.lstrip('\x1b[\x9b')+'>'
        return ansi_escape_sequence.subn(human_readable, string)[0]

    def add_tolerance(self, wantval, want):
        """
        Enlarge the real interval element ``wantval`` according to
        the tolerance options in ``want``.

        INPUT:

        - ``wantval`` -- a real interval element
        - ``want`` -- a :class:`MarkedOutput` describing the tolerance

        OUTPUT:

        - an interval element containing ``wantval``

        EXAMPLES::

            sage: from sage.doctest.parsing import MarkedOutput, SageOutputChecker
            sage: OC = SageOutputChecker()
            sage: want_tol = MarkedOutput().update(tol=0.0001)
            sage: want_abs = MarkedOutput().update(abs_tol=0.0001)
            sage: want_rel = MarkedOutput().update(rel_tol=0.0001)
            sage: OC.add_tolerance(RIF(pi.n(64)), want_tol).endpoints()
            (3.14127849432443, 3.14190681285516)
            sage: OC.add_tolerance(RIF(pi.n(64)), want_abs).endpoints()
            (3.14149265358979, 3.14169265358980)
            sage: OC.add_tolerance(RIF(pi.n(64)), want_rel).endpoints()
            (3.14127849432443, 3.14190681285516)
            sage: OC.add_tolerance(RIF(1e1000), want_tol)
            1.000?e1000
            sage: OC.add_tolerance(RIF(1e1000), want_abs)
            1.000000000000000?e1000
            sage: OC.add_tolerance(RIF(1e1000), want_rel)
            1.000?e1000
            sage: OC.add_tolerance(0, want_tol)
            0.000?
            sage: OC.add_tolerance(0, want_abs)
            0.000?
            sage: OC.add_tolerance(0, want_rel)
            0
        """
        if want.tol:
            if wantval == 0:
                return RIFtol(want.tol) * RIFtol(-1,1)
            else:
                return wantval * (1 + RIFtol(want.tol) * RIFtol(-1,1))
        elif want.abs_tol:
            return wantval + RIFtol(want.abs_tol) * RIFtol(-1,1)
        elif want.rel_tol:
            return wantval * (1 + RIFtol(want.rel_tol) * RIFtol(-1,1))
        else:
            return wantval

    def check_output(self, want, got, optionflags):
        """
        Checks to see if the output matches the desired output.

        If ``want`` is a :class:`MarkedOutput` instance, takes into account the desired tolerance.

        INPUT:

        - ``want`` -- a string or :class:`MarkedOutput`
        - ``got`` -- a string
        - ``optionflags`` -- an integer, passed down to :class:`doctest.OutputChecker`

        OUTPUT:

        - boolean, whether ``got`` matches ``want`` up to the specified tolerance.

        EXAMPLES::

            sage: from sage.doctest.parsing import MarkedOutput, SageOutputChecker
            sage: import doctest
            sage: optflag = doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS
            sage: rndstr = MarkedOutput("I'm wrong!").update(random=True)
            sage: tentol = MarkedOutput("10.0").update(tol=.1)
            sage: tenabs = MarkedOutput("10.0").update(abs_tol=.1)
            sage: tenrel = MarkedOutput("10.0").update(rel_tol=.1)
            sage: zerotol = MarkedOutput("0.0").update(tol=.1)
            sage: zeroabs = MarkedOutput("0.0").update(abs_tol=.1)
            sage: zerorel = MarkedOutput("0.0").update(rel_tol=.1)
            sage: zero = "0.0"
            sage: nf = "9.5"
            sage: ten = "10.05"
            sage: eps = "-0.05"
            sage: OC = SageOutputChecker()

        ::

            sage: OC.check_output(rndstr,nf,optflag)
            True

            sage: OC.check_output(tentol,nf,optflag)
            True
            sage: OC.check_output(tentol,ten,optflag)
            True
            sage: OC.check_output(tentol,zero,optflag)
            False

            sage: OC.check_output(tenabs,nf,optflag)
            False
            sage: OC.check_output(tenabs,ten,optflag)
            True
            sage: OC.check_output(tenabs,zero,optflag)
            False

            sage: OC.check_output(tenrel,nf,optflag)
            True
            sage: OC.check_output(tenrel,ten,optflag)
            True
            sage: OC.check_output(tenrel,zero,optflag)
            False

            sage: OC.check_output(zerotol,zero,optflag)
            True
            sage: OC.check_output(zerotol,eps,optflag)
            True
            sage: OC.check_output(zerotol,ten,optflag)
            False

            sage: OC.check_output(zeroabs,zero,optflag)
            True
            sage: OC.check_output(zeroabs,eps,optflag)
            True
            sage: OC.check_output(zeroabs,ten,optflag)
            False

            sage: OC.check_output(zerorel,zero,optflag)
            True
            sage: OC.check_output(zerorel,eps,optflag)
            False
            sage: OC.check_output(zerorel,ten,optflag)
            False

        More explicit tolerance checks::

            sage: _ = x  # rel tol 1e10
            sage: raise RuntimeError   # rel tol 1e10
            Traceback (most recent call last):
            ...
            RuntimeError
            sage: 1  # abs tol 2
            -0.5
            sage: print("0.9999")    # rel tol 1e-4
            1.0
            sage: print("1.00001")   # abs tol 1e-5
            1.0
            sage: 0  # rel tol 1
            1

        Spaces before numbers or between the sign and number are ignored::

            sage: print("[ - 1, 2]")  # abs tol 1e-10
            [-1,2]

        Tolerance on Python 3 for string results with unicode prefix::

            sage: a = 'Cyrano'; a
            'Cyrano'
            sage: b = ['Fermat', 'Euler']; b
            ['Fermat',  'Euler']
            sage: c = 'you'; c
            'you'
        """
        got = self.human_readable_escape_sequences(got)
        got = glpk_simplex_warning_regex.sub('', got)
        got = ld_warning_regex.sub('', got)
        if isinstance(want, MarkedOutput):
            if want.random:
                return True
            elif want.tol or want.rel_tol or want.abs_tol:
                # First check the doctest without the numbers
                want_str = [g[0] for g in float_regex.findall(want)]
                got_str = [g[0] for g in float_regex.findall(got)]
                if len(want_str) != len(got_str):
                    return False
                starwant = float_regex.sub('*', want)
                stargot = float_regex.sub('*', got)
                if not doctest.OutputChecker.check_output(self, starwant, stargot, optionflags):
                    return False

                # Now check the numbers
                want_values = [RIFtol(g) for g in want_str]
                want_intervals = [self.add_tolerance(v, want) for v in want_values]
                got_values = [RIFtol(g) for g in got_str]
                # The doctest is successful if the "want" and "got"
                # intervals have a non-empty intersection
                return all(a.overlaps(b) for a, b in zip(want_intervals, got_values))

        ok = doctest.OutputChecker.check_output(self, want, got, optionflags)

        if ok:
            return ok

        did_fixup = False
        for quick_check, fixup in _repr_fixups:
            do_fixup = quick_check(got, want)
            if do_fixup:
                got, want = fixup(got, want)
                did_fixup = True

        if not did_fixup:
            # Return the same result as before
            return ok

        return doctest.OutputChecker.check_output(self, want, got, optionflags)

    def output_difference(self, example, got, optionflags):
        r"""
        Report on the differences between the desired result and what
        was actually obtained.

        If ``want`` is a :class:`MarkedOutput` instance, takes into account the desired tolerance.

        INPUT:

        - ``example`` -- a :class:`doctest.Example` instance
        - ``got`` -- a string
        - ``optionflags`` -- an integer, passed down to :class:`doctest.OutputChecker`

        OUTPUT:

        - a string, describing how ``got`` fails to match ``example.want``

        EXAMPLES::

            sage: from sage.doctest.parsing import MarkedOutput, SageOutputChecker
            sage: import doctest
            sage: optflag = doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS
            sage: tentol = doctest.Example('',MarkedOutput("10.0\n").update(tol=.1))
            sage: tenabs = doctest.Example('',MarkedOutput("10.0\n").update(abs_tol=.1))
            sage: tenrel = doctest.Example('',MarkedOutput("10.0\n").update(rel_tol=.1))
            sage: zerotol = doctest.Example('',MarkedOutput("0.0\n").update(tol=.1))
            sage: zeroabs = doctest.Example('',MarkedOutput("0.0\n").update(abs_tol=.1))
            sage: zerorel = doctest.Example('',MarkedOutput("0.0\n").update(rel_tol=.1))
            sage: tlist = doctest.Example('',MarkedOutput("[10.0, 10.0, 10.0, 10.0, 10.0, 10.0]\n").update(abs_tol=0.987))
            sage: zero = "0.0"
            sage: nf = "9.5"
            sage: ten = "10.05"
            sage: eps = "-0.05"
            sage: L = "[9.9, 8.7, 10.3, 11.2, 10.8, 10.0]"
            sage: OC = SageOutputChecker()

        ::

            sage: print(OC.output_difference(tenabs,nf,optflag))
            Expected:
                10.0
            Got:
                9.5
            Tolerance exceeded:
                10.0 vs 9.5, tolerance 5e-1 > 1e-1

            sage: print(OC.output_difference(tentol,zero,optflag))
            Expected:
                10.0
            Got:
                0.0
            Tolerance exceeded:
                10.0 vs 0.0, tolerance 1e0 > 1e-1

            sage: print(OC.output_difference(tentol,eps,optflag))
            Expected:
                10.0
            Got:
                -0.05
            Tolerance exceeded:
                10.0 vs -0.05, tolerance 2e0 > 1e-1

            sage: print(OC.output_difference(tlist,L,optflag))
            Expected:
                [10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
            Got:
                [9.9, 8.7, 10.3, 11.2, 10.8, 10.0]
            Tolerance exceeded in 2 of 6:
                10.0 vs 8.7, tolerance 2e0 > 9.87e-1
                10.0 vs 11.2, tolerance 2e0 > 9.87e-1

        TESTS::

            sage: print(OC.output_difference(tenabs,zero,optflag))
            Expected:
                10.0
            Got:
                0.0
            Tolerance exceeded:
                10.0 vs 0.0, tolerance 1e1 > 1e-1

            sage: print(OC.output_difference(tenrel,zero,optflag))
            Expected:
                10.0
            Got:
                0.0
            Tolerance exceeded:
                10.0 vs 0.0, tolerance 1e0 > 1e-1

            sage: print(OC.output_difference(tenrel,eps,optflag))
            Expected:
                10.0
            Got:
                -0.05
            Tolerance exceeded:
                10.0 vs -0.05, tolerance 2e0 > 1e-1

            sage: print(OC.output_difference(zerotol,ten,optflag))
            Expected:
                0.0
            Got:
                10.05
            Tolerance exceeded:
                0.0 vs 10.05, tolerance 2e1 > 1e-1

            sage: print(OC.output_difference(zeroabs,ten,optflag))
            Expected:
                0.0
            Got:
                10.05
            Tolerance exceeded:
                0.0 vs 10.05, tolerance 2e1 > 1e-1

            sage: print(OC.output_difference(zerorel,eps,optflag))
            Expected:
                0.0
            Got:
                -0.05
            Tolerance exceeded:
                0.0 vs -0.05, tolerance +infinity > 1e-1

            sage: print(OC.output_difference(zerorel,ten,optflag))
            Expected:
                0.0
            Got:
                10.05
            Tolerance exceeded:
                0.0 vs 10.05, tolerance +infinity > 1e-1
        """
        got = self.human_readable_escape_sequences(got)
        want = example.want
        diff = doctest.OutputChecker.output_difference(self, example, got, optionflags)
        if isinstance(want, MarkedOutput) and (want.tol or want.abs_tol or want.rel_tol):
            if diff[-1] != "\n":
                diff += "\n"
            want_str = [g[0] for g in float_regex.findall(want)]
            got_str = [g[0] for g in float_regex.findall(got)]
            if len(want_str) == len(got_str):
                failures = []
                def fail(x, y, actual, desired):
                    failstr = "    {} vs {}, tolerance {} > {}".format(x, y,
                        RIFtol(actual).upper().str(digits=1, no_sci=False),
                        RIFtol(desired).center().str(digits=15, skip_zeroes=True, no_sci=False)
                    )
                    failures.append(failstr)

                for wstr, gstr in zip(want_str, got_str):
                    w = RIFtol(wstr)
                    g = RIFtol(gstr)
                    if not g.overlaps(self.add_tolerance(w, want)):
                        if want.tol:
                            if not w:
                                fail(wstr, gstr, abs(g), want.tol)
                            else:
                                fail(wstr, gstr, abs(1 - g/w), want.tol)
                        elif want.abs_tol:
                            fail(wstr, gstr, abs(g - w), want.abs_tol)
                        else:
                            fail(wstr, gstr, abs(1 - g/w), want.rel_tol)

                if failures:
                    if len(want_str) == 1:
                        diff += "Tolerance exceeded:\n"
                    else:
                        diff += "Tolerance exceeded in %s of %s:\n"%(len(failures), len(want_str))
                    diff += "\n".join(failures) + "\n"
            elif "..." in want:
                diff += "Note: combining tolerance (# tol) with ellipsis (...) is not supported\n"
        return diff
