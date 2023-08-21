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

import collections.abc
import doctest
import re

from collections import defaultdict
from functools import reduce

from sage.misc.cachefunc import cached_function
from sage.repl.preparse import preparse, strip_string_literals

from .external import available_software, external_software

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
            from sage.rings.real_mpfi import RealIntervalField
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
ansi_escape_sequence = re.compile(r'(\x1b[@-Z\\-~]|\x1b\[.*?[@-~]|\x9b.*?[@-~])')

special_optional_regex = 'arb216|arb218|py2|long time|not implemented|not tested|known bug'
tag_with_explanation_regex = fr'((?:\w|[.])+)\s*(?:\((.*?)\))?'
optional_regex = re.compile(fr'(?P<cmd>{special_optional_regex})\s*(?:\((?P<cmd_explanation>.*?)\))?|'
                            fr'[^ a-z]\s*(optional|needs)(?:\s|[:-])*(?P<tags>(?:(?:{tag_with_explanation_regex})\s*)*)',
                            re.IGNORECASE)
special_optional_regex = re.compile(special_optional_regex, re.IGNORECASE)
tag_with_explanation_regex = re.compile(tag_with_explanation_regex, re.IGNORECASE)

nodoctest_regex = re.compile(r'\s*(#+|%+|r"+|"+|\.\.)\s*nodoctest')
optionaltag_regex = re.compile(r'^(\w|[.])+$')
optionalfiledirective_regex = re.compile(r'\s*(#+|%+|r"+|"+|\.\.)\s*sage\.doctest: (.*)')


def parse_optional_tags(string, *, return_string_sans_tags=False):
    r"""
    Return a dictionary whose keys are optional tags from the following
    set that occur in a comment on the first line of the input string.

    - ``'long time'``
    - ``'not implemented'``
    - ``'not tested'``
    - ``'known bug'``
    - ``'py2'``
    - ``'arb216'``
    - ``'arb218'``
    - ``'optional - FEATURE...'`` or ``'needs FEATURE...'`` --
      the dictionary will just have the key ``'FEATURE'``

    The values, if non-``None``, are strings with optional explanations
    for a tag, which may appear in parentheses after the tag in ``string``.

    INPUT:

    - ``string`` -- a string

    - ``return_string_sans_tags`` -- (boolean, default ``False``); whether to
      additionally return ``string`` with the optional tags removed but other
      comments kept and a boolean ``is_persistent``

    EXAMPLES::

        sage: from sage.doctest.parsing import parse_optional_tags
        sage: parse_optional_tags("sage: magma('2 + 2')# optional: magma")
        {'magma': None}
        sage: parse_optional_tags("sage: #optional -- mypkg")
        {'mypkg': None}
        sage: parse_optional_tags("sage: print(1)  # parentheses are optional here")
        {}
        sage: parse_optional_tags("sage: print(1)  # optional")
        {}
        sage: sorted(list(parse_optional_tags("sage: #optional -- foo bar, baz")))
        ['bar', 'foo']
        sage: parse_optional_tags("sage: #optional -- foo.bar, baz")
        {'foo.bar': None}
        sage: parse_optional_tags("sage: #needs foo.bar, baz")
        {'foo.bar': None}
        sage: sorted(list(parse_optional_tags("    sage: factor(10^(10^10) + 1) # LoNg TiME, NoT TeSTED; OptioNAL -- P4cka9e")))
        ['long time', 'not tested', 'p4cka9e']
        sage: parse_optional_tags("    sage: raise RuntimeError # known bug")
        {'bug': None}
        sage: sorted(list(parse_optional_tags("    sage: determine_meaning_of_life() # long time, not implemented")))
        ['long time', 'not implemented']

    We don't parse inside strings::

        sage: parse_optional_tags("    sage: print('  # long time')")
        {}
        sage: parse_optional_tags("    sage: print('  # long time')  # not tested")
        {'not tested': None}

    UTF-8 works::

        sage: parse_optional_tags("'ěščřžýáíéďĎ'")
        {}

    Tags with parenthesized explanations::

        sage: parse_optional_tags("    sage: 1 + 1  # long time (1 year, 2 months??), optional - bliss (because)")
        {'bliss': 'because', 'long time': '1 year, 2 months??'}

    With ``return_string_sans_tags=True``::

        sage: parse_optional_tags("sage: print(1)  # very important 1  # optional - foo",
        ....:                     return_string_sans_tags=True)
        ({'foo': None}, 'sage: print(1)  # very important 1  ', False)
        sage: parse_optional_tags("sage: print(    # very important too  # optional - foo\n....:     2)",
        ....:                     return_string_sans_tags=True)
        ({'foo': None}, 'sage: print(    # very important too  \n....:     2)', False)
        sage: parse_optional_tags("sage: #this is persistent #needs scipy",
        ....:                     return_string_sans_tags=True)
        ({'scipy': None}, 'sage: #this is persistent ', True)
        sage: parse_optional_tags("sage: #this is not #needs scipy\n....: import scipy",
        ....:                     return_string_sans_tags=True)
        ({'scipy': None}, 'sage: #this is not \n....: import scipy', False)

    """
    safe, literals, state = strip_string_literals(string)
    split = safe.split('\n', 1)
    if len(split) > 1:
        first_line, rest = split
    else:
        first_line, rest = split[0], None

    sharp_index = first_line.find('#')
    if sharp_index < 0:                  # no comment
        if return_string_sans_tags:
            return {}, string, False
        else:
            return {}

    first_line_sans_comments, comment = first_line[:sharp_index] % literals, first_line[sharp_index:] % literals

    if return_string_sans_tags:
        # skip non-tag comments that precede the first tag comment
        if m := optional_regex.search(comment):
            sharp_index = comment[:m.start(0) + 1].rfind('#')
            if sharp_index >= 0:
                first_line = first_line_sans_comments + comment[:sharp_index]
                comment = comment[sharp_index:]
        else:
            # no tag comment
            return {}, string, False

    tags = {}
    for m in optional_regex.finditer(comment):
        cmd = m.group('cmd')
        if cmd and cmd.lower() == 'known bug':
            tags['bug'] = None  # so that such tests will be run by sage -t ... -only-optional=bug
        elif cmd:
            tags[cmd.lower()] = m.group('cmd_explanation') or None
        else:
            # optional/needs
            tags.update({m.group(1).lower(): m.group(2) or None
                         for m in tag_with_explanation_regex.finditer(m.group('tags'))})

    if return_string_sans_tags:
        is_persistent = tags and first_line_sans_comments.strip() == 'sage:' and not rest  # persistent (block-scoped) tag
        return tags, (first_line + '\n' + rest%literals if rest is not None
                      else first_line), is_persistent
    else:
        return tags


def parse_file_optional_tags(lines):
    r"""
    Scan the first few lines for file-level doctest directives.

    INPUT:

    - ``lines`` -- iterable of pairs ``(lineno, line)``.

    OUTPUT:

    a dictionary whose keys are strings (tags); see :func:`parse_optional_tags`

    EXAMPLES::

        sage: from sage.doctest.parsing import parse_file_optional_tags
        sage: filename = tmp_filename(ext=".pyx")
        sage: with open(filename, "r") as f:
        ....:     parse_file_optional_tags(enumerate(f))
        {}
        sage: with open(filename, "w") as f:
        ....:     _ = f.write("# nodoctest")
        sage: with open(filename, "r") as f:
        ....:     parse_file_optional_tags(enumerate(f))
        {'not tested': None}
        sage: with open(filename, "w") as f:
        ....:     _ = f.write("# sage.doctest: "    # broken in two source lines to avoid the pattern
        ....:                 "optional - xyz")     # of relint (multiline_doctest_comment)
        sage: with open(filename, "r") as f:
        ....:     parse_file_optional_tags(enumerate(f))
        {'xyz': None}
    """
    tags = {}
    for line_count, line in lines:
        if nodoctest_regex.match(line):
            tags['not tested'] = None
        if m := optionalfiledirective_regex.match(line):
            file_tag_string = m.group(2)
            tags.update(parse_optional_tags('#' + file_tag_string))
        if line_count >= 10:
            break
    return tags


@cached_function
def _standard_tags():
    r"""
    Return the set of the names of all standard features.

    EXAMPLES::

        sage: from sage.doctest.parsing import _standard_tags
        sage: sorted(_standard_tags())
        [..., 'numpy', ..., 'sage.rings.finite_rings', ...]
    """
    from sage.features.all import all_features
    return frozenset(feature.name for feature in all_features()
                     if feature._spkg_type() == 'standard')


def _tag_group(tag):
    r"""
    Classify a doctest tag as belonging to one of 4 groups.

    INPUT:

    - ``tag`` -- string

    OUTPUT:

    a string; one of ``'special'``, ``'optional'``, ``'standard'``, ``'sage'``

    EXAMPLES::

        sage: from sage.doctest.parsing import _tag_group
        sage: _tag_group('scipy')
        'standard'
        sage: _tag_group('sage.numerical.mip')
        'sage'
        sage: _tag_group('bliss')
        'optional'
        sage: _tag_group('not tested')
        'special'
    """
    if tag.startswith('sage.'):
        return 'sage'
    elif tag in _standard_tags():
        return 'standard'
    elif not special_optional_regex.fullmatch(tag):
        return 'optional'
    else:
        return 'special'


def unparse_optional_tags(tags, prefix='# '):
    r"""
    Return a comment string that sets ``tags``.

    INPUT:

    - ``tags`` -- dict or iterable of tags, as output by :func:`parse_optional_tags`

    - ``prefix`` -- to be put before a nonempty string

    EXAMPLES::

        sage: from sage.doctest.parsing import unparse_optional_tags
        sage: unparse_optional_tags({})
        ''
        sage: unparse_optional_tags({'magma': None})
        '# optional - magma'
        sage: unparse_optional_tags({'fictional_optional': None,
        ....:                        'sage.rings.number_field': None,
        ....:                        'scipy': 'just because',
        ....:                        'bliss': None})
        '# optional - bliss fictional_optional, needs scipy (just because) sage.rings.number_field'
        sage: unparse_optional_tags(['long time', 'not tested', 'p4cka9e'], prefix='')
        'long time, not tested, optional - p4cka9e'
    """
    group = defaultdict(set)
    if isinstance(tags, collections.abc.Mapping):
        for tag, explanation in tags.items():
            if tag == 'bug':
                tag = 'known bug'
            group[_tag_group(tag)].add(f'{tag} ({explanation})' if explanation else tag)
    else:
        for tag in tags:
            if tag == 'bug':
                tag = 'known bug'
            group[_tag_group(tag)].add(tag)

    tags = sorted(group.pop('special', []))
    if 'optional' in group:
        tags.append('optional - ' + " ".join(sorted(group.pop('optional'))))
    if 'standard' in group or 'sage' in group:
        tags.append('needs ' + " ".join(sorted(group.pop('standard', []))
                                        + sorted(group.pop('sage', []))))
    assert not group
    if tags:
        return prefix + ', '.join(tags)
    return ''


optional_tag_columns = [48, 56, 64, 72, 80, 84]
standard_tag_columns = [88, 100, 120, 160]


def update_optional_tags(line, tags=None, *, add_tags=None, remove_tags=None, force_rewrite=False):
    r"""
    Return the doctest ``line`` with tags changed.

    EXAMPLES::

        sage: from sage.doctest.parsing import update_optional_tags, optional_tag_columns, standard_tag_columns
        sage: ruler = ''
        sage: for column in optional_tag_columns:
        ....:     ruler += ' ' * (column - len(ruler)) + 'V'
        sage: for column in standard_tag_columns:
        ....:     ruler += ' ' * (column - len(ruler)) + 'v'
        sage: def print_with_ruler(lines):
        ....:     print('|' + ruler)
        ....:     for line in lines:
        ....:         print('|' + line)
        sage: print_with_ruler([  # the tags are obscured in the source file to avoid relint warnings
        ....:     update_optional_tags('    sage: something()  # opt' 'ional - latte_int',
        ....:                          remove_tags=['latte_int', 'wasnt_even_there']),
        ....:     update_optional_tags('    sage: nothing_to_be_seen_here()',
        ....:                          tags=['scipy', 'long time']),
        ....:     update_optional_tags('    sage: nothing_to_be_seen_here(honestly=True)',
        ....:                          add_tags=['scipy', 'long time']),
        ....:     update_optional_tags('    sage: nothing_to_be_seen_here(honestly=True, very=True)',
        ....:                          add_tags=['scipy', 'long time']),
        ....:     update_optional_tags('    sage: no_there_is_absolutely_nothing_to_be_seen_here_i_am_serious()#opt' 'ional:bliss',
        ....:                          add_tags=['scipy', 'long time']),
        ....:     update_optional_tags('    sage: ntbsh()  # abbrv for above#opt' 'ional:bliss',
        ....:                          add_tags={'scipy': None, 'long time': '30s on the highest setting'}),
        ....:     update_optional_tags('    sage: no_there_is_absolutely_nothing_to_be_seen_here_i_am_serious()  # really, you can trust me here',
        ....:                          add_tags=['scipy']),
        ....: ])
        |                                                V       V       V       V       V   V   v           v                   v                                       v
        |    sage: something()
        |    sage: nothing_to_be_seen_here()             # long time                             # needs scipy
        |    sage: nothing_to_be_seen_here(honestly=True)        # long time                     # needs scipy
        |    sage: nothing_to_be_seen_here(honestly=True, very=True)     # long time             # needs scipy
        |    sage: no_there_is_absolutely_nothing_to_be_seen_here_i_am_serious()         # long time, optional - bliss, needs scipy
        |    sage: ntbsh()  # abbrv for above            # long time (30s on the highest setting), optional - bliss, needs scipy
        |    sage: no_there_is_absolutely_nothing_to_be_seen_here_i_am_serious()  # really, you can trust me here                # needs scipy

    When no tags are changed, by default, the unchanged input is returned.
    We can force a rewrite; unconditionally or whenever standard tags are involved.
    But even when forced, if comments are already aligned at one of the standard alignment columns,
    this alignment is kept even if we would normally realign farther to the left::

        sage: print_with_ruler([
        ....:     update_optional_tags('    sage: unforced()       # opt' 'ional - latte_int'),
        ....:     update_optional_tags('    sage: unforced()  # opt' 'ional - latte_int',
        ....:                          add_tags=['latte_int']),
        ....:     update_optional_tags('    sage: forced()#opt' 'ional- latte_int',
        ....:                          force_rewrite=True),
        ....:     update_optional_tags('    sage: forced()  # opt' 'ional - scipy',
        ....:                          force_rewrite='standard'),
        ....:     update_optional_tags('    sage: aligned_with_below()                                  # opt' 'ional - 4ti2',
        ....:                          force_rewrite=True),
        ....:     update_optional_tags('    sage: aligned_with_above()                                  # opt' 'ional - 4ti2',
        ....:                          force_rewrite=True),
        ....:     update_optional_tags('    sage: also_already_aligned()                                                                                        # ne' 'eds scipy',
        ....:                          force_rewrite='standard'),
        ....:     update_optional_tags('    sage: two_columns_first_preserved()         # lo' 'ng time                             # ne' 'eds scipy',
        ....:                          force_rewrite='standard'),
        ....:     update_optional_tags('    sage: two_columns_first_preserved()                 # lo' 'ng time                                 # ne' 'eds scipy',
        ....:                          force_rewrite='standard'),
        ....: ])
        |                                                V       V       V       V       V   V   v           v                   v                                       v
        |    sage: unforced()       # optional - latte_int
        |    sage: unforced()  # optional - latte_int
        |    sage: forced()                              # optional - latte_int
        |    sage: forced()                                                                      # needs scipy
        |    sage: aligned_with_below()                                  # optional - 4ti2
        |    sage: aligned_with_above()                                  # optional - 4ti2
        |    sage: also_already_aligned()                                                                                        # needs scipy
        |    sage: two_columns_first_preserved()         # long time                             # needs scipy
        |    sage: two_columns_first_preserved()                 # long time                     # needs scipy

    Rewriting a persistent (block-scoped) tag::

        sage: print_with_ruler([
        ....:     update_optional_tags('    sage:    #opt' 'ional:magma sage.symbolic',
        ....:                          force_rewrite=True),
        ....: ])
        |                                                V       V       V       V       V   V   v           v                   v                                       v
        |    sage: # optional - magma, needs sage.symbolic
    """
    if not (m := re.match('( *sage: *)(.*)', line)):
        raise ValueError(f'line must start with a sage: prompt, got: {line}')

    current_tags, line_sans_tags, is_persistent = parse_optional_tags(line.rstrip(), return_string_sans_tags=True)

    if isinstance(tags, collections.abc.Mapping):
        new_tags = dict(tags)
    elif tags is not None:
        new_tags = {tag: None for tag in tags}
    else:
        new_tags = dict(current_tags)

    if add_tags is not None:
        if isinstance(add_tags, collections.abc.Mapping):
            new_tags.update(add_tags)
        else:
            new_tags.update({tag: None for tag in add_tags})

    if remove_tags is not None:
        for tag in remove_tags:
            new_tags.pop(tag, None)

    if not force_rewrite and new_tags == current_tags:
        return line

    if not new_tags:
        return line_sans_tags.rstrip()

    if (force_rewrite == 'standard'
            and new_tags == current_tags
            and not any(_tag_group(tag) in ['standard', 'sage']
                        for tag in new_tags)):
        return line

    if is_persistent:
        line = line_sans_tags.rstrip() + ' '
    else:
        group = defaultdict(set)
        for tag in new_tags:
            group[_tag_group(tag)].add(tag)
        tag_columns = (optional_tag_columns if group['optional'] or group['special']
                       else standard_tag_columns)

        if len(line_sans_tags) in tag_columns and line_sans_tags[-2:] == '  ':
            # keep alignment
            line = line_sans_tags
            pass
        else:
            # realign
            line = line_sans_tags.rstrip()
            for column in tag_columns:
                if len(line) <= column - 2:
                    line += ' ' * (column - 2 - len(line))
                    break
            line += '  '

        if (group['optional'] or group['special']) and (group['standard'] or group['sage']):
            # Try if two-column mode works better
            first_part = unparse_optional_tags({tag: explanation
                                                for tag, explanation in new_tags.items()
                                                if (tag in group['optional']
                                                    or tag in group['special'])})
            column = standard_tag_columns[0]
            if len(line + first_part) + 8 <= column:
                line += first_part
                line += ' ' * (column - len(line))
                line += unparse_optional_tags({tag: explanation
                                               for tag, explanation in new_tags.items()
                                               if not (tag in group['optional']
                                                       or tag in group['special'])})
                return line.rstrip()

    line += unparse_optional_tags(new_tags)
    return line


def parse_tolerance(source, want):
    r"""
    Return a version of ``want`` marked up with the tolerance tags
    specified in ``source``.

    INPUT:

    - ``source`` -- a string, the source of a doctest
    - ``want`` -- a string, the desired output of the doctest

    OUTPUT:

    ``want`` if there are no tolerance tags specified; a
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
    # regular expressions
    random_marker = re.compile('.*random', re.I)
    tolerance_pattern = re.compile(r'\b((?:abs(?:olute)?)|(?:rel(?:ative)?))? *?tol(?:erance)?\b( +[0-9.e+-]+)?')

    safe, literals, state = strip_string_literals(source)
    first_line = safe.split('\n', 1)[0]
    if '#' not in first_line:
        return want
    comment = first_line[first_line.find('#') + 1:]
    comment = comment[comment.index('(') + 1: comment.rindex(')')]
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


class OriginalSource():
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
    def __init__(self, optional_tags=(), long=False, *, probed_tags=(), file_optional_tags=()):
        r"""
        INPUT:

        - ``optional_tags`` -- a list or tuple of strings.
        - ``long`` -- boolean, whether to run doctests marked as taking a
          long time.
        - ``probed_tags`` -- a list or tuple of strings.
        - ``file_optional_tags`` -- an iterable of strings.

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
        self.probed_tags = probed_tags
        self.file_optional_tags = set(file_optional_tags)

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
          successive examples (exactly one unless long or optional
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
        guaranteed to persist::

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

        Optional tags at the start of an example block persist to the end of
        the block (delimited by a blank line)::

            sage: # long time, needs sage.rings.number_field
            sage: QQbar(I)^10000
            1
            sage: QQbar(I)^10000  # not tested
            I

            sage: # needs sage.rings.finite_rings
            sage: GF(7)
            Finite Field of size 7
            sage: GF(10)
            Traceback (most recent call last):
            ...
            ValueError: the order of a finite field must be a prime power

        Test that :trac:`26575` is resolved::

            sage: example3 = 'sage: Zp(5,4,print_mode="digits")(5)\n...00010'
            sage: parsed3 = DTP.parse(example3)
            sage: dte = parsed3[1]
            sage: dte.sage_source
            'Zp(5,4,print_mode="digits")(5)\n'
            sage: dte.want
            '...00010\n'

        Style warnings::

            sage: def parse(test_string):
            ....:     return [x if isinstance(x, str)
            ....:               else (getattr(x, 'warnings', None), x.sage_source, x.source)
            ....:             for x in DTP.parse(test_string)]

            sage: parse('sage: 1 # optional guava mango\nsage: 2 # optional guava\nsage: 3 # optional guava\nsage: 4 # optional guava\nsage: 5 # optional guava\n\nsage: 11 # optional guava')
            ['',
             (["Consider using a block-scoped tag by inserting the line 'sage: # optional - guava' just before this line to avoid repeating the tag 5 times"],
              '1 # optional guava mango\n',
              'None  # virtual doctest'),
             '',
             (None, '2 # optional guava\n', 'Integer(2) # optional guava\n'),
             '',
             (None, '3 # optional guava\n', 'Integer(3) # optional guava\n'),
             '',
             (None, '4 # optional guava\n', 'Integer(4) # optional guava\n'),
             '',
             (None, '5 # optional guava\n', 'Integer(5) # optional guava\n'),
             '\n',
             (None, '11 # optional guava\n', 'Integer(11) # optional guava\n'),
             '']

            sage: parse('sage: 1 # optional guava\nsage: 2 # optional guava mango\nsage: 3 # optional guava\nsage: 4 # optional guava\nsage: 5 # optional guava\n')
            ['',
             (["Consider using a block-scoped tag by inserting the line 'sage: # optional - guava' just before this line to avoid repeating the tag 5 times"],
              '1 # optional guava\n',
              'Integer(1) # optional guava\n'),
             '',
             '',
             (None, '3 # optional guava\n', 'Integer(3) # optional guava\n'),
             '',
             (None, '4 # optional guava\n', 'Integer(4) # optional guava\n'),
             '',
             (None, '5 # optional guava\n', 'Integer(5) # optional guava\n'),
             '']

            sage: parse('sage: # optional mango\nsage: 1 # optional guava\nsage: 2 # optional guava mango\nsage: 3 # optional guava\nsage: 4 # optional guava\n sage: 5 # optional guava\n')  # optional - guava mango
            ['',
             (["Consider updating this block-scoped tag to 'sage: # optional - guava mango' to avoid repeating the tag 5 times"],
              '# optional mango\n',
              'None  # virtual doctest'),
             '',
             '',
             '',
             '',
             '',
             '']

            sage: parse('::\n\n    sage: 1 # optional guava\n    sage: 2 # optional guava mango\n    sage: 3 # optional guava\n\n::\n\n    sage: 4 # optional guava\n     sage: 5 # optional guava\n')
            ['::\n\n',
            (None, '1 # optional guava\n', 'Integer(1) # optional guava\n'),
            '',
            '',
            (None, '3 # optional guava\n', 'Integer(3) # optional guava\n'),
            '\n::\n\n',
            (None, '4 # optional guava\n', 'Integer(4) # optional guava\n'),
            '',
            (None, '5 # optional guava\n', 'Integer(5) # optional guava\n'),
            '']

        TESTS::

            sage: parse("::\n\n    sage: # needs sage.combinat\n    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \\\n    ....:         import incidence_matrix_to_bit_rep_of_Vrep\n    sage: P = polytopes.associahedron(['A',3])\n\n")
            ['::\n\n',
            '',
            (None,
            'from sage.geometry.polyhedron.combinatorial_polyhedron.conversions import incidence_matrix_to_bit_rep_of_Vrep\n',
            'from sage.geometry.polyhedron.combinatorial_polyhedron.conversions import incidence_matrix_to_bit_rep_of_Vrep\n'),
            '',
            (None,
            "P = polytopes.associahedron(['A',3])\n",
            "P = polytopes.associahedron(['A',Integer(3)])\n"),
            '\n']

            sage: example4 = '::\n\n        sage: C.minimum_distance(algorithm="guava")  # optional - guava\n        ...\n        24\n\n'
            sage: parsed4 = DTP.parse(example4)
            sage: dte = parsed4[1]
            sage: dte.sage_source
            'C.minimum_distance(algorithm="guava")  # optional - guava\n'
            sage: dte.want
            '...\n24\n'
        """
        # Regular expressions
        find_sage_prompt = re.compile(r"^(\s*)sage: ", re.M)
        find_sage_continuation = re.compile(r"^(\s*)\.\.\.\.:", re.M)
        find_python_continuation = re.compile(r"^(\s*)\.\.\.([^\.])", re.M)
        python_prompt = re.compile(r"^(\s*)>>>", re.M)
        backslash_replacer = re.compile(r"""(\s*)sage:(.*)\\\ *
\ *((\.){4}:)?\ *""")

        # The following are used to allow ... at the beginning of output
        ellipsis_tag = "<TEMP_ELLIPSIS_TAG>"

        # Hack for non-standard backslash line escapes accepted by the current
        # doctest system.
        m = backslash_replacer.search(string)
        while m is not None:
            g = m.groups()
            string = string[:m.start()] + g[0] + "sage:" + g[1] + string[m.end():]
            m = backslash_replacer.search(string, m.start())

        replace_ellipsis = not python_prompt.search(string)
        if replace_ellipsis:
            # There are no >>> prompts, so we can allow ... to begin the output
            # We do so by replacing ellipses with a special tag, then putting them back after parsing
            string = find_python_continuation.sub(r"\1" + ellipsis_tag + r"\2", string)
        string = find_sage_prompt.sub(r"\1>>> sage: ", string)
        string = find_sage_continuation.sub(r"\1...", string)
        res = doctest.DocTestParser.parse(self, string, *args)
        filtered = []
        persistent_optional_tags = self.file_optional_tags
        persistent_optional_tag_setter = None
        persistent_optional_tag_setter_index = None
        first_example_in_block = None
        first_example_in_block_index = None
        tag_count_within_block = defaultdict(lambda: 0)

        def update_tag_counts(optional_tags):
            for tag in optional_tags:
                if tag not in persistent_optional_tags:
                    tag_count_within_block[tag] += 1
            tag_count_within_block[''] += 1

        def check_and_clear_tag_counts():
            if (num_examples := tag_count_within_block['']) >= 4:
                if overused_tags := set(tag for tag, count in tag_count_within_block.items()
                                        if tag and count >= num_examples):
                    overused_tags.update(persistent_optional_tags)
                    overused_tags.difference_update(self.file_optional_tags)
                    suggested = unparse_optional_tags(overused_tags, prefix='sage: # ')

                    if persistent_optional_tag_setter:
                        warning_example = persistent_optional_tag_setter
                        index = persistent_optional_tag_setter_index
                        warning = (f"Consider updating this block-scoped tag to '{suggested}' "
                                   f"to avoid repeating the tag {num_examples} times")
                    else:
                        warning_example = first_example_in_block
                        index = first_example_in_block_index
                        warning = (f"Consider using a block-scoped tag by "
                                   f"inserting the line '{suggested}' just before this line "
                                   f"to avoid repeating the tag {num_examples} times")

                    if not (index < len(filtered) and filtered[index] == warning_example):
                        # The example to which we want to attach our warning is
                        # not in ``filtered``. It is either the persistent tag line,
                        # or the first example of the block and not run because of unmet tags,
                        # or just a comment. Either way, we transform this example
                        # to a virtual example and attach the warning to it.
                        warning_example.sage_source = warning_example.source
                        if warning_example.sage_source.startswith("sage: "):
                            warning_example.sage_source = warning_example.source[6:]
                        warning_example.source = 'None  # virtual doctest'
                        warning_example.want = ''
                        filtered.insert(index, warning_example)

                    if not hasattr(warning_example, 'warnings'):
                        warning_example.warnings = []
                    warning_example.warnings.append(warning)
            tag_count_within_block.clear()

        for item in res:
            if isinstance(item, doctest.Example):
                optional_tags, source_sans_tags, is_persistent = parse_optional_tags(item.source, return_string_sans_tags=True)
                optional_tags = set(optional_tags)
                if is_persistent:
                    check_and_clear_tag_counts()
                    persistent_optional_tags = optional_tags
                    persistent_optional_tags.update(self.file_optional_tags)
                    persistent_optional_tag_setter = first_example_in_block = item
                    persistent_optional_tag_setter_index = len(filtered)
                    first_example_in_block_index = None
                    continue

                if not first_example_in_block:
                    first_example_in_block = item
                    first_example_in_block_index = len(filtered)
                update_tag_counts(optional_tags)
                optional_tags.update(persistent_optional_tags)
                item.optional_tags = frozenset(optional_tags)
                item.probed_tags = set()
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
                        extra = set(tag
                                    for tag in optional_tags
                                    if (tag not in self.optional_tags
                                        and tag not in available_software))
                        if extra:
                            if any(tag in external_software for tag in extra):
                                # never probe "external" software
                                continue
                            if all(tag in persistent_optional_tags for tag in extra):
                                # don't probe if test is only conditional
                                # on file-level or block-level tags
                                continue
                            if self.probed_tags is True:
                                item.probed_tags = extra
                            elif all(tag in self.probed_tags for tag in extra):
                                item.probed_tags = extra
                            else:
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
            else:
                if '\n' in item:
                    check_and_clear_tag_counts()
                    persistent_optional_tags = self.file_optional_tags
                    persistent_optional_tag_setter = first_example_in_block = None
                    persistent_optional_tag_setter_index = first_example_in_block_index = None
            filtered.append(item)

        check_and_clear_tag_counts()

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
                return '<ESC-' + ansi_escape[1] + '>'
            return '<CSI-' + ansi_escape.lstrip('\x1b[\x9b') + '>'
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
            sage: OC.add_tolerance(RIF(pi.n(64)), want_tol).endpoints()                 # needs sage.symbolic
            (3.14127849432443, 3.14190681285516)
            sage: OC.add_tolerance(RIF(pi.n(64)), want_abs).endpoints()                 # needs sage.symbolic
            (3.14149265358979, 3.14169265358980)
            sage: OC.add_tolerance(RIF(pi.n(64)), want_rel).endpoints()                 # needs sage.symbolic
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
                return RIFtol(want.tol) * RIFtol(-1, 1)
            else:
                return wantval * (1 + RIFtol(want.tol) * RIFtol(-1, 1))
        elif want.abs_tol:
            return wantval + RIFtol(want.abs_tol) * RIFtol(-1, 1)
        elif want.rel_tol:
            return wantval * (1 + RIFtol(want.rel_tol) * RIFtol(-1, 1))
        else:
            return wantval

    def check_output(self, want, got, optionflags):
        r"""
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

            sage: _ = x  # rel tol 1e10                                                 # needs sage.symbolic
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

        This illustrates that :trac:`33588` is fixed::

            sage: from sage.doctest.parsing import SageOutputChecker, SageDocTestParser
            sage: import doctest
            sage: optflag = doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS
            sage: DTP = SageDocTestParser(('sage','magma','guava'))
            sage: OC = SageOutputChecker()
            sage: example = "sage: 1.3090169943749475 # tol 1e-8\n1.3090169943749475"
            sage: ex = DTP.parse(example)[1]
            sage: OC.check_output(ex.want, '1.3090169943749475', optflag)
            True
            sage: OC.check_output(ex.want, 'ANYTHING1.3090169943749475', optflag)
            False
            sage: OC.check_output(ex.want, 'Long-step dual simplex will be used\n1.3090169943749475', optflag)
            True
        """
        # Regular expression for floats
        float_regex = re.compile(r'\s*([+-]?\s*((\d*\.?\d+)|(\d+\.?))([eE][+-]?\d+)?)')

        got = self.human_readable_escape_sequences(got)

        if isinstance(want, MarkedOutput):
            if want.random:
                return True
            elif want.tol or want.rel_tol or want.abs_tol:
                # First check that the number of occurrences of floats appearing match
                want_str = [g[0] for g in float_regex.findall(want)]
                got_str = [g[0] for g in float_regex.findall(got)]
                if len(want_str) != len(got_str):
                    return False

                # Then check the numbers
                want_values = [RIFtol(g) for g in want_str]
                want_intervals = [self.add_tolerance(v, want) for v in want_values]
                got_values = [RIFtol(g) for g in got_str]
                # The doctest is not successful if one of the "want" and "got"
                # intervals have an empty intersection
                if not all(a.overlaps(b) for a, b in zip(want_intervals, got_values)):
                    return False

                # Then check the part of the doctests without the numbers
                # Continue the check process with floats replaced by stars
                want = float_regex.sub('*', want)
                got = float_regex.sub('*', got)

        if doctest.OutputChecker.check_output(self, want, got, optionflags):
            return True
        else:
            # Last resort: try to fix-up the got string removing few typical warnings
            did_fixup, want, got = self.do_fixup(want, got)
            if did_fixup:
                return doctest.OutputChecker.check_output(self, want, got, optionflags)
            else:
                return False

    def do_fixup(self, want, got):
        r"""
        Performs few changes to the strings ``want`` and ``got``.

        For example, remove warnings to be ignored.

        INPUT:

        - ``want`` -- a string or :class:`MarkedOutput`
        - ``got`` -- a string

        OUTPUT:

        A tuple:

        - bool, ``True`` when some fixup were performed and ``False`` otherwise
        - string, edited wanted string
        - string, edited got string

        .. NOTE::

            Currently, the code only possibly changes the string ``got``
            while keeping ``want`` invariant. We keep open the possibility
            of adding a regular expression which would also change the
            ``want`` string. This is why ``want`` is an input and an output
            of the method even if currently kept invariant.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: OC = SageOutputChecker()
            sage: OC.do_fixup('1.3090169943749475','1.3090169943749475')
            (False, '1.3090169943749475', '1.3090169943749475')
            sage: OC.do_fixup('1.3090169943749475','ANYTHING1.3090169943749475')
            (False, '1.3090169943749475', 'ANYTHING1.3090169943749475')
            sage: OC.do_fixup('1.3090169943749475','Long-step dual simplex will be used\n1.3090169943749475')
            (True, '1.3090169943749475', '\n1.3090169943749475')

        When ``want`` is an instance of class :class:`MarkedOutput`::

            sage: from sage.doctest.parsing import SageOutputChecker, SageDocTestParser
            sage: import doctest
            sage: optflag = doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS
            sage: DTP = SageDocTestParser(('sage','magma','guava'))
            sage: OC = SageOutputChecker()
            sage: example = "sage: 1.3090169943749475\n1.3090169943749475"
            sage: ex = DTP.parse(example)[1]
            sage: ex.want
            '1.3090169943749475\n'
            sage: OC.do_fixup(ex.want,'1.3090169943749475')
            (False, '1.3090169943749475\n', '1.3090169943749475')
            sage: OC.do_fixup(ex.want,'ANYTHING1.3090169943749475')
            (False, '1.3090169943749475\n', 'ANYTHING1.3090169943749475')
            sage: OC.do_fixup(ex.want,'Long-step dual simplex will be used\n1.3090169943749475')
            (True, '1.3090169943749475\n', '\n1.3090169943749475')

        """
        did_fixup = False

        # The conditions in the below `if` are simple fast test on the expected
        # and/or actual output to determine if a fixup should be applied.

        if "Long-step" in got:
            # Version 4.65 of glpk prints the warning "Long-step dual
            # simplex will be used" frequently. When Sage uses a system
            # installation of glpk which has not been patched, we need to
            # ignore that message. See :trac:`29317`.
            glpk_simplex_warning_regex = re.compile(r'(Long-step dual simplex will be used)')
            got = glpk_simplex_warning_regex.sub('', got)
            did_fixup = True

        if "chained fixups" in got:
            # :trac:`34533` -- suppress warning on OS X 12.6 about chained fixups
            chained_fixup_warning_regex = re.compile(r'ld: warning: -undefined dynamic_lookup may not work with chained fixups')
            got = chained_fixup_warning_regex.sub('', got)
            did_fixup = True

        if "newer macOS version" in got:
            # :trac:`34741` -- suppress warning arising after
            # upgrading from macOS 12.X to 13.X.
            newer_macOS_version_regex = re.compile(r'.*dylib \(.*\) was built for newer macOS version \(.*\) than being linked \(.*\)')
            got = newer_macOS_version_regex.sub('', got)
            did_fixup = True

        if "insufficient permissions" in got:
            sympow_cache_warning_regex = re.compile(r'\*\*WARNING\*\* /var/cache/sympow/datafiles/le64 yields insufficient permissions')
            got = sympow_cache_warning_regex.sub('', got)
            did_fixup = True

        if "dylib" in got:
            # :trac:`31204` -- suppress warning about ld and OS version for
            # dylib files.
            ld_warning_regex = re.compile(r'^.*dylib.*was built for newer macOS version.*than being linked.*')
            got = ld_warning_regex.sub('', got)
            did_fixup = True

        if "pie being ignored" in got:
            # :trac:`30845` -- suppress warning on conda about ld
            ld_pie_warning_regex = re.compile(r'ld: warning: -pie being ignored. It is only used when linking a main executable')
            got = ld_pie_warning_regex.sub('', got)
            did_fixup = True

        return did_fixup, want, got

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
        # Regular expression for floats
        float_regex = re.compile(r'\s*([+-]?\s*((\d*\.?\d+)|(\d+\.?))([eE][+-]?\d+)?)')

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
                                fail(wstr, gstr, abs(1 - g / w), want.tol)
                        elif want.abs_tol:
                            fail(wstr, gstr, abs(g - w), want.abs_tol)
                        else:
                            fail(wstr, gstr, abs(1 - g / w), want.rel_tol)

                if failures:
                    if len(want_str) == 1:
                        diff += "Tolerance exceeded:\n"
                    else:
                        diff += "Tolerance exceeded in %s of %s:\n" % (len(failures), len(want_str))
                    diff += "\n".join(failures) + "\n"
            elif "..." in want:
                diff += "Note: combining tolerance (# tol) with ellipsis (...) is not supported\n"
        return diff
