r"""
Asymptotic Expansions --- Miscellaneous

AUTHORS:

- Daniel Krenn (2015)

ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

- Benjamin Hackl is supported by the Google Summer of Code 2015.


Functions, Classes and Methods
==============================
"""

#*****************************************************************************
# Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************


from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject


def repr_short_to_parent(s):
    r"""
    Helper method for the growth group factory, which converts a short
    representation string to a parent.

    INPUT:

    - ``s`` -- a string, short representation of a parent.

    OUTPUT:

    A parent.

    The possible short representations are shown in the examples below.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import repr_short_to_parent
        sage: repr_short_to_parent('ZZ')
        Integer Ring
        sage: repr_short_to_parent('QQ')
        Rational Field
        sage: repr_short_to_parent('SR')
        Symbolic Ring
        sage: repr_short_to_parent('NN')
        Non negative integer semiring
        sage: repr_short_to_parent('UU')
        Group of Roots of Unity

    TESTS::

        sage: repr_short_to_parent('abcdef')
        Traceback (most recent call last):
        ...
        ValueError: Cannot create a parent out of 'abcdef'.
        > *previous* ValueError: unknown specification abcdef
        > *and* NameError: name 'abcdef' is not defined
    """
    from sage.groups.misc_gps.argument_groups import ArgumentGroup
    from sage.misc.sage_eval import sage_eval

    def extract(s):
        try:
            return ArgumentGroup(specification=s)
        except Exception as e:
            e_ag = e
            e_ag.__traceback__ = None

        try:
            return sage_eval(s)
        except Exception as e:
            e_se = e
            e_se.__traceback__ = None

        raise combine_exceptions(
            ValueError("Cannot create a parent out of '%s'." % (s,)),
            e_ag, e_se)

    P = extract(s)

    from sage.misc.lazy_import import LazyImport
    if type(P) is LazyImport:
        P = P._get_object()

    from sage.structure.parent import is_Parent
    if not is_Parent(P):
        raise ValueError("'%s' does not describe a parent." % (s,))
    return P


def parent_to_repr_short(P):
    r"""
    Helper method which generates a short(er) representation string
    out of a parent.

    INPUT:

    - ``P`` -- a parent.

    OUTPUT:

    A string.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import parent_to_repr_short
        sage: parent_to_repr_short(ZZ)
        'ZZ'
        sage: parent_to_repr_short(QQ)
        'QQ'
        sage: parent_to_repr_short(SR)
        'SR'
        sage: parent_to_repr_short(RR)
        'RR'
        sage: parent_to_repr_short(CC)
        'CC'
        sage: parent_to_repr_short(ZZ['x'])
        'ZZ[x]'
        sage: parent_to_repr_short(QQ['d, k'])
        'QQ[d, k]'
        sage: parent_to_repr_short(QQ['e'])
        'QQ[e]'
        sage: parent_to_repr_short(SR[['a, r']])
        'SR[[a, r]]'
        sage: parent_to_repr_short(Zmod(3))
        'Ring of integers modulo 3'
        sage: parent_to_repr_short(Zmod(3)['g'])
        'Univariate Polynomial Ring in g over Ring of integers modulo 3'
    """
    from sage.rings.all import RR, CC, RIF, CIF, RBF, CBF
    from sage.rings.integer_ring import ZZ
    from sage.rings.rational_field import QQ
    from sage.symbolic.ring import SR
    from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
    from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing
    from sage.rings.power_series_ring import is_PowerSeriesRing
    def abbreviate(P):
        try:
            return P._repr_short_()
        except AttributeError:
            pass
        abbreviations = {ZZ: 'ZZ', QQ: 'QQ', SR: 'SR',
                         RR: 'RR', CC: 'CC',
                         RIF: 'RIF', CIF: 'CIF',
                         RBF: 'RBF', CBF: 'CBF'}
        try:
            return abbreviations[P]
        except KeyError:
            pass
        raise ValueError('Cannot abbreviate %s.' % (P,))

    poly = is_PolynomialRing(P) or is_MPolynomialRing(P)
    from sage.rings import multi_power_series_ring
    power = is_PowerSeriesRing(P) or \
            multi_power_series_ring.is_MPowerSeriesRing(P)

    if poly or power:
        if poly:
            op, cl = ('[', ']')
        else:
            op, cl = ('[[', ']]')
        try:
            s = abbreviate(P.base_ring()) + op + ', '.join(P.variable_names()) + cl
        except ValueError:
            s = str(P)
    else:
        try:
            s = abbreviate(P)
        except ValueError:
            s = str(P)

    return s


def split_str_by_op(string, op, strip_parentheses=True):
    r"""
    Split the given string into a tuple of substrings arising by
    splitting by ``op`` and taking care of parentheses.

    INPUT:

    - ``string`` -- a string.

    - ``op`` -- a string. This is used by
      :python:`str.split <library/stdtypes.html#str.split>`.
      Thus, if this is ``None``, then any whitespace string is a
      separator and empty strings are removed from the result.

    - ``strip_parentheses`` -- (default: ``True``) a boolean.

    OUTPUT:

    A tuple of strings.

    TESTS::

        sage: from sage.rings.asymptotic.misc import split_str_by_op
        sage: split_str_by_op('x^ZZ', '*')
        ('x^ZZ',)
        sage: split_str_by_op('log(x)^ZZ * y^QQ', '*')
        ('log(x)^ZZ', 'y^QQ')
        sage: split_str_by_op('log(x)**ZZ * y**QQ', '*')
        ('log(x)**ZZ', 'y**QQ')
        sage: split_str_by_op('a^b * * c^d', '*')
        Traceback (most recent call last):
        ...
        ValueError: 'a^b * * c^d' is invalid since a '*' follows a '*'.
        sage: split_str_by_op('a^b * (c*d^e)', '*')
        ('a^b', 'c*d^e')

    ::

        sage: split_str_by_op('(a^b)^c', '^')
        ('a^b', 'c')
        sage: split_str_by_op('a^(b^c)', '^')
        ('a', 'b^c')

    ::

        sage: split_str_by_op('(a) + (b)', op='+', strip_parentheses=True)
        ('a', 'b')
        sage: split_str_by_op('(a) + (b)', op='+', strip_parentheses=False)
        ('(a)', '(b)')
        sage: split_str_by_op(' ( t  ) ', op='+', strip_parentheses=False)
        ('( t  )',)

    ::

        sage: split_str_by_op(' ( t  ) ', op=None)
        ('t',)
        sage: split_str_by_op(' ( t  )s', op=None)
        ('(t)s',)
        sage: split_str_by_op(' ( t  ) s', op=None)
        ('t', 's')

    ::

        sage: split_str_by_op('(e^(n*log(n)))^SR.subring(no_variables=True)', '*')
        ('(e^(n*log(n)))^SR.subring(no_variables=True)',)
    """
    def is_balanced(s):
        open = 0
        for l in s:
            if l == '(':
                open += 1
            elif l == ')':
                open -= 1
            if open < 0:
                return False
        return bool(open == 0)

    factors = list()
    balanced = True
    if string and op is not None and string.startswith(op):
        raise ValueError("'%s' is invalid since it starts with a '%s'." %
                         (string, op))
    for s in string.split(op):
        if not s:
            factors[-1] += op
            balanced = False
            continue
        if not s.strip():
            raise ValueError("'%s' is invalid since a '%s' follows a '%s'." %
                             (string, op, op))
        if not balanced:
            s = factors.pop() + (op if op else '') + s
        balanced = is_balanced(s)
        factors.append(s)

    if not balanced:
        raise ValueError("Parentheses in '%s' are not balanced." % (string,))

    def strip(s):
        s = s.strip()
        if not s:
            return s
        if strip_parentheses and s[0] == '(' and s[-1] == ')':
            t = s[1:-1]
            if is_balanced(t):
                s = t
        return s.strip()

    return tuple(strip(f) for f in factors)


def repr_op(left, op, right=None, latex=False):
    r"""
    Create a string ``left op right`` with
    taking care of parentheses in its operands.

    INPUT:

    - ``left`` -- an element.

    - ``op`` -- a string.

    - ``right`` -- an element.

    - ``latex`` -- (default: ``False``) a boolean. If set, then
      LaTeX-output is returned.

    OUTPUT:

    A string.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import repr_op
        sage: repr_op('a^b', '^', 'c')
        '(a^b)^c'

    TESTS::

        sage: repr_op('a-b', '^', 'c')
        '(a-b)^c'
        sage: repr_op('a+b', '^', 'c')
        '(a+b)^c'

    ::

        sage: print(repr_op(r'\frac{1}{2}', '^', 'c', latex=True))
        \left(\frac{1}{2}\right)^c

    ::

        sage: repr_op('Arg', '_', 'Symbolic Ring')
        'Arg_(Symbolic Ring)'
    """
    left = str(left)
    right = str(right) if right is not None else ''

    def add_parentheses(s, op):
        if op in ('^', '_'):
            signals = ('^', '/', '*', '-', '+', ' ')
        else:
            return s
        if any(sig in s for sig in signals) or latex and s.startswith(r'\frac'):
            if latex:
                return r'\left({}\right)'.format(s)
            else:
                return '({})'.format(s)
        else:
            return s

    return add_parentheses(left, op) + op + add_parentheses(right, op)


def combine_exceptions(e, *f):
    r"""
    Helper function which combines the messages of the given exceptions.

    INPUT:

    - ``e`` -- an exception.

    - ``*f`` -- exceptions.

    OUTPUT:

    An exception.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import combine_exceptions
        sage: raise combine_exceptions(ValueError('Outer.'), TypeError('Inner.'))
        Traceback (most recent call last):
        ...
        ValueError: Outer.
        > *previous* TypeError: Inner.
        sage: raise combine_exceptions(ValueError('Outer.'),
        ....:                          TypeError('Inner1.'), TypeError('Inner2.'))
        Traceback (most recent call last):
        ...
        ValueError: Outer.
        > *previous* TypeError: Inner1.
        > *and* TypeError: Inner2.
        sage: raise combine_exceptions(ValueError('Outer.'),
        ....:                          combine_exceptions(TypeError('Middle.'),
        ....:                                             TypeError('Inner.')))
        Traceback (most recent call last):
        ...
        ValueError: Outer.
        > *previous* TypeError: Middle.
        >> *previous* TypeError: Inner.
    """
    import re
    msg = ('\n *previous* ' +
           '\n *and* '.join("%s: %s" % (ff.__class__.__name__, str(ff)) for ff in f))
    msg = re.sub(r'^([>]* \*previous\*)', r'>\1', msg, flags=re.MULTILINE)
    msg = re.sub(r'^([>]* \*and\*)', r'>\1', msg, flags=re.MULTILINE)
    msg = str(e.args if len(e.args) > 1 else e.args[0]) + msg
    e.args = (msg,)
    return e


def substitute_raise_exception(element, e):
    r"""
    Raise an error describing what went wrong with the substitution.

    INPUT:

    - ``element`` -- an element.

    - ``e`` -- an exception which is included in the raised error
      message.

    OUTPUT:

    Raise an exception of the same type as ``e``.

    TESTS::

        sage: from sage.rings.asymptotic.misc import substitute_raise_exception
        sage: substitute_raise_exception(x, Exception('blub'))
        Traceback (most recent call last):
        ...
        Exception: Cannot substitute in x in Symbolic Ring.
        > *previous* Exception: blub
    """
    raise combine_exceptions(
        type(e)('Cannot substitute in %s in %s.' %
                (element, element.parent())), e)


def bidirectional_merge_overlapping(A, B, key=None):
    r"""
    Merge the two overlapping tuples/lists.

    INPUT:

    - ``A`` -- a list or tuple (type has to coincide with type of ``B``).

    - ``B`` -- a list or tuple (type has to coincide with type of ``A``).

    - ``key`` -- (default: ``None``) a function. If ``None``, then the
      identity is used.  This ``key``-function applied on an element
      of the list/tuple is used for comparison. Thus elements with the
      same key are considered as equal.

    OUTPUT:

    A pair of lists or tuples (depending on the type of ``A`` and ``B``).

    .. NOTE::

        Suppose we can decompose the list `A=ac` and `B=cb` with
        lists `a`, `b`, `c`, where `c` is nonempty. Then
        :func:`bidirectional_merge_overlapping` returns the pair `(acb, acb)`.

        Suppose a ``key``-function is specified and `A=ac_A` and
        `B=c_Bb`, where the list of keys of the elements of `c_A`
        equals the list of keys of the elements of `c_B`. Then
        :func:`bidirectional_merge_overlapping` returns the pair `(ac_Ab, ac_Bb)`.

        After unsuccessfully merging `A=ac` and `B=cb`,
        a merge of `A=ca` and `B=bc` is tried.

    TESTS::

        sage: from sage.rings.asymptotic.misc import bidirectional_merge_overlapping
        sage: def f(L, s):
        ....:     return list((ell, s) for ell in L)
        sage: key = lambda k: k[0]
        sage: bidirectional_merge_overlapping(f([0..3], 'a'), f([5..7], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: bidirectional_merge_overlapping(f([0..2], 'a'), f([4..7], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: bidirectional_merge_overlapping(f([4..7], 'a'), f([0..2], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: bidirectional_merge_overlapping(f([0..3], 'a'), f([3..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'b')],
         [(0, 'a'), (1, 'a'), (2, 'a'), (3, 'b'), (4, 'b')])
        sage: bidirectional_merge_overlapping(f([3..4], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'b'), (2, 'b'), (3, 'a'), (4, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b'), (4, 'a')])
        sage: bidirectional_merge_overlapping(f([0..1], 'a'), f([0..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'b'), (3, 'b'), (4, 'b')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b'), (4, 'b')])
        sage: bidirectional_merge_overlapping(f([0..3], 'a'), f([0..1], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'a'), (3, 'a')])
        sage: bidirectional_merge_overlapping(f([0..3], 'a'), f([1..3], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'a'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: bidirectional_merge_overlapping(f([1..3], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: bidirectional_merge_overlapping(f([0..6], 'a'), f([3..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'a'), (5, 'a'), (6, 'a')],
         [(0, 'a'), (1, 'a'), (2, 'a'), (3, 'b'), (4, 'b'), (5, 'a'), (6, 'a')])
        sage: bidirectional_merge_overlapping(f([0..3], 'a'), f([1..2], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'a'), (1, 'b'), (2, 'b'), (3, 'a')])
        sage: bidirectional_merge_overlapping(f([1..2], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'a'), (2, 'a'), (3, 'b')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: bidirectional_merge_overlapping(f([1..3], 'a'), f([1..3], 'b'), key)
        ([(1, 'a'), (2, 'a'), (3, 'a')],
         [(1, 'b'), (2, 'b'), (3, 'b')])
    """
    if key is None:
        Akeys = A
        Bkeys = B
    else:
        Akeys = tuple(key(a) for a in A)
        Bkeys = tuple(key(b) for b in B)

    def find_overlapping_index(A, B):
        if len(B) > len(A) - 2:
            raise StopIteration
        matches = iter(i for i in range(1, len(A) - len(B))
                       if A[i:i+len(B)] == B)
        return next(matches)

    def find_mergedoverlapping_index(A, B):
        """
        Return in index i where to merge two overlapping tuples/lists ``A`` and ``B``.

        Then ``A + B[i:]`` or ``A[:-i] + B`` are the merged tuples/lists.

        Adapted from https://stackoverflow.com/a/30056066/1052778.
        """
        matches = iter(i for i in range(min(len(A), len(B)), 0, -1)
                       if A[-i:] == B[:i])
        return next(matches, 0)

    i = find_mergedoverlapping_index(Akeys, Bkeys)
    if i > 0:
        return A + B[i:], A[:-i] + B

    i = find_mergedoverlapping_index(Bkeys, Akeys)
    if i > 0:
        return B[:-i] + A, B + A[i:]

    try:
        i = find_overlapping_index(Akeys, Bkeys)
    except StopIteration:
        pass
    else:
        return A, A[:i] + B + A[i+len(B):]

    try:
        i = find_overlapping_index(Bkeys, Akeys)
    except StopIteration:
        pass
    else:
        return B[:i] + A + B[i+len(A):], B

    raise ValueError('Input does not have an overlap.')


def bidirectional_merge_sorted(A, B, key=None):
    r"""
    Merge the two tuples/lists, keeping the orders provided by them.

    INPUT:

    - ``A`` -- a list or tuple (type has to coincide with type of ``B``).

    - ``B`` -- a list or tuple (type has to coincide with type of ``A``).

    - ``key`` -- (default: ``None``) a function. If ``None``, then the
      identity is used.  This ``key``-function applied on an element
      of the list/tuple is used for comparison. Thus elements with the
      same key are considered as equal.

    .. NOTE::

        The two tuples/list need to overlap, i.e. need at least
        one key in common.

    OUTPUT:

    A pair of lists containing all elements totally ordered. (The first
    component uses ``A`` as a merge base, the second component ``B``.)

    If merging fails, then a
    :python:`RuntimeError<library/exceptions.html#exceptions.RuntimeError>`
    is raised.

    TESTS::

        sage: from sage.rings.asymptotic.misc import bidirectional_merge_sorted
        sage: def f(L, s):
        ....:     return list((ell, s) for ell in L)
        sage: key = lambda k: k[0]
        sage: bidirectional_merge_sorted(f([0..3], 'a'), f([5..7], 'b'), key)
        Traceback (most recent call last):
        ...
        RuntimeError: no common elements
        sage: bidirectional_merge_sorted(f([0..2], 'a'), f([4..7], 'b'), key)
        Traceback (most recent call last):
        ...
        RuntimeError: no common elements
        sage: bidirectional_merge_sorted(f([4..7], 'a'), f([0..2], 'b'), key)
        Traceback (most recent call last):
        ...
        RuntimeError: no common elements
        sage: bidirectional_merge_sorted(f([0..3], 'a'), f([3..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'b')],
         [(0, 'a'), (1, 'a'), (2, 'a'), (3, 'b'), (4, 'b')])
        sage: bidirectional_merge_sorted(f([3..4], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'b'), (2, 'b'), (3, 'a'), (4, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b'), (4, 'a')])
        sage: bidirectional_merge_sorted(f([0..1], 'a'), f([0..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'b'), (3, 'b'), (4, 'b')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b'), (4, 'b')])
        sage: bidirectional_merge_sorted(f([0..3], 'a'), f([0..1], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'a'), (3, 'a')])
        sage: bidirectional_merge_sorted(f([0..3], 'a'), f([1..3], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'a'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: bidirectional_merge_sorted(f([1..3], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: bidirectional_merge_sorted(f([0..6], 'a'), f([3..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'a'), (5, 'a'), (6, 'a')],
         [(0, 'a'), (1, 'a'), (2, 'a'), (3, 'b'), (4, 'b'), (5, 'a'), (6, 'a')])
        sage: bidirectional_merge_sorted(f([0..3], 'a'), f([1..2], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'a'), (1, 'b'), (2, 'b'), (3, 'a')])
        sage: bidirectional_merge_sorted(f([1..2], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'a'), (2, 'a'), (3, 'b')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: bidirectional_merge_sorted(f([1..3], 'a'), f([1..3], 'b'), key)
        ([(1, 'a'), (2, 'a'), (3, 'a')],
         [(1, 'b'), (2, 'b'), (3, 'b')])

    ::

        sage: bidirectional_merge_sorted(f([1, 2, 3], 'a'), f([1, 3], 'b'), key)
        ([(1, 'a'), (2, 'a'), (3, 'a')],
         [(1, 'b'), (2, 'a'), (3, 'b')])
        sage: bidirectional_merge_sorted(f([1, 4, 5, 6], 'a'), f([1, 2, 3, 4, 6], 'b'), key)
        ([(1, 'a'), (2, 'b'), (3, 'b'), (4, 'a'), (5, 'a'), (6, 'a')],
         [(1, 'b'), (2, 'b'), (3, 'b'), (4, 'b'), (5, 'a'), (6, 'b')])
        sage: bidirectional_merge_sorted(f([1, 2, 3, 4], 'a'), f([1, 3, 5], 'b'), key)
        Traceback (most recent call last):
        ...
        RuntimeError: sorting not unique
        sage: bidirectional_merge_sorted(f([1, 2], 'a'), f([2, 1], 'b'), key)
        Traceback (most recent call last):
        ...
        RuntimeError: sorting in lists not compatible
        sage: bidirectional_merge_sorted(f([1, 2, 4], 'a'), f([1, 3, 4], 'b'), key)
        Traceback (most recent call last):
        ...
        RuntimeError: sorting not unique
    """
    if key is None:
        Akeys = A
        Bkeys = B
    else:
        Akeys = tuple(key(a) for a in A)
        Bkeys = tuple(key(b) for b in B)

    matches = tuple((i, j)
                    for i, a in enumerate(Akeys)
                    for j, b in enumerate(Bkeys)
                    if a == b)
    if not matches:
        raise RuntimeError('no common elements')

    resultA = []
    resultB = []
    last = (0, 0)
    end = (len(A), len(B))
    for current in matches + (end,):
        if not all(a <= b for a, b in zip(last, current)):
            raise RuntimeError('sorting in lists not compatible')
        if last[0] == current[0]:
            resultA.extend(B[last[1]:current[1]])
            resultB.extend(B[last[1]:current[1]])
        elif last[1] == current[1]:
            resultA.extend(A[last[0]:current[0]])
            resultB.extend(A[last[0]:current[0]])
        else:
            raise RuntimeError('sorting not unique')
        if current != end:
            resultA.append(A[current[0]])
            resultB.append(B[current[1]])
            last = (current[0]+1, current[1]+1)

    return (resultA, resultB)


def log_string(element, base=None):
    r"""
    Return a representation of the log of the given element to the
    given base.

    INPUT:

    - ``element`` -- an object.

    - ``base`` -- an object or ``None``.

    OUTPUT:

    A string.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import log_string
        sage: log_string(3)
        'log(3)'
        sage: log_string(3, base=42)
        'log(3, base=42)'
    """
    basestr = ', base=' + str(base) if base else ''
    return 'log(%s%s)' % (element, basestr)


def strip_symbolic(expression):
    r"""
    Return, if possible, the underlying (numeric) object of
    the symbolic expression.

    If ``expression`` is not symbolic, then ``expression`` is returned.

    INPUT:

    - ``expression`` -- an object

    OUTPUT:

    An object.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import strip_symbolic
        sage: strip_symbolic(SR(2)); _.parent()
        2
        Integer Ring
        sage: strip_symbolic(SR(2/3)); _.parent()
        2/3
        Rational Field
        sage: strip_symbolic(SR('x')); _.parent()
        x
        Symbolic Ring
        sage: strip_symbolic(pi); _.parent()
        pi
        Symbolic Ring
    """
    from sage.structure.element import parent, Element
    from sage.symbolic.ring import SymbolicRing

    P = parent(expression)
    if isinstance(P, SymbolicRing):
        try:
            stripped = expression.pyobject()
            if isinstance(stripped, Element):
                return stripped
        except TypeError:
            pass
    return expression


class NotImplementedOZero(NotImplementedError):
    r"""
    A special :python:`NotImplementedError<library/exceptions.html#exceptions.NotImplementedError>`
    which is raised when the result is O(0) which means 0
    for sufficiently large values of the variable.
    """
    def __init__(self, asymptotic_ring=None, var=None, exact_part=0):
        r"""
        INPUT:

        - ``asymptotic_ring`` -- (default: ``None``) an :class:`AsymptoticRing` or ``None``.

        - ``var`` -- (default: ``None``) a string.

        Either ``asymptotic_ring`` or ``var`` has to be specified.

        - ``exact_part`` -- (default: ``0``) asymptotic expansion

        EXAMPLES::

            sage: A = AsymptoticRing('n^ZZ', ZZ)
            sage: from sage.rings.asymptotic.misc import NotImplementedOZero

            sage: raise NotImplementedOZero(A)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: got O(0)
            The error term O(0) means 0 for sufficiently large n.

            sage: raise NotImplementedOZero(var='m')
            Traceback (most recent call last):
            ...
            NotImplementedOZero: got O(0)
            The error term O(0) means 0 for sufficiently large m.

        TESTS::

            sage: raise NotImplementedOZero(A, var='m')
            Traceback (most recent call last):
            ...
            ValueError: specify either 'asymptotic_ring' or 'var'
            sage: raise NotImplementedOZero()
            Traceback (most recent call last):
            ...
            ValueError: specify either 'asymptotic_ring' or 'var'
        """
        if (asymptotic_ring is None) == (var is None):
            raise ValueError("specify either 'asymptotic_ring' or 'var'")

        if var is None:
            var = ', '.join(str(g) for g in asymptotic_ring.gens())
        message = ('got {}\n'.format(('{} + '.format(exact_part) if exact_part else '')
                                     + 'O(0)') +
                   'The error term O(0) '
                   'means 0 for sufficiently large {}.'.format(var))

        if asymptotic_ring is not None and isinstance(exact_part, int) and exact_part == 0:
            exact_part = asymptotic_ring.zero()
        self.exact_part = exact_part

        super(NotImplementedOZero, self).__init__(message)


class NotImplementedBZero(NotImplementedError):
    r"""
    A special :python:`NotImplementedError<library/exceptions.html#exceptions.NotImplementedError>`
    which is raised when the result is B(0) which means 0
    for sufficiently large values of the variable.
    """
    def __init__(self, asymptotic_ring=None, var=None, exact_part=0):
        r"""
        INPUT:

        - ``asymptotic_ring`` -- (default: ``None``) an :class:`AsymptoticRing` or ``None``.

        - ``var`` -- (default: ``None``) a string.

        Either ``asymptotic_ring`` or ``var`` has to be specified.

        - ``exact_part`` -- (default: ``0``) asymptotic expansion

        EXAMPLES::

            sage: A = AsymptoticRing('n^ZZ', ZZ)
            sage: from sage.rings.asymptotic.misc import NotImplementedBZero

            sage: raise NotImplementedBZero(A)
            Traceback (most recent call last):
            ...
            NotImplementedBZero: got B(0)
            The error term B(0) means 0 for sufficiently large n.

            sage: raise NotImplementedBZero(var='m')
            Traceback (most recent call last):
            ...
            NotImplementedBZero: got B(0)
            The error term B(0) means 0 for sufficiently large m.

            sage: AR.<n> = AsymptoticRing('n^QQ', QQ)
            sage: AR(0).B(42)
            Traceback (most recent call last):
            ...
            NotImplementedBZero: got B(0)
            The error term B(0) means 0 for sufficiently large n.

        TESTS::

            sage: raise NotImplementedBZero(A, var='m')
            Traceback (most recent call last):
            ...
            ValueError: specify either 'asymptotic_ring' or 'var'
            sage: raise NotImplementedBZero()
            Traceback (most recent call last):
            ...
            ValueError: specify either 'asymptotic_ring' or 'var'
        """
        if (asymptotic_ring is None) == (var is None):
            raise ValueError("specify either 'asymptotic_ring' or 'var'")

        if var is None:
            var = ', '.join(str(g) for g in asymptotic_ring.gens())
        message = ('got {}\n'.format(('{} + '.format(exact_part) if exact_part else '')
                                     + 'B(0)') +
                   'The error term B(0) '
                   'means 0 for sufficiently large {}.'.format(var))

        if asymptotic_ring is not None and isinstance(exact_part, int) and exact_part == 0:
            exact_part = asymptotic_ring.zero()
        self.exact_part = exact_part

        super(NotImplementedBZero, self).__init__(message)


def transform_category(category,
                       subcategory_mapping, axiom_mapping,
                       initial_category=None):
    r"""
    Transform ``category`` to a new category according to the given
    mappings.

    INPUT:

    - ``category`` -- a category.

    - ``subcategory_mapping`` -- a list (or other iterable) of triples
      ``(from, to, mandatory)``, where

      - ``from`` and ``to`` are categories and
      - ``mandatory`` is a boolean.

    - ``axiom_mapping`` -- a list (or other iterable) of triples
      ``(from, to, mandatory)``, where

      - ``from`` and ``to`` are strings describing axioms and
      - ``mandatory`` is a boolean.

    - ``initial_category`` -- (default: ``None``) a category. When
      transforming the given category, this ``initial_category`` is
      used as a starting point of the result. This means the resulting
      category will be a subcategory of ``initial_category``.
      If ``initial_category`` is ``None``, then the
      :class:`category of objects <sage.categories.objects.Objects>`
      is used.

    OUTPUT:

    A category.

    .. NOTE::

        Consider a subcategory mapping ``(from, to, mandatory)``. If
        ``category`` is a subcategory of ``from``, then the
        returned category will be a subcategory of ``to``. Otherwise and
        if ``mandatory`` is set, then an error is raised.

        Consider an axiom mapping ``(from, to, mandatory)``. If
        ``category`` is has axiom ``from``, then the
        returned category will have axiom ``to``. Otherwise and
        if ``mandatory`` is set, then an error is raised.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import transform_category
        sage: from sage.categories.additive_semigroups import AdditiveSemigroups
        sage: from sage.categories.additive_monoids import AdditiveMonoids
        sage: from sage.categories.additive_groups import AdditiveGroups
        sage: S = [
        ....:     (Sets(), Sets(), True),
        ....:     (Posets(), Posets(), False),
        ....:     (AdditiveMagmas(), Magmas(), False)]
        sage: A = [
        ....:     ('AdditiveAssociative', 'Associative', False),
        ....:     ('AdditiveUnital', 'Unital', False),
        ....:     ('AdditiveInverse', 'Inverse', False),
        ....:     ('AdditiveCommutative', 'Commutative', False)]
        sage: transform_category(Objects(), S, A)
        Traceback (most recent call last):
        ...
        ValueError: Category of objects is not
        a subcategory of Category of sets.
        sage: transform_category(Sets(), S, A)
        Category of sets
        sage: transform_category(Posets(), S, A)
        Category of posets
        sage: transform_category(AdditiveSemigroups(), S, A)
        Category of semigroups
        sage: transform_category(AdditiveMonoids(), S, A)
        Category of monoids
        sage: transform_category(AdditiveGroups(), S, A)
        Category of groups
        sage: transform_category(AdditiveGroups().AdditiveCommutative(), S, A)
        Category of commutative groups

    ::

        sage: transform_category(AdditiveGroups().AdditiveCommutative(), S, A,
        ....:     initial_category=Posets())
        Join of Category of commutative groups
            and Category of posets

    ::

        sage: transform_category(ZZ.category(), S, A)
        Category of commutative groups
        sage: transform_category(QQ.category(), S, A)
        Category of commutative groups
        sage: transform_category(SR.category(), S, A)
        Category of commutative groups
        sage: transform_category(Fields(), S, A)
        Category of commutative groups
        sage: transform_category(ZZ['t'].category(), S, A)
        Category of commutative groups

    ::

        sage: A[-1] = ('Commutative', 'AdditiveCommutative', True)
        sage: transform_category(Groups(), S, A)
        Traceback (most recent call last):
        ...
        ValueError: Category of groups does not have
        axiom Commutative.
    """
    if initial_category is None:
        from sage.categories.objects import Objects
        result = Objects()
    else:
        result = initial_category

    for A, B, mandatory in subcategory_mapping:
        if category.is_subcategory(A):
            result &= B
        elif mandatory:
            raise ValueError('%s is not a subcategory of %s.' %
                             (category, A))

    axioms = category.axioms()
    for A, B, mandatory in axiom_mapping:
        if A in axioms:
            result = result._with_axiom(B)
        elif mandatory:
            raise ValueError('%s does not have axiom %s.' %
                             (category, A))

    return result


class Locals(dict):
    r"""
    A frozen dictionary-like class for storing locals
    of an :class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticRing`.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import Locals
        sage: locals = Locals({'a': 42})
        sage: locals['a']
        42

    The object contains default values (see :meth:`default_locals`)
    for some keys::

        sage: locals['log']
        <function log at 0x...>
    """
    def __getitem__(self, key):
        r"""
        Return an item.

        TESTS::

            sage: from sage.rings.asymptotic.misc import Locals
            sage: locals = Locals()
            sage: locals
            {}
            sage: locals['log']  # indirect doctest
            <function log at 0x...>
        """
        try:
            return super(Locals, self).__getitem__(key)
        except KeyError as ke:
            try:
                return self.default_locals()[key]
            except KeyError:
                raise ke

    def __setitem__(self, key, value):
        r"""
        Set an item.

        This raises an error as the object is immutable.

        TESTS::

            sage: from sage.rings.asymptotic.misc import Locals
            sage: locals = Locals()
            sage: locals['a'] = 4  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: locals dictionary is frozen,
            therefore does not support item assignment
        """
        raise TypeError('locals dictionary is frozen, therefore does not support item assignment')

    @cached_method
    def _data_(self):
        r"""
        Return stored data as tuple sorted by their keys.

        TESTS::

            sage: from sage.rings.asymptotic.misc import Locals
            sage: locals = Locals({'a': 2, 'b': 1})
            sage: locals._data_()
            (('a', 2), ('b', 1))
        """
        return tuple(sorted(self.items(), key=lambda k_v: k_v[0]))

    def __hash__(self):
        r"""
        Return a hash value.

        TESTS::

            sage: from sage.rings.asymptotic.misc import Locals
            sage: locals = Locals({'a': 2, 'b': 1})
            sage: hash(locals)  # random
            -42
        """
        return hash(self._data_())

    @cached_method
    def default_locals(self):
        r"""
        Return the default locals used in
        the :class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticRing`.

        OUTPUT:

        A dictionary.

        EXAMPLES::

            sage: from sage.rings.asymptotic.misc import Locals
            sage: locals = Locals({'a': 2, 'b': 1})
            sage: locals
            {'a': 2, 'b': 1}
            sage: locals.default_locals()
            {'log': <function log at 0x...>}
            sage: locals['log']
            <function log at 0x...>
        """
        from sage.functions.log import log
        return {
            'log': log}


class WithLocals(SageObject):
    r"""
    A class extensions for handling local values; see also
    :class:`Locals`.

    This is used in the
    :class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticRing`.

    EXAMPLES::

        sage: A.<n> = AsymptoticRing('n^ZZ', QQ, locals={'a': 42})
        sage: A.locals()
        {'a': 42}
    """
    @staticmethod
    def _convert_locals_(locals):
        r"""
        This helper method return data converted to :class:`Locals`.

        TESTS::

            sage: A.<n> = AsymptoticRing('n^ZZ', QQ)
            sage: locals = A._convert_locals_({'a': 42}); locals
            {'a': 42}
            sage: type(locals)
            <class 'sage.rings.asymptotic.misc.Locals'>
        """
        if locals is None:
            return Locals()
        else:
            return Locals(locals)

    def locals(self, locals=None):
        r"""
        Return the actual :class:`Locals` object to be used.

        INPUT:

        - ``locals`` -- an object

          If ``locals`` is not ``None``, then a :class:`Locals` object
          is created and returned.
          If ``locals`` is ``None``, then a stored :class:`Locals` object,
          if any, is returned. Otherwise, an empty (i.e. no values except
          the default values)
          :class:`Locals` object is created and returned.

        OUTPUT:

        A :class:`Locals` object.

        TESTS::

            sage: A.<n> = AsymptoticRing('n^ZZ', QQ, locals={'a': 42})
            sage: A.locals()
            {'a': 42}
            sage: A.locals({'a': 41})
            {'a': 41}
            sage: A.locals({'b': -3})
            {'b': -3}
        """
        if locals is None:
            try:
                locals = self._locals_
            except AttributeError:
                pass
        return self._convert_locals_(locals)
