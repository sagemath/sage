"""
Miscellaneous functions

AUTHORS:

- William Stein

- William Stein (2006-04-26): added workaround for Windows where most
  users' home directory has a space in it.

- Robert Bradshaw (2007-09-20): Ellipsis range/iterator.

TESTS:

The following test, verifying that :trac:`16181` has been resolved, needs
to stay at the beginning of this file so that its context is not
poisoned by other tests::

    sage: sage.misc.misc.inject_variable('a', 0)
    sage: a
    0

Check the fix from :trac:`8323`::

    sage: 'name' in globals()
    False
    sage: 'func' in globals()
    False
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import time
import resource
import pdb
import warnings

import sage.misc.prandom as random
from .lazy_string import lazy_string
from sage.interfaces.quit import expect_objects
from sage.env import DOT_SAGE, HOSTNAME
from sage.misc.lazy_import import lazy_import

lazy_import("sage.misc.call", ["AttrCallObject", "attrcall", "call_method"],
            deprecation=29869)

lazy_import("sage.misc.verbose", ["verbose", "set_verbose", "set_verbose_files",
                                  "get_verbose_files", "unset_verbose_files", "get_verbose"],
            deprecation=17815)

lazy_import("sage.misc.repr", ["coeff_repr", "repr_lincomb"],
            deprecation=29892)

LOCAL_IDENTIFIER = '%s.%s' % (HOSTNAME, os.getpid())

#################################################################
# File and directory utilities
#################################################################


def sage_makedirs(dirname, mode=0o777):
    """
    Python version of ``mkdir -p``: try to create a directory, and also
    create all intermediate directories as necessary.  Succeed silently
    if the directory already exists (unlike ``os.makedirs()``).
    Raise other errors (like permission errors) normally.

    EXAMPLES::

        sage: from sage.misc.misc import sage_makedirs
        sage: sage_makedirs(DOT_SAGE) # no output

    The following fails because we are trying to create a directory in
    place of an ordinary file::

        sage: filename = tmp_filename()
        sage: sage_makedirs(filename)
        Traceback (most recent call last):
        ...
        FileExistsError: [Errno ...] File exists: ...
    """
    try:
        os.makedirs(dirname)
    except OSError:
        if not os.path.isdir(dirname):
            raise


# We create the DOT_SAGE directory (if it does not exist yet; note in particular
# that it may already have been created by the bin/sage script) with
# restrictive permissions, since otherwise possibly just anybody can easily see
# every command you type.

sage_makedirs(DOT_SAGE, mode=0o700)


def try_read(obj, splitlines=False):
    r"""
    Determine if a given object is a readable file-like object and if so
    read and return its contents.

    That is, the object has a callable method named ``read()`` which takes
    no arguments (except ``self``) then the method is executed and the
    contents are returned.

    Alternatively, if the ``splitlines=True`` is given, first ``splitlines()``
    is tried, then if that fails ``read().splitlines()``.

    If either method fails, ``None`` is returned.

    INPUT:

    - ``obj`` -- typically a `file` or `io.BaseIO` object, but any other
      object with a ``read()`` method is accepted.

    - ``splitlines`` -- `bool`, optional; if True, return a list of lines
      instead of a string.

    EXAMPLES::

        sage: import io
        sage: filename = tmp_filename()
        sage: from sage.misc.misc import try_read
        sage: with open(filename, 'w') as fobj:
        ....:     _ = fobj.write('a\nb\nc')
        sage: with open(filename) as fobj:
        ....:     print(try_read(fobj))
        a
        b
        c
        sage: with open(filename) as fobj:
        ....:     try_read(fobj, splitlines=True)
        ['a\n', 'b\n', 'c']

    The following example is identical to the above example on Python 3,
    but different on Python 2 where ``open != io.open``::

        sage: with io.open(filename) as fobj:
        ....:     print(try_read(fobj))
        a
        b
        c

    I/O buffers::

        sage: buf = io.StringIO('a\nb\nc')
        sage: print(try_read(buf))
        a
        b
        c
        sage: _ = buf.seek(0); try_read(buf, splitlines=True)
        ['a\n', 'b\n', 'c']
        sage: buf = io.BytesIO(b'a\nb\nc')
        sage: try_read(buf) == b'a\nb\nc'
        True
        sage: _ = buf.seek(0)
        sage: try_read(buf, splitlines=True) == [b'a\n', b'b\n', b'c']
        True

    Custom readable::

        sage: class MyFile(object):
        ....:     def read(self): return 'Hello world!'
        sage: try_read(MyFile())
        'Hello world!'
        sage: try_read(MyFile(), splitlines=True)
        ['Hello world!']

    Not readable::

        sage: try_read(1) is None
        True
    """

    if splitlines:
        try:
            return obj.readlines()
        except (AttributeError, TypeError):
            pass

    try:
        data = obj.read()
    except (AttributeError, TypeError):
        return

    if splitlines:
        try:
            data = data.splitlines()
        except (AttributeError, TypeError):
            # Not a string??
            data = [data]

    return data


#################################################
# Next we create the Sage temporary directory.
#################################################


@lazy_string
def SAGE_TMP():
    """
    EXAMPLES::

        sage: from sage.misc.misc import SAGE_TMP
        sage: SAGE_TMP
        l'.../temp/...'
    """
    d = os.path.join(DOT_SAGE, 'temp', HOSTNAME, str(os.getpid()))
    sage_makedirs(d)
    return d


@lazy_string
def ECL_TMP():
    """
    Temporary directory that should be used by ECL interfaces launched from
    Sage.

    EXAMPLES::

        sage: from sage.misc.misc import ECL_TMP
        sage: ECL_TMP
        l'.../temp/.../ecl'
    """
    d = os.path.join(str(SAGE_TMP), 'ecl')
    sage_makedirs(d)
    return d


@lazy_string
def SPYX_TMP():
    """
    EXAMPLES::

        sage: from sage.misc.misc import SPYX_TMP
        sage: SPYX_TMP
        l'.../temp/.../spyx'
    """
    return os.path.join(str(SAGE_TMP), 'spyx')


@lazy_string
def SAGE_TMP_INTERFACE():
    """
    EXAMPLES::

        sage: from sage.misc.misc import SAGE_TMP_INTERFACE
        sage: SAGE_TMP_INTERFACE
        l'.../temp/.../interface'
    """
    d = os.path.join(str(SAGE_TMP), 'interface')
    sage_makedirs(d)
    return d


SAGE_DB = os.path.join(DOT_SAGE, 'db')
sage_makedirs(SAGE_DB)

try:
    # Create the matplotlib config directory.
    sage_makedirs(os.environ["MPLCONFIGDIR"])
except KeyError:
    pass

#################################################################
# timing
#################################################################


def cputime(t=0, subprocesses=False):
    """
    Return the time in CPU seconds since Sage started, or with
    optional argument ``t``, return the time since ``t``. This is how
    much time Sage has spent using the CPU.  If ``subprocesses=False``
    this does not count time spent in subprocesses spawned by Sage
    (e.g., Gap, Singular, etc.). If ``subprocesses=True`` this
    function tries to take all subprocesses with a working
    ``cputime()`` implementation into account.

    The measurement for the main Sage process is done via a call to
    :func:`resource.getrusage()`, so it avoids the wraparound problems in
    :func:`time.clock()` on Cygwin.

    INPUT:

    - ``t`` - (optional) time in CPU seconds, if ``t`` is a result
      from an earlier call with ``subprocesses=True``, then
      ``subprocesses=True`` is assumed.

    - subprocesses -- (optional), include subprocesses (default:
      ``False``)

    OUTPUT:

    - ``float`` - time in CPU seconds if ``subprocesses=False``

    - :class:`GlobalCputime` - object which holds CPU times of
      subprocesses otherwise

    EXAMPLES::

        sage: t = cputime()
        sage: F = gp.factor(2^199-1)
        sage: cputime(t)          # somewhat random
        0.010999000000000092

        sage: t = cputime(subprocesses=True)
        sage: F = gp.factor(2^199-1)
        sage: cputime(t) # somewhat random
        0.091999

        sage: w = walltime()
        sage: F = gp.factor(2^199-1)
        sage: walltime(w)         # somewhat random
        0.58425593376159668

    .. NOTE::

        Even with ``subprocesses=True`` there is no guarantee that the
        CPU time is reported correctly because subprocesses can be
        started and terminated at any given time.
    """
    if isinstance(t, GlobalCputime):
        subprocesses = True

    if not subprocesses:
        try:
            t = float(t)
        except TypeError:
            t = 0.0
        u, s = resource.getrusage(resource.RUSAGE_SELF)[:2]
        return u + s - t
    else:
        if t == 0:
            ret = GlobalCputime(cputime())
            for s in expect_objects:
                S = s()
                if S and S.is_running():
                    try:
                        ct = S.cputime()
                        ret.total += ct
                        ret.interfaces[s] = ct
                    except NotImplementedError:
                        pass
            return ret
        else:
            if not isinstance(t, GlobalCputime):
                t = GlobalCputime(t)
            ret = GlobalCputime(cputime() - t.local)
            for s in expect_objects:
                S = s()
                if S and S.is_running():
                    try:
                        ct = S.cputime() - t.interfaces.get(s, 0.0)
                        ret.total += ct
                        ret.interfaces[s] = ct
                    except NotImplementedError:
                        pass
            return ret


class GlobalCputime:
    """
    Container for CPU times of subprocesses.

    AUTHOR:

    - Martin Albrecht - (2008-12): initial version

    EXAMPLES:

    Objects of this type are returned if ``subprocesses=True`` is
    passed to :func:`cputime`::

        sage: cputime(subprocesses=True) # indirect doctest, output random
        0.2347431

    We can use it to keep track of the CPU time spent in Singular for
    example::

        sage: t = cputime(subprocesses=True)
        sage: P = PolynomialRing(QQ,7,'x')
        sage: I = sage.rings.ideal.Katsura(P)
        sage: gb = I.groebner_basis() # calls Singular
        sage: cputime(subprocesses=True) - t # output random
        0.462987

    For further processing we can then convert this container to a
    float::

        sage: t = cputime(subprocesses=True)
        sage: float(t) #output somewhat random
        2.1088339999999999

    .. SEEALSO::

      :func:`cputime`
    """
    def __init__(self, t):
        """
        Create a new CPU time object which also keeps track of
        subprocesses.

        EXAMPLES::

            sage: from sage.misc.misc import GlobalCputime
            sage: ct = GlobalCputime(0.0); ct
            0.0...
        """
        self.total = t
        self.local = t
        self.interfaces = {}

    def __repr__(self):
        """
        EXAMPLES::

            sage: cputime(subprocesses=True) # indirect doctest, output random
            0.2347431
        """
        return str(self.total)

    def __add__(self, other):
        """
        EXAMPLES::

            sage: t = cputime(subprocesses=True)
            sage: P = PolynomialRing(QQ,7,'x')
            sage: I = sage.rings.ideal.Katsura(P)
            sage: gb = I.groebner_basis() # calls Singular
            sage: cputime(subprocesses=True) + t # output random
            2.798708
        """
        if not isinstance(other, GlobalCputime):
            other = GlobalCputime(other)
        ret = GlobalCputime(self.total + other.total)
        return ret

    def __sub__(self, other):
        """
        EXAMPLES::

            sage: t = cputime(subprocesses=True)
            sage: P = PolynomialRing(QQ,7,'x')
            sage: I = sage.rings.ideal.Katsura(P)
            sage: gb = I.groebner_basis() # calls Singular
            sage: cputime(subprocesses=True) - t # output random
            0.462987
        """
        if not isinstance(other, GlobalCputime):
            other = GlobalCputime(other)
        ret = GlobalCputime(self.total - other.total)
        return ret

    def __float__(self):
        """
        EXAMPLES::

            sage: t = cputime(subprocesses=True)
            sage: float(t) #output somewhat random
            2.1088339999999999
        """
        return float(self.total)


def walltime(t=0):
    """
    Return the wall time in second, or with optional argument t, return
    the wall time since time t. "Wall time" means the time on a wall
    clock, i.e., the actual time.

    INPUT:


    -  ``t`` - (optional) float, time in CPU seconds

    OUTPUT:

    -  ``float`` - time in seconds


    EXAMPLES::

        sage: w = walltime()
        sage: F = factor(2^199-1)
        sage: walltime(w)   # somewhat random
        0.8823847770690918
    """
    return time.time() - t


def union(x, y=None):
    """
    Return the union of x and y, as a list. The resulting list need not
    be sorted and can change from call to call.

    INPUT:


    -  ``x`` - iterable

    -  ``y`` - iterable (may optionally omitted)


    OUTPUT: list

    EXAMPLES::

        sage: answer = union([1,2,3,4], [5,6]); answer
        doctest:...: DeprecationWarning: sage.misc.misc.union is deprecated...
        See https://trac.sagemath.org/32096 for details.
        [1, 2, 3, 4, 5, 6]
        sage: union([1,2,3,4,5,6], [5,6]) == answer
        True
        sage: union((1,2,3,4,5,6), [5,6]) == answer
        True
        sage: union((1,2,3,4,5,6), set([5,6])) == answer
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(32096, "sage.misc.misc.union is deprecated, use 'list(set(x).union(y)' or a more suitable replacement")
    if y is None:
        return list(set(x))
    return list(set(x).union(y))


def uniq(x):
    """
    Return the sublist of all elements in the list x that is sorted and
    is such that the entries in the sublist are unique.

    EXAMPLES::

        sage: uniq([1, 1, 8, -5, 3, -5, -13, 13, -13])
        doctest:...: DeprecationWarning: the output of uniq(X) being sorted is deprecated; use sorted(set(X)) instead if you want sorted output
        See https://trac.sagemath.org/27014 for details.
        [-13, -5, 1, 3, 8, 13]
    """
    # After deprecation period, rename _stable_uniq -> uniq
    from sage.misc.superseded import deprecation
    deprecation(27014, "the output of uniq(X) being sorted is deprecated; use sorted(set(X)) instead if you want sorted output")
    return sorted(set(x))


def _stable_uniq(L):
    """
    Iterate over the elements of ``L``, yielding every element at most
    once: keep only the first occurance of any item.

    The items must be hashable.

    INPUT:

    - ``L`` -- iterable

    EXAMPLES::

        sage: from sage.misc.misc import _stable_uniq
        sage: L = [1, 1, 8, -5, 3, -5, 'a', 'x', 'a']
        sage: it = _stable_uniq(L)
        sage: it
        <generator object _stable_uniq at ...>
        sage: list(it)
        [1, 8, -5, 3, 'a', 'x']
    """
    seen = set()
    for x in L:
        if x in seen:
            continue
        yield x
        seen.add(x)


def exactly_one_is_true(iterable):
    r"""
    Return whether exactly one element of ``iterable`` evaluates ``True``.

    INPUT:

    - ``iterable`` -- an iterable object

    OUTPUT:

    A boolean.

    .. NOTE::

        The implementation is suggested by
        `stackoverflow entry <https://stackoverflow.com/a/16801605/1052778>`_.

    EXAMPLES::

        sage: from sage.misc.misc import exactly_one_is_true
        sage: exactly_one_is_true([])
        False
        sage: exactly_one_is_true([True])
        True
        sage: exactly_one_is_true([False])
        False
        sage: exactly_one_is_true([True, True])
        False
        sage: exactly_one_is_true([False, True])
        True
        sage: exactly_one_is_true([True, False, True])
        False
        sage: exactly_one_is_true([False, True, False])
        True
    """
    it = iter(iterable)
    return any(it) and not any(it)


def strunc(s, n=60):
    """
    Truncate at first space after position n, adding '...' if
    nontrivial truncation.
    """
    n = int(n)
    s = str(s)
    if len(s) > n:
        i = n
        while i < len(s) and s[i] != ' ':
            i += 1
        return s[:i] + " ..."
    return s


def newton_method_sizes(N):
    r"""
    Return a sequence of integers
    `1 = a_1 \leq a_2 \leq \cdots \leq a_n = N` such that
    `a_j = \lceil a_{j+1} / 2 \rceil` for all `j`.

    This is useful for Newton-style algorithms that double the
    precision at each stage. For example if you start at precision 1
    and want an answer to precision 17, then it's better to use the
    intermediate stages 1, 2, 3, 5, 9, 17 than to use 1, 2, 4, 8, 16,
    17.

    INPUT:


    -  ``N`` - positive integer


    EXAMPLES::

        sage: newton_method_sizes(17)
        [1, 2, 3, 5, 9, 17]
        sage: newton_method_sizes(16)
        [1, 2, 4, 8, 16]
        sage: newton_method_sizes(1)
        [1]

    AUTHORS:

    - David Harvey (2006-09-09)
    """

    N = int(N)
    if N < 1:
        raise ValueError("N (={}) must be a positive integer".format(N))

    output = []
    while N > 1:
        output.append(N)
        N = (N + 1) >> 1

    output.append(1)
    output.reverse()
    return output


#################################################################
# Generally useful
#################################################################


def compose(f, g):
    r"""
    Return the composition of one-variable functions: `f \circ g`

    See also :func:`nest()`

    INPUT:
        - `f` -- a function of one variable
        - `g` -- another function of one variable

    OUTPUT:
        A function, such that compose(f,g)(x) = f(g(x))

    EXAMPLES::

        sage: def g(x): return 3*x
        sage: def f(x): return x + 1
        sage: h1 = compose(f,g)
        sage: h2 = compose(g,f)
        sage: _ = var ('x')
        sage: h1(x)
        3*x + 1
        sage: h2(x)
        3*x + 3

    ::

        sage: _ = function('f g')
        sage: _ = var ('x')
        sage: compose(f,g)(x)
        f(g(x))

    """
    return lambda x: f(g(x))


def nest(f, n, x):
    """
    Return `f(f(...f(x)...))`, where the composition occurs n times.

    See also :func:`compose()` and :func:`self_compose()`

    INPUT:
        - `f` -- a function of one variable
        - `n` -- a nonnegative integer
        - `x` -- any input for `f`

    OUTPUT:
        `f(f(...f(x)...))`, where the composition occurs n times

    EXAMPLES::

        sage: def f(x): return x^2 + 1
        sage: x = var('x')
        sage: nest(f, 3, x)
        ((x^2 + 1)^2 + 1)^2 + 1

    ::

        sage: _ = function('f')
        sage: _ = var('x')
        sage: nest(f, 10, x)
        f(f(f(f(f(f(f(f(f(f(x))))))))))

    ::

        sage: _ = function('f')
        sage: _ = var('x')
        sage: nest(f, 0, x)
        x

    """
    from sage.rings.integer import Integer
    n = Integer(n)

    if n < 0:
        raise ValueError("n must be a nonnegative integer, not {}.".format(n))

    for i in range(n):
        x = f(x)
    return x


#################################################################
# The A \ b operator
#################################################################

class BackslashOperator:
    r"""
    Implements Matlab-style backslash operator for solving systems::

        A \ b

    The preparser converts this to multiplications using
    ``BackslashOperator()``.

    EXAMPLES::

        sage: preparse("A \\ matrix(QQ,2,1,[1/3,'2/3'])")
        "A  * BackslashOperator() * matrix(QQ,Integer(2),Integer(1),[Integer(1)/Integer(3),'2/3'])"
        sage: preparse("A \\ matrix(QQ,2,1,[1/3,2*3])")
        'A  * BackslashOperator() * matrix(QQ,Integer(2),Integer(1),[Integer(1)/Integer(3),Integer(2)*Integer(3)])'
        sage: preparse("A \\ B + C")
        'A  * BackslashOperator() * B + C'
        sage: preparse("A \\ eval('C+D')")
        "A  * BackslashOperator() * eval('C+D')"
        sage: preparse("A \\ x / 5")
        'A  * BackslashOperator() * x / Integer(5)'
        sage: preparse("A^3 \\ b")
        'A**Integer(3)  * BackslashOperator() * b'
    """
    def __rmul__(self, left):
        """
        EXAMPLES::

            sage: A = random_matrix(ZZ, 4)
            sage: while A.rank() != 4:
            ....:     A = random_matrix(ZZ, 4)
            sage: B = random_matrix(ZZ, 4)
            sage: temp = A * BackslashOperator()
            sage: temp.left is A
            True
            sage: X = temp * B
            sage: A * X == B
            True
        """
        self.left = left
        return self

    def __mul__(self, right):
        r"""
        EXAMPLES::

            sage: A = matrix(RDF, 5, 5, 2)
            sage: b = vector(RDF, 5, range(5))
            sage: v = A \ b
            sage: v.zero_at(1e-19)  # On at least one platform, we get a "negative zero"
            (0.0, 0.5, 1.0, 1.5, 2.0)
            sage: v = A._backslash_(b)
            sage: v.zero_at(1e-19)
            (0.0, 0.5, 1.0, 1.5, 2.0)
            sage: v = A * BackslashOperator() * b
            sage: v.zero_at(1e-19)
            (0.0, 0.5, 1.0, 1.5, 2.0)
        """
        return self.left._backslash_(right)


#################################################################
# is_iterator function
#################################################################
def is_iterator(it) -> bool:
    """
    Tests if it is an iterator.

    The mantra ``if hasattr(it, 'next')`` was used to tests if ``it`` is an
    iterator. This is not quite correct since ``it`` could have a ``next``
    methods with a different semantic.

    EXAMPLES::

        sage: it = iter([1,2,3])
        sage: is_iterator(it)
        True

        sage: class wrong():
        ....:    def __init__(self): self.n = 5
        ....:    def __next__(self):
        ....:        self.n -= 1
        ....:        if self.n == 0: raise StopIteration
        ....:        return self.n
        sage: x = wrong()
        sage: is_iterator(x)
        False
        sage: list(x)
        Traceback (most recent call last):
        ...
        TypeError: 'wrong' object is not iterable

        sage: class good(wrong):
        ....:    def __iter__(self): return self
        sage: x = good()
        sage: is_iterator(x)
        True
        sage: list(x)
        [4, 3, 2, 1]

        sage: P = Partitions(3)
        sage: is_iterator(P)
        False
        sage: is_iterator(iter(P))
        True
    """
    # see trac #7398 for a discussion
    try:
        return it is iter(it)
    except Exception:
        return False


#################################################################
# Useful but hard to classify
#################################################################


def random_sublist(X, s):
    """
    Return a pseudo-random sublist of the list X where the probability
    of including a particular element is s.

    INPUT:


    -  ``X`` - list

    -  ``s`` - floating point number between 0 and 1


    OUTPUT: list

    EXAMPLES::

        sage: from sage.misc.misc import is_sublist
        sage: S = [1,7,3,4,18]
        sage: sublist = random_sublist(S, 0.5); sublist  # random
        [1, 3, 4]
        sage: is_sublist(sublist, S)
        True
        sage: sublist = random_sublist(S, 0.5); sublist  # random
        [1, 3]
        sage: is_sublist(sublist, S)
        True
    """
    return [a for a in X if random.random() <= s]


def is_sublist(X, Y):
    """
    Test whether ``X`` is a sublist of ``Y``.

    EXAMPLES::

        sage: from sage.misc.misc import is_sublist
        sage: S = [1, 7, 3, 4, 18]
        sage: is_sublist([1, 7], S)
        True
        sage: is_sublist([1, 3, 4], S)
        True
        sage: is_sublist([1, 4, 3], S)
        False
        sage: is_sublist(S, S)
        True
    """
    X_i = 0
    for Y_i, y in enumerate(Y):
        if X_i == len(X):
            return True
        if y == X[X_i]:
            X_i += 1
    return X_i == len(X)


def some_tuples(elements, repeat, bound, max_samples=None):
    r"""
    Return an iterator over at most ``bound`` number of ``repeat``-tuples of
    ``elements``.

    INPUT:

    - ``elements`` -- an iterable
    - ``repeat`` -- integer (default ``None``), the length of the tuples to be returned.
      If ``None``, just returns entries from ``elements``.
    - ``bound`` -- the maximum number of tuples returned (ignored if ``max_samples`` given)
    - ``max_samples`` -- non-negative integer (default ``None``).  If given,
      then a sample of the possible tuples will be returned,
      instead of the first few in the standard order.

    OUTPUT:

    If ``max_samples`` is not provided, an iterator over the first
    ``bound`` tuples of length ``repeat``, in the standard nested-for-loop order.

    If ``max_samples`` is provided, a list of at most ``max_samples`` tuples,
    sampled uniformly from the possibilities.  In this case, ``elements``
    must be finite.

    TESTS::

        sage: from sage.misc.misc import some_tuples
        sage: l = some_tuples([0,1,2,3], 2, 3)
        sage: l
        <itertools.islice object at ...>
        sage: len(list(l))
        3

        sage: l = some_tuples(range(50), 3, 10)
        sage: len(list(l))
        10

        sage: l = some_tuples(range(3), 2, None, max_samples=10)
        sage: len(list(l))
        9
    """
    if max_samples is None:
        from itertools import islice, product
        P = elements if repeat is None else product(elements, repeat=repeat)
        return islice(P, int(bound))
    else:
        if not (hasattr(elements, '__len__') and hasattr(elements, '__getitem__')):
            elements = list(elements)
        n = len(elements)
        N = n if repeat is None else n**repeat
        if N <= max_samples:
            from itertools import product
            return elements if repeat is None else product(elements, repeat=repeat)
        return _some_tuples_sampling(elements, repeat, max_samples, n)


def _some_tuples_sampling(elements, repeat, max_samples, n):
    """
    Internal function for :func:`some_tuples`.

    TESTS::

        sage: from sage.misc.misc import _some_tuples_sampling
        sage: l = list(_some_tuples_sampling(range(3), 3, 2, 3))
        sage: len(l)
        2
        sage: all(len(tup) == 3 for tup in l)
        True
        sage: all(el in range(3) for tup in l for el in tup)
        True
        sage: l = list(_some_tuples_sampling(range(20), None, 4, 20))
        sage: len(l)
        4
        sage: all(el in range(20) for el in l)
        True
    """
    from sage.rings.integer import Integer
    N = n if repeat is None else n**repeat
    # We sample on range(N) and create tuples manually since we don't want to create the list of all possible tuples in memory
    for a in random.sample(range(N), max_samples):
        if repeat is None:
            yield elements[a]
        else:
            yield tuple(elements[j] for j in Integer(a).digits(n, padto=repeat))


def powerset(X):
    r"""
    Iterator over the *list* of all subsets of the iterable X, in no
    particular order. Each list appears exactly once, up to order.

    INPUT:

    -  ``X`` - an iterable

    OUTPUT: iterator of lists

    EXAMPLES::

        sage: list(powerset([1,2,3]))
        [[], [1], [2], [1, 2], [3], [1, 3], [2, 3], [1, 2, 3]]
        sage: [z for z in powerset([0,[1,2]])]
        [[], [0], [[1, 2]], [0, [1, 2]]]

    Iterating over the power set of an infinite set is also allowed::

        sage: i = 0
        sage: L = []
        sage: for x in powerset(ZZ):
        ....:     if i > 10:
        ....:         break
        ....:     else:
        ....:         i += 1
        ....:     L.append(x)
        sage: print(" ".join(str(x) for x in L))
        [] [0] [1] [0, 1] [-1] [0, -1] [1, -1] [0, 1, -1] [2] [0, 2] [1, 2]

    You may also use subsets as an alias for powerset::

        sage: subsets([1,2,3])
        <generator object ...powerset at 0x...>
        sage: list(subsets([1,2,3]))
        [[], [1], [2], [1, 2], [3], [1, 3], [2, 3], [1, 2, 3]]

        The reason we return lists instead of sets is that the elements of
        sets must be hashable and many structures on which one wants the
        powerset consist of non-hashable objects.

    AUTHORS:

    - William Stein

    - Nils Bruin (2006-12-19): rewrite to work for not-necessarily
      finite objects X.
    """
    yield []
    pairs = []
    power2 = 1
    for x in X:
        pairs.append((power2, x))
        next_power2 = power2 << 1
        for w in range(power2, next_power2):
            yield [x for m, x in pairs if m & w]
        power2 = next_power2


subsets = powerset


#################################################################
# Misc.
#################################################################

def exists(S, P):
    """
    If S contains an element x such that P(x) is True, this function
    returns True and the element x. Otherwise it returns False and
    None.

    Note that this function is NOT suitable to be used in an
    if-statement or in any place where a boolean expression is
    expected. For those situations, use the Python built-in

    any(P(x) for x in S)

    INPUT:


    -  ``S`` - object (that supports enumeration)

    -  ``P`` - function that returns True or False


    OUTPUT:


    -  ``bool`` - whether or not P is True for some element
       x of S

    -  ``object`` - x


    EXAMPLES: lambda functions are very useful when using the exists
    function::

        sage: exists([1,2,5], lambda x : x > 7)
        (False, None)
        sage: exists([1,2,5], lambda x : x > 3)
        (True, 5)

    The following example is similar to one in the MAGMA handbook. We
    check whether certain integers are a sum of two (small) cubes::

        sage: cubes = [t**3 for t in range(-10,11)]
        sage: exists([(x,y) for x in cubes for y in cubes], lambda v : v[0]+v[1] == 218)
        (True, (-125, 343))
        sage: exists([(x,y) for x in cubes for y in cubes], lambda v : v[0]+v[1] == 219)
        (False, None)
    """
    for x in S:
        if P(x):
            return True, x
    return False, None


def forall(S, P):
    """
    If P(x) is true every x in S, return True and None. If there is
    some element x in S such that P is not True, return False and x.

    Note that this function is NOT suitable to be used in an
    if-statement or in any place where a boolean expression is
    expected. For those situations, use the Python built-in

    all(P(x) for x in S)

    INPUT:

    -  ``S`` - object (that supports enumeration)

    -  ``P`` - function that returns True or False

    OUTPUT:


    -  ``bool`` - whether or not P is True for all elements
       of S

    -  ``object`` - x


    EXAMPLES: lambda functions are very useful when using the forall
    function. As a toy example we test whether certain integers are
    greater than 3.

    ::

        sage: forall([1,2,5], lambda x : x > 3)
        (False, 1)
        sage: forall([1,2,5], lambda x : x > 0)
        (True, None)

    Next we ask whether every positive integer less than 100 is a
    product of at most 2 prime factors::

        sage: forall(range(1,100),  lambda n : len(factor(n)) <= 2)
        (False, 30)

    The answer is no, and 30 is a counterexample. However, every
    positive integer 100 is a product of at most 3 primes.

    ::

        sage: forall(range(1,100),  lambda n : len(factor(n)) <= 3)
        (True, None)
    """
    for x in S:
        if not P(x):
            return False, x
    return True, None


#################################################################
# debug tracing
#################################################################
set_trace = pdb.set_trace


#################################################################
# Word wrap lines
#################################################################
def word_wrap(s, ncols=85):
    t = []
    if ncols == 0:
        return s
    for x in s.split('\n'):
        if not x or x.lstrip()[:5] == 'sage:':
            t.append(x)
            continue
        while len(x) > ncols:
            k = ncols
            while k > 0 and x[k] != ' ':
                k -= 1
            if k == 0:
                k = ncols
                end = '\\'
            else:
                end = ''
            t.append(x[:k] + end)
            x = x[k:]
            k = 0
            while k < len(x) and x[k] == ' ':
                k += 1
            x = x[k:]
        t.append(x)
    return '\n'.join(t)


def pad_zeros(s, size=3):
    """
    EXAMPLES::

        sage: pad_zeros(100)
        '100'
        sage: pad_zeros(10)
        '010'
        sage: pad_zeros(10, 5)
        '00010'
        sage: pad_zeros(389, 5)
        '00389'
        sage: pad_zeros(389, 10)
        '0000000389'
    """
    return "0" * (size - len(str(s))) + str(s)


def is_in_string(line, pos):
    r"""
    Return ``True`` if the character at position ``pos`` in ``line`` occurs
    within a string.

    EXAMPLES::

        sage: from sage.misc.misc import is_in_string
        sage: line = 'test(\'#\')'
        sage: is_in_string(line, line.rfind('#'))
        True
        sage: is_in_string(line, line.rfind(')'))
        False
    """
    i = 0
    in_single_quote = False
    in_double_quote = False
    in_triple_quote = False

    def in_quote():
        return in_single_quote or in_double_quote or in_triple_quote

    while i < pos:
        # Update quote parsing
        # We only do this if this quote isn't backquoted itself,
        # which is the case if the previous character isn't
        # a backslash, or it is but both previous characters
        # are backslashes.
        if line[i - 1: i] != '\\' or line[i - 2: i] == '\\\\':
            if line[i: i + 3] in ['"""', "'''"]:
                if not in_quote():
                    in_triple_quote = True
                elif in_triple_quote:
                    in_triple_quote = False
            elif line[i] == "'":
                if not in_quote():
                    in_single_quote = True
                elif in_single_quote:
                    in_single_quote = False
            elif line[i] == '"':
                if not in_quote():
                    in_double_quote = True
                elif in_double_quote:
                    in_double_quote = False
        i += 1
    return in_quote()


def get_main_globals():
    """
    Return the main global namespace.

    EXAMPLES::

        sage: from sage.misc.misc import get_main_globals
        sage: G = get_main_globals()
        sage: bla = 1
        sage: G['bla']
        1
        sage: bla = 2
        sage: G['bla']
        2
        sage: G['ble'] = 5
        sage: ble
        5

    This is analogous to :func:`globals`, except that it can be called
    from any function, even if it is in a Python module::

        sage: def f():
        ....:     G = get_main_globals()
        ....:     assert G['bli'] == 14
        ....:     G['blo'] = 42
        sage: bli = 14
        sage: f()
        sage: blo
        42

    ALGORITHM:

    The main global namespace is discovered by going up the frame
    stack until the frame for the :mod:`__main__` module is found.
    Should this frame not be found (this should not occur in normal
    operation), an exception "ValueError: call stack is not deep
    enough" will be raised by ``_getframe``.

    See :meth:`inject_variable_test` for a real test that this works
    within deeply nested calls in a function defined in a Python
    module.
    """
    import sys
    depth = 0
    while True:
        G = sys._getframe(depth).f_globals
        if G.get("__name__", None) == "__main__":
            break
        depth += 1
    return G


def inject_variable(name, value, warn=True):
    """
    Inject a variable into the main global namespace.

    INPUT:

    - ``name``  -- a string
    - ``value`` -- anything
    - ``warn`` -- a boolean (default: :obj:`False`)

    EXAMPLES::

        sage: from sage.misc.misc import inject_variable
        sage: inject_variable("a", 314)
        sage: a
        314

    A warning is issued the first time an existing value is overwritten::

        sage: inject_variable("a", 271)
        doctest:...: RuntimeWarning: redefining global value `a`
        sage: a
        271
        sage: inject_variable("a", 272)
        sage: a
        272

    That's because warn seem to not reissue twice the same warning:

        sage: from warnings import warn
        sage: warn("blah")
        doctest:...: UserWarning: blah
        sage: warn("blah")

    Warnings can be disabled::

        sage: b = 3
        sage: inject_variable("b", 42, warn=False)
        sage: b
        42

    Use with care!
    """
    assert isinstance(name, str)
    # Using globals() does not work, even in Cython, because
    # inject_variable is called not only from the interpreter, but
    # also from functions in various modules.
    G = get_main_globals()
    if name in G and warn:
        warnings.warn("redefining global value `%s`" % name,
                      RuntimeWarning, stacklevel=2)
    G[name] = value


def inject_variable_test(name, value, depth):
    """
    A function for testing deep calls to inject_variable

    EXAMPLES::

        sage: from sage.misc.misc import inject_variable_test
        sage: inject_variable_test("a0", 314, 0)
        sage: a0
        314
        sage: inject_variable_test("a1", 314, 1)
        sage: a1
        314
        sage: inject_variable_test("a2", 314, 2)
        sage: a2
        314
        sage: inject_variable_test("a2", 271, 2)
        doctest:...: RuntimeWarning: redefining global value `a2`
        sage: a2
        271
    """
    if depth == 0:
        inject_variable(name, value)
    else:
        inject_variable_test(name, value, depth - 1)
