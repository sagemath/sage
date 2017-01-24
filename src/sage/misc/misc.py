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

Test deprecation::

    sage: sage.misc.misc.srange(5)
    doctest:...: DeprecationWarning:
    Importing srange from here is deprecated. If you need to use it, please import it directly from sage.arith.srange
    See http://trac.sagemath.org/20094 for details.
    [0, 1, 2, 3, 4]
    sage: sage.misc.all.srange(5)
    doctest:...: DeprecationWarning:
    Importing srange from here is deprecated. If you need to use it, please import it directly from sage.arith.srange
    See http://trac.sagemath.org/20334 for details.
    [0, 1, 2, 3, 4]
    sage: sage.misc.misc.sxrange(5)
    doctest:...: DeprecationWarning:
    Importing sxrange from here is deprecated. If you need to use it, please import it directly from sage.arith.srange
    See http://trac.sagemath.org/20094 for details.
    <generator object at 0x...>
    sage: sage.misc.misc.cancel_alarm()
    doctest:...: DeprecationWarning:
    Importing cancel_alarm from here is deprecated. If you need to use it, please import it directly from cysignals.alarm
    See http://trac.sagemath.org/20002 for details.
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import
from six.moves import range

__doc_exclude=["cached_attribute", "cached_class_attribute", "lazy_prop",
               "generic_cmp", "to_gmp_hex", "todo",
               "typecheck", "prop", "strunc",
               "assert_attribute", "LOGFILE"]

from warnings import warn
import os
import stat
import sys
import time
import resource
import sage.misc.prandom as random
from .lazy_string import lazy_string

from sage.misc.lazy_import import lazy_import
lazy_import('sage.arith.srange', ('xsrange', 'srange', 'ellipsis_range', 'ellipsis_iter'), deprecation=20094)
lazy_import('sage.arith.srange', 'xsrange', 'sxrange', deprecation=20094)
lazy_import('cysignals.alarm', ('alarm', 'cancel_alarm'), deprecation=20002)


from sage.env import DOT_SAGE, HOSTNAME

LOCAL_IDENTIFIER = '%s.%s'%(HOSTNAME , os.getpid())

def sage_makedirs(dir):
    """
    Python version of ``mkdir -p``: try to create a directory, and also
    create all intermediate directories as necessary.  Succeed silently
    if the directory already exists (unlike ``os.makedirs()``).
    Raise other errors (like permission errors) normally.

    EXAMPLES::

        sage: from sage.misc.misc import sage_makedirs
        sage: sage_makedirs(DOT_SAGE) # no output

    The following fails because we are trying to create a directory in
    place of an ordinary file (the main Sage executable)::

        sage: sage_executable = os.path.join(SAGE_ROOT, 'sage')
        sage: sage_makedirs(sage_executable)
        Traceback (most recent call last):
        ...
        OSError: ...
    """
    try:
        os.makedirs(dir)
    except OSError:
        if not os.path.isdir(dir):
            raise


#################################################
# Now that the variable DOT_SAGE has been set,
# we make sure that the DOT_SAGE directory
# has restrictive permissions, since otherwise
# possibly just anybody can easily see every
# command you type, since it is in the history,
# and every worksheet you create, etc.
# We do the following:
#   1. If there is no DOT_SAGE, we create it.
#   2. Check to see if the permissions on DOT_SAGE are
#      sufficiently restrictive.  If not, we change them.

sage_makedirs(DOT_SAGE)

if hasattr(os, 'chmod'):
    _mode = os.stat(DOT_SAGE)[stat.ST_MODE]
    _desired_mode = 0o40700     # drwx------
    if _mode != _desired_mode:
        # On Cygwin, if the sage directory is not in a filesystem mounted with
        # 'acl' support, setting the permissions may fail silently, so only
        # print the message after we've changed the permissions and confirmed
        # that the change succeeded
        os.chmod(DOT_SAGE, _desired_mode)
        if os.stat(DOT_SAGE)[stat.ST_MODE] == _desired_mode:
            print("Setting permissions of DOT_SAGE directory so only you "
                  "can read and write it.")


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
def SPYX_TMP():
    """
    EXAMPLES::

        sage: from sage.misc.misc import SPYX_TMP
        sage: SPYX_TMP
        l'.../temp/.../spyx'
    """
    return os.path.join(SAGE_TMP, 'spyx')

@lazy_string
def SAGE_TMP_INTERFACE():
    """
    EXAMPLES::

        sage: from sage.misc.misc import SAGE_TMP_INTERFACE
        sage: SAGE_TMP_INTERFACE
        l'.../temp/.../interface'
    """
    d = os.path.join(SAGE_TMP, 'interface')
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
# Functions to help with interfacing with CXX code that
# uses the GMP library
#################################################################
def to_gmp_hex(n):
    from sage.misc.superseded import deprecation
    deprecation(21926, "to_gmp_hex() is deprecated")
    return hex(n).replace("L","").replace("0x","")

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

    .. note ::

      Even with ``subprocesses=True`` there is no guarantee that the
      CPU time is reported correctly because subprocesses can be
      started and terminated at any given time.
    """
    if isinstance(t, GlobalCputime):
        subprocesses=True

    if not subprocesses:
        try:
            t = float(t)
        except TypeError:
            t = 0.0
        u,s = resource.getrusage(resource.RUSAGE_SELF)[:2]
        return u+s - t
    else:
        if t == 0:
            ret = GlobalCputime(cputime())
            for s in sage.interfaces.quit.expect_objects:
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
            for s in sage.interfaces.quit.expect_objects:
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

    EXAMPLE:

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

        EXAMPLE::

            sage: from sage.misc.misc import GlobalCputime
            sage: ct = GlobalCputime(0.0); ct
            0.0...
        """
        self.total = t
        self.local = t
        self.interfaces = {}

    def __repr__(self):
        """
        EXAMPLE::

            sage: cputime(subprocesses=True) # indirect doctest, output random
            0.2347431
        """
        return str(self.total)

    def __add__(self, other):
        """
        EXAMPLE::

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
        EXAMPLE::

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
        EXAMPLE::

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


#################################################################
# simple verbosity system
#################################################################
LEVEL=0  # default

verbose_files = []

def verbose(mesg="", t=0, level=1, caller_name=None):
    """
    Print a message if the current verbosity is at least level.

    INPUT:


    -  ``mesg`` - str, a message to print

    -  ``t`` - int, optional, if included, will also print
       cputime(t), - which is the time since time t. Thus t should have
       been obtained with t=cputime()

    -  ``level`` - int, (default: 1) the verbosity level of
       what we are printing

    -  ``caller_name`` - string (default: None), the name
       of the calling function; in most cases Python can deduce this, so
       it need not be provided.


    OUTPUT: possibly prints a message to stdout; also returns
    cputime()

    EXAMPLE::

        sage: set_verbose(1)
        sage: t = cputime()
        sage: t = verbose("This is Sage.", t, level=1, caller_name="william")       # not tested
        VERBOSE1 (william): This is Sage. (time = 0.0)
        sage: set_verbose(0)
    """
    if level>LEVEL:
        return cputime()

    frame = sys._getframe(1).f_code
    file_name = frame.co_filename
    lineno = frame.co_firstlineno
    if 'all' in verbose_files or level<=0:
        show = True
    else:
        show = False
        for X in verbose_files:
            if file_name.find(X) != -1:
                show = True
                break

    if not show:
        return cputime()

    if t != 0 and mesg=="":
        mesg = "Finished."

    # see recipe 14.7 in Python Cookbook
    if caller_name is None:
        caller_name = frame.co_name
        if caller_name == "?: ":
            caller_name = ""
    short_file_name = os.path.split(frame.co_filename)[1]
    if '<' in short_file_name and '>' in short_file_name:
        s = "verbose %s (%s) %s"%(level, caller_name, mesg)
    else:
        s = "verbose %s (%s: %s, %s) %s"%(level, lineno, short_file_name, caller_name, mesg)
    if t!=0:
        s = s + " (time = %s)"%cputime(t)
    print(s)
    sys.stdout.flush()
    return cputime()

def todo(mesg=""):
    from sage.misc.superseded import deprecation
    deprecation(21926, "todo() is deprecated")
    caller_name = sys._getframe(1).f_code.co_name
    raise NotImplementedError("{}: todo -- {}".format(caller_name, mesg))

def set_verbose(level, files='all'):
    """
    Set the global Sage verbosity level.

    INPUT:

    - ``level`` - an integer between 0 and 2, inclusive.

    - ``files`` (default: 'all'): list of files to make verbose, or
       'all' to make ALL files verbose (the default).

    OUTPUT: changes the state of the verbosity flag and possibly
    appends to the list of files that are verbose.

    EXAMPLES::

        sage: set_verbose(2)
        sage: verbose("This is Sage.", level=1)  # not tested
        VERBOSE1 (?): This is Sage.
        sage: verbose("This is Sage.", level=2)  # not tested
        VERBOSE2 (?): This is Sage.
        sage: verbose("This is Sage.", level=3)  # not tested
        [no output]
        sage: set_verbose(0)
    """
    if isinstance(level, str):
        set_verbose_files([level])
    global LEVEL
    LEVEL = level
    if isinstance(files, str):
        files = [files]
    set_verbose_files(files)

def set_verbose_files(file_name):
    """

    """
    if not isinstance(file_name, list):
        file_name = [file_name]
    global verbose_files
    verbose_files = file_name

def get_verbose_files():
    """

    """
    return verbose_files

def unset_verbose_files(file_name):
    """

    """
    if not isinstance(file_name, list):
        file_name = [file_name]
    for X in file_name:
        verbose_files.remove(X)


def get_verbose():
    """
    Return the global Sage verbosity level.

    INPUT: int level: an integer between 0 and 2, inclusive.

    OUTPUT: changes the state of the verbosity flag.

    EXAMPLES::

        sage: get_verbose()
        0
        sage: set_verbose(2)
        sage: get_verbose()
        2
        sage: set_verbose(0)
    """
    global LEVEL
    return LEVEL



def generic_cmp(x,y):
    """
    Compare x and y and return -1, 0, or 1.

    This is similar to x.__cmp__(y), but works even in some cases
    when a .__cmp__ method isn't defined.
    """
    from sage.misc.superseded import deprecation
    deprecation(21926, "generic_cmp() is deprecated")
    if x<y:
        return -1
    elif x==y:
        return 0
    return 1

def cmp_props(left, right, props):
    for a in props:
        c = cmp(left.__getattribute__(a)(), right.__getattribute__(a)())
        if c: return c
    return 0

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
        [1, 2, 3, 4, 5, 6]
        sage: union([1,2,3,4,5,6], [5,6]) == answer
        True
        sage: union((1,2,3,4,5,6), [5,6]) == answer
        True
        sage: union((1,2,3,4,5,6), set([5,6])) == answer
        True
    """
    if y is None:
        return list(set(x))
    return list(set(x).union(y))

def uniq(x):
    """
    Return the sublist of all elements in the list x that is sorted and
    is such that the entries in the sublist are unique.

    EXAMPLES::

        sage: v = uniq([1,1,8,-5,3,-5,'a','x','a'])
        sage: v            # potentially random ordering of output
        ['a', 'x', -5, 1, 3, 8]
        sage: set(v) == set(['a', 'x', -5, 1, 3, 8])
        True
    """
    v = sorted(set(x))
    return v


def coeff_repr(c, is_latex=False):
    if not is_latex:
        try:
            return c._coeff_repr()
        except AttributeError:
            pass
    if isinstance(c, (int, long, float)):
        return str(c)
    if is_latex and hasattr(c, '_latex_'):
        s = c._latex_()
    else:
        s = str(c).replace(' ','')
    if s.find("+") != -1 or s.find("-") != -1:
        if is_latex:
            return "\\left(%s\\right)"%s
        else:
            return "(%s)"%s
    return s

def repr_lincomb(terms, is_latex=False, scalar_mult="*", strip_one=False, repr_monomial = None, latex_scalar_mult = None):
    """
    Compute a string representation of a linear combination of some
    formal symbols.

    INPUT:

    - ``terms`` -- list of terms, as pairs (support, coefficient)
    - ``is_latex`` -- whether to produce latex (default: ``False``)
    - ``scalar_mult`` -- string representing the multiplication (default:``'*'``)
    - ``latex_scalar_mult`` -- latex string representing the multiplication
      (default: ``''`` if ``scalar_mult`` is ``'*'``; otherwise ``scalar_mult``)
    - ``coeffs`` -- for backward compatibility

    OUTPUT:

    -  ``str`` - a string

    EXAMPLES::

        sage: repr_lincomb([('a',1), ('b',-2), ('c',3)])
        'a - 2*b + 3*c'
        sage: repr_lincomb([('a',0), ('b',-2), ('c',3)])
        '-2*b + 3*c'
        sage: repr_lincomb([('a',0), ('b',2), ('c',3)])
        '2*b + 3*c'
        sage: repr_lincomb([('a',1), ('b',0), ('c',3)])
        'a + 3*c'
        sage: repr_lincomb([('a',-1), ('b','2+3*x'), ('c',3)])
        '-a + (2+3*x)*b + 3*c'
        sage: repr_lincomb([('a', '1+x^2'), ('b', '2+3*x'), ('c', 3)])
        '(1+x^2)*a + (2+3*x)*b + 3*c'
        sage: repr_lincomb([('a', '1+x^2'), ('b', '-2+3*x'), ('c', 3)])
        '(1+x^2)*a + (-2+3*x)*b + 3*c'
        sage: repr_lincomb([('a', 1), ('b', -2), ('c', -3)])
        'a - 2*b - 3*c'
        sage: t = PolynomialRing(RationalField(),'t').gen()
        sage: repr_lincomb([('a', -t), ('s', t - 2), ('', t^2 + 2)])
        '-t*a + (t-2)*s + (t^2+2)'

    Examples for ``scalar_mult``::

        sage: repr_lincomb([('a',1), ('b',2), ('c',3)], scalar_mult='*')
        'a + 2*b + 3*c'
        sage: repr_lincomb([('a',2), ('b',0), ('c',-3)], scalar_mult='**')
        '2**a - 3**c'
        sage: repr_lincomb([('a',-1), ('b',2), ('c',3)], scalar_mult='**')
        '-a + 2**b + 3**c'

    Examples for ``scalar_mult`` and ``is_latex``::

        sage: repr_lincomb([('a',-1), ('b',2), ('c',3)], is_latex=True)
        '-a + 2b + 3c'
        sage: repr_lincomb([('a',-1), ('b',-1), ('c',3)], is_latex=True, scalar_mult='*')
        '-a - b + 3c'
        sage: repr_lincomb([('a',-1), ('b',2), ('c',-3)], is_latex=True, scalar_mult='**')
        '-a + 2**b - 3**c'
        sage: repr_lincomb([('a',-2), ('b',-1), ('c',-3)], is_latex=True, latex_scalar_mult='*')
        '-2*a - b - 3*c'

    Examples for ``strip_one``::

        sage: repr_lincomb([ ('a',1), (1,-2), ('3',3) ])
        'a - 2*1 + 3*3'
        sage: repr_lincomb([ ('a',-1), (1,1), ('3',3) ])
        '-a + 1 + 3*3'
        sage: repr_lincomb([ ('a',1), (1,-2), ('3',3) ], strip_one = True)
        'a - 2 + 3*3'
        sage: repr_lincomb([ ('a',-1), (1,1), ('3',3) ], strip_one = True)
        '-a + 1 + 3*3'
        sage: repr_lincomb([ ('a',1), (1,-1), ('3',3) ], strip_one = True)
        'a - 1 + 3*3'

    Examples for ``repr_monomial``::

        sage: repr_lincomb([('a',1), ('b',2), ('c',3)], repr_monomial = lambda s: s+"1")
        'a1 + 2*b1 + 3*c1'
    """
    # Setting scalar_mult: symbol used for scalar multiplication
    if is_latex:
        if latex_scalar_mult is not None:
            scalar_mult = latex_scalar_mult
        elif scalar_mult == "*":
            scalar_mult = ""

    if repr_monomial is None:
        if is_latex:
            repr_monomial = lambda monomial: monomial._latex_() if hasattr(monomial, '_latex_') else str(monomial)
        else:
            repr_monomial = str

    s = ""
    first = True

    if scalar_mult is None:
        scalar_mult = "" if is_latex else "*"

    for (monomial,c) in terms:
        if c != 0:
            coeff = coeff_repr(c)
            negative = False
            if len(coeff)>0 and coeff[0] == "-":
                negative = True
            try:
                if c < 0:
                    negative = True
            except NotImplementedError:
                # comparisons may not be implemented for some coefficients
                pass
            if negative:
                coeff = coeff_repr(-c, is_latex)
            else:
                coeff = coeff_repr(c, is_latex)
            if coeff == "1":
                coeff = ""
            if coeff != "0":
                if negative:
                    if first:
                        sign = "-" # add trailing space?
                    else:
                        sign = " - "
                else:
                    if first:
                        sign = ""
                    else:
                        sign= " + "
                b = repr_monomial(monomial)
                if len(b) > 0:
                    if  coeff != "":
                        if b =="1" and strip_one:
                            b = ""
                        else:
                            b = scalar_mult + b
                s += "%s%s%s"%(sign, coeff, b)
                first = False
    if first:
        return "0" # this can happen only if are only terms with coeff_repr(c) == "0"
    #elif s == "":
        #return "1" # is empty string representation invalid?
    else:
        return s


def strunc(s, n = 60):
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
        #return s[:n-4] + " ..."
    return s



def newton_method_sizes(N):
    r"""
    Returns a sequence of integers
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


def assert_attribute(x, attr, init=None):
    """
    If the object x has the attribute attr, do nothing. If not, set
    x.attr to init.
    """
    from sage.misc.superseded import deprecation
    deprecation(21926, "assert_attribute() is deprecated")
    if attr in x.__dict__: return
    if attr[:2] == "__":
        z = str(x.__class__).split("'")
        if len(z) > 1:
            z = z[1]
        else:
            z = z[0]
        attr = "_" + z[len(x.__module__)+1:] + attr
    x.__dict__[attr] = init


def compose(f, g):
    """
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


def self_compose(f, n):
    """
    Return the function `f` composed with itself `n` times.

    See :func:`nest()` if you want `f(f(...(f(x))...))` for
    known `x`.


    INPUT:
        - `f` -- a function of one variable
        - `n` -- a nonnegative integer

    OUTPUT:
        A function, the result of composing `f` with itself `n` times

    EXAMPLES::

        sage: def f(x): return x^2 + 1
        sage: g = self_compose(f, 3)
        doctest:... DeprecationWarning: self_compose() is deprecated, use nest() instead
        See http://trac.sagemath.org/21926 for details.
        sage: x = var('x')
        sage: g(x)
        ((x^2 + 1)^2 + 1)^2 + 1

    ::

        sage: def f(x): return x + 1
        sage: g = self_compose(f, 10000)
        sage: g(0)
        10000

    ::

        sage: x = var('x')
        sage: self_compose(sin, 0)(x)
        x

    """
    from sage.misc.superseded import deprecation
    deprecation(21926, "self_compose() is deprecated, use nest() instead")
    from sage.rings.all import Integer
    n = Integer(n)

    if n < 0:
        raise ValueError("n must be a nonnegative integer, not {}.".format(n))

    return lambda x: nest(f, n, x)


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
    from sage.rings.all import Integer
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
    """
    Implements Matlab-style backslash operator for solving systems::

        A \\ b

    The preparser converts this to multiplications using
    ``BackslashOperator()``.

    EXAMPLES::

        sage: preparse("A \ matrix(QQ,2,1,[1/3,'2/3'])")
        "A  * BackslashOperator() * matrix(QQ,Integer(2),Integer(1),[Integer(1)/Integer(3),'2/3'])"
        sage: preparse("A \ matrix(QQ,2,1,[1/3,2*3])")
        'A  * BackslashOperator() * matrix(QQ,Integer(2),Integer(1),[Integer(1)/Integer(3),Integer(2)*Integer(3)])'
        sage: preparse("A \ B + C")
        'A  * BackslashOperator() * B + C'
        sage: preparse("A \ eval('C+D')")
        "A  * BackslashOperator() * eval('C+D')"
        sage: preparse("A \ x / 5")
        'A  * BackslashOperator() * x / Integer(5)'
        sage: preparse("A^3 \ b")
        'A**Integer(3)  * BackslashOperator() * b'
    """
    def __rmul__(self, left):
        """
        EXAMPLES::

            sage: A = random_matrix(ZZ, 4)
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
        """
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
def is_iterator(it):
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
        ....:    def next(self):
        ....:        self.n -= 1
        ....:        if self.n == 0: raise StopIteration
        ....:        return self.n
        sage: x = wrong()
        sage: is_iterator(x)
        False
        sage: list(x)
        Traceback (most recent call last):
        ...
        TypeError: iteration over non-sequence

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

        sage: S = [1,7,3,4,18]
        sage: random_sublist(S, 0.5)
        [1, 3, 4]
        sage: random_sublist(S, 0.5)
        [1, 3]
    """
    return [a for a in X if random.random() <= s]


def some_tuples(elements, repeat, bound):
    r"""
    Return an iterator over at most ``bound`` number of ``repeat``-tuples of
    ``elements``.

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

    .. TODO::

        Currently, this only return an iterator over the first element of the
        Cartesian product. It would be smarter to return something more
        "random like" as it is used in tests. However, this should remain
        deterministic.
    """
    from itertools import islice, product
    return islice(product(elements, repeat=repeat), bound)

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
        <generator object powerset at 0x...>
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
    for x in X:
        pairs.append((2**len(pairs),x))
        for w in range(2**(len(pairs)-1), 2**(len(pairs))):
            yield [x for m, x in pairs if m & w]

subsets = powerset

#################################################################
# Type checking
#################################################################
def typecheck(x, C, var="x"):
    """
    Check that x is of instance C. If not raise a TypeError with an
    error message.
    """
    from sage.misc.superseded import deprecation
    deprecation(21926, "typecheck is deprecated, use isinstance instead")
    if not isinstance(x, C):
        raise TypeError("{} (={}) must be of type {}.".format(var, x, C))

#################################################################
# This will likely eventually be useful.
#################################################################

# From the Python Cookbook Ver 2, Recipe 20.4
class cached_attribute(object):
    """
    Computes attribute value and caches it in the instance.
    """
    def __init__(self, method, name=None):
        from sage.misc.superseded import deprecation
        deprecation(21926, "cached_attribute is deprecated")
        # record the unbound-method and the name
        self.method = method
        self.name = name or method.__name__
    def __get__(self, inst, cls):
        if inst is None:
            # instance attribute accessed on class, return self
            return self
        # compute, cache and return the instance's attribute value
        result = self.method(inst)
        setattr(inst, self.name, result)
        return result

class lazy_prop(object):
    def __init__(self, calculate_function):
        from sage.misc.superseded import deprecation
        deprecation(21926, "lazy_prop is deprecated")
        self._calculate = calculate_function
        self.__doc__ = calculate_function.__doc__

    def __call__(self, obj, _=None):
        if obj is None:
            return self
        value = self._calculate(obj)
        setattr(obj, self._calculate.__name__, value)
        return value

def prop(f):
    from sage.misc.superseded import deprecation
    deprecation(21926, "prop() is deprecated")
    return property(f, None, None, f.__doc__)


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
        if P(x): return True, x
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
        if not P(x): return False, x
    return True, None

#################################################################
# which source file?
#################################################################
import inspect
def sourcefile(object):
    """
    Work out which source or compiled file an object was defined in.
    """
    from sage.misc.superseded import deprecation
    deprecation(21926, "sourcefile(x) is deprecated, use inspect.getfile(x) instead")
    return inspect.getfile(object)


#################################################################
# debug tracing
#################################################################
import pdb
set_trace = pdb.set_trace



#################################################################
# Word wrap lines
#################################################################
def word_wrap(s, ncols=85):
    t = []
    if ncols == 0:
        return s
    for x in s.split('\n'):
        if len(x) == 0 or x.lstrip()[:5] == 'sage:':
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
            k=0
            while k < len(x) and x[k] == ' ':
                k += 1
            x = x[k:]
        t.append(x)
    return '\n'.join(t)


def getitem(v, n):
    r"""
    Variant of getitem that coerces to an int if a TypeError is
    raised.

    (This is not needed anymore - classes should define an
    __index__ method.)

    Thus, e.g., ``getitem(v,n)`` will work even if
    `v` is a Python list and `n` is a Sage integer.

    EXAMPLES::

        sage: v = [1,2,3]

    The following used to fail in Sage <= 1.3.7. Now it works fine::

        sage: v[ZZ(1)]
        2

    This always worked.

    ::

        sage: getitem(v, ZZ(1))
        doctest:... DeprecationWarning: getitem(v, n) is deprecated, use v[n] instead
        See http://trac.sagemath.org/21926 for details.
        2
    """
    from sage.misc.superseded import deprecation
    deprecation(21926, "getitem(v, n) is deprecated, use v[n] instead")
    try:
        return v[n]
    except TypeError:
        return v[int(n)]

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
    return "0"*(size-len(str(s))) + str(s)

import sage.server.support

def embedded():
    """
    Return True if this copy of Sage is running embedded in the Sage
    notebook.

    EXAMPLES::

        sage: sage.misc.misc.embedded()    # output True if in the notebook
        False
    """
    return sage.server.support.EMBEDDED_MODE


#############################################
# Operators
#############################################
class AttrCallObject(object):
    def __init__(self, name, args, kwds):
        """
        TESTS::

            sage: f = attrcall('core', 3); f
            *.core(3)
            sage: TestSuite(f).run()
        """
        self.name = name
        self.args = args
        self.kwds = kwds

    def __call__(self, x, *args):
        """
        Gets the ``self.name`` method from ``x``, calls it with
        ``self.args`` and ``args`` as positional parameters and
        ``self.kwds`` as keyword parameters, and returns the result.

        EXAMPLES::

            sage: core = attrcall('core', 3)
            sage: core(Partition([4,2]))
            [4, 2]

            sage: series = attrcall('series', x)
            sage: series(sin(x), 4)
            1*x + (-1/6)*x^3 + Order(x^4)
        """
        return getattr(x, self.name)(*(self.args+args), **self.kwds)

    def __repr__(self):
        """
        Returns a string representation of this object. The star in the
        output represents the object passed into self.

        EXAMPLES::

            sage: attrcall('core', 3)
            *.core(3)
            sage: attrcall('hooks', flatten=True)
            *.hooks(flatten=True)
            sage: attrcall('hooks', 3, flatten=True)
            *.hooks(3, flatten=True)
        """
        s =  "*.%s(%s"%(self.name, ", ".join(map(repr, self.args)))
        if self.kwds:
            if len(self.args) > 0:
                s += ", "
            s += ", ".join("%s=%s"%keyvalue for keyvalue in self.kwds.items())
        s += ")"
        return s

    def __eq__(self, other):
        """
        Equality testing

        EXAMPLES::

            sage: attrcall('core', 3, flatten = True) == attrcall('core', 3, flatten = True)
            True
            sage: attrcall('core', 2) == attrcall('core', 3)
            False
            sage: attrcall('core', 2) == 1
            False
        """
        return self.__class__ == other.__class__ and self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Equality testing

        EXAMPLES::

            sage: attrcall('core', 3, flatten = True) != attrcall('core', 3, flatten = True)
            False
            sage: attrcall('core', 2) != attrcall('core', 3)
            True
            sage: attrcall('core', 2) != 1
            True
        """
        return not self == other

    def __hash__(self):
        """
        Hash value

        This method tries to ensure that, when two ``attrcall``
        objects are equal, they have the same hash value.

        .. warning:: dicts are not hashable, so we instead hash their
        items; however the order of those items might differ. The
        proper fix would be to use a frozen dict for ``kwds``, when
        frozen dicts will be available in Python.

        EXAMPLES::

            sage: x = attrcall('core', 3, flatten = True, blah = 1)
            sage: hash(x)       # random # indirect doctest
            210434060
            sage: type(hash(x))
            <type 'int'>
            sage: y = attrcall('core', 3, blah = 1, flatten = True)
            sage: hash(y) == hash(x)
            True
            sage: y = attrcall('core', 3, flatten = True, blah = 2)
            sage: hash(y) != hash(x)
            True
            sage: hash(attrcall('core', 2)) != hash(attrcall('core', 3))
            True
            sage: hash(attrcall('core', 2)) != hash(1)
            True

        Note: a missing ``__hash__`` method here used to break the
        unique representation of parents taking ``attrcall`` objects
        as input; see :trac:`8911`.
        """
        return hash((self.args, tuple(self.kwds.items())))

def attrcall(name, *args, **kwds):
    """
    Returns a callable which takes in an object, gets the method named
    name from that object, and calls it with the specified arguments
    and keywords.

    INPUT:

     -  ``name`` - a string of the name of the method you
        want to call

     -  ``args, kwds`` - arguments and keywords to be passed
        to the method

    EXAMPLES::

        sage: f = attrcall('core', 3); f
        *.core(3)
        sage: [f(p) for p in Partitions(5)]
        [[2], [1, 1], [1, 1], [3, 1, 1], [2], [2], [1, 1]]
    """
    return AttrCallObject(name, args, kwds)

def call_method(obj, name, *args, **kwds):
    """
    Call the method ``name`` on ``obj``.

    This has to exist somewhere in Python!!!

    .. SEEALSO:: :func:`operator.methodcaller` :func:`attrcal`

    EXAMPLES::

        sage: from sage.misc.misc import call_method
        sage: call_method(1, "__add__", 2)
        3
    """
    return getattr(obj, name)(*args, **kwds)

def is_in_string(line, pos):
    r"""
    Returns True if the character at position pos in line occurs
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
        if line[i-1:i] != '\\' or line[i-2:i] == '\\\\':
            if line[i:i+3] in ['"""', "'''"]:
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


def inject_variable(name, value):
    """
    Inject a variable into the main global namespace.

    INPUT:

    - ``name``  -- a string
    - ``value`` -- anything

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

    Use with care!
    """
    assert isinstance(name, str)
    # Using globals() does not work, even in Cython, because
    # inject_variable is called not only from the interpreter, but
    # also from functions in various modules.
    G = get_main_globals()
    if name in G:
        warn("redefining global value `%s`"%name, RuntimeWarning, stacklevel = 2)
    G[name] = value


def inject_variable_test(name, value, depth):
    """
    A function for testing deep calls to inject_variable

    TESTS::

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
