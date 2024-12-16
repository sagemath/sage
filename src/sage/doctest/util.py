# sage_setup: distribution = sagemath-repl
"""
Utility functions

This module contains various utility functions and classes used in doctesting.

AUTHORS:

- David Roe (2012-03-27) -- initial version, based on Robert Bradshaw's code.
"""

# ****************************************************************************
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                     2012 Robert Bradshaw <robertwb@gmail.com>
#                     2012 William Stein <wstein@gmail.com>
#                     2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2014 Volker Braun
#                     2017 Frédéric Chapoton
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from time import time as walltime
from os import sysconf, times
from contextlib import contextmanager
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt


def count_noun(number, noun, plural=None, pad_number=False, pad_noun=False):
    """
    EXAMPLES::

        sage: from sage.doctest.util import count_noun
        sage: count_noun(1, "apple")
        '1 apple'
        sage: count_noun(1, "apple", pad_noun=True)
        '1 apple '
        sage: count_noun(1, "apple", pad_number=3)
        '  1 apple'
        sage: count_noun(2, "orange")
        '2 oranges'
        sage: count_noun(3, "peach", "peaches")
        '3 peaches'
        sage: count_noun(1, "peach", plural='peaches', pad_noun=True)
        '1 peach  '
    """
    if plural is None:
        plural = noun + "s"
    if pad_noun:
        # We assume that the plural is never shorter than the noun....
        pad_noun = " " * (len(plural) - len(noun))
    else:
        pad_noun = ""
    if pad_number:
        number_str = ("%%%sd" % pad_number) % number
    else:
        number_str = "%d" % number
    if number == 1:
        return "%s %s%s" % (number_str, noun, pad_noun)
    else:
        return "%s %s" % (number_str, plural)


def dict_difference(self, other):
    """
    Return a dict with all key-value pairs occurring in ``self`` but not
    in ``other``.

    EXAMPLES::

        sage: from sage.doctest.util import dict_difference
        sage: d1 = {1: 'a', 2: 'b', 3: 'c'}
        sage: d2 = {1: 'a', 2: 'x', 4: 'c'}
        sage: dict_difference(d2, d1)
        {2: 'x', 4: 'c'}

    ::

        sage: from sage.doctest.control import DocTestDefaults
        sage: D1 = DocTestDefaults()
        sage: D2 = DocTestDefaults(foobar='hello', timeout=100)
        sage: dict_difference(D2.__dict__, D1.__dict__)
        {'foobar': 'hello', 'timeout': 100}
    """
    D = {}
    for k, v in self.items():
        try:
            if other[k] == v:
                continue
        except KeyError:
            pass
        D[k] = v
    return D


class Timer:
    """
    A simple timer.

    EXAMPLES::

        sage: from sage.doctest.util import Timer
        sage: Timer()
        {}
        sage: TestSuite(Timer()).run()
    """

    def _proc_stat_cpu_seconds(self, path):
        r"""
        Parse a "stat" file from the ``/proc`` filesystem to get
        the cputime of a process.

        This also includes the times for child processes, but only
        those that have already terminated and for which ``wait()``
        was called. It is important to note that pexpect processes DO
        NOT fall into that category.

        The document ``Documentation/filesystems/proc.rst`` within the
        Linux kernel source tree defines a "stat" file.

        INPUT:

        - ``path`` -- string; the path to a "stat" file on the ``/proc``
          filesystem, typically "/proc/<pid>/stat", from which we will
          read cputime information

        OUTPUT:

        A nonnegative float representing the number of cpu-seconds
        used by the process associated with ``path``. An ``OSError`` is
        raised if anything goes wrong, which typically happens on
        platforms that don't store this information under ``/proc``.

        TESTS:

        About all we can say for certain is that this will return a
        nonnegative float or raise an ``OSError``::

            sage: from sage.doctest.util import Timer
            sage: cputime = float(0.0)
            sage: path = "/proc/1/stat"
            sage: try:
            ....:     cputime = Timer()._proc_stat_cpu_seconds(path)
            ....: except OSError:
            ....:     pass
            sage: cputime >= 0.0
            True
            sage: isinstance(cputime, float)
            True

        We can force an ``OSError`` with an invalid PID::

            sage: from sage.doctest.util import Timer
            sage: path = "/proc/-1/stat"
            sage: cputime = Timer()._proc_stat_cpu_seconds(path)
            Traceback (most recent call last):
            ...
            OSError: unable to access /proc/-1/stat

        Or with an unparseable file (wrong number of fields, non-float
        fields, et cetera)::

            sage: from tempfile import NamedTemporaryFile
            sage: from os import unlink
            sage: from sage.doctest.util import Timer
            sage: with NamedTemporaryFile(delete=False, mode="w") as f:
            ....:     _ = f.write("1 2 3 4 5")
            sage: cputime = Timer()._proc_stat_cpu_seconds(f.name)
            Traceback (most recent call last):
            ...
            OSError: unable to parse ...
            sage: os.unlink(f.name)
            sage: with NamedTemporaryFile(delete=False, mode="w") as f:
            ....:     _ = f.write("1 2 3 4 5 6 7 8 9 10 11 12 w x y z 17")
            sage: cputime = Timer()._proc_stat_cpu_seconds(f.name)
            Traceback (most recent call last):
            ...
            OSError: unable to parse ...
            sage: os.unlink(f.name)

        """
        try:
            with open(path, "r") as statfile:
                stats = statfile.read().split()
        except (FileNotFoundError, PermissionError) as e:
            # FileNotFoundError: bad PID, or no /proc support
            # PermissionError: can't read the stat file
            raise OSError(f"unable to access {path}") from e

        if len(stats) < 17:
            raise OSError(f"unable to parse {path}")

        try:
            # These fields used to be documented in the proc(5) man
            # page, but are now most easily found in the Linux kernel
            # documentation (Documentation/filesystems/proc.rst). The
            # intent is to sum the user- and kernel-mode "jiffies" for
            # both the given process and its children.
            cputicks = sum( float(s) for s in stats[13:17] )
        except (ArithmeticError, TypeError, ValueError) as e:
            # ArithmeticError: unexpected (non-numeric?) values in fields
            # TypeError/ValueError: fields can't be converted to float
            raise OSError(f"unable to parse {path}") from e

        try:
            hertz = sysconf("SC_CLK_TCK")
        except (ValueError) as e:
            # ValueError: SC_CLK_TCK doesn't exist
            raise OSError("SC_CLK_TCK sysconf not found") from e

        if hertz <= 0:
            # The python documentation for os.sysconf() says, "If the
            # configuration value specified by name isn’t defined, -1
            # is returned." Having tried this with a junk value, I
            # don't believe it: I got a ValueError that was handled
            # above. Nevertheless, we play it safe here and turn a -1
            # into an OSError. We check for zero, too, because we're
            # about to divide by it.
            raise OSError("SC_CLK_TCK sysconf is nonpositive")

        return (cputicks / hertz)

    def _quick_cputime(self, expect_objects):
        r"""
        A fast replacement for ``sage.misc.timing.cputime``.

        This is a "reliable" replacement (on Linux/BSD) that takes
        subprocesses (particularly pexpect interfaces) into
        account. The ``cputime()`` function from the ``misc`` module
        can be passed ``subprocesses=True``, but this has a few
        faults; mainly that it relies on each pexpect interface to
        implement its own ``cputime()`` function. And most of our
        pexpect interfaces either don't implement one, or implement
        one in a way that requires the subprocess (being pexpected) to
        be in perfect working condition -- that will often not be the
        case at the end of a doctest line.

        INPUT:

        - ``expect_objects`` -- list; a list of
          :class:`sage.interfaces.expect.Expect` instances whose CPU
          times will be included in the total

        OUTPUT:

        A float measuring the cputime in seconds of the sage process
        and all its subprocesses.

        TESTS:

        About all we can say for certain is that this will return a
        nonnegative float::

            sage: from sage.doctest.util import Timer
            sage: from sage.interfaces.quit import expect_objects
            sage: cputime = Timer()._quick_cputime(expect_objects)
            sage: cputime >= 0.0
            True
            sage: isinstance(cputime, float)
            True

        If an error occurs in :meth:`_proc_stat_cpu_seconds`, this
        function should still return a valid answer, albeit one that
        is missing timing information for the PID that failed::

            sage: class FakeExpect:
            ....:     def __call__(self):
            ....:         return self
            ....:     def is_running(self):
            ....:         return True
            ....:     def pid(self):
            ....:         return -1
            sage: e = FakeExpect()
            sage: from sage.doctest.util import Timer
            sage: cputime = Timer()._quick_cputime([e])
            sage: cputime >= 0.0
            True
            sage: isinstance(cputime, float)
            True
        """
        # Start by using os.times() to get the cputime for sage itself
        # and any subprocesses that have been wait()ed for and that
        # have terminated.
        cputime = sum( times()[:4] )

        # Now try to get the times for any pexpect interfaces, since
        # they do not fall into the category above.
        for s in expect_objects:
            S = s()
            if S and S.is_running():
                try:
                    # This will fail anywhere but linux/BSD, but
                    # there's no good cross-platform way to get the
                    # cputimes from pexpect interfaces without totally
                    # mucking up the doctests.
                    path = f"/proc/{S.pid()}/stat"
                    cputime += self._proc_stat_cpu_seconds(path)
                except OSError:
                    # If we're on macOS, we can fall back to using
                    # psutil, but only if it's installed. It's usually
                    # installed as a transitive dependency (ipython
                    # needs it), but it isn't explicitly listed as
                    # a dependency of sagelib.
                    try:
                        from psutil import (NoSuchProcess,
                                            Process,
                                            ZombieProcess)
                        try:
                            cputime += sum(Process(S.pid()).cpu_times()[0:2])
                        except (ValueError, NoSuchProcess, ZombieProcess):
                            # ValueError: invalid (e.g. negative) PID
                            # NoSuchProcess: it's gone
                            # ZombieProcess: PID refers to a zombie
                            pass
                    except ImportError:
                        pass

        return cputime

    def start(self):
        """
        Start the timer.

        Can be called multiple times to reset the timer.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: Timer().start()
            {'cputime': ..., 'walltime': ...}
        """
        from sage.interfaces.quit import expect_objects
        self.cputime = self._quick_cputime(expect_objects)
        self.walltime = walltime()
        return self

    def stop(self):
        """
        Stops the timer, recording the time that has passed since it
        was started.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: import time
            sage: timer = Timer().start()
            sage: time.sleep(float(0.5))
            sage: timer.stop()
            {'cputime': ..., 'walltime': ...}
        """
        from sage.interfaces.quit import expect_objects
        self.cputime = self._quick_cputime(expect_objects) - self.cputime
        self.walltime = walltime() - self.walltime
        return self

    def annotate(self, object):
        """
        Annotates the given object with the cputime and walltime
        stored in this timer.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: Timer().start().annotate(EllipticCurve)
            sage: EllipticCurve.cputime # random
            2.817255
            sage: EllipticCurve.walltime # random
            1332649288.410404
        """
        object.cputime = self.cputime
        object.walltime = self.walltime

    def __repr__(self):
        """
        String representation is from the dictionary.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: repr(Timer().start()) # indirect doctest
            "{'cputime': ..., 'walltime': ...}"
        """
        return str(self)

    def __str__(self):
        """
        String representation is from the dictionary.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: str(Timer().start()) # indirect doctest
            "{'cputime': ..., 'walltime': ...}"
        """
        return str(self.__dict__)

    def __eq__(self, other):
        """
        Comparison.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: Timer() == Timer()
            True
            sage: t = Timer().start()
            sage: loads(dumps(t)) == t
            True
        """
        if not isinstance(other, Timer):
            return False
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Test for non-equality.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: Timer() == Timer()
            True
            sage: t = Timer().start()
            sage: loads(dumps(t)) != t
            False
        """
        return not (self == other)


# Inheritance rather then delegation as globals() must be a dict
class RecordingDict(dict):
    """
    This dictionary is used for tracking the dependencies of an example.

    This feature allows examples in different doctests to be grouped
    for better timing data.  It's obtained by recording whenever
    anything is set or retrieved from this dictionary.

    EXAMPLES::

        sage: from sage.doctest.util import RecordingDict
        sage: D = RecordingDict(test=17)
        sage: D.got
        set()
        sage: D['test']
        17
        sage: D.got
        {'test'}
        sage: D.set
        set()
        sage: D['a'] = 1
        sage: D['a']
        1
        sage: D.set
        {'a'}
        sage: D.got
        {'test'}

    TESTS::

        sage: TestSuite(D).run()
    """
    def __init__(self, *args, **kwds):
        """
        Initialization arguments are the same as for a normal dictionary.

        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D.got
            set()
        """
        dict.__init__(self, *args, **kwds)
        self.start()

    def start(self):
        """
        We track which variables have been set or retrieved.
        This function initializes these lists to be empty.

        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D.set
            set()
            sage: D['a'] = 4
            sage: D.set
            {'a'}
            sage: D.start(); D.set
            set()
        """
        self.set = set()
        self.got = set()

    def __getitem__(self, name):
        """
        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D['a'] = 4
            sage: D.got
            set()
            sage: D['a'] # indirect doctest
            4
            sage: D.got
            set()
            sage: D['d']
            42
            sage: D.got
            {'d'}
        """
        if name not in self.set:
            self.got.add(name)
        return dict.__getitem__(self, name)

    def __setitem__(self, name, value):
        """
        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D['a'] = 4 # indirect doctest
            sage: D.set
            {'a'}
        """
        self.set.add(name)
        dict.__setitem__(self, name, value)

    def __delitem__(self, name):
        """
        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: del D['d'] # indirect doctest
            sage: D.set
            {'d'}
        """
        self.set.add(name)
        dict.__delitem__(self, name)

    def get(self, name, default=None):
        """
        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D.get('d')
            42
            sage: D.got
            {'d'}
            sage: D.get('not_here')
            sage: sorted(list(D.got))
            ['d', 'not_here']
        """
        if name not in self.set:
            self.got.add(name)
        return dict.get(self, name, default)

    def copy(self):
        """
        Note that set and got are not copied.

        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D['a'] = 4
            sage: D.set
            {'a'}
            sage: E = D.copy()
            sage: E.set
            set()
            sage: sorted(E.keys())
            ['a', 'd']
        """
        return RecordingDict(dict.copy(self))

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D['a'] = 4
            sage: D.get('not_here')
            sage: E = loads(dumps(D))
            sage: E.got
            {'not_here'}
        """
        return make_recording_dict, (dict(self), self.set, self.got)


def make_recording_dict(D, st, gt):
    """
    Auxiliary function for pickling.

    EXAMPLES::

        sage: from sage.doctest.util import make_recording_dict
        sage: D = make_recording_dict({'a':4,'d':42},set([]),set(['not_here']))
        sage: sorted(D.items())
        [('a', 4), ('d', 42)]
        sage: D.got
        {'not_here'}
    """
    ans = RecordingDict(D)
    ans.set = st
    ans.got = gt
    return ans


class NestedName:
    """
    Class used to construct fully qualified names based on indentation level.

    EXAMPLES::

        sage: from sage.doctest.util import NestedName
        sage: qname = NestedName('sage.categories.algebras')
        sage: qname[0] = 'Algebras'; qname
        sage.categories.algebras.Algebras
        sage: qname[4] = '__contains__'; qname
        sage.categories.algebras.Algebras.__contains__
        sage: qname[4] = 'ParentMethods'
        sage: qname[8] = 'from_base_ring'; qname
        sage.categories.algebras.Algebras.ParentMethods.from_base_ring

    TESTS::

        sage: TestSuite(qname).run()
    """
    def __init__(self, base):
        """
        INPUT:

        - ``base`` -- string; the name of the module

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname
            sage.categories.algebras
        """
        self.all = [base]

    def __setitem__(self, index, value):
        """
        Set the value at a given indentation level.

        INPUT:

        - ``index`` -- positive integer; the indentation level (often a
          multiple of 4, but not necessarily)
        - ``value`` -- string; the name of the class or function at that
          indentation level

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname[1] = 'Algebras' # indirect doctest
            sage: qname
            sage.categories.algebras.Algebras
            sage: qname.all
            ['sage.categories.algebras', None, 'Algebras']
        """
        if index < 0:
            raise ValueError
        while len(self.all) <= index:
            self.all.append(None)
        self.all[index+1:] = [value]

    def __str__(self):
        """
        Return a .-separated string giving the full name.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname[1] = 'Algebras'
            sage: qname[44] = 'at_the_end_of_the_universe'
            sage: str(qname) # indirect doctest
            'sage.categories.algebras.Algebras.at_the_end_of_the_universe'
        """
        return repr(self)

    def __repr__(self):
        """
        Return a .-separated string giving the full name.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname[1] = 'Algebras'
            sage: qname[44] = 'at_the_end_of_the_universe'
            sage: print(qname) # indirect doctest
            sage.categories.algebras.Algebras.at_the_end_of_the_universe
        """
        return '.'.join(a for a in self.all if a is not None)

    def __eq__(self, other):
        """
        Comparison is just comparison of the underlying lists.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname2 = NestedName('sage.categories.algebras')
            sage: qname == qname2
            True
            sage: qname[0] = 'Algebras'
            sage: qname2[2] = 'Algebras'
            sage: repr(qname) == repr(qname2)
            True
            sage: qname == qname2
            False
        """
        if not isinstance(other, NestedName):
            return False
        return self.all == other.all

    def __ne__(self, other):
        """
        Test for non-equality.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname2 = NestedName('sage.categories.algebras')
            sage: qname != qname2
            False
            sage: qname[0] = 'Algebras'
            sage: qname2[2] = 'Algebras'
            sage: repr(qname) == repr(qname2)
            True
            sage: qname != qname2
            True
        """
        return not (self == other)


@contextmanager
def ensure_interruptible_after(seconds: float, max_wait_after_interrupt: float = 0.2, inaccuracy_tolerance: float = 0.1):
    """
    Helper function for doctesting to ensure that the code is interruptible after a certain amount of time.
    This should only be used for internal doctesting purposes.

    EXAMPLES::

        sage: from sage.doctest.util import ensure_interruptible_after
        sage: with ensure_interruptible_after(1) as data: sleep(3)

    ``as data`` is optional, but if it is used, it will contain a few useful values::

        sage: data  # abs tol 0.2
        {'alarm_raised': True, 'elapsed': 1.0}

    ``max_wait_after_interrupt`` can be passed if the function may take longer than usual to be interrupted::

        sage: cython('''
        ....: from libc.time cimport clock_t, clock, CLOCKS_PER_SEC
        ....: from cysignals.signals cimport sig_check
        ....: cpdef void uninterruptible_sleep(double seconds):
        ....:     cdef clock_t target = clock() + <clock_t>(CLOCKS_PER_SEC * seconds)
        ....:     while clock() < target:
        ....:         pass
        ....: cpdef void check_interrupt_only_occasionally():
        ....:     for i in range(10):
        ....:         uninterruptible_sleep(0.8)
        ....:         sig_check()
        ....: ''')
        sage: with ensure_interruptible_after(1) as data:  # not passing max_wait_after_interrupt will raise an error
        ....:     check_interrupt_only_occasionally()
        Traceback (most recent call last):
        ...
        RuntimeError: Function is not interruptible within 1.0000 seconds, only after 1... seconds
        sage: with ensure_interruptible_after(1, max_wait_after_interrupt=0.9):
        ....:     check_interrupt_only_occasionally()

    TESTS::

        sage: data['elapsed']  # abs tol 0.3  # 1.6 = 0.8 * 2
        1.6

    ::

        sage: with ensure_interruptible_after(2) as data: sleep(1)
        Traceback (most recent call last):
        ...
        RuntimeError: Function terminates early after 1... < 2.0000 seconds
        sage: data  # abs tol 0.2
        {'alarm_raised': False, 'elapsed': 1.0}
        sage: with ensure_interruptible_after(1) as data: raise ValueError
        Traceback (most recent call last):
        ...
        ValueError
        sage: data  # abs tol 0.2
        {'alarm_raised': False, 'elapsed': 0.0}

    ::

        sage: # needs sage.misc.cython
        sage: with ensure_interruptible_after(1) as data: uninterruptible_sleep(2)
        Traceback (most recent call last):
        ...
        RuntimeError: Function is not interruptible within 1.0000 seconds, only after 2... seconds
        sage: data  # abs tol 0.2
        {'alarm_raised': True, 'elapsed': 2.0}
        sage: with ensure_interruptible_after(1): uninterruptible_sleep(2); raise RuntimeError
        Traceback (most recent call last):
        ...
        RuntimeError: Function is not interruptible within 1.0000 seconds, only after 2... seconds
        sage: data  # abs tol 0.2
        {'alarm_raised': True, 'elapsed': 2.0}
    """
    data = {}
    start_time = walltime()
    alarm(seconds)
    alarm_raised = False

    try:
        yield data
    except AlarmInterrupt:
        alarm_raised = True
    finally:
        cancel_alarm()
        elapsed = walltime() - start_time
        data["elapsed"] = elapsed
        data["alarm_raised"] = alarm_raised

    if elapsed > seconds + max_wait_after_interrupt:
        raise RuntimeError(
                f"Function is not interruptible within {seconds:.4f} seconds, only after {elapsed:.4f} seconds"
                + ("" if alarm_raised else " (__exit__ called before interrupt check)"))

    if alarm_raised:
        if elapsed < seconds - inaccuracy_tolerance:
            raise RuntimeError(f"Interrupted too early: {elapsed:.4f} < {seconds:.4f}, this should not happen")
    else:
        raise RuntimeError(f"Function terminates early after {elapsed:.4f} < {seconds:.4f} seconds")
