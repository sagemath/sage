# sage.doctest: needs sage.symbolic
"""
C Function Profiler Using Google Perftools

Note that the profiler samples 100x per second by default. In
particular, you cannot profile anything shorter than 10ms. You can
adjust the rate with the ``CPUPROFILE_FREQUENCY`` environment variable
if you want to change it.

EXAMPLES::

    sage: from sage.misc.gperftools import Profiler, run_100ms
    sage: prof = Profiler()
    sage: prof.start()       # optional - gperftools
    sage: run_100ms()
    sage: prof.stop()        # optional - gperftools
    PROFILE: interrupts/evictions/bytes = ...

REFERENCE:

Uses the `Google performance analysis tools
<https://github.com/gperftools/gperftools>`_. Note that they are not
included in Sage, you have to install them yourself on your system.

AUTHORS:

- Volker Braun (2014-03-31): initial version
"""

# ****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sys
import ctypes
import time
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.compat import find_library
from sage.cpython.string import bytes_to_str


libc = None
libprofiler = None


class Profiler(SageObject):

    def __init__(self, filename=None):
        """
        Interface to the gperftools profiler.

        INPUT:

        - ``filename`` -- string or ``None`` (default). The file name
          to log to. By default, a new temporary file is created.

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: Profiler()
            Profiler logging to ...
        """
        if filename is None:
            from sage.misc.temporary_file import tmp_filename
            self._filename = tmp_filename(ext='.perf')
        else:
            self._filename = filename

    def filename(self):
        """
        Return the file name.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: prof = Profiler()
            sage: prof.filename()
            '.../tmp_....perf'
        """
        return self._filename

    def _repr_(self):
        """
        Return string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: Profiler()
            Profiler logging to .../tmp....perf
        """
        return 'Profiler logging to {0}'.format(self.filename())

    def _libc(self):
        """
        Return libc.

        OUTPUT: a ctypes shared library handle

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: Profiler()._libc()
            <CDLL '...', handle ... at ...>
        """
        global libc
        if libc is not None:
            return libc
        name = find_library('c')
        if name:
            libc = ctypes.CDLL(name)
            return libc
        else:
            raise ImportError('failed to open libc')

    def _libprofiler(self):
        """
        Return libprofiler.

        OUTPUT: a ctypes shared library handle

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: Profiler()._libprofiler()    # optional - gperftools
            <CDLL 'libprofiler...', handle ... at ...>
        """
        global libprofiler
        if libprofiler is not None:
            return libprofiler
        import ctypes.util
        name = ctypes.util.find_library('profiler')
        if name:
            libprofiler = ctypes.CDLL(name)
            return libprofiler
        else:
            raise ImportError('failed to open libprofiler, make sure gperftools is installed')

    def start(self):
        """
        Start profiling.

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler, run_100ms
            sage: prof = Profiler()
            sage: prof.start()    # optional - gperftools
            sage: run_100ms()
            sage: prof.stop()     # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
        """
        from signal import SIGPROF, SIG_DFL
        self._previous_sigprof_handler = self._libc().signal(SIGPROF, SIG_DFL)
        profiler = self._libprofiler()
        self._t_start = time.time()
        rc = profiler.ProfilerStart(self.filename())
        if rc < 0:
            raise ValueError('profiler failed to start')

    def stop(self):
        """
        Stop the CPU profiler.

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler, run_100ms
            sage: prof = Profiler()
            sage: prof.start()    # optional - gperftools
            sage: run_100ms()
            sage: prof.stop()     # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
        """
        profiler = self._libprofiler()
        profiler.ProfilerStop()
        self._t_stop = time.time()
        if (self._t_stop - self._t_start) < 0.1:
            from warnings import warn
            warn('not enough samples, total runtime was '
                 'less than 100ms', RuntimeWarning)

    @cached_method
    def _pprof(self):
        """
        Return the name of the ``pprof`` binary.

        OUTPUT:

        String. The name of the gperftools ``pprof`` utility. A
        :exc:`OSError` is raised if it cannot be found.

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: prof = Profiler()
            sage: try:
            ....:     pp = prof._pprof()
            ....:     assert isinstance(pp, str)
            ....: except OSError:
            ....:     pass    # not installed
        """
        potential_names = ['google-pprof', 'pprof']
        from subprocess import check_output, CalledProcessError, STDOUT
        for name in potential_names:
            try:
                version = check_output([name, '--version'], stderr=STDOUT)
            except (CalledProcessError, OSError):
                continue
            version = bytes_to_str(version)
            if 'gperftools' not in version:
                from warnings import warn
                warn('the "{0}" utility does not appear to be the gperftools profiler'
                     .format(name), RuntimeWarning)
                continue
            return name
        raise OSError('unable to run pprof, please install gperftools')

    def _executable(self):
        """
        Return the name of the Sage Python interpreter.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: prof = Profiler()
            sage: prof._executable()
            '.../python...'
        """
        return sys.executable

    def _call_pprof(self, *args, **kwds):
        """
        Run the pprof binary.

        INPUT:

        - ``args`` -- list of strings; the arguments to ``pprof``

        - ``kwds`` -- keyword arguments passed to ``subprocess.check_call``

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: prof = Profiler()
            sage: prof._call_pprof('--help')   # optional - gperftools
            Usage:
            pprof [options] <program> <profiles>
            ...
        """
        from subprocess import check_call
        check_call([self._pprof()] + list(args), **kwds)

    def top(self, cumulative=True):
        """
        Print text report.

        OUTPUT: nothing; a textual report is printed to stdout

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler
            sage: prof = Profiler()
            sage: prof.start()    # optional - gperftools
            sage: # do something
            sage: prof.stop()     # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
            sage: prof.top()   # optional - gperftools
            Using local file ...
            Using local file ...
        """
        args = []
        if cumulative:
            args += ['--cum']
        args += ['--text', self._executable(), self.filename()]
        self._call_pprof(*args)

    def save(self, filename, cumulative=True, verbose=True):
        """
        Save report to disk.

        INPUT:

        - ``filename`` -- string; the filename to save at. Must end
          with one of ``.dot``, ``.ps``, ``.pdf``, ``.svg``, ``.gif``,
          or ``.txt`` to specify the output file format.

        - ``cumulative`` -- boolean (default: ``True``); whether to return
          cumulative timings

        - ``verbose`` -- boolean (default: ``True``); whether to print
          informational messages

        EXAMPLES::

            sage: from sage.misc.gperftools import Profiler, run_100ms
            sage: prof = Profiler()
            sage: prof.start()    # optional - gperftools
            sage: run_100ms()     # optional - gperftools
            sage: prof.stop()     # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
            sage: f = tmp_filename(ext='.txt')      # optional - gperftools
            sage: prof.save(f, verbose=False)       # optional - gperftools
        """
        args = []
        if cumulative:
            args += ['--cum']
        if filename.endswith('.dot'):
            args += ['--dot']
        elif filename.endswith('.ps'):
            args += ['--ps']
        elif filename.endswith('.pdf'):
            args += ['--pdf']
        elif filename.endswith('.svg'):
            args += ['--svg']
        elif filename.endswith('.txt'):
            args += ['--text']
        elif filename.endswith('.gif'):
            args += ['--gif']
        else:
            raise ValueError('unknown extension')
        args += [self._executable(), self.filename()]
        stderr = sys.stdout if verbose else False
        with open(filename, 'wb') as outfile:
            self._call_pprof(*args, stdout=outfile, stderr=stderr)


def crun(s, evaluator):
    """
    Profile single statement.

    - ``s`` -- string; Sage code to profile

    - ``evaluator`` -- callable to evaluate

    EXAMPLES::

        sage: import sage.misc.gperftools as gperf
        sage: ev = lambda ex:eval(ex, globals(), locals())
        sage: gperf.crun('gperf.run_100ms()', evaluator=ev)   # optional - gperftools
        PROFILE: interrupts/evictions/bytes = ...
        Using local file ...
        Using local file ...
    """
    prof = Profiler()
    from sage.repl.preparse import preparse
    py_s = preparse(s)
    prof.start()
    try:
        evaluator(py_s)
    finally:
        prof.stop()
    prof.top()


def run_100ms():
    """
    Used for doctesting.

    A function that performs some computation for more than (but not
    that much more than) 100ms.

    EXAMPLES::

        sage: from sage.misc.gperftools import run_100ms
        sage: run_100ms()
    """
    t0 = time.time()   # start
    t1 = t0 + 0.1      # end
    from sage.symbolic.ring import SR
    while time.time() < t1:
        sum(1 / (1 + SR(n) ** 2) for n in range(100))
