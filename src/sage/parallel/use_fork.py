"""
Parallel iterator built using the ``fork()`` system call
"""

# ****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from shutil import rmtree
from cysignals.alarm import AlarmInterrupt, alarm, cancel_alarm

from sage.interfaces.process import ContainChildren
from sage.misc.timing import walltime

from sage.misc.randstate import set_random_seed
from sage.misc.prandom import getrandbits

class WorkerData:
    """
    Simple class which stores data about a running ``p_iter_fork``
    worker.

    This just stores three attributes:

    - ``input`` -- the input value used by this worker

    - ``starttime`` -- the walltime when this worker started

    - ``failure`` -- an optional message indicating the kind of failure

    EXAMPLES::

        sage: from sage.parallel.use_fork import WorkerData
        sage: W = WorkerData(42); W
        <sage.parallel.use_fork.WorkerData object at ...>
        sage: W.starttime  # random
        1499330252.463206
    """
    def __init__(self, input_value, starttime=None, failure=""):
        r"""
        See the class documentation for description of the inputs.

        EXAMPLES::

            sage: from sage.parallel.use_fork import WorkerData
            sage: W = WorkerData(42)
        """
        self.input = input_value
        self.starttime = starttime or walltime()
        self.failure = failure


class p_iter_fork:
    """
    A parallel iterator implemented using ``fork()``.

    INPUT:

    - ``ncpus`` -- the maximal number of simultaneous subprocesses to spawn
    - ``timeout`` -- float (default: 0); wall time in seconds until a
      subprocess is automatically killed
    - ``verbose`` -- boolean (default: ``False``); whether to print anything
      about what the iterator does (e.g., killing subprocesses)
    - ``reset_interfaces`` -- boolean (default: ``True``); whether to reset all
      pexpect interfaces
    - ``reseed_rng`` -- boolean (default: ``False``); whether or not to reseed the rng in the subprocesses

    EXAMPLES::

        sage: X = sage.parallel.use_fork.p_iter_fork(2,3, False); X
        <sage.parallel.use_fork.p_iter_fork object at ...>
        sage: X.ncpus
        2
        sage: X.timeout
        3.0
        sage: X.verbose
        False
    """
    def __init__(self, ncpus, timeout=0, verbose=False, reset_interfaces=True, reseed_rng=False):
        """
        Create a ``fork()``-based parallel iterator.

        See the class documentation for description of the inputs.

        EXAMPLES::

            sage: X = sage.parallel.use_fork.p_iter_fork(2,3, False); X
            <sage.parallel.use_fork.p_iter_fork object at ...>
            sage: X.ncpus
            2
            sage: X.timeout
            3.0
            sage: X.verbose
            False
        """
        self.ncpus = int(ncpus)
        if self.ncpus != ncpus:  # check that there wasn't a roundoff
            raise TypeError("ncpus must be an integer")
        self.timeout = float(timeout)  # require a float
        self.verbose = verbose
        self.reset_interfaces = reset_interfaces
        self.reseed_rng = reseed_rng
        self.worker_seed = None

    def __call__(self, f, inputs):
        """
        Parallel iterator using ``fork()``.

        INPUT:

        - ``f`` -- a function (or more general, any callable)

        - ``inputs`` -- list of pairs ``(args, kwds)`` to be used as
          arguments to ``f``, where ``args`` is a tuple and ``kwds`` is
          a dictionary

        EXAMPLES::

            sage: F = sage.parallel.use_fork.p_iter_fork(2,3)
            sage: sorted(list( F( (lambda x: x^2), [([10],{}), ([20],{})])))
            [(([10], {}), 100), (([20], {}), 400)]
            sage: sorted(list( F( (lambda x, y: x^2+y), [([10],{'y':1}), ([20],{'y':2})])))
            [(([10], {'y': 1}), 101), (([20], {'y': 2}), 402)]

        TESTS:

        The output of functions decorated with :func:`parallel` is read
        as a pickle by the parent process. We intentionally break the
        unpickling and demonstrate that this failure is handled
        gracefully (the exception is put in the list instead of the
        answer)::

            sage: Polygen = parallel(polygen)
            sage: list(Polygen([QQ]))
            [(((Rational Field,), {}), x)]
            sage: from sage.misc.persist import unpickle_override, register_unpickle_override
            sage: register_unpickle_override('sage.rings.polynomial.polynomial_rational_flint', 'Polynomial_rational_flint', Integer)
            sage: L = list(Polygen([QQ]))
            sage: L
            [(((Rational Field,), {}),
              'INVALID DATA __init__() takes at most 2 positional arguments (4 given)')]

        Fix the unpickling::

            sage: del unpickle_override[('sage.rings.polynomial.polynomial_rational_flint', 'Polynomial_rational_flint')]
            sage: list(Polygen([QQ,QQ]))
            [(((Rational Field,), {}), x), (((Rational Field,), {}), x)]
        """
        import os
        import sys
        import signal
        from sage.misc.persist import loads
        from sage.misc.temporary_file import tmp_dir
        dir = tmp_dir()
        timeout = self.timeout

        n = self.ncpus
        inputs = list(inputs)
        if self.reseed_rng:
            seeds = [getrandbits(512) for _ in range(len(inputs))]
            vs = list(zip(inputs, seeds))
        else:
            vs = list(zip(inputs, [None]*len(inputs)))
        workers = {}
        try:
            while vs or workers:
                # Spawn up to n subprocesses
                while vs and len(workers) < n:
                    v0, seed0 = vs.pop(0)  # Input value and seed for the next subprocess
                    with ContainChildren():
                        pid = os.fork()
                        # The way fork works is that pid returns the
                        # nonzero pid of the subprocess for the master
                        # process and returns 0 for the subprocess.
                        if not pid:
                            # This is the subprocess.
                            if self.reseed_rng:
                                self.worker_seed = seed0
                            self._subprocess(f, dir, *v0)
                    workers[pid] = WorkerData(v0)

                if len(workers) > 0:
                    # Now wait for one subprocess to finish and report the result.
                    # However, wait at most the time since the oldest process started.
                    T = walltime()
                    if timeout:
                        oldest = min(W.starttime for W in workers.values())
                        alarm(max(timeout - (T - oldest), 0.1))

                    try:
                        pid = os.wait()[0]
                        cancel_alarm()
                        W = workers.pop(pid)
                    except AlarmInterrupt:
                        # Kill workers that are too old
                        for pid, W in workers.items():
                            if T - W.starttime > timeout:
                                if self.verbose:
                                    print(
                                        "Killing subprocess %s with input %s which took too long" % (pid, W.input))
                                os.kill(pid, signal.SIGKILL)
                                W.failure = " (timed out)"
                    except KeyError:
                        # Some other process exited, not our problem...
                        pass
                    else:
                        # collect data from process that successfully terminated
                        sobj = os.path.join(dir, '%s.sobj' % pid)
                        try:
                            with open(sobj, "rb") as file:
                                data = file.read()
                        except OSError:
                            answer = "NO DATA" + W.failure
                        else:
                            os.unlink(sobj)
                            try:
                                answer = loads(data, compress=False)
                            except Exception as E:
                                answer = "INVALID DATA {}".format(E)

                        out = os.path.join(dir, '%s.out' % pid)
                        try:
                            with open(out) as file:
                                sys.stdout.write(file.read())
                            os.unlink(out)
                        except OSError:
                            pass

                        yield (W.input, answer)
        finally:
            # Send SIGKILL signal to workers that are left.
            if workers:
                if self.verbose:
                    print("Killing any remaining workers...")
                sys.stdout.flush()
                for pid in workers:
                    try:
                        os.kill(pid, signal.SIGKILL)
                    except OSError:
                        # If kill() failed, it is most likely because
                        # the process already exited.
                        pass
                    else:
                        try:
                            os.waitpid(pid, 0)
                        except OSError as msg:
                            if self.verbose:
                                print(msg)

            # Clean up all temporary files.
            rmtree(dir)

    def _subprocess(self, f, dir, args, kwds={}):
        """
        Setup and run evaluation of ``f(*args, **kwds)``, storing the
        result in the given directory ``dir``.

        This method is called by each forked subprocess.

        INPUT:

        - ``f`` -- a function

        - ``dir`` -- name of a directory

        - ``args`` -- tuple with positional arguments for ``f``

        - ``kwds`` -- (optional) a dict with keyword arguments for ``f``

        TESTS:

        The method ``_subprocess`` is really meant to be run only in a
        subprocess. It doesn't print not return anything, the output is
        saved in pickles. It redirects stdout, so we save and later
        restore stdout in order not to break the doctester::

            sage: saved_stdout = sys.stdout
            sage: F = sage.parallel.use_fork.p_iter_fork(2,3)
            sage: F._subprocess(operator.add, tmp_dir(), (1, 2))
            sage: sys.stdout = saved_stdout
        """
        import os
        import sys
        try:
            from importlib import reload
        except ImportError:
            from imp import reload
        from sage.misc.persist import save

        # Make it so all stdout is sent to a file so it can
        # be displayed.
        out = os.path.join(dir, '%s.out' % os.getpid())
        sys.stdout = open(out, 'w')

        # Run some commands to tell Sage that its
        # pid has changed (forcing a reload of
        # misc).
        import sage.misc.misc
        reload(sage.misc.misc)

        # The pexpect interfaces (and objects defined in them) are
        # not valid.
        if self.reset_interfaces:
            try:
                from sage.interfaces.quit import invalidate_all
            except ImportError:
                pass
            else:
                invalidate_all()

        # Reseed rng, if requested.
        if self.reseed_rng:
            set_random_seed(self.worker_seed)

        # Now evaluate the function f.
        value = f(*args, **kwds)

        # And save the result to disk.
        sobj = os.path.join(dir, '%s.sobj' % os.getpid())
        save(value, sobj, compress=False)
