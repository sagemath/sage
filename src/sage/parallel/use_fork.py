"""
Parallel iterator built using the ``fork()`` system call
"""

#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function

from cysignals.alarm import AlarmInterrupt, alarm, cancel_alarm

class p_iter_fork:
    """
    A parallel iterator implemented using ``fork()``.

    EXAMPLES::

        sage: X = sage.parallel.use_fork.p_iter_fork(2,3, False); X
        <sage.parallel.use_fork.p_iter_fork instance at ...>
        sage: X.ncpus
        2
        sage: X.timeout
        3.0
        sage: X.verbose
        False
    """
    def __init__(self, ncpus, timeout=0, verbose=False, reset_interfaces=True):
        """
        Create a ``fork()``-based parallel iterator.

        INPUT:

            - ``ncpus`` -- the maximal number of simultaneous
              subprocesses to spawn
            - ``timeout`` -- (float, default: 0) wall time in seconds until
              a subprocess is automatically killed
            - ``verbose`` -- (default: False) whether to print
              anything about what the iterator does (e.g., killing
              subprocesses)
            - ``reset_interfaces`` -- (default: True) whether to reset
              all pexpect interfaces

        EXAMPLES::

            sage: X = sage.parallel.use_fork.p_iter_fork(2,3, False); X
            <sage.parallel.use_fork.p_iter_fork instance at ...>
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

    def __call__(self, f, inputs):
        """
        Parallel iterator using ``fork()``.

        INPUT:

        - ``f`` -- a Python function that need not be pickleable or anything else!

        - ``inputs`` -- a list of pickleable pairs ``(args, kwds)``, where
          ``args`` is a tuple and ``kwds`` is a dictionary.

        OUTPUT:

        EXAMPLES::

            sage: F = sage.parallel.use_fork.p_iter_fork(2,3)
            sage: sorted(list( F( (lambda x: x^2), [([10],{}), ([20],{})])))
            [(([10], {}), 100), (([20], {}), 400)]
            sage: sorted(list( F( (lambda x, y: x^2+y), [([10],{'y':1}), ([20],{'y':2})])))
            [(([10], {'y': 1}), 101), (([20], {'y': 2}), 402)]

        TESTS:

        The output of functions decorated with :func:parallel is read
        as a pickle by the parent process. We intentionally break the
        unpickling and demonstrate that this failure is handled
        gracefully (an exception is displayed and an empty list is
        returned)::

            sage: Polygen = parallel(polygen)
            sage: list(Polygen([QQ]))
            [(((Rational Field,), {}), x)]
            sage: from sage.structure.sage_object import unpickle_override, register_unpickle_override
            sage: register_unpickle_override('sage.rings.polynomial.polynomial_rational_flint', 'Polynomial_rational_flint', Integer)
            sage: L = list(Polygen([QQ]))
            ('__init__() takes at most 2 positional arguments (4 given)', <type 'sage.rings.integer.Integer'>, (Univariate Polynomial Ring in x over Rational Field, [0, 1], False, True))
            sage: L
            []

        Fix the unpickling::

            sage: del unpickle_override[('sage.rings.polynomial.polynomial_rational_flint', 'Polynomial_rational_flint')]
            sage: list(Polygen([QQ,QQ]))
            [(((Rational Field,), {}), x), (((Rational Field,), {}), x)]
        """
        n = self.ncpus
        v = list(inputs)
        import os, sys, signal
        from sage.structure.sage_object import load
        from sage.misc.all import tmp_dir, walltime
        dir = tmp_dir()
        timeout = self.timeout

        workers = {}
        try:
            while len(v) > 0 or len(workers) > 0:
                # Spawn up to n subprocesses
                while len(v) > 0 and len(workers) < n:
                    # Subprocesses shouldn't inherit unflushed buffers (cf. #11778):
                    sys.stdout.flush()
                    sys.stderr.flush()

                    pid = os.fork()
                    # The way fork works is that pid returns the
                    # nonzero pid of the subprocess for the master
                    # process and returns 0 for the subprocess.
                    if pid:
                        # This is the parent master process.
                        workers[pid] = [v[0], walltime(), '']
                        del v[0]
                    else:
                        # This is the subprocess.
                        self._subprocess(f, dir, v[0])

                if len(workers) > 0:
                    # Now wait for one subprocess to finish and report the result.
                    # However, wait at most the time since the oldest process started.
                    if timeout:
                        oldest = min([X[1] for X in workers.values()])
                        alarm(max(timeout - (walltime()-oldest), 0.1))

                    try:
                        pid = os.wait()[0]
                        cancel_alarm()
                        w = workers.pop(pid)
                    except AlarmInterrupt:
                        cancel_alarm()
                        # Kill workers that are too old
                        for pid, X in workers.iteritems():
                            if walltime() - X[1] > timeout:
                                if self.verbose:
                                    print(
                                        "Killing subprocess %s with input %s which took too long"
                                         % (pid, X[0]) )
                                os.kill(pid, signal.SIGKILL)
                                X[-1] = ' (timed out)'
                    except KeyError:
                        # Some other process exited, not our problem...
                        pass
                    else:
                        # collect data from process that successfully terminated
                        sobj = os.path.join(dir, '%s.sobj'%pid)
                        if not os.path.exists(sobj):
                            X = "NO DATA" + w[-1]  # the message field
                        else:
                            X = load(sobj, compress=False)
                            os.unlink(sobj)
                        out = os.path.join(dir, '%s.out'%pid)
                        if not os.path.exists(out):
                            output = "NO OUTPUT"
                        else:
                            output = open(out).read()
                            os.unlink(out)

                        if output.strip():
                            print(output, end="")

                        yield (w[0], X)

        except Exception as msg:
            print(msg)

        finally:
            # Clean up all temporary files.
            try:
                for X in os.listdir(dir):
                    os.unlink(os.path.join(dir, X))
                os.rmdir(dir)
            except OSError as msg:
                if self.verbose:
                    print(msg)

            # Send "kill -9" signal to workers that are left.
            if len(workers) > 0:
                if self.verbose:
                    print("Killing any remaining workers...")
                sys.stdout.flush()
                for pid in workers.keys():
                    try:
                        os.kill(pid, signal.SIGKILL)
                        os.waitpid(pid, 0)
                    except OSError as msg:
                        if self.verbose:
                            print(msg)

    def _subprocess(self, f, dir, x):
        """
        Setup and run evaluation of ``f(*x[0], **x[1])``, storing the
        result in the given directory ``dir``.  This method is called by each
        forked subprocess.

        INPUT:

            - ``f`` -- a function
            - ``dir`` -- name of a directory
            - ``x`` -- 2-tuple, with args and kwds

        EXAMPLES:

        We have only this indirect test, since a direct test would terminate the Sage session.

            sage: F = sage.parallel.use_fork.p_iter_fork(2,3)
            sage: sorted(list( F( (lambda x: x^2), [([10],{}), ([20],{})])))
            [(([10], {}), 100), (([20], {}), 400)]
        """
        import imp, os, sys
        from sage.structure.sage_object import save

        try:
            # Make it so all stdout is sent to a file so it can
            # be displayed.
            out = os.path.join(dir, '%s.out'%os.getpid())
            sys.stdout = open(out, 'w')

            # Run some commands to tell Sage that its
            # pid has changed (forcing a reload of
            # misc).
            import sage.misc.misc
            imp.reload(sage.misc.misc)

            # The pexpect interfaces (and objects defined in them) are
            # not valid.
            if self.reset_interfaces:
                sage.interfaces.quit.invalidate_all()

            # Now evaluate the function f.
            value = f(*x[0], **x[1])

            # And save the result to disk.
            sobj = os.path.join(dir, '%s.sobj'%os.getpid())
            save(value, sobj, compress=False)

        except Exception as msg:
            # Important to print this, so it is seen by the caller.
            print(msg)
        finally:
            sys.stdout.flush()
            os._exit(0)
