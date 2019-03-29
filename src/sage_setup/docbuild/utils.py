"""Miscellaneous utilities for running the docbuilder."""

import errno
import os


class WorkerDiedException(RuntimeError):
    """Raised if a worker process dies unexpected."""


def _build_many(target, args, processes=None):
    """
    Map a list of arguments in ``args`` to a single-argument target function
    ``target`` in parallel using ``NUM_THREADS`` (or ``processes`` if given)
    simultaneous processes.

    This is a simplified version of ``multiprocessing.Pool.map`` from the
    Python standard library which avoids a couple of its pitfalls.  In
    particular, it can abort (with a `RuntimeError`) without hanging if one of
    the worker processes unexpectedly dies.  It also avoids starting new
    processes from a pthread, which is known to result in bugs on versions of
    Cygwin prior to 3.0.0 (see
    https://trac.sagemath.org/ticket/27214#comment:25).

    On the other hand, unlike ``multiprocessing.Pool.map`` it does not return
    a result.  This is fine for the purpose of building multiple Sphinx
    documents in parallel.

    In the future this may be replaced by a generalized version of the more
    robust parallel processing implementation from ``sage.doctest.forker``.

    EXAMPLES::

        sage: from sage_setup.docbuild.utils import _build_many
        sage: def target(N):
        ....:     import time
        ....:     time.sleep(float(0.1))
        ....:     print('Processed task %s' % N)
        ....:
        sage: _build_many(target, range(8), processes=8)
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...

    If one of the worker processes errors out from an unhandled exception, or
    otherwise exits non-zero (e.g. killed by a signal) any in-progress tasks
    will be completed gracefully, but then a `RuntimeError` is raised and
    pending tasks are not started::

        sage: def target(N):
        ....:     import time
        ....:     if N == 4:
        ....:         # Task 4 is a poison pill
        ....:         1 / 0
        ....:     else:
        ....:         time.sleep(0.5)
        ....:         print('Processed task %s' % N)
        ....:

    Note: In practice this test might still show output from the other worker
    processes before the poison-pill is executed.  It may also display the
    traceback from the failing process on stderr.  However, due to how the
    doctest runner works, the doctest will only expect the final exception::

        sage: _build_many(target, range(8), processes=8)
        Traceback (most recent call last):
        ...
        WorkerDiedException: worker for 4 died with non-zero exit code 1
    """
    from multiprocessing import Process
    from .build_options import NUM_THREADS, ABORT_ON_ERROR

    if processes is None:
        processes = NUM_THREADS

    workers = [None] * processes
    queue = list(args)

    # Maps worker process PIDs to the name of the document it's working
    # on (the argument it was passed).  This is primarily used just for
    # debugging/information purposes.
    jobs = {}

    ### Utility functions ###

    def bring_out_yer_dead(w, exitcode):
        """
        Handle a dead / completed worker.  Raises WorkerDiedError if it
        returned with a non-zero exit code.
        """

        if w is None or exitcode is None:
            # I'm not dead yet! (or I haven't even been born yet)
            return w

        # Hack: If we wait()ed on this worker manually we have to tell it
        # it's dead:
        if w._popen.returncode is None:
            w._popen.returncode = exitcode

        if exitcode != 0 and ABORT_ON_ERROR:
            raise WorkerDiedException(
                "worker for {} died with non-zero exit code "
                "{}".format(jobs[w.pid], w.exitcode))

        jobs.pop(w.pid)
        # Helps multiprocessing with some internal bookkeeping
        w.join()

        return None

    def wait_for_one():
        """Wait for a single process and return its pid and exit code."""
        try:
            pid, sts = os.wait()
        except OSError as exc:
            # No more processes to wait on if ECHILD
            if exc.errno != errno.ECHILD:
                raise
            else:
                return None, None

        if os.WIFSIGNALED(sts):
            exitcode = -os.WTERMSIG(sts)
        else:
            exitcode = os.WEXITSTATUS(sts)

        return pid, exitcode

    def reap_workers(waited_pid=None, waited_exitcode=None):
        """
        This is the main worker handling loop.

        Checks if workers have completed their tasks and spawns new workers if
        there are more tasks on the queue.  Returns `False` if there is more
        work to be done or `True` if the work is complete.

        Raises a ``WorkerDiedException`` if a worker exits unexpectedly.
        """

        all_done = True

        for idx, w in enumerate(workers):
            if w is not None:
                if w.pid == waited_pid:
                    exitcode = waited_exitcode
                else:
                    exitcode = w.exitcode

                w = bring_out_yer_dead(w, exitcode)

            # Worker w is dead/not started, so start a new worker
            # in its place with the next document from the queue
            if w is None and queue:
                job = queue.pop(0)
                w = Process(target=target, args=(job,))
                w.start()
                jobs[w.pid] = job

            workers[idx] = w

            if w is not None:
                all_done = False

        # If all workers are dead and there are no more items to
        # process in the queue then we are done
        return all_done

    ### Main loop ###

    waited_pid = None  # Set along with waited_exitcode by calls to
                       # wait_for_one()
    waited_exitcode = None
    worker_exc = None  # Set to a WorkerDiedException if one occurs

    try:
        while True:
            # Check the status of each worker and break out of the loop if
            # all work is done.
            # We'll check each worker process against the returned
            # pid back at the top of the `while True` loop.  We also
            # check any other processes that may have exited in the
            # meantime
            try:
                if reap_workers(waited_pid, waited_exitcode):
                    break
            except WorkerDiedException as exc:
                worker_exc = exc
                break

            waited_pid, waited_exitcode = wait_for_one()
    finally:
        try:
            remaining_workers = [w for w in workers if w is not None]
            for w in remaining_workers:
                # Give any remaining workers a chance to shut down gracefully
                try:
                    w.terminate()
                except OSError as exc:
                    if exc.errno != errno.ESRCH:
                        # Otherwise it was already dead so this was expected
                        raise
            for w in remaining_workers:
                w.join()
        finally:
            if worker_exc is not None:
                # Re-raise the RuntimeError from bring_out_yer_dead set if a
                # worker died unexpectedly
                raise worker_exc
