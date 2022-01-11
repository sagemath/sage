"""Miscellaneous utilities for running the docbuilder."""

import errno
import os


class WorkerDiedException(RuntimeError):
    """Raised if a worker process dies unexpected."""

    def __init__(self, message, original_exception=None):
        super(WorkerDiedException, self).__init__(message)
        self.original_exception = original_exception


def build_many(target, args, processes=None):
    """
    Map a list of arguments in ``args`` to a single-argument target function
    ``target`` in parallel using ``multiprocessing.cpu_count()`` (or
    ``processes`` if given) simultaneous processes.

    This is a simplified version of ``multiprocessing.Pool.map`` from the
    Python standard library which avoids a couple of its pitfalls.  In
    particular, it can abort (with a `RuntimeError`) without hanging if one of
    the worker processes unexpectedly dies.  It also has semantics equivalent
    to ``maxtasksperchild=1``; that is, one process is started per argument.
    As such, this is inefficient for processing large numbers of fast tasks,
    but appropriate for running longer tasks (such as doc builds) which may
    also require significant cleanup.

    It also avoids starting new processes from a pthread, which results in at
    least two known issues:

        * On versions of Cygwin prior to 3.0.0 there were bugs in mmap handling
          on threads (see https://trac.sagemath.org/ticket/27214#comment:25).

        * When PARI is built with multi-threading support, forking a Sage
          process from a thread leaves the main Pari interface instance broken
          (see https://trac.sagemath.org/ticket/26608#comment:38).

    In the future this may be replaced by a generalized version of the more
    robust parallel processing implementation from ``sage.doctest.forker``.

    EXAMPLES::

        sage: from sage_docbuild.utils import build_many
        sage: def target(N):
        ....:     import time
        ....:     time.sleep(float(0.1))
        ....:     print('Processed task %s' % N)
        ....:
        sage: _ = build_many(target, range(8), processes=8)
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...
        Processed task ...

    Unlike the first version of `build_many` which was only intended to get
    around the Cygwin bug, this version can also return a result, and thus can
    be used as a replacement for `multiprocessing.Pool.map` (i.e. it still
    blocks until the result is ready)::

        sage: def square(N):
        ....:     return N * N
        sage: build_many(square, range(100))
        [0, 1, 4, 9, ..., 9604, 9801]

    If the target function raises an exception in any of the workers,
    `build_many` raises that exception and all other results are discarded.
    Any in-progress tasks may still be allowed to complete gracefully before
    the exception is raised::

        sage: def target(N):
        ....:     import time, os, signal
        ....:     if N == 4:
        ....:         # Task 4 is a poison pill
        ....:         1 / 0
        ....:     else:
        ....:         time.sleep(float(0.5))
        ....:         print('Processed task %s' % N)
        ....:

    Note: In practice this test might still show output from the other worker
    processes before the poison-pill is executed.  It may also display the
    traceback from the failing process on stderr.  However, due to how the
    doctest runner works, the doctest will only expect the final exception::

        sage: build_many(target, range(8), processes=8)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: rational division by zero

    Similarly, if one of the worker processes dies unexpectedly otherwise exits
    non-zero (e.g. killed by a signal) any in-progress tasks will be completed
    gracefully, but then a `RuntimeError` is raised and pending tasks are not
    started::

        sage: def target(N):
        ....:     import time, os, signal
        ....:     if N == 4:
        ....:         # Task 4 is a poison pill
        ....:         os.kill(os.getpid(), signal.SIGKILL)
        ....:     else:
        ....:         time.sleep(float(0.5))
        ....:         print('Processed task %s' % N)
        ....:
        sage: build_many(target, range(8), processes=8)
        Traceback (most recent call last):
        ...
        WorkerDiedException: worker for 4 died with non-zero exit code -9
    """
    from multiprocessing import Process, Queue, cpu_count, set_start_method
    # With OS X, Python 3.8 defaults to use 'spawn' instead of 'fork'
    # in multiprocessing, and Sage docbuilding doesn't work with
    # 'spawn'. See trac #27754.
    if os.uname().sysname == 'Darwin':
        set_start_method('fork', force=True)
    from queue import Empty

    if processes is None:
        processes = cpu_count()

    workers = [None] * processes
    tasks = enumerate(args)
    results = []
    result_queue = Queue()

    # Utility functions #
    def run_worker(target, queue, idx, task):
        try:
            result = target(task)
        except BaseException as exc:
            queue.put((None, exc))
        else:
            queue.put((idx, result))

    def bring_out_yer_dead(w, task, exitcode):
        """
        Handle a dead / completed worker.  Raises WorkerDiedException if it
        returned with a non-zero exit code.
        """

        if w is None or exitcode is None:
            # I'm not dead yet! (or I haven't even been born yet)
            return (w, task)

        # Hack: If we wait()ed on this worker manually we have to tell it
        # it's dead:
        if w._popen.returncode is None:
            w._popen.returncode = exitcode

        if exitcode != 0:
            raise WorkerDiedException(
                "worker for {} died with non-zero exit code "
                "{}".format(task[1], w.exitcode))

        # Get result from the queue; depending on ordering this may not be
        # *the* result for this worker, but for each completed worker there
        # should be *a* result so let's get it
        try:
            result = result_queue.get_nowait()
        except Empty:
            # Generally shouldn't happen but could in case of a race condition;
            # don't worry we'll collect any remaining results at the end.
            pass

        if result[0] is None:
            # Indicates that an exception occurred in the target function
            raise WorkerDiedException('', original_exception=result[1])
        else:
            results.append(result)

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
                w, task = w
                if w.pid == waited_pid:
                    exitcode = waited_exitcode
                else:
                    exitcode = w.exitcode

                w = bring_out_yer_dead(w, task, exitcode)

            # Worker w is dead/not started, so start a new worker
            # in its place with the next document from the queue
            if w is None:
                try:
                    task = next(tasks)
                except StopIteration:
                    pass
                else:
                    w = Process(target=run_worker,
                                args=((target, result_queue) + task))
                    w.start()
                    # Pair the new worker with the task it's performing (mostly
                    # for debugging purposes)
                    w = (w, task)

            workers[idx] = w

            if w is not None:
                all_done = False

        # If all workers are dead and there are no more items to
        # process in the queue then we are done
        return all_done

    # Main loop #

    waited_pid = None
    # Set along with waited_exitcode by calls to wait_for_one()

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
            for w, _ in remaining_workers:
                # Give any remaining workers a chance to shut down gracefully
                try:
                    w.terminate()
                except OSError as exc:
                    if exc.errno != errno.ESRCH:
                        # Otherwise it was already dead so this was expected
                        raise
            for w, _ in remaining_workers:
                w.join()
        finally:
            if worker_exc is not None:
                # Re-raise the RuntimeError from bring_out_yer_dead set if a
                # worker died unexpectedly, or the original exception if it's
                # wrapping one
                if worker_exc.original_exception:
                    raise worker_exc.original_exception
                else:
                    raise worker_exc

    # All workers should be shut down by now and should have completed without
    # error.  No new items will be added to the result queue, so we can get all
    # the remaining results, if any.
    while True:
        try:
            results.append(result_queue.get_nowait())
        except Empty:
            break

    # Return the results sorted according to their original task order
    return [r[1] for r in sorted(results, key=lambda r: r[0])]
