# sage_setup: distribution = sagemath-repl
r"""
Simple profiling tool

AUTHORS:

- David Harvey (August 2006)
- Martin Albrecht
"""

# *****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Harvey <dmharvey@math.harvard.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.timing import cputime
import inspect
import sys


class Profiler:
    r"""
    Keeps track of CPU time used between a series of user-defined checkpoints.

    It's probably not a good idea to use this class in an inner loop :-)

    EXAMPLES::

        from sage.misc.profiler import Profiler
        sage: def f():                        # not tested
        ....:     p = Profiler()

    Calling ``p(message)`` creates a checkpoint::

        sage: p("try factoring 15")           # not tested

    Do something time-consuming::

        sage: x = factor(15)                  # not tested

    You can create a checkpoints without a string; ``Profiler``
    will use the source code instead::

        sage: # not tested
        sage: p()
        sage: y = factor(25)
        sage: p("last step")
        sage: z = factor(35)
        sage: p()

    This will give a nice list of timings between checkpoints::

        sage: print(p)                        # not tested

    Let's try it out::

        sage: f()                             # not tested
            3.020s -- try factoring 15
           15.240s -- line 17: y = factor(25)
         5000.190s -- last step

    .. SEEALSO:: :func:`runsnake`

    .. TODO::

        - Add Pyrex source code inspection (I assume it doesn't
          currently do this)
        - Add ability to sort output by time
        - Add option to constructor to print timing immediately when
          checkpoint is reached
        - Migrate to Pyrex?
        - Add ability to return timings in a more machine-friendly
          format

    AUTHOR:

    - David Harvey (August 2006)
    """

    def __init__(self, systems=[], verbose=False):
        """
        INPUT:

        - ``systems`` -- list of interfaces to other system which implements a
          cputime method. The cputimes of all provided systems will be added
          to the cputime of Sage itself.
        """
        systems = [e.cputime for e in systems]
        self._cputime_functions = [cputime] + list(systems)
        self._verbose = bool(verbose)
        self.clear()

    def clear(self):
        # _checkpoints is list of pairs (details, time), where time is a float
        # and details is a triple (line_number, context, message)
        self._checkpoints = []
        self._active_details = None   # details from the last __call__() call
        self._last_cputime = [None]*len(self._cputime_functions)

    def __call__(self, message=None):
        """ Adds a checkpoint. """
        entry_times = [fn() for fn in self._cputime_functions]

        frame = inspect.currentframe().f_back
        try:
            frame_info = inspect.getframeinfo(frame, 9)
            line_number = frame_info[1]
            context = frame_info[3]
        finally:
            # This is here to prevent reference cycles
            # (according to the warning in the python documentation:
            # http://docs.python.org/lib/inspect-stack.html):
            del frame_info
            del frame

        if self._active_details is not None:
            _time = sum([entry_times[i]-self._last_cputime[i] for i in range(len(entry_times))])
            self._checkpoints.append((self._active_details, _time))

        self._active_details = (line_number, context, message)

        self._last_cputime = [fn() for fn in self._cputime_functions]
        if self._verbose:
            print(self.print_last())
            sys.stdout.flush()

    def __repr__(self):
        """ Returns a nicely formatted table of stored checkpoints and timings. """
        if not self._checkpoints:
            return "no checkpoints defined"

        output = []
        for ((line_number, context, message), time_used) in self._checkpoints:
            if message is None:
                # If the user hasn't given a message, we look for some
                # source code to print instead
                found = "(unknown)"
                for line in context[5:]:
                    line = line.strip()
                    if line != "":
                        found = line
                        break

                if len(found) > 60:
                    found = found[:60] + "..."   # in case the source line is really long
                message = "line %d: %s" % (line_number, found)

            output.append("%9.3fs -- %s" % (time_used, message))

        return "\n".join(output)

    def print_last(self):
        """
        Print the last profiler step.
        """
        if not self._checkpoints:
            return ""

        ((line_number, context, message), time_used) = self._checkpoints[-1]
        if message is None:
            # If the user hasn't given a message, we look for some
            # source code to print instead
            found = "(unknown)"
            for line in context[5:]:
                line = line.strip()
                if line != "":
                    found = line
                    break

            if len(found) > 60:
                found = found[:60] + "..."   # in case the source line is really long
            message = "line %d: %s" % (line_number, found)

        return "%9.3fs -- %s" % (time_used, message)

# end of file
