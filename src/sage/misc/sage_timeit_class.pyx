# sage_setup: distribution = sagemath-objects
# This is here in Cython so we can access the interpreter globals
# cython: old_style_globals=True
"""
The ``timeit`` command

This uses the function :func:`~sage.misc.sage_timeit.sage_timeit`.
"""

from sage.misc import sage_timeit


class SageTimeit:
    r"""nodetex
    Time execution of a command or block of commands.

    Displays the best WALL TIME for execution of the given code. This
    is based on the Python ``timeit`` module, which avoids a number of
    common traps for measuring execution times.  It is also based on
    IPython's ``%timeit`` command.

    TYPICAL INPUT FORMAT::

        timeit(statement, preparse=None, number=0, repeat=3, precision=3)

    EXAMPLES::

        sage: timeit('2^10000')
        625 loops, best of 3: ... per loop

    We illustrate some options::

        sage: timeit('2+2',precision=2,number=20,repeat=5)
        20 loops, best of 5: ... per loop

    The preparser is on by default (if it is on), but the preparse option
    allows us to override it::

        sage: timeit('2^10000', preparse=False, number=50)
        50 loops, best of 3: ... per loop

    The input can contain newlines::

        sage: timeit("a = 2\nb=131\nfactor(a^b-1)", number=25)                          # needs sage.libs.pari
        25 loops, best of 3: ... per loop

    .. SEEALSO:: :func:`runsnake`
    """
    def eval(self, code, globs=None, locals=None, **kwds):
        r"""
        This eval function is called when doing \%timeit in the notebook.

        INPUT:

        - ``code`` -- string of code to evaluate; may contain newlines

        - ``globs`` -- global variables; if not given, uses module scope
          globals

        - ``locals`` -- ignored completely

        - ``kwds`` -- passed onto sage_timeit. Common options are
          ``preparse``, ``number``, ``repeat``, ``precision``. See
          :func:`~sage.misc.sage_timeit.sage_timeit` for details.

        OUTPUT: string; timing information

        EXAMPLES::

            sage: timeit.eval("2+2")                     # random output
            '625 loops, best of 3: 1.47 us per loop'

        We emphasize that timeit times WALL TIME. This is good in the
        context of Sage where commands often call out to other
        subprocesses that do not appear in CPU time. ::

            sage: timeit('sleep(float(0.5))', number=3)  # long time (5s on sage.math, 2012)
            3 loops, best of 3: ... ms per loop
        """
        if globs is None:
            globs = globals()
        return sage_timeit.sage_timeit(code, globs, **kwds)

    def __call__(self, code, globals=None, **kwds):
        """
        INPUT:

        - ``code`` -- string; a line or block of code, which may
          contain newlines

        - ``globals`` -- (optional) global variables; if not given the
          globals of the calling module are used (e.g., if using this
          from the command line, the globals of the command line are
          used).

        - ``preparse`` -- boolean or ``None`` (default). Whether or
          not to preparse the input code using the Sage preparser.  If
          not specified, do the same thing as whatever was set by the
          preparser command.

        - ``**kwds`` -- passed onto self.eval(...). Common options are
          ``preparse``, ``number``, ``repeat``, ``precision``. See
          :func:`~sage.misc.sage_timeit.sage_timeit` for details.

        OUTPUT:

        Return the timing information, either as
        :class:`~sage.misc.sage_timeit.SageTimeitResult`
        or as a floating-point number containing the number of seconds
        (if ``seconds=True`` was passed).

        EXAMPLES::

            sage: timeit('2^10000', preparse=False, number=100)
            100 loops, best of 3: ... per loop
        """
        return self.eval(code, globals, **kwds)


timeit = SageTimeit()
