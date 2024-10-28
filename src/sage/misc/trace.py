"""
Interactively tracing execution of a command
"""


def trace(code, preparse=True):
    r"""
    Evaluate Sage code using the interactive tracer and return the
    result. The string ``code`` must be a valid expression
    enclosed in quotes (no assignments - the result of the expression
    is returned). In the Sage notebook this just raises a
    NotImplementedException.

    INPUT:

    - ``code`` -- string

    - ``preparse`` -- boolean (default: ``True``); if ``True``, run
      expression through the Sage preparser

    REMARKS: This function is extremely powerful! For example, if you
    want to step through each line of execution of, e.g.,
    ``factor(100)``, type

    ::

        sage: from sage.misc.trace import trace
        sage: trace("factor(100)")             # not tested

    then at the (Pdb) prompt type ``s`` (or ``step``), then press :kbd:`Return`
    over and over to step through every line of Python that is called
    in the course of the above computation. Type ``?`` at any time for
    help on how to use the debugger (e.g., ``l`` lists 11 lines around
    the current line; ``bt`` gives a back trace, etc.).

    Setting a break point: If you have some code in a file and would
    like to drop into the debugger at a given point, put the following
    code at that point in the file:

    ``import pdb; pdb.set_trace()``

    For an article on how to use the Python debugger, see
    http://www.onlamp.com/pub/a/python/2005/09/01/debugger.html

    TESTS:

    For tests we disable garbage collection, see :issue:`21258` ::

        sage: import gc
        sage: gc.disable()

    The only real way to test this is via pexpect spawning a
    sage subprocess that uses IPython::

        sage: # needs pexpect sage.all
        sage: import pexpect
        sage: s = pexpect.spawn('sage')
        sage: _ = s.sendline("from sage.misc.trace import trace; trace('print(factor(10))'); print(3+97)")
        sage: _ = s.expect('ipdb>', timeout=90)
        sage: _ = s.sendline("s"); _ = s.sendline("c")
        sage: _ = s.expect('100', timeout=90)

    Seeing the ipdb prompt and the 2 \* 5 in the output below is a
    strong indication that the trace command worked correctly::

        sage: print(s.before[s.before.find(b'--'):].decode())                           # needs pexpect sage.all
        --...
        ...ipdb> c
        ...2 * 5...

    Re-enable garbage collection::

        sage: gc.enable()
    """
    from IPython.core.debugger import Pdb
    pdb = Pdb()

    try:
        ipython = get_ipython()
    except NameError:
        raise NotImplementedError("the trace command can only be run from the Sage command-line")

    from sage.repl.preparse import preparse
    code = preparse(code)
    return pdb.run(code, ipython.user_ns)

    # this could also be useful; it drops
    # us into a debugger in an except block:
    #     import pdb; pdb.post_mortem(sys.exc_info()[2])
