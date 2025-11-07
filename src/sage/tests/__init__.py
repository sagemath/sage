import os
import select
import sys
from subprocess import PIPE, Popen


def check_executable(args, input='', timeout=100.0, pydebug_ignore_warnings=False, **kwds):
    r"""
    Run the program defined by ``args`` using the string ``input`` on
    the standard input.

    INPUT:

    - ``args`` -- list of program arguments, the first being the
      executable

    - ``input`` -- string serving as standard input; usually, this
      should end with a newline

    - ``timeout`` -- if the program produces no output for ``timeout``
      seconds, a :exc:`RuntimeError` is raised

    - ``pydebug_ignore_warnings`` -- boolean. Set the PYTHONWARNINGS environment variable to ignore
      Python warnings when on a Python debug build (`--with-pydebug`, e.g. from building with
      `SAGE_DEBUG=yes`). Debug builds do not install the default warning filters, which can break
      some doctests. Unfortunately the environment variable does not support regex message filters,
      so the filter will catch a bit more than the default filters. Hence we only enable it on debug
      builds.

    - ``**kwds`` -- additional keyword arguments passed to the
      :class:`Popen` constructor

    OUTPUT: a tuple ``(out, err, ret)`` with the standard output,
    standard error and exitcode of the program run.

    EXAMPLES::

        sage: from sage.tests import check_executable
        sage: (out, err, ret) = check_executable(["cat"], "Hello World!")
        sage: out
        'Hello World!'
        sage: err
        ''
        sage: ret
        0

    We test the timeout option::

        sage: (out, err, ret) = check_executable(["sleep", "1"], timeout=0.1)
        Traceback (most recent call last):
        ...
        RuntimeError: timeout in check_executable()
    """
    pexpect_env = dict(os.environ)
    try:
        del pexpect_env["TERM"]
    except KeyError:
        pass

    __with_pydebug = hasattr(sys, 'gettotalrefcount')   # This is a Python debug build (--with-pydebug)
    if __with_pydebug and pydebug_ignore_warnings:
        pexpect_env['PYTHONWARNINGS'] = ','.join([
            'ignore::DeprecationWarning',
        ])

    kwds['encoding'] = kwds.pop('encoding', 'utf-8')

    p = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, env=pexpect_env,
              **kwds)
    if input:
        p.stdin.write(input)

    p.stdin.close()
    fdout = p.stdout.fileno()
    fderr = p.stderr.fileno()
    out = []
    err = []

    while True:
        # Try reading from fdout and fderr
        rfd = []
        if fdout:
            rfd.append(fdout)
        if fderr:
            rfd.append(fderr)
        if len(rfd) == 0:
            break
        timeout = float(timeout)
        rlist = select.select(rfd, [], [], timeout)[0]

        if len(rlist) == 0:
            # Timeout!
            p.terminate()
            raise RuntimeError("timeout in check_executable()")
        if fdout in rlist:
            s = p.stdout.read(1024)
            if not s:
                fdout = None   # EOF
                p.stdout.close()
            out.append(s)
        if fderr in rlist:
            s = p.stderr.read(1024)
            if not s:
                fderr = None   # EOF
                p.stderr.close()
            err.append(s)

    # In case out or err contains a quoted string, force the use of
    # double quotes so that the output is enclosed in single
    # quotes. This avoids some doctest failures with some versions of
    # OS X and Xcode.
    out = ''.join(out)
    out = out.replace("'", '"')
    err = ''.join(err)
    err = err.replace("'", '"')

    return (out, err, p.wait())
