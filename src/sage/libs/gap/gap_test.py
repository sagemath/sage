import os
import select
import signal
import subprocess
import sys
import textwrap

import pytest
from sage.libs.gap.libgap import libgap


def test_initialization_preserves_sigchld(tmp_path):
    """Initialize in a fresh process with a host SIGCHLD handler installed."""
    script = textwrap.dedent('''
        import signal
        import subprocess
        import sys
        from sage.libs.gap.libgap import libgap
        from sage.libs.gap.saved_workspace import workspace

        assert workspace()[1] == bool(int(sys.argv[1]))

        received = []
        signal.signal(signal.SIGCHLD, lambda signum, frame: received.append(signum))
        assert libgap(1) == 1
        received.clear()

        # getsignal() only reads Python's cached handler, so send a real signal.
        signal.raise_signal(signal.SIGCHLD)
        assert received == [signal.SIGCHLD], received

        # The host must also retain ownership of its children's exit statuses.
        child = subprocess.run([sys.executable, '-c', 'raise SystemExit(23)'],
                               timeout=30, check=False)
        assert child.returncode == 23, child.returncode
    ''')
    # Exercise both a new workspace and loading the saved workspace.
    for cached in (0, 1):
        result = subprocess.run([sys.executable, '-c', script, str(cached)],
                                env={**os.environ, 'DOT_SAGE': str(tmp_path)},
                                capture_output=True, text=True, timeout=60,
                                check=False)
        assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('signum', [signal.SIGINT, signal.SIGALRM])
@pytest.mark.parametrize('interface', ['eval', 'function'])
def test_interrupt_restores_signal_handlers(signum, interface):
    """Interrupt GAP repeatedly and check recovery and the host's handlers."""
    script = textwrap.dedent(r'''
        import signal
        import sys
        from cysignals.signals import AlarmInterrupt
        from sage.libs.gap.libgap import libgap
        from sage.libs.gap.util import GAPError

        body = r'PrintTo("*stdout*", "ready\n"); while true do od;'
        spin = libgap.eval('function() ' + body + ' end')
        received = []
        for sig in (signal.SIGINT, signal.SIGALRM):
            signal.signal(sig, lambda signum, frame: received.append(signum))

        def check_handlers():
            received.clear()
            for sig in (signal.SIGINT, signal.SIGALRM):
                signal.raise_signal(sig)
            assert received == [signal.SIGINT, signal.SIGALRM], received
            received.clear()

        expected = KeyboardInterrupt if int(sys.argv[1]) == signal.SIGINT else AlarmInterrupt
        for _ in range(3):
            try:
                if sys.argv[2] == 'eval':
                    libgap.eval('if true then ' + body + ' fi')
                else:
                    spin()
            except expected:
                pass
            else:
                raise AssertionError('GAP was not interrupted')
            assert not received, received
            check_handlers()

            assert libgap.eval('1 + 1') == 2
            check_handlers()
            try:
                libgap.eval('Error("ordinary GAP error")')
            except GAPError as error:
                assert 'ordinary GAP error' in str(error), error
            else:
                raise AssertionError('GAPError was not raised')
            check_handlers()
    ''')
    # Wait until GAP itself announces that it is running before interrupting.
    # A Python thread cannot send the signal while libgap holds the GIL.
    process = subprocess.Popen([sys.executable, '-u', '-c', script,
                                str(int(signum)), interface],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True)
    try:
        for _ in range(3):
            ready, _, _ = select.select([process.stdout], [], [], 30)
            assert ready, 'GAP did not reach the interruptible computation'
            line = process.stdout.readline()
            assert line == 'ready\n', line or process.communicate(timeout=30)[1]
            process.send_signal(signum)
        stdout, stderr = process.communicate(timeout=30)
        assert process.returncode == 0, stdout + stderr
    finally:
        if process.poll() is None:
            process.kill()
        process.communicate()


def test_libgap_can_read_and_write_files(tmpfile):
    """
    Test that libgap can write to a temporary file and
    subsequently read it.

    See :issue:`16502`, :issue:`15833`.
    """
    message = "Ceci n'est pas une groupe"
    libgap.PrintTo(tmpfile.name, message)
    with open(tmpfile.name) as f:
        contents = f.read()
    assert contents == message
    SystemFile = libgap.function_factory('StringFile')
    assert SystemFile(tmpfile.name) == libgap(message)


@pytest.mark.long
def test_gc_loop_1():
    r"""
    Stress test for garbage collection in libgap.

    Manually runs the GAP garbage collector, and then creates 10,000
    instances of the cyclic group of order two in a loop. In each
    iteration, the python variable is overwritten, meaning that python
    is free to garbage collect the object.
    """
    libgap.collect()
    for _ in range(10000):
        libgap.CyclicGroup(2)
    assert True


@pytest.mark.long
def test_gc_loop_2():
    r"""
    Stress test for garbage collection in libgap.

    Create the free group on two elements (``a`` and ``b``) and then
    construct a quotient group of order two in a loop by specifying
    some relations. The python variables are overwritten in each
    iteration, meaning that python is free to garbage-collect them.
    (We also ensure that the quotient group has the expected order.)

    After that loop, we take one of the generators of the quotient
    group (from the final iteration), and compute its order in a
    loop. The python reference from the final iteration lives on,
    so this generator should not be collected.
    """
    G = libgap.FreeGroup(2)
    a, b = G.GeneratorsOfGroup()
    two = libgap(2)

    for _ in range(100):
        rel = libgap([a**2, b**2, a*b*a*b])
        H = G / rel
        H1 = H.GeneratorsOfGroup()[0]
        n = H1.Order()
        assert n == two

    result = True
    for i in range(300000):
        n = libgap.Order(H1)
        result &= (n == two)
    assert result


@pytest.mark.long
def test_gc_loop_3():
    r"""
    Stress test for garbage collection in libgap.

    Create the free group on two elements (``a`` and ``b``) and then
    add some of its elements to a list in a loop. A new, empty list is
    created at each iteration, so this mainly serves to guarantee that
    the generators of the group are not garbage-collected.
    """
    G = libgap.FreeGroup(2)
    a, b = G.GeneratorsOfGroup()
    for _ in range(300000):
        lis = libgap([])
        lis.Add(a ** 2)
        lis.Add(b ** 2)
        lis.Add(b * a)
    assert True
