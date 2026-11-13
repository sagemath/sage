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


@pytest.fixture
def no_pexpect_gap(monkeypatch):
    """Forbid both starting pexpect GAP and using its help interface."""
    from sage.interfaces.gap import Gap

    def no_pexpect(*args, **kwds):
        raise AssertionError('help must not start pexpect GAP')

    monkeypatch.setattr(Gap, '_start', no_pexpect)
    monkeypatch.setattr(Gap, 'help', no_pexpect)


@pytest.mark.parametrize('name', ['CyclicGroup', 'SymmetricGroup', 'Size', 'Socle'])
def test_function_help_uses_libgap(name, no_pexpect_gap):
    """Core and ambiguous help queries must not start pexpect GAP."""
    doc = getattr(libgap, name).__doc__
    assert name in doc
    assert '\x1b' not in doc
    assert f'{name}( ' in getattr(libgap, name).help(1, pager=False)


@pytest.mark.parametrize('name', ['Size', 'HasSize', 'SetIsFinite'])
def test_function_help_menu_matches_gap(name, no_pexpect_gap):
    """Use the same menu entries, order and line widths as GAP itself."""
    capture = libgap.eval('''function()
        local lines;
        lines := fail;
        return rec(pager := function(data) lines := data.lines; end,
                   result := function() return lines; end);
        end''')()
    screen = libgap.SizeScreen()
    try:
        libgap.SizeScreen([80, 24])
        with libgap.global_context('Pager', capture['pager']), \
             libgap.global_context('HELP_LAST', libgap.StructuralCopy(libgap.HELP_LAST)):
            libgap.HELP_SHOW_MATCHES(libgap.HELP_KNOWN_BOOKS[0],
                                    libgap.SIMPLE_STRING(name), True)
        expected = [libgap.StripEscapeSequences(line).sage()
                    for line in capture['result']()]
    finally:
        libgap.SizeScreen(screen)
    menu = getattr(libgap, name).__doc__
    assert 'several entries match' in menu
    assert 'help(n)' in menu
    assert menu.splitlines()[2:] == expected[1:]


def test_function_help_menu_selection(monkeypatch, no_pexpect_gap):
    """Selecting a menu entry must retrieve that entry, not always the first."""
    matches = libgap.eval('''function()
        local book;
        book := rec(bookname := "Example", handler := "SageHelpTest",
                    entries := [["First"], ["Second"], ["Third"]]);
        return [[[book, 1], [book, 2]], [[book, 3]]];
        end''')()
    monkeypatch.setattr(libgap, 'HELP_GET_MATCHES', lambda *args: matches)
    handler = {'SageHelpTest': {'HelpData': lambda book, entry, kind: libgap(f'Entry {entry}')}}
    monkeypatch.setattr(libgap, 'HELP_BOOK_HANDLER', handler)
    f = libgap.CyclicGroup
    assert f.__doc__.splitlines()[2:] == ['[1] Example: First',
                                         '[2] Example: Second',
                                         '[3] Example: Third']
    assert f.help(2, pager=False) == 'Entry 2'
    assert f.help(3, pager=False) == 'Entry 3'
    with pytest.raises(ValueError, match='help match number out of range'):
        f.help(4, pager=False)

    # GAP displays one exact match immediately, even with partial matches.
    monkeypatch.setattr(libgap, 'HELP_GET_MATCHES',
                        lambda *args: libgap([[matches[0][1]], matches[1]]))
    assert f.__doc__ == 'Entry 2'


def test_function_help_uses_sage_pager(monkeypatch, no_pexpect_gap):
    """Display the menu and the selected documentation through Sage's pager."""
    from IPython.core import page

    shown = []
    monkeypatch.setattr(page, 'page', shown.append)
    assert libgap.HasSize.help() is None
    assert 'several entries match' in shown.pop()
    assert libgap.HasSize.help(1) is None
    assert 'Size( ' in shown.pop()


def test_function_help_missing(capfd, no_pexpect_gap):
    """Missing documentation should return None without printing a message."""
    f = libgap.eval('x -> x')
    assert f.__doc__ is None
    libgap.SetNameFunction(f, 'SageUndocumentedFunctionHelpTest')
    assert f.__doc__ is None
    assert capfd.readouterr() == ('', '')


def test_function_help_package(no_pexpect_gap):
    """Use manuals from packages loaded in the current libgap session."""
    from sage.features.gap import GapPackage

    if not GapPackage('guava').is_present():
        pytest.skip('GUAVA is not available')

    assert 'RandomLinearCode' in libgap.RandomLinearCode.__doc__


@pytest.mark.parametrize('encoding', ['UTF-8', 'ISO-8859-1', 'ANSI_X3.4-1968'])
def test_function_help_encoding(encoding, no_pexpect_gap):
    """Decode help using GAP's terminal encoding."""
    original = libgap.GAPInfo['TermEncoding']
    set_encoding = libgap.eval('function(s) GAPInfo.TermEncoding := s; end')
    try:
        set_encoding(encoding)
        doc = libgap.SymmetricGroup.__doc__
        assert 'SymmetricGroup' in doc
        if encoding == 'UTF-8':
            assert '‣ SymmetricGroup' in doc
    finally:
        set_encoding(original)


@pytest.mark.parametrize('mode', ['single', 'menu', 'selection', 'error'])
def test_function_help_preserves_state(mode, monkeypatch, no_pexpect_gap):
    """Looking up documentation must preserve interactive help settings."""
    state = libgap.eval('''function()
        return rec(last := StructuralCopy(HELP_LAST), index := HELP_RING_IDX,
                   books := StructuralCopy(HELP_BOOK_RING),
                   topics := StructuralCopy(HELP_TOPIC_RING),
                   original_topics := StructuralCopy(HELP_ORIG_TOPIC_RING),
                   theme := StructuralCopy(GAPDocTextTheme),
                   viewers := UserPreference("HelpViewers"),
                   screen := SizeScreen(), encoding := GAPInfo.TermEncoding);
        end''')
    before = state()
    if mode == 'error':
        from sage.libs.gap.util import GAPError

        lookup_error = libgap.eval('function(a, b, c) Error("help lookup failed"); end')
        monkeypatch.setattr(libgap, 'HELP_GET_MATCHES', lookup_error)
        with pytest.raises(GAPError, match='help lookup failed'):
            libgap.CyclicGroup._instancedoc_()
    elif mode == 'single':
        assert libgap.CyclicGroup.__doc__
    elif mode == 'menu':
        assert 'several entries match' in libgap.HasSize.__doc__
    else:
        assert libgap.HasSize.help(1, pager=False)
    assert state() == before


@pytest.mark.parametrize('data, expected', [
    ('"first\\nsecond"', 'first\nsecond'),
    ('["first", "second"]', 'first\nsecond'),
    ('rec(lines := "first\\nsecond", start := 2)', 'second'),
    ('rec(lines := ["first", "second"], start := 2)', 'second'),
    ('rec(lines := "first\\nsecond")', 'first\nsecond'),
    ('fail', None),
])
def test_function_help_text_formats(data, expected, monkeypatch, no_pexpect_gap):
    """Help book handlers can return strings, lists of lines, or pager records."""
    matches = libgap.eval('[[], [[rec(handler := "SageHelpTest"), 1]]]')
    monkeypatch.setattr(libgap, 'HELP_GET_MATCHES', lambda *args: matches)
    handler = {'SageHelpTest': {'HelpData': lambda *args: libgap.eval(data)}}
    monkeypatch.setattr(libgap, 'HELP_BOOK_HANDLER', handler)
    assert libgap.CyclicGroup.__doc__ == expected


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
