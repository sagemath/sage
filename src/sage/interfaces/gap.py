r"""
Interface to GAP

Sage provides an interface to the GAP system. This system provides
extensive group theory, combinatorics, etc.

The GAP interface will only work if GAP is installed on your
computer; this should be the case, since GAP is included with Sage.
The interface offers three pieces of functionality:


#. ``gap_console()`` -- a function that dumps you into
   an interactive command-line GAP session.

#. ``gap(expr)`` -- evaluation of arbitrary GAP
   expressions, with the result returned as a string.

#. ``gap.new(expr)`` -- creation of a Sage object that
   wraps a GAP object. This provides a Pythonic interface to GAP. For
   example, if ``f=gap.new(10)``, then
   ``f.Factors()`` returns the prime factorization of
   `10` computed using GAP.


First Examples
--------------

We factor an integer using GAP::

    sage: n = gap(20062006); n
    20062006
    sage: n.parent()
    Gap
    sage: fac = n.Factors(); fac
    [ 2, 17, 59, 73, 137 ]
    sage: fac.parent()
    Gap
    sage: fac[1]
    2

GAP and Singular
----------------

This example illustrates conversion between Singular and GAP via
Sage as an intermediate step. First we create and factor a Singular
polynomial.

::

    sage: singular(389)
    389
    sage: R1 = singular.ring(0, '(x,y)', 'dp')
    sage: f = singular('9*x^16-18*x^13*y^2-9*x^12*y^3+9*x^10*y^4-18*x^11*y^2+36*x^8*y^4+18*x^7*y^5-18*x^5*y^6+9*x^6*y^4-18*x^3*y^6-9*x^2*y^7+9*y^8')
    sage: F = f.factorize()
    sage: print(F)
    [1]:
       _[1]=9
       _[2]=x^6-2*x^3*y^2-x^2*y^3+y^4
       _[3]=-x^5+y^2
    [2]:
       1,1,2

Next we convert the factor `-x^5+y^2` to a Sage
multivariate polynomial. Note that it is important to let
`x` and `y` be the generators of a polynomial ring,
so the eval command works.

::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: s = F[1][3].sage_polystring(); s
    '-x**5+y**2'
    sage: g = eval(s); g
    -x^5 + y^2

Next we create a polynomial ring in GAP and obtain its
indeterminates::

    sage: R = gap.PolynomialRing('Rationals', 2); R
    PolynomialRing( Rationals, ["x_1", "x_2"] )
    sage: I = R.IndeterminatesOfPolynomialRing(); I
    [ x_1, x_2 ]

In order to eval `g` in GAP, we need to tell GAP to view
the variables ``x0`` and ``x1`` as the two
generators of `R`. This is the one tricky part. In the GAP
interpreter the object ``I`` has its own name (which
isn't ``I``). We can access its name using
``I.name()``.

::

    sage: _ = gap.eval("x := %s[1];; y := %s[2];;"%(I.name(), I.name()))

Now `x_0` and `x_1` are defined, so we can
construct the GAP polynomial `f` corresponding to
`g`::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: f = gap(str(g)); f
    -x_1^5+x_2^2

We can call GAP functions on `f`. For example, we evaluate
the GAP ``Value`` function, which evaluates `f`
at the point `(1,2)`.

::

    sage: f.Value(I, [1,2])
    3
    sage: g(1,2)        # agrees
    3

Saving and loading objects
--------------------------

Saving and loading GAP objects (using the dumps method, etc.) is
*not* supported, since the output string representation of Gap
objects is sometimes not valid input to GAP. Creating classes that
wrap GAP objects *is* supported, via simply defining the a
_gap_init_ member function that returns a string that when
evaluated in GAP constructs the object. See
``groups/perm_gps/permgroup.py`` for a nontrivial
example of this.

Long Input
----------

The GAP interface reads in even very long input (using files) in a
robust manner, as long as you are creating a new object.

.. NOTE::

   Using ``gap.eval`` for long input is much less robust, and is not
   recommended.

::

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = gap(t)

Changing which GAP is used, and how
-----------------------------------

Set the environment variable :envvar:`SAGE_GAP_COMMAND` to specify
how GAP executable is called. E.g.  ::

  $ SAGE_GAP_COMMAND = "/usr/local/bin/gap -s 4G" ./sage

will use GAP installed in `/usr/local/bin`, with 4Gb RAM.


Set the environment variable :envvar:`SAGE_GAP_MEMORY` to specify the amount
of RAM allocated to :mod:`~sage.libs.gap.libgap` and to the GAP executable.
If :envvar:`SAGE_GAP_COMMAND` is set, as well, then
:envvar:`SAGE_GAP_MEMORY` is only used for `libgap`. ::

    sage: gap.eval('GAPInfo.CommandLineOptions.s') # not tested
    '"42m"'

After the GAP interface initialisation, setting :envvar:`SAGE_GAP_MEMORY` has no effect::

    sage: os.environ['SAGE_GAP_MEMORY'] = '24M'
    sage: gap.eval('GAPInfo.CommandLineOptions.s') # not tested
    '"4g"'


AUTHORS:

- David Joyner and William Stein: initial version(s)

- William Stein (2006-02-01): modified gap_console command so it uses
  exactly the same startup command as Gap.__init__.

- William Stein (2006-03-02): added tab completions: gap.[tab], x =
  gap(...), x.[tab], and docs, e.g., gap.function? and x.function?
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .expect import Expect, ExpectElement, FunctionElement, ExpectFunction
from .gap_workspace import gap_workspace_file, prepare_workspace_dir
from sage.cpython.string import bytes_to_str
from sage.env import SAGE_EXTCODE, SAGE_GAP_COMMAND, SAGE_GAP_MEMORY, GAP_ROOT_PATHS
from sage.misc.misc import is_in_string
from sage.misc.cachefunc import cached_method
from sage.misc.instancedoc import instancedoc
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.structure.element import ModuleElement

import sage.interfaces.abc

import re
import os
import pexpect
import time
import platform
import string
import warnings

WORKSPACE = gap_workspace_file()

first_try = True

if SAGE_GAP_COMMAND is None:
    # Passing -A allows us to use a minimal GAP installation without
    # producing errors at start-up. The files sage.g and sage.gaprc are
    # used to load any additional packages that may be available.
    gap_cmd = f'gap -A -l "{GAP_ROOT_PATHS}"'
    if SAGE_GAP_MEMORY is not None:
        gap_cmd += " -s " + SAGE_GAP_MEMORY + " -o " + SAGE_GAP_MEMORY
else:
    gap_cmd = SAGE_GAP_COMMAND


if platform.processor() == 'ia64' and os.path.exists('/usr/bin/prctl'):
    # suppress unaligned access to 0x..., ip=0x... warnings
    gap_cmd = 'prctl --unaligned=silent ' + gap_cmd


def gap_command(use_workspace_cache=True, local=True):
    if use_workspace_cache:
        if local:
            return "%s -L %s" % (gap_cmd, WORKSPACE), False
        else:
            # TO DO: Use remote workspace
            return gap_cmd, False
    else:
        return gap_cmd, True


# ########### Classes with methods for both the GAP3 and GAP4 interface

class Gap_generic(ExtraTabCompletion, Expect):
    r"""
    Generic interface to the GAP3/GAP4 interpreters.

    AUTHORS:

    - William Stein and David Joyner (interface for GAP4)

    - Franco Saliola (Feb 2010): refactored to separate out the generic
      code
    """
    _identical_function = "IsIdenticalObj"

    def _synchronize(self, timeout=0.5, cmd='%s;'):
        """
        Synchronize GAP pexpect interface.

        See the base method
        :meth:`~sage.interfaces.expect.Expect._synchronize` for more
        details.

        We override this method since we are looking at GAP package
        mode output, which is quite different from the normal
        (human-readable) interface.

        EXAMPLES::

            sage: gap('"ok"')
            ok
            sage: gap._expect.sendline()  # now we are out of sync
            1
            sage: gap._synchronize()
            sage: gap(123)
            123
        """
        if self._expect is None:
            return
        E = self._expect
        from sage.misc.prandom import randrange
        rnd = randrange(2147483647)
        cmd = str(rnd) + ';'
        try:
            E.sendline(cmd)
            E.expect(r'@[nf][@J\s>]*' + str(rnd), timeout=timeout)
            E.send(' ')
            E.expect('@i', timeout=timeout)
        except pexpect.TIMEOUT:
            self.interrupt()
        except pexpect.EOF:
            self._crash_msg()
            self.quit()

    def interrupt(self, tries=None, timeout=1, quit_on_fail=True):
        """
        Interrupt the GAP process.

        Gap installs a SIGINT handler, we call it directly instead of
        trying to sent Ctrl-C. Unlike
        :meth:`~sage.interfaces.expect.Expect.interrupt`, we only try
        once since we are knowing what we are doing.

        Sometimes GAP dies while interrupting.

        EXAMPLES::

            sage: gap._eval_line('while(1=1) do i:=1;; od;', wait_for_prompt=False)
            ''
            sage: rc = gap.interrupt(timeout=1)
            sage: [ gap(i) for i in range(10) ]   # check that it is still working
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        TESTS::

            sage: gap('"finished computation"'); gap.interrupt(); gap('"ok"')
            finished computation
            True
            ok
        """
        E = self._expect
        if E is None:
            return True
        # GAP oddity: If a computation is running and we send Ctrl-C,
        # it is stopped as expected. But if we are at the idle prompt,
        # nothing is happening UNTIL we run the next command (which is
        # then immediately interrupted).
        # There is apparently also a race in GAP between the signal
        # handler and input, if we don't wait a bit the result is
        # unpredictable.
        E.sendline(chr(3))
        time.sleep(0.1)
        E.sendline()
        try:
            # send a dummy command
            E.sendline('224433409;')
            # read everything up to the actual output of the command
            E.expect(r'@[nf][@J\s>]*224433409', timeout=timeout)
            E.send(' ')
            # the following input prompt should be the current input
            # prompt but GAP might be too confused to display it
            #    E.expect('@i', timeout=timeout)
            # Ideally, we would be finished here. But sometimes GAP
            # thinks it is still inside a do/od block. So we run some
            # more plain commands to get back into sync. These might
            # either complete successfully (output "@n+<number>") or
            # return a "Syntax error: od expected@J@f +<number>"
            E.sendline()
            time.sleep(0.1)
            E.sendline('224433437;')
            E.expect(r'@[nf][@J\s>]*224433437', timeout=timeout)
            E.sendline()
            time.sleep(0.1)
            E.sendline('224433479;')
            E.expect(r'@[nf][@J\s>]*224433479', timeout=timeout)
            E.send(' ')
            # the following input prompt is now the current input prompt
            E.expect('@i', timeout=timeout)
            success = True
        except (pexpect.TIMEOUT, pexpect.EOF):
            # GAP died or hangs indefinitely
            success = False

        if not success and quit_on_fail:
            self.quit()
        return success

    def _assign_symbol(self):
        r"""
        Return the assign symbol in GAP.

        TESTS::

            sage: gap = Gap()
            sage: print(gap._assign_symbol())
            :=
        """
        return ":="

    def _quit_string(self):
        """
        Return the string used to quit GAP.

        EXAMPLES::

            sage: gap._quit_string()
            'quit;'

        ::

            sage: g = Gap()
            sage: a = g(2); g.is_running()
            True
            sage: g.quit()
            sage: g.is_running()
            False
        """
        return 'quit;'

    def _read_in_file_command(self, filename):
        r"""
        Return the command use to read in a file in GAP.

        EXAMPLES::

            sage: gap._read_in_file_command('test')
            'Read("test");'

        ::

            sage: filename = tmp_filename()
            sage: with open(filename, 'w') as f:
            ....:     _ = f.write('xx := 22;\n')
            sage: gap.read(filename)
            sage: gap.get('xx').strip()
            '22'
        """
        return 'Read("%s");' % filename

    def _continuation_prompt(self):
        """
        Return the continuation prompt in GAP.

        EXAMPLES::

            sage: gap._continuation_prompt()
            '> '
        """
        return '> '

    def load_package(self, pkg, verbose=False):
        """
        Load the Gap package with the given name.

        If loading fails, raise a :exc:`RuntimeError` exception.

        TESTS::

            sage: gap.load_package("chevie")
            Traceback (most recent call last):
            ...
            RuntimeError: Error loading Gap package chevie. You may want to install gap_packages SPKG.
        """
        if verbose:
            print("Loading GAP package {}".format(pkg))
        x = self.eval('LoadPackage("{}")'.format(pkg))
        if x == 'fail':
            raise RuntimeError("Error loading Gap package " + str(pkg) + ". " +
                               "You may want to install gap_packages SPKG.")

    def eval(self, x, newlines=False, strip=True, split_lines=True, **kwds):
        r"""
        Send the code in the string s to the GAP interpreter and return the
        output as a string.

        INPUT:

        - ``s`` -- string containing GAP code

        - ``newlines`` -- boolean (default: ``True``); if ``False``,
          remove all backslash-newlines inserted by the GAP output formatter

        - ``strip`` -- ignored

        - ``split_lines`` -- boolean (default: ``True``); if ``True`` then each
          line is evaluated separately.  If ``False``, then the whole
          block of code is evaluated all at once.

        EXAMPLES::

            sage: gap.eval('2+2')
            '4'
            sage: gap.eval('Print(4); #test\n Print(6);')
            '46'
            sage: gap.eval('Print("#"); Print(6);')
            '#6'
            sage: gap.eval('4; \n 6;')
            '4\n6'
            sage: gap.eval('if 3>2 then\nPrint("hi");\nfi;')
            'hi'
            sage: gap.eval('## this is a test\nPrint("OK")')
            'OK'
            sage: gap.eval('Print("This is a test. Oh no, a #");# but this is a comment\nPrint("OK")')
            'This is a test. Oh no, a #OK'
            sage: gap.eval('if 4>3 then')
            ''
            sage: gap.eval('Print("Hi how are you?")')
            'Hi how are you?'
            sage: gap.eval('fi')
            ''

        TESTS:

        Whitespace is not stripped from the front of the result
        (:issue:`28439`)::

            sage: gap.eval(r'Print("  -\n\\\\-  ")')
            '  -\n\\\\-'
        """
        # '"
        # We remove all of the comments:  On each line, we try
        # to find a pound sign.  If we find it, we check to see if
        # it is occurring in a string.  If it is not in a string, we
        # strip off the comment.
        if not split_lines:
            input_line = str(x)
        else:
            input_line = ""
            for line in str(x).rstrip().split('\n'):
                pound_position = line.find('#')
                while pound_position != -1:
                    if not is_in_string(line, pound_position):
                        line = line[:pound_position]
                    pound_position = line.find('#', pound_position + 1)
                input_line += " " + line
            if not input_line.endswith(';'):
                input_line += ';'
        result = Expect.eval(self, input_line, **kwds)
        if not newlines:
            result = result.replace("\\\n", "")
        return result.rstrip()

    def _execute_line(self, line, wait_for_prompt=True, expect_eof=False):
        if self._expect is None:  # interface is down
            self._start()
        E = self._expect
        try:
            if len(line) > 4095:
                raise RuntimeError("Passing commands this long to gap would hang")
            E.sendline(line)
        except OSError:
            raise RuntimeError("Error evaluating %s in %s" % (line, self))
        if not wait_for_prompt:
            return (b'', b'')
        if len(line) == 0:
            return (b'', b'')
        try:
            terminal_echo = []   # to be discarded
            normal_outputs = []  # GAP stdout
            error_outputs = []   # GAP stderr
            current_outputs = terminal_echo
            while True:
                x = E.expect_list(self._compiled_full_pattern)
                current_outputs.append(E.before)
                if x == 0:   # @p
                    if E.after != b'@p1.':
                        warnings.warn(
                            "possibly wrong version of GAP package "
                            "interface. Crossing fingers and continuing.")
                elif x == 1:  # @@
                    current_outputs.append(b'@')
                elif x == 2:  # special char
                    c = ord(E.after[1:2]) - ord(b'A') + 1
                    s = bytes([c])
                    current_outputs.append(s)
                elif x == 3:  # garbage collection info, ignore
                    pass
                elif x == 4:  # @e -- break loop
                    E.sendline("quit;")
                elif x == 5:  # @c completion, doesn't seem to happen when -p is in use
                    warnings.warn("I didn't think GAP could do this")
                elif x == 6:  # @f GAP error message
                    current_outputs = error_outputs
                elif x == 7:  # @h help text, but this stopped happening with new help
                    warnings.warn("I didn't think GAP could do this")
                elif x == 8:  # @i awaiting normal input
                    break
                elif x == 9:  # @m finished running a child
                    pass   # there is no need to do anything
                elif x == 10:  # @n normal output line
                    current_outputs = normal_outputs
                elif x == 11:  # @r echoing input
                    current_outputs = terminal_echo
                elif x == 12:  # @sN shouldn't happen
                    warnings.warn("this should never happen")
                elif x == 13:  # @w GAP is trying to send a Window command
                    warnings.warn("this should never happen")
                elif x == 14:  # @x seems to be safely ignorable
                    pass
                elif x == 15:  # @z GAP starting a subprocess
                    pass  # there is no need to do anything
        except pexpect.EOF:
            if not expect_eof:
                raise RuntimeError("Unexpected EOF from %s executing %s" % (self, line))
        except OSError:
            raise RuntimeError("IO Error from %s executing %s" % (self, line))
        return (b"".join(normal_outputs), b"".join(error_outputs))

    def _keyboard_interrupt(self):
        """
        TESTS:

        We check that the gap interface behaves correctly after an
        interrupt::

            sage: gap(2)
            2
            sage: try:
            ....:     alarm(0.5)
            ....:     gap.eval('while(1=1) do i:=1;; od;', wait_for_prompt=True)
            ....: except KeyboardInterrupt:
            ....:     pass
            sage: gap(2)
            2
        """
        self.quit()
        raise KeyboardInterrupt("Ctrl-c pressed while running %s" % self)

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=True):
        r"""
        Evaluate a line of commands.

        REMARK:

        By default, a long command (length exceeding ``self._eval_using_file_cutoff``)
        is evaluated using :meth:`_eval_line_using_file`.

        If the command can not be evaluated since the interface
        has crashed, it is automatically restarted and tried
        again *once*.

        If the optional ``wait_for_prompt`` is ``False`` then even a very long line
        will not be evaluated by :meth:`_eval_line_using_file`, since this does not
        support the ``wait_for_prompt`` option.

        INPUT:

        - ``line`` -- string; a command
        - ``allow_use_file`` -- boolean (default: ``True``);
          allow to evaluate long commands using :meth:`_eval_line_using_file`
        - ``wait_for_prompt`` -- boolean (default: ``True``);
          wait until the prompt appears in the sub-process' output
        - ``restart_if_needed`` -- boolean (default: ``True``);
          if it is ``True``, the command evaluation is evaluated
          a second time after restarting the interface, if an
          :exc:`EOFError` occurred.

        TESTS::

            sage: from sage.interfaces.gap import gap
            sage: gap._eval_line('2+2;')
            '4'

        We test the ``wait_for_prompt`` option by sending a command that
        creates an infinite loop in the GAP sub-process. But if we don't
        wait for the prompt to appear in the output, we can interrupt
        the loop without raising a KeyboardInterrupt. At the same time,
        we test that the line is not forwarded to :meth:`_eval_line_using_file`,
        since that method would not support the ``wait_for_prompt`` option::

            sage: cutoff = gap._eval_using_file_cutoff
            sage: gap._eval_using_file_cutoff = 4
            sage: gap._eval_line('while(1=1) do i:=1;; od;', wait_for_prompt=False)
            ''
            sage: rc = gap.interrupt(timeout=1)
            sage: gap._eval_using_file_cutoff = cutoff

        The following tests against a bug fixed at :issue:`10296`::

            sage: gap(3)
            3
            sage: gap.eval('quit;')
            ''
            sage: a = gap(3)
            ** Gap crashed or quit executing '\$sage...:=3;;' **
            Restarting Gap and trying again
            sage: a
            3
        """
        expect_eof = self._quit_string() in line

        try:
            if self._expect is None:
                self._start()
            if allow_use_file and wait_for_prompt and len(line) > self._eval_using_file_cutoff:
                return self._eval_line_using_file(line)

            (normal, error) = self._execute_line(line, wait_for_prompt=wait_for_prompt,
                                                 expect_eof=expect_eof)

            # The internal method _execute_line returns bytes but the bytes it
            # returns should contain text (any terminal commands and other
            # garbage should be filtered out by this point); here we decode
            # them (on Python 3), currently just using the default encoding
            normal, error = bytes_to_str(normal), bytes_to_str(error)

            if error:
                if 'Error, Rebuild completion files!' in error:
                    error += "\nRunning gap_reset_workspace()..."
                    self.quit()
                    gap_reset_workspace()
                error = error.replace('\r', '')
                raise RuntimeError("%s produced error output\n%s\n   executing %s" % (self, error, line))
            if not normal:
                return ''

            if isinstance(wait_for_prompt, str) and normal.ends_with(wait_for_prompt):
                n = len(wait_for_prompt)
            elif normal.endswith(bytes_to_str(self._prompt)):
                n = len(self._prompt)
            elif normal.endswith(self._continuation_prompt()):
                n = len(self._continuation_prompt())
            else:
                n = 0
            out = normal[:-n]
            if out and out[-1] == "\n":
                out = out[:-1]
            return out

        except (RuntimeError, TypeError, pexpect.ExceptionPexpect) as exc:
            if not self._isalive():
                # We can't distinguish just EOF from an unexpectedly killed
                # process because pexpect catches EOF's and re-reraises them
                # But if we *were* expecting EOF then we should just let it
                # fail silently and return
                if expect_eof:
                    return ''

                print("** %s crashed or quit executing '%s' **" % (self, line))
                print("Restarting %s and trying again" % self)
                self._start()
                if line != '':
                    return self._eval_line(line, allow_use_file=allow_use_file)
                else:
                    return ''
            else:
                raise RuntimeError(exc)

        except KeyboardInterrupt:
            self._keyboard_interrupt()
            raise KeyboardInterrupt("Ctrl-c pressed while running %s" % self)

    def unbind(self, var):
        """
        Clear the variable named var.

        EXAMPLES::

            sage: gap.set('x', '2')
            sage: gap.get('x')
            '2'
            sage: gap.unbind('x')
            sage: gap.get('x')
            Traceback (most recent call last):
            ...
            RuntimeError: Gap produced error output
            Error, Variable: 'x' must have a value
            ...
        """
        self.eval('Unbind(%s)' % var)
        self.clear(var)

    def _contains(self, v1, v2):
        """
        EXAMPLES::

            sage: Integers = gap('Integers')
            sage: two = gap(2)
            sage: gap._contains(two.name(), Integers.name())
            True

        ::

            sage: 2 in gap('Integers')
            True
        """
        return self.eval('%s in %s' % (v1, v2)) == "true"

    def _true_symbol(self):
        """
        Return the symbol for truth in GAP.

        EXAMPLES::

            sage: gap._true_symbol()
            'true'
            sage: gap(2) == gap(2)
            True
        """
        return "true"

    def _false_symbol(self):
        """
        Return the symbol for falsity in GAP.

        EXAMPLES::

            sage: gap._false_symbol()
            'false'
            sage: gap(2) == gap(3)
            False
        """
        return "false"

    def _equality_symbol(self):
        """
        Return the symbol for equality in GAP.

        EXAMPLES::

            sage: gap._equality_symbol()
            '='
            sage: gap(2) == gap(3)
            False
            sage: gap(2) == gap(2)
            True
        """
        return "="

    def version(self):
        """
        Return the version of GAP being used.

        EXAMPLES::

            sage: print(gap.version())
            4...
        """
        return self.eval('GAPInfo.Version')[1:-1]

    def function_call(self, function, args=None, kwds=None):
        """
        Call the GAP function with ``args`` and ``kwds``.

        EXAMPLES::

            sage: gap.function_call('SymmetricGroup', [5])
            SymmetricGroup( [ 1 .. 5 ] )

        If the GAP function does not return a value, but prints something
        to the screen, then a string of the printed output is returned.

        ::

            sage: s = gap.function_call('Display', [gap.SymmetricGroup(5).CharacterTable()])
            sage: type(s)
            <class 'sage.interfaces.interface.AsciiArtString'>
            sage: s.startswith('CT')
            True

        TESTS:

        If the function call is too long, two ``gap.eval`` calls are made
        since returned values from commands in a file cannot be handled
        properly::

            sage: g = Gap()
            sage: g.function_call("ConjugacyClassesSubgroups", sage.interfaces.gap.GapElement(g, 'SymmetricGroup(2)', name = 'a_variable_with_a_very_very_very_long_name')) # random
            [ ConjugacyClassSubgroups(SymmetricGroup( [ 1 .. 2 ] ),Group( () )),
              ConjugacyClassSubgroups(SymmetricGroup( [ 1 .. 2 ] ),Group( [ (1,2) ] )) ]

        When the command itself is so long that it warrants use of a temporary
        file to be communicated to GAP, this does not cause problems since
        the file will contain a single command::

            sage: g.function_call("ConjugacyClassesSubgroups", sage.interfaces.gap.GapElement(g, 'SymmetricGroup(2)', name = 'a_variable_with_a_name_so_very_very_very_long_that_even_by_itself_will_make_expect_use_a_file')) # random
            [ ConjugacyClassSubgroups(SymmetricGroup( [ 1 .. 2 ] ),Group( () )),
              ConjugacyClassSubgroups(SymmetricGroup( [ 1 .. 2 ] ),Group( [ (1,2) ] )) ]
        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)

        # Here we have to do some magic because not all GAP
        # functions return a value.  If you try to store their
        # results to a variable, then GAP will complain.  Thus, before
        # we evaluate the function, we make it so that the marker string
        # is in the 'last' variable in GAP.  If the function returns a
        # value, then that value will be in 'last', otherwise it will
        # be the marker.
        marker = '__SAGE_LAST__:="__SAGE_LAST__";;'
        cmd = "%s(%s);;" % (function, ",".join([s.name() for s in args] +
                                               [f'{key}={value.name()}'
                                                for key, value in kwds.items()]))
        if len(marker) + len(cmd) <= self._eval_using_file_cutoff:
            # We combine the two commands so we only run eval() once and the
            #   only output would be from the second command
            res = self.eval(marker + cmd)
        else:
            self.eval(marker)
            res = self.eval(cmd)
        if self.eval(self._identical_function + '(last,__SAGE_LAST__)') != 'true':
            return self.new('last2;')
        else:
            if res.strip():
                from sage.interfaces.interface import AsciiArtString
                return AsciiArtString(res)

    def get_record_element(self, record, name):
        r"""
        Return the element of a GAP record identified by ``name``.

        INPUT:

        - ``record`` -- a GAP record
        - ``name`` -- string

        OUTPUT: :class:`GapElement`

        EXAMPLES::

            sage: rec = gap('rec( a := 1, b := "2" )')
            sage: gap.get_record_element(rec, 'a')
            1
            sage: gap.get_record_element(rec, 'b')
            2

        TESTS::

            sage: rec = gap('rec( a := 1, b := "2" )')
            sage: type(gap.get_record_element(rec, 'a'))
            <class 'sage.interfaces.gap.GapElement'>
        """
        return self('%s.%s' % (record.name(), name))


# We need to inherit from ModuleElement to support
# sage.structure.coerce_actions.ModuleAction and it needs to be first
# in the MRO because extension types should always come first.
@instancedoc
class GapElement_generic(ModuleElement, ExtraTabCompletion, ExpectElement):
    r"""
    Generic interface to the GAP3/GAP4 interpreters.

    AUTHORS:

    - William Stein and David Joyner (interface for GAP4)

    - Franco Saliola (Feb 2010): refactored to separate out the generic
      code
    """
    def _add_(self, other):
        """
        EXAMPLES::

            sage: a = gap(1)
            sage: a + a
            2
        """
        # This is just a copy of ExpectElement._add_ to fix the fact
        # that the abtract method ModuleElement._add_ comes first in
        # the MRO.
        return self._operation("+", other)

    def __bool__(self):
        """
        EXAMPLES::

            sage: bool(gap(2))
            True
            sage: gap(0).bool()
            False
            sage: gap('false').bool()
            False
        """
        P = self._check_valid()
        return self != P(0) and repr(self) != 'false'

    def __len__(self):
        """
        EXAMPLES::

            sage: v = gap('[1,2,3]'); v
            [ 1, 2, 3 ]
            sage: len(v)
            3

        len is also called implicitly by if::

            sage: if gap('1+1 = 2'):
            ....:     print("1 plus 1 does equal 2")
            1 plus 1 does equal 2

        ::

            sage: if gap('1+1 = 3'):
            ....:     print("it is true")
            ....: else:
            ....:     print("it is false")
            it is false
        """
        P = self.parent()
        if P.eval('%s = true' % self.name()) == 'true':
            return 1
        elif P.eval('%s = false' % self.name()) == 'true':
            return 0
        else:
            return int(self.Length())

    def is_string(self):
        """
        Tell whether this element is a string.

        EXAMPLES::

            sage: gap('"abc"').is_string()
            True
            sage: gap('[1,2,3]').is_string()
            False
        """
        return bool(self.IsString())

    def _matrix_(self, R):
        r"""
        Return matrix over the (Sage) ring R determined by self, where self
        should be a Gap matrix.

        EXAMPLES::

            sage: s = gap("(Z(7)^0)*[[1,2,3],[4,5,6]]"); s
            [ [ Z(7)^0, Z(7)^2, Z(7) ], [ Z(7)^4, Z(7)^5, Z(7)^3 ] ]
            sage: s._matrix_(GF(7))
            [1 2 3]
            [4 5 6]

        ::

            sage: s = gap("[[1,2], [3/4, 5/6]]"); s
            [ [ 1, 2 ], [ 3/4, 5/6 ] ]
            sage: m = s._matrix_(QQ); m
            [  1   2]
            [3/4 5/6]
            sage: parent(m)
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        ::

            sage: s = gap('[[Z(16),Z(16)^2],[Z(16)^3,Z(16)]]')
            sage: s._matrix_(GF(16,'a'))
            [  a a^2]
            [a^3   a]
        """
        v = self.DimensionsMat()
        n = int(v[1])
        m = int(v[2])

        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(R, n, m)
        entries = [[R(self[r, c]) for c in range(1, m + 1)]
                   for r in range(1, n + 1)]
        return M(entries)


class Gap(Gap_generic):
    r"""
    Interface to the GAP interpreter.

    AUTHORS:

    - William Stein and David Joyner
    """
    def __init__(self, max_workspace_size=None,
                 maxread=None, script_subdirectory=None,
                 use_workspace_cache=True,
                 server=None,
                 server_tmpdir=None,
                 logfile=None,
                 seed=None,
                 env={}):
        """
        EXAMPLES::

            sage: gap == loads(dumps(gap))
            True
        """
        self.__use_workspace_cache = use_workspace_cache
        cmd, self.__make_workspace = gap_command(use_workspace_cache, server is None)
        # -b: suppress banner
        # -p: enable "package output mode"; this confusingly named option
        #     causes GAP to output special control characters that are normally
        #     intended for communication with a window manager (i.e. for xgap)
        #     but that we also use to control GAP with pexepect
        # -T: disable interactive break loop when encountering errors
        # -E: disable readline support
        cmd += " -b -p -T -E"
        cmd += ' -m 64m '   # attempt at a workaround for http://tracker.gap-system.org/issues/224
        cmd += ' ' + os.path.join(SAGE_EXTCODE, 'gap', 'sage.g')
        Expect.__init__(self,
                        name='gap',
                        prompt='gap> ',
                        command=cmd,
                        maxread=maxread,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=True,
                        verbose_start=False,
                        logfile=logfile,
                        eval_using_file_cutoff=100,
                        env=env)
        self.__seq = 0
        self._seed = seed

    def set_seed(self, seed=None):
        """
        Set the seed for gap interpreter.

        The seed should be an integer.

        EXAMPLES::

            sage: g = Gap()
            sage: g.set_seed(0)
            0
            sage: [g.Random(1,10) for i in range(5)]
            [2, 3, 3, 4, 2]
        """
        if seed is None:
            seed = self.rand_seed()
        self.eval("Reset(GlobalMersenneTwister,%d);;" % seed)
        self.eval("Reset(GlobalRandomSource,%d);;" % seed)
        self._seed = seed
        return seed

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.gap import gap
            sage: gap.__reduce__()
            (<function reduce_load_GAP at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            Gap
        """
        return reduce_load_GAP, tuple([])

    def _next_var_name(self):
        r"""
        Return the next unused variable name.

        Note that names starting with dollar signs are valid GAP
        identifiers, but need to be escaped with a backslash starting
        with GAP-4.8.

        EXAMPLES::

            sage: g = Gap()
            sage: g._next_var_name()
            '\\$sage1'
            sage: g(2)^2
            4
            sage: g._next_var_name()
            '\\$sage...'
        """
        if len(self._available_vars) != 0:
            v = self._available_vars[0]
            del self._available_vars[0]
            return v
        self.__seq += 1
        return r'\$sage%s' % self.__seq

    def _start(self):
        """
        EXAMPLES::

            sage: g = Gap()
            sage: g.is_running()
            False
            sage: g._start()
            sage: g.is_running()
            True
            sage: g.quit()
        """
        if self.__use_workspace_cache:
            from sage.libs.gap.saved_workspace import timestamp
            try:
                # Check to see if we need to auto-regenerate the gap
                # workspace, i.e., if the gap script is more recent
                # than the saved workspace, which signals that gap has
                # been upgraded.
                if os.path.getmtime(WORKSPACE) < timestamp():
                    raise OSError("GAP workspace too old")
                # Set the modification time of the workspace to the
                # current time.  This ensures the workspace doesn't
                # get deleted too soon by gap_reset_workspace().
                os.utime(WORKSPACE, None)
            except OSError:
                gap_reset_workspace(verbose=False)

        global first_try
        n = self._session_number
        try:
            Expect._start(self, "Failed to start GAP.")
        except Exception:
            if self.__use_workspace_cache and first_try:
                first_try = False
                self.quit()
                gap_reset_workspace(verbose=False)
                Expect._start(self, "Failed to start GAP.")
                self._session_number = n
                self.__make_workspace = False
            else:
                raise

        if self.__use_workspace_cache and self.__make_workspace:
            self.save_workspace()
        # Now, as self._expect exists, we can compile some useful pattern:
        self._compiled_full_pattern = self._expect.compile_pattern_list([
            r'@p\d+\.', '@@', '@[A-Z]', r'@[123456!"#$%&][^+]*\+',
            '@e', '@c', '@f', '@h', '@i', '@m', '@n', '@r', r'@s\d', r'@w.*\+', '@x', '@z'])
        # read everything up to the first "ready" prompt
        self._expect.expect("@i")

        # set random seed
        self.set_seed(self._seed)

    def _function_class(self):
        """
        Return the GapFunction class.

        EXAMPLES::

            sage: gap._function_class()
            <class 'sage.interfaces.gap.GapFunction'>

        ::

            sage: type(gap.Order)
            <class 'sage.interfaces.gap.GapFunction'>
        """
        return GapFunction

    def cputime(self, t=None):
        r"""
        Return the amount of CPU time that the GAP session has used. If
        ``t`` is not None, then it returns the difference
        between the current CPU time and ``t``.

        EXAMPLES::

            sage: t = gap.cputime()
            sage: t  #random
            0.13600000000000001
            sage: gap.Order(gap.SymmetricGroup(5))
            120
            sage: gap.cputime(t)  #random
            0.059999999999999998
        """
        if t is not None:
            return self.cputime() - t
        self.eval('_r_ := Runtimes();')
        r = sum(eval(self.eval('[_r_.user_time, _r_.system_time, _r_.user_time_children, _r_.system_time_children]')))
        return r / 1000.0

    def save_workspace(self):
        r"""
        Save the GAP workspace.

        TESTS:

        We make sure that :issue:`9938` (GAP does not start if the path
        to the GAP workspace file contains more than 82 characters) is
        fixed::

            sage: ORIGINAL_WORKSPACE = sage.interfaces.gap.WORKSPACE
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(prefix="0"*80) as f:    # long time (4s on sage.math, 2013)
            ....:     sage.interfaces.gap.WORKSPACE = f.name
            ....:     gap = Gap()
            ....:     gap('3+2')
            5
            sage: sage.interfaces.gap.WORKSPACE = ORIGINAL_WORKSPACE
        """
        prepare_workspace_dir()

        # According to the GAP Reference Manual,
        # [https://www.gap-system.org/Manuals/doc/htm/ref/CHAP003.htm#SSEC011.1]
        # SaveWorkspace can only be used at the main gap> prompt. It cannot
        # be included in the body of a loop or function, or called from a
        # break loop.
        from sage.misc.temporary_file import atomic_write
        with atomic_write(WORKSPACE) as f:
            f.close()
            self.eval('SaveWorkspace("%s");' % (f.name), allow_use_file=False)

    # Todo -- this -- but there is a tricky "when does it end" issue!
    # Maybe do via a file somehow?
    def help(self, s, pager=True):
        """
        Print help on a given topic.

        EXAMPLES:

        Note: In order to ensure consistent unicode handling from GAP we
        start a GAP instance with a forced UTF-8 locale::

            sage: gap = Gap(env={'LC_CTYPE': 'en_US.UTF-8'})
            sage: print(gap.help('SymmetricGroup', pager=False))
            <BLANKLINE>
              50.1-... SymmetricGroup
            <BLANKLINE>
              ‣ SymmetricGroup( [filt, ]deg ) ─────────────────────────────────── function
            ...
            <BLANKLINE>
        """
        if self.is_remote():
            tmp_to_use = self._remote_tmpfile()
        else:
            tmp_to_use = self._local_tmpfile()
        self.eval('SetGAPDocTextTheme("none")')
        gap_encoding = str(self('GAPInfo.TermEncoding;'))
        self.eval(r'\$SAGE.tempfile := "%s";' % tmp_to_use)
        line = Expect.eval(self, "? %s" % s)
        Expect.eval(self, "? 1")
        match = re.search(r"Page from (\d+)", line)
        if match is None:
            print(line)
        else:
            (sline,) = match.groups()
            sline = int(sline) - 1
            if self.is_remote():
                self._get_tmpfile()
            with open(self._local_tmpfile(),
                      encoding=gap_encoding) as fobj:
                help = fobj.read()
            if pager:
                from IPython.core.page import page
                page(help, start=sline)
            else:
                # Find the n-th line and return from there
                idx = -1
                while sline:
                    try:
                        idx = help.find('\n', idx + 1)
                        sline -= 1
                    except ValueError:
                        # We ran out of lines early somehow; this shouldn't
                        # happen though
                        break

                return help[idx:]

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: gap.set('x', '2')
            sage: gap.get('x')
            '2'
        """
        cmd = ('%s:=%s;;' % (var, value)).replace('\n', '')
        self._eval_line(cmd, allow_use_file=True)

    def get(self, var, use_file=False):
        """
        Get the string representation of the variable var.

        EXAMPLES::

            sage: gap.set('x', '2')
            sage: gap.get('x')
            '2'
        """
        if use_file:
            tmp = self._local_tmpfile()
            if os.path.exists(tmp):
                os.unlink(tmp)
            self.eval('PrintTo("%s", %s);' % (tmp, var), strip=False)
            with open(tmp) as f:
                r = f.read()
            r = r.strip().replace("\\\n", "")
            os.unlink(tmp)
            return r
        else:
            return self.eval('Print(%s);' % var, newlines=False)

    def _pre_interact(self):
        """
        EXAMPLES::

            sage: gap._pre_interact()
            sage: gap._post_interact()
        """
        self._eval_line(r'\$SAGE.StartInteract();')

    def _post_interact(self):
        """
        EXAMPLES::

            sage: gap._pre_interact()
            sage: gap._post_interact()
        """
        self._eval_line(r'\$SAGE.StopInteract();')

    def _eval_line_using_file(self, line):
        i = line.find(':=')
        if i != -1:
            j = line.find('"')
            if j >= 0 and j < i:
                i = -1
        if i == -1:
            line0 = 'Print( %s );' % line.rstrip().rstrip(';')
            try:  # this is necessary, since Print requires something as input, and some functions (e.g., Read) return nothing.
                return Expect._eval_line_using_file(self, line0)
            except RuntimeError:
                return ''
        return Expect._eval_line_using_file(self, line)

    def console(self):
        """
        Spawn a new GAP command-line session.

        EXAMPLES::

            sage: gap.console()  # not tested
            *********   GAP, Version 4.5.7 of 14-Dec-2012 (free software, GPL)
            *  GAP  *   https://www.gap-system.org
            *********   Architecture: x86_64-unknown-linux-gnu-gcc-default64
            Libs used:  gmp, readline
            Loading the library and packages ...
            Packages:   GAPDoc 1.5.1
            Try '?help' for help. See also  '?copyright' and  '?authors'
            gap>
        """
        gap_console()

    def _object_class(self):
        """
        Return the GapElement class.

        EXAMPLES::

            sage: gap._object_class()
            <class 'sage.interfaces.gap.GapElement'>
            sage: type(gap(2))
            <class 'sage.interfaces.gap.GapElement'>
        """
        return GapElement

    def _an_element_(self):
        """
        EXAMPLES::

            sage: gap._an_element_()
            0
        """
        return self(0)

    def _function_element_class(self):
        """
        Return the GapFunctionElement class.

        EXAMPLES::

            sage: gap._function_element_class()
            <class 'sage.interfaces.gap.GapFunctionElement'>
            sage: type(gap.SymmetricGroup(4).Order)
            <class 'sage.interfaces.gap.GapFunctionElement'>
        """
        return GapFunctionElement

    @cached_method
    def _tab_completion(self):
        """
        Return additional tab completion entries.

        OUTPUT: list of strings

        EXAMPLES::

            sage: '{}' in gap._tab_completion()
            False
            sage: c = gap._tab_completion()
            sage: len(c) > 100
            True
            sage: 'Order' in c
            True
        """
        names = eval(self.eval('NamesSystemGVars()')) + \
            eval(self.eval('NamesUserGVars()'))
        return [n for n in names if n[0] in string.ascii_letters]


def gap_reset_workspace(max_workspace_size=None, verbose=False):
    r"""
    Call this to completely reset the GAP workspace, which is used by
    default when Sage first starts GAP.

    The first time you start GAP from Sage, it saves the startup state
    of GAP in a file ``$HOME/.sage/gap/workspace-gap-HASH``, where
    ``HASH`` is a hash of the directory where Sage is installed. This
    is useful because the subsequent startup of GAP is at least ten
    times as fast. But if you update GAP or any of its packages, those
    changes won't take effect until the workspace is reset.

    TESTS:

    Check that the race condition from :issue:`14242` has been fixed.
    We temporarily need to change the worksheet filename, and to set
    ``first_try=True`` to ensure that the new workspace is created::

        sage: # long time
        sage: ORIGINAL_WORKSPACE = sage.interfaces.gap.WORKSPACE
        sage: saved_first_try = sage.interfaces.gap.first_try
        sage: sage.interfaces.gap.first_try = True
        sage: sage.interfaces.gap.WORKSPACE = tmp_filename()
        sage: from multiprocessing import Process
        sage: import time
        sage: gap = Gap()  # reset GAP session
        sage: P = [Process(target=gap, args=("14242",)) for i in range(4)]
        sage: for p in P:  # indirect doctest
        ....:     p.start()
        ....:     time.sleep(float(0.2))
        sage: for p in P:
        ....:     p.join()
        sage: os.unlink(sage.interfaces.gap.WORKSPACE)
        sage: sage.interfaces.gap.WORKSPACE = ORIGINAL_WORKSPACE
        sage: sage.interfaces.gap.first_try = saved_first_try
    """
    # The use_workspace_cache=False causes a new workspace to
    # be created, and we save it immediately thereafter.
    g = Gap(use_workspace_cache=False, max_workspace_size=None)
    g.save_workspace()
    g.quit()


@instancedoc
class GapElement(GapElement_generic, sage.interfaces.abc.GapElement):
    def __getitem__(self, n):
        """
        EXAMPLES::

            sage: a = gap([1,2,3])
            sage: a[1]
            1
        """
        self._check_valid()
        if not isinstance(n, tuple):
            return self.parent().new('%s[%s]' % (self._name, n))
        return self.parent().new('%s%s' % (self._name,
                                           ''.join('[%s]' % x for x in n)))

    def str(self, use_file=False):
        """
        EXAMPLES::

            sage: print(gap(2))
            2
        """
        if use_file:
            P = self._check_valid()
            return P.get(self.name(), use_file=True)
        else:
            return repr(self)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: s = gap("[[1,2], [3/4, 5/6]]")
            sage: latex(s)
            \left[\left[1, 2\right], \left[\frac{3}{4}, \frac{5}{6}\right]\right]
        """
        from sage.misc.latex import latex
        return latex(self._sage_())

    @cached_method
    def _tab_completion(self):
        """
        Return additional tab completion entries.

        OUTPUT: list of strings

        EXAMPLES::

            sage: s5 = gap.SymmetricGroup(5)
            sage: 'Centralizer' in s5._tab_completion()
            True
        """
        P = self.parent()
        v = P.eval(r'\$SAGE.OperationsAdmittingFirstArgument(%s)' % self.name())
        v = v.replace('Tester(', '').replace('Setter(', '').replace(')', '').replace('\n', '')
        v = v.split(',')
        v = (oper.split('"')[1] for oper in v)
        v = [oper for oper in v
             if all(ch in string.ascii_letters for ch in oper)]
        return sorted(set(v))


@instancedoc
class GapFunctionElement(FunctionElement):
    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: gap = Gap(env={'LC_CTYPE': 'en_US.UTF-8'})
            sage: print(gap(4).SymmetricGroup.__doc__)
            <BLANKLINE>
              50.1-... SymmetricGroup
            <BLANKLINE>
              ‣ SymmetricGroup( [filt, ]deg ) ─────────────────────────────────── function
            ...
        """
        M = self._obj.parent()
        help = M.help(self._name, pager=False)
        return help


@instancedoc
class GapFunction(ExpectFunction):
    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: gap = Gap(env={'LC_CTYPE': 'en_US.UTF-8'})
            sage: print(gap.SymmetricGroup.__doc__)
            <BLANKLINE>
              50.1-... SymmetricGroup
            <BLANKLINE>
              ‣ SymmetricGroup( [filt, ]deg ) ─────────────────────────────────── function
            ...
        """
        M = self._parent
        help = M.help(self._name, pager=False)
        return help


def is_GapElement(x):
    """
    Return ``True`` if ``x`` is a :class:`GapElement`.

    This function is deprecated; use :func:`isinstance`
    (of :class:`sage.interfaces.abc.GapElement`) instead.

    EXAMPLES::

        sage: from sage.interfaces.gap import is_GapElement
        sage: is_GapElement(gap(2))
        doctest:...: DeprecationWarning: the function is_GapElement is deprecated; use isinstance(x, sage.interfaces.abc.GapElement) instead
        See https://github.com/sagemath/sage/issues/34823 for details.
        True
        sage: is_GapElement(2)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(34823, "the function is_GapElement is deprecated; use isinstance(x, sage.interfaces.abc.GapElement) instead")
    return isinstance(x, GapElement)


def gfq_gap_to_sage(x, F):
    """
    INPUT:

    - ``x`` -- GAP finite field element

    - ``F`` -- Sage finite field

    OUTPUT: element of ``F``

    EXAMPLES::

        sage: x = gap('Z(13)')
        sage: F = GF(13, 'a')
        sage: F(x)
        2
        sage: F(gap('0*Z(13)'))
        0
        sage: F = GF(13^2, 'a')
        sage: x = gap('Z(13)')
        sage: F(x)
        2
        sage: x = gap('Z(13^2)^3')
        sage: F(x)
        12*a + 11
        sage: F.multiplicative_generator()^3
        12*a + 11

    TESTS:

    Check that :issue:`18048` is fixed::

        sage: K.<a> = GF(16)
        sage: b = a^2 + a
        sage: K(b._gap_())
        a^2 + a

    AUTHOR:

    - David Joyner and William Stein
    """
    s = str(x)
    if s[:2] == '0*':
        return F(0)
    i1 = s.index("(")
    i2 = s.index(")")
    q = eval(s[i1 + 1:i2].replace('^', '**'))
    if not F.cardinality().is_power_of(q):
        raise ValueError('%r has no subfield of size %r' % (F, q))
    if s.find(')^') == -1:
        e = 1
    else:
        e = int(s[i2 + 2:])
    if F.degree() == 1:
        g = F(gap.eval('Int(Z(%s))' % q))
    elif F.is_conway():
        f = (F.cardinality() - 1) // (q - 1)
        g = F.multiplicative_generator() ** f
    else:
        raise ValueError('%r is not prime or defined by a Conway polynomial' % F)
    return g**e


def intmod_gap_to_sage(x):
    r"""
    INPUT:

    - ``x`` -- Gap integer mod ring element

    EXAMPLES::

        sage: a = gap(Mod(3, 18)); a
        ZmodnZObj( 3, 18 )
        sage: b = sage.interfaces.gap.intmod_gap_to_sage(a); b
        3
        sage: b.parent()
        Ring of integers modulo 18

        sage: a = gap(Mod(3, 17)); a
        Z(17)
        sage: b = sage.interfaces.gap.intmod_gap_to_sage(a); b
        3
        sage: b.parent()
        Finite Field of size 17

        sage: a = gap(Mod(0, 17)); a
        0*Z(17)
        sage: b = sage.interfaces.gap.intmod_gap_to_sage(a); b
        0
        sage: b.parent()
        Finite Field of size 17

        sage: a = gap(Mod(3, 65537)); a
        ZmodpZObj( 3, 65537 )
        sage: b = sage.interfaces.gap.intmod_gap_to_sage(a); b
        3
        sage: b.parent()
        Ring of integers modulo 65537
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.rings.finite_rings.integer_mod import Mod
    from sage.rings.integer import Integer
    s = str(x)
    m = re.search(r'Z\(([0-9]*)\)', s)
    if m:
        return gfq_gap_to_sage(x, FiniteField(Integer(m.group(1))))
    m = re.match(r'Zmod[np]ZObj\( ([0-9]*), ([0-9]*) \)', s)
    if m:
        return Mod(Integer(m.group(1)), Integer(m.group(2)))
    raise ValueError("Unable to convert Gap element '%s'" % s)


gap = Gap()


def reduce_load_GAP():
    """
    Return the GAP interface object defined in ``sage.interfaces.gap``.

    EXAMPLES::

        sage: from sage.interfaces.gap import reduce_load_GAP
        sage: reduce_load_GAP()
        Gap
    """
    return gap


def gap_console():
    """
    Spawn a new GAP command-line session.

    Note that in gap-4.5.7 you cannot use a workspace cache that had
    no commandline to restore a gap session with commandline.

    EXAMPLES::

        sage: gap_console()  # not tested
        *********   GAP, Version 4.5.7 of 14-Dec-2012 (free software, GPL)
        *  GAP  *   https://www.gap-system.org
        *********   Architecture: x86_64-unknown-linux-gnu-gcc-default64
        Libs used:  gmp, readline
        Loading the library and packages ...
        Packages:   GAPDoc 1.5.1
        Try '?help' for help. See also  '?copyright' and  '?authors'
        gap>

    TESTS::

        sage: import subprocess as sp
        sage: from sage.interfaces.gap import gap_command
        sage: cmd = 'echo "quit;" | ' + gap_command(use_workspace_cache=False)[0]
        sage: gap_startup = sp.check_output(cmd, shell=True,
        ....:                               stderr=sp.STDOUT,
        ....:                               encoding='latin1')
        sage: 'www.gap-system.org' in gap_startup
        True
        sage: 'Error' not in gap_startup
        True
        sage: 'sorry' not in gap_startup
        True
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%gap magics instead.')
    cmd, _ = gap_command(use_workspace_cache=False)
    cmd += ' ' + os.path.join(SAGE_EXTCODE, 'gap', 'console.g')
    os.system(cmd)
