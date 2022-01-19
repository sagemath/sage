r"""
Interface to Maple

AUTHORS:

- William Stein (2005): maple interface

- Gregg Musiker (2006-02-02): tutorial

- William Stein (2006-03-05): added tab completion, e.g., maple.[tab],
  and help, e.g, maple.sin?.

You must have the optional commercial Maple interpreter installed
and available as the command ``maple`` in your PATH in
order to use this interface. You do not have to install any
optional Sage packages.

Type ``maple.[tab]`` for a list of all the functions
available from your Maple install. Type
``maple.[tab]?`` for Maple's help about a given
function. Type ``maple(...)`` to create a new Maple
object, and ``maple.eval(...)`` to run a string using
Maple (and get the result back as a string).

EXAMPLES::

    sage: maple('3 * 5')                                 # optional - maple
    15
    sage: maple.eval('ifactor(2005)')                    # optional - maple
    '``(5)*``(401)'
    sage: maple.ifactor(2005)                            # optional - maple
    ``(5)*``(401)
    sage: maple.fsolve('x^2=cos(x)+4', 'x=0..5')         # optional - maple
    1.914020619
    sage: maple.factor('x^5 - y^5')                      # optional - maple
    (x-y)*(x^4+x^3*y+x^2*y^2+x*y^3+y^4)

If the string "error" (case insensitive) occurs in the output of
anything from Maple, a RuntimeError exception is raised.

Tutorial
--------

AUTHORS:

- Gregg Musiker (2006-02-02): initial version.

This tutorial is based on the Maple Tutorial for number theory from
http://www.math.mun.ca/~drideout/m3370/numtheory.html.

There are several ways to use the Maple Interface in Sage. We will
discuss two of those ways in this tutorial.


#. If you have a maple expression such as

   ::

       factor( (x^5-1));

   We can write that in sage as

   ::

       sage: maple('factor(x^5-1)')                 # optional - maple
       (x-1)*(x^4+x^3+x^2+x+1)

   Notice, there is no need to use a semicolon.

#. Since Sage is written in Python, we can also import maple
   commands and write our scripts in a Pythonic way. For example,
   ``factor()`` is a maple command, so we can also factor
   in Sage using

   ::

       sage: maple('(x^5-1)').factor()              # optional - maple
       (x-1)*(x^4+x^3+x^2+x+1)

   where ``expression.command()`` means the same thing as
   ``command(expression)`` in Maple. We will use this
   second type of syntax whenever possible, resorting to the first
   when needed.

   ::

       sage: maple('(x^12-1)/(x-1)').simplify()     # optional - maple
       x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1


The normal command will always reduce a rational function to the
lowest terms. The factor command will factor a polynomial with
rational coefficients into irreducible factors over the ring of
integers. So for example,

::

    sage: maple('(x^12-1)').factor( )           # optional - maple
    (x-1)*(x+1)*(x^2+x+1)*(x^2-x+1)*(x^2+1)*(x^4-x^2+1)

::

    sage: maple('(x^28-1)').factor( )           # optional - maple
    (x-1)*(x^6+x^5+x^4+x^3+x^2+x+1)*(x+1)*(x^6-x^5+x^4-x^3+x^2-x+1)*(x^2+1)*(x^12-x^10+x^8-x^6+x^4-x^2+1)

Another important feature of maple is its online help. We can
access this through sage as well. After reading the description of
the command, you can press q to immediately get back to your
original prompt.

Incidentally you can always get into a maple console by the
command ::

    sage: maple.console()          # not tested
    sage: !maple                   # not tested

Note that the above two commands are slightly different, and the
first is preferred.

For example, for help on the maple command fibonacci, we type

::

    sage: maple.help('fibonacci')  # not tested, since it uses a pager

We see there are two choices. Type

::

    sage: maple.help('combinat, fibonacci')   # not tested, since it uses a pager

We now see how the Maple command fibonacci works under the
combinatorics package. Try typing in

::

    sage: maple.fibonacci(10)                # optional - maple
    fibonacci(10)

You will get fibonacci(10) as output since Maple has not loaded the
combinatorics package yet. To rectify this type

::

    sage: maple('combinat[fibonacci]')(10)     # optional - maple
    55

instead.

If you want to load the combinatorics package for future
calculations, in Sage this can be done as

::

    sage: maple.with_package('combinat')       # optional - maple

or

::

    sage: maple.load('combinat')               # optional - maple

Now if we type ``maple.fibonacci(10)``, we get the
correct output::

    sage: maple.fibonacci(10)                  # optional - maple
    55

Some common maple packages include ``combinat``,
``linalg``, and ``numtheory``. To produce
the first 19 Fibonacci numbers, use the sequence command.

::

    sage: maple('seq(fibonacci(i),i=1..19)')     # optional - maple
    1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584,
    4181

Two other useful Maple commands are ifactor and isprime. For
example

::

    sage: maple.isprime(maple.fibonacci(27))     # optional - maple
    false
    sage: maple.ifactor(maple.fibonacci(27))     # optional - maple
    ``(2)*``(17)*``(53)*``(109)

Note that the isprime function that is included with Sage (which
uses PARI) is better than the Maple one (it is faster and gives a
provably correct answer, whereas Maple is sometimes wrong).

::

    sage: alpha = maple('(1+sqrt(5))/2')         # optional - maple
    sage: beta = maple('(1-sqrt(5))/2')          # optional - maple
    sage: f19  = alpha^19 - beta^19/maple('sqrt(5)')      # optional - maple
    sage: f19                                             # optional - maple
    (1/2+1/2*5^(1/2))^19-1/5*(1/2-1/2*5^(1/2))^19*5^(1/2)
    sage: f19.simplify()                # somewhat randomly ordered output; optional - maple
    6765+5778/5*5^(1/2)

Let's say we want to write a maple program now that squares a
number if it is positive and cubes it if it is negative. In maple,
that would look like

::

    mysqcu := proc(x)
    if x > 0 then x^2;
    else x^3; fi;
    end;

In Sage, we write

::

    sage: mysqcu = maple('proc(x) if x > 0 then x^2 else x^3 fi end')    # optional - maple
    sage: mysqcu(5)                                                      # optional - maple
    25
    sage: mysqcu(-5)                                                     # optional - maple
    -125

More complicated programs should be put in a separate file and
loaded.
"""

#############################################################################
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#############################################################################

import os

from .expect import Expect, ExpectElement, ExpectFunction, FunctionElement, gc_disabled

import pexpect

from sage.env import DOT_SAGE
from sage.misc.pager import pager
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.docs.instancedoc import instancedoc
from sage.structure.richcmp import rich_to_bool


COMMANDS_CACHE = '%s/maple_commandlist_cache.sobj' % DOT_SAGE


class Maple(ExtraTabCompletion, Expect):
    """
    Interface to the Maple interpreter.

    Type ``maple.[tab]`` for a list of all the functions
    available from your Maple install. Type
    ``maple.[tab]?`` for Maple's help about a given
    function. Type ``maple(...)`` to create a new Maple
    object, and ``maple.eval(...)`` to run a string using
    Maple (and get the result back as a string).
    """
    def __init__(self, maxread=None, script_subdirectory=None, server=None,
            server_tmpdir=None, logfile=None, ulimit=None):
        """
        Create an instance of the Maple interpreter.

        EXAMPLES::

            sage: maple == loads(dumps(maple))
            True
        """
        __maple_iface_opts = [
            'screenwidth=infinity',
            'errorcursor=false']
        __maple_command = 'maple -t -c "interface({})"'.format(
            ','.join(__maple_iface_opts))
        # errorcursor=false avoids maple command line interface to dump
        # into the editor when an error occurs. Thus pexpect interface
        # is not messed up if a maple error occurs.
        # screenwidth=infinity prevents maple command interface from cutting
        # your input lines. By doing this, file interface also works in the
        # event that  sage_user_home + sage_tmp_file_stuff exceeds the
        # length of 79 characters.
        Expect.__init__(self,
                        name='maple',
                        prompt='#-->',
                        command=__maple_command,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        ulimit=ulimit,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=False,
                        verbose_start=False,
                        logfile=logfile,
                        eval_using_file_cutoff=2048)  # 2048 is
        # a small enough value to avoid conflicts with the 4096 limit
        # hardcoded in Expect.

    def _function_class(self):
        """
        EXAMPLES::

            sage: maple._function_class()
            <class 'sage.interfaces.maple.MapleFunction'>

        ::

            sage: type(maple.diff)
            <class 'sage.interfaces.maple.MapleFunction'>
        """
        return MapleFunction

    def _keyboard_interrupt(self):
        """
        EXAMPLES::

            sage: maple._keyboard_interrupt()   # not tested
            Interrupting Maple...
            ...
            RuntimeError: Ctrl-c pressed while running Maple

        ::

            sage: maple('Matrix(8000,8000)')    # not tested
            #Press ctrl-c
            ^CInterrupting Maple...
            ...
            RuntimeError: Ctrl-c pressed while running Maple
        """
        print("Interrupting %s..." % self)
        self._expect.sendline(chr(3))  # send ctrl-c
        self._expect.expect(self._prompt)
        raise RuntimeError("Ctrl-c pressed while running %s" % self)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: maple.__reduce__()
            (<function reduce_load_Maple at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            Maple
        """
        return reduce_load_Maple, tuple([])

    def _read_in_file_command(self, filename):
        r"""
        Return the string used to read filename into Maple.

        EXAMPLES::

            sage: maple._read_in_file_command('test')
            'read "test"'

        ::

            sage: filename = tmp_filename()  # optional - maple
            sage: with open(filename, 'w') as f:   # optional - maple
            ....:     _ = f.write('xx := 22;\n')
            sage: maple.read(filename)    # optional - maple
            sage: maple.get('xx').strip() # optional - maple
            '22'
        """
        return 'read "%s"' % filename

    def _quit_string(self):
        """
        EXAMPLES::

            sage: maple._quit_string()
            'quit'

        ::

            sage: m = Maple()        # optional - maple
            sage: a = m(2)           # optional - maple
            sage: m.is_running()     # optional - maple
            True
            sage: m.quit()           # optional - maple
            sage: m.is_running()     # optional - maple
            False
        """
        return 'quit'

    def _install_hints(self):
        """
        Hints for installing Maple on your computer.

        AUTHORS:

        - William Stein and Justin Walker (2006-02-12).

        EXAMPLES::

            sage: print(maple._install_hints())
            In order...
        """
        return """

In order to use the Maple interface you need to have Maple installed
and have a script in your PATH called "maple" that runs the
command-line version of Maple.  Alternatively, you could use a remote
connection to a server running Maple; for hints, type
    print(maple._install_hints_ssh())

  (1) You might have to buy Maple (http://webstore.maplesoft.com/).

  (2) * LINUX: The maple script comes standard with your Maple install.

      * APPLE OS X:
          (a) create a file called maple (in your PATH), with the following contents:
             #!/bin/sh
             /Library/Frameworks/Maple.framework/Versions/Current/bin/maple $@
          (b) Save the file.
          (c) Make the file executable.
                chmod +x maple

      * WINDOWS:
        You must install Maple-for-Linux into the VMware machine (sorry, that's
        the only way at present).
"""

    def expect(self):
        """
        Return the pexpect object for this Maple session.

        EXAMPLES::

            sage: m = Maple()          # optional - maple
            sage: m.expect() is None   # optional - maple
            True
            sage: m._start()           # optional - maple
            sage: m.expect()           # optional - maple
            <pexpect.spawn instance at 0x...>
            sage: m.quit()             # optional - maple
        """
        return self._expect

    def console(self):
        r"""
        Spawn a new Maple command-line session.

        EXAMPLES::

            sage: maple.console() # not tested
                |^/|     Maple 11 (IBM INTEL LINUX)
            ._|\|   |/|_. Copyright (c) Maplesoft, a division of Waterloo Maple Inc. 2007
             \  MAPLE  /  All rights reserved. Maple is a trademark of
             <____ ____>  Waterloo Maple Inc.
                  |       Type ? for help.
            >
        """
        maple_console()

#     def killall(self):
#         """
#         Kill all running instances of the maple interpreter
#         on this system.

#         TODO: When Sage exits it doesn't correctly by default kill
#         all running Maple interpreters, for some strange reason.
#         Calling this function uses the kill and pidof operating system
#         programs to find all instances of cmaple and kill them.
#         """
#         import os
#         self._expect = None
#         while True:
#             pid = os.popen("pidof cmaple").read()[:-1]
#             if len(pid) > 0:
#                 os.system('kill -9 %s'%pid)
#             else:
#                 break

    def completions(self, s):
        """
        Return all commands that complete the command starting with the
        string s. This is like typing s[Ctrl-T] in the maple interpreter.

        EXAMPLES::

            sage: c = maple.completions('di')  # optional - maple
            sage: 'divide' in c                # optional - maple
            True
        """
        bs = chr(8) * len(s)
        if self._expect is None:
            self._start()
        E = self._expect
        E.sendline('%s%s%s' % (s, chr(20), bs))
        t = E.timeout
        E.timeout = 0.3  # since some things have no completion
        try:
            E.expect('----')
        except pexpect.TIMEOUT:
            E.timeout = t
            return []
        E.timeout = t
        v = E.before
        E.expect(self._prompt)
        E.expect(self._prompt)
        return v.split()[2:]

    def _commands(self):
        """
        Return list of all commands defined in Maple.

        EXAMPLES::

            sage: c = maple._commands() # optional - maple
            sage: len(c) > 100          # optional - maple
            True
            sage: 'dilog' in c          # optional - maple
            True
        """
        try:
            v = sum([self.completions(chr(65 + n)) for n in range(26)], []) + \
                sum([self.completions(chr(97 + n)) for n in range(26)], [])
        except RuntimeError:
            print("\n" * 3)
            print("*" * 70)
            print("WARNING: You do not have a working version of Maple installed!")
            print("*" * 70)
            v = []
        v.sort()
        return v

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        """
        Return a list of all the commands defined in Maple and optionally
        (per default) store them to disk.

        EXAMPLES::

            sage: c = maple._tab_completion(use_disk_cache=False, verbose=False) # optional - maple
            sage: len(c) > 100  # optional - maple
            True
            sage: 'dilog' in c  # optional - maple
            True
        """
        try:
            return self.__tab_completion
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__tab_completion = sage.misc.persist.load(COMMANDS_CACHE)
                    return self.__tab_completion
                except IOError:
                    pass
            if verbose:
                print("\nBuilding Maple command completion list (this takes")
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % COMMANDS_CACHE)
            v = self._commands()
            self.__tab_completion = v
            if len(v) > 200:
                # Maple is actually installed.
                sage.misc.persist.save(v, COMMANDS_CACHE)
            return v

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: maple._eval_line('2+2')  # optional - maple
            '4'
        """
        line += ';'
        with gc_disabled():
            z = Expect._eval_line(self, line, allow_use_file=allow_use_file,
                    wait_for_prompt=wait_for_prompt).replace('\\\n', '').strip()
            if z.lower().find("error") != -1:
                raise RuntimeError("An error occurred running a Maple command:\nINPUT:\n%s\nOUTPUT:\n%s" % (line, z))
        return z

    def _eval_line_using_file(self, line, *args, **kwargs):
        """
        EXAMPLES::

            sage: maple._eval_line_using_file('2+2')  # optional - maple
            '4'
        """
        line += ';'  # Adds the maple ";" thing like in self._eval_line
        return Expect._eval_line_using_file(self, line, *args, **kwargs)

    def cputime(self, t=None):
        r"""
        Return the amount of CPU time that the Maple session has used. If
        ``t`` is not None, then it returns the difference
        between the current CPU time and ``t``.

        EXAMPLES::

            sage: t = maple.cputime() # optional - maple
            sage: t                   # random; optional - maple
            0.02
            sage: x = maple('x')      # optional - maple
            sage: maple.diff(x^2, x)  # optional - maple
            2*x
            sage: maple.cputime(t)    # random; optional - maple
            0.0
        """
        if t is None:
            return float(self('time()'))
        else:
            return float(self('time() - %s' % float(t)))

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: maple.set('xx', '2') # optional - maple
            sage: maple.get('xx')      # optional - maple
            '2'
        """
        cmd = '%s:=%s:' % (var, value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError("Error executing code in Maple\nCODE:\n\t%s\nMaple ERROR:\n\t%s" % (cmd, out))

    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: maple.set('xx', '2') # optional - maple
            sage: maple.get('xx')      # optional - maple
            '2'
        """
        return self.eval('printf("%%q",%s)' % var)

    def _object_class(self):
        """
        Return the class of MapleElements.

        EXAMPLES::

            sage: maple._object_class()
            <class 'sage.interfaces.maple.MapleElement'>

        ::

            sage: m = maple(2)  # optional - maple
            sage: type(m)       # optional - maple
            <class 'sage.interfaces.maple.MapleElement'>
        """
        return MapleElement

    def _function_element_class(self):
        """
        Return the MapleFunctionElement class.

        EXAMPLES::

            sage: maple._function_element_class()
            <class 'sage.interfaces.maple.MapleFunctionElement'>

        ::

            sage: two = maple(2)  # optional - maple
            sage: type(two.gcd)   # optional - maple
            <class 'sage.interfaces.maple.MapleFunctionElement'>
        """
        return MapleFunctionElement

    def _equality_symbol(self):
        """
        Return the symbol used for equality testing in Maple.

        EXAMPLES::

            sage: maple._equality_symbol()
            '='

            sage: maple(2) == maple(2) # optional - maple
            True
        """
        return '='

    def _true_symbol(self):
        """
        Return the symbol used for truth in Maple.

        EXAMPLES::

            sage: maple._true_symbol()
            'true'

        ::

            sage: maple(2) == maple(2) # optional - maple
            True
        """
        return 'true'

    def _assign_symbol(self):
        """
        Return the symbol used for assignment in Maple.

        EXAMPLES::

            sage: maple._assign_symbol()
            ':='
        """
        return ":="

    def _source(self, s):
        """
        Tries to return the source code of a Maple function str as a
        string.

        EXAMPLES::

            sage: print(maple._source('curry').strip())  # optional - maple
            p -> subs('_X' = args[2 .. nargs], () -> p(_X, args))
            sage: maple._source('ZZZ')                  #not tested
            Traceback (most recent call last):
            ...
            Exception: no source code could be found
        """
        cmd = 'echo "interface(verboseproc=2): print(%s);" | maple -q' % s
        src = os.popen(cmd).read()
        if src.strip() == s:
            raise RuntimeError("no source code could be found")
        else:
            return src

    def source(self, s):
        """
        Display the Maple source (if possible) about s. This is the same as
        returning the output produced by the following Maple commands:

        interface(verboseproc=2): print(s)

        INPUT:


        -  ``s`` - a string representing the function whose
           source code you want


        EXAMPLES::

            sage: maple.source('curry')  #not tested
            p -> subs('_X' = args[2 .. nargs], () -> p(_X, args))
        """
        try:
            pager()(self._source(s))
        except Exception:
            pager()('No source code could be found.')

    def _help(self, string):
        r"""
        Return the Maple help on ``string``.

        EXAMPLES::

            sage: txt = maple._help('gcd')  # optional - maple
            sage: txt.find('gcd - greatest common divisor') > 0 # optional - maple
            True
        """
        return os.popen('echo "?%s" | maple -q' % string).read()

    def help(self, string):
        """
        Display Maple help about ``string``.

        This is the same as typing "?string" in the Maple console.

        INPUT:

        -  ``string`` - a string to search for in the maple help
           system

        EXAMPLES::

            sage: maple.help('Psi')  # not tested
            Psi - the Digamma and Polygamma functions
            ...
        """
        pager()(self._help(string))

    def with_package(self, package):
        """
        Make a package of Maple procedures available in the interpreter.

        INPUT:


        -  ``package`` - string


        EXAMPLES: Some functions are unknown to Maple until you use with to
        include the appropriate package.

        ::

            sage: maple.quit()   # reset maple; optional -- maple
            sage: maple('partition(10)')              # optional - maple
            partition(10)
            sage: maple('bell(10)')                   # optional - maple
            bell(10)
            sage: maple.with_package('combinat')      # optional - maple
            sage: maple('partition(10)')              # optional - maple
            [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 2], [1, 1, 1, 1, 1, 1, 2, 2], [1, 1, 1, 1, 2, 2, 2], [1, 1, 2, 2, 2, 2], [2, 2, 2, 2, 2], [1, 1, 1, 1, 1, 1, 1, 3], [1, 1, 1, 1, 1, 2, 3], [1, 1, 1, 2, 2, 3], [1, 2, 2, 2, 3], [1, 1, 1, 1, 3, 3], [1, 1, 2, 3, 3], [2, 2, 3, 3], [1, 3, 3, 3], [1, 1, 1, 1, 1, 1, 4], [1, 1, 1, 1, 2, 4], [1, 1, 2, 2, 4], [2, 2, 2, 4], [1, 1, 1, 3, 4], [1, 2, 3, 4], [3, 3, 4], [1, 1, 4, 4], [2, 4, 4], [1, 1, 1, 1, 1, 5], [1, 1, 1, 2, 5], [1, 2, 2, 5], [1, 1, 3, 5], [2, 3, 5], [1, 4, 5], [5, 5], [1, 1, 1, 1, 6], [1, 1, 2, 6], [2, 2, 6], [1, 3, 6], [4, 6], [1, 1, 1, 7], [1, 2, 7], [3, 7], [1, 1, 8], [2, 8], [1, 9], [10]]
            sage: maple('bell(10)')                   # optional - maple
            115975
            sage: maple('fibonacci(10)')              # optional - maple
            55
        """
        self.eval('with(%s)' % package)

    load = with_package

    def clear(self, var):
        """
        Clear the variable named ``var``.

        To clear a Maple variable, you must assign 'itself' to itself.
        In Maple 'expr' prevents expr to be evaluated.

        EXAMPLES::

            sage: maple.set('xx', '2')  # optional - maple
            sage: maple.get('xx')       # optional - maple
            '2'
            sage: maple.clear('xx')     # optional - maple
            sage: maple.get('xx')       # optional - maple
            'xx'
        """
        self.set(var, "'{}'".format(var))


@instancedoc
class MapleFunction(ExpectFunction):
    def _instancedoc(self):
        """
        Return the Maple help for this function. This gets called when
        doing "?" on self.

        EXAMPLES::

            sage: txt = maple.gcd.__doc__  # optional - maple
            sage: txt.find('gcd - greatest common divisor') > 0 # optional - maple
            True
        """
        M = self._parent
        return M._help(self._name)

    def _sage_src_(self):
        """
        Return the source code of ``self``.

        This is the function that eventually gets called when doing
        maple.gcd?? for example.

        EXAMPLES::

            sage: print(maple.curry._sage_src_().strip()) # optional - maple
            p -> subs('_X' = args[2 .. nargs], () -> p(_X, args))
            sage: maple.ZZZ._sage_src_()                 #not tested
            Traceback (most recent call last):
            ...
            Exception: no source code could be found
        """
        M = self._parent
        return M._source(self._name)


@instancedoc
class MapleFunctionElement(FunctionElement):
    def _instancedoc_(self):
        """
        Return the Maple help for this function.

        This gets called when doing "?" on ``self``.

        EXAMPLES::

            sage: two = maple(2)  # optional - maple
            sage: txt = two.gcd.__doc__  # optional - maple
            sage: txt.find('gcd - greatest common divisor') > 0 # optional - maple
            True
        """
        return self._obj.parent()._help(self._name)

    def _sage_src_(self):
        """
        Return the source code of self.

        EXAMPLES::

            sage: g = maple('gcd')                   # optional - maple
            sage: print(g.curry._sage_src_().strip()) # optional - maple
            p -> subs('_X' = args[2 .. nargs], () -> p(_X, args))
            sage: m = maple('2')                     # optional - maple
            sage: m.ZZZ._sage_src_()                 #not tested
            Traceback (most recent call last):
            ...
            Exception: no source code could be found
        """
        return self._obj.parent()._source(self._name)


@instancedoc
class MapleElement(ExtraTabCompletion, ExpectElement):

    def __float__(self):
        """
        Return a floating point version of self.

        EXAMPLES::

            sage: float(maple(1/2))  # optional - maple
            0.5
            sage: type(_)            # optional - maple
            <... 'float'>
        """
        return float(maple.eval('evalf(%s)' % self.name()))

    def __hash__(self):
        """
        Return a 64-bit integer representing the hash of self. Since
        Python uses 32-bit hashes, it will automatically convert the result
        of this to a 32-bit hash.

        These examples are optional, and require Maple to be installed. You
        don't need to install any Sage packages for this.

        EXAMPLES::

            sage: m = maple('x^2+y^2')                      # optional - maple
            sage: m.__hash__()                              # optional - maple
            188724254834261060184983038723355865733
            sage: hash(m)               # random            # optional - maple
            5035731711831192733
            sage: m = maple('x^2+y^3')                      # optional - maple
            sage: m.__hash__()          # random            # optional - maple
            264835029579301191531663246434344770556
            sage: hash(m)               # random            # optional - maple
            -2187277978252104690
        """
        return int(maple.eval('StringTools:-Hash(convert(%s, string))' % self.name())[1:-1], 16)

    def _richcmp_(self, other, op):
        """
        Compare equality between self and other, using maple.

        These examples are optional, and require Maple to be installed. You
        don't need to install any Sage packages for this.

        EXAMPLES::

            sage: a = maple(5)                             # optional - maple
            sage: b = maple(5)                             # optional - maple
            sage: a == b                                   # optional - maple
            True
            sage: a == 5                                   # optional - maple
            True

        ::

            sage: c = maple(3)                             # optional - maple
            sage: a == c                                   # optional - maple
            False
            sage: a < c                                    # optional - maple
            False
            sage: a < 6                                    # optional - maple
            True
            sage: c <= a                                   # optional - maple
            True

        ::

            sage: M = matrix(ZZ, 2, range(1,5))            # optional - maple
            sage: Mm = maple(M)                            # optional - maple
            sage: Mm == Mm                                 # optional - maple
            True

        TESTS::

            sage: x = var('x')
            sage: t = maple((x+1)^2)                       # optional - maple
            sage: u = maple(x^2+2*x+1)                     # optional - maple
            sage: u == t # todo: not implemented
            True         # returns False, should use 'testeq' in maple
            sage: maple.eval('testeq(%s = %s)' % (t.name(),u.name()))    # optional - maple
            'true'
        """
        P = self.parent()
        if P.eval("evalb(%s %s %s)" % (self.name(), P._equality_symbol(),
                                       other.name())) == P._true_symbol():
            return rich_to_bool(op, 0)
        # Maple does not allow comparing objects of different types and
        # it raises an error in this case.
        # We catch the error, and return True for <
        try:
            if P.eval("evalb(%s %s %s)" % (self.name(), P._lessthan_symbol(),
                                           other.name())) == P._true_symbol():
                return rich_to_bool(op, -1)
        except RuntimeError as e:
            msg = str(e)
            if 'is not valid' in msg and 'to < or <=' in msg:
                if (hash(str(self)) < hash(str(other))):
                    return rich_to_bool(op, -1)
                else:
                    return rich_to_bool(op, 1)
            else:
                raise RuntimeError(e)
        if P.eval("evalb(%s %s %s)" % (self.name(), P._greaterthan_symbol(),
                                       other.name())) == P._true_symbol():
            return rich_to_bool(op, 1)
        return NotImplemented

    def _mul_(self, right):
        """
        These examples are optional, and require Maple to be installed. You
        don't need to install any Sage packages for this.

        EXAMPLES::

            sage: t = maple(5); u = maple(3)                # optional - maple
            sage: t*u                                       # optional - maple
            15
            sage: t._mul_(u)                                # optional - maple
            15
            sage: M = matrix(ZZ,2,range(4))                 # optional - maple
            sage: Mm = maple(M)                             # optional - maple
            sage: Mm*Mm                                     # optional - maple
            Matrix(2, 2, [[2,3],[6,11]])

        ::

            sage: v = vector(ZZ,2,[2,3])                    # optional - maple
            sage: vm = maple(v)                             # optional - maple
            sage: vm*Mm                                     # optional - maple
            Vector[row](2, [6,11])

        ::

            sage: t*Mm                                      # optional - maple
            Matrix(2, 2, [[0,5],[10,15]])
        """
        P = self._check_valid()
        try:
            return P.new('%s . %s' % (self._name, right._name))
        except Exception as msg:
            raise TypeError(msg)

    def _tab_completion(self):
        """
        EXAMPLES::

            sage: a = maple(2) # optional - maple
            sage: 'sin' in a._tab_completion() # optional - maple
            True
        """
        return self.parent()._tab_completion()

    def _latex_(self):
        r"""
        You can output Maple expressions in latex.

        EXAMPLES::

            sage: print(latex(maple('(x^4 - y)/(y^2-3*x)')))   # optional - maple
            {\frac {{x}^{4}-y}{{y}^{2}-3\,x}}
            sage: print(latex(maple(pi - e^3)))                 # optional - maple
            \pi-{{\rm e}^{3}}
            sage: print(maple(pi - e^3)._latex_())              # optional - maple
            \pi-{{\rm e}^{3}}

        .. note::

           Some expressions might require the Maple style file
           ``maple2e.sty`` in order to latex correctly.
        """
        return self.parent().eval('latex(%s)' % self.name())

    def op(self, i):
        """
        Return the i-th operand of this expression.

        INPUT:

        - i -- an integer

        EXAMPLES::

            sage: V = maple(vector(QQ,[4,5,6]))      # optional - maple
            sage: V.op(1)                            # optional - maple
            3
            sage: V.op(2)                            # optional - maple
            {1 = 4, 2 = 5, 3 = 6}
        """
        return self.parent().op(i, self)

    def _sage_(self):
        r"""
        Convert a maple expression back to a Sage expression.

        This currently does not implement a parser for the Maple output language,
        therefore only very simple expressions will convert successfully.

        REFERENCE:

        https://www.asc.tuwien.ac.at/compmath/download/Monagan_Maple_Programming.pdf

        EXAMPLES::

            sage: m = maple('x^2 + 5*y')        # optional - maple
            sage: m.sage()                      # optional - maple
            x^2 + 5*y
            sage: m._sage_()                    # optional - maple
            x^2 + 5*y

            sage: m = maple('sin(sqrt(1-x^2)) * (1 - cos(1/x))^2')  # optional - maple
            sage: m.sage()                                          # optional - maple
            (cos(1/x) - 1)^2*sin(sqrt(-x^2 + 1))

        Some matrices can be converted back::

            sage: m = matrix(2, 2, [1, 2, x, 3])    # optional - maple
            sage: mm = maple(m)                     # optional - maple
            sage: mm.sage() == m                    # optional - maple
            True

        Some vectors can be converted back::

            sage: m = vector([1, x, 2, 3])          # optional - maple
            sage: mm = maple(m)                     # optional - maple
            sage: mm.sage() == m                    # optional - maple
            True

        Integers and rationals are converted as such::

            sage: maple(33).sage().parent()         # optional - maple
            Integer Ring
            sage: maple(191/5).sage().parent()      # optional - maple
            Rational Field

        Sets, lists, sequences::

            sage: maple("[4,5,6]").sage()           # optional - maple
            [4, 5, 6]
            sage: maple({14,33,6}).sage()           # optional - maple
            {6, 14, 33}
            sage: maple("seq(i**2,i=1..5)").sage()  # optional - maple
            (1, 4, 9, 16, 25)

        Strings::

            sage: maple('"banane"').sage()          # optional - maple
            '"banane"'

        Floats::

            sage: Z3 = maple('evalf(Zeta(3))')   # optional - maple
            sage: Z3.sage().parent()             # optional - maple
            Real Field with 53 bits of precision

            sage: sq5 = maple('evalf(sqrt(5),100)')   # optional - maple
            sage: sq5 = sq5.sage(); sq5               # optional - maple
            2.23606797749978969640...
            sage: sq5.parent()                        # optional - maple
            Real Field with 332 bits of precision

        Functions are not yet converted back correctly::

            sage: maple(hypergeometric([3,4],[5],x))  # optional - maple
            hypergeom([3, 4],[5],x)
            sage: _.sage()                # known bug # optional - maple
            hypergeometric((3, 4), (5,), x)
        """
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector
        from sage.rings.integer_ring import ZZ
        # The next few lines are a very crude excuse for a maple "parser"
        maple_type = repr(self.whattype())
        result = repr(self)
        result = result.replace("Pi", "pi")
        if maple_type == 'symbol':       # pi
            pass   # left to symbolic ring
        elif maple_type == 'string':     # "banane"
            return result
        elif maple_type == 'exprseq':   # 2, 2
            n = self.parent()(f"[{self._name}]").nops()._sage_()
            return tuple(self[i] for i in range(1, n + 1))
        elif maple_type == 'set':        # {1, 2}
            n = self.nops()._sage_()
            return set(self.op(i)._sage_() for i in range(1, n + 1))
        elif maple_type == 'list':        # [1, 2]
            n = self.nops()._sage_()
            return [self.op(i)._sage_() for i in range(1, n + 1)]
        elif maple_type == "Matrix":     # Matrix(2, 2, [[1,2],[3,4]])
            mn = self.op(1)
            m = mn[1]._sage_()
            n = mn[2]._sage_()
            coeffs = [self[i + 1, j + 1]._sage_()
                      for i in range(m) for j in range(n)]
            return matrix(m, n, coeffs)
        elif maple_type[:6] == "Vector":  # Vector[row](3, [4,5,6])
            n = self.op(1)._sage_()
            return vector([self[i + 1]._sage_() for i in range(n)])
        elif maple_type == 'integer':
            return ZZ(result)
        elif maple_type == 'fraction':
            return self.op(1)._sage_() / self.op(2)._sage_()
        elif maple_type == "function":
            pass  # TODO : here one should translate back function names
        elif maple_type == "float":
            from sage.rings.real_mpfr import RealField
            mantissa = len(repr(self.op(1)))
            prec = max(53, (mantissa * 13301) // 4004)
            R = RealField(prec)
            return R(result)
        elif maple_type == '`=`':        # (1, 1) = 2
            return (self.op(1)._sage_() == self.op(2)._sage())
        try:
            from sage.symbolic.ring import SR
            return SR(result)
        except Exception:
            raise NotImplementedError("Unable to parse Maple output: %s" % result)


maple = Maple()  # an instance


def reduce_load_Maple():
    """
    Return the maple object created in sage.interfaces.maple.

    EXAMPLES::

        sage: from sage.interfaces.maple import reduce_load_Maple
        sage: reduce_load_Maple()
        Maple
    """
    return maple


def maple_console():
    r"""
    Spawn a new Maple command-line session.

    EXAMPLES::

        sage: maple_console() #not tested
            |^/|     Maple 11 (IBM INTEL LINUX)
        ._|\|   |/|_. Copyright (c) Maplesoft, a division of Waterloo Maple Inc. 2007
         \  MAPLE  /  All rights reserved. Maple is a trademark of
         <____ ____>  Waterloo Maple Inc.
              |       Type ? for help.
        >
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%maple magics instead.')
    os.system('maple')


def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.maple import __doctest_cleanup
        sage: m = maple(2)         # optional - maple
        sage: maple.is_running()   # optional - maple
        True
        sage: __doctest_cleanup()
        sage: maple.is_running()
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()


"""
The following only works in Maple >= 9, I guess, but could
be useful.

From Jaap Spies:  In addition Maple has a nice feature the function

 > FunctionAdvisor();

 > FunctionAdvisor(topics, quiet);
      [DE, analytic_extension, asymptotic_expansion, branch_cuts,
      branch_points, calling_sequence, class_members,
      classify_function, definition, describe, differentiation_rule,
      function_classes, identities, integral_form,
      known_functions, relate, series, singularities, special_values,
      specialize, sum_form, synonyms]

 > FunctionAdvisor(syntax, hypergeom);
                                            hypergeom([a, b], [c], z)

Eventually this could be used to do an intelligent command
completion.
"""
