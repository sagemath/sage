r"""
Interface to MuPAD

AUTHOR:

- Mike Hansen
- William Stein

You must have the optional commercial MuPAD interpreter installed and
available as the command \code{mupkern} in your PATH in order to use
this interface.  You do not have to install any optional \sage
packages.

TESTS::

    sage: mupad.package('"MuPAD-Combinat"')                        # optional - mupad
    sage: combinat = mupad.combinat                                # optional - mupad
    sage: examples = mupad.examples                                # optional - mupad
    sage: S = examples.SymmetricFunctions()                        # optional - mupad
    sage: S.s[2,1]^2                                               # optional - mupad
    s[3, 3] + s[4, 2] + s[2, 2, 1, 1] + s[2, 2, 2] + 2 s[3, 2, 1] + s[4, 1, 1] +
    s[3, 1, 1, 1]
    sage: S.omega( S.s[3] )                                        # optional - mupad
    s[1, 1, 1]
    sage: s = S.s                                                  # optional - mupad
    sage: p = S.p                                                  # optional - mupad
    sage: s(s[2,1] + p[2,1])                                       # optional - mupad
    s[2, 1] + s[3] - s[1, 1, 1]
    sage: s(_)                                                     # optional - mupad
    s[2, 1] + s[3] - s[1, 1, 1]

    sage: combinat.tableaux.list(3)                                # optional - mupad # note: the order of the result seems to depend on the version of MuPAD / MuPAD-Combinat
                --                                      +---+ --
                |                                       | 3 |  |
                |                 +---+      +---+      +---+  |
                |                 | 3 |      | 2 |      | 2 |  |
                |  +---+---+---+  +---+---+  +---+---+  +---+  |
                |  | 1 | 2 | 3 |, | 1 | 2 |, | 1 | 3 |, | 1 |  |
                -- +---+---+---+  +---+---+  +---+---+  +---+ --
    sage: three = mupad(3)                                        # optional - mupad
    sage: three.combinat.tableaux.list()                          # optional - mupad
                --                                      +---+ --
                |                                       | 3 |  |
                |                 +---+      +---+      +---+  |
                |                 | 3 |      | 2 |      | 2 |  |
                |  +---+---+---+  +---+---+  +---+---+  +---+  |
                |  | 1 | 2 | 3 |, | 1 | 2 |, | 1 | 3 |, | 1 |  |
                -- +---+---+---+  +---+---+  +---+---+  +---+ --
    sage: t = _[1]                                                # optional - mupad
    sage: t                                                       # optional - mupad
                                 +---+---+---+
                                 | 1 | 2 | 3 |
                                 +---+---+---+
    sage: combinat.tableaux.conjugate(t)                          # optional - mupad
                                     +---+
                                     | 3 |
                                     +---+
                                     | 2 |
                                     +---+
                                     | 1 |
                                     +---+

    sage: combinat.ribbonsTableaux.list([2,2],[1,1],2)           # optional - mupad
                           -- +---+---+  +---+---+ --
                           |  |   | 2 |  |     2 |  |
                           |  +   +   +, +---+---+  |
                           |  | 1 |   |  | 1     |  |
                           -- +---+---+  +---+---+ --
    sage: combinat.tableaux.kAtom([2,1],3)                      # optional - mupad
                                  -- +---+     --
                                  |  | 2 |      |
                                  |  +---+---+  |
                                  |  | 1 | 1 |  |
                                  -- +---+---+ --
    sage: M = S.Macdonald()                                    # optional - mupad
    sage: a = M.P[1]^2                                         # optional - mupad
    sage: mupad.mapcoeffs(a, 'normal')                         # optional - mupad
                                 q - t + q t - 1
                          P[2] + --------------- P[1, 1]
                                     q t - 1

"""

#############################################################################
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#############################################################################

import os

from .expect import (Expect, ExpectElement, ExpectFunction,
                    FunctionElement)
from sage.interfaces.interface import AsciiArtString
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.env import DOT_SAGE
from sage.docs.instancedoc import instancedoc

COMMANDS_CACHE = '%s/mupad_commandlist_cache.sobj' % DOT_SAGE
PROMPT = ">>"
seq = 0


class Mupad(ExtraTabCompletion, Expect):
    """
    Interface to the MuPAD interpreter.
    """
    def __init__(self, maxread=None, script_subdirectory=None, server=None, server_tmpdir=None, logfile=None):
        """
        Create an instance of the MuPAD interpreter.

        EXAMPLES::

            sage: mupad == loads(dumps(mupad))                      # optional - mupad
            True
        """
        Expect.__init__(self,
                        name = 'MuPAD',
                        prompt = PROMPT,
                        # the -U SAGE=TRUE allows for MuPAD programs to test whether they are run from Sage
                        command = "mupkern -P e -U SAGE=TRUE",
                        script_subdirectory = script_subdirectory,
                        server = server,
                        server_tmpdir = server_tmpdir,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        logfile = None)




    def _function_class(self):
        """
        EXAMPLES::

            sage: mupad._function_class()
            <class 'sage.interfaces.mupad.MupadFunction'>

            sage: mdiff = mupad.diff; mdiff  # optional - mupad
            diff
            sage: type(mdiff)                # optional -- mupad
            <class 'sage.interfaces.mupad.MupadFunction'>
        """
        return MupadFunction

    def __reduce__(self):
        """
        EXAMPLES::

            sage: mupad.__reduce__()
            (<function reduce_load_mupad at 0x...>, ())

        """
        return reduce_load_mupad, tuple([])

    def _read_in_file_command(self, filename):
        """
        EXAMPLES::

            sage: mupad._read_in_file_command('test')
            'read("test")'

            sage: filename = tmp_filename()
            sage: f = open(filename, 'w')
            sage: _ = f.write('x := 2;\n')
            sage: f.close()
            sage: mupad.read(filename)   # optional - MuPAD
            sage: mupad.get('x').strip() # optional - mupad
            '2'

        """
        return 'read("%s")'%filename

    def _quit_string(self):
        """
        EXAMPLES::

            sage: mupad._quit_string()
            'quit'
        """
        return 'quit'

    def _install_hints(self):
        """
        Hints for installing MuPAD on your computer.

        EXAMPLES::

            sage: print(mupad._install_hints())
            <BLANKLINE>
            In order to use the MuPAD interface you need to have MuPAD installed
            ...

        """
        return """
In order to use the MuPAD interface you need to have MuPAD installed
and have a script in your PATH called "mupkern" that runs the
command-line version of MuPAD.

  (1) You might have to buy MuPAD.

  (2) * LINUX: The mupkern script comes standard with your Mupad install.

      * APPLE OS X:
         ???
"""

    def expect(self):
        """
        EXAMPLES::

            sage: a = mupad(1)   # optional - mupad
            sage: mupad.expect() # optional - mupad
            <pexpect.spawn instance at 0x...>
        """
        return self._expect

    def console(self):
        """
        Spawn a new MuPAD command-line session.

        EXAMPLES::

            sage: mupad.console() #not tested

               *----*    MuPAD Pro 4.0.2 -- The Open Computer Algebra System
              /|   /|
             *----* |    Copyright (c)  1997 - 2007  by SciFace Software
             | *--|-*                   All rights reserved.
             |/   |/
             *----*      Licensed to:   ...

        """
        mupad_console()

    def eval(self, code, strip=True, **kwds):
        """
        EXAMPLES::

            sage: mupad.eval('2+2')   # optional - mupad
                                                   4

        """
        s = Expect.eval(self, code, **kwds)
        return AsciiArtString(s)

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True,
                   need_output=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: mupad._eval_line('2+2')  # optional - mupad
            '                                       4'
            sage: mupad._eval_line('x::asdf') # optional - mupad
            Traceback (most recent call last):
            ...
            RuntimeError: Unknown slot "x::asdf" [slot]

        """
        if self._expect is None:
            self._start()
        if not need_output:
            E = self._expect
            E.sendline(line)
            return

        global seq
        seq += 1
        START = '__start__(%s+1)'%seq
        END = '__end__(%s+1)'%seq
        line = '%s; %s; %s;'%(START, line, END)
        START = '__start__(%s)'%(seq+1)
        END = '__end__(%s)'%(seq+1)

        E = self._expect
        E.sendline(line)
        E.expect(PROMPT)
        z = E.before
        i = z.find(START)
        if i == -1:
            raise RuntimeError("%s\nError evaluating code in MuPAD"%z)
        z = z[i+len(START)+2:]
        z = z.rstrip().rstrip(END).rstrip('"').rstrip().strip('\n').strip('\r').strip('\n').replace('\\\r\n','')
        i = z.find('Error: ')
        if i != -1:
            raise RuntimeError(z[i + 7:])
        return z

    def cputime(self, t=None):
        """
        EXAMPLES::

            sage: t = mupad.cputime() #random, optional - MuPAD
            0.11600000000000001
        """
        if t is None:
            return float(str(self('time()')))/1000
        else:
            return float(str(self('time() - %s'%float(t))))/1000

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: mupad.set('a', 4) # optional - mupad
            sage: mupad.get('a').strip() # optional - mupad
            '4'
        """
        cmd = '%s:=%s:'%(var,value)
        out = self.eval(cmd)
        i = out.find('Error: ')
        if i != -1:
            raise RuntimeError(out[i + 7:])

    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: mupad.set('a', 4) # optional - mupad
            sage: mupad.get('a').strip() # optional - mupad
            '4'

        """
        s = self.eval('%s'%var)
        i = s.find('=')
        return s[i+1:]

    def _object_class(self):
        """
        EXAMPLES::

            sage: mupad._object_class()
            <class 'sage.interfaces.mupad.MupadElement'>
        """
        return MupadElement

    def _equality_symbol(self):
        """
        EXAMPLES::

            sage: mupad._equality_symbol()
            '='
        """
        return '='

    def _assign_symbol(self):
        """
        EXAMPLES::

            sage: mupad._assign_symbol()
            ':='
        """
        return ":="

    def _continuation_prompt(self):
        """
        EXAMPLES::

            sage: mupad._continuation_prompt()
            '&>'
        """
        return "&>"

    def _commands(self):
        """
        Return list of all commands defined in MuPAD.

        EXAMPLES::

            sage: cmds = mupad._commands()  # optional - mupad
            sage: len(cmds) > 100 # optional - mupad
            True
            sage: 'diff' in cmds  # optional - mupad
            True
        """
        try:
            v = sum([self.completions(chr(65+n)) for n in range(26)], []) + \
                sum([self.completions(chr(97+n)) for n in range(26)], [])
        except RuntimeError:
            print("\n" * 3)
            print("*" * 70)
            print("WARNING: You do not have a working version of MuPAD installed!")
            print("*" * 70)
            v = []
        v.sort()
        return v

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        """
        EXAMPLES::

            sage: names = mupad._tab_completion() # optional - mupad
            sage: len(names) > 100 # optional - mupad
            True
            sage: 'combinat' in names # optional - mupad
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
                print("\nBuilding MuPAD command completion list (this takes")
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % COMMANDS_CACHE)
            v = self._commands()
            self.__tab_completion = v
            if len(v) > 200:
                # MuPAD is actually installed.
                sage.misc.persist.save(v, COMMANDS_CACHE)
            return v

    def completions(self, string, strip=False):
        """
        EXAMPLES::

            sage: mupad.completions('linal') # optional - mupad
            ['linalg']
        """
        res = self.eval('_pref(Complete)("%s")'%string).strip()
        res = res.replace('\n', '').split(',')
        res = [s.strip().strip('"') for s in res]
        res = [s for s in res if not s.endswith('::')]
        if strip:
            n = len(string)
            res = [s[n:] for s in res]

        return res if res != [''] else []


@instancedoc
class MupadFunction(ExtraTabCompletion, ExpectFunction):
    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: mupad.diff.__doc__
            No help on diff available
        """
        M = self._parent
        return M.help(self._name)

    def __getattr__(self, attrname):
        """
        EXAMPLES::

            sage: mupad.linalg.addRow
            linalg::addRow
        """
        if attrname[:1] == "_":
            raise AttributeError
        return MupadFunction(self._parent, self._name+"::"+attrname)

    def _tab_completion(self):
        """
        EXAMPLES::

            sage: mupad.linalg._tab_completion() # optional - mupad
            ['addCol',
             'addRow',
             ...
             'wiedemann']

        """
        res = self._parent.completions(self._name+"::", strip=True)
        return res if res != [] else self._parent._tab_completion()


@instancedoc
class MupadFunctionElement(ExtraTabCompletion, FunctionElement):
    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: x = mupad('x')  # optional - mupad
            sage: x.diff.__doc__  # optional - mupad
            No help on diff available

        """
        return self._obj.parent().help(self._name)

    def __getattr__(self, attrname):
        """
        EXAMPLES::

            sage: mupad.package('"MuPAD-Combinat"')  # optional - mupad-Combinat
            sage: combinat = mupad.combinat          # optional - mupad-Combinat
            sage: three = mupad(3)                   # optional - mupad-Combinat
            sage: type(three.combinat)               # optional - mupad-Combinat
            <class 'sage.interfaces.mupad.MupadFunctionElement'>
            sage: tableaux = three.combinat.tableaux # optional - mupad-Combinat
            sage: type(tableaux)                     # optional - mupad-Combinat
            <class 'sage.interfaces.mupad.MupadFunctionElement'>
        """
        P = self._obj.parent()
        if attrname[:1] == "_":
            if attrname not in self.__dict__:
                raise AttributeError
            else:
                return self.__dict__[attrname]
        name = self._name+"::"+attrname
        if P.eval('type(%s)'%name) == "DOM_DOMAIN":
            return MupadElement(P, name)
        else:
            return MupadFunctionElement(self._obj, name)

    def _tab_completion(self):
        """
        EXAMPLES::

            sage: three = mupad(3) # optional - mupad
            sage: 'list' in three.combinat.tableaux._tab_completion() # optional - mupad
            True
        """
        P = self._obj.parent()
        res = P.completions(self._name+"::", strip=True)
        return res if res != [] else P._tab_completion()


    def __call__(self, *args):
        """
        EXAMPLES::

            sage: mupad.package('"MuPAD-Combinat"') # optional - mupad-Combinat
            sage: combinat = mupad.combinat         # optional - mupad-Combinat
            sage: examples = mupad.examples         # optional - mupad-Combinat
            sage: S = examples.SymmetricFunctions() # optional - mupad-Combinat
            sage: type(S.omega)                     # optional - mupad-Combinat
            <class 'sage.interfaces.mupad.MupadFunctionElement'>
            sage: S.omega(S.s[3])                   # optional - mupad-Combinat
            s[1, 1, 1]
        """
        P = self._obj.parent()
        if P.eval('type(%s)'%(self._obj.name())).strip() == "DOM_DOMAIN":
            return P.function_call(self._name, list(args))
        else:
            return P.function_call(self._name, [self._obj] + list(args))


@instancedoc
class MupadElement(ExtraTabCompletion, ExpectElement):

    def __getattr__(self, attrname):
        """
        EXAMPLES::

            sage: mupad.package('"MuPAD-Combinat"') # optional - mupad-Combinat
            sage: S = mupad.examples.SymmetricFunctions() # optional - mupad-Combinat
            sage: type(S)                           # optional - mupad-Combinat
            <class 'sage.interfaces.mupad.MupadElement'>
            sage: S.s                               # optional - mupad-Combinat
            (examples::SymmetricFunctions(Dom::ExpressionField()))::s

            sage: x = mupad('x')                    # optional - mupad-Combinat
            sage: x.diff(x)                         # optional - mupad-Combinat
                                       1

        """
        if attrname[:1] == "_":
            if attrname not in self.__dict__:
                raise AttributeError
            else:
                return self.__dict__[attrname]
        P = self.parent()

        name = self._name + "::" + attrname
        try:
            if P.eval('type(%s::%s)'%(self.name(),attrname)).strip() == "DOM_DOMAIN":
                return P.new("%s::%s"%(self.name(),attrname))
            else:
                return MupadFunctionElement(self, name)
        except RuntimeError as err:
            if 'Unknown slot' in str(err):
                return MupadFunctionElement(self, attrname)
            else:
                raise err

    def _tab_completion(self):
        """
        EXAMPLES::

            sage: mupad.package('"MuPAD-Combinat"')       # optional - mupad-Combinat
            sage: S = mupad.examples.SymmetricFunctions() # optional - mupad-Combinat
            sage: 'HallLittlewood' in S._tab_completion()     # optional - mupad-Combinat
            True
        """
        res = self.parent().completions(self.name()+"::", strip=True)
        return res if res != [] else self.parent()._tab_completion()

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: mupad.package('"MuPAD-Combinat"') # optional - mupad-Combinat
            sage: S = mupad.examples.SymmetricFunctions() # optional - mupad-Combinat
            sage: latex(S) # optional - mupad-Combinat
            \mathrm{examples}{::}\mathrm{SymmetricFunctions}\left(\mathbb{E}\right)
        """
        self._check_valid()
        P = self.parent()
        s = P._eval_line('generate::TeX(%s)'%self.name())
        s = s.replace('\\\\','\\').strip().strip('"')
        return s

    def __len__(self):
        r"""
        The analogue in MuPAD of Python's len is the method nops

        EXAMPLES::

            sage: len(mupad([1,2,3])) # indirect doctest # optional - mupad
            3
            sage: type(len(mupad([1,2,3])))              # optional - mupad
            <... 'int'>

            sage: len(mupad(4))                          # optional - mupad
            1

        Implementing this is necessary for using MuPAD's lists as
        standard containers::

            sage: list(map(ZZ, list(mupad([1,2,3]))))    # optional - mupad
            [1, 2, 3]

            sage: [int(x) for x in mupad([1,2,3]) ]      # optional - mupad
            [1, 2, 3]

            sage: [int(x) for x in mupad("{1,2,3,5}") ]  # optional - mupad
            [1, 2, 3, 5]

        """
        return mupad.nops(self)

# An instance
mupad = Mupad()

def reduce_load_mupad():
    """
    EXAMPLES::

        sage: from sage.interfaces.mupad import reduce_load_mupad
        sage: reduce_load_mupad()
        Mupad
    """
    return mupad


def mupad_console():
    """
    Spawn a new MuPAD command-line session.

    EXAMPLES::

        sage: from sage.interfaces.mupad import mupad_console
        sage: mupad_console() #not tested

           *----*    MuPAD Pro 4.0.2 -- The Open Computer Algebra System
          /|   /|
         *----* |    Copyright (c)  1997 - 2007  by SciFace Software
         | *--|-*                   All rights reserved.
         |/   |/
         *----*      Licensed to:   ...

    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%mupad magics instead.')
    os.system('mupkern')


def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.mupad import __doctest_cleanup
        sage: m = mupad(2)         # optional - mupad
        sage: mupad.is_running()   # optional - mupad
        True
        sage: __doctest_cleanup()
        sage: mupad.is_running()   # optional - mupad
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()

