r"""
Interface to Sage

This is an expect interface to *another* copy of the Sage
interpreter.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import pickle
import re
import textwrap

from .expect import Expect, ExpectElement, FunctionElement
import sage.repl.preparse

from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.misc.persist import dumps, load
from sage.docs.instancedoc import instancedoc


class Sage(ExtraTabCompletion, Expect):
    r"""
    Expect interface to the Sage interpreter itself.

    INPUT:


    -  ``server`` - (optional); if specified runs Sage on a
       remote machine with address. You must have ssh keys setup so you
       can login to the remote machine by typing "ssh remote_machine" and
       no password, call _install_hints_ssh() for hints on how to do
       that.


    The version of Sage should be the same as on the local machine,
    since pickling is used to move data between the two Sage process.

    EXAMPLES: We create an interface to a copy of Sage. This copy of
    Sage runs as an external process with its own memory space, etc.

    ::

        sage: s = Sage()

    Create the element 2 in our new copy of Sage, and cube it.

    ::

        sage: a = s(2)
        sage: a^3
        8

    Create a vector space of dimension `4`, and compute its
    generators::

        sage: V = s('QQ^4')
        sage: V.gens()
        ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))

    Note that V is a not a vector space, it's a wrapper around an
    object (which happens to be a vector space), in another running
    instance of Sage.

    ::

        sage: type(V)
        <class 'sage.interfaces.sage0.SageElement'>
        sage: V.parent()
        Sage
        sage: g = V.0;  g
        (1, 0, 0, 0)
        sage: g.parent()
        Sage

    We can still get the actual parent by using the name attribute of
    g, which is the variable name of the object in the child process.

    ::

        sage: s('%s.parent()' % g.name())
        Vector space of dimension 4 over Rational Field

    Note that the memory space is completely different.

    ::

        sage: x = 10
        sage: s('x = 5')
        5
        sage: x
        10
        sage: s('x')
        5

    We can have the child interpreter itself make another child Sage
    process, so now three copies of Sage are running::

        sage: s3 = s('Sage()')
        sage: a = s3(10)
        sage: a
        10

    This `a=10` is in a subprocess of a subprocess of your
    original Sage.

    ::

        sage: _ = s.eval('%s.eval("x=8")' % s3.name())
        sage: s3('"x"')
        8
        sage: s('x')
        5
        sage: x
        10

    The double quotes are needed because the call to s3 first evaluates
    its arguments using the s interpreter, so the call to s3 is passed
    ``s('"x"')``, which is the string ``"x"``
    in the s interpreter.
    """
    def __init__(self,
                 logfile=None,
                 preparse=True,
                 python=False,
                 init_code=None,
                 server=None,
                 server_tmpdir=None,
                 remote_cleaner=True,
                 **kwds):
        """
        EXAMPLES::

            sage: sage0 == loads(dumps(sage0))
            True
        """

        if init_code is None:
            init_code = []
        elif isinstance(init_code, str):
            init_code = init_code.splitlines()
        else:
            try:
                init_code = list(init_code)
            except TypeError:
                raise TypeError(
                    'init_code should be a string or an iterable of lines '
                    'of code')

        if python:
            command = 'python -u'
            prompt = re.compile(b'>>> ')
            init_code.append('from sage.all import *')
        else:
            command = ' '.join([
                'sage-ipython',
                # Disable the IPython history (implemented as SQLite database)
                # to avoid problems with locking.
                '--HistoryManager.hist_file=:memory:',
                # Disable everything that prints ANSI codes
                '--colors=NoColor',
                '--no-term-title',
                '--simple-prompt',
            ])
            prompt = re.compile(b'sage: ')

        init_code.append('import pickle')
        init_code.append(textwrap.dedent("""
            def _sage0_load_local(filename):
                with open(filename, 'rb') as f:
                    return pickle.load(f)
        """))
        init_code.append(textwrap.dedent("""
            def _sage0_load_remote(filename):
                with open(filename, 'rb') as f:
                    return loads(f.read())
        """))

        Expect.__init__(self,
                        name='sage',
                        prompt=prompt,
                        command=command,
                        restart_on_ctrlc=False,
                        logfile=logfile,
                        init_code=init_code,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        remote_cleaner=remote_cleaner,
                        **kwds
                        )
        self._preparse = preparse

    def cputime(self, t=None):
        """
        Return cputime since this Sage subprocess was started.

        EXAMPLES::

            sage: sage0.cputime()     # random output
            1.3530439999999999
            sage: sage0('factor(2^157-1)')
            852133201 * 60726444167 * 1654058017289 * 2134387368610417
            sage: sage0.cputime()     # random output
            1.6462939999999999
        """
        s = self.eval('cputime(%s)' % t)
        i = s.rfind('m')
        if i != -1:
            s = s[i + 1:-1]
        return float(s)

    def _tab_completion(self):
        """
        EXAMPLES::

            sage: t = sage0._tab_completion()
            sage: len(t) > 100
            True
            sage: 'gcd' in t
            True
        """
        return eval(self.eval('print(repr(list(globals())))'))

    def __call__(self, x):
        """
        EXAMPLES::

            sage: a = sage0(4)
            sage: a.parent()
            Sage
            sage: a is sage0(a)
            True

        TESTS::

            sage: sage0(axiom(x^2+1)) #optional - axiom
            x^2 + 1

        """
        if isinstance(x, ExpectElement):
            if x.parent() is self:
                return x
            else:
                return self(x.sage())

        if isinstance(x, str):
            return SageElement(self, x)

        if self.is_local():
            with open(self._local_tmpfile(), 'wb') as fobj:
                fobj.write(pickle.dumps(x, 2))
            code = '_sage0_load_local({!r})'.format(self._local_tmpfile())
            return SageElement(self, code)
        else:
            with open(self._local_tmpfile(), 'wb') as fobj:
                fobj.write(dumps(x))   # my dumps is compressed by default
            self._send_tmpfile_to_server()
            code = '_sage0_load_remote({!r})'.format(self._remote_tmpfile())
            return SageElement(self, code)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: sage0.__reduce__()
            (<function reduce_load_Sage at 0x...>, ())
        """
        return reduce_load_Sage, tuple([])

    def _quit_string(self):
        """
        EXAMPLES::

            sage: sage0._quit_string()
            ''
        """
        return ""

    def preparse(self, x):
        """
        Returns the preparsed version of the string s.

        EXAMPLES::

            sage: sage0.preparse('2+2')
            'Integer(2)+Integer(2)'
        """
        return sage.repl.preparse.preparse(x)

    def eval(self, line, strip=True, **kwds):
        """
        Send the code x to a second instance of the Sage interpreter and
        return the output as a string.

        This allows you to run two completely independent copies of Sage at
        the same time in a unified way.

        INPUT:


        -  ``line`` - input line of code

        -  ``strip`` - ignored


        EXAMPLES::

            sage: sage0.eval('2+2')
            '4'
        """
        if self._preparse:
            line = self.preparse(line)
        line = self._wrap_multiline(line)
        return Expect.eval(self, line, **kwds).strip()

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: sage0.set('x', '2')
            sage: sage0.get('x')
            '2'
        """
        cmd = '%s=%s' % (var, value)
        out = self.eval(cmd)
        if 'Traceback' in out:
            raise TypeError("Error executing code in Sage\nCODE:\n\t%s\nSage ERROR:\n\t%s" % (cmd, out))

    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: sage0.set('x', '2')
            sage: sage0.get('x')
            '2'
        """
        return self.eval('print(%s)' % var).strip()

    def clear(self, var):
        """
        Clear the variable named var.

        Note that the exact format of the NameError for a cleared variable
        is slightly platform dependent, see :trac:`10539`.

        EXAMPLES::

            sage: sage0.set('x', '2')
            sage: sage0.get('x')
            '2'
            sage: sage0.clear('x')
            sage: 'NameError' in sage0.get('x')
            True
        """
        self.eval('del %s' % var)

    def _contains(self, v1, v2):
        """
        EXAMPLES::

            sage: sage0._contains('2', 'QQ')
            True
        """
        return self.eval('%s in %s' % (v1, v2)) == "True"

    def console(self):
        """
        Spawn a new Sage command-line session.

        EXAMPLES::

            sage: sage0.console() #not tested
            ----------------------------------------------------------------------
            | SageMath version ..., Release Date: ...                            |
            | Using Python ....   Type "help()" for help.                        |
            ----------------------------------------------------------------------
            ...
        """
        sage0_console()

    def version(self):
        """
        EXAMPLES::

            sage: sage0.version()
            'SageMath version ..., Release Date: ...'
            sage: sage0.version() == version()
            True
        """
        return sage0_version()

    def _object_class(self):
        """
        EXAMPLES::

            sage: sage0._object_class()
            <class 'sage.interfaces.sage0.SageElement'>
        """
        return SageElement

    def new(self, x):
        """
        EXAMPLES::

            sage: sage0.new(2)
            2
            sage: _.parent()
            Sage
        """
        return SageElement(self, x)

    @staticmethod
    def _wrap_multiline(s):
        r"""
        The Sage interface does not currently handle multi-line Python
        statements well.

        So given a multi-line Python statement, it is converted to an
        equivalent ``eval()`` call.

        EXAMPLES::

            sage: import textwrap
            sage: code = textwrap.dedent('''
            ....:     def foo():
            ....:         return 'foo'
            ....: ''')
            sage: print(sage0._wrap_multiline(code))
            eval(compile("def foo():\n    return 'foo'", '<stdin>', 'single'))
            sage: sage0.eval(code)
            ''
            sage: sage0.eval('foo()')
            "'foo'"
        """

        if '\n' in s:
            return "eval(compile({!r}, '<stdin>', 'single'))".format(s.strip())
        else:
            return s


@instancedoc
class SageElement(ExpectElement):

    def _rich_repr_(self, display_manager, **kwds):
        """
        Disable rich output

        This is necessary because otherwise our :meth:`__getattr__`
        would be called.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: m = sage0(4)
            sage: m._rich_repr_(get_display_manager()) is None
            True
        """
        return None

    def _repr_option(self, option):
        """
        Disable repr option.

        This is necessary because otherwise our :meth:`__getattr__`
        would be called.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: m = sage0(4)
            sage: m._repr_option('ascii_art')
            False
        """
        return False

    def __getattr__(self, attrname):
        """
        EXAMPLES::

            sage: m = sage0(4)
            sage: four_gcd = m.gcd
            sage: type(four_gcd)
            <class 'sage.interfaces.sage0.SageFunction'>
        """
        self._check_valid()
        return SageFunction(self, attrname)

    def _sage_(self):
        """
        Return local copy of self.

        EXAMPLES::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F == sage0(F)._sage_()
            True
        """
        P = self.parent()
        if P.is_remote():
            P.eval('save({}, {!r})'.format(self.name(), P._remote_tmpfile()))
            P._get_tmpfile_from_server(self)
            return load(P._local_tmp_file())
        else:
            P.eval('save({}, {!r})'.format(self.name(), P._local_tmpfile()))
            return load(P._local_tmpfile())


@instancedoc
class SageFunction(FunctionElement):
    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: four_gcd = sage0(4).gcd
            sage: four_gcd(6)
            2
        """
        P = self._obj.parent()

        # Important! Keep references to the argument values or else they may
        # get cleared from the interpreter before we complete the function
        # call.
        args = [P(x) for x in args]
        kwds = [(k, P(v)) for k, v in kwds.items()]

        arg_str = ','.join(x.name() for x in args)
        kwd_str = ','.join('%s=%s' % (k, v.name()) for k, v in kwds)

        if arg_str and kwd_str:
            args_str = '%s,%s' % (arg_str, kwd_str)
        else:
            args_str = arg_str + kwd_str  # At least one of these is empty

        call_str = '%s.%s(%s)' % (self._obj._name, self._name, args_str)
        return SageElement(P, call_str)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sage0(4).gcd
            <built-in method gcd of sage.rings.integer.Integer object at 0x...>
        """
        return str(self._obj.parent().eval('%s.%s' % (self._obj._name,
                                                      self._name)))


sage0 = Sage()


def reduce_load_Sage():
    """
    EXAMPLES::

        sage: from sage.interfaces.sage0 import reduce_load_Sage
        sage: reduce_load_Sage()
        Sage
    """
    return sage0


def reduce_load_element(s):
    """
    EXAMPLES::

        sage: from sage.interfaces.sage0 import reduce_load_element
        sage: s = dumps(1/2)
        sage: half = reduce_load_element(s); half
        1/2
        sage: half.parent()
        Sage
    """
    import base64
    s = base64.b32encode(s)
    sage0.eval('import base64')
    return sage0('loads(base64.b32decode({!r}))'.format(s))


def sage0_console():
    """
    Spawn a new Sage command-line session.

    EXAMPLES::

        sage: sage0_console() #not tested
        ----------------------------------------------------------------------
        | SageMath version ..., Release Date: ...                            |
        | Using Python ....   Type "help()" for help.                        |
        ----------------------------------------------------------------------
        ...
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%sage0 magics instead.')
    os.system('sage')


def sage0_version():
    """
    EXAMPLES::

        sage: from sage.interfaces.sage0 import sage0_version
        sage: sage0_version() == version()
        True
    """
    return str(sage0('version()'))
