# sage.doctest: needs sage.libs.pari
r"""
Interface to the GP calculator of PARI/GP

Type ``gp.[tab]`` for a list of all the functions
available from your Gp install. Type ``gp.[tab]?`` for
Gp's help about a given function. Type ``gp(...)`` to
create a new Gp object, and ``gp.eval(...)`` to evaluate a
string using Gp (and get the result back as a string).

EXAMPLES: We illustrate objects that wrap GP objects (gp is the
PARI interpreter)::

    sage: M = gp('[1,2;3,4]')
    sage: M
    [1, 2; 3, 4]
    sage: M * M
    [7, 10; 15, 22]
    sage: M + M
    [2, 4; 6, 8]
    sage: M.matdet()
    -2

::

    sage: E = gp.ellinit([1,2,3,4,5])
    sage: E.ellglobalred()
    [10351, [1, -1, 0, -1], 1, [11, 1; 941, 1], [[1, 5, 0, 1], [1, 5, 0, 1]]]
    sage: E.ellan(20)
    [1, 1, 0, -1, -3, 0, -1, -3, -3, -3, -1, 0, 1, -1, 0, -1, 5, -3, 4, 3]

::

    sage: primitive_root(7)
    3
    sage: x = gp("znlog( Mod(2,7), Mod(3,7))")
    sage: 3^x % 7
    2

::

    sage: print(gp("taylor(sin(x),x)"))
    x - 1/6*x^3 + 1/120*x^5 - 1/5040*x^7 + 1/362880*x^9 - 1/39916800*x^11 + 1/6227020800*x^13 - 1/1307674368000*x^15 + O(x^16)

GP has a powerful very efficient algorithm for numerical
computation of integrals.

::

    sage: gp("a = intnum(x=0,6,sin(x))")
    0.03982971334963397945434770208               # 32-bit
    0.039829713349633979454347702077075594548     # 64-bit
    sage: gp("a")
    0.03982971334963397945434770208               # 32-bit
    0.039829713349633979454347702077075594548     # 64-bit
    sage: gp.kill("a")
    sage: gp("a")
    a

Note that gp ASCII plots *do* work in Sage, as follows::

    sage: print(gp.eval("plot(x=0,6,sin(x))"))
    <BLANKLINE>
    0.9988963 |''''''''''''_x...x_''''''''''''''''''''''''''''''''''''''''''|
              |          x"        "x                                        |
              |        _"            "_                                      |
              |       x                x                                     |
              |      "                  "                                    |
              |     "                    "                                   |
              |   _"                      "_                                 |
              |  _                          _                                |
              | _                            _                               |
              |_                              _                              |
              _                                                              |
              `````````````````````````````````"``````````````````````````````
              |                                 "                            |
              |                                  "                           |
              |                                   "                          "
              |                                    "_                      _"|
              |                                      _                    _  |
              |                                       _                  _   |
              |                                        x                x    |
              |                                         "_            _"     |
              |                                           x_        _x       |
    -0.998955 |............................................."x____x".........|
              0                                                              6

The GP interface reads in even very long input (using files) in a
robust manner, as long as you are creating a new object.

::

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = gp.eval(t)
    sage: a = gp(t)

In Sage, the PARI large Galois groups datafiles should be installed
by default::

    sage: f = gp('x^9 - x - 2')
    sage: f.polgalois()
    [362880, -1, 34, "S9"]

TESTS:

Test error recovery::

    sage: x = gp('1/0')
    Traceback (most recent call last):
    ...
    TypeError: Error executing code in GP:
    CODE:
        sage[...]=1/0;
    PARI/GP ERROR:
      ***   at top-level: sage[...]=1/0
      ***                          ^--
      *** _/_: impossible inverse in gdiv: 0.

AUTHORS:

- William Stein

- David Joyner: some examples

- William Stein (2006-03-01): added tab completion for methods:
  gp.[tab] and x = gp(blah); x.[tab]

- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta

- William Stein (2006-05-17): updated to work with PARI 2.2.13-beta
"""

##########################################################################
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#
##########################################################################
import os

import sage.interfaces.abc

from sage.env import DOT_SAGE
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.libs.pari.all import pari
from sage.misc.instancedoc import instancedoc
from sage.misc.lazy_import import lazy_import
from sage.misc.verbose import verbose

lazy_import('sage.rings.cc', 'CC')

from .expect import Expect, ExpectElement, ExpectFunction, FunctionElement


class Gp(ExtraTabCompletion, Expect):
    """
    Interface to the PARI gp interpreter.

    Type ``gp.[tab]`` for a list of all the functions
    available from your Gp install. Type ``gp.[tab]?`` for
    Gp's help about a given function. Type ``gp(...)`` to
    create a new Gp object, and ``gp.eval(...)`` to evaluate a
    string using Gp (and get the result back as a string).

        INPUT:

        - ``stacksize`` -- integer (default: 10000000); the initial PARI
          stacksize in bytes (default: 10MB)
        - ``script_subdirectory`` -- string (default: ``None``); name of the
          subdirectory of ``SAGE_EXTCODE/pari`` from which to read scripts
        - ``logfile`` -- string (default: ``None``); log file for the pexpect
          interface
        - ``server`` -- name of remote server
        - ``server_tmpdir`` -- name of temporary directory on remote server
        - ``init_list_length`` -- integer (default: 1024); length of initial
          list of local variables
        - ``seed`` -- integer (default: random); random number generator seed
          for pari

        EXAMPLES::

            sage: Gp()
            PARI/GP interpreter
    """
    def __init__(self, stacksize=10000000,   # 10MB
                 maxread=None, script_subdirectory=None,
                 logfile=None,
                 server=None,
                 server_tmpdir=None,
                 init_list_length=1024,
                 seed=None):
        """
        Initialization of this PARI gp interpreter.

        INPUT:

        - ``stacksize`` -- integer (default: 10000000); the initial PARI
          stacksize in bytes (default: 10MB)
        - ``script_subdirectory`` -- string (default: ``None``); name of the
          subdirectory of SAGE_EXTCODE/pari from which to read scripts
        - ``logfile`` -- string (default: ``None``); log file for the pexpect
          interface
        - ``server`` -- name of remote server
        - ``server_tmpdir`` -- name of temporary directory on remote server
        - ``init_list_length`` -- integer (default: 1024); length of initial
          list of local variables.
        - ``seed`` -- integer (default random nonzero 31 bit integer); value of
          random seed

        EXAMPLES::

            sage: gp == loads(dumps(gp))
            True
        """
        Expect.__init__(self,
                        name='pari',
                        prompt='\\? ',
                        # --fast so the system gprc isn't read (we configure below)
                        command=f"gp --fast --emacs --quiet --stacksize {stacksize}",
                        maxread=maxread,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=False,
                        verbose_start=False,
                        logfile=logfile,
                        eval_using_file_cutoff=1024)
        self.__seq = 0
        self.__var_store_len = 0
        self.__init_list_length = init_list_length
        self._seed = seed

    def set_seed(self, seed=None):
        """
        Set the seed for gp interpreter.

        The seed should be an integer.

        EXAMPLES::

            sage: g = Gp()
            sage: g.set_seed(1)
            1
            sage: [g.random() for i in range(5)]
            [1546275796, 879788114, 1745191708, 771966234, 1247963869]
        """
        if seed is None:
            seed = self.rand_seed()
        self.eval("setrand(%d)" % seed)
        self._seed = seed
        return seed

    def _start(self, alt_message=None, block_during_init=True):
        Expect._start(self, alt_message, block_during_init)
        # disable memory debugging: those warnings can only confuse our
        # interface
        self._eval_line('default(debugmem,0);')
        # disable timer
        self._eval_line('default(timer,0);')
        # disable the break loop, otherwise gp will seem to hang on errors
        self._eval_line('default(breakloop,0);')
        # list of directories where gp will look for scripts (only current working directory)
        self._eval_line('default(path,".");')
        # executable for gp ?? help
        self._eval_line('default(help, "gphelp -detex");')
        # logfile disabled since Expect already logs
        self._eval_line('default(log,0);')
        self._eval_line("default(nbthreads,1);")
        # set random seed
        self.set_seed(self._seed)

    def _repr_(self):
        """
        String representation of this PARI gp interpreter.

        EXAMPLES::

            sage: gp  # indirect doctest
            PARI/GP interpreter
        """
        return 'PARI/GP interpreter'

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.gp import gp
            sage: gp.__reduce__()
            (<function reduce_load_GP at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            PARI/GP interpreter
        """
        return reduce_load_GP, tuple([])

    def _function_class(self):
        """
        Return the GpFunction class.

        EXAMPLES::

            sage: gp._function_class()
            <class 'sage.interfaces.expect.ExpectFunction'>
            sage: type(gp.gcd)
            <class 'sage.interfaces.expect.ExpectFunction'>
        """
        return GpFunction

    def _quit_string(self):
        """
        Return the string used to quit the GP interpreter.

        EXAMPLES::

            sage: gp._quit_string()
            '\\q'

        ::

            sage: g = Gp()
            sage: a = g(2)
            sage: g.is_running()
            True
            sage: g.quit()
            sage: g.is_running()
            False
        """
        return r"\q"

    def _read_in_file_command(self, filename):
        r"""
        Return the string used to read filename into GP.

        EXAMPLES::

            sage: gp._read_in_file_command('test')
            'read("test")'

        ::

            sage: filename = tmp_filename()
            sage: f = open(filename, 'w')
            sage: _ = f.write('x = 22;\n')
            sage: f.close()
            sage: gp.read(filename)
            sage: gp.get('x').strip()
            '22'
        """
        return 'read("%s")' % filename

    def _tab_completion(self):
        """
        EXAMPLES::

            sage: c = gp._tab_completion()
            sage: len(c) > 100
            True
            sage: 'gcd' in c
            True
        """
        try:
            b = self.__builtin
        except AttributeError:
            b = self.eval('?*').split()
            self.__builtin = b
        return b + self.eval('?0').split()

    def get_precision(self):
        """
        Return the current PARI precision for real number computations.

        EXAMPLES::

            sage: gp.get_precision()
            28              # 32-bit
            38              # 64-bit
        """
        return self.get_default('realprecision')

    get_real_precision = get_precision

    def set_precision(self, prec):
        """
        Set the PARI precision (in decimal digits) for real
        computations, and returns the old value.

        .. NOTE::

            PARI/GP rounds up precisions to the nearest machine word,
            so the result of :meth:`get_precision` is not always the
            same as the last value inputted to :meth:`set_precision`.

        EXAMPLES::

            sage: old_prec = gp.set_precision(53); old_prec
            28              # 32-bit
            38              # 64-bit
            sage: gp.get_precision()
            57
            sage: gp.set_precision(old_prec)
            57
            sage: gp.get_precision()
            28              # 32-bit
            38              # 64-bit
        """
        return self.set_default('realprecision', prec)

    set_real_precision = set_precision

    def get_series_precision(self):
        """
        Return the current PARI power series precision.

        EXAMPLES::

            sage: gp.get_series_precision()
            16
        """
        return self.get_default('seriesprecision')

    def set_series_precision(self, prec=None):
        """
        Set the PARI power series precision, and returns the old precision.

        EXAMPLES::

            sage: old_prec = gp.set_series_precision(50); old_prec
            16
            sage: gp.get_series_precision()
            50
            sage: gp.set_series_precision(old_prec)
            50
            sage: gp.get_series_precision()
            16
        """
        return self.set_default('seriesprecision', prec)

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: gp._eval_line('2+2')
            '4'

        TESTS:

        We verify that :issue:`11617` is fixed::

            sage: gp._eval_line('a='+str(list(range(2*10^5))))[:70]
            '[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,'
        """
        line = line.strip()
        if len(line) == 0:
            return ''
        a = Expect._eval_line(self, line,
                              allow_use_file=allow_use_file,
                              wait_for_prompt=wait_for_prompt)
        if a.find("the PARI stack overflows") != -1:
            verbose("automatically doubling the PARI stack and re-executing current input line")
            b = self.eval("allocatemem()")
            if b.find("Warning: not enough memory") != -1:
                raise RuntimeError(a)
            return self._eval_line(line, allow_use_file=allow_use_file,
                                   wait_for_prompt=wait_for_prompt)
        else:
            return a

    def cputime(self, t=None):
        """
        cputime for pari - cputime since the pari process was started.

        INPUT:

        - ``t`` -- (default: ``None``) if not None, then returns
          time since t


        .. warning::

           If you call gettime explicitly, e.g., gp.eval('gettime'),
           you will throw off this clock.

        EXAMPLES::

            sage: gp.cputime()          # random output
            0.0080000000000000002
            sage: gp.factor('2^157-1')
            [852133201, 1; 60726444167, 1; 1654058017289, 1; 2134387368610417, 1]
            sage: gp.cputime()          # random output
            0.26900000000000002
        """
        try:
            tm = self._last
        except AttributeError:
            tm = 0.0
        m = eval(self.eval('gettime()/1000.0')) + tm
        self._last = m
        if t:
            return m - t
        return m

    def set_default(self, var, value):
        """
        Set a PARI gp configuration variable, and return the old value.

        INPUT:

        - ``var`` -- string; the name of a PARI gp
          configuration variable (see ``gp.default()`` for a list)
        - ``value`` -- the value to set the variable to

        EXAMPLES::

            sage: old_prec = gp.set_default('realprecision', 110)
            sage: gp.get_default('realprecision')
            115
            sage: gp.set_default('realprecision', old_prec)
            115
            sage: gp.get_default('realprecision')
            28              # 32-bit
            38              # 64-bit
        """
        old = self.get_default(var)
        self._eval_line('default(%s,%s)' % (var, value))
        return old

    def get_default(self, var):
        """
        Return the current value of a PARI gp configuration variable.

        INPUT:

        - ``var`` -- string; the name of a PARI gp
          configuration variable (see ``gp.default()`` for a list)

        OUTPUT: string; the value of the variable

        EXAMPLES::

            sage: gp.get_default('log')
            0
            sage: gp.get_default('datadir')
            '.../share/pari'
            sage: gp.get_default('seriesprecision')
            16
            sage: gp.get_default('realprecision')
            28              # 32-bit
            38              # 64-bit
        """
        return eval(self._eval_line('default(%s)' % var))

    def set(self, var, value):
        """
        Set the GP variable var to the given value.

        INPUT:

        - ``var`` -- string; a valid GP variable identifier
        - ``value`` -- a value for the variable

        EXAMPLES::

            sage: gp.set('x', '2')
            sage: gp.get('x')
            '2'
        """
        cmd = '%s=%s;' % (var, value)
        out = self.eval(cmd)
        if out.find('***') != -1:
            raise TypeError("Error executing code in GP:\nCODE:\n\t%s\nPARI/GP ERROR:\n%s" % (cmd, out))

    def get(self, var):
        """
        Get the value of the GP variable var.

        INPUT:

        - ``var`` -- string; a valid GP variable identifier

        EXAMPLES::

            sage: gp.set('x', '2')
            sage: gp.get('x')
            '2'
        """
        return self.eval('print(%s)' % var)

    def kill(self, var):
        """
        Kill the value of the GP variable var.

        INPUT:

        - ``var`` -- string; a valid GP variable identifier

        EXAMPLES::

            sage: gp.set('xx', '22')
            sage: gp.get('xx')
            '22'
            sage: gp.kill('xx')
            sage: gp.get('xx')
            'xx'
        """
        self.eval('kill(%s)' % var)

    def _next_var_name(self):
        """
        Return the name of the next unused interface variable name.

        EXAMPLES::

            sage: g = Gp()
            sage: g._next_var_name()
            'sage[1]'
            sage: g(2)^2
            4
            sage: g._next_var_name()
            'sage[5]'

        TESTS:

        The vector of results is correctly resized when the stack has
        to be enlarged during this operation::

            sage: g = Gp(stacksize=3*10^6,init_list_length=12000)  # long time
            sage: for n in [1..13000]:  # long time
            ....:     a = g(n)
            sage: g('length(sage)')     # long time
            24000
        """
        self.__seq += 1
        if self.__seq >= self.__var_store_len:
            if self.__var_store_len == 0:
                self.eval('sage=vector(%s,k,0);' % self.__init_list_length)
                self.__var_store_len = self.__init_list_length
            else:
                self.eval('sage0=concat(sage, vector(%s,k,0));' % self.__var_store_len)
                self.eval('sage=sage0;')
                self.eval('kill(sage0);')
                self.__var_store_len *= 2
                verbose("doubling PARI/sage object vector: %s" % self.__var_store_len)
        return 'sage[%s]' % self.__seq

    def _reset_expect(self):
        """
        Reset state of the GP interface.

        EXAMPLES::

            sage: a = gp('10'); a
            10
            sage: gp.quit()  # indirect doctest
            sage: a
            (invalid PARI/GP interpreter object -- The pari session in which this object was defined is no longer running.)
            sage: gp("30!")
            265252859812191058636308480000000
        """
        self.__var_store_len = 0
        Expect._reset_expect(self)

    def console(self):
        """
        Spawn a new GP command-line session.

        EXAMPLES::

            sage: gp.console()  # not tested
            GP/PARI CALCULATOR Version 2.4.3 (development svn-12577)
            amd64 running linux (x86-64/GMP-4.2.1 kernel) 64-bit version
            compiled: Jul 21 2010, gcc-4.6.0 20100705 (experimental) (GCC)
            (readline v6.0 enabled, extended help enabled)
        """
        gp_console()

    def version(self):
        """
        Return the version of GP being used.

        EXAMPLES::

            sage: gp.version()  # not tested
            ((2, 4, 3), 'GP/PARI CALCULATOR Version 2.4.3 (development svn-12577)')
        """
        return gp_version()

    def _object_class(self):
        """
        Return the GpElement class.

        EXAMPLES::

            sage: gp._object_class()
            <class 'sage.interfaces.gp.GpElement'>
            sage: type(gp(2))
            <class 'sage.interfaces.gp.GpElement'>
        """
        return GpElement

    def _function_element_class(self):
        """
        Return the GpFunctionElement class.

        EXAMPLES::

            sage: gp._function_element_class()
            <class 'sage.interfaces.expect.FunctionElement'>
            sage: type(gp(2).gcd)
            <class 'sage.interfaces.expect.FunctionElement'>
        """
        return GpFunctionElement

    def _true_symbol(self):
        """
        Return the symbol used for truth in GP.

        EXAMPLES::

            sage: gp._true_symbol()
            '1'

        ::

            sage: gp(2) == gp(2)
            True
        """
        return '1'

    def _false_symbol(self):
        """
        Return the symbol used for falsity in GP.

        EXAMPLES::

            sage: gp._false_symbol()
            '0'

        ::

            sage: gp(2) == gp(3)
            False
        """
        return '0'

    def _equality_symbol(self):
        """
        Return the symbol used for equality in GP.

        EXAMPLES::

            sage: gp._equality_symbol()
            '=='

        ::

            sage: gp(2) == gp(2)
            True
        """
        return '=='

    def _exponent_symbol(self):
        """
        Return the symbol to denote the exponent of a number in GP.

        EXAMPLES::

            sage: gp._exponent_symbol()
            ' E'

        ::

            sage: repr(gp(10.^80)).replace(gp._exponent_symbol(), 'e')
            '1.0000000000000000000000000000000000000e80'    # 64-bit
            '1.000000000000000000000000000e80'              # 32-bit
        """
        return ' E'

    def help(self, command):
        r"""
        Return GP's help for ``command``.

        EXAMPLES::

            sage: gp.help('gcd')
            'gcd(x,{y}): greatest common divisor of x and y.'
        """
        return self.eval('?%s' % command).strip()

    def new_with_bits_prec(self, s, precision=0):
        r"""
        Create a GP object from s with ``precision`` bits of
        precision. GP actually automatically increases this precision to
        the nearest word (i.e. the next multiple of 32 on a 32-bit machine,
        or the next multiple of 64 on a 64-bit machine).

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: pi_def = gp(pi); pi_def
            3.141592653589793238462643383                  # 32-bit
            3.1415926535897932384626433832795028842        # 64-bit
            sage: pi_def.precision()
            28                                             # 32-bit
            38                                             # 64-bit
            sage: pi_150 = gp.new_with_bits_prec(pi, 150)
            sage: new_prec = pi_150.precision(); new_prec
            48                                             # 32-bit
            57                                             # 64-bit
            sage: old_prec = gp.set_precision(new_prec); old_prec
            28                                             # 32-bit
            38                                             # 64-bit
            sage: pi_150
            3.14159265358979323846264338327950288419716939938  # 32-bit
            3.14159265358979323846264338327950288419716939937510582098  # 64-bit
            sage: gp.set_precision(old_prec)
            48                                             # 32-bit
            57                                             # 64-bit
            sage: gp.get_precision()
            28                                             # 32-bit
            38                                             # 64-bit
        """
        if precision:
            old_prec = self.get_real_precision()
            prec = int(precision / 3.321928095)
            self.set_real_precision(prec)
            x = self(s)
            self.set_real_precision(old_prec)
        else:
            x = self(s)
        return x


@instancedoc
class GpElement(ExpectElement, sage.interfaces.abc.GpElement):
    """
    EXAMPLES: This example illustrates dumping and loading GP elements
    to compressed strings.

    ::

        sage: a = gp(39393)
        sage: loads(a.dumps()) == a
        True

    Since dumping and loading uses the string representation of the
    object, it need not result in an identical object from the point of
    view of PARI::

        sage: E = gp('ellinit([1,2,3,4,5])')
        sage: loads(dumps(E)) == E
        True
        sage: x = gp.Pi()/3
        sage: loads(dumps(x)) == x
        False
        sage: x
        1.047197551196597746154214461            # 32-bit
        1.0471975511965977461542144610931676281  # 64-bit
        sage: loads(dumps(x))
        1.047197551196597746154214461            # 32-bit
        1.0471975511965977461542144610931676281  # 64-bit

    The two elliptic curves look the same, but internally the floating
    point numbers are slightly different.
    """
    def _reduce(self):
        """
        Return the string representation of self, for pickling.

        Because the internal representation of a gp element is richer
        than the corresponding sage object, we use the string representation
        for pickling.

        EXAMPLES::

            sage: E = gp('ellinit([1,2,3,4,5])')
            sage: loads(dumps(E)) == E # indirect doctest
            True
            sage: gp(E.sage()) == E
            False
        """
        return repr(self)

    def _sage_(self):
        """
        Convert this GpElement into a Sage object, if possible.

        EXAMPLES::

            sage: gp(SR(I)).sage()                                                      # needs sage.symbolic
            i
            sage: gp(SR(I)).sage().parent()                                             # needs sage.symbolic
            Number Field in i with defining polynomial x^2 + 1 with i = 1*I

        ::

            sage: # needs sage.modules
            sage: M = Matrix(ZZ,2,2,[1,2,3,4]); M
            [1 2]
            [3 4]
            sage: gp(M)
            [1, 2; 3, 4]
            sage: gp(M).sage()
            [1 2]
            [3 4]
            sage: gp(M).sage() == M
            True

        Conversion of strings::

           sage: s = gp('"foo"')
           sage: s.sage()
           'foo'
           sage: type(s.sage())
           <class 'str'>
        """
        if self.is_string():
            return str(self)
        return pari(str(self)).sage()

    def is_string(self):
        """
        Tell whether this element is a string.

        EXAMPLES::

            sage: gp('"abc"').is_string()
            True
            sage: gp('[1,2,3]').is_string()
            False
        """
        return repr(self.type()) == 't_STR'

    def __float__(self):
        """
        Return Python float.

        EXAMPLES::

            sage: float(gp(10))
            10.0
        """
        return float(pari(str(self)))

    def __bool__(self):
        """
        EXAMPLES::

            sage: gp(2).bool()
            True
            sage: bool(gp(2))
            True
            sage: bool(gp(0))
            False
        """
        P = self._check_valid()
        return P.eval('%s != 0' % (self.name())) == '1'

    def _complex_mpfr_field_(self, CC):
        """
        Return ComplexField element of ``self``.

        INPUT:

        - ``CC`` -- a Complex or Real Field

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: z = gp(SR(1+15*I)); z
            1 + 15*I
            sage: z._complex_mpfr_field_(CC)
            1.00000000000000 + 15.0000000000000*I
            sage: CC(z) # CC(gp(1+15*I))
            1.00000000000000 + 15.0000000000000*I
            sage: CC(gp(11243.9812+15*I))
            11243.9812000000 + 15.0000000000000*I
            sage: ComplexField(10)(gp(11243.9812+15*I))
            11000. + 15.*I
        """
        # Multiplying by CC(1) is necessary here since
        # sage: pari(gp(1+I)).sage().parent()
        # Maximal Order in Number Field in i with defining polynomial x^2 + 1
        return CC((CC.one() * pari(self)).sage())

    def _complex_double_(self, CDF):
        """
        Return this value as a CDF element.

        EXAMPLES::

            sage: CDF(gp(pi+I*e))                                                       # needs sage.symbolic
            3.141592653589793 + 2.718281828459045*I
        """
        # Retrieving values from another computer algebra system is
        # slow anyway, right?
        cc_val = self._complex_mpfr_field_(CC)
        return CDF(cc_val)

    def __len__(self):
        """
        EXAMPLES::

            sage: len(gp([1,2,3]))
            3
        """
        return int(self.length())

    def __del__(self):
        """
        Note that clearing object is pointless since it wastes time and
        PARI/GP doesn't really free used memory.

        EXAMPLES::

            sage: a = gp(2)
            sage: a.__del__()
            sage: a
            2
            sage: del a
            sage: a
            Traceback (most recent call last):
            ...
            NameError: name 'a' is not defined
        """
        return

    # This is tempting -- but the (la)tex output is very very
    # out of date, e.g., for matrices it uses \pmatrix (which
    # causes an error if amsmath is loaded) and for rationals
    # it does nothing, etc.
    # def _latex_(self):
    #     P = self._check_valid()
    #     return P.eval('printtex(%s)'%self.name())

    def _tab_completion(self):
        """
        EXAMPLES::

            sage: 'gcd' in gp(2)._tab_completion()
            True
        """
        return self.parent()._tab_completion()


GpFunctionElement = FunctionElement
GpFunction = ExpectFunction


def is_GpElement(x):
    """
    Return ``True`` if ``x`` is of type :class:`GpElement`.

    This function is deprecated; use :func:`isinstance`
    (of :class:`sage.interfaces.abc.GpElement`) instead.

    EXAMPLES::

        sage: from sage.interfaces.gp import is_GpElement
        sage: is_GpElement(gp(2))
        doctest:...: DeprecationWarning: the function is_GpElement is deprecated; use isinstance(x, sage.interfaces.abc.GpElement) instead
        See https://github.com/sagemath/sage/issues/34804 for details.
        True
        sage: is_GpElement(2)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(34804, "the function is_GpElement is deprecated; use isinstance(x, sage.interfaces.abc.GpElement) instead")

    return isinstance(x, GpElement)


# An instance
gp = Gp(logfile=os.path.join(DOT_SAGE, 'gp-expect.log'))
# useful for debugging!


def reduce_load_GP():
    """
    Return the GP interface object defined in ``sage.interfaces.gp``.

    EXAMPLES::

        sage: from sage.interfaces.gp import reduce_load_GP
        sage: reduce_load_GP()
        PARI/GP interpreter
    """
    return gp


def gp_console():
    """
    Spawn a new GP command-line session.

    EXAMPLES::

        sage: gp.console()  # not tested
        GP/PARI CALCULATOR Version 2.4.3 (development svn-12577)
        amd64 running linux (x86-64/GMP-4.2.1 kernel) 64-bit version
        compiled: Jul 21 2010, gcc-4.6.0 20100705 (experimental) (GCC)
        (readline v6.0 enabled, extended help enabled)
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%gp magics instead.')
    os.system('gp')


def gp_version():
    """
    EXAMPLES::

        sage: gp.version()  # not tested
        ((2, 4, 3), 'GP/PARI CALCULATOR Version 2.4.3 (development svn-12577)')
    """
    v = gp.eval(r'\v')
    i = v.find("Version ")
    w = v[i + len("Version "):]
    i = w.find(' ')
    w = w[:i]
    t = tuple([int(n) for n in w.split('.')])
    k = v.find('\n')
    return t, v[:k].strip()
