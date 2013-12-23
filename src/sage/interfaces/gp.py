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
    [10351, [1, -1, 0, -1], 1]
    sage: E.ellan(20)
    [1, 1, 0, -1, -3, 0, -1, -3, -3, -3, -1, 0, 1, -1, 0, -1, 5, -3, 4, 3]

::

    sage: primitive_root(7)
    3
    sage: x = gp("znlog( Mod(2,7), Mod(3,7))")
    sage: 3^x % 7
    2

::

    sage: print gp("taylor(sin(x),x)")
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

    sage: print gp.eval("plot(x=0,6,sin(x))")
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

In Sage, the PARI large galois groups datafiles should be installed
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
      *** _/_: division by zero

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
#                  http://www.gnu.org/licenses/
#
##########################################################################

from expect import Expect, ExpectElement, ExpectFunction, FunctionElement
from sage.misc.misc import verbose
from sage.libs.pari.all import pari
import sage.rings.complex_field
## import sage.rings.all

class Gp(Expect):
    """
    Interface to the PARI gp interpreter.

    Type ``gp.[tab]`` for a list of all the functions
    available from your Gp install. Type ``gp.[tab]?`` for
    Gp's help about a given function. Type ``gp(...)`` to
    create a new Gp object, and ``gp.eval(...)`` to evaluate a
    string using Gp (and get the result back as a string).

        INPUT:

        - ``stacksize`` (int, default 10000000) -- the initial PARI
          stacksize in bytes (default 10MB)
        - ``maxread`` (int, default 100000) -- ??
        - ``script_subdirectory`` (string, default None) -- name of the subdirectory of SAGE_EXTCODE/pari from which to read scripts
        - ``logfile`` (string, default None) -- log file for the pexpect interface
        - ``server`` -- name of remote server
        - ``server_tmpdir`` -- name of temporary directory on remote server
        - ``init_list_length`` (int, default 1024) -- length of initial list of local variables.

        EXAMPLES::

            sage: Gp()
            PARI/GP interpreter
    """
    def __init__(self, stacksize=10000000,   # 10MB
                 maxread=100000, script_subdirectory=None,
                 logfile=None,
                 server=None,
                 server_tmpdir=None,
                 init_list_length=1024):
        """
        Initialization of this PARI gp interpreter.

        INPUT:

        - ``stacksize`` (int, default 10000000) -- the initial PARI
          stacksize in bytes (default 10MB)
        - ``maxread`` (int, default 100000) -- ??
        - ``script_subdirectory`` (string, default None) -- name of the subdirectory of SAGE_EXTCODE/pari from which to read scripts
        - ``logfile`` (string, default None) -- log file for the pexpect interface
        - ``server`` -- name of remote server
        - ``server_tmpdir`` -- name of temporary directory on remote server
        - ``init_list_length`` (int, default 1024) -- length of initial list of local variables.

        EXAMPLES::

            sage: gp == loads(dumps(gp))
            True
        """
        Expect.__init__(self,
                        name = 'pari',
                        prompt = '\\? ',
                        command = "gp --emacs --quiet --stacksize %s"%stacksize,
                        maxread = maxread,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        logfile=logfile,
                        eval_using_file_cutoff=1024)
        self.__seq = 0
        self.__var_store_len = 0
        self.__init_list_length = init_list_length


    def _repr_(self):
        """
        String representation of this PARI gp interpreter.

        EXAMPLES::

            sage: gp # indirect doctest
            PARI/GP interpreter
        """
        return 'PARI/GP interpreter'

    def __reduce__(self):
        """
        EXAMPLES::

            sage: gp.__reduce__()
            (<function reduce_load_GP at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            PARI/GP interpreter
        """
        return reduce_load_GP, tuple([])

    def _function_class(self):
        """
        Returns the GpFunction class.

        EXAMPLES::

            sage: gp._function_class()
            <class 'sage.interfaces.gp.GpFunction'>
            sage: type(gp.gcd)
            <class 'sage.interfaces.gp.GpFunction'>
        """
        return GpFunction

    def _quit_string(self):
        """
        Returns the string used to quit the GP interpreter.

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
        Returns the string used to read filename into GP.

        EXAMPLES::

            sage: gp._read_in_file_command('test')
            'read("test")'

        ::

            sage: filename = tmp_filename()
            sage: f = open(filename, 'w')
            sage: f.write('x = 22;\n')
            sage: f.close()
            sage: gp.read(filename)
            sage: gp.get('x').strip()
            '22'
        """
        return 'read("%s")'%filename

    def trait_names(self):
        """
        EXAMPLES::

            sage: c = gp.trait_names()
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

    def set_precision(self, prec=None):
        """
        Sets the PARI precision (in decimal digits) for real computations, and returns the old value.

        EXAMPLES::

            sage: old_prec = gp.set_precision(53); old_prec
            28              # 32-bit
            38              # 64-bit
            sage: gp.get_precision()
            53
            sage: gp.set_precision(old_prec)
            53
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
        Sets the PARI power series precision, and returns the old precision.

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

        We verify that trac 11617 is fixed::

            sage: gp._eval_line('a='+str(range(2*10^5)))[:70]
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
                raise RuntimeError, a
            return self._eval_line(line, allow_use_file=allow_use_file,
                                   wait_for_prompt=wait_for_prompt)
        else:
            return a

    def cputime(self, t=None):
        """
        cputime for pari - cputime since the pari process was started.

        INPUT:


        -  ``t`` - (default: None); if not None, then returns
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

    def set_default(self, var=None, value=None):
        """
        Set a PARI gp configuration variable, and return the old value.

        INPUT:

        - ``var`` (string, default None) -- the name of a PARI gp
          configuration variable.  (See ``gp.default()`` for a list.)
        - ``value`` -- the value to set the variable to.

        EXAMPLES::

            sage: old_prec = gp.set_default('realprecision',100); old_prec
            28              # 32-bit
            38              # 64-bit
            sage: gp.get_default('realprecision')
            100
            sage: gp.set_default('realprecision',old_prec)
            100
            sage: gp.get_default('realprecision')
            28              # 32-bit
            38              # 64-bit
        """
        old = self.get_default(var)
        self._eval_line('default(%s,%s)'%(var,value))
        return old

    def get_default(self, var=None):
        """
        Return the current value of a PARI gp configuration variable.

        INPUT:

        - ``var`` (string, default None) -- the name of a PARI gp
          configuration variable.  (See ``gp.default()`` for a list.)

        OUTPUT:

        (string) the value of the variable.

        EXAMPLES::

            sage: gp.get_default('log')
            0
            sage: gp.get_default('datadir')
            '.../local/share/pari'
            sage: gp.get_default('seriesprecision')
            16
            sage: gp.get_default('realprecision')
            28              # 32-bit
            38              # 64-bit
        """
        return eval(self._eval_line('default(%s)'%var))

    def set(self, var, value):
        """
        Set the GP variable var to the given value.

        INPUT:

        - ``var`` (string) -- a valid GP variable identifier
        - ``value`` -- a value for the variable

        EXAMPLES::

            sage: gp.set('x', '2')
            sage: gp.get('x')
            '2'
        """
        cmd = '%s=%s;'%(var,value)
        out = self.eval(cmd)
        if out.find('***') != -1:
            raise TypeError, "Error executing code in GP:\nCODE:\n\t%s\nPARI/GP ERROR:\n%s"%(cmd, out)


    def get(self, var):
        """
        Get the value of the GP variable var.

        INPUT:

        - ``var`` (string) -- a valid GP variable identifier

        EXAMPLES::

            sage: gp.set('x', '2')
            sage: gp.get('x')
            '2'
        """
        return self.eval('print(%s)'%var)

    def kill(self, var):
        """
        Kill the value of the GP variable var.

        INPUT:

        - ``var`` (string) -- a valid GP variable identifier

        EXAMPLES::

            sage: gp.set('xx', '22')
            sage: gp.get('xx')
            '22'
            sage: gp.kill('xx')
            sage: gp.get('xx')
            'xx'
        """
        self.eval('kill(%s)'%var)

    #def xclear(self, var):
        #"""
        #Clear the variable named var.
        #"""
        #for varname based memory -- only 65000 variables and then dead.
        #self.eval('kill(%s)'%var)
        # for array-based memory this is best:
        #self.eval('%s=0'%var)
        # However, I've commented it out, since PARI doesn't seem
        # to ever free any memory on its stack anyways.
        # Killing variables as above takes a lot of time in some
        # cases, also.

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
        """
        self.__seq += 1
        if self.__seq >= self.__var_store_len:
            if self.__var_store_len == 0:
                self.eval('sage=vector(%s,k,0);'%self.__init_list_length)
                self.__var_store_len = self.__init_list_length
            else:
                self.eval('sage=concat(sage, vector(%s,k,0));'%self.__var_store_len)
                self.__var_store_len *= 2
                verbose("doubling PARI/sage object vector: %s"%self.__var_store_len)
        return 'sage[%s]'%self.__seq

    def quit(self, verbose=False, timeout=0.25):
        """
        Terminate the GP process.

        EXAMPLES::

            sage: a = gp('10'); a
            10
            sage: gp.quit()
            sage: a
            Traceback (most recent call last):
            ...
            ValueError: The pari session in which this object was defined is no longer running.
            sage: gp(pi)
            3.1415926535897932384626433832795028842    # 64-bit
            3.141592653589793238462643383              # 32-bit
        """
        self.__var_store_len = 0
        Expect.quit(self, verbose=verbose, timeout=timeout)

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
        Returns the version of GP being used.

        EXAMPLES::

            sage: gp.version()  # not tested
            ((2, 4, 3), 'GP/PARI CALCULATOR Version 2.4.3 (development svn-12577)')
        """
        return gp_version()

    def _object_class(self):
        """
        Returns the GpElement class.

        EXAMPLES::

            sage: gp._object_class()
            <class 'sage.interfaces.gp.GpElement'>
            sage: type(gp(2))
            <class 'sage.interfaces.gp.GpElement'>
        """
        return GpElement

    def _function_element_class(self):
        """
        Returns the GpFunctionElement class.

        EXAMPLES::

            sage: gp._function_element_class()
            <class 'sage.interfaces.gp.GpFunctionElement'>

        ::

            sage: type(gp(2).gcd)
            <class 'sage.interfaces.gp.GpFunctionElement'>
        """
        return GpFunctionElement

    def _true_symbol(self):
        """
        Returns the symbol used for truth in GP.

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
        Returns the symbol used for falsity in GP.

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
        Returns the symbol used for equality in GP.

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
        Returns the symbol to denote the exponent of a number in GP.

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
        Returns GP's help for ``command``.

        EXAMPLES::

            sage: gp.help('gcd')
            'gcd(x,{y}): greatest common divisor of x and y.'
        """
        return self.eval('?%s'%command).strip()

    def new_with_bits_prec(self, s, precision = 0):
        r"""
        Creates a GP object from s with ``precision`` bits of
        precision. GP actually automatically increases this precision to
        the nearest word (i.e. the next multiple of 32 on a 32-bit machine,
        or the next multiple of 64 on a 64-bit machine).

        EXAMPLES::

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
            prec = int(precision/3.321928095)
            self.set_real_precision(prec)
            x = self(s)
            self.set_real_precision(old_prec)
        else:
            x = self(s)
        return x


class GpElement(ExpectElement):
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
        False
        sage: loads(E.dumps())
        [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351, [-1.618909932267371342378000940, -0.3155450338663143288109995302 - 2.092547096911958607981689447*I, -0.3155450338663143288109995302 + 2.092547096911958607981689447*I]~, 2.780740013766729771063197627, 1.390370006883364885531598814 - 1.068749776356193066159263548*I, 3.109648242324380328550149122 + 1.009741959000000000000000000 E-28*I, 1.554824121162190164275074561 + 1.064374745210273756943885994*I, 2.971915267817909670771647951] # 32-bit
        [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351, [-1.6189099322673713423780009396072169751, -0.31554503386631432881099953019639151248 - 2.0925470969119586079816894466366945829*I, -0.31554503386631432881099953019639151248 + 2.0925470969119586079816894466366945829*I]~, 2.7807400137667297710631976271813584994, 1.3903700068833648855315988135906792497 - 1.0687497763561930661592635474375038788*I, 3.1096482423243803285501491221965830079 + 2.3509887016445750160000000000000000000 E-38*I, 1.5548241211621901642750745610982915040 + 1.0643747452102737569438859937299427442*I, 2.9719152678179096707716479509361896060] # 64-bit
        sage: E
        [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351, [-1.618909932267371342378000940, -0.3155450338663143288109995302 - 2.092547096911958607981689447*I, -0.3155450338663143288109995302 + 2.092547096911958607981689447*I]~, 2.780740013766729771063197627, 1.390370006883364885531598814 - 1.068749776356193066159263548*I, 3.109648242324380328550149122 + 1.009741959 E-28*I, 1.554824121162190164275074561 + 1.064374745210273756943885994*I, 2.971915267817909670771647951] # 32-bit
        [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351, [-1.6189099322673713423780009396072169751, -0.31554503386631432881099953019639151248 - 2.0925470969119586079816894466366945829*I, -0.31554503386631432881099953019639151248 + 2.0925470969119586079816894466366945829*I]~, 2.7807400137667297710631976271813584994, 1.3903700068833648855315988135906792497 - 1.0687497763561930661592635474375038788*I, 3.1096482423243803285501491221965830079 + 2.350988701644575016 E-38*I, 1.5548241211621901642750745610982915040 + 1.0643747452102737569438859937299427442*I, 2.9719152678179096707716479509361896060]  # 64-bit

    The two elliptic curves look the same, but internally the floating
    point numbers are slightly different.
    """

    def _sage_(self):
        """
        Convert this GpElement into a Sage object, if possible.

        EXAMPLES::

            sage: gp(I).sage()
            i
            sage: gp(I).sage().parent()
            Maximal Order in Number Field in i with defining polynomial x^2 + 1

        ::

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
        """
        return pari(str(self)).python()

    def __long__(self):
        """
        Return Python long.

        EXAMPLES::

            sage: long(gp(10))
            10L
        """
        return long(str(self))

    def __float__(self):
        """
        Return Python float.

        EXAMPLES::

            sage: float(gp(10))
            10.0
        """
        return float(pari(str(self)))

    def bool(self):
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
        return P.eval('%s != 0'%(self.name())) == '1'

    def _complex_mpfr_field_(self, CC):
        """
        Return ComplexField element of self.

        INPUT:

        - ``CC`` -- a Complex or Real Field.

        EXAMPLES::

            sage: z = gp(1+15*I); z
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

        return CC((CC(1)*pari(self))._sage_())

    def _complex_double_(self, CDF):
        """
        Returns this value as a CDF element.

        EXAMPLES::

            sage: CDF(gp(pi+I*e))
            3.14159265359 + 2.71828182846*I
        """
        # Retrieving values from another computer algebra system is
        # slow anyway, right?
        cc_val = self._complex_mpfr_field_(sage.rings.complex_field.ComplexField())
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
    #def _latex_(self):
    #    P = self._check_valid()
    #    return P.eval('printtex(%s)'%self.name())

    def trait_names(self):
        """
        EXAMPLES::

            sage: 'gcd' in gp(2).trait_names()
            True
        """
        return self.parent().trait_names()


class GpFunctionElement(FunctionElement):
    pass

class GpFunction(ExpectFunction):
    pass



def is_GpElement(x):
    """
    Returns True of x is a GpElement.

    EXAMPLES::

        sage: from sage.interfaces.gp import is_GpElement
        sage: is_GpElement(gp(2))
        True
        sage: is_GpElement(2)
        False
    """
    return isinstance(x, GpElement)

from sage.env import DOT_SAGE, SAGE_ETC
import os

# Set GPRC environment variable to $SAGE_ETC/gprc.expect
os.environ["GPRC"] = os.path.join(SAGE_ETC, 'gprc.expect')

# An instance
gp = Gp(logfile=os.path.join(DOT_SAGE,'gp-expect.log')) # useful for debugging!

def reduce_load_GP():
    """
    Returns the GP interface object defined in sage.interfaces.gp.

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
    os.system('gp')


def gp_version():
    """
    EXAMPLES::

        sage: gp.version()  # not tested
        ((2, 4, 3), 'GP/PARI CALCULATOR Version 2.4.3 (development svn-12577)')
    """
    v = gp.eval(r'\v')
    i = v.find("Version ")
    w = v[i+len("Version "):]
    i = w.find(' ')
    w = w[:i]
    t = tuple([int(n) for n in w.split('.')])
    k = v.find('\n')
    return t, v[:k].strip()
