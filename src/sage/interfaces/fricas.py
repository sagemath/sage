r"""
Interface to FriCAS

.. TODO::

    - ``fricas(dilog(x))`` should be ``dilog(-(x-1))``, and some
      more conversions in ``sage.functions`` are missing

FriCAS is a free GPL-compatible (modified BSD license) general
purpose computer algebra system based on Axiom.  The FriCAS
website can be found at http://fricas.sourceforge.net/.

AUTHORS:

- Mike Hansen (2009-02): Split off the FriCAS interface from
  the Axiom interface.

- Martin Rubey, Bill Page (2016-08): Completely separate from Axiom,
  implement more complete translation from FriCAS to SageMath types.


EXAMPLES::

    sage: fricas('3 * 5')                                                       # optional - fricas
    15
    sage: a = fricas(3) * fricas(5); a                                          # optional - fricas
    15

The type of a is :class:`FriCASElement`, i.e., an element of the
FriCAS interpreter::

    sage: type(a)                                                               # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: a.parent()                                                            # optional - fricas
    FriCAS

The underlying FriCAS type of a is also available, via the type
method::

    sage: a.typeOf()                                                            # optional - fricas
    PositiveInteger

FriCAS objects are normally displayed using "ASCII art"::

    sage: fricas(2/3)                                                           # optional - fricas
      2
      -
      3
    sage: fricas('x^2 + 3/7')                                                   # optional - fricas
       2   3
      x  + -
           7

Functions defined in FriCAS are available as methods of the :class:`fricas<FriCAS>` object::

    sage: F = fricas.factor('x^5 - y^5'); F                                     # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    sage: type(F)                                                               # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: F.typeOf()                                                            # optional - fricas
    Factored(Polynomial(Integer))

We can also create a FriCAS polynomial and apply the function
``factor`` from FriCAS.  The notation ``f.factor()`` is consistent
with how the rest of SageMath works::

    sage: f = fricas('x^5 - y^5')                                               # optional - fricas
    sage: f^2                                                                   # optional - fricas
       10     5 5    10
      y   - 2x y  + x
    sage: f.factor()                                                            # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

For many FriCAS types, translation to an appropriate SageMath type is
available::

    sage: f.factor().sage()                                                     # optional - fricas
    (y - x) * (y^4 + y^3*x + y^2*x^2 + y*x^3 + x^4)

Control-C interruption works well with the FriCAS interface. For
example, try the following sum but with a much bigger range, and hit
control-C::

    sage:  f = fricas('(x^5 - y^5)^10000')                                      # not tested - fricas
    Interrupting FriCAS...
    ...
    KeyboardInterrupt: Ctrl-c pressed while running FriCAS

Let us demonstrate some features of FriCAS.  FriCAS can guess a
differential equation for the generating function for integer
partitions::

    sage: fricas("guessADE([partition n for n in 0..40], homogeneous==4)")      # optional - fricas
     [
       [
           n
         [x ]f(x):
              2    3 (iv)          2    2 ,             3  ,,,         2    2 ,,   2
             x f(x) f    (x) + (20x f(x) f (x) + 5x f(x) )f   (x) - 39x f(x) f  (x)
    <BLANKLINE>
           +
                 2     ,   2           2 ,           3  ,,        2 ,   4
             (12x f(x)f (x)  - 15x f(x) f (x) + 4f(x) )f  (x) + 6x f (x)
    <BLANKLINE>
           +
                      ,   3         2 ,   2
             10x f(x)f (x)  - 16f(x) f (x)
    <BLANKLINE>
             =
             0
         ,
                        2     3      4
        f(x)= 1 + x + 2x  + 3x  + O(x )]
       ]

FriCAS can solve linear ordinary differential equations::

    sage: fricas.set("y", "operator y")                                         # optional - fricas
    sage: fricas.set("deq", "x^3*D(y x, x, 3) + x^2*D(y x, x, 2) - 2*x*D(y x, x) + 2*y x - 2*x^4")  # optional - fricas
    sage: fricas.set("sol", "solve(deq, y, x)"); fricas("sol")                  # optional - fricas
                  5      3      2               3     2      3      3     2
                 x  - 10x  + 20x  + 4         2x  - 3x  + 1 x  - 1 x  - 3x  - 1
    [particular= --------------------,basis= [-------------,------,------------]]
                          15x                       x          x         x

    sage: fricas("sol.particular").sage()                                       # optional - fricas
    1/15*(x^5 - 10*x^3 + 20*x^2 + 4)/x
    sage: fricas("sol.basis").sage()                                            # optional - fricas
    [(2*x^3 - 3*x^2 + 1)/x, (x^3 - 1)/x, (x^3 - 3*x^2 - 1)/x]
    sage: fricas.eval(")clear values y deq sol")                                # optional - fricas
    ''

FriCAS can expand expressions into series::

    sage: x = var('x'); ex = sqrt(cos(x)); a = fricas(ex).series(x=0); a        # optional - fricas
        1  2    1  4    19   6     559   8     29161    10      11
    1 - - x  - -- x  - ---- x  - ------ x  - --------- x   + O(x  )
        4      96      5760      645120      116121600

    sage: a.coefficients()[38].sage()                                           # optional - fricas
    -29472026335337227150423659490832640468979/274214482066329363682430667508979749984665600000000

    sage: ex = sqrt(atan(x)); a = fricas(ex).series(x=0); a                     # optional - fricas
     1      5        9
     -      -        -
     2   1  2    31  2      6
    x  - - x  + --- x  + O(x )
         6      360

    sage: a.coefficient(9/2).sage()                                             # optional - fricas
    31/360

    sage: x = fricas("x::TaylorSeries Fraction Integer")                        # optional - fricas
    sage: y = fricas("y::TaylorSeries Fraction Integer")                        # optional - fricas
    sage: 2*(1+2*x+sqrt(1-4*x)-2*x*y).recip()                                   # optional - fricas
                   2      3     2 2     3      4       4       5
       1 + (x y + x ) + 2x  + (x y  + 2x y + 6x ) + (4x y + 18x )
     +
         3 3     4 2      5       6       5 2      6        7
       (x y  + 3x y  + 13x y + 57x ) + (6x y  + 40x y + 186x )
     +
         4 4     5 3      6 2       7        8       6 3      7 2       8         9
       (x y  + 4x y  + 21x y  + 130x y + 622x ) + (8x y  + 66x y  + 432x y + 2120x )
     +
         5 5     6 4      7 3       8 2        9         10
       (x y  + 5x y  + 30x y  + 220x y  + 1466x y + 7338x  ) + O(11)

FriCAS does some limits right::

    sage: x = var('x'); ex = x^2*exp(-x)*Ei(x) - x; fricas(ex).limit(x=oo)      # optional - fricas
    1

"""

###########################################################################
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#                     2007 Bill Page
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
###########################################################################
from __future__ import print_function
# from __future__ import absolute_import

from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.interfaces.expect import Expect, ExpectElement, FunctionElement, ExpectFunction
from sage.misc.misc import SAGE_TMP_INTERFACE
from sage.env import DOT_SAGE
import re
import six

FRICAS_SINGLE_LINE_START = 3 # where the output starts when it fits next to the line number
FRICAS_MULTI_LINE_START = 2  # and when it doesn't
FRICAS_LINE_LENGTH = 80      # length of a line, should match the line length in sage
# the following messages have, unfortunately, no markup.
FRICAS_WHAT_OPERATIONS_STRING = "Operations whose names satisfy the above pattern\(s\):"
FRICAS_ERROR_IN_LIBRARY_CODE = ">> Error detected within library code:"

# only the last command should be necessary to make the interface
# work, the other are optimizations.  Beware that lisp distinguishes
# between ' and ".
FRICAS_INIT_CODE = (
")set functions compile on",
")set message autoload off",
")set message type off",
")set output length " + str(FRICAS_LINE_LENGTH),
")lisp (setf |$ioHook|"
"            (lambda (x &optional args)"
"              (when (member x '(|startAlgebraOutput| |endOfAlgebraOutput|"
"                                |startKeyedMsg|      |endOfKeyedMsg|))"
"               (prin1 x)"
"               (princ #\\Newline))))")

FRICAS_LINENUMBER_OFF_CODE = ")lisp (setf |$IOindex| NIL)"
FRICAS_FIRST_PROMPT = "\(1\) -> "
FRICAS_LINENUMBER_OFF_PROMPT = "\(NIL\) -> "

class FriCAS(ExtraTabCompletion, Expect):
    """
    Interface to a FriCAS interpreter.
    """
    def __init__(self, name='fricas', command='fricas -nox -noclef',
                 script_subdirectory=None, logfile=None,
                 server=None, server_tmpdir=None):
        """
        Create an instance of the FriCAS interpreter.

        TESTS::

            sage: fricas == loads(dumps(fricas))                                # optional - fricas
            True
        """
        eval_using_file_cutoff = 4096-5 # magic number from Expect._eval_line (there might be a bug)
        assert max(len(c) for c in FRICAS_INIT_CODE) < eval_using_file_cutoff
        self.__eval_using_file_cutoff = eval_using_file_cutoff
        self._COMMANDS_CACHE = '%s/%s_commandlist_cache.sobj'%(DOT_SAGE, name)
        # we run the init code in _start to avoid spurious output
        Expect.__init__(self,
                        name = name,
                        prompt = FRICAS_FIRST_PROMPT,
                        command = command,
                        script_subdirectory = script_subdirectory,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        init_code = [],
                        logfile = logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)

    def _start(self):
        """
        Start the FriCAS interpreter and switch off the linenumbers.

        EXAMPLES::

            sage: a = FriCAS()                                                  # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False
            sage: a._start()                                                    # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            True
            sage: a.quit()                                                      # optional - fricas
        """
        # setting the prompt properly is necessary for restarting FriCAS
        self._prompt = FRICAS_FIRST_PROMPT
        Expect._start(self)
        for line in FRICAS_INIT_CODE:
            self.eval(line, reformat=False)
        # switching off the line numbers also modified the prompt
        self._prompt = FRICAS_LINENUMBER_OFF_PROMPT
        self.eval(FRICAS_LINENUMBER_OFF_CODE, reformat=False)

    def _quit_string(self):
        """
        Returns the string used to quit FriCAS.

        EXAMPLES::

            sage: fricas._quit_string()                                         # optional - fricas
            ')quit\r'
            sage: a = FriCAS()                                                  # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False
            sage: a._start()                                                    # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            True
            sage: a.quit()                                                      # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False

        TESTS::

            sage: import psutil                                                 # optional - fricas
            sage: p = fricas.pid(); pr = psutil.Process(p); pr                  # optional - fricas
            <psutil.Process(pid=..., name='sman') at ...>
            sage: pr.children()                                                 # optional - fricas
            [<psutil.Process(pid=..., name='AXIOMsys') at ...>,
             <psutil.Process(pid=..., name='session') at ...>,
             <psutil.Process(pid=..., name='spadclient') at ...>,
             <psutil.Process(pid=..., name='sman') at ...>]
            sage: fricas.quit()                                                 # optional - fricas
            sage: pr.is_running()                                               # optional - fricas, random
            False
        """
        return ')quit\r'

    def _commands(self):
        """
        Returns a list of commands available. This is done by parsing the
        result of the first section of the output of ')what things'.

        EXAMPLES::

            sage: cmds = fricas._commands()                                     # optional - fricas
            sage: len(cmds) > 100                                               # optional - fricas
            True
            sage: '<' in cmds                                                   # optional - fricas
            True
            sage: 'factor' in cmds                                              # optional - fricas
            True
        """
        output = self.eval(")what operations", reformat=False)
        m = re.search(FRICAS_WHAT_OPERATIONS_STRING + "\r\n(.*)\r\n\|startKeyedMsg\|", output, flags = re.DOTALL)
        l = m.groups()[0].split()
        return l

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        """
        Returns a list of all the commands defined in Fricas and optionally
        (per default) store them to disk.

        EXAMPLES::

            sage: c = fricas._tab_completion(use_disk_cache=False, verbose=False)         # optional - fricas
            sage: len(c) > 100                                                  # optional - fricas
            True
            sage: 'factor' in c                                                 # optional - fricas
            True
            sage: '**' in c                                                     # optional - fricas
            False
            sage: 'upperCase?' in c                                             # optional - fricas
            False
            sage: 'upperCase_q' in c                                            # optional - fricas
            True
            sage: 'upperCase_e' in c                                            # optional - fricas
            True
        """
        try:
            return self.__tab_completion
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__tab_completion = sage.misc.persist.load(self._COMMANDS_CACHE)
                    return self.__tab_completion
                except IOError:
                    pass
            if verbose:
                print("\nBuilding %s command completion list (this takes" % self)
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % self._COMMANDS_CACHE)
            v = self._commands()

            #Process we now need process the commands to strip out things which
            #are not valid Python identifiers.
            valid = re.compile('[^a-zA-Z0-9_]+')
            names = [x for x in v if valid.search(x) is None]

            #Change everything that ends with ? to _q and
            #everything that ends with ! to _e
            names += [x[:-1]+"_q" for x in v if x.endswith("?")]
            names += [x[:-1]+"_e" for x in v if x.endswith("!")]

            self.__tab_completion = names
            if len(v) > 200:
                # Fricas is actually installed.
                sage.misc.persist.save(v, self._COMMANDS_CACHE)
            return names

    def _read_in_file_command(self, filename):
        """
        Return the FriCAS command to read the file ``filename``.

        INPUT:

        - ``filename``, a string ending in '.input'.

        OUTPUT:

        - a string with the command for reading filename without output.

        TESTS:

        Evaluate a rather long line::

            sage: len(fricas([i for i in range(600)]))                          # optional - fricas, indirect doctest
            600

        """
        if not filename.endswith('.input'):
            raise ValueError("the filename must end with .input")

        return ')read %s )quiet'%filename

    def _local_tmpfile(self):
        """
        Return a local tmpfile ending with ".input" used to buffer long
        command lines sent to FriCAS.

        """
        try:
            return self.__local_tmpfile
        except AttributeError:
            self.__local_tmpfile = os.path.join(SAGE_TMP_INTERFACE, 'tmp' + str(self.pid()) + '.input')
            return self.__local_tmpfile

    def _remote_tmpfile(self):
        """
        Return a remote tmpfile ending with ".input" used to buffer long
        command lines sent to FriCAS.

        """
        try:
            return self.__remote_tmpfile
        except AttributeError:
            self.__remote_tmpfile = self._remote_tmpdir()+"/interface_%s:%s.input"%(LOCAL_IDENTIFIER,self.pid())
            return self.__remote_tmpfile

# what I expect from FriCAS:

# 1.) in set(self, var, value)
#
# no markers:
# there could be some "debugging" output, as in fricas("guessADE([1,1,1,1], debug==true)")
#
# startKeyedMsg: an error happened
#
# 2.) in get(self, var)
#
# |startAlgebraOutput\|...|endOfAlgebraOutput\|
#
# 3.) I also need a routine to send a system command and get its output.

    def _check_errors(self, line, output):
        """
        Check whether output contains an error and, if so, raise it.

        INPUT:

        - ``line``, a string that was sent to FriCAS.

        - ``output``, a string returned by FriCAS

        OUTPUT:

        None

        TESTS::

            sage: fricas.set("x", "[i fo83r i in 0..17]")                       # optional - fricas, indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: An error occurred when FriCAS evaluated '[i fo83r i in 0..17]':
              Line   1: x:=[i fo83r i in 0..17];
                       ...A..........B
              Error  A: Missing mate.
              Error  B: syntax error at top level
              Error  B: Possibly missing a ]
               3 error(s) parsing

            sage: fricas.set("x", "something stupid")                           # optional - fricas, indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: An error occurred when FriCAS evaluated 'something stupid':
               There are no library operations named something
                  Use HyperDoc Browse or issue
                                           )what op something
                  to learn if there is any operation containing " something " in its
                  name.
            <BLANKLINE>
               Cannot find a definition or applicable library operation named
                  something with argument type(s)
                                            Variable(stupid)
            <BLANKLINE>
                  Perhaps you should use "@" to indicate the required return type, or
                  "$" to specify which version of the function you need.

        """
        # otherwise there might be a message
        m = re.search("\|startKeyedMsg\|\r\n(.*)\r\n\|endOfKeyedMsg\|\r", output, flags = re.DOTALL)
        if m:
            replacements = [('|startKeyedMsg|\r\n', ''),
                            ('|endOfKeyedMsg|\r', '')]
            for old, new in replacements:
                output = output.replace(old, new)
            raise RuntimeError("An error occurred when FriCAS evaluated '%s':\n%s" % (line, output))

        # or even an error
        if FRICAS_ERROR_IN_LIBRARY_CODE in output:
            raise RuntimeError("An error occurred when FriCAS evaluated '%s':\n%s" % (line, output))

    def set(self, var, value):
        """
        Set a variable to a value in FriCAS.

        INPUT:

        - ``var``, ``value``: strings, the first representing a valid
          FriCAS variable identifier, the second a FriCAS expression.

        OUTPUT: None

        EXAMPLES::

            sage: fricas.set('xx', '2')                                         # optional - fricas
            sage: fricas.get('xx')                                              # optional - fricas
            '2'

        """
        cmd = '%s%s%s;'%(var,self._assign_symbol(), value)
        output = self.eval(cmd, reformat=False)
        self._check_errors(value, output)

    def get(self, var):
        r"""
 Get the string representation of the value (more precisely, the
        OutputForm) of a variable or expression in FriCAS.

        If FriCAS cannot evaluate `var` an error is raised.

        EXAMPLES::

            sage: fricas.set('xx', '2')                                         # optional - fricas
            sage: fricas.get('xx')                                              # optional - fricas
            '2'
            sage: a = fricas('(1 + sqrt(2))^5')                                 # optional - fricas
            sage: fricas.get(a.name())                                          # optional - fricas
            '   +-+\r\n29\\|2  + 41'
            sage: fricas.get('(1 + sqrt(2))^5')                                 # optional - fricas
            '   +-+\r\n29\\|2  + 41'
            sage: fricas.new('(1 + sqrt(2))^5')                                 # optional - fricas
               +-+
            29\|2  + 41
        """
        output = self.eval(str(var), reformat=False)
        # if there is AlgebraOutput we ask no more
        m = re.search("\|startAlgebraOutput\|\r\n(.*)\r\n\|endOfAlgebraOutput\|\r", output, flags = re.DOTALL)
        if m:
            lines = m.groups()[0].split("\r\n")
            if max(len(line) for line in lines) < FRICAS_LINE_LENGTH:
                return "\r\n".join(line[FRICAS_SINGLE_LINE_START:] for line in lines)
            else:
                return "\r\n".join(line[FRICAS_MULTI_LINE_START:] for line in lines)

        self._check_errors(var, output)

    def get_string(self, var):
        """
        Return the value of a FriCAS string as a string, without checking
        that it is a string.

        TESTS:

        We test that strings are returned properly::

            sage: r = fricas.get_string('concat([concat(string(i)," ") for i in 0..299])')   # optional - fricas
            sage: r == " ".join([str(i) for i in range(300)]) + ' '                          # optional - fricas
            True

            sage: fricas.get_string('concat([string(1) for i in 1..5])') == "1"*5            # optional - fricas
            True

            sage: fricas.get_string('concat([string(1) for i in 1..10000])') == "1"*10000    # optional - fricas
            True

        """
        return self.get(str(var)).replace("\r\n", "")[1:-1]

    def get_integer(self, var):
        """
        Return the value of a FriCAS integer as an integer, without
        checking that it is an integer.

        TESTS::

            sage: fricas.get_integer('factorial 1111') == factorial(1111)       # optional - fricas
            True

        """
        return int(self.get_unparsed_InputForm(str(var)))

    def get_boolean(self, var):
        """
        Return the value of a FriCAS boolean as a boolean, without checking
        that it is a boolean.

        TESTS::

            sage: fricas.get_boolean('(1=1)::Boolean') == True                  # optional - fricas
            True

            sage: fricas.get_boolean('(1=2)::Boolean') == False                 # optional - fricas
            True
        """
        return self.get(str(var)).replace("\r\n", "") == "true"

    def get_unparsed_InputForm(self, var):
        """
        Return the unparsed ``InputForm`` as a string.

        .. TODO::

            - catch errors, especially when InputForm is not available:

                - for example when integration returns ``"failed"``

                - ``UnivariatePolynomial``

            - should we provide workarounds, too?

        TESTS::

            sage: fricas.get_unparsed_InputForm('1..3')                         # optional - fricas
            '1..3$Segment(Integer())'

        """
        return self.get_string('unparse((%s)::InputForm)' %str(var))

    def _assign_symbol(self):
        """
        Return the symbol used for setting a variable in FriCAS.

        EXAMPLES::

            sage: fricas.set("x", "1");                                         # optional - fricas, indirect doctest
            sage: fricas.get("x")                                               # optional - fricas
            '1'
            sage: fricas.eval(")cl val x")                                      # optional - fricas
            ''
        """
        return ":="

    def _equality_symbol(self):
        """
        Return the equality testing logical symbol in FriCAS.

        EXAMPLES::

            sage: a = fricas(x==6); a                                           # optional - fricas, indirect doctest
            x= 6

        A warning:

            sage: fricas.set("x", 2);                                           # optional - fricas
            sage: a = fricas(x==6); a                                           # optional - fricas
            2= 6
            sage: fricas.eval(")cl val x")                                      # optional - fricas
            ''
        """
        return "="

    def _true_symbol(self):
        """
        Return the string used for True in FriCAS.

        EXAMPLES::

            sage: str(fricas("(1=1)@Boolean")) == fricas._true_symbol()         # optional - fricas
            True
        """
        return "true"

    def _false_symbol(self):
        """
        Return the string used for False in FriCAS.

        EXAMPLES::

            sage: str(fricas("(1~=1)@Boolean")) == fricas._false_symbol()       # optional - fricas
            True
        """
        return "false"

    def _inequality_symbol(self):
        """
        Return the string used for False in FriCAS.

        EXAMPLES::

            sage: fricas(x!=0)                                                  # optional - fricas, indirect doctest
            true
        """
        return '~='

    def __repr__(self):
        """
        EXAMPLES::

            sage: fricas                                                        # optional - fricas
            FriCAS
        """
        return "FriCAS"

    def __reduce__(self):
        """
        EXAMPLES::

            sage: fricas.__reduce__()                                           # optional - fricas
            (<function reduce_load_fricas at 0x...>, ())
            sage: f, args = _                                                   # optional - fricas
            sage: f(*args)                                                      # optional - fricas
            FriCAS
        """
        return reduce_load_fricas, tuple([])

    def eval(self, code, strip=True, synchronize=False, locals=None, allow_use_file=True,
             split_lines="nofile", reformat=True, **kwds):
        """
        Evaluate ``code`` using FriCAS.

        Except ``reformat``, all arguments are passed to
        :meth:`sage.interfaces.expect.Expect.eval`.

        INPUT:

        - ``reformat`` -- bool; remove the output markers when True.

        This can also be used to pass system commands to FriCAS.

        EXAMPLES::

            sage: fricas.set("x", "1783"); fricas("x")                               # optional - fricas
            1783
            sage: fricas.eval(")cl val x");                                          # optional - fricas
            ''
            sage: fricas("x")                                                        # optional - fricas
            x

        """
        output = Expect.eval(self, code, strip=strip,
                             synchronize=synchronize, locals=locals,
                             allow_use_file=allow_use_file, split_lines=split_lines,
                             **kwds)
        if reformat:
            replacements = [('|startAlgebraOutput|\r\n', ''),
                            ('|endOfAlgebraOutput|\r', ''),
                            ('|startKeyedMsg|\r\n', ''),
                            ('|endOfKeyedMsg|\r', '')]
            for old, new in replacements:
                output = output.replace(old, new)

        return output


    def _function_class(self):
        """
        Return the FriCASExpectFunction class.

        EXAMPLES::

            sage: fricas._function_class()                                      # optional - fricas
            <class 'sage.interfaces.fricas.FriCASExpectFunction'>
            sage: type(fricas.gcd)                                              # optional - fricas
            <class 'sage.interfaces.fricas.FriCASExpectFunction'>
        """
        return FriCASExpectFunction

    def _object_class(self):
        """
        EXAMPLES::

            sage: fricas._object_class()                                        # optional - fricas
            <class 'sage.interfaces.fricas.FriCASElement'>
            sage: type(fricas(2))                                               # optional - fricas
            <class 'sage.interfaces.fricas.FriCASElement'>
        """
        return FriCASElement

    def _function_element_class(self):
        """
        Returns the FriCAS function element class.

        EXAMPLES::

            sage: fricas._function_element_class()                              # optional - fricas
            <class 'sage.interfaces.fricas.FriCASFunctionElement'>
            sage: type(fricas(2).gcd)                                           # optional - fricas
            <class 'sage.interfaces.fricas.FriCASFunctionElement'>
        """
        return FriCASFunctionElement

    def console(self):
        """
        Spawn a new FriCAS command-line session.

        EXAMPLES::

            sage: fricas.console()                                              # not tested
                             FriCAS (AXIOM fork) Computer Algebra System
                                    Version: FriCAS 1.0.5
                     Timestamp: Thursday February 19, 2009 at 06:57:33
            -----------------------------------------------------------------------------
               Issue )copyright to view copyright notices.
               Issue )summary for a summary of useful system commands.
               Issue )quit to leave AXIOM and return to shell.
            -----------------------------------------------------------------------------
        """
        fricas_console()

class FriCASElement(ExpectElement):
    """
    Instances of this class represent objects in FriCAS.

    Using the method :meth:`sage` we can translate some of them to
    SageMath objects:

    .. automethod:: _sage_
    """
    def __len__(self):
        """
        Return the length of a list.

        EXAMPLES::

            sage: v = fricas('[x^i for i in 0..5]')                             # optional - fricas
            sage: len(v)                                                        # optional - fricas
            6
        """
        P = self._check_valid()
        l = P('#(%s)' %self._name)
        return l.sage()

    def __getitem__(self, n):
        """
        We implement the sage conventions here, translating to 0-based iterables.

        We do not check validity, since many objects in FriCAS are
        iterable, in particular Streams

        .. TODO::

            - can we somehow implement negative arguments?

        TEST:

            sage: fricas("[1,2,3]")[0]                                          # optional - fricas
            1

            sage: fricas("[1,2,3]")[3]                                          # optional - fricas
            Traceback (most recent call last):
            ...
            TypeError: An error occurred when FriCAS evaluated 'elt(...,...)':
            <BLANKLINE>
            >> Error detected within library code:
            index out of range
        """
        n = int(n)
        if n < 0:
            raise IndexError("index out of range")
        P = self._check_valid()
        # use "elt" instead of "." here because then the error
        # message is clearer
        return P.new("elt(%s,%s)" %(self._name, n+1))

    def __int__(self):
        """
        TEST::

            sage: int(fricas(2))                                                # optional - fricas
            2
        """
        return int(self.sage())

    def bool(self):
        """
        Coerce the expression into a boolean.

        EXAMPLES::

            sage: fricas("1=1").bool()                                          # optional - fricas
            True
            sage: fricas("1~=1").bool()                                         # optional - fricas
            False
        """
        P = self._check_valid()
        return P.new(self._name + "::Boolean").sage()

    def __bool__(self):
        """
        Check whether the expression is different from zero.

        EXAMPLES::

            sage: fricas(0).is_zero()                                           # optional - fricas, indirect doctest
            True
        """
        P = self._check_valid()
        return not P.new("zero?(%s)" % self._name).sage()

    __nonzero__ = __bool__

    def __long__(self):
        """
        TEST::

            sage: long(fricas('1'))                                             # optional - fricas
            1L
        """
        return long(self.sage())

    def __float__(self):
        """
        TEST::

            sage: float(fricas(2))                                              # optional - fricas
            2.0
        """
        return float(self.sage())

    def _integer_(self, ZZ=None):
        """
        EXAMPLES::

            sage: ZZ(fricas('1'))                                               # optional - fricas
            1
        """
        from sage.rings.all import ZZ
        return ZZ(self.sage())

    def _rational_(self):
        """
        EXAMPLES::

            sage: QQ(fricas('-1/2'))                                            # optional - fricas
            -1/2
        """
        from sage.rings.all import QQ
        return QQ(self.sage())

    def gen(self, n):
        """
        Return an error, since the n-th generator in FriCAS is not well defined.
        """
        raise NotImplementedError

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(fricas("sin(x+y)/exp(z)*log(1+%e)"))                    # optional - fricas
            {{\log \left( {{e+1}} \right)} \  {\sin \left( {{y+x}} \right)}} \over {{e} ^{z}}

            sage: latex(fricas("matrix([[1,2],[3,4]])"))                        # optional - fricas
            \left[ \begin{array}{cc} 1 & 2 \\ 3 & 4 \end{array}  \right]

            sage: latex(fricas("integrate(sin(x+1/x),x)"))                      # optional - fricas
            \int ^{\displaystyle x} {{\sin \left( {{{{{ \%A} ^{2}}+1} \over  \%A}} \right)} \  {d \%A}}
        """
        replacements = [('\sp ', '^'),
                        ('\sp{', '^{'),
                        ('\sb ', '_'),
                        ('\sb{', '_{')]
        P = self._check_valid()
        s = P.get_string("first tex(%s)" %self._name)
        for old, new in replacements:
            s = s.replace(old, new)
        return s

    def _get_sage_type(self, domain):
        """
        INPUT:

        - ``domain``, a FriCAS SExpression

        OUTPUT:

        - a corresponding Sage type

        EXAMPLES::

            sage: m = fricas("dom(1/2)::Any")                                   # optional - fricas
            sage: fricas(0)._get_sage_type(m)                                   # optional - fricas
            Rational Field
        """
        from sage.rings.all import ZZ, QQ, QQbar, PolynomialRing, RDF
        from sage.rings.fraction_field import FractionField
        from sage.rings.finite_rings.integer_mod_ring import Integers
        from sage.rings.real_mpfr import RealField
        from sage.symbolic.ring import SR
        from sage.matrix.constructor import matrix

        # first implement domains without arguments
        head = str(domain.car())
        if head in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ
        if head == "String":
            return str
        if head == "Float":
            P = self._check_valid()
            prec = max(P.new("length mantissa(%s)" %self._name).sage(), 53)
            return RealField(prec)
        if head == "DoubleFloat":
            return RDF
        if head == "AlgebraicNumber":
            return QQbar

        # now implement "functorial" types
        if head == "OrderedCompletion":
            # this is a workaround, I don't know how translate this
            return SR

        if head == "IntegerMod":
            return Integers(domain[1].integer().sage())

        if head == "Fraction":
            return FractionField(self._get_sage_type(domain[1]))

        if head == "Expression":
            return SR

        if head == "Polynomial":
            # this is a workaround, since in sage we always have to specify the variables
            return SR

        raise NotImplementedError("The translation of FriCAS type %s to sage is not yet implemented." %domain)

    def _sage_expression(self, unparsed_InputForm):
        """
        Convert an expression to an element of the Symbolic Ring.

        This does not depend on `self`.  Instead, for practical
        reasons of the implementation of `self._sage_`, it takes the
        unparsed InputForm as argument.

        .. TODO::

             We really should walk through the InputForm here.

        TESTS::

            sage: f = fricas('integrate(sin(x^2), x)'); f                       # optional - fricas
                       +---+
                       | 2
            fresnelS(x |--- )
                      \|%pi
            -----------------
                   +---+
                   | 2
                   |---
                  \|%pi
            sage: s = fricas.get_unparsed_InputForm(f._name); s                 # optional - fricas
            'fresnelS(x*(2/pi())^(1/2))/((2/pi())^(1/2))'
            sage: f._sage_expression(s)                                         # optional - fricas
            1/2*sqrt(2)*sqrt(pi)*fresnelS(sqrt(2)*x/sqrt(pi))

        """
        from sage.symbolic.ring import SR
        s = unparsed_InputForm
        replacements = [('pi()', 'pi '),
                        ('::Symbol', ' ')]
        for old, new in replacements:
            s = s.replace(old, new)
        try:
            return SR(s)
        except TypeError:
            raise NotImplementedError("The translation of the FriCAS Expression %s to sage is not yet implemented." %s)


    def _sage_(self):
        """
        Convert self to a Sage object.

        EXAMPLES:

        Floats::

            sage: fricas(2.1234).sage()                                         # optional - fricas
            2.12340000000000
            sage: _.parent()                                                    # optional - fricas
            Real Field with 53 bits of precision
            sage: a = RealField(100)(pi)                                        # optional - fricas
            sage: fricas(a).sage()                                              # optional - fricas
            3.1415926535897932384626433833
            sage: _.parent()                                                    # optional - fricas
            Real Field with 100 bits of precision
            sage: fricas(a).sage() == a                                         # optional - fricas
            True
            sage: fricas(2.0).sage()                                            # optional - fricas
            2.00000000000000
            sage: _.parent()                                                    # optional - fricas
            Real Field with 53 bits of precision

        Algebraic numbers::

            sage: a = fricas('(1 + sqrt(2))^5'); a                              # optional - fricas
               +-+
            29\|2  + 41
            sage: b = a.sage(); b                                               # optional - fricas
            82.0121933088198?
            sage: b.radical_expression()                                        # optional - fricas
            29*sqrt(2) + 41

        Integers modulo n::

            sage: fricas("((42^17)^1783)::IntegerMod(5^(5^5))").sage() == Integers(5^(5^5))((42^17)^1783) # optional - fricas
            True

        We can also convert FriCAS's polynomials to Sage polynomials::

            sage: a = fricas(x^2 + 1); a.typeOf()                               # optional - fricas
            Polynomial(Integer)
            sage: a.sage()                                                      # optional - fricas
            x^2 + 1
            sage: _.parent()                                                    # optional - fricas
            Univariate Polynomial Ring in x over Integer Ring
            sage: fricas('x^2 + y^2 + 1/2').sage()                              # optional - fricas
            y^2 + x^2 + 1/2
            sage: _.parent()                                                    # optional - fricas
            Multivariate Polynomial Ring in y, x over Rational Field

            sage: fricas("1$Polynomial Integer").sage()                         # optional - fricas
            1

            sage: fricas("x^2/2").sage()                                        # optional - fricas
            1/2*x^2

        Rational functions::

            sage: fricas("x^2 + 1/z").sage()                                    # optional - fricas
            x^2 + 1/z

        Expressions::

            sage: fricas("sin(x+y)/exp(z)*log(1+%e)").sage()                    # optional - fricas
            e^(-z)*log(e + 1)*sin(x + y)

            sage: fricas("factorial(n)").sage()                                 # optional - fricas
            factorial(n)

            sage: fricas("integrate(sin(x+y), x=0..1)").sage()                  # optional - fricas
            -cos(y + 1) + cos(y)

            sage: fricas("integrate(x*sin(1/x), x=0..1)").sage()                # optional - fricas
            'failed'

            sage: fricas("integrate(sin((x^2+1)/x),x)").sage()                  # optional - fricas
            integral(sin((x^2 + 1)/x), x)

        .. TODO::

            - Converting matrices and lists takes much too long.

        Matrices::

            sage: fricas("matrix [[x^n/2^m for n in 0..5] for m in 0..3]").sage()         # optional - fricas, long time
            [      1       x     x^2     x^3     x^4     x^5]
            [    1/2   1/2*x 1/2*x^2 1/2*x^3 1/2*x^4 1/2*x^5]
            [    1/4   1/4*x 1/4*x^2 1/4*x^3 1/4*x^4 1/4*x^5]
            [    1/8   1/8*x 1/8*x^2 1/8*x^3 1/8*x^4 1/8*x^5]

        Lists::

            sage: fricas("[2^n/x^n for n in 0..5]").sage()                      # optional - fricas, long time
            [1, 2/x, 4/x^2, 8/x^3, 16/x^4, 32/x^5]

            sage: fricas("[matrix [[i for i in 1..n]] for n in 0..5]").sage()   # optional - fricas, long time
            [[], [1], [1 2], [1 2 3], [1 2 3 4], [1 2 3 4 5]]

        Error handling::

            sage: s = fricas.guessPade("[fibonacci i for i in 0..10]"); s       # optional - fricas
                n        x
            [[[x ]- ----------]]
                     2
                    x  + x - 1
            sage: s.sage()                                                      # optional - fricas
            Traceback (most recent call last):
            ...
            NotImplementedError: The translation of the FriCAS Expression rootOfADE(n,...()) to sage is not yet implemented.

            sage: s = fricas("series(sqrt(1+x), x=0)"); s                       # optional - fricas
                  1     1  2    1  3    5   4    7   5    21   6    33   7    429   8
              1 + - x - - x  + -- x  - --- x  + --- x  - ---- x  + ---- x  - ----- x
                  2     8      16      128      256      1024      2048      32768
            +
               715   9    2431   10      11
              ----- x  - ------ x   + O(x  )
              65536      262144

            sage: s.sage()                                                      # optional - fricas
            Traceback (most recent call last):
            ...
            NotImplementedError: The translation of the FriCAS object
            <BLANKLINE>
                  1     1  2    1  3    5   4    7   5    21   6    33   7    429   8
              1 + - x - - x  + -- x  - --- x  + --- x  - ---- x  + ---- x  - ----- x
                  2     8      16      128      256      1024      2048      32768
            +
               715   9    2431   10      11
              ----- x  - ------ x   + O(x  )
              65536      262144
            <BLANKLINE>
            to sage is not yet implemented:
            An error occurred when FriCAS evaluated 'unparse(...::InputForm)':
            <BLANKLINE>
               Cannot convert the value from type Any to InputForm .
        """
        from sage.rings.all import ZZ, QQ, QQbar, PolynomialRing, RDF
        from sage.rings.fraction_field import FractionField
        from sage.rings.finite_rings.integer_mod_ring import Integers
        from sage.rings.real_mpfr import RealField
        from sage.symbolic.ring import SR
        from sage.matrix.constructor import matrix
        from sage.structure.factorization import Factorization
        from sage.misc.sage_eval import sage_eval

        # TODO: perhaps we should translate the type first?
        # TODO: perhaps we should get the InputForm as SExpression?

        # remember: fricas.new gives a FriCASElement

        # the coercion to Any gets rid of the Union domain
        P = self._check_valid()
        domain = P.new("dom((%s)::Any)" % self._name) # domain is now a fricas SExpression

        # first translate dummy domains such as "failed". we must not
        # recurse here!
        if P.get_boolean("string?(%s)" % domain._name):
            return P.get_string("string(%s)" % domain._name)

        # now translate domains which cannot be coerced to InputForm,
        # or where we do not need it.
        head = str(domain.car())
        if head == "List":
            n = P.get_integer('#(%s)' %self._name)
            return [P.new('elt(%s,%s)' %(self._name, k)).sage() for k in range(1, n+1)]

        if head == "Matrix":
            base_ring = self._get_sage_type(domain[1])
            rows = P.new('listOfLists(%s)' %self._name).sage()
            return matrix(base_ring, rows)

        if head == "Fraction":
            return P.new("numer(%s)" %self._name).sage()/P.new("denom(%s)" %self._name).sage()

        if head == "Factored":
            l = P.new('[[f.factor, f.exponent] for f in factors(%s)]' %self._name).sage()
            return Factorization([(p, e) for p,e in l])

        # finally translate domains with InputForm
        try:
            unparsed_InputForm = P.get_unparsed_InputForm(self._name)
        except RuntimeError as error:
            raise NotImplementedError("The translation of the FriCAS object\n\n%s\n\nto sage is not yet implemented:\n%s" %(self, error))

        if head == "Boolean":
            return unparsed_InputForm == "true"

        if head in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ(unparsed_InputForm)

        if head == "String":
            return unparsed_InputForm

        if head == "Float":
            # Warning: precision$Float gives the current precision,
            # whereas length(mantissa(self)) gives the precision of
            # self.
            prec = max(P.new("length mantissa(%s)" %self._name).sage(), 53)
            R = RealField(prec)
            x, e, b = unparsed_InputForm.lstrip('float(').rstrip(')').split(',')
            return R(ZZ(x)*ZZ(b)**ZZ(e))

        if head == "DoubleFloat":
            return RDF(unparsed_InputForm)

        if head == "AlgebraicNumber":
            s = unparsed_InputForm[:-len("::AlgebraicNumber()")]
            return sage_eval("QQbar(" + s + ")")

        if head == "IntegerMod":
            # one might be tempted not to go via InputForm here, but
            # it turns out to be safer to do it.
            n = unparsed_InputForm[len("index("):]
            n = n[:n.find(")")]
            return self._get_sage_type(domain)(n)

        if head == "Polynomial":
            base_ring = self._get_sage_type(domain[1])
            # the following is a bad hack, we should be getting a list here
            vars = P.get_unparsed_InputForm("variables(%s)" %self._name)[1:-1]
            if vars == "":
                return base_ring(unparsed_InputForm)
            else:
                R = PolynomialRing(base_ring, vars)
                return R(unparsed_InputForm)

        if head == "OrderedCompletion":
            # this is a workaround, I don't know how translate this
            if str(domain[1].car()) == "Expression":
                return self._sage_expression(unparsed_InputForm)

        if head == "Expression":
            # TODO: we also have Expression Complex Integer and the like
            if str(domain[1].car()) == "Integer":
                return self._sage_expression(unparsed_InputForm)

        if head == 'DistributedMultivariatePolynomial':
            base_ring = self._get_sage_type(domain[2])
            vars = domain[1].car()
            R = PolynomialRing(base_ring, vars)
            return R(unparsed_InputForm)

        raise NotImplementedError("The translation of the FriCAS object %s to sage is not yet implemented." %(unparsed_InputForm))

class FriCASFunctionElement(FunctionElement):
    def __init__(self, object, name):
        """
        Make FriCAS operation names valid python function identifiers.

        TESTS::

            sage: a = fricas('"Hello"')                                         # optional - fricas
            sage: a.upperCase_q                                                 # optional - fricas
            upperCase?
            sage: a.upperCase_e                                                 # optional - fricas
            upperCase!
            sage: a.upperCase_e()                                               # optional - fricas
            "HELLO"

        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        FunctionElement.__init__(self, object, name)

    pass

class FriCASExpectFunction(ExpectFunction):
    def __init__(self, parent, name):
        """
        Translate the pythonized function identifier back to a FriCAS
        operation name.

        TESTS::

            sage: fricas.upperCase_q                                            # optional - fricas
            upperCase?
            sage: fricas.upperCase_e                                            # optional - fricas
            upperCase!

        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        ExpectFunction.__init__(self, parent, name)

    pass

def is_FriCASElement(x):
    """
    Return ``True`` if ``x`` is of type :class:`FriCASElement`.

    EXAMPLES::

        sage: from sage.interfaces.fricas import is_FriCASElement               # optional - fricas
        sage: is_FriCASElement(fricas(2))                                       # optional - fricas
        True
        sage: is_FriCASElement(2)                                               # optional - fricas
        False
    """
    return isinstance(x, FriCASElement)

fricas = FriCAS()

def reduce_load_fricas():
    """
    Returns the FriCAS interface object defined in
    :sage.interfaces.fricas.

    EXAMPLES::

        sage: from sage.interfaces.fricas import reduce_load_fricas             # optional - fricas
        sage: reduce_load_fricas()                                              # optional - fricas
        FriCAS
    """
    return fricas

import os

def fricas_console():
    """
    Spawn a new FriCAS command-line session.

    EXAMPLES::

        sage: fricas_console()                                                  # not tested
                         FriCAS (AXIOM fork) Computer Algebra System
                                    Version: FriCAS 1.0.5
                     Timestamp: Thursday February 19, 2009 at 06:57:33
        -----------------------------------------------------------------------------
           Issue )copyright to view copyright notices.
           Issue )summary for a summary of useful system commands.
           Issue )quit to leave AXIOM and return to shell.
        -----------------------------------------------------------------------------
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%fricas magics instead.')
    os.system('fricas -nox')

def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.fricas import __doctest_cleanup              # optional - fricas
        sage: a = FriCAS()                                                      # optional - fricas
        sage: two = a(2)                                                        # optional - fricas
        sage: a.is_running()                                                    # optional - fricas
        True
        sage: __doctest_cleanup()                                               # optional - fricas
        sage: a.is_running()                                                    # optional - fricas
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
