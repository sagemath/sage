r"""
Interface to Axiom

.. TODO::

    - Evaluation using a file is not done. Any input line with more than a
      few thousand characters would hang the system, so currently it
      automatically raises an exception.

    - All completions of a given command.

    - Interactive help.

Axiom is a free GPL-compatible (modified BSD license) general
purpose computer algebra system whose development started in 1973
at IBM. It contains symbolic manipulation algorithms, as well as
implementations of special functions, including elliptic functions
and generalized hypergeometric functions. Moreover, Axiom has
implementations of many functions relating to the invariant theory
of the symmetric group `S_n.` For many links to Axiom
documentation see http://wiki.axiom-developer.org.

AUTHORS:

- Bill Page (2006-10): Created this (based on Maxima interface)


  .. NOTE::

     Bill Page put a huge amount of effort into the Sage Axiom
     interface over several days during the Sage Days 2 coding
     sprint. This is contribution is greatly appreciated.

- William Stein (2006-10): misc touchup.

- Bill Page (2007-08): Minor modifications to support axiom4sage-0.3

.. NOTE::

   The axiom4sage-0.3.spkg is based on an experimental version of the
   FriCAS fork of the Axiom project by Waldek Hebisch that uses
   pre-compiled cached Lisp code to build Axiom very quickly with
   clisp.

If the string "error" (case insensitive) occurs in the output of
anything from axiom, a :exc:`RuntimeError` exception is raised.

EXAMPLES: We evaluate a very simple expression in axiom.

::

    sage: axiom('3 * 5')                     #optional - axiom
    15
    sage: a = axiom(3) * axiom(5); a         #optional - axiom
    15

The type of a is AxiomElement, i.e., an element of the axiom
interpreter.

::

    sage: type(a)                            #optional - axiom
    <class 'sage.interfaces.axiom.AxiomElement'>
    sage: parent(a)                          #optional - axiom
    Axiom

The underlying Axiom type of a is also available, via the type
method::

    sage: a.type()                           #optional - axiom
    PositiveInteger

We factor `x^5 - y^5` in Axiom in several different ways.
The first way yields a Axiom object.

::

    sage: F = axiom.factor('x^5 - y^5'); F      #optional - axiom
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    sage: type(F)                               #optional - axiom
    <class 'sage.interfaces.axiom.AxiomElement'>
    sage: F.type()                              #optional - axiom
    Factored Polynomial Integer

Note that Axiom objects are normally displayed using "ASCII art".

::

    sage: a = axiom(2/3); a          #optional - axiom
      2
      -
      3
    sage: a = axiom('x^2 + 3/7'); a      #optional - axiom
       2   3
      x  + -
           7

The ``axiom.eval`` command evaluates an expression in
axiom and returns the result as a string. This is exact as if we
typed in the given line of code to axiom; the return value is what
Axiom would print out.

::

    sage: print(axiom.eval('factor(x^5 - y^5)'))   # optional - axiom
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    Type: Factored Polynomial Integer

We can create the polynomial `f` as a Axiom polynomial,
then call the factor method on it. Notice that the notation
``f.factor()`` is consistent with how the rest of Sage
works.

::

    sage: f = axiom('x^5 - y^5')                  #optional - axiom
    sage: f^2                                     #optional - axiom
       10     5 5    10
      y   - 2x y  + x
    sage: f.factor()                              #optional - axiom
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

Control-C interruption works well with the axiom interface, because
of the excellent implementation of axiom. For example, try the
following sum but with a much bigger range, and hit control-C.

::

    sage:  f = axiom('(x^5 - y^5)^10000')       # not tested
    Interrupting Axiom...
    ...
    <class 'exceptions.TypeError'>: Ctrl-c pressed while running Axiom

::

    sage: axiom('1/100 + 1/101')                  #optional - axiom
       201
      -----
      10100
    sage: a = axiom('(1 + sqrt(2))^5'); a         #optional - axiom
         +-+
      29\|2  + 41

TESTS:

We check to make sure the subst method works with keyword
arguments.

::

    sage: a = axiom(x+2); a  #optional - axiom
    x + 2
    sage: a.subst(x=3)       #optional - axiom
    5

We verify that Axiom floating point numbers can be converted to
Python floats.

::

    sage: float(axiom(2))     #optional - axiom
    2.0
"""

###########################################################################
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#                     2007 Bill Page
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
###########################################################################

import os
import re

import sage.interfaces.abc

from .expect import Expect, ExpectElement, FunctionElement, ExpectFunction
from sage.env import DOT_SAGE
from pexpect import EOF
from sage.misc.multireplace import multiple_replace
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.misc.instancedoc import instancedoc
from sage.structure.richcmp import rich_to_bool


# The Axiom commands ")what thing det" ")show Matrix" and ")display
# op det" commands, gives a list of all identifiers that begin in
# a certain way.  This could maybe be useful somehow... (?)  Also
# axiom has a lot a lot of ways for getting documentation from the
# system -- this could also be useful.

class PanAxiom(ExtraTabCompletion, Expect):
    """
    Interface to a PanAxiom interpreter.
    """
    def __init__(self, name='axiom', command='axiom -nox -noclef',
                 script_subdirectory=None, logfile=None,
                 server=None, server_tmpdir=None,
                 init_code=[')lisp (si::readline-off)']):
        """
        Create an instance of the Axiom interpreter.

        TESTS::

            sage: axiom == loads(dumps(axiom))
            True
        """
        eval_using_file_cutoff = 200
        self.__eval_using_file_cutoff = eval_using_file_cutoff
        self._COMMANDS_CACHE = '%s/%s_commandlist_cache.sobj' % (DOT_SAGE, name)
        Expect.__init__(self,
                        name=name,
                        prompt=r'\([0-9]+\) -> ',
                        command=command,
                        script_subdirectory=script_subdirectory,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        restart_on_ctrlc=False,
                        verbose_start=False,
                        init_code=init_code,
                        logfile=logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)
        self._prompt_wait = self._prompt

    def _start(self):
        """
        Start the Axiom interpreter.

        EXAMPLES::

            sage: a = Axiom()
            sage: a.is_running()
            False
            sage: a._start()     #optional - axiom
            sage: a.is_running() #optional - axiom
            True
            sage: a.quit()       #optional - axiom
        """
        Expect._start(self)
        self._eval_line(')set functions compile on', reformat=False)
        self._eval_line(')set output length 245', reformat=False)
        self._eval_line(')set message autoload off', reformat=False)

    def _read_in_file_command(self, filename):
        r"""
        EXAMPLES::

            sage: axiom._read_in_file_command('test.input')
            ')read test.input \n'
            sage: axiom._read_in_file_command('test')
            Traceback (most recent call last):
            ...
            ValueError: the filename must end with .input

        ::

            sage: filename = tmp_filename(ext='.input')
            sage: f = open(filename, 'w')
            sage: _ = f.write('xx := 22;\n')
            sage: f.close()
            sage: axiom.read(filename)    # optional - axiom
            sage: axiom.get('xx')         # optional - axiom
            '22'
        """
        if not filename.endswith('.input'):
            raise ValueError("the filename must end with .input")

        # For some reason this trivial comp
        # keeps certain random freezes from occurring.  Do not remove this.
        # The space before the \n is also important.
        return ')read %s \n' % filename

    def _quit_string(self):
        """
        Return the string used to quit Axiom.

        EXAMPLES::

            sage: axiom._quit_string()
            ')lisp (quit)'
            sage: a = Axiom()
            sage: a.is_running()
            False
            sage: a._start()     #optional - axiom
            sage: a.is_running() #optional - axiom
            True
            sage: a.quit()       #optional - axiom
            sage: a.is_running() #optional - axiom
            False
        """
        return ')lisp (quit)'

    def _commands(self):
        """
        Return a list of commands available. This is done by parsing the
        result of the first section of the output of ')what things'.

        EXAMPLES::

            sage: # optional - axiom
            sage: cmds = axiom._commands()
            sage: len(cmds) > 100
            True
            sage: '<' in cmds
            True
            sage: 'factor' in cmds
            True
        """
        s = self.eval(")what things")
        start = '\r\n\r\n#'
        i = s.find(start)
        end = "To get more information about"
        j = s.find(end)
        return s[i + len(start):j].split()

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        """
        Return a list of all the commands defined in Axiom and optionally
        (per default) store them to disk.

        EXAMPLES::

            sage: # optional - axiom
            sage: c = axiom._tab_completion(use_disk_cache=False, verbose=False)
            sage: len(c) > 100
            True
            sage: 'factor' in c
            True
            sage: '**' in c
            False
            sage: 'upperCase?' in c
            False
            sage: 'upperCase_q' in c
            True
            sage: 'upperCase_e' in c
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
                except OSError:
                    pass
            if verbose:
                print("\nBuilding %s command completion list (this takes" % self)
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % self._COMMANDS_CACHE)
            v = self._commands()

            # Process we now need process the commands to strip out things which
            # are not valid Python identifiers.
            valid = re.compile('[^a-zA-Z0-9_]+')
            names = [x for x in v if valid.search(x) is None]

            # Change everything that ends with ? to _q and
            # everything that ends with ! to _e
            names += [x[:-1] + "_q" for x in v if x.endswith("?")]
            names += [x[:-1] + "_e" for x in v if x.endswith("!")]

            self.__tab_completion = names
            if len(v) > 200:
                # Axiom is actually installed.
                sage.misc.persist.save(v, self._COMMANDS_CACHE)
            return names

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: axiom.set('xx', '2')    #optional - axiom
            sage: axiom.get('xx')         #optional - axiom
            '2'
        """
        cmd = '%s := %s' % (var, value)
        out = self._eval_line(cmd, reformat=False)

        if out.find("error") != -1:
            raise TypeError("Error executing code in Axiom\nCODE:\n\t%s\nAxiom ERROR:\n\t%s" % (cmd, out))

    def get(self, var):
        r"""
        Get the string value of the Axiom variable var.

        EXAMPLES::

            sage: # optional - axiom
            sage: axiom.set('xx', '2')
            sage: axiom.get('xx')
            '2'
            sage: a = axiom('(1 + sqrt(2))^5')
            sage: axiom.get(a.name())
            '     +-+\r\r\n  29\\|2  + 41'
        """
        s = self._eval_line(str(var))
        i = s.rfind('Type:')
        s = s[:i].rstrip().lstrip("\n")
        if '\n' not in s:
            s = s.strip()
        return s

    def _eval_line(self, line, reformat=True, allow_use_file=False,
                   wait_for_prompt=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: print(axiom._eval_line('2+2'))  # optional - axiom
              4
                                                       Type: PositiveInteger
        """
        from sage.misc.verbose import verbose
        if not wait_for_prompt:
            return Expect._eval_line(self, line)
        line = line.rstrip().rstrip(';')
        if line == '':
            return ''
        if len(line) > 3000:
            raise NotImplementedError("evaluation of long input lines (>3000 characters) in Axiom not yet implemented.")
        if self._expect is None:
            self._start()
        if allow_use_file and self.__eval_using_file_cutoff and \
           len(line) > self.__eval_using_file_cutoff:
            return self._eval_line_using_file(line)
        try:
            E = self._expect
            # debug
            # self._synchronize(cmd='1+%s\n')
            verbose("in = '%s'" % line, level=3)
            E.sendline(line)
            self._expect.expect(self._prompt)
            out = self._expect.before
            # debug
            verbose("out = '%s'" % out, level=3)
        except EOF:
            if self._quit_string() in line:
                return ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()

        if '>> Error detected within library code:' in out or \
           'Cannot find a definition or applicable library operation named' in out:
            raise RuntimeError(out)

        if not reformat:
            return out
        if 'error' in out:
            return out
        # out = out.lstrip()
        i = out.find('\n')
        out = out[i+1:]
        outs = out.split("\n")
        i = 0
        for line in outs:
            line = line.rstrip()
            if line[:4] == '   (':
                i = line.find('(')
                i += line[i:].find(')') + 1
                if line[i:] == "":
                    i = 0
                    outs = outs[1:]
                break
        return "\n".join(line[i:] for line in outs[1:])

    # define relational operators
    def _equality_symbol(self):
        """equality symbol

        EXAMPLES::

            sage: a = axiom(x==6); a    #optional axiom
            x= 6
        """
        return "="


class Axiom(PanAxiom):
    def __reduce__(self):
        """
        EXAMPLES::

            sage: Axiom().__reduce__()
            (<function reduce_load_Axiom at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            Axiom
        """
        return reduce_load_Axiom, tuple([])

    def _function_class(self):
        """
        Return the AxiomExpectFunction class.

        EXAMPLES::

            sage: axiom._function_class()
            <class 'sage.interfaces.axiom.PanAxiomExpectFunction'>
            sage: type(axiom.gcd)
            <class 'sage.interfaces.axiom.PanAxiomExpectFunction'>
        """
        return AxiomExpectFunction

    def _object_class(self):
        """
        EXAMPLES::

            sage: axiom._object_class()
            <class 'sage.interfaces.axiom.PanAxiomElement'>
            sage: type(axiom(2)) #optional - axiom
            <class 'sage.interfaces.axiom.PanAxiomElement'>
        """
        return AxiomElement

    def _function_element_class(self):
        """
        Return the Axiom function element class.

        EXAMPLES::

            sage: axiom._function_element_class()
            <class 'sage.interfaces.axiom.PanAxiomFunctionElement'>
            sage: type(axiom(2).gcd) #optional - axiom
            <class 'sage.interfaces.axiom.PanAxiomFunctionElement'>
        """
        return AxiomFunctionElement

    def console(self):
        """
        Spawn a new Axiom command-line session.

        EXAMPLES::

            sage: axiom.console() #not tested
                                    AXIOM Computer Algebra System
                                    Version: Axiom (January 2009)
                           Timestamp: Sunday January 25, 2009 at 07:08:54
            -----------------------------------------------------------------------------
               Issue )copyright to view copyright notices.
               Issue )summary for a summary of useful system commands.
               Issue )quit to leave AXIOM and return to shell.
            -----------------------------------------------------------------------------
        """
        axiom_console()


@instancedoc
class PanAxiomElement(ExpectElement, sage.interfaces.abc.AxiomElement):
    def __call__(self, x):
        """
        EXAMPLES::

            sage: f = axiom(x+2) #optional - axiom
            sage: f(2)           #optional - axiom
            4
        """
        self._check_valid()
        P = self.parent()
        return P('%s(%s)' % (self.name(), x))

    def _richcmp_(self, other, op):
        """
        EXAMPLES::

            sage: # optional - axiom
            sage: two = axiom(2)
            sage: two == 2
            True
            sage: two == 3
            False
            sage: two < 3
            True
            sage: two > 1
            True

            sage: # optional - axiom
            sage: a = axiom(1); b = axiom(2)
            sage: a == b
            False
            sage: a < b
            True
            sage: a > b
            False
            sage: b < a
            False
            sage: b > a
            True

        We can also compare more complicated object such as functions::

            sage: f = axiom('sin(x)'); g = axiom('cos(x)')    #optional - axiom
            sage: f == g                                      #optional - axiom
            False
        """
        P = self.parent()
        if 'true' in P.eval("(%s = %s) :: Boolean" % (self.name(), other.name())):
            return rich_to_bool(op, 0)
        elif 'true' in P.eval("(%s < %s) :: Boolean" % (self.name(), other.name())):
            return rich_to_bool(op, -1)
        elif 'true' in P.eval("(%s > %s) :: Boolean" % (self.name(), other.name())):
            return rich_to_bool(op, 1)

        return NotImplemented

    def type(self):
        """
        Return the type of an AxiomElement.

        EXAMPLES::

            sage: axiom(x+2).type()  #optional - axiom
            Polynomial Integer
        """
        P = self._check_valid()
        s = P._eval_line(self.name())
        i = s.rfind('Type:')
        return P(s[i+5:].strip())

    def __len__(self):
        """
        Return the length of a list.

        EXAMPLES::

            sage: v = axiom('[x^i for i in 0..5]')            # optional - axiom
            sage: len(v)                                      # optional - axiom
            6
        """
        P = self._check_valid()
        s = P.eval('# %s ' % self.name())
        i = s.rfind('Type')
        return int(s[:i - 1])

    def __getitem__(self, n):
        r"""
        Return the `n`-th element of this list.

        .. NOTE::

           Lists are 1-based.

        EXAMPLES::

            sage: # optional - axiom
            sage: v = axiom('[i*x^i for i in 0..5]'); v
                     2   3   4   5
              [0,x,2x ,3x ,4x ,5x ]
            sage: v[4]
                3
              3x
            sage: v[1]
            0
            sage: v[10]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        n = int(n)
        if n <= 0 or n > len(self):
            raise IndexError("index out of range")
        P = self._check_valid()
        if not isinstance(n, tuple):
            return P.new('%s(%s)' % (self._name, n))
        else:
            return P.new('%s(%s)' % (self._name, str(n)[1:-1]))

    def comma(self, *args):
        """
        Return an Axiom tuple from ``self`` and ``args``.

        EXAMPLES::

            sage: # optional - axiom
            sage: two = axiom(2)
            sage: two.comma(3)
            [2,3]
            sage: two.comma(3,4)
            [2,3,4]
            sage: _.type()
            Tuple PositiveInteger
        """
        P = self._check_valid()
        args = list(args)
        for i, arg in enumerate(args):
            if not isinstance(arg, AxiomElement) or arg.parent() is not P:
                args[i] = P(arg)
        cmd = "(" + ",".join(x.name() for x in [self] + args) + ")"
        return P(cmd)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: a = axiom(1/2) #optional - axiom
            sage: latex(a)       #optional - axiom
            \frac{1}{2}
        """
        self._check_valid()
        P = self.parent()
        s = P._eval_line('outputAsTex(%s)' % self.name(), reformat=False)
        if '$$' not in s:
            raise RuntimeError("Error texing axiom object.")
        i = s.find('$$')
        j = s.rfind('$$')
        s = s[i + 2:j]
        s = multiple_replace({'\r': '', '\n': ' ',
                              ' \\sp ': '^',
                              '\\arcsin ': '\\sin^{-1} ',
                              '\\arccos ': '\\cos^{-1} ',
                              '\\arctan ': '\\tan^{-1} '},
                             re.sub(r'\\leqno\(.*?\)', '', s))  # no eq number!
        return s

    def as_type(self, type):
        """
        Return ``self`` as type.

        EXAMPLES::

            sage: a = axiom(1.2); a            #optional - axiom
            1.2
            sage: a.as_type(axiom.DoubleFloat) #optional - axiom
            1.2
            sage: _.type()                     #optional - axiom
            DoubleFloat
        """
        P = self._check_valid()
        type = P(type)
        return P.new("%s :: %s" % (self.name(), type.name()))

    def unparsed_input_form(self):
        """
        Get the linear string representation of this object, if possible
        (often it isn't).

        EXAMPLES::

            sage: a = axiom(x^2+1); a     #optional - axiom
               2
              x  + 1
            sage: a.unparsed_input_form() #optional - axiom
            'x*x+1'
        """
        P = self._check_valid()
        s = P.eval('unparse(%s::InputForm)' % self._name)
        if 'translation error' in s or 'Cannot convert' in s:
            raise NotImplementedError
        s = multiple_replace({'\r\n': '',  # fix stupid Fortran-ish
                              'DSIN(': 'sin(',
                              'DCOS(': 'cos(',
                              'DTAN(': 'tan(',
                              'DSINH(': 'sinh('}, s)
        r = re.search(r'"(.*)"', s)
        return r.groups(0)[0] if r else s

    def _sage_(self):
        """
        Convert ``self`` to a Sage object.

        EXAMPLES::

            sage: a = axiom(1/2); a #optional - axiom
              1
              -
              2
            sage: a.sage()          #optional - axiom
            1/2
            sage: _.parent()        #optional - axiom
            Rational Field

            sage: gp(axiom(1/2))    #optional - axiom
            1/2

        DoubleFloat's in Axiom are converted to be in RDF in Sage.

        ::

            sage: axiom(2.0).as_type('DoubleFloat').sage()  #optional - axiom
            2.0
            sage: _.parent() #optional - axiom
            Real Double Field


            sage: axiom(2.1234)._sage_() #optional - axiom
            2.12340000000000
            sage: _.parent()             #optional - axiom
            Real Field with 53 bits of precision
            sage: a = RealField(100)(pi)                                                # needs sage.symbolic
            sage: axiom(a)._sage_()             # optional - axiom                      # needs sage.symbolic
            3.1415926535897932384626433833
            sage: _.parent()                    # optional - axiom                      # needs sage.symbolic
            Real Field with 100 bits of precision
            sage: axiom(a)._sage_() == a        # optional - axiom                      # needs sage.symbolic
            True
            sage: axiom(2.0)._sage_() #optional - axiom
            2.00000000000000
            sage: _.parent() #optional  - axiom
            Real Field with 53 bits of precision


        We can also convert Axiom's polynomials to Sage polynomials.
            sage: # optional - axiom, needs sage.symbolic
            sage: a = axiom(x^2 + 1)
            sage: a.type()
            Polynomial Integer
            sage: a.sage()
            x^2 + 1
            sage: _.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: axiom('x^2 + y^2 + 1/2').sage()
            y^2 + x^2 + 1/2
            sage: _.parent()
            Multivariate Polynomial Ring in y, x over Rational Field
        """
        P = self._check_valid()
        type = str(self.type())

        if type in ["Type", "Domain"]:
            return self._sage_domain()

        if type == "Float":
            from sage.rings.real_mpfr import RealField
            from sage.rings.integer_ring import ZZ
            prec = max(self.mantissa().length()._sage_(), 53)
            R = RealField(prec)
            x, e, b = self.unparsed_input_form().lstrip('float(').rstrip(')').split(',')
            return R(ZZ(x) * ZZ(b)**ZZ(e))
        elif type == "DoubleFloat":
            from sage.rings.real_double import RDF
            return RDF(repr(self))
        elif type in ["PositiveInteger", "Integer"]:
            from sage.rings.integer_ring import ZZ
            return ZZ(repr(self))
        elif type.startswith('Polynomial'):
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            base_ring = P(type.removeprefix('Polynomial '))._sage_domain()
            vars = str(self.variables())[1:-1]
            R = PolynomialRing(base_ring, vars)
            return R(self.unparsed_input_form())
        elif type.startswith('Fraction'):
            return self.numer().sage() / self.denom().sage()

        # If all else fails, try using the unparsed input form
        try:
            import sage.misc.sage_eval
            vars = sage.symbolic.ring.var(str(self.variables())[1:-1])
            if isinstance(vars, tuple):
                return sage.misc.sage_eval.sage_eval(self.unparsed_input_form(),
                                                     locals={str(x): x for x in vars})
            else:
                return sage.misc.sage_eval.sage_eval(self.unparsed_input_form(),
                                                     locals={str(vars): vars})
        except Exception:
            raise NotImplementedError

    def _sage_domain(self):
        """
        A helper function for converting Axiom domains to the corresponding
        Sage object.

        EXAMPLES::

            sage: axiom('Integer').sage()  #optional - axiom
            Integer Ring

            sage: axiom('Fraction Integer').sage()  #optional - axiom
            Rational Field

            sage: axiom('DoubleFloat').sage()  #optional - axiom
            Real Double Field
        """
        P = self._check_valid()
        name = str(self)
        if name == 'Integer':
            from sage.rings.integer_ring import ZZ
            return ZZ
        elif name == 'DoubleFloat':
            from sage.rings.real_double import RDF
            return RDF
        elif name.startswith('Fraction '):
            return P(name.lstrip('Fraction '))._sage_domain().fraction_field()

        raise NotImplementedError


AxiomElement = PanAxiomElement


@instancedoc
class PanAxiomFunctionElement(FunctionElement):
    def __init__(self, object, name):
        """
        TESTS::

            sage: # optional - axiom
            sage: a = axiom('"Hello"')
            sage: a.upperCase_q
            upperCase?
            sage: a.upperCase_e
            upperCase!
            sage: a.upperCase_e()
            "HELLO"
        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        FunctionElement.__init__(self, object, name)


AxiomFunctionElement = PanAxiomFunctionElement


@instancedoc
class PanAxiomExpectFunction(ExpectFunction):
    def __init__(self, parent, name):
        """
        TESTS::

            sage: axiom.upperCase_q
            upperCase?
            sage: axiom.upperCase_e
            upperCase!
        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        ExpectFunction.__init__(self, parent, name)


AxiomExpectFunction = PanAxiomExpectFunction


def is_AxiomElement(x):
    """
    Return ``True`` if ``x`` is of type :class:`AxiomElement`.

    EXAMPLES::

        sage: from sage.interfaces.axiom import is_AxiomElement
        sage: is_AxiomElement(2)
        doctest:...: DeprecationWarning: the function is_AxiomElement is deprecated; use isinstance(x, sage.interfaces.abc.AxiomElement) instead
        See https://github.com/sagemath/sage/issues/34804 for details.
        False
        sage: is_AxiomElement(axiom(2))  # optional - axiom
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(34804, "the function is_AxiomElement is deprecated; use isinstance(x, sage.interfaces.abc.AxiomElement) instead")

    return isinstance(x, AxiomElement)


# Instances
axiom = Axiom(name='axiom')


def reduce_load_Axiom():
    """
    Return the Axiom interface object defined in
    sage.interfaces.axiom.

    EXAMPLES::

        sage: from sage.interfaces.axiom import reduce_load_Axiom
        sage: reduce_load_Axiom()
        Axiom
    """
    return axiom


def axiom_console():
    """
    Spawn a new Axiom command-line session.

    EXAMPLES::

        sage: axiom_console() #not tested
                                AXIOM Computer Algebra System
                                Version: Axiom (January 2009)
                       Timestamp: Sunday January 25, 2009 at 07:08:54
        -----------------------------------------------------------------------------
           Issue )copyright to view copyright notices.
           Issue )summary for a summary of useful system commands.
           Issue )quit to leave AXIOM and return to shell.
        -----------------------------------------------------------------------------
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%axiom magics instead.')
    os.system('axiom -nox')
