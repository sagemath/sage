r"""
Pexpect Interface to Giac
(You should prefer the cython interface: giacpy_sage and its libgiac command)

(adapted by F. Han from William Stein and Gregg Musiker maple's interface)

You must have the  Giac interpreter installed
and available as the command ``giac`` in your PATH in
order to use this interface. You need a giac version
supporting "giac --sage" ( roughly after 0.9.1 ). In this case you do not have
to install any  optional Sage packages. If giac is not already installed, you can
download binaries or sources or spkg (follow the sources link) from the homepage:

Homepage <https://www-fourier.ujf-grenoble.fr/~parisse/giac.html>

Type ``giac.[tab]`` for a list of all the functions
available from your Giac install. Type
``giac.[tab]?`` for Giac's help about a given
function. Type ``giac(...)`` to create a new Giac
object, and ``giac.eval(...)`` to run a string using
Giac (and get the result back as a string).

If the giac spkg is installed, you should find the full html documentation there::

    $SAGE_LOCAL/share/giac/doc/en/cascmd_local/index.html


EXAMPLES::

    sage: giac('3 * 5')
    15
    sage: giac.eval('ifactor(2005)')
    '5*401'
    sage: giac.ifactor(2005)
    2005
    sage: l=giac.ifactors(2005) ; l; l[2]
    [5,1,401,1]
    401
    sage: giac.fsolve('x^2=cos(x)+4', 'x','0..5')
    [1.9140206190...
    sage: giac.factor('x^4 - y^4')
    (x-y)*(x+y)*(x^2+y^2)
    sage: R.<x,y>=QQ[];f=(x+y)^5;f2=giac(f);(f-f2).normal()
    0
    sage: x,y=giac('x,y'); giac.int(y/(cos(2*x)+cos(x)),x)     # random
    y*2*((-(tan(x/2)))/6+(-2*1/6/sqrt(3))*ln(abs(6*tan(x/2)-2*sqrt(3))/abs(6*tan(x/2)+2*sqrt(3))))


If the string "error" (case insensitive) occurs in the output of
anything from Giac, a RuntimeError exception is raised.

Tutorial
--------

AUTHORS:

- Gregg Musiker (2006-02-02): initial version.

- Frederic Han: adapted to giac.

- Marcelo Forets (2017-04-06): conversions and cleanup.

This tutorial is based on the Maple Tutorial for number theory from
http://www.math.mun.ca/~drideout/m3370/numtheory.html.

Syntax
~~~~~~~

There are several ways to use the Giac Interface in Sage. We will
discuss two of those ways in this tutorial.


#. If you have a giac expression such as

   ::

       factor( (x^4-1));

   We can write that in sage as

   ::

       sage: giac('factor(x^4-1)')
       (x-1)*(x+1)*(x^2+1)

   Notice, there is no need to use a semicolon.

#. Since Sage is written in Python, we can also import giac
   commands and write our scripts in a pythonic way. For example,
   ``factor()`` is a giac command, so we can also factor
   in Sage using

   ::

       sage: giac('(x^4-1)').factor()
       (x-1)*(x+1)*(x^2+1)

   where ``expression.command()`` means the same thing as
   ``command(expression)`` in Giac. We will use this
   second type of syntax whenever possible, resorting to the first
   when needed.

   ::

       sage: giac('(x^12-1)/(x-1)').normal()
       x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1

Some typical input
~~~~~~~~~~~~~~~~~~

The normal command will reduce a rational function to the
lowest terms. In giac, simplify is slower than normal because it
tries more sophisticated simplifications (ex algebraic extensions)
The factor command will factor a polynomial with
rational coefficients into irreducible factors over the ring of
integers (if your default configuration of giac (cf .xcasrc) has not
allowed square roots). So for example,


::

    sage: giac('(x^12-1)').factor( )
    (x-1)*(x+1)*(x^2+1)*(x^2-x+1)*(x^2+x+1)*(x^4-x^2+1)

::

    sage: giac('(x^28-1)').factor( )
    (x-1)*(x+1)*(x^2+1)*(x^6-x^5+x^4-x^3+x^2-x+1)*(x^6+x^5+x^4+x^3+x^2+x+1)*(x^12-x^10+x^8-x^6+x^4-x^2+1)

Giac console
~~~~~~~~~~~~~

Another important feature of giac is its online help. We can
access this through sage as well. After reading the description of
the command, you can press q to immediately get back to your
original prompt.

Incidentally you can always get into a giac console by the
command ::

    sage: giac.console()                       # not tested
    sage: !giac                                # not tested

Note that the above two commands are slightly different, and the
first is preferred.

For example, for help on the giac command factors, we type ::

    sage: giac.help('factors')                     # not tested

::

    sage: alpha = giac((1+sqrt(5))/2)
    sage: beta = giac(1-sqrt(5))/2
    sage: f19  = alpha^19 - beta^19/sqrt(5)
    sage: f19
    (sqrt(5)/2+1/2)^19-((-sqrt(5)+1)/2)^19/sqrt(5)
    sage: (f19-(5778*sqrt(5)+33825)/5).normal()
    0

Function definitions
~~~~~~~~~~~~~~~~~~~~

Let's say we want to write a giac program now that squares a
number if it is positive and cubes it if it is negative. In giac,
that would look like

::

    mysqcu := proc(x)
    if x > 0 then x^2;
    else x^3; fi;
    end;

In Sage, we write

::

    sage: mysqcu = giac('proc(x) if x > 0 then x^2 else x^3 fi end')
    sage: mysqcu(5)
    25
    sage: mysqcu(-5)
    -125

More complicated programs should be put in a separate file and
loaded.

Conversions
~~~~~~~~~~~~

The ``GiacElement.sage()`` method tries to convert a Giac object to a Sage
object. In many cases, it will just work. In particular, it should be able to
convert expressions entirely consisting of:

- numbers, i.e. integers, floats, complex numbers;
- functions and named constants also present in Sage, where Sage knows how to
  translate the function or constant's name from Giac's
- symbolic variables whose names don't pathologically overlap with
  objects already defined in Sage.

This method will not work when Giac's output includes functions unknown to Sage.

If you want to convert more complicated Giac expressions, you can
instead call ``GiacElement._sage_()`` and supply a translation dictionary::

    sage: g = giac('NewFn(x)')
    sage: g._sage_(locals={'NewFn': sin})
    sin(x)

Moreover, new conversions can be permanently added using Pynac's
``register_symbol``, and this is the recommended approach for library code.
For more details, see the documentation for ``._sage_()``.

TESTS:

Test that conversion of symbolic functions with latex names works (:trac:`31047`)::

    sage: var('phi')
    phi
    sage: function('Cp', latex_name='C_+')
    Cp
    sage: test = Cp(phi)._giac_()._sage_()
    sage: test.operator() == Cp
    True
    sage: test.operator()._latex_() == 'C_+'
    True
"""

#############################################################################
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#############################################################################

import os

from sage.interfaces.expect import Expect, ExpectElement, ExpectFunction, FunctionElement, gc_disabled

import pexpect

from sage.cpython.string import bytes_to_str
from sage.env import DOT_SAGE
from sage.misc.pager import pager
from sage.docs.instancedoc import instancedoc
from sage.structure.richcmp import rich_to_bool


COMMANDS_CACHE = '%s/giac_commandlist_cache.sobj'%DOT_SAGE

class Giac(Expect):
    r"""
    Interface to the Giac interpreter.

    You must have the optional  Giac interpreter installed and available as the command ``giac`` in your PATH in order to use this interface. Try the command: print(giac._install_hints()) for more informations on giac installation.

    Type ``giac.[tab]`` for a list of all the functions available from your Giac install.
    Type ``giac.[tab]?`` for Giac's help about a given function.
    Type ``giac(...)`` to create a new Giac object.

    Full html documentation for giac is available from your giac installation at ``$PREFIX``/share/giac/doc/en/cascmd_en/index.html

    EXAMPLES:

    Any Giac instruction can be evaluated as a string by the giac command. You can access the giac functions by adding the ``giac.`` prefix to the usual Giac name.

    ::

      sage: l=giac('normal((y+sqrt(2))^4)'); l
      y^4+4*sqrt(2)*y^3+12*y^2+8*sqrt(2)*y+4
      sage: f=giac('(u,v)->{ if (u<v){ [u,v] } else { [v,u] }}');f(1,2),f(3,1)
      ([1,2], [1,3])

    The output of the giac command is a Giac object, and it can be used for another giac command.

    ::

      sage: l.factors()
      [y+sqrt(2),4]
      sage: giac('(x^12-1)').factor( )
      (x-1)*(x+1)*(x^2+1)*(x^2-x+1)*(x^2+x+1)*(x^4-x^2+1)
      sage: giac('assume(y>0)'); giac('y^2=3').solve('y')
      y
      ...[sqrt(3)]

    You can create some Giac elements and avoid many quotes like this:

    ::

      sage: x,y,z=giac('x,y,z');type(y)
      <class 'sage.interfaces.giac.GiacElement'>
      sage: I1=(1/(cos(2*y)+cos(y))).integral(y,0,pi/4).simplify()
      sage: (I1-((-2*ln((sqrt(3)-3*tan(1/8*pi))/(sqrt(3)+3*tan(1/8*pi)))*sqrt(3)-3*tan(1/8*pi))/9)).normal()
      0
      sage: ((y+z*sqrt(5))*(y-sqrt(5)*z)).normal()
      y^2-5*z^2

    Polynomials or elements of SR can be evaluated directly by the giac interface.

    ::

      sage: R.<a,b> = QQ[]; f = (2+a+b)
      sage: p = giac.gcd(f^3+5*f^5,f^2+f^5); p; R(p.sage())
      sageVARa^2+2*sageVARa*sageVARb+4*sageVARa+sageVARb^2+4*sageVARb+4
      a^2 + 2*a*b + b^2 + 4*a + 4*b + 4

    Variable names in python and giac are independent::

        sage: a=sqrt(2);giac('Digits:=30;a:=5');a,giac('a'),giac(a),giac(a).evalf()
        30
        (sqrt(2), 5, sqrt(2), 1.41421356237309504880168872421)

    TESTS::

        sage: g = giac('euler_gamma').sage();g
        euler_gamma
        sage: g.n()
        0.577215664901533
    """
    def __init__(self, maxread=None, script_subdirectory=None, server=None, server_tmpdir=None, logfile=None):
        """
        Create an instance of the Giac interpreter.

        EXAMPLES::

            sage: giac == loads(dumps(giac))
            True
        """
        Expect.__init__(self,
                        name = 'giac',
                        prompt = '[0-9]*>> ',
                        command = "giac --sage",
                        env = {"LANG": "C"},
                        init_code= ['maple_mode(0);I:=i;'],      #  coercion could be broken in maple_mode
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,                        server = server,
                        server_tmpdir = server_tmpdir,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=1000)

    def _function_class(self):
        """
        EXAMPLES::

            sage: giac._function_class()
            <class 'sage.interfaces.giac.GiacFunction'>

        ::

            sage: type(giac.diff)
            <class 'sage.interfaces.giac.GiacFunction'>
        """
        return GiacFunction

    def _keyboard_interrupt(self):
        """
        The pexepect interface for giac has a very poor support of keyboard interruptions.
        """
        print("Interrupting %s..." % self)
        self._expect.sendline(chr(3))  # send ctrl-c
        self._expect.expect(self._prompt)
#        self._expect.expect(self._prompt)
        raise RuntimeError("Ctrl-c pressed while running %s"%self)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: giac.__reduce__()
            (<function reduce_load_Giac at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            Giac
        """
        return reduce_load_Giac, tuple([])

    def _read_in_file_command(self, filename):
        r"""
        Returns the string used to read filename into Giac.

        EXAMPLES::

            sage: giac._read_in_file_command('test')
            'read "test"'

        ::

            sage: filename = tmp_filename()
            sage: with open(filename,'w') as f:
            ....:     _ = f.write('xx := 22;\n')
            sage: giac.read(filename)
            sage: giac.get('xx').strip()
            '22'
        """
        return 'read "%s"' % filename

    def _quit_string(self):
        """
        EXAMPLES::

            sage: giac._quit_string()
            '@d'

        ::

            sage: m = Giac()
            sage: a = m(2)
            sage: m.is_running()
            True
            sage: m.quit()
            sage: m.is_running()
            False
        """
        return '@d'

    def _install_hints(self):
        """
        Hints for installing Giac on your computer.

        EXAMPLES::

            sage: print(giac._install_hints())
            In order...
        """
        return r"""

In order to use the Giac interface you need to have Giac installed
and have a program called "giac" in your PATH. You need a giac version
supporting "giac --sage" ( roughly after 0.9.1 of march 2011). Some giac
instructions  and the help's langage depend of you LANG variable. To obtain
inline help for  giac commands, you also need to have the program "cas_help"
in your PATH.


If giac is not already installed, you can download binaries or sources
or a spkg ( for the spkg follow the sources link) from the homepage:

Homepage http://www-fourier.ujf-grenoble.fr/~parisse/giac.html


Full html documentation for giac is available from your giac installation at:

    ``$PREFIX``/share/giac/doc/en/cascmd_en/index.html

If you got giac from the spkg then ``$PREFIX`` is ``$SAGE_LOCAL``

"""

    def expect(self):
        """
        Returns the pexpect object for this Giac session.

        EXAMPLES::

            sage: m = Giac()
            sage: m.expect() is None
            True
            sage: m._start()
            sage: m.expect()
            Giac with PID ... running .../giac --sage
            sage: m.quit()
        """
        return self._expect

    def console(self):
        """
        Spawn a new Giac command-line session.

        EXAMPLES::

            sage: giac_console()                   # not tested - giac
            ...
            Homepage http://www-fourier.ujf-grenoble.fr/~parisse/giac.html
            Released under the GPL license 3.0 or above
            See http://www.gnu.org for license details
            -------------------------------------------------
            Press CTRL and D simultaneously to finish session
            Type ?commandname for help
            0>>

        """
        giac_console()


    def completions(self, s):
        """
        Return all commands that complete the command starting with the
        string s.

        EXAMPLES::

            sage: c = giac.completions('cas')
            sage: 'cas_setup' in c
            True
        """
        if self._expect is None:
            self._start()
        E = self._expect
        E.sendline('%s%s%s' % (s, chr(63), chr(13)))
        t = E.timeout
        E.timeout=0.3  # since some things have no completion
        try:
            E.expect('----')
        except pexpect.TIMEOUT:
            E.timeout = t
            return []
        E.timeout = t
        v = bytes_to_str(E.before)
        E.expect(self._prompt)
        E.expect(self._prompt)
        return v.split()[1:]

    def _commands(self):
        """
        Return list of all commands defined in Giac.

        EXAMPLES::

            sage: c = giac._commands()
            sage: len(c) > 100
            True
            sage: 'Psi' in c
            True
        """
        try:
            v = sum([self.completions(chr(65+n)) for n in range(26)], []) + \
                sum([self.completions(chr(97+n)) for n in range(26)], [])
        except RuntimeError:
            print("\n" * 3)
            print("*" * 70)
            print("WARNING: You do not have a working version of Giac installed!")
            print("*" * 70)
            v = []
        v.sort()
        return v

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        """
        Returns a list of all the commands defined in Giac and optionally
        (per default) store them to disk.

        EXAMPLES::

            sage: c = giac._tab_completion(use_disk_cache=False, verbose=False)
            sage: len(c) > 100
            True
            sage: 'factors' in c
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
                print("\nBuilding Giac command completion list (this takes")
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % COMMANDS_CACHE)
            v = self._commands()
            self.__tab_completion = v
            if len(v) > 200:
                # Giac is actually installed.
                sage.misc.persist.save(v, COMMANDS_CACHE)
            return v


    def cputime(self, t=None):
        r"""
        Returns the amount of CPU time that the Giac session has used. If
        ``t`` is not None, then it returns the difference
        between the current CPU time and ``t``.

        EXAMPLES::

            sage: t = giac.cputime()
            sage: t                     # random
            0.02
            sage: x = giac('x')
            sage: giac.diff(x^2, x)
            2*x
            sage: giac.cputime(t)       # random
            0.0
        """
        if t is None:
            return float(self('time()'))
        else:
            return float(self('time() - %s'%float(t)))

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: giac._eval_line('2+2')
            '4'

            sage: A = matrix([range(280)])
            sage: GA = giac(A)

        TESTS::

            sage: h1 = 'int(sin(x)^2, x)'
            sage: h2 = 'int(cos(x)^2, x)'
            sage: giac_result = giac(h1) + giac(h2)
            sage: bool(giac_result.sage() == x)
            True

        """
        with gc_disabled():
            z = Expect._eval_line(self, line, allow_use_file=allow_use_file,
                    wait_for_prompt=wait_for_prompt)
            if z.lower().find("error") != -1:
                raise RuntimeError("An error occurred running a Giac command:\nINPUT:\n%s\nOUTPUT:\n%s"%(line, z))
        lines = (line for line in z.splitlines()
                 if not line.startswith('Evaluation time:'))
        return "\n".join(lines)

    def eval(self, code, strip=True, **kwds):
        r"""
        Send the code x to the Giac interpreter.
        Remark: To enable multi-lines codes in the notebook magic mode: ``%giac``,
        the ``\n`` are removed before sending the code to giac.

        INPUT:

        - code -- str
        - strip -- Default is True and removes ``\n``

        EXAMPLES::

            sage: giac.eval("2+2;\n3")
            '4,3'
            sage: giac.eval("2+2;\n3",False)
            '4\n3'
            sage: s='g(x):={\nx+1;\nx+2;\n}'
            sage: giac(s)
            ...x+1...x+2...
            sage: giac.g(5)
            7
        """
        #we remove \n to enable multiline code in the notebook magic mode %giac
        if strip:
             code = code.replace("\n","").strip()
        ans = Expect.eval(self, code, strip=strip, **kwds).strip()
        return ans



    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: giac.set('xx', '2')
            sage: giac.get('xx')
            '2'
        """
        cmd = '%s:=%s:;'%(var,value)   #if giac is not in maple mode ( maple_mode(0))
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError("Error executing code in Giac\nCODE:\n\t%s\nGiac ERROR:\n\t%s"%(cmd, out))


    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: giac.set('xx', '2')
            sage: giac.get('xx')
            '2'
        """
        s = self.eval('%s'%var)
        return s

    def _object_class(self):
        """
        Returns the class of GiacElements.

        EXAMPLES::

            sage: giac._object_class()
            <class 'sage.interfaces.giac.GiacElement'>

        ::

            sage: m = giac(2)
            sage: type(m)
            <class 'sage.interfaces.giac.GiacElement'>
        """
        return GiacElement

    def _function_element_class(self):
        """
        Returns the GiacFunctionElement class.

        EXAMPLES::

            sage: giac._function_element_class()
            <class 'sage.interfaces.giac.GiacFunctionElement'>

        ::

            sage: two = giac(2)
            sage: type(two.gcd)
            <class 'sage.interfaces.giac.GiacFunctionElement'>
        """
        return GiacFunctionElement

    def _equality_symbol(self):
        """
        Returns the symbol used for equality testing in Giac.

        EXAMPLES::

            sage: giac._equality_symbol()
            '=='

            sage: giac(2) == giac(2)
            True
        """
        return '=='

    def _true_symbol(self):
        """
        Returns the symbol used for truth in Giac.

        EXAMPLES::

            sage: giac._true_symbol()
            'true'

        ::

            sage: giac(2) == giac(2)
            True
        """
        return 'true'

    def _assign_symbol(self):
        """
        Returns the symbol used for assignment in Giac.

        EXAMPLES::

            sage: giac._assign_symbol()
            ':='
        """
        return ":="

    def _help(self, string):
        r"""
        Return the Giac help on ``string``.

        EXAMPLES::

            sage: giac._help('gcd')  # not tested ; output may vary (LANG)
            "...gcd - greatest common divisor of polynomials...
        """
        return os.popen('cas_help %s' % string).read()
        # return os.popen('echo "?%s" | giac' % string).read()

    def help(self, string):
        """
        Display Giac help about string.

        This is the same as typing "?string" in the Giac console.

        INPUT:

        -  ``string`` -- a string to search for in the giac help system

        EXAMPLES::

            sage: giac.help('Psi')         # not tested - depends of giac and $LANG
            Psi(a,n)=nth-derivative of the function DiGamma (=ln@Gamma) at point a (Psi(a,0)=Psi(a))...
        """
        pager()(self._help(string))

    def clear(self, var):
        """
        Clear the variable named var.

        EXAMPLES::

            sage: giac.set('xx', '2')
            sage: giac.get('xx')
            '2'
            sage: giac.clear('xx')
            sage: giac.get('xx')
            'xx'
        """
        self.eval('purge(%s)'%var)

    def version(self):
        """
        Wrapper for giac's version().

        EXAMPLES::

            sage: giac.version()
            "giac...

        """
        return giac('version()')


@instancedoc
class GiacFunction(ExpectFunction):
    def _instancedoc_(self):
        """
        Returns the Giac help for this function. This gets called when
        doing "?" on self.

        EXAMPLES::

            sage: giac.gcd.__doc__  # random
            "gcd - greatest common divisor of polynomials...
        """
        M = self._parent
        return M._help(self._name)


@instancedoc
class GiacFunctionElement(FunctionElement):
    def _instancedoc_(self):
        """
        Returns the Giac help for this function. This gets called when
        doing "?" on self.

        EXAMPLES::

            sage: two = giac(2)
            sage: two.gcd.__doc__  # random
            "...gcd - greatest common divisor of polynomials...
        """
        return self._obj.parent()._help(self._name)


@instancedoc
class GiacElement(ExpectElement):
    def __float__(self):
        """
        Returns a floating point version of self.

        EXAMPLES::

            sage: float(giac(1/2))
            0.5
            sage: type(_)
            <class 'float'>
        """
        return float(giac.eval('evalf(%s)' % self.name()))

    def unapply(self, var):
        """
        Creates a Giac function in the given arguments from a Giac symbol.

        EXAMPLES::

            sage: f=giac('y^3+1+t')
            sage: g=(f.unapply('y,t'))
            sage: g
            (y,t)->y^3+1+t
            sage: g(1,2)
            4
        """
        return giac('unapply(%s,%s)'%(self,var))


    def __hash__(self):
        """
        Returns a  integer representing the hash of self.

        These examples are optional, and require Giac to be installed. You
        don't need to install any Sage packages for this.

        EXAMPLES::

            sage: m = giac('x^2+y^2')
            sage: hash(m)              # random
            4614285348919569149
        """
        return hash(giac.eval('string(%s);'%self.name()))

    def _richcmp_(self, other, op):
        """
        Compare equality between self and other, using giac.

        These examples are optional, and require Giac to be installed. You
        do not need to install any Sage packages for this.

        EXAMPLES::

            sage: a = giac(5)
            sage: b = giac(5)
            sage: a == b
            True
            sage: a == 5
            True

        ::

            sage: c = giac(3)
            sage: a == c
            False
            sage: a < c
            False
            sage: a < 6
            True
            sage: c <= a
            True

        ::

        TESTS::

            sage: x = var('x')
            sage: t = giac((x+1)^2)
            sage: u = giac(x^2+2*x+1)
            sage: u == t
            False
        """
        P = self.parent()
        if P.eval("evalb(%s %s %s)"%(self.name(), P._equality_symbol(),
                                 other.name())) == P._true_symbol():
            return rich_to_bool(op, 0)
        # (to be tested with giac). Maple  does not allow comparing objects of different types and
        # it raises an error in this case.
        # We catch the error, and return True for <
        try:
            if P.eval("evalb(%s %s %s)"%(self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
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
        if P.eval("evalb(%s %s %s)"%(self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
            return rich_to_bool(op, 1)

        return NotImplemented

    def _tab_completion(self):
        """
        EXAMPLES::

            sage: a = giac(2)
            sage: 'sin' in a._tab_completion()
            True
        """
        return self.parent()._tab_completion()

    def __len__(self):
        """
        EXAMPLES::

            sage: len(giac([1,2,3]))
            3
        """
        return int(self.size())

    def __iter__(self):
        """
        EXAMPLES::

            sage: l = giac([1,2,3])
            sage: list(iter(l))
            [1, 2, 3]
        """
        for i in range(len(self)):  # zero-indexed if giac is maple_mode(0)
            yield self[i]

    def __del__(self):
        """
        Note that clearing object is pointless since it wastes time.
        (Ex: otherwise doing a=0 after a = (giac('x+y+z')^40).normal() is very slow )

        EXAMPLES::

            sage: a = giac(2)
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

    def _latex_(self):
        r"""
        You can output Giac expressions in latex.

        EXAMPLES::

            sage: M = matrix(QQ, [[1, 2], [3, 4]])
            sage: latex(M)
            \left(\begin{array}{rr}
            1 & 2 \\
            3 & 4
            \end{array}\right)
            sage: gM = giac(M)
            sage: latex(gM)
            \left...\begin{array}{cc}...1...&...2...\\...3...&...4...\end{array}\right...
            sage: gf = giac('(x^4 - y)/(y^2-3*x)')
            sage: latex(gf)          # output changed slightly from 1.5.0-63 to 1.5.0-87
            \frac{...x^{4}...-...y...}{...y^{2}-3...x...}

        """
        s = self.parent().eval('latex(%s)'%self.name())
        if s.startswith('"'):
            s = s[1:]
        if s.endswith('"'):
            s = s[:-1]
        s = s.strip()
        return s

    def _matrix_(self, R):
        r"""
        Return matrix over the (Sage) ring R determined by self, where self
        should be a  Giac matrix.

        .. WARNING:: It is slow, do not convert big matrices.

        EXAMPLES::

            sage: R.<x,y>=QQ[]
            sage: M=giac('matrix(4,4,(k,l)->(x^k-y^l))'); M
            matrix[[0,1-y,1-y^2,1-y^3],[x-1,x-y,x-y^2,x-y^3],[x^2-1,x^2-y,x^2-y^2,x^2-y^3],[x^3-1,x^3-y,x^3-y^2,x^3-y^3]]
            sage: M.eigenvals()             # random
            0,0,(x^3+x^2+x-y^3-y^2-y+sqrt(x^6+2*x^5+3*x^4-14*x^3*y^3+2*x^3*y^2+2*x^3*y+6*x^3+2*x^2*y^3-14*x^2*y^2+2*x^2*y+5*x^2+2*x*y^3+2*x*y^2-14*x*y+4*x+y^6+2*y^5+3*y^4+6*y^3+5*y^2+4*y-12))/2,(x^3+x^2+x-y^3-y^2-y-sqrt(x^6+2*x^5+3*x^4-14*x^3*y^3+2*x^3*y^2+2*x^3*y+6*x^3+2*x^2*y^3-14*x^2*y^2+2*x^2*y+5*x^2+2*x*y^3+2*x*y^2-14*x*y+4*x+y^6+2*y^5+3*y^4+6*y^3+5*y^2+4*y-12))/2
            sage: Z=matrix(R,M);Z
            [         0     -y + 1   -y^2 + 1   -y^3 + 1]
            [     x - 1      x - y   -y^2 + x   -y^3 + x]
            [   x^2 - 1    x^2 - y  x^2 - y^2 -y^3 + x^2]
            [   x^3 - 1    x^3 - y  x^3 - y^2  x^3 - y^3]
            sage: parent(Z)
            Full MatrixSpace of 4 by 4 dense matrices over Multivariate Polynomial Ring in x, y over Rational Field
        """
        v = self.dim()
        n = int(v[0])
        m = int(v[1])

        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(R, n, m)
        entries = [[R(self[r, c]) for c in range(m)] for r in range(n)]
        return M(entries)

    def _sage_(self, locals={}):
        r"""
        Convert a giac expression back to a Sage expression, if possible.

        .. NOTE::

            This method works successfully when Giac returns a result
            or list of results that consist only of:
            - numbers, i.e. integers, floats, complex numbers;
            - functions and named constants also present in Sage, where:
                - Sage knows how to translate the function or constant's name
                from Giac's naming scheme through the ``symbol_table``, or
                - you provide a translation dictionary ``locals``.

        New conversions can be added using Pynac's ``register_symbol``.
        This is the recommended approach for library code.

        .. WARNING:: List conversion is slow.

        EXAMPLES::

            sage: m = giac('x^2 + 5*y')
            sage: m.sage()
            x^2 + 5*y

        ::

            sage: m = giac('sin(2*sqrt(1-x^2)) * (1 - cos(1/x))^2')
            sage: m.trigexpand().sage()
            2*cos(sqrt(-x^2 + 1))*cos(1/x)^2*sin(sqrt(-x^2 + 1)) - 4*cos(sqrt(-x^2 + 1))*cos(1/x)*sin(sqrt(-x^2 + 1)) + 2*cos(sqrt(-x^2 + 1))*sin(sqrt(-x^2 + 1))

        Converting a custom name using the ``locals`` dictionary::

            sage: ex = giac('myFun(x)')
            sage: ex._sage_({'myFun': sin})
            sin(x)

        Same but by adding a new entry to the ``symbol_table``::

            sage: ex = giac('myFun(x)')
            sage: sage.symbolic.expression.register_symbol(sin, {'giac':'myFun'})
            sage: ex._sage_()
            sin(x)

        Conversion of lists::

            sage: L = giac('solve((2/3)^x-2, x)'); L
            list[ln(2)/(ln(2)-ln(3))]
            sage: L.sage()
            [-log(2)/(log(3) - log(2))]

        TESTS:

        Check conversion of Booleans (:trac:`28705`)::

            sage: giac('true')._sage_(), giac('false')._sage_()
            (True, False)

        Check that variables and constants are not mixed up (:trac:`30133`)::

            sage: ee, ii, pp = SR.var('e,i,pi')
            sage: giac(ee * ii * pp).sage().variables()
            (e, i, pi)
            sage: giac(e * i * pi).sage().variables()
            ()
        """
        from sage.symbolic.expression import symbol_table
        from sage.calculus.calculus import symbolic_expression_from_string, SR_parser_giac

        result = repr(self) # string representation

        if str(self.type()) not in ['DOM_LIST', 'vector', 'vecteur']:

            # Merge the user-specified locals dictionary and the symbol_table
            # (locals takes priority)
            lsymbols = symbol_table['giac'].copy()
            lsymbols.update(locals)

            try:
                return symbolic_expression_from_string(result, lsymbols,
                    accept_sequence=True, parser=SR_parser_giac)

            except Exception:
                raise NotImplementedError("Unable to parse Giac output: %s" % result)
        else:
            return [entry.sage() for entry in self]

    def integral(self, var='x', min=None, max=None):
        r"""
        Return the integral of self with respect to the variable x.

        INPUT:


        -  ``var`` - variable

        -  ``min`` - default: None

        -  ``max`` - default: None


        Returns the definite integral if xmin is not None, otherwise
        returns an indefinite integral.

        EXAMPLES::

            sage: y=giac('y');f=(sin(2*y)/y).integral(y).simplify(); f
            Si(2*y)
            sage: f.diff(y).simplify()
            sin(2*y)/y

        ::

            sage: f = giac('exp(x^2)').integral('x',0,1) ; f
            1.46265174...
            sage: x,y=giac('x'),giac('y');integrate(cos(x+y),'x=0..pi').simplify()
            -2*sin(y)
        """
        if min is None:
            return giac('int(%s,%s)' % (self.name(), var))
        else:
            if max is None:
                raise ValueError("neither or both of min/max must be specified.")
        return giac('int(%s,%s,%s,%s)' % (self.name(), var,
                                          giac(min), giac(max)))

    integrate = integral

    def sum(self, var, min=None, max=None):
        r"""
        Return the sum of self with respect to the variable x.

        INPUT:

        -  ``var`` - variable

        -  ``min`` - default: None

        -  ``max`` - default: None

        Returns the definite integral if xmin is not None, otherwise
        returns an indefinite integral.

        EXAMPLES::

            sage: giac('1/(1+k^2)').sum('k',-oo,+infinity).simplify()
            (pi*exp(pi)^2+pi)/(exp(pi)^2-1)
        """
        if min is None:
            return giac('sum(%s,%s)' % (self.name(), var))
        else:
            if max is None:
                raise ValueError("neither or both of min/max must be specified.")
            return giac('sum(%s,%s,%s,%s)' % (self.name(), var,
                                              giac(min), giac(max)))


# An instance
giac = Giac()

def reduce_load_Giac():
    """
    Returns the giac object created in sage.interfaces.giac.

    EXAMPLES::

        sage: from sage.interfaces.giac import reduce_load_Giac
        sage: reduce_load_Giac()
        Giac
    """
    return giac


def giac_console():
    """
    Spawn a new Giac command-line session.

    EXAMPLES::

        sage: giac.console()  # not tested - giac
        ...
        Homepage http://www-fourier.ujf-grenoble.fr/~parisse/giac.html
        Released under the GPL license 3.0 or above
        See http://www.gnu.org for license details
        -------------------------------------------------
        Press CTRL and D simultaneously to finish session
        Type ?commandname for help
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%giac magics instead.')
    os.system('giac')


def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.giac import __doctest_cleanup
        sage: m = giac(2)
        sage: giac.is_running()
        True
        sage: __doctest_cleanup()
        sage: giac.is_running()
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
