r"""
Interface to Mathematica

The Mathematica interface will only work if Mathematica is installed on your
computer with a command line interface that runs when you give the ``math``
command. The interface lets you send certain Sage objects to Mathematica,
run Mathematica functions, import certain Mathematica expressions to Sage,
or any combination of the above.
The Sage command::

    sage: print(mathematica._install_hints())
    ...

prints more information on Mathematica installation.

To send a Sage object ``sobj`` to Mathematica, call ``mathematica(sobj)``.
This exports the Sage object to Mathematica and returns a new Sage object
wrapping the Mathematica expression/variable, so that you can use the
Mathematica variable from within Sage. You can then call Mathematica
functions on the new object; for example::

    sage: mobj = mathematica(x^2-1)             # optional - mathematica
    sage: mobj.Factor()                         # optional - mathematica
    (-1 + x)*(1 + x)

In the above example the factorization is done using Mathematica's
``Factor[]`` function.

To see Mathematica's output you can simply print the Mathematica wrapper
object. However if you want to import Mathematica's output back to Sage,
call the Mathematica wrapper object's ``sage()`` method. This method returns
a native Sage object::

    sage: # optional - mathematica
    sage: mobj = mathematica(x^2-1)
    sage: mobj2 = mobj.Factor(); mobj2
    (-1 + x)*(1 + x)
    sage: mobj2.parent()
    Mathematica
    sage: sobj = mobj2.sage(); sobj
    (x + 1)*(x - 1)
    sage: sobj.parent()
    Symbolic Ring


If you want to run a Mathematica function and don't already have the input
in the form of a Sage object, then it might be simpler to input a string to
``mathematica(expr)``. This string will be evaluated as if you had typed it
into Mathematica::

    sage: mathematica('Factor[x^2-1]')          # optional - mathematica
    (-1 + x)*(1 + x)
    sage: mathematica('Range[3]')               # optional - mathematica
    {1, 2, 3}

If you don't want Sage to go to the trouble of creating a wrapper for the
Mathematica expression, then you can call ``mathematica.eval(expr)``, which
returns the result as a Mathematica AsciiArtString formatted string. If you
want the result to be a string formatted like Mathematica's InputForm, call
``repr(mobj)`` on the wrapper object ``mobj``. If you want a string
formatted in Sage style, call ``mobj._sage_repr()``::

    sage: mathematica.eval('x^2 - 1')           # optional - mathematica
                   2
             -1 + x
    sage: repr(mathematica('Range[3]'))         # optional - mathematica
    '{1, 2, 3}'
    sage: mathematica('Range[3]')._sage_repr()  # optional - mathematica
    '[1, 2, 3]'

Finally, if you just want to use a Mathematica command line from within
Sage, the function ``mathematica_console()`` dumps you into an interactive
command-line Mathematica session. This is an enhanced version of the usual
Mathematica command-line, in that it provides readline editing and history
(the usual one doesn't!)

Tutorial
--------

We follow some of the tutorial from
http://library.wolfram.com/conferences/devconf99/withoff/Basic1.html/.

For any of this to work you must buy and install the Mathematica
program, and it must be available as the command
``math`` in your PATH.

Syntax
~~~~~~

Now make 1 and add it to itself. The result is a Mathematica
object.

::

    sage: m = mathematica
    sage: a = m(1) + m(1); a                # optional - mathematica
    2
    sage: a.parent()                        # optional - mathematica
    Mathematica
    sage: m('1+1')                          # optional - mathematica
    2
    sage: m(3)**m(50)                       # optional - mathematica
    717897987691852588770249

The following is equivalent to ``Plus[2, 3]`` in
Mathematica::

    sage: m = mathematica
    sage: m(2).Plus(m(3))                   # optional - mathematica
    5

We can also compute `7(2+3)`.

::

    sage: m(7).Times(m(2).Plus(m(3)))       # optional - mathematica
    35
    sage: m('7(2+3)')                       # optional - mathematica
    35

Some typical input
~~~~~~~~~~~~~~~~~~

We solve an equation and a system of two equations::

    sage: # optional - mathematica
    sage: eqn = mathematica('3x + 5 == 14')
    sage: eqn
    5 + 3*x == 14
    sage: eqn.Solve('x')
    {{x -> 3}}
    sage: sys = mathematica('{x^2 - 3y == 3, 2x - y == 1}')
    sage: print(sys)
               2
             {x  - 3 y == 3, 2 x - y == 1}
    sage: sys.Solve('{x, y}')
    {{x -> 0, y -> -1}, {x -> 6, y -> 11}}

Assignments and definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you assign the mathematica `5` to a variable `c`
in Sage, this does not affect the `c` in Mathematica.

::

    sage: c = m(5)                          # optional - mathematica
    sage: print(m('b + c x'))               # optional - mathematica
                 b + c x
    sage: print(m('b') + c*m('x'))          # optional - mathematica
             b + 5 x

The Sage interfaces changes Sage lists into Mathematica lists::

    sage: m = mathematica
    sage: eq1 = m('x^2 - 3y == 3')          # optional - mathematica
    sage: eq2 = m('2x - y == 1')            # optional - mathematica
    sage: v = m([eq1, eq2]); v              # optional - mathematica
    {x^2 - 3*y == 3, 2*x - y == 1}
    sage: v.Solve(['x', 'y'])               # optional - mathematica
    {{x -> 0, y -> -1}, {x -> 6, y -> 11}}

Function definitions
~~~~~~~~~~~~~~~~~~~~

Define mathematica functions by simply sending the definition to
the interpreter.

::

    sage: m = mathematica
    sage: _ = mathematica('f[p_] = p^2');   # optional - mathematica
    sage: m('f[9]')                         # optional - mathematica
    81

Numerical Calculations
~~~~~~~~~~~~~~~~~~~~~~

We find the `x` such that `e^x - 3x = 0`.

::

    sage: eqn = mathematica('Exp[x] - 3x == 0') # optional - mathematica
    sage: eqn.FindRoot(['x', 2])                # optional - mathematica
    {x -> 1.512134551657842}

Note that this agrees with what the PARI interpreter gp produces::

    sage: gp('solve(x=1,2,exp(x)-3*x)')
    1.512134551657842473896739678              # 32-bit
    1.5121345516578424738967396780720387046    # 64-bit

Next we find the minimum of a polynomial using the two different
ways of accessing Mathematica::

    sage: mathematica('FindMinimum[x^3 - 6x^2 + 11x - 5, {x,3}]')  # optional - mathematica
    {0.6150998205402516, {x -> 2.5773502699629733}}
    sage: f = mathematica('x^3 - 6x^2 + 11x - 5')  # optional - mathematica
    sage: f.FindMinimum(['x', 3])                  # optional - mathematica
    {0.6150998205402516, {x -> 2.5773502699629733}}

Polynomial and Integer Factorization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We factor a polynomial of degree 200 over the integers.

::

    sage: R.<x> = PolynomialRing(ZZ)
    sage: f = (x**100+17*x+5)*(x**100-5*x+20)
    sage: f
    x^200 + 12*x^101 + 25*x^100 - 85*x^2 + 315*x + 100
    sage: g = mathematica(str(f))            # optional - mathematica
    sage: print(g)                           # optional - mathematica
                               2       100       101    200
             100 + 315 x - 85 x  + 25 x    + 12 x    + x
    sage: g                                  # optional - mathematica
    100 + 315*x - 85*x^2 + 25*x^100 + 12*x^101 + x^200
    sage: print(g.Factor())                  # optional - mathematica
                          100               100
             (20 - 5 x + x   ) (5 + 17 x + x   )

We can also factor a multivariate polynomial::

    sage: f = mathematica('x^6 + (-y - 2)*x^5 + (y^3 + 2*y)*x^4 - y^4*x^3')  # optional - mathematica
    sage: print(f.Factor())                  # optional - mathematica
              3                  2    3
             x  (x - y) (-2 x + x  + y )

We factor an integer::

    sage: # optional - mathematica
    sage: n = mathematica(2434500)
    sage: n.FactorInteger()
    {{2, 2}, {3, 2}, {5, 3}, {541, 1}}
    sage: n = mathematica(2434500)
    sage: F = n.FactorInteger(); F
    {{2, 2}, {3, 2}, {5, 3}, {541, 1}}
    sage: F[1]
    {2, 2}
    sage: F[4]
    {541, 1}

Mathematica's ECM package is no longer available.

Long Input
----------

The Mathematica interface reads in even very long input (using
files) in a robust manner.

::

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = mathematica(t)        # optional - mathematica
    sage: a = mathematica.eval(t)   # optional - mathematica

Loading and saving
------------------

Mathematica has an excellent ``InputForm`` function,
which makes saving and loading Mathematica objects possible. The
first examples test saving and loading to strings.

::

    sage: # optional - mathematica
    sage: x = mathematica(pi/2)
    sage: print(x)
             Pi
             --
             2
    sage: loads(dumps(x)) == x
    True
    sage: n = x.N(50)
    sage: print(n)
                  1.5707963267948966192313216916397514420985846996876
    sage: loads(dumps(n)) == n
    True

Complicated translations
------------------------

The ``mobj.sage()`` method tries to convert a Mathematica object to a Sage
object. In many cases, it will just work. In particular, it should be able to
convert expressions entirely consisting of:

- numbers, i.e. integers, floats, complex numbers;
- functions and named constants also present in Sage, where:

    - Sage knows how to translate the function or constant's name from
      Mathematica's, or
    - the Sage name for the function or constant is trivially related to
      Mathematica's;

- symbolic variables whose names don't pathologically overlap with
  objects already defined in Sage.

This method will not work when Mathematica's output includes:

- strings;
- functions unknown to Sage;
- Mathematica functions with different parameters/parameter order to
  the Sage equivalent.

If you want to convert more complicated Mathematica expressions, you can
instead call ``mobj._sage_()`` and supply a translation dictionary::

    sage: m = mathematica('NewFn[x]')       # optional - mathematica
    sage: m._sage_(locals={('NewFn', 1): sin})   # optional - mathematica
    sin(x)

For more details, see the documentation for ``._sage_()``.


OTHER Examples::

    sage: def math_bessel_K(nu,x):
    ....:     return mathematica(nu).BesselK(x).N(20)
    sage: math_bessel_K(2,I)                      # optional - mathematica
    -2.59288617549119697817 + 0.18048997206696202663*I

::

    sage: slist = [[1, 2], 3., 4 + I]
    sage: mlist = mathematica(slist); mlist     # optional - mathematica
    {{1, 2}, 3., 4 + I}
    sage: slist2 = list(mlist); slist2          # optional - mathematica
    [{1, 2}, 3., 4 + I]
    sage: slist2[0]                             # optional - mathematica
    {1, 2}
    sage: slist2[0].parent()                    # optional - mathematica
    Mathematica
    sage: slist3 = mlist.sage(); slist3         # optional - mathematica
    [[1, 2], 3.00000000000000, I + 4]

::

    sage: mathematica('10.^80')     # optional - mathematica
    1.*^80
    sage: mathematica('10.^80').sage()  # optional - mathematica
    1.00000000000000e80

AUTHORS:

- William Stein (2005): first version

- Doug Cutrell (2006-03-01): Instructions for use under Cygwin/Windows.

- Felix Lawrence (2009-08-21): Added support for importing Mathematica lists
  and floats with exponents.

TESTS:

Check that numerical approximations via Mathematica's `N[]` function work
correctly (:issue:`18888`, :issue:`28907`)::

    sage: # optional - mathematica
    sage: mathematica('Pi/2').N(10)
    1.5707963268
    sage: mathematica('Pi').N(10)
    3.1415926536
    sage: mathematica('Pi').N(50)
    3.14159265358979323846264338327950288419716939937511
    sage: str(mathematica('Pi*x^2-1/2').N())
                    2
    -0.5 + 3.14159 x

Check that Mathematica's `E` exponential symbol is correctly backtranslated
as Sage's `e` (:issue:`29833`)::

    sage: x = var('x')
    sage: (e^x)._mathematica_().sage()  # optional -- mathematica
    e^x
    sage: exp(x)._mathematica_().sage() # optional -- mathematica
    e^x

Check that all trig/hyperbolic functions and their reciprocals are correctly
translated to Mathematica (:issue:`34087`)::

    sage: # optional - mathematica
    sage: x=var('x')
    sage: FL=[sin, cos, tan, csc, sec, cot,
    ....:     sinh, cosh, tanh, csch, sech, coth]
    sage: IFL=[arcsin, arccos, arctan, arccsc,
    ....:      arcsec, arccot, arcsinh, arccosh,
    ....:      arctanh, arccsch, arcsech, arccoth]
    sage: [mathematica.TrigToExp(u(x)).sage()
    ....:  for u in FL]
    [-1/2*I*e^(I*x) + 1/2*I*e^(-I*x),
     1/2*e^(I*x) + 1/2*e^(-I*x),
     (-I*e^(I*x) + I*e^(-I*x))/(e^(I*x) + e^(-I*x)),
     2*I/(e^(I*x) - e^(-I*x)),
     2/(e^(I*x) + e^(-I*x)),
     -(-I*e^(I*x) - I*e^(-I*x))/(e^(I*x) - e^(-I*x)),
     -1/2*e^(-x) + 1/2*e^x,
     1/2*e^(-x) + 1/2*e^x,
     -e^(-x)/(e^(-x) + e^x) + e^x/(e^(-x) + e^x),
     -2/(e^(-x) - e^x),
     2/(e^(-x) + e^x),
     -(e^(-x) + e^x)/(e^(-x) - e^x)]
    sage: [mathematica.TrigToExp(u(x)).sage()
    ....:  for u in IFL]
    [-I*log(I*x + sqrt(-x^2 + 1)),
     1/2*pi + I*log(I*x + sqrt(-x^2 + 1)),
     -1/2*I*log(I*x + 1) + 1/2*I*log(-I*x + 1),
     -I*log(sqrt(-1/x^2 + 1) + I/x),
     1/2*pi + I*log(sqrt(-1/x^2 + 1) + I/x),
     -1/2*I*log(I/x + 1) + 1/2*I*log(-I/x + 1),
     log(x + sqrt(x^2 + 1)),
     log(sqrt(x + 1)*sqrt(x - 1) + x),
     1/2*log(x + 1) - 1/2*log(-x + 1),
     log(sqrt(1/x^2 + 1) + 1/x),
     log(sqrt(1/x + 1)*sqrt(1/x - 1) + 1/x),
     1/2*log(1/x + 1) - 1/2*log(-1/x + 1)]
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

import os
import re

from sage.misc.cachefunc import cached_method
from sage.interfaces.expect import (Expect, ExpectElement, ExpectFunction,
                                    FunctionElement)
from sage.interfaces.interface import AsciiArtString
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.misc.instancedoc import instancedoc
from sage.structure.richcmp import rich_to_bool


def clean_output(s):
    if s is None:
        return ''
    i = s.find('Out[')
    j = i + s[i:].find('=')
    s = s[:i] + ' ' * (j + 1 - i) + s[j + 1:]
    s = s.replace('\\\n', '')
    return s.strip('\n')


def _un_camel(name):
    """
    Convert `CamelCase` to `camel_case`.

    EXAMPLES::

        sage: sage.interfaces.mathematica._un_camel('CamelCase')
        'camel_case'
        sage: sage.interfaces.mathematica._un_camel('EllipticE')
        'elliptic_e'
        sage: sage.interfaces.mathematica._un_camel('FindRoot')
        'find_root'
        sage: sage.interfaces.mathematica._un_camel('GCD')
        'gcd'
    """
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


class Mathematica(ExtraTabCompletion, Expect):
    """
    Interface to the Mathematica interpreter.
    """
    def __init__(self, maxread=None, script_subdirectory=None, logfile=None, server=None,
                 server_tmpdir=None, command=None, verbose_start=False):
        r"""
        TESTS:

        Test that :issue:`28075` is fixed::

            sage: repr(mathematica.eval("Print[1]; Print[2]; Print[3]"))  # optional - mathematica
            '1\n2\n3'
        """
        # We use -rawterm to get a raw text interface in Mathematica 9 or later.
        # This works around the following issues of Mathematica 9 or later
        # (tested with Mathematica 11.0.1 for Mac OS X x86 (64-bit))
        #
        # 1) If TERM is unset and input is a pseudoterminal, Mathematica shows no
        # prompts, so pexpect will not work.
        #
        # 2) If TERM is set (to dumb, lpr, vt100, or xterm), there will be
        # prompts; but there is bizarre echoing behavior by Mathematica (not
        # the terminal driver).  For example, with TERM=dumb, many spaces and
        # \r's are echoed.  With TERM=vt100 or better, in addition, many escape
        # sequences are printed.
        #
        if command is None:
            command = os.getenv('SAGE_MATHEMATICA_COMMAND') or 'math -rawterm'
        eval_using_file_cutoff = 1024
        # Removing terminal echo using "stty -echo" is not essential but it slightly
        # improves performance (system time) and eliminates races of the terminal echo
        # as a possible source of error.
        if server:
            command = 'stty -echo; {}'.format(command)
        else:
            command = 'sh -c "stty -echo; {}"'.format(command)
        Expect.__init__(self,
                        name='mathematica',
                        terminal_echo=False,
                        command=command,
                        prompt=r'In\[[0-9]+\]:= ',
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        verbose_start=verbose_start,
                        logfile=logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)

    def _read_in_file_command(self, filename):
        return '<<"%s"' % filename

    def _keyboard_interrupt(self):
        print("Interrupting %s..." % self)
        e = self._expect
        e.sendline(chr(3))  # send ctrl-c
        e.expect('Interrupt> ')
        e.sendline("a")  # a -- abort
        e.expect(self._prompt)
        return e.before

    def _install_hints(self):
        """
        Hints for installing mathematica on your computer.

        AUTHORS:

        - William Stein and Justin Walker (2006-02-12)
        """
        return """
In order to use the Mathematica interface you need to have Mathematica
installed and have a script in your PATH called "math" that runs the
command-line version of Mathematica. Alternatively, you could use a
remote connection to a server running Mathematica -- for hints, type
    print(mathematica._install_hints_ssh())


  (1) You might have to buy Mathematica (https://www.wolfram.com/), or
  install a currently (Feb 2022) free for personal use Wolfram Engine
  (https://www.wolfram.com/engine/).

  (2) * LINUX: The math script usually comes standard with your Mathematica install.
        However, on some systems it may be called wolfram,
        or, in case of Wolfram Engine, wolframengine, while math is absent.
        In this case, assuming wolfram, respectively, wolframengine,
        is in your PATH,
          (a) create a file called math (in your PATH):
              #!/bin/sh
              /usr/bin/env wolfram $@

        respectively,
          (a') create a file called math (in your PATH):
              #!/bin/sh
              /usr/bin/env wolframengine $@

          (b) Make the file executable.
                chmod +x math

      * Apple macOS: for Mathematica,
          (a) create a file called math (in your PATH):
              #!/bin/sh
              /Applications/Mathematica.app/Contents/MacOS/MathKernel $@

          (a') for Wolfram Engine, follow the Linux step (a') above.

          The path in the above script must be modified if you installed
          Mathematica elsewhere or installed an old version of
          Mathematica that has the version in the .app name.

          (b) Make the file executable.
                chmod +x math

      * WINDOWS:

        Install Mathematica for Linux into the VMware virtual machine, or in
        a WSL/WSL2 Linux installation with Sage installed there (sorry,
        that's the only ways at present).
"""

    def eval(self, code, strip=True, **kwds):
        s = Expect.eval(self, code, **kwds)
        if strip:
            return AsciiArtString(clean_output(s))
        else:
            return AsciiArtString(s)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s;' % (var, value)
        out = self._eval_line(cmd, allow_use_file=True)
        if len(out) > 8:
            raise TypeError("Error executing code in Mathematica\nCODE:\n\t%s\nMathematica ERROR:\n\t%s" % (cmd, out))

    def get(self, var, ascii_art=False):
        """
        Get the value of the variable var.

        AUTHORS:

        - William Stein

        - Kiran Kedlaya (2006-02-04): suggested using InputForm
        """
        if ascii_art:
            return self.eval(var, strip=True)
        return self.eval('InputForm[%s, NumberMarks->False]' % var, strip=True)

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=False):
        s = Expect._eval_line(self, line,
             allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt)
        return str(s).strip('\n')

    def _function_call_string(self, function, args, kwds):
        """
        Return the string used to make function calls.

        EXAMPLES::

            sage: mathematica._function_call_string('Sin', ['x'], [])
            'Sin[x]'
        """
        return "%s[%s]" % (function, ",".join(args))

    def _left_list_delim(self):
        return "{"

    def _right_list_delim(self):
        return "}"

    def _left_func_delim(self):
        return "["

    def _right_func_delim(self):
        return "]"

    ###########################################
    # System -- change directory, etc.
    ###########################################
    def chdir(self, dir):
        """
        Change Mathematica's current working directory.

        EXAMPLES::

            sage: mathematica.chdir('/')          # optional - mathematica
            sage: mathematica('Directory[]')      # optional - mathematica
            "/"
        """
        self.eval('SetDirectory["%s"]' % dir)

    def _true_symbol(self):
        return 'True'

    def _false_symbol(self):
        return 'False'

    def _equality_symbol(self):
        return '=='

    def _assign_symbol(self):
        return ":="

    def _exponent_symbol(self):
        """
        Return the symbol used to denote the exponent of a number in
        Mathematica.

        EXAMPLES::

            sage: mathematica._exponent_symbol()        # optional - mathematica
            '*^'

        ::

            sage: bignum = mathematica('10.^80')  # optional - mathematica
            sage: repr(bignum)                          # optional - mathematica
            '1.*^80'
            sage: repr(bignum).replace(mathematica._exponent_symbol(), 'e').strip() # optional - mathematica
            '1.e80'
        """
        return "*^"

    def _object_class(self):
        return MathematicaElement

    def console(self, readline=True):
        mathematica_console(readline=readline)

    def _tab_completion(self):
        a = self.eval('Names["*"]')
        return a.replace('$', '').replace('\n \n>', '').replace(',', '').replace('}', '').replace('{', '').split()

    def help(self, cmd):
        return self.eval('? %s' % cmd)

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MathematicaFunction(self, attrname)


@instancedoc
class MathematicaElement(ExpectElement):
    def __getitem__(self, n):
        return self.parent().new('%s[[%s]]' % (self._name, n))

    def __getattr__(self, attrname):
        self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        return MathematicaFunctionElement(self, attrname)

    def __float__(self, precision=16):
        P = self.parent()
        return float(P.eval('N[%s,%s]' % (self.name(), precision)))

    def _reduce(self):
        return self.parent().eval('InputForm[%s]' % self.name()).strip()

    def __reduce__(self):
        return reduce_load, (self._reduce(), )

    def _latex_(self):
        z = self.parent().eval('TeXForm[%s]' % self.name())
        i = z.find('=')
        return z[i + 1:].strip()

    def _repr_(self):
        P = self.parent()
        return P.get(self._name, ascii_art=False).strip()

    def _sage_(self, locals={}):
        r"""
        Attempt to return a Sage version of this object.

        This method works successfully when Mathematica returns a result
        or list of results that consist only of:

        - numbers, i.e. integers, floats, complex numbers;
        - functions and named constants also present in Sage, where:
            - Sage knows how to translate the function or constant's name
              from Mathematica's naming scheme, or
            - you provide a translation dictionary `locals`, or
            - the Sage name for the function or constant is simply the
              Mathematica name in lower case;

        - symbolic variables whose names do not pathologically overlap with
          objects already defined in Sage.

        This method will not work when Mathematica's output includes:

        - strings;
        - functions unknown to Sage that are not specified in `locals`;
        - Mathematica functions with different parameters/parameter order to
          the Sage equivalent. In this case, define a function to do the
          parameter conversion, and pass it in via the locals dictionary.

        EXAMPLES:

        Mathematica lists of numbers/constants become Sage lists of
        numbers/constants::

            sage: # optional - mathematica
            sage: m = mathematica('{{1., 4}, Pi, 3.2e100, I}')
            sage: s = m.sage(); s
            [[1.00000000000000, 4], pi, 3.20000000000000*e100, I]
            sage: s[1].n()
            3.14159265358979
            sage: s[3]^2
            -1

        ::

            sage: m = mathematica('x^2 + 5*y')      # optional - mathematica
            sage: m.sage()                          # optional - mathematica
            x^2 + 5*y

        ::

            sage: m = mathematica('Sin[Sqrt[1-x^2]] * (1 - Cos[1/x])^2')  # optional - mathematica
            sage: m.sage()                          # optional - mathematica
            (cos(1/x) - 1)^2*sin(sqrt(-x^2 + 1))

        ::

            sage: m = mathematica('NewFn[x]')       # optional - mathematica
            sage: m._sage_(locals={('NewFn', 1): sin})   # optional - mathematica
            sin(x)

        ::

            sage: var('bla')                        # optional - mathematica
            bla
            sage: m = mathematica('bla^2')          # optional - mathematica
            sage: bla^2 - m.sage()                  # optional - mathematica
            0

        ::

            sage: # optional - mathematica
            sage: m = mathematica('bla^2')
            sage: mb = m.sage()
            sage: var('bla')
            bla
            sage: bla^2 - mb
            0


        AUTHORS:

        - Felix Lawrence (2010-11-03): Major rewrite to use ._sage_repr() and
          sage.calculus.calculus.symbolic_expression_from_string() for greater
          compatibility, while still supporting conversion of symbolic
          expressions.

        TESTS:

        Check that :issue:`28814` is fixed::

            sage: mathematica('Exp[1000.0]').sage()  # optional - mathematica
            1.97007111401700e434
            sage: mathematica('1/Exp[1000.0]').sage()  # optional - mathematica
            5.07595889754950e-435
            sage: mathematica(RealField(100)(1/3)).sage()  # optional - mathematica
            0.3333333333333333333333333333335
        """
        from sage.symbolic.expression import symbol_table
        from sage.symbolic.constants import constants_name_table as constants
        from sage.calculus.calculus import symbolic_expression_from_string
        from sage.calculus.calculus import _find_func as find_func

        # Get Mathematica's output and perform preliminary formatting
        res = self._sage_repr()
        if '"' in res:
            raise NotImplementedError("String conversion from Mathematica \
                does not work.  Mathematica's output was: %s" % res)

        # Find all the mathematica functions, constants and symbolic variables
        # present in `res`.  Convert MMA functions and constants to their
        # Sage equivalents (if possible), using `locals` and
        # `sage.symbolic.pynac.symbol_table['mathematica']` as translation
        # dictionaries.  If a MMA function or constant is not in either
        # dictionary, then we use a variety of tactics listed in `autotrans`.
        # If a MMA variable is not in any dictionary, then create an
        # identically named Sage equivalent.

        # Merge the user-specified locals dictionary and the symbol_table
        # (locals takes priority)
        lsymbols = symbol_table['mathematica'].copy()
        lsymbols.update(locals)

        # Strategies for translating unknown functions/constants:
        autotrans = [str.lower,      # Try it in lower case
                     _un_camel,    # Convert `CamelCase` to `camel_case`
                     lambda x: x]     # Try the original name

        # Find the MMA funcs/vars/constants - they start with a letter.
        # Exclude exponents (e.g. 'e8' from 4.e8)
        p = re.compile(r'(?<!\.)[a-zA-Z]\w*')
        for m in p.finditer(res):
            # If the function, variable or constant is already in the
            # translation dictionary, then just move on.
            if m.group() in lsymbols:
                pass
            # Now try to translate all other functions -- try each strategy
            # in `autotrans` and check if the function exists in Sage
            elif m.end() < len(res) and res[m.end()] == '(':
                for t in autotrans:
                    f = find_func(t(m.group()), create_when_missing=False)
                    if f is not None:
                        lsymbols[m.group()] = f
                        break
                else:
                    raise NotImplementedError("Don't know a Sage equivalent \
                        for Mathematica function '%s'.  Please specify one \
                        manually using the 'locals' dictionary" % m.group())
            # Check if Sage has an equivalent constant
            else:
                for t in autotrans:
                    if t(m.group()) in constants:
                        lsymbols[m.group()] = constants[t(m.group())]
                        break
            # If Sage has never heard of the variable, then
            # symbolic_expression_from_string will automatically create it
        try:
            return symbolic_expression_from_string(res, lsymbols,
                accept_sequence=True)
        except Exception:
            raise NotImplementedError("Unable to parse Mathematica \
                output: %s" % res)

    def __str__(self):
        P = self._check_valid()
        return P.get(self._name, ascii_art=True)

    def __len__(self):
        """
        Return the object's length, evaluated by mathematica.

        EXAMPLES::

            sage: len(mathematica([1,1.,2]))    # optional - mathematica
            3

        AUTHORS:
        - Felix Lawrence (2009-08-21)
        """
        return int(self.Length())

    @cached_method
    def _is_graphics(self):
        """
        Test whether the mathematica expression is graphics.

        OUTPUT: boolean

        EXAMPLES::

            sage: P = mathematica('Plot[Sin[x],{x,-2Pi,4Pi}]')   # optional - mathematica
            sage: P._is_graphics()                               # optional - mathematica
            True
        """
        P = self._check_valid()
        return P.eval('InputForm[%s]' % self.name()).strip().startswith('Graphics[')

    def save_image(self, filename, ImageSize=600):
        r"""
        Save a mathematica graphics.

        INPUT:

        - ``filename`` -- string; the filename to save as. The
          extension determines the image file format

        - ``ImageSize`` -- integer; the size of the resulting image

        EXAMPLES::

            sage: P = mathematica('Plot[Sin[x],{x,-2Pi,4Pi}]')   # optional - mathematica
            sage: filename = tmp_filename()                      # optional - mathematica
            sage: P.save_image(filename, ImageSize=800)                # optional - mathematica
        """
        P = self._check_valid()
        if not self._is_graphics():
            raise ValueError('mathematica expression is not graphics')
        filename = os.path.abspath(filename)
        s = 'Export["%s", %s, ImageSize->%s]' % (filename, self.name(),
                                                 ImageSize)
        P.eval(s)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: P = mathematica('Plot[Sin[x],{x,-2Pi,4Pi}]')   # optional - mathematica

        The following test requires a working X display on Linux so that the
        Mathematica frontend can do the rendering (:issue:`23112`)::

            sage: P._rich_repr_(dm)                              # optional - mathematica mathematicafrontend
            OutputImagePng container
        """
        if self._is_graphics():
            OutputImagePng = display_manager.types.OutputImagePng
            if display_manager.preferences.graphics == 'disable':
                return
            if OutputImagePng in display_manager.supported_output():
                return display_manager.graphics_from_save(
                    self.save_image, kwds, '.png', OutputImagePng)
        else:
            OutputLatex = display_manager.types.OutputLatex
            dmp = display_manager.preferences.text
            if dmp is None or dmp == 'plain':
                return
            if dmp == 'latex' and OutputLatex in display_manager.supported_output():
                return OutputLatex(self._latex_())

    def show(self, ImageSize=600):
        r"""
        Show a mathematica expression immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        INPUT:

        - ``ImageSize`` -- integer; the size of the resulting image

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image.

        EXAMPLES::

            sage: Q = mathematica('Sin[x Cos[y]]/Sqrt[1-x^2]')   # optional - mathematica
            sage: show(Q)                                        # optional - mathematica
            Sin[x*Cos[y]]/Sqrt[1 - x^2]

        The following example starts a Mathematica frontend to do the rendering
        (:issue:`28819`)::

            sage: P = mathematica('Plot[Sin[x],{x,-2Pi,4Pi}]')   # optional - mathematica
            sage: show(P)                                        # optional - mathematica mathematicafrontend
            sage: P.show(ImageSize=800)                          # optional - mathematica mathematicafrontend
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, ImageSize=ImageSize)

    def str(self):
        return str(self)

    def _richcmp_(self, other, op):
        P = self.parent()
        if P.eval("%s < %s" % (self.name(), other.name())).strip() == 'True':
            return rich_to_bool(op, -1)
        elif P.eval("%s > %s" % (self.name(), other.name())).strip() == 'True':
            return rich_to_bool(op, 1)
        elif P.eval("%s == %s" % (self.name(), other.name())).strip() == 'True':
            return rich_to_bool(op, 0)
        return NotImplemented

    def __bool__(self):
        """
        Return whether this Mathematica element is not identical to ``False``.

        EXAMPLES::

            sage: bool(mathematica(True))  # optional - mathematica
            True
            sage: bool(mathematica(False))  # optional - mathematica
            False

        In Mathematica, `0` cannot be used to express falsity::

            sage: bool(mathematica(0))  # optional - mathematica
            True
        """
        P = self._check_valid()
        cmd = '%s===%s' % (self._name, P._false_symbol())
        return P.eval(cmd).strip() != P._true_symbol()

    def n(self, *args, **kwargs):
        r"""
        Numerical approximation by converting to Sage object first.

        Convert the object into a Sage object and return its numerical
        approximation. See documentation of the function
        :func:`sage.misc.functional.n` for details.

        EXAMPLES::

            sage: mathematica('Pi').n(10)    # optional -- mathematica
            3.1
            sage: mathematica('Pi').n()      # optional -- mathematica
            3.14159265358979
            sage: mathematica('Pi').n(digits=10)   # optional -- mathematica
            3.141592654
        """
        return self._sage_().n(*args, **kwargs)


@instancedoc
class MathematicaFunction(ExpectFunction):
    def _instancedoc_(self):
        M = self._parent
        return M.help(self._name)


@instancedoc
class MathematicaFunctionElement(FunctionElement):
    def _instancedoc_(self):
        M = self._obj.parent()
        return M.help(self._name)


# An instance
mathematica = Mathematica()


def reduce_load(X):
    return mathematica(X)


def mathematica_console(readline=True):
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%mathematica magics instead.')
    if not readline:
        os.system('math')
        return
    else:
        os.system('math-readline')
        return


# some tools for online interface

def request_wolfram_alpha(input, verbose=False):
    r"""
    Request Wolfram Alpha website.

    INPUT:

    - ``input`` -- string
    - ``verbose`` -- boolean (default: ``False``)

    OUTPUT: json

    EXAMPLES::

        sage: from sage.interfaces.mathematica import request_wolfram_alpha
        sage: page_data = request_wolfram_alpha('integrate Sin[x]')      # optional internet
        sage: [str(a) for a in sorted(page_data.keys())]                 # optional internet
        ['queryresult']
        sage: [str(a) for a in sorted(page_data['queryresult'].keys())]  # optional internet
        ['datatypes',
         'encryptedEvaluatedExpression',
         'encryptedParsedExpression',
         'error',
         'host',
         'id',
         'inputstring',
         'numpods',
         'parsetimedout',
         'parsetiming',
         'pods',
         'recalculate',
         'related',
         'server',
         'sponsorCategories',
         'success',
         'timedout',
         'timedoutpods',
         'timing',
         'version']
    """
    from urllib.parse import urlencode
    from urllib.request import Request, build_opener, HTTPCookieProcessor, HTTPSHandler
    import json
    from http.cookiejar import CookieJar
    from ssl import create_default_context as default_context

    # we need cookies for this...
    cj = CookieJar()
    opener = build_opener(HTTPCookieProcessor(cj),
                          HTTPSHandler(context=default_context()))
    # build initial query for code
    req = Request("https://www.wolframalpha.com/input/api/v1/code")
    resp = opener.open(req)
    # the website returns JSON containing the code
    page_data = json.loads(resp.read().decode("utf-8"))
    if "code" not in page_data:
        raise ValueError("Wolfram did not return a code")
    proxy_code = page_data['code']
    if verbose:
        print("Code: {}".format(proxy_code))
        print("Cookies: {}".format(cj))
    # now we can make a request
    # some parameters documented here:
    #   https://products.wolframalpha.com/api/documentation/#parameter-reference
    # the following are the parameters used by the website
    params = {
        'assumptionsversion': '2',
        'async': 'true',
        'banners': 'raw',
        'debuggingdata': 'false',
        'format': 'image,plaintext,imagemap,sound,minput,moutput',
        'formattimeout': '8',
        'input': input,
        'output': 'JSON',
        'parsetimeout': '5',
        'podinfosasync': 'true',
        'proxycode': proxy_code,
        'recalcscheme': 'parallel',
        'sbsdetails': 'true',
        'scantimeout': '0.5',
        'sponsorcategories': 'true',
        'statemethod': 'deploybutton',
        'storesubpodexprs': 'true'}
    # # we can also change some parameters
    # params = {
    #     'assumptionsversion': '2',
    #     'banners': 'raw',
    #     'format': 'minput,moutput',
    #     'formattimeout': '8',
    #     'input': input,
    #     'output': 'JSON',
    #     'parsetimeout': '5',
    #     'proxycode': proxy_code,
    #     'scantimeout': '0.5',
    #     'storesubpodexprs': 'true'
    # }
    params = urlencode(params)
    url = "https://www.wolframalpha.com/input/json.jsp?%s" % params
    req = Request(url)
    req.add_header('Referer', "https://www.wolframalpha.com/input/")  # seems important
    resp = opener.open(req)
    # the website returns JSON containing the code
    return json.loads(resp.read().decode("utf-8"))


def parse_moutput_from_json(page_data, verbose=False):
    r"""
    Return the list of outputs found in the json (with key ``'moutput'``).

    INPUT:

    - ``page_data`` -- json obtained from Wolfram Alpha
    - ``verbose`` -- boolean (default: ``False``)

    OUTPUT: list of unicode strings

    EXAMPLES::

        sage: from sage.interfaces.mathematica import request_wolfram_alpha
        sage: from sage.interfaces.mathematica import parse_moutput_from_json
        sage: page_data = request_wolfram_alpha('integrate Sin[x]') # optional internet
        sage: parse_moutput_from_json(page_data)                    # optional internet
        ['-Cos[x]']

    ::

        sage: page_data = request_wolfram_alpha('Sin[x]')           # optional internet
        sage: L = parse_moutput_from_json(page_data)                # optional internet
        sage: sorted(L)                                             # optional internet
        ['-Cos[x]', '{x == 0}', '{x == Pi C[1], Element[C[1], Integers]}']

    TESTS::

        sage: page_data = request_wolfram_alpha('Integrate(Sin[z], y)')  # optional internet
        sage: parse_moutput_from_json(page_data)                         # optional internet
        Traceback (most recent call last):
        ...
        ValueError: asking wolframalpha.com was not successful
    """
    queryresult = page_data['queryresult']
    if not queryresult['success']:
        raise ValueError('asking wolframalpha.com was not successful')
    if 'pods' not in queryresult:
        raise ValueError('json object contains no pods')
    pods = queryresult['pods']
    if verbose:
        print("  Query successful: {}".format(queryresult['success']))
        print("  Number of results: {}".format(len(pods)))
    L = []
    for i, result in enumerate(pods):
        if verbose:
            print("  Result #{}".format(i))
            print("    Title: {}".format(result['title']))
        if 'subpods' not in result:
            continue
        subpods = result['subpods']
        for j, subpod in enumerate(subpods):
            if verbose:
                print("    Subpod #{}".format(j))
                if 'minput' in subpod.keys():
                    print("      MInput: {}".format(subpod['minput']))
                if 'moutput' in subpod.keys():
                    print("      MOutput: {}".format(subpod['moutput']))
            if 'moutput' in subpod.keys():
                L.append(subpod['moutput'])
    return L


def symbolic_expression_from_mathematica_string(mexpr):
    r"""
    Translate a mathematica string into a symbolic expression.

    INPUT:

    - ``mexpr`` -- string

    OUTPUT: symbolic expression

    EXAMPLES::

        sage: from sage.interfaces.mathematica import symbolic_expression_from_mathematica_string
        sage: symbolic_expression_from_mathematica_string('-Cos[x]')
        -cos(x)
    """
    from sage.symbolic.expression import symbol_table
    from sage.symbolic.constants import constants_name_table as constants
    from sage.calculus.calculus import symbolic_expression_from_string
    from sage.calculus.calculus import _find_func as find_func

    expr = mexpr.replace('\n', ' ').replace('\r', '')
    expr = expr.replace('[', '(').replace(']', ')')
    expr = expr.replace('{', '[').replace('}', ']')
    lsymbols = symbol_table['mathematica'].copy()
    lsymbols_names_only = [s[0] for s in lsymbols]
    autotrans = [lambda x:x.lower(),      # Try it in lower case
                 _un_camel,      # Convert `CamelCase` to `camel_case`
                 lambda x: x]     # Try the original name
    # Find the MMA funcs/vars/constants - they start with a letter.
    # Exclude exponents (e.g. 'e8' from 4.e8)
    p = re.compile(r'(?<!\.)[a-zA-Z]\w*')

    for m in p.finditer(expr):
        # If the function, variable or constant is already in the
        # translation dictionary, then just move on.
        if m.group() in lsymbols_names_only:
            pass
        # Now try to translate all other functions -- try each strategy
        # in `autotrans` and check if the function exists in Sage
        elif m.end() < len(expr) and expr[m.end()] == '(':
            for t in autotrans:
                f = find_func(t(m.group()), create_when_missing=False)
                if f is not None:
                    lsymbols[(m.group(), f.number_of_arguments())] = f
                    break
            else:
                raise NotImplementedError("Don't know a Sage equivalent for Mathematica function '%s'." % m.group())
        # Check if Sage has an equivalent constant
        else:
            for t in autotrans:
                if t(m.group()) in constants:
                    lsymbols[(m.group(), 0)] = constants[t(m.group())]
                    break
    return symbolic_expression_from_string(expr, lsymbols, accept_sequence=True)
