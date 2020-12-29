r"""
Interface to KASH

Sage provides an interface to the KASH computer algebra system,
which is a *free* (as in beer!) but *closed source* program for
algebraic number theory that shares much common code with Magma. To
use KASH, you must first install it. Visit its web page:
http://page.math.tu-berlin.de/~kant/kash.html

.. TODO::

    Update the following sentence.

It is not enough to just have KASH installed on your computer.

The KASH interface offers three pieces of functionality:


#. ``kash_console()`` - A function that dumps you into
   an interactive command-line KASH session. Alternatively,

   type ``!kash`` from the Sage prompt.

#. ``kash(expr)`` - Creation of a Sage object that
   wraps a KASH object. This provides a Pythonic interface to KASH.
   For example, if ``f=kash.new(10)``, then
   ``f.Factors()`` returns the prime factorization of
   `10` computed using KASH.

#. ``kash.function_name(args ...)`` - Call the
   indicated KASH function with the given arguments are return the
   result as a KASH object.

#. ``kash.eval(expr)`` - Evaluation of arbitrary KASH
   expressions, with the result returned as a string.


Issues
------

For some reason hitting Control-C to interrupt a calculation
doesn't work correctly. (TODO)

Tutorial
--------

The examples in this tutorial require that kash
be installed.

Basics
~~~~~~

Basic arithmetic is straightforward. First, we obtain the result as
a string.

::

    sage: kash.eval('(9 - 7) * (5 + 6)')                # optional -- kash
    '22'

Next we obtain the result as a new KASH object.

::

    sage: a = kash('(9 - 7) * (5 + 6)'); a              # optional -- kash
    22
    sage: a.parent()                                    # optional -- kash
    Kash

We can do arithmetic and call functions on KASH objects::

    sage: a*a                                           # optional -- kash
    484
    sage: a.Factorial()                                 # optional -- kash
    1124000727777607680000

Integrated Help
~~~~~~~~~~~~~~~

Use the ``kash.help(name)`` command to get help about a
given command. This returns a list of help for each of the
definitions of ``name``. Use ``print
kash.help(name)`` to nicely print out all signatures.

Arithmetic
~~~~~~~~~~

Using the ``kash.new`` command we create Kash objects
on which one can do arithmetic.

::

    sage: a = kash(12345)                          # optional -- kash
    sage: b = kash(25)                             # optional -- kash
    sage: a/b                                      # optional -- kash
    2469/5
    sage: a**b                                     # optional -- kash
    1937659030411463935651167391656422626577614411586152317674869233464019922771432158872187137603759765625

Variable assignment
~~~~~~~~~~~~~~~~~~~

Variable assignment using ``kash`` is takes place in
Sage.

::

    sage: a = kash('32233')                        # optional -- kash
    sage: a                                        # optional -- kash
    32233

In particular, ``a`` is not defined as part of the KASH
session itself.

::

    sage: kash.eval('a')                           # optional -- kash
    "Error, the variable 'a' must have a value"

Use ``a.name()`` to get the name of the KASH variable::

    sage: a.name()                                 # somewhat random; optional - kash
    'sage0'
    sage: kash(a.name())                           # optional -- kash
    32233

Integers and Rationals
~~~~~~~~~~~~~~~~~~~~~~

We illustrate arithmetic with integers and rationals in KASH.

::

    sage: F = kash.Factorization(4352)             # optional -- kash
    sage: F[1]                                     # optional -- kash
    <2, 8>
    sage: F[2]                                     # optional -- kash
    <17, 1>
    sage: F                                        # optional -- kash
    [ <2, 8>, <17, 1> ], extended by:
      ext1 := 1,
      ext2 := Unassign

.. note::

   For some very large numbers KASH's integer factorization seems much
   faster than PARI's (which is the default in Sage).

::

    sage: kash.GCD(15,25)                          # optional -- kash
    5
    sage: kash.LCM(15,25)                          # optional -- kash
    75
    sage: kash.Div(25,15)                          # optional -- kash
    1
    sage: kash(17) % kash(5)                       # optional -- kash
    2
    sage: kash.IsPrime(10007)                      # optional -- kash
    TRUE
    sage: kash.IsPrime(2005)                       # optional -- kash
    FALSE

    sage: kash.NextPrime(10007)                    # optional -- kash
    10009

Real and Complex Numbers
~~~~~~~~~~~~~~~~~~~~~~~~

::

    sage: kash.Precision()                         # optional -- kash
    30
    sage: kash('R')                                # optional -- kash
    Real field of precision 30
    sage: kash.Precision(40)                       # optional -- kash
    40
    sage: kash('R')                                # optional -- kash
    Real field of precision 40
    sage: z = kash('1 + 2*I')                      # optional -- kash
    sage: z                                        # optional -- kash
    1.000000000000000000000000000000000000000 + 2.000000000000000000000000000000000000000*I
    sage: z*z                                      # optional -- kash
    -3.000000000000000000000000000000000000000 + 4.000000000000000000000000000000000000000*I

    sage: kash.Cos('1.24')                         # optional -- kash
    0.3247962844387762365776934156973803996992
    sage: kash('1.24').Cos()                       # optional -- kash
    0.3247962844387762365776934156973803996992

    sage: kash.Exp('1.24')                         # optional -- kash
    3.455613464762675598057615494121998175400

    sage: kash.Precision(30)                       # optional -- kash
    30
    sage: kash.Log('3+4*I')                        # optional -- kash
    1.60943791243410037460075933323 + 0.927295218001612232428512462922*I
    sage: kash.Log('I')                            # optional -- kash
    1.57079632679489661923132169164*I

    sage: kash.Sqrt(4)                             # optional -- kash
    2.00000000000000000000000000000
    sage: kash.Sqrt(2)                             # optional -- kash
    1.41421356237309504880168872421

    sage: kash.Floor('9/5')                        # optional -- kash
    1
    sage: kash.Floor('3/5')                        # optional -- kash
    0

    sage: x_c = kash('3+I')                        # optional -- kash
    sage: x_c.Argument()                           # optional -- kash
    0.321750554396642193401404614359
    sage: x_c.Imaginary()                          # optional -- kash
    1.00000000000000000000000000000

Lists
~~~~~

Note that list appends are completely different in KASH than in
Python. Use underscore after the function name for the mutation
version.

::

    sage: v = kash([1,2,3]); v                    # optional -- kash
    [ 1, 2, 3 ]
    sage: v[1]                                    # optional -- kash
    1
    sage: v[3]                                    # optional -- kash
    3
    sage: v.Append([5])                           # optional -- kash
    [ 1, 2, 3, 5 ]
    sage: v                                       # optional -- kash
    [ 1, 2, 3 ]
    sage: v.Append_([5, 6])                       # optional -- kash
    SUCCESS
    sage: v                                       # optional -- kash
    [ 1, 2, 3, 5, 6 ]
    sage: v.Add(5)                                # optional -- kash
    [ 1, 2, 3, 5, 6, 5 ]
    sage: v                                       # optional -- kash
    [ 1, 2, 3, 5, 6 ]
    sage: v.Add_(5)                               # optional -- kash
    SUCCESS
    sage: v                                       # optional -- kash
    [ 1, 2, 3, 5, 6, 5 ]

The ``Apply`` command applies a function to each
element of a list::

    sage: L = kash([1,2,3,4])                    # optional -- kash
    sage: L.Apply('i -> 3*i')                    # optional -- kash
    [ 3, 6, 9, 12 ]
    sage: L                                      # optional -- kash
    [ 1, 2, 3, 4 ]
    sage: L.Apply('IsEven')                      # optional -- kash
    [ FALSE, TRUE, FALSE, TRUE ]
    sage: L                                      # optional -- kash
    [ 1, 2, 3, 4 ]

Ranges
~~~~~~

the following are examples of ranges.

::

    sage: L = kash('[1..10]')                    # optional -- kash
    sage: L                                      # optional -- kash
    [ 1 .. 10 ]
    sage: L = kash('[2,4..100]')                 # optional -- kash
    sage: L                                      # optional -- kash
    [ 2, 4 .. 100 ]

Sequences
~~~~~~~~~

Tuples
~~~~~~

Polynomials
~~~~~~~~~~~

::

    sage: f = kash('X^3 + X + 1')                # optional -- kash
    sage: f + f                                  # optional -- kash
    2*X^3 + 2*X + 2
    sage: f * f                                  # optional -- kash
    X^6 + 2*X^4 + 2*X^3 + X^2 + 2*X + 1
    sage: f.Evaluate(10)                         # optional -- kash
    1011
    sage: Qx = kash.PolynomialAlgebra('Q')       # optional -- kash
    sage: Qx.gen(1)**5 + kash('7/3')   # sage1 below somewhat random; optional -- kash
    sage1.1^5 + 7/3

Number Fields
~~~~~~~~~~~~~

We create an equation order.

::

    sage: f = kash('X^5 + 4*X^4 - 56*X^2 -16*X + 192')    # optional -- kash
    sage: OK = f.EquationOrder()                          # optional -- kash
    sage: OK                                              # optional -- kash
    Equation Order with defining polynomial X^5 + 4*X^4 - 56*X^2 - 16*X + 192 over Z

::

    sage: f = kash('X^5 + 4*X^4 - 56*X^2 -16*X + 192')    # optional -- kash
    sage: O = f.EquationOrder()                           # optional -- kash
    sage: a = O.gen(2)                                    # optional -- kash
    sage: a                                               # optional -- kash
    [0, 1, 0, 0, 0]
    sage: O.Basis()        # output somewhat random; optional -- kash
    [
    _NG.1,
    _NG.2,
    _NG.3,
    _NG.4,
    _NG.5
    ]
    sage: O.Discriminant()              # optional -- kash
    1364202618880
    sage: O.MaximalOrder()    # name sage2 below somewhat random; optional -- kash
    Maximal Order of sage2

    sage: O = kash.MaximalOrder('X^3 - 77')                  # optional -- kash
    sage: I = O.Ideal(5,[2, 1, 0])                           # optional -- kash
    sage: I                    # name sage14 below random; optional -- kash
    Ideal of sage14
    Two element generators:
    [5, 0, 0]
    [2, 1, 0]

    sage: F = I.Factorisation()                  # optional -- kash
    sage: F                    # name sage14 random; optional -- kash
    [
    <Prime Ideal of sage14
    Two element generators:
    [5, 0, 0]
    [2, 1, 0], 1>
    ]

Determining whether an ideal is principal.

::

    sage: I.IsPrincipal()                      # optional -- kash
    FALSE, extended by:
    ext1 := Unassign

Computation of class groups and unit groups::

    sage: f = kash('X^5 + 4*X^4 - 56*X^2 -16*X + 192')         # optional -- kash
    sage: O = kash.EquationOrder(f)                            # optional -- kash
    sage: OK = O.MaximalOrder()                                # optional -- kash
    sage: OK.ClassGroup()       # name sage32 below random; optional -- kash
    Abelian Group isomorphic to Z/6
      Defined on 1 generator
      Relations:
      6*sage32.1 = 0, extended by:
      ext1 := Mapping from: grp^abl: sage32 to ids/ord^num: _AA

::

    sage: U = OK.UnitGroup()                                  # optional -- kash
    sage: U        # name sage34 below random; optional -- kash
    Abelian Group isomorphic to Z/2 + Z + Z
      Defined on 3 generators
      Relations:
      2*sage34.1 = 0, extended by:
      ext1 := Mapping from: grp^abl: sage34 to ord^num: sage30

    sage: kash.Apply('x->%s.ext1(x)'%U.name(), U.Generators().List())     # optional -- kash
    [ [1, -1, 0, 0, 0], [1, 1, 0, 0, 0], [-1, 0, 0, 0, 0] ]

Function Fields
~~~~~~~~~~~~~~~

::

    sage: k = kash.FiniteField(25)                                 # optional -- kash
    sage: kT = k.RationalFunctionField()                           # optional -- kash
    sage: kTy = kT.PolynomialAlgebra()                             # optional -- kash
    sage: T = kT.gen(1)                                            # optional -- kash
    sage: y = kTy.gen(1)                                           # optional -- kash
    sage: f = y**3 + T**4 + 1                                      # optional -- kash

Long Input
----------

The KASH interface reads in even very long input (using files) in a
robust manner, as long as you are creating a new object.

.. note::

   Using ``kash.eval`` for long input is much less robust, and is not
   recommended.

::

    sage: a = kash(range(10000))                                  # optional -- kash

Note that KASH seems to not support string or integer literals with
more than 1024 digits, which is why the above example uses a list
unlike for the other interfaces.
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

from .expect import Expect, ExpectElement
from sage.docs.instancedoc import instancedoc
import os

from sage.misc.sage_eval import sage_eval


class Kash(Expect):
    r"""
    Interface to the Kash interpreter.

    AUTHORS:

    - William Stein and David Joyner
    """
    def __init__(self,
                 max_workspace_size=None,
                 maxread=None,
                 script_subdirectory=None,
                 restart_on_ctrlc=True,
                 logfile=None,
                 server=None,
                 server_tmpdir=None):
        """
        INPUT:
            max_workspace_size -- (default: None)
                    set maximal workspace memory usage to <mem>
                    <mem> stands for byte-wise allocation
                    <mem>k stands for kilobyte-wise allocation
                    <mem>m stands for megabyte-wise allocation
        """
        cmd = "kash3 -b -c -d  "
        if max_workspace_size is not None:
            cmd += " -a %s" % int(max_workspace_size)
        Expect.__init__(self,
                        name='kash',
                        prompt='kash% ',
                        command=cmd,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=True,
                        verbose_start=False,
                        logfile=logfile,
                        eval_using_file_cutoff=100,
                        init_code=['X:=ZX.1;']
                        )
        # The above init_code programs around a bug reported by Jack Schmidt

        self.__seq = 0

    def clear(self, var):
        """
        Clear the variable named ``var``.

        Kash variables have a record structure, so if sage1 is a
        polynomial ring, sage1.1 will be its indeterminate.  This
        prevents us from easily reusing variables, since sage1.1
        might still have references even if sage1 does not.

        For now, we don't implement variable clearing to avoid these
        problems, and instead implement this method with a noop.
        """
        pass

    def _read_in_file_command(self, filename):
        return 'Read("%s");' % filename

    def _eval_line_using_file(self, line):
        F = open(self._local_tmpfile(), 'w')
        F.write(line)
        F.close()
        tmp_to_use = self._local_tmpfile()
        if self.is_remote():
            self._send_tmpfile_to_server()
            tmp_to_use = self._remote_tmpfile()
        return self._eval_line(self._read_in_file_command(tmp_to_use),
                               allow_use_file=False)

    # Change the default for KASH, since eval using a file doesn't
    # work except for setting variables.
    def _eval_line(self, line, allow_use_file=False, wait_for_prompt=True, restart_if_needed=False):
        return Expect._eval_line(self, line, allow_use_file=allow_use_file,
                                 wait_for_prompt=wait_for_prompt)

    def __reduce__(self):
        return reduce_load_Kash, tuple([])

    def _quit_string(self):
        return 'quit;'

    def _start(self):
        try:
            Expect._start(self)
        except RuntimeError:
            # TODO: replace this error with something more accurate.
            from sage.misc.package import PackageNotFoundError
            raise PackageNotFoundError("kash")
        # Turn off the annoying timer.
        self.eval('Time(false);')

    def _object_class(self):
        return KashElement

    def eval(self, x, newlines=False, strip=True, **kwds):
        r"""
        Send the code in the string s to the Kash interpreter and return
        the output as a string.

        INPUT:


        -  ``s`` - string containing Kash code.

        -  ``newlines`` - bool (default: True); if False,
           remove all backslash-newlines inserted by the Kash output
           formatter.

        -  ``strip`` - ignored
        """
        x = str(x)
        x = x.rstrip()
        if len(x) == 0 or x[len(x) - 1] != ';':
            x += ';'
        s = Expect.eval(self, x, **kwds)
        i = s.find('\r\n')
        if i != -1:
            s = s[i + 2:]
        if newlines:
            return s
        else:
            return s.replace("\\\n", "")

#     def help(self, name=None):
#         """
#         Return help on KASH commands.

#         EXAMPLES::

#             sage: X = kash.help('IntegerRing')   # optional - kash

#         """
#         if name is None:
#           print '\nTo use KASH help enter kash.help(s). '
#           print 'The syntax of the string s is given below.\n'
#           print self.eval('?')
#         elif name[0] == '?':
#           print self.eval(name)
#         else:
#           print self.eval('?%s'%name)

    def help(self, name=None):
        """
        Return help on KASH commands.

        This returns help on all commands with a given name.  If name
        is ``None``, return the location of the installed Kash HTML
        documentation.

        EXAMPLES::

            sage: X = kash.help('IntegerRing')   # random; optional -- kash
            1439: IntegerRing() -> <ord^rat>
            1440: IntegerRing(<elt-ord^rat> m) -> <res^rat>
            1441: IntegerRing(<seq()> Q) -> <res^rat>
            1442: IntegerRing(<fld^rat> K) -> <ord^rat>
            1443: IntegerRing(<fld^fra> K) -> <ord^num>
            1444: IntegerRing(<rng> K) -> <rng>
            1445: IntegerRing(<fld^pad> L) -> <ord^pad>

        There is one entry in X for each item found in the documentation
        for this function: If you type ``print(X[0])`` you will
        get help on about the first one, printed nicely to the screen.

        AUTHORS:

        - Sebastion Pauli (2006-02-04): during Sage coding sprint
        """
        if name is None:
            print('\nTo use KASH help enter kash.help(s). ')
            print('The syntax of the string s is given below.\n')
            print(self.eval('?'))
            return
        name = str(name)
        if name[0] == '?':
            print(self.eval(name))
        else:
            print(self.eval('?%s' % name))

    def _doc(self, V):
        if V.lstrip()[:11] == 'No matches.':
            return KashDocumentation([])
        V = V.split('\n')[1:-1]
        X = []
        for C in V:
            i = C.find('m')
            j = C.find(':')
            try:
                n = int(C[i + 1:j])
            except ValueError:
                full = C
            else:
                full = self.eval('?%s' % n)
            X.append(full)
        return KashDocumentation(X)

    def help_search(self, name):
        return self._doc(self.eval('?*%s' % name))

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s:=%s;;' % (var, value)
        out = self._eval_line(cmd, allow_use_file=True)
        if out.lower().find('error') != -1:
            raise TypeError("Error executing code in Kash\nCODE:\n\t%s\nKash ERROR:\n\t%s" % (cmd, out))

    def get(self, var):
        """
        Get the value of the variable var.
        """
        return self.eval('%s;' % var, newlines=False)

    # def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
    #    self.eval('Unbind(%s)'%var)

    def _contains(self, v1, v2):
        return self.eval('%s in %s' % (v1, v2)) == "TRUE"

    def _assign_symbol(self):
        return ":="

    def _equality_symbol(self):
        return "="

    def _true_symbol(self):
        return "TRUE"

    def _false_symbol(self):
        return "FALSE"

    def function_call(self, function, args=None, kwds=None):
        """
        EXAMPLES::

            sage: kash.function_call('ComplexToPolar', [1+I], {'Results' : 1})   # optional -- kash
            1.41421356237309504880168872421
        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)
        s = self._function_call_string(function,
                                       [s.name() for s in args],
                                       ['%s:=%s' % (key, value.name())
                                        for key, value in kwds.items()])
        return self.new(s)

    def _function_call_string(self, function, args, kwds):
        """
        Return the string used to make function calls.

        EXAMPLES::

            sage: kash._function_call_string('Expand', ['x', 'y'], ['Prec:=10'])
            'Expand(x,y,rec(Prec:=10))'
        """
        if not kwds:
            return "%s(%s)" % (function, ",".join(args))
        return "%s(%s,rec(%s))" % (function, ",".join(args), ",".join(kwds))

    def console(self):
        kash_console()

    def version(self):
        return kash_version()


@instancedoc
class KashElement(ExpectElement):
    def __mod__(self, other):
        self._check_valid()
        if not isinstance(other, KashElement):
            other = self.parent()(other)
        other._check_valid()
        return self.parent()('%s mod %s' % (self._name, other._name))

    def __len__(self):
        self._check_valid()
        return int(self.parent().eval('Length(%s)' % self.name()))

    def __bool__(self):
        """
        Return ``True`` if this Kash element is not 0 or FALSE.

        EXAMPLES::

            sage: bool(kash('FALSE'))                   # optional -- kash
            False
            sage: bool(kash('TRUE'))                    # optional -- kash
            True

            sage: bool(kash(0))                         # optional -- kash
            False
            sage: bool(kash(1))                         # optional -- kash
            True
        """

        # Kash has separate integer and boolean types, and FALSE does not
        # compare equal to 0 (i.e, FALSE = 0 is FALSE)

        # Python 2.x uses __nonzero__ for type conversion to 'bool', so we
        # have to test against FALSE, and sage.structure.element.Element's
        # default implementation of is_zero() is to return 'not self', so
        # our boolean conversion also has to test against 0.

        P = self.parent()
        return (P.eval('%s = FALSE' % self.name()) == 'FALSE' and
                P.eval('%s = 0' % self.name()) == 'FALSE')

    __nonzero__ = __bool__

    def _sage_(self, locals={}, *args):
        """
        Convert this object to Sage.

        A translation dictionary `locals` can be provided to map Kash
        names and objects to Sage objects.

        EXAMPLES::

            sage: kash('1234').sage()                   # optional -- kash
            1234

            sage: kash('X^2+X').sage({'X': x})          # optional -- kash
            x^2 + x

            sage: kQ = kash.RationalField()             # optional -- kash
            sage: kR = kQ.PolynomialAlgebra()           # optional -- kash

            sage: R.<x> = QQ[]                          # optional -- kash
            sage: ka = (x^2+x).subs({x : kR.1})         # random; optional -- kash
            sage541.1^2 + sage541.1
            sage: ka.sage({kR.1: x})                    # optional -- kash
            x^2 + x

            sage: R.<x,y> = QQ[]                        # optional -- kash
            sage: ka = (x^2+x).subs({x : kR.1})         # random; optional -- kash
            sage541.1^2 + sage541.1
            sage: ka.sage({kR.1: x})                    # optional -- kash
            x^2 + x

        """

        string = self._sage_repr()

        i = 1
        parsedict = {}
        for key, val in locals.items():
            name = 'sage' + str(i)
            string = string.replace(str(key), name)
            parsedict[name] = val
            i = i + 1

        try:
            return sage_eval(string, locals=parsedict)
        except Exception:
            raise NotImplementedError("Unable to parse output: %s" % string)


class KashDocumentation(list):
    def __repr__(self):
        if len(self) == 0:
            return "No matches."
        return '\n'.join(self)


def is_KashElement(x):
    return isinstance(x, KashElement)

######

######


kash = Kash()


def reduce_load_Kash():
    return kash


def kash_console():
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%kash magics instead.')
    os.system("kash3 ")


def kash_version():
    return kash.eval('VERSION')


def __doctest_cleanup():
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
