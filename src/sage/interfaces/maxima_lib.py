r"""
Library interface to Maxima

Maxima is a free GPL'd general purpose computer algebra system whose
development started in 1968 at MIT. It contains symbolic manipulation
algorithms, as well as implementations of special functions, including
elliptic functions and generalized hypergeometric functions. Moreover,
Maxima has implementations of many functions relating to the invariant
theory of the symmetric group `S_n`. (However, the commands for group
invariants, and the corresponding Maxima documentation, are in
French.) For many links to Maxima documentation, see
http://maxima.sourceforge.net/documentation.html.

AUTHORS:

- William Stein (2005-12): Initial version

- David Joyner: Improved documentation

- William Stein (2006-01-08): Fixed bug in parsing

- William Stein (2006-02-22): comparisons (following suggestion of
  David Joyner)

- William Stein (2006-02-24): *greatly* improved robustness by adding
  sequence numbers to IO bracketing in _eval_line

- Robert Bradshaw, Nils Bruin, Jean-Pierre Flori (2010,2011): Binary library
  interface

For this interface, Maxima is loaded into ECL which is itself loaded
as a C library in Sage. Translations between Sage and Maxima objects
(which are nothing but wrappers to ECL objects) is made as much as possible
directly, but falls back to the string based conversion used by the
classical Maxima Pexpect interface in case no new implementation has been made.

This interface is the one used for calculus by Sage
and is accessible as `maxima_calculus`::

    sage: maxima_calculus
    Maxima_lib

Only one instance of this interface can be instantiated,
so the user should not try to instantiate another one,
which is anyway set to raise an error::

    sage: from sage.interfaces.maxima_lib import MaximaLib
    sage: MaximaLib()
    Traceback (most recent call last):
    ...
    RuntimeError: Maxima interface in library mode can only be instantiated once

Changed besselexpand to true in init_code -- automatically simplify Bessel functions to trig functions when appropriate when true. Examples:

For some infinite sums, a closed expression can be found. By default, "maxima" is used for that::

    sage: x,n,k = var("x","n","k")
    sage: sum((-x)^n/(factorial(n)*factorial(n+3/2)),n,0,oo)
    -1/2*(2*x*cos(2*sqrt(x)) - sqrt(x)*sin(2*sqrt(x)))/(sqrt(pi)*x^2)

Maxima has some flags that affect how the result gets simplified (By default, besselexpand is false in Maxima; however in 5.39 this test does not show any difference, as, apparently, another expansion path is used)::

    sage: maxima_calculus("besselexpand:false")
    false
    sage: x,n,k = var("x","n","k")
    sage: sum((-x)^n/(factorial(n)*factorial(n+3/2)),n,0,oo)
    -1/2*(2*x*cos(2*sqrt(x)) - sqrt(x)*sin(2*sqrt(x)))/(sqrt(pi)*x^2)
    sage: maxima_calculus("besselexpand:true")
    true

The output is parseable (i. e. :issue:`31796` is fixed)::

    sage: foo = maxima_calculus('a and (b or c)') ; foo
    a and (b or c)
    sage: bar = maxima_calculus(foo) ; bar
    a and (b or c)
    sage: bar == foo
    True

TESTS:

Check our workaround for a race in ecl works, see :issue:`26968`.
We use a temporary `MAXIMA_USERDIR` so it's empty; we place it
in `DOT_SAGE` since we expect it to have more latency than `/tmp`.

    sage: import tempfile, subprocess
    sage: tmpdir = tempfile.TemporaryDirectory(dir=DOT_SAGE)
    sage: _ = subprocess.run(['sage', '-c',  # long time
    ....: f'''
    ....: import os
    ....: os.environ["MAXIMA_USERDIR"] = "{tmpdir.name}"
    ....: if not os.fork():
    ....:     import sage.interfaces.maxima_lib
    ....: else:
    ....:     import sage.interfaces.maxima_lib
    ....:     os.wait()
    ....: '''])
    sage: tmpdir.cleanup()
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

from sage.structure.element import Expression
from sage.symbolic.ring import SR

from sage.libs.ecl import EclObject, ecl_eval

from .maxima_abstract import (MaximaAbstract, MaximaAbstractFunction,
                              MaximaAbstractElement, MaximaAbstractFunctionElement,
                              MaximaAbstractElementFunction)
from sage.misc.instancedoc import instancedoc
from sage.env import MAXIMA_FAS

import sage.rings.real_double
import sage.symbolic.expression
import sage.symbolic.integration.integral

from sage.rings.number_field.number_field_element_base import NumberFieldElement_base
from sage.symbolic.operators import FDerivativeOperator, add_vararg, mul_vararg


# We begin here by initializing Maxima in library mode
# i.e. loading it into ECL
ecl_eval("(setf *load-verbose* NIL)")
if MAXIMA_FAS:
    ecl_eval("(require 'maxima \"{}\")".format(MAXIMA_FAS))
else:
    ecl_eval("(require 'maxima)")
ecl_eval("(in-package :maxima)")
ecl_eval("(set-locale-subdir)")

# This workaround has to happen before any call to (set-pathnames).
# To be safe please do not call anything other than
# (set-locale-subdir) before this block.
try:
    ecl_eval("(set-pathnames)")
except RuntimeError:
    # Recover from :issue:`26968` by creating `*maxima-objdir*` here.
    # This cannot be done before calling `(set-pathnames)` since
    # `*maxima-objdir*` is computed there.
    # We use python `os.makedirs()` which is immune to the race.
    # Using `(ensure-directories-exist ...)` in lisp would be
    # subject to the same race condition and since `*maxima-objdir*`
    # has multiple components this is quite plausible to happen.
    maxima_objdir = ecl_eval("*maxima-objdir*").python()[1:-1]
    import os
    os.makedirs(maxima_objdir, exist_ok=True)
    # Call `(set-pathnames)` again to complete its job.
    ecl_eval("(set-pathnames)")

ecl_eval("(initialize-runtime-globals)")
ecl_eval("(setq $nolabels t))")
ecl_eval("(defun add-lineinfo (x) x)")
ecl_eval('(defun principal nil (cond ($noprincipal (diverg)) ((not pcprntd) (merror "Divergent Integral"))))')
ecl_eval("(remprop 'mfactorial 'grind)")  # don't use ! for factorials (#11539)
ecl_eval("(setf $errormsg nil)")

# The following is an adaptation of the "retrieve" function in maxima
# itself. This routine is normally responsible for displaying a
# question and returning the answer. Our version throws an error in
# which the text of the question is included. This is accomplished by
# redirecting *standard-output* to a string.
#
# After an update in Issue 31553, this routine also preprocesses the
# text to replace space symbols with strings. This prevents those
# symbols from being turned into ugly newlines -- a problem that we
# used to avoid with a custom patch.
ecl_eval(r"""
(defun retrieve (msg flag &aux (print? nil))
  (declare (special msg flag print?))
  (setq msg (mapcar #'(lambda (x) (if (eq x '| |) " " x)) msg))
  (or (eq flag 'noprint) (setq print? t))
  (error
    (concatenate 'string
      "Maxima asks: "
      (string-trim '(#\Newline)
                   (with-output-to-string (*standard-output*)
                     (cond ((not print?)
                            (setq print? t)
                            (format-prompt t ""))
                           ((null msg)
                            (format-prompt t ""))
                           ((atom msg)
                            (format-prompt t "~A" msg)
                            (mterpri))
                           ((eq flag t)
                            (format-prompt t "~{~A~}" (cdr msg))
                            (mterpri))
                           (t
                            (format-prompt t "~M" msg)
                            (mterpri))))))))
""")

# Redirection of ECL and Maxima stdout to /dev/null
ecl_eval(r"""(defparameter *dev-null* (make-two-way-stream
              (make-concatenated-stream) (make-broadcast-stream)))""")
ecl_eval("(setf original-standard-output *standard-output*)")
ecl_eval("(setf *standard-output* *dev-null*)")
# ecl_eval("(setf *error-output* *dev-null*)")

# Default options set in Maxima
# display2d -- no ascii art output
# keepfloat -- don't automatically convert floats to rationals

init_code = ['besselexpand : true', 'display2d : false', 'domain : complex', 'keepfloat : true',
             'load(to_poly_solve)', 'load(simplify_sum)',
             'load(diag)']


# Turn off the prompt labels, since computing them *very
# dramatically* slows down the maxima interpret after a while.
# See the function makelabel in suprv1.lisp.
# Many thanks to andrej.vodopivec@gmail.com and also
# Robert Dodier for figuring this out!
# See trac # 6818.
init_code.append('nolabels : true')
for l in init_code:
    ecl_eval("#$%s$" % l)
# To get more debug information uncomment the next line
# should allow to do this through a method
# ecl_eval("(setf *standard-output* original-standard-output)")

# This is the main function (ECL object) used for evaluation
# This returns an EclObject
maxima_eval = ecl_eval("""
(defun maxima-eval( form )
    (with-$error (meval form)))
""")

# Number of instances of this interface
maxima_lib_instances = 0

# Here we define several useful ECL/Maxima objects
# The Maxima string function can change the structure of its input
# maxprint=EclObject("$STRING")
maxprint = EclObject(r"""(defun mstring-for-sage (form)
                         (coerce (mstring form) 'string))""").eval()
meval = EclObject("MEVAL")
msetq = EclObject("MSETQ")
mlist = EclObject("MLIST")
mequal = EclObject("MEQUAL")
mnotequal = EclObject("MNOTEQUAL")
cadadr = EclObject("CADADR")

max_integrate = EclObject("$INTEGRATE")
max_sum = EclObject("$SUM")
max_simplify_sum = EclObject("$SIMPLIFY_SUM")
max_prod = EclObject("$PRODUCT")
max_simplify_prod = EclObject("$SIMPLIFY_PRODUCT")
max_ratsimp = EclObject("$RATSIMP")
max_limit = EclObject("$LIMIT")
max_tlimit = EclObject("$TLIMIT")
max_plus = EclObject("$PLUS")
max_minus = EclObject("$MINUS")
max_use_grobner = EclObject("$USE_GROBNER")
max_to_poly_solve = EclObject("$TO_POLY_SOLVE")
max_at = EclObject("%AT")


def stdout_to_string(s):
    r"""
    Evaluate command ``s`` and catch Maxima stdout
    (not the result of the command!) into a string.

    INPUT:

    - ``s`` -- string; command to evaluate

    OUTPUT: string

    This is currently used to implement :meth:`~MaximaLibElement.display2d`.

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import stdout_to_string
        sage: stdout_to_string('1+1')
        ''
        sage: stdout_to_string('disp(1+1)')
        '2\n\n'
    """
    return ecl_eval(r"""(with-output-to-string (*standard-output*)
                          (maxima-eval #$%s$))""" % s).python()[1:-1]


def max_to_string(s):
    r"""
    Return the Maxima string corresponding to this ECL object.

    INPUT:

    - ``s`` -- ECL object

    OUTPUT: string

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, max_to_string
        sage: ecl = maxima_lib(cos(x)).ecl()
        sage: max_to_string(ecl)
        'cos(_SAGE_VAR_x)'
    """
    return maxprint(s).python()[1:-1]


my_mread = ecl_eval("""
(defun my-mread (cmd)
  (caddr (mread (make-string-input-stream cmd))))
""")


def parse_max_string(s):
    r"""
    Evaluate string in Maxima without *any* further simplification.

    INPUT:

    - ``s`` -- string

    OUTPUT: ECL object

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import parse_max_string
        sage: parse_max_string('1+1')
        <ECL: ((MPLUS) 1 1)>
    """
    return my_mread('"%s;"' % s)


class MaximaLib(MaximaAbstract):
    """
    Interface to Maxima as a Library.

    OUTPUT: Maxima interface as a Library

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import MaximaLib, maxima_lib
        sage: isinstance(maxima_lib,MaximaLib)
        True

    Only one such interface can be instantiated::

        sage: MaximaLib()
        Traceback (most recent call last):
        ...
        RuntimeError: Maxima interface in library mode can only
        be instantiated once
    """
    def __init__(self):
        """
        Create an instance of the Maxima interpreter.
        See ``MaximaLib`` for full documentation.

        TESTS::

            sage: from sage.interfaces.maxima_lib import MaximaLib, maxima_lib
            sage: MaximaLib == loads(dumps(MaximaLib))
            True
            sage: maxima_lib == loads(dumps(maxima_lib))
            True

        We make sure labels are turned off (see :issue:`6816`)::

            sage: 'nolabels : true' in maxima_lib._MaximaLib__init_code
            True
        """
        global maxima_lib_instances
        if maxima_lib_instances > 0:
            raise RuntimeError("Maxima interface in library mode can only be instantiated once")
        maxima_lib_instances += 1

        global init_code
        self.__init_code = init_code

        MaximaAbstract.__init__(self, "maxima_lib")
        self.__seq = 0

    def _coerce_from_special_method(self, x):
        r"""
        Coerce ``x`` into ``self`` trying to call a special underscore method.

        INPUT:

        - ``x`` -- object to coerce into self

        OUTPUT: Maxima element equivalent to ``x``

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: xmax = maxima_lib._coerce_from_special_method(x)
            sage: type(xmax)
            <class 'sage.interfaces.maxima_lib.MaximaLibElement'>
        """
        if isinstance(x, EclObject):
            return MaximaLibElement(self, self._create(x))
        return MaximaAbstract._coerce_from_special_method(self, x)

    def __reduce__(self):
        r"""
        Implement __reduce__ for ``MaximaLib``.

        OUTPUT: a couple consisting of:

        - the function to call for unpickling

        - a tuple of arguments for the function

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.__reduce__()
            (<function reduce_load_MaximaLib at 0x...>, ())
        """
        return reduce_load_MaximaLib, tuple([])

    # This outputs a string
    def _eval_line(self, line, locals=None, reformat=True, **kwds):
        r"""
        Evaluate the line in Maxima.

        INPUT:

        - ``line`` -- string; text to evaluate

        - ``locals`` -- ``None`` (ignored); this is used for compatibility with the
          Sage notebook's generic system interface

        - ``reformat`` -- boolean; whether to strip output or not

        - ``**kwds`` -- all other arguments are currently ignored

        OUTPUT: string representing Maxima output

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._eval_line('1+1')
            '2'
            sage: maxima_lib._eval_line('1+1;')
            '2'
            sage: maxima_lib._eval_line('1+1$')
            ''
            sage: maxima_lib._eval_line('randvar : cos(x)+sin(y)$')
            ''
            sage: maxima_lib._eval_line('randvar')
            'sin(y)+cos(x)'
        """
        result = ''
        while line:
            ind_dollar = line.find("$")
            ind_semi = line.find(";")
            if ind_dollar == -1 or (ind_semi >= 0 and ind_dollar > ind_semi):
                if ind_semi == -1:
                    statement = line
                    line = ''
                else:
                    statement = line[:ind_semi]
                    line = line[ind_semi + 1:]
                if statement:
                    result = ((result + '\n') if result else '') + max_to_string(maxima_eval("#$%s$" % statement))
            else:
                statement = line[:ind_dollar]
                line = line[ind_dollar + 1:]
                if statement:
                    maxima_eval("#$%s$" % statement)
        if not reformat:
            return result
        return ' '.join(x.strip() for x in result.split())

    eval = _eval_line

    ###########################################
    # Direct access to underlying lisp interpreter.
    ###########################################
    def lisp(self, cmd):
        """
        Send a lisp command to maxima.

        INPUT:

        - ``cmd`` -- string

        OUTPUT: ECL object

        .. NOTE::

           The output of this command is very raw - not pretty.

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.lisp("(+ 2 17)")
            <ECL: 19>
        """
        return ecl_eval(cmd)

    def set(self, var, value):
        """
        Set the variable var to the given value.

        INPUT:

        - ``var`` -- string

        - ``value`` -- string

        OUTPUT: none

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.set('xxxxx', '2')
            sage: maxima_lib.get('xxxxx')
            '2'
        """
        if not isinstance(value, str):
            raise TypeError
        cmd = '%s : %s$' % (var, value.rstrip(';'))
        self.eval(cmd)

    def clear(self, var):
        """
        Clear the variable named var.

        INPUT:

        - ``var`` -- string

        OUTPUT: none

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.set('xxxxx', '2')
            sage: maxima_lib.get('xxxxx')
            '2'
            sage: maxima_lib.clear('xxxxx')
            sage: maxima_lib.get('xxxxx')
            'xxxxx'
        """
        try:
            self.eval('kill(%s)$' % var)
            ecl_eval("(unintern '$%s)" % var)
        except (TypeError, AttributeError):
            pass

    def get(self, var):
        """
        Get the string value of the variable ``var``.

        INPUT:

        - ``var`` -- string

        OUTPUT: string

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.set('xxxxx', '2')
            sage: maxima_lib.get('xxxxx')
            '2'
        """
        s = self.eval('%s;' % var)
        return s

    def _create(self, value, name=None):
        r"""
        Create a variable with given value and name.

        INPUT:

        - ``value`` -- string or ECL object

        - ``name`` -- string (default: ``None``); name to use for the variable,
          an automatically generated name is used if this is none

        OUTPUT: string; the name of the created variable

        EXAMPLES:

        Creation from strings::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._create('3','var3')
            'var3'
            sage: maxima_lib.get('var3')
            '3'
            sage: s = maxima_lib._create('3')
            sage: s # random output
            'sage9'
            sage: s[:4] == 'sage'
            True

        And from ECL object::

            sage: c = maxima_lib(x+cos(19)).ecl()
            sage: maxima_lib._create(c,'m')
            'm'
            sage: maxima_lib.get('m')
            '_SAGE_VAR_x+cos(19)'
            sage: maxima_lib.clear('m')
        """
        name = self._next_var_name() if name is None else name
        try:
            if isinstance(value, EclObject):
                maxima_eval([[msetq], cadadr("#$%s$#$" % name), value])
            else:
                self.set(name, value)
        except RuntimeError as error:
            s = str(error)
            if "Is" in s:  # Maxima asked for a condition
                self._missing_assumption(s)
            else:
                raise
        return name

    def _function_class(self):
        r"""
        Return the Python class of Maxima functions.

        OUTPUT: type

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._function_class()
            <class 'sage.interfaces.interface.InterfaceFunction'>
        """
        return MaximaLibFunction

    def _object_class(self):
        r"""
        Return the Python class of Maxima elements.

        OUTPUT: type

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._object_class()
            <class 'sage.interfaces.maxima_lib.MaximaLibElement'>
        """
        return MaximaLibElement

    def _function_element_class(self):
        r"""
        Return the Python class of Maxima functions of elements.

        OUTPUT: type

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._function_element_class()
            <class 'sage.interfaces.interface.InterfaceFunctionElement'>
        """
        return MaximaLibFunctionElement

    def _object_function_class(self):
        r"""
        Return the Python class of Maxima user-defined functions.

        OUTPUT: type

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._object_function_class()
            <class 'sage.interfaces.maxima_lib.MaximaLibElementFunction'>
        """
        return MaximaLibElementFunction

    # some helper functions to wrap the calculus use of the maxima interface.
    # these routines expect arguments living in the symbolic ring
    # and return something that is hopefully coercible into the symbolic
    # ring again.

    def sr_integral(self, *args):
        """
        Helper function to wrap calculus use of Maxima's integration.

        TESTS::

            sage: a,b=var('a,b')
            sage: integrate(1/(x^3 *(a+b*x)^(1/3)),x)
            Traceback (most recent call last):
            ...
            ValueError: Computation failed since Maxima requested additional
            constraints; using the 'assume' command before evaluation
            *may* help (example of legal syntax is 'assume(a>0)', see
            `assume?` for more details)
            Is a positive or negative?
            sage: assume(a>0)
            sage: integrate(1/(x^3 *(a+b*x)^(1/3)),x)
            2/9*sqrt(3)*b^2*arctan(1/3*sqrt(3)*(2*(b*x + a)^(1/3) + a^(1/3))/a^(1/3))/a^(7/3) - 1/9*b^2*log((b*x + a)^(2/3) + (b*x + a)^(1/3)*a^(1/3) + a^(2/3))/a^(7/3) + 2/9*b^2*log((b*x + a)^(1/3) - a^(1/3))/a^(7/3) + 1/6*(4*(b*x + a)^(5/3)*b^2 - 7*(b*x + a)^(2/3)*a*b^2)/((b*x + a)^2*a^2 - 2*(b*x + a)*a^3 + a^4)
            sage: var('x, n')
            (x, n)
            sage: integral(x^n,x)
            Traceback (most recent call last):
            ...
            ValueError: Computation failed since Maxima requested additional
            constraints; using the 'assume' command before evaluation
            *may* help (example of legal syntax is 'assume(n>0)',
            see `assume?` for more details)
            Is n equal to -1?
            sage: assume(n+1>0)
            sage: integral(x^n,x)
            x^(n + 1)/(n + 1)
            sage: forget()
            sage: assumptions()  # Check the assumptions really were forgotten
            []

        Make sure the abs_integrate package is being used,
        :issue:`11483`. The following are examples from the Maxima
        abs_integrate documentation::

            sage: integrate(abs(x), x)
            1/2*x*abs(x)

        ::

            sage: integrate(sgn(x) - sgn(1-x), x)  # known bug
            abs(x - 1) + abs(x)

        This is a known bug in Sage symbolic limits code, see
        :issue:`17892` and https://sourceforge.net/p/maxima/bugs/3237/ ::

            sage: integrate(1 / (1 + abs(x-5)), x, -5, 6) # not tested -- known bug
            log(11) + log(2)

        ::

            sage: integrate(1/(1 + abs(x)), x)  # known bug
            1/2*(log(x + 1) + log(-x + 1))*sgn(x) + 1/2*log(x + 1) - 1/2*log(-x + 1)

        ::

            sage: integrate(cos(x + abs(x)), x)  # known bug
            -1/2*x*sgn(x) + 1/4*(sgn(x) + 1)*sin(2*x) + 1/2*x

        The last example relies on the following simplification::

            sage: maxima("realpart(signum(x))")
            signum(x)

        An example from sage-support thread e641001f8b8d1129::

            sage: f = e^(-x^2/2)/sqrt(2*pi) * sgn(x-1)
            sage: integrate(f, x, -Infinity, Infinity)  # known bug
            -erf(1/2*sqrt(2))

        From :issue:`8624`::

            sage: integral(abs(cos(x))*sin(x),(x,pi/2,pi))
            1/2

        ::

            sage: integrate(sqrt(x + sqrt(x)), x).canonicalize_radical()  # known bug
            1/12*((8*x - 3)*x^(1/4) + 2*x^(3/4))*sqrt(sqrt(x) + 1) + 1/8*log(sqrt(sqrt(x) + 1) + x^(1/4)) - 1/8*log(sqrt(sqrt(x) + 1) - x^(1/4))

        And :issue:`11594`::

            sage: integrate(abs(x^2 - 1), x, -2, 2)  # known bug
            4

        This definite integral returned zero (incorrectly) in at least
        Maxima 5.23. The correct answer is now given (:issue:`11591`)::

            sage: f = (x^2)*exp(x) / (1+exp(x))^2
            sage: integrate(f, (x, -infinity, infinity))
            1/3*pi^2

        The following integral was computed incorrectly in versions of
        Maxima before 5.27 (see :issue:`12947`)::

            sage: a = integrate(x*cos(x^3),(x,0,1/2)).n()
            sage: a.real()
            0.124756040961038
            sage: a.imag().abs() < 3e-17
            True
        """
        try:
            return max_to_sr(maxima_eval(([max_integrate],
                                          [sr_to_max(SR(a)) for a in args])))
        except RuntimeError as error:
            s = str(error)
            if "Divergent" in s or "divergent" in s:
                # in pexpect interface, one looks for this
                # - e.g. integrate(1/x^3,x,-1,3) gives a principal value
                # if "divergent" in s or 'Principal Value' in s:
                raise ValueError("Integral is divergent.")
            elif "Is" in s:  # Maxima asked for a condition
                self._missing_assumption(s)
            else:
                raise

    def sr_sum(self, *args):
        """
        Helper function to wrap calculus use of Maxima's summation.

        TESTS:

        Check that :issue:`16224` is fixed::

            sage: k = var('k')
            sage: sum(x^(2*k)/factorial(2*k), k, 0, oo).canonicalize_radical()
            cosh(x)

        ::

            sage: x, y, k, n = var('x, y, k, n')
            sage: sum(binomial(n,k) * x^k * y^(n-k), k, 0, n)
            (x + y)^n
            sage: q, a = var('q, a')
            sage: sum(a*q^k, k, 0, oo)
            Traceback (most recent call last):
            ...
            ValueError: Computation failed since Maxima requested additional
            constraints; using the 'assume' command before evaluation *may* help
            (example of legal syntax is 'assume(abs(q)-1>0)', see `assume?`
            for more details)
            Is abs(q)-1 positive, negative or zero?
            sage: assume(q > 1)
            sage: sum(a*q^k, k, 0, oo)
            Traceback (most recent call last):
            ...
            ValueError: Sum is divergent.
            sage: forget()
            sage: assume(abs(q) < 1)
            sage: sum(a*q^k, k, 0, oo)
            -a/(q - 1)
            sage: forget()
            sage: assumptions() # check the assumptions were really forgotten
            []

        Taking the sum of all natural numbers informs us that the sum
        is divergent.  Maxima (before 5.29.1) used to ask questions
        about `m`, leading to a different error (see :issue:`11990`)::

            sage: m = var('m')
            sage: sum(m, m, 0, infinity)
            Traceback (most recent call last):
            ...
            ValueError: Sum is divergent.

        An error with an infinite sum in Maxima (before 5.30.0,
        see :issue:`13712`)::

            sage: n = var('n')
            sage: sum(1/((2*n-1)^2*(2*n+1)^2*(2*n+3)^2), n, 0, oo)
            3/256*pi^2

        Maxima correctly detects division by zero in a symbolic sum
        (see :issue:`11894`)::

            sage: sum(1/(m^4 + 2*m^3 + 3*m^2 + 2*m)^2, m, 0, infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: ECL says: Zero to negative power computed.

        Similar situation for :issue:`12410`::

            sage: x = var('x')
            sage: sum(1/x*(-1)^x, x, 0, oo)
            Traceback (most recent call last):
            ...
            RuntimeError: ECL says: Zero to negative power computed.
        """
        try:
            return max_to_sr(maxima_eval([[max_ratsimp],
                                          [[max_simplify_sum],
                                           ([max_sum],
                                            [sr_to_max(SR(a)) for a in args])]]))
        except RuntimeError as error:
            s = str(error)
            if "divergent" in s:
                # in pexpect interface, one looks for this;
                # could not find an example where 'Pole encountered' occurred, though
                # if "divergent" in s or 'Pole encountered' in s:
                raise ValueError("Sum is divergent.")
            elif "Is" in s:  # Maxima asked for a condition
                self._missing_assumption(s)
            else:
                raise

    def sr_prod(self, *args):
        """
        Helper function to wrap calculus use of Maxima's product.

        TESTS::

            sage: from sage.calculus.calculus import symbolic_product
            sage: _ = var('n')
            sage: symbolic_product(x,x,1,n)
            factorial(n)
            sage: symbolic_product(2*x,x,1,n)
            2^n*factorial(n)
        """
        try:
            return max_to_sr(maxima_eval([[max_ratsimp],
                                          [[max_simplify_prod],
                                           ([max_prod],
                                            [sr_to_max(SR(a)) for a in args])]]))
        except RuntimeError as error:
            s = str(error)
            if "divergent" in s:
                raise ValueError("Product is divergent.")
            elif "Is" in s:  # Maxima asked for a condition
                self._missing_assumption(s)
            else:
                raise

    def sr_limit(self, expr, v, a, dir=None):
        """
        Helper function to wrap calculus use of Maxima's limits.

        TESTS::

            sage: f = (1+1/x)^x
            sage: limit(f,x = oo)
            e
            sage: limit(f,x = 5)
            7776/3125

        Domain to real, a regression in 5.46.0, see https://sf.net/p/maxima/bugs/4138 ::

            sage: maxima_calculus.eval("domain:real")
            ...
            sage: limit(f,x = 1.2).n()
            2.06961575467...
            sage: maxima_calculus.eval("domain:complex");
            ...
            sage: var('a')
            a
            sage: limit(x^a,x=0)
            Traceback (most recent call last):
            ...
            ValueError: Computation failed since Maxima requested additional
            constraints; using the 'assume' command before evaluation
            *may* help (example of legal syntax is 'assume(a>0)', see `assume?`
            for more details)
            Is a positive, negative or zero?
            sage: assume(a>0)
            sage: limit(x^a,x=0)  # random - not needed for maxima 5.46.0
            Traceback (most recent call last):
            ...
            ValueError: Computation failed ...
            Is a an integer?
            sage: assume(a,'integer')
            sage: assume(a,'even')  # Yes, Maxima will ask this too
            sage: limit(x^a,x=0)
            0
            sage: forget()
            sage: assumptions() # check the assumptions were really forgotten
            []

        The second limit below was computed incorrectly prior to
        Maxima 5.24 (:issue:`10868`)::

            sage: f(n) = 2 + 1/factorial(n)
            sage: limit(f(n), n=infinity)
            2
            sage: limit(1/f(n), n=infinity)
            1/2

        The limit below was computed incorrectly prior to Maxima 5.30
        (see :issue:`13526`)::

            sage: n = var('n')
            sage: l = (3^n + (-2)^n) / (3^(n+1) + (-2)^(n+1))
            sage: l.limit(n=oo)
            1/3

        The following limit computation used to incorrectly return 0
        or infinity, depending on the domain (see :issue:`15033`)::

            sage: m = sage.calculus.calculus.maxima
            sage: _ = m.eval('domain: real')   # much faster than 'domain: complex'
            sage: limit(gamma(x + 1/2)/(sqrt(x)*gamma(x)), x=infinity)
            1
            sage: _ = m.eval('domain: complex')
        """
        try:
            L = [sr_to_max(SR(aa)) for aa in [expr, v, a]]
            if dir == "plus":
                L.append(max_plus)
            elif dir == "minus":
                L.append(max_minus)
            return max_to_sr(maxima_eval(([max_limit], L)))
        except RuntimeError as error:
            s = str(error)
            if "Is" in s:  # Maxima asked for a condition
                self._missing_assumption(s)
            else:
                raise

    def sr_tlimit(self, expr, v, a, dir=None):
        """
        Helper function to wrap calculus use of Maxima's Taylor series limits.

        TESTS::

            sage: f = (1+1/x)^x
            sage: limit(f, x = I, taylor=True)
            (-I + 1)^I
        """
        L = [sr_to_max(SR(aa)) for aa in [expr, v, a]]
        if dir == "plus":
            L.append(max_plus)
        elif dir == "minus":
            L.append(max_minus)
        return max_to_sr(maxima_eval(([max_tlimit], L)))

    def _missing_assumption(self, errstr):
        """
        Helper function for unified handling of failed computation because an
        assumption was missing.

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._missing_assumption('Is xyz a thing?')
            Traceback (most recent call last):
            ...
            ValueError: Computation failed ...
            Is xyz a thing?
        """
        j = errstr.find('Is ')
        errstr = errstr[j:]
        jj = 2
        if errstr[3] == ' ':
            jj = 3
        k = errstr.find(' ', jj + 1)

        outstr = "Computation failed since Maxima requested additional constraints; using the 'assume' command before evaluation *may* help (example of legal syntax is 'assume("\
            + errstr[jj + 1:k] + ">0)', see `assume?` for more details)\n" + errstr
        outstr = outstr.replace('_SAGE_VAR_', '')
        raise ValueError(outstr)


def is_MaximaLibElement(x):
    r"""
    Return ``True`` if ``x`` is of type :class:`MaximaLibElement`.

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, is_MaximaLibElement
        sage: is_MaximaLibElement(1)
        doctest:...: DeprecationWarning: the function is_MaximaLibElement is deprecated; use isinstance(x, sage.interfaces.abc.MaximaLibElement) instead
        See https://github.com/sagemath/sage/issues/34804 for details.
        False
        sage: m = maxima_lib(1)
        sage: is_MaximaLibElement(m)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(34804, "the function is_MaximaLibElement is deprecated; use isinstance(x, sage.interfaces.abc.MaximaLibElement) instead")

    return isinstance(x, MaximaLibElement)


@instancedoc
class MaximaLibElement(MaximaAbstractElement):
    r"""
    Element of Maxima through library interface.

    EXAMPLES:

    Elements of this class should not be created directly.
    The targeted parent should be used instead::

        sage: from sage.interfaces.maxima_lib import maxima_lib
        sage: maxima_lib(4)
        4
        sage: maxima_lib(log(x))
        log(_SAGE_VAR_x)
    """

    def ecl(self):
        r"""
        Return the underlying ECL object of this MaximaLib object.

        OUTPUT: ECL object

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib(x+cos(19)).ecl()
            <ECL: ((MPLUS SIMP) ((%COS SIMP) 19) |$_SAGE_VAR_x|)>
        """
        try:
            return self._ecl
        except AttributeError:
            self._ecl = maxima_eval("#$%s$" % self._name)
            return self._ecl

    def to_poly_solve(self, vars, options=""):
        r"""
        Use Maxima's to_poly_solver package.

        INPUT:

        - ``vars`` -- symbolic expressions

        - ``options`` -- string (default="")

        OUTPUT: Maxima object

        EXAMPLES:

        The zXXX below are names for arbitrary integers and
        subject to change::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: sol = maxima_lib(sin(x) == 0).to_poly_solve(x)
            sage: sol.sage()
            [[x == pi*z...]]
        """
        if options.find("use_grobner=true") != -1:
            cmd = EclObject([[max_to_poly_solve], self.ecl(), sr_to_max(vars),
                             [[mequal], max_use_grobner, True]])
        else:
            cmd = EclObject([[max_to_poly_solve], self.ecl(), sr_to_max(vars)])
        return self.parent()(maxima_eval(cmd))

    def display2d(self, onscreen=True):
        r"""
        Return the 2d representation of this Maxima object.

        INPUT:

        - ``onscreen`` -- boolean (default: ``True``); whether to print or return

        OUTPUT:

        The representation is printed if onscreen is set to True
        and returned as a string otherwise.

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: F = maxima_lib('x^5 - y^5').factor()
            sage: F.display2d()
                                   4      3    2  2    3      4
                       - (y - x) (y  + x y  + x  y  + x  y + x )
        """
        self._check_valid()
        P = self.parent()
        P._eval_line('display2d : true$')
        s = stdout_to_string('disp(%s)' % self.name())
        # s = P._eval_line('disp(%s)$'%self.name())
        P._eval_line('display2d : false$')
        s = s.strip('\r\n')

        # if ever want to dedent, see
        # http://mail.python.org/pipermail/python-list/2006-December/420033.html
        if onscreen:
            print(s)
        else:
            return s


MaximaLibFunctionElement = MaximaAbstractFunctionElement
MaximaLibFunction = MaximaAbstractFunction


@instancedoc
class MaximaLibElementFunction(MaximaLibElement, MaximaAbstractElementFunction):
    pass


# The (unique) instance
maxima_lib = MaximaLib()
maxima = maxima_lib


def reduce_load_MaximaLib():
    r"""
    Unpickle the (unique) Maxima library interface.

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import reduce_load_MaximaLib
        sage: reduce_load_MaximaLib()
        Maxima_lib
    """
    return maxima_lib


#############################################
# Smart translations between SR and Maxima
#############################################

car = EclObject("car")
cdr = EclObject("cdr")
caar = EclObject("caar")
cadr = EclObject("cadr")
cddr = EclObject("cddr")
caddr = EclObject("caddr")
caaadr = EclObject("caaadr")
cadadr = EclObject("cadadr")
meval = EclObject("meval")
NIL = EclObject("NIL")
lisp_length = EclObject("length")

# Dictionaries for standard operators
sage_op_dict = {
    sage.functions.other.abs: "MABS",
    add_vararg: "MPLUS",
    sage.symbolic.expression.operator.truediv: "MQUOTIENT",
    sage.symbolic.expression.operator.eq: "MEQUAL",
    sage.symbolic.expression.operator.ge: "MGEQP",
    sage.symbolic.expression.operator.gt: "MGREATERP",
    sage.symbolic.expression.operator.le: "MLEQP",
    sage.symbolic.expression.operator.lt: "MLESSP",
    mul_vararg: "MTIMES",
    sage.symbolic.expression.operator.ne: "MNOTEQUAL",
    sage.symbolic.expression.operator.neg: "MMINUS",
    sage.symbolic.expression.operator.pow: "MEXPT",
    sage.symbolic.expression.operator.or_: "MOR",
    sage.symbolic.expression.operator.and_: "MAND",
    sage.functions.log.ln: "%LOG",
    sage.functions.log.log: "%LOG",
    sage.functions.log.lambert_w: "%LAMBERT_W",
    sage.functions.other.factorial: "MFACTORIAL",
    sage.functions.error.erf: "%ERF",
    sage.functions.gamma.gamma_inc: "%GAMMA_INCOMPLETE",
    sage.functions.other.conjugate: "$CONJUGATE",
}
# we compile the dictionary
sage_op_dict = dict([(k, EclObject(sage_op_dict[k])) for k in sage_op_dict])
max_op_dict = dict([(sage_op_dict[k], k) for k in sage_op_dict])


# Here we correct the dictionaries for some simple operators

max_op_dict[mequal] = sage.symbolic.expression.operator.eq
max_op_dict[mnotequal] = sage.symbolic.expression.operator.ne


def sage_rat(x, y):
    r"""
    Return quotient x/y.

    INPUT:

    - ``x`` -- integer

    - ``y`` -- integer

    OUTPUT: rational

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import sage_rat
        sage: sage_rat(1,7)
        1/7
    """
    return x / y


mplus = EclObject("MPLUS")
mtimes = EclObject("MTIMES")
rat = EclObject("RAT")
max_op_dict[mplus] = add_vararg
max_op_dict[mtimes] = mul_vararg
max_op_dict[rat] = sage_rat


# Here we build dictionaries for operators needing special conversions.
ratdisrep = EclObject("ratdisrep")
mrat = EclObject("MRAT")
mqapply = EclObject("MQAPPLY")
max_li = EclObject("$LI")
max_psi = EclObject("$PSI")
max_hyper = EclObject("$%F")
max_array = EclObject("ARRAY")
mdiff = EclObject("%DERIVATIVE")
max_lambert_w = sage_op_dict[sage.functions.log.lambert_w]
max_harmo = EclObject("$GEN_HARMONIC_NUMBER")
max_pochhammer = EclObject("$POCHHAMMER")


def mrat_to_sage(expr):
    r"""
    Convert a Maxima MRAT expression to Sage SR.

    INPUT:

    - ``expr`` -- ECL object; a Maxima MRAT expression

    OUTPUT: symbolic expression

    Maxima has an optimised representation for multivariate
    rational expressions. The easiest way to translate those
    to SR is by first asking Maxima to give the generic representation
    of the object. That is what RATDISREP does in Maxima.

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, mrat_to_sage
        sage: var('x y z')
        (x, y, z)
        sage: c = maxima_lib((x+y^2+z^9)/x^6+z^8/y).rat()
        sage: c
        (_SAGE_VAR_y*_SAGE_VAR_z^9+_SAGE_VAR_x^6*_SAGE_VAR_z^8+_SAGE_VAR_y^3+_SAGE_VAR_x*_SAGE_VAR_y)/(_SAGE_VAR_x^6*_SAGE_VAR_y)
        sage: c.ecl()
        <ECL: ((MRAT SIMP (|$_SAGE_VAR_x| |$_SAGE_VAR_y| |$_SAGE_VAR_z|)
        ...>
        sage: mrat_to_sage(c.ecl())
        (x^6*z^8 + y*z^9 + y^3 + x*y)/(x^6*y)
    """
    return max_to_sr(meval(EclObject([[ratdisrep], expr])))


def mqapply_to_sage(expr):
    r"""
    Special conversion rule for MQAPPLY expressions.

    INPUT:

    - ``expr`` -- ECL object; a Maxima MQAPPLY expression

    OUTPUT: symbolic expression

    MQAPPLY is used for function as li[x](y) and psi[x](y).

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, mqapply_to_sage
        sage: c = maxima_lib('li[2](3)')
        sage: c.ecl()
        <ECL: ((MQAPPLY SIMP) (($LI SIMP ARRAY) 2) 3)>
        sage: mqapply_to_sage(c.ecl())
        dilog(3)
    """
    if caaadr(expr) == max_li:
        return sage.functions.log.polylog(max_to_sr(cadadr(expr)),
                                          max_to_sr(caddr(expr)))
    if caaadr(expr) == max_psi:
        return sage.functions.gamma.psi(max_to_sr(cadadr(expr)),
                                        max_to_sr(caddr(expr)))
    if caaadr(expr) == max_hyper:
        return sage.functions.hypergeometric.hypergeometric(mlist_to_sage(car(cdr(cdr(expr)))),
                                                            mlist_to_sage(car(cdr(cdr(cdr(expr))))),
                                                            max_to_sr(car(cdr(cdr(cdr(cdr(expr)))))))
    else:
        op = max_to_sr(cadr(expr))
        max_args = cddr(expr)
        args = [max_to_sr(a) for a in max_args]
        return op(*args)


def mdiff_to_sage(expr):
    r"""
    Special conversion rule for %DERIVATIVE expressions.

    INPUT:

    - ``expr`` -- ECL object; a Maxima %DERIVATIVE expression

    OUTPUT: symbolic expression

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, mdiff_to_sage
        sage: f = maxima_lib('f(x)').diff('x',4)
        sage: f.ecl()
        <ECL: ((%DERIVATIVE SIMP) (($F SIMP) $X) $X 4)>
        sage: mdiff_to_sage(f.ecl())
        diff(f(x), x, x, x, x)
    """
    return max_to_sr(expr.cadr()).diff(*[max_to_sr(e) for e in expr.cddr()])


def mlist_to_sage(expr):
    r"""
    Special conversion rule for MLIST expressions.

    INPUT:

    - ``expr`` -- ECL object; a Maxima MLIST expression (i.e., a list)

    OUTPUT: a Python list of converted expressions

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, mlist_to_sage
        sage: L=maxima_lib("[1,2,3]")
        sage: L.ecl()
        <ECL: ((MLIST SIMP) 1 2 3)>
        sage: mlist_to_sage(L.ecl())
        [1, 2, 3]
    """
    return [max_to_sr(x) for x in expr.cdr()]


def max_at_to_sage(expr):
    r"""
    Special conversion rule for AT expressions.

    INPUT:

    - ``expr`` -- ECL object; a Maxima AT expression

    OUTPUT: symbolic expression

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, max_at_to_sage
        sage: a=maxima_lib("'at(f(x,y,z),[x=1,y=2,z=3])")
        sage: a
        'at(f(x,y,z),[x = 1,y = 2,z = 3])
        sage: max_at_to_sage(a.ecl())
        f(1, 2, 3)
        sage: a=maxima_lib("'at(f(x,y,z),x=1)")
        sage: a
        'at(f(x,y,z),x = 1)
        sage: max_at_to_sage(a.ecl())
        f(1, y, z)
    """
    arg = max_to_sr(expr.cadr())
    subsarg = caddr(expr)
    if caar(subsarg) == mlist:
        subsvalues = {v.lhs(): v.rhs() for v in max_to_sr(subsarg)}
    else:
        v = max_to_sr(subsarg)
        subsvalues = dict([(v.lhs(), v.rhs())])
    return SR(arg).subs(subsvalues)


def dummy_integrate(expr):
    r"""
    We would like to simply tie Maxima's integrate to
    sage.calculus.calculus.dummy_integrate, but we're being
    imported there so to avoid circularity we define it here.

    INPUT:

    - ``expr`` -- ECL object; a Maxima %INTEGRATE expression

    OUTPUT: symbolic expression

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, dummy_integrate
        sage: f = maxima_lib('f(x)').integrate('x')
        sage: f.ecl()
        <ECL: ((%INTEGRATE SIMP) (($F SIMP) $X) $X)>
        sage: dummy_integrate(f.ecl())
        integrate(f(x), x)

    ::

        sage: f = maxima_lib('f(x)').integrate('x',0,10)
        sage: f.ecl()
        <ECL: ((%INTEGRATE SIMP) (($F SIMP) $X) $X 0 10)>
        sage: dummy_integrate(f.ecl())
        integrate(f(x), x, 0, 10)
    """
    args = [max_to_sr(a) for a in cdr(expr)]
    if len(args) == 4:
        return sage.symbolic.integration.integral.definite_integral(*args,
                                                                    hold=True)
    return sage.symbolic.integration.integral.indefinite_integral(*args,
                                                                  hold=True)


def max_harmonic_to_sage(expr):
    """
    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, max_to_sr
        sage: c=maxima_lib(harmonic_number(x,2))
        sage: c.ecl()
        <ECL: (($GEN_HARMONIC_NUMBER SIMP) 2 |$_SAGE_VAR_x|)>
        sage: max_to_sr(c.ecl())
        harmonic_number(x, 2)
    """
    return sage.functions.log.harmonic_number(max_to_sr(caddr(expr)),
                                              max_to_sr(cadr(expr)))


def max_pochhammer_to_sage(expr):
    """
    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, max_to_sr
        sage: c = maxima_lib('pochhammer(x,n)')
        sage: c.ecl()
        <ECL: (($POCHHAMMER SIMP) $X $N)>
        sage: max_to_sr(c.ecl())
        gamma(n + x)/gamma(x)
    """
    from sage.functions.gamma import gamma
    x = max_to_sr(cadr(expr))
    y = max_to_sr(caddr(expr))
    return gamma(x + y) / gamma(x)


# The dictionaries
special_max_to_sage = {
    mrat: mrat_to_sage,
    mqapply: mqapply_to_sage,
    mdiff: mdiff_to_sage,
    EclObject("%INTEGRATE"): dummy_integrate,
    max_at: max_at_to_sage,
    mlist: mlist_to_sage,
    max_harmo: max_harmonic_to_sage,
    max_pochhammer: max_pochhammer_to_sage
}

special_sage_to_max = {
    sage.functions.log.polylog: lambda N, X: [[mqapply], [[max_li, max_array], N], X],
    sage.functions.gamma.psi1: lambda X: [[mqapply], [[max_psi, max_array], 0], X],
    sage.functions.gamma.psi2: lambda N, X: [[mqapply], [[max_psi, max_array], N], X],
    sage.functions.log.lambert_w: lambda N, X: [[max_lambert_w], X] if N == EclObject(0) else [[mqapply], [[max_lambert_w, max_array], N], X],
    sage.functions.log.harmonic_number: lambda N, X: [[max_harmo], X, N],
    sage.functions.hypergeometric.hypergeometric: lambda A, B, X: [[mqapply], [[max_hyper, max_array], lisp_length(A.cdr()), lisp_length(B.cdr())], A, B, X]
}


# Dictionaries for symbols
sage_sym_dict = {}
max_sym_dict = {}


# Generic conversion functions

max_i = EclObject("$%I")


def pyobject_to_max(obj):
    r"""
    Convert a (simple) Python object into a Maxima object.

    INPUT:

    - ``expr`` -- Python object

    OUTPUT: ECL object

    .. NOTE::

       This uses functions defined in sage.libs.ecl.

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import pyobject_to_max
        sage: pyobject_to_max(4)
        <ECL: 4>
        sage: pyobject_to_max('z')
        <ECL: Z>
        sage: var('x')
        x
        sage: pyobject_to_max(x)
        Traceback (most recent call last):
        ...
        TypeError: Unimplemented type for python_to_ecl
    """
    if isinstance(obj, sage.rings.rational.Rational):
        return EclObject(obj) if (obj.denom().is_one()) else EclObject([[rat], obj.numer(), obj.denom()])
    elif isinstance(obj, NumberFieldElement_base):
        from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic
        if isinstance(obj, NumberFieldElement_quadratic) and obj.parent().defining_polynomial().list() == [1, 0, 1]:
            re, im = obj.list()
            return EclObject([[mplus], pyobject_to_max(re), [[mtimes], pyobject_to_max(im), max_i]])
    return EclObject(obj)


# This goes from SR to EclObject
def sr_to_max(expr):
    r"""
    Convert a symbolic expression into a Maxima object.

    INPUT:

    - ``expr`` -- symbolic expression

    OUTPUT: ECL object

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import sr_to_max
        sage: var('x')
        x
        sage: sr_to_max(x)
        <ECL: $X>
        sage: sr_to_max(cos(x))
        <ECL: ((%COS) $X)>
        sage: f = function('f')(x)
        sage: sr_to_max(f.diff())
        <ECL: ((%DERIVATIVE) (($F) $X) $X 1)>

    TESTS:

    We should be able to convert derivatives evaluated at a point,
    :issue:`12796`::

        sage: from sage.interfaces.maxima_lib import sr_to_max, max_to_sr
        sage: f = function('f')
        sage: f_prime = f(x).diff(x)
        sage: max_to_sr(sr_to_max(f_prime(x = 1)))
        D[0](f)(1)
    """
    global sage_op_dict, max_op_dict
    global sage_sym_dict, max_sym_dict
    if isinstance(expr, (list, tuple)):
        return EclObject(([mlist], [sr_to_max(e) for e in expr]))
    op = expr.operator()
    if op:
        # Stolen from sage.symbolic.expression_conversion
        # Should be defined in a function and then put in special_sage_to_max
        # For that, we should change the API of the functions there
        # (we need to have access to op, not only to expr.operands()
        if isinstance(op, FDerivativeOperator):
            args = expr.operands()
            if (not all(isinstance(v, Expression) and v.is_symbol() for v in args)
                    or len(args) != len(set(args))):
                # An evaluated derivative of the form f'(1) is not a
                # symbolic variable, yet we would like to treat it
                # like one. So, we replace the argument `1` with a
                # temporary variable e.g. `_symbol0` and then evaluate
                # the derivative f'(_symbol0) symbolically at
                # _symbol0=1. See trac #12796. Note that we cannot use
                # SR.temp_var here since two conversions of the same
                # expression have to be equal.
                temp_args = [SR.symbol("_symbol%s" % i) for i in range(len(args))]
                f = sr_to_max(op.function()(*temp_args))
                params = op.parameter_set()
                deriv_max = [[mdiff], f]
                for i in set(params):
                    deriv_max.extend([sr_to_max(temp_args[i]), EclObject(params.count(i))])
                at_eval = sr_to_max([temp_args[i] == args[i] for i in range(len(args))])
                return EclObject([[max_at], deriv_max, at_eval])

            f = sr_to_max(op.function()(*args))
            params = op.parameter_set()
            deriv_max = []
            [deriv_max.extend([sr_to_max(args[i]), EclObject(params.count(i))]) for i in set(params)]
            l = [[mdiff], f]
            l.extend(deriv_max)
            return EclObject(l)
        elif (op in special_sage_to_max):
            return EclObject(special_sage_to_max[op](*[sr_to_max(o) for o in expr.operands()]))
        elif op == tuple:
            return EclObject(([mlist], list(sr_to_max(op) for op in expr.operands())))
        elif not (op in sage_op_dict):
            # Maxima does some simplifications automatically by default
            # so calling maxima(expr) can change the structure of expr
            # op_max=caar(maxima(expr).ecl())
            # This should be safe if we treated all special operators above
            # furthermore, this should already use any _maxima_ methods on op, so use any
            # conversion methods that are registered in pynac.
            op_max = maxima(op).ecl()
            if op_max in max_op_dict:
                raise RuntimeError("Encountered operator mismatch in sr-to-maxima translation")
            sage_op_dict[op] = op_max
            max_op_dict[op_max] = op
        return EclObject(([sage_op_dict[op]],
                          [sr_to_max(o) for o in expr.operands()]))
    elif expr.is_symbol() or expr._is_registered_constant_():
        if expr not in sage_sym_dict:
            sym_max = maxima(expr).ecl()
            sage_sym_dict[expr] = sym_max
            max_sym_dict[sym_max] = expr
        return sage_sym_dict[expr]
    else:
        try:
            return pyobject_to_max(expr.pyobject())
        except TypeError:
            return maxima(expr).ecl()


# This goes from EclObject to SR
from sage.symbolic.expression import symbol_table
max_to_pynac_table = symbol_table['maxima']


def max_to_sr(expr):
    r"""
    Convert a Maxima object into a symbolic expression.

    INPUT:

    - ``expr`` -- ECL object

    OUTPUT: symbolic expression

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, max_to_sr
        sage: f = maxima_lib('f(x)')
        sage: f.ecl()
        <ECL: (($F SIMP) $X)>
        sage: max_to_sr(f.ecl())
        f(x)

    TESTS::

        sage: from sage.interfaces.maxima_lib import sr_to_max, max_to_sr
        sage: f = function('f')(x).diff()
        sage: bool(max_to_sr(sr_to_max(f)) == f)
        True
    """
    if expr.consp():
        op_max = caar(expr)
        if op_max in special_max_to_sage:
            return special_max_to_sage[op_max](expr)
        if op_max not in max_op_dict:
            op_max_str = maxprint(op_max).python()[1:-1]
            if op_max_str in max_to_pynac_table:   # nargs ?
                op = max_to_pynac_table[op_max_str]
            else:
                # This could be unsafe if the conversion to SR
                # changes the structure of expr
                sage_expr = SR(maxima(expr))
                op = sage_expr.operator()
            if op in sage_op_dict:
                raise RuntimeError("Encountered operator mismatch in maxima-to-sr translation")
            max_op_dict[op_max] = op
            sage_op_dict[op] = op_max
        else:
            op = max_op_dict[op_max]
        max_args = cdr(expr)
        args = [max_to_sr(a) for a in max_args]
        return op(*args)
    elif expr.symbolp():
        if expr not in max_sym_dict:
            sage_symbol = SR(maxima(expr))
            sage_sym_dict[sage_symbol] = expr
            max_sym_dict[expr] = sage_symbol
        return max_sym_dict[expr]
    else:
        e = expr.python()
        if isinstance(e, float):
            return sage.rings.real_double.RealDoubleElement(e)
        return e
