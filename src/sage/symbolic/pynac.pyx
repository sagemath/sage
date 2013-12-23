###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

# distutils: libraries = gsl

cdef extern from "pynac_cc.h":
    long double sage_logl(long double)
    long double sage_sqrtl(long double)
    long double sage_tgammal(long double)
    long double sage_lgammal(long double)

include "sage/ext/cdefs.pxi"
include "sage/ext/stdsage.pxi"
include "sage/ext/python.pxi"

from sage.libs.ginac cimport *

# for complex log and log gamma
include "sage/gsl/gsl_complex.pxi"
include "sage/gsl/gsl_sf_result.pxi"
include "sage/gsl/gsl_gamma.pxi"


from sage.structure.element import Element
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.real_mpfr import RR, RealField
from sage.rings.complex_field import ComplexField
from sage.rings.all import CC

from sage.symbolic.expression cimport Expression, new_Expression_from_GEx
from sage.symbolic.function import get_sfunction_from_serial
from sage.symbolic.function cimport Function, parent_c
from sage.symbolic.constants_c cimport PynacConstant

import ring

from sage.rings.integer cimport Integer

import math
from sage.libs.mpmath import utils as mpmath_utils

#################################################################
# Symbolic function helpers
#################################################################


cdef public object ex_to_pyExpression(GEx juice):
    """
    Convert given GiNaC::ex object to a python Expression instance.

    Used to pass parameters to custom power and series functions.
    """
    cdef Expression nex
    nex = <Expression>PY_NEW(Expression)
    GEx_construct_ex(&nex._gobj, juice)
    nex._parent = ring.SR
    return nex

cdef public object exvector_to_PyTuple(GExVector seq):
    """
    Converts arguments list given to a function to a PyTuple.

    Used to pass arguments to python methods assigned to custom
    evaluation, derivative, etc. functions of symbolic functions.

    We convert Python objects wrapped in symbolic expressions back to regular
    Python objects.
    """
    from sage.symbolic.ring import SR
    res = []
    for i in range(seq.size()):
        if is_a_numeric(seq.at(i)):
            res.append(py_object_from_numeric(seq.at(i)))
        else:
            res.append(new_Expression_from_GEx(SR, seq.at(i)))
    return tuple(res)

cdef public GEx pyExpression_to_ex(object res) except *:
    """
    Converts an Expression object to a GiNaC::ex.

    Used to pass return values of custom python evaluation, derivation
    functions back to C++ level.
    """
    if res is None:
        raise TypeError, "function returned None, expected return value of type sage.symbolic.expression.Expression"
    try:
        t = ring.SR.coerce(res)
    except TypeError, err:
        raise TypeError, "function did not return a symbolic expression or an element that can be coerced into a symbolic expression"
    return (<Expression>t)._gobj

cdef public object paramset_to_PyTuple(const_paramset_ref s):
    """
    Converts a std::multiset<unsigned> to a PyTuple.

    Used to pass a list of parameter numbers with respect to which a function
    is differentiated to the printing functions py_print_fderivative and
    py_latex_fderivative.
    """
    cdef GParamSetIter itr = s.begin()
    res = []
    while itr.is_not_equal(s.end()):
        res.append(itr.obj())
        itr.inc()
    return res

def paramset_from_Expression(Expression e):
    """
    EXAMPLES::

        sage: from sage.symbolic.pynac import paramset_from_Expression
        sage: f = function('f')
        sage: paramset_from_Expression(f(x).diff(x))
        [0L]
    """
    return paramset_to_PyTuple(ex_to_fderivative(e._gobj).get_parameter_set())

cdef int GINAC_FN_SERIAL = 0

cdef set_ginac_fn_serial():
    """
    Initialize the GINAC_FN_SERIAL variable to the number of functions
    defined by GiNaC. This allows us to prevent collisions with C++ level
    special functions when a user asks to construct a symbolic function
    with the same name.
    """
    global GINAC_FN_SERIAL
    GINAC_FN_SERIAL = g_registered_functions().size()

cdef public int py_get_ginac_serial():
    """
    Returns the number of C++ level functions defined by GiNaC.

    EXAMPLES::

        sage: from sage.symbolic.pynac import get_ginac_serial
        sage: get_ginac_serial() >= 40
        True
    """
    global GINAC_FN_SERIAL
    return GINAC_FN_SERIAL

def get_ginac_serial():
    """
    Number of C++ level functions defined by GiNaC. (Defined mainly for testing.)

    EXAMPLES::

        sage: sage.symbolic.pynac.get_ginac_serial() >= 40
        True
    """
    return py_get_ginac_serial()

#################################################################
# Printing helpers
#################################################################

##########################################################################
# Pynac's precedence levels, as extracted from the raw source code on
# 2009-05-15.  If this changes in Pynac it could cause a bug in
# printing.  But it's hardcoded in Pynac now, so there's not much to
# be done (at present).
#    Container: 10
#    Expairseq: 10
#    Relational: 20
#    Numeric: 30
#    Pseries: 38
#    Addition: 40
#    Integral: 45
#    Multiplication: 50
#    Noncummative mult: 50
#    Index: 55
#    Power: 60
#    Clifford: 65
#    Function: 70
#    Structure: 70
##########################################################################

cdef public stdstring* py_repr(object o, int level) except +:
    """
    Return string representation of o.  If level > 0, possibly put
    parentheses around the string.
    """
    s = repr(o)
    if level >= 20:
        # s may need parens (e.g., is in an exponent), so decide if we
        # have to put parentheses around s:
        # A regexp might seem better, but I don't think it's really faster.
        # It would be more readable. Python does the below (with in) very quickly.
        if level <= 50:
            t = s[1:]   # ignore leading minus
        else:
            t = s
        # Python complexes are always printed with parentheses
        # we try to avoid double parantheses
        if not PY_TYPE_CHECK_EXACT(o, complex) and \
                (' ' in t or '/' in t or '+' in t or '-' in t or '*' in t \
                or '^' in t):
            s = '(%s)'%s
    return string_from_pystr(s)

cdef public stdstring* py_latex(object o, int level) except +:
    """
    Return latex string representation of o.  If level > 0, possibly
    put parentheses around the string.
    """
    from sage.misc.latex import latex
    s = latex(o)
    if level >= 20:
        if ' ' in s or '/' in s or '+' in s or '-' in s or '*' in s or '^' in s or '\\frac' in s:
            s = '\\left(%s\\right)'%s
    return string_from_pystr(s)

cdef stdstring* string_from_pystr(object py_str):
    """
    Creates a C++ string with the same contents as the given python string.

    Used when passing string output to Pynac for printing, since we don't want
    to mess with reference counts of the python objects and we cannot guarantee
    they won't be garbage collected before the output is printed.

    WARNING: You *must* call this with py_str a str type, or it will segfault.
    """
    cdef char *t_str = PyString_AsString(py_str)
    cdef Py_ssize_t slen = len(py_str)
    cdef stdstring* sout = stdstring_construct_cstr(t_str, slen)
    return sout

cdef public stdstring* py_latex_variable(char* var_name) except +:
    """
    Returns a c++ string containing the latex representation of the given
    variable name.

    Real work is done by the function sage.misc.latex.latex_variable_name.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_latex_variable_for_doctests
        sage: py_latex_variable = py_latex_variable_for_doctests

        sage: py_latex_variable('a')
        a
        sage: py_latex_variable('abc')
        \mbox{abc}
        sage: py_latex_variable('a_00')
        a_{00}
        sage: py_latex_variable('sigma_k')
        \sigma_{k}
        sage: py_latex_variable('sigma389')
        \sigma_{389}
        sage: py_latex_variable('beta_00')
        \beta_{00}
    """
    cdef Py_ssize_t slen
    from sage.misc.latex import latex_variable_name
    py_vlatex = latex_variable_name(var_name)
    return string_from_pystr(py_vlatex)

def py_latex_variable_for_doctests(x):
    """
    Internal function used so we can doctest a certain cdef'd method.

    EXAMPLES::

        sage: sage.symbolic.pynac.py_latex_variable_for_doctests('x')
        x
        sage: sage.symbolic.pynac.py_latex_variable_for_doctests('sigma')
        \sigma
    """
    assert isinstance(x, str)
    cdef stdstring* ostr = py_latex_variable(PyString_AsString(x))
    print(ostr.c_str())
    stdstring_delete(ostr)

def py_print_function_pystring(id, args, fname_paren=False):
    """
    Return a string with the representation of the symbolic function specified
    by the given id applied to args.

    INPUT:

        id --   serial number of the corresponding symbolic function
        params -- Set of parameter numbers with respect to which to take
                    the derivative.
        args -- arguments of the function.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_print_function_pystring, get_ginac_serial
        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: var('x,y,z')
        (x, y, z)
        sage: foo = function('foo', nargs=2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_function_pystring(i, (x,y))
        'foo(x, y)'
        sage: py_print_function_pystring(i, (x,y), True)
        '(foo)(x, y)'
        sage: def my_print(self, *args): return "my args are: " + ', '.join(map(repr, args))
        sage: foo = function('foo', nargs=2, print_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_function_pystring(i, (x,y))
        'my args are: x, y'
    """
    cdef Function func = get_sfunction_from_serial(id)
    # This function is called from two places, from function::print in Pynac
    # and from py_print_fderivative. function::print checks if the serial
    # belongs to a function defined at the C++ level. There are no C++ level
    # functions that return an instance of fderivative when derivated. Hence,
    # func will never be None.
    assert(func is not None)

    # if function has a custom print function call it
    if hasattr(func,'_print_'):
        res = func._print_(*args)
        # make sure the output is a string
        if res is None:
            return ""
        if not isinstance(res, str):
            return str(res)
        return res

    # otherwise use default output
    if fname_paren:
        olist = ['(', func._name, ')']
    else:
        olist = [func._name]
    olist.extend(['(', ', '.join(map(repr, args)), ')'])
    return ''.join(olist)

cdef public stdstring* py_print_function(unsigned id, object args) except +:
    return string_from_pystr(py_print_function_pystring(id, args))

def py_latex_function_pystring(id, args, fname_paren=False):
    """
    Return a string with the latex representation of the symbolic function
    specified by the given id applied to args.

    See documentation of py_print_function_pystring for more information.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_latex_function_pystring, get_ginac_serial
        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: var('x,y,z')
        (x, y, z)
        sage: foo = function('foo', nargs=2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        '{\\rm foo}\\left(x, y^{z}\\right)'
        sage: py_latex_function_pystring(i, (x,y^z), True)
        '\\left({\\rm foo}\\right)\\left(x, y^{z}\\right)'
        sage: py_latex_function_pystring(i, (int(0),x))
        '{\\rm foo}\\left(0, x\\right)'

    Test latex_name::

        sage: foo = function('foo', nargs=2, latex_name=r'\mathrm{bar}')
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        '\\mathrm{bar}\\left(x, y^{z}\\right)'

    Test custom func::

        sage: def my_print(self, *args): return "my args are: " + ', '.join(map(repr, args))
        sage: foo = function('foo', nargs=2, print_latex_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        'my args are: x, y^z'


    """
    cdef Function func = get_sfunction_from_serial(id)
    # This function is called from two places, from function::print in Pynac
    # and from py_latex_fderivative. function::print checks if the serial
    # belongs to a function defined at the C++ level. There are no C++ level
    # functions that return an instance of fderivative when derivated. Hence,
    # func will never be None.
    assert(func is not None)

    # if function has a custom print method call it
    if hasattr(func, '_print_latex_'):
        res = func._print_latex_(*args)
        # make sure the output is a string
        if res is None:
            return ""
        if not isinstance(res, str):
            return str(res)
        return res

    # otherwise, use the latex name if defined
    if func._latex_name:
        name = func._latex_name
    else:
        # if latex_name is not defined, then call
        # latex_variable_name with "is_fname=True" flag
        from sage.misc.latex import latex_variable_name
        name = latex_variable_name(func._name, is_fname=True)
    if fname_paren:
        olist = [r'\left(', name, r'\right)']
    else:
        olist = [name]
    # print the arguments
    from sage.misc.latex import latex
    olist.extend([r'\left(', ', '.join([latex(x) for x in args]),
        r'\right)'] )
    return ''.join(olist)

cdef public stdstring* py_latex_function(unsigned id, object args) except +:
    return string_from_pystr(py_latex_function_pystring(id, args))

cdef public stdstring* py_print_fderivative(unsigned id, object params,
        object args) except +:
    """
    Return a string with the representation of the derivative of the symbolic
    function specified by the given id, lists of params and args.

    INPUT:

        id --   serial number of the corresponding symbolic function
        params -- Set of parameter numbers with respect to which to take
                    the derivative.
        args -- arguments of the function.


    """
    ostr = ''.join(['D[', ', '.join([repr(int(x)) for x in params]), ']'])
    fstr = py_print_function_pystring(id, args, True)
    py_res = ostr + fstr
    return string_from_pystr(py_res)

def py_print_fderivative_for_doctests(id, params, args):
    """
    Used for testing a cdef'd function.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_print_fderivative_for_doctests as py_print_fderivative, get_ginac_serial

        sage: var('x,y,z')
        (x, y, z)
        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: foo = function('foo', nargs=2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1](foo)(x, y^z)

    Test custom print function::

        sage: def my_print(self, *args): return "func_with_args(" + ', '.join(map(repr, args)) +')'
        sage: foo = function('foo', nargs=2, print_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1]func_with_args(x, y^z)

    """
    cdef stdstring* ostr = py_print_fderivative(id, params, args)
    print(ostr.c_str())
    stdstring_delete(ostr)

cdef public stdstring* py_latex_fderivative(unsigned id, object params,
        object args) except +:
    """
    Return a string with the latex representation of the derivative of the
    symbolic function specified by the given id, lists of params and args.

    See documentation of py_print_fderivative for more information.

    """
    ostr = ''.join(['D[', ', '.join([repr(int(x)) for x in params]), ']'])
    fstr = py_latex_function_pystring(id, args, True)
    py_res = ostr + fstr
    return string_from_pystr(py_res)

def py_latex_fderivative_for_doctests(id, params, args):
    """
    Used internally for writing doctests for certain cdef'd functions.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_latex_fderivative_for_doctests as py_latex_fderivative, get_ginac_serial

        sage: var('x,y,z')
        (x, y, z)
        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: foo = function('foo', nargs=2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1]\left({\rm foo}\right)\left(x, y^{z}\right)

    Test latex_name::

        sage: foo = function('foo', nargs=2, latex_name=r'\mathrm{bar}')
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1]\left(\mathrm{bar}\right)\left(x, y^{z}\right)

    Test custom func::

        sage: def my_print(self, *args): return "func_with_args(" + ', '.join(map(repr, args)) +')'
        sage: foo = function('foo', nargs=2, print_latex_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1]func_with_args(x, y^z)
    """
    cdef stdstring* ostr = py_latex_fderivative(id, params, args)
    print(ostr.c_str())
    stdstring_delete(ostr)

#################################################################
# Archive helpers
#################################################################

from sage.structure.sage_object import loads, dumps
cdef public stdstring* py_dumps(object o) except +:
    s = dumps(o, compress=False)
    # pynac archive format terminates atoms with zeroes.
    # since pickle output can break the archive format
    # we use the base64 data encoding
    import base64
    s = base64.b64encode(s)
    return string_from_pystr(s)

cdef public object py_loads(object s) except +:
    import base64
    s = base64.b64decode(s)
    return loads(s)

cdef public object py_get_sfunction_from_serial(unsigned s) except +:
    """
    Return the Python object associated with a serial.
    """
    return get_sfunction_from_serial(s)

cdef public unsigned py_get_serial_from_sfunction(object f) except +:
    """
    Given a Function object return its serial.

    Python's unpickling mechanism is used to unarchive a symbolic function with
    custom methods (evaluation, differentiation, etc.). Pynac extracts a string
    representation from the archive, calls loads() to recreate the stored
    function. This allows us to extract the serial from the Python object to
    set the corresponding member variable of the C++ object representing this
    function.
    """
    return (<Function>f)._serial

cdef public unsigned py_get_serial_for_new_sfunction(stdstring &s,
        unsigned nargs) except +:
    """
    Return a symbolic function with the given name and number of arguments.

    When unarchiving a user defined symbolic function, Pynac goes through
    the registry of existing functions. If there is no function already defined
    with the archived name and number of arguments, this method is called to
    create one and set up the function tables properly.
    """
    from sage.symbolic.function_factory import function_factory
    cdef Function fn = function_factory(s.c_str(), nargs)
    return fn._serial


#################################################################
# Modular helpers
#################################################################

from sage.structure.element cimport Element

cdef public int py_get_parent_char(object o) except -1:
    if isinstance(o, Element):
        return (<Element>o)._parent.characteristic()
    else:
        return 0

#################################################################
# power helpers
#################################################################


from sage.rings.rational cimport rational_power_parts
cdef public object py_rational_power_parts(object base, object exp) except +:
    if not PY_TYPE_CHECK_EXACT(base, Rational):
        base = Rational(base)
    if not PY_TYPE_CHECK_EXACT(exp, Rational):
        exp = Rational(exp)
    res= rational_power_parts(base, exp)
    return res + (bool(res[0] == 1),)

#################################################################
# Binomial Coefficients
#################################################################


cdef public object py_binomial_int(int n, unsigned int k) except +:
    cdef bint sign
    if n < 0:
        n = -n + (k-1)
        sign = k%2
    else:
        sign = 0
    cdef Integer ans = PY_NEW(Integer)
    # Compute the binomial coefficient using GMP.
    mpz_bin_uiui(ans.value, n, k)
    # Return the answer or the negative of it (only if k is odd and n is negative).
    if sign:
        return -ans
    else:
        return ans

cdef public object py_binomial(object n, object k) except +:
    # Keep track of the sign we should use.
    cdef bint sign
    if n < 0:
        n = k-n-1
        sign = k%2
    else:
        sign = 0
    # Convert n and k to unsigned ints.
    cdef unsigned int n_ = n, k_ = k
    cdef Integer ans = PY_NEW(Integer)
    # Compute the binomial coefficient using GMP.
    mpz_bin_uiui(ans.value, n_, k_)
    # Return the answer or the negative of it (only if k is odd and n is negative).
    if sign:
        return -ans
    else:
        return ans

def test_binomial(n, k):
    """
    The Binomial coefficients.  It computes the binomial coefficients.  For
    integer n and k and positive n this is the number of ways of choosing k
    objects from n distinct objects.  If n is negative, the formula
    binomial(n,k) == (-1)^k*binomial(k-n-1,k) is used to compute the result.

    INPUT:
        n, k -- integers, with k >= 0.

    OUTPUT:
        integer

    EXAMPLES::

        sage: import sage.symbolic.pynac
        sage: sage.symbolic.pynac.test_binomial(5,2)
        10
        sage: sage.symbolic.pynac.test_binomial(-5,3)
        -35
        sage: -sage.symbolic.pynac.test_binomial(3-(-5)-1, 3)
        -35
    """
    return py_binomial(n, k)

#################################################################
# GCD
#################################################################
import sage.rings.arith
cdef public object py_gcd(object n, object k) except +:
    if PY_TYPE_CHECK(n, Integer) and PY_TYPE_CHECK(k, Integer):
        if mpz_cmp_si((<Integer>n).value,1) == 0:
            return n
        elif mpz_cmp_si((<Integer>k).value,1) == 0:
            return k
        return n.gcd(k)

    if PY_TYPE_CHECK_EXACT(n, Rational) and PY_TYPE_CHECK_EXACT(k, Rational):
        return n.content(k)
    try:
        return sage.rings.arith.gcd(n,k)
    except (TypeError, ValueError, AttributeError):
        # some strange meaning in case of weird things with no usual lcm.
        return 1


#################################################################
# LCM
#################################################################
cdef public object py_lcm(object n, object k) except +:
    if PY_TYPE_CHECK(n, Integer) and PY_TYPE_CHECK(k, Integer):
        if mpz_cmp_si((<Integer>n).value,1) == 0:
            return k
        elif mpz_cmp_si((<Integer>k).value,1) == 0:
            return n
        return n.lcm(k)
    try:
        return sage.rings.arith.lcm(n,k)
    except (TypeError, ValueError, AttributeError):
        # some strange meaning in case of weird things with no usual lcm, e.g.,
        # elements of finite fields.
        return 1


#################################################################
# Real Part
#################################################################
cdef public object py_real(object x) except +:
    """
    Returns the real part of x.

    TESTS::

        sage: from sage.symbolic.pynac import py_real_for_doctests as py_real
        sage: py_real(I)
        0
        sage: py_real(CC(1,5))
        1.00000000000000
        sage: py_real(CC(1))
        1.00000000000000
        sage: py_real(RR(1))
        1.00000000000000

        sage: py_real(Mod(2,7))
        2

        sage: py_real(QQ['x'].gen())
        x
        sage: py_real(float(2))
        2.0
        sage: py_real(complex(2,2))
        2.0
    """
    if PY_TYPE_CHECK_EXACT(x, float) or PY_TYPE_CHECK_EXACT(x, int) or \
            PY_TYPE_CHECK_EXACT(x, long):
        return x
    elif PY_TYPE_CHECK_EXACT(x, complex):
        return x.real

    try:
        return x.real()
    except AttributeError:
        pass
    try:
        return x.real_part()
    except AttributeError:
        pass

    return x # assume x is real

def py_real_for_doctests(x):
    """
    Used for doctesting py_real.

    TESTS::

        sage: from sage.symbolic.pynac import py_real_for_doctests
        sage: py_real_for_doctests(I)
        0
    """
    return py_real(x)

#################################################################
# Imaginary Part
#################################################################
cdef public object py_imag(object x) except +:
    """
    Return the imaginary part of x.

    TESTS::

        sage: from sage.symbolic.pynac import py_imag_for_doctests as py_imag
        sage: py_imag(I)
        1
        sage: py_imag(CC(1,5))
        5.00000000000000
        sage: py_imag(CC(1))
        0.000000000000000
        sage: py_imag(RR(1))
        0
        sage: py_imag(Mod(2,7))
        0

        sage: py_imag(QQ['x'].gen())
        0
        sage: py_imag(float(2))
        0.0
        sage: py_imag(complex(2,2))
        2.0
    """
    if PY_TYPE_CHECK_EXACT(x, float):
        return float(0)
    if PY_TYPE_CHECK_EXACT(x, complex):
        return x.imag
    try:
        return x.imag()
    except AttributeError:
        pass
    try:
        return x.imag_part()
    except AttributeError:
        pass


    return 0 # assume x is real

def py_imag_for_doctests(x):
    """
    Used for doctesting py_imag.

    TESTS::

        sage: from sage.symbolic.pynac import py_imag_for_doctests
        sage: py_imag_for_doctests(I)
        1
    """
    return py_imag(x)


#################################################################
# Conjugate
#################################################################
cdef public object py_conjugate(object x) except +:
    try:
        return x.conjugate()
    except AttributeError:
        return x # assume is real since it doesn't have an imag attribute.

cdef public bint py_is_rational(object x) except +:
    return PY_TYPE_CHECK_EXACT(x, Rational) or \
           PY_TYPE_CHECK_EXACT(x, Integer) or\
           IS_INSTANCE(x, int) or IS_INSTANCE(x, long)

cdef public bint py_is_equal(object x, object y) except +:
    """
    Return True precisely if x and y are equal.
    """
    return bool(x==y)

cdef public bint py_is_integer(object x) except +:
    r"""
    Returns True if pynac should treat this object as an integer.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_is_integer_for_doctests
        sage: py_is_integer = py_is_integer_for_doctests

        sage: py_is_integer(1r)
        True
        sage: py_is_integer(long(1))
        True
        sage: py_is_integer(3^57)
        True
        sage: py_is_integer(SR(5))
        True
        sage: py_is_integer(4/2)
        True
        sage: py_is_integer(QQbar(sqrt(2))^2)
        True
        sage: py_is_integer(3.0)
        False
        sage: py_is_integer(3.0r)
        False
    """
    return IS_INSTANCE(x, int) or IS_INSTANCE(x, long) or PY_TYPE_CHECK(x, Integer) or \
           (IS_INSTANCE(x, Element) and
            ((<Element>x)._parent.is_exact() or (<Element>x)._parent == ring.SR) and
            (x in ZZ))

def py_is_integer_for_doctests(x):
    """
    Used internally for doctesting purposes.

    TESTS::

        sage: sage.symbolic.pynac.py_is_integer_for_doctests(1r)
        True
        sage: sage.symbolic.pynac.py_is_integer_for_doctests(1/3)
        False
        sage: sage.symbolic.pynac.py_is_integer_for_doctests(2)
        True
    """
    return py_is_integer(x)

cdef public bint py_is_even(object x) except +:
    try:
        return not(x%2)
    except StandardError:
        try:
            return not(ZZ(x)%2)
        except StandardError:
            pass
    return 0


cdef public bint py_is_crational(object x) except +:
    if py_is_rational(x):
        return True
    elif isinstance(x, Element) and (<Element>x)._parent is pynac_I._parent:
        return True
    else:
        return False

def py_is_crational_for_doctest(x):
    """
    Returns True if pynac should treat this object as an element of `\QQ(i)`.

    TESTS::

        sage: from sage.symbolic.pynac import py_is_crational_for_doctest
        sage: py_is_crational_for_doctest(1)
        True
        sage: py_is_crational_for_doctest(-2r)
        True
        sage: py_is_crational_for_doctest(1.5)
        False
        sage: py_is_crational_for_doctest(I.pyobject())
        True
        sage: py_is_crational_for_doctest(I.pyobject()+1/2)
        True
    """
    return py_is_crational(x)

cdef public bint py_is_real(object a) except +:
    if PyInt_CheckExact(a) or PY_TYPE_CHECK(a, Integer) or\
            PyLong_CheckExact(a) or PY_TYPE_CHECK_EXACT(a, float):
        return True
    return py_imag(a) == 0

import sage.rings.arith
cdef public bint py_is_prime(object n) except +:
    try:
        return n.is_prime()
    except StandardError:  # yes, I'm doing this on purpose.
        pass
    try:
        return sage.rings.arith.is_prime(n)
    except StandardError:
        pass
    return False

cdef public object py_numer(object n) except +:
    """
    Return the numerator of the given object. This is called for
    typesetting coefficients.

    TESTS::

        sage: from sage.symbolic.pynac import py_numer_for_doctests as py_numer
        sage: py_numer(2r)
        2
        sage: py_numer(3)
        3
        sage: py_numer(2/3)
        2
        sage: C.<i> = NumberField(x^2+1)
        sage: py_numer(2/3*i)
        2*i
        sage: class no_numer:
        ...     def denominator(self):
        ...         return 5
        ...     def __mul__(left, right):
        ...         return 42
        ...
        sage: py_numer(no_numer())
        42
    """
    if isinstance(n, (int, long, Integer)):
        return n
    try:
        return n.numerator()
    except AttributeError:
        try:
            return n*n.denominator()
        except AttributeError:
            return n

def py_numer_for_doctests(n):
    """
    This function is used to test py_numer().

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_numer_for_doctests
        sage: py_numer_for_doctests(2/3)
        2
    """
    return py_numer(n)

cdef public object py_denom(object n) except +:
    """
    Return the denominator of the given object. This is called for
    typesetting coefficients.

    TESTS::

        sage: from sage.symbolic.pynac import py_denom_for_doctests as py_denom
        sage: py_denom(5)
        1
        sage: py_denom(2/3)
        3
        sage: C.<i> = NumberField(x^2+1)
        sage: py_denom(2/3*i)
        3
    """
    if isinstance(n, (int, long, Integer)):
        return 1
    try:
        return n.denominator()
    except AttributeError:
        return 1

def py_denom_for_doctests(n):
    """
    This function is used to test py_denom().

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_denom_for_doctests
        sage: py_denom_for_doctests(2/3)
        3
    """
    return py_denom(n)

cdef public bint py_is_cinteger(object x) except +:
    return py_is_integer(x) or (py_is_crational(x) and py_denom(x) == 1)

def py_is_cinteger_for_doctest(x):
    """
    Returns True if pynac should treat this object as an element of `\ZZ(i)`.

    TESTS::

        sage: from sage.symbolic.pynac import py_is_cinteger_for_doctest
        sage: py_is_cinteger_for_doctest(1)
        True
        sage: py_is_cinteger_for_doctest(long(-3))
        True
        sage: py_is_cinteger_for_doctest(I.pyobject())
        True
        sage: py_is_cinteger_for_doctest(I.pyobject() - 3)
        True
        sage: py_is_cinteger_for_doctest(I.pyobject() + 1/2)
        False
    """
    return py_is_cinteger(x)

cdef public object py_float(object n, PyObject* parent) except +:
    """
    Evaluate pynac numeric objects numerically.

    TESTS::

        sage: from sage.symbolic.pynac import py_float_for_doctests as py_float
        sage: py_float(I, ComplexField(10))
        1.0*I
        sage: py_float(pi, RealField(100))
        3.1415926535897932384626433833
        sage: py_float(10, CDF)
        10.0
        sage: type(py_float(10, CDF))
        <type 'sage.rings.complex_double.ComplexDoubleElement'>
        sage: py_float(1/2, CC)
        0.500000000000000
        sage: type(py_float(1/2, CC))
        <type 'sage.rings.complex_number.ComplexNumber'>
    """
    if parent is not NULL:
        return (<object>parent)(n)
    else:
        try:
            return RR(n)
        except TypeError:
            return CC(n)

def py_float_for_doctests(n, prec):
    """
    This function is for testing py_float.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_float_for_doctests
        sage: py_float_for_doctests(pi, RealField(80))
        3.1415926535897932384626
    """
    return py_float(n, <PyObject*>prec)

# TODO: Optimize this
from sage.rings.real_double import RDF
cdef public object py_RDF_from_double(double x) except +:
    return RDF(x)

#################################################################
# SPECIAL FUNCTIONS
#################################################################
cdef public object py_tgamma(object x) except +:
    """
    The gamma function exported to pynac.

    TESTS::

        sage: from sage.symbolic.pynac import py_tgamma_for_doctests as py_tgamma
        sage: py_tgamma(4)
        6
        sage: py_tgamma(1/2)
        1.77245385090552
    """
    if PY_TYPE_CHECK_EXACT(x, int) or PY_TYPE_CHECK_EXACT(x, long):
        x = float(x)
    if PY_TYPE_CHECK_EXACT(x, float):
        return sage_tgammal(x)

    # try / except blocks are faster than
    # if hasattr(x, 'gamma')
    try:
        res = x.gamma()
    except AttributeError:
        return CC(x).gamma()

    # the result should be numeric, however the gamma method of rationals may
    # return symbolic expressions. for example (1/2).gamma() -> sqrt(pi).
    if isinstance(res, Expression):
        try:
            return RR(res)
        except ValueError:
            return CC(res)
    return res

def py_tgamma_for_doctests(x):
    """
    This function is for testing py_tgamma().

    TESTS::

        sage: from sage.symbolic.pynac import py_tgamma_for_doctests
        sage: py_tgamma_for_doctests(3)
        2
    """
    return py_tgamma(x)

from sage.rings.arith import factorial
cdef public object py_factorial(object x) except +:
    """
    The factorial function exported to pynac.

    TESTS::

        sage: from sage.symbolic.pynac import py_factorial_py as py_factorial
        sage: py_factorial(4)
        24
        sage: py_factorial(-2/3)
        2.67893853470775
    """
    # factorial(x) is only defined for non-negative integers x
    # so we first test if x can be coerced into ZZ and is non-negative.
    # If this is not the case then we return the symbolic expression gamma(x+1)
    # This fixes Trac 9240
    try:
        x_in_ZZ = ZZ(x)
        coercion_success = True
    except TypeError:
        coercion_success = False

    if coercion_success and x_in_ZZ >= 0:
        return factorial(x)
    else:
        return py_tgamma(x+1)

def py_factorial_py(x):
    """
    This function is a python wrapper around py_factorial(). This wrapper
    is needed when we override the eval() method for GiNaC's factorial
    function in sage.functions.other.Function_factorial.

    TESTS::

        sage: from sage.symbolic.pynac import py_factorial_py
        sage: py_factorial_py(3)
        6
    """
    return py_factorial(x)

cdef public object py_doublefactorial(object x) except +:
    n = Integer(x)
    if n < -1:
        raise ValueError, "argument must be >= -1"
    from sage.misc.misc_c import prod  # fast balanced product
    return prod([n - 2*i for i in range(n//2)])

def doublefactorial(n):
    """
    The double factorial combinatorial function:
        n!! == n * (n-2) * (n-4) * ... * ({1|2}) with 0!! == (-1)!! == 1.

    INPUT:
        n -- an integer > = 1

    EXAMPLES::

        sage: from sage.symbolic.pynac import doublefactorial
        sage: doublefactorial(-1)
        1
        sage: doublefactorial(0)
        1
        sage: doublefactorial(1)
        1
        sage: doublefactorial(5)
        15
        sage: doublefactorial(20)
        3715891200
        sage: prod( [20,18,..,2] )
        3715891200
    """
    return py_doublefactorial(n)


from sage.libs.pari.all import pari
cdef public object py_fibonacci(object n) except +:
    return Integer(pari(n).fibonacci())

cdef public object py_step(object n) except +:
    """
    Return step function of n.
    """
    cdef int c = cmp(n, 0)
    if c < 0:
        return ZERO
    elif c > 0:
        return ONE
    return ONE_HALF

from sage.rings.arith import bernoulli
cdef public object py_bernoulli(object x) except +:
    return bernoulli(x)

cdef public object py_sin(object x) except +:
    """
    TESTS::

        sage: sin(float(2)) #indirect doctest
        0.9092974268256817
        sage: sin(2.)
        0.909297426825682
        sage: sin(2.*I)
        3.62686040784702*I
        sage: sin(QQbar(I))
        1.17520119364380*I
    """
    try:
        return x.sin()
    except AttributeError:
        pass
    try:
        return RR(x).sin()
    except (TypeError, ValueError):
        return CC(x).sin()

cdef public object py_cos(object x) except +:
    """
    TESTS::
        sage: cos(float(2)) #indirect doctest
        -0.4161468365471424
        sage: cos(2.)
        -0.416146836547142
        sage: cos(2.*I)
        3.76219569108363
        sage: cos(QQbar(I))
        1.54308063481524
    """
    try:
        return x.cos()
    except AttributeError:
        pass
    try:
        return RR(x).cos()
    except (TypeError, ValueError):
        return CC(x).cos()

cdef public object py_zeta(object x) except +:
    """
    Return the value of the zeta function at the given value.

    The value is expected to be a numerical object, in RR, CC, RDF or CDF,
    different from 1.

    TESTS::

        sage: from sage.symbolic.pynac import py_zeta_for_doctests as py_zeta
        sage: py_zeta(CC.0)
        0.00330022368532410 - 0.418155449141322*I
        sage: py_zeta(CDF(5))
        1.03692775514
        sage: py_zeta(RealField(100)(5))
        1.0369277551433699263313654865
    """
    try:
        return x.zeta()
    except AttributeError:
        return CC(x).zeta()

def py_zeta_for_doctests(x):
    """
    This function is for testing py_zeta().

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_zeta_for_doctests
        sage: py_zeta_for_doctests(CC.0)
        0.00330022368532410 - 0.418155449141322*I
    """
    return py_zeta(x)

# py_exp(float(1)) should return the value of e
# Note, the nearest IEEE 754 value of e is NOT the same as
# the correctly rounded decimal value.
# The numerical value of e                = 2.71828182845904523536...
# Correctly rounded decimal number        = 2.7182818284590452
# Nearest IEEE 754 format number          = 2.7182818284590451
# On Sun Blade 1000 with SPARC processors = 2.7182818284590455
cdef public object py_exp(object x) except +:
    """
    Return the value of the exp function at the given value.

    The value is expected to be a numerical object, in RR, CC, RDF or CDF.

    TESTS::

        sage: from sage.symbolic.pynac import py_exp_for_doctests as py_exp
        sage: py_exp(CC(1))
        2.71828182845905
        sage: py_exp(CC(.5*I))
        0.877582561890373 + 0.479425538604203*I
        sage: py_exp(float(1))
        2.718281828459045...
        sage: py_exp(QQbar(I))
        0.540302305868140 + 0.841470984807897*I
    """
    if PY_TYPE_CHECK_EXACT(x, float):
        return math.exp(x)
    try:
        return x.exp()
    except AttributeError:
        pass
    try:
        return RR(x).exp()
    except (TypeError, ValueError):
        return CC(x).exp()

def py_exp_for_doctests(x):
    """
    This function tests py_exp.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_exp_for_doctests
        sage: py_exp_for_doctests(CC(2))
        7.38905609893065
    """
    return py_exp(x)

cdef public object py_log(object x) except +:
    """
    Return the value of the log function at the given value.

    The value is expected to be a numerical object, in RR, CC, RDF or CDF.

    TESTS::

        sage: from sage.symbolic.pynac import py_log_for_doctests as py_log
        sage: py_log(CC(e))
        1.00000000000000
        sage: py_log(CC.0)
        1.57079632679490*I
        sage: py_log(float(e))
        1.0
        sage: py_log(float(0))
        -inf
        sage: py_log(float(-1))
        3.141592653589793j
        sage: py_log(int(1))
        0.0
        sage: py_log(long(1))
        0.0
        sage: py_log(int(0))
        -inf
        sage: py_log(long(0))
        -inf
        sage: py_log(complex(0))
        -inf
        sage: py_log(2)
        0.693147180559945
    """
    cdef gsl_complex res
    cdef double real, imag
    if PY_TYPE_CHECK_EXACT(x, int) or PY_TYPE_CHECK_EXACT(x, long):
        x = float(x)
    if PY_TYPE_CHECK_EXACT(x, float):
        if (<float>x) > 0:
            return sage_logl(x)
        elif x < 0:
            res = gsl_complex_log(gsl_complex_rect(PyFloat_AsDouble(x), 0))
            return PyComplex_FromDoubles(res.dat[0], res.dat[1])
        else:
            return float('-inf')
    elif PY_TYPE_CHECK_EXACT(x, complex):
        real = PyComplex_RealAsDouble(x)
        imag = PyComplex_ImagAsDouble(x)
        if real == 0 and imag == 0:
            return float('-inf')
        res = gsl_complex_log(gsl_complex_rect(real, imag))
        return PyComplex_FromDoubles(res.dat[0], res.dat[1])
    elif isinstance(x, Integer):
        return x.log().n()
    elif hasattr(x, 'log'):
        return x.log()
    try:
        return RR(x).log()
    except (TypeError, ValueError):
        return CC(x).log()

def py_log_for_doctests(x):
    """
    This function tests py_log.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_log_for_doctests
        sage: py_log_for_doctests(CC(e))
        1.00000000000000
    """
    return py_log(x)

cdef public object py_tan(object x) except +:
    try:
        return x.tan()
    except AttributeError:
        pass
    try:
        return RR(x).tan()
    except TypeError:
        return CC(x).tan()

cdef public object py_asin(object x) except +:
    try:
        return x.arcsin()
    except AttributeError:
        return RR(x).arcsin()

cdef public object py_acos(object x) except +:
    try:
        return x.arccos()
    except AttributeError:
        return RR(x).arccos()

cdef public object py_atan(object x) except +:
    try:
        return x.arctan()
    except AttributeError:
        return RR(x).arctan()

cdef public object py_atan2(object x, object y) except +:
    """
    Return the value of the two argument arctan function at the given values.

    The values are expected to be numerical objects, for example in RR, CC,
    RDF or CDF.

    Note that the usual call signature of this function has the arguments
    reversed.

    TESTS::

        sage: from sage.symbolic.pynac import py_atan2_for_doctests as py_atan2
        sage: py_atan2(0, 1)
        1.57079632679490
        sage: py_atan2(0.r, 1.r)
        1.5707963267948966
        sage: CC100 = ComplexField(100)
        sage: py_atan2(CC100(0), CC100(1))
        1.5707963267948966192313216916
        sage: RR100 = RealField(100)
        sage: py_atan2(RR100(0), RR100(1))
        1.5707963267948966192313216916
    """
    from sage.symbolic.constants import pi
    parent = parent_c(x)
    assert parent is parent_c(y)
    if parent is ZZ:
        parent = RR
    pi_n = parent(pi)
    cdef int sgn_y = cmp(y, 0)
    cdef int sgn_x = cmp(x, 0)
    if sgn_y:
        if sgn_x > 0:
            return py_atan(abs(y/x)) * sgn_y
        elif sgn_x == 0:
            return pi_n/2 * sgn_y
        else:
            return (pi_n - py_atan(abs(y/x))) * sgn_y
    else:
        if sgn_x > 0:
            return 0
        elif x == 0:
            raise ValueError, "arctan2(0,0) undefined"
        else:
            return pi_n

def py_atan2_for_doctests(x, y):
    """
    Wrapper function to test py_atan2.

    TESTS::

        sage: from sage.symbolic.pynac import py_atan2_for_doctests
        sage: py_atan2_for_doctests(0., 1.)
        1.57079632679490
    """
    return py_atan2(x, y)

cdef public object py_sinh(object x) except +:
    try:
        return x.sinh()
    except AttributeError:
        return RR(x).sinh()


cdef public object py_cosh(object x) except +:
    if PY_TYPE_CHECK_EXACT(x, float):
        return math.cosh(x)
    try:
        return x.cosh()
    except AttributeError:
        return RR(x).cosh()

cdef public object py_tanh(object x) except +:
    try:
        return x.tanh()
    except AttributeError:
        return RR(x).tanh()


cdef public object py_asinh(object x) except +:
    try:
        return x.arcsinh()
    except AttributeError:
        return CC(x).arcsinh()

cdef public object py_acosh(object x) except +:
    try:
        return x.arccosh()
    except AttributeError:
        pass
    try:
        return RR(x).arccosh()
    except TypeError:
        return CC(x).arccosh()


cdef public object py_atanh(object x) except +:
    try:
        return x.arctanh()
    except AttributeError:
        return CC(x).arctanh()

cdef public object py_lgamma(object x) except +:
    """
    Return the value of the log gamma function at the given value.

    The value is expected to be a numerical object, in RR, CC, RDF or CDF.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_lgamma_for_doctests as py_lgamma
        sage: py_lgamma(4)
        1.79175946922805
        sage: py_lgamma(4.r)
        1.791759469228055
        sage: py_lgamma(4r)
        1.791759469228055
        sage: py_lgamma(CC.0)
        -0.650923199301856 - 1.87243664726243*I
        sage: py_lgamma(ComplexField(100).0)
        -0.65092319930185633888521683150 - 1.8724366472624298171188533494*I
    """
    cdef gsl_sf_result lnr, arg
    cdef gsl_complex res
    if PY_TYPE_CHECK_EXACT(x, int) or PY_TYPE_CHECK_EXACT(x, long):
        x = float(x)
    if PY_TYPE_CHECK_EXACT(x, float):
         return sage_lgammal(x)
    elif PY_TYPE_CHECK_EXACT(x, complex):
        gsl_sf_lngamma_complex_e(PyComplex_RealAsDouble(x),PyComplex_ImagAsDouble(x), &lnr, &arg)
        res = gsl_complex_polar(lnr.val, arg.val)
        return PyComplex_FromDoubles(res.dat[0], res.dat[1])
    elif isinstance(x, Integer):
        return x.gamma().log().n()

    # try / except blocks are faster than
    # if hasattr(x, 'log_gamma')
    try:
        return x.log_gamma()
    except AttributeError:
        pass

    try:
        return x.gamma().log()
    except AttributeError:
        pass

    return CC(x).gamma().log()

def py_lgamma_for_doctests(x):
    """
    This function tests py_lgamma.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_lgamma_for_doctests
        sage: py_lgamma_for_doctests(CC(I))
        -0.650923199301856 - 1.87243664726243*I
    """
    return py_lgamma(x)

cdef public object py_isqrt(object x) except +:
    return Integer(x).isqrt()

cdef public object py_sqrt(object x) except +:
    try:
        # WORRY: What if Integer's sqrt calls symbolic one and we go in circle?
        return x.sqrt()
    except AttributeError, msg:
        return sage_sqrtl(float(x))

cdef public object py_abs(object x) except +:
    return abs(x)

cdef public object py_mod(object x, object n) except +:
    """
    Return x mod n. Both x and n are assumed to be integers.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_mod_for_doctests as py_mod
        sage: py_mod(I.parent(5), 4)
        1
        sage: py_mod(3, -2)
        -1
        sage: py_mod(3/2, 2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer


    Note: The original code for this function in GiNaC checks if the arguments
    are integers, and returns 0 otherwise. We omit this check, since all the
    calls to py_mod are preceded by an integer check. We also raise an error
    if converting the arguments to integers fails, since silently returning 0
    would hide possible misuses of this function.

    Regarding the sign of the return value, the CLN reference manual says:

        If x and y are both >= 0, mod(x,y) = rem(x,y) >= 0. In general,
        mod(x,y) has the sign of y or is zero, and rem(x,y) has the sign of
        x or is zero.

    This matches the behavior of the % operator for integers in Sage.
    """
    return Integer(x) % Integer(n)

def py_mod_for_doctests(x, n):
    """
    This function is a python wrapper so py_mod can be tested. The real tests
    are in the docstring for py_mod.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_mod_for_doctests
        sage: py_mod_for_doctests(5, 2)
        1
    """
    return py_mod(x, n)

cdef public object py_smod(object a, object b) except +:
    # Modulus (in symmetric representation).
    # Equivalent to Maple's mods.
    # returns a mod b in the range [-iquo(abs(b)-1,2), iquo(abs(b),2)]
    a = Integer(a); b = Integer(b)
    b = abs(b)
    c = a % b
    if c > b//2:
        c -= b
    return c

cdef public object py_irem(object x, object n) except +:
    return Integer(x) % Integer(n)

cdef public object py_iquo(object x, object n) except +:
    return Integer(x)//Integer(n)

cdef public object py_iquo2(object x, object n) except +:
    x = Integer(x); n = Integer(n)
    try:
        q = x//n
        r = x - q*n
        return q, r
    except (TypeError, ValueError):
        return 0, 0

cdef public int py_int_length(object x) except -1:
    # Size in binary notation.  For integers, this is the smallest n >= 0 such
    # that -2^n <= x < 2^n. If x > 0, this is the unique n > 0 such that
    # 2^(n-1) <= x < 2^n.  This returns 0 if x is not an integer.
    return Integer(x).nbits()

from sage.structure.parent import Parent
cdef public object py_li(object x, object n, object parent) except +:
    """
    Returns a numerical approximation of polylog(n, x) with precision given
    by the ``parent`` argument.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_li_for_doctests as py_li
        sage: py_li(0,2,RR)
        0.000000000000000
        sage: py_li(-1,2,RR)
        -0.822467033424113
        sage: py_li(0, 1, float)
        0.000000000000000
    """
    import mpmath
    if isinstance(parent, Parent) and hasattr(parent, 'prec'):
        prec = parent.prec()
    else:
        prec = 53
    return mpmath_utils.call(mpmath.polylog, n, x, prec=prec)

def py_li_for_doctests(x, n, parent):
    """
    This function is a python wrapper so py_li can be tested. The real tests
    are in the docstring for py_li.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_li_for_doctests
        sage: py_li_for_doctests(0,2,float)
        0.000000000000000
    """
    return py_li(x, n, parent)

cdef public object py_psi(object x) except +:
    """
    EXAMPLES::

        sage: from sage.symbolic.pynac import py_psi_for_doctests as py_psi
        sage: py_psi(0)
        Traceback (most recent call last):
        ...
        ValueError: polygamma pole
        sage: py_psi(1)
        -0.577215664901533
        sage: euler_gamma.n()
        0.577215664901533
    """
    import mpmath
    if PY_TYPE_CHECK(x, Element) and hasattr((<Element>x)._parent, 'prec'):
        prec = (<Element>x)._parent.prec()
    else:
        prec = 53
    return mpmath_utils.call(mpmath.psi, 0, x, prec=prec)

def py_psi_for_doctests(x):
    """
    This function is a python wrapper so py_psi can be tested. The real tests
    are in the docstring for py_psi.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_psi_for_doctests
        sage: py_psi_for_doctests(2)
        0.422784335098467
    """
    return py_psi(x)

cdef public object py_psi2(object n, object x) except +:
    """
    EXAMPLES::

        sage: from sage.symbolic.pynac import py_psi2_for_doctests as py_psi2
        sage: py_psi2(2, 1)
        -2.40411380631919
    """
    import mpmath
    if PY_TYPE_CHECK(x, Element) and hasattr((<Element>x)._parent, 'prec'):
        prec = (<Element>x)._parent.prec()
    else:
        prec = 53
    return mpmath_utils.call(mpmath.psi, n, x, prec=prec)

def py_psi2_for_doctests(n, x):
    """
    This function is a python wrapper so py_psi2 can be tested. The real tests
    are in the docstring for py_psi2.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_psi2_for_doctests
        sage: py_psi2_for_doctests(1, 2)
        0.644934066848226
    """
    return py_psi2(n, x)

cdef public object py_li2(object x) except +:
    """
    EXAMPLES::

        sage: from sage.symbolic.pynac import py_li2_for_doctests as py_li2
        sage: py_li2(-1.1)
        -0.890838090262283
    """
    import mpmath
    if PY_TYPE_CHECK(x, Element) and hasattr((<Element>x)._parent, 'prec'):
        prec = (<Element>x)._parent.prec()
    else:
        prec = 53
    return mpmath_utils.call(mpmath.polylog, 2, x, prec=prec)


def py_li2_for_doctests(x):
    """
    This function is a python wrapper so py_psi2 can be tested. The real tests
    are in the docstring for py_psi2.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_li2_for_doctests
        sage: py_li2_for_doctests(-1.1)
        -0.890838090262283
    """
    return py_li2(x)

##################################################################
# Constants
##################################################################

cdef public GConstant py_get_constant(const_char_ptr name) except +:
    """
    Returns a constant given its name. This is called by
    constant::unarchive in constant.cpp in Pynac and is used for
    pickling.
    """
    from sage.symbolic.constants import constants_name_table
    cdef PynacConstant pc
    c = constants_name_table.get(name, None)
    if c is None:
        raise RuntimeError
    else:
        pc = c._pynac
        return pc.object

cdef public object py_eval_constant(unsigned serial, object parent) except +:
    from sage.symbolic.constants import constants_table
    constant = constants_table[serial]
    return parent(constant)

cdef public object py_eval_unsigned_infinity() except +:
    """
    Returns unsigned_infinity.
    """
    from sage.rings.infinity import unsigned_infinity
    return unsigned_infinity

def py_eval_unsigned_infinity_for_doctests():
    """
    This function tests py_eval_unsigned_infinity.

    TESTS::
        sage: from sage.symbolic.pynac import py_eval_unsigned_infinity_for_doctests as py_eval_unsigned_infinity
        sage: py_eval_unsigned_infinity()
        Infinity
    """
    return py_eval_unsigned_infinity()

cdef public object py_eval_infinity() except +:
    """
    Returns positive infinity, i.e., oo.
    """
    from sage.rings.infinity import infinity
    return infinity

def py_eval_infinity_for_doctests():
    """
    This function tests py_eval_infinity.

    TESTS::
        sage: from sage.symbolic.pynac import py_eval_infinity_for_doctests as py_eval_infinity
        sage: py_eval_infinity()
        +Infinity
    """
    return py_eval_infinity()

cdef public object py_eval_neg_infinity() except +:
    """
    Returns minus_infinity.
    """
    from sage.rings.infinity import minus_infinity
    return minus_infinity

def py_eval_neg_infinity_for_doctests():
    """
    This function tests py_eval_neg_infinity.

    TESTS::
        sage: from sage.symbolic.pynac import py_eval_neg_infinity_for_doctests as py_eval_neg_infinity
        sage: py_eval_neg_infinity()
        -Infinity
    """
    return py_eval_neg_infinity()

##################################################################
# Constructors
##################################################################
cdef Integer z = Integer(0)
cdef public object py_integer_from_long(long x) except +:
    cdef Integer z = PY_NEW(Integer)
    mpz_init_set_si(z.value, x)
    return z

cdef public object py_integer_from_python_obj(object x) except +:
    return Integer(x)


ZERO = ring.SR(0)
ONE = ring.SR(1)
ONE_HALF = ring.SR(Rational((1,2)))

symbol_table = {'functions':{}}
def register_symbol(obj, conversions):
    """
    Add an object to the symbol table, along with how to convert it to
    other systems such as Maxima, Mathematica, etc.  This table is used
    to convert *from* other systems back to Sage.

    INPUTS:

        - `obj` -- a symbolic object or function.

        - `conversions` -- a dictionary of conversions, where the keys
                           are the names of interfaces (e.g.,
                           'maxima'), and the values are the string
                           representation of obj in that system.



    EXAMPLES::

        sage: sage.symbolic.pynac.register_symbol(SR(5),{'maxima':'five'})
        sage: SR(maxima_calculus('five'))
        5
    """
    conversions = dict(conversions)
    try:
        conversions['sage'] = obj.name()
    except AttributeError:
        pass
    for system, value in conversions.iteritems():
        system_table = symbol_table.get(system, None)
        if system_table is None:
            symbol_table[system] = system_table = {}
        system_table[value] = obj



import sage.rings.integer
ginac_pyinit_Integer(sage.rings.integer.Integer)

import sage.rings.real_double
ginac_pyinit_Float(sage.rings.real_double.RDF)

cdef Element pynac_I
I = None

def init_pynac_I():
    """
    Initialize the numeric I object in pynac. We use the generator of QQ(i).

    EXAMPLES::

        sage: sage.symbolic.pynac.init_pynac_I()
        sage: type(sage.symbolic.pynac.I)
        <type 'sage.symbolic.expression.Expression'>
        sage: type(sage.symbolic.pynac.I.pyobject())
        <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>

    TESTS:

    Check that :trac:`10064` is fixed::

        sage: y = I*I*x / x # so y is the expression -1
        sage: y.is_positive()
        False
        sage: z = -x / x
        sage: z.is_positive()
        False
        sage: bool(z == y)
        True
    """
    global pynac_I, I
    from sage.rings.number_field.number_field import QuadraticField
    K = QuadraticField(-1, 'I', embedding=CC.gen(), latex_name='i')
    pynac_I = K.gen()
    ginac_pyinit_I(pynac_I)
    I = new_Expression_from_GEx(ring.SR, g_I)


def init_function_table():
    """
    Initializes the function pointer table in Pynac.  This must be
    called before Pynac is used; otherwise, there will be segfaults.
    """

    py_funcs.py_binomial_int = &py_binomial_int
    py_funcs.py_binomial = &py_binomial
    py_funcs.py_gcd = &py_gcd
    py_funcs.py_lcm = &py_lcm
    py_funcs.py_real = &py_real
    py_funcs.py_imag = &py_imag
    py_funcs.py_numer = &py_numer
    py_funcs.py_denom = &py_denom
    py_funcs.py_conjugate = &py_conjugate

    py_funcs.py_is_rational = &py_is_rational
    py_funcs.py_is_crational = &py_is_crational
    py_funcs.py_is_real = &py_is_real
    py_funcs.py_is_integer = &py_is_integer
    py_funcs.py_is_equal = &py_is_equal
    py_funcs.py_is_even = &py_is_even
    py_funcs.py_is_cinteger = &py_is_cinteger
    py_funcs.py_is_prime = &py_is_prime

    py_funcs.py_integer_from_long = &py_integer_from_long
    py_funcs.py_integer_from_python_obj = &py_integer_from_python_obj

    py_funcs.py_float = &py_float
    py_funcs.py_RDF_from_double = &py_RDF_from_double

    py_funcs.py_factorial = &py_factorial
    py_funcs.py_doublefactorial = &py_doublefactorial
    py_funcs.py_fibonacci = &py_fibonacci
    py_funcs.py_step = &py_step
    py_funcs.py_bernoulli = &py_bernoulli
    py_funcs.py_sin = &py_sin
    py_funcs.py_cos = &py_cos
    py_funcs.py_zeta = &py_zeta
    py_funcs.py_exp = &py_exp
    py_funcs.py_log = &py_log
    py_funcs.py_tan = &py_tan
    py_funcs.py_asin = &py_asin
    py_funcs.py_acos = &py_acos
    py_funcs.py_atan = &py_atan
    py_funcs.py_atan2 = &py_atan2
    py_funcs.py_sinh = &py_sinh
    py_funcs.py_cosh = &py_cosh
    py_funcs.py_tanh = &py_tanh
    py_funcs.py_asinh = &py_asinh
    py_funcs.py_acosh = &py_acosh
    py_funcs.py_atanh = &py_atanh
    py_funcs.py_tgamma = &py_tgamma
    py_funcs.py_lgamma = &py_lgamma
    py_funcs.py_isqrt = &py_isqrt
    py_funcs.py_sqrt = &py_sqrt
    py_funcs.py_abs = &py_abs
    py_funcs.py_mod = &py_mod
    py_funcs.py_smod = &py_smod
    py_funcs.py_irem = &py_irem
    py_funcs.py_iquo = &py_iquo
    py_funcs.py_iquo2 = &py_iquo2
    py_funcs.py_li = &py_li
    py_funcs.py_li2 = &py_li2
    py_funcs.py_psi = &py_psi
    py_funcs.py_psi2 = &py_psi2

    py_funcs.py_int_length = &py_int_length

    py_funcs.py_eval_constant = &py_eval_constant
    py_funcs.py_eval_unsigned_infinity = &py_eval_unsigned_infinity
    py_funcs.py_eval_infinity = &py_eval_infinity
    py_funcs.py_eval_neg_infinity = &py_eval_neg_infinity

    py_funcs.py_get_parent_char = &py_get_parent_char

    py_funcs.py_latex = &py_latex
    py_funcs.py_repr = &py_repr

    py_funcs.py_dumps = &py_dumps
    py_funcs.py_loads = &py_loads

    py_funcs.exvector_to_PyTuple = &exvector_to_PyTuple
    py_funcs.pyExpression_to_ex = &pyExpression_to_ex
    py_funcs.ex_to_pyExpression = &ex_to_pyExpression
    py_funcs.py_print_function = &py_print_function
    py_funcs.py_latex_function = &py_latex_function
    py_funcs.py_get_ginac_serial = &py_get_ginac_serial
    py_funcs.py_get_sfunction_from_serial = &py_get_sfunction_from_serial
    py_funcs.py_get_serial_from_sfunction = &py_get_serial_from_sfunction
    py_funcs.py_get_serial_for_new_sfunction = &py_get_serial_for_new_sfunction

    py_funcs.py_get_constant = &py_get_constant
    py_funcs.py_print_fderivative =  &py_print_fderivative
    py_funcs.py_latex_fderivative =  &py_latex_fderivative
    py_funcs.paramset_to_PyTuple = &paramset_to_PyTuple
    py_funcs.py_rational_power_parts = &py_rational_power_parts

init_function_table()
init_pynac_I()

"""
Some tests for the formal square root of -1.

EXAMPLES::

    sage: I
    I
    sage: I^2
    -1

Note that conversions to real fields will give TypeErrors::

    sage: float(I)
    Traceback (most recent call last):
    ...
    TypeError: unable to simplify to float approximation
    sage: gp(I)
    I
    sage: RR(I)
    Traceback (most recent call last):
    ...
    TypeError: Unable to convert x (='1.00000000000000*I') to real number.

We can convert to complex fields::

    sage: C = ComplexField(200); C
    Complex Field with 200 bits of precision
    sage: C(I)
    1.0000000000000000000000000000000000000000000000000000000000*I
    sage: I._complex_mpfr_field_(ComplexField(53))
    1.00000000000000*I

    sage: I._complex_double_(CDF)
    1.0*I
    sage: CDF(I)
    1.0*I

    sage: z = I + I; z
    2*I
    sage: C(z)
    2.0000000000000000000000000000000000000000000000000000000000*I
    sage: 1e8*I
    1.00000000000000e8*I

    sage: complex(I)
    1j

    sage: QQbar(I)
    1*I

    sage: abs(I)
    1

    sage: I.minpoly()
    x^2 + 1
    sage: maxima(2*I)
    2*%i

TESTS::

    sage: repr(I)
    'I'
    sage: latex(I)
    i
"""


set_ginac_fn_serial()
