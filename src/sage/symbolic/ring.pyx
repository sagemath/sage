# -*- coding: utf-8 -*-
"""
The symbolic ring
"""

# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.ext.cplusplus cimport ccrepr

from sage.libs.pynac.pynac cimport *

from sage.rings.integer cimport Integer
from sage.rings.real_mpfr cimport RealNumber

from sage.symbolic.expression cimport Expression, new_Expression_from_GEx, new_Expression_from_pyobject, is_Expression

from sage.misc.latex import latex_variable_name
from sage.cpython.string cimport str_to_bytes, bytes_to_str, char_to_str
from sage.structure.element cimport RingElement, Element, Matrix
from sage.categories.morphism cimport Morphism
from sage.structure.coerce cimport is_numpy_type

from sage.rings.all import RR, CC, ZZ

import keyword
import operator

# Do not allow any of these keywords as identifiers for symbolic variables
KEYWORDS = set(keyword.kwlist).union(['exec', 'print', 'None', 'True',
                                      'False', 'nonlocal'])


cdef class SymbolicRing(CommutativeRing):
    """
    Symbolic Ring, parent object for all symbolic expressions.
    """
    def __init__(self, base_ring = None):
        """
        Initialize the Symbolic Ring.

        EXAMPLES::

            sage: SR
            Symbolic Ring

        TESTS::

            sage: isinstance(SR, sage.symbolic.ring.SymbolicRing)
            True
            sage: TestSuite(SR).run(skip=['_test_divides'])

        """
        if base_ring is None:
            base_ring = self
        CommutativeRing.__init__(self, base_ring)
        self._populate_coercion_lists_(convert_method_name='_symbolic_')
        self.symbols = {}

    def __reduce__(self):
        """
        EXAMPLES::

           sage: loads(dumps(SR)) == SR           # indirect doctest
           True
        """
        return the_SymbolicRing, tuple([])

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLES::

            sage: repr(SR)
            'Symbolic Ring'
        """
        return "Symbolic Ring"

    def _latex_(self):
        """
        Return latex representation of the symbolic ring.

        EXAMPLES::

            sage: latex(SR)
            \text{SR}
            sage: M = MatrixSpace(SR, 2); latex(M)
            \mathrm{Mat}_{2\times 2}(\text{SR})
        """
        return r'\text{SR}'

    cpdef _coerce_map_from_(self, R):
        """
        EXAMPLES::

            sage: SR.coerce(int(2))
            2
            sage: SR.coerce(-infinity)
            -Infinity
            sage: SR.coerce(unsigned_infinity)
            Infinity
            sage: SR.has_coerce_map_from(ZZ['t'])
            True
            sage: SR.has_coerce_map_from(ZZ['t,u,v'])
            True
            sage: SR.has_coerce_map_from(Frac(ZZ['t,u,v']))
            True
            sage: SR.has_coerce_map_from(GF(5)['t'])
            True
            sage: SR.has_coerce_map_from(SR['t'])
            False
            sage: SR.has_coerce_map_from(Integers(8))
            True
            sage: SR.has_coerce_map_from(GF(9, 'a'))
            True
            sage: SR.has_coerce_map_from(RealBallField())
            True
            sage: SR.has_coerce_map_from(ComplexBallField())
            True
            sage: SR.has_coerce_map_from(UnsignedInfinityRing)
            True

        TESTS::

            sage: SR.has_coerce_map_from(pari)
            False

        Check if arithmetic with bools works (see :trac:`9560`)::

            sage: SR.has_coerce_map_from(bool)
            True
            sage: SR(5)*True; True*SR(5)
            5
            5
            sage: SR(5)+True; True+SR(5)
            6
            6
            sage: SR(5)-True
            4

        TESTS::

            sage: SR.has_coerce_map_from(SR.subring(accepting_variables=('a',)))
            True
            sage: SR.has_coerce_map_from(SR.subring(rejecting_variables=('r',)))
            True
            sage: SR.has_coerce_map_from(SR.subring(no_variables=True))
            True

            sage: SR.has_coerce_map_from(AA)
            True
            sage: SR.has_coerce_map_from(QQbar)
            True
        """
        if isinstance(R, type):
            if R in [int, float, long, complex, bool]:
                return True

            if is_numpy_type(R):
                import numpy
                if (issubclass(R, numpy.integer) or
                    issubclass(R, numpy.floating) or
                    issubclass(R, numpy.complexfloating)):
                    return NumpyToSRMorphism(R)
                else:
                    return None

            if 'sympy' in R.__module__:
                from sympy.core.basic import Basic
                if issubclass(R, Basic):
                    return UnderscoreSageMorphism(R, self)

            return False
        else:
            from sage.rings.real_mpfr import mpfr_prec_min

            from sage.rings.fraction_field import is_FractionField
            from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing
            from sage.rings.real_mpfi import is_RealIntervalField
            from sage.rings.complex_interval_field import is_ComplexIntervalField
            from sage.rings.real_arb import RealBallField
            from sage.rings.complex_arb import ComplexBallField
            from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
            from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
            from sage.rings.polynomial.laurent_polynomial_ring import is_LaurentPolynomialRing

            from sage.rings.all import (ComplexField,
                                        RLF, CLF,
                                        InfinityRing,
                                        UnsignedInfinityRing)
            from sage.rings.finite_rings.finite_field_base import is_FiniteField

            from sage.interfaces.maxima import Maxima

            from .subring import GenericSymbolicSubring

            if R._is_numerical():
                # Almost anything with a coercion into any precision of CC
                return R not in (RLF, CLF)
            elif is_PolynomialRing(R) or is_MPolynomialRing(R) or is_FractionField(R) or is_LaurentPolynomialRing(R):
                base = R.base_ring()
                return base is not self and self.has_coerce_map_from(base)
            elif (R is InfinityRing or R is UnsignedInfinityRing
                  or is_RealIntervalField(R) or is_ComplexIntervalField(R)
                  or isinstance(R, RealBallField)
                  or isinstance(R, ComplexBallField)
                  or is_IntegerModRing(R) or is_FiniteField(R)):
                return True
            elif isinstance(R, GenericSymbolicSubring):
                return True

    def _element_constructor_(self, x):
        r"""
        Coerce `x` into the symbolic expression ring SR.

        EXAMPLES::

            sage: a = SR(-3/4); a
            -3/4
            sage: type(a)
            <type 'sage.symbolic.expression.Expression'>
            sage: a.parent()
            Symbolic Ring
            sage: K.<a> = QuadraticField(-3)
            sage: a + sin(x)
            I*sqrt(3) + sin(x)
            sage: x=var('x'); y0,y1=PolynomialRing(ZZ,2,'y').gens()
            sage: x+y0/y1
            x + y0/y1
            sage: x.subs(x=y0/y1)
            y0/y1
            sage: x + int(1)
            x + 1

        If `a` is already in the symbolic expression ring, coercing returns
        `a` itself (not a copy)::

            sage: a = SR(-3/4); a
            -3/4
            sage: SR(a) is a
            True

        A Python complex number::

            sage: SR(complex(2,-3))
            (2-3j)

        TESTS::

            sage: SR._coerce_(int(5))
            5
            sage: SR._coerce_(5)
            5
            sage: SR._coerce_(float(5))
            5.0
            sage: SR._coerce_(5.0)
            5.00000000000000

        An interval arithmetic number::

            sage: SR._coerce_(RIF(pi))
            3.141592653589794?

        The complex number `I`::

            sage: si = SR.coerce(I)
            sage: si^2
            -1
            sage: bool(si == CC.0)
            True

        Polynomial ring element factorizations::

            sage: R.<x> = QQ[]
            sage: SR(factor(5*x^2 - 5))
            5*(x + 1)*(x - 1)
            sage: R.<x,y> = QQ[]
            sage: SR(factor(x^2 - y^2))
            (x + y)*(x - y)
            sage: R.<x,y,z> = QQ[]
            sage: SR(factor(x^2*y^3 + x^2*y^2*z - x*y^3 - x*y^2*z - 2*x*y*z - 2*x*z^2 + 2*y*z + 2*z^2))
            (x*y^2 - 2*z)*(x - 1)*(y + z)

        Asymptotic expansions::

            sage: A.<x, y> = AsymptoticRing(growth_group='x^ZZ * y^QQ * log(y)^ZZ', coefficient_ring=ZZ)
            sage: s = SR(3*x^5 * log(y) + 4*y^(3/7) + O(x*log(y))); s
            3*x^5*log(y) + 4*y^(3/7) + Order(x*log(y))
            sage: s.operator(), s.operands()
            (<function add_vararg at 0x...>,
             [3*x^5*log(y), 4*y^(3/7), Order(x*log(y))])
            sage: t = s.operands()[0]; t
            3*x^5*log(y)
            sage: t.operator(), t.operands()
            (<function mul_vararg at 0x...>, [x^5, log(y), 3])

        We get a sensible error message if conversion fails::

            sage: SR(int)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert <... 'int'> to a symbolic expression
            sage: r^(1/2)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for ** or pow(): 'R' and 'sage.rings.rational.Rational'

        Check that :trac:`22068` is fixed::

            sage: _ = var('x')
            sage: sin(x).subs(x=RR('NaN'))
            sin(NaN)
            sage: SR(RR('NaN')).is_real()
            False
            sage: sin(x).subs(x=float('NaN'))
            sin(NaN)
            sage: SR(float('NaN')).is_real()
            False
            sage: sin(x).subs(x=complex('NaN'))
            sin(NaN)

        Check that :trac:`24072` is solved::

            sage: x = polygen(GF(3))
            sage: a = SR.var('a')
            sage: (2*x + 1) * a
            Traceback (most recent call last):
            ...
            TypeError: positive characteristic not allowed in symbolic computations

        Check support for unicode characters (:trac:`29280`)::

            sage: SR('λ + 2λ')
            3*λ
            sage: SR('μ') is var('μ')
            True
            sage: SR('λ + * 1')
            Traceback (most recent call last):
            ...
            TypeError: Malformed expression: λ + * !!!  1
        """
        cdef GEx exp
        if is_Expression(x):
            return new_Expression_from_GEx(self, (<Expression>x)._gobj)
        if hasattr(x, '_symbolic_'):
            return x._symbolic_(self)
        elif isinstance(x, str):
            try:
                from sage.calculus.calculus import symbolic_expression_from_string
                return self(symbolic_expression_from_string(x))
            except SyntaxError as err:
                msg, s, pos = err.args
                raise TypeError("%s: %s !!! %s" % (msg, s[:pos], s[pos:]))

        from sage.rings.infinity import (infinity, minus_infinity,
                                         unsigned_infinity)
        from sage.structure.factorization import Factorization

        if isinstance(x, RealNumber):
            if x.is_NaN():
                from sage.symbolic.constants import NaN
                return NaN
            exp = x
        elif isinstance(x, (float, complex)):
            if not (x == x):
                from sage.symbolic.constants import NaN
                return NaN
            exp = x
        elif isinstance(x, long):
            exp = x
        elif isinstance(x, int):
            exp = GEx(<long>x)
        elif x is infinity:
            return new_Expression_from_GEx(self, g_Infinity)
        elif x is minus_infinity:
            return new_Expression_from_GEx(self, g_mInfinity)
        elif x is unsigned_infinity:
            return new_Expression_from_GEx(self, g_UnsignedInfinity)
        elif isinstance(x, (RingElement, Matrix)):
            if x.parent().characteristic():
                raise TypeError('positive characteristic not allowed in symbolic computations')
            exp = x
        elif isinstance(x, Factorization):
            from sage.misc.all import prod
            return prod([SR(p)**e for p,e in x], SR(x.unit()))
        else:
            raise TypeError(f"unable to convert {x!r} to a symbolic expression")

        return new_Expression_from_GEx(self, exp)

    def _force_pyobject(self, x, bint force=False, bint recursive=True):
        """
        Wrap the given Python object in a symbolic expression even if it
        cannot be coerced to the Symbolic Ring.

        INPUT:

        - ``x`` - a Python object.

        - ``force`` - bool, default ``False``, if True, the Python object
          is taken as is without attempting coercion or list traversal.

        - ``recursive`` - bool, default ``True``, disables recursive
          traversal of lists.

        EXAMPLES::

            sage: t = SR._force_pyobject(QQ); t
            Rational Field
            sage: type(t)
            <type 'sage.symbolic.expression.Expression'>

        Testing tuples::

            sage: t = SR._force_pyobject((1, 2, x, x+1, x+2)); t
            (1, 2, x, x + 1, x + 2)
            sage: t.subs(x = 2*x^2)
            (1, 2, 2*x^2, 2*x^2 + 1, 2*x^2 + 2)
            sage: t.op[0]
            1
            sage: t.op[2]
            x

        It also works if the argument is a ``list``::

            sage: t = SR._force_pyobject([1, 2, x, x+1, x+2]); t
            (1, 2, x, x + 1, x + 2)
            sage: t.subs(x = 2*x^2)
            (1, 2, 2*x^2, 2*x^2 + 1, 2*x^2 + 2)
            sage: SR._force_pyobject((QQ, RR, CC))
            (Rational Field, Real Field with 53 bits of precision, Complex Field with 53 bits of precision)
            sage: t = SR._force_pyobject((QQ, (x, x + 1, x + 2), CC)); t
            (Rational Field, (x, x + 1, x + 2), Complex Field with 53 bits of precision)
            sage: t.subs(x=x^2)
            (Rational Field, (x^2, x^2 + 1, x^2 + 2), Complex Field with 53 bits of precision)

        If ``recursive`` is ``False`` the inner tuple is taken as a Python
        object. This prevents substitution as above::

            sage: t = SR._force_pyobject((QQ, (x, x + 1, x + 2), CC), recursive=False)
            sage: t
            (Rational Field, (x, x + 1, x + 2), Complex Field with 53 bits
            of precision)
            sage: t.subs(x=x^2)
            (Rational Field, (x, x + 1, x + 2), Complex Field with 53 bits
            of precision)
        """
        cdef GEx exp
        cdef GExprSeq ex_seq
        cdef GExVector ex_v
        if force:
            exp = x

        else:
            # first check if we can do it the nice way
            if isinstance(x, Expression):
                return x
            try:
                return self._coerce_(x)
            except TypeError:
                pass

            # tuples can be packed into exprseq
            if isinstance(x, (tuple, list)):
                for e in x:
                    obj = SR._force_pyobject(e, force=(not recursive))
                    ex_v.push_back( (<Expression>obj)._gobj )

                ex_seq = GExprSeq(ex_v)

                exp = GEx(ex_seq)
            else:
                exp = x

        return new_Expression_from_GEx(self, exp)

    def wild(self, unsigned int n=0):
        """
        Return the n-th wild-card for pattern matching and substitution.

        INPUT:

        - ``n`` - a nonnegative integer

        OUTPUT:

        - `n^{th}` wildcard expression

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: w0 = SR.wild(0); w1 = SR.wild(1)
            sage: pattern = sin(x)*w0*w1^2; pattern
            $1^2*$0*sin(x)
            sage: f = atan(sin(x)*3*x^2); f
            arctan(3*x^2*sin(x))
            sage: f.has(pattern)
            True
            sage: f.subs(pattern == x^2)
            arctan(x^2)

        TESTS:

        Check that :trac:`15047` is fixed::

            sage: latex(SR.wild(0))
            \$0

        Check that :trac:`21455` is fixed::

            sage: coth(SR.wild(0))
            coth($0)
        """
        return new_Expression_from_GEx(self, g_wild(n))

    def __contains__(self, x):
        r"""
        True if there is an element of the symbolic ring that is equal to x
        under ``==``.

        EXAMPLES:

        The symbolic variable x is in the symbolic ring.::

            sage: x.parent()
            Symbolic Ring
            sage: x in SR
            True

        2 is also in the symbolic ring since it is equal to something in
        SR, even though 2's parent is not SR.

        ::

            sage: 2 in SR
            True
            sage: parent(2)
            Integer Ring
            sage: 1/3 in SR
            True
        """
        try:
            x2 = self(x)
            return bool(x2 == x)
        except TypeError:
            return False

    def characteristic(self):
        """
        Return the characteristic of the symbolic ring, which is 0.

        OUTPUT:

        - a Sage integer

        EXAMPLES::

            sage: c = SR.characteristic(); c
            0
            sage: type(c)
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(0)

    def _an_element_(self):
        """
        Return an element of the symbolic ring, which is used by the
        coercion model.

        EXAMPLES::

            sage: SR._an_element_()
            some_variable
        """
        return self.symbol('some_variable')

    def is_field(self, proof = True):
        """
        Returns True, since the symbolic expression ring is (for the most
        part) a field.

        EXAMPLES::

            sage: SR.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return False, since the Symbolic Ring is infinite.

        EXAMPLES::

            sage: SR.is_finite()
            False
        """
        return False

    cpdef bint is_exact(self) except -2:
        """
        Return False, because there are approximate elements in the
        symbolic ring.

        EXAMPLES::

            sage: SR.is_exact()
            False

        Here is an inexact element.

        ::

            sage: SR(1.9393)
            1.93930000000000
        """
        return False

    def pi(self):
        """
        EXAMPLES::

            sage: SR.pi() is pi
            True
        """
        from sage.symbolic.constants import pi
        return self(pi)

    def I(self):
        r"""
        The imaginary unit, viewed as an element of the symbolic ring.

        EXAMPLES::

            sage: SR.I()^2
            -1
            sage: SR.I().parent()
            Symbolic Ring
        """
        from sage.symbolic.constants import I
        return I

    cpdef Expression symbol(self, name=None, latex_name=None, domain=None):
        """
        EXAMPLES::

            sage: t0 = SR.symbol("t0")
            sage: t0.conjugate()
            conjugate(t0)

            sage: t1 = SR.symbol("t1", domain='real')
            sage: t1.conjugate()
            t1

            sage: t0.abs()
            abs(t0)

            sage: t0_2 = SR.symbol("t0", domain='positive')
            sage: t0_2.abs()
            t0
            sage: bool(t0_2 == t0)
            True
            sage: t0.conjugate()
            t0

            sage: SR.symbol() # temporary variable
            symbol...

        We propagate the domain to the assumptions database::

            sage: n = var('n', domain='integer')
            sage: solve([n^2 == 3],n)
            []

        TESTS:

        Test that the parent is set correctly (inheritance)::

            sage: from sage.symbolic.ring import SymbolicRing
            sage: class MySymbolicRing(SymbolicRing):
            ....:     def _repr_(self):
            ....:         return 'My Symbolic Ring'
            sage: MySR = MySymbolicRing()
            sage: MySR.symbol('x').parent()
            My Symbolic Ring
            sage: MySR.var('x').parent()  # indirect doctest
            My Symbolic Ring
            sage: MySR.var('blub').parent()  # indirect doctest
            My Symbolic Ring
            sage: MySR.an_element().parent()
            My Symbolic Ring
        """
        cdef GSymbol symb
        cdef Expression e

        # check if there is already a symbol with same name
        e = self.symbols.get(name)

        # fast path to get an already existing variable
        if e is not None:
            if domain is None:
                if latex_name is None:
                    return e

            # get symbol
            symb = ex_to_symbol(e._gobj)
            if latex_name is not None:
                symb.set_texname(str_to_bytes(latex_name))
            if domain is not None:
                symb.set_domain(sage_domain_to_ginac_domain(domain))
            e._gobj = GEx(symb)
            if domain is not None:
                send_sage_domain_to_maxima(e, domain)

            return e

        else: # initialize a new symbol
            # Construct expression
            e = <Expression>Expression.__new__(Expression)
            e._parent = self

            if name is None: # Check if we need a temporary anonymous new symbol
                symb = ginac_new_symbol()
                if domain is not None:
                    symb.set_domain(sage_domain_to_ginac_domain(domain))
            else:
                if latex_name is None:
                    latex_name = latex_variable_name(name)
                if domain is not None:
                    ginac_domain = sage_domain_to_ginac_domain(domain)
                else:
                    ginac_domain = domain_complex
                symb = ginac_symbol(str_to_bytes(name),
                                    str_to_bytes(latex_name), ginac_domain)
                self.symbols[name] = e

            e._gobj = GEx(symb)
            if domain is not None:
                send_sage_domain_to_maxima(e, domain)

        return e

    def var(self, name, latex_name=None, n=None, domain=None):
        """
        Return a symbolic variable as an element of the symbolic ring.

        INPUT:

        - ``name`` -- string or list of strings with the name(s) of the symbolic variable(s)

        - ``latex_name`` -- (optional) string used when printing in latex mode, if not specified use ``'name'``

        - ``n`` -- (optional) positive integer; number of symbolic variables, indexed from `0` to `n-1`

        - ``domain`` -- (optional) specify the domain of the variable(s); it is the complex plane
          by default, and possible options are (non-exhaustive list, see note below):
          ``'real'``, ``'complex'``, ``'positive'``, ``'integer'`` and ``'noninteger'``

        OUTPUT:

        Symbolic expression or tuple of symbolic expressions.

        .. SEEALSO::

            This function does not inject the variable(s) into the global namespace.
            For that purpose see :meth:`var()<sage.calculus.var.var>`.

        .. NOTE::

            For a comprehensive list of acceptable features type ``'maxima('features')'``,
            and see also the documentation of :ref:`sage.symbolic.assumptions`.

        EXAMPLES:

        Create a variable `zz` (complex by default)::

            sage: zz = SR.var('zz'); zz
            zz

        The return type is a symbolic expression::

            sage: type(zz)
            <type 'sage.symbolic.expression.Expression'>

        We can specify the domain as well::

            sage: zz = SR.var('zz', domain='real')
            sage: zz.is_real()
            True

        The real domain is also set with the integer domain::

            sage: SR.var('x', domain='integer').is_real()
            True

        The ``name`` argument does not have to match the left-hand side variable::

            sage: t = SR.var('theta2'); t
            theta2

        Automatic indexing is available as well::

            sage: x = SR.var('x', 4)
            sage: x[0], x[3]
            (x0, x3)
            sage: sum(x)
            x0 + x1 + x2 + x3

        TESTS::

            sage: var(' x y  z    ')
            (x, y, z)
            sage: var(' x  ,  y ,  z    ')
            (x, y, z)
            sage: var(' ')
            Traceback (most recent call last):
            ...
            ValueError: You need to specify the name of the new variable.

            var(['x', 'y ', ' z '])
            (x, y, z)
            var(['x,y'])
            Traceback (most recent call last):
            ...
            ValueError: The name "x,y" is not a valid Python identifier.

        Check that :trac:`17206` is fixed::

            sage: var1 = var('var1', latex_name=r'\sigma^2_1'); latex(var1)
            {\sigma^2_1}

        The number of variables should be an integer greater or equal than 1::

            sage: SR.var('K', -273)
            Traceback (most recent call last):
            ...
            ValueError: the number of variables should be a positive integer

        The argument ``n`` can only handle a single variable::

            sage: SR.var('x y', 4)
            Traceback (most recent call last):
            ...
            ValueError: cannot specify n for multiple symbol names

        Check that :trac:`28353` is fixed: Constructions that suggest multiple
        variables but actually only give one variable name return a 1-tuple::

            sage: SR.var(['x'])
            (x,)
            sage: SR.var('x,')
            (x,)
            sage: SR.var(['x'], n=4)
            Traceback (most recent call last):
            ...
            ValueError: cannot specify n for multiple symbol names
        """
        if is_Expression(name):
            return name
        if not isinstance(name, (basestring, list, tuple)):
            name = repr(name)

        is_multiple = False

        if isinstance(name, (list, tuple)):
            names_list = [s.strip() for s in name]
            is_multiple = True
        else:
            name = name.strip()
            if ',' in name:
                names_list = [s.strip() for s in name.split(',') if s.strip()]
                is_multiple = True
            elif ' ' in name:
                names_list = [s.strip() for s in name.split()]
                is_multiple = True
            else:
                names_list = [name] if name else []

        for s in names_list:
            if not isidentifier(s):
                raise ValueError(f'The name "{s}" is not a valid Python identifier.')

        formatted_latex_name = None
        if latex_name is not None and n is None:
            try:
                n = operator.index(latex_name)
                latex_name = None
            except TypeError:
                formatted_latex_name = '{{{0}}}'.format(latex_name)

        if not names_list:
            raise ValueError('You need to specify the name of the new variable.')

        if is_multiple:
            if latex_name is not None:
                raise ValueError("cannot specify latex_name for multiple symbol names")
            if n is not None:
                raise ValueError("cannot specify n for multiple symbol names")
            return tuple([self.symbol(s, domain=domain) for s in names_list])
        else:
            if n is not None:
                if n > 0:
                    name = [name + str(i) for i in range(n)]
                    if latex_name is None:
                        return tuple([self.symbol(name[i], domain=domain) for i in range(n)])
                    else:
                        formatted_latex_name = ['{{{}}}_{{{}}}'.format(latex_name, str(i)) for i in range(n)]
                        return tuple([self.symbol(name[i], latex_name=formatted_latex_name[i], domain=domain) for i in range(n)])
                else:
                    raise ValueError("the number of variables should be a positive integer")
            else:
                return self.symbol(name, latex_name=formatted_latex_name, domain=domain)

    def _repr_element_(self, Expression x):
        """
        Returns the string representation of the element x.  This is
        used so that subclasses of the SymbolicRing (such the a
        CallableSymbolicExpressionRing) can provide their own
        implementations of how to print Expressions.

        EXAMPLES::

            sage: SR._repr_element_(x+2)
            'x + 2'
        """
        return ccrepr(x._gobj)

    def _latex_element_(self, Expression x):
        r"""
        Returns the standard LaTeX version of the expression `x`.

        EXAMPLES::

            sage: latex(sin(x+2))
            \sin\left(x + 2\right)
            sage: latex(var('theta') + 2)
            \theta + 2
        """
        return char_to_str(GEx_to_str_latex(&x._gobj))

    def _call_element_(self, _the_element, *args, **kwds):
        """
        EXAMPLES::

            sage: x,y=var('x,y')
            sage: f = x+y
            sage: f.variables()
            (x, y)
            sage: f()
            x + y
            sage: f(3)
            doctest:...: DeprecationWarning: Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)
            See http://trac.sagemath.org/5930 for details.
            y + 3
            sage: f(x=3)
            y + 3
            sage: f(3,4)
            7
            sage: f(x=3,y=4)
            7
            sage: f(2,3,4)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 2
            sage: f(x=2,y=3,z=4)
            5

        ::

            sage: f({x:3})
            y + 3
            sage: f({x:3,y:4})
            7
            sage: f(x=3)
            y + 3
            sage: f(x=3,y=4)
            7

        ::

            sage: a = (2^(8/9))
            sage: a(4)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 0


        Note that you make get unexpected results when calling
        symbolic expressions and not explicitly giving the variables::

            sage: f = function('Gamma')(var('z'), var('w')); f
            Gamma(z, w)
            sage: f(2)
            Gamma(z, 2)
            sage: f(2,5)
            Gamma(5, 2)

        Thus, it is better to be explicit::

            sage: f(z=2)
            Gamma(2, w)
        """
        if not args:
            d = None
        elif len(args) == 1 and isinstance(args[0], dict):
            d = args[0]
        else:
            from sage.misc.superseded import deprecation
            deprecation(5930, "Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)")
            d = {}

            vars = _the_element.variables()
            for i, arg in enumerate(args):
                try:
                    d[ vars[i] ] = arg
                except IndexError:
                    raise ValueError("the number of arguments must be less than or equal to %s"%len(vars))

        return _the_element.subs(d, **kwds)

    def subring(self, *args, **kwds):
        r"""
        Create a subring of this symbolic ring.

        INPUT:

        Choose one of the following keywords to create a subring.

        - ``accepting_variables`` (default: ``None``) -- a tuple or other
          iterable of variables. If specified, then a symbolic subring of
          expressions in only these variables is created.

        - ``rejecting_variables`` (default: ``None``) -- a tuple or other
          iterable of variables. If specified, then a symbolic subring of
          expressions in variables distinct to these variables is
          created.

        - ``no_variables`` (default: ``False``) -- a boolean. If set,
          then a symbolic subring of constant expressions (i.e.,
          expressions without a variable) is created.

        OUTPUT:

        A ring.

        EXAMPLES:

        Let us create a couple of symbolic variables first::

            sage: V = var('a, b, r, s, x, y')

        Now we create a symbolic subring only accepting expressions in
        the variables `a` and `b`::

            sage: A = SR.subring(accepting_variables=(a, b)); A
            Symbolic Subring accepting the variables a, b

        An element is
        ::

            sage: A.an_element()
            a

        From our variables in `V` the following are valid in `A`::

            sage: tuple(v for v in V if v in A)
            (a, b)

        Next, we create a symbolic subring rejecting expressions with
        given variables::

            sage: R = SR.subring(rejecting_variables=(r, s)); R
            Symbolic Subring rejecting the variables r, s

        An element is
        ::

            sage: R.an_element()
            some_variable

        From our variables in `V` the following are valid in `R`::

            sage: tuple(v for v in V if v in R)
            (a, b, x, y)

        We have a third kind of subring, namely the subring of
        symbolic constants::

            sage: C = SR.subring(no_variables=True); C
            Symbolic Constants Subring

        Note that this subring can be considered as a special accepting
        subring; one without any variables.

        An element is
        ::

            sage: C.an_element()
            I*pi*e

        None of our variables in `V` is valid in `C`::

            sage: tuple(v for v in V if v in C)
            ()

        .. SEEALSO::

            :doc:`subring`
        """
        if self is not SR:
            raise NotImplementedError('Cannot create subring of %s.' % (self,))
        from .subring import SymbolicSubring
        return SymbolicSubring(*args, **kwds)

    def _fricas_init_(self):
        """
        Return a FriCAS representation of ``self``.

        EXAMPLES::

            sage: fricas(SR)          # indirect doctest, optional - fricas
            Expression(Integer)
        """
        return 'Expression Integer'


SR = SymbolicRing()


cdef unsigned sage_domain_to_ginac_domain(object domain) except? 3474701533:
    """
    TESTS::

        sage: var('x', domain='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo': domain must be one of 'complex', 'real', 'positive' or 'integer'
    """
    # convert the domain argument to something easy to parse
    if domain is RR or domain == 'real':
        return domain_real
    elif domain == 'positive':
        return domain_positive
    elif domain is CC or domain == 'complex':
        return domain_complex
    elif domain is ZZ or domain == 'integer':
        return domain_integer
    else:
        raise ValueError(repr(domain)+": domain must be one of 'complex', 'real', 'positive' or 'integer'")

cdef void send_sage_domain_to_maxima(Expression v, object domain) except +:
    from sage.symbolic.assumptions import assume
    # convert the domain argument to something easy to parse
    if domain is RR or domain == 'real':
        assume(v, 'real')
    elif domain == 'positive':
        assume(v>0)
    elif domain is CC or domain == 'complex':
        assume(v, 'complex')
    elif domain is ZZ or domain == 'integer':
        assume(v, 'integer')
    else:
        raise ValueError(repr(domain)+": domain must be one of 'complex', 'real', 'positive' or 'integer'")

cdef class NumpyToSRMorphism(Morphism):
    r"""
    A morphism from numpy types to the symbolic ring.

    TESTS:

    We check that :trac:`8949` and :trac:`9769` are fixed (see also :trac:`18076`)::

        sage: import numpy
        sage: f(x) = x^2
        sage: f(numpy.int8('2'))
        4
        sage: f(numpy.int32('3'))
        9

    Note that the answer is a Sage integer and not a numpy type::

        sage: a = f(numpy.int8('2')).pyobject()
        sage: type(a)
        <type 'sage.rings.integer.Integer'>

    This behavior also applies to standard functions::

        sage: cos(int('2'))
        cos(2)
        sage: numpy.cos(int('2'))
        -0.4161468365471424
    """
    cdef _intermediate_ring

    def __init__(self, numpy_type):
        """
        A Morphism which constructs Expressions from NumPy floats and
        complexes by converting them to elements of either RDF or CDF.

        INPUT:

        - ``numpy_type`` - a numpy number type

        EXAMPLES::

            sage: import numpy
            sage: from sage.symbolic.ring import NumpyToSRMorphism
            sage: f = NumpyToSRMorphism(numpy.float64)
            sage: f(numpy.float64('2.0'))
            2.0
            sage: _.parent()
            Symbolic Ring

            sage: NumpyToSRMorphism(str)
            Traceback (most recent call last):
            ...
            TypeError: <... 'str'> is not a numpy number type
        """
        Morphism.__init__(self, numpy_type, SR)

        import numpy
        if issubclass(numpy_type, numpy.integer):
            from sage.rings.all import ZZ
            self._intermediate_ring = ZZ
        elif issubclass(numpy_type, numpy.floating):
            from sage.rings.all import RDF
            self._intermediate_ring = RDF
        elif issubclass(numpy_type, numpy.complexfloating):
            from sage.rings.all import CDF
            self._intermediate_ring = CDF
        else:
            raise TypeError("{} is not a numpy number type".format(numpy_type))

    cpdef Element _call_(self, a):
        """
        EXAMPLES:

        This should be called when coercing or converting a NumPy
        float or complex to the Symbolic Ring::

            sage: import numpy
            sage: SR(numpy.int32('1')).pyobject().parent()
            Integer Ring
            sage: SR(numpy.int64('-2')).pyobject().parent()
            Integer Ring

            sage: SR(numpy.float16('1')).pyobject().parent()
            Real Double Field
            sage: SR(numpy.float64('2.0')).pyobject().parent()
            Real Double Field

            sage: SR(numpy.complex64(1jr)).pyobject().parent()
            Complex Double Field
        """
        return new_Expression_from_pyobject(self.codomain(), self._intermediate_ring(a))

cdef class UnderscoreSageMorphism(Morphism):
    def __init__(self, t, R):
        """
        A Morphism which constructs Expressions from an arbitrary Python
        object by calling the :meth:`_sage_` method on the object.

        EXAMPLES::

            sage: import sympy
            sage: from sage.symbolic.ring import UnderscoreSageMorphism
            sage: b = sympy.var('b')
            sage: f = UnderscoreSageMorphism(type(b), SR)
            sage: f(b)
            b
            sage: _.parent()
            Symbolic Ring
        """
        import sage.categories.homset
        from sage.sets.pythonclass import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(t), R))

    cpdef Element _call_(self, a):
        """
        EXAMPLES:

        This should be called when coercing or converting a SymPy
        object to the Symbolic Ring::

            sage: import sympy
            sage: b = sympy.var('b')
            sage: bool(SR(b) == SR(b._sage_()))
            True
        """
        return self.codomain()(a._sage_())


def the_SymbolicRing():
    """
    Return the unique symbolic ring object.

    (This is mainly used for unpickling.)

    EXAMPLES::

        sage: sage.symbolic.ring.the_SymbolicRing()
        Symbolic Ring
        sage: sage.symbolic.ring.the_SymbolicRing() is sage.symbolic.ring.the_SymbolicRing()
        True
        sage: sage.symbolic.ring.the_SymbolicRing() is SR
        True
    """
    return SR

def is_SymbolicExpressionRing(R):
    """
    Returns True if *R* is the symbolic expression ring.

    EXAMPLES::

        sage: from sage.symbolic.ring import is_SymbolicExpressionRing
        sage: is_SymbolicExpressionRing(ZZ)
        False
        sage: is_SymbolicExpressionRing(SR)
        True
    """
    return R is SR

def var(name, **kwds):
    """
    EXAMPLES::

        sage: from sage.symbolic.ring import var
        sage: var("x y z")
        (x, y, z)
        sage: var("x,y,z")
        (x, y, z)
        sage: var("x , y , z")
        (x, y, z)
        sage: var("z")
        z

    TESTS:

    These examples test that variables can only be made from valid
    identifiers.  See :trac:`7496` (and :trac:`9724`) for details::

        sage: var(' ')
        Traceback (most recent call last):
        ...
        ValueError: You need to specify the name of the new variable.
        sage: var('3')
        Traceback (most recent call last):
        ...
        ValueError: The name "3" is not a valid Python identifier.
    """
    return SR.var(name, **kwds)


def is_SymbolicVariable(x):
    """
    Return ``True`` if ``x`` is a variable.

    EXAMPLES::

        sage: from sage.symbolic.ring import is_SymbolicVariable
        sage: is_SymbolicVariable(x)
        True
        sage: is_SymbolicVariable(x+2)
        False

    TESTS::

        sage: ZZ['x']
        Univariate Polynomial Ring in x over Integer Ring
    """
    return is_Expression(x) and is_a_symbol((<Expression>x)._gobj)


def isidentifier(x):
    """
    Return whether ``x`` is a valid identifier.

    INPUT:

    - ``x`` -- a string

    OUTPUT:

    Boolean. Whether the string ``x`` can be used as a variable name.

    This function should return ``False`` for keywords, so we can not
    just use the ``isidentifier`` method of strings,
    because, for example, it returns ``True`` for "def" and for "None".

    EXAMPLES::

        sage: from sage.symbolic.ring import isidentifier
        sage: isidentifier('x')
        True
        sage: isidentifier(' x')   # can't start with space
        False
        sage: isidentifier('ceci_n_est_pas_une_pipe')
        True
        sage: isidentifier('1 + x')
        False
        sage: isidentifier('2good')
        False
        sage: isidentifier('good2')
        True
        sage: isidentifier('lambda s:s+1')
        False
        sage: isidentifier('None')
        False
        sage: isidentifier('lambda')
        False
        sage: isidentifier('def')
        False
    """
    if x in KEYWORDS:
        return False
    return x.isidentifier()
