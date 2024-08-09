"""
The symbolic ring
"""

# ****************************************************************************
#       Copyright (C) 2008      William Stein <wstein@gmail.com>
#       Copyright (C) 2008-2013 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2009      Mike Hansen
#       Copyright (C) 2011      Karl-Dieter Crisman
#       Copyright (C) 2011-2012 Volker Braun
#       Copyright (C) 2013-2019 Frédéric Chapoton
#       Copyright (C) 2014-2020 Marc Mezzarobba
#       Copyright (C) 2015      Bruno Grenet
#       Copyright (C) 2015-2016 Daniel Krenn
#       Copyright (C) 2015-2016 Jeroen Demeyer
#       Copyright (C) 2015-2017 Vincent Delecroix
#       Copyright (C) 2015-2018 Ralf Stephan
#       Copyright (C) 2016      Julian Rüth
#       Copyright (C) 2017      Marcelo Forets
#       Copyright (C) 2018      Martin Rubey
#       Copyright (C) 2019      E. Madison Bray
#       Copyright (C) 2019      Markus Wageringel
#       Copyright (C) 2021      Marius Gerbershagen
#       Copyright (C) 2021      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer cimport Integer

import sage.rings.abc

from sage.symbolic.expression cimport (
    _latex_Expression,
    _repr_Expression,
    new_Expression,
    new_Expression_from_pyobject,
    new_Expression_wild,
    new_Expression_symbol,
)

from sage.categories.commutative_rings import CommutativeRings
from sage.structure.element cimport Element, Expression
from sage.structure.parent cimport Parent
from sage.categories.morphism cimport Morphism
from sage.structure.coerce cimport is_numpy_type

import sage.rings.abc
from sage.rings.integer_ring import ZZ

# is_SymbolicVariable used to be defined here; re-export it
from sage.symbolic.expression import _is_SymbolicVariable as is_SymbolicVariable

import keyword
import operator

# Do not allow any of these keywords as identifiers for symbolic variables
KEYWORDS = set(keyword.kwlist).union(['exec', 'print', 'None', 'True',
                                      'False', 'nonlocal'])


cdef class SymbolicRing(sage.rings.abc.SymbolicRing):
    """
    Symbolic Ring, parent object for all symbolic expressions.
    """
    def __init__(self, base_ring=None):
        """
        Initialize the Symbolic Ring.

        This is a commutative ring of symbolic expressions and functions.

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
        Parent.__init__(self, base_ring, category=CommutativeRings())
        self._populate_coercion_lists_(convert_method_name='_symbolic_')
        self.symbols = {}

    def __reduce__(self):
        """
        EXAMPLES::

           sage: loads(dumps(SR)) == SR           # indirect doctest
           True
        """
        return the_SymbolicRing, tuple()

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: repr(SR)
            'Symbolic Ring'
        """
        return "Symbolic Ring"

    def _latex_(self):
        r"""
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

        Check if arithmetic with bools works (see :issue:`9560`)::

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
            from sage.rings.fraction_field import FractionField_generic
            from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
            from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
            from sage.rings.polynomial.laurent_polynomial_ring_base import LaurentPolynomialRing_generic
            from sage.rings.infinity import InfinityRing, UnsignedInfinityRing
            from sage.rings.real_lazy import RLF, CLF
            from sage.rings.finite_rings.finite_field_base import FiniteField

            from sage.symbolic.subring import GenericSymbolicSubring

            if R._is_numerical():
                # Almost anything with a coercion into any precision of CC
                return R not in (RLF, CLF)
            elif isinstance(R, PolynomialRing_general) or isinstance(R, MPolynomialRing_base) or isinstance(R, FractionField_generic) or isinstance(R, LaurentPolynomialRing_generic):
                base = R.base_ring()
                return base is not self and self.has_coerce_map_from(base)
            elif (R is InfinityRing or R is UnsignedInfinityRing
                  or isinstance(R, (sage.rings.abc.RealIntervalField,
                                    sage.rings.abc.ComplexIntervalField,
                                    sage.rings.abc.RealBallField,
                                    sage.rings.abc.ComplexBallField,
                                    sage.rings.abc.IntegerModRing,
                                    FiniteField))):
                return True
            elif isinstance(R, GenericSymbolicSubring):
                return True

    def _element_constructor_(self, x):
        r"""
        Convert `x` into the symbolic expression ring SR.

        EXAMPLES::

            sage: a = SR(-3/4); a
            -3/4
            sage: type(a)
            <class 'sage.symbolic.expression.Expression'>
            sage: a.parent()
            Symbolic Ring
            sage: K.<a> = QuadraticField(-3)
            sage: a + sin(x)
            I*sqrt(3) + sin(x)
            sage: x = var('x'); y0,y1 = PolynomialRing(ZZ,2,'y').gens()
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

        Any proper subset of the complex numbers::

            sage: SR(NN)
            Non negative integer semiring
            sage: SR(ZZ)
            Integer Ring
            sage: SR(Set([1/2, 2/3, 3/4]))
            {3/4, 2/3, 1/2}
            sage: SR(RealSet(0, 1))
            (0, 1)

        TESTS::

            sage: SR.coerce(int(5))
            5
            sage: SR.coerce(5)
            5
            sage: SR.coerce(float(5))
            5.0
            sage: SR.coerce(5.0)
            5.00000000000000

        An interval arithmetic number::

            sage: SR.coerce(RIF(pi))
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

        Check that :issue:`22068` is fixed::

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

        Check that :issue:`24072` is solved::

            sage: x = polygen(GF(3))
            sage: a = SR.var('a')
            sage: (2*x + 1) * a
            Traceback (most recent call last):
            ...
            TypeError: positive characteristic not allowed in symbolic computations

        Check support for unicode characters (:issue:`29280`)::

            sage: SR('λ + 2λ')
            3*λ
            sage: SR('μ') is var('μ')
            True
            sage: SR('λ + * 1')
            Traceback (most recent call last):
            ...
            TypeError: Malformed expression: λ + * !!!  1
        """
        return new_Expression(self, x)

    def _force_pyobject(self, x, bint force=False, bint recursive=True):
        r"""
        Wrap the given Python object in a symbolic expression even if it
        cannot be coerced to the Symbolic Ring.

        INPUT:

        - ``x`` -- a Python object

        - ``force`` -- boolean (default: ``False``); if ``True``, the Python object
          is taken as is without attempting coercion or list traversal

        - ``recursive`` -- boolean (default: ``True``); disables recursive
          traversal of lists

        EXAMPLES::

            sage: t = SR._force_pyobject(QQ); t
            Rational Field
            sage: type(t)
            <class 'sage.symbolic.expression.Expression'>

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
        return new_Expression_from_pyobject(self, x, force, recursive)

    def wild(self, unsigned int n=0):
        r"""
        Return the n-th wild-card for pattern matching and substitution.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: n-th wildcard expression

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

        Check that :issue:`15047` is fixed::

            sage: latex(SR.wild(0))
            \$0

        Check that :issue:`21455` is fixed::

            sage: coth(SR.wild(0))
            coth($0)
        """
        return new_Expression_wild(self, n)

    def __contains__(self, x):
        r"""
        ``True`` if there is an element of the symbolic ring that is equal to x
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

        OUTPUT: a Sage integer

        EXAMPLES::

            sage: c = SR.characteristic(); c
            0
            sage: type(c)
            <class 'sage.rings.integer.Integer'>
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

    def is_field(self, proof=True):
        """
        Return ``True``, since the symbolic expression ring is (for the most
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

        TESTS:

        Test that :issue:`32404` is fixed::

            sage: SR0 = SR.subring(no_variables=True)
            sage: SR0.I().parent()
            Symbolic Constants Subring
        """
        from sage.symbolic.constants import I
        return self(I)

    def symbol(self, name=None, latex_name=None, domain=None):
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
        return new_Expression_symbol(self, name, latex_name, domain)

    def temp_var(self, n=None, domain=None):
        """
        Return one or multiple new unique symbolic variables as an element
        of the symbolic ring. Use this instead of SR.var() if there is a
        possibility of name clashes occuring. Call SR.cleanup_var() once
        the variables are no longer needed or use a `with SR.temp_var()
        as ...` construct.

        INPUT:

        - ``n`` -- (optional) positive integer; number of symbolic variables

        - ``domain`` -- (optional) specify the domain of the variable(s);

        EXAMPLES:

        Simple definition of a functional derivative::

            sage: def functional_derivative(expr,f,x):
            ....:     with SR.temp_var() as a:
            ....:         return expr.subs({f(x):a}).diff(a).subs({a:f(x)})
            sage: f = function('f')
            sage: a = var('a')
            sage: functional_derivative(f(a)^2+a,f,a)
            2*f(a)

        Contrast this to a similar implementation using SR.var(),
        which gives a wrong result in our example::

            sage: def functional_derivative(expr,f,x):
            ....:     a = SR.var('a')
            ....:     return expr.subs({f(x):a}).diff(a).subs({a:f(x)})
            sage: f = function('f')
            sage: a = var('a')
            sage: functional_derivative(f(a)^2+a,f,a)
            2*f(a) + 1

        TESTS:

            sage: x = SR.temp_var()
            sage: y = SR.temp_var()
            sage: bool(x == x)
            True
            sage: bool(x == y)
            False
            sage: bool(x.parent()(x._maxima_()) == x)
            True
        """
        if n is None:
            return self.symbol(None, domain=domain)
        return TemporaryVariables([self.temp_var(domain=domain)
                                   for i in range(n)])

    def cleanup_var(self, symbol):
        """
        Cleans up a variable, removing assumptions about the
        variable and allowing for it to be garbage collected

        INPUT:

        - ``symbol`` -- a variable or a list of variables

        TESTS:

            sage: from sage.symbolic.assumptions import assumptions
            sage: symbols_copy = SR.symbols.copy()
            sage: assumptions_copy = assumptions().copy()
            sage: x = SR.temp_var(domain='real')
            sage: SR.cleanup_var(x)
            sage: symbols_copy == SR.symbols
            True
            sage: assumptions_copy == assumptions()
            True
        """
        from sage.symbolic.assumptions import assumptions
        if isinstance(symbol, (list, tuple)):
            for s in symbol:
                self.cleanup_var(s)
        else:
            try:
                name = self._repr_element_(symbol)
                del self.symbols[name]
            except KeyError:
                pass
            for asm in assumptions():
                if asm.has(symbol):
                    asm.forget()

    def var(self, name, latex_name=None, n=None, domain=None):
        r"""
        Return a symbolic variable as an element of the symbolic ring.

        INPUT:

        - ``name`` -- string or list of strings with the name(s) of the symbolic variable(s)

        - ``latex_name`` -- (optional) string used when printing in latex mode, if not specified use ``'name'``

        - ``n`` -- (optional) positive integer; number of symbolic variables, indexed from `0` to `n-1`

        - ``domain`` -- (optional) specify the domain of the variable(s); it is None
          by default, and possible options are (non-exhaustive list, see note below):
          ``'real'``, ``'complex'``, ``'positive'``, ``'integer'`` and ``'noninteger'``

        OUTPUT: symbolic expression or tuple of symbolic expressions

        .. SEEALSO::

            This function does not inject the variable(s) into the global namespace.
            For that purpose see :meth:`var()<sage.calculus.var.var>`.

        .. NOTE::

            For a comprehensive list of acceptable features type ``'maxima('features')'``,
            and see also the documentation of :ref:`sage.symbolic.assumptions`.

        EXAMPLES:

        Create a variable `zz`::

            sage: zz = SR.var('zz'); zz
            zz

        The return type is a symbolic expression::

            sage: type(zz)
            <class 'sage.symbolic.expression.Expression'>

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

        Check that :issue:`17206` is fixed::

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

        Check that :issue:`28353` is fixed: Constructions that suggest multiple
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
        if isinstance(name, Expression):
            return name
        if not isinstance(name, (str, list, tuple)):
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
            # warn on bad symbol names, but only once
            # symbol... names are temporary variables created with
            #   SR.temp_var
            # _symbol... names are used in the conversion of
            #   derivatives of symbolic functions to maxima and other
            #   external libraries
            if self.symbols.get(s) is None and ((s.startswith('symbol') and s[6:].isdigit()) or (s.startswith('_symbol') and s[7:].isdigit())):
                import warnings
                warnings.warn(f'The name "{name}" may clash with names used internally in sagemath. It is recommended to choose a different name for your variable.')

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

    def _repr_element_(self, x):
        """
        Return the string representation of the expression ``x``.

        This is used so that subclasses of :class:`SymbolicRing` (such as a
        :class:`~sage.symbolic.callable.CallableSymbolicExpressionRing`)
        can provide their own implementations of how to print expressions.

        EXAMPLES::

            sage: SR._repr_element_(x+2)
            'x + 2'
        """
        return _repr_Expression(x)

    def _latex_element_(self, x):
        r"""
        Return the standard LaTeX version of the expression ``x``.

        EXAMPLES::

            sage: latex(sin(x+2))
            \sin\left(x + 2\right)
            sage: latex(var('theta') + 2)
            \theta + 2
        """
        return _latex_Expression(x)

    def _call_element_(self, _the_element, *args, **kwds):
        """
        EXAMPLES::

            sage: x, y = var('x,y')
            sage: f = x+y
            sage: f.variables()
            (x, y)
            sage: f()
            x + y
            sage: f(3)
            Traceback (most recent call last):
            ...
            TypeError: Substitution using function-call syntax and unnamed arguments
            has been removed. You can use named arguments instead, like EXPR(x=..., y=...)
            sage: f(x=3)
            y + 3
            sage: f(3, 4)
            Traceback (most recent call last):
            ...
            TypeError: Substitution using function-call syntax and unnamed arguments
            has been removed. You can use named arguments instead, like EXPR(x=..., y=...)
            sage: f(x=3, y=4)
            7
            sage: f(2, 3, 4)
            Traceback (most recent call last):
            ...
            TypeError: Substitution using function-call syntax and unnamed arguments
            has been removed. You can use named arguments instead, like EXPR(x=..., y=...)
            sage: f(x=2, y=3, z=4)
            5

        ::

            sage: f({x: 3})
            y + 3
            sage: f({x: 3, y: 4})
            7
            sage: f(x=3)
            y + 3
            sage: f(x=3, y=4)
            7

        ::

            sage: a = (2^(8/9))
            sage: a(4)
            Traceback (most recent call last):
            ...
            TypeError: Substitution using function-call syntax and unnamed arguments
            has been removed. You can use named arguments instead, like EXPR(x=..., y=...)

        Note that the application of arguments to a function defined using `function`
        creates an ordinary expression, not a callable symbolic expression.  Hence,
        calling this expression using function-call syntax and unnamed arguments
        leads to an error::

            sage: f = function('Gamma')(var('z'), var('w')); f
            Gamma(z, w)
            sage: f(2)
            Traceback (most recent call last):
            ...
            TypeError: Substitution using function-call syntax and unnamed arguments
            has been removed. You can use named arguments instead, like EXPR(x=..., y=...)
            sage: f(2,5)
            Traceback (most recent call last):
            ...
            TypeError: Substitution using function-call syntax and unnamed arguments
            has been removed. You can use named arguments instead, like EXPR(x=..., y=...)

        Thus, it is better to be explicit::

            sage: f(z=2)
            Gamma(2, w)
        """
        if not args:
            d = None
        elif len(args) == 1 and isinstance(args[0], dict):
            d = args[0]
        else:
            raise TypeError("Substitution using function-call syntax "
                            "and unnamed arguments has been removed. You "
                            "can use named arguments instead, like "
                            "EXPR(x=..., y=...)")
        return _the_element.subs(d, **kwds)

    def subring(self, *args, **kwds):
        r"""
        Create a subring of this symbolic ring.

        INPUT:

        Choose one of the following keywords to create a subring.

        - ``accepting_variables`` -- (default: ``None``) a tuple or other
          iterable of variables. If specified, then a symbolic subring of
          expressions in only these variables is created.

        - ``rejecting_variables`` -- (default: ``None``) a tuple or other
          iterable of variables. If specified, then a symbolic subring of
          expressions in variables distinct to these variables is
          created.

        - ``no_variables`` -- boolean (default: ``False``); if set,
          then a symbolic subring of constant expressions (i.e.,
          expressions without a variable) is created.

        OUTPUT: a ring

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
            raise NotImplementedError('cannot create subring of %s' % (self,))
        from sage.symbolic.subring import SymbolicSubring
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


cdef class NumpyToSRMorphism(Morphism):
    r"""
    A morphism from numpy types to the symbolic ring.

    TESTS:

    We check that :issue:`8949` and :issue:`9769` are fixed (see also :issue:`18076`)::

        sage: import numpy                                                              # needs numpy
        sage: if int(numpy.version.short_version[0]) > 1:                               # needs numpy
        ....:     numpy.set_printoptions(legacy="1.25")                                 # needs numpy
        sage: f(x) = x^2
        sage: f(numpy.int8('2'))                                                        # needs numpy
        4
        sage: f(numpy.int32('3'))                                                       # needs numpy
        9

    Note that the answer is a Sage integer and not a numpy type::

        sage: a = f(numpy.int8('2')).pyobject()                                         # needs numpy
        sage: type(a)                                                                   # needs numpy
        <class 'sage.rings.integer.Integer'>

    This behavior also applies to standard functions::

        sage: cos(int('2'))
        cos(2)
        sage: numpy.cos(int('2'))                                                       # needs numpy
        -0.4161468365471424
    """
    cdef _intermediate_ring

    def __init__(self, numpy_type):
        """
        A Morphism which constructs Expressions from NumPy floats and
        complexes by converting them to elements of either RDF or CDF.

        INPUT:

        - ``numpy_type`` -- a numpy number type

        EXAMPLES::

            sage: # needs numpy
            sage: import numpy
            sage: from sage.symbolic.ring import NumpyToSRMorphism
            sage: f = NumpyToSRMorphism(numpy.float64)
            sage: f(numpy.float64('2.0'))
            2.0
            sage: _.parent()
            Symbolic Ring

            sage: NumpyToSRMorphism(str)                                                # needs numpy
            Traceback (most recent call last):
            ...
            TypeError: <... 'str'> is not a numpy number type
        """
        Morphism.__init__(self, numpy_type, SR)

        import numpy
        if issubclass(numpy_type, numpy.integer):
            from sage.rings.integer_ring import ZZ
            self._intermediate_ring = ZZ
        elif issubclass(numpy_type, numpy.floating):
            from sage.rings.real_double import RDF
            self._intermediate_ring = RDF
        elif issubclass(numpy_type, numpy.complexfloating):
            from sage.rings.complex_double import CDF
            self._intermediate_ring = CDF
        else:
            raise TypeError("{} is not a numpy number type".format(numpy_type))

    cpdef Element _call_(self, a):
        """
        EXAMPLES:

        This should be called when coercing or converting a NumPy
        float or complex to the Symbolic Ring::

            sage: # needs numpy
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
        return new_Expression_from_pyobject(self.codomain(), self._intermediate_ring(a), True)


cdef class UnderscoreSageMorphism(Morphism):
    def __init__(self, t, R):
        """
        A Morphism which constructs Expressions from an arbitrary Python
        object by calling the :meth:`_sage_` method on the object.

        EXAMPLES::

            sage: # needs sympy
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
        from sage.interfaces.sympy import sympy_init
        sympy_init()

    cpdef Element _call_(self, a):
        """
        EXAMPLES:

        This should be called when coercing or converting a SymPy
        object to the Symbolic Ring::

            sage: import sympy                                                          # needs sympy
            sage: b = sympy.var('b')                                                    # needs sympy
            sage: bool(SR(b) == SR(b._sage_()))                                         # needs sympy
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
    identifiers.  See :issue:`7496` (and :issue:`9724`) for details::

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


def isidentifier(x):
    """
    Return whether ``x`` is a valid identifier.

    INPUT:

    - ``x`` -- string

    OUTPUT: boolean; whether the string ``x`` can be used as a variable name

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


class TemporaryVariables(tuple):
    """
    Instances of this class can be used with Python `with` to
    automatically clean up after themselves.
    """
    def __enter__(self):
        return self

    def __exit__(self, *args):
        """
        TESTS::

            sage: symbols_copy = SR.symbols.copy()
            sage: with SR.temp_var(n=2) as temp_vars: pass
            sage: symbols_copy == SR.symbols
            True
        """
        SR.cleanup_var(self)
        return False
