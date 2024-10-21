# sage_setup: distribution = sagemath-objects
r"""
The coercion model

The coercion model manages how elements of one parent get related to elements
of another. For example, the integer 2 can canonically be viewed as an element
of the rational numbers. (The parent of a non-element is its Python type.)

::

    sage: ZZ(2).parent()
    Integer Ring
    sage: QQ(2).parent()
    Rational Field

The most prominent role of the coercion model is to make sense of binary
operations between elements that have distinct parents. It does this by
finding a parent where both elements make sense, and doing the operation
there. For example::

    sage: a = 1/2; a.parent()
    Rational Field
    sage: b = ZZ['x'].gen(); b.parent()
    Univariate Polynomial Ring in x over Integer Ring
    sage: a + b
    x + 1/2
    sage: (a + b).parent()
    Univariate Polynomial Ring in x over Rational Field

If there is a coercion (see below) from one of the parents to the other,
the operation is always performed in the codomain of that coercion. Otherwise
a reasonable attempt to create a new parent with coercion maps from both
original parents is made. The results of these discoveries are cached.
On failure, a :exc:`TypeError` is always raised.

Some arithmetic operations (such as multiplication) can indicate an action
rather than arithmetic in a common parent. For example::

    sage: E = EllipticCurve('37a')                                                      # needs sage.schemes
    sage: P = E(0,0)                                                                    # needs sage.schemes
    sage: 5*P                                                                           # needs sage.schemes
    (1/4 : -5/8 : 1)

where there is action of `\ZZ` on the points of `E` given by the additive
group law. Parents can specify how they act on or are acted upon by other
parents.

There are two kinds of ways to get from one parent to another, coercions and
conversions.

Coercions are canonical (possibly modulo a finite number of
deterministic choices) morphisms, and the set of all coercions between
all parents forms a commuting diagram (modulo possibly rounding
issues). `\ZZ \rightarrow \QQ` is an example of a
coercion. These are invoked implicitly by the coercion model.

Conversions try to construct an element out of their input if at all possible.
Examples include sections of coercions, creating an element from a string or
list, etc. and may fail on some inputs of a given type while succeeding on
others (i.e. they may not be defined on the whole domain). Conversions are
always explicitly invoked, and never used by the coercion model to resolve
binary operations.

For more information on how to specify coercions, conversions, and actions,
see the documentation for :class:`Parent`.
"""

# ****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport (PyTypeObject, PyObject_CallObject,
                             PyObject_RichCompare, Py_TYPE,
                             Py_EQ, Py_NE, Py_LT, Py_LE, Py_GT)
from cpython.weakref cimport PyWeakref_GET_OBJECT, PyWeakref_NewRef
from libc.string cimport strncmp
cimport gmpy2

cdef mul, truediv
from operator import mul, truediv

from sage.structure.richcmp cimport rich_to_bool, revop
from sage.structure.sage_object cimport SageObject
from sage.structure.parent cimport Parent_richcmp_element_without_coercion
from sage.structure.element cimport bin_op_exception, parent, Element
from sage.structure.coerce_exceptions import CoercionException
from sage.rings.integer_fake cimport is_Integer
from sage.categories.map cimport Map
from sage.categories.morphism import IdentityMorphism
from sage.categories.action cimport Action, PrecomposedAction
from sage.sets.pythonclass cimport Set_PythonType

import traceback

from fractions import Fraction
cdef type FractionType = <type>Fraction

cpdef py_scalar_parent(py_type):
    """
    Return the Sage equivalent of the given python type, if one exists.
    If there is no equivalent, return ``None``.

    EXAMPLES::

        sage: from sage.structure.coerce import py_scalar_parent
        sage: py_scalar_parent(int)
        Integer Ring
        sage: py_scalar_parent(float)
        Real Double Field
        sage: py_scalar_parent(complex)                                                 # needs sage.rings.complex_double
        Complex Double Field
        sage: py_scalar_parent(bool)
        Integer Ring
        sage: py_scalar_parent(dict),
        (None,)

        sage: import fractions
        sage: py_scalar_parent(fractions.Fraction)
        Rational Field

        sage: # needs numpy
        sage: import numpy
        sage: py_scalar_parent(numpy.int16)
        Integer Ring
        sage: py_scalar_parent(numpy.int32)
        Integer Ring
        sage: py_scalar_parent(numpy.uint64)
        Integer Ring
        sage: py_scalar_parent(numpy.double)
        Real Double Field

        sage: import gmpy2
        sage: py_scalar_parent(gmpy2.mpz)
        Integer Ring
        sage: py_scalar_parent(gmpy2.mpq)
        Rational Field
        sage: py_scalar_parent(gmpy2.mpfr)
        Real Double Field
        sage: py_scalar_parent(gmpy2.mpc)                                               # needs sage.rings.complex_double
        Complex Double Field
    """
    if issubclass(py_type, int):
        import sage.rings.integer_ring
        return sage.rings.integer_ring.ZZ
    if py_type is FractionType:
        import sage.rings.rational_field
        return sage.rings.rational_field.QQ
    elif issubclass(py_type, float):
        import sage.rings.real_double
        return sage.rings.real_double.RDF
    elif issubclass(py_type, complex):
        import sage.rings.complex_double
        return sage.rings.complex_double.CDF
    elif is_numpy_type(py_type):
        import numpy
        if issubclass(py_type, numpy.integer):
            import sage.rings.integer_ring
            return sage.rings.integer_ring.ZZ
        elif issubclass(py_type, numpy.floating):
            import sage.rings.real_double
            return sage.rings.real_double.RDF
        elif issubclass(py_type, numpy.complexfloating):
            import sage.rings.complex_double
            return sage.rings.complex_double.CDF
        else:
            return None
    elif issubclass(py_type, gmpy2.mpz):
        import sage.rings.integer_ring
        return sage.rings.integer_ring.ZZ
    elif issubclass(py_type, gmpy2.mpq):
        import sage.rings.rational_field
        return sage.rings.rational_field.QQ
    elif issubclass(py_type, gmpy2.mpfr):
        import sage.rings.real_double
        return sage.rings.real_double.RDF
    elif issubclass(py_type, gmpy2.mpc):
        import sage.rings.complex_double
        return sage.rings.complex_double.CDF
    else:
        return None

cpdef py_scalar_to_element(x):
    """
    Convert ``x`` to a Sage :class:`~sage.structure.element.Element` if possible.

    If ``x`` was already an :class:`~sage.structure.element.Element` or if there is no obvious
    conversion possible, just return ``x`` itself.

    EXAMPLES::

        sage: from sage.structure.coerce import py_scalar_to_element
        sage: x = py_scalar_to_element(42)
        sage: x, parent(x)
        (42, Integer Ring)
        sage: x = py_scalar_to_element(int(42))
        sage: x, parent(x)
        (42, Integer Ring)
        sage: x = py_scalar_to_element(float(42))
        sage: x, parent(x)
        (42.0, Real Double Field)
        sage: x = py_scalar_to_element(complex(42))                                     # needs sage.rings.complex_double
        sage: x, parent(x)                                                              # needs sage.rings.complex_double
        (42.0, Complex Double Field)
        sage: py_scalar_to_element('hello')
        'hello'

        sage: from fractions import Fraction
        sage: f = Fraction(int(2^100), int(3^100))
        sage: py_scalar_to_element(f)
        1267650600228229401496703205376/515377520732011331036461129765621272702107522001

    Note that bools are converted to 0 or 1::

        sage: py_scalar_to_element(False), py_scalar_to_element(True)
        (0, 1)

    Test gmpy2's types::

        sage: import gmpy2
        sage: x = py_scalar_to_element(gmpy2.mpz(42))
        sage: x, parent(x)
        (42, Integer Ring)
        sage: x = py_scalar_to_element(gmpy2.mpq('3/4'))
        sage: x, parent(x)
        (3/4, Rational Field)
        sage: x = py_scalar_to_element(gmpy2.mpfr(42.57))
        sage: x, parent(x)
        (42.57, Real Double Field)
        sage: x = py_scalar_to_element(gmpy2.mpc(int(42), int(42)))                     # needs sage.rings.complex_double
        sage: x, parent(x)                                                              # needs sage.rings.complex_double
        (42.0 + 42.0*I, Complex Double Field)

    Test compatibility with :func:`py_scalar_parent`::

        sage: from sage.structure.coerce import py_scalar_parent
        sage: elt = [True, int(42), float(42), complex(42)]
        sage: for x in elt:                                                             # needs sage.rings.complex_double
        ....:     assert py_scalar_parent(type(x)) == py_scalar_to_element(x).parent()

        sage: import numpy                                                              # needs numpy
        sage: elt = [numpy.int8('-12'),  numpy.uint8('143'),                            # needs numpy
        ....:        numpy.int16('-33'), numpy.uint16('122'),
        ....:        numpy.int32('-19'), numpy.uint32('44'),
        ....:        numpy.int64('-3'),  numpy.uint64('552'),
        ....:        numpy.float16('-1.23'), numpy.float32('-2.22'),
        ....:        numpy.float64('-3.412'), numpy.complex64(1.2+I),
        ....:        numpy.complex128(-2+I)]
        sage: for x in elt:                                                             # needs numpy
        ....:     assert py_scalar_parent(type(x)) == py_scalar_to_element(x).parent()

        sage: elt = [gmpy2.mpz(42), gmpy2.mpq('3/4'),
        ....:        gmpy2.mpfr(42.57), gmpy2.mpc(int(42), int(42))]
        sage: for x in elt:                                                             # needs sage.rings.complex_double
        ....:     assert py_scalar_parent(type(x)) == py_scalar_to_element(x).parent()
    """
    if isinstance(x, Element):
        return x
    elif isinstance(x, int):
        from sage.rings.integer import Integer
        return Integer(x)
    elif type(x) is FractionType:
        from sage.rings.rational import Rational
        return Rational(x)
    elif isinstance(x, float):
        from sage.rings.real_double import RDF
        return RDF(x)
    elif isinstance(x, complex):
        from sage.rings.complex_double import CDF
        return CDF(x)
    elif is_numpy_type(type(x)):
        import numpy
        if isinstance(x, numpy.integer):
            from sage.rings.integer import Integer
            return Integer(x)
        elif isinstance(x, numpy.floating):
            from sage.rings.real_double import RDF
            return RDF(x)
        elif isinstance(x, numpy.complexfloating):
            from sage.rings.complex_double import CDF
            return CDF(x)
        else:
            return x
    elif type(x) is gmpy2.mpz:
            from sage.rings.integer import Integer
            return Integer(x)
    elif type(x) is gmpy2.mpq:
        from sage.rings.rational import Rational
        return Rational(x)
    elif type(x) is gmpy2.mpfr:
        from sage.rings.real_double import RDF
        return RDF(x)
    elif type(x) is gmpy2.mpc:
        from sage.rings.complex_double import CDF
        return CDF(x)
    else:
        return x


cpdef bint parent_is_integers(P) except -1:
    """
    Check whether the type or parent represents the ring of integers.

    EXAMPLES::

        sage: from sage.structure.coerce import parent_is_integers
        sage: parent_is_integers(int)
        True
        sage: parent_is_integers(float)
        False
        sage: parent_is_integers(bool)
        True
        sage: parent_is_integers(dict)
        False

        sage: import numpy                                                              # needs numpy
        sage: parent_is_integers(numpy.int16)                                           # needs numpy
        True
        sage: parent_is_integers(numpy.uint64)                                          # needs numpy
        True
        sage: parent_is_integers(float)
        False

        sage: import gmpy2
        sage: parent_is_integers(gmpy2.mpz)
        True
        sage: parent_is_integers(gmpy2.mpq)
        False

    Ensure (:issue:`27893`) is fixed::

        sage: K.<f> = QQ[]
        sage: gmpy2.mpz(2) * f
        2*f
    """
    if isinstance(P, type):
        if issubclass(P, int):
            return True
        elif is_numpy_type(P):
            from numpy import integer
            return issubclass(P, integer)
        elif issubclass(P, gmpy2.mpz):
            return True
        else:
            return False
    else:
        from sage.rings.integer_ring import ZZ
        return P is ZZ


def parent_is_numerical(P):
    r"""
    Test if elements of the parent or type ``P`` can be numerically evaluated
    as complex numbers (in a canonical way).

    EXAMPLES::

        sage: from sage.structure.coerce import parent_is_numerical
        sage: import gmpy2
        sage: [parent_is_numerical(R) for R in [QQ, int, complex, gmpy2.mpc]]           # needs sage.rings.complex_double
        [True, True, True, True]
        sage: [parent_is_numerical(R) for R in [RR, CC]]                                # needs sage.rings.real_mpfr
        [True, True]
        sage: parent_is_numerical(QuadraticField(-1))                                   # needs sage.rings.number_field
        True
        sage: import numpy; parent_is_numerical(numpy.complexfloating)                  # needs numpy
        True
        sage: parent_is_numerical(SR)                                                   # needs sage.symbolic
        False
        sage: [parent_is_numerical(R) for R in [QQ['x'], QQ[['x']], str]]
        [False, False, False]
        sage: [parent_is_numerical(R) for R in [RIF, RBF, CIF, CBF]]                    # needs sage.libs.flint
        [False, False, False, False]
    """
    if not isinstance(P, Parent):
        P = py_scalar_parent(P)
        if P is None:
            return False
    return P._is_numerical()


def parent_is_real_numerical(P):
    r"""
    Test if elements of the parent or type ``P`` can be numerically evaluated
    as real numbers (in a canonical way).

    EXAMPLES::

        sage: from sage.structure.coerce import parent_is_real_numerical
        sage: import gmpy2
        sage: [parent_is_real_numerical(R) for R in [RR, QQ, ZZ, RLF, int, float, gmpy2.mpq]]
        [True, True, True, True, True, True, True]
        sage: parent_is_real_numerical(QuadraticField(2))                               # needs sage.rings.number_field
        True
        sage: import numpy; parent_is_real_numerical(numpy.integer)                     # needs numpy
        True
        sage: parent_is_real_numerical(QuadraticField(-1))                              # needs sage.rings.number_field
        False
        sage: [parent_is_real_numerical(R)                                              # needs numpy
        ....:  for R in [CC, complex, gmpy2.mpc, numpy.complexfloating]]
        [False, False, False, False]
        sage: [parent_is_real_numerical(R) for R in [QQ['x'], QQ[['x']], str]]
        [False, False, False]
        sage: parent_is_real_numerical(SR)                                              # needs sage.symbolic
        False
        sage: [parent_is_real_numerical(R) for R in [RIF, RBF, CIF, CBF]]               # needs sage.libs.flint
        [False, False, False, False]
    """
    if not isinstance(P, Parent):
        P = py_scalar_parent(P)
        if P is None:
            return False
    return P._is_real_numerical()


cpdef bint is_numpy_type(t) noexcept:
    """
    Return ``True`` if and only if `t` is a type whose name starts
    with ``numpy.``

    EXAMPLES::

        sage: from sage.structure.coerce import is_numpy_type

        sage: # needs numpy
        sage: import numpy
        sage: is_numpy_type(numpy.int16)
        True
        sage: is_numpy_type(numpy.floating)
        True
        sage: is_numpy_type(numpy.ndarray)
        True
        sage: is_numpy_type(numpy.matrix)
        True

        sage: is_numpy_type(int)
        False
        sage: is_numpy_type(Integer)
        False
        sage: is_numpy_type(Sudoku)                                                     # needs sage.combinat
        False
        sage: is_numpy_type(None)
        False

    TESTS:

    This used to crash Sage (:issue:`20715`)::

        sage: is_numpy_type(object)
        False
        sage: 1 + object()
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for +: 'Integer Ring' and
        '<class 'object'>'
    """
    if not isinstance(t, type):
        return False
    cdef PyTypeObject* T = <PyTypeObject*>t
    if strncmp(T.tp_name, "numpy.", 6) == 0:
        return True
    # Check base type. This is needed to detect numpy.matrix.
    if T.tp_base != NULL and strncmp(T.tp_base.tp_name, "numpy.", 6) == 0:
        return True
    return False

cpdef bint is_mpmath_type(t) noexcept:
    r"""
    Check whether the type ``t`` is a type whose name starts with
    ``sage.libs.mpmath.``.

    EXAMPLES::

        sage: # needs mpmath
        sage: from sage.structure.coerce import is_mpmath_type
        sage: is_mpmath_type(int)
        False
        sage: import sage.libs.mpmath.all
        sage: is_mpmath_type(sage.libs.mpmath.all.mpc(2))
        False
        sage: is_mpmath_type(type(sage.libs.mpmath.all.mpc(2)))
        True
        sage: is_mpmath_type(type(sage.libs.mpmath.all.mpf(2)))
        True
    """
    return isinstance(t, type) and \
           strncmp((<PyTypeObject*>t).tp_name, "sage.libs.mpmath.", 17) == 0


cdef class CoercionModel:
    """
    See also :mod:`sage.categories.pushout`

    EXAMPLES::

        sage: f = ZZ['t', 'x'].0 + QQ['x'].0 + CyclotomicField(13).gen(); f             # needs sage.rings.number_field
        t + x + zeta13
        sage: f.parent()                                                                # needs sage.rings.number_field
        Multivariate Polynomial Ring in t, x
         over Cyclotomic Field of order 13 and degree 12
        sage: ZZ['x','y'].0 + ~Frac(QQ['y']).0
        (x*y + 1)/y
        sage: MatrixSpace(ZZ['x'], 2, 2)(2) + ~Frac(QQ['x']).0                          # needs sage.modules
        [(2*x + 1)/x           0]
        [          0 (2*x + 1)/x]
        sage: f = ZZ['x,y,z'].0 + QQ['w,x,z,a'].0; f
        w + x
        sage: f.parent()
        Multivariate Polynomial Ring in w, x, y, z, a over Rational Field
        sage: ZZ['x,y,z'].0 + ZZ['w,x,z,a'].1
        2*x

    TESTS:

    Check that :issue:`8426` is fixed (see also :issue:`18076`)::

        sage: import numpy                                                              # needs numpy
        sage: if int(numpy.version.short_version[0]) > 1:                               # needs numpy
        ....:     numpy.set_printoptions(legacy="1.25")                                 # needs numpy

        sage: # needs sage.rings.real_mpfr
        sage: x = polygen(RR)
        sage: numpy.float32('1.5') * x                                                  # needs numpy
        1.50000000000000*x
        sage: x * numpy.float32('1.5')                                                  # needs numpy
        1.50000000000000*x
        sage: p = x**3 + 2*x - 1
        sage: p(float('1.2'))
        3.12800000000000
        sage: p(int('2'))
        11.0000000000000

    This used to fail (see :issue:`18076`)::

        sage: 1/3 + numpy.int8('12')                                                    # needs numpy
        37/3
        sage: -2/3 + numpy.int16('-2')                                                  # needs numpy
        -8/3
        sage: 2/5 + numpy.uint8('2')                                                    # needs numpy
        12/5

    The numpy types do not interact well with the Sage coercion framework. More
    precisely, if a numpy type is the first operand in a binary operation then
    this operation is done in numpy. The result is hence a numpy type::

        sage: numpy.uint8('2') + 3                                                      # needs numpy
        5
        sage: type(_)                                                                   # needs numpy
        <class 'numpy.int32'>  # 32-bit
        <class 'numpy.int64'>  # 64-bit

        sage: numpy.int8('12') + 1/3                                                    # needs numpy
        12.333333333333334
        sage: type(_)                                                                   # needs numpy
        <class 'numpy.float64'>

    AUTHOR:

    - Robert Bradshaw
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce import CoercionModel
            sage: cm = CoercionModel()
            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x^2 - 2, 'a')                                         # needs sage.rings.number_field
            sage: A = cm.get_action(ZZ, K, operator.mul)                                # needs sage.rings.number_field
            sage: f, g = cm.coercion_maps(QQ, int)
            sage: f, g = cm.coercion_maps(ZZ, int)
        """
        self.reset_cache()

    def reset_cache(self):
        """
        Clear the coercion cache.

        This should have no impact on the result of arithmetic operations, as
        the exact same coercions and actions will be re-discovered when needed.

        It may be useful for debugging, and may also free some memory.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: len(cm.get_cache()[0])    # random
            42
            sage: cm.reset_cache()
            sage: cm.get_cache()
            ({}, {})
        """
        # This MUST be a mapping of tuples, where each
        # tuple contains at least two elements that are either
        # None or of type Map.
        self._coercion_maps = TripleDict()
        # This MUST be a mapping to actions.
        self._action_maps = TripleDict()
        # This is a mapping from Parents to Parents, storing the result of division in the given parent.
        self._division_parents = TripleDict()

    def get_cache(self):
        """
        This returns the current cache of coercion maps and actions, primarily
        useful for debugging and introspection.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.canonical_coercion(1, 2/3)
            (1, 2/3)
            sage: maps, actions = cm.get_cache()

        Now let us see what happens when we do a binary operations with
        an integer and a rational::

            sage: left_morphism_ref, right_morphism_ref = maps[ZZ, QQ]

        Note that by :issue:`14058` the coercion model only stores a weak
        reference to the coercion maps in this case::

            sage: left_morphism_ref
            <weakref at ...; to 'sage.rings.rational.Z_to_Q' at ...>

        Moreover, the weakly referenced coercion map uses only a weak
        reference to the codomain::

            sage: left_morphism_ref()
            (map internal to coercion system -- copy before use)
            Natural morphism:
              From: Integer Ring
              To:   Rational Field

        To get an actual valid map, we simply copy the weakly referenced
        coercion map::

            sage: print(copy(left_morphism_ref()))
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: print(right_morphism_ref)
            None

        We can see that it coerces the left operand from an integer to a
        rational, and doesn't do anything to the right.

        Now for some actions::

            sage: R.<x> = ZZ['x']
            sage: 1/2 * x
            1/2*x
            sage: maps, actions = cm.get_cache()
            sage: act = actions[QQ, R, operator.mul]; act
            Left scalar multiplication by Rational Field
             on Univariate Polynomial Ring in x over Integer Ring
            sage: act.actor()
            Rational Field
            sage: act.domain()
            Univariate Polynomial Ring in x over Integer Ring
            sage: act.codomain()
            Univariate Polynomial Ring in x over Rational Field
            sage: act(1/5, x+10)
            1/5*x + 2
        """
        d1 = {(S, R): mors for (S, R, op), mors in self._coercion_maps.items()}
        d2 = self._action_maps.copy()
        return d1, d2

    def record_exceptions(self, bint value=True):
        r"""
        Enables (or disables) recording of the exceptions suppressed during
        arithmetic.

        Each time that record_exceptions is called (either enabling or disabling
        the record), the exception_stack is cleared.

        TESTS::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.record_exceptions()
            sage: cm._test_exception_stack()
            sage: cm.exception_stack()
            ['Traceback (most recent call last):\n  File "sage/structure/coerce.pyx", line ...TypeError: just a test']
            sage: cm.record_exceptions(False)
            sage: cm._test_exception_stack()
            sage: cm.exception_stack()
            []
        """
        self._record_exceptions = value
        self._exceptions_cleared = True
        self._exception_stack = []

    cpdef _record_exception(self):
        r"""
        Push the last exception that occurred onto the stack for later reference,
        for internal use.

        If the stack has not yet been flagged as cleared, we clear it now (rather
        than wasting time to do so for successful operations).

        TESTS::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.record_exceptions()
            sage: 1 + 1/2 + 2  # make sure there aren't any errors hanging around
            7/2
            sage: cm.exception_stack()
            []
            sage: cm._test_exception_stack()
            sage: cm.exception_stack()
            ['Traceback (most recent call last):\n  File "sage/structure/coerce.pyx", line ...TypeError: just a test']

        The function _test_exception_stack is executing the following code::

            try:
                raise TypeError("just a test")
            except TypeError:
                cm._record_exception()
        """
        if not self._record_exceptions:
            return
        if not self._exceptions_cleared:
            self._exception_stack = []
            self._exceptions_cleared = True
        self._exception_stack.append(traceback.format_exc().strip())

    def _test_exception_stack(self):
        r"""
        A function to test the exception stack.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.record_exceptions()
            sage: 1 + 1/11 # make sure there aren't any errors hanging around
            12/11
            sage: cm.exception_stack()
            []
            sage: cm._test_exception_stack()
            sage: cm.exception_stack()
            ['Traceback (most recent call last):\n  File "sage/structure/coerce.pyx", line ...TypeError: just a test']
        """
        try:
            raise TypeError("just a test")
        except TypeError:
            self._record_exception()

    def exception_stack(self):
        r"""
        Return the list of exceptions that were caught in the course of
        executing the last binary operation. Useful for diagnosis when
        user-defined maps or actions raise exceptions that are caught in
        the course of coercion detection.

        If all went well, this should be the empty list. If things aren't
        happening as you expect, this is a good place to check. See also
        :func:`coercion_traceback`.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.record_exceptions()
            sage: 1/2 + 2
            5/2
            sage: cm.exception_stack()
            []
            sage: 1/2 + GF(3)(2)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +:
            'Rational Field' and 'Finite Field of size 3'

        Now see what the actual problem was::

            sage: import traceback
            sage: cm.exception_stack()
            ['Traceback (most recent call last):...', 'Traceback (most recent call last):...']
            sage: print(cm.exception_stack()[-1])
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Rational Field' and 'Finite Field of size 3'

        This is typically accessed via the :func:`coercion_traceback` function.

        ::

            sage: coercion_traceback()
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Rational Field' and 'Finite Field of size 3'
        """
        if not self._exceptions_cleared:
            self._exception_stack = []
            self._exceptions_cleared = True
        return self._exception_stack

    def explain(self, xp, yp, op=mul, int verbosity=2):
        """
        This function can be used to understand what coercions will happen
        for an arithmetic operation between xp and yp (which may be either
        elements or parents). If the parent of the result can be determined
        then it will be returned.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()

            sage: cm.explain(ZZ, ZZ)
            Identical parents, arithmetic performed immediately.
            Result lives in Integer Ring
            Integer Ring

            sage: cm.explain(QQ, int)
            Coercion on right operand via
                Native morphism:
                  From: Set of Python objects of class 'int'
                  To:   Rational Field
            Arithmetic performed after coercions.
            Result lives in Rational Field
            Rational Field

            sage: R = ZZ['x']
            sage: cm.explain(R, QQ)
            Action discovered.
                Right scalar multiplication by Rational Field
                 on Univariate Polynomial Ring in x over Integer Ring
            Result lives in Univariate Polynomial Ring in x over Rational Field
            Univariate Polynomial Ring in x over Rational Field

            sage: cm.explain(ZZ['x'], QQ, operator.add)
            Coercion on left operand via
                Ring morphism:
                  From: Univariate Polynomial Ring in x over Integer Ring
                  To:   Univariate Polynomial Ring in x over Rational Field
                  Defn: Induced from base ring by
                        Natural morphism:
                          From: Integer Ring
                          To:   Rational Field
            Coercion on right operand via
                Polynomial base injection morphism:
                  From: Rational Field
                  To:   Univariate Polynomial Ring in x over Rational Field
            Arithmetic performed after coercions.
            Result lives in Univariate Polynomial Ring in x over Rational Field
            Univariate Polynomial Ring in x over Rational Field

        Sometimes with non-sage types there is not enough information to deduce
        what will actually happen::

            sage: R100 = RealField(100)                                                 # needs sage.rings.real_mpfr
            sage: cm.explain(R100, float, operator.add)                                 # needs sage.rings.real_mpfr
            Right operand is numeric, will attempt coercion in both directions.
            Unknown result parent.
            sage: parent(R100(1) + float(1))                                            # needs sage.rings.real_mpfr
            <class 'float'>
            sage: cm.explain(QQ, float, operator.add)
            Right operand is numeric, will attempt coercion in both directions.
            Unknown result parent.
            sage: parent(QQ(1) + float(1))
            <class 'float'>

        Special care is taken to deal with division::

            sage: cm.explain(ZZ, ZZ, operator.truediv)
            Identical parents, arithmetic performed immediately.
            Result lives in Rational Field
            Rational Field

            sage: ZZx = ZZ['x']
            sage: QQx = QQ['x']
            sage: cm.explain(ZZx, QQx, operator.truediv)
            Coercion on left operand via
                Ring morphism:
                  From: Univariate Polynomial Ring in x over Integer Ring
                  To:   Univariate Polynomial Ring in x over Rational Field
                  Defn: Induced from base ring by
                        Natural morphism:
                          From: Integer Ring
                          To:   Rational Field
            Arithmetic performed after coercions.
            Result lives in Fraction Field of Univariate Polynomial Ring in x over Rational Field
            Fraction Field of Univariate Polynomial Ring in x over Rational Field

            sage: cm.explain(int, ZZ, operator.truediv)
            Coercion on left operand via
                Native morphism:
                  From: Set of Python objects of class 'int'
                  To:   Integer Ring
            Arithmetic performed after coercions.
            Result lives in Rational Field
            Rational Field

            sage: cm.explain(ZZx, ZZ, operator.truediv)
            Action discovered.
                Right inverse action by Rational Field
                 on Univariate Polynomial Ring in x over Integer Ring
                with precomposition on right by Natural morphism:
                  From: Integer Ring
                  To:   Rational Field
            Result lives in Univariate Polynomial Ring in x over Rational Field
            Univariate Polynomial Ring in x over Rational Field

        .. NOTE::

           This function is accurate only in so far as :meth:`analyse` is kept
           in sync with the :meth:`bin_op` and
           :meth:`canonical_coercion` which are kept separate for
           maximal efficiency.
        """
        all, res = self.analyse(xp, yp, op)
        indent = " " * 4
        if verbosity >= 2:
            print("\n".join(s if isinstance(s, str)
                            else (indent + repr(s).replace("\n", "\n" + indent))
                            for s in all))
        elif verbosity >= 1:
            print("\n".join(s for s in all if isinstance(s, str)))
        if verbosity >= 1:
            if res is None:
                print("Unknown result parent.")
            else:
                print("Result lives in {}".format(res))
        return res

    cpdef analyse(self, xp, yp, op=mul):
        """
        Emulate the process of doing arithmetic between xp and yp, returning
        a list of steps and the parent that the result will live in.

        The :meth:`explain` method is easier to use, but if one wants access to
        the actual morphism and action objects (rather than their string
        representations), then this is the function to use.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: GF7 = GF(7)
            sage: steps, res = cm.analyse(GF7, ZZ)
            sage: steps
            ['Coercion on right operand via',
             Natural morphism:
              From: Integer Ring
              To:   Finite Field of size 7,
             'Arithmetic performed after coercions.']
            sage: res
            Finite Field of size 7
            sage: f = steps[1]; type(f)
            <class 'sage.rings.finite_rings.integer_mod.Integer_to_IntegerMod'>
            sage: f(100)
            2
        """
        self._exceptions_cleared = False
        res = None
        if not isinstance(xp, type) and not isinstance(xp, Parent):
            xp = parent(xp)
        if not isinstance(yp, type) and not isinstance(yp, Parent):
            yp = parent(yp)

        all = []
        if xp is yp:
            all.append("Identical parents, arithmetic performed immediately." % xp)
            if op is truediv and isinstance(xp, Parent):
                xp = self.division_parent(xp)
            return all, xp
        if xp == yp:
            all.append("Equal but distinct parents.")

        action = self.get_action(xp, yp, op)
        if action is not None:
            all.append("Action discovered.")
            all.append(action)
            return all, action.codomain()

        homs = self.discover_coercion(xp, yp)
        if homs is not None:
            x_mor, y_mor = homs
            if x_mor is not None:
                x_mor = x_mor.__copy__()
                all.append("Coercion on left operand via")
                all.append(x_mor)
                res = x_mor.codomain()
            if y_mor is not None:
                y_mor = y_mor.__copy__()
                all.append("Coercion on right operand via")
                all.append(y_mor)
                if res is not None and res is not y_mor.codomain():
                    raise RuntimeError("BUG in coercion model: codomains not equal", x_mor, y_mor)
                res = y_mor.codomain()
            all.append("Arithmetic performed after coercions.")
            if op is truediv and isinstance(res, Parent):
                res = self.division_parent(res)
            return all, res

        if isinstance(yp, Parent) and xp in [int, float, complex, bool]:
            mor = yp._internal_coerce_map_from(xp)
            if mor is not None:
                mor = mor.__copy__()
                all.append("Coercion on numeric left operand via")
                all.append(mor)
                if op is truediv and isinstance(yp, Parent):
                    yp = self.division_parent(yp)
                return all, yp
            all.append("Left operand is numeric, will attempt coercion in both directions.")
        elif type(xp) is type:
            all.append("Left operand is not Sage element, will try _sage_.")

        if isinstance(xp, Parent) and yp in [int, float, complex, bool]:
            mor = xp._internal_coerce_map_from(yp)
            if mor is not None:
                mor = mor.__copy__()
                all.append("Coercion on numeric right operand via")
                all.append(mor)
                if op is truediv and isinstance(xp, Parent):
                    xp = self.division_parent(xp)
                return all, xp
            all.append("Right operand is numeric, will attempt coercion in both directions.")
        elif type(yp) is type:
            all.append("Right operand is not Sage element, will try _sage_.")

        return all, None

    def common_parent(self, *args):
        """
        Compute a common parent for all the inputs.

        It's essentially an `n`-ary canonical coercion except it can
        operate on parents rather than just elements.

        INPUT:

        - ``args`` -- set of elements and/or parents

        OUTPUT:

        A :class:`Parent` into which each input should coerce, or raises a
        :exc:`TypeError` if no such :class:`Parent` can be found.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.common_parent(ZZ, QQ)
            Rational Field
            sage: cm.common_parent(ZZ, QQ, RR)                                          # needs sage.rings.real_mpfr
            Real Field with 53 bits of precision
            sage: ZZT = ZZ[['T']]
            sage: QQT = QQ['T']
            sage: cm.common_parent(ZZT, QQT, RDF)
            Power Series Ring in T over Real Double Field
            sage: cm.common_parent(4r, 5r)
            <class 'int'>
            sage: cm.common_parent(int, float, ZZ)
            <class 'float'>
            sage: real_fields = [RealField(prec) for prec in [10,20..100]]              # needs sage.rings.real_mpfr
            sage: cm.common_parent(*real_fields)                                        # needs sage.rings.real_mpfr
            Real Field with 10 bits of precision

        There are some cases where the ordering does matter, but if a parent
        can be found it is always the same::

            sage: QQxy = QQ['x,y']
            sage: QQyz = QQ['y,z']
            sage: cm.common_parent(QQxy, QQyz) == cm.common_parent(QQyz, QQxy)
            True
            sage: QQzt = QQ['z,t']
            sage: cm.common_parent(QQxy, QQyz, QQzt)
            Multivariate Polynomial Ring in x, y, z, t over Rational Field
            sage: cm.common_parent(QQxy, QQzt, QQyz)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Multivariate Polynomial Ring in x, y over Rational Field' and
            'Multivariate Polynomial Ring in z, t over Rational Field'
        """
        base = None
        for x in args:
            if not isinstance(x, Parent) and not isinstance(x, type):
                x = parent(x)
            if base is None:
                base = x
            if isinstance(base, Parent) and (<Parent>base).has_coerce_map_from(x):
                continue
            elif isinstance(x, Parent) and (<Parent>x).has_coerce_map_from(base):
                base = x
            else:
                a = base.an_element() if isinstance(base, Parent) else base(1)
                b = x.an_element() if isinstance(x, Parent) else x(1)
                base = parent(self.canonical_coercion(a, b)[0])
        return base

    cpdef division_parent(self, Parent P):
        r"""
        Deduces where the result of division in ``P`` lies by
        calculating the inverse of ``P.one()`` or ``P.an_element()``.

        The result is cached.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.division_parent(ZZ)
            Rational Field
            sage: cm.division_parent(QQ)
            Rational Field
            sage: ZZx = ZZ['x']
            sage: cm.division_parent(ZZx)
            Fraction Field of Univariate Polynomial Ring in x over Integer Ring
            sage: K = GF(41)
            sage: cm.division_parent(K)
            Finite Field of size 41
            sage: Zmod100 = Integers(100)
            sage: cm.division_parent(Zmod100)
            Ring of integers modulo 100
            sage: S5 = SymmetricGroup(5)                                                # needs sage.groups
            sage: cm.division_parent(S5)                                                # needs sage.groups
            Symmetric group of order 5! as a permutation group
        """
        try:
            return self._division_parents.get(P, None, None)
        except KeyError:
            pass
        try:
            ret = parent(~P.one())
        except Exception:
            self._record_exception()
            ret = parent(~P.an_element())
        self._division_parents.set(P, None, None, ret)
        return ret

    cpdef bin_op(self, x, y, op):
        """
        Execute the operation ``op`` on `x` and `y`.

        It first looks for an action
        corresponding to ``op``, and failing that, it tries to coerce `x` and `y`
        into a common parent and calls ``op`` on them.

        If it cannot make sense of the operation, a :exc:`TypeError` is raised.

        INPUT:

        - ``x`` -- the left operand

        - ``y`` -- the right operand

        - ``op`` -- a python function taking 2 arguments

          .. NOTE::

             ``op`` is often an arithmetic operation, but need not be so.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.bin_op(1/2, 5, operator.mul)
            5/2

        The operator can be any callable::

            sage: R.<x> = ZZ['x']
            sage: cm.bin_op(x^2 - 1, x + 1, gcd)
            x + 1

        Actions are detected and performed::

            sage: M = matrix(ZZ, 2, 2, range(4))                                        # needs sage.modules
            sage: V = vector(ZZ, [5,7])                                                 # needs sage.modules
            sage: cm.bin_op(M, V, operator.mul)                                         # needs sage.modules
            (7, 31)

        TESTS::

            sage: class Foo():
            ....:     def __rmul__(self, left):
            ....:         return 'hello'
            sage: H = Foo()
            sage: print(int(3)*H)
            hello
            sage: print(Integer(3)*H)
            hello
            sage: print(H*3)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: '<class '__main__.Foo'>' and 'Integer Ring'

            sage: class Nonsense():
            ....:     def __init__(self, s):
            ....:         self.s = s
            ....:     def __repr__(self):
            ....:         return self.s
            ....:     def __mul__(self, x):
            ....:         return Nonsense(self.s + chr(x%256))
            ....:     __add__ = __mul__
            ....:     def __rmul__(self, x):
            ....:         return Nonsense(chr(x%256) + self.s)
            ....:     __radd__ = __rmul__
            sage: a = Nonsense('blahblah')
            sage: a*80
            blahblahP
            sage: 80*a
            Pblahblah
            sage: a+80
            blahblahP
            sage: 80+a
            Pblahblah
        """
        self._exceptions_cleared = False

        # If parents are equal, we can just call op()
        xp = parent(x)
        yp = parent(y)
        if xp is yp:
            return op(x,y)

        # Actions take preference over common-parent coercions
        try:
            action = self._action_maps.get(xp, yp, op)
        except KeyError:
            action = self.get_action(xp, yp, op, x, y)
        if action is not None:
            if (<Action>action)._is_left:
                return (<Action>action)._act_(x, y)
            else:
                return (<Action>action)._act_(y, x)

        # Now coerce to a common parent and do the operation there
        try:
            xy = self.canonical_coercion(x, y)
        except TypeError:
            self._record_exception()
        else:
            return PyObject_CallObject(op, xy)

        if op is mul:
            # elements may also act on non-elements
            # (e.g. sequences or parents)
            if not isinstance(y, Element) or not isinstance(x, Element):
                try:
                    if hasattr(x, '_act_on_'):
                        res = x._act_on_(y, True)
                        if res is not None: return res
                except CoercionException:
                    self._record_exception()

                try:
                    if hasattr(x, '_acted_upon_'):
                        res = x._acted_upon_(y, True)
                        if res is not None: return res
                except CoercionException:
                    self._record_exception()

                try:
                    if hasattr(y, '_act_on_'):
                        res = y._act_on_(x, False)
                        if res is not None: return res
                except CoercionException:
                    self._record_exception()

                try:
                    if hasattr(y, '_acted_upon_'):
                        res = y._acted_upon_(x, False)
                        if res is not None: return res
                except CoercionException:
                    self._record_exception()

        if not isinstance(y, Element):
            op_name = op.__name__
            mul_method = getattr(y, '__r%s__'%op_name, None)
            if mul_method is not None:
                res = mul_method(x)
                if res is not None and res is not NotImplemented:
                    return res

        # We should really include the underlying error.
        # This causes so much headache.
        raise bin_op_exception(op, x, y)

    cpdef canonical_coercion(self, x, y):
        r"""
        Given two elements `x` and `y`, with parents `S` and `R` respectively,
        find a common parent `Z` such that there are coercions
        `f: S \to Z` and `g: R \to Z` and return `f(x), g(y)`,
        which will have the same parent.

        Raises a type error if no such `Z` can be found.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.canonical_coercion(mod(2, 10), 17)
            (2, 7)

            sage: # needs sage.modules
            sage: x, y = cm.canonical_coercion(1/2, matrix(ZZ, 2, 2, range(4)))
            sage: x
            [1/2   0]
            [  0 1/2]
            sage: y
            [0 1]
            [2 3]
            sage: parent(x) is parent(y)
            True

        There is some support for non-Sage datatypes as well::

            sage: x, y = cm.canonical_coercion(int(5), 10)
            sage: type(x), type(y)
            (<class 'sage.rings.integer.Integer'>, <class 'sage.rings.integer.Integer'>)

            sage: x, y = cm.canonical_coercion(int(5), complex(3))
            sage: type(x), type(y)
            (<class 'complex'>, <class 'complex'>)

            sage: class MyClass:
            ....:     def _sage_(self):
            ....:         return 13
            sage: a, b = cm.canonical_coercion(MyClass(), 1/3)
            sage: a, b
            (13, 1/3)
            sage: type(a)
            <class 'sage.rings.rational.Rational'>

        We also make an exception for 0, even if `\ZZ` does not map in::

            sage: canonical_coercion(vector([1, 2, 3]), 0)                              # needs sage.modules
            ((1, 2, 3), (0, 0, 0))
            sage: canonical_coercion(GF(5)(0), float(0))
            (0, 0)
        """
        xp = parent(x)
        yp = parent(y)
        if xp is yp:
            return x,y

        cdef Element x_elt, y_elt
        coercions = self.coercion_maps(xp, yp)
        if coercions is not None:
            x_map, y_map = coercions
            if x_map is not None:
                x_elt = (<Map>x_map)._call_(x)
            else:
                x_elt = x
            if y_map is not None:
                y_elt = (<Map>y_map)._call_(y)
            else:
                y_elt = y
            if x_elt is None:
                raise RuntimeError("BUG in map, returned None %s %s %s" % (x, type(x_map), x_map))
            elif y_elt is None:
                raise RuntimeError("BUG in map, returned None %s %s %s" % (y, type(y_map), y_map))
            if x_elt._parent is y_elt._parent:
                # We must verify this as otherwise we are prone to
                # getting into an infinite loop in c, and the above
                # maps may be written by (imperfect) users.
                return x_elt,y_elt
            elif x_elt._parent == y_elt._parent:
                # TODO: Non-uniqueness of parents strikes again!
                y_elt = parent(x_elt)(y_elt)
                if x_elt._parent is y_elt._parent:
                    return x_elt,y_elt
            self._coercion_error(x, x_map, x_elt, y, y_map, y_elt)

        cdef bint x_numeric = isinstance(x, (int, float, complex))
        cdef bint y_numeric = isinstance(y, (int, float, complex))

        if not x_numeric and is_numpy_type(type(x)):
            import numpy
            x_numeric = isinstance(x, numpy.number)
        if not y_numeric and is_numpy_type(type(y)):
            import numpy
            y_numeric = isinstance(y, numpy.number)

        if x_numeric and y_numeric:
            ty = type(x + y)
            return ty(x), ty(y)

        # Now handle the native python + sage object cases
        # that were not taken care of above.
        elif x_numeric:
            try:
                sage_parent = py_scalar_parent(type(x))
                if sage_parent is None or sage_parent.has_coerce_map_from(yp):
                    return x, x.__class__(y)
                else:
                    return self.canonical_coercion(sage_parent(x), y)
            except (TypeError, ValueError):
                self._record_exception()

        elif y_numeric:
            try:
                sage_parent = py_scalar_parent(type(y))
                if sage_parent is None or sage_parent.has_coerce_map_from(xp):
                    return y.__class__(x), y
                else:
                    return self.canonical_coercion(x, sage_parent(y))
            except (TypeError, ValueError):
                self._record_exception()

        # See if the non-objects define a _sage_ method.
        if not isinstance(x, SageObject) or not isinstance(y, SageObject):
            try:
                x = x._sage_()
                y = y._sage_()
            except AttributeError:
                self._record_exception()
            else:
                return self.canonical_coercion(x, y)

        # Allow coercion of 0 even if no coercion from Z
        if (x_numeric or is_Integer(x)) and not x and type(yp) is not type:
            try:
                return yp(0), y
            except Exception:
                self._record_exception()

        if (y_numeric or is_Integer(y)) and not y and type(xp) is not type:
            try:
                return x, xp(0)
            except Exception:
                self._record_exception()

        raise TypeError("no common canonical parent for objects with parents: '%s' and '%s'"%(xp, yp))

    cpdef coercion_maps(self, R, S):
        r"""
        Give two parents `R` and `S`, return a pair of coercion maps
        `f: R \rightarrow Z` and `g: S \rightarrow Z` , if such a `Z`
        can be found.

        In the (common) case that `R=Z` or `S=Z` then ``None`` is returned
        for `f` or `g` respectively rather than constructing (and subsequently
        calling) the identity morphism.

        If no suitable `f, g` can be found, a single ``None`` is returned.
        This result is cached.

        .. NOTE::

            By :issue:`14711`, coerce maps should be copied when using them
            outside of the coercion system, because they may become defunct
            by garbage collection.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: f, g = cm.coercion_maps(ZZ, QQ)
            sage: print(copy(f))
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: print(g)
            None

            sage: ZZx = ZZ['x']
            sage: f, g = cm.coercion_maps(ZZx, QQ)
            sage: print(f)
            (map internal to coercion system -- copy before use)
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in x over Rational Field
            sage: print(g)
            (map internal to coercion system -- copy before use)
            Polynomial base injection morphism:
              From: Rational Field
              To:   Univariate Polynomial Ring in x over Rational Field

            sage: K = GF(7)
            sage: cm.coercion_maps(QQ, K) is None
            True

        Note that to break symmetry, if there is a coercion map in both
        directions, the parent on the left is used::

            sage: # needs sage.modules
            sage: V = QQ^3
            sage: W = V.__class__(QQ, 3)
            sage: V == W
            True
            sage: V is W
            False
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.coercion_maps(V, W)
            (None,
             (map internal to coercion system -- copy before use)
             Coercion map:
               From: Vector space of dimension 3 over Rational Field
               To:   Vector space of dimension 3 over Rational Field)
            sage: cm.coercion_maps(W, V)
            (None,
             (map internal to coercion system -- copy before use)
             Coercion map:
               From: Vector space of dimension 3 over Rational Field
               To:   Vector space of dimension 3 over Rational Field)
            sage: v = V([1,2,3])
            sage: w = W([1,2,3])
            sage: parent(v + w) is V
            True
            sage: parent(w + v) is W
            True

        TESTS:

        We check that with :issue:`14058`, parents are still eligible for
        garbage collection after being involved in binary operations::

            sage: # needs sage.libs.pari
            sage: import gc
            sage: T = type(GF(2))
            sage: gc.collect()  # random
            852
            sage: N0 = len(list(o for o in gc.get_objects() if type(o) is T))
            sage: L = [ZZ(1) + GF(p)(1) for p in prime_range(2, 50)]
            sage: N1 = len(list(o for o in gc.get_objects() if type(o) is T))
            sage: N1 > N0
            True
            sage: del L
            sage: gc.collect()  # random
            3939
            sage: N2 = len(list(o for o in gc.get_objects() if type(o) is T))
            sage: N2 - N0
            0
        """
        try:
            refs = self._coercion_maps.get(R, S, None)
            if refs is None:
                return None
            R_map_ref, S_map_ref = refs
            if R_map_ref is None:
                S_map = <object>PyWeakref_GET_OBJECT(S_map_ref)
                if S_map is not None:
                    return None, S_map
            elif S_map_ref is None:
                R_map = <object>PyWeakref_GET_OBJECT(R_map_ref)
                if R_map is not None:
                    return R_map, None
            else:
                R_map = <object>PyWeakref_GET_OBJECT(R_map_ref)
                S_map = <object>PyWeakref_GET_OBJECT(S_map_ref)
                if R_map is not None and S_map is not None:
                    return R_map, S_map
        except KeyError:
            pass
        homs = self.discover_coercion(R, S)
        if 0:
            # This breaks too many things that are going to change
            # in the new coercion model anyways.
            # COERCE TODO: Enable it then.
            homs = self.verify_coercion_maps(R, S, homs)
        else:
            if homs is not None:
                x_map, y_map = homs
                if x_map is not None and not isinstance(x_map, Map):
                    raise RuntimeError("BUG in coercion model: coerce_map_from must return a Map")
                if y_map is not None and not isinstance(y_map, Map):
                    raise RuntimeError("BUG in coercion model: coerce_map_from must return a Map")
        if homs is None:
            refs = None
            swap = None
        else:
            R_map, S_map = homs
            R_map_ref = None if R_map is None else PyWeakref_NewRef(R_map, None)
            S_map_ref = None if S_map is None else PyWeakref_NewRef(S_map, None)
            refs = R_map_ref, S_map_ref
            if R_map is None and isinstance(S, Parent) and (<Parent>S).has_coerce_map_from(R):
                swap = None, PyWeakref_NewRef((<Parent>S).coerce_map_from(R), None)
            else:
                swap = S_map_ref, R_map_ref
        self._coercion_maps.set(R, S, None, refs)
        self._coercion_maps.set(S, R, None, swap)
        return homs

    cpdef verify_coercion_maps(self, R, S, homs, bint fix=False):
        """
        Make sure this is a valid pair of homomorphisms from `R` and `S` to a common parent.
        This function is used to protect the user against buggy parents.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: homs = QQ.coerce_map_from(ZZ), None
            sage: cm.verify_coercion_maps(ZZ, QQ, homs) == homs
            True
            sage: homs = QQ.coerce_map_from(ZZ), RR.coerce_map_from(QQ)
            sage: cm.verify_coercion_maps(ZZ, QQ, homs) == homs
            Traceback (most recent call last):
            ...
            RuntimeError: ('BUG in coercion model, codomains must be identical',
            Natural morphism:
              From: Integer Ring
              To:   Rational Field,
            Generic map:
              From: Rational Field
              To:   Real Field with 53 bits of precision)
        """
        if homs is None:
            return None
        cdef Map x_map, y_map
        R_map, S_map = homs
        if isinstance(R, type):
            R = Set_PythonType(R)
        elif isinstance(S, type):
            S = Set_PythonType(S)
        if R_map is None:
            R_map = IdentityMorphism(R)
        elif S_map is None:
            S_map = IdentityMorphism(S)
        # Make sure the domains are correct
        if R_map.domain() is not R:
            if fix:
                connecting = R_map.domain()._internal_coerce_map_from(R)
                if connecting is not None:
                    R_map = R_map * connecting
            if R_map.domain() is not R:
                raise RuntimeError("BUG in coercion model, left domain must be original parent", R, R_map)
        if S_map is not None and S_map.domain() is not S:
            if fix:
                connecting = S_map.domain()._internal_coerce_map_from(S)
                if connecting is not None:
                    S_map = S_map * connecting
            if S_map.domain() is not S:
                raise RuntimeError("BUG in coercion model, right domain must be original parent", S, S_map)
        # Make sure the codomains are correct
        if R_map.codomain() is not S_map.codomain():
            if fix:
                connecting = R_map.codomain()._internal_coerce_map_from(S_map.codomain())
                if connecting is not None:
                    S_map = connecting * S_map
                else:
                    connecting = S_map.codomain()._internal_coerce_map_from(R_map.codomain())
                    if connecting is not None:
                        R_map = connecting * R_map
            if R_map.codomain() is not S_map.codomain():
                raise RuntimeError("BUG in coercion model, codomains must be identical", R_map, S_map)
        if isinstance(R_map, IdentityMorphism):
            R_map = None
        elif isinstance(S_map, IdentityMorphism):
            S_map = None
        return R_map, S_map

    cpdef discover_coercion(self, R, S):
        """
        This actually implements the finding of coercion maps as described in
        the :meth:`coercion_maps` method.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()

        If R is S, then two identity morphisms suffice::

            sage: cm.discover_coercion(SR, SR)                                          # needs sage.symbolic
            (None, None)

        If there is a coercion map either direction, use that::

            sage: cm.discover_coercion(ZZ, QQ)
            ((map internal to coercion system -- copy before use)
            Natural morphism:
              From: Integer Ring
              To:   Rational Field, None)
            sage: cm.discover_coercion(RR, QQ)                                          # needs sage.rings.real_mpfr
            (None, (map internal to coercion system -- copy before use)
             Generic map:
              From: Rational Field
              To:   Real Field with 53 bits of precision)

        Otherwise, try and compute an appropriate cover::

            sage: ZZxy = ZZ['x,y']
            sage: cm.discover_coercion(ZZxy, RDF)
            ((map internal to coercion system -- copy before use)
            Call morphism:
              From: Multivariate Polynomial Ring in x, y over Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Real Double Field,
             (map internal to coercion system -- copy before use)
             Polynomial base injection morphism:
              From: Real Double Field
              To:   Multivariate Polynomial Ring in x, y over Real Double Field)

        Sometimes there is a reasonable "cover," but no canonical coercion::

            sage: sage.categories.pushout.pushout(QQ, QQ^3)                             # needs sage.modules
            Vector space of dimension 3 over Rational Field
            sage: print(cm.discover_coercion(QQ, QQ^3))                                 # needs sage.modules
            None
        """
        if R is S:
            return None, None

        # See if there is a natural coercion from R to S
        if isinstance(R, Parent):
            mor = (<Parent>R)._internal_coerce_map_from(S)
            if mor is not None:
                return None, mor

        # See if there is a natural coercion from S to R
        if isinstance(S, Parent):
            mor = (<Parent>S)._internal_coerce_map_from(R)
            if mor is not None:
                return mor, None

        # Try base extending
        if isinstance(R, Parent) and isinstance(S, Parent):
            from sage.categories.pushout import pushout
            try:
                Z = pushout(R, S)
                coerce_R = Z._internal_coerce_map_from(R)
                coerce_S = Z._internal_coerce_map_from(S)
                if coerce_R is None:
                    raise TypeError("No coercion from %s to pushout %s" % (R, Z))
                if coerce_S is None:
                    raise TypeError("No coercion from %s to pushout %s" % (S, Z))
                return coerce_R, coerce_S
            except Exception:
                self._record_exception()

        return None

    cpdef get_action(self, R, S, op=mul, r=None, s=None):
        """
        Get the action of R on S or S on R associated to the operation op.

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: ZZx = ZZ['x']
            sage: cm.get_action(ZZx, ZZ, operator.mul)
            Right scalar multiplication by Integer Ring
             on Univariate Polynomial Ring in x over Integer Ring
            sage: cm.get_action(ZZx, QQ, operator.mul)
            Right scalar multiplication by Rational Field
             on Univariate Polynomial Ring in x over Integer Ring
            sage: QQx = QQ['x']
            sage: cm.get_action(QQx, int, operator.mul)
            Right scalar multiplication by Integer Ring
             on Univariate Polynomial Ring in x over Rational Field
            with precomposition on right by Native morphism:
              From: Set of Python objects of class 'int'
              To:   Integer Ring

            sage: A = cm.get_action(QQx, ZZ, operator.truediv); A
            Right inverse action by Rational Field
             on Univariate Polynomial Ring in x over Rational Field
            with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: x = QQx.gen()
            sage: A(x+10, 5)
            1/5*x + 2
        """
        try:
            return self._action_maps.get(R, S, op)
        except KeyError:
            pass
        action = self.discover_action(R, S, op, r, s)
        action = self.verify_action(action, R, S, op)
        self._action_maps.set(R, S, op, action)
        return action

    cpdef verify_action(self, action, R, S, op, bint fix=True):
        r"""
        Verify that ``action`` takes an element of R on the left and S
        on the right, raising an error if not.

        This is used for consistency checking in the coercion model.

        EXAMPLES::

            sage: R.<x> = ZZ['x']
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.verify_action(R.get_action(QQ), R, QQ, operator.mul)
            Right scalar multiplication by Rational Field
             on Univariate Polynomial Ring in x over Integer Ring
            sage: cm.verify_action(R.get_action(QQ), RDF, R, operator.mul)
            Traceback (most recent call last):
            ...
            RuntimeError: There is a BUG in the coercion model:
                Action found for R <built-in function mul> S does not have the correct domains
                R = Real Double Field
                S = Univariate Polynomial Ring in x over Integer Ring
                (should be Univariate Polynomial Ring in x over Integer Ring, Rational Field)
                action = Right scalar multiplication by Rational Field
                         on Univariate Polynomial Ring in x over Integer Ring
                (<class 'sage.structure.coerce_actions.RightModuleAction'>)
        """
        if action is None:
            return action
        cdef bint ok = True
        try:
            if action.left_domain() is not R:
                ok &= isinstance(R, type) and action.left_domain()._type is R
            if action.right_domain() is not S:
                ok &= isinstance(S, type) and action.right_domain()._type is S
        except AttributeError:
            ok = False
        if not ok:
            if isinstance(R, type):
                R = Set_PythonType(R)
            if isinstance(S, type):
                S = Set_PythonType(S)

            # Non-unique parents
            if fix and action.left_domain() is not R and action.left_domain() == R:
                action = PrecomposedAction(action, action.left_domain()._internal_coerce_map_from(R), None)
            if fix and action.right_domain() is not S and action.right_domain() == S:
                action = PrecomposedAction(action, None, action.right_domain()._internal_coerce_map_from(S))

            if action.left_domain() is not R or action.right_domain() is not S:
                raise RuntimeError("""There is a BUG in the coercion model:
                Action found for R %s S does not have the correct domains
                R = %s
                S = %s
                (should be %s, %s)
                action = %s (%s)
                """ % (op, R, S, action.left_domain(), action.right_domain(), action, type(action)))

        return action

    cpdef discover_action(self, R, S, op, r=None, s=None):
        """
        INPUT:

        - ``R`` -- the left :class:`Parent` (or type)
        - ``S`` -- the right :class:`Parent` (or type)
        - ``op`` -- the operand, typically an element of the :mod:`operator` module
        - ``r`` -- (optional) element of `R`
        - ``s`` -- (optional) element of `S`

        OUTPUT: an action `A` such that `s` ``op`` `r` is given by `A(s,r)`

        The steps taken are illustrated below.

        EXAMPLES::

            sage: P.<x> = ZZ['x']
            sage: P.get_action(ZZ)
            Right scalar multiplication by Integer Ring on
             Univariate Polynomial Ring in x over Integer Ring
            sage: ZZ.get_action(P) is None
            True
            sage: cm = sage.structure.element.get_coercion_model()

        If `R` or `S` is a :class:`Parent`, ask it for an action by/on `R`::

            sage: cm.discover_action(ZZ, P, operator.mul)
            Left scalar multiplication by Integer Ring on
             Univariate Polynomial Ring in x over Integer Ring

        If `R` or `S` a type, recursively call :meth:`get_action`
        with the Sage versions of `R` and/or `S`::

            sage: cm.discover_action(P, int, operator.mul)
            Right scalar multiplication by Integer Ring on
             Univariate Polynomial Ring in x over Integer Ring
             with precomposition on right by Native morphism:
              From: Set of Python objects of class 'int'
              To:   Integer Ring

        If ``op`` is division, look for action on ``right`` by inverse::

            sage: cm.discover_action(P, ZZ, operator.truediv)
            Right inverse action by Rational Field on
             Univariate Polynomial Ring in x over Integer Ring
            with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Rational Field

        Check that :issue:`17740` is fixed::

            sage: R = GF(5)['x']
            sage: cm.discover_action(R, ZZ, operator.truediv)
            Right inverse action by Finite Field of size 5
             on Univariate Polynomial Ring in x over Finite Field of size 5
            with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Finite Field of size 5
            sage: cm.bin_op(R.gen(), 7, operator.truediv).parent()
            Univariate Polynomial Ring in x over Finite Field of size 5

        Check that :issue:`18221` is fixed::

            sage: # needs sage.combinat sage.modules
            sage: F.<x> = FreeAlgebra(QQ)
            sage: x / 2
            1/2*x
            sage: cm.discover_action(F, ZZ, operator.truediv)
            Right inverse action by Rational Field on
             Free Algebra on 1 generators (x,) over Rational Field
             with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        if isinstance(R, Parent):
            action = (<Parent>R).get_action(S, op, True, r, s)
            if action is not None:
                return action

        if isinstance(S, Parent):
            action = (<Parent>S).get_action(R, op, False, s, r)
            if action is not None:
                return action

        if type(R) is type:
            sageR = py_scalar_parent(R)
            if sageR is not None:
                action = self.get_action(sageR, S, op, s=s)
                if action is not None:
                    return PrecomposedAction(action, sageR._internal_coerce_map_from(R), None)

        if type(S) is type:
            sageS = py_scalar_parent(S)
            if sageS is not None:
                action = self.get_action(R, sageS, op, r=r)
                if action is not None:
                    return PrecomposedAction(action, None, sageS._internal_coerce_map_from(S))

        if op is truediv:
            # Division on right is the same acting on right by inverse, if it is so defined.
            right_mul = None
            try:
                right_mul = self.get_action(R, S, mul)
            except NotImplementedError:
                self._record_exception()

            if right_mul is not None and not right_mul.is_left():
                try:
                    action = ~right_mul
                    if action.right_domain() != S:
                        action = PrecomposedAction(action, None,
                                                   action.right_domain()._internal_coerce_map_from(S))
                    return action
                except TypeError: # action may not be invertible
                    self._record_exception()

            # It's possible an action is defined on the fraction field itself.
            try:
                K = S._pseudo_fraction_field()
            except AttributeError:
                pass
            else:
                if K is not S:
                    try:
                        right_mul = self.get_action(R, K, mul)
                    except NotImplementedError:
                        self._record_exception()

                    if right_mul is not None and not right_mul.is_left():
                        try:
                            return PrecomposedAction(~right_mul, None, K.coerce_map_from(S))
                        except TypeError: # action may not be invertible
                            self._record_exception()

        return None

    cpdef richcmp(self, x, y, int op):
        """
        Given two arbitrary objects ``x`` and ``y``, coerce them to
        a common parent and compare them using rich comparison operator
        ``op``.

        EXAMPLES::

            sage: from sage.structure.element import get_coercion_model
            sage: from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE
            sage: richcmp = get_coercion_model().richcmp
            sage: richcmp(None, None, op_EQ)
            True
            sage: richcmp(None, 1, op_LT)
            True
            sage: richcmp("hello", None, op_LE)
            False
            sage: richcmp(-1, 1, op_GE)
            False
            sage: richcmp(int(1), float(2), op_GE)
            False

        If there is no coercion, we only support ``==`` and ``!=``::

            sage: x = QQ.one(); y = GF(2).one()
            sage: richcmp(x, y, op_EQ)
            False
            sage: richcmp(x, y, op_NE)
            True
            sage: richcmp(x, y, op_GT)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for >:
            'Rational Field' and 'Finite Field of size 2'

        We support non-Sage types with the usual Python convention::

            sage: class AlwaysEqual():
            ....:     def __eq__(self, other):
            ....:         return True
            sage: x = AlwaysEqual()
            sage: x == 1
            True
            sage: 1 == x
            True
        """
        # Some very special cases
        if x is None or x is Ellipsis:
            return rich_to_bool(op, 0 if x is y else -1)
        if y is None or y is Ellipsis:
            return rich_to_bool(op, 0 if x is y else 1)

        cdef bint y_is_Element = isinstance(y, Element)

        # Check for manual __richcmp__ override (only on y since
        # x.__richcmp__ would already have been called)
        if y_is_Element:
            if (<Element>y)._parent.get_flag(Parent_richcmp_element_without_coercion):
                return Py_TYPE(y).tp_richcompare(y, x, revop(op))

        # Coerce to a common parent
        try:
            x, y = self.canonical_coercion(x, y)
        except (TypeError, NotImplementedError):
            pass
        else:
            # The common parent should not be one which explicitly
            # asked to *not* use coercion for comparisons.
            assert not (isinstance(x, Element) and
                (<Element>x)._parent.get_flag(Parent_richcmp_element_without_coercion))
            return PyObject_RichCompare(x, y, op)

        # Comparing with coercion didn't work, try something else.

        # Try y.__richcmp__(x, revop) where revop is the reversed
        # operation (<= becomes >=).
        # This only makes sense when y is not a Sage Element, otherwise
        # we would end up trying the same coercion again.
        if not y_is_Element and Py_TYPE(y).tp_richcompare:
            res = Py_TYPE(y).tp_richcompare(y, x, revop(op))
            if res is not NotImplemented:
                return res

        # At this point, we have 2 objects which cannot be coerced to
        # a common parent. So we assume that they are not equal.
        if op == Py_EQ:
            return False
        if op == Py_NE:
            return True

        # It does not make sense to compare x and y with an inequality,
        # so we raise an exception.
        if op == Py_LT:
            raise bin_op_exception('<', x, y)
        elif op == Py_LE:
            raise bin_op_exception('<=', x, y)
        elif op == Py_GT:
            raise bin_op_exception('>', x, y)
        else:
            raise bin_op_exception('>=', x, y)

    def _coercion_error(self, x, x_map, x_elt, y, y_map, y_elt):
        """
        This function is only called when someone has incorrectly implemented
        a user-defined part of the coercion system (usually, a morphism).

        EXAMPLES::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm._coercion_error('a', 'f', 'f(a)', 'b', 'g', 'g(b)')
            Traceback (most recent call last):
            ...
            RuntimeError: There is a bug in the coercion code in Sage.
            Both x (='f(a)') and y (='g(b)') are supposed to have identical parents but they don't.
            In fact, x has parent '<class 'str'>'
            whereas y has parent '<class 'str'>'
            Original elements 'a' (parent <class 'str'>) and 'b' (parent <class 'str'>) and maps
            <class 'str'> 'f'
            <class 'str'> 'g'
        """
        raise RuntimeError("""There is a bug in the coercion code in Sage.
Both x (=%r) and y (=%r) are supposed to have identical parents but they don't.
In fact, x has parent '%s'
whereas y has parent '%s'
Original elements %r (parent %s) and %r (parent %s) and maps
%s %r
%s %r""" % (x_elt, y_elt, parent(x_elt), parent(y_elt),
            x, parent(x), y, parent(y),
            type(x_map), x_map, type(y_map), y_map))


coercion_model = CoercionModel()
