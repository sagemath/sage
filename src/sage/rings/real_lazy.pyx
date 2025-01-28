# sage.doctest: needs sage.rings.real_interval_field sage.rings.real_mpfr
"""
Lazy real and complex numbers

These classes are very lazy, in the sense that it doesn't really do anything
but simply sits between exact rings of characteristic 0 and the real numbers.
The values are actually computed when they are cast into a field of fixed
precision.

The main purpose of these classes is to provide a place for exact rings (e.g.
number fields) to embed for the coercion model (as only one embedding can be
specified in the forward direction).

TESTS:

Bug :issue:`21991`::

    sage: a = QuadraticField(5).gen()                                                   # needs sage.rings.number_field
    sage: u = -573147844013817084101/2*a + 1281597540372340914251/2                     # needs sage.rings.number_field
    sage: RealIntervalField(128)(RLF(u)).is_exact()                                     # needs sage.rings.number_field
    False
"""

# ****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import math
import cmath

cdef add, sub, mul, truediv, pow, neg, inv
from operator import add, sub, mul, pow, neg, inv, truediv

cdef canonical_coercion
from sage.structure.element import canonical_coercion
from sage.structure.element cimport parent
from sage.structure.richcmp cimport richcmp

import sage.categories.map
from sage.categories.morphism cimport Morphism
from sage.rings.ring cimport Field
import sage.rings.infinity
from sage.rings.integer import Integer

cdef QQ, RR, CC, RealField, ComplexField
from sage.rings.rational_field import QQ

cdef late_import():
    global RR, CC, RealField, ComplexField
    if CC is not None:
        return
    try:
        from sage.rings.real_mpfr import RR, RealField
        from sage.rings.complex_mpfr import ComplexField
        from sage.rings.cc import CC
    except ImportError:
        pass

cdef _QQx = None

cdef QQx():
    global _QQx
    if _QQx is None:
        _QQx = QQ['x']
    return _QQx

cdef named_unops = [ 'sqrt', 'erf', 'gamma', 'abs',
                     'floor', 'ceil', 'trunc',
                     'exp', 'log', 'log10', 'log2',
                     'sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan',
                     'csc', 'sec', 'cot',
                     'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh' ]

cdef named_constants = [ 'pi', 'e',
                         'euler_constant', 'catalan_constant' ]

cdef class LazyField(Field):
    """
    The base class for lazy real fields.

    .. WARNING::

        LazyField uses :meth:`__getattr__`, to implement::

            sage: CLF.pi
            3.141592653589794?

        I (NT, 20/04/2012) did not manage to have ``__getattr__`` call
        :meth:`Parent.__getattr__` in case of failure; hence we can't
        use this ``__getattr__`` trick for extension types to recover
        the methods from categories. Therefore, at this point, no
        concrete subclass of this class should be an extension type
        (which is probably just fine)::

            sage: RLF.__class__
            <class 'sage.rings.real_lazy.RealLazyField_class_with_category'>
            sage: CLF.__class__
            <class 'sage.rings.real_lazy.ComplexLazyField_class_with_category'>
    """
    def __init__(self, base=None, names=None, normalize=True, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RLF # indirect doctest
            Real Lazy Field
        """
        Field.__init__(self,base or self, names=names, normalize=normalize, category=category)

    Element = LazyWrapper

    def __getattr__(self, name):
        """
        Simulates a list of methods found on the real/complex rings.

        EXAMPLES::

            sage: a = CLF.pi() * CLF.I(); a
            3.141592653589794?*I
            sage: CDF(a)
            3.141592653589793*I
        """
        if name in named_constants:
            return LazyConstant(self, name)
        elif name == 'I' and self == CLF:
            return LazyConstant(self, name)
        else:
            raise AttributeError(name)

    cpdef _coerce_map_from_(self, R):
        r"""
        The only things that coerce into this ring are exact rings that
        embed into `\RR` or `\CC` (depending on whether this field
        is real or complex), that is, exact rings that coerce into all
        rings into which this ring coerces.

        .. NOTE::

            The rings into which this ring coerces are currently the
            corresponding floating-point fields (RealField(p) or
            ComplexField(p)), machine-precision floating-point fields (RDF,
            CDF), and interval fields (RealIntervalField(p),
            ComplexIntervalField(p)). This method should be updated if a new
            parent is added that declares a coercion from RLF/CLF but not from
            one of these, otherwise coercions of elements of type LazyWrapper
            into the new ring might fail.

        EXAMPLES::

            sage: RLF.has_coerce_map_from(ZZ) # indirect doctest
            True
            sage: RLF.has_coerce_map_from(QQ)
            True
            sage: RLF.has_coerce_map_from(AA)                                           # needs sage.rings.number_field
            True
            sage: CLF.has_coerce_map_from(QQbar)                                        # needs sage.rings.number_field
            True
            sage: RLF.has_coerce_map_from(RDF)
            False

            sage: CLF.has_coerce_map_from(QQ)
            True
            sage: CLF.has_coerce_map_from(QQbar)                                        # needs sage.rings.number_field
            True
            sage: CLF.has_coerce_map_from(CC)
            False
            sage: CLF.has_coerce_map_from(RLF)
            True
        """
        if isinstance(R, type):
            if R is int:
                from sage.sets.pythonclass import Set_PythonType
                return LazyWrapperMorphism(Set_PythonType(R), self)
        elif R.is_exact():
            ivf = self.interval_field()
            mor = ivf.coerce_map_from(R)
            # Indirect coercions might lead to loops both in the coercion
            # discovery algorithm and when trying to convert LazyWrappers,
            # so we only consider direct coercions.
            if mor is not None and not isinstance(mor, sage.categories.map.FormalCompositeMap):
                mor = ivf.middle_field().coerce_map_from(R)
                if mor is not None and not isinstance(mor, sage.categories.map.FormalCompositeMap):
                    return LazyWrapperMorphism(R, self)
            # We can skip the test for a coercion to RDF/CDF since RR/CC
            # already coerce into it.

    def algebraic_closure(self):
        """
        Return the algebraic closure of ``self``, i.e., the complex lazy
        field.

        EXAMPLES::

            sage: RLF.algebraic_closure()
            Complex Lazy Field

            sage: CLF.algebraic_closure()
            Complex Lazy Field
        """
        return CLF

    cpdef interval_field(self, prec=None):
        """
        Abstract method to create the corresponding interval field.

        TESTS::

            sage: RLF.interval_field() # indirect doctest
            Real Interval Field with 53 bits of precision
        """
        raise NotImplementedError("subclasses must override this method")


class RealLazyField_class(LazyField):
    r"""
    This class represents the set of real numbers to unspecified precision.
    For the most part it simply wraps exact elements and defers evaluation
    until a specified precision is requested.

    Its primary use is to connect the exact rings (such as number fields) to
    fixed precision real numbers. For example, to specify an embedding of a
    number field `K` into `\RR` one can map into this field and the
    coercion will then be able to carry the mapping to real fields of any
    precision.

    EXAMPLES::

        sage: a = RLF(1/3)
        sage: a
        0.3333333333333334?
        sage: a + 1/5
        0.5333333333333334?
        sage: a = RLF(1/3)
        sage: a
        0.3333333333333334?
        sage: a + 5
        5.333333333333334?
        sage: RealField(100)(a+5)
        5.3333333333333333333333333333

    ::

        sage: CC.0 + RLF(1/3)
        0.333333333333333 + 1.00000000000000*I
        sage: ComplexField(200).0 + RLF(1/3)
        0.33333333333333333333333333333333333333333333333333333333333 + 1.0000000000000000000000000000000000000000000000000000000000*I

    TESTS::

        sage: TestSuite(RLF).run()
    """
    def interval_field(self, prec=None):
        """
        Return the interval field that represents the same mathematical
        field as ``self``.

        EXAMPLES::

            sage: RLF.interval_field()
            Real Interval Field with 53 bits of precision
            sage: RLF.interval_field(200)
            Real Interval Field with 200 bits of precision
        """
        from sage.rings.real_mpfi import RIF, RealIntervalField
        if prec is None:
            return RIF
        else:
            return RealIntervalField(prec)

    def construction(self):
        """
        Return the functorial construction of ``self``, namely, the
        completion of the rationals at infinity to infinite precision.

        EXAMPLES::

            sage: c, S = RLF.construction(); S
            Rational Field
            sage: RLF == c(S)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return CompletionFunctor(sage.rings.infinity.Infinity,
                                 sage.rings.infinity.Infinity,
                                 {'type': 'RLF'}), QQ

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(RLF) # indirect doctest
            \Bold{R}
        """
        return "\\Bold{R}"

    def gen(self, i=0):
        """
        Return the `i`-th generator of ``self``.

        EXAMPLES::

            sage: RLF.gen()
            1
        """
        if i == 0:
            return self(Integer(1))
        else:
            raise ValueError("RLF has only one generator.")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RealLazyField()
            Real Lazy Field
        """
        return "Real Lazy Field"

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: hash(RLF) == hash(RealLazyField())
            True
        """
        return 1501555429

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: RLF == loads(dumps(RLF))
            True
            sage: RLF is loads(dumps(RLF))
            True
        """
        return RealLazyField, ()


RLF = RealLazyField_class()


def RealLazyField():
    """
    Return the lazy real field.

    EXAMPLES:

    There is only one lazy real field::

        sage: RealLazyField() is RealLazyField()
        True
    """
    return RLF


class ComplexLazyField_class(LazyField):
    """
    This class represents the set of complex numbers to unspecified precision.
    For the most part it simply wraps exact elements and defers evaluation
    until a specified precision is requested.

    For more information, see the documentation of the
    :class:`RLF <sage.rings.real_lazy.RealLazyField_class>`.

    EXAMPLES::

        sage: a = CLF(-1).sqrt()
        sage: a
        1*I
        sage: CDF(a)
        1.0*I
        sage: ComplexField(200)(a)
        1.0000000000000000000000000000000000000000000000000000000000*I

    TESTS::

        sage: TestSuite(CLF).run()
    """
    def __init__(self):
        """
        This lazy field does not evaluate its elements until they are cast into
        a field of fixed precision.

        EXAMPLES::

            sage: a = RLF(1/3); a
            0.3333333333333334?
            sage: Reals(200)(a)
            0.33333333333333333333333333333333333333333333333333333333333
        """
        LazyField.__init__(self, base=RLF)
        self._populate_coercion_lists_(coerce_list=[LazyWrapperMorphism(RLF, self)])

    def interval_field(self, prec=None):
        """
        Return the interval field that represents the same mathematical
        field as ``self``.

        EXAMPLES::

            sage: CLF.interval_field()
            Complex Interval Field with 53 bits of precision
            sage: CLF.interval_field(333)                                               # needs sage.rings.complex_interval_field
            Complex Interval Field with 333 bits of precision
            sage: CLF.interval_field() is CIF
            True
        """
        if prec is None:
            from sage.rings.cif import CIF
            return CIF
        else:
            from sage.rings.complex_interval_field import ComplexIntervalField
            return ComplexIntervalField(prec)

    def gen(self, i=0):
        """
        Return the `i`-th generator of ``self``.

        EXAMPLES::

            sage: CLF.gen()
            1*I
            sage: ComplexField(100)(CLF.gen())                                          # needs sage.rings.number_field
            1.0000000000000000000000000000*I
        """
        if i == 0:
            from sage.rings.complex_double import CDF
            return LazyAlgebraic(self, [1, 0, 1], CDF.gen())
        else:
            raise ValueError("CLF has only one generator.")

    def construction(self):
        """
        Return the functorial construction of ``self``, namely,
        algebraic closure of the real lazy field.

        EXAMPLES::

            sage: c, S = CLF.construction(); S
            Real Lazy Field
            sage: CLF == c(S)
            True
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), RLF)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CLF) # indirect doctest
            \Bold{C}
        """
        return "\\Bold{C}"

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CLF
            Complex Lazy Field
        """
        return "Complex Lazy Field"

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_lazy import ComplexLazyField_class
            sage: hash(CLF) == hash(ComplexLazyField_class())
            True
        """
        return -1382606040

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: CLF == loads(dumps(CLF))
            True
            sage: CLF is loads(dumps(CLF))
            True
        """
        return ComplexLazyField, ()


CLF = ComplexLazyField_class()


def ComplexLazyField():
    """
    Return the lazy complex field.

    EXAMPLES:

    There is only one lazy complex field::

        sage: ComplexLazyField() is ComplexLazyField()
        True
    """
    return CLF


cdef int get_new_prec(R, int depth) except -1:
    """
    There are depth operations, so we want at least that many more digits of
    precision.

    Field creation may be expensive, so we want to avoid incrementing by 1 so
    that it is more likely for cached fields to be used.
    """
    cdef int needed_prec = R.prec()
    needed_prec += depth
    if needed_prec % 10 != 0:
        needed_prec += 10 - needed_prec % 10
    return needed_prec


cdef class LazyFieldElement(FieldElement):

    cpdef _add_(left, right):
        """
        Add ``left`` with ``right``.

        EXAMPLES::

            sage: RLF(5) + RLF(1/2) # indirect doctest
            5.5000000000000000?
        """
        if isinstance(left, LazyWrapper) and isinstance(right, LazyWrapper):
            try:
                return left._new_wrapper((<LazyWrapper?>left)._value + (<LazyWrapper?>right)._value)
            except TypeError:
                pass
        return left._new_binop(left, right, add)

    cpdef _sub_(left, right):
        """
        Subtract ``right`` from ``left``.

        EXAMPLES::

            sage: CLF(5) - 2 # indirect doctest
            3
        """
        if isinstance(left, LazyWrapper) and isinstance(right, LazyWrapper):
            try:
                return left._new_wrapper((<LazyWrapper?>left)._value - (<LazyWrapper?>right)._value)
            except TypeError:
                pass
        return left._new_binop(left, right, sub)

    cpdef _mul_(left, right):
        """
        Multiply ``left`` with ``right``.

        EXAMPLES::

            sage: CLF(10) * RLF(5) # indirect doctest
            50
        """
        if isinstance(left, LazyWrapper) and isinstance(right, LazyWrapper):
            try:
                return left._new_wrapper((<LazyWrapper?>left)._value * (<LazyWrapper?>right)._value)
            except TypeError:
                pass
        return left._new_binop(left, right, mul)

    cpdef _div_(left, right):
        """
        Divide ``left`` by ``right``.

        EXAMPLES::

            sage: a = RLF(1) / RLF(6); a # indirect doctest
            0.1666666666666667?
            sage: Reals(300)(a)
            0.166666666666666666666666666666666666666666666666666666666666666666666666666666666666666667
        """
        if isinstance(left, LazyWrapper) and isinstance(right, LazyWrapper):
            try:
                return left._new_wrapper((<LazyWrapper?>left)._value / (<LazyWrapper?>right)._value)
            except TypeError:
                pass
        return left._new_binop(left, right, truediv)

    def __pow__(left, right, dummy):
        """
        Raise ``left`` to the ``right`` power.

        EXAMPLES::

            sage: a = RLF(2) ^ (1/2); a
            1.414213562373095?
            sage: Reals(300)(a)
            1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753
        """
        if isinstance(left, LazyWrapper) and isinstance(right, LazyWrapper):
            try:
                return left._new_wrapper((<LazyWrapper>left)._value ** (<LazyWrapper>right)._value)
            except TypeError:
                pass
        if not isinstance(left, LazyFieldElement):
            left = (<LazyFieldElement>right)._new_wrapper(left)
        elif not isinstance(right, LazyFieldElement):
            right = (<LazyFieldElement>left)._new_wrapper(right)
        return (<LazyFieldElement>left)._new_binop(left, right, pow)

    def __neg__(self):
        """
        Return the negation of ``self``.

        EXAMPLES::

            sage: -RLF(7)
            -7
        """
        return self._new_unop(self, neg)

    def __invert__(self):
        """
        Take the reciprocal of ``self``.

        EXAMPLES::

          sage: a = ~RLF(6); a
          0.1666666666666667?
          sage: Reals(90)(a)
          0.16666666666666666666666667
        """
        return self._new_unop(self, inv)

    cpdef _richcmp_(self, other, int op):
        """
        If things are being wrapped, tries to compare values. That failing, it
        tries to compare intervals, which may return a false negative.

        EXAMPLES::

            sage: RLF(3) == RLF(9/3)
            True
            sage: RLF(3) == RLF(4)
            False
            sage: RLF(3) < RLF(5/3)
            False

        TESTS::

            sage: from sage.rings.real_lazy import LazyBinop
            sage: RLF(3) < LazyBinop(RLF, 5, 3, operator.truediv)
            False
            sage: from sage.rings.real_lazy import LazyWrapper
            sage: LazyWrapper(RLF, 3) < LazyWrapper(RLF, 5/3)
            False
            sage: from sage.rings.real_lazy import LazyUnop
            sage: RLF(3) < LazyUnop(RLF, 2, sqrt)
            False
            sage: from sage.rings.real_lazy import LazyNamedUnop
            sage: RLF(3) < LazyNamedUnop(RLF, 0, 'sin')
            False
            sage: from sage.rings.real_lazy import LazyConstant
            sage: RLF(3) < LazyConstant(RLF, 'e')
            False
        """
        left = self
        try:
            if isinstance(self, LazyWrapper) and isinstance(other, LazyWrapper):
                left, right = canonical_coercion((<LazyWrapper>self)._value, (<LazyWrapper>other)._value)
                return richcmp(left, right, op)
        except TypeError:
            pass
        left, right = self.approx(), other.approx()
        return richcmp(left.endpoints(), right.endpoints(), op)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: a = RLF(3)
            sage: hash(a)
            3
        """
        return hash(complex(self))

    cdef LazyFieldElement _new_wrapper(self, value):
        cdef LazyWrapper e = <LazyWrapper>LazyWrapper.__new__(LazyWrapper)
        e._parent = self._parent
        e._value = value
        return e

    cdef LazyFieldElement _new_binop(self, LazyFieldElement left, LazyFieldElement right, op):
        cdef LazyBinop e = <LazyBinop>LazyBinop.__new__(LazyBinop)
        e._parent = self._parent
        e._left = left
        e._right = right
        e._op = op
        return e

    cdef LazyFieldElement _new_unop(self, LazyFieldElement arg, op):
        cdef LazyUnop e = <LazyUnop>LazyUnop.__new__(LazyUnop)
        e._parent = self._parent
        e._op = op
        e._arg = arg
        return e

    def _repr_(self):
        """
        The string representation of ``self`` is an interval in which
        ``self`` is contained.

        EXAMPLES::

            sage: RLF(3) # indirect doctest
            3
            sage: RLF(1/3)
            0.3333333333333334?
        """
        return str(self.approx())

    def approx(self):
        """
        Return ``self`` as an element of an interval field.

        EXAMPLES::

            sage: CLF(1/6).approx()
            0.1666666666666667?
            sage: CLF(1/6).approx().parent()
            Complex Interval Field with 53 bits of precision

        When the absolute value is involved, the result might be real::

            sage: # needs sage.symbolic
            sage: z = exp(CLF(1 + I/2)); z
            2.38551673095914? + 1.303213729686996?*I
            sage: r = z.abs(); r
            2.71828182845905?
            sage: parent(z.approx())
            Complex Interval Field with 53 bits of precision
            sage: parent(r.approx())
            Real Interval Field with 53 bits of precision
        """
        return self.eval(self._parent.interval_field())

    def _real_double_(self, R):
        """
        Return ``self`` as a real double.

        EXAMPLES::

            sage: a = RLF(3)
            sage: RDF(a) # indirect doctest
            3.0
        """
        return self.eval(R)

    def _complex_double_(self, R):
        """
        Return ``self`` as a complex double.

        EXAMPLES::

            sage: a = RLF(5)
            sage: CDF(a) # indirect doctest
            5.0
            sage: a = CLF(-1)^(1/4)
            sage: CDF(a)
            0.7071067811865476 + 0.7071067811865475*I
        """
        return self.eval(R)

    def _generic_(self, R):
        """
        Return ``self`` in a generic ring ``R``.

        EXAMPLES::

            sage: a = RLF(2/3)
            sage: RR(a) # indirect doctest
            0.666666666666667
            sage: RR(a^2)
            0.444444444444444
        """
        return self.eval(R)

    _real_mpfi_ = _complex_mpfi_ =_mpfr_ = _complex_mpfr_field_ = _generic_

    def __complex__(self):
        """
        Return ``self`` as a complex.

        EXAMPLES::

            sage: complex(CLF(-1)^(1/4))
            (0.707106781186547...+0.707106781186547...j)
        """
        try:
            return self.eval(complex)
        except Exception:
            return complex(self.eval(ComplexField(53)))

    cpdef eval(self, R):
        """
        Abstract method for converting ``self`` into an element of ``R``.

        EXAMPLES::

            sage: a = RLF(12)
            sage: a.eval(ZZ)
            12
        """
        raise NotImplementedError("Subclasses must override this method.")

    cpdef int depth(self) noexcept:
        """
        Abstract method for returning the depth of ``self`` as an arithmetic
        expression.

        This is the  maximum number of dependent intermediate expressions when
        evaluating ``self``, and is used to determine the precision needed to
        get the final result to the desired number of bits.

        It is equal to the maximum of the right and left depths, plus one.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyBinop
            sage: a = LazyBinop(RLF, 6, 8, operator.mul)
            sage: a.depth()
            1
        """
        raise NotImplementedError("Subclasses must override this method.")

    def __dir__(self):
        """
        Add the named_unops to ``__dir__`` so that tab completion works.

        TESTS::

            sage: "log" in RLF(sqrt(8)).__dir__()                                       # needs sage.symbolic
            True
        """
        return FieldElement.__dir__(self) + named_unops

    def __getattribute__(self, name):
        """
        Simulates a list of methods found on the real/complex mpfr classes.

        EXAMPLES::

            sage: a = RLF(3)
            sage: a.sqrt()
            1.732050807568878?
            sage: sin(a)                                                                # needs sage.symbolic
            0.1411200080598673?
            sage: RealField(160)(tanh(RLF(3)))
            0.99505475368673045133188018525548847509781385470
        """
        if name in named_unops:
            return LazyNamedUnop(self._parent, self, name)
        else:
            return FieldElement.__getattribute__(self, name)

    def continued_fraction(self):
        r"""
        Return the continued fraction of ``self``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: a = RLF(sqrt(2)) + RLF(sqrt(3))
            sage: cf = a.continued_fraction()
            sage: cf
            [3; 6, 1, 5, 7, 1, 1, 4, 1, 38, 43, 1, 3, 2, 1, 1, 1, 1, 2, 4, ...]
            sage: cf.convergent(100)
            444927297812646558239761867973501208151173610180916865469/141414466649174973335183571854340329919207428365474086063
        """
        from sage.rings.continued_fraction import ContinuedFraction_real
        return ContinuedFraction_real(self)


def make_element(parent, *args):
    """
    Create an element of ``parent``.

    EXAMPLES::

        sage: a = RLF(pi) + RLF(sqrt(1/2))  # indirect doctest                          # needs sage.symbolic
        sage: bool(loads(dumps(a)) == a)                                                # needs sage.symbolic
        True
    """
    return parent(*args)


cdef class LazyWrapper(LazyFieldElement):

    cpdef int depth(self) noexcept:
        """
        Return the depth of ``self`` as an expression, which is always 0.

        EXAMPLES::

            sage: RLF(4).depth()
            0
        """
        return 0

    def __init__(self, LazyField parent, value, check=True):
        """
        A lazy element that simply wraps an element of another ring.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyWrapper
            sage: a = LazyWrapper(RLF, 3)
            sage: a._value
            3
        """
        FieldElement.__init__(self, parent)
        self._value = value
        if check:
            self._parent.interval_field()(value)

    def __neg__(self):
        """
        Return the negation of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyWrapper
            sage: a = LazyWrapper(RLF, 3)
            sage: (-a)._value
            -3
        """
        return self._new_wrapper(-self._value)

    def __invert__(self):
        """
        Return the reciprocal of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyWrapper
            sage: a = LazyWrapper(RLF, 23)
            sage: ~a
            0.04347826086956522?
            sage: (~a)._value
            1/23
        """
        return self._new_wrapper(~self._value)

    def __float__(self):
        """
        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyWrapper
            sage: a = LazyWrapper(CLF, 19)
            sage: float(a)
            19.0
        """
        return <double>self._value

    def __bool__(self):
        """
        Check to see if ``self`` is not zero.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyWrapper
            sage: not LazyWrapper(RLF, 1)
            False
            sage: not LazyWrapper(RLF, 0)
            True
        """
        return not not self._value

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: hash(CLF(-1))
            -2
            sage: hash(RLF(9/4)) == hash(9/4)
            True
        """
        return hash(self._value)

    cpdef eval(self, R):
        """
        Convert ``self`` into an element of ``R``.

        EXAMPLES::

            sage: a = RLF(12)
            sage: a.eval(ZZ)
            12
            sage: a.eval(ZZ).parent()
            Integer Ring
        """
        try:
            mor = R.convert_map_from(parent(self._value))
        except AttributeError:
            return R(self._value)
        if mor is not None and self.parent() not in mor.domains():
            return mor(self._value)
        else:
            raise TypeError("unable to convert {!r} to an element of {}".format(self._value, R))

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: a = RLF(2)
            sage: loads(dumps(a)) == a
            True
        """
        return make_element, (self._parent, self._value)

    def continued_fraction(self):
        r"""
        Return the continued fraction of ``self``.

        EXAMPLES::

            sage: a = RLF(sqrt(2))                                                      # needs sage.symbolic
            sage: a.continued_fraction()                                                # needs sage.symbolic
            [1; 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, ...]
        """
        from sage.rings.continued_fraction import ContinuedFraction_real, ContinuedFraction_infinite
        if isinstance(self._value, (ContinuedFraction_infinite, ContinuedFraction_real)):
            return self._value
        return ContinuedFraction_real(self)


cdef class LazyBinop(LazyFieldElement):

    def __init__(self, LazyField parent, left, right, op):
        """
        A lazy element representing a binary (usually arithmetic) operation
        between two other lazy elements.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyBinop
            sage: a = LazyBinop(RLF, 2, 1/3, operator.add)
            sage: a
            2.333333333333334?
            sage: Reals(200)(a)
            2.3333333333333333333333333333333333333333333333333333333333
        """
        FieldElement.__init__(self, parent)
        if not isinstance(left, LazyFieldElement):
            left = self._new_wrapper(left)
        if not isinstance(right, LazyFieldElement):
            right = self._new_wrapper(right)
        self._left = left
        self._right = right
        self._op = op

    cpdef int depth(self) noexcept:
        """
        Return the depth of ``self`` as an arithmetic expression.

        This is the  maximum number of dependent intermediate expressions when
        evaluating ``self``, and is used to determine the precision needed to
        get the final result to the desired number of bits.

        It is equal to the maximum of the right and left depths, plus one.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyBinop
            sage: a = LazyBinop(RLF, 6, 8, operator.mul)
            sage: a.depth()
            1
            sage: b = LazyBinop(RLF, 2, a, operator.sub)
            sage: b.depth()
            2
        """
        cdef int left = self._left.depth()
        cdef int right = self._right.depth()
        return 1 + (left if left > right else right)

    cpdef eval(self, R):
        """
        Convert the operands to elements of ``R``, then perform the operation
        on them.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyBinop
            sage: a = LazyBinop(RLF, 6, 8, operator.add)
            sage: a.eval(RR)
            14.0000000000000

        A bit absurd::

            sage: a.eval(str)
            '68'
        """
        left = self._left.eval(R)
        right = self._right.eval(R)
        if self._op is add:
            return left + right
        elif self._op is mul:
            return left * right
        elif self._op is sub:
            return left - right
        elif self._op is truediv:
            return left / right
        elif self._op is pow:
            return left ** right
        else:
            # We only do a call after testing the above because it is a python call.
            return self._op(left, right)

    def __float__(self):
        """
        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyBinop
            sage: a = LazyBinop(RLF, 3, 1/2, operator.sub)
            sage: float(a)
            2.5
            sage: type(float(a))
            <... 'float'>
        """
        cdef double left = self._left
        cdef double right = self._right
        if self._op is add:
            return left + right
        elif self._op is mul:
            return left * right
        elif self._op is sub:
            return left - right
        elif self._op is truediv:
            return left / right
        elif self._op is pow:
            return left ** right
        else:
            # We only do a call here because it is a python call.
            return self._op(left, right)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyBinop
            sage: a = LazyBinop(RLF, 5, 1/2, operator.sub)
            sage: b = LazyBinop(RLF, 4, 1/2, operator.add)
            sage: hash(a) == hash(b)
            False
        """
        return hash(self._op(hash(self._left), hash(self._right)))

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: from sage.rings.real_lazy import LazyBinop
            sage: a = LazyBinop(CLF, 3, 2, operator.truediv)
            sage: loads(dumps(a)) == a
            True
        """
        return make_element, (LazyBinop, self._parent, self._left, self._right, self._op)


cdef class LazyUnop(LazyFieldElement):

    def __init__(self, LazyField parent, arg, op):
        """
        Represent a unevaluated single function of one variable.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyUnop
            sage: a = LazyUnop(RLF, 3, sqrt); a
            1.732050807568878?
            sage: a._arg
            3
            sage: a._op
            <function sqrt at ...>
            sage: Reals(100)(a)
            1.7320508075688772935274463415
            sage: Reals(100)(a)^2
            3.0000000000000000000000000000
        """
        FieldElement.__init__(self, parent)
        if not isinstance(arg, LazyFieldElement):
            arg = self._new_wrapper(arg)
        self._op = op
        self._arg = arg

    cpdef int depth(self) noexcept:
        """
        Return the depth of ``self`` as an arithmetic expression.

        This is the  maximum number of dependent intermediate expressions when
        evaluating ``self``, and is used to determine the precision needed to
        get the final result to the desired number of bits.

        It is equal to one more than the depth of its operand.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyUnop
            sage: a = LazyUnop(RLF, 3, sqrt)
            sage: a.depth()
            1
            sage: b = LazyUnop(RLF, a, sin)
            sage: b.depth()
            2
        """
        return 1 + self._arg.depth()

    cpdef eval(self, R):
        """
        Convert ``self`` into an element of ``R``.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyUnop
            sage: a = LazyUnop(RLF, 3, sqrt)
            sage: a.eval(ZZ)                                                            # needs sage.symbolic
            sqrt(3)
        """
        arg = self._arg.eval(R)
        if self._op is neg:
            return -arg
        elif self._op is inv:
            return ~arg
        return self._op(self._arg.eval(R))

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: hash(RLF(sin(1))) == hash(RLF(sin(1)))                                # needs sage.symbolic
            True
        """
        return hash(self._op(hash(self._arg)))

    def __float__(self):
        """
        Convert ``self`` into a floating point.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyUnop
            sage: a = LazyUnop(RLF, 3, sqrt)
            sage: float(a)
            1.7320508075688772
        """
        return self._op(<double>self._arg)

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: from sage.rings.real_lazy import LazyUnop
            sage: a = LazyUnop(RLF, 7, sqrt)
            sage: float(loads(dumps(a))) == float(a)
            True
        """
        return make_element, (LazyUnop, self._parent, self._arg, self._op)


cdef class LazyNamedUnop(LazyUnop):

    def __init__(self, LazyField parent, arg, op, extra_args=None):
        """
        This class is used to represent the many named methods attached to real
        numbers, and is instantiated by the ``__getattr__`` method of
        :class:`LazyElements`.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyNamedUnop
            sage: a = LazyNamedUnop(RLF, 1, 'arcsin')
            sage: RR(a)
            1.57079632679490
            sage: a = LazyNamedUnop(RLF, 9, 'log', extra_args=(3,))
            sage: RR(a)
            2.00000000000000
        """
        LazyUnop.__init__(self, parent, arg, op)
        if extra_args is not None and not isinstance(extra_args, tuple):
            raise TypeError("extra args must be a tuple")
        self._extra_args = extra_args

    cpdef eval(self, R):
        """
        Convert ``self`` into an element of ``R``.

        TESTS::

            sage: from sage.rings.real_lazy import LazyNamedUnop
            sage: a = LazyNamedUnop(RLF, 4, 'sqrt')
            sage: RR(a) # indirect doctest
            2.00000000000000
            sage: a.sqrt()
            1.414213562373095?
            sage: RealField(212)(a)
            2.00000000000000000000000000000000000000000000000000000000000000
            sage: float(a)
            2.0

        Now for some extra arguments::

            sage: a = RLF(100)
            sage: a.log(10)
            2
            sage: float(a.log(10))
            2.0
        """
        arg = self._arg.eval(R)
        cdef bint has_extra_args = self._extra_args is not None and len(self._extra_args) > 0
        if type(R) is type:
            f = getattr(math, self._op)
            if has_extra_args:
                return f(arg, *self._extra_args)
            else:
                return f(arg)
        else:
            f = getattr(arg, self._op)
            if has_extra_args:
                return f(*self._extra_args)
            else:
                return f()

    def approx(self):
        """
        Do something reasonable with functions that are not defined on the
        interval fields.

        TESTS::

            sage: from sage.rings.real_lazy import LazyNamedUnop
            sage: LazyNamedUnop(RLF, 8, 'sqrt') # indirect doctest
            2.828427124746190?
        """
        try:
            return LazyUnop.approx(self)
        except AttributeError:
            # not everything defined on interval fields
            # this is less info though, but mostly just want to print it
            interval_field = self._parent.interval_field()
            return self.eval(interval_field.middle_field())

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyNamedUnop
            sage: a = LazyNamedUnop(RLF, 1, 'sin')
            sage: b = LazyNamedUnop(RLF, 1, 'cos')
            sage: hash(a) == hash(b)
            False
        """
        return hash(complex(self))

    def __float__(self):
        """
        TESTS::

            sage: from sage.rings.real_lazy import LazyNamedUnop
            sage: a = LazyNamedUnop(RLF, 1, 'sin')
            sage: float(a)
            0.8414709848078965
        """
        return self.eval(float)

    def __call__(self, *args):
        """
        TESTS::

            sage: a = RLF(32)
            sage: a.log(2)
            5
            sage: float(a.log(2))
            5.0

        What is going on here in the background is::

            sage: from sage.rings.real_lazy import LazyNamedUnop
            sage: b = LazyNamedUnop(RLF, a, 'log')
            sage: b(2)
            5
            sage: b(2)._extra_args
            (2,)
        """
        self._extra_args = args
        return self

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.rings.real_lazy import LazyNamedUnop
            sage: a = LazyNamedUnop(RLF, 1, 'sin')
            sage: float(loads(dumps(a))) == float(a)
            True
        """
        return make_element, (LazyNamedUnop, self._parent, self._arg, self._op, self._extra_args)

cdef class LazyConstant(LazyFieldElement):

    cdef readonly _name
    cdef readonly _extra_args

    def __init__(self, LazyField parent, name, extra_args=None):
        """
        This class represents a real or complex constant (such as ``pi``
        or ``I``).

        TESTS::

            sage: a = RLF.pi(); a
            3.141592653589794?
            sage: RealField(300)(a)
            3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482

            sage: from sage.rings.real_lazy import LazyConstant
            sage: a = LazyConstant(RLF, 'euler_constant')
            sage: RealField(200)(a)
            0.57721566490153286060651209008240243104215933593992359880577
        """
        LazyFieldElement.__init__(self, parent)
        self._name = name
        self._extra_args = extra_args

    cpdef eval(self, R):
        """
        Convert ``self`` into an element of ``R``.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyConstant
            sage: a = LazyConstant(RLF, 'e')
            sage: RDF(a) # indirect doctest
            2.718281828459045
            sage: a = LazyConstant(CLF, 'I')
            sage: CC(a)
            1.00000000000000*I

        TESTS:

        Check that :issue:`26839` is fixed::

            sage: RLF.pi().eval(float)
            3.141592653589793
            sage: type(RLF.pi().eval(float)) is float
            True

            sage: RLF.pi().eval(complex)
            (3.141592653589793+0j)
            sage: type(RLF.pi().eval(complex)) is complex
            True

            sage: RLF.pi().eval(RealBallField(128))
            [3.1415926535897932384626433832795028842 +/- 1.06e-38]

            sage: float(sin(RLF.pi()))
            1.2246467991473532e-16
        """
        # special handling of e and I
        if self._name == 'I':
            if R is float:
                raise ValueError('I is not a real number')
            elif R is complex:
                return 1j
            else:
                I = R.gen()
                if I*I != -R.one():
                    raise TypeError("The complex constant I is not in this complex field.")
                return I

        elif self._name == 'e':
            if R is float:
                return math.e
            elif R is complex:
                return complex(math.e)
            else:
                return R(1).exp()

        elif R is float:
            # generic float
            return getattr(math, self._name)
        elif R is complex:
            # generic complex
            return complex(getattr(cmath, self._name))
        else:
            # generic Sage parent
            f = getattr(R, self._name)
            if self._extra_args is None:
                return f()
            else:
                return f(*self._extra_args)

    def __call__(self, *args):
        """
        TESTS::

            sage: CLF.I()
            1*I
            sage: CDF(CLF.I())
            1.0*I
        """
        self._extra_args = args
        return self

    def __hash__(self):
        """
        Return the hash value of ``self``.

        TESTS::

            sage: from sage.rings.real_lazy import LazyConstant
            sage: a = LazyConstant(RLF, 'e')
            sage: hash(a) == hash(1)
            False
        """
        return hash(complex(self))

    def __float__(self):
        """
        TESTS::

            sage: from sage.rings.real_lazy import LazyConstant
            sage: a = LazyConstant(RLF, 'pi')
            sage: float(a)
            3.141592653589793
        """
        interval_field = self._parent.interval_field()
        return <double>self.eval(interval_field.middle_field())

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.rings.real_lazy import LazyConstant
            sage: a = LazyConstant(RLF, 'pi')
            sage: float(loads(dumps(a))) == float(a)
            True
        """
        return make_element, (LazyConstant, self._parent, self._name, self._extra_args)


cdef class LazyAlgebraic(LazyFieldElement):

    cdef readonly _poly
    cdef readonly _root_approx
    cdef readonly int _prec
    cdef readonly _quadratic_disc
    cdef readonly _root

    def __init__(self, parent, poly, approx, int prec=0):
        r"""
        This represents an algebraic number, specified by a polynomial over
        `\QQ` and a real or complex approximation.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: from sage.rings.real_lazy import LazyAlgebraic
            sage: a = LazyAlgebraic(RLF, x^2-2, 1.5)
            sage: a
            1.414213562373095?
        """
        LazyFieldElement.__init__(self, parent)
        self._poly = QQx()(poly)
        self._root = None
        if prec is None:
            prec = approx.parent().prec()
        self._prec = prec
        if self._poly.degree() == 2:
            c, b, a = self._poly.list()
            self._quadratic_disc = b*b - 4*a*c

        late_import()

        if isinstance(parent, RealLazyField_class):
            if not self._poly.number_of_real_roots():
                raise ValueError("%s has no real roots" % self._poly)
            approx = (RR if prec == 0 else RealField(prec))(approx)
        else:
            approx = (CC if prec == 0 else ComplexField(prec))(approx)
        self._root_approx = approx

    cpdef eval(self, R):
        """
        Convert ``self`` into an element of ``R``.

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyAlgebraic
            sage: a = LazyAlgebraic(CLF, QQ['x'].cyclotomic_polynomial(7), 0.6+0.8*CC.0)
            sage: a
            0.6234898018587335? + 0.7818314824680299?*I
            sage: ComplexField(150)(a)  # indirect doctest                              # needs sage.rings.number_field
            0.62348980185873353052500488400423981063227473 + 0.78183148246802980870844452667405775023233452*I

            sage: a = LazyAlgebraic(CLF, QQ['x'].0^2-7, -2.0)
            sage: RR(a)                                                                 # needs sage.rings.number_field
            -2.64575131106459
            sage: RR(a)^2                                                               # needs sage.rings.number_field
            7.00000000000000
        """
        if isinstance(R, type):
            if self._prec < 53:
                self.eval(self.parent().interval_field(64)) # up the prec
        elif R.is_exact() or self._prec < R.prec():
            # Carl Witty said:
            # Quadratic equation faster and more accurate than roots(),
            # but the current code doesn't do the right thing with interval
            # arithmetic (it returns a point interval) so it's being disabled
            # for now
#             if self._quadratic_disc is not None:
#                 c, b, a = self._poly.list()
#                 if self._root_approx.real() < -b/2*a:
#                     z = (-b - R(self._quadratic_disc).sqrt()) / (2*a)
#                 else:
#                     z = (-b + R(self._quadratic_disc).sqrt()) / (2*a)
#                 if z.parent() is not R:
#                     z = R(z)
#                 self._root_approx = z
#                 from sage.rings.complex_interval_field import is_IntervalField
#                 if is_IntervalField(R):
#                     self._root_approx = (self._root_approx.upper() + self._root_approx.lower()) / 2
#                 self._prec = R.prec()
#                 return R(self._root_approx)
            if self._root is None:
                # This could be done much more efficiently with Newton iteration,
                # but will require some care to make sure we get the right root, and
                # to the correct precision.
                from sage.rings.qqbar import AA, QQbar
                roots = self._poly.roots(ring = AA if isinstance(self._parent, RealLazyField_class) else QQbar)
                best_root = roots[0][0]
                min_dist = abs(self._root_approx - best_root)
                for r, e in roots[1:]:
                    dist = abs(self._root_approx - r)
                    if dist < min_dist:
                        best_root = r
                        min_dist = dist
                self._root = best_root
        if self._root is not None:
            return R(self._root)

    def __float__(self):
        """
        TESTS::

            sage: x = polygen(QQ)
            sage: from sage.rings.real_lazy import LazyAlgebraic
            sage: a = LazyAlgebraic(RLF, x^3-10, 1.5)
            sage: float(a)
            2.154434690031883...
        """
        return self.eval(float)

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.rings.real_lazy import LazyAlgebraic
            sage: x = polygen(QQ)
            sage: a = LazyAlgebraic(RLF, x^2 - 2, 1.5)
            sage: float(loads(dumps(a))) == float(a)
            True
        """
        return make_element, (LazyAlgebraic, self._parent, self._poly, self._root_approx, self._prec)


cdef class LazyWrapperMorphism(Morphism):

    def __init__(self, domain, LazyField codomain):
        """
        This morphism coerces elements from anywhere into lazy rings
        by creating a wrapper element (as fast as possible).

        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyWrapperMorphism
            sage: f = LazyWrapperMorphism(QQ, RLF)
            sage: a = f(3); a
            3
            sage: type(a)
            <class 'sage.rings.real_lazy.LazyWrapper'>
            sage: a._value
            3
            sage: a._value.parent()
            Rational Field
        """
        from sage.categories.homset import Hom
        Morphism.__init__(self, Hom(domain, codomain))

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.real_lazy import LazyWrapperMorphism
            sage: f = LazyWrapperMorphism(QQ, CLF)
            sage: a = f(1/3); a # indirect doctest
            0.3333333333333334?
            sage: type(a)
            <class 'sage.rings.real_lazy.LazyWrapper'>
            sage: Reals(100)(a)
            0.33333333333333333333333333333

        Note that it doesn't double-wrap lazy elements::

            sage: f = LazyWrapperMorphism(RLF, CLF)
            sage: x = RLF(20)
            sage: f(x)
            20
            sage: f(x)._value
            20
        """
        cdef LazyWrapper e = <LazyWrapper>LazyWrapper.__new__(LazyWrapper)
        e._parent = self._codomain
        if type(x) is LazyWrapper:
            e._value = (<LazyWrapper>x)._value
        else:
            e._value = x
        return e
