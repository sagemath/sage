# sage_setup: distribution = sagemath-objects
"""
Coerce maps
"""
import types

from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.sets.pythonclass cimport Set_PythonType

cdef object BuiltinMethodType = type(repr)

# COERCE TODO: remove or integrate better (as this bit is only checked on an error)
cdef bint print_warnings = 0


cdef class DefaultConvertMap(Map):
    """
    This morphism simply calls the codomain's element_constructor method,
    passing in the codomain as the first argument.

    EXAMPLES::

        sage: QQ[['x']].coerce_map_from(QQ)
        Coercion map:
          From: Rational Field
          To:   Power Series Ring in x over Rational Field
    """
    def __init__(self, domain, codomain, category=None):
        """
        TESTS:

        Maps of this type are morphisms in the category of sets with
        partial maps (see :issue:`15618`)::

            sage: f = GF(11).convert_map_from(GF(7)); f                                 # needs sage.rings.finite_rings
            Conversion map:
              From: Finite Field of size 7
              To:   Finite Field of size 11
            sage: f.parent()                                                            # needs sage.rings.finite_rings
            Set of Morphisms
             from Finite Field of size 7
             to Finite Field of size 11
             in Category of sets with partial maps

        Test that :issue:`23211` is resolved::

            sage: f._is_coercion                                                        # needs sage.rings.finite_rings
            False
            sage: QQ[['x']].coerce_map_from(QQ)._is_coercion
            True

        This class is deprecated when used directly::

            sage: from sage.structure.coerce_maps import DefaultConvertMap
            sage: DefaultConvertMap(ZZ, ZZ)
            doctest:...: DeprecationWarning: DefaultConvertMap is deprecated, use DefaultConvertMap_unique instead.
            This probably means that _element_constructor_ should be a method and not some other kind of callable
            See https://github.com/sagemath/sage/issues/26879 for details.
            Conversion map:
              From: Integer Ring
              To:   Integer Ring
        """
        # The base class DefaultConvertMap is deprecated, only the
        # derived class DefaultConvertMap_unique should be used.
        # When removing this deprecation, this class should be merged
        # into DefaultConvertMap_unique.
        if not isinstance(self, DefaultConvertMap_unique):
            from sage.misc.superseded import deprecation_cython as deprecation
            deprecation(26879, "DefaultConvertMap is deprecated, use DefaultConvertMap_unique instead. This probably means that _element_constructor_ should be a method and not some other kind of callable")

        if not isinstance(domain, Parent):
            domain = Set_PythonType(domain)
        if category is None:
            from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
            category = SetsWithPartialMaps()
        parent = domain.Hom(codomain, category=category)
        Map.__init__(self, parent)
        self._coerce_cost = 100
        if (<Parent>codomain)._element_constructor is None:
            raise RuntimeError("BUG in coercion model, no element constructor for {}".format(type(codomain)))

    def _repr_type(self):
        r"""
        Return a printable type for this morphism.

        EXAMPLES::

            sage: f = GF(11).convert_map_from(GF(7))                                    # needs sage.rings.finite_rings
            sage: f._repr_type()                                                        # needs sage.rings.finite_rings
            'Conversion'
        """
        return self._repr_type_str or ("Coercion" if self._is_coercion else "Conversion")

    cpdef Element _call_(self, x):
        """
        Create an element of the codomain from a single element of the domain.

        EXAMPLES::

            sage: f = QQ[['x']].coerce_map_from(QQ)
            sage: f(2/3).parent()
            Power Series Ring in x over Rational Field
        """
        cdef Parent C = self._codomain
        try:
            return C._element_constructor(C, x)
        except Exception:
            if print_warnings:
                print(type(C), C)
                print(type(C._element_constructor), C._element_constructor)
            raise

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        Create an element of the codomain from an element of the domain, with extra arguments.

        EXAMPLES::

            sage: f = QQ[['x']].coerce_map_from(QQ)
            sage: f(2/3, 4)
            2/3 + O(x^4)
        """
        cdef Parent C = self._codomain
        try:
            if len(args) == 0:
                if len(kwds) == 0:
                    # This line is apparently never used in any tests (hivert, 2009-04-28)
                    return C._element_constructor(C, x)
                else:
                    return C._element_constructor(C, x, **kwds)
            else:
                if len(kwds) == 0:
                    return C._element_constructor(C, x, *args)
                else:
                    return C._element_constructor(C, x, *args, **kwds)
        except Exception:
            if print_warnings:
                print(type(C), C)
                print(type(C._element_constructor), C._element_constructor)
            raise


cdef class DefaultConvertMap_unique(DefaultConvertMap):
    """
    This morphism simply defers action to the codomain's
    element_constructor method, WITHOUT passing in the codomain as the
    first argument.

    This is used for creating elements that don't take a parent as the
    first argument to their __init__ method, for example, Integers,
    Rationals, Algebraic Reals... all have a unique parent. It is also
    used when the element_constructor is a bound method (whose self
    argument is assumed to be bound to the codomain).
    """
    cpdef Element _call_(self, x):
        cdef Parent C = self._codomain
        try:
            return C._element_constructor(x)
        except Exception:
            if print_warnings:
                print(type(C), C)
                print(type(C._element_constructor), C._element_constructor)
            raise

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        cdef Parent C = self._codomain
        try:
            if len(args) == 0:
                if len(kwds) == 0:
                    return C._element_constructor(x)
                else:
                    return C._element_constructor(x, **kwds)
            else:
                if len(kwds) == 0:
                    return C._element_constructor(x, *args)
                else:
                    return C._element_constructor(x, *args, **kwds)
        except Exception:
            if print_warnings:
                print(type(C), C)
                print(type(C._element_constructor), C._element_constructor)
            raise


cdef class NamedConvertMap(Map):
    """
    This is used for creating elements via the _xxx_ methods.

    For example, many elements implement an _integer_ method to
    convert to ZZ, or a _rational_ method to convert to QQ.
    """

    def __init__(self, domain, codomain, method_name):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.structure.coerce_maps import NamedConvertMap
            sage: var('t')
            t
            sage: mor = NamedConvertMap(SR, QQ['t'], '_polynomial_')
            sage: mor(t^2/4 + 1)
            1/4*t^2 + 1
            sage: mor = NamedConvertMap(SR, GF(7)[['t']], '_polynomial_')
            sage: mor(t^2/4 + 1)
            1 + 2*t^2
        """
        if isinstance(domain, type):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, codomain)
        self._coerce_cost = 400
        self.method_name = method_name
        self._repr_type_str = "Conversion via %s method" % self.method_name

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.structure.coerce_maps import NamedConvertMap
            sage: var('t')
            t
            sage: phi = NamedConvertMap(SR, QQ['t'], '_polynomial_')
            sage: psi = copy(phi)    # indirect doctest
            sage: psi
            Conversion via _polynomial_ method map:
              From: Symbolic Ring
              To:   Univariate Polynomial Ring in t over Rational Field
            sage: phi == psi                    # not implemented
            True
            sage: psi(t^2/4 + 1)
            1/4*t^2 + 1
            sage: psi(t^2/4 + 1) == phi(t^2/4 + 1)
            True
        """
        slots = Map._extra_slots(self)
        slots['method_name'] = self.method_name
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.structure.coerce_maps import NamedConvertMap
            sage: var('t')
            t
            sage: phi = NamedConvertMap(SR, QQ['t'], '_polynomial_')
            sage: psi = copy(phi)    # indirect doctest
            sage: psi
            Conversion via _polynomial_ method map:
              From: Symbolic Ring
              To:   Univariate Polynomial Ring in t over Rational Field
            sage: phi == psi                    # not implemented
            True
            sage: psi(t^2/4 + 1)
            1/4*t^2 + 1
            sage: psi(t^2/4 + 1) == phi(t^2/4 + 1)
            True
        """
        self.method_name = _slots['method_name']
        Map._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_maps import NamedConvertMap
            sage: f = NamedConvertMap(GF(5), QQ, '_integer_'); f
            Conversion via _integer_ method map:
              From: Finite Field of size 5
              To:   Rational Field
            sage: f(19)
            4
            sage: f(19).parent()
            Rational Field
        """
        cdef Parent C = self._codomain
        try:
            method = getattr(x, self.method_name)
        except AttributeError:
            if print_warnings:
                print(type(x), x)
                print(type(C), C)
                print(self.method_name)
            raise TypeError(f"cannot coerce {x} to {C}")
        cdef Map m
        cdef Element e = method(C)
        if e is None:
            raise RuntimeError("BUG in coercion model: {} method of {} returned None".format(self.method_name, type(x)))
        if e._parent is not C:
            m = C._internal_convert_map_from(e._parent)
            if m is None or m is self:
                raise TypeError
            e = m._call_(e)
        return e

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_maps import NamedConvertMap
            sage: f = NamedConvertMap(SR, ZZ['x'], '_polynomial_')                      # needs sage.symbolic
            sage: f(x^2 + 1, check=True)                                                # needs sage.symbolic
            x^2 + 1
        """
        cdef Parent C = self._codomain
        return C._element_constructor(self._call_(x), *args, **kwds)


# Perhaps this could be a method, extracting (<PyMethodDescrObject *>(<object>Parent).coerce_map_from).d_method.ml_meth and/or PyCFunction_GET_FUNCTION(method)
# and constructing a CCallableConvertMap_class if it is bound to the codomain.

cdef class CallableConvertMap(Map):
    def __init__(self, domain, codomain, func, parent_as_first_arg=None):
        r"""
        This lets one easily create maps from any callable object.

        This is especially useful to create maps from bound methods.

        EXAMPLES::

            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x): return x/2
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(3)
            3/2
            sage: f
            Conversion via foo map:
              From: Integer Ring
              To:   Rational Field

        Create a homomorphism from `\RR` to `\RR^+` viewed as additive groups.

        ::

            sage: # needs sage.symbolic
            sage: f = CallableConvertMap(RR, RR, exp, parent_as_first_arg=False)
            sage: f(0)
            1.00000000000000
            sage: f(1)
            2.71828182845905
            sage: f(-3)
            0.0497870683678639
        """
        if isinstance(domain, type):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, codomain)
        self._coerce_cost = 100
        self._func = func
        if parent_as_first_arg is None:
            if isinstance(func, types.MethodType):
                # can't easily access self
                parent_as_first_arg = False
            elif isinstance(func, BuiltinMethodType):
                parent_as_first_arg = codomain is func.__self__
            else:
                parent_as_first_arg = True
        self._parent_as_first_arg = parent_as_first_arg
        try:
            self._repr_type_str = "Conversion via %s" % self._func.__name__
        except AttributeError:
            self._repr_type_str = "Conversion via %s" % self._func

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x): return x^2
            sage: f = CallableConvertMap(ZZ, ZZ, foo)
            sage: g = copy(f)     # indirect doctest
            sage: f == g          # not implemented (todo: implement comparison)
            True
            sage: f(3) == g(3)
            True
        """
        slots = Map._extra_slots(self)
        slots['_func'] = self._func
        slots['_parent_as_first_arg'] = self._parent_as_first_arg
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x): return x^2
            sage: f = CallableConvertMap(ZZ, ZZ, foo)
            sage: g = copy(f)     # indirect doctest
            sage: f == g          # not implemented (todo: implement comparison)
            True
            sage: f(3) == g(3)
            True
        """
        self._func = _slots['_func']
        self._parent_as_first_arg = _slots['_parent_as_first_arg']
        Map._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Because self._func may be anything we do a little bit of sanity
        checking (the return value must be an element with the correct parent).

        TESTS::

            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x): return x
            sage: f = CallableConvertMap(ZZ, ZZ, foo)
            sage: f(0)
            0
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(0)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned element with wrong parent (expected Rational Field got Integer Ring)
            sage: def foo(P, x): return None
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(0)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned None
        """
        cdef Element y
        cdef Parent C = self._codomain
        try:
            if self._parent_as_first_arg:
                y = self._func(C, x)
            else:
                y = self._func(x)
        except Exception:
            if print_warnings:
                print(self._func)
                print(C)
            raise
        if y is None:
            raise RuntimeError("BUG in coercion model: {} returned None".format(self._func))
        elif y._parent is not C:
            raise RuntimeError("BUG in coercion model: {} returned element with wrong parent (expected {} got {})".format(self._func, C, y._parent))
        return y

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        TESTS::

            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x, y): return x or y
            sage: f = CallableConvertMap(ZZ, ZZ, foo)
            sage: f(0, 3)
            3
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(0, 3)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned element with wrong parent (expected Rational Field got Integer Ring)
            sage: f(None, None)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned None
        """
        cdef Element y
        try:
            if self._parent_as_first_arg:
                y = self._func(self._codomain, x, *args, **kwds)
            else:
                y = self._func(x, *args, **kwds)
        except Exception:
            if print_warnings:
                print(self._func)
                print(self._codomain)
            raise
        if y is None:
            raise RuntimeError("BUG in coercion model: %s returned None" % (self._func))
        elif y._parent is not self._codomain:
            raise RuntimeError("BUG in coercion model: %s returned element with wrong parent (expected %s got %s)" % (self._func, self._codomain, y._parent))
        return y


cdef class CCallableConvertMap_class(Map):
    cdef Element (*_func)(Parent, object)
    cdef public _name

    def __init__(self, domain, codomain, name):
        if isinstance(domain, type):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, codomain)
        self._coerce_cost = 10
        self._name = name

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: from sage.structure.coerce_maps import test_CCallableConvertMap
            sage: f = test_CCallableConvertMap(QQ, 'test')
            sage: f(1/3)
            -8/27
        """
        return self._func(self._codomain, x)

    def _repr_type(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_maps import test_CCallableConvertMap
            sage: test_CCallableConvertMap(ZZ, 'any name')
            Conversion via c call 'any name' map:
              From: Integer Ring
              To:   Integer Ring
            sage: test_CCallableConvertMap(ZZ, None)
            Conversion via c call at 0x... map:
              From: Integer Ring
              To:   Integer Ring
        """
        if self._name is None:
            return "Conversion via c call at 0x%x" % <long>self._func
        else:
            return "Conversion via c call '%s'" % self._name


cdef Map CCallableConvertMap(domain, codomain, void* func, name):
    """
    Use this to create a map from domain to codomain by calling func
    (which must be a function pointer taking a Parent and object, and
    returning an Element in the given Parent).

    This is the c analogue of CallableConvertMap.
    """
    # Cython doesn't yet accept function pointers as arguments,
    # change this when it does.
    cdef CCallableConvertMap_class map = CCallableConvertMap_class(domain, codomain, name)
    map._func = <Element (*)(Parent, object)>func
    return map

cpdef Element _ccall_test_function(codomain, x):
    """
    For testing CCallableConvertMap_class. Returns x*x*x-x in the codomain.

    TESTS::

        sage: from sage.structure.coerce_maps import _ccall_test_function
        sage: _ccall_test_function(ZZ, 1)
        0
        sage: _ccall_test_function(ZZ, 2)
        6
        sage: _ccall_test_function(ZZ, -3)
        -24
    """
    return codomain(x*x*x-x)


def test_CCallableConvertMap(domain, name=None):
    """
    For testing CCallableConvertMap_class.

    TESTS::

        sage: from sage.structure.coerce_maps import test_CCallableConvertMap
        sage: f = test_CCallableConvertMap(ZZ, 'test'); f
        Conversion via c call 'test' map:
          From: Integer Ring
          To:   Integer Ring
        sage: f(3)
        24
        sage: f(9)
        720
    """
    return CCallableConvertMap(domain, domain, <void*>&_ccall_test_function, name)


cdef class ListMorphism(Map):

    cdef Map _real_morphism

    def __init__(self, domain, Map real_morphism):
        if not isinstance(domain, Parent):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, real_morphism.codomain())
        self._coerce_cost = real_morphism._coerce_cost + 3
        self._real_morphism = real_morphism
        self._repr_type_str = "List"

    cdef dict _extra_slots(self):
        slots = Map._extra_slots(self)
        slots['_real_morphism'] = self._real_morphism
        return slots

    cdef _update_slots(self, dict _slots):
        self._real_morphism = _slots['_real_morphism']
        Map._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        try:
            x = x._data
        except AttributeError:
            x = list(x)
        return self._real_morphism._call_(x)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        try:
            x = x._data
        except AttributeError:
            x = list(x)
        return self._real_morphism._call_with_args(x, args, kwds)


cdef class TryMap(Map):
    def __init__(self, Map morphism_preferred, Map morphism_backup, error_types=None):
        """
        TESTS::

            sage: sage.structure.coerce_maps.TryMap(RDF.coerce_map_from(QQ), RDF.coerce_map_from(ZZ))
            Traceback (most recent call last):
            ...
            TypeError: incorrectly matching parent
        """
        if (morphism_preferred.domain() is not morphism_backup.domain()
             or morphism_preferred.codomain() is not morphism_backup.codomain()):
            raise TypeError("incorrectly matching parent")
        Map.__init__(self, morphism_preferred.parent())
        self._map_p = morphism_preferred
        self._map_b = morphism_backup
        if error_types is None:
            self._error_types = (ValueError, TypeError, AttributeError)
        else:
            self._error_types = error_types

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: map1 = sage.structure.coerce_maps.CallableConvertMap(ZZ, QQ, lambda parent, x: 1/x)
            sage: map2 = QQ.coerce_map_from(ZZ)
            sage: map = sage.structure.coerce_maps.TryMap(map1, map2, error_types=(ZeroDivisionError,))
            sage: cmap = copy(map)     # indirect doctest
            sage: cmap == map          # not implemented (todo: implement comparison)
            True
            sage: map(3) == cmap(3)
            True
            sage: map(0) == cmap(0)
            True
        """
        slots = Map._extra_slots(self)
        slots['_map_p'] = self._map_p
        slots['_map_b'] = self._map_b
        slots['_error_types'] = self._error_types
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: map1 = sage.structure.coerce_maps.CallableConvertMap(ZZ, QQ, lambda parent, x: 1/x)
            sage: map2 = QQ.coerce_map_from(ZZ)
            sage: map = sage.structure.coerce_maps.TryMap(map1, map2, error_types=(ZeroDivisionError,))
            sage: cmap = copy(map)     # indirect doctest
            sage: cmap == map          # not implemented (todo: implement comparison)
            True
            sage: map(3) == cmap(3)
            True
            sage: map(0) == cmap(0)
            True
        """
        self._map_p = _slots['_map_p']
        self._map_b = _slots['_map_b']
        self._error_types = _slots['_error_types']
        Map._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: map1 = sage.structure.coerce_maps.CallableConvertMap(ZZ, QQ, lambda parent, x: 1/x)
            sage: map2 = QQ.coerce_map_from(ZZ)
            sage: map = sage.structure.coerce_maps.TryMap(map1, map2, error_types=(ZeroDivisionError,))
            sage: map(3)
            1/3
            sage: map(-7)
            -1/7
            sage: map(0)
            0
        """
        try:
            return self._map_p._call_(x)
        except self._error_types:
            return self._map_b._call_(x)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        EXAMPLES::

            sage: map1 = sage.structure.coerce_maps.CallableConvertMap(ZZ, QQ, lambda parent, x, y:  y/x)
            sage: map2 = sage.structure.coerce_maps.CallableConvertMap(ZZ, QQ, lambda parent, x, y: 23/1)
            sage: map = sage.structure.coerce_maps.TryMap(map1, map2, error_types=(ZeroDivisionError,))
            sage: map._call_with_args(3, (2,))
            2/3
            sage: map._call_with_args(-7, (5,))
            -5/7
            sage: map._call_with_args(0, (1,))
            23
        """
        try:
            return self._map_p._call_with_args(x, args, kwds)
        except self._error_types:
            return self._map_b._call_with_args(x, args, kwds)
