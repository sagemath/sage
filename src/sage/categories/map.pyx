# sage_setup: distribution = sagemath-objects
r"""
Base class for maps

AUTHORS:

- Robert Bradshaw: initial implementation

- Sebastien Besnier (2014-05-5): :class:`FormalCompositeMap` contains
  a list of Map instead of only two Map. See :issue:`16291`.

- Sebastian Oehms   (2019-01-19): :meth:`section` added to :class:`FormalCompositeMap`.
  See :issue:`27081`.
"""
# ****************************************************************************
#       Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories import homset
import weakref
from sage.ext.stdsage cimport HAS_DICTIONARY
from sage.arith.power cimport generic_power
from sage.sets.pythonclass cimport Set_PythonType
from sage.misc.constant_function import ConstantFunction
from sage.structure.element cimport parent
from cpython.object cimport PyObject_RichCompare


def unpickle_map(_class, parent, _dict, _slots):
    """
    Auxiliary function for unpickling a map.

    TESTS::

        sage: R.<x,y> = QQ[]
        sage: f = R.hom([x+y, x-y], R)
        sage: f == loads(dumps(f))  # indirect doctest
        True
    """
    # should we use slots?
    # from element.pyx
    cdef Map mor = _class.__new__(_class)
    mor._set_parent(parent)
    mor._update_slots(_slots)
    if HAS_DICTIONARY(mor):
        mor.__dict__ = _dict
    return mor


def is_Map(x):
    """
    Auxiliary function: Is the argument a map?

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: f = R.hom([x+y, x-y], R)
        sage: from sage.categories.map import is_Map
        sage: is_Map(f)
        doctest:warning...
        DeprecationWarning: The function is_Map is deprecated; use 'isinstance(..., Map)' instead.
        See https://github.com/sagemath/sage/issues/38103 for details.
        True
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38103, "The function is_Map is deprecated; use 'isinstance(..., Map)' instead.")
    return isinstance(x, Map)


cdef class Map(Element):
    """
    Basic class for all maps.

    .. NOTE::

        The call method is of course not implemented in this base class. This must
        be done in the sub classes, by overloading ``_call_`` and possibly also
        ``_call_with_args``.

    EXAMPLES:

    Usually, instances of this class will not be constructed directly, but
    for example like this::

        sage: from sage.categories.morphism import SetMorphism
        sage: X.<x> = ZZ[]
        sage: Y = ZZ
        sage: phi = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
        sage: phi(x^2+2*x-1)
        -1
        sage: R.<x,y> = QQ[]
        sage: f = R.hom([x+y, x-y], R)
        sage: f(x^2+2*x-1)
        x^2 + 2*x*y + y^2 + 2*x + 2*y - 1
    """

    def __init__(self, parent, codomain=None):
        """
        INPUT:

        There can be one or two arguments of this init method. If it is one argument,
        it must be a hom space. If it is two arguments, it must be two parent structures
        that will be domain and codomain of the map-to-be-created.

        TESTS::

            sage: from sage.categories.map import Map

        Using a hom space::

            sage: Map(Hom(QQ, ZZ, Rings()))
            Generic map:
              From: Rational Field
              To:   Integer Ring

        Using domain and codomain::

            sage: Map(QQ['x'], SymmetricGroup(6))                                       # needs sage.groups
            Generic map:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Symmetric group of order 6! as a permutation group
        """
        if codomain is not None:
            if isinstance(parent, type):
                parent = Set_PythonType(parent)
            parent = homset.Hom(parent, codomain)
        elif not isinstance(parent, homset.Homset):
            raise TypeError("parent (=%s) must be a Homspace" % parent)
        Element.__init__(self, parent)
        D = parent.domain()
        C = parent.codomain()
        self._category_for = parent.homset_category()
        self._codomain = C
        self.domain = ConstantFunction(D)
        self.codomain = ConstantFunction(C)
        self._is_coercion = False
        if D.is_exact() and C.is_exact():
            self._coerce_cost = 10  # default value.
        else:
            self._coerce_cost = 10000  # inexact morphisms are bad.

    def __copy__(self):
        """
        Return copy, with strong references to domain and codomain.

        .. NOTE::

            To implement copying on sub-classes, do not override this method, but
            implement cdef methods ``_extra_slots()`` returning a dictionary and
            ``_update_slots()`` using this dictionary to fill the cdef or cpdef
            slots of the subclass.

        EXAMPLES::

            sage: phi = QQ['x']._internal_coerce_map_from(ZZ)
            sage: phi.domain
            <weakref at ...; to 'sage.rings.integer_ring.IntegerRing_class' at ...>
            sage: type(phi)
            <class 'sage.categories.map.FormalCompositeMap'>
            sage: psi = copy(phi)   # indirect doctest
            sage: psi
            Composite map:
              From: Integer Ring
              To:   Univariate Polynomial Ring in x over Rational Field
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
                    then
                      Polynomial base injection morphism:
                      From: Rational Field
                      To:   Univariate Polynomial Ring in x over Rational Field
            sage: psi.domain
            The constant function (...) -> Integer Ring
            sage: psi(3)
            3
        """
        cdef Map out = Element.__copy__(self)
        # Element.__copy__ updates the __dict__, but not the slots.
        # Let's do this now, but with strong references.
        out._parent = self.parent()  # self._parent might be None
        out._update_slots(self._extra_slots())
        return out

    def parent(self):
        r"""
        Return the homset containing this map.

        .. NOTE::

            The method :meth:`_make_weak_references`, that is used for the maps
            found by the coercion system, needs to remove the usual strong
            reference from the coercion map to the homset containing it. As long
            as the user keeps strong references to domain and codomain of the map,
            we will be able to reconstruct the homset. However, a strong reference
            to the coercion map does not prevent the domain from garbage collection!

        EXAMPLES::

            sage: Q = QuadraticField(-5)                                                # needs sage.rings.number_field
            sage: phi = CDF._internal_convert_map_from(Q)                               # needs sage.rings.number_field
            sage: print(phi.parent())                                                   # needs sage.rings.number_field
            Set of field embeddings
             from Number Field in a with defining polynomial x^2 + 5
                  with a = 2.236067977499790?*I
               to Complex Double Field

        We now demonstrate that the reference to the coercion map `\phi` does
        not prevent `Q` from being garbage collected::

            sage: import gc
            sage: del Q                                                                 # needs sage.rings.number_field
            sage: _ = gc.collect()
            sage: phi.parent()                                                          # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            ValueError: This map is in an invalid state,
            the domain has been garbage collected

        You can still obtain copies of the maps used by the coercion system with
        strong references::

            sage: # needs sage.rings.number_field
            sage: Q = QuadraticField(-5)
            sage: phi = CDF.convert_map_from(Q)
            sage: print(phi.parent())
            Set of field embeddings
             from Number Field in a with defining polynomial x^2 + 5
                  with a = 2.236067977499790?*I
               to Complex Double Field
            sage: import gc
            sage: del Q
            sage: _ = gc.collect()
            sage: phi.parent()
            Set of field embeddings
             from Number Field in a with defining polynomial x^2 + 5
                  with a = 2.236067977499790?*I
               to Complex Double Field
        """
        if self._parent is None:
            D = self.domain()
            C = self._codomain
            if C is None or D is None:
                raise ValueError("This map is in an invalid state, the domain has been garbage collected")
            return homset.Hom(D, C, self._category_for)
        return self._parent

    def _make_weak_references(self):
        """
        Only store weak references to domain and codomain of this map.

        .. NOTE::

            This method is internally used on maps that are used for coercions
            or conversions between parents. Without using this method, some objects
            would stay alive indefinitely as soon as they are involved in a coercion
            or conversion.

        .. SEEALSO::

            :meth:`_make_strong_references`

        EXAMPLES::

            sage: Q = QuadraticField(-5)                                                # needs sage.rings.number_field
            sage: phi = CDF._internal_convert_map_from(Q)                               # needs sage.rings.number_field

        By :issue:`14711`, maps used in the coercion and conversion system
        use *weak* references to domain and codomain, in contrast to other
        maps::

            sage: phi.domain                                                            # needs sage.rings.number_field
            <weakref at ...; to 'NumberField_quadratic_with_category' at ...>
            sage: phi._make_strong_references()                                         # needs sage.rings.number_field
            sage: print(phi.domain)                                                     # needs sage.rings.number_field
            The constant function (...) -> Number Field in a
            with defining polynomial x^2 + 5 with a = 2.236067977499790?*I

        Now, as there is a strong reference, `Q` cannot be garbage collected::

            sage: # needs sage.rings.number_field
            sage: import gc
            sage: _ = gc.collect()
            sage: C = Q.__class__.__base__
            sage: x = None
            sage: numberQuadFields = len([x for x in gc.get_objects()
            ....:                         if isinstance(x, C)])
            sage: del Q, x
            sage: _ = gc.collect()
            sage: numberQuadFields == len([x for x in gc.get_objects()
            ....:                          if isinstance(x, C)])
            True

        However, if we now make the references weak again, the number field can
        be garbage collected, which of course makes the map and its parents
        invalid. This is why :meth:`_make_weak_references` should only be used
        if one really knows what one is doing::

            sage: # needs sage.rings.number_field
            sage: phi._make_weak_references()
            sage: _ = gc.collect()
            sage: numberQuadFields == len([x for x in gc.get_objects()
            ....:                          if isinstance(x, C)]) + 1
            True
            sage: phi
            Defunct map
        """
        if not isinstance(self.domain, ConstantFunction):
            return
        self.domain = weakref.ref(self.domain())
        # Save the category before clearing the parent.
        self._category_for = self._parent.homset_category()
        self._parent = None

    def _make_strong_references(self):
        """
        Store strong references to domain and codomain of this map.

        .. NOTE::

            By default, maps keep strong references to domain and codomain,
            preventing them thus from garbage collection. However, in Sage's
            coercion system, these strong references are replaced by weak
            references, since otherwise some objects would stay alive indefinitely
            as soon as they are involved in a coercion or conversion.

        .. SEEALSO::

            :meth:`_make_weak_references`

        EXAMPLES::

            sage: Q = QuadraticField(-5)                                                # needs sage.rings.number_field
            sage: phi = CDF._internal_convert_map_from(Q)                               # needs sage.rings.number_field

        By :issue:`14711`, maps used in the coercion and conversion system
        use *weak* references to domain and codomain, in contrast to other
        maps::

            sage: phi.domain                                                            # needs sage.rings.number_field
            <weakref at ...; to 'NumberField_quadratic_with_category' at ...>
            sage: phi._make_strong_references()                                         # needs sage.rings.number_field
            sage: print(phi.domain)                                                     # needs sage.rings.number_field
            The constant function (...) -> Number Field in a
            with defining polynomial x^2 + 5 with a = 2.236067977499790?*I

        Now, as there is a strong reference, `Q` cannot be garbage collected::

            sage: # needs sage.rings.number_field
            sage: import gc
            sage: _ = gc.collect()
            sage: C = Q.__class__.__base__
            sage: x = None
            sage: numberQuadFields = len([x for x in gc.get_objects()
            ....:                         if isinstance(x, C)])
            sage: del Q, x
            sage: _ = gc.collect()
            sage: numberQuadFields == len([x for x in gc.get_objects()
            ....:                          if isinstance(x, C)])
            True

        However, if we now make the references weak again, the number field can
        be garbage collected, which of course makes the map and its parents
        invalid. This is why :meth:`_make_weak_references` should only be used
        if one really knows what one is doing::

            sage: # needs sage.rings.number_field
            sage: phi._make_weak_references()
            sage: _ = gc.collect()
            sage: numberQuadFields == len([x for x in gc.get_objects()
            ....:                          if isinstance(x, C)]) + 1
            True
            sage: phi
            Defunct map
            sage: phi._make_strong_references()
            Traceback (most recent call last):
            ...
            RuntimeError: The domain of this map became garbage collected
            sage: phi.parent()
            Traceback (most recent call last):
            ...
            ValueError: This map is in an invalid state, the domain has been garbage collected
        """
        if isinstance(self.domain, ConstantFunction):
            return
        D = self.domain()
        C = self._codomain
        if D is None or C is None:
            raise RuntimeError("The domain of this map became garbage collected")
        self.domain = ConstantFunction(D)
        self._parent = homset.Hom(D, C, self._category_for)

    cdef _update_slots(self, dict slots):
        """
        Set various attributes of this map to implement unpickling.

        INPUT:

        - ``slots`` -- dictionary of slots to be updated;
          the dictionary must have the keys ``'_domain'`` and
          ``'_codomain'``, and may have the keys ``'_repr_type_str'``
          and ``'_is_coercion'``

        TESTS:

        Since it is a ``cdef``d method, it is tested using a dummy python method.
        ::

            sage: # needs sage.rings.real_mpfr
            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._update_slots_test({"_domain": RR, "_codomain": QQ}) # indirect doctest
            sage: f.domain()
            Real Field with 53 bits of precision
            sage: f.codomain()
            Rational Field
            sage: f._repr_type_str
            sage: f._update_slots_test({"_repr_type_str": "bla", "_domain": RR, "_codomain": QQ})
            sage: f._repr_type_str
            'bla'
        """
        # todo: the following can break during unpickling of complex
        # objects with circular references! In that case, _slots might
        # contain incomplete objects.
        self.domain = ConstantFunction(slots['_domain'])
        self._codomain = slots['_codomain']
        self.codomain = ConstantFunction(self._codomain)

        # Several pickles exist without the following, so these are
        # optional
        self._repr_type_str = slots.get('_repr_type_str')
        self._is_coercion = slots.get('_is_coercion')

    def _update_slots_test(self, _slots):
        """
        A Python method to test the cdef _update_slots method.

        TESTS::

            sage: # needs sage.rings.real_mpfr
            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._update_slots_test({"_domain": RR, "_codomain": QQ})
            sage: f.domain()
            Real Field with 53 bits of precision
            sage: f.codomain()
            Rational Field
            sage: f._repr_type_str
            sage: f._update_slots_test({"_repr_type_str": "bla", "_domain": RR, "_codomain": QQ})
            sage: f._repr_type_str
            'bla'
        """
        self._update_slots(_slots)

    cdef dict _extra_slots(self):
        """
        Return a dict with attributes to pickle and copy this map.
        """
        return dict(
            _domain=self.domain(),
            _codomain=self._codomain,
            _is_coercion=self._is_coercion,
            _repr_type_str=self._repr_type_str)

    def _extra_slots_test(self):
        """
        A Python method to test the cdef _extra_slots method.

        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._extra_slots_test()
            {'_codomain': Integer Ring,
             '_domain': Rational Field,
             '_is_coercion': False,
             '_repr_type_str': None}
        """
        return self._extra_slots()

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings())); f
            Generic map:
              From: Rational Field
              To:   Integer Ring
            sage: loads(dumps(f))  # indirect doctest
            Generic map:
              From: Rational Field
              To:   Integer Ring
        """
        if HAS_DICTIONARY(self):
            _dict = self.__dict__
        else:
            _dict = {}
        return unpickle_map, (type(self), self.parent(), _dict, self._extra_slots())

    def _repr_type(self):
        """
        Return a string describing the specific type of this map, to be used when printing ``self``.

        .. NOTE::

            By default, the string ``'Generic'`` is returned. Subclasses may overload this method.

        EXAMPLES::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: print(f._repr_type())
            Generic
            sage: R.<x,y> = QQ[]
            sage: phi = R.hom([x+y, x-y], R)
            sage: print(phi._repr_type())
            Ring
        """
        if self._repr_type_str is None:
            return "Generic"
        else:
            return self._repr_type_str

    def _repr_defn(self):
        """
        Return a string describing the definition of ``self``, to be used when printing ``self``.

        .. NOTE::

            By default, the empty string is returned. Subclasses may overload this method.

        EXAMPLES::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._repr_defn() == ''
            True
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: print(f._repr_defn())
            x |--> x + y
            y |--> x - y
        """
        return ""

    def _repr_(self):
        """
        .. NOTE::

            The string representation is based on the strings returned by
            :meth:`_repr_defn` and :meth:`_repr_type`, as well as the domain
            and the codomain.

            A map that has been subject to :meth:`_make_weak_references` has
            probably been used internally in the coercion system. Hence, it
            may become defunct by garbage collection of the domain. In this
            case, a warning is printed accordingly.

        EXAMPLES::

            sage: from sage.categories.map import Map
            sage: Map(Hom(QQ, ZZ, Rings()))    # indirect doctest
            Generic map:
              From: Rational Field
              To:   Integer Ring
            sage: R.<x,y> = QQ[]
            sage: R.hom([x+y, x-y], R)
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x + y
                    y |--> x - y

        TESTS::

            sage: # needs sage.rings.number_field
            sage: Q = QuadraticField(-5)
            sage: phi = CDF._internal_coerce_map_from(Q); phi
            (map internal to coercion system -- copy before use)
            Composite map:
              From: Number Field in a with defining polynomial x^2 + 5 with a = 2.236067977499790?*I
              To:   Complex Double Field
            sage: del Q
            sage: import gc
            sage: _ = gc.collect()
            sage: phi
            Defunct map
        """
        D = self.domain()
        if D is None:
            return "Defunct map"
        s = "%s map:" % self._repr_type()
        s += "\n  From: %s" % D
        s += "\n  To:   %s" % self._codomain
        if isinstance(self.domain, ConstantFunction):
            d = self._repr_defn()
            if d != '':
                s += "\n  Defn: %s" % ('\n        '.join(d.split('\n')))
        else:
            d = "(map internal to coercion system -- copy before use)"
            s = d + "\n" + s
        return s

    def _default_repr_(self):
        D = self.domain()
        if D is None:
            return "Defunct map"
        s = "%s map:" % self._repr_type()
        s += "\n  From: %s" % D
        s += "\n  To:   %s" % self._codomain
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s" % ('\n        '.join(d.split('\n')))
        return s

    def domains(self):
        """
        Iterate over the domains of the factors of a (composite) map.

        This default implementation simply yields the domain of this map.

        .. SEEALSO:: :meth:`FormalCompositeMap.domains`

        EXAMPLES::

            sage: list(QQ.coerce_map_from(ZZ).domains())
            [Integer Ring]
        """
        yield self.domain()

    def category_for(self):
        """
        Return the category ``self`` is a morphism for.

        .. NOTE::

            This is different from the category of maps to which this
            map belongs *as an object*.

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: phi = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi.category_for()
            Category of rings
            sage: phi.category()
            Category of homsets of unital magmas and additive unital additive magmas
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: f.category_for()
            Join of Category of unique factorization domains
             and Category of algebras with basis over
              (number fields and quotient fields and metric spaces)
             and Category of commutative algebras over
              (number fields and quotient fields and metric spaces)
             and Category of infinite sets
            sage: f.category()
            Category of endsets of unital magmas
             and right modules over (number fields and quotient fields and metric spaces)
             and left modules over (number fields and quotient fields and metric spaces)


        FIXME: find a better name for this method
        """
        if self._category_for is None:
            # This can happen if the map is the result of unpickling.
            # We have initialised self._parent, but could not set
            # self._category_for at that moment, because it could
            # happen that the parent was not fully constructed and
            # did not know its category yet.
            self._category_for = self._parent.homset_category()
        return self._category_for

    def __call__(self, x, *args, **kwds):
        """
        Apply this map to ``x``.

        IMPLEMENTATION:

        - To implement the call method in a subclass of Map, implement
          :meth:`_call_` and possibly also :meth:`_call_with_args` and
          :meth:`pushforward`.
        - If the parent of ``x`` cannot be coerced into the domain of
          ``self``, then the method ``pushforward`` is called with ``x``
          and the other given arguments, provided it is implemented.
          In that way, ``self`` could be applied to objects like ideals
          or sub-modules.
        - If there is no coercion and if ``pushforward`` is not implemented
          or fails, ``_call_`` is  called after conversion into the domain
          (which may be possible even when there is no coercion); if there
          are additional arguments (or keyword arguments),
          :meth:`_call_with_args` is called instead. Note that the
          positional arguments after ``x`` are passed as a tuple to
          :meth:`_call_with_args` and the named arguments are passed
          as a dictionary.

        INPUT:

        - ``x`` -- an element coercible to the domain of ``self``; also objects
          like ideals are supported in some cases

        OUTPUT:

        an element (or ideal, etc.)

        EXAMPLES::

            sage: R.<x,y> = QQ[]; phi = R.hom([y, x])
            sage: phi(y)          # indirect doctest
            x

        We take the image of an ideal::

            sage: I = ideal(x, y); I
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: phi(I)
            Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field

        TESTS:

        We test that the map can be applied to something that converts
        (but not coerces) into the domain and can *not* be dealt with
        by :meth:`pushforward` (see :issue:`10496`)::

            sage: D = {(0, 2): -1, (0, 0): -1, (1, 1): 7, (2, 0): 1/3}
            sage: phi(D)
            -x^2 + 7*x*y + 1/3*y^2 - 1

        We test what happens if the argument can't be converted into
        the domain::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(ZZ, QQ, Rings()))
            sage: f(1/2)
            Traceback (most recent call last):
            ...
            TypeError: 1/2 fails to convert into the map's domain Integer Ring,
            but a `pushforward` method is not properly implemented

        We test that the default call method really works as described
        above (that was fixed in :issue:`10496`)::

            sage: class FOO(Map):
            ....:   def _call_(self, x):
            ....:       print("_call_ {}".format(parent(x)))
            ....:       return self.codomain()(x)
            ....:   def _call_with_args(self, x, args=(), kwds={}):
            ....:       print("_call_with_args {}".format(parent(x)))
            ....:       return self.codomain()(x)^kwds.get('exponent', 1)
            ....:   def pushforward(self, x, exponent=1):
            ....:       print("pushforward {}".format(parent(x)))
            ....:       return self.codomain()(1/x)^exponent
            sage: f = FOO(ZZ, QQ)
            sage: f(1/1)   #indirect doctest
            pushforward Rational Field
            1

        ``_call_`` and ``_call_with_args_`` are used *after* coercion::

            sage: f(int(1))
            _call_ Integer Ring
            1
            sage: f(int(2), exponent=2)
            _call_with_args Integer Ring
            4

        ``pushforward`` is called without conversion::

            sage: f(1/2)
            pushforward Rational Field
            2
            sage: f(1/2, exponent=2)
            pushforward Rational Field
            4

        If the argument does not coerce into the domain, and if
        ``pushforward`` fails, ``_call_`` is tried after conversion::

            sage: g = FOO(QQ, ZZ)
            sage: g(SR(3))                                                              # needs sage.symbolic
            pushforward Symbolic Ring
            _call_ Rational Field
            3
            sage: g(SR(3), exponent=2)                                                  # needs sage.symbolic
            pushforward Symbolic Ring
            _call_with_args Rational Field
            9

        If conversion fails as well, an error is raised::

            sage: h = FOO(ZZ, ZZ)
            sage: h(2/3)
            Traceback (most recent call last):
            ...
            TypeError: 2/3 fails to convert into the map's domain Integer Ring,
            but a `pushforward` method is not properly implemented
        """
        P = parent(x)
        cdef Parent D = self.domain()
        if P is D:  # we certainly want to call _call_/with_args
            if not args and not kwds:
                return self._call_(x)
            return self._call_with_args(x, args, kwds)
        # Is there coercion?
        converter = D._internal_coerce_map_from(P)
        if converter is None:
            try:
                return self.pushforward(x, *args, **kwds)
            except (AttributeError, TypeError, NotImplementedError):
                pass
            try:
                x = D(x)
            except (TypeError, NotImplementedError):
                raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented" % (x, D))
        else:
            x = converter(x)
        if not args and not kwds:
            return self._call_(x)
        return self._call_with_args(x, args, kwds)

    cpdef Element _call_(self, x):
        """
        Call method with a single argument, not implemented in the base class.

        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f(1/2)             # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: <class 'sage.categories.map.Map'>
        """
        raise NotImplementedError(type(self))

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        Call method with multiple arguments, not implemented in the base class.

        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f(1/2, 2, foo='bar')      # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: _call_with_args not overridden to accept arguments for <class 'sage.categories.map.Map'>
        """
        if not args and not kwds:
            return self(x)
        raise NotImplementedError("_call_with_args not overridden "
                                  f"to accept arguments for {type(self)}")

    def __mul__(self, right):
        r"""
        The multiplication * operator is operator composition.

        IMPLEMENTATION:

        If you want to change the behaviour of composition for
        derived classes, please overload :meth:`_composition_`
        (but not :meth:`_composition`!) of the left factor.

        INPUT:

        - ``self`` -- Map
        - ``right`` -- Map

        OUTPUT:

        The map `x \mapsto self(right(x))`.

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: Z = QQ
            sage: phi_xy = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi_yz = SetMorphism(Hom(Y, Z, CommutativeAdditiveMonoids()), lambda y: QQ(y)/2)
            sage: phi_yz * phi_xy
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn:   Generic morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Integer Ring
                    then
                      Generic morphism:
                      From: Integer Ring
                      To:   Rational Field

        If ``right`` is a ring homomorphism given by the images of
        generators, then it is attempted to form the composition
        accordingly. Only if this fails, or if the result does not
        belong to the given homset, a formal composite map is
        returned (as above).
        ::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: f = R.hom([a+b, a-b])
            sage: g = S.hom([x+y, x-y])
            sage: f*g
            Ring endomorphism of Multivariate Polynomial Ring in a, b over Rational Field
              Defn: a |--> 2*a
                    b |--> 2*b
            sage: h = SetMorphism(Hom(S, QQ, Rings()), lambda p: p.lc())
            sage: h*f
            Composite map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Rational Field
              Defn:   Ring morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Multivariate Polynomial Ring in a, b over Rational Field
                      Defn: x |--> a + b
                            y |--> a - b
                    then
                      Generic morphism:
                      From: Multivariate Polynomial Ring in a, b over Rational Field
                      To:   Rational Field
        """
        if not isinstance(right, Map):
            raise TypeError("right (=%s) must be a map to multiply it by %s" % (right, self))
        if right.codomain() != self.domain():
            raise TypeError("self (=%s) domain must equal right (=%s) codomain" % (self, right))
        return self._composition(right)

    def _composition(self, right):
        """
        Composition of maps, which generically returns a :class:`CompositeMap`.

        INPUT:

        - ``self`` -- a Map in some ``Hom(Y, Z, category_left)``
        - ``right`` -- a Map in some ``Hom(X, Y, category_right)``

        OUTPUT:

        Returns the composition of ``self`` and ``right`` as a
        morphism in ``Hom(X, Z, category)`` where ``category`` is the
        meet of ``category_left`` and ``category_right``.

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: Z = QQ
            sage: phi_xy = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi_yz = SetMorphism(Hom(Y, Z, CommutativeAdditiveMonoids()), lambda y: QQ(y)/2)
            sage: phi_yz._composition(phi_xy)
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn:   Generic morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Integer Ring
                    then
                      Generic morphism:
                      From: Integer Ring
                      To:   Rational Field
            sage: phi_yz.category_for()
            Category of commutative additive monoids
        """
        category = self.category_for()._meet_(right.category_for())
        H = homset.Hom(right.domain(), self._codomain, category)
        return self._composition_(right, H)

    def _composition_(self, right, homset):
        """
        INPUT:

        - ``self``, ``right`` -- maps
        - ``homset`` -- a homset

        ASSUMPTION:

        The codomain of ``right`` is contained in the domain of ``self``.
        This assumption is not verified.

        OUTPUT:

        Returns a formal composite map, the composition of ``right``
        followed by ``self``, as a morphism in ``homset``.

        Classes deriving from :class:`Map` are encouraged to override
        this whenever meaningful. This is the case, e.g., for ring
        homomorphisms.

        EXAMPLES::

            sage: Rx.<x> = ZZ['x']
            sage: Ry.<y> = ZZ['y']
            sage: Rz.<z> = ZZ['z']
            sage: phi_xy = Rx.hom([y+1]); phi_xy
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in y over Integer Ring
              Defn: x |--> y + 1
            sage: phi_yz = Ry.hom([z+1]); phi_yz
            Ring morphism:
              From: Univariate Polynomial Ring in y over Integer Ring
              To:   Univariate Polynomial Ring in z over Integer Ring
              Defn: y |--> z + 1
            sage: phi_xz = phi_yz._composition_(phi_xy, Hom(Rx, Rz, Monoids()))
            sage: phi_xz
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in z over Integer Ring
              Defn:   Ring morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Univariate Polynomial Ring in y over Integer Ring
                      Defn: x |--> y + 1
                    then
                      Ring morphism:
                      From: Univariate Polynomial Ring in y over Integer Ring
                      To:   Univariate Polynomial Ring in z over Integer Ring
                      Defn: y |--> z + 1
            sage: phi_xz.category_for()
            Category of monoids

        TESTS:

        This illustrates that it is not tested whether the maps can actually
        be composed, i.e., whether codomain and domain match.
        ::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f_R = R.hom([x+y, x-y], R)
            sage: f_S = S.hom([a+b, a-b], S)
            sage: foo_bar = f_R._composition_(f_S, Hom(S, R, Monoids()))
            sage: foo_bar(a)
            2*x

        However, it is tested when attempting to compose the maps in
        the usual multiplicative notation::

            sage: f_R*f_S
            Traceback (most recent call last):
            ...
            TypeError: self (=Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x + y
                    y |--> x - y) domain must equal right (=Ring endomorphism of Multivariate Polynomial Ring in a, b over Rational Field
              Defn: a |--> a + b
                    b |--> a - b) codomain
        """
        return FormalCompositeMap(homset, right, self)

    def pre_compose(self, right):
        """
        INPUT:

        - ``self`` -- a Map in some ``Hom(Y, Z, category_left)``
        - ``left`` -- a Map in some ``Hom(X, Y, category_right)``

        Returns the composition of ``right`` followed by ``self`` as a
        morphism in ``Hom(X, Z, category)`` where ``category`` is the
        meet of ``category_left`` and ``category_right``.

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: Z = QQ
            sage: phi_xy = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi_yz = SetMorphism(Hom(Y, Z, Monoids()), lambda y: QQ(y**2))
            sage: phi_xz = phi_yz.pre_compose(phi_xy); phi_xz
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn:   Generic morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Integer Ring
                    then
                      Generic morphism:
                      From: Integer Ring
                      To:   Rational Field
            sage: phi_xz.category_for()
            Category of monoids
        """
        D = self.domain()
        if D is not right.codomain():
            right = right.extend_codomain(D)
        return self._composition(right)

    def post_compose(self, left):
        """
        INPUT:

        - ``self`` -- a Map in some ``Hom(X, Y, category_right)``
        - ``left`` -- a Map in some ``Hom(Y, Z, category_left)``

        Returns the composition of ``self`` followed by ``left`` as a
        morphism in ``Hom(X, Z, category)`` where ``category`` is the
        meet of ``category_left`` and ``category_right``.

        Caveat: see the current restrictions on :meth:`Category.meet`

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: Z = QQ
            sage: phi_xy = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi_yz = SetMorphism(Hom(Y, Z, Monoids()), lambda y: QQ(y**2))
            sage: phi_xz = phi_xy.post_compose(phi_yz); phi_xz
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn:   Generic morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Integer Ring
                    then
                      Generic morphism:
                      From: Integer Ring
                      To:   Rational Field
            sage: phi_xz.category_for()
            Category of monoids
        """
        return left._composition(self)

    def extend_domain(self, new_domain):
        r"""
        INPUT:

        - ``self`` -- a member of Hom(Y, Z)
        - ``new_codomain`` -- an object X such that there is a canonical coercion
          `\phi` in Hom(X, Y)

        OUTPUT:

        An element of Hom(X, Z) obtained by composing self with `\phi`.  If
        no canonical `\phi` exists, a :exc:`TypeError` is raised.

        EXAMPLES::

            sage: # needs sage.rings.complex_double
            sage: mor = CDF.coerce_map_from(RDF)
            sage: mor.extend_domain(QQ)
            Composite map:
              From: Rational Field
              To:   Complex Double Field
              Defn:   Native morphism:
                      From: Rational Field
                      To:   Real Double Field
                    then
                      Native morphism:
                      From: Real Double Field
                      To:   Complex Double Field
            sage: mor.extend_domain(ZZ['x'])
            Traceback (most recent call last):
            ...
            TypeError: No coercion from Univariate Polynomial Ring in x over Integer Ring
            to Real Double Field
        """
        D = self.domain()
        if D is None:
            raise ValueError("This map became defunct by garbage collection")
        cdef Map connecting = D._internal_coerce_map_from(new_domain)
        if connecting is None:
            raise TypeError("No coercion from %s to %s" % (new_domain, D))
        elif connecting.codomain() is not D:
            raise RuntimeError("BUG: coerce_map_from should always return a map to self (%s)" % D)
        else:
            return self.pre_compose(connecting.__copy__())

    def extend_codomain(self, new_codomain):
        r"""
        INPUT:

        - ``self`` -- a member of Hom(X, Y)
        - ``new_codomain`` -- an object Z such that there is a canonical coercion
          `\phi` in Hom(Y, Z)

        OUTPUT:

        An element of Hom(X, Z) obtained by composing ``self`` with `\phi`.  If
        no canonical `\phi` exists, a :exc:`TypeError` is raised.

        EXAMPLES::

            sage: mor = QQ.coerce_map_from(ZZ)
            sage: mor.extend_codomain(RDF)
            Composite map:
              From: Integer Ring
              To:   Real Double Field
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
                    then
                      Native morphism:
                      From: Rational Field
                      To:   Real Double Field
            sage: mor.extend_codomain(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: No coercion from Rational Field to Finite Field of size 7
        """
        cdef Map connecting = new_codomain._internal_coerce_map_from(self._codomain)
        if connecting is None:
            raise TypeError("No coercion from %s to %s" % (self._codomain, new_codomain))
        elif connecting.domain() is not self._codomain:
            raise RuntimeError("BUG: coerce_map_from should always return a map from its input (%s)" % new_codomain)
        else:
            return self.post_compose(connecting.__copy__())

    def is_surjective(self):
        """
        Tell whether the map is surjective (not implemented in the base class).

        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f.is_surjective()
            Traceback (most recent call last):
            ...
            NotImplementedError: <class 'sage.categories.map.Map'>
        """
        raise NotImplementedError(type(self))

    cpdef _pow_int(self, n):
        """
        TESTS::

            sage: R.<x> = ZZ['x']
            sage: phi = R.hom([x+1]); phi
            Ring endomorphism of Univariate Polynomial Ring in x over Integer Ring
              Defn: x |--> x + 1

            sage: phi^0
            Identity endomorphism of Univariate Polynomial Ring in x over Integer Ring

            sage: phi^2 == phi*phi
            True

            sage: S.<y> = QQ[]
            sage: psi = R.hom([y^2])
            sage: psi^1
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in y over Rational Field
              Defn: x |--> y^2
            sage: psi^2
            Traceback (most recent call last):
            ...
            TypeError: self must be an endomorphism

            sage: # needs sage.rings.number_field
            sage: K.<a> = NumberField(x^4 - 5*x + 5)
            sage: C5.<z> = CyclotomicField(5)
            sage: tau = K.hom([z - z^2]); tau
            Ring morphism:
              From: Number Field in a with defining polynomial x^4 - 5*x + 5
              To:   Cyclotomic Field of order 5 and degree 4
              Defn: a |--> -z^2 + z
            sage: tau^-1
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Number Field in a with defining polynomial x^4 - 5*x + 5
              Defn: z |--> 3/11*a^3 + 4/11*a^2 + 9/11*a - 14/11
        """
        if self.domain() is not self._codomain and n != 1 and n != -1:
            raise TypeError("self must be an endomorphism")
        if n == 0:
            from sage.categories.morphism import IdentityMorphism
            return IdentityMorphism(self._parent)
        return generic_power(self, n)

    def section(self):
        """
        Return a section of ``self``.

        .. NOTE::

            By default, it returns ``None``. You may override it in subclasses.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: print(f.section())
            None

            sage: f = QQ.coerce_map_from(ZZ); f
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: ff = f.section(); ff
            Generic map:
              From: Rational Field
              To:   Integer Ring
            sage: ff(4/2)
            2
            sage: parent(ff(4/2)) is ZZ
            True
            sage: ff(1/2)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
        """
        return None

    def __hash__(self):
        """
        Return the hash of this map.

        TESTS::

            sage: f = sage.rings.morphism.RingMap(ZZ.Hom(ZZ))
            sage: type(f)
            <class 'sage.rings.morphism.RingMap'>
            sage: hash(f) == hash(f)
            True
            sage: {f: 1}[f]
            1
        """
        D = self.domain()
        if D is None:
            raise ValueError("This map became defunct by garbage collection")
        return hash((self.domain(), self._codomain))

cdef class Section(Map):
    """
    A formal section of a map.

    .. NOTE::

        Call methods are not implemented for the base class ``Section``.

    EXAMPLES::

        sage: from sage.categories.map import Section
        sage: R.<x,y> = ZZ[]
        sage: S.<a,b> = QQ[]
        sage: f = R.hom([a+b, a-b])
        sage: sf = Section(f); sf
        Section map:
          From: Multivariate Polynomial Ring in a, b over Rational Field
          To:   Multivariate Polynomial Ring in x, y over Integer Ring
        sage: sf(a)
        Traceback (most recent call last):
        ...
        NotImplementedError: <class 'sage.categories.map.Section'>
    """

    def __init__(self, map):
        """
        INPUT:

        - ``map`` -- a map

        TESTS::

            sage: from sage.categories.map import Section
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: sf = Section(f); sf
            Section map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
        """
        from sage.categories.homset import Hom
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        Map.__init__(self, Hom(map.codomain(), map.domain(), SetsWithPartialMaps()))
        self._inverse = map    # TODO: Use this attribute somewhere!

    cdef dict _extra_slots(self):
        """
        Helper for pickling and copying.

        TESTS::

            sage: from sage.categories.map import Section
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: sf = Section(f)
            sage: copy(sf)   # indirect doctest
            Section map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
        """
        slots = Map._extra_slots(self)
        slots['_inverse'] = self._inverse
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for pickling and copying.

        TESTS::

            sage: from sage.categories.map import Section
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: sf = Section(f)
            sage: copy(sf)   # indirect doctest
            Section map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
        """
        Map._update_slots(self, _slots)
        self._inverse = _slots['_inverse']

    def _repr_type(self):
        """
        Return a string describing the type of this map (which is "Section").

        TESTS::

            sage: from sage.categories.map import Section
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: sf = Section(f)
            sage: sf         # indirect doctest
            Section map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
        """
        return "Section"

    def inverse(self):
        """
        Return inverse of ``self``.

        TESTS::

            sage: from sage.categories.map import Section
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y, x-y], R)
            sage: sf = Section(f)
            sage: sf.inverse()
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x + y
                    y |--> x - y
        """
        return self._inverse


cdef class FormalCompositeMap(Map):
    """
    Formal composite maps.

    A formal composite map is formed by two maps, so that the codomain of the
    first map is contained in the domain of the second map.

    .. NOTE::

        When calling a composite with additional arguments, these arguments are
        *only* passed to the second underlying map.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: S.<a> = QQ[]
        sage: from sage.categories.morphism import SetMorphism
        sage: f = SetMorphism(Hom(R, S, Rings()), lambda p: p[0]*a^p.degree())
        sage: g = S.hom([2*x])
        sage: f*g
        Composite map:
          From: Univariate Polynomial Ring in a over Rational Field
          To:   Univariate Polynomial Ring in a over Rational Field
          Defn:   Ring morphism:
                  From: Univariate Polynomial Ring in a over Rational Field
                  To:   Univariate Polynomial Ring in x over Rational Field
                  Defn: a |--> 2*x
                then
                  Generic morphism:
                  From: Univariate Polynomial Ring in x over Rational Field
                  To:   Univariate Polynomial Ring in a over Rational Field
        sage: g*f
        Composite map:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Univariate Polynomial Ring in x over Rational Field
          Defn:   Generic morphism:
                  From: Univariate Polynomial Ring in x over Rational Field
                  To:   Univariate Polynomial Ring in a over Rational Field
                then
                  Ring morphism:
                  From: Univariate Polynomial Ring in a over Rational Field
                  To:   Univariate Polynomial Ring in x over Rational Field
                  Defn: a |--> 2*x
        sage: (f*g)(2*a^2+5)
        5*a^2
        sage: (g*f)(2*x^2+5)
        20*x^2
    """

    def __init__(self, parent, first, second=None):
        """
        INPUT:

        - ``parent`` -- a homset
        - ``first`` -- a map or a list of maps
        - ``second`` -- a map or None

        .. NOTE::

            The intended use is of course that the codomain of the
            first map is contained in the domain of the second map,
            so that the two maps can be composed, and that the
            composition belongs to ``parent``. However, none of
            these conditions is verified in the init method.

            The user is advised to compose two maps ``f`` and ``g``
            in multiplicative notation, ``g*f``, since this will in
            some cases return a more efficient map object than a
            formal composite map.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b, a-b])
            sage: g = S.hom([x+y, x-y])
            sage: H = Hom(R, R, Rings())
            sage: from sage.categories.map import FormalCompositeMap
            sage: m = FormalCompositeMap(H, f, g); m
            Composite map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn:   Ring morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Multivariate Polynomial Ring in a, b over Rational Field
                      Defn: x |--> a + b
                            y |--> a - b
                    then
                      Ring morphism:
                      From: Multivariate Polynomial Ring in a, b over Rational Field
                      To:   Multivariate Polynomial Ring in x, y over Rational Field
                      Defn: a |--> x + y
                            b |--> x - y
            sage: m(x), m(y)
            (2*x, 2*y)
        """
        Map.__init__(self, parent)

        if isinstance(first, (list, tuple)):
            self.__list = first
            self._coerce_cost = sum((<Map>f)._coerce_cost for f in first)
            return

        self.__list = []
        if isinstance(first, FormalCompositeMap):
            self.__list += (<FormalCompositeMap>first).__list
        else:
            self.__list += [first]

        if isinstance(second, FormalCompositeMap):
            self.__list += (<FormalCompositeMap>second).__list
        else:
            self.__list += [second]
        self._coerce_cost = (<Map>first)._coerce_cost + (<Map>second)._coerce_cost

    def __copy__(self):
        """
        Since :meth:`_extra_slots` would return the uncopied constituents
        of this composite map, we cannot rely on the default copying method
        of maps.

        TESTS::

            sage: copy(QQ['q,t'].coerce_map_from(int))   # indirect doctest
            Composite map:
              From: Set of Python objects of class 'int'
              To:   Multivariate Polynomial Ring in q, t over Rational Field
              Defn:   Native morphism:
                      From: Set of Python objects of class 'int'
                      To:   Rational Field
                    then
                      Polynomial base injection morphism:
                      From: Rational Field
                      To:   Multivariate Polynomial Ring in q, t over Rational Field
        """
        return FormalCompositeMap(self.parent(), [f.__copy__() for f in self.__list])

    cdef _update_slots(self, dict _slots):
        """
        Used in pickling and copying.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b, a-b])
            sage: g = S.hom([x+y, x-y])
            sage: from sage.categories.map import FormalCompositeMap
            sage: H = Hom(R, R, Rings())
            sage: m = FormalCompositeMap(H, f, g)
            sage: m == loads(dumps(m))    # indirect doctest
            True
        """
        self.__list = _slots['__list']
        Map._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Used in pickling and copying.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b, a-b])
            sage: g = S.hom([x+y, x-y])
            sage: from sage.categories.map import FormalCompositeMap
            sage: H = Hom(R, R, Rings())
            sage: m = FormalCompositeMap(H, f, g)
            sage: m == loads(dumps(m))    # indirect doctest
            True
        """
        slots = Map._extra_slots(self)
        slots['__list'] = self.__list
        return slots

    def __richcmp__(self, other, int op):
        """
        TESTS::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b, a-b])
            sage: g = S.hom([x+y, x-y])
            sage: from sage.categories.map import FormalCompositeMap
            sage: H = Hom(R, R, Rings())
            sage: m = FormalCompositeMap(H, f, g)
            sage: m == loads(dumps(m))
            True

            sage: m == None
            False
            sage: m == 2
            False
        """
        if type(self) is not type(other):
            return NotImplemented
        left = (<FormalCompositeMap>self).__list
        right = (<FormalCompositeMap>other).__list
        return PyObject_RichCompare(left, right, op)

    def __hash__(self):
        """
        Return the hash of this map.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b, a-b])
            sage: g = S.hom([x+y, x-y])
            sage: from sage.categories.map import FormalCompositeMap
            sage: H = Hom(R, R, Rings())
            sage: m = FormalCompositeMap(H, f, g)
            sage: hash(m) == hash(m)
            True
            sage: {m: 1}[m]
            1
            sage: n = FormalCompositeMap(Hom(S, S, Rings()), g, f)
            sage: hash(m) == hash(n)
            False
            sage: len({m: 1, n: 2}.keys())
            2
        """
        return hash(tuple(self.__list))

    def __getitem__(self, i):
        r"""
        Return the `i`-th map of the formal composition.

        If ``self`` represents `f_n \circ f_{n-1} \circ \cdots \circ
        f_1 \circ f_0`, then ``self[i]`` gives `f_i`.  Support
        negative indices as ``list.__getitem__``.  Raise an error if
        the index does not match, in the same way as
        ``list.__getitem__``.

        EXAMPLES::

            sage: from sage.categories.map import Map
            sage: f = Map(ZZ, QQ)
            sage: g = Map(QQ, ZZ)
            sage: (f*g)[0]
            Generic map:
              From: Rational Field
              To:   Integer Ring
            sage: (f*g)[1]
            Generic map:
              From: Integer Ring
              To:   Rational Field
            sage: (f*g)[-1]
            Generic map:
              From: Integer Ring
              To:   Rational Field
            sage: (f*g)[-2]
            Generic map:
              From: Rational Field
              To:   Integer Ring
            sage: (f*g)[-3]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: (f*g)[2]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self.__list[i]

    cpdef Element _call_(self, x):
        """
        Call with a single argument.

        TESTS::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R, S, Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: (g*f)((x+1)^2), (f*g)((a+1)^2)     # indirect doctest
            (4*x^2, a^2)
        """
        for f in self.__list:
            x = f._call_(x)
        return x

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        Additional arguments are only passed to the last applied map.

        TESTS::

            sage: from sage.categories.morphism import SetMorphism
            sage: R.<x> = QQ[]
            sage: def foo(x, *args, **kwds):
            ....:     print('foo called with {} {}'.format(args, kwds))
            ....:     return x
            sage: def bar(x, *args, **kwds):
            ....:     print('bar called with {} {}'.format(args, kwds))
            ....:     return x
            sage: f = SetMorphism(Hom(R, R, Rings()), foo)
            sage: b = SetMorphism(Hom(R, R, Rings()), bar)
            sage: c = b*f
            sage: c(2, 'hello world', test=1)     # indirect doctest
            foo called with () {}
            bar called with ('hello world',) {'test': 1}
            2
            sage: c = f*b
            sage: c(2, 'hello world', test=1)
            bar called with () {}
            foo called with ('hello world',) {'test': 1}
            2
        """
        for f in self.__list[:-1]:
            x = f._call_(x)
        return self.__list[-1]._call_with_args(x, args, kwds)

    def _repr_type(self):
        """
        Return a string describing the type of ``self``, namely "Composite".

        TESTS::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R, S, Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: f*g            # indirect doctest
            Composite map:
              From: Univariate Polynomial Ring in a over Rational Field
              To:   Univariate Polynomial Ring in a over Rational Field
              Defn:   Ring morphism:
                      From: Univariate Polynomial Ring in a over Rational Field
                      To:   Univariate Polynomial Ring in x over Rational Field
                      Defn: a |--> 2*x
                    then
                      Generic morphism:
                      From: Univariate Polynomial Ring in x over Rational Field
                      To:   Univariate Polynomial Ring in a over Rational Field
        """
        return "Composite"

    def _repr_defn(self):
        """
        Return a string describing the definition of ``self``.

        The return value is obtained from the string representations
        of the two constituents.

        TESTS::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R, S, Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: f*g            # indirect doctest
            Composite map:
              From: Univariate Polynomial Ring in a over Rational Field
              To:   Univariate Polynomial Ring in a over Rational Field
              Defn:   Ring morphism:
                      From: Univariate Polynomial Ring in a over Rational Field
                      To:   Univariate Polynomial Ring in x over Rational Field
                      Defn: a |--> 2*x
                    then
                      Generic morphism:
                      From: Univariate Polynomial Ring in x over Rational Field
                      To:   Univariate Polynomial Ring in a over Rational Field
        """
        s = "  %s" % (self.__list[0])
        for f in self.__list[1:]:
            s += "\nthen\n  %s" % f
        return s

    def first(self):
        r"""
        Return the first map in the formal composition.

        If ``self`` represents `f_n \circ f_{n-1} \circ \cdots \circ
        f_1 \circ f_0`, then ``self.first()`` returns `f_0`.  We have
        ``self == self.then() * self.first()``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R, S, Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: fg = f * g
            sage: fg.first() == g
            True
            sage: fg == fg.then() * fg.first()
            True
        """
        return self.__list[0]

    def then(self):
        r"""
        Return the tail of the list of maps.

        If ``self`` represents `f_n \circ f_{n-1} \circ \cdots \circ
        f_1 \circ f_0`, then ``self.first()`` returns `f_n \circ
        f_{n-1} \circ \cdots \circ f_1`.  We have ``self ==
        self.then() * self.first()``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R, S, Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: (f*g).then() == f
            True

            sage: f = QQ.coerce_map_from(ZZ)
            sage: f = f.extend_domain(ZZ).extend_codomain(QQ)
            sage: f.then()
            Composite map:
            From: Integer Ring
            To:   Rational Field
            Defn:   Natural morphism:
            From: Integer Ring
            To:   Rational Field
            then
            Identity endomorphism of Rational Field
        """
        if len(self.__list) == 2:
            return self.__list[1]
        domain = self.__list[0].codomain()
        codomain = self.codomain()
        H = homset.Hom(domain, codomain, category=self._category_for)
        return FormalCompositeMap(H, self.__list[1:])

    def is_injective(self):
        """
        Tell whether ``self`` is injective.

        It raises :exc:`NotImplementedError` if it cannot be determined.

        EXAMPLES::

            sage: # needs sage.modules
            sage: V1 = QQ^2
            sage: V2 = QQ^3
            sage: phi1 = (QQ^1).hom(Matrix([[1, 1]]), V1)
            sage: phi2 = V1.hom(Matrix([[1, 2, 3], [4, 5, 6]]), V2)

        If both constituents are injective, the composition is injective::

            sage: from sage.categories.map import FormalCompositeMap
            sage: c1 = FormalCompositeMap(Hom(QQ^1, V2, phi1.category_for()),           # needs sage.modules
            ....:                         phi1, phi2)
            sage: c1.is_injective()                                                     # needs sage.modules
            True

        If it cannot be determined whether the composition is injective,
        an error is raised::

            sage: psi1 = V2.hom(Matrix([[1, 2], [3, 4], [5, 6]]), V1)                   # needs sage.modules
            sage: c2 = FormalCompositeMap(Hom(V1, V1, phi2.category_for()),             # needs sage.modules
            ....:                         phi2, psi1)
            sage: c2.is_injective()                                                     # needs sage.modules
            Traceback (most recent call last):
            ...
            NotImplementedError: not enough information to deduce injectivity

        If the first map is surjective and the second map is not injective,
        then the composition is not injective::

            sage: psi2 = V1.hom([[1], [1]], QQ^1)                                       # needs sage.modules
            sage: c3 = FormalCompositeMap(Hom(V2, QQ^1, phi2.category_for()),           # needs sage.modules
            ....:                         psi2, psi1)
            sage: c3.is_injective()                                                     # needs sage.modules
            False

        TESTS:

        Check that :issue:`23205` has been resolved::

            sage: f = QQ.hom(QQbar) * ZZ.hom(QQ)                                        # needs sage.rings.number_field
            sage: f.is_injective()                                                      # needs sage.rings.number_field
            True
        """
        try:
            # we try the category first
            # as of 2017-06, the MRO of this class does not get patched to
            # include the category's MorphismMethods (because it is a Cython
            # class); therefore, we cannot simply call "super" but need to
            # invoke the category method explicitly
            return self.getattr_from_category('is_injective')()
        except (AttributeError, NotImplementedError):
            pass

        injectives = []
        for f in self.__list:
            if f.is_injective():
                injectives.append(f)
            else:
                break
        else:
            return True

        if all(f.is_surjective() for f in injectives):
            return False

        raise NotImplementedError("not enough information to deduce injectivity")

    def is_surjective(self):
        """
        Tell whether ``self`` is surjective.

        It raises :exc:`NotImplementedError` if it cannot be determined.

        EXAMPLES::

            sage: from sage.categories.map import FormalCompositeMap
            sage: V3 = QQ^3                                                             # needs sage.modules
            sage: V2 = QQ^2                                                             # needs sage.modules
            sage: V1 = QQ^1                                                             # needs sage.modules

        If both maps are surjective, the composition is surjective::

            sage: # needs sage.modules
            sage: phi32 = V3.hom(Matrix([[1, 2], [3, 4], [5, 6]]), V2)
            sage: phi21 = V2.hom(Matrix([[1], [1]]), V1)
            sage: c_phi = FormalCompositeMap(Hom(V3, V1, phi32.category_for()),
            ....:                            phi32, phi21)
            sage: c_phi.is_surjective()
            True

        If the second map is not surjective, the composition is not
        surjective::

            sage: FormalCompositeMap(Hom(V3, V1, phi32.category_for()),                 # needs sage.modules
            ....:                    phi32,
            ....:                    V2.hom(Matrix([[0], [0]]), V1)).is_surjective()
            False

        If the second map is an isomorphism and the first map is not
        surjective, then the composition is not surjective::

            sage: FormalCompositeMap(Hom(V2, V1, phi32.category_for()),                 # needs sage.modules
            ....:                    V2.hom(Matrix([[0], [0]]), V1),
            ....:                    V1.hom(Matrix([[1]]), V1)).is_surjective()
            False

        Otherwise, surjectivity of the composition cannot be determined::

            sage: FormalCompositeMap(Hom(V2, V1, phi32.category_for()),                 # needs sage.modules
            ....:     V2.hom(Matrix([[1, 1], [1, 1]]), V2),
            ....:     V2.hom(Matrix([[1], [1]]), V1)).is_surjective()
            Traceback (most recent call last):
            ...
            NotImplementedError: not enough information to deduce surjectivity
        """
        try:
            # we try the category first
            # as of 2017-06, the MRO of this class does not get patched to
            # include the category's MorphismMethods (because it is a Cython
            # class); therefore, we cannot simply call "super" but need to
            # invoke the category method explicitly
            return self.getattr_from_category('is_surjective')()
        except (AttributeError, NotImplementedError):
            pass

        surjectives = []
        for f in self.__list[::-1]:
            if f.is_surjective():
                surjectives.append(f)
            else:
                break
        else:
            return True

        if all(f.is_injective() for f in surjectives):
            return False

        raise NotImplementedError("not enough information to deduce surjectivity")

    def domains(self):
        """
        Iterate over the domains of the factors of this map.

        (This is useful in particular to check for loops in coercion maps.)

        .. SEEALSO:: :meth:`Map.domains`

        EXAMPLES::

            sage: f = QQ.coerce_map_from(ZZ)
            sage: g = MatrixSpace(QQ, 2, 2).coerce_map_from(QQ)                         # needs sage.modules
            sage: list((g * f).domains())                                               # needs sage.modules
            [Integer Ring, Rational Field]
        """
        for f in self.__list:
            yield f.domain()

    def section(self):
        """
        Compute a section map from sections of the factors of
        ``self`` if they have been implemented.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: incl = P.coerce_map_from(ZZ)
            sage: sect = incl.section(); sect
            Composite map:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Integer Ring
              Defn:   Generic map:
                      From: Univariate Polynomial Ring in x over Rational Field
                      To:   Rational Field
                    then
                      Generic map:
                      From: Rational Field
                      To:   Integer Ring
            sage: p = x + 5; q = x + 2
            sage: sect(p-q)
            3

        the following example has been attached to :meth:`_integer_`
        of :class:`sage.rings.polynomial.polynomial_element.Polynomial`
        before (see comment there)::

            sage: k = GF(47)
            sage: R.<x> = PolynomialRing(k)
            sage: R.coerce_map_from(ZZ).section()
            Composite map:
              From: Univariate Polynomial Ring in x over Finite Field of size 47
              To:   Integer Ring
              Defn:   Generic map:
                      From: Univariate Polynomial Ring in x over Finite Field of size 47
                      To:   Finite Field of size 47
                    then
                      Lifting map:
                      From: Finite Field of size 47
                      To:   Integer Ring
            sage: ZZ(R(45))                 # indirect doctest
            45
            sage: ZZ(3*x + 45)              # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: 3*x + 45 is not a constant polynomial
        """
        sections = []
        for m in reversed(list(self)):
            try:
                sec = m.section()
            except TypeError:
                return None
            if sec is None:
                return None
            sections.append(sec)

        from sage.categories.homset import Hom
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        H = Hom(self.codomain(), self.domain(), category=SetsWithPartialMaps())
        return FormalCompositeMap(H, sections)
