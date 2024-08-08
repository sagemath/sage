# sage_setup: distribution = sagemath-objects
r"""
Base class for old-style parent objects

CLASS HIERARCHY:

SageObject
    Parent
        ParentWithBase
            ParentWithGens


TESTS:

This came up in some subtle bug once::

    sage: gp(2) + gap(3)                                                                # needs sage.libs.gap sage.libs.pari
    5
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.superseded import deprecation
from sage.structure.coerce cimport py_scalar_parent
from sage.ext.stdsage cimport HAS_DICTIONARY
from sage.sets.pythonclass cimport Set_PythonType, Set_PythonType_class

from cpython.object cimport *
from cpython.bool cimport *

cdef inline check_old_coerce(Parent p):
    if p._element_constructor is not None:
        raise RuntimeError("%s still using old coercion framework" % p)


cdef class Parent(parent.Parent):
    """
    Parents are the Sage / mathematical analogues of container objects
    in computer science.

    TESTS::

        sage: # needs sage.modules
        sage: V = VectorSpace(GF(2,'a'), 2)
        sage: V.list()
        [(0, 0), (1, 0), (0, 1), (1, 1)]
        sage: MatrixSpace(GF(3), 1, 1).list()
        [[0], [1], [2]]
        sage: DirichletGroup(3).list()                                                  # needs sage.libs.pari sage.modular
        [Dirichlet character modulo 3 of conductor 1 mapping 2 |--> 1,
         Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1]

        sage: # needs sage.rings.finite_rings
        sage: K = GF(7^6,'a')
        sage: K.list()[:10]                     # long time
        [0, 1, 2, 3, 4, 5, 6, a, a + 1, a + 2]
        sage: K.<a> = GF(4)
        sage: K.list()
        [0, a, a + 1, 1]
    """

    def __cinit__(self):
        self._has_coerce_map_from = MonoDict()

    def __init__(self, *, category=None):
        self.init_coerce(False)
        self._coerce_from_hash = MonoDict()
        self._set_element_constructor()
        if category is not None:
            self._init_category_(category)

    ##########################################################
    # New Coercion support functionality
    ##########################################################

    cdef __coerce_map_from_c(self, S):
        """
        EXAMPLES:

        Check to make sure that we handle coerce maps from Python
        native types correctly::

            sage: QQ['q,t'].coerce_map_from(int)
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
        check_old_coerce(self)
        if S is self:
            from sage.categories.homset import Hom
            return Hom(self, self).identity()
        elif S == self:
            # non-unique parents
            from sage.categories.homset import Hom
            from sage.categories.morphism import CallMorphism
            return CallMorphism(Hom(S, self))
        elif isinstance(S, Set_PythonType_class):
            return self.__coerce_map_from_c(S._type)
        if self._coerce_from_hash is None: # this is because parent.__init__() does not always get called
            self.init_coerce()
        cdef object ret
        try:
            ret = self._coerce_from_hash.get(S)
            return ret
        except KeyError:
            pass

        mor = self.__coerce_map_from_c_impl(S)
        import sage.categories.morphism
        import sage.categories.map
        if mor is True:
            mor = sage.categories.morphism.CallMorphism(S, self)
        elif mor is False:
            mor = None
        elif mor is not None and not isinstance(mor, sage.categories.map.Map):
            raise AssertionError("__coerce_map_from_c_impl must return a boolean, None, or an explicit Map")

        if mor is None and isinstance(S, type):
            #Convert Python types to native Sage types
            sage_type = py_scalar_parent(S)
            if sage_type is None:
                self._coerce_from_hash[S] = None
                return None
            mor = self.__coerce_map_from_c(sage_type)
            if mor is not None:
                mor = mor * sage_type._internal_coerce_map_from(S)

        if mor is not None:
            self._coerce_from_hash.set(S, mor) # TODO: if this is None, could it be non-None in the future?

        return mor

    cdef __coerce_map_from_c_impl(self, S):
        check_old_coerce(self)
        import sage.categories.morphism
        from sage.categories.map import Map
        from sage.categories.homset import Hom
        cdef parent.Parent R
        for mor in self._coerce_from_list:
            if isinstance(mor, Map):
                R = mor.domain()
            else:
                R = mor
                mor = sage.categories.morphism.CallMorphism(Hom(R, self))
                i = self._coerce_from_list.index(R)
                self._coerce_from_list[i] = mor # cache in case we need it again
            if R is S:
                return mor
            else:
                connecting = R._internal_coerce_map_from(S)
                if connecting is not None:
                    return mor * connecting

        if self.__has_coerce_map_from_c(S):
            if isinstance(S, type):
                S = Set_PythonType(S)
            return sage.categories.morphism.CallMorphism(Hom(S, self))
        else:
            return None

    ##############################################
    # Coercion support functionality
    ##############################################

    def _coerce_(self, x):            # Call this from Python (do not override!)
        if self._element_constructor is not None:
            from sage.misc.superseded import deprecation
            deprecation(33497, "_coerce_ is deprecated, use coerce instead")
            return self.coerce(x)
        check_old_coerce(self)
        return self._coerce_c(x)

    cpdef _coerce_c(self, x):          # DO NOT OVERRIDE THIS (call it)
        if self._element_constructor is not None:
            from sage.misc.superseded import deprecation
            deprecation(33497, "_coerce_c is deprecated, use coerce instead")
            return self.coerce(x)
        check_old_coerce(self)
        try:
            P = x.parent()   # todo -- optimize
            if P is self:
                return x
        except AttributeError as msg:
            pass
        if HAS_DICTIONARY(self):
            return self._coerce_impl(x)
        else:
            return self._coerce_c_impl(x)

    cdef _coerce_c_impl(self, x):     # OVERRIDE THIS FOR CYTHON CLASSES
        """
        Canonically coerce ``x`` in assuming that the parent of ``x`` is not
        equal to ``self``.
        """
        check_old_coerce(self)
        raise TypeError

    def _coerce_impl(self, x):        # OVERRIDE THIS FOR PYTHON CLASSES
        """
        Canonically coerce ``x`` in assuming that the parent of ``x`` is not
        equal to ``self``.
        """
        check_old_coerce(self)
        return self._coerce_c_impl(x)

    cdef __has_coerce_map_from_c(self, S):
        check_old_coerce(self)
        if self == S:
            return True
        try:
            return self._has_coerce_map_from.get(S)
        except KeyError:
            pass
        try:
            self._coerce_c((<parent.Parent?>S).an_element())
        except TypeError:
            ans = False
        else:
            ans = True
        self._has_coerce_map_from.set(S, ans)
        return ans

    def _an_element_impl(self):     # override this in Python
        """
        Return an element of ``self``.

        Want it in sufficient generality
        that poorly-written functions will not work when they are not
        supposed to. This is cached so does not have to be super fast.
        """
        check_old_coerce(self)
        try:
            return self.gen(0)
        except Exception:
            pass

        try:
            return self.gen()
        except Exception:
            pass

        from sage.rings.infinity import infinity
        for x in ['_an_element_', 'pi', 1.2, 2, 1, 0, infinity]:
            try:
                return self(x)
            except Exception:
                pass

        raise NotImplementedError(f"_an_element_ is not implemented for {self}")

    ###############################################################
    # Coercion Compatibility Layer
    ###############################################################
    cpdef _coerce_map_from_(self, S):
        if self._element_constructor is None:
            return self.__coerce_map_from_c(S)
        else:
            return parent.Parent._coerce_map_from_(self, S)

    def _an_element_(self):
        if self._element_constructor is not None:
            return parent.Parent._an_element_(self)
        if self._cache_an_element is not None:
            return self._cache_an_element
        self._cache_an_element = self._an_element_impl()
        return self._cache_an_element

    cpdef _generic_convert_map(self, S, category=None):
        r"""
        Return a default conversion from ``S``.

        EXAMPLES::

           sage: R.<x,y> = QQ[]
           sage: R._generic_convert_map(QQ).category_for()
           Category of sets with partial maps
        """
        if self._element_constructor is None:
            if hasattr(self, '_element_constructor_'):
                assert callable(self._element_constructor_)
                self._element_constructor = self._element_constructor_
            else:
                from sage.categories.morphism import CallMorphism
                from sage.categories.homset import Hom
                if isinstance(S, type):
                    S = Set_PythonType(S)
                return CallMorphism(Hom(S, self))
        return parent.Parent._generic_convert_map(self, S, category)
