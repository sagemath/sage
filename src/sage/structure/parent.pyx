r"""
Base class for parent objects

CLASS HIERARCHY::

    SageObject
        CategoryObject
            Parent

A simple example of registering coercions::

    sage: class A_class(Parent):
    ....:   def __init__(self, name):
    ....:       Parent.__init__(self, name=name)
    ....:       self._populate_coercion_lists_()
    ....:       self.rename(name)
    ....:   #
    ....:   def category(self):
    ....:       return Sets()
    ....:   #
    ....:   def _element_constructor_(self, i):
    ....:       assert(isinstance(i, (int, Integer)))
    ....:       return ElementWrapper(self, i)
    ....:
    sage: A = A_class("A")
    sage: B = A_class("B")
    sage: C = A_class("C")

    sage: def f(a):
    ....:   return B(a.value+1)
    ....:
    sage: class MyMorphism(Morphism):
    ....:   def __init__(self, domain, codomain):
    ....:       Morphism.__init__(self, Hom(domain, codomain))
    ....:   #
    ....:   def _call_(self, x):
    ....:       return self.codomain()(x.value)
    ....:
    sage: f = MyMorphism(A,B)
    sage: f
        Generic morphism:
          From: A
          To:   B
    sage: B.register_coercion(f)
    sage: C.register_coercion(MyMorphism(B,C))
    sage: A(A(1)) == A(1)
    True
    sage: B(A(1)) == B(1)
    True
    sage: C(A(1)) == C(1)
    True

    sage: A(B(1))
    Traceback (most recent call last):
    ...
    AssertionError

When implementing an element of a ring, one would typically provide the
element class with ``_rmul_`` and/or ``_lmul_`` methods for the action of a
base ring, and with ``_mul_`` for the ring multiplication. However, prior to
:trac:`14249`, it would have been necessary to additionally define a method
``_an_element_()`` for the parent. But now, the following example works::

    sage: from sage.structure.element import RingElement
    sage: class MyElement(RingElement):
    ....:      def __init__(self, parent, x, y):
    ....:          RingElement.__init__(self, parent)
    ....:      def _mul_(self, other):
    ....:          return self
    ....:      def _rmul_(self, other):
    ....:          return self
    ....:      def _lmul_(self, other):
    ....:          return self
    sage: class MyParent(Parent):
    ....:      Element = MyElement

Now, we define ::

    sage: P = MyParent(base=ZZ, category=Rings())
    sage: a = P(1,2)
    sage: a*a is a
    True
    sage: a*2 is a
    True
    sage: 2*a is a
    True

TESTS:

This came up in some subtle bug once::

    sage: gp(2) + gap(3)
    5

"""

cimport element
cimport sage.categories.morphism as morphism
cimport sage.categories.map as map
from sage.structure.debug_options import debug
from sage.structure.sage_object import SageObject
from sage.structure.misc import (dir_with_other_class, getattr_from_other_class,
                                 is_extension_type)
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.sets_cat import Sets, EmptySetError
from copy import copy
from sage.misc.sage_itertools import unique_merge
from sage.misc.lazy_format import LazyFormat

# Create a dummy attribute error, using some kind of lazy error message,
# so that neither the error itself not the message need to be created
# repeatedly, which would cost time.

from sage.structure.misc cimport AttributeErrorMessage
cdef AttributeErrorMessage dummy_error_message = AttributeErrorMessage(None, '')
dummy_attribute_error = AttributeError(dummy_error_message)

# TODO: define this once?

cdef object elt_parent = None

cdef inline parent_c(x):
    if PY_TYPE_CHECK(x, element.Element):
        return (<element.Element>x)._parent
#    elif hasattr(x, 'parent'):
#        return x.parent()
#    else:
#        return <object>PY_TYPE(x)
    else:
        try:
            return x.parent()
        except AttributeError:
            return <object>PY_TYPE(x)

cdef _record_exception():
    from element import get_coercion_model
    get_coercion_model()._record_exception()

cdef object _Integer
cdef bint is_Integer(x):
    global _Integer
    if _Integer is None:
        from sage.rings.integer import Integer as _Integer
    return PY_TYPE_CHECK_EXACT(x, _Integer) or PY_TYPE_CHECK_EXACT(x, int)

# for override testing
cdef extern from "descrobject.h":
    ctypedef struct PyMethodDef:
        void *ml_meth
    ctypedef struct PyMethodDescrObject:
        PyMethodDef *d_method
    void* PyCFunction_GET_FUNCTION(object)
    bint PyCFunction_Check(object)

cdef extern from *:
    Py_ssize_t PyDict_Size(object)
    Py_ssize_t PyTuple_GET_SIZE(object)

    ctypedef class __builtin__.dict [object PyDictObject]:
        cdef Py_ssize_t ma_fill
        cdef Py_ssize_t ma_used

    void* PyDict_GetItem(object, object)

cdef inline Py_ssize_t PyDict_GET_SIZE(o):
    return (<dict>o).ma_used

###############################################################################
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2008 Burcin Erocal   <burcin@erocal.org>
#       Copyright (C) 2008 Mike Hansen     <mhansen@gmail.com>
#       Copyright (C) 2008 David Roe       <roed@math.harvard.edu>
#       Copyright (C) 2007 William Stein   <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

import operator
import weakref

from category_object import CategoryObject
from generators import Generators
from coerce_exceptions import CoercionException

# TODO: Theses should probably go in the printer module (but lots of stuff imports them from parent)
from category_object import localvars

cdef object BuiltinMethodType = type(repr)

from cpython.object cimport *
from cpython.bool cimport *
include 'sage/ext/stdsage.pxi'


def is_Parent(x):
    """
    Return True if x is a parent object, i.e., derives from
    sage.structure.parent.Parent and False otherwise.

    EXAMPLES::

        sage: from sage.structure.parent import is_Parent
        sage: is_Parent(2/3)
        False
        sage: is_Parent(ZZ)
        True
        sage: is_Parent(Primes())
        True
    """
    return PY_TYPE_CHECK(x, Parent)

cdef bint guess_pass_parent(parent, element_constructor):
    if PY_TYPE_CHECK(element_constructor, types.MethodType):
        return False
    elif PY_TYPE_CHECK(element_constructor, BuiltinMethodType):
        return element_constructor.__self__ is not parent
    else:
        return True

from sage.categories.category import Category
from sage.structure.dynamic_class import dynamic_class
Sets_parent_class = Sets().parent_class

cdef inline bint good_as_coerce_domain(S):
    """
    Determine whether the input can be the domain of a map.

    NOTE:

    This is the same as being an object in a category, or
    being a type. Namely, in Sage, we do consider coercion maps
    from the type ``<int>`` to, say, `ZZ`.

    TESTS:

    If an instance `S` is not suitable as domain of a map, then
    the non-existence of a coercion or conversion map from `S`
    to some other parent is not cached, by :trac:`13378`::

        sage: P.<x,y> = QQ[]
        sage: P.is_coercion_cached(x)
        False
        sage: P.coerce_map_from(x)
        sage: P.is_coercion_cached(x)  # indirect doctest
        False

    """
    return isinstance(S,CategoryObject) or isinstance(S,type)

cdef inline bint good_as_convert_domain(S):
    return isinstance(S,SageObject) or isinstance(S,type)

cdef class Parent(category_object.CategoryObject):

    def __init__(self, base=None, *, category=None, element_constructor=None,
                 gens=None, names=None, normalize=True, facade=None, **kwds):
        """
        Base class for all parents.

        Parents are the Sage/mathematical analogues of container
        objects in computer science.

        INPUT:

        - ``base`` -- An algebraic structure considered to be the
          "base" of this parent (e.g. the base field for a vector
          space).

        - ``category`` -- a category or list/tuple of categories. The
          category in which this parent lies (or list or tuple
          thereof).  Since categories support more general
          super-categories, this should be the most specific category
          possible. If category is a list or tuple, a JoinCategory is
          created out of them.  If category is not specified, the
          category will be guessed (see
          :class:`~sage.structure.category_object.CategoryObject`),
          but won't be used to inherit parent's or element's code from
          this category.

        - ``element_constructor`` -- A class or function that creates
          elements of this Parent given appropriate input (can also be
          filled in later with :meth:`_populate_coercion_lists_`)

        - ``gens`` -- Generators for this object (can also be filled in
          later with :meth:`_populate_generators_`)

        - ``names`` -- Names of generators.

        - ``normalize`` -- Whether to standardize the names (remove
          punctuation, etc)

        - ``facade`` -- a parent, or tuple thereof, or ``True``

        If ``facade`` is specified, then ``Sets().Facades()`` is added
        to the categories of the parent. Furthermore, if ``facade`` is
        not ``True``, the internal attribute ``_facade_for`` is set
        accordingly for use by
        :meth:`Sets.Facades.ParentMethods.facade_for`.

        Internal invariants:

        - ``self._element_init_pass_parent == guess_pass_parent(self,
          self._element_constructor)`` Ensures that :meth:`__call__`
          passes down the parent properly to
          :meth:`_element_constructor`.  See :trac:`5979`.

        .. TODO::

            Eventually, category should be
            :class:`~sage.categories.sets_cat.Sets` by default.

        .. automethod:: __call__
        .. automethod:: _populate_coercion_lists_
        .. automethod:: __mul__
        .. automethod:: __contains__
        .. automethod:: _coerce_map_from_
        .. automethod:: _convert_map_from_
        .. automethod:: _get_action_
        .. automethod:: _an_element_
        .. automethod:: _repr_option
        .. automethod:: _init_category_
        """
        # TODO: in the long run, we want to get rid of the element_constructor = argument
        # (element_constructor would always be set by inheritance)
        # Then, all the element_constructor logic should be moved to init_coerce.
        CategoryObject.__init__(self, base = base)
        if facade is not None and facade is not False:
            if isinstance(facade, Parent):
                facade = (facade,)
            if facade is not True:
                assert isinstance(facade, tuple)
                self._facade_for = facade
            if category is None:
                category = Sets().Facades()
            else:
                if isinstance(category, (tuple,list)):
                    category = Category.join(tuple(category) + (Sets().Facades(),))
                else:
                    category = Category.join((category,Sets().Facades()))
        # Setting the categories is currently done in a separate
        # method to let some subclasses (like ParentsWithBase)
        # call it without calling the full constructor
        self._init_category_(category)

        if len(kwds) > 0:
            if debug.bad_parent_warnings:
                print "Illegal keywords for %s: %s" % (type(self), kwds)
        # TODO: many classes don't call this at all, but __new__ crashes Sage
        if debug.bad_parent_warnings:
            if element_constructor is not None and not callable(element_constructor):
                print "coerce BUG: Bad element_constructor provided", type(self), type(element_constructor), element_constructor
        if gens is not None:
            self._populate_generators_(gens, names, normalize)
        elif names is not None:
            self._assign_names(names, normalize)
        if element_constructor is None:
            self._set_element_constructor()
        else:
            self._element_constructor = element_constructor
            self._element_init_pass_parent = guess_pass_parent(self, element_constructor)
        self.init_coerce(False)

        for cls in self.__class__.mro():
            # this calls __init_extra__ if it is *defined* in cls (not in a super class)
            if "__init_extra__" in cls.__dict__:
                cls.__init_extra__(self)

    def _init_category_(self, category):
        """
        Initialize the category framework

        Most parents initialize their category upon construction, and
        this is the recommended behavior. For example, this happens
        when the constructor calls :meth:`Parent.__init__` directly or
        indirectly. However, some parents defer this for performance
        reasons. For example,
        :mod:`sage.matrix.matrix_space.MatrixSpace` does not.

        EXAMPLES::

            sage: P = Parent()
            sage: P.category()
            Category of sets
            sage: class MyParent(Parent):
            ....:     def __init__(self):
            ....:         self._init_category_(Groups())
            sage: MyParent().category()
            Category of groups
        """
        CategoryObject._init_category_(self, category)

        # This substitutes the class of this parent to a subclass
        # which also subclasses the parent_class of the category

        if category is not None: #isinstance(self._category, Category) and not isinstance(self, Set_generic):
            category = self._category # CategoryObject may have done some argument processing
            # Some parent class may readily have their category classes attached
            # TODO: assert that the category is consistent
            if not issubclass(self.__class__, Sets_parent_class) and not is_extension_type(self.__class__):
                #documentation transfer is handled by dynamic_class
                self.__class__ = dynamic_class(
                    '{0}_with_category'.format(self.__class__.__name__),
                    (self.__class__, category.parent_class, ),
                    doccls=self.__class__)

    def _refine_category_(self, category):
        """
        Change the category of ``self`` into a subcategory.

        INPUT:

        - ``category`` -- a category or list or tuple thereof

        The new category is obtained by adjoining ``category`` to the
        current one.

        .. NOTE::

            The class of ``self`` might be replaced by a sub-class.

        .. seealso::

            :meth:`CategoryObject._refine_category`

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: Q = P.quotient(x^2+2)
            sage: Q.category()
            Join of Category of commutative rings and Category of subquotients of monoids and Category of quotients of semigroups
            sage: first_class = Q.__class__
            sage: Q._refine_category_(Fields())
            sage: Q.category()
            Join of Category of fields and Category of subquotients of monoids and Category of quotients of semigroups
            sage: first_class == Q.__class__
            False
            sage: TestSuite(Q).run()


        TESTS:

        Here is a test against :trac:`14471`. Refining the category will issue
        a warning, if this change affects the hash value (note that this will
        only be seen in doctest mode)::

            sage: class MyParent(Parent):
            ....:     def __hash__(self):
            ....:         return hash(type(self))   # subtle mistake
            sage: a = MyParent()
            sage: h_a = hash(a)
            sage: a._refine_category_(Algebras(QQ))
            hash of <class '__main__.MyParent_with_category'> changed in
            Parent._refine_category_ during initialisation

            sage: b = MyParent(category=Rings())
            sage: h_b = hash(b)
            sage: h_a == h_b
            False
            sage: b._refine_category_(Algebras(QQ))
            hash of <class '__main__.MyParent_with_category'> changed in
            Parent._refine_category_ during refinement
            sage: hash(a) == hash(b)
            True
            sage: hash(a) != h_a
            True

        """
        if debug.refine_category_hash_check:
            # check that the hash stays the same after refinement
            hash_old = hash(self)
        if self._category is None:
            self._init_category_(category)
            if debug.refine_category_hash_check and hash_old != hash(self):
                print 'hash of {0} changed in Parent._refine_category_ during initialisation' \
                    .format(str(self.__class__))
            return
        if category is self._category:
            return
        CategoryObject._refine_category_(self, category)
        category = self._category

        # This substitutes the class of this parent to a subclass
        # which also subclasses the parent_class of the category.
        # However, we only do so if we don't have an extension class.
        if not is_extension_type(self.__class__):
            # We tested in the very beginning that this parent
            # had its category initialised. Hence, the class
            # is already a dynamic class.
            base = self.__class__.__base__
            #documentation transfer is handled by dynamic_class
            self.__class__     = dynamic_class("%s_with_category"%base.__name__,
                                               (base, category.parent_class, ),
                                               doccls=base)
        # If the element class has already been assigned, it
        # needs to be erased now.
        try:
            self.__dict__.__delitem__('element_class')
        except (AttributeError, KeyError):
            pass
        if debug.refine_category_hash_check and hash_old != hash(self):
            print 'hash of {0} changed in Parent._refine_category_ during refinement' \
                .format(str(self.__class__))


    # This probably should go into Sets().Parent
    @lazy_attribute
    def element_class(self):
        """
        The (default) class for the elements of this parent

        FIXME's and design issues:

        - If self.Element is "trivial enough", should we optimize it away with:
          self.element_class = dynamic_class("%s.element_class"%self.__class__.__name__, (category.element_class,), self.Element)
        - This should lookup for Element classes in all super classes
        """
        try: #if hasattr(self, 'Element'):
            return self.__make_element_class__(self.Element, "%s.element_class"%self.__class__.__name__)
        except AttributeError: #else:
            return NotImplemented


    def __make_element_class__(self, cls, name = None, inherit = None):
        """
        A utility to construct classes for the elements of this
        parent, with appropriate inheritance from the element class of
        the category (only for pure python types so far).
        """
        if name is None:
            name = "%s_with_category"%cls.__name__
        # By default, don't fiddle with extension types yet; inheritance from
        # categories will probably be achieved in a different way
        if inherit is None:
            inherit = not is_extension_type(cls)
        if inherit:
            return dynamic_class(name, (cls, self.category().element_class))
        else:
            return cls

    def _set_element_constructor(self):
        """
        This function is used in translating from the old to the new coercion model.

        It is called from sage.structure.parent_old.Parent.__init__
        when an old style parent provides a _element_constructor_ method.

        It just asserts that this _element_constructor_ is callable and
        also sets self._element_init_pass_parent

        EXAMPLES::

            sage: k = GF(5); k._element_constructor # indirect doctest
            <bound method FiniteField_prime_modn_with_category._element_constructor_ of Finite Field of size 5>
        """
        try: #if hasattr(self, '_element_constructor_'):
            _element_constructor_ = self._element_constructor_
        except (AttributeError, TypeError):
            # Remark: A TypeError can actually occur;
            # it is a possible reason for "hasattr" to return False
            return
        assert callable(_element_constructor_)
        self._element_constructor = _element_constructor_
        self._element_init_pass_parent = guess_pass_parent(self, self._element_constructor)

    def category(self):
        """
        EXAMPLES::

            sage: P = Parent()
            sage: P.category()
            Category of sets
            sage: class MyParent(Parent):
            ....:     def __init__(self): pass
            sage: MyParent().category()
            Category of sets
        """
        if self._category is None:
            # COERCE TODO: we shouldn't need this
            self._category = Sets()
        return self._category

    def _test_category(self, **options):
        """
        Run generic tests on the method :meth:`.category`.

        See also: :class:`TestSuite`.

        EXAMPLES::

            sage: C = Sets().example()
            sage: C._test_category()

        Let us now write a parent with broken categories:

            sage: class MyParent(Parent):
            ....:     def __init__(self):
            ....:         pass
            sage: P = MyParent()
            sage: P._test_category()
            Traceback (most recent call last):
            ...
            AssertionError: category of self improperly initialized

        To fix this, :meth:`MyParent.__init__` should initialize the
        category of ``self`` by calling :meth:`._init_category` or
        ``Parent.__init__(self, category = ...)``.
        """
        tester = self._tester(**options)
        SageObject._test_category(self, tester = tester)
        category = self.category()
        tester.assert_(category.is_subcategory(Sets()))
        # Tests that self inherits methods from the categories
        if not is_extension_type(self.__class__):
            # For usual Python classes, that should be done with
            # standard inheritance
            tester.assertTrue(isinstance(self, category.parent_class),
                LazyFormat("category of self improperly initialized")%self)
        else:
            # For extension types we just check that inheritance
            # occurs on one specific method.
            # _test_an_element from Sets().ParentMethods is a good
            # candidate because it's unlikely to be overriden in self.
            tester.assertTrue(hasattr(self, "_test_an_element"),
                LazyFormat("category of self improperly initialized")%self)

    def _test_eq(self, **options):
        """
        Test that ``self`` is equal to ``self`` and different to ``None``.

        See also: :class:`TestSuite`.

        TESTS::

            sage: O = Parent()
            sage: O._test_eq()

        Let us now write a broken class method::

            sage: class CCls(Parent):
            ....:     def __eq__(self, other):
            ....:         return True
            sage: CCls()._test_eq()
            Traceback (most recent call last):
            ...
            AssertionError: broken equality: <class '__main__.CCls'> == None

        Let us now break inequality::

            sage: class CCls(Parent):
            ....:     def __ne__(self, other):
            ....:         return True
            sage: CCls()._test_eq()
            Traceback (most recent call last):
            ...
            AssertionError: broken non-equality: <class '__main__.CCls'> != itself
        """
        tester = self._tester(**options)

        # We don't use assertEqual / assertNonEqual in order to be
        # 100% sure we indeed call the operators == and !=, whatever
        # the version of Python is (see #11236)
        tester.assertTrue(self == self,
                   LazyFormat("broken equality: %s == itself is False")%self)
        tester.assertFalse(self == None,
                   LazyFormat("broken equality: %s == None")%self)
        tester.assertFalse(self != self,
                   LazyFormat("broken non-equality: %s != itself")%self)
        tester.assertTrue(self != None,
                   LazyFormat("broken non-equality: %s != None is False")%self)

    cdef int init_coerce(self, bint warn=True) except -1:
        if self._coerce_from_hash is None:
            if warn:
                print "init_coerce() for ", type(self)
                raise ZeroDivisionError("hello")
            self._initial_coerce_list = []
            self._initial_action_list = []
            self._initial_convert_list = []
            self._coerce_from_list = []
            self._coerce_from_hash = MonoDict(23)
            self._action_list = []
            self._action_hash = TripleDict(23)
            self._convert_from_list = []
            self._convert_from_hash = MonoDict(53)
            self._embedding = None

    def __getattr__(self, str name):
        """
        Let cat be the category of ``self``. This method emulates
        ``self`` being an instance of both ``Parent`` and
        ``cat.parent_class``, in that order, for attribute lookup.

        NOTE:

        Attributes beginning with two underscores but not ending with
        an unnderscore are considered private and are thus exempted
        from the lookup in ``cat.parent_class``.

        EXAMPLES:

        We test that ZZ (an extension type) inherits the methods from
        its categories, that is from ``EuclideanDomains().parent_class``::

            sage: ZZ._test_associativity
            <bound method EuclideanDomains.parent_class._test_associativity of Integer Ring>
            sage: ZZ._test_associativity(verbose = True)
            sage: TestSuite(ZZ).run(verbose = True)
            running ._test_additive_associativity() . . . pass
            running ._test_an_element() . . . pass
            running ._test_associativity() . . . pass
            running ._test_category() . . . pass
            running ._test_characteristic() . . . pass
            running ._test_distributivity() . . . pass
            running ._test_elements() . . .
              Running the test suite of self.an_element()
              running ._test_category() . . . pass
              running ._test_eq() . . . pass
              running ._test_nonzero_equal() . . . pass
              running ._test_not_implemented_methods() . . . pass
              running ._test_pickling() . . . pass
              pass
            running ._test_elements_eq_reflexive() . . . pass
            running ._test_elements_eq_symmetric() . . . pass
            running ._test_elements_eq_transitive() . . . pass
            running ._test_elements_neq() . . . pass
            running ._test_eq() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_one() . . . pass
            running ._test_pickling() . . . pass
            running ._test_prod() . . . pass
            running ._test_some_elements() . . . pass
            running ._test_zero() . . . pass

            sage: Sets().example().sadfasdf
            Traceback (most recent call last):
            ...
            AttributeError: 'PrimeNumbers_with_category' object has no attribute 'sadfasdf'

        TESTS:

        We test that "private" attributes are not requested from the element class
        of the category (:trac:`10467`)::

            sage: P = Parent(QQ, category=CommutativeRings())
            sage: P.category().parent_class.__foo = 'bar'
            sage: P.__foo
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.structure.parent.Parent' object has no attribute '__foo'

        """
        if (name.startswith('__') and not name.endswith('_')) or self._category is None:
            dummy_error_message.cls = type(self)
            dummy_error_message.name = name
            raise dummy_attribute_error
        try:
            return self.__cached_methods[name]
        except KeyError:
            attr = getattr_from_other_class(self, self._category.parent_class, name)
            self.__cached_methods[name] = attr
            return attr
        except TypeError:
            attr = getattr_from_other_class(self, self._category.parent_class, name)
            self.__cached_methods = {name:attr}
            return attr

    def __dir__(self):
        """
        Let cat be the category of ``self``. This method emulates
        ``self`` being an instance of both ``Parent`` and
        ``cat.parent_class``, in that order, for attribute directory.

        EXAMPLES::

            sage: for s in dir(ZZ):
            ....:     if s[:6] == "_test_": print s
            _test_additive_associativity
            _test_an_element
            _test_associativity
            _test_category
            _test_characteristic
            _test_distributivity
            _test_elements
            _test_elements_eq_reflexive
            _test_elements_eq_symmetric
            _test_elements_eq_transitive
            _test_elements_neq
            _test_eq
            _test_not_implemented_methods
            _test_one
            _test_pickling
            _test_prod
            _test_some_elements
            _test_zero
            sage: F = GF(9,'a')
            sage: dir(F)
            [..., '__class__', ..., '_test_pickling', ..., 'extension', ...]

        """
        return dir_with_other_class(self, self.category().parent_class)

    def _introspect_coerce(self):
        """
        Used for debugging the coercion model.

        EXAMPLES::

            sage: sorted(QQ._introspect_coerce().items())
            [('_action_hash', <sage.structure.coerce_dict.TripleDict object at ...>),
             ('_action_list', []),
             ('_coerce_from_hash', <sage.structure.coerce_dict.MonoDict object at ...>),
             ('_coerce_from_list', []),
             ('_convert_from_hash', <sage.structure.coerce_dict.MonoDict object at ...>),
             ('_convert_from_list', [...]),
             ('_element_init_pass_parent', False),
             ('_embedding', None),
             ('_initial_action_list', []),
             ('_initial_coerce_list', []),
             ('_initial_convert_list', [])]
        """
        return {
            '_coerce_from_list': self._coerce_from_list,
            '_coerce_from_hash': self._coerce_from_hash,
            '_action_list': self._action_list,
            '_action_hash': self._action_hash,
            '_convert_from_list': self._convert_from_list,
            '_convert_from_hash': self._convert_from_hash,
            '_embedding': self._embedding,
            '_initial_coerce_list': self._initial_coerce_list,
            '_initial_action_list': self._initial_action_list,
            '_initial_convert_list': self._initial_convert_list,
            '_element_init_pass_parent': self._element_init_pass_parent,
        }

    def __getstate__(self):
        """
        Used for pickling.

        TESTS::

            sage: loads(dumps(RR['x'])) == RR['x']
            True
        """
        d = CategoryObject.__getstate__(self)
        d['_embedding'] = self._embedding
        d['_element_constructor'] = self._element_constructor
        d['_convert_method_name'] = self._convert_method_name
        d['_element_init_pass_parent'] = self._element_init_pass_parent
        d['_initial_coerce_list'] = self._initial_coerce_list
        d['_initial_action_list'] = self._initial_action_list
        d['_initial_convert_list'] = self._initial_convert_list
        return d

    def __setstate__(self, d):
        """
        Used for pickling.

        TESTS::

            sage: loads(dumps(CDF['x'])) == CDF['x']
            True
        """
        CategoryObject.__setstate__(self, d)
        try:
            version = d['_pickle_version']
        except KeyError:
            version = 0
        if version == 1:
            self.init_coerce(False) # Really, do we want to init this with the same initial data as before?
            self._populate_coercion_lists_(coerce_list=d['_initial_coerce_list'] or [],
                                           action_list=d['_initial_action_list'] or [],
                                           convert_list=d['_initial_convert_list'] or [],
                                           embedding=d['_embedding'],
                                           convert_method_name=d['_convert_method_name'],
                                           element_constructor = d['_element_constructor'],
                                           init_no_parent=not d['_element_init_pass_parent'],
                                           unpickling=True)

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        INPUT:

        - ``key`` -- string. A key for different metadata informations
          that can be inquired about.

        Valid ``key`` arguments are:

        - ``'ascii_art'``: The :meth:`_repr_` output is multi-line
          ascii art and each line must be printed starting at the same
          column, or the meaning is lost.

        - ``'element_ascii_art'``: same but for the output of the
          elements. Used in :mod:`sage.misc.displayhook`.

        - ``'element_is_atomic'``: the elements print atomically, that
          is, parenthesis are not required when *printing* out any of
          `x - y`, `x + y`, `x^y` and `x/y`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: ZZ._repr_option('ascii_art')
            False
            sage: MatrixSpace(ZZ, 2)._repr_option('element_ascii_art')
            True
        """
        if not isinstance(key, basestring):
            raise ValueError('key must be a string')
        defaults = {
            'ascii_art': False,
            'element_ascii_art': False,
            'element_is_atomic': False,
            }
        return defaults[key]

    def is_atomic_repr(self):
        """
        The old way to signal atomic string reps.

        True if the elements have atomic string representations, in the
        sense that if they print at s, then -s means the negative of s. For
        example, integers are atomic but polynomials are not.

        EXAMPLES::

            sage: Parent().is_atomic_repr()
            doctest:...: DeprecationWarning: Use _repr_option to return metadata about string rep
            See http://trac.sagemath.org/14040 for details.
            False
        """
        from sage.misc.superseded import deprecation
        deprecation(14040, 'Use _repr_option to return metadata about string rep')
        return False

    def __call__(self, x=0, *args, **kwds):
        """
        This is the generic call method for all parents.

        When called, it will find a map based on the Parent (or type) of x.
        If a coercion exists, it will always be chosen. This map will
        then be called (with the arguments and keywords if any).

        By default this will dispatch as quickly as possible to
        :meth:`_element_constructor_` though faster pathways are
        possible if so desired.

        TESTS:

        We check that the invariant::

                self._element_init_pass_parent == guess_pass_parent(self, self._element_constructor)

        is preserved (see :trac:`5979`)::

            sage: class MyParent(Parent):
            ....:     def _element_constructor_(self, x):
            ....:         print self, x
            ....:         return sage.structure.element.Element(parent = self)
            ....:     def _repr_(self):
            ....:         return "my_parent"
            ....:
            sage: my_parent = MyParent()
            sage: x = my_parent("bla")
            my_parent bla
            sage: x.parent()         # indirect doctest
            my_parent

            sage: x = my_parent()    # shouldn't this one raise an error?
            my_parent 0
            sage: x = my_parent(3)   # todo: not implemented  why does this one fail???
            my_parent 3


        """
        if self._element_constructor is None:
            # Neither __init__ nor _populate_coercion_lists_ have been called...
            try:
                assert callable(self._element_constructor_)
                self._element_constructor = self._element_constructor_
                self._element_init_pass_parent = guess_pass_parent(self, self._element_constructor)
            except (AttributeError, AssertionError):
                raise NotImplementedError
        cdef Py_ssize_t i
        cdef R = parent_c(x)
        cdef bint no_extra_args = PyTuple_GET_SIZE(args) == 0 and PyDict_GET_SIZE(kwds) == 0
        if R is self and no_extra_args:
            return x

        # Here we inline the first part of convert_map_from for speed.
        # (Yes, the virtual function overhead can matter.)
        if self._convert_from_hash is None: # this is because parent.__init__() does not always get called
            self.init_coerce()
        cdef map.Map mor
        try:
            mor = <map.Map> self._convert_from_hash.get(R)
        except KeyError:
            mor = <map.Map> self.convert_map_from(R)

        if mor is not None:
            if no_extra_args:
                return mor._call_(x)
            else:
                return mor._call_with_args(x, args, kwds)

        raise TypeError("No conversion defined from %s to %s"%(R, self))

    def __mul__(self,x):
        """
        This is a multiplication method that more or less directly
        calls another attribute ``_mul_`` (single underscore). This
        is because ``__mul__`` can not be implemented via inheritance
        from the parent methods of the category, but ``_mul_`` can
        be inherited. This is, e.g., used when creating twosided
        ideals of matrix algebras. See :trac:`7797`.

        EXAMPLE::

            sage: MS = MatrixSpace(QQ,2,2)

        This matrix space is in fact an algebra, and in particular
        it is a ring, from the point of view of categories::

            sage: MS.category()
            Category of algebras over Rational Field
            sage: MS in Rings()
            True

        However, its class does not inherit from the base class
        ``Ring``::

            sage: isinstance(MS,Ring)
            False

        Its ``_mul_`` method is inherited from the category, and
        can be used to create a left or right ideal::

            sage: MS._mul_.__module__
            'sage.categories.rings'
            sage: MS*MS.1      # indirect doctest
            Left Ideal
            (
              [0 1]
              [0 0]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: MS*[MS.1,2]
            Left Ideal
            (
              [0 1]
              [0 0],
            <BLANKLINE>
              [2 0]
              [0 2]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: MS.1*MS
            Right Ideal
            (
              [0 1]
              [0 0]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: [MS.1,2]*MS
            Right Ideal
            (
              [0 1]
              [0 0],
            <BLANKLINE>
              [2 0]
              [0 2]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        """
        # generic multiplication method. It defers to
        # _mul_, which may be defined via categories.
        _mul_ = None
        switch = False
        try:
            if isinstance(self,Parent):
                _mul_ = self._mul_
        except AttributeError:
            pass
        if _mul_ is None:
            try:
                if isinstance(x,Parent):
                    _mul_ = x._mul_
                    switch = True
            except AttributeError:
                pass
        if _mul_ is None:
            raise TypeError("For implementing multiplication, provide the method '_mul_' for %s resp. %s"%(self,x))
        if switch:
            return _mul_(self,switch_sides=True)
        return _mul_(x)

    #############################################################################
    # Containment testing
    #############################################################################
    def __contains__(self, x):
        r"""
        True if there is an element of self that is equal to x under
        ==, or if x is already an element of self.  Also, True in other
        cases involving the Symbolic Ring, which is handled specially.

        For many structures we test this by using :meth:`__call__` and
        then testing equality between x and the result.

        The Symbolic Ring is treated differently because it is
        ultra-permissive about letting other rings coerce in, but
        ultra-strict about doing comparisons.

        EXAMPLES::

            sage: 2 in Integers(7)
            True
            sage: 2 in ZZ
            True
            sage: Integers(7)(3) in ZZ
            True
            sage: 3/1 in ZZ
            True
            sage: 5 in QQ
            True
            sage: I in RR
            False
            sage: SR(2) in ZZ
            True
            sage: RIF(1, 2) in RIF
            True
            sage: pi in RIF # there is no element of RIF equal to pi
            False
            sage: sqrt(2) in CC
            True
            sage: pi in RR
            True
            sage: pi in CC
            True
            sage: pi in RDF
            True
            sage: pi in CDF
            True

        TESTS:

        Check that :trac:`13824` is fixed::

            sage: 4/3 in GF(3)
            False
            sage: 15/50 in GF(25, 'a')
            False
            sage: 7/4 in Integers(4)
            False
            sage: 15/36 in Integers(6)
            False
        """
        P = parent_c(x)
        if P is self or P == self:
            return True
        try:
            x2 = self(x)
            EQ = (x2 == x)
            if EQ is True:
                return True
            elif EQ is False:
                return False
            elif EQ:
                return True
            else:
                from sage.symbolic.expression import is_Expression
                if is_Expression(EQ):  # if comparing gives an Expression, then it must be an equation.
                    # We return *true* here, even though the equation
                    # EQ must have evaluated to False for us to get to
                    # this point. The reason is because... in practice
                    # SR is ultra-permissive about letting other rings
                    # coerce in, but ultra-strict about doing
                    # comparisons.
                    return True
                return False
        except (TypeError, ValueError, ZeroDivisionError):
            return False

    cpdef coerce(self, x):
        """
        Return x as an element of self, if and only if there is a canonical
        coercion from the parent of x to self.

        EXAMPLES::

            sage: QQ.coerce(ZZ(2))
            2
            sage: ZZ.coerce(QQ(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Integer Ring

        We make an exception for zero::

            sage: V = GF(7)^7
            sage: V.coerce(0)
            (0, 0, 0, 0, 0, 0, 0)
        """
        mor = self.coerce_map_from(parent_c(x))
        if mor is None:
            if is_Integer(x) and not x:
                try:
                    return self(0)
                except StandardError:
                    _record_exception()
            raise TypeError("no canonical coercion from %s to %s" % (parent_c(x), self))
        else:
            return (<map.Map>mor)._call_(x)

    def _list_from_iterator_cached(self):
        r"""
        Return a list of the elements of ``self``.

        OUTPUT:

        A list of all the elements produced by the iterator
        defined for the object.  The result is cached. An
        infinite set may define an iterator, allowing one
        to search through the elements, but a request by
        this method for the entire list should fail.

        NOTE:

        Some objects ``X`` do not know if they are finite or not. If
        ``X.is_finite()`` fails with a ``NotImplementedError``, then
        ``X.list()`` will simply try. In that case, it may run without
        stopping.

        However, if ``X`` knows that it is infinite, then running
        ``X.list()`` will raise an appropriate error, while running
        ``list(X)`` will run indefinitely.  For many Sage objects
        ``X``, using ``X.list()`` is preferable to using ``list(X)``.

        Nevertheless, since the whole list of elements is created and
        cached by ``X.list()``, it may be better to do ``for x in
        X:``, not ``for x in X.list():``.

        EXAMPLES::

            sage: R = Integers(11)
            sage: R.list()    # indirect doctest
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

            sage: ZZ.list()
            Traceback (most recent call last):
            ...
            ValueError: since it is infinite, cannot list Integer Ring

        This is the motivation for :trac:`10470` ::

            sage: (QQ^2).list()
            Traceback (most recent call last):
            ...
            ValueError: since it is infinite, cannot list Vector space of dimension 2 over Rational Field

        TESTS:

        The following tests the caching by adjusting the cached version::

            sage: R = Integers(3)
            sage: R.list()
            [0, 1, 2]
            sage: R._list[0] = 'junk'
            sage: R.list()
            ['junk', 1, 2]

        Here we test that for an object that does not know whether it
        is finite or not.  Calling ``X.list()`` simply tries to create
        the list (but here it fails, since the object is not
        iterable). This was fixed :trac:`11350` ::

            sage: R.<t,p> = QQ[]
            sage: Q = R.quotient(t^2-t+1)
            sage: Q.is_finite()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: Q.list()
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support iteration

        Here is another example. We artificially create a version of
        the ring of integers that does not know whether it is finite
        or not::

            sage: from sage.rings.integer_ring import IntegerRing_class
            sage: class MyIntegers_class(IntegerRing_class):
            ....:      def is_finite(self):
            ....:          raise NotImplementedError
            sage: MyIntegers = MyIntegers_class()
            sage: MyIntegers.is_finite()
            Traceback (most recent call last):
            ...
            NotImplementedError

        Asking for ``list(MyIntegers)`` below will never finish without
        pressing Ctrl-C.  We let it run for 1 second and then interrupt::

            sage: alarm(1.0); list(MyIntegers)
            Traceback (most recent call last):
            ...
            AlarmInterrupt

        """
        try:
            if self._list is not None:
                return self._list
        except AttributeError:
            pass
        # if known to be infinite, give up
        # if unsure, proceed (which will hang for an infinite set)
        try:
            if not self.is_finite():
                raise ValueError('since it is infinite, cannot list %s' % self )
        except (AttributeError, NotImplementedError):
            pass
        the_list = list(self.__iter__())
        try:
            self._list = the_list
        except AttributeError:
            pass
        return the_list

    def __nonzero__(self):
        """
        By default, all Parents are treated as True when used in an if
        statement. Override this method if other behavior is desired
        (for example, for empty sets).

        EXAMPLES::

            sage: if ZZ: print "Yes"
            Yes
        """
        return True

    cpdef bint _richcmp_helper(left, right, int op) except -2:
        """
        Compare left and right.

        EXAMPLES::

            sage: ZZ < QQ
            True
        """
        cdef int r

        if not PY_TYPE_CHECK(right, Parent) or not PY_TYPE_CHECK(left, Parent):
            # One is not a parent -- use arbitrary ordering
            if (<PyObject*>left) < (<PyObject*>right):
                r = -1
            elif (<PyObject*>left) > (<PyObject*>right):
                r = 1
            else:
                r = 0

        else:
            # Both are parents -- but need *not* have the same type.
            if HAS_DICTIONARY(left):
                r = left.__cmp__(right)
            else:
                r = left._cmp_(right)

        if op == 0:  #<
            return r  < 0
        elif op == 2: #==
            return r == 0
        elif op == 4: #>
            return r  > 0
        elif op == 1: #<=
            return r <= 0
        elif op == 3: #!=
            return r != 0
        elif op == 5: #>=
            return r >= 0

    # Should be moved and merged into the EnumeratedSets() category (#12955)
    def __len__(self):
        """
        Returns the number of elements in self. This is the naive algorithm
        of listing self and counting the elements.

        EXAMPLES::

            sage: len(GF(5))
            5
            sage: len(MatrixSpace(GF(2), 3, 3))
            512
        """
        return len(self.list())

    # Should be moved and merged into the EnumeratedSets() category (#12955)
    def __getitem__(self, n):
        """
        Returns the `n^{th}` item or slice `n` of self,
        by getting self as a list.

        EXAMPLES::

            sage: MatrixSpace(GF(3), 2, 2)[9]
            [0 2]
            [0 0]
            sage: MatrixSpace(GF(3), 2, 2)[0]
            [0 0]
            [0 0]
            sage: VectorSpace(GF(7), 3)[:10]
            [(0, 0, 0),
             (1, 0, 0),
             (2, 0, 0),
             (3, 0, 0),
             (4, 0, 0),
             (5, 0, 0),
             (6, 0, 0),
             (0, 1, 0),
             (1, 1, 0),
             (2, 1, 0)]

        TESTS:

        We test the workaround described in #12956 to let categories
        override this default implementation::

            sage: class As(Category):
            ....:     def super_categories(self): return [Sets()]
            ....:     class ParentMethods:
            ....:         def __getitem__(self, n):
            ....:             return 'coucou'
            sage: class A(Parent):
            ....:     def __init__(self):
            ....:         Parent.__init__(self, category=As())
            sage: a = A()
            sage: a[1]
            'coucou'
        """
        try:
            return super(Parent, self).__getitem__(n)
        except AttributeError:
            return self.list()[n]


    #################################################################################
    # Generators and Homomorphisms
    #################################################################################

    def _is_valid_homomorphism_(self, codomain, im_gens):
       r"""
       Return True if ``im_gens`` defines a valid homomorphism
       from self to codomain; otherwise return False.

       If determining whether or not a homomorphism is valid has not
       been implemented for this ring, then a NotImplementedError exception
       is raised.
       """
       raise NotImplementedError("Verification of correctness of homomorphisms from %s not yet implemented."%self)

    def Hom(self, codomain, category=None):
        r"""
        Return the homspace ``Hom(self, codomain, category)``.

        INPUT:

        - ``codomain`` -- a parent
        - ``category`` -- a category or ``None`` (default: ``None``)
          If ``None``, the meet of the category of ``self`` and
          ``codomain`` is used.

        OUTPUT:

        The homspace of all homomorphisms from ``self`` to
        ``codomain`` in the category ``category``.

        .. SEEALSO:: :func:`~sage.categories.homset.Hom`

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: R.Hom(QQ)
            Set of Homomorphisms from Multivariate Polynomial Ring in x, y over Rational Field to Rational Field

        Homspaces are defined for very general Sage objects, even elements of familiar rings::

            sage: n = 5; Hom(n,7)
            Set of Morphisms from 5 to 7 in Category of elements of Integer Ring
            sage: z=(2/3); Hom(z,8/1)
            Set of Morphisms from 2/3 to 8 in Category of elements of Rational Field

        This example illustrates the optional third argument::

            sage: QQ.Hom(ZZ, Sets())
            Set of Morphisms from Rational Field to Integer Ring in Category of sets

        A parent may specify how to construct certain homsets by
        implementing a method :meth:`_Hom_`(codomain, category).
        See :func:`~sage.categories.homset.Hom` for details.
        """
        from sage.categories.homset import Hom
        return Hom(self, codomain, category)

    def hom(self, im_gens, codomain=None, check=None):
       r"""
       Return the unique homomorphism from self to codomain that
       sends ``self.gens()`` to the entries of ``im_gens``.
       Raises a TypeError if there is no such homomorphism.

       INPUT:

       - ``im_gens`` -- the images in the codomain of the generators
         of this object under the homomorphism

       - ``codomain`` -- the codomain of the homomorphism

       - ``check`` -- whether to verify that the images of generators
         extend to define a map (using only canonical coercions).

       OUTPUT:

       A homomorphism self --> codomain

       .. NOTE::

          As a shortcut, one can also give an object X instead of
          ``im_gens``, in which case return the (if it exists)
          natural map to X.

       EXAMPLES:

       Polynomial Ring: We first illustrate construction of a few
       homomorphisms involving a polynomial ring::

           sage: R.<x> = PolynomialRing(ZZ)
           sage: f = R.hom([5], QQ)
           sage: f(x^2 - 19)
           6

           sage: R.<x> = PolynomialRing(QQ)
           sage: f = R.hom([5], GF(7))
           Traceback (most recent call last):
           ...
           TypeError: images do not define a valid homomorphism

           sage: R.<x> = PolynomialRing(GF(7))
           sage: f = R.hom([3], GF(49,'a'))
           sage: f
           Ring morphism:
             From: Univariate Polynomial Ring in x over Finite Field of size 7
             To:   Finite Field in a of size 7^2
             Defn: x |--> 3
           sage: f(x+6)
           2
           sage: f(x^2+1)
           3

       Natural morphism::

           sage: f = ZZ.hom(GF(5))
           sage: f(7)
           2
           sage: f
           Ring Coercion morphism:
             From: Integer Ring
             To:   Finite Field of size 5

       There might not be a natural morphism, in which case a
       ``TypeError`` is raised::

           sage: QQ.hom(ZZ)
           Traceback (most recent call last):
           ...
           TypeError: Natural coercion morphism from Rational Field to Integer Ring not defined.
       """
       if isinstance(im_gens, Parent):
           return self.Hom(im_gens).natural_map()
       from sage.structure.sequence import Sequence_generic, Sequence
       if codomain is None:
           im_gens = Sequence(im_gens)
           codomain = im_gens.universe()
       if isinstance(im_gens, (Sequence_generic, Generators)):
            im_gens = list(im_gens)
       if check is None:
           return self.Hom(codomain)(im_gens)
       else:
           return self.Hom(codomain)(im_gens, check=check)

    #################################################################################
    # New Coercion support functionality
    #################################################################################

    def _populate_coercion_lists_(self,
                                  coerce_list=[],
                                  action_list=[],
                                  convert_list=[],
                                  embedding=None,
                                  convert_method_name=None,
                                  element_constructor=None,
                                  init_no_parent=None,
                                  bint unpickling=False):
        """
        This function allows one to specify coercions, actions, conversions
        and embeddings involving this parent.

        IT SHOULD ONLY BE CALLED DURING THE __INIT__ method, often at the end.

        INPUT:

        - ``coerce_list`` -- a list of coercion Morphisms to self and
          parents with canonical coercions to self

        - ``action_list`` -- a list of actions on and by self

        - ``convert_list`` -- a list of conversion Maps to self and
           parents with conversions to self

        - ``embedding`` -- a single Morphism from self

        - ``convert_method_name`` -- a name to look for that other elements
          can implement to create elements of self (e.g. _integer_)

        - ``element_constructor`` -- A callable object used by the
          __call__ method to construct new elements. Typically the
          element class or a bound method (defaults to
          self._element_constructor_).

        - ``init_no_parent`` -- if True omit passing self in as the
          first argument of element_constructor for conversion. This
          is useful if parents are unique, or element_constructor is a
          bound method (this latter case can be detected
          automatically).
        """
        self.init_coerce(False)

        if element_constructor is None and not unpickling:
            try:
                element_constructor = self._element_constructor_
            except AttributeError:
                raise RuntimeError("element_constructor must be provided, either as an _element_constructor_ method or via the _populate_coercion_lists_ call")
        self._element_constructor = element_constructor
        self._element_init_pass_parent = guess_pass_parent(self, element_constructor)

        if not PY_TYPE_CHECK(coerce_list, list):
            raise ValueError("%s_populate_coercion_lists_: coerce_list is type %s, must be list" % (type(coerce_list), type(self)))
        if not PY_TYPE_CHECK(action_list, list):
            raise ValueError("%s_populate_coercion_lists_: action_list is type %s, must be list" % (type(action_list), type(self)))
        if not PY_TYPE_CHECK(convert_list, list):
            raise ValueError("%s_populate_coercion_lists_: convert_list is type %s, must be list" % (type(convert_list), type(self)))

        self._initial_coerce_list = copy(coerce_list)
        self._initial_action_list = copy(action_list)
        self._initial_convert_list = copy(convert_list)

        self._convert_method_name = convert_method_name
        if init_no_parent is not None:
            self._element_init_pass_parent = not init_no_parent

        for mor in coerce_list:
            self.register_coercion(mor)
        for action in action_list:
            self.register_action(action)
        for mor in convert_list:
            self.register_conversion(mor)
        if embedding is not None:
            self.register_embedding(embedding)


    def _unset_coercions_used(self):
        r"""
        Pretend that this parent has never been interrogated by the coercion
        model, so that it is possible to add coercions, conversions, and
        actions.  Does not remove any existing embedding.

        WARNING::

            For internal use only!
        """
        self._coercions_used = False
        import sage.structure.element
        sage.structure.element.get_coercion_model().reset_cache()

    def _unset_embedding(self):
        r"""
        Pretend that this parent has never been interrogated by the
        coercion model, and remove any existing embedding.

        WARNING::

            This does *not* make it safe to add an entirely new embedding!  It
            is possible that a `Parent` has cached information about the
            existing embedding; that cached information *is not* removed by
            this call.

            For internal use only!
        """
        self._embedding = None
        self._unset_coercions_used()

    cpdef bint is_coercion_cached(self, domain):
        """

        """
        return domain in self._coerce_from_hash

    cpdef bint is_conversion_cached(self, domain):
        """
        """
        return domain in self._convert_from_hash

    cpdef register_coercion(self, mor):
        r"""
        Update the coercion model to use `mor : P \to \text{self}` to coerce
        from a parent ``P`` into ``self``.

        For safety, an error is raised if another coercion has already
        been registered or discovered between ``P`` and ``self``.

        EXAMPLES::

            sage: K.<a> = ZZ['a']
            sage: L.<b> = ZZ['b']
            sage: L_into_K = L.hom([-a]) # non-trivial automorphism
            sage: K.register_coercion(L_into_K)

            sage: K(0) + b
            -a
            sage: a + b
            0
            sage: K(b) # check that convert calls coerce first; normally this is just a
            -a

            sage: L(0) + a in K # this goes through the coercion mechanism of K
            True
            sage: L(a) in L # this still goes through the convert mechanism of L
            True

            sage: K.register_coercion(L_into_K)
            Traceback (most recent call last):
            ...
            AssertionError: coercion from Univariate Polynomial Ring in b over Integer Ring to Univariate Polynomial Ring in a over Integer Ring already registered or discovered
        """
        if PY_TYPE_CHECK(mor, map.Map):
            if mor.codomain() is not self:
                raise ValueError("Map's codomain must be self (%s) is not (%s)" % (self, mor.codomain()))
        elif PY_TYPE_CHECK(mor, Parent) or PY_TYPE_CHECK(mor, type):
            mor = self._generic_convert_map(mor)
        else:
            raise TypeError("coercions must be parents or maps (got %s)" % type(mor))

        assert not (self._coercions_used and mor.domain() in self._coerce_from_hash), "coercion from %s to %s already registered or discovered"%(mor.domain(), self)
        self._coerce_from_list.append(mor)
        self._coerce_from_hash.set(mor.domain(),mor)

    cpdef register_action(self, action):
        r"""
        Update the coercion model to use ``action`` to act on self.

        ``action`` should be of type ``sage.categories.action.Action``.

        EXAMPLES::

            sage: import sage.categories.action
            sage: import operator

            sage: class SymmetricGroupAction(sage.categories.action.Action):
            ....:     "Act on a multivariate polynomial ring by permuting the generators."
            ....:     def __init__(self, G, M, is_left=True):
            ....:         sage.categories.action.Action.__init__(self, G, M, is_left, operator.mul)
            ....:
            ....:     def _call_(self, g, a):
            ....:         if not self.is_left():
            ....:             g, a = a, g
            ....:         D = {}
            ....:         for k, v in a.dict().items():
            ....:             nk = [0]*len(k)
            ....:             for i in range(len(k)):
            ....:                 nk[g(i+1)-1] = k[i]
            ....:             D[tuple(nk)] = v
            ....:         return a.parent()(D)

            sage: R.<x, y, z> = QQ['x, y, z']
            sage: G = SymmetricGroup(3)
            sage: act = SymmetricGroupAction(G, R)
            sage: t = x + 2*y + 3*z

            sage: act(G((1, 2)), t)
            2*x + y + 3*z
            sage: act(G((2, 3)), t)
            x + 3*y + 2*z
            sage: act(G((1, 2, 3)), t)
            3*x + y + 2*z

        This should fail, since we haven't registered the left
        action::

            sage: G((1,2)) * t
            Traceback (most recent call last):
            ...
            TypeError: ...

        Now let's make it work::

            sage: R._unset_coercions_used()
            sage: R.register_action(act)
            sage: G((1, 2)) * t
            2*x + y + 3*z
        """
        assert not self._coercions_used, "coercions must all be registered up before use"
        from sage.categories.action import Action
        if isinstance(action, Action):
            if action.actor() is self:
                self._action_list.append(action)
                self._action_hash.set(action.domain(), action.operation(), action.is_left(), action)
            elif action.domain() is self:
                self._action_list.append(action)
                self._action_hash.set(action.actor(), action.operation(), not action.is_left(), action)
            else:
                raise ValueError("Action must involve self")
        else:
            raise TypeError("actions must be actions")

    cpdef register_conversion(self, mor):
        r"""
        Update the coercion model to use `\text{mor} : P \to \text{self}` to convert
        from ``P`` into ``self``.

        EXAMPLES::

            sage: K.<a> = ZZ['a']
            sage: M.<c> = ZZ['c']
            sage: M_into_K = M.hom([a]) # trivial automorphism
            sage: K._unset_coercions_used()
            sage: K.register_conversion(M_into_K)

            sage: K(c)
            a
            sage: K(0) + c
            Traceback (most recent call last):
            ...
            TypeError: ...
        """
        assert not (self._coercions_used and mor.domain() in self._convert_from_hash), "conversion from %s to %s already registered or discovered"%(mor.domain(), self)
        if isinstance(mor, map.Map):
            if mor.codomain() is not self:
                raise ValueError("Map's codomain must be self")
            self._convert_from_list.append(mor)
            self._convert_from_hash.set(mor.domain(),mor)
        elif PY_TYPE_CHECK(mor, Parent) or PY_TYPE_CHECK(mor, type):
            t = mor
            mor = self._generic_convert_map(mor)
            self._convert_from_list.append(mor)
            self._convert_from_hash.set(t, mor)
            self._convert_from_hash.set(mor.domain(), mor)
        else:
            raise TypeError("conversions must be parents or maps")

    cpdef register_embedding(self, embedding):
        r"""
        Add embedding to coercion model.

        This method updates the coercion model to use
        `\text{embedding} : \text{self} \to P` to embed ``self`` into
        the parent ``P``.

        There can only be one embedding registered; it can only be registered
        once; and it must be registered before using this parent in the
        coercion model.

        EXAMPLES::

            sage: S3 = AlternatingGroup(3)
            sage: G = SL(3, QQ)
            sage: p = S3[2]; p.matrix()
            [0 0 1]
            [1 0 0]
            [0 1 0]

        By default, one can't mix matrices and permutations::

            sage: G(p)
            Traceback (most recent call last):
            ...
            TypeError: entries must be coercible to a list or integer
            sage: G(1) * p
            Traceback (most recent call last):
            ...
            TypeError: ...
            sage: phi = S3.hom(lambda p: G(p.matrix()), codomain = G)
            sage: phi(p)
            [0 0 1]
            [1 0 0]
            [0 1 0]
            sage: S3._unset_coercions_used()
            sage: S3.register_embedding(phi)
            sage: S3.coerce_embedding()
            Generic morphism:
              From: Alternating group of order 3!/2 as a permutation group
              To:   Special Linear Group of degree 3 over Rational Field
            sage: S3.coerce_embedding()(p)
            [0 0 1]
            [1 0 0]
            [0 1 0]

        This does not work since matrix groups are still old-style
        parents (see :trac:`14014`)::

            sage: G(p)                               # todo: not implemented
            sage: G(1) * p                           # todo: not implemented

        Some more advanced examples::

            sage: x = QQ['x'].0
            sage: t = abs(ZZ.random_element(10^6))
            sage: K = NumberField(x^2 + 2*3*7*11, "a"+str(t))
            sage: a = K.gen()
            sage: K_into_MS = K.hom([a.matrix()])
            sage: K._unset_coercions_used()
            sage: K.register_embedding(K_into_MS)

            sage: L = NumberField(x^2 + 2*3*7*11*19*31, "b"+str(abs(ZZ.random_element(10^6))))
            sage: b = L.gen()
            sage: L_into_MS = L.hom([b.matrix()])
            sage: L._unset_coercions_used()
            sage: L.register_embedding(L_into_MS)

            sage: K.coerce_embedding()(a)
            [   0    1]
            [-462    0]
            sage: L.coerce_embedding()(b)
            [      0       1]
            [-272118       0]

            sage: a.matrix() * b
            [-272118       0]
            [      0    -462]
            sage: a * b.matrix()
            [-272118       0]
            [      0    -462]
        """
        assert not self._coercions_used, "coercions must all be registered up before use"
        assert self._embedding is None, "only one embedding allowed"

        if isinstance(embedding, map.Map):
            if embedding.domain() is not self:
                raise ValueError("embedding's domain must be self")
            self._embedding = embedding
        elif isinstance(embedding, Parent):
            self._embedding = embedding._generic_convert_map(self)
        elif embedding is not None:
            raise TypeError("embedding must be a parent or map")

    def coerce_embedding(self):
        """
        Returns the embedding of self into some other parent, if such a parent
        exists.

        This does not mean that there are no coercion maps from self into other
        fields, this is simply a specific morphism specified out of self and
        usually denotes a special relationship (e.g. sub-objects, choice of
        completion, etc.)

        EXAMPLES::

            sage: K.<a>=NumberField(x^3+x^2+1,embedding=1)
            sage: K.coerce_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 + x^2 + 1
              To:   Real Lazy Field
              Defn: a -> -1.465571231876768?
            sage: K.<a>=NumberField(x^3+x^2+1,embedding=CC.gen())
            sage: K.coerce_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 + x^2 + 1
              To:   Complex Lazy Field
              Defn: a -> 0.2327856159383841? + 0.7925519925154479?*I
        """
        return self._embedding

    cpdef _generic_convert_map(self, S):
        r"""
        Returns the default conversion map based on the data provided to
        :meth:`_populate_coercion_lists_`.

        This called when :meth:`_coerce_map_from_` returns ``True``.

        If a ``convert_method_name`` is provided, it creates a
        ``NamedConvertMap``, otherwise it creates a
        ``DefaultConvertMap`` or ``DefaultConvertMap_unique``
        depending on whether or not init_no_parent is set.

        EXAMPLES:
        """
        import coerce_maps
        if self._convert_method_name is not None:
            # handle methods like _integer_
            if PY_TYPE_CHECK(S, type):
                element_constructor = S
            elif PY_TYPE_CHECK(S, Parent):
                element_constructor = (<Parent>S)._element_constructor
                if not PY_TYPE_CHECK(element_constructor, type):
                    # if element_constructor is not an actual class, get the element class
                    element_constructor = type(S.an_element())
            else:
                element_constructor = None
            if element_constructor is not None and hasattr(element_constructor, self._convert_method_name):
                return coerce_maps.NamedConvertMap(S, self, self._convert_method_name)

        if self._element_init_pass_parent:
            return coerce_maps.DefaultConvertMap(S, self)
        else:
            return coerce_maps.DefaultConvertMap_unique(S, self)

    def _coerce_map_via(self, v, S):
        """
        This attempts to construct a morphism from S to self by passing through
        one of the items in v (tried in order).

        S may appear in the list, in which case algorithm will never progress
        beyond that point.

        This is similar in spirit to the old {{{_coerce_try}}}, and useful when
        defining _coerce_map_from_

        INPUT:

        - ``v`` - A list (iterator) of parents with coercions into self. There
          MUST be maps provided from each item in the list to self.

        - ``S`` - the starting parent

        EXAMPLES::

            sage: CDF._coerce_map_via([ZZ, RR, CC], int)
            Composite map:
              From: Set of Python objects of type 'int'
              To:   Complex Double Field
              Defn:   Native morphism:
                      From: Set of Python objects of type 'int'
                      To:   Integer Ring
                    then
                      Native morphism:
                      From: Integer Ring
                      To:   Complex Double Field

            sage: CDF._coerce_map_via([ZZ, RR, CC], QQ)
            Composite map:
              From: Rational Field
              To:   Complex Double Field
              Defn:   Generic map:
                      From: Rational Field
                      To:   Real Field with 53 bits of precision
                    then
                      Native morphism:
                      From: Real Field with 53 bits of precision
                      To:   Complex Double Field

            sage: CDF._coerce_map_via([ZZ, RR, CC], CC)
            Generic map:
              From: Complex Field with 53 bits of precision
              To:   Complex Double Field
        """
        cdef Parent R
        for R in v:
            if R is None:
                continue
            if R is S:
                return self.coerce_map_from(R)
            connecting = R.coerce_map_from(S)
            if connecting is not None:
                return self.coerce_map_from(R) * connecting

    cpdef bint has_coerce_map_from(self, S) except -2:
        """
        Return True if there is a natural map from S to self.
        Otherwise, return False.

        EXAMPLES::

            sage: RDF.has_coerce_map_from(QQ)
            True
            sage: RDF.has_coerce_map_from(QQ['x'])
            False
            sage: RDF['x'].has_coerce_map_from(QQ['x'])
            True
            sage: RDF['x,y'].has_coerce_map_from(QQ['x'])
            True
        """
        if S is self:
            return True
        elif S == self:
            if debug.unique_parent_warnings:
                print "Warning: non-unique parents %s"%(type(S))
            return True
        return self.coerce_map_from(S) is not None

    cpdef _coerce_map_from_(self, S):
        """
        Override this method to specify coercions beyond those specified
        in coerce_list.

        If no such coercion exists, return None or False. Otherwise, it may
        return either an actual Map to use for the coercion, a callable
        (in which case it will be wrapped in a Map), or True (in which case
        a generic map will be provided).
        """
        return None

    cpdef coerce_map_from(self, S):
        """
        This returns a Map object to coerce from S to self if one exists,
        or None if no such coercion exists.

        EXAMPLES::

            sage: ZZ.coerce_map_from(int)
            Native morphism:
              From: Set of Python objects of type 'int'
              To:   Integer Ring
            sage: QQ.coerce_map_from(ZZ)
            Natural morphism:
              From: Integer Ring
              To:   Rational Field

        By :trac:`12313`, a special kind of weak key dictionary is used to
        store coercion and conversion maps, namely
        :class:`~sage.structure.coerce_dict.MonoDict`. In that way, a memory
        leak was fixed that would occur in the following test::

            sage: import gc
            sage: _ = gc.collect()
            sage: K = GF(1<<55,'t')
            sage: for i in range(50):
            ....:   a = K.random_element()
            ....:   E = EllipticCurve(j=a)
            ....:   b = K.has_coerce_map_from(E)
            sage: _ = gc.collect()
            sage: len([x for x in gc.get_objects() if isinstance(x,type(E))])
            1

        TESTS:

        The following was fixed in :trac:`12969`::

            sage: R = QQ['q,t'].fraction_field()
            sage: Sym = sage.combinat.sf.sf.SymmetricFunctions(R)
            sage: H = Sym.macdonald().H()
            sage: P = Sym.macdonald().P()
            sage: m = Sym.monomial()
            sage: Ht = Sym.macdonald().Ht()
            sage: phi = m.coerce_map_from(P)
            sage: Ht.coerce_map_from(P)
            Composite map:
              From: Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald P basis
              To:   Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald Ht basis
              Defn:   Composite map:
                      From: Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald P basis
                      To:   Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Schur basis
                      Defn:   Generic morphism:
                              From: Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald P basis
                              To:   Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald J basis
                            then
                              Generic morphism:
                              From: Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald J basis
                              To:   Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Schur basis
                    then
                      Generic morphism:
                      From: Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Schur basis
                      To:   Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald Ht basis

        The following was fixed in :trac:`4740`::

            sage: F = GF(13)
            sage: F.coerce_map_from(F) is F.coerce_map_from(F)
            True

        """
        if not good_as_coerce_domain(S):
            return None
        self._coercions_used = True
        cdef map.Map mor

        if isinstance(S, Set_PythonType_class):
            return self.coerce_map_from(S._type)
        if self._coerce_from_hash is None: # this is because parent.__init__() does not always get called
            self.init_coerce(False)

        try:
            return self._coerce_from_hash.get(S)
        except KeyError:
            pass

        if S is self:
            from sage.categories.homset import Hom
            mor = Hom(self, self).identity()
            self._coerce_from_hash.set(S, mor)
            return mor

        if S == self:
            # non-unique parents
            if debug.unique_parent_warnings:
                print "Warning: non-unique parents %s"%(type(S))
            mor = self._generic_convert_map(S)
            self._coerce_from_hash.set(S, mor)
            return mor

        try:
            _register_pair(self, S, "coerce")
            mor = self.discover_coerce_map_from(S)
            #if mor is not None:
            #    # Need to check that this morphism doesn't connect previously unconnected parts of the coercion diagram
            #    if self._embedding is not None and not self._embedding.codomain().has_coerce_map_from(S):
            #        # The following if statement may call this function with self and S.  If so, we want to return None,
            #        # so that it doesn't use this path for the existence of a coercion path.
            #        # We disable this for now because it is too strict
            #        pass
            #        # print "embed problem: the following morphisms connect unconnected portions of the coercion graph\n%s\n%s"%(self._embedding, mor)
            #        # mor = None
            if mor is not None:
                # NOTE: this line is what makes the coercion detection stateful
                # self._coerce_from_list.append(mor)
                pass
            # It may be that the only coercion from S to self is
            # via another parent X. But if the pair (S,X) is temporarily
            # disregarded (using _register_pair, to avoid infinite recursion)
            # then we are not allowed to cache the absence of a coercion
            # from S to self. See #12969
            if (mor is not None) or _may_cache_none(self, S, "coerce"):
                self._coerce_from_hash.set(S,mor)
            return mor
        except CoercionException, ex:
            _record_exception()
            return None
        finally:
            _unregister_pair(self, S, "coerce")


    cdef discover_coerce_map_from(self, S):
        """
        Precedence for discovering a coercion S -> self goes as follows:

        1. If S has an embedding into T, look for T -> self and return composition

        2. If self._coerce_map_from_(S) is NOT exactly one of

           - DefaultConvertMap
           - DefaultConvertMap_unique
           - NamedConvertMap

           return this map

        3. Traverse the coercion lists looking for another map
           returning the map from step (2) if none is found.

        In the future, multiple paths may be discovered and compared.

        TESTS:

        Regression test for :trac:`12919` (probably not 100% robust)::

            sage: class P(Parent):
            ....:     def __init__(self):
            ....:         Parent.__init__(self, category=Sets())
            ....:     Element=ElementWrapper
            sage: A = P(); a = A('a')
            sage: B = P(); b = B('b')
            sage: C = P(); c = C('c')
            sage: D = P(); d = D('d')
            sage: Hom(A, B)(lambda x: b).register_as_coercion()
            sage: Hom(B, A)(lambda x: a).register_as_coercion()
            sage: Hom(C, B)(lambda x: b).register_as_coercion()
            sage: Hom(D, C)(lambda x: c).register_as_coercion()
            sage: A(d)
            'a'

        Another test::

            sage: K = NumberField([x^2-2, x^2-3], 'a,b')
            sage: M = K.absolute_field('c')
            sage: M_to_K, K_to_M = M.structure()
            sage: M.register_coercion(K_to_M)
            sage: K.register_coercion(M_to_K)
            sage: phi = M.coerce_map_from(QQ)
            sage: p = QQ.random_element()
            sage: c = phi(p) - p; c
            0
            sage: c.parent() is M
            True
            sage: K.coerce_map_from(QQ)
            Conversion map:
            From: Rational Field
            To:   Number Field in a with defining polynomial x^2 - 2 over its base field
        """
        best_mor = None
        if PY_TYPE_CHECK(S, Parent) and (<Parent>S)._embedding is not None:
            if (<Parent>S)._embedding.codomain() is self:
                return (<Parent>S)._embedding
            connecting = self.coerce_map_from((<Parent>S)._embedding.codomain())
            if connecting is not None:
                return (<Parent>S)._embedding.post_compose(connecting)

        cdef map.Map mor
        user_provided_mor = self._coerce_map_from_(S)

        if user_provided_mor is False:
            user_provided_mor = None

        elif user_provided_mor is not None:

            from sage.categories.map import Map
            from coerce_maps import DefaultConvertMap, DefaultConvertMap_unique, NamedConvertMap, CallableConvertMap

            if user_provided_mor is True:
                mor = self._generic_convert_map(S)
            elif PY_TYPE_CHECK(user_provided_mor, Map):
                mor = <map.Map>user_provided_mor
            elif callable(user_provided_mor):
                mor = CallableConvertMap(user_provided_mor)
            else:
                raise TypeError("_coerce_map_from_ must return None, a boolean, a callable, or an explicit Map (called on %s, got %s)" % (type(self), type(user_provided_mor)))

            if (PY_TYPE_CHECK_EXACT(mor, DefaultConvertMap) or
                  PY_TYPE_CHECK_EXACT(mor, DefaultConvertMap_unique) or
                  PY_TYPE_CHECK_EXACT(mor, NamedConvertMap)) and not mor._force_use:
                # If there is something better in the list, try to return that instead
                # This is so, for example, _coerce_map_from_ can return True but still
                # take advantage of the _populate_coercion_lists_ data.
                best_mor = mor
            else:
                return mor

        from sage.categories.homset import Hom

        cdef int num_paths = 1 # this is the number of paths we find before settling on the best (the one with lowest coerce_cost).
                               # setting this to 1 will make it return the first path found.
        cdef int mor_found = 0
        cdef Parent R
        # Recurse.  Note that if S is the domain of one of the maps in self._coerce_from_list,
        # we will have stuck the map into _coerce_map_hash and thus returned it already.
        for mor in self._coerce_from_list:
            if mor._domain is self:
                continue
            if mor._domain is S:
                if best_mor is None or mor._coerce_cost < best_mor._coerce_cost:
                    best_mor = mor
                mor_found += 1
                if mor_found  >= num_paths:
                    return best_mor
            else:
                connecting = None
                if EltPair(mor._domain, S, "coerce") not in _coerce_test_dict:
                    connecting = mor._domain.coerce_map_from(S)
                if connecting is not None:
                    mor = mor * connecting
                    if best_mor is None or mor._coerce_cost < best_mor._coerce_cost:
                        best_mor = mor
                    mor_found += 1
                    if mor_found  >= num_paths:
                        return best_mor

        return best_mor


    cpdef convert_map_from(self, S):
        """
        This function returns a Map from S to self, which may or may not
        succeed on all inputs. If a coercion map from S to self exists,
        then the it will be returned. If a coercion from self to S exists,
        then it will attempt to return a section of that map.

        Under the new coercion model, this is the fastest way to convert
        elements of S to elements of self (short of manually constructing
        the elements) and is used by __call__.

        EXAMPLES::

            sage: m = ZZ.convert_map_from(QQ)
            sage: m(-35/7)
            -5
            sage: parent(m(-35/7))
            Integer Ring
        """
        if not good_as_convert_domain(S):
            return None
        if self._convert_from_hash is None: # this is because parent.__init__() does not always get called
            self.init_coerce()
        try:
            return self._convert_from_hash.get(S)
        except KeyError:
            mor = self.discover_convert_map_from(S)
            self._convert_from_list.append(mor)
            self._convert_from_hash.set(S, mor)
            return mor

    cdef discover_convert_map_from(self, S):

        cdef map.Map mor = self.coerce_map_from(S)
        if mor is not None:
            return mor

        if PY_TYPE_CHECK(S, Parent):
            mor = S.coerce_map_from(self)
            if mor is not None:
                mor = mor.section()
                if mor is not None:
                    return mor

        user_provided_mor = self._convert_map_from_(S)

        if user_provided_mor is not None:
            if PY_TYPE_CHECK(user_provided_mor, map.Map):
                return user_provided_mor
            elif callable(user_provided_mor):
                from coerce_maps import CallableConvertMap
                return CallableConvertMap(user_provided_mor)
            else:
                raise TypeError("_convert_map_from_ must return a map or callable (called on %s, got %s)" % (type(self), type(user_provided_mor)))

        if not PY_TYPE_CHECK(S, type) and not PY_TYPE_CHECK(S, Parent):
            # Sequences is used as a category and a "Parent"
            from sage.categories.category_types import Sequences
            from sage.structure.coerce_maps import ListMorphism
            if isinstance(S, Sequences):
                return ListMorphism(S, self.convert_map_from(list))

        mor = self._generic_convert_map(S)
        return mor

    cpdef _convert_map_from_(self, S):
        """
        Override this method to provide additional conversions beyond those
        given in convert_list.

        This function is called after coercions are attempted. If there is a
        coercion morphism in the opposite direction, one should consider
        adding a section method to that.

        This MUST return a Map from S to self, or None. If None is returned
        then a generic map will be provided.
        """
        return None

    cpdef get_action(self, S, op=operator.mul, bint self_on_left=True, self_el=None, S_el=None):
        """
        Returns an action of self on S or S on self.

        To provide additional actions, override :meth:`_get_action_`.

        TESTS::

            sage: M = QQ['y']^3
            sage: M.get_action(ZZ['x']['y'])
            Right scalar multiplication by Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integer Ring on Ambient free module of rank 3 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field
            sage: M.get_action(ZZ['x']) # should be None
        """
        try:
            if self._action_hash is None: # this is because parent.__init__() does not always get called
                self.init_coerce()
            return self._action_hash.get(S, op, self_on_left)
        except KeyError:
            pass

        action = self._get_action_(S, op, self_on_left)
        if action is None:
            action = self.discover_action(S, op, self_on_left, self_el, S_el)

        if action is not None:
            from sage.categories.action import Action
            if not isinstance(action, Action):
                raise TypeError("get_action_impl must return None or an Action")
            # We do NOT add to the list, as this would lead to errors as in
            # the example above.

        self._action_hash.set(S, op, self_on_left, action)
        return action


    cdef discover_action(self, S, op, bint self_on_left, self_el=None, S_el=None):
        # G acts on S, G -> G', R -> S => G' acts on R (?)
        # NO! ZZ[x,y] acts on Matrices(ZZ[x]) but ZZ[y] does not.
        # What may be true is that if the action's destination is S, then this can be allowed.
        # Note: a is either None or a sample elements of self.
        # If needed, it will be passed to Left/RightModuleAction.
        from sage.categories.action import Action, PrecomposedAction
        from sage.categories.homset import Hom
        from coerce_actions import LeftModuleAction, RightModuleAction
        cdef Parent R
        for action in self._action_list:
            if PY_TYPE_CHECK(action, Action) and action.operation() is op:
                if self_on_left:
                    if action.left_domain() is not self: continue
                    R = action.right_domain()
                else:
                    if action.right_domain() is not self: continue
                    R = action.left_domain()
            elif op is operator.mul and isinstance(action, Parent):
                try:
                    R = action
                    _register_pair(self, R, "action") # to kill circular recursion
                    if self_on_left:
                        action = LeftModuleAction(R, self, a=S_el, g=self_el) # self is acted on from right
                    else:
                        action = RightModuleAction(R, self, a=S_el, g=self_el) # self is acted on from left
                    ## The following two lines are disabled to prevent the following from working:
                    ## sage: x, y = var('x,y')
                    ## sage: parent(ZZ[x][y](1)*vector(QQ[y],[1,2]))
                    ## sage: parent(ZZ[x](1)*vector(QQ[y],[1,2]))
                    ## We will hopefully come up with a way to reinsert them, because they increase the scope
                    ## of discovered actions.
                    #i = self._action_list.index(R)
                    #self._action_list[i] = action
                except CoercionException:
                    _record_exception()
                    continue
                finally:
                    _unregister_pair(self, R, "action")
            else:
                continue # only try mul if not specified
            if R is S:
                return action
            else:
                connecting = R.coerce_map_from(S) # S -> R
                if connecting is not None:
                    if self_on_left:
                        return PrecomposedAction(action, None, connecting)
                    else:
                        return PrecomposedAction(action, connecting, None)

        # We didn't find an action in the list, but maybe the elements
        # define special action methods
        if op is operator.mul:
            # TODO: if _xmul_/_x_action_ code does stuff like
            # if self == 0:
            #    return self
            # then an_element() == 0 could be very bad.
            try:
                _register_pair(self, S, "action") # this is to avoid possible infinite loops

                # detect actions defined by _rmul_, _lmul_, _act_on_, and _acted_upon_ methods
                from coerce_actions import detect_element_action
                action = detect_element_action(self, S, self_on_left, self_el, S_el)
                if action is not None:
                    return action

                try:
                    # maybe there is a more clever way of detecting ZZ than importing here...
                    from sage.rings.integer_ring import ZZ
                    if S is ZZ and not self.has_coerce_map_from(ZZ):
                        from sage.structure.coerce_actions import IntegerMulAction
                        action = IntegerMulAction(S, self, not self_on_left, self_el)
                        return action
                except (CoercionException, TypeError):
                    _record_exception()

            finally:
                _unregister_pair(self, S, "action")


    cpdef _get_action_(self, S, op, bint self_on_left):
        """
        Override this method to provide an action of self on S or S on self
        beyond what was specified in action_list.

        This must return an action which accepts an element of self and an
        element of S (in the order specified by self_on_left).
        """
        return None

    def construction(self):
        """
        Returns a pair (functor, parent) such that functor(parent) return self.
        If this ring does not have a functorial construction, return None.

        EXAMPLES::

            sage: QQ.construction()
            (FractionField, Integer Ring)
            sage: f, R = QQ['x'].construction()
            sage: f
            Poly[x]
            sage: R
            Rational Field
            sage: f(R)
            Univariate Polynomial Ring in x over Rational Field
        """
        return None

    # TODO: remove once all parents in Sage will inherit properly from
    # Sets().ParentMethods.an_element
    cpdef an_element(self):
        r"""
        Returns a (preferably typical) element of this parent.

        This is used both for illustration and testing purposes. If
        the set ``self`` is empty, :meth:`an_element` raises the
        exception :class:`EmptySetError`.

        This calls :meth:`_an_element_` (which see), and caches the
        result. Parent are thus encouraged to override :meth:`_an_element_`.

        EXAMPLES::

            sage: CDF.an_element()
            1.0*I
            sage: ZZ[['t']].an_element()
            t

        In case the set is empty, an :class:`EmptySetError` is raised::

            sage: Set([]).an_element()
            Traceback (most recent call last):
            ...
            EmptySetError
        """
        if self.__an_element is None:
            self.__an_element = self._an_element_()
        return self.__an_element

    def _an_element_(self):
        """
        Returns an element of self. Want it in sufficient generality
        that poorly-written functions won't work when they're not
        supposed to. This is cached so doesn't have to be super fast.

        EXAMPLES::

            sage: QQ._an_element_()
            1/2
            sage: ZZ['x,y,z']._an_element_()
            x

        TESTS:

        Since ``Parent`` comes before the parent classes provided by
        categories in the hierarchy of classes, we make sure that this
        default implementation of :meth:`_an_element_` does not
        override some provided by the categories.  Eventually, this
        default implementation should be moved into the categories to
        avoid this workaround::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.category()
            Category of facade finite enumerated sets
            sage: super(Parent, S)._an_element_
            Cached version of <function _an_element_from_iterator at ...>
            sage: S._an_element_()
            1
            sage: S = FiniteEnumeratedSet([])
            sage: S._an_element_()
            Traceback (most recent call last):
            ...
            EmptySetError

        """
        try:
            return super(Parent, self)._an_element_()
        except EmptySetError:
            raise
        except StandardError:
            _record_exception()
            pass

        try:
            return self.gen(0)
        except StandardError:
            _record_exception()
            pass

        try:
            return self.gen()
        except StandardError:
            _record_exception()
            pass

        from sage.rings.infinity import infinity
        for x in ['_an_element_', 'pi', 1.2, 2, 1, 0, infinity]:
            # This weird looking list is to try to get an element
            # which doesn't coerce other places.
            try:
                return self(x)
            except (TypeError, NameError, NotImplementedError, AttributeError, ValueError):
                _record_exception()

        raise NotImplementedError("please implement _an_element_ for %s" % self)

    cpdef bint is_exact(self) except -2:
        """
        Test whether the ring is exact.

        .. NOTE::

            This defaults to true, so even if it does return ``True``
            you have no guarantee (unless the ring has properly
            overloaded this).

        OUTPUT:

        Return True if elements of this ring are represented exactly, i.e.,
        there is no precision loss when doing arithmetic.

        EXAMPLES::

            sage: QQ.is_exact()
            True
            sage: ZZ.is_exact()
            True
            sage: Qp(7).is_exact()
            False
            sage: Zp(7, type='capped-abs').is_exact()
            False
        """
        return True

#    cpdef base_extend(self, other, category=None):
#        """
#        EXAMPLES:
#            sage: QQ.base_extend(GF(7))
#            Traceback (most recent call last):
#            ...
#            TypeError: base extension not defined for Rational Field
#            sage: ZZ.base_extend(GF(7))
#            Finite Field of size 7
#        """
#        # Use the coerce map if a base extend is not defined in the category.
#        # this is the default implementation.
#        try:
#            if category is None:
#                method = self._categories[0].get_object_method("base_extend") # , self._categories[1:])
#            else:
#                method = category.get_object_method("base_extend")
#            if method is not None:
#                return method(self)
#            elif other.has_coerce_map_from(self):
#                return other
#            else:
#                raise TypeError, "base extension not defined for %s" % self
#        except AttributeError:
#            raise TypeError, "base extension not defined for %s" % self

############################################################################
# Set base class --
############################################################################


cdef class Set_generic(Parent): # Cannot use Parent because Element._parent is Parent
    """
    Abstract base class for sets.

    TESTS::

        sage: Set(QQ).category()
        Category of sets

    """
#     def category(self):
#         # TODO: remove once all subclasses specify their category, or
#         # the constructor of Parent sets it to Sets() by default
#         """
#         The category that this set belongs to, which is the category
#         of all sets.

#         EXAMPLES:
#             sage: Set(QQ).category()
#             Category of sets
#         """
#         from sage.categories.sets_cat import Sets
#         return Sets()

    def object(self):
        return self

    def __nonzero__(self):
        """
        A set is considered True unless it is empty, in which case it is
        considered to be False.

        EXAMPLES::

            sage: bool(Set(QQ))
            True
            sage: bool(Set(GF(3)))
            True
        """
        return not (self.is_finite() and len(self) == 0)


import types
cdef _type_set_cache = {}

cpdef Parent Set_PythonType(theType):
    """
    Return the (unique) Parent that represents the set of Python objects
    of a specified type.

    EXAMPLES::

        sage: from sage.structure.parent import Set_PythonType
        sage: Set_PythonType(list)
        Set of Python objects of type 'list'
        sage: Set_PythonType(list) is Set_PythonType(list)
        True
        sage: S = Set_PythonType(tuple)
        sage: S([1,2,3])
        (1, 2, 3)

      S is a parent which models the set of all lists:
        sage: S.category()
        Category of sets

    EXAMPLES::

        sage: R = sage.structure.parent.Set_PythonType(int)
        sage: S = sage.structure.parent.Set_PythonType(float)
        sage: Hom(R, S)
        Set of Morphisms from Set of Python objects of type 'int' to Set of Python objects of type 'float' in Category of sets

    """
    try:
        return _type_set_cache[theType]
    except KeyError:
        _type_set_cache[theType] = theSet = Set_PythonType_class(theType)
        return theSet

cdef class Set_PythonType_class(Set_generic):

    cdef _type

    def __init__(self, theType):
        """
        EXAMPLES::

            sage: S = sage.structure.parent.Set_PythonType(float)
            sage: S.category()
            Category of sets
        """
        Set_generic.__init__(self, element_constructor=theType, category=Sets())
        self._type = theType

    def __call__(self, x):
        """
        This doesn't return Elements, but actual objects of the given type.

        EXAMPLES::

            sage: S = sage.structure.parent.Set_PythonType(float)
            sage: S(5)
            5.0
            sage: S(9/3)
            3.0
            sage: S(1/3)
            0.333333333333333...

        """
        if PY_TYPE_CHECK(x,self._type):
            return x
        return self._type(x)

    def __hash__(self):
        """
        TESTS::

            sage: S = sage.structure.parent.Set_PythonType(int)
            sage: hash(S) == -hash(int)
            True
        """
        return -hash(self._type)

    cpdef int _cmp_(self, other) except -100:
        """
        Two Python type sets are considered the same if they contain the same
        type.

        EXAMPLES::

            sage: S = sage.structure.parent.Set_PythonType(int)
            sage: S == S
            True
            sage: S == sage.structure.parent.Set_PythonType(float)
            False
        """
        if self is other:
            return 0
        if isinstance(other, Set_PythonType_class):
            return cmp(self._type, other._type)
        else:
            return cmp(self._type, other)

    def __contains__(self, x):
        """
        Only things of the right type (or subtypes thereof) are considered to
        belong to the set.

        EXAMPLES::

            sage: S = sage.structure.parent.Set_PythonType(tuple)
            sage: (1,2,3) in S
            True
            sage: () in S
            True
            sage: [1,2] in S
            False
        """
        return isinstance(x, self._type)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sage.structure.parent.Set_PythonType(tuple)
            Set of Python objects of type 'tuple'
            sage: sage.structure.parent.Set_PythonType(Integer)
            Set of Python objects of type 'sage.rings.integer.Integer'
            sage: sage.structure.parent.Set_PythonType(Parent)
            Set of Python objects of type 'sage.structure.parent.Parent'
        """
        return "Set of Python objects of %s"%(str(self._type)[1:-1])

    def object(self):
        """
        EXAMPLES::

            sage: S = sage.structure.parent.Set_PythonType(tuple)
            sage: S.object()
            <type 'tuple'>
        """
        return self._type

    def cardinality(self):
        """
        EXAMPLES::

            sage: S = sage.structure.parent.Set_PythonType(bool)
            sage: S.cardinality()
            2
            sage: S = sage.structure.parent.Set_PythonType(int)
            sage: S.cardinality()
            4294967296                        # 32-bit
            18446744073709551616              # 64-bit
            sage: S = sage.structure.parent.Set_PythonType(float)
            sage: S.cardinality()
            18437736874454810627
            sage: S = sage.structure.parent.Set_PythonType(long)
            sage: S.cardinality()
            +Infinity
        """
        from sage.rings.integer import Integer
        two = Integer(2)
        if self._type is bool:
            return two
        elif self._type is int:
            import sys
            return two*sys.maxint + 2
        elif self._type is float:
            return 2 * two**52 * (two**11 - 1) + 3 # all NaN's are the same from Python's point of view
        else:
            # probably
            import sage.rings.infinity
            return sage.rings.infinity.infinity

# These functions are to guarantee that user defined _lmul_, _rmul_,
# _act_on_, _acted_upon_ do not in turn call __mul__ on their
# arguments, leading to an infinite loop.

cdef dict _coerce_test_dict = {}

cdef class EltPair:
    cdef x, y, tag
    def __init__(self, x, y, tag):
        self.x = x
        self.y = y
        self.tag = tag

    def __richcmp__(EltPair self, EltPair other, int op):
        cdef bint eq = self.x is other.x and self.y is other.y and self.tag is other.tag
        if op in [Py_EQ, Py_GE, Py_LE]:
            return eq
        else:
            return not eq

    def __hash__(self):
        """
        EXAMPLES::

            sage: from sage.structure.parent import EltPair
            sage: a = EltPair(ZZ, QQ, "coerce")
            sage: hash(a) == hash((ZZ, QQ, "coerce"))
            True
        """
        return hash((self.x, self.y, self.tag))

    def short_repr(self):
        return self.tag, hex(<long><void*>self.x), hex(<long><void*>self.y)

    def __repr__(self):
        return "%r: %r (%r), %r (%r)" % (self.tag, self.x, type(self.x), self.y, type(self.y))

cdef bint _may_cache_none(x, y, tag) except -1:
    # Are we allowed to cache the absence of a coercion
    # from y to x? We are only allowed, if y is *not*
    # part of any coerce path that is temporarily disregarded,
    # with the only exception of the path from y to x.
    # See #12969.
    cdef EltPair P
    for P in _coerce_test_dict.iterkeys():
        if (P.y is y) and (P.x is not x) and (P.tag is tag):
            return 0
    return 1

cdef bint _register_pair(x, y, tag) except -1:
    # Means: We will temporarily disregard coercions from
    # y to x when looking for a coercion path by depth first
    # search. This is to avoid infinite recursion.
    both = EltPair(x,y,tag)

    if both in _coerce_test_dict:
        xp = type(x) if PY_TYPE_CHECK(x, Parent) else parent_c(x)
        yp = type(y) if PY_TYPE_CHECK(y, Parent) else parent_c(y)
        raise CoercionException("Infinite loop in action of %s (parent %s) and %s (parent %s)!" % (x, xp, y, yp))
    _coerce_test_dict[both] = True
    return 0

cdef bint _unregister_pair(x, y, tag) except -1:
    try:
        _coerce_test_dict.pop(EltPair(x,y,tag), None)
    except (ValueError, CoercionException):
        pass

empty_set = Set_generic()

def normalize_names(ngens, names):
    """
    TESTS::

        sage: sage.structure.parent.normalize_names(5, 'x')
        ('x0', 'x1', 'x2', 'x3', 'x4')
        sage: sage.structure.parent.normalize_names(2, ['x','y'])
        ('x', 'y')
    """
    return empty_set.normalize_names(ngens, names)
