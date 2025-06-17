.. -*- coding: utf-8 -*-

.. _coercion_and_categories:

=================================================
How to implement new algebraic structures in Sage
=================================================

.. contents::
   :depth: 3
   :class: this-will-duplicate-information-and-it-is-still-useful-here

--------------------------------------
Sage's category and coercion framework
--------------------------------------

.. MODULEAUTHOR::
    Simon King,
    Friedrich\-Schiller\-Universität Jena,
    <simon.king@uni-jena.de>
    © 2011/2013

.. linkall

The aim of this tutorial is to explain how one can benefit from Sage's
category framework and coercion model when implementing new algebraic
structures. It is based on a worksheet created in 2011.

We illustrate the concepts of Sage's category framework and coercion model by
means of a detailed example, namely a toy implementation of fraction fields.
The code is developed step by step, so that the reader can focus on one detail
in each part of this tutorial. The complete code can be found in the appendix.

Outline
=======

- Use existing base classes

  For using Sage's coercion system, it is essential to work with sub\--classes
  of :class:`sage.structure.parent.Parent` or
  :class:`sage.structure.element.Element`, respectively. They provide default
  implementations of many "magical" double-underscore Python methods, which
  must not be overridden. Instead, the actual implementation should be in
  *single underscore* methods, such as ``_add_`` or ``_mul_``.

- Turn your parent structure into an object of a category

  Declare the category during initialisation\---Your parent structure will
  inherit further useful methods and consistency tests.

- Provide your parent structure with an element class

  Assign to it an attribute called ``Element``\---The elements will inherit
  further useful methods from the category. In addition, some basic
  conversions will immediately work.

- Implement further conversions

  Never override a parent's ``__call__`` method! Provide the method
  ``_element_constructor_`` instead.

- Declare coercions

  If a conversion happens to be a morphism, you may consider to turn it into a
  coercion. It will then *implicitly* be used in arithmetic operations.

- Advanced coercion:  Define construction functors for your parent structure

  Sage will automatically create new parents for you when needed, by the
  so\--called :func:`sage.categories.pushout.pushout` construction.

- Run the automatic test suites

  Each method should be documented and provide a doc test (we are not giving
  examples here). In addition, any method defined for the objects or elements
  of a category should be supported by a test method, that is executed when
  running the test suite.

Base classes
============

In Sage, a "Parent" is an object of a category and contains elements.  Parents
should inherit from :class:`sage.structure.parent.Parent` and their elements
from :class:`sage.structure.element.Element`.

Sage provides appropriate sub\--classes of
:class:`~sage.structure.parent.Parent` and
:class:`~sage.structure.element.Element` for a variety of more concrete
algebraic structures, such as groups, rings, or fields, and of their
elements. But some old stuff in Sage doesn't use it.  **Volunteers for
refactoring are welcome!**


The parent
----------

Since we wish to implement a special kind of fields, namely fraction fields,
it makes sense to build on top of the base class
:class:`sage.rings.ring.Field` provided by Sage.  ::

    sage: from sage.rings.ring import Field


This base class provides a lot more methods than a general parent::

    sage: [p for p in dir(Field) if p not in dir(Parent)]
    ['_CommutativeRing__fraction_field',
     '__iter__',
     '__len__',
     '__rxor__',
     '__xor__',
     '_coerce_',
     '_coerce_c',
     '_coerce_impl',
     '_default_category',
     '_gens',
     '_latex_names',
     '_list',
     '_one_element',
     '_pseudo_fraction_field',
     '_zero_element',
     'an_embedding',
     'base_extend',
     'epsilon',
     'extension',
     'fraction_field',
     'gen',
     'gens',
     'ngens',
     'one',
     'order',
     'zero']

The following is a very basic implementation of fraction fields, that needs to
be complemented later.
::

    sage: from sage.structure.unique_representation import UniqueRepresentation
    sage: class MyFrac(UniqueRepresentation, Field):
    ....:     def __init__(self, base):
    ....:         if base not in IntegralDomains():
    ....:             raise ValueError("%s is no integral domain" % base)
    ....:         Field.__init__(self, base)
    ....:     def _repr_(self):
    ....:         return "NewFrac(%s)"%repr(self.base())
    ....:     def base_ring(self):
    ....:         return self.base().base_ring()
    ....:     def characteristic(self):
    ....:         return self.base().characteristic()

.. end ouf output

This basic implementation is formed by the following steps:

- Any ring in Sage has a **base** and a **base ring**. The "usual" fraction
  field of a ring `R` has the base `R` and the base ring ``R.base_ring()``::

      sage: Frac(QQ['x']).base(), Frac(QQ['x']).base_ring()
      (Univariate Polynomial Ring in x over Rational Field, Rational Field)


  Declaring the base is easy: We just pass it as an argument to the field
  constructor.
  ::

      sage: Field(ZZ['x']).base()
      Univariate Polynomial Ring in x over Integer Ring

  .. end of output

  We are implementing a separate method returning the base ring.

- Python uses double\--underscore methods for arithmetic methods and string
  representations. Sage's base classes often have a default implementation,
  and it is requested to **implement SINGLE underscore methods _repr_, and
  similarly _add_, _mul_ etc.**

- You are encouraged to **make your parent "unique"**. That's to say, parents
  should only evaluate equal if they are identical. Sage provides frameworks
  to create unique parents. We use here the most easy one: Inheriting from the
  class :class:`sage.structure.unique_representation.UniqueRepresentation` is
  enough. Making parents unique can be quite important for an efficient
  implementation, because the repeated creation of "the same" parent would
  take a lot of time.

- Fraction fields are only defined for integral domains. Hence, we raise an
  error if the given ring does not belong to the category of integral
  domains. This is our first use case of categories.

- Last, we add a method that returns the characteristic of the field. We don't
  go into details, but some automated tests that we study below implicitly
  rely on this method.

We see that our basic implementation correctly refuses a ring that is not an
integral domain::

    sage: MyFrac(ZZ['x'])
    NewFrac(Univariate Polynomial Ring in x over Integer Ring)
    sage: MyFrac(Integers(15))
    Traceback (most recent call last):
    ...
    ValueError: Ring of integers modulo 15 is no integral domain

.. NOTE::

    Inheritance from :class:`~sage.structure.unique_representation.UniqueRepresentation`
    automatically provides our class with pickling, preserving the unique
    parent condition. If we had defined the class in some external module or
    in an interactive session, pickling would work immediately.

    However, for making the following example work in Sage's doctesting
    framework, we need to assign our class as an attribute of the ``__main__``
    module, so that the class can be looked up during unpickling.

::

    sage: import __main__
    sage: __main__.MyFrac = MyFrac
    sage: loads(dumps(MyFrac(ZZ))) is MyFrac(ZZ)
    True

.. NOTE::

    In the following sections, we will successively add or change details of
    ``MyFrac``. Rather than giving a full class definition in each step, we
    define new versions of ``MyFrac`` by inheriting from the previously
    defined version of ``MyFrac``. We believe this will help the reader to
    focus on the single detail that is relevant in each section.

    The complete code can be found in the appendix.

The elements
------------

We use the base class :class:`sage.structure.element.FieldElement`. Note that
in the creation of field elements it is not tested that the given parent is a
field::

    sage: from sage.structure.element import FieldElement
    sage: FieldElement(ZZ)
    Generic element of a structure

Our toy implementation of fraction field elements is based on the following
considerations:

- A fraction field element is defined by numerator and denominator, which both
  need to be elements of the base. There should be methods returning numerator
  resp. denominator.

- The denominator must not be zero, and (provided that the base is an ordered
  ring) we can make it non-negative, without loss of generality. By default,
  the denominator is one.

- The string representation is returned by the single\--underscore method
  ``_repr_``. In order to make our fraction field elements distinguishable
  from those already present in Sage, we use a different string representation.

- Arithmetic is implemented in single\--underscore method ``_add_``, ``_mul_``,
  etc. **We do not override the default double underscore __add__, __mul__**,
  since otherwise, we could not use Sage's coercion model.

- Comparisons can be implemented using ``_richcmp_``.
  This automatically makes the relational operators like
  ``==`` and ``<`` work. Inside this method, you can use
  the ``richcmp`` functions and related tools provided by sage.

  Note that ``_richcmp_`` should be provided,
  since otherwise comparison does not work::

      sage: class Foo(sage.structure.element.Element):
      ....:  def __init__(self, parent, x):
      ....:      self.x = x
      ....:  def _repr_(self):
      ....:      return "<%s>" % self.x
      sage: a = Foo(ZZ, 1)
      sage: b = Foo(ZZ, 2)
      sage: a <= b
      Traceback (most recent call last):
      ...
      TypeError: '<=' not supported between instances of 'Foo' and 'Foo'

- In the single underscore methods, we can assume that
  *both arguments belong to the same parent*.
  This is one benefit of the coercion model.

- When constructing new elements as the result of arithmetic operations, we do
  not directly name our class, but we use ``self.__class__``. Later, this will
  come in handy.

This gives rise to the following code::

    sage: class MyElement(FieldElement):
    ....:     def __init__(self, parent,n,d=None):
    ....:         B = parent.base()
    ....:         if d is None:
    ....:             d = B.one()
    ....:         if n not in B or d not in B:
    ....:             raise ValueError("Numerator and denominator must be elements of %s"%B)
    ....:         # Numerator and denominator should not just be "in" B,
    ....:         # but should be defined as elements of B
    ....:         d = B(d)
    ....:         n = B(n)
    ....:         if d==0:
    ....:             raise ZeroDivisionError("The denominator must not be zero")
    ....:         if d<0:
    ....:             self.n = -n
    ....:             self.d = -d
    ....:         else:
    ....:             self.n = n
    ....:             self.d = d
    ....:         FieldElement.__init__(self,parent)
    ....:     def numerator(self):
    ....:         return self.n
    ....:     def denominator(self):
    ....:         return self.d
    ....:     def _repr_(self):
    ....:         return "(%s):(%s)"%(self.n,self.d)
    ....:     def _richcmp_(self, other, op):
    ....:         from sage.structure.richcmp import richcmp
    ....:         return richcmp(self.n*other.denominator(), other.numerator()*self.d, op)
    ....:     def _add_(self, other):
    ....:         C = self.__class__
    ....:         D = self.d*other.denominator()
    ....:         return C(self.parent(), self.n*other.denominator()+self.d*other.numerator(), D)
    ....:     def _sub_(self, other):
    ....:         C = self.__class__
    ....:         D = self.d*other.denominator()
    ....:         return C(self.parent(), self.n*other.denominator()-self.d*other.numerator(),D)
    ....:     def _mul_(self, other):
    ....:         C = self.__class__
    ....:         return C(self.parent(), self.n*other.numerator(), self.d*other.denominator())
    ....:     def _div_(self, other):
    ....:         C = self.__class__
    ....:         return C(self.parent(), self.n*other.denominator(), self.d*other.numerator())

.. end of output


Features and limitations of the basic implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Thanks to the single underscore methods, some basic arithmetics works, **if**
we stay inside a single parent structure::

    sage: P = MyFrac(ZZ)
    sage: a = MyElement(P, 3, 4)
    sage: b = MyElement(P, 1, 2)
    sage: a+b, a-b, a*b, a/b
    ((10):(8), (2):(8), (3):(8), (6):(4))
    sage: a-b == MyElement(P, 1, 4)
    True

There is a default implementation of element tests. We can already do
::

    sage: a in P
    True

.. end of output

since `a` is defined as an element of `P`. However, we cannot verify yet that
the integers are contained in the fraction field of the ring of integers. It
does not even give a wrong answer, but results in an error::

    sage: 1 in P
    Traceback (most recent call last):
    ...
    NotImplementedError: cannot construct elements of NewFrac(Integer Ring)

.. end of output

We will take care of this later.

Categories in Sage
==================

Sometimes the base classes do not reflect the mathematics: The set of `m\times
n` matrices over a field forms, in general, not more than a vector
space. Hence, this set (called :class:`~sage.matrix.matrix_space.MatrixSpace`)
is not implemented on top of :class:`sage.rings.ring.Ring`.  However, if
`m=n`, then the matrix space is an algebra, thus, is a ring.

From the point of view of Python base classes, both cases are the same::

    sage: MS1 = MatrixSpace(QQ,2,3)
    sage: isinstance(MS1, Ring)
    False
    sage: MS2 = MatrixSpace(QQ,2)
    sage: isinstance(MS2, Ring)
    False

.. end of output

Sage's category framework can differentiate the two cases::

    sage: Rings()
    Category of rings
    sage: MS1 in Rings()
    False
    sage: MS2 in Rings()
    True

.. end of output

And indeed, ``MS2`` has *more* methods than ``MS1``::

    sage: import inspect
    sage: L1 = len([s for s in dir(MS1) if inspect.ismethod(getattr(MS1,s,None))])
    sage: L2 = len([s for s in dir(MS2) if inspect.ismethod(getattr(MS2,s,None))])
    sage: L1 < L2
    True

This is because the class of ``MS2`` also inherits from the parent
class for algebras::

    sage: MS1.__class__.__bases__
    (<class 'sage.matrix.matrix_space.MatrixSpace'>,
    <class 'sage.categories.category.JoinCategory.parent_class'>)
    sage: MS2.__class__.__bases__
    (<class 'sage.matrix.matrix_space.MatrixSpace'>,
    <class 'sage.categories.category.JoinCategory.parent_class'>)

.. end of output

Below, we will explain how this can be taken advantage of.

It is no surprise that our parent `P` defined above knows that it belongs to
the category of fields, as it is derived from the base class of fields.

::

    sage: P.category()
    Category of fields

.. end of output

However, we could choose a smaller category, namely the category of quotient fields.

Why should one choose a category?
---------------------------------

One can provide **default methods** *for all objects* of a category, and
*for all elements* of such objects. Hence, the category framework is a way
to inherit useful stuff that is not present in the base classes.  These
default methods do not rely on implementation details, but on mathematical
concepts.

In addition, the categories define **test suites** for their objects and
elements\---see the last section. Hence, one also gets basic sanity tests for
free.


How does the  *category framework* work?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Abstract base classes for the objects ("parent_class") and the elements of
objects ("element_class") are provided by attributes of the category. During
initialisation of a parent, the class of the parent is *dynamically changed*
into a sub\--class of the category's parent class. Likewise, sub\--classes of
the category's element class are available for the creation of elements of the
parent, as explained below.

A dynamic change of classes does not work in Cython. Nevertheless, method
inheritance still works, by virtue of a ``__getattr__`` method.

.. NOTE::

    It is strongly recommended to use the category framework both in Python
    and in Cython.

Let us see whether there is any gain in choosing the category of quotient
fields instead of the category of fields::

    sage: QuotientFields().parent_class, QuotientFields().element_class
    (<class 'sage.categories.quotient_fields.QuotientFields.parent_class'>,
     <class 'sage.categories.quotient_fields.QuotientFields.element_class'>)
    sage: [p for p in dir(QuotientFields().parent_class) if p not in dir(Fields().parent_class)]
    []
    sage: [p for p in dir(QuotientFields().element_class) if p not in dir(Fields().element_class)]
    ['_derivative',
     'denominator',
     'derivative',
     'numerator',
     'partial_fraction_decomposition']

.. end of output

So, there is no immediate gain for our fraction fields, but additional methods
become available to our fraction field elements. Note that some of these
methods are place-holders: There is no default implementation, but it is
*required* (respectively is *optional*) to implement these methods::

    sage: QuotientFields().element_class.denominator
    <abstract method denominator at ...>
    sage: from sage.misc.abstract_method import abstract_methods_of_class
    sage: abstract_methods_of_class(QuotientFields().element_class)['optional']
    ['_add_', '_mul_']
    sage: abstract_methods_of_class(QuotientFields().element_class)['required']
    ['__bool__', 'denominator', 'numerator']

Hence, when implementing elements of a quotient field, it is *required* to
implement methods returning the denominator and the numerator, and a method
that tells whether the element is nonzero, and in addition, it is *optional*
(but certainly recommended) to provide some arithmetic methods. If one forgets
to implement the required methods, the test suites of the category framework
will complain\---see below.


Implementing the category framework for the parent
--------------------------------------------------

We simply need to declare the correct category by an optional argument of the
field constructor, where we provide the possibility to override the default
category::

    sage: from sage.categories.quotient_fields import QuotientFields
    sage: class MyFrac(MyFrac):
    ....:     def __init__(self, base, category=None):
    ....:         if base not in IntegralDomains():
    ....:             raise ValueError("%s is no integral domain" % base)
    ....:         Field.__init__(self, base, category=category or QuotientFields())

When constructing instances of ``MyFrac``, their class is dynamically changed
into a new class called ``MyFrac_with_category``. It is a common sub\--class of
``MyFrac`` and of the category's parent class::

    sage: P = MyFrac(ZZ)
    sage: type(P)
    <class '__main__.MyFrac_with_category'>
    sage: isinstance(P,MyFrac)
    True
    sage: isinstance(P,QuotientFields().parent_class)
    True

The fraction field `P` inherits additional methods. For example, the base
class :class:`~sage.rings.fields.Field` does not have a method ``sum``. But
`P` inherits such method from the category of commutative additive
monoids\---see
:meth:`~sage.categories.commutative_additive_monoids.CommutativeAdditiveMonoids.ParentMethods.sum`::

    sage: P.sum.__module__
    'sage.categories.additive_monoids'

.. end of output

We have seen above that we can add elements. Nevertheless, the ``sum`` method
does not work, yet::

    sage: a = MyElement(P, 3, 4)
    sage: b = MyElement(P, 1, 2)
    sage: c = MyElement(P, -1, 2)
    sage: P.sum([a, b, c])
    Traceback (most recent call last):
    ...
    NotImplementedError: cannot construct elements of NewFrac(Integer Ring)

.. end of output

The reason is that the ``sum`` method starts with the return value of
``P.zero()``, which defaults to ``P(0)``\---but the conversion of
integers into ``P`` is not implemented, yet.

Implementing the category framework for the elements
----------------------------------------------------

Similar to what we have seen for parents, a new class is dynamically created
that combines the element class of the parent's category with the class that
we have implemented above. However, the category framework is implemented in a
different way for elements than for parents:

- We provide the parent `P` (or its class) with an attribute called
  "``Element``", whose value is a class.
- The parent *automatically* obtains an attribute ``P.element_class``, that
  subclasses both ``P.Element`` and ``P.category().element_class``.

Hence, for providing our fraction fields with their own element classes, **we
just need to add a single line to our class**::

    sage: class MyFrac(MyFrac):
    ....:     Element = MyElement


.. end of output

This little change provides several benefits:

- We can now create elements by simply calling the parent::

      sage: P = MyFrac(ZZ)
      sage: P(1), P(2,3)
      ((1):(1), (2):(3))

- There is a method ``zero`` returning the expected result::

      sage: P.zero()
      (0):(1)

- The ``sum`` method mentioned above suddenly works::

      sage: a = MyElement(P, 9, 4)
      sage: b = MyElement(P, 1, 2)
      sage: c = MyElement(P, -1, 2)
      sage: P.sum([a,b,c])
      (36):(16)

- Exponentiation now works out of the box using the multiplication
  that we defined::

    sage: a^3
    (729):(64)


What did happen behind the scenes to make this work?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provided ``P.Element``, and thus obtain ``P.element_class``, which is a
*lazy attribute*.  It provides a *dynamic* class, which is a sub\--class of
both ``MyElement`` defined above and of ``P.category().element_class``::

    sage: P.__class__.element_class
    <sage.misc.lazy_attribute.lazy_attribute object at ...>
    sage: P.element_class
    <class '__main__.MyFrac_with_category.element_class'>
    sage: type(P.element_class)
    <class 'sage.structure.dynamic_class.DynamicInheritComparisonMetaclass'>
    sage: issubclass(P.element_class, MyElement)
    True
    sage: issubclass(P.element_class,P.category().element_class)
    True

.. end of output

The *default* ``__call__`` method of `P` passes the given arguments to
``P.element_class``, adding the argument ``parent=P``. This is why we are now
able to create elements by calling the parent.

In particular, these elements are instances of that new dynamic class::

    sage: type(P(2,3))
    <class '__main__.MyFrac_with_category.element_class'>

.. end of output

.. NOTE::

    *All* elements of `P` should use the element class. In order to make sure
    that this also holds for the result of arithmetic operations, we created
    them as instances of ``self.__class__`` in the arithmetic methods of
    ``MyElement``.

``P.zero()`` defaults to returning ``P(0)`` and thus returns an
instance of ``P.element_class``. Since ``P.sum([...])`` starts the summation with
``P.zero()`` and the class of the sum only depends on the first
summand, by our implementation, we have::

    sage: type(a)
    <class '__main__.MyElement'>
    sage: isinstance(a,P.element_class)
    False
    sage: type(P.sum([a,b,c]))
    <class '__main__.MyFrac_with_category.element_class'>

.. end of output

The method ``factor`` provided by ``P.category().element_class`` (see above)
simply works::

    sage: a; a.factor(); P(6,4).factor()
    (9):(4)
    2^-2 * 3^2
    2^-1 * 3

.. end of output

But that's surprising: The element `a` is just an instance of ``MyElement``,
but not of ``P.element_class``, and its class does not know about the factor
method.  In fact, this is due to a ``__getattr__`` method defined for
:class:`sage.structure.element.Element`.
::

    sage: hasattr(type(a), 'factor')
    False
    sage: hasattr(P.element_class, 'factor')
    True
    sage: hasattr(a, 'factor')
    True

.. end of output

A first note on performance
---------------------------

The category framework is sometimes blamed for speed regressions, as in
:issue:`9138` and :issue:`11900`. But if the category framework is *used
properly*, then it is fast. For illustration, we determine the time needed to
access an attribute inherited from the element class. First, we consider an
element that uses the class that we implemented above, but does not use the
category framework properly::

    sage: type(a)
    <class '__main__.MyElement'>
    sage: timeit('a.factor',number=1000)     # random
    1000 loops, best of 3: 2 us per loop

.. end of output

Now, we consider an element that is equal to `a`, but uses the category
framework properly::

    sage: a2 = P(9,4)
    sage: a2 == a
    True
    sage: type(a2)
    <class '__main__.MyFrac_with_category.element_class'>
    sage: timeit('a2.factor',number=1000)    # random
    1000 loops, best of 3: 365 ns per loop

.. end of output

So,  *don't be afraid of using categories!*


Coercion\---the basics
======================

Theoretical background
----------------------

Coercion is not just *type conversion*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"Coercion" in the C programming language means "automatic type
conversion". However, in Sage, coercion is involved if one wants to be able to
do arithmetic, comparisons, etc. between elements of distinct parents. Hence,
**coercion is not about a change of types, but about a change of parents.**

As an illustration, we show that elements of the same type may very well belong
to rather different parents::

    sage: P1 = QQ['v,w']; P2 = ZZ['w,v']
    sage: type(P1.gen()) == type(P2.gen())
    True
    sage: P1 == P2
    False

.. end of output

`P_2` naturally is a sub\--ring of `P_1`. So, it makes sense to be able to add
elements of the two rings\---the result should then live in `P_1`, and indeed
it does::

    sage: (P1.gen()+P2.gen()).parent() is P1
    True

.. end of output

It would be rather inconvenient if one needed to *manually* convert an element
of `P_2` into `P_1` before adding. The coercion system does that conversion
automatically.

Not every conversion is a coercion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A coercion happens implicitly, without being explicitly requested by the
user. Hence, coercion must be based on mathematical rigour. In our example,
any element of `P_2` can be naturally interpreted as an element of `P_1`. We
thus have::

    sage: P1.has_coerce_map_from(P2)
    True
    sage: P1.coerce_map_from(P2)
    Coercion map:
      From: Multivariate Polynomial Ring in w, v over Integer Ring
      To:   Multivariate Polynomial Ring in v, w over Rational Field

While there is a conversion from `P_1` to `P_2` (namely restricted to
polynomials with integral coefficients), this conversion is not a coercion::

    sage: P2.convert_map_from(P1)
    Conversion map:
      From: Multivariate Polynomial Ring in v, w over Rational Field
      To:   Multivariate Polynomial Ring in w, v over Integer Ring
    sage: P2.has_coerce_map_from(P1)
    False
    sage: P2.coerce_map_from(P1) is None
    True

.. end of output

The four axioms requested for coercions
.......................................

1. A coercion is a morphism in an appropriate category.

   This first axiom has two implications:

   A. A coercion is defined on all elements of a parent.

      A polynomial of degree zero over the integers can be interpreted as an
      integer\---but the attempt to convert a polynomial of non-zero degree
      would result in an error::

          sage: ZZ(P2.one())
          1
          sage: ZZ(P2.gen(1))
          Traceback (most recent call last):
          ...
          TypeError: v is not a constant polynomial

      Hence, we only have a *partial* map. This is fine for a *conversion*,
      but a partial map does not qualify as a *coercion*.

   B. Coercions are structure preserving.

      Any real number can be converted to an integer, namely by
      rounding. However, such a conversion is not useful in arithmetic
      operations, since the underlying algebraic structure is not preserved::

          sage: int(1.6)+int(2.7) == int(1.6+2.7)
          False

      .. end of output

      The structure that is to be preserved depends on the category of the
      involved parents. For example, the coercion from the integers into the
      rational field is a homomorphism of euclidean domains::

          sage: QQ.coerce_map_from(ZZ).category_for()
          Join of Category of euclidean domains
          and Category of noetherian rings
          and Category of infinite sets
          and Category of metric spaces

      .. end of output

2. There is at most one coercion from one parent to another

   In addition, if there is a *coercion* from `P_2` to `P_1`, then a
   *conversion* from `P_2` to `P_1` is defined for all elements of `P_2` and
   coincides with the coercion.
   Nonetheless, user-exposed maps are copies of the internally used maps whence
   the lack of identity between different instantiations::

       sage: P1.coerce_map_from(P2) is P1.convert_map_from(P2)
       False

   For internally used maps, the maps are identical::

       sage: P1._internal_coerce_map_from(P2) is P1._internal_convert_map_from(P2)
       True

   .. end of output

3. Coercions can be composed

   If there is a coercion `\varphi: P_1 \to P_2` and another coercion `\psi:
   P_2 \to P_3`, then the composition of `\varphi` followed by `\psi` must
   yield the unique coercion from `P_1` to `P_3`.

4. The identity is a coercion

   Together with the two preceding axioms, it follows: If there are coercions
   from `P_1` to `P_2` and from `P_2` to `P_1`, then they are mutually
   inverse.


Implementing a conversion
-------------------------

We have seen above that some conversions into our fraction fields became
available after providing the attribute ``Element``.  However, we cannot
convert elements of a fraction field into elements of another fraction field,
yet::

    sage: P(2/3)
    Traceback (most recent call last):
    ...
    ValueError: Numerator and denominator must be elements of Integer Ring

.. end of output

For implementing a conversion, **the default __call__ method should (almost)
never be overridden.** Instead, **we implement the method
_element_constructor_**, that should return an instance of the parent's
element class.  Some old parent classes violate that rule\---please help to
refactor them!
::

    sage: class MyFrac(MyFrac):
    ....:     def _element_constructor_(self, *args, **kwds):
    ....:         if len(args)!=1:
    ....:             return self.element_class(self, *args, **kwds)
    ....:         x = args[0]
    ....:         try:
    ....:             P = x.parent()
    ....:         except AttributeError:
    ....:             return self.element_class(self, x, **kwds)
    ....:         if P in QuotientFields() and P != self.base():
    ....:             return self.element_class(self, x.numerator(), x.denominator(), **kwds)
    ....:         return self.element_class(self, x, **kwds)


.. end of output

In addition to the conversion from the base ring and from pairs of base ring
elements, we now also have a conversion from the rationals to our fraction
field of `\ZZ`:


::

    sage: P = MyFrac(ZZ)
    sage: P(2); P(2,3); P(3/4)
    (2):(1)
    (2):(3)
    (3):(4)

.. end of output

Recall that above, the test `1 \in P` failed with an error. We try again and
find that the error has disappeared. This is because we are now able to
convert the integer `1` into `P`. But the containment test still yields a
wrong answer::

    sage: 1 in P
    False

.. end of output

The technical reason: We have a conversion `P(1)` of `1` into `P`, but this is
not known as a coercion\---yet!
::

    sage: P.has_coerce_map_from(ZZ), P.has_coerce_map_from(QQ)
    (False, False)

.. end of output

Establishing a coercion
-----------------------

There are two main ways to make Sage use a particular conversion as a
coercion:

- One can use :meth:`sage.structure.parent.Parent.register_coercion`, normally
  during initialisation of the parent (see documentation of the method).
- A more flexible way is to provide a method ``_coerce_map_from_`` for the
  parent.

Let `P` and `R` be parents. If ``P._coerce_map_from_(R)`` returns ``False``
or ``None``, then there is no coercion from `R` to `P`. If it returns a map
with domain `R` and codomain `P`, then this map is used for coercion. If it
returns ``True``, then the conversion from `R` to `P` is used as coercion.

Note that in the following implementation, we need a special case for the
rational field, since ``QQ.base()`` is not the ring of integers.
::

    sage: class MyFrac(MyFrac):
    ....:     def _coerce_map_from_(self, S):
    ....:         if self.base().has_coerce_map_from(S):
    ....:             return True
    ....:         if S in QuotientFields():
    ....:             if self.base().has_coerce_map_from(S.base()):
    ....:                 return True
    ....:             if hasattr(S,'ring_of_integers') and self.base().has_coerce_map_from(S.ring_of_integers()):
    ....:                 return True


.. end of output

By the method above, a parent coercing into the base ring will also coerce
into the fraction field, and a fraction field coerces into another fraction
field if there is a coercion of the corresponding base rings. Now, we have::

    sage: P = MyFrac(QQ['x'])
    sage: P.has_coerce_map_from(ZZ['x']), P.has_coerce_map_from(Frac(ZZ['x'])), P.has_coerce_map_from(QQ)
    (True, True, True)

.. end of output

We can now use coercion from `\ZZ[x]` and from `\QQ` into `P` for arithmetic
operations between the two rings::

    sage: 3/4+P(2)+ZZ['x'].gen(), (P(2)+ZZ['x'].gen()).parent() is P
    ((4*x + 11):(4), True)

.. end of output

Equality and element containment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Recall that above, the test `1 \in P` gave a wrong answer. Let us repeat the
test now::

    sage: 1 in P
    True

.. end of output

Why is that?

The default element containment test `x \in P` is based on the interplay of
three building blocks: conversion, coercion, and equality test.

#. Clearly, if the conversion `P(x)` raises an error, then `x` cannot be seen as an element of `P`. On the other hand, a conversion `P(x)` can generally do very nasty things. So, the fact that `P(x)` works without error is necessary, but not sufficient for `x \in P`.
#. If `P` is the parent of `x`, then the conversion `P(x)` will not change `x` (at least, that's the default). Hence, we will have `x=P(x)`.
#. Sage uses coercion not only for arithmetic operations, but also for comparison: *If* there is a coercion from the parent of `x` to `P`, then the equality test ``x==P(x)`` reduces to ``P(x)==P(x)``. Otherwise, ``x==P(x)`` will evaluate as false.

That leads to the following default implementation of element containment testing:

.. NOTE::

    `x \in P` holds if and only if the test ``x==P(x)`` does not
    raise an error and evaluates as true.

If the user is not happy with that behaviour, the "magical" Python method
``__contains__`` can be overridden.

Coercion\---the advanced parts
==============================

So far, we are able to add integers and rational numbers to elements of our
new implementation of the fraction field of `\ZZ`.

::

    sage: P = MyFrac(ZZ)


.. end of output

::

    sage: 1/2+P(2,3)+1
    (13):(6)

.. end of output

Surprisingly, we can even add a polynomial over the integers to an element of
`P`, even though the *result lives in a new parent*, namely in a polynomial
ring over `P`::

    sage: P(1/2) + ZZ['x'].gen(), (P(1/2) + ZZ['x'].gen()).parent() is P['x']
    ((1):(1)*x + (1):(2), True)

.. end of output

In the next, seemingly more easy example, there "obviously" is a coercion from
the fraction field of `\ZZ` to the fraction field of `\ZZ[x]`.  However, Sage
does not know enough about our new implementation of fraction fields. Hence,
it does not recognise the coercion::

    sage: Frac(ZZ['x']).has_coerce_map_from(P)
    False

.. end of output

Two obvious questions arise:

#. How / why has the new ring been constructed in the example above?
#. How can we establish a coercion from `P`  to  `\mathrm{Frac}(\ZZ[x])`?

The key to answering both question is the construction of parents from simpler
pieces, that we are studying now. Note that we will answer the second question
*not* by providing a coercion from `P`  to  `\mathrm{Frac}(\ZZ[x])`, but by
teaching Sage to automatically construct `\mathrm{MyFrac}(\ZZ[x])` and coerce
both `P` and `\mathrm{Frac}(\ZZ[x])` into it.

If we are lucky, a parent can tell how it has been constructed::

    sage: Poly,R = QQ['x'].construction()
    sage: Poly,R
    (Poly[x], Rational Field)
    sage: Fract,R = QQ.construction()
    sage: Fract,R
    (FractionField, Integer Ring)

In both cases, the first value returned by
:meth:`~sage.structure.parent.Parent.construction` is a mathematical
construction, called *construction functor*\---see
:class:`~sage.categories.pushout.ConstructionFunctor`. The second return value
is a simpler parent to which the construction functor is applied.

Being functors, the same construction can be applied to different objects of a
category::

    sage: Poly(QQ) is QQ['x']
    True
    sage: Poly(ZZ) is ZZ['x']
    True
    sage: Poly(P) is P['x']
    True
    sage: Fract(QQ['x'])
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

Let us see on which categories these construction functors are defined::

    sage: Poly.domain()
    Category of rings
    sage: Poly.codomain()
    Category of rings
    sage: Fract.domain()
    Category of integral domains
    sage: Fract.codomain()
    Category of fields

In particular, the construction functors can be composed::

    sage: Poly*Fract
    Poly[x](FractionField(...))
    sage: (Poly*Fract)(ZZ) is QQ['x']
    True

.. end of output

In addition, it is often assumed that we have a coercion from input to output of the
construction functor::

    sage: ((Poly*Fract)(ZZ)).coerce_map_from(ZZ)
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

.. end of output

Construction functors do not necessarily commute::

    sage: (Fract*Poly)(ZZ)
    Fraction Field of Univariate Polynomial Ring in x over Integer Ring

.. end of output


The pushout of construction functors
------------------------------------

We can now formulate our problem. We have parents `P_1`, `P_2` and `R`, and
construction functors `F_1`, `F_2`, such that `P_1 = F_1(R)` and `P_2 =
F_2(R)`. We want to find a new construction functor `F_3`, such that both
`P_1` and `P_2` coerce into `P_3 = F_3(R)`.

In analogy to a notion of category theory, `P_3` is called the
:func:`~sage.categories.pushout.pushout` of `P_1` and `P_2`; and similarly
`F_3` is called the pushout of `F_1` and `F_2`.
::

    sage: from sage.categories.pushout import pushout
    sage: pushout(Fract(ZZ),Poly(ZZ))
    Univariate Polynomial Ring in x over Rational Field

.. end of output

`F_1\circ F_2` and `F_2\circ F_1` are natural candidates for the pushout of
`F_1` and `F_2`. However, the order of the functors must rely on a canonical
choice. "Indecomposable" construction functors have a *rank*, and this allows
to order them canonically:

.. NOTE::

    If ``F1.rank`` is smaller than ``F2.rank``, then the pushout is `F_2\circ
    F_1` (hence, `F_1` is applied first).

We have
::

    sage: Fract.rank, Poly.rank
    (5, 9)

.. end of output

and thus the pushout is
::

    sage: Fract.pushout(Poly), Poly.pushout(Fract)
    (Poly[x](FractionField(...)), Poly[x](FractionField(...)))

.. end of output

This is why the example above has worked.

However, only "elementary" construction functors have a rank::

    sage: (Fract*Poly).rank
    Traceback (most recent call last):
    ...
    AttributeError: 'CompositeConstructionFunctor' object has no attribute 'rank'...

.. end of output

Shuffling composite construction functors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If composed construction fuctors `...\circ F_2\circ F_1` and `...\circ
G_2\circ G_1` are given, then Sage determines their pushout by *shuffling* the
constituents:

- If ``F1.rank < G1.rank`` then we apply `F_1` first, and continue with `...\circ F_3\circ F_2` and `...\circ G_2\circ G_1`.
- If ``F1.rank > G1.rank`` then we apply `G_1` first, and continue with `...\circ F_2\circ F_1` and `...\circ G_3\circ G_2`.

If ``F1.rank == G1.rank``, then the tie needs to be broken by other techniques
(see below).

As an illustration, we first get us some functors and then see how chains of
functors are shuffled.
::

    sage: AlgClos, R = CC.construction(); AlgClos
    AlgebraicClosureFunctor

.. end of output

::

    sage: Compl, R = RR.construction(); Compl
    Completion[+Infinity, prec=53]

.. end of output

::

    sage: Matr, R = (MatrixSpace(ZZ,3)).construction(); Matr
    MatrixFunctor

.. end of output

::

    sage: AlgClos.rank, Compl.rank, Fract.rank, Poly.rank, Matr.rank
    (3, 4, 5, 9, 10)

.. end of output

When we apply ``Fract``, ``AlgClos``, ``Poly`` and ``Fract`` to the ring of
integers, we obtain::

    sage: (Fract*Poly*AlgClos*Fract)(ZZ)
    Fraction Field of Univariate Polynomial Ring in x over Algebraic Field

.. end of output

When we apply ``Compl``, ``Matr`` and ``Poly`` to the ring of integers, we
obtain::

    sage: (Poly*Matr*Compl)(ZZ)
    Univariate Polynomial Ring in x over Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

.. end of output

Applying the shuffling procedure yields
::

    sage: (Poly*Matr*Fract*Poly*AlgClos*Fract*Compl)(ZZ)
    Univariate Polynomial Ring in x over Full MatrixSpace of 3 by 3 dense matrices over Fraction Field of Univariate Polynomial Ring in x over Algebraic Field

.. end of output

and this is indeed equal to the pushout found by Sage::

    sage: pushout((Fract*Poly*AlgClos*Fract)(ZZ), (Poly*Matr*Compl)(ZZ))
    Univariate Polynomial Ring in x over Full MatrixSpace of 3 by 3 dense matrices over Fraction Field of Univariate Polynomial Ring in x over Algebraic Field

.. end of output

Breaking the tie
^^^^^^^^^^^^^^^^

If ``F1.rank==G1.rank`` then Sage's pushout constructions offers two ways to
proceed:

#. Construction functors have a method :meth:`~sage.categories.pushout.ConstructionFunctor.merge` that either returns ``None`` or returns a construction functor\---see below. If either ``F1.merge(G1)`` or ``G1.merge(F1)`` returns a construction functor `H_1`, then we apply `H_1` and continue with `...\circ F_3\circ F_2` and `...\circ G_3\circ G_2`.
#. Construction functors have a method :meth:`~sage.categories.pushout.ConstructionFunctor.commutes`. If either ``F1.commutes(G1)`` or ``G1.commutes(F1)`` returns ``True``, then we apply both `F_1` and `G_1` in any order, and continue with `...\circ F_3\circ F_2` and `...\circ G_3\circ G_2`.

By default, ``F1.merge(G1)`` returns ``F1`` if ``F1==G1``, and returns
``None`` otherwise. The ``commutes()`` method exists, but it seems that so far
nobody has implemented two functors of the same rank that commute.

Establishing a default implementation
-------------------------------------

The typical application of
:meth:`~sage.categories.pushout.ConstructionFunctor.merge` is to provide a
coercion between *different implementations* of the *same algebraic
structure*.

.. NOTE::

    If ``F1(P)`` and ``F2(P)`` are different implementations of the same thing, then ``F1.merge(F2)(P)`` should return the default implementation.

We want to boldly turn our toy implementation of fraction fields into the new
default implementation. Hence:

- Next, we implement a new version of the "usual" fraction field functor, having the same rank, but returning our new implementation.
- We make our new implementation the default, by virtue of a merge method.
- Since our fraction fields accept an optional argument ``category``, we pass
  the optional arguments to the construction functor, which will in turn use
  it to create a fraction field.

.. WARNING::

  - Do not override the default ``__call__`` method of :class:`~sage.categories.pushout.ConstructionFunctor`\---implement ``_apply_functor`` instead.
  - Declare domain and codomain of the functor during initialisation.

::

    sage: from sage.categories.pushout import ConstructionFunctor
    sage: class MyFracFunctor(ConstructionFunctor):
    ....:     rank = 5
    ....:     def __init__(self, args=None, kwds=None):
    ....:         self.args = args or ()
    ....:         self.kwds = kwds or {}
    ....:         ConstructionFunctor.__init__(self, IntegralDomains(), Fields())
    ....:     def _apply_functor(self, R):
    ....:         return MyFrac(R,*self.args,**self.kwds)
    ....:     def merge(self, other):
    ....:         if isinstance(other, (type(self), sage.categories.pushout.FractionField)):
    ....:             return self


.. end of output

::

    sage: MyFracFunctor()
    MyFracFunctor

.. end of output

We verify that our functor can really be used to construct our implementation of fraction fields, and that it can be merged with either itself or the usual fraction field constructor:


::

    sage: MyFracFunctor()(ZZ)
    NewFrac(Integer Ring)

.. end of output

::

    sage: MyFracFunctor().merge(MyFracFunctor())
    MyFracFunctor

.. end of output

::

    sage: MyFracFunctor().merge(Fract)
    MyFracFunctor

.. end of output

There remains to let our new fraction fields know about the new construction
functor. The arguments that were used when creating the fraction field are
stored as an attribute---this is a feature provided by
:class:`~sage.structure.unique_representation.CachedRepresentation`. We pass
all but the first of these arguments to the construction functor, such that
the construction functor is able to reconstruct the fraction field.

::

    sage: class MyFrac(MyFrac):
    ....:     def construction(self):
    ....:         return MyFracFunctor(self._reduction[1][1:], self._reduction[2]), self.base()


.. end of output

::

    sage: MyFrac(ZZ['x']).construction()
    (MyFracFunctor, Univariate Polynomial Ring in x over Integer Ring)

.. end of output

Due to merging, we have:


::

    sage: pushout(MyFrac(ZZ['x']), Frac(QQ['x']))
    NewFrac(Univariate Polynomial Ring in x over Rational Field)

.. end of output

A second note on performance
----------------------------

Being able to do arithmetics involving elements of different parents, with the
automatic creation of a pushout to contain the result, is certainly
convenient\---but one should not rely on it, if speed matters. Simply the
conversion of elements into different parents takes time. Moreover, by
:issue:`14058`, the pushout may be subject to Python's cyclic garbage
collection. Hence, if one does not keep a strong reference to it, the same
parent may be created repeatedly, which is a waste of time. In the following
example, we illustrate the slow\--down resulting from blindly relying on
coercion::

    sage: ZZxy = ZZ['x','y']
    sage: a = ZZxy('x')
    sage: b = 1/2
    sage: timeit("c = a+b")    # random
    10000 loops, best of 3: 172 us per loop
    sage: QQxy = QQ['x','y']
    sage: timeit("c2 = QQxy(a)+QQxy(b)") # random
    10000 loops, best of 3: 168 us per loop
    sage: a2 = QQxy(a)
    sage: b2 = QQxy(b)
    sage: timeit("c2 = a2+b2") # random
    100000 loops, best of 3: 10.5 us per loop

Hence, if one avoids the explicit or implicit conversion into the pushout, but
works in the pushout right away, one can get a more than 10\--fold speed\--up.

The test suites of the category framework
=========================================

The category framework does not only provide functionality but also a test
framework. This section logically belongs to the section on categories, but
without the bits that we have implemented in the section on coercion, our
implementation of fraction fields would not have passed the tests yet.

"Abstract" methods
------------------

We have already seen above that a category can require/suggest certain parent
or element methods, that the user must/should implement. This is in order to
smoothly blend with the methods that already exist in Sage.

The methods that ought to be provided are called
:func:`~sage.misc.abstract_method.abstract_method`. Let us see what methods
are needed for quotient fields and their elements::

    sage: from sage.misc.abstract_method import abstract_methods_of_class

.. end of output

::

    sage: abstract_methods_of_class(QuotientFields().parent_class)['optional']
    []
    sage: abstract_methods_of_class(QuotientFields().parent_class)['required']
    ['__contains__']

.. end of output

Hence, the only required method (that is actually required for all parents
that belong to the category of sets) is an element containment test. That's
fine, because the base class :class:`~sage.structure.parent.Parent` provides a
default containment test.

The elements have to provide more::

    sage: abstract_methods_of_class(QuotientFields().element_class)['optional']
    ['_add_', '_mul_']
    sage: abstract_methods_of_class(QuotientFields().element_class)['required']
    ['__bool__', 'denominator', 'numerator']

.. end of output

Hence, the elements must provide ``denominator()`` and ``numerator()``
methods, and must be able to tell whether they are zero or not. The base class
:class:`~sage.structure.element.Element` provides a default ``__bool__()``
method. In addition, the elements may provide Sage's single underscore
arithmetic methods (actually any ring element *should* provide them).

The ``_test_...`` methods
-------------------------

If a parent or element method's name start with "_test_", it gives rise to a
test in the automatic test suite. For example, it is tested

- whether a parent `P` actually is an instance of the parent class of the category of `P`,
- whether the user has implemented the required abstract methods,
- whether some defining structural properties (e.g., commutativity) hold.

For example, if one forgets to implement required methods, one obtains the
following error::

    sage: class Foo(Parent):
    ....:  Element = sage.structure.element.Element
    ....:  def __init__(self):
    ....:      Parent.__init__(self, category=QuotientFields())
    sage: Bar = Foo()
    sage: bar = Bar.element_class(Bar)
    sage: bar._test_not_implemented_methods()
    Traceback (most recent call last):
    ...
    AssertionError: Not implemented method: denominator

Here are the tests that form the test suite of quotient fields::

    sage: [t for t in dir(QuotientFields().parent_class) if t.startswith('_test_')]
    ['_test_additive_associativity',
     '_test_an_element',
     '_test_associativity',
     '_test_cardinality',
     '_test_characteristic',
     '_test_characteristic_fields',
     '_test_construction',
     '_test_distributivity',
     '_test_divides',
     '_test_elements',
     '_test_elements_eq_reflexive',
     '_test_elements_eq_symmetric',
     '_test_elements_eq_transitive',
     '_test_elements_neq',
     '_test_euclidean_degree',
     '_test_fraction_field',
     '_test_gcd_vs_xgcd',
     '_test_one',
     '_test_prod',
     '_test_quo_rem',
     '_test_some_elements',
     '_test_zero',
     '_test_zero_divisors']

.. end of output

We have implemented all abstract methods (or inherit them from base classes),
we use the category framework, and we have implemented coercions. So, we are
confident that the test suite runs without an error. In fact, it does!

.. NOTE::

    The following trick with the ``__main__`` module is only needed in
    doctests, not in an interactive session or when defining the classes
    externally.

::

    sage: __main__.MyFrac = MyFrac
    sage: __main__.MyElement = MyElement
    sage: P = MyFrac(ZZ['x'])
    sage: TestSuite(P).run()

.. end of output

Let us see what tests are actually performed::

    sage: TestSuite(P).run(verbose=True)
    running ._test_additive_associativity() . . . pass
    running ._test_an_element() . . . pass
    running ._test_associativity() . . . pass
    running ._test_cardinality() . . . pass
    running ._test_category() . . . pass
    running ._test_characteristic() . . . pass
    running ._test_characteristic_fields() . . . pass
    running ._test_construction() . . . pass
    running ._test_distributivity() . . . pass
    running ._test_divides() . . . pass
    running ._test_elements() . . .
      Running the test suite of self.an_element()
      running ._test_category() . . . pass
      running ._test_eq() . . . pass
      running ._test_new() . . . pass
      running ._test_nonzero_equal() . . . pass
      running ._test_not_implemented_methods() . . . pass
      running ._test_pickling() . . . pass
      pass
    running ._test_elements_eq_reflexive() . . . pass
    running ._test_elements_eq_symmetric() . . . pass
    running ._test_elements_eq_transitive() . . . pass
    running ._test_elements_neq() . . . pass
    running ._test_eq() . . . pass
    running ._test_euclidean_degree() . . . pass
    running ._test_fraction_field() . . . pass
    running ._test_gcd_vs_xgcd() . . . pass
    running ._test_new() . . . pass
    running ._test_not_implemented_methods() . . . pass
    running ._test_one() . . . pass
    running ._test_pickling() . . . pass
    running ._test_prod() . . . pass
    running ._test_quo_rem() . . . pass
    running ._test_some_elements() . . . pass
    running ._test_zero() . . . pass
    running ._test_zero_divisors() . . . pass

.. end of output

Implementing a new category with additional tests
-------------------------------------------------

As one can see, tests are also performed on elements. There are methods that
return one element or a list of some elements, relying on "typical" elements
that can be found in most algebraic structures.
::

    sage: P.an_element(); P.some_elements()
    (2):(1)
    [(2):(1)]

.. end of output

Unfortunately, the list of elements that is returned by the default method is
of length one, and that single element could also be a bit more interesting.
The method an_element relies on a method ``_an_element_()``, so, we implement
that. We also override the some_elements method.
::

    sage: class MyFrac(MyFrac):
    ....:     def _an_element_(self):
    ....:         a = self.base().an_element()
    ....:         b = self.base_ring().an_element()
    ....:         if (a+b)!=0:
    ....:             return self(a)**2/(self(a+b)**3)
    ....:         if b != 0:
    ....:             return self(a)/self(b)**2
    ....:         return self(a)**2*self(b)**3
    ....:     def some_elements(self):
    ....:         return [self.an_element(),self(self.base().an_element()),self(self.base_ring().an_element())]


.. end of output

::

    sage: P = MyFrac(ZZ['x'])
    sage: P.an_element(); P.some_elements()
    (x^2):(x^3 + 3*x^2 + 3*x + 1)
    [(x^2):(x^3 + 3*x^2 + 3*x + 1), (x):(1), (1):(1)]

.. end of output

Now, as we have more interesting elements, we may also add a test for the
"factor" method. Recall that the method was inherited from the category, but
it appears that it is not tested.

Normally, a test for a method defined by a category should be provided by the
same category. Hence, since ``factor`` is defined in the category of quotient
fields, a test should be added there. But we won't change source code here and
will instead create a sub\--category.

Apparently, If `e` is an element of a quotient field, the product of the
factors returned by ``e.factor()`` should be equal to `e`. For forming the
product, we use the ``prod`` method, that, no surprise, is inherited from
another category::

    sage: P.prod.__module__
    'sage.categories.monoids'

.. end of output

When we want to create a sub\--category, we need to provide a method
:meth:`~sage.categories.category.Category.super_categories`, that returns a
list of all immediate super categories (here: category of quotient fields).

.. WARNING::

    A sub\--category `S` of a category `C` is *not* implemented as a
    sub\--class of ``C.__class__``! `S` becomes a sub\--category of `C` only
    if ``S.super_categories()`` returns (a sub\--category of) `C`!

The parent and element methods of a category are provided as methods of
classes that are the attributes ``ParentMethods`` and ``Element Methods`` of
the category, as follows::

    sage: from sage.categories.category import Category
    sage: class QuotientFieldsWithTest(Category): # do *not* inherit from QuotientFields, but ...
    ....:     def super_categories(self):
    ....:         return [QuotientFields()]       # ... declare QuotientFields as a super category!
    ....:     class ParentMethods:
    ....:         pass
    ....:     class ElementMethods:
    ....:         def _test_factorisation(self, **options):
    ....:             P = self.parent()
    ....:             assert self == P.prod([P(b)**e for b,e in self.factor()])


.. end of output

We provide an instance of our quotient field implementation with that new
category. Note that categories have a default ``_repr_`` method, that guesses
a good string representation from the name of the class:
``QuotientFieldsWithTest`` becomes "quotient fields with test".

.. NOTE::

    The following trick with the ``__main__`` module is only needed in
    doctests, not in an interactive session or when defining the classes
    externally.

::

    sage: __main__.MyFrac = MyFrac
    sage: __main__.MyElement = MyElement
    sage: __main__.QuotientFieldsWithTest = QuotientFieldsWithTest
    sage: P = MyFrac(ZZ['x'], category=QuotientFieldsWithTest())
    sage: P.category()
    Category of quotient fields with test

.. end of output

The new test is inherited from the category. Since ``an_element()`` is returning a
complicated element, ``_test_factorisation`` is a serious test::

    sage: P.an_element()._test_factorisation
    <bound method QuotientFieldsWithTest.ElementMethods._test_factorisation of (x^2):(x^3 + 3*x^2 + 3*x + 1)>

.. end of output

::

    sage: P.an_element().factor()
    (x + 1)^-3 * x^2

.. end of output

Last, we observe that the new test has automatically become part of the test
suite. We remark that the existing tests became more serious as well, since we
made :meth:`sage.structure.parent.Parent.an_element` return something more
interesting.
::

    sage: TestSuite(P).run(verbose=True)
    running ._test_additive_associativity() . . . pass
    running ._test_an_element() . . . pass
    running ._test_associativity() . . . pass
    running ._test_cardinality() . . . pass
    running ._test_category() . . . pass
    running ._test_characteristic() . . . pass
    running ._test_characteristic_fields() . . . pass
    running ._test_construction() . . . pass
    running ._test_distributivity() . . . pass
    running ._test_divides() . . . pass
    running ._test_elements() . . .
      Running the test suite of self.an_element()
      running ._test_category() . . . pass
      running ._test_eq() . . . pass
      running ._test_factorisation() . . . pass
      running ._test_new() . . . pass
      running ._test_nonzero_equal() . . . pass
      running ._test_not_implemented_methods() . . . pass
      running ._test_pickling() . . . pass
      pass
    running ._test_elements_eq_reflexive() . . . pass
    running ._test_elements_eq_symmetric() . . . pass
    running ._test_elements_eq_transitive() . . . pass
    running ._test_elements_neq() . . . pass
    running ._test_eq() . . . pass
    running ._test_euclidean_degree() . . . pass
    running ._test_fraction_field() . . . pass
    running ._test_gcd_vs_xgcd() . . . pass
    running ._test_new() . . . pass
    running ._test_not_implemented_methods() . . . pass
    running ._test_one() . . . pass
    running ._test_pickling() . . . pass
    running ._test_prod() . . . pass
    running ._test_quo_rem() . . . pass
    running ._test_some_elements() . . . pass
    running ._test_zero() . . . pass
    running ._test_zero_divisors() . . . pass

.. end of output

Appendix: The complete code
===========================

.. highlight:: python
   :linenothreshold: 2

::

    # Importing base classes, ...
    import sage
    from sage.rings.ring import Field
    from sage.structure.element import FieldElement
    from sage.categories.category import Category
    # ... the UniqueRepresentation tool,
    from sage.structure.unique_representation import UniqueRepresentation
    # ... some categories, and ...
    from sage.categories.fields import Fields
    from sage.categories.quotient_fields import QuotientFields
    from sage.categories.integral_domains import IntegralDomains
    # construction functors
    from sage.categories.pushout import ConstructionFunctor

    # Fraction field elements
    class MyElement(FieldElement):
        def __init__(self, parent, n, d=None):
            if parent is None:
                raise ValueError("The parent must be provided")
            B = parent.base()
            if d is None:
                # The default denominator is one
                d = B.one()
            # verify that both numerator and denominator belong to the base
            if n not in B or d not in B:
                raise ValueError("Numerator and denominator must be elements of %s"%B)
            # Numerator and denominator should not just be "in" B,
            # but should be defined as elements of B
            d = B(d)
            n = B(n)
            # the denominator must not be zero
            if d==0:
                raise ZeroDivisionError("The denominator must not be zero")
            # normalize the denominator: WLOG, it shall be non-negative.
            if d<0:
                self.n = -n
                self.d = -d
            else:
                self.n = n
                self.d = d
            FieldElement.__init__(self,parent)

        # Methods required by the category of fraction fields:
        def numerator(self):
            return self.n
        def denominator(self):
            return self.d

        # String representation (single underscore!)
        def _repr_(self):
            return "(%s):(%s)"%(self.n,self.d)

        # Comparison: We can assume that both arguments are coerced
        # into the same parent, which is a fraction field. Hence, we
        # are allowed to use the denominator() and numerator() methods
        # on the second argument.
        def _richcmp_(self, other, op):
            from sage.structure.richcmp import richcmp
            return richcmp(self.n*other.denominator(), other.numerator()*self.d, op)

        # Arithmetic methods, single underscore. We can assume that both
        # arguments are coerced into the same parent.
        # We return instances of self.__class__, because self.__class__ will
        # eventually be a sub-class of MyElement.
        def _add_(self, other):
            C = self.__class__
            D = self.d*other.denominator()
            return C(self.parent(), self.n*other.denominator()+self.d*other.numerator(),D)
        def _sub_(self, other):
            C = self.__class__
            D = self.d*other.denominator()
            return C(self.parent(), self.n*other.denominator()-self.d*other.numerator(),D)
        def _mul_(self, other):
            C = self.__class__
            return C(self.parent(), self.n*other.numerator(), self.d*other.denominator())
        def _div_(self, other):
            C = self.__class__
            return C(self.parent(), self.n*other.denominator(), self.d*other.numerator())

    # Inheritance from UniqueRepresentation implements the unique parent
    # behaviour. Moreover, it implements pickling (provided that Python
    # succeeds to look up the class definition).
    class MyFrac(UniqueRepresentation, Field):
        # Implement the category framework for elements, which also
        # makes some basic conversions work.
        Element = MyElement

        # Allow to pass to a different category, by an optional argument
        def __init__(self, base, category=None):
            # Fraction fields only exist for integral domains
            if base not in IntegralDomains():
                raise ValueError("%s is no integral domain" % base)
            # Implement the category framework for the parent
            Field.__init__(self, base, category=category or QuotientFields())

        # Single-underscore method for string representation
        def _repr_(self):
            return "NewFrac(%s)"%repr(self.base())

        # Two methods that are implicitly used in some tests
        def base_ring(self):
            return self.base().base_ring()
        def characteristic(self):
            return self.base().characteristic()

        # Implement conversions. Do not override __call__!
        def _element_constructor_(self, *args, **kwds):
            if len(args)!=1:
               return self.element_class(self, *args, **kwds)
            x = args[0]
            try:
                P = x.parent()
            except AttributeError:
                return self.element_class(self, x, **kwds)
            if P in QuotientFields() and P != self.base():
                return self.element_class(self, x.numerator(), x.denominator(), **kwds)
            return self.element_class(self, x, **kwds)

        # Implement coercion from the base and from fraction fields
        # over a ring that coerces into the base
        def _coerce_map_from_(self, S):
            if self.base().has_coerce_map_from(S):
                return True
            if S in QuotientFields():
                if self.base().has_coerce_map_from(S.base()):
                    return True
                if hasattr(S,'ring_of_integers') and self.base().has_coerce_map_from(S.ring_of_integers()):
                    return True
        # Tell how this parent was constructed, in order to enable pushout constructions
        def construction(self):
            return MyFracFunctor(), self.base()

        # return some elements of this parent
        def _an_element_(self):
            a = self.base().an_element()
            b = self.base_ring().an_element()
            if (a+b)!=0:
                return self(a)**2/(self(a+b)**3)
            if b != 0:
                return self(a)/self(b)**2
            return self(a)**2*self(b)**3
        def some_elements(self):
            return [self.an_element(),self(self.base().an_element()),self(self.base_ring().an_element())]


    # A construction functor for our implementation of fraction fields
    class MyFracFunctor(ConstructionFunctor):
        # The rank is the same for Sage's original fraction field functor
        rank = 5
        def __init__(self):
            # The fraction field construction is a functor
            # from the category of integral domains into the category of
            # fields
            # NOTE: We could actually narrow the codomain and use the
            # category QuotientFields()
            ConstructionFunctor.__init__(self, IntegralDomains(), Fields())
        # Applying the functor to an object. Do not override __call__!
        def _apply_functor(self, R):
            return MyFrac(R)
        # Note: To apply the functor to morphisms, implement
        #       _apply_functor_to_morphism

        # Make sure that arithmetic involving elements of Frac(R) and
        # MyFrac(R) works and yields elements of MyFrac(R)
        def merge(self, other):
            if isinstance(other, (type(self), sage.categories.pushout.FractionField)):
                return self

    # A quotient field category with additional tests.
    # Notes:
    # - Category inherits from UniqueRepresentation. Hence, there
    #   is only one category for given arguments.
    # - Since QuotientFieldsWithTest is a singleton (there is only
    #   one instance of this class), we could inherit from
    #   sage.categories.category_singleton.Category_singleton
    #   rather than from sage.categories.category.Category
    class QuotientFieldsWithTest(Category):
        # Our category is a sub-category of the category of quotient fields,
        # by means of the following method.
        def super_categories(self):
            return [QuotientFields()]

        # Here, we could implement methods that are available for
        # all objects in this category.
        class ParentMethods:
            pass

        # Here, we add a new test that is available for all elements
        # of any object in this category.
        class ElementMethods:
            def _test_factorisation(self, **options):
                P = self.parent()
                # The methods prod() and factor() are inherited from
                # some other categories.
                assert self == P.prod([P(b)**e for b,e in self.factor()])


.. highlight:: python
   :linenothreshold: 22222
