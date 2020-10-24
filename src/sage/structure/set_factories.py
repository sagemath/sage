r"""
Set factories
=============


A *set factory* `F` is a device for constructing some :class:`Parent`
`P` that models subsets of a big set `S`. Typically, each such parent
is constructed as the subset of `S` of all elements satisfying a
certain collection of constraints `cons`. In such a hierarchy of
subsets, one needs an easy and flexible control on how elements are
constructed. For example, one may want to construct the elements of
`P` in some subclass of the class of the elements of `S`. On other
occasions, one also often needs `P` to be a facade parent, whose
elements are represented as elements of `S` (see
:class:`~sage.categories.facade_sets.FacadeSets`).

The role of a set factory is twofold:

- *Manage a database* of constructors for the different parents `P = F(cons)`
  depending on the various kinds of constraints `cons`. Note: currently there
  is no real support for that. We are gathering use cases before fixing the
  interface.

- Ensure that the elements `e = P(...)` created by the different parents
  follows a consistent policy concerning their *class and parent*.

.. RUBRIC:: Basic usage: constructing parents through a factory

The file :mod:`sage.structure.set_factories_example` shows an example of a
:class:`SetFactory` together with typical implementation. Note that the
written code is intentionally kept minimal, many things and in particular
several iterators could be written in a more efficient way.

Consider the set `S` of couples `(x,y)` with `x` and `y` in `I:=\{0,1,2,3,4\}`.
We represent an element of `S` as a 2-elements tuple, wrapped in a class
:class:`~.set_factories_example.XYPair` deriving from :class:`ElementWrapper`.
You can create a :class:`~.set_factories_example.XYPair` with any
:class:`Parent`::

    sage: from sage.structure.set_factories import *
    sage: from sage.structure.set_factories_example import *
    sage: p = XYPair(Parent(), (0,1)); p
    (0, 1)

Now, given `(a, b)\in S` we want to consider the following subsets of
`S`

.. MATH::

    S_a := \{(x,y) \in S \mid x = a\},

    S^b := \{(x,y) \in S \mid y = b\},

    S_a^b := \{(x,y) \in S \mid x = a, y = b\}.

The constraints considered here are admittedly trivial. In a realistic
example, there would be much more of them. And for some sets of constraints
no good enumeration algorithms would be known.

In Sage, those sets are constructed by passing the constraints to the
factory. We first create the set with no constraints at all::

    sage: XYPairs
    Factory for XY pairs
    sage: S = XYPairs(); S.list()
    [(0, 0), (1, 0), ..., (4, 0), (0, 1), (1, 1), ..., (3, 4), (4, 4)]
    sage: S.cardinality()
    25

Let us construct `S_2`, `S^3` and `S_2^3`::

    sage: Sx2 = XYPairs(x=2); Sx2.list()
    [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4)]
    sage: Sy3 = XYPairs(y=3); Sy3.list()
    [(0, 3), (1, 3), (2, 3), (3, 3), (4, 3)]
    sage: S23 = XYPairs(x=2, y=3); S23.list()
    [(2, 3)]

Set factories provide an alternative way to build subsets of an
already constructed set: each set constructed by a factory has a
method :meth:`~ParentWithSetFactory.subset` which accept new
constraints. Sets constructed by the factory or the
:meth:`~ParentWithSetFactory.subset` methods are identical::

    sage: Sx2s = S.subset(x=2); Sx2 is Sx2s
    True
    sage: Sx2.subset(y=3) is S23
    True

It is not possible to change an already given constraint::

    sage: S23.subset(y=5)
    Traceback (most recent call last):
    ...
    ValueError: Duplicate value for constraints 'y': was 3 now 5

.. RUBRIC:: Constructing custom elements: policies

We now come to the point of factories: constructing custom elements. The
writer of :func:`~.set_factories_example.XYPairs` decided that, by default,
the parents ``Sx2``, ``Sy3`` and ``S23`` are facade for parent ``S``. This
means that each element constructed by those subsets behaves as if they where
directly constructed by ``S`` itself::

    sage: Sx2.an_element().parent()
    AllPairs
    sage: el = Sx2.an_element()
    sage: el.parent() is S
    True
    sage: type(el) is S.element_class
    True

This is not always desirable. The device which decides how to construct an
element is called a *policy* (see :class:`SetFactoryPolicy`). Each factory
should have a default policy. Here is the policy of
:func:`~.set_factories_example.XYPairs`::

    sage: XYPairs._default_policy
    Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent AllPairs[=Factory for XY pairs(())]

This means that with the current policy, the parent builds elements with class
``XYPair`` and parent ``AllPairs`` which is itself constructed by calling the
factory :func:`~.set_factories_example.XYPairs` with constraints ``()``. There
is a lot of flexibility to change that. We now illustrate how to make a few
different choices.

1 - In a first use case, we want to add some methods to the constructed
elements. As illustration, we add here a new method ``sum`` which returns
`x+y`. We therefore create a new class for the elements which inherits from
:class:`~.set_factories_example.XYPair`::

    sage: class NewXYPair(XYPair):
    ....:     def sum(self):
    ....:         return sum(self.value)

Here is an instance of this class (with a dummy parent)::

    sage: el = NewXYPair(Parent(), (2,3))
    sage: el.sum()
    5

We now want to have subsets generating those new elements while still having a
single real parent (the one with no constraint) for each element. The
corresponding policy is called :class:`TopMostParentPolicy`. It takes three
parameters:

- the factory;
- the parameters for void constraint;
- the class used for elements.

Calling the factory with this policy returns a new set which builds its
elements with the new policy::

    sage: new_policy = TopMostParentPolicy(XYPairs, (), NewXYPair)
    sage: NewS = XYPairs(policy=new_policy)
    sage: el = NewS.an_element(); el
    (0, 0)
    sage: el.sum()
    0
    sage: el.parent() is NewS
    True
    sage: isinstance(el, NewXYPair)
    True

Newly constructed subsets inherit the policy::

    sage: NewS2 = NewS.subset(x=2)
    sage: el2 = NewS2.an_element(); el2
    (2, 0)
    sage: el2.sum()
    2
    sage: el2.parent() is NewS
    True

2 - In a second use case, we want the elements to remember which
parent created them. The corresponding policy is called
:class:`SelfParentPolicy`. It takes only two parameters:

- the factory;
- the class used for elements.

Here is an example::

    sage: selfpolicy = SelfParentPolicy(XYPairs, NewXYPair)
    sage: SelfS = XYPairs(policy=selfpolicy)
    sage: el = SelfS.an_element()
    sage: el.parent() is SelfS
    True

Now all subsets are the parent of the elements that they create::

    sage: SelfS2 = SelfS.subset(x=2)
    sage: el2 = SelfS2.an_element()
    sage: el2.parent() is NewS
    False
    sage: el2.parent() is SelfS2
    True

3 - Finally, a common use case is to construct simple python object which are
not Sage :class:`sage.structure.Element`. As an example, we show how to build
a parent ``TupleS`` which construct pairs as tuple. The corresponding policy
is called :class:`BareFunctionPolicy`. It takes two parameters:

- the factory;
- the function called to construct the elements.

Here is how to do it::

    sage: cons = lambda t, check: tuple(t) # ignore the check parameter
    sage: tuplepolicy = BareFunctionPolicy(XYPairs, cons)
    sage: P = XYPairs(x=2, policy=tuplepolicy)
    sage: P.list()
    [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4)]
    sage: el = P.an_element()
    sage: type(el)
    <... 'tuple'>

Here are the currently implemented policies:

- :class:`FacadeParentPolicy`: reuse an existing parent together with
  its element_class

- :class:`TopMostParentPolicy`: use a parent created by the factory
  itself and provide a class ``Element`` for it. In this case, we need
  to specify the set of constraints which build this parent taking the
  ownership of all elements and the class which will be used for the
  ``Element``.

- :class:`SelfParentPolicy`: provide systematically ``Element`` and
  element_class and ensure that the parent is ``self``.

- :class:`BareFunctionPolicy`: instead of calling a class constructor element
  are passed to a function provided to the policy.

.. TODO::

    Generalize :class:`TopMostParentPolicy` to be able to have several
    topmost parents.

.. RUBRIC:: Technicalities: how policies inform parents

Parents built from factories should inherit from
:class:`ParentWithSetFactory`. This class provide a few methods
related to factories and policies. The ``__init__`` method of
:class:`ParentWithSetFactory` must be provided with the set of
constraints and the policy. A parent built from a factory must create
elements through a call to the method ``_element_constructor_``. The
current way in which policies inform parents how to builds their elements is
set by a few attributes. So the class must accept attribute
adding. The precise details of which attributes are set may be subject
to change in the future.

.. RUBRIC:: How to write a set factory

.. SEEALSO:: :mod:`.set_factories_example` for an example of a factory.

Here are the specifications:

A parent built from a factory should

- *inherit* from :class:`ParentWithSetFactory`. It should accept a
  ``policy`` argument and pass it verbatim to the ``__init__`` method
  of :class:`ParentWithSetFactory` together with the set of
  constraints;

- *create its elements* through calls to the method
  ``_element_constructor_``;

- *define a method* :class:`ParentWithSetFactory.check_element` which
  checks if a built element indeed belongs to it. The method should
  accept an extra keyword parameter called ``check`` specifying which
  level of check should be performed. It will only be called when
  ``bool(check)`` evaluates to ``True``.

The constructor of the elements of a parent from a factory should:

- receive the parent as first mandatory argument;

- accept an extra optional keyword parameter called ``check`` which is meant
  to tell if the input must be checked or not. The precise meaning of
  ``check`` is intentionally left vague. The only intent is that if
  ``bool(check)`` evaluates to ``False``, no check is performed at all.

A factory should

- *define a method* ``__call__`` which is responsible for calling the
  appropriate parent constructor given the constraints;

- *define a method* overloading :meth:`SetFactory.add_constraints` which is
  responsible of computing the union of two sets of constraints;

- *optionally* define a method or an attribute ``_default_policy`` passed to
  the :class:`ParentWithSetFactory` if no policy is given to the factory.


.. TODO::

    There is currently no support for dealing with sets of
    constraints. The set factory and the parents must cooperate to
    consistently handle them. More support, together with a generic mechanism
    to select the appropriate parent class from the constraints, will be added
    as soon as we have gathered sufficiently enough use-cases.

AUTHORS:

- Florent Hivert (2011-2012): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2012 Florent Hivert <florent.hivert at lri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.abstract_method import abstract_method

####################################################
#                   Factories                      #
####################################################


class SetFactory(UniqueRepresentation, SageObject):
    r"""
    This class is currently just a stub that we will be using to add
    more structures on factories.

    TESTS::

        sage: from sage.structure.set_factories import SetFactory
        sage: S = SetFactory()
        sage: S.__call__("foo")
        Traceback (most recent call last):
        ...
        NotImplementedError: <abstract method __call__ at ...>
        sage: S.add_constraints("foo")
        Traceback (most recent call last):
        ...
        NotImplementedError: <abstract method add_constraints at ...>
    """
    @abstract_method
    def __call__(self, *constraints, **consdict):
        r"""
        Construct the parent associated with the constraints in
        argument. This should return a :class:`Parent`.

        .. NOTE::

            Currently there is no specification on how constraints are
            passed as arguments.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs()
            AllPairs
            sage: XYPairs(3)
            {(3, b) | b in range(5)}

            sage: XYPairs(x=3)
            {(3, b) | b in range(5)}

            sage: XYPairs(y=2)
            {(a, 2) | a in range(5)}

        TESTS::

            sage: from sage.structure.set_factories import SetFactory
            sage: F = SetFactory()
            sage: F()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method __call__ at 0x...>
        """

    @abstract_method
    def add_constraints(self, cons, *args, **opts):
        r"""
        Add constraints to the set of constraints `cons`.

        Should return a set of constraints.

        .. NOTE::

            Currently there is no specification on how constraints are
            passed as arguments.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs.add_constraints((3,),((None, 2), {}))
            (3, 2)

            sage: XYPairs.add_constraints((3,),((None, None), {'y': 2}))
            (3, 2)

        TESTS::

            sage: from sage.structure.set_factories import SetFactory
            sage: F = SetFactory()
            sage: F.add_constraints(())
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method add_constraints at 0x...>
        """

    # TODO : default policy ?

####################################################
#                    Policies                      #
####################################################


class SetFactoryPolicy(UniqueRepresentation, SageObject):
    r"""
    Abstract base class for policies.

    A policy is a device which is passed to a parent inheriting from
    :class:`ParentWithSetFactory` in order to set-up the element
    construction framework.

    INPUT:

    - ``factory`` -- a :class:`SetFactory`

    .. WARNING::

        This class is a base class for policies, one should not try
        to create instances.
    """
    def __init__(self, factory):
        r"""
        TESTS::

            sage: from sage.structure.set_factories import SetFactoryPolicy
            sage: from sage.structure.set_factories_example import XYPairs
            sage: S = SetFactoryPolicy(XYPairs); S
            <sage.structure.set_factories.SetFactoryPolicy object at ...>
        """
        assert isinstance(factory, SetFactory)
        self._factory = factory

    def factory(self):
        r"""
        Return the factory for ``self``.

        EXAMPLES::

            sage: from sage.structure.set_factories import SetFactoryPolicy, SelfParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: XYPairs._default_policy.factory()
            Factory for XY pairs
            sage: XYPairs._default_policy.factory() is XYPairs
            True

        TESTS::

            sage: policy = SetFactoryPolicy(XYPairs)
            sage: policy.factory()
            Factory for XY pairs
            sage: SelfParentPolicy(XYPairs, XYPair).factory()
            Factory for XY pairs
        """
        return self._factory

    def self_element_constructor_attributes(self, Element):
        r"""
        Element Constructor Attributes for non facade parent.

        The list of attributes which must be set during the init of a
        non facade parent with factory.

        INPUT:

        - ``Element`` -- the class used for the elements

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: pol = XYPairs._default_policy
            sage: pol.self_element_constructor_attributes(XYPair)
            {'Element': <class 'sage.structure.set_factories_example.XYPair'>,
             '_parent_for': 'self'}
        """
        return {'_parent_for': "self", 'Element': Element}

    def facade_element_constructor_attributes(self, parent):
        r"""
        Element Constructor Attributes for facade parent.

        The list of attributes which must be set during the init of a
        facade parent with factory.

        INPUT:

        - ``parent`` -- the actual parent for the elements

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: pol = XYPairs._default_policy
            sage: pol.facade_element_constructor_attributes(XYPairs())
            {'_facade_for': AllPairs,
             '_parent_for': AllPairs,
             'element_class': <class 'sage.structure.set_factories_example.AllPairs_with_category.element_class'>}
        """
        return {'_parent_for': parent,
                '_facade_for': parent,
                'element_class': parent.element_class}

    @abstract_method
    def element_constructor_attributes(self, constraints):
        r"""
        Element constructor attributes.

        INPUT:

        - ``constraints`` -- a bunch of constraints

        Should return the attributes that are prerequisite for element
        construction. This is coordinated with
        :meth:`ParentWithSetFactory._element_constructor_`. Currently two standard
        attributes are provided in
        :meth:`facade_element_constructor_attributes` and
        :meth:`self_element_constructor_attributes`. You should return the
        one needed depending on the given constraints.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: pol = XYPairs._default_policy
            sage: pol.element_constructor_attributes(())
            {'Element': <class 'sage.structure.set_factories_example.XYPair'>,
             '_parent_for': 'self'}
            sage: pol.element_constructor_attributes((1))
            {'_facade_for': AllPairs,
             '_parent_for': AllPairs,
             'element_class': <class 'sage.structure.set_factories_example.AllPairs_with_category.element_class'>}
        """


class SelfParentPolicy(SetFactoryPolicy):
    r"""
    Policy where each parent is a standard parent.

    INPUT:

    - ``factory`` -- an instance of :class:`SetFactory`
    - ``Element`` -- a subclass of :class:`~.element.Element`

    Given a factory ``F`` and a class ``E``, returns a policy for
    parent ``P`` creating elements in class ``E`` and parent ``P``
    itself.

    EXAMPLES::

        sage: from sage.structure.set_factories import SelfParentPolicy
        sage: from sage.structure.set_factories_example import XYPairs, XYPair, Pairs_Y
        sage: pol = SelfParentPolicy(XYPairs, XYPair)
        sage: S = Pairs_Y(3, pol)
        sage: el = S.an_element()
        sage: el.parent() is S
        True

        sage: class Foo(XYPair): pass
        sage: pol = SelfParentPolicy(XYPairs, Foo)
        sage: S = Pairs_Y(3, pol)
        sage: el = S.an_element()
        sage: isinstance(el, Foo)
        True
    """
    def __init__(self, factory, Element):
        r"""
        TESTS::

            sage: from sage.structure.set_factories import SelfParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: S = SelfParentPolicy(XYPairs, XYPair); S
            Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent ``self``
            sage: TestSuite(S).run(skip='_test_category')
        """
        self._Element = Element
        SetFactoryPolicy.__init__(self, factory)

    def element_constructor_attributes(self, constraints):
        r"""
        Return the element constructor attributes as per
        :meth:`SetFactoryPolicy.element_constructor_attributes`

        INPUT:

        - ``constraints`` -- a bunch of constraints

        TESTS::

            sage: from sage.structure.set_factories import SelfParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: pol = SelfParentPolicy(XYPairs, XYPair)
            sage: pol.element_constructor_attributes(())
            {'Element': <class 'sage.structure.set_factories_example.XYPair'>,
             '_parent_for': 'self'}
        """
        return self.self_element_constructor_attributes(self._Element)

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories import SelfParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: SelfParentPolicy(XYPairs, XYPair)    # indirect doctest
            Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent ``self``
        """
        return "Set factory policy for {} with parent ``self``".format(self._Element)


class TopMostParentPolicy(SetFactoryPolicy):
    r"""
    Policy where the parent of the elements is the topmost parent.

    INPUT:

    - ``factory`` -- an instance of :class:`SetFactory`
    - ``top_constraints`` -- the empty set of constraints.
    - ``Element`` -- a subclass of :class:`~.element.Element`

    Given a factory ``F`` and a class ``E``, returns a policy for
    parent ``P`` creating element in class ``E`` and parent
    ``factory(*top_constraints, policy)``.

    EXAMPLES::

        sage: from sage.structure.set_factories_example import XYPairs, XYPair
        sage: P = XYPairs(); P.policy()
        Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent AllPairs[=Factory for XY pairs(())]
    """
    def __init__(self, factory, top_constraints, Element):
        """
        TESTS::

            sage: from sage.structure.set_factories import TopMostParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: T = TopMostParentPolicy(XYPairs, (), XYPair); T
            Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent AllPairs[=Factory for XY pairs(())]
            sage: TestSuite(T).run(skip='_test_category')
        """
        # assert(isinstance(top_constraints, tuple))
        self._top_constraints = top_constraints
        self._Element = Element
        SetFactoryPolicy.__init__(self, factory)

    def element_constructor_attributes(self, constraints):
        r"""
        Return the element constructor attributes as per
        :meth:`SetFactoryPolicy.element_constructor_attributes`.

        INPUT:

        - ``constraints`` -- a bunch of constraints

        TESTS::

            sage: from sage.structure.set_factories import TopMostParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: pol = TopMostParentPolicy(XYPairs, (), XYPair)
            sage: pol.element_constructor_attributes(())
            {'Element': <class 'sage.structure.set_factories_example.XYPair'>,
             '_parent_for': 'self'}
            sage: pol.element_constructor_attributes((1))
            {'_facade_for': AllPairs,
             '_parent_for': AllPairs,
             'element_class': <class 'sage.structure.set_factories_example.AllPairs_with_category.element_class'>}
        """
        factory = self._factory
        if constraints == self._top_constraints:
            return self.self_element_constructor_attributes(self._Element)
        else:
            return self.facade_element_constructor_attributes(
                factory(*self._top_constraints, policy=self))

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories import TopMostParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: TopMostParentPolicy(XYPairs, (), XYPair)  # indirect doctest
            Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent AllPairs[=Factory for XY pairs(())]
        """
        return "Set factory policy for {} with parent {}[={}({})]".format(
            self._Element, self._factory(*self._top_constraints, policy=self),
            self._factory, self._top_constraints)


class FacadeParentPolicy(SetFactoryPolicy):
    r"""
    Policy for facade parent.

    INPUT:

    - ``factory`` -- an instance of :class:`SetFactory`
    - ``parent`` -- an instance of :class:`Parent`

    Given a factory ``F`` and a class ``E``, returns a policy for
    parent ``P`` creating elements as if they were created by
    ``parent``.

    EXAMPLES::

        sage: from sage.structure.set_factories import SelfParentPolicy, FacadeParentPolicy
        sage: from sage.structure.set_factories_example import XYPairs, XYPair

    We create a custom standard parent ``P``::

        sage: selfpolicy = SelfParentPolicy(XYPairs, XYPair)
        sage: P = XYPairs(x=2, policy=selfpolicy)
        sage: pol = FacadeParentPolicy(XYPairs, P)
        sage: P2 = XYPairs(x=2, y=3, policy=pol)
        sage: el = P2.an_element()
        sage: el.parent() is P
        True
        sage: type(el) is P.element_class
        True

    If ``parent`` is itself a facade parent, then transitivity is
    correctly applied::

        sage: P =  XYPairs()
        sage: P2 = XYPairs(x=2)
        sage: P2.category()
        Category of facade finite enumerated sets
        sage: pol = FacadeParentPolicy(XYPairs, P)
        sage: P23 = XYPairs(x=2, y=3, policy=pol)
        sage: el = P2.an_element()
        sage: el.parent() is P
        True
        sage: type(el) is P.element_class
        True
    """
    def __init__(self, factory, parent):
        r"""
        TESTS::

            sage: from sage.structure.set_factories import FacadeParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: F = FacadeParentPolicy(XYPairs, XYPairs()); F
            Set factory policy for facade parent AllPairs
            sage: TestSuite(F).run(skip='_test_category')
        """
        self._parent_for = parent
        SetFactoryPolicy.__init__(self, factory)

    def element_constructor_attributes(self, constraints):
        r"""
        Return the element constructor attributes as per
        :meth:`SetFactoryPolicy.element_constructor_attributes`.

        INPUT:

        - ``constraints`` -- a bunch of constraints

        TESTS::

            sage: from sage.structure.set_factories import FacadeParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: pol = FacadeParentPolicy(XYPairs, XYPairs())
            sage: pol.element_constructor_attributes(())
            {'_facade_for': AllPairs,
             '_parent_for': AllPairs,
             'element_class': <class 'sage.structure.set_factories_example.AllPairs_with_category.element_class'>}
            sage: pol.element_constructor_attributes((1))
            {'_facade_for': AllPairs,
             '_parent_for': AllPairs,
             'element_class': <class 'sage.structure.set_factories_example.AllPairs_with_category.element_class'>}
        """
        return self.facade_element_constructor_attributes(
            self._parent_for._parent_for)

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories import FacadeParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair
            sage: FacadeParentPolicy(XYPairs, XYPairs())  # indirect doctest
            Set factory policy for facade parent AllPairs
        """
        return "Set factory policy for facade parent {}".format(
            self._parent_for)




class BareFunctionPolicy(SetFactoryPolicy):
    r"""
    Policy where element are constructed using a bare function.

    INPUT:

    - ``factory`` -- an instance of :class:`SetFactory`
    - ``constructor`` -- a function

    Given a factory ``F`` and a function ``c``, returns a policy for
    parent ``P`` creating element using the function ``f``.

    EXAMPLES::

        sage: from sage.structure.set_factories import BareFunctionPolicy
        sage: from sage.structure.set_factories_example import XYPairs
        sage: cons = lambda t, check: tuple(t) # ignore the check parameter
        sage: tuplepolicy = BareFunctionPolicy(XYPairs, cons)
        sage: P = XYPairs(x=2, policy=tuplepolicy)
        sage: el = P.an_element()
        sage: type(el)
        <... 'tuple'>
    """
    def __init__(self, factory, constructor):
        """
        TESTS::

            sage: from sage.structure.set_factories import BareFunctionPolicy
            sage: from sage.structure.set_factories_example import XYPairs
            sage: pol = BareFunctionPolicy(XYPairs, tuple)
            sage: TestSuite(pol).run(skip='_test_category')
        """
        # assert(isinstance(top_constraints, tuple))
        self._constructor = constructor
        SetFactoryPolicy.__init__(self, factory)

    def element_constructor_attributes(self, constraints):
        r"""
        Return the element constructor attributes as per
        :meth:`SetFactoryPolicy.element_constructor_attributes`.

        INPUT:

        - ``constraints`` -- a bunch of constraints

        TESTS::

            sage: from sage.structure.set_factories import BareFunctionPolicy
            sage: from sage.structure.set_factories_example import XYPairs
            sage: pol = BareFunctionPolicy(XYPairs, tuple)
            sage: pol.element_constructor_attributes(())
            {'_element_constructor_': <... 'tuple'>, '_parent_for': None}
        """
        return {'_element_constructor_' : self._constructor, '_parent_for' : None}

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories import BareFunctionPolicy
            sage: from sage.structure.set_factories_example import XYPairs
            sage: BareFunctionPolicy(XYPairs, tuple)
            Set factory policy for bare function <... 'tuple'>
        """
        return "Set factory policy for bare function {}".format(self._constructor)


####################################################
#                     Parent                       #
####################################################


class ParentWithSetFactory(Parent):
    r"""
    Abstract class for parent belonging to a set factory.

    INPUT:

    - ``constraints`` -- a set of constraints
    - ``policy`` -- the policy for element construction
    - ``category`` -- the category of the parent (default to ``None``)

    Depending on the constraints and the policy, initialize the parent
    in a proper category to set up element construction.

    EXAMPLES::

        sage: from sage.structure.set_factories_example import XYPairs, PairsX_
        sage: P = PairsX_(3, XYPairs._default_policy)
        sage: P is XYPairs(3)
        True
        sage: P.category()
        Category of facade finite enumerated sets
    """
    def __init__(self, constraints, policy, category=None):
        r"""
        TESTS::

            sage: from sage.structure.set_factories import ParentWithSetFactory
            sage: from sage.structure.set_factories_example import XYPairs
            sage: isinstance(XYPairs(3), ParentWithSetFactory)  # indirect doctest
            True
        """
        self._constraints = constraints
        assert isinstance(policy, SetFactoryPolicy)
        self._policy = policy
        policy_attributes = policy.element_constructor_attributes(constraints)
        for attrname, attr in policy_attributes.items():
            if attr == "self":
                setattr(self, attrname, self)
            else:
                setattr(self, attrname, attr)
        assert self._parent_for is None or isinstance(self._parent_for, Parent)
        Parent.__init__(self,
                        category=category,
                        facade=policy_attributes.get('_facade_for', None))

    def constraints(self):
        r"""
        Return the constraints defining ``self``.

        .. NOTE::

            Currently there is no specification on how constraints are
            passed as arguments.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs().constraints()
            ()
            sage: XYPairs(x=3).constraints()
            (3, None)
            sage: XYPairs(y=2).constraints()
            (None, 2)
        """
        return self._constraints

    def policy(self):
        r"""
        Return the policy used when ``self`` was created.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs().policy()
            Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent AllPairs[=Factory for XY pairs(())]
            sage: XYPairs(x=3).policy()
            Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent AllPairs[=Factory for XY pairs(())]
        """
        return self._policy

    def facade_policy(self):
        r"""
        Return the policy for parent facade for ``self``.

        EXAMPLES::

            sage: from sage.structure.set_factories import SelfParentPolicy
            sage: from sage.structure.set_factories_example import XYPairs, XYPair

        We create a custom standard parent ``P``::

            sage: selfpolicy = SelfParentPolicy(XYPairs, XYPair)
            sage: P = XYPairs(x=2, policy=selfpolicy)
            sage: P.facade_policy()
            Set factory policy for facade parent {(2, b) | b in range(5)}

        Now passing ``P.facade_policy()`` creates parent which are facade for
        ``P``::

            sage: P3 = XYPairs(x=2, y=3, policy=P.facade_policy())
            sage: P3.facade_for() == (P,)
            True
            sage: el = P3.an_element()
            sage: el.parent() is P
            True
        """
        return FacadeParentPolicy(self.factory(), self)

    def factory(self):
        r"""
        Return the factory which built ``self``.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs().factory() is XYPairs
            True
            sage: XYPairs(x=3).factory() is XYPairs
            True
        """
        return self._policy.factory()

    def subset(self, *args, **options):
        r"""
        Return a subset of ``self`` by adding more constraints.

        .. NOTE::

            Currently there is no specification on how constraints are
            passed as arguments.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: S = XYPairs()
            sage: S3 = S.subset(x=3)
            sage: S3.list()
            [(3, 0), (3, 1), (3, 2), (3, 3), (3, 4)]

        TESTS::

            sage: S3 is XYPairs(3)
            True
            sage: S3 is XYPairs(x=3)
            True
        """
        factory = self.factory()
        constr = factory.add_constraints(self._constraints,
                                         (args, options))
        return factory(*constr, policy=self._policy)

    def _test_subset(self, **options):
        r"""
        Tests that subsets with no extra parameters returns
        ``self``.

        Currently, only the test that one gets the same parent when no
        more constraints are given, is performed.

        .. TODO::

            Straighten the test when handling of constraints will be
            specified.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: S = XYPairs()
            sage: S._test_subset()
        """
        tester = self._tester(**options)
        tester.assertIs(self.subset(), self)

    @abstract_method
    def check_element(self, x, check):
        r"""
        Check that ``x`` verifies the constraints of ``self``.

        INPUT:

        - ``x`` -- an instance of ``self.element_class``.

        - ``check`` -- the level of checking to be performed (usually a
          boolean).

        This method may assume that ``x`` was properly constructed by
        ``self`` or a possible super-set of ``self`` for which
        ``self`` is a facade. It should return nothing if ``x``
        verifies the constraints and raise a
        :exc:`~exceptions.ValueError` explaining which constraints
        ``x`` fails otherwise.

        The method should accept an extra parameter check specifying
        which level of check should be performed. It will only be
        called when ``bool(check)`` evaluates to ``True``.

        .. TODO::

            Should we always call check element and let it decide
            which check has to be performed ?

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: S = XYPairs()
            sage: el = S((2,3))
            sage: S.check_element(el, True)
            sage: XYPairs(x=2).check_element(el, True)
            sage: XYPairs(x=3).check_element(el, True)
            Traceback (most recent call last):
            ...
            ValueError: Wrong first coordinate
            sage: XYPairs(y=4).check_element(el, True)
            Traceback (most recent call last):
            ...
            ValueError: Wrong second coordinate
        """

    def __contains__(self, x):
        r"""
        Default implementation for ``__contains__``.

        INPUT:

        - ``x`` -- any object

        Check for class, parent and calls ``self.check_element(x)``.

        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: S = XYPairs()
            sage: el = S((2,3))
            sage: el in S
            True
            sage: el in XYPairs(x=2)
            True
            sage: el in XYPairs(x=3)
            False
            sage: el in XYPairs(y=4)
            False
        """
        if (isinstance(x, self.element_class) and
            x.parent() == self._parent_for):  # TODO: is_parent_of ???
            try:
                self.check_element(x, True)
            except ValueError:
                return False
            else:
                return True
        return False

    def __call__(self, *args, **keywords):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: S = XYPairs()
            sage: el = S((2,3)); el
            (2, 3)
            sage: S(el) is el
            True

            sage: XYPairs(x=3)((2,3))
            Traceback (most recent call last):
            ...
            ValueError: Wrong first coordinate

            sage: XYPairs(x=3)(el)
            Traceback (most recent call last):
            ...
            ValueError: Wrong first coordinate
        """
        # Ensure idempotence of element construction
        if (len(args) == 1 and
            isinstance(args[0], self.element_class) and
            args[0].parent() == self._parent_for):
            check = keywords.get("check", True)
            if check:
                self.check_element(args[0], check)
            return args[0]
        else:
            return Parent.__call__(self, *args, **keywords)

    # QUESTION: Should we call:
    #     self._parent_for._element_constructor_
    # currently we do not call it directly because:
    #  - it may do some extra check we dont want to perform ?
    #  - calling directly element_class should be faster
    def _element_constructor_(self, *args, **keywords):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs()((2,3)) # indirect doctest
            (2, 3)
            sage: XYPairs(x=3)((3,3)) # indirect doctest
            (3, 3)
            sage: XYPairs(x=3)((2,3)) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Wrong first coordinate

            sage: XYPairs(x=3)((2,3), check=False) # Don't do this at home, kids
            (2, 3)
        """
        check = keywords.get("check", True)
        res = self.element_class(self._parent_for, *args, **keywords)
        if check:
            self.check_element(res, check)
        return res
