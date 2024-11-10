r"""
Cartesian Products of Growth Groups

See :doc:`growth_group` for a description.

AUTHORS:

- Benjamin Hackl (2015)
- Daniel Krenn (2015)
- Clemens Heuberger (2016)

ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

- Benjamin Hackl is supported by the Google Summer of Code 2015.


TESTS::

    sage: from sage.rings.asymptotic.growth_group import GrowthGroup
    sage: A = GrowthGroup('(QQ_+)^x * x^ZZ'); A
    Growth Group QQ^x * x^ZZ
    sage: A.construction()
    (The cartesian_product functorial construction,
     (Growth Group QQ^x, Growth Group x^ZZ))
    sage: A.construction()[1][0].construction()
    (ExponentialGrowthGroup[x], Rational Field)
    sage: A.construction()[1][1].construction()
    (MonomialGrowthGroup[x], Integer Ring)
    sage: B = GrowthGroup('x^ZZ * y^ZZ'); B
    Growth Group x^ZZ * y^ZZ
    sage: B.construction()
    (The cartesian_product functorial construction,
     (Growth Group x^ZZ, Growth Group y^ZZ))
    sage: C = GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ'); C
    Growth Group x^ZZ * log(x)^ZZ * y^ZZ
    sage: C.construction()
    (The cartesian_product functorial construction,
     (Growth Group x^ZZ * log(x)^ZZ, Growth Group y^ZZ))
    sage: C.construction()[1][0].construction()
    (The cartesian_product functorial construction,
     (Growth Group x^ZZ, Growth Group log(x)^ZZ))
    sage: C.construction()[1][1].construction()
    (MonomialGrowthGroup[y], Integer Ring)

::

    sage: cm = sage.structure.element.get_coercion_model()
    sage: D = GrowthGroup('(QQ_+)^x * x^QQ')
    sage: cm.common_parent(A, D)
    Growth Group QQ^x * x^QQ
    sage: E = GrowthGroup('(ZZ_+)^x * x^QQ')
    sage: cm.record_exceptions()  # not tested, see #19411
    sage: cm.common_parent(A, E)
    Growth Group QQ^x * x^QQ
    sage: for t in cm.exception_stack():  # not tested, see #19411
    ....:     print(t)

::

    sage: assume(SR.an_element() > 0)
    sage: F = GrowthGroup('(SR_+)^n * n^ZZ * UU^n'); F
    Growth Group SR^n * n^ZZ * UU^n
    sage: G = GrowthGroup('QQ^n * n^QQ'); G
    Growth Group QQ^n * n^QQ * Signs^n
    sage: cm.common_parent(F, G)
    Growth Group SR^n * n^QQ * UU^n
    sage: forget()

::

    sage: A.an_element()
    (1/2)^x*x
    sage: tuple(E.an_element())
    (1, x^(1/2))

Classes and Methods
===================
"""

#*****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2014--2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.factory import UniqueFactory


class CartesianProductFactory(UniqueFactory):
    r"""
    Create various types of Cartesian products depending on its input.

    INPUT:

    - ``growth_groups`` -- tuple (or other iterable) of growth groups

    - ``order`` -- (default: ``None``) if specified, then this order
      is taken for comparing two Cartesian product elements. If ``order`` is
      ``None`` this is determined automatically.

    .. NOTE::

        The Cartesian product of growth groups is again a growth
        group. In particular, the resulting structure is partially
        ordered.

        The order on the product is determined as follows:

        - Cartesian factors with respect to the same variable are
          ordered lexicographically. This causes
          ``GrowthGroup('x^ZZ * log(x)^ZZ')`` and
          ``GrowthGroup('log(x)^ZZ * x^ZZ')`` to produce two
          different growth groups.

        - Factors over different variables are equipped with the
          product order (i.e. the comparison is component-wise).

        Also, note that the sets of variables of the Cartesian
        factors have to be either equal or disjoint.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: A = GrowthGroup('x^ZZ'); A
        Growth Group x^ZZ
        sage: B = GrowthGroup('log(x)^ZZ'); B
        Growth Group log(x)^ZZ
        sage: C = cartesian_product([A, B]); C  # indirect doctest
        Growth Group x^ZZ * log(x)^ZZ
        sage: C._le_ == C.le_lex
        True
        sage: D = GrowthGroup('y^ZZ'); D
        Growth Group y^ZZ
        sage: E = cartesian_product([A, D]); E  # indirect doctest
        Growth Group x^ZZ * y^ZZ
        sage: E._le_ == E.le_product
        True
        sage: F = cartesian_product([C, D]); F  # indirect doctest
        Growth Group x^ZZ * log(x)^ZZ * y^ZZ
        sage: F._le_ == F.le_product
        True
        sage: cartesian_product([A, E]); G  # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: The growth groups (Growth Group x^ZZ, Growth Group x^ZZ * y^ZZ)
        need to have pairwise disjoint or equal variables.
        sage: cartesian_product([A, B, D])  # indirect doctest
        Growth Group x^ZZ * log(x)^ZZ * y^ZZ

    TESTS::

        sage: from sage.rings.asymptotic.growth_group_cartesian import CartesianProductFactory
        sage: CartesianProductFactory('factory')([A, B], category=Groups() & Posets())
        Growth Group x^ZZ * log(x)^ZZ
        sage: CartesianProductFactory('factory')([], category=Sets())
        Traceback (most recent call last):
        ...
        TypeError: Cannot create Cartesian product without factors.

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G1 = GrowthGroup('x^QQ')
        sage: G2 = GrowthGroup('log(x)^ZZ')
        sage: G = cartesian_product([G1, G2])
        sage: cartesian_product([G1, G2], category=G.category()) is G
        True
    """
    def create_key_and_extra_args(self, growth_groups, category, **kwds):
        r"""
        Given the arguments and keywords, create a key that uniquely
        determines this object.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group_cartesian import CartesianProductFactory
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: A = GrowthGroup('x^ZZ')
            sage: CartesianProductFactory('factory').create_key_and_extra_args(
            ....:     [A], category=Sets(), order='blub')
            (((Growth Group x^ZZ,), Category of posets), {'order': 'blub'})
        """

        # CartesianProductPosets automatically add Posets() to their categories
        from sage.categories.category import Category
        from sage.categories.posets import Posets
        if not isinstance(category, tuple):
            category = (category,)
        category = Category.join(category + (Posets(),))

        return (tuple(growth_groups), category), kwds

    def create_object(self, version, args, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: cartesian_product([GrowthGroup('x^ZZ')])  # indirect doctest
            Growth Group x^ZZ
        """
        growth_groups, category = args
        if not growth_groups:
            raise TypeError('Cannot create Cartesian product without factors.')
        order = kwds.pop('order', None)
        if order is not None:
            return GenericProduct(growth_groups, category, order=order, **kwds)

        vg = tuple((g.variable_names(), g) for g in growth_groups)

        # check if all groups have a variable
        if not all(v for v, _ in vg):
            raise NotImplementedError('Growth groups %s have no variable.' %
                                      tuple(g for g in growth_groups
                                            if not g.variable_names()))

        # sort by variables
        from itertools import groupby, product
        vgs = tuple((v, tuple(gs)) for v, gs in
                    groupby(sorted(vg, key=lambda k: k[0]), key=lambda k: k[0]))

        # check whether variables are pairwise disjoint
        for u, w in product(iter(v for v, _ in vgs), repeat=2):
            if u != w and not set(u).isdisjoint(set(w)):
                raise ValueError('The growth groups %s need to have pairwise '
                                 'disjoint or equal variables.' % (growth_groups,))

        # build Cartesian products
        u_groups = list()
        for _, gs in vgs:
            gs = tuple(g for _, g in gs)
            if len(gs) > 1:
                u_groups.append(UnivariateProduct(gs, category, **kwds))
            else:
                u_groups.append(gs[0])

        if len(u_groups) > 1:
            m_group = MultivariateProduct(tuple(u_groups), category, **kwds)
        else:
            m_group = u_groups[0]
        return m_group


CartesianProductGrowthGroups = CartesianProductFactory('CartesianProductGrowthGroups')


from sage.combinat.posets.cartesian_product import CartesianProductPoset
from .growth_group import GenericGrowthGroup


class GenericProduct(CartesianProductPoset, GenericGrowthGroup):
    r"""
    A Cartesian product of growth groups.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: P = GrowthGroup('x^QQ')
        sage: L = GrowthGroup('log(x)^ZZ')
        sage: C = cartesian_product([P, L], order='lex'); C  # indirect doctest
        Growth Group x^QQ * log(x)^ZZ
        sage: C.an_element()
        x^(1/2)*log(x)

    ::

        sage: Px = GrowthGroup('x^QQ')
        sage: Lx = GrowthGroup('log(x)^ZZ')
        sage: Cx = cartesian_product([Px, Lx], order='lex')  # indirect doctest
        sage: Py = GrowthGroup('y^QQ')
        sage: C = cartesian_product([Cx, Py], order='product'); C  # indirect doctest
        Growth Group x^QQ * log(x)^ZZ * y^QQ
        sage: C.an_element()
        x^(1/2)*log(x)*y^(1/2)

    .. SEEALSO::

        :class:`~sage.sets.cartesian_product.CartesianProduct`,
        :class:`~sage.combinat.posets.cartesian_product.CartesianProductPoset`.
    """

    __classcall__ = CartesianProductPoset.__classcall__

    def __init__(self, sets, category, **kwds):
        r"""
        See :class:`GenericProduct` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ * y^ZZ')  # indirect doctest
            Growth Group x^ZZ * y^ZZ

        Check :issue:`26452`::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
            sage: R = QQ.extension(x^2+1, 'i')
            sage: P = MonomialGrowthGroup(R, 'w')
            sage: L = MonomialGrowthGroup(ZZ, 'log(w)')
            sage: cartesian_product([P, L])
            Growth Group w^(Number Field in i with defining polynomial x^2 + 1) * log(w)^ZZ
        """
        order = kwds.pop('order')
        CartesianProductPoset.__init__(self, sets, category, order, **kwds)

        vars = sum(iter(factor.variable_names()
                        for factor in self.cartesian_factors()),
                   tuple())
        from itertools import groupby
        from .growth_group import Variable
        Vars = Variable(tuple(v for v, _ in groupby(vars)), repr=self._repr_short_())

        GenericGrowthGroup.__init__(self, sets[0], Vars, self.category(), **kwds)

    __hash__ = CartesianProductPoset.__hash__

    def some_elements(self):
        r"""
        Return some elements of this Cartesian product of growth groups.

        See :class:`TestSuite` for a typical use case.

        OUTPUT: an iterator

        EXAMPLES::

            sage: from itertools import islice
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('(QQ_+)^y * x^QQ * log(x)^ZZ')
            sage: tuple(islice(G.some_elements(), 10r))
            (x^(1/2)*(1/2)^y,
             x^(-1/2)*log(x)*2^y,
             x^2*log(x)^(-1),
             x^(-2)*log(x)^2*42^y,
             log(x)^(-2)*(2/3)^y,
             x*log(x)^3*(3/2)^y,
             x^(-1)*log(x)^(-3)*(4/5)^y,
             x^42*log(x)^4*(5/4)^y,
             x^(2/3)*log(x)^(-4)*(6/7)^y,
             x^(-2/3)*log(x)^5*(7/6)^y)
        """
        return iter(
            self(c) for c in
            zip(*tuple(F.some_elements() for F in self.cartesian_factors())))

    def _create_element_in_extension_(self, element):
        r"""
        Create an element in an extension of this Cartesian product of
        growth groups which is chosen according to the input ``element``.

        INPUT:

        - ``element`` -- tuple

        OUTPUT: an element

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('z^ZZ * log(z)^ZZ')
            sage: z = G('z')[0]
            sage: lz = G('log(z)')[1]
            sage: G._create_element_in_extension_((z^3, lz)).parent()
            Growth Group z^ZZ * log(z)^ZZ
            sage: G._create_element_in_extension_((z^(1/2), lz)).parent()
            Growth Group z^QQ * log(z)^ZZ

        ::

            sage: G._create_element_in_extension_((3, 3, 3))
            Traceback (most recent call last):
            ...
            ValueError: Cannot create (3, 3, 3) as a Cartesian product like
            Growth Group z^ZZ * log(z)^ZZ.
        """
        factors = self.cartesian_factors()
        if len(element) != len(factors):
            raise ValueError('Cannot create %s as a Cartesian product like %s.' %
                             (element, self))

        if all(n.parent() is f for n, f in zip(element, factors)):
            parent = self
        else:
            parent = self._underlying_class()(tuple(n.parent() for n in element),
                                            category=self.category())
        return parent(element)

    def _element_constructor_(self, data):
        r"""
        Convert the given object to an element of this Cartesian product.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G_log = GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ')

        Conversion from the symbolic ring works::

            sage: x,y = var('x y')
            sage: G(x^-3*y^2)
            x^(-3)*y^2
            sage: G(x^4), G(y^2)
            (x^4, y^2)
            sage: G(1)
            1

        Even more complex expressions can be parsed::

            sage: G_log(x^42*log(x)^-42*y^42)
            x^42*log(x)^(-42)*y^42

        TESTS::

            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G('x'), G('y')
            (x, y)

        ::

            sage: G_log(log(x))
            log(x)

        ::

            sage: G(G.cartesian_factors()[0].gen())
            x

        ::

            sage: GrowthGroup('QQ^x * x^QQ')(['x^(1/2)'])
            x^(1/2)
            sage: l = GrowthGroup('x^ZZ * log(x)^ZZ')(['x', 'log(x)']); l
            x*log(x)
            sage: type(l)
            <class 'sage.rings.asymptotic.growth_group_cartesian.UnivariateProduct_with_category.element_class'>
            sage: GrowthGroup('(QQ_+)^x * x^QQ')(['2^log(x)'])
            Traceback (most recent call last):
            ...
            ValueError: ['2^log(x)'] is not in Growth Group QQ^x * x^QQ.
            > *previous* ValueError: 2^log(x) is not in any of the factors of
            Growth Group QQ^x * x^QQ
            >> *previous* ValueError: 2^log(x) is not in Growth Group QQ^x.
            >> *and* ValueError: 2^log(x) is not in Growth Group x^QQ.
            sage: GrowthGroup('(QQ_+)^x * x^QQ')(['2^log(x)', 'x^55'])
            Traceback (most recent call last):
            ...
            ValueError: ['2^log(x)', 'x^55'] is not in Growth Group QQ^x * x^QQ.
            > *previous* ValueError: 2^log(x) is not in any of the factors of
            Growth Group QQ^x * x^QQ
            >> *previous* ValueError: 2^log(x) is not in Growth Group QQ^x.
            >> *and* ValueError: 2^log(x) is not in Growth Group x^QQ.

        ::

            sage: n = GrowthGroup('n^ZZ * log(n)^ZZ')('n')
            sage: G = GrowthGroup('(QQ_+)^n * n^ZZ * log(n)^ZZ')
            sage: G(n).value
            (1, n, 1)
        """
        from sage.sets.cartesian_product import CartesianProduct
        from sage.symbolic.ring import SR

        def convert_factors(data, raw_data):
            try:
                return self._convert_factors_(data)
            except ValueError as e:
                from .misc import combine_exceptions
                raise combine_exceptions(
                    ValueError('%s is not in %s.' % (raw_data, self)), e)

        if data == 1:
            return self.one()

        elif data is None:
            raise ValueError('%s cannot be converted.' % (data,))

        elif type(data) is self.element_class and data.parent() == self:
            return data

        elif isinstance(data, str):
            from .misc import split_str_by_op
            return convert_factors(split_str_by_op(data, '*'), data)

        elif hasattr(data, 'parent'):
            P = data.parent()

            if P is self:
                return data

            elif P is SR:
                from sage.symbolic.operators import mul_vararg
                if data.operator() == mul_vararg:
                    return convert_factors(data.operands(), data)

            # room for other parents (e.g. polynomial ring et al.)

        try:
            return super()._element_constructor_(data)
        except (TypeError, ValueError):
            pass
        if isinstance(data, (tuple, list, CartesianProduct.Element)):
            return convert_factors(tuple(data), data)

        return convert_factors((data,), data)

    _repr_ = GenericGrowthGroup._repr_

    def _repr_short_(self):
        r"""
        A short (shorter than :meth:`._repr_`) representation string
        for this Cartesian product of growth groups.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^QQ')
            sage: L = GrowthGroup('log(x)^ZZ')
            sage: cartesian_product([P, L], order='lex')._repr_short_()
            'x^QQ * log(x)^ZZ'
        """
        return ' * '.join(S._repr_short_() for S in self.cartesian_factors())

    def _convert_factors_(self, factors):
        r"""
        Helper method. Try to convert some ``factors`` to an
        element of one of the Cartesian factors and return the product of
        all these factors.

        INPUT:

        - ``factors`` -- tuple or other iterable

        OUTPUT: an element of this Cartesian product

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * log(x)^QQ * y^QQ')
            sage: e1 = G._convert_factors_([x^2])
            sage: (e1, e1.parent())
            (x^2, Growth Group x^ZZ * log(x)^QQ * y^QQ)

        ::

            sage: G = GrowthGroup('(QQ_+)^n * n^ZZ * UU^n')
            sage: n = SR.var('n')
            sage: G((-2)^n)
            2^n*(-1)^n
        """
        from sage.misc.misc_c import prod
        from .growth_group import PartialConversionValueError
        from .misc import combine_exceptions

        def get_factors(data):
            result = []
            errors = []
            for factor in self.cartesian_factors():
                try:
                    try:
                        result.append((factor, factor(data)))
                        break
                    except PartialConversionValueError as e:
                        try:
                            element, todo = e.element.split()
                        except NotImplementedError as nie:
                            raise combine_exceptions(
                                ValueError('cannot split {}: no splitting '
                                           'implemented'.format(e.element)),
                                nie)
                        except ValueError as ve:
                            raise combine_exceptions(
                                ValueError('cannot split {} after failed '
                                           'conversion into element of '
                                           '{}'.format(e.element, factor)),
                                ve)
                        assert todo is not None
                        result.append((factor, element))
                        data = todo
                except (ValueError, TypeError) as error:
                    errors.append(error)
            if not result:
                raise combine_exceptions(
                    ValueError('%s is not in any of the factors of %s' % (data, self)),
                    *errors)
            return result

        return prod(self.cartesian_injection(*fs)
                    for f in factors for fs in get_factors(f))

    def cartesian_injection(self, factor, element):
        r"""
        Inject the given element into this Cartesian product at the given factor.

        INPUT:

        - ``factor`` -- a growth group (a factor of this Cartesian product)

        - ``element`` -- an element of ``factor``

        OUTPUT: an element of this Cartesian product

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * y^QQ')
            sage: G.cartesian_injection(G.cartesian_factors()[1], 'y^7')
            y^7
        """
        return self(tuple((f.one() if f != factor else element)
                          for f in self.cartesian_factors()))

    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` coerces into this growth group.

        INPUT:

        - ``S`` -- a parent

        OUTPUT: boolean

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: A = GrowthGroup('QQ^x * x^QQ')
            sage: B = GrowthGroup('QQ^x * x^ZZ')
            sage: A.has_coerce_map_from(B) # indirect doctest
            True
            sage: B.has_coerce_map_from(A) # indirect doctest
            False
        """
        if CartesianProductPoset.has_coerce_map_from(self, S):
            return True

        elif isinstance(S, GenericProduct):
            factors = S.cartesian_factors()
        else:
            factors = (S,)

        if all(any(g.has_coerce_map_from(f) for g in self.cartesian_factors())
               for f in factors):
            return True

    def _pushout_(self, other):
        r"""
        Construct the pushout of this and the other growth group. This is called by
        :func:`sage.categories.pushout.pushout`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.categories.pushout import pushout
            sage: cm = sage.structure.element.get_coercion_model()
            sage: A = GrowthGroup('(QQ_+)^x * x^ZZ')
            sage: B = GrowthGroup('x^ZZ * log(x)^ZZ')
            sage: A._pushout_(B)
            Growth Group QQ^x * x^ZZ * log(x)^ZZ
            sage: pushout(A, B)
            Growth Group QQ^x * x^ZZ * log(x)^ZZ
            sage: cm.discover_coercion(A, B)
            ((map internal to coercion system -- copy before use)
             Coercion map:
               From: Growth Group QQ^x * x^ZZ
               To:   Growth Group QQ^x * x^ZZ * log(x)^ZZ,
             (map internal to coercion system -- copy before use)
             Coercion map:
               From: Growth Group x^ZZ * log(x)^ZZ
               To:   Growth Group QQ^x * x^ZZ * log(x)^ZZ)
            sage: cm.common_parent(A, B)
            Growth Group QQ^x * x^ZZ * log(x)^ZZ

        ::

            sage: C = GrowthGroup('(QQ_+)^x * x^QQ * y^ZZ')
            sage: D = GrowthGroup('x^ZZ * log(x)^QQ * (QQ_+)^z')
            sage: C._pushout_(D)
            Growth Group QQ^x * x^QQ * log(x)^QQ * y^ZZ * QQ^z
            sage: cm.common_parent(C, D)
            Growth Group QQ^x * x^QQ * log(x)^QQ * y^ZZ * QQ^z
            sage: A._pushout_(D)
            Growth Group QQ^x * x^ZZ * log(x)^QQ * QQ^z
            sage: cm.common_parent(A, D)
            Growth Group QQ^x * x^ZZ * log(x)^QQ * QQ^z
            sage: cm.common_parent(B, D)
            Growth Group x^ZZ * log(x)^QQ * QQ^z
            sage: cm.common_parent(A, C)
            Growth Group QQ^x * x^QQ * y^ZZ
            sage: E = GrowthGroup('log(x)^ZZ * y^ZZ')
            sage: cm.common_parent(A, E)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Growth Group QQ^x * x^ZZ' and 'Growth Group log(x)^ZZ * y^ZZ'

        ::

            sage: F = GrowthGroup('z^QQ')
            sage: pushout(C, F)
            Growth Group QQ^x * x^QQ * y^ZZ * z^QQ

        ::

            sage: pushout(GrowthGroup('(QQ_+)^x * x^ZZ'), GrowthGroup('(ZZ_+)^x * x^QQ'))
            Growth Group QQ^x * x^QQ
            sage: cm.common_parent(GrowthGroup('(QQ_+)^x * x^ZZ'), GrowthGroup('(ZZ_+)^x * x^QQ'))
            Growth Group QQ^x * x^QQ

        ::

            sage: pushout(GrowthGroup('(QQ_+)^n * n^QQ'), GrowthGroup('(SR_+)^n'))
            Growth Group SR^n * n^QQ

        ::

            sage: cm.common_parent(GrowthGroup('n^ZZ * log(n)^ZZ * UU^n'),
            ....:                  GrowthGroup('n^QQ * UU^n'))
            Growth Group n^QQ * log(n)^ZZ * UU^n
        """
        from .growth_group import GenericGrowthGroup, AbstractGrowthGroupFunctor
        from .misc import bidirectional_merge_sorted
        from sage.structure.element import get_coercion_model

        Sfactors = self.cartesian_factors()
        if isinstance(other, GenericProduct):
            Ofactors = other.cartesian_factors()
        elif isinstance(other, GenericGrowthGroup):
            Ofactors = (other,)
        elif (other.construction() is not None and
              isinstance(other.construction()[0], AbstractGrowthGroupFunctor)):
            Ofactors = (other,)
        else:
            return

        def pushout_univariate_factors(self, other, var, Sfactors, Ofactors):
            try:
                return bidirectional_merge_sorted(
                    Sfactors, Ofactors,
                    lambda f: (f._underlying_class(), f._var_.var_repr))
            except RuntimeError:
                pass

            cm = get_coercion_model()
            try:
                Z = cm.common_parent(*Sfactors+Ofactors)
                return (Z,), (Z,)
            except TypeError:
                pass

            def subfactors(F):
                for f in F:
                    if isinstance(f, GenericProduct):
                        yield from subfactors(f.cartesian_factors())
                    else:
                        yield f

            try:
                return bidirectional_merge_sorted(
                    tuple(subfactors(Sfactors)), tuple(subfactors(Ofactors)),
                    lambda f: (f._underlying_class(), f._var_.var_repr))
            except RuntimeError:
                pass

            from sage.structure.coerce_exceptions import CoercionException
            raise CoercionException(
                'Cannot construct the pushout of %s and %s: The factors '
                'with variables %s are not overlapping, '
                'no common parent was found, and '
                'splitting the factors was unsuccessful.' % (self, other, var))

        # A wrapper around an iterator that stores additional intermediate data.
        # This deviates slightly from the iterator protocol:
        # At the end of the iteration the data is reset to None instead
        # of raising a StopIteration.
        class it:
            def __init__(self, it):
                self.it = it
                self.var = None
                self.factors = None

            def next_custom(self):
                try:
                    self.var, factors = next(self.it)
                    self.factors = tuple(factors)
                except StopIteration:
                    self.var = None
                    self.factors = tuple()

        from itertools import groupby
        S = it(groupby(Sfactors, key=lambda k: k.variable_names()))
        O = it(groupby(Ofactors, key=lambda k: k.variable_names()))

        newS = []
        newO = []

        S.next_custom()
        O.next_custom()
        while S.var is not None or O.var is not None:
            if S.var is not None and O.var is not None and S.var < O.var:
                newS.extend(S.factors)
                newO.extend(S.factors)
                S.next_custom()
            elif S.var is not None and O.var is not None and S.var > O.var:
                newS.extend(O.factors)
                newO.extend(O.factors)
                O.next_custom()
            else:
                SL, OL = pushout_univariate_factors(self, other, S.var,
                                                    S.factors, O.factors)
                newS.extend(SL)
                newO.extend(OL)
                S.next_custom()
                O.next_custom()

        assert (len(newS) == len(newO))

        if (len(Sfactors) == len(newS) and
            len(Ofactors) == len(newO)):
            # We had already all factors in each of self and
            # other, thus splitting it in subproblems (one for
            # each factor) is the strategy to use. If a pushout is
            # possible :func:`sage.categories.pushout.pushout`
            # will manage this by itself.
            return

        from sage.categories.pushout import pushout
        from sage.categories.cartesian_product import cartesian_product
        return pushout(cartesian_product(newS), cartesian_product(newO))

    def gens_monomial(self):
        r"""
        Return a tuple containing monomial generators of this growth group.

        OUTPUT: a tuple containing elements of this growth group

        .. NOTE::

            This method calls the ``gens_monomial()`` method on the
            individual factors of this Cartesian product and
            concatenates the respective outputs.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ * log(z)^ZZ')
            sage: G.gens_monomial()
            (x, y)

        TESTS::

            sage: all(g.parent() == G for g in G.gens_monomial())
            True
        """
        return sum(iter(
            tuple(self.cartesian_injection(factor, g) for g in factor.gens_monomial())
            for factor in self.cartesian_factors()),
                   tuple())

    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT: a tuple of strings

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ * log(z)^ZZ').variable_names()
            ('x', 'y', 'z')
        """
        vars = sum(iter(factor.variable_names()
                        for factor in self.cartesian_factors()),
                   tuple())
        from itertools import groupby
        return tuple(v for v, _ in groupby(vars))

    class Element(CartesianProductPoset.Element):

        from .growth_group import _is_lt_one_
        is_lt_one = _is_lt_one_

        def _repr_(self, latex=False):
            r"""
            A representation string for this Cartesian product element.

            INPUT:

            - ``latex`` -- boolean (default: ``False``); if set, then
              LaTeX-output is returned

            OUTPUT: string

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: P = GrowthGroup('x^QQ')
                sage: L = GrowthGroup('log(x)^ZZ')
                sage: cartesian_product([P, L], order='lex').an_element()._repr_()
                'x^(1/2)*log(x)'
            """
            if latex:
                from sage.misc.latex import latex as latex_repr
                f = latex_repr
            else:
                f = repr

            mul = ' ' if latex else '*'
            s = mul.join(f(v) for v in self.value if not v.is_one())
            if s == '':
                return '1'
            return s

        def _latex_(self):
            r"""
            A representation string for this Cartesian product element.

            OUTPUT: string

            TESTS::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: P = GrowthGroup('x^QQ')
                sage: L = GrowthGroup('log(x)^ZZ')
                sage: latex(cartesian_product([P, L], order='lex').an_element())  # indirect doctest
                x^{\frac{1}{2}} \log\left(x\right)
                sage: latex(GrowthGroup('(QQ_+)^n * n^QQ').an_element())  # indirect doctest
                \left(\frac{1}{2}\right)^{n} n^{\frac{1}{2}}
            """
            return self._repr_(latex=True)

        def __pow__(self, exponent):
            r"""
            Calculate the power of this growth element to the given
            ``exponent``.

            INPUT:

            - ``exponent`` -- a number

            OUTPUT: a growth element

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('x^ZZ * y^QQ * z^ZZ')
                sage: x, y, z = G.gens_monomial()
                sage: (x^5 * y * z^5)^(1/5)  # indirect doctest
                x*y^(1/5)*z

            ::

                sage: G = GrowthGroup('x^QQ * log(x)^QQ'); x = G('x')
                sage: (x^(21/5) * log(x)^7)^(1/42)  # indirect doctest
                x^(1/10)*log(x)^(1/6)
            """
            return self.parent()._create_element_in_extension_(
                tuple(x ** exponent for x in self.cartesian_factors()))

        def factors(self):
            r"""
            Return the atomic factors of this growth element. An atomic factor
            cannot be split further and is not the identity (`1`).

            OUTPUT: a tuple of growth elements

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ')
                sage: x, y = G.gens_monomial()
                sage: x.factors()
                (x,)
                sage: f = (x * y).factors(); f
                (x, y)
                sage: tuple(factor.parent() for factor in f)
                (Growth Group x^ZZ, Growth Group y^ZZ)
                sage: f = (x * log(x)).factors(); f
                (x, log(x))
                sage: tuple(factor.parent() for factor in f)
                (Growth Group x^ZZ, Growth Group log(x)^ZZ)

            ::

                sage: G = GrowthGroup('x^ZZ * log(x)^ZZ * log(log(x))^ZZ * y^QQ')
                sage: x, y = G.gens_monomial()
                sage: f = (x * log(x) * y).factors(); f
                (x, log(x), y)
                sage: tuple(factor.parent() for factor in f)
                (Growth Group x^ZZ, Growth Group log(x)^ZZ, Growth Group y^QQ)

            ::

                sage: G.one().factors()
                ()
            """
            return sum(iter(f.factors()
                            for f in self.cartesian_factors()
                            if not f.is_one()),
                       tuple())

        from .growth_group import _log_factor_, _log_
        log = _log_
        log_factor = _log_factor_

        def _log_factor_(self, base=None, locals=None):
            r"""
            Helper method for calculating the logarithm of the factorization
            of this element.

            INPUT:

            - ``base`` -- the base of the logarithm; if ``None``
              (default value) is used, the natural logarithm is taken

            - ``locals`` -- dictionary which may contain the following keys and
              values:

              - ``'log'`` -- value: a function; if not used, then the usual
                :class:`log <sage.functions.log.Function_log>` is taken

            OUTPUT:

            A tuple of pairs, where the first entry is either a growth
            element or something out of which we can construct a growth element
            and the second a multiplicative coefficient.

            TESTS::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('QQ^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ')
                sage: x, y = G.gens_monomial()
                sage: (x * y).log_factor()  # indirect doctest
                ((log(x), 1), (log(y), 1))
            """
            if self.is_one():
                return tuple()

            def try_create_growth(g):
                try:
                    return self.parent()(g)
                except (TypeError, ValueError):
                    return g

            try:
                return sum(iter(tuple((try_create_growth(g), c)
                                      for g, c in
                                      factor._log_factor_(base=base,
                                                          locals=locals))
                                for factor in self.cartesian_factors()
                                if factor != factor.parent().one()),
                           tuple())
            except (ArithmeticError, TypeError, ValueError) as e:
                from .misc import combine_exceptions
                raise combine_exceptions(
                    ArithmeticError('Cannot build log(%s) in %s.' %
                                    (self, self.parent())), e)

        from .growth_group import _rpow_
        rpow = _rpow_

        def _rpow_element_(self, base):
            r"""
            Return an element which is the power of ``base`` to this element.

            INPUT:

            - ``base`` -- an element

            OUTPUT: a growth element

            .. NOTE::

                The parent of the result can be different from the parent
                of this element.

            A :exc:`ValueError` is raised if the calculation is not possible
            within this method. (Then the calling method should take care
            of the calculation.)

            TESTS::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('QQ^x * x^ZZ * log(x)^ZZ')
                sage: lx = log(G('x'))
                sage: rp = lx._rpow_element_('e'); rp
                x
                sage: rp.parent()
                Growth Group x^ZZ
            """
            factors = self.factors()
            if len(factors) != 1:
                raise ValueError  # calling method has to deal with it...
            from .growth_group import MonomialGrowthGroup
            factor = factors[0]
            if not isinstance(factor.parent(), MonomialGrowthGroup):
                raise ValueError  # calling method has to deal with it...
            return factor._rpow_element_(base)

        def exp(self):
            r"""
            The exponential of this element.

            OUTPUT: a growth element

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('x^ZZ * log(x)^ZZ * log(log(x))^ZZ')
                sage: x = G('x')
                sage: exp(log(x))
                x
                sage: exp(log(log(x)))
                log(x)

            ::

                sage: exp(x)
                Traceback (most recent call last):
                ...
                ArithmeticError: Cannot construct e^x in
                Growth Group x^ZZ * log(x)^ZZ * log(log(x))^ZZ
                > *previous* TypeError: unsupported operand parent(s) for *:
                'Growth Group x^ZZ * log(x)^ZZ * log(log(x))^ZZ' and
                'Growth Group (e^x)^ZZ'

            TESTS::

                sage: E = GrowthGroup("(e^y)^QQ * y^QQ * log(y)^QQ")
                sage: y = E('y')
                sage: log(exp(y))
                y
                sage: exp(log(y))
                y
            """
            return self.rpow('e')

        def __invert__(self):
            r"""
            Return the multiplicative inverse of this Cartesian product.

            OUTPUT: a growth element

            .. NOTE::

                The result may live in a larger parent than we started with.

            TESTS::

                 sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                 sage: G = GrowthGroup('(ZZ_+)^x * x^ZZ')
                 sage: g = G('2^x * x^3')
                 sage: (~g).parent()
                 Growth Group QQ^x * x^ZZ
            """
            return self.parent()._create_element_in_extension_(
                tuple(~x for x in self.cartesian_factors()))

        def _substitute_(self, rules):
            r"""
            Substitute the given ``rules`` in this
            Cartesian product growth element.

            INPUT:

            - ``rules`` -- dictionary; the neutral element of the group is
              replaced by the value to key ``'_one_'``

            OUTPUT: an object

            TESTS::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('x^QQ * log(x)^QQ')
                sage: G(x^3 * log(x)^5)._substitute_({'x': SR.var('z')})
                z^3*log(z)^5
                sage: _.parent()
                Symbolic Ring
                sage: G(x^3 * log(x)^5)._substitute_({'x': 2.2})  # rel tol 1e-6
                3.24458458945
                sage: _.parent()
                Real Field with 53 bits of precision
                sage: G(1 / x)._substitute_({'x': 0})
                Traceback (most recent call last):
                ...
                ZeroDivisionError: Cannot substitute in x^(-1) in
                Growth Group x^QQ * log(x)^QQ.
                > *previous* ZeroDivisionError: Cannot substitute in x^(-1) in
                Growth Group x^QQ.
                >> *previous* ZeroDivisionError: rational division by zero
                sage: G(1)._substitute_({'_one_': 'one'})
                'one'
            """
            if self.is_one():
                return rules['_one_']
            from sage.symbolic.operators import mul_vararg
            try:
                return mul_vararg(
                    *tuple(x._substitute_(rules)
                           for x in self.cartesian_factors()))
            except (ArithmeticError, TypeError, ValueError) as e:
                from .misc import substitute_raise_exception
                substitute_raise_exception(self, e)

        def _singularity_analysis_(self, var, zeta, precision):
            r"""
            Perform singularity analysis on this growth element.

            INPUT:

            - ``var`` -- string denoting the variable

            - ``zeta`` -- a number

            - ``precision`` -- integer

            OUTPUT:

            An asymptotic expansion for `[z^n] f` where `n` is ``var``
            and `f` has this growth element as a singular expansion
            in `T=\frac{1}{1-\frac{z}{\zeta}}\to \infty` where this
            element is a growth element in `T`.

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('exp(x)^QQ * x^QQ * log(x)^QQ')
                sage: G(x^(1/2))._singularity_analysis_('n', 2, precision=2)
                1/sqrt(pi)*(1/2)^n*n^(-1/2) - 1/8/sqrt(pi)*(1/2)^n*n^(-3/2)
                + O((1/2)^n*n^(-5/2))
                sage: G(log(x))._singularity_analysis_('n', 1, precision=5)
                n^(-1) + O(n^(-3))
                sage: G(x*log(x))._singularity_analysis_('n', 1, precision=5)
                log(n) + euler_gamma + 1/2*n^(-1) + O(n^(-2))

            TESTS::

                sage: G('exp(x)*log(x)')._singularity_analysis_('n', 1, precision=5)
                Traceback (most recent call last):
                ...
                NotImplementedError: singularity analysis of exp(x)*log(x)
                not implemented
                sage: G('exp(x)*x*log(x)')._singularity_analysis_('n', 1, precision=5)
                Traceback (most recent call last):
                ...
                NotImplementedError: singularity analysis of exp(x)*x*log(x)
                not yet implemented since it has more than two factors
                sage: G(1)._singularity_analysis_('n', 2, precision=3)
                Traceback (most recent call last):
                ...
                NotImplementedOZero: got O(0)
                The error term O(0) means 0 for sufficiently large n.
                sage: G('exp(x)')._singularity_analysis_('n', 2, precision=3)
                Traceback (most recent call last):
                ...
                NotImplementedError: singularity analysis of exp(x)
                not implemented
            """
            factors = self.factors()
            if len(factors) == 0:
                from .asymptotic_expansion_generators import asymptotic_expansions
                from .misc import NotImplementedOZero
                raise NotImplementedOZero(var=var, exact_part=0)
            elif len(factors) == 1:
                return factors[0]._singularity_analysis_(
                    var=var, zeta=zeta, precision=precision)
            elif len(factors) == 2:
                from .growth_group import MonomialGrowthGroup
                from sage.rings.integer_ring import ZZ

                a, b = factors
                if all(isinstance(f.parent(), MonomialGrowthGroup)
                       for f in factors) \
                        and a.parent().gens_monomial() \
                        and b.parent().gens_logarithmic() \
                        and a.parent().variable_name() == \
                            b.parent().variable_name():
                    if b.exponent not in ZZ:
                        raise NotImplementedError(
                            'singularity analysis of {} not implemented '
                            'since exponent {} of {} is not an integer'.format(
                                self, b.exponent, b.parent().gen()))

                    from sage.rings.asymptotic.asymptotic_expansion_generators import \
                        asymptotic_expansions
                    return asymptotic_expansions.SingularityAnalysis(
                        var=var, zeta=zeta, alpha=a.exponent,
                        beta=ZZ(b.exponent), delta=0,
                        precision=precision, normalized=False)
                else:
                    raise NotImplementedError(
                        'singularity analysis of {} not implemented'.format(self))
            else:
                raise NotImplementedError(
                    'singularity analysis of {} not yet implemented '
                    'since it has more than two factors'.format(self))

        def variable_names(self):
            r"""
            Return the names of the variables of this growth element.

            OUTPUT: a tuple of strings

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('QQ^m * m^QQ * log(n)^ZZ')
                sage: G('2^m * m^4 * log(n)').variable_names()
                ('m', 'n')
                sage: G('2^m * m^4').variable_names()
                ('m',)
                sage: G('log(n)').variable_names()
                ('n',)
                sage: G('m^3').variable_names()
                ('m',)
                sage: G('m^0').variable_names()
                ()
            """
            vars = sum(iter(factor.variable_names()
                            for factor in self.factors()),
                       tuple())
            from itertools import groupby
            return tuple(v for v, _ in groupby(vars))

    CartesianProduct = CartesianProductGrowthGroups


class UnivariateProduct(GenericProduct):
    r"""
    A Cartesian product of growth groups with the same variables.

    .. NOTE::

        A univariate product of growth groups is ordered
        lexicographically. This is motivated by the assumption
        that univariate growth groups can be ordered in a chain
        with respect to the growth they model (e.g.
        ``x^ZZ * log(x)^ZZ``: polynomial growth dominates
        logarithmic growth).

    .. SEEALSO::

        :class:`MultivariateProduct`,
        :class:`GenericProduct`.
    """

    def __init__(self, sets, category, **kwargs):
        r"""
        See :class:`UnivariateProduct` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: type(GrowthGroup('x^ZZ * log(x)^ZZ'))  # indirect doctest
            <class 'sage.rings.asymptotic.growth_group_cartesian.UnivariateProduct_with_category'>
        """
        super().__init__(sets, category, order='lex', **kwargs)

    CartesianProduct = CartesianProductGrowthGroups


class MultivariateProduct(GenericProduct):
    r"""
    A Cartesian product of growth groups with pairwise disjoint
    (or equal) variable sets.

    .. NOTE::

        A multivariate product of growth groups is ordered by
        means of the product order, i.e. component-wise. This is
        motivated by the assumption that different variables are
        considered to be independent (e.g. ``x^ZZ * y^ZZ``).

    .. SEEALSO::

        :class:`UnivariateProduct`,
        :class:`GenericProduct`.
    """
    def __init__(self, sets, category, **kwargs):
        r"""

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: type(GrowthGroup('x^ZZ * y^ZZ'))  # indirect doctest
            <class 'sage.rings.asymptotic.growth_group_cartesian.MultivariateProduct_with_category'>
        """
        super().__init__(sets, category, order='product', **kwargs)

    CartesianProduct = CartesianProductGrowthGroups
