r"""
Multiplicative Abelian Groups With Values

Often, one ends up with a set that forms an Abelian group. It would be
nice if one could return an Abelian group class to encapsulate the
data. However,
:func:`~sage.groups.abelian_gps.abelian_group.AbelianGroup` is an
abstract Abelian group defined by generators and relations. This
module implements :class:`AbelianGroupWithValues` that allows the
group elements to be decorated with values.

An example where this module is used is the unit group of a number
field, see :mod:`sage.rings.number_field.unit_group`. The units form a
finitely generated Abelian group. We can think of the elements either
as abstract Abelian group elements or as particular numbers in the
number field. The :func:`AbelianGroupWithValues` keeps track of these
associated values.

.. warning::

    Really, this requires a group homomorphism from the abstract
    Abelian group to the set of values. This is only checked if you
    pass the ``check=True`` option to :func:`AbelianGroupWithValues`.

EXAMPLES:

Here is `\ZZ_6` with value `-1` assigned to the generator::

    sage: Z6 = AbelianGroupWithValues([-1], [6], names='g')
    sage: g = Z6.gen(0)
    sage: g.value()
    -1
    sage: g*g
    g^2
    sage: (g*g).value()
    1
    sage: for i in range(7):
    ....:     print((i, g^i, (g^i).value()))
    (0, 1, 1)
    (1, g, -1)
    (2, g^2, 1)
    (3, g^3, -1)
    (4, g^4, 1)
    (5, g^5, -1)
    (6, 1, 1)

The elements come with a coercion embedding into the
:meth:`~AbelianGroupWithValues_class.values_group`, so you can use the
group elements instead of the values::

    sage: CF3.<zeta> = CyclotomicField(3)
    sage: Z3.<g> = AbelianGroupWithValues([zeta], [3])
    sage: Z3.values_group()
    Cyclotomic Field of order 3 and degree 2
    sage: g.value()
    zeta
    sage: CF3(g)
    zeta
    sage: g + zeta
    2*zeta
    sage: zeta + g
    2*zeta
"""

##########################################################################
#  Copyright (C) 2012 Volker Braun  <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL):
#
#                  http://www.gnu.org/licenses/
##########################################################################

from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.categories.morphism import Morphism
from sage.groups.abelian_gps.abelian_group import AbelianGroup_class, _normalize
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement


def AbelianGroupWithValues(values, n, gens_orders=None, names='f', check=False, values_group=None):
    r"""
    Construct an Abelian group with values associated to the generators.

    INPUT:

    - ``values`` -- a list/tuple/iterable of values that you want to
      associate to the generators.

    - ``n`` -- integer (optional). If not specified, will be derived
       from ``gens_orders``.

    - ``gens_orders`` -- a list of non-negative integers in the form
       `[a_0, a_1, \dots, a_{n-1}]`, typically written in increasing
       order. This list is padded with zeros if it has length less
       than n. The orders of the commuting generators, with `0`
       denoting an infinite cyclic factor.

    -  ``names`` -- (optional) names of generators

    - ``values_group`` -- a parent or ``None`` (default). The common
      parent of the values. This might be a group, but can also just
      contain the values. For example, if the values are units in a
      ring then the ``values_group`` would be the whole ring. If
      ``None`` it will be derived from the values.

    EXAMPLES::

        sage: G = AbelianGroupWithValues([-1], [6])
        sage: g = G.gen(0)
        sage: for i in range(7):
        ....:     print((i, g^i, (g^i).value()))
        (0, 1, 1)
        (1, f, -1)
        (2, f^2, 1)
        (3, f^3, -1)
        (4, f^4, 1)
        (5, f^5, -1)
        (6, 1, 1)
        sage: G.values_group()
        Integer Ring

    The group elements come with a coercion embedding into the
    :meth:`values_group`, so you can use them like their
    :meth:`~sage.groups.abelian_gps.value.AbelianGroupWithValuesElement.value`
    ::

        sage: G.values_embedding()
        Generic morphism:
          From: Multiplicative Abelian group isomorphic to C6
          To:   Integer Ring
        sage: g.value()
        -1
        sage: 0 + g
        -1
        sage: 1 + 2*g
        -1
    """
    if check:
        raise NotImplementedError('checking that the values are a homomorphism is not implemented')
    gens_orders, names = _normalize(n, gens_orders, names)
    if values_group is None:
        from sage.structure.sequence import Sequence
        values_group = Sequence(values).universe()
    values = tuple( values_group(val) for val in values )
    M = AbelianGroupWithValues_class(gens_orders, names, values, values_group)
    return M


class AbelianGroupWithValuesEmbedding(Morphism):
    """
    The morphism embedding the Abelian group with values in its values group.

    INPUT:

    - ``domain`` -- a :class:`AbelianGroupWithValues_class`

    - ``codomain`` -- the values group (need not be in the category of
      groups, e.g. symbolic ring).

    EXAMPLES::

        sage: Z4.<g> = AbelianGroupWithValues([I], [4])
        sage: embedding = Z4.values_embedding();  embedding
        Generic morphism:
          From: Multiplicative Abelian group isomorphic to C4
          To:   Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        sage: embedding(1)
        1
        sage: embedding(g)
        I
        sage: embedding(g^2)
        -1
    """

    def __init__(self, domain, codomain):
        """
        Construct the morphism

        TESTS::

            sage: Z4 = AbelianGroupWithValues([I], [4])
            sage: from sage.groups.abelian_gps.values import AbelianGroupWithValuesEmbedding
            sage: AbelianGroupWithValuesEmbedding(Z4, Z4.values_group())
            Generic morphism:
              From: Multiplicative Abelian group isomorphic to C4
              To:   Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        assert domain.values_group() is codomain
        from sage.categories.homset import Hom
        Morphism.__init__(self, Hom(domain, codomain))

    def _call_(self, x):
        """
        Return the value associated to ``x``

        INPUT:

        - ``x`` -- a group element

        OUTPUT:

        Its value.

        EXAMPLES::

            sage: Z4.<g> = AbelianGroupWithValues([I], [4])
            sage: embedding = Z4.values_embedding()
            sage: embedding(g)
            I
            sage: embedding._call_(g)
            I
        """
        return x.value()


class AbelianGroupWithValuesElement(AbelianGroupElement):
    """
    An element of an Abelian group with values assigned to generators.

    INPUT:

    - ``exponents`` -- tuple of integers. The exponent vector defining
      the group element.

    - ``parent`` -- the parent.

    - ``value`` -- the value assigned to the group element or ``None``
      (default). In the latter case, the value is computed as needed.

    EXAMPLES::

        sage: F = AbelianGroupWithValues([1,-1], [2,4])
        sage: a,b = F.gens()
        sage: TestSuite(a*b).run()
    """

    def __init__(self, parent, exponents, value=None):
        """
        Create an element

        EXAMPLES::

            sage: F = AbelianGroupWithValues([1,-1], [2,4])
            sage: a,b = F.gens()
            sage: a*b^-1 in F
            True
            sage: (a*b^-1).value()
            -1
        """
        self._value = value
        AbelianGroupElement.__init__(self, parent, exponents)

    def value(self):
        """
        Return the value of the group element.

        OUTPUT:

        The value according to the values for generators, see
        :meth:`~AbelianGroupWithValues.gens_values`.

        EXAMPLES::

            sage: G = AbelianGroupWithValues([5], 1)
            sage: G.0.value()
            5
        """
        if self._value is None:
            values = self.parent().gens_values()
            self._value = prod( v**e for v,e in zip(values, self.exponents()) )
        return self._value

    def _div_(left, right):
        """
        Divide ``left`` by ``right``

        TESTS::

            sage: G.<a,b> = AbelianGroupWithValues([5,2], 2)
            sage: a._div_(b)
            a*b^-1
            sage: a/b
            a*b^-1
            sage: (a/b).value()
            5/2
        """
        m = AbelianGroupElement._div_(left, right)
        m._value = left.value() / right.value()
        return m

    def _mul_(left, right):
        """
        Multiply ``left`` and ``right``

        TESTS::

            sage: G.<a,b> = AbelianGroupWithValues([5,2], 2)
            sage: a._mul_(b)
            a*b
            sage: a*b
            a*b
            sage: (a*b).value()
            10
        """
        m = AbelianGroupElement._mul_(left, right)
        m._value = left.value() * right.value()
        return m

    def __pow__(self, n):
        """
        Exponentiate ``self``

        INPUT:

        - ``n`` -- integer. The exponent.

        TESTS::

            sage: G.<a,b> = AbelianGroupWithValues([5,2], 2)
            sage: a^3
            a^3
            sage: (a^3).value()
            125
        """
        m = Integer(n)
        if n != m:
            raise TypeError('argument n (= '+str(n)+') must be an integer.')
        pow_self = AbelianGroupElement.__pow__(self, m)
        pow_self._value = pow(self.value(), m)
        return pow_self

    def inverse(self):
        """
        Return the inverse element.

        EXAMPLES::

            sage: G.<a,b> = AbelianGroupWithValues([2,-1], [0,4])
            sage: a.inverse()
            a^-1
            sage: a.inverse().value()
            1/2
            sage: a.__invert__().value()
            1/2
            sage: (~a).value()
            1/2
            sage: (a*b).value()
            -2
            sage: (a*b).inverse().value()
            -1/2
        """
        m = AbelianGroupElement.inverse(self)
        m._value = ~self.value()
        return m

    __invert__ = inverse



class AbelianGroupWithValues_class(AbelianGroup_class):
    """
    The class of an Abelian group with values associated to the generator.

    INPUT:

    - ``generator_orders`` -- tuple of integers. The orders of the
      generators.

    - ``names`` -- string or list of strings. The names for the generators.

    - ``values`` -- Tuple the same length as the number of
      generators. The values assigned to the generators.

    - ``values_group`` -- the common parent of the values.

    EXAMPLES::

        sage: G.<a,b> = AbelianGroupWithValues([2,-1], [0,4])
        sage: TestSuite(G).run()
    """
    Element = AbelianGroupWithValuesElement

    def __init__(self, generator_orders, names, values, values_group):
        """
        The Python constructor

        TESTS::

            sage: G = AbelianGroupWithValues([2,-1], [0,4]); G
            Multiplicative Abelian group isomorphic to Z x C4

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.explain(G, ZZ, operator.add)
            Coercion on left operand via
                Generic morphism:
                  From: Multiplicative Abelian group isomorphic to Z x C4
                  To:   Integer Ring
            Arithmetic performed after coercions.
            Result lives in Integer Ring
            Integer Ring
        """
        self._values = values
        self._values_group = values_group
        AbelianGroup_class.__init__(self, generator_orders, names)
        self._populate_coercion_lists_(embedding=self.values_embedding())
        if self.ngens() != len(self._values):
            raise ValueError('need one value per generator')

    def gen(self, i=0):
        """
        The `i`-th generator of the abelian group.

        INPUT:

        - ``i`` -- integer (default: 0). The index of the generator.

        OUTPUT:

        A group element.

        EXAMPLES::

            sage: F = AbelianGroupWithValues([1,2,3,4,5], 5,[],names='a')
            sage: F.0
            a0
            sage: F.0.value()
            1
            sage: F.2
            a2
            sage: F.2.value()
            3

            sage: G = AbelianGroupWithValues([-1,0,1], [2,1,3])
            sage: G.gens()
            (f0, 1, f2)
        """
        g = AbelianGroup_class.gen(self, i)
        g._value = self._values[i]
        return g

    def gens_values(self):
        """
        Return the values associated to the generators.

        OUTPUT:

        A tuple.

        EXAMPLES::

            sage: G = AbelianGroupWithValues([-1,0,1], [2,1,3])
            sage: G.gens()
            (f0, 1, f2)
            sage: G.gens_values()
            (-1, 0, 1)
        """
        return self._values

    def values_group(self):
        """
        The common parent of the values.

        The values need to form a multiplicative group, but can be
        embedded in a larger structure. For example, if the values are
        units in a ring then the :meth:`values_group` would be the
        whole ring.

        OUTPUT:

        The common parent of the values, containing the group
        generated by all values.

        EXAMPLES::

            sage: G = AbelianGroupWithValues([-1,0,1], [2,1,3])
            sage: G.values_group()
            Integer Ring

            sage: Z4 = AbelianGroupWithValues([I], [4])
            sage: Z4.values_group()
            Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        return self._values_group

    def values_embedding(self):
        """
        Return the embedding of ``self`` in :meth:`values_group`.

        OUTPUT:

        A morphism.

        EXAMPLES::

            sage: Z4 = AbelianGroupWithValues([I], [4])
            sage: Z4.values_embedding()
            Generic morphism:
              From: Multiplicative Abelian group isomorphic to C4
              To:   Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        return AbelianGroupWithValuesEmbedding(self, self.values_group())
