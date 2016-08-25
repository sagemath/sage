r"""
Finite simple continued fractions

Sage implements the field ``ContinuedFractionField`` (or ``CFF``
for short) of finite simple continued fractions. This is really
isomorphic to the field `\QQ` of rational numbers, but with different
printing and semantics.  It should be possible to use this field in
most cases where one could use `\QQ`, except arithmetic is *much* slower.

EXAMPLES:

We can create matrices, polynomials, vectors, etc., over the continued fraction
field::

    sage: a = random_matrix(CFF, 4)
    doctest:...: DeprecationWarning: CFF (ContinuedFractionField) is deprecated, use QQ instead
    See http://trac.sagemath.org/20012 for details.
    sage: a
    [    [-1; 2] [-1; 1, 94]      [0; 2]       [-12]]
    [       [-1]      [0; 2]  [-1; 1, 3]   [0; 1, 2]]
    [    [-3; 2]         [0]   [0; 1, 2]        [-1]]
    [        [1]        [-1]      [0; 3]         [1]]
    sage: f = a.charpoly()
    doctest:...: DeprecationWarning: CFF (ContinuedFractionField) is deprecated, use QQ instead
    See http://trac.sagemath.org/20012 for details.
    sage: f
    [1]*x^4 + ([-2; 3])*x^3 + [14; 1, 1, 1, 9, 1, 8]*x^2 + ([-13; 4, 1, 2, 1, 1, 1, 1, 1, 2, 2])*x + [-6; 1, 5, 9, 1, 5]
    sage: f(a)
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    sage: vector(CFF, [1/2, 2/3, 3/4, 4/5])
    ([0; 2], [0; 1, 2], [0; 1, 3], [0; 1, 4])

AUTHORS:

- Niles Johnson (2010-08): ``random_element()`` should pass on ``*args`` and
  ``**kwds`` (:trac:`3893`).
"""
from __future__ import print_function
from __future__ import absolute_import

from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import Field
from sage.structure.element import FieldElement

from .continued_fraction import ContinuedFraction_periodic, ZZ_0
from sage.misc.superseded import deprecation

class ContinuedFractionField(UniqueRepresentation,Field):
    """
    The field of rational implemented as continued fraction.

    The code here is deprecated since in all situations it is better to use
    ``QQ``.

    .. SEEALSO::

        :func:`continued_fraction`

    EXAMPLES::

        sage: CFF
        QQ as continued fractions
        sage: CFF([0,1,3,2])
        [0; 1, 3, 2]
        sage: CFF(133/25)
        [5; 3, 8]

        sage: CFF.category()
        Category of fields

    The continued fraction field inherits from the base class
    :class:`sage.rings.ring.Field`. However it was initialised as such only
    since trac ticket :trac:`11900`::

        sage: CFF.category()
        Category of fields
    """
    class Element(ContinuedFraction_periodic,FieldElement):
        r"""
        A continued fraction of a rational number.

        EXAMPLES::

            sage: CFF(1/3)
            [0; 3]
            sage: CFF([1,2,3])
            [1; 2, 3]
        """
        def __init__(self, x1, x2=None):
            r"""
            INPUT:

            - ``parent`` - the parent

            - ``x`` - the quotients of the continued fraction

            TESTS::

                sage: TestSuite(CFF.an_element()).run()
            """
            deprecation(20012, "CFF (ContinuedFractionField) is deprecated, use QQ instead")
            ContinuedFraction_periodic.__init__(self, x1)
            FieldElement.__init__(self, parent=CFF)

        def _add_(self, other):
            r"""
            Add two continued fractions.

            EXAMPLES::

                sage: CFF(1/3) + CFF([0,1,2,3])
                [1; 30]
            """
            return self.parent()(self.value() + other.value())

        def _mul_(self, other):
            r"""
            Multiply two continued fractions.

            EXAMPLES::

                sage: CFF(1/3) * CFF([0,1,2,3])
                [0; 4, 3, 2]
            """
            return self.parent()(self.value() * other.value())

        def _div_(self, other):
            r"""
            Divides two continued fractions.

            EXAMPLES::

                sage: CFF(1/3) / CFF(4/5)
                [0; 2, 2, 2]
            """
            return self.parent()(self.value() / other.value())

        def __reduce__(self):
            r"""
            Pickling helper.

            EXAMPLES::

                sage: x = CFF(1/3)
                sage: loads(dumps(x)) == x
                True
            """
            return (self.parent(), (self.value(),))

    def __init__(self):
        r"""
        TESTS::

            sage: TestSuite(CFF(1/3)).run()
            sage: TestSuite(CFF([1,2,3])).run()
        """
        Field.__init__(self, self)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(CFF)
            \Bold{CFF}
        """
        return "\\Bold{CFF}"

    def _repr_(self):
        """
        EXAMPLES::

            sage: CFF
            QQ as continued fractions
        """
        return "QQ as continued fractions"

    def an_element(self):
        r"""
        Returns a continued fraction.

        EXAMPLES::

            sage: CFF.an_element()
            [-1; 2, 3]
        """
        return self([-1,2,3])

    def some_elements(self):
        r"""
        Return some continued fractions.

        EXAMPLES::

            sage: CFF.some_elements()
            ([0], [1], [1], [-1; 2], [3; 1, 2, 3])
        """
        return (self([0]), self([1]), self([0,1]), self([-1,2]), self([3,1,2,3]))

    def is_field(self, proof=True):
        """
        Return True.

        EXAMPLES::

            sage: CFF.is_field()
            True
        """
        return True

    def is_exact(self):
        r"""
        Return True.

        EXAMPLES::

            sage: CFF.is_exact()
            True
        """
        return True

    def is_finite(self):
        """
        Return False, since the continued fraction field is not finite.

        EXAMPLES::

            sage: CFF.is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Return 0, since the continued fraction field has characteristic 0.

        EXAMPLES::

            sage: c = CFF.characteristic(); c
            0
            sage: parent(c)
            Integer Ring
        """
        return ZZ_0

    def order(self):
        """
        EXAMPLES::

            sage: CFF.order()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    def random_element(self, *args, **kwds):
        """
        Return a somewhat random continued fraction (the result is either
        finite or ultimately periodic).

        INPUT:

        - ``args``, ``kwds`` - arguments passed to ``QQ.random_element``

        EXAMPLES::

            sage: CFF.random_element() # random
            [0; 4, 7]
        """
        from sage.rings.rational_field import QQ
        return self(QQ.random_element())

    def _element_constructor_(self, data, *extra_args):
        r"""
        Build an element of that field.

        TESTS::

            sage: CFF(1/3)
            [0; 3]
            sage: CFF([1,3,2])
            [1; 3, 2]
            sage: CFF(CFF(1/3))
            [0; 3]
        """
        if isinstance(data, FieldElement) and data.parent() is self:
            data = list(data)
        if extra_args:
            print("data", data, type(data))
            print("extra_args", extra_args, type(extra_args[0]))
            data = list(extra_args[0])
        if not isinstance(data, (tuple,list)):
            from sage.rings.rational_field import QQ
            data = QQ(data).continued_fraction_list()
        else:
            from .continued_fraction import check_and_reduce_pair
            data,_ = check_and_reduce_pair(data, [])
        return self.element_class(data)

    def _coerce_map_from_(self, R):
        r"""
        Return True for ZZ and QQ.

        EXAMPLES::

            sage: CFF.has_coerce_map_from(ZZ) # indirect doctest
            True
            sage: CFF.has_coerce_map_from(QQ)
            True
            sage: CFF.has_coerce_map_from(RR)
            False
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        return R is ZZ or R is QQ

CFF = ContinuedFractionField()

# Unpickling support is needed as the class ContinuedFractionField_class has
# been renamed into ContinuedFractionField in the ticket 14567
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.rings.contfrac', 'ContinuedFractionField_class',ContinuedFractionField)

