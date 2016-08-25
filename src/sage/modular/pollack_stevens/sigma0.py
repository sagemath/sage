# -*- coding: utf-8 -*-
r"""
The matrix monoid `\Sigma_0(N)`.

This stands for a monoid of matrices over `\ZZ`, `\QQ`, `\ZZ_p`, or `\QQ_p`,
depending on an integer `N \ge 1`. This class exists in order to act on p-adic
distribution spaces.

Over `\QQ` or `\ZZ`, it is the monoid of matrices `2\times2` matrices
`\begin{pmatrix} a & b \\ c & d \end{pmatrix}`
such that
- `ad - bc \ne 0`,
- `a` is integral and invertible at the primes dividing `N`,
- `c` has valuation at least `v_p(N)` for each `p` dividing `N` (but may be
  non-integral at other places).

The value `N=1` is allowed, in which case the second and third conditions are vacuous.

EXAMPLES::

    sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
    sage: S1 = Sigma0(1); S3 = Sigma0(3)
    sage: S1([3, 0, 0, 1])
    [3 0]
    [0 1]
    sage: S3([3, 0, 0, 1]) # boom
    Traceback (most recent call last):
    ...
    TypeError: 3 is not a unit at 3
    sage: S3([5,0,0,1])
    [5 0]
    [0 1]
    sage: S3([1, 0, 0, 3])
    [1 0]
    [0 3]
    sage: matrix(ZZ, 2, [1,0,0,1]) in S1
    True

AUTHORS:

    - David Pollack (2012): initial version
"""

# Warning to developers: when working with Sigma0 elements it is generally a
# good idea to avoid using the entries of x.matrix() directly; rather, use the
# "adjuster" mechanism. The purpose of this is to allow us to seamlessly change
# conventions for matrix actions (since there are several in use in the
# literature and no natural "best" choice).
from __future__ import print_function
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.abstract_method import abstract_method
from sage.structure.factory import UniqueFactory
from sage.structure.element import MonoidElement
from sage.categories.monoids import Monoids
from sage.categories.morphism import Morphism
from sage.structure.parent import Parent
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.unique_representation import UniqueRepresentation


class Sigma0ActionAdjuster(UniqueRepresentation):

    # Can one make an abstract class in Sage?

    @abstract_method
    def __call__(self, x):
        r"""
        Given a :class:`Sigma0element` ``x``, return four integers.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import _default_adjuster
            sage: A = _default_adjuster()
            sage: A(matrix(ZZ, 2, [1,2,3,4])) # indirect doctest
            (1, 2, 3, 4)
        """# mm TODO
        pass


class _default_adjuster(Sigma0ActionAdjuster):
    """
    A callable object that does nothing to a matrix, returning its entries in the natural order.

    INPUT:

    - ``g`` -- a `2 \times 2` matrix

    OUTPUT:

    - a 4-tuple consisting of the entries of the matrix

    EXAMPLES::

        sage: A = sage.modular.pollack_stevens.sigma0._default_adjuster(); A
        <sage.modular.pollack_stevens.sigma0._default_adjuster object at 0x...>
        sage: TestSuite(A).run()
    """
    def __call__(self, g):
        """
        EXAMPLES::

            sage: T = sage.modular.pollack_stevens.sigma0._default_adjuster()
            sage: T(matrix(ZZ,2,[1..4])) # indirect doctest
            (1, 2, 3, 4)
        """
        return tuple(g.list())


class Sigma0_factory(UniqueFactory):
    r"""
    Create the monoid of non-singular matrices, upper triangular mod `N`.

    INPUT:

    - ``N`` (integer) -- the level (should be strictly positive)
    - ``base_ring`` (commutative ring, default `\ZZ`) -- the base
      ring (normally `\ZZ` or a `p`-adic ring)
    - ``adjuster`` -- None, or a callable which takes a `2 \times 2` matrix and returns
      a 4-tuple of integers. This is supplied in order to support differing
      conventions for the action of `2 \times 2` matrices on distributions.

    EXAMPLE::

        sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
        sage: Sigma0(3)
        Monoid Sigma0(3) with coefficients in Integer Ring
    """

    def create_key(self, N, base_ring=ZZ, adjuster=None):
        r"""
        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: Sigma0.create_key(3)
            (3, Integer Ring, <sage.modular.pollack_stevens.sigma0._default_adjuster object at 0x...>)
            sage: TestSuite(Sigma0).run()
        """
        N = ZZ(N)
        if N <= 0:
            raise ValueError("Modulus should be > 0")
        if adjuster is None:
            adjuster = _default_adjuster()

        if base_ring not in (QQ, ZZ):
            try:
                if not N.is_power_of(base_ring.prime()):
                    raise ValueError("Modulus must equal base ring prime")
            except AttributeError:
                raise ValueError("Base ring must be QQ, ZZ or a p-adic field")
        return (N, base_ring, adjuster)

    def create_object(self, version, key):
        r"""
        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: Sigma0(3) # indirect doctest
            Monoid Sigma0(3) with coefficients in Integer Ring
        """
        return Sigma0_class(*key)

Sigma0 = Sigma0_factory('sage.modular.pollack_stevens.sigma0.Sigma0')


class Sigma0Element(MonoidElement):
    r"""
    An element of the monoid Sigma0. This is a wrapper around a `2 \times 2` matrix.
    """##mm TODO
    def __init__(self, parent, mat):
        r"""
        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: s = Sigma0(3)([1,4,3,3]) # indirect doctest
            sage: TestSuite(s).run()
        """
        self._mat = mat
        MonoidElement.__init__(self, parent)

    def __hash__(self):
        r"""
        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: s = Sigma0(3)([1,4,3,3])
            sage: hash(s) # indirect doctest
            11

        # TODO: the doctest is probably wrong on 32-bit machines
        """
        return hash(self.matrix())

    def det(self):
        r"""
        Return the determinant of this matrix, which is (by assumption) non-zero.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: s = Sigma0(3)([1,4,3,3])
            sage: s.det()
            -9
        """
        return self.matrix().det()

    def _mul_(self, other):
        r"""
        Return the product of two Sigma0 elements.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: s = Sigma0(3)([1,4,3,3])
            sage: t = Sigma0(15)([4,0,0,1])
            sage: u = s*t; u # indirect doctest
            [ 4  4]
            [12  3]
            sage: type(u)
            <class 'sage.modular.pollack_stevens.sigma0.Sigma0_class_with_category.element_class'>
            sage: u.parent()
            Monoid Sigma0(3) with coefficients in Integer Ring
        """
        return self.parent()(self._mat * other._mat, check=False)

    def __cmp__(self, other):
        r"""
        Compare two elements (of a common Sigma0 object).

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: s = Sigma0(3)([1,4,3,3])
            sage: t = Sigma0(3)([4,0,0,1])
            sage: s == t
            False
            sage: s == Sigma0(1)([1,4,3,3])
            True

        This uses the coercion model to find a common parent, with occasionally surprising results:

            sage: t == Sigma0(5)([4, 0, 0, 1])
            False
        """
        return cmp(self._mat, other._mat)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: s = Sigma0(3)([1,4,3,3])
            sage: s._repr_()
            '[1 4]\n[3 3]'
        """
        return self.matrix().__repr__()

    def matrix(self):
        r"""
        Return self as a matrix (forgetting the additional data that it is in Sigma0(N)).

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: s = Sigma0(3)([1,4,3,3])
            sage: sm = s.matrix()
            sage: type(s)
            <class 'sage.modular.pollack_stevens.sigma0.Sigma0_class_with_category.element_class'>
            sage: type(sm)
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: s == sm
            True
        """
        return self._mat

    def inverse(self):
        r"""
        Return the inverse of self. This will raise an error if the result is not in the monoid.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: s = Sigma0(3)([1,4,3,13])
            sage: s.inverse()
            [13 -4]
            [-3  1]
            sage: Sigma0(3)([1, 0, 0, 3]).inverse()
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        .. todo::

            In an ideal world this would silently extend scalars to `\QQ` if
            the inverse has non-integer entries but is still in `\Sigma_0(N)`
            locally at `N`. But we do not use such functionality, anyway.
        """
        return self.parent()(~self._mat)


class _Sigma0Embedding(Morphism):
    r"""
    A Morphism object giving the natural inclusion of `\Sigma_0` into the
    appropriate matrix space. This snippet of code is fed to the coercion
    framework so that "x * y" will work if ``x`` is a matrix and ``y`` is a `\Sigma_0`
    element (returning a matrix, *not* a Sigma0 element).
    """
    def __init__(self, domain):
        r"""
        TESTS::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0, _Sigma0Embedding
            sage: x = _Sigma0Embedding(Sigma0(3))
            sage: TestSuite(x).run(skip=['_test_category'])

        # TODO: The category test breaks because _Sigma0Embedding is not an instance of
        # the element class of its parent (a homset in the category of
        # monoids). I have no idea how to fix this.
        """
        Morphism.__init__(self, domain.Hom(domain._matrix_space, category=Monoids()))

    def _call_(self, x):
        r"""
        Return a matrix.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0, _Sigma0Embedding
            sage: S = Sigma0(3)
            sage: x = _Sigma0Embedding(S)
            sage: x(S([1,0,0,3])).parent() # indirect doctest
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        return x.matrix()

    def _cmp_(self, other):
        r"""
        Required for pickling.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0, _Sigma0Embedding
            sage: S = Sigma0(3)
            sage: x = _Sigma0Embedding(S)
            sage: x == loads(dumps(x))
            True
        """
        return cmp(type(self), type(other)) or cmp(self.domain(), other.domain())


class Sigma0_class(Parent):
    r"""

    """ #mm TODO
    Element = Sigma0Element

    def __init__(self, N, base_ring, adjuster):
        r"""
        Standard init function. For args documentation see the factory
        function.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: S = Sigma0(3) # indirect doctest
            sage: TestSuite(S).run()
        """
        self._N = N
        self._primes = list(N.factor())
        self._base_ring = base_ring
        self._adjuster = adjuster
        self._matrix_space = MatrixSpace(base_ring, 2)
        Parent.__init__(self, category=Monoids())
        self.register_embedding(_Sigma0Embedding(self))

    def _an_element_(self):
        r"""
        Return an element of self. This is implemented in a rather dumb way.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: S = Sigma0(3)
            sage: S.an_element() # indirect doctest
            [1 0]
            [0 1]
        """
        return self([1, 0, 0, 1])

    def level(self):
        r"""
        If this monoid is `\Sigma_0(N)`, return `N`.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: S = Sigma0(3)
            sage: S.level()
            3
        """
        return self._N

    def base_ring(self):
        r"""
        Return the base ring.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: S = Sigma0(3)
            sage: S.base_ring()
            Integer Ring
        """
        return self._base_ring

    def _coerce_map_from_(self, other):
        r"""
        Find out wheter other coerces into self.

        The *only* thing that coerces canonically into `\Sigma_0` is another
        `\Sigma_0`. It is *very bad* if integers are allowed to coerce in, as
        this leads to a noncommutative coercion diagram whenever we let
        `\Sigma_0` act on anything..

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: Sigma0(1, QQ).has_coerce_map_from(Sigma0(3, ZZ)) # indirect doctest
            True
            sage: Sigma0(1, ZZ).has_coerce_map_from(ZZ)
            False

        (If something changes that causes the last doctest above to return
        True, then the entire purpose of this class is violated, and all sorts
        of nasty things will go wrong with scalar multiplication of
        distributions. Do not let this happen!)
        """
        return (isinstance(other, Sigma0_class)
                and self.level().divides(other.level())
                and self.base_ring().has_coerce_map_from(other.base_ring()))

    def _element_constructor_(self, x, check=True):
        r"""
        Construct an element of self from x.

        INPUT:

        - ``x`` -- something that one can make into a matrix over the
          appropriate base ring
        - ``check`` (boolean, default True) -- if True, then check that this
          matrix actually satisfies the conditions.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: S = Sigma0(3)
            sage: S([1,0,0,3]) # indirect doctest
            [1 0]
            [0 3]
            sage: S([3,0,0,1]) # boom
            Traceback (most recent call last):
            ...
            TypeError: 3 is not a unit at 3
            sage: S(Sigma0(1)([3,0,0,1]), check=False) # don't do this
            [3 0]
            [0 1]
        """
        if isinstance(x, Sigma0Element):
            x = x.matrix()
        if check:
            x = self._matrix_space(x)
            a, b, c, d = self._adjuster(x)
            for (p, e) in self._primes:
                if c.valuation(p) < e:
                    raise TypeError("level %s^%s does not divide %s" % (p, e, c))
                if a.valuation(p) != 0:
                    raise TypeError("%s is not a unit at %s" % (a, p))
            if x.det() == 0:
                raise TypeError("matrix must be nonsingular")
        x.set_immutable()
        return self.element_class(self, x)

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: S = Sigma0(3)
            sage: S._repr_()
            'Monoid Sigma0(3) with coefficients in Integer Ring'
        """
        return 'Monoid Sigma0(%s) with coefficients in %s' % (self.level(),
                                                              self.base_ring())
