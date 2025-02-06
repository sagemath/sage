# -*- coding: utf-8 -*-
r"""
The Chow Ring of a Scheme

In general, the Chow ring of an algebraic variety is a graded ring, an algebraic
analogue to the cohomology ring of its underlying topological space. Its
elements correspond to formal linear combinations of algebraic subvarieties
modulo rational equivalence. The product structure comes from the intersection
of subvarieties.

See http://en.wikipedia.org/wiki/Chow_ring for an introduction to Chow Rings.

Calculating the Chow ring of an algebraic (projective, smooth)
variety is a non trivial task. Here we suppose that the Chow ring is already
given by generators and relations and will first provide some standard methods
as computing a linear basis, the intersection matrix or the dual basis.
Note also that we always work over `\QQ`.

Later we will provide methods to compute the Chow ring for related algebraic
varieties as a Grassmann bundle or a blowup (grass, blowup).

For example the Chow ring of the projective line is given by its generator `h`
in degree `1` subject to the relation `h^2`::

    sage: A = ChowRing('h', 1, 'h^2')

Another example is the Chow ring of the projective plane blown up in a point
given by two classes `E` and `H` in degree one subject to the relations
`E.H=0` and `E^2+H^2=0`. Note that the latter relations are not in a standard
basis. However when retrieving the relations, a standard basis is returned::

    sage: A = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
    sage: A.rels()
    [E^3, H^2 + E^2, H*E]

An important method for a ChowRing is the computation of a linear basis
of algebraic cycles :meth:`ChowRing_generic.basis`::

    sage: A = ChowRing('h', 1, 'h^3')  # P2
    sage: A.basis()
    [h^2, h, 1]

Note however we do not ask for a ChowRing to be Artin, even though this is the
case, at least for projective smooth algebraic varieties. So the following is
perfectly valid for us::

    sage: A = ChowRing(['c1', 'c2'], [1, 2])

The reason for this is that it is sometimes very useful to compute in these
more general rings, especially while computing the relations for a particular
variety. Note that if `A` is not Artin, it does not make sens to compute a
linear basis and an error message is thrown::

    sage: A.basis()
    Traceback (most recent call last):
    ...
    ValueError: Ring is not of dimension 0.

If A is the Chow ring of a projective smooth algebraic variety, a *point_class*
is an element in top degree that integrates to one. For example,
for the projective space of dimension `n` with generator `h`, the
point_class is `h^n`.
Once we have a point class, we can compute the intersection matrix
and the basis dual to those given by :meth:`ChowRing_generic.basis` with
respect to the intersection matrix::

    sage: A.<h> = ChowRing('h', 1, 'h^4')  # P3
    sage: A.basis()
    [h^3, h^2, h, 1]
    sage: A.set_point_class(h^3)
    sage: A.intersection_matrix()
    [0 0 0 1]
    [0 0 1 0]
    [0 1 0 0]
    [1 0 0 0]
    sage: A.dual_basis_slow()
    [1, h, h^2, h^3]

Note that :meth:`dual_basis_slow` is, as its name indicates, slow at least for
large examples (e.g. with a large linear basis of A). This is due to the fact
that :meth:`ChowRing_generic.dual_basis_slow` simply computes the intersection
matrix, then takes its inverse. A much faster method is
:meth:`ChowRing_generic.dual_basis` which computes only the relevant (i.e. in
complementary degrees) part of the intersection matrix.
Setting the optional parameter verbose=True allows to follow the computation::

    sage: A.<H, E> = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
    sage: A.basis()
    [E^2, E, H, 1]
    sage: A.set_point_class(H^2)
    sage: A.dual_basis(verbose=True)
    Computing block 0 of size 1 ...
    Computing block 1 of size 2 ...
    [-1, -E, H, -E^2]

Finally, there is a particular PointChowRing which is defined as the ChowRing
with no generators and relations::

    sage: A = ChowRing()
    sage: A.set_dimension(0)
    sage: A.set_point_class(1)
    sage: B = PointChowRing
    sage: A == B
    True

In this case the underlying ring is `\QQ[]/(O)` and not simply `\QQ` as one
might suspect, in order to get uniformly quotient rings::

    sage: A = ChowRing(); A
    Quotient of Multivariate Polynomial Ring in no variables over Rational Field by the ideal (0)
    sage: A.gens()
    ()
    sage: A.ngens()
    0
    sage: A.gen()
    Traceback (most recent call last):
    ...
    ValueError: Generator not defined.
    sage: A.variable_names()
    ()
    sage: A.variable_name()
    Traceback (most recent call last):
    ...
    IndexError: tuple index out of range

In the general case, when get as expected::

    sage: A = ChowRing('h', 1, 'h^2')  # ChowRing of P1
    sage: A.gens()
    (h,)
    sage: A.variable_names()
    ('h',)
    sage: A.variable_name()
    'h'
    sage: A = ChowRing(['c1', 'c2'], [1, 2])
    sage: A.variable_names()
    ('c1', 'c2')
    sage: A.variable_name()
    'c1'

TESTS::

    sage: A = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
    sage: TestSuite(A).run()

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""

# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import QQ
from sage.modules.all import vector
from sage.matrix.all import matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.quotient_ring import QuotientRing_generic
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.libs.singular.function_factory import singular_function, lib as singular_lib
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer


def is_chowRing(R):
    """
    Test whether ``R`` is a ChowRing.

    INPUT:

    - ``R`` -- anything.

    OUTPUT:

    Boolean. Whether ``R`` derives from :class:`ChowRing_generic`.

    TESTS::

        sage: A = PointChowRing
        sage: is_chowRing(A)
        True
        sage: S = PolynomialRing(QQ, ['x', 'y']).quotient('x*y')
        sage: is_chowRing(S)
        False
    """
    return isinstance(R, ChowRing_generic)


def ChowRingFromRing(S, names=None, name=None, latex_name=None):
    """
    Returns a ChowRing, given a ring S.

    INPUT:

    - ``S`` -- A (quotient of) a multivariate Polynomial Ring

    OUTPUT:

    - The ChowRing defined by ``S``

    EXAMPLES::

        sage: from sage.schemes.chow.ring import ChowRingFromRing
        sage: variables, degrees, relations = ['x', 'y'], [2, 3], ['x*y']
        sage: AX = ChowRing(variables, degrees, relations)
        sage: OT = TermOrder('wdegrevlex', (2, 3))
        sage: S = PolynomialRing(QQ, ['x', 'y'], order=OT).quotient('x*y'); S
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y)
        sage: AS = ChowRingFromRing(S)
        sage: AS == AX
        True

    Note that we expect the ring ``S`` to be (quotient of) a multivariate
    Polynomial Ring::

        sage: S = PolynomialRing(QQ, 'h').quotient('h^2'); S
        Univariate Quotient Polynomial Ring in hbar over Rational Field with modulus h^2
        sage: AS = ChowRingFromRing(S)
        Traceback (most recent call last):
        ...
        TypeError: Expected (Quotient of) Multivariate Polynomial Ring.

    """

    # Ring Check
    if isinstance(S, MPolynomialRing_base):
        R, I = S, S.ideal(0)
    elif isinstance(S, QuotientRing_generic):
        R, I = S.cover_ring(), S.defining_ideal()
    else:
        err = "Expected (Quotient of) Multivariate Polynomial Ring."
        raise TypeError(err)
    if not(isinstance(R, MPolynomialRing_base)):
        err = "Expected (Quotient of) Multivariate Polynomial Ring."
        raise TypeError(err)
    # Get names and check that numbers correspond
    names = R.variable_names() if names is None else names
    if len(names) != R.ngens():
        raise ValueError("Number of names different from number of variables.")
    return ChowRing_generic(R, I, names, name, latex_name)


def ChowRing(generators=None, degrees=None, relations=None,
             names=None, name=None, latex_name=None):
    r"""
    Return a ChowRing, given its dimension, generators, degrees, relations
    and point_class.

    INPUT:

    - ``generators``-- An optional (list of) strings, the generators

    - ``degrees``-- An optional (list of) integers, the degrees of the generators

    - ``relations``-- An optional (list of) strings, the relations between the generators

    - ``names``-- An optional (list of) strings, names of the generators in the ChowRing

    - ``name``-- An optional string, the name of the ChowRing

    - ``latex_name``-- An optional string, a latex representation of the ChowRing

    OUTPUT:

    - A ChowRing

    EXAMPLES:

    The ChowRing of the projective line is given by a generator `h` in degree
    `1`, subject to the relation `h^2`::

        sage: A = ChowRing('h', 1, 'h^2')

    It can equally be defined more verbosely by::

        sage: B = ChowRing(generators='h', degrees=1, relations='h^2')
        sage: A == B
        True

    If there is more than one generator (or more than one relation) they have
    to be specified as lists of strings::

        sage: A = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])

    Note that by default the generators in the ChowRing have the same names
    as in the cover_ring as opposed to a quotient by a PolynomialRing::

        sage: A = ChowRing('h', 1, 'h^2')
        sage: A.gens()
        (h,)
        sage: P = PolynomialRing(QQ, 'h').quotient('h^2')
        sage: P.gens()
        (hbar,)

    As for QuotientRings, the generators can be named explicitly::

        sage: A.<x> = ChowRing('h', 1, 'h^2')
        sage: A.gens()
        (x,)

    Only the dimension is mandatory. For example, one may define a point
    Chow Ring as follows::

        sage: A = ChowRing(0); A
        Quotient of Multivariate Polynomial Ring in no variables over Rational Field by the ideal (0)

    Please note that the internal structure of A is not the field `\QQ` but
    rather the Quotient of `\QQ` by the ideal `(0)`. This is due to the fact
    that ChowRings are uniformly represented as quotient rings for us.

    TESTS::

        sage: A = ChowRing('c1', 'c2')
        Traceback (most recent call last):
        ...
        TypeError: Expect list of integers for degrees.

        sage: A = ChowRing('h', -1)
        Traceback (most recent call last):
        ...
        ValueError: Expect positive integers for degrees.

        sage: A = ChowRing('h', 0)
        Traceback (most recent call last):
        ...
        ValueError: Except a positive degree for every generator.

        sage: A = ChowRing('h', [1,2])
        Traceback (most recent call last):
        ...
        ValueError: Number of generators and degrees differ.

        sage: A = ChowRing(degrees=1, relations='h^2')
        No generators, degrees ignored.
        No generators, relations ignored.

        sage: A = ChowRing([], [], [])
        sage: A.set_point_class(1)
        sage: A(3).integral()
        3
    """
    # Construct a quotient ring by generators and relations.
    if generators:  # Not None neither []
        if isinstance(generators, str):
            generators = [generators]
        for gen in generators:
            if not isinstance(gen, str):
                raise TypeError("Expect list of strings for generators.")
        if degrees:
            if isinstance(degrees, Integer) or isinstance(degrees, int):
                degrees = [int(degrees)]
            for deg in degrees:
                deg = int(deg) if isinstance(deg, Integer) else deg
                if not isinstance(deg, int):
                    raise TypeError("Expect list of integers for degrees.")
                if deg <= 0:
                    raise ValueError("Expect positive integers for degrees.")
            if len(generators) != len(degrees):
                raise ValueError("Number of generators and degrees differ.")
        else:
            raise ValueError("Except a positive degree for every generator.")
        R = PolynomialRing(QQ, len(generators), names=tuple(generators),
                           order=TermOrder('wdegrevlex', tuple(degrees))) # type: ignore
        I = R.ideal(0)
        if relations:
            if isinstance(relations, str):
                relations = [relations]
            if not all(isinstance(r, str) for r in relations):
                raise TypeError("Expect list of strings as relations")
            I = R.ideal([R(rel) for rel in relations]) # type: ignore
    else:
        # The ChowRing of the point: QQ[]/(0). This might seem artificial but
        # permits us not to distinguish between QQ as a rational field
        # and QQ seen of type QuotientRing_generic.
        R = PolynomialRing(QQ, 0, names=[])
        I = R.ideal(0)
        if degrees:
            print(Warning("No generators, degrees ignored."))
        if relations:
            print(Warning("No generators, relations ignored."))
    return ChowRingFromRing(R.quotient(I), names, name, latex_name)


class ChowRing_generic(QuotientRing_generic):

    from sage.schemes.chow.ring_element import ChowRingElement
    Element = ChowRingElement

    def __init__(self, R, I, names=None, name=None, latex_name=None):
        r"""
        Construct a :class:`ChowRing_generic`.

        .. WARNING::

            This class does not perform any checks of correctness of input.
            Use :func:`ChowRing` or :func:`ChowRingFromRing` to construct a
            ChowRing.

        INPUT:

        - ``R`` -- `\QQ` or a multivariate polynomial ring;

        - ``I`` -- an ideal of ``R``;

        - ``names`` -- list of variable names;

        - ``name``-- An optional string, the name of the ChowRing

        - ``latex_name``-- An optional string, a latex representation of the ChowRing

        OUTPUT:

        - :class:`ChowRing <ChowRing_generic>`.

        TESTS ::

            sage: A = ChowRing()

        """
        if R.ngens() != 0:
            # Ensure standard basis before taking the quotient.
            I = R.ideal(I.groebner_basis('libsingular:std'))
        QuotientRing_generic.__init__(self, R, I, names)
        self._name = name
        self._latex_name = latex_name
        # Options
        self._dimension = None
        self._point_class = None

    def _cache_key(self):
        return(self.parent(),str(self))

    def __hash__(self):
        return hash(type(self))

    def __eq__(self, other):
        """
        Two ChowRings are considered to be "equal" their dimensions are equal,
        the map sending the i-th generator to the i-th generator of the
        underlying rings induces and isomorphism, and the point_classes are the
        same under this identification.

        TESTS::

            sage: from sage.schemes.chow.ring import ChowRingFromRing
            sage: AX = ChowRing('h', 1, 'h^2')
            sage: S = PolynomialRing(QQ, 1, ['h']).quotient('h^2')
            sage: AS = ChowRingFromRing(S)
            sage: AS == AX
            True

        """
        # Check instances
        if not isinstance(self, ChowRing_generic):
            return False
        if not isinstance(other, ChowRing_generic):
            return False
        # Check dimension
        if self.dimension() != other.dimension():
            return False
        # Check number of generators
        if self.ngens() != other.ngens():
            return False
        # Either there are or not. In the latter cas, the underlying cover_ring
        # is the rational field.
        if self.ngens() != 0:
            # Check degrees
            if self.degs() != other.degs():
                return False
            A, I = self.cover_ring(), self.defining_ideal()
            B, J = other.cover_ring(), other.defining_ideal()
            f = A.hom(other.gens(), B)
            JJ = f.pushforward(I)
            if not ideals_are_equal(J, JJ, B):
                return False
                # Finally check that other.point_class() equals self.point_class()
            sopc = None
            if self.point_class() is not None:
                spc = A.ideal(self.point_class())
                sopc = f.pushforward(spc).gen(0)
        else:
            # Rational field. Just check point_classes
            sopc = None if self.point_class() is None else self.point_class()
            # Check PointChowScheme Classes:
        opc = other.point_class()
        if not(opc is None and sopc is None):  # At least one is not None
            if opc is None or sopc is None:    # One is None
                return False
            if str(opc) != str(sopc):          # Need string comparison here
                return False
        return True

    def __ne__(self, other):
        """
        The opposite of equal.

        TESTS::

            sage: A = ChowRing()
            sage: B = ChowRing('h', 1)
            sage: A != B
            True
        """
        return not self.__eq__(other)

    def __call__(self, x):
        r"""
        Create an element of this group ChowGroup.

        INPUT:

            -- x: anything that coerces to a ChowRingElement

        OUTPUT:

            -- a ChowRingElement instance whose parent is self.

        EXAMPLE::

            sage: A = ChowRing('h', 1, 'h^3')
            sage: A('3+2*h+2*h^2')
            2*h^2 + 2*h + 3
        """
        return self._element_constructor_(x)

    def _repr_(self):
        r"""
        Return a string representation of this ChowRing.

        OUTPUT:

        - string.

        TESTS::

            sage: PointChowRing._repr_()
            'Quotient of Multivariate Polynomial Ring in no variables over Rational Field by the ideal (0)'
        """
        return QuotientRing_generic._repr_(self)

    def __str__(self):
        return self._name if self._name else self._repr_()

    def _latex_(self):
        r"""
        Return a latex representation of this ChowRing.

        OUTPUT:

        - string.

        TESTS::

            sage: PointChowRing._latex_()
            '\\Bold{Q}[]/\\left(0\\right)\\Bold{Q}[]'
        """
        return self._name if self._name else QuotientRing_generic._latex_(self)

    def dimension(self):
        r"""
        Return the variety dimension (see `:meth:set_dimension`).

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^2')
            sage: A.dimension()  # None as not set yet.
            sage: A.set_dimension(1)
            sage: A.dimension()
            1

        """
        return self._dimension

    def set_dimension(self, dimension):
        r"""
        Set the dimension. If this ChowRing is the Chow ring of a
        variety `X`, set the dimension to `\dim(X)`

        INPUT:

        - an integer

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^3')   # P2
            sage: A.set_dimension(2)    # dim(P2)

        TESTS::

            sage: A = ChowRing('h', 1, 'h^3')   # P2
            sage: A.set_dimension(1)
            Traceback (most recent call last):
            ...
            ValueError: Expect the dimension to be at least max-degree.
        """
        # Dimension check
        dimension = int(dimension) if isinstance(dimension, Integer) else dimension
        if not isinstance(dimension, int):
            raise TypeError("Expect an integer for the dimension.")
        if dimension < 0:
            raise ValueError("Expect non negative integer for the dimension.")
        if dimension < self.max_degree():
            raise ValueError("Expect the dimension to be at least max-degree.")
        self._dimension = dimension

    def point_class(self):
        r"""
        Return the point class.

        OUTPUT:

        - a ring element

        EXAMPLES::

            sage: A.<h> = ChowRing('h', 1, 'h^3')   # P2
            sage: A.set_point_class(h^2)
            sage: A.point_class()
            h^2

        """
        return self._point_class

    def set_point_class(self, point_class):
        r"""
        Set the point class.

        EXAMPLES::

            sage: A.<h> = ChowRing('h', 1, 'h^2')   # P2
            sage: A.set_point_class(h)

        TESTS::

            sage: A = ChowRing('h', 1, 'h^3')   # P2
            sage: A.set_point_class('h^2')
            sage: A.point_class()
            h^2
            sage: A.set_point_class('k^2')
            Traceback (most recent call last):
            ...
            TypeError: Can't coerce pointclass to ring element.
            sage: A.set_point_class('h')
            Traceback (most recent call last):
            ...
            ValueError: Expect the point class in max-degree.
        """
        # Check point_class
        try:
            self(point_class)
        except Exception as e:
            raise TypeError("Can't coerce pointclass to ring element.")
        pc = self(point_class)
        # Check in max degree
        if pc.lift().degree() != self.max_degree():
            raise ValueError("Expect the point class in max-degree.")
        self._point_class = pc

    @cached_method
    def degs(self):
        """
        Return the degrees of the generators of this ChowRing.

        OUTPUT:

        A tuple of integers representing the degrees of the generators.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^2')  # ChowRing of P1
            sage: A.degs()
            (1,)
            sage: A = ChowRing(['c1', 'c2'], [1, 2])
            sage: A.degs()
            (1, 2)

        CORNER CASE::

            sage: A = PointChowRing
            sage: A.degs()
            ()
        """
        if self.ngens() == 0:
            return ()
        R = self.cover_ring()
        deg = singular_function('deg')
        return tuple([deg(x, ring=R) for x in R.gens()])

    def deg(self, i=None):
        """
        Return the degree of the i-th generator.  If i is unspecified the
        degree of the first generator (eg i=0) is returned if exists.

        INPUT:

        - ``i`` -- an (optional) integer.

        OUTPUT:

        The degree of the i-th generator.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^2')  # ChowRing of P1
            sage: A.deg()
            1
            sage: A = ChowRing(['c1', 'c2'], [1, 2])
            sage: A.deg()
            1
            sage: A.deg(1)
            2

        CORNER CASE::

            sage: A = PointChowRing
            sage: A.deg()
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        i = 0 if i is None else i  # Default value
        return self.degs()[i]

    @cached_method
    def rels(self):
        """
        Return the relations of this ChowRing.

        OUTPUT:

        A list of ring elements representing the relations.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^2')  # ChowRing of P1
            sage: A.rels()
            [h^2]
            sage: A = ChowRing(['c1', 'c2'], [1, 2])
            sage: A.rels()
            []

        TESTS::

            sage: A = PointChowRing
            sage: A.rels()
            []
        """
        relations = self.defining_ideal().gens()
        if relations == [0] or relations == (0,):
            return []
        return relations

    @cached_method
    def nrels(self):
        """
        Return the number of relations of this ChowRing

        OUTPUT:

        An integer, the number of relations.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^2')  # ChowRing of P1
            sage: A.nrels()
            1
            sage: A = ChowRing(['c1', 'c2'], [1, 2])
            sage: A.nrels()
            0
        """
        return len(self.rels())

    def rel(self, i=None):
        """
        Return the i-th relation.  If i is unspecified the first relation
        (eg i=0) is returned if exists.

        INPUT:

        - ``i`` -- an (optional) integer.

        OUTPUT:

        The i-th relation.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^2')  # ChowRing of P1
            sage: A.rel()
            h^2

        TESTS::

            sage: A = ChowRing(['c1', 'c2'], [1, 2])
            sage: A.rel()
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        i = 0 if i is None else i  # Default value
        return self.rels()[i]

    def is_point(self):
        """
        Returns True if this ChowRing corresponds to a point, e.g. is of
        dimension 0 and has no generators.

        OUTPUT:

        - True or False

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^4')  # ChowRing of P3
            sage: A.is_point()
            False

            sage: A = ChowRing()
            sage: A.is_point()
            True
        """
        return self.ngens() == 0

    def is_point_chowring(self):
        """
        Returns True if this ChowRing is the instance PointChowRing.

        OUTPUT:

        - True or False

        EXAMPLES::

            sage: A = ChowRing()
            sage: A.is_point_chowring()
            False

            sage: A = PointChowRing
            sage: A.is_point_chowring()
            True

        """
        return self is PointChowRing

    # noinspection PyUnusedLocal
    def _Hom_(self, Y, category=None, check=True):
        from sage.schemes.chow.ring_homset import ChowRingHomSet
        return ChowRingHomSet(self, Y)

    def identity_morphism(self):
        """
        Returns the identity morphism of this ChowRing.

        OUTPUT:

        - The identity morphism

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^2')
            sage: f = A.identity_morphism()
            sage: h = A.gen()
            sage: f(h)
            h
            sage: f(1)
            1

            sage: A = PointChowRing
            sage: f = A.identity_morphism()
            sage: f(1)
            1

        """
        return self.hom([self(x) for x in self.cover_ring().gens()], self)

    @cached_method
    def krull_dimension(self):
        """
        Returns the Krull dimension of the underlying ring of this ChowRing.
        Remark: this is *not* the dimension of this ChowRing specified while
        defining this ChowRing.

        OUTPUT:

        - An integer, the Krull dimension of this ring.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^4')  # P3
            sage: A.krull_dimension()
            0
            sage: A = ChowRing(['c1', 'c2'], [1, 2])  # No relations!
            sage: A.krull_dimension()
            2
        """
        if self.ngens() == 0:
            return 0
        return self.defining_ideal().dimension()

    @cached_method
    def basis(self):
        """
        Returns a basis of the underlying ring of this ChowRing if the
        underlying ring is Artin.

        OUTPUT:

        - A list of elements of the underlying ring.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^4')  # P3
            sage: A.basis()
            [h^3, h^2, h, 1]

            sage: A = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: A.basis()
            [E^2, E, H, 1]

            sage: A = ChowRing(['c1', 'c2'], [1, 2])  # No relations!
            sage: A.basis()
            Traceback (most recent call last):
            ...
            ValueError: Ring is not of dimension 0.

        TESTS::

            sage: A = PointChowRing
            sage: A.basis()
            [1]

        """
        if self.ngens() == 0:  # Special case of no generators
            return [1]
        if self.krull_dimension() != 0:
            raise ValueError("Ring is not of dimension 0.")
        return [self(b) for b in self.defining_ideal().normal_basis()]

    @cached_method
    def max_degree(self):
        r"""
        Return the maximal degree of this ChowRing if Artin else 0.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^4')  # P3
            sage: A.max_degree()
            3

            sage: A = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: A.max_degree()
            2

            sage: A = ChowRing('h', 1)  # not Artin
            sage: A.max_degree()
            0

        TESTS::

            sage: A = PointChowRing
            sage: A.max_degree()
            0
        """
        # Corner case:
        if self.ngens() == 0:
            return 0
        # Artin ?
        if self.krull_dimension() != 0:
            return 0
        # Compute max degree:
        return max([self(x).lift().degree() for x in self.basis()])

    @cached_method
    def basis_by_degree(self):
        """
        Return a basis of the underlying ring of this ChowRing ordered
        by the degrees of the elements.

        OUTPUT:

        - A list of list of elements of the underlying ring.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^4')  # P3
            sage: A.basis_by_degree()
            [[1], [h], [h^2], [h^3]]

            sage: A = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: A.basis_by_degree()
            [[1], [E, H], [E^2]]

        TESTS::

            sage: A = PointChowRing
            sage: A.basis_by_degree()
            [[1]]
        """
        if self.ngens() == 0:  # Special case of no generators
            return [[1]]
        bbd = [[] for _ in range(self.max_degree() + 1)]
        for e in self.basis():
            d = self(e).lift().degree()
            bbd[d].append(e)
        return bbd

    @cached_method
    def intersection_matrix(self):
        """
        Returns the intersection matrix of this ChowRing.

        OUTPUT:

        - A matrix, the intersection matrix of this ChowRing

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^3')  # P2
            sage: A.set_point_class('h^2')
            sage: A.intersection_matrix()
            [0 0 1]
            [0 1 0]
            [1 0 0]

            sage: A.<H, E> = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: A.basis()
            [E^2, E, H, 1]
            sage: A.set_point_class(H^2)
            sage: A.intersection_matrix()
            [ 0  0  0 -1]
            [ 0 -1  0  0]
            [ 0  0  1  0]
            [-1  0  0  0]

        TESTS::

            sage: A = ChowRing()
            sage: A.set_point_class(2)
            sage: A.intersection_matrix()
            [1/2]
        """
        rb, n = self.basis(), len(self.basis())
        if self.point_class() is None:
            raise ValueError("Need a point_class before calling integral.")
        m = [self(rb[i] * rb[j]).integral()
             for i in range(n) for j in range(n)]
        return matrix(QQ, n, n, m)

    def dual_basis_slow(self):
        """
        Returns the basis dual to self.basis() with respect to the intersection
        matrix. It is computed by simply inverting the intersection matrix.
        For a large ring basis, this method is somehow slow and has been
        optimized in self.dual_basis().

        OUTPUT:

        - A list of ring elements representing the basis dual to self.basis()

        EXAMPLES::

            sage: A = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: A.set_point_class('H^2')
            sage: A.dual_basis_slow()
            [-1, -E, H, -E^2]

        TESTS::

            sage: A = ChowRing()
            sage: A.set_point_class(2)
            sage: A.intersection_matrix()
            [1/2]
        """
        rb = self.basis()
        QMinv = self.intersection_matrix().inverse()
        return list(QMinv * vector(rb))

    @cached_method
    def dual_basis(self, verbose=False):
        r"""
        Returns the basis dual to :meth:`basis` with respect to the intersection
        matrix. We carefully compute only the relevant part of the intersection
        matrix in order to invert several smaller matrices instead of one large
        matrix as in :meth:`dual_basis_slow`. This method is mainly needed to
        compute `f_*` for morphisms `f:X\rightarrow Y` hence
        especially while computing a blowup and its tangent bundle.

        Remark: This requires quite some computational time for large examples,
        so we print out which block we compute if verbose is True.

        INPUT:

        - ``verbose`` -- an (optional) flag for printing intermediate steps

        OUTPUT:

        - A list of ring elements representing the basis dual to :meth:`basis`.

        EXAMPLES::

            sage: A.<H, E> = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: A.set_point_class(H^2)
            sage: A.dual_basis(verbose=True)
            Computing block 0 of size 1 ...
            Computing block 1 of size 2 ...
            [-1, -E, H, -E^2]

        TESTS::

            sage: A = ChowRing()
            sage: A.set_point_class(2)
            sage: A.dual_basis()
            [1/2]
        """
        dual_basis_dict = {}
        if self.point_class() is None:
            raise ValueError("Need a point_class before calling dual_basis.")
        R, pc = self.cover_ring(), self.point_class().lift() # type: ignore
        if pc.degree() < self.max_degree():
            raise ValueError("Expect the point_class in at least in max degree")
        if self.ngens() == 0:
            return [self.basis()[0] / pc]
        division = singular_function('division')
        #  We will compute the dictionary degree by degree.
        rbbd, maxdeg = self.basis_by_degree(), len(self.basis_by_degree())
        for d in range(1 + maxdeg // 2):
            # Ring basis in degree d and complementary degree maxdeg - d - 1
            rbd, rbdc, n = rbbd[d], rbbd[maxdeg - d - 1], len(rbbd[d])
            if verbose:
                print("Computing block %s of size %s ..." % (d, n))
            # Compute the intersection matrix relative to degree d
            m = [QQ(division(self(rbd[i] * rbdc[j]).lift(),
                             pc, ring=R)[0][0, 0])
                 for i in range(n) for j in range(n)]
            QMd = matrix(QQ, n, n, m)
            QMd_inverse = QMd.inverse()
            dbd = list(QMd_inverse.transpose() * vector(rbdc))
            db_dict = dict([(rbd[i], dbd[i]) for i in range(n)])
            # Add to dictionary
            dual_basis_dict.update(db_dict)
            if not ((d == maxdeg / 2) and (maxdeg % 2 == 1)):
                # Compute the dual basis for elements of rbdc
                # This is not needed if maxdeg is odd and d is maximal
                # as then rdb = rbdc.
                dbdc = list(QMd_inverse * vector(rbd))
                dbc_dict = dict([(rbdc[i], dbdc[i]) for i in range(n)])
                dual_basis_dict.update(dbc_dict)
        return [dual_basis_dict[e] for e in self.basis()]

    def dual_basis_dict(self, verbose=False):
        r"""
        Returns the dictionary `e_i\rightarrow e_i^*` where `e_i` is the basis
        given by :meth:`basis` and `e_i^*` is the dual element.
        If verbose is set to True, intermediate steps are printing during the
        calculation of the dual basis.

        INPUT:

        - ``verbose`` -- an (optional) flag for printing intermediate steps

        OUTPUT:

        - A dictionary associating `e_i\rightarrow e_i^*`.

        EXAMPLES::

            sage: A = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: A.set_point_class('H^2')
            sage: DualBasis = A.dual_basis_dict()
            sage: [(key, DualBasis[key]) for key in sorted(DualBasis)]
            [(1, -E^2), (E, -E), (H, H), (E^2, -1)]
        """
        rb, db = self.basis(), self.dual_basis(verbose=verbose)
        return dict([(rb[i], db[i]) for i in range(len(rb))])


PointChowRing = ChowRing()
PointChowRing.set_dimension(0)
PointChowRing.set_point_class(1)


# Helpers
def ideals_are_equal(I, J, A):
    r"""
    Internal Helper. Fast check that `I\subset J` and `J\subset I` in A using
    Singular.

    INPUT:

    - ``I`` -- an ideal in the ring ``A``
    - ``J`` -- an ideal in the ring ``A``
    - ``A`` -- a ring

    OUTPUT:
        True or False depending whether I = J in A or not.

    TESTS::

        sage: from sage.schemes.chow.ring import ideals_are_equal
        sage: R = PolynomialRing(QQ, 3, names=['a', 'b', 'c'])
        sage: a, b, c = R.gens()
        sage: I = R.ideal(a^2 - b, a^3)
        sage: J = R.ideal(b^2 , a * b, a^2 - b)
        sage: ideals_are_equal(I, J, R)
        True

        sage: from sage.schemes.chow.ring import ideals_are_equal
        sage: K = R.ideal(a^2 - b, a^2)
        sage: ideals_are_equal(I, K, R)
        False

    """
    singular_reduce = singular_function('reduce')
    std = singular_function('std')
    # Ensure I and J are in standard form.
    I, J = std(I, ring=A), std(J, ring=A)
    # Check I \subset J
    L = singular_reduce(I, J, ring=A, attributes={J: {'isSB': 1}})
    for x in L:  # Could call simplify, but the following is equally fast.
        if str(x) != '0':
            return False
    # Check J \subset I
    L = singular_reduce(J, I, ring=A, attributes={I: {'isSB': 1}})
    for x in L:  # Could call simplify, but the following is equally fast.
        if str(x) != '0':
            return False
    return True


def Kernel(f):
    r"""
    Given a morphism `f:A\rightarrow B` between (quotients of) multivariate
    polynomial rings, the kernel of `f` is returned.

    EXAMPLE (http://trac.sagemath.org/sage_trac/ticket/9792)::

        sage: from sage.schemes.chow.ring import Kernel
        sage: R = PolynomialRing(QQ, 4, names=['x', 'y', 'z', 'w'])
        sage: S = PolynomialRing(QQ, 2, names=['s', 't'])
        sage: s, t = S.gens()
        sage: f = R.hom([s ** 4, s ** 3 * t, s * t ** 3, t ** 4], S)
        sage: Kernel(f).gens()
        [y*z - x*w, z^3 - y*w^2, x*z^2 - y^2*w, y^3 - x^2*z]

    TESTS (taken from the Singular documentation of alg_kernel and kernel)::

        sage: from sage.schemes.chow.ring import Kernel
        sage: R = PolynomialRing(QQ, 3, names=['a', 'b', 'c'])
        sage: S = PolynomialRing(QQ, 6, names=['x', 'y', 'z', 'u', 'v', 'w'])
        sage: x, y, z, u, v, w = S.gens()
        sage: f = R.hom([x - w, u^2 * w + 1, y * z - v], S)
        sage: Kernel(f).gens()
        [0]

        sage: from sage.schemes.chow.ring import Kernel
        sage: R = PolynomialRing(QQ, 3, names=['a', 'b', 'c'])
        sage: S = PolynomialRing(QQ, 3, names=['x', 'y', 'z']).quotient(x - y)
        sage: x, y, z = S.gens()
        sage: f = R.hom([x, y, x^2 - y^3], S)
        sage: Kernel(f).gens()
        [a - b, b^3 - b^2 + c]

        sage: from sage.schemes.chow.ring import Kernel
        sage: R = PolynomialRing(QQ, 4, names=['a', 'b', 'c', 'd'])
        sage: S = PolynomialRing(QQ, 3, names=['a', 'b', 'c'])
        sage: a, b, c = S.gens()
        sage: f = R.hom([a, b, c, 0], S)
        sage: Kernel(f).gens()
        [d]
    """
    #######################################################################
    # Get rings A, B from f such that f: A -> B and lift to cover rings:
    # ff: AA -> BB with A = AA/AI and B = BB/BI
    #######################################################################
    singular_lib('algebra.lib')
    sing_alg_kernel = singular_function('alg_kernel')
    A, B = f.domain(), f.codomain()
    AA = A
    BB, BI, = B, B.ideal(0)
    ff = f
    if isinstance(A, QuotientRing_generic):
        AA = A.cover_ring()
        ff = f.morphism_from_cover()
    if isinstance(B, QuotientRing_generic):
        BB, BI = B.cover_ring(), B.defining_ideal()
        ff = AA.hom([BB(str(x.lift())) for x in ff.im_gens()], BB)
    s = 0 if BI.is_zero() else BI.ngens()

    # Two choices: either BI is trivial or not.
    if s:
        # The codomain is a quotient ring.
        # We consider AA[u_1, ..., u_s] -> BB defined by ff and BI.gens()
        # then take its kernel and finally set x_i=0 to get the kernel of ff.
        a_vars = [str(v) for v in AA.gens()]
        # Find a new variable name not yet used in AA:
        var_name = 'u'
        t_vars = ['u%d' % i for i in range(1, s + 1)]
        while [c for c in t_vars if c in a_vars]:
            var_name += var_name
            t_vars = ['%s%d' % (var_name, i) for i in range(1, s + 1)]
        # Construct R = AA[u_1, ..., u_s]
        r_vars = t_vars + a_vars  # We always place the new vars "on the left"!
        r_term_order = TermOrder('wdegrevlex', (1,) * s) + AA.term_order()
        R = PolynomialRing(QQ, len(r_vars), names=r_vars, order=r_term_order)
        # Construct g: R -> BB
        g = R.hom(BI.gens() + ff.im_gens(), BB)
        kern = sing_alg_kernel(BB.ideal(g.im_gens()), R, ring=BB).split(',')
        kern_g = R.ideal(kern)
        subs_dict = dict([(R(t), R(0)) for t in t_vars]) # type: ignore
        kern_ff = AA.ideal(kern_g.subs(subs_dict).gens())
    else:
        # The codomain is a Multivariate Polynomial ring.
        # Hence it is enough to take the kernel of ff:AA->BB
        kern = sing_alg_kernel(BB.ideal(ff.im_gens()), AA, ring=BB).split(',')
        kern_ff = AA.ideal(kern)

    # Return the kernel in standard basis as an ideal of A.
    sing_std = singular_function('std')
    return A.ideal(sing_std(AA.ideal(kern_ff), ring=AA))


def Kernel2(f):
    r"""
    Given a morphism `f:A\rightarrow B` between (quotients of) multivariate
    polynomial rings, the kernel of `f` is returned.

    EXAMPLE (http://trac.sagemath.org/sage_trac/ticket/9792)::

        sage: from sage.schemes.chow.ring import Kernel
        sage: R = PolynomialRing(QQ, 4, names=['x', 'y', 'z', 'w'])
        sage: S = PolynomialRing(QQ, 2, names=['s', 't'])
        sage: s, t = S.gens()
        sage: f = R.hom([s ** 4, s ** 3 * t, s * t ** 3, t ** 4], S)
        sage: Kernel(f).gens()
        [y*z - x*w, z^3 - y*w^2, x*z^2 - y^2*w, y^3 - x^2*z]

    TESTS (taken from the Singular documentation of alg_kernel and kernel)::

        sage: from sage.schemes.chow.ring import Kernel
        sage: R = PolynomialRing(QQ, 3, names=['a', 'b', 'c'])
        sage: S = PolynomialRing(QQ, 6, names=['x', 'y', 'z', 'u', 'v', 'w'])
        sage: x, y, z, u, v, w = S.gens()
        sage: f = R.hom([x - w, u^2 * w + 1, y * z - v], S)
        sage: Kernel(f).gens()
        [0]

        sage: from sage.schemes.chow.ring import Kernel
        sage: R = PolynomialRing(QQ, 3, names=['a', 'b', 'c'])
        sage: S = PolynomialRing(QQ, 3, names=['x', 'y', 'z']).quotient(x - y)
        sage: x, y, z = S.gens()
        sage: f = R.hom([x, y, x^2 - y^3], S)
        sage: Kernel(f).gens()
        [a - b, b^3 - b^2 + c]

        sage: from sage.schemes.chow.ring import Kernel
        sage: R = PolynomialRing(QQ, 4, names=['a', 'b', 'c', 'd'])
        sage: S = PolynomialRing(QQ, 3, names=['a', 'b', 'c'])
        sage: a, b, c = S.gens()
        sage: f = R.hom([a, b, c, 0], S)
        sage: Kernel(f).gens()
        [d]
    """
    A = f.domain()
    B = f.codomain()
    g = f.morphism_from_cover() if isinstance(A, QuotientRing_generic) else f
    f_sing = ",".join(str(gen) for gen in g.im_gens())
    from sage.all import singular
    singular.set_ring(singular(A))
    A_sing_name = singular.current_ring_name()
    singular.set_ring(singular(B))
    B_sing_name = singular.current_ring_name()
    singular.eval('''
        setring %s;
        ideal f=%s;
        setring %s;
        option(redSB);
        def kern = std(kernel(%s,f));
    ''' % (B_sing_name, f_sing, A_sing_name, B_sing_name))
    return A.ideal([A(k) for k in singular.get('kern').split(',')])
