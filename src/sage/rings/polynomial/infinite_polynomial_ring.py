r"""
Infinite Polynomial Rings

By Infinite Polynomial Rings, we mean polynomial rings in a countably
infinite number of variables. The implementation consists of a wrapper
around the current *finite* polynomial rings in Sage.

AUTHORS:

- Simon King <simon.king@nuigalway.ie>
- Mike Hansen <mhansen@gmail.com>

An Infinite Polynomial Ring has finitely many generators `x_\ast,
y_\ast,...` and infinitely many variables of the form `x_0, x_1, x_2,
..., y_0, y_1, y_2,...,...`.  We refer to the natural number `n` as
the *index* of the variable `x_n`.

INPUT:

- ``R`` -- the base ring; it has to be a commutative ring, and in some
  applications it must even be a field
- ``names`` -- a finite list of generator names; generator names must be alpha-numeric
- ``order`` -- (optional) string; the default order is ``'lex'`` (lexicographic).
  ``'deglex'`` is degree lexicographic, and ``'degrevlex'`` (degree reverse
  lexicographic) is possible but discouraged.

Each generator ``x`` produces an infinite sequence of variables
``x[1], x[2], ...`` which are printed on screen as ``x_1, x_2, ...``
and are latex typeset as `x_{1}, x_{2}`.  Then, the Infinite
Polynomial Ring is formed by polynomials in these variables.

By default, the monomials are ordered lexicographically. Alternatively,
degree (reverse) lexicographic ordering is possible as well. However, we
do not guarantee that the computation of Groebner bases will terminate
in this case.

In either case, the variables of a Infinite Polynomial Ring X are ordered
according to the following rule:

  ``X.gen(i)[m] > X.gen(j)[n]`` if and only if ``i<j or (i==j and m>n)``

We provide a 'dense' and a 'sparse' implementation. In the dense
implementation, the Infinite Polynomial Ring carries a finite
polynomial ring that comprises *all* variables up to the maximal index
that has been used so far. This is potentially a very big ring and may
also comprise many variables that are not used.

In the sparse implementation, we try to keep the underlying finite
polynomial rings small, using only those variables that are really
needed. By default, we use the dense implementation, since it usually
is much faster.

EXAMPLES::

    sage: X.<x,y> = InfinitePolynomialRing(ZZ, implementation='sparse')
    sage: A.<alpha,beta> = InfinitePolynomialRing(QQ, order='deglex')

    sage: f = x[5] + 2; f
    x_5 + 2
    sage: g = 3*y[1]; g
    3*y_1

It has some advantages to have an underlying ring that is not
univariate.  Hence, we always have at least two variables::

    sage: g._p.parent()
    Multivariate Polynomial Ring in y_1, y_0 over Integer Ring

    sage: f2 = alpha[5] + 2; f2
    alpha_5 + 2
    sage: g2 = 3*beta[1]; g2
    3*beta_1
    sage: A.polynomial_ring()
    Multivariate Polynomial Ring in alpha_5, alpha_4, alpha_3, alpha_2, alpha_1, alpha_0,
     beta_5, beta_4, beta_3, beta_2, beta_1, beta_0 over Rational Field

Of course, we provide the usual polynomial arithmetic::

    sage: f + g
    x_5 + 3*y_1 + 2
    sage: p = x[10]^2*(f+g); p
    x_10^2*x_5 + 3*x_10^2*y_1 + 2*x_10^2
    sage: p2 = alpha[10]^2*(f2+g2); p2
    alpha_10^2*alpha_5 + 3*alpha_10^2*beta_1 + 2*alpha_10^2

There is a permutation action on the variables, by permuting positive
variable indices::

    sage: P = Permutation(((10,1)))
    sage: p^P
    x_5*x_1^2 + 3*x_1^2*y_10 + 2*x_1^2
    sage: p2^P
    alpha_5*alpha_1^2 + 3*alpha_1^2*beta_10 + 2*alpha_1^2

Note that `x_0^P = x_0`, since the permutations only change *positive*
variable indices.

We also implemented ideals of Infinite Polynomial Rings. Here, it is
thoroughly assumed that the ideals are set-wise invariant under the
permutation action. We therefore refer to these ideals as *Symmetric
Ideals*. Symmetric Ideals are finitely generated modulo addition,
multiplication by ring elements and permutation of variables. If the
base ring is a field, one can compute Symmetric Groebner Bases::

    sage: J = A * (alpha[1]*beta[2])
    sage: J.groebner_basis()                                                            # needs sage.combinat sage.libs.singular
    [alpha_1*beta_2, alpha_2*beta_1]

For more details, see :class:`~sage.rings.polynomial.symmetric_ideal.SymmetricIdeal`.

Infinite Polynomial Rings can have any commutative base ring. If the
base ring of an Infinite Polynomial Ring is a (classical or infinite)
Polynomial Ring, then our implementation tries to merge everything
into *one* ring. The basic requirement is that the monomial orders
match. In the case of two Infinite Polynomial Rings, the
implementations must match. Moreover, name conflicts should be
avoided. An overlap is only accepted if the order of variables can be
uniquely inferred, as in the following example::

    sage: A.<a,b,c> = InfinitePolynomialRing(ZZ)
    sage: B.<b,c,d> = InfinitePolynomialRing(A)
    sage: B
    Infinite polynomial ring in a, b, c, d over Integer Ring

This is also allowed if finite polynomial rings are involved::

    sage: A.<a_3,a_1,b_1,c_2,c_0> = ZZ[]
    sage: B.<b,c,d> = InfinitePolynomialRing(A, order='degrevlex')
    sage: B
    Infinite polynomial ring in b, c, d over
     Multivariate Polynomial Ring in a_3, a_1 over Integer Ring

It is no problem if one generator of the Infinite Polynomial Ring is
called ``x`` and one variable of the base ring is also called
``x``. This is since no *variable* of the Infinite Polynomial Ring
will be called ``x``. However, a problem arises if the underlying
classical Polynomial Ring has a variable ``x_1``, since this can be
confused with a variable of the Infinite Polynomial Ring. In this
case, an error will be raised::

    sage: X.<x,y_1> = ZZ[]
    sage: Y.<x,z> = InfinitePolynomialRing(X)

Note that ``X`` is not merged into ``Y``; this is since the monomial
order of ``X`` is 'degrevlex', but of ``Y`` is 'lex'.
::

    sage: Y
    Infinite polynomial ring in x, z over
     Multivariate Polynomial Ring in x, y_1 over Integer Ring

The variable ``x`` of ``X`` can still be interpreted in ``Y``,
although the first generator of ``Y`` is called ``x`` as well::

    sage: x
    x_*
    sage: X('x')
    x
    sage: Y(X('x'))
    x
    sage: Y('x')
    x

But there is only merging if the resulting monomial order is uniquely
determined. This is not the case in the following examples, and thus
an error is raised::

    sage: X.<y_1,x> = ZZ[]
    sage: Y.<y,z> = InfinitePolynomialRing(X)
    Traceback (most recent call last):
    ...
    CoercionException: Overlapping variables (('y', 'z'),['y_1']) are incompatible
    sage: Y.<z,y> = InfinitePolynomialRing(X)
    Traceback (most recent call last):
    ...
    CoercionException: Overlapping variables (('z', 'y'),['y_1']) are incompatible
    sage: X.<x_3,y_1,y_2> = PolynomialRing(ZZ, order='lex')
    sage: # y_1 and y_2 would be in opposite order in an Infinite Polynomial Ring
    sage: Y.<y> = InfinitePolynomialRing(X)
    Traceback (most recent call last):
    ...
    CoercionException: Overlapping variables (('y',),['y_1', 'y_2']) are incompatible


If the type of monomial orderings (e.g., 'degrevlex' versus 'lex') or
if the implementations do not match, there is no simplified
construction available::

    sage: X.<x,y> = InfinitePolynomialRing(ZZ)
    sage: Y.<z> = InfinitePolynomialRing(X, order='degrevlex')
    sage: Y
    Infinite polynomial ring in z over Infinite polynomial ring in x, y over Integer Ring
    sage: Y.<z> = InfinitePolynomialRing(X, implementation='sparse')
    sage: Y
    Infinite polynomial ring in z over Infinite polynomial ring in x, y over Integer Ring

TESTS:

Infinite Polynomial Rings are part of Sage's coercion system. Hence,
we can do arithmetic, so that the result lives in a ring into which
all constituents coerce.
::

    sage: R.<a,b> = InfinitePolynomialRing(ZZ)
    sage: X.<x> = InfinitePolynomialRing(R)
    sage: x[2]/2+(5/3)*a[3]*x[4] + 1
    5/3*a_3*x_4 + 1/2*x_2 + 1

    sage: R.<a,b> = InfinitePolynomialRing(ZZ, implementation='sparse')
    sage: X.<x> = InfinitePolynomialRing(R)
    sage: x[2]/2+(5/3)*a[3]*x[4] + 1
    5/3*a_3*x_4 + 1/2*x_2 + 1

    sage: R.<a,b> = InfinitePolynomialRing(ZZ, implementation='sparse')
    sage: X.<x> = InfinitePolynomialRing(R, implementation='sparse')
    sage: x[2]/2+(5/3)*a[3]*x[4] + 1
    5/3*a_3*x_4 + 1/2*x_2 + 1

    sage: R.<a,b> = InfinitePolynomialRing(ZZ)
    sage: X.<x> = InfinitePolynomialRing(R, implementation='sparse')
    sage: x[2]/2+(5/3)*a[3]*x[4] + 1
    5/3*a_3*x_4 + 1/2*x_2 + 1

Check that :issue:`22514` is fixed::

    sage: R.<x> = InfinitePolynomialRing(ZZ)
    sage: a = R(3)
    sage: a.is_constant()
    True
    sage: a.constant_coefficient()
    3
    sage: a.degree()
    0
    sage: b = R("2")
    sage: b.parent() is R
    True
    sage: S.<y> = ZZ[]
    sage: Q.<z> = InfinitePolynomialRing(S)
    sage: a = Q(1+y)
    sage: a.is_constant()
    True
    sage: a.constant_coefficient()
    y + 1
"""
# ****************************************************************************
#       Copyright (C) 2009 Simon King <simon.king@nuigalway.ie> and
#                          Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import operator
import re
from functools import reduce

from sage.rings.ring import CommutativeRing
from sage.categories.rings import Rings
from sage.structure.all import SageObject, parent
from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method

###################################################
#  The Construction Functor

from sage.categories.pushout import InfinitePolynomialFunctor


###############################################################
# Ring Factory framework

class InfinitePolynomialRingFactory(UniqueFactory):
    """
    A factory for creating infinite polynomial ring elements.  It
    makes sure that they are unique as well as handling pickling.
    For more details, see
    :class:`~sage.structure.factory.UniqueFactory` and
    :mod:`~sage.rings.polynomial.infinite_polynomial_ring`.

    EXAMPLES::

        sage: A.<a> = InfinitePolynomialRing(QQ)
        sage: B.<b> = InfinitePolynomialRing(A)
        sage: B.construction()
        [InfPoly{[a,b], "lex", "dense"}, Rational Field]
        sage: R.<a,b> = InfinitePolynomialRing(QQ)
        sage: R is B
        True
        sage: X.<x> = InfinitePolynomialRing(QQ)
        sage: X2.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
        sage: X is X2
        False

        sage: X is loads(dumps(X))
        True
    """
    def create_key(self, R, names=('x',), order='lex', implementation='dense'):
        """
        Create a key which uniquely defines the infinite polynomial ring.

        TESTS::

            sage: InfinitePolynomialRing.create_key(QQ, ('y1',))
            (InfPoly{[y1], "lex", "dense"}(FractionField(...)), Integer Ring)
            sage: _[0].all
            [FractionField, InfPoly{[y1], "lex", "dense"}]
            sage: InfinitePolynomialRing.create_key(QQ, names=['beta'], order='deglex',
            ....:                                   implementation='sparse')
            (InfPoly{[beta], "deglex", "sparse"}(FractionField(...)), Integer Ring)
            sage: _[0].all
            [FractionField, InfPoly{[beta], "deglex", "sparse"}]
            sage: InfinitePolynomialRing.create_key(QQ, names=['x','y'],
            ....:                                   implementation='dense')
            (InfPoly{[x,y], "lex", "dense"}(FractionField(...)), Integer Ring)
            sage: _[0].all
            [FractionField, InfPoly{[x,y], "lex", "dense"}]

        If no generator name is provided, a generator named 'x',
        lexicographic order and the dense implementation are assumed::

            sage: InfinitePolynomialRing.create_key(QQ)
            (InfPoly{[x], "lex", "dense"}(FractionField(...)), Integer Ring)
            sage: _[0].all
            [FractionField, InfPoly{[x], "lex", "dense"}]

        If it is attempted to use no generator, a :exc:`ValueError` is raised::

            sage: InfinitePolynomialRing.create_key(ZZ, names=[])
            Traceback (most recent call last):
            ...
            ValueError: Infinite Polynomial Rings must have at least one generator
        """
        if isinstance(names, list):
            names = tuple(names)
        if not names:
            raise ValueError("Infinite Polynomial Rings must have at least one generator")
        if len(names) > len(set(names)):
            raise ValueError("the variable names must be distinct")
        F = InfinitePolynomialFunctor(names, order, implementation)
        while hasattr(R, 'construction'):
            C = R.construction()
            if C is None:
                break
            F = F * C[0]
            R = C[1]
        return (F, R)

    def create_object(self, version, key):
        """
        Return the infinite polynomial ring corresponding to the key ``key``.

        TESTS::

            sage: InfinitePolynomialRing.create_object('1.0', InfinitePolynomialRing.create_key(ZZ, ('x3',)))
            Infinite polynomial ring in x3 over Integer Ring
        """
        if len(key) > 2:
            # We got an old pickle. By calling the ring constructor, it will automatically
            # be transformed into the new scheme
            return InfinitePolynomialRing(*key)
        # By now, we have different unique keys, based on construction functors
        C, R = key
        from sage.categories.pushout import CompositeConstructionFunctor, InfinitePolynomialFunctor
        if isinstance(C, CompositeConstructionFunctor):
            F = C.all[-1]
            if len(C.all) > 1:
                R = CompositeConstructionFunctor(*C.all[:-1])(R)
        else:
            F = C
        if not isinstance(F, InfinitePolynomialFunctor):
            raise TypeError("we expected an InfinitePolynomialFunctor, not %s" % type(F))
        if F._imple == 'sparse':
            return InfinitePolynomialRing_sparse(R, F._gens, order=F._order)
        return InfinitePolynomialRing_dense(R, F._gens, order=F._order)


InfinitePolynomialRing = InfinitePolynomialRingFactory('InfinitePolynomialRing')


##############################################################
#  An auxiliary dictionary-like class that returns variables

class InfiniteGenDict:
    """
    A dictionary-like class that is suitable for usage in ``sage_eval``.

    The generators of an Infinite Polynomial Ring are not
    variables. Variables of an Infinite Polynomial Ring are returned
    by indexing a generator. The purpose of this class is to return a
    variable of an Infinite Polynomial Ring, given its string
    representation.

    EXAMPLES::

        sage: R.<a,b> = InfinitePolynomialRing(ZZ)
        sage: D = R.gens_dict() # indirect doctest
        sage: D._D
        [InfiniteGenDict defined by ['a', 'b'], {'1': 1}]
        sage: D._D[0]['a_15']
        a_15
        sage: type(_)
        <class 'sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_dense'>
        sage: sage_eval('3*a_3*b_5-1/2*a_7', D._D[0])
        -1/2*a_7 + 3*a_3*b_5
    """
    def __init__(self, Gens):
        """
        INPUT:

        - ``Gens`` -- list of generators of an infinite polynomial ring

        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: D = R.gens_dict() # indirect doctest
            sage: D._D
            [InfiniteGenDict defined by ['a', 'b'], {'1': 1}]
            sage: D._D == loads(dumps(D._D)) # indirect doctest
            True
        """
        self._D = dict(zip(((hasattr(X, '_name') and X._name) or repr(X) for X in Gens), Gens))

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: D = R.gens_dict() # indirect doctest
            sage: D._D
            [InfiniteGenDict defined by ['a', 'b'], {'1': 1}]
            sage: D._D == loads(dumps(D._D)) # indirect doctest
            True
        """
        if isinstance(other, InfiniteGenDict):
            return self._D == other._D
        return False

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: D = R.gens_dict() # indirect doctest
            sage: D._D
            [InfiniteGenDict defined by ['a', 'b'], {'1': 1}]
            sage: D._D != loads(dumps(D._D)) # indirect doctest
            False
        """
        return not (self == other)

    def __repr__(self):
        """
        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: D = R.gens_dict()
            sage: D._D # indirect doctest
            [InfiniteGenDict defined by ['a', 'b'], {'1': 1}]
        """
        return "InfiniteGenDict defined by %s" % sorted(self._D)

    def __getitem__(self, k):
        """
        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: D = R.gens_dict() # indirect doctest
            sage: D._D
            [InfiniteGenDict defined by ['a', 'b'], {'1': 1}]
            sage: D._D[0]['a_15']
            a_15
            sage: type(_)
            <class 'sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_dense'>
        """
        if not isinstance(k, str):
            raise KeyError("string expected")
        L = k.split('_')
        try:
            if len(L) == 2:
                return self._D[L[0]][int(L[1])]
        except (TypeError, ValueError):
            pass
        raise KeyError("%s is not a variable name" % k)


class GenDictWithBasering:
    """
    A dictionary-like class that is suitable for usage in ``sage_eval``.

    This pseudo-dictionary accepts strings as index, and then walks down
    a chain of base rings of (infinite) polynomial rings until it finds
    one ring that has the given string as variable name, which is then
    returned.

    EXAMPLES::

        sage: R.<a,b> = InfinitePolynomialRing(ZZ)
        sage: D = R.gens_dict() # indirect doctest
        sage: D
        GenDict of Infinite polynomial ring in a, b over Integer Ring
        sage: D['a_15']
        a_15
        sage: type(_)
        <class 'sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_dense'>
        sage: sage_eval('3*a_3*b_5-1/2*a_7', D)
        -1/2*a_7 + 3*a_3*b_5
    """
    def __init__(self, parent, start):
        """
        INPUT:

        - ``parent`` -- a ring
        - ``start`` -- some dictionary, usually the dictionary of variables of ``parent``

        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: D = R.gens_dict() # indirect doctest
            sage: D
            GenDict of Infinite polynomial ring in a, b over Integer Ring
            sage: D['a_15']
            a_15
            sage: type(_)
            <class 'sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_dense'>
            sage: sage_eval('3*a_3*b_5-1/2*a_7', D)
            -1/2*a_7 + 3*a_3*b_5

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring import GenDictWithBasering
            sage: R = ZZ['x']['y']['a','b']['c']
            sage: D = GenDictWithBasering(R,R.gens_dict())
            sage: R.gens_dict()['a']
            Traceback (most recent call last):
            ...
            KeyError: 'a'
            sage: D['a']
            a
        """
        P = self._P = parent
        if isinstance(start, list):
            self._D = start
            return
        self._D = [start]
        while hasattr(P, 'base_ring') and (P.base_ring() is not P):
            P = P.base_ring()
            D = P.gens_dict()
            if isinstance(D, GenDictWithBasering):
                self._D.extend(D._D)
                break
            else:
                self._D.append(D)

    def __next__(self):
        """
        Return a dictionary that can be used to interpret strings in the base ring of ``self``.

        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(QQ['t'])
            sage: D = R.gens_dict()
            sage: D
            GenDict of Infinite polynomial ring in a, b over Univariate Polynomial Ring in t over Rational Field
            sage: next(D)
            GenDict of Univariate Polynomial Ring in t over Rational Field
            sage: sage_eval('t^2', next(D))
            t^2
        """
        if len(self._D) <= 1:
            raise ValueError("no next term for %s available" % self)
        return GenDictWithBasering(self._P.base_ring(), self._D[1:])

    next = __next__

    def __repr__(self):
        """
        TESTS::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: D = R.gens_dict() # indirect doctest
            sage: D
            GenDict of Infinite polynomial ring in a, b over Integer Ring
        """
        return "GenDict of " + repr(self._P)

    def __getitem__(self, k):
        """
        TESTS::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: D = R.gens_dict() # indirect doctest
            sage: D
            GenDict of Infinite polynomial ring in a, b over Integer Ring
            sage: D['a_15']
            a_15
            sage: type(_)
            <class 'sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_dense'>
        """
        for D in self._D:
            try:
                return D[k]
            except KeyError:
                pass
        raise KeyError("{} is not a variable name of {} or its iterated base rings".format(k, self._P))


##############################################################
#  The sparse implementation

class InfinitePolynomialRing_sparse(CommutativeRing):
    r"""
    Sparse implementation of Infinite Polynomial Rings.

    An Infinite Polynomial Ring with generators `x_\ast, y_\ast,
    ...` over a field `F` is a free commutative `F`-algebra generated
    by `x_0, x_1, x_2, ..., y_0, y_1, y_2, ..., ...` and is equipped
    with a permutation action on the generators, namely `x_n^P =
    x_{P(n)}, y_{n}^P=y_{P(n)}, ...` for any permutation `P` (note
    that variables of index zero are invariant under such
    permutation).

    It is known that any permutation invariant ideal in an Infinite
    Polynomial Ring is finitely generated modulo the permutation
    action -- see :class:`~sage.rings.polynomial.symmetric_ideal.SymmetricIdeal`
    for more details.

    Usually, an instance of this class is created using
    ``InfinitePolynomialRing`` with the optional parameter
    ``implementation='sparse'``. This takes care of uniqueness of
    parent structures. However, a direct construction is possible, in
    principle::

        sage: X.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
        sage: Y.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
        sage: X is Y
        True
        sage: from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing_sparse
        sage: Z = InfinitePolynomialRing_sparse(QQ, ['x','y'], 'lex')

    Nevertheless, since infinite polynomial rings are supposed to be unique
    parent structures, they do not evaluate equal. ::

        sage: Z == X
        False

    The last parameter ('lex' in the above example) can also be
    'deglex' or 'degrevlex'; this would result in an Infinite
    Polynomial Ring in degree lexicographic or degree reverse
    lexicographic order.

    See :mod:`~sage.rings.polynomial.infinite_polynomial_ring` for
    more details.
    """
    def __init__(self, R, names, order):
        """
        INPUT:

        - ``R`` -- base ring
        - ``names`` -- list of generator names
        - ``order`` -- string determining the monomial order of the infinite polynomial ring

        EXAMPLES::

            sage: X.<alpha,beta> = InfinitePolynomialRing(ZZ, implementation='sparse')

        Infinite Polynomial Rings are unique parent structures::

            sage: X is loads(dumps(X))
            True
            sage: p=alpha[10]*beta[2]^3+2*alpha[1]*beta[3]
            sage: p
            alpha_10*beta_2^3 + 2*alpha_1*beta_3

        We define another Infinite Polynomial Ring with same generator
        names but a different order. These rings are different, but
        allow for coercion::

            sage: Y.<alpha,beta> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse')
            sage: Y is X
            False
            sage: q=beta[2]^3*alpha[10]+beta[3]*alpha[1]*2
            sage: q
            alpha_10*beta_2^3 + 2*alpha_1*beta_3
            sage: p==q
            True
            sage: X.gen(1)[2]*Y.gen(0)[1]
            alpha_1*beta_2
        """
        if not names:
            names = ['x']
        for n in names:
            if not (isinstance(n, str) and n.isalnum() and (not n[0].isdigit())):
                raise ValueError("generator names must be alphanumeric strings not starting with a digit, but %s is not" % n)
        if len(names) != len(set(names)):
            raise ValueError("generator names must be pairwise different")
        self._names = tuple(names)
        if not isinstance(order, str):
            raise TypeError("the monomial order must be given as a string")
        if R not in Rings().Commutative():
            raise TypeError("the given 'base ring' (= %s) must be a commutative ring" % R)

        # now, the input is accepted
        if hasattr(R, '_underlying_ring'):
            self._underlying_ring = R._underlying_ring
        else:
            self._underlying_ring = R.base_ring()

        # some tools to analyse polynomial string representations.
        self._identify_variable = lambda x, y: (-self._names.index(x), int(y))
        self._find_maxshift = re.compile('_([0-9]+)')  # findall yields stringrep of the shifts
        self._find_variables = re.compile('[a-zA-Z0-9]+_[0-9]+')
        self._find_varpowers = re.compile(r'([a-zA-Z0-9]+)_([0-9]+)\^?([0-9]*)')  # findall yields triple "generator_name", "index", "exponent"

        # Create some small underlying polynomial ring.
        # It is used to ensure that the parent of the underlying
        # polynomial of an element of self is actually a *multi*variate
        # polynomial ring.
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        if len(names) == 1:
            VarList = [names[0] + '_0', names[0] + '_1']
        else:
            VarList = [X + '_0' for X in names]
        VarList.sort(key=self.varname_key, reverse=True)
        self._minP = PolynomialRing(R, len(VarList), VarList)

        # some basic data
        self._order = order
        self._name_dict = {name: i for i, name in enumerate(names)}
        from sage.categories.commutative_algebras import CommutativeAlgebras
        CommutativeRing.__init__(self, R, category=CommutativeAlgebras(R))

        self._populate_coercion_lists_()

    def __repr__(self):
        """
        EXAMPLES::

            sage: InfinitePolynomialRing(QQ)         # indirect doctest
            Infinite polynomial ring in x over Rational Field

            sage: X.<alpha,beta> = InfinitePolynomialRing(ZZ, order='deglex'); X
            Infinite polynomial ring in alpha, beta over Integer Ring
        """
        return "Infinite polynomial ring in {} over {}".format(", ".join(self._names), self._base)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.misc.latex import latex
            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: latex(X)        # indirect doctest
            \Bold{Q}[x_{\ast}, y_{\ast}]
        """
        from sage.misc.latex import latex
        vars = ', '.join(latex(X) for X in self.gens())
        return f"{latex(self.base_ring())}[{vars}]"

    @cached_method
    def _an_element_(self):
        """
        Return an element of this ring.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(QQ)
            sage: R.an_element() # indirect doctest
            x_1
        """
        x = self.gen(0)
        return x[1]

    @cached_method
    def one(self):
        """
        TESTS::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: X.one()
            1
        """
        from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
        return InfinitePolynomial(self, self._base(1))

    #####################
    # coercion

    def construction(self):
        """
        Return the construction of ``self``.

        OUTPUT:

        A pair ``F,R``, where ``F`` is a construction functor and ``R`` is a ring,
        so that ``F(R) is self``.

        EXAMPLES::

            sage: R.<x,y> = InfinitePolynomialRing(GF(5))
            sage: R.construction()
            [InfPoly{[x,y], "lex", "dense"}, Finite Field of size 5]
        """
        return [InfinitePolynomialFunctor(self._names, self._order, 'sparse'), self._base]

    def _coerce_map_from_(self, S):
        r"""
        Coerce things into ``self``.

        NOTE:

        Any coercion will preserve variable names.

        EXAMPLES:

        Here, we check to see that elements of a *finitely* generated
        polynomial ring with appropriate variable names coerce
        correctly into the Infinite Polynomial Ring::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: px0 = PolynomialRing(QQ,'x_0').gen(0)
            sage: px0 + x[0]  # indirect doctest
            2*x_0
            sage: px0==x[0]
            True

        It is possible to construct an Infinite Polynomial Ring whose
        base ring is another Infinite Polynomial Ring::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: X.<x> = InfinitePolynomialRing(R)
            sage: a[2]*x[3]+x[1]*a[4]^2
            a_4^2*x_1 + a_2*x_3
        """
        # Use Construction Functors!
        from sage.categories.pushout import pushout
        try:
            # the following line should not test "pushout is self", but
            # only "pushout == self", since we also allow coercion from
            # dense to sparse implementation!
            P = pushout(self, S)
            # We do not care about the orders. But base ring and generators
            # of the pushout should remain the same as in self.
            return P._names == self._names and P._base == self._base
        except (TypeError, ValueError):
            return False

    def _element_constructor_(self, x):
        """
        Return an element of ``self``.

        INPUT:

        - ``x`` -- any object that can be interpreted in ``self``

        TESTS::

            sage: X = InfinitePolynomialRing(QQ)
            sage: a = X(2); a     # indirect doctest
            2
            sage: a.parent()
            Infinite polynomial ring in x over Rational Field

            sage: R = PolynomialRing(ZZ, ['x_3'])
            sage: b = X(R.gen()); b
            x_3
            sage: b.parent()
            Infinite polynomial ring in x over Rational Field
            sage: X('(x_1^2+2/3*x_4)*(x_2+x_5)')
            2/3*x_5*x_4 + x_5*x_1^2 + 2/3*x_4*x_2 + x_2*x_1^2

            sage: Y = InfinitePolynomialRing(ZZ)
            sage: Y('1/3')
            Traceback (most recent call last):
            ...
            ValueError: cannot convert 1/3 into an element of Infinite polynomial ring in x over Integer Ring

        .. WARNING::

            The :issue:`37756` is not yet fixed::

                sage: L.<x, y> = QQ[]
                sage: R.<a> = InfinitePolynomialRing(QQ)
                sage: M = InfinitePolynomialRing(L, names=["a"])
                sage: c = a[0]
                sage: M(c)  # known bug
                a_0
        """
        from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
        # In many cases, the easiest solution is to "simply" evaluate
        # the string representation.
        from sage.misc.sage_eval import sage_eval
        if isinstance(x, str):
            try:
                x = sage_eval(x, self.gens_dict())
            except (TypeError, ValueError, SyntaxError):
                raise ValueError(f"cannot convert {x} into an element of {self}")
            P = parent(x)
            if P is self:
                return x
            elif self._base.has_coerce_map_from(P):
                return InfinitePolynomial(self, self._base(x))
            else:
                raise ValueError(f"cannot convert {x} into an element of {self}")

        if isinstance(parent(x), InfinitePolynomialRing_sparse):
            # the easy case - parent == self - is already past
            if x.parent() is self._base:  # another easy case
                return InfinitePolynomial(self, x)
            xmaxind = x.max_index()  # save for later
            x = x._p
        else:
            xmaxind = -1

        # Now, we focus on the underlying classical polynomial ring.
        # First, try interpretation in the base ring.
        try:
            from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict
            if isinstance(self._base, MPolynomialRing_polydict):
                x = sage_eval(repr(), next(self.gens_dict()))
            else:
                x = self._base(x)
            # remark: Conversion to self._P (if applicable)
            # is done in InfinitePolynomial()
            return InfinitePolynomial(self, x)
        except (TypeError, ValueError):
            pass

        # By now, we can assume that x has a parent, because
        # types like int have already been done in the previous step;
        # and also it is not an InfinitePolynomial.
        # If it is not a polynomial (duck typing: we need
        # the variables attribute), we fall back to using strings
        if not hasattr(x, 'variables'):
            try:
                return sage_eval(repr(x), self.gens_dict())
            except (TypeError, ValueError, SyntaxError):
                raise ValueError(f"cannot convert {x} into an element of {self}")

        # direct conversion will only be used if the underlying polynomials are libsingular.
        from sage.rings.polynomial.multi_polynomial import MPolynomial_libsingular
        # try interpretation in self._P, if we have a dense implementation
        if hasattr(self, '_P'):
            if x.parent() is self._P:
                return InfinitePolynomial(self, x)
            # It's a shame to use sage_eval. However, it's even more of a shame
            # that MPolynomialRing_polydict does not work in complicated settings.
            # So, if self._P is libsingular (and this will be the case in many
            # applications!), we do it "nicely". Otherwise, we have to use sage_eval.
            if isinstance(x, MPolynomial_libsingular):
                from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
                if isinstance(self._P, MPolynomialRing_libsingular):
                    if xmaxind == -1:  # Otherwise, x has been an InfinitePolynomial
                        # We infer the correct variable shift.
                        # Note: Since we are in the "libsingular" case, there are
                        # no further "variables" hidden in the base ring of x.parent()
                        try:
                            VarList = [repr(v) for v in x.variables()]
                            # since interpretation in base ring
                            # was impossible, it *must* have
                            # variables
                            # This tests admissibility on the fly:
                            VarList.sort(key=self.varname_key, reverse=True)
                        except ValueError:
                            raise ValueError("cannot convert {} into an element of {} - variables are not admissible".format(x, self))
                        xmaxind = max([int(v.split('_')[1]) for v in VarList])
                    try:
                        # Apparently, in libsingular, the polynomial conversion is not done by
                        # name but by position, if the number of variables in the parents coincide.
                        # So, we shift self._P to achieve xmaxind, and if the number of variables is
                        # the same then we shift further. We then *must* be
                        # able to convert x into self._P, or conversion to self is
                        # impossible (and will be done in InfinitePolynomial(...)
                        if self._max < xmaxind:
                            self.gen()[xmaxind]
                        if self._P.ngens() == x.parent().ngens():
                            self.gen()[self._max + 1]
                        # conversion to self._P will be done in InfinitePolynomial.__init__
                        return InfinitePolynomial(self, x)
                    except (ValueError, TypeError, NameError):
                        raise ValueError("cannot convert {} (from {}, but variables {}) into an element of {} - no conversion into underlying polynomial ring {}".format(x, x.parent(), x.variables(), self, self._P))
            # By now, x or self._P are not libsingular. Since MPolynomialRing_polydict
            # is too buggy, we use string evaluation
            try:
                return sage_eval(repr(x), self.gens_dict())
            except (ValueError, TypeError, NameError):
                raise ValueError("cannot convert {} into an element of {} - no conversion into underlying polynomial ring".format(x, self))

        # By now, we are in the sparse case.
        try:
            VarList = [repr(v) for v in x.variables()]
            # since interpretation in base ring
            # was impossible, it *must* have
            # variables
            # This tests admissibility on the fly:
            VarList.sort(key=self.varname_key, reverse=True)
        except ValueError:
            raise ValueError("cannot convert {} into an element of {} - variables are not admissible".format(x, self))

        if len(VarList) == 1:
            # univariate polynomial rings are crab. So, make up another variable.
            if VarList[0] == self._names[0] + '_0':
                VarList.append(self._names[0] + '_1')
            else:
                VarList.append(self._names[0] + '_0')
        # We ensure that polynomial conversion is done by names;
        # the problem is that it is done by names if the number of variables coincides.
        if len(VarList) == x.parent().ngens():
            BigList = x.parent().variable_names()
            ind = 2
            while self._names[0] + '_' + str(ind) in BigList:
                ind += 1
            VarList.append(self._names[0] + '_' + str(ind))
        try:
            VarList.sort(key=self.varname_key, reverse=True)
        except ValueError:
            raise ValueError("cannot convert {} into an element of {}; the variables are not admissible".format(x, self))

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(self._base, VarList, order=self._order)
        if isinstance(x, MPolynomial_libsingular):  # everything else is so buggy that it's even not worth to try.
            from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            if isinstance(R, MPolynomialRing_libsingular):
                try:
                    # Problem: If there is only a partial overlap in the variables
                    # of x.parent() and R, then R(x) raises an error (which, I think,
                    # is a bug, since we talk here about conversion, not coercion).
                    # Hence, for being on the safe side, we coerce into a pushout ring:
                    x = R(1) * x
                    return InfinitePolynomial(self, x)
                except (TypeError, ValueError):
                    # OK, last resort, to be on the safe side
                    try:
                        return sage_eval(repr(x), self.gens_dict())
                    except (ValueError, TypeError, NameError):
                        raise ValueError("cannot convert {} into an element of {}; conversion of the underlying polynomial failed".format(x, self))
        try:
            return sage_eval(repr(x), self.gens_dict())
        except (ValueError, TypeError, NameError):
            raise ValueError(f"cannot convert {x} into an element of {self}")

    def tensor_with_ring(self, R):
        """
        Return the tensor product of ``self`` with another ring.

        INPUT:

        - ``R`` -- a ring

        OUTPUT:

        An infinite polynomial ring that, mathematically, can be seen as the
        tensor product of ``self`` with ``R``.

        NOTE:

        It is required that the underlying ring of ``self`` coerces into ``R``.
        Hence, the tensor product is in fact merely an extension of the base
        ring.

        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: R.tensor_with_ring(QQ)
            Infinite polynomial ring in a, b over Rational Field
            sage: R
            Infinite polynomial ring in a, b over Integer Ring

        The following tests against a bug that was fixed at :issue:`10468`::

            sage: R.<x,y> = InfinitePolynomialRing(QQ)
            sage: R.tensor_with_ring(QQ) is R
            True
        """
        if not R.has_coerce_map_from(self._underlying_ring):
            raise TypeError("we cannot tensor with " + repr(R))
        B = self.base_ring()
        if hasattr(B, 'tensor_with_ring'):
            return InfinitePolynomialRing(B.tensor_with_ring(R), self._names, self._order, implementation='sparse')
        if hasattr(B, 'change_ring'):  # e.g., polynomial rings
            return InfinitePolynomialRing(B.change_ring(R), self._names, self._order, implementation='sparse')
        # try to find the correct base ring in other ways:
        try:
            o = B.one() * R.one()
        except (TypeError, ValueError):
            raise TypeError("we cannot tensor with " + repr(R))
        return InfinitePolynomialRing(o.parent(), self._names, self._order, implementation='sparse')

    # Basic Ring Properties
    # -- some stuff that is useful for quotient rings etc.
    def is_noetherian(self):
        """
        Return ``False``, since polynomial rings in infinitely many
        variables are never Noetherian rings.

        Since Infinite Polynomial Rings must have at least one
        generator, they have infinitely many variables and are thus
        not Noetherian, as a ring.

        .. NOTE::

            Infinite Polynomial Rings over a field `F` are Noetherian as
            `F(G)` modules, where `G` is the symmetric group of the
            natural numbers. But this is not what the method
            ``is_noetherian()`` is answering.

        TESTS::

            sage: R = InfinitePolynomialRing(GF(2))
            sage: R
            Infinite polynomial ring in x over Finite Field of size 2
            sage: R.is_noetherian()
            False

            sage: R.<x> = InfinitePolynomialRing(QQ)
            sage: R.is_noetherian()
            False
        """
        return False

    def is_field(self, *args, **kwds):
        """
        Return ``False`` since Infinite Polynomial Rings are never fields.

        Since Infinite Polynomial Rings must have at least one generator,
        they have infinitely many variables and thus never are fields.

        EXAMPLES::

            sage: R.<x, y> = InfinitePolynomialRing(QQ)
            sage: R.is_field()
            False

        TESTS::

            sage: R = InfinitePolynomialRing(GF(2))
            sage: R
            Infinite polynomial ring in x over Finite Field of size 2
            sage: R.is_field()
            False

        :issue:`9443`::

            sage: W = PowerSeriesRing(InfinitePolynomialRing(QQ,'a'),'x')
            sage: W.is_field()
            False
        """
        return False

    def varname_key(self, x):
        """
        Key for comparison of variable names.

        INPUT:

        - ``x`` -- string of the form ``a+'_'+str(n)``, where a is the
          name of a generator, and n is an integer

        OUTPUT: a key used to sort the variables

        THEORY:

        The order is defined as follows:

        x<y `\\iff` the string ``x.split('_')[0]`` is later in the list of
        generator names of ``self`` than ``y.split('_')[0]``, or
        (``x.split('_')[0]==y.split('_')[0]`` and
        ``int(x.split('_')[1])<int(y.split('_')[1])``)

        EXAMPLES::

            sage: X.<alpha,beta> = InfinitePolynomialRing(ZZ)
            sage: X.varname_key('alpha_1')
            (0, 1)
            sage: X.varname_key('beta_10')
            (-1, 10)
            sage: X.varname_key('beta_1')
            (-1, 1)
            sage: X.varname_key('alpha_10')
            (0, 10)
            sage: X.varname_key('alpha_1')
            (0, 1)
            sage: X.varname_key('alpha_10')
            (0, 10)
        """
        try:
            return self._identify_variable(*x.split('_', 1))
        except (KeyError, ValueError, TypeError):
            raise ValueError("%s is not a valid variable name" % x)

    def ngens(self):
        """
        Return the number of generators for this ring.

        Since there
        are countably infinitely many variables in this polynomial
        ring, by 'generators' we mean the number of infinite families
        of variables. See :mod:`~sage.rings.polynomial.infinite_polynomial_ring`
        for more details.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(ZZ)
            sage: X.ngens()
            1

            sage: X.<x1,x2> = InfinitePolynomialRing(QQ)
            sage: X.ngens()
            2
        """
        return len(self._names)

    @cached_method
    def gen(self, i=None):
        """
        Return the `i`-th 'generator' (see the description in :meth:`.ngens`)
        of this infinite polynomial ring.

        EXAMPLES::

            sage: X = InfinitePolynomialRing(QQ)
            sage: x = X.gen()
            sage: x[1]
            x_1
            sage: X.gen() is X.gen(0)
            True
            sage: XX = InfinitePolynomialRing(GF(5))
            sage: XX.gen(0) is XX.gen()
            True
        """
        if i is not None and i > len(self._names):
            raise ValueError
        j = i if i is not None else 0
        res = InfinitePolynomialGen(self, self._names[j])
        if i is None:
            key = ((0,), ())
            if key in self._cache__gen:
                return self._cache__gen[key]
            else:
                self._cache__gen[key] = res
        return res

    def _first_ngens(self, n):
        """
        Used by the preparser for R.<x> = ...

        EXAMPLES::

            sage: InfinitePolynomialRing(ZZ, 'a')._first_ngens(1)
            (a_*,)
        """
        # It may be that we merge variables. If this is the case,
        # the new variables (as used by R.<x>  = ...) come *last*,
        # but in order.
        return self.gens()[-n:]

    @cached_method
    def gens_dict(self) -> GenDictWithBasering:
        """
        Return a dictionary-like object containing the infinitely many
        ``{var_name:variable}`` pairs.

        EXAMPLES::

            sage: R = InfinitePolynomialRing(ZZ, 'a')
            sage: D = R.gens_dict()
            sage: D
            GenDict of Infinite polynomial ring in a over Integer Ring
            sage: D['a_5']
            a_5
        """
        return GenDictWithBasering(self, InfiniteGenDict(self.gens()))

    def _ideal_class_(self, n=0):
        """
        Return :class:`SymmetricIdeals` (see there for further details).

        TESTS::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ)
            sage: R._ideal_class_()
            <class 'sage.rings.polynomial.symmetric_ideal.SymmetricIdeal'>
        """
        import sage.rings.polynomial.symmetric_ideal
        return sage.rings.polynomial.symmetric_ideal.SymmetricIdeal

    def characteristic(self):
        """
        Return the characteristic of the base field.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(GF(25,'a'))                          # needs sage.rings.finite_rings
            sage: X                                                                     # needs sage.rings.finite_rings
            Infinite polynomial ring in x, y over Finite Field in a of size 5^2
            sage: X.characteristic()                                                    # needs sage.rings.finite_rings
            5
        """
        return self._base.characteristic()

    def is_integral_domain(self, *args, **kwds):
        """
        An infinite polynomial ring is an integral domain if and only if
        the base ring is.  Arguments are passed to is_integral_domain
        method of base ring.

        EXAMPLES::

            sage: R.<x, y> = InfinitePolynomialRing(QQ)
            sage: R.is_integral_domain()
            True

        TESTS:

        :issue:`9443`::

            sage: W = PolynomialRing(InfinitePolynomialRing(QQ,'a'),2,'x,y')
            sage: W.is_integral_domain()
            True
        """
        return self._base.is_integral_domain(*args, **kwds)

    def krull_dimension(self, *args, **kwds):
        """
        Return ``Infinity``, since polynomial rings in infinitely many
        variables have infinite Krull dimension.

        EXAMPLES::

            sage: R.<x, y> = InfinitePolynomialRing(QQ)
            sage: R.krull_dimension()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    def order(self):
        """
        Return ``Infinity``, since polynomial rings have infinitely
        many elements.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(GF(2))
            sage: R.order()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    # Other bases
    def key_basis(self):
        r"""
        Return the basis of ``self`` given by key polynomials.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(GF(2))
            sage: R.key_basis()                                                         # needs sage.combinat sage.modules
            Key polynomial basis over Finite Field of size 2
        """
        from sage.combinat.key_polynomial import KeyPolynomialBasis
        return KeyPolynomialBasis(self)


class InfinitePolynomialGen(SageObject):
    """
    This class provides the object which is responsible for returning
    variables in an infinite polynomial ring (implemented in
    :meth:`.__getitem__`).

    EXAMPLES::

        sage: X.<x1,x2> = InfinitePolynomialRing(RR)
        sage: x1
        x1_*
        sage: x1[5]
        x1_5
        sage: x1 == loads(dumps(x1))
        True
    """

    def __init__(self, parent, name):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: loads(dumps(x))
            x_*
        """
        self._name = name
        self._parent = parent
        self._output = {}

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialGen
            sage: x2 = InfinitePolynomialGen(X, 'x')
            sage: x2 == x
            True
        """
        if not isinstance(other, InfinitePolynomialGen):
            return False
        return (self._name, self._parent) == (other._name, other._parent)

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialGen
            sage: x2 = InfinitePolynomialGen(X, 'x')
            sage: x2 != x
            False
        """
        return not (self == other)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.misc.latex import latex
            sage: X.<x,x1,xx> = InfinitePolynomialRing(QQ)
            sage: latex(x) # indirect doctest
            x_{\ast}
            sage: latex(x1) # indirect doctest
            \mathit{x1}_{\ast}
            sage: latex(xx) # indirect doctest
            \mathit{xx}_{\ast}
            sage: latex(x[2]) # indirect doctest
            x_{2}
            sage: latex(x1[3]) # indirect doctest
            \mathit{x1}_{3}
        """
        from sage.misc.latex import latex_variable_name
        return latex_variable_name(self._name + '_ast')

    def __getitem__(self, i):
        """
        Return the variable ``x[i]`` where ``x`` is this
        :class:`sage.rings.polynomial.infinite_polynomial_ring.InfinitePolynomialGen`,
        and i is a nonnegative integer.

        EXAMPLES::

            sage: X.<alpha> = InfinitePolynomialRing(QQ)
            sage: alpha[1]
            alpha_1
        """
        if int(i) != i:
            raise ValueError("the index (= %s) must be an integer" % i)
        i = int(i)
        if i < 0:
            raise ValueError("the index (= %s) must be nonnegative" % i)
        P = self._parent
        from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial_dense, InfinitePolynomial_sparse
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        OUT = self._output.get(i)
        if hasattr(P, '_P'):
            if i <= P._max:
                # return InfinitePolynomial_dense(P, P._P.gen(P._P.variable_names().index(self._name+'_'+str(i))))
                if OUT is None:
                    self._output[i] = InfinitePolynomial_dense(P, P._P.gen(P._P.variable_names().index(self._name + '_' + str(i))))
                else:
                    if OUT._p.parent() is not P._P:
                        OUT._p = P._P(OUT._p)
                return self._output[i]
            # Calculate all of the new names needed
            try:
                names = [[name + '_' + str(j) for name in P._names]
                         for j in range(i + 1)]
            except OverflowError:
                raise IndexError("variable index is too big - consider using the sparse implementation")
            names = reduce(operator.add, names)
            names.sort(key=P.varname_key, reverse=True)
            # Create the new polynomial ring
            P._P = PolynomialRing(P.base_ring(), names, order=P._order)
            # Get the generators
            P._max = i
            # return InfinitePolynomial_dense(P, P._P.gen(P._P.variable_names().index(self._name+'_'+str(i))))
            self._output[i] = InfinitePolynomial_dense(P, P._P.gen(P._P.variable_names().index(self._name + '_' + str(i))))
            return self._output[i]
        # Now, we are in the sparse implementation
        if OUT is not None:  # in the sparse implementation, this is ok
            return OUT
        if i == 0:
            names = [self._name + '_0', self._name + '_1']
        else:
            names = [self._name + '_0', self._name + '_' + str(i)]
        names.sort(key=P.varname_key, reverse=True)
        Pol = PolynomialRing(P.base_ring(), names, order=P._order)
        # return InfinitePolynomial_sparse(P, Pol.gen(names.index(self._name+'_'+str(i))))
        self._output[i] = InfinitePolynomial_sparse(P, Pol.gen(names.index(self._name + '_' + str(i))))
        return self._output[i]

    def _repr_(self):
        """
        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: x  # indirect doctest
            x_*
        """
        return self._name + '_*'

    def __str__(self):
        """
        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: print(x) # indirect doctest
            Generator for the x's in Infinite polynomial ring in x, y over Rational Field
        """
        return f"Generator for the {self._name}'s in {self._parent}"


##############################################################
#  The dense implementation

class InfinitePolynomialRing_dense(InfinitePolynomialRing_sparse):
    """
    Dense implementation of Infinite Polynomial Rings.

    Compared with :class:`~sage.rings.polynomial.infinite_polynomial_ring.InfinitePolynomialRing_sparse`,
    from which this class inherits, it keeps a polynomial ring that comprises all elements that have
    been created so far.
    """
    def __init__(self, R, names, order):
        """
        EXAMPLES::

            sage: X.<x2,alpha,y4> = InfinitePolynomialRing(ZZ, implementation='dense')
            sage: X == loads(dumps(X))
            True
        """
        if not names:
            names = ['x']
        # Generate the initial polynomial ring
        self._max = 0
        InfinitePolynomialRing_sparse.__init__(self, R, names, order)
        self._P = self._minP
        # self._pgens = self._P.gens()

    #####################
    # Coercion

    def construction(self):
        """
        Return the construction of ``self``.

        OUTPUT:

        A pair ``F,R``, where ``F`` is a construction functor and ``R`` is a ring,
        so that ``F(R) is self``.

        EXAMPLES::

            sage: R.<x,y> = InfinitePolynomialRing(GF(5))
            sage: R.construction()
            [InfPoly{[x,y], "lex", "dense"}, Finite Field of size 5]
        """
        return [InfinitePolynomialFunctor(self._names, self._order, 'dense'), self._base]

    def tensor_with_ring(self, R):
        """
        Return the tensor product of ``self`` with another ring.

        INPUT:

        - ``R`` -- a ring

        OUTPUT:

        An infinite polynomial ring that, mathematically, can be seen as the
        tensor product of ``self`` with ``R``.

        NOTE:

        It is required that the underlying ring of ``self`` coerces into ``R``.
        Hence, the tensor product is in fact merely an extension of the base
        ring.

        EXAMPLES::

            sage: R.<a,b> = InfinitePolynomialRing(ZZ, implementation='sparse')
            sage: R.tensor_with_ring(QQ)
            Infinite polynomial ring in a, b over Rational Field
            sage: R
            Infinite polynomial ring in a, b over Integer Ring

        The following tests against a bug that was fixed at :issue:`10468`::

            sage: R.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: R.tensor_with_ring(QQ) is R
            True
        """
        if not R.has_coerce_map_from(self._underlying_ring):
            raise TypeError("we cannot tensor with " + repr(R))
        B = self.base_ring()
        if hasattr(B, 'tensor_with_ring'):
            return InfinitePolynomialRing(B.tensor_with_ring(R), self._names, self._order, implementation='dense')
        if hasattr(B, 'change_ring'):  # e.g., polynomial rings
            return InfinitePolynomialRing(B.change_ring(R), self._names, self._order, implementation='dense')
        # try to find the correct base ring in other ways:
        try:
            o = B.one() * R.one()
        except (TypeError, ValueError):
            raise TypeError("we cannot tensor with " + repr(R))
        return InfinitePolynomialRing(o.parent(), self._names, self._order, implementation='dense')

    def polynomial_ring(self):
        """
        Return the underlying *finite* polynomial ring.

        .. NOTE::

           The ring returned can change over time as more variables
           are used.

           Since the rings are cached, we create here a ring with variable
           names that do not occur in other doc tests, so that we avoid
           side effects.

        EXAMPLES::

            sage: X.<xx, yy> = InfinitePolynomialRing(ZZ)
            sage: X.polynomial_ring()
            Multivariate Polynomial Ring in xx_0, yy_0 over Integer Ring
            sage: a = yy[3]
            sage: X.polynomial_ring()
            Multivariate Polynomial Ring in xx_3, xx_2, xx_1, xx_0, yy_3, yy_2, yy_1, yy_0
             over Integer Ring
        """
        return self._P
