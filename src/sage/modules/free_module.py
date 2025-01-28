r"""
Free modules

Sage supports computation with free modules over an arbitrary commutative ring.
Nontrivial functionality is available over `\ZZ`, fields, and some principal
ideal domains (e.g. `\QQ[x]` and rings of integers of number fields). All free
modules over an integral domain are equipped with an embedding in an ambient
vector space and an inner product, which you can specify and change.

Create the free module of rank `n` over an arbitrary commutative ring `R` using
the command ``FreeModule(R,n)``. Equivalently, ``R^n`` also creates that free
module.

The following example illustrates the creation of both a vector space and a
free module over the integers and a submodule of it.  Use the functions
``FreeModule``, ``span`` and member functions of free modules to create free
modules.  *Do not use the FreeModule_xxx constructors directly.*

EXAMPLES::

    sage: V = VectorSpace(QQ, 3)
    sage: W = V.subspace([[1,2,7], [1,1,0]])
    sage: W
    Vector space of degree 3 and dimension 2 over Rational Field
    Basis matrix:
    [ 1  0 -7]
    [ 0  1  7]
    sage: C = VectorSpaces(FiniteField(7))
    sage: C
    Category of vector spaces over Finite Field of size 7
    sage: C(W)
    Vector space of degree 3 and dimension 2 over Finite Field of size 7
    Basis matrix:
    [1 0 0]
    [0 1 0]

::

    sage: M = ZZ^3
    sage: C = VectorSpaces(FiniteField(7))
    sage: C(M)
    Vector space of dimension 3 over Finite Field of size 7
    sage: W = M.submodule([[1,2,7], [8,8,0]])
    sage: C(W)
    Vector space of degree 3 and dimension 2 over Finite Field of size 7
    Basis matrix:
    [1 0 0]
    [0 1 0]

We illustrate the exponent notation for creation of free modules.

::

    sage: ZZ^4
    Ambient free module of rank 4 over the principal ideal domain Integer Ring
    sage: QQ^2
    Vector space of dimension 2 over Rational Field
    sage: RR^3
    Vector space of dimension 3 over Real Field with 53 bits of precision

Base ring::

    sage: R.<x,y> = QQ[]
    sage: M = FreeModule(R,2)
    sage: M.base_ring()
    Multivariate Polynomial Ring in x, y over Rational Field

::

    sage: VectorSpace(QQ, 10).base_ring()
    Rational Field

Enumeration of `\ZZ^n` happens in order of increasing `1`-norm
primarily and increasing `\infty`-norm secondarily::

    sage: print([v for _,v in zip(range(31), ZZ^3)])
    [(0, 0, 0),
     (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
     (1, 1, 0), (-1, 1, 0), (1, -1, 0), (-1, -1, 0), (1, 0, 1), (-1, 0, 1), (1, 0, -1), (-1, 0, -1), (0, 1, 1), (0, -1, 1), (0, 1, -1), (0, -1, -1),
     (2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2),
     (1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1), (1, 1, -1), ...]

For other infinite enumerated base rings (i.e., rings which
are objects of the category :class:`InfiniteEnumeratedSets`),
a free module of rank `r` is enumerated by applying
:meth:`FreeModule_ambient.linear_combination_of_basis`
to all vectors in `\ZZ^r`, enumerated in the way shown above.

TESTS:

We intersect a zero-dimensional vector space with a
1-dimension submodule.

::

    sage: V = (QQ^1).span([])
    sage: W = ZZ^1
    sage: V.intersection(W)
    Free module of degree 1 and rank 0 over Integer Ring
    Echelon basis matrix:
    []

We construct subspaces of real and complex double vector spaces and
verify that the element types are correct::

    sage: V = FreeModule(RDF, 3); V
    Vector space of dimension 3 over Real Double Field
    sage: V.0
    (1.0, 0.0, 0.0)
    sage: type(V.0)
    <class 'sage.modules.vector_real_double_dense.Vector_real_double_dense'>
    sage: W = V.span([V.0]); W
    Vector space of degree 3 and dimension 1 over Real Double Field
    Basis matrix:
    [1.0 0.0 0.0]
    sage: type(W.0)
    <class 'sage.modules.vector_real_double_dense.Vector_real_double_dense'>
    sage: V = FreeModule(CDF, 3); V
    Vector space of dimension 3 over Complex Double Field
    sage: type(V.0)
    <class 'sage.modules.vector_complex_double_dense.Vector_complex_double_dense'>
    sage: W = V.span_of_basis([CDF.0 * V.1]); W
    Vector space of degree 3 and dimension 1 over Complex Double Field
    User basis matrix:
    [  0.0 1.0*I   0.0]
    sage: type(W.0)
    <class 'sage.modules.vector_complex_double_dense.Vector_complex_double_dense'>

Basis vectors are immutable::

    sage: A = span([[1,2,3], [4,5,6]], ZZ)
    sage: A.0
    (1, 2, 3)
    sage: A.0[0] = 5
    Traceback (most recent call last):
    ...
    ValueError: vector is immutable; please change a copy instead (use copy())

Among other things, this tests that we can save and load submodules
and elements::

    sage: M = ZZ^3
    sage: TestSuite(M).run()
    sage: W = M.span_of_basis([[1,2,3],[4,5,19]])
    sage: TestSuite(W).run()
    sage: v = W.0 + W.1
    sage: TestSuite(v).run()

AUTHORS:

- William Stein (2005, 2007)

- David Kohel (2007, 2008)

- Niles Johnson (2010-08): (:issue:`3893`) ``random_element()`` should pass on ``*args`` and ``**kwds``.

- Simon King (2010-12): (:issue:`8800`) fixed a bug in ``denominator()``.

- Simon King (2010-12), Peter Bruin (June 2014): (:issue:`10513`) new coercion model and category framework.
"""

###########################################################################
#       Copyright (C) 2005, 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2007, 2008 David Kohel <kohel@iml.univ-mrs.fr>
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
###########################################################################

import itertools
from warnings import warn

import sage.matrix.matrix_space
import sage.misc.latex as latex
import sage.rings.abc
import sage.rings.infinity
import sage.rings.integer
import sage.rings.integer_ring
import sage.rings.rational_field
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.integral_domains import IntegralDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.misc.randstate import current_randstate
from sage.modules import free_module_element
from sage.modules.module import Module
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.structure.factory import UniqueFactory
from sage.structure.richcmp import (
    op_EQ,
    op_GE,
    op_GT,
    op_LE,
    op_LT,
    op_NE,
    revop,
    rich_to_bool,
    richcmp,
    richcmp_method,
    richcmp_not_equal,
)
from sage.structure.sequence import Sequence


###############################################################################
#
# Constructor functions
#
###############################################################################

class FreeModuleFactory(UniqueFactory):
    r"""
    Factory class for the finite-dimensional free modules with standard basis
    """
    def create_key(self, base_ring, rank, sparse=False, inner_product_matrix=None):
        """
        TESTS::

            sage: loads(dumps(ZZ^6)) is ZZ^6
            True
            sage: loads(dumps(RDF^3)) is RDF^3
            True
        """
        rank = int(sage.rings.integer.Integer(rank))

        if inner_product_matrix is not None:
            inner_product_matrix = sage.matrix.matrix_space.MatrixSpace(base_ring, rank)(inner_product_matrix)
            inner_product_matrix.set_immutable()

        return (base_ring, rank, sparse, inner_product_matrix)

    def create_object(self, version, key):
        """
        TESTS::

            sage: TestSuite(ZZ^6).run()
            sage: TestSuite(RDF^3).run()

        Check that :issue:`34380` is fixed::

            sage: R.<x,y> = QQ[]
            sage: Q = R.quo(R.ideal([x^2 - y^2 - 1]))
            sage: Q.is_integral_domain()                                                # needs sage.libs.singular
            True
            sage: Q2 = FreeModule(Q, 2)                                                 # needs sage.libs.singular
            sage: from sage.modules.free_module import FreeModule_ambient_domain
            sage: isinstance(Q2, FreeModule_ambient_domain)                             # needs sage.libs.singular
            True
        """
        base_ring, rank, sparse, inner_product_matrix = key

        if inner_product_matrix is not None:
            from sage.modules.free_quadratic_module import FreeQuadraticModule
            return FreeQuadraticModule(base_ring, rank, inner_product_matrix=inner_product_matrix, sparse=sparse)

        if not isinstance(sparse, bool):
            raise TypeError("Argument sparse (= %s) must be True or False" % sparse)

        if base_ring not in CommutativeRings():
            warn("You are constructing a free module\n"
                 "over a noncommutative ring. Sage does not have a concept\n"
                 "of left/right and both sided modules, so be careful.\n"
                 "It's also not guaranteed that all multiplications are\n"
                 "done from the right side.")
            # raise TypeError("the base_ring must be a commutative ring")

        if not sparse and isinstance(base_ring, sage.rings.abc.RealDoubleField):
            return RealDoubleVectorSpace_class(rank)

        if not sparse and isinstance(base_ring, sage.rings.abc.ComplexDoubleField):
            return ComplexDoubleVectorSpace_class(rank)

        try:
            if base_ring.is_field():
                return FreeModule_ambient_field(base_ring, rank, sparse=sparse)
        except NotImplementedError:
            pass

        if base_ring in PrincipalIdealDomains():
            return FreeModule_ambient_pid(base_ring, rank, sparse=sparse)

        if (isinstance(base_ring, sage.rings.abc.Order)
                and base_ring.is_maximal() and base_ring.class_number() == 1):
            return FreeModule_ambient_pid(base_ring, rank, sparse=sparse)

        if base_ring in IntegralDomains():
            return FreeModule_ambient_domain(base_ring, rank, sparse=sparse)

        return FreeModule_ambient(base_ring, rank, sparse=sparse)


FreeModuleFactory_with_standard_basis = FreeModuleFactory("FreeModule")


def FreeModule(base_ring, rank_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
               with_basis='standard', rank=None, basis_keys=None, **args):
    r"""
    Create a free module over the given commutative ``base_ring``.

    ``FreeModule`` can be called with the following positional arguments:

    - ``FreeModule(base_ring, rank, ...)``

    - ``FreeModule(base_ring, basis_keys, ...)``

    INPUT:

    - ``base_ring`` -- a commutative ring

    - ``rank`` -- nonnegative integer

    - ``basis_keys`` -- a finite or enumerated family of arbitrary objects

    - ``sparse`` -- boolean (default: ``False``)

    - ``inner_product_matrix`` -- the inner product matrix (default: ``None``)

    - ``with_basis`` -- either ``'standard'`` (the default), in which case
      a free module with the standard basis as the distinguished basis is created;
      or ``None``, in which case a free module without distinguished basis is
      created.

    - further options may be accepted by various implementation classes

    OUTPUT: a free module

    This factory function creates instances of various specialized classes
    depending on the input.  Not all combinations of options are
    implemented.

    - If the parameter ``basis_keys`` is provided, it must be a finite
      or enumerated family of objects, and an instance of
      :class:`CombinatorialFreeModule` is created.

       EXAMPLES::

           sage: CombinatorialFreeModule(QQ, ['a','b','c'])
           Free module generated by {'a', 'b', 'c'} over Rational Field

       It has a distinguished standard basis that is indexed by the provided
       ``basis_keys``. See the documentation of :class:`CombinatorialFreeModule`
       for more examples and details, including its :class:`UniqueRepresentation`
       semantics.

    - If the parameter ``with_basis`` is set to ``None``, then a free module
      of the given ``rank`` without distinguished basis is created.  It is
      represented by an instance of :class:`FiniteRankFreeModule`.

       EXAMPLES::

           sage: FiniteRankFreeModule(ZZ, 3, name='M')
           Rank-3 free module M over the Integer Ring

       See the documentation of :class:`FiniteRankFreeModule` for more
       options, examples, and details.

    - If ``rank`` is provided and the option ``with_basis`` is left at its
      default value, ``'standard'``, then a free ambient module with
      distinguished standard basis indexed by ``range(rank)`` is created.
      There is only one dense and one sparse free ambient module of
      given ``rank`` over ``base_ring``.

      EXAMPLES::

          sage: FreeModule(Integers(8), 10)
          Ambient free module of rank 10 over Ring of integers modulo 8

      The remainder of this documentation discusses this case of
      free ambient modules.

    EXAMPLES:

    First we illustrate creating free modules over various base fields.
    The base field affects the free module that is created. For
    example, free modules over a field are vector spaces, and free
    modules over a principal ideal domain are special in that more
    functionality is available for them than for completely general
    free modules.

    ::

        sage: FreeModule(QQ,10)
        Vector space of dimension 10 over Rational Field
        sage: FreeModule(ZZ,10)
        Ambient free module of rank 10 over the principal ideal domain Integer Ring
        sage: FreeModule(FiniteField(5), 10)
        Vector space of dimension 10 over Finite Field of size 5
        sage: FreeModule(Integers(7), 10)
        Vector space of dimension 10 over Ring of integers modulo 7
        sage: FreeModule(PolynomialRing(QQ,'x'), 5)
        Ambient free module of rank 5 over the principal ideal domain
         Univariate Polynomial Ring in x over Rational Field
        sage: FreeModule(PolynomialRing(ZZ,'x'), 5)
        Ambient free module of rank 5 over the integral domain
         Univariate Polynomial Ring in x over Integer Ring

    Of course we can make rank 0 free modules::

        sage: FreeModule(RealField(100),0)
        Vector space of dimension 0 over Real Field with 100 bits of precision

    Next we create a free module with sparse representation of
    elements. Functionality with sparse modules is *identical* to dense
    modules, but they may use less memory and arithmetic may be faster
    (or slower!).

    ::

        sage: M = FreeModule(ZZ,200,sparse=True)
        sage: M.is_sparse()
        True
        sage: type(M.0)
        <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>

    The default is dense.

    ::

        sage: M = ZZ^200
        sage: type(M.0)
        <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

    Note that matrices associated in some way to sparse free modules
    are sparse by default::

        sage: M = FreeModule(Integers(8), 2)
        sage: A = M.basis_matrix()
        sage: A.is_sparse()
        False
        sage: Ms = FreeModule(Integers(8), 2, sparse=True)
        sage: M == Ms  # as mathematical objects they are equal
        True
        sage: Ms.basis_matrix().is_sparse()
        True

    We can also specify an inner product matrix, which is used when
    computing inner products of elements.

    ::

        sage: A = MatrixSpace(ZZ, 2)([[1,0], [0,-1]])
        sage: M = FreeModule(ZZ, 2, inner_product_matrix=A)
        sage: v, w = M.gens()
        sage: v.inner_product(w)
        0
        sage: v.inner_product(v)
        1
        sage: w.inner_product(w)
        -1
        sage: (v+2*w).inner_product(w)
        -2

    You can also specify the inner product matrix by giving anything
    that coerces to an appropriate matrix. This is only useful if the
    inner product matrix takes values in the base ring.

    ::

        sage: FreeModule(ZZ, 2, inner_product_matrix=1).inner_product_matrix()
        [1 0]
        [0 1]
        sage: FreeModule(ZZ, 2, inner_product_matrix=[1,2,3,4]).inner_product_matrix()
        [1 2]
        [3 4]
        sage: FreeModule(ZZ, 2, inner_product_matrix=[[1,2], [3,4]]).inner_product_matrix()
        [1 2]
        [3 4]

    .. TODO::

        Refactor modules such that it only counts what category the base
        ring belongs to, but not what is its Python class.

    EXAMPLES::

        sage: FreeModule(QQ, ['a', 'b', 'c'])
        Free module generated by {'a', 'b', 'c'} over Rational Field
        sage: _.category()
        Category of finite dimensional vector spaces with basis over Rational Field

        sage: FreeModule(QQ, 3, with_basis=None)
        3-dimensional vector space over the Rational Field
        sage: _.category()
        Category of finite dimensional vector spaces over Rational Field

        sage: FreeModule(QQ, [1, 2, 3, 4], with_basis=None)
        4-dimensional vector space over the Rational Field
        sage: _.category()
        Category of finite dimensional vector spaces over Rational Field

    TESTS::

        sage: FreeModule(QQ, ['a', 2, 3, 4], with_basis=None)
        Traceback (most recent call last):
        ...
        NotImplementedError: FiniteRankFreeModule only supports integer ranges
        as basis_keys, got ['a', 2, 3, 4]
        sage: FreeModule(QQ, [1, 3, 5], with_basis=None)
        Traceback (most recent call last):
        ...
        NotImplementedError: FiniteRankFreeModule only supports integer ranges
        as basis_keys, got [1, 3, 5]

        sage: FreeModule(ZZ, rank=3, basis_keys=['c','d'])
        Traceback (most recent call last):
        ...
        ValueError: inconsistent rank: should be cardinality of ['c', 'd'] but got 3
        sage: FreeModule(QQ, ['a', 'b'], rank=2)
        Free module generated by {'a', 'b'} over Rational Field
        sage: FreeModule(QQ, 2, basis_keys=['a', 'b'])
        Free module generated by {'a', 'b'} over Rational Field
        sage: FreeModule(QQ, ['a', 'b'], rank=3)
        Traceback (most recent call last):
        ...
        ValueError: inconsistent rank: should be cardinality of ['a', 'b'] but got 3
    """
    if rank_or_basis_keys is not None:
        try:
            n = sage.rings.integer_ring.ZZ(rank_or_basis_keys)
        except (TypeError, ValueError):
            if basis_keys is not None:
                raise ValueError("duplicate values for basis_keys")
            basis_keys = rank_or_basis_keys
        else:
            if rank is not None:
                raise ValueError("duplicate values for rank")
            rank = n

    if rank is not None and basis_keys is not None and rank != len(basis_keys):
        raise ValueError(f"inconsistent rank: should be cardinality of {basis_keys} "
                         f"but got {rank}")

    if not with_basis:
        if inner_product_matrix is not None:
            raise NotImplementedError
        from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
        if basis_keys:
            if not all(key in sage.rings.integer_ring.ZZ for key in basis_keys):
                raise NotImplementedError(f'FiniteRankFreeModule only supports integer ranges as basis_keys, got {basis_keys}')
            start_index = min(basis_keys)
            end_index = max(basis_keys)
            rank = end_index - start_index + 1
            # Check that the ordered list of basis_keys is the range from start_index to end_index
            if (len(basis_keys) != rank
                or not all(key == index
                           for key, index in zip(basis_keys,
                                                 range(start_index, end_index + 1)))):
                raise NotImplementedError(f'FiniteRankFreeModule only supports integer ranges as basis_keys, got {basis_keys}')
            return FiniteRankFreeModule(base_ring, rank, start_index=start_index, **args)
        return FiniteRankFreeModule(base_ring, rank, **args)
    elif with_basis == 'standard':
        if rank is not None and basis_keys is None:
            return FreeModuleFactory_with_standard_basis(base_ring, rank, sparse,
                                                         inner_product_matrix, **args)
        else:
            if inner_product_matrix is not None:
                raise NotImplementedError
            if rank is not None and rank != len(basis_keys):
                raise ValueError(f'inconsistent basis_keys: should be of cardinality {rank}, '
                                 f'got {basis_keys}')
            from sage.combinat.free_module import CombinatorialFreeModule
            return CombinatorialFreeModule(base_ring, basis_keys, **args)
    else:
        raise NotImplementedError


def VectorSpace(K, dimension_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
                with_basis='standard', dimension=None, basis_keys=None, **args):
    """
    EXAMPLES:

    The base can be complicated, as long as it is a field.

    ::

        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x over Integer Ring
        sage: V.basis()
        [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

    The base must be a field or a :exc:`TypeError` is raised.

    ::

        sage: VectorSpace(ZZ,5)
        Traceback (most recent call last):
        ...
        TypeError: Argument K (= Integer Ring) must be a field.
    """
    if not K.is_field():
        raise TypeError("Argument K (= %s) must be a field." % K)
    if sparse not in (True, False):
        raise TypeError("Argument sparse (= %s) must be a boolean." % sparse)
    return FreeModule(K, dimension_or_basis_keys, sparse, inner_product_matrix,
                      with_basis=with_basis, rank=dimension, basis_keys=basis_keys,
                      **args)


def span(gens, base_ring=None, check=True, already_echelonized=False):
    r"""
    Return the span of the vectors in ``gens`` using scalars from ``base_ring``.

    INPUT:

    - ``gens`` -- list of either vectors or lists of ring elements
      used to generate the span

    - ``base_ring`` -- (default: ``None``) a principal ideal domain
      for the ring of scalars

    - ``check`` -- (default: ``True``) passed to the ``span()`` method
      of the ambient module

    - ``already_echelonized`` -- (default: ``False``) set to ``True``
      if the vectors form the rows of a matrix in echelon form, in
      order to skip the computation of an echelonized basis for the
      span.

    OUTPUT:

    A module (or vector space) that is all the linear combinations of the
    free module elements (or vectors) with scalars from the
    ring (or field) given by ``base_ring``.  See the examples below
    describing behavior when the base ring is not specified and/or
    the module elements are given as lists that do not carry
    explicit base ring information.

    EXAMPLES:

    The vectors in the list of generators can be given as
    lists, provided a base ring is specified and the elements of the list
    are in the ring (or the fraction field of the ring).  If the
    base ring is a field, the span is a vector space.  ::

        sage: V = span([[1,2,5], [2,2,2]], QQ); V
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -3]
        [ 0  1  4]

        sage: span([V.gen(0)], QuadraticField(-7,'a'))                                  # needs sage.rings.number_field
        Vector space of degree 3 and dimension 1 over Number Field in a
         with defining polynomial x^2 + 7 with a = 2.645751311064591?*I
        Basis matrix:
        [ 1  0 -3]

        sage: span([[1,2,3], [2,2,2], [1,2,5]], GF(2))
        Vector space of degree 3 and dimension 1 over Finite Field of size 2
        Basis matrix:
        [1 0 1]

    If the base ring is not a field, then a module is created.
    The entries of the vectors can lie outside the ring, if they
    are in the fraction field of the ring.  ::

        sage: span([[1,2,5], [2,2,2]], ZZ)
        Free module of degree 3 and rank 2 over Integer Ring
        Echelon basis matrix:
        [ 1  0 -3]
        [ 0  2  8]

        sage: span([[1,1,1], [1,1/2,1]], ZZ)
        Free module of degree 3 and rank 2 over Integer Ring
        Echelon basis matrix:
        [  1   0   1]
        [  0 1/2   0]

        sage: R.<x> = QQ[]
        sage: M= span( [[x, x^2+1], [1/x, x^3]], R); M
        Free module of degree 2 and rank 2 over
        Univariate Polynomial Ring in x over Rational Field
        Echelon basis matrix:
        [          1/x           x^3]
        [            0 x^5 - x^2 - 1]
        sage: M.basis()[0][0].parent()
        Fraction Field of Univariate Polynomial Ring in x over Rational Field

    A base ring can be inferred if the generators are given as a
    list of vectors. ::

        sage: span([vector(QQ, [1,2,3]), vector(QQ, [4,5,6])])
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -1]
        [ 0  1  2]
        sage: span([vector(QQ, [1,2,3]), vector(ZZ, [4,5,6])])
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -1]
        [ 0  1  2]
        sage: span([vector(ZZ, [1,2,3]), vector(ZZ, [4,5,6])])
        Free module of degree 3 and rank 2 over Integer Ring
        Echelon basis matrix:
        [1 2 3]
        [0 3 6]

    TESTS::

        sage: span([[1,2,3], [2,2,2], [1,2/3,5]], ZZ)
        Free module of degree 3 and rank 3 over Integer Ring
        Echelon basis matrix:
        [  1   0  13]
        [  0 2/3   6]
        [  0   0  14]
        sage: span([[1,2,3], [2,2,2], [1,2,QQ['x'].gen()]], ZZ)
        Traceback (most recent call last):
        ...
        ValueError: The elements of gens (= [[1, 2, 3], [2, 2, 2], [1, 2, x]]) must be defined over base_ring (= Integer Ring) or its field of fractions.

    For backwards compatibility one can also give the base ring as the
    first argument.  ::

        sage: span(QQ,[[1,2],[3,4]])
        Vector space of degree 2 and dimension 2 over Rational Field
        Basis matrix:
        [1 0]
        [0 1]

    The base ring must be a principal ideal domain (PID).  ::

        sage: span([[1,2,3]], Integers(6))
        Traceback (most recent call last):
        ...
        TypeError: The base_ring (= Ring of integers modulo 6)
        must be a principal ideal domain.

    Fix :issue:`5575`::

        sage: V = QQ^3
        sage: span([V.0, V.1])
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [1 0 0]
        [0 1 0]

    Improve error message from :issue:`12541`::

        sage: span({0:vector([0,1])}, QQ)
        Traceback (most recent call last):
        ...
        TypeError: generators must be lists of ring elements
        or free module elements!
    """
    if gens in CommutativeRings():
        # we allow the old input format with first input the base_ring.
        # Do we want to deprecate it?..
        base_ring, gens = gens, base_ring

    try:
        if base_ring is None:
            gens = Sequence(gens)
            R = gens.universe().base_ring()
        else:
            gens = list(gens)
            R = base_ring
    except TypeError:
        raise TypeError("generators must be given as an iterable structure")

    if R not in PrincipalIdealDomains():
        raise TypeError("The base_ring (= %s) must be a principal ideal "
                        "domain." % R)
    if not gens:
        return FreeModule(R, 0)
    else:
        x = gens[0]
        if isinstance(x, free_module_element.FreeModuleElement):
            M = x.parent()
        else:
            try:
                x = list(x)
            except TypeError:
                raise TypeError("generators must be lists of ring elements or "
                                "free module elements!")
            M = FreeModule(R, len(x))
            try:
                gens = [M(_) for _ in gens]
            except TypeError:
                R = R.fraction_field()
                M = FreeModule(R, len(x))
                try:
                    gens = [M(_) for _ in gens]
                except TypeError:
                    raise ValueError("The elements of gens (= %s) must be "
                                     "defined over base_ring (= %s) or its "
                                     "field of fractions." % (gens, base_ring))
        return M.span(gens=gens, base_ring=base_ring, check=check,
                      already_echelonized=already_echelonized)


def basis_seq(V, vecs):
    """
    This converts a list vecs of vectors in V to a Sequence of
    immutable vectors.

    Should it? I.e. in most ``other`` parts of the system the return type
    of basis or generators is a tuple.

    EXAMPLES::

        sage: V = VectorSpace(QQ,2)
        sage: B = V.gens()
        sage: B
        ((1, 0), (0, 1))
        sage: v = B[0]
        sage: v[0] = 0 # immutable
        Traceback (most recent call last):
        ...
        ValueError: vector is immutable; please change a copy instead (use copy())
        sage: sage.modules.free_module.basis_seq(V, V.gens())
        [(1, 0), (0, 1)]
    """
    for z in vecs:
        z.set_immutable()
    return Sequence(vecs, universe=V, check=False, immutable=True, cr=True)


###############################################################################
#
# Base class for all free modules
#
###############################################################################

def is_FreeModule(M):
    """
    Return ``True`` if M inherits from ``FreeModule_generic``.

    EXAMPLES::

        sage: from sage.modules.free_module import is_FreeModule
        sage: V = ZZ^3
        sage: is_FreeModule(V)
        doctest:warning...
        DeprecationWarning: the function is_FreeModule is deprecated;
        use 'isinstance(..., FreeModule_generic)' instead
        See https://github.com/sagemath/sage/issues/37924 for details.
        True
        sage: W = V.span([ V.random_element() for i in range(2) ])
        sage: is_FreeModule(W)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37924, "the function is_FreeModule is deprecated; use 'isinstance(..., FreeModule_generic)' instead")
    return isinstance(M, FreeModule_generic)


@richcmp_method
class Module_free_ambient(Module):
    """
    Base class for modules with elements represented by elements of a free
    module.

    Modules whose elements are represented by elements of a free module (such
    as submodules, quotients, and subquotients of a free module) should be
    either a subclass of this class or :class:`FreeModule_generic`, which
    itself is a subclass of this class. If the modules have bases and ranks,
    then use :class:`FreeModule_generic`. Otherwise, use this class.

    INPUT:

    - ``base_ring`` -- a commutative ring

    - ``degree`` -- nonnegative integer; degree of the ambient free module

    - ``sparse`` -- boolean (default: ``False``)

    - ``category`` -- category (default: ``None``)

    If ``base_ring`` is a field, then the default category is the category of
    finite-dimensional vector spaces over that field; otherwise it is the
    category of finite-dimensional free modules over that ring.  In addition,
    the category is intersected with the category of finite enumerated sets if
    the ring is finite or the rank is 0.

    EXAMPLES::

        sage: S.<x,y,z> = PolynomialRing(QQ)
        sage: M = S**2
        sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
        sage: N.gens()
        [(x - y, z), (y*z, x*z)]
        sage: N.degree()
        2
    """
    def __init__(self, base_ring, degree, sparse=False, category=None):
        """
        Initialize.

        TESTS::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
            sage: TestSuite(N).run(skip=['_test_elements', '_test_pickling'])
        """
        degree = sage.rings.integer.Integer(degree)
        if degree < 0:
            raise ValueError("degree (=%s) must be nonnegative" % degree)

        from sage.categories.modules_with_basis import ModulesWithBasis
        modules_category = ModulesWithBasis(base_ring.category()).FiniteDimensional()
        try:
            if base_ring.is_finite() or degree == 0:
                modules_category = modules_category.Enumerated().Finite()
        except (ValueError, TypeError, AttributeError, NotImplementedError):
            pass
        category = modules_category.or_subcategory(category, join=True)

        if not hasattr(self, 'Element'):
            self.Element = element_class(base_ring, sparse)

        super().__init__(base_ring, category=category)
        self.__degree = degree
        self.__is_sparse = sparse

    def _element_constructor_(self, x, coerce=True, copy=True, check=True):
        r"""
        Create an element of this module from ``x``.

        The ``coerce`` and ``copy`` arguments are passed on to the underlying
        element constructor.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y*z , x*z])])
            sage: Q = M.quotient_module(N)
            sage: Q(0)
            (0, 0)
            sage: Q([x, x + y])
            (x, x + y)
            sage: phi = Q.coerce_map_from(M)
            sage: phi(M.gen(1))
            (0, 1)
            sage: _ in Q
            True
        """
        if isinstance(x, (int, sage.rings.integer.Integer)) and x == 0:
            return self.zero_vector()
        elif isinstance(x, free_module_element.FreeModuleElement):
            if x.parent() is self:
                if copy:
                    return x.__copy__()
                else:
                    return x
            x = x.list()
        if check and self.coordinate_ring().is_exact():
            # No check if x belongs to this module as there is no algorithm.
            try:
                R = self.base_ring()
                for d in x:
                    if d not in R:
                        raise ArithmeticError
            except ArithmeticError:
                raise TypeError("element {!r} is not in free module".format(x))
        return self.element_class(self, x, coerce, copy)

    def degree(self):
        """
        Return the degree of this free module. This is the dimension of the
        ambient vector space in which it is embedded.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 10)
            sage: W = M.submodule([M.gen(0), 2*M.gen(3) - M.gen(0), M.gen(0) + M.gen(3)])
            sage: W.degree()
            10
            sage: W.rank()
            2
        """
        return self.__degree

    def is_sparse(self):
        """
        Return ``True`` if the underlying representation of this module uses
        sparse vectors, and ``False`` otherwise.

        EXAMPLES::

            sage: FreeModule(ZZ, 2).is_sparse()
            False
            sage: FreeModule(ZZ, 2, sparse=True).is_sparse()
            True
        """
        return self.__is_sparse

    def _an_element_(self):
        """
        Return an arbitrary element of a free module.

        EXAMPLES::

            sage: V = VectorSpace(QQ,2)
            sage: V._an_element_()
            (1, 0)
            sage: U = V.submodule([[1,0]])
            sage: U._an_element_()
            (1, 0)
            sage: W = V.submodule([])
            sage: W._an_element_()
            (0, 0)
        """
        try:
            return self.gen(0)
        except ValueError:
            return self.zero()

    def some_elements(self):
        r"""
        Return some elements of this free module.

        See :class:`TestSuite` for a typical use case.

        OUTPUT: an iterator

        EXAMPLES::

            sage: F = FreeModule(ZZ, 2)
            sage: tuple(F.some_elements())
            ((1, 0),
             (1, 1),
             (0, 1),
             (-1, 2),
             (-2, 3),
             ...
             (-49, 50))

            sage: F = FreeModule(QQ, 3)
            sage: tuple(F.some_elements())
            ((1, 0, 0),
             (1/2, 1/2, 1/2),
             (1/2, -1/2, 2),
             (-2, 0, 1),
             (-1, 42, 2/3),
             (-2/3, 3/2, -3/2),
             (4/5, -4/5, 5/4),
             ...
             (46/103823, -46/103823, 103823/46))

            sage: F = FreeModule(SR, 2)                                                 # needs sage.symbolic
            sage: tuple(F.some_elements())                                              # needs sage.symbolic
            ((1, 0), (some_variable, some_variable))
        """
        yield self.an_element()
        yield self.base().an_element() * sum(self.gens())
        some_elements_base = iter(self.base().some_elements())
        n = self.degree()
        while True:
            L = list(itertools.islice(some_elements_base, int(n)))
            if len(L) != n:
                return
            try:
                yield self(L)
            except (TypeError, ValueError):
                pass

    def coordinate_ring(self):
        """
        Return the ring over which the entries of the vectors are
        defined.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
            sage: N.coordinate_ring()
            Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.base_ring()

    def zero_vector(self):
        """
        Return the zero vector in this module.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 2)
            sage: M.zero_vector()
            (0, 0)
            sage: M(0)
            (0, 0)
            sage: M.span([[1,1]]).zero_vector()
            (0, 0)
            sage: M.zero_submodule().zero_vector()
            (0, 0)
        """
        # Do *not* cache this -- it must be computed fresh each time, since
        # it is used by __call__ to make a new copy of the 0 element.

        return self.element_class(self, 0)

    @cached_method
    def zero(self):
        """
        Return the zero vector in this module.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 2)
            sage: M.zero()
            (0, 0)
            sage: M.span([[1,1]]).zero()
            (0, 0)
            sage: M.zero_submodule().zero()
            (0, 0)
            sage: M.zero_submodule().zero().is_mutable()
            False
        """
        res = self.element_class(self, 0)
        res.set_immutable()
        return res

    def zero_submodule(self):
        """
        Return the zero submodule of this module.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: M.zero_submodule()
            Submodule of Ambient free module of rank 2 over the integral domain Multivariate Polynomial Ring in x, y, z over Rational Field
            Generated by the rows of the matrix:
            []
        """
        return self.submodule([], check=False)

    def relations_matrix(self):
        r"""
        Return the matrix of relations of ``self``.

        EXAMPLES::

            sage: V = GF(2)^2
            sage: V.relations_matrix()
            []
            sage: W = V.subspace([[1, 0]])
            sage: W.relations_matrix()
            []

            sage: Q = V / W
            sage: Q.relations_matrix()
            [1 0]

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: M.relations_matrix()
            []

            sage: N = M.submodule([vector([x - y, z]), vector([y*z, x*z])])
            sage: Q = M.quotient_module(N)
            sage: Q.relations_matrix()
            [x - y     z]
            [  y*z   x*z]
        """
        return self.relations().matrix()

    def __richcmp__(self, other, op):
        r"""
        Rich comparison via containment in the same ambient space.

        Two modules compare if their ambient module/space is equal.
        Ambient spaces are equal if they have the same
        base ring, rank and inner product matrix.

        EXAMPLES:

        We compare rank three free modules over the integers,
        rationals, and complex numbers. Note the free modules
        ``QQ^3`` (and hence ``ZZ^3``) and ``CC^3`` are incomparable
        because of the different ambient vector spaces::

            sage: QQ^3 <= CC^3
            False
            sage: CC^3 <= QQ^3
            False
            sage: QQ^3 <= QQ^3
            True

            sage: QQ^3 <= ZZ^3
            False
            sage: ZZ^3 <= QQ^3
            True
            sage: ZZ^3 <= CC^3
            False
            sage: CC^3 <= ZZ^3
            False

        Comparison with a submodule::

            sage: A = QQ^3
            sage: V = span([[1,2,3], [5,6,7], [8,9,10]], QQ)
            sage: V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: V <= A
            True
            sage: A <= V
            False
            sage: L1 = span([[1,2,3], [5,6,7], [8,9,10]], ZZ)
            sage: L2 = span([[2,4,6], [10,12,14], [16,18,20]], ZZ)
            sage: 2*L1 <= L2
            True
            sage: L2 <= L1
            True
            sage: L1 <= L2
            False

        More exotic comparisons::

            sage: # needs sage.symbolic
            sage: R1 = ZZ[sqrt(2)]
            sage: F1 = R1^3
            sage: V1 = F1.span([[sqrt(2), sqrt(2), 0]])
            sage: F2 = ZZ^3
            sage: V2 = F2.span([[2,2,0]])
            sage: V2 <= V1  # Different ambient vector spaces
            False
            sage: V1 <= V2
            False

            sage: R2.<x> = GF(5)[]
            sage: F3 = R2^3
            sage: V3 = F3.span([[x^5 - 1, 1 + x + x^2 + x^3 + x^4, 0]])
            sage: W3 = F3.span([[1,1,0], [0,4,0]])
            sage: V3 <= W3
            True
            sage: W3 <= V3
            False

        We compare a one dimensional space to a two dimensional space::

            sage: V = span([[1,2,3], [5,6,7], [8,9,10]], QQ)
            sage: M = span([[5,6,7]], QQ)
            sage: M <= V
            True
            sage: V <= M
            False

        We test that :issue:`5525` is fixed::

            sage: A = (QQ^1).span([[1/3]], ZZ); B = (QQ^1).span([[1]], ZZ)
            sage: A.intersection(B)
            Free module of degree 1 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1]

        We create the module `\ZZ^3`, and the submodule generated by
        one vector `(1,1,0)`, and check whether certain elements are
        in the submodule::

            sage: R = FreeModule(ZZ, 3)
            sage: V = R.submodule([R.gen(0) + R.gen(1)])
            sage: R.gen(0) + R.gen(1) in V
            True
            sage: R.gen(0) + 2*R.gen(1) in V
            False

            sage: w = (1/2)*(R.gen(0) + R.gen(1))
            sage: w
            (1/2, 1/2, 0)
            sage: w.parent()
            Vector space of dimension 3 over Rational Field
            sage: w in V
            False
            sage: V.coordinates(w)
            [1/2]

        We check an example over a more general ring::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y*z, x*z])])
            sage: Q = M.quotient_module(N)
            sage: Q == Q
            True
            sage: Q == M
            False
            sage: Q == N
            False
            sage: Q != M
            True
            sage: Q != N
            True

        When equality cannot be checked, we get both ``==`` and ``!=``
        returning ``False``::

            sage: P.<x,y> = QQ[]
            sage: M = P**2
            sage: S1 = M.submodule([(x, y), (y, x)])
            sage: S2 = M.submodule([(x, y), (y - x, x - y)])
            sage: S1 == S2
            False
            sage: S1 != S2
            False

        TESTS::

            sage: QQ^3 < ZZ^3
            False
            sage: ZZ^3 < QQ^3
            True
            sage: ZZ^3 < CC^3
            False
            sage: CC^3 < ZZ^3
            False

        Comparison with a submodule::

            sage: A = QQ^3
            sage: V = span([[1,2,3], [5,6,7], [8,9,10]], QQ)
            sage: V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: V < A
            True
            sage: A < V
            False
            sage: L1 = span([[1,2,3], [5,6,7], [8,9,10]], ZZ)
            sage: L2 = span([[2,4,6], [10,12,14], [16,18,20]], ZZ)
            sage: 2*L1 < L2
            False
            sage: L2 < L1
            True
            sage: L1 < L2
            False

        We compare a `\ZZ`-module to a one-dimensional space::

            sage: V = span([[5,6,7]], ZZ).scale(1/11);  V
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [5/11 6/11 7/11]
            sage: M = span([[5,6,7]], QQ)
            sage: V < M
            True
            sage: M < V
            False

        We compare rank three free modules over the rationals and
        complex numbers::

            sage: QQ^3 >= CC^3
            False
            sage: CC^3 >= QQ^3
            False
            sage: QQ^3 >= QQ^3
            True

        Comparison with a submodule::

            sage: A = QQ^3
            sage: V = span([[1,2,3], [5,6,7], [8,9,10]], QQ)
            sage: V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: V >= A
            False
            sage: A >= V
            True
            sage: L1 = span([[1,2,3], [5,6,7], [8,9,10]], ZZ)
            sage: L2 = span([[2,4,6], [10,12,14], [16,18,20]], ZZ)
            sage: 2*L1 >= L2
            True
            sage: L2 >= L1
            False
            sage: L1 >= L2
            True

        We compare rank three free modules over the rationals and
        complex numbers::

            sage: QQ^3 > CC^3
            False
            sage: CC^3 > QQ^3
            False
            sage: QQ^3 > QQ^3
            False

        Comparison with a submodule::

            sage: A = QQ^3
            sage: V = span([[1,2,3], [5,6,7], [8,9,10]], QQ)
            sage: V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: V > A
            False
            sage: A > V
            True
            sage: L1 = span([[1,2,3], [5,6,7], [8,9,10]], ZZ)
            sage: L2 = span([[2,4,6], [10,12,14], [16,18,20]], ZZ)
            sage: 2*L1 > L2
            False
            sage: L2 > L1
            False
            sage: L1 > L2
            True
        """
        if self is other:
            return rich_to_bool(op, 0)
        if not isinstance(other, Module_free_ambient):
            return NotImplemented

        # Check equality first if needed
        try:
            if op == op_EQ:
                return self._eq(other)
            if op == op_NE:
                return not self._eq(other)
        except NotImplementedError:
            return False

        try:
            if op == op_LE:
                return self.is_submodule(other)
            if op == op_GE:
                return other.is_submodule(self)
            if op == op_LT:
                return (not self._eq(other)) and self.is_submodule(other)
            if op == op_GT:
                return (not self._eq(other)) and other.is_submodule(self)
        except NotImplementedError:
            return NotImplemented

    def _eq(self, other):
        r"""
        Return if ``self`` is equal to ``other``.

        Ambient spaces are considered equal if they have the same
        rank, basering and inner product matrix.

        Modules in the same ambient space are partially ordered by inclusion.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y*z, x*z])])
            sage: Q = M.quotient_module(N)
            sage: M._eq(Q)
            False
            sage: M.zero_submodule()._eq(Q.zero_submodule())
            False
            sage: Q.zero_submodule()._eq(M.zero_submodule())
            False
            sage: M.zero_submodule()._eq(N.zero_submodule())
            True
        """
        if self.degree() != other.degree():
            return False
        if self.base_ring() != other.base_ring():
            return False

        from sage.modules.quotient_module import QuotientModule_free_ambient
        lq = isinstance(self, QuotientModule_free_ambient)
        rq = isinstance(other, QuotientModule_free_ambient)
        if lq or rq:
            # if the relations agree we continue with the covers
            if lq:
                lx = self.relations()
                self = self.cover()
            else:
                lx = self.zero_submodule()
            if rq:
                rx = other.relations()
                other = other.cover()
            else:
                rx = other.zero_submodule()
            if lx != rx:
                return False

        # This method is overridden for free modules in FreeModule_generic.
        # Hence we know self and other are not ambient, but they are contained
        # in the same ambient space
        return self.is_submodule(other) and other.is_submodule(self)

    def is_submodule(self, other):
        r"""
        Return ``True`` if ``self`` is a submodule of ``other``.

        EXAMPLES:

        Submodule testing over general rings is not guaranteed to work in
        all cases. However, it will raise an error when it is unable to
        determine containment.

        The zero module can always be tested::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y*z, x*z])])
            sage: N.zero_submodule().is_submodule(M)
            True
            sage: N.zero_submodule().is_submodule(N)
            True
            sage: M.zero_submodule().is_submodule(N)
            True

        It also respects which module it is constructed from::

            sage: Q = M.quotient_module(N)
            sage: Q.zero_submodule().is_submodule(M)
            False
            sage: Q.zero_submodule().is_submodule(N)
            False
            sage: M.zero_submodule().is_submodule(Q)
            False
            sage: N.zero_submodule().is_submodule(Q)
            False
        """
        if self is other:
            return True
        if not isinstance(other, Module_free_ambient):
            return False
        if self.base_ring() != other.base_ring():
            return False

        if not (self.ambient_module() == other.ambient_module()):
            if not (self.ambient_module() != other.ambient_module()):
                raise NotImplementedError("could not determine containment")
            return False
        if other is self.ambient_module():
            return True

        from sage.modules.quotient_module import QuotientModule_free_ambient
        if isinstance(other, QuotientModule_free_ambient):
            # if the relations agree we continue with the covers
            if isinstance(self, QuotientModule_free_ambient):
                if other.relations() != self.relations():
                    return False
                self = self.cover()
            else:
                if other.relations() != 0:
                    return False
            other = other.cover()

        if self.degree() != other.degree():
            return False
        R = self.base_ring()
        S = other.base_ring()
        if R != S:
            try:
                if not R.is_subring(S):
                    return False
            except NotImplementedError:
                if not R.fraction_field().is_subring(S):
                    raise NotImplementedError("could not determine if %s is a "
                                              "subring of %s" % (R, S))
        if not self.gens():
            # self is the zero module
            return True
        if not other.gens():
            # other is the zero module
            return False
        if self.ambient_module().is_submodule(other):
            return True

        raise NotImplementedError("could not determine containment")

    _submodule_class = LazyImport("sage.modules.submodule", "Submodule_free_ambient")

    def span(self, gens, base_ring=None, check=True, already_echelonized=False):
        r"""
        Return the `R`-span of ``gens``, where `R` is the ``base_ring``.

        The default `R` is the base ring of ``self``. Note that this span need
        not be a submodule of ``self``, nor even of the ambient space. It must,
        however, be contained in the ambient vector space, i.e., the
        ambient space tensored with the fraction field of `R`.

        INPUT:

        - ``gens`` -- list of vectors

        - ``base_ring`` -- (optional) a ring

        - ``check`` -- boolean (default: ``True``); whether or not to
          coerce entries of gens into base field

        - ``already_echelonized`` -- boolean (default: ``False``);
          set this if you know the gens are already in echelon form

        EXAMPLES::

            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace([[2, 3, 4]]); W
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 5 2]
            sage: W.span([[1, 1, 1]])
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 1 1]

        Over a general ring::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: M.span([vector([x - y, z]), vector([y*z, x*z])])
            Submodule of Ambient free module of rank 2 over the integral domain
             Multivariate Polynomial Ring in x, y, z over Rational Field
            Generated by the rows of the matrix:
            [x - y     z]
            [  y*z   x*z]

        Over a PID::

            sage: V = FreeModule(ZZ, 3)
            sage: W = V.submodule([V.gen(0)])
            sage: W.span([V.gen(1)])
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1 0]
            sage: W.submodule([V.gen(1)])
            Traceback (most recent call last):
            ...
            ArithmeticError: argument gens (= [(0, 1, 0)]) does not generate a submodule of self
            sage: V.span([[1,0,0], [1/5,4,0], [6,3/4,0]])
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1/5   0   0]
            [  0 1/4   0]

        It also works with other things than integers::

            sage: R.<x> = QQ[]
            sage: L = R^1
            sage: a = L.span([(1/x,)])
            sage: a
            Free module of degree 1 and rank 1
             over Univariate Polynomial Ring in x over Rational Field
            Echelon basis matrix:
            [1/x]
            sage: b=L.span([(1/x,)])
            sage: a(b.gens()[0])
            (1/x)
            sage: L2 = R^2
            sage: L2.span([[(x^2+x)/(x^2-3*x+2), 1/5], [(x^2+2*x)/(x^2-4*x+3), x]])
            Free module of degree 2 and rank 2
             over Univariate Polynomial Ring in x over Rational Field
            Echelon basis matrix:
            [x/(x^3 - 6*x^2 + 11*x - 6)  2/15*x^2 - 17/75*x - 1/75]
            [                         0 x^3 - 11/5*x^2 - 3*x + 4/5]

        Note that the ``base_ring`` can make a huge difference. We
        repeat the previous example over the fraction field of R and
        get a simpler vector space. ::

            sage: L2.span([[(x^2+x)/(x^2-3*x+2), 1/5], [(x^2+2*x)/(x^2-4*x+3), x]],
            ....:         base_ring=R.fraction_field())
            Vector space of degree 2 and dimension 2 over
             Fraction Field of Univariate Polynomial Ring in x over Rational Field
            Basis matrix:
            [1 0]
            [0 1]

        TESTS::

            sage: V = FreeModule(RDF, 3)
            sage: W = V.submodule([V.gen(0)])
            sage: W.span([V.gen(1)], base_ring=GF(7))
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [0 1 0]
            sage: v = V((1, pi, log(2))); v                                             # needs sage.symbolic
            (1.0, 3.141592653589793, 0.6931471805599453)
            sage: W.span([v], base_ring=GF(7))                                          # needs sage.rings.finite_rings sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: argument gens (= [(1.0, 3.141592653589793, 0.6931471805599453)]) is not compatible with base_ring (= Finite Field of size 7)
            sage: W = V.submodule([v])                                                  # needs sage.symbolic
            sage: W.span([V.gen(2)], base_ring=GF(7))
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [0 0 1]
        """
        if isinstance(gens, FreeModule_generic):
            gens = gens.gens()
        if base_ring is None or base_ring is self.base_ring():
            return self._submodule_class(self.ambient_module(), gens, check=check, already_echelonized=already_echelonized)

        # The base ring has changed
        try:
            M = self.ambient_module().change_ring(base_ring)
        except TypeError:
            raise ValueError("argument base_ring (= %s) is not compatible " % base_ring +
                             "with the base ring (= %s)" % self.base_ring())
        try:
            return M.span(gens)
        except TypeError:
            raise ValueError("argument gens (= %s) is not compatible " % gens +
                             "with base_ring (= %s)" % base_ring)

    def submodule(self, gens, check=True, already_echelonized=False):
        r"""
        Create the `R`-submodule of the ambient module with given generators,
        where `R` is the base ring of ``self``.

        INPUT:

        - ``gens`` -- list of free module elements or a free module

        - ``check`` -- boolean (default: ``True``); whether or not to verify
          that the gens are in ``self``

        OUTPUT:

        The submodule spanned by the vectors in the list ``gens``. The basis for
        the subspace is always put in reduced row echelon form (if possible).

        EXAMPLES:

        We create a submodule of `\ZZ^3`::

            sage: M = FreeModule(ZZ, 3)
            sage: B = M.basis()
            sage: W = M.submodule([B[0] + B[1], 2*B[1] - B[2]])
            sage: W
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  1  0]
            [ 0  2 -1]

        We create a submodule of a submodule::

            sage: W.submodule([3*B[0] + 3*B[1]])
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [3 3 0]

        We try to create a submodule that is not really a submodule,
        which results in an :exc:`ArithmeticError` exception::

            sage: W.submodule([B[0] - B[1]])
            Traceback (most recent call last):
            ...
            ArithmeticError: argument gens (= [(1, -1, 0)]) does not generate a submodule of self

        Next we create a submodule of a free module over the principal ideal
        domain `\QQ[x]`, which uses the general Hermite normal form functionality::

            sage: R = PolynomialRing(QQ, 'x'); x = R.gen()
            sage: M = FreeModule(R, 3)
            sage: B = M.basis()
            sage: W = M.submodule([x*B[0], 2*B[1] - x*B[2]]); W
            Free module of degree 3 and rank 2
             over Univariate Polynomial Ring in x over Rational Field
            Echelon basis matrix:
            [ x  0  0]
            [ 0  2 -x]
            sage: W.ambient_module()
            Ambient free module of rank 3 over the principal ideal domain
             Univariate Polynomial Ring in x over Rational Field

        Over a generic ring::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: A = S**2
            sage: A.submodule([vector([x - y,z]), vector([y*z, x*z])])
            Submodule of Ambient free module of rank 2 over the integral domain
             Multivariate Polynomial Ring in x, y, z over Rational Field
            Generated by the rows of the matrix:
            [x - y     z]
            [  y*z   x*z]
        """
        if isinstance(gens, Module_free_ambient):
            gens = gens.gens()
        V = self.span(gens, check=check, already_echelonized=already_echelonized)
        if check:
            if not V.is_submodule(self):
                raise ArithmeticError("argument gens (= %s) does not generate "
                                      "a submodule of self" % gens)
        return V

    def quotient_module(self, sub, check=True):
        r"""
        Return the quotient of ``self`` by the given subspace ``sub``.

        INPUT:

        - ``sub`` -- a submodule of ``self`` or something that can
          be turned into one via ``self.submodule(sub)``

        - ``check`` -- boolean (default: ``True``); whether or not to check that
          ``sub`` is a submodule

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
            sage: M.quotient(N)
            Quotient module by Submodule of Ambient free module of rank 2 over
             the integral domain Multivariate Polynomial Ring in x, y, z over Rational Field
            Generated by the rows of the matrix:
            [x - y     z]
            [  y*z   x*z]
        """
        if isinstance(sub, Module_free_ambient) and self.base_ring() != sub.base_ring():
            raise ValueError("base rings must be the same")
        if check and (not isinstance(sub, Module_free_ambient) or not sub.is_submodule(self)):
            try:
                sub = self.submodule(sub)
            except (TypeError, ArithmeticError):
                raise ArithmeticError("sub must be a subspace of self")
        from sage.modules.quotient_module import QuotientModule_free_ambient
        return QuotientModule_free_ambient(self, sub)

    def __truediv__(self, sub):
        """
        Return the quotient of ``self`` by the given submodule sub.

        EXAMPLES::

            sage: V1 = ZZ^2; W1 = V1.span([[1,2],[3,4]])
            sage: V1/W1
            Finitely generated module V/W over Integer Ring with invariants (2)
            sage: V2 = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W2 = V2.span([2*V2.0+4*V2.1, 9*V2.0+12*V2.1, 4*V2.2])
            sage: V2/W2
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
        """
        return self.quotient(sub, check=True)

    def free_resolution(self, *args, **kwds):
        r"""
        Return a free resolution of ``self``.

        For input options, see
        :class:`~sage.homology.free_resolution.FreeResolution`.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
            sage: res = N.free_resolution(); res                                        # needs sage.libs.singular
            S^2 <-- S^2 <-- 0
            sage: ascii_art(res.chain_complex())                                        # needs sage.libs.singular
                        [x - y   y*z]
                        [    z   x*z]
             0 <-- C_0 <-------------- C_1 <-- 0
        """
        from sage.rings.polynomial.multi_polynomial_libsingular import (
            MPolynomialRing_libsingular,
        )
        if isinstance(self.base_ring(), MPolynomialRing_libsingular):
            from sage.homology.free_resolution import FiniteFreeResolution_singular
            return FiniteFreeResolution_singular(self, *args, **kwds)

        if isinstance(self, FreeModule_generic):
            from sage.homology.free_resolution import FiniteFreeResolution_free_module
            return FiniteFreeResolution_free_module(self, *args, **kwds)

        raise NotImplementedError("the module must be a free module or "
                                  "have the base ring be a polynomial ring using Singular")

    def graded_free_resolution(self, *args, **kwds):
        r"""
        Return a graded free resolution of ``self``.

        For input options, see
        :class:`~sage.homology.graded_resolution.GradedFiniteFreeResolution`.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
            sage: N.graded_free_resolution(shifts=[1, -1])                              # needs sage.libs.singular
            S(-1)⊕S(1) <-- S(-2)⊕S(-3) <-- 0
            sage: N.graded_free_resolution(shifts=[2, 3])                               # needs sage.libs.singular
            S(-2)⊕S(-3) <-- S(-3)⊕S(-4) <-- 0

            sage: N = M.submodule([vector([x^3 - y^6, z^2]), vector([y * z, x])])
            sage: N.graded_free_resolution(degrees=[2, 1, 3], shifts=[2, 3])            # needs sage.libs.singular
            S(-2)⊕S(-3) <-- S(-6)⊕S(-8) <-- 0
        """
        from sage.rings.polynomial.multi_polynomial_libsingular import (
            MPolynomialRing_libsingular,
        )
        if isinstance(self.base_ring(), MPolynomialRing_libsingular):
            from sage.homology.graded_resolution import (
                GradedFiniteFreeResolution_singular,
            )
            return GradedFiniteFreeResolution_singular(self, *args, **kwds)

        if isinstance(self, FreeModule_generic):
            from sage.homology.graded_resolution import (
                GradedFiniteFreeResolution_free_module,
            )
            return GradedFiniteFreeResolution_free_module(self, *args, **kwds)

        raise NotImplementedError("the module must be a free module or "
                                  "have the base ring be a polynomial ring using Singular")


class FreeModule_generic(Module_free_ambient):
    """
    Base class for all free modules.

    INPUT:

    - ``base_ring`` -- a commutative ring

    - ``rank`` -- nonnegative integer

    - ``degree`` -- nonnegative integer

    - ``sparse`` -- boolean (default: ``False``)

    - ``coordinate_ring`` -- a ring containing ``base_ring``
      (default: equal to ``base_ring``)

    - ``category`` -- category (default: ``None``)

    If ``base_ring`` is a field, then the default category is the
    category of finite-dimensional vector spaces over that field;
    otherwise it is the category of finite-dimensional free modules
    over that ring.  In addition, the category is intersected with the
    category of finite enumerated sets if the ring is finite or the
    rank is 0.

    EXAMPLES::

        sage: PolynomialRing(QQ,3,'x')^3
        Ambient free module of rank 3 over the integral domain
         Multivariate Polynomial Ring in x0, x1, x2 over Rational Field

        sage: FreeModule(GF(7), 3).category()
        Category of enumerated finite dimensional vector spaces with basis over
         (finite enumerated fields and subquotients of monoids and quotients of semigroups)
        sage: V = QQ^4; V.category()
        Category of finite dimensional vector spaces with basis over
         (number fields and quotient fields and metric spaces)
        sage: V = GF(5)**20; V.category()
        Category of enumerated finite dimensional vector spaces with basis over
         (finite enumerated fields and subquotients of monoids and quotients of semigroups)
        sage: FreeModule(ZZ,3).category()
        Category of finite dimensional modules with basis over
         (Dedekind domains and euclidean domains
          and noetherian rings and infinite enumerated sets
          and metric spaces)
        sage: (QQ^0).category()
        Category of finite enumerated finite dimensional vector spaces with basis
         over (number fields and quotient fields and metric spaces)

    TESTS:

    Check that :issue:`17576` is fixed::

        sage: V = VectorSpace(RDF, 3)
        sage: v = vector(RDF, [1, 2, 3, 4])
        sage: v in V
        False
    """
    def __init__(self, base_ring, rank, degree, sparse=False,
                 coordinate_ring=None, category=None):
        """
        Create the free module of given rank ``rank`` over the given base
        ring ``base_ring``.

        TESTS::

            sage: M = FreeModule(ZZ, 20, sparse=False)
            sage: x = M.random_element()
            sage: type(x)
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>
            sage: M.element_class
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

            sage: N = FreeModule(ZZ, 20, sparse=True)
            sage: y = N.random_element()
            sage: type(y)
            <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
            sage: N.element_class
            <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
        """
        if base_ring not in CommutativeRings():
            warn("You are constructing a free module\n"
                 "over a noncommutative ring. Sage does not have a concept\n"
                 "of left/right and both sided modules, so be careful.\n"
                 "It's also not guaranteed that all multiplications are\n"
                 "done from the right side.")

        if coordinate_ring is None:
            coordinate_ring = base_ring

        if not hasattr(self, 'Element'):
            self.Element = element_class(coordinate_ring, sparse)

        rank = sage.rings.integer.Integer(rank)
        if rank < 0:
            raise ValueError("rank (=%s) must be nonnegative" % rank)

        Module_free_ambient.__init__(self, base_ring, degree=degree, sparse=sparse, category=category)
        self.__coordinate_ring = coordinate_ring
        self.__uses_ambient_inner_product = True
        self.__rank = rank
        self._gram_matrix = None

    def construction(self):
        """
        The construction functor and base ring for ``self``.

        EXAMPLES::

            sage: R = PolynomialRing(QQ, 3, 'x')
            sage: V = R^5
            sage: V.construction()
            (VectorFunctor, Multivariate Polynomial Ring in x0, x1, x2 over Rational Field)
        """
        from sage.categories.pushout import VectorFunctor
        if hasattr(self, '_inner_product_matrix'):
            return VectorFunctor(self.rank(), self.is_sparse(),
                                 self.inner_product_matrix()), self.base_ring()
        return VectorFunctor(self.rank(), self.is_sparse()), self.base_ring()

    # FIXME: what's the level of generality of FreeModuleHomspace?
    # Should there be a category for free modules accepting it as hom space?
    # See similar method for FreeModule_generic_field class
    def _Hom_(self, Y, category):
        from sage.modules.free_module_homspace import FreeModuleHomspace
        return FreeModuleHomspace(self, Y, category)

    def dense_module(self):
        """
        Return corresponding dense module.

        EXAMPLES:

        We first illustrate conversion with ambient spaces::

            sage: M = FreeModule(QQ, 3)
            sage: S = FreeModule(QQ, 3, sparse=True)
            sage: M.sparse_module()
            Sparse vector space of dimension 3 over Rational Field
            sage: S.dense_module()
            Vector space of dimension 3 over Rational Field
            sage: M.sparse_module() == S
            True
            sage: S.dense_module() == M
            True
            sage: M.dense_module() == M
            True
            sage: S.sparse_module() == S
            True

        Next we create a subspace::

            sage: M = FreeModule(QQ, 3, sparse=True)
            sage: V = M.span([ [1,2,3] ] ); V
            Sparse vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
            sage: V.sparse_module()
            Sparse vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
        """
        if self.is_sparse():
            return self._dense_module()
        return self

    def _dense_module(self):
        """
        Create a dense module with the same defining data as ``self``.

        N.B. This function is for internal use only! See ``dense_module`` for
        use.

        EXAMPLES::

            sage: M = FreeModule(Integers(8), 3)
            sage: S = FreeModule(Integers(8), 3, sparse=True)
            sage: M is S._dense_module()
            True
        """
        A = self.ambient_module().dense_module()
        return A.span(self.basis())

    def sparse_module(self):
        """
        Return the corresponding sparse module with the same defining
        data.

        EXAMPLES:

        We first illustrate conversion with ambient spaces::

            sage: M = FreeModule(Integers(8), 3)
            sage: S = FreeModule(Integers(8), 3, sparse=True)
            sage: M.sparse_module()
            Ambient sparse free module of rank 3 over Ring of integers modulo 8
            sage: S.dense_module()
            Ambient free module of rank 3 over Ring of integers modulo 8
            sage: M.sparse_module() is S
            True
            sage: S.dense_module() is M
            True
            sage: M.dense_module() is M
            True
            sage: S.sparse_module() is S
            True

        Next we convert a subspace::

            sage: M = FreeModule(QQ,3)
            sage: V = M.span([ [1,2,3] ] ); V
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
            sage: V.sparse_module()
            Sparse vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
        """
        if self.is_sparse():
            return self
        return self._sparse_module()

    def _sparse_module(self):
        """
        Create a sparse module with the same defining data as ``self``.

        N.B. This function is for internal use only! See ``sparse_module`` for
        use.

        EXAMPLES::

            sage: M = FreeModule(Integers(8), 3)
            sage: S = FreeModule(Integers(8), 3, sparse=True)
            sage: M._sparse_module() is S
            True
        """
        A = self.ambient_module().sparse_module()
        return A.span(self.basis())

    def _element_constructor_(self, x, coerce=True, copy=True, check=True):
        r"""
        Create an element of this free module from ``x``.

        The ``coerce`` and ``copy`` arguments are
        passed on to the underlying element constructor. If
        ``check`` is ``True``, confirm that the
        element specified by x does in fact lie in ``self``.

        .. NOTE::

           In the case of an inexact base ring (i.e. RDF), we don't
           verify that the element is in the subspace, even when
           ``check=True``, to account for numerical instability
           issues.

        EXAMPLES::

            sage: M = ZZ^4
            sage: M([1,-1,0,1])  #indirect doctest
            (1, -1, 0, 1)
            sage: M(0)
            (0, 0, 0, 0)

        ::

            sage: N = M.submodule([[1,0,0,0], [0,1,1,0]])
            sage: N([1,1,1,0])
            (1, 1, 1, 0)
            sage: N((3,-2,-2,0))
            (3, -2, -2, 0)
            sage: N((0,0,0,1))
            Traceback (most recent call last):
            ...
            TypeError: element (0, 0, 0, 1) is not in free module

        Beware that using check=False can create invalid results::

            sage: N((0,0,0,1), check=False)
            (0, 0, 0, 1)
            sage: N((0,0,0,1), check=False) in N
            True
        """
        if (isinstance(x, (int, sage.rings.integer.Integer)) and x == 0):
            return self.zero_vector()
        if isinstance(x, free_module_element.FreeModuleElement):
            if x.parent() is self:
                if copy:
                    return x.__copy__()
                else:
                    return x
            x = x.list()
        if check and self.coordinate_ring().is_exact():
            if isinstance(self, FreeModule_ambient):
                return self.element_class(self, x, coerce, copy)
            try:
                c = self.coordinates(x)
                R = self.base_ring()
                for d in c:
                    if d not in R:
                        raise ArithmeticError
            except ArithmeticError:
                raise TypeError("element {!r} is not in free module".format(x))
        return self.element_class(self, x, coerce, copy)

    def _eq(self, other):
        r"""
        Return if ``self`` is equal to ``other``.

        Ambient spaces are considered equal if they have the same
        rank, basering and inner product matrix.

        Modules in the same ambient space are partially ordered by inclusion.

        EXAMPLES::

            sage: L = IntegralLattice("U")
            sage: L is IntegralLattice("U")
            False
            sage: L._eq(IntegralLattice("U"))
            True
        """
        if not isinstance(other, FreeModule_generic):
            return False
        if self.rank() != other.rank():
            return False
        if self.base_ring() != other.base_ring():
            return False
        # We do not want to create an inner product matrix in memory if
        # self and other use the dot product
        if not (self._inner_product_is_dot_product()
                and other._inner_product_is_dot_product()):
            # This only affects free_quadratic_modules
            if self.inner_product_matrix() != other.inner_product_matrix():
                return False
        from sage.modules.quotient_module import FreeModule_ambient_field_quotient
        lq = isinstance(self, FreeModule_ambient_field_quotient)
        rq = isinstance(other, FreeModule_ambient_field_quotient)
        if lq or rq:
            # if the relations agree we continue with the covers.
            if lq:
                lx = self.relations()
                self = self.cover()
            else:
                lx = self.zero_submodule()
            if rq:
                rx = other.relations()
                other = other.cover()
            else:
                rx = other.zero_submodule()
            if lx != rx:
                return False
        if isinstance(self, FreeModule_ambient) and isinstance(other, FreeModule_ambient):
            return True
        # self and other are not ambient.
        # but they are contained in the same ambient space

        # We use self.echelonized_basis_matrix() == other.echelonized_basis_matrix()
        # with the matrix to avoid a circular reference.
        from sage.rings.integer_ring import IntegerRing
        if self.base_ring().is_field() or self.base_ring() is IntegerRing():
            # We know that the Hermite normal form is unique here.
            return self.echelonized_basis_matrix() == other.echelonized_basis_matrix()
        return self.is_submodule(other) and other.is_submodule(self)

    def is_submodule(self, other):
        r"""
        Return ``True`` if ``self`` is a submodule of ``other``.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3)
            sage: V = M.ambient_vector_space()
            sage: X = V.span([[1/2,1/2,0], [1/2,0,1/2]], ZZ)
            sage: Y = V.span([[1,1,1]], ZZ)
            sage: N = X + Y
            sage: M.is_submodule(X)
            False
            sage: M.is_submodule(Y)
            False
            sage: Y.is_submodule(M)
            True
            sage: N.is_submodule(M)
            False
            sage: M.is_submodule(N)
            True

            sage: M = FreeModule(ZZ, 2)
            sage: M.is_submodule(M)
            True
            sage: N = M.scale(2)
            sage: N.is_submodule(M)
            True
            sage: M.is_submodule(N)
            False
            sage: N = M.scale(1/2)
            sage: N.is_submodule(M)
            False
            sage: M.is_submodule(N)
            True

        Since :meth:`basis` is not implemented in general, submodule
        testing does not work for all PID's. However, trivial cases are
        already used (and useful) for coercion, e.g.::

            sage: QQ(1/2) * vector(ZZ['x']['y'], [1,2,3,4])
            (1/2, 1, 3/2, 2)
            sage: vector(ZZ['x']['y'], [1,2,3,4]) * QQ(1/2)
            (1/2, 1, 3/2, 2)

        TESTS::

            M = QQ^3 / [[1,2,3]]
            V = QQ^2
            V.is_submodule(M)
            False

            sage: M1 = QQ^3 / [[1,2,3]]
            sage: V1 = span(QQ, [(1,0,0)]) + M1.relations()
            sage: M2 = V1 / M1.relations()
            sage: M2.is_submodule(M1)  # Different ambient vector spaces
            False
        """
        if self is other:
            return True
        if not isinstance(other, FreeModule_generic):
            return False
        # Removing this try-except block changes the behavior of
        #   is_submodule for (QQ^2).is_submodule(CC^2)
        try:
            if self.ambient_vector_space() != other.ambient_vector_space():
                return False
            if other is other.ambient_vector_space():
                return True
        except AttributeError:
            # Not all free modules have an ambient_vector_space.
            pass
        from sage.modules.quotient_module import FreeModule_ambient_field_quotient
        if isinstance(other, FreeModule_ambient_field_quotient):
            # if the relations agree we continue with the covers.
            if isinstance(self, FreeModule_ambient_field_quotient):
                if other.relations() != self.relations():
                    return False
                self = self.cover()
            else:
                if other.relations() != 0:
                    return False
            other = other.cover()

        if other.rank() < self.rank():
            return False
        if other.degree() != self.degree():
            return False
        if self._inner_product_is_dot_product() and other._inner_product_is_dot_product():
            pass
        else:
            if self.inner_product_matrix() != other.inner_product_matrix():
                return False
        R = self.base_ring()
        S = other.base_ring()
        if R != S:
            try:
                if not R.is_subring(S):
                    return False
            except NotImplementedError:
                if not R.fraction_field().is_subring(S):
                    raise NotImplementedError("could not determine if %s is a "
                                              "subring of %s" % (R, S))
        # now R is a subring of S
        if other.is_ambient() and S.is_field():
            return True
        try:
            M = other.basis_matrix().solve_left(self.basis_matrix())
        except ValueError:
            return False
        except TypeError:
            # only if solve_left does not eat a matrix
            # else this is far to inefficient
            try:
                M = [list(other.basis_matrix().solve_left(self.basis_matrix()[i])) for i in range(self.basis_matrix().nrows())]
            except ValueError:
                return False
            from sage.misc.flatten import flatten
            return all(x in S for x in flatten(M))
        return all(x in S for x in M.list())

    def __iter__(self):
        r"""
        Return iterator over the elements of this free module.

        EXAMPLES::

            sage: V = VectorSpace(GF(4, 'a'), 2)                                        # needs sage.rings.finite_rings
            sage: [x for x in V]                                                        # needs sage.rings.finite_rings
            [(0, 0), (a, 0), (a + 1, 0), (1, 0), (0, a), (a, a), (a + 1, a), (1, a),
             (0, a + 1), (a, a + 1), (a + 1, a + 1), (1, a + 1), (0, 1), (a, 1),
             (a + 1, 1), (1, 1)]

        ::

            sage: W = V.subspace([V([1, 1])])                                           # needs sage.rings.finite_rings
            sage: [x for x in W]                                                        # needs sage.rings.finite_rings
            [(0, 0), (a, a), (a + 1, a + 1), (1, 1)]

        Free modules over enumerated infinite rings (i.e., those in the
        category :class:`InfiniteEnumeratedSets`) iterate over module
        elements ordered by (primarily) the 1-norm and (secondarily) the
        `\infty`-norm of their coordinate vectors::

            sage: it = iter(ZZ^3)
            sage: vs = [next(it) for _ in range(1000)]
            sage: vs[:50]
            [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1), (1, 1, 0), (-1, 1, 0), (1, -1, 0), (-1, -1, 0), (1, 0, 1), (-1, 0, 1), (1, 0, -1), (-1, 0, -1), (0, 1, 1), (0, -1, 1), (0, 1, -1), (0, -1, -1), (2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2), (1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1), (1, 1, -1), (-1, 1, -1), (1, -1, -1), (-1, -1, -1), (2, 1, 0), (-2, 1, 0), (2, -1, 0), (-2, -1, 0), (2, 0, 1), (-2, 0, 1), (2, 0, -1), (-2, 0, -1), (0, 2, 1), (0, -2, 1), (0, 2, -1), (0, -2, -1), (1, 2, 0), (1, -2, 0), (-1, 2, 0), (-1, -2, 0), (1, 0, 2)]
            sage: vs[775:825]
            [(1, 6, 1), (1, -6, 1), (1, 6, -1), (1, -6, -1), (-1, 6, 1), (-1, -6, 1), (-1, 6, -1), (-1, -6, -1), (2, 6, 0), (2, -6, 0), (-2, 6, 0), (-2, -6, 0), (1, 1, 6), (1, 1, -6), (-1, 1, 6), (-1, 1, -6), (1, -1, 6), (1, -1, -6), (-1, -1, 6), (-1, -1, -6), (2, 0, 6), (2, 0, -6), (-2, 0, 6), (-2, 0, -6), (0, 2, 6), (0, 2, -6), (0, -2, 6), (0, -2, -6), (7, 1, 0), (-7, 1, 0), (7, -1, 0), (-7, -1, 0), (7, 0, 1), (-7, 0, 1), (7, 0, -1), (-7, 0, -1), (0, 7, 1), (0, -7, 1), (0, 7, -1), (0, -7, -1), (1, 7, 0), (1, -7, 0), (-1, 7, 0), (-1, -7, 0), (1, 0, 7), (1, 0, -7), (-1, 0, 7), (-1, 0, -7), (0, 1, 7), (0, 1, -7)]

        TESTS::

            sage: V = VectorSpace(GF(2, 'a'), 2)
            sage: V.list()
            [(0, 0), (1, 0), (0, 1), (1, 1)]

        Test ``iter(ZZ^n)`` and the like::

            sage: it = iter(ZZ^5)
            sage: vs = [next(it) for _ in range(10000)]
            sage: _ = [v.set_immutable() for v in vs]
            sage: len(set(vs))  # no duplicates
            10000
            sage: onenorm = lambda v: sum(map(abs, v))
            sage: maxnorm = lambda v: max(map(abs, v))
            sage: norms = [(onenorm(v), maxnorm(v)) for v in vs]
            sage: all(n <= m for n,m in zip(norms, norms[1:]))  # ordered by (1-norm, max-norm)
            True
            sage: [sum(o <= b and m <= b for o,m in norms) for b in range(8)]  # did not miss any
            [1, 11, 61, 231, 681, 1683, 3653, 7183]
        """
        G = self.gens()
        if not G:
            yield self(0)
            return

        R = self.base_ring()

        if R in InfiniteEnumeratedSets():
            # This makes iter(ZZ^n) produce vectors in a "natural" order,
            # rather than only vectors with the first component nonzero.
            # Algorithm: Initial version which ordered by max-norm due to
            # Aleksei Udovenko, adapted by Lorenz Panny to order by 1-norm
            # primarily and by max-norm secondarily.
            def aux(length, norm, max_):
                if not 0 <= norm <= length * max_:
                    return  # there are no such vectors
                if norm == max_ == 0:
                    yield (0,) * length
                    return
                for pos in range(length):
                    for lnorm in range(norm - max_ + 1):
                        for lmax in range(max_):
                            for left in aux(pos, lnorm, lmax):
                                for rmax in range(max_ + 1):
                                    for right in aux(length - 1 - pos,
                                                     norm - max_ - lnorm, rmax):
                                        for mid in (+max_, -max_):
                                            yield left + (mid,) + right
            n = len(G)
            for norm in itertools.count(0):
                mm = (norm + n - 1) // n
                for max_ in range(mm, norm + 1):
                    for vec in aux(n, norm, max_):
                        yield self.linear_combination_of_basis(vec)
            assert False  # should loop forever

        iters = [iter(R) for _ in range(len(G))]
        for x in iters:
            next(x)     # put at 0
        zero = R.zero()
        v = [zero for _ in range(len(G))]
        n = 0
        z = self(0)
        yield z
        while n < len(G):
            try:
                v[n] = next(iters[n])
                yield self.linear_combination_of_basis(v)
                n = 0
            except StopIteration:
                iters[n] = iter(R)  # reset
                next(iters[n])     # put at 0
                v[n] = zero
                n += 1

    def cardinality(self):
        r"""
        Return the cardinality of the free module.

        OUTPUT: either an integer or ``+Infinity``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = FiniteField(9)
            sage: V = VectorSpace(k, 3)
            sage: V.cardinality()
            729
            sage: W = V.span([[1,2,1], [0,1,1]])
            sage: W.cardinality()
            81

            sage: R = IntegerModRing(12)
            sage: M = FreeModule(R, 2)
            sage: M.cardinality()
            144

            sage: (QQ^3).cardinality()
            +Infinity

        TESTS:

        Check that :issue:`22987` is fixed::

            sage: VectorSpace(QQ, 0).cardinality()
            1
        """
        if not self.rank():
            return sage.rings.integer.Integer(1)
        return self.base_ring().cardinality() ** self.rank()

    __len__ = cardinality  # for backward compatibility

    def basis(self):
        """
        Return the basis of this module.

        EXAMPLES::

            sage: FreeModule(Integers(12),3).basis()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        raise NotImplementedError

    def gens(self) -> tuple:
        """
        Return a tuple of basis elements of ``self``.

        EXAMPLES::

            sage: FreeModule(Integers(12),3).gens()
            ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        """
        return tuple(self.basis())

    def basis_matrix(self, ring=None):
        """
        Return the matrix whose rows are the basis for this free module.

        INPUT:

        - ``ring`` -- (default: ``self.coordinate_ring()``) a ring over
          which the matrix is defined

        EXAMPLES::

            sage: FreeModule(Integers(12), 3).basis_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

        ::

            sage: M = FreeModule(GF(7), 3).span([[2,3,4], [1,1,1]]); M
            Vector space of degree 3 and dimension 2 over Finite Field of size 7
            Basis matrix:
            [1 0 6]
            [0 1 2]
            sage: M.basis_matrix()
            [1 0 6]
            [0 1 2]

        ::

            sage: M = FreeModule(GF(7), 3).span_of_basis([[2,3,4], [1,1,1]])
            sage: M.basis_matrix()
            [2 3 4]
            [1 1 1]

        ::

            sage: M = FreeModule(QQ, 2).span_of_basis([[1,-1], [1,0]]); M
            Vector space of degree 2 and dimension 2 over Rational Field
            User basis matrix:
            [ 1 -1]
            [ 1  0]
            sage: M.basis_matrix()
            [ 1 -1]
            [ 1  0]

        TESTS:

        See :issue:`3699`::

            sage: K = FreeModule(ZZ, 2000)
            sage: I = K.basis_matrix()

        See :issue:`17585`::

            sage: ((ZZ^2)*2).basis_matrix().parent()
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: ((ZZ^2)*2).basis_matrix(RDF).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Real Double Field

            sage: M = (ZZ^2)*(1/2)
            sage: M.basis_matrix()
            [1/2   0]
            [  0 1/2]
            sage: M.basis_matrix().parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: M.basis_matrix(QQ).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: M.basis_matrix(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: matrix has denominators so can...t change to ZZ
        """
        try:
            A = self.__basis_matrix
        except AttributeError:
            MAT = sage.matrix.matrix_space.MatrixSpace(self.coordinate_ring(),
                                                       len(self.basis()),
                                                       self.degree(),
                                                       sparse=self.is_sparse())
            if self.is_ambient():
                A = MAT.identity_matrix()
            else:
                A = MAT(self.basis())
            A.set_immutable()
            self.__basis_matrix = A
        if ring is None or ring is A.base_ring():
            return A
        else:
            return A.change_ring(ring)

    def echelonized_basis_matrix(self):
        """
        The echelonized basis matrix (not implemented for this module).

        This example works because M is an ambient module. Submodule
        creation should exist for generic modules.

        EXAMPLES::

            sage: R = IntegerModRing(12)
            sage: S.<x,y> = R[]
            sage: M = FreeModule(S, 3)
            sage: M.echelonized_basis_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

        TESTS::

            sage: from sage.modules.free_module import FreeModule_generic
            sage: FreeModule_generic.echelonized_basis_matrix(M)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def matrix(self):
        """
        Return the basis matrix of this module, which is the matrix whose
        rows are a basis for this module.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 2)
            sage: M.matrix()
            [1 0]
            [0 1]
            sage: M.submodule([M.gen(0) + M.gen(1), M.gen(0) - 2*M.gen(1)]).matrix()
            [1 1]
            [0 3]
        """
        return self.basis_matrix()

    def direct_sum(self, other):
        """
        Return the direct sum of ``self`` and ``other`` as a free module.

        EXAMPLES::

            sage: V = (ZZ^3).span([[1/2,3,5], [0,1,-3]]); V
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1/2   0  14]
            [  0   1  -3]
            sage: W = (ZZ^3).span([[1/2,4,2]]); W
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1/2   4   2]
            sage: V.direct_sum(W)
            Free module of degree 6 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1/2   0  14   0   0   0]
            [  0   1  -3   0   0   0]
            [  0   0   0 1/2   4   2]
        """
        if not isinstance(other, FreeModule_generic):
            raise TypeError("other must be a free module")
        if other.base_ring() != self.base_ring():
            raise TypeError("base rings of self and other must be the same")
        return self.basis_matrix().block_sum(other.basis_matrix()).row_module(self.base_ring())

    def coordinates(self, v, check=True):
        """
        Write `v` in terms of the basis for ``self``.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also verify
          that `v` is really in ``self``

        OUTPUT: list

        Returns a list `c` such that if `B` is the basis for ``self``, then

        .. MATH::

            \\sum c_i B_i = v.

        If `v` is not in self, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 2); M0, M1 = M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0 - M1)
            [2, -1]
        """
        return self.coordinate_vector(v, check=check).list()

    def coordinate_vector(self, v, check=True):
        """
        Return the vector whose coefficients give `v` as a linear
        combination of the basis for ``self``.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also verify that
          `v` is really in ``self``

        OUTPUT: list

        EXAMPLES::

            sage: M = FreeModule(ZZ, 2); M0, M1 = M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinate_vector(2*M0 - M1)
            (2, -1)
        """
        raise NotImplementedError

    def coordinate_module(self, V):
        r"""
        Suppose ``V`` is a submodule of ``self`` (or a module commensurable
        with ``self``), and that ``self`` is a free module over `R` of rank
        `n`. Let `\phi` be the map from ``self`` to
        `R^n` that sends the basis vectors of ``self`` in order to the
        standard basis of `R^n`. This function returns the image
        `\phi(V)`.

        .. WARNING::

           If there is no integer `d` such that `dV` is a submodule
           of ``self``, then this function will give total nonsense.

        EXAMPLES:

        We illustrate this function with some
        `\ZZ`-submodules of `\QQ^3`::

            sage: V = (ZZ^3).span([[1/2,3,5], [0,1,-3]])
            sage: W = (ZZ^3).span([[1/2,4,2]])
            sage: V.coordinate_module(W)
            Free module of degree 2 and rank 1 over Integer Ring
            User basis matrix:
            [1 4]
            sage: V.0 + 4*V.1
            (1/2, 4, 2)

        In this example, the coordinate module isn't even in `\ZZ^3`::

            sage: W = (ZZ^3).span([[1/4,2,1]])
            sage: V.coordinate_module(W)
            Free module of degree 2 and rank 1 over Integer Ring
            User basis matrix:
            [1/2   2]

        The following more elaborate example illustrates using this
        function to write a submodule in terms of integral cuspidal modular
        symbols::

            sage: # needs sage.modular
            sage: M = ModularSymbols(54)
            sage: S = M.cuspidal_subspace()
            sage: K = S.integral_structure(); K
            Free module of degree 19 and rank 8 over Integer Ring
            Echelon basis matrix:
            [ 0  1  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            ...
            sage: L = M[0].integral_structure(); L
            Free module of degree 19 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 0  1  1  0 -2  1 -1  1 -1 -2  2  0  0  0  0  0  0  0  0]
            [ 0  0  3  0 -3  2 -1  2 -1 -4  2 -1 -2  1  2  0  0 -1  1]
            sage: K.coordinate_module(L)
            Free module of degree 8 and rank 2 over Integer Ring
            User basis matrix:
            [ 1  1  1 -1  1 -1  0  0]
            [ 0  3  2 -1  2 -1 -1 -2]
            sage: K.coordinate_module(L).basis_matrix() * K.basis_matrix()
            [ 0  1  1  0 -2  1 -1  1 -1 -2  2  0  0  0  0  0  0  0  0]
            [ 0  0  3  0 -3  2 -1  2 -1 -4  2 -1 -2  1  2  0  0 -1  1]
        """
        if not isinstance(V, FreeModule_generic):
            raise ValueError("V must be a free module")
        A = self.basis_matrix()
        A = A.matrix_from_columns(A.pivots()).transpose()
        B = V.basis_matrix()
        B = B.matrix_from_columns(self.basis_matrix().pivots()).transpose()
        S = A.solve_right(B).transpose()
        return (self.base_ring()**S.ncols()).span_of_basis(S.rows())

    def dimension(self):
        """
        Return the dimension of this free module.

        EXAMPLES::

            sage: M = FreeModule(FiniteField(19), 100)
            sage: W = M.submodule([M.gen(50)])
            sage: W.dimension()
            1
        """
        return self.rank()

    def codimension(self):
        """
        Return the codimension of this free module, which is the
        dimension of the ambient space minus the dimension of this
        free module.

        EXAMPLES::

            sage: M = Matrix(3, 4, range(12))
            sage: V = M.left_kernel(); V
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 1 -2  1]
            sage: V.dimension()
            1
            sage: V.codimension()
            2

        The codimension of an ambient space is always zero::

            sage: (QQ^10).codimension()
            0
        """
        return self.degree() - self.rank()

    def discriminant(self):
        """
        Return the discriminant of this free module.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3)
            sage: M.discriminant()
            1
            sage: W = M.span([[1,2,3]])
            sage: W.discriminant()
            14
            sage: W2 = M.span([[1,2,3], [1,1,1]])
            sage: W2.discriminant()
            6
        """
        return self.gram_matrix().determinant()

    def base_field(self):
        """
        Return the base field, which is the fraction field of the base ring
        of this module.

        EXAMPLES::

            sage: FreeModule(GF(3), 2).base_field()
            Finite Field of size 3
            sage: FreeModule(ZZ, 2).base_field()
            Rational Field
            sage: FreeModule(PolynomialRing(GF(7), 'x'), 2).base_field()
            Fraction Field of Univariate Polynomial Ring in x
             over Finite Field of size 7
        """
        return self.base_ring().fraction_field()

    def coordinate_ring(self):
        """
        Return the ring over which the entries of the vectors are
        defined.

        This is the same as :meth:`base_ring` unless an explicit basis
        was given over the fraction field.

        EXAMPLES::

            sage: M = ZZ^2
            sage: M.coordinate_ring()
            Integer Ring

        ::

            sage: M = (ZZ^2) * (1/2)
            sage: M.base_ring()
            Integer Ring
            sage: M.coordinate_ring()
            Rational Field

        ::

            sage: R.<x> = QQ[]
            sage: L = R^2
            sage: L.coordinate_ring()
            Univariate Polynomial Ring in x over Rational Field
            sage: L.span([(x,0), (1,x)]).coordinate_ring()
            Univariate Polynomial Ring in x over Rational Field
            sage: L.span([(x,0), (1,1/x)]).coordinate_ring()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
            sage: L.span([]).coordinate_ring()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self.__coordinate_ring

    def free_module(self):
        """
        Return this free module. (This is used by the
        ``FreeModule`` functor, and simply returns self.)

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3)
            sage: M.free_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        return self

    def gen(self, i=0):
        """
        Return the `i`-th generator for ``self``.

        Here `i` is between 0 and rank - 1, inclusive.

        INPUT:

        - ``i`` -- integer (default: 0)

        OUTPUT: `i`-th basis vector for ``self``

        EXAMPLES::

            sage: n = 5
            sage: V = QQ^n
            sage: B = [V.gen(i) for i in range(n)]
            sage: B
            [(1, 0, 0, 0, 0),
            (0, 1, 0, 0, 0),
            (0, 0, 1, 0, 0),
            (0, 0, 0, 1, 0),
            (0, 0, 0, 0, 1)]
            sage: V.gens() == tuple(B)
            True

        TESTS::

            sage: (QQ^3).gen(4/3)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert rational 4/3 to an integer
        """
        if i < 0 or i >= self.rank():
            raise ValueError("Generator %s not defined." % i)
        return self.basis()[i]

    def gram_matrix(self):
        r"""
        Return the gram matrix associated to this free module, defined to
        be `G = B*A*B.transpose()`, where A is the inner product matrix
        (induced from the ambient space), and B the basis matrix.

        EXAMPLES::

            sage: V = VectorSpace(QQ,4)
            sage: u = V([1/2,1/2,1/2,1/2])
            sage: v = V([0,1,1,0])
            sage: w = V([0,0,1,1])
            sage: M = span([u,v,w], ZZ)
            sage: M.inner_product_matrix() == V.inner_product_matrix()
            True
            sage: L = M.submodule_with_basis([u,v,w])
            sage: L.inner_product_matrix() == M.inner_product_matrix()
            True
            sage: L.gram_matrix()
            [1 1 1]
            [1 2 1]
            [1 1 2]
        """
        if self.is_ambient():
            return sage.matrix.matrix_space.MatrixSpace(self.base_ring(), self.degree(), sparse=True)(1)
        else:
            if self._gram_matrix is None:
                B = self.basis_matrix()
                self._gram_matrix = B * B.transpose()
            return self._gram_matrix

    def has_user_basis(self) -> bool:
        """
        Return ``True`` if the basis of this free module is
        specified by the user, as opposed to being the default echelon
        form.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.subspace([[2,'1/2', 1]])
            sage: W.has_user_basis()
            False
            sage: W = V.subspace_with_basis([[2,'1/2',1]])
            sage: W.has_user_basis()
            True
        """
        return False

    def hom(self, im_gens, codomain=None, **kwds):
        """
        Override the hom method to handle the case of morphisms given by left-multiplication
        of a matrix and the codomain is not given.

        EXAMPLES::

            sage: W = ZZ^2; W.hom(matrix(1, [1, 2]), side='right')
            Free module morphism defined as left-multiplication by the matrix
            [1 2]
            Domain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 1 over the principal ideal domain Integer Ring
            sage: V = QQ^2; V.hom(identity_matrix(2), side='right')
            Vector space morphism represented as left-multiplication by the matrix:
            [1 0]
            [0 1]
            Domain: Vector space of dimension 2 over Rational Field
            Codomain: Vector space of dimension 2 over Rational Field
        """
        from sage.structure.element import Matrix
        if codomain is None and isinstance(im_gens, Matrix):
            side = kwds.get("side", "left")
            n = im_gens.nrows() if side == "right" else im_gens.ncols()
            from sage.categories.pushout import pushout
            R = pushout(self.base_ring(), im_gens.base_ring())
            codomain = R**n
        return super().hom(im_gens, codomain, **kwds)

    def inner_product_matrix(self):
        """
        Return the default identity inner product matrix associated to this
        module.

        By definition this is the inner product matrix of the ambient
        space, hence may be of degree greater than the rank of the module.

        TODO: Differentiate the image ring of the inner product from the
        base ring of the module and/or ambient space. E.g. On an integral
        module over ZZ the inner product pairing could naturally take
        values in ZZ, QQ, RR, or CC.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3)
            sage: M.inner_product_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return sage.matrix.matrix_space.MatrixSpace(self.base_ring(), self.degree(), sparse=True)(1)

    def _inner_product_is_dot_product(self):
        """
        Return whether or not the inner product on this module is induced
        by the dot product on the ambient vector space. This is used
        internally by the inner_product function for optimization.

        EXAMPLES::

            sage: FreeModule(ZZ, 3)._inner_product_is_dot_product()
            True
            sage: FreeModule(ZZ, 3, inner_product_matrix=1)._inner_product_is_dot_product()
            True
            sage: FreeModule(ZZ, 2, inner_product_matrix=[1,0,-1,0])._inner_product_is_dot_product()
            False

        ::

            sage: M = FreeModule(QQ, 3)
            sage: M2 = M.span([[1,2,3]])
            sage: M2._inner_product_is_dot_product()
            True
        """
        return True

    def is_ambient(self):
        """
        Return ``False`` since this is not an ambient free module.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3).span([[1,2,3]]); M
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            sage: M.is_ambient()
            False
            sage: M = (ZZ^2).span([[1,0], [0,1]])
            sage: M
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 1]
            sage: M.is_ambient()
            False
            sage: M == M.ambient_module()
            True
        """
        return False

    def is_dense(self):
        """
        Return ``True`` if the underlying representation of
        this module uses dense vectors, and ``False`` otherwise.

        EXAMPLES::

            sage: FreeModule(ZZ, 2).is_dense()
            True
            sage: FreeModule(ZZ, 2, sparse=True).is_dense()
            False
        """
        return not self.is_sparse()

    def is_full(self):
        """
        Return ``True`` if the rank of this module equals its
        degree.

        EXAMPLES::

            sage: FreeModule(ZZ, 2).is_full()
            True
            sage: M = FreeModule(ZZ, 2).span([[1,2]])
            sage: M.is_full()
            False
        """
        return self.rank() == self.degree()

    def is_finite(self):
        """
        Return ``True`` if the underlying set of this free module is finite.

        EXAMPLES::

            sage: FreeModule(ZZ, 2).is_finite()
            False
            sage: FreeModule(Integers(8), 2).is_finite()
            True
            sage: FreeModule(ZZ, 0).is_finite()
            True
        """
        return self.base_ring().is_finite() or self.rank() == 0

    def ngens(self):
        """
        Return the number of basis elements of this free module.

        EXAMPLES::

            sage: FreeModule(ZZ, 2).ngens()
            2
            sage: FreeModule(ZZ, 0).ngens()
            0
            sage: FreeModule(ZZ, 2).span([[1,1]]).ngens()
            1
        """
        try:
            return self.__ngens
        except AttributeError:
            self.__ngens = self.rank()
        return self.__ngens

    def nonembedded_free_module(self):
        """
        Return an ambient free module that is isomorphic to this free
        module.

        Thus if this free module is of rank `n` over a ring
        `R`, then this function returns `R^n`, as an
        ambient free module.

        EXAMPLES::

            sage: FreeModule(ZZ, 2).span([[1,1]]).nonembedded_free_module()
            Ambient free module of rank 1 over the principal ideal domain Integer Ring
        """
        return FreeModule(self.base_ring(), self.rank())

    def random_element(self, prob=1.0, *args, **kwds):
        """
        Return a random element of ``self``.

        INPUT:

        - ``prob`` -- float. Each coefficient will be set to zero with
          probability `1-prob`. Otherwise coefficients will be chosen
          randomly from base ring (and may be zero).

        - ``*args``, ``**kwds`` -- passed on to the :func:`random_element`
          function of the base ring

        EXAMPLES::

            sage: M = FreeModule(ZZ, 2).span([[1, 1]])
            sage: v = M.random_element()
            sage: v.parent() is M
            True
            sage: v in M
            True

        Small entries are likely::

            sage: for i in [-2, -1, 0, 1, 2]:
            ....:     while vector([i, i]) != M.random_element():
            ....:         pass

        Large entries appear as well::

            sage: while abs(M.random_element()[0]) < 100:
            ....:     pass

        Passes extra positional or keyword arguments through::

            sage: all(i in range(5, 10) for i in M.random_element(1.0, 5, 10))
            True
        """
        rand = current_randstate().python_random().random
        R = self.base_ring()
        prob = float(prob)
        c = [0 if rand() > prob else R.random_element(*args, **kwds) for _ in range(self.rank())]
        return self.linear_combination_of_basis(c)

    def rank(self):
        """
        Return the rank of this free module.

        EXAMPLES::

            sage: FreeModule(Integers(6), 10000000).rank()
            10000000
            sage: FreeModule(ZZ, 2).span([[1,1], [2,2], [3,4]]).rank()
            2
        """
        return self.__rank

    def __bool__(self):
        """
        Return ``True`` if and only if the rank of this module is
        nonzero. In other words, this returns ``False`` for the zero
        module and ``True`` otherwise (apart from the exceptional case
        where the base ring is the zero ring).

        EXAMPLES::

            sage: bool(QQ^0)
            False
            sage: bool(QQ^1)
            True
            sage: M = Matrix(2, 3, range(6))
            sage: bool(M.right_kernel())
            True
            sage: bool(M.left_kernel())
            False

        When the base ring is the zero ring, we still look at the
        "rank" (which may not be mathematically meaningful)::

            sage: M = Integers(1)^4; M
            Ambient free module of rank 4 over Ring of integers modulo 1
            sage: M.rank()
            4
            sage: bool(M)
            True
            sage: M.cardinality()
            1
        """
        return bool(self.rank())

    def uses_ambient_inner_product(self):
        r"""
        Return ``True`` if the inner product on this module is
        the one induced by the ambient inner product.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 2)
            sage: W = M.submodule([[1,2]])
            sage: W.uses_ambient_inner_product()
            True
            sage: W.inner_product_matrix()
            [1 0]
            [0 1]

        ::

            sage: W.gram_matrix()
            [5]
        """
        return self.__uses_ambient_inner_product

    def are_linearly_dependent(self, vecs):
        """
        Return ``True`` if the vectors ``vecs`` are linearly dependent and
        ``False`` otherwise.

        EXAMPLES::

            sage: M = QQ^3
            sage: vecs = [M([1,2,3]), M([4,5,6])]
            sage: M.are_linearly_dependent(vecs)
            False
            sage: vecs.append(M([3,3,3]))
            sage: M.are_linearly_dependent(vecs)
            True

            sage: R.<x> = QQ[]
            sage: M = FreeModule(R, 2)
            sage: vecs = [M([x^2+1, x+1]), M([x+2, 2*x+1])]
            sage: M.are_linearly_dependent(vecs)
            False
            sage: vecs.append(M([-2*x+1, -2*x^2+1]))
            sage: M.are_linearly_dependent(vecs)
            True
        """
        from sage.matrix.constructor import matrix
        A = matrix(vecs)
        A.echelonize()
        return any(row.is_zero() for row in A.rows())

    def _magma_init_(self, magma):
        """
        EXAMPLES::

            sage: magma(QQ^9)                                   # optional - magma
            Full Vector space of degree 9 over Rational Field
            sage: (QQ^9)._magma_init_(magma)                    # optional - magma
            'RSpace(_sage_[...],9)'

        ::

            sage: magma(Integers(8)^2)                          # optional - magma
            Full RSpace of degree 2 over IntegerRing(8)
            sage: magma(FreeModule(QQ['x'], 2))                 # optional - magma
            Full RSpace of degree 2 over Univariate Polynomial Ring in x over Rational Field

        ::

            sage: # optional - magma
            sage: A = matrix([[1,0],[0,-1]])
            sage: M = FreeModule(ZZ,2,inner_product_matrix=A); M
            Ambient free quadratic module of rank 2 over the principal ideal domain Integer Ring
            Inner product matrix:
            [ 1  0]
            [ 0 -1]
            sage: M._magma_init_(magma)
            'RSpace(_sage_[...],2,_sage_ref...)'
            sage: m = magma(M); m
            Full RSpace of degree 2 over Integer Ring
            Inner Product Matrix:
            [ 1  0]
            [ 0 -1]
            sage: m.Type()
            ModTupRng
            sage: m.sage()
            Ambient free quadratic module of rank 2 over the principal ideal domain Integer Ring
            Inner product matrix:
            [ 1  0]
            [ 0 -1]
            sage: m.sage() is M
            True

        Now over a field::

            sage: # optional - magma
            sage: N = FreeModule(QQ,2,inner_product_matrix=A); N
            Ambient quadratic space of dimension 2 over Rational Field
            Inner product matrix:
            [ 1  0]
            [ 0 -1]
            sage: n = magma(N); n
            Full Vector space of degree 2 over Rational Field
            Inner Product Matrix:
            [ 1  0]
            [ 0 -1]
            sage: n.Type()
            ModTupFld
            sage: n.sage()
            Ambient quadratic space of dimension 2 over Rational Field
            Inner product matrix:
            [ 1  0]
            [ 0 -1]
            sage: n.sage() is N
            True

        How about some inexact fields::

            sage: # optional - magma, needs sage.symbolic
            sage: v = vector(RR, [1, pi, 5/6])
            sage: F = v.parent()
            sage: M = magma(F); M
            Full Vector space of degree 3 over Real field of precision 15
            sage: M.Type()
            ModTupFld
            sage: m = M.sage(); m
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: m is F
            True

        For interval fields, we can convert to Magma but there is no
        interval field in Magma so we cannot convert back::

            sage: # optional - magma, needs sage.symbolic
            sage: v = vector(RealIntervalField(100), [1, pi, 0.125])
            sage: F = v.parent()
            sage: M = magma(v.parent()); M
            Full Vector space of degree 3 over Real field of precision 30
            sage: M.Type()
            ModTupFld
            sage: m = M.sage(); m
            Vector space of dimension 3 over Real Field with 100 bits of precision
            sage: m is F
            False
        """
        K = magma(self.base_ring())
        if not self._inner_product_is_dot_product():
            M = magma(self.inner_product_matrix())
            return "RSpace(%s,%s,%s)" % (K.name(), self.rank(), M._ref())
        else:
            return "RSpace(%s,%s)" % (K.name(), self.rank())

    def _macaulay2_(self, macaulay2=None):
        r"""
        EXAMPLES::

            sage: R = QQ^2
            sage: macaulay2(R)          # optional - macaulay2
              2
            QQ
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2
        if hasattr(self, '_inner_product_matrix'):
            raise NotImplementedError
        else:
            return macaulay2(self.base_ring())**self.rank()

    def scale(self, other):
        """
        Return the product of this module by the number other, which is the
        module spanned by other times each basis vector.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3)
            sage: M.scale(2)
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [2 0 0]
            [0 2 0]
            [0 0 2]

        ::

            sage: a = QQ('1/3')
            sage: M.scale(a)
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1/3   0   0]
            [  0 1/3   0]
            [  0   0 1/3]
        """
        if other == 0:
            return self.zero_submodule()
        if other == 1 or other == -1:
            return self
        return self.span([v * other for v in self.basis()])

    def __radd__(self, other):
        """
        EXAMPLES::

            sage: int(0) + QQ^3
            Vector space of dimension 3 over Rational Field
            sage: sum([QQ^3, QQ^3])
            Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        if other == 0:
            return self
        else:
            raise TypeError

    def _mul_(self, other, switch_sides=False):
        r"""
        Multiplication of the basis by ``other``.

        EXAMPLES::

            sage: A = ZZ^3
            sage: A * 3
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [3 0 0]
            [0 3 0]
            [0 0 3]

            sage: V = A.span([A([1,2,2]), A([-1,0,2])])
            sage: 2 * V
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 2  0 -4]
            [ 0  4  8]

            sage: m = matrix(3, range(9))
            sage: A * m
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 3  0 -3]
            [ 0  1  2]
            sage: m * A
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 3  0 -3]
            [ 0  1  2]

        TESTS:

        Check that :issue:`17705` is fixed::

            sage: V = GF(2)^2
            sage: W = V.subspace([[1, 0]])
            sage: x = matrix(GF(2), [[1, 1], [0, 1]])
            sage: W*x
            Vector space of degree 2 and dimension 1
             over Finite Field of size 2
            Basis matrix:
            [1 1]
        """
        B = self.basis_matrix()
        B = other * B if switch_sides else B * other
        return self.span(B.rows())

    def relations(self):
        """
        Return the module of relations of ``self``.

        EXAMPLES::

            sage: V = GF(2)^2
            sage: V.relations() == V.zero_submodule()
            True
            sage: W = V.subspace([[1, 0]])
            sage: W.relations() == V.zero_submodule()
            True

            sage: Q = V / W
            sage: Q.relations() == W
            True
        """
        return self.zero_submodule()


class FreeModule_generic_domain(FreeModule_generic):
    """
    Base class for free modules over an integral domain.
    """
    def __init__(self, base_ring, rank, degree, sparse=False, coordinate_ring=None, category=None):
        """
        Create a free module over an integral domain.

        EXAMPLES::

            sage: FreeModule(ZZ, 2)
            Ambient free module of rank 2
             over the principal ideal domain Integer Ring
            sage: FreeModule(PolynomialRing(GF(7), 'x'), 2)
            Ambient free module of rank 2
             over the principal ideal domain Univariate Polynomial Ring in x
              over Finite Field of size 7
        """
        FreeModule_generic.__init__(self, base_ring, rank, degree, sparse, coordinate_ring, category=category)

    def __add__(self, other):
        r"""
        Return the sum of ``self`` and ``other``, where both ``self`` and
        ``other`` must be submodules of the ambient vector space.

        EXAMPLES:

        We add two vector spaces::

            sage: V  = VectorSpace(QQ, 3)
            sage: W  = V.subspace([V([1,1,0])])
            sage: W2 = V.subspace([V([1,-1,0])])
            sage: W + W2
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]

        We add two free `\ZZ`-modules.

        ::

            sage: M = FreeModule(ZZ, 3)
            sage: W = M.submodule([M([1,0,2])])
            sage: W2 = M.submodule([M([2,0,-4])])
            sage: W + W2
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0 2]
            [0 0 8]

        We can also add free `\ZZ`-modules embedded
        non-integrally into an ambient space.

        ::

            sage: V = VectorSpace(QQ, 3)
            sage: W = M.span([1/2*V.0 - 1/3*V.1])

        Here the command ``M.span(...)`` creates the span of
        the indicated vectors over the base ring of `M`.

        ::

            sage: W2 = M.span([1/3*V.0 + V.1])
            sage: W + W2
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1/6  7/3    0]
            [   0 11/3    0]

        We add two modules over `\ZZ`::

            sage: A = Matrix(ZZ, 3, 3, [3, 0, -1, 0, -2, 0, 0, 0, -2])
            sage: V = (A+2).kernel()
            sage: W = (A-3).kernel()
            sage: V+W
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [5 0 0]
            [0 1 0]
            [0 0 1]

        We add a module to 0::

            sage: ZZ^3 + 0
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        if not isinstance(other, FreeModule_generic):
            if other == 0:
                return self
            raise TypeError("other (=%s) must be a free module" % other)
        if not (self.ambient_vector_space() == other.ambient_vector_space()):
            raise TypeError("ambient vector spaces must be equal")
        return self.span(self.basis() + other.basis())


class FreeModule_generic_pid(FreeModule_generic_domain):
    """
    Base class for all free modules over a PID.
    """
    def __init__(self, base_ring, rank, degree, sparse=False, coordinate_ring=None, category=None):
        """
        Create a free module over a PID.

        EXAMPLES::

            sage: FreeModule(ZZ, 2)
            Ambient free module of rank 2
             over the principal ideal domain Integer Ring
            sage: FreeModule(PolynomialRing(GF(7), 'x'), 2)
            Ambient free module of rank 2
             over the principal ideal domain Univariate Polynomial Ring in x
              over Finite Field of size 7
        """
        super().__init__(base_ring, rank, degree, sparse, coordinate_ring, category=category)

    def index_in(self, other):
        """
        Return the lattice index [other:self] of ``self`` in other, as an
        element of the base field. When ``self`` is contained in other, the
        lattice index is the usual index. If the index is infinite, then
        this function returns infinity.

        EXAMPLES::

            sage: L1 = span([[1,2]], ZZ)
            sage: L2 = span([[3,6]], ZZ)
            sage: L2.index_in(L1)
            3

        Note that the free modules being compared need not be integral.

        ::

            sage: L1 = span([['1/2','1/3'], [4,5]], ZZ)
            sage: L2 = span([[1,2], [3,4]], ZZ)
            sage: L2.index_in(L1)
            12/7
            sage: L1.index_in(L2)
            7/12
            sage: L1.discriminant() / L2.discriminant()
            49/144

        The index of a lattice of infinite index is infinite.

        ::

            sage: L1 = FreeModule(ZZ, 2)
            sage: L2 = span([[1,2]], ZZ)
            sage: L2.index_in(L1)
            +Infinity
        """
        if not isinstance(other, FreeModule_generic):
            raise TypeError("other must be a free module")

        if self.ambient_vector_space() != other.ambient_vector_space():
            raise ArithmeticError("self and other must be embedded in the same ambient space.")

        if self.base_ring() != other.base_ring():
            raise NotImplementedError("lattice index only defined for modules over the same base ring.")

        if other.base_ring().is_field():
            if self == other:
                return sage.rings.integer.Integer(1)
            else:
                if self.is_subspace(other):
                    return sage.rings.infinity.infinity
            raise ArithmeticError("self must be contained in the vector space spanned by other.")

        C = [other.coordinates(b) for b in self.basis()]

        if self.rank() < other.rank():
            return sage.rings.infinity.infinity

        a = sage.matrix.matrix_space.MatrixSpace(self.base_field(), self.rank())(C).determinant()
        if isinstance(self.base_ring(), sage.rings.integer_ring.IntegerRing_class):
            return a.abs()
        elif isinstance(self.base_ring, sage.rings.abc.Order):
            return self.base_ring().ideal(a).norm()
        else:
            raise NotImplementedError

    def intersection(self, other):
        r"""
        Return the intersection of ``self`` and ``other``.

        EXAMPLES:

        We intersect two submodules one of which is clearly
        contained in the other::

            sage: A = ZZ^2
            sage: M1 = A.span([[1,1]])
            sage: M2 = A.span([[3,3]])
            sage: M1.intersection(M2)
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [3 3]
            sage: M1.intersection(M2) is M2
            True

        We intersection two submodules of `\ZZ^3` of rank
        `2`, whose intersection has rank `1`::

            sage: A = ZZ^3
            sage: M1 = A.span([[1,1,1], [1,2,3]])
            sage: M2 = A.span([[2,2,2], [1,0,0]])
            sage: M1.intersection(M2)
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [2 2 2]

        We compute an intersection of two `\ZZ`-modules that
        are not submodules of `\ZZ^2`::

            sage: A = ZZ^2
            sage: M1 = A.span([[1,2]]).scale(1/6)
            sage: M2 = A.span([[1,2]]).scale(1/15)
            sage: M1.intersection(M2)
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1/3 2/3]

        We intersect a `\ZZ`-module with a `\QQ`-vector space::

            sage: A = ZZ^3
            sage: L = ZZ^3
            sage: V = QQ^3
            sage: W = L.span([[1/2,0,1/2]])
            sage: K = V.span([[1,0,1], [0,0,1]])
            sage: W.intersection(K)
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1/2   0 1/2]
            sage: K.intersection(W)
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1/2   0 1/2]

        We intersect two modules over the ring of integers of a number field::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: L.<w> = NumberField(x^2 - x + 2)
            sage: OL = L.ring_of_integers()
            sage: V = L**3
            sage: W1 = V.span([[0,w/5,0], [1,0,-1/17]], OL)
            sage: W2 = V.span([[0,(1-w)/5,0]], OL)
            sage: W1.intersection(W2)
            Free module of degree 3 and rank 1 over Maximal Order generated by w
             in Number Field in w with defining polynomial x^2 - x + 2
            Echelon basis matrix:
            [  0 2/5   0]

        TESTS:

        Check that :issue:`24702` is fixed::

            sage: L = FreeQuadraticModule(ZZ,2,matrix.identity(2))
            sage: S1 = L.submodule([(1,0)])
            sage: S2 = L.submodule([(0,1)])
            sage: S1.intersection(S2).ambient_module() == S1.ambient_module()
            True
        """
        if not isinstance(other, FreeModule_generic):
            raise TypeError("other must be a free module")

        if self.ambient_vector_space() != other.ambient_vector_space():
            raise ArithmeticError("self and other must be embedded in the same ambient space.")

        if self.base_ring() != other.base_ring():
            if other.base_ring().is_field():
                return other.intersection(self)
            raise NotImplementedError("intersection of modules over different base rings (neither a field) is not implemented.")

        # dispense with the three easy cases
        if self == self.ambient_vector_space() or other.is_submodule(self):
            return other
        elif other == other.ambient_vector_space() or self.is_submodule(other):
            return self
        elif self.rank() == 0 or other.rank() == 0:
            if self.base_ring().is_field():
                return other.zero_submodule()
            else:
                return self.zero_submodule()

        # standard algorithm for computing intersection of general submodule
        if self.dimension() <= other.dimension():
            V1 = self
            V2 = other
        else:
            V1 = other
            V2 = self
        A1 = V1.basis_matrix()
        A2 = V2.basis_matrix()
        S = A1.stack(A2)
        K = S.integer_kernel(self.base_ring()).basis_matrix()
        n = int(V1.dimension())
        K = K.matrix_from_columns(range(n))
        B = K * A1
        return self.span(B)

    def __and__(self, other):
        r"""
        Return the intersection of ``self`` and ``other``.

        See :meth:`intersection`.

        EXAMPLES:

        We intersect two submodules one of which is clearly
        contained in the other::

            sage: A = ZZ^2
            sage: M1 = A.span([[1,1]])
            sage: M2 = A.span([[3,3]])
            sage: M1 & M2
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [3 3]
            sage: M1 & M2 is M2
            True

        We intersection two submodules of `\ZZ^3` of rank
        `2`, whose intersection has rank `1`::

            sage: A = ZZ^3
            sage: M1 = A.span([[1,1,1], [1,2,3]])
            sage: M2 = A.span([[2,2,2], [1,0,0]])
            sage: M1 & M2
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [2 2 2]
        """
        return self.intersection(other)

    def zero_submodule(self):
        """
        Return the zero submodule of this module.

        EXAMPLES::

            sage: V = FreeModule(ZZ,2)
            sage: V.zero_submodule()
            Free module of degree 2 and rank 0 over Integer Ring
            Echelon basis matrix:
            []
        """
        return self.submodule([], check=False, already_echelonized=True)

    def denominator(self):
        """
        The denominator of the basis matrix of ``self`` (i.e. the LCM of the
        coordinate entries with respect to the basis of the ambient
        space).

        EXAMPLES::

            sage: V = QQ^3
            sage: L = V.span([[1,1/2,1/3], [-1/5,2/3,3]],ZZ)
            sage: L
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1/5 19/6 37/3]
            [   0 23/6 46/3]
            sage: L.denominator()
            30
        """
        return self.basis_matrix().denominator()

    def index_in_saturation(self):
        r"""
        Return the index of this module in its saturation, i.e., its
        intersection with `R^n`.

        EXAMPLES::

            sage: W = span([[2,4,6]], ZZ)
            sage: W.index_in_saturation()
            2
            sage: W = span([[1/2,1/3]], ZZ)
            sage: W.index_in_saturation()
            1/6
        """
        # TODO: There is probably a much faster algorithm in this case.
        return self.index_in(self.saturation())

    def saturation(self):
        r"""
        Return the saturated submodule of `R^n` that spans the same
        vector space as ``self``.

        EXAMPLES:

        We create a 1-dimensional lattice that is obviously not
        saturated and saturate it.

        ::

            sage: L = span([[9,9,6]], ZZ); L
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [9 9 6]
            sage: L.saturation()
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [3 3 2]

        We create a lattice spanned by two vectors, and saturate.
        Computation of discriminants shows that the index of lattice in its
        saturation is `3`, which is a prime of congruence between
        the two generating vectors.

        ::

            sage: L = span([[1,2,3], [4,5,6]], ZZ)
            sage: L.saturation()
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: L.discriminant()
            54
            sage: L.saturation().discriminant()
            6

        Notice that the saturation of a non-integral lattice `L` is
        defined, but the result is integral hence does not contain
        `L`::

            sage: L = span([['1/2',1,3]], ZZ)
            sage: L.saturation()
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1 2 6]

        TESTS:

        We check that :issue:`24702` is fixed::

            sage: L = FreeQuadraticModule(ZZ,1,matrix.identity(1))
            sage: S = 2*L
            sage: S.saturation().ambient_module() == L
            True
        """
        R = self.base_ring()
        if R.is_field():
            return self
        try:
            A, _ = self.basis_matrix()._clear_denom()
            S = self.span(A.saturation())
        except AttributeError:
            # fallback in case _clear_denom isn't written
            V = self.vector_space()
            A = self.ambient_module()
            S = V.intersection(A)
        # Return exactly self if it is already saturated.
        return self if self == S else S

    def span_of_basis(self, basis, base_ring=None, check=True, already_echelonized=False):
        r"""
        Return the free R-module with the given basis, where R is the base
        ring of ``self`` or user specified base_ring.

        Note that this R-module need not be a submodule of self, nor even
        of the ambient space. It must, however, be contained in the ambient
        vector space, i.e., the ambient space tensored with the fraction
        field of R.

        EXAMPLES::

            sage: M = FreeModule(ZZ,3)
            sage: W = M.span_of_basis([M([1,2,3])])

        Next we create two free `\ZZ`-modules, neither of
        which is a submodule of `W`.

        ::

            sage: W.span_of_basis([M([2,4,0])])
            Free module of degree 3 and rank 1 over Integer Ring
            User basis matrix:
            [2 4 0]

        The following module isn't in the ambient module `\ZZ^3`
        but is contained in the ambient vector space `\QQ^3`::

            sage: V = M.ambient_vector_space()
            sage: W.span_of_basis([ V([1/5,2/5,0]), V([1/7,1/7,0]) ])
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1/5 2/5   0]
            [1/7 1/7   0]

        Of course the input basis vectors must be linearly independent::

            sage: W.span_of_basis([ [1,2,0], [2,4,0] ])
            Traceback (most recent call last):
            ...
            ValueError: The given basis vectors must be linearly independent.
        """
        if isinstance(basis, FreeModule_generic):
            basis = basis.gens()
        if base_ring is None or base_ring == self.base_ring():
            try:
                if self.is_dense():
                    from sage.modules.free_module_integer import (
                        FreeModule_submodule_with_basis_integer,
                    )
                    return FreeModule_submodule_with_basis_integer(self.ambient_module(),
                                                                   basis=basis, check=check,
                                                                   already_echelonized=already_echelonized,
                                                                   lll_reduce=False)
            except TypeError:
                pass

            return FreeModule_submodule_with_basis_pid(
                self.ambient_module(), basis=basis, check=check,
                already_echelonized=already_echelonized)
        else:
            try:
                M = self.change_ring(base_ring)
            except TypeError:
                raise ValueError("Argument base_ring (= %s) is not compatible " % base_ring +
                                 "with the base ring (= %s)." % self.base_ring())
            try:
                return M.span_of_basis(basis)
            except TypeError:
                raise ValueError("Argument gens (= %s) is not compatible " % basis +
                                 "with base_ring (= %s)." % base_ring)

    def submodule_with_basis(self, basis, check=True, already_echelonized=False):
        r"""
        Create the `R`-submodule of the ambient vector space with given
        basis, where `R` is the base ring of ``self``.

        INPUT:

        - ``basis`` -- list of linearly independent vectors

        - ``check`` -- whether or not to verify that each gen is in
          the ambient vector space

        OUTPUT: ``FreeModule``; the `R`-submodule with given basis

        EXAMPLES:

        First we create a submodule of `\\ZZ^3`::

            sage: M = FreeModule(ZZ, 3)
            sage: B = M.basis()
            sage: N = M.submodule_with_basis([B[0]+B[1], 2*B[1]-B[2]])
            sage: N
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [ 1  1  0]
            [ 0  2 -1]

        A list of vectors in the ambient vector space may fail to generate
        a submodule.

        ::

            sage: V = M.ambient_vector_space()
            sage: X = M.submodule_with_basis([ V(B[0]+B[1])/2, V(B[1]-B[2])/2])
            Traceback (most recent call last):
            ...
            ArithmeticError: The given basis does not generate a submodule of self.

        However, we can still determine the R-span of vectors in the
        ambient space, or over-ride the submodule check by setting check to
        False.

        ::

            sage: X = V.span([ V(B[0]+B[1])/2, V(B[1]-B[2])/2 ], ZZ)
            sage: X
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1/2    0   1/2]
            [   0  1/2  -1/2]
            sage: Y = M.submodule([ V(B[0]+B[1])/2, V(B[1]-B[2])/2 ], check=False)
            sage: X == Y
            True

        Next we try to create a submodule of a free module over the
        principal ideal domain `\QQ[x]`, using our general Hermite normal form implementation::

            sage: R = PolynomialRing(QQ, 'x'); x = R.gen()
            sage: M = FreeModule(R, 3)
            sage: B = M.basis()
            sage: W = M.submodule_with_basis([x*B[0], 2*B[0]- x*B[2]]); W
            Free module of degree 3 and rank 2 over Univariate Polynomial Ring in x over Rational Field
            User basis matrix:
            [ x  0  0]
            [ 2  0 -x]
        """
        V = self.span_of_basis(basis=basis, check=check, already_echelonized=already_echelonized)
        if check:
            if not V.is_submodule(self):
                raise ArithmeticError("The given basis does not generate a submodule of self.")
        return V

    def vector_space_span(self, gens, check=True):
        r"""
        Create the vector subspace of the ambient vector space with given
        generators.

        INPUT:

        - ``gens`` -- list of vector in ``self``

        - ``check`` -- whether or not to verify that each gen
          is in the ambient vector space

        OUTPUT: a vector subspace

        EXAMPLES:

        We create a `2`-dimensional subspace of `\QQ^3`.

        ::

            sage: V = VectorSpace(QQ, 3)
            sage: B = V.basis()
            sage: W = V.vector_space_span([B[0]+B[1], 2*B[1]-B[2]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0  1/2]
            [   0    1 -1/2]

        We create a subspace of a vector space over
        `\QQ(i)`.

        ::

            sage: R.<x> = QQ[]
            sage: K = NumberField(x^2 + 1, 'a'); a = K.gen()                            # needs sage.rings.number_field
            sage: V = VectorSpace(K, 3)                                                 # needs sage.rings.number_field
            sage: V.vector_space_span([2*V.gen(0) + 3*V.gen(2)])                        # needs sage.rings.number_field
            Vector space of degree 3 and dimension 1
             over Number Field in a with defining polynomial x^2 + 1
            Basis matrix:
            [  1   0 3/2]

        We use the ``vector_space_span`` command to create a
        vector subspace of the ambient vector space of a submodule of
        `\ZZ^3`.

        ::

            sage: M = FreeModule(ZZ, 3)
            sage: W = M.submodule([M([1,2,3])])
            sage: W.vector_space_span([M([2,3,4])])
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [  1 3/2   2]
        """
        if isinstance(gens, FreeModule_generic):
            gens = gens.gens()
        return FreeModule_submodule_field(self.ambient_vector_space(), gens, check=check)

    def vector_space_span_of_basis(self, basis, check=True):
        """
        Create the vector subspace of the ambient vector space with given
        basis.

        INPUT:

        - ``basis`` -- list of linearly independent vectors

        - ``check`` -- whether or not to verify that each gen is in
          the ambient vector space

        OUTPUT: a vector subspace with user-specified basis

        EXAMPLES::

            sage: V = VectorSpace(QQ, 3)
            sage: B = V.basis()
            sage: W = V.vector_space_span_of_basis([B[0] + B[1], 2*B[1] - B[2]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [ 1  1  0]
            [ 0  2 -1]
        """
        return FreeModule_submodule_with_basis_field(self.ambient_vector_space(), basis, check=check)

    def quotient_module(self, sub, check=True, **kwds):
        """
        Return the quotient of ``self`` by the given submodule sub.

        INPUT:

        - ``sub`` -- a submodule of ``self``, or something that can
          be turned into one via ``self.submodule(sub)``

        - ``check`` -- boolean (default: ``True``); whether or not to check
          that ``sub`` is a submodule

        - further named arguments, that are passed to the constructor
          of the quotient space

        EXAMPLES::

            sage: A = ZZ^3; V = A.span([[1,2,3], [4,5,6]])
            sage: Q = V.quotient( [V.0 + V.1] ); Q
            Finitely generated module V/W over Integer Ring with invariants (0)
        """
        # Calling is_subspace may be way too slow and repeat work done below.
        # It will be very desirable to somehow do this step better.
        if check and (not isinstance(sub, FreeModule_generic) or not sub.is_submodule(self)):
            try:
                sub = self.submodule(sub)
            except (TypeError, ArithmeticError):
                raise ArithmeticError("sub must be a subspace of self")
        if self.base_ring() == sage.rings.integer_ring.ZZ:
            from sage.modules.fg_pid.fgp_module import FGP_Module
            return FGP_Module(self, sub, check=False, **kwds)

        raise NotImplementedError("quotients of modules over rings other than fields or ZZ is not fully implemented")


class FreeModule_generic_field(FreeModule_generic_pid):
    """
    Base class for all free modules over fields.
    """
    def __init__(self, base_field, dimension, degree, sparse=False, category=None):
        """
        Create a vector space over a field.

        EXAMPLES::

            sage: FreeModule(QQ, 2)
            Vector space of dimension 2 over Rational Field
            sage: FreeModule(FiniteField(2), 7)
            Vector space of dimension 7 over Finite Field of size 2

        We test that objects of this type are initialised correctly;
        see :issue:`11166` (the failing ``repr`` is fine because this
        is an abstract base class)::

            sage: from sage.modules.free_module import FreeModule_generic_field
            sage: FreeModule_generic_field(QQ, 5, 5)
            <repr(<sage.modules.free_module.FreeModule_generic_field_with_category at 0x...>) failed: NotImplementedError>
        """
        if base_field not in Fields():
            raise TypeError("The base_field (=%s) must be a field" % base_field)
        super().__init__(base_field, dimension, degree, sparse=sparse, category=category)

    def _Hom_(self, Y, category):
        r"""
        Return a homspace whose morphisms have this vector space as domain.

        This is called by the general methods such as
        :meth:`sage.structure.parent.Parent.Hom`.

        INPUT:

        - ``Y`` -- a free module (or vector space) that will
          be the codomain of the morphisms in returned homspace
        - ``category`` -- the category for the homspace

        OUTPUT:

        If ``Y`` is a free module over a field, in other words, a vector space,
        then this returns a space of homomorphisms between vector spaces,
        in other words a space of linear transformations.

        If ``Y`` is a free module that is not a vector space, then
        the returned space contains homomorphisms between free modules.

        EXAMPLES::

            sage: V = QQ^2
            sage: W = QQ^3
            sage: H = V._Hom_(W, category=None)
            sage: type(H)
            <class 'sage.modules.vector_space_homspace.VectorSpaceHomspace_with_category'>
            sage: H
            Set of Morphisms (Linear Transformations)
             from Vector space of dimension 2 over Rational Field
               to Vector space of dimension 3 over Rational Field

            sage: V = QQ^2
            sage: W = ZZ^3
            sage: H = V._Hom_(W, category=None)
            sage: type(H)
            <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
            sage: H
            Set of Morphisms from Vector space of dimension 2 over Rational
             Field to Ambient free module of rank 3 over the principal ideal
             domain Integer Ring in Category of commutative additive groups
        """
        if Y.base_ring().is_field():
            from sage.modules import vector_space_homspace
            return vector_space_homspace.VectorSpaceHomspace(self, Y, category)
        from sage.modules import free_module_homspace
        return free_module_homspace.FreeModuleHomspace(self, Y, category)

    def scale(self, other):
        """
        Return the product of ``self`` by the number other, which is the module
        spanned by ``other`` times each basis vector. Since ``self`` is a vector
        space this product equals ``self`` if ``other`` is nonzero, and is the zero
        vector space if ``other`` is 0.

        EXAMPLES::

            sage: V = QQ^4
            sage: V.scale(5)
            Vector space of dimension 4 over Rational Field
            sage: V.scale(0)
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []

        ::

            sage: W = V.span([[1,1,1,1]])
            sage: W.scale(2)
            Vector space of degree 4 and dimension 1 over Rational Field
            Basis matrix:
            [1 1 1 1]
            sage: W.scale(0)
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []

        ::

            sage: V = QQ^4; V
            Vector space of dimension 4 over Rational Field
            sage: V.scale(3)
            Vector space of dimension 4 over Rational Field
            sage: V.scale(0)
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        if other == 0:
            return self.zero_submodule()
        return self

    def __add__(self, other):
        """
        Return the sum of ``self`` and ``other``.

        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: V0 = V.span([V.gen(0)])
            sage: V2 = V.span([V.gen(2)])
            sage: V0 + V2
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 0 1]
            sage: QQ^3 + 0
            Vector space of dimension 3 over Rational Field
        """
        if not isinstance(other, FreeModule_generic_field):
            if other == 0:
                return self
            raise TypeError("other must be a Vector Space")
        V = self.ambient_vector_space()
        if V != other.ambient_vector_space():
            raise ArithmeticError("self and other must have the same ambient space")
        return V.span(self.basis() + other.basis())

    def echelonized_basis_matrix(self):
        """
        Return basis matrix for ``self`` in row echelon form.

        EXAMPLES::

            sage: V = FreeModule(QQ, 3).span_of_basis([[1,2,3],[4,5,6]])
            sage: V.basis_matrix()
            [1 2 3]
            [4 5 6]
            sage: V.echelonized_basis_matrix()
            [ 1  0 -1]
            [ 0  1  2]
        """
        return self.basis_matrix().echelon_form()

    def intersection(self, other):
        """
        Return the intersection of ``self`` and ``other``, which must be
        `R`-submodules of a common ambient vector space.

        EXAMPLES::

            sage: V  = VectorSpace(QQ,3)
            sage: W1 = V.submodule([V.gen(0), V.gen(0) + V.gen(1)])
            sage: W2 = V.submodule([V.gen(1), V.gen(2)])
            sage: W1.intersection(W2)
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
            sage: W2.intersection(W1)
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
            sage: V.intersection(W1)
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            sage: W1.intersection(V)
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            sage: Z = V.submodule([])
            sage: W1.intersection(Z)
            Vector space of degree 3 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        if not isinstance(other, FreeModule_generic):
            raise TypeError("other must be a free module")

        if self.ambient_vector_space() != other.ambient_vector_space():
            raise ArithmeticError("self and other must have the same ambient space.")

        if self.rank() == 0 or other.rank() == 0:
            if self.base_ring().is_field():
                return other.zero_submodule()
            else:
                return self.zero_submodule()

        if self.base_ring() != other.base_ring():
            # Now other is over a ring R whose fraction field K is the base field of V = self.
            # We compute the intersection using the following algorithm:
            # 1. By explicitly computing the nullspace of the matrix whose rows
            #    are a basis for self, we obtain the matrix over a linear map
            #         phi:  K^n ----> W
            #    with kernel equal to V = self.
            # 2. Compute the kernel over R of Phi restricted to other.  Do this
            #    by clearing denominators, computing the kernel of a matrix with
            #    entries in R, then restoring denominators to the answer.
            K = self.base_ring()
            B = self.basis_matrix().transpose()
            W = B.kernel()
            phi = W.basis_matrix().transpose()

            # To restrict phi to other, we multiply the basis matrix for other
            # by phi, thus computing the image of each basis vector.
            X = other.basis_matrix()
            psi = X * phi

            # Now psi is a matrix that defines an R-module morphism from other to some
            # R-module, whose kernel defines the long sought for intersection of self and other.
            L = psi.integer_kernel()

            # Finally the kernel of the intersection has basis the linear combinations of
            # the basis of other given by a basis for L.
            G = L.basis_matrix() * other.basis_matrix()
            return other.span(G.rows())

        # dispense with the three easy cases
        if self == self.ambient_vector_space():
            return other
        elif other == other.ambient_vector_space():
            return self
        elif self.dimension() == 0 or other.dimension() == 0:
            return self.zero_submodule()

        # standard algorithm for computing intersection of general subspaces
        if self.dimension() <= other.dimension():
            V1 = self
            V2 = other
        else:
            V1 = other
            V2 = self
        A1 = V1.basis_matrix()
        A2 = V2.basis_matrix()
        S = A1.stack(A2)
        K = S.kernel()
        n = int(V1.dimension())
        B = [A1.linear_combination_of_rows(v.list()[:n]) for v in K.basis()]
        return self.ambient_vector_space().submodule(B, check=False)

    def is_subspace(self, other):
        """
        ``True`` if this vector space is a subspace of ``other``.

        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: W = V.subspace([V.gen(0), V.gen(0) + V.gen(1)])
            sage: W2 = V.subspace([V.gen(1)])
            sage: W.is_subspace(V)
            True
            sage: W2.is_subspace(V)
            True
            sage: W.is_subspace(W2)
            False
            sage: W2.is_subspace(W)
            True
        """
        return self.is_submodule(other)

    def span_of_basis(self, basis, base_ring=None, check=True, already_echelonized=False):
        r"""
        Return the free K-module with the given basis, where K is the base
        field of ``self`` or user specified base_ring.

        Note that this span is a subspace of the ambient vector space, but
        need not be a subspace of ``self``.

        INPUT:

        - ``basis`` -- list of vectors

        - ``check`` -- boolean (default: ``True``); whether or not to
          coerce entries of gens into base field

        - ``already_echelonized`` -- boolean (default: ``False``);
          set this if you know the gens are already in echelon form

        EXAMPLES::

            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace([[2,3,4]]); W
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 5 2]
            sage: W.span_of_basis([[2,2,2], [3,3,0]])
            Vector space of degree 3 and dimension 2 over Finite Field of size 7
            User basis matrix:
            [2 2 2]
            [3 3 0]

        The basis vectors must be linearly independent or a
        :exc:`ValueError` exception is raised::

            sage: W.span_of_basis([[2,2,2], [3,3,3]])
            Traceback (most recent call last):
            ...
            ValueError: The given basis vectors must be linearly independent.
        """
        if isinstance(basis, FreeModule_generic):
            basis = basis.gens()
        if base_ring is None:
            return FreeModule_submodule_with_basis_field(
                self.ambient_module(), basis=basis, check=check, already_echelonized=already_echelonized)
        else:
            try:
                M = self.change_ring(base_ring)
            except TypeError:
                raise ValueError("Argument base_ring (= %s) is not compatible with the base field (= %s)." % (
                    base_ring, self.base_field()))
            try:
                return M.span_of_basis(basis)
            except TypeError:
                raise ValueError("Argument basis (= %s) is not compatible with base_ring (= %s)." % (basis, base_ring))

    def subspace(self, gens, check=True, already_echelonized=False):
        """
        Return the subspace of ``self`` spanned by the elements of gens.

        INPUT:

        - ``gens`` -- list of vectors

        - ``check`` -- boolean (default: ``True``); verify that gens
          are all in ``self``

        - ``already_echelonized`` -- boolean (default: ``False``); set
          to ``True`` if you know the gens are in Echelon form

        EXAMPLES:

        First we create a 1-dimensional vector subspace of an
        ambient `3`-dimensional space over the finite field of
        order `7`::

            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace([[2,3,4]]); W
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 5 2]

        Next we create an invalid subspace, but it's allowed since
        ``check=False``. This is just equivalent to computing
        the span of the element::

            sage: W.subspace([[1,1,0]], check=False)
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 1 0]

        With ``check=True`` (the default) the mistake is correctly
        detected and reported with an :exc:`ArithmeticError` exception::

            sage: W.subspace([[1,1,0]], check=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: argument gens (= [[1, 1, 0]]) does not generate a submodule of self
        """
        return self.submodule(gens, check=check, already_echelonized=already_echelonized)

    def subspaces(self, dim):
        """
        Iterate over all subspaces of dimension dim.

        INPUT:

        - ``dim`` -- integer; dimension of subspaces to be generated

        EXAMPLES::

            sage: V = VectorSpace(GF(3), 5)
            sage: len(list(V.subspaces(0)))
            1
            sage: len(list(V.subspaces(1)))
            121
            sage: len(list(V.subspaces(2)))
            1210
            sage: len(list(V.subspaces(3)))
            1210
            sage: len(list(V.subspaces(4)))
            121
            sage: len(list(V.subspaces(5)))
            1

        ::

            sage: V = VectorSpace(GF(3), 5)
            sage: V = V.subspace([V([1,1,0,0,0]), V([0,0,1,1,0])])
            sage: list(V.subspaces(1))
            [Vector space of degree 5 and dimension 1 over Finite Field of size 3
            Basis matrix:
            [1 1 0 0 0],
             Vector space of degree 5 and dimension 1 over Finite Field of size 3
            Basis matrix:
            [1 1 1 1 0],
             Vector space of degree 5 and dimension 1 over Finite Field of size 3
            Basis matrix:
            [1 1 2 2 0],
             Vector space of degree 5 and dimension 1 over Finite Field of size 3
            Basis matrix:
            [0 0 1 1 0]]
        """
        if not self.base_ring().is_finite():
            raise RuntimeError("Base ring must be finite.")
        b = self.basis_matrix()
        from sage.matrix.echelon_matrix import reduced_echelon_matrix_iterator
        for m in reduced_echelon_matrix_iterator(self.base_ring(), dim, self.dimension(), self.is_sparse(), copy=False):
            yield self.subspace((m * b).rows())

    def subspace_with_basis(self, gens, check=True, already_echelonized=False):
        """
        Same as ``self.submodule_with_basis(...)``.

        EXAMPLES:

        We create a subspace with a user-defined basis.

        ::

            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace_with_basis([[2,2,2], [1,2,3]]); W
            Vector space of degree 3 and dimension 2 over Finite Field of size 7
            User basis matrix:
            [2 2 2]
            [1 2 3]

        We then create a subspace of the subspace with user-defined basis.

        ::

            sage: W1 = W.subspace_with_basis([[3,4,5]]); W1
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            User basis matrix:
            [3 4 5]

        Notice how the basis for the same subspace is different if we
        merely use the ``subspace`` command.

        ::

            sage: W2 = W.subspace([[3,4,5]]); W2
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 6 4]

        Nonetheless the two subspaces are equal (as mathematical objects)::

            sage: W1 == W2
            True
        """
        return self.submodule_with_basis(gens, check=check, already_echelonized=already_echelonized)

    def complement(self):
        r"""
        Return the complement of ``self`` in the
        :meth:`~sage.modules.free_module.FreeModule_ambient_field.ambient_vector_space`.

        EXAMPLES::

            sage: V = QQ^3
            sage: V.complement()
            Vector space of degree 3 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: V == V.complement().complement()
            True
            sage: W = V.span([[1, 0, 1]])
            sage: X = W.complement(); X
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  0]
            sage: X.complement() == W
            True
            sage: X + W == V
            True

        Even though we construct a subspace of a subspace, the
        orthogonal complement is still done in the ambient vector
        space `\QQ^3`::

            sage: V = QQ^3
            sage: W = V.subspace_with_basis([[1,0,1],[-1,1,0]])
            sage: X = W.subspace_with_basis([[1,0,1]])
            sage: X.complement()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  0]

        All these complements are only done with respect to the inner
        product in the usual basis.  Over finite fields, this means
        we can get complements which are only isomorphic to a vector
        space decomposition complement. ::

            sage: F2 = GF(2, 'x')
            sage: V = F2^6
            sage: W = V.span([[1,1,0,0,0,0]]); W
            Vector space of degree 6 and dimension 1 over Finite Field of size 2
            Basis matrix:
            [1 1 0 0 0 0]
            sage: W.complement()
            Vector space of degree 6 and dimension 5 over Finite Field of size 2
            Basis matrix:
            [1 1 0 0 0 0]
            [0 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
            sage: W.intersection(W.complement())
            Vector space of degree 6 and dimension 1 over Finite Field of size 2
            Basis matrix:
            [1 1 0 0 0 0]
        """
        # Check simple cases
        if self.dimension() == 0:
            return self.ambient_vector_space()
        if self.dimension() == self.ambient_vector_space().dimension():
            return self.submodule([])
        return self.basis_matrix().right_kernel()

    def vector_space(self, base_field=None):
        """
        Return the vector space associated to ``self``. Since ``self`` is a vector
        space this function simply returns ``self``, unless the base field is
        different.

        EXAMPLES::

            sage: V = span([[1,2,3]],QQ); V
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
            sage: V.vector_space()
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
        """
        if base_field is None:
            return self
        return self.change_ring(base_field)

    def zero_submodule(self):
        """
        Return the zero submodule of ``self``.

        EXAMPLES::

            sage: (QQ^4).zero_submodule()
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        return self.zero_subspace()

    def zero_subspace(self):
        """
        Return the zero subspace of ``self``.

        EXAMPLES::

            sage: (QQ^4).zero_subspace()
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        return self.submodule([], check=False, already_echelonized=True)

    def linear_dependence(self, vectors, zeros='left', check=True):
        r"""
        Return a list of vectors giving relations of linear dependence for the
        input list of vectors.

        Can be used to check linear independence of a set of vectors.

        INPUT:

        - ``vectors`` -- list of vectors; all from the same vector
          space

        - ``zeros`` -- (default: ``'left'``) ``'left'`` or ``'right'``
          as a general preference for where zeros are located in the
          returned coefficients

        - ``check`` -- (default: ``True``) if ``True`` each item in
          the list ``vectors`` is checked for membership in ``self``.
          Set to ``False`` if you can be certain the vectors come from
          the vector space.

        OUTPUT:

        Returns a list of vectors.  The scalar entries of each vector provide
        the coefficients for a linear combination of the input vectors that
        will equal the zero vector in ``self``.  Furthermore, the returned list
        is linearly independent in the vector space over the same base field
        with degree equal to the length of the list ``vectors``.

        The linear independence of ``vectors`` is equivalent to the returned
        list being empty, so this provides a test - see the examples below.

        The returned vectors are always independent, and with ``zeros`` set to
        ``'left'`` they have 1s in their first nonzero entries and a qualitative
        disposition to having zeros in the low-index entries.  With ``zeros`` set
        to ``'right'`` the situation is reversed with a qualitative disposition
        for zeros in the high-index entries.

        If the vectors in ``vectors`` are made the rows of a matrix `V` and
        the returned vectors are made the rows of a matrix `R`, then the
        matrix product `RV` is a zero matrix of the proper size.  And
        `R` is a matrix of full rank.  This routine uses kernels of
        matrices to compute these relations of linear dependence,
        but handles all the conversions between sets of vectors
        and matrices. If speed is important, consider working with
        the appropriate matrices and kernels instead.

        EXAMPLES:

        We begin with two linearly independent vectors, and add three
        non-trivial linear combinations to the set.  We illustrate
        both types of output and check a selected relation of linear
        dependence. ::

            sage: v1 = vector(QQ, [2, 1, -4, 3])
            sage: v2 = vector(QQ, [1, 5, 2, -2])
            sage: V = QQ^4
            sage: V.linear_dependence([v1,v2])
            []

            sage: v3 = v1 + v2
            sage: v4 = 3*v1 - 4*v2
            sage: v5 = -v1 + 2*v2
            sage: L = [v1, v2, v3, v4, v5]

            sage: relations = V.linear_dependence(L, zeros='left')
            sage: relations
            [(1, 0, 0, -1, -2), (0, 1, 0, -1/2, -3/2), (0, 0, 1, -3/2, -7/2)]
            sage: v2 + (-1/2)*v4 + (-3/2)*v5
            (0, 0, 0, 0)

            sage: relations = V.linear_dependence(L, zeros='right')
            sage: relations
            [(-1, -1, 1, 0, 0), (-3, 4, 0, 1, 0), (1, -2, 0, 0, 1)]
            sage: z = sum([relations[2][i]*L[i] for i in range(len(L))])
            sage: z == zero_vector(QQ, 4)
            True

        A linearly independent set returns an empty list,
        a result that can be tested. ::

            sage: v1 = vector(QQ, [0,1,-3])
            sage: v2 = vector(QQ, [4,1,0])
            sage: V = QQ^3
            sage: relations = V.linear_dependence([v1, v2]); relations
            []
            sage: relations == []
            True

        Exact results result from exact fields. We start with three
        linearly independent vectors and add in two linear combinations
        to make a linearly dependent set of five vectors. ::

            sage: F = FiniteField(17)
            sage: v1 = vector(F, [1, 2, 3, 4, 5])
            sage: v2 = vector(F, [2, 4, 8, 16, 15])
            sage: v3 = vector(F, [1, 0, 0, 0, 1])
            sage: (F^5).linear_dependence([v1, v2, v3]) == []
            True
            sage: L = [v1, v2, v3, 2*v1+v2, 3*v2+6*v3]
            sage: (F^5).linear_dependence(L)
            [(1, 0, 16, 8, 3), (0, 1, 2, 0, 11)]
            sage: v1 + 16*v3 + 8*(2*v1+v2) + 3*(3*v2+6*v3)
            (0, 0, 0, 0, 0)
            sage: v2 + 2*v3 + 11*(3*v2+6*v3)
            (0, 0, 0, 0, 0)
            sage: (F^5).linear_dependence(L, zeros='right')
            [(15, 16, 0, 1, 0), (0, 14, 11, 0, 1)]

        TESTS:

        With ``check=True`` (the default) a mismatch between vectors
        and the vector space is caught. ::

            sage: v1 = vector(RR, [1,2,3])
            sage: v2 = vector(RR, [1,2,3,4])
            sage: (RR^3).linear_dependence([v1,v2], check=True)
            Traceback (most recent call last):
            ...
            ValueError: vector (1.00000000000000, 2.00000000000000, 3.00000000000000, 4.00000000000000)
            is not an element of Vector space of dimension 3 over Real Field with 53 bits of precision

        The ``zeros`` keyword is checked. ::

            sage: (QQ^3).linear_dependence([vector(QQ,[1,2,3])], zeros='bogus')
            Traceback (most recent call last):
            ...
            ValueError: 'zeros' keyword must be 'left' or 'right', not 'bogus'

        An empty input set is linearly independent, vacuously. ::

            sage: (QQ^3).linear_dependence([]) == []
            True
        """
        if check:
            for v in vectors:
                if v not in self:
                    raise ValueError('vector %s is not an element of %s' % (v, self))
        if zeros == 'left':
            basis = 'echelon'
        elif zeros == 'right':
            basis = 'pivot'
        else:
            raise ValueError("'zeros' keyword must be 'left' or 'right', not '%s'" % zeros)
        import sage.matrix.constructor
        A = sage.matrix.constructor.matrix(vectors)  # as rows, so get left kernel
        return A.left_kernel(basis=basis).basis()

    def __truediv__(self, sub):
        """
        Return the quotient of ``self`` by the given subspace sub.

        EXAMPLES::

            sage: V = RDF^3; W = V.span([[1,0,-1], [1,-1,0]])
            sage: Q = V/W; Q
            Vector space quotient V/W of dimension 1 over Real Double Field where
            V: Vector space of dimension 3 over Real Double Field
            W: Vector space of degree 3 and dimension 2 over Real Double Field
            Basis matrix:
            [ 1.0  0.0 -1.0]
            [ 0.0  1.0 -1.0]
            sage: type(Q)
            <class 'sage.modules.quotient_module.FreeModule_ambient_field_quotient_with_category'>
            sage: V([1,2,3])
            (1.0, 2.0, 3.0)
            sage: Q == V.quotient(W)
            True
            sage: Q(W.0)
            (0.0)
        """
        return self.quotient(sub, check=True)

    def quotient_module(self, sub, check=True):
        """
        Return the quotient of ``self`` by the given subspace sub.

        INPUT:

        - ``sub`` -- a submodule of ``self``, or something that can
          be turned into one via ``self.submodule(sub)``

        - ``check`` -- boolean (default: ``True``); whether or not to check
          that ``sub`` is a submodule

        EXAMPLES::

            sage: A = QQ^3; V = A.span([[1,2,3], [4,5,6]])
            sage: Q = V.quotient( [V.0 + V.1] ); Q
            Vector space quotient V/W of dimension 1 over Rational Field where
            V: Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 1 1]
            sage: Q(V.0 + V.1)
            (0)

        We illustrate that the base rings must be the same::

            sage: (QQ^2)/(ZZ^2)
            Traceback (most recent call last):
            ...
            ValueError: base rings must be the same
        """
        # Calling is_submodule may be way too slow and repeat work done below.
        # It will be very desirable to somehow do this step better.
        if isinstance(sub, FreeModule_generic) and self.base_ring() != sub.base_ring():
            raise ValueError("base rings must be the same")
        if check and (not isinstance(sub, FreeModule_generic) or not sub.is_subspace(self)):
            try:
                sub = self.subspace(sub)
            except (TypeError, ArithmeticError):
                raise ArithmeticError("sub must be a subspace of self")
        A, L = self.__quotient_matrices(sub)
        from sage.modules import quotient_module
        return quotient_module.FreeModule_ambient_field_quotient(self, sub, A, L)

    def __quotient_matrices(self, sub):
        r"""
        This internal function is used by
        ``self.quotient(...)``.

        EXAMPLES::

            sage: V = QQ^3; W = V.span([[1,0,-1], [1,-1,0]])
            sage: A, L = V._FreeModule_generic_field__quotient_matrices(W)
            sage: A
            [1]
            [1]
            [1]
            sage: L
            [1 0 0]

        The quotient and lift maps are used to compute in the quotient and
        to lift::

            sage: Q = V/W
            sage: Q(W.0)
            (0)
            sage: Q.lift_map()(Q.0)
            (1, 0, 0)
            sage: Q(Q.lift_map()(Q.0))
            (1)

        An example in characteristic 5::

            sage: A = GF(5)^2; B = A.span([[1,3]]); A / B
            Vector space quotient V/W of dimension 1 over Finite Field of size 5 where
            V: Vector space of dimension 2 over Finite Field of size 5
            W: Vector space of degree 2 and dimension 1 over Finite Field of size 5
            Basis matrix:
            [1 3]
        """
        # 2. Find a basis C for a another submodule of self, so that
        #    B + C is a basis for self.
        # 3. Then the quotient map is:
        #     x |---> 'write in terms of basis for C and take the last m = #C-#B components.
        # 4. And a section of this map is:
        #     x |---> corresponding linear combination of entries of last m entries
        #    of the basis C.

        # Step 1: Find bases for spaces
        B = sub.basis_matrix()
        S = self.basis_matrix()

        n = self.dimension()
        m = n - sub.dimension()

        # Step 2: Extend basis B to a basis for self.
        # We do this by simply finding the pivot rows of the matrix
        # whose rows are a basis for sub concatenated with a basis for
        # self.
        C = B.stack(S).transpose()
        A = C.matrix_from_columns(C.pivots()).transpose()

        # Step 3: Compute quotient map
        # The quotient map is given by writing in terms of the above basis,
        # then taking the last #C columns

        # Compute the matrix D "change of basis from S to A"
        # that writes each element of the basis
        # for ``self`` in terms of the basis of rows of A, i.e.,
        # want to find D such that
        #                D * A = S
        # where D is a square n x n matrix.
        # Our algorithm is to note that D is determined if we just
        # replace both A and S by the submatrix got from their pivot
        # columns.
        P = A.pivots()
        AA = A.matrix_from_columns(P)
        SS = S.matrix_from_columns(P)
        D = SS * AA**(-1)

        # Compute the image of each basis vector for ``self`` under the
        # map "write an element of ``self`` in terms of the basis A" then
        # take the last n-m components.
        Q = D.matrix_from_columns(range(n - m, n))

        # Step 4. Section map
        # The lifting or section map
        Dinv = D**(-1)
        L = Dinv.matrix_from_rows(range(n - m, n))

        return Q, L

    def quotient_abstract(self, sub, check=True, **kwds):
        r"""
        Return an ambient free module isomorphic to the quotient space
        of ``self`` modulo ``sub``, together with maps from ``self`` to
        the quotient, and a lifting map in the other direction.

        Use ``self.quotient(sub)`` to obtain the quotient
        module as an object equipped with natural maps in both directions,
        and a canonical coercion.

        INPUT:

        - ``sub`` -- a submodule of ``self`` or something that can
          be turned into one via ``self.submodule(sub)``

        - ``check`` -- boolean (default: ``True``); whether or not to check
          that sub is a submodule

        - further named arguments, that are currently ignored.

        OUTPUT:

        - ``U``; the quotient as an abstract *ambient* free module

        - ``pi``; projection map to the quotient

        - ``lift``; lifting map back from quotient

        EXAMPLES::

            sage: V = GF(19)^3
            sage: W = V.span_of_basis([[1,2,3], [1,0,1]])
            sage: U, pi, lift = V.quotient_abstract(W)
            sage: pi(V.2)
            (18)
            sage: pi(V.0)
            (1)
            sage: pi(V.0 + V.2)
            (0)

        Another example involving a quotient of one subspace by another::

            sage: A = matrix(QQ, 4,4, [0,1,0,0, 0,0,1,0, 0,0,0,1, 0,0,0,0])
            sage: V = (A^3).kernel()
            sage: W = A.kernel()
            sage: U, pi, lift = V.quotient_abstract(W)
            sage: [pi(v) == 0 for v in W.gens()]
            [True]
            sage: [pi(lift(b)) == b for b in U.basis()]
            [True, True]
        """
        # Calling is_subspace may be way too slow and repeat work done below.
        # It will be very desirable to somehow do this step better.
        if check and (not isinstance(sub, FreeModule_generic) or not sub.is_subspace(self)):
            try:
                sub = self.subspace(sub)
            except (TypeError, ArithmeticError):
                raise ArithmeticError("sub must be a subspace of self")

        A, L = self.__quotient_matrices(sub)
        quomap = self.hom(A)
        quo = quomap.codomain()
        liftmap = quo.Hom(self)(L)

        return quomap.codomain(), quomap, liftmap


###############################################################################
#
# Generic ambient free module R^n for some commutative ring R
#
###############################################################################

class FreeModule_ambient(FreeModule_generic):
    """
    Ambient free module over a commutative ring.
    """
    def __init__(self, base_ring, rank, sparse=False, coordinate_ring=None, category=None):
        """
        The free module of given rank over the given base_ring.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``rank`` -- nonnegative integer

        - ``sparse`` -- boolean (default: ``False``)

        - ``coordinate_ring`` -- a ring containing ``base_ring``
          (default: equal to ``base_ring``)

        EXAMPLES::

            sage: FreeModule(ZZ, 4)
            Ambient free module of rank 4 over the principal ideal domain Integer Ring

        TESTS:

        We check that the creation of a submodule does not trigger
        the construction of a basis of the ambient space. See :issue:`15953`::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: V = VectorSpace(F, 1)
            sage: v = V.random_element()
            sage: _ = V.subspace([v])
            sage: hasattr(V, '_FreeModule_ambient__basis')
            False
            sage: _ = V.basis()
            sage: hasattr(V, '_FreeModule_ambient__basis')
            True
        """
        FreeModule_generic.__init__(self, base_ring, rank=rank,
                                    degree=rank, sparse=sparse,
                                    coordinate_ring=coordinate_ring,
                                    category=category)

    def __hash__(self):
        """
        The hash is obtained from the rank and the base ring.

        .. TODO::

            Make pickling so that the hash is available early enough.

        EXAMPLES::

            sage: V = QQ^7
            sage: hash(V) == hash((V.rank(), V.base_ring()))
            True
        """
        try:
            return hash((self.rank(), self.base_ring()))
        except AttributeError:
            # This is a fallback because sometimes hash is called during object
            # reconstruction (unpickle), and the above fields haven't been
            # filled in yet.
            return 0

    def _coerce_map_from_(self, M):
        """
        Return a coercion map from `M` to ``self``, or ``None``.

        TESTS:

        Make sure :issue:`10513` is fixed (no coercion from a quotient
        vector space to an isomorphic abstract vector space)::

            sage: M = QQ^3 / [[1,2,3]]
            sage: V = QQ^2
            sage: V.coerce_map_from(M)
        """
        from sage.modules.quotient_module import FreeModule_ambient_field_quotient
        from sage.modules.submodule import Submodule_free_ambient

        if isinstance(M, FreeModule_ambient_field_quotient):
            # No forgetful map.
            return None
        if isinstance(M, FreeModule_ambient):
            if (self.base_ring().has_coerce_map_from(M.base_ring()) and
                    self.rank() == M.rank()):
                # We could return M.hom(self.basis(), self), but the
                # complexity of this is quadratic in space and time,
                # since it constructs a matrix.
                return True
        elif isinstance(M, Submodule_free_ambient):
            if (self.base_ring().has_coerce_map_from(M.base_ring()) and
                    self.rank() == M.degree()):
                return True
        return super()._coerce_map_from_(M)

    def _dense_module(self):
        """
        Create a dense module with the same defining data as ``self``.

        N.B. This function is for internal use only! See dense_module for
        use.

        EXAMPLES::

            sage: M = FreeModule(Integers(8),3)
            sage: S = FreeModule(Integers(8),3, sparse=True)
            sage: M is S._dense_module()
            True
        """
        return FreeModule(base_ring=self.base_ring(), rank=self.rank(), sparse=False)

    def _sparse_module(self):
        """
        Create a sparse module with the same defining data as ``self``.

        N.B. This function is for internal use only! See sparse_module for
        use.

        EXAMPLES::

            sage: M = FreeModule(Integers(8),3)
            sage: S = FreeModule(Integers(8),3, sparse=True)
            sage: M._sparse_module() is S
            True
        """
        return FreeModule(base_ring=self.base_ring(), rank=self.rank(), sparse=True)

    def echelonized_basis_matrix(self):
        """
        The echelonized basis matrix of ``self``.

        EXAMPLES::

            sage: V = ZZ^4
            sage: W = V.submodule([ V.gen(i)-V.gen(0) for i in range(1,4) ])
            sage: W.basis_matrix()
            [ 1  0  0 -1]
            [ 0  1  0 -1]
            [ 0  0  1 -1]
            sage: W.echelonized_basis_matrix()
            [ 1  0  0 -1]
            [ 0  1  0 -1]
            [ 0  0  1 -1]
            sage: U = V.submodule_with_basis([ V.gen(i)-V.gen(0) for i in range(1,4) ])
            sage: U.basis_matrix()
            [-1  1  0  0]
            [-1  0  1  0]
            [-1  0  0  1]
            sage: U.echelonized_basis_matrix()
            [ 1  0  0 -1]
            [ 0  1  0 -1]
            [ 0  0  1 -1]
        """
        return self.basis_matrix()

    def _echelon_matrix_richcmp(self, other, op):
        r"""
        Compare the free module ``self`` with ``other``.

        This compares modules by their ambient spaces, then by dimension,
        then in order by their echelon matrices. However, if
        ``other`` is a sub-module or is a quotient module then its
        total comparison method is used instead of generic comparison.

        EXAMPLES:

        We compare rank three free modules over the integers and
        rationals::

            sage: from sage.structure.richcmp import op_LT,op_LE,op_EQ,op_NE,op_GT,op_GE
            sage: (QQ^3)._echelon_matrix_richcmp(CC^3, op_LT)
            True
            sage: (CC^3)._echelon_matrix_richcmp(QQ^3, op_LT)
            False
            sage: (CC^3)._echelon_matrix_richcmp(QQ^3, op_GT)
            True

        ::

            sage: from sage.structure.richcmp import op_LT,op_LE,op_EQ,op_NE,op_GT,op_GE
            sage: Q = QQ; Z = ZZ
            sage: (Q^3)._echelon_matrix_richcmp(Z^3, op_GT)
            True
            sage: (Q^3)._echelon_matrix_richcmp(Z^3, op_LT)
            False
            sage: (Z^3)._echelon_matrix_richcmp(Q^3, op_LT)
            True
            sage: (Z^3)._echelon_matrix_richcmp(Q^3, op_GT)
            False
            sage: (Q^3)._echelon_matrix_richcmp(Z^3, op_EQ)
            False
            sage: (Q^3)._echelon_matrix_richcmp(Q^3, op_EQ)
            True

        Comparison with a sub-module::

            sage: from sage.structure.richcmp import op_LT,op_LE,op_EQ,op_NE,op_GT,op_GE
            sage: V = span([[1,2,3], [5,6,7], [8,9,10]], QQ)
            sage: V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: A = QQ^3
            sage: V._echelon_matrix_richcmp(A, op_LT)
            True
            sage: A._echelon_matrix_richcmp(V, op_LT)
            False

        Comparison with a quotient module (see :issue:`10513`)::

            sage: from sage.structure.richcmp import op_LT,op_LE,op_EQ,op_NE,op_GT,op_GE
            sage: M = QQ^3 / [[1,2,3]]
            sage: V = QQ^2
            sage: V._echelon_matrix_richcmp(M, op_EQ)
            False
            sage: M._echelon_matrix_richcmp(V, op_EQ)
            False
        """
        if self is other:
            return rich_to_bool(op, 0)

        if not isinstance(other, FreeModule_generic):
            return NotImplemented

        from sage.modules.quotient_module import FreeModule_ambient_field_quotient
        if isinstance(other, FreeModule_ambient):
            if (isinstance(other, FreeModule_ambient_field_quotient) or
                    isinstance(self, FreeModule_ambient_field_quotient)):
                return richcmp(self, other, op)

            lx = self.rank()
            rx = other.rank()
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

            lx = self.base_ring()
            rx = other.base_ring()
            if lx == rx:
                # We do not want to create an inner product matrix in memory if
                # self and other use the dot product
                if self._inner_product_is_dot_product() and other._inner_product_is_dot_product():
                    return rich_to_bool(op, 0)
                else:
                    # this only affects free_quadratic_modules
                    lx = self.inner_product_matrix()
                    rx = other.inner_product_matrix()
                    return richcmp(lx, rx, op)

            try:
                if lx.is_subring(rx):
                    return rich_to_bool(op, -1)
                elif rx.is_subring(lx):
                    return rich_to_bool(op, 1)
            except NotImplementedError:
                pass
            return richcmp_not_equal(lx, rx, op)
        else:
            # now other is not ambient or is a quotient;
            # it knows how to do the comparison.
            return other._echelon_matrix_richcmp(self, revop(op))

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: R = ZZ.quo(12)
            sage: M = R^12
            sage: M
            Ambient free module of rank 12 over Ring of integers modulo 12
            sage: print(M._repr_())
            Ambient free module of rank 12 over Ring of integers modulo 12

        The system representation can be overwritten, but leaves _repr_
        unmodified.

        ::

            sage: M.rename('M')
            sage: M
            M
            sage: print(M._repr_())
            Ambient free module of rank 12 over Ring of integers modulo 12

        Sparse modules print this fact.

        ::

            sage: N = FreeModule(R,12,sparse=True)
            sage: N
            Ambient sparse free module of rank 12 over Ring of integers modulo 12

        (Now clean up again.)

        ::

            sage: M.reset_name()
            sage: M
            Ambient free module of rank 12 over Ring of integers modulo 12
        """
        if self.is_sparse():
            return "Ambient sparse free module of rank %s over %s" % (self.rank(), self.base_ring())
        else:
            return "Ambient free module of rank %s over %s" % (self.rank(), self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of this ambient free module.

        EXAMPLES::

            sage: latex(QQ^3)   # indirect doctest
            \Bold{Q}^{3}

        ::

            sage: A = GF(5)^20
            sage: latex(A)      # indirect doctest
            \Bold{F}_{5}^{20}

        ::

            sage: A = PolynomialRing(QQ, 3, 'x')^20
            sage: latex(A)                              # indirect doctest
            (\Bold{Q}[x_{0}, x_{1}, x_{2}])^{20}
        """
        t = "%s" % latex.latex(self.base_ring())
        if t.find(" ") != -1:
            t = "(%s)" % t
        return "%s^{%s}" % (t, self.rank())

    def is_ambient(self):
        """
        Return ``True`` since this module is an ambient
        module.

        EXAMPLES::

            sage: A = QQ^5; A.is_ambient()
            True
            sage: A = (QQ^5).span([[1,2,3,4,5]]); A.is_ambient()
            False
        """
        return True

    def ambient_module(self):
        """
        Return ``self``, since ``self`` is ambient.

        EXAMPLES::

            sage: A = QQ^5; A.ambient_module()
            Vector space of dimension 5 over Rational Field
            sage: A = ZZ^5; A.ambient_module()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
        """
        return self

    def basis(self):
        """
        Return a basis for this ambient free module.

        OUTPUT:

        - ``Sequence`` -- an immutable sequence with universe
          this ambient free module

        EXAMPLES::

            sage: A = ZZ^3; B = A.basis(); B
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
            sage: B.universe()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        try:
            return self.__basis
        except AttributeError:
            ZERO = self(0)
            one = self.coordinate_ring().one()
            w = []
            for n in range(self.rank()):
                v = ZERO.__copy__()
                v.set(n, one)
                w.append(v)
            self.__basis = basis_seq(self, w)
            return self.__basis

    def echelonized_basis(self):
        """
        Return a basis for this ambient free module in echelon form.

        EXAMPLES::

            sage: A = ZZ^3; A.echelonized_basis()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        return self.basis()

    def change_ring(self, R):
        """
        Return the ambient free module over ``R`` of the same rank as ``self``.

        This also preserves the sparsity.

        EXAMPLES::

            sage: A = ZZ^3; A.change_ring(QQ)
            Vector space of dimension 3 over Rational Field
            sage: A = ZZ^3; A.change_ring(GF(5))
            Vector space of dimension 3 over Finite Field of size 5

        For ambient modules any change of rings is defined::

            sage: A = GF(5)**3; A.change_ring(QQ)
            Vector space of dimension 3 over Rational Field

        TESTS:

        Check for :issue:`29630`::

            sage: V = VectorSpace(QQ, 2, sparse=True)
            sage: V.change_ring(RR).is_sparse()
            True
        """
        if self.base_ring() is R:
            return self
        from sage.modules.free_quadratic_module import FreeQuadraticModule_generic
        if isinstance(self, FreeQuadraticModule_generic):
            return FreeModule(R, self.rank(),
                              inner_product_matrix=self.inner_product_matrix(),
                              sparse=self.is_sparse())
        return FreeModule(R, self.rank(), sparse=self.is_sparse())

    def linear_combination_of_basis(self, v):
        """
        Return the linear combination of the basis for ``self`` obtained from
        the elements of the list v.

        INPUT:

        - ``v`` -- list

        EXAMPLES::

            sage: V = span([[1,2,3], [4,5,6]], ZZ)
            sage: V
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            [0 3 6]
            sage: V.linear_combination_of_basis([1,1])
            (1, 5, 9)

        This should raise an error if the resulting element is not in self::

            sage: W = span([[2,4]], ZZ)
            sage: W.linear_combination_of_basis([1/2])
            Traceback (most recent call last):
            ...
            TypeError: element [1, 2] is not in free module
        """
        return self(v)

    def coordinate_vector(self, v, check=True):
        """
        Write `v` in terms of the standard basis for ``self`` and
        return the resulting coefficients in a vector over the fraction
        field of the base ring.

        Returns a vector `c` such that if `B` is the basis for self, then

        .. MATH::

            \\sum c_i B_i = v.

        If `v` is not in self, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: V = Integers(16)^3
            sage: v = V.coordinate_vector([1,5,9]); v
            (1, 5, 9)
            sage: v.parent()
            Ambient free module of rank 3 over Ring of integers modulo 16
        """
        return self(v)

    def echelon_coordinate_vector(self, v, check=True):
        r"""
        Same as ``self.coordinate_vector(v)``, since ``self`` is
        an ambient free module.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also
          verify that `v` is really in ``self``

        OUTPUT: list

        EXAMPLES::

            sage: V = QQ^4
            sage: v = V([-1/2,1/2,-1/2,1/2])
            sage: v
            (-1/2, 1/2, -1/2, 1/2)
            sage: V.coordinate_vector(v)
            (-1/2, 1/2, -1/2, 1/2)
            sage: V.echelon_coordinate_vector(v)
            (-1/2, 1/2, -1/2, 1/2)
            sage: W = V.submodule_with_basis([[1/2,1/2,1/2,1/2],[1,0,1,0]])
            sage: W.coordinate_vector(v)
            (1, -1)
            sage: W.echelon_coordinate_vector(v)
            (-1/2, 1/2)
        """
        return self.coordinate_vector(v, check=check)

    def echelon_coordinates(self, v, check=True):
        """
        Return the coordinate vector of v in terms of the echelon basis
        for ``self``.

        EXAMPLES::

            sage: U = VectorSpace(QQ,3)
            sage: [ U.coordinates(v) for v in U.basis() ]
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            sage: [ U.echelon_coordinates(v) for v in U.basis() ]
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            sage: V = U.submodule([[1,1,0],[0,1,1]])
            sage: V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  1]
            sage: [ V.coordinates(v) for v in V.basis() ]
            [[1, 0], [0, 1]]
            sage: [ V.echelon_coordinates(v) for v in V.basis() ]
            [[1, 0], [0, 1]]
            sage: W = U.submodule_with_basis([[1,1,0],[0,1,1]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [0 1 1]
            sage: [ W.coordinates(v) for v in W.basis() ]
            [[1, 0], [0, 1]]
            sage: [ W.echelon_coordinates(v) for v in W.basis() ]
            [[1, 1], [0, 1]]
        """
        return self.coordinates(v, check=check)

    def random_element(self, prob=1.0, *args, **kwds):
        """
        Return a random element of ``self``.

        INPUT:

        - ``prob`` -- float. Each coefficient will be set to zero with
          probability `1-prob`. Otherwise coefficients will be chosen
          randomly from base ring (and may be zero).

        - ``*args``, ``**kwds`` -- passed on to random_element function of base
          ring

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3)
            sage: M.random_element().parent() is M
            True

        Passes extra positional or keyword arguments through::

            sage: all(i in range(5, 10) for i in M.random_element(1.0, 5, 10))
            True

        ::

            sage: M = FreeModule(ZZ, 16)
            sage: M.random_element().parent() is M
            True

            sage: def add_sample(**kwds):
            ....:     global total, zeros
            ....:     v = M.random_element(**kwds)
            ....:     total += M.rank()
            ....:     zeros += sum(i == 0 for i in v)

            sage: total = 0
            sage: zeros = 0
            sage: add_sample()
            sage: expected = 1/5
            sage: while abs(zeros/total - expected) > 0.01:
            ....:     add_sample()

            sage: total = 0
            sage: zeros = 0
            sage: add_sample(prob=0.3)
            sage: expected = 1/5 * 3/10 + 7/10
            sage: while abs(zeros/total - expected) > 0.01:
            ....:     add_sample(prob=0.3)

            sage: total = 0
            sage: zeros = 0
            sage: add_sample(prob=0.7)
            sage: expected = 1/5 * 7/10 + 3/10
            sage: while abs(zeros/total - expected) > 0.01:
            ....:     add_sample(prob=0.7)
        """
        rand = current_randstate().python_random().random
        R = self.base_ring()
        v = self(0)
        prob = float(prob)
        for i in range(self.rank()):
            if rand() <= prob:
                v[i] = R.random_element(*args, **kwds)
        return v

    def gen(self, i=0):
        """
        Return the `i`-th generator for ``self``.

        Here `i` is between 0 and rank - 1, inclusive.

        INPUT:

        - ``i`` -- integer (default: 0)

        OUTPUT: `i`-th basis vector for ``self``

        EXAMPLES::

            sage: n = 5
            sage: V = QQ^n
            sage: B = [V.gen(i) for i in range(n)]
            sage: B
            [(1, 0, 0, 0, 0),
            (0, 1, 0, 0, 0),
            (0, 0, 1, 0, 0),
            (0, 0, 0, 1, 0),
            (0, 0, 0, 0, 1)]
            sage: V.gens() == tuple(B)
            True

        TESTS::

            sage: (QQ^3).gen(4/3)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert rational 4/3 to an integer

        Check that :issue:`10262` and :issue:`13304` are fixed
        (coercions involving :class:`FreeModule_ambient` used to take
        quadratic time and space in the rank of the module)::

            sage: vector([0]*50000)/1
            (0, 0, 0, ..., 0)
        """
        if i < 0 or i >= self.rank():
            raise ValueError("Generator %s not defined." % i)
        try:
            return self.__basis[i]
        except AttributeError:
            v = self(0)
            v[i] = self.base_ring().one()
            v.set_immutable()
            return v

    def _sympy_(self):
        """
        Return a SymPy ``ProductSet`` corresponding to ``self``.

        EXAMPLES::

            sage: sZZ3 = (ZZ^3)._sympy_(); sZZ3                                         # needs sympy
            ProductSet(Integers, Integers, Integers)
            sage: (1, 2, 3) in sZZ3                                                     # needs sympy
            True
        """
        from sage.interfaces.sympy import sympy_init
        from sympy import ProductSet
        sympy_init()
        return ProductSet(*([self.coordinate_ring()] * self.rank()))


###############################################################################
#
# Ambient free modules over an integral domain
#
###############################################################################

class FreeModule_ambient_domain(FreeModule_generic_domain, FreeModule_ambient):
    """
    Ambient free module over an integral domain.

    EXAMPLES::

        sage: FreeModule(PolynomialRing(GF(5), 'x'), 3)
        Ambient free module of rank 3 over the principal ideal domain
         Univariate Polynomial Ring in x over Finite Field of size 5
    """
    def __init__(self, base_ring, rank, sparse=False, coordinate_ring=None, category=None):
        """
        Create the ambient free module of given rank over the given integral
        domain.

        TESTS::

            sage: A = FreeModule(PolynomialRing(GF(5),'x'), 3)
            sage: TestSuite(A).run()
        """
        FreeModule_ambient.__init__(self, base_ring, rank, sparse, coordinate_ring, category=category)

    def _repr_(self):
        """
        Return the string representation of this free module.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ,'x')
            sage: M = FreeModule(R, 7)
            sage: M
            Ambient free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring
            sage: print(M._repr_())
            Ambient free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring

        The system representation can be overwritten, but leaves ``_repr_`` unmodified.

        ::

            sage: M.rename('M')
            sage: M
            M
            sage: print(M._repr_())
            Ambient free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring

        Sparse modules print this fact.

        ::

            sage: N = FreeModule(R, 7, sparse=True)
            sage: N
            Ambient sparse free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring

        (Now clean up again.)

        ::

            sage: M.reset_name()
            sage: M
            Ambient free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring
        """
        if self.is_sparse():
            return "Ambient sparse free module of rank %s over the integral domain %s" % (
                self.rank(), self.base_ring())
        else:
            return "Ambient free module of rank %s over the integral domain %s" % (
                self.rank(), self.base_ring())

    def ambient_vector_space(self):
        """
        Return the ambient vector space, which is this free module tensored
        with its fraction field.

        EXAMPLES::

            sage: M = ZZ^3
            sage: V = M.ambient_vector_space(); V
            Vector space of dimension 3 over Rational Field

        If an inner product on the module is specified, then this is preserved
        on the ambient vector space.

        ::

            sage: N = FreeModule(ZZ,4,inner_product_matrix=1)
            sage: U = N.ambient_vector_space()
            sage: U
            Ambient quadratic space of dimension 4 over Rational Field
            Inner product matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: P = N.submodule_with_basis([[1,-1,0,0],[0,1,-1,0],[0,0,1,-1]])
            sage: P.gram_matrix()
            [ 2 -1  0]
            [-1  2 -1]
            [ 0 -1  2]
            sage: U == N.ambient_vector_space()
            True
            sage: U == V
            False
        """
        try:
            return self.__ambient_vector_space
        except AttributeError:
            self.__ambient_vector_space = FreeModule(self.base_field(), self.rank(), sparse=self.is_sparse())
            return self.__ambient_vector_space

    def coordinate_vector(self, v, check=True):
        """
        Write `v` in terms of the standard basis for ``self`` and
        return the resulting coefficients in a vector over the fraction
        field of the base ring.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also verify
          that `v` is really in ``self``

        OUTPUT: list

        The output is a vector `c` such that if `B` is the basis for ``self``,
        then

        .. MATH::

            \\sum c_i B_i = v.

        If `v` is not in self, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: V = ZZ^3
            sage: v = V.coordinate_vector([1,5,9]); v
            (1, 5, 9)
            sage: v.parent()
            Vector space of dimension 3 over Rational Field
        """
        # Calling the element constructor directly, since the
        # usual call method indirectly relies on coordinate_vector,
        # hence, an infinite recursion would result
        try:
            out = self.ambient_vector_space()._element_constructor_(v)
        except TypeError:
            raise ArithmeticError("Error transforming the given vector into the ambient vector space")
        if check and out not in self:
            raise ArithmeticError("The given vector does not belong to this free module")
        return out

    def vector_space(self, base_field=None):
        """
        Return the vector space obtained from ``self`` by tensoring with the
        fraction field of the base ring and extending to the field.

        EXAMPLES::

            sage: M = ZZ^3;  M.vector_space()
            Vector space of dimension 3 over Rational Field
        """
        if base_field is None:
            R = self.base_ring()
            return self.change_ring(R.fraction_field())
        else:
            return self.change_ring(base_field)


###############################################################################
#
# Ambient free modules over a principal ideal domain
#
###############################################################################

class FreeModule_ambient_pid(FreeModule_generic_pid, FreeModule_ambient_domain):
    """
    Ambient free module over a principal ideal domain.
    """
    def __init__(self, base_ring, rank, sparse=False, coordinate_ring=None, category=None):
        """
        Create the ambient free module of given rank over the given
        principal ideal domain.

        INPUT:

        - ``base_ring`` -- a principal ideal domain

        - ``rank`` -- nonnegative integer

        - ``sparse`` -- boolean (default: ``False``)

        - ``coordinate_ring`` -- a ring containing ``base_ring``
          (default: equal to ``base_ring``)

        EXAMPLES::

            sage: ZZ^3
            Ambient free module of rank 3 over the principal ideal domain Integer Ring

        We create the same module with coordinates in ``QQ``::

            sage: from sage.modules.free_module import FreeModule_ambient_pid
            sage: M = FreeModule_ambient_pid(ZZ, 3, coordinate_ring=QQ)
            sage: M
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: v = M.basis()[0]; v
            (1, 0, 0)
            sage: type(v)
            <class 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        """
        FreeModule_ambient_domain.__init__(self, base_ring=base_ring,
                                           rank=rank, sparse=sparse,
                                           coordinate_ring=coordinate_ring,
                                           category=category)

    def _repr_(self) -> str:
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: M = FreeModule(ZZ,7)
            sage: M
            Ambient free module of rank 7 over the principal ideal domain Integer Ring
            sage: print(M._repr_())
            Ambient free module of rank 7 over the principal ideal domain Integer Ring

        The system representation can be overwritten, but leaves _repr_
        unmodified.

        ::

            sage: M.rename('M')
            sage: M
            M
            sage: print(M._repr_())
            Ambient free module of rank 7 over the principal ideal domain Integer Ring

        Sparse modules print this fact.

        ::

            sage: N = FreeModule(ZZ,7,sparse=True)
            sage: N
            Ambient sparse free module of rank 7 over the principal ideal domain Integer Ring

        (Now clean up again.)

        ::

            sage: M.reset_name()
            sage: M
            Ambient free module of rank 7 over the principal ideal domain Integer Ring
        """
        if self.is_sparse():
            return "Ambient sparse free module of rank %s over the principal ideal domain %s" % (
                self.rank(), self.base_ring())
        else:
            return "Ambient free module of rank %s over the principal ideal domain %s" % (
                self.rank(), self.base_ring())


###############################################################################
#
# Ambient free modules over a field (vector spaces)
#
###############################################################################

class FreeModule_ambient_field(FreeModule_generic_field, FreeModule_ambient_pid):
    def __init__(self, base_field, dimension, sparse=False, category=None):
        """
        Create the ambient vector space of given dimension over the given
        field.

        INPUT:

        - ``base_field`` -- a field

        - ``dimension`` -- nonnegative integer

        - ``sparse`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: QQ^3
            Vector space of dimension 3 over Rational Field
        """
        FreeModule_ambient_pid.__init__(self, base_field, dimension, sparse=sparse, category=category)

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: V = FreeModule(QQ,7)
            sage: V
            Vector space of dimension 7 over Rational Field
            sage: print(V._repr_())
            Vector space of dimension 7 over Rational Field

        The system representation can be overwritten, but leaves _repr_
        unmodified.

        ::

            sage: V.rename('V')
            sage: V
            V
            sage: print(V._repr_())
            Vector space of dimension 7 over Rational Field

        Sparse modules print this fact.

        ::

            sage: U = FreeModule(QQ,7,sparse=True)
            sage: U
            Sparse vector space of dimension 7 over Rational Field

        (Now clean up again.)

        ::

            sage: V.reset_name()
            sage: V
            Vector space of dimension 7 over Rational Field
        """
        if self.is_sparse():
            return "Sparse vector space of dimension %s over %s" % (self.dimension(), self.base_ring())
        else:
            return "Vector space of dimension %s over %s" % (self.dimension(), self.base_ring())

    def ambient_vector_space(self):
        """
        Return ``self`` as the ambient vector space.

        EXAMPLES::

            sage: M = QQ^3
            sage: M.ambient_vector_space()
            Vector space of dimension 3 over Rational Field
        """
        return self

    def base_field(self):
        """
        Return the base field of this vector space.

        EXAMPLES::

            sage: M = QQ^3
            sage: M.base_field()
            Rational Field
        """
        return self.base_ring()

    def _element_constructor_(self, e, *args, **kwds):
        """
        Create an element of this vector space.

        EXAMPLES::

            sage: k.<a> = GF(3^4)                                                       # needs sage.rings.finite_rings
            sage: VS = k.vector_space(map=False)                                        # needs sage.rings.finite_rings
            sage: VS(a)                                                                 # needs sage.rings.finite_rings
            (0, 1, 0, 0)
        """
        try:
            k = e.parent()
            if isinstance(k, FiniteField) and k.base_ring() == self.base_ring() and k.degree() == self.degree():
                return self(e._vector_())
        except AttributeError:
            pass
        return FreeModule_generic_field._element_constructor_(self, e, *args, **kwds)


class RealDoubleVectorSpace_class(FreeModule_ambient_field):
    def __init__(self, n):
        FreeModule_ambient_field.__init__(self, sage.rings.real_double.RDF, n)

    def coordinates(self, v):
        return v


class ComplexDoubleVectorSpace_class(FreeModule_ambient_field):
    def __init__(self, n):
        FreeModule_ambient_field.__init__(self, sage.rings.complex_double.CDF, n)

    def coordinates(self, v):
        return v


###############################################################################
#
# R-Submodule of K^n where K is the fraction field of a principal ideal domain R
#
###############################################################################

class FreeModule_submodule_with_basis_pid(FreeModule_generic_pid):
    r"""
    Construct a submodule of a free module over PID with a distinguished basis.

    INPUT:

    - ``ambient`` -- ambient free module over a principal ideal domain `R`,
      i.e. `R^n`

    - ``basis`` -- list of elements of `K^n`, where `K` is the fraction field
      of `R`; these elements must be linearly independent and will be used as
      the default basis of the constructed submodule

    - ``check`` -- boolean (default: ``True``); if ``False``, correctness of
      the input will not be checked and type conversion may be omitted, use
      with care

    - ``echelonize`` -- (default: ``False``) if ``True``, ``basis`` will be
      echelonized and the result will be used as the default basis of the
      constructed submodule;

    - ``echelonized_basis`` -- (default: ``None``) if not ``None``, must be
      the echelonized basis spanning the same submodule as ``basis``

    - ``already_echelonized`` -- boolean (default: ``False``); if ``True``,
      ``basis`` must be already given in the echelonized form

    OUTPUT: `R`-submodule of `K^n` with the user-specified ``basis``

    EXAMPLES::

        sage: M = ZZ^3
        sage: W = M.span_of_basis([[1,2,3], [4,5,6]]); W
        Free module of degree 3 and rank 2 over Integer Ring
        User basis matrix:
        [1 2 3]
        [4 5 6]

    Now we create a submodule of the ambient vector space, rather than
    ``M`` itself::

        sage: W = M.span_of_basis([[1,2,3/2], [4,5,6]]); W
        Free module of degree 3 and rank 2 over Integer Ring
        User basis matrix:
        [  1   2 3/2]
        [  4   5   6]
    """
    def __init__(self, ambient, basis, check=True,
                 echelonize=False, echelonized_basis=None,
                 already_echelonized=False,
                 category=None):
        r"""
        See :class:`FreeModule_submodule_with_basis_pid` for documentation.

        TESTS::

            sage: M = ZZ^3
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]])
            sage: TestSuite(W).run()

        We test that the issue at :issue:`9502` is solved::

            sage: parent(W.basis()[0])
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
            sage: parent(W.echelonized_basis()[0])
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]

        Now we test that the issue introduced at :issue:`9502` and reported at
        :issue:`10250` is solved as well::

            sage: V = (QQ^2).span_of_basis([[1,1]])
            sage: w = sqrt(2) * V([1,1])                                                # needs sage.symbolic
            sage: 3 * w                                                                 # needs sage.symbolic
            (3*sqrt(2), 3*sqrt(2))

        TESTS:

        Test that the category is determined as intended::

            sage: from sage.modules.free_module import FreeModule_ambient_pid, FreeModule_submodule_with_basis_pid
            sage: V = FreeModule_ambient_pid(QQ, 3, category=Algebras(QQ))
            sage: V.category()
            Category of finite dimensional algebras with basis over Rational Field
            sage: W = FreeModule_submodule_with_basis_pid(V, [[1,2,3]])
            sage: W.category()
            Join of
             Category of finite dimensional vector spaces with basis over (number fields and quotient fields and metric spaces) and
             Category of subobjects of sets
            sage: W = FreeModule_submodule_with_basis_pid(V, [[1,2,3]], category=Algebras(QQ))
            sage: W.category()
            Join of
             Category of finite dimensional algebras with basis over Rational Field and
             Category of subobjects of sets
        """
        if not isinstance(ambient, FreeModule_ambient_pid):
            raise TypeError("ambient (=%s) must be ambient." % ambient)
        self.__ambient_module = ambient
        R = ambient.base_ring()
        R_coord = R

        # Convert all basis elements to the ambient module
        try:
            basis = [ambient(x) for x in basis]
        except TypeError:
            # That failed, try the ambient vector space instead
            V = ambient.ambient_vector_space()
            R_coord = V.base_ring()
            try:
                basis = [V(x) for x in basis]
            except TypeError:
                raise TypeError("each element of basis must be in "
                                "the ambient vector space")

        if echelonize and not already_echelonized:
            basis = self._echelonized_basis(ambient, basis)

        # Adapted from Module_free_ambient.__init__
        from sage.categories.modules_with_basis import ModulesWithBasis
        modules_category = ModulesWithBasis(R.category()).FiniteDimensional()
        try:
            if R.is_finite() or len(basis) == 0:
                modules_category = modules_category.Enumerated().Finite()
        except (ValueError, TypeError, AttributeError, NotImplementedError):
            pass
        modules_category = modules_category.Subobjects()
        category = modules_category.or_subcategory(category, join=True)

        FreeModule_generic_pid.__init__(self, base_ring=R, coordinate_ring=R_coord,
                                        rank=len(basis), degree=ambient.degree(),
                                        sparse=ambient.is_sparse(), category=category)
        C = self.element_class
        w = [C(self, x.list(), coerce=False, copy=False) for x in basis]
        self.__basis = basis_seq(self, w)

        if echelonize or already_echelonized:
            self.__echelonized_basis = self.__basis
        else:
            if echelonized_basis is None:
                echelonized_basis = self._echelonized_basis(ambient, basis)
            w = [C(self, x.list(), coerce=False, copy=True)
                 for x in echelonized_basis]
            self.__echelonized_basis = basis_seq(self, w)
        if check and len(basis) != len(self.__echelonized_basis):
            raise ValueError("The given basis vectors must be linearly "
                             "independent.")

    def __hash__(self):
        """
        The hash is given by the basis.

        EXAMPLES::

            sage: M = ZZ^3
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]])
            sage: hash(W) == hash(W.basis())
            True
        """
        return hash(self.__basis)

    def _echelon_matrix_richcmp(self, other, op):
        r"""
        Compare the free module ``self`` with ``other``.

        Modules are ordered by their ambient spaces, then by dimension,
        then in order by their echelon matrices.

        .. NOTE::

           Use :meth:`is_submodule` to determine if one
           module is a submodule of another.

        EXAMPLES:

        First we compare two equal vector spaces.

        ::

            sage: from sage.structure.richcmp import op_LT,op_LE,op_EQ,op_NE,op_GT,op_GE
            sage: V = span([[1,2,3], [5,6,7], [8,9,10]], QQ)
            sage: W = span([[5,6,7], [8,9,10]], QQ)
            sage: V._echelon_matrix_richcmp(W,op_EQ)
            True

        Next we compare a one dimensional space to the two dimensional
        space defined above.

        ::

            sage: from sage.structure.richcmp import op_LT,op_LE,op_EQ,op_NE,op_GT,op_GE
            sage: M = span([[5,6,7]], QQ)
            sage: V._echelon_matrix_richcmp(M,op_EQ)
            False
            sage: M._echelon_matrix_richcmp(V, op_LT)
            True
            sage: V._echelon_matrix_richcmp(M, op_LT)
            False

        We compare a `\ZZ`-module to the one-dimensional
        space above.::

            sage: from sage.structure.richcmp import op_LT,op_LE,op_EQ,op_NE,op_GT,op_GE
            sage: V = span([[5,6,7]], ZZ).scale(1/11);  V
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [5/11 6/11 7/11]
            sage: V._echelon_matrix_richcmp(M, op_LT)
            True
            sage: M._echelon_matrix_richcmp(V, op_LT)
            False
        """
        if self is other:
            return rich_to_bool(op, 0)
        if not isinstance(other, FreeModule_generic):
            return NotImplemented
        lx = self.ambient_vector_space()
        rx = other.ambient_vector_space()
        if lx != rx:
            return lx._echelon_matrix_richcmp(rx, op)

        lx = self.dimension()
        rx = other.dimension()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.base_ring()
        rx = other.base_ring()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        # We use self.echelonized_basis_matrix() == other.echelonized_basis_matrix()
        # with the matrix to avoid a circular reference.
        return richcmp(self.echelonized_basis_matrix(),
                       other.echelonized_basis_matrix(), op)

    def construction(self):
        """
        Return the functorial construction of ``self``, namely, the subspace
        of the ambient module spanned by the given basis.

        EXAMPLES::

            sage: M = ZZ^3
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]]); W
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
            sage: c, V = W.construction()
            sage: c(V) == W
            True
        """
        from sage.categories.pushout import SubspaceFunctor
        return SubspaceFunctor(self.basis()), self.ambient_module()

    def echelonized_basis_matrix(self):
        """
        Return basis matrix for ``self`` in row echelon form.

        EXAMPLES::

            sage: V = FreeModule(ZZ, 3).span_of_basis([[1,2,3],[4,5,6]])
            sage: V.basis_matrix()
            [1 2 3]
            [4 5 6]
            sage: V.echelonized_basis_matrix()
            [1 2 3]
            [0 3 6]
        """
        try:
            return self.__echelonized_basis_matrix
        except AttributeError:
            pass
        self._echelonized_basis(self.ambient_module(), self.__basis)
        return self.__echelonized_basis_matrix

    def _echelonized_basis(self, ambient, basis):
        """
        Given the ambient space and a basis, construct and cache the
        echelonized basis matrix and returns its rows.

        N.B. This function is for internal use only!

        EXAMPLES::

            sage: M = ZZ^3
            sage: N = M.submodule_with_basis([[1,1,0],[0,2,1]])
            sage: N._echelonized_basis(M,N.basis())
            [(1, 1, 0), (0, 2, 1)]
            sage: V = QQ^3
            sage: W = V.submodule_with_basis([[1,1,0],[0,2,1]])
            sage: W._echelonized_basis(V,W.basis())
            [(1, 0, -1/2), (0, 1, 1/2)]
            sage: V = SR^3                                                              # needs sage.symbolic
            sage: W = V.submodule_with_basis([[1,0,1]])
            sage: W._echelonized_basis(V, W.basis())
            [(1, 0, 1)]
        """
        # Return the first rank rows (i.e., the nonzero rows).
        d = self._denominator(basis)
        MAT = sage.matrix.matrix_space.MatrixSpace(
            ambient.base_ring(), len(basis), ambient.degree(), sparse=ambient.is_sparse())
        if d != 1:
            basis = [x * d for x in basis]
        A = MAT(basis)
        E = A.echelon_form()
        if d != 1:
            E = E.matrix_over_field() * (~d)   # divide out denominator
        r = E.rank()
        if r < E.nrows():
            E = E.matrix_from_rows(range(r))
        self.__echelonized_basis_matrix = E
        return E.rows()

    def _denominator(self, B):
        """
        The LCM of the denominators of the given list B.

        N.B.: This function is for internal use only!

        EXAMPLES::

            sage: V = QQ^3
            sage: L = V.span([[1,1/2,1/3], [-1/5,2/3,3]], ZZ)
            sage: L
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1/5 19/6 37/3]
            [   0 23/6 46/3]
            sage: L._denominator(L.echelonized_basis_matrix().list())
            30
        """
        if not B:
            return 1
        d = B[0].denominator()
        from sage.arith.functions import lcm
        for x in B[1:]:
            d = lcm(d, x.denominator())
        return d

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: L = ZZ^8
            sage: E = L.submodule_with_basis([ L.gen(i) - L.gen(0) for i in range(1,8) ])
            sage: E # indirect doctest
            Free module of degree 8 and rank 7 over Integer Ring
            User basis matrix:
            [-1  1  0  0  0  0  0  0]
            [-1  0  1  0  0  0  0  0]
            [-1  0  0  1  0  0  0  0]
            [-1  0  0  0  1  0  0  0]
            [-1  0  0  0  0  1  0  0]
            [-1  0  0  0  0  0  1  0]
            [-1  0  0  0  0  0  0  1]

        ::

            sage: M = FreeModule(ZZ,8,sparse=True)
            sage: N = M.submodule_with_basis([ M.gen(i) - M.gen(0) for i in range(1,8) ])
            sage: N # indirect doctest
            Sparse free module of degree 8 and rank 7 over Integer Ring
            User basis matrix:
            [-1  1  0  0  0  0  0  0]
            [-1  0  1  0  0  0  0  0]
            [-1  0  0  1  0  0  0  0]
            [-1  0  0  0  1  0  0  0]
            [-1  0  0  0  0  1  0  0]
            [-1  0  0  0  0  0  1  0]
            [-1  0  0  0  0  0  0  1]
        """
        if self.is_sparse():
            s = "Sparse free module of degree %s and rank %s over %s\n" % (
                self.degree(), self.rank(), self.base_ring()) + \
                "User basis matrix:\n%r" % self.basis_matrix()
        else:
            s = "Free module of degree %s and rank %s over %s\n" % (
                self.degree(), self.rank(), self.base_ring()) + \
                "User basis matrix:\n%r" % self.basis_matrix()
        return s

    def _latex_(self):
        r"""
        Return latex representation of this free module.

        EXAMPLES::

            sage: A = ZZ^3
            sage: M = A.span_of_basis([[1,2,3],[4,5,6]])
            sage: M._latex_()
            '\\mathrm{RowSpan}_{\\Bold{Z}}\\left(\\begin{array}{rrr}\n1 & 2 & 3 \\\\\n4 & 5 & 6\n\\end{array}\\right)'
        """
        return "\\mathrm{RowSpan}_{%s}%s" % (latex.latex(self.base_ring()), latex.latex(self.basis_matrix()))

    def ambient_module(self):
        """
        Return the ambient module related to the `R`-module self,
        which was used when creating this module, and is of the form
        `R^n`. Note that ``self`` need not be contained in the ambient
        module, though ``self`` will be contained in the ambient vector space.

        EXAMPLES::

            sage: A = ZZ^3
            sage: M = A.span_of_basis([[1,2,'3/7'],[4,5,6]])
            sage: M
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [  1   2 3/7]
            [  4   5   6]
            sage: M.ambient_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: M.is_submodule(M.ambient_module())
            False
        """
        return self.__ambient_module

    # Sets.Subquotients.ParentMethods
    def ambient(self):
        """
        Return the ambient module or space for ``self``.

        EXAMPLES::

            sage: M = ZZ^3
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]]); W
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
            sage: W.ambient()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring

        Now we create a submodule of the ambient vector space, rather than
        ``M`` itself::

            sage: W = M.span_of_basis([[1,2,3/2],[4,5,6]]); W
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [  1   2 3/2]
            [  4   5   6]
            sage: W.ambient()
            Vector space of dimension 3 over Rational Field

        A submodule of a submodule::

            sage: M = ZZ^3
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]]); W
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
            sage: U = W.span_of_basis([[5,7,9]]); U
            Free module of degree 3 and rank 1 over Integer Ring
            User basis matrix:
            [5 7 9]
            sage: U.ambient()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        if self.base_ring() == self.coordinate_ring():
            return self.ambient_module()
        else:
            return self.ambient_vector_space()

    # Sets.Subquotients.ParentMethods
    @lazy_attribute
    def lift(self):
        r"""
        The lift (embedding) map from ``self`` to the ambient module or space.

        EXAMPLES::

            sage: M = ZZ^3
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]]); W
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
            sage: W.lift
            Generic morphism:
            From: Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
            To:   Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: w = W([5,7,9])
            sage: m = W.lift(w); m
            (5, 7, 9)
            sage: m.parent()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        ambient = self.ambient()
        return self.module_morphism(function=ambient, codomain=ambient)

    # Sets.Subquotients.ParentMethods
    @lazy_attribute
    def retract(self):
        r"""
        The retract map from the ambient space.

        This is a partial map, which gives an error for elements not in the subspace.

        Calling this map on elements of the ambient space is the same as calling the
        element constructor of ``self``.

        EXAMPLES::

            sage: M = ZZ^3
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]]); W
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
            sage: W.retract
            Generic morphism:
            From: Ambient free module of rank 3 over the principal ideal domain Integer Ring
            To:   Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
            sage: m = M([5, 7, 9])
            sage: w = W.retract(m); w
            (5, 7, 9)
            sage: w.parent()
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
        """
        return self.ambient().module_morphism(function=self, codomain=self)

    def relations(self):
        r"""
        Return the submodule defining the relations of ``self`` as a
        subquotient (considering the ambient module as a quotient module).

        EXAMPLES::

            sage: V = GF(2)^2
            sage: W = V.subspace([[1, 0]])
            sage: W.relations() == V.zero_submodule()
            True

            sage: Q = V / W
            sage: Q.relations() == W
            True
            sage: Q.zero_submodule().relations() == W
            True
        """
        return self.__ambient_module.relations()

    def echelon_coordinates(self, v, check=True):
        r"""
        Write `v` in terms of the echelonized basis for ``self``.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also
          verify that `v` is really in ``self``

        OUTPUT: list

        Returns a list `c` such that if `B` is the basis
        for self, then

        .. MATH::

            \sum c_i B_i = v.

        If `v` is not in self, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: A = ZZ^3
            sage: M = A.span_of_basis([[1,2,'3/7'],[4,5,6]])
            sage: M.coordinates([8,10,12])
            [0, 2]
            sage: M.echelon_coordinates([8,10,12])
            [8, -2]
            sage: B = M.echelonized_basis(); B
            [(1, 2, 3/7), (0, 3, -30/7)]
            sage: 8*B[0] - 2*B[1]
            (8, 10, 12)

        We do an example with a sparse vector space::

            sage: V = VectorSpace(QQ,5, sparse=True)
            sage: W = V.subspace_with_basis([[0,1,2,0,0], [0,-1,0,0,-1/2]])
            sage: W.echelonized_basis()
            [(0, 1, 0, 0, 1/2), (0, 0, 1, 0, -1/4)]
            sage: W.echelon_coordinates([0,0,2,0,-1/2])
            [0, 2]
        """
        if not isinstance(v, free_module_element.FreeModuleElement):
            v = self.ambient_vector_space()(v)
        elif v.degree() != self.degree():
            raise ArithmeticError("vector is not in free module")
        # Find coordinates of v with respect to rref basis.
        E = self.echelonized_basis_matrix()
        P = E.pivots()
        w = v.list_from_positions(P)
        # Next use the transformation matrix from the rref basis
        # to the echelon basis.
        T = self._rref_to_echelon_matrix()
        x = T.linear_combination_of_rows(w).list(copy=False)
        if not check:
            return x
        if v.parent() is self:
            return x
        lc = E.linear_combination_of_rows(x)
        if list(lc) != list(v):
            raise ArithmeticError("vector is not in free module")
        return x

    def user_to_echelon_matrix(self):
        """
        Return matrix that transforms a vector written with respect to the
        user basis of ``self`` to one written with respect to the echelon
        basis. The matrix acts from the right, as is usual in Sage.

        EXAMPLES::

            sage: A = ZZ^3
            sage: M = A.span_of_basis([[1,2,3],[4,5,6]])
            sage: M.echelonized_basis()
            [(1, 2, 3), (0, 3, 6)]
            sage: M.user_to_echelon_matrix()
            [ 1  0]
            [ 4 -1]

        The vector `v=(5,7,9)` in `M` is `(1,1)`
        with respect to the user basis. Multiplying the above matrix on the
        right by this vector yields `(5,-1)`, which has components
        the coordinates of `v` with respect to the echelon basis.

        ::

            sage: v0,v1 = M.basis(); v = v0+v1
            sage: e0,e1 = M.echelonized_basis()
            sage: v
            (5, 7, 9)
            sage: 5*e0 + (-1)*e1
            (5, 7, 9)
        """
        try:
            return self.__user_to_echelon_matrix
        except AttributeError:
            if self.base_ring().is_field():
                self.__user_to_echelon_matrix = self._user_to_rref_matrix()
            else:
                rows = sum([self.echelon_coordinates(b, check=False)
                            for b in self.basis()], [])
                M = sage.matrix.matrix_space.MatrixSpace(self.base_ring().fraction_field(),
                                                         self.dimension(),
                                                         sparse=self.is_sparse())
                self.__user_to_echelon_matrix = M(rows)
        return self.__user_to_echelon_matrix

    def echelon_to_user_matrix(self):
        """
        Return matrix that transforms the echelon basis to the user basis
        of ``self``. This is a matrix `A` such that if `v` is a
        vector written with respect to the echelon basis for ``self`` then
        `vA` is that vector written with respect to the user basis
        of ``self``.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.span_of_basis([[1,2,3],[4,5,6]])
            sage: W.echelonized_basis()
            [(1, 0, -1), (0, 1, 2)]
            sage: A = W.echelon_to_user_matrix(); A
            [-5/3  2/3]
            [ 4/3 -1/3]

        The vector `(1,1,1)` has coordinates `v=(1,1)` with
        respect to the echelonized basis for ``self``. Multiplying `vA`
        we find the coordinates of this vector with respect to the user
        basis.

        ::

            sage: v = vector(QQ, [1,1]); v
            (1, 1)
            sage: v * A
            (-1/3, 1/3)
            sage: u0, u1 = W.basis()
            sage: (-u0 + u1)/3
            (1, 1, 1)
        """
        try:
            return self.__echelon_to_user_matrix
        except AttributeError:
            self.__echelon_to_user_matrix = ~self.user_to_echelon_matrix()
            return self.__echelon_to_user_matrix

    def _user_to_rref_matrix(self):
        """
        Return a transformation matrix from the user specified basis to row
        reduced echelon form, for this module over a PID.

        Note: For internal use only! See ``user_to_echelon_matrix``.

        EXAMPLES::

            sage: M = ZZ^3
            sage: N = M.submodule_with_basis([[1,1,0],[0,1,1]])
            sage: T = N.user_to_echelon_matrix(); T # indirect doctest
            [1 1]
            [0 1]
            sage: N.basis_matrix()
            [1 1 0]
            [0 1 1]
            sage: N.echelonized_basis_matrix()
            [ 1  0 -1]
            [ 0  1  1]
            sage: T * N.echelonized_basis_matrix() == N.basis_matrix()
            True
        """
        try:
            return self.__user_to_rref_matrix
        except AttributeError:
            A = self.basis_matrix()
            P = self.echelonized_basis_matrix().pivots()
            T = A.matrix_from_columns(P)
            self.__user_to_rref_matrix = T
        return self.__user_to_rref_matrix

    def _rref_to_user_matrix(self):
        """
        Return a transformation matrix from row reduced echelon form to
        the user specified basis, for this module over a PID.

        Note: For internal use only! See ``user_to_echelon_matrix``.

        EXAMPLES::

            sage: M = ZZ^3
            sage: N = M.submodule_with_basis([[1,1,0],[0,1,1]])
            sage: U = N.echelon_to_user_matrix(); U # indirect doctest
            [ 1 -1]
            [ 0  1]
            sage: N.echelonized_basis_matrix()
            [ 1  0 -1]
            [ 0  1  1]
            sage: N.basis_matrix()
            [1 1 0]
            [0 1 1]
            sage: U * N.basis_matrix() == N.echelonized_basis_matrix()
            True
        """
        try:
            return self.__rref_to_user_matrix
        except AttributeError:
            self.__rref_to_user_matrix = ~self._user_to_rref_matrix()
            return self.__rref_to_user_matrix

    def _echelon_to_rref_matrix(self):
        """
        Return a transformation matrix from the some matrix to the row
        reduced echelon form for this module over a PID.

        Note: For internal use only! and not used!

        EXAMPLES::

            sage: M = ZZ^3
            sage: N = M.submodule_with_basis([[1,1,0],[1,1,2]])
            sage: N
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 1 0]
            [1 1 2]
            sage: T = N._echelon_to_rref_matrix(); T
            [1 0]
            [0 2]
            sage: type(T)
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: U = N._rref_to_echelon_matrix(); U
            [  1   0]
            [  0 1/2]
            sage: type(U)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
        """
        try:
            return self.__echelon_to_rref_matrix
        except AttributeError:
            A = self.echelonized_basis_matrix()
            T = A.matrix_from_columns(A.pivots())
            self.__echelon_to_rref_matrix = T
        return self.__echelon_to_rref_matrix

    def _rref_to_echelon_matrix(self):
        """
        Return a transformation matrix from row reduced echelon form to
        some matrix for this module over a PID.

        Note: For internal use only!

        EXAMPLES::

            sage: M = ZZ^3
            sage: N = M.submodule_with_basis([[1,1,0],[1,1,2]])
            sage: N
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 1 0]
            [1 1 2]
            sage: T = N._echelon_to_rref_matrix(); T
            [1 0]
            [0 2]
            sage: type(T)
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: U = N._rref_to_echelon_matrix(); U
            [  1   0]
            [  0 1/2]
            sage: type(U)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
        """
        try:
            return self.__rref_to_echelon_matrix
        except AttributeError:
            self.__rref_to_echelon_matrix = ~self._echelon_to_rref_matrix()
            return self.__rref_to_echelon_matrix

    def vector_space(self, base_field=None):
        """
        Return the vector space associated to this free module via tensor
        product with the fraction field of the base ring.

        EXAMPLES::

            sage: A = ZZ^3; A
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: A.vector_space()
            Vector space of dimension 3 over Rational Field
            sage: M = A.span_of_basis([['1/3',2,'3/7'],[4,5,6]]); M
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1/3   2 3/7]
            [  4   5   6]
            sage: M.vector_space()
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1/3   2 3/7]
            [  4   5   6]
        """
        if base_field is None:
            V = self.ambient_vector_space()
            return V.submodule_with_basis(self.basis())
        return self.change_ring(base_field)

    def ambient_vector_space(self):
        """
        Return the ambient vector space in which this free module is
        embedded.

        EXAMPLES::

            sage: M = ZZ^3;  M.ambient_vector_space()
            Vector space of dimension 3 over Rational Field

        ::

            sage: N = M.span_of_basis([[1,2,'1/5']])
            sage: N
            Free module of degree 3 and rank 1 over Integer Ring
            User basis matrix:
            [  1   2 1/5]
            sage: M.ambient_vector_space()
            Vector space of dimension 3 over Rational Field
            sage: M.ambient_vector_space() is N.ambient_vector_space()
            True

        If an inner product on the module is specified, then this
        is preserved on the ambient vector space.

        ::

            sage: M = FreeModule(ZZ,4,inner_product_matrix=1)
            sage: V = M.ambient_vector_space()
            sage: V
            Ambient quadratic space of dimension 4 over Rational Field
            Inner product matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: N = M.submodule([[1,-1,0,0],[0,1,-1,0],[0,0,1,-1]])
            sage: N.gram_matrix()
            [2 1 1]
            [1 2 1]
            [1 1 2]
            sage: V == N.ambient_vector_space()
            True
        """
        return self.ambient_module().ambient_vector_space()

    def basis(self):
        """
        Return the user basis for this free module.

        EXAMPLES::

            sage: V = ZZ^3
            sage: V.basis()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
            sage: M = V.span_of_basis([['1/8',2,1]])
            sage: M.basis()
            [(1/8, 2, 1)]
        """
        return self.__basis

    def change_ring(self, R):
        """
        Return the free module over `R` obtained by coercing each
        element of the basis of ``self`` into a vector over the
        fraction field of `R`, then taking the resulting `R`-module.

        INPUT:

        - ``R`` -- a principal ideal domain

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.subspace([[2, 1/2, 1]])
            sage: W.change_ring(GF(7))
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 2 4]

        ::

            sage: M = (ZZ^2) * (1/2)
            sage: N = M.change_ring(QQ)
            sage: N
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            sage: N = M.change_ring(QQ['x'])
            sage: N
            Free module of degree 2 and rank 2
             over Univariate Polynomial Ring in x over Rational Field
            Echelon basis matrix:
            [1/2   0]
            [  0 1/2]
            sage: N.coordinate_ring()
            Univariate Polynomial Ring in x over Rational Field

        The ring must be a principal ideal domain::

            sage: M.change_ring(ZZ['x'])
            Traceback (most recent call last):
            ...
            TypeError: the new ring Univariate Polynomial Ring in x over Integer Ring
            should be a principal ideal domain
        """
        if self.base_ring() is R:
            return self
        if R not in PrincipalIdealDomains():
            raise TypeError("the new ring %r should be a principal ideal domain" % R)

        K = R.fraction_field()
        V = VectorSpace(K, self.degree())
        B = [V(b) for b in self.basis()]
        M = self.ambient_module().change_ring(R)
        if self.has_user_basis():
            return M.span_of_basis(B)
        else:
            return M.span(B)

    def coordinate_vector(self, v, check=True):
        """
        Write `v` in terms of the user basis for ``self``.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also verify that
          `v` is really in ``self``

        OUTPUT: list

        The output is a vector `c` such that if `B` is the basis for ``self``,
        then

        .. MATH::

            \\sum c_i B_i = v.

        If `v` is not in self, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: V = ZZ^3
            sage: M = V.span_of_basis([['1/8',2,1]])
            sage: M.coordinate_vector([1,16,8])
            (8)
        """
        # First find the coordinates of v wrt echelon basis.
        w = self.echelon_coordinate_vector(v, check=check)
        # Next use transformation matrix from echelon basis to
        # user basis.
        T = self.echelon_to_user_matrix()
        return T.linear_combination_of_rows(w)

    def echelonized_basis(self):
        """
        Return the basis for ``self`` in echelon form.

        EXAMPLES::

            sage: V = ZZ^3
            sage: M = V.span_of_basis([['1/2',3,1], [0,'1/6',0]])
            sage: M.basis()
            [(1/2, 3, 1), (0, 1/6, 0)]
            sage: B = M.echelonized_basis(); B
            [(1/2, 0, 1), (0, 1/6, 0)]
            sage: V.span(B) == M
            True
        """
        return self.__echelonized_basis

    def echelon_coordinate_vector(self, v, check=True):
        """
        Write `v` in terms of the echelonized basis for ``self``.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also
          verify that `v` is really in ``self``

        Returns a list `c` such that if `B` is the echelonized basis
        for ``self``, then

        .. MATH::

            \\sum c_i B_i = v.

        If `v` is not in self, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: V = ZZ^3
            sage: M = V.span_of_basis([['1/2',3,1], [0,'1/6',0]])
            sage: B = M.echelonized_basis(); B
            [(1/2, 0, 1), (0, 1/6, 0)]
            sage: M.echelon_coordinate_vector(['1/2', 3, 1])
            (1, 18)
        """
        return FreeModule(self.base_ring().fraction_field(), self.rank())(self.echelon_coordinates(v, check=check))

    def has_user_basis(self):
        """
        Return ``True`` if the basis of this free module is
        specified by the user, as opposed to being the default echelon
        form.

        EXAMPLES::

            sage: V = ZZ^3; V.has_user_basis()
            False
            sage: M = V.span_of_basis([[1,3,1]]); M.has_user_basis()
            True
            sage: M = V.span([[1,3,1]]); M.has_user_basis()
            False
        """
        return True

    def linear_combination_of_basis(self, v):
        """
        Return the linear combination of the basis for ``self`` obtained from
        the coordinates of v.

        INPUT:

        - ``v`` -- list

        EXAMPLES::

            sage: V = span([[1,2,3], [4,5,6]], ZZ); V
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            [0 3 6]
            sage: V.linear_combination_of_basis([1,1])
            (1, 5, 9)

        This should raise an error if the resulting element is not in self::

            sage: W = (QQ**2).span([[2, 0], [0, 8]], ZZ)
            sage: W.linear_combination_of_basis([1, -1/2])
            Traceback (most recent call last):
            ...
            TypeError: element [2, -4] is not in free module
        """
        R = self.base_ring()
        check = (not R.is_field()) and any(a not in R for a in list(v))
        return self(self.basis_matrix().linear_combination_of_rows(v),
                    check=check, copy=False, coerce=False)


class FreeModule_submodule_pid(FreeModule_submodule_with_basis_pid):
    """
    An `R`-submodule of `K^n` where `K` is the
    fraction field of a principal ideal domain `R`.

    EXAMPLES::

        sage: M = ZZ^3
        sage: W = M.span_of_basis([[1,2,3],[4,5,19]]); W
        Free module of degree 3 and rank 2 over Integer Ring
        User basis matrix:
        [ 1  2  3]
        [ 4  5 19]

    Generic tests, including saving and loading submodules and elements::

        sage: TestSuite(W).run()
        sage: v = W.0 + W.1
        sage: TestSuite(v).run()
    """
    def __init__(self, ambient, gens, check=True, already_echelonized=False,
                 category=None):
        """
        Create an embedded free module over a PID.

        EXAMPLES::

            sage: V = ZZ^3
            sage: W = V.span([[1,2,3],[4,5,6]])
            sage: W
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            [0 3 6]
        """
        FreeModule_submodule_with_basis_pid.__init__(self, ambient, basis=gens,
                                                     echelonize=True,
                                                     already_echelonized=already_echelonized,
                                                     category=category)

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: M = ZZ^8
            sage: L = M.submodule([ M.gen(i) - M.gen(0) for i in range(1,8) ])
            sage: L # indirect doctest
            Free module of degree 8 and rank 7 over Integer Ring
            Echelon basis matrix:
            [ 1  0  0  0  0  0  0 -1]
            [ 0  1  0  0  0  0  0 -1]
            [ 0  0  1  0  0  0  0 -1]
            [ 0  0  0  1  0  0  0 -1]
            [ 0  0  0  0  1  0  0 -1]
            [ 0  0  0  0  0  1  0 -1]
            [ 0  0  0  0  0  0  1 -1]
        """
        if self.is_sparse():
            s = "Sparse free module of degree %s and rank %s over %s\n" % (
                self.degree(), self.rank(), self.base_ring()) + \
                "Echelon basis matrix:\n%s" % self.basis_matrix()
        else:
            s = "Free module of degree %s and rank %s over %s\n" % (
                self.degree(), self.rank(), self.base_ring()) + \
                "Echelon basis matrix:\n%s" % self.basis_matrix()
        return s

    def coordinate_vector(self, v, check=True):
        """
        Write `v` in terms of the user basis for ``self``.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also verify that
           `v` is really in ``self``

        OUTPUT: list

        The output is a list `c` such that if `B` is the basis for ``self``,
        then

        .. MATH::

            \\sum c_i B_i = v.

        If `v` is not in self, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: V = ZZ^3
            sage: W = V.span_of_basis([[1,2,3],[4,5,6]])
            sage: W.coordinate_vector([1,5,9])
            (5, -1)
        """
        return self.echelon_coordinate_vector(v, check=check)

    def has_user_basis(self):
        r"""
        Return ``True`` if the basis of this free module is
        specified by the user, as opposed to being the default echelon
        form.

        EXAMPLES::

            sage: A = ZZ^3; A
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: A.has_user_basis()
            False
            sage: W = A.span_of_basis([[2,'1/2',1]])
            sage: W.has_user_basis()
            True
            sage: W = A.span([[2,'1/2',1]])
            sage: W.has_user_basis()
            False
        """
        return False


FreeModule_generic_pid._submodule_class = FreeModule_submodule_pid


class FreeModule_submodule_with_basis_field(FreeModule_generic_field, FreeModule_submodule_with_basis_pid):
    """
    An embedded vector subspace with a distinguished user basis.

    EXAMPLES::

        sage: M = QQ^3; W = M.submodule_with_basis([[1,2,3], [4,5,19]]); W
        Vector space of degree 3 and dimension 2 over Rational Field
        User basis matrix:
        [ 1  2  3]
        [ 4  5 19]

    Since this is an embedded vector subspace with a distinguished user
    basis possibly different than the echelonized basis, the
    echelon_coordinates() and user coordinates() do not agree::

        sage: V = QQ^3

    ::

        sage: W = V.submodule_with_basis([[1,2,3], [4,5,6]])
        sage: W
        Vector space of degree 3 and dimension 2 over Rational Field
        User basis matrix:
        [1 2 3]
        [4 5 6]

    ::

        sage: v = V([1,5,9])
        sage: W.echelon_coordinates(v)
        [1, 5]
        sage: vector(QQ, W.echelon_coordinates(v)) * W.echelonized_basis_matrix()
        (1, 5, 9)

    ::

        sage: v = V([1,5,9])
        sage: W.coordinates(v)
        [5, -1]
        sage: vector(QQ, W.coordinates(v)) * W.basis_matrix()
        (1, 5, 9)

    Generic tests, including saving and loading submodules and elements::

        sage: TestSuite(W).run()

        sage: K.<x> = FractionField(PolynomialRing(QQ,'x'))
        sage: M = K^3; W = M.span_of_basis([[1,1,x]])
        sage: TestSuite(W).run()
    """
    def __init__(self, ambient, basis, check=True,
                 echelonize=False, echelonized_basis=None,
                 already_echelonized=False,
                 category=None):
        """
        Create a vector space with given basis.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.span_of_basis([[1,2,3],[4,5,6]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 2 3]
            [4 5 6]
        """
        FreeModule_submodule_with_basis_pid.__init__(
            self, ambient, basis=basis, check=check, echelonize=echelonize,
            echelonized_basis=echelonized_basis, already_echelonized=already_echelonized,
            category=category)

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: V = VectorSpace(QQ,5)
            sage: U = V.submodule([ V.gen(i) - V.gen(0) for i in range(1,5) ])
            sage: U # indirect doctest
            Vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
            sage: print(U._repr_())
            Vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]

        The system representation can be overwritten, but leaves _repr_
        unmodified.

        ::

            sage: U.rename('U')
            sage: U
            U
            sage: print(U._repr_())
            Vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]

        Sparse vector spaces print this fact.

        ::

            sage: VV = VectorSpace(QQ,5,sparse=True)
            sage: UU = VV.submodule([ VV.gen(i) - VV.gen(0) for i in range(1,5) ])
            sage: UU # indirect doctest
            Sparse vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]

        (Now clean up again.)

        ::

            sage: U.reset_name()
            sage: U
            Vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
        """
        if self.is_sparse():
            return "Sparse vector space of degree %s and dimension %s over %s\n" % (
                self.degree(), self.dimension(), self.base_field()) + \
                "User basis matrix:\n%r" % self.basis_matrix()

        return "Vector space of degree %s and dimension %s over %s\n" % (
            self.degree(), self.dimension(), self.base_field()) + \
            "User basis matrix:\n%r" % self.basis_matrix()

    def _denominator(self, B):
        """
        Given a list (of field elements) returns 1 as the common
        denominator.

        N.B.: This function is for internal use only!

        EXAMPLES::

            sage: U = QQ^3
            sage: U
            Vector space of dimension 3 over Rational Field
            sage: U.denominator()
            1
            sage: V = U.span([[1,1/2,1/3], [-1/5,2/3,3]])
            sage: V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0 -5/3]
            [   0    1    4]
            sage: W = U.submodule_with_basis([[1,1/2,1/3], [-1/5,2/3,3]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [   1  1/2  1/3]
            [-1/5  2/3    3]
            sage: W._denominator(W.echelonized_basis_matrix().list())
            1
        """
        return 1

    def _echelonized_basis(self, ambient, basis):
        """
        Given the ambient space and a basis, construct and cache the
        echelonized basis matrix and returns its rows.

        N.B. This function is for internal use only!

        EXAMPLES::

            sage: M = ZZ^3
            sage: N = M.submodule_with_basis([[1,1,0],[0,2,1]])
            sage: N._echelonized_basis(M,N.basis())
            [(1, 1, 0), (0, 2, 1)]
            sage: V = QQ^3
            sage: W = V.submodule_with_basis([[1,1,0],[0,2,1]])
            sage: W._echelonized_basis(V,W.basis())
            [(1, 0, -1/2), (0, 1, 1/2)]
        """
        MAT = sage.matrix.matrix_space.MatrixSpace(
            base_ring=ambient.base_ring(),
            nrows=len(basis), ncols=ambient.degree(),
            sparse=ambient.is_sparse())
        A = MAT(basis)
        E = A.echelon_form()
        # Return the first rank rows (i.e., the nonzero rows).
        return E.rows()[:E.rank()]

    def is_ambient(self):
        """
        Return False since this is not an ambient module.

        EXAMPLES::

            sage: V = QQ^3
            sage: V.is_ambient()
            True
            sage: W = V.span_of_basis([[1,2,3],[4,5,6]])
            sage: W.is_ambient()
            False
        """
        return False


class FreeModule_submodule_field(FreeModule_submodule_with_basis_field):
    """
    An embedded vector subspace with echelonized basis.

    EXAMPLES:

    Since this is an embedded vector subspace with echelonized basis,
    the echelon_coordinates() and user coordinates() agree::

        sage: V = QQ^3
        sage: W = V.span([[1,2,3],[4,5,6]])
        sage: W
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -1]
        [ 0  1  2]

    ::

        sage: v = V([1,5,9])
        sage: W.echelon_coordinates(v)
        [1, 5]
        sage: vector(QQ, W.echelon_coordinates(v)) * W.basis_matrix()
        (1, 5, 9)
        sage: v = V([1,5,9])
        sage: W.coordinates(v)
        [1, 5]
        sage: vector(QQ, W.coordinates(v)) * W.basis_matrix()
        (1, 5, 9)
    """
    def __init__(self, ambient, gens, check=True, already_echelonized=False, category=None):
        """
        Create an embedded vector subspace with echelonized basis.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.span([[1,2,3],[4,5,6]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
        """
        if isinstance(gens, FreeModule_generic):
            gens = gens.gens()
        FreeModule_submodule_with_basis_field.__init__(self, ambient,
                                                       basis=gens, check=check,
                                                       echelonize=not already_echelonized,
                                                       already_echelonized=already_echelonized,
                                                       category=category)

    def _repr_(self) -> str:
        """
        The default printing representation of ``self``.

        EXAMPLES::

            sage: V = VectorSpace(QQ,5)
            sage: U = V.submodule([ V.gen(i) - V.gen(0) for i in range(1,5) ])
            sage: U # indirect doctest
            Vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
            sage: print(U._repr_())
            Vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]

        The system representation can be overwritten, but leaves _repr_
        unmodified.

        ::

            sage: U.rename('U')
            sage: U
            U
            sage: print(U._repr_())
            Vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]

        Sparse vector spaces print this fact.

        ::

            sage: VV = VectorSpace(QQ,5,sparse=True)
            sage: UU = VV.submodule([ VV.gen(i) - VV.gen(0) for i in range(1,5) ])
            sage: UU # indirect doctest
            Sparse vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]

        (Now clean up again.)

        ::

            sage: U.reset_name()
            sage: U
            Vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
        """
        if self.is_sparse():
            return "Sparse vector space of degree %s and dimension %s over %s\n" % (
                self.degree(), self.dimension(), self.base_field()) + \
                "Basis matrix:\n%r" % self.basis_matrix()
        else:
            return "Vector space of degree %s and dimension %s over %s\n" % (
                self.degree(), self.dimension(), self.base_field()) + \
                "Basis matrix:\n%r" % self.basis_matrix()

    def echelon_coordinates(self, v, check=True):
        """
        Write `v` in terms of the echelonized basis of ``self``.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also
          verify that `v` is really in ``self``

        OUTPUT: list

        The output is a list `c` such that if `B` is the basis for ``self``,
        then

        .. MATH::

            \\sum c_i B_i = v.

        If `v` is not in ``self``, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.span([[1,2,3],[4,5,6]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]

        ::

            sage: v = V([1,5,9])
            sage: W.echelon_coordinates(v)
            [1, 5]
            sage: vector(QQ, W.echelon_coordinates(v)) * W.basis_matrix()
            (1, 5, 9)
        """
        if not isinstance(v, free_module_element.FreeModuleElement):
            v = self.ambient_vector_space()(v)
        if v.degree() != self.degree():
            raise ArithmeticError("v (=%s) is not in self" % v)
        E = self.echelonized_basis_matrix()
        P = E.pivots()
        if not P:
            if check and v != 0:
                raise ArithmeticError("vector is not in free module")
            return []
        w = v.list_from_positions(P)
        if not check:
            # It's really really easy.
            return w
        if v.parent() is self:   # obvious that v is really in here.
            return w
        # the "linear_combination_of_rows" call dominates the runtime
        # of this function, in the check==False case when the parent
        # of v is not self.
        lc = E.linear_combination_of_rows(w)
        if lc.list() != v.list():
            raise ArithmeticError("vector is not in free module")
        return w

    def coordinate_vector(self, v, check=True):
        """
        Write `v` in terms of the user basis for ``self``.

        INPUT:

        - ``v`` -- vector

        - ``check`` -- boolean (default: ``True``); if ``True``, also verify that
           `v` is really in ``self``

        OUTPUT: list

        The output is a list `c` such that if `B` is the basis for ``self``, then

        .. MATH::

            \\sum c_i B_i = v.

        If `v` is not in ``self``, raise an :exc:`ArithmeticError` exception.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.span([[1,2,3],[4,5,6]]); W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: v = V([1,5,9])
            sage: W.coordinate_vector(v)
            (1, 5)
            sage: W.coordinates(v)
            [1, 5]
            sage: vector(QQ, W.coordinates(v)) * W.basis_matrix()
            (1, 5, 9)

        ::

            sage: V = VectorSpace(QQ,5, sparse=True)
            sage: W = V.subspace([[0,1,2,0,0], [0,-1,0,0,-1/2]])
            sage: W.coordinate_vector([0,0,2,0,-1/2])
            (0, 2)
        """
        return self.echelon_coordinate_vector(v, check=check)

    def has_user_basis(self):
        """
        Return ``True`` if the basis of this free module is
        specified by the user, as opposed to being the default echelon
        form.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.subspace([[2,'1/2', 1]])
            sage: W.has_user_basis()
            False
            sage: W = V.subspace_with_basis([[2,'1/2',1]])
            sage: W.has_user_basis()
            True
        """
        return False


FreeModule_generic_field._submodule_class = FreeModule_submodule_field


###############################################################################

def element_class(R, is_sparse):
    """
    The class of the vectors (elements of a free module) with base ring
    R and boolean is_sparse.

    EXAMPLES::

        sage: FF = FiniteField(2)
        sage: P = PolynomialRing(FF,'x')
        sage: sage.modules.free_module.element_class(QQ, is_sparse=True)
        <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
        sage: sage.modules.free_module.element_class(QQ, is_sparse=False)
        <class 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        sage: sage.modules.free_module.element_class(ZZ, is_sparse=True)
        <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
        sage: sage.modules.free_module.element_class(ZZ, is_sparse=False)
        <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>
        sage: sage.modules.free_module.element_class(FF, is_sparse=True)
        <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
        sage: sage.modules.free_module.element_class(FF, is_sparse=False)               # needs sage.rings.finite_rings
        <class 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>
        sage: sage.modules.free_module.element_class(GF(7), is_sparse=False)
        <class 'sage.modules.vector_modn_dense.Vector_modn_dense'>
        sage: sage.modules.free_module.element_class(P, is_sparse=True)
        <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
        sage: sage.modules.free_module.element_class(P, is_sparse=False)
        <class 'sage.modules.free_module_element.FreeModuleElement_generic_dense'>
    """
    import sage.rings.integer_ring
    if isinstance(R, sage.rings.integer_ring.IntegerRing_class) and not is_sparse:
        from sage.modules.vector_integer_dense import Vector_integer_dense
        return Vector_integer_dense
    elif isinstance(R, sage.rings.rational_field.RationalField) and not is_sparse:
        from sage.modules.vector_rational_dense import Vector_rational_dense
        return Vector_rational_dense
    elif isinstance(R, sage.rings.abc.IntegerModRing) and not is_sparse:
        if R.order() == 2:
            try:
                from sage.modules.vector_mod2_dense import Vector_mod2_dense
            except ImportError:
                pass
            else:
                return Vector_mod2_dense
        try:
            from sage.modules.vector_modn_dense import MAX_MODULUS, Vector_modn_dense
        except ImportError:
            pass
        else:
            if R.order() < MAX_MODULUS:
                return Vector_modn_dense
        return free_module_element.FreeModuleElement_generic_dense
    elif isinstance(R, sage.rings.abc.RealDoubleField) and not is_sparse:
        try:
            from sage.modules.vector_real_double_dense import Vector_real_double_dense
        except ImportError:
            pass
        else:
            return Vector_real_double_dense
    elif isinstance(R, sage.rings.abc.ComplexDoubleField) and not is_sparse:
        try:
            from sage.modules.vector_complex_double_dense import (
                Vector_complex_double_dense,
            )
        except ImportError:
            pass
        else:
            return Vector_complex_double_dense
    elif isinstance(R, sage.rings.abc.CallableSymbolicExpressionRing) and not is_sparse:
        import sage.modules.vector_callable_symbolic_dense
        return sage.modules.vector_callable_symbolic_dense.Vector_callable_symbolic_dense
    elif isinstance(R, sage.rings.abc.SymbolicRing):
        if not is_sparse:
            import sage.modules.vector_symbolic_dense
            return sage.modules.vector_symbolic_dense.Vector_symbolic_dense
        else:
            import sage.modules.vector_symbolic_sparse
            return sage.modules.vector_symbolic_sparse.Vector_symbolic_sparse

    if is_sparse:
        return free_module_element.FreeModuleElement_generic_sparse
    else:
        return free_module_element.FreeModuleElement_generic_dense


@richcmp_method
class EchelonMatrixKey:
    r"""
    A total ordering on free modules for sorting.

    This class orders modules by their ambient spaces, then by dimension,
    then in order by their echelon matrices. If a function returns a list
    of free modules, this can be used to sort the output and thus render
    it deterministic.

    INPUT:

    - ``obj`` -- a free module

    EXAMPLES::

        sage: V = span([[1,2,3], [5,6,7], [8,9,10]], QQ)
        sage: W = span([[5,6,7], [8,9,10]], QQ)
        sage: X = span([[5,6,7]], ZZ).scale(1/11)
        sage: Y = CC^3
        sage: Z = ZZ^2
        sage: modules = [V,W,X,Y,Z]
        sage: modules_sorted = [Z,X,V,W,Y]
        sage: from sage.modules.free_module import EchelonMatrixKey
        sage: modules.sort(key=EchelonMatrixKey)
        sage: modules == modules_sorted
        True
    """
    def __init__(self, obj):
        r"""
        Create a container for a free module with a total ordering.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: V = span(R, [[x, 1 + x], [x^2, 2 + x]])
            sage: W = R^2
            sage: from sage.modules.free_module import EchelonMatrixKey
            sage: V = EchelonMatrixKey(V)
            sage: W = EchelonMatrixKey(W)
            sage: V < W
            False
        """
        self.obj = obj

    def __richcmp__(self, other, op):
        r"""
        A total ordering on free modules.

        TESTS::

            sage: from sage.modules.free_module import EchelonMatrixKey
            sage: Y = EchelonMatrixKey(CC^3)
            sage: Z = EchelonMatrixKey(ZZ^2)
            sage: Z < Y
            True
        """
        return self.obj._echelon_matrix_richcmp(other.obj, op)
