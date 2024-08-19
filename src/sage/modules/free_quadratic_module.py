r"""
Free quadratic modules

Sage supports computation with free quadratic modules over an arbitrary
commutative ring.  Nontrivial functionality is available over `\ZZ` and
fields.  All free modules over an integral domain are equipped with an
embedding in an ambient vector space and an inner product, which you
can specify and change.

Create the free module of rank `n` over an arbitrary commutative ring `R`
using the command ``FreeModule(R,n)`` with a given ``inner_product_matrix``.

The following example illustrates the creation of both a vector space
and a free module over the integers and a submodule of it.  Use the functions
``FreeModule``, ``span`` and member functions of free modules
to create free modules.  ''Do not use the ``FreeModule_xxx`` constructors
directly.''

EXAMPLES::

    sage: M = Matrix(QQ, [[2,1,0], [1,2,1], [0,1,2]])
    sage: V = VectorSpace(QQ, 3, inner_product_matrix=M)
    sage: type(V)
    <class 'sage.modules.free_quadratic_module.FreeQuadraticModule_ambient_field_with_category'>
    sage: V.inner_product_matrix()
    [2 1 0]
    [1 2 1]
    [0 1 2]
    sage: W = V.subspace([[1,2,7], [1,1,0]])
    sage: type(W)
    <class 'sage.modules.free_quadratic_module.FreeQuadraticModule_submodule_field_with_category'>
    sage: W
    Quadratic space of degree 3 and dimension 2 over Rational Field
    Basis matrix:
    [ 1  0 -7]
    [ 0  1  7]
    Inner product matrix:
    [2 1 0]
    [1 2 1]
    [0 1 2]
    sage: W.gram_matrix()
    [ 100 -104]
    [-104  114]

TESTS::

    sage: M = Matrix(QQ, [[2,1,0], [1,2,1], [0,1,2]])
    sage: V = VectorSpace(QQ, 3, inner_product_matrix=M)
    sage: V == loads(dumps(V))
    True
    sage: W = QuadraticSpace(QQ, 3, M)
    sage: W == V
    True

AUTHORS:

- David Kohel (2008-06): First created (based on free_module.py)
"""
# ****************************************************************************
#       Copyright (C) 2008 David Kohel <kohel@iml.univ-mrs.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import weakref

from sage.categories.commutative_rings import CommutativeRings
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.modules import free_module
from sage.rings.ring import Field, IntegralDomain
import sage.matrix.matrix_space
import sage.misc.latex as latex

# #############################################################################
#
# Constructor functions
#
# #############################################################################
_cache = {}


def FreeQuadraticModule(base_ring, rank, inner_product_matrix,
                        sparse=False, inner_product_ring=None):
    r"""
    Create the free quadratic module over the given commutative ring of the given rank.

    INPUT:

    - ``base_ring`` -- a commutative ring

    - ``rank`` -- nonnegative integer

    - ``inner_product_matrix`` -- the inner product matrix

    - ``sparse`` -- boolean (default: ``False``)

    - ``inner_product_ring`` -- the inner product codomain ring (default: ``None``)

    OUTPUT:

    A free quadratic module (with given inner product matrix).

    .. NOTE::

        In Sage, it is the case that there is only one dense and one
        sparse free ambient quadratic module of rank `n` over `R` and
        given inner product matrix.

    EXAMPLES::

        sage: M2 = FreeQuadraticModule(ZZ, 2, inner_product_matrix=[1,2,3,4])
        sage: M2 is FreeQuadraticModule(ZZ, 2, inner_product_matrix=[1,2,3,4])
        True
        sage: M2.inner_product_matrix()
        [1 2]
        [3 4]
        sage: M3 = FreeModule(ZZ, 2, inner_product_matrix=[[1,2],[3,4]])
        sage: M3 is M2
        True

    TESTS:

    Check for :issue:`10577`::

        sage: m = matrix.diagonal(GF(2), [1,1])
        sage: V2 = VectorSpace(GF(2), 2, inner_product_matrix=m)
        sage: deepcopy(V2)
        Ambient quadratic space of dimension 2 over Finite Field of size 2
        Inner product matrix:
        [1 0]
        [0 1]
    """
    global _cache
    rank = int(rank)

    # In order to use coercion into the inner_product_ring we need to pass
    # this ring into the vector classes.
    if inner_product_ring is not None:
        raise NotImplementedError("an inner_product_ring cannot currently be defined")

    # We intentionally create a new matrix instead of using the given
    # inner_product_matrix. This ensures that the matrix has the correct
    # parent space. It also gets rid of subdivisions which is good
    # because matrices with and without subdivisions compare equal.
    # Because of uniqueness, we need a canonical matrix, which is the one
    # without subdivisions.
    MS = sage.matrix.matrix_space.MatrixSpace(base_ring, rank)
    inner_product_matrix = MS(list(inner_product_matrix))
    inner_product_matrix.set_immutable()

    key = (base_ring, rank, inner_product_matrix, sparse)

    if key in _cache:
        M = _cache[key]()
        if M is not None:
            return M

    if base_ring not in CommutativeRings():
        raise TypeError("base_ring must be a commutative ring")

    # elif not sparse and isinstance(base_ring,sage.rings.real_double.RealDoubleField_class):
    #     M = RealDoubleQuadraticSpace_class(rank, inner_product_matrix=inner_product_matrix, sparse=False)

    # elif not sparse and isinstance(base_ring,sage.rings.complex_double.ComplexDoubleField_class):
    #     M = ComplexDoubleQuadraticSpace_class(rank, inner_product_matrix=inner_product_matrix, sparse=False)

    elif base_ring.is_field():
        M = FreeQuadraticModule_ambient_field(
            base_ring, rank, sparse=sparse, inner_product_matrix=inner_product_matrix)

    elif base_ring in PrincipalIdealDomains():
        M = FreeQuadraticModule_ambient_pid(
            base_ring, rank, sparse=sparse, inner_product_matrix=inner_product_matrix)

    elif isinstance(base_ring, IntegralDomain) or base_ring.is_integral_domain():
        M = FreeQuadraticModule_ambient_domain(
            base_ring, rank, sparse=sparse, inner_product_matrix=inner_product_matrix)
    else:
        M = FreeQuadraticModule_ambient(
            base_ring, rank, sparse=sparse, inner_product_matrix=inner_product_matrix)

    _cache[key] = weakref.ref(M)
    return M


def QuadraticSpace(K, dimension, inner_product_matrix, sparse=False):
    """
    EXAMPLES:

    The base can be complicated, as long as it is a field::

        sage: F.<x> = FractionField(PolynomialRing(ZZ,'x'))
        sage: D = diagonal_matrix([x, x - 1, x + 1])
        sage: V = QuadraticSpace(F, 3, D)
        sage: V
        Ambient quadratic space of dimension 3 over
         Fraction Field of Univariate Polynomial Ring in x over Integer Ring
        Inner product matrix:
        [    x     0     0]
        [    0 x - 1     0]
        [    0     0 x + 1]
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]

    The base must be a field or a :exc:`TypeError` is raised::

        sage: QuadraticSpace(ZZ, 5, identity_matrix(ZZ,2))
        Traceback (most recent call last):
        ...
        TypeError: argument K (= Integer Ring) must be a field
    """
    if not K.is_field():
        raise TypeError(f"argument K (= {K}) must be a field")
    if sparse not in (True, False):
        raise TypeError("Argument sparse (= %s) must be a boolean." % sparse)
    return FreeQuadraticModule(K, rank=dimension, inner_product_matrix=inner_product_matrix, sparse=sparse)


InnerProductSpace = QuadraticSpace


# #############################################################################
#
# Base class for all free modules
#
# #############################################################################

def is_FreeQuadraticModule(M):
    """
    Return ``True`` if `M` is a free quadratic module.

    EXAMPLES::

        sage: from sage.modules.free_quadratic_module import is_FreeQuadraticModule
        sage: U = FreeModule(QQ,3)
        sage: is_FreeQuadraticModule(U)
        doctest:warning...
        DeprecationWarning: the function is_FreeQuadraticModule is deprecated;
        use 'isinstance(..., FreeQuadraticModule_generic)' instead
        See https://github.com/sagemath/sage/issues/37924 for details.
        False
        sage: V = FreeModule(QQ,3,inner_product_matrix=diagonal_matrix([1,1,1]))
        sage: is_FreeQuadraticModule(V)
        True
        sage: W = FreeModule(QQ,3,inner_product_matrix=diagonal_matrix([2,3,3]))
        sage: is_FreeQuadraticModule(W)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37924, "the function is_FreeQuadraticModule is deprecated; use 'isinstance(..., FreeQuadraticModule_generic)' instead")
    return isinstance(M, FreeQuadraticModule_generic)


class FreeQuadraticModule_generic(free_module.FreeModule_generic):
    """
    Base class for all free quadratic modules.

    Modules are ordered by inclusion in the same ambient space.

    TESTS:

    We compare rank three free modules over the integers,
    rationals, and complex numbers::

        sage: Q3 = FreeQuadraticModule(QQ,3,matrix.identity(3))
        sage: C3 = FreeQuadraticModule(CC,3,matrix.identity(3))
        sage: Z3 = FreeQuadraticModule(ZZ,3,matrix.identity(3))
        sage: Q3 < C3
        False
        sage: C3 < Q3
        False
        sage: C3 > Q3
        False
        sage: Q3 > Z3
        True
        sage: Q3 < Z3
        False
        sage: Z3 < Q3
        True
        sage: Z3 > Q3
        False
        sage: Q3 == Z3
        False
        sage: Q3 == Q3
        True

        sage: V = Q3.span([[1,2,3], [5,6,7], [8,9,10]])
        sage: V < Q3
        True
        sage: Q3 < V
        False

    The :meth:`inner_product_matrix` is part of the comparison::

        sage: Q3zero = FreeQuadraticModule(QQ,3,matrix.zero(3))
        sage: Q3zero == Q3
        False

    We test that :issue:`23915` is fixed::

        sage: M1 = FreeQuadraticModule(ZZ,1,matrix.identity(1))
        sage: M2 = FreeQuadraticModule(ZZ,1,matrix.identity(1)*2)
        sage: M1 == M2
        False
    """
    def __init__(self, base_ring, rank, degree, inner_product_matrix, sparse=False):
        """
        Create the free module of given rank over the given ``base_ring``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``rank`` -- nonnegative integer

        EXAMPLES::

            sage: R = PolynomialRing(QQ,3,'x')
            sage: FreeModule(R,3,inner_product_matrix=diagonal_matrix(list(R.gens())))
            Ambient free quadratic module of rank 3 over the integral domain Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
            Inner product matrix:
            [x0  0  0]
            [ 0 x1  0]
            [ 0  0 x2]
        """
        free_module.FreeModule_generic.__init__(
            self, base_ring=base_ring, rank=rank, degree=degree, sparse=sparse)
        self._inner_product_matrix = inner_product_matrix

    def _dense_module(self):
        """
        Create a dense module with the same defining data as ``self``.

        ..  NOTE:: This function is for internal use only! See ``dense_module`` for use.

        EXAMPLES::

            sage: A = diagonal_matrix([1,2,2])
            sage: M = FreeModule(Integers(8),3,inner_product_matrix=A)
            sage: S = FreeModule(Integers(8),3,inner_product_matrix=A,sparse=True)
            sage: M is S._dense_module()
            True
        """
        A = self.ambient_module().dense_module()
        return A.span(self.basis())

    def _sparse_module(self):
        """
        Create a sparse module with the same defining data as ``self``.

        .. NOTE:: This function is for internal use only! See ``sparse_module`` for use.

        EXAMPLES::

            sage: A = diagonal_matrix([1,2,2])
            sage: M = FreeModule(Integers(8),3,inner_product_matrix=A)
            sage: S = FreeModule(Integers(8),3,inner_product_matrix=A,sparse=True)
            sage: M._sparse_module() is S
            True
        """
        A = self.ambient_module().sparse_module()
        return A.span(self.basis())

    def ambient_module(self):
        """
        Return the ambient module associated to this module.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: M = FreeModule(R,2)
            sage: M.ambient_module()
            Ambient free module of rank 2 over the integral domain Multivariate Polynomial Ring in x, y over Rational Field

            sage: V = FreeModule(QQ, 4).span([[1,2,3,4], [1,0,0,0]]); V
            Vector space of degree 4 and dimension 2 over Rational Field
            Basis matrix:
            [  1   0   0   0]
            [  0   1 3/2   2]
            sage: V.ambient_module()
            Vector space of dimension 4 over Rational Field
        """
        return FreeQuadraticModule(self.base_ring(), self.degree(), self.inner_product_matrix())

    def determinant(self):
        """
        Return the determinant of this free module.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3, inner_product_matrix=1)
            sage: M.determinant()
            1
            sage: N = M.span([[1,2,3]])
            sage: N.determinant()
            14
            sage: P = M.span([[1,2,3], [1,1,1]])
            sage: P.determinant()
            6
        """
        return self.gram_matrix().determinant()

    def discriminant(self):
        """
        Return the discriminant of this free module.

        This is defined to be `(-1)^r` of the determinant, where `r = n/2`
        (`n` even) or `(n-1)/2` (`n` odd) for a module of rank `n`.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3)
            sage: M.discriminant()
            1
            sage: N = M.span([[1,2,3]])
            sage: N.discriminant()
            14
            sage: P = M.span([[1,2,3], [1,1,1]])
            sage: P.discriminant()
            6

        TESTS::

            sage: M = FreeQuadraticModule(ZZ, 2, matrix.identity(2))
            sage: M.discriminant()
            -1
            sage: M = FreeQuadraticModule(QQ, 3, matrix.identity(3))
            sage: M.discriminant()
            -1
        """
        r = self.rank() // 2
        return (-1)**r * self.gram_matrix().determinant()

    def gram_matrix(self):
        """
        Return the Gram matrix associated to this free module.

        This is defined to be ``B*A*B.transpose()``, where ``A`` is the
        inner product matrix (induced from the ambient space), and ``B``
        the basis matrix.

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
            return self.inner_product_matrix()
        if self._gram_matrix is None:
            A = self.inner_product_matrix()
            B = self.basis_matrix()
            self._gram_matrix = B * A * B.transpose()
        return self._gram_matrix

    def inner_product_matrix(self):
        """
        Return the inner product matrix associated to this module.

        By definition, this is the inner product matrix of the ambient
        space, hence may be of degree greater than the rank of the
        module.

        .. NOTE:: The inner product does not have to be symmetric (see examples).

        .. TODO::

            Differentiate the image ring of the inner product from the base ring of
            the module and/or ambient space.  E.g. On an integral module over ZZ the inner
            product pairing could naturally take values in ZZ, QQ, RR, or CC.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 3)
            sage: M.inner_product_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

        The inner product does not have to be symmetric or definite::

            sage: N = FreeModule(ZZ,2,inner_product_matrix=[[1,-1],[2,5]])
            sage: N.inner_product_matrix()
            [ 1 -1]
            [ 2  5]
            sage: u, v = N.basis()
            sage: u.inner_product(v)
            -1
            sage: v.inner_product(u)
            2

        The inner product matrix is defined with respect to the ambient space::

            sage: V = QQ^3
            sage: u = V([1/2,1,1])
            sage: v = V([1,1,1/2])
            sage: M = span([u,v], ZZ)
            sage: M.inner_product_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: M.inner_product_matrix() == V.inner_product_matrix()
            True
            sage: M.gram_matrix()
            [ 1/2 -3/4]
            [-3/4 13/4]
        """
        return self._inner_product_matrix

    def _inner_product_is_dot_product(self):
        """
        Return whether or not the inner product on this module is induced by
        the dot product on the ambient vector space.

        This is used internally by the ``inner_product`` function for
        optimization.

        EXAMPLES::

            sage: FreeModule(ZZ, 3)._inner_product_is_dot_product()
            True
            sage: FreeModule(ZZ, 3, inner_product_matrix=1)._inner_product_is_dot_product()
            True
            sage: FreeModule(ZZ, 2, inner_product_matrix=[1,0,-1,0])._inner_product_is_dot_product()
            False

            sage: M = FreeModule(QQ, 3)
            sage: M2 = M.span([[1,2,3]])
            sage: M2._inner_product_is_dot_product()
            True
        """
        return self.inner_product_matrix() == 1

    def _inner_product_is_diagonal(self):
        """
        Return whether or not the inner product on this module is induced by
        the dot product on the ambient vector space.

        This is used internally by the ``inner_product`` function for
        optimization.

        .. NOTE::

            The ``FreeModule`` classes have the identity inner product matrix,
            while ``FreeQuadraticModules`` must have an ``inner_product_matrix``, although
            it can be diagonal.

        EXAMPLES::

            sage: M0 = FreeModule(ZZ, 3, inner_product_matrix=1)
            sage: M0._inner_product_is_diagonal()
            True
            sage: D = diagonal_matrix([3,5,7])
            sage: M1 = FreeModule(ZZ, 3, inner_product_matrix=D)
            sage: M1._inner_product_is_diagonal()
            True
            sage: A = Matrix([[2,1,0],[1,2,1],[0,1,2]])
            sage: M2 = FreeModule(ZZ, 3, inner_product_matrix=A)
            sage: M2._inner_product_is_diagonal()
            False
            sage: M3 = FreeModule(ZZ, 2, inner_product_matrix=[1,0,-1,0])
            sage: M3._inner_product_is_diagonal()
            False

        .. TODO:: Actually use the diagonal form of the inner product.
        """
        A = self.inner_product_matrix()
        D = sage.matrix.constructor.diagonal_matrix(A.diagonal())
        return A == D


class FreeQuadraticModule_generic_pid(free_module.FreeModule_generic_pid,
                                      FreeQuadraticModule_generic):
    """
    Class of all free modules over a PID.
    """
    def __init__(self, base_ring, rank, degree, inner_product_matrix, sparse=False):
        """
        Create a free module over a PID.

        EXAMPLES::

            sage: FreeModule(ZZ, 2, inner_product_matrix=Matrix([[2,1],[1,2]]))
            Ambient free quadratic module of rank 2 over the principal ideal domain Integer Ring
            Inner product matrix:
            [2 1]
            [1 2]
        """
        free_module.FreeModule_generic_pid.__init__(
            self, base_ring=base_ring, rank=rank, degree=degree, sparse=sparse)
        self._inner_product_matrix = inner_product_matrix

    def span(self, gens, check=True, already_echelonized=False):
        """
        Return the `R`-span of the given list of gens, where `R`
        is the base ring of ``self``.

        Note that this span need not be a submodule of ``self``, nor even
        of the ambient space.  It must, however, be contained in the
        ambient vector space, i.e., the ambient space tensored with
        the fraction field of `R`.

        EXAMPLES::

            sage: V = FreeModule(ZZ,3)
            sage: W = V.submodule([V.gen(0)])
            sage: W.span([V.gen(1)])
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1 0]
            sage: W.submodule([V.gen(1)])
            Traceback (most recent call last):
            ...
            ArithmeticError: argument gens (= [(0, 1, 0)]) does not generate a submodule of self
        """
        return FreeQuadraticModule_submodule_pid(
            self.ambient_module(), gens, inner_product_matrix=self.inner_product_matrix(),
            check=check, already_echelonized=already_echelonized)

    def span_of_basis(self, basis, check=True, already_echelonized=False):
        r"""
        Return the free `R`-module with the given basis, where `R`
        is the base ring of ``self``.

        Note that this `R`-module need not be a submodule of ``self``, nor
        even of the ambient space.  It must, however, be contained in
        the ambient vector space, i.e., the ambient space tensored
        with the fraction field of `R`.

        EXAMPLES::

            sage: M = FreeModule(ZZ,3)
            sage: W = M.span_of_basis([M([1,2,3])])

        Next we create two free `\ZZ`-modules, neither of which is a
        submodule of `W`::

            sage: W.span_of_basis([M([2,4,0])])
            Free module of degree 3 and rank 1 over Integer Ring
            User basis matrix:
            [2 4 0]

        The following module is not even in the ambient space::

            sage: Q = QQ
            sage: W.span_of_basis([ Q('1/5')*M([1,2,0]), Q('1/7')*M([1,1,0]) ])
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
        return FreeQuadraticModule_submodule_with_basis_pid(
            self.ambient_module(), basis=basis, inner_product_matrix=self.inner_product_matrix(),
            check=check, already_echelonized=already_echelonized)

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
        return FreeQuadraticModule_submodule_pid(
            self.ambient_module(), [], self.inner_product_matrix(), check=False)


class FreeQuadraticModule_generic_field(free_module.FreeModule_generic_field,
                                        FreeQuadraticModule_generic_pid):
    """
    Base class for all free modules over fields.
    """
    def __init__(self, base_field, dimension, degree, inner_product_matrix, sparse=False):
        """
        Create a vector space over a field.

        EXAMPLES::

            sage: FreeModule(QQ, 2, inner_product_matrix=[[2,1],[1,2]])
            Ambient quadratic space of dimension 2 over Rational Field
            Inner product matrix:
            [2 1]
            [1 2]
            sage: FreeModule(FiniteField(2), 7, inner_product_matrix=1)
            Ambient quadratic space of dimension 7 over Finite Field of size 2
            Inner product matrix:
            [1 0 0 0 0 0 0]
            [0 1 0 0 0 0 0]
            [0 0 1 0 0 0 0]
            [0 0 0 1 0 0 0]
            [0 0 0 0 1 0 0]
            [0 0 0 0 0 1 0]
            [0 0 0 0 0 0 1]
        """
        if not isinstance(base_field, Field):
            raise TypeError("the base_field (=%s) must be a field" % base_field)
        free_module.FreeModule_generic_field.__init__(
            self, base_field=base_field, dimension=dimension, degree=degree, sparse=sparse)
        self._inner_product_matrix = inner_product_matrix

    def span(self, gens, check=True, already_echelonized=False):
        """
        Return the `K`-span of the given list of gens, where `K` is the
        base field of ``self``.

        Note that this span is a subspace of the ambient vector space,
        but need not be a subspace of ``self``.

        INPUT:

        - ``gens`` -- list of vectors

        - ``check`` -- boolean (default: ``True``); whether or not to coerce
          entries of gens into base field

        - ``already_echelonized`` -- boolean (default: ``False``); set this if
          you know the gens are already in echelon form

        EXAMPLES::

            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace([[2,3,4]]); W
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 5 2]
            sage: W.span([[1,1,1]])
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 1 1]
        """
        if isinstance(gens, free_module.FreeModule_generic):
            gens = gens.gens()
        if not isinstance(gens, (list, tuple)):
            raise TypeError("gens (=%s) must be a list or tuple" % gens)

        return FreeQuadraticModule_submodule_field(
            self.ambient_module(), gens,
            inner_product_matrix=self.inner_product_matrix(),
            check=check, already_echelonized=already_echelonized)

    def span_of_basis(self, basis, check=True, already_echelonized=False):
        r"""
        Return the free `K`-module with the given basis, where `K`
        is the base field of ``self``.

        Note that this span is a subspace of the ambient vector space,
        but need not be a subspace of ``self``.

        INPUT:

        - ``basis`` -- list of vectors

        - ``check`` -- boolean (default: ``True``); whether or not to coerce
          entries of gens into base field

        - ``already_echelonized`` -- boolean (default: ``False``); set this if
          you know the gens are already in echelon form

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
        return FreeQuadraticModule_submodule_with_basis_field(
            self.ambient_module(), basis=basis,
            inner_product_matrix=self.inner_product_matrix(),
            check=check, already_echelonized=already_echelonized)


# #############################################################################
#
# Generic ambient free modules, i.e., of the form R^n for some commutative ring R.
#
# #############################################################################

class FreeQuadraticModule_ambient(free_module.FreeModule_ambient,
                                  FreeQuadraticModule_generic):
    """
    Ambient free module over a commutative ring.
    """
    def __init__(self, base_ring, rank, inner_product_matrix, sparse=False):
        """
        The free module of given rank over the given ``base_ring``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``rank`` -- nonnegative integer

        EXAMPLES::

            sage: FreeModule(ZZ, 4)
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
        """
        free_module.FreeModule_ambient.__init__(self, base_ring=base_ring, rank=rank, sparse=sparse)
        self._inner_product_matrix = inner_product_matrix

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

        The system representation can be overwritten, but leaves
        :meth:`_repr_` unmodified::

            sage: M.rename('M')
            sage: M
            M
            sage: print(M._repr_())
            Ambient free module of rank 12 over Ring of integers modulo 12

        Sparse modules print this fact::

            sage: N = FreeModule(R,12,sparse=True)
            sage: N
            Ambient sparse free module of rank 12 over Ring of integers modulo 12
        """
        if self.is_sparse():
            return "Ambient sparse free quadratic module of rank %s over %s\n" % (self.rank(), self.base_ring()) + \
                "Inner product matrix:\n%s" % self.inner_product_matrix()
        return "Ambient free quadratic module of rank %s over %s\n" % (self.rank(), self.base_ring()) + \
            "Inner product matrix:\n%s" % self.inner_product_matrix()

    def _latex_(self):
        r"""
        Return a latex representation of this ambient free quadratic module.

        EXAMPLES::

            sage: latex(QQ^3) # indirect doctest
            \Bold{Q}^{3}

            sage: A = GF(5)^20; latex(A)
            \Bold{F}_{5}^{20}

            sage: A = PolynomialRing(QQ,3,'x')^20; latex(A)
            (\Bold{Q}[x_{0}, x_{1}, x_{2}])^{20}

            sage: V = QuadraticSpace(QQ,3,inner_product_matrix=[[2,1,0],[1,4,1],[0,1,8]])
            sage: latex(V)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # How do we want to represent this object?
        raise NotImplementedError

    def _dense_module(self):
        """
        Create a dense module with the same defining data as ``self``.

        .. NOTE:: This function is for internal use only! See dense_module for use.

        EXAMPLES::

            sage: A = diagonal_matrix([1,2,2])
            sage: M = FreeModule(Integers(8),3,inner_product_matrix=A)
            sage: S = FreeModule(Integers(8),3,inner_product_matrix=A,sparse=True)
            sage: M is S._dense_module()
            True
        """
        return FreeQuadraticModule(base_ring=self.base_ring(), rank=self.rank(),
                                   inner_product_matrix=self.inner_product_matrix(),
                                   sparse=False)

    def _sparse_module(self):
        """
        Create a sparse module with the same defining data as ``self``.

        .. NOTE:: This function is for internal use only! See sparse_module for use.

        EXAMPLES::

            sage: A = diagonal_matrix([1,2,2])
            sage: M = FreeModule(Integers(8),3,inner_product_matrix=A)
            sage: S = FreeModule(Integers(8),3,inner_product_matrix=A,sparse=True)
            sage: M._sparse_module() is S
            True
        """
        return FreeQuadraticModule(base_ring=self.base_ring(), rank=self.rank(),
                                   inner_product_matrix=self.inner_product_matrix(),
                                   sparse=True)


# #############################################################################
#
# Ambient free modules over an integral domain.
#
# #############################################################################

class FreeQuadraticModule_ambient_domain(free_module.FreeModule_ambient_domain,
                                         FreeQuadraticModule_ambient):
    """
    Ambient free quadratic module over an integral domain.
    """
    def __init__(self, base_ring, rank, inner_product_matrix, sparse=False):
        """
        EXAMPLES::

            sage: FreeModule(PolynomialRing(GF(5),'x'), 3)
            Ambient free module of rank 3 over the principal ideal domain
            Univariate Polynomial Ring in x over Finite Field of size 5
        """
        free_module.FreeModule_ambient.__init__(self, base_ring=base_ring, rank=rank, sparse=sparse)
        self._inner_product_matrix = inner_product_matrix

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ,'x')
            sage: M = FreeModule(R,7)
            sage: M
            Ambient free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring
            sage: print(M._repr_())
            Ambient free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring

        The system representation can be overwritten, but leaves
        :meth:`_repr_` unmodified::

            sage: M.rename('M')
            sage: M
            M
            sage: print(M._repr_())
            Ambient free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring

        Sparse modules print this fact::

            sage: N = FreeModule(R,7,sparse=True)
            sage: N
            Ambient sparse free module of rank 7 over the integral domain Univariate Polynomial Ring in x over Integer Ring

        Here is a construction of a free quadratic module with generic
        symmetric inner product matrix::

            sage: R.<a,b,c> = PolynomialRing(QQ,3)
            sage: M = FreeModule(R, 2, inner_product_matrix=[[2*a,b],[b,2*c]])
            sage: M
            Ambient free quadratic module of rank 2 over the integral domain Multivariate Polynomial Ring in a, b, c over Rational Field
            Inner product matrix:
            [2*a   b]
            [  b 2*c]
            sage: M.determinant()
            -b^2 + 4*a*c
        """
        if self.is_sparse():
            return "Ambient sparse free quadratic module of rank %s over the integral domain %s\n" % (
                self.rank(), self.base_ring()) + \
                "Inner product matrix:\n%s" % self.inner_product_matrix()
        return "Ambient free quadratic module of rank %s over the integral domain %s\n" % (
            self.rank(), self.base_ring()) + \
            "Inner product matrix:\n%s" % self.inner_product_matrix()

    def ambient_vector_space(self):
        """
        Return the ambient vector space, which is this free module tensored
        with its fraction field.

        EXAMPLES::

            sage: M = ZZ^3;  M.ambient_vector_space()
            Vector space of dimension 3 over Rational Field
        """
        try:
            return self.__ambient_vector_space
        except AttributeError:
            self.__ambient_vector_space = FreeQuadraticModule(
                self.base_field(), self.rank(),
                inner_product_matrix=self.inner_product_matrix(), sparse=self.is_sparse())
            return self.__ambient_vector_space


# #############################################################################
#
# Ambient free modules over a principal ideal domain.
#
# #############################################################################

class FreeQuadraticModule_ambient_pid(free_module.FreeModule_ambient_pid,
                                      FreeQuadraticModule_generic_pid,
                                      FreeQuadraticModule_ambient_domain):
    """
    Ambient free quadratic module over a principal ideal domain.
    """
    def __init__(self, base_ring, rank, inner_product_matrix, sparse=False):
        """
        Create the ambient free module of given rank over the given
        principal ideal domain.

        INPUT:

        - ``base_ring`` -- a principal ideal domain

        - ``rank`` -- nonnegative integer

        - ``sparse`` -- boolean (default: ``False``)

        - ``inner_product_matrix`` -- boolean (default: ``None``)

        EXAMPLES::

            sage: ZZ^3
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: FreeModule(ZZ,3,inner_product_matrix=Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]]))
            Ambient free quadratic module of rank 3 over the principal ideal domain Integer Ring
            Inner product matrix:
            [ 2 -1  0]
            [-1  2 -1]
            [ 0 -1  2]
        """
        free_module.FreeModule_ambient_pid.__init__(self, base_ring=base_ring, rank=rank, sparse=sparse)
        self._inner_product_matrix = inner_product_matrix

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: M = FreeModule(ZZ, 2, inner_product_matrix=[[2,-1],[-1,2]])
            sage: M
            Ambient free quadratic module of rank 2 over the principal ideal domain Integer Ring
            Inner product matrix:
            [ 2 -1]
            [-1  2]

        Without a user specified inner product the class and printing
        is simpler::

            sage: M = FreeModule(ZZ,7)
            sage: M
            Ambient free module of rank 7 over the principal ideal domain Integer Ring
            sage: print(M._repr_())
            Ambient free module of rank 7 over the principal ideal domain Integer Ring

        The system representation can be overwritten, but leaves
        :meth:`_repr_` unmodified::

            sage: M.rename('M')
            sage: M
            M
            sage: print(M._repr_())
            Ambient free module of rank 7 over the principal ideal domain Integer Ring

        Sparse modules print this fact::

            sage: N = FreeModule(ZZ,7,sparse=True)
            sage: N
            Ambient sparse free module of rank 7 over the principal ideal domain Integer Ring
        """
        if self.is_sparse():
            return "Ambient sparse free quadratic module of rank %s over the principal ideal domain %s\n" % (
                self.rank(), self.base_ring()) + \
                "Inner product matrix:\n%s" % self.inner_product_matrix()
        return "Ambient free quadratic module of rank %s over the principal ideal domain %s\n" % (
            self.rank(), self.base_ring()) + \
            "Inner product matrix:\n%s" % self.inner_product_matrix()


# #############################################################################
#
# Ambient free modules over a field (i.e., a vector space).
#
# #############################################################################

class FreeQuadraticModule_ambient_field(free_module.FreeModule_ambient_field,
                                        FreeQuadraticModule_generic_field,
                                        FreeQuadraticModule_ambient_pid):

    def __init__(self, base_field, dimension, inner_product_matrix, sparse=False):
        """
        Create the ambient vector space of given dimension over the given field.

        INPUT:

        - ``base_field`` -- a field

        - ``dimension`` -- nonnegative integer

        - ``sparse`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: VectorSpace(QQ,3,inner_product_matrix=[[2,1,0],[1,2,0],[0,1,2]])
            Ambient quadratic space of dimension 3 over Rational Field
            Inner product matrix:
            [2 1 0]
            [1 2 0]
            [0 1 2]

        TESTS:

        Check for :issue:`10606`::

            sage: D = matrix.diagonal(ZZ, [1,1])
            sage: V = VectorSpace(GF(46349), 2, inner_product_matrix=D)                 # needs sage.rings.finite_rings
            sage: deepcopy(V)                                                           # needs sage.rings.finite_rings
            Ambient quadratic space of dimension 2 over Finite Field
            of size 46349
            Inner product matrix:
            [1 0]
            [0 1]
        """
        free_module.FreeModule_ambient_field.__init__(
            self, base_field=base_field, dimension=dimension, sparse=sparse)
        self._inner_product_matrix = inner_product_matrix

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: V = FreeModule(QQ,7)
            sage: V
            Vector space of dimension 7 over Rational Field
            sage: print(V._repr_())
            Vector space of dimension 7 over Rational Field

        The system representation can be overwritten, but leaves
        :meth:`_repr_` unmodified::

            sage: V.rename('V')
            sage: V
            V
            sage: print(V._repr_())
            Vector space of dimension 7 over Rational Field

        Sparse modules print this fact::

            sage: U = FreeModule(QQ,7,sparse=True)
            sage: U
            Sparse vector space of dimension 7 over Rational Field
        """
        if self.is_sparse():
            return "Ambient sparse free quadratic space of dimension %s over %s\n" % (self.rank(), self.base_ring()) + \
                "Inner product matrix:\n%s" % self.inner_product_matrix()
        return "Ambient quadratic space of dimension %s over %s\n" % (self.rank(), self.base_ring()) + \
            "Inner product matrix:\n%s" % self.inner_product_matrix()


# #############################################################################
#
# R-Submodule of K^n where K is the fraction field of a principal ideal domain R.
#
# #############################################################################


class FreeQuadraticModule_submodule_with_basis_pid(free_module.FreeModule_submodule_with_basis_pid,
                                                   FreeQuadraticModule_generic_pid):
    r"""
    An `R`-submodule of `K^n` with distinguished basis, where `K` is
    the fraction field of a principal ideal domain `R`.

    Modules are ordered by inclusion.

    EXAMPLES:

    First we compare two equal vector spaces::

        sage: A = FreeQuadraticModule(QQ,3,2*matrix.identity(3))
        sage: V = A.span([[1,2,3], [5,6,7], [8,9,10]])
        sage: W = A.span([[5,6,7], [8,9,10]])
        sage: V == W
        True

    Next we compare a one dimensional space to the two dimensional
    space defined above::

        sage: M = A.span([[5,6,7]])
        sage: V == M
        False
        sage: M < V
        True
        sage: V < M
        False

    We compare a `\ZZ`-module to the one-dimensional space above::

        sage: V = A.span([[5,6,7]])
        sage: V = V.change_ring(ZZ).scale(1/11)
        sage: V < M
        True
        sage: M < V
        False
    """
    def __init__(self, ambient, basis, inner_product_matrix,
                 check=True, echelonize=False, echelonized_basis=None,
                 already_echelonized=False):
        """
        Create a free module with basis over a PID.

        EXAMPLES::

            sage: A = diagonal_matrix([1,2,2])
            sage: M = FreeQuadraticModule(ZZ,3,inner_product_matrix=A)
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]]); W
            Free quadratic module of degree 3 and rank 2 over Integer Ring
            Basis matrix:
            [1 2 3]
            [4 5 6]
            Inner product matrix:
            [1 0 0]
            [0 2 0]
            [0 0 2]

            sage: W = M.span_of_basis([[1,2,3/2],[4,5,6]]); W
            Free quadratic module of degree 3 and rank 2 over Integer Ring
            Basis matrix:
            [  1   2 3/2]
            [  4   5   6]
            Inner product matrix:
            [1 0 0]
            [0 2 0]
            [0 0 2]

        TESTS:

        We test that :issue:`23703` is fixed::

            sage: A = FreeQuadraticModule(ZZ, 1, matrix.identity(1))
            sage: B = A.span([[1/2]])
            sage: C = B.span([[1]])
            sage: B.intersection(C) == C.intersection(B)
            True
        """
        free_module.FreeModule_submodule_with_basis_pid.__init__(
            self, ambient=ambient, basis=basis, check=check,
            echelonize=echelonize, echelonized_basis=echelonized_basis, already_echelonized=already_echelonized)
        self._inner_product_matrix = inner_product_matrix

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
            s = "Sparse free quadratic module of degree %s and rank %s over %s\n" % (
                self.degree(), self.rank(), self.base_ring()) + \
                "Basis matrix:\n%r\n" % self.basis_matrix() + \
                "Inner product matrix:\n%r" % self.inner_product_matrix()
        else:
            s = "Free quadratic module of degree %s and rank %s over %s\n" % (
                self.degree(), self.rank(), self.base_ring()) + \
                "Basis matrix:\n%r\n" % self.basis_matrix() + \
                "Inner product matrix:\n%r" % self.inner_product_matrix()
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
        return "\\mathrm{RowSpan}_{%s}%s" % (latex.latex(self.base_ring()),
                                             latex.latex(self.basis_matrix()))

    def change_ring(self, R):
        """
        Return the free module over `R` obtained by coercing each
        element of ``self`` into a vector over the fraction field of `R`,
        then taking the resulting `R`-module.

        This raises a :exc:`TypeError` if coercion is not possible.

        INPUT:

        - ``R`` -- a principal ideal domain

        EXAMPLES:

        Changing rings preserves the inner product and the user basis::

            sage: V = QQ^3
            sage: W = V.subspace([[2, '1/2', 1]]); W
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [  1 1/4 1/2]
            sage: W.change_ring(GF(7))
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 2 4]

            sage: N = FreeModule(ZZ, 2, inner_product_matrix=[[1,-1], [2,5]])
            sage: N.inner_product_matrix()
            [ 1 -1]
            [ 2  5]
            sage: Np = N.change_ring(RDF)
            sage: Np.inner_product_matrix()
            [ 1.0 -1.0]
            [ 2.0  5.0]
        """
        if self.base_ring() is R:
            return self
        K = R.fraction_field()
        A = self.inner_product_matrix()
        V = QuadraticSpace(K, self.degree(), inner_product_matrix=A)
        B = [V(b) for b in self.basis()]
        M = FreeQuadraticModule(R, self.degree(), inner_product_matrix=A)
        if self.has_user_basis():
            return M.span_of_basis(B)
        return M.span(B)


class FreeQuadraticModule_submodule_pid(free_module.FreeModule_submodule_pid,
                                        FreeQuadraticModule_submodule_with_basis_pid):
    """
    An `R`-submodule of `K^n` where `K` is the fraction field of a
    principal ideal domain `R`.

    EXAMPLES::

        sage: M = ZZ^3
        sage: W = M.span_of_basis([[1,2,3], [4,5,19]]); W
        Free module of degree 3 and rank 2 over Integer Ring
        User basis matrix:
        [ 1  2  3]
        [ 4  5 19]

    We can save and load submodules and elements::

        sage: loads(W.dumps()) == W
        True
        sage: v = W.0 + W.1
        sage: loads(v.dumps()) == v
        True
    """
    def __init__(self, ambient, gens, inner_product_matrix, check=True, already_echelonized=False):
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
        free_module.FreeModule_submodule_pid.__init__(
            self, ambient=ambient, gens=gens, check=check, already_echelonized=already_echelonized)
        self._inner_product_matrix = inner_product_matrix

    def _repr_(self):
        """
        The printing representation of ``self``.

        EXAMPLES::

            sage: M = FreeModule(ZZ,8,inner_product_matrix=1)
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


class FreeQuadraticModule_submodule_with_basis_field(free_module.FreeModule_submodule_with_basis_field,
                                                     FreeQuadraticModule_generic_field,
                                                     FreeQuadraticModule_submodule_with_basis_pid):
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
    ``echelon_coordinates()`` and user ``coordinates()`` do not agree::

        sage: V = QQ^3
        sage: W = V.submodule_with_basis([[1,2,3], [4,5,6]])
        sage: W
        Vector space of degree 3 and dimension 2 over Rational Field
        User basis matrix:
        [1 2 3]
        [4 5 6]

        sage: v = V([1,5,9])
        sage: W.echelon_coordinates(v)
        [1, 5]
        sage: vector(QQ, W.echelon_coordinates(v)) * W.echelonized_basis_matrix()
        (1, 5, 9)

        sage: v = V([1,5,9])
        sage: W.coordinates(v)
        [5, -1]
        sage: vector(QQ, W.coordinates(v)) * W.basis_matrix()
        (1, 5, 9)

    We can load and save submodules::

        sage: loads(W.dumps()) == W
        True

        sage: K.<x> = FractionField(PolynomialRing(QQ,'x'))
        sage: M = K^3; W = M.span_of_basis([[1,1,x]])
        sage: loads(W.dumps()) == W
        True
    """
    def __init__(self, ambient, basis, inner_product_matrix,
                 check=True, echelonize=False, echelonized_basis=None, already_echelonized=False):
        """
        Create a vector space with given basis.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.span_of_basis([[1,2,3], [4,5,6]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 2 3]
            [4 5 6]
            sage: V = VectorSpace(QQ, 3, inner_product_matrix=1)
            sage: V.span_of_basis([[1,2,3], [4,5,6]])
            Quadratic space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 2 3]
            [4 5 6]
            Inner product matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        free_module.FreeModule_submodule_with_basis_field.__init__(
            self, ambient=ambient, basis=basis, check=check,
            echelonize=echelonize, echelonized_basis=echelonized_basis, already_echelonized=already_echelonized)
        self._inner_product_matrix = inner_product_matrix

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

        The system representation can be overwritten, but leaves
        :meth:`_repr_` unmodified::

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

        Sparse vector spaces print this fact::

            sage: V = VectorSpace(QQ,5,sparse=True)
            sage: U = V.submodule([ V.gen(i) - V.gen(0) for i in range(1,5) ])
            sage: U # indirect doctest
            Sparse vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
        """
        if self.is_sparse():
            return "Sparse quadratic space of degree %s and dimension %s over %s\n" % (
                self.degree(), self.dimension(), self.base_field()) + \
                "Basis matrix:\n%r" % self.basis_matrix() + \
                "Inner product matrix:\n%r" % self.inner_product_matrix()
        return "Quadratic space of degree %s and dimension %s over %s\n" % (
            self.degree(), self.dimension(), self.base_field()) + \
            "Basis matrix:\n%r\n" % self.basis_matrix() + \
            "Inner product matrix:\n%r" % self.inner_product_matrix()


class FreeQuadraticModule_submodule_field(free_module.FreeModule_submodule_field,
                                          FreeQuadraticModule_submodule_with_basis_field):
    """
    An embedded vector subspace with echelonized basis.

    EXAMPLES:

    Since this is an embedded vector subspace with echelonized basis,
    the methods :meth:`echelon_coordinates` and :meth:`coordinates` return the same
    coordinates::

        sage: V = QQ^3
        sage: W = V.span([[1,2,3], [4,5,6]])
        sage: W
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -1]
        [ 0  1  2]

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
    def __init__(self, ambient, gens, inner_product_matrix, check=True, already_echelonized=False):
        """
        Create an embedded vector subspace with echelonized basis.

        EXAMPLES::

            sage: V = QQ^3
            sage: W = V.span([[1,2,3], [4,5,6]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
        """
        free_module.FreeModule_submodule_field.__init__(
            self, ambient=ambient, gens=gens, check=check, already_echelonized=already_echelonized)
        self._inner_product_matrix = inner_product_matrix

    def _repr_(self):
        """
        The default printing representation of ``self``.

        EXAMPLES::

            sage: V = VectorSpace(QQ, 5)
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

        The system representation can be overwritten, but leaves
        :meth:`_repr_` unmodified::

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

        Sparse vector spaces print this fact::

            sage: V = VectorSpace(QQ, 5, sparse=True)
            sage: U = V.submodule([ V.gen(i) - V.gen(0) for i in range(1,5) ])
            sage: U # indirect doctest
            Sparse vector space of degree 5 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
        """
        if self.is_sparse():
            return "Sparse quadratic space of degree %s and dimension %s over %s\n" % (
                self.degree(), self.dimension(), self.base_field()) + \
                "Basis matrix:\n%r\n" % self.basis_matrix() + \
                "Inner product matrix:\n%r" % self.inner_product_matrix()
        return "Quadratic space of degree %s and dimension %s over %s\n" % (
            self.degree(), self.dimension(), self.base_field()) + \
            "Basis matrix:\n%r\n" % self.basis_matrix() + \
            "Inner product matrix:\n%r" % self.inner_product_matrix()
