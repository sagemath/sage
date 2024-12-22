r"""
Projective `n` space over a ring

EXAMPLES:

We construct projective space over various rings of various dimensions.

The simplest projective space::

    sage: ProjectiveSpace(0)
    Projective Space of dimension 0 over Integer Ring

A slightly bigger projective space over `\QQ`::

    sage: X = ProjectiveSpace(1000, QQ); X
    Projective Space of dimension 1000 over Rational Field
    sage: X.dimension()
    1000

We can use "over" notation to create projective spaces over various
base rings.

::

    sage: X = ProjectiveSpace(5)/QQ; X
    Projective Space of dimension 5 over Rational Field
    sage: X/CC                                                                          # needs sage.rings.real_mpfr
    Projective Space of dimension 5 over Complex Field with 53 bits of precision

The third argument specifies the printing names of the generators of the
homogeneous coordinate ring. Using the method :meth:`objgens` you can obtain both
the space and the generators as ready to use variables. ::

    sage: P2, vars = ProjectiveSpace(10, QQ, 't').objgens()
    sage: vars
    (t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)

You can alternatively use the special syntax with ``<`` and ``>``.

::

    sage: P2.<x,y,z> = ProjectiveSpace(2, QQ)
    sage: P2
    Projective Space of dimension 2 over Rational Field
    sage: P2.coordinate_ring()
    Multivariate Polynomial Ring in x, y, z over Rational Field

The first of the three lines above is just equivalent to the two lines::

    sage: P2 = ProjectiveSpace(2, QQ, 'xyz')
    sage: x,y,z = P2.gens()

For example, we use `x,y,z` to define the intersection of
two lines.

::

    sage: V = P2.subscheme([x + y + z, x + y - z]); V
    Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
      x + y + z,
      x + y - z
    sage: V.dimension()                                                                 # needs sage.libs.singular
    0

AUTHORS:

- Ben Hutz: (June 2012): support for rings

- Ben Hutz (9/2014): added support for Cartesian products

- Rebecca Lauren Miller (March 2016) : added point_transformation_matrix
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import product

import sage.rings.abc

from sage.arith.misc import gcd, binomial
from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.categories.map import Map
from sage.categories.number_fields import NumberFields
from sage.categories.rings import Rings
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.permutation import Permutation
from sage.combinat.subset import Subsets
from sage.misc.latex import latex
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.misc.persist import register_unpickle_override
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.fraction_field import FractionField
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ, RationalField
from sage.schemes.generic.ambient_space import AmbientSpace
from sage.structure.category_object import normalize_names
from sage.structure.unique_representation import UniqueRepresentation
from sage.schemes.projective.projective_homset import (SchemeHomset_points_projective_ring,
                                                       SchemeHomset_points_projective_field,
                                                       SchemeHomset_polynomial_projective_space)
from sage.schemes.projective.projective_morphism import (SchemeMorphism_polynomial_projective_space,
                                                         SchemeMorphism_polynomial_projective_space_field,
                                                         SchemeMorphism_polynomial_projective_space_finite_field)
from sage.schemes.projective.projective_point import (SchemeMorphism_point_projective_ring,
                                                      SchemeMorphism_point_projective_field,
                                                      SchemeMorphism_point_projective_finite_field)

lazy_import('sage.combinat.integer_vector_weighted', 'WeightedIntegerVectors')
lazy_import('sage.combinat.tuple', ['Tuples', 'UnorderedTuples'])
lazy_import('sage.dynamics.arithmetic_dynamics.projective_ds', 'DynamicalSystem_projective')
lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.modules.free_module_element', 'prepare')
lazy_import('sage.schemes.generic.algebraic_scheme', 'AlgebraicScheme_subscheme')
lazy_import('sage.schemes.product_projective.space',
            ['ProductProjectiveSpaces', 'ProductProjectiveSpaces_ring'])
lazy_import('sage.schemes.projective.projective_subscheme',
            ['AlgebraicScheme_subscheme_projective', 'AlgebraicScheme_subscheme_projective_field'])


# for better efficiency
_Fields = Fields()
_Rings = Rings()
_CommRings = _Rings.Commutative()


def is_ProjectiveSpace(x):
    r"""
    Return ``True`` if ``x`` is a projective space.

    In other words, if ``x`` is an ambient space `\mathbb{P}^n_R`,
    where `R` is a ring and `n\geq 0` is an integer.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_space import is_ProjectiveSpace
        sage: is_ProjectiveSpace(ProjectiveSpace(5, names='x'))
        doctest:warning...
        DeprecationWarning: The function is_ProjectiveSpace is deprecated; use 'isinstance(..., ProjectiveSpace_ring)' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        True
        sage: is_ProjectiveSpace(ProjectiveSpace(5, GF(9, 'alpha'), names='x'))         # needs sage.rings.finite_rings
        True
        sage: is_ProjectiveSpace(Spec(ZZ))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_ProjectiveSpace is deprecated; use 'isinstance(..., ProjectiveSpace_ring)' instead.")
    return isinstance(x, ProjectiveSpace_ring)


def ProjectiveSpace(n, R=None, names=None):
    r"""
    Return projective space of dimension ``n`` over the ring ``R``.

    EXAMPLES: The dimension and ring can be given in either order.

    ::

        sage: ProjectiveSpace(3, QQ)
        Projective Space of dimension 3 over Rational Field
        sage: ProjectiveSpace(5, QQ)
        Projective Space of dimension 5 over Rational Field
        sage: P = ProjectiveSpace(2, QQ, names='XYZ'); P
        Projective Space of dimension 2 over Rational Field
        sage: P.coordinate_ring()
        Multivariate Polynomial Ring in X, Y, Z over Rational Field

    The divide operator does base extension.

    ::

        sage: ProjectiveSpace(5)/GF(17)
        Projective Space of dimension 5 over Finite Field of size 17

    The default base ring is `\ZZ`.

    ::

        sage: ProjectiveSpace(5)
        Projective Space of dimension 5 over Integer Ring

    There is also a projective space associated each polynomial ring.

    ::

        sage: R = GF(7)['x,y,z']
        sage: P = ProjectiveSpace(R); P
        Projective Space of dimension 2 over Finite Field of size 7
        sage: P.coordinate_ring()
        Multivariate Polynomial Ring in x, y, z over Finite Field of size 7
        sage: P.coordinate_ring() is R
        True

    ::

        sage: ProjectiveSpace(3, Zp(5), 'y')                                            # needs sage.rings.padics
        Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20

    ::

        sage: ProjectiveSpace(2, QQ, 'x,y,z')
        Projective Space of dimension 2 over Rational Field

    ::

        sage: PS.<x,y> = ProjectiveSpace(1, CC); PS                                     # needs sage.rings.real_mpfr
        Projective Space of dimension 1 over Complex Field with 53 bits of precision

    ::

        sage: R.<x,y,z> = QQ[]
        sage: ProjectiveSpace(R).variable_names()
        ('x', 'y', 'z')

    Projective spaces are not cached, i.e., there can be several with
    the same base ring and dimension (to facilitate gluing
    constructions).

    ::

        sage: R.<x> = QQ[]
        sage: ProjectiveSpace(R)
        Projective Space of dimension 0 over Rational Field

    TESTS::

        sage: R.<x,y> = QQ[]
        sage: P.<z,w> = ProjectiveSpace(R)
        Traceback (most recent call last):
        ...
        NameError: variable names passed to ProjectiveSpace conflict with names in ring

    ::

        sage: R.<x,y> = QQ[]
        sage: P.<x,y> = ProjectiveSpace(R)
        sage: P.gens() == R.gens()
        True
    """
    if isinstance(n, (MPolynomialRing_base, PolynomialRing_general)) and R is None:
        if names is not None:
            # Check for the case that the user provided a variable name
            # That does not match what we wanted to use from R
            names = normalize_names(n.ngens(), names)
            if n.variable_names() != names:
                # The provided name doesn't match the name of R's variables
                raise NameError("variable names passed to ProjectiveSpace conflict with names in ring")
        A = ProjectiveSpace(n.ngens() - 1, n.base_ring(),
                            names=n.variable_names())
        A._coordinate_ring = n
        return A
    if names is None:
        names = 'x'
    if isinstance(R, (Integer, int)):
        n, R = R, n
    if R is None:
        R = ZZ  # default is the integers
    if R in _Fields:
        if isinstance(R, FiniteField):
            return ProjectiveSpace_finite_field(n, R, names)
        if isinstance(R, RationalField):
            return ProjectiveSpace_rational_field(n, R, names)
        else:
            return ProjectiveSpace_field(n, R, names)
    elif R in _CommRings:
        return ProjectiveSpace_ring(n, R, names)
    raise TypeError("R (=%s) must be a commutative ring" % R)


class ProjectiveSpace_ring(UniqueRepresentation, AmbientSpace):
    """
    Projective space of dimension `n` over the ring
    `R`.

    EXAMPLES::

        sage: X.<x,y,z,w> = ProjectiveSpace(3, QQ)
        sage: X.base_scheme()
        Spectrum of Rational Field
        sage: X.base_ring()
        Rational Field
        sage: X.structure_morphism()
        Scheme morphism:
          From: Projective Space of dimension 3 over Rational Field
          To:   Spectrum of Rational Field
          Defn: Structure map
        sage: X.coordinate_ring()
        Multivariate Polynomial Ring in x, y, z, w over Rational Field

    Loading and saving::

        sage: loads(X.dumps()) == X
        True
        sage: P = ProjectiveSpace(ZZ, 1, 'x')
        sage: loads(P.dumps()) is P
        True

    Equality and hashing::

        sage: ProjectiveSpace(QQ, 3, 'a') == ProjectiveSpace(ZZ, 3, 'a')
        False
        sage: ProjectiveSpace(ZZ, 1, 'a') == ProjectiveSpace(ZZ, 0, 'a')
        False
        sage: ProjectiveSpace(ZZ, 2, 'a') == AffineSpace(ZZ, 2, 'a')
        False

        sage: ProjectiveSpace(QQ, 3, 'a') != ProjectiveSpace(ZZ, 3, 'a')
        True
        sage: ProjectiveSpace(ZZ, 1, 'a') != ProjectiveSpace(ZZ, 0, 'a')
        True
        sage: ProjectiveSpace(ZZ, 2, 'a') != AffineSpace(ZZ, 2, 'a')
        True

        sage: hash(ProjectiveSpace(QQ, 3, 'a')) == hash(ProjectiveSpace(ZZ, 3, 'a'))
        False
        sage: hash(ProjectiveSpace(ZZ, 1, 'a')) == hash(ProjectiveSpace(ZZ, 0, 'a'))
        False
        sage: hash(ProjectiveSpace(ZZ, 2, 'a')) == hash(AffineSpace(ZZ, 2, 'a'))
        False
    """
    @staticmethod
    def __classcall__(cls, n, R=ZZ, names=None):
        """
        EXAMPLES::

            sage: ProjectiveSpace(QQ, 2, names='XYZ') is ProjectiveSpace(QQ, 2, names='XYZ')
            True
        """
        normalized_names = normalize_names(n + 1, names)
        return super().__classcall__(cls, n, R, normalized_names)

    def __init__(self, n, R=ZZ, names=None):
        """
        Initialization function.

        EXAMPLES::

            sage: ProjectiveSpace(3, Zp(5), 'y')                                        # needs sage.rings.padics
            Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20
        """
        AmbientSpace.__init__(self, n, R)
        self._assign_names(names)

    def ngens(self):
        """
        Return the number of generators of this projective space.

        This is the number of variables in the coordinate ring of ``self``.

        EXAMPLES::

            sage: ProjectiveSpace(3, QQ).ngens()
            4
            sage: ProjectiveSpace(7, ZZ).ngens()
            8
        """
        return self.dimension_relative() + 1

    def _check_satisfies_equations(self, v):
        """
        Return ``True`` if ``v`` defines a point on the scheme; raise a
        :exc:`TypeError` otherwise.

        EXAMPLES::

            sage: P = ProjectiveSpace(2, ZZ)
            sage: P._check_satisfies_equations([1, 1, 0])
            True

        ::

            sage: P = ProjectiveSpace(1, QQ)
            sage: P._check_satisfies_equations((1/2, 0))
            True

        ::

            sage: P = ProjectiveSpace(2, ZZ)
            sage: P._check_satisfies_equations([0, 0, 0])
            Traceback (most recent call last):
            ...
            TypeError: the zero vector is not a point in projective space

        ::

            sage: P = ProjectiveSpace(2, ZZ)
            sage: P._check_satisfies_equations((1, 0))
            Traceback (most recent call last):
            ...
            TypeError: the list v=(1, 0) must have 3 components

        ::

            sage: P = ProjectiveSpace(2, ZZ)
            sage: P._check_satisfies_equations([1/2, 0, 1])
            Traceback (most recent call last):
            ...
            TypeError: the components of v=[1/2, 0, 1] must be elements of Integer Ring
        """
        if not isinstance(v, (list, tuple)):
            raise TypeError('the argument v=%s must be a list or tuple' % v)
        n = self.ngens()
        if not len(v) == n:
            raise TypeError('the list v=%s must have %s components' % (v, n))
        R = self.base_ring()
        for coord in v:
            if coord not in R:
                raise TypeError('the components of v=%s must be elements of %s' % (v, R))
        zero = [R(0)] * n
        if v == zero:
            raise TypeError('the zero vector is not a point in projective space')
        return True

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme.

        EXAMPLES::

            sage: ProjectiveSpace(3, GF(19^2,'alpha'), 'abcd').coordinate_ring()        # needs sage.rings.finite_rings
            Multivariate Polynomial Ring in a, b, c, d over Finite Field in alpha of size 19^2

        ::

            sage: ProjectiveSpace(3).coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring

        ::

            sage: ProjectiveSpace(2, QQ, ['alpha', 'beta', 'gamma']).coordinate_ring()
            Multivariate Polynomial Ring in alpha, beta, gamma over Rational Field
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            self._coordinate_ring = PolynomialRing(self.base_ring(),
                                        self.variable_names(),
                                        self.dimension_relative() + 1)
            return self._coordinate_ring

    def _validate(self, polynomials):
        """
        If ``polynomials`` is a tuple of valid polynomial functions on
        ``self``, return ``polynomials``, otherwise raise :exc:`TypeError`.

        Since this is a projective space, polynomials must be homogeneous.

        INPUT:

        - ``polynomials`` -- tuple of polynomials in the coordinate ring of
          this space

        OUTPUT: tuple of polynomials in the coordinate ring of this space

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._validate([x*y - z^2, x])
            [x*y - z^2, x]

        ::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._validate((x*y - z, x))
            Traceback (most recent call last):
            ...
            TypeError: x*y - z is not a homogeneous polynomial

        ::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._validate(x*y - z)
            Traceback (most recent call last):
            ...
            TypeError: the argument polynomials=x*y - z must be a list or tuple
        """
        if not isinstance(polynomials, (list, tuple)):
            raise TypeError('the argument polynomials=%s must be a list or tuple' % polynomials)
        for f in polynomials:
            if not f.is_homogeneous():
                raise TypeError("%s is not a homogeneous polynomial" % f)
        return polynomials

    def __pow__(self, m):
        """
        Return the Cartesian power of this space.

        INPUT:

        - ``m`` -- integer

        OUTPUT: product of projective spaces

        EXAMPLES::

            sage: P = ProjectiveSpace(1, QQ, 'x')
            sage: P3 = P^3; P3
            Product of projective spaces P^1 x P^1 x P^1 over Rational Field
            sage: P3.variable_names()
            ('x0', 'x1', 'x2', 'x3', 'x4', 'x5')

        As you see, custom variable names are not preserved by power operator,
        since there is no natural way to make new ones in general.
        """
        mm = int(m)
        if mm != m:
            raise ValueError("m must be an integer")
        return ProductProjectiveSpaces([self.dimension_relative()] * mm, self.base_ring())

    def __mul__(self, right):
        r"""
        Create the product of projective spaces.

        INPUT:

        - ``right`` -- a projective space, product of projective spaces, or subscheme

        OUTPUT: a product of projective spaces or subscheme

        EXAMPLES::

            sage: P1 = ProjectiveSpace(QQ, 1, 'x')
            sage: P2 = ProjectiveSpace(QQ, 2, 'y')
            sage: P1*P2
            Product of projective spaces P^1 x P^2 over Rational Field

            ::

            sage: S.<t,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: T.<a,b> = ProjectiveSpace(QQ, 1)
            sage: T*S
            Product of projective spaces P^1 x P^3 x P^2 over Rational Field

        ::

            sage: S = ProjectiveSpace(ZZ, 2, 't')
            sage: T = ProjectiveSpace(ZZ, 3, 'x')
            sage: T.inject_variables()
            Defining x0, x1, x2, x3
            sage: X = T.subscheme([x0*x2 - x1*x3])
            sage: S*X
            Closed subscheme of Product of projective spaces P^2 x P^3 over Integer Ring defined by:
              x0*x2 - x1*x3

        ::

            sage: S = ProjectiveSpace(QQ, 3, 'x')
            sage: T = AffineSpace(2, QQ, 'y')
            sage: S*T
            Traceback (most recent call last):
            ...
            TypeError: Affine Space of dimension 2 over Rational Field must be a
            projective space, product of projective spaces, or subscheme
        """
        if self.base_ring() != right.base_ring():
            raise ValueError('Must have the same base ring')

        if isinstance(right, ProductProjectiveSpaces_ring):
            return ProductProjectiveSpaces([self] + right.components())
        elif isinstance(right, ProjectiveSpace_ring):
            if self is right:
                return self.__pow__(2)
            return ProductProjectiveSpaces([self, right])
        elif isinstance(right, AlgebraicScheme_subscheme):
            AS = self * right.ambient_space()
            CR = AS.coordinate_ring()
            n = self.ambient_space().coordinate_ring().ngens()

            phi = self.ambient_space().coordinate_ring().hom(list(CR.gens()[:n]), CR)
            psi = right.ambient_space().coordinate_ring().hom(list(CR.gens()[n:]), CR)
            return AS.subscheme([phi(t) for t in self.defining_polynomials()] + [psi(t) for t in right.defining_polynomials()])
        else:
            raise TypeError('%s must be a projective space, product of projective spaces, or subscheme' % right)

    def _latex_(self):
        r"""
        Return a LaTeX representation of this projective space.

        EXAMPLES::

            sage: print(latex(ProjectiveSpace(1, ZZ, 'x')))
            {\mathbf P}_{\Bold{Z}}^{1}

        TESTS::

            sage: ProjectiveSpace(11, Zp(5), 'y')._latex_()                             # needs sage.rings.padics
            '{\\mathbf P}_{\\Bold{Z}_{5}}^{11}'
        """
        return "{\\mathbf P}_{%s}^{%s}" % (latex(self.base_ring()), self.dimension_relative())

    def _linear_system_as_kernel(self, d, pt, m):
        """
        Return a matrix whose kernel consists of the coefficient vectors
        of the degree ``d`` hypersurfaces (wrt lexicographic ordering of its
        monomials) with multiplicity at least ``m`` at ``pt``.

        INPUT:

        - ``d`` -- nonnegative integer

        - ``pt`` -- a point of ``self`` (possibly represented by a list with at
          least one component equal to 1)

        - ``m`` -- nonnegative integer

        OUTPUT:

        A matrix of size `\binom{m-1+n}{n}` x `\binom{d+n}{n}` where n is the
        relative dimension of ``self``. The base ring of the matrix is a ring that
        contains the base ring of ``self`` and the coefficients of the given point.

        EXAMPLES:

        If the degree `d` is 0, then a matrix consisting of the first unit vector
        is returned::

            sage: P = ProjectiveSpace(GF(5), 2, names='x')
            sage: pt = P([1, 1, 1])
            sage: P._linear_system_as_kernel(0, pt, 3)                                  # needs sage.modules
            [1]
            [0]
            [0]
            [0]
            [0]
            [0]

        If the multiplicity `m` is 0, then a matrix with zero rows
        is returned::

            sage: P = ProjectiveSpace(GF(5), 2, names='x')
            sage: pt = P([1, 1, 1])
            sage: M = P._linear_system_as_kernel(2, pt, 0)                              # needs sage.modules
            sage: [M.nrows(), M.ncols()]                                                # needs sage.modules
            [0, 6]

        The base ring does not need to be a field or even an integral domain.
        In this case, the point can be given by a list::

            sage: R = Zmod(4)
            sage: P = ProjectiveSpace(R, 2, names='x')
            sage: pt = [R(1), R(3), R(0)]
            sage: P._linear_system_as_kernel(3, pt, 2)                                  # needs sage.modules
            [1 3 0 1 0 0 3 0 0 0]
            [0 1 0 2 0 0 3 0 0 0]
            [0 0 1 0 3 0 0 1 0 0]

        When representing a point by a list at least one component must be 1
        (even when the base ring is a field and the list gives a well-defined
        point in projective space)::

            sage: R = GF(5)
            sage: P = ProjectiveSpace(R, 2, names='x')
            sage: pt = [R(3), R(3), R(0)]
            sage: P._linear_system_as_kernel(3, pt, 2)                                  # needs sage.modules
            Traceback (most recent call last):
            ...
            TypeError: at least one component of pt=[3, 3, 0] must be equal to 1

        The components of the list do not have to be elements of the base ring
        of the projective space. It suffices if there exists a common parent.
        For example, the kernel of the following matrix corresponds to
        hypersurfaces of degree 2 in 3-space with multiplicity at least 2 at a
        general point in the third affine patch::

            sage: P = ProjectiveSpace(QQ, 3, names='x')
            sage: RPol.<t0,t1,t2,t3> = PolynomialRing(QQ, 4)
            sage: pt = [t0,t1,1,t3]
            sage: P._linear_system_as_kernel(2, pt, 2)                                  # needs sage.modules
            [ 2*t0    t1     1    t3     0     0     0     0     0     0]
            [    0    t0     0     0  2*t1     1    t3     0     0     0]
            [ t0^2 t0*t1    t0 t0*t3  t1^2    t1 t1*t3     1    t3  t3^2]
            [    0     0     0    t0     0     0    t1     0     1  2*t3]

        .. TODO::

            Use this method as starting point to implement a class
            LinearSystem for linear systems of hypersurfaces.
        """
        if not isinstance(d, (int, Integer)):
            raise TypeError('the argument d=%s must be an integer' % d)
        if d < 0:
            raise ValueError('the integer d=%s must be nonnegative' % d)
        if not isinstance(pt, (list, tuple,
                               SchemeMorphism_point_projective_ring)):
            raise TypeError('the argument pt=%s must be a list, tuple, or '
                            'point on a projective space' % pt)
        pt, R = prepare(pt, None)
        n = self.dimension_relative()
        if not len(pt) == n + 1:
            raise TypeError('the sequence pt=%s must have %s '
                            'components' % (pt, n + 1))
        if not R.has_coerce_map_from(self.base_ring()):
            raise TypeError('unable to find a common ring for all elements')
        try:
            i = pt.index(1)
        except Exception:
            raise TypeError('at least one component of pt=%s must be equal '
                            'to 1' % pt)
        pt = pt[:i] + pt[i + 1:]
        if not isinstance(m, (int, Integer)):
            raise TypeError('the argument m=%s must be an integer' % m)
        if m < 0:
            raise ValueError('the integer m=%s must be nonnegative' % m)
        # the components of partials correspond to partial derivatives
        # of order at most m-1 with respect to n variables
        partials = IntegerVectors(m - 1, n + 1).list()
        # the components of monoms correspond to monomials of degree
        # at most d in n variables
        monoms = IntegerVectors(d, n + 1).list()
        M = matrix(R, len(partials), len(monoms))
        for row in range(M.nrows()):
            e = partials[row][:i] + partials[row][i + 1:]
            for col in range(M.ncols()):
                f = monoms[col][:i] + monoms[col][i + 1:]
                if all(f[j] >= e[j] for j in range(n)):
                    M[row, col] = prod(binomial(fj, ej) * ptj**(fj - ej)
                                       for ptj, fj, ej in zip(pt, f, e)
                                       if fj > ej)
        return M

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: P2._morphism(P2.Hom(P2), [x,y,z])
            Scheme endomorphism of Projective Space of dimension 2 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x : y : z)
        """
        return SchemeMorphism_polynomial_projective_space(*args, **kwds)

    def _homset(self, *args, **kwds):
        """                                          ii
        Construct the Hom-set

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(2, QQ)
            sage: Hom(P, P)
            Set of morphisms
              From: Projective Space of dimension 2 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
        """
        return SchemeHomset_polynomial_projective_space(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: P2._point_homset(Spec(GF(3)), P2)
            Set of rational points of Projective Space of dimension 2 over Finite Field of size 3
        """
        return SchemeHomset_points_projective_ring(*args, **kwds)

    def point(self, v, check=True):
        """
        Create a point on this projective space.

        INPUT:

        - ``v`` -- anything that defines a point

        - ``check`` -- boolean (default: ``True``); whether
          to check the defining data for consistency

        OUTPUT: a point of this projective space

        EXAMPLES::

            sage: P2 = ProjectiveSpace(QQ, 2)
            sage: P2.point([4,5])
            (4 : 5 : 1)

        ::

            sage: P = ProjectiveSpace(QQ, 1)
            sage: P.point(infinity)
            (1 : 0)

        ::

            sage: P = ProjectiveSpace(QQ, 2)
            sage: P.point(infinity)
            Traceback (most recent call last):
            ...
            ValueError: +Infinity not well defined in dimension > 1

        ::

            sage: P = ProjectiveSpace(ZZ, 2)
            sage: P.point([infinity])
            Traceback (most recent call last):
             ...
            ValueError: [+Infinity] not well defined in dimension > 1
        """
        from sage.rings.infinity import infinity
        if v is infinity or (isinstance(v, (list, tuple)) and
                             len(v) == 1 and v[0] is infinity):
            if self.dimension_relative() > 1:
                raise ValueError("%s not well defined in dimension > 1" % v)
            v = [1, 0]

        return self.point_homset()(v, check=check)

    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: point_homset = P2._point_homset(Spec(GF(3)), P2)
            sage: P2._point(point_homset, [1,2,3])
            (2 : 1 : 0)
        """
        return SchemeMorphism_point_projective_ring(*args, **kwds)

    def _repr_(self):
        """
        Return a string representation of this projective space.

        EXAMPLES::

            sage: ProjectiveSpace(1, ZZ, 'x')
            Projective Space of dimension 1 over Integer Ring

        TESTS::

            sage: ProjectiveSpace(3, Zp(5), 'y')._repr_()                               # needs sage.rings.padics
            'Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20'
        """
        return "Projective Space of dimension %s over %s" % (self.dimension_relative(), self.base_ring())

    def _repr_generic_point(self, v=None):
        """
        Return a string representation of the generic point
        corresponding to the list of polys ``v`` on this projective space.

        If ``v`` is None, the representation of the generic point of
        the projective space is returned.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._repr_generic_point([z*y - x^2])
            '(-x^2 + y*z)'
            sage: P._repr_generic_point()
            '(x : y : z)'
        """
        if v is None:
            v = self.gens()
        return '(%s)' % (" : ".join(repr(f) for f in v))

    def _latex_generic_point(self, v=None):
        """
        Return a LaTeX representation of the generic point
        corresponding to the list of polys ``v`` on this projective space.

        If ``v`` is None, the representation of the generic point of
        the projective space is returned.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._latex_generic_point([z*y - x^2])
            '\\left(-x^{2} + y z\\right)'
            sage: P._latex_generic_point()
            '\\left(x : y : z\\right)'
        """
        if v is None:
            v = self.gens()
        return '\\left(%s\\right)' % (" : ".join(str(latex(f)) for f in v))

    def change_ring(self, R):
        r"""
        Return a projective space over ring ``R``.

        INPUT:

        - ``R`` -- commutative ring or morphism

        OUTPUT: projective space over ``R``

        .. NOTE::

            There is no need to have any relation between ``R`` and the base ring
            of this space, if you want to have such a relation, use
            ``self.base_extend(R)`` instead.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: PQ = P.change_ring(QQ); PQ
            Projective Space of dimension 2 over Rational Field
            sage: PQ.change_ring(GF(5))
            Projective Space of dimension 2 over Finite Field of size 5

        ::

            sage: K.<w> = QuadraticField(2)                                             # needs sage.rings.number_field
            sage: P = ProjectiveSpace(K, 2, 't')                                        # needs sage.rings.number_field
            sage: P.change_ring(K.embeddings(QQbar)[0])                                 # needs sage.rings.number_field
            Projective Space of dimension 2 over Algebraic Field
        """
        if isinstance(R, Map):
            return ProjectiveSpace(self.dimension_relative(), R.codomain(),
                               self.variable_names())
        else:
            return ProjectiveSpace(self.dimension_relative(), R,
                               self.variable_names())

    def is_projective(self):
        """
        Return that this ambient space is projective `n`-space.

        EXAMPLES::

            sage: ProjectiveSpace(3,QQ).is_projective()
            True
        """
        return True

    def subscheme(self, X):
        """
        Return the closed subscheme defined by ``X``.

        INPUT:

        - ``X`` -- list or tuple of equations

        EXAMPLES::

            sage: A.<x,y,z> = ProjectiveSpace(2, QQ)
            sage: X = A.subscheme([x*z^2, y^2*z, x*y^2]); X
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*z^2,
              y^2*z,
              x*y^2
            sage: X.defining_polynomials ()
            (x*z^2, y^2*z, x*y^2)
            sage: I = X.defining_ideal(); I
            Ideal (x*z^2, y^2*z, x*y^2) of Multivariate Polynomial Ring in x, y, z
             over Rational Field
            sage: I.groebner_basis()                                                    # needs sage.libs.singular
            [x*y^2, y^2*z,  x*z^2]
            sage: X.dimension()                                                         # needs sage.libs.singular
            0
            sage: X.base_ring()
            Rational Field
            sage: X.base_scheme()
            Spectrum of Rational Field
            sage: X.structure_morphism()
            Scheme morphism:
              From: Closed subscheme of Projective Space of dimension 2
                    over Rational Field defined by: x*z^2, y^2*z, x*y^2
              To:   Spectrum of Rational Field
              Defn: Structure map

        TESTS::

            sage: TestSuite(X).run(skip=["_test_an_element", "_test_elements",\
            ....: "_test_elements_eq", "_test_some_elements", "_test_elements_eq_reflexive",\
            ....: "_test_elements_eq_symmetric", "_test_elements_eq_transitive",\
            ....: "_test_elements_neq"])
        """
        R = self.base_ring()
        if R.is_field() and R.is_exact():
            return AlgebraicScheme_subscheme_projective_field(self, X)

        return AlgebraicScheme_subscheme_projective(self, X)

    def points_of_bounded_height(self, **kwds):
        r"""
        Return an iterator of the points in ``self`` of absolute multiplicative
        height of at most the given bound.

        ALGORITHM:

        This is an implementation of Algorithm 6 in [Krumm2016]_.

        INPUT: keyword arguments:

        - ``bound`` -- a real number

        - ``precision`` -- (default: 53) a positive integer

        OUTPUT: an iterator of points of bounded height

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: sorted(list(P.points_of_bounded_height(bound=2)))
            [(-2 : 1), (-1 : 1), (-1/2 : 1), (0 : 1),
             (1/2 : 1), (1 : 0), (1 : 1), (2 : 1)]

        ::

            sage: u = QQ['u'].0
            sage: P.<x,y,z> = ProjectiveSpace(NumberField(u^2 - 2, 'v'), 2)             # needs sage.rings.number_field
            sage: len(list(P.points_of_bounded_height(bound=2)))                        # needs sage.rings.number_field
            265

        ::

            sage: # needs sage.rings.number_field
            sage: CF.<a> = CyclotomicField(3)
            sage: R.<x> = CF[]
            sage: L.<l> = CF.extension(x^3 + 2)
            sage: Q.<x,y> = ProjectiveSpace(L, 1)
            sage: sorted(list(Q.points_of_bounded_height(bound=1)))
            [(0 : 1), (1 : 0), (a + 1 : 1), (a : 1),
             (-1 : 1), (-a - 1 : 1), (-a : 1), (1 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: F.<a> = NumberField(x^4 - 8*x^2 + 3)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2)
            sage: all(exp(p.global_height()) <= 1                                       # needs sage.symbolic
            ....:     for p in P.points_of_bounded_height(bound=1))
            True

        ::

            sage: K.<a> = CyclotomicField(3)                                            # needs sage.rings.number_field
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)                                     # needs sage.rings.number_field
            sage: len(list(P.points_of_bounded_height(bound=1)))                        # needs sage.rings.number_field
            57

        ::

            sage: u = QQ['u'].0
            sage: K.<k> = NumberField(u^2 - 2)                                          # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(K, 1)                                       # needs sage.rings.number_field
            sage: len(list(P.points_of_bounded_height(bound=2)))                        # needs sage.rings.number_field
            24

        ::

            sage: R.<x> = QQ[]
            sage: K.<k> = NumberField(x^4 - 8*x^2 + 3)                                  # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(K, 1)                                       # needs sage.rings.number_field
            sage: len(list(P.points_of_bounded_height(bound=2)))                        # needs sage.rings.number_field
            108

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<v> = NumberField(x^5 + x^3 + 1)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: L = P.points_of_bounded_height(bound=1.2)
            sage: len(list(L))
            109

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = QuadraticField(2)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: sorted(list(P.points_of_bounded_height(bound=2)))
            [(-v - 2 : 1), (-v - 1 : 1), (-2 : 1), (-1/2*v - 1 : 1), (-v : 1), (-1 : 1),
             (-1/2*v : 1), (v - 2 : 1), (-1/2 : 1), (-v + 1 : 1), (1/2*v - 1 : 1), (0 : 1),
             (-1/2*v + 1 : 1), (v - 1 : 1), (1/2 : 1), (-v + 2 : 1), (1/2*v : 1), (1 : 0),
             (1 : 1), (v : 1), (1/2*v + 1 : 1), (2 : 1), (v + 1 : 1), (v + 2 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(3*x^2 + 1)
            sage: P.<z,w> = ProjectiveSpace(K, 1)
            sage: sorted(list(P.points_of_bounded_height(bound=1)))
            [(-1 : 1), (-3/2*a - 1/2 : 1), (3/2*a - 1/2 : 1), (0 : 1),
             (-3/2*a + 1/2 : 1), (3/2*a + 1/2 : 1), (1 : 0), (1 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(3*x^2 + 1)
            sage: O = K.maximal_order()
            sage: P.<z,w> = ProjectiveSpace(O, 1)
            sage: len(sorted(list(P.points_of_bounded_height(bound=2))))
            44

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3 - 7)
            sage: O = K.maximal_order()
            sage: P.<z,w> = ProjectiveSpace(O, 1)
            sage: len(sorted(list(P.points_of_bounded_height(bound=2))))
            28

        ::

            sage: P.<w,z> = ProjectiveSpace(ZZ, 1)
            sage: sorted(list(P.points_of_bounded_height(bound=2)))
            [(-2 : -1), (-2 : 1), (-1 : -2), (-1 : -1),
             (-1 : 0), (-1 : 1), (-1 : 2), (0 : -1)]

        ::

            sage: R.<x> = QQ[]
            sage: P.<z,w> = ProjectiveSpace(R, 1)
            sage: P.points_of_bounded_height(bound=2)
            Traceback (most recent call last):
            ...
            NotImplementedError: self must be a projective space over
            a number field or a ring of integers

        ::

            sage: # needs sage.rings.number_field
            sage: K.<i> = NumberField(x^2 + 1)
            sage: PK.<t> = K[]
            sage: L.<a> = K.extension(t^4  - i)
            sage: P.<z,w> = ProjectiveSpace(L, 1)
            sage: sorted(list(P.points_of_bounded_height(bound=1)))
            [(0 : 1), (1 : 0), (a : 1), (a^2 : 1), (a^3 : 1), (i : 1),
             (i*a : 1), (i*a^2 : 1), (i*a^3 : 1), (-1 : 1), (-a : 1), (-a^2 : 1),
             (-a^3 : 1), (-i : 1), (-i*a : 1), (-i*a^2 : 1), (-i*a^3 : 1), (1 : 1)]
        """
        from sage.schemes.projective.proj_bdd_height import (
            ZZ_points_of_bounded_height,
            QQ_points_of_bounded_height,
            IQ_points_of_bounded_height,
            points_of_bounded_height
        )

        R = self.base_ring()

        # Check the base ring is the rational field, a number field,
        # or the ring of integers
        is_ring_of_ints = False

        if isinstance(R, RationalField):
            field_type = False
        elif R in NumberFields():
            # True for the rational field as well, so check RationalField first
            field_type = True
        elif R is ZZ or (isinstance(R, sage.rings.abc.Order) and R.is_integrally_closed()):  # Ensure ring of integers / maximal order
            is_ring_of_ints = True
        else:
            raise NotImplementedError("self must be a projective space over a number field or a ring of integers")

        bound = kwds.pop('bound')
        prec = kwds.pop('precision', 53)

        # Convert between absolute and relative height for calling Krumm's algorithm
        bound = bound**R.absolute_degree()

        dim = self.dimension_relative()

        # When R is the ring of integers
        if is_ring_of_ints:
            fraction_field = FractionField(R)

            # Field of fraction is the rational field
            if fraction_field == QQ:
                return ZZ_points_of_bounded_height(self, dim, bound)

            # Field of fraction is a number field
            r1, r2 = fraction_field.signature()
            r = r1 + r2 - 1

            if fraction_field.is_relative():
                deg = fraction_field.relative_degree()
            else:
                deg = fraction_field.degree()

            if deg == 2 and r == 0:
                return IQ_points_of_bounded_height(self, fraction_field, dim, bound)

            return points_of_bounded_height(self, fraction_field, dim, bound, prec)

        # When R is a field
        if field_type:
            # For checking whether R is imaginary quadratic field
            r1, r2 = R.signature()
            r = r1 + r2 - 1

            if R.is_relative():
                deg = R.relative_degree()
            else:
                deg = R.degree()

            if deg == 2 and r == 0:
                return IQ_points_of_bounded_height(self, R, dim, bound)

            return points_of_bounded_height(self, R, dim, bound, prec)
        else:
            return QQ_points_of_bounded_height(self, dim, bound)

    def affine_patch(self, i, AA=None):
        r"""
        Return the `i`-th affine patch of this projective space.

        This is an ambient affine space `\mathbb{A}^n_R,` where
        `R` is the base ring of ``self``, whose "projective embedding"
        map is `1` in the `i`-th factor.

        INPUT:

        - ``i`` -- integer between 0 and dimension of ``self``, inclusive

        - ``AA`` -- (default: ``None``) ambient affine space, this is constructed
          if it is not given

        OUTPUT: an ambient affine space with fixed projective_embedding map

        EXAMPLES::

            sage: PP = ProjectiveSpace(5) / QQ
            sage: AA = PP.affine_patch(2)
            sage: AA
            Affine Space of dimension 5 over Rational Field
            sage: AA.projective_embedding()
            Scheme morphism:
              From: Affine Space of dimension 5 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1, x3, x4, x5) to
                    (x0 : x1 : 1 : x3 : x4 : x5)
            sage: AA.projective_embedding(0)
            Scheme morphism:
              From: Affine Space of dimension 5 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1, x3, x4, x5) to
                    (1 : x0 : x1 : x3 : x4 : x5)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.affine_patch(0).projective_embedding(0).codomain() == P
            True
        """
        i = int(i)   # implicit type checking
        n = self.dimension_relative()
        if i < 0 or i > n:
            raise ValueError("argument i (= %s) must be between 0 and %s" % (i, n))
        try:
            A = self.__affine_patches[i]
            # assume that if you've passed in a new affine space you
            # want to override the existing patch
            if AA is None or A == AA:
                return A
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        # if no ith patch exists, we may still be here with AA==None
        if AA is None:
            from sage.schemes.affine.affine_space import AffineSpace
            g = self.gens()
            gens = g[:i] + g[i + 1:]
            AA = AffineSpace(n, self.base_ring(), names=gens,
                             ambient_projective_space=self,
                             default_embedding_index=i)
        elif AA.dimension_relative() != n:
            raise ValueError("affine space must be of the dimension %s" % (n))
        self.__affine_patches[i] = AA
        return AA

    def _an_element_(self):
        r"""
        Return a (preferably typical) element of this space.

        This is used both for illustration and testing purposes.

        OUTPUT: a point in this projective space

        EXAMPLES::

            sage: ProjectiveSpace(ZZ, 3, 'x').an_element()
            (7 : 6 : 5 : 1)

            sage: ProjectiveSpace(PolynomialRing(ZZ,'y'), 3, 'x').an_element()
            (7*y : 6*y : 5*y : 1)
        """
        n = self.dimension_relative()
        R = self.base_ring()
        return self([(7 - i) * R.an_element() for i in range(n)] + [R.one()])

    def Lattes_map(self, E, m):
        r"""
        Given an elliptic curve ``E`` and an integer ``m`` return
        the Lattes map associated to multiplication by `m`.

        In other words, the rational map on the quotient
        `E/\{\pm 1\} \cong \mathbb{P}^1` associated to `[m]:E \to E`.

        INPUT:

        - ``E`` -- an elliptic curve

        - ``m`` -- integer

        OUTPUT: a dynamical system on this projective space

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: E = EllipticCurve(QQ,[-1, 0])                                         # needs sage.schemes
            sage: P.Lattes_map(E, 2)                                                    # needs sage.schemes
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (1/4*x^4 + 1/2*x^2*y^2 + 1/4*y^4 : x^3*y - x*y^3)

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(GF(37), 1)
            sage: E = EllipticCurve([1, 1])                                             # needs sage.rings.finite_rings sage.schemes
            sage: f = P.Lattes_map(E, 2); f                                             # needs sage.rings.finite_rings sage.schemes
            Dynamical System of Projective Space of dimension 1 over Finite Field of size 37
              Defn: Defined on coordinates by sending (x : y) to
                    (-9*x^4 + 18*x^2*y^2 - 2*x*y^3 - 9*y^4 : x^3*y + x*y^3 + y^4)
        """
        if self.dimension_relative() != 1:
            raise TypeError("must be dimension 1")
        if self.base_ring() != E.base_ring():
            E = E.change_ring(self.base_ring())

        L = E.multiplication_by_m(m, x_only=True)
        F = [L.numerator(), L.denominator()]
        R = self.coordinate_ring()
        x, y = R.gens()
        phi = F[0].parent().hom([x], R)
        F = [phi(F[0]).homogenize(y), phi(F[1]).homogenize(y) * y]
        return DynamicalSystem_projective(F, domain=self)

    def cartesian_product(self, other):
        r"""
        Return the Cartesian product of this projective space and
        ``other``.

        INPUT:

        - ``other`` -- a projective space with the same base ring as this space

        OUTPUT: a Cartesian product of projective spaces

        EXAMPLES::

            sage: P1 = ProjectiveSpace(QQ, 1, 'x')
            sage: P2 = ProjectiveSpace(QQ, 2, 'y')
            sage: PP = P1.cartesian_product(P2); PP
            Product of projective spaces P^1 x P^2 over Rational Field
            sage: PP.gens()
            (x0, x1, y0, y1, y2)
        """
        return ProductProjectiveSpaces([self, other])

    def chebyshev_polynomial(self, n, kind='first', monic=False):
        """
        Generates an endomorphism of this projective line by a Chebyshev polynomial.

        Chebyshev polynomials are a sequence of recursively defined orthogonal
        polynomials. Chebyshev of the first kind are defined as `T_0(x) = 1`,
        `T_1(x) = x`, and `T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)`. Chebyshev of
        the second kind are defined as `U_0(x) = 1`,
        `U_1(x) = 2x`, and `U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x)`.

        INPUT:

        - ``n`` -- nonnegative integer

        - ``kind`` -- ``'first'`` (default) or ``'second'`` specifying which
          kind of Chebyshev the user would like to generate

        - ``monic`` -- boolean (default: ``False``) specifying if the
          polynomial defining the system should be monic or not

        OUTPUT: :class:`DynamicalSystem_projective`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(5, 'first')                                    # needs sage.symbolic
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (16*x^5 - 20*x^3*y^2 + 5*x*y^4 : y^5)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(3, 'second')                                   # needs sage.symbolic
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (8*x^3 - 4*x*y^2 : y^3)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(3, 2)                                          # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: keyword 'kind' must have a value of either 'first' or 'second'

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(-4, 'second')
            Traceback (most recent call last):
            ...
            ValueError: first parameter 'n' must be a nonnegative integer

        ::

            sage: P = ProjectiveSpace(QQ, 2, 'x')
            sage: P.chebyshev_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: projective space must be of dimension 1

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(3, monic=True)                                 # needs sage.symbolic
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 - 3*x*y^2 : y^3)

        ::

            sage: F.<t> = FunctionField(QQ)
            sage: P.<y,z> = ProjectiveSpace(F, 1)
            sage: P.chebyshev_polynomial(4, monic=True)                                 # needs sage.symbolic
            Dynamical System of Projective Space of dimension 1
             over Rational function field in t over Rational Field
              Defn: Defined on coordinates by sending (y : z) to
                    (y^4 + (-4)*y^2*z^2 + 2*z^4 : z^4)
        """
        if self.dimension_relative() != 1:
            raise TypeError("projective space must be of dimension 1")
        n = ZZ(n)
        if (n < 0):
            raise ValueError("first parameter 'n' must be a nonnegative integer")
        # use the affine version and then homogenize.
        A = self.affine_patch(1)
        f = A.chebyshev_polynomial(n, kind)
        if monic and self.base().characteristic() != 2:
            f = f.homogenize(1)
            return f.conjugate(matrix([[~ZZ(2), 0], [0, 1]]))
        return f.homogenize(1)

    def veronese_embedding(self, d, CS=None, order='lex'):
        r"""
        Return the degree ``d`` Veronese embedding from this projective space.

        INPUT:

        - ``d`` -- positive integer

        - ``CS`` -- (default: ``None``) a projective ambient space to embed
          into. If this projective space has dimension `N`, the dimension of
          ``CS`` must be `\binom{N + d}{d} - 1`. This is constructed if not
          specified.

        - ``order`` -- string (default: ``'lex'``); a monomial order to use to
          arrange the monomials defining the embedding. The monomials will be
          arranged from greatest to least with respect to this order.

        OUTPUT: a scheme morphism from this projective space to ``CS``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: vd = P.veronese_embedding(4, order='invlex')                          # needs sage.combinat
            sage: vd                                                                    # needs sage.combinat
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 4 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (y^4 : x*y^3 : x^2*y^2 : x^3*y : x^4)

        Veronese surface::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: Q.<q,r,s,t,u,v> = ProjectiveSpace(QQ, 5)
            sage: vd = P.veronese_embedding(2, Q)                                       # needs sage.combinat
            sage: vd                                                                    # needs sage.combinat
            Scheme morphism:
              From: Projective Space of dimension 2 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 : x*y : x*z : y^2 : y*z : z^2)
            sage: vd(P.subscheme([]))                                                   # needs sage.combinat sage.libs.singular
            Closed subscheme of Projective Space of dimension 5 over Rational Field
             defined by:
              -u^2 + t*v,
              -s*u + r*v,
              -s*t + r*u,
              -s^2 + q*v,
              -r*s + q*u,
              -r^2 + q*t
        """
        d = ZZ(d)
        if d <= 0:
            raise ValueError("(=%s) must be a positive integer" % d)
        N = self.dimension()
        # construct codomain space if not given
        if CS is None:
            CS = ProjectiveSpace(self.base_ring(), binomial(N + d, d) - 1)
        else:
            if not isinstance(CS, ProjectiveSpace_ring):
                raise TypeError("(=%s) must be a projective space" % CS)
            if CS.dimension() != binomial(N + d, d) - 1:
                raise TypeError("(=%s) has the wrong dimension to serve as the codomain space" % CS)

        R = self.coordinate_ring().change_ring(order=order)
        monomials = sorted([R({tuple(v): 1}) for v in WeightedIntegerVectors(d, [1] * (N + 1))])
        monomials.reverse()  # order the monomials greatest to least via the given monomial order
        return Hom(self, CS)(monomials)

    def point_transformation_matrix(self, points_source, points_target, normalize=True):
        r"""
        Returns a unique element of PGL that transforms one set of points to another.

        Given a projective space of dimension n and a set of n+2 source points and a set of n+2 target
        points in the same projective space, such that no n+1 points of each set are linearly dependent
        find the unique element of PGL that translates the source points to the target points.

        .. warning::
            over non-exact rings such as the ComplexField, the returned matrix could
            be very far from correct.

        INPUT:

        - ``points_source`` -- points in source projective space

        - ``points_target`` -- points in target projective space

        - ``normalize`` -- boolean (default: ``True``); if the returned matrix
          should be normalized. Only works over exact rings. If the base ring
          is a field, the matrix is normalized so that the last nonzero entry
          in the last row is 1. If the base ring is a ring, then the matrix is
          normalized so that the entries are elements of the base ring.

        OUTPUT: transformation matrix - element of PGL

        ALGORITHM:

        See [Hutz2007]_, Proposition 2.16 for details.

        EXAMPLES::

            sage: P1.<a,b,c> = ProjectiveSpace(QQ, 2)
            sage: points_source = [P1([1, 4, 1]), P1([1, 2, 2]), P1([3, 5, 1]), P1([1, -1, 1])]
            sage: points_target = [P1([5, -2, 7]), P1([3, -2, 3]), P1([6, -5, 9]), P1([3, 6, 7])]
            sage: m = P1.point_transformation_matrix(points_source, points_target); m   # needs sage.modules
            [ -13/59 -128/59  -25/59]
            [538/177    8/59  26/177]
            [ -45/59 -196/59       1]
            sage: [m*points_source[i] == points_target[i] for i in range(4)]            # needs sage.modules
            [True, True, True, True]

        ::

            sage: P.<a,b> = ProjectiveSpace(GF(13),  1)
            sage: points_source = [P([-6, 7]), P([1, 4]), P([3, 2])]
            sage: points_target = [P([-1, 2]), P([0, 2]), P([-1, 6])]
            sage: P.point_transformation_matrix(points_source, points_target)           # needs sage.modules
            [10  4]
            [10  1]

        ::

            sage: P.<a,b> = ProjectiveSpace(QQ, 1)
            sage: points_source = [P([-6, -4]), P([1, 4]), P([3, 2])]
            sage: points_target = [P([-1, 2]), P([0, 2]), P([-7, -3])]
            sage: P.point_transformation_matrix(points_source, points_target)           # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: source points not independent

        ::

            sage: R.<t> = FunctionField(QQ)
            sage: P.<a,b> = ProjectiveSpace(R, 1)
            sage: points_source = [P([-6*t, 7]), P([1, 4]), P([3, 2])]
            sage: points_target = [P([-1, 2*t]), P([0, 2]), P([-1, 6])]
            sage: P.point_transformation_matrix(points_source, points_target)           # needs sage.modules
            [             (1/3*t + 7/12)/(t^2 - 53/24*t)       (-1/12*t - 7/48)/(t^2 - 53/24*t)]
            [(-2/3*t^2 - 7/36*t - 35/12)/(t^2 - 53/24*t)                                      1]

        ::

            sage: P1.<a,b,c> = ProjectiveSpace(RR, 2)
            sage: points_source = [P1([1, 4, 1]), P1([1, 2, 2]), P1([3, 5, 1]), P1([1, -1, 1])]
            sage: points_target = [P1([5, -2, 7]), P1([3, -2, 3]), P1([6, -5, 9]), P1([3, 6, 7])]
            sage: P1.point_transformation_matrix(points_source,        # abs tol 1e-13  # needs sage.modules
            ....:                                points_target)
            [-0.0619047619047597  -0.609523809523810  -0.119047619047621]
            [  0.853968253968253  0.0380952380952380  0.0412698412698421]
            [ -0.214285714285712  -0.933333333333333   0.280952380952379]

        ::

            sage: P1.<a,b,c> = ProjectiveSpace(ZZ, 2)
            sage: points_source = [P1([1, 4, 1]), P1([1, 2, 2]), P1([3, 5, 1]), P1([1, -1, 1])]
            sage: points_target = [P1([5, -2, 7]), P1([3, -2, 3]), P1([6, -5, 9]), P1([3, 6, 7])]
            sage: P1.point_transformation_matrix(points_source, points_target)          # needs sage.modules
            [ -39 -384  -75]
            [ 538   24   26]
            [-135 -588  177]

        ::

            sage: P1.<a,b,c> = ProjectiveSpace(ZZ, 2)
            sage: points_source = [P1([1, 4, 1]), P1([1, 2, 2]), P1([3, 5, 1]), P1([1, -1, 1])]
            sage: points_target = [P1([5, -2, 7]), P1([3, -2, 3]), P1([6, -5, 9]), P1([3, 6, 7])]
            sage: P1.point_transformation_matrix(points_source, points_target,          # needs sage.modules
            ....:                                normalize=False)
            [-13/30 -64/15   -5/6]
            [269/45   4/15  13/45]
            [  -3/2 -98/15  59/30]

        ::

            sage: R.<t> = ZZ[]
            sage: P.<a,b> = ProjectiveSpace(R, 1)
            sage: points_source = [P([-6*t, 7]), P([1, 4]), P([3, 2])]
            sage: points_target = [P([-1, 2*t]), P([0, 2]), P([-1, 6])]
            sage: P.point_transformation_matrix(points_source, points_target)           # needs sage.modules
            [         -48*t - 84           12*t + 21]
            [96*t^2 + 28*t + 420    -144*t^2 + 318*t]

        TESTS::

            sage: P.<a,b> = ProjectiveSpace(QQ, 1)
            sage: points_source = [P([-6, -1]), P([1, 4]), P([3, 2])]
            sage: points_target = [P([-1, 2]), P([0, 2]), P([-2, 4])]
            sage: P.point_transformation_matrix(points_source, points_target)           # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: target points not independent

        ::

            sage: P.<a,b,c> = ProjectiveSpace(QQ, 2)
            sage: points_source = [P([1, 4, 1]), P([2, -7, 9]), P([3, 5, 1])]
            sage: points_target = [P([5, -2, 7]), P([3, -2, 3]), P([6, -5, 9]), P([6, -1, 1])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: incorrect number of points in source, need 4 points

        ::

            sage: P.<a,b,c> = ProjectiveSpace(QQ, 2)
            sage: points_source = [P([1, 4, 1]), P([2, -7, 9]), P([3, 5, 1]), P([1, -1, 1])]
            sage: points_target = [P([5, -2, 7]), P([3, -2, 3]), P([6, -5, 9]), P([6, -1, 1]), P([7, 8, -9])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: incorrect number of points in target, need 4 points

        ::

            sage: P.<a,b,c> = ProjectiveSpace(QQ, 2)
            sage: P1.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: points_source = [P([1, 4, 1]), P([2, -7, 9]), P([3, 5, 1]), P1([1, -1, 1])]
            sage: points_target=[P([5, -2, 7]), P([3, -2, 3]), P([6, -5, 9]), P([6, -1, 1])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: source points not in self

        ::

            sage: P.<a,b,c> = ProjectiveSpace(QQ, 2)
            sage: P1.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: points_source = [P([1, 4, 1]), P([2, -7, 9]), P([3, 5, 1]), P([1, -1, 1])]
            sage: points_target = [P([5, -2, 7]), P([3, -2, 3]), P([6, -5, 9]), P1([6, -1, 1])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: target points not in self

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: points_source = [P(1, 0, 0), P(0, 1, 0), P(0, 0, 1), P(1, -1, -1)]
            sage: points_target = [P(0, 1, 0), P(-2, 0, 1), P(0, 0, 1), P(1, -1, -1)]
            sage: P.point_transformation_matrix(points_source, points_target,           # needs sage.modules
            ....:                               normalize=True)
            [ 0 -2  0]
            [-2  0  0]
            [ 0  1  1]
        """
        r = self.base_ring()
        n = self.dimension_relative()
        # makes sure there aren't to few or two many points
        if len(points_source) != n + 2:
            raise ValueError("incorrect number of points in source, need %d points" % (n + 2))
        if len(points_target) != n + 2:
            raise ValueError("incorrect number of points in target, need %d points" % (n + 2))
        if any(x.codomain() != self for x in points_source):
            raise ValueError("source points not in self")
        if any(x.codomain() != self for x in points_target):
            raise ValueError("target points not in self")
        Ms = matrix(r, [list(s) for s in points_source])
        if any(m == 0 for m in Ms.minors(n + 1)):
            raise ValueError("source points not independent")
        Mt = matrix(r, [list(t) for t in points_target])
        if any(l == 0 for l in Mt.minors(n + 1)):
            raise ValueError("target points not independent")

        # get_matrix calculates the transform from the list of points
        # [ [1 : 0 : 0 : ... ]
        #   [0 : 1 : 0 : ... ]
        #   [0 : 0 : 1 : ... ]
        #   ...
        #   [1 : 1 : 1 : ... ] ]
        # to the list of points S
        def get_matrix(S, N):
            a = matrix(N+1, N+1, [S[j][i] for i in range(N+1) for j in range(N+1)])
            b = matrix(N+1, 1, list(S[N+1]))
            X = a.solve_right(b)
            m = matrix(N+1, N+1, [X[i,0]*S[i][j] for i in range(N+1) for j in range(N+1)])
            m = m.transpose()
            return m

        m_source = get_matrix(points_source, n)
        m_target = get_matrix(points_target, n)
        return_mat = m_target*m_source.inverse()
        if normalize:
            R = self.base_ring()
            if R.is_exact():
                if R.is_field():
                    last_row = list(return_mat.rows()[-1])[:]
                    last_ele = last_row.pop()
                    while last_ele == 0:
                        last_ele = last_row.pop()
                    return_mat *= ZZ(1)/last_ele
                else:
                    lcm = return_mat[0][0].denominator()
                    for row in return_mat.rows():
                        for ele in row:
                            lcm = lcm.lcm(ele.denominator())
                    return_mat *= lcm
        return return_mat

    def hyperplane_transformation_matrix(self, plane_1, plane_2):
        r"""
        Return a PGL element sending ``plane_1`` to ``plane_2``.

        ``plane_1`` and ``plane_2`` must be hyperplanes (subschemes of
        codimension 1, each defined by a single linear homogeneous equation).

        INPUT:

        - ``plane_1``, ``plane_2`` -- hyperplanes of this projective space

        OUTPUT: an element of PGL

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: plane1 = P.subscheme(x)
            sage: plane2 = P.subscheme(y)
            sage: m = P.hyperplane_transformation_matrix(plane1, plane2); m             # needs sage.modules
            [0 1]
            [1 0]
            sage: plane2(m*P((0,1)))                                                    # needs sage.modules
            (1 : 0)

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: plane1 = P.subscheme(x + 2*y + z)
            sage: plane2 = P.subscheme(2*x + y + z)
            sage: P.hyperplane_transformation_matrix(plane1, plane2)                    # needs sage.modules
            [1 0 0 0]
            [0 4 0 0]
            [0 0 2 0]
            [0 0 0 1]

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: plane1 = P.subscheme(x + y)
            sage: plane2 = P.subscheme(y)
            sage: P.hyperplane_transformation_matrix(plane1, plane2)                    # needs sage.modules
            [-1  0]
            [ 1  1]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = CyclotomicField(3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: plane1 = P.subscheme(x - 2*v*y + z)
            sage: plane2 = P.subscheme(x + v*y + v*z)
            sage: m = P.hyperplane_transformation_matrix(plane1, plane2); m             # needs sage.modules
            [   v    0    0]
            [   0 -2*v    0]
            [   0    0    1]

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<k> = NumberField(x^2 + 1)
            sage: P.<x,y,z,w> = ProjectiveSpace(K, 3)
            sage: plane1 = P.subscheme(k*x + 2*k*y + z)
            sage: plane2 = P.subscheme(7*k*x + y + 9*z)
            sage: m = P.hyperplane_transformation_matrix(plane1, plane2); m             # needs sage.modules
            [   1    0    0    0]
            [   0 14*k    0    0]
            [   0    0  7/9    0]
            [   0    0    0    1]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = CyclotomicField(3)
            sage: R.<t> = K[]
            sage: F.<w> = K.extension(t^5 + 2)
            sage: G.<u> = F.absolute_field()
            sage: P.<x,y,z> = ProjectiveSpace(G, 2)
            sage: plane1 = P.subscheme(x - 2*u*y + z)
            sage: plane2 = P.subscheme(x + u*y + z)
            sage: m = P.hyperplane_transformation_matrix(plane1, plane2)                # needs sage.modules
            sage: plane2(m*P((2*u, 1, 0)))                                              # needs sage.modules
            (-u : 1 : 0)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(FiniteField(2), 2)
            sage: plane1 = P.subscheme(x + y + z)
            sage: plane2 = P.subscheme(z)
            sage: P.hyperplane_transformation_matrix(plane1, plane2)                    # needs sage.modules
            [1 0 0]
            [1 1 0]
            [1 1 1]

        ::

            sage: R.<t> = QQ[]
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: plane1 = P.subscheme(x + 9*t*y + z)
            sage: plane2 = P.subscheme(x + z)
            sage: P.hyperplane_transformation_matrix(plane1, plane2)                    # needs sage.modules
            [  1 9*t   0]
            [  1   0   0]
            [  0   0   1]

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: plane1 = P.subscheme(x^2)
            sage: plane2 = P.subscheme(y)
            sage: P.hyperplane_transformation_matrix(plane1, plane2)
            Traceback (most recent call last):
            ...
            ValueError: plane_1 must be defined by a single degree 1 equation
        """
        if not isinstance(plane_1, AlgebraicScheme_subscheme_projective):
            raise TypeError('plane_1 must be a subscheme')
        if not isinstance(plane_2, AlgebraicScheme_subscheme_projective):
            raise TypeError('plane_2 must be a subscheme')
        if plane_1.ambient_space() != self:
            raise ValueError('plane_1 must be a subscheme of this projective space')
        if plane_2.ambient_space() != self:
            raise ValueError('plane_2 must be a subscheme of this projective space')
        if len(plane_1.defining_polynomials()) > 1 or plane_1.defining_polynomials()[0].degree() != 1:
            raise ValueError('plane_1 must be defined by a single degree 1 equation')
        if len(plane_2.defining_polynomials()) > 1 or plane_2.defining_polynomials()[0].degree() != 1:
            raise ValueError('plane_2 must be defined by a single degree 1 equation')
        N = self.dimension_relative()
        CR = self.coordinate_ring()
        points = []
        from sage.rings.rational_field import QQ
        P_QQ = ProjectiveSpace(QQ, N)
        # to determine the PGL transform, we need N+2 points source points and N+2 target points,
        # of which no N+1 are co-planar. Additionally, in order to map plane_1 to plane_2, N source
        # points must lie on plane_1, and N target points must lie on plane_2
        for plane in [plane_1, plane_2]:
            source_points = []
            nonzero_places = []
            height_1 = P_QQ.points_of_bounded_height(bound=1)

            # first we find N planar points
            # we have a single linear equation with N+1 variables
            # first we add a point for each variable with coefficient 0
            # giving us J points added in this loop
            for i in range(N+1):
                if plane.defining_polynomials()[0].coefficient(CR.gens()[i]) == 0:
                    L = [0]*(N+1)
                    L[i] = 1
                    source_points.append(self(L))
                else:
                    nonzero_places.append(i)
            # next we add a point for each variable with nonzero coefficient, except the last
            # giving us a total of (N+1) - J - 1 = N - J points added in this loop
            # resulting in exactly J + (N-J) = N points on the plane
            for i in range(len(nonzero_places)-1):
                nonzero_place1 = nonzero_places[i]
                nonzero_place2 = nonzero_places[i+1]
                L = [0]*(N+1)
                L[nonzero_place1] = -1*plane.defining_polynomials()[0].coefficient(CR.gens()[nonzero_place2])
                L[nonzero_place2] = plane.defining_polynomials()[0].coefficient(CR.gens()[nonzero_place1])
                source_points.append(self(L))

            # next we add independent points until we have N+2 points total
            for point in height_1:
                if len(source_points) == N:
                    try:
                        plane(point)
                    except (ValueError, TypeError):
                        source_points.append(self(point))
                        base_list = [list(s) for s in source_points]
                elif len(source_points) == N + 1:
                    Ms = matrix(base_list + [point.change_ring(self.base_ring())])
                    if not any(m == 0 for m in Ms.minors(N + 1)):
                        source_points.append(self(point))
                        break
            if len(source_points) != N+2:
                raise NotImplementedError('Failed to automatically find sufficient independent points.' +
                    ' Please find the necessary independent points manually, then use point transformation matrix.')
            points.append(source_points)
        return self.point_transformation_matrix(points[0], points[1])

    def is_linearly_independent(self, points, n=None):
        r"""
        Return whether the set of points is linearly independent.

        Alternatively, specify ``n`` to check if every subset of
        size ``n`` is linearly independent.

        INPUT:

        - ``points`` -- list of points in this projective space

        - ``n`` -- (optional) positive integer less than or equal to the length
          of ``points``. Specifies the size of the subsets to check for
          linear independence.

        OUTPUT: ``True`` if ``points`` is linearly independent, ``False`` otherwise

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: points = [P((1, 0, 1)), P((1, 2, 1)), P((1, 3, 4))]
            sage: P.is_linearly_independent(points)                                     # needs sage.modules
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: points = [P((1, 0, 1)), P((1, 2, 1)), P((1, 3, 4)), P((0, 0, 1))]
            sage: P.is_linearly_independent(points, 2)                                  # needs sage.modules
            True

        ::

            sage: R.<c> = QQ[]
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: points = [P((c, 0, 1)), P((0, c, 1)), P((1, 0, 4)), P((0, 0, 1))]
            sage: P.is_linearly_independent(points, 3)                                  # needs sage.modules
            False

        ::

            sage: R.<c> = QQ[]
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R), 2)
            sage: points = [P((c, 0, 1)), P((0, c, 1)), P((1, 3, 4)), P((0, 0, 1))]
            sage: P.is_linearly_independent(points, 3)                                  # needs sage.modules
            True

        ::

            sage: # needs sage.rings.number_field
            sage: K.<k> = CyclotomicField(3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: points = [P((k, k^2, 1)), P((0, k, 1)), P((1, 0, 4)), P((0, 0, 1))]
            sage: P.is_linearly_independent(points, 3)                                  # needs sage.modules
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: points = [P((1, 0)), P((1, 1))]
            sage: P.is_linearly_independent(points)                                     # needs sage.modules
            True

        TESTS::

            sage: points = [P(1, 0), P(1, 1), P(2, 1)]
            sage: P.is_linearly_independent(points, 5)
            Traceback (most recent call last):
            ...
            ValueError: n must be a nonnegative integer not greater than the length of points
        """
        if not isinstance(points, list):
            raise TypeError("points must be a list")
        if any(not isinstance(point, SchemeMorphism_point_projective_ring) for point in points):
            raise TypeError("points must be a list of projective points")
        if any(x.codomain() != self for x in points):
            raise ValueError("points not in this projective space")
        if n is None:
            M = matrix([list(t) for t in points])
            return M.rank() == len(points)
        n = Integer(n)
        if n < 1 or n > len(points):
            raise ValueError('n must be a nonnegative integer not greater than the length of points')
        all_subsets = Subsets(range(len(points)), n)
        linearly_independent = True
        for subset in all_subsets:
            point_list = []
            for index in subset:
                point_list.append(list(points[index]))
            M = matrix(point_list)
            if M.rank() != n:
                linearly_independent = False
                break
        return linearly_independent


class ProjectiveSpace_field(ProjectiveSpace_ring):
    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: P2._point_homset(Spec(GF(3)), P2)
            Set of rational points of Projective Space of dimension 2 over Finite Field of size 3
        """
        return SchemeHomset_points_projective_field(*args, **kwds)

    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: point_homset = P2._point_homset(Spec(GF(3)), P2)
            sage: P2._point(point_homset, [1,2,3])
            (2 : 1 : 0)
        """
        return SchemeMorphism_point_projective_field(*args, **kwds)

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: P2._morphism(P2.Hom(P2), [x,y,z])
            Scheme endomorphism of Projective Space of dimension 2 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x : y : z)
        """
        return SchemeMorphism_polynomial_projective_space_field(*args, **kwds)

    def subscheme_from_Chow_form(self, Ch, dim):
        r"""
        Return the subscheme defined by the Chow equations associated to the Chow form ``Ch``.

        These equations define the subscheme set-theoretically, but only for smooth
        subschemes and hypersurfaces do they define the subscheme as a scheme.

        ALGORITHM:

        The Chow form is a polynomial in the Plucker coordinates. The Plucker coordinates
        are the bracket polynomials. We first re-write the Chow form in terms of the dual
        Plucker coordinates. Then we expand `Ch(span(p,L)` for a generic point `p` and a
        generic linear subspace `L`. The coefficients as polynomials in the coordinates
        of `p` are the equations defining the subscheme. [DalbecSturmfels].

        INPUT:

        - ``Ch`` -- a homogeneous polynomial

        - ``dim`` -- the dimension of the associated scheme

        OUTPUT: a projective subscheme

        EXAMPLES::

            sage: P = ProjectiveSpace(QQ, 4, 'z')
            sage: R.<x0,x1,x2,x3,x4> = PolynomialRing(QQ)
            sage: H = x1^2 + x2^2 + 5*x3*x4
            sage: P.subscheme_from_Chow_form(H, 3)                                      # needs sage.modules
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
              -5*z0*z1 + z2^2 + z3^2

        ::

            sage: P = ProjectiveSpace(QQ, 3, 'z')
            sage: R.<x0,x1,x2,x3,x4,x5> = PolynomialRing(QQ)
            sage: H = x1 - x2 - x3 + x5 + 2*x0
            sage: P.subscheme_from_Chow_form(H, 1)                                      # needs sage.modules
            Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
              -z1 + z3,
              z0 + z2 + z3,
              -z1 - 2*z3,
              -z0 - z1 + 2*z2

        ::

            sage: # needs sage.libs.singular
            sage: P.<x0,x1,x2,x3> = ProjectiveSpace(GF(7), 3)
            sage: X = P.subscheme([x3^2 + x1*x2, x2 - x0])
            sage: Ch = X.Chow_form(); Ch
            t0^2 - 2*t0*t3 + t3^2 - t2*t4 - t4*t5
            sage: Y = P.subscheme_from_Chow_form(Ch, 1); Y
            Closed subscheme of Projective Space of dimension 3
             over Finite Field of size 7 defined by:
              x1*x2 + x3^2,
              -x0*x2 + x2^2,
              -x0*x1 - x1*x2 - 2*x3^2,
              x0^2 - x0*x2,
              x0*x1 + x3^2,
              -2*x0*x3 + 2*x2*x3,
              2*x0*x3 - 2*x2*x3,
              x0^2 - 2*x0*x2 + x2^2
            sage: I = Y.defining_ideal()
            sage: I.saturation(I.ring().ideal(list(I.ring().gens())))[0]
            Ideal (x0 - x2, x1*x2 + x3^2) of Multivariate Polynomial Ring
             in x0, x1, x2, x3 over Finite Field of size 7
        """
        if not Ch.is_homogeneous():
            raise ValueError("Chow form must be a homogeneous polynomial")
        n = self.dimension_relative()
        R = Ch.parent()
        if binomial(n + 1, n - dim) != R.ngens():
            raise ValueError("for given dimension, there should be %d variables in the Chow form" % binomial(n + 1, n - dim))
        # create the brackets associated to variables
        L1 = []
        for t in UnorderedTuples(list(range(n + 1)), dim + 1):
            if all(t[i] < t[i + 1] for i in range(dim)):
                L1.append(list(t))
        # create the dual brackets
        L2 = []
        signs = []
        for l in L1:
            s = []
            for v in range(n + 1):
                if v not in l:
                    s.append(v)
            t1 = [b + 1 for b in l]
            t2 = [b + 1 for b in s]
            perm = Permutation(t1 + t2)
            signs.append(perm.sign())
            L2.append(s)
        # create the polys associated to dual brackets
        if n - dim - 1 > 0:
            S = PolynomialRing(R.base_ring(), n + 1, 'z')
            T = PolynomialRing(S, (n + 1) * (n - dim - 1), 's')
            M = matrix(T, n - dim, n + 1, list(S.gens()) + list(T.gens()))
        else:
            T = PolynomialRing(R.base_ring(), n + 1, 'z')
            M = matrix(T, n - dim, n + 1, list(T.gens()))
        coords = []
        for i in range(len(L2)):
            coords.append(signs[i] * M.matrix_from_columns(L2[i]).det())
        # substitute in dual brackets to chow form
        phi = R.hom(coords, T)
        ch = phi(Ch)
        # coefficients are polys in zs which are the chow equations for the chow form
        if n - dim - 1 > 0:
            return self.subscheme(ch.coefficients())
        else:
            return self.subscheme(ch)

    def curve(self, F):
        r"""
        Return a curve defined by ``F`` in this projective space.

        INPUT:

        - ``F`` -- a polynomial, or a list or tuple of polynomials in
          the coordinate ring of this projective space

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P.curve([y^2 - x*z])                                                  # needs sage.schemes
            Projective Plane Curve over Rational Field defined by y^2 - x*z
        """
        from sage.schemes.curves.constructor import Curve
        return Curve(F, self)

    def line_through(self, p, q):
        """
        Return the line through ``p`` and ``q``.

        INPUT:

        - ``p``, ``q`` -- distinct rational points of the projective space

        EXAMPLES::

            sage: P3.<x0,x1,x2,x3> = ProjectiveSpace(3, QQ)
            sage: p1 = P3(1, 2, 3, 4)
            sage: p2 = P3(4, 3, 2, 1)
            sage: P3.line_through(p1, p2)                                               # needs sage.libs.singular sage.schemes
            Projective Curve over Rational Field defined by
              -5/4*x0 + 5/2*x1 - 5/4*x2,        -5/2*x0 + 15/4*x1 - 5/4*x3,
              -5/4*x0 + 15/4*x2 - 5/2*x3,       -5/4*x1 + 5/2*x2 - 5/4*x3
            sage: p3 = P3(2,4,6,8)
            sage: P3.line_through(p1, p3)
            Traceback (most recent call last):
            ...
            ValueError: not distinct points
        """
        if p == q:
            raise ValueError("not distinct points")

        from sage.schemes.curves.constructor import Curve

        m = matrix(3, list(self.gens()) + list(p) + list(q))
        return Curve([f for f in m.minors(3) if f])


class ProjectiveSpace_finite_field(ProjectiveSpace_field):
    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: point_homset = P2._point_homset(Spec(GF(3)), P2)
            sage: P2._point(point_homset, [1,2,3])
            (2 : 1 : 0)
        """
        return SchemeMorphism_point_projective_finite_field(*args, **kwds)

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: P2._morphism(P2.Hom(P2), [x,y,z])
            Scheme endomorphism of Projective Space of dimension 2 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x : y : z)
        """
        return SchemeMorphism_polynomial_projective_space_finite_field(*args, **kwds)

    def __iter__(self):
        r"""
        Return iterator over the elements of this projective space.

        Note that iteration is over the decomposition
        `\mathbb{P}^n = \mathbb{A}^n \cup \mathbb{P}^n-1`, where
        `\mathbb{A}^n` is the `n`-th affine patch and
        `\mathbb{P}^n-1` is the hyperplane at infinity
        `x_n = 0`.

        EXAMPLES::

            sage: FF = FiniteField(3)
            sage: PP = ProjectiveSpace(0, FF)
            sage: [ x for x in PP ]
            [(1)]
            sage: PP = ProjectiveSpace(1, FF)
            sage: [ x for x in PP ]
            [(0 : 1), (1 : 1), (2 : 1), (1 : 0)]
            sage: PP = ProjectiveSpace(2, FF)
            sage: [ x for x in PP ]
            [(0 : 0 : 1),
             (0 : 1 : 1),
             (0 : 2 : 1),
             (1 : 0 : 1),
             (1 : 1 : 1),
             (1 : 2 : 1),
             (2 : 0 : 1),
             (2 : 1 : 1),
             (2 : 2 : 1),
             (0 : 1 : 0),
             (1 : 1 : 0),
             (2 : 1 : 0),
             (1 : 0 : 0)]

        AUTHORS:

        - David Kohel, John Cremona

        .. TODO::

            Iteration for point sets over finite fields, and return of
            iter of point set over base field. Note that the point set does not
            know whether this is a projective space or subscheme.
        """
        n = self.dimension_relative()
        R = self.base_ring()
        zero = (R.zero(), )
        one = (R.one(), )
        PHom = self.point_homset()
        C = PHom.codomain()

        for k in range(n + 1): # position of last 1 before the 0's
            for v in product(*[R for _ in range(n - k)]):
                yield C._point(PHom, v + one + zero * k, check=False)

    def rational_points(self, F=None):
        """
        Return the list of ``F``-rational points on this projective space,
        where ``F`` is a given finite field, or the base ring of this space.

        EXAMPLES::

            sage: P = ProjectiveSpace(1, GF(3))
            sage: P.rational_points()
            [(0 : 1), (1 : 1), (2 : 1), (1 : 0)]
            sage: sorted(P.rational_points(GF(3^2, 'b')), key=str)                      # needs sage.rings.finite_rings
            [(0 : 1), (1 : 0), (1 : 1), (2 : 1),
             (2*b + 1 : 1), (2*b + 2 : 1), (2*b : 1),
             (b + 1 : 1), (b + 2 : 1), (b : 1)]
        """
        if F is None:
            return list(self)
        elif not isinstance(F, FiniteField):
            raise TypeError("second argument (= %s) must be a finite field" % F)
        return list(self.base_extend(F))

    def rational_points_dictionary(self):
        r"""
        Return dictionary of points.

        OUTPUT: dictionary

        EXAMPLES::

            sage: P1 = ProjectiveSpace(GF(7), 1, 'x')
            sage: P1.rational_points_dictionary()
            {(0 : 1): 0,
             (1 : 0): 7,
             (1 : 1): 1,
             (2 : 1): 2,
             (3 : 1): 3,
             (4 : 1): 4,
             (5 : 1): 5,
             (6 : 1): 6}
        """
        n = self.dimension_relative()
        R = self.base_ring()
        D = {}
        zero = R.zero()
        i = n
        index = 0
        while not i < 0:
            P = [zero for _ in range(i)] + [R.one()]
            P += [zero for _ in range(n - i)]
            D.update({self(P): index})
            index += 1
            iters = [iter(R) for _ in range(i)]
            for x in iters:
                next(x)  # put at zero
            j = 0
            while j < i:
                try:
                    P[j] = next(iters[j])
                    D.update({self(P): index})
                    index += 1
                    j = 0
                except StopIteration:
                    iters[j] = iter(R)  # reset
                    next(iters[j])  # put at zero
                    P[j] = zero
                    j += 1
            i -= 1
        return D


class ProjectiveSpace_rational_field(ProjectiveSpace_field):
    def rational_points(self, bound=0):
        r"""
        Return the projective points `(x_0:\cdots:x_n)` over
        `\QQ` with `|x_i| \leq` bound.

        ALGORITHM:

        The very simple algorithm works as follows: every point
        `(x_0:\cdots:x_n)` in projective space has a unique
        largest index `i` for which `x_i` is not
        zero. The algorithm then iterates downward on this
        index. We normalize by choosing `x_i` positive. Then,
        the points `x_0,\ldots,x_{i-1}` are the points of
        affine `i`-space that are relatively prime to
        `x_i`. We access these by using the Tuples method.

        INPUT:

        - ``bound`` -- integer

        EXAMPLES::

            sage: PP = ProjectiveSpace(0, QQ)
            sage: PP.rational_points(1)
            [(1)]
            sage: PP = ProjectiveSpace(1, QQ)
            sage: PP.rational_points(2)
            [(-2 : 1), (-1 : 1), (0 : 1), (1 : 1), (2 : 1), (-1/2 : 1), (1/2 : 1), (1 : 0)]
            sage: PP = ProjectiveSpace(2, QQ)
            sage: PP.rational_points(2)
            [(-2 : -2 : 1), (-1 : -2 : 1), (0 : -2 : 1), (1 : -2 : 1), (2 : -2 : 1),
             (-2 : -1 : 1), (-1 : -1 : 1), (0 : -1 : 1), (1 : -1 : 1), (2 : -1 : 1),
             (-2 : 0 : 1), (-1 : 0 : 1), (0 : 0 : 1), (1 : 0 : 1), (2 : 0 : 1), (-2 : 1 : 1),
             (-1 : 1 : 1), (0 : 1 : 1), (1 : 1 : 1), (2 : 1 : 1), (-2 : 2 : 1),
             (-1 : 2 : 1), (0 : 2 : 1), (1 : 2 : 1), (2 : 2 : 1), (-1/2 : -1 : 1),
             (1/2 : -1 : 1), (-1 : -1/2 : 1), (-1/2 : -1/2 : 1), (0 : -1/2 : 1),
             (1/2 : -1/2 : 1), (1 : -1/2 : 1), (-1/2 : 0 : 1), (1/2 : 0 : 1), (-1 : 1/2 : 1),
             (-1/2 : 1/2 : 1), (0 : 1/2 : 1), (1/2 : 1/2 : 1), (1 : 1/2 : 1), (-1/2 : 1 : 1),
             (1/2 : 1 : 1), (-2 : 1 : 0), (-1 : 1 : 0), (0 : 1 : 0), (1 : 1 : 0),
             (2 : 1 : 0), (-1/2 : 1 : 0), (1/2 : 1 : 0), (1 : 0 : 0)]

        AUTHORS:

        - Benjamin Antieau (2008-01-12)
        """
        if not bound > 0:
            raise ValueError("argument bound (= %s) must be a positive integer")

        n = self.dimension_relative()

        Q = [k - bound for k in range(2 * bound + 1)]  # the affine coordinates
        R = [(k + 1) for k in range(bound)]         # the projective coordinate
        S = [Tuples(Q, (k + 1)) for k in range(n)]
        pts = []

        i = n
        while i > 0:
            P = [0 for _ in range(n + 1)]
            for ai in R:
                P[i] = ai
                for tup in S[i - 1]:
                    if gcd((ai,) + tup) == 1:
                        for j in range(i):
                            P[j] = tup[j]
                        pts.append(self(P))
            i -= 1

        # now do i=0; this is treated as a special case so that
        # we don't have all points (1:0),(2,0),(3,0),etc.
        P = [0 for _ in range(n + 1)]
        P[0] = 1
        pts.append(self(P))
        return pts


# fix the pickles from moving projective_space.py
register_unpickle_override('sage.schemes.generic.projective_space',
                           'ProjectiveSpace_field',
                           ProjectiveSpace_field)

register_unpickle_override('sage.schemes.generic.projective_space',
                           'ProjectiveSpace_rational_field',
                           ProjectiveSpace_rational_field)
