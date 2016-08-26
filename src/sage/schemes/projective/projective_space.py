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
    sage: X/CC
    Projective Space of dimension 5 over Complex Field with 53 bits of precision

The third argument specifies the printing names of the generators of the
homogenous coordinate ring. Using the method `.objgens()` you can obtain both
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

    sage: V = P2.subscheme([x+y+z, x+y-z]); V
    Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
     x + y + z,
     x + y - z
    sage: V.dimension()
    0

AUTHORS:

- Ben Hutz: (June 2012): support for rings

- Ben Hutz (9/2014): added support for Cartesian products

- Rebecca Lauren Miller (March 2016) : added point_transformation_matrix
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.arith.misc import binomial

from sage.rings.all import (PolynomialRing,
                            Integer,
                            ZZ)

from sage.rings.ring import CommutativeRing
from sage.rings.rational_field import is_RationalField
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField

from sage.categories.fields import Fields
_Fields = Fields()

from sage.categories.homset import Hom, End
from sage.categories.number_fields import NumberFields
from sage.categories.map import Map

from sage.misc.all import (latex,
                           prod)
from sage.structure.category_object import normalize_names
from sage.arith.all import gcd, binomial
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.permutation import Permutation
from sage.combinat.tuple import Tuples
from sage.combinat.tuple import UnorderedTuples
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import prepare

from sage.schemes.generic.ambient_space import AmbientSpace
from sage.schemes.projective.projective_homset import (SchemeHomset_points_projective_ring,
                                                       SchemeHomset_points_projective_field)
from sage.schemes.projective.projective_point import (SchemeMorphism_point_projective_ring,
                                                      SchemeMorphism_point_projective_field,
                                                      SchemeMorphism_point_projective_finite_field)
from sage.schemes.projective.projective_morphism import (SchemeMorphism_polynomial_projective_space,
                                                         SchemeMorphism_polynomial_projective_space_field,
                                                         SchemeMorphism_polynomial_projective_space_finite_field)

def is_ProjectiveSpace(x):
    r"""
    Return True if ``x`` is a projective space.

    In other words, if ``x`` is an ambient space `\mathbb{P}^n_R`,
    where `R` is a ring and `n\geq 0` is an integer.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_space import is_ProjectiveSpace
        sage: is_ProjectiveSpace(ProjectiveSpace(5, names='x'))
        True
        sage: is_ProjectiveSpace(ProjectiveSpace(5, GF(9,'alpha'), names='x'))
        True
        sage: is_ProjectiveSpace(Spec(ZZ))
        False
    """
    return isinstance(x, ProjectiveSpace_ring)

def ProjectiveSpace(n, R=None, names='x'):
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

    There is also an projective space associated each polynomial ring.

    ::

        sage: R = GF(7)['x,y,z']
        sage: P = ProjectiveSpace(R); P
        Projective Space of dimension 2 over Finite Field of size 7
        sage: P.coordinate_ring()
        Multivariate Polynomial Ring in x, y, z over Finite Field of size 7
        sage: P.coordinate_ring() is R
        True

    ::

        sage: ProjectiveSpace(3, Zp(5), 'y')
        Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20

    ::

        sage: ProjectiveSpace(2,QQ,'x,y,z')
        Projective Space of dimension 2 over Rational Field

    ::

        sage: PS.<x,y>=ProjectiveSpace(1,CC)
        sage: PS
        Projective Space of dimension 1 over Complex Field with 53 bits of precision

    Projective spaces are not cached, i.e., there can be several with
    the same base ring and dimension (to facilitate gluing
    constructions).
    """
    if is_MPolynomialRing(n) and R is None:
        A = ProjectiveSpace(n.ngens()-1, n.base_ring())
        A._coordinate_ring = n
        return A
    if isinstance(R, (int, long, Integer)):
        n, R = R, n
    if R is None:
        R = ZZ  # default is the integers
    if R in _Fields:
        if is_FiniteField(R):
            return ProjectiveSpace_finite_field(n, R, names)
        if is_RationalField(R):
            return ProjectiveSpace_rational_field(n, R, names)
        else:
            return ProjectiveSpace_field(n, R, names)
    elif isinstance(R, CommutativeRing):
        return ProjectiveSpace_ring(n, R, names)
    else:
        raise TypeError("R (=%s) must be a commutative ring"%R)

class ProjectiveSpace_ring(AmbientSpace):
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
    """
    def __init__(self, n, R=ZZ, names=None):
        """
        Initialization function.

        EXAMPLES::

            sage: ProjectiveSpace(3, Zp(5), 'y')
            Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20
        """
        names = normalize_names(n+1, names)
        AmbientSpace.__init__(self, n, R)
        self._assign_names(names)

    def ngens(self):
        """
        Return the number of generators of this projective space.

        This is the number of variables in the coordinate ring of self.

        EXAMPLES::

            sage: ProjectiveSpace(3, QQ).ngens()
            4
            sage: ProjectiveSpace(7, ZZ).ngens()
            8
        """
        return self.dimension_relative() + 1

    def _check_satisfies_equations(self, v):
        """
        Return True if ``v`` defines a point on the scheme; raise a
        TypeError otherwise.

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
            raise TypeError('the argument v=%s must be a list or tuple'%v)
        n = self.ngens()
        if not len(v) == n:
            raise TypeError('the list v=%s must have %s components'%(v, n))
        R = self.base_ring()
        for coord in v:
            if not coord in R:
                raise TypeError('the components of v=%s must be elements of %s'%(v, R))
        zero = [R(0)]*n
        if v == zero:
            raise TypeError('the zero vector is not a point in projective space')
        return True

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme.

        EXAMPLES::

            sage: ProjectiveSpace(3, GF(19^2,'alpha'), 'abcd').coordinate_ring()
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
                               self.variable_names(), self.dimension_relative()+1)
            return self._coordinate_ring

    def _validate(self, polynomials):
        """
        If ``polynomials`` is a tuple of valid polynomial functions on self,
        return ``polynomials``, otherwise raise TypeError.

        Since this is a projective space, polynomials must be homogeneous.

        INPUT:

        - ``polynomials`` -- tuple of polynomials in the coordinate ring of
            this space.

        OUTPUT:

        - tuple of polynomials in the coordinate ring of this space.

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
            raise TypeError('the argument polynomials=%s must be a list or tuple'%polynomials)
        for f in polynomials:
            if not f.is_homogeneous():
                raise TypeError("%s is not a homogeneous polynomial" % f)
        return polynomials

    def __cmp__(self, right):
        """
        Compare equality of two projective spaces.

        EXAMPLES::

            sage: ProjectiveSpace(QQ, 3, 'a') == ProjectiveSpace(ZZ, 3, 'a')
            False
            sage: ProjectiveSpace(ZZ, 1, 'a') == ProjectiveSpace(ZZ, 0, 'a')
            False
            sage: ProjectiveSpace(ZZ, 2, 'a') == AffineSpace(ZZ, 2, 'a')
            False
            sage: loads(AffineSpace(ZZ, 1, 'x').dumps()) == AffineSpace(ZZ, 1, 'x')
            True
        """
        if not isinstance(right, ProjectiveSpace_ring):
            return -1
        return cmp([self.dimension_relative(), self.coordinate_ring()],
                   [right.dimension_relative(), right.coordinate_ring()])

    def _latex_(self):
        r"""
        Return a LaTeX representation of this projective space.

        EXAMPLES::

            sage: print(latex(ProjectiveSpace(1, ZZ, 'x')))
            {\mathbf P}_{\Bold{Z}}^1

        TESTS::

            sage: ProjectiveSpace(3, Zp(5), 'y')._latex_()
            '{\\mathbf P}_{\\ZZ_{5}}^3'
        """
        return "{\\mathbf P}_{%s}^%s"%(latex(self.base_ring()), self.dimension_relative())

    def _linear_system_as_kernel(self, d, pt, m):
        """
        Return a matrix whose kernel consists of the coefficient vectors
        of the degree ``d`` hypersurfaces (wrt lexicographic ordering of its
        monomials) with multiplicity at least ``m`` at ``pt``.

        INPUT:

        -  ``d`` -- a nonnegative integer.

        -  ``pt`` -- a point of self (possibly represented by a list with at \
                     least one component equal to 1).

        -  ``m`` -- a nonnegative integer.

        OUTPUT:

        A matrix of size `\binom{m-1+n}{n}` x `\binom{d+n}{n}` where n is the
        relative dimension of self. The base ring of the matrix is a ring that
        contains the base ring of self and the coefficients of the given point.

        EXAMPLES:

        If the degree `d` is 0, then a matrix consisting of the first unit vector
        is returned::

            sage: P = ProjectiveSpace(GF(5), 2, names='x')
            sage: pt = P([1, 1, 1])
            sage: P._linear_system_as_kernel(0, pt, 3)
            [1]
            [0]
            [0]
            [0]
            [0]
            [0]

        If the multiplcity ``m`` is 0, then the a matrix with zero rows is returned::

            sage: P = ProjectiveSpace(GF(5), 2, names='x')
            sage: pt = P([1, 1, 1])
            sage: M = P._linear_system_as_kernel(2, pt, 0)
            sage: [M.nrows(), M.ncols()]
            [0, 6]

        The base ring does not need to be a field or even an integral domain.
        In this case, the point can be given by a list::

            sage: R = Zmod(4)
            sage: P = ProjectiveSpace(R, 2, names='x')
            sage: pt = [R(1), R(3), R(0)]
            sage: P._linear_system_as_kernel(3, pt, 2)
            [1 3 0 1 0 0 3 0 0 0]
            [0 1 0 2 0 0 3 0 0 0]
            [0 0 1 0 3 0 0 1 0 0]

        When representing a point by a list at least one component must be 1
        (even when the base ring is a field and the list gives a well-defined
        point in projective space)::

            sage: R = GF(5)
            sage: P = ProjectiveSpace(R, 2, names='x')
            sage: pt = [R(3), R(3), R(0)]
            sage: P._linear_system_as_kernel(3, pt, 2)
            Traceback (most recent call last):
            ...
            TypeError: at least one component of pt=[3, 3, 0] must be equal
                          to 1

        The components of the list do not have to be elements of the base ring
        of the projective space. It suffices if there exists a common parent.
        For example, the kernel of the following matrix corresponds to
        hypersurfaces of degree 2 in 3-space with multiplicity at least 2 at a
        general point in the third affine patch::

            sage: P = ProjectiveSpace(QQ,3,names='x')
            sage: RPol.<t0,t1,t2,t3> = PolynomialRing(QQ,4)
            sage: pt = [t0,t1,1,t3]
            sage: P._linear_system_as_kernel(2,pt,2)
            [ 2*t0    t1     1    t3     0     0     0     0     0     0]
            [    0    t0     0     0  2*t1     1    t3     0     0     0]
            [ t0^2 t0*t1    t0 t0*t3  t1^2    t1 t1*t3     1    t3  t3^2]
            [    0     0     0    t0     0     0    t1     0     1  2*t3]

        .. TODO::

            Use this method as starting point to implement a class
            LinearSystem for linear systems of hypersurfaces.

        """
        if not isinstance(d, (int, Integer)):
            raise TypeError('the argument d=%s must be an integer'%d)
        if d < 0:
            raise ValueError('the integer d=%s must be nonnegative'%d)
        if not isinstance(pt, (list, tuple, \
                               SchemeMorphism_point_projective_ring)):
            raise TypeError('the argument pt=%s must be a list, tuple, or '
                            'point on a projective space'%pt)
        pt, R = prepare(pt, None)
        n = self.dimension_relative()
        if not len(pt) == n+1:
            raise TypeError('the sequence pt=%s must have %s '
                            'components'%(pt, n + 1))
        if not R.has_coerce_map_from(self.base_ring()):
            raise TypeError('unable to find a common ring for all elements')
        try:
            i = pt.index(1)
        except Exception:
            raise TypeError('at least one component of pt=%s must be equal '
                            'to 1'%pt)
        pt = pt[:i] + pt[i+1:]
        if not isinstance(m, (int, Integer)):
            raise TypeError('the argument m=%s must be an integer'%m)
        if m < 0:
            raise ValueError('the integer m=%s must be nonnegative'%m)
        # the components of partials correspond to partial derivatives
        # of order at most m-1 with respect to n variables
        partials = IntegerVectors(m-1, n+1).list()
        # the components of monoms correspond to monomials of degree
        # at most d in n variables
        monoms = IntegerVectors(d, n+1).list()
        M = matrix(R,len(partials),len(monoms))
        for row in range(M.nrows()):
            e = partials[row][:i] + partials[row][i+1:]
            for col in range(M.ncols()):
                f = monoms[col][:i] + monoms[col][i+1:]
                if min([f[j]-e[j] for j in range(n)]) >= 0:
                    M[row,col] = prod([ binomial(f[j],e[j]) * pt[j]**(f[j]-e[j])
                                        for j in (k for k in range(n) if f[k] > e[k]) ])
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

            sage: ProjectiveSpace(3, Zp(5), 'y')._repr_()
            'Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20'
        """
        return "Projective Space of dimension %s over %s"%(self.dimension_relative(), self.base_ring())

    def _repr_generic_point(self, v=None):
        """
        Return a string representation of the generic point
        corresponding to the list of polys ``v`` on this projective space.

        If ``v`` is None, the representation of the generic point of
        the projective space is returned.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._repr_generic_point([z*y-x^2])
            '(-x^2 + y*z)'
            sage: P._repr_generic_point()
            '(x : y : z)'
        """
        if v is None:
            v = self.gens()
        return '(%s)'%(" : ".join([repr(f) for f in v]))

    def _latex_generic_point(self, v=None):
        """
        Return a LaTeX representation of the generic point
        corresponding to the list of polys ``v`` on this projective space.

        If ``v`` is None, the representation of the generic point of
        the projective space is returned.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._latex_generic_point([z*y-x^2])
            '\\left(- x^{2} + y z\\right)'
            sage: P._latex_generic_point()
            '\\left(x : y : z\\right)'
        """
        if v is None:
            v = self.gens()
        return '\\left(%s\\right)'%(" : ".join([str(latex(f)) for f in v]))

    def change_ring(self, R):
        r"""
        Return a projective space over ring ``R``.

        INPUT:

        - ``R`` -- commutative ring or morphism.

        OUTPUT:

        - projective space over ``R``.

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

            sage: K.<w> = QuadraticField(2)
            sage: P = ProjectiveSpace(K,2,'t')
            sage: P.change_ring(K.embeddings(QQbar)[0])
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

        -  ``X`` - a list or tuple of equations.

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
            Ideal (x*z^2, y^2*z, x*y^2) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: I.groebner_basis()
            [x*y^2, y^2*z,  x*z^2]
            sage: X.dimension()
            0
            sage: X.base_ring()
            Rational Field
            sage: X.base_scheme()
            Spectrum of Rational Field
            sage: X.structure_morphism()
            Scheme morphism:
              From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*z^2,
              y^2*z,
              x*y^2
              To:   Spectrum of Rational Field
              Defn: Structure map

            sage: TestSuite(X).run(skip=["_test_an_element", "_test_elements",\
            "_test_elements_eq", "_test_some_elements", "_test_elements_eq_reflexive",\
            "_test_elements_eq_symmetric", "_test_elements_eq_transitive",\
            "_test_elements_neq"])
        """
        from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
        return AlgebraicScheme_subscheme_projective(self, X)

    def affine_patch(self, i, AA=None):
        r"""
        Return the `i^{th}` affine patch of this projective space.

        This is an ambient affine space `\mathbb{A}^n_R,` where
        `R` is the base ring of self, whose "projective embedding"
        map is `1` in the `i^{th}` factor.

        INPUT:

        - ``i`` -- integer between 0 and dimension of self, inclusive.

        - ``AA`` -- (default: None) ambient affine space, this is constructed
                if it is not given.

        OUTPUT:

        - An ambient affine space with fixed projective_embedding map.

        EXAMPLES::

            sage: PP = ProjectiveSpace(5) / QQ
            sage: AA = PP.affine_patch(2)
            sage: AA
            Affine Space of dimension 5 over Rational Field
            sage: AA.projective_embedding()
            Scheme morphism:
              From: Affine Space of dimension 5 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1, x2, x3, x4) to
                    (x0 : x1 : 1 : x2 : x3 : x4)
            sage: AA.projective_embedding(0)
            Scheme morphism:
              From: Affine Space of dimension 5 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1, x2, x3, x4) to
                    (1 : x0 : x1 : x2 : x3 : x4)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: P.affine_patch(0).projective_embedding(0).codomain() == P
            True
        """
        i = int(i)   # implicit type checking
        n = self.dimension_relative()
        if i < 0 or i > n:
            raise ValueError("argument i (= %s) must be between 0 and %s"%(i, n))
        try:
            A = self.__affine_patches[i]
            #assume that if you've passed in a new affine space you want to override
            #the existing patch
            if AA is None or A == AA:
                return(A)
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        #if no ith patch exists, we may still be here with AA==None
        if AA == None:
            from sage.schemes.affine.affine_space import AffineSpace
            AA = AffineSpace(n, self.base_ring(), names = 'x')
        elif AA.dimension_relative() != n:
                raise ValueError("affine space must be of the dimension %s"%(n))
        AA._default_embedding_index = i
        phi = AA.projective_embedding(i, self)
        self.__affine_patches[i] = AA
        return AA

    def _an_element_(self):
        r"""
        Returns a (preferably typical) element of this space.

        This is used both for illustration and testing purposes.

        OUTPUT: a point in this projective space.

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

        - ``E`` -- an elliptic curve.

        - ``m`` -- an integer.

        OUTPUT: an endomorphism of this projective space.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: E = EllipticCurve(QQ,[-1, 0])
            sage: P.Lattes_map(E, 2)
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 + 2*x^2*y^2 + y^4 : 4*x^3*y - 4*x*y^3)
        """
        if self.dimension_relative() != 1:
            raise TypeError("must be dimension 1")

        L = E.multiplication_by_m(m, x_only = True)
        F = [L.numerator(), L.denominator()]
        R = self.coordinate_ring()
        x, y = R.gens()
        phi = F[0].parent().hom([x],R)
        F = [phi(F[0]).homogenize(y), phi(F[1]).homogenize(y)*y]
        H = Hom(self,self)
        return(H(F))

    def cartesian_product(self, other):
        r"""
        Return the Cartesian product of this projective space and
        ``other``.

        INPUT:

        - ``other`` - A projective space with the same base ring as this space.

        OUTPUT:

        - A Cartesian product of projective spaces.

        EXAMPLES::

            sage: P1 = ProjectiveSpace(QQ, 1, 'x')
            sage: P2 = ProjectiveSpace(QQ, 2, 'y')
            sage: PP = P1.cartesian_product(P2); PP
            Product of projective spaces P^1 x P^2 over Rational Field
            sage: PP.gens()
            (x0, x1, y0, y1, y2)
        """
        from sage.schemes.product_projective.space import ProductProjectiveSpaces
        return ProductProjectiveSpaces([self, other])

    def chebyshev_polynomial(self, n, kind='first'):
        """
        Generates an endomorphism of this projective line by a Chebyshev polynomial.
    
        Chebyshev polynomials are a sequence of recursively defined orthogonal
        polynomials. Chebyshev of the first kind are defined as `T_0(x) = 1`,
        `T_1(x) = x`, and `T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)`. Chebyshev of
        the second kind are defined as `U_0(x) = 1`,
        `U_1(x) = 2x`, and `U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x)`.
    
        INPUT:
    
        - ``n`` -- a non-negative integer.
    
        - ``kind`` -- ``first`` or ``second`` specifying which kind of chebyshev the user would like
          to generate. Defaults to ``first``.
    
        OUTPUT: :class:`SchemeMorphism_polynomial_projective_space`
    
        EXAMPLES::
    
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(5, 'first')
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to
            (16*x^5 - 20*x^3*y^2 + 5*x*y^4 : y^5)
    
        ::
    
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(3, 'second')
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to
            (8*x^3 - 4*x*y^2 : y^3)
    
        ::
    
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(3, 2)
            Traceback (most recent call last):
            ...
            ValueError: keyword 'kind' must have a value of either 'first' or 'second'
    
        ::
    
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P.chebyshev_polynomial(-4, 'second')
            Traceback (most recent call last):
            ...
            ValueError: first parameter 'n' must be a non-negative integer

        ::

            sage: P = ProjectiveSpace(QQ, 2, 'x')
            sage: P.chebyshev_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: projective space must be of dimension 1
        """
        if self.dimension_relative() != 1:
            raise TypeError("projective space must be of dimension 1")
        n = ZZ(n)
        if (n < 0):
            raise ValueError("first parameter 'n' must be a non-negative integer")
        #use the affine version and then homogenize.        
        A = self.affine_patch(1)
        f = A.chebyshev_polynomial(n, kind)
        return f.homogenize(1)

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

    def points_of_bounded_height(self, bound, prec=53):
        r"""
        Returns an iterator of the points in self of absolute height of at most the given bound.

        Bound check is strict for the rational field. Requires self to be projective space
        over a number field. Uses the Doyle-Krumm algorithm for computing algebraic numbers
        up to a given height [Doyle-Krumm]_.

        INPUT:

        - ``bound`` - a real number.

        - ``prec`` - the precision to use to compute the elements of bounded height for number fields.

        OUTPUT:

        - an iterator of points in this space.

        .. WARNING::

           In the current implementation, the output of the [Doyle-Krumm]_ algorithm
           cannot be guaranteed to be correct due to the necessity of floating point
           computations. In some cases, the default 53-bit precision is
           considerably lower than would be required for the algorithm to
           generate correct output.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: list(P.points_of_bounded_height(5))
            [(0 : 1), (1 : 1), (-1 : 1), (1/2 : 1), (-1/2 : 1), (2 : 1), (-2 : 1), (1/3 : 1),
            (-1/3 : 1), (3 : 1), (-3 : 1), (2/3 : 1), (-2/3 : 1), (3/2 : 1), (-3/2 : 1), (1/4 : 1),
            (-1/4 : 1), (4 : 1), (-4 : 1), (3/4 : 1), (-3/4 : 1), (4/3 : 1), (-4/3 : 1), (1 : 0)]

        ::

            sage: u = QQ['u'].0
            sage: P.<x,y,z> = ProjectiveSpace(NumberField(u^2 - 2, 'v'), 2)
            sage: len(list(P.points_of_bounded_height(1.5)))
            57
        """
        if (is_RationalField(self.base_ring())):
            ftype = False # stores whether the field is a number field or the rational field
        elif (self.base_ring() in NumberFields()): # true for rational field as well, so check is_RationalField first
            ftype = True
        else:
            raise NotImplementedError("self must be projective space over a number field")

        bound = bound**(self.base_ring().absolute_degree()) # convert to relative height

        n = self.dimension_relative()
        R = self.base_ring()
        zero = R(0)
        i = n
        while not i < 0:
            P = [ zero for _ in range(i) ] + [ R(1) ] + [ zero for _ in range(n-i) ]
            yield self(P)
            if (ftype == False): # if rational field
                iters = [ R.range_by_height(bound) for _ in range(i) ]
            else: # if number field
                iters = [ R.elements_of_bounded_height(bound, precision=prec) for _ in range(i) ]
            for x in iters: next(x) # put at zero
            j = 0
            while j < i:
                try:
                    P[j] = next(iters[j])
                    yield self(P)
                    j = 0
                except StopIteration:
                    if (ftype == False): # if rational field
                        iters[j] = R.range_by_height(bound) # reset
                    else: # if number field
                        iters[j] = R.elements_of_bounded_height(bound, precision=prec) # reset
                    next(iters[j]) # put at zero
                    P[j] = zero
                    j += 1
            i -= 1

    def subscheme_from_Chow_form(self, Ch, dim):
        r"""
        Returns the subscheme defined by the Chow equations associated to the Chow form ``Ch``.

        These equations define the subscheme set-theoretically, but only for smooth
        subschemes and hypersurfaces do they define the subscheme as a scheme.

        ALGORITHM:

        The Chow form is a polynomial in the Plucker coordinates. The Plucker coordinates
        are the bracket polynomials. We first re-write the Chow form in terms of the dual
        Plucker coordinates. Then we expand `Ch(span(p,L)` for a generic point `p` and a
        generic linear subspace `L`. The coefficients as polynomials in the coordinates
        of `p` are the equations defining the subscheme. [DalbecSturmfels].

        INPUT:

        - ``Ch`` - a homogeneous polynomial.

        - ``dim`` - the dimension of the associated scheme.

        OUTPUT: a projective subscheme.

        EXAMPLES::

            sage: P = ProjectiveSpace(QQ, 4, 'z')
            sage: R.<x0,x1,x2,x3,x4> = PolynomialRing(QQ)
            sage: H = x1^2 + x2^2 + 5*x3*x4
            sage: P.subscheme_from_Chow_form(H,3)
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
              -5*z0*z1 + z2^2 + z3^2

        ::

            sage: P = ProjectiveSpace(QQ, 3, 'z')
            sage: R.<x0,x1,x2,x3,x4,x5> = PolynomialRing(QQ)
            sage: H = x1-x2-x3+x5+2*x0
            sage: P.subscheme_from_Chow_form(H, 1)
            Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
              -z1 + z3,
              z0 + z2 + z3,
              -z1 - 2*z3,
              -z0 - z1 + 2*z2

        ::

            sage: P.<x0,x1,x2,x3> = ProjectiveSpace(GF(7), 3)
            sage: X = P.subscheme([x3^2+x1*x2,x2-x0])
            sage: Ch = X.Chow_form();Ch
            t0^2 - 2*t0*t3 + t3^2 - t2*t4 - t4*t5
            sage: Y = P.subscheme_from_Chow_form(Ch, 1); Y
            Closed subscheme of Projective Space of dimension 3 over Finite Field of
            size 7 defined by:
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
            Ideal (x0 - x2, x1*x2 + x3^2) of Multivariate Polynomial Ring in x0, x1,
            x2, x3 over Finite Field of size 7
        """
        if not Ch.is_homogeneous():
            raise ValueError("Chow form must be a homogeneous polynomial")
        n = self.dimension_relative()
        R = Ch.parent()
        if binomial(n+1,n-dim) != R.ngens():
            raise ValueError("for given dimension, there should be %d variables in the Chow form" %binomial(n+1,n-dim))
        vars = list(R.gens())
        #create the brackets associated to variables
        L1 = []
        for t in UnorderedTuples(range(n+1), dim+1):
            if all([t[i]<t[i+1] for i in range(dim)]):
                L1.append(t)
        #create the dual brackets
        L2 = []
        signs = []
        for l in L1:
            s = []
            for v in range(n+1):
                if not v in l:
                    s.append(v)
            t1 = [b+1 for b in l]
            t2 = [b+1 for b in s]
            perm = Permutation(t1+t2)
            signs.append(perm.sign())
            L2.append(s)
        #create the polys associated to dual brackets
        if n-dim-1 > 0:
            S = PolynomialRing(R.base_ring(),n+1,'z')
            T = PolynomialRing(S,(n+1)*(n-dim-1),'s')
            M = matrix(T,n-dim,n+1,list(S.gens())+list(T.gens()))
        else:
            T = PolynomialRing(R.base_ring(),n+1,'z')
            M = matrix(T,n-dim,n+1,list(T.gens()))
        coords=[]
        for i in range(len(L2)):
            coords.append(signs[i]*M.matrix_from_columns(L2[i]).det())
        #substitute in dual brackets to chow form
        phi = R.hom(coords,T)
        ch = phi(Ch)
        #coefficients are polys in zs which are the chow equations for the chow form
        if n-dim-1 > 0:
            X = self.subscheme(ch.coefficients())
        else:
            X = self.subscheme(ch)
        return X

    def point_transformation_matrix(self, points_source, points_target):
        r"""

        Returns a unique element of PGL that transforms one set of points to another.

        Given a projective space of degree n and a set of n+2 source points and a set of n+2 target
        points in the same projective space, such that no n+1 points of each set are linearly dependent
        finds the unique element of PGL that translates the source points to the target points.


        Warning :: will not work over precision fields

        INPUT:

            - ``points_source`` - points in source projective space.

            - ``points_target`` - points in target projective space.

        OUTPUT: Transformation matrix - element of PGL.

        EXAMPLES::

            sage: P1.<a,b,c>=ProjectiveSpace(QQ, 2)
            sage: points_source=[P1([1,4,1]),P1([1,2,2]),P1([3,5,1]),P1([1,-1,1])]
            sage: points_target=[P1([5,-2,7]),P1([3,-2,3]),P1([6,-5,9]), P1([3,6,7])]
            sage: m = P1.point_transformation_matrix(points_source, points_target); m
            [ -13/59 -128/59  -25/59]
            [538/177    8/59  26/177]
            [ -45/59 -196/59       1]
            sage: [P1(list(m*vector(list(points_source[i])))) == points_target[i] for i in range(4)]
            [True, True, True, True]

        ::

            sage: P.<a,b> = ProjectiveSpace(GF(13),1)
            sage: points_source = [P([-6,7]), P([1,4]), P([3,2])]
            sage: points_target = [P([-1,2]), P([0,2]), P([-1,6])]
            sage: P.point_transformation_matrix(points_source, points_target)
            [10  4]
            [10  1]

        ::

            sage: P.<a,b> = ProjectiveSpace(QQ,1)
            sage: points_source = [P([-6,-4]), P([1,4]), P([3,2])]
            sage: points_target = [P([-1,2]), P([0,2]), P([-7,-3])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: source points not independent

        ::

            sage: P.<a,b> = ProjectiveSpace(QQ,1)
            sage: points_source = [P([-6,-1]), P([1,4]), P([3,2])]
            sage: points_target = [P([-1,2]), P([0,2]), P([-2,4])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: target points not independent

        ::

            sage: P.<a,b,c>=ProjectiveSpace(QQ, 2)
            sage: points_source=[P([1,4,1]),P([2,-7,9]),P([3,5,1])]
            sage: points_target=[P([5,-2,7]),P([3,-2,3]),P([6,-5,9]),P([6,-1,1])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: incorrect number of points in source, need 4 points

        ::

            sage: P.<a,b,c>=ProjectiveSpace(QQ, 2)
            sage: points_source=[P([1,4,1]),P([2,-7,9]),P([3,5,1]),P([1,-1,1])]
            sage: points_target=[P([5,-2,7]),P([3,-2,3]),P([6,-5,9]),P([6,-1,1]),P([7,8,-9])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: incorrect number of points in target, need 4 points

        ::

            sage: P.<a,b,c>=ProjectiveSpace(QQ, 2)
            sage: P1.<x,y,z>=ProjectiveSpace(QQ, 2)
            sage: points_source=[P([1,4,1]),P([2,-7,9]),P([3,5,1]),P1([1,-1,1])]
            sage: points_target=[P([5,-2,7]),P([3,-2,3]),P([6,-5,9]),P([6,-1,1])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: source points not in self

        ::

            sage: P.<a,b,c>=ProjectiveSpace(QQ, 2)
            sage: P1.<x,y,z>=ProjectiveSpace(QQ, 2)
            sage: points_source=[P([1,4,1]),P([2,-7,9]),P([3,5,1]),P([1,-1,1])]
            sage: points_target=[P([5,-2,7]),P([3,-2,3]),P([6,-5,9]),P1([6,-1,1])]
            sage: P.point_transformation_matrix(points_source, points_target)
            Traceback (most recent call last):
            ...
            ValueError: target points not in self
        """
        r = self.base_ring()
        n = self.dimension_relative()
        P = ProjectiveSpace(r, n**2+2*n,'p')
        # makes sure there aren't to few or two many points
        if len(points_source)!= n + 2:
            raise ValueError ("incorrect number of points in source, need %d points"%(n+2))
        if len(points_target)!= n + 2:
            raise ValueError ("incorrect number of points in target, need %d points"%(n+2))
        if any([x.codomain()!=self for x in points_source]):
            raise ValueError ("source points not in self")
        if any([x.codomain()!=self for x in points_target]):
            raise ValueError ("target points not in self")
        # putting points as the rows of the matrix
        Ms = matrix(r, [list(s) for s in points_source])
        if any([m == 0 for m in Ms.minors(n+1)]):
            raise ValueError("source points not independent")
        Mt = matrix(r, [list(t) for t in points_target])
        if any([l == 0 for l in Mt.minors(n+1)]):
            raise ValueError("target points not independent")
        A = matrix(P.coordinate_ring(), n+1, n+1, P.gens())
        #transpose to get image points and then get the list of image points with columns
        funct = (A*Ms.transpose()).columns()
        eq = []
        for k in range(n+2):# n+2 num f point and n is size of pts
            eq = eq+ [funct[k][i]*points_target[k][j] - funct[k][j]*points_target[k][i]\
                for i in range(0,n+1) for j in range(i+1, n+1)]
        v = P.subscheme(eq)
        w = v.rational_points()
        return matrix(r, n+1, n+1, list(w[0]))

    def curve(self,F):
        r"""
        Return a curve defined by ``F`` in this projective space.

        INPUT:

        - ``F`` -- a polynomial, or a list or tuple of polynomials in the coorinate ring
          of this projective space.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P.curve([y^2 - x*z])
            Projective Plane Curve over Rational Field defined by y^2 - x*z
        """
        from sage.schemes.curves.constructor import Curve
        return Curve(F, self)

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
        `\mathbb{P}^n = \mathbb{A}A^n \cup \mathbb{P}^n-1`, where
        `\mathbb{A}A^n` is the `n`-th affine patch and
        `\mathbb{P}^n-1` is the hyperplane at infinity
        `x_n = 0`.

        EXAMPLES::

            sage: FF = FiniteField(3)
            sage: PP = ProjectiveSpace(0,FF)
            sage: [ x for x in PP ]
            [(1)]
            sage: PP = ProjectiveSpace(1,FF)
            sage: [ x for x in PP ]
            [(0 : 1), (1 : 1), (2 : 1), (1 : 0)]
            sage: PP = ProjectiveSpace(2,FF)
            sage: [ x for x in PP ]
            [(0 : 0 : 1),
             (1 : 0 : 1),
             (2 : 0 : 1),
             (0 : 1 : 1),
             (1 : 1 : 1),
             (2 : 1 : 1),
             (0 : 2 : 1),
             (1 : 2 : 1),
             (2 : 2 : 1),
             (0 : 1 : 0),
             (1 : 1 : 0),
             (2 : 1 : 0),
             (1 : 0 : 0)]

        AUTHORS:

        - David Kohel

        .. TODO::

            Iteration for point sets over finite fields, and return of
            iter of point set over base field. Note that the point set does not
            know whether this is a projective space or subscheme.
        """
        n = self.dimension_relative()
        R = self.base_ring()
        zero = R(0)
        i = n
        while not i < 0:
            P = [ zero for _ in range(i) ] + [ R(1) ] + [ zero for _ in range(n-i) ]
            yield self(P)
            iters = [ iter(R) for _ in range(i) ]
            for x in iters: next(x) # put at zero
            j = 0
            while j < i:
                try:
                    P[j] = next(iters[j])
                    yield self(P)
                    j = 0
                except StopIteration:
                    iters[j] = iter(R)  # reset
                    next(iters[j]) # put at zero
                    P[j] = zero
                    j += 1
            i -= 1

    def rational_points(self, F=None):
        """
        Return the list of ``F``-rational points on this projective space,
        where ``F`` is a given finite field, or the base ring of this space.

        EXAMPLES::

            sage: P = ProjectiveSpace(1, GF(3))
            sage: P.rational_points()
            [(0 : 1), (1 : 1), (2 : 1), (1 : 0)]
            sage: P.rational_points(GF(3^2, 'b'))
            [(0 : 1), (b : 1), (b + 1 : 1), (2*b + 1 : 1), (2 : 1), (2*b : 1), (2*b + 2 : 1), (b + 2 : 1), (1 : 1), (1 : 0)]
        """
        if F is None:
            return [ P for P in self ]
        elif not is_FiniteField(F):
            raise TypeError("second argument (= %s) must be a finite field"%F)
        return [ P for P in self.base_extend(F) ]

    def rational_points_dictionary(self):
        r"""
        Return dictionary of points.

        OUTPUT:

        - dictionary

        EXAMPLES::

            sage: P1 = ProjectiveSpace(GF(7),1,'x')
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
        D={}
        zero = R(0)
        i = n
        index=0
        while not i < 0:
            P = [ zero for _ in range(i) ] + [ R(1) ] + [ zero for _ in range(n-i) ]
            D.update({self(P):index})
            index+=1
            iters = [ iter(R) for _ in range(i) ]
            for x in iters: next(x) # put at zero
            j = 0
            while j < i:
                try:
                    P[j] = next(iters[j])
                    D.update({self(P):index})
                    index+=1
                    j = 0
                except StopIteration:
                    iters[j] = iter(R)  # reset
                    next(iters[j]) # put at zero
                    P[j] = zero
                    j += 1
            i -= 1
        return(D)

class ProjectiveSpace_rational_field(ProjectiveSpace_field):
    def rational_points(self, bound=0):
        r"""
        Returns the projective points `(x_0:\cdots:x_n)` over
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

        -  ``bound`` - integer.

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
            (-2 : 0 : 1), (-1 : 0 : 1), (0 : 0 : 1), (1 : 0 : 1), (2 : 0 : 1), (-2 :
            1 : 1), (-1 : 1 : 1), (0 : 1 : 1), (1 : 1 : 1), (2 : 1 : 1), (-2 : 2 :
            1), (-1 : 2 : 1), (0 : 2 : 1), (1 : 2 : 1), (2 : 2 : 1), (-1/2 : -1 :
            1), (1/2 : -1 : 1), (-1 : -1/2 : 1), (-1/2 : -1/2 : 1), (0 : -1/2 : 1),
            (1/2 : -1/2 : 1), (1 : -1/2 : 1), (-1/2 : 0 : 1), (1/2 : 0 : 1), (-1 :
            1/2 : 1), (-1/2 : 1/2 : 1), (0 : 1/2 : 1), (1/2 : 1/2 : 1), (1 : 1/2 :
            1), (-1/2 : 1 : 1), (1/2 : 1 : 1), (-2 : 1 : 0), (-1 : 1 : 0), (0 : 1 :
            0), (1 : 1 : 0), (2 : 1 : 0), (-1/2 : 1 : 0), (1/2 : 1 : 0), (1 : 0 :
            0)]

        AUTHORS:

        - Benjamin Antieau (2008-01-12)
        """
        if not bound > 0:
            raise ValueError("argument bound (= %s) must be a positive integer")

        n = self.dimension_relative()

        Q = [k-bound for k in range(2*bound+1)]      # the affine coordinates
        R = [(k+1) for k in range(bound)]            # the projective coordinate
        S = [Tuples(Q, (k+1)) for k in range(n)]
        pts = []

        i = n
        while i > 0:
            P = [ 0 for _ in range(n+1) ]
            for ai in R:
                P[i] = ai
                for tup in S[i-1]:
                    if gcd([ai] + tup) == 1:
                        for j in range(i):
                            P[j] = tup[j]
                        pts.append(self(P))
            i -= 1

        # now do i=0; this is treated as a special case so that
        # we don't have all points (1:0),(2,0),(3,0),etc.
        P = [ 0 for _ in range(n+1) ]; P[0] = 1
        pts.append(self(P))
        return pts


#fix the pickles from moving projective_space.py
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.schemes.generic.projective_space',
                           'ProjectiveSpace_field',
                           ProjectiveSpace_field)

register_unpickle_override('sage.schemes.generic.projective_space',
                           'ProjectiveSpace_rational_field',
                           ProjectiveSpace_rational_field)

