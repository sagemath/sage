r"""
Algebraic schemes

An algebraic scheme is defined by a set of polynomials in some
suitable affine or projective coordinates. Possible ambient spaces are

  * Affine spaces (:class:`AffineSpace
    <sage.schemes.affine.affine_space.AffineSpace_generic>`),

  * Projective spaces (:class:`ProjectiveSpace
    <sage.schemes.projective.projective_space.ProjectiveSpace_ring>`), or

  * Toric varieties (:class:`ToricVariety
    <sage.schemes.toric.variety.ToricVariety_field>`).

Note that while projective spaces are of course toric varieties themselves,
they are implemented differently in Sage due to efficiency considerations.
You still can create a projective space as a toric variety if you wish.

In the following, we call the corresponding subschemes affine
algebraic schemes, projective algebraic schemes, or toric algebraic
schemes. In the future other ambient spaces, perhaps by means of
gluing relations, may be intoduced.

Generally, polynomials `p_0, p_1, \dots, p_n` define an ideal
`I=\left<p_0, p_1, \dots, p_n\right>`. In the projective and toric case, the
polynomials (and, therefore, the ideal) must be homogeneous. The
associated subscheme `V(I)` of the ambient space is, roughly speaking,
the subset of the ambient space on which all polynomials vanish simultaneously.

.. WARNING::

    You should not construct algebraic scheme objects directly. Instead, use
    ``.subscheme()`` methods of ambient spaces. See below for examples.

EXAMPLES:

We first construct the ambient space, here the affine space `\QQ^2`::

    sage: A2 = AffineSpace(2, QQ, 'x, y')
    sage: A2.coordinate_ring().inject_variables()
    Defining x, y

Now we can write polynomial equations in the variables `x` and `y`. For
example, one equation cuts out a curve (a one-dimensional subscheme)::

    sage: V = A2.subscheme([x^2+y^2-1]); V
    Closed subscheme of Affine Space of dimension 2
    over Rational Field defined by:
      x^2 + y^2 - 1
    sage: V.dimension()
    1

Here is a more complicated example in a projective space::

    sage: P3 = ProjectiveSpace(3, QQ, 'x')
    sage: P3.inject_variables()
    Defining x0, x1, x2, x3
    sage: Q = matrix([[x0, x1, x2], [x1, x2, x3]]).minors(2); Q
    [-x1^2 + x0*x2, -x1*x2 + x0*x3, -x2^2 + x1*x3]
    sage: twisted_cubic = P3.subscheme(Q)
    sage: twisted_cubic
    Closed subscheme of Projective Space of dimension 3
    over Rational Field defined by:
      -x1^2 + x0*x2,
      -x1*x2 + x0*x3,
      -x2^2 + x1*x3
    sage: twisted_cubic.dimension()
    1

Note that there are 3 equations in the 3-dimensional ambient space,
yet the subscheme is 1-dimensional. One can show that it is not
possible to eliminate any of the equations, that is, the twisted cubic
is **not** a complete intersection of two polynomial equations.

Let us look at one affine patch, for example the one where `x_0=1` ::

    sage: patch = twisted_cubic.affine_patch(0)
    sage: patch
    Closed subscheme of Affine Space of dimension 3
    over Rational Field defined by:
      -x0^2 + x1,
      -x0*x1 + x2,
      -x1^2 + x0*x2
    sage: patch.embedding_morphism()
    Scheme morphism:
      From: Closed subscheme of Affine Space of dimension 3
      over Rational Field defined by:
      -x0^2 + x1,
      -x0*x1 + x2,
      -x1^2 + x0*x2
      To:   Closed subscheme of Projective Space of dimension 3
      over Rational Field defined by:
      x1^2 - x0*x2,
      x1*x2 - x0*x3,
      x2^2 - x1*x3
      Defn: Defined on coordinates by sending (x0, x1, x2) to
            (1 : x0 : x1 : x2)


AUTHORS:

- David Kohel (2005): initial version.
- William Stein (2005): initial version.
- Andrey Novoseltsev (2010-05-17): subschemes of toric varieties.
- Volker Braun (2010-12-24): documentation of schemes and
  refactoring. Added coordinate neighborhoods and is_smooth()
- Ben Hutz (2014): subschemes of Cartesian products of projective space
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


#*** A quick overview over the class hierarchy:
# class AlgebraicScheme(scheme.Scheme)
#    class AlgebraicScheme_subscheme
#       class AlgebraicScheme_subscheme_affine
#       class AlgebraicScheme_subscheme_projective
#       class AlgebraicScheme_subscheme_toric
#          class AlgebraicScheme_subscheme_affine_toric
#    class AlgebraicScheme_quasi

from sage.arith.misc import binomial
from sage.matrix.constructor import matrix
from sage.combinat.tuple import UnorderedTuples

from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.categories.morphism import Morphism

from sage.interfaces.all import singular

from sage.rings.all import ZZ
from sage.rings.ideal import is_Ideal
from sage.rings.rational_field import is_RationalField
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.fraction_field import FractionField

from sage.misc.all import prod
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.misc.latex import latex
from sage.misc.misc import is_iterator
from sage.structure.all import Sequence
from sage.calculus.functions import jacobian

import sage.schemes.affine
from . import ambient_space
from . import scheme



#*******************************************************************
def is_AlgebraicScheme(x):
    """
    Test whether ``x`` is an algebraic scheme.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean. Whether ``x`` is an an algebraic scheme, that is, a
    subscheme of an ambient space over a ring defined by polynomial
    equations.

    EXAMPLES::

        sage: A2 = AffineSpace(2, QQ, 'x, y')
        sage: A2.coordinate_ring().inject_variables()
        Defining x, y
        sage: V = A2.subscheme([x^2+y^2]); V
        Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
          x^2 + y^2
        sage: from sage.schemes.generic.algebraic_scheme import is_AlgebraicScheme
        sage: is_AlgebraicScheme(V)
        True

    Affine space is itself not an algebraic scheme, though the closed
    subscheme defined by no equations is::

        sage: from sage.schemes.generic.algebraic_scheme import is_AlgebraicScheme
        sage: is_AlgebraicScheme(AffineSpace(10, QQ))
        False
        sage: V = AffineSpace(10, QQ).subscheme([]); V
        Closed subscheme of Affine Space of dimension 10 over Rational Field defined by:
          (no polynomials)
        sage: is_AlgebraicScheme(V)
        True

    We create a more complicated closed subscheme::

        sage: A,x = AffineSpace(10, QQ).objgens()
        sage: X = A.subscheme([sum(x)]); X
        Closed subscheme of Affine Space of dimension 10 over Rational Field defined by:
        x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
        sage: is_AlgebraicScheme(X)
        True

    ::

        sage: is_AlgebraicScheme(QQ)
        False
        sage: S = Spec(QQ)
        sage: is_AlgebraicScheme(S)
        False
    """
    return isinstance(x, AlgebraicScheme)



#*******************************************************************
class AlgebraicScheme(scheme.Scheme):
    """
    An algebraic scheme presented as a subscheme in an ambient space.

    This is the base class for all algebraic schemes, that is, schemes
    defined by equations in affine, projective, or toric ambient
    spaces.
    """

    def __init__(self, A):
        """
        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme
            sage: P = ProjectiveSpace(3, ZZ)
            sage: P.category()
            Category of schemes over Integer Ring
            sage: S = AlgebraicScheme(P); S
            Subscheme of Projective Space of dimension 3 over Integer Ring
            sage: S.category()
            Category of schemes over Integer Ring
        """
        if not ambient_space.is_AmbientSpace(A):
            raise TypeError("A (=%s) must be an ambient space")
        self.__A = A
        self.__divisor_group = {}
        scheme.Scheme.__init__(self, A.base_scheme())

    def _latex_(self):
        """
        Return a LaTeX representation of this algebraic scheme.

        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme
            sage: P = ProjectiveSpace(3, ZZ)
            sage: S = AlgebraicScheme(P); S
            Subscheme of Projective Space of dimension 3 over Integer Ring
            sage: S._latex_()
            '\text{Subscheme of } {\\mathbf P}_{\\Bold{Z}}^3'
        """
        return "\text{Subscheme of } %s" % latex(self.__A)

    def is_projective(self):
        """
        Return True if self is presented as a subscheme of an ambient
        projective space.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: PP.<x,y,z,w> = ProjectiveSpace(3,QQ)
            sage: f = x^3 + y^3 + z^3 + w^3
            sage: R = f.parent()
            sage: I = [f] + [f.derivative(zz) for zz in PP.gens()]
            sage: V = PP.subscheme(I)
            sage: V.is_projective()
            True
            sage: AA.<x,y,z,w> = AffineSpace(4,QQ)
            sage: V = AA.subscheme(I)
            sage: V.is_projective()
            False

        Note that toric varieties are implemented differently than
        projective spaces. This is why this method returns ``False``
        for toric varieties::

            sage: PP.<x,y,z,w> = toric_varieties.P(3)
            sage: V = PP.subscheme(x^3 + y^3 + z^3 + w^3)
            sage: V.is_projective()
            False
        """
        return self.ambient_space().is_projective()

    def coordinate_ring(self):
        """
        Return the coordinate ring of this algebraic scheme.  The
        result is cached.

        OUTPUT:

        The coordinate ring. Usually a polynomial ring, or a quotient
        thereof.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x-y, x-z])
            sage: S.coordinate_ring()
            Quotient of Multivariate Polynomial Ring in x, y, z over Integer Ring by the ideal (x - y, x - z)
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            R = self.__A.coordinate_ring()
            I = self.defining_ideal()
            Q = R.quotient(I)
            self._coordinate_ring = Q
            return Q

    def ambient_space(self):
        """
        Return the ambient space of this algebraic scheme.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, GF(5))
            sage: S = A.subscheme([])
            sage: S.ambient_space()
            Affine Space of dimension 2 over Finite Field of size 5

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x-y, x-z])
            sage: S.ambient_space() is P
            True
        """
        return self.__A

    def embedding_morphism(self):
        r"""
        Return the default embedding morphism of ``self``.

        If the scheme `Y` was constructed as a neighbourhood of a
        point `p \in X`, then :meth:`embedding_morphism` returns a
        local isomorphism `f:Y\to X` around the preimage point
        `f^{-1}(p)`. The latter is returned by
        :meth:`embedding_center`.

        If the algebraic scheme `Y` was not constructed as a
        neighbourhood of a point, then the embedding in its
        :meth:`ambient_space` is returned.

        OUTPUT:

        A scheme morphism whose
        :meth:`~morphism.SchemeMorphism.domain` is ``self``.

        * By default, it is the tautological embedding into its own
          ambient space :meth:`ambient_space`.

        * If the algebraic scheme (which itself is a subscheme of an
          auxiliary :meth:`ambient_space`) was constructed as a patch
          or neighborhood of a point then the embedding is the
          embedding into the original scheme.

        * A ``NotImplementedError`` is raised if the construction of
          the embedding morphism is not implemented yet.

        EXAMPLES::

            sage: A2.<x,y> = AffineSpace(QQ,2)
            sage: C = A2.subscheme(x^2+y^2-1)
            sage: C.embedding_morphism()
              Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^2 + y^2 - 1
              To:   Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (x, y)
            sage: P1xP1.<x,y,u,v> = toric_varieties.P1xP1()
            sage: P1 = P1xP1.subscheme(x-y)
            sage: P1.embedding_morphism()
            Scheme morphism:
            From: Closed subscheme of 2-d CPR-Fano toric variety covered
                  by 4 affine patches defined by:
            x - y
            To:   2-d CPR-Fano toric variety covered by 4 affine patches
            Defn: Defined on coordinates by sending [x : y : u : v] to
                  [y : y : u : v]

        So far, the embedding was just in the own ambient space. Now a
        bit more interesting examples::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P2.subscheme((x^2-y^2)*z)
            sage: p = (1,1,0)
            sage: nbhd = X.neighborhood(p)
            sage: nbhd
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -x0^2*x1 - 2*x0*x1

        Note that `p=(1,1,0)` is a singular point of `X`. So the
        neighborhood of `p` is not just affine space. The
        :meth:neighborhood` method returns a presentation of
        the neighborhood as a subscheme of an auxiliary 2-dimensional
        affine space::

            sage: nbhd.ambient_space()
            Affine Space of dimension 2 over Rational Field

        But its :meth:`embedding_morphism` is not into this auxiliary
        affine space, but the original subscheme `X`::

            sage: nbhd.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -x0^2*x1 - 2*x0*x1
              To:   Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2*z - y^2*z
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1 : x0 + 1 : x1)

        A couple more examples::

            sage: patch1 = P1xP1.affine_patch(1)
            sage: patch1
            2-d affine toric variety
            sage: patch1.embedding_morphism()
              Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [y : u] to
                    [1 : y : u : 1]
            sage: subpatch = P1.affine_patch(1)
            sage: subpatch
            Closed subscheme of 2-d affine toric variety defined by:
              -y + 1
            sage: subpatch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of 2-d affine toric variety defined by:
              -y + 1
              To:   Closed subscheme of 2-d CPR-Fano toric variety covered
                    by 4 affine patches defined by:
              x - y
              Defn: Defined on coordinates by sending [y : u] to
                    [1 : y : u : 1]
        """
        if '_embedding_morphism' in self.__dict__:
            hom = self._embedding_morphism
            if isinstance(hom, tuple):
                raise hom[0]
            return hom
        ambient = self.ambient_space()
        return self.hom(self.coordinate_ring().gens(), ambient)

    def embedding_center(self):
        r"""
        Return the distinguished point, if there is any.

        If the scheme `Y` was constructed as a neighbourhood of a
        point `p \in X`, then :meth:`embedding_morphism` returns a
        local isomorphism `f:Y\to X` around the preimage point
        `f^{-1}(p)`. The latter is returned by
        :meth:`embedding_center`.

        OUTPUT:

        A point of ``self``. Raises ``AttributeError`` if there is no
        distinguished point, depending on how ``self`` was
        constructed.

        EXAMPLES::

            sage: P3.<w,x,y,z> = ProjectiveSpace(QQ,3)
            sage: X = P3.subscheme( (w^2-x^2)*(y^2-z^2) )
            sage: p = [1,-1,3,4]
            sage: nbhd = X.neighborhood(p); nbhd
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              x0^2*x2^2 - x1^2*x2^2 + 6*x0^2*x2 - 6*x1^2*x2 + 2*x0*x2^2 +
              2*x1*x2^2 - 7*x0^2 + 7*x1^2 + 12*x0*x2 + 12*x1*x2 - 14*x0 - 14*x1
            sage: nbhd.embedding_center()
            (0, 0, 0)
            sage: nbhd.embedding_morphism()(nbhd.embedding_center())
            (1/4 : -1/4 : 3/4 : 1)
            sage: nbhd.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              x0^2*x2^2 - x1^2*x2^2 + 6*x0^2*x2 - 6*x1^2*x2 + 2*x0*x2^2 +
              2*x1*x2^2 - 7*x0^2 + 7*x1^2 + 12*x0*x2 + 12*x1*x2 - 14*x0 - 14*x1
              To:   Closed subscheme of Projective Space of dimension 3 over Rational Field defined by:
              w^2*y^2 - x^2*y^2 - w^2*z^2 + x^2*z^2
              Defn: Defined on coordinates by sending (x0, x1, x2) to
                    (x0 + 1 : x1 - 1 : x2 + 3 : 4)
        """
        if '_embedding_center' in self.__dict__:
            return self._embedding_center
        raise AttributeError('This algebraic scheme does not have a designated point.')

    def ngens(self):
        """
        Return the number of generators of the ambient space of this
        algebraic scheme.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, GF(5))
            sage: S = A.subscheme([])
            sage: S.ngens()
            2
            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x-y, x-z])
            sage: P.ngens()
            3
        """
        return self.__A.ngens()

    def _repr_(self):
        """
        Return a string representation of this algebraic scheme.

        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme
            sage: P = ProjectiveSpace(3, ZZ)
            sage: S = AlgebraicScheme(P); S
            Subscheme of Projective Space of dimension 3 over Integer Ring
            sage: S._repr_()
            'Subscheme of Projective Space of dimension 3 over Integer Ring'
        """
        return "Subscheme of %s"%self.__A

    def _homset(self, *args, **kwds):
        """
        Construct the Hom-set

        INPUT:

        Same as :class:`sage.schemes.generic.homset.SchemeHomset_generic`.

        OUTPUT:

        The Hom-set of the ambient space.

        EXAMPLES::

            sage: P1.<x,y> = toric_varieties.P1()
            sage: type(P1.Hom(P1))
            <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>
            sage: X = P1.subscheme(x-y)
            sage: type(X.Hom(X))
            <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>

        ::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: P1xP1._homset(P1xP1,P1)
            Set of morphisms
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   1-d CPR-Fano toric variety covered by 2 affine patches
        """
        return self.__A._homset(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set. For internal use only.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, ZZ)
            sage: P2._point_homset(Spec(ZZ), P2)
            Set of rational points of Projective Space of dimension 2 over Integer Ring
        """
        return self.__A._point_homset(*args, **kwds)

    def _point(self, *args, **kwds):
        r"""
        Construct a point of ``self``. For internal use only.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, QQ)
            sage: point_homset = P2._point_homset(Spec(QQ), P2)
            sage: P2._point(point_homset, [1,2,1])
            (1 : 2 : 1)
        """
        return self.__A._point(*args, **kwds)



#*******************************************************************
class AlgebraicScheme_quasi(AlgebraicScheme):
    """
    The quasi-affine or quasi-projective scheme `X - Y`, where `X` and `Y`
    are both closed subschemes of a common ambient affine or projective
    space.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`complement` method of algebraic schemes.

    OUTPUT:

    An instance of :class:`AlgebraicScheme_quasi`.

    EXAMPLES::

        sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
        sage: S = P.subscheme([])
        sage: T = P.subscheme([x-y])
        sage: T.complement(S)
        Quasi-projective subscheme X - Y of Projective Space of dimension 2 over
        Integer Ring, where X is defined by:
          (no polynomials)
        and Y is defined by:
          x - y
    """

    def __init__(self, X, Y):
        """
        The constructor.

        INPUT:

        - ``X``, ``Y`` -- two subschemes of the same ambient space.

        TESTS::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_quasi
            sage: AlgebraicScheme_quasi(S, T)
            Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring, where X is defined by:
              (no polynomials)
            and Y is defined by:
              x - y
        """
        self.__X = X
        self.__Y = Y
        if not isinstance(X, AlgebraicScheme_subscheme):
            raise TypeError("X must be a closed subscheme of an ambient space.")
        if not isinstance(Y, AlgebraicScheme_subscheme):
            raise TypeError("Y must be a closed subscheme of an ambient space.")
        if X.ambient_space() != Y.ambient_space():
            raise ValueError("X and Y must be embedded in the same ambient space.")
        # _latex_ and _repr_ assume all of the above conditions and should be
        # probably changed if they are relaxed!
        A = X.ambient_space()
        self._base_ring = A.base_ring()
        AlgebraicScheme.__init__(self, A)

    def _latex_(self):
        """
        Return a LaTeX representation of this algebraic scheme.

        EXAMPLES::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_quasi
            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = AlgebraicScheme_quasi(S, T); U
            Quasi-projective subscheme X - Y of Projective Space of dimension 2
            over Integer Ring, where X is defined by:
              (no polynomials)
            and Y is defined by:
              x - y
            sage: U._latex_()
            '\\text{Quasi-projective subscheme }
             (X\\setminus Y)\\subset {\\mathbf P}_{\\Bold{Z}}^2,\\text{ where }
             X \\text{ is defined by }\\text{no polynomials},\\text{ and }
             Y \\text{ is defined by } x - y.'
        """
        if sage.schemes.affine.affine_space.is_AffineSpace(self.ambient_space()):
            t = "affine"
        else:
            t = "projective"
        X = ', '.join(latex(f) for f in self.__X.defining_polynomials())
        if not X:
            X = r"\text{no polynomials}"
        Y = ', '.join(latex(f) for f in self.__Y.defining_polynomials())
        if not Y:
            Y = r"\text{no polynomials}"
        return (r"\text{Quasi-%s subscheme } (X\setminus Y)\subset %s,"
                r"\text{ where } X \text{ is defined by }%s,"
                r"\text{ and } Y \text{ is defined by } %s."
                % (t, latex(self.ambient_space()), X, Y))

    def _repr_(self):
        r"""
        Return a string representation of this algebraic scheme.

        EXAMPLES::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_quasi
            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = AlgebraicScheme_quasi(S, T); U
            Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring, where X is defined by:
              (no polynomials)
            and Y is defined by:
              x - y
            sage: U._repr_()
            'Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring, where X is defined by:\n  (no polynomials)\nand Y is defined by:\n  x - y'
        """
        if sage.schemes.affine.affine_space.is_AffineSpace(self.ambient_space()):
            t = "affine"
        else:
            t = "projective"
        return ("Quasi-%s subscheme X - Y of %s, where X is defined by:\n%s\n"
                "and Y is defined by:\n%s"
                % (t, self.ambient_space(), str(self.__X).split("\n", 1)[1],
                   str(self.__Y).split("\n", 1)[1]))

    def X(self):
        """
        Return the scheme `X` such that self is represented as `X - Y`.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U.X() is S
            True
        """
        return self.__X

    def Y(self):
        """
        Return the scheme `Y` such that self is represented as `X - Y`.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U.Y() is T
            True
        """
        return self.__Y

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme, or
        raise a TypeError.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U._check_satisfies_equations([1, 2, 0])
            True
            sage: U._check_satisfies_equations([1, 1, 0])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 1, 0] do not define a point on
            Quasi-projective subscheme X - Y of Projective Space of dimension 2
            over Integer Ring, where X is defined by:
              (no polynomials)
            and Y is defined by:
              x - y

            sage: U._check_satisfies_equations([1, 4])
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match number of variables in parent

            sage: A.<x, y> = AffineSpace(2, GF(7))
            sage: S = A.subscheme([x^2-y])
            sage: T = A.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U._check_satisfies_equations([2, 4])
            True
            sage: U.point([2,4])
            (2, 4)
            sage: U._check_satisfies_equations(_)
            True
            sage: U._check_satisfies_equations([1, 1])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 1] do not define a point on Quasi-affine
            subscheme X - Y of Affine Space of dimension 2 over Finite
            Field of size 7, where X is defined by:
              x^2 - y
            and Y is defined by:
              x - y
            sage: U._check_satisfies_equations([1, 0])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 0] do not define a point on Quasi-affine
            subscheme X - Y of Affine Space of dimension 2 over Finite
            Field of size 7, where X is defined by:
              x^2 - y
            and Y is defined by:
              x - y

        TESTS:

        The bug reported at #12211 has been fixed::

            sage: P.<x, y, z, w> = ProjectiveSpace(3, QQ)
            sage: S = P.subscheme([x])
            sage: T = P.subscheme([y, z])
            sage: U = T.complement(S)
            sage: U._check_satisfies_equations([0, 0, 1, 1])
            True
        """
        coords = list(v)
        for f in self.__X.defining_polynomials():
            if f(coords) != 0:
                raise TypeError("Coordinates %s do not define a point on %s"%(v,self))
        for f in self.__Y.defining_polynomials():
            if f(coords) != 0:
                return True
        raise TypeError("Coordinates %s do not define a point on %s"%(v,self))

    def rational_points(self, F=None, bound=0):
        """
        Return the set of rational points on this algebraic scheme
        over the field `F`.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, GF(7))
            sage: S = A.subscheme([x^2-y])
            sage: T = A.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U.rational_points()
            [(2, 4), (3, 2), (4, 2), (5, 4), (6, 1)]
            sage: U.rational_points(GF(7^2, 'b'))
            [(2, 4), (3, 2), (4, 2), (5, 4), (6, 1), (b, b + 4), (b + 1, 3*b + 5), (b + 2, 5*b + 1),
            (b + 3, 6), (b + 4, 2*b + 6), (b + 5, 4*b + 1), (b + 6, 6*b + 5), (2*b, 4*b + 2),
            (2*b + 1, b + 3), (2*b + 2, 5*b + 6), (2*b + 3, 2*b + 4), (2*b + 4, 6*b + 4),
            (2*b + 5, 3*b + 6), (2*b + 6, 3), (3*b, 2*b + 1), (3*b + 1, b + 2), (3*b + 2, 5),
            (3*b + 3, 6*b + 3), (3*b + 4, 5*b + 3), (3*b + 5, 4*b + 5), (3*b + 6, 3*b + 2),
            (4*b, 2*b + 1), (4*b + 1, 3*b + 2), (4*b + 2, 4*b + 5), (4*b + 3, 5*b + 3),
            (4*b + 4, 6*b + 3), (4*b + 5, 5), (4*b + 6, b + 2), (5*b, 4*b + 2), (5*b + 1, 3),
            (5*b + 2, 3*b + 6), (5*b + 3, 6*b + 4), (5*b + 4, 2*b + 4), (5*b + 5, 5*b + 6),
            (5*b + 6, b + 3), (6*b, b + 4), (6*b + 1, 6*b + 5), (6*b + 2, 4*b + 1), (6*b + 3, 2*b + 6),
            (6*b + 4, 6), (6*b + 5, 5*b + 1), (6*b + 6, 3*b + 5)]
        """
        if F is None:
            F = self.base_ring()

        if bound == 0:
            if is_RationalField(F):
                raise TypeError("A positive bound (= %s) must be specified."%bound)
            if not is_FiniteField(F):
                raise TypeError("Argument F (= %s) must be a finite field."%F)
        pts = []
        for P in self.ambient_space().rational_points(F):
            try:
                if self._check_satisfies_equations(list(P)):
                    pts.append(P)
            except TypeError:
                pass
        pts.sort()
        return pts



#*******************************************************************
class AlgebraicScheme_subscheme(AlgebraicScheme):
    """
    An algebraic scheme presented as a closed subscheme is defined by
    explicit polynomial equations. This is as opposed to a general
    scheme, which could, e.g., be the Neron model of some object, and
    for which we do not want to give explicit equations.

    INPUT:

    -  ``A`` - ambient space (e.g. affine or projective `n`-space)

    -  ``polynomials`` - single polynomial, ideal or iterable of defining
        polynomials; in any case polynomials must belong to the coordinate
        ring of the ambient space and define valid polynomial functions (e.g.
        they should be homogeneous in the case of a projective space)

    OUTPUT:

    - algebraic scheme

    EXAMPLES::

        sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
        sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
        sage: P.subscheme([x^2-y*z])
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x^2 - y*z
        sage: AlgebraicScheme_subscheme(P, [x^2-y*z])
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x^2 - y*z
    """

    def __init__(self, A, polynomials):
        """
        See ``AlgebraicScheme_subscheme`` for documentation.

        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: P.subscheme([x^2-y*z])
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
            sage: AlgebraicScheme_subscheme(P, [x^2-y*z])
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
        """
        from sage.rings.polynomial.multi_polynomial_sequence import is_PolynomialSequence

        AlgebraicScheme.__init__(self, A)
        self._base_ring = A.base_ring()
        R = A.coordinate_ring()
        if is_Ideal(polynomials):
            I = polynomials
            polynomials = I.gens()
            if I.ring() is R: # Otherwise we will recompute I later after
                self.__I = I  # converting generators to the correct ring
        if isinstance(polynomials, tuple) or is_PolynomialSequence(polynomials) or is_iterator(polynomials):
            polynomials = list(polynomials)
        elif not isinstance(polynomials, list):
            # Looks like we got a single polynomial
            polynomials = [polynomials]
        for n, f in enumerate(polynomials):
            try:
                polynomials[n] = R(f)
            except TypeError:
                raise TypeError("%s cannot be converted to a polynomial in "
                                "the coordinate ring of this %s!" % (f, A))
        polynomials = tuple(polynomials)
        self.__polys = A._validate(polynomials)

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme, or
        raise a TypeError.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: S = P.subscheme([x^2-y*z])
            sage: S._check_satisfies_equations([1, 1, 1])
            True
            sage: S._check_satisfies_equations([1, 0, 1])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 0, 1] do not define a point on Closed subscheme
            of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
            sage: S._check_satisfies_equations([0, 0, 0])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [0, 0, 0] do not define a point on Closed subscheme
            of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
        """
        coords = list(v)
        for f in self.defining_polynomials():
            if f(coords) != 0:   # it must be "!=0" instead of "if f(v)", e.g.,
                                 # because of p-adic base rings.
                raise TypeError("Coordinates %s do not define a point on %s"%(coords,self))
        try:
            return self.ambient_space()._check_satisfies_equations(coords)
        except TypeError:
            raise TypeError("Coordinates %s do not define a point on %s"%(coords,self))

    def base_extend(self, R):
        """
        Return the base change to the ring `R` of this scheme.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(11))
            sage: S = P.subscheme([x^2-y*z])
            sage: S.base_extend(GF(11^2, 'b'))
            Closed subscheme of Projective Space of dimension 2 over Finite Field in b of size 11^2 defined by:
              x^2 - y*z
            sage: S.base_extend(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: no natural map from the base ring (=Finite Field of size 11) to R (=Integer Ring)!
        """
        A = self.ambient_space().base_extend(R)
        return A.subscheme(self.__polys)

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, QQ)
            sage: X = A.subscheme([x*y, z])
            sage: X == A.subscheme([z, x*y])
            True
            sage: X == A.subscheme([x*y, z^2])
            False
            sage: B.<u, v, t> = AffineSpace(3, QQ)
            sage: X == B.subscheme([u*v, t])
            False
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            return -1
        A = self.ambient_space()
        if other.ambient_space() != A:
            return -1
        return cmp(self.defining_ideal(), other.defining_ideal())

    def _latex_(self):
        """
        Return a LaTeX representation of this scheme.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(11))
            sage: S = P.subscheme([x^2-y*z])
            sage: S
            Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:
              x^2 - y*z
            sage: S._latex_()
            '\\text{Closed subscheme of } {\\mathbf P}_{\\Bold{F}_{11}}^2 \\text{ defined by } x^{2} - y z'
            sage: S = P.subscheme([x^2-y*z, x^5])
            sage: S
            Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:
              x^2 - y*z,
              x^5
            sage: S._latex_()
            '\\text{Closed subscheme of } {\\mathbf P}_{\\Bold{F}_{11}}^2 \\text{ defined by } x^{2} - y z, x^{5}'
        """
        polynomials = ', '.join(latex(f) for f in self.defining_polynomials())
        if not polynomials:
            polynomials = r"\text{no polynomials}"
        return (r"\text{Closed subscheme of } %s \text{ defined by } %s"
                % (latex(self.ambient_space()), polynomials))

    def _repr_(self):
        r"""
        Return a string representation of this scheme.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(11))
            sage: S = P.subscheme([x^2-y*z])
            sage: S
            Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:
              x^2 - y*z
            sage: S._repr_()
            'Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:\n  x^2 - y*z'
            sage: S = P.subscheme([x^2-y*z, x^5])
            sage: S
            Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:
              x^2 - y*z,
              x^5
            sage: S._repr_()
            'Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:\n  x^2 - y*z,\n  x^5'
        """
        polynomials = ',\n  '.join(str(f) for f in self.defining_polynomials())
        if not polynomials:
            polynomials = '(no polynomials)'
        return ("Closed subscheme of %s defined by:\n  %s"
                % (self.ambient_space(), polynomials))

    def defining_polynomials(self):
        """
        Return the polynomials that define this scheme as a subscheme
        of its ambient space.

        OUTPUT:

        A tuple of polynomials in the coordinate ring of the ambient
        space.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x^2-y*z, x^3+z^3])
            sage: S.defining_polynomials()
            (x^2 - y*z, x^3 + z^3)
        """
        return self.__polys

    def defining_ideal(self):
        """
        Return the ideal that defines this scheme as a subscheme
        of its ambient space.

        OUTPUT:

        An ideal in the coordinate ring of the ambient space.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x^2-y*z, x^3+z^3])
            sage: S.defining_ideal()
            Ideal (x^2 - y*z, x^3 + z^3) of Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        try:
            return self.__I
        except AttributeError:
            R = self.ambient_space().coordinate_ring()
            self.__I = R.ideal(self.defining_polynomials())
            return self.__I

    # Note: dimension must be implemented by the derived classes
    def codimension(self):
        r"""
        Return the codimension of the algebraic subscheme.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: PP.<x,y,z,w,v> = ProjectiveSpace(4,QQ)
            sage: V = PP.subscheme(x*y)
            sage: V.codimension()
            1
            sage: V.dimension()
            3
        """
        return self.ambient_space().dimension() - self.dimension()

    def irreducible_components(self):
        r"""
        Return the irreducible components of this algebraic scheme, as
        subschemes of the same ambient space.

        OUTPUT:

        an immutable sequence of irreducible subschemes of the ambient
        space of this scheme

        The components are cached.

        EXAMPLES:

        We define what is clearly a union of four hypersurfaces in
        `\P^4_{\QQ}` then find the irreducible components::

            sage: PP.<x,y,z,w,v> = ProjectiveSpace(4,QQ)
            sage: V = PP.subscheme( (x^2 - y^2 - z^2)*(w^5 -  2*v^2*z^3)* w * (v^3 - x^2*z) )
            sage: V.irreducible_components()
            [
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
            w,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
            x^2 - y^2 - z^2,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
            x^2*z - v^3,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
            w^5 - 2*z^3*v^2
            ]

        We verify that the irrelevant ideal isn't accidently returned
        (see :trac:`6920`)::

            sage: PP.<x,y,z,w> = ProjectiveSpace(3,QQ)
            sage: f = x^3 + y^3 + z^3 + w^3
            sage: R = f.parent()
            sage: I = [f] + [f.derivative(zz) for zz in PP.gens()]
            sage: V = PP.subscheme(I)
            sage: V.irreducible_components()
            [
            <BLANKLINE>
            ]

        The same polynomial as above defines a scheme with a
        nontrivial irreducible component in affine space (instead of
        the empty scheme as above)::

            sage: AA.<x,y,z,w> = AffineSpace(4,QQ)
            sage: V = AA.subscheme(I)
            sage: V.irreducible_components()
            [
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              w,
              z,
              y,
              x
            ]
        """
        try:
            return self.__irreducible_components
        except AttributeError:
            pass
        I = self.defining_ideal()
        P = I.associated_primes()
        if self.is_projective():
            # In the projective case, we must exclude the prime ideals
            # that contain the irrelevant ideal, which is the ideal
            # generated by the variables, which are the gens of the
            # base ring.
            G = I.ring().gens()
            # We make a list of ideals with the property that "any"
            # of the elements of G are not in the ideal.
            P = [J for J in P if any(g not in J for g in G)]

        A = self.ambient_space()
        C = Sequence([A.subscheme(X) for X in P], check=False, cr=True)
        C.sort()
        C.set_immutable()
        self.__irreducible_components = C
        return C

    def is_irreducible(self):
        r"""
        Return whether this subscheme is or is not irreducible.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: K = QuadraticField(-3)
            sage: P.<x,y,z,w,t,u> = ProjectiveSpace(K, 5)
            sage: X = P.subscheme([x*y - z^2 - K.0*t^2, t*w*x + y*z^2 - u^3])
            sage: X.is_irreducible()
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme([(y + x - z)^2])
            sage: X.is_irreducible()
            False

        ::

            sage: A.<x,y,z,w> = AffineSpace(GF(17), 4)
            sage: X = A.subscheme([x*y*z^2 - x*y*z*w - z*w^2 + w^3, x^3*y*z*w - x*y^3*z - x^2*y*z*w \
            - x^2*w^3 + y^2*w^2 + x*w^3])
            sage: X.is_irreducible()
            False
        """
        return self.defining_ideal().is_prime()

    def Jacobian_matrix(self):
        r"""
        Return the matrix `\frac{\partial f_i}{\partial x_j}` of
        (formal) partial derivatives.

        OUTPUT:

        A matrix of polynomials.

        EXAMPLES::

            sage: P3.<w,x,y,z> = ProjectiveSpace(3, QQ)
            sage: twisted_cubic = P3.subscheme(matrix([[w, x, y],[x, y, z]]).minors(2))
            sage: twisted_cubic.Jacobian_matrix()
            [   y -2*x    w    0]
            [   z   -y   -x    w]
            [   0    z -2*y    x]
            
        This example addresses ticket :trac:`20512`::
        
            sage: X = P3.subscheme([])
            sage: X.Jacobian_matrix().base_ring() == P3.coordinate_ring()
            True
        """
        R = self.ambient_space().coordinate_ring()
        l = self.defining_polynomials()
        if len(l) == 0:
            return sage.matrix.constructor.Matrix(R, 0)
        return jacobian(l, R.gens())

    def Jacobian(self):
        r"""
        Return the Jacobian ideal.

        This is the ideal generated by

        * the `d\times d` minors of the Jacobian matrix, where `d` is
          the :meth:`codimension` of the algebraic scheme, and

        * the defining polynomials of the algebraic scheme. Note that
          some authors do not include these in the definition of the
          Jacobian ideal. An example of a reference that does include
          the defining equations is [LazarsfeldJacobian].

        OUTPUT:

        An ideal in the coordinate ring of the ambient space.

        REFERENCES:

        ..  [LazarsfeldJacobian]
            Robert Lazarsfeld:
            Positivity in algebraic geometry II;
            Positivity for Vector Bundles, and Multiplier Ideals,
            page 181.

        EXAMPLES::

            sage: P3.<w,x,y,z> = ProjectiveSpace(3, QQ)
            sage: twisted_cubic = P3.subscheme(matrix([[w, x, y],[x, y, z]]).minors(2))
            sage: twisted_cubic.Jacobian()
            Ideal (-x^2 + w*y, -x*y + w*z, -y^2 + x*z, x*z, -2*w*z, w*y, 3*w*y, -2*w*x,
            w^2, y*z, -2*x*z, w*z, 3*w*z, -2*w*y, w*x, z^2, -2*y*z, x*z, 3*x*z, -2*w*z,
            w*y) of Multivariate Polynomial Ring in w, x, y, z over Rational Field
            sage: twisted_cubic.defining_ideal()
            Ideal (-x^2 + w*y, -x*y + w*z, -y^2 + x*z) of Multivariate Polynomial Ring
            in w, x, y, z over Rational Field
        
        This example addresses ticket :trac:`20512`::
        
            sage: X = P3.subscheme([])
            sage: X.Jacobian() == P3.coordinate_ring().unit_ideal()
            True
        """
        d = self.codimension()
        minors = self.Jacobian_matrix().minors(d)
        I = self.defining_ideal()
        minors = tuple([ I.reduce(m) for m in minors ])
        return I.ring().ideal(I.gens() + minors)

    def reduce(self):
        r"""
        Return the corresponding reduced algebraic space associated to this
        scheme.

        EXAMPLES: First we construct the union of a doubled and tripled
        line in the affine plane over `\QQ` ::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: X = A.subscheme([(x-1)^2*(x-y)^3]); X
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^5 - 3*x^4*y + 3*x^3*y^2 - x^2*y^3 - 2*x^4 + 6*x^3*y
              - 6*x^2*y^2 + 2*x*y^3 + x^3 - 3*x^2*y + 3*x*y^2 - y^3
            sage: X.dimension()
            1

        Then we compute the corresponding reduced scheme::

            sage: Y = X.reduce(); Y
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^2 - x*y - x + y

        Finally, we verify that the reduced scheme `Y` is the union
        of those two lines::

            sage: L1 = A.subscheme([x-1]); L1
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x - 1
            sage: L2 = A.subscheme([x-y]); L2
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x - y
            sage: W = L1.union(L2); W             # taken in ambient space
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^2 - x*y - x + y
            sage: Y == W
            True
        """
        try:
            return self._reduce
        except AttributeError:
            r = self.defining_ideal().radical()
            A = self.ambient_space()
            V = A.subscheme(r)
            V._reduce = V       # so knows it is already reduced!
            self._reduce = V
            return V

    def union(self, other):
        """
        Return the scheme-theoretic union of self and other in their common
        ambient space.

        EXAMPLES: We construct the union of a line and a tripled-point on
        the line.

        ::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: I = ideal([x,y])^3
            sage: P = A.subscheme(I)
            sage: L = A.subscheme([y-1])
            sage: S = L.union(P); S
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
            y^4 - y^3,
            x*y^3 - x*y^2,
            x^2*y^2 - x^2*y,
            x^3*y - x^3
            sage: S.dimension()
            1
            sage: S.reduce()
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
            y^2 - y,
            x*y - x

        We can also use the notation "+" for the union::

            sage: A.subscheme([x]) + A.subscheme([y^2 - (x^3+1)])
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
            x^4 - x*y^2 + x

        Saving and loading::

            sage: loads(S.dumps()) == S
            True
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            raise TypeError("other (=%s) must be a closed algebraic subscheme of an ambient space"%other)
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError("other (=%s) must be in the same ambient space as self"%other)
        return A.subscheme(self.defining_ideal().intersection(other.defining_ideal()))

    __add__ = union

    def intersection(self, other):
        """
        Return the scheme-theoretic intersection of self and other in their
        common ambient space.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: X = A.subscheme([x^2-y])
            sage: Y = A.subscheme([y])
            sage: X.intersection(Y)
            Closed subscheme of Affine Space of dimension 2 over Integer Ring defined by:
              x^2 - y,
              y
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            raise TypeError("other (=%s) must be a closed algebraic subscheme of an ambient space"%other)
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError("other (=%s) must be in the same ambient space as self"%other)
        return A.subscheme(self.defining_ideal() + other.defining_ideal())

    def complement(self, other=None):
        """
        Return the scheme-theoretic complement other - self, where
        self and other are both closed algebraic subschemes of the
        same ambient space.

        If other is unspecified, it is taken to be the ambient space
        of self.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: X = A.subscheme([x+y-z])
            sage: Y = A.subscheme([x-y+z])
            sage: Y.complement(X)
            Quasi-affine subscheme X - Y of Affine Space of
            dimension 3 over Integer Ring, where X is defined by:
              x + y - z
            and Y is defined by:
              x - y + z
            sage: Y.complement()
            Quasi-affine subscheme X - Y of Affine Space of
            dimension 3 over Integer Ring, where X is defined by:
              (no polynomials)
            and Y is defined by:
              x - y + z
            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: X = P.subscheme([x^2+y^2+z^2])
            sage: Y = P.subscheme([x*y+y*z+z*x])
            sage: Y.complement(X)
            Quasi-projective subscheme X - Y of Projective Space of
            dimension 2 over Rational Field, where X is defined by:
              x^2 + y^2 + z^2
            and Y is defined by:
              x*y + x*z + y*z
            sage: Y.complement(P)
            Quasi-projective subscheme X - Y of Projective Space of
            dimension 2 over Rational Field, where X is defined by:
              (no polynomials)
            and Y is defined by:
              x*y + x*z + y*z
        """
        A = self.ambient_space()
        if other is None:
            other = A.subscheme([])
        elif not isinstance(other, AlgebraicScheme_subscheme):
            if other == A:
                other = A.subscheme([])
            else:
                raise TypeError("Argument other (=%s) must be a closed algebraic subscheme of an ambient space"%other)
        if other.ambient_space() != A:
            raise ValueError("other (=%s) must be in the same ambient space as self"%other)
        return AlgebraicScheme_quasi(other, self)

    def rational_points(self, bound=0, F=None):
        """
        Return the rational points on the algebraic subscheme.

        EXAMPLES:

        Enumerate over a projective scheme over a number field::

            sage: u = QQ['u'].0
            sage: K.<v> = NumberField(u^2 + 3)
            sage: A.<x,y> = ProjectiveSpace(K,1)
            sage: X=A.subscheme(x^2 - y^2)
            sage: X.rational_points(3)
            [(-1 : 1), (1 : 1)]

        One can enumerate points up to a given bound on a projective scheme
        over the rationals::

            sage: E = EllipticCurve('37a')
            sage: E.rational_points(bound=8)
            [(-1 : -1 : 1), (-1 : 0 : 1), (0 : -1 : 1), (0 : 0 : 1), (0 : 1 : 0), (1/4 : -5/8 : 1),
            (1/4 : -3/8 : 1), (1 : -1 : 1), (1 : 0 : 1), (2 : -3 : 1), (2 : 2 : 1)]

        For a small finite field, the complete set of points can be
        enumerated. ::

            sage: Etilde = E.base_extend(GF(3))
            sage: Etilde.rational_points()
            [(0 : 0 : 1), (0 : 1 : 0), (0 : 2 : 1), (1 : 0 : 1),
             (1 : 2 : 1), (2 : 0 : 1), (2 : 2 : 1)]

        The class of hyperelliptic curves does not (yet) support
        desingularization of the places at infinity into two points::

            sage: FF = FiniteField(7)
            sage: P.<x> = PolynomialRing(FiniteField(7))
            sage: C = HyperellipticCurve(x^8+x+1)
            sage: C.rational_points()
            [(0 : 1 : 0), (0 : 1 : 1), (0 : 6 : 1), (2 : 0 : 1),
             (4 : 0 : 1), (6 : 1 : 1), (6 : 6 : 1)]

        TODO:

        1. The above algorithms enumerate all projective points and
           test whether they lie on the scheme; Implement a more naive
           sieve at least for covers of the projective line.

        2. Implement Stoll's model in weighted projective space to
           resolve singularities and find two points (1 : 1 : 0) and
           (-1 : 1 : 0) at infinity.
        """
        if F is None:
            F = self.base_ring()
        X = self.base_extend(F)(F)
        if F in NumberFields() or F == ZZ:
            try:
                return X.points(bound) # checks for proper bound done in points functions
            except TypeError:
                raise TypeError("Unable to enumerate points over %s."%F)
        try:
            return X.points()
        except TypeError:
            raise TypeError("Unable to enumerate points over %s."%F)

    def change_ring(self, R):
        r"""
        Returns a new algebraic subscheme which is this subscheme coerced to ``R``.

        INPUT:

        - ``R`` -- ring or morphism.

        OUTPUT:

        - A new algebraic subscheme which is this subscheme coerced to ``R``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: X = P.subscheme([3*x^2-y^2])
            sage: H = Hom(X,X)
            sage: X.change_ring(GF(3))
            Closed subscheme of Projective Space of dimension 1 over Finite Field of size 3 defined by:
            -y^2

        ::

            sage: K.<w> = QuadraticField(2)
            sage: R.<z> = K[]
            sage: L.<v> = K.extension(z^3-5)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: X = P.subscheme(x - w*y)
            sage: X.change_ring(L)
            Closed subscheme of Projective Space of dimension 1 over Number Field in v with
            defining polynomial z^3 - 5 over its base field defined by:
              x + (-w)*y

        ::

            sage: K.<w> = QuadraticField(2)
            sage: R.<z> = K[]
            sage: L.<v> = K.extension(z^3-5)
            sage: P.<x,y,z> = AffineSpace(L,3)
            sage: X = P.subscheme([x-w*y, z^2-v*x])
            sage: emb = L.embeddings(QQbar)
            sage: X.change_ring(emb[0])
            Closed subscheme of Affine Space of dimension 3 over Algebraic Field
            defined by:
              x + (-1.414213562373095? + 0.?e-16*I)*y,
              z^2 + (0.8549879733383485? + 1.480882609682365?*I)*x

        ::

            sage: K.<w> = QuadraticField(2)
            sage: R.<z> = K[]
            sage: L.<v> = K.extension(z^3-5)
            sage: P.<x,y,z> = AffineSpace(L,3)
            sage: X = P.subscheme([x-w*y, z^2-v*x])
            sage: emb = L.embeddings(QQbar)
            sage: X.change_ring(emb[1])
            Closed subscheme of Affine Space of dimension 3 over Algebraic Field
            defined by:
              x + (-1.414213562373095? + 0.?e-16*I)*y,
              z^2 + (0.8549879733383485? - 1.480882609682365?*I)*x

        ::

            sage: K.<w> = QuadraticField(-3)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: X = P.subscheme(x-w*y)
            sage: X.change_ring(CC)
            Closed subscheme of Projective Space of dimension 1 over Complex Field
            with 53 bits of precision defined by:
              x + (-1.73205080756888*I)*y

        ::

            sage: K.<w> = QuadraticField(3)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: X = P.subscheme(x-w*y)
            sage: X.change_ring(RR)
            Closed subscheme of Projective Space of dimension 1 over Real Field
            with 53 bits of precision defined by:
              x - 1.73205080756888*y

        ::

            sage: K.<v> = CyclotomicField(7)
            sage: O = K.maximal_order()
            sage: P.<x,y> = ProjectiveSpace(O, 1)
            sage: X = P.subscheme([x^2+O(v)*y^2])
            sage: X.change_ring(CC)
            Closed subscheme of Projective Space of dimension 1 over Complex Field
            with 53 bits of precision defined by:
              x^2 + (0.623489801858734 + 0.781831482468030*I)*y^2
            sage: X.change_ring(K).change_ring(K.embeddings(QQbar)[0])
            Closed subscheme of Projective Space of dimension 1 over Algebraic Field defined by:
              x^2 + (-0.9009688679024191? - 0.4338837391175581?*I)*y^2

        ::

            sage: R.<x> = QQ[]
            sage: f = x^6-2
            sage: L.<b> = NumberField(f, embedding=f.roots(CC)[2][0])
            sage: A.<x,y> = AffineSpace(L, 2)
            sage: H = Hom(A,A)
            sage: X = A.subscheme([b*x^2, y^2])
            sage: X.change_ring(CC)
            Closed subscheme of Affine Space of dimension 2 over Complex Field with
            53 bits of precision defined by:
              (-0.561231024154687 - 0.972080648619833*I)*x^2,
              y^2
        """
        K = self.base_ring()
        AS = self.ambient_space()
        new_AS = AS.change_ring(R)
        I = [f.change_ring(R) for f in self.defining_polynomials()]
        return(new_AS.subscheme(I))

    def weil_restriction(self):
        r"""
        Compute the Weil restriction of this variety over some extension
        field. If the field is a finite field, then this computes
        the Weil restriction to the prime subfield.

        A Weil restriction of scalars - denoted `Res_{L/k}` - is a
        functor which, for any finite extension of fields `L/k` and
        any algebraic variety `X` over `L`, produces another
        corresponding variety `Res_{L/k}(X)`, defined over `k`. It is
        useful for reducing questions about varieties over large
        fields to questions about more complicated varieties over
        smaller fields.

        This function does not compute this Weil restriction directly
        but computes on generating sets of polynomial ideals:

        Let `d` be the degree of the field extension `L/k`, let `a` a
        generator of `L/k` and `p` the minimal polynomial of
        `L/k`. Denote this ideal by `I`.

        Specifically, this function first maps each variable `x` to
        its representation over `k`: `\sum_{i=0}^{d-1} a^i x_i`. Then
        each generator of `I` is evaluated over these representations
        and reduced modulo the minimal polynomial `p`. The result is
        interpreted as a univariate polynomial in `a` and its
        coefficients are the new generators of the returned ideal.

        If the input and the output ideals are radical, this is
        equivalent to the statement about algebraic varieties above.

        OUTPUT: Affine subscheme - the Weil restriction of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^5-2)
            sage: R.<x> = K[]
            sage: L.<v> = K.extension(x^2+1)
            sage: A.<x,y> = AffineSpace(L,2)
            sage: X = A.subscheme([y^2-L(w)*x^3-v])
            sage: X.weil_restriction()
            Closed subscheme of Affine Space of dimension 4 over Number Field in w
            with defining polynomial x^5 - 2 defined by:
              (-w)*z0^3 + (3*w)*z0*z1^2 + z2^2 - z3^2,
              (-3*w)*z0^2*z1 + (w)*z1^3 + 2*z2*z3 - 1
            sage: X.weil_restriction().ambient_space() is A.weil_restriction()
            True

        ::

            sage: A.<x,y,z> = AffineSpace(GF(5^2,'t'),3)
            sage: X = A.subscheme([y^2-x*z, z^2+2*y])
            sage: X.weil_restriction()
            Closed subscheme of Affine Space of dimension 6 over Finite Field of
            size 5 defined by:
              z2^2 - 2*z3^2 - z0*z4 + 2*z1*z5,
              2*z2*z3 + z3^2 - z1*z4 - z0*z5 - z1*z5,
              z4^2 - 2*z5^2 + 2*z2,
              2*z4*z5 + z5^2 + 2*z3
        """
        try:
            X = self.__weil_restriction
        except AttributeError:
            L = self.base_ring()
            if L.is_finite():
                d = L.degree()
            else:
                d = L.relative_degree()

            if d == 1:
                X = self
            else:
                A = self.ambient_space().weil_restriction()
                I = self.defining_ideal().weil_restriction()
                X = A.subscheme(I)
            self.__weil_restriction = X
        return X

#*******************************************************************
# Affine varieties
#*******************************************************************
class AlgebraicScheme_subscheme_affine(AlgebraicScheme_subscheme):
    """
    Construct an algebraic subscheme of affine space.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`~sage.schemes.affine.affine_space.AffineSpace_generic.subscheme`
        method of :class:`affine space
        <sage.schemes.affine.affine_space.AffineSpace_generic>`.

    INPUT:

    - ``A`` -- ambient :class:`affine space
      <sage.schemes.affine.affine_space.AffineSpace_generic>`

    - ``polynomials`` -- single polynomial, ideal or iterable of
      defining polynomials.

    EXAMPLES::

        sage: A3.<x, y, z> = AffineSpace(3, QQ)
        sage: A3.subscheme([x^2-y*z])
        Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
          x^2 - y*z

    TESTS::

        sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_affine
        sage: AlgebraicScheme_subscheme_affine(A3, [x^2-y*z])
        Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
          x^2 - y*z
    """

    def _morphism(self, *args, **kwds):
        """
        A morphism between two schemes in your category, usually defined via
        polynomials. Your morphism class should derive from
        :class:`SchemeMorphism_polynomial`. These morphisms will usually be
        elements of the Hom-set
        :class:`~sage.schemes.generic.homset.SchemeHomset_generic`.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(2, ZZ)
            sage: P2._morphism(P2.Hom(P2), [x,y,z])
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x : y : z)

        """
        return self.ambient_space()._morphism(*args, **kwds)

    def dimension(self):
        """
        Return the dimension of the affine algebraic subscheme.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: A.subscheme([]).dimension()
            2
            sage: A.subscheme([x]).dimension()
            1
            sage: A.subscheme([x^5]).dimension()
            1
            sage: A.subscheme([x^2 + y^2 - 1]).dimension()
            1
            sage: A.subscheme([x*(x-1), y*(y-1)]).dimension()
            0

        Something less obvious::

            sage: A.<x,y,z,w> = AffineSpace(4, QQ)
            sage: X = A.subscheme([x^2, x^2*y^2 + z^2, z^2 - w^2, 10*x^2 + w^2 - z^2])
            sage: X
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              x^2,
              x^2*y^2 + z^2,
              z^2 - w^2,
              10*x^2 - z^2 + w^2
            sage: X.dimension()
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.defining_ideal().dimension()
            return self.__dimension

    def projective_embedding(self, i=None, PP=None):
        """
        Returns a morphism from this affine scheme into an ambient
        projective space of the same dimension.

        The codomain of this morphism is the projective closure of this affine scheme in ``PP``,
        if given, or otherwise in a new projective space that is constructed.

        INPUT:

        -  ``i`` -- integer (default: dimension of self = last
           coordinate) determines which projective embedding to compute. The
           embedding is that which has a 1 in the i-th coordinate, numbered
           from 0.

        -  ``PP`` -- (default: None) ambient projective space, i.e., ambient space
            of codomain of morphism; this is constructed if it is not given.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: S = A.subscheme([x*y-z])
            sage: S.projective_embedding()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Integer Ring defined by:
              x*y - z
              To:   Closed subscheme of Projective Space of dimension 3 over Integer Ring defined by:
              x0*x1 - x2*x3
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x : y : z : 1)

        ::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: P = ProjectiveSpace(3,ZZ,'u')
            sage: S = A.subscheme([x^2-y*z])
            sage: S.projective_embedding(1,P)
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Integer
            Ring defined by:
              x^2 - y*z
              To:   Closed subscheme of Projective Space of dimension 3 over Integer
            Ring defined by:
              u0^2 - u2*u3
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x : 1 : y : z)

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([y - x^2, z - x^3])
            sage: X.projective_embedding()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Rational
            Field defined by:
              -x^2 + y,
              -x^3 + z
              To:   Closed subscheme of Projective Space of dimension 3 over
            Rational Field defined by:
              x0^2 - x1*x3,
              x0*x1 - x2*x3,
              x1^2 - x0*x2
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x : y : z : 1)

        """
        AA = self.ambient_space()
        n = AA.dimension_relative()
        if i is None:
            try:
                i = self._default_embedding_index
            except AttributeError:
                i = int(n)
        else:
            i = int(i)
        if i < 0 or i > n:
            raise ValueError("Argument i (=%s) must be between 0 and %s, inclusive"%(i, n))
        try:
            phi = self.__projective_embedding[i]
            #assume that if you've passed in a new ambient projective space
            #you want to override the existing embedding
            if PP is None or phi.codomain().ambient_space() == PP:
                return(phi)
        except AttributeError:
            self.__projective_embedding = {}
        except KeyError:
            pass
        if PP is None:
            PP = AA.projective_embedding(i).codomain()
        elif PP.dimension_relative() != n:
            raise ValueError("Projective Space must be of dimension %s"%(n))
        PR = PP.coordinate_ring()
        # Groebner basis w.r.t. a graded monomial order computed here to ensure
        # after homogenization, the basis elements will generate the defining
        # ideal of the projective closure of this affine subscheme
        R = AA.coordinate_ring()
        G = self.defining_ideal().groebner_basis()
        v = list(PP.gens())
        z = v.pop(i)
        phi = R.hom(v,PR)
        v.append(z)
        X = PP.subscheme([phi(f).homogenize(i) for f in G])
        v = list(R.gens())
        v.insert(i, R(1))
        phi = self.hom(v, X)
        self.__projective_embedding[i] = phi
        return phi

    def projective_closure(self, i=None, PP=None):
        r"""
        Return the projective closure of this affine subscheme.

        INPUT:

        - ``i`` -- (default: None) determines the embedding to use to compute the projective
          closure of this affine subscheme. The embedding used is the one which has a 1 in the
          i-th coordinate, numbered from 0.

        -  ``PP`` -- (default: None) ambient projective space, i.e., ambient space
           of codomain of morphism; this is constructed if it is not given.

        OUTPUT:

        - a projective subscheme.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ,4)
            sage: X = A.subscheme([x^2 - y, x*y - z, y^2 - w, x*z - w, y*z - x*w, z^2 - y*w])
            sage: X.projective_closure()
            Closed subscheme of Projective Space of dimension 4 over Rational Field
            defined by:
              x0^2 - x1*x4,
              x0*x1 - x2*x4,
              x1^2 - x3*x4,
              x0*x2 - x3*x4,
              x1*x2 - x0*x3,
              x2^2 - x1*x3

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: P.<a,b,c,d> = ProjectiveSpace(QQ, 3)
            sage: X = A.subscheme([z - x^2 - y^2])
            sage: X.projective_closure(1, P).ambient_space() == P
            True
        """
        return self.projective_embedding(i, PP).codomain()

    def is_smooth(self, point=None):
        r"""
        Test whether the algebraic subscheme is smooth.

        INPUT:

        - ``point`` -- A point or ``None`` (default). The point to
          test smoothness at.

        OUTPUT:

        Boolean. If no point was specified, returns whether the
        algebraic subscheme is smooth everywhere. Otherwise,
        smoothness at the specified point is tested.

        EXAMPLES::

            sage: A2.<x,y> = AffineSpace(2,QQ)
            sage: cuspidal_curve = A2.subscheme([y^2-x^3])
            sage: cuspidal_curve
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -x^3 + y^2
            sage: smooth_point = cuspidal_curve.point([1,1])
            sage: smooth_point in cuspidal_curve
            True
            sage: singular_point = cuspidal_curve.point([0,0])
            sage: singular_point in cuspidal_curve
            True
            sage: cuspidal_curve.is_smooth(smooth_point)
            True
            sage: cuspidal_curve.is_smooth(singular_point)
            False
            sage: cuspidal_curve.is_smooth()
            False
        """
        R = self.ambient_space().coordinate_ring()
        if not point is None:
            self._check_satisfies_equations(point)
            point_subs = dict(zip(R.gens(), point))
            Jac = self.Jacobian().subs(point_subs)
            return not Jac.is_zero()

        # testing smoothness everywhere tends to be expensive
        try:
            return self._smooth
        except AttributeError:
            pass
        sing_dim = self.Jacobian().dimension()
        self._smooth = (sing_dim == -1)
        return self._smooth

    def intersection_multiplicity(self, X, P):
        r"""
        Return the intersection multiplicity of this subscheme and the subscheme ``X`` at the point ``P``.

        The intersection of this subscheme with ``X`` must be proper, that is `\mathrm{codim}(self\cap
        X) = \mathrm{codim}(self) + \mathrm{codim}(X)`, and must also be finite. We use Serre's Tor
        formula to compute the intersection multiplicity. If `I`, `J` are the defining ideals of ``self``, ``X``,
        respectively, then this is `\sum_{i=0}^{\infty}(-1)^i\mathrm{length}(\mathrm{Tor}_{\mathcal{O}_{A,p}}^{i}
        (\mathcal{O}_{A,p}/I,\mathcal{O}_{A,p}/J))` where `A` is the affine ambient space of these subschemes.

        INPUT:

        - ``X`` -- subscheme in the same ambient space as this subscheme.

        - ``P`` -- a point in the intersection of this subscheme with ``X``.

        OUTPUT: An integer.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3 - x^2], A)
            sage: D = Curve([y^2 + x^3], A)
            sage: Q = A([0,0])
            sage: C.intersection_multiplicity(D, Q)
            4

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^6 - 3*a^5 + 5*a^4 - 5*a^3 + 5*a^2 - 3*a + 1)
            sage: A.<x,y,z,w> = AffineSpace(K, 4)
            sage: X = A.subscheme([x*y, y*z + 7, w^3 - x^3])
            sage: Y = A.subscheme([x - z^3 + z + 1])
            sage: Q = A([0, -7*b^5 + 21*b^4 - 28*b^3 + 21*b^2 - 21*b + 14, -b^5 + 2*b^4 - 3*b^3 \
            + 2*b^2 - 2*b, 0])
            sage: X.intersection_multiplicity(Y, Q)
            3

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([z^2 - 1])
            sage: Y = A.subscheme([z - 1, y - x^2])
            sage: Q = A([1,1,1])
            sage: X.intersection_multiplicity(Y, Q)
            Traceback (most recent call last):
            ...
            TypeError: the intersection of this subscheme and (=Closed subscheme of Affine Space of dimension 3
            over Rational Field defined by: z - 1, -x^2 + y) must be proper and finite

        ::

            sage: A.<x,y,z,w,t> = AffineSpace(QQ, 5)
            sage: X = A.subscheme([x*y, t^2*w, w^3*z])
            sage: Y = A.subscheme([y*w + z])
            sage: Q = A([0,0,0,0,0])
            sage: X.intersection_multiplicity(Y, Q)
            Traceback (most recent call last):
            ...
            TypeError: the intersection of this subscheme and (=Closed subscheme of Affine Space of dimension 5
            over Rational Field defined by: y*w + z) must be proper and finite
        """
        AA = self.ambient_space()
        if AA != X.ambient_space():
            raise TypeError("this subscheme and (=%s) must be defined in the same ambient space"%X)
        W = self.intersection(X)
        try:
            W._check_satisfies_equations(P)
        except TypeError:
            raise TypeError("(=%s) must be a point in the intersection of this subscheme and (=%s)"%(P,X))
        if AA.dimension() != self.dimension() + X.dimension() or W.dimension() != 0:
            raise TypeError("the intersection of this subscheme and (=%s) must be proper and finite"%X)
        I = self.defining_ideal()
        J = X.defining_ideal()
        # move P to the origin and localize
        chng_coords = [AA.gens()[i] + P[i] for i in range(AA.dimension_relative())]
        R = AA.coordinate_ring().change_ring(order="negdegrevlex")
        Iloc = R.ideal([f(chng_coords) for f in I.gens()])
        Jloc = R.ideal([f(chng_coords) for f in J.gens()])
        # compute the intersection multiplicity with Serre's Tor formula using Singular
        singular.lib("homolog.lib")
        i = 0
        s = 0
        t = sum(singular.Tor(i, Iloc, Jloc).std().hilb(2).sage())
        while t != 0:
            s = s + ((-1)**i)*t
            i = i + 1
            t = sum(singular.Tor(i, Iloc, Jloc).std().hilb(2).sage())
        return s

    def multiplicity(self, P):
        r"""
        Return the multiplicity of ``P`` on this subscheme.

        This is computed as the multiplicity of the local ring of this subscheme corresponding to ``P``. This
        subscheme must be defined over a field. An error is raised if ``P`` is not a point on this subscheme.

        INPUT:

        - ``P`` -- a point on this subscheme.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: X = A.subscheme([z*y - x^7, w - 2*z])
            sage: Q1 = A([1,1/3,3,6])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = A([0,0,0,0])
            sage: X.multiplicity(Q2)
            2

        ::

            sage: A.<x,y,z,w,v> = AffineSpace(GF(23), 5)
            sage: C = A.curve([x^8 - y, y^7 - z, z^3 - 1, w^5 - v^3])
            sage: Q = A([22,1,1,0,0])
            sage: C.multiplicity(Q)
            3

        ::

            sage: K.<a> = QuadraticField(-1)
            sage: A.<x,y,z,w,t> = AffineSpace(K, 5)
            sage: X = A.subscheme([y^7 - x^2*z^5 + z^3*t^8 - x^2*y^4*z - t^8])
            sage: Q1 = A([1,1,0,1,-1])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = A([0,0,0,-a,0])
            sage: X.multiplicity(Q2)
            7
        """
        if not self.base_ring() in Fields():
            raise TypeError("subscheme must be defined over a field")

        # Check whether P is a point on this subscheme
        try:
            P = self(P)
        except TypeError:
            raise TypeError("(=%s) is not a point on (=%s)"%(P,self))

        # Apply a linear change of coordinates to self so that P is sent to the origin
        # and then compute the multiplicity of the local ring of the translated subscheme 
        # corresponding to the point (0,...,0)
        AA = self.ambient_space()
        chng_coords = [AA.gens()[i] + P[i] for i in range(AA.dimension_relative())]
        R = AA.coordinate_ring().change_ring(order='negdegrevlex')
        I = R.ideal([f(chng_coords) for f in self.defining_polynomials()])
        return singular.mult(singular.std(I)).sage()


#*******************************************************************
# Projective varieties
#*******************************************************************
class AlgebraicScheme_subscheme_projective(AlgebraicScheme_subscheme):
    """
    Construct an algebraic subscheme of projective space.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`~sage.schemes.projective.projective_space.ProjectiveSpace_field.subscheme`
        method of :class:`projective space
        <sage.schemes.projective.projective_space.ProjectiveSpace_field>`.

    INPUT:

    - ``A`` -- ambient :class:`projective space
      <sage.schemes.projective.projective_space.ProjectiveSpace_field>`.

    - ``polynomials`` -- single polynomial, ideal or iterable of
      defining homogeneous polynomials.

    EXAMPLES::

        sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
        sage: P.subscheme([x^2-y*z])
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x^2 - y*z

    TESTS::

        sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
        sage: AlgebraicScheme_subscheme_projective(P, [x^2-y*z])
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x^2 - y*z
    """

    def _morphism(self, *args, **kwds):
        r"""
        Construct a morphism determined by action on points of ``self``.

        For internal use only.

        INPUT:

        - same as for
          :class:`~sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space`.

        OUTPUT:

        - :class:`~sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space`.

        TESTS::

            sage: P1.<x,y> = ProjectiveSpace(1,QQ)
            sage: P2 = ProjectiveSpace(2,QQ)
            sage: H12 = P1.Hom(P2)
            sage: H12([x^2,x*y, y^2])    # indirect doctest
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                (x^2 : x*y : y^2)
            sage: P1._morphism(H12, [x^2,x*y, y^2])
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                (x^2 : x*y : y^2)
        """
        return self.ambient_space()._morphism(*args, **kwds)

    def dimension(self):
        """
        Return the dimension of the projective algebraic subscheme.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(2, QQ)
            sage: P2.subscheme([]).dimension()
            2
            sage: P2.subscheme([x]).dimension()
            1
            sage: P2.subscheme([x^5]).dimension()
            1
            sage: P2.subscheme([x^2 + y^2 - z^2]).dimension()
            1
            sage: P2.subscheme([x*(x-z), y*(y-z)]).dimension()
            0

        Something less obvious::

            sage: P3.<x,y,z,w,t> = ProjectiveSpace(4, QQ)
            sage: X = P3.subscheme([x^2, x^2*y^2 + z^2*t^2, z^2 - w^2, 10*x^2 + w^2 - z^2])
            sage: X
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
              x^2,
              x^2*y^2 + z^2*t^2,
              z^2 - w^2,
              10*x^2 - z^2 + w^2
            sage: X.dimension()
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.defining_ideal().dimension() - 1
            return self.__dimension

    def affine_patch(self, i, AA = None):
        r"""
        Return the `i^{th}` affine patch of this projective scheme.
        This is the intersection with this `i^{th}` affine patch of
        its ambient space.

        INPUT:

        - ``i`` -- integer between 0 and dimension of self, inclusive.

        - ``AA`` -- (default: None) ambient affine space, this is constructed
            if it is not given.

        OUTPUT:

        An affine algebraic scheme with fixed
        :meth:`embedding_morphism` equal to the default
        :meth:`projective_embedding` map`.

        EXAMPLES::

            sage: PP = ProjectiveSpace(2, QQ, names='X,Y,Z')
            sage: X,Y,Z = PP.gens()
            sage: C = PP.subscheme(X^3*Y + Y^3*Z + Z^3*X)
            sage: U = C.affine_patch(0)
            sage: U
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x0^3*x1 + x1^3 + x0
            sage: U.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x0^3*x1 + x1^3 + x0
              To:   Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              X^3*Y + Y^3*Z + X*Z^3
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1 : x0 : x1)
            sage: U.projective_embedding() is U.embedding_morphism()
            True

        ::

            sage: A.<x,y,z> = AffineSpace(QQ,3)
            sage: X = A.subscheme([x-y*z])
            sage: Y = X.projective_embedding(1).codomain()
            sage: Y.affine_patch(1,A).ambient_space() == A
            True

        ::

            sage: P.<u,v,w> = ProjectiveSpace(2,ZZ)
            sage: S = P.subscheme([u^2-v*w])
            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: S.affine_patch(1, A)
            Closed subscheme of Affine Space of dimension 2 over Integer Ring
            defined by:
              x^2 - y
        """
        i = int(i)   # implicit type checking
        PP = self.ambient_space()
        n = PP.dimension_relative()
        if i < 0 or i > n:
            raise ValueError("Argument i (= %s) must be between 0 and %s."%(i, n))
        try:
            A = self.__affine_patches[i]
            #assume that if you've passed in a new ambient affine space
            #you want to override the existing patch
            if AA is None or A.ambient_space() == AA:
                return self.__affine_patches[i]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        if AA is None:
            AA = PP.affine_patch(i)
        elif AA.dimension_relative() != n:
            raise ValueError("Affine Space must be of the dimension %s"%(n))
        phi = AA.projective_embedding(i, PP)
        polys = self.defining_polynomials()
        xi = phi.defining_polynomials()
        U = AA.subscheme([ f(xi) for f in polys ])
        U._default_embedding_index = i
        phi = U.projective_embedding(i, PP)
        self.__affine_patches[i] = U
        U._embedding_morphism = phi
        return U

    def _best_affine_patch(self, point):
        r"""
        Return the best affine patch of the ambient projective space.

        The "best" affine patch is where you end up dividing by the
        homogeneous coordinate with the largest absolutue
        value. Division by small numbers is numerically unstable.

        INPUT:

        - ``point`` -- a point of the algebraic subscheme.

        OUTPUT:

        Integer. The index of the patch. See :meth:`affine_patch`.

        EXAMPLES::

            sage: P.<x,y,z>= ProjectiveSpace(QQ,2)
            sage: S = P.subscheme(x+2*y+3*z)
            sage: S._best_affine_patch(P.point([0,-3,2]))
            1
            sage: S._best_affine_patch([0,-3,2])
            1

        TESTS::

            sage: F = GF(3)
            sage: P.<x,y,z>= ProjectiveSpace(F,2)
            sage: S._best_affine_patch([0,1,2])
            2
        """
        point = list(point)
        try:
            abs_point = [abs(_) for _ in point]
        except ArithmeticError:
            # our base ring does not know abs
            abs_point = point
        # find best patch
        i_max = 0
        p_max = abs_point[i_max]
        for i in range(1,len(point)):
            if abs_point[i]>p_max:
                i_max = i
                p_max = abs_point[i_max]
        return i_max

    def neighborhood(self, point):
        r"""
        Return an affine algebraic subscheme isomorphic to a
        neighborhood of the ``point``.

        INPUT:

        - ``point`` -- a point of the projective subscheme.

        OUTPUT:

        An affine algebraic scheme (polynomial equations in affine
        space) ``result`` such that

        * :meth:`embedding_morphism
          <AlgebraicScheme.embedding_morphism>` is an isomorphism to a
          neighborhood of ``point``

        * :meth:`embedding_center <AlgebraicScheme.embedding_center>`
          is mapped to ``point``.

        EXAMPLES::

            sage: P.<x,y,z>= ProjectiveSpace(QQ,2)
            sage: S = P.subscheme(x+2*y+3*z)
            sage: s = S.point([0,-3,2]); s
            (0 : -3/2 : 1)
            sage: patch = S.neighborhood(s); patch
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x0 + 3*x1
            sage: patch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x0 + 3*x1
              To:   Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x + 2*y + 3*z
              Defn: Defined on coordinates by sending (x0, x1) to
                    (x0 : -3/2 : x1 + 1)
            sage: patch.embedding_center()
            (0, 0)
            sage: patch.embedding_morphism()([0,0])
            (0 : -3/2 : 1)
            sage: patch.embedding_morphism()(patch.embedding_center())
            (0 : -3/2 : 1)
        """
        point = list(point)
        self._check_satisfies_equations(point)
        PP = self.ambient_space()
        n = PP.dimension()
        i = self._best_affine_patch(point)

        patch_cover = PP.affine_patch(i)
        R = patch_cover.coordinate_ring()

        phi = list(point)
        for j in range(0,i):
            phi[j] = phi[j] + R.gen(j)
        for j in range(i,n):
            phi[j+1] = phi[j+1] + R.gen(j)

        pullback_polys = [f(phi) for f in self.defining_polynomials()]
        patch = patch_cover.subscheme(pullback_polys)
        patch_hom = patch.hom(phi,self)
        patch._embedding_center = patch.point([0]*n)
        patch._embedding_morphism = patch_hom
        return patch

    def is_smooth(self, point=None):
        r"""
        Test whether the algebraic subscheme is smooth.

        INPUT:

        - ``point`` -- A point or ``None`` (default). The point to
          test smoothness at.

        OUTPUT:

        Boolean. If no point was specified, returns whether the
        algebraic subscheme is smooth everywhere. Otherwise,
        smoothness at the specified point is tested.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: cuspidal_curve = P2.subscheme([y^2*z-x^3])
            sage: cuspidal_curve
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              -x^3 + y^2*z
            sage: cuspidal_curve.is_smooth([1,1,1])
            True
            sage: cuspidal_curve.is_smooth([0,0,1])
            False
            sage: cuspidal_curve.is_smooth()
            False
            sage: P2.subscheme([y^2*z-x^3+z^3+1/10*x*y*z]).is_smooth()
            True

        TESTS::

            sage: H = P2.subscheme(x)
            sage: H.is_smooth()  # one of the few cases where the cone over the subvariety is smooth
            True
        """
        if not point is None:
            self._check_satisfies_equations(point)
            R = self.ambient_space().coordinate_ring()
            point_subs = dict(zip(R.gens(), point))
            Jac = self.Jacobian().subs(point_subs)
            return not Jac.is_zero()

        # testing smoothness everywhere tends to be expensive
        try:
            return self._smooth
        except AttributeError:
            pass
        sing_dim = self.Jacobian().dimension()
        # We really test the affine cone here; the origin is always a
        # singular point:
        self._smooth = (sing_dim <= 0)
        return self._smooth

    def orbit(self, f, N):
        r"""
        Returns the orbit of this scheme by ``f``.

        If `N` is an integer it returns `[self,f(self),\ldots,f^N(self)]`.
        If `N` is a list or tuple `N=[m,k]` it returns `[f^m(self),\ldots,f^k(self)`].

        INPUT:

        - ``f`` -- a :class:`SchemeMorphism_polynomial` with ``self`` in ``f.domain()``

        - ``N`` -- a non-negative integer or list or tuple of two non-negative integers

        OUTPUT:

        - a list of projective subschemes

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: H = End(P)
            sage: f = H([(x-2*y)^2,(x-2*z)^2,(x-2*w)^2,x^2])
            sage: f.orbit(P.subscheme([x]),5)
            [Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
               x,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
               w,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
               z - w,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
               y - z,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
               x - y,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
               x - w]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(PS, P1)
            sage: f = H([x^2, y^2])
            sage: X = PS.subscheme([x-y])
            sage: X.orbit(f,2)
            Traceback (most recent call last):
            ...
            TypeError: map must be an endomorphism for iteration

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: X = PS.subscheme([x-y])
            sage: X.orbit(f,[-1,2])
            Traceback (most recent call last):
            ...
            TypeError: orbit bounds must be non-negative
        """
        if not f.is_endomorphism():
            raise TypeError("map must be an endomorphism for iteration")
        if not isinstance(N,(list,tuple)):
            N = [0,N]
        N[0] = ZZ(N[0])
        N[1] = ZZ(N[1])
        if N[0] < 0 or N[1] < 0:
            raise TypeError("orbit bounds must be non-negative")
        if N[0] > N[1]:
            return([])

        Q = self
        for i in range(1, N[0]+1):
            Q = f(Q)
        Orb = [Q]

        for i in range(N[0]+1, N[1]+1):
            Q = f(Q)
            Orb.append(Q)
        return(Orb)

    def nth_iterate(self, f, n):
        r"""
        The nth forward image of this scheme by the map ``f``.

        INPUT:

        - ``f`` -- a SchmemMorphism_polynomial with ``self`` in ``f.domain()``

        - ``n`` -- a positive integer.

        OUTPUT:

        - A subscheme in ``f.codomain()``

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: H = End(P)
            sage: f = H([y^2, z^2, x^2, w^2])
            sage: f.nth_iterate(P.subscheme([x-w,y-z]), 3)
            Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
              y - z,
              x - w

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: X = PS.subscheme([x-y])
            sage: X.nth_iterate(f,-2)
            Traceback (most recent call last):
            ...
            TypeError: must be a forward orbit

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: P2.<u,v,w>=ProjectiveSpace(QQ, 2)
            sage: H = Hom(PS, P2)
            sage: f = H([x^2, y^2, z^2])
            sage: X = PS.subscheme([x-y])
            sage: X.nth_iterate(f,2)
            Traceback (most recent call last):
            ...
            TypeError: map must be an endomorphism for iteration

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: X = PS.subscheme([x-y])
            sage: X.nth_iterate(f,2.5)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
        """
        n = ZZ(n)
        if n < 0:
            raise TypeError("must be a forward orbit")
        return self.orbit(f,[n,n+1])[0]

    def _forward_image(self, f, check = True):
        """
        Compute the forward image of this subscheme by the morphism ``f``.

        The forward image is computed through elimination and ``f`` must be
        a morphism for this to be well defined.
        In particular, let $X = V(h_1,\ldots, h_t)$ and define the ideal
        $I = (h_1,\ldots,h_t,y_0-f_0(\bar{x}), \ldots, y_n-f_n(\bar{x}))$.
        Then the elimination ideal $I_{n+1} = I \cap K[y_0,\ldots,y_n]$ is a homogeneous
        ideal and $self(X) = V(I_{n+1})$.

        INPUT:

        - ``f`` -- a map whose domain contains ``self``

        - ``check`` -- Boolean, if `False` no input checking is done

        OUTPUT:

         - a subscheme in the codomain of ``f``.

        EXAMPLES::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(PS)
            sage: f = H([x^2, y^2-2*z^2, z^2])
            sage: X = PS.subscheme(y-2*z)
            sage: X._forward_image(f)
            Closed subscheme of Projective Space of dimension 2 over Rational Field
            defined by:
              y - 2*z

        ::

            sage: set_verbose(None)
            sage: PS.<x,y,z,w> = ProjectiveSpace(ZZ, 3)
            sage: H = End(PS)
            sage: f = H([y^2, x^2, w^2, z^2])
            sage: X = PS.subscheme([z^2+y*w, x-w])
            sage: f(X)
            Closed subscheme of Projective Space of dimension 3 over Integer Ring
            defined by:
              y - z,
              x*z - w^2

        ::

            sage: PS.<x,y,z,w> = ProjectiveSpace(CC, 3)
            sage: H = End(PS)
            sage: f = H([x^2 + y^2, y^2, z^2-y^2, w^2])
            sage: X = PS.subscheme([z-2*w])
            sage: f(X)
            Closed subscheme of Projective Space of dimension 3 over Complex Field
            with 53 bits of precision defined by:
              y + z + (-4.00000000000000)*w

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R), 2)
            sage: H = End(P)
            sage: f = H([x^2 + 2*y*z, t^2*y^2, z^2])
            sage: f([t^2*y-z])
            Closed subscheme of Projective Space of dimension 2 over Fraction Field
            of Univariate Polynomial Ring in t over Rational Field defined by:
              y + (-1/t^2)*z

        ::

            sage: set_verbose(-1)
            sage: PS.<x,y,z> = ProjectiveSpace(Qp(3), 2)
            sage: H = End(PS)
            sage: f = H([x^2,2*y^2,z^2])
            sage: X = PS.subscheme([2*x-y,z])
            sage: f(X)
            Closed subscheme of Projective Space of dimension 2 over 3-adic Field
            with capped relative precision 20 defined by:
              z,
              x + (1 + 3^2 + 3^4 + 3^6 + 3^8 + 3^10 + 3^12 + 3^14 + 3^16 + 3^18 +
            O(3^20))*y

        ::

            sage: R.<y0,y1,y2,y3> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R), 2)
            sage: H = End(P)
            sage: f = H([y0*x^2+y1*z^2, y2*y^2+y3*z^2, z^2])
            sage: X = P.subscheme(x*z)
            sage: X._forward_image(f)
            Closed subscheme of Projective Space of dimension 2 over Fraction Field
            of Multivariate Polynomial Ring in y0, y1, y2, y3 over Rational Field
            defined by:
              x*z + (-y1)*z^2

            ::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P5.<z0,z1,z2,z3,z4,z5> = ProjectiveSpace(QQ, 5)
            sage: H = Hom(P2, P5)
            sage: f = H([x^2,x*y,x*z,y^2,y*z,z^2]) #Veronese map
            sage: X = P2.subscheme([])
            sage: f(X)
            Closed subscheme of Projective Space of dimension 5 over Rational Field
            defined by:
              -z4^2 + z3*z5,
              -z2*z4 + z1*z5,
              -z2*z3 + z1*z4,
              -z2^2 + z0*z5,
              -z1*z2 + z0*z4,
              -z1^2 + z0*z3

            ::

            sage: P2.<x,y,z>=ProjectiveSpace(QQ, 2)
            sage: P3.<u,v,w,t>=ProjectiveSpace(QQ, 3)
            sage: H = Hom(P2, P3)
            sage: X = P2.subscheme([x-y,x-z])
            sage: f = H([x^2,y^2,z^2,x*y])
            sage: f(X)
            Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
              w - t,
              v - t,
              u - t

            ::

            sage: P1.<u,v> = ProjectiveSpace(QQ, 1)
            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P2,P1)
            sage: f = H([x^2,y*z])
            sage: X = P2.subscheme([x-y])
            sage: f(X)
            Traceback (most recent call last):
            ...
            TypeError: map must be a morphism

            ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = End(PS)
            sage: f = H([x^3, x*y^2, x*z^2])
            sage: X = PS.subscheme([x-y])
            sage: X._forward_image(f)
            Traceback (most recent call last):
            ...
            TypeError: map must be a morphism

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(QQ, 1)
            sage: Y = P1.subscheme([u-v])
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: Y._forward_image(f)
            Traceback (most recent call last):
            ...
            TypeError: subscheme must be in ambient space of domain of map
        """
        dom = f.domain()
        codom = f.codomain()
        if check:
            if not f.is_morphism():
                raise TypeError("map must be a morphism")
            if self.ambient_space() != dom:
                raise TypeError("subscheme must be in ambient space of domain of map")
        CR_dom = dom.coordinate_ring()
        CR_codom = codom.coordinate_ring()
        n = CR_dom.ngens()
        m = CR_codom.ngens()
        #can't call eliminate if the base ring is polynomial so we do it ourselves
        #with a lex ordering
        R = PolynomialRing(f.base_ring(), n+m, 'tempvar', order = 'lex')
        Rvars = R.gens()[0 : n]
        phi = CR_dom.hom(Rvars,R)
        zero = n*[0]
        psi = R.hom(zero + list(CR_codom.gens()),CR_codom)
        #set up ideal
        L = R.ideal([phi(t) for t in self.defining_polynomials()] + [R.gen(n+i) - phi(f[i]) for i in range(m)])
        G = L.groebner_basis() #eliminate
        newL = []
        #get only the elimination ideal portion
        for i in range (len(G)-1,0,-1):
            v = G[i].variables()
            if all([Rvars[j] not in v for j in range(n)]):
                newL.append(psi(G[i]))
        return(codom.subscheme(newL))

    def preimage(self, f, k=1, check=True):
        r"""
        The subscheme that maps to this scheme by the map `f^k`.

        In particular, `f^{-k}(V(h_1,\ldots,h_t)) = V(h_1 \circ f^k, \ldots, h_t \circ f^k)`.
        Map must be a morphism and also must be an endomorphism for `k > 1`.

        INPUT:

        - ``f`` - a map whose codomain contains this scheme

        - ``k`` - a positive integer

        - ``check`` -- Boolean, if ``False`` no input checking is done

        OUTPUT:

        - a subscheme in the domain of ``f``.

        Examples::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = End(PS)
            sage: f = H([y^2, x^2, z^2])
            sage: X = PS.subscheme([x-y])
            sage: X.preimage(f)
            Closed subscheme of Projective Space of dimension 2 over Integer Ring
            defined by:
              -x^2 + y^2

        ::

            sage: P.<x,y,z,w,t> = ProjectiveSpace(QQ, 4)
            sage: H = End(P)
            sage: f = H([x^2-y^2, y^2, z^2, w^2, t^2+w^2])
            sage: f.rational_preimages(P.subscheme([x-z, t^2, w-t]))
            Closed subscheme of Projective Space of dimension 4 over Rational Field
            defined by:
              x^2 - y^2 - z^2,
              w^4 + 2*w^2*t^2 + t^4,
              -t^2

        ::

            sage: P1.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P3.<u,v,w,t> = ProjectiveSpace(QQ, 3)
            sage: H = Hom(P1, P3)
            sage: X = P3.subscheme([u-v, 2*u-w, u+t])
            sage: f = H([x^2,y^2, x^2+y^2, x*y])
            sage: X.preimage(f)
            Closed subscheme of Projective Space of dimension 1 over Rational Field
            defined by:
              x^2 - y^2,
              x^2 - y^2,
              x^2 + x*y

        ::

            sage: P1.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P3.<u,v,w,t> = ProjectiveSpace(QQ, 3)
            sage: H = Hom(P3, P1)
            sage: X = P1.subscheme([x-y])
            sage: f = H([u^2, v^2])
            sage: X.preimage(f)
            Traceback (most recent call last):
            ...
            TypeError: map must be a morphism

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = End(PS)
            sage: f = H([x^2, x^2, x^2])
            sage: X = PS.subscheme([x-y])
            sage: X.preimage(f)
            Traceback (most recent call last):
            ...
            TypeError: map must be a morphism

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: P1.<u,v> = ProjectiveSpace(ZZ, 1)
            sage: Y = P1.subscheme([u^2-v^2])
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: Y.preimage(f)
            Traceback (most recent call last):
            ...
            TypeError: subscheme must be in ambient space of codomain

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: Y = P.subscheme([x-y])
            sage: H = End(P)
            sage: f = H([x^2, y^2, z^2])
            sage: Y.preimage(f, k=2)
            Closed subscheme of Projective Space of dimension 2 over Rational Field
            defined by:
              x^4 - y^4
        """
        dom = f.domain()
        codom = f.codomain()
        if check:
            if not f.is_morphism():
                raise TypeError("map must be a morphism")
            if self.ambient_space() != codom:
                raise TypeError("subscheme must be in ambient space of codomain")
            k = ZZ(k)
            if k <= 0:
                raise ValueError("k (=%s) must be a positive integer"%(k))
            if k > 1 and not f.is_endomorphism():
                raise TypeError("map must be an endomorphism")
        R = codom.coordinate_ring()
        F = f.nth_iterate_map(k)
        dict = {R.gen(i): F[i] for i in range(codom.dimension_relative()+1)}
        return(dom.subscheme([t.subs(dict) for t in self.defining_polynomials()]))

    def dual(self):
        r"""
        Return the projective dual of the given subscheme of projective space.

        INPUT:

        - ``X`` -- A subscheme of projective space. At present, ``X`` is
          required to be an irreducible and reduced hypersurface defined
          over `\QQ` or a finite field.

        OUTPUT:

        - The dual of ``X`` as a subscheme of the dual projective space.

        EXAMPLES:

        The dual of a smooth conic in the plane is also a smooth conic::

            sage: R.<x, y, z> = QQ[]
            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: I = R.ideal(x^2 + y^2 + z^2)
            sage: X = P.subscheme(I)
            sage: X.dual()
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              y0^2 + y1^2 + y2^2

        The dual of the twisted cubic curve in projective 3-space is a singular
        quartic surface. In the following example, we compute the dual of this
        surface, which by double duality is equal to the twisted cubic itself.
        The output is the twisted cubic as an intersection of three quadrics::

            sage: R.<x, y, z, w> = QQ[]
            sage: P.<x, y, z, w> = ProjectiveSpace(3, QQ)
            sage: I = R.ideal(y^2*z^2 - 4*x*z^3 - 4*y^3*w + 18*x*y*z*w - 27*x^2*w^2)
            sage: X = P.subscheme(I)
            sage: X.dual()
            Closed subscheme of Projective Space of dimension 3 over
            Rational Field defined by:
              y2^2 - y1*y3,
              y1*y2 - y0*y3,
              y1^2 - y0*y2

        The singular locus of the quartic surface in the last example
        is itself supported on a twisted cubic::

            sage: X.Jacobian().radical()
            Ideal (z^2 - 3*y*w, y*z - 9*x*w, y^2 - 3*x*z) of Multivariate
            Polynomial Ring in x, y, z, w over Rational Field

        An example over a finite field::

            sage: R = PolynomialRing(GF(61), 'a,b,c')
            sage: P.<a, b, c> = ProjectiveSpace(2, R.base_ring())
            sage: X = P.subscheme(R.ideal(a*a+2*b*b+3*c*c))
            sage: X.dual()
            Closed subscheme of Projective Space of dimension 2 over
            Finite Field of size 61 defined by:
            y0^2 - 30*y1^2 - 20*y2^2

        TESTS::

            sage: R = PolynomialRing(Qp(3), 'a,b,c')
            sage: P.<a, b, c> = ProjectiveSpace(2, R.base_ring())
            sage: X = P.subscheme(R.ideal(a*a+2*b*b+3*c*c))
            sage: X.dual()
            Traceback (most recent call last):
            ...
            NotImplementedError: base ring must be QQ or a finite field
        """
        from sage.libs.singular.function_factory import ff

        K = self.base_ring()
        if not(is_RationalField(K) or is_FiniteField(K)):
            raise NotImplementedError("base ring must be QQ or a finite field")
        I = self.defining_ideal()
        m = I.ngens()
        n = I.ring().ngens() - 1
        if (m != 1 or (n < 1) or I.is_zero()
            or I.is_trivial() or not I.is_prime()):
            raise NotImplementedError("At the present, the method is only"
                                      " implemented for irreducible and"
                                      " reduced hypersurfaces and the given"
                                      " list of generators for the ideal must"
                                      " have exactly one element.")
        R = PolynomialRing(K, 'x', n + 1)
        Pd = sage.schemes.projective.projective_space.ProjectiveSpace(n, K, 'y')
        Rd = Pd.coordinate_ring()
        x = R.gens()
        y = Rd.gens()
        S = PolynomialRing(K, x + y + ('t',))
        if S.has_coerce_map_from(I.ring()):
            T = PolynomialRing(K, 'w', n + 1)
            I_S = (I.change_ring(T)).change_ring(S)
        else:
            I_S = I.change_ring(S)
        f_S = I_S.gens()[0]
        z = S.gens()
        J = I_S
        for i in range(n + 1):
            J = J + S.ideal(z[-1] * f_S.derivative(z[i]) - z[i + n + 1])

        sat = ff.elim__lib.sat

        max_ideal = S.ideal(z[n + 1: 2 * n + 2])
        J_sat_gens = sat(J, max_ideal)[0]
        J_sat = S.ideal(J_sat_gens)
        L = J_sat.elimination_ideal(z[0: n + 1] + (z[-1],))
        return Pd.subscheme(L.change_ring(Rd))

    def Chow_form(self):
        r"""
        Returns the Chow form associated to this subscheme.

        For a `k`-dimensional subvariety of `\mathbb{P}^N` of degree `D`.
        The `(N-k-1)`-dimensional projective linear subspaces of `\mathbb{P}^N`
        meeting `X` form a hypersurface in the Grassmannian `G(N-k-1,N)`.
        The homogeneous form of degree `D` defining this hypersurface in Plucker
        coordinates is called the Chow form of `X`.

        The base ring needs to be a number field, finite field, or `\QQbar`.

        ALGORITHM:

        For a `k`-dimension subscheme `X` consider the `k+1` linear forms
        `l_i = u_{i0}x_0 + \cdots + u_{in}x_n`. Let `J` be the ideal in the
        polynomial ring `K[x_i,u_{ij}]` defined by the equations of `X` and the `l_i`.
        Let `J'` be the saturation of `J` with respect to the irrelevant ideal of
        the ambient projective space of `X`. The elimination ideal `I = J' \cap K[u_{ij}]`
        is a principal ideal, let `R` be its generator. The Chow form is obtained by
        writing `R` as a polynomial in Plucker coordinates (i.e. bracket polynomials).
        [DalbecSturmfels].

        OUTPUT: a homogeneous polynomial.

        REFERENCES:

        .. [DalbecSturmfels] J. Dalbec and B. Sturmfels. Invariant methods in discrete and computational geometry,
           chapter Introduction to Chow forms, pages 37-58. Springer Netherlands, 1994.

        EXAMPLES::

            sage: P.<x0,x1,x2,x3> = ProjectiveSpace(GF(17), 3)
            sage: X = P.subscheme([x3+x1,x2-x0,x2-x3])
            sage: X.Chow_form()
            t0 - t1 + t2 + t3

        ::

            sage: P.<x0,x1,x2,x3> = ProjectiveSpace(QQ,3)
            sage: X = P.subscheme([x3^2 -101*x1^2 - 3*x2*x0])
            sage: X.Chow_form()
            t0^2 - 101*t2^2 - 3*t1*t3

        ::

            sage: P.<x0,x1,x2,x3>=ProjectiveSpace(QQ,3)
            sage: X = P.subscheme([x0*x2-x1^2, x0*x3-x1*x2, x1*x3-x2^2])
            sage: Ch = X.Chow_form(); Ch
            t2^3 + 2*t2^2*t3 + t2*t3^2 - 3*t1*t2*t4 - t1*t3*t4 + t0*t4^2 + t1^2*t5
            sage: Y = P.subscheme_from_Chow_form(Ch, 1); Y
            Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
              x2^2*x3 - x1*x3^2,
              -x2^3 + x0*x3^2,
              -x2^2*x3 + x1*x3^2,
              x1*x2*x3 - x0*x3^2,
              3*x1*x2^2 - 3*x0*x2*x3,
              -2*x1^2*x3 + 2*x0*x2*x3,
              -3*x1^2*x2 + 3*x0*x1*x3,
              x1^3 - x0^2*x3,
              x2^3 - x1*x2*x3,
              -3*x1*x2^2 + 2*x1^2*x3 + x0*x2*x3,
              2*x0*x2^2 - 2*x0*x1*x3,
              3*x1^2*x2 - 2*x0*x2^2 - x0*x1*x3,
              -x0*x1*x2 + x0^2*x3,
              -x0*x1^2 + x0^2*x2,
              -x1^3 + x0*x1*x2,
              x0*x1^2 - x0^2*x2
            sage: I = Y.defining_ideal()
            sage: I.saturation(I.ring().ideal(list(I.ring().gens())))[0]
            Ideal (x2^2 - x1*x3, x1*x2 - x0*x3, x1^2 - x0*x2) of Multivariate
            Polynomial Ring in x0, x1, x2, x3 over Rational Field
        """
        I = self.defining_ideal()
        P = self.ambient_space()
        R = P.coordinate_ring()
        N = P.dimension()+1
        d = self.dimension()
        #create the ring for the generic linear hyperplanes
        # u0x0 + u1x1 + ...
        SS = PolynomialRing(R.base_ring(), 'u', N*(d+1), order='lex')
        vars = SS.variable_names() + R.variable_names()
        S = PolynomialRing(R.base_ring(), vars, order='lex')
        n = S.ngens()
        newcoords = [S.gen(n-N+t) for t in range(N)]
        #map the generators of the subscheme into the ring with the hyperplane variables
        phi = R.hom(newcoords,S)
        phi(self.defining_polynomials()[0])
        #create the dim(X)+1 linear hyperplanes
        l = []
        for i in range(d+1):
            t = 0
            for j in range(N):
                t += S.gen(N*i + j)*newcoords[j]
            l.append(t)
        #intersect the hyperplanes with X
        J = phi(I) + S.ideal(l)
        #saturate the ideal with respect to the irrelevant ideal
        J2 = J.saturation(S.ideal([phi(t) for t in R.gens()]))[0]
        #elimante the original variables to be left with the hyperplane coefficients 'u'
        E = J2.elimination_ideal(newcoords)
        #create the plucker coordinates
        D = binomial(N,N-d-1) #number of plucker coordinates
        tvars = [str('t') + str(i) for i in range(D)] #plucker coordinates
        T = PolynomialRing(R.base_ring(), tvars+list(S.variable_names()), order='lex')
        L = []
        coeffs = [T.gen(i) for i in range(0+len(tvars), N*(d+1)+len(tvars))]
        M = matrix(T,d+1,N,coeffs)
        i = 0
        for c in M.minors(d+1):
            L.append(T.gen(i)-c)
            i += 1
        #create the ideal that we can use for eliminating to get a polynomial
        #in the plucker coordinates (brackets)
        br = T.ideal(L)
        #create a mapping into a polynomial ring over the plucker coordinates
        #and the hyperplane coefficients
        psi = S.hom(coeffs + [0 for i in range(N)],T)
        E2 = T.ideal([psi(u) for u in E.gens()] +br)
        #eliminate the hyperplane coefficients
        CH = E2.elimination_ideal(coeffs)
        #CH should be a principal ideal, but because of the relations among
        #the plucker coordinates, the elimination will probably have several generators

        #get the relations among the plucker coordinates
        rel = br.elimination_ideal(coeffs)
        #reduce CH with respect to the relations
        reduced = []
        for f in CH.gens():
            reduced.append(f.reduce(rel))
        #find the principal generator

        #polynomial ring in just the plucker coordinates
        T2 = PolynomialRing(R.base_ring(), tvars)
        alp = T.hom(tvars + (N*(d+1) +N)*[0], T2)
        #get the degress of the reduced generators of CH
        degs = [u.degree() for u in reduced]
        mind = max(degs)
        #need the smallest degree form that did not reduce to 0
        for d in degs:
            if d < mind and d >0:
                mind = d
        ind = degs.index(mind)
        CF = reduced[ind] #this should be the Chow form of X
        #check that it is correct (i.e., it is a principal generator for CH + the relations)
        rel2 = rel + [CF]
        assert all([f in rel2 for f in CH.gens()]), "did not find a principal generator"
        return(alp(CF))

    def degree(self):
        r"""
        Return the degree of this projective subscheme.

        If `P(t) = a_{m}t^m + \ldots + a_{0}` is the the Hilbert polynomial of this subscheme, then
        the degree is `a_{m}m!`.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y,z,w,t,u> = ProjectiveSpace(QQ, 5)
            sage: X = P.subscheme([x^7 + x*y*z*t^4 - u^7])
            sage: X.degree()
            7

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(GF(13), 3)
            sage: X = P.subscheme([y^3 - w^3, x + 7*z])
            sage: X.degree()
            3

        ::

            sage: P.<x,y,z,w,u> = ProjectiveSpace(QQ, 4)
            sage: C = P.curve([x^7 - y*z^3*w^2*u, w*x^2 - y*u^2, z^3 + y^3])
            sage: C.degree()
            63
        """
        P = self.defining_ideal().hilbert_polynomial()
        return P.leading_coefficient()*P.degree().factorial()

    def intersection_multiplicity(self, X, P):
        r"""
        Return the intersection multiplicity of this subscheme and the subscheme ``X`` at the point ``P``.

        This uses the intersection_multiplicity function for affine subschemes on affine patches of this subscheme
        and ``X`` that contain ``P``.

        INPUT:

        - ``X`` -- subscheme in the same ambient space as this subscheme.

        - ``P`` -- a point in the intersection of this subscheme with ``X``.

        OUTPUT: An integer.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: C = Curve([x^4 - z^2*y^2], P)
            sage: D = Curve([y^4*z - x^5 - x^3*z^2], P)
            sage: Q1 = P([0,1,0])
            sage: C.intersection_multiplicity(D, Q1)
            4
            sage: Q2 = P([0,0,1])
            sage: C.intersection_multiplicity(D, Q2)
            6

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^4 + 1)
            sage: P.<x,y,z,w> = ProjectiveSpace(K, 3)
            sage: X = P.subscheme([x^2 + y^2 - z*w])
            sage: Y = P.subscheme([y*z - x*w, z - w])
            sage: Q1 = P([b^2,1,0,0])
            sage: X.intersection_multiplicity(Y, Q1)
            1
            sage: Q2 = P([1/2*b^3-1/2*b,1/2*b^3-1/2*b,1,1])
            sage: X.intersection_multiplicity(Y, Q2)
            1

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: X = P.subscheme([x^2 - z^2, y^3 - w*x^2])
            sage: Y = P.subscheme([w^2 - 2*x*y + z^2, y^2 - w^2])
            sage: Q = P([1,1,-1,1])
            sage: X.intersection_multiplicity(Y, Q)
            Traceback (most recent call last):
            ...
            TypeError: the intersection of this subscheme and (=Closed subscheme of Affine Space of dimension 3
            over Rational Field defined by:
              x1^2 + x2^2 - 2*x0,
              x0^2 - x2^2) must be proper and finite
        """
        try:
            self.ambient_space()(P)
        except TypeError:
            raise TypeError("(=%s) must be a point in the ambient space of this subscheme and (=%s)"%(P,X))
        # find an affine chart of the ambient space of this curve that contains P
        n = self.ambient_space().dimension_relative()
        for i in range(n + 1):
            if P[i] != 0:
                break
        X1 = self.affine_patch(i)
        X2 = X.affine_patch(i)
        return X1.intersection_multiplicity(X2, X1(P.dehomogenize(i)))

    def multiplicity(self, P):
        r"""
        Return the multiplicity of ``P`` on this subscheme.

        This is computed as the multiplicity of the corresponding point on an affine patch of this subscheme
        that contains ``P``. This subscheme must be defined over a field. An error is returned if ``P``
        not a point on this subscheme.

        INPUT:

        - ``P`` -- a point on this subscheme.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: P.<x,y,z,w,t> = ProjectiveSpace(QQ, 4)
            sage: X = P.subscheme([y^2 - x*t, w^7 - t*w*x^5 - z^7])
            sage: Q1 = P([0,0,1,1,1])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = P([1,0,0,0,0])
            sage: X.multiplicity(Q2)
            3
            sage: Q3 = P([0,0,0,0,1])
            sage: X.multiplicity(Q3)
            7

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(CC, 3)
            sage: X = P.subscheme([z^5*x^2*w - y^8])
            sage: Q = P([2,0,0,1])
            sage: X.multiplicity(Q)
            5

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(GF(29), 3)
            sage: C = Curve([y^17 - x^5*w^4*z^8, x*y - z^2], P)
            sage: Q = P([3,0,0,1])
            sage: C.multiplicity(Q)
            8
        """
        if not self.base_ring() in Fields():
            raise TypeError("subscheme must be defined over a field")

        # Check whether P is a point on this subscheme
        try:
            P = self(P)
        except TypeError:
            raise TypeError("(=%s) is not a point on (=%s)"%(P,self))

        # Find an affine chart of the ambient space of self that contains P
        i = 0
        while(P[i] == 0):
            i = i + 1
        X = self.affine_patch(i)
        return X.multiplicity(X(P.dehomogenize(i)))

class AlgebraicScheme_subscheme_product_projective(AlgebraicScheme_subscheme_projective):

    @cached_method
    def segre_embedding(self, PP=None):
        r"""
        Return the Segre embedding of this subscheme into the appropriate projective
        space.

        INPUT:

        - ``PP`` -- (default: ``None``) ambient image projective space;
          this is constructed if it is not given.

        OUTPUT:

        Hom from this subscheme to the appropriate subscheme of projective space

        EXAMPLES::

            sage: X.<x,y,z,w,u,v> = ProductProjectiveSpaces([2,2], QQ)
            sage: P = ProjectiveSpace(QQ,8,'t')
            sage: L = (-w - v)*x + (-w*y - u*z)
            sage: Q = (-u*w - v^2)*x^2 + ((-w^2 - u*w + (-u*v - u^2))*y + (-w^2 - u*v)*z)*x + \
            ((-w^2 - u*w - u^2)*y^2 + (-u*w - v^2)*z*y + (-w^2 + (-v - u)*w)*z^2)
            sage: W = X.subscheme([L,Q])
            sage: phi = W.segre_embedding(P)
            sage: phi.codomain().ambient_space() == P
            True

        ::

            sage: PP.<x,y,u,v,s,t> = ProductProjectiveSpaces([1,1,1], CC)
            sage: PP.subscheme([]).segre_embedding()
            Scheme morphism:
              From: Closed subscheme of Product of projective spaces P^1 x P^1 x P^1
            over Complex Field with 53 bits of precision defined by:
              (no polynomials)
              To:   Closed subscheme of Projective Space of dimension 7 over Complex
            Field with 53 bits of precision defined by:
              -u5*u6 + u4*u7,
              -u3*u6 + u2*u7,
              -u3*u4 + u2*u5,
              -u3*u5 + u1*u7,
              -u3*u4 + u1*u6,
              -u3*u4 + u0*u7,
              -u2*u4 + u0*u6,
              -u1*u4 + u0*u5,
              -u1*u2 + u0*u3
              Defn: Defined by sending (x : y , u : v , s : t) to
                    (x*u*s : x*u*t : x*v*s : x*v*t : y*u*s : y*u*t : y*v*s : y*v*t).

        ::

            sage: PP.<x,y,z,u,v,s,t> = ProductProjectiveSpaces([2,1,1], ZZ)
            sage: PP.subscheme([x^3, u-v, s^2-t^2]).segre_embedding()
            Scheme morphism:
              From: Closed subscheme of Product of projective spaces P^2 x P^1 x P^1
            over Integer Ring defined by:
              x^3,
              u - v,
              s^2 - t^2
              To:   Closed subscheme of Projective Space of dimension 11 over
            Integer Ring defined by:
              u10^2 - u11^2,
              u9 - u11,
              u8 - u10,
              -u7*u10 + u6*u11,
              u6*u10 - u7*u11,
              u6^2 - u7^2,
              u5 - u7,
              u4 - u6,
              u3^3,
              -u3*u10 + u2*u11,
              u2*u10 - u3*u11,
              -u3*u6 + u2*u7,
              u2*u6 - u3*u7,
              u2*u3^2,
              u2^2 - u3^2,
              u1 - u3,
              u0 - u2
              Defn: Defined by sending (x : y : z , u : v , s : t) to
                    (x*u*s : x*u*t : x*v*s : x*v*t : y*u*s : y*u*t : y*v*s : y*v*t :
            z*u*s : z*u*t : z*v*s : z*v*t).
        """
        AS = self.ambient_space()
        CR = AS.coordinate_ring()
        N = AS.dimension_relative_components()
        M = prod([n+1 for n in N]) - 1

        vars = list(AS.coordinate_ring().variable_names()) + ['u' + str(i) for i in range(M+1)]
        from sage.rings.all import PolynomialRing
        R = PolynomialRing(AS.base_ring(), AS.ngens()+M+1, vars, order='lex')

        #set-up the elimination for the segre embedding
        mapping = []
        k = AS.ngens()
        index = AS.num_components()*[0]
        for count in range(M + 1):
            mapping.append(R.gen(k+count)-prod([CR(AS[i].gen(index[i])) for i in range(len(index))]))
            for i in range(len(index)-1, -1, -1):
                if index[i] == N[i]:
                    index[i] = 0
                else:
                    index[i] += 1
                    break #only increment once

        #change the defining ideal of the subscheme into the variables
        I = R.ideal(list(self.defining_polynomials()) + mapping)
        J = I.groebner_basis()
        s = set(R.gens()[:AS.ngens()])
        n = len(J)-1
        L = []
        while s.isdisjoint(J[n].variables()):
            L.append(J[n])
            n = n-1

        #create new subscheme
        if PP is None:
            from sage.schemes.projective.projective_space import ProjectiveSpace
            PS = ProjectiveSpace(self.base_ring(), M, R.gens()[AS.ngens():])
            Y = PS.subscheme(L)
        else:
            if PP.dimension_relative() != M:
                raise ValueError("projective space %s must be dimension %s")%(PP, M)
            S = PP.coordinate_ring()
            psi = R.hom([0]*k + list(S.gens()), S)
            L = [psi(l) for l in L]
            Y = PP.subscheme(L)

        #create embedding for points
        mapping = []
        index = AS.num_components()*[0]
        for count in range(M + 1):
            mapping.append(prod([CR(AS[i].gen(index[i])) for i in range(len(index))]))
            for i in range(len(index)-1, -1, -1):
                if index[i] == N[i]:
                    index[i] = 0
                else:
                    index[i] += 1
                    break #only increment once
        phi = self.hom(mapping, Y)

        return phi

    def dimension(self):
        r"""
        Return the dimension of the algebraic subscheme.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: X.<x,y,z,w,u,v> = ProductProjectiveSpaces([2,2],QQ)
            sage: L = (-w - v)*x + (-w*y - u*z)
            sage: Q = (-u*w - v^2)*x^2 + ((-w^2 - u*w + (-u*v - u^2))*y + (-w^2 - u*v)*z)*x + \
            ((-w^2 - u*w - u^2)*y^2 + (-u*w - v^2)*z*y + (-w^2 + (-v - u)*w)*z^2)
            sage: W = X.subscheme([L,Q])
            sage: W.dimension()
            2

        ::

            sage: PP.<x,y,z,u,v,s,t> = ProductProjectiveSpaces([2,1,1], QQ)
            sage: X = PP.subscheme([x^3, x^5+y^5, z^6, x*u-v*y, s^2-t^2])
            sage: X.dimension()
            -1

        ::

            sage: PP = ProductProjectiveSpaces([2,1,3], CC, 't')
            sage: PP.subscheme([]).dimension()
            6

        ::

            sage: PP = ProductProjectiveSpaces([1,3,1], ZZ, 't')
            sage: PP.subscheme([]).dimension()
            5

        ::

            sage: PP.<x,y,u,v,s,t> = ProductProjectiveSpaces([1,1,1], CC)
            sage: X = PP.subscheme([x^2-y^2, u-v, s^2-t^2])
            sage: X.dimension()
            0
        """
        try:
            return self.__dimension
        except AttributeError:
            try:
                #move to field to compute radical
                X = self.change_ring(FractionField(self.base_ring()))
                PP = X.ambient_space()
                I = X.defining_ideal().radical()
                #check if the irrelevant ideal of any component is in the radical
                if any([all([t in I for t in PS.gens()]) for PS in PP.components()]):
                    self.__dimension = -1
                else:
                    self.__dimension = I.dimension() - PP.num_components()
            except TypeError:  #cannot compute radical for this base ring
                phi = self.segre_embedding()
                self.__dimension = phi.codomain().defining_ideal().dimension() - 1
            return self.__dimension

    def is_smooth(self, point=None):
        r"""
        Test whether the algebraic subscheme is smooth.

        EXAMPLES::

            sage: X.<x,y,z,w,u,v> = ProductProjectiveSpaces([2,2],QQ)
            sage: L = (-w - v)*x + (-w*y - u*z)
            sage: Q = (-u*w - v^2)*x^2 + ((-w^2 - u*w + (-u*v - u^2))*y + (-w^2 - u*v)*z)*x + \
            ((-w^2 - u*w - u^2)*y^2 + (-u*w - v^2)*z*y + (-w^2 + (-v - u)*w)*z^2)
            sage: W = X.subscheme([L,Q])
            sage: W.is_smooth()
            Traceback (most recent call last):
            ...
            NotImplementedError: Not Implemented
        """
        raise NotImplementedError("Not Implemented")

    def affine_patch(self, I, return_embedding = False):
        r"""
        Return the `I^{th}` affine patch of this projective scheme
        where 'I' is a multi-index.

        INPUT:

        - ``I`` -- a list or tuple of positive integers

        - ``return_embedding`` -- Boolean, if true the projective embedding is also returned

        OUTPUT:

        - An affine algebraic scheme

        - An embedding into a product of projective space (optional)

        EXAMPLES::

            sage: PP.<x,y,z,w,u,v> = ProductProjectiveSpaces([3,1],QQ)
            sage: W = PP.subscheme([y^2*z-x^3,z^2-w^2,u^3-v^3])
            sage: W.affine_patch([0,1],True)
            (Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              x0^2*x1 - 1,
              x1^2 - x2^2,
              x3^3 - 1, Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              x0^2*x1 - 1,
              x1^2 - x2^2,
              x3^3 - 1
              To:   Closed subscheme of Product of projective spaces P^3 x P^1 over Rational Field defined by:
              -x^3 + y^2*z,
              z^2 - w^2,
              u^3 - v^3
              Defn: Defined on coordinates by sending (x0, x1, x2, x3) to
                    (1 : x0 : x1 : x2 , x3 : 1))
        """
        if not isinstance(I, (list, tuple)):
            raise TypeError('The argument I=%s must be a list or tuple of positice integers'%I)
        PP = self.ambient_space()
        N = PP.dimension_relative_components()
        if len(I) != len(N):
            raise ValueError('The argument I=%s must have %s entries'%(I,len(N)))
        I = tuple([int(i) for i in I])   # implicit type checking
        for i in range(len(I)):
            if I[i] < 0 or I[i] > N[i]:
                raise ValueError("Argument i (= %s) must be between 0 and %s."%(I[i], N[i]))
        #see if we've already created this affine patch
        try:
            if return_embedding:
                return self.__affine_patches[I]
            else:
                return self.__affine_patches[I][0]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        from sage.schemes.affine.affine_space import AffineSpace
        AA = AffineSpace(PP.base_ring(),sum(N),'x')
        v = list(AA.gens())
        # create the projective embedding
        index = 0
        for i in range(len(I)):
            v.insert(index+I[i],1)
            index += N[i]+1
        phi = AA.hom(v,self)
        #find the image of the subscheme
        polys = self.defining_polynomials()
        xi = phi.defining_polynomials()
        U = AA.subscheme([ f(xi) for f in polys ])
        phi = U.hom(v,self)
        self.__affine_patches.update({I:(U,phi)})
        if return_embedding:
            return U,phi
        else:
            return U

    def intersection_multiplicity(self, X, P):
        r"""
        Return the intersection multiplicity of this subscheme and the subscheme ``X`` at the point ``P``.

        This uses the intersection_multiplicity function for affine subschemes on affine patches of this subscheme
        and ``X`` that contain ``P``.

        INPUT:

        - ``X`` -- subscheme in the same ambient space as this subscheme.

        - ``P`` -- a point in the intersection of this subscheme with ``X``.

        OUTPUT: An integer.

        EXAMPLES:

        Multiplicity of a fixed point of the map `z^2 + \frac{1}{4}`::

            sage: PP.<x,y,u,v> = ProductProjectiveSpaces(QQ, [1,1])
            sage: G = PP.subscheme([(x^2 + 1/4*y^2)*v - y^2*u])
            sage: D = PP.subscheme([x*v - y*u])
            sage: G.intersection(D).rational_points()
            [(1 : 0 , 1 : 0), (1/2 : 1 , 1/2 : 1)]
            sage: Q = PP([1/2,1,1/2,1])
            sage: G.intersection_multiplicity(D, Q)
            2

        ::

            sage: F.<a> = GF(4)
            sage: PP.<x,y,z,u,v,w> = ProductProjectiveSpaces(F, [2,2])
            sage: X = PP.subscheme([z^5 + 3*x*y^4 + 8*y^5, u^2 - v^2])
            sage: Y = PP.subscheme([x^6 + z^6, w*z - v*y])
            sage: Q = PP([a,a+1,1,a,a,1])
            sage: X.intersection_multiplicity(Y, Q)
            16

        ::

            sage: PP.<x,y,z,u,v,w> = ProductProjectiveSpaces(QQ, [2,2])
            sage: X = PP.subscheme([x^2*u^3 + y*z*u*v^2, x - y])
            sage: Y = PP.subscheme([u^3 - w^3, x*v - y*w, z^3*w^2 - y^3*u*v])
            sage: Q = PP([0,0,1,0,1,0])
            sage: X.intersection_multiplicity(Y, Q)
            Traceback (most recent call last):
            ...
            TypeError: the intersection of this subscheme and (=Closed subscheme of Affine Space of dimension 4
            over Rational Field defined by: x2^3 - x3^3, -x1*x3 + x0, -x1^3*x2 + x3^2) must be proper and finite
        """
        PP = self.ambient_space()
        try:
            PP(P)
        except TypeError:
            raise TypeError("(=%s) must be a point in the ambient space of this subscheme and (=%s)"%(P,X))
        # find an affine chart of the ambient space of this subscheme that contains P
        indices = []
        aff_pt = []
        for i in range(PP.num_components()):
            Q = P[i]
            j = 0
            while Q[j] == 0:
                j = j + 1
            indices.append(j)
            T = list(Q)
            t = T.pop(j)
            aff_pt.extend([1/t*T[k] for k in range(PP.components()[i].dimension_relative())])
        X1 = self.affine_patch(indices)
        X2 = X.affine_patch(indices)
        return X1.intersection_multiplicity(X2, X1.ambient_space()(aff_pt))

    def multiplicity(self, P):
        r"""
        Return the multiplicity of ``P`` on this subscheme.

        This is computed as the multiplicity of the corresponding point on an affine patch of this subscheme
        that contains ``P``. This subscheme must be defined over a field. An error is returned if ``P``
        not a point on this subscheme.

        INPUT:

        - ``P`` -- a point on this subscheme.

        OUTPUT: an integer.

        EXAMPLES::

            sage: PP.<x,y,z,w> = ProductProjectiveSpaces(QQ, [1,1])
            sage: X = PP.subscheme([x^4*z^3 - y^4*w^3])
            sage: Q1 = PP([1,1,1,1])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = PP([0,1,1,0])
            sage: X.multiplicity(Q2)
            3

        ::

            sage: PP.<x,y,z,w,u> = ProductProjectiveSpaces(GF(11), [1,2])
            sage: X = PP.subscheme([x^7*u - y^7*z, u^6*x^2 - w^3*z^3*x*y - w^6*y^2])
            sage: Q1 = PP([1,0,10,1,0])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = PP([1,0,1,0,0])
            sage: X.multiplicity(Q2)
            4
        """
        PP = self.ambient_space()
        try:
            PP(P)
        except TypeError:
            raise TypeError("(=%s) must be a point in the ambient space of this subscheme and (=%s)"%(P,X))
        # find an affine chart of the ambient space of this subscheme that contains P
        indices = []
        aff_pt = []
        for i in range(PP.num_components()):
            Q = P[i]
            j = 0
            while Q[j] == 0:
                j = j + 1
            indices.append(j)
            T = list(Q)
            t = T.pop(j)
            aff_pt.extend([1/t*T[k] for k in range(PP.components()[i].dimension_relative())])
        X = self.affine_patch(indices)
        return X.multiplicity(X.ambient_space()(aff_pt))

#*******************************************************************
# Toric varieties
#*******************************************************************
class AlgebraicScheme_subscheme_toric(AlgebraicScheme_subscheme):
    r"""
    Construct an algebraic subscheme of a toric variety.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`~ToricVariety_field.subscheme` method of :class:`toric
        varieties
        <sage.schemes.toric.variety.ToricVariety_field>`.

    INPUT:

    - ``toric_variety`` -- ambient :class:`toric variety
      <ToricVariety_field>`;

    - ``polynomials`` -- single polynomial, list, or ideal of defining
      polynomials in the coordinate ring of ``toric_variety``.

    OUTPUT:

    - :class:`algebraic subscheme of a toric variety
      <AlgebraicScheme_subscheme_toric>`.

    TESTS::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1.inject_variables()
        Defining s, t, x, y
        sage: import sage.schemes.generic.algebraic_scheme as SCM
        sage: X = SCM.AlgebraicScheme_subscheme_toric(
        ....:       P1xP1, [x*s + y*t, x^3+y^3])
        sage: X
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          s*x + t*y,
          x^3 + y^3

    A better way to construct the same scheme as above::

        sage: P1xP1.subscheme([x*s + y*t, x^3+y^3])
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          s*x + t*y,
          x^3 + y^3
    """

    # Implementation note: if the toric variety is affine you should
    # construct instances of the derived class
    # AlgebraicScheme_subscheme_affine_toric instead.

    def __init__(self, toric_variety, polynomials):
        r"""
        See :class:`AlgebraicScheme_subscheme_toric` for documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: import sage.schemes.generic.algebraic_scheme as SCM
            sage: X = SCM.AlgebraicScheme_subscheme_toric(
            ....:       P1xP1, [x*s + y*t, x^3+y^3])
            sage: X
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 4 affine patches defined by:
              s*x + t*y,
              x^3 + y^3
        """
        # Just to make sure that keyword arguments will be passed correctly
        super(AlgebraicScheme_subscheme_toric, self).__init__(toric_variety,
                                                              polynomials)

    def _morphism(self, *args, **kwds):
        r"""
        Construct a morphism determined by action on points of ``self``.

        INPUT:

        - same as for
          :class:`~sage.schemes.toric.morphism.SchemeMorphism_polynomial_toric_variety`.

        OUTPUT:

        - :class:`~sage.schemes.toric.morphism.SchemeMorphism_polynomial_toric_variety`.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: P1 = P1xP1.subscheme(s - t)
            sage: H = P1.Hom(P1xP1)
            sage: H([s, s, x, y])
            Scheme morphism:
              From: Closed subscheme of 2-d CPR-Fano toric variety
              covered by 4 affine patches defined by:
              s - t
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [s : t : x : y] to
                    [s : s : x : y]

            sage: sbar, tbar, xbar, ybar = P1.coordinate_ring().gens()
            sage: P1._morphism(H, [sbar, sbar, xbar, ybar])
            Scheme morphism:
              From: Closed subscheme of 2-d CPR-Fano toric variety
              covered by 4 affine patches defined by:
              s - t
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [s : t : x : y] to
                    [t : t : x : y]
        """
        from sage.schemes.toric.morphism import SchemeMorphism_polynomial_toric_variety
        return SchemeMorphism_polynomial_toric_variety(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        r"""
        Construct a Hom-set for ``self``.

        INPUT:

        - same as for
          :class:`~sage.schemes.generic.homset.SchemeHomset_points_toric_field`.

        OUTPUT:

        :class:`~sage.schemes.toric.homset.SchemeHomset_points_subscheme_toric_field`.

        TESTS::

            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: quadric = P2.subscheme([x^2 + y^2 + z^2])
            sage: quadric._point_homset(Spec(QQ), quadric)
            Set of rational points of Closed subscheme of 2-d CPR-Fano
            toric variety covered by 3 affine patches defined by:
              x^2 + y^2 + z^2
            sage: type(quadric.point_set())
            <class 'sage.schemes.toric.homset.SchemeHomset_points_subscheme_toric_field_with_category'>
        """
        from sage.schemes.toric.homset import SchemeHomset_points_subscheme_toric_field
        return SchemeHomset_points_subscheme_toric_field(*args, **kwds)

    def fan(self):
        """
        Return the fan of the ambient space.

        OUTPUT:

        A fan.

        EXAMPLES::

            sage: P2.<x,y,z> = toric_varieties.P(2)
            sage: E = P2.subscheme([x^2+y^2+z^2])
            sage: E.fan()
            Rational polyhedral fan in 2-d lattice N
        """
        return self.ambient_space().fan()

    def affine_patch(self, i):
        r"""
        Return the ``i``-th affine patch of ``self`` as an affine
        toric algebraic scheme.

        INPUT:

        - ``i`` -- integer, index of a generating cone of the fan of the
          ambient space of ``self``.

        OUTPUT:

        - subscheme of an affine :class:`toric variety
          <sage.schemes.toric.variety.ToricVariety_field>`
          corresponding to the pull-back of ``self`` by the embedding
          morphism of the ``i``-th :meth:`affine patch of the ambient
          space
          <sage.schemes.toric.variety.ToricVariety_field.affine_patch>`
          of ``self``.

        The result is cached, so the ``i``-th patch is always the same object
        in memory.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: patch1 = P1xP1.affine_patch(1)
            sage: patch1.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [t : x] to
                    [1 : t : x : 1]
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: P1 = P1xP1.subscheme(x-y)
            sage: subpatch = P1.affine_patch(1)
            sage: subpatch
            Closed subscheme of 2-d affine toric variety defined by:
              x - 1
        """
        i = int(i)   # implicit type checking
        try:
            return self._affine_patches[i]
        except AttributeError:
            self._affine_patches = dict()
        except KeyError:
            pass
        ambient_patch = self.ambient_space().affine_patch(i)
        phi_p = ambient_patch.embedding_morphism().defining_polynomials()
        patch = ambient_patch.subscheme(
                            [p(phi_p) for p in self.defining_polynomials()])
        patch._embedding_morphism = patch.hom(phi_p, self, check=False)
        self._affine_patches[i] = patch
        return patch

    def affine_algebraic_patch(self, cone=None, names=None):
        r"""
        Return the affine patch corresponding to ``cone`` as an affine
        algebraic scheme.

        INPUT:

        - ``cone`` -- a :class:`Cone
          <sage.geometry.cone.ConvexRationalPolyhedralCone>` `\sigma`
          of the fan. It can be omitted for an affine toric variety,
          in which case the single generating cone is used.

        OUTPUT:

        An :class:`affine algebraic subscheme
        <sage.schemes.generic.algebraic_scheme.AlgebraicScheme_subscheme_affine>`
        corresponding to the patch `\mathop{Spec}(\sigma^\vee \cap M)`
        associated to the cone `\sigma`.

        See also :meth:`affine_patch`, which expresses the patches as
        subvarieties of affine toric varieties instead.

        REFERENCES:

        ..

            David A. Cox, "The Homogeneous Coordinate Ring of a Toric
            Variety", Lemma 2.2.
            http://www.arxiv.org/abs/alg-geom/9210008v2

        EXAMPLES::

            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: cone = P2.fan().generating_cone(0)
            sage: V = P2.subscheme(x^3+y^3+z^3)
            sage: V.affine_algebraic_patch(cone)
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              z0^3 + z1^3 + 1

            sage: cone = Cone([(0,1),(2,1)])
            sage: A2Z2.<x,y> = AffineToricVariety(cone)
            sage: A2Z2.affine_algebraic_patch()
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              -z0*z1 + z2^2
            sage: V = A2Z2.subscheme(x^2+y^2-1)
            sage: patch = V.affine_algebraic_patch();  patch
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              -z0*z1 + z2^2,
              z0 + z1 - 1
            sage: nbhd_patch = V.neighborhood([1,0]).affine_algebraic_patch();  nbhd_patch
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              -z0*z1 + z2^2,
              z0 + z1 - 1
            sage: nbhd_patch.embedding_center()
            (0, 1, 0)

        Here we got two defining equations. The first one describes
        the singularity of the ambient space and the second is the
        pull-back of `x^2+y^2-1` ::

            sage: lp = LatticePolytope([(1,0,0),(1,1,0),(1,1,1),(1,0,1),(-2,-1,-1)],
            ....:                      lattice=ToricLattice(3))
            sage: X.<x,y,u,v,t> = CPRFanoToricVariety(Delta_polar=lp)
            sage: Y = X.subscheme(x*v+y*u+t)
            sage: cone = Cone([(1,0,0),(1,1,0),(1,1,1),(1,0,1)])
            sage: Y.affine_algebraic_patch(cone)
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              z0*z2 - z1*z3,
              z1 + z3 + 1
        """
        from sage.modules.all import vector
        from sage.misc.all import prod
        ambient = self.ambient_space()
        fan = ambient.fan()
        if cone is None:
            assert ambient.is_affine()
            cone = fan.generating_cone(0)
        else:
            cone = fan.embed(cone)
        # R/I = C[sigma^dual cap M]
        R, I, dualcone = ambient._semigroup_ring(cone, names)

        # inhomogenize the Cox homogeneous polynomial with respect to the given cone
        inhomogenize = dict( (ambient.coordinate_ring().gen(i), 1)
                             for i in range(0,fan.nrays())
                             if not i in cone.ambient_ray_indices() )
        polynomials = [ p.subs(inhomogenize) for p in self.defining_polynomials() ]

        # map the monomial x^{D_m} to m, see reference.
        n_rho_matrix = cone.rays().matrix()
        def pullback_polynomial(p):
            result = R.zero()
            for coefficient, monomial in p:
                exponent = monomial.exponents()[0]
                exponent = [ exponent[i] for i in cone.ambient_ray_indices() ]
                exponent = vector(ZZ,exponent)
                m = n_rho_matrix.solve_right(exponent)
                assert all(x in ZZ for x in m), \
                    'The polynomial '+str(p)+' does not define a ZZ-divisor!'
                m_coeffs = dualcone.Hilbert_coefficients(m)
                result += coefficient * prod(R.gen(i)**m_coeffs[i]
                                             for i in range(0,R.ngens()))
            return result

        # construct the affine algebraic scheme to use as patch
        polynomials = [pullback_polynomial(_) for _ in polynomials]
        patch_cover = sage.schemes.affine.affine_space.AffineSpace(R)
        polynomials = list(I.gens()) + polynomials
        polynomials = [x for x in polynomials if not x.is_zero()]
        patch = patch_cover.subscheme(polynomials)

        # TODO: If the cone is not smooth, then the coordinate_ring()
        # of the affine toric variety is wrong; it should be the
        # G-invariant part. So we can't construct the embedding
        # morphism in that case.
        if cone.is_smooth():
            x = ambient.coordinate_ring().gens()
            phi = []
            for i in range(0,fan.nrays()):
                if i in cone.ambient_ray_indices():
                    phi.append(pullback_polynomial(x[i]))
                else:
                    phi.append(1)
            patch._embedding_morphism = patch.hom(phi, self)
        else:
            patch._embedding_morphism = (NotImplementedError,
               'I only know how to construct embedding morphisms for smooth patches')

        try:
            point = self.embedding_center()
        except AttributeError:
            return patch

        # it remains to find the preimage of point
        # map m to the monomial x^{D_m}, see reference.
        F = ambient.coordinate_ring().fraction_field()
        image = []
        for m in dualcone.Hilbert_basis():
            x_Dm = prod([ F.gen(i)**(m*n) for i,n in enumerate(fan.rays()) ])
            image.append(x_Dm)
        patch._embedding_center = tuple( f(list(point)) for f in image )
        return patch

    def _best_affine_patch(self, point):
        r"""
        Return the best affine patch of the ambient toric variety.

        INPUT:

        - ``point`` -- a point of the algebraic subscheme.

        OUTPUT:

        Integer. The index of the patch. See :meth:`affine_patch`.

        EXAMPLES::

            sage: P.<x,y,z>= toric_varieties.P2()
            sage: S = P.subscheme(x+2*y+3*z)
            sage: S._best_affine_patch(P.point([2,-3,0]))
            1
            sage: S._best_affine_patch([2,-3,0])
            1
        """
        # TODO: this method should pick a "best" patch in the sense
        # that it is numerically stable to dehomogenize, see the
        # corresponding method for projective varieties.
        point = list(point)
        zeros = set(i for i, coord in enumerate(point) if coord == 0)
        for cone_idx, cone in enumerate(self.ambient_space().fan().generating_cones()):
            if zeros.issubset(cone.ambient_ray_indices()):
                return cone_idx
        assert False, 'The point must not have been a point of the toric variety.'

    def neighborhood(self, point):
        r"""
        Return an toric algebraic scheme isomorphic to neighborhood of
        the ``point``.

        INPUT:

        - ``point`` -- a point of the toric algebraic scheme.

        OUTPUT:

        An affine toric algebraic scheme (polynomial equations in an
        affine toric variety) with fixed
        :meth:`~AlgebraicScheme.embedding_morphism` and
        :meth:`~AlgebraicScheme.embedding_center`.

        EXAMPLES::

            sage: P.<x,y,z>= toric_varieties.P2()
            sage: S = P.subscheme(x+2*y+3*z)
            sage: s = S.point([0,-3,2]); s
            [0 : -3 : 2]
            sage: patch = S.neighborhood(s); patch
            Closed subscheme of 2-d affine toric variety defined by:
              x + 2*y + 6
            sage: patch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of 2-d affine toric variety defined by:
              x + 2*y + 6
              To:   Closed subscheme of 2-d CPR-Fano toric variety covered by 3 affine patches defined by:
              x + 2*y + 3*z
              Defn: Defined on coordinates by sending [x : y] to
                    [-2*y - 6 : y : 2]
            sage: patch.embedding_center()
            [0 : -3]
            sage: patch.embedding_morphism()(patch.embedding_center())
            [0 : -3 : 2]

        A more complicated example::

            sage: dP6.<x0,x1,x2,x3,x4,x5> = toric_varieties.dP6()
            sage: twoP1 = dP6.subscheme(x0*x3)
            sage: patch = twoP1.neighborhood([0,1,2, 3,4,5]); patch
            Closed subscheme of 2-d affine toric variety defined by:
              3*x0
            sage: patch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of 2-d affine toric variety defined by:
              3*x0
              To:   Closed subscheme of 2-d CPR-Fano toric variety covered by 6 affine patches defined by:
              x0*x3
              Defn: Defined on coordinates by sending [x0 : x1] to
                    [0 : x1 : 2 : 3 : 4 : 5]
            sage: patch.embedding_center()
            [0 : 1]
            sage: patch.embedding_morphism()(patch.embedding_center())
            [0 : 1 : 2 : 3 : 4 : 5]
        """
        point = list(point)
        self._check_satisfies_equations(point)
        PP = self.ambient_space()
        n = PP.dimension()
        fan = PP.fan()
        cone_idx = self._best_affine_patch(point)
        cone = fan.generating_cone(cone_idx)

        patch_cover = PP.affine_patch(cone_idx)
        R = patch_cover.coordinate_ring()
        phi = []
        point_preimage = []
        for i in range(0,fan.nrays()):
            try:
                ray_index = cone.ambient_ray_indices().index(i)
                phi.append(R.gen(ray_index))
                point_preimage.append(point[i])
            except ValueError:
                phi.append(point[i])
        pullback_polys = [f(phi) for f in self.defining_polynomials()]
        patch = patch_cover.subscheme(pullback_polys)
        S = patch.coordinate_ring()
        phi_reduced = [S(t) for t in phi]

        patch._embedding_center = patch(point_preimage)
        patch._embedding_morphism = patch.hom(phi_reduced,self)
        return patch

    def dimension(self):
        """
        Return the dimension of ``self``.

        OUTPUT:

        Integer. If ``self`` is empty, `-1` is returned.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: P1 = P1xP1.subscheme(s-t)
            sage: P1.dimension()
            1
            sage: P1xP1.subscheme([s-t, (s-t)^2]).dimension()
            1
            sage: P1xP1.subscheme([s, t]).dimension()
            -1
        """
        if '_dimension' in self.__dict__:
            return self._dimension
        npatches = self.ambient_space().fan().ngenerating_cones()
        dims = [ self.affine_patch(i).dimension() for i in range(0,npatches) ]
        self._dimension = max(dims)
        return self._dimension

    def is_smooth(self, point=None):
        r"""
        Test whether the algebraic subscheme is smooth.

        INPUT:

        - ``point`` -- A point or ``None`` (default). The point to
          test smoothness at.

        OUTPUT:

        Boolean. If no point was specified, returns whether the
        algebraic subscheme is smooth everywhere. Otherwise,
        smoothness at the specified point is tested.

        EXAMPLES::

            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: cuspidal_curve = P2.subscheme([y^2*z-x^3])
            sage: cuspidal_curve
            Closed subscheme of 2-d CPR-Fano toric variety covered by 3 affine patches defined by:
              -x^3 + y^2*z
            sage: cuspidal_curve.is_smooth([1,1,1])
            True
            sage: cuspidal_curve.is_smooth([0,0,1])
            False
            sage: cuspidal_curve.is_smooth()
            False

        Any sufficiently generic cubic hypersurface is smooth::

            sage: P2.subscheme([y^2*z-x^3+z^3+1/10*x*y*z]).is_smooth()
            True

        A more complicated example::

            sage: dP6.<x0,x1,x2,x3,x4,x5> = toric_varieties.dP6()
            sage: disjointP1s = dP6.subscheme(x0*x3)
            sage: disjointP1s.is_smooth()
            True
            sage: intersectingP1s = dP6.subscheme(x0*x1)
            sage: intersectingP1s.is_smooth()
            False

        A smooth hypersurface in a compact singular toric variety::

            sage: lp = LatticePolytope([(1,0,0),(1,1,0),(1,1,1),(1,0,1),(-2,-1,-1)],
            ....:                      lattice=ToricLattice(3))
            sage: X.<x,y,u,v,t> = CPRFanoToricVariety(Delta_polar=lp)
            sage: Y = X.subscheme(x*v+y*u+t)
            sage: cone = Cone([(1,0,0),(1,1,0),(1,1,1),(1,0,1)])
            sage: Y.is_smooth()
            True
        """
        if not point is None:
            toric_patch = self.neighborhood(point)
            return toric_patch.is_smooth(toric_patch.embedding_center())

        # testing smoothness everywhere tends to be expensive
        if '_smooth' in self.__dict__:
            return self._smooth
        npatches = self.ambient_space().fan().ngenerating_cones()
        self._smooth = all(self.affine_patch(i).is_smooth() for i in range(0,npatches))
        return self._smooth


class AlgebraicScheme_subscheme_affine_toric(AlgebraicScheme_subscheme_toric):
    r"""
    Construct an algebraic subscheme of an affine toric variety.

    .. WARNING::

        You should not create objects of this class directly. The preferred
        method to construct such subschemes is to use
        :meth:`~ToricVariety_field.subscheme` method of
        :class:`toric varieties <ToricVariety_field>`.

    INPUT:

    - ``toric_variety`` -- ambient :class:`affine toric variety
      <ToricVariety_field>`;

    - ``polynomials`` -- single polynomial, list, or ideal of defining
      polynomials in the coordinate ring of ``toric_variety``.

    OUTPUT:

    A :class:`algebraic subscheme of an affine toric variety
    <AlgebraicScheme_subscheme_affine_toric>`.

    TESTS::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1.inject_variables()
        Defining s, t, x, y
        sage: import sage.schemes.generic.algebraic_scheme as SCM
        sage: X = SCM.AlgebraicScheme_subscheme_toric(
        ....:       P1xP1, [x*s + y*t, x^3+y^3])
        sage: X
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          s*x + t*y,
          x^3 + y^3

    A better way to construct the same scheme as above::

        sage: P1xP1.subscheme([x*s + y*t, x^3+y^3])
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          s*x + t*y,
          x^3 + y^3
    """

    def __init__(self, toric_variety, polynomials):
        r"""
        See :class:`AlgebraicScheme_subscheme_toric` for documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: import sage.schemes.generic.algebraic_scheme as SCM
            sage: X = SCM.AlgebraicScheme_subscheme_toric(
            ....:       P1xP1, [x*s + y*t, x^3+y^3])
            sage: X
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 4 affine patches defined by:
              s*x + t*y,
              x^3 + y^3
        """
        assert toric_variety.is_affine(), 'The toric variety must be affine!'
        # Just to make sure that keyword arguments will be passed correctly
        super(AlgebraicScheme_subscheme_affine_toric, self).__init__(toric_variety,
                                                                     polynomials)

    def dimension(self):
        """
        Return the dimension of ``self``.

        OUTPUT:

        - integer.

        EXAMPLES::

            sage: P1xP1.<s0,s1,t0,t1> = toric_varieties.P1xP1()
            sage: P1 = P1xP1.subscheme(s0-s1)
            sage: P1.dimension()
            1

        A more complicated example where the ambient toric variety is
        not smooth::

            sage: X.<x,y> = toric_varieties.A2_Z2()
            sage: X.is_smooth()
            False
            sage: Y = X.subscheme([x*y, x^2])
            sage: Y
            Closed subscheme of 2-d affine toric variety defined by:
              x*y,
              x^2
            sage: Y.dimension()
            1
        """
        if '_dimension' in self.__dict__:
            return self._dimension

        if self.ambient_space().is_smooth():
            self._dimension = self.defining_ideal().dimension()
        else:
            self._dimension = self.affine_algebraic_patch().dimension()
        return self._dimension

    def is_smooth(self, point=None):
        r"""
        Test whether the algebraic subscheme is smooth.

        INPUT:

        - ``point`` -- A point or ``None`` (default). The point to
          test smoothness at.

        OUTPUT:

        Boolean. If no point was specified, returns whether the
        algebraic subscheme is smooth everywhere. Otherwise,
        smoothness at the specified point is tested.

        EXAMPLES::

            sage: A2.<x,y> = toric_varieties.A2()
            sage: cuspidal_curve = A2.subscheme([y^2-x^3])
            sage: cuspidal_curve
            Closed subscheme of 2-d affine toric variety defined by:
              -x^3 + y^2
            sage: cuspidal_curve.is_smooth([1,1])
            True
            sage: cuspidal_curve.is_smooth([0,0])
            False
            sage: cuspidal_curve.is_smooth()
            False
            sage: circle = A2.subscheme(x^2+y^2-1)
            sage: circle.is_smooth([1,0])
            True
            sage: circle.is_smooth()
            True

        A more complicated example where the ambient toric variety is
        not smooth::

            sage: X.<x,y> = toric_varieties.A2_Z2()    # 2-d affine space mod Z/2
            sage: X.is_smooth()
            False
            sage: Y = X.subscheme([x*y, x^2])   # (twice the x=0 curve) mod Z/2
            sage: Y
            Closed subscheme of 2-d affine toric variety defined by:
              x*y,
              x^2
            sage: Y.dimension()   # Y is a Weil divisor but not Cartier
            1
            sage: Y.is_smooth()
            True
            sage: Y.is_smooth([0,0])
            True
        """
        if not point is None:
            self._check_satisfies_equations(point)
            if self.ambient_space().is_smooth():
                R = self.ambient_space().coordinate_ring()
                point_subs = dict(zip(R.gens(), point))
                Jac = self.Jacobian().subs(point_subs)
                return not Jac.is_zero()
            else:
                self._embedding_center = self.point(point)
                affine = self.affine_algebraic_patch()
                return affine.is_smooth(affine.embedding_center())

        # testing smoothness everywhere tends to be expensive
        if '_smooth' in self.__dict__:
            return self._smooth

        if self.ambient_space().is_smooth():
            sing_dim = self.Jacobian().dimension()
            self._smooth = (sing_dim == -1)
        else:
            self._smooth = self.affine_algebraic_patch().is_smooth()

        return self._smooth



