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
      -x1^2 + x0*x2,
      -x1*x2 + x0*x3,
      -x2^2 + x1*x3
      Defn: Defined on coordinates by sending (x0, x1, x2) to
            (1 : x0 : x1 : x2)


AUTHORS:

- David Kohel (2005): initial version.
- William Stein (2005): initial version.
- Andrey Novoseltsev (2010-05-17): subschemes of toric varieties.
- Volker Braun (2010-12-24): documentation of schemes and
  refactoring. Added coordinate neighborhoods and is_smooth()
"""

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



from sage.rings.all import ZZ

from sage.rings.ideal import is_Ideal
from sage.rings.rational_field import is_RationalField
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.finite_rings.constructor import is_FiniteField

from sage.misc.latex import latex
from sage.misc.misc import is_iterator
from sage.structure.all import Sequence
from sage.calculus.functions import jacobian

import sage.schemes.projective
import sage.schemes.affine
import ambient_space
import scheme



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

        sage: A, x = AffineSpace(10, QQ).objgens()
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
            raise TypeError, "A (=%s) must be an ambient space"
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
                raise(hom[0], hom[1])
            return hom
        ambient = self.ambient_space()
        return self.hom(ambient.coordinate_ring().gens(), ambient)

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
        raise AttributeError, 'This algebraic scheme does not have a designated point.'

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
            raise TypeError, "X must be a closed subscheme of an ambient space."
        if not isinstance(Y, AlgebraicScheme_subscheme):
            raise TypeError, "Y must be a closed subscheme of an ambient space."
        if X.ambient_space() != Y.ambient_space():
            raise ValueError, "X and Y must be embedded in the same ambient space."
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
        """
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
                raise TypeError, "Coordinates %s do not define a point on %s"%(v,self)
        for f in self.__Y.defining_polynomials():
            if f(coords) != 0:
                return True
        raise TypeError, "Coordinates %s do not define a point on %s"%(v,self)

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
                raise TypeError, "A positive bound (= %s) must be specified."%bound
            if not is_FiniteField(F):
                raise TypeError, "Argument F (= %s) must be a finite field."%F
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
                raise TypeError, "Coordinates %s do not define a point on %s"%(coords,self)
        try:
            return self.ambient_space()._check_satisfies_equations(coords)
        except TypeError:
            raise TypeError, "Coordinates %s do not define a point on %s"%(coords,self)

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
        """
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

        OUTPUT: an immutable sequence of irreducible subschemes of the
        ambient space of this scheme

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
        (see trac 6920)::

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
        """
        R = self.ambient_space().coordinate_ring()
        return jacobian(self.defining_polynomials(), R.gens())

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
            raise TypeError, "other (=%s) must be a closed algebraic subscheme of an ambient space"%other
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError, "other (=%s) must be in the same ambient space as self"%other
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
            raise TypeError, "other (=%s) must be a closed algebraic subscheme of an ambient space"%other
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError, "other (=%s) must be in the same ambient space as self"%other
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
                raise TypeError, \
                      "Argument other (=%s) must be a closed algebraic subscheme of an ambient space"%other
        if other.ambient_space() != A:
            raise ValueError, "other (=%s) must be in the same ambient space as self"%other
        return AlgebraicScheme_quasi(other, self)

    def rational_points(self, F=None, bound=0):
        """
        Return the rational points on the algebraic subscheme.

        EXAMPLES:

        One can enumerate points up to a given bound on a projective scheme
        over the rationals::

            sage: E = EllipticCurve('37a')
            sage: E.rational_points(bound=8)
            [(-1 : -1 : 1), (-1 : 0 : 1), (0 : -1 : 1), (0 : 0 : 1), (0 : 1 : 0), (1/4 : -5/8 : 1), (1/4 : -3/8 : 1), (1 : -1 : 1), (1 : 0 : 1), (2 : -3 : 1), (2 : 2 : 1)]

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
        if F == None:
            F = self.base_ring()
        X = self(F)
        if is_RationalField(F) or F == ZZ:
            if not bound > 0:
                raise TypeError, "A positive bound (= %s) must be specified."%bound
            try:
                return X.points(bound)
            except TypeError:
                raise TypeError, "Unable to enumerate points over %s."%F
        try:
            return X.points()
        except TypeError:
            raise TypeError, "Unable to enumerate points over %s."%F

    def change_ring(self,R):
        r"""
        Returns a new projective subscheme whose base ring is self coerced to R.

        EXAMPLES::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: X=P.subscheme([3*x^2-y^2])
            sage: H=Hom(X,X)
            sage: X.change_ring(GF(3))
            Closed subscheme of Projective Space of dimension 1 over Finite Field of size 3 defined by:
            -y^2
        """
        A=self.ambient_space().change_ring(R)
        I=self.defining_ideal().change_ring(A.coordinate_ring())
        return(A.subscheme(I))

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

    def projective_embedding(self, i=None, X=None):
        """
        Returns a morphism from this affine scheme into an ambient
        projective space of the same dimension.

        INPUT:

        -  ``i`` -- integer (default: dimension of self = last
           coordinate) determines which projective embedding to compute. The
           embedding is that which has a 1 in the i-th coordinate, numbered
           from 0.


        -  ``X`` -- (default: None) projective scheme, i.e., codomain of
           morphism; this is constructed if it is not given.

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
            raise ValueError, \
                  "Argument i (=%s) must be between 0 and %s, inclusive"%(i, n)
        try:
            return self.__projective_embedding[i]
        except AttributeError:
            self.__projective_embedding = {}
        except KeyError:
            pass
        if X is None:
            PP = sage.schemes.projective.projective_space.ProjectiveSpace(n, AA.base_ring())
            v = list(PP.gens())
            z = v.pop(i)
            v.append(z)
            polys = self.defining_polynomials()
            X = PP.subscheme([ f.homogenize()(v) for f in polys ])
        R = AA.coordinate_ring()
        v = list(R.gens())
        v.insert(i, R(1))
        phi = self.hom(v, X)
        self.__projective_embedding[i] = phi
        return phi

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

        OUPUT:

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

    def affine_patch(self, i):
        r"""
        Return the `i^{th}` affine patch of this projective scheme.
        This is the intersection with this `i^{th}` affine patch of
        its ambient space.

        INPUT:

        - ``i`` -- integer between 0 and dimension of self, inclusive.

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
        """
        i = int(i)   # implicit type checking
        PP = self.ambient_space()
        n = PP.dimension()
        if i < 0 or i > n:
            raise ValueError, "Argument i (= %s) must be between 0 and %s."%(i, n)
        try:
            return self.__affine_patches[i]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        AA = PP.affine_patch(i)
        phi = AA.projective_embedding()
        polys = self.defining_polynomials()
        xi = phi.defining_polynomials()
        U = AA.subscheme([ f(xi) for f in polys ])
        U._default_embedding_index = i
        phi = U.projective_embedding(i, self)
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
            abs_point = map(abs, point)
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

        OUTPUT

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
        # We really test the affine cone here; the origin is always a singular point:
        self._smooth = (sing_dim <= 0)
        return self._smooth


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

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan, "x s y t")
        sage: P1xP1.inject_variables()
        Defining x, s, y, t
        sage: import sage.schemes.generic.algebraic_scheme as SCM
        sage: X = SCM.AlgebraicScheme_subscheme_toric(
        ...         P1xP1, [x*s + y*t, x^3+y^3])
        sage: X
        Closed subscheme of 2-d toric variety
        covered by 4 affine patches defined by:
          x*s + y*t,
          x^3 + y^3

    A better way to construct the same scheme as above::

        sage: P1xP1.subscheme([x*s + y*t, x^3+y^3])
        Closed subscheme of 2-d toric variety
        covered by 4 affine patches defined by:
          x*s + y*t,
          x^3 + y^3
    """

    # Implementation note: if the toric variety is affine you should
    # construct instances of the derived class
    # AlgebraicScheme_subscheme_affine_toric instead.

    def __init__(self, toric_variety, polynomials):
        r"""
        See :class:`AlgebraicScheme_subscheme_toric` for documentation.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: import sage.schemes.generic.algebraic_scheme as SCM
            sage: X = SCM.AlgebraicScheme_subscheme_toric(
            ...         P1xP1, [x*s + y*t, x^3+y^3])
            sage: X
            Closed subscheme of 2-d toric variety
            covered by 4 affine patches defined by:
              x*s + y*t,
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

        OUPUT:

        - :class:`~sage.schemes.toric.morphism.SchemeMorphism_polynomial_toric_variety`.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.inject_variables()
            Defining z0, z1, z2, z3
            sage: P1 = P1xP1.subscheme(z0-z2)
            sage: H = P1.Hom(P1xP1)
            sage: H([z0,z1,z0,z3])
            Scheme morphism:
              From: Closed subscheme of 2-d toric variety
              covered by 4 affine patches defined by:
              z0 - z2
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [z0 : z1 : z2 : z3] to
                    [z2 : z1 : z2 : z3]

            sage: P1._morphism(H, [z0,z1,z0,z3])
            Scheme morphism:
              From: Closed subscheme of 2-d toric variety
              covered by 4 affine patches defined by:
              z0 - z2
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [z0 : z1 : z2 : z3] to
                    [z2 : z1 : z2 : z3]
        """
        from sage.schemes.toric.morphism import SchemeMorphism_polynomial_toric_variety
        return SchemeMorphism_polynomial_toric_variety(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        r"""
        Construct a Hom-set for ``self``.

        INPUT:

        - same as for
          :class:`~sage.schemes.generic.homset.SchemeHomset_points_toric_field`.

        OUPUT:

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

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: patch1 = P1xP1.affine_patch(1)
            sage: patch1.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [y : t] to
                    [1 : 1 : y : t]
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: P1 = P1xP1.subscheme(x-y)
            sage: subpatch = P1.affine_patch(1)
            sage: subpatch
            Closed subscheme of 2-d affine toric variety defined by:
              -y + 1
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

            sage: lp = LatticePolytope([(1,0,0),(1,1,0),(1,1,1),(1,0,1),(-2,-1,-1)])
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
        polynomials = map(pullback_polynomial, polynomials)
        patch_cover = sage.schemes.affine.affine_space.AffineSpace(R)
        polynomials = list(I.gens()) + polynomials
        polynomials = filter( lambda x:not x.is_zero(), polynomials)
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

        OUTPUT

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

        patch._embedding_center = patch(point_preimage)
        patch._embedding_morphism = patch.hom(phi,self)
        return patch

    def dimension(self):
        """
        Return the dimension of ``self``.

        OUTPUT:

        Integer. If ``self`` is empty, `-1` is returned.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.inject_variables()
            Defining z0, z1, z2, z3
            sage: P1 = P1xP1.subscheme(z0-z2)
            sage: P1.dimension()
            1
            sage: P1xP1.subscheme([z0-z2, (z0-z2)^2]).dimension()
            1
            sage: P1xP1.subscheme([z0,z2]).dimension()
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

            sage: lp = LatticePolytope(matrix([(1,0,0),(1,1,0),(1,1,1),(1,0,1),(-2,-1,-1)]).transpose())
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

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan, "x s y t")
        sage: P1xP1.inject_variables()
        Defining x, s, y, t
        sage: import sage.schemes.generic.algebraic_scheme as SCM
        sage: X = SCM.AlgebraicScheme_subscheme_toric(
        ...         P1xP1, [x*s + y*t, x^3+y^3])
        sage: X
        Closed subscheme of 2-d toric variety
        covered by 4 affine patches defined by:
          x*s + y*t,
          x^3 + y^3

    A better way to construct the same scheme as above::

        sage: P1xP1.subscheme([x*s + y*t, x^3+y^3])
        Closed subscheme of 2-d toric variety
        covered by 4 affine patches defined by:
          x*s + y*t,
          x^3 + y^3
    """

    def __init__(self, toric_variety, polynomials):
        r"""
        See :class:`AlgebraicScheme_subscheme_toric` for documentation.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: import sage.schemes.generic.algebraic_scheme as SCM
            sage: X = SCM.AlgebraicScheme_subscheme_toric(
            ...         P1xP1, [x*s + y*t, x^3+y^3])
            sage: X
            Closed subscheme of 2-d toric variety
            covered by 4 affine patches defined by:
              x*s + y*t,
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



