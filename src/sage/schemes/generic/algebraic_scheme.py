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
gluing relations, may be introduced.

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

    sage: V = A2.subscheme([x^2 + y^2 - 1]); V
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x^2 + y^2 - 1
    sage: V.dimension()                                                                 # needs sage.libs.singular
    1

Here is a more complicated example in a projective space::

    sage: P3 = ProjectiveSpace(3, QQ, 'x')
    sage: P3.inject_variables()
    Defining x0, x1, x2, x3
    sage: Q = matrix([[x0, x1, x2], [x1, x2, x3]]).minors(2); Q                         # needs sage.modules
    [-x1^2 + x0*x2, -x1*x2 + x0*x3, -x2^2 + x1*x3]
    sage: twisted_cubic = P3.subscheme(Q); twisted_cubic                                # needs sage.libs.singular sage.modules
    Closed subscheme of Projective Space of dimension 3 over Rational Field defined by:
      -x1^2 + x0*x2,
      -x1*x2 + x0*x3,
      -x2^2 + x1*x3
    sage: twisted_cubic.dimension()                                                     # needs sage.libs.singular sage.modules
    1

Note that there are 3 equations in the 3-dimensional ambient space,
yet the subscheme is 1-dimensional. One can show that it is not
possible to eliminate any of the equations, that is, the twisted cubic
is **not** a complete intersection of two polynomial equations.

Let us look at one affine patch, for example the one where `x_0=1` ::

    sage: patch = twisted_cubic.affine_patch(0); patch                                  # needs sage.libs.singular sage.modules
    Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
      -x1^2 + x2,
      -x1*x2 + x3,
      -x2^2 + x1*x3
    sage: patch.embedding_morphism()                                                    # needs sage.libs.singular sage.modules
    Scheme morphism:
      From: Closed subscheme of Affine Space of dimension 3 over Rational Field
            defined by: -x1^2 + x2, -x1*x2 + x3, -x2^2 + x1*x3
      To:   Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by: x1^2 - x0*x2, x1*x2 - x0*x3, x2^2 - x1*x3
      Defn: Defined on coordinates by sending (x1, x2, x3) to (1 : x1 : x2 : x3)


AUTHORS:

- David Kohel, William Stein (2005): initial version

- Andrey Novoseltsev (2010-05-17): subschemes of toric varieties

- Volker Braun (2010-12-24): documentation of schemes and refactoring; added
  coordinate neighborhoods and is_smooth()

- Ben Hutz (2014): subschemes of Cartesian products of projective space

- Ben Hutz (2017): split subschemes types into respective folders
"""

# ****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.abc

from sage.categories.number_fields import NumberFields

from sage.rings.ideal import Ideal_generic
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import RationalField
from sage.rings.finite_rings.finite_field_base import FiniteField

from sage.misc.latex import latex
from sage.misc.misc import is_iterator

from sage.structure.all import Sequence
from sage.structure.richcmp import richcmp, richcmp_method

from sage.arith.functions import lcm
from sage.arith.misc import gcd

import sage.schemes.affine
from . import ambient_space
from . import scheme

def is_AlgebraicScheme(x):
    """
    Test whether ``x`` is an algebraic scheme.

    INPUT:

    - ``x`` -- anything

    OUTPUT:

    boolean; whether ``x`` is an algebraic scheme, that is, a
    subscheme of an ambient space over a ring defined by polynomial
    equations.

    EXAMPLES::

        sage: A2 = AffineSpace(2, QQ, 'x, y')
        sage: A2.coordinate_ring().inject_variables()
        Defining x, y
        sage: V = A2.subscheme([x^2 + y^2]); V
        Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
          x^2 + y^2
        sage: from sage.schemes.generic.algebraic_scheme import is_AlgebraicScheme
        sage: is_AlgebraicScheme(V)
        doctest:warning...
        DeprecationWarning: The function is_AlgebraicScheme is deprecated; use 'isinstance(..., AlgebraicScheme)' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
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
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_AlgebraicScheme is deprecated; use 'isinstance(..., AlgebraicScheme)' instead.")
    return isinstance(x, AlgebraicScheme)


# ****************************************************************************
# A quick overview over the class hierarchy:
#
# class AlgebraicScheme(scheme.Scheme)
#    class AlgebraicScheme_subscheme
#       class AlgebraicScheme_subscheme_affine
#       class AlgebraicScheme_subscheme_projective
#       class AlgebraicScheme_subscheme_toric
#          class AlgebraicScheme_subscheme_affine_toric
#    class AlgebraicScheme_quasi
# ****************************************************************************


class AlgebraicScheme(scheme.Scheme):
    """
    An algebraic scheme presented as a subscheme in an ambient space.

    This is the base class for all algebraic schemes, that is, schemes
    defined by equations in affine, projective, or toric ambient
    spaces.
    """
    def __init__(self, A, category=None):
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
        if not isinstance(A, ambient_space.AmbientSpace):
            raise TypeError("A (=%s) must be an ambient space")
        self.__A = A
        self.__divisor_group = {}
        scheme.Scheme.__init__(self, A.base_scheme(), category=category)

    def _latex_(self):
        r"""
        Return a LaTeX representation of this algebraic scheme.

        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme
            sage: P = ProjectiveSpace(3, ZZ)
            sage: S = AlgebraicScheme(P); S
            Subscheme of Projective Space of dimension 3 over Integer Ring
            sage: S._latex_()
            '\\text{Subscheme of ${\\mathbf P}_{\\Bold{Z}}^{3}$}'
        """
        return r"\text{{Subscheme of ${}$}}".format(latex(self.__A))

    def is_projective(self):
        """
        Return ``True`` if ``self`` is presented as a subscheme of an ambient
        projective space.

        OUTPUT: boolean

        EXAMPLES::

            sage: PP.<x,y,z,w> = ProjectiveSpace(3, QQ)
            sage: f = x^3 + y^3 + z^3 + w^3
            sage: R = f.parent()
            sage: I = [f] + [f.derivative(zz) for zz in PP.gens()]
            sage: V = PP.subscheme(I)
            sage: V.is_projective()
            True
            sage: AA.<x,y,z,w> = AffineSpace(4, QQ)
            sage: V = AA.subscheme(I)
            sage: V.is_projective()
            False

        Note that toric varieties are implemented differently than
        projective spaces. This is why this method returns ``False``
        for toric varieties::

            sage: # needs sage.geometry.polyhedron sage.graphs
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
            sage: S = P.subscheme([x - y, x - z])
            sage: S.coordinate_ring()
            Quotient of Multivariate Polynomial Ring in x, y, z over Integer Ring
             by the ideal (x - y, x - z)
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
            sage: S = P.subscheme([x - y, x - z])
            sage: S.ambient_space() is P
            True
        """
        return self.__A

    def identity_morphism(self):
        """
        Return the identity morphism.

        OUTPUT: the identity morphism of the scheme ``self``

        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X.identity_morphism()
            Scheme endomorphism of Spectrum of Rational Field
              Defn: Identity map
        """
        from sage.schemes.generic.morphism import SchemeMorphism_polynomial_id
        return SchemeMorphism_polynomial_id(self)

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

        * A :exc:`NotImplementedError` is raised if the construction of
          the embedding morphism is not implemented yet.

        EXAMPLES::

            sage: A2.<x,y> = AffineSpace(QQ, 2)
            sage: C = A2.subscheme(x^2 + y^2 - 1)
            sage: C.embedding_morphism()                                                # needs sage.libs.singular
              Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: x^2 + y^2 - 1
              To:   Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to (x, y)

            sage: # needs sage.geometry.polyhedron sage.graphs sage.libs.singular
            sage: P1xP1.<x,y,u,v> = toric_varieties.P1xP1()
            sage: P1 = P1xP1.subscheme(x - y)
            sage: P1.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of 2-d CPR-Fano toric variety covered
                    by 4 affine patches defined by: x - y
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [x : y : u : v] to [y : y : u : v]

        So far, the embedding was just in the own ambient space. Now a
        bit more interesting examples::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme((x^2-y^2)*z)
            sage: p = (1,1,0)
            sage: nbhd = X.neighborhood(p); nbhd
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -y^2*z - 2*y*z

        Note that `p=(1,1,0)` is a singular point of `X`. So the
        neighborhood of `p` is not just affine space. The
        :meth:`neighborhood` method returns a presentation of
        the neighborhood as a subscheme of an auxiliary 2-dimensional
        affine space::

            sage: nbhd.ambient_space()
            Affine Space of dimension 2 over Rational Field

        But its :meth:`embedding_morphism` is not into this auxiliary
        affine space, but the original subscheme `X`::

            sage: nbhd.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: -y^2*z - 2*y*z
              To:   Closed subscheme of Projective Space of dimension 2 over Rational Field
                    defined by: x^2*z - y^2*z
              Defn: Defined on coordinates by sending (y, z) to (1 : y + 1 : z)

        A couple more examples::

            sage: # needs sage.geometry.polyhedron sage.graphs sage.libs.singular
            sage: patch1 = P1xP1.affine_patch(1); patch1
            2-d affine toric variety
            sage: patch1.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [y : u] to [1 : y : u : 1]
            sage: subpatch = P1.affine_patch(1); subpatch
            Closed subscheme of 2-d affine toric variety defined by:
              -y + 1
            sage: subpatch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of 2-d affine toric variety defined by: -y + 1
              To:   Closed subscheme of 2-d CPR-Fano toric variety covered
                    by 4 affine patches defined by: x - y
              Defn: Defined on coordinates by sending [y : u] to [1 : y : u : 1]
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

        A point of ``self``. This raises :exc:`AttributeError` if there
        is no distinguished point, depending on how ``self`` was constructed.

        EXAMPLES::

            sage: P3.<w,x,y,z> = ProjectiveSpace(QQ, 3)
            sage: X = P3.subscheme( (w^2-x^2)*(y^2-z^2) )
            sage: p = [1,-1,3,4]
            sage: nbhd = X.neighborhood(p); nbhd
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              w^2*y^2 - x^2*y^2 + 6*w^2*y - 6*x^2*y + 2*w*y^2 +
              2*x*y^2 - 7*w^2 + 7*x^2 + 12*w*y + 12*x*y - 14*w - 14*x
            sage: nbhd.embedding_center()
            (0, 0, 0)
            sage: nbhd.embedding_morphism()(nbhd.embedding_center())
            (1/4 : -1/4 : 3/4 : 1)
            sage: nbhd.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Rational Field
                    defined by: w^2*y^2 - x^2*y^2 + 6*w^2*y - 6*x^2*y + 2*w*y^2 + 2*x*y^2
                                - 7*w^2 + 7*x^2 + 12*w*y + 12*x*y - 14*w - 14*x
              To:   Closed subscheme of Projective Space of dimension 3 over Rational Field
                    defined by: w^2*y^2 - x^2*y^2 - w^2*z^2 + x^2*z^2
              Defn: Defined on coordinates by sending (w, x, y) to
                    (w + 1 : x - 1 : y + 3 : 4)
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
            sage: S = P.subscheme([x - y, x - z])
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
        return "Subscheme of %s" % self.__A

    def _homset(self, *args, **kwds):
        """
        Construct the Hom-set.

        INPUT:

        Same as :class:`sage.schemes.generic.homset.SchemeHomset_generic`.

        OUTPUT: the Hom-set of the ambient space

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron sage.graphs
            sage: P1.<x,y> = toric_varieties.P1()
            sage: type(P1.Hom(P1))
            <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>
            sage: X = P1.subscheme(x - y)
            sage: type(X.Hom(X))
            <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>

        ::

            sage: # needs sage.geometry.polyhedron sage.graphs
            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: P1xP1._homset(P1xP1, P1)
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


class AlgebraicScheme_quasi(AlgebraicScheme):
    """
    The quasi-affine or quasi-projective scheme `X - Y`, where `X` and `Y`
    are both closed subschemes of a common ambient affine or projective
    space.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`complement` method of algebraic schemes.

    OUTPUT: an instance of :class:`AlgebraicScheme_quasi`

    EXAMPLES::

        sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
        sage: S = P.subscheme([])
        sage: T = P.subscheme([x - y])
        sage: T.complement(S)
        Quasi-projective subscheme X - Y
         of Projective Space of dimension 2 over Integer Ring,
         where X is defined by:
           (no polynomials)
         and Y is defined by:
           x - y
    """

    def __init__(self, X, Y):
        """
        The constructor.

        INPUT:

        - ``X``, ``Y`` -- two subschemes of the same ambient space

        TESTS::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x - y])
            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_quasi
            sage: AlgebraicScheme_quasi(S, T)
            Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring,
             where X is defined by:
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
            sage: T = P.subscheme([x - y])
            sage: U = AlgebraicScheme_quasi(S, T); U
            Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring,
             where X is defined by:
               (no polynomials)
             and Y is defined by:
               x - y
            sage: U._latex_()
            '\\text{Quasi-projective subscheme }
             (X\\setminus Y)\\subset {\\mathbf P}_{\\Bold{Z}}^{2},\\text{ where }
             X \\text{ is defined by }\\text{no polynomials},\\text{ and }
             Y \\text{ is defined by } x - y.'
        """
        if isinstance(self.ambient_space(), sage.schemes.affine.affine_space.AffineSpace_generic):
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
            sage: T = P.subscheme([x - y])
            sage: U = AlgebraicScheme_quasi(S, T); U
            Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring,
             where X is defined by:
               (no polynomials)
             and Y is defined by:
               x - y
            sage: U._repr_()
            'Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring, where X is defined by:\n  (no polynomials)\nand Y is defined by:\n  x - y'
        """
        if isinstance(self.ambient_space(), sage.schemes.affine.affine_space.AffineSpace_generic):
            t = "affine"
        else:
            t = "projective"
        return ("Quasi-%s subscheme X - Y of %s, where X is defined by:\n%s\n"
                "and Y is defined by:\n%s"
                % (t, self.ambient_space(), str(self.__X).split("\n", 1)[1],
                   str(self.__Y).split("\n", 1)[1]))

    def X(self):
        """
        Return the scheme `X` such that ``self`` is represented as `X - Y`.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x - y])
            sage: U = T.complement(S)
            sage: U.X() is S
            True
        """
        return self.__X

    def Y(self):
        """
        Return the scheme `Y` such that ``self`` is represented as `X - Y`.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x - y])
            sage: U = T.complement(S)
            sage: U.Y() is T
            True
        """
        return self.__Y

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme, or
        raise a :exc:`TypeError`.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x - y])
            sage: U = T.complement(S)
            sage: U._check_satisfies_equations([1, 2, 0])
            True
            sage: U._check_satisfies_equations([1, 1, 0])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 1, 0] do not define a point on
            Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring,
             where X is defined by:
               (no polynomials)
             and Y is defined by:
               x - y

            sage: U._check_satisfies_equations([1, 4])
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match number of variables in parent

            sage: A.<x, y> = AffineSpace(2, GF(7))
            sage: S = A.subscheme([x^2 - y])
            sage: T = A.subscheme([x - y])
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

        The bug reported at :issue:`12211` has been fixed::

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
                raise TypeError("Coordinates %s do not define a point on %s" % (v,self))
        for f in self.__Y.defining_polynomials():
            if f(coords) != 0:
                return True
        raise TypeError("Coordinates %s do not define a point on %s" % (v,self))

    def rational_points(self, **kwds):
        """
        Return the set of rational points on this algebraic scheme
        over the field `F`.

        INPUT: keyword arguments:

        - ``bound`` -- integer (default: 0); the bound for the coordinates for
          subschemes with dimension at least 1

        - ``F`` -- field (default: base ring); the field to compute
          the rational points over

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, GF(7))
            sage: S = A.subscheme([x^2 - y])
            sage: T = A.subscheme([x - y])
            sage: U = T.complement(S)
            sage: U.rational_points()
            [(2, 4), (3, 2), (4, 2), (5, 4), (6, 1)]
            sage: U.rational_points(F=GF(7^2, 'b'))                                     # needs sage.rings.finite_rings
            [(2, 4), (3, 2), (4, 2), (5, 4), (6, 1), (b, b + 4), (b + 1, 3*b + 5),
             (b + 2, 5*b + 1), (b + 3, 6), (b + 4, 2*b + 6), (b + 5, 4*b + 1),
             (b + 6, 6*b + 5), (2*b, 4*b + 2), (2*b + 1, b + 3), (2*b + 2, 5*b + 6),
             (2*b + 3, 2*b + 4), (2*b + 4, 6*b + 4), (2*b + 5, 3*b + 6), (2*b + 6, 3),
             (3*b, 2*b + 1), (3*b + 1, b + 2), (3*b + 2, 5), (3*b + 3, 6*b + 3),
             (3*b + 4, 5*b + 3), (3*b + 5, 4*b + 5), (3*b + 6, 3*b + 2),
             (4*b, 2*b + 1), (4*b + 1, 3*b + 2), (4*b + 2, 4*b + 5),
             (4*b + 3, 5*b + 3), (4*b + 4, 6*b + 3), (4*b + 5, 5), (4*b + 6, b + 2),
             (5*b, 4*b + 2), (5*b + 1, 3), (5*b + 2, 3*b + 6), (5*b + 3, 6*b + 4),
             (5*b + 4, 2*b + 4), (5*b + 5, 5*b + 6), (5*b + 6, b + 3), (6*b, b + 4),
             (6*b + 1, 6*b + 5), (6*b + 2, 4*b + 1), (6*b + 3, 2*b + 6), (6*b + 4, 6),
             (6*b + 5, 5*b + 1), (6*b + 6, 3*b + 5)]
        """
        F = kwds.get('F', None)
        bound = kwds.get('bound', 0)
        if F is None:
            F = self.base_ring()

        if bound == 0:
            if isinstance(F, RationalField):
                raise TypeError("A positive bound (= %s) must be specified." % bound)
            if not isinstance(F, FiniteField):
                raise TypeError("Argument F (= %s) must be a finite field." % F)
        pts = []
        for P in self.ambient_space().rational_points(F):
            try:
                if self._check_satisfies_equations(list(P)):
                    pts.append(P)
            except TypeError:
                pass
        pts.sort()
        return pts


@richcmp_method
class AlgebraicScheme_subscheme(AlgebraicScheme):
    """
    An algebraic scheme presented as a closed subscheme is defined by
    explicit polynomial equations. This is as opposed to a general
    scheme, which could, e.g., be the Neron model of some object, and
    for which we do not want to give explicit equations.

    INPUT:

    - ``A`` -- ambient space (e.g. affine or projective `n`-space)

    - ``polynomials`` -- single polynomial, ideal or iterable of defining
      polynomials; in any case polynomials must belong to the coordinate
      ring of the ambient space and define valid polynomial functions (e.g.
      they should be homogeneous in the case of a projective space)

    OUTPUT: algebraic scheme

    EXAMPLES::

        sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
        sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
        sage: P.subscheme([x^2 - y*z])
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x^2 - y*z
        sage: AlgebraicScheme_subscheme(P, [x^2 - y*z])
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x^2 - y*z
    """

    def __init__(self, A, polynomials, category=None):
        """
        See ``AlgebraicScheme_subscheme`` for documentation.

        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: P.subscheme([x^2 - y*z])
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
            sage: AlgebraicScheme_subscheme(P, [x^2 - y*z])
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
        """
        from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence_generic

        AlgebraicScheme.__init__(self, A, category=category)
        self._base_ring = A.base_ring()
        R = A.coordinate_ring()
        if isinstance(polynomials, Ideal_generic):
            I = polynomials
            polynomials = I.gens()
            if I.ring() is R: # Otherwise we will recompute I later after
                self.__I = I  # converting generators to the correct ring
        if isinstance(polynomials, tuple) or isinstance(polynomials, PolynomialSequence_generic) or is_iterator(polynomials):
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
        raise a :exc:`TypeError`.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: S = P.subscheme([x^2 - y*z])
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
                raise TypeError("Coordinates %s do not define a point on %s" % (coords,self))
        try:
            return self.ambient_space()._check_satisfies_equations(coords)
        except TypeError:
            raise TypeError("Coordinates %s do not define a point on %s" % (coords,self))

    def base_extend(self, R):
        """
        Return the base change to the ring `R` of this scheme.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(11))
            sage: S = P.subscheme([x^2 - y*z])
            sage: S.base_extend(GF(11^2, 'b'))                                          # needs sage.rings.finite_rings
            Closed subscheme of Projective Space of dimension 2
             over Finite Field in b of size 11^2
             defined by: x^2 - y*z
            sage: S.base_extend(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: no natural map from the base ring (=Finite Field of size 11)
            to R (=Integer Ring)!
        """
        A = self.ambient_space().base_extend(R)
        return A.subscheme(self.__polys)

    def __richcmp__(self, other, op):
        """
        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, QQ)
            sage: X = A.subscheme([x*y, z])
            sage: X == A.subscheme([z, x*y])
            True
            sage: X == A.subscheme([x*y, z^2])                                          # needs sage.libs.singular
            False
            sage: B.<u, v, t> = AffineSpace(3, QQ)
            sage: X == B.subscheme([u*v, t])
            False
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            return NotImplemented
        A = self.ambient_space()
        if other.ambient_space() != A:
            return NotImplemented
        return richcmp(self.defining_ideal(), other.defining_ideal(), op)

    def _latex_(self):
        """
        Return a LaTeX representation of this scheme.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(11))
            sage: S = P.subscheme([x^2 - y*z]); S
            Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:
              x^2 - y*z
            sage: S._latex_()
            '\\text{Closed subscheme of } {\\mathbf P}_{\\Bold{F}_{11}}^{2} \\text{ defined by } x^{2} - y z'
            sage: S = P.subscheme([x^2 - y*z, x^5]); S
            Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:
              x^2 - y*z,
              x^5
            sage: S._latex_()
            '\\text{Closed subscheme of } {\\mathbf P}_{\\Bold{F}_{11}}^{2} \\text{ defined by } x^{2} - y z, x^{5}'
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
            sage: S = P.subscheme([x^2 - y*z]); S
            Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:
              x^2 - y*z
            sage: S._repr_()
            'Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:\n  x^2 - y*z'
            sage: S = P.subscheme([x^2 - y*z, x^5]); S
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
            sage: S = P.subscheme([x^2 - y*z, x^3 + z^3])
            sage: S.defining_polynomials()
            (x^2 - y*z, x^3 + z^3)
        """
        return self.__polys

    def normalize_defining_polynomials(self):
        r"""
        Function to normalize the coefficients of defining polynomials
        of given subscheme.

        Normalization as in removing denominator from all the coefficients,
        and then removing any common factor between the coefficients.
        It takes LCM of denominators and then removes common factor among
        coefficients, if any.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: S = A.subscheme([2*x^2 + 4*x*y, 1/8*x + 1/3*y])
            sage: S.normalize_defining_polynomials()
            sage: S.defining_polynomials()
            (x^2 + 2*x*y, 3*x + 8*y)
        """
        BR = self.base_ring()
        if (BR == ZZ
                or isinstance(BR, (sage.rings.abc.AlgebraicField, sage.rings.abc.Order))
                or BR in NumberFields()):
            normalized_polys = []
            initial_polys = list(self.__polys)

            for P in initial_polys:
                # stores value which need to be mutliplied to make all coefficient integers
                mult = lcm([c.denominator() for c in P.coefficients()])
                P = mult*P
                # stores the common factor from all coefficients
                div = gcd(list(P.coefficients()))
                poly_ring = P.parent() # need to coerce, since division might change base ring
                P = poly_ring((BR.one()/div)*P)
                normalized_polys.append(P)

            self.__polys = tuple(normalized_polys)

        else:
                raise NotImplementedError("currently normalization is implemented "
                    "only for QQbar, number fields and number field orders")

    def defining_ideal(self):
        """
        Return the ideal that defines this scheme as a subscheme
        of its ambient space.

        OUTPUT: an ideal in the coordinate ring of the ambient space

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x^2 - y*z, x^3 + z^3])
            sage: S.defining_ideal()
            Ideal (x^2 - y*z, x^3 + z^3) of Multivariate Polynomial Ring in x, y, z
             over Integer Ring
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

        OUTPUT: integer

        EXAMPLES::

            sage: PP.<x,y,z,w,v> = ProjectiveSpace(4, QQ)
            sage: V = PP.subscheme(x*y)
            sage: V.codimension()                                                       # needs sage.libs.singular
            1
            sage: V.dimension()                                                         # needs sage.libs.singular
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

            sage: PP.<x,y,z,w,v> = ProjectiveSpace(4, QQ)
            sage: V = PP.subscheme((x^2 - y^2 - z^2) * (w^5 - 2*v^2*z^3) * w * (v^3 - x^2*z))
            sage: V.irreducible_components()                                            # needs sage.libs.singular
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

        We verify that the irrelevant ideal is not accidentally returned
        (see :issue:`6920`)::

            sage: PP.<x,y,z,w> = ProjectiveSpace(3, QQ)
            sage: f = x^3 + y^3 + z^3 + w^3
            sage: R = f.parent()
            sage: I = [f] + [f.derivative(zz) for zz in PP.gens()]
            sage: V = PP.subscheme(I)
            sage: V.irreducible_components()                                            # needs sage.libs.singular
            [
            <BLANKLINE>
            ]

        The same polynomial as above defines a scheme with a
        nontrivial irreducible component in affine space (instead of
        the empty scheme as above)::

            sage: AA.<x,y,z,w> = AffineSpace(4, QQ)
            sage: V = AA.subscheme(I)
            sage: V.irreducible_components()                                            # needs sage.libs.singular
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
        C.sort(key=lambda scheme: scheme.defining_ideal().gens())
        C.set_immutable()
        self.__irreducible_components = C
        return C

    def is_irreducible(self):
        r"""
        Return whether this subscheme is or is not irreducible.

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QuadraticField(-3)
            sage: P.<x,y,z,w,t,u> = ProjectiveSpace(K, 5)
            sage: X = P.subscheme([x*y - z^2 - K.0*t^2, t*w*x + y*z^2 - u^3])
            sage: X.is_irreducible()                                                    # needs sage.libs.singular
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme([(y + x - z)^2])
            sage: X.is_irreducible()                                                    # needs sage.libs.singular
            False

        ::

            sage: A.<x,y,z,w> = AffineSpace(GF(17), 4)
            sage: X = A.subscheme([
            ....:         x*y*z^2 - x*y*z*w - z*w^2 + w^3,
            ....:         x^3*y*z*w - x*y^3*z - x^2*y*z*w - x^2*w^3 + y^2*w^2 + x*w^3
            ....:     ])
            sage: X.is_irreducible()                                                    # needs sage.libs.singular
            False
        """
        return self.defining_ideal().is_prime()

    def Jacobian_matrix(self):
        r"""
        Return the matrix `\frac{\partial f_i}{\partial x_j}` of
        (formal) partial derivatives.

        OUTPUT: a matrix of polynomials

        EXAMPLES::

            sage: P3.<w,x,y,z> = ProjectiveSpace(3, QQ)
            sage: twisted_cubic = P3.subscheme(matrix([[w, x, y],                       # needs sage.libs.singular
            ....:                                      [x, y, z]]).minors(2))
            sage: twisted_cubic.Jacobian_matrix()                                       # needs sage.libs.singular
            [   y -2*x    w    0]
            [   z   -y   -x    w]
            [   0    z -2*y    x]

        This example addresses issue :issue:`20512`::

            sage: X = P3.subscheme([])
            sage: X.Jacobian_matrix().base_ring() == P3.coordinate_ring()               # needs sage.libs.singular
            True
        """
        R = self.ambient_space().coordinate_ring()
        l = self.defining_polynomials()
        if len(l) == 0:
            return sage.matrix.constructor.Matrix(R, 0)

        from sage.calculus.functions import jacobian

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
          the defining equations is [Laz2004]_, p. 181.

        OUTPUT: an ideal in the coordinate ring of the ambient space

        EXAMPLES::

            sage: P3.<w,x,y,z> = ProjectiveSpace(3, QQ)
            sage: twisted_cubic = P3.subscheme(matrix([[w, x, y],                       # needs sage.libs.singular
            ....:                                      [x, y, z]]).minors(2))
            sage: twisted_cubic.Jacobian()                                              # needs sage.libs.singular
            Ideal (-x^2 + w*y, -x*y + w*z, -y^2 + x*z, x*z, -2*w*z, w*y, 3*w*y,
                   -2*w*x, w^2, y*z, -2*x*z, w*z, 3*w*z, -2*w*y, w*x, z^2, -2*y*z,
                   x*z, 3*x*z, -2*w*z, w*y)
             of Multivariate Polynomial Ring in w, x, y, z over Rational Field
            sage: twisted_cubic.defining_ideal()                                        # needs sage.libs.singular
            Ideal (-x^2 + w*y, -x*y + w*z, -y^2 + x*z)
             of Multivariate Polynomial Ring in w, x, y, z over Rational Field

        This example addresses issue :issue:`20512`::

            sage: X = P3.subscheme([])
            sage: X.Jacobian() == P3.coordinate_ring().unit_ideal()                     # needs sage.libs.singular
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
            sage: X.dimension()                                                         # needs sage.libs.singular
            1

        Then we compute the corresponding reduced scheme::

            sage: Y = X.reduce(); Y                                                     # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^2 - x*y - x + y

        Finally, we verify that the reduced scheme `Y` is the union
        of those two lines::

            sage: # needs sage.libs.singular
            sage: L1 = A.subscheme([x - 1]); L1
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x - 1
            sage: L2 = A.subscheme([x - y]); L2
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
        Return the scheme-theoretic union of ``self`` and ``other`` in their common
        ambient space.

        EXAMPLES: We construct the union of a line and a tripled-point on
        the line.

        ::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: I = ideal([x, y])^3
            sage: P = A.subscheme(I)
            sage: L = A.subscheme([y - 1])
            sage: S = L.union(P); S                                                     # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              y^4 - y^3,
              x*y^3 - x*y^2,
              x^2*y^2 - x^2*y,
              x^3*y - x^3
            sage: S.dimension()                                                         # needs sage.libs.singular
            1
            sage: S.reduce()                                                            # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              y^2 - y,
              x*y - x

        We can also use the notation "+" for the union::

            sage: A.subscheme([x]) + A.subscheme([y^2 - (x^3+1)])                       # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^4 - x*y^2 + x

        Saving and loading::

            sage: loads(S.dumps()) == S                                                 # needs sage.libs.singular
            True
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            raise TypeError("other (=%s) must be a closed algebraic subscheme of an ambient space" % other)
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError("other (=%s) must be in the same ambient space as self" % other)
        return A.subscheme(self.defining_ideal().intersection(other.defining_ideal()))

    def __pow__(self, m):
        """
        Return the Cartesian power of this space.

        INPUT:

        - ``m`` -- integer

        OUTPUT: subscheme of product of ambient spaces

        EXAMPLES::

            sage: P2.<y0,y1,y2> = ProjectiveSpace(ZZ, 2)
            sage: Z = P2.subscheme([y0^2 - y1*y2, y2])
            sage: Z**3
            Closed subscheme of Product of projective spaces P^2 x P^2 x P^2 over
             Integer Ring defined by:
              x0^2 - x1*x2,
              x2,
              x3^2 - x4*x5,
              x5,
              x6^2 - x7*x8,
              x8

            sage: A2.<x,y> = AffineSpace(QQ, 2)
            sage: V = A2.subscheme([x^2 - y, x - 1])
            sage: V**4
            Closed subscheme of Affine Space of dimension 8 over Rational Field
             defined by:
              x0^2 - x1,
              x0 - 1,
              x2^2 - x3,
              x2 - 1,
              x4^2 - x5,
              x4 - 1,
              x6^2 - x7,
              x6 - 1

            sage: T.<x0,x1,x2,x3,x4,x5> = ProductProjectiveSpaces([2,2], ZZ)
            sage: X = T.subscheme([x0*x4 - x1*x3])
            sage: X^2
            Closed subscheme of Product of projective spaces P^2 x P^2 x P^2 x P^2
             over Integer Ring defined by:
              -x1*x3 + x0*x4,
              -x7*x9 + x6*x10

            sage: E = EllipticCurve([0,0,0,0,1])                                        # needs sage.schemes
            sage: E^2                                                                   # needs sage.schemes
            Closed subscheme of Product of projective spaces P^2 x P^2
             over Rational Field defined by:
              -x0^3 + x1^2*x2 - x2^3,
              -x3^3 + x4^2*x5 - x5^3
        """
        AS = self.ambient_space().__pow__(m)
        CR = AS.coordinate_ring()
        n = self.ambient_space().coordinate_ring().ngens()

        polys = []
        for i in range(m):
            phi = self.ambient_space().coordinate_ring().hom(list(CR.gens()[n*i : n*(i+1)]), CR)
            polys.extend([phi(t) for t in self.defining_polynomials()])
        return AS.subscheme(polys)

    def __mul__(self, right):
        r"""
        Create the product of subschemes.

        INPUT:

        - ``right`` -- a subscheme of similar type

        OUTPUT: a subscheme of a the product of the ambient spaces

        EXAMPLES::

            sage: S = ProductProjectiveSpaces([1,2,1], ZZ, 't')
            sage: T = ProductProjectiveSpaces([2,2], ZZ, 'x')
            sage: T.inject_variables()
            Defining x0, x1, x2, x3, x4, x5
            sage: X = T.subscheme([x0*x4 - x1*x3])
            sage: X*S
            Closed subscheme of
             Product of projective spaces P^2 x P^2 x P^1 x P^2 x P^1 over Integer Ring
             defined by:
              -x1*x3 + x0*x4

        ::

            sage: S = ProjectiveSpace(ZZ, 2, 't')
            sage: T.<x0,x1,x2,x3> = ProjectiveSpace(ZZ, 3)
            sage: X = T.subscheme([x0*x2 - x1*x3])
            sage: X*S
            Closed subscheme of Product of projective spaces P^3 x P^2 over Integer Ring
             defined by:
              x0*x2 - x1*x3

        ::

            sage: A2 = AffineSpace(ZZ, 2, 't')
            sage: A3.<x0,x1,x2> = AffineSpace(ZZ, 3)
            sage: X = A3.subscheme([x0*x2 - x1])
            sage: X*A2
            Closed subscheme of Affine Space of dimension 5 over Integer Ring
             defined by:
              x0*x2 - x1

        ::

            sage: T.<x0,x1,x2,x3,x4,x5> = ProductProjectiveSpaces([2,2], ZZ)
            sage: X = T.subscheme([x0*x4 - x1*x3])
            sage: X*X
            Closed subscheme of Product of projective spaces P^2 x P^2 x P^2 x P^2
             over Integer Ring defined by:
              -x1*x3 + x0*x4,
              -x7*x9 + x6*x10

        ::

            sage: P1.<z0,z1> = ProjectiveSpace(ZZ, 1)
            sage: Y = P1.subscheme([z0 - z1])
            sage: T.<x0,x1,x2,x3,x4,x5> = ProductProjectiveSpaces([2,2], ZZ)
            sage: X = T.subscheme([x0*x4 - x1*x3])
            sage: X*Y
            Closed subscheme of Product of projective spaces P^2 x P^2 x P^1
             over Integer Ring defined by:
              -x1*x3 + x0*x4,
              z0 - z1

        ::

            sage: A3.<x0,x1,x2> = AffineSpace(ZZ, 3)
            sage: X = A3.subscheme([x0*x2 - x1])
            sage: P1.<u,v> = ProjectiveSpace(ZZ, 1)
            sage: Y = P1.subscheme([u - v])
            sage: X*Y
            Traceback (most recent call last):
            ...
            TypeError: Projective Space of dimension 1 over Integer Ring must be an affine space or affine subscheme
            sage: Y*X
            Traceback (most recent call last):
            ...
            TypeError: Affine Space of dimension 3 over Integer Ring must be a projective space,
            product of projective spaces, or subscheme
            sage: PP.<a,b,c,d> = ProductProjectiveSpaces(ZZ, [1,1])
            sage: Z = PP.subscheme([a*d - b*c])
            sage: X*Z
            Traceback (most recent call last):
            ...
            TypeError: Product of projective spaces P^1 x P^1 over Integer Ring must be an affine space or affine subscheme
            sage: Z*X
            Traceback (most recent call last):
            ...
            TypeError: Affine Space of dimension 3 over Integer Ring must be a projective space,
            product of projective spaces, or subscheme
        """
        #This will catch any ambient space mismatches
        AS = self.ambient_space()*right.ambient_space()
        CR = AS.coordinate_ring()
        n = self.ambient_space().coordinate_ring().ngens()

        phi = self.ambient_space().coordinate_ring().hom(list(CR.gens()[:n]), CR)
        psi = right.ambient_space().coordinate_ring().hom(list(CR.gens()[n:]), CR)
        return AS.subscheme([phi(t) for t in self.defining_polynomials()] + [psi(t) for t in right.defining_polynomials()])

    __add__ = union

    def intersection(self, other):
        """
        Return the scheme-theoretic intersection of ``self`` and ``other`` in their
        common ambient space.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: X = A.subscheme([x^2 - y])
            sage: Y = A.subscheme([y])
            sage: X.intersection(Y)
            Closed subscheme of Affine Space of dimension 2 over Integer Ring defined by:
              x^2 - y,
              y
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            raise TypeError("other (=%s) must be a closed algebraic subscheme of an ambient space" % other)
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError("other (=%s) must be in the same ambient space as self" % other)
        return A.subscheme(self.defining_ideal() + other.defining_ideal())

    def complement(self, other=None):
        """
        Return the scheme-theoretic complement ``other - self``, where
        ``self`` and ``other`` are both closed algebraic subschemes of the
        same ambient space.

        If ``other`` is unspecified, it is taken to be the ambient space
        of ``self``.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: X = A.subscheme([x + y - z])
            sage: Y = A.subscheme([x - y + z])
            sage: Y.complement(X)
            Quasi-affine subscheme X - Y of
             Affine Space of dimension 3 over Integer Ring,
             where X is defined by:
               x + y - z
             and Y is defined by:
               x - y + z
            sage: Y.complement()
            Quasi-affine subscheme X - Y of
             Affine Space of dimension 3 over Integer Ring,
             where X is defined by:
               (no polynomials)
             and Y is defined by:
               x - y + z
            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: X = P.subscheme([x^2 + y^2 + z^2])
            sage: Y = P.subscheme([x*y + y*z + z*x])
            sage: Y.complement(X)
            Quasi-projective subscheme X - Y of
             Projective Space of dimension 2 over Rational Field,
             where X is defined by:
               x^2 + y^2 + z^2
             and Y is defined by:
               x*y + x*z + y*z
            sage: Y.complement(P)
            Quasi-projective subscheme X - Y of
             Projective Space of dimension 2 over Rational Field,
             where X is defined by:
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
                raise TypeError("Argument other (=%s) must be a closed algebraic subscheme of an ambient space" % other)
        if other.ambient_space() != A:
            raise ValueError("other (=%s) must be in the same ambient space as self" % other)
        return AlgebraicScheme_quasi(other, self)

    def rational_points(self, **kwds):
        """
        Return the rational points on the algebraic subscheme.

        For a dimension 0 subscheme, if the base ring is a numerical field
        such as the ComplexField the results returned could be very far from correct.
        If the polynomials defining the subscheme are defined over a number field, you
        will get better results calling rational points with `F` defined as the number
        field and the base ring as the field of definition. If the base ring
        is a number field, the embedding into ``F`` must be known.

        In the case of numerically approximated points, the points are returned over as
        points of the ambient space.

        For a dimension greater than 0 scheme, depending on bound size, either the
        points in the ambient space are enumerated or a sieving algorithm lifting points
        modulo primes is used. See the documentation in homset for the details of the
        sieving algorithm.

        INPUT: keyword arguments:

        - ``bound`` -- integer (default: 0); the bound for the coordinates for
          subschemes with dimension at least 1

        - ``prec`` -- integer (default: 53); the precision to use to
          compute the elements of bounded height for number fields

        - ``F`` -- field (default: base ring). The field to compute
          the rational points over

        - ``point_tolerance`` -- positive real number (default: 10^(-10)).
          For numerically inexact fields, two points are considered the same
          if their coordinates are within tolerance.

        - ``zero_tolerance`` -- positive real number (default: 10^(-10)).
          For numerically inexact fields, points are on the subscheme if they
          satisfy the equations to within tolerance.

        - ``tolerance`` -- a rational number in (0,1] used in Doyle-Krumm
          algorithm-4

        OUTPUT: list of points in subscheme or ambient space

        .. WARNING::

           For numerically inexact fields such as ComplexField or RealField the
           list of points returned is very likely to be incomplete at best.

        EXAMPLES:

        Enumerate over a projective scheme over a number field::

            sage: # needs sage.rings.number_field
            sage: u = QQ['u'].0
            sage: K.<v> = NumberField(u^2 + 3)
            sage: A.<x,y> = ProjectiveSpace(K, 1)
            sage: X = A.subscheme(x^2 - y^2)
            sage: X.rational_points(bound=3)                                            # needs sage.libs.singular
            [(-1 : 1), (1 : 1)]

        One can enumerate points up to a given bound on a projective scheme
        over the rationals::

            sage: E = EllipticCurve('37a')                                              # needs sage.schemes
            sage: E.rational_points(bound=8)                                            # needs sage.libs.singular sage.schemes
            [(-1 : -1 : 1), (-1 : 0 : 1), (0 : -1 : 1), (0 : 0 : 1), (0 : 1 : 0),
             (1/4 : -5/8 : 1), (1/4 : -3/8 : 1), (1 : -1 : 1), (1 : 0 : 1),
             (2 : -3 : 1), (2 : 2 : 1)]

        For a small finite field, the complete set of points can be
        enumerated. ::

            sage: Etilde = E.base_extend(GF(3))                                         # needs sage.schemes
            sage: Etilde.rational_points()                                              # needs sage.libs.singular sage.schemes
            [(0 : 1 : 0), (0 : 0 : 1), (0 : 2 : 1), (1 : 0 : 1),
             (1 : 2 : 1), (2 : 0 : 1), (2 : 2 : 1)]

        The class of hyperelliptic curves does not (yet) support
        desingularization of the places at infinity into two points::

            sage: FF = FiniteField(7)
            sage: P.<x> = PolynomialRing(FiniteField(7))
            sage: C = HyperellipticCurve(x^8 + x + 1)                                   # needs sage.schemes
            sage: C.rational_points()                                                   # needs sage.libs.singular sage.schemes
            [(0 : 1 : 0), (0 : 1 : 1), (0 : 6 : 1), (2 : 0 : 1),
             (4 : 0 : 1), (6 : 1 : 1), (6 : 6 : 1)]

        ::

            sage: # needs sage.rings.number_field sage.rings.real_mpfr
            sage: K.<v> = QuadraticField(-3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: X = P.subscheme([x^2 - v^2*x*z, y*x - v*z^2])
            sage: X.rational_points(F=CC)                                               # needs sage.libs.singular
            [(-3.00000000000000 : -0.577350269189626*I : 1.00000000000000),
             (0.000000000000000 : 1.00000000000000 : 0.000000000000000)]

        ::

            sage: # needs sage.rings.number_field sage.rings.real_mpfr
            sage: K.<v> = QuadraticField(3)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: X = A.subscheme([x^2 - v^2*y, y*x - v])
            sage: X.rational_points(F=RR)                                               # needs sage.libs.singular
            [(1.73205080756888, 1.00000000000000)]

        .. TODO::

            Implement Stoll's model in weighted projective space to
            resolve singularities and find two points (1 : 1 : 0) and
            (-1 : 1 : 0) at infinity.
        """
        F = kwds.pop('F', None)
        if F is None: #sometimes None is passed in
            F = self.base_ring()
        if F in NumberFields() or F == ZZ:
            X = self.base_extend(F)(F)
            try:
                return X.points(**kwds) # checks for proper bound done in points functions
            except TypeError:
                raise TypeError("Unable to enumerate points over %s." % F)
        elif (self.base_ring() in NumberFields() or self.base_ring() == ZZ)\
          and hasattr(F, 'precision'):
            #we are numerically approximating number field points
            return self(self.base_ring()).numerical_points(F=F, **kwds)
        try:
            X = self.base_extend(F)(F)
            return X.points()
        except TypeError:
            raise TypeError("Unable to enumerate points over %s." % F)

    def change_ring(self, R):
        r"""
        Return a new algebraic subscheme which is this subscheme coerced to ``R``.

        INPUT:

        - ``R`` -- ring or morphism

        OUTPUT: a new algebraic subscheme which is this subscheme coerced to ``R``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: X = P.subscheme([3*x^2 - y^2])
            sage: H = Hom(X, X)
            sage: X.change_ring(GF(3))
            Closed subscheme of Projective Space of dimension 1
             over Finite Field of size 3 defined by: -y^2

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(2)
            sage: R.<z> = K[]
            sage: L.<v> = K.extension(z^3 - 5)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: X = P.subscheme(x - w*y)                                              # needs sage.libs.singular
            sage: X.change_ring(L)                                                      # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1 over
             Number Field in v with defining polynomial z^3 - 5 over its base field
             defined by: x + (-w)*y

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(2)
            sage: R.<z> = K[]
            sage: L.<v> = K.extension(z^3 - 5)
            sage: P.<x,y,z> = AffineSpace(L, 3)
            sage: X = P.subscheme([x - w*y, z^2 - v*x])                                 # needs sage.libs.singular
            sage: emb = L.embeddings(QQbar)                                             # needs sage.libs.singular
            sage: X.change_ring(emb[0])                                                 # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 3 over Algebraic Field
             defined by:
              x + (-1.414213562373095? + 0.?e-16*I)*y,
              z^2 + (0.8549879733383485? + 1.480882609682365?*I)*x

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(2)
            sage: R.<z> = K[]
            sage: L.<v> = K.extension(z^3 - 5)
            sage: P.<x,y,z> = AffineSpace(L, 3)
            sage: X = P.subscheme([x - w*y, z^2 - v*x])                                 # needs sage.libs.singular
            sage: emb = L.embeddings(QQbar)                                             # needs sage.libs.singular
            sage: X.change_ring(emb[1])                                                 # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 3 over Algebraic Field
             defined by:
              x + (-1.414213562373095? + 0.?e-16*I)*y,
              z^2 + (0.8549879733383485? - 1.480882609682365?*I)*x

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(-3)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: X = P.subscheme(x - w*y)                                              # needs sage.libs.singular
            sage: X.change_ring(CC)                                                     # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1
             over Complex Field with 53 bits of precision defined by:
              x + (-1.73205080756888*I)*y

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(3)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: X = P.subscheme(x - w*y)                                              # needs sage.libs.singular
            sage: X.change_ring(RR)                                                     # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1
             over Real Field with 53 bits of precision defined by:
              x - 1.73205080756888*y

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = CyclotomicField(7)
            sage: O = K.maximal_order()
            sage: P.<x,y> = ProjectiveSpace(O, 1)
            sage: X = P.subscheme([x^2 + O(v)*y^2])                                     # needs sage.libs.singular
            sage: X.change_ring(CC)                                                     # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1
             over Complex Field with 53 bits of precision defined by:
              x^2 + (0.623489801858734 + 0.781831482468030*I)*y^2
            sage: X.change_ring(K).change_ring(K.embeddings(QQbar)[3])                  # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1
             over Algebraic Field defined by:
              x^2 + (-0.9009688679024191? - 0.4338837391175581?*I)*y^2

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: f = x^6 - 2
            sage: L.<b> = NumberField(f, embedding=f.roots(CC)[2][0])
            sage: A.<x,y> = AffineSpace(L, 2)
            sage: H = Hom(A, A)
            sage: X = A.subscheme([b*x^2, y^2])                                         # needs sage.libs.singular
            sage: X.change_ring(CC)                                                     # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 2
             over Complex Field with 53 bits of precision defined by:
              (-0.561231024154687 - 0.972080648619833*I)*x^2,
              y^2
        """
        AS = self.ambient_space()
        new_AS = AS.change_ring(R)
        I = [f.change_ring(R) for f in self.defining_polynomials()]
        return new_AS.subscheme(I)

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

        OUTPUT: affine subscheme; the Weil restriction of ``self``

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^5 - 2)
            sage: R.<x> = K[]
            sage: L.<v> = K.extension(x^2 + 1)
            sage: A.<x,y> = AffineSpace(L, 2)
            sage: X = A.subscheme([y^2 - L(w)*x^3 - v])                                 # needs sage.libs.singular
            sage: X.weil_restriction()                                                  # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 4
             over Number Field in w with defining polynomial x^5 - 2 defined by:
              (-w)*z0^3 + (3*w)*z0*z1^2 + z2^2 - z3^2,
              (-3*w)*z0^2*z1 + w*z1^3 + 2*z2*z3 - 1
            sage: X.weil_restriction().ambient_space() is A.weil_restriction()          # needs sage.libs.singular
            True

        ::

            sage: A.<x,y,z> = AffineSpace(GF(5^2, 't'), 3)                              # needs sage.rings.finite_rings
            sage: X = A.subscheme([y^2 - x*z, z^2 + 2*y])                               # needs sage.libs.singular sage.rings.finite_rings
            sage: X.weil_restriction()                                                  # needs sage.libs.singular sage.rings.finite_rings
            Closed subscheme of Affine Space of dimension 6
             over Finite Field of size 5 defined by:
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

    def specialization(self, D=None, phi=None):
        r"""
        Specialization of this subscheme.

        Given a family of maps defined over a polynomial ring. A specialization
        is a particular member of that family. The specialization can be specified either
        by a dictionary or a :class:`SpecializationMorphism`.

        INPUT:

        - ``D`` -- dictionary (optional)

        - ``phi`` -- :class:`SpecializationMorphism` (optional)

        OUTPUT: :class:`SchemeMorphism_polynomial`

        EXAMPLES::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: X = P.subscheme([x^2 + c*y^2])
            sage: X.specialization(dict({c:2}))
            Closed subscheme of Projective Space of dimension 1 over Rational Field defined by:
              x^2 + 2*y^2

        ::

            sage: R.<c> = PolynomialRing(QQ)
            sage: S.<a,b> = R[]
            sage: P.<x,y,z> = AffineSpace(S, 3)
            sage: X = P.subscheme([x^2 + a*c*y^2 - b*z^2])
            sage: from sage.rings.polynomial.flatten import SpecializationMorphism
            sage: phi = SpecializationMorphism(P.coordinate_ring(),
            ....:                              dict({c: 2, a: 1}))
            sage: X.specialization(phi=phi)                                             # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 3
             over Univariate Polynomial Ring in b over Rational Field defined by:
              x^2 + 2*y^2 + (-b)*z^2
        """
        if D is None:
            if phi is None:
                raise ValueError("either the dictionary or the specialization must be provided")
        else:
            from sage.rings.polynomial.flatten import SpecializationMorphism
            phi = SpecializationMorphism(self.ambient_space().coordinate_ring(),D)
        amb = self.ambient_space().change_ring(phi.codomain().base_ring())
        return amb.subscheme([phi(g) for g in self.defining_polynomials()])
