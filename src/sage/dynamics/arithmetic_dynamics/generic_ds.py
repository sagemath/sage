r"""
Generic dynamical systems on schemes

This is the generic class for dynamical systems and contains the exported
constructor functions. The constructor functions can take either polynomials
(or rational functions in the affine case) or morphisms from which to
construct a dynamical system. If the domain is not specified, it is
constructed. However, if you plan on working with points or subvarieties
in the domain, it recommended to specify the domain. For products of
projective spaces the domain must be specified.

The initialization checks are always performed by the constructor functions.
It is possible, but not recommended, to skip these checks by calling the
class initialization directly.

AUTHORS:

- Ben Hutz (July 2017): initial version
"""

#*****************************************************************************
#       Copyright (C) 2017 Ben Hutz <bn4941@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy

from sage.categories.homset import End
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_import import lazy_import
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.rational_field import QQ
from sage.schemes.affine.affine_space import AffineSpace_generic
from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine
from sage.schemes.berkovich.berkovich_space import Berkovich_Cp
from sage.schemes.generic.morphism import SchemeMorphism_polynomial

lazy_import('sage.rings.algebraic_closure_finite_field', 'AlgebraicClosureFiniteField_generic')
lazy_import('sage.rings.qqbar', 'AlgebraicField_common')


class DynamicalSystem(SchemeMorphism_polynomial,
                      metaclass=InheritComparisonClasscallMetaclass):
    r"""
    Base class for dynamical systems of schemes.

    INPUT:

    - ``polys_or_rat_fncts`` -- list of polynomials or rational functions,
      all of which should have the same parent

    - ``domain`` -- an affine or projective scheme, or product of
      projective schemes, on which ``polys`` defines an endomorphism
      (Subschemes are also ok)

    - ``names`` -- (default: ``('X', 'Y')``) tuple of strings to be used
      as coordinate names for a projective space that is constructed

      The following combinations of ``morphism_or_polys`` and
      ``domain`` are meaningful:

      * ``morphism_or_polys`` is a SchemeMorphism; ``domain`` is
        ignored in this case

      * ``morphism_or_polys`` is a list of homogeneous polynomials
        that define a rational endomorphism of ``domain``

      * ``morphism_or_polys`` is a list of homogeneous polynomials and
        ``domain`` is unspecified; ``domain`` is then taken to be the
        projective space of appropriate dimension over the common parent
        of the elements in ``morphism_or_polys``

      * ``morphism_or_polys`` is a single polynomial or rational
        function; ``domain`` is ignored and taken to be a
        1-dimensional projective space over the base ring of
        ``morphism_or_polys`` with coordinate names given by ``names``

    EXAMPLES::

        sage: A.<x> = AffineSpace(QQ,1)
        sage: f = DynamicalSystem_affine([x^2 + 1])
        sage: type(f)
        <class 'sage.dynamics.arithmetic_dynamics.affine_ds.DynamicalSystem_affine_field'>

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
        sage: type(f)
        <class 'sage.dynamics.arithmetic_dynamics.projective_ds.DynamicalSystem_projective_field'>

    ::

        sage: P1.<x,y> = ProjectiveSpace(CC,1)
        sage: H = End(P1)
        sage: DynamicalSystem(H([y, x]))
        Dynamical System of Projective Space of dimension 1 over Complex Field
        with 53 bits of precision
          Defn: Defined on coordinates by sending (x : y) to (y : x)

    :class:`DynamicalSystem` defaults to projective::

        sage: R.<x,y,z> = QQ[]
        sage: DynamicalSystem([x^2, y^2, z^2])
        Dynamical System of Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to
                (x^2 : y^2 : z^2)

    ::

        sage: A.<x,y> = AffineSpace(QQ, 2)
        sage: DynamicalSystem([y, x], domain=A)
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to (y, x)
        sage: H = End(A)
        sage: DynamicalSystem(H([y, x]))
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to (y, x)

    Note that ``domain`` is ignored if an endomorphism is passed in::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: P2.<x,y> = ProjectiveSpace(CC, 1)
        sage: H = End(P2)
        sage: f = H([CC.0*x^2, y^2])
        sage: g = DynamicalSystem(f, domain=P)
        sage: g.domain()
        Projective Space of dimension 1 over Complex Field with 53 bits of precision

    Constructing a common parent::

        sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
        sage: DynamicalSystem([CC.0*x^2, 4/5*y^2])
        Dynamical System of
         Projective Space of dimension 1 over Complex Field with 53 bits of precision
          Defn: Defined on coordinates by sending (x : y) to
                (1.00000000000000*I*x^2 : 0.800000000000000*y^2)

        sage: # needs sage.rings.finite_rings
        sage: P.<x,y> = ProjectiveSpace(GF(5), 1)
        sage: K.<t> = GF(25)
        sage: DynamicalSystem([GF(5)(3)*x^2, K(t)*y^2])
        Dynamical System of Projective Space of dimension 1 over Finite Field in t of size 5^2
          Defn: Defined on coordinates by sending (x : y) to
                (-2*x^2 : t*y^2)
    """

    @staticmethod
    def __classcall_private__(cls, morphism_or_polys, domain=None, names=None):
        r"""
        Return the appropriate dynamical system on a scheme.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: DynamicalSystem(t^2 - 3)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (X : Y) to
                    (X^2 - 3*Y^2 : Y^2)
        """
        if isinstance(morphism_or_polys, SchemeMorphism_polynomial):
            domain = morphism_or_polys.domain()
        if domain is not None:
            if isinstance(domain, AffineSpace_generic) or isinstance(domain, AlgebraicScheme_subscheme_affine):
                from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine
                return DynamicalSystem_affine(morphism_or_polys, domain)
            if isinstance(domain, Berkovich_Cp):
                from sage.dynamics.arithmetic_dynamics.berkovich_ds import DynamicalSystem_Berkovich
                return DynamicalSystem_Berkovich(morphism_or_polys,domain)

        from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective
        return DynamicalSystem_projective(morphism_or_polys, domain, names)

    def __init__(self, polys_or_rat_fncts, domain):
        r"""
        The Python constructor.

        EXAMPLES::

            sage: from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: isinstance(f, DynamicalSystem)
            True
        """
        H = End(domain)
        # All consistency checks are done by the public class constructors,
        # so we can set check=False here.
        SchemeMorphism_polynomial.__init__(self, H, polys_or_rat_fncts, check=False)

    def _repr_type(self):
        r"""
        Return a string representation of the type of a dynamical system.

        OUTPUT: string

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2])
            sage: f._repr_type()
            'Dynamical System'
        """
        return "Dynamical System"

    def _repr_(self):
        r"""
        Return a string representation of a dynamical system.

        OUTPUT: string

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2])
            sage: f._repr_()
            'Dynamical System of Projective Space of dimension 1 over Rational Field\n
              Defn: Defined on coordinates by sending (x : y) to\n        (x^3 : x*y^2)'
        """
        s = "%s of %s" % (self._repr_type(), self.domain())
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s" % ('\n        '.join(self._repr_defn().split('\n')))
        return s

    def as_scheme_morphism(self):
        """
        Return this dynamical system as :class:`SchemeMorphism_polynomial`.

        OUTPUT: :class:`SchemeMorphism_polynomial`

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space'>

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space_field'>

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(5), 1)
            sage: f = DynamicalSystem_projective([x^2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space_finite_field'>

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: f = DynamicalSystem_affine([x^2 - 2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space'>

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([x^2 - 2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space_field'>

        ::

            sage: A.<x,y> = AffineSpace(GF(3), 2)
            sage: f = DynamicalSystem_affine([x^2 - 2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space_finite_field'>
        """
        H = End(self.domain())
        return H(list(self))

    def change_ring(self, R, check=True):
        r"""
        Return a new dynamical system which is this map coerced to ``R``.

        If ``check`` is ``True``, then the initialization checks are performed.

        INPUT:

        - ``R`` -- ring or morphism

        OUTPUT:

        A new :class:`DynamicalSystem_projective` that is this map
        coerced to ``R``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: f = DynamicalSystem_projective([3*x^2, y^2])
            sage: f.change_ring(GF(5))
            Dynamical System of Projective Space of dimension 1 over Finite Field of size 5
              Defn: Defined on coordinates by sending (x : y) to (-2*x^2 : y^2)
        """
        f = self.as_scheme_morphism()
        F = f.change_ring(R)
        return F.as_dynamical_system()

    def specialization(self, D=None, phi=None, homset=None):
        r"""
        Specialization of this dynamical system.

        Given a family of maps defined over a polynomial ring. A
        specialization is a particular member of that family. The
        specialization can be specified either by a dictionary or
        a :class:`SpecializationMorphism`.

        INPUT:

        - ``D`` -- (optional) dictionary

        - ``phi`` -- (optional) SpecializationMorphism

        - ``homset`` -- (optional) homset of specialized map

        OUTPUT: :class:`DynamicalSystem`

        EXAMPLES::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2], domain=P)
            sage: f.specialization({c:1})
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 + y^2 : y^2)
        """
        F = self.as_scheme_morphism().specialization(D, phi, homset)
        return F.as_dynamical_system()

    def field_of_definition_critical(self, return_embedding=False, simplify_all=False, names='a'):
        r"""
        Return smallest extension of the base field which contains the critical points.

        Ambient space of dynamical system must be either the affine line or projective
        line over a number field or finite field.

        INPUT:

        - ``return_embedding`` -- boolean (default: ``False``); if ``True``, return an
          embedding of base field of dynamical system into the returned number field or
          finite field. Note that computing this embedding might be expensive.

        - ``simplify_all`` -- boolean (default: ``False``); if ``True``, simplify
          intermediate fields and also the resulting number field. Note that this
          is not implemented for finite fields and has no effect.

        - ``names`` -- (optional) string to be used as generator for returned number field
          or finite field

        OUTPUT:

        If ``return_embedding`` is ``False``, the field of definition as an absolute number
        field or finite field.  If ``return_embedding`` is ``True``, a tuple
        ``(K, phi)`` where ``phi`` is an embedding of the base field in ``K``.

        EXAMPLES:

        Note that the number of critical points is `2d-2`, but `(1:0)` has multiplicity 2 in this case::

            sage: # needs sage.libs.singular sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([1/3*x^3 + x*y^2, y^3], domain=P)
            sage: f.critical_points()
            [(1 : 0)]
            sage: N.<a> = f.field_of_definition_critical(); N
            Number Field in a with defining polynomial x^2 + 1
            sage: g = f.change_ring(N)
            sage: g.critical_points()
            [(-a : 1), (a : 1), (1 : 0)]

        ::

            sage: # needs sage.libs.singular sage.rings.number_field
            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem([z^4 + 2*z^2 + 2], domain=A)
            sage: K.<a> = f.field_of_definition_critical(); K
            Number Field in a with defining polynomial z^2 + 1

        ::

            sage: # needs sage.libs.singular sage.rings.finite_rings
            sage: G.<a> = GF(9)
            sage: R.<z> = G[]
            sage: R.irreducible_element(3, algorithm='first_lexicographic')
            z^3 + (a + 1)*z + a
            sage: A.<x> = AffineSpace(G,1)
            sage: f = DynamicalSystem([x^4 + (2*a+2)*x^2 + a*x], domain=A)
            sage: f[0].derivative(x).univariate_polynomial().is_irreducible()
            True
            sage: f.field_of_definition_critical(return_embedding=True, names='b')
            (Finite Field in b of size 3^6,
             Ring morphism:
               From: Finite Field in a of size 3^2
               To:   Finite Field in b of size 3^6
               Defn: a |--> 2*b^5 + 2*b^3 + b^2 + 2*b + 2)
        """
        ds = copy(self)
        space = ds.domain().ambient_space()
        K = ds.base_ring()
        if space.dimension() != 1:
            raise ValueError('ambient space of dynamical system must be either the affine line or projective line')
        if isinstance(K, (AlgebraicClosureFiniteField_generic, AlgebraicField_common)):
            if return_embedding:
                return (K, K.hom(K))
            else:
                return K
        if space.is_projective():
            ds = ds.dehomogenize(1)
        f,g = ds[0].numerator(), ds[0].denominator()
        CR = space.coordinate_ring()
        if CR.is_field():
            #want the polynomial ring not the fraction field
            CR = CR.ring()
        x = CR.gen(0)
        poly = (g*CR(f).derivative(x) - f*CR(g).derivative(x)).univariate_polynomial()
        if isinstance(ds.base_ring(), FiniteField):
            return poly.splitting_field(names, map=return_embedding)
        else:
            K = poly.splitting_field(names, map=return_embedding, simplify_all=simplify_all)
            if return_embedding:
                N = K[0]
            else:
                N = K
            if N.absolute_degree() == 1:
                if return_embedding:
                    return (QQ,ds.base_ring().embeddings(QQ)[0])
                else:
                    return QQ
            else:
                return K

    def field_of_definition_periodic(self, n, formal=False, return_embedding=False, simplify_all=False, names='a'):
        r"""
        Return smallest extension of the base field which contains all fixed points
        of the ``n``-th iterate

        Ambient space of dynamical system must be either the affine line
        or projective line over a number field or finite field.

        INPUT:

        - ``n`` -- positive integer

        - ``formal`` -- boolean (default: ``False``); ``True`` signals to return number
          field or finite field over which the formal periodic points are defined, where a
          formal periodic point is a root of the ``n``-th dynatomic polynomial.
          ``False`` specifies to find number field or finite field over which all periodic
          points of the ``n``-th iterate are defined.

        - ``return_embedding`` -- boolean (default: ``False``); if ``True``, return
          an embedding of base field of dynamical system into the returned number
          field or finite field. Note that computing this embedding might be expensive.

        - ``simplify_all`` -- boolean (default: ``False``); if ``True``, simplify
          intermediate fields and also the resulting number field. Note that this
          is not implemented for finite fields and has no effect.

        - ``names`` -- (optional) string to be used as generator for returned number
          field or finite field

        OUTPUT:

        If ``return_embedding`` is ``False``, the field of definition as an absolute
        number field or finite field.  If ``return_embedding`` is ``True``, a tuple
        ``(K, phi)`` where ``phi`` is an embedding of the base field in ``K``.

        EXAMPLES::

            sage: # needs sage.libs.singular sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x^2, y^2], domain=P)
            sage: f.periodic_points(3, minimal=False)
            [(0 : 1), (1 : 0), (1 : 1)]
            sage: N.<a> = f.field_of_definition_periodic(3); N
            Number Field in a with defining polynomial x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: sorted(f.periodic_points(3, minimal=False, R=N), key=str)
            [(-a^5 - a^4 - a^3 - a^2 - a - 1 : 1),
             (0 : 1),
             (1 : 0),
             (1 : 1),
             (a : 1),
             (a^2 : 1),
             (a^3 : 1),
             (a^4 : 1),
             (a^5 : 1)]

        ::

            sage: # needs sage.libs.singular sage.rings.number_field
            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem([(z^2 + 1)/(2*z + 1)], domain=A)
            sage: K.<a> = f.field_of_definition_periodic(2); K
            Number Field in a with defining polynomial z^4 + 12*z^3 + 39*z^2 + 18*z + 171
            sage: F.<b> = f.field_of_definition_periodic(2, formal=True); F
            Number Field in b with defining polynomial z^2 + 3*z + 6

        ::

            sage: # needs sage.rings.finite_rings
            sage: G.<a> = GF(4)
            sage: A.<x> = AffineSpace(G, 1)
            sage: f = DynamicalSystem([x^2 + (a+1)*x + 1], domain=A)
            sage: g = f.nth_iterate_map(2)[0]
            sage: (g-x).univariate_polynomial().factor()
            (x + 1) * (x + a + 1) * (x^2 + a*x + 1)
            sage: f.field_of_definition_periodic(2, return_embedding=True, names='b')
            (Finite Field in b of size 2^4,
             Ring morphism:
               From: Finite Field in a of size 2^2
               To:   Finite Field in b of size 2^4
               Defn: a |--> b^2 + b)
        """
        ds = copy(self)
        n = int(n)
        K = ds.base_ring()
        if n < 1:
            raise ValueError('`n` must be >= 1')
        space = ds.domain().ambient_space()
        if space.dimension() != 1:
            raise NotImplementedError("not implemented for affine or projective spaces of dimension >1")
        if isinstance(K, (AlgebraicClosureFiniteField_generic, AlgebraicField_common)):
            if return_embedding:
                return (K, K.hom(K))
            else:
                return K
        if space.is_projective():
            ds = ds.dehomogenize(1)
        CR = space.coordinate_ring()
        if CR.is_field():
            #want the polynomial ring not the fraction field
            CR = CR.ring()
        x = CR.gen(0)
        if formal:
            poly = ds.dynatomic_polynomial(n)
            poly = CR(poly).univariate_polynomial()
        else:
            fn = ds.nth_iterate_map(n)
            f,g = fn[0].numerator(), fn[0].denominator()
            poly = (f - g*x).univariate_polynomial()
        if isinstance(ds.base_ring(), FiniteField):
            return poly.splitting_field(names, map=return_embedding)
        else:
            K = poly.splitting_field(names, map=return_embedding, simplify_all=simplify_all)
            if return_embedding:
                N = K[0]
            else:
                N = K
            if N.absolute_degree() == 1:
                if return_embedding:
                    return (QQ,ds.base_ring().embeddings(QQ)[0])
                else:
                    return QQ
            else:
                return K

    def field_of_definition_preimage(self, point, n, return_embedding=False, simplify_all=False, names='a'):
        r"""
        Return smallest extension of the base field which contains the
        ``n``-th preimages of ``point``

        Ambient space of dynamical system must be either the affine line or
        projective line over a number field or finite field.

        INPUT:

        - ``point`` -- a point in this map's domain

        - ``n`` -- positive integer

        - ``return_embedding`` -- boolean (default: ``False``); if ``True``, return
          an embedding of base field of dynamical system into the returned number
          field or finite field. Note that computing this embedding might be expensive.

        - ``simplify_all`` -- boolean (default: ``False``); if ``True``, simplify
          intermediate fields and also the resulting number field. Note that this
          is not implemented for finite fields and has no effect.

        - ``names`` -- (optional) string to be used as generator for returned
          number field or finite field

        OUTPUT:

        If ``return_embedding`` is ``False``, the field of definition as an absolute
        number field or finite field.  If ``return_embedding`` is ``True``, a tuple
        ``(K, phi)`` where ``phi`` is an embedding of the base field in ``K``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([1/3*x^2 + 2/3*x*y, x^2 - 2*y^2], domain=P)
            sage: N.<a> = f.field_of_definition_preimage(P(1,1), 2,                     # needs sage.rings.number_field
            ....:                                        simplify_all=True); N
            Number Field in a with defining polynomial
             x^8 - 4*x^7 - 128*x^6 + 398*x^5 + 3913*x^4 - 8494*x^3 - 26250*x^2 + 30564*x - 2916

        ::

            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem([z^2], domain=A)
            sage: K.<a> = f.field_of_definition_preimage(A(1), 3); K                    # needs sage.rings.number_field
            Number Field in a with defining polynomial z^4 + 1

        ::

            sage: G = GF(5)
            sage: P.<x,y> = ProjectiveSpace(G, 1)
            sage: f = DynamicalSystem([x^2 + 2*y^2, y^2], domain=P)
            sage: f.field_of_definition_preimage(P(2,1), 2, return_embedding=True,      # needs sage.rings.number_field
            ....:                                names='a')
            (Finite Field in a of size 5^2,
             Ring morphism:
               From: Finite Field of size 5
               To:   Finite Field in a of size 5^2
               Defn: 1 |--> 1)
        """
        ds = copy(self)
        n = int(n)
        if n < 1:
            raise ValueError('`n` must be >= 1')
        space = ds.domain().ambient_space()
        if space.dimension() != 1:
            raise NotImplementedError("not implemented for affine or projective spaces of dimension >1")
        try:
            point = space(point)
        except TypeError:
            raise TypeError('`point` must be in {}'.format(ds.domain()))
        if space.is_projective():
            ds = ds.dehomogenize(1)
        else:
            point = (point[0],1)
        fn = ds.nth_iterate_map(n)
        f, g = fn[0].numerator(), fn[0].denominator()
        CR = space.coordinate_ring()
        if CR.is_field():
            #want the polynomial ring not the fraction field
            CR = CR.ring()
        poly = (f*point[1] - g*CR(point[0])).univariate_polynomial()
        if isinstance(ds.base_ring(), FiniteField):
            return poly.splitting_field(names, map=return_embedding)
        else:
            K = poly.splitting_field(names, map=return_embedding, simplify_all=simplify_all)
            if return_embedding:
                N = K[0]
            else:
                N = K
            if N.absolute_degree() == 1:
                if return_embedding:
                    return (QQ, ds.base_ring().embeddings(QQ)[0])
                else:
                    return QQ
            else:
                return K
