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

from __future__ import absolute_import, print_function
from sage.categories.homset import End
from six import add_metaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.schemes.affine.affine_space import is_AffineSpace
from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine

@add_metaclass(InheritComparisonClasscallMetaclass)
class DynamicalSystem(SchemeMorphism_polynomial):
    r"""
    Base class for dynamical systems of schemes.

    INPUT:

    - ``polys_or_rat_fncts`` -- a list of polynomials or rational functions,
      all of which should have the same parent

    - ``domain`` -- an affine or projective scheme, or product of
      projective schemes, on which ``polys`` defines an endomorphism.
      Subschemes are also ok

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
        projective space of appropriate dimension over the base ring of
        the first element of ``morphism_or_polys``

      * ``morphism_or_polys`` is a single polynomial or rational
        function; ``domain`` is ignored and taken to be a
        1-dimensional projective space over the base ring of
        ``morphism_or_polys`` with coordinate names given by ``names``

    EXAMPLES::

        sage: A.<x> = AffineSpace(QQ,1)
        sage: f = DynamicalSystem_affine([x^2+1])
        sage: type(f)
        <class 'sage.dynamics.arithmetic_dynamics.affine_ds.DynamicalSystem_affine_field'>

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem_projective([x^2+y^2, y^2])
        sage: type(f)
        <class 'sage.dynamics.arithmetic_dynamics.projective_ds.DynamicalSystem_projective_field'>

    ::

        sage: P1.<x,y> = ProjectiveSpace(CC,1)
        sage: H = End(P1)
        sage: DynamicalSystem(H([y, x]))
        Dynamical System of Projective Space of dimension 1 over Complex Field
        with 53 bits of precision
          Defn: Defined on coordinates by sending (x : y) to
                (y : x)

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
              Defn: Defined on coordinates by sending (x, y) to
                    (y, x)
            sage: H = End(A)
            sage: DynamicalSystem(H([y, x]))
            Dynamical System of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (y, x)
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
        if not domain is None:
            if is_AffineSpace(domain) or isinstance(domain, AlgebraicScheme_subscheme_affine):
                from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine
                return DynamicalSystem_affine(morphism_or_polys, domain)

        from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective
        return DynamicalSystem_projective(morphism_or_polys, domain, names)

    def __init__(self, polys_or_rat_fncts, domain):
        r"""
        The Python constructor.

        EXAMPLES::

            sage: from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2+y^2, y^2])
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
        s = "%s of %s"%(self._repr_type(), self.domain())
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s"%('\n        '.join(self._repr_defn().split('\n')))
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
            sage: f = DynamicalSystem_projective([x^2-y^2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space_field'>

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(5), 1)
            sage: f = DynamicalSystem_projective([x^2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space_finite_field'>

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space'>

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space_field'>

        ::

            sage: A.<x,y> = AffineSpace(GF(3), 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
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
              Defn: Defined on coordinates by sending (x : y) to
                    (-2*x^2 : y^2)
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
            sage: f = DynamicalSystem_projective([x^2 + c*y^2,y^2], domain=P)
            sage: f.specialization({c:1})
            Dynamical System of Projective Space of dimension 1 over Rational Field
                  Defn: Defined on coordinates by sending (x : y) to
                        (x^2 + y^2 : y^2)
        """
        F = self.as_scheme_morphism().specialization(D, phi, homset)
        return F.as_dynamical_system()

