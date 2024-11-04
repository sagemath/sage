# sage.doctest: needs sage.rings.padics
r"""
Dynamical systems on Berkovich space over `\CC_p`.

A dynamical system on Berkovich space over `\CC_p` is
determined by a dynamical system on `A^1(\CC_p)` or `P^1(\CC_p)`,
which naturally induces a dynamical system on affine or
projective Berkovich space.

For an exposition of dynamical systems on Berkovich space, see chapter
7 of [Ben2019]_, or for a more involved exposition, chapter 2 of [BR2010]_.

AUTHORS:

 - Alexander Galarraga (August 14th, 2020): initial implementation
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.number_fields import NumberFields
from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective
from sage.matrix.constructor import Matrix
from sage.misc.classcall_metaclass import typecall
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_import import lazy_import
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.schemes.affine.affine_space import AffineSpace_generic
from sage.schemes.berkovich.berkovich_space import (Berkovich_Cp_Affine,
                                Berkovich_Cp_Projective, Berkovich_Cp,
                                Berkovich_Element_Cp_Affine)
from sage.schemes.projective.projective_space import ProjectiveSpace_ring
from sage.structure.element import Element

lazy_import('sage.rings.padics.padic_base_generic', 'pAdicBaseGeneric')


class DynamicalSystem_Berkovich(Element, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A dynamical system on Berkovich space over `\CC_p`.

    A dynamical system on Berkovich space over `\CC_p` is
    determined by a dynamical system on `A^1(\CC_p)` or `P^1(\CC_p)`,
    which naturally induces a dynamical system on affine or
    projective Berkovich space.

    INPUT:

    - ``dynamical_system`` -- a :class:`DynamicalSystem`
      over affine or projective space. If this input is not defined
      over a `p`-adic field, then ``domain`` MUST be specified.

    - ``domain`` -- (optional) affine or projective Berkovich space
      over `\CC_p`. ``domain`` must be specified if ``dynamical_system``
      is defined over a number field.

    - ``ideal`` -- (optional) an ideal of the ``base_ring`` of the domain
      of ``dynamical_system``. Used to create ``domain`` as a Berkovich
      space backed by a number field more efficiently, see examples.

    EXAMPLES:

    We can easily create a dynamical system on Berkovich space
    using a dynamical system on projective space over `\QQ_p`::

        sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
        sage: f = DynamicalSystem_projective([2*x^2 + 4*y^2, 3*x^2 + 9*y^2])
        sage: DynamicalSystem_Berkovich(f)
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20
         induced by the map
          Defn: Defined on coordinates by sending (x : y) to
                ((2 + O(3^20))*x^2 + (1 + 3 + O(3^20))*y^2
                 : (3 + O(3^21))*x^2 + (3^2 + O(3^22))*y^2)

    Or directly from polynomials::

        sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: DynamicalSystem_Berkovich([x^2 + y^2, y^2])
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20
         induced by the map
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 + y^2 : y^2)

    :class:`DynamicalSystem_Berkovich` defaults to projective::

        sage: R.<x,y> = Qp(3)[]
        sage: DynamicalSystem_Berkovich([x^2, y^2])
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20
         induced by the map
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

    To create an affine dynamical system on Berkovich space, pass an
    affine dynamical system to :class:`DynamicalSystem_Berkovich`::

        sage: A.<z> = AffineSpace(Qp(3), 1)
        sage: f = DynamicalSystem_affine(z^2 + 1)
        sage: DynamicalSystem_Berkovich(f)
        Dynamical system of Affine Berkovich line over Cp(3) of precision 20 induced by the map
          Defn: Defined on coordinates by sending (z) to
                (z^2 + 1 + O(3^20))

    ``domain`` can be used to specify the type of dynamical system::

        sage: A.<z> = AffineSpace(Qp(3), 1)
        sage: C = Berkovich_Cp_Affine(3)
        sage: DynamicalSystem_Berkovich([z^2 + 1], C)
        Dynamical system of Affine Berkovich line over Cp(3) of precision 20 induced by the map
          Defn: Defined on coordinates by sending (z) to
                (z^2 + 1 + O(3^20))

    We can create dynamical systems which act on Berkovich spaces backed by number fields::

        sage: # needs sage.rings.number_field
        sage: R.<z> = QQ[]
        sage: A.<a> = NumberField(z^2 + 1)
        sage: ideal = A.prime_above(2)
        sage: P.<x,y> = ProjectiveSpace(A, 1)
        sage: B = Berkovich_Cp_Projective(P, ideal)
        sage: DynamicalSystem_Berkovich([x^2 + y^2, 2*a*x*y], B)
        Dynamical system of Projective Berkovich line over Cp(2), with base Number Field
        in a with defining polynomial z^2 + 1 induced by the map
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 + y^2 : (2*a)*x*y)

    We can use the optional parameter ``ideal`` to create the
    same dynamical system more efficiently::

        sage: # needs sage.rings.number_field
        sage: R.<z> = QQ[]
        sage: A.<a> = NumberField(z^2 + 1)
        sage: prime_ideal = A.prime_above(2)
        sage: P.<x,y> = ProjectiveSpace(A, 1)
        sage: DynamicalSystem_Berkovich([x^2 + y^2, 2*a*x*y], ideal=prime_ideal)
        Dynamical system of Projective Berkovich line over Cp(2), with base Number Field
        in a with defining polynomial z^2 + 1 induced by the map
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 + y^2 : (2*a)*x*y)

    Creating a map on Berkovich space
    creates the Berkovich space it acts on::

        sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
        sage: f = DynamicalSystem_projective([x^2, y^2])
        sage: g = DynamicalSystem_Berkovich(f)
        sage: B = g.domain(); B
        Projective Berkovich line over Cp(3) of precision 20

    The image of type I point is the image of the center::

        sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
        sage: F = DynamicalSystem_Berkovich([x^2, y^2])
        sage: B = F.domain()
        sage: Q1 = B(2)
        sage: F(Q1)
        Type I point centered at (1 + 3 + O(3^20) : 1 + O(3^20))

    For type II/III points with no poles in the corresponding disk,
    the image is the type II/III point corresponding to the image
    of the disk::

        sage: Q2 = B(0, 3)
        sage: F(Q2)
        Type II point centered at (0 : 1 + O(3^20)) of radius 3^2

    The image of any type II point can be computed::

        sage: g = DynamicalSystem_projective([x^2 + y^2, x*y])
        sage: G = DynamicalSystem_Berkovich(g)
        sage: Q3 = B(0, 1)
        sage: G(Q3)
        Type II point centered at (0 : 1 + O(3^20)) of radius 3^0

    The image of type III points can be computed has long as the
    corresponding disk contains no poles of the dynamical system::

        sage: Q4 = B(1/9, 1.5)
        sage: G(Q4)
        Type III point centered at (3^-2 + 3^2 + O(3^18) : 1 + O(3^20))
        of radius 1.50000000000000

    Sometimes, however, the poles are contained in an extension of
    `\QQ_p` that Sage does not support::

        sage: H = DynamicalSystem_Berkovich([x*y^2, x^3 + 20*y^3])
        sage: H(Q4) # not tested
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot check if poles lie in type III disk

    ``Q4``, however, does not contain any poles of ``H`` (this
    can be checked using pencil and paper or the number field functionality
    in Sage). There are two ways around this error: the first and simplest is
    to have ``H`` act on a Berkovich space backed by a number field::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: B = Berkovich_Cp_Projective(P, 3)
        sage: H = DynamicalSystem_Berkovich([x*y^2, x^3 + 20*y^3], B)
        sage: Q4 = B(1/9, 1.5)
        sage: H(Q4)                                                                     # needs sage.rings.number_field
        Type III point centered at (81/14581 : 1) of radius 0.00205761316872428

    Alternatively, if checking for poles in the disk has been done already,
    ``type_3_pole_check`` can be set to ``False``::

        sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
        sage: H = DynamicalSystem_Berkovich([x*y^2, x^3 + 20*y^3])
        sage: B = H.domain()
        sage: Q4 = B(1/9, 1.5)
        sage: H(Q4, False)
        Type III point centered at (3^4 + 3^10 + 2*3^11 + 2*3^13 + 2*3^14 +
        2*3^15 + 3^17 + 2*3^18 + 2*3^19 + 3^20 + 3^21 + 3^22 + O(3^24) : 1 +
        O(3^20)) of radius 0.00205761316872428

    WARNING: setting ``type_3_pole_check`` to ``False`` can lead to
    mathematically incorrect answers.
    """

    @staticmethod
    def __classcall_private__(cls, dynamical_system, domain=None, ideal=None):
        """
        Return the appropriate dynamical system on Berkovich space.

        EXAMPLES::

            sage: R.<t> = Qp(3)[]
            sage: f = DynamicalSystem_affine(t^2 - 3)
            sage: DynamicalSystem_Berkovich(f)
            Dynamical system of Affine Berkovich line over Cp(3) of precision 20 induced by the map
              Defn: Defined on coordinates by sending ((1 + O(3^20))*t) to
                    ((1 + O(3^20))*t^2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 +
                    2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + 2*3^12 + 2*3^13 + 2*3^14 +
                    2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19 + 2*3^20 + O(3^21))
        """
        if not (isinstance(domain, Berkovich_Cp) or domain is None):
            raise TypeError('domain must be a Berkovich space over Cp, not %s' % domain)

        if isinstance(domain, Berkovich_Cp_Affine):
            if not isinstance(dynamical_system, DynamicalSystem_affine):
                try:
                    dynamical_system = DynamicalSystem_affine(dynamical_system)
                except (TypeError, ValueError):
                    raise TypeError('domain was affine Berkovich space, but dynamical_system did not '
                        'convert to an affine dynamical system')
        if isinstance(domain, Berkovich_Cp_Projective):
            if not isinstance(dynamical_system, DynamicalSystem_projective):
                try:
                    dynamical_system = DynamicalSystem_projective(dynamical_system)
                except (TypeError, ValueError):
                    raise TypeError('domain was projective Berkovich space, but dynamical_system did not convert '
                        'to a projective dynamical system')

        if not isinstance(dynamical_system, DynamicalSystem):
            try:
                dynamical_system = DynamicalSystem(dynamical_system)
            except (TypeError, ValueError):
                raise TypeError('dynamical_system did not convert to a dynamical system')
        morphism_domain = dynamical_system.domain()

        if not isinstance(morphism_domain.base_ring(), pAdicBaseGeneric):
            if morphism_domain.base_ring() in NumberFields():
                if domain is None and ideal is not None:
                    if isinstance(morphism_domain, AffineSpace_generic):
                        domain = Berkovich_Cp_Affine(morphism_domain.base_ring(), ideal)
                    else:
                        domain = Berkovich_Cp_Projective(morphism_domain, ideal)
                else:
                    if ideal is not None:
                        if ideal != domain.ideal():
                            raise ValueError('conflicting inputs for ideal and domain')
            else:
                raise ValueError('base ring of domain of dynamical_system must be p-adic or a number field '
                    'not %s' % morphism_domain.base_ring())

        if isinstance(morphism_domain, AffineSpace_generic):
            return DynamicalSystem_Berkovich_affine(dynamical_system, domain)

        return DynamicalSystem_Berkovich_projective(dynamical_system, domain)

    def __init__(self, dynamical_system, domain):
        r"""
        The Python constructor.

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_projective([2*x^2 + 4*y^2, 3*x^2 + 9*y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: isinstance(g, DynamicalSystem_Berkovich)
            True
        """
        self._system = dynamical_system
        self._domain = domain

    def __eq__(self, other):
        """
        Equality operator.

        EXAMPLES::

            sage: R.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_Berkovich([x^2, y^2])
            sage: f == f
            True

        ::

            sage: g = DynamicalSystem_Berkovich([x^3, x*y^2])
            sage: f == g
            True
        """
        if not isinstance(other, type(self)):
            return False
        return self._system == other._system

    def __neq__(self, other):
        """
        Inequality operator.

        EXAMPLES::

            sage: R.<x, y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_Berkovich([x^2, y^2])
            sage: f != f
            False

        ::

            sage: g = DynamicalSystem_Berkovich([x^2 + y^2, y^2])
            sage: f != g
            True
        """
        return not (self == other)

    def domain(self):
        """
        Return the domain of this dynamical system.

        OUTPUT: a Berkovich space over ``Cp``

        EXAMPLES::

            sage: Q.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_projective([3*x^2, 2*y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: g.domain()
            Projective Berkovich line over Cp(3) of precision 20
        """
        return self._domain

    def as_scheme_dynamical_system(self):
        r"""
        Return this dynamical system as :class:`DynamicalSystem`.

        OUTPUT: an affine or projective :class:`DynamicalSystem`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_Berkovich([x^2 + y^2, x*y])
            sage: f.as_scheme_dynamical_system()
            Dynamical System of Projective Space of dimension 1 over
             3-adic Field with capped relative precision 20
              Defn: Defined on coordinates by sending (x : y) to (x^2 + y^2 : x*y)
        """
        return self._system

    def __getitem__(self, i):
        """
        Return the `i`-th polynomial.

        INPUT:

        - ``i`` -- integer

        OUTPUT: an element of polynomial ring or a
        fraction field of a polynomial ring

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, 2*y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: g[0]
            x^2 + y^2
        """
        return self._system._polys[i]

    def defining_polynomials(self):
        """
        Return the defining polynomials.

        OUTPUT: a tuple of polynomials that defines the
        dynamical system.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_projective([2*x^2 + 4*y^2, 3*x^2 + 9*y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: g.defining_polynomials()
            ((2 + O(3^20))*x^2 + (1 + 3 + O(3^20))*y^2,
             (3 + O(3^21))*x^2 + (3^2 + O(3^22))*y^2)
        """
        return self._system._polys

    def base_ring(self):
        """
        The base ring of this dynamical system, that is, the field of definition of the coefficients.

        OUTPUT: a field

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_Berkovich([x^2 + y^2, y^2])
            sage: f.base_ring()
            3-adic Field with capped relative precision 20

        ::

            sage: # needs sage.rings.number_field
            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^3 + 20)
            sage: P.<x,y> = ProjectiveSpace(A, 1)
            sage: f = DynamicalSystem_Berkovich([x^2, x^2 + y^2], ideal=A.prime_above(2))
            sage: f.base_ring()
            Number Field in a with defining polynomial z^3 + 20
        """
        return self.domain().base_ring()

    def _repr_(self):
        r"""
        Return a string representation of this dynamical system.

        OUTPUT: string

        EXAMPLES::

            sage: Q.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_projective([3*x^2, 2*y^2])
            sage: f = DynamicalSystem_Berkovich(f)
            sage: f._repr_()
            'Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map\n
              Defn: Defined on coordinates by sending (x : y) to\n        ((3 + O(3^21))*x^2 : (2 + O(3^20))*y^2)'
        """
        domain_str = self._domain._repr_()
        return "Dynamical system of " + domain_str + " induced by the map" + \
            "\n  Defn: %s" % ('\n        '.join(self._system._repr_defn().split('\n')))


class DynamicalSystem_Berkovich_projective(DynamicalSystem_Berkovich):
    r"""
    A dynamical system on projective Berkovich space over `\CC_p`.

    A dynamical system on projective Berkovich space over `\CC_p` is
    determined by a dynamical system on `A^1(\CC_p)` or `P^1(\CC_p)`,
    which naturally induces a dynamical system on affine or
    projective Berkovich space.

    INPUT:

    - ``dynamical_system`` -- a :class:`DynamicalSystem_Projective`
      of relative dimension 1. If this input is not defined
      over a `p`-adic field, then ``domain`` MUST be specified.

    - ``domain`` -- (optional) projective Berkovich space
      over `\CC_p`. If the input to ``dynamical_system`` is
      not defined over a `p`-adic field, ``domain``
      must be specified.

    EXAMPLES:

    We can easily create a dynamical system on Berkovich space
    using a dynamical system on projective space over `\QQ_p`::

        sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
        sage: f = DynamicalSystem_projective([1/2*x^2 + x*y + 3*y^2, 3*x^2 + 9*y^2])
        sage: DynamicalSystem_Berkovich(f)
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20
        induced by the map
          Defn: Defined on coordinates by sending (x : y) to
                ((2 + 3 + 3^2 + 3^3 + 3^4 + 3^5 + 3^6 + 3^7 + 3^8 + 3^9 + 3^10 + 3^11
                + 3^12 + 3^13 + 3^14 + 3^15 + 3^16 + 3^17 + 3^18 + 3^19 + O(3^20))*x^2
                + x*y + (3 + O(3^21))*y^2 : (3 + O(3^21))*x^2 + (3^2 + O(3^22))*y^2)

    Or from a morphism::

        sage: P1.<x,y> = ProjectiveSpace(Qp(3), 1)
        sage: H = End(P1)
        sage: DynamicalSystem_Berkovich(H([y, x]))
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20
         induced by the map
          Defn: Defined on coordinates by sending (x : y) to (y : x)

    Or from polynomials::

        sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
        sage: DynamicalSystem_Berkovich([x^2+y^2, y^2])
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20
         induced by the map
          Defn: Defined on coordinates by sending (x : y) to (x^2 + y^2 : y^2)
    """
    @staticmethod
    def __classcall_private__(cls, dynamical_system, domain=None):
        """
        Return the approapriate dynamical system on projective Berkovich space over ``Cp``.

        EXAMPLES::

            sage: P1.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: from sage.dynamics.arithmetic_dynamics.berkovich_ds import DynamicalSystem_Berkovich_projective
            sage: DynamicalSystem_Berkovich_projective([y, x])
            Dynamical system of Projective Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x : y) to (y : x)
        """
        if not isinstance(dynamical_system, DynamicalSystem):
            if not isinstance(dynamical_system, DynamicalSystem_projective):
                dynamical_system = DynamicalSystem_projective(dynamical_system)
            else:
                raise TypeError('affine dynamical system passed to projective constructor')
        R = dynamical_system.base_ring()
        morphism_domain = dynamical_system.domain()
        if not isinstance(morphism_domain, ProjectiveSpace_ring):
            raise TypeError('the domain of dynamical_system must be projective space, not %s' % morphism_domain)
        if morphism_domain.dimension_relative() != 1:
            raise ValueError('domain was not relative dimension 1')
        if not isinstance(R, pAdicBaseGeneric):
            if domain is None:
                raise TypeError('dynamical system defined over %s, not p-adic, ' % morphism_domain.base_ring() +
                    'and domain is None')
            if not isinstance(domain, Berkovich_Cp_Projective):
                raise TypeError('domain was %s, not a projective Berkovich space over Cp' % domain)
            if domain.base() != morphism_domain:
                raise ValueError('base of domain was %s, with coordinate ring %s ' % (domain.base(),
                    domain.base().coordinate_ring()) + 'while dynamical_system acts on %s, ' % morphism_domain +
                        'with coordinate ring %s' % morphism_domain.coordinate_ring())
        else:
            domain = Berkovich_Cp_Projective(morphism_domain)
        return typecall(cls, dynamical_system, domain)

    def __init__(self, dynamical_system, domain=None):
        """
        Python constructor.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: DynamicalSystem_Berkovich([x^2 + x*y + 2*y^2, 2*x*y])
            Dynamical system of Projective Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + x*y + (2 + O(3^20))*y^2 : (2 + O(3^20))*x*y)
        """
        DynamicalSystem_Berkovich.__init__(self, dynamical_system, domain)

    def scale_by(self, t):
        """
        Scale each coordinate of this dynamical system by a factor of `t`.

        INPUT:

        - ``t`` -- a ring element

        OUTPUT: none

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_Berkovich([x^2, y^2])
            sage: f.scale_by(x); f
            Dynamical system of Projective Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x : y) to (x^3 : x*y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: Q.<z> = QQ[]
            sage: A.<a> = NumberField(z^3 + 20)
            sage: ideal = A.prime_above(3)
            sage: P.<x,y> = ProjectiveSpace(A, 1)
            sage: B = Berkovich_Cp_Projective(P, ideal)
            sage: f = DynamicalSystem_Berkovich([x^2 + y^2, 2*x*y], B)
            sage: f.scale_by(2); f
            Dynamical system of Projective Berkovich line over Cp(3), with base Number
            Field in a with defining polynomial z^3 + 20 induced by the map
              Defn: Defined on coordinates by sending (x : y) to
                    (2*x^2 + 2*y^2 : 4*x*y)
        """
        self._system.scale_by(t)

    def normalize_coordinates(self):
        r"""
        Normalize the coordinates of the inducing map.

        OUTPUT: none

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_Berkovich([2*x^2, 2*y^2])
            sage: f.normalize_coordinates(); f
            Dynamical system of Projective Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)


        Normalize_coordinates may sometimes fail over `p`-adic fields::

            sage: g = DynamicalSystem_Berkovich([2*x^2, x*y])
            sage: g.normalize_coordinates()  # not tested
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce since the denominator is not 1

        To fix this issue, create a system on Berkovich space backed
        by a number field::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: B = Berkovich_Cp_Projective(P, 3)
            sage: g = DynamicalSystem_Berkovich([2*x^2, x*y], B)
            sage: g.normalize_coordinates(); g
            Dynamical system of Projective Berkovich line over Cp(3), with base Rational Field
             induced by the map
              Defn: Defined on coordinates by sending (x : y) to (2*x : y)
        """
        self._system.normalize_coordinates()

    def conjugate(self, M, adjugate=False, new_ideal=None):
        r"""
        Conjugate this dynamical system by ``M``, i.e. `M^{-1} \circ f \circ M`.

        If possible the new map will be defined over the same space.
        Otherwise, will try to coerce to the base ring of ``M``.

        INPUT:

        - ``M`` -- a square invertible matrix

        - ``adjugate`` -- boolean (default: ``False``); also classically
          called adjoint, takes a square matrix ``M`` and finds the transpose
          of its cofactor matrix. Used for conjugation in place of inverse
          when specified ``'True'``. Functionality is the same in projective space.

        - ``new_ideal`` -- (optional) an ideal of the ``base_ring`` of ``M``.
          Used to specify an extension in the case where ``M`` is not defined
          over the same number field as this dynamical system.

        OUTPUT: a dynamical system

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, 2*y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: g.conjugate(Matrix([[1, 1], [0, 1]]))
            Dynamical system of Projective Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + (2 + O(3^20))*x*y : (2 + O(3^20))*y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_Berkovich([x^2 + y^2, y^2], ideal=5)
            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^2 + 1)
            sage: conj = Matrix([[1, a], [0, 1]])
            sage: f.conjugate(conj)
            Dynamical system of Projective Berkovich line over Cp(5), with base Number Field
            in a with defining polynomial z^2 + 1 induced by the map
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + (2*a)*x*y + (-a)*y^2 : y^2)

        We can use ``new_ideal`` to specify a new domain when
        the base ring of ``M`` and of this dynamical system are not the
        same::

            sage: # needs sage.rings.number_field
            sage: ideal = A.ideal(5).factor()[1][0]; ideal
            Fractional ideal (2*a + 1)
            sage: g = f.conjugate(conj, new_ideal=ideal)
            sage: g.domain().ideal()
            Fractional ideal (2*a + 1)
        """
        if self.domain().is_padic_base():
            return DynamicalSystem_Berkovich(self._system.conjugate(M, adjugate=adjugate))
        from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
        if not (isinstance(new_ideal, NumberFieldFractionalIdeal) or new_ideal is None or new_ideal in ZZ):
            raise TypeError('new_ideal must be an ideal of a number field, not %s' % new_ideal)
        new_system = self._system.conjugate(M, adjugate=adjugate)
        system_domain = new_system.domain()
        if new_ideal is None:
            if not system_domain.base_ring() is QQ:
                new_ideal = system_domain.base_ring().prime_above(self.domain().ideal())
            else:
                new_ideal = self.domain().ideal()
        return DynamicalSystem_Berkovich(new_system, ideal=new_ideal)

    def resultant(self, normalize=False):
        r"""
        Compute the resultant of the defining polynomials of
        this dynamical system.

        If ``normalize`` is ``True``, then first normalize the coordinate
        functions with :meth:`normalize_coordinates`.

        INPUT:

        - ``normalize`` -- boolean (default: ``False``)

        OUTPUT: an element of the base ring of this map

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_Berkovich([x^2 + y^2, y^2])
            sage: f.resultant()
            1 + O(3^20)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^3 + 20)
            sage: P.<x,y> = ProjectiveSpace(A, 1)
            sage: f = DynamicalSystem_Berkovich([2*x^2, x^2 + y^2], ideal=A.prime_above(2))
            sage: f.resultant()
            4
        """
        return self._system.resultant(normalize=normalize)

    def dehomogenize(self, n):
        """
        Return the map induced by the standard dehomogenization.

        The dehomogenization is done at the ``n[0]`` coordinate
        of the domain and the ``n[1]`` coordinate of the codomain.

        INPUT:

        - ``n`` -- tuple of nonnegative integers; if `n` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: a dynamical system on affine Berkovich space

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y + y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: g.dehomogenize(1)
            Dynamical system of Affine Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x) to
                    ((x^2 + 1 + O(3^20))/(x + 1 + O(3^20)))
        """
        new_system = self._system.dehomogenize(n)
        base_ring = self.domain().base_ring()
        ideal = self.domain().ideal()
        new_domain = Berkovich_Cp_Affine(base_ring, ideal)
        return DynamicalSystem_Berkovich_affine(new_system, new_domain)

    def __call__(self, x, type_3_pole_check=True):
        """
        Makes dynamical systems on Berkovich space over ``Cp`` callable.

        INPUT:

        - ``x`` -- a point of projective Berkovich space over ``Cp``

        - ``type_3_pole_check`` -- boolean (default: ``True``); WARNING:
          changing the value of ``type_3_pole_check`` can lead to mathematically
          incorrect answers. Only set to ``False`` if there are NO
          poles of the dynamical system in the disk corresponding
          to the type III point ``x``. See Examples.

        OUTPUT: a point of projective Berkovich space over ``Cp``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: g = DynamicalSystem_projective([x^2 + y^2, x*y])
            sage: G = DynamicalSystem_Berkovich(g)
            sage: B = G.domain()
            sage: Q3 = B(0, 1)
            sage: G(Q3)
            Type II point centered at (0 : 1 + O(3^20)) of radius 3^0

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)
            sage: H = DynamicalSystem_Berkovich([x*y^2, x^3 + 20*y^3])
            sage: B = H.domain()
            sage: Q4 = B(1/9, 1.5)
            sage: H(Q4, False)
            Type III point centered at (3^4 + 3^10 + 2*3^11 + 2*3^13 + 2*3^14 +
            2*3^15 + 3^17 + 2*3^18 + 2*3^19 + 3^20 + 3^21 + 3^22 + O(3^24) : 1 +
            O(3^20)) of radius 0.00205761316872428

        ALGORITHM:

        - For type II points, we use the approach outlined in Example
          7.37 of [Ben2019]_
        - For type III points, we use Proposition 7.6 of [Ben2019]_
        """
        if not isinstance(x.parent(), Berkovich_Cp_Projective):
            try:
                x = self.domain()(x)
            except (TypeError, ValueError):
                raise TypeError('action of dynamical system not defined on %s' % x.parent())
        if x.parent().is_padic_base() != self.domain().is_padic_base():
            raise ValueError('x was not backed by the same type of field as f')
        if x.prime() != self.domain().prime():
            raise ValueError('x and f are defined over Berkovich spaces over Cp for different p')
        if x.type_of_point() == 1:
            return self.domain()(self._system(x.center()))
        if x.type_of_point() == 4:
            raise NotImplementedError('action on Type IV points not implemented')
        f = self._system
        if x.type_of_point() == 2:
            if self.domain().is_number_field_base():
                ideal = self.domain().ideal()
                ring_of_integers = self.domain().base_ring().ring_of_integers()
            field = f.domain().base_ring()
            M = Matrix([[field(x.prime()**(-1 * x.power())), x.center()[0]], [field(0), field(1)]])
            F = list(f*M)
            R = field['z']
            S = f.domain().coordinate_ring()
            z = R.gen(0)
            dehomogenize_hom = S.hom([z, 1])
            for i in range(len(F)):
                F[i] = dehomogenize_hom(F[i])
            lcm = field(1)
            for poly in F:
                for i in poly:
                    if i != 0:
                        lcm = i.denominator().lcm(lcm)
            for i in range(len(F)):
                F[i] *= lcm
            gcd = [i for i in F[0] if i != 0][0]
            for poly in F:
                for i in poly:
                    if i != 0:
                        gcd = gcd * i * gcd.lcm(i).inverse_of_unit()
            for i in range(len(F)):
                F[i] *= gcd.inverse_of_unit()
            gcd = F[0].gcd(F[1])
            F[0] = F[0].quo_rem(gcd)[0]
            F[1] = F[1].quo_rem(gcd)[0]
            fraction = []
            for poly in F:
                new_poly = []
                for i in poly:
                    if self.domain().is_padic_base():
                        new_poly.append(i.residue())
                    else:
                        new_poly.append(ring_of_integers(i).mod(ideal))
                new_poly = R(new_poly)
                fraction.append(new_poly)
            gcd = fraction[0].gcd(fraction[1])
            num = fraction[0].quo_rem(gcd)[0]
            dem = fraction[1].quo_rem(gcd)[0]
            if dem.is_zero():
                f = DynamicalSystem_affine(F[0] / F[1]).homogenize(1)
                f = f.conjugate(Matrix([[0, 1], [1, 0]]))
                g = DynamicalSystem_Berkovich(f)
                return g(self.domain()(QQ(0), QQ(1))).involution_map()
            # if the reduction is not constant, the image is the Gauss point
            if not (num.is_constant() and dem.is_constant()):
                return self.domain()(QQ(0), QQ(1))
            if self.domain().is_padic_base():
                reduced_value = field(num * dem.inverse_of_unit()).lift_to_precision(field.precision_cap())
            else:
                reduced_value = field(num * dem.inverse_of_unit())
            new_num = F[0]-reduced_value*F[1]
            if self.domain().is_padic_base():
                power_of_p = min([i.valuation() for i in new_num])
            else:
                power_of_p = min([i.valuation(ideal) for i in new_num])
            inverse_map = field(x.prime()**power_of_p) * z + reduced_value
            if self.domain().is_padic_base():
                return self.domain()(inverse_map(0), (inverse_map(1) - inverse_map(0)).abs())
            else:
                val = (inverse_map(1) - inverse_map(0)).valuation(ideal)
                if val == Infinity:
                    return self.domain()(inverse_map(0), 0)
                return self.domain()(inverse_map(0), x.prime()**(-1 * val))
        # point is now type III, so we compute using Proposition 7.6 [of Benedetto]
        affine_system = f.dehomogenize(1)
        dem = affine_system.defining_polynomials()[0].denominator().univariate_polynomial()
        if type_3_pole_check:
            if self.domain().is_padic_base():
                factorization = [i[0] for i in dem.factor()]
                for factor in factorization:
                    if factor.degree() >= 2:
                        try:
                            factor_root_field = factor.root_field('a')
                            factor = factor.change_ring(factor_root_field)
                        except (TypeError, ValueError):
                            raise NotImplementedError('cannot check if poles lie in type III disk')
                    else:
                        factor_root_field = factor.base_ring()
                    center = factor_root_field(x.center()[0])
                    for pole in [i[0] for i in factor.roots()]:
                        if (center - pole).abs() <= x.radius():
                            raise NotImplementedError('image of type III point not implemented when poles in disk')
            else:
                dem_splitting_field, embedding = dem.splitting_field('a', True)
                poles = [i[0] for i in dem.roots(dem_splitting_field)]
                primes_above = dem_splitting_field.primes_above(self.domain().ideal())
                # check if all primes of the extension map the roots to outside
                # the disk corresponding to the type III point
                for prime in primes_above:
                    no_poles = True
                    for pole in poles:
                        valuation = (embedding(x.center()[0]) - pole).valuation(prime)
                        if valuation == Infinity:
                            no_poles = False
                            break
                        elif x.prime()**(-1 * valuation/prime.absolute_ramification_index()) <= x.radius():
                            no_poles = False
                            break
                    if not no_poles:
                        break
                if not no_poles:
                    raise NotImplementedError('image of type III not implemented when poles in disk')
        nth_derivative = f.dehomogenize(1).defining_polynomials()[0]
        variable = nth_derivative.parent().gens()[0]
        a = x.center()[0]
        Taylor_expansion = []
        from sage.arith.misc import factorial
        for i in range(f.degree() + 1):
            Taylor_expansion.append(nth_derivative(a) * 1/factorial(i))
            nth_derivative = nth_derivative.derivative(variable)
        r = x.radius()
        new_center = f(a)
        if self.domain().is_padic_base():
            new_radius = max([Taylor_expansion[i].abs()*r**i for i in range(1, len(Taylor_expansion))])
        else:
            if prime is None:
                prime = x.parent().ideal()
                dem_splitting_field = x.parent().base_ring()
            p = x.prime()
            new_radius = 0
            for i in range(1, len(Taylor_expansion)):
                valuation = dem_splitting_field(Taylor_expansion[i]).valuation(prime)
                new_radius = max(new_radius, p**(-valuation/prime.absolute_ramification_index())*r**i)
        return self.domain()(new_center, new_radius)


class DynamicalSystem_Berkovich_affine(DynamicalSystem_Berkovich):
    r"""
    A dynamical system of the affine Berkovich line over `\CC_p`.

    INPUT:

    - ``dynamical_system`` -- a :class:`DynamicalSystem_affine`
      of relative dimension 1

    - ``domain`` -- (optional) affine or projective Berkovich space
      over `\CC_p`. If the input to ``dynamical_system`` is
      not defined over `\QQ_p` or a finite extension, ``domain``
      must be specified.

    EXAMPLES:

    A dynamical system of the affine Berkovich line is
    induced by a dynamical system on `\QQ_p` or an extension
    of `\QQ_p`::

        sage: A.<x> = AffineSpace(Qp(5), 1)
        sage: f = DynamicalSystem_affine([(x^2 + 1)/x])
        sage: DynamicalSystem_Berkovich(f)
        Dynamical system of Affine Berkovich line over Cp(5) of precision 20 induced by the map
          Defn: Defined on coordinates by sending (x) to ((x^2 + 1 + O(5^20))/x)

    Dynamical system can be created from a morphism::

        sage: H = End(A)
        sage: phi = H([x + 3])
        sage: DynamicalSystem_Berkovich(phi)
        Dynamical system of Affine Berkovich line over Cp(5) of precision 20 induced by the map
          Defn: Defined on coordinates by sending (x) to (x + 3 + O(5^20))
    """
    @staticmethod
    def __classcall_private__(cls, dynamical_system, domain=None):
        """
        Return the appropriate dynamical system on affine Berkovich space over ``Cp``.

        EXAMPLES::

            sage: A.<x> = AffineSpace(Qp(3), 1)
            sage: from sage.dynamics.arithmetic_dynamics.berkovich_ds import DynamicalSystem_Berkovich_affine
            sage: DynamicalSystem_Berkovich_affine(DynamicalSystem_affine(x^2))
            Dynamical system of Affine Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x) to (x^2)
        """
        if not isinstance(dynamical_system, DynamicalSystem):
            if not isinstance(dynamical_system, DynamicalSystem_affine):
                dynamical_system = DynamicalSystem_projective(dynamical_system)
            else:
                raise TypeError('projective dynamical system passed to affine constructor')
        R = dynamical_system.base_ring()
        morphism_domain = dynamical_system.domain()
        if not isinstance(morphism_domain, AffineSpace_generic):
            raise TypeError('the domain of dynamical_system must be affine space, not %s' % morphism_domain)
        if morphism_domain.dimension_relative() != 1:
            raise ValueError('domain not relative dimension 1')
        if not isinstance(R, pAdicBaseGeneric):
            if domain is None:
                raise TypeError('dynamical system defined over %s, not padic, ' % morphism_domain.base_ring() +
                    'and domain was not specified')
            if not isinstance(domain, Berkovich_Cp_Affine):
                raise TypeError('domain was %s, not an affine Berkovich space over Cp' % domain)
        else:
            domain = Berkovich_Cp_Affine(morphism_domain.base_ring())
        return typecall(cls, dynamical_system, domain)

    def __init__(self, dynamical_system, domain):
        """
        Python constructor.

        EXAMPLES::

            sage: A.<x> = AffineSpace(Qp(3), 1)
            sage: from sage.dynamics.arithmetic_dynamics.berkovich_ds import DynamicalSystem_Berkovich_affine
            sage: DynamicalSystem_Berkovich_affine(DynamicalSystem_affine(x^3))
            Dynamical system of Affine Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x) to (x^3)
        """
        DynamicalSystem_Berkovich.__init__(self, dynamical_system, domain)

    def homogenize(self, n):
        """
        Return the homogenization of this dynamical system.

        For dynamical systems of Berkovich space, this is the dynamical
        system of projective Berkovich space induced by the homogenization of
        the dynamical system.

        INPUT:

        - ``n`` -- tuple of nonnegative integers. If `n` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: a dynamical system on projective Berkovich space

        EXAMPLES::

            sage: A.<x> = AffineSpace(Qp(3), 1)
            sage: f = DynamicalSystem_affine(1/x)
            sage: f = DynamicalSystem_Berkovich(f)
            sage: f.homogenize(1)
            Dynamical system of Projective Berkovich line over Cp(3) of precision 20
             induced by the map
              Defn: Defined on coordinates by sending (x0 : x1) to (x1 : x0)
        """
        new_system = self._system.homogenize(n)
        ideal = self.domain().ideal()
        base_space = new_system.domain()
        new_domain = Berkovich_Cp_Projective(base_space, ideal)
        return DynamicalSystem_Berkovich_projective(new_system, new_domain)

    def __call__(self, x):
        """
        Makes this dynamical system callable.

        EXAMPLES::

            sage: P.<x> = AffineSpace(Qp(3), 1)
            sage: f = DynamicalSystem_affine(x^2)
            sage: g = DynamicalSystem_Berkovich(f)
            sage: B = g.domain()
            sage: Q1 = B(2)
            sage: g(Q1)
            Type I point centered at 1 + 3 + O(3^20)
        """
        if not isinstance(x, Berkovich_Element_Cp_Affine):
            try:
                x = self.domain()(x)
            except (TypeError, ValueError):
                raise ValueError('action of dynamical system not defined on %s' % x)
        proj_system = self.homogenize(1)
        return proj_system(x.as_projective_point()).as_affine_point()
