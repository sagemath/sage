# sage_setup: distribution = sagemath-categories
"""
Ambient spaces
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.categories.commutative_rings import CommutativeRings
from sage.schemes.generic.scheme import Scheme


def is_AmbientSpace(x):
    """
    Return ``True`` if `x` is an ambient space.

    EXAMPLES::

        sage: from sage.schemes.generic.ambient_space import is_AmbientSpace
        sage: is_AmbientSpace(ProjectiveSpace(3, ZZ))
        doctest:warning...
        DeprecationWarning: The function is_AmbientSpace is deprecated; use 'isinstance(..., AmbientSpace)' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        True
        sage: is_AmbientSpace(AffineSpace(2, QQ))
        True
        sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
        sage: is_AmbientSpace(P.subscheme([x + y + z]))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_AmbientSpace is deprecated; use 'isinstance(..., AmbientSpace)' instead.")
    return isinstance(x, AmbientSpace)


class AmbientSpace(Scheme):
    """
    Base class for ambient spaces over a ring.

    INPUT:

    - ``n`` -- dimension

    - ``R`` -- ring
    """
    def __init__(self, n, R=ZZ):
        """
        TESTS::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: TestSuite(A).run() # not tested (abstract scheme with no elements?)
        """
        if R not in CommutativeRings():
            raise TypeError("R (={}) must be a commutative ring".format(R))
        if n < 0:
            raise ValueError("n (={}) must be nonnegative".format(n))
        self._dimension_relative = Integer(n)
        Scheme.__init__(self, R)

    #######################################################################
    # Derived classes must overload all of the following functions
    #######################################################################
    def _latex_(self):
        """
        TESTS::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._latex_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_(self):
        """
        TESTS::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._repr_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_generic_point(self, coords=None):
        """
        TESTS::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._repr_generic_point([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _latex_generic_point(self, coords=None):
        """
        TESTS::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._latex_generic_point([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme, or
        raise a :exc:`TypeError`.

        TESTS::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._check_satisfies_equations([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _validate(self, polynomials):
        """
        If ``polynomials`` is a tuple of valid polynomial functions on
        ``self``, return ``polynomials``, otherwise raise :exc:`TypeError`.

        INPUT:

        - ``polynomials`` -- tuple of polynomials in the coordinate ring of
          ``self``

        OUTPUT: tuple of polynomials in the coordinate ring of ``self``

        TESTS::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(3, ZZ)
            sage: A._validate((x + 1, 1))                                               # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: ambient spaces must override "_validate" method!
        """
        raise NotImplementedError('ambient spaces must override "_validate" '
                                  'method!')

    def change_ring(self, R):
        r"""
        Return an ambient space over ring `R` and otherwise the same as ``self``.

        INPUT:

        - ``R`` -- commutative ring

        OUTPUT: ambient space over ``R``

        .. NOTE::

            There is no need to have any relation between `R` and the base ring
            of  self, if you want to have such a relation, use
            ``self.base_extend(R)`` instead.

        TESTS::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5)
            sage: A.change_ring(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: ambient spaces must override "change_ring" method!
        """
        raise NotImplementedError(
                        'ambient spaces must override "change_ring" method!')

    #######################################################################
    # End overloads
    #######################################################################

    def is_projective(self):
        """
        Return whether this ambient space is projective n-space.

        EXAMPLES::

            sage: AffineSpace(3, QQ).is_projective()
            False
            sage: ProjectiveSpace(3, QQ).is_projective()
            True
        """
        # overloaded in the projective space derived class
        return False

    def base_extend(self, R):
        """
        Return the natural extension of ``self`` over ``R``.

        INPUT:

        - ``R`` -- a commutative ring, such that there is a natural map from
          the base ring of ``self`` to ``R``

        OUTPUT: an ambient space over ``R`` of the same structure as ``self``

        .. NOTE::

            A :exc:`ValueError` is raised if there is no such natural map.
            If you need to drop this condition, use ``self.change_ring(R)``.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: PQ = P.base_extend(QQ); PQ
            Projective Space of dimension 2 over Rational Field
            sage: PQ.base_extend(GF(5))
            Traceback (most recent call last):
            ...
            ValueError: no natural map from the base ring (=Rational Field)
            to R (=Finite Field of size 5)!
        """
        if R in CommutativeRings():
            if self.base_ring() == R:
                return self
            if not R.has_coerce_map_from(self.base_ring()):
                raise ValueError(
                    "no natural map from the base ring (=%s) to R (=%s)!"
                    % (self.base_ring(), R))
            return self.change_ring(R)
        else:
            raise NotImplementedError(
                        "extension of spaces over %s to %s is not implemented!"
                        % (self.base_ring(), R))

    def ambient_space(self):
        """
        Return the ambient space of the scheme self, in this case self
        itself.

        EXAMPLES::

            sage: P = ProjectiveSpace(4, ZZ)
            sage: P.ambient_space() is P
            True

            sage: A = AffineSpace(2, GF(3))
            sage: A.ambient_space()
            Affine Space of dimension 2 over Finite Field of size 3
        """
        return self

    def defining_polynomials(self):
        """
        Return the defining polynomials of the scheme ``self``.  Since
        ``self`` is an ambient space, this is an empty list.

        EXAMPLES::

            sage: ProjectiveSpace(2, QQ).defining_polynomials()
            ()
            sage: AffineSpace(0, ZZ).defining_polynomials()
            ()
        """
        return ()

    def identity_morphism(self):
        """
        Return the identity morphism.

        OUTPUT: the identity morphism of the scheme ``self``

        EXAMPLES::

            sage: A = AffineSpace(2, GF(3))
            sage: A.identity_morphism()
            Scheme endomorphism of Affine Space of dimension 2 over Finite Field of size 3
              Defn: Identity map

            sage: P = ProjectiveSpace(3, ZZ)
            sage: P.identity_morphism()
            Scheme endomorphism of Projective Space of dimension 3 over Integer Ring
              Defn: Identity map
        """
        from sage.schemes.generic.morphism import SchemeMorphism_polynomial_id
        return SchemeMorphism_polynomial_id(self)

    ######################################################################
    # Associated MPolynomial ring generators
    ######################################################################

    def gen(self, n=0):
        """
        Return the `n`-th generator of the coordinate ring of the
        scheme ``self``.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P.gen(1)
            y
        """
        return self.coordinate_ring().gen(n)

    def gens(self):
        """
        Return the generators of the coordinate ring of the scheme
        ``self``.

        EXAMPLES::

            sage: AffineSpace(0, QQ).gens()
            ()

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(5))
            sage: P.gens()
            (x, y, z)
        """
        return self.coordinate_ring().gens()

    def ngens(self):
        """
        Return the number of generators of the coordinate ring of the
        scheme ``self``.

        EXAMPLES::

            sage: AffineSpace(0, QQ).ngens()
            0

            sage: ProjectiveSpace(50, ZZ).ngens()
            51
        """
        return len(self.gens())

##     def assign_names(self, names=None):
##         """
##         EXAMPLES::

##             sage: A = AffineSpace(2, QQ, 'ab'); A
##             Affine Space of dimension 2 over Rational Field
##             sage: A.coordinate_ring()
##             Polynomial Ring in a, b over Rational Field
##             sage: A._assign_names('xy'); A.coordinate_ring()
##             Polynomial Ring in x, y over Rational Field
##         """
##         self.coordinate_ring()._assign_names(names)

    def dimension_absolute(self):
        """
        Return the absolute dimension of this scheme.

        EXAMPLES::

            sage: A2Q = AffineSpace(2, QQ)
            sage: A2Q.dimension_absolute()
            2
            sage: A2Q.dimension()
            2
            sage: A2Z = AffineSpace(2, ZZ)
            sage: A2Z.dimension_absolute()
            3
            sage: A2Z.dimension()
            3
        """
        base = self.base_scheme()
        if base.is_noetherian():
            return self.dimension_relative() + base.dimension()
        raise NotImplementedError("cannot compute the dimension of this scheme")

    dimension = dimension_absolute

    def dimension_relative(self):
        """
        Return the relative dimension of this scheme over its base.

        EXAMPLES::

            sage: A2Q = AffineSpace(2, QQ)
            sage: A2Q.dimension_relative()
            2
            sage: A2Z = AffineSpace(2, ZZ)
            sage: A2Z.dimension_relative()
            2
        """
        return self._dimension_relative
