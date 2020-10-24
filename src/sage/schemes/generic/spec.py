"""
The Spec functor

AUTHORS:

- William Stein (2006): initial implementation

- Peter Bruin (2014): rewrite Spec as a functor
"""

# ******************************************************************************
#  Copyright (C) 2006 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ******************************************************************************

from sage.categories.functor import Functor
from sage.rings.integer_ring import ZZ
from sage.schemes.generic.scheme import AffineScheme
from sage.structure.unique_representation import UniqueRepresentation


def Spec(R, S=None):
    r"""
    Apply the Spec functor to `R`.

    INPUT:

    - ``R`` -- either a commutative ring or a ring homomorphism

    - ``S`` -- a commutative ring (optional), the base ring

    OUTPUT:

    - ``AffineScheme`` -- the affine scheme `\mathrm{Spec}(R)`

    EXAMPLES::

        sage: Spec(QQ)
        Spectrum of Rational Field
        sage: Spec(PolynomialRing(QQ, 'x'))
        Spectrum of Univariate Polynomial Ring in x over Rational Field
        sage: Spec(PolynomialRing(QQ, 'x', 3))
        Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        sage: X = Spec(PolynomialRing(GF(49,'a'), 3, 'x')); X
        Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Finite Field in a of size 7^2
        sage: TestSuite(X).run()

    Applying ``Spec`` twice to the same ring gives identical output
    (see :trac:`17008`)::

        sage: A = Spec(ZZ); B = Spec(ZZ)
        sage: A is B
        True

    A ``TypeError`` is raised if the input is not a commutative ring::

        sage: Spec(5)
        Traceback (most recent call last):
        ...
        TypeError: x (=5) is not in Category of commutative rings
        sage: Spec(FreeAlgebra(QQ,2, 'x'))
        Traceback (most recent call last):
        ...
        TypeError: x (=Free Algebra on 2 generators (x0, x1) over Rational Field) is not in Category of commutative rings

    TESTS::

        sage: X = Spec(ZZ)
        sage: X
        Spectrum of Integer Ring
        sage: X.base_scheme()
        Spectrum of Integer Ring
        sage: X.base_ring()
        Integer Ring
        sage: X.dimension()
        1
        sage: Spec(QQ,QQ).base_scheme()
        Spectrum of Rational Field
        sage: Spec(RDF,QQ).base_scheme()
        Spectrum of Rational Field
    """
    return SpecFunctor(S)(R)


class SpecFunctor(Functor, UniqueRepresentation):
    """
    The Spec functor.
    """
    def __init__(self, base_ring=None):
        """
        EXAMPLES::

            sage: from sage.schemes.generic.spec import SpecFunctor
            sage: SpecFunctor()
            Spec functor from Category of commutative rings to Category of schemes
            sage: SpecFunctor(QQ)
            Spec functor from Category of commutative rings to
             Category of schemes over Rational Field
        """
        from sage.categories.all import CommutativeRings, Schemes

        if base_ring is None:
            domain = CommutativeRings()
            codomain = Schemes()
        elif base_ring in CommutativeRings():
            # We would like to use CommutativeAlgebras(base_ring) as
            # the domain; we use CommutativeRings() instead because
            # currently many algebras are not yet considered to be in
            # CommutativeAlgebras(base_ring) by the category framework.
            domain = CommutativeRings()
            codomain = Schemes(AffineScheme(base_ring))
        else:
            raise TypeError('base (= {}) must be a commutative ring'.format(base_ring))
        self._base_ring = base_ring
        super(SpecFunctor, self).__init__(domain, codomain)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.schemes.generic.spec import SpecFunctor
            sage: SpecFunctor(QQ)
            Spec functor from Category of commutative rings to
             Category of schemes over Rational Field
        """
        return 'Spec functor from {} to {}'.format(self.domain(), self.codomain())

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.schemes.generic.spec import SpecFunctor
            sage: latex(SpecFunctor())
            \mathrm{Spec}\colon \mathbf{CommutativeRings} \longrightarrow \mathbf{Schemes}
        """
        return r'\mathrm{{Spec}}\colon {} \longrightarrow {}'.format(
               self.domain()._latex_(), self.codomain()._latex_())

    def _apply_functor(self, A):
        """
        Apply the Spec functor to the commutative ring ``A``.

        EXAMPLES::

            sage: from sage.schemes.generic.spec import SpecFunctor
            sage: F = SpecFunctor()
            sage: F(RR) # indirect doctest
            Spectrum of Real Field with 53 bits of precision
        """
        # The second argument of AffineScheme defaults to None.
        # However, AffineScheme has unique representation, so there is
        # a difference between calling it with or without explicitly
        # giving this argument.
        if self._base_ring is None:
            return AffineScheme(A)
        return AffineScheme(A, self._base_ring)

    def _apply_functor_to_morphism(self, f):
        """
        Apply the Spec functor to the ring homomorphism ``f``.

        EXAMPLES::

            sage: from sage.schemes.generic.spec import SpecFunctor
            sage: F = SpecFunctor(GF(7))
            sage: A.<x, y> = GF(7)[]
            sage: B.<t> = GF(7)[]
            sage: f = A.hom((t^2, t^3))
            sage: Spec(f) # indirect doctest
            Affine Scheme morphism:
              From: Spectrum of Univariate Polynomial Ring in t over Finite Field of size 7
              To:   Spectrum of Multivariate Polynomial Ring in x, y over Finite Field of size 7
              Defn: Ring morphism:
                     From: Multivariate Polynomial Ring in x, y over Finite Field of size 7
                     To:   Univariate Polynomial Ring in t over Finite Field of size 7
                     Defn: x |--> t^2
                           y |--> t^3
        """
        A = f.domain()
        B = f.codomain()
        return self(B).hom(f, self(A))


SpecZ = Spec(ZZ)


# Compatibility with older versions of this module

from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.schemes.generic.spec', 'Spec', AffineScheme)
