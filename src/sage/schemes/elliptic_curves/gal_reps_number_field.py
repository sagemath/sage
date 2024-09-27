# sage.doctest: needs sage.rings.number_field
r"""
Galois representations for elliptic curves over number fields

This file contains the code to compute for which primes the Galois
representation attached to an elliptic curve (over an arbitrary number field)
is surjective. The functions in this file are called by the ``is_surjective``
and ``non_surjective`` methods of an elliptic curve over a number field.

EXAMPLES::

    sage: x = polygen(ZZ, 'x')
    sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
    sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
    sage: rho = E.galois_representation()
    sage: rho.is_surjective(29) # Cyclotomic character not surjective.
    False
    sage: rho.is_surjective(31) # See Section 5.10 of [Ser1972].
    True
    sage: rho.non_surjective()  # long time (4s on sage.math, 2014)
    [3, 5, 29]

    sage: E = EllipticCurve_from_j(1728).change_ring(K) # CM
    sage: E.galois_representation().non_surjective()  # long time (2s on sage.math, 2014)
    [0]

AUTHORS:

- Eric Larson (2012-05-28): initial version.
- Eric Larson (2014-08-13): added isogeny_bound function.
- John Cremona (2016, 2017): various efficiency improvements to _semistable_reducible_primes
- John Cremona (2017): implementation of Billerey's algorithm to find all reducible primes

REFERENCES:

- [Ser1972]_
- [Sut2012]_
"""
# ****************************************************************************
#       Copyright (C) 2012 Eric Larson <elarson3@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import legendre_symbol, primes
from sage.misc.functional import cyclotomic_polynomial
from sage.misc.lazy_import import lazy_import
from sage.modules.free_module import VectorSpace
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.set import Set
from sage.structure.sage_object import SageObject

lazy_import('sage.rings.number_field.number_field', 'NumberField')


class GaloisRepresentation(SageObject):
    r"""
    The compatible family of Galois representation
    attached to an elliptic curve over a number field.

    Given an elliptic curve `E` over a number field `K`
    and a rational prime number `p`, the `p^n`-torsion
    `E[p^n]` points of `E` is a representation of the
    absolute Galois group `G_K` of `K`. As `n` varies
    we obtain the Tate module `T_p E` which is
    a representation of `G_K` on a free `\ZZ_p`-module
    of rank `2`. As `p` varies the representations
    are compatible.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 + 1, 'a')
        sage: E = EllipticCurve('11a1').change_ring(K)
        sage: rho = E.galois_representation()
        sage: rho
        Compatible family of Galois representations associated to the Elliptic Curve
         defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20)
         over Number Field in a with defining polynomial x^2 + 1
    """

    def __init__(self, E):
        r"""
        See ``GaloisRepresentation`` for documentation.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: rho = E.galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve
             defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20)
             over Number Field in a with defining polynomial x^2 + 1
            sage: loads(rho.dumps()) == rho
            True
        """
        self.E = E

    def __repr__(self):
        r"""
        Return a string representation of the class.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: rho = E.galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve
             defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20)
             over Number Field in a with defining polynomial x^2 + 1

            sage: K.<a> = NumberField(x^2-x+1)
            sage: E = EllipticCurve([0,0,0,a,0])
            sage: E.galois_representation()
            Compatible family of Galois representations associated to the
             CM Elliptic Curve defined by y^2 = x^3 + a*x
             over Number Field in a with defining polynomial x^2 - x + 1
        """
        if self.E.has_cm():
            return "Compatible family of Galois representations associated to the CM " + repr(self.E)
        else:
            return "Compatible family of Galois representations associated to the " + repr(self.E)

    def __eq__(self,other):
        r"""
        Compare two Galois representations.

        We define two compatible families of representations attached
        to elliptic curves to be equal if the curves are isomorphic.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 + 1, 'a'); a = K.gen()
            sage: rho1 = EllipticCurve_from_j(1 + a).galois_representation()
            sage: rho2 = EllipticCurve_from_j(2 + a).galois_representation()
            sage: rho1 == rho1
            True
            sage: rho1 == rho2
            False
            sage: rho1 == 42
            False
        """
        if type(self) is not type(other):
            return False
        return self.E.is_isomorphic(other.E)

    def elliptic_curve(self):
        r"""
        Return the elliptic curve associated to this representation.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 + 1, 'a'); a = K.gen()
            sage: E = EllipticCurve_from_j(a)
            sage: rho = E.galois_representation()
            sage: rho.elliptic_curve() == E
            True
        """
        return self.E

    def non_surjective(self, A=100):
        r"""
        Return a list of primes `p` including all primes for which the mod-`p`
        representation might not be surjective.

        INPUT:

        - ``A`` -- integer; a bound on the number of traces of Frobenius to use
          while trying to prove surjectivity

        OUTPUT:

        A list of primes where mod-`p` representation is very likely not
        surjective. At any prime not in this list, the representation is
        definitely surjective. If `E` has CM, the list ``[0]`` is returned.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.non_surjective() # See Section 5.10 of [Ser1972].  # long time
            [3, 5, 29]
            sage: K = NumberField(x**2 + 3, 'a'); a = K.gen()
            sage: E = EllipticCurve([0, -1, 1, -10, -20]).change_ring(K) # X_0(11)
            sage: rho = E.galois_representation()
            sage: rho.non_surjective()  # long time (4s on sage.math, 2014)
            [3, 5]
            sage: K = NumberField(x**2 + 1, 'a'); a = K.gen()
            sage: E = EllipticCurve_from_j(1728).change_ring(K) # CM
            sage: rho = E.galois_representation()
            sage: rho.non_surjective()
            [0]
            sage: K = NumberField(x**2 - 5, 'a'); a = K.gen()
            sage: E = EllipticCurve_from_j(146329141248*a - 327201914880) # CM
            sage: rho = E.galois_representation()
            sage: rho.non_surjective() # long time (3s on sage.math, 2014)
            [0]

        TESTS:

        An example which failed until fixed at :issue:`19229`::

            sage: K.<a> = NumberField(x^2-x+1)
            sage: E = EllipticCurve([a+1,1,1,0,0])
            sage: rho = E.galois_representation()
            sage: rho.non_surjective()
            [2, 3]
        """
        if self.E.has_cm():
            return [0]
        return _non_surjective(self.E, A)

    def is_surjective(self, p, A=100):
        r"""
        Return ``True`` if the mod-p representation is (provably)
        surjective onto `Aut(E[p]) = GL_2(\GF{p})`.  Return
        ``False`` if it is (probably) not.

        INPUT:

        - ``p`` -- prime number

        - ``A`` -- integer; a bound on the number of traces of Frobenius to use
          while trying to prove surjectivity

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.is_surjective(29) # Cyclotomic character not surjective.
            False
            sage: rho.is_surjective(7) # See Section 5.10 of [Ser1972].
            True

        If `E` is defined over `\QQ`, then the exceptional primes for `E_{/K}`
        are the same as the exceptional primes for `E`, except for those primes
        that are ramified in `K/\QQ` or are less than `[K:\QQ]`::

            sage: K = NumberField(x**2 + 11, 'a')
            sage: E = EllipticCurve([2, 14])
            sage: rhoQQ = E.galois_representation()
            sage: rhoK = E.change_ring(K).galois_representation()
            sage: rhoQQ.is_surjective(2) == rhoK.is_surjective(2)
            False
            sage: rhoQQ.is_surjective(3) == rhoK.is_surjective(3)
            True
            sage: rhoQQ.is_surjective(5) == rhoK.is_surjective(5)
            True

        For CM curves, the mod-p representation is never surjective::

            sage: K.<a> = NumberField(x^2 - x + 1)
            sage: E = EllipticCurve([0,0,0,0,a])
            sage: E.has_cm()
            True
            sage: rho = E.galois_representation()
            sage: any(rho.is_surjective(p) for p in [2,3,5,7])
            False
        """
        if self.E.has_cm():
            return False
        return not _exceptionals(self.E, [p], A)

    def isogeny_bound(self, A=100):
        r"""
        Return a list of primes `p` including all primes for which
        the image of the mod-`p` representation is contained in a
        Borel.

        .. NOTE::

            For the actual list of primes `p` at which the
            representation is reducible see :meth:`reducible_primes()`.

        INPUT:

        - ``A`` -- integer; a bound on the number of traces of Frobenius to
          use while trying to prove the mod-`p` representation is not contained
          in a Borel

        OUTPUT:

        A list of primes which contains (but may not be equal to) all `p` for
        which the image of the mod-`p` representation is contained in a Borel
        subgroup.  At any prime not in this list, the image is definitely not
        contained in a Borel. If E has `CM` defined over `K`, the list ``[0]``
        is returned.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound() # See Section 5.10 of [Ser1972].  # long time
            [3, 5]
            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve_from_j(K(1728))             # CM over K
            sage: E.galois_representation().isogeny_bound()
            [0]
            sage: E = EllipticCurve_from_j(K(0))                # CM NOT over K
            sage: E.galois_representation().isogeny_bound()     # long time
            [2, 3]
            sage: E = EllipticCurve_from_j(K(2268945/128))      # c.f. [Sut2012]
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound()  # No 7-isogeny, but...   # long time
            [7]

        For curves with rational CM, there are infinitely many primes
        `p` for which the mod-`p` representation is reducible, and [0]
        is returned::

            sage: K.<a> = NumberField(x^2 - x + 1)
            sage: E = EllipticCurve([0,0,0,0,a])
            sage: E.has_rational_cm()
            True
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound()
            [0]

        An example (an elliptic curve with everywhere good reduction
        over an imaginary quadratic field with quite large
        discriminant), which failed until fixed at :issue:`21776`::

            sage: K.<a> = NumberField(x^2 - x + 112941801)
            sage: E = EllipticCurve([a+1,a-1,a,-23163076*a + 266044005933275,57560769602038*a - 836483958630700313803])
            sage: E.conductor().norm()
            1
            sage: GR = E.galois_representation()
            sage: GR.isogeny_bound()
            []
        """
        if self.E.has_rational_cm():
            return [0]

        E = _over_numberfield(self.E)
        K = E.base_field()

        char = lambda P: P.smallest_integer() # cheaper than constructing the residue field

        # semistable reducible primes (we are now not in the CM case)
        bad_primes = _semistable_reducible_primes(E)

        # primes of additive reduction
        bad_primesK = (K.ideal(E.c4()) + K.ideal(E.discriminant())).prime_factors()
        bad_primes += [char(P) for P in bad_primesK]

        # ramified primes
        bad_primes += K.absolute_discriminant().prime_factors()

        # remove repeats:
        bad_primes = list(Set(bad_primes))

        return Frobenius_filter(E, bad_primes, A)

    def reducible_primes(self):
        r"""
        Return a list of primes `p` for which the mod-`p`
        representation is reducible, or [0] for CM curves.

        OUTPUT:

        A list of those primes `p` for which the mod-`p` representation is
        contained in a Borel subgroup, i.e. is reducible.  If E has CM
        *defined over K*, the list ``[0]`` is returned (in this case the
        representation is reducible for infinitely many primes).

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.reducible_primes() # See Section 5.10 of [Ser1972].  # long time
            [3, 5]

            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve_from_j(K(1728))             # CM over K
            sage: E.galois_representation().reducible_primes()
            [0]
            sage: E = EllipticCurve_from_j(K(0))                # CM but NOT over K
            sage: E.galois_representation().reducible_primes()  # long time
            [2, 3]
            sage: E = EllipticCurve_from_j(K(2268945/128))      # c.f. [Sut2012]
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound()  # No 7-isogeny, but...   # long time
            [7]
            sage: rho.reducible_primes()                        # long time
            []

        For curves with rational CM, there are infinitely many primes
        `p` for which the mod-`p` representation is reducible, and [0]
        is returned::

            sage: K.<a> = NumberField(x^2 - x + 1)
            sage: E = EllipticCurve([0,0,0,0,a])
            sage: E.has_rational_cm()
            True
            sage: rho = E.galois_representation()
            sage: rho.reducible_primes()
            [0]
        """
        if self.E.has_rational_cm():
            return [0]

        return [l for l in self.isogeny_bound() if self.E.isogenies_prime_degree(l)]


def _non_surjective(E, patience=100):
    r"""
    Return a list of primes `p` including all primes for which the mod-`p`
    representation might not be surjective.

    INPUT:

    - ``E`` -- :class:`EllipticCurve` (over a number field)

    - ``A`` -- integer; a bound on the number of traces of Frobenius to use
      while trying to prove surjectivity

    OUTPUT:

    A list of primes where mod-`p` representation is very likely not
    surjective. At any prime not in this list, the representation is
    definitely surjective. If E has CM, a :exc:`ValueError` is raised.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._non_surjective(E) # See Section 5.10 of [Ser1972].  # long time
        [3, 5, 29]
        sage: E = EllipticCurve_from_j(1728).change_ring(K) # CM
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._non_surjective(E)
        Traceback (most recent call last):
        ...
        ValueError: The curve E should not have CM.
    """
    if E.has_cm():
        raise ValueError("The curve E should not have CM.")

    E = _over_numberfield(E)
    K = E.base_field()

    exceptional_primes = [2, 3, 5, 7, 11, 13, 17, 19]
    # The possible primes l unramified in K/QQ for which the image of the mod l
    # Galois representation could be contained in an exceptional subgroup.

    # Find the places of additive reduction.
    SA = []
    for P, n in E.conductor().factor():
        if n > 1:
            SA.append(P)
    # TODO: All we really require is a list of primes that include all
    # primes at which E has additive reduction. Perhaps we can speed
    # up things by doing something less time-consuming here that produces
    # a list with some extra terms? (Of course, the longer this list is,
    # the slower the rest of the computation is, so it is not clear that
    # this would help...)

    char = lambda P: P.smallest_integer() # cheaper than constructing the residue field

    bad_primes = exceptional_primes
    bad_primes += [char(P) for P in SA]
    bad_primes += K.discriminant().prime_factors()
    bad_primes += _semistable_reducible_primes(E)
    bad_primes += _possible_normalizers(E, SA)

    bad_primes = list(Set(bad_primes))

    return _exceptionals(E, bad_primes, patience)


def Frobenius_filter(E, L, patience=100):
    r"""
    Determine which primes in L might have an image contained in a
    Borel subgroup, by checking of traces of Frobenius.

    .. NOTE::

        This function will sometimes return primes for which the image
        is not contained in a Borel subgroup.  This issue cannot always
        be fixed by increasing patience as it may be a result of a
        failure of a local-global principle for isogenies.

    INPUT:

    - ``E`` -- :class:`EllipticCurve` over a number field

    - ``L`` -- list of prime numbers

    - ``patience`` -- positive integer (default: 100) bounding the number of
      traces of Frobenius to use while trying to prove irreducibility

    OUTPUT:

    The list of all primes `\ell` in L for which the mod `\ell` image might be
    contained in a Borel subgroup of `GL_2(\mathbf{F}_{\ell})`.

    EXAMPLES::

        sage: E = EllipticCurve('11a1') # has a 5-isogeny
        sage: sage.schemes.elliptic_curves.gal_reps_number_field.Frobenius_filter(E,primes(40))  # long time
        [5]

    Example to show that the output may contain primes where the
    representation is in fact reducible.  Over `\QQ` the following is
    essentially the unique such example by [Sut2012]_::

        sage: E = EllipticCurve_from_j(2268945/128)
        sage: sage.schemes.elliptic_curves.gal_reps_number_field.Frobenius_filter(E, [7, 11])  # long time
        [7]

    This curve does possess a 7-isogeny modulo every prime of good
    reduction, but has no rational 7-isogeny::

        sage: E.isogenies_prime_degree(7)
        []

    A number field example::

        sage: K.<i> = QuadraticField(-1)
        sage: E = EllipticCurve([1+i, -i, i, -399-240*i,  2627+2869*i])
        sage: sage.schemes.elliptic_curves.gal_reps_number_field.Frobenius_filter(E, primes(20))  # long time
        [2, 3]

    Here the curve really does possess isogenies of degrees 2 and 3::

        sage: [len(E.isogenies_prime_degree(l)) for l in [2,3]]
        [1, 1]
    """
    E = _over_numberfield(E).global_integral_model()
    K = E.base_field()

    L = list(set(L)) # Remove duplicates from L and makes a copy for output
    L.sort()

    include_2 = False
    if 2 in L: # c.f. Section 5.3(a) of [Ser1972].
        L.remove(2)
        include_2 = not E.division_polynomial(2).is_irreducible()

    K_is_Q = (K == QQ)
    from sage.arith.misc import primes
    from sage.rings.infinity import infinity

    def primes_iter():
        for p in primes(start=2, stop=infinity):
            if K_is_Q:
                if E.has_good_reduction(p):
                    yield ZZ.ideal(p)
            else:
                for P in K.primes_above(p):
                    if E.has_good_reduction(P):
                        yield P
    numP = 0
    for P in primes_iter():
        if not L or numP == patience:  # stop if no primes are left, or patience is exhausted
            break

        numP += 1

        # Discard any l for which the Frobenius polynomial at P is
        # irreducible modulo l

        disc = E.reduction(P).frobenius_polynomial().discriminant()

        L = [l for l in L if legendre_symbol(disc,l) != -1]

        #print("After using {} primes P, {}  primes l remain".format(numP,len(L)))

    if include_2:
        L = [2] + L
    return L


def _exceptionals(E, L, patience=1000):
    r"""
    Determine which primes in L are exceptional for E.

    This is done using Proposition 19 of Section 2.8 of Serre's
    *Propriétés Galoisiennes des Points d'Ordre Fini des Courbes Elliptiques*
    [Ser1972]_.

    INPUT:

    - ``E`` -- :class:`EllipticCurve` over a number field

    - ``L`` -- list of prime numbers

    - ``patience`` -- integer; a bound on the number of traces of Frobenius to
      use while trying to prove surjectivity

    OUTPUT:

    The list of all primes l in L for which the mod l image might fail to be
    surjective.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._exceptionals(E, [29, 31])
        [29]

    For CM curves an error is raised::

        sage: E = EllipticCurve_from_j(1728).change_ring(K) # CM
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._exceptionals(E,[2,3,5])
        Traceback (most recent call last):
        ...
        ValueError: The curve E should not have CM.
    """
    if E.has_cm():
        raise ValueError("The curve E should not have CM.")

    E = _over_numberfield(E)
    K = E.base_field()

    output = []

    L = list(set(L)) # Remove duplicates from L.

    for l in L:
        if l == 2: # c.f. Section 5.3(a) of [Ser1972].
            if (E.j_invariant() - 1728).is_square():
                output.append(2)
            elif not E.division_polynomial(2).is_irreducible():
                output.append(2)

        elif l == 3: # c.f. Section 5.3(b) of [Ser1972].
            if K(-3).is_square():
                output.append(3)
            elif not (K['x'].gen()**3 - E.j_invariant()).is_irreducible():
                output.append(3)
            elif not E.division_polynomial(3).is_irreducible():
                output.append(3)

        elif (K.discriminant() % l) == 0:
            if not K['x'](cyclotomic_polynomial(l)).is_irreducible():
                # I.E. if the action on lth roots of unity is not surjective
                # (We want this since as a Galois module, \wedge^2 E[l]
                # is isomorphic to the lth roots of unity.)
                output.append(l)

    for l in output:
        L.remove(l)
    if 2 in L:
        L.remove(2)
    if 3 in L:
        L.remove(3)

    # If the image is not surjective, then it is contained in one of the
    # maximal subgroups. So, we start by creating a dictionary between primes
    # l in L and possible maximal subgroups in which the mod l image could
    # be contained. This information is stored as a triple whose elements
    # are True/False according to whether the mod l image could be contained
    # in:
    #      0. A Borel or normalizer of split Cartan subgroup.
    #      1. A nonsplit Cartan subgroup or its normalizer.
    #      2. An exceptional subgroup of GL_2.

    D = {}
    for l in L:
        D[l] = [True, True, True]

    for P in deg_one_primes_iter(K):
        try:
            trace = E.change_ring(P.residue_field()).trace_of_frobenius()
        except ArithmeticError: # Bad reduction at P.
            continue

        patience -= 1

        determinant = P.norm()
        discriminant = trace**2 - 4 * determinant

        unexc = [] # Primes we discover are unexceptional go here.

        for l in D:
            tr = GF(l)(trace)
            det = GF(l)(determinant)
            disc = GF(l)(discriminant)

            if tr == 0:
                # I.E. if Frob_P could be contained in the normalizer of
                # a Cartan subgroup, but not in the Cartan subgroup.
                continue

            if disc == 0:
                # I.E. If the matrix might be non-diagonalizable over F_{p^2}.
                continue

            if legendre_symbol(disc, l) == 1:
                # If the matrix is diagonalizable over F_p, it can't be
                # contained in a non-split Cartan subgroup. Since we've
                # gotten rid of the case where it is contained in the
                # of a nonsplit Cartan subgroup but not the Cartan subgroup,
                D[l][1] = False
            else:
                # If the matrix is not diagonalizable over F_p, it can't
                # be contained Borel subgroup.
                D[l][0] = False

            if det != 0:  # c.f. [Ser1972], Section 2.8, Prop. 19
                u = trace**2 / det
                if u not in (1, 2, 4) and u**2 - 3 * u + 1 != 0:
                    D[l][2] = False

            if D[l] == [False, False, False]:
                unexc.append(l)

        for l in unexc:
            D.pop(l)
        unexc = []

        if (not D) or (patience == 0):
            break

    output.extend(D)
    output.sort()
    return output


def _over_numberfield(E):
    r"""
    Return `E`, defined over an absolute ``NumberField``.

    This is necessary since if `E` is defined over `\QQ`, then we
    cannot use Sage commands available for number fields, and if the
    base field is relative then other problems result.

    INPUT:

    - ``E`` -- :class:`EllipticCurve` over a number field

    OUTPUT:

    - If `E` is defined over an absolute number field, returns `E`.

    - If `E` is defined over a relative number field, returns `E` defined over the absolute base field.

    - If `E` is defined over `\QQ`, returns `E` defined over the number field `\QQ`.

    EXAMPLES::

        sage: E = EllipticCurve([1, 2])
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._over_numberfield(E)
        Elliptic Curve defined by y^2 = x^3 + x + 2 over Number Field in a with defining polynomial x
    """
    K = E.base_field()

    if K == QQ:
        from sage.rings.polynomial.polynomial_ring import polygen
        K = NumberField(polygen(QQ), 'a')
    else:
        K = K.absolute_field('a')
    return E.change_ring(K)


def deg_one_primes_iter(K, principal_only=False):
    r"""
    Return an iterator over degree 1 primes of ``K``.

    INPUT:

    - ``K`` -- a number field
    - ``principal_only`` -- boolean; if ``True``, only yield principal primes

    OUTPUT:

    An iterator over degree 1 primes of `K` up to the given norm,
    optionally yielding only principal primes.

    EXAMPLES::

        sage: K.<a> = QuadraticField(-5)
        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import deg_one_primes_iter
        sage: it = deg_one_primes_iter(K)
        sage: [next(it) for _ in range(6)]
        [Fractional ideal (2, a + 1),
         Fractional ideal (3, a + 1),
         Fractional ideal (3, a + 2),
         Fractional ideal (a),
         Fractional ideal (7, a + 3),
         Fractional ideal (7, a + 4)]
        sage: it = deg_one_primes_iter(K, True)
        sage: [next(it) for _ in range(6)]
        [Fractional ideal (a),
         Fractional ideal (-2*a + 3),
         Fractional ideal (2*a + 3),
         Fractional ideal (a + 6),
         Fractional ideal (a - 6),
         Fractional ideal (-3*a + 4)]
    """
    # imaginary quadratic fields have no principal primes of norm < disc / 4
    start = K.discriminant().abs() // 4 if principal_only and K.signature() == (0,1) else 2

    K_is_Q = (K == QQ)

    for p in primes(start=start, stop=Infinity):
        if K_is_Q:
            yield ZZ.ideal(p)
        else:
            for P in K.primes_above(p, degree=1):
                if not principal_only or P.is_principal():
                    yield P


def _semistable_reducible_primes(E, verbose=False):
    r"""Find a list containing all semistable primes l unramified in K/QQ
    for which the Galois image for E could be reducible.

    INPUT:

    - ``E`` -- :class:`EllipticCurve` over a number field

    OUTPUT:

    A list of primes, which contains all primes `l` unramified in
    `K/\mathbb{QQ}`, such that `E` is semistable at all primes lying
    over `l`, and the Galois image at `l` is reducible. If `E` has CM
    defined over its ground field, a :exc:`ValueError` is raised.

    EXAMPLES::

        sage: E = EllipticCurve([0, -1, 1, -10, -20]) # X_0(11)
        sage: 5 in sage.schemes.elliptic_curves.gal_reps_number_field._semistable_reducible_primes(E)
        True

    This example, over a quintic field with Galois group `S_5`, took a
    very long time before :issue:`22343`::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^5 - 6*x^3 + 8*x - 1)
        sage: E = EllipticCurve(K, [a^3 - 2*a, a^4 - 2*a^3 - 4*a^2 + 6*a + 1,
        ....:                       a + 1, -a^3 + a + 1, -a])
        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import _semistable_reducible_primes
        sage: _semistable_reducible_primes(E)
        [2, 5, 53, 1117]
    """
    if verbose:
        print("In _semistable_reducible_primes with E={}".format(E.ainvs()))
    K = E.base_field()
    d = K.degree()

    deg_one_primes = deg_one_primes_iter(K, principal_only=True)

    bad_primes = set()  # This will store the output.

    # We find two primes (of distinct residue characteristics) which are
    # of degree 1, unramified in K/Q, and at which E has good reduction.
    # Each of these primes will give us a nontrivial divisibility constraint
    # on the exceptional primes l. For both of these primes P, we precompute
    # a generator and the characteristic polynomial of Frob_P^12.

    precomp = []
    last_p = 0 # The residue characteristic of the most recent prime.

    while len(precomp) < 2:
        P = next(deg_one_primes)
        p = P.norm()
        if p != last_p and (d == 1 or P.ramification_index() == 1) and E.has_good_reduction(P):
            precomp.append(P)
            last_p = p

    Px, Py = precomp
    x, y = (PP.gens_reduced()[0] for PP in precomp)
    EmodPx = E.reduction(Px) if d > 1 else E.reduction(x)
    EmodPy = E.reduction(Py) if d > 1 else E.reduction(y)
    fxpol = EmodPx.frobenius_polynomial()
    fypol = EmodPy.frobenius_polynomial()
    fx12pol = fxpol.adams_operator_on_roots(12) # roots are 12th powers of those of fxpol
    fy12pol = fypol.adams_operator_on_roots(12)
    px = x.norm() if d > 1 else x
    py = y.norm() if d > 1 else x
    Zx = fxpol.parent()
    xpol = x.charpoly() if d > 1 else Zx([-x,1])
    ypol = y.charpoly() if d > 1 else Zx([-y,1])

    if verbose:
        print("Finished precomp, x={} (p={}), y={} (p={})".format(x,px,y,py))

    for w in range(1 + d // 2):
        if verbose:
            print("w = {}".format(w))
        gx = xpol.symmetric_power(w).adams_operator_on_roots(12).resultant(fx12pol)
        gy = ypol.symmetric_power(w).adams_operator_on_roots(12).resultant(fy12pol)
        if verbose:
            print("computed gx and gy")

        gxn = Integer(gx.absolute_norm()) if d > 1 else gx
        gyn = Integer(gy.absolute_norm()) if d > 1 else gy
        gxyn = gxn.gcd(gyn)
        if gxyn:
            xprimes = gxyn.prime_factors()
            if verbose:
                print("adding prime factors {} of {} to {}".format(xprimes, gxyn, sorted(bad_primes)))
            bad_primes.update(xprimes)
            if verbose:
                print("...done, bad_primes now {}".format(sorted(bad_primes)))
            continue
        else:
            if verbose:
                print("gx and gy both 0!")

        ## It is possible that our curve has CM. ##

        # Our character must be of the form Nm^K_F for an imaginary
        # quadratic subfield F of K (which is the CM field if E has CM).

        # Note that this can only happen when d is even, w=d/2, and K
        # contains (or the Galois closure of K contains?) the
        # imaginary quadratic field F = Q(sqrt(a)) which is the
        # splitting field of both fx12pol and fy12pol.  We compute a
        # and relativise K over F:

        a = fx12pol.discriminant().squarefree_part()

        # Construct a field isomorphic to K but a relative extension over QQ(sqrt(a)).

        # See #19229: the names given here, which are not used, should
        # not be the name of the generator of the base field.

        rootsa = K(a).sqrt(all=True) # otherwise if a is not a square the
                                     # returned result is in the symbolic ring!
        try:
            roota = rootsa[0]
        except IndexError:
            raise RuntimeError("error in _semistable_reducible_primes: K={} does not contain sqrt({})".format(K,a))
        K_rel = K.relativize(roota, ['name1','name2'])
        iso = K_rel.structure()[1] # an isomorphism from K to K_rel

        ## We try again to find a nontrivial divisibility condition. ##

        div = 0
        patience = 5 * K.absolute_degree()
        # Number of Frobenius elements to check before suspecting that E
        # has CM and computing the set of CM j-invariants of K to check.
        # TODO: Is this the best value for this parameter?

        while div == 0 and patience > 0:
            P = next(deg_one_primes) # a prime of K not K_rel
            while E.has_bad_reduction(P):
                P = next(deg_one_primes)

            if verbose:
                print("trying P = {}...".format(P))
            EmodP = E.reduction(P)
            fpol = EmodP.frobenius_polynomial()
            if verbose:
                print("...good reduction, frobenius poly = {}".format(fpol))
            x = iso(P.gens_reduced()[0]).relative_norm()
            xpol = x.charpoly().adams_operator_on_roots(12)
            div2 = Integer(xpol.resultant(fpol.adams_operator_on_roots(12)) // x.norm()**12)
            if div2:
                div = div2.isqrt()
                assert div2 == div**2
                if verbose:
                    print("...div = {}".format(div))
            else:
                if verbose:
                    print("...div = 0, continuing")
                patience -= 1

        if patience == 0:
            # We suspect that E has CM, so we check:
            if E.has_cm():
                raise ValueError("In _semistable_reducible_primes, the curve E should not have CM.")

        assert div != 0
        # We found our divisibility constraint.

        xprimes = div.prime_factors()
        if verbose:
            print("...adding prime factors {} of {} to {}...".format(xprimes,div, sorted(bad_primes)))
        bad_primes.update(xprimes)
        if verbose:
            print("...done, bad_primes now {}".format(sorted(bad_primes)))

    L = sorted(bad_primes)
    return L


def _possible_normalizers(E, SA):
    r"""Find a list containing all primes `l` such that the Galois image at `l`
    is contained in the normalizer of a Cartan subgroup, such that the
    corresponding quadratic character is ramified only at the given primes.

    INPUT:

    - ``E`` -- :class:`EllipticCurve` over a number field K

    - ``SA`` -- list of primes of K

    OUTPUT:

    - a list of primes, which contains all primes `l` such that the Galois
      image at `l` is contained in the normalizer of a Cartan subgroup, such
      that the corresponding quadratic character is ramified only at primes in SA.

    - if `E` has geometric CM that is not defined over its ground field, a
      :exc:`ValueError` is raised.

    EXAMPLES::

        sage: E = EllipticCurve([0,0,0,-56,4848])
        sage: 5 in sage.schemes.elliptic_curves.gal_reps_number_field._possible_normalizers(E, [ZZ.ideal(2)])
        True

    For CM curves, an error is raised::

        sage: K.<i> = QuadraticField(-1)
        sage: E = EllipticCurve_from_j(1728).change_ring(K) # CM
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._possible_normalizers(E, [])
        Traceback (most recent call last):
        ...
        ValueError: The curve E should not have CM.
    """
    if E.has_cm():
        raise ValueError("The curve E should not have CM.")

    E = _over_numberfield(E)
    K = E.base_field()
    SA = [K.ideal(I.gens()) for I in SA]

    selmer_gens = K.selmer_generators(SA, 2) # Generators of the selmer group.

    if not selmer_gens:
        return []

    V = VectorSpace(GF(2), len(selmer_gens))
    # We think of this as the character group of the selmer group.

    traces_list = []
    W = V.zero_subspace()

    deg_one_primes = deg_one_primes_iter(K)

    while W.dimension() < V.dimension() - 1:
        P = next(deg_one_primes)

        k = P.residue_field()

        defines_valid_character = True
        # A prime P defines a quadratic residue character
        # on the Selmer group. This variable will be set
        # to zero if any elements of the selmer group are
        # zero mod P (i.e. the character is ramified).

        splitting_vector = [] # This will be the values of this
        # character on the generators of the Selmer group.

        for a in selmer_gens:
            abar = k(a)
            if abar == 0:
                # Ramification.
                defines_valid_character = False
                break

            if abar.is_square():
                splitting_vector.append(GF(2)(0))
            else:
                splitting_vector.append(GF(2)(1))

        if not defines_valid_character:
            continue

        if splitting_vector in W:
            continue

        try:
            Etilde = E.change_ring(k)
        except ArithmeticError: # Bad reduction.
            continue

        tr = Etilde.trace_of_frobenius()

        if tr == 0:
            continue

        traces_list.append(tr)

        W = W + V.span([splitting_vector])

    bad_primes = set()

    for i in traces_list:
        for p in i.prime_factors():
            bad_primes.add(p)

    # We find the unique vector v in V orthogonal to W:
    v = W.matrix().transpose().kernel().basis()[0]

    # We find the element a of the selmer group corresponding to v:
    a = 1
    for i in range(len(selmer_gens)):
        if v[i] == 1:
            a *= selmer_gens[i]

    # Since we've already included the above bad primes, we can assume
    # that the quadratic character corresponding to the exceptional primes
    # we're looking for is given by mapping into Gal(K[\sqrt{a}]/K).

    patience = 5 * K.degree()
    # Number of Frobenius elements to check before suspecting that E
    # has CM and computing the set of CM j-invariants of K to check.
    # TODO: Is this the best value for this parameter?

    while True:
        P = next(deg_one_primes)

        k = P.residue_field()

        if not k(a).is_square():
            try:
                tr = E.change_ring(k).trace_of_frobenius()
            except ArithmeticError: # Bad reduction.
                continue

            if tr == 0:
                patience -= 1

            else:
                for p in tr.prime_factors():
                    bad_primes.add(p)

                bad_primes = sorted(bad_primes)
                return bad_primes

#
# Code for Billerey's algorithm to find reducible primes
#
# See "Critères d'irréductibilité pour les représentations des courbes
# elliptiques", Nicolas Billerey, https://arxiv.org/abs/0908.1084
#


def Billerey_P_l(E, l):
    r"""
    Return Billerey's `P_l^*` as defined in [Bil2011]_, equation (9).

    INPUT:

    - ``E`` -- an elliptic curve over a number field `K`, given by a
      global integral model

    - ``l`` -- a rational prime

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import Billerey_P_l
        sage: [Billerey_P_l(E,l) for l in primes(10)]
        [x^2 + 8143*x + 16777216,
         x^2 + 451358*x + 282429536481,
         x^4 - 664299076*x^3 + 205155493652343750*x^2 - 39595310449600219726562500*x + 3552713678800500929355621337890625,
         x^4 - 207302404*x^3 - 377423798538689366394*x^2 - 39715249826471656586987520004*x + 36703368217294125441230211032033660188801]
    """
    K = E.base_field()
    qq = K.primes_above(l)
    # if len(qq) == K.degree():
    #     return None
    from sage.rings.polynomial.polynomial_ring import polygen
    from operator import mul
    P = polygen(ZZ)-1
    for q in qq:
        e = K(l).valuation(q)
        P = P.composed_op(E.reduction(q).frobenius_polynomial().adams_operator_on_roots(12*e), mul, monic=True)
    return P


def Billerey_B_l(E,l,B=0):
    r"""
    Return Billerey's `B_l`, adapted from the definition in [Bil2011]_, after (9).

    INPUT:

    - ``E`` -- an elliptic curve over a number field `K`, given by a
      global integral model

    - ``l`` -- integer; a rational prime

    - ``B`` -- integer; 0 or LCM of previous `B_l` the prime-to-B part of this
      `B_l` is ignored

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import Billerey_B_l
        sage: [Billerey_B_l(E,l) for l in primes(15)]
        [1123077552537600,
         227279663773903886745600,
         0,
         0,
         269247154818492941287713746693964214802283882086400,
         0]
    """
    d = E.base_field().absolute_degree()
    P = Billerey_P_l(E, l)
    if P is None:
        return ZZ.zero()
    # We compute the factors one at a time since if any is 0 we quit:
    B_l = ZZ(1)
    for k in range(1 + d // 2):
        factor = ZZ(P(l**(12*k)))
        if factor:
            B_l *= factor.gcd(B)
        else:
            return ZZ(0)
    return B_l


def Billerey_R_q(E, q, B=0):
    r"""
    Return Billerey's `R_q`, adapted from the definition in [Bil2011]_, Theorem 2.8.

    INPUT:

    - ``E`` -- an elliptic curve over a number field `K`, given by a
      global integral model

    - ``q`` -- a prime ideal of `K`

    - ``B`` -- integer; 0 or LCM of previous `R_q` the prime-to-B part of this
      `R_q` is ignored

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import Billerey_R_q
        sage: [Billerey_R_q(E,K.prime_above(l)) for l in primes(10)]
        [1123077552537600,
        227279663773903886745600,
        51956919562116960000000000000000,
        252485933820556361829926400000000]
    """
    K = E.base_field()
    d = K.absolute_degree()
    h = K.class_number()
    P = E.reduction(q).frobenius_polynomial().adams_operator_on_roots(12*h)
    Q = ((q**h).gens_reduced()[0]).absolute_minpoly().adams_operator_on_roots(12)

    # We compute the factors one at a time since if any is 0 we quit:
    R_q = ZZ(1)
    for k in range(1 + d // 2):
        # the following would be in QQ if we did not coerce
        factor = ZZ(P.resultant(Q.compose_power(k)))
        if factor:
            R_q *= factor.gcd(B)
        else:
            return ZZ.zero()
    return R_q


def Billerey_B_bound(E, max_l=200, num_l=8, small_prime_bound=0, debug=False):
    r"""
    Compute Billerey's bound `B`.

    We compute `B_l` for `l` up to ``max_l`` (at most) until ``num_l``
    nonzero values are found (at most).  Return the list of primes
    dividing all `B_l` computed, excluding those dividing 6 or
    ramified or of bad reduction or less than small_prime_bound.  If
    no nonzero values are found return ``[0]``.

    INPUT:

    - ``E`` -- an elliptic curve over a number field `K`, given by a
      global integral model

    - ``max_l`` -- integer (default: 200); maximum size of primes l to check

    - ``num_l`` integer (default: 8); maximum number of primes l to check

    - ``small_prime_bound`` -- integer (default: 0); remove primes less
      than this from the output

    - ``debug`` -- boolean (default: ``False``); if ``True`` prints details

    .. NOTE::

        The purpose of the small_prime_bound is that it is faster to
        deal with these using the local test; by ignoring them here,
        we enable the algorithm to terminate sooner when there are no
        large reducible primes, which is always the case in practice.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import Billerey_B_bound
        sage: Billerey_B_bound(E)
        [5]

    If we do not use enough primes `l`, extraneous primes will be
    included which are not reducible primes::

        sage: Billerey_B_bound(E, num_l=6)
        [5, 7]

    Similarly if we do not use large enough primes `l`::

        sage: Billerey_B_bound(E, max_l=50, num_l=8)
        [5, 7]
        sage: Billerey_B_bound(E, max_l=100, num_l=8)
        [5]

    This curve does have a rational 5-isogeny::

        sage: len(E.isogenies_prime_degree(5))
        1
    """
    if debug:
        print("Computing B-bound for {} with max_l={}, num_l={}".format(E.ainvs(),max_l,num_l) + " (ignoring primes under {})".format(small_prime_bound) if small_prime_bound else "")
    B = ZZ.zero()
    ells = []
    K = E.base_field()
    DK = K.discriminant()
    ED = E.discriminant().norm()
    B0 = ZZ(6*DK*ED)

    def remove_primes(B):
        B1 = B.prime_to_m_part(B0)
        for p in primes(small_prime_bound):
            B1 = B1.prime_to_m_part(p)
        return B1
    ll = primes(5,max_l) # iterator
    while B != 1 and len(ells) < num_l:
        try:
            l = next(ll)
            while B0.valuation(l):
                l = next(ll)
        except StopIteration:
            break
        if debug:
            print("..trying l={}".format(l))
        b = Billerey_B_l(E,l,B)
        if b:
            if debug:
                print("..ok, B_l = {}".format(b))
            if B:
                B = B.gcd(b)
            else:
                B = remove_primes(b)
            ells.append(l)
            if debug:
                print("..so far, B = {} using l in {}".format(B,ells))
        else:
            if debug:
                print("..B_l=0 for l={}".format(l))

    if B:
        res = [p for p,e in B.factor()]
        if debug:
            print("..returning {}".format(res))
        return res
    # or we failed to find any nonzero values...
    if debug:
        print("..failed to find a bound")
    return [0]


def Billerey_R_bound(E, max_l=200, num_l=8, small_prime_bound=None, debug=False):
    r"""
    Compute Billerey's bound `R`.

    We compute `R_q` for `q` dividing primes `\ell` up to ``max_l``
    (at most) until ``num_l`` nonzero values are found (at most).
    Return the list of primes dividing all ``R_q`` computed, excluding
    those dividing 6 or ramified or of bad reduction or less than
    small_prime_bound.  If no nonzero values are found return ``[0]``.

    INPUT:

    - ``E`` -- an elliptic curve over a number field `K`, given by a
      global integral model

    - ``max_l`` -- integer (default: 200); maximum size of rational primes
      l for which the primes q above l are checked

    - ``num_l`` -- integer (default: 8); maximum number of rational primes
      l for which the primes q above l are checked

    - ``small_prime_bound`` -- integer (default: 0); remove primes less
      than this from the output

    - ``debug`` -- boolean (default: ``False``); if ``True`` prints details

    .. NOTE::

        The purpose of the ``small_prime_bound`` is that it is faster to
        deal with these using the local test; by ignoring them here,
        we enable the algorithm to terminate sooner when there are no
        large reducible primes, which is always the case in practice.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import Billerey_R_bound
        sage: Billerey_R_bound(E)
        [5]

    We may get no bound at all if we do not use enough primes::

        sage: Billerey_R_bound(E, max_l=2, debug=False)
        [0]

    Or we may get a bound but not a good one if we do not use enough primes::

        sage: Billerey_R_bound(E, num_l=1, debug=False)
        [5, 17, 67, 157]

    In this case two primes is enough to restrict the set of possible
    reducible primes to just `\{5\}`.  This curve does have a rational 5-isogeny::

        sage: Billerey_R_bound(E, num_l=2, debug=False)
        [5]
        sage: len(E.isogenies_prime_degree(5))
        1
    """
    if debug:
        print("Computing R-bound for {} with max_l={}, num_l={}".format(E.ainvs(),max_l,num_l) + " (ignoring primes under {})".format(small_prime_bound) if small_prime_bound else "")
    B = ZZ.zero()
    ells = []
    K = E.base_field()
    DK = K.discriminant()
    ED = E.discriminant().norm()
    B0 = ZZ(6*DK*ED)

    def remove_primes(B):
        B1 = B.prime_to_m_part(B0)
        for p in primes(small_prime_bound):
            B1 = B1.prime_to_m_part(p)
        return B1
    ll = primes(5, max_l) # iterator
    while len(ells) < num_l and B != 1:
        try:
            l = next(ll)
            while B0.valuation(l):
                l = next(ll)
        except StopIteration:
            break
        q = K.prime_above(l)
        if debug:
            print("..trying q={} above l={}".format(q,l))
        b = Billerey_R_q(E,q,B)
        if b:
            if debug:
                print("..ok, R_q = {}, type={}".format(b,type(b)))
            if B:
                B = B.gcd(b)
            else:
                B = remove_primes(b)
            ells.append(l)
            if debug:
                print("..so far, B = {} using l in {}".format(B,ells))

    if B:
        res = B.support()
        if debug:
            print("..returning {}".format(res))
        return res
    # or we failed to find any nonzero values...
    if debug:
        print("..failed to find a bound")
    return [0]


def reducible_primes_Billerey(E, num_l=None, max_l=None, verbose=False):
    r"""
    Return a finite set of primes `\ell` containing all those for which
    `E` has a `K`-rational ell-isogeny, where `K` is the base field of
    `E`: i.e., the mod-`\ell` representation is irreducible for all
    `\ell` outside the set returned.

    INPUT:

    - ``E`` -- an elliptic curve defined over a number field `K`

    - ``max_l`` -- integer or ``None`` (default); the maximum prime
      `\ell` to use for the B-bound and R-bound.  If ``None``, a
      default value will be used.

    - ``num_l`` -- integer or ``None`` (default); the number of primes
      `\ell` to use for the B-bound and R-bound.  If ``None``, a
      default value will be used.

    .. NOTE::

        If ``E`` has CM then [0] is returned.  In this case use the
        function
        sage.schemes.elliptic_curves.isogeny_class.possible_isogeny_degrees

    We first compute Billeray's B_bound using at most ``num_l`` primes
    of size up to ``max_l``.  If that fails we compute Billeray's
    R_bound using at most ``num_q`` primes of size up to ``max_q``.

    Provided that one of these methods succeeds in producing a finite
    list of primes we check these using a local condition, and finally
    test that the primes returned actually are reducible.  Otherwise
    we return ``[0]``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import reducible_primes_Billerey
        sage: x = polygen(ZZ, 'x')
        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: reducible_primes_Billerey(E)  # long time
        [3, 5]
        sage: K = NumberField(x**2 + 1, 'a')
        sage: E = EllipticCurve_from_j(K(1728)) # CM over K
        sage: reducible_primes_Billerey(E)  # long time
        [0]
        sage: E = EllipticCurve_from_j(K(0)) # CM but NOT over K
        sage: reducible_primes_Billerey(E)  # long time
        [2, 3]

    An example where a prime is not reducible but passes the test::

        sage: E = EllipticCurve_from_j(K(2268945/128)).global_minimal_model() # c.f. [Sut2012]
        sage: reducible_primes_Billerey(E)  # long time
        [7]

    TESTS:

    Test that this function works with non-integral models (see :issue:`34174`)::

        sage: K.<a> = QuadraticField(4569)
        sage: j = 46969655/32768
        sage: E = EllipticCurve(j=K(j))
        sage: EK = E.change_ring(K)
        sage: C = EK.isogeny_class(minimal_models=False)    # long time
        sage: len(C)                                        # long time
        4
    """
    #verbose=True
    if verbose:
        print("E = {}, finding reducible primes using Billerey's algorithm".format(E.ainvs()))

    # Set parameters to default values if not given:
    if max_l is None:
        max_l = 200
    if num_l is None:
        num_l = 8

    K = E.base_field()
    DK = K.discriminant()

    # We replace E by an integral model if necessary, since this
    # function and the helper functions need this:
    E1 = E.global_integral_model()
    ED = E1.discriminant().norm()
    B0 = ZZ(6*DK*ED).prime_divisors()

    # Billeray's algorithm will be faster if we tell it to ignore
    # small primes; these can be tested using the naive algorithm.

    if verbose:
        print("First doing naive test of primes up to {}...".format(max_l))

    max_small_prime = 200
    OK_small_primes = reducible_primes_naive(E1, max_l=max_small_prime, num_P=200, verbose=verbose)
    if verbose:
        print("Naive test of primes up to {} returns {}.".format(max_small_prime, OK_small_primes))

    B1 = Billerey_B_bound(E1, max_l, num_l, max_small_prime, verbose)
    if B1 == [0]:
        if verbose:
            print("...  B_bound ineffective using max_l={}, moving on to R-bound".format(max_l))

        B1 = Billerey_R_bound(E1,max_l, num_l, max_small_prime, verbose)
        if B1 == [0]:
            if verbose:
                print("... R_bound ineffective using max_l={}",format(max_l))
            return [0]
        if verbose:
            print("... R_bound = {}".format(B1))
    else:
        if verbose:
            print("... B_bound = {}".format(B1))
    B = sorted(set(B0 + B1 + OK_small_primes))
    if verbose:
        print("... combined bound = {}".format(B))

    num_p = 100
    B = Frobenius_filter(E1, B, num_p)
    if verbose:
        print("... after Frobenius filter = {}".format(B))
    return B


def reducible_primes_naive(E, max_l=None, num_P=None, verbose=False):
    r"""
    Return locally reducible primes `\ell` up to ``max_l``.

    The list of primes `\ell` returned consists of all those up to
    ``max_l`` such that `E` mod `P` has an `\ell`-isogeny, where `K`
    is the base field of `E`, for ``num_P`` primes `P` of `K`.  In
    most cases `E` then has a `K`-rational `\ell`-isogeny, but there
    are rare exceptions.

    INPUT:

    - ``E`` -- an elliptic curve defined over a number field `K`

    - ``max_l`` -- integer or ``None`` (default); the maximum prime
      `\ell` to test

    - ``num_P`` -- integer or ``None`` (default); the number of primes
      `P` of `K` to use in testing each `\ell`

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import reducible_primes_naive
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^4 - 5*x^2 + 3)
        sage: E = EllipticCurve(K, [a^2 - 2, -a^2 + 3, a^2 - 2, -50*a^2 + 35, 95*a^2 - 67])
        sage: reducible_primes_naive(E,num_P=10)
        [2, 5, 53, 173, 197, 241, 293, 317, 409, 557, 601, 653, 677, 769, 773, 797]
        sage: reducible_primes_naive(E,num_P=15)
        [2, 5, 197, 557, 653, 769]
        sage: reducible_primes_naive(E,num_P=20)
        [2, 5]
        sage: reducible_primes_naive(E)  # long time
        [2, 5]
        sage: [phi.degree() for phi in E.isogenies_prime_degree()]  # long time
        [2, 2, 2, 5]
    """
    if max_l is None:
        max_l = 1000
    if num_P is None:
        num_P = 100
    if verbose:
        print("E = {}, finding reducible primes up to {} using Frobenius filter with {} primes".format(E.ainvs(), max_l, num_P))

    B = Frobenius_filter(E, primes(max_l), num_P)
    if verbose:
        print("... returning {}".format(B))
    return B
