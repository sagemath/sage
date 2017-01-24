# -*- coding: utf-8 -*-
r"""
Galois representations for elliptic curves over number fields.

This file contains the code to compute for which primes the Galois
representation attached to an elliptic curve (over an arbitrary
number field) is surjective. The functions in this file are called by
the `is_surjective` and `non_surjective` methods of an elliptic curve
over a number field.

EXAMPLES::

    sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
    sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
    sage: rho = E.galois_representation()
    sage: rho.is_surjective(29) # Cyclotomic character not surjective.
    False
    sage: rho.is_surjective(31) # See Section 5.10 of [Serre72].
    True
    sage: rho.non_surjective()  # long time (4s on sage.math, 2014)
    [3, 5, 29]

    sage: E = EllipticCurve_from_j(1728).change_ring(K) # CM
    sage: E.galois_representation().non_surjective()  # long time (2s on sage.math, 2014)
    [0]

AUTHORS:

- Eric Larson (2012-05-28): initial version.
- Eric Larson (2014-08-13): added isogeny_bound function.

REFERENCES:

.. [Serre72] Serre. Propriétés galoisiennes des points d'ordre fini des
    courbes elliptiques. Inventiones mathematicae, 1972.

.. [Sutherland12] Sutherland. A local-global principle for rational
    isogenies of prime degree. Journal de Theorie des Nombres de Bordeaux,
    2012.

"""

#*****************************************************************************
#       Copyright (C) 2012 Eric Larson <elarson3@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from six.moves import range

from sage.structure.sage_object import SageObject
from sage.rings.number_field.number_field import NumberField, QuadraticField
from sage.schemes.elliptic_curves.cm import cm_j_invariants
from sage.rings.rational_field import QQ
from sage.modules.free_module import VectorSpace
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.misc.functional import cyclotomic_polynomial
from sage.arith.all import legendre_symbol
from sage.sets.set import Set

class GaloisRepresentation(SageObject):
    r"""
    The compatible family of Galois representation
    attached to an elliptic curve over a number field.

    Given an elliptic curve `E` over a number field `K`
    and a rational prime number `p`, the `p^n`-torsion
    `E[p^n]` points of `E` is a representation of the
    absolute Galois group `G_K` of `K`. As `n` varies
    we obtain the Tate module `T_p E` which is a
    a representation of `G_K` on a free `\ZZ_p`-module
    of rank `2`. As `p` varies the representations
    are compatible.

    EXAMPLES::

        sage: K = NumberField(x**2 + 1, 'a')
        sage: E = EllipticCurve('11a1').change_ring(K)
        sage: rho = E.galois_representation()
        sage: rho
        Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in a with defining polynomial x^2 + 1
    """


    def __init__(self, E):
        r"""
        See ``GaloisRepresentation`` for documentation.

        EXAMPLES::

            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: rho = E.galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in a with defining polynomial x^2 + 1
            sage: loads(rho.dumps()) == rho
            True
        """
        self.E = E

    def __repr__(self):
        r"""
        Return a string representation of the class.

        EXAMPLES::

            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: rho = E.galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in a with defining polynomial x^2 + 1

            sage: K.<a> = NumberField(x^2-x+1)
            sage: E = EllipticCurve([0,0,0,a,0])
            sage: E.galois_representation()
            Compatible family of Galois representations associated to the CM Elliptic Curve defined by y^2 = x^3 + a*x over Number Field in a with defining polynomial x^2 - x + 1
        """
        if self.E.has_cm():
            return "Compatible family of Galois representations associated to the CM " + repr(self.E)
        else:
            return "Compatible family of Galois representations associated to the " + repr(self.E)


    def __eq__(self,other):
        r"""
        Compares two Galois representations.

        We define two compatible families of representations attached
        to elliptic curves to be equal if the curves are isomorphic.

        EXAMPLES::

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

        * ``A`` - int (a bound on the number of traces of Frobenius to use
                     while trying to prove surjectivity).

        OUTPUT:

        - ``list`` - A list of primes where mod-`p` representation is
          very likely not surjective. At any prime not in this list,
          the representation is definitely surjective. If `E` has CM,
          the list [0] is returned.

        EXAMPLES::

            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.non_surjective() # See Section 5.10 of [Serre72].
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

        An example which failed until fixed at :trac:`19229`::

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
        surjective onto `Aut(E[p]) = GL_2(\mathbb{F}_p)`.  Return
        ``False`` if it is (probably) not.

        INPUT:

        * ``p`` - int - a prime number.

        * ``A`` - int - a bound on the number of traces of Frobenius to use
                     while trying to prove surjectivity.

        EXAMPLES::

            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.is_surjective(29) # Cyclotomic character not surjective.
            False
            sage: rho.is_surjective(7) # See Section 5.10 of [Serre72].
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

            sage: K.<a> = NumberField(x^2-x+1)
            sage: E = EllipticCurve([0,0,0,0,a])
            sage: E.has_cm()
            True
            sage: rho = E.galois_representation()
            sage: any(rho.is_surjective(p) for p in [2,3,5,7])
            False
        """
        if self.E.has_cm():
            return False
        return (_exceptionals(self.E, [p], A) == [])

    def isogeny_bound(self, A=100):
        r"""
        Returns a list of primes `p` including all primes for which
        the image of the mod-`p` representation is contained in a
        Borel.

        .. NOTE::

           For the actual list of primes `p` at which the
           representation is reducible see :meth:`reducible_primes()`.

        INPUT:

        - ``A`` - int (a bound on the number of traces of Frobenius to
                     use while trying to prove the mod-`p`
                     representation is not contained in a Borel).

        OUTPUT:

        - ``list`` - A list of primes which contains (but may not be
          equal to) all `p` for which the image of the mod-`p`
          representation is contained in a Borel subgroup.  At any
          prime not in this list, the image is definitely not
          contained in a Borel. If E has `CM` defined over `K`, the list
          [0] is returned.

        EXAMPLES::

            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound() # See Section 5.10 of [Serre72].
            [3, 5]
            sage: K = NumberField(x**2 + 1, 'a')
            sage: EllipticCurve_from_j(K(1728)).galois_representation().isogeny_bound() # CM over K
            [0]
            sage: EllipticCurve_from_j(K(0)).galois_representation().isogeny_bound() # CM NOT over K
            [2, 3]
            sage: E = EllipticCurve_from_j(K(2268945/128)) # c.f. [Sutherland12]
            sage: E.galois_representation().isogeny_bound() # No 7-isogeny, but...
            [7]

        For curves with rational CM, there are infinitely many primes
        `p` for which the mod-`p` representation is reducible, and [0]
        is returned::

            sage: K.<a> = NumberField(x^2-x+1)
            sage: E = EllipticCurve([0,0,0,0,a])
            sage: E.has_rational_cm()
            True
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound()
            [0]

        An example (an elliptic curve with everywhere good reduction
        over an imaginary quadratic field with quite large
        discriminant), which failed until fixed at :trac:`21776`::

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

        return _maybe_borels(E, bad_primes, A)

    def reducible_primes(self):
        r"""
        Returns a list of primes `p` for which the mod-`p`
        representation is reducible, or [0] for CM curves.

        OUTPUT:

        - ``list`` - A list of those primes `p` for which the mod-`p`
          representation is contained in a Borel subgroup, i.e. is
          reducible.  If E has CM *defined over K*, the list [0] is
          returned (in this case the representation is reducible for
          infinitely many primes).

        EXAMPLES::

            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound() # See Section 5.10 of [Serre72].
            [3, 5]
            sage: rho.reducible_primes()
            [3, 5]

            sage: K = NumberField(x**2 + 1, 'a')
            sage: EllipticCurve_from_j(K(1728)).galois_representation().isogeny_bound() # CM over K
            [0]
            sage: EllipticCurve_from_j(K(0)).galois_representation().reducible_primes() # CM but NOT over K
            [2, 3]
            sage: E = EllipticCurve_from_j(K(2268945/128)) # c.f. [Sutherland12]
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound() # ... but there is no 7-isogeny ...
            [7]
            sage: rho.reducible_primes()
            []

        For curves with rational CM, there are infinitely many primes
        `p` for which the mod-`p` representation is reducible, and [0]
        is returned::

            sage: K.<a> = NumberField(x^2-x+1)
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
    Returns a list of primes `p` including all primes for which the mod-`p`
    representation might not be surjective.

    INPUT:

    - ``E`` - EllipticCurve (over a number field).

    - ``A`` - int (a bound on the number of traces of Frobenius to use
                 while trying to prove surjectivity).

    OUTPUT:

    - ``list`` - A list of primes where mod-`p` representation is very likely
      not surjective. At any prime not in this list, the representation is
      definitely surjective. If E has CM, a ValueError is raised.

    EXAMPLES::

        sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
        sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._non_surjective(E) # See Section 5.10 of [Serre72].
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


def _maybe_borels(E, L, patience=100):
    r"""
    Determine which primes in L might have an image contained in a
    Borel subgroup, using straight-forward checking of traces of
    Frobenius.

    .. NOTE:

       This function will sometimes return primes for which the image
       is not contained in a Borel subgroup.  This issue cannot always
       be fixed by increasing patience as it may be a result of a
       failure of a local-global principle for isogenies.

    INPUT:

    - ``E`` - EllipticCurve - over a number field.

    - ``L`` - list - a list of prime numbers.

    - ``patience`` - int (a positive integer bounding the number of
                          traces of Frobenius to use while trying to
                          prove irreducibility).

    OUTPUT: list - The list of all primes `\ell` in L for which the
                   mod `\ell` image might be contained in a Borel
                   subgroup of `GL_2(\mathbf{F}_{\ell})`.

    EXAMPLES::

        sage: E = EllipticCurve('11a1') # has a 5-isogeny
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._maybe_borels(E,primes(40))
        [5]

    Example to show that the output may contain primes where the
    representation is in fact reducible.  Over `\QQ` the following is
    essentially the unique such example by [Sutherland12]_::

        sage: E = EllipticCurve_from_j(2268945/128)
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._maybe_borels(E, [7, 11])
        [7]

    This curve does possess a 7-isogeny modulo every prime of good
    reduction, but has no rational 7-isogeny::

        sage: E.isogenies_prime_degree(7)
        []

    A number field example:

        sage: K.<i> = QuadraticField(-1)
        sage: E = EllipticCurve([1+i, -i, i, -399-240*i,  2627+2869*i])
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._maybe_borels(E, primes(20))
        [2, 3]

    Here the curve really does possess isognies of degrees 2 and 3::

        sage: [len(E.isogenies_prime_degree(l)) for l in [2,3]]
        [1, 1]

    """
    E = _over_numberfield(E)
    K = E.base_field()

    L = sorted(set(L)) # Remove duplicates from L and makes a copy for output

    include_2 = False
    if 2 in L: # c.f. Section 5.3(a) of [Serre72].
        L.remove(2)
        include_2 = not E.division_polynomial(2).is_irreducible()

    for P in deg_one_primes_iter(K):
        if not (L and patience): # stop if no primes are left, or
                                 # patience is exhausted
            break

        patience -= 1

        # Check whether the Frobenius polynomial at P is irreducible
        # modulo each l, dropping l from the list if so.

        try:
            trace = E.change_ring(P.residue_field()).trace_of_frobenius()
        except ArithmeticError: # Bad reduction at P.
            continue

        determinant = P.norm()
        discriminant = trace**2 - 4 * determinant

        for l in L:
            if legendre_symbol(discriminant,l)==-1:
                L.remove(l)

    if include_2:
        L = [2] + L
    return L

def _exceptionals(E, L, patience=1000):
    r"""
    Determine which primes in L are exceptional for E, using Proposition 19
    of Section 2.8 of Serre's ``Proprietes Galoisiennes des Points d'Ordre
    Fini des Courbes Elliptiques'' [Serre72].

    INPUT:

    - ``E`` - EllipticCurve - over a number field.

    - ``L`` - list - a list of prime numbers.

    - ``patience`` - int (a bound on the number of traces of Frobenius to
                          use while trying to prove surjectivity).

    OUTPUT: list - The list of all primes l in L for which the mod l image
                   might fail to be surjective.

    EXAMPLES::

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
        if l == 2: # c.f. Section 5.3(a) of [Serre72].
            if (E.j_invariant() - 1728).is_square():
                output.append(2)
            elif not E.division_polynomial(2).is_irreducible():
                output.append(2)

        elif l == 3: # c.f. Section 5.3(b) of [Serre72].
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

            if det != 0: # c.f. [Serre72], Section 2.8, Prop. 19
                u = trace**2 / det
                if u not in (1, 2, 4) and u**2 - 3 * u + 1 != 0:
                    D[l][2] = False


            if D[l] == [False, False, False]:
                unexc.append(l)

        for l in unexc:
            D.pop(l)
        unexc = []

        if (D == {}) or (patience == 0):
            break

    for l in D:
        output.append(l)

    output.sort()
    return output


def _over_numberfield(E):
    r"""Return `E`, defined over a NumberField object. This is necessary
    since if `E` is defined over `\QQ`, then we cannot use SAGE commands
    available for number fields.

    INPUT:

    - ``E`` - EllipticCurve - over a number field.

    OUTPUT:

    - If `E` is defined over a NumberField, returns E.

    - If `E` is defined over QQ, returns E defined over the NumberField QQ.

    EXAMPLES::

        sage: E = EllipticCurve([1, 2])
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._over_numberfield(E)
        Elliptic Curve defined by y^2 = x^3 + x + 2 over Number Field in a with defining polynomial x
    """

    K = E.base_field()

    if K == QQ:
        x = QQ['x'].gen()
        K = NumberField(x, 'a')
        E = E.change_ring(K)

    return E


def _tr12(tr, det):
    r"""Compute `X^{12} + Y^{12}` given `X + Y` and `X * Y`.

    INPUT:

    - ``tr`` - The value of `X + Y`.

    - ``det`` - The value of `X * Y`.

    OUTPUT: The value of `X^{12} + Y^{12}`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.gal_reps_number_field import *
        sage: X, Y = QQ['X, Y'].gens()
        sage: sage.schemes.elliptic_curves.gal_reps_number_field._tr12(X + Y, X * Y)
        X^12 + Y^12
    """

    det3 = det**3
    return ((tr * (tr**2 - 3 * det))**2 - 2 * det3)**2 - 2 * det3**2

def deg_one_primes_iter(K, principal_only=False):
    r"""
    Return an iterator over degree 1 primes of ``K``.

    INPUT:

    - ``K`` -- a number field
    - ``principal_only`` -- bool; if ``True``, only yield principal primes

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
         Fractional ideal (-a),
         Fractional ideal (7, a + 3),
         Fractional ideal (7, a + 4)]
        sage: it = deg_one_primes_iter(K, True)
        sage: [next(it) for _ in range(6)]
        [Fractional ideal (-a),
         Fractional ideal (-2*a + 3),
         Fractional ideal (2*a + 3),
         Fractional ideal (a + 6),
         Fractional ideal (a - 6),
         Fractional ideal (-3*a + 4)]
    """
    # imaginary quadratic fields have no principal primes of norm < disc / 4
    start = K.discriminant().abs() // 4 if principal_only and K.signature() == (0,1) else 2

    from sage.arith.misc import primes
    from sage.rings.infinity import infinity
    for p in primes(start=start, stop=infinity):
        for P in K.primes_above(p, degree=1):
            if not principal_only or P.is_principal():
                yield P

def _semistable_reducible_primes(E):
    r"""Find a list containing all semistable primes l unramified in K/QQ
    for which the Galois image for E could be reducible.

    INPUT:

    - ``E`` - EllipticCurve - over a number field.

    OUTPUT:

    A list of primes, which contains all primes `l` unramified in
    `K/\mathbb{QQ}`, such that `E` is semistable at all primes lying
    over `l`, and the Galois image at `l` is reducible. If `E` has CM
    defined over its ground field, a ``ValueError`` is raised.

    EXAMPLES::

        sage: E = EllipticCurve([0, -1, 1, -10, -20]) # X_0(11)
        sage: 5 in sage.schemes.elliptic_curves.gal_reps_number_field._semistable_reducible_primes(E)
        True
    """

    E = _over_numberfield(E)
    K = E.base_field()

    deg_one_primes = deg_one_primes_iter(K, principal_only=True)

    bad_primes = set([]) # This will store the output.

    # We find two primes (of distinct residue characteristics) which are
    # of degree 1, unramified in K/Q, and at which E has good reduction.
    # Both of these primes will give us a nontrivial divisibility constraint
    # on the exceptional primes l. For both of these primes P, we precompute
    # a generator and the trace of Frob_P^12.

    precomp = []
    last_char = 0 # The residue characteristic of the most recent prime.

    while len(precomp) < 2:
        P = next(deg_one_primes)

        # the iterator tests this already
        # if not P.is_principal():
        #     continue

        det = P.norm()
        if det == last_char:
            continue

        if P.ramification_index() != 1:
            continue

        if E.has_bad_reduction(P):
            continue

        tr = E.reduction(P).trace_of_frobenius()
        x = P.gens_reduced()[0]

        precomp.append((x, _tr12(tr, det)))
        last_char = det

    x, tx = precomp[0]
    y, ty = precomp[1]

    Kgal = K.galois_closure('b')
    maps = K.embeddings(Kgal)

    for i in range(2 ** (K.degree() - 1)):
        ## We iterate through all possible characters. ##

        # Here, if i = i_{l-1} i_{l-2} cdots i_1 i_0 in binary, then i
        # corresponds to the character prod sigma_j^{i_j}.

        phi1x = 1
        phi2x = 1
        phi1y = 1
        phi2y = 1

        # We compute the two algebraic characters at x and y:
        for j in range(K.degree()):
            if i % 2 == 1:
                phi1x *= maps[j](x)
                phi1y *= maps[j](y)
            else:
                phi2x *= maps[j](x)
                phi2y *= maps[j](y)
            i = int(i/2)

        # Any prime with reducible image must divide both of:
        gx = phi1x**12 + phi2x**12 - tx
        gy = phi1y**12 + phi2y**12 - ty

        if (gx != 0) or (gy != 0):
            for prime in Integer(Kgal.ideal([gx, gy]).norm()).prime_factors():
                bad_primes.add(prime)

            continue

        ## It is possible that our curve has CM. ##

        # Our character must be of the form Nm^K_F for an imaginary
        # quadratic subfield F of K (which is the CM field if E has CM).
        # We compute F:

        a = (Integer(phi1x + phi2x)**2 - 4 * x.norm()).squarefree_part()

        # See #19229: the name given here, which is not used, should
        # not be the name of the generator of the base field.
        F = QuadraticField(a, 'gal_rep_nf_sqrt_a')

        # Next, we turn K into relative number field over F.

        K = K.relativize(F.embeddings(K)[0], K.variable_name()+'0')
        E = E.change_ring(K.structure()[1])

        ## We try to find a nontrivial divisibility condition. ##

        patience = 5 * K.absolute_degree()
        # Number of Frobenius elements to check before suspecting that E
        # has CM and computing the set of CM j-invariants of K to check.
        # TODO: Is this the best value for this parameter?

        while True:
            P = next(deg_one_primes)

            # the iterator tests this already
            # if not P.is_principal():
            #     continue

            try:
                tr = E.change_ring(P.residue_field()).trace_of_frobenius()
            except ArithmeticError: # Bad reduction at P.
                continue

            x = P.gens_reduced()[0].norm(F)
            div = (x**12).trace() - _tr12(tr, x.norm())

            patience -= 1

            if div != 0:
                # We found our divisibility constraint.

                for prime in Integer(div).prime_factors():
                    bad_primes.add(prime)

                # Turn K back into an absolute number field.

                E = E.change_ring(K.structure()[0])
                K = K.structure()[0].codomain()

                break

            if patience == 0:
                # We suspect that E has CM, so we check:
                f = K.structure()[0]
                if f(E.j_invariant()) in cm_j_invariants(f.codomain()):
                    raise ValueError("The curve E should not have CM.")

    L = sorted(bad_primes)
    return L


def _possible_normalizers(E, SA):
    r"""Find a list containing all primes `l` such that the Galois image at `l`
    is contained in the normalizer of a Cartan subgroup, such that the
    corresponding quadratic character is ramified only at the given primes.

    INPUT:

    - ``E`` - EllipticCurve - over a number field K.

    - ``SA`` - list - a list of primes of K.

    OUTPUT:

    - list - A list of primes, which contains all primes `l` such that the
             Galois image at `l` is contained in the normalizer of a Cartan
             subgroup, such that the corresponding quadratic character is
             ramified only at primes in SA.

    - If `E` has geometric CM that is not defined over its ground field, a
      ValueError is raised.

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

    x = K['x'].gen()
    selmer_group = K.selmer_group(SA, 2) # Generators of the selmer group.

    if selmer_group == []:
        return []

    V = VectorSpace(GF(2), len(selmer_group))
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

        for a in selmer_group:
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

    bad_primes = set([])

    for i in traces_list:
        for p in i.prime_factors():
            bad_primes.add(p)

    # We find the unique vector v in V orthogonal to W:
    v = W.matrix().transpose().kernel().basis()[0]

    # We find the element a of the selmer group corresponding to v:
    a = 1
    for i in range(len(selmer_group)):
        if v[i] == 1:
            a *= selmer_group[i]

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

