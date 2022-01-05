# -*- coding: utf-8 -*-
r"""
Galois representations attached to elliptic curves

Given an elliptic curve `E` over `\QQ`
and a rational prime number `p`, the `p^n`-torsion
`E[p^n]` points of `E` is a representation of the
absolute Galois group `G_{\QQ}` of `\QQ`. As `n` varies
we obtain the Tate module `T_p E` which is a
a representation of `G_\QQ` on a free `\ZZ_p`-module
of rank `2`. As `p` varies the representations
are compatible.

Currently sage can decide whether the Galois module
`E[p]` is reducible, i.e., if `E` admits an isogeny
of degree `p`, and whether the image of
the representation on `E[p]` is surjective onto
`\text{Aut}(E[p]) = GL_2(\GF{p})`.

The following are the most useful functions for the class ``GaloisRepresentation``.

For the reducibility:

- ``is_reducible(p)``
- ``is_irreducible(p)``
- ``reducible_primes()``

For the image:

- ``is_surjective(p)``
- ``non_surjective()``
- ``image_type(p)``

For the classification of the representation

- ``is_semistable(p)``
- ``is_unramified(p, ell)``
- ``is_crystalline(p)``

EXAMPLES::

    sage: E = EllipticCurve('196a1')
    sage: rho = E.galois_representation()
    sage: rho.is_irreducible(7)
    True
    sage: rho.is_reducible(3)
    True
    sage: rho.is_irreducible(2)
    True
    sage: rho.is_surjective(2)
    False
    sage: rho.is_surjective(3)
    False
    sage: rho.is_surjective(5)
    True
    sage: rho.reducible_primes()
    [3]
    sage: rho.non_surjective()
    [2, 3]
    sage: rho.image_type(2)
    'The image is cyclic of order 3.'
    sage: rho.image_type(3)
    'The image is contained in a Borel subgroup as there is a 3-isogeny.'
    sage: rho.image_type(5)
    'The image is all of GL_2(F_5).'

For semi-stable curve it is known that the representation is
surjective if and only if it is irreducible::

    sage: E = EllipticCurve('11a1')
    sage: rho = E.galois_representation()
    sage: rho.non_surjective()
    [5]
    sage: rho.reducible_primes()
    [5]

For cm curves it is not true that there are only finitely many primes for which the
Galois representation mod p is surjective onto `GL_2(\GF{p})`::

    sage: E = EllipticCurve('27a1')
    sage: rho = E.galois_representation()
    sage: rho.non_surjective()
    [0]
    sage: rho.reducible_primes()
    [3]
    sage: E.has_cm()
    True
    sage: rho.image_type(11)
    'The image is contained in the normalizer of a non-split Cartan group. (cm)'

REFERENCES:

- [Ser1972]_
- [Ser1987]_
- [Coj2005]_

AUTHORS:

- chris wuthrich (02/10) - moved from ell_rational_field.py.

"""

######################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
######################################################################

from sage.structure.sage_object import SageObject
import sage.arith.all as arith
from sage.rings.fast_arith import prime_range
import sage.misc.all as misc
from sage.misc.verbose import verbose
import sage.rings.all as rings
from sage.rings.all import RealField, GF

from math import sqrt
from sage.libs.pari.all import pari


def _ex_set(p):
    """
    Gives the set of the only values of trace^2/det
    that appear in a exceptional subgroup in PGL_2(F_p).

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.gal_reps import _ex_set
        sage: for p in prime_range(3,30): print("{} {}".format(p, _ex_set(p)))
        3 [0, 1, 2, 1]
        5 [0, 1, 2, 4]
        7 [0, 1, 2, 4]
        11 [0, 1, 2, 4, 9, 5]
        13 [0, 1, 2, 4]
        17 [0, 1, 2, 4]
        19 [0, 1, 2, 4, 16, 6]
        23 [0, 1, 2, 4]
        29 [0, 1, 2, 4, 25, 7]
    """
    k = GF(p)
    res = [ k(0), k(1), k(2), k(4) ]
    R = k['X']
    f = R([1,-3,1]) #(X**2 - 3*X+1)
    ro = f.roots()
    for a in ro:
        if a[0] not in res:
            res.append(a[0])
    return res


class GaloisRepresentation(SageObject):
    r"""
    The compatible family of Galois representation
    attached to an elliptic curve over the rational numbers.

    Given an elliptic curve `E` over `\QQ` and a rational
    prime number `p`, the `p^n`-torsion `E[p^n]` points of
    `E` is a representation of the absolute Galois group.
    As `n` varies we obtain the Tate module `T_p E` which is
    a representation of the absolute Galois group on a free
    `\ZZ_p`-module of rank `2`. As `p` varies the
    representations are compatible.

    EXAMPLES::

        sage: rho = EllipticCurve('11a1').galois_representation()
        sage: rho
        Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

    """

    def __init__(self, E):
        r"""
        see ``GaloisRepresentation`` for documentation

        EXAMPLES::

            sage: rho = EllipticCurve('11a1').galois_representation()
            sage: loads(rho.dumps()) == rho
            True

        """
        self.__image_type = {}
        self._E = E

    def __repr__(self):
        r"""
        string representation of the class

        EXAMPLES::

            sage: rho = EllipticCurve([0,1]).galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field

        """
        return "Compatible family of Galois representations associated to the " + repr(self._E)

    def __eq__(self,other):
        r"""
        Compares two Galois representations.
        We define tho compatible families of representations
        attached to elliptic curves to be isomorphic if the curves are equal.

        EXAMPLES::

            sage: rho = EllipticCurve('11a1').galois_representation()
            sage: rho2 = EllipticCurve('11a2').galois_representation()
            sage: rho == rho
            True
            sage: rho == rho2
            False
            sage: rho == 34
            False

        """
        # if rho_E = rho_E' then the L-functions agree,
        # so E and E' are isogenous
        # except for p=2
        # the mod p representations will be different
        # for p dividing the degree of the isogeny
        # anyway, there should not be a _compatible_
        # isomorphism between rho and rho' unless E
        # is isomorphic to E'
        # Note that rho can not depend on the Weierstrass model
        if type(self) is not type(other):
            return False
        return self._E.is_isomorphic(other._E)

    def elliptic_curve(self):
        r"""
        The elliptic curve associated to this representation.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: rho = E.galois_representation()
            sage: rho.elliptic_curve() == E
            True

        """
        from copy import copy
        return copy(self._E)

#####################################################################
# reducibility
#####################################################################

    def is_reducible(self, p):
        r"""
        Return True if the mod-p representation is
        reducible. This is equivalent to the existence of an
        isogeny defined over `\QQ` of degree `p` from the
        elliptic curve.

        INPUT:

        -  ``p`` - a prime number

        OUTPUT:

        - a boolean

        The answer is cached.

        EXAMPLES::

            sage: rho = EllipticCurve('121a').galois_representation()
            sage: rho.is_reducible(7)
            False
            sage: rho.is_reducible(11)
            True
            sage: EllipticCurve('11a').galois_representation().is_reducible(5)
            True
            sage: rho = EllipticCurve('11a2').galois_representation()
            sage: rho.is_reducible(5)
            True
            sage: EllipticCurve('11a2').torsion_order()
            1

        """
        try:
            return self.__is_reducible[p]
        except AttributeError:
            self.__is_reducible = {}
        except KeyError:
            pass

        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        # we do is_surjective first, since this is
        # much easier than computing isogeny_class
        t = self._is_surjective(p, A=-1)
        if t:
            self.__is_reducible[p] = False
            return False  # definitely not reducible
        isogeny_matrix = self._E.isogeny_class().matrix(fill=True)
        v = isogeny_matrix.row(0) # first row
        for a in v:
            if a != 0 and a % p == 0:
                self.__is_reducible[p] = True
                return True
        self.__is_reducible[p] = False
        return False

    def is_irreducible(self, p):
        r"""
        Return True if the mod p representation is irreducible.

        INPUT:

        -  ``p`` - a prime number

        OUTPUT:

        - a boolean

        EXAMPLES::

            sage: rho = EllipticCurve('37b').galois_representation()
            sage: rho.is_irreducible(2)
            True
            sage: rho.is_irreducible(3)
            False
            sage: rho.is_reducible(2)
            False
            sage: rho.is_reducible(3)
            True
        """
        return not self.is_reducible(p)

    def reducible_primes(self):
        r"""
        Return a list of the primes `p` such that the mod-`p`
        representation is reducible. For all other primes the
        representation is irreducible.

        EXAMPLES::

            sage: rho = EllipticCurve('225a').galois_representation()
            sage: rho.reducible_primes()
            [3]
        """
        try:
            return self.__reducible_primes
        except AttributeError:
            pass

        E = self._E
        j = E.j_invariant()
        from .isogeny_small_degree import sporadic_j
        if j in sporadic_j: # includes all CM j-invariants
            R = [sporadic_j[j]]
        else:
            R = [l for l in [2,3,5,7,13] if len(E.isogenies_prime_degree(l))>0]
        self.__reducible_primes = R
        return R

#####################################################################
# image
#####################################################################

    def is_surjective(self, p, A=1000):
        r"""
        Return True if the mod-p representation is
        surjective onto `Aut(E[p]) = GL_2(\GF{p})`.

        False if it is not, or None if we were unable to
        determine whether it is or not.

        INPUT:

        -  ``p`` - int (a prime number)

        -  ``A`` - int (a bound on the number of a_p to use)

        OUTPUT:

        - boolean. True if the mod-p representation is surjective
          and False if not.

        The answer is cached.

        EXAMPLES::

            sage: rho = EllipticCurve('37b').galois_representation()
            sage: rho.is_surjective(2)
            True
            sage: rho.is_surjective(3)
            False

        ::

            sage: rho = EllipticCurve('121a1').galois_representation()
            sage: rho.non_surjective()
            [11]
            sage: rho.is_surjective(5)
            True
            sage: rho.is_surjective(11)
            False

            sage: rho = EllipticCurve('121d1').galois_representation()
            sage: rho.is_surjective(5)
            False
            sage: rho.is_surjective(11)
            True

        Here is a case, in which the algorithm does not return an answer::

            sage: rho = EllipticCurve([0,0,1,2580,549326]).galois_representation()
            sage: rho.is_surjective(7)

        In these cases, one can use image_type to get more information about the image::

            sage: rho.image_type(7)
            'The image is contained in the normalizer of a split Cartan group.'


        REMARKS:

        1. If `p \geq 5` then the mod-p representation is
           surjective if and only if the p-adic representation is
           surjective. When `p = 2, 3` there are
           counterexamples. See papers of Dokchitsers and Elkies
           for more details.

        2. For the primes `p=2` and 3, this will always answer either
           True or False. For larger primes it might give None.

        """
        if not arith.is_prime(p):
            raise TypeError("p (=%s) must be prime." % p)
        A = int(A)
        key = (p, A)
        try:
            return self.__is_surjective[key]
        except KeyError:
            pass
        except AttributeError:
            self.__is_surjective = {}

        ans = self._is_surjective(p, A)
        self.__is_surjective[key] = ans
        return ans


    def _is_surjective(self, p, A):
        r"""
        helper function for ``is_surjective``

        The value of `A` is as before, except that
        `A=-1` is a special code to stop before
        testing reducibility (to avoid an infinite loop).

        EXAMPLES::

            sage: rho = EllipticCurve('37b').galois_representation()
            sage: rho._is_surjective(7,100)
            True

        TEST for :trac:`8451`::

            sage: E = EllipticCurve('648a1')
            sage: rho = E.galois_representation()
            sage: rho._is_surjective(5,1000)

        """
        T = self._E.torsion_subgroup().order()
        if T % p == 0 and p != 2 :
            # we could probably determine the group structure directly
            self.__image_type[p] = "The image is meta-cyclic inside a Borel subgroup as there is a %s-torsion point on the curve."%p
            return False

        R = rings.PolynomialRing(self._E.base_ring(), 'x')
        x = R.gen()

        if p == 2:
            # E is isomorphic to  [0,b2,0,8*b4,16*b6]
            b2,b4,b6,b8=self._E.b_invariants()
            f = x**3 + b2*x**2 + 8*b4*x + 16*b6
            if not f.is_irreducible():
                if len(f.roots()) > 2:
                    self.__image_type[p] = "The image is trivial as all 2-torsion points are rational."
                else:
                    self.__image_type[p] = "The image is cyclic of order 2 as there is exactly one rational 2-torsion point."
                return False    #, '2-torsion'
            if arith.is_square(f.discriminant()):
                self.__image_type[p] = "The image is cyclic of order 3."
                return False    #, "A3"
            self.__image_type[p] = "The image is all of GL_2(F_2), i.e. a symmetric group of order 6."
            return True  #, None

        if p == 3:
            # Algorithm: Let f be the 3-division polynomial, which is
            # a polynomial of degree 4.  Then I claim that this
            # polynomial has Galois group S_4 if and only if the
            # representation rhobar_{E,3} is surjective.  If the group
            # is S_4, then S_4 is a quotient of the image of
            # rhobar_{E,3}.  Since S_4 has order 24 and GL_2(F_3)
            # has order 48, the only possibility we have to consider
            # is that the image of rhobar is isomorphic to S_4.
            # But this is not the case because S_4 is not a subgroup
            # of GL_2(F_3).    If it were, it would be normal, since
            # it would have index 2.  But there is a *unique* normal
            # subgroup of GL_2(F_3) of index 2, namely SL_2(F_3),
            # and SL_2(F_3) is not isomorphic to S_4 (S_4 has a normal
            # subgroup of index 2 and SL_2(F_3) does not.)
            # (What's a simple way to see that SL_2(F_3) is the
            # unique index-2 normal subgroup?  I didn't see an obvious
            # reason, so just used the NormalSubgroups command in MAGMA
            # and it output exactly one of index 2.)

            #sage: G = SymmetricGroup(4)
            #sage: [H.group_id() for H in G.conjugacy_classes_subgroups()]
            #[[1, 1], [2, 1], [2, 1], [3, 1], [4, 2], [4, 2], [4, 1], [6, 1], [8, 3], [12, 3], [24, 12]]
            #sage: G = GL(2,GF(3)).as_matrix_group().as_permutation_group()
            #sage: [H.group_id() for H in G.conjugacy_classes_subgroups()]
            #[[1, 1], [2, 1], [2, 1], [3, 1], [4, 2], [4, 1], [6, 2], [6, 1], [6, 1], [8, 4], [8, 1], [8, 3], [12, 4], [16, 8], [24, 3], [48, 29]]

            # Here's Noam Elkies proof for the other direction:

            #> Let E be an elliptic curve over Q.  Is the mod-3
            #> representation E[3]  surjective if and only if the
            #> (degree 4) division polynomial has Galois group S_4?  I
            #> can see why the group being S_4 implies the
            #> representation is surjective, but the converse is not
            #> clear to me.
            # I would have thought that this is the easier part: to
            # say that E[3] is surjective is to say the 3-torsion
            # field Q(E[3]) has Galois group GL_2(Z/3) over Q.  Let
            # E[3]+ be the subfield fixed by the element -1 of
            # GL_2(Z/3).  Then E[3] has Galois group PGL_2(Z/3), which
            # is identified with S_4 by its action on the four
            # 3-element subgroups of E[3].  Each such subgroup is in
            # turn determined by the x-coordinate shared by its two
            # nonzero points.  So, if E[3] is surjective then any
            # permutation of those x-coordinates is realized by some
            # element of Gal(E[3]+/Q).  Thus the Galois group of the
            # division polynomial (whose roots are those
            # x-coordinates) maps surjectively to S_4, which means it
            # equals S_4.

            f = self._E.division_polynomial(3)
            if not f.is_irreducible():
                return False   #, "reducible_3-divpoly"
            n = pari(f).polgalois()[0]
            if n == 24:
                self.__image_type[p] = "The image is all of GL_2(F_3)."
                return True   #, None
            else:
                return False   #, "3-divpoly_galgroup_order_%s"%n

        if self._E.has_cm():
            return False   #, "CM"

        # Now we try to prove that the rep IS surjective.

        Np = self._E.conductor() * p
        signs = []
        # this follows Proposition 19 in Serre.
        # there was a bug in the original implementation,
        # counter-examples were 324b1 and 324d1, 648a1 and 648c1 for p=5
        exclude_exceptional_image = False
        ex_setp = _ex_set(p)
        ell = 4
        k = GF(p)

        if A == -1:
            Am = 1000
        else:
            Am = A
        while ell < Am:
            ell = arith.next_prime(ell)
            if Np % ell != 0:
                a_ell = self._E.ap(ell)
                if a_ell % p != 0:
                    if not exclude_exceptional_image:
                        u = k(a_ell)**2 * k(ell)**(-1)
                        if u not in ex_setp:
                            exclude_exceptional_image = True
                    s = arith.kronecker(a_ell**2 - 4*ell, p)
                    if s != 0 and s not in signs:
                        signs.append(s)
                    if len(signs) == 2 and exclude_exceptional_image:
                        self.__image_type[p] = "The image is all of GL_2(F_%s)."%p
                        return True   #,None

        if A == -1: # we came in from is reducible. Now go out with False
            return False

        if self.is_reducible(p):
            return False  #, Borel

        # if we reach this, then we do not know if it is surjective. Most likely
        # not but we can't be certain. See trac 11271.
        verbose("We can not conclude if the representation is surjective or not. Increasing the parameter A may help.")
        return None

    def non_surjective(self, A=1000):
        r"""
        Return a list of primes p such that the mod-p representation
        *might* not be surjective. If `p` is not in the returned list,
        then the mod-p representation is provably surjective.

        By a theorem of Serre, there are only finitely
        many primes in this list, except when the curve has
        complex multiplication.

        If the curve has CM, we simply return the
        sequence [0] and do no further computation.

        INPUT:

        - ``A`` - an integer (default 1000). By increasing this parameter
          the resulting set might get smaller.

        OUTPUT:

        -  ``list`` - if the curve has CM, returns [0].
           Otherwise, returns a list of primes where mod-p representation is
           very likely not surjective. At any prime not in this list, the
           representation is definitely surjective.

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -38, 90])  # 361A
            sage: E.galois_representation().non_surjective()   # CM curve
            [0]

        ::

            sage: E = EllipticCurve([0, -1, 1, 0, 0]) # X_1(11)
            sage: E.galois_representation().non_surjective()
            [5]

            sage: E = EllipticCurve([0, 0, 1, -1, 0]) # 37A
            sage: E.galois_representation().non_surjective()
            []

            sage: E = EllipticCurve([0,-1,1,-2,-1])   # 141C
            sage: E.galois_representation().non_surjective()
            [13]

        ::

            sage: E = EllipticCurve([1,-1,1,-9965,385220]) # 9999a1
            sage: rho = E.galois_representation()
            sage: rho.non_surjective()
            [2]

            sage: E = EllipticCurve('324b1')
            sage: rho = E.galois_representation()
            sage: rho.non_surjective()
            [3, 5]

        ALGORITHM:
        We first find an upper bound `B` on the possible primes. If `E`
        is semi-stable, we can take `B=11` by a result of Mazur. There is
        a bound by Serre in the case that the `j`-invariant is not integral
        in terms of the smallest prime of good reduction. Finally
        there is an unconditional bound by Cojocaru, but which depends
        on the conductor of `E`.
        For the prime below that bound we call ``is_surjective``.

        """
        if self._E.has_cm():
            verbose("cm curve")
            return [0]
        N = self._E.conductor()
        if self._E.is_semistable():
            # Mazur's bound
            C = 11
            verbose("semistable -- so bound is 11")
        elif not self._E.j_invariant().is_integral():
            # prop 24 in Serre
            vs = self._E.j_invariant().denominator().prime_factors()
            C1 = arith.gcd([-arith.valuation(self._E.j_invariant(),v) for v in vs])
            p0 = 2
            while self._E.has_bad_reduction(p0):
                p0 = arith.next_prime(p0+1)
            C2 = (sqrt(p0)+1)**8
            C = max(C1,C2)
            verbose("j is not integral -- Serre's bound is %s"%C)
            C3 = 1 + 4*sqrt(6)*int(N)/3 * sqrt(misc.mul([1+1.0/int(p) for p,_ in arith.factor(N)]))
            C = min(C,C3)
            verbose("conductor = %s, and bound is %s"%(N,C))
        else:
            # Cojocaru's bound (depends on the conductor)
            C = 1 + 4*sqrt(6)*int(N)/3 * sqrt(misc.mul([1+1.0/int(p) for p,_ in arith.factor(N)]))
            verbose("conductor = %s, and bound is %s"%(N,C))
        B = []
        p = 2
        while p <= C:
            t = self.is_surjective(p, A=A)
            verbose("(%s,%s)"%(p,t))
            # both False and None will be appended here.
            if not t:
                B.append(p)
            p = arith.next_prime(p)
        return B

    def image_type(self, p):
        r"""
        Return a string describing the image of the
        mod-p representation.
        The result is provably correct, but only
        indicates what sort of an image we have. If
        one wishes to determine the exact group one
        needs to work a bit harder. The probabilistic
        method of image_classes or Sutherland's galrep
        package can give a very good guess what the
        image should be.

        INPUT:

        - ``p``  a prime number

        OUTPUT:

        - a string.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: rho = E.galois_representation()
            sage: rho.image_type(5)
            'The image is all of GL_2(F_5).'

            sage: E = EllipticCurve('11a1')
            sage: rho = E.galois_representation()
            sage: rho.image_type(5)
            'The image is meta-cyclic inside a Borel subgroup as there is a 5-torsion point on the curve.'

            sage: EllipticCurve('27a1').galois_representation().image_type(5)
            'The image is contained in the normalizer of a non-split Cartan group. (cm)'
            sage: EllipticCurve('30a1').galois_representation().image_type(5)
            'The image is all of GL_2(F_5).'
            sage: EllipticCurve("324b1").galois_representation().image_type(5)
            'The image in PGL_2(F_5) is the exceptional group S_4.'

            sage: E = EllipticCurve([0,0,0,-56,4848])
            sage: rho = E.galois_representation()

            sage: rho.image_type(5)
            'The image is contained in the normalizer of a split Cartan group.'

            sage: EllipticCurve('49a1').galois_representation().image_type(7)
            'The image is contained in a Borel subgroup as there is a 7-isogeny.'

            sage: EllipticCurve('121c1').galois_representation().image_type(11)
            'The image is contained in a Borel subgroup as there is a 11-isogeny.'
            sage: EllipticCurve('121d1').galois_representation().image_type(11)
            'The image is all of GL_2(F_11).'
            sage: EllipticCurve('441f1').galois_representation().image_type(13)
            'The image is contained in a Borel subgroup as there is a 13-isogeny.'

            sage: EllipticCurve([1,-1,1,-5,2]).galois_representation().image_type(5)
            'The image is contained in the normalizer of a non-split Cartan group.'
            sage: EllipticCurve([0,0,1,-25650,1570826]).galois_representation().image_type(5)
            'The image is contained in the normalizer of a split Cartan group.'
            sage: EllipticCurve([1,-1,1,-2680,-50053]).galois_representation().image_type(7)    # the dots (...) in the output fix #11937 (installed 'Kash' may give additional output); long time (2s on sage.math, 2014)
            'The image is a... group of order 18.'
            sage: EllipticCurve([1,-1,0,-107,-379]).galois_representation().image_type(7)       # the dots (...) in the output fix #11937 (installed 'Kash' may give additional output); long time (1s on sage.math, 2014)
            'The image is a... group of order 36.'
            sage: EllipticCurve([0,0,1,2580,549326]).galois_representation().image_type(7)
            'The image is contained in the normalizer of a split Cartan group.'

        Test :trac:`14577`::

            sage: EllipticCurve([0, 1, 0, -4788, 109188]).galois_representation().image_type(13)
            'The image in PGL_2(F_13) is the exceptional group S_4.'

        Test :trac:`14752`::

            sage: EllipticCurve([0, 0, 0, -1129345880,-86028258620304]).galois_representation().image_type(11)
            'The image is contained in the normalizer of a non-split Cartan group.'

        For `p=2`::

            sage: E = EllipticCurve('11a1')
            sage: rho = E.galois_representation()
            sage: rho.image_type(2)
            'The image is all of GL_2(F_2), i.e. a symmetric group of order 6.'

            sage: rho = EllipticCurve('14a1').galois_representation()
            sage: rho.image_type(2)
            'The image is cyclic of order 2 as there is exactly one rational 2-torsion point.'

            sage: rho = EllipticCurve('15a1').galois_representation()
            sage: rho.image_type(2)
            'The image is trivial as all 2-torsion points are rational.'

            sage: rho = EllipticCurve('196a1').galois_representation()
            sage: rho.image_type(2)
            'The image is cyclic of order 3.'

        `p=3`::

            sage: rho = EllipticCurve('33a1').galois_representation()
            sage: rho.image_type(3)
            'The image is all of GL_2(F_3).'

            sage: rho = EllipticCurve('30a1').galois_representation()
            sage: rho.image_type(3)
            'The image is meta-cyclic inside a Borel subgroup as there is a 3-torsion point on the curve.'

            sage: rho = EllipticCurve('50b1').galois_representation()
            sage: rho.image_type(3)
            'The image is contained in a Borel subgroup as there is a 3-isogeny.'

            sage: rho = EllipticCurve('3840h1').galois_representation()
            sage: rho.image_type(3)
            'The image is contained in a dihedral group of order 8.'

            sage: rho = EllipticCurve('32a1').galois_representation()
            sage: rho.image_type(3)
            'The image is a semi-dihedral group of order 16, gap.SmallGroup([16,8]).'


        ALGORITHM: Mainly based on Serre's paper.
        """
        if not arith.is_prime(p):
            raise TypeError("p (=%s) must be prime." % p)
        try:
            return self.__image_type[p]
        except KeyError:
            pass
        except AttributeError:
            self.__image_type = {}

        # we eliminate step by step possibilities.
        # checks if the image is all of GL_2, while doing so, it may detect certain other classes already
        # see _is_surjective. It will have set __image_type

        self._is_surjective(p, A=1000)  # this sets __image_type
        try:
            return self.__image_type[p]
        except KeyError:
            pass

        # check if the rep is reducible

        if self.is_reducible(p):
            self.__image_type[p] = "The image is contained in a Borel subgroup as there is a %s-isogeny."%p
            return self.__image_type[p]

        # if we are then the image of rho is not surjective and not contained in a Borel subgroup
        # there are three cases left:
        # it could be in a normalizer of a split Cartan,
        #                  normalizer of a non-split Cartan,
        # or the image in PGL_2 is one of the three exceptional groups A_4 S_4 A_5

        non_split_str = "The image is contained in the normalizer of a non-split Cartan group."
        split_str =     "The image is contained in the normalizer of a split Cartan group."
        s4_str =        "The image in PGL_2(F_%s) is the exceptional group S_4."%p
        a4_str =        "The image in PGL_2(F_%s) is the exceptional group A_4."%p
        a5_str =        "The image in PGL_2(F_%s) is the exceptional group A_5."%p

        # we first treat p=3 and 5 separately. p=2 has already been done.

        if p == 3:
            # this implies that the image of rhobar in PGL_2 = S_4
            # determines completely the image of rho
            f = self._E.division_polynomial(3)
            if not f.is_irreducible():
                # must be a product of two polynomials of degree 2
                self.__image_type[p] = "The image is contained in a dihedral group of order 8."
                return self.__image_type[p]
            n = pari(f).polgalois()[0]
            # the following is due to a simple classification of all subgroups of GL_2(F_3)
            if n == 2:
                self.__image_type[p] = "The image is a cyclic group of order 4."
            elif n == 4:
                for ell in prime_range(5,1000):
                    if ell % 3 == 2 and self._E.ap(ell) % 3 != 0:
                        # there is an element of order 8 in the image
                        self.__image_type[p] = "The image is a cyclic group of order 8."
                        return self.__image_type[p]
                self.__image_type[p] = "The image is a group of order 8, most likely a quaternion group."
            elif n == 6:
                self.__image_type[p] = "The image is a dihedral group of order 12."
            elif n == 8:
                self.__image_type[p] = "The image is a semi-dihedral group of order 16, gap.SmallGroup([16,8])."
            elif n == 12:
                self.__image_type[p] = "The image is SL_2(F_3)."
            else:
                raise RuntimeError("Bug in image_type for p = 3.")
            return self.__image_type[p]

        # we also eliminate cm curves

        if self._E.has_cm():
            if self._E.is_good(p) and self._E.is_ordinary(p):
                self.__image_type[p] = split_str + " (cm)"
                return self.__image_type[p]
            if self._E.is_supersingular(p):
                self.__image_type[p] = non_split_str + " (cm)"
                return self.__image_type[p]
            else:
                # if the reduction is bad (additive nec.) then
                # the image is in a Borel subgroup and we should have found this before
                raise NotImplementedError("image_type is not implemented for cm-curves at bad prime.")

        # now to p=5 where a lot of non-standard thing happen
        # we run through primes if we hit an element of order 6 in PGL_2, we know that it is the normaliser of a NON-split Cartan
        # if we find both an element of order 3 and one of order 4, we know that we have a S_4 in PGL_2

        if p == 5:
            # we filter here a few cases and leave the rest to the computation of the Galois group later
            ell = 1
            k = GF(p)
            Np = self._E.conductor()*p
            has_an_el_order_4 = False
            has_an_el_order_3 = False
            while ell < 10000:
                ell = arith.next_prime(ell)
                if Np % ell != 0:
                    a_ell = self._E.ap(ell)
                    u = k(a_ell)**2 * k(ell)**(-1)
                    if u == 3:
                        verbose("found an element of order 6",2)
                        # found an element of order 6:
                        self.__image_type[p] = non_split_str
                        return self.__image_type[p]

                    if u == 2 and not has_an_el_order_4:
                        # found an element of order 4
                        verbose("found an element of order 4",2)
                        has_an_el_order_4 = True
                        if has_an_el_order_3:
                            self.__image_type[p] = s4_str
                            return self.__image_type[p]

                    if u == 1 and not has_an_el_order_3:
                        # found an element of order 3
                        verbose("found an element of order 3",2)
                        has_an_el_order_3 = True
                        if has_an_el_order_4:
                            self.__image_type[p] = s4_str
                            return self.__image_type[p]

            verbose("p=5 and we could not determine the image, yet", 2)
            # we have not yet determined the image, there are only the following possible subgroups of PGL_2
            # (unless we were unlucky and none of the elements of order 6 showed up above, for instance)
            # A_4       of order 12 with elements of order 2 and 3
            # D_4       of order  8 with elements of order 2 and 4  (normalizer of the SPLIT)
            # Z/2 x Z/2  of order 4 with elements of order 2        (in the normalizer of the SPLIT)
            # S_3        of order 6 with elements of order 2 and 3  (inside the normalizer of the NON-split)

            # we compute the splitting field of the 5-division polynomial. Its degree is equal to the above order or the double of it.
            # That allows us to determine almost all cases.

            f = self._E.division_polynomial(5)
            from sage.rings.number_field.splitting_field import SplittingFieldAbort
            try:
                K = f.splitting_field('x', degree_multiple=240, abort_degree=24)
            except SplittingFieldAbort:
                pass
            else:
                if K.degree() in (4,8,16):
                    self.__image_type[p] = split_str
                    return self.__image_type[p]
                if K.degree() == 24:
                    self.__image_type[p] = a4_str
                    return self.__image_type[p]
                if K.degree() == 6:
                    self.__image_type[p] = non_split_str
                    return self.__image_type[p]
                if K.degree() == 12:
                    # PGL - image could be a S_3 in the normalizer of the split or A4
                    self.__image_type[p] = "The image is of order 6 or 12. Probably contained in the normalizer of the split Cartan g."
                    return self.__image_type[p]

        ## now E has no cm, is not semi-stable,
        ## p > 5,
        ## rho_p it is (most probably) not surjective , it is not contained in a Borel subgroup.
        # trying to detect that the image is not exceptional in PGL_2
        # this uses Serre 2.6.iii and Prop 19
        # the existence of certain classes in the image rules out certain cases.
        # we run over the small prime until we are left with only one case
        # for p = 5 this could never distinguish any from an exceptional S_4 or A_4,
        # that is why the case 5 is treated a part before

        else:
            ex_setp = _ex_set(p)
            ell = 1
            k = GF(p)
            Np = self._E.conductor()*p
            could_be_exc = 1
            could_be_split = 1
            could_be_non_split = 1
            # loops over primes as long as we still have two options left
            while ell < 10000 and (could_be_exc + could_be_split + could_be_non_split  > 1):
                ell = arith.next_prime(ell)
                if Np % ell != 0:
                    a_ell = self._E.ap(ell)
                    u = k(a_ell)**2 * k(ell)**(-1)
                    if (u not in ex_setp) and could_be_exc == 1:
                        # it can not be in the exceptional
                        verbose("the image cannot be exceptional, found u=%s"%u,2)
                        could_be_exc = 0
                    if a_ell != 0 and arith.kronecker(a_ell**2 - 4*ell,p) == 1 and could_be_non_split == 1:
                        # it can not be in the normalizer of the non-split Cartan
                        verbose("the image cannot be non-split, found u=%s"%u,2)
                        could_be_non_split = 0
                    if a_ell != 0 and arith.kronecker(a_ell**2 - 4*ell,p) == -1 and could_be_split == 1:
                        # it can not be in the normalizer of the split Cartan
                        verbose("the image cannot be split, found u=%s"%u,2)
                        could_be_split = 0

            assert could_be_exc + could_be_split + could_be_non_split  > 0, "bug in image_type."

            if could_be_exc + could_be_split + could_be_non_split == 1:
                # it is only one of the three cases:
                if could_be_split == 1 :
                    self.__image_type[p] = split_str
                    return self.__image_type[p]
                if could_be_non_split == 1 :
                    self.__image_type[p] = non_split_str
                    return self.__image_type[p]
                if could_be_exc == 1:
                    # here we can distinguish further
                    could_be_a4 = 1
                    could_be_s4 = 1
                    could_be_a5 = 1
                    if p % 5 != 1 and p % 5 != 4 :
                        could_be_a5 = 0
                    # elements of order 5 # bug corrected see trac 14577
                    R = k['X']
                    f = R([1,-3,1]) #(X**2 - 3*X+1)
                    el5 = f.roots()
                    # loops over primes as long as we still have two options left
                    while ell < 10000 and (could_be_s4 + could_be_a4 + could_be_a5  > 1):
                        ell = arith.next_prime(ell)
                        if Np % ell != 0:
                            a_ell = self._E.ap(ell)
                            u = k(a_ell)**2 * k(ell)**(-1)
                            if u == 2:
                                # it can not be A4 not A5 as they have no elements of order 4
                                could_be_a4 = 0
                                could_be_a5 = 0
                            if u in el5 :
                                # it can not be A4 or S4 as they have no elements of order 5
                                could_be_a4 = 0
                                could_be_s4 = 0

                    assert (could_be_s4 + could_be_a4 + could_be_a5  > 0), "bug in image_type."

                    if could_be_s4 + could_be_a4 + could_be_a5 == 1:
                        if could_be_s4 == 1:
                            self.__image_type[p] = s4_str
                            return self.__image_type[p]
                        if could_be_a4 == 1:
                            self.__image_type[p] = a4_str
                            return self.__image_type[p]
                        if could_be_a5 == 1:
                            self.__image_type[p] = a5_str
                            return self.__image_type[p]

                    else:
                        self.__image_type[p] = "The image in PGL_2(F_%s) is an exceptional group A_4, S_4 or A_5, but we could not determine which one."%p
                        return self.__image_type[p]

        # If all fails, we probably have a fairly small group and we can try to detect it using the galois_group
        if p <= 13:
            K = self._E.division_field(p, 'z')
            d = K.absolute_degree()

            verbose("field of degree %s.  try to compute Galois group"%(d),2)
            # If the degree is too big, we have no chance at the Galois
            # group.  K.galois_group calls is_galois which used to rely on
            # pari's Galois group computations, so degree < 12
            try:
                if d > 15:
                    raise Exception()
                G = K.galois_group()
            except Exception:
                self.__image_type[p] = "The image is a group of order %s."%d
                return self.__image_type[p]

            else:
                if G.is_abelian():
                    ab = ""
                else:
                    ab = "non-"
                self.__image_type[p] = "The image is a " + ab + "abelian group of order %s."%G.order()
                return self.__image_type[p]

        ## everything failed :

        self.__image_type[p] = "The image could not be determined. Sorry."
        return self.__image_type[p]


    def image_classes(self,p,bound=10000):
        r"""
        This function returns, given the representation `\rho`
        a list of `p` values that add up to 1, representing the
        frequency of the conjugacy classes of the projective image
        of `\rho` in `PGL_2(\GF{p})`.

        Let `M` be a matrix in `GL_2(\GF{p})`, then define
        `u(M) = \text{tr}(M)^2/\det(M)`, which only depends on the
        conjugacy class of `M` in `PGL_2(\GF{p})`. Hence this defines
        a map `u: PGL_2(\GF{p}) \to \GF{p}`, which is almost
        a bijection between conjugacy classes of the source
        and `\GF{p}` (the elements of order `p` and the identity
        map to `4` and both classes of elements of order 2 map to 0).

        This function returns the frequency with which the values of
        `u` appeared among the images of the Frobenius elements
        `a_{\ell}`at `\ell` for good primes `\ell\neq p` below a
        given ``bound``.

        INPUT:

        -  a prime ``p``

        -  a natural number ``bound`` (optional, default=10000)

        OUTPUT:

        - a list of `p` real numbers in the interval `[0,1]` adding up to 1

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: rho = E.galois_representation()
            sage: rho.image_classes(5)
            [0.2095, 0.1516, 0.2445, 0.1728, 0.2217]

            sage: E = EllipticCurve('11a1')
            sage: rho = E.galois_representation()
            sage: rho.image_classes(5)
            [0.2467, 0.0000, 0.5049, 0.0000, 0.2484]

        ::

            sage: EllipticCurve('27a1').galois_representation().image_classes(5)
            [0.5839, 0.1645, 0.0000, 0.1702, 0.08143]
            sage: EllipticCurve('30a1').galois_representation().image_classes(5)
            [0.1956, 0.1801, 0.2543, 0.1728, 0.1972]
            sage: EllipticCurve('32a1').galois_representation().image_classes(5)
            [0.6319, 0.0000, 0.2492, 0.0000, 0.1189]
            sage: EllipticCurve('900a1').galois_representation().image_classes(5)
            [0.5852, 0.1679, 0.0000, 0.1687, 0.07824]
            sage: EllipticCurve('441a1').galois_representation().image_classes(5)
            [0.5860, 0.1646, 0.0000, 0.1679, 0.08150]
            sage: EllipticCurve('648a1').galois_representation().image_classes(5)
            [0.3945, 0.3293, 0.2388, 0.0000, 0.03749]

        ::

            sage: EllipticCurve('784h1').galois_representation().image_classes(7)
            [0.5049, 0.0000, 0.0000, 0.0000, 0.4951, 0.0000, 0.0000]
            sage: EllipticCurve('49a1').galois_representation().image_classes(7)
            [0.5045, 0.0000, 0.0000, 0.0000, 0.4955, 0.0000, 0.0000]

            sage: EllipticCurve('121c1').galois_representation().image_classes(11)
            [0.1001, 0.0000, 0.0000, 0.0000, 0.1017, 0.1953, 0.1993, 0.0000, 0.0000, 0.2010, 0.2026]
            sage: EllipticCurve('121d1').galois_representation().image_classes(11)
            [0.08869, 0.07974, 0.08706, 0.08137, 0.1001, 0.09439, 0.09764, 0.08218, 0.08625, 0.1017, 0.1009]

            sage: EllipticCurve('441f1').galois_representation().image_classes(13)
            [0.08232, 0.1663, 0.1663, 0.1663, 0.08232, 0.0000, 0.1549, 0.0000, 0.0000, 0.0000, 0.0000, 0.1817, 0.0000]

        REMARKS:

        Conjugacy classes of subgroups of `PGL_2(\GF{5})`

        For the case `p=5`, the order of an element determines almost the value of `u`:

        +-------+---+---+---+---+--------+
        |`u`    | 0 | 1 | 2 | 3 | 4      |
        +-------+---+---+---+---+--------+
        |orders | 2 | 3 | 4 | 6 | 1 or 5 |
        +-------+---+---+---+---+--------+

        Here we give here the full table of all conjugacy classes of subgroups with the values
        that ``image_classes`` should give (as ``bound`` tends to `\infty`). Comparing with the output
        of the above examples, it is now easy to guess what the image is.

        +---------+-----+------------------------------------------+
        |subgroup |order| frequencies of values of `u`             |
        +=========+=====+==========================================+
        | trivial | 1   | [0.0000, 0.0000, 0.0000, 0.0000, 1.000]  |
        +---------+-----+------------------------------------------+
        | cyclic  | 2   | [0.5000, 0.0000, 0.0000, 0.0000, 0.5000] |
        +---------+-----+------------------------------------------+
        | cyclic  |  2  |  [0.5000, 0.0000, 0.0000, 0.0000, 0.5000]|
        +---------+-----+------------------------------------------+
        | cyclic  |  3  |  [0.0000, 0.6667, 0.0000, 0.0000, 0.3333]|
        +---------+-----+------------------------------------------+
        | Klein   |  4  |  [0.7500, 0.0000, 0.0000, 0.0000, 0.2500]|
        +---------+-----+------------------------------------------+
        | cyclic  |  4  |  [0.2500, 0.0000, 0.5000, 0.0000, 0.2500]|
        +---------+-----+------------------------------------------+
        | Klein   |  4  |  [0.7500, 0.0000, 0.0000, 0.0000, 0.2500]|
        +---------+-----+------------------------------------------+
        | cyclic  |  5  |   [0.0000, 0.0000, 0.0000, 0.0000, 1.000]|
        +---------+-----+------------------------------------------+
        | cyclic  |  6  |  [0.1667, 0.3333, 0.0000, 0.3333, 0.1667]|
        +---------+-----+------------------------------------------+
        | `S_3`   |  6  |  [0.5000, 0.3333, 0.0000, 0.0000, 0.1667]|
        +---------+-----+------------------------------------------+
        | `S_3`   |  6  | [0.5000, 0.3333, 0.0000, 0.0000, 0.1667] |
        +---------+-----+------------------------------------------+
        | `D_4`   |  8  |  [0.6250, 0.0000, 0.2500, 0.0000, 0.1250]|
        +---------+-----+------------------------------------------+
        | `D_5`   |  10 |  [0.5000, 0.0000, 0.0000, 0.0000, 0.5000]|
        +---------+-----+------------------------------------------+
        | `A_4`   |  12 | [0.2500, 0.6667, 0.0000, 0.0000, 0.08333]|
        +---------+-----+------------------------------------------+
        | `D_6`   |  12 | [0.5833, 0.1667, 0.0000, 0.1667, 0.08333]|
        +---------+-----+------------------------------------------+
        | Borel   |  20 |  [0.2500, 0.0000, 0.5000, 0.0000, 0.2500]|
        +---------+-----+------------------------------------------+
        | `S_4`   |  24 | [0.3750, 0.3333, 0.2500, 0.0000, 0.04167]|
        +---------+-----+------------------------------------------+
        | `PSL_2` |  60 |  [0.2500, 0.3333, 0.0000, 0.0000, 0.4167]|
        +---------+-----+------------------------------------------+
        | `PGL_2` |  120| [0.2083, 0.1667, 0.2500, 0.1667, 0.2083] |
        +---------+-----+------------------------------------------+
        """
        res = [0 for i in range(p)]
        co = 0
        ell = 1
        while ell <= bound:
            ell = arith.next_prime(ell)
            if ell != p and self._E.is_good(ell):
                d = (self._E.ap(ell)**2 * ell.inverse_mod(p)) % p
                res[d] += 1
                co += 1
        Rt = RealField(16)
        res = [Rt(x)/Rt(co) for x in res]
        return res

#####################################################################
# classification of ell and p-adic reps
#####################################################################

# ell-adic reps

    def is_unramified(self,p,ell):
        r"""
        Return true if the Galois representation to `GL_2(\ZZ_p)` is unramified at `\ell`, i.e.
        if the inertia group at a place above `\ell` in `\text{Gal}(\bar\QQ/\QQ)` has trivial image in
        `GL_2(\ZZ_p)`.

        For a Galois representation attached to an elliptic curve `E`, this returns True if `\ell\neq p`
        and `E` has good reduction at `\ell`.

        INPUT:

        - ``p``   a prime
        - ``ell`` another prime

        OUTPUT:

        - Boolean

        EXAMPLES::

            sage: rho = EllipticCurve('20a3').galois_representation()
            sage: rho.is_unramified(5,7)
            True
            sage: rho.is_unramified(5,5)
            False
            sage: rho.is_unramified(7,5)
            False

        This says that the 5-adic representation is unramified at 7, but the 7-adic representation is ramified at 5.
        """
        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        if not arith.is_prime(ell):
            raise ValueError('ell (=%s) must be prime' % ell)
        return (ell != p) and self._E.has_good_reduction(ell)

    def is_unipotent(self,p,ell):
        r"""
        Return true if the Galois representation to `GL_2(\ZZ_p)` is unipotent at `\ell\neq p`, i.e.
        if the inertia group at a place above `\ell` in `\text{Gal}(\bar\QQ/\QQ)` maps into a Borel subgroup.

        For a Galois representation attached to an elliptic curve `E`, this returns True if
        `E` has semi-stable reduction at `\ell`.

        INPUT:

        - ``p``   a prime
        - ``ell`` a different prime

        OUTPUT:

        - Boolean

        EXAMPLES::

            sage: rho = EllipticCurve('120a1').galois_representation()
            sage: rho.is_unipotent(2,5)
            True
            sage: rho.is_unipotent(5,2)
            False
            sage: rho.is_unipotent(5,7)
            True
            sage: rho.is_unipotent(5,3)
            True
            sage: rho.is_unipotent(5,5)
            Traceback (most recent call last):
            ...
            ValueError: unipotent is not defined for l = p, use semistable instead.
        """
        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        if not arith.is_prime(ell):
            raise ValueError('ell (=%s) must be prime' % ell)
        if ell == p:
            raise ValueError("unipotent is not defined for l = p, use semistable instead.")
        return not self._E.has_additive_reduction(ell)

    def is_quasi_unipotent(self,p,ell):
        r"""
        Return true if the Galois representation to `GL_2(\ZZ_p)` is quasi-unipotent at `\ell\neq p`, i.e. if there is a finite extension `K/\QQ` such that the inertia group at a place above `\ell` in `\text{Gal}(\bar\QQ/K)` maps into a Borel subgroup.

        For a Galois representation attached to an elliptic curve `E`, this returns always True.

        INPUT:

        - ``p``   a prime
        - ``ell`` a different prime

        OUTPUT:

        - Boolean

        EXAMPLES::

            sage: rho = EllipticCurve('11a3').galois_representation()
            sage: rho.is_quasi_unipotent(11,13)
            True
        """
        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        if not arith.is_prime(ell):
            raise ValueError('ell (=%s) must be prime' % ell)
        if ell == p:
            raise ValueError("quasi unipotent is not defined for l = p, use semistable instead.")
        return True

# p-adic reps

    def is_ordinary(self,p):
        r"""
        Return true if the `p`-adic Galois representation to `GL_2(\ZZ_p)` is ordinary, i.e.
        if the image of the decomposition group in `\text{Gal}(\bar\QQ/\QQ)` above he prime `p` maps into
        a Borel subgroup.

        For an elliptic curve `E`, this is to ask whether `E` is ordinary at `p`, i.e. good ordinary or multiplicative.

        INPUT:

        - ``p`` a prime

        OUTPUT:

        - a Boolean

        EXAMPLES::

            sage: rho = EllipticCurve('11a3').galois_representation()
            sage: rho.is_ordinary(11)
            True
            sage: rho.is_ordinary(5)
            True
            sage: rho.is_ordinary(19)
            False
        """
        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        if self._E.has_additive_reduction(p):
            raise NotImplementedError('is_ordinary is only implemented for semi-stable representations')
        return self._E.has_multiplicative_reduction(p) or (self._E.has_good_reduction(p) and self._E.ap(p) % p != 0)

    def is_crystalline(self,p):
        r"""
        Return true is the `p`-adic Galois representation to `GL_2(\ZZ_p)` is crystalline.

        For an elliptic curve `E`, this is to ask whether `E` has good reduction at `p`.

        INPUT:

        - ``p`` a prime

        OUTPUT:

        - a Boolean

        EXAMPLES::

            sage: rho = EllipticCurve('64a1').galois_representation()
            sage: rho.is_crystalline(5)
            True
            sage: rho.is_crystalline(2)
            False
        """
        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        return self._E.has_good_reduction(p)

    def is_potentially_crystalline(self,p):
        r"""
        Return true is the `p`-adic Galois representation to `GL_2(\ZZ_p)` is potentially crystalline, i.e.
        if there is a finite extension `K/\QQ_p` such that the `p`-adic representation becomes crystalline.

        For an elliptic curve `E`, this is to ask whether `E` has potentially good reduction at `p`.

        INPUT:

        - ``p`` a prime

        OUTPUT:

        - a Boolean

        EXAMPLES::

            sage: rho = EllipticCurve('37b1').galois_representation()
            sage: rho.is_potentially_crystalline(37)
            False
            sage: rho.is_potentially_crystalline(7)
            True
        """
        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        return self._E.j_invariant().valuation(p) >= 0


    def is_semistable(self,p):
        r"""
        Return true if the `p`-adic Galois representation to `GL_2(\ZZ_p)` is semistable.

        For an elliptic curve `E`, this is to ask whether `E` has semistable reduction at `p`.

        INPUT:

        - ``p`` a prime

        OUTPUT:

        - a Boolean

        EXAMPLES::

            sage: rho = EllipticCurve('20a3').galois_representation()
            sage: rho.is_semistable(2)
            False
            sage: rho.is_semistable(3)
            True
            sage: rho.is_semistable(5)
            True
        """
        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        return not self._E.has_additive_reduction(p)

    def is_potentially_semistable(self,p):
        r"""
        Return true if the `p`-adic Galois representation to `GL_2(\ZZ_p)` is potentially semistable.

        For an elliptic curve `E`, this returns True always

        INPUT:

        - ``p`` a prime

        OUTPUT:

        - a Boolean

        EXAMPLES::

            sage: rho = EllipticCurve('27a2').galois_representation()
            sage: rho.is_potentially_semistable(3)
            True
        """
        if not arith.is_prime(p):
            raise ValueError('p (=%s) must be prime' % p)
        return True
