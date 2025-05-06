r"""
Conductor and reduction types for genus 2 curves

AUTHORS:

- Qing Liu and Henri Cohen (1994-1998): wrote genus2reduction C
  program

- William Stein (2006-03-05): wrote Sage interface to genus2reduction

- Jeroen Demeyer (2014-09-17): replace genus2reduction program by PARI
  library call (:issue:`15808`)

ACKNOWLEDGMENT: (From Liu's website:) Many thanks to Henri Cohen who started
writing this program. After this program is available, many people pointed out
to me (mathematical as well as programming) bugs : B. Poonen, E. Schaefer, C.
Stahlke, M. Stoll, F. Villegas.  So thanks to all of them. Thanks also go to
Ph. Depouilly who help me to compile the program.

Also Liu has given me explicit permission to include genus2reduction with Sage
and for people to modify the C source code however they want.
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2014 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.libs.pari import pari

roman_numeral = ["", "I", "II", "III", "IV", "V", "VI", "VII"]


class ReductionData(SageObject):
    r"""
    Reduction data for a genus 2 curve.

    How to read ``local_data`` attribute, i.e., if this
    class is R, then the following is the meaning of
    ``R.local_data[p]``.

    For each prime number `p` dividing the discriminant of
    `y^2+Q(x)y=P(x)`, there are two lines.

    The first line contains information about the stable reduction
    after field extension. Here are the meanings of the symbols of
    stable reduction:

    (I) The stable reduction is smooth (i.e. the curve has potentially
    good reduction).

    (II) The stable reduction is an elliptic curve `E` with an
    ordinary double point. `j` mod `p` is the modular
    invariant of `E`.

    (III) The stable reduction is a projective line with two ordinary
    double points.

    (IV) The stable reduction is two projective lines crossing
    transversally at three points.

    (V) The stable reduction is the union of two elliptic curves
    `E_1` and `E_2` intersecting transversally at one
    point. Let `j_1`, `j_2` be their modular
    invariants, then `j_1+j_2` and `j_1 j_2` are
    computed (they are numbers mod `p`).

    (VI) The stable reduction is the union of an elliptic curve
    `E` and a projective line which has an ordinary double
    point. These two components intersect transversally at one point.
    `j` mod `p` is the modular invariant of
    `E`.

    (VII) The stable reduction is as above, but the two components are
    both singular.

    In the cases (I) and (V), the Jacobian `J(C)` has
    potentially good reduction. In the cases (III), (IV) and (VII),
    `J(C)` has potentially multiplicative reduction. In the two
    remaining cases, the (potential) semi-abelian reduction of
    `J(C)` is extension of an elliptic curve (with modular
    invariant `j` mod `p`) by a torus.

    The second line contains three data concerning the reduction at
    `p` without any field extension.


    #. The first symbol describes the REDUCTION AT `p` of
       `C`. We use the symbols of Namikawa-Ueno for the type of
       the reduction (Namikawa, Ueno:"The complete classification of
       fibers in pencils of curves of genus two", Manuscripta Math., vol.
       9, (1973), pages 143-186.) The reduction symbol is followed by the
       corresponding page number (or just an indication) in the above
       article. The lower index is printed by , for instance, [I2-II-5]
       means [I_2-II-5]. Note that if `K` and `K'` are
       Kodaira symbols for singular fibers of elliptic curves, [K-K'-m]
       and [K'-K-m] are the same type. Finally, [K-K'-1] (not the same as
       [K-K'-1]) is [K'-K-alpha] in the notation of Namikawa-Ueno. The
       figure [2I_0-m] in Namikawa-Ueno, page 159 must be denoted by
       [2I_0-(m+1)].

    #. The second datum is the GROUP OF CONNECTED COMPONENTS (over an
       ALGEBRAIC CLOSURE (!) of `\GF{p}`) of the Neron
       model of J(C). The symbol (n) means the cyclic group with n
       elements. When n=0, (0) is the trivial group (1).
       ``Hn`` is isomorphic to (2)x(2) if n is even and to (4)
       otherwise.

       Note - The set of rational points of `\Phi` can be computed
       using Theorem 1.17 in S. Bosch and Q. Liu "Rational points of the
       group of components of a Neron model", Manuscripta Math. 98 (1999),
       275-293.

    #. Finally, `f` is the exponent of the conductor of
       `J(C)` at `p`.


    .. warning::

       Be careful regarding the formula:

       .. MATH::

          \text{valuation of the naive minimal discriminant} = f + n - 1 + 11c(X).

       (Q. Liu : "Conducteur et discriminant minimal de courbes de genre
       2", Compositio Math. 94 (1994) 51-79, Theoreme 2) is valid only if
       the residual field is algebraically closed as stated in the paper.
       So this equality does not hold in general over
       `\QQ_p`. The fact is that the minimal discriminant
       may change after unramified extension. One can show however that,
       at worst, the change will stabilize after a quadratic unramified
       extension (Q. Liu : "Modèles entiers de courbes hyperelliptiques
       sur un corps de valuation discrète", Trans. AMS 348 (1996),
       4577-4610, Section 7.2, Proposition 4).
    """
    def __init__(self, pari_result, P, Q, Pmin, Qmin, minimal_disc,
                 local_data, conductor):
        self.pari_result = pari_result
        self.P = P
        self.Q = Q
        self.Pmin = Pmin
        self.Qmin = Qmin
        self.minimal_disc = minimal_disc
        self.local_data = local_data
        self.conductor = conductor

    def _repr_(self):
        if self.Q == 0:
            yterm = ''
        else:
            yterm = '+ (%s)*y ' % self.Q

        s = 'Reduction data about this proper smooth genus 2 curve:\n'
        s += '\ty^2 %s= %s\n' % (yterm, self.P)
        if self.Qmin:
            s += 'A Minimal Equation:\n\ty^2 + (%s)y = %s\n' % (self.Qmin, self.Pmin)
        else:
            s += 'A Minimal Equation:\n\ty^2 = %s\n' % self.Pmin
        s += 'Minimal Discriminant: %s\n' % self.minimal_disc
        s += 'Conductor: %s\n' % self.conductor
        s += 'Local Data:\n%s' % self._local_data_str()
        return s

    def _local_data_str(self):
        s = ''
        D = self.local_data
        K = sorted(D.keys())
        for p in K:
            s += 'p=%s\n%s\n' % (p, D[p])
        s = '\t' + '\n\t'.join(s.strip().split('\n'))
        return s


def divisors_to_string(divs):
    """
    Convert a list of numbers (representing the orders of cyclic groups
    in the factorization of a finite abelian group) to a string
    according to the format shown in the examples.

    INPUT:

    - ``divs`` -- a (possibly empty) list of numbers

    OUTPUT: string representation of these numbers

    EXAMPLES::

        sage: from sage.interfaces.genus2reduction import divisors_to_string
        sage: print(divisors_to_string([]))
        (1)
        sage: print(divisors_to_string([5]))
        (5)
        sage: print(divisors_to_string([5]*6))
        (5)^6
        sage: print(divisors_to_string([2,3,4]))
        (2)x(3)x(4)
        sage: print(divisors_to_string([6,2,2]))
        (6)x(2)^2
    """
    s = ""
    n = 0  # How many times have we seen the current divisor?
    for i in range(len(divs)):
        n += 1
        if i+1 == len(divs) or divs[i+1] != divs[i]:
            # Next divisor is different or we are done? Print current one
            if s:
                s += "x"
            s += "(%s)" % divs[i]
            if n > 1:
                s += "^%s" % n
            n = 0

    return s or "(1)"


class Genus2reduction(SageObject):
    r"""
    Conductor and Reduction Types for Genus 2 Curves.

    Use ``R = genus2reduction(Q, P)`` to obtain reduction
    information about the Jacobian of the projective smooth curve
    defined by `y^2 + Q(x)y = P(x)`. Type ``R?``
    for further documentation and a description of how to interpret the
    local reduction data.

    EXAMPLES::

        sage: x = QQ['x'].0
        sage: R = genus2reduction(x^3 - 2*x^2 - 2*x + 1, -5*x^5)
        sage: R.conductor
        1416875
        sage: factor(R.conductor)
        5^4 * 2267

    The discriminant is always minimal::

        sage: factor(R.minimal_disc)
        2^3 * 5^5 * 2267

    Printing R summarizes all the information computed about the curve

    ::

        sage: R
        Reduction data about this proper smooth genus 2 curve:
            y^2 + (x^3 - 2*x^2 - 2*x + 1)*y = -5*x^5
        A Minimal Equation:
            y^2 ...
        Minimal Discriminant: 56675000
        Conductor: 1416875
        Local Data:
            p=2
            (potential) stable reduction:  (II), j=1
            p=5
            (potential) stable reduction:  (I)
            reduction at p: [V] page 156, (3), f=4
            p=2267
            (potential) stable reduction:  (II), j=432
            reduction at p: [I{1-0-0}] page 170, (1), f=1

    Here are some examples of curves with modular Jacobians::

        sage: R = genus2reduction(x^3 + x + 1, -2*x^5 - 3*x^2 + 2*x - 2)
        sage: factor(R.conductor)
        23^2
        sage: factor(genus2reduction(x^3 + 1, -x^5 - 3*x^4 + 2*x^2 + 2*x - 2).conductor)
        29^2
        sage: factor(genus2reduction(x^3 + x + 1, x^5 + 2*x^4 + 2*x^3 + x^2 - x - 1).conductor)
        5^6

    EXAMPLES::

        sage: genus2reduction(0, x^6 + 3*x^3 + 63)
        Reduction data about this proper smooth genus 2 curve:
                y^2 = x^6 + 3*x^3 + 63
        A Minimal Equation:
                y^2 ...
        Minimal Discriminant: -10628388316852992
        Conductor: 2893401
        Local Data:
                p=2
                (potential) stable reduction:  (V), j1+j2=0, j1*j2=0
                p=3
                (potential) stable reduction:  (I)
                reduction at p: [III{9}] page 184, (3)^2, f=10
                p=7
                (potential) stable reduction:  (V), j1+j2=0, j1*j2=0
                reduction at p: [I{0}-II-0] page 159, (1), f=2

    In the above example, Liu remarks that in fact at `p=2`,
    the reduction is [II-II-0] page 163, (1), `f=8`. So the
    conductor of J(C) is actually `2 \cdot 2893401=5786802`.

    A MODULAR CURVE:

    Consider the modular curve `X_1(13)` defined by an
    equation

    .. MATH::

                   y^2 + (x^3-x^2-1)y = x^2 - x.

    We have::

        sage: genus2reduction(x^3-x^2-1, x^2 - x)
        Reduction data about this proper smooth genus 2 curve:
                y^2 + (x^3 - x^2 - 1)*y = x^2 - x
        A Minimal Equation:
                y^2 ...
        Minimal Discriminant: -169
        Conductor: 169
        Local Data:
                p=13
                (potential) stable reduction:  (V), j1+j2=0, j1*j2=0
                reduction at p: [I{0}-II-0] page 159, (1), f=2

    So the curve has good reduction at 2. At `p=13`, the stable
    reduction is union of two elliptic curves, and both of them have 0
    as modular invariant. The reduction at 13 is of type [I_0-II-0]
    (see Namikawa-Ueno, page 159). It is an elliptic curve with a cusp.
    The group of connected components of the Neron model of
    `J(C)` is trivial, and the exponent of the conductor of
    `J(C)` at `13` is `f=2`. The conductor of
    `J(C)` is `13^2`. (Note: It is a theorem of
    Conrad-Edixhoven-Stein that the component group of
    `J(X_1(p))` is trivial for all primes `p`.)
    """
    def __init__(self):
        pass

    def _repr_(self):
        """
        EXAMPLES::

            sage: genus2reduction
            Genus 2 reduction PARI interface
        """
        return "Genus 2 reduction PARI interface"

    def __call__(self, Q, P):
        """
        Compute and return the :class:`ReductionData` corresponding to
        the genus 2 curve `y^2 + Q(x) y = P(x)`.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: genus2reduction(x^3 - 2*x^2 - 2*x + 1, -5*x^5)
            Reduction data about this proper smooth genus 2 curve:
                    y^2 + (x^3 - 2*x^2 - 2*x + 1)*y = -5*x^5
            A Minimal Equation:
                    y^2 ...
            Minimal Discriminant: 56675000
            Conductor: 1416875
            Local Data:
                    p=2
                    (potential) stable reduction:  (II), j=1
                    p=5
                    (potential) stable reduction:  (I)
                    reduction at p: [V] page 156, (3), f=4
                    p=2267
                    (potential) stable reduction:  (II), j=432
                    reduction at p: [I{1-0-0}] page 170, (1), f=1

        ::

            sage: genus2reduction(x^2 + 1, -5*x^5)
            Reduction data about this proper smooth genus 2 curve:
                    y^2 + (x^2 + 1)*y = -5*x^5
            A Minimal Equation:
                    y^2 ...
            Minimal Discriminant: 48838125
            Conductor: 32025
            Local Data:
                    p=3
                    (potential) stable reduction:  (II), j=1
                    reduction at p: [I{1-0-0}] page 170, (1), f=1
                    p=5
                    (potential) stable reduction:  (IV)
                    reduction at p: [I{1-1-2}] page 182, (5), f=2
                    p=7
                    (potential) stable reduction:  (II), j=4
                    reduction at p: [I{1-0-0}] page 170, (1), f=1
                    p=61
                    (potential) stable reduction:  (II), j=57
                    reduction at p: [I{2-0-0}] page 170, (2), f=1

        Verify that we fix :issue:`5573`::

            sage: genus2reduction(x^3 + x^2 + x,-2*x^5 + 3*x^4 - x^3 - x^2 - 6*x - 2)
            Reduction data about this proper smooth genus 2 curve:
                    y^2 + (x^3 + x^2 + x)*y = -2*x^5 + 3*x^4 - x^3 - x^2 - 6*x - 2
            A Minimal Equation:
                    y^2 ...
            Minimal Discriminant: 1520984142
            Conductor: 954
            Local Data:
                    p=2
                    (potential) stable reduction:  (II), j=1
                    reduction at p: [I{1-0-0}] page 170, (1), f=1
                    p=3
                    (potential) stable reduction:  (I)
                    reduction at p: [II] page 155, (1), f=2
                    p=53
                    (potential) stable reduction:  (II), j=12
                    reduction at p: [I{1-0-0}] page 170, (1), f=1
        """
        R = PolynomialRing(QQ, 'x')
        P = R(P)
        Q = R(Q)
        if P.degree() > 6:
            raise ValueError("P (=%s) must have degree at most 6" % P)
        if Q.degree() >= 4:
            raise ValueError("Q (=%s) must have degree at most 3" % Q)

        res = pari.genus2red([P, Q])
        conductor = ZZ(res[0])
        Pmin = R(res[2][0])
        Qmin = R(res[2][1])
        minimal_disc = ZZ(pari.hyperelldisc(res[2]))

        local_data = {}
        for red in res[3]:
            p = ZZ(red[0])

            t = red[1]
            data = "(potential) stable reduction:  (%s)" % roman_numeral[int(t[0])]
            t = t[1]
            if len(t) == 1:
                data += ", j=%s" % t[0].lift()
            elif len(t) == 2:
                data += ", j1+j2=%s, j1*j2=%s" % (t[0].lift(), t[1].lift())

            t = red[2]
            if t:
                data += "\nreduction at p: %s, " % str(t[0]).replace('"', '').replace("(tame) ", "")
                data += divisors_to_string(t[1]) + ", f=" + str(res[0].valuation(red[0]))

            local_data[p] = data

        return ReductionData(res, P, Q, Pmin, Qmin, minimal_disc, local_data, conductor)

    def __reduce__(self):
        return _reduce_load_genus2reduction, tuple([])


# An instance
genus2reduction = Genus2reduction()


def _reduce_load_genus2reduction():
    return genus2reduction
