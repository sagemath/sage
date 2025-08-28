# -*- coding: utf-8 -*-
r"""
Finite Ring Extensions

Let `f: A\rightarrow  B` be a finite ring extension, given as a ring map
between (quotients of) multivariate polynomial rings.

Let `m_1, \dots, m_s` be elements of `B` that generate `B` as an `A`-module.
Consider the surjective ring morphism `\psi:A[z_1, ..., z_s]\rightarrow B`
defined by `f` and evaluation on `m_1, \dots, m_s`.

Then the class FiniteRingExtension returns the quotient ring `C/I` with
`C=A[z_1, ..., z_s]` and `I=Ker(\psi)` together with a ring morphism
`\phi:B\rightarrow C/I`, inverse to the morphism `\psi:C/I\rightarrow B`
induced and still denoted by `\psi`.

Moreover, the annihilator `Ann_A(B)` of `B` as an `A`-module is calculated
as well as a presentation matrix `P` of `B` as an `A`-module:

.. math::

    A^p\ \xrightarrow{P}\ A^s\ \xrightarrow{m}\ B\ \longrightarrow\ 0

EXAMPLE:

Consider the twisted cubic in `\mathbb{P}^3`. On the level of ChowRings the
inclusion `i` induces the ring map

.. math::

    A^{*}\mathbb{P}^3\ \xrightarrow{\ j\ }\ A^{*}\mathbb{P}^1.

::

    sage: A.<h> = ChowRing('h', 1, 'h^4')  # P3
    sage: B.<w> = ChowRing('w', 1, 'w^2')  # P1
    sage: j = A.hom([3*w], B)  # j=i^{*} sends h to 3*w.
    sage: F = FiniteRingExtension(j); F
    Quotient of Multivariate Polynomial Ring in z, h over Rational Field by the ideal (h^2, z - 1)
    sage: F.psi()
    Ring morphism:
      From: Quotient of Multivariate Polynomial Ring in z, h over Rational Field by the ideal (h^2, z - 1)
      To:   Quotient of Multivariate Polynomial Ring in w over Rational Field by the ideal (w^2)
      Defn: 1 |--> 1
            h |--> 3*w
    sage: F.phi()
    Ring morphism:
      From: Quotient of Multivariate Polynomial Ring in w over Rational Field by the ideal (w^2)
      To:   Quotient of Multivariate Polynomial Ring in z, h over Rational Field by the ideal (h^2, z - 1)
      Defn: w |--> 1/3*h
    sage: F.ann()
    Ideal (h^2) of Quotient of Multivariate Polynomial Ring in h over Rational Field by the ideal (h^4)
    sage: F.prm()  # presentation matrix
    [0]
    sage: F.mgs()  # module generators
    (1,)
    sage: F.mds()  # module generator degrees
    (0,)
    sage: F.nvs()  # new variables that have been introduced
    ('z',)

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""


# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import QQ
from sage.matrix.all import matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.quotient_ring import QuotientRing_generic
from sage.rings.polynomial.term_order import TermOrder
from sage.libs.singular.function_factory import singular_function, lib as singular_lib
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base

from sage.schemes.chow.ring import Kernel2


class FiniteRingExtension(QuotientRing_generic):

    def __init__(self, f, var_name='z'):
        r"""
        Construct a :class:`FiniteRingExtension`.

        INPUT:

        - ``f`` -- a finite ring morphism between (quotients of) multivariate polynomial rings.

        OUTPUT:

        - :class:`FiniteRingExtension <QuotientRing_generic>`.

        EXAMPLES:

        Consider the morphism on the level of Chow rings induced by the Veronese
        embedding (http://en.wikipedia.org/wiki/Veronese_surface)
        `i:\mathbb{P}^2\longrightarrow\mathbb{P}^5`::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f)
            sage: F.cover_ring().gens()
            (z, y)
            sage: F.defining_ideal().gens()
            [y^3, z - 1]

        Another examples is the morphism on the level of Chow rings induced by
        the Segre embedding
        `\mathbb{P}^2\times\mathbb{P}^2\longrightarrow\mathbb{P}^8`::

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f);
            sage: F.cover_ring().gens()
            (z1, z2, z3, l)
            sage: F.defining_ideal().gens()
            [l^5, z3 - 1, 2*z2*l^3 - l^4, 3*z1*l - 3*z2*l^2 + l^3, z2^2 - z1, z1*z2, z1^2]

        TESTS:

        We test a corner case with codomain a quotient of a multivariate
        polynomial ring in no variables (used while blowing up `\mathbb{P}^2`
        in a point)::

            sage: A.<h> = ChowRing('h', 1, 'h^3')
            sage: B = ChowRing()
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.cover_ring().gens()
            (z, h)
            sage: F.defining_ideal().gens()
            [h, z - 1]

        Finally run the test suite on F::

            sage: TestSuite(F).run()
        """

        #######################################################################
        # Options
        #######################################################################

        if not isinstance(var_name, str):
            m = "String as var_name in FiniteRingExtension expected."
            raise ValueError(m, var_name)

        #######################################################################
        # Singular functions we need...
        #######################################################################

        singular_lib('algebra.lib')
        deg = singular_function('deg')
        sing_reduce = singular_function('reduce')
        division = singular_function('division')
        std = singular_function('std')
        kbase = singular_function('kbase')

        #######################################################################
        # Get rings A, B from f such that f: A -> B and lift to cover rings:
        # ff: AA -> BB with A = AA/AI and B = BB/BI
        #######################################################################

        A, B = f.domain(), f.codomain()
        ff = f  # A priori A, B are MPolynomial.
        # Check the ring A and get the rings AA and AI:
        if isinstance(A, MPolynomialRing_base):
            AA = A
        elif isinstance(A, QuotientRing_generic):
            AA = A.cover_ring()
            ff = f.morphism_from_cover()
        else:
            err = "Domain is not (quotient of) Multivariate Polynomial Ring."
            raise TypeError(err)
        # Check the ring B and get the rings BB and BI:
        if isinstance(B, MPolynomialRing_base):
            BB, BI = B, B.ideal(0)
        elif isinstance(B, QuotientRing_generic):
            BB, BI = B.cover_ring(), B.defining_ideal()
            ff = AA.hom([BB(str(x.lift())) for x in ff.im_gens()], BB)
        else:
            err = "Codomain is not (quotient of) Multivariate Polynomial Ring."
            raise TypeError(err)

        #######################################################################
        # Get the module generators m_1,..., m_s of B as an A-module.
        # In order to do this, we calculate a kbase of the ideal given by the
        # image of A under f.
        # Remark: as sage.libs.singular.function only works for polynomial
        # rings (and not quotients of) we use the above lifts then reduce to B.
        #######################################################################

        if B.ngens() == 0:
            # As B = Q in this case, we know the module_gens and module_degs.
            module_gens = [B(1)]
            module_degs = [0]
        else:
            # Compute using singular: kbase
            II = std(BI.gens() + ff.im_gens(), ring=BB)
            kbasis = kbase(II, attributes={II: {'isSB': 1}}, ring=BB)
            module_gens = [B(b) for b in kbasis]
            if module_gens == [B(0)]:
                raise ValueError("Extension is not finite.")
            module_degs = [deg(b, ring=BB) for b in kbasis]

        self._mgs = tuple(module_gens)
        self._mds = tuple(module_degs)
        s = len(module_gens)

        #######################################################################
        # Get the variable names from A, B and their number.
        # Remark: in Singular the variable names in a quotient are the same
        # as opposed to sage which by default adds an 'bar'.
        # To ensure singulars behavior we get the variable names of AA.
        #######################################################################

        a_vars, a = [str(v) for v in AA.gens()], len(AA.gens())
        b_vars, b = [str(v) for v in BB.gens()], len(BB.gens())

        #######################################################################
        # Construct new variables
        #               z_1, ..., z_s corresponding to module_gens
        # as well as    z_{s+1}, ..., z_{s+b} corresponding to b_vars
        #######################################################################

        t_vars = ['%s%d' % (var_name, i) for i in range(1, s + 1)]
        t_vars = [var_name] if s == 1 else t_vars
        n_vars = ['%s%d' % (var_name, i) for i in range(s + 1, s + b + 1)]
        # Ensure that the variable 'z' does not occur in the variables of A.
        while [c for c in n_vars if c in a_vars]:
            var_name += var_name
            t_vars = ['%s%d' % (var_name, i) for i in range(1, s + 1)]
            t_vars = [var_name] if s == 1 else t_vars
            n_vars = ['%s%d' % (var_name, i) for i in range(s + 1, s + b + 1)]
        self._nvs = tuple(t_vars)

        #######################################################################
        # Construct an intermediate ring R defined as follows:
        #   R = Q[z_{s+1}, ..., z_{b+s}, z_{1}, ..., z_{s}, a_1, ..., a_a]
        # and a map chi: R -> B given by B.gens(), module_gens, f.im_gens()
        # then compute kern = kernel(chi)
        #######################################################################

        r_vars = n_vars + t_vars + a_vars
        r_term_order = TermOrder('wdegrevlex', (1,) * s) + AA.term_order()
        if B.ngens():
            # Add the order for B if B has generators
            # noinspection PyAugmentAssignment
            r_term_order = TermOrder('wdegrevlex', (1,) * b) + r_term_order

        R = PolynomialRing(QQ, len(r_vars), names=r_vars, order=r_term_order)
        chi_im_gens = list(B.gens()) + module_gens + ff.im_gens()

        if B.ngens():
            # Compute the kernel using alg_kernel
            chi = R.hom(chi_im_gens, BB.quo(BI, names=b_vars))
            kern = Kernel2(chi)
        else:
            # The kernel is simply given by those elements with trivial image.
            # Keeping attention on the order though
            z1 = [R(AA.gen(i)) for i in range(a) if B(ff.im_gens()[i]) == B(0)]
            z2 = [R.gen(i) - R(module_gens[i]) for i in range(s)]
            kern = R.ideal(z1 + z2)
        kern_gens, l = kern.gens(), len(kern.gens())

        #######################################################################
        # Now we construct the ring C = A[z_1, ..., z_s]/I and the morphisms
        # psi: A[z_1,...,z_s]/I -> B given by f and m and its inverse
        # phi: B -> A[z_1,...,z_s]/I.
        # Remark: we keep z_1, ..., z_s as variable names in the quotient.
        #######################################################################
        # The ring C = A[z_1,...,z_s]/I
        # TODO: discuss about module_degs as degrees for t_vars

        c_vars = t_vars + a_vars
        c_term_order = TermOrder('wdegrevlex', (1,) * s) + A.term_order()
        CC = PolynomialRing(QQ, len(c_vars), names=c_vars, order=c_term_order)
        relations = CC.ideal([CC(str(x)) for x in kern_gens[0:l - b]])
        CI = CC.ideal(std(relations))
        # The morphism psi: A[z_1,...,z_s]/I -> B
        self._B = B
        self._psi_gens = module_gens + ff.im_gens()
        # The morphism phi: B -> A[z_1,...,z_s]/I
        if b:
            self._phi_gens = sing_reduce(R.ideal([R.gen(i) for i in range(b)]),
                                         kern, ring=R,
                                         attributes={kern: {'isSB': 1}})
        else:
            self._phi_gens = []

        #######################################################################
        # Finally we calculate the Annihilator and the presentation matrix
        #######################################################################

        weights = (1,) * (b + s) + (0,) * a
        kern_degs = [deg(x, weights, ring=R) for x in kern_gens]
        i = 0
        l = len(kern_gens)
        for d in kern_degs:
            if d != 0:
                break
            else:
                i += 1
        # Ann_A(B) is given by the first i kernel generators
        if i:
            self._ann = A.ideal([A(str(x)) for x in kern_gens[0:i]])
        else:
            self._ann = A.ideal(0)
        j = 0
        for d in kern_degs[i:l - b]:
            if d != 1:
                break
            else:
                j += 1
        # The presentation matrix prm
        IX = R.ideal(kern_gens[i + 1:i + j])
        L = division(R.ideal(kern_gens[i + 1:i + j]), R.ideal(t_vars), ring=R)
        P = matrix(R, IX.ngens(), s, (L[0]).transpose())
        Q = matrix(R, IX.ngens(), 1, L[1])
        M = matrix(R, 1, s)
        M[0, s - 1] = 1
        self._prm = (P + (Q * M)).transpose().apply_map(lambda x: A(str(x)), A)
        QuotientRing_generic.__init__(self, CC, CI, c_vars)

    def ann(self):
        r"""
        Return the annihilator `Ann_A(B)` of `B` seen as an `A`-module via `f`.

        EXAMPLES::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f)
            sage: F.ann().gens()
            [y^3]

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f)
            sage: F.ann().gens()
            [l^5]

        TESTS::

            sage: A.<h> = ChowRing('h', 1, 'h^3')
            sage: B = ChowRing()
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.ann().gens()
            [h]
        """
        return self._ann

    def nvs(self):
        r"""
        Return the tuple of 'new' variables `(z_1, \dots, z_s)` corresponding
        to the module generators of the `A`-module `B` (via `f`).

        EXAMPLES::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f, var_name='t')
            sage: F.nvs()
            ('t',)

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f)
            sage: F.nvs()
            ('z1', 'z2', 'z3')

        TESTS::

            sage: A.<h> = ChowRing('h', 1, 'h^3')
            sage: B = ChowRing()
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.nvs()
            ('z',)
        """
        return self._nvs

    def prm(self):
        r"""
        Return the presentation matrix.

        .. math::

            A^{*}\mathbb{P}^3\ \xrightarrow{\ j\ }\ A^{*}\mathbb{P}^1.

        EXAMPLES::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f)
            sage: F.prm()
            [0]

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f)
            sage: F.prm()
            [     0    3*l]
            [ 2*l^3 -3*l^2]
            [  -l^4    l^3]

        TESTS::

            sage: A.<h> = ChowRing('h', 1, 'h^3')
            sage: B = ChowRing()
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.prm()
            [0]
        """
        return self._prm

    def mgs(self):
        """
        Return the tuple of module generators of `B` seen as an `A`-module
        via `f`.

        EXAMPLES::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f)
            sage: F.mgs()
            (1,)

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f)
            sage: F.mgs()
            (k^2, k, 1)

        TESTS::

            sage: A.<h> = ChowRing('h', 1, 'h^3')
            sage: B = ChowRing()
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.mgs()
            (1,)
        """
        return self._mgs

    def mds(self):
        """
        Return the tuple of the degrees of the module generators of `B`
        seen as an `A`-module via `f`.

        EXAMPLES::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f)
            sage: F.mds()
            (0,)

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f)
            sage: F.mds()
            (2, 1, 0)

        TESTS::

            sage: B = ChowRing()
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.mds()
            (0,)
        """
        return self._mds

    def push_down(self, v):
        """
        Return the push_down of a sequence of elements in B.

        INPUT:

        - ``v``-- a list or tuple of elements in B

        EXAMPLES::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f)
            sage: F.push_down([x^2])
            [1/4*y^2]

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f)
            sage: F.push_down((h, k))
            [ 0  0]
            [-1  1]
            [ l  0]

        TESTS::

            sage: A.<h> = ChowRing('h', 1, 'h^3')
            sage: B = ChowRing()
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.push_down([1])
            [1]
        """
        # Singular functions we need
        sing_division = singular_function('division')
        sing_reduce = singular_function('reduce')
        # Get ideal in self and sizes
        w = [self.phi()(y) for y in v]
        s, t = len(self.nvs()), len(w)
        # Lift
        CC = self.cover_ring()
        CI = self.defining_ideal()
        lft = [x.lift() for x in w]
        # Compute division
        if lft == [CC(0)] * t:
            P = matrix(self, t, s, lft * s)
            Q = matrix(self, t, 1, lft)
        else:
            JJ = CC.ideal(lft)
            JJ = sing_reduce(JJ, CI, ring=CC, attributes={CI: {'isSB': 1}})
            L = sing_division(CC.ideal(JJ), CC.ideal(self.nvs()), ring=CC)
            P = matrix(self, t, s, (L[0]).transpose())
            Q = matrix(self, t, 1, L[1])
        M = matrix(self, 1, s)
        M[0, s - 1] = 1

        return ((P + (Q * M)).transpose()).apply_map(lambda x: self(str(x)))

    def psi(self):
        r"""
        Return he isomorphism `\psi:C/I\longrightarrow B` given by `f` and the
        module generators.

        EXAMPLES::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f)
            sage: F.psi()
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in z, y over Rational Field by the ideal (y^3, z - 1)
              To:   Quotient of Multivariate Polynomial Ring in x over Rational Field by the ideal (x^3)
              Defn: 1 |--> 1
                    y |--> 2*x

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f)
            sage: F.psi()
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in z1, z2, z3, l over Rational Field by the ideal (l^5, z3 - 1, 2*z2*l^3 - l^4, 3*z1*l - 3*z2*l^2 + l^3, z2^2 - z1, z1*z2, z1^2)
              To:   Quotient of Multivariate Polynomial Ring in h, k over Rational Field by the ideal (h^3, k^3)
              Defn: z1 |--> k^2
                    z2 |--> k
                    1 |--> 1
                    l |--> h + k

        TESTS::

            sage: A.<h> = ChowRing('h', 1, 'h^3')
            sage: B = ChowRing(0)
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.psi()
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in z, h over Rational Field by the ideal (h, z - 1)
              To:   Quotient of Multivariate Polynomial Ring in no variables over Rational Field by the ideal (0)
              Defn: 1 |--> 1
                    0 |--> 0
        """
        return self.hom(self._psi_gens, self._B)

    def phi(self):
        r"""
        Return the isomorphism `\phi:B\longrightarrow C/I` inverse to `\psi`.

        EXAMPLES::

            sage: A.<y> = ChowRing('y', 1, 'y^6')  # P5
            sage: B.<x> = ChowRing('x', 1, 'x^3')  # P2
            sage: f = A.hom([2*x], B)
            sage: F = FiniteRingExtension(f)
            sage: F.phi()
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in x over Rational Field by the ideal (x^3)
              To:   Quotient of Multivariate Polynomial Ring in z, y over Rational Field by the ideal (y^3, z - 1)
              Defn: x |--> 1/2*y

            sage: A.<l> = ChowRing('l', 1, 'l^9')
            sage: B.<h,k> = ChowRing(['h', 'k'], [1, 1], ['h^3', 'k^3'])
            sage: f = A.hom([h+k], B)  # Segre embedding P2 x P2 to P8
            sage: F = FiniteRingExtension(f)
            sage: F.phi()
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in h, k over Rational Field by the ideal (h^3, k^3)
              To:   Quotient of Multivariate Polynomial Ring in z1, z2, z3, l over Rational Field by the ideal (l^5, z3 - 1, 2*z2*l^3 - l^4, 3*z1*l - 3*z2*l^2 + l^3, z2^2 - z1, z1*z2, z1^2)
              Defn: h |--> -z2 + l
                    k |--> z2

        TESTS::

            sage: A.<h> = ChowRing('h', 1, 'h^3')
            sage: B = ChowRing()
            sage: f = A.hom([B(0)], B)
            sage: F = FiniteRingExtension(f)
            sage: F.phi()
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in no variables over Rational Field by the ideal (0)
              To:   Quotient of Multivariate Polynomial Ring in z, h over Rational Field by the ideal (h, z - 1)
        """
        return self._B.hom([self(str(x)) for x in self._phi_gens], self)
