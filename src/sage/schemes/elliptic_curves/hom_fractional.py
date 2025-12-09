r"""
Fractional morphisms of elliptic curves

Algorithms involving advanced computations with endomorphisms
or isogenies of elliptic curves sometimes involve divisions of
elliptic-curve morphisms by integers.
This operation yields another well-defined morphism whenever the
kernel of the numerator contains the kernel of the denominator.

The class :class:`EllipticCurveHom_fractional` represents symbolic
fractions `\varphi / n` where `\varphi\colon E\to E'` is any
:class:`~sage.schemes.elliptic_curves.hom.EllipticCurveHom` whose
kernel contains the `n`-torsion subgroup of `E`.
Functionality for converting this fraction to a more explicit form
is provided (:meth:`~EllipticCurveHom_fractional.to_isogeny_chain`).

EXAMPLES:

Division by an integer::

    sage: E = EllipticCurve(GF(419), [-1, 0])
    sage: phi = (E.frobenius_isogeny() + 1) / 2
    sage: phi
    Fractional elliptic-curve morphism of degree 105:
      Numerator:   Sum morphism:
        From: Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419
        To:   Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419
        Via:  (Frobenius endomorphism of degree 419:
            From: Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419
            To:   Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419,
          Scalar-multiplication endomorphism [1]
            of Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419)
      Denominator: 2
    sage: phi.degree()
    105
    sage: phi.to_isogeny_chain()
    Composite morphism of degree 105 = 1*3*5*7:
      From: Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419
      To:   Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419

Right division of isogenies::

    sage: E = EllipticCurve(GF(419), [1, 0])
    sage: ker = E(125, 70)
    sage: phi = E.isogeny(ker); phi
    Isogeny of degree 35
      from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 419
      to Elliptic Curve defined by y^2 = x^3 + 289*x + 323 over Finite Field of size 419
    sage: psi = E.isogeny(5*ker); psi
    Isogeny of degree 7
      from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 419
      to Elliptic Curve defined by y^2 = x^3 + 285*x + 87 over Finite Field of size 419
    sage: chi = phi.divide_right(psi); chi
    Fractional elliptic-curve morphism of degree 5:
      Numerator:   Composite morphism of degree 245 = 7*35:
      From: Elliptic Curve defined by y^2 = x^3 + 285*x + 87 over Finite Field of size 419
      To:   Elliptic Curve defined by y^2 = x^3 + 289*x + 323 over Finite Field of size 419
      Denominator: 7
    sage: phi == chi * psi
    True

Left division of isogenies::

    sage: E = EllipticCurve(GF(419), [1, 0])
    sage: ker = E(125, 70)
    sage: phi = E.isogeny(ker); phi
    Isogeny of degree 35
      from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 419
      to Elliptic Curve defined by y^2 = x^3 + 289*x + 323 over Finite Field of size 419
    sage: tmp = E.isogeny(7*ker)
    sage: psi = tmp.codomain().isogeny(tmp(ker)); psi
    Isogeny of degree 7
      from Elliptic Curve defined by y^2 = x^3 + 269*x + 82 over Finite Field of size 419
      to Elliptic Curve defined by y^2 = x^3 + 289*x + 323 over Finite Field of size 419
    sage: chi = phi.divide_left(psi); chi
    Fractional elliptic-curve morphism of degree 5:
      Numerator:   Composite morphism of degree 245 = 35*7:
      From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 419
      To:   Elliptic Curve defined by y^2 = x^3 + 269*x + 82 over Finite Field of size 419
      Denominator: 7
    sage: phi == psi * chi
    True

AUTHORS:

- Lorenz Panny (2024)
"""
from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

from sage.rings.integer_ring import ZZ

from sage.schemes.elliptic_curves.hom import EllipticCurveHom


class EllipticCurveHom_fractional(EllipticCurveHom):
    r"""
    This class represents a (symbolic) quotient of an isogeny divided by an integer.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_fractional import EllipticCurveHom_fractional
        sage: phi = EllipticCurve([1,1]).scalar_multiplication(-2)
        sage: EllipticCurveHom_fractional(phi, 2)
        Fractional elliptic-curve morphism of degree 1:
          Numerator:   Scalar-multiplication endomorphism [-2]
                         of Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
          Denominator: 2
        sage: EllipticCurveHom_fractional(phi, 3)
        Traceback (most recent call last):
        ...
        ValueError: Scalar-multiplication endomorphism [-2]
          of Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
          is not divisible by 3
    """

    def __init__(self, phi, d, *, check=True) -> None:
        r"""
        Construct a (symbolic) quotient of an isogeny divided by an integer.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_fractional import EllipticCurveHom_fractional
            sage: phi = EllipticCurve(GF(11), [1,1]).scalar_multiplication(-3)
            sage: EllipticCurveHom_fractional(phi, 3)
            Fractional elliptic-curve morphism of degree 1:
              Numerator:   Scalar-multiplication endomorphism [-3]
                             of Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 11
              Denominator: 3
        """
        if not isinstance(phi, EllipticCurveHom):
            raise TypeError('not an elliptic-curve morphism')

        d = ZZ(d)
        if d < 0:
            phi = -phi
            d = -d
        if not d:
            raise ZeroDivisionError('cannot divide isogeny by zero')

        if check:
            if not (d**2).divides(phi.degree()):
                raise ValueError(f'{phi} is not divisible by {d}')

            E = phi.domain()
            for l, e in d.factor():
                F = E.division_field(l**e)
                EE = E.change_ring(F)

                psi = E.division_polynomial(l**e).radical()
                psi //= psi.gcd(E.division_polynomial(l**(e-1)))
                if psi.degree() < l**2//2:
                    assert l == E.base_field().characteristic()
                    if psi.is_one():
                        # supersingular, hence [p] is the only p^2-isogeny up to isomorphism
                        continue
                    # ordinary, hence [p] is Frobenius times its dual
                    insep = phi.inseparable_degree().valuation(l)
                    sep = phi.degree().valuation(l) - insep
                    if sep < e or insep < e:
                        raise ValueError(f'{phi} is not divisible by {l**e}')
                    continue

                P, Q = _torsion_gens(E, EE, l, e, psi=psi)

                if phi._eval(P) or phi._eval(Q):
                    raise ValueError(f'{phi} is not divisible by {l**e}')

        self._phi = phi
        self._d = d
        self._degree = self._phi.degree() // d**2
        self._domain = phi.domain()
        self._codomain = phi.codomain()
        EllipticCurveHom.__init__(self, self._domain, self._codomain)

    def _call_(self, P):
        r"""
        Evaluate this fractional elliptic-curve morphism at a point.

        EXAMPLES::

            sage: E = EllipticCurve(GF(13), [9,0])
            sage: pi = E.frobenius_isogeny()
            sage: phi = (1 + pi) / 2
            sage: phi.degree()
            2
            sage: P, Q = E(6,7), E(2, 0)
            sage: phi(P)
            (4 : 3 : 1)
            sage: phi(Q)
            (0 : 0 : 1)
        """
        return self._eval(P)

    def _eval(self, P):
        r"""
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self`` at the given point

        EXAMPLES::

            sage: E = EllipticCurve(GF(13), [9,0])
            sage: pi = E.frobenius_isogeny()
            sage: phi = (1 + pi) / 2
            sage: phi.degree()
            2
            sage: EE = E.change_ring(GF(13^3))
            sage: EE.base_field().inject_variables()
            Defining z3
            sage: P = EE(4*z3^2 + 11*z3, 8*z3^2 + 11*z3 + 8)
            sage: Q = EE(z3^2 + 7*z3 + 6, 3*z3^2 + 7*z3 + 10)
            sage: phi._eval(P)
            (9*z3^2 + 6*z3 + 6 : 4*z3^2 + 9*z3 + 3 : 1)
            sage: phi._eval(Q)
            (z3^2 + 9*z3 : 10*z3^2 + 6*z3 + 10 : 1)
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')
        k = Sequence(P).universe()

        if not P:
            return self._codomain.base_extend(k).zero()

        # TODO this should really be a "divide point by possibly
        # extending the base field" method
        F = k
        n = P.order()
        m = self._d.prime_to_m_part(n)
        P *= m.inverse_mod(n)
        for q, e in (self._d//m).factor():
            for _ in range(e):
                f = P.division_points(q, poly_only=True)
                try:
                    f.any_root(assume_squarefree=True)
                except ValueError:
                    g = f.factor()[0][0]
                    F = F.extension(g.degree())
                    g.any_root(ring=F)
                P = P.change_ring(F).division_points(q)[0]

        Q = self._phi._eval(P).change_ring(k)

        return self._codomain.base_extend(k)(*Q)

    def _repr_(self) -> str:
        r"""
        Return a description of this fractional elliptic-curve morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(13), [9,0])
            sage: pi = E.frobenius_isogeny()
            sage: (1 + pi) / 2
            Fractional elliptic-curve morphism of degree 2:
              Numerator:   Sum morphism:
                From: Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 13
                To:   Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 13
                Via:  (Scalar-multiplication endomorphism [1]
                         of Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 13,
                       Frobenius endomorphism of degree 13:
                         From: Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 13
                         To:   Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 13)
              Denominator: 2
        """
        return f'Fractional elliptic-curve morphism of degree {self._degree}:' \
            f'\n  Numerator:   {self._phi}' \
            f'\n  Denominator: {self._d}'

    @cached_method
    def to_isogeny_chain(self):
        r"""
        Convert this fractional elliptic-curve morphism into a (non-fractional)
        :class:`~sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite`
        object representing the same morphism.

        EXAMPLES::

            sage: p = 419
            sage: E = EllipticCurve(GF(p^2), [1,0])
            sage: iota = E.automorphisms()[2]   # sqrt(-1)
            sage: pi = E.frobenius_isogeny()    # sqrt(-p)
            sage: endo = (iota + pi) / 2
            sage: endo.degree()
            105
            sage: endo.to_isogeny_chain()
            Composite morphism of degree 105 = 1*3*5*7:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 419^2
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 419^2
            sage: endo.to_isogeny_chain() == endo
            True
        """
        E = self._domain

        ker = []
        insep = 0
        for l, e in self._phi.degree().factor():
            F = E.division_field(l**e)
            EE = E.change_ring(F)

            psi = E.division_polynomial(l**e).radical()
            psi //= psi.gcd(E.division_polynomial(l**(e-1)))
            if psi.degree() < l**2//2:
                assert l == E.base_field().characteristic()
                if psi.is_one():
                    # supersingular, hence [p] is the only p^2-isogeny up to isomorphism
                    insep += 2*self._d.valuation(l)
                else:
                    # ordinary, hence [p] is Frobenius times its dual
                    insep += self._d.valuation(l)
                    ker.append(self._d.p_primary_part(l) * EE.lift_x(psi.any_root(ring=F)))
                continue

            P, Q = _torsion_gens(E, EE, l, e, psi=psi)
            if self.is_endomorphism():
                RS = None
            else:
                RS = _torsion_gens(self._codomain, self._codomain.change_ring(F), l, e)

            mat = self._phi.matrix_on_subgroup((P, Q), RS)
            for row in filter(bool, self._d.p_primary_part(l) * mat.left_kernel_matrix()):
                K = sum(ZZ(c)*T for c, T in zip(row, (P, Q)))
                K.set_order(multiple=l**e)
                assert self._eval(K) == 0
                ker.append(K)

        from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
        chain = EllipticCurveHom_composite(E, [])
        ker = ker[::-1]
        while ker:
            if not (P := ker.pop()):
                continue
            (l, e), = P.order().factor()
            K = l**(e-1)*P
            if e > 1:
                ker.append(P)
            poly = E.kernel_polynomial_from_point(K, algorithm='basic')  # FIXME algorithm='basic' is a workaround for #34907
            step = E.isogeny(poly)
            chain = step * chain
            ker = [step._eval(T) for T in ker]
            E = chain.codomain()

        for iso in E.isomorphisms(self._codomain):
            if self.scaling_factor() == iso.scaling_factor() * chain.scaling_factor():
                break
        else:
            assert False, 'bug in converting fractional isogeny to isogeny chain'

        return iso * chain

    # EllipticCurveHom methods

    @staticmethod
    def _composition_impl(left, right):
        r"""
        Specialized composition method for fractional elliptic-curve morphisms.

        TESTS::

            sage: E = EllipticCurve(GF(419), [-1,0])
            sage: pi = E.frobenius_isogeny()
            sage: endo = (1 + pi) / 2
            sage: 12 * endo
            Composite morphism of degree 15120 = 420*36:
              From: Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419
              To:   Elliptic Curve defined by y^2 = x^3 + 418*x over Finite Field of size 419
        """
        from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
        if isinstance(left, EllipticCurveHom_scalar):
            left, right = right, left
        if isinstance(left, EllipticCurveHom_fractional) and isinstance(right, EllipticCurveHom_scalar):
            r = right._m / left._d
            num, den = r.numerator(), r.denominator()
            if num.is_one():
                f = left._phi
            elif (-num).is_one():
                f = -left._phi
            else:
                f = num * left._phi
            if den.is_one():
                return f
            return EllipticCurveHom_fractional(f, den)
        return NotImplemented

    @staticmethod
    def _comparison_impl(left, right, op):
        r"""
        Compare a fractional elliptic-curve morphism to another elliptic-curve morphism.

        Called by :meth:`EllipticCurveHom._richcmp_`.

        TESTS::

            sage: E = EllipticCurve(GF(419), [-1,0])
            sage: pi = E.frobenius_isogeny()
            sage: endo = (1 + pi) / 2
            sage: 2 * endo - 1 == pi
            True
            sage: 1 - 2 * endo == pi
            False
        """
        assert isinstance(left, EllipticCurveHom_fractional) or isinstance(right, EllipticCurveHom_fractional)
        if isinstance(left, EllipticCurveHom_fractional):
            right = left._d * right
            left = left._phi
        if isinstance(right, EllipticCurveHom_fractional):
            left = right._d * left
            right = right._phi
        return EllipticCurveHom._richcmp_(left, right, op)

    def rational_maps(self):
        r"""
        Return the pair of explicit rational maps defining this fractional isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(GF(419), [1,0])
            sage: phi = E.isogeny(E(185, 73)); phi.rational_maps()
            ((x^5 + 189*x^4 + 9*x^3 + 114*x^2 + 11*x + 206)/(x^4 + 189*x^3 - 105*x^2 - 171*x - 155),
             (x^6*y + 74*x^5*y - 127*x^4*y + 148*x^3*y + 182*x^2*y + 115*x*y + 43*y)/(x^6 + 74*x^5 - 13*x^4 + x^3 - 88*x^2 - 157*x - 179))
            sage: ((phi + phi) / 2).rational_maps()
            ((x^5 + 189*x^4 + 9*x^3 + 114*x^2 + 11*x + 206)/(x^4 + 189*x^3 - 105*x^2 - 171*x - 155),
             (x^6*y + 74*x^5*y - 127*x^4*y + 148*x^3*y + 182*x^2*y + 115*x*y + 43*y)/(x^6 + 74*x^5 - 13*x^4 + x^3 - 88*x^2 - 157*x - 179))
        """
        return self.to_isogeny_chain().rational_maps()

    def x_rational_map(self):
        r"""
        Return the `x`-coordinate rational map of this fractional isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(GF(419), [1,0])
            sage: phi = E.isogeny(E(185, 73)); phi.x_rational_map()
            (x^5 + 189*x^4 + 9*x^3 + 114*x^2 + 11*x + 206)/(x^4 + 189*x^3 + 314*x^2 + 248*x + 264)
            sage: ((phi + phi) / 2).x_rational_map()
            (x^5 + 189*x^4 + 9*x^3 + 114*x^2 + 11*x + 206)/(x^4 + 189*x^3 + 314*x^2 + 248*x + 264)
        """
        return self.to_isogeny_chain().x_rational_map()

    def kernel_polynomial(self):
        r"""
        Return the kernel polynomial of this fractional isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(GF(419), [1,0])
            sage: phi = E.isogeny(E(185, 73)); phi.kernel_polynomial()
            x^2 + 304*x + 39
            sage: ((phi + phi) / 2).kernel_polynomial()
            x^2 + 304*x + 39
        """
        return self.to_isogeny_chain().kernel_polynomial()

    @cached_method
    def dual(self):
        r"""
        Return the dual of this fractional isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(GF(419), [1,0])
            sage: phi = E.isogeny(E(185, 73))
            sage: ((phi + phi) / 2).dual() == phi.dual()
            True
        """
        psi = EllipticCurveHom_fractional(self._phi.dual(), self._d)
        psi.dual.set_cache(self)
        return psi

    def formal(self, *args):
        r"""
        Return the formal isogeny corresponding to this fractional
        isogeny as a power series in the variable `t=-x/y` on the
        domain curve.

        EXAMPLES::

            sage: E = EllipticCurve(GF(419), [-1,0])
            sage: pi = E.frobenius_isogeny()
            sage: ((1 + pi) / 2).formal()
            210*t + 26*t^5 + 254*t^9 + 227*t^13 + 36*t^17 + 74*t^21 + O(t^23)
        """
        return self.to_isogeny_chain().formal(*args)

    def scaling_factor(self):
        r"""
        Return the Weierstrass scaling factor associated to this
        fractional isogeny.

        The scaling factor is the constant `u` (in the base field)
        such that `\varphi^* \omega_2 = u \omega_1`, where
        `\varphi: E_1\to E_2` is this morphism and `\omega_i` are
        the standard Weierstrass differentials on `E_i` defined by
        `\mathrm dx/(2y+a_1x+a_3)`.

        EXAMPLES::

            sage: E = EllipticCurve(GF(419), [-1,0])
            sage: pi = E.frobenius_isogeny()
            sage: ((1 + pi) / 2).scaling_factor()
            210
        """
        # FIXME this can crash when p | d
        return self._phi.scaling_factor() / self._d

    def inseparable_degree(self):
        r"""
        Return the inseparable degree of this morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(419), [-1,0])
            sage: pi = E.frobenius_isogeny()
            sage: ((1 + pi) / 2).inseparable_degree()
            1
            sage: E = EllipticCurve(GF(419), [-1,0])
            sage: pi = E.frobenius_isogeny()
            sage: ((3*pi - pi) / 2).inseparable_degree()
            419
        """
        return self._phi.inseparable_degree() / self._domain.scalar_multiplication(self._d).inseparable_degree()


def _torsion_gens(E, EE, l, e, psi=None):
    if psi is None:
        psi = E.division_polynomial(l**e).radical()
        psi //= psi.gcd(E.division_polynomial(l**(e-1)))

    xs = iter(psi.roots(ring=EE.base_field(), multiplicities=False))
    P = EE.lift_x(next(xs))
    while True:
        Q = EE.lift_x(next(xs))
        if not (P.weil_pairing(Q, l**e)**(l**(e-1))).is_one():
            break
    else:
        assert False, f'bug in finding {l**e}-torsion basis'
    return P, Q
