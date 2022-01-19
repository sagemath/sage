r"""
Composite morphisms of elliptic curves

It is often computationally convenient (for example, in cryptography)
to factor an isogeny of large degree into a composition of isogenies
of smaller (prime) degree. This class implements such a decomposition
while exposing (close to) the same interface as "normal", unfactored
elliptic-curve isogenies.

.. WARNING::

    This module is currently considered experimental.
    It may change in a future release without prior warning, or even
    be removed altogether if things turn out to be unfixably broken.

EXAMPLES:

The following example would take quite literally forever with the
straightforward :class:`EllipticCurveIsogeny` implementation, but
decomposing into prime steps is exponentially faster::

    sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
    doctest:warning
    ...
    sage: p = 3 * 2^143 - 1
    sage: GF(p^2).inject_variables()
    Defining z2
    sage: E = EllipticCurve(GF(p^2), [1,0])
    sage: P = E.lift_x(31415926535897932384626433832795028841971 - z2)
    sage: P.order().factor()
    2^143
    sage: EllipticCurveHom_composite(E, P)
    Composite morphism of degree 11150372599265311570767859136324180752990208 = 2^143:
      From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 33451117797795934712303577408972542258970623^2
      To:   Elliptic Curve defined by y^2 = x^3 + (18676616716352953484576727486205473216172067*z2+32690199585974925193292786311814241821808308)*x + (3369702436351367403910078877591946300201903*z2+15227558615699041241851978605002704626689722) over Finite Field in z2 of size 33451117797795934712303577408972542258970623^2

Yet, the interface provided by :class:`EllipticCurveHom_composite`
is identical to :class:`EllipticCurveIsogeny` and other instantiations
of :class:`EllipticCurveHom`::

    sage: E = EllipticCurve(GF(419), [0,1])
    sage: P = E.lift_x(33); P.order()
    35
    sage: psi = EllipticCurveHom_composite(E, P); psi
    Composite morphism of degree 35 = 5*7:
      From: Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 419
      To:   Elliptic Curve defined by y^2 = x^3 + 101*x + 285 over Finite Field of size 419
    sage: psi(E.lift_x(11))
    (352 : 73 : 1)
    sage: psi.rational_maps()
    ((x^35 + 162*x^34 + 186*x^33 + 92*x^32 - ... + 44*x^3 + 190*x^2 + 80*x - 72)/(x^34 + 162*x^33 - 129*x^32 + 41*x^31 + ... + 66*x^3 - 191*x^2 + 119*x + 21),
     (x^51*y - 176*x^50*y + 115*x^49*y - 120*x^48*y + ... + 72*x^3*y + 129*x^2*y + 163*x*y + 178*y)/(x^51 - 176*x^50 + 11*x^49 + 26*x^48 - ... - 77*x^3 + 185*x^2 + 169*x - 128))
    sage: psi.kernel_polynomial()
    x^17 + 81*x^16 + 7*x^15 + 82*x^14 + 49*x^13 + 68*x^12 + 109*x^11 + 326*x^10 + 117*x^9 + 136*x^8 + 111*x^7 + 292*x^6 + 55*x^5 + 389*x^4 + 175*x^3 + 43*x^2 + 149*x + 373
    sage: psi.dual()
    Composite morphism of degree 35 = 7*5:
      From: Elliptic Curve defined by y^2 = x^3 + 101*x + 285 over Finite Field of size 419
      To:   Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 419
    sage: psi.formal()
    t + 211*t^5 + 417*t^7 + 159*t^9 + 360*t^11 + 259*t^13 + 224*t^15 + 296*t^17 + 139*t^19 + 222*t^21 + O(t^23)

Equality is decided correctly (and, in some cases, much faster than
comparing :meth:`EllipticCurveHom.rational_maps`) even when distinct
factorizations of the same isogeny are compared::

    sage: psi == EllipticCurveIsogeny(E, P)
    True

We can easily obtain the individual factors of the composite map::

    sage: psi.factors()
    (Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 419 to Elliptic Curve defined by y^2 = x^3 + 140*x + 214 over Finite Field of size 419,
     Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 140*x + 214 over Finite Field of size 419 to Elliptic Curve defined by y^2 = x^3 + 101*x + 285 over Finite Field of size 419)

AUTHORS:

- Lukas Zobernig (2020): initial proof-of-concept version
- Lorenz Panny (2021): :class:`EllipticCurveHom` interface,
  documentation and tests, equality testing
"""

from sage.structure.richcmp import op_EQ, op_NE
from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

from sage.arith.misc import prod, integer_floor
from sage.functions.log import log
from sage.functions.other import sqrt

from sage.rings.finite_rings import finite_field_base
from sage.rings.number_field import number_field_base

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.hom import EllipticCurveHom
from sage.schemes.elliptic_curves.ell_curve_isogeny import EllipticCurveIsogeny
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

from sage.misc.superseded import experimental_warning
experimental_warning(32744, 'EllipticCurveHom_composite is experimental code.')

#TODO: implement sparse strategies? (cf. the SIKE cryptosystem)

def _eval_factored_isogeny(phis, P):
    """
    This method pushes a point `P` through a given sequence ``phis``
    of compatible isogenies.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves import hom_composite
        sage: E = EllipticCurve(GF(419), [1,0])
        sage: Q = E(21,8)
        sage: phis = []
        sage: while len(phis) < 10:
        ....:     P = list(sorted(E(0).division_points(7)))[1]
        ....:     phis.append(E.isogeny(P))
        ....:     E = phis[-1].codomain()
        sage: R = hom_composite._eval_factored_isogeny(phis, Q); R
        (290 : 183 : 1)
        sage: R in E
        True
    """
    for phi in phis:
        P = phi(P)
    return P

def _compute_factored_isogeny_prime_power(P, l, e):
    """
    This method takes a point `P` of order `l^e` and returns
    a sequence of degree-`l` isogenies whose composition has
    the subgroup generated by `P` as its kernel.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves import hom_composite
        sage: E = EllipticCurve(GF(8191), [1,0])
        sage: P = E.random_point()
        sage: (l,e), = P.order().factor()
        sage: phis = hom_composite._compute_factored_isogeny_prime_power(P,l,e)
        sage: hom_composite._eval_factored_isogeny(phis, P)
        (0 : 1 : 0)
        sage: [phi.degree() for phi in phis] == [l]*e
        True
    """
    E = P.curve()
    phis = []
    for i in range(e):
        K = l**(e-1-i) * P
        phi = EllipticCurveIsogeny(E, K)
        E = phi.codomain()
        P = phi(P)
        phis.append(phi)
    return phis

def _compute_factored_isogeny_single_generator(P):
    """
    This method takes a point `P` and returns a sequence of
    prime-degree isogenies whose composition has the subgroup
    generated by `P` as its kernel.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves import hom_composite
        sage: E = EllipticCurve(GF(419), [1,0])
        sage: P = E(42,321)
        sage: phis = hom_composite._compute_factored_isogeny_single_generator(P)
        sage: list(sorted(phi.degree() for phi in phis))
        [2, 2, 3, 5, 7]
        sage: hom_composite._eval_factored_isogeny(phis, P)
        (0 : 1 : 0)
    """
    phis = []
    h = P.order()
    for l,e in P.order().factor():
        h //= l**e
        psis = _compute_factored_isogeny_prime_power(h*P, l, e)
        P = _eval_factored_isogeny(psis, P)
        phis += psis
    return phis

def _compute_factored_isogeny(kernel):
    """
    This method takes a set of points on an elliptic curve
    and returns a sequence of isogenies whose composition
    has the subgroup generated by that subset as its kernel.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves import hom_composite
        sage: E = EllipticCurve(GF(419), [-1,0])
        sage: Ps = [E(41,99), E(41,-99), E(51,14), E(21,21), E(33,17)]
        sage: phis = hom_composite._compute_factored_isogeny(Ps)
        sage: [phi.degree() for phi in phis]
        [2, 3, 5, 7, 2]
        sage: {hom_composite._eval_factored_isogeny(phis, P) for P in Ps}
        {(0 : 1 : 0)}
    """
    phis = []
    ker = list(kernel)
    while ker:
        K, ker = ker[0], ker[1:]
        psis = _compute_factored_isogeny_single_generator(K)
        ker = [_eval_factored_isogeny(psis, P) for P in ker]
        phis += psis
    return phis

class EllipticCurveHom_composite(EllipticCurveHom):

    _degree = None
    _phis = None

    def __init__(self, E, kernel, codomain=None):
        """
        Construct a composite isogeny with given kernel (and optionally,
        prescribed codomain curve). The isogeny is decomposed into steps
        of prime degree.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(419), [1,0])
            sage: EllipticCurveHom_composite(E, E.lift_x(23))
            Composite morphism of degree 105 = 3*5*7:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 419
              To:   Elliptic Curve defined by y^2 = x^3 + 373*x + 126 over Finite Field of size 419


        The given kernel generators need not be independent::

            sage: K.<a> = NumberField(x^2 - x - 5)
            sage: E = EllipticCurve('210.b6').change_ring(K)
            sage: E.torsion_subgroup()
            Torsion Subgroup isomorphic to Z/12 + Z/2 associated to the Elliptic Curve defined by y^2 + x*y + y = x^3 + (-578)*x + 2756 over Number Field in a with defining polynomial x^2 - x - 5
            sage: EllipticCurveHom_composite(E, E.torsion_points())
            Composite morphism of degree 24 = 2^3*3:
              From: Elliptic Curve defined by y^2 + x*y + y = x^3 + (-578)*x + 2756 over Number Field in a with defining polynomial x^2 - x - 5
              To:   Elliptic Curve defined by y^2 + x*y + y = x^3 + (-89915533/16)*x + (-328200928141/64) over Number Field in a with defining polynomial x^2 - x - 5

        TESTS::

            sage: E = EllipticCurve(GF(19), [1,0])
            sage: P = E.random_point()
            sage: psi = EllipticCurveHom_composite(E, P)
            sage: psi   # random
            Composite morphism of degree 10 = 2*5:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19
              To:   Elliptic Curve defined by y^2 = x^3 + 14*x over Finite Field of size 19

        ::

            sage: EllipticCurveHom_composite(E, E.lift_x(3), codomain=E)
            Composite morphism of degree 20 = 2^2*5:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19
        """
        if not isinstance(E, EllipticCurve_generic):
            raise ValueError(f'not an elliptic curve: {E}')

        if not isinstance(kernel, list) and not isinstance(kernel, tuple):
            kernel = [kernel]

        for P in kernel:
            if not P in E:
                raise ValueError(f'given point {P} does not lie on {E}')

        self._phis = _compute_factored_isogeny(kernel)

        if not self._phis:
            self._phis = [WeierstrassIsomorphism(E, (1,0,0,0))]

        if codomain is not None:
            if not isinstance(codomain, EllipticCurve_generic):
                raise ValueError(f'not an elliptic curve: {codomain}')
            iso = self._phis[-1].codomain().isomorphism_to(codomain)
            if hasattr(self._phis[-1], '_set_post_isomorphism'):
                self._phis[-1]._set_post_isomorphism(iso)
            else:
                self._phis.append(iso)

        self._phis = tuple(self._phis)  # make immutable
        self.__perform_inheritance_housekeeping()

    def __perform_inheritance_housekeeping(self):
        """
        Internal helper function to set values of the superclass.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve([1,0])
            sage: phi = EllipticCurveHom_composite(E, E(0,0))   # implicit doctest
            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: print(EllipticCurveHom._repr_(phi))
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 = x^3 + x over Rational Field
              To:   Elliptic Curve defined by y^2 = x^3 - 4*x over Rational Field
            sage: phi.domain()
            Elliptic Curve defined by y^2 = x^3 + x over Rational Field
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 - 4*x over Rational Field
        """
        self._degree = prod(phi.degree() for phi in self._phis)
        self._domain = self._phis[0].domain()
        self._codomain = self._phis[-1].codomain()
        EllipticCurveHom.__init__(self, self._domain, self._codomain)

    @classmethod
    def from_factors(cls, maps, E=None, strict=True):
        r"""
        This method constructs a :class:`EllipticCurveHom_composite`
        object encapsulating a given sequence of compatible isogenies.

        The isogenies are composed in left-to-right order, i.e., the
        resulting composite map equals `f_{n-1} \circ \dots \circ f_0`
        where `f_i` denotes ``maps[i]``.

        INPUT:

        - ``maps`` -- sequence of :class:`EllipticCurveHom` objects
        - ``E`` (optional) -- the domain elliptic curve
        - ``strict`` (optional, default: ``True``) -- if ``True``,
          always return an :class:`EllipticCurveHom_composite` object;
          else may return another :class:`EllipticCurveHom` type

        OUTPUT: the composite of ``maps``

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(43), [1,0])
            sage: P, = E.gens()
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: psi = EllipticCurveHom_composite.from_factors(phi.factors())
            sage: psi == phi
            True

        TESTS::

            sage: E = EllipticCurve('4730k1')
            sage: EllipticCurveHom_composite.from_factors([], E) == E.multiplication_by_m_isogeny(1)
            True

        ::

            sage: E = EllipticCurve(GF(419), [1,0])
            sage: P, = E.gens()
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: EllipticCurveHom_composite.from_factors(phi.factors()) == phi
            True
        """
        maps = tuple(maps)
        if not maps and E is None:
            raise ValueError('need either factors or domain')
        if E is None:
            E = maps[0].domain()

        for phi in maps:
            if not isinstance(phi, EllipticCurveHom):
                raise TypeError(f'not an elliptic-curve isogeny: {phi}')
            if phi.domain() != E:
                raise ValueError(f'isogeny has incorrect domain: {phi}')
            E = phi.codomain()

        if not maps:
            maps = (WeierstrassIsomorphism(E, (1,0,0,0)),)

        if len(maps) == 1 and not strict:
            return maps[0]

        result = cls.__new__(cls)
        result._phis = maps
        result.__perform_inheritance_housekeeping()
        return result

    def _call_(self, P):
        """
        Evaluate this composite isogeny at a point.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: K.<a> = NumberField(x^2 - x - 5)
            sage: E = EllipticCurve('210.b6').change_ring(K)
            sage: psi = EllipticCurveHom_composite(E, E.torsion_points())
            sage: R = E.lift_x(15/4 * (a+3))
            sage: psi(R)    # indirect doctest
            (1033648757/303450 : 58397496786187/1083316500*a - 62088706165177/2166633000 : 1)
        """
        return _eval_factored_isogeny(self._phis, P)

    def _eval(self, P):
        r"""
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self'' at the given point

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(j=Mod(1728,419))
            sage: K, = E.gens()
            sage: psi = EllipticCurveHom_composite(E, 4*K)
            sage: Ps = E.change_ring(GF(419**2))(0).division_points(5)
            sage: {psi._eval(P).curve() for P in Ps}
            {Elliptic Curve defined by y^2 = x^3 + 373*x + 126 over Finite Field in z2 of size 419^2}
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')
        Q = P
        for phi in self._phis:
            Q = phi._eval(Q)
        k = Sequence(tuple(P) + tuple(Q)).universe()
        return self._codomain.base_extend(k)(*Q)

    def _repr_(self):
        """
        Return basic facts about this composite isogeny as a string.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(43), [1,0])
            sage: P, = E.gens()
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: phi   # indirect doctest
            Composite morphism of degree 44 = 2^2*11:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43
            sage: phi * phi * phi * phi * phi * phi * phi   # indirect doctest
            Composite morphism of degree 319277809664 = 2^2*11*2^2*11*2^2*11*2^2*11*2^2*11*2^2*11*2^2*11:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43
        """
        from itertools import groupby
        degs = [phi.degree() for phi in self._phis]
        grouped = [(d, sum(1 for _ in g)) for d,g in groupby(degs)]
        degs_str = '*'.join(str(d) + (f'^{e}' if e>1 else '') for d,e in grouped)
        return f'Composite morphism of degree {self._degree} = {degs_str}:' \
                f'\n  From: {self._domain}' \
                f'\n  To:   {self._codomain}'

    def factors(self):
        r"""
        Return the factors of this composite isogeny as a tuple.

        The isogenies are returned in left-to-right order, i.e.,
        the returned tuple `(f_1,...,f_n)` corresponds to the map
        `f_n \circ \dots \circ f_1`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(43), [1,0])
            sage: P, = E.gens()
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: phi.factors()
            (Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43 to Elliptic Curve defined by y^2 = x^3 + 39*x over Finite Field of size 43,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 39*x over Finite Field of size 43 to Elliptic Curve defined by y^2 = x^3 + 42*x + 26 over Finite Field of size 43,
             Isogeny of degree 11 from Elliptic Curve defined by y^2 = x^3 + 42*x + 26 over Finite Field of size 43 to Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43)
        """
        return self._phis


    @staticmethod
    def _composition_impl(left, right):
        """
        Helper method to compose other elliptic-curve morphisms with
        :class:`EllipticCurveHom_composite` objects. Called by
        :meth:`EllipticCurveHom._composition_`.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve([i+1, i, 0, -4, -6*i])
            sage: P,Q = E.lift_x(i-5), E.lift_x(-4*i)
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: psi = phi.codomain().isogeny(phi(Q))
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: iso1 = WeierstrassIsomorphism(E, (-1, 0, -i-1, 0))
            sage: iso2 = psi.codomain().isomorphism_to(E)
            sage: psi * phi                 # indirect doctest
            Composite morphism of degree 16 = 2^2*4:
              From: Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I) over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-3331/4)*x + (-142593/8*I) over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: iso2 * EllipticCurveHom_composite.from_factors([phi, psi]) # indirect doctest
            Composite morphism of degree 16 = 4^2:
              From: Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I) over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I) over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: phi * iso1                # indirect doctest
            Composite morphism of degree 4 = 2^2:
              From: Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I) over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (480*I-694)*x + (-7778*I+5556) over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: iso2 * psi * phi * iso1   # indirect doctest
            Composite morphism of degree 16 = 2^2*4:
              From: Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I) over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I) over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        if isinstance(left, EllipticCurveHom_composite):
            if isinstance(right, WeierstrassIsomorphism) and hasattr(left.factors()[0], '_set_pre_isomorphism'):    #XXX bit of a hack
                return EllipticCurveHom_composite.from_factors((left.factors()[0] * right,) + left.factors()[1:], strict=False)
            if isinstance(right, EllipticCurveHom_composite):
                return EllipticCurveHom_composite.from_factors(right.factors() + left.factors())
            if isinstance(right, EllipticCurveHom):
                return EllipticCurveHom_composite.from_factors((right,) + left.factors())
        if isinstance(right, EllipticCurveHom_composite):
            if isinstance(left, WeierstrassIsomorphism) and hasattr(right.factors()[-1], '_set_post_isomorphism'):  #XXX bit of a hack
                return EllipticCurveHom_composite.from_factors(right.factors()[:-1] + (left * right.factors()[-1],), strict=False)
            if isinstance(left, EllipticCurveHom):
                return EllipticCurveHom_composite.from_factors(right.factors() + (left,))
        return NotImplemented


    # EllipticCurveHom methods

    def degree(self):
        """
        Return the degree of this composite isogeny.

        Degrees are multiplicative, so this is the product of the
        degrees of the individual factors.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(419), [1,0])
            sage: P, = E.gens()
            sage: phi = EllipticCurveHom_composite(E, P+P)
            sage: phi.degree()
            210

        TESTS::

            sage: phi.degree() == prod(f.degree() for f in phi.factors())
            True
        """
        return self._degree

    def _richcmp_(self, other, op):
        r"""
        Compare this composite isogeny to another elliptic-curve morphism.

        ALGORITHM:

        Over finite fields and number fields, we use the fact that two
        isogenies of equal degree `d` are the same if and only if they
        act identically on more than `4d` points. (It suffices to check
        this on a few points that generate a large enough subgroup.)

        If the domain curve does not have sufficiently many rational
        points, the base field is extended first. Since an extension of
        degree `O(\log(d))` suffices, the complexity is polynomial in
        the representation size of this morphism.

        Over more general base fields, we fall back to comparing the
        results of :meth:`rational_maps`, which takes time at least
        linear in the degree.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(QuadraticField(-3), [0,16])
            sage: P,Q = E.lift_x(0), E.lift_x(-4)
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: psi = phi.codomain().isogeny(phi(Q))
            sage: psi = psi.codomain().isomorphism_to(E) * psi
            sage: comp = psi * phi
            sage: mu = E.multiplication_by_m_isogeny(phi.degree())
            sage: sum(a*comp == mu for a in E.automorphisms())
            1

        ::

            sage: E = EllipticCurve(GF(431**2), [1,0])
            sage: P,Q = E.gens()
            sage: phi1 = EllipticCurveHom_composite(E, P)
            sage: phi2 = EllipticCurveHom_composite(phi1.codomain(), phi1(Q))
            sage: psi1 = EllipticCurveHom_composite(E, Q)
            sage: psi2 = EllipticCurveHom_composite(psi1.codomain(), psi1(P))
            sage: phi2 * phi1 == psi2 * psi1
            True
        """
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        if op != op_EQ:
            return NotImplemented

        if self.domain() != other.domain():
            return False
        if self.codomain() != other.codomain():
            return False
        if self.degree() != other.degree():
            return False

        E = self.domain()
        F = E.base_ring()

        if isinstance(F, finite_field_base.FiniteField):
            q = F.cardinality()
            d = self.degree()
            e = integer_floor(1 + 2 * log(2*sqrt(d) + 1, q))    # from Hasse bound
            e = next(i for i,n in enumerate(E.count_points(e+1), 1) if n > 4*d)
            EE = E.base_extend(F.extension(e))
            Ps = EE.gens()
            return all(self._eval(P) == other._eval(P) for P in Ps)

        elif isinstance(F, number_field_base.NumberField):
            for _ in range(100):
                P = E.lift_x(F.random_element(), extend=True)
                if not P.has_finite_order():
                    return self._eval(P) == other._eval(P)
            else:
                assert False, "_richcmp_() couldn't find a point of infinite order"

        else:
            # fall back to generic method
            return self.rational_maps() == other.rational_maps()

    def rational_maps(self):
        """
        Return the pair of explicit rational maps defining this composite
        isogeny.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])
            sage: P = E.lift_x(7321)
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: phi.rational_maps()
            ((x^9 + 27463*x^8 + 21204*x^7 - 5750*x^6 + 1610*x^5 + 14440*x^4 + 26605*x^3 - 15569*x^2 - 3341*x + 1267)/(x^8 + 27463*x^7 + 26871*x^6 + 5999*x^5 - 20194*x^4 - 6310*x^3 + 24366*x^2 - 20905*x - 13867),
             (x^12*y + 8426*x^11*y + 5667*x^11 + 27612*x^10*y + 26124*x^10 + 9688*x^9*y - 22715*x^9 + 19864*x^8*y + 498*x^8 + 22466*x^7*y - 14036*x^7 + 8070*x^6*y + 19955*x^6 - 20765*x^5*y - 12481*x^5 + 12672*x^4*y + 24142*x^4 - 23695*x^3*y + 26667*x^3 + 23780*x^2*y + 17864*x^2 + 15053*x*y - 30118*x + 17539*y - 23609)/(x^12 + 8426*x^11 + 21945*x^10 - 22587*x^9 + 22094*x^8 + 14603*x^7 - 26255*x^6 + 11171*x^5 - 16508*x^4 - 14435*x^3 - 2170*x^2 + 29081*x - 19009))

        TESTS::

            sage: f = phi.codomain().defining_polynomial()
            sage: g = E.defining_polynomial().subs({2:1})
            sage: f(*phi.rational_maps(), 1) % g
            0

        ::

            sage: phi.rational_maps()[0].parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Finite Field of size 65537
            sage: phi.rational_maps()[1].parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Finite Field of size 65537
        """
        fx, fy = self._phis[-1].rational_maps()
        for phi in self._phis[:-1][::-1]:
            gx, gy = phi.rational_maps()
            fx, fy = fx(gx, gy), fy(gx, gy)
        return (fx, fy)

    def x_rational_map(self):
        """
        Return the `x`-coordinate rational map of this composite isogeny.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])
            sage: P = E.lift_x(7321)
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: phi.x_rational_map() == phi.rational_maps()[0]
            True

        TESTS::

            sage: phi.x_rational_map().parent()
            Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 65537
        """
        fx = self._phis[-1].x_rational_map()
        for phi in self._phis[:-1][::-1]:
            fx = fx(phi.x_rational_map())
        return fx

    def kernel_polynomial(self):
        """
        Return the kernel polynomial of this composite isogeny.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])
            sage: P = E.lift_x(7321)
            sage: phi = EllipticCurveHom_composite(E, P); phi
            Composite morphism of degree 9 = 3^2:
              From: Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 65537
              To:   Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 28339*x + 59518 over Finite Field of size 65537
            sage: phi.kernel_polynomial()
            x^4 + 46500*x^3 + 19556*x^2 + 7643*x + 15952
        """
        # shouldn't there be a better algorithm for this?
        return self.x_rational_map().denominator().radical()

    @cached_method
    def dual(self):
        """
        Return the dual of this composite isogeny.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])
            sage: P = E.lift_x(7321)
            sage: phi = EllipticCurveHom_composite(E, P); phi
            Composite morphism of degree 9 = 3^2:
              From: Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 65537
              To:   Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 28339*x + 59518 over Finite Field of size 65537
            sage: psi = phi.dual(); psi
            Composite morphism of degree 9 = 3^2:
              From: Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 28339*x + 59518 over Finite Field of size 65537
              To:   Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 65537
            sage: psi * phi == phi.domain().multiplication_by_m_isogeny(phi.degree())
            True
            sage: phi * psi == psi.domain().multiplication_by_m_isogeny(psi.degree())
            True
        """
        phis = (phi.dual() for phi in self._phis[::-1])
        return EllipticCurveHom_composite.from_factors(phis)

    def is_separable(self):
        """
        Determine whether this composite isogeny is separable.

        A composition of isogenies is separable if and only if
        all factors are.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(7^2), [3,2])
            sage: P = E.lift_x(1)
            sage: phi = EllipticCurveHom_composite(E, P); phi
            Composite morphism of degree 7 = 7:
              From: Elliptic Curve defined by y^2 = x^3 + 3*x + 2 over Finite Field in z2 of size 7^2
              To:   Elliptic Curve defined by y^2 = x^3 + 3*x + 2 over Finite Field in z2 of size 7^2
            sage: phi.is_separable()
            True
        """
        return all(phi.is_separable() for phi in self._phis)

    def formal(self, prec=20):
        """
        Return the formal isogeny corresponding to this composite
        isogeny as a power series in the variable `t=-x/y` on the
        domain curve.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])
            sage: P = E.lift_x(7321)
            sage: phi = EllipticCurveHom_composite(E, P)
            sage: phi.formal()
            t + 54203*t^5 + 48536*t^6 + 40698*t^7 + 37808*t^8 + 21111*t^9 + 42381*t^10 + 46688*t^11 + 657*t^12 + 38916*t^13 + 62261*t^14 + 59707*t^15 + 30767*t^16 + 7248*t^17 + 60287*t^18 + 50451*t^19 + 38305*t^20 + 12312*t^21 + 31329*t^22 + O(t^23)
            sage: (phi.dual() * phi).formal(prec=5)
            9*t + 65501*t^2 + 65141*t^3 + 59183*t^4 + 21491*t^5 + 8957*t^6 + 999*t^7 + O(t^8)
        """
        res = self._phis[-1].formal(prec=prec)
        for phi in self._phis[:-1][::-1]:
            res = res(phi.formal(prec=prec))
        return res

    def is_injective(self):
        """
        Determine whether this composite morphism has trivial kernel.

        In other words, return ``True`` if and only if ``self`` is a
        purely inseparable isogeny.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve([1,0])
            sage: phi = EllipticCurveHom_composite(E, E(0,0))
            sage: phi.is_injective()
            False
            sage: E = EllipticCurve_from_j(GF(3).algebraic_closure()(0))
            sage: nu = EllipticCurveHom_composite.from_factors(E.automorphisms())
            sage: nu
            Composite morphism of degree 1 = 1^12:
              From: Elliptic Curve defined by y^2 = x^3 + x over Algebraic closure of Finite Field of size 3
              To:   Elliptic Curve defined by y^2 = x^3 + x over Algebraic closure of Finite Field of size 3
            sage: nu.is_injective()
            True
        """
        return all(phi.is_injective() for phi in self._phis)


    @staticmethod
    def make_default():
        r"""
        Calling this method will override the composition method
        of :class:`EllipticCurveHom` such that it constructs a
        :class:`EllipticCurveHom_composite` object by default,
        rather than a :class:`sage.categories.map.FormalCompositeMap`.

        .. WARNING::

            This method exists only temporarily to make testing more
            convenient while :class:`EllipticCurveHom_composite` is
            experimental.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(587), [1,0])
            sage: P = E(3,404)
            sage: phi = E.isogeny(7*P)
            sage: psi = phi.codomain().isogeny(phi(P))
            sage: psi * phi
            Composite map:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 587
              To:   Elliptic Curve defined by y^2 = x^3 + 296*x + 164 over Finite Field of size 587
              Defn:   Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 587 to Elliptic Curve defined by y^2 = x^3 + 126*x + 500 over Finite Field of size 587
                    then
                      Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 126*x + 500 over Finite Field of size 587 to Elliptic Curve defined by y^2 = x^3 + 296*x + 164 over Finite Field of size 587
            sage: EllipticCurveHom_composite.make_default()
            sage: psi * phi
            Composite morphism of degree 49 = 7^2:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 587
              To:   Elliptic Curve defined by y^2 = x^3 + 296*x + 164 over Finite Field of size 587
            sage: (psi * phi).factors()
            (Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 587 to Elliptic Curve defined by y^2 = x^3 + 126*x + 500 over Finite Field of size 587,
             Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 126*x + 500 over Finite Field of size 587 to Elliptic Curve defined by y^2 = x^3 + 296*x + 164 over Finite Field of size 587)
        """
        def _composition_(self, other, homset):
            if not isinstance(self, EllipticCurveHom) or not isinstance(other, EllipticCurveHom):
                raise TypeError(f'cannot compose {type(self)} with {type(other)}')
            ret = self._composition_impl(self, other)
            if ret is not NotImplemented:
                return ret
            ret = other._composition_impl(self, other)
            if ret is not NotImplemented:
                return ret
            return EllipticCurveHom_composite.from_factors([other, self])
        EllipticCurveHom._composition_ = _composition_

