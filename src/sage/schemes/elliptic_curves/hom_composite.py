r"""
Composite morphisms of elliptic curves

It is often computationally convenient (for example, in cryptography)
to factor an isogeny of large degree into a composition of isogenies
of smaller (prime) degree. This class implements such a decomposition
while exposing (close to) the same interface as "normal", unfactored
elliptic-curve isogenies.

EXAMPLES:

The following example would take quite literally forever with the
straightforward :class:`EllipticCurveIsogeny` implementation, but
decomposing into prime steps is exponentially faster::

    sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
    sage: p = 3 * 2^143 - 1
    sage: GF(p^2).inject_variables()                                                    # optional - sage.rings.finite_rings
    Defining z2
    sage: E = EllipticCurve(GF(p^2), [1,0])                                             # optional - sage.rings.finite_rings
    sage: P = E.lift_x(31415926535897932384626433832795028841971 - z2)                  # optional - sage.rings.finite_rings
    sage: P.order().factor()                                                            # optional - sage.rings.finite_rings
    2^143
    sage: EllipticCurveHom_composite(E, P)                                              # optional - sage.rings.finite_rings
    Composite morphism of degree 11150372599265311570767859136324180752990208 = 2^143:
      From: Elliptic Curve defined by y^2 = x^3 + x
            over Finite Field in z2 of size 33451117797795934712303577408972542258970623^2
      To:   Elliptic Curve defined by y^2 = x^3 + (18676616716352953484576727486205473216172067*z2+32690199585974925193292786311814241821808308)*x
            + (3369702436351367403910078877591946300201903*z2+15227558615699041241851978605002704626689722)
            over Finite Field in z2 of size 33451117797795934712303577408972542258970623^2

Yet, the interface provided by :class:`EllipticCurveHom_composite`
is identical to :class:`EllipticCurveIsogeny` and other instantiations
of :class:`EllipticCurveHom`::

    sage: E = EllipticCurve(GF(419), [0,1])                                             # optional - sage.rings.finite_rings
    sage: P = E.lift_x(33); P.order()                                                   # optional - sage.rings.finite_rings
    35
    sage: psi = EllipticCurveHom_composite(E, P); psi                                   # optional - sage.rings.finite_rings
    Composite morphism of degree 35 = 5*7:
      From: Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 419
      To:   Elliptic Curve defined by y^2 = x^3 + 101*x + 285 over Finite Field of size 419
    sage: psi(E.lift_x(11))                                                             # optional - sage.rings.finite_rings
    (352 : 73 : 1)
    sage: psi.rational_maps()                                                           # optional - sage.rings.finite_rings
    ((x^35 + 162*x^34 + 186*x^33 + 92*x^32 - ... + 44*x^3 + 190*x^2 + 80*x
      - 72)/(x^34 + 162*x^33 - 129*x^32 + 41*x^31 + ... + 66*x^3 - 191*x^2 + 119*x + 21),
     (x^51*y - 176*x^50*y + 115*x^49*y - 120*x^48*y + ... + 72*x^3*y + 129*x^2*y + 163*x*y
      + 178*y)/(x^51 - 176*x^50 + 11*x^49 + 26*x^48 - ... - 77*x^3 + 185*x^2 + 169*x - 128))
    sage: psi.kernel_polynomial()                                                       # optional - sage.rings.finite_rings
    x^17 + 81*x^16 + 7*x^15 + 82*x^14 + 49*x^13 + 68*x^12 + 109*x^11 + 326*x^10
     + 117*x^9 + 136*x^8 + 111*x^7 + 292*x^6 + 55*x^5 + 389*x^4 + 175*x^3 + 43*x^2 + 149*x + 373
    sage: psi.dual()                                                                    # optional - sage.rings.finite_rings
    Composite morphism of degree 35 = 7*5:
      From: Elliptic Curve defined by y^2 = x^3 + 101*x + 285 over Finite Field of size 419
      To:   Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 419
    sage: psi.formal()                                                                  # optional - sage.rings.finite_rings
    t + 211*t^5 + 417*t^7 + 159*t^9 + 360*t^11 + 259*t^13 + 224*t^15 + 296*t^17 + 139*t^19 + 222*t^21 + O(t^23)

Equality is decided correctly (and, in some cases, much faster than
comparing :meth:`EllipticCurveHom.rational_maps`) even when distinct
factorizations of the same isogeny are compared::

    sage: psi == EllipticCurveIsogeny(E, P)                                             # optional - sage.rings.finite_rings
    True

We can easily obtain the individual factors of the composite map::

    sage: psi.factors()                                                                 # optional - sage.rings.finite_rings
    (Isogeny of degree 5
      from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 419
        to Elliptic Curve defined by y^2 = x^3 + 140*x + 214 over Finite Field of size 419,
     Isogeny of degree 7
      from Elliptic Curve defined by y^2 = x^3 + 140*x + 214 over Finite Field of size 419
        to Elliptic Curve defined by y^2 = x^3 + 101*x + 285 over Finite Field of size 419)

AUTHORS:

- Lukas Zobernig (2020): initial proof-of-concept version
- Lorenz Panny (2021): :class:`EllipticCurveHom` interface,
  documentation and tests, equality testing
"""

from sage.structure.richcmp import op_EQ
from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

from sage.arith.misc import prod

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.hom import EllipticCurveHom, compare_via_evaluation
from sage.schemes.elliptic_curves.ell_curve_isogeny import EllipticCurveIsogeny
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism, identity_morphism

# TODO: Implement sparse strategies? (cf. the SIKE cryptosystem)


def _eval_factored_isogeny(phis, P):
    """
    This method pushes a point `P` through a given sequence ``phis``
    of compatible isogenies.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves import hom_composite
        sage: E = EllipticCurve(GF(419), [1,0])                                         # optional - sage.rings.finite_rings
        sage: Q = E(21, 8)                                                              # optional - sage.rings.finite_rings
        sage: phis = []                                                                 # optional - sage.rings.finite_rings
        sage: while len(phis) < 10:                                                     # optional - sage.rings.finite_rings
        ....:     P = list(sorted(E(0).division_points(7)))[1]
        ....:     phis.append(E.isogeny(P))
        ....:     E = phis[-1].codomain()
        sage: R = hom_composite._eval_factored_isogeny(phis, Q); R                      # optional - sage.rings.finite_rings
        (290 : 183 : 1)
        sage: R in E                                                                    # optional - sage.rings.finite_rings
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
        sage: E = EllipticCurve(GF(8191), [1,0])                                        # optional - sage.rings.finite_rings
        sage: P = E.random_point()                                                      # optional - sage.rings.finite_rings
        sage: (l,e), = P.order().factor()                                               # optional - sage.rings.finite_rings
        sage: phis = hom_composite._compute_factored_isogeny_prime_power(P, l, e)       # optional - sage.rings.finite_rings
        sage: hom_composite._eval_factored_isogeny(phis, P)                             # optional - sage.rings.finite_rings
        (0 : 1 : 0)
        sage: [phi.degree() for phi in phis] == [l]*e                                   # optional - sage.rings.finite_rings
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
        sage: E = EllipticCurve(GF(419), [1,0])                                         # optional - sage.rings.finite_rings
        sage: P = E(42, 321)                                                            # optional - sage.rings.finite_rings
        sage: phis = hom_composite._compute_factored_isogeny_single_generator(P)        # optional - sage.rings.finite_rings
        sage: list(sorted(phi.degree() for phi in phis))                                # optional - sage.rings.finite_rings
        [2, 2, 3, 5, 7]
        sage: hom_composite._eval_factored_isogeny(phis, P)                             # optional - sage.rings.finite_rings
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
        sage: E = EllipticCurve(GF(419), [-1,0])                                        # optional - sage.rings.finite_rings
        sage: Ps = [E(41,99), E(41,-99), E(51,14), E(21,21), E(33,17)]                  # optional - sage.rings.finite_rings
        sage: phis = hom_composite._compute_factored_isogeny(Ps)                        # optional - sage.rings.finite_rings
        sage: [phi.degree() for phi in phis]                                            # optional - sage.rings.finite_rings
        [2, 3, 5, 7, 2]
        sage: {hom_composite._eval_factored_isogeny(phis, P) for P in Ps}               # optional - sage.rings.finite_rings
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

    def __init__(self, E, kernel, codomain=None, model=None):
        """
        Construct a composite isogeny with given kernel (and optionally,
        prescribed codomain curve). The isogeny is decomposed into steps
        of prime degree.

        The ``codomain`` and ``model`` parameters have the same meaning
        as for :class:`EllipticCurveIsogeny`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(419), [1,0])                                     # optional - sage.rings.finite_rings
            sage: EllipticCurveHom_composite(E, E.lift_x(23))                           # optional - sage.rings.finite_rings
            Composite morphism of degree 105 = 3*5*7:
              From: Elliptic Curve defined by y^2 = x^3 + x
                    over Finite Field of size 419
              To:   Elliptic Curve defined by y^2 = x^3 + 373*x + 126
                    over Finite Field of size 419

        The given kernel generators need not be independent::

            sage: K.<a> = NumberField(x^2 - x - 5)                                      # optional - sage.rings.number_field
            sage: E = EllipticCurve('210.b6').change_ring(K)                            # optional - sage.rings.number_field
            sage: E.torsion_subgroup()                                                  # optional - sage.rings.number_field
            Torsion Subgroup isomorphic to Z/12 + Z/2 associated to the Elliptic Curve
             defined by y^2 + x*y + y = x^3 + (-578)*x + 2756
              over Number Field in a with defining polynomial x^2 - x - 5
            sage: EllipticCurveHom_composite(E, E.torsion_points())                     # optional - sage.rings.number_field
            Composite morphism of degree 24 = 2^3*3:
              From: Elliptic Curve defined by y^2 + x*y + y = x^3 + (-578)*x + 2756
                    over Number Field in a with defining polynomial x^2 - x - 5
              To:   Elliptic Curve defined by
                    y^2 + x*y + y = x^3 + (-89915533/16)*x + (-328200928141/64)
                    over Number Field in a with defining polynomial x^2 - x - 5

        TESTS::

            sage: E = EllipticCurve(GF(19), [1,0])                                      # optional - sage.rings.finite_rings
            sage: P = E.random_point()                                                  # optional - sage.rings.finite_rings
            sage: psi = EllipticCurveHom_composite(E, P)                                # optional - sage.rings.finite_rings
            sage: psi   # random                                                        # optional - sage.rings.finite_rings
            Composite morphism of degree 10 = 2*5:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19
              To:   Elliptic Curve defined by y^2 = x^3 + 14*x over Finite Field of size 19

        ::

            sage: EllipticCurveHom_composite(E, E.lift_x(3), codomain=E)                # optional - sage.rings.finite_rings
            Composite morphism of degree 20 = 2^2*5:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19

        ::

            sage: E = EllipticCurve(GF((2^127-1)^2), [1,0])                             # optional - sage.rings.finite_rings
            sage: K = 2^30 * E.random_point()                                           # optional - sage.rings.finite_rings
            sage: psi = EllipticCurveHom_composite(E, K, model='montgomery')            # optional - sage.rings.finite_rings
            sage: psi.codomain().a_invariants()                                         # optional - sage.rings.finite_rings
            (0, ..., 0, 1, 0)
        """
        if not isinstance(E, EllipticCurve_generic):
            raise ValueError(f'not an elliptic curve: {E}')

        if not isinstance(kernel, list) and not isinstance(kernel, tuple):
            kernel = [kernel]

        for P in kernel:
            if P not in E:
                raise ValueError(f'given point {P} does not lie on {E}')

        self._phis = _compute_factored_isogeny(kernel)

        if not self._phis:
            self._phis = [identity_morphism(E)]

        if model is not None:
            if codomain is not None:
                raise ValueError("cannot specify a codomain curve and model name simultaneously")

            from sage.schemes.elliptic_curves.ell_field import compute_model
            codomain = compute_model(self._phis[-1].codomain(), model)

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
            sage: E = EllipticCurve(GF(43), [1,0])                                      # optional - sage.rings.finite_rings
            sage: P, = E.gens()                                                         # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P)                                # optional - sage.rings.finite_rings
            sage: psi = EllipticCurveHom_composite.from_factors(phi.factors())          # optional - sage.rings.finite_rings
            sage: psi == phi                                                            # optional - sage.rings.finite_rings
            True

        TESTS::

            sage: E = EllipticCurve('4730k1')
            sage: EllipticCurveHom_composite.from_factors([], E) == E.scalar_multiplication(1)
            True

        ::

            sage: E = EllipticCurve(GF(419), [1,0])                                     # optional - sage.rings.finite_rings
            sage: P, = E.gens()                                                         # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P)                                # optional - sage.rings.finite_rings
            sage: EllipticCurveHom_composite.from_factors(phi.factors()) == phi         # optional - sage.rings.finite_rings
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
            maps = (identity_morphism(E),)

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
            sage: K.<a> = NumberField(x^2 - x - 5)                                                  # optional - sage.rings.number_field
            sage: E = EllipticCurve('210.b6').change_ring(K)                                        # optional - sage.rings.number_field
            sage: psi = EllipticCurveHom_composite(E, E.torsion_points())                           # optional - sage.rings.number_field
            sage: R = E.lift_x(15/4 * (a+3))                                                        # optional - sage.rings.number_field
            sage: psi(R)    # indirect doctest                                                      # optional - sage.rings.number_field
            (1033648757/303450 : 58397496786187/1083316500*a - 62088706165177/2166633000 : 1)

        Check that copying the order over works::

            sage: E = EllipticCurve(GF(431), [1,0])                                                 # optional - sage.rings.finite_rings
            sage: P, = E.gens()                                                                     # optional - sage.rings.finite_rings
            sage: Q = 2^99*P; Q.order()                                                             # optional - sage.rings.finite_rings
            27
            sage: phi = E.isogeny(3^99*P, algorithm='factored')                                     # optional - sage.rings.finite_rings
            sage: phi(Q)._order                                                                     # optional - sage.rings.finite_rings
            27
        """
        return _eval_factored_isogeny(self._phis, P)

    def _eval(self, P):
        r"""
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self`` at the given point

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(j=Mod(1728,419))
            sage: K, = E.gens()
            sage: psi = EllipticCurveHom_composite(E, 4*K)
            sage: Ps = E.change_ring(GF(419**2))(0).division_points(5)                              # optional - sage.rings.finite_rings
            sage: {psi._eval(P).curve() for P in Ps}                                                # optional - sage.rings.finite_rings
            {Elliptic Curve defined by y^2 = x^3 + 373*x + 126 over Finite Field in z2 of size 419^2}
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')
        k = Sequence(P).universe()

        Q = P
        for phi in self._phis:
            Q = phi._eval(Q)

        return self._codomain.base_extend(k)(*Q)

    def _repr_(self):
        """
        Return basic facts about this composite isogeny as a string.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(43), [1,0])                                                  # optional - sage.rings.finite_rings
            sage: P, = E.gens()                                                                     # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P)                                            # optional - sage.rings.finite_rings
            sage: phi   # indirect doctest                                                          # optional - sage.rings.finite_rings
            Composite morphism of degree 44 = 2^2*11:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43
            sage: phi * phi * phi * phi * phi * phi * phi   # indirect doctest                      # optional - sage.rings.finite_rings
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
            sage: E = EllipticCurve(GF(43), [1,0])                                      # optional - sage.rings.finite_rings
            sage: P, = E.gens()                                                         # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P)                                # optional - sage.rings.finite_rings
            sage: phi.factors()                                                         # optional - sage.rings.finite_rings
            (Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43
                to Elliptic Curve defined by y^2 = x^3 + 39*x over Finite Field of size 43,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 39*x over Finite Field of size 43
                to Elliptic Curve defined by y^2 = x^3 + 42*x + 26 over Finite Field of size 43,
             Isogeny of degree 11
              from Elliptic Curve defined by y^2 = x^3 + 42*x + 26 over Finite Field of size 43
                to Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 43)
        """
        return self._phis

    # EllipticCurveHom methods

    @staticmethod
    def _composition_impl(left, right):
        """
        Helper method to compose other elliptic-curve morphisms with
        :class:`EllipticCurveHom_composite` objects. Called by
        :meth:`EllipticCurveHom._composition_`.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve([i + 1, i, 0, -4, -6*i])                                        # optional - sage.rings.number_field
            sage: P,Q = E.lift_x(i - 5), E.lift_x(-4*i)                                             # optional - sage.rings.number_field
            sage: phi = EllipticCurveHom_composite(E, P)                                            # optional - sage.rings.number_field
            sage: psi = phi.codomain().isogeny(phi(Q))                                              # optional - sage.rings.number_field
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: iso1 = WeierstrassIsomorphism(E, (-1, 0, -i - 1, 0))                              # optional - sage.rings.number_field
            sage: iso2 = psi.codomain().isomorphism_to(E)                                           # optional - sage.rings.number_field
            sage: psi * phi                 # indirect doctest                                      # optional - sage.rings.number_field
            Composite morphism of degree 16 = 2^2*4:
              From: Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I)
                    over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-3331/4)*x + (-142593/8*I)
                    over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: iso2 * EllipticCurveHom_composite.from_factors([phi, psi]) # indirect doctest     # optional - sage.rings.number_field
            Composite morphism of degree 16 = 4^2:
              From: Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I)
                    over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I)
                    over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: phi * iso1                # indirect doctest                                      # optional - sage.rings.number_field
            Composite morphism of degree 4 = 2^2:
              From: Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I)
                    over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (480*I-694)*x + (-7778*I+5556)
                    over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: iso2 * psi * phi * iso1   # indirect doctest                                      # optional - sage.rings.number_field
            Composite morphism of degree 16 = 2^2*4:
              From: Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I)
                    over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Elliptic Curve defined by y^2 + (I+1)*x*y = x^3 + I*x^2 + (-4)*x + (-6*I)
                    over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        if isinstance(left, EllipticCurveHom_composite):
            if isinstance(right, WeierstrassIsomorphism) and hasattr(left.factors()[0], '_set_pre_isomorphism'):    # XXX bit of a hack
                return EllipticCurveHom_composite.from_factors((left.factors()[0] * right,) + left.factors()[1:], strict=False)
            if isinstance(right, EllipticCurveHom_composite):
                return EllipticCurveHom_composite.from_factors(right.factors() + left.factors())
            if isinstance(right, EllipticCurveHom):
                return EllipticCurveHom_composite.from_factors((right,) + left.factors())
        if isinstance(right, EllipticCurveHom_composite):
            if isinstance(left, WeierstrassIsomorphism) and hasattr(right.factors()[-1], '_set_post_isomorphism'):  # XXX bit of a hack
                return EllipticCurveHom_composite.from_factors(right.factors()[:-1] + (left * right.factors()[-1],), strict=False)
            if isinstance(left, EllipticCurveHom):
                return EllipticCurveHom_composite.from_factors(right.factors() + (left,))
        return NotImplemented

    @staticmethod
    def _comparison_impl(left, right, op):
        r"""
        Compare a composite isogeny to another elliptic-curve morphism.

        Called by :meth:`EllipticCurveHom._richcmp_`.

        ALGORITHM:

        If possible, we use
        :func:`~sage.schemes.elliptic_curves.hom.compare_via_evaluation`.
        The complexity in that case is polynomial in the representation
        size of this morphism.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(QuadraticField(-3), [0,16])                                     # optional - sage.rings.number_field
            sage: P,Q = E.lift_x(0), E.lift_x(-4)                                                   # optional - sage.rings.number_field
            sage: phi = EllipticCurveHom_composite(E, P)                                            # optional - sage.rings.number_field
            sage: psi = phi.codomain().isogeny(phi(Q))                                              # optional - sage.rings.number_field
            sage: psi = psi.codomain().isomorphism_to(E) * psi                                      # optional - sage.rings.number_field
            sage: comp = psi * phi                                                                  # optional - sage.rings.number_field
            sage: mu = E.scalar_multiplication(phi.degree())                                        # optional - sage.rings.number_field
            sage: sum(a*comp == mu for a in E.automorphisms())                                      # optional - sage.rings.number_field
            1

        ::

            sage: E = EllipticCurve(GF(431**2), [1,0])                                              # optional - sage.rings.finite_rings
            sage: P,Q = E.gens()                                                                    # optional - sage.rings.finite_rings
            sage: phi1 = EllipticCurveHom_composite(E, P)                                           # optional - sage.rings.finite_rings
            sage: phi2 = EllipticCurveHom_composite(phi1.codomain(), phi1(Q))                       # optional - sage.rings.finite_rings
            sage: psi1 = EllipticCurveHom_composite(E, Q)                                           # optional - sage.rings.finite_rings
            sage: psi2 = EllipticCurveHom_composite(psi1.codomain(), psi1(P))                       # optional - sage.rings.finite_rings
            sage: phi2 * phi1 == psi2 * psi1                                                        # optional - sage.rings.finite_rings
            True
        """
        if op != op_EQ:
            return NotImplemented
        try:
            return compare_via_evaluation(left, right)
        except NotImplementedError:
            return NotImplemented

    def rational_maps(self):
        """
        Return the pair of explicit rational maps defining this composite
        isogeny.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])                             # optional - sage.rings.finite_rings
            sage: P = E.lift_x(7321)                                                    # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P)                                # optional - sage.rings.finite_rings
            sage: phi.rational_maps()                                                   # optional - sage.rings.finite_rings
            ((x^9 + 27463*x^8 + 21204*x^7 - 5750*x^6 + 1610*x^5 + 14440*x^4
              + 26605*x^3 - 15569*x^2 - 3341*x + 1267)/(x^8 + 27463*x^7 + 26871*x^6
              + 5999*x^5 - 20194*x^4 - 6310*x^3 + 24366*x^2 - 20905*x - 13867),
             (x^12*y + 8426*x^11*y + 5667*x^11 + 27612*x^10*y + 26124*x^10 + 9688*x^9*y
              - 22715*x^9 + 19864*x^8*y + 498*x^8 + 22466*x^7*y - 14036*x^7 + 8070*x^6*y
              + 19955*x^6 - 20765*x^5*y - 12481*x^5 + 12672*x^4*y + 24142*x^4 - 23695*x^3*y
              + 26667*x^3 + 23780*x^2*y + 17864*x^2 + 15053*x*y - 30118*x + 17539*y
              - 23609)/(x^12 + 8426*x^11 + 21945*x^10 - 22587*x^9 + 22094*x^8 + 14603*x^7
              - 26255*x^6 + 11171*x^5 - 16508*x^4 - 14435*x^3 - 2170*x^2 + 29081*x - 19009))

        TESTS::

            sage: f = phi.codomain().defining_polynomial()                              # optional - sage.rings.finite_rings
            sage: g = E.defining_polynomial().subs({2:1})                               # optional - sage.rings.finite_rings
            sage: f(*phi.rational_maps(), 1) % g                                        # optional - sage.rings.finite_rings
            0

        ::

            sage: phi.rational_maps()[0].parent()                                       # optional - sage.rings.finite_rings
            Fraction Field of
             Multivariate Polynomial Ring in x, y over Finite Field of size 65537
            sage: phi.rational_maps()[1].parent()                                       # optional - sage.rings.finite_rings
            Fraction Field of
             Multivariate Polynomial Ring in x, y over Finite Field of size 65537
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
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])                             # optional - sage.rings.finite_rings
            sage: P = E.lift_x(7321)                                                    # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P)                                # optional - sage.rings.finite_rings
            sage: phi.x_rational_map() == phi.rational_maps()[0]                        # optional - sage.rings.finite_rings
            True

        TESTS::

            sage: phi.x_rational_map().parent()                                         # optional - sage.rings.finite_rings
            Fraction Field of Univariate Polynomial Ring in x
             over Finite Field of size 65537
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
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])                             # optional - sage.rings.finite_rings
            sage: P = E.lift_x(7321)                                                    # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P); phi                           # optional - sage.rings.finite_rings
            Composite morphism of degree 9 = 3^2:
              From: Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
                    over Finite Field of size 65537
              To:   Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 28339*x + 59518
                    over Finite Field of size 65537
            sage: phi.kernel_polynomial()                                               # optional - sage.rings.finite_rings
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
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])                             # optional - sage.rings.finite_rings
            sage: P = E.lift_x(7321)                                                    # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P); phi                           # optional - sage.rings.finite_rings
            Composite morphism of degree 9 = 3^2:
              From: Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
                    over Finite Field of size 65537
              To:   Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 28339*x + 59518
                    over Finite Field of size 65537
            sage: psi = phi.dual(); psi                                                 # optional - sage.rings.finite_rings
            Composite morphism of degree 9 = 3^2:
              From: Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 28339*x + 59518
                    over Finite Field of size 65537
              To:   Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
                    over Finite Field of size 65537
            sage: psi * phi == phi.domain().scalar_multiplication(phi.degree())         # optional - sage.rings.finite_rings
            True
            sage: phi * psi == psi.domain().scalar_multiplication(psi.degree())         # optional - sage.rings.finite_rings
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
            sage: E = EllipticCurve(GF(7^2), [3,2])                                     # optional - sage.rings.finite_rings
            sage: P = E.lift_x(1)                                                       # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P); phi                           # optional - sage.rings.finite_rings
            Composite morphism of degree 7 = 7:
              From: Elliptic Curve defined by y^2 = x^3 + 3*x + 2
                    over Finite Field in z2 of size 7^2
              To:   Elliptic Curve defined by y^2 = x^3 + 3*x + 2
                    over Finite Field in z2 of size 7^2
            sage: phi.is_separable()                                                    # optional - sage.rings.finite_rings
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
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])                             # optional - sage.rings.finite_rings
            sage: P = E.lift_x(7321)                                                    # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P)                                # optional - sage.rings.finite_rings
            sage: phi.formal()                                                          # optional - sage.rings.finite_rings
            t + 54203*t^5 + 48536*t^6 + 40698*t^7 + 37808*t^8 + 21111*t^9 + 42381*t^10
             + 46688*t^11 + 657*t^12 + 38916*t^13 + 62261*t^14 + 59707*t^15
             + 30767*t^16 + 7248*t^17 + 60287*t^18 + 50451*t^19 + 38305*t^20
             + 12312*t^21 + 31329*t^22 + O(t^23)
            sage: (phi.dual() * phi).formal(prec=5)                                     # optional - sage.rings.finite_rings
            9*t + 65501*t^2 + 65141*t^3 + 59183*t^4 + 21491*t^5 + 8957*t^6
             + 999*t^7 + O(t^8)
        """
        res = self._phis[-1].formal(prec=prec)
        for phi in self._phis[:-1][::-1]:
            res = res(phi.formal(prec=prec))
        return res

    def scaling_factor(self):
        r"""
        Return the Weierstrass scaling factor associated to this
        composite morphism.

        The scaling factor is the constant `u` (in the base field)
        such that `\varphi^* \omega_2 = u \omega_1`, where
        `\varphi: E_1\to E_2` is this morphism and `\omega_i` are
        the standard Weierstrass differentials on `E_i` defined by
        `\mathrm dx/(2y+a_1x+a_3)`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])                             # optional - sage.rings.finite_rings
            sage: P = E.lift_x(7321)                                                    # optional - sage.rings.finite_rings
            sage: phi = EllipticCurveHom_composite(E, P)                                # optional - sage.rings.finite_rings
            sage: phi = WeierstrassIsomorphism(phi.codomain(), [7,8,9,10]) * phi        # optional - sage.rings.finite_rings
            sage: phi.formal()                                                          # optional - sage.rings.finite_rings
            7*t + 65474*t^2 + 511*t^3 + 61316*t^4 + 20548*t^5 + 45511*t^6 + 37285*t^7
             + 48414*t^8 + 9022*t^9 + 24025*t^10 + 35986*t^11 + 55397*t^12 + 25199*t^13
             + 18744*t^14 + 46142*t^15 + 9078*t^16 + 18030*t^17 + 47599*t^18
             + 12158*t^19 + 50630*t^20 + 56449*t^21 + 43320*t^22 + O(t^23)
            sage: phi.scaling_factor()                                                  # optional - sage.rings.finite_rings
            7

        ALGORITHM: The scaling factor is multiplicative under
        composition, so we return the product of the individual
        scaling factors associated to each factor.
        """
        return prod(phi.scaling_factor() for phi in self._phis)

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
            sage: E = EllipticCurve_from_j(GF(3).algebraic_closure()(0))                # optional - sage.rings.finite_rings
            sage: nu = EllipticCurveHom_composite.from_factors(E.automorphisms())       # optional - sage.rings.finite_rings
            sage: nu                                                                    # optional - sage.rings.finite_rings
            Composite morphism of degree 1 = 1^12:
              From: Elliptic Curve defined by y^2 = x^3 + x
                    over Algebraic closure of Finite Field of size 3
              To:   Elliptic Curve defined by y^2 = x^3 + x
                    over Algebraic closure of Finite Field of size 3
            sage: nu.is_injective()                                                     # optional - sage.rings.finite_rings
            True
        """
        return all(phi.is_injective() for phi in self._phis)

    @staticmethod
    def make_default():
        r"""
        This method does nothing and will be removed.

        (It is a leftover from the time when :class:`EllipticCurveHom_composite`
        wasn't the default yet.)

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            sage: EllipticCurveHom_composite.make_default()
            doctest:warning ...
        """
        from sage.misc.superseded import deprecation
        deprecation(34410, 'calling EllipticCurveHom_composite.make_default() is no longer necessary')
