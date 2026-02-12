r"""
Supersingular Isogeny Diffie-Hellman

Toy implementation of the historical SIDH scheme.
This implementation serves as an example of a class that implements :class:`KeyExchangeBase`
rather than :class:`CommutativeKeyExchangeBase`.

.. WARNING::

    This is a toy implementation of a broken cryptographic scheme for educational
    use only! Do not use this implementation, or any cryptographic features of
    Sage, in any setting where security is needed!

AUTHORS:

- Taha Hedayat (2025-12-09): initial version
- Vincent Macri (2026-01-13): parameter set generation and cleanup
"""

# ****************************************************************************
#       Copyright (C) 2025 Taha Hedayat <tahah22121001@gmail.com>
#       Copyright (C) 2026 Vincent Macri <vincent.macri@ucalgary.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import annotations

import random
from typing import TYPE_CHECKING, Self

from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.constructor import EllipticCurve

from .key_exchange_base import KeyExchangeBase

if TYPE_CHECKING:
    from sage.rings.finite_rings.element_base import FinitePolyExtElement
    from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_finite_field
    from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

    PublicKeySIDH = tuple[
        EllipticCurve_finite_field,
        EllipticCurvePoint_finite_field,
        EllipticCurvePoint_finite_field,
    ]
    SecretKeySIDH = Integer | int


class SIDH(KeyExchangeBase):
    r"""
    Supersingular isogeny Diffie-Hellman key exchange.

    This implementation uses the notation of [Cos2020]_.

    For an easier-to-use constructor, see :meth:`parameter_set_from_prime`
    which implements parameter set generation as described in [Jao2022]_.

    EXAMPLES:

    This example comes from [Cos2020]_.::

        sage: e_A = 4
        sage: e_B = 3
        sage: p = 2^e_A * 3^e_B - 1
        sage: K.<i> = GF(p^2, modulus=x^2 + 1)
        sage: a0 = 329 * i + 423
        sage: E = EllipticCurve(K, [0, a0, 0, 1, 0])
        sage: PA = E(100 * i + 248, 304 * i + 199)
        sage: QA = E(426 * i + 394, 51 * i + 79)
        sage: PB = E(358 * i + 275, 410 * i + 104)
        sage: QB = E(20 * i + 185, 281 * i + 239)
        sage: toy_sidh = key_exchange.SIDH(E, PA, QA, PB, QB)
        doctest:...: FutureWarning: SageMath's key exchange functionality is experimental and might change in the future.
                     See https://github.com/sagemath/sage/issues/41218 for details.
        sage: TestSuite(toy_sidh).run()
    """

    def __init__(
        self,
        E: EllipticCurve_finite_field,
        PA: EllipticCurvePoint_finite_field,
        QA: EllipticCurvePoint_finite_field,
        PB: EllipticCurvePoint_finite_field,
        QB: EllipticCurvePoint_finite_field,
    ) -> None:
        K = E.base_field()
        self._p: Integer = K.characteristic()
        n: Integer = self._p + 1
        self._e_A: Integer = n.valuation(2)
        self._e_B: Integer = n.valuation(3)

        if 2**self._e_A * 3**self._e_B - 1 != self._p:
            raise ValueError('curve must be over characteristic 2^e_A * 3^e_B - 1')

        if K.degree() != 2:
            raise ValueError('curve must be defined over a degree two field extension')

        self._E = E

        def validate_point(P):
            if P not in self._E:
                raise ValueError(f'{P} is not on {self._E}')

        validate_point(PA)
        validate_point(QA)
        validate_point(PB)
        validate_point(QB)
        self._PA = PA
        self._QA = QA
        self._PB = PB
        self._QB = QB

    @classmethod
    def parameter_set_from_prime(cls, prime: Integer | int) -> Self:
        r"""
        Generate an SIDH instance from the specified prime, generating the parameter
        set as in [Jao2022]_.

        EXAMPLES::

            sage: key_exchange.SIDH.parameter_set_from_prime(431)
            SIDH with parameter set: (431, Elliptic Curve defined by y^2 = x^3 + 6*x^2 + x over Finite Field in i of size 431^2, ...
        """
        p: Integer = Integer(prime)
        if not p.is_prime():
            raise ValueError('p must be prime')

        n: Integer = p + ZZ.one()
        e_A: Integer = n.valuation(2)
        e_B: Integer = n.valuation(3)

        x = ZZ['x'].gen()
        f = x**3 + 6 * x**2 + x
        Fp = FiniteField(p)

        Fp2 = Fp.extension(modulus=x**2 + 1, name='i')
        i = Fp2.gen()

        E = EllipticCurve(Fp2, [0, 6, 0, 1, 0])

        # As defined in [Jao2022]_
        def canonical_square_root(s) -> FinitePolyExtElement:
            r: FinitePolyExtElement = Fp2(s).sqrt()
            alpha, beta = r.list()
            alpha = alpha.lift()
            beta = beta.lift()
            if (alpha != 0 and alpha % 2 == 0) or (alpha == 0 and beta % 2 == 0):
                assert r**2 == s
                return r
            assert (-r) ** 2 == s
            return -r

        def trial_points():
            for c in range(p):
                x0 = Fp2(i + c)
                y0 = Fp2(f.subs(x=x0))
                if not y0.is_square():
                    continue
                y0 = canonical_square_root(y0)
                if E.is_on_curve(x0, y0):
                    yield E(x0, y0)

        PA_mults = (
            E(-3 + 2 * canonical_square_root(2), 0),
            E(-3 - 2 * canonical_square_root(2), 0),
        )
        for R in trial_points():
            PA = (3**e_B) * R
            if 2 ** (e_A - 1) * PA in PA_mults:
                break
        else:
            raise ValueError(
                'failed to find PA, try using SIDH constructor directly and specify all parameters'
            )

        E_00 = E(0, 0)
        for R in trial_points():
            QA = (3**e_B) * R
            if 2 ** (e_A - 1) * QA == E_00:
                break
        else:
            raise ValueError(
                'failed to find QA, try using SIDH constructor directly and specify all parameters'
            )

        for c in range(p):
            fc = f.subs(x=c)
            if not Fp(fc).is_square():
                continue
            PB = 2 ** (e_A - 1) * E(c, canonical_square_root(fc))
            if PB.order() == 3**e_B:
                break
        else:
            raise ValueError(
                'failed to find PB, try using SIDH constructor directly and specify all parameters'
            )

        for c in range(p):
            fc = f.subs(x=c)
            if Fp(fc).is_square() or not Fp2(fc).is_square():
                continue
            QB = 2 ** (e_A - 1) * E(c, canonical_square_root(fc))
            if QB.order() == 3**e_B:
                break
        else:
            raise ValueError(
                'failed to find QB, try using SIDH constructor directly and specify all parameters'
            )
        return cls(E, PA, QA, PB, QB)

    @classmethod
    def named_parameter_set(cls, name: str) -> Self:
        r"""
        Return an SIDH instance corresponding to a named parameter set.

        INPUT:

        - ``name`` -- one of the following:

            - ``"toy"``: a small SIDH instance over `\GF{431^2}`,
                         this toy example is from [Cos2020]_.
            - ``"p434"``: SIDH instance with the same parameters
                          as SIKEp434 from [Jao2022]_.
            - ``"p503"``: SIDH instance with the same parameters
                          as SIKEp503 from [Jao2022]_.
            - ``"p610"``: SIDH instance with the same parameters
                          as SIKEp610 from [Jao2022]_.
            - ``"p751"``: SIDH instance with the same parameters
                          as SIKEp751 from [Jao2022]_.

        TESTS:

        Test that p434 gives the same public basis points as the SIKE specification::

            sage: p434 = key_exchange.SIDH.named_parameter_set('p434')
            sage: p, E, P2, Q2, P3, Q3 = p434.parameters()
            sage: xQ20, xQ21 = Q2[0].list()
            sage: xQ20 == 0x0000C7461738340EFCF09CE388F666EB38F7F3AFD42DC0B664D9F461F31AA2EDC6B4AB71BD42F4D7C058E13F64B237EF7DDD2ABC0DEB0C6C
            True
            sage: xQ21 == 0x000025DE37157F50D75D320DD0682AB4A67E471586FBC2D31AA32E6957FA2B2614C4CD40A1E27283EAAF4272AE517847197432E2D61C85F5
            True
            sage: yQ20, yQ21 = Q2[1].list()
            sage: yQ20 == 0x0001D407B70B01E4AEE172EDF491F4EF32144F03F5E054CEF9FDE5A35EFA3642A11817905ED0D4F193F31124264924A5F64EFE14B6EC97E5
            True
            sage: yQ21 == 0x0000E7DEC8C32F50A4E735A839DCDB89FE0763A184C525F7B7D0EBC0E84E9D83E9AC53A572A25D19E1464B509D97272AE761657B4765B3D6
            True
            sage: xP20, xP21 = P2[0].list()
            sage: xP20 == 0x00003CCFC5E1F050030363E6920A0F7A4C6C71E63DE63A0E6475AF621995705F7C84500CB2BB61E950E19EAB8661D25C4A50ED279646CB48
            True
            sage: xP21 == 0x0001AD1C1CAE7840EDDA6D8A924520F60E573D3B9DFAC6D189941CB22326D284A8816CC4249410FE80D68047D823C97D705246F869E3EA50
            True
            sage: yP20, yP21 = P2[1].list()
            sage: yP20 == 0x0001AB066B84949582E3F66688452B9255E72A017C45B148D719D9A63CDB7BE6F48C812E33B68161D5AB3A0A36906F04A6A6957E6F4FB2E0
            True
            sage: yP21 == 0x0000FD87F67EA576CE97FF65BF9F4F7688C4C752DCE9F8BD2B36AD66E04249AAF8337C01E6E4E1A844267BA1A1887B433729E1DD90C7DD2F
            True
            sage: xQ30, xQ31 = Q3[0].list()
            sage: xQ30 == 0x00012E84D7652558E694BF84C1FBDAAF99B83B4266C32EC65B10457BCAF94C63EB063681E8B1E7398C0B241C19B9665FDB9E1406DA3D3846
            True
            sage: xQ31 == 0
            True
            sage: yQ30, yQ31 = Q3[1].list()
            sage: yQ30 == 0
            True
            sage: yQ31 == 0x0000EBAAA6C731271673BEECE467FD5ED9CC29AB564BDED7BDEAA86DD1E0FDDF399EDCC9B49C829EF53C7D7A35C3A0745D73C424FB4A5FD2
            True
            sage: xP30, xP31 = P3[0].list()
            sage: xP30 == 0x00008664865EA7D816F03B31E223C26D406A2C6CD0C3D667466056AAE85895EC37368BFC009DFAFCB3D97E639F65E9E45F46573B0637B7A9
            True
            sage: xP31 == 0
            True
            sage: yP30, yP31 = P3[1].list()
            sage: yP30 == 0x00006AE515593E73976091978DFBD70BDA0DD6BCAEEBFDD4FB1E748DDD9ED3FDCF679726C67A3B2CC12B39805B32B612E058A4280764443B
            True
            sage: yP31 == 0
            True
        """
        if name == 'toy':
            R = ZZ['x']
            x = R.gen()
            p = 2**4 * 3**3 - 1
            K = FiniteField(p**2, name='i', modulus=x**2 + 1)
            i = K.gen()
            E = EllipticCurve(K, [0, 329 * i + 423, 0, 1, 0])
            PA = E(100 * i + 248, 304 * i + 199)
            QA = E(426 * i + 394, 51 * i + 79)
            PB = E(358 * i + 275, 410 * i + 104)
            QB = E(20 * i + 185, 281 * i + 239)
            toy = cls(E, PA, QA, PB, QB)
            toy.rename('sidh-toy')
            return toy
        elif name == 'p434':
            e_A = 0xD8
            e_B = 0x89
        elif name == 'p503':
            e_A = 0xFA
            e_B = 0x9F
        elif name == 'p610':
            e_A = 0x131
            e_B = 0xC0
        elif name == 'p751':
            e_A = 0x174
            e_B = 0xEF
        else:
            return super().named_parameter_set()

        p = 2**e_A * 3**e_B - 1
        sidh = cls.parameter_set_from_prime(p)
        sidh.rename(f'sidh-{name}')
        return sidh

    def parameters(
        self,
    ) -> tuple[
        Integer,
        EllipticCurve_finite_field,
        EllipticCurvePoint_finite_field,
        EllipticCurvePoint_finite_field,
        EllipticCurvePoint_finite_field,
        EllipticCurvePoint_finite_field,
    ]:
        r"""
        Return the parameter set of the SIDH instance.

        OUTPUT:

        A tuple (`E`, `PA`, `QA`, `PB`, `QB`) where:

        - `p` is the characteristic of the finite field `\GF{p^2}`
        - `E` is the starting curve
        - `PA` and `QA` are the generators for Alice's secret key
        - `PB` and `QB` are the generators for Bob's secret key
        """
        return (self._p, self._E, self._PA, self._QA, self._PB, self._QB)

    def alice_secret_key(self) -> Integer:
        r"""
        Generate Alice's secret key.

        EXAMPLES::

            sage: toy_sidh = key_exchange.SIDH.named_parameter_set('toy')
            sage: 0 <= toy_sidh.alice_secret_key() <= 2^4 - 1
            True
        """
        return Integer(random.randint(0, 2**self._e_A - 1))

    def bob_secret_key(self) -> Integer:
        r"""
        Generate Bob's secret key.

        EXAMPLES::

            sage: toy_sidh = key_exchange.SIDH.named_parameter_set('toy')
            sage: 0 <= toy_sidh.bob_secret_key() <= 3^3 - 1
            True
        """
        return Integer(random.randint(0, 3**self._e_B - 1))

    def alice_public_key(self, alice_secret_key: SecretKeySIDH) -> PublicKeySIDH:
        r"""
        Generate a valid public key for Alice.

        INPUT:

        - ``alice_secret_key`` -- Alice's secret key that will be used to generate
            the public key

        OUTPUT:

        Alice's public key as a tuple `(E_A, P'_B, Q'_B)`.
        """
        phi_A, _ = self.secret_isogeny_path(
            self._E, alice_secret_key, self._PA, self._QA
        )
        return phi_A.codomain(), phi_A(self._PB), phi_A(self._QB)

    def bob_public_key(self, bob_secret_key: SecretKeySIDH) -> PublicKeySIDH:
        r"""
        Generate a valid public key for Alice.

        INPUT:

        - ``alice_secret_key`` -- Alice's secret key that will be used to generate
            the public key

        OUTPUT:

        Bob's public key as a tuple `(E_B, P'_A, Q'_A)`.
        """
        phi_B, _ = self.secret_isogeny_path(self._E, bob_secret_key, self._PB, self._QB)
        return phi_B.codomain(), phi_B(self._PA), phi_B(self._QA)

    def alice_compute_shared_secret(
        self, alice_secret_key: SecretKeySIDH, bob_public_key: PublicKeySIDH
    ) -> Integer:
        r"""
        Compute the shared secret using Alice's secret key and Bob's public key.
        """
        E_B, P_A1, Q_A1 = bob_public_key
        phi_A1, _ = self.secret_isogeny_path(E_B, alice_secret_key, P_A1, Q_A1)
        E_AB = phi_A1.codomain()
        return E_AB.j_invariant()

    def bob_compute_shared_secret(
        self, bob_secret_key: SecretKeySIDH, alice_public_key: PublicKeySIDH
    ) -> Integer:
        r"""
        Compute the shared secret using Bob's secret key and Alice's public key.
        """
        E_A, P_B1, Q_B1 = alice_public_key
        phi_B1, _ = self.secret_isogeny_path(E_A, bob_secret_key, P_B1, Q_B1)
        E_BA = phi_B1.codomain()
        return E_BA.j_invariant()

    def secret_isogeny_path(
        self,
        start_curve: EllipticCurve_finite_field,
        secret_key: SecretKeySIDH,
        P: EllipticCurvePoint_finite_field,
        Q: EllipticCurvePoint_finite_field,
    ) -> tuple[EllipticCurveHom_composite, list[Integer]]:
        r"""
        Compute the secret isogeny and path of j-invariants.
        The kernel of returned isogeny is generated by ``P + secret_key * Q``.

        The purpose of this method is to demonstrate how implementations of
        :class:`~sage.crypto.public_key.key_exchange_base.KeyExchangeBase` may
        implement additional methods besides those defined in
        :class:`~sage.crypto.public_key.key_exchange_base.KeyExchangeBase`
        that compute some additional interesting information, in this case
        the path of j-invariants.

        INPUT:

        - ``start_curve`` -- the starting curve
        - ``secret_key`` -- the secret key
        - ``P`` -- first basis point
        - ``Q`` -- second basis point

        OUTPUT:

        A tuple ``(iso, path)`` where ``iso`` is an elliptic curve isogeny
        with kernel generated by ``P + secret_key * Q`` and ``path`` is
        a list of the j-invariants of the isomorphism classes of the curves on
        the isogeny path.

        EXAMPLES:

        We work through the key generation example from [Cos2020]_::

            sage: toy_sidh = key_exchange.SIDH.named_parameter_set('toy')
            sage: p, E0, PA, QA, PB, QB = toy_sidh.parameters()
            sage: alice_secret_key = 11
            sage: phi_A, alice_keygen_path = toy_sidh.secret_isogeny_path(E0, alice_secret_key, PA, QA)
            sage: alice_keygen_path
            [87*i + 190, 107, 344*i + 190, 350*i + 65, 222*i + 118]
            sage: bob_secret_key = 2
            sage: phi_B, bob_keygen_path = toy_sidh.secret_isogeny_path(E0, bob_secret_key, PB, QB)
            sage: bob_keygen_path
            [87*i + 190, 106*i + 379, 325*i + 379, 344*i + 190]
        """
        iso = start_curve.isogeny(P + secret_key * Q, algorithm='factored')
        j_invariant_path = [start_curve.j_invariant()]
        j_invariant_path.extend(E.codomain().j_invariant() for E in iso.factors())
        return iso, j_invariant_path
