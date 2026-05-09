"""
Gross-Zagier L-series
"""
from sage.lfunctions.pari import lfun_generic, LFunction
from sage.modular.dirichlet import kronecker_character
from sage.modular.modform.l_series_gross_zagier_coeffs import gross_zagier_L_series
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject


class GrossZagierLseries(SageObject):

    def __init__(self, E, A, prec=53, max_imaginary_part=0) -> None:
        r"""
        Class for the Gross-Zagier `L`-series.

        This is attached to a pair `(E,A)` where `E` is an elliptic curve over
        `\QQ` and `A` is an ideal class in an imaginary quadratic number field.

        For the exact definition, in the more general setting of modular forms
        instead of elliptic curves, see section IV of [GZ1986]_.

        INPUT:

        - ``E`` -- an elliptic curve over `\QQ`

        - ``A`` -- an ideal class in an imaginary quadratic number field

        - ``prec`` -- integer (default: 53); giving the required precision

        - ``max_imaginary_part`` -- real number (default: 0)

        EXAMPLES::

            sage: from sage.modular.modform.l_series_gross_zagier import GrossZagierLseries
            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: G = GrossZagierLseries(e, A); G
            Gross Zagier L-series attached to Elliptic Curve defined by
            y^2 + y = x^3 - x over Rational Field with ideal class
            Fractional ideal class (2, 1/2*a)

        TESTS::

            sage: K.<b> = QuadraticField(131)
            sage: A = K.class_group().one()
            sage: G = GrossZagierLseries(e, A)
            Traceback (most recent call last):
            ...
            ValueError: A is not an ideal class in an imaginary quadratic field
        """
        self._N = N = E.conductor()
        ideal = A.ideal()
        K = A.gens()[0].parent()
        D = K.disc()
        if not (K.degree() == 2 and D < 0):
            raise ValueError("A is not an ideal class in an"
                             " imaginary quadratic field")
        Q = ideal.quadratic_form().reduced_form()
        epsilon = - kronecker_character(D)(N)

        # first compute the number of required terms
        Lpari = lfun_generic(N**2 * D**2,
                             [0, 0, 1, 1],
                             weight=2, eps=epsilon)
        L = LFunction(Lpari, prec=prec, max_im=max_imaginary_part)
        nterms = Integer(L.cost())
        if nterms > 1e6:
            # just takes way to long
            raise ValueError(f"Too many terms: {nterms}")

        zeta_ord = ideal.number_field().zeta_order()
        an_list = gross_zagier_L_series(E.anlist(nterms + 1), Q, N, zeta_ord)
        Lpari.init_coeffs(an_list[1:])
        self._lfunction = LFunction(Lpari, prec=prec)

        msg = f"Gross Zagier L-series attached to {E} with ideal class {A}"
        self.rename(msg)

    def __call__(self, s, der=0):
        r"""
        Return the value at `s`.

        INPUT:

        - ``s`` -- complex number

        - ``der`` -- order of derivative (default: 0)

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series_gross_zagier import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)
            sage: G(3)
            -0.272946890617449
            sage: G(3, 1)
            0.212442670030197
        """
        return self._lfunction.derivative(s, der)

    def taylor_series(self, s=1, series_prec=6, var='z'):
        r"""
        Return the Taylor series at `s`.

        INPUT:

        - ``s`` -- complex number (default: 1)
        - ``series_prec`` -- number of terms (default: 6) in the Taylor series
        - ``var`` -- variable (default: ``'z'``)

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series_gross_zagier import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)
            sage: G.taylor_series(2,3)
            -0.613002046122888 + 0.490374999263489*z - 0.122903033710347*z^2 + O(z^3)
        """
        return self._lfunction.taylor_series(s, series_prec, var)
