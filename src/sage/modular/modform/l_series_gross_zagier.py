from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.lfunctions.dokchitser import Dokchitser
from .l_series_gross_zagier_coeffs import gross_zagier_L_series
from sage.modular.dirichlet import kronecker_character


class GrossZagierLseries(SageObject):

    def __init__(self, E, A, prec=53):
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

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series_gross_zagier import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)

        TESTS::

            sage: K.<b> = QuadraticField(131)
            sage: A = K.class_group().one()
            sage: G = GrossZagierLseries(e, A)
            Traceback (most recent call last):
            ...
            ValueError: A is not an ideal class in an imaginary quadratic field
        """
        self._E = E
        self._N = N = E.conductor()
        self._A = A
        ideal = A.ideal()
        K = A.gens()[0].parent()
        D = K.disc()
        if not (K.degree() == 2 and D < 0):
            raise ValueError("A is not an ideal class in an"
                             " imaginary quadratic field")
        Q = ideal.quadratic_form().reduced_form()
        epsilon = - kronecker_character(D)(N)
        self._dokchister = Dokchitser(N ** 2 * D ** 2,
                                      [0, 0, 1, 1],
                                      weight=2, eps=epsilon, prec=prec)
        self._nterms = nterms = Integer(self._dokchister.num_coeffs())
        if nterms > 1e6:
            # just takes way to long
            raise ValueError("Too many terms: {}".format(nterms))
        zeta_ord = ideal.number_field().zeta_order()
        an_list = gross_zagier_L_series(E.anlist(nterms + 1), Q, N, zeta_ord)
        self._dokchister.gp().set('a', an_list[1:])
        self._dokchister.init_coeffs('a[k]', 1)

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
            -0.272946890617590
            sage: G(3, 1)
            0.212442670030741
        """
        return self._dokchister.derivative(s, der)

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
            -0.613002046122894 + 0.490374999263514*z - 0.122903033710382*z^2 + O(z^3)
        """
        return self._dokchister.taylor_series(s, series_prec, var)

    def _repr_(self):
        """
        Return the string representation.

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series_gross_zagier import GrossZagierLseries
            sage: GrossZagierLseries(e, A)
            Gross Zagier L-series attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field with ideal class Fractional ideal class (2, 1/2*a)
        """
        msg = "Gross Zagier L-series attached to {} with ideal class {}"
        return msg.format(self._E, self._A)
