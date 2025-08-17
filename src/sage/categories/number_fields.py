# sage_setup: distribution = sagemath-categories
r"""
Number fields
"""
# ****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.basic import Fields


class NumberFields(Category_singleton):
    r"""
    The category of number fields.

    EXAMPLES:

    We create the category of number fields::

        sage: C = NumberFields()
        sage: C
        Category of number fields

    By definition, it is infinite::

        sage: NumberFields().Infinite() is NumberFields()
        True

    Notice that the rational numbers `\QQ` *are* considered as
    an object in this category::

        sage: RationalField() in C
        True

    However, we can define a degree 1 extension of `\QQ`, which is of
    course also in this category::

        sage: x = PolynomialRing(RationalField(), 'x').gen()
        sage: K = NumberField(x - 1, 'a'); K                                            # needs sage.rings.number_field
        Number Field in a with defining polynomial x - 1
        sage: K in C                                                                    # needs sage.rings.number_field
        True

    Number fields all lie in this category, regardless of the name
    of the variable::

        sage: K = NumberField(x^2 + 1, 'a')                                             # needs sage.rings.number_field
        sage: K in C                                                                    # needs sage.rings.number_field
        True

    TESTS::

        sage: TestSuite(NumberFields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: NumberFields().super_categories()
            [Category of infinite fields]
        """
        return [Fields().Infinite()]

    def __contains__(self, x):
        r"""
        Return ``True`` if ``x`` is a number field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 1, 'a') in NumberFields()                           # needs sage.rings.number_field
            True
            sage: QuadraticField(-97, 'theta') in NumberFields()                        # needs sage.rings.number_field
            True
            sage: CyclotomicField(97) in NumberFields()                                 # needs sage.rings.number_field
            True

        Note that the rational numbers QQ are a number field::

            sage: QQ in NumberFields()
            True
            sage: ZZ in NumberFields()
            False
        """
        from sage.rings.number_field.number_field_base import NumberField
        return isinstance(x, NumberField)

    def _call_(self, x):
        r"""
        Construct an object in this category from the data in ``x``,
        or raise a :exc:`TypeError`.

        EXAMPLES::

            sage: C = NumberFields()
            sage: x = polygen(QQ, 'x')

            sage: C(QQ)
            Rational Field

            sage: C(NumberField(x^2 + 1, 'a'))                                          # needs sage.rings.number_field
            Number Field in a with defining polynomial x^2 + 1

            sage: C(UnitGroup(NumberField(x^2 + 1, 'a')))  # indirect doctest           # needs sage.rings.number_field
            Number Field in a with defining polynomial x^2 + 1

            sage: C(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: unable to canonically associate a number field to Integer Ring
        """
        try:
            return x.number_field()
        except AttributeError:
            raise TypeError("unable to canonically associate a number field to %s" % x)

    class ParentMethods:
        def zeta_function(self, prec=53,
                          max_imaginary_part=0,
                          max_asymp_coeffs=40, algorithm='pari'):
            r"""
            Return the Dedekind zeta function of this number field.

            Actually, this returns an interface for computing with the
            Dedekind zeta function `\zeta_F(s)` of the number field `F`.

            INPUT:

            - ``prec`` -- integer (default: 53); bits of precision

            - ``max_imaginary_part`` -- real (default: 0)

            - ``max_asymp_coeffs`` -- integer (default: 40)

            - ``algorithm`` -- (default: ``'pari'``) either ``'gp'`` or
              ``'pari'``

            OUTPUT: the zeta function of this number field

            If algorithm is ``'gp'``, this returns an interface to Tim
            Dokchitser's gp script for computing with `L`-functions.

            If algorithm is ``'pari'``, this returns instead an interface to Pari's
            own general implementation of `L`-functions.

            EXAMPLES::

                sage: K.<a> = NumberField(ZZ['x'].0^2 + ZZ['x'].0 - 1)                  # needs sage.rings.number_field
                sage: Z = K.zeta_function(); Z                                          # needs sage.rings.number_field sage.symbolic
                PARI zeta function associated to Number Field in a
                 with defining polynomial x^2 + x - 1
                sage: Z(-1)                                                             # needs sage.rings.number_field sage.symbolic
                0.0333333333333333

                sage: x = polygen(QQ, 'x')
                sage: L.<a, b, c> = NumberField([x^2 - 5, x^2 + 3, x^2 + 1])            # needs sage.rings.number_field
                sage: Z = L.zeta_function()                                             # needs sage.rings.number_field sage.symbolic
                sage: Z(5)                                                              # needs sage.rings.number_field sage.symbolic
                1.00199015670185

            Using the algorithm "pari"::

                sage: K.<a> = NumberField(ZZ['x'].0^2 + ZZ['x'].0 - 1)                  # needs sage.rings.number_field
                sage: Z = K.zeta_function(algorithm='pari')                             # needs sage.rings.number_field sage.symbolic
                sage: Z(-1)                                                             # needs sage.rings.number_field sage.symbolic
                0.0333333333333333

                sage: x = polygen(QQ, 'x')
                sage: L.<a, b, c> = NumberField([x^2 - 5, x^2 + 3, x^2 + 1])            # needs sage.rings.number_field
                sage: Z = L.zeta_function(algorithm='pari')                             # needs sage.rings.number_field sage.symbolic
                sage: Z(5)                                                              # needs sage.rings.number_field sage.symbolic
                1.00199015670185

            TESTS::

                sage: QQ.zeta_function()                                                # needs sage.symbolic
                PARI zeta function associated to Rational Field
            """
            if algorithm == 'gp':
                from sage.lfunctions.dokchitser import Dokchitser
                r1, r2 = self.signature()
                zero = [0]
                one = [1]
                Z = Dokchitser(conductor=abs(self.absolute_discriminant()),
                               gammaV=(r1 + r2) * zero + r2 * one,
                               weight=1,
                               eps=1,
                               poles=[1],
                               prec=prec)
                s = 'nf = nfinit(%s);' % self.absolute_polynomial()
                s += 'dzk = dirzetak(nf,cflength());'
                Z.init_coeffs('dzk[k]', pari_precode=s,
                              max_imaginary_part=max_imaginary_part,
                              max_asymp_coeffs=max_asymp_coeffs)
                Z.check_functional_equation()
                Z.rename('Dokchitser Zeta function associated to %s' % self)
                return Z

            if algorithm == 'pari':
                from sage.lfunctions.pari import lfun_number_field, LFunction
                Z = LFunction(lfun_number_field(self), prec=prec)
                Z.rename('PARI zeta function associated to %s' % self)
                return Z

            raise ValueError('algorithm must be "gp" or "pari"')

        def _test_absolute_disc(self, **options):
            r"""
            Run basic tests for the method :meth:`absolute_discriminant` of ``self``.

            See the documentation for :class:`TestSuite` for information on
            further options.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`

            EXAMPLES::

                sage: x = polygen(ZZ, 'x')
                sage: S = NumberField(x**3 - x - 1, 'a')                                # needs sage.rings.number_field
                sage: S._test_absolute_disc()                                           # needs sage.rings.number_field
            """
            from sage.rings.integer import Integer
            tester = self._tester(**options)
            tester.assertIsInstance(self.absolute_discriminant(), Integer)

    class ElementMethods:
        pass
