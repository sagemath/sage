"""
Local Generic

Superclass for `p`-adic and power series rings.

AUTHORS:

- David Roe
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from sage.rings.ring import CommutativeRing
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationRings, CompleteDiscreteValuationFields
from sage.structure.category_object import check_default_category
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

class LocalGeneric(CommutativeRing):
    def __init__(self, base, prec, names, element_class, category=None):
        """
        Initializes self.

        EXAMPLES::

            sage: R = Zp(5) #indirect doctest
            sage: R.precision_cap()
            20

        In :trac:`14084`, the category framework has been implemented for p-adic rings::

            sage: TestSuite(R).run()
            sage: K = Qp(7)
            sage: TestSuite(K).run()

        TESTS::

            sage: R = Zp(5, 5, 'fixed-mod')
            sage: R._repr_option('element_is_atomic')
            False
        """
        self._prec = prec
        self.Element = element_class
        default_category = getattr(self, '_default_category', None)
        if self.is_field():
            category = CompleteDiscreteValuationFields()
        else:
            category = CompleteDiscreteValuationRings()
        category = category.Metric().Complete()
        if default_category is not None:
            category = check_default_category(default_category, category)
        Parent.__init__(self, base, names=(names,), normalize=False, category=category)

    def is_capped_relative(self):
        """
        Returns whether this `p`-adic ring bounds precision in a capped
        relative fashion.

        The relative precision of an element is the power of `p`
        modulo which the unit part of that element is defined.  In a
        capped relative ring, the relative precision of elements are
        bounded by a constant depending on the ring.

        EXAMPLES::

            sage: R = ZpCA(5, 15)
            sage: R.is_capped_relative()
            False
            sage: R(5^7)
            5^7 + O(5^15)
            sage: S = Zp(5, 15)
            sage: S.is_capped_relative()
            True
            sage: S(5^7)
            5^7 + O(5^22)
        """
        return False

    def is_capped_absolute(self):
        """
        Returns whether this `p`-adic ring bounds precision in a
        capped absolute fashion.

        The absolute precision of an element is the power of `p`
        modulo which that element is defined.  In a capped absolute
        ring, the absolute precision of elements are bounded by a
        constant depending on the ring.

        EXAMPLES::

            sage: R = ZpCA(5, 15)
            sage: R.is_capped_absolute()
            True
            sage: R(5^7)
            5^7 + O(5^15)
            sage: S = Zp(5, 15)
            sage: S.is_capped_absolute()
            False
            sage: S(5^7)
            5^7 + O(5^22)
        """
        return False

    def is_fixed_mod(self):
        """
        Returns whether this `p`-adic ring bounds precision in a fixed
        modulus fashion.

        The absolute precision of an element is the power of `p`
        modulo which that element is defined.  In a fixed modulus
        ring, the absolute precision of every element is defined to be
        the precision cap of the parent.  This means that some
        operations, such as division by `p`, don't return a well defined
        answer.

        EXAMPLES::

            sage: R = ZpFM(5,15)
            sage: R.is_fixed_mod()
            True
            sage: R(5^7,absprec=9)
            5^7 + O(5^15)
            sage: S = ZpCA(5, 15)
            sage: S.is_fixed_mod()
            False
            sage: S(5^7,absprec=9)
            5^7 + O(5^9)
        """
        return False

    def is_floating_point(self):
        """
        Returns whether this `p`-adic ring bounds precision in a floating
        point fashion.

        The relative precision of an element is the power of `p`
        modulo which the unit part of that element is defined.  In a
        floating point ring, elements do not store precision, but arithmetic
        operations truncate to a relative precision depending on the ring.

        EXAMPLES::

            sage: R = ZpCR(5, 15)
            sage: R.is_floating_point()
            False
            sage: R(5^7)
            5^7 + O(5^22)
            sage: S = ZpFP(5, 15)
            sage: S.is_floating_point()
            True
            sage: S(5^7)
            5^7
        """
        return False

    def is_lazy(self):
        """
        Returns whether this `p`-adic ring bounds precision in a lazy
        fashion.

        In a lazy ring, elements have mechanisms for computing
        themselves to greater precision.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.is_lazy()
            False
        """
        return False

    def _latex_(self):
        r"""
        Latex.

        EXAMPLES::

            sage: latex(Zq(27,names='a')) #indirect doctest
            \mathbf{Z}_{3^{3}}
        """
        return self._repr_(do_latex = True)

    def change(self, **kwds):
        r"""
        Return a new ring with changed attributes.

        INPUT:

        The following arguments are applied to every ring in the tower:

        - ``type`` -- string, the precision type
        - ``p`` -- the prime of the ground ring.  Defining polynomials
                   will be converted to the new base rings.
        - ``print_mode`` -- string
        - ``print_pos`` -- bool
        - ``print_sep`` -- string
        - ``print_alphabet`` -- dict
        - ``show_prec`` -- bool
        - ``check`` -- bool

        The following arguments are only applied to the top ring in the tower:

        - ``var_name`` -- string
        - ``res_name`` -- string
        - ``unram_name`` -- string
        - ``ram_name`` -- string
        - ``names`` -- string
        - ``modulus`` -- polynomial

        The following arguments have special behavior:

        - ``prec`` -- integer.  If the precision is increased on an extension ring,
                       the precision on the base is increased as necessary (respecting ramification).
                       If the precision is decreased, the precision of the base is unchanged.

        - ``field`` -- bool.  If ``True``, switch to a tower of fields via the fraction field.
                        If False, switch to a tower of rings of integers.

        - ``q`` -- prime power.  Replace the initial unramified extension of `\QQ_p` or `\ZZ_p`
                    with an unramified extension of residue cardinality `q`.
                    If the initial extension is ramified, add in an unramified extension.

        - ``base`` -- ring or field. Use a specific base ring instead of recursively
                       calling :meth:`change` down the tower.

        See the :mod:`constructors <sage.rings.padics.factory>` for more details on the
        meaning of these arguments.

        EXAMPLES:

        We can use this method to change the precision::

            sage: Zp(5).change(prec=40)
            5-adic Ring with capped relative precision 40

        or the precision type::

            sage: Zp(5).change(type="capped-abs")
            5-adic Ring with capped absolute precision 20

        or even the prime::

            sage: ZpCA(3).change(p=17)
            17-adic Ring with capped absolute precision 20

        You can switch between the ring of integers and its fraction field::

            sage: ZpCA(3).change(field=True)
            3-adic Field with capped relative precision 20

        You can also change print modes::

            sage: R = Zp(5).change(prec=5, print_mode='digits')
            sage: repr(~R(17))
            '...13403'

        Changing print mode to 'digits' works for Eisenstein extensions::

            sage: S.<x> = ZZ[]
            sage: W.<w> = Zp(3).extension(x^4 + 9*x^2 + 3*x - 3)
            sage: W.print_mode()
            'series'
            sage: W.change(print_mode='digits').print_mode()
            'digits'

        You can change extensions::

            sage: K.<a> = QqFP(125, prec=4)
            sage: K.change(q=64)
            Unramified Extension in a defined by x^6 + x^4 + x^3 + x + 1 with floating precision 4 over 2-adic Field
            sage: R.<x> = QQ[]
            sage: K.change(modulus = x^2 - x + 2, print_pos=False)
            Unramified Extension in a defined by x^2 - x + 2 with floating precision 4 over 5-adic Field

        and variable names::

            sage: K.change(names='b')
            Unramified Extension in b defined by x^3 + 3*x + 3 with floating precision 4 over 5-adic Field

        and precision::

            sage: Kup = K.change(prec=8); Kup
            Unramified Extension in a defined by x^3 + 3*x + 3 with floating precision 8 over 5-adic Field
            sage: Kup.base_ring()
            5-adic Field with floating precision 8

        If you decrease the precision, the precision of the base stays the same::

            sage: Kdown = K.change(prec=2); Kdown
            Unramified Extension in a defined by x^3 + 3*x + 3 with floating precision 2 over 5-adic Field
            sage: Kdown.precision_cap()
            2
            sage: Kdown.base_ring()
            5-adic Field with floating precision 4

        Changing the prime works for extensions::

            sage: x = polygen(ZZ)
            sage: R.<a> = Zp(5).extension(x^2 + 2)
            sage: S = R.change(p=7)
            sage: S.defining_polynomial(exact=True)
            x^2 + 2
            sage: A.<y> = Zp(5)[]
            sage: R.<a> = Zp(5).extension(y^2 + 2)
            sage: S = R.change(p=7)
            sage: S.defining_polynomial(exact=True)
            y^2 + 2

        ::

            sage: R.<a> = Zq(5^3)
            sage: S = R.change(prec=50)
            sage: S.defining_polynomial(exact=True)
            x^3 + 3*x + 3
        """
        # We support both print_* and * for *=mode, pos, sep, alphabet
        for atr in ('print_mode', 'print_pos', 'print_sep', 'print_alphabet'):
            if atr in kwds:
                kwds[atr[6:]] = kwds.pop(atr)
        def get_unramified_modulus(q, res_name):
            from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
            return GF(q, res_name).modulus().change_ring(ZZ)
        n = None
        q = None
        from .padic_base_generic import pAdicBaseGeneric
        if 'q' in kwds and isinstance(self.base_ring(), pAdicBaseGeneric):
            q = kwds.pop('q')
            if not isinstance(q, Integer):
                raise TypeError("q must be an integer")
            p, n = q.is_prime_power(get_data=True)
            if n == 0:
                raise ValueError("q must be a prime power")
            if 'p' in kwds and kwds['p'] != p:
                raise ValueError("q does not match p")
            kwds['p'] = p
        functor, ring = self.construction()
        functor = copy(functor)
        if 'mode' in kwds and 'show_prec' not in kwds:
            new_type = kwds.get('type', self._prec_type())
            cur_type = self._prec_type()
            cur_mode = self._printer._print_mode()
            cur_show_prec = self._printer._show_prec()
            from .factory import _default_show_prec
            if cur_show_prec == _default_show_prec(cur_type, cur_mode):
                kwds['show_prec'] = _default_show_prec(new_type, kwds['mode'])
            else:
                raise RuntimeError
        p = kwds.get('p', functor.p if hasattr(functor, 'p') else self.prime())
        curpstr = str(self.prime())
        functor_dict = getattr(functor, "extras", getattr(functor, "kwds", None))
        # If we are switching to 'digits', or changing p, need to ensure a large enough alphabet.
        if 'alphabet' not in kwds and (kwds.get('mode') == 'digits' or
           (functor_dict['print_mode'].get('mode') == 'digits' and p > getattr(functor, "p", p))):
            from .padic_printing import _printer_defaults
            kwds['alphabet'] = _printer_defaults.alphabet()[:p]
        # For fraction fields of fixed-mod rings, we need to explicitly set show_prec = False
        if 'field' in kwds and 'type' not in kwds:
            if self._prec_type() == 'capped-abs':
                kwds['type'] = 'capped-rel'
            elif self._prec_type() == 'fixed-mod':
                kwds['type'] = 'floating-point'
                kwds['show_prec'] = False # This can be removed once printing of fixed mod elements is changed.

        # There are two kinds of functors possible:
        # CompletionFunctor and AlgebraicExtensionFunctor
        # We distinguish them by the presence of ``prec``,
        if hasattr(functor, "prec"):
            functor.extras = copy(functor.extras)
            functor.extras['print_mode'] = copy(functor.extras['print_mode'])
            if 'type' in kwds and kwds['type'] not in functor._dvr_types:
                raise ValueError("completion type must be one of %s"%(", ".join(functor._dvr_types[1:])))
            if 'field' in kwds:
                field = kwds.pop('field')
                if field:
                    ring = ring.fraction_field()
                elif ring.is_field():
                    ring = ring.ring_of_integers()
            for atr in ('p', 'prec', 'type'):
                if atr in kwds:
                    setattr(functor, atr, kwds.pop(atr))
            if q is not None:
                if 'names' in kwds:
                    names = kwds.pop('names')
                elif 'unram_name' in kwds:
                    names = kwds.pop('unram_name')
                else:
                    raise TypeError("You must specify the name of the generator")
                res_name = kwds.pop('res_name', names + '0')
                modulus = kwds.pop('modulus', get_unramified_modulus(q, res_name))
                implementation = kwds.pop('implementation', 'FLINT')
            # We have to change the way p prints in the default case
            if 'names' in kwds:
                functor.extras['names'] = kwds.pop('names')
            elif functor.extras['names'][0] == curpstr:
                functor.extras['names'] = (str(p),)
            for atr in ('ram_name', 'var_name'):
                if atr in kwds:
                    functor.extras['print_mode'][atr] = kwds.pop(atr)
                elif functor.extras['print_mode'].get(atr) == curpstr:
                    functor.extras['print_mode'][atr] = str(p)
            if 'check' in kwds:
                functor.extras['check'] = kwds.pop('check')
            for atr in ('mode', 'pos', 'unram_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet', 'show_prec'):
                if atr in kwds:
                    functor.extras['print_mode'][atr] = kwds.pop(atr)
            if kwds:
                raise ValueError("Extra arguments received: %s"%(", ".join(kwds.keys())))
            if q is not None:
                # Create an unramified extension
                base = functor(ring)
                from .factory import ExtensionFactory
                modulus = modulus.change_ring(base)
                return ExtensionFactory(base=base, premodulus=modulus, names=names, res_name=res_name, unram=True, implementation=implementation)
        else:
            functor.kwds = copy(functor.kwds)
            functor.kwds['print_mode'] = copy(functor.kwds['print_mode'])
            if 'prec' in kwds:
                prec = kwds.pop('prec')
                baseprec = (prec - 1) // self.e() + 1
                if baseprec > self.base_ring().precision_cap():
                    kwds['prec'] = baseprec
                functor.kwds['prec'] = prec
            from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
            if 'names' in kwds:
                functor.names = [kwds.pop('names')]
            modulus = None
            if 'modulus' in kwds:
                modulus = kwds.pop('modulus')
                if n is not None and modulus.degree() != n:
                    raise ValueError("modulus must have degree matching q")
            elif q is not None and self.f() != 1:
                # If q is specified, replace the modulus with one from q.
                modulus = get_unramified_modulus(q, functor.kwds.get('res_name', functor.names[0] + '0'))
            for atr in ('var_name', 'res_name', 'unram_name', 'ram_name'):
                if atr in kwds:
                    functor.kwds[atr] = kwds.pop(atr)
            if 'check' in kwds:
                functor.kwds['check'] = kwds['check']
            for atr in ('mode', 'pos', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet', 'show_prec'):
                if atr in kwds:
                    functor.kwds['print_mode'][atr] = kwds[atr]
            if 'base' in kwds:
                ring = kwds['base']
            else:
                if q is not None and self.f() == 1:
                    kwds['q'] = q
                ring = ring.change(**kwds)
            if modulus is None:
                if len(functor.polys) != 1:
                    raise RuntimeError("Unexpected number of defining polynomials")
                modulus = functor.polys[0]
            if isinstance(modulus.base_ring(), pAdicBaseGeneric):
                modulus.change_ring(ring)
            functor.polys = [modulus]
        return functor(ring)

    def precision_cap(self):
        r"""
        Returns the precision cap for ``self``.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - integer -- ``self``'s precision cap

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.precision_cap()
            10
            sage: R = Zp(3, 10,'capped-rel'); R.precision_cap()
            10
            sage: R = Zp(3, 10,'capped-abs'); R.precision_cap()
            10

        .. NOTE::

            This will have different meanings depending on the type of
            local ring.  For fixed modulus rings, all elements are
            considered modulo ``self.prime()^self.precision_cap()``.
            For rings with an absolute cap (i.e. the class
            ``pAdicRingCappedAbsolute``), each element has a precision
            that is tracked and is bounded above by
            ``self.precision_cap()``.  Rings with relative caps
            (e.g. the class ``pAdicRingCappedRelative``) are the same
            except that the precision is the precision of the unit
            part of each element.  For lazy rings, this gives the
            initial precision to which elements are computed.
        """
        return self._prec

    def is_exact(self):
        r"""
        Returns whether this p-adic ring is exact, i.e. False.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            boolean -- whether self is exact, i.e. False.

        EXAMPLES::

            #sage: R = Zp(5, 3, 'lazy'); R.is_exact()
            #False
            sage: R = Zp(5, 3, 'fixed-mod'); R.is_exact()
            False
        """
        return False

    def residue_characteristic(self):
        r"""
        Returns the characteristic of ``self``'s residue field.

        INPUT:

        - ``self`` -- a p-adic ring.

        OUTPUT:

        - integer -- the characteristic of the residue field.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.residue_characteristic()
            3
        """
        return self.residue_class_field().characteristic()

    def defining_polynomial(self, var = 'x'):
        r"""
        Returns the defining polynomial of this local ring, i.e. just ``x``.

        INPUT:

        - ``self`` -- a local ring
        - ``var`` -- string (default: ``'x'``) the name of the variable

        OUTPUT:

        - polynomial -- the defining polynomial of this ring as an extension over its ground ring

        EXAMPLES::

            sage: R = Zp(3, 3, 'fixed-mod'); R.defining_polynomial('foo')
            (1 + O(3^3))*foo + (O(3^3))
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(self, var).gen()

    def ground_ring(self):
        r"""
        Returns ``self``.

        Will be overridden by extensions.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - the ground ring of ``self``, i.e., itself

        EXAMPLES::

            sage: R = Zp(3, 5, 'fixed-mod')
            sage: S = Zp(3, 4, 'fixed-mod')
            sage: R.ground_ring() is R
            True
            sage: S.ground_ring() is R
            False
        """
        return self

    def ground_ring_of_tower(self):
        r"""
        Returns ``self``.

        Will be overridden by extensions.

        INPUT:

        - ``self`` -- a `p`-adic ring

        OUTPUT:

        - the ground ring of the tower for ``self``, i.e., itself

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.ground_ring_of_tower()
            5-adic Ring with capped relative precision 20
        """
        return self

    def degree(self):
        r"""
        Returns the degree of ``self`` over the ground ring, i.e. 1.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - integer -- the degree of this ring, i.e., 1

        EXAMPLES::

            sage: R = Zp(3, 10, 'capped-rel'); R.degree()
            1
        """
        return Integer(1)

    def ramification_index(self, K = None):
        r"""
        Returns the ramification index over the ground ring: 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - integer -- the ramification index of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.ramification_index()
            1
        """
        if K is None or K is self:
            return Integer(1)
        else:
            raise ValueError("K should be a subring of self")

    def e(self, K = None):
        r"""
        Returns the ramification index over the ground ring: 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring
        - ``K`` -- a subring of ``self`` (default ``None``)

        OUTPUT:

        - integer -- the ramification index of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.e()
            1
        """
        return self.ramification_index(K)

    def inertia_degree(self, K=None):
        r"""
        Returns the inertia degree over ``K`` (defaults to the ground ring): 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring
        - ``K`` -- a subring of ``self`` (default None)

        OUTPUT:

        - integer -- the inertia degree of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.inertia_degree()
            1
        """
        return Integer(1)

    def residue_class_degree(self, K=None):
        r"""
        Returns the inertia degree over the ground ring: 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring
        - ``K`` -- a subring (default ``None``)

        OUTPUT:

        - integer -- the inertia degree of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.residue_class_degree()
            1
        """
        return self.inertia_degree(K)

    def f(self, K=None):
        r"""
        Returns the inertia degree over the ground ring: 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring
        - ``K`` -- a subring (default ``None``)

        OUTPUT:

        - integer -- the inertia degree of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.f()
            1
        """
        return self.inertia_degree(K)

    def inertia_subring(self):
        r"""
        Returns the inertia subring, i.e. ``self``.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - the inertia subring of self, i.e., itself

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.inertia_subring()
            5-adic Ring with capped relative precision 20
        """
        return self

    def maximal_unramified_subextension(self):
        r"""
        Returns the maximal unramified subextension.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - the maximal unramified subextension of ``self``

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.maximal_unramified_subextension()
            5-adic Ring with capped relative precision 20
        """
        return self.inertia_subring()

#    def get_extension(self):
#        r"""
#        Returns the trivial extension of self.
#        """
#        raise NotImplementedError

    def uniformiser(self):
        """
        Returns a uniformiser for ``self``, ie a generator for the unique maximal ideal.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.uniformiser()
            5 + O(5^21)
            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.uniformiser()
            t + O(t^21)
        """
        return self.uniformizer()

    def uniformiser_pow(self, n):
        """
        Returns the `n`th power of the uniformiser of ``self`` (as an element of ``self``).

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.uniformiser_pow(5)
            5^5 + O(5^25)
        """
        return self.uniformizer_pow(n)

    def is_finite(self):
        r"""
        Returns whether this ring is finite, i.e. ``False``.

        INPUT:

        - ``self`` -- a `p`-adic ring

        OUTPUT:

        - boolean -- whether self is finite, i.e., ``False``

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.is_finite()
            False
        """
        return False

    def ext(self, *args, **kwds):
        """
        Constructs an extension of self.  See ``extension`` for more details.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.uniformiser()
            t + O(t^21)
        """
        return self.extension(*args, **kwds)

    def _test_add_bigoh(self, **options):
        r"""
        Perform tests on ``add_bigoh``.

        EXAMPLES::

            sage: K = Qp(3)
            sage: K._test_add_bigoh()

        """
        tester = self._tester(**options)
        for x in tester.some_elements():
            tester.assertEqual(x.add_bigoh(x.precision_absolute()), x)
            from sage.rings.all import infinity
            tester.assertEqual(x.add_bigoh(infinity), x)
            tester.assertEqual(x.add_bigoh(x.precision_absolute()+1), x)

            y = x.add_bigoh(0)
            tester.assertIs(y.parent(), self)
            if self.is_capped_absolute():
                tester.assertEqual(y.precision_absolute(), 0)
                tester.assertEqual(y, self.zero())
            elif self.is_capped_relative():
                tester.assertLessEqual(y.precision_absolute(), 0)
            elif self.is_fixed_mod() or self.is_floating_point():
                tester.assertGreaterEqual((x-y).valuation(), 0)
            else:
                raise NotImplementedError

            # if absprec < 0, then the result is in the fraction field (see #13591)
            y = x.add_bigoh(-1)
            tester.assertIs(y.parent(), self.fraction_field())
            if not self.is_floating_point() and not self.is_fixed_mod():
                tester.assertLessEqual(y.precision_absolute(), -1)

            # make sure that we handle very large values correctly
            absprec = Integer(2)**1000
            tester.assertEqual(x.add_bigoh(absprec), x)

    def _test_residue(self, **options):
        r"""
        Perform some tests on the residue field of this ring.

        EXAMPLES::

            sage: R = Zp(2)
            sage: R._test_residue()

        """
        tester = self._tester(**options)
        tester.assertEqual(self.residue_field().characteristic(), self.residue_characteristic())

        for x in tester.some_elements():
            errors = []
            if x.precision_absolute() <= 0:
                from .precision_error import PrecisionError
                errors.append(PrecisionError)
            if x.valuation() < 0:
                errors.append(ValueError)
            if errors:
                with tester.assertRaises(tuple(errors)):
                    x.residue()
                continue
            y = x.residue()
            # residue() is in `Z/pZ` which is not identical to the residue field `F_p`
            tester.assertEqual(y.parent().cardinality(), self.residue_field().cardinality())
            z = self(y)
            tester.assertGreater((x-z).valuation(), 0)

        for x in self.residue_field().some_elements():
            y = self(x)
            if x.is_zero():
                tester.assertGreater(y.valuation(), 0)
            else:
                tester.assertEqual(y.valuation(), 0)
            z = y.residue()
            tester.assertEqual(x, z)
