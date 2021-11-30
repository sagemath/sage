r"""
Local Generic

Superclass for `p`-adic and power series rings.

AUTHORS:

- David Roe
"""

# *****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from copy import copy
from sage.rings.ring import CommutativeRing
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationRings, CompleteDiscreteValuationFields
from sage.structure.category_object import check_default_category
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

class LocalGeneric(CommutativeRing):
    def __init__(self, base, prec, names, element_class, category=None):
        r"""
        Initialize ``self``.

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

            sage: R = Zp(3, 10,'fixed-mod')
            sage: R.is_finite()
            False
            sage: R.cardinality()
            +Infinity

            sage: Qp(11).is_finite()
            False
            sage: Qp(11).cardinality()
            +Infinity
        """
        self._prec = prec
        self.Element = element_class
        default_category = getattr(self, '_default_category', None)
        if self.is_field():
            category = CompleteDiscreteValuationFields()
        else:
            category = CompleteDiscreteValuationRings()
        category = category.Metric().Complete().Infinite()
        if default_category is not None:
            category = check_default_category(default_category, category)
        Parent.__init__(self, base, names=(names,), normalize=False, category=category)

    def is_capped_relative(self):
        r"""
        Return whether this `p`-adic ring bounds precision in a capped
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
        r"""
        Return whether this `p`-adic ring bounds precision in a
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
        r"""
        Return whether this `p`-adic ring bounds precision in a fixed
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
            5^7
            sage: S = ZpCA(5, 15)
            sage: S.is_fixed_mod()
            False
            sage: S(5^7,absprec=9)
            5^7 + O(5^9)
        """
        return False

    def is_floating_point(self):
        r"""
        Return whether this `p`-adic ring bounds precision in a floating
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

    def is_lattice_prec(self):
        r"""
        Return whether this `p`-adic ring bounds precision using
        a lattice model.

        In lattice precision, relationships between elements
        are stored in a precision object of the parent, which
        allows for optimal precision tracking at the cost of
        increased memory usage and runtime.

        EXAMPLES::

            sage: R = ZpCR(5, 15)
            sage: R.is_lattice_prec()
            False
            sage: x = R(25, 8)
            sage: x - x
            O(5^8)
            sage: S = ZpLC(5, 15)
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See http://trac.sagemath.org/23505 for details.
            sage: S.is_lattice_prec()
            True
            sage: x = S(25, 8)
            sage: x - x
            O(5^30)
        """
        return False

    def is_relaxed(self):
        r"""
        Return whether this `p`-adic ring bounds precision in a relaxed
        fashion.

        In a relaxed ring, elements have mechanisms for computing
        themselves to greater precision.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.is_relaxed()
            False
        """
        return False

    def _latex_(self):
        r"""
        Latex.

        EXAMPLES::

            sage: latex(Zq(27,names='a')) #indirect doctest
            \Bold{Z}_{3^{3}}
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
        - ``label`` -- string (only for lattice precision)

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
            2-adic Unramified Extension Field in a defined by x^6 + x^4 + x^3 + x + 1
            sage: R.<x> = QQ[]
            sage: K.change(modulus = x^2 - x + 2, print_pos=False)
            5-adic Unramified Extension Field in a defined by x^2 - x + 2

        and variable names::

            sage: K.change(names='b')
            5-adic Unramified Extension Field in b defined by x^3 + 3*x + 3

        and precision::

            sage: Kup = K.change(prec=8); Kup
            5-adic Unramified Extension Field in a defined by x^3 + 3*x + 3
            sage: Kup.precision_cap()
            8
            sage: Kup.base_ring()
            5-adic Field with floating precision 8

        If you decrease the precision, the precision of the base stays the same::

            sage: Kdown = K.change(prec=2); Kdown
            5-adic Unramified Extension Field in a defined by x^3 + 3*x + 3
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

        Changing label for lattice precision (the precision lattice is not copied)::

            sage: R = ZpLC(37, (8,11))
            sage: S = R.change(label = "change"); S
            37-adic Ring with lattice-cap precision (label: change)
            sage: S.change(label = "new")
            37-adic Ring with lattice-cap precision (label: new)
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
        functor, ring = self.construction(forbid_frac_field=True)
        functor = copy(functor)
        if 'mode' in kwds and 'show_prec' not in kwds:
            new_type = kwds.get('type', self._prec_type())
            cur_type = self._prec_type()
            cur_mode = self._printer._print_mode()
            cur_show_prec = self._printer._show_prec()
            from .factory import _canonicalize_show_prec
            if cur_show_prec == _canonicalize_show_prec(cur_type, cur_mode):
                kwds['show_prec'] = _canonicalize_show_prec(new_type, kwds['mode'])
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
            # Labels for lattice precision
            if 'label' in kwds:
                functor.extras['label'] = kwds.pop('label')
            elif 'label' in functor.extras and functor.type not in ['lattice-cap','lattice-float']:
                del functor.extras['label']
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
                # This will need to be modified once lattice precision supports extensions
                prec = kwds.pop('prec')
                baseprec = (prec - 1) // self.relative_e() + 1
                if baseprec > self.base_ring().precision_cap():
                    kwds['prec'] = baseprec
                functor.precs = [prec]
            from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
            if 'names' in kwds:
                functor.names = [kwds.pop('names')]
            modulus = None
            if 'modulus' in kwds:
                modulus = kwds.pop('modulus')
                if n is not None and modulus.degree() != n:
                    raise ValueError("modulus must have degree matching q")
            elif q is not None:
                if self.relative_e() == 1:
                    # If q is specified, replace the modulus with one from q.
                    modulus = get_unramified_modulus(q, functor.kwds.get('res_name', functor.names[0] + '0'))
                elif self.relative_f() != 1:
                    raise ValueError("Cannot change q in mixed extensions")
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
                if q is not None and self.relative_f() == 1:
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
        Return the precision cap for this ring.

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
            part of each element.
        """
        return self._prec

    def _precision_cap(self):
        r"""
        Return the precision cap for this ring, in the format
        used by the factory methods to create the ring.

        For most `p`-adic types, this is the same as :meth:`precision_cap`,
        but there is a difference for lattice precision.

        EXAMPLES::

            sage: Zp(17,34)._precision_cap()
            34
        """
        return self._prec

    def is_exact(self):
        r"""
        Return whether this `p`-adic ring is exact, i.e. ``False``.

        EXAMPLES::

            sage: R = Zp(5, 3, 'fixed-mod'); R.is_exact()
            False
        """
        return False

    def residue_characteristic(self):
        r"""
        Return the characteristic of ``self``'s residue field.

        INPUT:

        - ``self`` -- a p-adic ring.

        OUTPUT:

        The characteristic of the residue field.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.residue_characteristic()
            3
        """
        return self.residue_class_field().characteristic()

    def defining_polynomial(self, var='x', exact=False):
        r"""
        Return the defining polynomial of this local ring

        INPUT:

        - ``var`` -- string (default: ``'x'``), the name of the variable

        - ``exact`` -- a boolean (default: ``False``), whether to return the
          underlying exact  defining polynomial rather than the one with coefficients
          in the base ring.

        OUTPUT:

        The defining polynomial of this ring as an extension over its ground ring

        EXAMPLES::

            sage: R = Zp(3, 3, 'fixed-mod')

            sage: R.defining_polynomial().parent()
            Univariate Polynomial Ring in x over 3-adic Ring of fixed modulus 3^3
            sage: R.defining_polynomial('foo')
            foo

            sage: R.defining_polynomial(exact=True).parent()
            Univariate Polynomial Ring in x over Integer Ring
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        if exact:
            from sage.rings.integer_ring import ZZ
            return PolynomialRing(ZZ,var).gen()
        else:
            return PolynomialRing(self,var).gen()

    def ground_ring(self):
        r"""
        Return ``self``.

        Will be overridden by extensions.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        The ground ring of ``self``, i.e., itself.

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
        Return ``self``.

        Will be overridden by extensions.

        INPUT:

        - ``self`` -- a `p`-adic ring

        OUTPUT:

        The ground ring of the tower for ``self``, i.e., itself.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.ground_ring_of_tower()
            5-adic Ring with capped relative precision 20
        """
        return self


    def absolute_degree(self):
        r"""
        Return the degree of this extension over the prime p-adic field/ring.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.absolute_degree()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.absolute_degree()
            2
        """
        return self.absolute_e() * self.absolute_f()

    def relative_degree(self):
        r"""
        Return the degree of this extension over its base field/ring.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.relative_degree()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.relative_degree()
            2
        """
        return self.absolute_degree() // self.base_ring().absolute_degree()

    def degree(self):
        r"""
        Return the degree of this extension.

        Raise an error if the base ring/field is itself an extension.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.degree()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.degree()
            2
        """
        if self.base_ring().absolute_degree() == 1:
            return self.absolute_degree()
        else:
            raise NotImplementedError("For a relative p-adic ring or field you must use relative_degree or absolute_degree as appropriate")


    def absolute_e(self):
        r"""
        Return the absolute ramification index of this ring/field.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.absolute_e()
            1

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.absolute_e()
            2
        """
        # Override this in subclasses (if appropriate)
        if self is self.base_ring():
            return ZZ(1)
        else:
            return self.base_ring().absolute_e()

    def absolute_ramification_index(self):
        r"""
        Return the absolute ramification index of this ring/field.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.absolute_ramification_index()
            1

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.absolute_ramification_index()
            2
        """
        return self.absolute_e()

    def relative_e(self):
        r"""
        Return the ramification index of this extension over its base ring/field.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.relative_e()
            1

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.relative_e()
            2
        """
        return self.absolute_e() // self.base_ring().absolute_e()

    def relative_ramification_index(self):
        r"""
        Return the ramification index of this extension over its base ring/field.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.relative_ramification_index()
            1

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.relative_ramification_index()
            2
        """
        return self.relative_e()

    def e(self):
        r"""
        Return the ramification index of this extension.

        Raise an error if the base ring/field is itself an extension.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.e()
            1

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.e()
            2
        """
        if self.base_ring().absolute_degree() == 1:
            return self.absolute_e()
        else:
            raise NotImplementedError("For a relative p-adic ring or field you must use relative_e or absolute_e as appropriate")

    def ramification_index(self):
        r"""
        Return the ramification index of this extension.

        Raise an error if the base ring/field is itself an extension.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.ramification_index()
            1

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.ramification_index()
            2
        """
        return self.e()


    def absolute_f(self):
        r"""
        Return the degree of the residue field of this ring/field
        over its prime subfield.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.absolute_f()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.absolute_f()
            1
        """
        # Override this in subclasses (if appropriate)
        if self is self.base_ring():
            return ZZ(1)
        else:
            return self.base_ring().absolute_f()

    def absolute_inertia_degree(self):
        r"""
        Return the degree of the residue field of this ring/field
        over its prime subfield.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.absolute_inertia_degree()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.absolute_inertia_degree()
            1
        """
        return self.absolute_f()

    def relative_f(self):
        r"""
        Return the degree of the residual extension over its base ring/field.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.relative_f()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.relative_f()
            1
        """
        return self.absolute_f() // self.base_ring().absolute_f()

    def relative_inertia_degree(self):
        r"""
        Return the degree of the residual extension over its base ring/field.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.relative_inertia_degree()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.relative_inertia_degree()
            1
        """
        return self.relative_f()

    def f(self):
        r"""
        Return the degree of the residual extension.

        Raise an error if the base ring/field is itself an extension.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.f()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.f()
            1
        """
        if self.base_ring().absolute_degree() == 1:
            return self.absolute_f()
        else:
            raise NotImplementedError("For a relative p-adic ring or field you must use relative_f or absolute_f as appropriate")

    def inertia_degree(self):
        r"""
        Return the degree of the residual extension.

        Raise an error if the base ring/field is itself an extension.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.inertia_degree()
            5

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.inertia_degree()
            1
        """
        return self.f()

    def inertia_subring(self):
        r"""
        Return the inertia subring, i.e. ``self``.

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
        Return the maximal unramified subextension.

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
#        Return the trivial extension of self.
#        """
#        raise NotImplementedError

    def uniformiser(self):
        r"""
        Return a uniformiser for ``self``, ie a generator for the unique maximal ideal.

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
        r"""
        Return the `n`th power of the uniformiser of ``self`` (as an element of ``self``).

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.uniformiser_pow(5)
            5^5 + O(5^25)
        """
        return self.uniformizer_pow(n)

    def ext(self, *args, **kwds):
        r"""
        Construct an extension of self.  See :meth:`extension` for more details.

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
            from sage.rings.infinity import infinity
            tester.assertEqual(x.add_bigoh(infinity), x)
            tester.assertEqual(x.add_bigoh(x.precision_absolute()+1), x)

            y = x.add_bigoh(0)
            tester.assertIs(y.parent(), self)
            if self.is_capped_absolute():
                tester.assertEqual(y.precision_absolute(), 0)
                tester.assertEqual(y, self.zero())
            elif self.is_capped_relative() or self.is_lattice_prec():
                tester.assertLessEqual(y.precision_absolute(), 0)
            elif self.is_fixed_mod() or self.is_floating_point():
                tester.assertGreaterEqual((x-y).valuation(), 0)

            # if absprec < 0, then the result is in the fraction field (see #13591)
            y = x.add_bigoh(-1)
            tester.assertIs(y.parent(), self.fraction_field())
            if not self.is_floating_point() and not self.is_fixed_mod():
                tester.assertLessEqual(y.precision_absolute(), -1)

            # make sure that we handle very large values correctly
            if self._prec_type() not in [ 'lattice-float', 'relaxed' ]:   # no cap in these models
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

    def _matrix_flatten_precision(self, M):
        r"""
        Rescale rows and columns of ``M`` so that the minimal
        absolute precision of each row and column is equal to
        the cap.

        This method is useful for increasing the numerical
        stability. It is called by :meth:`_matrix_smith_form`
        and :meth:`_matrix_determinant`

        Only for internal use.

        OUTPUT:

        The lists of valuations by which rows and columns,
        respectively, have been shifted.

        EXAMPLES::

            sage: K = Qp(2, print_mode='digits', prec=10)
            sage: M = matrix(K, 2, 2, [K(1,5),K(2,7),K(3,3),K(5,8)])
            sage: M
            [   ...00001  ...0000010]
            [     ...011 ...00000101]
            sage: K._matrix_flatten_precision(M)
            ([5, 7], [0, -2])
            sage: M
            [   ...0000100000    ...0000010000]
            [   ...0110000000 ...0000010100000]
        """
        parent = M.base_ring()
        cap = parent.precision_cap()
        n = M.nrows()
        m = M.ncols()
        shift_rows = n * [ ZZ(0) ]
        shift_cols = m * [ ZZ(0) ]
        for i in range(n):
            prec = min(M[i,j].precision_absolute() for j in range(m))
            if prec is Infinity or prec == cap:
                continue
            shift_rows[i] = s = cap - prec
            for j in range(m):
                M[i,j] <<= s
        for j in range(m):
            prec = min(M[i,j].precision_absolute() for i in range(n))
            if prec is Infinity or prec == cap:
                continue
            shift_cols[j] = s = cap - prec
            for i in range(n):
                M[i,j] <<= s
        return shift_rows, shift_cols


    def _matrix_smith_form(self, M, transformation, integral, exact):
        r"""
        Return the Smith normal form of the matrix `M`.

        This method gets called by
        :meth:`sage.matrix.matrix2.Matrix.smith_form` to compute the Smith
        normal form over local rings and fields.

        The entries of the Smith normal form are normalized such that non-zero
        entries of the diagonal are powers of the distinguished uniformizer.

        INPUT:

        - ``M`` -- a matrix over this ring

        - ``transformation`` -- a boolean; whether the transformation matrices
          are returned

        - ``integral`` -- a subring of the base ring or ``True``; the entries
          of the transformation matrices are in this ring.  If ``True``, the
          entries are in the ring of integers of the base ring.

        - ``exact`` -- boolean.  If ``True``, the diagonal smith form will
          be exact, or raise a ``PrecisionError`` if this is not possible.
          If ``False``, the diagonal entries will be inexact, but the
          transformation matrices will be exact.

        EXAMPLES::

            sage: A = Zp(5, prec=10, print_mode="digits")
            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

            sage: S, L, R = M.smith_form()  # indirect doctest
            sage: S
            [ ...1     0]
            [    0 ...10]
            sage: L
            [...222222223          ...]
            [...444444444         ...2]
            sage: R
            [...0000000001 ...2222222214]
            [            0 ...0000000001]

        If not needed, it is possible to avoid the computation of
        the transformations matrices `L` and `R`::

            sage: M.smith_form(transformation=False)  # indirect doctest
            [ ...1     0]
            [    0 ...10]

        This method works for rectangular matrices as well::

            sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
            sage: S, L, R = M.smith_form()  # indirect doctest
            sage: S
            [ ...1     0]
            [    0 ...10]
            [    0     0]
            sage: L
            [...222222223          ...          ...]
            [...444444444         ...2          ...]
            [...444444443         ...1         ...1]
            sage: R
            [...0000000001 ...2222222214]
            [            0 ...0000000001]

        If some of the elementary divisors have valuation larger than the
        minimum precision of any entry in the matrix, then they are
        reported as an inexact zero::

            sage: A = ZpCA(5, prec=10)
            sage: M = matrix(A, 2, 2, [5, 5, 5, 5])
            sage: M.smith_form(transformation=False, exact=False)  # indirect doctest
            [5 + O(5^10)     O(5^10)]
            [    O(5^10)     O(5^10)]

        However, an error is raised if the precision on the entries is
        not enough to determine which column to use as a pivot at some point::

            sage: M = matrix(A, 2, 2, [A(0,5), A(5^6,10), A(0,8), A(5^7,10)]); M
            [       O(5^5) 5^6 + O(5^10)]
            [       O(5^8) 5^7 + O(5^10)]
            sage: M.smith_form(transformation=False, exact=False)  # indirect doctest
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to compute Smith normal form

        TESTS::

            sage: A = ZpCR(5, prec=10)
            sage: M = zero_matrix(A, 2)
            sage: M.smith_form(transformation=False)  # indirect doctest
            [0 0]
            [0 0]

            sage: M = matrix(2, 2, [ A(0,10), 0, 0, 0] )
            sage: M.smith_form(transformation=False)  # indirect doctest
            Traceback (most recent call last):
            ...
            PrecisionError: some elementary divisors indistinguishable from zero (try exact=False)
            sage: M.smith_form(transformation=False, exact=False)  # indirect doctest
            [O(5^10) O(5^10)]
            [O(5^10) O(5^10)]
        """
        from sage.rings.infinity import infinity
        from .precision_error import PrecisionError
        from copy import copy
        n = M.nrows()
        m = M.ncols()
        if m > n:
            ## It's easier below if we can always deal with precision on left.
            if transformation:
                d, u, v = self._matrix_smith_form(M.transpose(), True, integral, exact)
                return d.transpose(), v.transpose(), u.transpose()
            else:
                return self._matrix_smith_form(M.transpose(), False, integral, exact).transpose()
        smith = M.parent()(0)
        S = copy(M)
        Z = self.integer_ring()
        if integral is None or integral is self or integral is (not self.is_field()):
            integral = not self.is_field()
            R = self
        elif integral is True or integral is Z:
            # This is a field, but we want the integral smith form
            # The diagonal matrix may not be integral, but the transformations should be
            R = Z
            integral = True
        elif integral is False or integral is self.fraction_field():
            # This is a ring, but we want the field smith form
            # The diagonal matrix should be over this ring, but the transformations should not
            R = self.fraction_field()
            integral = False
        else:
            raise NotImplementedError("Smith normal form over this subring")
        ## the difference between ball_prec and inexact_ring is just for lattice precision.
        ball_prec = R._prec_type() in ['capped-rel','capped-abs']
        inexact_ring = R._prec_type() not in ['fixed-mod','floating-point']

        if not integral:
            shift_rows, shift_cols = self._matrix_flatten_precision(S)

        precS = min(x.precision_absolute() for x in S.list())
        if transformation:
            from sage.matrix.special import identity_matrix
            left = identity_matrix(R,n)
            right = identity_matrix(R,m)

        if ball_prec and precS is infinity: # capped-rel and M = 0 exactly
            return (smith, left, right) if transformation else smith

        val = -infinity
        for piv in range(m): # m <= n
            curval = infinity
            pivi = pivj = piv
            # allzero tracks whether every possible pivot is zero.
            # if so, we can stop.  allzero is also used in detecting some
            # precision problems: if we can't determine what pivot has
            # the smallest valuation, or if exact=True and some elementary
            # divisor is zero modulo the working precision
            allzero = True
            # allexact is tracked because there is one case where we can correctly
            # deduce the exact smith form even with some elementary divisors zero:
            # if the bottom right block consists entirely of exact zeros.
            allexact = True
            for i in range(piv,n):
                for j in range(piv,m):
                    Sij = S[i,j]
                    v = Sij.valuation()
                    allzero = allzero and Sij.is_zero()
                    if exact: # we only care in this case
                        allexact = allexact and Sij.precision_absolute() is infinity
                    if v < curval:
                        pivi = i
                        pivj = j
                        curval = v
                        if v == val:
                            break
                else:
                    continue
                break
            val = curval

            if inexact_ring and not allzero and val >= precS:
                if ball_prec:
                    raise PrecisionError("not enough precision to compute Smith normal form")
                precS = min([ S[i,j].precision_absolute() for i in range(piv,n) for j in range(piv,m) ])
                if val >= precS:
                    raise PrecisionError("not enough precision to compute Smith normal form")

            if allzero:
                if exact:
                    if allexact:
                        # We need to finish checking allexact since we broke out of the loop early
                        for i in range(i,n):
                            for j in range(piv,m):
                                allexact = allexact and S[i,j].precision_absolute() is infinity
                                if not allexact:
                                    break
                            else:
                                continue
                            break
                    if not allexact:
                        raise PrecisionError("some elementary divisors indistinguishable from zero (try exact=False)")
                break

            # We swap the lowest valuation pivot into position
            S.swap_rows(pivi,piv)
            S.swap_columns(pivj,piv)
            if transformation:
                left.swap_rows(pivi,piv)
                right.swap_columns(pivj,piv)

            # ... and clear out this row and column.  Note that we
            # will deal with precision later, thus the call to lift_to_precision
            smith[piv,piv] = self(1) << val
            inv = (S[piv,piv] >> val).inverse_of_unit()
            if ball_prec:
                inv = inv.lift_to_precision()
            for i in range(piv+1,n):
                scalar = -inv * Z(S[i,piv] >> val)
                if ball_prec:
                    scalar = scalar.lift_to_precision()
                S.add_multiple_of_row(i,piv,scalar,piv+1)
                if transformation:
                    left.add_multiple_of_row(i,piv,scalar)
            if transformation:
                left.rescale_row(piv,inv)
                for j in range(piv+1,m):
                    scalar = -inv * Z(S[piv,j] >> val)
                    if ball_prec:
                        scalar = scalar.lift_to_precision()
                    right.add_multiple_of_column(j,piv,scalar)
        else:
            # We use piv as an upper bound on a range below, and need to set it correctly
            # in the case that we didn't break out of the loop
            piv = m
        # We update the precision on left
        # The bigoh measures the effect of multiplying by row operations
        # on the left in order to clear out the digits in the smith form
        # with valuation at least precS
        if ball_prec and exact and transformation:
            for j in range(n):
                delta = min(left[i,j].valuation() - smith[i,i].valuation() for i in range(piv))
                if delta is not infinity:
                    for i in range(n):
                        left[i,j] = left[i,j].add_bigoh(precS + delta)
        ## Otherwise, we update the precision on smith
        if ball_prec and not exact:
            smith = smith.apply_map(lambda x: x.add_bigoh(precS))
        ## We now have to adjust the elementary divisors (and precision) in the non-integral case
        if not integral:
            for i in range(piv):
                v = smith[i,i].valuation()
                if transformation:
                    for j in range(n):
                        left[i,j] >>= v
                if exact:
                    smith[i,i] = self(1)
                else:
                    for j in range(n):
                        smith[i,j] = smith[i,j] >> v
            if transformation:
                for i in range(n):
                    for j in range(n):
                        left[i,j] <<= shift_rows[j]
                for i in range(m):
                    for j in range(m):
                        right[i,j] <<= shift_cols[i]
        if transformation:
            return smith, left, right
        else:
            return smith

    def _test_matrix_smith(self, **options):
        r"""
        Test that :meth:`_matrix_smith_form` works correctly.

        EXAMPLES::

            sage: ZpCA(5, 15)._test_matrix_smith()

        """
        tester = self._tester(**options)
        tester.assertEqual(self.residue_field().characteristic(), self.residue_characteristic())

        from itertools import chain
        from sage.all import MatrixSpace
        from .precision_error import PrecisionError
        matrices = chain(*[MatrixSpace(self, n, m).some_elements() for n in (1,3,7) for m in (1,4,7)])
        for M in tester.some_elements(matrices):
            bases = [self]
            if self is not self.integer_ring():
                bases.append(self.integer_ring())
            for base in bases:
                try:
                    S,U,V = M.smith_form(integral=base)
                except PrecisionError:
                    continue

                if self.is_exact() or self._prec_type() not in ['fixed-mod','floating-point']:
                    tester.assertEqual(U*M*V, S)

                tester.assertEqual(U.nrows(), U.ncols())
                tester.assertEqual(U.base_ring(), base)

                tester.assertEqual(V.nrows(), V.ncols())
                tester.assertEqual(V.base_ring(), base)

                for d in S.diagonal():
                    if not d.is_zero():
                        tester.assertTrue(d.unit_part().is_one())

                for (d,dd) in zip(S.diagonal(), S.diagonal()[1:]):
                    tester.assertTrue(d.divides(dd))

    def _matrix_determinant(self, M):
        r"""
        Return the determinant of the matrix `M`.

        This method gets called by
        :meth:`sage.matrix.matrix2.Matrix.determinant`.

        INPUT:

        - ``M`` -- a matrix over this ring

        ALGORITHM:

        We flatten the absolute precision in order to increase
        the numerical stability.

        We row-echelonize the matrix by always choosing the
        pivot of smallest valuation and allowing permutations
        of columns.

        Then we compute separately the value of the determinant
        (as the product of the diagonal entries of the row-echelon
        form) and a bound on the precision on it.

        EXAMPLES::

            sage: R = Qp(5,10)
            sage: M = matrix(R, 2, 2, [1, 6, 2, 7])
            sage: M.determinant()  # indirect doctest
            4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)

            sage: (5*M).determinant()  # indirect doctest
            4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12)

        Sometimes, we gain precision on the determinant::

            sage: M = matrix(R, 3, 3,
            ....:             [R(16820,7), R(73642,7), R( 3281,7),
            ....:              R(67830,7), R(63768,7), R(76424,7),
            ....:              R(37790,7), R(38784,7), R(69287,7)])
            sage: M.determinant()  # indirect doctest
            4*5^5 + 4*5^6 + 3*5^7 + 2*5^8 + O(5^9)

        TESTS:

        We check the stability of our algorithm::

            sage: for dim in range(3,10):
            ....:     M = matrix(dim, dim, [ R(1) for _ in range(dim^2) ])
            ....:     print(M.determinant())
            O(5^20)
            O(5^30)
            O(5^40)
            O(5^50)
            O(5^60)
            O(5^70)
            O(5^80)

            sage: A = random_matrix(Qp(5),4)
            sage: B = random_matrix(Qp(5),4)
            sage: (A*B).det() == A.det()*B.det()
            True
            sage: A.change_ring(QQ).det() == A.det()
            True
            sage: matrix(Qp(37),[0]).determinant()
            0
            sage: matrix(Qp(37),[O(37)]).determinant()
            O(37)
        """
        n = M.nrows()

        # For 2x2 matrices, we use the formula
        if n == 2:
            return M[0,0]*M[1,1] - M[0,1]*M[1,0]

        R = M.base_ring()
        track_precision = R._prec_type() in ['capped-rel','capped-abs']

        S = copy(M)
        shift_rows, shift_cols = self._matrix_flatten_precision(S)
        shift = sum(shift_rows) + sum(shift_cols)
        det = R(1)

        sign = 1
        valdet = 0
        val = -Infinity
        for piv in range(n):
            pivi = pivj = piv
            curval = S[pivi, pivj].valuation()
            for i in range(piv,n):
                for j in range(piv,n):
                    v = S[i,j].valuation()
                    if v < curval:
                        pivi = i
                        pivj = j
                        curval = v
                        if v == val:
                            break
                else:
                    continue
                break
            val = curval
            if S[pivi,pivj] == 0:
                if track_precision:
                    return R(0, valdet + (n-piv)*val - shift)
                else:
                    return R(0)

            valdet += val
            S.swap_rows(pivi,piv)
            if pivi > piv:
                sign = -sign
            S.swap_columns(pivj,piv)
            if pivj > piv:
                sign = -sign

            det *= S[piv,piv]
            inv = ~(S[piv,piv] >> val)
            for i in range(piv+1,n):
                scalar = -inv * (S[i,piv] >> val)
                if track_precision:
                    scalar = scalar.lift_to_precision()
                S.add_multiple_of_row(i,piv,scalar)

        if track_precision:
            relprec = +Infinity
            relprec_neg = 0
            for i in range(n):
                prec = Infinity
                for j in range(n):
                    prec = min(prec, S[i,j].precision_absolute())
                prec -= S[i,i].valuation()
                if prec < relprec:
                    relprec = prec
                if prec < 0:
                    relprec_neg += prec
            if relprec_neg < 0:
                relprec = relprec_neg
            det = (sign*det).add_bigoh(valdet+relprec)
        else:
            det = sign*det
        return det >> shift
