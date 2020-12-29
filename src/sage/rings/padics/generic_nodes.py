"""
`p`-Adic Generic Nodes

This file contains a bunch of intermediate classes for the `p`-adic
parents, allowing a function to be implemented at the right level of
generality.

AUTHORS:

- David Roe
"""

# ****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.padics.local_generic import LocalGeneric
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.ring import EuclideanDomain, Field
from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import infinity, SignError
from .lattice_precision import PrecisionLattice, PrecisionModule
from sage.rings.padics.precision_error import PrecisionError
from .padic_lattice_element import pAdicLatticeElement, pAdicLatticeCapElement, pAdicLatticeFloatElement


class CappedAbsoluteGeneric(LocalGeneric):
    def is_capped_absolute(self):
        """
        Returns whether this `p`-adic ring bounds precision in a
        capped absolute fashion.

        The absolute precision of an element is the power of `p` modulo
        which that element is defined.  In a capped absolute ring, the
        absolute precision of elements are bounded by a constant
        depending on the ring.

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
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: ZpCA(5)._prec_type()
            'capped-abs'
        """
        return 'capped-abs'

class CappedRelativeGeneric(LocalGeneric):
    def is_capped_relative(self):
        """
        Returns whether this `p`-adic ring bounds precision in a capped
        relative fashion.

        The relative precision of an element is the power of p modulo
        which the unit part of that element is defined.  In a capped
        relative ring, the relative precision of elements are bounded
        by a constant depending on the ring.

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
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: Zp(5)._prec_type()
            'capped-rel'
        """
        return 'capped-rel'

class FixedModGeneric(LocalGeneric):
    def is_fixed_mod(self):
        """
        Returns whether this `p`-adic ring bounds precision in a fixed
        modulus fashion.

        The absolute precision of an element is the power of p modulo
        which that element is defined.  In a fixed modulus ring, the
        absolute precision of every element is defined to be the
        precision cap of the parent.  This means that some operations,
        such as division by `p`, don't return a well defined answer.

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
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: ZpFM(5)._prec_type()
            'fixed-mod'
        """
        return 'fixed-mod'

class FloatingPointGeneric(LocalGeneric):
    def is_floating_point(self):
        """
        Returns whether this `p`-adic ring uses a floating point precision model.

        Elements in the floating point model are stored by giving a
        valuation and a unit part.  Arithmetic is done where the unit
        part is truncated modulo a fixed power of the uniformizer,
        stored in the precision cap of the parent.

        EXAMPLES::

            sage: R = ZpFP(5,15)
            sage: R.is_floating_point()
            True
            sage: R(5^7,absprec=9)
            5^7
            sage: S = ZpCR(5,15)
            sage: S.is_floating_point()
            False
            sage: S(5^7,absprec=9)
            5^7 + O(5^9)
        """
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: ZpFP(5)._prec_type()
            'floating-point'
        """
        return 'floating-point'

    def _test_distributivity(self, **options):
        r"""
        Test the distributivity of `*` on `+` on (not necessarily
        all) elements of this set.

        p-adic floating point rings only satisfy distributivity
        up to a precision that depends on the elements.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

        EXAMPLES:

        By default, this method runs the tests only on the
        elements returned by ``self.some_elements()``::

            sage: R = ZpFP(5,3)
            sage: R.some_elements()
            [0, 1, 5, 1 + 3*5 + 3*5^2, 5 + 4*5^2 + 4*5^3]
            sage: R._test_distributivity()

        However, the elements tested can be customized with the
        ``elements`` keyword argument::

            sage: R._test_distributivity(elements=[R(0),~R(0),R(42)])

        See the documentation for :class:`TestSuite` for more information.
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        from sage.misc.misc import some_tuples
        for x,y,z in some_tuples(S, 3, tester._max_runs):
            yz_prec = min(y.precision_absolute(), z.precision_absolute())
            yz_val = (y + z).valuation()
            try:
                prec = min(x.valuation() + yz_val + min(x.precision_relative(), yz_prec - yz_val),
                           x.valuation() + y.valuation() + (x * y).precision_relative(),
                           x.valuation() + z.valuation() + (x * z).precision_relative())
            except SignError:
                pass
            else:
                if prec > -infinity:
                    # only check left distributivity, since multiplication commutative
                    tester.assertTrue((x * (y + z)).is_equal_to((x * y) + (x * z),prec))

    def _test_additive_associativity(self, **options):
        r"""
        Test associativity for (not necessarily all) elements of this
        additive semigroup.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

        EXAMPLES:

        By default, this method tests only the elements returned by
        ``self.some_elements()``::

            sage: R = QpFP(7,3)
            sage: R._test_additive_associativity()

        However, the elements tested can be customized with the
        ``elements`` keyword argument::

            sage: R._test_additive_associativity(elements = [R(0), ~R(0), R(42)])

        See the documentation for :class:`TestSuite` for more information.
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        from sage.misc.misc import some_tuples
        for x,y,z in some_tuples(S, 3, tester._max_runs):
            tester.assertTrue(((x + y) + z).is_equal_to(x + (y + z), min(x.precision_absolute(), y.precision_absolute(), z.precision_absolute())))

class FloatingPointRingGeneric(FloatingPointGeneric):
    pass
class FloatingPointFieldGeneric(FloatingPointGeneric):#, sage.rings.ring.Field):
    pass
class CappedRelativeRingGeneric(CappedRelativeGeneric):
    pass
class CappedRelativeFieldGeneric(CappedRelativeGeneric):#, sage.rings.ring.Field):
    pass

class pAdicLatticeGeneric(pAdicGeneric):
    r"""
    An implementation of the `p`-adic rationals with lattice precision.

    INPUT:

    - `p` -- the underlying prime number

    - ``prec`` -- the precision

    - ``subtype`` -- either ``"cap"`` or ``"float"``,
      specifying the precision model used for tracking precision

    - ``label`` -- a string or ``None`` (default: ``None``)

    TESTS::

        sage: R = ZpLC(17)   # indirect doctest
        doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/23505 for details.
        sage: R._prec_type()
        'lattice-cap'

        sage: R = ZpLF(17)   # indirect doctest
        sage: R._prec_type()
        'lattice-float'

        sage: R = QpLC(17)   # indirect doctest
        sage: R._prec_type()
        'lattice-cap'

        sage: R = QpLF(17)   # indirect doctest
        sage: R._prec_type()
        'lattice-float'
    """
    def __init__(self, p, prec, print_mode, names, label=None):
        """
        Initialization.

        TESTS::

            sage: R = ZpLC(17)   # indirect doctest
            sage: R._prec_type()
            'lattice-cap'
            sage: R._subtype
            'cap'

            sage: R = ZpLF(17)   # indirect doctest
            sage: R._prec_type()
            'lattice-float'
            sage: R._subtype
            'float'
        """
        from sage.rings.padics.lattice_precision import pRational
        self._approx_zero = pRational(p, 0)
        self._approx_one = pRational(p, 1)
        self._approx_minusone = pRational(p, -1)
        if label is None:
            self._label = None
        else:
            self._label = str(label)
        # We do not use the standard attribute element_class
        # because we need to be careful with precision
        # Instead we implement _element_constructor_ (cf below)
        if self._subtype == 'cap':
            (self._prec_cap_relative, self._prec_cap_absolute) = prec
            self._zero_cap = None
            self._precision = PrecisionLattice(p, label)
            element_class = pAdicLatticeCapElement
        elif self._subtype == 'float':
            self._prec_cap_relative = prec
            self._prec_cap_absolute = infinity
            self._zero_cap = prec
            self._precision = PrecisionModule(p, label, prec)
            element_class = pAdicLatticeFloatElement
        else:
            raise ValueError("subtype must be either 'cap' or 'float'")
        self._element_class = self.__make_element_class__(element_class)
        pAdicGeneric.__init__(self, self, p, prec, print_mode, names, None)

    def _prec_type(self):
        """
        Return the precision handling type.

        EXAMPLES::

            sage: ZpLC(5)._prec_type()
            'lattice-cap'
        """
        return 'lattice-' + self._subtype

    def is_lattice_prec(self):
        """
        Returns whether this `p`-adic ring bounds precision using
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
            sage: S.is_lattice_prec()
            True
            sage: x = S(25, 8)
            sage: x - x
            O(5^30)
        """
        return True

    def precision_cap(self):
        """
        Return the relative precision cap for this ring if it is finite.
        Otherwise return the absolute precision cap.

        EXAMPLES::

            sage: R = ZpLC(3)
            sage: R.precision_cap()
            20
            sage: R.precision_cap_relative()
            20

            sage: R = ZpLC(3, prec=(infinity,20))
            sage: R.precision_cap()
            20
            sage: R.precision_cap_relative()
            +Infinity
            sage: R.precision_cap_absolute()
            20

        .. SEEALSO::

            :meth:`precision_cap_relative`, :meth:`precision_cap_absolute`
        """
        if self._prec_cap_relative is not infinity:
            return self._prec_cap_relative
        else:
            return self._prec_cap_absolute

    def _precision_cap(self):
        """
        Return the pair of precisions (for ``lattice-cap``)
        or the relative precision cap (for ``lattice-float``).

        EXAMPLES::

            sage: R = ZpLC(11, (27,37))
            sage: R._precision_cap()
            (27, 37)
            sage: R = ZpLF(11, 14)
            sage: R._precision_cap()
            14
        """
        if self._subtype == 'cap':
            return (self._prec_cap_relative, self._prec_cap_absolute)
        else:
            return self._prec_cap_relative

    def precision_cap_relative(self):
        """
        Return the relative precision cap for this ring.

        EXAMPLES::

            sage: R = ZpLC(3)
            sage: R.precision_cap_relative()
            20

            sage: R = ZpLC(3, prec=(infinity,20))
            sage: R.precision_cap_relative()
            +Infinity

        .. SEEALSO::

            :meth:`precision_cap`, :meth:`precision_cap_absolute`
        """
        return self._prec_cap_relative

    def precision_cap_absolute(self):
        """
        Return the absolute precision cap for this ring.

        EXAMPLES::

            sage: R = ZpLC(3)
            sage: R.precision_cap_absolute()
            40

            sage: R = ZpLC(3, prec=(infinity,20))
            sage: R.precision_cap_absolute()
            20

        .. SEEALSO::

            :meth:`precision_cap`, :meth:`precision_cap_relative`
        """
        return self._prec_cap_absolute

    def precision(self):
        """
        Return the lattice precision object attached to this parent.

        EXAMPLES::

            sage: R = ZpLC(5, label='precision')
            sage: R.precision()
            Precision lattice on 0 objects (label: precision)

            sage: x = R(1, 10); y = R(1, 5)
            sage: R.precision()
            Precision lattice on 2 objects (label: precision)

        .. SEEALSO::

            :class:`sage.rings.padics.lattice_precision.PrecisionLattice`
        """
        return self._precision

    def label(self):
        """
        Return the label of this parent.

        NOTE:

        Labels can be used to distinguish between parents with
        the same defining data.

        They are useful in the lattice precision framework in order
        to limit the size of the lattice modeling the precision (which
        is roughly the number of elements having this parent).

        Elements of a parent with some label do not coerce to a parent
        with a different label. However conversions are allowed.

        EXAMPLES::

            sage: R = ZpLC(5)
            sage: R.label()  # no label by default

            sage: R = ZpLC(5, label='mylabel')
            sage: R.label()
            'mylabel'

        Labels are typically useful to isolate computations.
        For example, assume that we first want to do some calculations
        with matrices::

            sage: R = ZpLC(5, label='matrices')
            sage: M = random_matrix(R, 4, 4)
            sage: d = M.determinant()

        Now, if we want to do another unrelated computation, we can
        use a different label::

            sage: R = ZpLC(5, label='polynomials')
            sage: S.<x> = PolynomialRing(R)
            sage: P = (x-1)*(x-2)*(x-3)*(x-4)*(x-5)

        Without labels, the software would have modeled the
        precision on the matrices and on the polynomials using the same
        lattice (manipulating a lattice of higher
        dimension can have a significant impact on performance).
        """
        return self._label

    def _element_constructor_(self, x, prec=None):
        """
        Create an element of this parent.

        INPUT:

        - ``x``: the datum from which the element is created

        - ``prec`` -- an integer or ``None`` (the default); the
          absolute precision of the created element

        NOTE:

        This function tries to be sharp on precision as much as
        possible.
        For instance, if the datum ``x`` is itself an element of the
        same parent, the software remembers that the created element
        is actually equal to ``x`` (at infinite precision)::

            sage: R = ZpLC(2, prec=(infinity,50))
            sage: x = R(1, 10); x
            1 + O(2^10)
            sage: y = R(x)   # indirect doctest
            sage: y
            1 + O(2^10)
            sage: x - y
            O(2^50)

        TESTS::

            sage: R(x, prec=5)
            1 + O(2^5)
        """
        # We first try the _copy method which is sharp on precision
        try:
            if prec is None:
                return x._copy(parent=self)
            elif x.parent() is self:
                return x.add_bigoh(prec)
            else:
                return x._copy(parent=self).add_bigoh(prec)
        except (TypeError, ValueError, AttributeError):
            pass
        return self._element_class(self, x, prec)

    def convert_multiple(self, *elts):
        """
        Convert a list of elements to this parent.

        NOTE:

        This function tries to be sharp on precision as much as
        possible.
        In particular, if the precision of the input elements are
        handled by a lattice, diffused digits of precision are
        preserved during the conversion.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: x = R(1, 10); y = R(1, 5)
            sage: x,y = x+y, x-y

        Remark that the pair `(x,y)` has diffused digits of precision::

            sage: x
            2 + O(2^5)
            sage: y
            O(2^5)
            sage: x + y
            2 + O(2^11)

            sage: R.precision().diffused_digits([x,y])
            6

        As a consequence, if we convert ``x`` and ``y`` separately, we
        loose some precision::

            sage: R2 = ZpLC(2, label='copy')
            sage: x2 = R2(x); y2 = R2(y)
            sage: x2
            2 + O(2^5)
            sage: y2
            O(2^5)
            sage: x2 + y2
            2 + O(2^5)

            sage: R2.precision().diffused_digits([x2,y2])
            0

        On the other hand, this issue disappears when we use multiple
        conversion::

            sage: x2,y2 = R2.convert_multiple(x,y)
            sage: x2 + y2
            2 + O(2^11)

            sage: R2.precision().diffused_digits([x2,y2])
            6
        """
        p = self.prime()

        # We sort elements by precision lattice
        elt_by_prec = { }
        elt_other = [ ]
        indices = { }
        for i in range(len(elts)):
            x = elts[i]
            idx = id(x)
            if idx in indices:
                indices[idx].append(i)
            else:
                indices[idx] = [i]
            if isinstance(x, pAdicLatticeElement):
                prec = x.parent().precision()
                if prec.prime() != p:
                    raise TypeError("conversion between different p-adic rings not supported")
                if prec in elt_by_prec:
                    elt_by_prec[prec].append(x)
                else:
                    elt_by_prec[prec] = [x]
            else:
                elt_other.append(x)

        # We create the elements
        ans = len(elts)*[None]
        selfprec = self._precision
        # First the elements with precision lattice
        for (prec, L) in elt_by_prec.items():
            if prec is selfprec:
                # Here, we use the _copy method in order
                # to be sharp on precision
                for x in L:
                    y = x._copy(parent=self)
                    for i in indices[id(x)]:
                        ans[i] = y
            else:
                try:
                    lattice = prec.precision_lattice(L)
                except PrecisionError:
                    raise NotImplementedError("multiple conversion of a set of variables for which the module precision is not a lattice is not implemented yet")
                for j in range(len(L)):
                    x = L[j]; dx = [ ]
                    for i in range(j):
                        dx.append([L[i], lattice[i,j]])
                    prec = lattice[j,j].valuation(p)
                    y = self._element_class(self, x.value(), prec, dx=dx, dx_mode='values', check=False, reduce=False)
                    for i in indices[id(x)]:
                        ans[i] = y
                    L[j] = y
        # Now the other elements
        for x in elt_other:
            y = self._element_class(self, x)
            for i in indices[id(x)]:
                ans[i] = y

        # We return the created elements
        return ans

def is_pAdicRing(R):
    """
    Returns ``True`` if and only if ``R`` is a `p`-adic ring (not a
    field).

    EXAMPLES::

        sage: is_pAdicRing(Zp(5))
        True
        sage: is_pAdicRing(RR)
        False
    """
    return isinstance(R, pAdicRingGeneric)

class pAdicRingGeneric(pAdicGeneric, EuclideanDomain):
    def is_field(self, proof = True):
        """
        Returns whether this ring is actually a field, ie ``False``.

        EXAMPLES::

            sage: Zp(5).is_field()
            False
        """
        return False


    def krull_dimension(self):
        r"""
        Returns the Krull dimension of self, i.e. 1

        INPUT:

        - self -- a `p`-adic ring

        OUTPUT:

        - the Krull dimension of self.  Since self is a `p`-adic ring,
          this is 1.

        EXAMPLES::

            sage: Zp(5).krull_dimension()
            1
        """
        return 1

    def _xgcd_univariate_polynomial(self, f, g):
        """
        Extended gcd for univariate polynomial rings over ``self``.

        Should not be called directly. Use ``f.xgcd(g)`` instead.

        INPUT:

         - ``f``, ``g`` - the polynomials of which to take the xgcd

        OUTPUT:

         - A tuple (a, b, c) which satisfies ``a = b*f + c*g``. There
           is no guarantee that a, b, and c are minimal.

        .. WARNING::

            The computations are performed using the standard Euclidean
            algorithm which might produce mathematically incorrect results in
            some cases. See :trac:`13439`.

        EXAMPLES::

            sage: R.<x> = Zp(3,3)[]
            sage: f = x + 1
            sage: f.xgcd(f^2)
            ((1 + O(3^3))*x + 1 + O(3^3), 1 + O(3^3), 0)

        We check that :trac:`13439` has been fixed::

            sage: R.<x> = Zp(3,3)[]
            sage: f = 3*x + 7
            sage: g = 5*x + 9
            sage: f.xgcd(f*g)
            ((3 + O(3^4))*x + 1 + 2*3 + O(3^3), 1 + O(3^3), 0)

            sage: R.<x> = Zp(3)[]
            sage: f = 357555295953*x + 257392844
            sage: g = 225227399*x - 511940255230575
            sage: f.xgcd(f*g)
            ((3^9 + O(3^29))*x + 2 + 2*3 + 3^2 + 2*3^5 + 3^6 + 3^7
             + 3^8 + 3^10 + 3^11 + 2*3^13 + 3^14 + 3^16 + 2*3^19 +
            O(3^20), 1 + 2*3^2 + 3^4 + 2*3^5 + 3^6 + 3^7 +
             2*3^8 + 2*3^10 + 2*3^12 + 3^13 + 3^14 + 3^15 + 2*3^17
             + 3^18 + O(3^20), 0)

        We check low precision computations::

            sage: R.<x> = Zp(3,1)[]
            sage: h = 3*x + 7
            sage: i = 4*x + 9
            sage: h.xgcd(h*i)
            ((3 + O(3^2))*x + 1 + O(3), 1 + O(3), 0)
        """
        from sage.misc.stopgap import stopgap
        stopgap("Extended gcd computations over p-adic fields are performed using the standard Euclidean algorithm which might produce mathematically incorrect results in some cases.", 13439)

        base_ring = f.base_ring()
        fracfield = base_ring.fraction_field()
        f_field = f.change_ring(fracfield)
        g_field = g.change_ring(fracfield)
        xgcd = fracfield._xgcd_univariate_polynomial(f_field,g_field)
        lcm = base_ring(1)
        for f in xgcd:
            for i in f:
                lcm = (i.denominator()).lcm(lcm)
        returnlst = []
        for f in xgcd:
            f *= lcm
            returnlst.append(f.change_ring(base_ring))
        return tuple(returnlst)

    def _gcd_univariate_polynomial(self, f, g):
        """
        gcd for univariate polynomial rings over ``self``.

        INPUT:

         - ``f``, ``g`` - the polynomials of which to take the gcd

        OUTPUT: A polynomial

        EXAMPLES::

            sage: R.<a> = Zq(27)
            sage: K.<x> = R[]
            sage: h = 3*x + a
            sage: i = 4*x + 2
            sage: h.gcd(h*i)
            (3 + O(3^21))*x + a + O(3^20)
        """
        return self._xgcd_univariate_polynomial(f, g)[0]


def is_pAdicField(R):
    """
    Returns ``True`` if and only if ``R`` is a `p`-adic field.

    EXAMPLES::

        sage: is_pAdicField(Zp(17))
        False
        sage: is_pAdicField(Qp(17))
        True
    """
    return isinstance(R, pAdicFieldGeneric)


class pAdicFieldGeneric(pAdicGeneric, Field):
    pass

    #def class_field(self, group=None, map=None, generators=None):
    #    raise NotImplementedError

    #def composite(self, subfield1, subfield2):
    #    raise NotImplementedError

    #def norm_equation(self):
    #    raise NotImplementedError

    #def norm_group(self):
    #    raise NotImplementedError

    #def norm_group_discriminant(self, group=None, map=None, generators=None):
    #    raise NotImplementedError

    #def number_of_extensions(self, degree, discriminant=None, e=None, f=None):
    #    raise NotImplementedError

    #def list_of_extensions(self, degree, discriminant=None, e=None, f=None):
    #    raise NotImplementedError

    #def subfield(self, list):
    #    raise NotImplementedError

    #def subfield_lattice(self):
    #    raise NotImplementedError

    #def subfields_of_degree(self, n):
    #    raise NotImplementedError

class pAdicFixedModRingGeneric(pAdicRingGeneric, FixedModGeneric):
    pass
class pAdicCappedAbsoluteRingGeneric(pAdicRingGeneric, CappedAbsoluteGeneric):
    pass
class pAdicCappedRelativeRingGeneric(pAdicRingGeneric, CappedRelativeRingGeneric):
    pass
class pAdicCappedRelativeFieldGeneric(pAdicFieldGeneric, CappedRelativeFieldGeneric):
    pass
class pAdicFloatingPointRingGeneric(pAdicRingGeneric, FloatingPointRingGeneric):
    pass
class pAdicFloatingPointFieldGeneric(pAdicFieldGeneric, FloatingPointFieldGeneric):
    pass

class pAdicRingBaseGeneric(pAdicBaseGeneric, pAdicRingGeneric):
    def construction(self, forbid_frac_field=False):
        """
        Returns the functorial construction of self, namely,
        completion of the rational numbers with respect a given prime.

        Also preserves other information that makes this field unique
        (e.g. precision, rounding, print mode).

        INPUT:

        - ``forbid_frac_field`` -- ignored, for compatibility with other p-adic types.

        EXAMPLES::

            sage: K = Zp(17, 8, print_mode='val-unit', print_sep='&')
            sage: c, L = K.construction(); L
            Integer Ring
            sage: c(L)
            17-adic Ring with capped relative precision 8
            sage: K == c(L)
            True

        TESTS::

            sage: R = ZpLC(13,(31,41))
            sage: R._precision_cap()
            (31, 41)
            sage: F, Z = R.construction()
            sage: S = F(Z)
            sage: S._precision_cap()
            (31, 41)
        """
        from sage.categories.pushout import CompletionFunctor
        extras = {'print_mode':self._printer.dict(), 'type':self._prec_type(), 'names':self._names}
        if hasattr(self, '_label'):
            extras['label'] = self._label
        return (CompletionFunctor(self.prime(), self._precision_cap(), extras), ZZ)

    def random_element(self, algorithm='default'):
        r"""
        Returns a random element of self, optionally using the
        algorithm argument to decide how it generates the
        element. Algorithms currently implemented:

        - default: Choose `a_i`, `i >= 0`, randomly between `0` and
          `p-1` until a nonzero choice is made. Then continue choosing
          `a_i` randomly between `0` and `p-1` until we reach
          precision_cap, and return `\sum a_i p^i`.

        EXAMPLES::

            sage: Zp(5,6).random_element()
            3 + 3*5 + 2*5^2 + 3*5^3 + 2*5^4 + 5^5 + O(5^6)
            sage: ZpCA(5,6).random_element()
            4*5^2 + 5^3 + O(5^6)
            sage: ZpFM(5,6).random_element()
            2 + 4*5^2 + 2*5^4 + 5^5
        """
        if (algorithm == 'default'):
            if self.is_capped_relative():
                i = 0
                a_i = ZZ.random_element(self.prime())
                while a_i.is_zero():
                    i += 1
                    a_i = ZZ.random_element(self.prime())
                return self((self.prime()**i)*(a_i + self.prime()*ZZ.random_element(self.prime_pow.pow_Integer_Integer(self.precision_cap()-1))))
            else:
                return self(ZZ.random_element(self.prime_pow.pow_Integer_Integer(self.precision_cap())))
        else:
            raise NotImplementedError("Don't know %s algorithm"%algorithm)

    #def unit_group(self):
    #    raise NotImplementedError

    #def unit_group_gens(self):
    #    raise NotImplementedError

    #def principal_unit_group(self):
    #    raise NotImplementedError

class pAdicFieldBaseGeneric(pAdicBaseGeneric, pAdicFieldGeneric):
    def composite(self, subfield1, subfield2):
        r"""
        Returns the composite of two subfields of self, i.e., the
        largest subfield containing both

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``subfield1`` -- a subfield
        - ``subfield2`` -- a subfield

        OUTPUT:

        - the composite of subfield1 and subfield2

        EXAMPLES::

            sage: K = Qp(17); K.composite(K, K) is K
            True
        """
        #should be overridden for extension fields
        if (subfield1 is self) and (subfield2 is self):
            return self
        raise ValueError("Arguments must be subfields of self.")

    def subfields_of_degree(self, n):
        r"""
        Returns the number of subfields of self of degree `n`

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``n`` -- an integer

        OUTPUT:

        - integer -- the number of subfields of degree ``n`` over self.base_ring()

        EXAMPLES::

            sage: K = Qp(17)
            sage: K.subfields_of_degree(1)
            1
        """
        if n == 1:
            return 1
        else:
            return 0

    def subfield(self, list):
        r"""
        Returns the subfield generated by the elements in list

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``list`` -- a list of elements of ``self``

        OUTPUT:

        - the subfield of ``self`` generated by the elements of list

        EXAMPLES::

            sage: K = Qp(17); K.subfield([K(17), K(1827)]) is K
            True
        """
        for x in list:
            if x not in self:
                raise TypeError("Members of the list of generators must be elements of self.")
        return self

    def construction(self, forbid_frac_field=False):
        """
        Returns the functorial construction of ``self``, namely,
        completion of the rational numbers with respect a given prime.

        Also preserves other information that makes this field unique
        (e.g. precision, rounding, print mode).

        INPUT:

        - ``forbid_frac_field`` -- require a completion functor rather
          than a fraction field functor.  This is used in the
          :meth:`sage.rings.padics.local_generic.LocalGeneric.change` method.

        EXAMPLES::

            sage: K = Qp(17, 8, print_mode='val-unit', print_sep='&')
            sage: c, L = K.construction(); L
            17-adic Ring with capped relative precision 8
            sage: c
            FractionField
            sage: c(L)
            17-adic Field with capped relative precision 8
            sage: K == c(L)
            True

        We can get a completion functor by forbidding the fraction field::

            sage: c, L = K.construction(forbid_frac_field=True); L
            Rational Field
            sage: c
            Completion[17, prec=8]
            sage: c(L)
            17-adic Field with capped relative precision 8
            sage: K == c(L)
            True

        TESTS::

            sage: R = QpLC(13,(31,41))
            sage: R._precision_cap()
            (31, 41)
            sage: F, Z = R.construction()
            sage: S = F(Z)
            sage: S._precision_cap()
            (31, 41)
        """
        from sage.categories.pushout import FractionField, CompletionFunctor
        if forbid_frac_field:
            extras = {'print_mode':self._printer.dict(), 'type':self._prec_type(), 'names':self._names}
            if hasattr(self, '_label'):
                extras['label'] = self._label
            return (CompletionFunctor(self.prime(), self._precision_cap(), extras), QQ)
        else:
            return FractionField(), self.integer_ring()
