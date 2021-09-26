# -*- coding: utf-8 -*-
r"""
Puiseux Series Ring

The ring of Puiseux series.

AUTHORS:

- Chris Swierczewski 2016: initial version on https://github.com/abelfunctions/abelfunctions/tree/master/abelfunctions
- Frédéric Chapoton 2016: integration of code
- Travis Scrimshaw, Sebastian Oehms 2019-2020: basic improvements and completions

REFERENCES:

- :wikipedia:`Puiseux_series`
"""


# ****************************************************************************
#       Copyright (c) 2016 Chris Swierczewski
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.cachefunc import cached_method
from sage.rings.infinity import infinity
from sage.rings.puiseux_series_ring_element import PuiseuxSeries
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeRing
from sage.structure.element import parent
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.laurent_series_ring_element import LaurentSeries
from sage.rings.power_series_ring import is_PowerSeriesRing
from sage.rings.power_series_ring_element import PowerSeries


class PuiseuxSeriesRing(UniqueRepresentation, CommutativeRing):
    """
    Rings of Puiseux series.

    EXAMPLES::

        sage: P = PuiseuxSeriesRing(QQ, 'y')
        sage: y = P.gen()
        sage: f = y**(4/3) + y**(-5/6); f
        y^(-5/6) + y^(4/3)
        sage: f.add_bigoh(2)
        y^(-5/6) + y^(4/3) + O(y^2)
        sage: f.add_bigoh(1)
        y^(-5/6) + O(y)
    """
    @staticmethod
    def __classcall__(cls, *args, **kwds):
        r"""
        TESTS::

            sage: L = PuiseuxSeriesRing(QQ, 'q')
            sage: L is PuiseuxSeriesRing(QQ, name='q')
            True
            sage: Lp.<q> = PuiseuxSeriesRing(QQ)
            sage: L is Lp
            True
            sage: loads(dumps(L)) is L
            True

            sage: L.variable_names()
            ('q',)
            sage: L.variable_name()
            'q'
        """
        if not kwds and len(args) == 1 and isinstance(args[0], LaurentSeriesRing):
            laurent_series = args[0]
        else:
            laurent_series = LaurentSeriesRing(*args, **kwds)

        return super(PuiseuxSeriesRing, cls).__classcall__(cls, laurent_series)

    def __init__(self, laurent_series):
        """
        Generic class for Puiseux series rings.

        EXAMPLES::

            sage: P = PuiseuxSeriesRing(QQ, 'y')
            sage: TestSuite(P).run()

            sage: P = PuiseuxSeriesRing(ZZ, 'x')
            sage: TestSuite(P).run()

            sage: P = PuiseuxSeriesRing(ZZ['a,b'], 'x')
            sage: TestSuite(P).run()
        """
        base_ring = laurent_series.base_ring()

        # If self is R(( x^(1/e) )) then the corresponding Laurent series
        # ring will be R(( x ))
        self._laurent_series_ring = laurent_series

        CommutativeRing.__init__(self, base_ring,
                                 names=laurent_series.variable_names(),
                                 category=laurent_series.category())

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: PuiseuxSeriesRing(AA, 'y')
            Puiseux Series Ring in y over Algebraic Real Field
        """
        s = "Puiseux Series Ring in {} over {}".format(self.variable_name(),
                                                       self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    def base_extend(self, R):
        """
        Extend the coefficients.

        INPUT:

        - ``R`` -- a ring

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.base_extend(QQ)
            Puiseux Series Ring in y over Rational Field
        """
        return PuiseuxSeriesRing(self._laurent_series_ring.base_extend(R))

    def change_ring(self, R):
        """
        Return a Puiseux series ring over another ring.

        INPUT:

        - ``R`` -- a ring

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.change_ring(QQ)
            Puiseux Series Ring in y over Rational Field
        """
        return PuiseuxSeriesRing(self._laurent_series_ring.change_ring(R))

    def is_sparse(self):
        """
        Return whether ``self`` is sparse.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.is_sparse()
            False
        """
        return self.laurent_series_ring().is_sparse()

    def is_dense(self):
        """
        Return whether ``self`` is dense.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.is_dense()
            True
        """
        return self.laurent_series_ring().is_dense()

    def is_field(self, proof=True):
        r"""
        Return whether ``self`` is a field.

        A Puiseux series ring is a field if and only
        its base ring is a field.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.is_field()
            False
            sage: A.change_ring(QQ).is_field()
            True
        """
        return self.base_ring().is_field(proof=proof)

    def fraction_field(self):
        r"""
        Return the fraction field of this ring of Laurent series.

        If the base ring is a field, then Puiseux series are already a field.
        If the base ring is a domain, then the Puiseux series over its fraction
        field is returned. Otherwise, raise a ``ValueError``.

        EXAMPLES::

            sage: R = PuiseuxSeriesRing(ZZ, 't', 30).fraction_field()
            sage: R
            Puiseux Series Ring in t over Rational Field
            sage: R.default_prec()
            30

            sage: PuiseuxSeriesRing(Zmod(4), 't').fraction_field()
            Traceback (most recent call last):
            ...
            ValueError: must be an integral domain
        """
        from sage.categories.integral_domains import IntegralDomains
        from sage.categories.fields import Fields
        if self in Fields():
            return self
        elif self in IntegralDomains():
            return PuiseuxSeriesRing(self._laurent_series_ring.fraction_field())
        else:
            raise ValueError('must be an integral domain')

    def residue_field(self):
        r"""
        Return the residue field of this Puiseux series field
        if it is a complete discrete valuation field (i.e. if
        the base ring is a field, in which case it is also the
        residue field).

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(GF(17))
            sage: R.residue_field()
            Finite Field of size 17

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: R.residue_field()
            Traceback (most recent call last):
            ...
            TypeError: the base ring is not a field
        """
        if not self.base_ring().is_field():
            raise TypeError("the base ring is not a field")
        return self.base_ring()

    def uniformizer(self):
        r"""
        Return a uniformizer of this Puiseux series field if it is
        a discrete valuation field (i.e. if the base ring is actually
        a field). Otherwise, an error is raised.

        EXAMPLES::

            sage: R.<t> = PuiseuxSeriesRing(QQ)
            sage: R.uniformizer()
            t

            sage: R.<t> = PuiseuxSeriesRing(ZZ)
            sage: R.uniformizer()
            Traceback (most recent call last):
            ...
            TypeError: the base ring is not a field
        """
        if not self.base_ring().is_field():
            raise TypeError("the base ring is not a field")
        return self.gen()

    Element = PuiseuxSeries

    def _element_constructor_(self, x, e=1, prec=infinity):
        r"""
        Construct a Puiseux series from ``x``.

        INPUT:

        - ``x`` -- an object that can be converted into a Puiseux series
        - ``e`` -- (default: ``1``) the ramification index of the series
        - ``prec`` -- (default: ``infinity``) the precision of the series
          as a rational number

        EXAMPLES::

            sage: P = PuiseuxSeriesRing(QQ, 'y')
            sage: y = P.gen()
            sage: P([1,3,5,7])
            1 + 3*y + 5*y^2 + 7*y^3
            sage: P(33/14)
            33/14

            sage: Q = PowerSeriesRing(QQ, 'y')
            sage: z = Q([1,2,4,5]).O(6); z
            1 + 2*y + 4*y^2 + 5*y^3 + O(y^6)
            sage: P(z) + y**(1/2)
            1 + y^(1/2) + 2*y + 4*y^2 + 5*y^3 + O(y^6)

            sage: Q = LaurentSeriesRing(QQ, 'y')
            sage: z = Q([3,2,1,2]).add_bigoh(5); z
            3 + 2*y + y^2 + 2*y^3 + O(y^5)
            sage: P(z) + y**(1/2)
            3 + y^(1/2) + 2*y + y^2 + 2*y^3 + O(y^5)

            sage: from sage.modular.etaproducts import qexp_eta
            sage: y^(1/24)*qexp_eta(P, prec=30)
            y^(1/24) - y^(25/24) - y^(49/24) + y^(121/24) + y^(169/24) - y^(289/24) - y^(361/24) + y^(529/24) + y^(625/24) + O(y^(721/24))
        """
        P = parent(x)

        # 1. x is a Puiseux series belonging to this ring.
        #    This is short-circuited by the coercion framework.
        #if isinstance(x, self.element_class) and P is self:
        #    return x
        # 2. x is a Puiseux series but not an element of this ring. the laurent
        #    part should be coercible to the laurent series ring of self
        if isinstance(x, self.element_class):
            l = self._laurent_series_ring(x.laurent_part())
            e = x.ramification_index()
        # 3. x is a member of the base ring then convert x to a laurent series
        #    and set the ramification index of the Puiseux series to 1.
        elif P is self.base_ring():
            l = self._laurent_series_ring(x)
            e = 1
        # 4. x is a Laurent or power series with the same base ring
        elif (isinstance(x, (LaurentSeries, PowerSeries))
              and P is self.base_ring()):
            l = self._laurent_series_ring(x)
        # 5. everything else: try to coerce to laurent series ring
        else:
            l = self._laurent_series_ring(x)

        # finally, construct an instance of the element class and adding
        # the precision afterwards (see also :trac:`28993`).
        return self.element_class(self, l, e=e).add_bigoh(prec)

    def _coerce_map_from_(self, P):
        r"""
        Return a coercion map from `P` to `self`, or `True`, or `None`.

        The following rings admit a coercion map to the Puiseux series ring
        `A((x-a)^(1/e))`:

        - any ring that admits a coercion map to `A`

        - any Laurent series ring, power series ring, or polynomial ring in the
          variable `(x-a)` over a ring admitting a coercion map to `A`

        - any Puiseux series ring with the same center `a` and ramification
          index equal to a multiple of `self`'s ramification index. For
          example, Puiseux series in (x-a)^(1/2) can be interpreted as Puiseux
          series in (x-a)^(1/4).

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: 5 in R, 1/5 in R              # indirect doctests
            (True, False)

            sage: p = x^(1/2) + x**3-x**(-1/4)
            sage: p.laurent_part() in R         # indirect doctests
            True

            sage: Q.<x> = PuiseuxSeriesRing(QQ) # indirect doctests
            sage: p in Q
            True
        """
        # any ring that has a coercion map to A
        A = self.base_ring()
        if A is P:
            return True
        f = A.coerce_map_from(P)
        if f is not None:
            return self.coerce_map_from(A) * f

        # Laurent series rings, power series rings, and polynomial rings with
        # the same variable name and the base rings are coercible
        if ((isinstance(P, PuiseuxSeriesRing) or isinstance(P, LaurentSeriesRing) or
             is_PowerSeriesRing(P)) and
                P.variable_name() == self.variable_name() and
                A.has_coerce_map_from(P.base_ring())):
            return True

        # # other Puiseux series rings with the same variable name and
        # # center. Puiseux series rings with difference ramification indices are
        # # coercible to each other.
        # if (is_PuiseuxSeriesRing(P) and
        #     P.variable_name() == self.variable_name()):
        #     return True

    @cached_method
    def gen(self, n=0):
        r"""
        Return the generator of ``self``.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(AA, 'z')
            sage: A.gen()
            z
        """
        if n != 0:
            raise IndexError("generator {} not defined".format(n))
        return self.element_class(self, [0, 1], e=1)

    def ngens(self):
        r"""
        Return the number of generators of ``self``, namely 1.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(AA, 'z')
            sage: A.ngens()
            1
        """
        return 1

    def laurent_series_ring(self):
        r"""
        Return the underlying Laurent series ring.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(AA, 'z')
            sage: A.laurent_series_ring()
            Laurent Series Ring in z over Algebraic Real Field
        """
        return self._laurent_series_ring

    def default_prec(self):
        """
        Return the default precision of ``self``.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(AA, 'z')
            sage: A.default_prec()
            20
        """
        return self.laurent_series_ring().default_prec()

