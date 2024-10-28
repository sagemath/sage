# sage.doctest: needs sage.modules
"""
Clifford algebra elements

AUTHORS:

- Travis Scrimshaw (2013-09-06): Initial version
- Trevor Karn (2022-07-10): Rewrite multiplication using bitsets
"""
# ****************************************************************************
#       Copyright (C) 2022 Trevor K. Karn <karnx018 at umn.edu>
#                 (C) 2022 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.parent cimport Parent
from sage.data_structures.bitset cimport Bitset
from sage.algebras.weyl_algebra import repr_from_monomials
from sage.data_structures.blas_dict cimport scal
from copy import copy

cdef class CliffordAlgebraElement(IndexedFreeModuleElement):
    """
    An element in a Clifford algebra.

    TESTS::

        sage: Q = QuadraticForm(ZZ, 3, [1, 2, 3, 4, 5, 6])
        sage: Cl.<x,y,z> = CliffordAlgebra(Q)
        sage: elt = ((x^3-z)*x + y)^2
        sage: TestSuite(elt).run()
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: ((x^3-z)*x + y)^2
            -2*x*y*z - x*z + 5*x - 4*y + 2*z + 2
            sage: Cl.zero()
            0
        """
        return repr_from_monomials(self.list(), self._parent._repr_term)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        TESTS::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: latex( ((x^3-z)*x + y)^2 )
            -2  x y z -  x z + 5  x - 4  y + 2  z + 2
            sage: Cl.<x0,x1,x2> = CliffordAlgebra(Q)
            sage: latex(  (x1 - x2)*x0 + 5*x0*x1*x2 )
            5  x_{0} x_{1} x_{2} -  x_{0} x_{1} +  x_{0} x_{2} - 1
        """
        return repr_from_monomials(self.list(), self._parent._latex_term, True)

    cdef _mul_(self, other):
        """
        Return ``self`` multiplied by ``other``.

        INPUT:

        - ``other`` -- element of the same Clifford algebra as ``self``

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: (x^3 - z*y)*x*(y*z + x*y*z)
            x*y*z + y*z - 24*x + 12*y + 2*z - 24
            sage: y*x
            -x*y + 2
            sage: z*x
            -x*z + 3
            sage: z*z
            6
            sage: x*0
            0
            sage: 0*x
            0

        :issue:`34707`::

            sage: Q = QuadraticForm(QQ, 2, [0,5,0])
            sage: C.<p,q> = CliffordAlgebra(Q)
            sage: (q * p) * q
            5*q
            sage: q * (p * q)
            5*q
        """
        Q = self._parent._quadratic_form
        zero = self._parent._base.zero()
        cdef dict next_level, cur, d = {}
        cdef FrozenBitset ml, mr, t
        cdef Py_ssize_t i, j
        cdef CliffordAlgebraElement rhs = <CliffordAlgebraElement> other

        # Special case when multiplying by 0
        if not self._monomial_coefficients:
            return self
        if not rhs._monomial_coefficients:
            return rhs

        # Special case when multiplying by an element of the base ring
        if len(self._monomial_coefficients) == 1:
            ml, cl = next(iter(self._monomial_coefficients.items()))
            if ml.isempty():
                return rhs._mul_term_self(ml, cl)
        if len(rhs._monomial_coefficients) == 1:
            mr, cr = next(iter(rhs._monomial_coefficients.items()))
            if mr.isempty():
                return self._mul_self_term(mr, cr)

        for ml, cl in self:
            # Distribute the current term ``cl`` * ``ml`` over ``other``.
            cur = copy(other._monomial_coefficients)  # The current distribution of the term
            for i in reversed(ml):
                # Distribute the current factor ``e[i]`` (the ``i``-th
                # element of the standard basis).
                next_level = {}
                # At the end of the following for-loop, ``next`` will be
                # the dictionary describing the element
                # ``e[i]`` * (the element described by the dictionary ``cur``)
                # (where ``e[i]`` is the ``i``-th standard basis vector).
                for mr, cr in cur.items():

                    # Commute the factor as necessary until we are in order
                    for j in mr:
                        if i <= j:
                            break
                        # Add the additional term from the commutation
                        # get a non-frozen bitset to manipulate
                        t = Bitset(mr)  # a mutable copy
                        t.discard(j)
                        t = FrozenBitset(t)
                        next_level[t] = next_level.get(t, zero) + cr * Q[i, j]
                        # Note: ``Q[i,j] == Q(e[i]+e[j]) - Q(e[i]) - Q(e[j])`` for
                        # ``i != j``, where ``e[k]`` is the ``k``-th standard
                        # basis vector.
                        cr = -cr
                        if next_level[t] == zero:
                            del next_level[t]

                    # Check to see if we have a squared term or not
                    mr = Bitset(mr)  # temporarily mutable
                    if i in mr:
                        mr.discard(i)
                        cr *= Q[i, i]
                        # Note: ``Q[i,i] == Q(e[i])`` where ``e[i]`` is the
                        # ``i``-th standard basis vector.
                    else:
                        # mr is implicitly sorted
                        mr.add(i)
                    mr = FrozenBitset(mr)  # refreeze it
                    next_level[mr] = next_level.get(mr, zero) + cr
                    if next_level[mr] == zero:
                        del next_level[mr]
                cur = next_level

            # Add the distributed terms to the total
            for index, coeff in cur.items():
                d[index] = d.get(index, zero) + cl * coeff
                if d[index] == zero:
                    del d[index]

        return self.__class__(self.parent(), d)

    cdef CliffordAlgebraElement _mul_self_term(self, FrozenBitset supp, coeff):
        r"""
        Multiply ``self * term`` with the ``term`` having support ``supp``
        and coefficient ``coeff``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: r = sum(E.basis())
            sage: x * y  # indirect doctest
            x*y
            sage: y * x  # indirect doctest
            -x*y
            sage: r * x  # indirect doctest
            x*y*z - x*y - x*z + x
            sage: r * -x  # indirect doctest
            -x*y*z + x*y + x*z - x
            sage: r * (2*x)  # indirect doctest
            2*x*y*z - 2*x*y - 2*x*z + 2*x
            sage: r * y  # indirect doctest
            -x*y*z + x*y - y*z + y
            sage: r * z  # indirect doctest
            x*y*z + x*z + y*z + z
            sage: r * (x*y)  # indirect doctest
            x*y*z + x*y
            sage: r * (-x*y)  # indirect doctest
            -x*y*z - x*y
            sage: r * (x*y*z)  # indirect doctest
            x*y*z
            sage: r * 1 == r  # indirect doctest
            True
            sage: r * -1 == -r  # indirect doctest
            True
            sage: r * 2  # indirect doctest
            2*x*y*z + 2*x*y + 2*x*z + 2*y*z + 2*x + 2*y + 2*z + 2
        """
        if supp.isempty():  # Multiplication by a base ring element
            if coeff == self._parent._base.one():
                return self
            if coeff == -self._parent._base.one():
                return self._neg_()

            return type(self)(self._parent,
                              scal(coeff, self._monomial_coefficients,
                                   factor_on_left=False))

        return type(self)(self._parent, {supp: coeff}) * self

    cdef CliffordAlgebraElement _mul_term_self(self, FrozenBitset supp, coeff):
        r"""
        Multiply ``term * self`` with the ``term`` having support ``supp``
        and coefficient ``coeff``.
        """
        if supp.isempty():  # Multiplication by a base ring element
            if coeff == self._parent._base.one():
                return self
            if coeff == -self._parent._base.one():
                return self._neg_()

            return type(self)(self._parent,
                              scal(coeff, self._monomial_coefficients,
                                   factor_on_left=True))

        return type(self)(self._parent, {supp: coeff}) * self

    def list(self):
        """
        Return the list of monomials and their coefficients in ``self``
        (as a list of `2`-tuples, each of which has the form
        ``(monomial, coefficient)``).

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y
            sage: elt.list()
            [(1, 5), (01, 1)]
        """
        return sorted(self._monomial_coefficients.items(), key=lambda m: (-len(m[0]), list(m[0])))

    def support(self):
        """
        Return the support of ``self``.

        This is the list of all monomials which appear with nonzero
        coefficient in ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y
            sage: elt.support()
            [1, 01]
        """
        return sorted(self._monomial_coefficients, key=lambda x: (-len(x), list(x)))

    def reflection(self):
        r"""
        Return the image of the reflection automorphism on ``self``.

        The *reflection automorphism* of a Clifford algebra is defined
        as the linear endomorphism of this algebra which maps

        .. MATH::

            x_1 \wedge x_2 \wedge \cdots \wedge x_m \mapsto
            (-1)^m x_1 \wedge x_2 \wedge \cdots \wedge x_m.

        It is an algebra automorphism of the Clifford algebra.

        :meth:`degree_negation` is an alias for :meth:`reflection`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y + x*z
            sage: r = elt.reflection(); r
            x*z - 5*x - y
            sage: r.reflection() == elt
            True

        TESTS:

        We check that the reflection is an involution::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: all(x.reflection().reflection() == x for x in Cl.basis())
            True
        """
        return self.__class__(self._parent,
                              {m: (-1)**len(m) * c for m, c in self})

    degree_negation = reflection

    def transpose(self):
        r"""
        Return the transpose of ``self``.

        The transpose is an anti-algebra involution of a Clifford algebra
        and is defined (using linearity) by

        .. MATH::

            x_1 \wedge x_2 \wedge \cdots \wedge x_m \mapsto
            x_m \wedge \cdots \wedge x_2 \wedge x_1.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y + x*z
            sage: t = elt.transpose(); t
            -x*z + 5*x + y + 3
            sage: t.transpose() == elt
            True
            sage: Cl.one().transpose()
            1

        TESTS:

        We check that the transpose is an involution::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: all(x.transpose().transpose() == x for x in Cl.basis())
            True

        Zero is sent to zero::

            sage: Cl.zero().transpose() == Cl.zero()
            True
        """
        P = self._parent
        if not self._monomial_coefficients:
            return P.zero()
        g = P.gens()
        return P.sum(c * P.prod(g[i] for i in reversed(m)) for m, c in self)

    def conjugate(self):
        r"""
        Return the Clifford conjugate of ``self``.

        The Clifford conjugate of an element `x` of a Clifford algebra is
        defined as

        .. MATH::

            \bar{x} := \alpha(x^t) = \alpha(x)^t

        where `\alpha` denotes the :meth:`reflection <reflection>`
        automorphism and `t` the :meth:`transposition <transpose>`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y + x*z
            sage: c = elt.conjugate(); c
            -x*z - 5*x - y + 3
            sage: c.conjugate() == elt
            True

        TESTS:

        We check that the conjugate is an involution::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: all(x.conjugate().conjugate() == x for x in Cl.basis())
            True
        """
        return self.reflection().transpose()

    clifford_conjugate = conjugate


cdef class ExteriorAlgebraElement(CliffordAlgebraElement):
    """
    An element of an exterior algebra.
    """
    cdef _mul_(self, other):
        """
        Return ``self`` multiplied by ``other``.

        INPUT:

        - ``other`` -- element of the same exterior algebra as ``self``

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: x*y
            x*y
            sage: y*x
            -x*y
            sage: z*y*x
            -x*y*z
            sage: (x*z)*y
            -x*y*z
            sage: (3*x + y)^2
            0
            sage: (x - 3*y + z/3)^2
            0
            sage: (x+y) * (y+z)
            x*y + x*z + y*z

            sage: E.<x,y,z,w> = ExteriorAlgebra(QQ)
            sage: (x * y) * (w * z)
            -x*y*z*w
            sage: x * y * w * z
            -x*y*z*w
            sage: (z * w) * (x * y)
            x*y*z*w

            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: r = sum(E.basis())
            sage: r*r
            4*a*b*c*d + 4*a*b*c + 4*a*b*d + 4*a*c*d + 4*b*c*d
             + 2*a*b + 2*a*c + 2*a*d + 2*b*c + 2*b*d + 2*c*d
             + 2*a + 2*b + 2*c + 2*d + 1
        """
        cdef Parent P = self._parent
        zero = P._base.zero()
        cdef dict d
        cdef ExteriorAlgebraElement rhs = <ExteriorAlgebraElement> other
        cdef list to_remove

        cdef FrozenBitset ml, mr, t
        cdef Py_ssize_t n, num_cross, tot_cross, i, j

        # Special case: one of them is zero
        if not self._monomial_coefficients:
            return self
        if not rhs._monomial_coefficients:
            return rhs

        # Special case: other is a single term
        if len(rhs._monomial_coefficients) == 1:
            mr, cr = next(iter(rhs._monomial_coefficients.items()))
            return self._mul_self_term(mr, cr)

        # Special case: self is a single term
        if len(self._monomial_coefficients) == 1:
            ml, cl = next(iter(self._monomial_coefficients.items()))
            return rhs._mul_term_self(ml, cl)

        # Do some special processing for the constant monomial in ml
        ml = FrozenBitset()
        if ml in self._monomial_coefficients:
            const_coeff = self._monomial_coefficients[ml]
            d = dict(rhs._monomial_coefficients)  # Make a shallow copy
            to_remove = []
            if const_coeff != P._base.one():
                for k in d:
                    d[k] *= const_coeff
                    if not d[k]:  # there might be zero divisors
                        to_remove.append(k)
                for k in to_remove:
                    del d[k]
        else:
            d = {}

        n = P.ngens()
        for ml, cl in self._monomial_coefficients.items():
            # ml for "monomial on the left"
            if ml.isempty():  # We already handled the trivial element
                continue
            for mr, cr in rhs._monomial_coefficients.items():
                # mr for "monomial on the right"
                if mr.isempty():
                    t = ml
                else:
                    if not ml.isdisjoint(mr):
                        # if they intersect nontrivially, move along.
                        continue
                    t = <FrozenBitset> ml._union(mr)
                    it = iter(mr)
                    j = next(it)

                    num_cross = 0  # keep track of the number of signs
                    tot_cross = 0
                    for i in ml:
                        while i > j:
                            num_cross += 1
                            try:
                                j = next(it)
                            except StopIteration:
                                j = n + 1
                        tot_cross += num_cross
                    if tot_cross % 2:
                        cr = -cr

                val = d.get(t, zero) + cl * cr
                if not val:
                    del d[t]
                else:
                    d[t] = val

        return self.__class__(P, d)

    cdef CliffordAlgebraElement _mul_self_term(self, FrozenBitset supp, coeff):
        r"""
        Multiply ``self * term`` with the ``term`` having support ``supp``
        and coefficient ``coeff``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: r = sum(E.basis())
            sage: x * y  # indirect doctest
            x*y
            sage: y * x  # indirect doctest
            -x*y
            sage: r * x  # indirect doctest
            x*y*z - x*y - x*z + x
            sage: r * -x  # indirect doctest
            -x*y*z + x*y + x*z - x
            sage: r * (2*x)  # indirect doctest
            2*x*y*z - 2*x*y - 2*x*z + 2*x
            sage: r * y  # indirect doctest
            -x*y*z + x*y - y*z + y
            sage: r * z  # indirect doctest
            x*y*z + x*z + y*z + z
            sage: r * (x*y)  # indirect doctest
            x*y*z + x*y
            sage: r * (-x*y)  # indirect doctest
            -x*y*z - x*y
            sage: r * (x*y*z)  # indirect doctest
            x*y*z
            sage: r * 1 == r  # indirect doctest
            True
            sage: r * -1 == -r  # indirect doctest
            True
            sage: r * 2  # indirect doctest
            2*x*y*z + 2*x*y + 2*x*z + 2*y*z + 2*x + 2*y + 2*z + 2
        """
        cdef dict d
        cdef list to_remove
        cdef Py_ssize_t num_cross, tot_cross, i, j
        cdef FrozenBitset ml

        if supp.isempty():  # Multiplication by a base ring element
            if coeff == self._parent._base.one():
                return self
            if coeff == -self._parent._base.one():
                return self._neg_()

            return type(self)(self._parent,
                              scal(coeff, self._monomial_coefficients,
                                   factor_on_left=False))

        n = self._parent.ngens()
        d = {}
        for ml, cl in self._monomial_coefficients.items():
            # ml for "monomial on the left"
            if not ml.isdisjoint(supp):
                # if they intersect nontrivially, move along.
                continue
            t = <FrozenBitset> ml._union(supp)
            it = iter(supp)
            j = next(it)

            num_cross = 0  # keep track of the number of signs
            tot_cross = 0
            for i in ml:
                while i > j:
                    num_cross += 1
                    try:
                        j = next(it)
                    except StopIteration:
                        j = n + 1
                tot_cross += num_cross
            if tot_cross % 2:
                d[t] = -cl
            else:
                d[t] = cl

        if coeff == -self._parent._base.one():
            for k in d:
                d[k] = -d[k]
        elif coeff != self._parent._base.one():
            to_remove = []
            for k in d:
                d[k] *= coeff
                if not d[k]:  # there might be zero divisors
                    to_remove.append(k)
            for k in to_remove:
                del d[k]
        return type(self)(self._parent, d)

    cdef CliffordAlgebraElement _mul_term_self(self, FrozenBitset supp, coeff):
        r"""
        Multiply ``term * self`` with the ``term`` having support ``supp``
        and coefficient ``coeff``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: r = sum(E.basis())
            sage: x * r  # indirect doctest
            x*y*z + x*y + x*z + x
            sage: (-x) * r  # indirect doctest
            -x*y*z - x*y - x*z - x
            sage: (2*x) * r  # indirect doctest
            2*x*y*z + 2*x*y + 2*x*z + 2*x
            sage: y * r  # indirect doctest
            -x*y*z - x*y + y*z + y
            sage: z * r  # indirect doctest
            x*y*z - x*z - y*z + z
            sage: (x*y) * r  # indirect doctest
            x*y*z + x*y
            sage: (-x*y) * r  # indirect doctest
            -x*y*z - x*y
            sage: (x*y*z) * r  # indirect doctest
            x*y*z
            sage: 1 * r == r  # indirect doctest
            True
            sage: -1 * r == -r  # indirect doctest
            True
            sage: 2 * r  # indirect doctest
            2*x*y*z + 2*x*y + 2*x*z + 2*y*z + 2*x + 2*y + 2*z + 2
        """
        cdef dict d
        cdef list to_remove
        cdef Py_ssize_t n, num_cross, tot_cross, i, j
        cdef FrozenBitset mr, t

        if supp.isempty():  # Multiplication by a base ring element
            if coeff == self._parent._base.one():
                return self
            if coeff == -self._parent._base.one():
                return self._neg_()

            return type(self)(self._parent,
                              scal(coeff, self._monomial_coefficients,
                                   factor_on_left=True))

        n = self._parent.ngens()
        d = {}
        mr = FrozenBitset()
        # We need to special case the constant coefficient
        const_coeff = None
        if mr in self._monomial_coefficients:
            const_coeff = self._monomial_coefficients.pop(mr)
            d[supp] = const_coeff

        for mr, cr in self._monomial_coefficients.items():
            # mr for "monomial on the right"
            if not supp.isdisjoint(mr):
                # if they intersect nontrivially, move along.
                continue
            t = <FrozenBitset> supp._union(mr)
            it = iter(mr)
            j = next(it)  # We assume mr is non-empty here

            num_cross = 0  # keep track of the number of signs
            tot_cross = 0
            for i in supp:
                while i > j:
                    num_cross += 1
                    try:
                        j = next(it)
                    except StopIteration:
                        j = n + 1
                tot_cross += num_cross
            if tot_cross % 2:
                d[t] = -cr
            else:
                d[t] = cr

        if coeff == -self._parent._base.one():
            for k in d:
                d[k] = -d[k]
        elif coeff != self._parent._base.one():
            to_remove = []
            for k in d:
                d[k] = coeff * d[k]  # This will work for non-commutative base rings
                if not d[k]:  # there might be zero divisors
                    to_remove.append(k)
            for k in to_remove:
                del d[k]

        # Add back the constant coefficient since we removed it for the special case
        if const_coeff is not None:
            self._monomial_coefficients[FrozenBitset()] = const_coeff
        return type(self)(self._parent, d)

    def reduce(self, I, left=True):
        r"""
        Reduce ``self`` with respect to the elements in ``I``.

        INPUT:

        - ``I`` -- list of exterior algebra elements or an ideal
        - ``left`` -- boolean; if reduce as a left ideal (``True``)
          or right ideal (``False``), ignored if ``I`` is an ideal

        EXAMPLES::

            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: f = (a + b*c) * d
            sage: f.reduce([a + b*c], True)
            2*a*d
            sage: f.reduce([a + b*c], False)
            0

            sage: I = E.ideal([a + b*c])
            sage: f.reduce(I)
            0
        """
        from sage.algebras.clifford_algebra import ExteriorAlgebraIdeal
        if isinstance(I, ExteriorAlgebraIdeal):
            return I.reduce(self)

        f = self
        E = self._parent

        cdef FrozenBitset lm, s
        for g in I:
            lm = g.leading_support()
            reduction = True
            while reduction:
                supp = f.support()
                reduction = False
                for s in supp:
                    if lm <= s:
                        reduction = True
                        mon = E.monomial(s - lm)
                        if left:
                            gp = mon * g
                            f = f - f[s] / gp[s] * gp
                        else:
                            gp = g * mon
                            f = f - f[s] / gp[s] * gp
                        break
        return f

    def interior_product(self, x):
        r"""
        Return the interior product (also known as antiderivation) of
        ``self`` with respect to ``x`` (that is, the element
        `\iota_{x}(\text{self})` of the exterior algebra).

        If `V` is an `R`-module, and if `\alpha` is a fixed element of
        `V^*`, then the *interior product* with respect to `\alpha` is
        an `R`-linear map
        `i_{\alpha} \colon \Lambda(V) \to \Lambda(V)`, determined by
        the following requirements:

        - `i_{\alpha}(v) = \alpha(v)` for all `v \in V = \Lambda^1(V)`,
        - it is a graded derivation of degree `-1`: all `x` and `y`
          in `\Lambda(V)` satisfy

        .. MATH::

            i_{\alpha}(x \wedge y) = (i_{\alpha} x) \wedge y
            + (-1)^{\deg x} x \wedge (i_{\alpha} y).

        It can be shown that this map `i_{\alpha}` is graded of
        degree `-1` (that is, sends `\Lambda^k(V)` into
        `\Lambda^{k-1}(V)` for every `k`).

        When `V` is a finite free `R`-module, the interior product can
        also be defined by

        .. MATH::

            (i_{\alpha} \omega)(u_1, \ldots, u_k)
            = \omega(\alpha, u_1, \ldots, u_k),

        where `\omega \in \Lambda^k(V)` is thought of as an
        alternating multilinear mapping from
        `V^* \times \cdots \times V^*` to `R`.

        Since Sage is only dealing with exterior powers of modules
        of the form `R^d` for some nonnegative integer `d`, the
        element `\alpha \in V^*` can be thought of as an element of
        `V` (by identifying the standard basis of `V = R^d` with its
        dual basis). This is how `\alpha` should be passed to this
        method.

        We then extend the interior product to all
        `\alpha \in \Lambda (V^*)` by

        .. MATH::

            i_{\beta \wedge \gamma} = i_{\gamma} \circ i_{\beta}.

        INPUT:

        - ``x`` -- element of (or coercing into) `\Lambda^1(V)`
          (for example, an element of `V`); this plays the role of
          `\alpha` in the above definition

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: x.interior_product(x)
            1
            sage: (x + x*y).interior_product(2*y)
            -2*x
            sage: (x*z + x*y*z).interior_product(2*y - x)
            -2*x*z - y*z - z
            sage: x.interior_product(E.one())
            x
            sage: E.one().interior_product(x)
            0
            sage: x.interior_product(E.zero())
            0
            sage: E.zero().interior_product(x)
            0

        REFERENCES:

        - :wikipedia:`Exterior_algebra#Interior_product`
        """
        P = self._parent
        return P.sum([c * cx * P.interior_product_on_basis(m, mx)
                      for m, c in self for mx, cx in x])

    antiderivation = interior_product

    def hodge_dual(self):
        r"""
        Return the Hodge dual of ``self``.

        The Hodge dual of an element `\alpha` of the exterior algebra is
        defined as `i_{\alpha} \sigma`, where `\sigma` is the volume
        form
        (:meth:`~sage.algebras.clifford_algebra.ExteriorAlgebra.volume_form`)
        and `i_{\alpha}` denotes the antiderivation function with
        respect to `\alpha` (see :meth:`interior_product` for the
        definition of this).

        .. NOTE::

            The Hodge dual of the Hodge dual of a homogeneous element
            `p` of `\Lambda(V)` equals `(-1)^{k(n-k)} p`, where
            `n = \dim V` and `k = \deg(p) = |p|`.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: x.hodge_dual()
            y*z
            sage: (x*z).hodge_dual()
            -y
            sage: (x*y*z).hodge_dual()
            1
            sage: [a.hodge_dual().hodge_dual() for a in E.basis()]
            [1, x, y, z, x*y, x*z, y*z, x*y*z]
            sage: (x + x*y).hodge_dual()
            y*z + z
            sage: (x*z + x*y*z).hodge_dual()
            -y + 1
            sage: E = ExteriorAlgebra(QQ, 'wxyz')
            sage: [a.hodge_dual().hodge_dual() for a in E.basis()]
            [1, -w, -x, -y, -z, w*x, w*y, w*z, x*y, x*z, y*z,
             -w*x*y, -w*x*z, -w*y*z, -x*y*z, w*x*y*z]
        """
        volume_form = self._parent.volume_form()
        return volume_form.interior_product(self)

    def constant_coefficient(self):
        """
        Return the constant coefficient of ``self``.

        .. TODO::

            Define a similar method for general Clifford algebras once
            the morphism to exterior algebras is implemented.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: elt = 5*x + y + x*z + 10
            sage: elt.constant_coefficient()
            10
            sage: x.constant_coefficient()
            0
        """
        return self._monomial_coefficients.get(self._parent.one_basis(),
                                               self._parent._base.zero())

    def scalar(self, other):
        r"""
        Return the standard scalar product of ``self`` with ``other``.

        The standard scalar product of `x, y \in \Lambda(V)` is
        defined by `\langle x, y \rangle = \langle x^t y \rangle`, where
        `\langle a \rangle` denotes the degree-0 term of `a`, and where
        `x^t` denotes the transpose
        (:meth:`~sage.algebras.clifford_algebra.CliffordAlgebraElement.transpose`)
        of `x`.

        .. TODO::

            Define a similar method for general Clifford algebras once
            the morphism to exterior algebras is implemented.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: elt = 5*x + y + x*z
            sage: elt.scalar(z + 2*x)
            0
            sage: elt.transpose() * (z + 2*x)
            -2*x*y + 5*x*z + y*z
        """
        return (self.transpose() * other).constant_coefficient()

cdef class CohomologyRAAGElement(CliffordAlgebraElement):
    """
    An element in the cohomology of a right-angled Artin group.

    .. SEEALSO::

        :class:`~sage.groups.raag.CohomologyRAAG`
    """
    cdef _mul_(self, other):
        """
        Return ``self`` multiplied by ``other``.

        EXAMPLES::

            sage: # needs sage.graphs sage.groups
            sage: C4 = graphs.CycleGraph(4)
            sage: A = groups.misc.RightAngledArtin(C4)
            sage: H = A.cohomology()
            sage: b = sum(H.basis())
            sage: b * b
            2*e0*e2 + 2*e1*e3 + 2*e0 + 2*e1 + 2*e2 + 2*e3 + 1
        """
        zero = self._parent._base.zero()
        cdef frozenset I = frozenset(self._parent._indices)
        cdef dict d = {}
        cdef list t
        cdef tuple tp
        cdef tuple ml, mr
        cdef Py_ssize_t pos, i, j
        cdef bint negate

        for ml, cl in self._monomial_coefficients.items():
            for mr, cr in other._monomial_coefficients.items():
                # Create the next term
                tp = tuple(sorted(mr + ml))
                if any(tp[i] == tp[i+1] for i in range(len(tp)-1)):  # e_i ^ e_i = 0
                    continue
                if tp not in I:
                    # not an independent set, so this term is also 0
                    continue

                t = list(mr)
                negate = False
                for i in reversed(ml):
                    pos = 0
                    for j in t:
                        assert i != j
                        if i < j:
                            break
                        pos += 1
                        negate = not negate
                    t.insert(pos, i)

                if negate:
                    cr = -cr

                d[tp] = d.get(tp, zero) + cl * cr
                if d[tp] == zero:
                    del d[tp]

        return self.__class__(self._parent, d)
