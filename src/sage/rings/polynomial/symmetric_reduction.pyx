# sage.doctest: needs sage.libs.singular
r"""
Symmetric Reduction of Infinite Polynomials

:class:`~sage.rings.polynomial.symmetric_reduction.SymmetricReductionStrategy`
provides a framework for efficient symmetric reduction of Infinite
Polynomials, see :mod:`~sage.rings.polynomial.infinite_polynomial_element`.

AUTHORS:

- Simon King <simon.king@nuigalway.ie>

THEORY:

According to M. Aschenbrenner and C. Hillar [AB2007]_, Symmetric
Reduction of an element `p` of an Infinite Polynomial Ring `X` by some
other element `q` means the following:

    1. Let `M` and `N` be the leading terms of `p` and `q`.
    2. Test whether there is a permutation `P` that does not
       diminish the variable indices occurring in `N`
       and preserves their order, so that there is some term
       `T\in X` with `T N^P = M`. If there is no such permutation,
       return `p`.
    3. Replace `p` by `p-T q^P` and continue with step 1.


When reducing one polynomial `p` with respect to a list `L` of other
polynomials, there usually is a choice of order on which the
efficiency crucially depends. Also it helps to modify the polynomials
on the list in order to simplify the basic reduction steps.

The preparation of `L` may be expensive. Hence, if the same list is
used many times then it is reasonable to perform the preparation only
once. This is the background of
:class:`~sage.rings.polynomial.symmetric_reduction.SymmetricReductionStrategy`.

Our current strategy is to keep the number of terms in the polynomials
as small as possible. For this, we sort `L` by increasing number of
terms. If several elements of `L` allow for a reduction of `p`, we
choose the one with the smallest number of terms. Later on, it should
be possible to implement further strategies for choice.

When adding a new polynomial `q` to `L`, we first reduce `q` with
respect to `L`. Then, we test heuristically whether it is possible to
reduce the number of terms of the elements of `L` by reduction modulo
`q`.  That way, we see best chances to keep the number of terms in
intermediate reduction steps relatively small.

EXAMPLES:

First, we create an infinite polynomial ring and one of its elements::

    sage: X.<x,y> = InfinitePolynomialRing(QQ)
    sage: p = y[1]*y[3] + y[1]^2*x[3]

We want to symmetrically reduce it by another polynomial. So, we put
this other polynomial into a list and create a Symmetric Reduction
Strategy object::

    sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
    sage: S = SymmetricReductionStrategy(X, [y[2]^2*x[1]])
    sage: S
    Symmetric Reduction Strategy in
     Infinite polynomial ring in x, y over Rational Field, modulo
        x_1*y_2^2
    sage: S.reduce(p)
    x_3*y_1^2 + y_3*y_1

The preceding is correct, since any permutation that turns
``y[2]^2*x[1]`` into a factor of ``y[1]^2*x[3]`` interchanges the
variable indices 1 and 2 -- which is not allowed in a symmetric
reduction. However, reduction by ``y[1]^2*x[2]`` works, since one can
change variable index 1 into 2 and 2 into 3. So, we add this to
``S``::

    sage: S.add_generator(y[1]^2*x[2])
    sage: S
    Symmetric Reduction Strategy in
     Infinite polynomial ring in x, y over Rational Field, modulo
        x_2*y_1^2,
        x_1*y_2^2
    sage: S.reduce(p)                                                                   # needs sage.combinat
    y_3*y_1

The next example shows that tail reduction is not done, unless it is
explicitly advised::

    sage: S.reduce(x[3] + 2*x[2]*y[1]^2 + 3*y[2]^2*x[1])                                # needs sage.combinat
    x_3 + 2*x_2*y_1^2 + 3*x_1*y_2^2
    sage: S.tailreduce(x[3] + 2*x[2]*y[1]^2 + 3*y[2]^2*x[1])                            # needs sage.combinat
    x_3

However, it is possible to ask for tailreduction already when the
Symmetric Reduction Strategy is created::

    sage: S2 = SymmetricReductionStrategy(X, [y[2]^2*x[1],y[1]^2*x[2]], tailreduce=True)
    sage: S2
    Symmetric Reduction Strategy in
     Infinite polynomial ring in x, y over Rational Field, modulo
        x_2*y_1^2,
        x_1*y_2^2
    with tailreduction
    sage: S2.reduce(x[3] + 2*x[2]*y[1]^2 + 3*y[2]^2*x[1])                               # needs sage.combinat
    x_3

"""

# ****************************************************************************
#       Copyright (C) 2009 Simon King <king@mathematik.nuigalway.ie>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import sys
from sage.structure.richcmp cimport richcmp, Py_NE, Py_EQ


cdef class SymmetricReductionStrategy:
    """
    A framework for efficient symmetric reduction of InfinitePolynomial, see
    :mod:`~sage.rings.polynomial.infinite_polynomial_element`.

    INPUT:

    - ``Parent`` -- an Infinite Polynomial Ring, see
      :mod:`~sage.rings.polynomial.infinite_polynomial_element`.
    - ``L`` -- (list, default the empty list) List of elements of ``Parent``
      with respect to which will be reduced.
    - ``good_input`` -- (bool, default ``None``) If this optional parameter
      is true, it is assumed that each element of ``L`` is symmetrically
      reduced with respect to the previous elements of ``L``.

    EXAMPLES::

        sage: X.<y> = InfinitePolynomialRing(QQ)
        sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
        sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]], good_input=True)
        sage: S.reduce(y[3] + 2*y[2]*y[1]^2 + 3*y[2]^2*y[1])
        y_3 + 3*y_2^2*y_1 + 2*y_2*y_1^2
        sage: S.tailreduce(y[3] + 2*y[2]*y[1]^2 + 3*y[2]^2*y[1])                        # needs sage.combinat
        y_3

    """
    def __init__(self, Parent, L=None, tailreduce=False, good_input=None):
        """
        EXAMPLES::

            sage: X.<y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]], good_input=True)
            sage: S == loads(dumps(S))
            True

        """
        self._parent = Parent
        if hasattr(Parent, '_P'):
            self._R = Parent._P
        else:
            self._R = None
        self._lm = []
        self._lengths = []
        self._min_lm = None
        self._tail = int(tailreduce)
        if not (L is None):
            for p in L:
                self.add_generator(p, good_input=good_input)

    def __getinitargs__(self):
        r"""
        Used for pickling.

        EXAMPLES::

            sage: X.<y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]], good_input=True)
            sage: S.__getinitargs__()
            (Infinite polynomial ring in y over Rational Field, [], 0, None)
        """
        return (self._parent, [], self._tail, None)

    def __getstate__(self):
        r"""
        Used for pickling.

        EXAMPLES::

            sage: X.<y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]], good_input=True)
            sage: S.__getstate__()
            ([y_2*y_1^2, y_2^2*y_1], [1, 1], y_2*y_1^2, 0, Infinite polynomial ring in y over Rational Field)
        """
        # Apparently, for pickling it is needed to update self._lm and
        # self._min_lm before calling dumps...
        R = self._parent
        self._lm = [R(x) for x in self._lm]  # I have no idea why -- but it seems needed

        self._min_lm = R(self._min_lm)
        return (self._lm, self._lengths, self._min_lm,
                self._tail, self._parent)

    def __setstate__(self, L):  # (lm, lengths, min_lm, tail)
        r"""
        Used for pickling.

        EXAMPLES::

            sage: X.<y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]], good_input=True)
            sage: S == loads(dumps(S)) # indirect doctest
            True
        """
        self._lm = L[0]
        self._lengths = L[1]
        self._min_lm = L[2]
        self._tail = L[3]
        self._parent = L[4]
        if hasattr(self._parent, '_P'):
            self._R = self._parent._P
        else:
            self._R = None

    def __richcmp__(self, other, op):
        r"""
        Standard comparison function.

        EXAMPLES::

            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]], tailreduce=True)
            sage: S == 17
            False
            sage: S == SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]], tailreduce=False)
            False
            sage: S == SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]], tailreduce=True)
            True
        """
        if not isinstance(other, SymmetricReductionStrategy):
            if op in [Py_NE, Py_EQ]:
                return (op == Py_NE)
            else:
                return NotImplemented
        cdef SymmetricReductionStrategy left = self
        cdef SymmetricReductionStrategy right = other
        return richcmp((left._parent, left._lm, left._tail),
                       (right._parent, right._lm, right._tail), op)

    def gens(self):
        """
        Return the list of Infinite Polynomials modulo which self reduces.

        EXAMPLES::

            sage: X.<y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]])
            sage: S
            Symmetric Reduction Strategy in
             Infinite polynomial ring in y over Rational Field, modulo
                y_2*y_1^2,
                y_2^2*y_1
            sage: S.gens()
            [y_2*y_1^2, y_2^2*y_1]

        """
        return self._lm

    def setgens(self, L):
        """
        Define the list of Infinite Polynomials modulo which self reduces.

        INPUT:

        ``L`` -- a list of elements of the underlying infinite polynomial ring.

        .. NOTE::

            It is not tested if ``L`` is a good input. That method simply
            assigns a *copy* of ``L`` to the generators of self.

        EXAMPLES::

            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: X.<y> = InfinitePolynomialRing(QQ)
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]])
            sage: R = SymmetricReductionStrategy(X)
            sage: R.setgens(S.gens())
            sage: R
            Symmetric Reduction Strategy in
             Infinite polynomial ring in y over Rational Field, modulo
                y_2*y_1^2,
                y_2^2*y_1
            sage: R.gens() is S.gens()
            False
            sage: R.gens() == S.gens()
            True

        """
        self._lm = [X for X in L]

    def reset(self):
        """
        Remove all polynomials from ``self``.

        EXAMPLES::

            sage: X.<y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1],y[1]^2*y[2]])
            sage: S
            Symmetric Reduction Strategy in
             Infinite polynomial ring in y over Rational Field, modulo
                y_2*y_1^2,
                y_2^2*y_1
            sage: S.reset()
            sage: S
            Symmetric Reduction Strategy in Infinite polynomial ring in y over Rational Field

        """
        self._lm = []
        self._lengths = []
        self._min_lm = None

    def __repr__(self):
        """
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: S = SymmetricReductionStrategy(X, [y[2]^2*y[1], y[1]^2*y[2]], tailreduce=True)
            sage: S  # indirect doctest
            Symmetric Reduction Strategy in Infinite polynomial ring in x, y over Rational Field, modulo
                y_2*y_1^2,
                y_2^2*y_1
            with tailreduction

        """
        s = "Symmetric Reduction Strategy in %s" % self._parent
        if self._lm:
            s += ", modulo\n    %s" % (',\n    '.join(str(X) for X in self._lm))
        if self._tail:
            s += '\nwith tailreduction'
        return s

    def __call__(self, p):
        """
        INPUT:

        A polynomial or an infinite polynomial

        OUTPUT:

        A polynomial whose parent ring allows for coercion of any generator of self

        EXAMPLES::

            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: X.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: a, b = y[2]^2*y[1], y[1]^2*y[2]
            sage: p = y[3]*x[2]*x[1]
            sage: S = SymmetricReductionStrategy(X, [a,b])
            sage: p._p.parent().has_coerce_map_from(a._p.parent())
            False
            sage: q = S(p)
            sage: q.parent().has_coerce_map_from(a._p.parent())
            True
            sage: S(p) == S(p._p)
            True

        """
        if hasattr(p, '_p'):
            p = p._p
        if self._R is None:
            self._R = p.parent()
            if hasattr(self._parent, '_P'):
                self._parent._P = self._R
            return p
        if self._R.has_coerce_map_from(p.parent()):
            return self._R(p)
        if p.parent().has_coerce_map_from(self._R):
            self._R = p.parent()
            if hasattr(self._parent, '_P'):
                self._parent._P = self._R
            return p
        # now we really need to work...
        R = self._R
        VarList = list(set(list(R.variable_names()) + list(p.parent().variable_names())))
        VarList.sort(key=self._parent.varname_key, reverse=True)
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        self._R = PolynomialRing(self._parent.base_ring(), VarList,
                                 order=self._parent._order)
        if hasattr(self._parent, '_P'):
            self._parent._P = self._R
        return self._R(p)

    def add_generator(self, p, good_input=None):
        """
        Add another polynomial to ``self``.

        INPUT:

        - ``p`` -- An element of the underlying infinite polynomial ring.
        - ``good_input`` -- (bool, default ``None``) If ``True``, it is
          assumed that ``p`` is reduced with respect to ``self``. Otherwise,
          this reduction will be done first (which may cost some time).

        .. NOTE::

            Previously added polynomials may be modified. All input is
            prepared in view of an efficient symmetric reduction.

        EXAMPLES::

            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: S = SymmetricReductionStrategy(X)
            sage: S
            Symmetric Reduction Strategy in
             Infinite polynomial ring in x, y over Rational Field
            sage: S.add_generator(y[3] + y[1]*(x[3]+x[1]))
            sage: S
            Symmetric Reduction Strategy in
             Infinite polynomial ring in x, y over Rational Field, modulo
                x_3*y_1 + x_1*y_1 + y_3

        Note that the first added polynomial will be simplified when
        adding a suitable second polynomial::

            sage: S.add_generator(x[2] + x[1])                                          # needs sage.combinat
            sage: S                                                                     # needs sage.combinat
            Symmetric Reduction Strategy in
             Infinite polynomial ring in x, y over Rational Field, modulo
                y_3,
                x_2 + x_1

        By default, reduction is applied to any newly added
        polynomial. This can be avoided by specifying the optional
        parameter 'good_input'::

            sage: # needs sage.combinat
            sage: S.add_generator(y[2] + y[1]*x[2])
            sage: S
            Symmetric Reduction Strategy in
             Infinite polynomial ring in x, y over Rational Field, modulo
                y_3,
                x_1*y_1 - y_2,
                x_2 + x_1
            sage: S.reduce(x[3] + x[2])
            -2*x_1
            sage: S.add_generator(x[3] + x[2], good_input=True)
            sage: S
            Symmetric Reduction Strategy in
             Infinite polynomial ring in x, y over Rational Field, modulo
                y_3,
                x_3 + x_2,
                x_1*y_1 - y_2,
                x_2 + x_1

        In the previous example, ``x[3] + x[2]`` is added without
        being reduced to zero.

        """
        from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
        p = InfinitePolynomial(self._parent, self(p))
        cdef SymmetricReductionStrategy tmpStrategy
        if good_input is None:
            p = self.reduce(p)
        if p._p == 0:
            return
        cdef int i = 0
        cdef int l = len(self._lm)
        cdef int newLength = len(p._p.coefficients())
        p = p / p.lc()
        if (self._min_lm is None) or (p.lm() < self._min_lm):
            self._min_lm = p.lm()
        while (i < l) and (self._lengths[i] < newLength):
            i += 1
        self._lm.insert(i, p)
        self._lengths.insert(i, newLength)
        # return
        i += 1
        l += 1
        if i < l:
            tmpStrategy = SymmetricReductionStrategy(self._parent, [p],
                                                     tailreduce=False,
                                                     good_input=True)
        else:
            return
        cdef int j
        while i < l:
            q = tmpStrategy.reduce(self._lm[i].lm()) + tmpStrategy.reduce(self._lm[i].tail())
            if q._p == 0:
                self._lm.pop(i)
                self._lengths.pop(i)
                l -= 1
                i -= 1
            else:
                q_len = len(q._p.coefficients())
                if q_len < self._lengths[i]:
                    self._lm.pop(i)
                    self._lengths.pop(i)
                    j = 0
                    while (j < i) and (self._lengths[j] < q_len):
                        j += 1
                    self._lm.insert(j, q)
                    self._lengths.insert(j, q_len)
            i += 1

    def reduce(self, p, notail=False, report=None):
        """
        Symmetric reduction of an infinite polynomial.

        INPUT:

        - ``p`` -- an element of the underlying infinite polynomial ring.
        - ``notail`` -- (bool, default ``False``) If ``True``, tail reduction
          is avoided (but there is no guarantee that there will be no tail
          reduction at all).
        - ``report`` -- (object, default ``None``) If not ``None``, print
          information on the progress of the computation.

        OUTPUT:

        Reduction of ``p`` with respect to ``self``.

        .. NOTE::

            If tail reduction shall be forced, use :meth:`tailreduce`.

        EXAMPLES::

            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: S = SymmetricReductionStrategy(X, [y[3]], tailreduce=True)
            sage: S.reduce(y[4]*x[1] + y[1]*x[4])
            x_4*y_1
            sage: S.reduce(y[4]*x[1] + y[1]*x[4], notail=True)
            x_4*y_1 + x_1*y_4

        Last, we demonstrate the ``report`` option::

            sage: S = SymmetricReductionStrategy(X, [x[2] + y[1],
            ....:                                    x[2]*y[3] + x[1]*y[2] + y[4],
            ....:                                    y[3] + y[2]])
            sage: S
            Symmetric Reduction Strategy in
             Infinite polynomial ring in x, y over Rational Field, modulo
                y_3 + y_2,
                x_2 + y_1,
                x_1*y_2 + y_4 - y_3*y_1
            sage: S.reduce(x[3] + x[1]*y[3] + x[1]*y[1], report=True)
            :::>
            x_1*y_1 + y_4 - y_3*y_1 - y_1

        Each ':' indicates that one reduction of the leading monomial
        was performed. Eventually, the '>' indicates that the
        computation is finished.

        """
        from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
        cdef list lml = self._lm
        if not lml:
            if report is not None:
                print('>')
            return p
        if p.lm() < self._min_lm:
            if report is not None:
                print('>')
            return p
        cdef list REDUCTOR
        while True:
            REDUCTOR = []
            for q in lml:
                c, P, _ = q.symmetric_cancellation_order(p)
                if c is not None and c <= 0:
                    REDUCTOR = [self(q ** P)]
                    break
            if not REDUCTOR:
                new_p = p
                break
            p = self(p)  # now this is a usual polynomial
            R = self._R
            if hasattr(p, 'reduce'):
                new_p = InfinitePolynomial(self._parent,
                                           p.reduce([R(X) for X in REDUCTOR]))
            else:
                new_p = InfinitePolynomial(self._parent, p % (REDUCTOR * R))
            if report is not None:
                sys.stdout.write(':')
                sys.stdout.flush()
            if (new_p._p == p) or (new_p._p == 0):
                break
            p = new_p  # now this is an infinite polynomial
        p = new_p
        if (not self._tail) or notail or (p._p == 0):
            if report is not None:
                print('>')
            return p
        # there remains to perform tail reduction
        REM = p.lt()
        p = p.tail()
        p = self.tailreduce(p, report=report)
        if report is not None:
            print('>')
        return p + REM

    def tailreduce(self, p, report=None):
        """
        Symmetric reduction of an infinite polynomial, with forced tail reduction.

        INPUT:

        - ``p`` -- an element of the underlying infinite polynomial ring.
        - ``report`` -- (object, default ``None``) If not ``None``, print
          information on the progress of the computation.

        OUTPUT:

        Reduction (including the non-leading elements) of ``p`` with respect to ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: S = SymmetricReductionStrategy(X, [y[3]])
            sage: S.reduce(y[4]*x[1] + y[1]*x[4])
            x_4*y_1 + x_1*y_4
            sage: S.tailreduce(y[4]*x[1] + y[1]*x[4])                                   # needs sage.combinat
            x_4*y_1

        Last, we demonstrate the 'report' option::

            sage: S = SymmetricReductionStrategy(X, [x[2] + y[1],
            ....:                                    x[2]*x[3] + x[1]*y[2] + y[4],
            ....:                                    y[3] + y[2]])
            sage: S
            Symmetric Reduction Strategy in
             Infinite polynomial ring in x, y over Rational Field, modulo
                y_3 + y_2,
                x_2 + y_1,
                x_1*y_2 + y_4 + y_1^2
            sage: S.tailreduce(x[3] + x[1]*y[3] + x[1]*y[1], report=True)               # needs sage.combinat
            T[3]:::>
            T[3]:>
            x_1*y_1 - y_2 + y_1^2 - y_1

        The protocol means the following.
         * 'T[3]' means that we currently do tail reduction for a polynomial
           with three terms.
         * ':::>' means that there were three reductions of leading terms.
         * The tail of the result of the preceding reduction still has three
           terms. One reduction of leading terms was possible, and then the
           final result was obtained.
        """
        if not self._lm:
            return p
        OUT = p.parent()(0)
        while p._p != 0:
            if report is not None:
                sys.stdout.write('T[%d]' % len(p._p.coefficients()))
                sys.stdout.flush()
            p = self.reduce(p, notail=True, report=report)
            OUT = OUT + p.lt()
            p = p.tail()
            if p.lm() < self._min_lm:
                return OUT + p
        return OUT
