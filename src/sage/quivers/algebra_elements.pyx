"""
Path algebra elements

AUTHORS:

- Simon King (2015-08)
"""

# ****************************************************************************
#     Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
include "algebra_elements.pxi"

from sage.misc.repr import repr_lincomb
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool


cdef class PathAlgebraElement(RingElement):
    """
    Elements of a :class:`~sage.quivers.algebra.PathAlgebra`.

    .. NOTE::

        Upon creation of a path algebra, one can choose among several
        monomial orders, which are all positive or negative degree
        orders. Monomial orders that are not degree orders are not
        supported.

    EXAMPLES:

    After creating a path algebra and getting hold of its generators, one can
    create elements just as usual::

        sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
        sage: A.inject_variables()
        Defining e_0, e_1, e_2, a, b, c, d, e, f
        sage: x = a+2*b+3*c+5*e_0+3*e_2
        sage: x
        5*e_0 + a + 2*b + 3*c + 3*e_2

    The path algebra decomposes as a direct sum according to start- and endpoints::

        sage: x.sort_by_vertices()
        [(5*e_0, 0, 0),
         (a, 0, 1),
         (2*b, 0, 2),
         (3*c, 1, 0),
         (3*e_2, 2, 2)]
        sage: (x^3+x^2).sort_by_vertices()
        [(150*e_0 + 33*a*c, 0, 0),
         (30*a + 3*a*c*a, 0, 1),
         (114*b + 6*a*c*b, 0, 2),
         (90*c + 9*c*a*c, 1, 0),
         (18*c*a, 1, 1),
         (54*c*b, 1, 2),
         (36*e_2, 2, 2)]

    For a consistency test, we create a path algebra that is isomorphic to a
    free associative algebra, and compare arithmetic with two other
    implementations of free algebras (note that the letterplace implementation
    only allows weighted homogeneous elements)::

        sage: F.<x,y,z> = FreeAlgebra(GF(25,'t'))
        sage: pF = x+y*z*x+2*y-z+1
        sage: pF2 = x^4+x*y*x*z+2*z^2*x*y
        sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'))
        sage: pP = sage_eval('x+y*z*x+2*y-z+1', P.gens_dict())
        sage: pP^5+3*pP^3 == sage_eval(repr(pF^5+3*pF^3), P.gens_dict())
        True
        sage: L.<x,y,z> = FreeAlgebra(GF(25,'t'), implementation='letterplace')
        sage: pL2 = x^4+x*y*x*z+2*z^2*x*y
        sage: pP2 = sage_eval('x^4+x*y*x*z+2*z^2*x*y', P.gens_dict())
        sage: pP2^7 == sage_eval(repr(pF2^7), P.gens_dict())
        True
        sage: pP2^7 == sage_eval(repr(pL2^7), P.gens_dict())
        True

    When the Cython implementation of path algebra elements was
    introduced, it was faster than both the default implementation and
    the letterplace implementation of free algebras. The following
    timings where obtained with a 32-bit operating system; using 64-bit
    on the same machine, the letterplace implementation has not become
    faster, but the timing for path algebra elements has improved by
    about 20%::

        sage: # not tested
        sage: timeit('pF^5+3*pF^3')
        1 loops, best of 3: 338 ms per loop
        sage: timeit('pP^5+3*pP^3')
        100 loops, best of 3: 2.55 ms per loop
        sage: timeit('pF2^7')
        10000 loops, best of 3: 513 ms per loop
        sage: timeit('pL2^7')
        125 loops, best of 3: 1.99 ms per loop
        sage: timeit('pP2^7')
        10000 loops, best of 3: 1.54 ms per loop

    So, if one is merely interested in basic arithmetic operations for
    free associative algebras, it could make sense to model the free
    associative algebra as a path algebra. However, standard basis
    computations are not available for path algebras, yet. Hence, to
    implement computations in graded quotients of free algebras, the
    letterplace implementation currently is the only option.
    """
    def __cinit__(self):
        """
        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = a+2*b+3*c+5*e_0+3*e_2   # indirect doctest
            sage: x
            5*e_0 + a + 2*b + 3*c + 3*e_2
        """
        self.data = NULL

    def __dealloc__(self):
        """
        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = a+2*b+3*c+5*e_0+3*e_2
            sage: del x       # indirect doctest
        """
        homog_poly_free(self.data)

    def __init__(self, S, data):
        """
        Do not call directly.

        INPUT:

        - ``S`` -- path algebra

        - ``data`` -- dictionary; most of its keys are
          :class:`~sage.quivers.paths.QuiverPath`, the value giving its
          coefficient

        .. NOTE::

            Monomial orders that are not degree orders are not supported.

        EXAMPLES::

            sage: P1 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'))
            sage: P1.inject_variables()     # indirect doctest
            Defining e_1, x, y, z
            sage: (x+2*z+1)^2
            e_1 + 4*z + 2*x + 4*z*z + 2*x*z + 2*z*x + x*x
            sage: P2 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P2.inject_variables()
            Defining e_1, x, y, z
            sage: (x+2*z+1)^2
            4*z*z + 2*x*z + 2*z*x + x*x + 4*z + 2*x + e_1
            sage: P3 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='negdeglex')
            sage: P3.inject_variables()
            Defining e_1, x, y, z
            sage: (x+2*z+1)^2
            e_1 + 4*z + 2*x + 4*z*z + 2*z*x + 2*x*z + x*x
            sage: P4 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='deglex')
            sage: P4.inject_variables()
            Defining e_1, x, y, z
            sage: (x+2*z+1)^2
            4*z*z + 2*z*x + 2*x*z + x*x + 4*z + 2*x + e_1
        """
        self._hash = -1
        order = S.order_string()
        if order == "negdegrevlex":
            self.cmp_terms = negdegrevlex
        elif order == "degrevlex":
            self.cmp_terms = degrevlex
        elif order == "negdeglex":
            self.cmp_terms = negdeglex
        elif order == "deglex":
            self.cmp_terms = deglex
        else:
            raise ValueError("unknown term order '{}'".format(order))
        cdef QuiverPath tmp = None
        RingElement.__init__(self, S)
        cdef dict homog = {}
        cdef list L
        for tmp, c in data.iteritems():
            sig_check()
            homog.setdefault((tmp.initial_vertex(),
                              tmp.terminal_vertex()), []).append((tmp, c))
        cdef path_homog_poly_t *HP
        for (s, e), L in sorted(homog.iteritems(), reverse=True):
            sig_check()
            HP = homog_poly_init_list(s, e, L, self.cmp_terms, -1)
            HP.nxt = self.data
            self.data = HP

    def __reduce__(self):
        """
        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'))
            sage: p = sage_eval('(x+2*z+1)^3', P.gens_dict())
            sage: loads(dumps(p)) == p     # indirect doctest
            True
        """
        return path_algebra_element_unpickle, (self._parent, homog_poly_pickle(self.data))

    cdef list _sorted_items_for_printing(self):
        """
        Return list of pairs ``(M,c)``, where ``c`` is a coefficient and ``M``
        will be passed to ``self.parent()._repr_monomial`` resp. to
        ``self.parent()._latex_monomial``, providing the indices of the
        algebra generators occurring in the monomial.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: X         # indirect doctest
            5*e_0 + a + 2*b + 3*c + 3*e_2
            sage: latex(X)  # indirect doctest
            5 e_0 + a + 2 b + 3 c + 3 e_2
        """
        cdef path_homog_poly_t *H = self.data
        cdef list L, L_total
        cdef size_t i
        cdef path_term_t * T
        L_total = []
        cdef list vertices = self._parent.quiver().vertices(sort=True)
        cdef mp_size_t offset = len(vertices)
        while H != NULL:
            L = []  # data for a single component (given by start- and endpoints)
            T = H.poly.lead
            while T != NULL:
                sig_check()
                if T.mon.path.length:
                    L.append(([offset + biseq_getitem(T.mon.path, i) for i in range(<size_t>T.mon.path.length)],
                              <object>(T.coef)))
                else:
                    L.append(([vertices.index(H.start)], <object>(T.coef)))
                T = T.nxt
            if <size_t>len(L) != H.poly.nterms:
                print("Term count of polynomial is wrong, got", len(L),
                      "expected", H.poly.nterms)
            L_total.extend(L)
            H = H.nxt
        return L_total

    def _repr_(self):
        """
        String representation.

        .. NOTE::

            The terms are first sorted by initial and terminal
            vertices, and only then by the given monomial order.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: X         # indirect doctest
            5*e_0 + a + 2*b + 3*c + 3*e_2
        """
        return repr_lincomb(self._sorted_items_for_printing(), strip_one=True,
            scalar_mult=self.parent()._print_options['scalar_mult'],
            repr_monomial=self._parent._repr_monomial)

    def _latex_(self):
        r"""
        Latex string representation.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: latex(X)  # indirect doctest
            5 e_0 + a + 2 b + 3 c + 3 e_2
            sage: latex(X*X)
            10 e_0 + 3 a\cdot c + 5 a + b + 3 c\cdot a + 6 c\cdot b + 9 e_2
        """
        return repr_lincomb(self._sorted_items_for_printing(),
            scalar_mult=self.parent()._print_options['scalar_mult'],
            latex_scalar_mult=self.parent()._print_options['latex_scalar_mult'],
            repr_monomial=self._parent._latex_monomial,
            is_latex=True, strip_one=True)

    # Basic properties

    def __bool__(self):
        """
        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: bool(a+b+c+d)   # indirect doctest
            True
            sage: bool(((a+b+c+d)-(a+b))-(c+d))
            False
        """
        return self.data != NULL

    def __len__(self):
        """
        Return the number of terms appearing in this element.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: X = a+2*b+3*c+5*e_0+3*e_2
            sage: len(X)
            5
            sage: len(X^5)
            17
        """
        cdef size_t l = 0
        cdef path_homog_poly_t *H = self.data
        while H != NULL:
            sig_check()
            l += H.poly.nterms
            H = H.nxt
        return l

    cpdef ssize_t degree(self) except -2:
        """
        Return the degree, provided the element is homogeneous.

        An error is raised if the element is not homogeneous.

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'))
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: q = (x+y+2*z)^3
            sage: q.degree()
            3
            sage: p = (x+2*z+1)^3
            sage: p.degree()
            Traceback (most recent call last):
            ...
            ValueError: element is not homogeneous
        """
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef mp_size_t deg = 0
        cdef bint zero = True
        while H != NULL:
            sig_check()
            T = H.poly.lead
            while T != NULL:
                if zero:
                    deg = term_total_degree(T)
                elif deg != term_total_degree(T):
                    raise ValueError("element is not homogeneous")
                zero = False
                T = T.nxt
            H = H.nxt
        if zero:
            return -1
        return deg

    def is_homogeneous(self):
        """
        Tells whether this element is homogeneous.

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'))
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: q = (x+y+2*z)^3
            sage: q.is_homogeneous()
            True
            sage: p = (x+2*z+1)^3
            sage: p.is_homogeneous()
            False
        """
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef mp_size_t deg = 0
        cdef bint zero = True
        while H != NULL:
            T = H.poly.lead
            while T != NULL:
                sig_check()
                if zero:
                    deg = term_total_degree(T)
                elif deg != term_total_degree(T):
                    return False
                zero = False
                T = T.nxt
            H = H.nxt
        return True

    cpdef dict monomial_coefficients(self):
        """
        Return the dictionary keyed by the monomials appearing
        in this element, the values being the coefficients.

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: p = (x+2*z+1)^3
            sage: sorted(p.monomial_coefficients().items())
            [(x*x*x, 1),
             (z*x*x, 2),
             (x*z*x, 2),
             (z*z*x, 4),
             (x*x*z, 2),
             (z*x*z, 4),
             (x*z*z, 4),
             (z*z*z, 3),
             (x*x, 3),
             (z*x, 1),
             (x*z, 1),
             (z*z, 2),
             (x, 3),
             (z, 1),
             (e_1, 1)]

        Note that the dictionary can be fed to the algebra, to reconstruct the
        element::

            sage: P(p.monomial_coefficients()) == p
            True
        """
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        cdef dict D = {}
        while H != NULL:
            T = H.poly.lead
            while T != NULL:
                tmp = sample._new_(H.start, H.end)
                biseq_init_copy(tmp._path, T.mon.path)
                D[tmp] = <object>T.coef
                T = T.nxt
            H = H.nxt
        return D

    cpdef list coefficients(self):
        """
        Return the list of coefficients.

        .. NOTE::

            The order in which the coefficients are returned corresponds to the
            order in which the terms are printed. That is *not* the same as the
            order given by the monomial order, since the terms are first ordered
            according to initial and terminal vertices, before applying the
            monomial order of the path algebra.

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: p = (x+2*z+1)^3
            sage: p
            3*z*z*z + 4*x*z*z + 4*z*x*z + 2*x*x*z + 4*z*z*x + 2*x*z*x + 2*z*x*x + x*x*x + 2*z*z + x*z + z*x + 3*x*x + z + 3*x + e_1
            sage: p.coefficients()
            [3, 4, 4, 2, 4, 2, 2, 1, 2, 1, 1, 3, 1, 3, 1]
        """
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef list L = []
        while H != NULL:
            T = H.poly.lead
            while T != NULL:
                L.append(<object>T.coef)
                T = T.nxt
            H = H.nxt
        return L

    cpdef list monomials(self):
        """
        Return the list of monomials appearing in this element.

        .. NOTE::

            The order in which the monomials are returned corresponds to the
            order in which the element's terms are printed. That is *not* the
            same as the order given by the monomial order, since the terms are
            first ordered according to initial and terminal vertices, before
            applying the monomial order of the path algebra.

            The monomials are not elements of the underlying partial
            semigroup, but of the algebra.

        .. SEEALSO:: :meth:`support`

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: p = (x+2*z+1)^3
            sage: p
            3*z*z*z + 4*x*z*z + 4*z*x*z + 2*x*x*z + 4*z*z*x + 2*x*z*x + 2*z*x*x + x*x*x + 2*z*z + x*z + z*x + 3*x*x + z + 3*x + e_1
            sage: p.monomials()
            [z*z*z,
             x*z*z,
             z*x*z,
             x*x*z,
             z*z*x,
             x*z*x,
             z*x*x,
             x*x*x,
             z*z,
             x*z,
             z*x,
             x*x,
             z,
             x,
             e_1]
            sage: p.monomials()[1].parent() is P
            True
        """
        cdef path_homog_poly_t *H = self.data
        cdef path_homog_poly_t *out
        cdef path_term_t *T
        cdef object one = self.base_ring().one()
        cdef list L = []
        while H != NULL:
            T = H.poly.lead
            while T != NULL:
                out = homog_poly_create(H.start, H.end)
                out.poly.lead = term_create_blank(one)
                mon_copy(out.poly.lead.mon, T.mon)
                out.poly.lead.nxt = NULL
                out.poly.nterms = 1
                L.append(self._new_(out))
                T = T.nxt
            H = H.nxt
        return L

    cpdef list terms(self):
        """
        Return the list of terms.

        .. NOTE::

            The order in which the terms are returned corresponds to the order
            in which they are printed. That is *not* the same as the
            order given by the monomial order, since the terms are first
            ordered according to initial and terminal vertices, before
            applying the monomial order of the path algebra.

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: p = (x+2*z+1)^3
            sage: p
            3*z*z*z + 4*x*z*z + 4*z*x*z + 2*x*x*z + 4*z*z*x + 2*x*z*x + 2*z*x*x + x*x*x + 2*z*z + x*z + z*x + 3*x*x + z + 3*x + e_1
            sage: p.terms()
            [3*z*z*z,
             4*x*z*z,
             4*z*x*z,
             2*x*x*z,
             4*z*z*x,
             2*x*z*x,
             2*z*x*x,
             x*x*x,
             2*z*z,
             x*z,
             z*x,
             3*x*x,
             z,
             3*x,
             e_1]
        """
        cdef path_homog_poly_t *H = self.data
        cdef path_homog_poly_t *out
        cdef path_term_t *T
        cdef list L = []
        while H != NULL:
            T = H.poly.lead
            while T != NULL:
                out = homog_poly_create(H.start, H.end)
                out.poly.lead = term_copy(T)
                out.poly.lead.nxt = NULL
                out.poly.nterms = 1
                L.append(self._new_(out))
                T = T.nxt
            H = H.nxt
        return L

    cpdef list support(self):
        """
        Return the list of monomials, as elements of the underlying partial semigroup.

        .. NOTE::

            The order in which the monomials are returned corresponds to the
            order in which the element's terms are printed. That is *not* the
            same as the order given by the monomial order, since the terms are
            first ordered according to initial and terminal vertices, before
            applying the monomial order of the path algebra.

        .. SEEALSO:: :meth:`monomials`

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: p = (x+2*z+1)^3
            sage: p
            3*z*z*z + 4*x*z*z + 4*z*x*z + 2*x*x*z + 4*z*z*x + 2*x*z*x + 2*z*x*x + x*x*x + 2*z*z + x*z + z*x + 3*x*x + z + 3*x + e_1
            sage: p.support()
            [z*z*z,
             x*z*z,
             z*x*z,
             x*x*z,
             z*z*x,
             x*z*x,
             z*x*x,
             x*x*x,
             z*z,
             x*z,
             z*x,
             x*x,
             z,
             x,
             e_1]
            sage: p.support()[1].parent() is P.semigroup()
            True
        """
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        cdef list L = []
        while H != NULL:
            T = H.poly.lead
            while T != NULL:
                tmp = sample._new_(H.start, H.end)
                biseq_init_copy(tmp._path, T.mon.path)
                L.append(tmp)
                T = T.nxt
            H = H.nxt
        return L

    def support_of_term(self):
        """
        If ``self`` consists of a single term, return the corresponding
        element of the underlying path semigroup.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = 4*a*d*c*b*e
            sage: x.support_of_term()
            a*d*c*b*e
            sage: x.support_of_term().parent() is A.semigroup()
            True
            sage: (x + f).support_of_term()
            Traceback (most recent call last):
            ...
            ValueError: 4*a*d*c*b*e + f is not a single term
        """
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        if self.data != NULL and self.data.nxt == NULL:
            if self.data.poly.lead != NULL:
                tmp = sample._new_(self.data.start, self.data.end)
                biseq_init_copy(tmp._path, self.data.poly.lead.mon.path)
                return tmp
        raise ValueError("{} is not a single term".format(self))

    cpdef object coefficient(self, QuiverPath P):
        """
        Return the coefficient of a monomial.

        INPUT:

        - ``P`` -- an element of the underlying partial semigroup

        OUTPUT:

        The coefficient of the given semigroup element in ``self``, or zero if
        it does not appear.

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: p = (x+2*z+1)^3
            sage: p
            3*z*z*z + 4*x*z*z + 4*z*x*z + 2*x*x*z + 4*z*z*x + 2*x*z*x + 2*z*x*x + x*x*x + 2*z*z + x*z + z*x + 3*x*x + z + 3*x + e_1
            sage: p.coefficient(sage_eval('x*x*z', P.semigroup().gens_dict()))
            2
            sage: p.coefficient(sage_eval('z*x*x*x', P.semigroup().gens_dict()))
            0
        """
        if self.data == NULL:
            return self.base_ring().zero()
        H = homog_poly_get_predecessor_of_component(self.data, P._start, P._end)
        if H == NULL:
            if self.data.start != P._start or self.data.end != P._end:
                return self.base_ring().zero()
            H = self.data
        else:
            H = H.nxt
        if H == NULL:
            return self.base_ring().zero()
        # Now, H points to the component that belongs to K
        cdef path_mon_t pM
        mon_create_keep(pM, P._path, -1, 0, 0)
        T = H.poly.lead
        while T != NULL:
            if self.cmp_terms(T.mon, pM) == 0:
                return <object>T.coef
            T = T.nxt
        return self.base_ring().zero()

    def __iter__(self):
        """
        Iterate over the pairs (monomial, coefficient) appearing in ``self``.

        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: p = (x+2*z+1)^3
            sage: p
            3*z*z*z + 4*x*z*z + 4*z*x*z + 2*x*x*z + 4*z*z*x + 2*x*z*x + 2*z*x*x + x*x*x + 2*z*z + x*z + z*x + 3*x*x + z + 3*x + e_1
            sage: list(p)   # indirect doctest
            [(z*z*z, 3),
             (x*z*z, 4),
             (z*x*z, 4),
             (x*x*z, 2),
             (z*z*x, 4),
             (x*z*x, 2),
             (z*x*x, 2),
             (x*x*x, 1),
             (z*z, 2),
             (x*z, 1),
             (z*x, 1),
             (x*x, 3),
             (z, 1),
             (x, 3),
             (e_1, 1)]
        """
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        while H != NULL:
            T = H.poly.lead
            while T != NULL:
                sig_check()
                tmp = sample._new_(H.start, H.end)
                biseq_init_copy(tmp._path, T.mon.path)
                yield (tmp, <object>T.coef)
                T = T.nxt
            H = H.nxt

    cdef PathAlgebraElement _new_(self, path_homog_poly_t *h):
        """
        Create a new path algebra element from C interface data.
        """
        cdef PathAlgebraElement out = type(self).__new__(type(self))
        out._parent = self._parent
        out.cmp_terms = self.cmp_terms
        out.data = h
        out._hash = -1
        return out

    def __copy__(self):
        """
        EXAMPLES::

            sage: P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P.inject_variables()
            Defining e_1, x, y, z
            sage: p = (x+2*z+1)^3
            sage: copy(p) is p
            False
            sage: copy(p) == p   # indirect doctest
            True
        """
        return self._new_(homog_poly_copy(self.data))

    def __getitem__(self, k):
        """
        Either return the coefficient in ``self`` of an element of the
        underlying partial semigroup, or the sum of terms of ``self`` whose
        monomials have a given initial and terminal vertex.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: X = (a+2*b+3*c+5*e_0+3*e_2)^3
            sage: X[A.semigroup()('c')]
            75
            sage: X.sort_by_vertices()
            [(125*e_0 + 30*a*c, 0, 0),
             (25*a + 3*a*c*a, 0, 1),
             (98*b + 6*a*c*b, 0, 2),
             (75*c + 9*c*a*c, 1, 0),
             (15*c*a, 1, 1),
             (48*c*b, 1, 2),
             (27*e_2, 2, 2)]
            sage: X.sort_by_vertices()
            [(125*e_0 + 30*a*c, 0, 0),
             (25*a + 3*a*c*a, 0, 1),
             (98*b + 6*a*c*b, 0, 2),
             (75*c + 9*c*a*c, 1, 0),
             (15*c*a, 1, 1),
             (48*c*b, 1, 2),
             (27*e_2, 2, 2)]
            sage: X[0,2]
            98*b + 6*a*c*b
        """
        cdef path_homog_poly_t *H
        cdef path_term_t *T
        cdef path_mon_t kM
        cdef PathAlgebraElement out
        cdef QuiverPath K
        if isinstance(k, tuple):
            H = homog_poly_get_predecessor_of_component(self.data, k[0], k[1])
            if H == NULL:
                if self.data.start == k[0] and self.data.end == k[1]:
                    out = self._new_(homog_poly_create(self.data.start, self.data.end))
                    out.data.nxt = NULL
                    poly_icopy(out.data.poly, self.data.poly)
                else:
                    return self._new_(NULL)
            else:
                if H.nxt == NULL or H.nxt.start != k[0] or H.nxt.end != k[1]:
                    return self._new_(NULL)
                out = self._new_(homog_poly_create(H.nxt.start, H.nxt.end))
                out.data.nxt = NULL
                poly_icopy(out.data.poly, H.nxt.poly)
            return out
        elif isinstance(k, QuiverPath):
            if self.data == NULL:
                return self.base_ring().zero()
            K = k
            H = homog_poly_get_predecessor_of_component(self.data, K._start, K._end)
            if H == NULL:
                if self.data.start != K._start or self.data.end != K._end:
                    return self.base_ring().zero()
                H = self.data
            else:
                H = H.nxt
            if H == NULL:
                return self.base_ring().zero()
            # Now, H points to the component that belongs to K
            mon_create_keep(kM, K._path, -1, 0, 0)
            T = H.poly.lead
            while T != NULL:
                sig_check()
                if self.cmp_terms(T.mon, kM) == 0:
                    return <object>T.coef
                T = T.nxt
        return self.base_ring().zero()

    def sort_by_vertices(self):
        """
        Return a list of triples ``(element, v1, v2)``, where ``element`` is
        an element whose monomials all have initial vertex ``v1`` and terminal
        vertex ``v2``, so that the sum of elements is ``self``.

        EXAMPLES::

            sage: A1 = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: A1.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = (b*e*b*e+4*b+e_0)^2
            sage: y = (a*c*b+1)^3
            sage: x.sort_by_vertices()
            [(e_0 + 2*b*e*b*e + b*e*b*e*b*e*b*e, 0, 0), (4*b + 4*b*e*b*e*b, 0, 2)]
            sage: sum(c[0] for c in x.sort_by_vertices()) == x
            True
            sage: y.sort_by_vertices()
            [(e_0, 0, 0), (3*a*c*b, 0, 2), (e_1, 1, 1), (e_2, 2, 2)]
            sage: sum(c[0] for c in y.sort_by_vertices()) == y
            True
        """
        cdef path_homog_poly_t * H = self.data
        cdef PathAlgebraElement out
        cdef list C = []
        while H != NULL:
            out = self._new_(homog_poly_create(H.start, H.end))
            out.data.nxt = NULL
            sig_check()
            poly_icopy(out.data.poly, H.poly)
            C.append((out, H.start, H.end))
            H = H.nxt
        return C

    ####
    # Arithmetics
    # Hash and Comparison
    def __hash__(self):
        """
        The hash is cached, to make it faster.

        EXAMPLES::

            sage: P1 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(3,'t'))
            sage: P2 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(3,'t'), order='deglex')
            sage: P1.inject_variables()
            Defining e_1, x, y, z
            sage: p = x+y
            sage: P2.inject_variables()
            Defining e_1, x, y, z
            sage: q = x+y
            sage: D = dict([(p^i,i) for i in range(1,8)])
            sage: len(D)
            7
            sage: hash(q^5) == hash(p^5)    # indirect doctest
            True
            sage: D[q^6]
            6
        """
        if self._hash == -1:
            self._hash = hash(frozenset(self.monomial_coefficients().items()))
        return self._hash

    cpdef _richcmp_(left, right, int op):
        """
        Helper for comparison of path algebra elements.

        .. NOTE::

            First, the comparison is by initial vertices of
            monomials. Then, the terminal vertices are compared. Last,
            the given monomial order is applied for monomials that
            have the same initial and terminal vertices.

        EXAMPLES::

            sage: A1 = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: A1.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = (b*e*b*e+4*b+e_0)^2
            sage: y = (a*c*b+1)^3
            sage: x.sort_by_vertices()
            [(e_0 + 2*b*e*b*e + b*e*b*e*b*e*b*e, 0, 0), (4*b + 4*b*e*b*e*b, 0, 2)]
            sage: y.sort_by_vertices()
            [(e_0, 0, 0), (3*a*c*b, 0, 2), (e_1, 1, 1), (e_2, 2, 2)]

        The two elements are distinguished by monomials with initial and
        terminal vertex `0`. Hence, `x` should evaluate bigger than `y`::

            sage: x > y    # indirect doctest
            True
        """
        cdef PathAlgebraElement other = right
        cdef PathAlgebraElement self = left
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2 = other.data
        while H1 != NULL and H2 != NULL:
            v1 = H1.start
            v2 = H2.start
            if v1 != v2:
                return richcmp_not_equal(v1, v2, op)

            v1 = H1.end
            v2 = H2.end
            if v1 != v2:
                return richcmp_not_equal(v1, v2, op)

            w1 = H1.poly
            w2 = H2.poly
            if w1 != w2:
                return poly_richcmp(H1.poly, H2.poly, self.cmp_terms, op)

            H1 = H1.nxt
            H2 = H2.nxt
        if H1 == NULL:
            if H2 == NULL:
                return rich_to_bool(op, 0)
            return rich_to_bool(op, -1)
        return rich_to_bool(op, 1)

    # negation
    cpdef _neg_(self):
        """
        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(GF(3))
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = b*e*b*e+4*b*e+e_0
            sage: -x    # indirect doctest
            2*e_0 + 2*b*e + 2*b*e*b*e
        """
        return self._new_(homog_poly_neg(self.data))

    # addition
    cpdef _add_(self, other):
        """
        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(GF(3))
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = b*e*b*e+4*b*e+e_0
            sage: y = a*c+1
            sage: x+y    # indirect doctest
            2*e_0 + b*e + a*c + b*e*b*e + e_1 + e_2
        """
        cdef PathAlgebraElement right = other
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2 = right.data
        cdef path_poly_t *P
        cdef path_homog_poly_t *out = NULL
        cdef path_homog_poly_t *tmp
        while True:
            sig_check()
            if H1 == NULL:
                if out == NULL:
                    if H2 == NULL:
                        return self._new_(NULL)
                    return self._new_(homog_poly_copy(H2))
                else:
                    if H2 != NULL:
                        # If out is not NULL then tmp isn't either
                        tmp.nxt = homog_poly_copy(H2)
                    return self._new_(out)
            elif H2 == NULL:
                if out == NULL:
                    if H1 == NULL:
                        return self._new_(NULL)
                    return self._new_(homog_poly_copy(H1))
                else:
                    if H1 != NULL:
                        # If out is not NULL then tmp isn't either
                        tmp.nxt = homog_poly_copy(H1)
                    return self._new_(out)
            else:
                if (H1.start > H2.start) or (H1.start == H2.start and H1.end > H2.end):
                    if out == NULL:
                        out = homog_poly_create(H2.start, H2.end)
                        poly_icopy(out.poly, H2.poly)
                        tmp = out
                    else:
                        tmp.nxt = homog_poly_create(H2.start, H2.end)
                        tmp = tmp.nxt
                        poly_icopy(tmp.poly, H2.poly)
                    H2 = H2.nxt
                elif (H1.start < H2.start) or (H1.end < H2.end):
                    if out == NULL:
                        out = homog_poly_create(H1.start, H1.end)
                        poly_icopy(out.poly, H1.poly)
                        tmp = out
                    else:
                        tmp.nxt = homog_poly_create(H1.start, H1.end)
                        tmp = tmp.nxt
                        poly_icopy(tmp.poly, H1.poly)
                    H1 = H1.nxt
                else:
                    # start- and endpoints match
                    P = poly_add(H1.poly, H2.poly, self.cmp_terms)
                    if P.lead != NULL:
                        if out == NULL:
                            out = homog_poly_init_poly(H1.start, H1.end, P)
                            tmp = out
                        else:
                            tmp.nxt = homog_poly_init_poly(H1.start, H1.end, P)
                            tmp = tmp.nxt
                    else:
                        poly_free(P)
                    H1 = H1.nxt
                    H2 = H2.nxt

    cpdef _sub_(self, other):
        """
        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(GF(3))
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = b*e*b*e+4*b*e+1
            sage: y = a*c-1  # indirect doctest
            sage: x-y        # indirect doctest
            2*e_0 + b*e + 2*a*c + b*e*b*e
        """
        cdef PathAlgebraElement right = other
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2 = right.data
        cdef path_poly_t *P
        cdef path_homog_poly_t *out = NULL
        cdef path_homog_poly_t *tmp
        while True:
            sig_check()
            if H1 == NULL:
                if out == NULL:
                    if H2 == NULL:
                        return self._new_(NULL)
                    sig_check()
                    return self._new_(homog_poly_copy(H2))
                else:
                    if H2 != NULL:
                        # If out is not NULL then tmp isn't either
                        tmp.nxt = homog_poly_copy(H2)
                    return self._new_(out)
            elif H2 == NULL:
                if out == NULL:
                    if H1 == NULL:
                        return self._new_(NULL)
                    return self._new_(homog_poly_copy(H1))
                else:
                    if H1 != NULL:
                        # If out is not NULL then tmp isn't either
                        tmp.nxt = homog_poly_copy(H1)
                    return self._new_(out)
            else:
                if (H1.start > H2.start) or (H1.start == H2.start and H1.end > H2.end):
                    if out == NULL:
                        sig_on()
                        out = homog_poly_create(H2.start, H2.end)
                        poly_icopy(out.poly, H2.poly)
                        sig_off()
                        tmp = out
                    else:
                        sig_on()
                        tmp.nxt = homog_poly_create(H2.start, H2.end)
                        tmp = tmp.nxt
                        poly_icopy(tmp.poly, H2.poly)
                        sig_off()
                    H2 = H2.nxt
                elif (H1.start < H2.start) or (H1.end < H2.end):
                    if out == NULL:
                        sig_on()
                        out = homog_poly_create(H1.start, H1.end)
                        poly_icopy(out.poly, H1.poly)
                        sig_off()
                        tmp = out
                    else:
                        sig_on()
                        tmp.nxt = homog_poly_create(H1.start, H1.end)
                        tmp = tmp.nxt
                        poly_icopy(tmp.poly, H1.poly)
                        sig_off()
                    H1 = H1.nxt
                else:
                    # start- and endpoints match
                    sig_on()
                    P = poly_sub(H1.poly, H2.poly, self.cmp_terms)
                    if P.lead != NULL:
                        if out == NULL:
                            out = homog_poly_init_poly(H1.start, H1.end, P)
                            tmp = out
                        else:
                            tmp.nxt = homog_poly_init_poly(H1.start, H1.end, P)
                            tmp = tmp.nxt
                    else:
                        poly_free(P)
                    sig_off()
                    H1 = H1.nxt
                    H2 = H2.nxt

# (scalar) multiplication

    cpdef _lmul_(self, Element right):
        """
        EXAMPLES::

            sage: from sage.quivers.algebra_elements import PathAlgebraElement
            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: x = sage_eval('3*a+3*b+3*c+3*e_0+3*e_2', A.gens_dict())
            sage: x*2   # indirect doctest
            6*e_0 + 6*a + 6*b + 6*c + 6*e_2

        ::

            sage: z = sage_eval('a+2*b+5*c+5*e_0+3*e_2', A.gens_dict())
            sage: z
            5*e_0 + a + 2*b + 5*c + 3*e_2
            sage: z*3
            3*a + 6*b + 9*e_2
        """
        if self.data == NULL:
            return self
        cdef path_homog_poly_t * out = homog_poly_scale(self.data, right)
        cdef path_homog_poly_t * outnxt
        if out.poly.nterms == 0:
            # homog_poly_scale will remove zero components, except the first.
            # Thus, we can return self._new_(out.nxt), but need to free the
            # memory occupied by out first.
            outnxt = out.nxt
            poly_free(out.poly)
            sig_free(out)
            return self._new_(outnxt)
        return self._new_(out)

    cpdef _rmul_(self, Element left):
        """
        EXAMPLES::

            sage: from sage.quivers.algebra_elements import PathAlgebraElement
            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: x = sage_eval('3*a+3*b+3*c+3*e_0+3*e_2', A.gens_dict())
            sage: 2*x   # indirect doctest
            6*e_0 + 6*a + 6*b + 6*c + 6*e_2

        ::

            sage: z = sage_eval('a+2*b+5*c+5*e_0+3*e_2', A.gens_dict())
            sage: z
            5*e_0 + a + 2*b + 5*c + 3*e_2
            sage: 3*z
            3*a + 6*b + 9*e_2
        """
        if self.data == NULL:
            return self
        cdef path_homog_poly_t * out = homog_poly_scale(self.data, left)
        cdef path_homog_poly_t * outnxt
        if out.poly.nterms == 0:
            # homog_poly_scale will remove zero components, except the first.
            # Thus, we can return self._new_(out.nxt), but need to free the
            # memory occupied by out first.
            outnxt = out.nxt
            poly_free(out.poly)
            sig_free(out)
            return self._new_(outnxt)
        return self._new_(out)

    def __truediv__(self, x):
        """
        Division by coefficients.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: X/2
            10*e_0 + 8*a + b + 9*c + 9*e_2
            sage: (X/2)*2 == X    # indirect doctest
            True

        ::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: X = a+2*b+3*c+5*e_0+3*e_2
            sage: X/4
            5/4*e_0 + 1/4*a + 1/2*b + 3/4*c + 3/4*e_2
            sage: (X/4).parent()
            Path algebra of Looped multi-digraph on 3 vertices over Rational Field
            sage: (X/4)*4 == X
            True
        """
        cdef PathAlgebraElement sample
        if isinstance(self, PathAlgebraElement):
            sample = self
            x = ~(sample._parent._base(x))
            if x.parent() is not sample._parent._base:
                sample = sample._parent._semigroup.algebra(x.parent())(0)
            return sample._new_(homog_poly_scale((<PathAlgebraElement>self).data, x))
        raise TypeError("do not know how to divide {} by {}".format(x, self))

# Multiplication in the algebra

    cpdef _mul_(self, other):
        """
        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = b*e*b*e+4*b*e+e_0
            sage: y = a*c+5*f*e
            sage: x*y
            a*c + 4*b*e*a*c + b*e*b*e*a*c
            sage: y*x
            a*c + 4*a*c*b*e + a*c*b*e*b*e + 5*f*e + 5*f*e*b*e + 5*f*e*b*e*b*e
            sage: y*y
            a*c*a*c + 5*f*e*a*c
            sage: x*x
            e_0 + 8*b*e + 3*b*e*b*e + 8*b*e*b*e*b*e + b*e*b*e*b*e*b*e

        ::

            sage: x = b*e*b*e+4*b*e+e_0
            sage: y = a*c+d*c*b*f
            sage: x*(y+x) == x*y+x*x
            True

        TESTS:

        We compare against the multiplication in free algebras, which is
        implemented independently::

            sage: F.<x,y,z> = FreeAlgebra(GF(25,'t'))
            sage: A = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'))
            sage: pF = x+2*y-z+1
            sage: pA = sage_eval('x+2*y-z+1', A.gens_dict())
            sage: pA^5 == sage_eval(repr(pF^5), A.gens_dict())
            True
        """
        if self.data == NULL:
            return self
        cdef PathAlgebraElement right = other
        if right.data == NULL:
            return right
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2
        cdef path_term_t *T2
        cdef path_homog_poly_t *out_orig = NULL
        cdef path_homog_poly_t *out = NULL
        cdef path_homog_poly_t *nxt
        cdef path_term_t *P1start
        while H1 != NULL:
            H2 = right.data
            while H2 != NULL:
                sig_check()
                if H2.start == H1.end:
                    out = homog_poly_get_predecessor_of_component(out_orig, H1.start, H2.end)
                    if out == NULL:
                        if out_orig == NULL:
                            out_orig = homog_poly_create(H1.start, H2.end)
                        else:
                            if out_orig.start != H1.start or out_orig.end != H2.end:
                                nxt = out_orig
                                out_orig = homog_poly_create(H1.start, H2.end)
                                out_orig.nxt = nxt
                    else:
                        if out.nxt==NULL or out.nxt.start != H1.start or out.nxt.end != H2.end:
                            nxt = out.nxt
                            out.nxt = homog_poly_create(H1.start, H2.end)
                            out.nxt.nxt = nxt
                    T2 = H2.poly.lead
                    # now, either out==NULL, and we need to put the product
                    # into out_orig; or out!=NULL, and we need to put the
                    # product into out.nxt
                    if out == NULL:
                        P1start = out_orig.poly.lead
                        while T2 != NULL:
                            P1start = poly_iadd_lmul(out_orig.poly, <object>T2.coef, H1.poly,
                                                     T2.mon.path, self.cmp_terms, -1, 0, 0, P1start)
                            if P1start == H1.poly.lead:
                                P1start = out_orig.poly.lead
                            T2 = T2.nxt
                    else:
                        P1start = out.nxt.poly.lead
                        while T2 != NULL:
                            P1start = poly_iadd_lmul(out.nxt.poly, <object>T2.coef, H1.poly,
                                                     T2.mon.path, self.cmp_terms, -1, 0, 0, P1start)
                            if P1start == H1.poly.lead:
                                P1start = out.nxt.poly.lead
                            T2 = T2.nxt
                H2 = H2.nxt
            H1 = H1.nxt
        while out_orig != NULL and out_orig.poly.lead == NULL:
            tmp = out_orig.nxt
            sig_check()
            sig_free(out_orig.poly)
            sig_free(out_orig)
            out_orig = tmp
        if out_orig == NULL:
            return self._new_(NULL)
        tmp = out_orig
        while tmp.nxt != NULL:
            if tmp.nxt.poly.lead == NULL:
                sig_check()
                nxt = tmp.nxt.nxt
                sig_free(tmp.nxt.poly)
                sig_free(tmp.nxt)
                tmp.nxt = nxt
            else:
                tmp = tmp.nxt
        return self._new_(out_orig)

cpdef PathAlgebraElement path_algebra_element_unpickle(P, list data):
    """
    Auxiliary function for unpickling.

    EXAMPLES::

        sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15), order='negdeglex')
        sage: A.inject_variables()
        Defining e_0, e_1, e_2, a, b, c, d, e, f
        sage: X = a+2*b+3*c+5*e_0+3*e_2
        sage: loads(dumps(X)) == X    # indirect doctest
        True
    """
    cdef PathAlgebraElement out = P.element_class.__new__(P.element_class)
    out._parent = P
    order = P.order_string()
    if order == "negdegrevlex":
        out.cmp_terms = negdegrevlex
    elif order == "degrevlex":
        out.cmp_terms = degrevlex
    elif order == "negdeglex":
        out.cmp_terms = negdeglex
    elif order == "deglex":
        out.cmp_terms = deglex
    else:
        raise ValueError("unknown term order '{}'".format(order))
    out.data = homog_poly_unpickle(data)
    out._hash = -1
    return out
