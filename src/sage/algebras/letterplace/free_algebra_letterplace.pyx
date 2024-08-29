###############################################################################
#
#       Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
#
###############################################################################

"""
Free associative unital algebras, implemented via Singular's letterplace rings

AUTHOR:

- Simon King (2011-03-21): :issue:`7797`

With this implementation, Groebner bases out to a degree bound and
normal forms can be computed for twosided weighted homogeneous ideals
of free algebras. For now, all computations are restricted to weighted
homogeneous elements, i.e., other elements cannot be created by
arithmetic operations.

EXAMPLES::

    sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
    sage: F
    Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
    sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
    sage: I
    Twosided Ideal (x*y + y*z, x*x + x*y - y*x - y*y) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
    sage: x*(x*I.0-I.1*y+I.0*y)-I.1*y*z
    x*y*x*y + x*y*y*y - x*y*y*z + x*y*z*y + y*x*y*z + y*y*y*z
    sage: x^2*I.0-x*I.1*y+x*I.0*y-I.1*y*z in I
    True

The preceding containment test is based on the computation of Groebner
bases with degree bound::

    sage: I.groebner_basis(degbound=4)
    Twosided Ideal (x*y + y*z,
        x*x - y*x - y*y - y*z,
        y*y*y - y*y*z + y*z*y - y*z*z,
        y*y*x + y*y*z + y*z*x + y*z*z,
        y*y*z*y - y*y*z*z + y*z*z*y - y*z*z*z,
        y*z*y*y - y*z*y*z + y*z*z*y - y*z*z*z,
        y*y*z*x + y*y*z*z + y*z*z*x + y*z*z*z,
        y*z*y*x + y*z*y*z + y*z*z*x + y*z*z*z) of Free Associative Unital
        Algebra on 3 generators (x, y, z) over Rational Field

When reducing an element by `I`, the original generators are chosen::

    sage: (y*z*y*y).reduce(I)
    y*z*y*y

However, there is a method for computing the normal form of an
element, which is the same as reduction by the Groebner basis out to
the degree of that element::

    sage: (y*z*y*y).normal_form(I)
    y*z*y*z - y*z*z*y + y*z*z*z
    sage: (y*z*y*y).reduce(I.groebner_basis(4))
    y*z*y*z - y*z*z*y + y*z*z*z

The default term order derives from the degree reverse lexicographic
order on the commutative version of the free algebra::

    sage: F.commutative_ring().term_order()
    Degree reverse lexicographic term order

A different term order can be chosen, and of course may yield a
different normal form::

    sage: L.<a,b,c> = FreeAlgebra(QQ, implementation='letterplace', order='lex')
    sage: L.commutative_ring().term_order()
    Lexicographic term order
    sage: J = L*[a*b+b*c,a^2+a*b-b*c-c^2]*L
    sage: J.groebner_basis(4)
    Twosided Ideal (2*b*c*b - b*c*c + c*c*b,
        a*b + b*c,
        -a*c*c + 2*b*c*a + 2*b*c*c + c*c*a,
        a*c*c*b - 2*b*c*c*b + b*c*c*c,
        a*a - 2*b*c - c*c,
        a*c*c*a - 2*b*c*c*a - 4*b*c*c*c - c*c*c*c) of Free Associative Unital
        Algebra on 3 generators (a, b, c) over Rational Field
    sage: (b*c*b*b).normal_form(J)
    1/2*b*c*c*b - 1/2*c*c*b*b

Here is an example with degree weights::

    sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[1,2,3])
    sage: (x*y+z).degree()
    3

TESTS::

    sage: TestSuite(F).run()
    sage: TestSuite(L).run()
    sage: loads(dumps(F)) is F
    True

    sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
    sage: F.is_commutative()
    False
    sage: FreeAlgebra(QQ, implementation='letterplace', names=['x']).is_commutative()
    True

.. TODO::

    The computation of Groebner bases only works for global term
    orderings, and all elements must be weighted homogeneous with respect
    to positive integral degree weights. It is ongoing work in Singular to
    lift these restrictions.

    We support coercion from the letterplace wrapper to the corresponding
    generic implementation of a free algebra
    (:class:`~sage.algebras.free_algebra.FreeAlgebra_generic`), but there
    is no coercion in the opposite direction, since the generic
    implementation also comprises non-homogeneous elements.

    We also do not support coercion from a subalgebra, or between free
    algebras with different term orderings, yet.
"""
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.libs.singular.function import lib, singular_function
from sage.libs.singular.function cimport RingWrap
from sage.libs.singular.ring cimport singular_ring_delete, singular_ring_reference
from sage.categories.algebras import Algebras
from sage.rings.noncommutative_ideals import IdealMonoid_nc
from sage.rings.polynomial.plural cimport new_CRing
from sage.misc.cachefunc import cached_method

#####################
# Define some singular functions
lib("freegb.lib")
freeAlgebra = singular_function("freeAlgebra")

# unfortunately we cannot set Singular attributes for MPolynomialRing_libsingular
# Hence, we must constantly work around Letterplace's sanity checks,
# and cannot use the following library functions:
# set_letterplace_attributes = singular_function("setLetterplaceAttributes")
# lpMult = singular_function("lpMult")

#####################
# Auxiliary functions

cdef MPolynomialRing_libsingular make_letterplace_ring(base_ring, blocks):
    """
    Create a polynomial ring in block order.

    INPUT:

    - ``base_ring`` -- a multivariate polynomial ring
    - ``blocks`` -- the number of blocks to be formed

    OUTPUT:

    A multivariate polynomial ring in block order, all blocks
    isomorphic (as ordered rings) with the given ring, and the
    variable names of the `n`-th block (`n>0`) ending with
    ``"_%d"%n``.

    TESTS:

    Note that, since the algebras are cached, we need to choose
    a different base ring, since other doctests could have a
    side effect on the attained degree bound::

        sage: F.<x,y,z> = FreeAlgebra(GF(17), implementation='letterplace')
        sage: L.<a,b,c> = FreeAlgebra(GF(17), implementation='letterplace', order='lex')
        sage: F.set_degbound(4)
        sage: F.current_ring()  # indirect doctest
        Multivariate Polynomial Ring in x, y, z, x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3 over Finite Field of size 17
        sage: F.current_ring().term_order()
        Block term order with blocks:
        (Degree reverse lexicographic term order of length 3,
         Degree reverse lexicographic term order of length 3,
         Degree reverse lexicographic term order of length 3,
         Degree reverse lexicographic term order of length 3)
        sage: L.set_degbound(2)
        sage: L.current_ring().term_order()
        Block term order with blocks:
        (Lexicographic term order of length 3,
         Lexicographic term order of length 3)
    """
    T0 = base_ring.term_order()
    T = T0
    cdef i
    cdef tuple names0 = base_ring.variable_names()
    cdef list names = list(names0)
    for i in range(1, blocks):
        T += T0
        names.extend([x + '_' + str(i) for x in names0])
    return PolynomialRing(base_ring.base_ring(), names, order=T,
                          implementation='singular')


#####################
# The free algebra

cdef class FreeAlgebra_letterplace(Parent):
    """
    Finitely generated free algebra, with arithmetic restricted to weighted homogeneous elements.

    .. NOTE::

        The restriction to weighted homogeneous elements should be lifted
        as soon as the restriction to homogeneous elements is lifted in
        Singular's "Letterplace algebras".

    EXAMPLES::

        sage: K.<z> = GF(25)
        sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
        sage: F
        Free Associative Unital Algebra on 3 generators (a, b, c) over Finite Field in z of size 5^2
        sage: P = F.commutative_ring()
        sage: P
        Multivariate Polynomial Ring in a, b, c over Finite Field in z of size 5^2

    We can do arithmetic as usual, as long as we stay (weighted) homogeneous::

        sage: (z*a+(z+1)*b+2*c)^2
        (z + 3)*a*a + (2*z + 3)*a*b + (2*z)*a*c + (2*z + 3)*b*a + (3*z + 4)*b*b + (2*z + 2)*b*c + (2*z)*c*a + (2*z + 2)*c*b - c*c
        sage: a+1
        Traceback (most recent call last):
        ...
        ArithmeticError: can only add elements of the same weighted degree
    """
    # It is not really a free algebra over the given generators. Rather,
    # it is a free algebra over the commutative monoid generated by the given generators.
    def __init__(self, R, degrees=None):
        """
        INPUT:

        A multivariate polynomial ring of type :class:`~sage.rings.polynomial.multipolynomial_libsingular.MPolynomialRing_libsingular`.

        OUTPUT: the free associative version of the given commutative ring

        .. NOTE::

            One is supposed to use the :class:`FreeAlgebra` constructor,
            in order to use the cache.

        TESTS::

            sage: from sage.algebras.letterplace.free_algebra_letterplace import FreeAlgebra_letterplace
            sage: FreeAlgebra_letterplace(QQ['x','y'])
            Free Associative Unital Algebra on 2 generators (x, y) over Rational Field
            sage: FreeAlgebra_letterplace(QQ['x'])
            Traceback (most recent call last):
            ...
            TypeError: a letterplace algebra must be provided by a polynomial ring of type <... 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

        ::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
            sage: TestSuite(F).run()
        """
        if not isinstance(R, MPolynomialRing_libsingular):
            raise TypeError("a letterplace algebra must be provided by a polynomial ring of type %s" % MPolynomialRing_libsingular)

        self._ngens = R.ngens()
        if degrees is None:
            varnames = R.variable_names()
            self._nb_slackvars = 0
        else:
            varnames = R.variable_names()[:-1]
            self._nb_slackvars = 1
        number_of_vars = self._ngens - self._nb_slackvars

        base_ring = R.base_ring()
        if number_of_vars <= 1:
            cat = Algebras(base_ring).Commutative()
        else:
            cat = Algebras(base_ring)
        Parent.__init__(self, base=base_ring, names=varnames,
                        normalize=False, category=cat)

        self._commutative_ring = R
        self._current_ring = make_letterplace_ring(R, 1)
        self._degbound = 1
        if degrees is None:
            self._degrees = tuple([int(1)] * self._ngens)
        else:
            if (not isinstance(degrees, (tuple, list))) \
                    or len(degrees) != number_of_vars \
                    or any(i <= 0 for i in degrees):
                raise TypeError("the generator degrees must be given by a list or tuple of %d positive integers" % (self._ngens - 1))
            self._degrees = tuple([int(i) for i in degrees])
            self.set_degbound(max(self._degrees))
        self._populate_coercion_lists_(coerce_list=[base_ring])

    def __reduce__(self):
        """
        TESTS::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
            sage: loads(dumps(F)) is F    # indirect doctest
            True
        """
        from sage.algebras.free_algebra import FreeAlgebra
        if self._nb_slackvars == 0:
            return FreeAlgebra, (self._commutative_ring,)
        return FreeAlgebra, (self._commutative_ring, None, None, None,
                             None, None, None, None, self._degrees)

    # Small methods
    def ngens(self):
        """
        Return the number of generators.

        EXAMPLES::

            sage: F.<a,b,c> = FreeAlgebra(QQ, implementation='letterplace')
            sage: F.ngens()
            3
        """
        return self._ngens - self._nb_slackvars

    @cached_method
    def gen(self, i):
        """
        Return the `i`-th generator.

        INPUT:

        - ``i`` -- integer

        OUTPUT: the generator with index `i`

        EXAMPLES::

            sage: F.<a,b,c> = FreeAlgebra(QQ, implementation='letterplace')
            sage: F.1 is F.1  # indirect doctest
            True
            sage: F.gen(2)
            c
        """
        if i >= self._ngens - self._nb_slackvars:
            raise ValueError("this free algebra only has %d generators" % (self._ngens - self._nb_slackvars))
        deg = self._degrees[i]
        p = self._current_ring.gen(i)
        cdef int n
        cdef int j = self._ngens - 1
        for n in range(1, deg):
            j += self._ngens
            p *= self._current_ring.gen(j)
        return FreeAlgebraElement_letterplace(self, p)

    def gens(self):
        """
        Return the tuple of generators.

        EXAMPLES::

            sage: F.<a,b,c> = FreeAlgebra(QQ, implementation='letterplace')
            sage: F.gens()
            (a, b, c)
        """
        return tuple(self.gen(i) for i in range(self.ngens()))

    def current_ring(self):
        """
        Return the commutative ring that is used to emulate
        the non-commutative multiplication out to the current degree.

        EXAMPLES::

            sage: F.<a,b,c> = FreeAlgebra(QQ, implementation='letterplace')
            sage: F.current_ring()
            Multivariate Polynomial Ring in a, b, c over Rational Field
            sage: a*b
            a*b
            sage: F.current_ring()
            Multivariate Polynomial Ring in a, b, c, a_1, b_1, c_1 over Rational Field
            sage: F.set_degbound(3)
            sage: F.current_ring()
            Multivariate Polynomial Ring in a, b, c, a_1, b_1, c_1, a_2, b_2, c_2 over Rational Field
        """
        return self._current_ring

    def commutative_ring(self):
        """
        Return the commutative version of this free algebra.

        .. NOTE::

            This commutative ring is used as a unique key of the free algebra.

        EXAMPLES::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace')
            sage: F
            Free Associative Unital Algebra on 3 generators (a, b, c) over Finite Field in z of size 5^2
            sage: F.commutative_ring()
            Multivariate Polynomial Ring in a, b, c over Finite Field in z of size 5^2
            sage: FreeAlgebra(F.commutative_ring()) is F
            True
        """
        return self._commutative_ring

    def term_order_of_block(self):
        """
        Return the term order that is used for the commutative version of this free algebra.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: F.term_order_of_block()
            Degree reverse lexicographic term order
            sage: L.<a,b,c> = FreeAlgebra(QQ, implementation='letterplace',order='lex')
            sage: L.term_order_of_block()
            Lexicographic term order
        """
        return self._commutative_ring.term_order()

    def generator_degrees(self):
        return self._degrees

    # Some basic properties of this ring
    def is_field(self, proof=True):
        """
        Tell whether this free algebra is a field.

        .. NOTE::

            This would only be the case in the degenerate case of no generators.
            But such an example cannot be constructed in this implementation.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: F.is_field()
            False
        """
        return (not (self._ngens - self._nb_slackvars)) and self._base.is_field(proof=proof)

    def _repr_(self):
        """
        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: F     # indirect doctest
            Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field

        The degree weights are not part of the string representation::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
            sage: F
            Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
        """
        return "Free Associative Unital Algebra on %d generators %s over %s" % (self._ngens - self._nb_slackvars, self.gens(), self._base)

    def _latex_(self):
        r"""
        Representation of this free algebra in LaTeX.

        EXAMPLES::

            sage: F.<bla,alpha,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[1,2,3])
            sage: latex(F)
            \Bold{Q}\langle \mathit{bla}, \alpha, z\rangle
        """
        from sage.misc.latex import latex
        return "%s\\langle %s\\rangle" % (latex(self.base_ring()), ', '.join(self.latex_variable_names()))

    def degbound(self):
        """
        Return the degree bound that is currently used.

        .. NOTE::

            When multiplying two elements of this free algebra, the degree
            bound will be dynamically adapted. It can also be set by
            :meth:`set_degbound`.

        EXAMPLES:

        In order to avoid we get a free algebras from the cache that
        was created in another doctest and has a different degree
        bound, we choose a base ring that does not appear in other tests::

            sage: F.<x,y,z> = FreeAlgebra(ZZ, implementation='letterplace')
            sage: F.degbound()
            1
            sage: x*y
            x*y
            sage: F.degbound()
            2
            sage: F.set_degbound(4)
            sage: F.degbound()
            4
        """
        return self._degbound

    def set_degbound(self, d):
        """
        Increase the degree bound that is currently in place.

        .. NOTE::

            The degree bound cannot be decreased.

        EXAMPLES:

        In order to avoid we get a free algebras from the cache that
        was created in another doctest and has a different degree
        bound, we choose a base ring that does not appear in other tests::

            sage: F.<x,y,z> = FreeAlgebra(GF(251), implementation='letterplace')
            sage: F.degbound()
            1
            sage: x*y
            x*y
            sage: F.degbound()
            2
            sage: F.set_degbound(4)
            sage: F.degbound()
            4
            sage: F.set_degbound(2)
            sage: F.degbound()
            4
        """
        if d <= self._degbound:
            return
        self._degbound = d
        self._current_ring = make_letterplace_ring(self._commutative_ring, d)

#    def base_extend(self, R):
#        if self._base.has_coerce_map_from(R):
#            return self

    ################################################
    # Ideals

    def _ideal_class_(self, n=0):
        """
        Return the class :class:`~sage.algebras.letterplace.letterplace_ideal.LetterplaceIdeal`.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = [x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: I
            Right Ideal (x*y + y*z, x*x + x*y - y*x - y*y) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
            sage: type(I) is F._ideal_class_()
            True
        """
        from sage.algebras.letterplace.letterplace_ideal import LetterplaceIdeal
        return LetterplaceIdeal

    def ideal_monoid(self):
        """
        Return the monoid of ideals of this free algebra.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(GF(2), implementation='letterplace')
            sage: F.ideal_monoid()
            Monoid of ideals of Free Associative Unital Algebra on 2 generators (x, y) over Finite Field of size 2
            sage: F.ideal_monoid() is F.ideal_monoid()
            True
        """
        if self.__monoid is None:
            self.__monoid = IdealMonoid_nc(self)
        return self.__monoid

    # Auxiliary methods
    cdef str exponents_to_string(self, E):
        """
        This auxiliary method is used for the string representation of elements of this free algebra.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(GF(2), implementation='letterplace')
            sage: x*y*x*z   # indirect doctest
            x*y*x*z

        It should be possible to use the letterplace algebra to implement the
        free algebra generated by the elements of a finitely generated free abelian
        monoid. However, we cannot use it, yet. So, for now, we raise an error::

            sage: from sage.algebras.letterplace.free_algebra_element_letterplace import FreeAlgebraElement_letterplace
            sage: P = F.commutative_ring()
            sage: FreeAlgebraElement_letterplace(F, P.0*P.1^2+P.1^3) # indirect doctest
            <repr(<sage.algebras.letterplace.free_algebra_element_letterplace.FreeAlgebraElement_letterplace at 0x...>) failed: NotImplementedError:
              Apparently you tried to view the letterplace algebra with
              shift-multiplication as the free algebra over a finitely
              generated free abelian monoid.
              In principle, this is correct, but it is not implemented, yet.>
        """
        cdef int ngens = self._ngens
        cdef int nblocks = len(E) // ngens
        cdef int i, j, base, exp, var_ind
        cdef list out = []
        cdef list tmp
        for i from 0<=i<nblocks:
            # warning : the loop index is changed during the loop
            base = i*ngens
            tmp = [(j, E[base + j]) for j in range(ngens) if E[base + j]]
            if not tmp:
                continue
            var_ind, exp = tmp[0]
            if len(tmp) > 1 or exp > 1:
                raise NotImplementedError("\n  Apparently you tried to view the letterplace algebra with\n  shift-multiplication as the free algebra over a finitely\n  generated free abelian monoid.\n  In principle, this is correct, but it is not implemented, yet.")
            out.append(self._names[var_ind])
            i += self._degrees[var_ind] - 1
        return '*'.join(out)

    # Auxiliary methods
    cdef str exponents_to_latex(self, E):
        r"""
        This auxiliary method is used for the representation of elements of this free algebra as a latex string.

        EXAMPLES::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K, implementation='letterplace', degrees=[1,2,3])
            sage: -(a*b*(z+1)-c)^2
            (2*z + 1)*a*b*a*b + (z + 1)*a*b*c + (z + 1)*c*a*b - c*c
            sage: latex(-(a*b*(z+1)-c)^2)     # indirect doctest
            \left(2 z + 1\right) a b a b + \left(z + 1\right) a b c + \left(z + 1\right) c a b - c c
        """
        cdef int ngens = self._ngens
        cdef int nblocks = len(E) // ngens
        cdef int i, j, base, exp, var_ind
        cdef list out = []
        cdef list tmp
        cdef list names = self.latex_variable_names()
        for i from 0<=i<nblocks:
            # warning : the loop index is changed during the loop
            base = i*ngens
            tmp = [(j, E[base+j]) for j in range(ngens) if E[base+j]]
            if not tmp:
                continue
            var_ind, exp = tmp[0]
            if len(tmp) > 1 or exp > 1:
                raise NotImplementedError("\n  Apparently you tried to view the letterplace algebra with\n  shift-multiplication as the free algebra over a finitely\n  generated free abelian monoid.\n  In principle, this is correct, but it is not implemented, yet.")
            out.append(names[var_ind])
            i += (self._degrees[var_ind] - 1)
        return ' '.join(out)

    def _reductor_(self, g, d):
        """
        Return a commutative ideal that can be used to compute the normal
        form of a free algebra element of a given degree.

        INPUT:

        - ``g`` -- list of elements of this free algebra
        - ``d`` -- integer

        OUTPUT:

        An ideal such that reduction of a letterplace polynomial by that ideal corresponds
        to reduction of an element of degree at most ``d`` by ``g``.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: p = y*x*y + y*y*y + y*z*y - y*z*z
            sage: p.reduce(I)
            y*y*y - y*y*z + y*z*y - y*z*z
            sage: G = F._reductor_(I.gens(),3); G
            Ideal (x*y_1 + y*z_1, x_1*y_2 + y_1*z_2, x*x_1 + x*y_1 - y*x_1 - y*y_1, x_1*x_2 + x_1*y_2 - y_1*x_2 - y_1*y_2) of Multivariate Polynomial Ring in x, y, z, x_1, y_1, z_1, x_2, y_2, z_2... over Rational Field

        We do not use the usual reduction method for polynomials in
        Sage, since it does the reductions in a different order
        compared to Singular. Therefore, we call the original Singular
        reduction method, and prevent a warning message by asserting
        that `G` is a Groebner basis. ::

            sage: from sage.libs.singular.function import singular_function
            sage: poly_reduce = singular_function("NF")
            sage: q = poly_reduce(p.letterplace_polynomial(), G, ring=F.current_ring(), attributes={G:{"isSB":1}}); q
            y*y_1*y_2 - y*y_1*z_2 + y*z_1*y_2 - y*z_1*z_2
            sage: p.reduce(I).letterplace_polynomial() == q
            True
        """
        cdef list out = []
        C = self.current_ring()
        cdef FreeAlgebraElement_letterplace x
        ngens = self._ngens
        cdef list G = [C(x._poly) for x in g]
        from sage.groups.perm_gps.permgroup_named import CyclicPermutationGroup
        CG = CyclicPermutationGroup(C.ngens())
        for y in G:
            out.extend([y] + [y * CG[ngens * (n + 1)]
                              for n in range(d - y.degree())])
        return C.ideal(out)

    ###########################
    # Coercion
    cpdef _coerce_map_from_(self, S):
        """
        A ring ``R`` coerces into ``self``, if:

        - it coerces into the current polynomial ring, or
        - it is a free graded algebra in letterplace implementation,
          the generator names of ``R`` are a proper subset of the
          generator names of ``self``, the degrees of equally named
          generators are equal, and the base ring of ``R`` coerces
          into the base ring of ``self``.

        TESTS:

        Coercion from the base ring::

            sage: F.<x,y,z> = FreeAlgebra(GF(5), implementation='letterplace')
            sage: 5 == F.zero()    # indirect doctest
            True

        Coercion from another free graded algebra::

            sage: F.<t,y,z> = FreeAlgebra(ZZ, implementation='letterplace', degrees=[4,2,3])
            sage: G = FreeAlgebra(GF(5), implementation='letterplace', names=['x','y','z','t'], degrees=[1,2,3,4])
            sage: t*G.0       # indirect doctest
            t*x
        """
        if self == S or self._current_ring.has_coerce_map_from(S):
            return True
        cdef int i
        # Do we have another letterplace algebra?
        if not isinstance(S, FreeAlgebra_letterplace):
            return False
        # Do the base rings coerce?
        if not self.base_ring().has_coerce_map_from(S.base_ring()):
            return False
        # Do the names match?
        cdef tuple degs, Sdegs, names, Snames
        names = self.variable_names()
        Snames = S.variable_names()
        if not set(names).issuperset(Snames):
            return False
        # Do the degrees match
        degs = self._degrees
        Sdegs = (<FreeAlgebra_letterplace>S)._degrees
        for i in range(S.ngens()):
            if degs[names.index(Snames[i])] != Sdegs[i]:
                return False
        return True

    def _an_element_(self):
        """
        Return an element.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: F.an_element()   # indirect doctest
            x
        """
        return FreeAlgebraElement_letterplace(self, self._current_ring.an_element(), check=False)

#    def random_element(self, degree=2, terms=5):
#        """
#        Return a random element of a given degree and with a given number of terms.
#
#        INPUT:
#
#        - ``degree`` -- the maximal degree of the output (default: 2)
#        - ``terms`` -- the maximal number of terms of the output (default: 5)
#
#        NOTE:
#
#        This method is currently not useful at all.
#
#        Not tested.
#        """
#        self.set_degbound(degree)
#        while(1):
#            p = self._current_ring.random_element(degree=degree,terms=terms)
#            if p.is_homogeneous():
#                break
#        return FreeAlgebraElement_letterplace(self, p, check=False)

    def _from_dict_(self, dict D, check=True):
        """
        Create an element from a dictionary.

        INPUT:

        - ``D`` -- dictionary; keys: tuples of exponents, values:
          the coefficients of the corresponding monomial
          in the to-be-created element
        - ``check`` -- boolean (default: ``True``);
          this is forwarded to the initialisation of
          :class:`~sage.algebras.letterplace.free_algebra_element_letterplace.FreeAlgebraElement_letterplace`

        TESTS:

        This method applied to the dictionary of any element must
        return the same element. This must hold true even if the
        underlying letterplace ring has been extended in the meantime.
        ::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: p = 3*x*y+2*z^2
            sage: F.set_degbound(10)
            sage: p == F._from_dict_(dict(p))
            True

        For the empty dictionary, zero is returned::

            sage: F._from_dict_({})
            0
        """
        if not D:
            return self.zero()
        cdef Py_ssize_t l
        for e in D:
            l = len(e)
            break
        cdef dict out = {}
        self.set_degbound(l // self._ngens)
        cdef Py_ssize_t n = self._current_ring.ngens()
        for e, c in D.iteritems():
            out[tuple(e) + (0,) * (n - l)] = c
        return FreeAlgebraElement_letterplace(self, self._current_ring(out),
                                              check=check)

    def _element_constructor_(self, x):
        """
        Return an element of this free algebra.

        INPUT:

        An element of a free algebra with a proper subset of generator
        names, or anything that can be interpreted in the polynomial
        ring that is used to implement the letterplace algebra out to
        the current degree bound, or a string that can be interpreted
        as an expression in the algebra (provided that the
        coefficients are numerical).

        EXAMPLES::

            sage: F.<t,y,z> = FreeAlgebra(ZZ, implementation='letterplace', degrees=[4,2,3])

        Conversion of a number::

            sage: F(3)
            3

        Interpretation of a string as an algebra element::

            sage: F('t*y+3*z^2')
            t*y + 3*z*z

        Conversion from the currently underlying polynomial ring::

            sage: F.set_degbound(3)
            sage: P = F.current_ring()
            sage: F(P.0*P.7*P.11*P.15*P.17*P.23 - 2*P.2*P.7*P.11*P.14*P.19*P.23)
            t*y - 2*z*z

        Conversion from a graded sub-algebra::

            sage: G = FreeAlgebra(GF(5), implementation='letterplace', names=['x','y','z','t'], degrees=[1,2,3,4])
            sage: G(t*y + 2*y^3 - 4*z^2)   # indirect doctest
            (2)*y*y*y + z*z + t*y
        """
        if isinstance(x, str):
            from sage.misc.sage_eval import sage_eval
            return sage_eval(x, locals=self.gens_dict())
        try:
            P = x.parent()
        except AttributeError:
            P = None
        if P is self:
            (<FreeAlgebraElement_letterplace>x)._poly = self._current_ring((<FreeAlgebraElement_letterplace>x)._poly)
            return x
        if isinstance(P, FreeAlgebra_letterplace):
            self.set_degbound(P.degbound())
            Ppoly = (<FreeAlgebra_letterplace>P)._current_ring
            Gens = self._current_ring.gens()
            Names = self._current_ring.variable_names()
            PNames = list(Ppoly.variable_names())
            # translate the slack variables
            PNames[P.ngens(): len(PNames): P.ngens() + 1] = list(Names[self.ngens(): len(Names): self.ngens() + 1])[:P.degbound()]
            x = Ppoly.hom([Gens[Names.index(asdf)] for asdf in PNames])(x.letterplace_polynomial())
        return FreeAlgebraElement_letterplace(self, self._current_ring(x))


cdef class FreeAlgebra_letterplace_libsingular():
    """
    Internally used wrapper around a Singular Letterplace polynomial ring.
    """

    def __cinit__(self, MPolynomialRing_libsingular commutative_ring,
                  int degbound):
        cdef RingWrap rw = freeAlgebra(commutative_ring, degbound)
        self._lp_ring = singular_ring_reference(rw._ring)
        # `_lp_ring` viewed as `MPolynomialRing_libsingular` with additional
        # letterplace attributes set (for internal use only)
        self._lp_ring_internal = new_CRing(rw, commutative_ring.base_ring())
        self._commutative_ring = commutative_ring

    def __init__(self, commutative_ring, degbound):
        self._ngens = commutative_ring.ngens() * degbound

    def __dealloc__(self):
        r"""
        Carefully deallocate the ring, without changing ``currRing``
        (since this method can be at unpredictable times due to garbage
        collection).
        """
        singular_ring_delete(self._lp_ring)
