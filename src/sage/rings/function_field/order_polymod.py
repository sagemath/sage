# sage.doctest: needs sage.rings.function_field
r"""
Orders of function fields: extension
"""

# ****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2011      Maarten Derickx <m.derickx.student@gmail.com>
#                     2011      Julian Rueth <julian.rueth@gmail.com>
#                     2017-2020 Kwankyu Lee
#                     2019      Brent Baccala
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.functions import lcm
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from .ideal import FunctionFieldIdeal
from .ideal_polymod import (
    FunctionFieldIdeal_polymod,
    FunctionFieldIdeal_global,
    FunctionFieldIdealInfinite_polymod
)
from .order import FunctionFieldMaximalOrder, FunctionFieldMaximalOrderInfinite


class FunctionFieldMaximalOrder_polymod(FunctionFieldMaximalOrder):
    """
    Maximal orders of extensions of function fields.
    """

    def __init__(self, field, ideal_class=FunctionFieldIdeal_polymod):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: TestSuite(O).run()
        """
        FunctionFieldMaximalOrder.__init__(self, field, ideal_class)

        from sage.modules.free_module_element import vector
        from .function_field_polymod import FunctionField_integral

        if isinstance(field, FunctionField_integral):
            basis = field._maximal_order_basis()
        else:
            model, from_model, to_model = field.monic_integral_model('z')
            basis = [from_model(g) for g in model._maximal_order_basis()]

        V, fr, to = field.vector_space()
        R = field.base_field().maximal_order()

        # This module is over R, but linear algebra over R (MaximalOrder)
        # is not well supported in Sage. So we keep it as a vector space
        # over rational function field.
        self._module = V.span_of_basis([to(b) for b in basis])
        self._module_base_ring = R
        self._basis = tuple(basis)
        self._from_module = fr
        self._to_module = to

        # multiplication table (lower triangular)
        n = len(basis)
        self._mtable = []
        for i in range(n):
            row = [self._coordinate_vector(basis[i] * basis[j])
                   for j in range(n)]
            self._mtable.append(row)

        zero = vector(R._ring, n * [0])

        def mul_vecs(f, g):
            s = zero
            for i in range(n):
                if f[i].is_zero():
                    continue
                for j in range(n):
                    if g[j].is_zero():
                        continue
                    s += f[i] * g[j] * self._mtable[i][j]
            return s
        self._mul_vecs = mul_vecs

        # We prepare for using Kummer's theorem to decompose primes. Note
        # that Kummer's theorem applies to most places. Here we find
        # places for which the theorem does not apply.

        # this element is integral over k[x] and a generator of the field.
        for gen in basis:
            phi = gen.minimal_polynomial()
            if phi.degree() == n:
                break

        assert phi.degree() == n

        gen_vec = self._coordinate_vector(gen)
        g = gen_vec.parent().gen(0) # x
        gen_vec_pow = [g]
        for i in range(n):
            g = mul_vecs(g, gen_vec)
            gen_vec_pow.append(g)

        # find places where {1,gen,...,gen^(n-1)} is not integral basis
        W = V.span_of_basis([to(gen ** i) for i in range(phi.degree())])

        supp = []
        for g in basis:
            for c in W.coordinate_vector(to(g), check=False):
                if not c.is_zero():
                    supp += [f for f,_ in c.denominator().factor()]
        supp = set(supp)

        self._kummer_gen = gen
        self._kummer_gen_vec_pow = gen_vec_pow
        self._kummer_polynomial = phi
        self._kummer_places = supp

    def _element_constructor_(self, f):
        """
        Construct an element of this order from ``f``.

        INPUT:

        - ``f`` -- element convertible to the function field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x*Y + x^2 + 1)
            sage: O = L.maximal_order()
            sage: y in O
            True
            sage: 1/y in O
            False
            sage: x in O
            True
            sage: 1/x in O
            False
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: 1 in O
            True
            sage: y in O
            False
            sage: x*y in O
            True
            sage: x^2*y in O
            True
        """
        F = self.function_field()
        f = F(f)
        # check if f is in this order
        if not all(e in self._module_base_ring for e in self.coordinate_vector(f)):
            raise TypeError( "{} is not an element of {}".format(f, self) )

        return f

    def ideal_with_gens_over_base(self, gens):
        """
        Return the fractional ideal with basis ``gens`` over the
        maximal order of the base field.

        INPUT:

        - ``gens`` -- list of elements that generates the ideal over the
          maximal order of the base field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order(); O
            Maximal order of Function field in y defined by y^2 + 6*x^3 + 6
            sage: I = O.ideal_with_gens_over_base([1, y]);  I
            Ideal (1) of Maximal order of Function field in y defined by y^2 + 6*x^3 + 6
            sage: I.module()
            Free module of degree 2 and rank 2 over
             Maximal order of Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [1 0]
            [0 1]

        There is no check if the resulting object is really an ideal::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([y]); I
            Ideal (y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: y in I
            True
            sage: y^2 in I
            False
        """
        return self._ideal_from_vectors([self.coordinate_vector(g) for g in gens])

    def _ideal_from_vectors(self, vecs):
        """
        Return an ideal generated as a module by vectors over rational function
        field.

        INPUT:

        - ``vec`` -- list of vectors

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: v1 = O.coordinate_vector(x^3 + 1)
            sage: v2 = O.coordinate_vector(y)
            sage: v1
            (x^3 + 1, 0)
            sage: v2
            (0, 1)
            sage: O._ideal_from_vectors([v1, v2])
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + 6*x^3 + 6
        """
        d = lcm([v.denominator() for v in vecs])
        vecs = [[(d*c).numerator() for c in v] for v in vecs]
        return self._ideal_from_vectors_and_denominator(vecs, d, check=False)

    def _ideal_from_vectors_and_denominator(self, vecs, d=1, check=True):
        """
        Return an ideal generated as a module by vectors divided by ``d`` over
        the polynomial ring underlying the rational function field.

        INPUT:

        - ``vec`` -- list of vectors over the polynomial ring

        - ``d`` -- (default: 1) a nonzero element of the polynomial ring

        - ``check`` -- boolean (default: ``True``); if ``True``, compute the
          real denominator of the vectors, possibly different from ``d``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y^2)
            sage: m = I.basis_matrix()
            sage: v1 = m[0]
            sage: v2 = m[1]
            sage: v1
            (x^3 + 1, 0)
            sage: v2
            (0, x^3 + 1)
            sage: O._ideal_from_vectors([v1, v2])  # indirect doctest
            Ideal (x^3 + 1) of Maximal order of Function field in y
            defined by y^2 + 6*x^3 + 6
        """
        from sage.matrix.constructor import matrix
        from .hermite_form_polynomial import reversed_hermite_form

        R = self._module_base_ring._ring

        d = R(d) # make it sure that d is in the polynomial ring

        if check and not d.is_one(): # check if d is true denominator
            M = []
            g = d
            for v in vecs:
                for c in v:
                    g = g.gcd(c)
                    if g.is_one():
                        break
                else:
                    M += list(v)
                    continue # for v in vecs
                mat = matrix(R, vecs)
                break
            else:
                d = d // g
                mat = matrix(R, len(vecs), [c // g for c in M])
        else:
            mat = matrix(R, vecs)

        # IMPORTANT: make it sure that pivot polynomials monic
        # so that we get a unique hnf. Here the hermite form
        # algorithm also makes the pivots monic.

        # compute the reversed hermite form with zero rows deleted
        reversed_hermite_form(mat)
        i = 0
        while i < mat.nrows() and mat.row(i).is_zero():
            i += 1
        hnf = mat[i:] # remove zero rows

        return self.ideal_monoid().element_class(self, hnf, d)

    def ideal(self, *gens, **kwargs):
        """
        Return the fractional ideal generated by the elements in ``gens``.

        INPUT:

        - ``gens`` -- list of generators

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2 - 4)
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: S = L.maximal_order()
            sage: S.ideal(1/y)
            Ideal ((1/(x^3 + 1))*y) of Maximal order of Function field
            in y defined by y^2 + 6*x^3 + 6
            sage: I2 = S.ideal(x^2 - 4); I2
            Ideal (x^2 + 3) of Maximal order of Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 == S.ideal(I)
            True

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2 - 4)
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: S = L.maximal_order()
            sage: S.ideal(1/y)
            Ideal ((1/(x^3 + 1))*y) of
             Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: I2 = S.ideal(x^2-4); I2
            Ideal (x^2 - 4) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: I2 == S.ideal(I)
            True
        """
        if len(gens) == 1:
            gens = gens[0]
            if not isinstance(gens, (list, tuple)):
                if isinstance(gens, FunctionFieldIdeal):
                    gens = gens.gens()
                else:
                    gens = (gens,)
        F = self.function_field()
        mgens = [b*F(g) for g in gens for b in self.basis()]
        return self.ideal_with_gens_over_base(mgens)

    def polynomial(self):
        """
        Return the defining polynomial of the function field of which this is an order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.polynomial()
            y^4 + x*y + 4*x + 1

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.polynomial()
            y^4 + x*y + 4*x + 1
        """
        return self._field.polynomial()

    def basis(self):
        """
        Return a basis of the order over the maximal order of the base function
        field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.basis()
            (1, y, y^2, y^3)

            sage: K.<x> = FunctionField(QQ)
            sage: R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18)
            sage: O = F.maximal_order()
            sage: O.basis()
            (1, 1/x^4*y, 1/x^9*y^2, 1/x^13*y^3)
        """
        return self._basis

    def gen(self, n=0):
        """
        Return the ``n``-th generator of the order.

        The basis elements of the order are generators.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: O = L.maximal_order()
            sage: O.gen()
            1
            sage: O.gen(1)
            y
            sage: O.gen(2)
            (1/(x^3 + x^2 + x))*y^2
            sage: O.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: there are only 3 generators
        """
        if not ( n >= 0 and n < self.ngens() ):
            raise IndexError("there are only {} generators".format(self.ngens()))

        return self._basis[n]

    def ngens(self):
        """
        Return the number of generators of the order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order()
            sage: Oinf.ngens()
            3
        """
        return len(self._basis)

    def free_module(self):
        """
        Return the free module formed by the basis over the maximal order of the base field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O.free_module()
            Free module of degree 4 and rank 4 over
             Maximal order of Rational function field in x over Finite Field of size 7
            User basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        return self._module.change_ring(self._module_base_ring)

    def coordinate_vector(self, e):
        """
        Return the coordinates of ``e`` with respect to the basis of this order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O.coordinate_vector(y)
            (0, 1, 0, 0)
            sage: O.coordinate_vector(x*y)
            (0, x, 0, 0)

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: f = (x + y)^3
            sage: O.coordinate_vector(f)
            (x^3, 3*x^2, 3*x, 1)
        """
        return self._module.coordinate_vector(self._to_module(e))

    def _coordinate_vector(self, e):
        """
        Return the coordinate vector of ``e`` with respect to the basis
        of the order.

        Assuming ``e`` is in the maximal order, the coordinates are given
        as univariate polynomials in the underlying ring of the maximal
        order of the rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O._coordinate_vector(y)
            (0, 1, 0, 0)
            sage: O._coordinate_vector(x*y)
            (0, x, 0, 0)
        """
        from sage.modules.free_module_element import vector

        v = self._module.coordinate_vector(self._to_module(e), check=False)
        return vector([c.numerator() for c in v])

    @cached_method
    def different(self):
        """
        Return the different ideal of the function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O.different()
            Ideal (y^3 + 2*x)
            of Maximal order of Function field in y defined by y^4 + x*y + 4*x + 1
        """
        return ~self.codifferent()

    @cached_method
    def codifferent(self):
        """
        Return the codifferent ideal of the function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O.codifferent()
            Ideal (1, (1/(x^4 + 4*x^3 + 3*x^2 + 6*x + 4))*y^3
            + ((5*x^3 + 6*x^2 + x + 6)/(x^4 + 4*x^3 + 3*x^2 + 6*x + 4))*y^2
            + ((x^3 + 2*x^2 + 2*x + 2)/(x^4 + 4*x^3 + 3*x^2 + 6*x + 4))*y
            + 6*x/(x^4 + 4*x^3 + 3*x^2 + 6*x + 4)) of Maximal order of Function field
            in y defined by y^4 + x*y + 4*x + 1
        """
        T = self._codifferent_matrix()
        return self._ideal_from_vectors(T.inverse().columns())

    @cached_method
    def _codifferent_matrix(self):
        """
        Return the matrix `T` defined in Proposition 4.8.19 of [Coh1993]_.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O._codifferent_matrix()
            [      4       0       0     4*x]
            [      0       0     4*x 5*x + 3]
            [      0     4*x 5*x + 3       0]
            [    4*x 5*x + 3       0   3*x^2]
        """
        from sage.matrix.constructor import matrix

        rows = []
        for u in self.basis():
            row = []
            for v in self.basis():
                row.append((u*v).trace())
            rows.append(row)
        T = matrix(rows)
        return T

    @cached_method
    def decomposition(self, ideal):
        """
        Return the decomposition of the prime ideal.

        INPUT:

        - ``ideal`` -- prime ideal of the base maximal order

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: o = K.maximal_order()
            sage: O = F.maximal_order()
            sage: p = o.ideal(x + 1)
            sage: O.decomposition(p)
            [(Ideal (x + 1, y + 1) of Maximal order
             of Function field in y defined by y^3 + x^6 + x^4 + x^2, 1, 1),
             (Ideal (x + 1, (1/(x^3 + x^2 + x))*y^2 + y + 1) of Maximal order
             of Function field in y defined by y^3 + x^6 + x^4 + x^2, 2, 1)]

        ALGORITHM:

        In principle, we're trying to compute a primary decomposition
        of the extension of ``ideal`` in ``self`` (an order, and therefore
        a ring). However, while we have primary decomposition methods
        for polynomial rings, we lack any such method for an order.
        Therefore, we construct ``self`` mod ``ideal`` as a
        finite-dimensional algebra, a construct for which we do
        support primary decomposition.

        See :issue:`28094` and https://github.com/sagemath/sage/files/10659303/decomposition.pdf.gz

        .. TODO::

            Use Kummer's theorem to shortcut this code if possible, like as
            done in :meth:`FunctionFieldMaximalOrder_global.decomposition()`
        """
        from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra import FiniteDimensionalAlgebra
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector

        F = self.function_field()
        n = F.degree()

        # Base rational function field
        K = self.function_field().base_field()

        # Univariate polynomial ring isomorphic to the maximal order of K
        o = PolynomialRing(K.constant_field(), K.gen())

        # Prime ideal in o defined by the generator of ideal in the maximal
        # order of K
        p = o(ideal.gen().numerator())

        # Residue field k = o mod p
        k = o.quo(p)

        # Given an element of the function field expressed as a K-vector times
        # the basis of this order, construct the n n-by-n matrices that show
        # how to multiply by each of the basis elements.
        matrices = [matrix(o, [self.coordinate_vector(b1*b2) for b1 in self.basis()])
                    for b2 in self.basis()]

        # Let O denote the maximal order self. When reduced modulo p,
        # matrices_reduced give the multiplication matrices used to form the
        # algebra O mod pO.
        matrices_reduced = [M.mod(p) for M in matrices]
        cat = CommutativeAlgebras(k).FiniteDimensional().WithBasis()
        A = FiniteDimensionalAlgebra(k, matrices_reduced, category=cat)

        # Each prime ideal of the algebra A corresponds to a prime ideal of O,
        # and since the algebra is an Artinian ring, all of its prime ideals
        # are maximal [stacks 00JA]. Thus, we find all of our factors by
        # iterating over the algebra's maximal ideals.
        factors = []
        for q in A.maximal_ideals():
            if q == A.zero_ideal():
                # The zero ideal is the unique maximal ideal, which means that
                # A is a field, and the ideal itself is a prime ideal.
                P = self.ideal(p)

                P.is_prime.set_cache(True)
                P._prime_below = ideal
                P._relative_degree = n
                P._ramification_index = 1
                P._beta = [1] + [0]*(n-1)
            else:
                Q = q.basis_matrix().apply_map(lambda e: e.lift())
                P = self.ideal(p, *Q*vector(self.basis()))

                # Now we compute an element beta in O but not in pO such that
                # beta*P in pO.

                # Since beta is in k[x]-module O, we keep beta as a vector
                # in k[x] with respect to the basis of O. As long as at least
                # one element in this vector is not divisible by p, beta will
                # not be in pO. To ensure that beta*P is in pO, multiplying
                # beta by each of P's generators must produce a vector whose
                # elements are multiples of p. We can ensure that all this
                # occurs by constructing a matrix in k, and finding a nonzero
                # vector in the kernel of the matrix.

                m = []
                for g in q.basis_matrix():
                    m.extend(matrix([g * mr for mr in matrices_reduced]).columns())
                beta = [c.lift() for c in matrix(m).right_kernel().basis()[0]]

                r = q
                index = 1
                while True:
                    rq = r*q
                    if rq == r:
                        break
                    r = rq
                    index = index + 1

                P.is_prime.set_cache(True)
                P._prime_below = ideal
                P._relative_degree = n - q.basis_matrix().nrows()
                P._ramification_index = index
                P._beta = beta

            factors.append((P, P._relative_degree, P._ramification_index))

        return factors


class FunctionFieldMaximalOrderInfinite_polymod(FunctionFieldMaximalOrderInfinite):
    """
    Maximal infinite orders of function fields.

    INPUT:

    - ``field`` -- function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)                   # needs sage.rings.finite_rings
        sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)                                # needs sage.rings.finite_rings
        sage: F.maximal_order_infinite()                                                # needs sage.rings.finite_rings
        Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2

        sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                                 # needs sage.rings.finite_rings
        sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                    # needs sage.rings.finite_rings
        sage: L.maximal_order_infinite()                                                # needs sage.rings.finite_rings
        Maximal infinite order of Function field in y defined by y^2 + y + (x^2 + 1)/x
    """
    def __init__(self, field, category=None):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order_infinite()
            sage: TestSuite(O).run()
        """
        FunctionFieldMaximalOrderInfinite.__init__(self, field, ideal_class=FunctionFieldIdealInfinite_polymod)

        M, from_M, to_M = field._inversion_isomorphism()
        basis = [from_M(g) for g in M.maximal_order().basis()]

        V, from_V, to_V = field.vector_space()
        R = field.base_field().maximal_order_infinite()

        self._basis = tuple(basis)
        self._module = V.span_of_basis([to_V(v) for v in basis])
        self._module_base_ring = R
        self._to_module = to_V

    def _element_constructor_(self, f):
        """
        Make ``f`` an element of this order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)
            sage: 1 in Oinf
            True
            sage: 1/x*y in Oinf
            True
            sage: x*y in Oinf
            False
            sage: 1/x in Oinf
            True
        """
        F = self.function_field()

        try:
            f = F(f)
        except TypeError:
            raise TypeError("unable to convert to an element of {}".format(F))

        O = F.base_field().maximal_order_infinite()
        coordinates = self.coordinate_vector(f)
        if not all(c in O for c in coordinates):
            raise TypeError("%r is not an element of %r" % (f,self))

        return f

    def basis(self):
        """
        Return a basis of this order as a module over the maximal order
        of the base function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x^2*y, (1/(x^4 + x^3 + x^2))*y^2)

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)
        """
        return self._basis

    def gen(self, n=0):
        """
        Return the ``n``-th generator of the order.

        The basis elements of the order are generators.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.gen()
            1
            sage: Oinf.gen(1)
            1/x^2*y
            sage: Oinf.gen(2)
            (1/(x^4 + x^3 + x^2))*y^2
            sage: Oinf.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: there are only 3 generators
        """
        if not ( n >= 0 and n < self.ngens() ):
            raise IndexError("there are only {} generators".format(self.ngens()))

        return self._basis[n]

    def ngens(self):
        """
        Return the number of generators of the order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.ngens()
            3
        """
        return len(self._basis)

    def ideal(self, *gens):
        """
        Return the ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- tuple of elements of the function field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x, y); I
            Ideal (y) of Maximal infinite order of Function field
            in y defined by y^3 + x^6 + x^4 + x^2

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(x, y); I
            Ideal (x) of Maximal infinite order of Function field
            in y defined by y^2 + y + (x^2 + 1)/x
        """
        if len(gens) == 1:
            gens = gens[0]
            if not type(gens) in (list,tuple):
                gens = (gens,)
        mgens = [g * b for g in gens for b in self._basis]
        return self.ideal_with_gens_over_base(mgens)

    def ideal_with_gens_over_base(self, gens):
        """
        Return the ideal generated by ``gens`` as a module.

        INPUT:

        - ``gens`` -- tuple of elements of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: Oinf.ideal_with_gens_over_base((x^2, y, (1/(x^2 + x + 1))*y^2))
            Ideal (y) of Maximal infinite order of Function field in y
            defined by y^3 + x^6 + x^4 + x^2
        """
        F = self.function_field()
        iF, from_iF, to_iF = F._inversion_isomorphism()
        iO = iF.maximal_order()

        ideal = iO.ideal_with_gens_over_base([to_iF(g) for g in gens])

        if not ideal.is_zero():
            # Now the ideal does not correspond exactly to the ideal in the
            # maximal infinite order through the inversion isomorphism. The
            # reason is that the ideal also has factors not lying over x.
            # The following procedure removes the spurious factors. The idea
            # is that for an integral ideal I, J_n = I + (xO)^n stabilizes
            # if n is large enough, and then J_n is the I with the spurious
            # factors removed. For a fractional ideal, we also need to find
            # the largest factor x^m that divides the denominator.
            from sage.matrix.special import block_matrix
            from .hermite_form_polynomial import reversed_hermite_form

            d = ideal.denominator()
            h = ideal.hnf()
            x = d.parent().gen()

            # find the largest factor x^m that divides the denominator
            i = 0
            while d[i].is_zero():
                i += 1
            d = x ** i

            # find the largest n such that I + (xO)^n stabilizes
            h1 = h
            MS = h1.matrix_space()
            k = MS.identity_matrix()
            while True:
                k = x * k

                h2 = block_matrix([[h],[k]])
                reversed_hermite_form(h2)
                i = 0
                while i < h2.nrows() and h2.row(i).is_zero():
                    i += 1
                h2 = h2[i:] # remove zero rows

                if h2 == h1:
                    break
                h1 = h2

            # reconstruct ideal
            ideal = iO._ideal_from_vectors_and_denominator(list(h1), d)

        return self.ideal_monoid().element_class(self, ideal)

    def _to_iF(self, I):
        """
        Return the ideal in the inverted function field from ``I``.

        INPUT:

        - ``I`` -- ideal of the function field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(y)
            sage: Oinf._to_iF(I)
            Ideal (1, 1/x*s) of Maximal order of Function field in s
            defined by s^2 + x*s + x^3 + x
        """
        F = self.function_field()
        iF,from_iF,to_iF = F._inversion_isomorphism()
        iO = iF.maximal_order()
        iI = iO.ideal_with_gens_over_base([to_iF(b) for b in I.gens_over_base()])
        return iI

    def decomposition(self):
        r"""
        Return prime ideal decomposition of `pO_\infty` where `p` is the unique
        prime ideal of the maximal infinite order of the rational function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: Oinf.decomposition()
            [(Ideal ((1/(x^4 + x^3 + x^2))*y^2 + 1) of Maximal infinite order
             of Function field in y defined by y^3 + x^6 + x^4 + x^2, 1, 1),
             (Ideal ((1/(x^4 + x^3 + x^2))*y^2 + 1/x^2*y + 1) of Maximal infinite order
             of Function field in y defined by y^3 + x^6 + x^4 + x^2, 2, 1)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.decomposition()
            [(Ideal (1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x, 1, 2)]

        ::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: Oinf.decomposition()
            [(Ideal (1/x^2*y - 1) of Maximal infinite order
             of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2, 1, 1),
             (Ideal ((1/(x^4 + x^3 + x^2))*y^2 + 1/x^2*y + 1) of Maximal infinite order
             of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2, 2, 1)]

        ::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.decomposition()
            [(Ideal (1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x, 1, 2)]
        """
        F = self.function_field()
        iF,from_iF,to_iF = F._inversion_isomorphism()

        x = iF.base_field().gen()
        iO = iF.maximal_order()
        io = iF.base_field().maximal_order()
        ip = io.ideal(x)

        dec = []
        for iprime, deg, exp in iO.decomposition(ip):
            prime = self.ideal_monoid().element_class(self, iprime)
            dec.append((prime, deg, exp))
        return dec

    def different(self):
        """
        Return the different ideal of the maximal infinite order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.different()
            Ideal (1/x) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        T = self._codifferent_matrix()
        codiff_gens = []
        for c in T.inverse().columns():
            codiff_gens.append(sum([ci*bi for ci,bi in zip(c,self.basis())]))
        codiff = self.ideal_with_gens_over_base(codiff_gens)
        return ~codiff

    @cached_method
    def _codifferent_matrix(self):
        """
        Return the codifferent matrix of the maximal infinite order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf._codifferent_matrix()
            [    0   1/x]
            [  1/x 1/x^2]
        """
        from sage.matrix.constructor import matrix

        rows = []
        for u in self.basis():
            row = []
            for v in self.basis():
                row.append((u*v).trace())
            rows.append(row)
        T = matrix(rows)
        return T

    def coordinate_vector(self, e):
        """
        Return the coordinates of ``e`` with respect to the basis of the order.

        INPUT:

        - ``e`` -- element of the function field

        The returned coordinates are in the base maximal infinite order if and only
        if the element is in the order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: f = 1/y^2
            sage: f in Oinf
            True
            sage: Oinf.coordinate_vector(f)
            ((x^3 + x^2 + x)/(x^4 + 1), x^3/(x^4 + 1))
        """
        return self._module.coordinate_vector(self._to_module(e))


class FunctionFieldMaximalOrder_global(FunctionFieldMaximalOrder_polymod):
    """
    Maximal orders of global function fields.

    INPUT:

    - ``field`` -- function field to which this maximal order belongs

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]                                 # needs sage.rings.finite_rings
        sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)                                  # needs sage.rings.finite_rings
        sage: L.maximal_order()                                                         # needs sage.rings.finite_rings
        Maximal order of Function field in y defined by y^4 + x*y + 4*x + 1
    """

    def __init__(self, field):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: TestSuite(O).run()
        """
        FunctionFieldMaximalOrder_polymod.__init__(self, field, ideal_class=FunctionFieldIdeal_global)

    @cached_method
    def p_radical(self, prime):
        """
        Return the ``prime``-radical of the maximal order.

        INPUT:

        - ``prime`` -- prime ideal of the maximal order of the base
          rational function field

        The algorithm is outlined in Section 6.1.3 of [Coh1993]_.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2 * (x^2 + x + 1)^2)
            sage: o = K.maximal_order()
            sage: O = F.maximal_order()
            sage: p = o.ideal(x + 1)
            sage: O.p_radical(p)
            Ideal (x + 1) of Maximal order of Function field in y
            defined by y^3 + x^6 + x^4 + x^2
        """
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector

        g = prime.gens()[0]

        if not (g.denominator() == 1 and g.numerator().is_irreducible()):
            raise ValueError('not a prime ideal')

        F = self.function_field()
        n = F.degree()
        o = prime.ring()
        p = g.numerator()

        # Fp is isomorphic to the residue field o/p
        Fp, fr_Fp, to_Fp = o._residue_field_global(p)

        # exp = q^j should be at least extension degree where q is
        # the order of the residue field o/p
        q = F.constant_base_field().order()**p.degree()
        exp = q
        while exp <= F.degree():
            exp = exp**q

        # radical equals to the kernel of the map x |-> x^exp
        mat = []
        for g in self.basis():
            v = [to_Fp(c) for c in self._coordinate_vector(g**exp)]
            mat.append(v)
        mat = matrix(Fp, mat)
        ker = mat.kernel()

        # construct module generators of the p-radical
        vecs = []
        for i in range(n):
            v = vector([p if j == i else 0 for j in range(n)])
            vecs.append(v)
        for b in ker.basis():
            v = vector([fr_Fp(c) for c in b])
            vecs.append(v)

        return self._ideal_from_vectors(vecs)

    @cached_method
    def decomposition(self, ideal):
        """
        Return the decomposition of the prime ideal.

        INPUT:

        - ``ideal`` -- prime ideal of the base maximal order

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: o = K.maximal_order()
            sage: O = F.maximal_order()
            sage: p = o.ideal(x + 1)
            sage: O.decomposition(p)
            [(Ideal (x + 1, y + 1) of Maximal order
             of Function field in y defined by y^3 + x^6 + x^4 + x^2, 1, 1),
             (Ideal (x + 1, (1/(x^3 + x^2 + x))*y^2 + y + 1) of Maximal order
             of Function field in y defined by y^3 + x^6 + x^4 + x^2, 2, 1)]
        """
        from sage.matrix.constructor import matrix

        F = self.function_field()
        n = F.degree()

        p = ideal.gen().numerator()
        o = ideal.ring()

        # Fp is isomorphic to the residue field o/p
        Fp, fr, to = o._residue_field_global(p)
        P,X = Fp['X'].objgen()

        V = Fp**n # Ob = O/pO

        mtable = []
        for i in range(n):
            row = []
            for j in range(n):
                row.append( V([to(e) for e in self._mtable[i][j]]) )
            mtable.append(row)

        if p not in self._kummer_places:
            #####################################
            # Decomposition by Kummer's theorem #
            #####################################
            # gen is self._kummer_gen
            gen_vec_pow = self._kummer_gen_vec_pow
            mul_vecs = self._mul_vecs

            f = self._kummer_polynomial
            fp = P([to(c.numerator()) for c in f.list()])
            decomposition = []
            for q, exp in fp.factor():
                # construct O.ideal([p,q(gen)])
                gen_vecs = list(matrix.diagonal(n * [p]))
                c = q.list()

                # q(gen) in vector form
                qgen = sum(fr(c[i]) * gen_vec_pow[i] for i in range(len(c)))

                I = matrix.identity(o._ring, n)
                for i in range(n):
                    gen_vecs.append(mul_vecs(qgen,I[i]))
                prime = self._ideal_from_vectors_and_denominator(gen_vecs)

                # Compute an element beta in O but not in pO. How to find beta
                # is explained in Section 4.8.3 of [Coh1993]. We keep beta
                # as a vector over k[x] with respect to the basis of O.

                # p and qgen generates the prime; modulo pO, qgenb generates the prime
                qgenb = [to(qgen[i]) for i in range(n)]
                m = []
                for i in range(n):
                    m.append(sum(qgenb[j] * mtable[i][j] for j in range(n)))
                beta = [fr(coeff) for coeff in matrix(m).left_kernel().basis()[0]]

                prime.is_prime.set_cache(True)
                prime._prime_below = ideal
                prime._relative_degree = q.degree()
                prime._ramification_index = exp
                prime._beta = beta

                prime._kummer_form = (p, qgen)

                decomposition.append((prime, q.degree(), exp))
        else:
            #############################
            # Buchman-Lenstra algorithm #
            #############################
            from sage.matrix.special import block_matrix
            from sage.modules.free_module_element import vector

            pO = self.ideal(p)
            Ip = self.p_radical(ideal)
            Ob = matrix.identity(Fp, n)

            def bar(I):  # transfer to O/pO
                m = []
                for v in I._hnf:
                    m.append([to(e) for e in v])
                h = matrix(m).echelon_form()
                return cut_last_zero_rows(h)

            def liftb(Ib):
                m = [vector([fr(e) for e in v]) for v in Ib]
                m.extend(pO._hnf)
                return self._ideal_from_vectors_and_denominator(m, 1)

            def cut_last_zero_rows(h):
                i = h.nrows()
                while i > 0 and h.row(i-1).is_zero():
                    i -= 1
                return h[:i]

            def mul_vec(v1, v2):
                s = 0
                for i in range(n):
                    for j in range(n):
                        s += v1[i] * v2[j] * mtable[i][j]
                return s

            def pow(v, r):  # r > 0
                m = v
                while r > 1:
                    m = mul_vec(m, v)
                    r -= 1
                return m

            # Algorithm 6.2.7 of [Coh1993]
            def div(Ib, Jb):
                # compute a basis of Jb/Ib
                sJb = Jb.row_space()
                sIb = Ib.row_space()
                sJbsIb,proj_sJbsIb,lift_sJbsIb = sJb.quotient_abstract(sIb)
                supplement_basis = [lift_sJbsIb(v) for v in sJbsIb.basis()]

                m = []
                for b in V.gens(): # basis of Ob = O/pO
                    b_row = [] # row vector representation of the map a -> a*b
                    for a in supplement_basis:
                        b_row += lift_sJbsIb(proj_sJbsIb( mul_vec(a,b) ))
                    m.append(b_row)
                return matrix(Fp,n,m).left_kernel().basis_matrix()

            # Algorithm 6.2.5 of [Coh1993]
            def mul(Ib, Jb):
                m = []
                for v1 in Ib:
                    for v2 in Jb:
                        m.append(mul_vec(v1,v2))
                h = matrix(m).echelon_form()
                return cut_last_zero_rows(h)

            def add(Ib,Jb):
                m = block_matrix([[Ib], [Jb]])
                h = m.echelon_form()
                return cut_last_zero_rows(h)

            # K_1, K_2, ...
            Lb = IpOb = bar(Ip+pO)
            Kb = [Lb]
            while not Lb.is_zero():
                Lb = mul(Lb,IpOb)
                Kb.append(Lb)

            # J_1, J_2, ...
            Jb = [Kb[0]] + [div(Kb[j],Kb[j-1]) for j in range(1,len(Kb))]

            # H_1, H_2, ...
            Hb = [div(Jb[j],Jb[j+1]) for j in range(len(Jb)-1)] + [Jb[-1]]

            q = Fp.order()

            def split(h):
                # VsW represents O/H as a vector space
                W = h.row_space() # H/pO
                VsW,to_VsW,lift_to_V = V.quotient_abstract(W)

                # compute the space K of elements in O/H that satisfy a^q-a=0
                l = [lift_to_V(b) for b in VsW.basis()]

                images = [to_VsW(pow(x, q) - x) for x in l]
                K = VsW.hom(images, VsW).kernel()

                if K.dimension() == 0:
                    return []
                if K.dimension() == 1: # h is prime
                    return [(liftb(h),VsW.dimension())] # relative degree

                # choose a such that a^q - a is 0 but a is not in Fp
                for a in K.basis():
                    # IMPORTANT: This criterion is based on the assumption
                    # that O.basis() starts with 1.
                    if a.support() != [0]:
                        break
                else:
                    raise AssertionError("no appropriate value found")

                a = lift_to_V(a)
                # compute the minimal polynomial of a
                m = [to_VsW(Ob[0])] # 1 in VsW
                apow = a
                while True:
                    v = to_VsW(apow)
                    try:
                        sol = matrix(m).solve_left(v)
                    except ValueError:
                        m.append(v)
                        apow = mul_vec(apow, a)
                        continue
                    break

                minpol = X**len(sol) - P(list(sol))

                # The minimal polynomial of a has only linear factors and at least two
                # of them. We set f to the first factor and g to the product of the rest.
                fac = minpol.factor()
                f = fac[0][0]
                g = (fac/f).expand()
                d,u,v = f.xgcd(g)

                assert d == 1, "Not relatively prime {} and {}".format(f,g)

                # finally, idempotent!
                e = lift_to_V(sum([c1*c2 for c1,c2 in zip(u*f,m)]))

                h1 = add(h, matrix([mul_vec(e,Ob[i]) for i in range(n)]))
                h2 = add(h, matrix([mul_vec(Ob[0]-e,Ob[i]) for i in range(n)]))

                return split(h1) + split(h2)

            decomposition = []
            for i in range(len(Hb)):
                index = i + 1 # Hb starts with H_1
                for prime, degree in split(Hb[i]):
                    # Compute an element beta in O but not in pO. How to find beta
                    # is explained in Section 4.8.3 of [Coh1993]. We keep beta
                    # as a vector over k[x] with respect to the basis of O.
                    m = []
                    for i in range(n):
                        r = []
                        for g in prime._hnf:
                            r += sum(to(g[j]) * mtable[i][j] for j in range(n))
                        m.append(r)
                    beta = [fr(e) for e in matrix(m).left_kernel().basis()[0]]

                    prime.is_prime.set_cache(True)
                    prime._prime_below = ideal
                    prime._relative_degree = degree
                    prime._ramification_index = index
                    prime._beta = beta

                    decomposition.append((prime, degree, index))

        return decomposition
