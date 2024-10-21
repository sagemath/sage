# sage_setup: distribution = sagemath-singular
# sage.doctest: needs sage.rings.function_field
"""
Places of function fields: extension
"""

# ****************************************************************************
#       Copyright (C) 2016-2022 Kwankyu Lee <ekwankyu@gmail.com>
#                     2019      Brent Baccala
#                     2021      Jonathan Kliem
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import sage
from sage.arith.functions import lcm
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.rings.number_field.number_field_base import NumberField

from .place import FunctionFieldPlace


class FunctionFieldPlace_polymod(FunctionFieldPlace):
    """
    Places of extensions of function fields.
    """
    def place_below(self):
        """
        Return the place lying below the place.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: OK = K.maximal_order()
            sage: OL = L.maximal_order()
            sage: p = OK.ideal(x^2 + x + 1)
            sage: dec = OL.decomposition(p)
            sage: q = dec[0][0].place()
            sage: q.place_below()
            Place (x^2 + x + 1)
        """
        return self.prime_ideal().prime_below().place()

    def relative_degree(self):
        """
        Return the relative degree of the place.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: OK = K.maximal_order()
            sage: OL = L.maximal_order()
            sage: p = OK.ideal(x^2 + x + 1)
            sage: dec = OL.decomposition(p)
            sage: q = dec[0][0].place()
            sage: q.relative_degree()
            1
        """
        return self._prime._relative_degree

    def degree(self):
        """
        Return the degree of the place.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: OK = K.maximal_order()
            sage: OL = L.maximal_order()
            sage: p = OK.ideal(x^2 + x + 1)
            sage: dec = OL.decomposition(p)
            sage: q = dec[0][0].place()
            sage: q.degree()
            2
        """
        return self.relative_degree() * self.place_below().degree()

    def is_infinite_place(self):
        """
        Return ``True`` if the place is above the unique infinite place
        of the underlying rational function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: pls = L.places()
            sage: [p.is_infinite_place() for p in pls]
            [True, True, False]
            sage: [p.place_below() for p in pls]
            [Place (1/x), Place (1/x), Place (x)]
        """
        return self.place_below().is_infinite_place()

    def local_uniformizer(self):
        """
        Return an element of the function field that has a simple zero
        at the place.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: pls = L.places()
            sage: [p.local_uniformizer().valuation(p) for p in pls]
            [1, 1, 1, 1, 1]
        """
        gens = self._prime.gens()
        for g in gens:
            if g.valuation(self) == 1:
                return g
        assert False, "Internal error"

    def gaps(self):
        """
        Return the gap sequence for the place.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x,y).place()
            sage: p.gaps() # a Weierstrass place
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)                                # needs sage.rings.finite_rings
            sage: [p.gaps() for p in L.places()]                                        # needs sage.rings.finite_rings
            [[1, 2, 4], [1, 2, 4], [1, 2, 4]]
        """
        if self.degree() == 1:
            return self._gaps_rational() # faster for rational places
        else:
            return self._gaps_wronskian()

    def _gaps_rational(self):
        """
        Return the gap sequence for the rational place.

        This method computes the gap numbers using the definition of gap
        numbers. The dimension of the multiple of the prime divisor
        supported at the place is computed by Hess' algorithm.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x, y).place()
            sage: p.gaps()  # indirect doctest
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                  # needs sage.rings.finite_rings
            sage: [p.gaps() for p in L.places()]  # indirect doctest                    # needs sage.rings.finite_rings
            [[1, 2, 4], [1, 2, 4], [1, 2, 4]]
        """
        from sage.matrix.constructor import matrix

        F = self.function_field()
        n = F.degree()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        R = O._module_base_ring._ring
        one = R.one()

        # Hess' Riemann-Roch basis algorithm stripped down for gaps computation
        def dim_RR(M):
            den = lcm([e.denominator() for e in M.list()])
            mat = matrix(R, M.nrows(), [(den*e).numerator() for e in M.list()])

            # initialise pivot_row and conflicts list
            pivot_row = [[] for i in range(n)]
            conflicts = []
            for i in range(n):
                bestp = -1
                best = -1
                for c in range(n):
                    d = mat[i,c].degree()
                    if d >= best:
                        bestp = c
                        best = d

                if best >= 0:
                    pivot_row[bestp].append((i,best))
                    if len(pivot_row[bestp]) > 1:
                        conflicts.append(bestp)

            # while there is a conflict, do a simple transformation
            while conflicts:
                c = conflicts.pop()
                row = pivot_row[c]
                i,ideg = row.pop()
                j,jdeg = row.pop()

                if jdeg > ideg:
                    i,j = j,i
                    ideg,jdeg = jdeg,ideg

                coeff = - mat[i,c].lc() / mat[j,c].lc()
                s = coeff * one.shift(ideg - jdeg)

                mat.add_multiple_of_row(i, j, s)

                row.append((j,jdeg))

                bestp = -1
                best = -1
                for c in range(n):
                    d = mat[i,c].degree()
                    if d >= best:
                        bestp = c
                        best = d

                if best >= 0:
                    pivot_row[bestp].append((i,best))
                    if len(pivot_row[bestp]) > 1:
                        conflicts.append(bestp)

            dim = 0
            for j in range(n):
                i,ideg = pivot_row[j][0]
                k = den.degree() - ideg + 1
                if k > 0:
                    dim += k
            return dim

        V,fr,to = F.vector_space()

        prime_inv = ~ self.prime_ideal()
        I = O.ideal(1)
        J = Oinf.ideal(1)

        B = matrix([to(b) for b in J.gens_over_base()])
        C = matrix([to(v) for v in I.gens_over_base()])

        prev = dim_RR(C * B.inverse())
        gaps = []
        g = F.genus()
        i = 1
        if self.is_infinite_place():
            while g:
                J = J * prime_inv
                B = matrix([to(b) for b in J.gens_over_base()])
                dim = dim_RR(C * B.inverse())
                if dim == prev:
                    gaps.append(i)
                    g -= 1
                else:
                    prev = dim
                i += 1
        else: # self is a finite place
            Binv = B.inverse()
            while g:
                I = I * prime_inv
                C = matrix([to(v) for v in I.gens_over_base()])
                dim = dim_RR(C * Binv)
                if dim == prev:
                    gaps.append(i)
                    g -= 1
                else:
                    prev = dim
                i += 1

        return gaps

    def _gaps_wronskian(self):
        """
        Return the gap sequence for the place.

        This method implements the local version of Hess' Algorithm 30 of [Hes2002b]_
        based on the Wronskian determinant.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x, y).place()
            sage: p._gaps_wronskian()  # a Weierstrass place
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)                                # needs sage.rings.finite_rings
            sage: [p._gaps_wronskian() for p in L.places()]                             # needs sage.rings.finite_rings
            [[1, 2, 4], [1, 2, 4], [1, 2, 4]]
        """
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector

        F = self.function_field()
        R,fr_R,to_R = self._residue_field()
        der = F.higher_derivation()

        sep = self.local_uniformizer()

        # a differential divisor satisfying
        # v_p(W) = 0 for the place p
        W = sep.differential().divisor()

        # Step 3:
        basis = W._basis()
        d = len(basis)
        M = matrix([to_R(b) for b in basis])
        if M.rank() == 0:
            return []

        # Steps 4, 5, 6, 7:
        e = 1
        gaps = [1]
        while M.nrows() < d:
            row = vector([to_R(der._derive(basis[i], e, sep)) for i in range(d)])
            if row not in M.row_space():
                M = matrix(M.rows() + [row])
                M.echelonize()
                gaps.append(e + 1)
            e += 1

        return gaps

    def residue_field(self, name=None):
        """
        Return the residue field of the place.

        INPUT:

        - ``name`` -- string; name of the generator of the residue field

        OUTPUT:

        - a field isomorphic to the residue field

        - a ring homomorphism from the valuation ring to the field

        - a ring homomorphism from the field to the valuation ring

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: k, fr_k, to_k = p.residue_field()
            sage: k
            Finite Field of size 2
            sage: fr_k
            Ring morphism:
              From: Finite Field of size 2
              To:   Valuation ring at Place (x, x*y)
            sage: to_k
            Ring morphism:
              From: Valuation ring at Place (x, x*y)
              To:   Finite Field of size 2
            sage: to_k(y)
            Traceback (most recent call last):
            ...
            TypeError: y fails to convert into the map's domain
            Valuation ring at Place (x, x*y)...
            sage: to_k(1/y)
            0
            sage: to_k(y/(1+y))
            1
        """
        return self.valuation_ring().residue_field(name=name)

    @cached_method
    def _residue_field(self, name=None):
        """
        Return the residue field of the place along with the functions
        mapping from and to it.

        INPUT:

        - ``name`` -- string (default: ``None``); name of the generator
          of the residue field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: k,fr_k,to_k = p._residue_field()
            sage: k
            Finite Field of size 2
            sage: [fr_k(e) for e in k]
            [0, 1]

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(9)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + Y - x^4)
            sage: p = L.places()[-1]
            sage: p.residue_field()
            (Finite Field in z2 of size 3^2, Ring morphism:
               From: Finite Field in z2 of size 3^2
               To:   Valuation ring at Place (x + 1, y + 2*z2), Ring morphism:
               From: Valuation ring at Place (x + 1, y + 2*z2)
               To:   Finite Field in z2 of size 3^2)

        ::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + Y - x^4)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x)
            sage: [p.residue_field() for p in L.places_above(I.place())]
            [(Rational Field, Ring morphism:
                From: Rational Field
                To:   Valuation ring at Place (x, y, y^2), Ring morphism:
                From: Valuation ring at Place (x, y, y^2)
                To:   Rational Field),
             (Number Field in s with defining polynomial x^2 - 2*x + 2, Ring morphism:
                From: Number Field in s with defining polynomial x^2 - 2*x + 2
                To:   Valuation ring at Place (x, x*y, y^2 + 1), Ring morphism:
                From: Valuation ring at Place (x, x*y, y^2 + 1)
                To:   Number Field in s with defining polynomial x^2 - 2*x + 2)]
            sage: for p in L.places_above(I.place()):
            ....:    k, fr_k, to_k = p.residue_field()
            ....:    assert all(fr_k(k(e)) == e for e in range(10))
            ....:    assert all(to_k(fr_k(e)) == e for e in [k.random_element() for i in [1..10]])

        ::

            sage: # needs sage.rings.number_field
            sage: K.<x> = FunctionField(QQbar); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + Y - x^4)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x)
            sage: [p.residue_field() for p in L.places_above(I.place())]
            [(Algebraic Field, Ring morphism:
                From: Algebraic Field
                To:   Valuation ring at Place (x, y - I, y^2 + 1), Ring morphism:
                From: Valuation ring at Place (x, y - I, y^2 + 1)
                To:   Algebraic Field), (Algebraic Field, Ring morphism:
                From: Algebraic Field
                To:   Valuation ring at Place (x, y, y^2), Ring morphism:
                From: Valuation ring at Place (x, y, y^2)
                To:   Algebraic Field), (Algebraic Field, Ring morphism:
                From: Algebraic Field
                To:   Valuation ring at Place (x, y + I, y^2 + 1), Ring morphism:
                From: Valuation ring at Place (x, y + I, y^2 + 1)
                To:   Algebraic Field)]
        """
        F = self.function_field()
        prime = self.prime_ideal()  # Let P be this prime ideal

        if self.is_infinite_place():
            _F, from_F, to_F = F._inversion_isomorphism()
            _prime = prime._ideal
            _place = _prime.place()

            K, _from_K, _to_K = _place._residue_field(name=name)

            from_K = lambda e: from_F(_from_K(e))
            to_K = lambda f: _to_K(to_F(f))
            return K, from_K, to_K

        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector

        O = F.maximal_order()
        Obasis = O.basis()

        M = prime.hnf()
        R = M.base_ring() # univariate polynomial ring
        n = M.nrows() # extension degree of the function field

        # Step 1: construct a vector space representing the residue field
        #
        # Given an (reversed) HNF basis M for a prime ideal P of O, every
        # element of O mod P can be represented by a vector of polynomials of
        # degrees less than those of the (anti)diagonal elements of M. In turn,
        # the vector of polynomials can be represented by the vector of the
        # coefficients of the polynomials. V is the space of these vectors.

        k = F.constant_base_field()
        degs = [M[i,i].degree() for i in range(n)]
        deg = sum(degs) # degree of the place

        # Let V = k**deg

        def to_V(e):
            """
            An example to show the idea: Suppose that::

                    [x 0 0]
                M = [0 1 0] and v = (x^10, x^7 + x^3, x^7 + x^4 + x^3 + 1)
                    [1 0 1]

            Then to_V(e) = [1]
            """
            v = O._coordinate_vector(e)
            vec = []
            for i in reversed(range(n)):
                q,r = v[i].quo_rem(M[i,i])
                v -= q * M[i]
                for j in range(degs[i]):
                    vec.append(r[j])
            return vector(vec)

        def fr_V(vec): # to_O
            vec = vec.list()
            pos = 0
            e = F(0)
            for i in reversed(range(n)):
                if degs[i] == 0:
                    continue
                else:
                    end = pos + degs[i]
                    e += R(vec[pos:end]) * Obasis[i]
                    pos = end
            return e

        # Step 2: find a primitive element of the residue field

        def candidates():
            # Trial 1: this suffices for places obtained from Kummers' theorem
            # and for places of function fields over number fields or QQbar

            # Note that a = O._kummer_gen is a simple generator of O/prime over
            # o/p. If b is a simple generator of o/p over the constant base field
            # k, then the set a + k * b contains a simple generator of O/prime
            # over k (as there are finite number of intermediate fields).
            a = O._kummer_gen
            if a is not None:
                K,fr_K,_ = self.place_below().residue_field()
                b = fr_K(K.gen())
                if isinstance(k, (NumberField, sage.rings.abc.AlgebraicField)):
                    kk = ZZ
                else:
                    kk = k
                for c in kk:
                    if c != 0:
                        yield a + c * b

            # Trial 2: basis elements of the maximal order
            for gen in reversed(Obasis):
                yield gen

            import itertools

            # Trial 3: exhaustive search in O using only polynomials
            # with coefficients 0 or 1
            for d in range(deg):
                G = itertools.product(itertools.product([0,1],repeat=d+1), repeat=n)
                for g in G:
                    gen = sum([R(c1)*c2 for c1,c2 in zip(g, Obasis)])
                    yield gen

            # Trial 4: exhaustive search in O using all polynomials
            for d in range(deg):
                G = itertools.product(R.polynomials(max_degree=d), repeat=n)
                for g in G:
                    # discard duplicate cases
                    if max(c.degree() for c in g) != d:
                        continue
                    for j in range(n):
                        if g[j] != 0:
                            break
                    if g[j].leading_coefficient() != 1:
                        continue

                    gen = sum([c1*c2 for c1,c2 in zip(g, Obasis)])
                    yield gen

        # Search for a primitive element. It is such an element g of O
        # whose powers span the vector space V.
        for gen in candidates():
            g = F.one()
            m = []
            for i in range(deg):
                m.append(to_V(g))
                g *= gen
            mat = matrix(m)
            if mat.rank() == deg:
                break

        # Step 3: compute the minimal polynomial of g
        min_poly = R((-mat.solve_left(to_V(g))).list() + [1])

        # Step 4: construct the residue field K as an extension of the base
        # constant field using the minimal polynomial and compute vector space
        # representation W of K along with maps between them
        if deg > 1:
            if isinstance(k, NumberField):
                if name is None:
                    name = 's'
                K = k.extension(min_poly, names=name)

                def from_W(e):
                    return K(list(e))

                def to_W(e):
                    return vector(K(e))
            else:
                K = k.extension(deg, name=name)

                # primitive element in K corresponding to g in O mod P
                prim = min_poly.roots(K)[0][0]

                W, from_W, to_W = K.vector_space(k, basis=[prim**i for i in range(deg)], map=True)
        else: # deg == 1
            K = k

            def from_W(e):
                return K(e[0])

            def to_W(e):
                return vector([e])

        # Step 5: compute the matrix of change of basis, from V to W via K
        C = mat.inverse()

        # Step 6: construct the maps between the residue field of the valuation
        # ring at P and K, via O and V and W

        def from_K(e):
            return fr_V(to_W(e) * mat)

        # As explained in Section 4.8.3 of [Coh1993]_, alpha has a simple pole
        # at this place and no other poles at finite places.
        p = prime.prime_below().gen().numerator()
        beta = prime._beta
        alpha = ~p * sum(c1*c2 for c1,c2 in zip(beta, Obasis))
        alpha_powered_by_ramification_index = alpha ** prime._ramification_index

        def to_K(f):
            if f not in O:
                den = O.coordinate_vector(f).denominator()
                num = den * f

                # s powered by the valuation of den at the prime
                alpha_power = alpha_powered_by_ramification_index ** den.valuation(p)
                rn = num * alpha_power # in O
                rd = den * alpha_power # in O but not in prime

                # Note that rn is not in O if and only if f is
                # not in the valuation ring. Hence f is in the
                # valuation ring if and only if this procedure
                # does not fall into an infinite loop.
                return to_K(rn) / to_K(rd)

            return from_W(to_V(f) * C)

        return K, from_K, to_K

    def valuation_ring(self):
        """
        Return the valuation ring at the place.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: p.valuation_ring()
            Valuation ring at Place (x, x*y)
        """
        from .valuation_ring import FunctionFieldValuationRing

        return FunctionFieldValuationRing(self.function_field(), self)
