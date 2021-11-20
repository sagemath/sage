"""
Places of function fields

The places of a function field correspond, one-to-one, to valuation rings of
the function field, each of which defines a discrete valuation for the elements
of the function field. "Finite" places are in one-to-one correspondence with
the prime ideals of the finite maximal order while places "at infinity" are in
one-to-one correspondence with the prime ideals of the infinite maximal order.

EXAMPLES:

All rational places of a function field can be computed::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
    sage: L.places()
    [Place (1/x, 1/x^3*y^2 + 1/x),
     Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1),
     Place (x, y)]

The residue field associated with a place is given as an extension of the
constant field::

    sage: F.<x> = FunctionField(GF(2))
    sage: O = F.maximal_order()
    sage: p = O.ideal(x^2 + x + 1).place()
    sage: k, fr_k, to_k = p.residue_field()
    sage: k
    Finite Field in z2 of size 2^2

The homomorphisms are between the valuation ring and the residue field::

    sage: fr_k
    Ring morphism:
      From: Finite Field in z2 of size 2^2
      To:   Valuation ring at Place (x^2 + x + 1)
    sage: to_k
    Ring morphism:
      From: Valuation ring at Place (x^2 + x + 1)
      To:   Finite Field in z2 of size 2^2

AUTHORS:

- Kwankyu Lee (2017-04-30): initial version

- Brent Baccala (2019-12-20): function fields of characteristic zero

"""
#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method

from sage.arith.all import lcm

from sage.rings.integer_ring import ZZ
from sage.rings.qqbar import QQbar
from sage.rings.number_field.number_field_base import NumberField

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp

from sage.categories.sets_cat import Sets

from sage.modules.free_module_element import vector

from sage.matrix.constructor import matrix


class FunctionFieldPlace(Element):
    """
    Places of function fields.

    INPUT:

    - ``parent`` -- place set of a function field

    - ``prime`` -- prime ideal associated with the place

    EXAMPLES::

        sage: K.<x>=FunctionField(GF(2)); _.<Y> = K[]
        sage: L.<y>=K.extension(Y^3 + x + x^3*Y)
        sage: L.places_finite()[0]
        Place (x, y)
    """
    def __init__(self, parent, prime):
        """
        Initialize the place.

        TESTS::

            sage: K.<x>=FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y>=K.extension(Y^3 + x + x^3*Y)
            sage: p = L.places_finite()[0]
            sage: TestSuite(p).run()
        """
        Element.__init__(self, parent)

        self._prime = prime

    def __hash__(self):
        """
        Return the hash of the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y>=K.extension(Y^3 + x + x^3*Y)
            sage: p = L.places_finite()[0]
            sage: {p: 1}
            {Place (x, y): 1}
        """
        return hash((self.function_field(), self._prime))

    def _repr_(self):
        """
        Return the string representation of the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p = L.places_finite()[0]
            sage: p
            Place (x, y)
        """
        try:
            gens = self._prime.gens_two()
        except AttributeError:
            gens = self._prime.gens()
        gens_str = ', '.join(repr(g) for g in gens)
        return "Place ({})".format(gens_str)

    def _latex_(self):
        r"""
        Return the LaTeX representation of the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: p = L.places_finite()[0]
            sage: latex(p)
            \left(y\right)
        """
        return self._prime._latex_()

    def _richcmp_(self, other, op):
        """
        Compare the place with ``other`` place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: p1, p2, p3 = L.places()[:3]
            sage: p1 < p2
            True
            sage: p2 < p1
            False
            sage: p1 == p3
            False
        """
        from sage.rings.function_field.order import FunctionFieldOrderInfinite

        # effect that places at infinity are ordered first
        s = not isinstance(self._prime.ring(), FunctionFieldOrderInfinite)
        o = not isinstance(other._prime.ring(), FunctionFieldOrderInfinite)
        return richcmp((s, self._prime), (o, other._prime), op)

    def _acted_upon_(self, other, self_on_left):
        """
        Define integer multiplication upon the prime divisor
        of the place on the left.

        The output is a divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x+1,y)
            sage: P = I.place()
            sage: -3*P + 5*P
            2*Place (x + 1, y)
        """
        if self_on_left:
            raise TypeError("only left multiplication by integers is allowed")
        return other * self.divisor()

    def _neg_(self):
        """
        Return the negative of the prime divisor of this place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p1, p2, p3 = L.places()[:3]
            sage: -p1 + p2
            - Place (1/x, 1/x^3*y^2 + 1/x)
             + Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
        """
        from .divisor import divisor
        return divisor(self.function_field(), {self: -1})

    def _add_(self, other):
        """
        Return the divisor that is the sum of the place and ``other``.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p1, p2, p3 = L.places()[:3]
            sage: p1 + p2 + p3
            Place (1/x, 1/x^3*y^2 + 1/x)
             + Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             + Place (x, y)
        """
        from .divisor import prime_divisor
        return prime_divisor(self.function_field(), self) + other

    def __radd__(self, other):
        """
        Return the prime divisor of the place if ``other`` is zero.

        This is only to support the ``sum`` function, that adds
        the argument to initial (int) zero.

        EXAMPLES::

            sage: k.<a>=GF(2)
            sage: K.<x>=FunctionField(k)
            sage: sum(K.places_finite())
            Place (x) + Place (x + 1)

        Note that this does not work, as wanted::

            sage: 0 + K.place_infinite()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +: ...

        The reason is that the ``0`` is a Sage integer, for which
        the coercion system applies.
        """
        if other == 0:
            from .divisor import prime_divisor
            return prime_divisor(self.function_field(), self)
        raise NotImplementedError

    def function_field(self):
        """
        Return the function field to which the place belongs.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p = L.places()[0]
            sage: p.function_field() == L
            True
        """
        return self.parent()._field

    def prime_ideal(self):
        """
        Return the prime ideal associated with the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p = L.places()[0]
            sage: p.prime_ideal()
            Ideal (1/x^3*y^2 + 1/x) of Maximal infinite order of Function field
            in y defined by y^3 + x^3*y + x
        """
        return self._prime

    def divisor(self, multiplicity=1):
        """
        Return the prime divisor corresponding to the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x+1,y)
            sage: P = I.place()
            sage: P.divisor()
            Place (x + 1, y)
        """
        from .divisor import prime_divisor
        return prime_divisor(self.function_field(), self, multiplicity)


class FunctionFieldPlace_rational(FunctionFieldPlace):
    """
    Places of rational function fields.
    """
    def degree(self):
        """
        Return the degree of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: i = O.ideal(x^2+x+1)
            sage: p = i.place()
            sage: p.degree()
            2
        """
        if self.is_infinite_place():
            return 1
        else:
            return self._prime.gen().numerator().degree()

    def is_infinite_place(self):
        """
        Return ``True`` if the place is at infinite.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: F.places()
            [Place (1/x), Place (x), Place (x + 1)]
            sage: [p.is_infinite_place() for p in F.places()]
            [True, False, False]
        """
        F = self.function_field()
        return self.prime_ideal().ring() == F.maximal_order_infinite()

    def local_uniformizer(self):
        """
        Return a local uniformizer of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: F.places()
            [Place (1/x), Place (x), Place (x + 1)]
            sage: [p.local_uniformizer() for p in F.places()]
            [1/x, x, x + 1]
        """
        return self.prime_ideal().gen()

    def residue_field(self, name=None):
        """
        Return the residue field of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: p = O.ideal(x^2 + x + 1).place()
            sage: k, fr_k, to_k = p.residue_field()
            sage: k
            Finite Field in z2 of size 2^2
            sage: fr_k
            Ring morphism:
              From: Finite Field in z2 of size 2^2
              To:   Valuation ring at Place (x^2 + x + 1)
            sage: to_k
            Ring morphism:
              From: Valuation ring at Place (x^2 + x + 1)
              To:   Finite Field in z2 of size 2^2
        """
        return self.valuation_ring().residue_field(name=name)

    def _residue_field(self, name=None):
        """
        Return the residue field of the place along with the maps from
        and to it.

        INPUT:

        - ``name`` -- string; name of the generator of the residue field

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: i = O.ideal(x^2+x+1)
            sage: p = i.place()
            sage: R, fr, to = p._residue_field()
            sage: R
            Finite Field in z2 of size 2^2
            sage: [fr(e) for e in R.list()]
            [0, x, x + 1, 1]
            sage: to(x*(x+1)) == to(x) * to(x+1)
            True
        """
        F = self.function_field()
        prime = self.prime_ideal()

        if self.is_infinite_place():
            K = F.constant_base_field()

            def from_K(e):
                return F(e)

            def to_K(f):
                n = f.numerator()
                d = f.denominator()

                n_deg = n.degree()
                d_deg =d.degree()

                if n_deg < d_deg:
                    return K(0)
                elif n_deg == d_deg:
                    return n.lc() / d.lc()
                else:
                    raise TypeError("not in the valuation ring")
        else:
            O = F.maximal_order()
            K, from_K, _to_K = O._residue_field(prime, name=name)

            def to_K(f):
                if f in O: # f.denominator() is 1
                    return _to_K(f.numerator())
                else:
                    d = F(f.denominator())
                    n = d * f

                    nv = prime.valuation(O.ideal(n))
                    dv = prime.valuation(O.ideal(d))

                    if nv > dv:
                        return K(0)
                    elif dv > nv:
                        raise TypeError("not in the valuation ring")

                    s = ~prime.gen()
                    rd = d * s**dv # in O but not in prime
                    rn = n * s**nv # in O but not in prime
                    return to_K(rn) / to_K(rd)

        return K, from_K, to_K

    def valuation_ring(self):
        """
        Return the valuation ring at the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: p.valuation_ring()
            Valuation ring at Place (x, x*y)
        """
        from .valuation_ring import FunctionFieldValuationRing

        return FunctionFieldValuationRing(self.function_field(), self)


class FunctionFieldPlace_polymod(FunctionFieldPlace):
    """
    Places of extensions of function fields.
    """
    def place_below(self):
        """
        Return the place lying below the place.

        EXAMPLES::

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

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x,y).place()
            sage: p.gaps() # a Weierstrass place
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: [p.gaps() for p in L.places()]
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

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x,y).place()
            sage: p.gaps()  # indirect doctest
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: [p.gaps() for p in L.places()]  # indirect doctest
            [[1, 2, 4], [1, 2, 4], [1, 2, 4]]
        """
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

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x,y).place()
            sage: p._gaps_wronskian() # a Weierstrass place
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: [p._gaps_wronskian() for p in L.places()]
            [[1, 2, 4], [1, 2, 4], [1, 2, 4]]
        """
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

        - ``name`` -- string (default: `None`); name of the generator
          of the residue field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: k,fr_k,to_k = p._residue_field()
            sage: k
            Finite Field of size 2
            sage: [fr_k(e) for e in k]
            [0, 1]

        ::

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
            _F, from_F, to_F  = F._inversion_isomorphism()
            _prime = prime._ideal
            _place = _prime.place()

            K, _from_K, _to_K = _place._residue_field(name=name)

            from_K = lambda e: from_F(_from_K(e))
            to_K = lambda f: _to_K(to_F(f))
            return K, from_K, to_K

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
            An example to show the idea: Suppose that

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
                if isinstance(k, NumberField) or k is QQbar:
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
                    name='s'
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

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: p.valuation_ring()
            Valuation ring at Place (x, x*y)
        """
        from .valuation_ring import FunctionFieldValuationRing

        return FunctionFieldValuationRing(self.function_field(), self)


class PlaceSet(UniqueRepresentation, Parent):
    """
    Sets of Places of function fields.

    INPUT:

    - ``field`` -- function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
        sage: L.place_set()
        Set of places of Function field in y defined by y^3 + x^3*y + x
    """
    Element = FunctionFieldPlace

    def __init__(self, field):
        """
        Initialize the set of places of the function ``field``.

        TESTS::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: places = L.place_set()
            sage: TestSuite(places).run()
        """
        self.Element = field._place_class
        Parent.__init__(self, category = Sets().Infinite())

        self._field = field

    def _repr_(self):
        """
        Return the string representation of the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: L.place_set()
            Set of places of Function field in y defined by y^3 + x^3*y + x
        """
        return "Set of places of {}".format(self._field)

    def _element_constructor_(self, ideal):
        """
        Create a place from the prime ``ideal``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: places = L.place_set()
            sage: O = L.maximal_order()
            sage: places(O.ideal(x,y))
            Place (x, y)
        """
        if not ideal.is_prime():
            raise TypeError("not a prime ideal")

        return self.element_class(self, ideal)

    def _an_element_(self):
        """
        Return a place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: places = L.place_set()
            sage: places.an_element()  # random
            Ideal (x) of Maximal order of Rational function field in x
            over Finite Field of size 2
        """
        d = 1
        while True:
            try:
                p = self._field.places(d).pop()
            except IndexError:
                d = d + 1
            else:
                break
        return p

    def function_field(self):
        """
        Return the function field to which this place set belongs.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: PS = L.place_set()
            sage: PS.function_field() == L
            True
        """
        return self._field
