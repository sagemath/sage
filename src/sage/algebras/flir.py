r"""
FLIR (Finite Intersection Laurent Rings): a class for working with FLIRs and their elements.

A FLIR is represented by a base chart (a Laurent polynomial ring over a field)
and a list of additional charts, each defined by explicit substitutions from the base chart.
The class supports operations on FLIR elements, including factorization into prime divisors.
It also supports computations of divisors and class groups.

This file provides two main classes: 

- :class: `FLIR`

- :class: `FLIRElement`

:class: `FLIR` represents a FLIR, which is a finite intersection of Laurent polynomial rings over a field. 
It provides, besides all the algebraic features,methods for working with the divisor group, class group, and prime divisors. 
A FLIR is represented by a base chart and a list of additional charts, each defined by explicit substitutions from the base chart.

:class: `FLIRElement` represents an element of a FLIR. It provides methods for factorization into irreducibles and for computing the divisor of an element.


Auxiliary classes:
- :class: `FLIRChart` -- a Laurent polynomial ring with a specified fraction field and birational maps.

- :class: `FLIRPrimeDivisor` -- a height-one prime ideal of a FLIR, represented by a chart and an irreducible polynomial.

- :class: `FLIRDivisor` -- an element of the divisor group of a FLIR, represented as a formal sum of prime divisors.

- :class: `FLIRDivisorGroup` -- the free abelian group on the prime divisors of a FLIR.

- :class: `ClassGroupData` -- a data structure for storing class group information, including the group itself, the map from divisors to classes, and the list of prime divisors.

- :class: `FLIRFactorization` -- a data structure for storing the factorization of an element of a FLIR into irreducibles.

        
      
REFERENCES:
- Mara Pompili and Daniel Smertnig, "Factoriality and Class Groups of Upper Cluster Algebras
  and Finite Laurent Intersection Rings: A Computational Approach", 2026. arXiv:2601.07520.

AUTHORS: 
- Mara Pompili (2026-06-25): initial version
- Daniel Smertnig (2026-06-25): initial version

EXAMPLES:

We begin by creating a base chart and a simple birational chart defined by `y_1 = x_1` and `y_2 = x_2/x_1`::

    sage: K = QQ
    sage: C0 = FLIRChart.base(K, ("x1","x2"))
    sage: x1, x2 = C0.F.gens()
    sage: C1 = FLIRChart(K, ("y1","y2"), base_fraction_field=C0.F,
    ....:              this_to_base=[x1, x2/x1])
    sage: C1.compute_base_to_this(C0)
    sage: C1.base_to_this
    [y1, y1*y2]

Then a FLIR can be created from the base chart and the additional chart::

    sage: A = FLIR(C0, [C1])
    sage: A
    FLIR over Rational Field 
      rank n = 2
      #charts = 1
      base vars = ('x1', 'x2')


Simple operations on FLIR behave as expected::

    sage: f = x1 + x2
    sage: g = x1*x2
    sage: A(f + g)
    x1*x2 + x1 + x2
    sage: A(f * g)
    x1^2*x2 + x1*x2^2

More attention must be paid to division, as the result may not be in the FLIR::

    sage: A(f / g)
    (x1 + x2)/(x1*x2)
    sage: A(g / f)
    Traceback (most recent call last):
    ...
    ValueError: Not in FLIR: not Laurent in chart ('y1', 'y2').
    Substituted: y1*y2/(y2 + 1)

We can compute the divisor class group of a FLIR::
    sage: A.class_group()
    Trivial Abelian group
    sage: A3 = example_A3()
    sage: A3.class_group()
    Multiplicative Abelian group isomorphic to Z

We can also compute the divisor of an element of the FLIR::

    sage: z1, z2, z3 = A3._base_gens()
    sage: z = A3(z1*z3)
    sage: A3.divisor(z)
    1*PrimeDivisor(chart=('x_21', 'x_22', 'x_23'), p=x_22 + 1) +
     1*PrimeDivisor(chart=('x_31', 'x_32', 'x_33'), p=x_32 + 1) +
     2*PrimeDivisor(chart=('x_41', 'x_42', 'x_43'), p=x_42 + 1)
    

and the image of the divisor in the class group::

    sage: A3.divisor_class(z)
    (0)

    

One can also compute the factorization of an element of the FLIR into irreducibles. 

    sage: w = A3(z2 + 1)
    sage: w.factor()
    [[((x2 + 1)/x1, 1), (x1, 1)], [((x2 + 1)/x3, 1), (x3, 1)]]
    sage: w.atoms()
    [x3, x1, (x2 + 1)/x1, (x2 + 1)/x3]

Moreover, the factorization can be represented in a verbose way::

    sage: w.factor(verbose=True)
    Factorizations of x2 + 1:
    Atoms: [x3, x1, (x2 + 1)/x1, (x2 + 1)/x3]
    Number of factorizations: 2
      1: 1 * ((x2 + 1)/x1) * (x1)
      2: 1 * ((x2 + 1)/x3) * (x3)

"""

from typing import Any, Dict,List
from dataclasses import dataclass

from sage.categories.additive_groups import AdditiveGroups
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.matrix.constructor import matrix
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.misc.latex import latex 
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element import Element, CommutativeAlgebraElement
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp, op_EQ, op_NE
from sage.structure.unique_representation import UniqueRepresentation

# -----------------------------------------------------------------------------
# Auxiliary Classes
# -----------------------------------------------------------------------------

class FLIRChart:
    r"""
    A birational chart for a fixed base fraction field.

    A chart consists of a Laurent polynomial ring
    `L = K[y_1^{\\pm 1}, \\dots, y_n^{\\pm 1}]` with fraction field `F = \\mathrm{Frac}(L)`,
    together with (optionally) birational maps to/from a chosen *base* fraction field.

    DATA:

    - ``this_to_base`` -- list of length ``n``; the images of the chart generators
      `y_i` in the base fraction field (typically `K(x_1,\\dots,x_n)`)

    - ``base_to_this`` -- list of length ``n``; the images of the base generators
      `x_i` in the chart fraction field `F`

    INPUT:

    - ``K`` -- a ring/field; the coefficient ring

    - ``var_names`` -- iterable of strings; names of the chart variables

    - ``term_order`` -- string (default: ``"lex"``); term order used for polynomial
      and Groebner computations

    - ``base_fraction_field`` -- a fraction field (optional); the base field in which
      ``this_to_base`` lives

    - ``base_to_this`` -- list (optional); images of base generators in the chart field

    - ``this_to_base`` -- list (optional); images of chart generators in the base field


    EXAMPLES::

        sage: K = QQ
        sage: C0 = FLIRChart.base(K, ("x1","x2"))
        sage: C0
        FLIRChart(('x1', 'x2'))over Rational Field with substitutions:
          x1 -> x1
          x2 -> x2

    A simple birational chart given by `y_1 = x_1` and `y_2 = x_2/x_1`::

        sage: x1, x2 = C0.F.gens()
        sage: C1 = FLIRChart(K, ("y1","y2"), base_fraction_field=C0.F,
        ....:              this_to_base=[x1, x2/x1])
        sage: C1.compute_base_to_this(C0)
        sage: C1.base_to_this
        [y1, y1*y2]

        sage: C1._substitute_from_base(x2/x1 + x1)
        y1 + y2

    TESTS::

        sage: K = QQ
        sage: C0 = FLIRChart.base(K, ("x1","x2"))
        sage: x1, x2 = C0.F.gens()
        sage: C1 = FLIRChart(K, ("y1","y2"), base_fraction_field=C0.F,
        ....:              this_to_base=[x1, x2/x1])
        sage: C1.compute_base_to_this(C0)
        sage: phi = C0.F.hom(C1.base_to_this, codomain=C1.F)      # x -> g(y)
        sage: psi = C1.F.hom(C1.this_to_base, codomain=C0.F)      # y -> f(x)
        sage: [psi(phi(z)) for z in C0.F.gens()] == list(C0.F.gens())
        True
        sage: [phi(psi(z)) for z in C1.F.gens()] == list(C1.F.gens())
        True
    """



    def __init__(
        self,
        K,
        var_names,
        term_order="lex",
        base_fraction_field=None,
        base_to_this=None,
        this_to_base=None,
    ):
        
        r"""
        Initialize a chart.

        INPUT:

        - ``K`` -- a field; the coefficient ring

        - ``var_names`` -- iterable of strings; names for chart variables

        - ``term_order`` -- string (default: ``"lex"``); monomial order

        - ``base_fraction_field`` -- (optional) base fraction field

        - ``base_to_this`` -- (optional) list of images of base generators in this chart

        - ``this_to_base`` -- (optional) list of images of chart generators in the base field

        TESTS::

            sage: K = QQ
            sage: C = FLIRChart(K, ("y1","y2"))
            sage: C.n
            2
            sage: C.base_fraction_field is None
            True
        """
        self.base_ring = K
        self.var_names = tuple(var_names)
        self.n = len(self.var_names)
        self.term_order = term_order

        self.L = LaurentPolynomialRing(K, self.n, names=self.var_names, order=term_order)
        self.F = self.L.fraction_field()
        self.P = PolynomialRing(K, self.n, names=self.var_names, order=term_order)

        self.base_fraction_field = base_fraction_field
        self.base_to_this = list(base_to_this) if base_to_this is not None else None
        self.this_to_base = list(this_to_base) if this_to_base is not None else None

    @classmethod
    def base(cls, K, var_names, term_order="lex"):
        """
        Create the base chart. Its "base_fraction_field" is its own fraction field,
        and both maps are the identity.
        INPUT:

        - ``K`` -- a ring/field

        - ``var_names`` -- iterable of strings

        - ``term_order`` -- string (default: ``"lex"``)


        EXAMPLES::

            sage: C = FLIRChart.base(QQ, ("x","y"))
            sage: C.base_fraction_field == C.F
            True
            sage: C.base_to_this
            [x, y]
            sage: C.this_to_base
            [x, y]
        """
        chart = cls(K, var_names, term_order=term_order)
        chart.base_fraction_field = chart.F
        chart.base_to_this = list(chart.F.gens())
        chart.this_to_base = list(chart.F.gens())
        return chart

    def __eq__(self, other):
        if not isinstance(other, FLIRChart):
            return NotImplemented
        return (
            self.base_ring == other.base_ring
            and self.var_names == other.var_names
            and self.term_order == other.term_order
            and self.this_to_base == other.this_to_base
            and self.base_to_this == other.base_to_this
            and self.base_fraction_field == other.base_fraction_field
        )

    def __hash__(self):
        return hash((
            self.base_ring,
            self.var_names,
            self.term_order,
            tuple(self.this_to_base) if self.this_to_base is not None else None,
            tuple(self.base_to_this) if self.base_to_this is not None else None,
            self.base_fraction_field,
        ))
    
    def _substitute_from_base(self, f_in_base):
        """
        Substitute a base expression into this chart using ``base_to_this``.
        
        INPUT:

        - ``f_in_base`` -- an element of ``self.base_fraction_field``

        OUTPUT:

        The image of ``f_in_base`` in ``self.F`` obtained by substituting base generators
        `x_i` with ``self.base_to_this[i]``.

        Raises an error if no base substitution data is stored.

        EXAMPLES::

            sage: K = QQ
            sage: C0 = FLIRChart.base(K, ("x1","x2"))
            sage: x1, x2 = C0.F.gens()
            sage: C1 = FLIRChart(K, ("y1","y2"),
            ....:              base_fraction_field=C0.F,
            ....:              this_to_base=[x1, x2/x1])
            sage: C1.compute_base_to_this(C0)
            sage: C1._substitute_from_base(x2/x1 + x1)
            y1 + y2

        TESTS::
            sage: C = FLIRChart(QQ, ("y1","y2"))
            sage: C._substitute_from_base(1)
            Traceback (most recent call last):
            ...
            ValueError: No base substitution data stored for this chart.
        """
        
        if self.base_fraction_field is None or self.base_to_this is None:
            raise ValueError("No base substitution data stored for this chart.")
        phi = self.base_fraction_field.hom(self.base_to_this, codomain=self.F)
        return phi(f_in_base)

    def _laurent_poly_to_poly_up_to_unit(self, h):
        r"""
        Clear negative exponents in a Laurent polynomial up to a monomial unit.

        INPUT:

        - ``h`` -- a nonzero element of ``self.L`` (a Laurent polynomial)

        OUTPUT:

        A pair ``(g, u)`` where:

        - ``g`` is a polynomial in ``self.P``,
        - ``u`` is a Laurent monomial in ``self.L``,

        such that `g = u \\cdot h` in the Laurent ring.

        EXAMPLES::

            sage: C = FLIRChart.base(QQ, ("y1","y2"))
            sage: y1, y2 = C.L.gens()
            sage: h = y1^-2 * y2 + y1^-1
            sage: g, u = C._laurent_poly_to_poly_up_to_unit(h)
            sage: C.L(u) * C.L(h) == C.L(g)
            True
            sage: g
            y1 + y2
        
        TESTS::
            sage: C = FLIRChart.base(QQ, ("y1","y2"))
            sage: C._laurent_poly_to_poly_up_to_unit(0)
            Traceback (most recent call last):
            ...
            ValueError: Divisor of 0 is undefined.
        """
        

        hL = self.L(h)
        if hL == 0:
            raise ValueError("Divisor of 0 is undefined.")

        mons = hL.monomials()
        exps = [m.exponents()[0] for m in mons]
        min_exp = [min(e[i] for e in exps) for i in range(self.n)]

        shift = self.L.monomial(tuple([-e for e in min_exp]))
        gL = shift * hL
        gP = self.P(gL)  # exponents now >= 0
        return gP, shift

    def compute_base_to_this(self, base_chart: "FLIRChart"):
        r"""
        Compute ``base_to_this`` as the inverse of ``this_to_base`` by elimination.

        INPUT:

        - ``base_chart`` -- a :class:`FLIRChart`; the base chart. It must satisfy
        ``base_chart.F == self.base_fraction_field`` and have the same dimension.

        OUTPUT:

        This method sets ``self.base_to_this`` to a list of length ``n`` of elements of
        ``self.F``, giving the base generators expressed as rational functions in the
        chart generators.

        ALGORITHM:

        Let `y_j = f_j(x)` be the map stored in ``this_to_base``.
        The method clears denominators and forms polynomial equations
        in `K[x_1,\\dots,x_n,y_1,\\dots,y_n]`:
        `y_j * \\mathrm{den}(f_j) - \\mathrm{num}(f_j) = 0`.
        A Groebner basis with lex order is computed, and linear relations are extracted
        to solve for each `x_i` in terms of the `y`'s. The result is verified by
        composing the two maps on generators.

        EXAMPLES::
            sage: K = QQ
            sage: C0 = FLIRChart.base(K, ("x1","x2"))
            sage: x1, x2 = C0.F.gens()
            sage: C1 = FLIRChart(K, ("y1","y2"), base_fraction_field=C0.F,
            ....:              this_to_base=[x1, x2/x1])
            sage: C1.compute_base_to_this(C0)
            sage: C1.base_to_this
            [y1, y1*y2]

        TESTS::
            sage: K = QQ
            sage: C0 = FLIRChart.base(K, ("x1","x2"))
            sage: x1, x2 = C0.F.gens()
            sage: Cbad = FLIRChart(K, ("y1","y2"), base_fraction_field=C0.F)
            sage: Cbad.compute_base_to_this(C0)
            Traceback (most recent call last):
            ...
            ValueError: Cannot invert: this_to_base is not set.

            sage: C0b = FLIRChart.base(K, ("x1","x2","x3"))
            sage: C1 = FLIRChart(K, ("y1","y2"), base_fraction_field=C0.F,
            ....:              this_to_base=[x1, x2/x1])
            sage: C1.compute_base_to_this(C0b)
            Traceback (most recent call last):
            ...
            ValueError: Dimension mismatch between base and chart.
        """

        if self.this_to_base is None:
            raise ValueError("Cannot invert: this_to_base is not set.")
        if base_chart.n != self.n:
            raise ValueError("Dimension mismatch between base and chart.")
        if self.base_fraction_field != base_chart.F:
            raise ValueError("Base fraction field mismatch.")

        K = base_chart.base_ring
        n = self.n
        xnames = list(base_chart.var_names)
        ynames = list(self.var_names)

        R = PolynomialRing(K, xnames + ynames, order="lex")
        xs = R.gens()[:n]
        ys = R.gens()[n:]

        phi_x = base_chart.P.hom(list(xs), R)

        eqs = []
        for j in range(n):
            fj = self.this_to_base[j]  # in base fraction field
            num = fj.numerator()
            den = fj.denominator()

            numL = base_chart.L(num)
            denL = base_chart.L(den)

            mons = list(numL.monomials()) + list(denL.monomials())
            exps = [m.exponents()[0] for m in mons]
            min_exp = [min(e[i] for e in exps) for i in range(n)]

            shift = base_chart.L.monomial(tuple([-e for e in min_exp]))
            numP = base_chart.P(shift * numL)
            denP = base_chart.P(shift * denL)

            eqs.append(ys[j] * phi_x(denP) - phi_x(numP))

        G = R.ideal(eqs).groebner_basis()

        inv_x = [None] * n
        xs_list = list(xs)

        for g in G:
            vars_in = g.variables()
            x_in = [v for v in vars_in if v in xs_list]
            if len(x_in) != 1:
                continue
            x = x_in[0]
            if g.degree(x) != 1:
                continue
            if any(v in xs_list and v != x for v in vars_in):
                continue

            a = g.coefficient({x: 1})
            if a == 0:
                continue
            b = g.subs({x: 0})
            i = xs_list.index(x)
            inv_x[i] = -b / a

        if any(v is None for v in inv_x):
            missing = [xnames[i] for i, v in enumerate(inv_x) if v is None]
            raise ValueError(
                f"Failed to invert chart embedding: could not solve for base vars {missing}."
            )

        yF = self.F.gens()
        psi_poly = R.hom([self.F(0)] * n + list(yF), self.F)
        RF = R.fraction_field()

        def psi_frac(expr):
            expr = RF(expr)
            return psi_poly(expr.numerator()) / psi_poly(expr.denominator())

        self.base_to_this = [psi_frac(expr) for expr in inv_x]

        # sanity checks
        phi_base_to_chart = base_chart.F.hom(self.base_to_this, codomain=self.F)  # x -> g(y)
        phi_chart_to_base = self.F.hom(self.this_to_base, codomain=base_chart.F)  # y -> f(x)

        for xi in base_chart.F.gens():
            if phi_chart_to_base(phi_base_to_chart(xi)) != xi:
                raise ValueError("Inversion check failed on base generators.")
        for yj in self.F.gens():
            if phi_base_to_chart(phi_chart_to_base(yj)) != yj:
                raise ValueError("Inversion check failed on chart generators.")
    
    def _repr_(self):

        if self.this_to_base is not None:
            lines = [f"FLIRChart({self.var_names})over {self.base_ring} with substitutions:"]
            for name, img in zip(self.var_names, self.this_to_base):
                lines.append(f"  {name} -> {img}")
            return "\n".join(lines)

        else:
            return f"FLIRChart({self.var_names}) over {self.base_ring}. No substitution data."



    __repr__ = _repr_


@dataclass(frozen=True) 
class FLIRPrimeDivisor:
    r"""
    A prime divisor (i.e., height-one prime) used in the divisor ideals group of a FLIR.

    A prime divisor is represented by:

    - ``chart`` -- a :class:`FLIRChart`; the chart in which the prime is defined

    - ``irreducible`` -- an irreducible polynomial in ``chart.P`` (typically normalized
      to have leading coefficient 1)


    EXAMPLES::
        sage: K = QQ
        sage: C = FLIRChart.base(K, ("x1","x2"))
        sage: x1, x2 = C.P.gens()
        sage: P = FLIRPrimeDivisor(C, x1 + 1)
        sage: P
        PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1)

    TESTS::
        sage: K = QQ
        sage: C = FLIRChart.base(K, ("x1","x2"))
        sage: x1, x2 = C.P.gens()
        sage: P1 = FLIRPrimeDivisor(C, x1 + 1)
        sage: P2 = FLIRPrimeDivisor(C, x1 + 1)
        sage: P1 == P2
        True
        sage: d = {P1: 3}
        sage: d[P2]
        3
    """
    chart: FLIRChart
    irreducible: Any

    def __repr__(self):
        return f"PrimeDivisor(chart={self.chart.var_names}, p={self.irreducible})"



class FLIRDivisor(Element):
    r"""
    An element of the divisor group ``Div(A)``.

    A divisor is a finitely supported formal `\\ZZ`-linear combination of prime divisors.
    It is stored as a dictionary ``{P: e}`` where:

    - ``P`` is a :class:`FLIRPrimeDivisor`
    - ``e`` is an integer in ``ZZ``

    INPUT:

    - ``parent`` -- a :class:`FLIRDivisorGroup`

    - ``data`` -- (optional) either a dict ``{prime: exponent}`` or an iterable of
        pairs ``(prime, exponent)``

    Zero coefficients are discarded automatically.

    EXAMPLES::
        sage: K = QQ
        sage: C = FLIRChart.base(K, ("x1","x2"))
        sage: x1, x2 = C.P.gens()
        sage: P = FLIRPrimeDivisor(C, x1 + 1)
        sage: Q = FLIRPrimeDivisor(C, x2 + 1)
        sage: # Build a tiny fake divisor group just for examples:
        sage: class _DummyA: pass
        sage: G = FLIRDivisorGroup(_DummyA())
        sage: D = G({P: 2, Q: -1}); D
        2*PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1) +
         -1*PrimeDivisor(chart=('x1', 'x2'), p=x2 + 1)
        sage: D.support() == {P, Q}
        True
        sage: D.coeff(P), D.coeff(Q)
        (2, -1)

    Arithmetic in ``Div(A)``::
        sage: E = G({P: 1})
        sage: D + E
        3*PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1) +
         -1*PrimeDivisor(chart=('x1', 'x2'), p=x2 + 1)
        sage: D - E
        1*PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1) + 
         -1*PrimeDivisor(chart=('x1', 'x2'), p=x2 + 1)
        sage: -E
        -1*PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1)

    TESTS::
        sage: K = QQ
        sage: C = FLIRChart.base(K, ("x1","x2"))
        sage: x1, x2 = C.P.gens()
        sage: P = FLIRPrimeDivisor(C, x1 + 1)
        sage: class _DummyA: pass
        sage: G = FLIRDivisorGroup(_DummyA())
        sage: (G({P: 1}) + G({P: -1}))
        0
    """

    def __init__(self, parent, data=None):
        super().__init__(parent)
        d = {}
        if data is None:
            self._data = d
            return

        items = data.items() if isinstance(data, dict) else data
        for P, e in items:
            e = ZZ(e)
            if e == 0:
                continue
            d[P] = d.get(P, ZZ(0)) + e
            if d[P] == 0:
                del d[P]
        self._data = d

    def _repr_(self):
        if not self._data:
            return "0"
        parts = [f"{e}*{P}" for P, e in self._data.items()]
        return " +\n ".join(parts)

    def __iter__(self):
        return iter(self._data.items())

    def _richcmp_(self, other, op):
        r"""
        Compare two divisors.

        Divisors are equal exactly when they have the same parent and the same
        prime coefficients. The insertion order of the dictionary does not
        matter.

        EXAMPLES::

            sage: K = QQ
            sage: C = FLIRChart.base(K, ("x1","x2"))
            sage: x1, x2 = C.P.gens()
            sage: P = FLIRPrimeDivisor(C, x1 + 1)
            sage: Q = FLIRPrimeDivisor(C, x2 + 1)
            sage: class _DummyA: pass
            sage: G = FLIRDivisorGroup(_DummyA())
            sage: G({P: 1}) == G({P: 1})
            True
            sage: G({P: 1, Q: 2}) == G({Q: 2, P: 1})
            True
            sage: G({P: 1}) == G({P: 2})
            False
        """
        if op == op_EQ:
            return self.parent() is other.parent() and self._data == other._data
        if op == op_NE:
            return self.parent() is not other.parent() or self._data != other._data
        return NotImplemented

    def support(self):
        r"""
        Return the support of this divisor as a set of primes.

        EXAMPLES::
            sage: K = QQ
            sage: C = FLIRChart.base(K, ("x1","x2"))
            sage: x1, x2 = C.P.gens()
            sage: P = FLIRPrimeDivisor(C, x1 + 1)
            sage: class _DummyA: pass
            sage: G = FLIRDivisorGroup(_DummyA())
            sage: D = G({P: 3})
            sage: D.support()
            {PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1)}
        """
        return set(self._data)

    def coeff(self, P):
        r"""
        Return the coefficient of the prime ``P`` in this divisor.

        EXAMPLES::
            sage: C = FLIRChart.base(QQ, ("x1","x2"))
            sage: x1, x2 = C.P.gens()
            sage: P = FLIRPrimeDivisor(C, x1 + 1)
            sage: class _DummyA: pass
            sage: G = FLIRDivisorGroup(_DummyA())
            sage: D = G({P: 3})
            sage: D.coeff(P)
            3
        """
        return self._data.get(P, ZZ(0))
    
    def coeffs(self):
        r"""
        Return the list of coefficients of this divisor.

        EXAMPLES::
            sage: C = FLIRChart.base(QQ, ("x1","x2"))
            sage: x1, x2 = C.P.gens()
            sage: x1, x2 = C.P.gens()
            sage: P = FLIRPrimeDivisor(C, x1 + 1)
            sage: Q = FLIRPrimeDivisor(C, x2 + 1)
            sage: class _DummyA: pass
            sage: G = FLIRDivisorGroup(_DummyA())
            sage: D = G({P: 3, Q: 5})
            sage: D.coeffs()
            [3, 5]

        """
        return list(self._data.values())

    def _add_(self, other):
        out = dict(self._data)
        for P, e in other._data.items():
            out[P] = out.get(P, ZZ(0)) + e
            if out[P] == 0:
                del out[P]
        return self.parent()(out)

    def _neg_(self):
        return self.parent()({P: -e for P, e in self._data.items()})

    def _sub_(self, other):
        return self._add_(other._neg_())
    
    def _scalar_mul(self, n):
        n = ZZ(n)
        if n == 0:
            return self.parent().zero()
        return self.parent()({P: n*e for P, e in self._data.items()})

    def _lmul_(self, n):
        return self._scalar_mul(n)

    def _rmul_(self, n):
        return self._scalar_mul(n)
    
    def is_effective(self):
        r"""
        Return ``True`` if all coefficients are nonnegative.

        EXAMPLES::
            sage: K = QQ
            sage: C = FLIRChart.base(K, ("x1","x2"))
            sage: x1, x2 = C.P.gens()
            sage: P = FLIRPrimeDivisor(C, x1 + 1)
            sage: class _DummyA: pass
            sage: G = FLIRDivisorGroup(_DummyA())
            sage: G({P: 2}).is_effective()
            True
            sage: G({P: -1}).is_effective()
            False
        """
        return all(e >= 0 for e in self._data.values())
    
    def is_principal(self):
        r"""
        Return whether this divisor is principal in the ambient FLIR.

        EXAMPLES::
            sage: A = example_A3()
            sage: C = A.charts
            sage: c = C[0]
            sage: x1, x2, x3 = c.P.gens()
            sage: P = FLIRPrimeDivisor(c, x1 + 1)
            sage: Q = A.extra_primes()
            sage: q = Q[3]
            sage: D = A.divisor_group()({P: 1, q: 2})
            sage: D.is_principal()
            False
            sage: D2 = A.divisor_group()({P: 1,})
            sage: D2.is_principal()
            True

        .. WARNING::

            This method requires that the parent divisor group was created by a
            :class:`FLIR` instance, i.e. ``self.parent().flir()`` must exist.
        """
        return self.parent().flir().is_principal_divisor(self)
    
    def gen(self):
        r"""
        Return a principal generator if this divisor is principal.

        Raises ``ValueError`` if the divisor is not principal.

        EXAMPLES::
            sage: A = example_A3()
            sage: C = A.charts
            sage: c = C[0]
            sage: x1, x2, x3 = c.P.gens()
            sage: P = FLIRPrimeDivisor(c, x2 + 1)
            sage: Q = A.extra_primes()[2]
            sage: D = A.divisor_group()({P: 1, Q: 1})
            sage: D.is_principal()
            True
            sage: D.gen()
            (x2 + 1)/x1

        .. WARNING::

            This method requires that the parent divisor group was created by a
            :class:`FLIR` instance, i.e. ``self.parent().flir()`` must exist.
        """
        if not self.is_principal():
            raise ValueError("Not a principal divisor.")
        return self.parent().flir().principal_generator(self)

class FLIRDivisorGroup(Parent):
    r"""
    The divisor group ``Div(A)`` of a FLIR.

    This is the free abelian group on the height-1 prime ideals of `A`,
    implemented as finitely supported formal sums.

    Elements are instances of :class:`FLIRDivisor`, represented internally by a dict
    ``{prime: exponent}`` with coefficients in ``ZZ``.

    INPUT:

    - ``flir`` -- the ambient :class:`FLIR` (or an object providing the same interface)

    - ``base_ring`` -- (default: ``ZZ``) the coefficient ring for divisor coefficients

    EXAMPLES::
        sage: class _DummyA: pass
        sage: G = FLIRDivisorGroup(_DummyA()); G
        Divisor group Div(A) of <..._DummyA...> with coefficients in Integer Ring
        sage: G.zero()
        0

    Creating divisors::
        sage: K = QQ
        sage: C = FLIRChart.base(K, ("x1","x2"))
        sage: x1, x2 = C.P.gens()
        sage: P = FLIRPrimeDivisor(C, x1 + 1)
        sage: D = G.prime(P); D
        1*PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1)

        sage: G({P: 2})
        2*PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1)
    """

    Element = FLIRDivisor

    def __init__(self, flir, base_ring=ZZ):
        self._A = flir
        self._base_ring = base_ring
        Parent.__init__(self, category=AdditiveGroups()) # type: ignore[call-arg]

    def _repr_(self):
        return f"Divisor group Div(A) of {self._A} with coefficients in {self._base_ring}"

    def _element_constructor_(self, x=None):
        r"""
        Construct an element of ``self``.

        INPUT:

        - ``x`` -- ``None`` (zero divisor), a dict, an iterable of pairs, or an existing
          :class:`FLIRDivisor` in this parent.
        """
        if x is None or x == 0:
            return self.element_class(self, {})
        if isinstance(x, FLIRDivisor) and x.parent() is self:
            return x
        # allow dict / list of pairs
        return self.element_class(self, x)

    def zero(self):
        r"""
        Return the zero divisor.

        EXAMPLES::
            sage: class _DummyA: pass
            sage: G = FLIRDivisorGroup(_DummyA())
            sage: G.zero()
            0
        """
        return self()

    def prime(self, P):
        r"""
        Return the divisor `1\\cdot P`.

        INPUT:

        - ``P`` -- a :class:`FLIRPrimeDivisor`

        EXAMPLES::
            sage: K = QQ
            sage: C = FLIRChart.base(K, ("x1","x2"))
            sage: x1, x2 = C.P.gens()
            sage: P = FLIRPrimeDivisor(C, x1 + 1)
            sage: class _DummyA: pass
            sage: G = FLIRDivisorGroup(_DummyA())
            sage: G.prime(P)
            1*PrimeDivisor(chart=('x1', 'x2'), p=x1 + 1)
        """
        return self({P: 1})

    def flir(self):
        """Return the ambient FLIR object associated to this divisor group."""
        return self._A
    
class FLIRFactorization:
    r"""
    Formal factorizations of an element of a FLIR.

    Conventions:
      - self._element: a FLIRElement
      - self._atoms: list of FLIRElement that are atoms
      - self._factorizations: list of factorizations, each factorization is a list of (atom, exponent)

    EXAMPLES:

    A non-trivial example with two distinct factorizations, illustrating that
    a FLIR need not have unique factorization::

        sage: A = example_A2_generalized()
        sage: x1, x2 = A._base_gens()
        sage: f = A((x1 + 1)**2)
        sage: F = FLIRFactorization.from_element(f, verbose=False)
        sage: F.atoms()
        [x2, x1 + 1, (x1^2 + 2*x1 + 1)/x2]
        sage: F.num_atoms()
        3
        sage: F.num_factorizations()
        2
        sage: F.is_irreducible()
        False
        sage: F.sort_all_factorizations()
        [[((x1^2 + 2*x1 + 1)/x2, 1), (x2, 1)], [(x1 + 1, 2)]]

    The verbose representation displays both factorizations::

        sage: Fv = FLIRFactorization.from_element(f, verbose=True)
        sage: print(Fv)
        Factorizations of x1^2 + 2*x1 + 1:
        Atoms: [x2, x1 + 1, (x1^2 + 2*x1 + 1)/x2]
        Number of factorizations: 2
          1: 1 * ((x1^2 + 2*x1 + 1)/x2) * (x2)
          2: 1 * (x1 + 1)^2

    Other queries on this factorization::

        sage: Fv.units()
        [1, 1]
        sage: Fv.set_of_lengths()
        {2}
        sage: Fv.number_atoms_in_factorizations()
        [1, 2]
        sage: Fv.particular_factorization(0)
        [(x1 + 1, 2)]
        sage: Fv.particular_factorization(1)
        [(x2, 1), ((x1^2 + 2*x1 + 1)/x2, 1)]
        sage: Fv.sorted_factorizations_with_units()
        [(1, [((x1^2 + 2*x1 + 1)/x2, 1), (x2, 1)]), (1, [(x1 + 1, 2)])]

    An irreducible element has a unique factorization consisting of a single
    atom with exponent one::

        sage: g = A(x1 + x2)
        sage: G = FLIRFactorization.from_element(g, verbose=False)
        sage: G.atoms()
        [x1 + x2]
        sage: G.num_atoms()
        1
        sage: G.num_factorizations()
        1
        sage: G.is_irreducible()
        True
        sage: G.is_unit()
        False
        sage: G
        [[(x1 + x2, 1)]]
        sage: Gv = FLIRFactorization.from_element(g, verbose=True)
        sage: print(Gv)
        x1 + x2 is irreducible

    A unit has no atoms at all::

        sage: u = A.one()
        sage: U = FLIRFactorization.from_element(u, verbose=False)
        sage: U.num_atoms()
        0
        sage: U.is_unit()
        True
        sage: U
        [[]]

    """

    def __init__(self, element, atoms, factorizations, units=None, verbose=False):
        self._element = element
        self._atoms = list(atoms)
        self._factorizations = [
            [(a, int(e)) for (a, e) in fac if int(e) != 0]
            for fac in factorizations
        ]

        one = element.parent().fraction_field(1)

        if units is None:
            self._units = [one] * len(self._factorizations)
        elif isinstance(units, (list, tuple)):
            if len(units) != len(self._factorizations):
                raise ValueError("number of units must match number of factorizations")
            self._units = list(units)
        else:
            self._units = [units] * len(self._factorizations)

        self._verbose = verbose

    @classmethod
    def from_element(cls, f, verbose=False):
        """
        Build FLIRFactorization from a FLIRElement 
        """
        A = f.parent()
        atoms, factorizations, units = A.factorizations(f)
        return cls(f, atoms=atoms, factorizations=factorizations, units=units, verbose=verbose)
    

    def element(self):
        return self._element

    def atoms(self):
        return self._atoms
    
    def units(self):
        return self._units

    def has_nontrivial_units(self):
        one = self._element.parent().fraction_field(1)
        return any(u != one for u in self._units)

    def factorizations(self):
        return self._factorizations
    
    def factorization(self):
        """
        Return the factorization if it is unique, otherwise raise an error.
        """
        if self.num_factorizations() == 1:
            return self._factorizations[0]
        raise ValueError("Not a unique factorization.")

    def num_atoms(self):
        return len(self._atoms)

    def num_factorizations(self):
        return len(self._factorizations)

    def is_unit(self):
        return self.num_atoms() == 0

    def is_irreducible(self):
        # irreducible if exactly one atom and exactly one factorization with exponent 1
        return (
            self.num_atoms() == 1
            and self.num_factorizations() == 1
            and len(self._factorizations[0]) == 1
            and self._factorizations[0][0][1] == 1
        )

    def set_representation(self, verbose=True):
        self._verbose = bool(verbose)

    def set_of_lengths(self):
        """
        Length = sum of exponents in a factorization.
        """
        return {sum(exp for _, exp in fac) for fac in self._factorizations}

    def number_atoms_in_factorizations(self):
        """
        Number of distinct atoms appearing in each factorization.
        """
        return [len(fac) for fac in self._factorizations]

    def particular_factorization(self, i):
        """
        Return the i-th factorization as a list of (atom, exp).
        """
        return self._factorizations[i]

    # ---------------- sorting ----------------

    def sort_factors(self, key=None):
        """
        Sort the factors inside each factorization.
        Default: by repr(atom), then exponent.
        """
        sorted_fact_list = []
        for fac in self._factorizations:
            fac2 = list(fac)
            if key is None:
                fac2.sort(key=lambda t: (repr(t[0]), t[1]))
            else:
                fac2.sort(key=key)
            sorted_fact_list.append(fac2)
        return sorted_fact_list

    def sorted_factorizations_with_units(self):
        """
        Sort factorizations globally, preserving their units.

        Sorting convention:
          - more distinct factors first
          - then lexicographically by repr(atom), exponent
        """
        pairs = [
            (u, fac)
            for u, fac in zip(self._units, self.sort_factors())
        ]

        def factorization_key(pair):
            _, fac = pair
            length = -len(fac)  # more distinct factors first
            lex = [(repr(a), e) for a, e in fac]
            return (length, lex)

        return sorted(pairs, key=factorization_key)

    def sort_all_factorizations(self):
        """
        Sort factorizations globally:
          - more distinct factors first
          - then lex by repr of factors.
        """
        return [fac for _, fac in self.sorted_factorizations_with_units()]

    # ---------------- display ----------------

    def __repr__(self):
        if not self._verbose:
            return repr(self.sort_all_factorizations())

        f = self._element

        if self.is_irreducible():
            return f"{f} is irreducible"
        if self.is_unit():
            return f"{f} is a unit"

        sorted_pairs = self.sorted_factorizations_with_units()

        one = f.parent().fraction_field(1)
        show_units = any(u != one for u in self._units)

        rep = f"Factorizations of {f}:\n"
        rep += f"Atoms: {self._atoms}\n"
        rep += f"Number of factorizations: {len(self._factorizations)}\n"

        for i, (u, fac) in enumerate(sorted_pairs, start=1):
            rhs = " * ".join(
                f"({a})^{e}" if e != 1 else f"({a})"
                for a, e in fac
            )

            if rhs == "":
                rhs = "1"

            if show_units:
                rep += f"  {i}: {u} * {rhs}\n"
            else:
                rep += f"  {i}: {rhs}\n"

        return rep

    # ---------------- LaTeX ----------------

    def latex(self, env="aligned", include_display_math=True):
        """
        LaTeX block for all factorizations.
        Uses Sage's latex() on the objects.
        """
        f_ltx = latex(self._element)

        if self.is_irreducible():
            out = rf"{f_ltx}\ \\text{{is irreducible}}"
            return rf"\[{out}\]" if include_display_math else out

        if self.is_unit():
            out = rf"{f_ltx}\ \\text{{is a unit}}"
            return rf"\[{out}\]" if include_display_math else out

        sorted_pairs = self.sorted_factorizations_with_units()
        show_units = self.has_nontrivial_units()

        lines = []
        for idx, (u, fac) in enumerate(sorted_pairs):
            pieces = []

            if show_units:
                pieces.append(latex(u))

            for a, e in fac:
                a_ltx = latex(a)
                if e == 1:
                    pieces.append(rf"\\left({a_ltx}\\right)")
                else:
                    pieces.append(rf"\\left({a_ltx}\\right)^{{{e}}}")

            rhs = r" \\cdot ".join(pieces) if pieces else latex(u)

            if idx == 0:
                lines.append(rf"{f_ltx} &= {rhs}")
            else:
                lines.append(rf" &= {rhs}")

        body = r"\\ ".join(lines)
        block = rf"\\begin{{{env}}}" + body + rf"\\end{{{env}}}"
        return r"\[" + block + r"\]" if include_display_math else block

class ClassGroupData:
    r"""
    Data and algorithms for `Cl(A)` presented as cokernel of `M: ZZ^n -> ZZ^r`,
    where M is the matrix of valuations of the extra primes `P_i` on the base generators `x_j`. 

    Internally, we have:
        -  Z^r is the free abelian group on the `P_i`, represented by `H = FreeModule(ZZ, r)`
        -  R = im(M) is the subgroup generated by the columns of M, represented by `R = H.submodule(gens)`
        -  Q = H/R is the quotient module representing `Cl(A)`, with quotient map `qmap: H -> Q`.

    The Smith normal form of `M` is used to:
        - compute the invariants of `Cl(A)` (free rank + torsion)
        - solve `M*x = rhs` over `ZZ` when needed for principality tests and generator construction.

    INPUT:

    - ``M``: an `r x n` integer matrix, where `r` = number of extra primes `P_i` and `n` = number of base generators `x_j`.
    - ``primes``: a list of the extra primes `P_i` corresponding to the rows of `M`.

    EXAMPLES::
        sage: M = matrix(ZZ, [[2, 0], [0, 3]])
        sage: primes = ['P1', 'P2']
        sage: data = ClassGroupData(M, primes)
        sage: data.Cl()
        Multiplicative Abelian group isomorphic to C6 
        
        sage: M = matrix(ZZ, [[2, 0], [0, 6]])
        sage: data = ClassGroupData(M, primes)
        sage: data.r, data.n
        (2, 2)
        sage: data.diag
        [2, 6]
        sage: data.Cl()
        Multiplicative Abelian group isomorphic to C2 x C6
    
    Mapping a vector `v \\in \\ZZ^r` to its class in the quotient module::
        sage: data.class_in_quotient([2, 0]) == 0 # (2,0) is in the image 
        True
        sage: data.class_in_quotient([1, 0]) == 0 
        False
    
    Solving `M x = b` over `\\ZZ` (strict by default)::
        sage: # A solvable, uniquely-solvable system (square full rank)
        sage: M = matrix(ZZ, [[1, 2],
        ....:                 [0, 3]])
        sage: data = ClassGroupData(M, primes=("P1","P2"))
        sage: b = vector(ZZ, [5, 6])
        sage: x = data.solve_Mx_eq_rhs(b)   # require_solution=True, require_unique=True by default
        sage: (M * x.column()).column(0) == b.column()
        False

    An inconsistent system ::
        sage: M = matrix(ZZ, [[2, 0],
        ....:                 [0, 6]])
        sage: data = ClassGroupData(M, primes=("P1","P2"))
        sage: data.solve_Mx_eq_rhs([1, 0])
        Traceback (most recent call last):
        ...
        ValueError: No integer solution to M*x = rhs.

    The same inconsistent system, but in permissive mode we get ``None``::

        sage: data.solve_Mx_eq_rhs([1, 0], require_solution=False) is None
        True
    
    TESTS::
        sage: M = matrix(ZZ, [[0, 0], [0, 0]])
        sage: data = ClassGroupData(M, primes=("P1","P2"))
        sage: data.diag
        []
        sage: data.Cl().invariants()
        (0, 0)
        sage: M = matrix(ZZ, [[2, 4]])      
        sage: data = ClassGroupData(M, primes=("P1",))
        sage: data.solve_Mx_eq_rhs([2])
        Traceback (most recent call last):
        ...
        ValueError: Infinitely many integer solutions (nontrivial kernel).

        sage: # In permissive mode we still get a valid solution
        sage: x = data.solve_Mx_eq_rhs([2], require_unique=False)
        sage: (M * x.column())[0,0]
        2

        sage: M = matrix(ZZ, [[0, 0],
        ....:                 [0, 0]])
        sage: data = ClassGroupData(M, primes=("P1","P2"))
        sage: data.solve_Mx_eq_rhs([0, 0])
        Traceback (most recent call last):
        ...
        ValueError: Infinitely many integer solutions (nontrivial kernel).

        sage: # But permissive mode returns one solution
        sage: data.solve_Mx_eq_rhs([0, 0], require_unique=False)
        (0, 0)
    """
    def __init__(self, M, primes):
        self.M = matrix(ZZ, M)                     # r x n integer matrix
        self.primes = tuple(primes)      
        self.r = M.nrows()
        self.n = M.ncols()

        self.H = FreeModule(ZZ, self.r) #H = Z^r

        # R = im(M) as a submodule of H, generated by columns of M
        # Each column is in ZZ^r -> convert to element of H
        
        gens = [list(self.M.column(j)) for j in range(self.n)]

        if gens:
            G = matrix(ZZ, gens)          
            E = G.echelon_form(algorithm='flint')

            nonzero = [
                self.H(list(E.row(i)))
                for i in range(E.nrows())
                if not E.row(i).is_zero()
            ]

            if nonzero:
                self.R = self.H.span(nonzero, already_echelonized=True)
            else:
                self.R = self.H.zero_submodule()
        else:
            self.R = self.H.zero_submodule()

        # The quotient module representing Cl(A)
        self.Q = self.H.quotient(self.R)
        # Map H -> Q
        self.qmap = self.Q.quotient_map()
    
        # SNF:
        self.S, self.U, self.V = self.M.smith_form(
           transformation=True,
        )

        self.diag = []
        m = min(self.S.nrows(), self.S.ncols())
        for k in range(m):
            d = ZZ(self.S[k, k])
            if d == 0:
                break
            self.diag.append(int(d))

        # invariants for AbelianGroup
        self._m = len(self.diag)
        self.torsion_indices = [k for k, d in enumerate(self.diag) if d != 1]
        self.torsion_moduli = [self.diag[k] for k in self.torsion_indices]

        invariant_factors = [d for d in self.diag if d != 1]
        free_rank = int(self.r - self._m)

        self._Cl = AbelianGroup([0] * free_rank + invariant_factors)

    def Cl(self):
        r"""
        Return the class group as a Sage :class:`~sage.groups.abelian_gps.abelian_group.AbelianGroup`.

        The returned group has invariant factors determined from the Smith normal form
        of ``self.M``.

        OUTPUT:

        A Sage :class:`~sage.groups.abelian_gps.abelian_group.AbelianGroup`.

        EXAMPLES::
            sage: M = matrix(ZZ, [[4, 0],
            ....:                 [0, 2]])
            sage: data = ClassGroupData(M, primes=("P1","P2"))
            sage: data.Cl().invariants()
            (2, 4)

        TESTS::
            sage: M = matrix(ZZ, [[1, 0],
            ....:                 [0, 1]])
            sage: data = ClassGroupData(M, primes=("P1","P2"))
            sage: # cokernel is trivial
            sage: data.Cl().order()
            1
        """
        return self._Cl


    def class_in_quotient(self, v):
        r"""
        Return the class of a vector in the cokernel module ``Q = H/R``.

        INPUT:

        - ``v`` -- iterable or vector; an element of `\\ZZ^r` (coordinates with respect
          to the ordering of ``self.primes``).

        OUTPUT:

        An element of the quotient module ``self.Q``.

        EXAMPLES::
            sage: M = matrix(ZZ, [[2],
            ....:                 [0]])
            sage: data = ClassGroupData(M, primes=("P1","P2"))
            sage: q = data.class_in_quotient([2,0])
            sage: q
            (0, 0)
            sage: data.class_in_quotient([1,0]) == 0
            False

        TESTS::
            sage: M = matrix(ZZ, [[3, 0],
            ....:                 [0, 0]])
            sage: data = ClassGroupData(M, primes=("P1","P2"))
            sage: data.class_in_quotient([3,0]) == 0
            True
        """
        hv = self.H(vector(ZZ, v))
        return self.qmap(hv)

    def is_zero(self, v):
        r"""
        Test whether a vector maps to zero in the quotient ``Q = H/R``.

        EXAMPLES::
            sage: M = matrix(ZZ, [[2, 0], [0, 2]])
            sage: data = ClassGroupData(M, primes=("P1","P2"))
            sage: data.is_zero([2,0])
            True
            sage: data.is_zero([1,0])
            False

        TESTS::
            sage: M = matrix(ZZ, [[0], [0]])
            sage: data = ClassGroupData(M, primes=("P1","P2"))
            sage: data.is_zero([0,0])
            True
            sage: data.is_zero([1,0])
            False
        """
        q = self.class_in_quotient(v)
        return (q == 0)

    # ---------------------------------------
    # Solve M*x = rhs over ZZ 
    # ---------------------------------------

    def solve_Mx_eq_rhs(self, rhs, *, require_solution=True, require_unique=True):
        r"""
        Solve the integral linear system ``M*x = rhs``.

        This method uses the Smith normal form `U M V = S` to decide solvability.


        INPUT:

        - ``rhs`` --  vector; an element of `\\ZZ^r`.

        - ``require_solution`` -- boolean (default: ``True``); if ``True`` then raise
        a ``ValueError`` when no integer solution exists. If ``False`` then return ``None`` in this case.

        - ``require_unique`` -- boolean (default: ``True``); if ``True`` then raise a
        ``ValueError`` when the set of integer solutions is infinite (equivalently,
        when `\\ker(M)\\neq 0`, i.e. when `n > \\mathrm{rank}(M)`), even if a solution exists. If ``False`` then return one solution in this case.

        OUTPUT:

        - a vector in `\\ZZ^n` giving a solution `x`, if a unique solution exists;
        - if no solution exists, then either raise a ``ValueError`` (if ``require_solution=True``) or return ``None`` (if ``require_solution=False``);
        - if infinitely many solutions exist, then either raise a ``ValueError`` (if ``require_unique=True``) or return one solution (if ``require_unique=False``).

        EXAMPLES::
            sage: # A full-rank square system (unique solution)
            sage: M = matrix(ZZ, [[3, 1],
            ....:                 [2, 1]])
            sage: data = ClassGroupData(M, primes=("P1","P2"))
            sage: rhs = vector(ZZ, [10, 7])
            sage: data.solve_Mx_eq_rhs(rhs)
            (3, 1)

        No solution::
            sage: M = matrix(ZZ, [[4, 0],
            ....:                 [0, 6]])
            sage: data = ClassGroupData(M, primes=("P1","P2"))
            sage: data.solve_Mx_eq_rhs([2, 3])   
            Traceback (most recent call last):
            ...
            ValueError: No integer solution to M*x = rhs.

        Infinitely many solutions::
            sage: M = matrix(ZZ, [[1, 0]])
            sage: data = ClassGroupData(M, primes=("P1",))
            sage: data.solve_Mx_eq_rhs([1])
            Traceback (most recent call last):
            ...
            ValueError: Infinitely many integer solutions (nontrivial kernel).

        A less trivial example::
            sage: M = matrix(ZZ, [[4, 6, 2],
            ....:                 [2, 8, 4],
            ....:                 [6, 2, 10]])
            sage: data = ClassGroupData(M, primes=("P1","P2","P3"))
            sage: x0 = vector(ZZ, [3, -2, 5])
            sage: rhs = M * x0
            sage: x = data.solve_Mx_eq_rhs(rhs)
            sage: M * x == rhs
            True
            sage: x == x0
            True

        """
        if  hasattr(rhs, "list"):
            rhs = rhs.list()
        b = vector(ZZ, rhs).column()
        # We already have U*M*V = S, so solve S*y = U*b, then x = V*y.
        b1 = self.U * b

        diag_len = min(self.S.nrows(), self.S.ncols())

        # Determine rank = number of nonzero diagonal entries of S
        rank = 0
        for i in range(diag_len):
            if ZZ(self.S[i, i]) != 0:
                rank += 1

        # Solvability checks in SNF coordinates:
        for i in range(diag_len):
            d = ZZ(self.S[i, i])
            if d == 0:
                if b1[i, 0] != 0:
                    if require_solution:
                        raise ValueError("No integer solution to M*x = rhs.")
                    return None
            else:
                if b1[i, 0] % d != 0:
                    if require_solution:
                        raise ValueError("No integer solution to M*x = rhs.")
                    return None

        for i in range(diag_len, self.S.nrows()):
            if b1[i, 0] != 0:
                if require_solution:
                    raise ValueError("No integer solution to M*x = rhs.")
                return None

        # Uniqueness check: if solvable and ker(M) != 0, there are infinitely many solutions
        if require_unique and self.n > rank:
            raise ValueError("Infinitely many integer solutions (nontrivial kernel).")

        # Build one solution y, then x = V*y
        y = vector(ZZ, [0] * self.S.ncols()).column()
        for i in range(diag_len):
            d = ZZ(self.S[i, i])
            if d != 0:
                y[i, 0] = b1[i, 0] // d

        xcol = self.V * y
     
        sol = vector(ZZ, [xcol[i, 0] for i in range(xcol.nrows())])

        assert (self.M * sol).column() == b

        return sol

# --------------------------------------------------------------------------------
# FLIRElement
# --------------------------------------------------------------------------------

class FLIRElement(CommutativeAlgebraElement):
    r"""
    An element of a FLIR given by an expression in the base fraction field.

    An element is represented by an element ``f`` of the base fraction field
    ``A.base.F``. By default we verify membership in the FLIR:

    .. MATH::

        f \\in A \\iff \text{for every chart } i,\\text{the substitution of } f
        \\text{ is a Laurent polynomial in chart } i.
    
    INPUT:

    - ``parent`` -- a :class:`FLIR`

    - ``f_in_base`` -- something coercible to ``parent._base_chart.F``

    - ``check`` -- boolean (default: ``True``); whether to verify membership

    EXAMPLES:

    A FLIR with a single chart (the base chart)::
        sage: K = QQ
        sage: base = FLIRChart.base(K, ("x1","x2"))
        sage: A = FLIR(base, (base,), compute_base_to_charts=False)
        sage: x1, x2 = base.F.gens()
        sage: f = A(x1/x2); f
        x1/x2
        sage: (f + 1) * (f - 1)
        (x1^2 - x2^2)/x2^2

    A FLIR with two charts.

    Define a second chart by ``y1 = (x2+1)/x1`` and ``y2 = x2``::
        sage: K = QQ
        sage: base = FLIRChart.base(K, ("x1","x2"))
        sage: x1, x2 = base.F.gens()
        sage: ch1 = FLIRChart(K, ("y1","y2"),
        ....:              base_fraction_field=base.F,
        ....:              this_to_base=[(x2 + 1)/x1, x2])
        sage: A = FLIR(base, (base, ch1), compute_base_to_charts=True)

    The element ``(x2+1)/x1`` lies in the FLIR (it becomes ``y1`` in the second chart)::
        sage: f = A((x2+1)/x1); f
        (x2 + 1)/x1
        sage: # Check substitution is Laurent in chart1:
        sage: ch1._substitute_from_base((x2+1)/x1)
        y1

    Membership checking rejects elements that are not Laurent after substitution.
    For example, ``1/(x1)`` becomes ``y1/(y2+1)`` in chart1, which is not Laurent::
        
        sage: A(1/(x1))
        Traceback (most recent call last):
        ...
        ValueError: Not in FLIR: not Laurent in chart ('y1', 'y2').
        ...
        

    One can bypass membership checking with ``check=False``::
        
        sage: A(1/(x1), check=False)
        1/x1

    TESTS::
        sage: K = QQ
        sage: base = FLIRChart.base(K, ("x1","x2"))
        sage: x1, x2 = base.F.gens()
        sage: ch1 = FLIRChart(K, ("y1","y2"),
        ....:              base_fraction_field=base.F,
        ....:              this_to_base=[(x2 + 1)/x1, x2])
        sage: A = FLIR(base, (base, ch1), compute_base_to_charts=True)
        sage: f = A((x2+1)/x1)
        sage: g = A((x2+1)/x1)
        sage: f == g
        True
        sage: (f + A(0)) == f
        True
        sage: (A(1) * f) == f
        True

    """

    
    def __init__(self, parent, f, check=True):
        self._f = f
        if check:
            parent._check_membership(f)
        CommutativeAlgebraElement.__init__(self, parent)

    def _repr_(self):
        return repr(self._f)
    
    def _richcmp_(self, other, op):
        # mathematical comparison
        if not isinstance(other, FLIRElement):
            return NotImplemented
        return richcmp((self.parent(), self._f), (other.parent(), other._f), op)


    @property
    def f(self):
        return self._f

    def _add_(self, other):
        return self.parent()(self._f + other._f, check=False)

    def _sub_(self, other):
        return self.parent()(self._f - other._f, check=False)

    def _mul_(self, other):
        return self.parent()(self._f * other._f, check=False)

    def _neg_(self):
        return self.parent()(-self._f, check=False)

    def __pow__(self, n):
        return self.parent()(self._f ** int(n), check=False)
    
    def __truediv__(self, other):
        """Division in the ambient fraction field, with membership check."""
        A = self.parent()
        if not isinstance(other, FLIRElement) or other.parent() is not A:
            other = A(other)  # try coercion into the FLIR

        if other._f == 0:
            raise ZeroDivisionError("Division by zero.")

        q = self._f / other._f  # computed in the ambient fraction field

        A._check_membership(q)

        return A(q, check=False)

    def divisor(self) -> FLIRDivisor:
        """
        Compute div_A(f) as an element of Div(A).
        """
        if self.f == 0:
            raise ValueError("Divisor of 0 is undefined.")

        A = self.parent()
        DivA = A.Div()

        # P_i-part via valuations 
        Pis = A.extra_primes()
        data: Dict[FLIRPrimeDivisor, int] = {}
        for Pj in Pis:
            e = A._valuation_of_base_element_at_prime(self.f, Pj)
            if e != 0:
                data[Pj] = int(e)

        # base primes from factorization in base Laurent ring 
        # We compute numerator/denominator in base chart, clear monomial, factor polynomial part in base.P
        poly = A._base_chart._laurent_poly_to_poly_up_to_unit(self.f)[0]


        if poly == 0:
            raise ValueError("Divisor of 0 is undefined.")
        for q, e in poly.factor():
            if q.is_unit():
                continue
            q0 = A._normalize(q)

            if q0 in list(A._base_chart.P.gens()):
                continue
            Pbase = FLIRPrimeDivisor(A._base_chart, q0)

            data[Pbase] = data.get(Pbase, 0) + int(e)
            if data[Pbase] == 0:
                del data[Pbase]

        return DivA(data)

    def factor(self, verbose=False):
        return FLIRFactorization.from_element(self, verbose=verbose)
    
    def atoms(self):
        return self.factor().atoms()



# -----------------------------------------------------------------------------
# FLIR
# -----------------------------------------------------------------------------


class FLIR(Parent, UniqueRepresentation):
    r"""
    A finite Laurent intersection ring (FLIR).

    Let `K` be a field and let

    .. MATH::

        A = \\bigcap_{i=0}^m K[y_{i,1}^{\\pm 1}, \\dots, y_{i,n}^{\\pm 1}]
        \\subseteq K(x_1,\\dots,x_n),

    where chart `0` is a chosen base Laurent polynomial ring

    .. MATH::

        K[x_1^{\\pm 1},\\dots,x_n^{\\pm 1}],

    and each additional chart is given by an explicit birational substitution from
    its Laurent coordinates into the base fraction field.

    This class represents the ring `A` as a subring of the base fraction field.
    An element of `A` is stored by its expression in the base fraction field, and
    membership is checked by verifying that it becomes Laurent in every chart.

    In addition to ring operations, the class provides:

    - the distinguished height-1 prime divisors containing `x_1 \\cdots x_n`,
    - the divisor group `\\mathrm{Div}(A)`,
    - the class group `\\mathrm{Cl}(A)`,
    - divisor-class computations,
    - principality tests and principal generators,
    - atom search and factorization by divisor methods.

    INPUT:

    - ``base_chart`` -- a :class:`FLIRChart`; the chosen base Laurent chart

    - ``charts`` -- iterable of :class:`FLIRChart`; the charts defining the intersection.
      Usually this contains the base chart together with finitely many birational charts.

    - ``compute_base_to_charts`` -- boolean (default: ``True``); if ``True``, compute
      the inverse substitutions ``base_to_this`` for charts where only
      ``this_to_base`` is given.

    - ``category`` -- optional Sage category. By default this is taken to be a
      commutative algebra over the coefficient field.

    EXAMPLES::

    The simplest FLIR is just a Laurent polynomial ring itself, viewed as the
    intersection of one chart::
        sage: K = QQ
        sage: C0 = FLIRChart.base(K, ("x1","x2"))
        sage: A = FLIR(C0, (C0,), compute_base_to_charts=False)
        sage: A
        FLIR over Rational Field
          rank n = 2
          #charts = 1
          base vars = ('x1', 'x2')

    Its fraction field is the fraction field of the base Laurent ring::
        sage: A.fraction_field is C0.F
        True

    Elements are represented by expressions in the base fraction field::
        sage: x1, x2 = C0.F.gens()
        sage: f = A((x1 + 1)/x2)
        sage: f
        (x1 + 1)/x2

    Arithmetic is performed in the ambient fraction field::
        sage: g = A(x1)
        sage: f + g
        (x1*x2 + x1 + 1)/x2
        sage: f * g
        (x1^2 + x1)/x2
        sage: -f
        (-x1 - 1)/x2

    Since there is only one chart, every Laurent expression in the base chart is in `A`::
        sage: A(x1^-3 * x2 + x2^-1)
        (x1^3 + x2^2)/(x1^3*x2)

    We now build a nontrivial example with a second chart. Let

    .. MATH::

        y_1 = \\frac{x_2+1}{x_1}, \\qquad y_2 = x_2.

    Then the corresponding FLIR consists of rational functions that are Laurent in
    both the `x`-chart and the `y`-chart::
        sage: base = FLIRChart.base(QQ, ("x1","x2"))
        sage: x1, x2 = base.F.gens()
        sage: ch1 = FLIRChart(QQ, ("y1","y2"),
        ....:                 base_fraction_field=base.F,
        ....:                 this_to_base=[(x2 + 1)/x1, x2])
        sage: A = FLIR(base, (base, ch1), compute_base_to_charts=True)

    The chart inverse is computed automatically::
        sage: ch1.base_to_this
        [(y2 + 1)/y1, y2]

    The element `(x_2+1)/x_1` belongs to `A`, because it becomes the Laurent monomial
    `y_1` in the second chart::
        sage: A((x2 + 1)/x1)
        (x2 + 1)/x1

    But `1/x_1` does not belong to `A`, since after substitution it becomes
    `y_1/(y_2+1)`, which is not Laurent::
        sage: A(1/x1)
        Traceback (most recent call last):
        ...
        ValueError: Not in FLIR: not Laurent in chart ('y1', 'y2').
        ...

    Membership checks can be bypassed when needed for internal constructions::
        sage: A(1/x1, check=False)
        1/x1

    The distinguished extra primes are the height-1 primes containing `x_1 \\cdots x_n`
    that are not already visible as coordinate primes in the base chart::
        sage: A.extra_primes()
        [PrimeDivisor(chart=('y1', 'y2'), p=y2 + 1)]

    The divisor group is implemented as a free abelian group on prime divisors::
        sage: DivA = A.Div()
        sage: DivA
        Divisor group Div(A) of FLIR over Rational Field
          rank n = 2
          #charts = 2
          base vars = ('x1', 'x2') with coefficients in Integer Ring

    An element determines a divisor in `\\mathrm{Div}(A)`::
        sage: f = A((x2 + 1)/x1)
        sage: D = f.divisor()
        sage: D
        1*PrimeDivisor(chart=('x1', 'x2'), p=x2 + 1)

    Since `f` is itself an element of the ring, its divisor is principal and therefore
    trivial in the class group::
        sage: A.divisor_class(D) == 0
        True
        sage: D.is_principal()
        True

    A principal generator can be reconstructed from a principal divisor::
        sage: A.principal_generator(D)
        (x2 + 1)/x1

    A larger example is given by the helper function ``example_A3()``::
        sage: A3 = example_A3()
        sage: A3.n
        3
        sage: len(A3.charts)
        5
        sage: len(A3.extra_primes()) == 4
        True

    In that example one can compute divisors and atoms::
        sage: x1, x2, x3 = A3._base_chart.L.gens()
        sage: a = A3(x1**2 * (x1 + x3) / x2)
        sage: D = a.divisor()
        sage: D
        2*PrimeDivisor(chart=('x_21', 'x_22', 'x_23'), p=x_22 + 1) +
        3*PrimeDivisor(chart=('x_41', 'x_42', 'x_43'), p=x_42 + 1) +
        1*PrimeDivisor(chart=('x1', 'x2', 'x3'), p=x1 + x3)
        sage: D.is_effective()
        True
        sage: a.atoms()
        [x1, (x1 + x3)/x2]

        
    Factorizations are returned as a list of atoms and a list of exponent patterns::
        sage: a.factor()
        [[((x1 + x3)/x2, 1), (x1, 2)]]
    
    An example of non-half-factorial element::
        sage: x1, x2, x3 = A3._base_chart.L.gens()
        sage: P1, P2, P3, P4 = A3.extra_primes()
        sage: P2 = A3.divisor({P2: 1})
        sage: P4 = A3.divisor({P4: 1})
        sage: g = A3((x1*x2*x3)**2 + ((x1-x3)**2)*(x3+1))
        sage: f = A3((x1*x2*x3)**2 + (((x2+1)**2)/(x1*x3))*(x3+1))
        sage: Q1 = g.divisor()
        sage: Q2 = f.divisor()
        sage: div = P2 + P2 + Q1 + Q2 + P4 + P4 
        sage: h = A3.principal_generator(div)
        sage: h.factor().set_of_lengths()
        {4, 5}


    TESTS:

    Construction is unique-representation compatible::
        sage: C0 = FLIRChart.base(QQ, ("x","y"))
        sage: A1 = FLIR(C0, (C0,), compute_base_to_charts=False)
        sage: A2 = FLIR(C0, (C0,), compute_base_to_charts=False)
        sage: A1 is A2
        True

    Coercion of scalars from the base field works::
        sage: A = FLIR(C0, (C0,), compute_base_to_charts=False)
        sage: A(1)
        1
        sage: A(QQ(3)/2)
        3/2

    The zero and unit elements are always accepted without membership checks::
        sage: A(0)
        0
        sage: A(1)
        1

    Incomplete charts are rejected::
        sage: Cbad = FLIRChart(QQ, ("u","v"), base_fraction_field=C0.F)
        sage: FLIR(C0, (C0, Cbad), compute_base_to_charts=True)
        Traceback (most recent call last):
        ...
        ValueError: Chart ('u', 'v') incomplete: need this_to_base at least.

    Base-field mismatch in charts is caught when inversion is requested::
        sage: C0a = FLIRChart.base(QQ, ("x1","x2"))
        sage: C0b = FLIRChart.base(QQ, ("u1","u2"))
        sage: x1, x2 = C0a.F.gens()
        sage: C1 = FLIRChart(QQ, ("y1","y2"),
        ....:                base_fraction_field=C0a.F,
        ....:                this_to_base=[x1, x2/x1])
        sage: C1.compute_base_to_this(C0b)
        Traceback (most recent call last):
        ...
        ValueError: Base fraction field mismatch.
    """
    Element = FLIRElement

    @staticmethod
    def __classcall__(cls, base_chart, charts, compute_base_to_charts=True, category=None):
        # normalize inputs so caching works (UniqueRepresentation)
        charts = tuple(charts)
        compute_base_to_charts = bool(compute_base_to_charts)
        return super().__classcall__(
            cls,
            base_chart,
            charts,
            compute_base_to_charts=compute_base_to_charts,
            category=category,
        )

    def _init_flir_structure(self, base_chart: FLIRChart, charts, compute_base_to_charts=True):
        self._base_chart = base_chart
        self.charts = list(charts)
        self.n = self._base_chart.n

        self._extra_primes_cache = None
        self._class_data_cache = None
        self._div_group_cache = None

        if compute_base_to_charts:
            for ch in self.charts:
                if ch is self._base_chart:
                    continue
                if ch.base_fraction_field is None:
                    ch.base_fraction_field = self._base_chart.F
                if ch.base_to_this is None and ch.this_to_base is not None:
                    ch.compute_base_to_this(self._base_chart)
                if ch.base_to_this is None or ch.this_to_base is None:
                    raise ValueError(f"Chart {ch.var_names} incomplete: need this_to_base at least.")
                
    def __init__(self, base_chart: FLIRChart, charts, compute_base_to_charts=True, category=None):
        self._init_flir_structure(base_chart, charts, compute_base_to_charts=compute_base_to_charts)

        if category is None:
            category = CommutativeAlgebras(base_chart.base_ring.category())
        
        Parent.__init__(self, base=self._base_chart.base_ring, category=category)
        self.element_class = self.Element # type: ignore[assignment]

    @property
    def fraction_field(self):
        r"""
        The fraction field of the FLIR is the fraction field of the base Laurent ring.
        EXAMPLES::
            sage: C = FLIRChart.base(QQ, ("x","y"))
            sage: A = FLIR(C, (C,), compute_base_to_charts=False)
            sage: A.fraction_field is C.F
            True
        """
        return self._base_chart.F

    def _repr_(self) -> str:
        lines = []
        lines.append(f"FLIR over {self._base_chart.base_ring}")
        lines.append(f"  rank n = {self.n}")
        lines.append(f"  #charts = {len(self.charts)}")
        lines.append(f"  base vars = {self._base_chart.var_names}")
        return "\n".join(lines)

    __repr__ = _repr_

    def _coerce_map_from_(self, S):
        # allow coercion from the base ring 
        if self.base_ring().has_coerce_map_from(S):
            return True
        return super()._coerce_map_from_(S)
    
    def _element_constructor_(self, x, check=True):
        r"""
        We allow coercion from the base ring, and we also allow direct 
        construction from the fraction field of the base chart (which is the fraction field of the FLIR). 
        In either case, we check membership by verifying that the element becomes Laurent in every chart.

        TESTS::
            sage: C = FLIRChart.base(QQ, ("x","y"))
            sage: A = FLIR(C, (C,), compute_base_to_charts=False)
            sage: x, y = C.F.gens()
            sage: A(x/y)
            x/y
            sage: f = A(x)
            sage: A(f) is f
            True
        """
        if isinstance(x, FLIRElement) and x.parent() is self:
            return x
        # make sure scalars coerce into the base fraction field
        x = self._base_chart.F(x)
        # avoid membership checks for 0,1 (used during Parent init)
        if x == 0 or x == 1:
            check = False
        return self.element_class(self, x, check=check)

    
    def _check_membership(self, f):
        r"""
        Check that the given expression in the base fraction field 
        is actually an element of the FLIR by verifying that it 
        becomes Laurent in every chart.

        TESTS::
            sage: base = FLIRChart.base(QQ, ("x1","x2"))
            sage: x1, x2 = base.F.gens()
            sage: ch1 = FLIRChart(QQ, ("y1","y2"),
            ....:                 base_fraction_field=base.F,
            ....:                 this_to_base=[(x2 + 1)/x1, x2])
            sage: A = FLIR(base, (base, ch1), compute_base_to_charts=True)
            sage: A._check_membership((x2 + 1)/x1)
            sage: A._check_membership(1/x1)
            Traceback (most recent call last):
            ...
            ValueError: Not in FLIR: not Laurent in chart ('y1', 'y2').
            ...
        """
        for ch in self.charts:
            expr = ch._substitute_from_base(f)
            try:
                _ = ch.L(expr)
            except Exception as e:
                raise ValueError(
                    f"Not in FLIR: not Laurent in chart {ch.var_names}.\nSubstituted: {expr}"
                ) from e
    
    def _normalize(self, p):
        r"""
        Normalize a nonzero polynomial by dividing by its leading coefficient.
        The leading term depends on the term order of the parent polynomial ring.
        """
        if p == 0:
            return p
        return p / p.leading_coefficient()

    def _base_gens(self):
        return self._base_chart.F.gens()

    def _xprod_in_base(self):
        xs = self._base_gens()
        out = xs[0]
        for i in range(1, self.n):
            out *= xs[i]
        return out

    def _prod_of_chart_gens_as_base_expr(self, chart: FLIRChart):
        if chart.this_to_base is None:
            raise ValueError("Chart missing this_to_base substitution.")
        out = chart.this_to_base[0]
        for i in range(1, chart.n):
            out *= chart.this_to_base[i]
        return out

    def _gcd_polys(self, P, polys):
        if not polys:
            return P(0)
        g = polys[0]
        for h in polys[1:]:
            g = g.gcd(h)
            if g == 1:
                return g
        return g

    # -------------------------------------------------------
    # extra primes: primes over x1*...*xn, i.e. primes missing in FLIR._base_chart
    # -------------------------------------------------------

    def _extra_primes_in_chart(self, i: int) -> List[FLIRPrimeDivisor]:
        """
        For chart i:
          f_i = image of x1*...*xn in chart i
          g_{j,i} = image of (prod of chart j variables) in chart i for j<i
        Return irreducible factors of gcd(f_i, g_{0,i},...,g_{i-1,i}) (after clearing Laurent units).
        """
        chart_i = self.charts[i]

        xprod_base = self._xprod_in_base()
        f_i = chart_i._substitute_from_base(xprod_base)
        fP, _ = chart_i._laurent_poly_to_poly_up_to_unit(f_i)

        gcd_list = [fP]
        for j in range(0, i):
            chart_j = self.charts[j]
            yjprod_base = self._prod_of_chart_gens_as_base_expr(chart_j)
            gij = chart_i._substitute_from_base(yjprod_base)
            gP, _ = chart_i._laurent_poly_to_poly_up_to_unit(gij)
            gcd_list.append(gP)

        G = self._gcd_polys(chart_i.P, gcd_list)
        if G == 0 or G == 1:
            return []

        primes = []
        for f, _e in G.factor():
            if f.is_unit():
                continue
            primes.append(FLIRPrimeDivisor(chart_i, self._normalize(f)))
        return primes

    def extra_primes(self, recompute: bool = False) -> List[FLIRPrimeDivisor]:
        r"""
        Return the list P1,...,Pr of prime divisors containing x1*...*xn.

        EXAMPLES::
            sage: base = FLIRChart.base(QQ, ("x1","x2"))
            sage: x1, x2 = base.F.gens()
            sage: ch1 = FLIRChart(QQ, ("y1","y2"),
            ....:                 base_fraction_field=base.F,
            ....:                 this_to_base=[(x2 + 1)/x1, x2])
            sage: A = FLIR(base, (base, ch1), compute_base_to_charts=True)
            sage: A.extra_primes()
            [PrimeDivisor(chart=('y1', 'y2'), p=y2 + 1)]

        TESTS::
            sage: A.extra_primes() is A.extra_primes()
            True
            sage: A.extra_primes(recompute=True) == A.extra_primes()
            True
        """
        if (not recompute) and self._extra_primes_cache is not None:
            return self._extra_primes_cache

        primes = []
        for i in range(len(self.charts)):
            primes.extend(self._extra_primes_in_chart(i))

        # deduplicate 
        primes = list(dict.fromkeys(primes))

        self._extra_primes_cache = primes
        return primes

    # ------------------------
    # valuation machinery
    # ------------------------

    def _valuation_of_base_element_at_prime(self, base_expr, P: FLIRPrimeDivisor) -> int:
        """
        v_P(base_expr): multiplicity of P.irreducible in the factorization after substituting into chart.
        """
        chart = P.chart
        expr_in_chart = chart._substitute_from_base(base_expr)
        gP, _ = chart._laurent_poly_to_poly_up_to_unit(expr_in_chart)

        e = 0
        for f, k in gP.factor():
            if f.is_unit():
                continue
            if self._normalize(f) == P.irreducible:
                e += int(k)
        return e

    def _coeffs_on_extra_primes(self, base_expr) -> Dict[FLIRPrimeDivisor, int]:
        """
        Return dict Pi -> v_{Pi}(base_expr) for Pi in extra_primes.
        """
        coeffs = {}
        for Pi in self.extra_primes():
            ei = int(self._valuation_of_base_element_at_prime(base_expr, Pi))
            if ei != 0:
                coeffs[Pi] = ei
        return coeffs

    # ------------------------
    # class group
    # ------------------------

    def divisor_group(self):
        """Return the divisor group associated to this object.

        Constructs (on first call) and returns an FLIRDivisorGroup representing the
        divisor group of this FLIR instance. The constructed group is created with
        base_ring=ZZ and cached on the instance as `_div_group_cache` so subsequent
        calls return the same object.

        The cache attribute name is `_div_group_cache`.
        """

        if getattr(self, "_div_group_cache", None) is None:
            self._div_group_cache = FLIRDivisorGroup(self, base_ring=ZZ)
        return self._div_group_cache
    
    def Div(self):
        """
        Return the divisor group Div(A) of this FLIR.

        This is the free abelian group on all height-1 prime divisors of A.
        Elements are instances of :class:`FLIRDivisor`. The group instance is
        created once and cached by :meth:`divisor_group`; use ``A.Div()`` as
        the canonical accessor for the ambient divisor group.
        """
        return self.divisor_group()

    def divisor(self, data=None):
        """
        Create a divisor in Div(A) from `data`.

        Accepted forms:
         - None or 0 -> the zero divisor
         - FLIRDivisor (in the same Div(A)) -> returned unchanged
         - FLIRPrimeDivisor -> the prime with coefficient 1
         - dict {prime: exponent} or iterable of (prime, exponent) pairs -> passed to Div(A)
         - FLIRElement or anything coercible to the base fraction field -> principal divisor of that element

        Examples::
            sage: A = example_A3()
            sage: C1, C2, C3, C4, C5 = A.charts
            sage: x1, x2, x3 = C1.P.gens()
            sage: y1, y2, y3 = C3.P.gens()
            sage: A.divisor() == A.Div().zero()
            True
            sage: P = FLIRPrimeDivisor(A._base_chart, x1 + x3)
            sage: A.divisor(P)
            1*PrimeDivisor(chart=('x1', 'x2', 'x3'), p=x1 + x3)
            
            sage: Q = FLIRPrimeDivisor(C3, y1*y3 + y2)
            sage: A.divisor({P: 2, Q: 3})  
            2*PrimeDivisor(chart=('x1', 'x2', 'x3'), p=x1 + x3) +
             3*PrimeDivisor(chart=('x_21', 'x_22', 'x_23'), p=x_21*x_23 + x_22)

            sage: x = A._base_chart.F.gens()[0]
            sage: A.divisor(x) 
            1*PrimeDivisor(chart=('x_21', 'x_22', 'x_23'), p=x_22 + 1) +
             1*PrimeDivisor(chart=('x_41', 'x_42', 'x_43'), p=x_42 + 1) 
        """
        DivA = self.Div()

        # zero
        if data is None or data == 0:
            return DivA.zero()

        # already a divisor
        if isinstance(data, FLIRDivisor):
            if data.parent() is DivA:
                return data
            raise ValueError("Divisor belongs to a different ambient divisor group.")

        # single prime -> 1*P
        if isinstance(data, FLIRPrimeDivisor):
            return DivA({data: 1})

        # dict or iterable-of-pairs -> let FLIRDivisorGroup handle construction/validation
        if isinstance(data, dict):
            return DivA(data)
        if isinstance(data, (list, tuple)):
            # could be list of pairs or a sequence representing a divisor
            try:
                return DivA(data)
            except Exception:
                # fall through to attempt element coercion
                pass

        # FLIRElement -> use its divisor
        if isinstance(data, FLIRElement):
            return data.divisor()

        # finally: try to coerce to a base-field element and form principal divisor (no membership check)
        try:
            f = self._base_chart.F(data)
        except Exception as ex:
            raise ValueError("Cannot construct a divisor from the provided data.") from ex

        elt = self.element_class(self, f, check=False)
        return elt.divisor()
        return self.divisor_group()(data)

    def _compute_relation_matrix(self):
        """
        Build M = C^T, where C_{i,j} = v_{P_j}(x_i)
        """
        primes = self.extra_primes()
        r = len(primes)
        xs = self._base_gens()

        C = matrix(ZZ, self.n, r)
        for i in range(self.n):
            for j, Pj in enumerate(primes):
                C[i,j] = ZZ(self._valuation_of_base_element_at_prime(xs[i], Pj))
        return C.transpose(), primes
    
    def class_data(self, recompute=False):
        """ Return the ClassGroupData object (cached)"""

        if (self._class_data_cache is None) or recompute:
            M, primes = self._compute_relation_matrix()
            self._class_data_cache = ClassGroupData(M, primes)
        return self._class_data_cache

    def class_group(self):
        """
        Return Cl(A) as the quotient module Q = Z^r /im(M).

        EXAMPLES::
            sage: # trivial FLIR (single base chart) has trivial class group
            sage: C = FLIRChart.base(QQ, ("x1","x2"))
            sage: A = FLIR(C, (C,), compute_base_to_charts=False)
            sage: A.class_group().order()
            1
            sage: A3 = example_A3()
            sage: A3.class_group()
            Multiplicative Abelian group isomorphic to Z
            sage: A3.class_group().invariants()  
            (0,)

        """
        return self.class_data().Cl()
    

    # ------------------------
    # class of a divisor
    # ------------------------


    def _to_H(self, D: FLIRDivisor):
        """
        Project ``D \\in Div(A)`` to ``hv \\in H = \\ZZ^r`` (FreeModule element)
        with respect to the ordering of extra_primes().
        """

        if isinstance(D, FLIRPrimeDivisor):
            D = self.divisor({D: 1})
        elif isinstance(D, dict):
            D = self.divisor(D)

        data = self.class_data()
        Pis = data.primes
        coords = [ZZ(0)] * len(Pis)

        for P, e in D:
            e = ZZ(e)
            if P in Pis:
                coords[Pis.index(P)] += e
                continue
            if P.chart == self._base_chart:
                p = self._base_chart.F(P.irreducible)
                for j, Pj in enumerate(Pis):
                    coords[j] -= e * ZZ(self._valuation_of_base_element_at_prime(p, Pj))
                continue
            raise ValueError(
                "Prime in divisor is neither a distinguished P_i nor a base prime (base, r)"
            )
        return data.H(coords)

    def divisor_class(self, D):
        r"""Compute the class of a divisor in the divisor class group.

            Accept:
            - FLIRElement (the divisor of the element is used)
            - FLIRDivisor
            - FLIRPrimeDivisor (interpreted as 1*P)
            - dict {prime: exponent}

            Return: class [D] in Cl(A) = Z^r / im(M).

            EXAMPLES::

                sage: A3 = example_A3()
                sage: z1, z2, z3 = A3._base_gens()
                sage: z = A3(z1*z3)
                sage: D = A3.divisor(z)
                sage: A3.class_group()
                Multiplicative Abelian group isomorphic to Z
                sage: A3.divisor_class(D)
                (0)

            The class of the divisor of an element is trivial, since it is principal.
            ``divisor_class`` accepts a :class:`FLIRElement` directly, giving the
            same result as passing its divisor::

                sage: A3.divisor_class(z) == A3.divisor_class(D)
                True
            """
        if isinstance(D, FLIRElement):
            D = self.divisor(D)
        elif isinstance(D, FLIRPrimeDivisor):
            D = self.divisor({D: 1})
        elif isinstance(D, dict):
            D = self.divisor(D)

        hv = self._to_H(D)
        return self.class_data().class_in_quotient(hv)

    def is_principal_divisor(self, D: FLIRDivisor) -> bool:
        """
        Check if the divisor D is principal, i.e. if [D] = 0 in Cl(A).

        EXAMPLES:

        The divisor of an element of the FLIR is always principal::

            sage: A3 = example_A3()
            sage: z1, z2, z3 = A3._base_gens()
            sage: z = A3(z1*z3)
            sage: D = A3.divisor(z)
            sage: A3.is_principal_divisor(D)
            True

        A single prime appearing in ``D``, taken on its own, need not be
        principal::

            sage: P = list(D.support())[0]
            sage: A3.is_principal_divisor(A3.divisor({P: 1}))
            False
        """
        return self.class_data().class_in_quotient(self._to_H(D)) == 0
    
    def principal_generator(self, D: FLIRDivisor):
        """
        If D is principal, attempt to produce ``f \\in K(x)`` such that ``div(f) = D``.

        f = t * x^m
        where t is the product of base primes (base,r) appearing in D,
        and m solves M*m = rhs over ZZ with rhs = D - div(t).

        Returns ``None`` if ``D`` is not principal.

        EXAMPLES:

        Recovering a generator for the divisor of a known element::

            sage: A3 = example_A3()
            sage: z1, z2, z3 = A3._base_gens()
            sage: z = A3(z1*z3)
            sage: D = A3.divisor(z)
            sage: g = A3.principal_generator(D)
            sage: g
            x1*x3
            sage: A3.divisor(g) == D
            True

        If the divisor is not principal, ``None`` is returned::

            sage: P = list(D.support())[0]
            sage: A3.principal_generator(A3.divisor({P: 1})) is None
            True
        """
        if not self.is_principal_divisor(D):
            return None
        
        data = self.class_data()
        Pis = data.primes

        t = self._base_chart.F(1)
        for P, e in D:
            if (P.chart is self._base_chart):
                t *= self._base_chart.F(P.irreducible) ** ZZ(e)
        Dt = self(t, check=False).divisor() # div(t) in Div(A)
        R = D - Dt

        rhs = [ZZ(R.coeff(Pi)) for Pi in Pis]

        m = data.solve_Mx_eq_rhs(rhs, require_unique=False)
        if m is None:
            return None
        
        xs = self._base_gens()
        mon = self._base_chart.F(1)
        for i in range(self.n):
            mon *= xs[i] ** m[i]

        gen = t * mon

        return self(gen)

    # ----------------------------------------------------------------
    # helpers factorization
    # ----------------------------------------------------------------

    def _divisor_basis(self, D: FLIRDivisor):
        """
        Return a stable ordered list of prime divisors occuring in the support of D
        """
        Ps = sorted(D.support(), key=lambda P: repr(P))
        return Ps
    
    def _divisor_to_vector(self, D: FLIRDivisor, basis):
        """
        Given a divisor D and a list of prime divisors as basis, return the vector of coefficients.
        """
        return vector(ZZ, [D.coeff(P) for P in basis])
    
    def _vector_to_divisor(self, v, basis):
        """
        Given a vector of coefficients and a list of prime divisors as basis, return the corresponding divisor.
        """
        data = {P: int(e) for P, e in zip(basis, v) if e != 0}
        return self.divisor(data)
    
    def _is_effective_on_basis(self, v):
        return all(int(c) >= 0 for c in v)
    
    def _componentwise_leq(self, a, b):
        return all(int(x) <= int(y) for x, y in zip(a, b))
    

    def _principal_support_lattice(self, D):
        """
        Return (support, K), where support = [Q1,...,Qt] is the support of D
        in a fixed order, and K is the lattice of c in Z^t such that
        sum c_i Q_i is principal.

        This computes K = { c : B*c in Im_Z(M) } using the Smith normal form of M.
        """
        support = sorted(D.support(), key=lambda Q: repr(Q))
        data = self.class_data()
        M = matrix(ZZ, data.M)
        r = M.nrows()
        t = len(support)

        # B = [ h(Q1) ... h(Qt) ]  as an r x t integer matrix
        cols = []
        for Q in support:
            h = self._to_H(self.divisor({Q: 1}))
            cols.append(vector(ZZ, h))
        B = matrix(ZZ, r, t, lambda i, j: ZZ(cols[j][i]))

        # Smith normal form: U*M*V = S
        S, U = self.class_data().S, self.class_data().U

        # rank = number of nonzero diagonal entries
        diag = [ZZ(S[i, i]) for i in range(min(S.nrows(), S.ncols()))]
        rk = sum(1 for d in diag if d != 0)

        UB = U * B

        Aeq_rows = []

        # divisibility rows
        for i in range(rk):
            row = [ZZ(UB[i, j]) for j in range(t)] + [ZZ(0)] * rk
            row[t + i] = -ZZ(diag[i])
            Aeq_rows.append(row)

        # zero rows
        for i in range(rk, r):
            row = [ZZ(UB[i, j]) for j in range(t)] + [ZZ(0)] * rk
            Aeq_rows.append(row)

        if Aeq_rows:
            Aeq = matrix(ZZ, len(Aeq_rows), t + rk, Aeq_rows)
            Kbig = Aeq.right_kernel()
            proj = [vector(ZZ, g[:t]) for g in Kbig.basis()]
        else:
            # Degenerate case: no constraints, so K = Z^t
            proj = [vector(ZZ, [1 if i == j else 0 for i in range(t)]) for j in range(t)]

        Zt = ZZ**t  # type: ignore[operator]
        K = Zt.submodule(proj)

        return support, K


    def _principal_effective_subdivisors_via_kernel(self, D):
        """
        Enumerate all effective principal subdivisors E <= D
        using the SNF-based lattice and polyhedral enumeration.
        """
        support, K = self._principal_support_lattice(D)
        bounds = vector(ZZ, [ZZ(D.coeff(Q)) for Q in support])

        basis = [vector(ZZ, b) for b in K.basis()]
        t = len(support)
        r = len(basis)

        if r == 0:
            return [self.divisor()]

        # A has the basis vectors of K as columns, so c = A*m
        A = matrix(ZZ, t, r, lambda i, j: ZZ(basis[j][i]))

        # Polyhedron: 0 <= A*m <= bounds
        # Sage inequalities are b + a_1 x_1 + ... + a_r x_r >= 0
        ieqs = []

        for i in range(t):
            # (A*m)_i >= 0
            ieqs.append([0] + [ZZ(A[i, j]) for j in range(r)])

        for i in range(t):
            # bounds[i] - (A*m)_i >= 0
            ieqs.append([ZZ(bounds[i])] + [-ZZ(A[i, j]) for j in range(r)])

        P = Polyhedron(ieqs=ieqs, base_ring=QQ)

        out = []
        seen = set()

        for m in P.integral_points():
            m = vector(ZZ, m)
            c = vector(ZZ, A * m)

            # safety check
            if not all(0 <= c[i] <= bounds[i] for i in range(t)):
                continue

            tup = tuple(int(c[i]) for i in range(t))
            if tup in seen:
                continue
            seen.add(tup)

            E = self.divisor({Q: tup[i] for i, Q in enumerate(support) if tup[i] != 0})
            out.append(E)

        return out
    
    
    def _find_atoms(self, a: FLIRElement):
        if a.f == 0:
            raise ValueError("Zero has no factorizations.")

        Da = a.divisor()
        basis = self._divisor_basis(Da)

        principal_subdivisors = self._principal_effective_subdivisors_via_kernel(Da)
        vectors_of_principal_subdivisors = [
            self._divisor_to_vector(E, basis) for E in principal_subdivisors
        ]

        def minimal_vectors(vectors):
            vecs = [tuple(v) for v in vectors if any(v)]  # remove zero

            minimals = []
            for u in vecs:
                is_minimal = True
                for v in vecs:
                    if v == u:
                        continue
                    if all(vi <= ui for vi, ui in zip(v, u)) and any(vi < ui for vi, ui in zip(v, u)):
                        is_minimal = False
                        break
                if is_minimal:
                    minimals.append(u)

            return minimals

        atoms_valuations = minimal_vectors(vectors_of_principal_subdivisors)

        atoms = []
        seen = set()
        for kv in atoms_valuations:
            E = self._vector_to_divisor(kv, basis)
            gen = self.principal_generator(E)
            elt = self(gen)
            key = elt.f
            if key not in seen:
                seen.add(key)
                atoms.append(elt)

        return atoms
    
    def _all_solutions_bounded (self, target,  cols, bounds):
        """
        Solve target = sum ci * cols[i] with ci in [0, bounds[i]], ci\\in ZZ_{\\ge 0}.
        """
        n = len(cols)
        m = len(target)
        sols = []

        def feasible(rem, i):
            if any(x < 0 for x in rem):
                return False
            for j in range(m):
                max_cover = 0
                for k in range(i,n):
                    max_cover += bounds[k] * cols[k][j]
                if rem[j] > max_cover:
                    return False
            return True
        
        def rec(i, rem, tvec):
            if i == n:
                if all(x == 0 for x in rem):
                    sols.append(tuple(tvec))
                return
            if not feasible(rem, i):
                return
            col = cols[i]
            ub = bounds[i]
            for j in range(m):
                cj = col[j]
                if cj > 0:
                    ub = min(ub, rem[j] // cj)
            for t in range(ub + 1):
                rec(i + 1, rem - t * col, tvec + [t])
        
        rec(0, target, [])
        return sols
    
    def factorizations(self, a: FLIRElement):
        """
        Return  all factorizations of a as follows:
            - list of atoms that divide a
            - list of factorizations, where each factorization is a list of (atom, exponent) pairs.
            - units 
        """

        if a.f == 0:
            raise ValueError("Zero has no factorizations.")
        Da = a.divisor()

        basis = self._divisor_basis(Da)
        va = self._divisor_to_vector(Da, basis)
        if not self._is_effective_on_basis(va):
            raise ValueError(f"Negative coefficients in divisor. The element is not an element of {self}.")
        
        atoms = self._find_atoms(a)
        if len(atoms) == 0:
            return atoms, [[]], [a]
        
        
        atom_divs = [self._divisor_to_vector(b.divisor(), basis) for b in atoms]
        bounds =[]
        for col in atom_divs:
            ub = None
            for j in range(len(va)):
                cj = col[j]
                if cj > 0:
                    q = va[j] // cj
                    ub = q if ub is None else min(ub, q)
            if ub is None:
                ub = 0
            bounds.append(ub)
        sol_tuples = self._all_solutions_bounded(va, atom_divs, bounds)
        factorizations = []
        units = []

        for t in sol_tuples:
            prod = self(1, check=False)
            fac = []
            for b, e in zip(atoms, t):
                if e:
                    prod *= b ** e
                    fac.append((b, e))
            eps = a / prod
            factorizations.append(fac)
            units.append(eps)
        
        return atoms, factorizations, units
    

def example_A3(K=None):
    r"""
    Return a standard rank-`3` FLIR example with five charts.

    This example is modeled on the Cluster Algebra of type `A_3`.
        
    """
    if K is None:
        K = QQ

    base = FLIRChart.base(K, ["x1", "x2", "x3"])
    x1, x2, x3 = base.F.gens()

    chart1 = FLIRChart(K, ["x_11", "x_12", "x_13"],
                       term_order=base.term_order,
                       base_fraction_field=base.F)
    chart1.this_to_base = [x1, (x1 + x3) / x2, x3]

    chart2 = FLIRChart(K, ["x_21", "x_22", "x_23"],
                       term_order=base.term_order,
                       base_fraction_field=base.F)
    chart2.this_to_base = [(x2 + 1) / x1, x2, x3]

    chart3 = FLIRChart(K, ["x_31", "x_32", "x_33"],
                       term_order=base.term_order,
                       base_fraction_field=base.F)
    chart3.this_to_base = [x1, x2, (1 + x2) / x3]

    chart4 = FLIRChart(K, ["x_41", "x_42", "x_43"],
                       term_order=base.term_order,
                       base_fraction_field=base.F)
    chart4.this_to_base = [(x2 + 1) / x1, x2, (x2 + 1) / x3]

    A = FLIR(base, [base, chart1, chart2, chart3, chart4], compute_base_to_charts=True)
    return A

def example_A2_generalized(K=None):

    if K is None:
        K = QQ

    base = FLIRChart.base(K, ["x1", "x2"])
    x1, x2 = base.F.gens()

    chart1 = FLIRChart(K, ["x_11", "x_12"],
                       term_order=base.term_order,
                       base_fraction_field=base.F)
    chart1.this_to_base = [(x2 + 1) / x1, x2]

    chart2 = FLIRChart(K, ["x_21", "x_22"],
                       term_order=base.term_order,
                       base_fraction_field=base.F)
    chart2.this_to_base = [x1, (x1 + 1) ** 2 / x2]

    return FLIR(base, [base, chart1, chart2], compute_base_to_charts=True)
