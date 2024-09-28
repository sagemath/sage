r"""
Special extensions of function fields

This module currently implements only constant field extension.

Constant field extensions
-------------------------

EXAMPLES:

Constant field extension of the rational function field over rational numbers::

    sage: K.<x> = FunctionField(QQ)
    sage: N.<a> = QuadraticField(2)                                                     # needs sage.rings.number_field
    sage: L = K.extension_constant_field(N)                                             # needs sage.rings.number_field
    sage: L                                                                             # needs sage.rings.number_field
    Rational function field in x over Number Field in a with defining
    polynomial x^2 - 2 with a = 1.4142... over its base
    sage: d = (x^2 - 2).divisor()                                                       # needs sage.libs.pari sage.modules
    sage: d                                                                             # needs sage.libs.pari sage.modules
    -2*Place (1/x)
     + Place (x^2 - 2)
    sage: L.conorm_divisor(d)                                                           # needs sage.libs.pari sage.modules sage.rings.number_field
    -2*Place (1/x)
     + Place (x - a)
     + Place (x + a)

Constant field extension of a function field over a finite field::

    sage: # needs sage.rings.finite_rings sage.rings.function_field
    sage: K.<x> = FunctionField(GF(2)); R.<Y> = K[]
    sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
    sage: E = F.extension_constant_field(GF(2^3))
    sage: E
    Function field in y defined by y^3 + x^6 + x^4 + x^2 over its base
    sage: p = F.get_place(3)
    sage: E.conorm_place(p)  # random
    Place (x + z3, y + z3^2 + z3)
     + Place (x + z3^2, y + z3)
     + Place (x + z3^2 + z3, y + z3^2)
    sage: q = F.get_place(2)
    sage: E.conorm_place(q)  # random
    Place (x + 1, y^2 + y + 1)
    sage: E.conorm_divisor(p + q)  # random
    Place (x + 1, y^2 + y + 1)
     + Place (x + z3, y + z3^2 + z3)
     + Place (x + z3^2, y + z3)
     + Place (x + z3^2 + z3, y + z3^2)

AUTHORS:

- Kwankyu Lee (2021-12-24): added constant field extension
"""

# ****************************************************************************
#       Copyright (C) 2021-2022 Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.ring_extension import RingExtension_generic

from .constructor import FunctionField


class FunctionFieldExtension(RingExtension_generic):
    """
    Abstract base class of function field extensions.
    """
    pass


class ConstantFieldExtension(FunctionFieldExtension):
    """
    Constant field extension.

    INPUT:

    - ``F`` -- a function field whose constant field is `k`

    - ``k_ext`` -- an extension of `k`
    """
    def __init__(self, F, k_ext):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); R.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: E = F.extension_constant_field(GF(2^3))
            sage: TestSuite(E).run(skip=['_test_elements', '_test_pickling'])
        """
        k = F.constant_base_field()
        F_base = F.base_field()

        F_ext_base = FunctionField(k_ext, F_base.variable_name())

        if F.degree() > 1:
            # construct constant field extension F_ext of F
            def_poly = F.polynomial().base_extend(F_ext_base)
            F_ext = F_ext_base.extension(def_poly, names=def_poly.variable_name())
        else: # rational function field
            F_ext = F_ext_base

        # embedding of F into F_ext
        embedk = k_ext.coerce_map_from(k)
        embedF_base = F_base.hom(F_ext_base.gen(), embedk)

        if F.degree() > 1:
            embedF = F.hom(F_ext.gen(), embedF_base)
        else:
            embedF = embedF_base

        self._embedk = embedk
        self._embedF = embedF
        self._F_ext = F_ext
        self._k = k

        super().__init__(embedF, is_backend_exposed=True)

    def top(self):
        """
        Return the top function field of this extension.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); R.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: E = F.extension_constant_field(GF(2^3))
            sage: E.top()
            Function field in y defined by y^3 + x^6 + x^4 + x^2
        """
        return self._F_ext

    def defining_morphism(self):
        """
        Return the defining morphism of this extension.

        This is the morphism from the base to the top.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); R.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: E = F.extension_constant_field(GF(2^3))
            sage: E.defining_morphism()
            Function Field morphism:
              From: Function field in y defined by y^3 + x^6 + x^4 + x^2
              To:   Function field in y defined by y^3 + x^6 + x^4 + x^2
              Defn: y |--> y
                    x |--> x
                    1 |--> 1
        """
        return self._embedF

    def conorm_place(self, p):
        """
        Return the conorm of the place `p` in this extension.

        INPUT:

        - ``p`` -- place of the base function field

        OUTPUT: divisor of the top function field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); R.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: E = F.extension_constant_field(GF(2^3))
            sage: p = F.get_place(3)
            sage: d = E.conorm_place(p)
            sage: [pl.degree() for pl in d.support()]
            [1, 1, 1]
            sage: p = F.get_place(2)
            sage: d = E.conorm_place(p)
            sage: [pl.degree() for pl in d.support()]
            [2]
        """
        embedF = self.defining_morphism()

        O_ext = self.maximal_order()
        Oinf_ext = self.maximal_order_infinite()

        if p.is_infinite_place():
            ideal = Oinf_ext.ideal([embedF(g) for g in p.prime_ideal().gens()])
        else:
            ideal = O_ext.ideal([embedF(g) for g in p.prime_ideal().gens()])

        return ideal.divisor()

    def conorm_divisor(self, d):
        """
        Return the conorm of the divisor ``d`` in this extension.

        INPUT:

        - ``d`` -- divisor of the base function field

        OUTPUT: a divisor of the top function field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); R.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: E = F.extension_constant_field(GF(2^3))
            sage: p1 = F.get_place(3)
            sage: p2 = F.get_place(2)
            sage: c = E.conorm_divisor(2*p1 + 3*p2)
            sage: c1 = E.conorm_place(p1)
            sage: c2 = E.conorm_place(p2)
            sage: c == 2*c1 + 3*c2
            True
        """
        div_top = self.divisor_group()

        c = div_top.zero()
        for pl, mul in d.list():
            c += mul * self.conorm_place(pl)
        return c
