r"""
PARI Groups
"""
from __future__ import absolute_import

from sage.groups.old import Group
from sage.libs.all import pari_gen
from sage.rings.all import Integer
from sage.structure.parent import Parent


class PariGroup(Group):
    def __init__(self, x, degree=None):
        """
        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(pari_group=True); G
            PARI group [24, -1, 5, "S4"] of degree 4
            sage: G.category()
            Category of finite groups

        Caveat: fix those tests and/or document precisely that this is
        an abstract group without explicit elements::

            sage: TestSuite(G).run(skip = ["_test_an_element",\
                                           "_test_associativity",\
                                           "_test_elements",\
                                           "_test_elements_eq_reflexive",\
                                           "_test_elements_eq_symmetric",\
                                           "_test_elements_eq_transitive",\
                                           "_test_elements_neq",\
                                           "_test_enumerated_set_contains",\
                                           "_test_enumerated_set_iter_cardinality",\
                                           "_test_enumerated_set_iter_list",\
                                           "_test_inverse",\
                                           "_test_one",\
                                           "_test_prod",\
                                           "_test_some_elements"])
        """
        if not isinstance(x, pari_gen):
            raise TypeError("x (=%s) must be a PARI gen" % x)
        self.__x = x
        self.__degree = degree
        from sage.categories.finite_groups import FiniteGroups
        Parent.__init__(self, category=FiniteGroups())

    def __repr__(self):
        return "PARI group %s of degree %s" % (self.__x, self.__degree)

    def __eq__(self, other):
        """
        Test equality.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f1 = x^4 - 17*x^3 - 2*x + 1
            sage: f2 = x^3 - x - 1
            sage: G1 = f1.galois_group(pari_group=True)
            sage: G2 = f2.galois_group(pari_group=True)
            sage: G1 == G1
            True
            sage: G1 == G2
            False
        """
        return (isinstance(other, PariGroup) and
            (self.__x, self.__degree) == (other.__x, other.__degree))

    def __ne__(self, other):
        """
        Test inequality.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f1 = x^4 - 17*x^3 - 2*x + 1
            sage: f2 = x^3 - x - 1
            sage: G1 = f1.galois_group(pari_group=True)
            sage: G2 = f2.galois_group(pari_group=True)
            sage: G1 != G1
            False
            sage: G1 != G2
            True
        """
        return not (self == other)

    def __pari__(self):
        return self.__x

    def degree(self):
        """
        Return the degree of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f1 = x^4 - 17*x^3 - 2*x + 1
            sage: G1 = f1.galois_group(pari_group=True)
            sage: G1.degree()
            4
        """
        return self.__degree

    def order(self):
        """
        Return the order of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f1 = x^4 - 17*x^3 - 2*x + 1
            sage: G1 = f1.galois_group(pari_group=True)
            sage: G1.order()
            24
        """
        return Integer(self.__x[0])

    def permutation_group(self):
        if self.__degree is None:
            raise NotImplementedError
        from .perm_gps import permgroup_named
        return permgroup_named.TransitiveGroup(self.__degree, self.__x[2])

    _permgroup_ = permutation_group
