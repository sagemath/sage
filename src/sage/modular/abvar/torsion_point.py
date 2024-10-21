# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.rings.number_field
"""
Torsion points on modular abelian varieties

AUTHORS:

- William Stein (2007-03)

- Peter Bruin (2014-12): move TorsionPoint to a separate file
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2014 Peter Bruin <P.J.Bruin@math.leidenuniv.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp, rich_to_bool


class TorsionPoint(ModuleElement):
    r"""
    An element of a finite subgroup of a modular abelian variety.

    INPUT:

    - ``parent`` -- a finite subgroup of a modular abelian variety

    - ``element`` -- a `\QQ`-vector space element that represents
      this element in terms of the ambient rational homology

    - ``check`` -- boolean (default: ``True``); whether to check that
      element is in the appropriate vector space

    EXAMPLES:

    The following calls the :class:`TorsionPoint` constructor implicitly::

        sage: J = J0(11)
        sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
        Finite subgroup with invariants [15] over QQbar of Abelian variety J0(11) of dimension 1
        sage: type(G.0)
        <class 'sage.modular.abvar.finite_subgroup.FiniteSubgroup_lattice_with_category.element_class'>
    """
    def __init__(self, parent, element, check=True):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/2,0], [0,1/2]])
            sage: TestSuite(G).run() # long time
        """
        ModuleElement.__init__(self, parent)
        if check:
            if element not in parent.abelian_variety().vector_space():
                raise TypeError("element must be a vector in the abelian variety's rational homology (embedded in the ambient Jacobian product)")
        if element.denominator() == 1:
            element = element.parent().zero_vector()
        self.__element = element

    def element(self):
        r"""
        Return a vector over `\QQ` defining ``self``.

        OUTPUT:

        - A vector in the rational homology of the ambient modular
          Jacobian variety.

        EXAMPLES:

        We create some elements of `J_0(11)`::

            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of
             Abelian variety J0(11) of dimension 1
            sage: G.0.element()
            (1/3, 0)

        The underlying element is a vector over the rational numbers::

            sage: v = (G.0-G.1).element(); v
            (1/3, -1/5)
            sage: type(v)
            <class 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        """
        return self.__element

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        .. NOTE::

            Since they are represented as equivalences classes of
            rational homology modulo integral homology, we represent
            an element corresponding to `v` in the rational homology
            by ``[v]``.

        EXAMPLES::

            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of Abelian variety J0(11) of dimension 1
            sage: G.0._repr_()
            '[(1/3, 0)]'
        """
        return '[%s]' % self.__element

    def _add_(self, other):
        """
        Add two finite subgroup elements with the same parent. This is
        called implicitly by +.

        EXAMPLES::

            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._add_(G.1)
            [(1/3, 1/5)]
            sage: G.0 + G.1
            [(1/3, 1/5)]
        """
        P = self.parent()
        return P.element_class(P, self.__element + other.__element, check=False)

    def _sub_(self, other):
        """
        Subtract two finite subgroup elements with the same parent. This is
        called implicitly by +.

        EXAMPLES::

            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._sub_(G.1)
            [(1/3, -1/5)]
            sage: G.0 - G.1
            [(1/3, -1/5)]
        """
        P = self.parent()
        return P.element_class(P, self.__element - other.__element, check=False)

    def _neg_(self):
        """
        Negate a finite subgroup element.

        EXAMPLES::

            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._neg_()
            [(-1/3, 0)]
        """
        P = self.parent()
        return P.element_class(P, -self.__element, check=False)

    def _rmul_(self, left):
        """
        Left multiply a finite subgroup element by an integer.

        EXAMPLES::

            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._rmul_(2)
            [(2/3, 0)]
            sage: 2*G.0
            [(2/3, 0)]
        """
        P = self.parent()
        return P.element_class(P, left * self.__element, check=False)

    def _lmul_(self, right):
        """
        Right multiply a finite subgroup element by an integer.

        EXAMPLES::

            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._lmul_(2)
            [(2/3, 0)]
            sage: G.0 * 2
            [(2/3, 0)]
        """
        P = self.parent()
        return P.element_class(P, self.__element * right, check=False)

    def _richcmp_(self, right, op):
        """
        Compare ``self`` and ``right``.

        INPUT:

        - ``self``, ``right`` -- elements of the same finite abelian
          variety subgroup

        - ``op`` -- comparison operator (see :mod:`sage.structure.richcmp`)

        OUTPUT: boolean

        EXAMPLES::

            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0 > G.1
            True
            sage: G.0 == G.0
            True
            sage: 3*G.0 == 0
            True
            sage: 3*G.0 == 5*G.1
            True

        We make sure things that should not be equal are not::

            sage: H = J0(14).finite_subgroup([[1/3,0]])
            sage: G.0 == H.0
            False
            sage: G.0
            [(1/3, 0)]
            sage: H.0
            [(1/3, 0)]
        """
        A = self.parent().abelian_variety()
        from sage.rings.rational_field import QQ
        if self.__element.change_ring(QQ) - right.__element.change_ring(QQ) in A.lattice():
            return rich_to_bool(op, 0)
        return richcmp(self.__element, right.__element, op)

    def additive_order(self):
        """
        Return the additive order of ``self``.

        EXAMPLES::

            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0.additive_order()
            3
            sage: G.1.additive_order()
            5
            sage: (G.0 + G.1).additive_order()
            15
            sage: (3*G.0).additive_order()
            1
        """
        return self._relative_element().denominator()

    def _relative_element(self):
        """
        Return coordinates of ``self`` on a basis for the integral
        homology of the containing abelian variety.

        OUTPUT: vector

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: A = J0(43)[1]; A
            Simple abelian subvariety 43b(1,43) of dimension 2 of J0(43)
            sage: C = A.cuspidal_subgroup(); C
            Finite subgroup with invariants [7] over QQ of
             Simple abelian subvariety 43b(1,43) of dimension 2 of J0(43)
            sage: x = C.0; x
            [(0, 1/7, 0, 6/7, 0, 5/7)]
            sage: x._relative_element()
            (0, 1/7, 6/7, 5/7)
        """
        # check=False prevents testing that the element is really in
        # the lattice, not just in the corresponding QQ-vector space.
        return self.parent().abelian_variety().lattice().coordinate_vector(self.__element, check=False)
