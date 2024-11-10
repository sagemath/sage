# sage.doctest: optional - gap_package_quagroup sage.combinat sage.libs.gap sage.modules
"""
Quantum Groups Using GAP's QuaGroup Package

AUTHORS:

- Travis Scrimshaw (03-2017): initial version

See the :gap_package:`documentation for GAP's QuaGroup package <quagroup/doc/chap0_mj.html>`,
originally authored by Willem Adriaan de Graaf.
"""

# ****************************************************************************
#  Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import re

from copy import copy

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject
from sage.structure.richcmp import op_EQ, op_NE, richcmp
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.family import Family
from sage.combinat.root_system.cartan_type import CartanType
from sage.libs.gap.libgap import libgap
from sage.features.gap import GapPackage
from sage.rings.rational_field import QQ
from sage.categories.algebras import Algebras
from sage.categories.cartesian_product import cartesian_product
from sage.categories.fields import Fields
from sage.categories.homset import HomsetWithBase, Hom
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.modules import Modules
from sage.categories.morphism import Morphism
from sage.categories.rings import Rings

lazy_import('sage.graphs.digraph', 'DiGraph')


class QuaGroupModuleElement(Element):
    """
    Base class for elements created using QuaGroup.
    """
    def __init__(self, parent, libgap_elt):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: TestSuite(Q.an_element()).run()
        """
        self._libgap = libgap(libgap_elt)
        Element.__init__(self, parent)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: Q.an_element()
            1 + (q)*F[a1] + E[a1] + (q^2-1-q^-2 + q^-4)*[ K1 ; 2 ] + K1
             + (-q^-1 + q^-3)*K1[ K1 ; 1 ]

            sage: Q = QuantumGroup(['D',4])
            sage: Q.F_simple()
            Finite family {1: F[a1], 2: F[a2], 3: F[a3], 4: F[a4]}
        """
        # We add some space between the terms
        # FIXME: This doesn't work to avoid within the () for the coeff's
        c = re.compile(r"\+(?! [^(]* \))")
        ret = re.sub(c, ' + ', repr(self._libgap))
        # Replace Ei and Fi with the corresponding root in short form.
        # Do the largest index first so, e.g., F12 gets replaced as 12
        #   instead of as 1.
        for i, al in reversed(list(enumerate(self.parent()._pos_roots))):
            short = '+'.join('%s*a%s' % (coeff, index)
                             if coeff != 1 else 'a%s' % index
                             for index, coeff in al)
            ret = ret.replace('F%s' % (i + 1), 'F[%s]' % short)
            ret = ret.replace('E%s' % (i + 1), 'E[%s]' % short)
        return ret

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: latex(Q.an_element())
            1+{(q)} F_{\alpha_{1}}+E_{\alpha_{1}}+{(q^{2}-1-q^{-2}+q^{-4})}
             [ K_{1} ; 2 ]+K_{1}+{(-q^{-1}+q^{-3})} K_{1}[ K_{1} ; 1 ]

            sage: Q = QuantumGroup(['D',4])
            sage: latex(list(Q.F_simple()))
            \left[F_{\alpha_{1}}, F_{\alpha_{2}},
             F_{\alpha_{3}}, F_{\alpha_{4}}\right]
        """
        from sage.misc.latex import latex
        ret = repr(self._libgap)
        # Do the largest index first so, e.g., F12 gets replaced as 12
        #   instead of as 1.
        for i, al in reversed(list(enumerate(self.parent()._pos_roots))):
            ret = ret.replace('F%s' % (i + 1), 'F_{%s}' % latex(al))
            ret = ret.replace('E%s' % (i + 1), 'E_{%s}' % latex(al))
        for i, ii in reversed(list(enumerate(self.parent()._cartan_type.index_set()))):
            ret = ret.replace('K%s' % (i + 1), 'K_{%s}' % ii)
        # Fugly string parsing to get good looking latex
        # TODO: Find a better way
        ret = ret.replace('(', '{(')
        ret = ret.replace(')', ')}')
        ret = ret.replace('v0', 'v_0')
        ret = ret.replace('*', ' ')
        c = re.compile(r"q\^-?[0-9]*")
        for m in reversed(list(c.finditer(ret))):
            ret = ret[:m.start() + 2] + '{' + ret[m.start() + 2:m.end()] + '}' + ret[m.end():]
        return ret

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: x = Q.an_element()
            sage: loads(dumps(x)) == x
            True
        """
        data = self._libgap.ExtRepOfObj()
        R = self.base_ring()
        ret = []
        for i in range(len(data) // 2):
            ret.append(data[2 * i].sage())
            ret.append(R(str(data[2 * i + 1])))
        return (_unpickle_generic_element, (self.parent(), ret))

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',3])
            sage: x = Q.an_element()
            sage: hash(x) == hash(x.gap())
            True
        """
        return hash(self._libgap)

    def _richcmp_(self, other, op):
        """
        Rich comparison of ``self`` and ``other`` by ``op``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: x = Q.an_element()
            sage: F1, F12, F2 = Q.F()
            sage: q = Q.q()
            sage: x == F1
            False
            sage: x != F1
            True
            sage: F2 * F1
            (q)*F[a1]*F[a2] + F[a1+a2]
            sage: F2 * F1 == q * F1 * F2 + F12
            True
        """
        return richcmp(self._libgap, other._libgap, op)

    def gap(self):
        r"""
        Return the gap representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',3])
            sage: x = Q.an_element()
            sage: x.gap()
            1+(q)*F1+E1+(q^4-1-q^-4+q^-8)*[ K1 ; 2 ]+K1+(-q^-2+q^-6)*K1[ K1 ; 1 ]
        """
        return self._libgap

    _libgap_ = _gap_ = gap

    def _add_(self, other):
        r"""
        Add ``self`` and ``other``.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: F1, F2 = Q.F_simple()
            sage: F1 * F2 + F2 * F1
            (q^3 + 1)*F[a1]*F[a2] + F[a1+a2]
        """
        return self.__class__(self.parent(), self._libgap + other._libgap)

    def _sub_(self, other):
        r"""
        Subtract ``self`` and ``other``.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: F1, F2 = Q.F_simple()
            sage: F1 * F2 - F2 * F1
            (-q^3 + 1)*F[a1]*F[a2] + (-1)*F[a1+a2]
        """
        return self.__class__(self.parent(), self._libgap - other._libgap)

    def _acted_upon_(self, scalar, self_on_left=True):
        r"""
        Return the action of ``scalar`` on ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: q = Q.q()
            sage: x = Q.one().f_tilde([1,2,1,1,2,2]); x
            F[a1+a2]^(3)
            sage: 3 * x
            (3)*F[a1+a2]^(3)
            sage: x * (5/3)
            (5/3)*F[a1+a2]^(3)
            sage: q^-10 * x
            (q^-10)*F[a1+a2]^(3)
            sage: (1 + q^2 - q^-1) * x
            (q^2 + 1-q^-1)*F[a1+a2]^(3)
        """
        try:
            scalar = self.parent().base_ring()(scalar)
            scalar = scalar.subs(q=self.parent()._libgap_q)
        except (TypeError, ValueError):
            return None
        return self.__class__(self.parent(), self._libgap * libgap(scalar))

    def e_tilde(self, i):
        r"""
        Return the action of the Kashiwara operator
        `\widetilde{e}_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set or a list to
          perform a string of operators

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: x = Q.one().f_tilde([1,2,1,1,2,2])
            sage: x.e_tilde([2,2,1,2])
            F[a1]^(2)
        """
        # Do not override this method, instead implement _et
        if isinstance(i, (list, tuple)):
            ret = self
            for j in i:
                if not ret:  # ret == 0
                    return ret
                ret = ret._et(j)
            return ret
        return self._et(i)

    def f_tilde(self, i):
        r"""
        Return the action of the Kashiwara operator
        `\widetilde{f}_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set or a list to
          perform a string of operators

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: Q.one().f_tilde(1)
            F[a1]
            sage: Q.one().f_tilde(2)
            F[a2]
            sage: Q.one().f_tilde([1,2,1,1,2])
            F[a1]*F[a1+a2]^(2)
        """
        # Do not override this method, instead implement _ft
        if isinstance(i, (list, tuple)):
            ret = self
            for j in i:
                if not ret:  # ret == 0
                    return ret
                ret = ret._ft(j)
            return ret
        return self._ft(i)


class QuantumGroup(UniqueRepresentation, Parent):
    r"""
    A Drinfel'd-Jimbo quantum group (implemented using the optional GAP
    package ``QuaGroup``).

    EXAMPLES:

    We check the quantum Serre relations. We first we import the
    `q`-binomial using the `q`-int for quantum groups::

        sage: from sage.algebras.quantum_groups.q_numbers import q_binomial

    We verify the Serre relations for type `A_2`::

        sage: Q = algebras.QuantumGroup(['A',2])
        sage: F1,F12,F2 = Q.F()
        sage: q = Q.q()
        sage: F1^2*F2 - q_binomial(2,1,q) * F1*F2*F1 + F2*F1^2
        0

    We verify the Serre relations for type `B_2`::

        sage: Q = algebras.QuantumGroup(['B',2])
        sage: F1, F12, F122, F2 = Q.F()
        sage: F1^2*F2 - q_binomial(2,1,q^2) * F1*F2*F1 + F2*F1^2
        0
        sage: (F2^3*F1 - q_binomial(3,1,q) * F2^2*F1*F2
        ....:  + q_binomial(3,2,q) * F2*F1*F2^2 - F1*F2^3)
        0

    REFERENCES:

    - :wikipedia:`Quantum_group`
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, q=None):
        """
        Initialize ``self``.

        TESTS::

            sage: Q = QuantumGroup(['A',2])
            sage: Q is QuantumGroup('A2', None)
            True
        """
        cartan_type = CartanType(cartan_type)
        return super().__classcall__(cls, cartan_type, q)

    def __init__(self, cartan_type, q):
        """
        Initialize ``self``.

        TESTS::

            sage: Q = QuantumGroup(['A',2])
            sage: TestSuite(Q).run()            # long time

            sage: Q = QuantumGroup(['G',2])
            sage: TestSuite(Q).run()            # long time
        """
        self._cartan_type = cartan_type
        GapPackage("QuaGroup", spkg='gap_package_quagroup').require()
        libgap.LoadPackage('QuaGroup')
        R = libgap.eval('RootSystem("%s",%s)' % (cartan_type.type(), cartan_type.rank()))
        Q = self._cartan_type.root_system().root_lattice()
        I = cartan_type.index_set()
        self._pos_roots = [Q.sum_of_terms([(ii, root[i])
                                           for i, ii in enumerate(I)
                                           if root[i] != 0])
                           for root in R.PositiveRootsInConvexOrder().sage()]
        if q is None:
            self._libgap = R.QuantizedUEA()
            self._libgap_q = libgap.eval('_q')
            self._libgap_base = libgap.eval('QuantumField')
            base_field = QQ['q'].fraction_field()
            q = base_field.gen()
        else:
            base_field = q.parent()
            self._libgap = R.QuantizedUEA(base_field, q)
            self._libgap_base = libgap(base_field)
            self._libgap_q = libgap(q)
        self._q = q
        Parent.__init__(self, base=base_field, category=HopfAlgebras(Fields()))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: QuantumGroup(['A',2])
            Quantum Group of type ['A', 2] with q=q
        """
        return "Quantum Group of type {} with q={}".format(self._cartan_type, self._q)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(QuantumGroup(['A',3]))
            U_{q}(A_{3})
            sage: zeta3 = CyclotomicField(3).gen()
            sage: latex(QuantumGroup(['G',2], q=zeta3))
            U_{\zeta_{3}}(G_2)
        """
        from sage.misc.latex import latex
        return "U_{%s}(%s)" % (latex(self._q), latex(self._cartan_type))

    def gap(self):
        """
        Return the gap representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.gap()
            QuantumUEA( <root system of type A2>, Qpar = q )
        """
        return self._libgap

    _libgap_ = _gap_ = gap

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def _element_constructor_(self, elt):
        """
        Construct an element of ``self`` from ``elt``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q(0)
            0
            sage: Q(4)
            (4)*1
            sage: Q(4).parent() is Q
            True
            sage: Q(Q.q()).parent() is Q
            True
            sage: Q(Q.an_element()) == Q.an_element()
            True
        """
        if not elt:
            return self.zero()
        if elt in self.base_ring():
            return elt * self.one()
        return self.element_class(self, elt)

    # Special elements
    # ----------------

    @cached_method
    def one(self):
        """
        Return the multiplicative identity of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.one()
            1
        """
        return self.element_class(self, self._libgap.One())

    @cached_method
    def zero(self):
        """
        Return the multiplicative identity of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.zero()
            0
        """
        return self.element_class(self, self._libgap.ZeroImmutable())

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.gens()
            (F[a1], F[a1+a2], F[a2],
             K1, (-q + q^-1)*[ K1 ; 1 ] + K1,
             K2, (-q + q^-1)*[ K2 ; 1 ] + K2,
             E[a1], E[a1+a2], E[a2])
        """
        return tuple([self.element_class(self, gen)
                      for gen in self._libgap.GeneratorsOfAlgebra()])

    def E(self):
        r"""
        Return the family of generators `\{E_{\alpha}\}_{\alpha \in \Phi}`,
        where `\Phi` is the root system of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: list(Q.E())
            [E[a1], E[a1+a2], E[a1+2*a2], E[a2]]
        """
        N = len(self._pos_roots) + len(self._cartan_type.index_set()) * 2
        d = {al: self.gens()[N + i] for i, al in enumerate(self._pos_roots)}
        return Family(self._pos_roots, d.__getitem__)

    def E_simple(self):
        r"""
        Return the family of generators `\{E_i := E_{\alpha_i}\}_{i \in I}`.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: Q.E_simple()
            Finite family {1: E[a1], 2: E[a2]}
        """
        I = self._cartan_type.index_set()
        gens = self.algebra_generators()
        d = {i: gens['E%s' % i] for i in I}
        return Family(I, d.__getitem__)

    def F(self):
        r"""
        Return the family of generators `\{F_{\alpha}\}_{\alpha \in \Phi}`,
        where `\Phi` is the root system of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: list(Q.F())
            [F[a1], F[3*a1+a2], F[2*a1+a2], F[3*a1+2*a2], F[a1+a2], F[a2]]
        """
        d = {al: self.gens()[i] for i, al in enumerate(self._pos_roots)}
        return Family(self._pos_roots, d.__getitem__)

    def F_simple(self):
        r"""
        Return the family of generators `\{F_i := F_{\alpha_i}\}_{i \in I}`.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: Q.F_simple()
            Finite family {1: F[a1], 2: F[a2]}
        """
        I = self._cartan_type.index_set()
        gens = self.algebra_generators()
        d = {i: gens['F%s' % i] for i in I}
        return Family(I, d.__getitem__)

    def K(self):
        r"""
        Return the family of generators `\{K_i\}_{i \in I}`.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',3])
            sage: Q.K()
            Finite family {1: K1, 2: K2, 3: K3}
            sage: Q.K_inverse()
            Finite family {1: (-q + q^-1)*[ K1 ; 1 ] + K1,
                           2: (-q + q^-1)*[ K2 ; 1 ] + K2,
                           3: (-q + q^-1)*[ K3 ; 1 ] + K3}
        """
        N = len(self._pos_roots)
        I = self._cartan_type.index_set()
        d = {ii: self.gens()[N + 2 * i] for i, ii in enumerate(I)}
        return Family(I, d.__getitem__)

    def K_inverse(self):
        r"""
        Return the family of generators `\{K_i^{-1}\}_{i \in I}`.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',3])
            sage: Q.K_inverse()
            Finite family {1: (-q + q^-1)*[ K1 ; 1 ] + K1,
                           2: (-q + q^-1)*[ K2 ; 1 ] + K2,
                           3: (-q + q^-1)*[ K3 ; 1 ] + K3}
        """
        N = len(self._pos_roots)
        I = self._cartan_type.index_set()
        d = {ii: self.gens()[N + 2 * i + 1] for i, ii in enumerate(I)}
        return Family(I, d.__getitem__)

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: list(Q.algebra_generators())
            [F[a1], F[a2],
             K1, K2,
             (-q + q^-1)*[ K1 ; 1 ] + K1, (-q + q^-1)*[ K2 ; 1 ] + K2,
             E[a1], E[a2]]
        """
        I = self._cartan_type.index_set()
        simples = self._cartan_type.root_system().root_lattice().simple_roots()
        ret = {}
        for i, al in enumerate(simples):
            ii = I[i]
            ret['F%s' % ii] = self.F()[al]
            ret['K%s' % ii] = self.K()[ii]
            ret['Ki%s' % ii] = self.K_inverse()[ii]
            ret['E%s' % ii] = self.E()[al]
        keys = (['F%s' % i for i in I] + ['K%s' % i for i in I]
                + ['Ki%s' % i for i in I] + ['E%s' % i for i in I])
        return Family(keys, ret.__getitem__)

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.an_element()
            1 + (q)*F[a1] + E[a1] + (q^2-1-q^-2 + q^-4)*[ K1 ; 2 ]
             + K1 + (-q^-1 + q^-3)*K1[ K1 ; 1 ]
        """
        i = self._cartan_type.index_set()[0]
        al = self._cartan_type.root_system().root_lattice().simple_root(i)
        return self.E()[al] + self.K()[i] + self.K_inverse()[i]**2 + self.q()*self.F()[al]

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: Q.some_elements()
            [1 + (q)*F[a1] + E[a1] + (q^2-1-q^-2 + q^-4)*[ K1 ; 2 ]
              + K1 + (-q^-1 + q^-3)*K1[ K1 ; 1 ],
             K1, F[a1], E[a1]]
        """
        return ([self.an_element()] + list(self.K())
                + list(self.F_simple()) + list(self.E_simple()))

    def q(self):
        """
        Return the parameter `q`.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',3])
            sage: Q.q()
            q
            sage: zeta3 = CyclotomicField(3).gen()
            sage: Q = QuantumGroup(['B',2], q=zeta3)
            sage: Q.q()
            zeta3
        """
        return self._q

    # Misc
    # ----

    def _Hom_(self, Y, category):
        """
        Return the highest weight module of weight ``weight`` of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: H = Hom(Q, B); H
            Set of Morphisms from Quantum Group of type ['A', 2] with q=q to
             Lower Half of Quantum Group of type ['A', 2] with q=q in Category of rings
            sage: type(H)
            <class '...QuantumGroupHomset_with_category_with_equality_by_id'>
        """
        if category is not None and not category.is_subcategory(Rings()):
            raise TypeError("%s is not a subcategory of Rings()" % category)
        if Y not in Rings():
            raise TypeError("%s is not a ring" % Y)
        return QuantumGroupHomset(self, Y, category=category)

    def highest_weight_module(self, weight):
        """
        Return the highest weight module of weight ``weight`` of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.highest_weight_module([1,3])
            Highest weight module of weight Lambda[1] + 3*Lambda[2] of
             Quantum Group of type ['A', 2] with q=q
        """
        return HighestWeightModule(self, weight)

    def lower_half(self):
        """
        Return the lower half of the quantum group ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.lower_half()
            Lower Half of Quantum Group of type ['A', 2] with q=q
        """
        return LowerHalfQuantumGroup(self)

    # Hopf structure
    # --------------

    def coproduct(self, elt, n=1):
        r"""
        Return the coproduct of ``elt`` (iterated ``n`` times).

        The comultiplication `\Delta \colon U_q(\mathfrak{g}) \to
        U_q(\mathfrak{g}) \otimes U_q(\mathfrak{g})` is defined by

        .. MATH::

            \begin{aligned}
            \Delta(E_i) &= E_i \otimes 1 + K_i \otimes E_i, \\
            \Delta(F_i) &= F_i \otimes K_i^{-1} + 1 \otimes F_i, \\
            \Delta(K_i) &= K_i \otimes K_i.
            \end{aligned}

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: [Q.coproduct(e) for e in Q.E()]
            [1*(E[a1]<x>1) + 1*(K1<x>E[a1]),
             1*(E[a1+a2]<x>1) + 1*(K1*K2<x>E[a1+a2]) + q^2-q^-2*(K2*E[a1]<x>E[a2]),
             q^4-q^2-1 + q^-2*(E[a1]<x>E[a2]^(2)) + 1*(E[a1+2*a2]<x>1)
              + 1*(K1<x>E[a1+2*a2]) + q-q^-1*(K1*K2[ K2 ; 1 ]<x>E[a1+2*a2])
              + q-q^-1*(K2*E[a1+a2]<x>E[a2]) + q^5-2*q^3
              + 2*q^-1-q^-3*(K2[ K2 ; 1 ]*E[a1]<x>E[a2]^(2)),
             1*(E[a2]<x>1) + 1*(K2<x>E[a2])]
            sage: [Q.coproduct(f, 2) for f in Q.F_simple()]
            [1*(1<x>1<x>F[a1]) + -q^2 + q^-2*(1<x>F[a1]<x>[ K1 ; 1 ])
              + 1*(1<x>F[a1]<x>K1) + q^4-2 + q^-4*(F[a1]<x>[ K1 ; 1 ]<x>[ K1 ; 1 ])
              + -q^2 + q^-2*(F[a1]<x>[ K1 ; 1 ]<x>K1) + -q^2
              + q^-2*(F[a1]<x>K1<x>[ K1 ; 1 ]) + 1*(F[a1]<x>K1<x>K1),
             1*(1<x>1<x>F[a2]) + -q + q^-1*(1<x>F[a2]<x>[ K2 ; 1 ])
              + 1*(1<x>F[a2]<x>K2) + q^2-2 + q^-2*(F[a2]<x>[ K2 ; 1 ]<x>[ K2 ; 1 ])
              + -q + q^-1*(F[a2]<x>[ K2 ; 1 ]<x>K2) + -q
              + q^-1*(F[a2]<x>K2<x>[ K2 ; 1 ]) + 1*(F[a2]<x>K2<x>K2)]
        """
        D = self._libgap.ComultiplicationMap(n+1)
        # TODO: This is not the correct parent. Need to create it.
        return self.element_class(self, libgap.Image(D, elt._libgap))

    def antipode(self, elt):
        r"""
        Return the antipode of ``elt``.

        The antipode `S \colon U_q(\mathfrak{g}) \to U_q(\mathfrak{g})`
        is the anti-automorphism defined by

        .. MATH::

            S(E_i) = -K_i^{-1}E_i, \qquad
            S(F_i) = -F_iK_i, \qquad
            S(K_i) = K_i^{-1}.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: [Q.antipode(f) for f in Q.F()]
            [(-1)*F[a1]*K1,
             (-q^6 + q^2)*F[a1]*F[a2]*K1*K2 + (-q^4)*F[a1+a2]*K1*K2,
             (-q^8 + q^6 + q^4-q^2)*F[a1]*F[a2]^(2)*K1
              + (-q^9 + 2*q^7-2*q^3 + q)*F[a1]*F[a2]^(2)*K1*K2[ K2 ; 1 ]
              + (-q^5 + q^3)*F[a1+a2]*F[a2]*K1
              + (-q^6 + 2*q^4-q^2)*F[a1+a2]*F[a2]*K1*K2[ K2 ; 1 ]
              + (-q^4)*F[a1+2*a2]*K1 + (-q^5 + q^3)*F[a1+2*a2]*K1*K2[ K2 ; 1 ],
             (-1)*F[a2]*K2]
        """
        S = self._libgap.AntipodeMap()
        return self.element_class(self, libgap.Image(S, elt._libgap))

    def counit(self, elt):
        r"""
        Return the counit of ``elt``.

        The counit `\varepsilon \colon U_q(\mathfrak{g}) \to \QQ(q)` is
        defined by

        .. MATH::

            \varepsilon(E_i) = \varepsilon(F_i) = 0, \qquad
            \varepsilon(K_i) = 1.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: x = Q.an_element()^2
            sage: Q.counit(x)
            4
            sage: Q.counit(Q.one())
            1
            sage: Q.counit(Q.zero())
            0
        """
        # We need to extract the constant coefficient because the
        #   counit in QuaGroup doesn't support it
        R = self.base_ring()
        ext_rep = list(elt._libgap.ExtRepOfObj())
        constant = R.zero()
        for i in range(len(ext_rep) // 2):
            if ext_rep[2 * i].Length() == 0:
                ext_rep.pop(2 * i)  # Pop the key
                constant = R(str(ext_rep.pop(2 * i)))  # Pop the coefficient
                break
        # To reconstruct, we need the following
        F = self._libgap.FamilyObj().ElementsFamily()
        elt = F.ObjByExtRep(ext_rep)
        co = self._libgap.CounitMap()
        return R(str(co(elt))) + constant

    class Element(QuaGroupModuleElement):
        def _mul_(self, other):
            r"""
            Multiply ``self`` and ``other``.

            EXAMPLES::

                sage: Q = QuantumGroup(['G',2])
                sage: F1, F2 = Q.F_simple()
                sage: F1 * F2 * F1 * F2
                F[a1]*F[a1+a2]*F[a2] + (q^7 + q^5 + q + q^-1)*F[a1]^(2)*F[a2]^(2)
                sage: E1, E2 = Q.E_simple()
                sage: F1 * E1
                F[a1]*E[a1]
                sage: E1 * F1
                F[a1]*E[a1] + [ K1 ; 1 ]
            """
            return self.__class__(self.parent(), self._libgap * other._libgap)

        def bar(self):
            r"""
            Return the bar involution on ``self``.

            The bar involution is defined by

            .. MATH::

                \overline{E_i} = E_i, \qquad\qquad
                \overline{F_i} = F_i, \qquad\qquad
                \overline{K_i} = K_i^{-1}.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: [gen.bar() for gen in Q.gens()]
                [F[a1],
                 (q-q^-1)*F[a1]*F[a2] + F[a1+a2],
                 F[a2],
                 (-q + q^-1)*[ K1 ; 1 ] + K1, K1,
                 (-q + q^-1)*[ K2 ; 1 ] + K2, K2,
                 E[a1],
                 (-q^2 + 1)*E[a1]*E[a2] + (q^2)*E[a1+a2],
                 E[a2]]
            """
            bar = self.parent()._libgap.BarAutomorphism()
            return self.__class__(self.parent(), libgap.Image(bar, self._libgap))

        def omega(self):
            r"""
            Return the action of the `\omega` automorphism on ``self``.

            The `\omega` automorphism is defined by

            .. MATH::

                \omega(E_i) = F_i, \qquad\qquad
                \omega(F_i) = E_i, \qquad\qquad
                \omega(K_i) = K_i^{-1}.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: [gen.omega() for gen in Q.gens()]
                [E[a1],
                 (-q)*E[a1+a2],
                 E[a2],
                 (-q + q^-1)*[ K1 ; 1 ] + K1,
                 K1,
                 (-q + q^-1)*[ K2 ; 1 ] + K2,
                 K2,
                 F[a1],
                 (-q^-1)*F[a1+a2],
                 F[a2]]
            """
            omega = self.parent()._libgap.AutomorphismOmega()
            return self.__class__(self.parent(), libgap.Image(omega, self._libgap))

        def tau(self):
            r"""
            Return the action of the `\tau` anti-automorphism on ``self``.

            The `\tau` anti-automorphism is defined by

            .. MATH::

                \tau(E_i) = E_i, \qquad\qquad
                \tau(F_i) = F_i, \qquad\qquad
                \tau(K_i) = K_i^{-1}.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: [gen.tau() for gen in Q.gens()]
                [F[a1],
                 (-q^2 + 1)*F[a1]*F[a2] + (-q)*F[a1+a2],
                 F[a2],
                 (-q + q^-1)*[ K1 ; 1 ] + K1,
                 K1,
                 (-q + q^-1)*[ K2 ; 1 ] + K2,
                 K2,
                 E[a1],
                 (q-q^-1)*E[a1]*E[a2] + (-q)*E[a1+a2],
                 E[a2]]
            """
            tau = self.parent()._libgap.AntiAutomorphismTau()
            return self.__class__(self.parent(), libgap.Image(tau, self._libgap))

        def braid_group_action(self, braid):
            r"""
            Return the action of the braid group element ``braid``.

            The braid group operator `T_i \colon U_q(\mathfrak{g}) \to
            U_q(\mathfrak{g})` is defined by

            .. MATH::

                \begin{aligned}
                T_i(E_i) &= -F_iK_i, \\
                T_i(E_j) &= \sum_{k=0}^{-a_{ij}} (-1)^k q_i^{-k} E_i^{(-a_{ij}-k)} E_j E_i^{(k)} \text{ if } i \neq j,\\
                T_i(K_j) &= K_jK_i^{a_{ij}}, \\
                T_i(F_i) &= -K_i^{-1}E_i, \\
                T_i(F_j) &= \sum_{k=0}^{-a_{ij}} (-1)^k q_i^{-k} F_i^{(k)} F_j F_i^{(-a_{ij}-k)} \text{ if } i \neq j,
                \end{aligned}

            where `a_{ij} = \langle \alpha_j, \alpha_i^\vee \rangle` is the
            `(i,j)`-entry of the Cartan matrix associated to `\mathfrak{g}`.

            INPUT:

            - ``braid`` -- a reduced word of a braid group element

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: F1 = Q.F_simple()[1]
                sage: F1.braid_group_action([1])
                (q-q^-1)*[ K1 ; 1 ]*E[a1] + (-1)*K1*E[a1]
                sage: F1.braid_group_action([1,2])
                F[a2]
                sage: F1.braid_group_action([2,1])
                (-q^3 + 3*q-3*q^-1 + q^-3)*[ K1 ; 1 ]*[ K2 ; 1 ]*E[a1]*E[a2]
                 + (q^3-2*q + q^-1)*[ K1 ; 1 ]*[ K2 ; 1 ]*E[a1+a2]
                 + (q^2-2 + q^-2)*[ K1 ; 1 ]*K2*E[a1]*E[a2]
                 + (-q^2 + 1)*[ K1 ; 1 ]*K2*E[a1+a2]
                 + (q^2-2 + q^-2)*K1*[ K2 ; 1 ]*E[a1]*E[a2]
                 + (-q^2 + 1)*K1*[ K2 ; 1 ]*E[a1+a2]
                 + (-q + q^-1)*K1*K2*E[a1]*E[a2] + (q)*K1*K2*E[a1+a2]
                sage: F1.braid_group_action([1,2,1]) == F1.braid_group_action([2,1,2])
                True
                sage: F1.braid_group_action([]) == F1
                True
            """
            if not braid:
                return self
            QU = self.parent()._libgap
            tau = QU.AntiAutomorphismTau()
            ret = QU.IdentityMapping()
            for i in braid:
                if i < 0:
                    i = -i
                    T = QU.AutomorphismTalpha(i)
                    ret *= tau * T * tau
                else:
                    ret *= QU.AutomorphismTalpha(i)
            return self.__class__(self.parent(), libgap.Image(ret, self._libgap))

        def _et(self, i):
            r"""
            Return the action of the Kashiwara operator `\widetilde{e}_i`
            on ``self``.

            EXAMPLES::

                sage: Q = QuantumGroup(['G',2])
                sage: [(g.e_tilde(1), g.e_tilde(2)) for g in Q.F()]
                [(1, 0), (0, F[a1]^(3)), (0, F[a1]^(2)),
                 (0, F[3*a1+a2]), (0, F[a1]), (0, 1)]

            TESTS::

                sage: Q = QuantumGroup(['A',2])
                sage: Q.one()._et(1)
                0
                sage: Q.zero().e_tilde(1)
                0
            """
            if not self:  # self == 0
                return self
            ret = self._libgap.Ealpha(i)
            if not ret:
                return self.parent().zero()
            return self.__class__(self.parent(), ret)

        def _ft(self, i):
            r"""
            Return the action of the Kashiwara operator `\widetilde{f}_i`
            on ``self``.

            EXAMPLES::

                sage: Q = QuantumGroup(['G',2])
                sage: [(g._ft(1), g._ft(2)) for g in Q.F()]
                [(F[a1]^(2), F[a1+a2]),
                 (F[a1]*F[3*a1+a2], F[3*a1+2*a2]),
                 (F[a1]*F[2*a1+a2], F[a1+a2]^(2)),
                 (F[a1]*F[3*a1+2*a2], F[a1+a2]^(3)),
                 (F[a1]*F[a1+a2], F[a1+a2]*F[a2]),
                 (F[a1]*F[a2], F[a2]^(2))]
                sage: Q.one().f_tilde([1,2,1,1,2,2])
                F[2*a1+a2]*F[a1+a2]*F[a2]

            TESTS::

                sage: Q = QuantumGroup(['A',2])
                sage: Q.zero().f_tilde(1)
                0
            """
            if not self:  # self == 0
                return self
            ret = self._libgap.Falpha(i)
            if not ret:
                return self.parent().zero()
            return self.__class__(self.parent(), ret)


#####################################################################
# Morphisms

class QuantumGroupMorphism(Morphism):
    r"""
    A morphism whose domain is a quantum group.
    """
    def __init__(self, parent, im_gens, check=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: F, K, Ki, E = Q.gens()
            sage: phi = Q.hom([E, Ki, K, F])
            sage: TestSuite(phi).run(skip='_test_category')
        """
        self._repr_type_str = "Quantum group homomorphism"
        Morphism.__init__(self, parent)
        Q = parent.domain()
        self._im_gens = tuple(im_gens)
        if check and len(im_gens) != len(Q.algebra_generators()):
            raise ValueError("number of images must equal the number of generators")
        self._libgap = Q._libgap.QEAHomomorphism(parent.codomain(), im_gens)

    def __reduce__(self):
        r"""
        For pickling.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: F, K, Ki, E = Q.gens()
            sage: phi = Q.hom([E, Ki, K, F])
            sage: loads(dumps(phi)) == phi
            True
        """
        return (self.parent(), (self._im_gens,))

    def _call_(self, val):
        r"""
        Return the image of ``val`` under ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: F, K, Ki, E = Q.gens()
            sage: phi = Q.hom([E, Ki, K, F])
            sage: phi(F)
            E[a1]
            sage: phi(E*F)
            F[a1]*E[a1]
            sage: phi(F*E)
            F[a1]*E[a1] + [ K1 ; 1 ]
            sage: phi(E*K)
            (-q + q^-1)*F[a1]*[ K1 ; 1 ] + F[a1]*K1
            sage: phi(F*E) == phi(F) * phi(E)
            True
        """
        try:
            return self.codomain()(self._libgap.ImageElm(val))
        except TypeError:
            return self.codomain()(str(self._libgap.ImageElm(val)))

    def __richcmp__(self, other, op):
        r"""
        Rich comparison of ``self`` and ``other`` by ``op``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: F, K, Ki, E = Q.gens()
            sage: phi = Q.hom([E, Ki, K, F])
            sage: psi = Q.hom([F, K, Ki, E])
            sage: phi == Q.hom([E, Ki, K, F])
            True
            sage: phi == psi
            False
            sage: psi != Q.hom([F, K, Ki, E])
            False
            sage: phi != psi
            True

            sage: QB = QuantumGroup(['B',3])
            sage: QC = QuantumGroup(['C',3])
            sage: x = ZZ.one()
            sage: phi = QB.hom([x]*len(QB.algebra_generators()))
            sage: psi = QC.hom([x]*len(QC.algebra_generators()))
            sage: phi.im_gens() == psi.im_gens()
            True
            sage: phi == psi
            False
        """
        if op == op_EQ:
            return (type(self) is type(other)
                    and self.domain() is other.domain()
                    and self._im_gens == other._im_gens)
        if op == op_NE:
            return not (self == other)
        return NotImplemented

    def im_gens(self):
        r"""
        Return the image of the generators under ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: F, K, Ki, E = Q.gens()
            sage: phi = Q.hom([E, Ki, K, F])
            sage: phi.im_gens()
            (E[a1], (-q + q^-1)*[ K1 ; 1 ] + K1, K1, F[a1])
        """
        return self._im_gens

    def _repr_defn(self):
        r"""
        Used in constructing the string representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: F, K, Ki, E = Q.gens()
            sage: phi = Q.hom([E, Ki, K, F])
            sage: print(phi._repr_defn())
            F[a1] |--> E[a1]
            K1 |--> (-q + q^-1)*[ K1 ; 1 ] + K1
            (-q + q^-1)*[ K1 ; 1 ] + K1 |--> K1
            E[a1] |--> F[a1]
        """
        return '\n'.join('%s |--> %s' % (gen, self._im_gens[i])
                         for i, gen in enumerate(self.domain().algebra_generators()))


class QuantumGroupHomset(HomsetWithBase):
    r"""
    The homset whose domain is a quantum group.
    """
    def __call__(self, im_gens, check=True):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: H = Hom(Q, Q)
            sage: F, K, Ki, E = Q.gens()
            sage: phi = H([E, Ki, K, F]); phi
            Quantum group homomorphism endomorphism of Quantum Group of type ['A', 1] with q=q
              Defn: F[a1] |--> E[a1]
                    K1 |--> (-q + q^-1)*[ K1 ; 1 ] + K1
                    (-q + q^-1)*[ K1 ; 1 ] + K1 |--> K1
                    E[a1] |--> F[a1]
            sage: H(phi) == phi
            True
            sage: H2 = Hom(Q, Q, Modules(Fields()))
            sage: H == H2
            False
            sage: H2(phi)
            Quantum group homomorphism endomorphism of Quantum Group of type ['A', 1] with q=q
              Defn: F[a1] |--> E[a1]
                    K1 |--> (-q + q^-1)*[ K1 ; 1 ] + K1
                    (-q + q^-1)*[ K1 ; 1 ] + K1 |--> K1
                    E[a1] |--> F[a1]
        """
        if isinstance(im_gens, QuantumGroupMorphism):
            if im_gens.parent() is self:
                return im_gens
            if im_gens.parent() != self:
                return QuantumGroupMorphism(self, im_gens.im_gens())
            raise TypeError("unable to coerce {}".format(im_gens))
        return QuantumGroupMorphism(self, im_gens)


def projection_lower_half(Q):
    r"""
    Return the projection onto the lower half of the quantum group.

    EXAMPLES::

        sage: from sage.algebras.quantum_groups.quantum_group_gap import projection_lower_half
        sage: Q = QuantumGroup(['G',2])
        sage: phi = projection_lower_half(Q); phi
        Quantum group homomorphism endomorphism of Quantum Group of type ['G', 2] with q=q
          Defn: F[a1] |--> F[a1]
                F[a2] |--> F[a2]
                K1 |--> 0
                K2 |--> 0
                (-q + q^-1)*[ K1 ; 1 ] + K1 |--> 0
                (-q^3 + q^-3)*[ K2 ; 1 ] + K2 |--> 0
                E[a1] |--> 0
                E[a2] |--> 0
        sage: all(phi(f) == f for f in Q.F())
        True
        sage: all(phi(e) == Q.zero() for e in Q.E())
        True
        sage: all(phi(K) == Q.zero() for K in Q.K())
        True
    """
    I = Q._cartan_type.index_set()
    return Hom(Q, Q)(list(Q.F_simple()) + [Q.zero()] * (len(I) * 3))


#####################################################################
# Representations

class QuaGroupRepresentationElement(QuaGroupModuleElement):
    """
    Element of a quantum group representation.
    """
    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: F1, F2 = Q.F_simple()
            sage: q = Q.q()
            sage: V = Q.highest_weight_module([2,1])
            sage: v = V.highest_weight_vector()
            sage: x = (2 - q) * v + F1*v + q*F2*F1*v
            sage: loads(dumps(x)) == x
            True
        """
        return (self.parent(), (self.monomial_coefficients(),))

    def _acted_upon_(self, scalar, self_on_left=False):
        r"""
        Return the action of ``scalar`` on ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['B',2])
            sage: F1, F2 = Q.F_simple()
            sage: q = Q.q()
            sage: V = Q.highest_weight_module([2,1])
            sage: v = V.highest_weight_vector()
            sage: F1 * v
            F[a1]*v0
            sage: F2 * v
            F[a2]*v0
            sage: F1^2 * v
            (q^2 + q^-2)*F[a1]^(2)*v0
            sage: F2^2 * v
            0*v0
            sage: (F1 * F2) * v
            F[a1]*F[a2]*v0
            sage: F1 * (F2 * v)
            F[a1]*F[a2]*v0
            sage: (2 - q) * v + F1*v + q*F2*F1*v
            (-q + 2)*1*v0 + F[a1]*v0 + (q^3)*F[a1]*F[a2]*v0 + (q)*F[a1+a2]*v0
        """
        try:
            if scalar.parent() is self.parent()._Q:
                if self_on_left:  # Only act: scalar * v
                    return None
                return self.__class__(self.parent(), scalar._libgap ** self._libgap)
        except AttributeError:
            pass
        return QuaGroupModuleElement._acted_upon_(self, scalar, self_on_left)

    _lmul_ = _acted_upon_

    def _et(self, i):
        r"""
        Return the action of `\widetilde{e}_i` on ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: v = V.highest_weight_vector()
            sage: v._et(1)
            0*v0
            sage: V.zero().e_tilde(1)
            0*v0
        """
        if not self:  # self == 0
            return self
        V = self.parent()
        ret = V._libgap.Ealpha(self._libgap, i)
        return self.__class__(V, ret)

    def _ft(self, i):
        r"""
        Return the action of `\widetilde{e}_i` on ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['C',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: v = V.highest_weight_vector()
            sage: v._ft(1)
            F[a1]*v0
            sage: v._ft(2)
            F[a2]*v0
            sage: v.f_tilde([1,1])
            0*v0
            sage: v.f_tilde([2,2])
            0*v0
            sage: v.f_tilde([2,1,1])
            (-q^-3)*F[a1]*F[a1+a2]*v0 + (-q^-4)*F[2*a1+a2]*v0
            sage: v.f_tilde([1,2,2])
            F[a1+a2]*F[a2]*v0
            sage: V.zero().f_tilde(1)
            0*v0
        """
        if not self:  # self == 0
            return self
        V = self.parent()
        ret = V._libgap.Falpha(self._libgap, i)
        return self.__class__(V, ret)

    def monomial_coefficients(self, copy=True):
        r"""
        Return the dictionary of ``self`` whose keys are the basis indices
        and the values are coefficients.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: v = V.highest_weight_vector()
            sage: F1, F2 = Q.F_simple()
            sage: q = Q.q()
            sage: x = v + F1*v + q*F2*F1*v; x
            1*v0 + F[a1]*v0 + (q^2)*F[a1]*F[a2]*v0 + (q)*F[a1+a2]*v0
            sage: sorted(x.monomial_coefficients().items(), key=str)
            [(0, 1), (1, 1), (3, q^2), (4, q)]
        """
        R = self.parent()._Q.base_ring()
        B = self.parent()._libgap.Basis()
        data = [R(str(c)) for c in libgap.Coefficients(B, self._libgap)]
        return {i: c for i, c in enumerate(data) if c != 0}

    def _vector_(self, R=None, order=None, sparse=False):
        """
        Return ``self`` as a vector.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: v = V.highest_weight_vector()
            sage: vector(v)
            (1, 0, 0, 0, 0, 0, 0, 0)
            sage: F1, F2 = Q.F_simple()
            sage: q = Q.q()
            sage: x = v + F1*v + q*F2*F1*v; x
            1*v0 + F[a1]*v0 + (q^2)*F[a1]*F[a2]*v0 + (q)*F[a1+a2]*v0
            sage: vector(x)
            (1, 1, 0, q^2, q, 0, 0, 0)

            sage: v._vector_(sparse=True)
            (1, 0, 0, 0, 0, 0, 0, 0)
            sage: x._vector_(sparse=True)
            (1, 1, 0, q^2, q, 0, 0, 0)

            sage: M = V.submodule([V.an_element()])
            sage: M
            Free module generated by {0} over Fraction Field of Univariate Polynomial Ring in q over Rational Field
        """
        V = self.parent()._dense_free_module(R)
        if sparse:
            V = V.sparse_module()
            if order is None:
                return V(self.monomial_coefficients())

        v = copy(V.zero())
        if order is None:
            for i, c in self.monomial_coefficients().items():
                v[i] = c
        else:
            for i, c in self.monomial_coefficients().items():
                v[order[i]] = c
        return v


class CrystalGraphVertex(SageObject):
    r"""
    Helper class used as the vertices of a crystal graph.
    """
    def __init__(self, V, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.quantum_group_gap import CrystalGraphVertex
            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: v = CrystalGraphVertex(V, '<F2*v0>')
            sage: TestSuite(v).run()
        """
        self.V = V
        self.s = s

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.quantum_group_gap import CrystalGraphVertex
            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: v = CrystalGraphVertex(V, '<F2*v0>')
            sage: hash(v) == hash('<F2*v0>')
            True
        """
        return hash(self.s)

    def __eq__(self, other):
        """
        Check equality of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.quantum_group_gap import CrystalGraphVertex
            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: v = CrystalGraphVertex(V, '<F2*v0>')
            sage: vp = CrystalGraphVertex(V, '<F2*v0>')
            sage: v == vp
            True
            sage: vpp = CrystalGraphVertex(V, '<1*v0>')
            sage: v == vpp
            False
        """
        return isinstance(other, CrystalGraphVertex) and self.s == other.s

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.quantum_group_gap import CrystalGraphVertex
            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: CrystalGraphVertex(V, '<F2*v0>')
            <F2*v0>
        """
        return self.s

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.quantum_group_gap import CrystalGraphVertex
            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: v = CrystalGraphVertex(V, '<F2*v0>')
            sage: latex(v)
            \langle F_{\alpha_{1} + \alpha_{2}} v_0 \rangle
        """
        # Essentially same as QuaGroupModuleElement._latex_
        from sage.misc.latex import latex
        ret = self.s[1:-1]  # Strip leading '<' and trailing '>'
        for i, al in enumerate(self.V._pos_roots):
            ret = ret.replace('F%s' % (i + 1), 'F_{%s}' % latex(al))
            ret = ret.replace('E%s' % (i + 1), 'E_{%s}' % latex(al))
        for i, ii in enumerate(self.V._cartan_type.index_set()):
            ret = ret.replace('K%s' % (i + 1), 'K_{%s}' % ii)
        # Fugly string parsing to get good looking latex
        # TODO: Find a better way
        ret = ret.replace('(', '{(')
        ret = ret.replace(')', ')}')
        ret = ret.replace('v0', 'v_0')
        ret = ret.replace('*', ' ')
        ret = ret.replace('<x>', ' \\otimes ')
        c = re.compile(r"q\^-?[0-9]*")
        for m in reversed(list(c.finditer(ret))):
            ret = ret[:m.start()+2]+'{'+ret[m.start()+2:m.end()]+'}'+ret[m.end():]
        return '\\langle {} \\rangle'.format(ret)


class QuantumGroupModule(Parent, UniqueRepresentation):
    r"""
    Abstract base class for quantum group representations.
    """
    def __init__(self, Q, category):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['G',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: TestSuite(V).run()
        """
        self._Q = Q
        self._libgap_q = Q._libgap_q
        self._libgap_base = Q._libgap_base
        self._cartan_type = Q._cartan_type
        self._pos_roots = Q._pos_roots
        Parent.__init__(self, base=Q.base_ring(), category=category)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: S = T.highest_weight_decomposition()[0]
            sage: latex(S)
            \begin{tikzpicture}...
            ...
            \end{tikzpicture}
        """
        from sage.misc.latex import latex
        return latex(self.crystal_graph())

    def gap(self):
        r"""
        Return the gap representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: V.gap()
            <8-dimensional left-module over QuantumUEA( <root system of type A2>,
             Qpar = q )>
        """
        return self._libgap

    _libgap_ = _gap_ = gap

    def _element_constructor_(self, elt):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: q = Q.q()
            sage: V(0)
            0*v0
            sage: V({1: q^2 - q^-2, 3: 2})
            (q^2-q^-2)*F[a1]*v0 + (2)*F[a1]*F[a2]*v0
        """
        if not elt:
            return self.zero()
        if isinstance(elt, dict):
            return self._from_dict(elt)
        return self.element_class(self, elt)

    @cached_method
    def basis(self):
        r"""
        Return a basis of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: V.basis()
            Family (1*v0, F[a1]*v0, F[a2]*v0, F[a1]*F[a2]*v0, F[a1+a2]*v0,
                    F[a1]*F[a1+a2]*v0, F[a1+a2]*F[a2]*v0, F[a1+a2]^(2)*v0)
        """
        return Family([self.element_class(self, b) for b in self._libgap.Basis()])

    @cached_method
    def crystal_basis(self):
        r"""
        Return the crystal basis of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: V.crystal_basis()
            Family (1*v0, F[a1]*v0, F[a2]*v0, F[a1]*F[a2]*v0,
                    (q)*F[a1]*F[a2]*v0 + F[a1+a2]*v0, F[a1+a2]*F[a2]*v0,
                    (-q^-2)*F[a1]*F[a1+a2]*v0, (-q^-1)*F[a1+a2]^(2)*v0)
        """
        return Family([self.element_class(self, b) for b in self._libgap.CrystalBasis()])

    @cached_method
    def R_matrix(self):
        """
        Return the `R`-matrix of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',1])
            sage: V = Q.highest_weight_module([1])
            sage: V.R_matrix()
            [       1        0        0        0]
            [       0        q -q^2 + 1        0]
            [       0        0        q        0]
            [       0        0        0        1]
        """
        R = self._libgap.RMatrix()
        F = self._Q.base_ring()
        from sage.matrix.constructor import matrix
        M = matrix(F, [[F(str(elt)) for elt in row] for row in R])
        M.set_immutable()
        return M

    def crystal_graph(self):
        r"""
        Return the crystal graph of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: G = V.crystal_graph(); G
            Digraph on 8 vertices

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: G.is_isomorphic(B.digraph(), edge_labels=True)
            True
        """
        G = self._libgap.CrystalGraph()
        vertices = [CrystalGraphVertex(self, repr(p)) for p in G['points']]
        edges = [[vertices[e[0][0]-1], vertices[e[0][1]-1], e[1]]
                 for e in G['edges'].sage()]
        G = DiGraph([vertices, edges], format='vertices_and_edges')
        from sage.graphs.dot2tex_utils import have_dot2tex
        if have_dot2tex():
            G.set_latex_options(format='dot2tex',
                                edge_labels=True,
                                color_by_label=self._cartan_type._index_set_coloring)
        return G

    @cached_method
    def zero(self):
        r"""
        Return the zero element of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: V.zero()
            0*v0
        """
        return self.element_class(self, self._libgap.ZeroImmutable())


class HighestWeightModule(QuantumGroupModule):
    """
    A highest weight module of a quantum group.
    """
    @staticmethod
    def __classcall_private__(cls, Q, weight):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: La = Q.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: V = Q.highest_weight_module([1,3])
            sage: V is Q.highest_weight_module(La[1]+3*La[2])
            True
        """
        P = Q._cartan_type.root_system().weight_lattice()
        if isinstance(weight, (list, tuple)):
            La = P.fundamental_weights()
            weight = P.sum(la * weight[i] for i, la in enumerate(La))
        else:
            weight = P(weight)
        return super().__classcall__(cls, Q, weight)

    def __init__(self, Q, weight):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: TestSuite(V).run()
        """
        self._libgap = Q._libgap.HighestWeightModule(list(weight.to_vector()))
        self._weight = weight
        cat = Modules(Q.base_ring()).FiniteDimensional().WithBasis()
        QuantumGroupModule.__init__(self, Q, cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.highest_weight_module([1,1])
            Highest weight module of weight Lambda[1] + Lambda[2] of
             Quantum Group of type ['A', 2] with q=q
        """
        return "Highest weight module of weight {} of {}".format(self._weight, self._Q)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,2])
            sage: latex(V)
            V(\Lambda_{1} + 2 \Lambda_{2})
        """
        from sage.misc.latex import latex
        return "V({})".format(latex(self._weight))

    @cached_method
    def highest_weight_vector(self):
        """
        Return the highest weight vector of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: V.highest_weight_vector()
            1*v0
        """
        return self.element_class(self, self._libgap.HighestWeightsAndVectors()[1][0][0])

    an_element = highest_weight_vector

    def tensor(self, *V, **options):
        """
        Return the tensor product of ``self`` with ``V``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: Vp = Q.highest_weight_module([1,0])
            sage: Vp.tensor(V)
            Highest weight module of weight Lambda[1] of Quantum Group of type ['A', 2] with q=q
             # Highest weight module of weight Lambda[1] + Lambda[2] of Quantum Group of type ['A', 2] with q=q
        """
        return TensorProductOfHighestWeightModules(self, *V, **options)

    Element = QuaGroupRepresentationElement


class TensorProductOfHighestWeightModules(QuantumGroupModule):
    def __init__(self, *modules, **options):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,1])
            sage: T = tensor([V,V])
            sage: TestSuite(T).run()
        """
        Q = modules[0]._Q
        self._modules = tuple(modules)
        self._libgap = libgap.TensorProductOfAlgebraModules([m._libgap for m in modules])
        cat = Modules(Q.base_ring()).TensorProducts().FiniteDimensional().WithBasis()
        QuantumGroupModule.__init__(self, Q, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: T
            Highest weight module of weight Lambda[1] of Quantum Group of type ['A', 2] with q=q
             # Highest weight module of weight Lambda[1] of Quantum Group of type ['A', 2] with q=q
        """
        return " # ".join(repr(M) for M in self._modules)

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: latex(T)
            V(\Lambda_{1}) \otimes V(\Lambda_{1})
        """
        from sage.misc.latex import latex
        return " \\otimes ".join(latex(M) for M in self._modules)

    @lazy_attribute
    def _highest_weights_and_vectors(self):
        """
        Return the highest weights and the corresponding vectors.

        .. NOTE::

            The resulting objects are GAP objects.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([0,1])
            sage: T = tensor([V,V])
            sage: T._highest_weights_and_vectors
            [ [ [ 0, 2 ], [ 1, 0 ] ],
             [ [ 1*(1*v0<x>1*v0) ], [ -q^-1*(1*v0<x>F3*v0)+1*(F3*v0<x>1*v0) ] ] ]
        """
        return self._libgap.HighestWeightsAndVectors()

    def highest_weight_vectors(self):
        r"""
        Return the highest weight vectors of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: T.highest_weight_vectors()
            [1*(1*v0<x>1*v0), -q^-1*(1*v0<x>F[a1]*v0) + 1*(F[a1]*v0<x>1*v0)]
        """
        return [self.element_class(self, v)
                for vecs in self._highest_weights_and_vectors[1]
                for v in vecs]

    some_elements = highest_weight_vectors

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: T.an_element()
            1*(1*v0<x>1*v0)
        """
        return self.highest_weight_vectors()[0]

    @cached_method
    def highest_weight_decomposition(self):
        """
        Return the highest weight decomposition of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: T.highest_weight_decomposition()
            [Highest weight submodule with weight 2*Lambda[1] generated by 1*(1*v0<x>1*v0),
             Highest weight submodule with weight Lambda[2] generated by -q^-1*(1*v0<x>F[a1]*v0) + 1*(F[a1]*v0<x>1*v0)]
        """
        return [HighestWeightSubmodule(self, self.element_class(self, v), tuple(wt.sage()))
                for wt, vecs in zip(*self._highest_weights_and_vectors)
                for v in vecs]

    def tensor_factors(self):
        r"""
        Return the factors of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: T.tensor_factors()
            (Highest weight module of weight Lambda[1] of Quantum Group of type ['A', 2] with q=q,
             Highest weight module of weight Lambda[1] of Quantum Group of type ['A', 2] with q=q)
        """
        return self._modules

    Element = QuaGroupRepresentationElement


class HighestWeightSubmodule(QuantumGroupModule):
    def __init__(self, ambient, gen, weight):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: S = T.highest_weight_decomposition()[0]
            sage: TestSuite(S).run()
        """
        self._ambient = ambient
        # We do not use the generic ambient category since submodules of tensor
        #   products are considered to be tensor products.
        # This should be reverted after this has changed.
        #cat = ambient.category()
        cat = Modules(ambient.base_ring()).FiniteDimensional().WithBasis()
        QuantumGroupModule.__init__(self, ambient._Q, cat.Subobjects())

        self._gen = gen

        self._libgap = self._ambient._libgap.HWModuleByGenerator(gen, weight)

        # Convert the weight to an element of the weight lattice
        P = self._Q._cartan_type.root_system().weight_lattice()
        La = P.fundamental_weights()
        self._weight = P.sum(la * weight[i] for i, la in enumerate(La))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: T.highest_weight_decomposition()
            [Highest weight submodule with weight 2*Lambda[1]
                generated by 1*(1*v0<x>1*v0),
             Highest weight submodule with weight Lambda[2]
                generated by -q^-1*(1*v0<x>F[a1]*v0) + 1*(F[a1]*v0<x>1*v0)]
        """
        return "Highest weight submodule with weight {} generated by {}".format(self._weight, self._gen)

    @lazy_attribute
    def _ambient_basis_map(self):
        """
        A dict that maps the basis of ``self`` to the ambient module.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: S = T.highest_weight_decomposition()[0]
            sage: S._ambient_basis_map
            {0: 1*(1*v0<x>1*v0),
             1: 1*(1*v0<x>F[a1]*v0) + q^-1*(F[a1]*v0<x>1*v0),
             2: 1*(F[a1]*v0<x>F[a1]*v0),
             3: 1*(1*v0<x>F[a1+a2]*v0) + q^-1*(F[a1+a2]*v0<x>1*v0),
             4: 1*(F[a1]*v0<x>F[a1+a2]*v0) + q^-1*(F[a1+a2]*v0<x>F[a1]*v0),
             5: 1*(F[a1+a2]*v0<x>F[a1+a2]*v0)}
        """
        B = list(self.basis())
        d = {self.highest_weight_vector(): self._gen}
        todo = {self.highest_weight_vector()}
        I = self._cartan_type.index_set()
        while todo:
            x = todo.pop()
            for i in I:
                y = x.f_tilde(i)
                if y and y not in d:
                    d[y] = d[x].f_tilde(i)
                    todo.add(y)
        return {B.index(k): d[k] for k in d}

    def ambient(self):
        """
        Return the ambient module of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: S = T.highest_weight_decomposition()[0]
            sage: S.ambient() is T
            True
        """
        return self._ambient

    @lazy_attribute
    def lift(self):
        """
        The lift morphism from ``self`` to the ambient space.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: S = T.highest_weight_decomposition()[0]
            sage: S.lift
            Generic morphism:
              From: Highest weight submodule with weight 2*Lambda[1] generated by 1*(1*v0<x>1*v0)
              To:   Highest weight module ... # Highest weight module ...
            sage: x = sum(S.basis())
            sage: x.lift()
            1*(1*v0<x>1*v0) + 1*(1*v0<x>F[a1]*v0) + 1*(1*v0<x>F[a1+a2]*v0)
             + q^-1*(F[a1]*v0<x>1*v0) + 1*(F[a1]*v0<x>F[a1]*v0)
             + 1*(F[a1]*v0<x>F[a1+a2]*v0) + q^-1*(F[a1+a2]*v0<x>1*v0)
             + q^-1*(F[a1+a2]*v0<x>F[a1]*v0) + 1*(F[a1+a2]*v0<x>F[a1+a2]*v0)
        """
        return self.module_morphism(self._ambient_basis_map.__getitem__,
                                    codomain=self._ambient, unitriangular='lower')

    def retract(self, elt):
        """
        The retract map from the ambient space to ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: all(S.retract(S.lift(x)) == x
            ....:     for S in T.highest_weight_decomposition()
            ....:     for x in S.basis())
            True
        """
        c = self.lift.matrix().solve_right(elt._vector_())
        return self._from_dict(c.dict(), coerce=False, remove_zeros=False)

    def highest_weight_vector(self):
        """
        Return the highest weight vector of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: S = T.highest_weight_decomposition()[1]
            sage: u = S.highest_weight_vector(); u
            (1)*e.1
            sage: u.lift()
            -q^-1*(1*v0<x>F[a1]*v0) + 1*(F[a1]*v0<x>1*v0)
        """
        I = self._cartan_type.index_set()
        zero = self._libgap.ZeroImmutable()
        for v in self.basis():
            if all(self._libgap.Ealpha(v._libgap, i) == zero for i in I):
                return v
        return self.zero()

    an_element = highest_weight_vector

    def crystal_graph(self, use_ambient=True):
        """
        Return the crystal graph of ``self``.

        INPUT:

        - ``use_ambient`` -- boolean (default: ``True``); if ``True``,
          the vertices are given in terms of the ambient module

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: V = Q.highest_weight_module([1,0])
            sage: T = tensor([V,V])
            sage: S = T.highest_weight_decomposition()[1]
            sage: G = S.crystal_graph()
            sage: sorted(G.vertices(sort=False), key=str)
            [<-q^-1*(1*v0<x>F[a1+a2]*v0) + 1*(F[a1+a2]*v0<x>1*v0)>,
             <-q^-1*(1*v0<x>F[a1]*v0) + 1*(F[a1]*v0<x>1*v0)>,
             <-q^-1*(F[a1]*v0<x>F[a1+a2]*v0) + 1*(F[a1+a2]*v0<x>F[a1]*v0)>]
            sage: sorted(S.crystal_graph(False).vertices(sort=False), key=str)
            [<(1)*e.1>, <(1)*e.2>, <(1)*e.3>]
        """
        G = self._libgap.CrystalGraph()
        if not use_ambient:
            return QuantumGroupModule.crystal_graph(self)
        # Mostly a copy; there is likely a better way with a helper function
        B = self.basis()
        d = {repr(B[k]._libgap): '<{!r}>'.format(self._ambient_basis_map[k])
             for k in self._ambient_basis_map}
        vertices = [CrystalGraphVertex(self, d[repr(p)[1:-1]])
                    for p in G['points']]
        edges = [[vertices[e[0][0]-1], vertices[e[0][1]-1], e[1]]
                 for e in G['edges'].sage()]
        G = DiGraph([vertices, edges], format='vertices_and_edges')
        from sage.graphs.dot2tex_utils import have_dot2tex
        if have_dot2tex():
            G.set_latex_options(format='dot2tex',
                                edge_labels=True,
                                color_by_label=self._cartan_type._index_set_coloring)
        return G

    Element = QuaGroupRepresentationElement


# TODO: Generalized this to Verma modules
class LowerHalfQuantumGroup(Parent, UniqueRepresentation):
    """
    The lower half of the quantum group.
    """
    @staticmethod
    def __classcall_private__(cls, Q):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.quantum_group_gap import LowerHalfQuantumGroup
            sage: Q = QuantumGroup(['A',2])
            sage: Q.lower_half() is LowerHalfQuantumGroup(Q)
            True
        """
        from sage.combinat.root_system.cartan_type import CartanType_abstract
        if isinstance(Q, CartanType_abstract):
            Q = QuantumGroup(Q)
        return super().__classcall__(cls, Q)

    def __init__(self, Q):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: TestSuite(B).run()
        """
        self._Q = Q
        self._libgap = Q._libgap
        self._libgap_q = Q._libgap_q
        self._libgap_base = Q._libgap_base
        self._cartan_type = Q._cartan_type
        self._pos_roots = Q._pos_roots
        self._proj = projection_lower_half(Q)
        B = Q.base_ring()
        Parent.__init__(self, base=B, category=Algebras(B).WithBasis().Subobjects())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: Q.lower_half()
            Lower Half of Quantum Group of type ['A', 2] with q=q
        """
        return "Lower Half of {}".format(self._Q)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: latex(Q.lower_half())
            U^-_{q}(A_{2})
        """
        from sage.misc.latex import latex
        return "U^-_{%s}(%s)" % (latex(self._Q._q), latex(self._cartan_type))

    def _element_constructor_(self, elt):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: q = Q.q()
            sage: B(0)
            0
            sage: B(1 + q^2)
            (q^2 + 1)*1
            sage: B({(1,2,0): q, (0,0,2): q^2 - 2})
            (q)*F[a1]*F[a1+a2]^(2) + (q^2-2)*F[a2]^(2)
        """
        if not elt:
            return self.zero()
        if isinstance(elt, dict):
            return self._from_dict(elt)
        if elt in self.base_ring():
            return elt * self.one()
        if elt.parent() is self._Q:
            return self.element_class(self, self._proj(elt)._libgap)
        return self.element_class(self, elt)

    def ambient(self):
        r"""
        Return the ambient quantum group of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: B.ambient() is Q
            True
        """
        return self._Q

    @cached_method
    def highest_weight_vector(self):
        """
        Return the highest weight vector of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: B.highest_weight_vector()
            1
        """
        return self.element_class(self, self._Q.one()._libgap)

    one = highest_weight_vector
    an_element = highest_weight_vector

    @cached_method
    def zero(self):
        """
        Return the zero element of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: B.zero()
            0
        """
        return self.element_class(self, self._Q._libgap.ZeroImmutable())

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: B.algebra_generators()
            Finite family {1: F[a1], 2: F[a2]}
        """
        F = self._Q.F_simple()
        keys = F.keys()
        d = {i: self.element_class(self, F[i]._libgap) for i in keys}
        return Family(keys, d.__getitem__)

    gens = algebra_generators

    def _construct_monomial(self, k):
        """
        Construct a monomial of ``self`` indexed by ``k``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: B._construct_monomial((1,2,1))
            F[a1]*F[a1+a2]^(2)*F[a2]
            sage: B._construct_monomial((3,0,1))
            F[a1]^(3)*F[a2]
        """
        F = self._libgap.FamilyObj().ElementsFamily()
        one = self._libgap_base.One()
        data = []
        for i, val in enumerate(k):
            if val == 0:
                continue
            data.append(i + 1)
            data.append(val)
        return self.element_class(self, F.ObjByExtRep([data, one]))

    @cached_method
    def basis(self):
        r"""
        Return the basis of ``self``.

        This returns the PBW basis of ``self``, which is given by
        monomials in `\{F_{\alpha}\}`, where `\alpha` runs over all
        positive roots.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: basis = B.basis(); basis
            Lazy family (monomial(i))_{i in The Cartesian product of
             (Non negative integers, Non negative integers, Non negative integers)}
            sage: basis[1,2,1]
            F[a1]*F[a1+a2]^(2)*F[a2]
            sage: basis[1,2,4]
            F[a1]*F[a1+a2]^(2)*F[a2]^(4)
            sage: basis[1,0,4]
            F[a1]*F[a2]^(4)
        """
        I = cartesian_product([NonNegativeIntegers()]*len(self._pos_roots))
        return Family(I, self._construct_monomial, name='monomial')

    def _construct_canonical_basis_elts(self, k):
        r"""
        Construct the monomial elements of ``self`` indexed by ``k``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: B._construct_canonical_basis_elts((1,2))
            [F[a1]*F[a2]^(2), (q^2)*F[a1]*F[a2]^(2) + F[a1+a2]*F[a2]]
        """
        B = self._libgap.CanonicalBasis()
        return [self.element_class(self, v) for v in B.PBWElements(k)]

    @cached_method
    def canonical_basis_elements(self):
        r"""
        Construct the monomial elements of ``self`` indexed by ``k``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: C = B.canonical_basis_elements(); C
            Lazy family (Canonical basis(i))_{i in The Cartesian product of
             (Non negative integers, Non negative integers)}
            sage: C[2,1]
            [F[a1]^(2)*F[a2], F[a1]*F[a1+a2] + (q^2)*F[a1]^(2)*F[a2]]
            sage: C[1,2]
            [F[a1]*F[a2]^(2), (q^2)*F[a1]*F[a2]^(2) + F[a1+a2]*F[a2]]
        """
        I = cartesian_product([NonNegativeIntegers()]*len(self._cartan_type.index_set()))
        return Family(I, self._construct_canonical_basis_elts, name='Canonical basis')

    def lift(self, elt):
        r"""
        Lift ``elt`` to the ambient quantum group of ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: x = B.lift(B.an_element()); x
            1
            sage: x.parent() is Q
            True
        """
        return self._Q.element_class(self._Q, elt._libgap)

    def retract(self, elt):
        r"""
        Retract ``elt`` from the ambient quantum group to ``self``.

        EXAMPLES::

            sage: Q = QuantumGroup(['A',2])
            sage: B = Q.lower_half()
            sage: x = Q.an_element(); x
            1 + (q)*F[a1] + E[a1] + (q^2-1-q^-2 + q^-4)*[ K1 ; 2 ]
             + K1 + (-q^-1 + q^-3)*K1[ K1 ; 1 ]
            sage: B.retract(x)
            1 + (q)*F[a1]
        """
        return self.element_class(self, self._proj(elt)._libgap)

    class Element(QuaGroupModuleElement):
        """
        An element of the lower half of the quantum group.
        """
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: B = Q.lower_half()
                sage: F1, F2 = Q.F_simple()
                sage: v = B.highest_weight_vector(); v
                1
                sage: 2 * v
                (2)*1
                sage: v * (3/2)
                (3/2)*1
                sage: F1 * v
                F[a1]
                sage: F2 * (F1 * v)
                (q)*F[a1]*F[a2] + F[a1+a2]
                sage: (F1 * v) * F2
                F[a1]*F[a2]
            """
            try:
                if scalar.parent() is self.parent()._Q:
                    if self_on_left:
                        ret = self._libgap * scalar._libgap
                    else:
                        ret = scalar._libgap * self._libgap
                    return self.__class__(self.parent(), self.parent()._proj(ret)._libgap)
            except AttributeError:
                pass
            return QuaGroupModuleElement._acted_upon_(self, scalar, self_on_left)

        _lmul_ = _acted_upon_

        def _mul_(self, other):
            r"""
            Multiply ``self`` and ``other``.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: B = Q.lower_half()
                sage: F1, F2 = Q.F_simple()
                sage: v = B.highest_weight_vector()
                sage: f1, f2 = F1 * v, F2 * v
                sage: f1 * f2
                F[a1]*F[a2]
                sage: f1^2 * f2
                (q + q^-1)*F[a1]^(2)*F[a2]
                sage: f2 * f1^2 * f2
                (q + q^-1)*F[a1]*F[a1+a2]*F[a2]
                 + (q^4 + 2*q^2 + 1)*F[a1]^(2)*F[a2]^(2)
            """
            ret = self.parent()._proj(self._libgap * other._libgap)
            return self.__class__(self.parent(), ret._libgap)

        def monomial_coefficients(self, copy=True):
            r"""
            Return the dictionary of ``self`` whose keys are the basis
            indices and the values are coefficients.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: B = Q.lower_half()
                sage: x = B.retract(Q.an_element()); x
                1 + (q)*F[a1]
                sage: sorted(x.monomial_coefficients().items(), key=str)
                [((0, 0, 0), 1), ((1, 0, 0), q)]
            """
            ext_rep = self._libgap.ExtRepOfObj()
            num_pos_roots = len(self.parent()._pos_roots)
            R = self.parent().base_ring()
            d = {}
            for i in range(len(ext_rep)//2):
                exp = [0] * num_pos_roots
                mon = ext_rep[2*i].sage()
                for j in range(len(mon)//2):
                    exp[mon[2*j]-1] = mon[2*j+1]
                d[tuple(exp)] = R(str(ext_rep[2*i+1]))
            return d

        def bar(self):
            r"""
            Return the bar involution on ``self``.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: F1, F2 = Q.F_simple()
                sage: B = Q.lower_half()
                sage: x = B(Q.an_element()); x
                1 + (q)*F[a1]
                sage: x.bar()
                1 + (q^-1)*F[a1]
                sage: (F1*x).bar() == F1 * x.bar()
                True
                sage: (F2*x).bar() == F2 * x.bar()
                True

                sage: Q = QuantumGroup(['G',2])
                sage: F1, F2 = Q.F_simple()
                sage: q = Q.q()
                sage: B = Q.lower_half()
                sage: x = B(q^-2*F1*F2^2*F1)
                sage: x
                (q + q^-5)*F[a1]*F[a1+a2]*F[a2]
                 + (q^8 + q^6 + q^2 + 1)*F[a1]^(2)*F[a2]^(2)
                sage: x.bar()
                (q^5 + q^-1)*F[a1]*F[a1+a2]*F[a2]
                 + (q^12 + q^10 + q^6 + q^4)*F[a1]^(2)*F[a2]^(2)
            """
            bar = self.parent()._libgap.BarAutomorphism()
            # bar does not introduce E/K/Ki's
            return self.__class__(self.parent(), libgap.Image(bar, self._libgap))

        def tau(self):
            r"""
            Return the action of the `\tau` anti-automorphism on ``self``.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: F1, F2 = Q.F_simple()
                sage: B = Q.lower_half()
                sage: x = B(Q.an_element()); x
                1 + (q)*F[a1]
                sage: x.tau()
                1 + (q)*F[a1]
                sage: (F1*x).tau() == x.tau() * F1.tau()
                True
                sage: (F2*x).tau() == x.tau() * F2.tau()
                True

                sage: Q = QuantumGroup(['G',2])
                sage: F1, F2 = Q.F_simple()
                sage: q = Q.q()
                sage: B = Q.lower_half()
                sage: x = B(q^-2*F1*F2^2*F1)
                sage: x
                (q + q^-5)*F[a1]*F[a1+a2]*F[a2]
                 + (q^8 + q^6 + q^2 + 1)*F[a1]^(2)*F[a2]^(2)
                sage: x.tau()
                (q + q^-5)*F[a1]*F[a1+a2]*F[a2]
                 + (q^8 + q^6 + q^2 + 1)*F[a1]^(2)*F[a2]^(2)
            """
            tau = self.parent()._libgap.AntiAutomorphismTau()
            # tau does not introduce E/K/Ki's
            return self.__class__(self.parent(), libgap.Image(tau, self._libgap))

        def braid_group_action(self, braid):
            r"""
            Return the action of the braid group element ``braid``
            projected into ``self``.

            INPUT:

            - ``braid`` -- a reduced word of a braid group element

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: L = Q.lower_half()
                sage: v = L.highest_weight_vector().f_tilde([1,2,2,1]); v
                F[a1]*F[a1+a2]*F[a2]
                sage: v.braid_group_action([1])
                (-q^3-q)*F[a2]^(2)
                sage: v.braid_group_action([]) == v
                True
            """
            if not braid:
                return self
            Q = self.parent()
            QU = Q._libgap
            tau = QU.AntiAutomorphismTau()
            f = QU.IdentityMapping()
            for i in braid:
                if i < 0:
                    i = -i
                    T = QU.AutomorphismTalpha(i)
                    f *= tau * T * tau
                else:
                    f *= QU.AutomorphismTalpha(i)
            ret = libgap.Image(f, self._libgap)
            return self.__class__(Q, Q._proj(ret)._libgap)

        def _et(self, i):
            r"""
            Return the action of `\widetilde{e}_i` on ``self``.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: L = Q.lower_half()
                sage: v = L.highest_weight_vector()
                sage: v._et(1)
                0
                sage: w = v.f_tilde([1,2,1]); w
                F[a1]*F[a1+a2]
                sage: w._et(1)
                F[a1+a2]
                sage: w._et(2)
                F[a1]^(2)
                sage: L.zero().e_tilde(1)
                0
            """
            if not self:  # self == 0
                return self
            Q = self.parent()
            ret = self._libgap.Ealpha(i)
            if not ret:
                return self.parent().zero()
            return self.__class__(Q, Q._proj(ret)._libgap)

        def _ft(self, i):
            r"""
            Return the action of `\widetilde{e}_i` on ``self``.

            EXAMPLES::

                sage: Q = QuantumGroup(['A',2])
                sage: L = Q.lower_half()
                sage: v = L.highest_weight_vector()
                sage: v._ft(1)
                F[a1]
                sage: L.zero().f_tilde(1)
                0
            """
            if not self:  # self == 0
                return self
            Q = self.parent()
            ret = self._libgap.Falpha(i)
            if not ret:
                return self.parent().zero()
            return self.__class__(Q, Q._proj(ret)._libgap)


def _unpickle_generic_element(parent, data):
    """
    Used to unpickle an element of ``parent`` using ``data``.

    EXAMPLES::

        sage: Q = QuantumGroup(['D',4])
        sage: x = Q.an_element()
        sage: loads(dumps(x)) == x  # indirect doctest
        True
    """
    F = parent._libgap.FamilyObj().ElementsFamily()
    ret = []
    # We need to multiply by this to get the right type in GAP
    one = parent._libgap_base.One()
    for i in range(len(data) // 2):
        ret.append(libgap(data[2 * i]))
        ret.append(one * libgap(data[2 * i + 1].subs(q=parent._libgap_q)))
    return parent.element_class(parent, F.ObjByExtRep(ret))
