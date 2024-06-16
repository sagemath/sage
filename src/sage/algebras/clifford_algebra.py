# sage.doctest: needs sage.modules
r"""
Clifford Algebras

AUTHORS:

- Travis Scrimshaw (2013-09-06): Initial version
- Trevor K. Karn (2022-07-27): Rewrite basis indexing using FrozenBitset
"""
# ****************************************************************************
#       Copyright (C) 2013-2022 Travis Scrimshaw <tcscrims at gmail.com>
#                 (C) 2022 Trevor Karn <karnx018 at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import (richcmp_method, op_EQ, op_NE,
                                    op_LT, op_GT, op_LE, op_GE, rich_to_bool)
from sage.data_structures.bitset import Bitset, FrozenBitset

from sage.algebras.clifford_algebra_element import CliffordAlgebraElement, ExteriorAlgebraElement
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.categories.fields import Fields
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.modules.with_basis.morphism import ModuleMorphismByLinearity
from sage.categories.poor_man_map import PoorManMap
from sage.rings.integer_ring import ZZ
from sage.rings.noncommutative_ideals import Ideal_nc
from sage.modules.free_module import FreeModule, FreeModule_generic
from sage.matrix.constructor import Matrix
from sage.matrix.args import MatrixArgs
from sage.sets.family import Family
from sage.combinat.free_module import CombinatorialFreeModule
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.typeset.ascii_art import ascii_art
from sage.typeset.unicode_art import unicode_art
import unicodedata


class CliffordAlgebraIndices(UniqueRepresentation, Parent):
    r"""
    A facade parent for the indices of Clifford algebra.
    Users should not create instances of this class directly.
    """
    def __init__(self, Qdim, degree=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: idx = CliffordAlgebraIndices(7)
            sage: idx._nbits
            7
            sage: idx._cardinality
            128
            sage: i = idx.an_element(); i
            1111
            sage: type(i)
            <class 'sage.data_structures.bitset.FrozenBitset'>

            sage: idx = CliffordAlgebraIndices(7, 3)
            sage: idx._nbits
            7
            sage: idx._degree
            3
            sage: idx._cardinality
            35

            sage: idx = CliffordAlgebraIndices(7, 0)
            sage: idx._nbits
            7
            sage: idx._degree
            0
            sage: idx._cardinality
            1
        """
        self._nbits = Qdim
        if degree is None:
            self._cardinality = 2 ** Qdim
        else:
            from sage.arith.misc import binomial
            self._cardinality = binomial(Qdim, degree)
        self._degree = degree
        # the if statement here is in case Qdim is 0.
        category = FiniteEnumeratedSets().Facade()
        Parent.__init__(self, category=category, facade=True)

    def _element_constructor_(self, x):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: idx = CliffordAlgebraIndices(7)
            sage: idx([1,3,6])
            0101001
            sage: for i in range(7): print(idx(i))
            1
            01
            001
            0001
            00001
            000001
            0000001

            sage: idx = CliffordAlgebraIndices(0)
            sage: idx([])
            0
        """
        if isinstance(x, (list, tuple, set, frozenset)):
            if len(x) > self._nbits:
                raise ValueError(f"{x=} is too long")
            if not x:
                return FrozenBitset()
            return FrozenBitset(x)

        if isinstance(x, int):
            return FrozenBitset((x,))

    def __call__(self, el):
        r"""
        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: idx = CliffordAlgebraIndices(7)
            sage: idx([1,3,6])
            0101001
            sage: E = ExteriorAlgebra(QQ, 7)
            sage: B = E.basis()
        """
        if not isinstance(el, Element):
            return self._element_constructor_(el)
        else:
            return Parent.__call__(self, el)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: idx = CliffordAlgebraIndices(7)
            sage: idx.cardinality() == 2^7
            True
            sage: len(idx) == 2^7
            True

            sage: idx = CliffordAlgebraIndices(7, 3)
            sage: idx.cardinality() == binomial(7, 3)
            True
            sage: len(idx) == binomial(7, 3)
            True
        """
        return self._cardinality

    __len__ = cardinality

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: CliffordAlgebraIndices(7)
            Subsets of {0,1,...,6}
            sage: CliffordAlgebraIndices(0)
            Subsets of {}
            sage: CliffordAlgebraIndices(1)
            Subsets of {0}
            sage: CliffordAlgebraIndices(2)
            Subsets of {0,1}
            sage: CliffordAlgebraIndices(5, 3)
            Subsets of {0,1,...,4} of size 3
        """
        if self._degree is not None:
            extra = f" of size {self._degree}"
        else:
            extra = ""
        if self._nbits == 0:
            return "Subsets of {}" + extra
        if self._nbits == 1:
            return "Subsets of {0}" + extra
        if self._nbits == 2:
            return "Subsets of {0,1}" + extra
        return f"Subsets of {{0,1,...,{self._nbits-1}}}" + extra

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: latex(CliffordAlgebraIndices(7))
            \mathcal{P}(\{0,1,\ldots,6\})
            sage: latex(CliffordAlgebraIndices(0))
            \mathcal{P}(\emptyset)
            sage: latex(CliffordAlgebraIndices(1))
            \mathcal{P}(\{0\})
            sage: latex(CliffordAlgebraIndices(2))
            \mathcal{P}(\{0,1\})
            sage: latex(CliffordAlgebraIndices(2, 1))
            \mathcal{P}(\{0,1\}, 1)
        """
        if self._degree is not None:
            extra = f", {self._degree}"
        else:
            extra = ""
        if self._nbits == 0:
            return f"\\mathcal{{P}}(\\emptyset{extra})"
        if self._nbits == 1:
            return f"\\mathcal{{P}}(\\{{0\\}}{extra})"
        if self._nbits == 2:
            return f"\\mathcal{{P}}(\\{{0,1\\}}{extra})"
        return f"\\mathcal{{P}}(\\{{0,1,\\ldots,{self._nbits-1}\\}}{extra})"

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: idx = CliffordAlgebraIndices(3)
            sage: for i in idx:
            ....:     print(i)
            0
            1
            01
            001
            11
            101
            011
            111

            sage: idx = CliffordAlgebraIndices(5, 3)
            sage: list(idx)
            [111, 1101, 11001, 1011, 10101, 10011, 0111, 01101, 01011, 00111]

            sage: idx = CliffordAlgebraIndices(7, 0)
            sage: list(idx)
            [0]
        """
        import itertools
        n = self._nbits
        if self._degree is not None:
            if self._degree == 0:  # special corner case
                yield FrozenBitset()
                return
            for C in itertools.combinations(range(n), self._degree):
                yield FrozenBitset(C)
            return

        yield FrozenBitset()
        k = 1
        while k <= n:
            for C in itertools.combinations(range(n), k):
                yield FrozenBitset(C)
            k += 1

    def __contains__(self, elt):
        r"""
        Check containment of ``elt`` in ``self``.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: idx = CliffordAlgebraIndices(3)
            sage: int(8) in idx  # representing the set {4}
            False
            sage: int(5) in idx  # representing the set {1,3}
            True
            sage: FrozenBitset('1') in idx
            True
            sage: FrozenBitset('000001') in idx
            False

            sage: idx = CliffordAlgebraIndices(6, 3)
            sage: FrozenBitset('01011') in idx
            True
            sage: FrozenBitset('00011') in idx
            False
            sage: int(7) in idx
            True
            sage: int(8) in idx
            False

            sage: idx = CliffordAlgebraIndices(7, 0)
            sage: FrozenBitset() in idx
            True
            sage: FrozenBitset('01') in idx
            False
            sage: int(0) in idx
            True
            sage: int(5) in idx
            False
        """
        if isinstance(elt, int):
            if self._degree is not None and sum(ZZ(elt).bits()) != self._degree:
                return False
            return elt < self._cardinality and elt >= 0
        if not isinstance(elt, FrozenBitset):
            return False
        if self._degree is not None and len(elt) != self._degree:
            return False
        return elt.capacity() <= self._nbits

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import CliffordAlgebraIndices
            sage: idx = CliffordAlgebraIndices(0)
            sage: idx._an_element_()
            0
            sage: idx = CliffordAlgebraIndices(1)
            sage: idx._an_element_()
            1
            sage: idx = CliffordAlgebraIndices(2)
            sage: idx._an_element_()
            01
            sage: idx = CliffordAlgebraIndices(3)
            sage: idx._an_element_()
            11
            sage: idx = CliffordAlgebraIndices(5, 3)
            sage: idx._an_element_()
            111
            sage: idx = CliffordAlgebraIndices(7, 0)
            sage: idx._an_element_()
            0
        """
        if not self._nbits:
            return FrozenBitset()

        if self._degree is not None:
            if self._degree == 0:  # special corner case
                return FrozenBitset()
            return FrozenBitset(range(self._degree))

        from sage.combinat.subset import SubsetsSorted
        X = SubsetsSorted(range(self._nbits))
        return FrozenBitset(X.an_element())


class CliffordAlgebra(CombinatorialFreeModule):
    r"""
    The Clifford algebra of a quadratic form.

    Let `Q : V \to \mathbf{k}` denote a quadratic form on a vector space `V`
    over a field `\mathbf{k}`. The Clifford algebra `Cl(V, Q)` is defined as
    `T(V) / I_Q` where `T(V)` is the tensor algebra of `V` and `I_Q` is the
    two-sided ideal generated by all elements of the form `v \otimes v - Q(v)`
    for all `v \in V`.

    We abuse notation to denote the projection of a pure tensor
    `x_1 \otimes x_2 \otimes \cdots \otimes x_m \in T(V)` onto
    `T(V) / I_Q = Cl(V, Q)` by `x_1 \wedge x_2 \wedge \cdots \wedge x_m`.
    This is motivated by the fact that `Cl(V, Q)` is the exterior algebra
    `\wedge V` when `Q = 0` (one can also think of a Clifford algebra as
    a quantization of the exterior algebra). See :class:`ExteriorAlgebra`
    for the concept of an exterior algebra.

    From the definition, a basis of `Cl(V, Q)` is given by monomials of
    the form

    .. MATH::

        \{ e_{i_1} \wedge \cdots \wedge e_{i_k} \mid 1 \leq i_1 < \cdots <
        i_k \leq n \},

    where `n = \dim(V)` and where `\{ e_1, e_2, \cdots, e_n \}` is any
    fixed basis of `V`. Hence

    .. MATH::

        \dim(Cl(V, Q)) = \sum_{k=0}^n \binom{n}{k} = 2^n.

    .. NOTE::

        The algebra `Cl(V, Q)` is a `\ZZ / 2\ZZ`-graded algebra, but not
        (in general) `\ZZ`-graded (in a reasonable way).

    This construction satisfies the following universal property. Let
    `i : V \to Cl(V, Q)` denote the natural inclusion (which is an
    embedding). Then for every associative `\mathbf{k}`-algebra `A`
    and any `\mathbf{k}`-linear map `j : V \to A` satisfying

    .. MATH::

        j(v)^2 = Q(v) \cdot 1_A

    for all `v \in V`, there exists a unique `\mathbf{k}`-algebra
    homomorphism `f : Cl(V, Q) \to A` such that `f \circ i = j`.
    This property determines the Clifford algebra uniquely up to
    canonical isomorphism. The inclusion `i` is commonly used to
    identify `V` with a vector subspace of `Cl(V)`.

    The Clifford algebra `Cl(V, Q)` is a `\ZZ_2`-graded algebra
    (where `\ZZ_2 = \ZZ / 2 \ZZ`); this grading is determined by
    placing all elements of `V` in degree `1`. It is also an
    `\NN`-filtered algebra, with the filtration too being defined
    by placing all elements of `V` in degree `1`. The :meth:`degree` gives
    the `\NN`-*filtration* degree, and to get the super degree use instead
    :meth:`~sage.categories.super_modules.SuperModules.ElementMethods.is_even_odd`.

    The Clifford algebra also can be considered as a covariant functor
    from the category of vector spaces equipped with quadratic forms
    to the category of algebras. In fact, if `(V, Q)` and `(W, R)`
    are two vector spaces endowed with quadratic forms, and if
    `g : W \to V` is a linear map preserving the quadratic form,
    then we can define an algebra morphism
    `Cl(g) : Cl(W, R) \to Cl(V, Q)` by requiring that it send every
    `w \in W` to `g(w) \in V`. Since the quadratic form `R` on `W`
    is uniquely determined by the quadratic form `Q` on `V` (due to
    the assumption that `g` preserves the quadratic form), this fact
    can be rewritten as follows: If `(V, Q)` is a vector space with a
    quadratic form, and `W` is another vector space, and
    `\phi : W \to V` is any linear map, then we obtain an algebra
    morphism `Cl(\phi) : Cl(W, \phi(Q)) \to Cl(V, Q)` where
    `\phi(Q) = \phi^T \cdot Q \cdot \phi` (we consider `\phi` as a
    matrix) is the quadratic form `Q` pulled back to `W`. In fact, the
    map `\phi` preserves the quadratic form because of

    .. MATH::

        \phi(Q)(x) = x^T \cdot \phi^T \cdot Q \cdot \phi \cdot x
        = (\phi \cdot x)^T \cdot Q \cdot (\phi \cdot x) = Q(\phi(x)).

    Hence we have `\phi(w)^2 = Q(\phi(w)) = \phi(Q)(w)` for all `w \in W`.

    REFERENCES:

    - :wikipedia:`Clifford_algebra`

    INPUT:

    - ``Q`` -- a quadratic form
    - ``names`` -- (default: ``'e'``) the generator names

    EXAMPLES:

    To create a Clifford algebra, all one needs to do is specify a
    quadratic form::

        sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
        sage: Cl = CliffordAlgebra(Q)
        sage: Cl
        The Clifford algebra of the Quadratic form in 3 variables
         over Integer Ring with coefficients:
        [ 1 2 3 ]
        [ * 4 5 ]
        [ * * 6 ]

    We can also explicitly name the generators. In this example, the
    Clifford algebra we construct is an exterior algebra (since we
    choose the quadratic form to be zero)::

        sage: Q = QuadraticForm(ZZ, 4, [0]*10)
        sage: Cl.<a,b,c,d> = CliffordAlgebra(Q)
        sage: a*d
        a*d
        sage: d*c*b*a + a + 4*b*c
        a*b*c*d + 4*b*c + a
    """
    @staticmethod
    def __classcall_private__(cls, Q, names=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl1.<e0,e1,e2> = CliffordAlgebra(Q)
            sage: Cl2 = CliffordAlgebra(Q)
            sage: Cl3 = CliffordAlgebra(Q, ['e0','e1','e2'])
            sage: Cl1 is Cl2 and Cl2 is Cl3
            True
        """
        if not isinstance(Q, QuadraticForm):
            raise ValueError("{} is not a quadratic form".format(Q))
        if names is None:
            names = 'e'
        names = tuple(names)
        if len(names) != Q.dim():
            if len(names) == 1:
                names = tuple('{}{}'.format(names[0], i) for i in range(Q.dim()))
            else:
                raise ValueError("the number of variables does not match the number of generators")
        return super().__classcall__(cls, Q, names)

    def __init__(self, Q, names, category=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl = CliffordAlgebra(Q)
            sage: Cl.category()
            Category of finite dimensional super algebras with basis over
             (Dedekind domains and euclidean domains
              and noetherian rings
              and infinite enumerated sets and metric spaces)
            sage: TestSuite(Cl).run()

        TESTS:

        We check that the basis elements are indeed indexed by
        *strictly increasing* tuples::

            sage: Q = QuadraticForm(ZZ, 9)
            sage: Cl = CliffordAlgebra(Q)
            sage: ba = Cl.basis().keys()
            sage: all(FrozenBitset(format(i,'b')[::-1]) in ba for i in range(2**9))
            True
        """
        self._quadratic_form = Q
        R = Q.base_ring()
        category = AlgebrasWithBasis(R.category()).Super().Filtered().FiniteDimensional().or_subcategory(category)
        indices = CliffordAlgebraIndices(Q.dim())
        CombinatorialFreeModule.__init__(self, R, indices, category=category, sorting_key=tuple)
        self._assign_names(names)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: CliffordAlgebra(Q)
            The Clifford algebra of the Quadratic form in 3 variables
             over Integer Ring with coefficients:
            [ 1 2 3 ]
            [ * 4 5 ]
            [ * * 6 ]
        """
        return "The Clifford algebra of the {}".format(self._quadratic_form)

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl._repr_term((0,2))
            'x*z'
            sage: Cl._repr_term(FrozenBitset('101'))
            'x*z'
            sage: Cl._repr_term(())
            '1'
            sage: Cl._repr_term((1,))
            'y'
        """
        if not m:
            return '1'
        term = ''
        for i in m:
            if term:
                term += '*'
            term += self.variable_names()[i]
        return term

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the basis element indexed
        by ``m``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl._latex_term((0,2))
            ' x z'
        """
        if not m:
            return '1'
        term = ''
        for i in m:
            term += ' ' + self.latex_variable_names()[i]
        return term

    def _coerce_map_from_(self, V):
        """
        Return if there is a coerce map from ``V`` into ``self``.

        The things which coerce into ``self`` are:

        - Clifford algebras with the same generator names and an equal
          quadratic form over a ring which coerces into the base
          ring of ``self``.
        - The underlying free module of ``self``.
        - The base ring of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Qp = QuadraticForm(QQ, 3, [1,2,3,4,5,6])
            sage: Cl = CliffordAlgebra(Q)
            sage: Clp = CliffordAlgebra(Qp)
            sage: Cl.has_coerce_map_from(Clp)
            False
            sage: Clp.has_coerce_map_from(Cl)
            True

        Check that we preserve the multiplicative structure::

            sage: all(Clp(b)*Clp(b) == Clp(b*b) for b in Cl.basis())
            True

        Check from the underlying free module::

            sage: M = ZZ^3
            sage: Mp = QQ^3
            sage: Cl.has_coerce_map_from(M)
            True
            sage: Cl.has_coerce_map_from(Mp)
            False
            sage: Clp.has_coerce_map_from(M)
            True
            sage: Clp.has_coerce_map_from(Mp)
            True

        Names matter::

            sage: Cln = CliffordAlgebra(Q, names=['x','y','z'])
            sage: Cln.has_coerce_map_from(Cl)
            False
            sage: Cl.has_coerce_map_from(Cln)
            False

        Non-injective homomorphisms of base rings don't cause zero
        values in the coordinate dictionary (this had to be manually
        ensured)::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Qp = QuadraticForm(Integers(3), 3, [1,2,3,4,5,6])
            sage: Cl = CliffordAlgebra(Q)
            sage: Clp = CliffordAlgebra(Qp)
            sage: a = Cl.basis()[(1,2)]
            sage: a
            e1*e2
            sage: Clp(a) # so far so good
            e1*e2
            sage: Clp(3*a) # but now
            0
            sage: Clp(3*a) == 0
            True
            sage: b = Cl.basis()[(0,2)]
            sage: Clp(3*a-4*b)
            2*e0*e2
        """
        if isinstance(V, CliffordAlgebra):
            Q = self._quadratic_form
            try:
                return (V.variable_names() == self.variable_names() and
                        V._quadratic_form.change_ring(self.base_ring()) == Q)
            except (TypeError, AttributeError):
                return False

        if self.free_module().has_coerce_map_from(V):
            return True

        return super()._coerce_map_from_(V)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Qp = QuadraticForm(QQ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Clp = CliffordAlgebra(Qp, names=['x','y','z'])
            sage: M = ZZ^3
            sage: Mp = QQ^3
            sage: Cl(2/3)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x=2/3 an element of self
            sage: Clp(2/3)
            2/3
            sage: Clp(x)
            x
            sage: M = ZZ^3
            sage: Clp( M((1,-3,2)) )
            x - 3*y + 2*z

        Zero coordinates are handled appropriately::

            sage: Q3 = QuadraticForm(Integers(3), 3, [1,2,3,4,5,6])
            sage: Cl3 = CliffordAlgebra(Q3, names='xyz')  # different syntax for a change
            sage: Cl3( M((1,-3,2)) )
            x + 2*z
        """
        # This is the natural lift morphism of the underlying free module
        if x in self.free_module():
            R = self.base_ring()
            if x.parent().base_ring() is R:
                return self.element_class(self, {FrozenBitset((i,)): c for i, c in x.items()})
            # if the base ring is different, attempt to coerce it into R
            return self.element_class(self, {FrozenBitset((i,)): R(c) for i, c in x.items() if R(c) != R.zero()})

        if (isinstance(x, CliffordAlgebraElement)
                and self.has_coerce_map_from(x.parent())):
            R = self.base_ring()
            return self.element_class(self, {i: R(c) for i, c in x if R(c) != R.zero()})

        if isinstance(x, tuple):
            R = self.base_ring()
            return self.element_class(self, {FrozenBitset((i,)): R.one() for i in x})

        try:
            return super()._element_constructor_(x)
        except TypeError:
            raise TypeError(f'do not know how to make {x=} an element of self')

    def _basis_index_function(self, x):
        """
        Given an integer indexing the basis, return the correct
        bitset.

        For backwards compatibility, tuples are also accepted.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl = CliffordAlgebra(Q)
            sage: Cl._basis_index_function(7)
            111
            sage: Cl._basis_index_function(5)
            101
            sage: Cl._basis_index_function(4)
            001

            sage: Cl._basis_index_function((0, 1, 2))
            111
            sage: Cl._basis_index_function((0, 2))
            101
            sage: Cl._basis_index_function((2,))
            001
        """
        Q = self._quadratic_form
        format_style = f"0{Q.dim()}b"

        # if the input is a tuple, assume that it has
        # entries in {0, ..., 2**Q.dim()-1}
        if isinstance(x, tuple):
            return FrozenBitset(x, capacity=Q.dim())

        # slice the output of format in order to make conventions
        # of format and FrozenBitset agree.
        return FrozenBitset(format(x, format_style)[::-1], capacity=Q.dim())

    def gen(self, i):
        """
        Return the ``i``-th standard generator of the algebra ``self``.

        This is the ``i``-th basis vector of the vector space on which
        the quadratic form defining ``self`` is defined, regarded as an
        element of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: [Cl.gen(i) for i in range(3)]
            [x, y, z]
        """
        return self._from_dict({FrozenBitset((i,)): self.base_ring().one()}, remove_zeros=False)

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.algebra_generators()
            Finite family {'x': x, 'y': y, 'z': z}
        """
        d = {x: self.gen(i) for i, x in enumerate(self.variable_names())}
        return Family(self.variable_names(), lambda x: d[x])

    def gens(self):
        r"""
        Return the generators of ``self`` (as an algebra).

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.gens()
            (x, y, z)
        """
        return tuple(self.algebra_generators())

    @cached_method
    def ngens(self):
        """
        Return the number of algebra generators of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.ngens()
            3
        """
        return self._quadratic_form.dim()

    @cached_method
    def one_basis(self):
        """
        Return the basis index of the element ``1``. The element ``1``
        is indexed by the emptyset, which is represented by the
        :class:`sage.data_structures.bitset.Bitset` ``0``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.one_basis()
            0
        """
        return FrozenBitset()

    def is_commutative(self):
        """
        Check if ``self`` is a commutative algebra.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.is_commutative()
            False
        """
        return self._quadratic_form.dim() < 2

    def quadratic_form(self):
        """
        Return the quadratic form of ``self``.

        This is the quadratic form used to define ``self``. The
        quadratic form on ``self`` is yet to be implemented.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.quadratic_form()
            Quadratic form in 3 variables over Integer Ring with coefficients:
            [ 1 2 3 ]
            [ * 4 5 ]
            [ * * 6 ]
        """
        return self._quadratic_form

    def degree_on_basis(self, m):
        r"""
        Return the degree of the monomial indexed by ``m``.

        We are considering the Clifford algebra to be `\NN`-filtered,
        and the degree of the monomial ``m`` is the length of ``m``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.degree_on_basis((0,))
            1
            sage: Cl.degree_on_basis((0,1))
            2
        """
        return ZZ(len(m))

    def graded_algebra(self):
        """
        Return the associated graded algebra of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.graded_algebra()
            The exterior algebra of rank 3 over Integer Ring
        """
        return ExteriorAlgebra(self.base_ring(), self.variable_names())

    @cached_method
    def free_module(self):
        """
        Return the underlying free module `V` of ``self``.

        This is the free module on which the quadratic form that was
        used to construct ``self`` is defined.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.free_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        return FreeModule(self.base_ring(), self._quadratic_form.dim())

    def dimension(self):
        """
        Return the rank of ``self`` as a free module.

        Let `V` be a free `R`-module of rank `n`; then, `Cl(V, Q)` is a
        free `R`-module of rank `2^n`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.dimension()
            8
        """
        return ZZ(2)**self._quadratic_form.dim()

    def pseudoscalar(self):
        r"""
        Return the unit pseudoscalar of ``self``.

        Given the basis `e_1, e_2, \ldots, e_n` of the underlying
        `R`-module, the unit pseudoscalar is defined as
        `e_1 \cdot e_2 \cdots e_n`.

        This depends on the choice of basis.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.pseudoscalar()
            x*y*z

            sage: Q = QuadraticForm(ZZ, 0, [])
            sage: Cl = CliffordAlgebra(Q)
            sage: Cl.pseudoscalar()
            1

        REFERENCES:

        - :wikipedia:`Classification_of_Clifford_algebras#Unit_pseudoscalar`
        """
        d = self._quadratic_form.dim()
        return self.element_class(self, {tuple(range(d)): self.base_ring().one()})

    def lift_module_morphism(self, m, names=None):
        r"""
        Lift the matrix ``m`` to an algebra morphism of Clifford algebras.

        Given a linear map `m : W \to V` (here represented by a matrix
        acting on column vectors), this method returns the algebra
        morphism `Cl(m) : Cl(W, m(Q)) \to Cl(V, Q)`, where `Cl(V, Q)`
        is the Clifford algebra ``self`` and where `m(Q)` is the pullback
        of the quadratic form `Q` to `W`. See the documentation
        of :class:`CliffordAlgebra` for how this pullback and the
        morphism `Cl(m)` are defined.

        .. NOTE::

            This is a map into ``self``.

        INPUT:

        - ``m`` -- a matrix
        - ``names`` -- (default: ``'e'``) the names of the generators of the
          Clifford algebra of the domain of (the map represented by) ``m``

        OUTPUT: the algebra morphism `Cl(m)` from `Cl(W, m(Q))` to ``self``

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: m = matrix([[1,-1,-1],[0,1,-1],[1,1,1]])
            sage: phi = Cl.lift_module_morphism(m, 'abc')
            sage: phi
            Generic morphism:
              From: The Clifford algebra of the Quadratic form in 3 variables over Integer Ring with coefficients:
            [ 10 17 3 ]
            [ * 11 0 ]
            [ * * 5 ]
              To:   The Clifford algebra of the Quadratic form in 3 variables over Integer Ring with coefficients:
            [ 1 2 3 ]
            [ * 4 5 ]
            [ * * 6 ]
            sage: a,b,c = phi.domain().gens()
            sage: phi(a)
            x + z
            sage: phi(b)
            -x + y + z
            sage: phi(c)
            -x - y + z
            sage: phi(a + 3*b)
            -2*x + 3*y + 4*z
            sage: phi(a) + 3*phi(b)
            -2*x + 3*y + 4*z
            sage: phi(a*b)
            x*y + 2*x*z - y*z + 7
            sage: phi(b*a)
            -x*y - 2*x*z + y*z + 10
            sage: phi(a*b + c)
            x*y + 2*x*z - y*z - x - y + z + 7
            sage: phi(a*b) + phi(c)
            x*y + 2*x*z - y*z - x - y + z + 7

        We check that the map is an algebra morphism::

            sage: phi(a)*phi(b)
            x*y + 2*x*z - y*z + 7
            sage: phi(a*b)
            x*y + 2*x*z - y*z + 7
            sage: phi(a*a)
            10
            sage: phi(a)*phi(a)
            10
            sage: phi(b*a)
            -x*y - 2*x*z + y*z + 10
            sage: phi(b) * phi(a)
            -x*y - 2*x*z + y*z + 10
            sage: phi((a + b)*(a + c)) == phi(a + b) * phi(a + c)
            True

        We can also lift arbitrary linear maps::

            sage: m = matrix([[1,1],[0,1],[1,1]])
            sage: phi = Cl.lift_module_morphism(m, 'ab')
            sage: a,b = phi.domain().gens()
            sage: phi(a)
            x + z
            sage: phi(b)
            x + y + z
            sage: phi(a*b)
            x*y - y*z + 15
            sage: phi(a)*phi(b)
            x*y - y*z + 15
            sage: phi(b*a)
            -x*y + y*z + 12
            sage: phi(b)*phi(a)
            -x*y + y*z + 12

            sage: m = matrix([[1,1,1,2], [0,1,1,1], [0,1,1,1]])
            sage: phi = Cl.lift_module_morphism(m, 'abcd')
            sage: a,b,c,d = phi.domain().gens()
            sage: phi(a)
            x
            sage: phi(b)
            x + y + z
            sage: phi(c)
            x + y + z
            sage: phi(d)
            2*x + y + z
            sage: phi(a*b*c + d*a)
            -x*y - x*z + 21*x + 7
            sage: phi(a*b*c*d)
            21*x*y + 21*x*z + 42

        TESTS:

        Check that the resulting morphism knows it is for
        finite-dimensional algebras (:issue:`25339`)::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: m = matrix([[1,-1,-1],[0,1,-1],[1,1,1]])
            sage: phi = Cl.lift_module_morphism(m, 'abc')
            sage: phi.category_for()
            Category of finite dimensional super algebras with basis over
             (Dedekind domains and euclidean domains
              and noetherian rings
              and infinite enumerated sets and metric spaces)
            sage: phi.matrix()
            [  1   0   0   0   7  -3  -7   0]
            [  0   1  -1  -1   0   0   0 -17]
            [  0   0   1  -1   0   0   0  -4]
            [  0   1   1   1   0   0   0   3]
            [  0   0   0   0   1  -1   2   0]
            [  0   0   0   0   2   2   0   0]
            [  0   0   0   0  -1   1   2   0]
            [  0   0   0   0   0   0   0   4]
        """
        Q = self._quadratic_form(m)
        # If R is a quadratic form and m is a matrix, then R(m) returns
        # the quadratic form m^t R m.

        if Q == self._quadratic_form and names is None:
            Cl = self
        else:
            Cl = CliffordAlgebra(Q, names)

        n = self._quadratic_form.dim()
        f = lambda x: self.prod(self._from_dict({FrozenBitset((j, )): m[j, i] for j in range(n)},
                                remove_zeros=True) for i in x)
        cat = AlgebrasWithBasis(self.category().base_ring()).Super().FiniteDimensional()
        return Cl.module_morphism(on_basis=f, codomain=self, category=cat)

    def lift_isometry(self, m, names=None):
        r"""
        Lift an invertible isometry ``m`` of the quadratic form of
        ``self`` to a Clifford algebra morphism.

        Given an invertible linear map `m : V \to W` (here represented by
        a matrix acting on column vectors), this method returns the
        algebra morphism `Cl(m)` from `Cl(V, Q)` to `Cl(W, m^{-1}(Q))`,
        where `Cl(V, Q)` is the Clifford algebra ``self`` and where
        `m^{-1}(Q)` is the pullback of the quadratic form `Q` to `W` along
        the inverse map `m^{-1} : W \to V`. See the documentation of
        :class:`CliffordAlgebra` for how this pullback and the morphism
        `Cl(m)` are defined.

        INPUT:

        - ``m`` -- an isometry of the quadratic form of ``self``
        - ``names`` -- (default: ``'e'``) the names of the generators of
          the Clifford algebra of the codomain of (the map represented by)
          ``m``

        OUTPUT: the algebra morphism `Cl(m)` from ``self`` to `Cl(W, m^{-1}(Q))`

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: m = matrix([[1,1,2],[0,1,1],[0,0,1]])
            sage: phi = Cl.lift_isometry(m, 'abc')
            sage: phi(x)
            a
            sage: phi(y)
            a + b
            sage: phi(x*y)
            a*b + 1
            sage: phi(x) * phi(y)
            a*b + 1
            sage: phi(z*y)
            a*b - a*c - b*c
            sage: phi(z) * phi(y)
            a*b - a*c - b*c
            sage: phi(x + z) * phi(y + z) == phi((x + z) * (y + z))
            True

        TESTS:

        Check that the resulting morphism knows it is for
        finite-dimensional algebras (:issue:`25339`)::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: m = matrix([[1,1,2],[0,1,1],[0,0,1]])
            sage: phi = Cl.lift_isometry(m, 'abc')
            sage: phi.category_for()
            Category of finite dimensional super algebras with basis over
             (Dedekind domains and euclidean domains
              and noetherian rings
              and infinite enumerated sets and metric spaces)
            sage: phi.matrix()
            [ 1  0  0  0  1  2  5  0]
            [ 0  1  1  2  0  0  0  5]
            [ 0  0  1  1  0  0  0 -1]
            [ 0  0  0  1  0  0  0  1]
            [ 0  0  0  0  1  1 -1  0]
            [ 0  0  0  0  0  1  1  0]
            [ 0  0  0  0  0  0  1  0]
            [ 0  0  0  0  0  0  0  1]
        """
        MS = m.parent()
        if not m.is_invertible():
            raise ValueError('{} is not invertible')
        Q = self._quadratic_form(MS(m.inverse()))

        if Q == self._quadratic_form and names is None:
            Cl = self
        else:
            if names is None:
                names = 'e'
            Cl = CliffordAlgebra(Q, names)

        n = Q.dim()

        f = lambda x: Cl.prod(Cl._from_dict({FrozenBitset((j, )): m[j, i] for j in range(n)},
                              remove_zeros=True) for i in x)
        cat = AlgebrasWithBasis(self.category().base_ring()).Super().FiniteDimensional()
        return self.module_morphism(on_basis=f, codomain=Cl, category=cat)

    # This is a general method for finite dimensional algebras with bases
    #   and should be moved to the corresponding category once there is
    #   a category level method for getting the indexing set of the basis;
    #   similar to #15289 but on a category level.
    @cached_method
    def center_basis(self):
        """
        Return a list of elements which correspond to a basis for the center
        of ``self``.

        This assumes that the ground ring can be used to compute the
        kernel of a matrix.

        .. SEEALSO::

            :meth:`supercenter_basis`,
            http://math.stackexchange.com/questions/129183/center-of-clifford-algebra-depending-on-the-parity-of-dim-v

        .. TODO::

            Deprecate this in favor of a method called `center()` once
            subalgebras are properly implemented in Sage.

        EXAMPLES::

            sage: Q = QuadraticForm(QQ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Z = Cl.center_basis(); Z
            (1, -2/5*x*y*z + x - 3/5*y + 2/5*z)
            sage: all(z*b - b*z == 0 for z in Z for b in Cl.basis())
            True

            sage: Q = QuadraticForm(QQ, 3, [1,-2,-3, 4, 2, 1])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Z = Cl.center_basis(); Z
            (1, -x*y*z + x + 3/2*y - z)
            sage: all(z*b - b*z == 0 for z in Z for b in Cl.basis())
            True

            sage: Q = QuadraticForm(QQ, 2, [1,-2,-3])
            sage: Cl.<x,y> = CliffordAlgebra(Q)
            sage: Cl.center_basis()
            (1,)

            sage: Q = QuadraticForm(QQ, 2, [-1,1,-3])
            sage: Cl.<x,y> = CliffordAlgebra(Q)
            sage: Cl.center_basis()
            (1,)

        A degenerate case::

            sage: Q = QuadraticForm(QQ, 3, [4,4,-4,1,-2,1])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.center_basis()
            (1, x*y*z + x - 2*y - 2*z, x*y + x*z - 2*y*z)

        The most degenerate case (the exterior algebra)::

            sage: Q = QuadraticForm(QQ, 3)
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.center_basis()
            (1, x*y, x*z, y*z, x*y*z)
        """
        R = self.base_ring()
        B = self.basis()
        K = list(B.keys())
        k = len(K)
        d = {}
        for a, i in enumerate(K):
            Bi = B[i]
            for b, j in enumerate(K):
                Bj = B[j]
                for m, c in (Bi*Bj - Bj*Bi):
                    d[(a, K.index(m)+k*b)] = c
        m = Matrix(R, d, nrows=k, ncols=k*k, sparse=True)
        from_vector = lambda x: self.sum_of_terms(((K[i], c) for i, c in x.items()),
                                                  distinct=True)
        return tuple(map(from_vector, m.kernel().basis()))

    # Same as center except for superalgebras
    @cached_method
    def supercenter_basis(self):
        """
        Return a list of elements which correspond to a basis for the
        supercenter of ``self``.

        This assumes that the ground ring can be used to compute the
        kernel of a matrix.

        .. SEEALSO::

            :meth:`center_basis`,
            http://math.stackexchange.com/questions/129183/center-of-clifford-algebra-depending-on-the-parity-of-dim-v

        .. TODO::

            Deprecate this in favor of a method called `supercenter()` once
            subalgebras are properly implemented in Sage.

        EXAMPLES::

            sage: Q = QuadraticForm(QQ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: SZ = Cl.supercenter_basis(); SZ
            (1,)
            sage: all(z.supercommutator(b) == 0 for z in SZ for b in Cl.basis())
            True

            sage: Q = QuadraticForm(QQ, 3, [1,-2,-3, 4, 2, 1])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.supercenter_basis()
            (1,)

            sage: Q = QuadraticForm(QQ, 2, [1,-2,-3])
            sage: Cl.<x,y> = CliffordAlgebra(Q)
            sage: Cl.supercenter_basis()
            (1,)

            sage: Q = QuadraticForm(QQ, 2, [-1,1,-3])
            sage: Cl.<x,y> = CliffordAlgebra(Q)
            sage: Cl.supercenter_basis()
            (1,)

        Singular vectors of a quadratic form generate in the supercenter::

            sage: Q = QuadraticForm(QQ, 3, [1/2,-2,4,256/249,3,-185/8])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.supercenter_basis()
            (1, x + 249/322*y + 22/161*z)

            sage: Q = QuadraticForm(QQ, 3, [4,4,-4,1,-2,1])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.supercenter_basis()
            (1, x + 2*z, y + z, x*y + x*z - 2*y*z)

        The most degenerate case::

            sage: Q = QuadraticForm(QQ, 3)
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.supercenter_basis()
            (1, x, y, z, x*y, x*z, y*z, x*y*z)
        """
        R = self.base_ring()
        B = self.basis()
        K = list(B.keys())
        k = len(K)
        d = {}
        for a, i in enumerate(K):
            Bi = B[i]
            for b, j in enumerate(K):
                Bj = B[j]
                if len(i) % 2 and len(j) % 2:
                    supercommutator = Bi * Bj + Bj * Bi
                else:
                    supercommutator = Bi * Bj - Bj * Bi
                for m, c in supercommutator:
                    d[(a, K.index(m) + k * b)] = c
        m = Matrix(R, d, nrows=k, ncols=k * k, sparse=True)
        from_vector = lambda x: self.sum_of_terms(((K[i], c) for i, c in x.items()),
                                                  distinct=True)
        return tuple(map(from_vector, m.kernel().basis()))

    Element = CliffordAlgebraElement


class ExteriorAlgebra(CliffordAlgebra):
    r"""
    An exterior algebra of a free module over a commutative ring.

    Let `V` be a module over a commutative ring `R`. The exterior algebra
    (or Grassmann algebra) `\Lambda(V)` of `V` is defined as the quotient
    of the tensor algebra `T(V)` of `V` modulo the two-sided ideal
    generated by all tensors of the form `x \otimes x` with `x \in V`. The
    multiplication on `\Lambda(V)` is denoted by `\wedge` (so
    `v_1 \wedge v_2 \wedge \cdots \wedge v_n` is the projection of
    `v_1 \otimes v_2 \otimes \cdots \otimes v_n` onto `\Lambda(V)`) and
    called the "exterior product" or "wedge product".

    If `V` is a rank-`n` free `R`-module with a basis
    `\{e_1, \ldots, e_n\}`, then `\Lambda(V)` is the `R`-algebra
    noncommutatively generated by the `n` generators `e_1, \ldots, e_n`
    subject to the relations `e_i^2 = 0` for all `i`, and
    `e_i e_j = - e_j e_i` for all `i < j`. As an `R`-module,
    `\Lambda(V)` then has a basis `(\bigwedge_{i \in I} e_i)` with `I`
    ranging over the subsets of `\{1, 2, \ldots, n\}` (where
    `\bigwedge_{i \in I} e_i` is the wedge product of `e_i` for `i`
    running through all elements of `I` from smallest to largest), and
    hence is free of rank `2^n`.

    The exterior algebra of an `R`-module `V` can also be realized
    as the Clifford algebra of `V` for the quadratic form `Q` given by
    `Q(v) = 0` for all vectors `v \in V`. See :class:`CliffordAlgebra`
    for the notion of a Clifford algebra.

    The exterior algebra of an `R`-module `V` is a connected `\ZZ`-graded
    Hopf superalgebra. It is commutative in the super sense (i.e., the
    odd elements anticommute and square to `0`).

    This class implements the exterior algebra `\Lambda(R^n)` for
    `n` a nonnegative integer.

    INPUT:

    - ``R`` -- the base ring, *or* the free module whose exterior algebra
      is to be computed

    - ``names`` -- list of strings to name the generators of the
      exterior algebra; this list can either have one entry only (in which
      case the generators will be called ``e + '0'``, ``e + '1'``, ...,
      ``e + 'n-1'``, with ``e`` being said entry), or have ``n`` entries
      (in which case these entries will be used directly as names for the
      generators)

    - ``n`` -- the number of generators, i.e., the rank of the free
      module whose exterior algebra is to be computed (this doesn't have
      to be provided if it can be inferred from the rest of the input)

    REFERENCES:

    - :wikipedia:`Exterior_algebra`
    """
    @staticmethod
    def __classcall_private__(cls, R, names=None, n=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: E1.<e0,e1,e2> = ExteriorAlgebra(QQ)
            sage: E2 = ExteriorAlgebra(QQ, 3)
            sage: E3 = ExteriorAlgebra(QQ, ['e0','e1','e2'])
            sage: E1 is E2 and E2 is E3
            True
        """
        if names is None:
            names = 'e'
        elif names in ZZ:
            n = names
            names = 'e'

        if isinstance(R, FreeModule_generic):
            if n is not None and n != R.dimension():
                raise ValueError("the number of variables does not match the dimension")
            n = R.dimension()
            R = R.base_ring()

        names = tuple(names)
        if n is not None and len(names) != n:
            if len(names) == 1:
                names = tuple('{}{}'.format(names[0], i) for i in range(n))
            else:
                raise ValueError("the number of variables does not match the number of generators")
        return super().__classcall__(cls, R, names)

    def __init__(self, R, names):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.category()
            Category of finite dimensional supercommutative supercocommutative
             super Hopf algebras with basis over Rational Field
            sage: TestSuite(E).run()

            sage: TestSuite(ExteriorAlgebra(GF(3), ['a', 'b'])).run()
        """
        cat = HopfAlgebrasWithBasis(R).FiniteDimensional().Supercommutative().Supercocommutative()
        CliffordAlgebra.__init__(self, QuadraticForm(R, len(names)), names, category=cat)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ExteriorAlgebra(QQ, 3)
            The exterior algebra of rank 3 over Rational Field
        """
        return "The exterior algebra of rank {} over {}".format(self.ngens(), self.base_ring())

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by
        ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E._repr_term((0,1,2))
            'x*y*z'
            sage: y*x + x*z
            -x*y + x*z
        """
        if len(m) == 0:
            return '1'
        term = ''
        for i in m:
            if len(term) != 0:
                term += '*'
            term += self.variable_names()[i]
        return term

    def _ascii_art_term(self, m):
        r"""
        Return ascii art for the basis element indexed by ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E._ascii_art_term((0,1,2))
            x/\y/\z
            sage: ascii_art(y*x + 2*x*z)
            -x/\y + 2*x/\z
        """
        if len(m) == 0:
            return ascii_art('1')
        wedge = '/\\'
        return ascii_art(*[repr(self.basis()[FrozenBitset((i, ))]) for i in m], sep=wedge)

    def _unicode_art_term(self, m):
        """
        Return unicode art for the basis element indexed by ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E._unicode_art_term((0,1,2))
            x∧y∧z
            sage: unicode_art(y*x + x*z)
            -x∧y + x∧z
        """
        if len(m) == 0:
            return unicode_art('1')
        wedge = unicodedata.lookup('LOGICAL AND')
        return unicode_art(*[self.variable_names()[i] for i in m], sep=wedge)

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the basis element indexed
        by ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E._latex_term((0,1,2))
            ' x \\wedge y \\wedge z'
            sage: E.<x0,x1,x2> = ExteriorAlgebra(QQ)
            sage: E._latex_term((0,1,2))
            ' x_{0} \\wedge x_{1} \\wedge x_{2}'
            sage: E._latex_term(())
            '1'
            sage: E._latex_term((0,))
            ' x_{0}'
        """
        if len(m) == 0:
            return '1'
        term = ''
        for i in m:
            if len(term) != 0:
                term += ' \\wedge'
            term += ' ' + self.latex_variable_names()[i]
        return term

    def lift_morphism(self, phi, names=None):
        r"""
        Lift the matrix ``m`` to an algebra morphism of exterior algebras.

        Given a linear map `\phi : V \to W` (here represented by a matrix
        acting on column vectors over the base ring of `V`), this method
        returns the algebra morphism
        `\Lambda(\phi) : \Lambda(V) \to \Lambda(W)`. This morphism is defined
        on generators `v_i \in \Lambda(V)` by `v_i \mapsto \phi(v_i)`.

        .. NOTE::

            This is the map going out of ``self`` as opposed to
            :meth:`~sage.algebras.clifford_algebra.CliffordAlgebraElement.lift_module_morphism()`
            for general Clifford algebras.

        INPUT:

        - ``phi`` -- a linear map `\phi` from `V` to `W`, encoded as a
          matrix
        - ``names`` -- (default: ``'e'``) the names of the generators of
          the Clifford algebra of the domain of (the map represented by)
          ``phi``

        OUTPUT: the algebra morphism `\Lambda(\phi)` from ``self`` to
        `\Lambda(W)`

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: phi = matrix([[0,1],[1,1],[1,2]]); phi
            [0 1]
            [1 1]
            [1 2]
            sage: L = E.lift_morphism(phi, ['a','b','c']); L
            Generic morphism:
              From: The exterior algebra of rank 2 over Rational Field
              To:   The exterior algebra of rank 3 over Rational Field
            sage: L(x)
            b + c
            sage: L(y)
            a + b + 2*c
            sage: L.on_basis()((1,))
            a + b + 2*c
            sage: p = L(E.one()); p
            1
            sage: p.parent()
            The exterior algebra of rank 3 over Rational Field
            sage: L(x*y)
            -a*b - a*c + b*c
            sage: L(x)*L(y)
            -a*b - a*c + b*c
            sage: L(x + y)
            a + 2*b + 3*c
            sage: L(x) + L(y)
            a + 2*b + 3*c
            sage: L(1/2*x + 2)
            1/2*b + 1/2*c + 2
            sage: L(E(3))
            3

            sage: psi = matrix([[1, -3/2]]); psi
            [   1 -3/2]
            sage: Lp = E.lift_morphism(psi, ['a']); Lp
            Generic morphism:
              From: The exterior algebra of rank 2 over Rational Field
              To:   The exterior algebra of rank 1 over Rational Field
            sage: Lp(x)
            a
            sage: Lp(y)
            -3/2*a
            sage: Lp(x + 2*y + 3)
            -2*a + 3

        TESTS:

        Check that the resulting morphism knows it is for
        finite-dimensional algebras (:issue:`25339`)::

            sage: E = ExteriorAlgebra(ZZ, 'e', 3)
            sage: T = jordan_block(0, 2).block_sum(jordan_block(0, 1))
            sage: phi = E.lift_morphism(T)
            sage: phi.category_for()
            Category of finite dimensional super algebras with basis over Integer Ring
            sage: phi.matrix()
            [1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
        """
        n = phi.nrows()
        R = self.base_ring()
        E = ExteriorAlgebra(R, names, n)
        f = lambda x: E.prod(E._from_dict({FrozenBitset((j, )): phi[j, i] for j in range(n)},
                             remove_zeros=True) for i in x)
        cat = AlgebrasWithBasis(R).Super().FiniteDimensional()
        return self.module_morphism(on_basis=f, codomain=E, category=cat)

    def volume_form(self):
        r"""
        Return the volume form of ``self``.

        Given the basis `e_1, e_2, \ldots, e_n` of the underlying
        `R`-module, the volume form is defined as `e_1 \wedge e_2
        \wedge \cdots \wedge e_n`.

        This depends on the choice of basis.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.volume_form()
            x*y*z
        """
        d = self._quadratic_form.dim()
        return self.element_class(self, {tuple(range(d)): self.base_ring().one()})

    def boundary(self, s_coeff):
        r"""
        Return the boundary operator `\partial` defined by the structure
        coefficients ``s_coeff`` of a Lie algebra.

        For more on the boundary operator, see
        :class:`ExteriorAlgebraBoundary`.

        INPUT:

        - ``s_coeff`` -- dictionary whose keys are in `I \times I`, where
          `I` is the index set of the underlying vector space `V`, and whose
          values can be coerced into 1-forms (degree 1 elements) in ``E``
          (usually, these values will just be elements of `V`)

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.boundary({(0,1): z, (1,2): x, (2,0): y})
            Boundary endomorphism of The exterior algebra of rank 3 over Rational Field
        """
        return ExteriorAlgebraBoundary(self, s_coeff)

    def coboundary(self, s_coeff):
        r"""
        Return the coboundary operator `d` defined by the structure
        coefficients ``s_coeff`` of a Lie algebra.

        For more on the coboundary operator, see
        :class:`ExteriorAlgebraCoboundary`.

        INPUT:

        - ``s_coeff`` -- dictionary whose keys are in `I \times I`, where
          `I` is the index set of the underlying vector space `V`, and whose
          values can be coerced into 1-forms (degree 1 elements) in ``E``
          (usually, these values will just be elements of `V`)

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            Coboundary endomorphism of The exterior algebra of rank 3 over Rational Field
        """
        return ExteriorAlgebraCoboundary(self, s_coeff)

    def degree_on_basis(self, m):
        r"""
        Return the degree of the monomial indexed by ``m``.

        The degree of ``m`` in the `\ZZ`-grading of ``self`` is defined
        to be the length of ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.degree_on_basis(())
            0
            sage: E.degree_on_basis((0,))
            1
            sage: E.degree_on_basis((0,1))
            2
        """
        return ZZ(len(m))

    def coproduct_on_basis(self, a):
        r"""
        Return the coproduct on the basis element indexed by ``a``.

        The coproduct is defined by

        .. MATH::

            \Delta(e_{i_1} \wedge \cdots \wedge e_{i_m}) = \sum_{k=0}^m
            \sum_{\sigma \in Ush_{k,m-k}} (-1)^{\sigma}
            (e_{i_{\sigma(1)}} \wedge \cdots \wedge e_{i_{\sigma(k)}}) \otimes
            (e_{i_{\sigma(k+1)}} \wedge \cdots \wedge e_{i_{\sigma(m)}}),

        where `Ush_{k,m-k}` denotes the set of all `(k,m-k)`-unshuffles
        (i.e., permutations in `S_m` which are increasing on the interval
        `\{1, 2, \ldots, k\}` and on the interval
        `\{k+1, k+2, \ldots, k+m\}`).

        .. WARNING::

            This coproduct is a homomorphism of superalgebras, not a
            homomorphism of algebras!

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.coproduct_on_basis((0,))
            1 # x + x # 1
            sage: E.coproduct_on_basis((0,1))
            1 # x*y + x # y - y # x + x*y # 1
            sage: E.coproduct_on_basis((0,1,2))
            1 # x*y*z + x # y*z - y # x*z + x*y # z
             + z # x*y - x*z # y + y*z # x + x*y*z # 1
        """
        from sage.combinat.combinat import unshuffle_iterator
        one = self.base_ring().one()
        L = unshuffle_iterator(tuple(a), one)
        return self.tensor_square()._from_dict(
            {tuple(FrozenBitset(e) if e else FrozenBitset() for e in t): c for t, c in L if c},
            coerce=False,
            remove_zeros=False)

    def antipode_on_basis(self, m):
        r"""
        Return the antipode on the basis element indexed by ``m``.

        Given a basis element `\omega`, the antipode is defined by
        `S(\omega) = (-1)^{\deg(\omega)} \omega`.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.antipode_on_basis(())
            1
            sage: E.antipode_on_basis((1,))
            -y
            sage: E.antipode_on_basis((1,2))
            y*z
        """
        return self.term(m, (-self.base_ring().one())**len(m))

    def counit(self, x):
        r"""
        Return the counit of ``x``.

        The counit of an element `\omega` of the exterior algebra
        is its constant coefficient.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: elt = x*y - 2*x + 3
            sage: E.counit(elt)
            3
        """
        return x.constant_coefficient()

    def interior_product_on_basis(self, a, b):
        r"""
        Return the interior product `\iota_b a` of ``a`` with respect to
        ``b``.

        See :meth:`~sage.algebras.clifford_algebra.CliffordAlgebra.Element.interior_product`
        for more information.

        In this method, ``a`` and ``b`` are supposed to be
        basis elements (see
        :meth:`~sage.algebras.clifford_algebra.CliffordAlgebra.Element.interior_product`
        for a method that computes interior product of arbitrary
        elements), and to be input as their keys.

        This depends on the choice of basis of the vector space
        whose exterior algebra is ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: k = list(E.basis().keys())
            sage: E.interior_product_on_basis(k[1], k[1])
            1
            sage: E.interior_product_on_basis(k[5], k[1])
            z
            sage: E.interior_product_on_basis(k[2], k[5])
            0
            sage: E.interior_product_on_basis(k[5], k[2])
            0
            sage: E.interior_product_on_basis(k[7], k[5])
            -y

        Check :issue:`34694`::

            sage: # needs sage.symbolic
            sage: E = ExteriorAlgebra(SR,'e',3)
            sage: E.inject_variables()
            Defining e0, e1, e2
            sage: a = (e0*e1).interior_product(e0)
            sage: a * e0
            -e0*e1
        """
        sgn = True
        t = list(a)
        for i in b:
            if i not in t:
                return self.zero()
            if t.index(i) % 2:
                sgn = not sgn
            t.remove(i)
        R = self.base_ring()
        if not t:  # catch empty sets
            t = None
        return self.term(FrozenBitset(t), (R.one() if sgn else - R.one()))

    def lifted_bilinear_form(self, M):
        r"""
        Return the bilinear form on the exterior algebra ``self``
        `= \Lambda(V)` which is obtained by lifting the bilinear
        form `f` on `V` given by the matrix ``M``.

        Let `V` be a module over a commutative ring `R`, and let
        `f : V \times V \to R` be a bilinear form on `V`. Then,
        a bilinear form `\Lambda(f) : \Lambda(V) \times
        \Lambda(V) \to R` on `\Lambda(V)` can be canonically
        defined as follows: For every `n \in \NN`, `m \in \NN`,
        `v_1, v_2, \ldots, v_n, w_1, w_2, \ldots, w_m \in V`,
        we define

        .. MATH::

            \Lambda(f)
            ( v_1 \wedge v_2 \wedge \cdots \wedge v_n ,
              w_1 \wedge w_2 \wedge \cdots \wedge w_m )
            := \begin{cases}
              0, &\mbox{if } n \neq m ; \\
              \det G, & \mbox{if } n = m \end{cases} ,

        where `G` is the `n \times m`-matrix whose
        `(i, j)`-th entry is `f(v_i, w_j)`. This bilinear form
        `\Lambda(f)` is known as the bilinear form on
        `\Lambda(V)` obtained by lifting the bilinear form `f`.
        Its restriction to the `1`-st homogeneous component
        `V` of `\Lambda(V)` is `f`.

        The bilinear form `\Lambda(f)` is symmetric if `f` is.

        INPUT:

        - ``M`` -- a matrix over the same base ring as ``self``,
          whose `(i, j)`-th entry is `f(e_i, e_j)`, where
          `(e_1, e_2, \ldots, e_N)` is the standard basis of the
          module `V` for which ``self`` `= \Lambda(V)` (so that
          `N = \dim(V)`), and where `f` is the bilinear form
          which is to be lifted.

        OUTPUT:

        A bivariate function which takes two elements `p` and
        `q` of ``self`` to `\Lambda(f)(p, q)`.

        .. NOTE::

            This takes a bilinear form on `V` as matrix, and
            returns a bilinear form on ``self`` as a function in
            two arguments. We do not return the bilinear form as
            a matrix since this matrix can be huge and one often
            needs just a particular value.

        .. TODO::

            Implement a class for bilinear forms and rewrite this
            method to use that class.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: M = Matrix(QQ, [[1, 2, 3], [2, 3, 4], [3, 4, 5]])
            sage: Eform = E.lifted_bilinear_form(M)
            sage: Eform
            Bilinear Form from The exterior algebra of rank 3 over Rational
            Field (+) The exterior algebra of rank 3 over Rational Field to
            Rational Field
            sage: Eform(x*y, y*z)
            -1
            sage: Eform(x*y, y)
            0
            sage: Eform(x*(y+z), y*z)
            -3
            sage: Eform(x*(y+z), y*(z+x))
            0
            sage: N = Matrix(QQ, [[3, 1, 7], [2, 0, 4], [-1, -3, -1]])
            sage: N.determinant()
            -8
            sage: Eform = E.lifted_bilinear_form(N)
            sage: Eform(x, E.one())
            0
            sage: Eform(x, x*z*y)
            0
            sage: Eform(E.one(), E.one())
            1
            sage: Eform(E.zero(), E.one())
            0
            sage: Eform(x, y)
            1
            sage: Eform(z, y)
            -3
            sage: Eform(x*z, y*z)
            20
            sage: Eform(x+x*y+x*y*z, z+z*y+z*y*x)
            11

        TESTS:

        Exterior algebra over a zero space (a border case)::

            sage: E = ExteriorAlgebra(QQ, 0)
            sage: M = Matrix(QQ, [])
            sage: Eform = E.lifted_bilinear_form(M)
            sage: Eform(E.one(), E.one())
            1
            sage: Eform(E.zero(), E.one())
            0

        .. TODO::

            Another way to compute this bilinear form seems to be to
            map `x` and `y` to the appropriate Clifford algebra and
            there compute `x^t y`, then send the result back to the
            exterior algebra and return its constant coefficient. Or
            something like this. Once the maps to the Clifford and
            back are implemented, check if this is faster.
        """
        R = self.base_ring()

        def lifted_form(x, y):
            result = R.zero()
            for mx, cx in x:
                for my, cy in y:
                    n = len(mx)
                    m = len(my)
                    if m != n:
                        continue
                    matrix_list = [M[i, j] for i in mx for j in my]
                    MA = MatrixArgs(R, n, matrix_list)
                    del matrix_list
                    result += cx * cy * MA.matrix(False).determinant()
            return result
        from sage.categories.cartesian_product import cartesian_product
        return PoorManMap(lifted_form, domain=cartesian_product([self, self]),
                          codomain=self.base_ring(),
                          name="Bilinear Form")

    def _ideal_class_(self, n=0):
        """
        Return the class that is used to implement ideals of ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: type(E.ideal(x*y - z))
            <class 'sage.algebras.clifford_algebra.ExteriorAlgebraIdeal'>

        TESTS::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E._ideal_class_()
            <class 'sage.algebras.clifford_algebra.ExteriorAlgebraIdeal'>
        """
        return ExteriorAlgebraIdeal

    Element = ExteriorAlgebraElement


#####################################################################
# Differentials


class ExteriorAlgebraDifferential(ModuleMorphismByLinearity,
                                  UniqueRepresentation,
                                  metaclass=InheritComparisonClasscallMetaclass):
    r"""
    Internal class to store the data of a boundary or coboundary of
    an exterior algebra `\Lambda(L)` defined by the structure
    coefficients of a Lie algebra `L`.

    See :class:`ExteriorAlgebraBoundary` and
    :class:`ExteriorAlgebraCoboundary` for the actual classes, which
    inherit from this.

    .. WARNING::

        This is not a general class for differentials on the exterior
        algebra.
    """
    @staticmethod
    def __classcall__(cls, E, s_coeff):
        """
        Standardize the structure coefficients to ensure a unique
        representation.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import ExteriorAlgebraDifferential
            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par1 = ExteriorAlgebraDifferential(E, {(0,1): z, (1,2): x, (2,0): y})
            sage: par2 = ExteriorAlgebraDifferential(E, {(0,1): z, (1,2): x, (0,2): -y})
            sage: par3 = ExteriorAlgebraDifferential(E, {(1,0): {2:-1}, (1,2): {0:1}, (2,0):{1:1}})
            sage: par1 is par2
            True
            sage: par1 is par3
            True
            sage: par2 is par3
            True

            sage: par4 = ExteriorAlgebraDifferential(E, {})
            sage: par5 = ExteriorAlgebraDifferential(E, {(1,0): 0, (1,2): {}, (0,2): E.zero()})
            sage: par6 = ExteriorAlgebraDifferential(E, {(1,0): 0, (1,2): 0, (0,2): 0})
            sage: par4 is par5 and par5 is par6
            True
        """
        d = {}

        for k, v in dict(s_coeff).items():
            if not v:  # Strip terms with 0
                continue

            if isinstance(v, dict):
                R = E.base_ring()
                v = E._from_dict({FrozenBitset((i,)): R(c) for i, c in v.items()})
            else:
                # Make sure v is in ``E``
                v = E(v)
                # It's okay if v.degree results in an error
                #   (we'd throw a similar error) unless v == 0 (which
                #   is what v.list() is testing for)
                if v.list() and v.degree() != 1:
                    raise ValueError("elements must be degree 1")

            if k[0] < k[1]:
                d[tuple(k)] = v
            else:
                d[(k[1], k[0])] = -v

        from sage.sets.family import Family
        return super().__classcall__(cls, E, Family(d))

    def __init__(self, E, s_coeff):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2):x, (2,0):y})

        We skip the pickling test as there is an infinite recursion when
        doing equality checks::

            sage: TestSuite(par).run(skip='_test_pickling')

        Check that it knows it is a finite-dimensional algebra
        morphism (:issue:`25339`):;

            sage: par.category_for()
            Category of finite dimensional algebras with basis over Rational Field
            sage: par.matrix()
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0]
            [ 0  0  0  0  0 -1  0  0]
            [ 0  0  0  0  1  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
        """
        self._s_coeff = s_coeff

        # Technically this preserves the grading but with a shift of -1
        cat = AlgebrasWithBasis(E.base_ring()).FiniteDimensional()
        ModuleMorphismByLinearity.__init__(self, domain=E, codomain=E, category=cat)

    def homology(self, deg=None, **kwds):
        """
        Return the homology determined by ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
            sage: par.homology()
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 0 over Rational Field,
             2: Vector space of dimension 0 over Rational Field,
             3: Vector space of dimension 1 over Rational Field}
            sage: d = E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            sage: d.homology()
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 0 over Rational Field,
             2: Vector space of dimension 0 over Rational Field,
             3: Vector space of dimension 1 over Rational Field}
        """
        return self.chain_complex().homology(deg, **kwds)


class ExteriorAlgebraBoundary(ExteriorAlgebraDifferential):
    r"""
    The boundary `\partial` of an exterior algebra `\Lambda(L)` defined
    by the structure coefficients of `L`.

    Let `L` be a Lie algebra. We give the exterior algebra
    `E = \Lambda(L)` a chain complex structure by considering a
    differential `\partial : \Lambda^{k+1}(L) \to \Lambda^k(L)` defined by

    .. MATH::

        \partial(x_1 \wedge x_2 \wedge \cdots \wedge x_{k+1})
        = \sum_{i < j} (-1)^{i+j+1}
        [x_i, x_j] \wedge x_1 \wedge \cdots \wedge \hat{x}_i \wedge \cdots
        \wedge \hat{x}_j \wedge \cdots \wedge x_{k+1}

    where `\hat{x}_i` denotes a missing index. The corresponding homology is
    the Lie algebra homology.

    INPUT:

    - ``E`` -- an exterior algebra of a vector space `L`
    - ``s_coeff`` -- dictionary whose keys are in `I \times I`, where
      `I` is the index set of the basis of the vector space `L`, and whose
      values can be coerced into 1-forms (degree 1 elements) in ``E``;
      this dictionary will be used to define the Lie algebra structure
      on `L` (indeed, the `i`-th coordinate of the Lie bracket of the
      `j`-th and `k`-th basis vectors of `L` for `j < k` is set to be
      the value at the key `(j, k)` if this key appears in ``s_coeff``,
      or otherwise the negated of the value at the key `(k, j)`)

    .. WARNING::

        The values of ``s_coeff`` are supposed to be coercible into
        1-forms in ``E``; but they can also be dictionaries themselves
        (in which case they are interpreted as giving the coordinates of
        vectors in ``L``). In the interest of speed, these dictionaries
        are not sanitized or checked.

    .. WARNING::

        For any two distinct elements `i` and `j` of `I`, the dictionary
        ``s_coeff`` must have only one of the pairs `(i, j)` and
        `(j, i)` as a key. This is not checked.

    EXAMPLES:

    We consider the differential given by Lie algebra given by the cross
    product `\times` of `\RR^3`::

        sage: E.<x,y,z> = ExteriorAlgebra(QQ)
        sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
        sage: par(x)
        0
        sage: par(x*y)
        z
        sage: par(x*y*z)
        0
        sage: par(x+y-y*z+x*y)
        -x + z
        sage: par(E.zero())
        0

    We check that `\partial \circ \partial = 0`::

        sage: p2 = par * par
        sage: all(p2(b) == 0 for b in E.basis())
        True

    Another example: the Lie algebra `\mathfrak{sl}_2`, which has a
    basis `e,f,h` satisfying `[h,e] = 2e`, `[h,f] = -2f`, and `[e,f] = h`::

        sage: E.<e,f,h> = ExteriorAlgebra(QQ)
        sage: par = E.boundary({(0,1): h, (2,1): -2*f, (2,0): 2*e})
        sage: par(E.zero())
        0
        sage: par(e)
        0
        sage: par(e*f)
        h
        sage: par(f*h)
        2*f
        sage: par(h*f)
        -2*f
        sage: C = par.chain_complex(); C
        Chain complex with at most 4 nonzero terms over Rational Field
        sage: ascii_art(C)
                                  [ 0 -2  0]       [0]
                                  [ 0  0  2]       [0]
                    [0 0 0]       [ 1  0  0]       [0]
         0 <-- C_0 <-------- C_1 <----------- C_2 <---- C_3 <-- 0
        sage: C.homology()
        {0: Vector space of dimension 1 over Rational Field,
         1: Vector space of dimension 0 over Rational Field,
         2: Vector space of dimension 0 over Rational Field,
         3: Vector space of dimension 1 over Rational Field}

    Over the integers::

        sage: C = par.chain_complex(R=ZZ); C
        Chain complex with at most 4 nonzero terms over Integer Ring
        sage: ascii_art(C)
                                  [ 0 -2  0]       [0]
                                  [ 0  0  2]       [0]
                    [0 0 0]       [ 1  0  0]       [0]
         0 <-- C_0 <-------- C_1 <----------- C_2 <---- C_3 <-- 0
        sage: C.homology()
        {0: Z, 1: C2 x C2, 2: 0, 3: Z}

    REFERENCES:

    - :wikipedia:`Exterior_algebra#Lie_algebra_homology`
    """
    def _repr_type(self):
        """
        TESTS::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
            sage: par._repr_type()
            'Boundary'
        """
        return "Boundary"

    def _on_basis(self, m):
        """
        Return the differential on the basis element indexed by ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
            sage: par._on_basis(FrozenBitset())
            0
            sage: par._on_basis((0,))
            0
            sage: par._on_basis((0,1))
            z
            sage: par._on_basis((0,2))
            -y
            sage: par._on_basis((0,1,2))
            0
        """
        from itertools import combinations
        E = self.domain()
        sc = self._s_coeff
        keys = sc.keys()

        s = E.zero()

        for b, (i, j) in enumerate(combinations(m, 2)):
            if (i, j) not in keys:
                continue
            t = Bitset(m)
            t.discard(i)
            t.discard(j)
            s += sc[i, j] * E.term(FrozenBitset(t), (-1)**b)

        return s

    @cached_method
    def chain_complex(self, R=None):
        """
        Return the chain complex over ``R`` determined by ``self``.

        INPUT:

        - ``R`` -- the base ring; the default is the base ring of
          the exterior algebra

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
            sage: C = par.chain_complex(); C
            Chain complex with at most 4 nonzero terms over Rational Field
            sage: ascii_art(C)
                                      [ 0  0  1]       [0]
                                      [ 0 -1  0]       [0]
                        [0 0 0]       [ 1  0  0]       [0]
             0 <-- C_0 <-------- C_1 <----------- C_2 <---- C_3 <-- 0

        TESTS:

        This still works in degree `1`::

            sage: E.<x> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({})
            sage: C = par.chain_complex(); C
            Chain complex with at most 2 nonzero terms over Rational Field
            sage: ascii_art(C)
                        [0]
             0 <-- C_0 <---- C_1 <-- 0

        Also in degree `0`::

            sage: E = ExteriorAlgebra(QQ, 0)
            sage: par = E.boundary({})
            sage: C = par.chain_complex(); C
            Chain complex with at most 1 nonzero terms over Rational Field
            sage: ascii_art(C)
             0 <-- C_0 <-- 0
        """
        from sage.homology.chain_complex import ChainComplex
        from sage.matrix.constructor import Matrix
        E = self.domain()
        n = E.ngens()
        if R is None:
            R = E.base_ring()

        if n == 0:
            # Special case because there are no matrices and thus the
            # ChainComplex constructor needs the dimension of the
            # 0th degree space explicitly given.
            return ChainComplex({1: Matrix(R, [[]])}, degree=-1)
            # If you are reading this because you changed something about
            # the ChainComplex constructor and the doctests are failing:
            # This should return a chain complex with degree -1 and
            # only one nontrivial module, namely a free module of rank 1,
            # situated in degree 0.

        # Group the basis into degrees
        basis_by_deg = {deg: [] for deg in range(n+1)}
        for b in E.basis().keys():
            basis_by_deg[len(b)].append(b)

        # Construct the transition matrices
        data = {}
        prev_basis = basis_by_deg[0]
        for deg in range(1, n+1):
            # Make sure within each basis we're sorted by lex
            basis = sorted(basis_by_deg[deg])
            mat = []
            for b in basis:
                ret = self._on_basis(b)
                mat.append([ret.coefficient(p) for p in prev_basis])
            data[deg] = Matrix(mat).transpose().change_ring(R)
            prev_basis = basis

        return ChainComplex(data, degree=-1)


class ExteriorAlgebraCoboundary(ExteriorAlgebraDifferential):
    r"""
    The coboundary `d` of an exterior algebra `\Lambda(L)` defined
    by the structure coefficients of a Lie algebra `L`.

    Let `L` be a Lie algebra. We endow its exterior algebra
    `E = \Lambda(L)` with a cochain complex structure by considering a
    differential `d : \Lambda^k(L) \to \Lambda^{k+1}(L)` defined by

    .. MATH::

        d x_i = \sum_{j < k} s_{jk}^i x_j x_k,

    where `(x_1, x_2, \ldots, x_n)` is a basis of `L`, and where
    `s_{jk}^i` is the `x_i`-coordinate of the Lie bracket `[x_j, x_k]`.

    The corresponding cohomology is the Lie algebra cohomology of `L`.

    This can also be thought of as the exterior derivative, in which case
    the resulting cohomology is the de Rham cohomology of a manifold whose
    exterior algebra of differential forms is ``E``.

    INPUT:

    - ``E`` -- an exterior algebra of a vector space `L`
    - ``s_coeff`` -- dictionary whose keys are in `I \times I`, where
      `I` is the index set of the basis of the vector space `L`, and whose
      values can be coerced into 1-forms (degree 1 elements) in ``E``;
      this dictionary will be used to define the Lie algebra structure
      on `L` (indeed, the `i`-th coordinate of the Lie bracket of the
      `j`-th and `k`-th basis vectors of `L` for `j < k` is set to be
      the value at the key `(j, k)` if this key appears in ``s_coeff``,
      or otherwise the negated of the value at the key `(k, j)`)

    .. WARNING::

        For any two distinct elements `i` and `j` of `I`, the dictionary
        ``s_coeff`` must have only one of the pairs `(i, j)` and
        `(j, i)` as a key. This is not checked.

    EXAMPLES:

    We consider the differential coming from the Lie algebra given by the
    cross product `\times` of `\RR^3`::

        sage: E.<x,y,z> = ExteriorAlgebra(QQ)
        sage: d = E.coboundary({(0,1): z, (1,2): x, (0, 2): -y})
        sage: d(x)
        y*z
        sage: d(y)
        -x*z
        sage: d(x+y-y*z)
        -x*z + y*z
        sage: d(x*y)
        0
        sage: d(E.one())
        0
        sage: d(E.zero())
        0

    We check that `d \circ d = 0`::

        sage: d2 = d * d
        sage: all(d2(b) == 0 for b in E.basis())
        True

    Another example: the Lie algebra `\mathfrak{sl}_2`, which has a
    basis `e,f,h` satisfying `[h,e] = 2e`, `[h,f] = -2f`, and `[e,f] = h`::

        sage: E.<e,f,h> = ExteriorAlgebra(QQ)
        sage: d = E.coboundary({(0,1): h, (2,1): -2*f, (2,0): 2*e})
        sage: d(E.zero())
        0
        sage: d(e)
        -2*e*h
        sage: d(f)
        2*f*h
        sage: d(h)
        e*f
        sage: d(e*f)
        0
        sage: d(f*h)
        0
        sage: d(e*h)
        0
        sage: C = d.chain_complex(); C
        Chain complex with at most 4 nonzero terms over Rational Field
        sage: ascii_art(C)
                                  [ 0  0  1]       [0]
                                  [-2  0  0]       [0]
                    [0 0 0]       [ 0  2  0]       [0]
         0 <-- C_3 <-------- C_2 <----------- C_1 <---- C_0 <-- 0
        sage: C.homology()
        {0: Vector space of dimension 1 over Rational Field,
         1: Vector space of dimension 0 over Rational Field,
         2: Vector space of dimension 0 over Rational Field,
         3: Vector space of dimension 1 over Rational Field}

    Over the integers::

        sage: C = d.chain_complex(R=ZZ); C
        Chain complex with at most 4 nonzero terms over Integer Ring
        sage: ascii_art(C)
                                  [ 0  0  1]       [0]
                                  [-2  0  0]       [0]
                    [0 0 0]       [ 0  2  0]       [0]
         0 <-- C_3 <-------- C_2 <----------- C_1 <---- C_0 <-- 0
        sage: C.homology()
        {0: Z, 1: 0, 2: C2 x C2, 3: Z}

    REFERENCES:

    - :wikipedia:`Exterior_algebra#Differential_geometry`
    """
    def __init__(self, E, s_coeff):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({(0,1): z, (1,2):x, (2,0):y})
            sage: TestSuite(d).run() # known bug - morphisms are properly in a category
        """
        # Construct the dictionary of costructure coefficients, i.e. given
        # [x_j, x_k] = \sum_i s_{jk}^i x_i, we get x^i |-> \sum_{j<k} s_{jk}^i x^j x^k.
        # This dictionary might contain 0 values and might also be missing
        # some keys (both times meaning that the respective `s_{jk}^i` are
        # zero for all `j` and `k`).
        self._cos_coeff = {}
        zero = E.zero()
        B = E.basis()
        for k, v in dict(s_coeff).items():
            if k[0] > k[1]:  # k will have length 2
                k = sorted(k)
                v = -v

            k = B[FrozenBitset(k)]
            for m, c in v:
                self._cos_coeff[m] = self._cos_coeff.get(m, zero) + c * k
        ExteriorAlgebraDifferential.__init__(self, E, s_coeff)

    def _repr_type(self):
        """
        TESTS::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            sage: d._repr_type()
            'Coboundary'
        """
        return "Coboundary"

    def _on_basis(self, m):
        r"""
        Return the differential on the basis element indexed by ``m``.

        EXAMPLES:

        The vector space `\RR^3` made into a Lie algebra using the
        cross product::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({(0,1): z, (1,2): x, (0,2): -y})
            sage: d._on_basis(())
            0
            sage: d._on_basis((0,))
            y*z
            sage: d._on_basis((1,))
            -x*z
            sage: d._on_basis((2,))
            x*y
            sage: d._on_basis((0,1))
            0
            sage: d._on_basis((0,2))
            0
            sage: d._on_basis((0,1,2))
            0
        """
        E = self.domain()
        cc = self._cos_coeff

        tot = E.zero()

        for sgn, i in enumerate(m):
            k = FrozenBitset((i,))
            if k in cc:
                below = tuple([j for j in m if j < i])
                above = tuple([j for j in m if j > i])

                # a hack to deal with empty bitsets
                if not below:
                    below = E.one()
                else:
                    below = E.monomial(FrozenBitset(below))
                if not above:
                    above = E.one()
                else:
                    above = E.monomial(FrozenBitset(above))

                tot += (-1)**sgn * below * cc[k] * above

        return tot

    @cached_method
    def chain_complex(self, R=None):
        """
        Return the chain complex over ``R`` determined by ``self``.

        INPUT:

        - ``R`` -- the base ring; the default is the base ring of
          the exterior algebra

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            sage: C = d.chain_complex(); C
            Chain complex with at most 4 nonzero terms over Rational Field
            sage: ascii_art(C)
                                      [ 0  0  1]       [0]
                                      [ 0 -1  0]       [0]
                        [0 0 0]       [ 1  0  0]       [0]
             0 <-- C_3 <-------- C_2 <----------- C_1 <---- C_0 <-- 0

        TESTS:

        This still works in degree `1`::

            sage: E.<x> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({})
            sage: C = d.chain_complex(); C
            Chain complex with at most 2 nonzero terms over Rational Field
            sage: ascii_art(C)
                        [0]
             0 <-- C_1 <---- C_0 <-- 0

        Also in degree `0`::

            sage: E = ExteriorAlgebra(QQ, 0)
            sage: d = E.coboundary({})
            sage: C = d.chain_complex(); C
            Chain complex with at most 1 nonzero terms over Rational Field
            sage: ascii_art(C)
             0 <-- C_0 <-- 0
        """
        from sage.homology.chain_complex import ChainComplex
        from sage.matrix.constructor import Matrix
        E = self.domain()
        n = E.ngens()
        if R is None:
            R = E.base_ring()

        if n == 0:
            # Special case because there are no matrices and thus the
            # ChainComplex constructor needs the dimension of the
            # 0th degree space explicitly given.
            return ChainComplex({-1: Matrix(R, [[]])}, degree=1)
            # If you are reading this because you changed something about
            # the ChainComplex constructor and the doctests are failing:
            # This should return a chain complex with degree 1 and
            # only one nontrivial module, namely a free module of rank 1,
            # situated in degree 0.

        # Group the basis into degrees
        basis_by_deg = {deg: [] for deg in range(n+1)}
        for b in E.basis().keys():
            basis_by_deg[len(b)].append(b)

        # Construct the transition matrices
        data = {}
        basis = basis_by_deg[0]
        for deg in range(n):
            # Make sure within each basis we're sorted by lex
            next_basis = sorted(basis_by_deg[deg+1])
            mat = []
            for b in basis:
                ret = self._on_basis(b)
                try:
                    mat.append([ret.coefficient(p) for p in next_basis])
                except AttributeError:  # if ret is in E.base_ring()
                    mat.append([E.base_ring()(ret)]*len(next_basis))
            data[deg] = Matrix(mat).transpose().change_ring(R)
            basis = next_basis

        return ChainComplex(data, degree=1)


@richcmp_method
class ExteriorAlgebraIdeal(Ideal_nc):
    """
    An ideal of the exterior algebra.

    EXAMPLES::

        sage: E.<x,y,z> = ExteriorAlgebra(QQ)
        sage: I = E.ideal(x*y); I
        Twosided Ideal (x*y) of The exterior algebra of rank 3 over Rational Field

    We can also use it to build a quotient::

        sage: Q = E.quotient(I); Q
        Quotient of The exterior algebra of rank 3 over Rational Field by the ideal (x*y)
        sage: Q.inject_variables()
        Defining xbar, ybar, zbar
        sage: xbar * ybar
        0
    """
    def __init__(self, ring, gens, coerce=True, side='twosided'):
        """
        Initialize ``self``.

        EXAMPLES:

        We skip the category test because the ideals are not a proper
        element class of the monoid of all ideals::

            sage: E.<y, x> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([x*y - x, x*y - 1])
            sage: TestSuite(I).run(skip='_test_category')

            sage: I = E.ideal([x*y - 3, 0, 2*3])
            sage: TestSuite(I).run(skip='_test_category')

            sage: I = E.ideal([])
            sage: TestSuite(I).run(skip='_test_category')
        """
        self._groebner_strategy = None
        self._reduced = False
        self._homogeneous = all(x.is_super_homogeneous() for x in gens if x)
        if self._homogeneous:
            side = "twosided"
        Ideal_nc.__init__(self, ring, gens, coerce, side)

    def reduce(self, f):
        """
        Reduce ``f`` modulo ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: I = E.ideal(x*y);
            sage: I.reduce(x*y + x*y*z + z)
            z
            sage: I.reduce(x*y + x + y)
            x + y
            sage: I.reduce(x*y + x*y*z)
            0

            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([a+b*c])
            sage: I.reduce(I.gen(0) * d)
            0
        """
        if self._groebner_strategy is None:
            self.groebner_basis()
        R = self.ring()
        return self._groebner_strategy.reduce(R(f))

    def _contains_(self, f):
        r"""
        Return ``True`` if ``f`` is in this ideal,
        ``False`` otherwise.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([x, x*y*z + 2*x*z + 3*y*z], side='left')
            sage: I.groebner_basis()
            (x, y*z)
            sage: x in I
            True
            sage: y*z in I
            True
            sage: x + 3*y*z in I
            True
            sage: x + 3*y in I
            False
            sage: x*y in I
            True
            sage: x + x*y + y*z + x*z in I
            True

        .. NOTE::

            Requires computation of a Groebner basis, which can be a very
            expensive operation.
        """
        return not self.reduce(f)

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([x, x*y*z + 2*x*z + 3*y*z])
            sage: I == I
            True
            sage: Ip = E.ideal([x, y*z])
            sage: Ip == I
            True
            sage: Ip <= I
            True
            sage: Ip < I
            False
            sage: Ip >= I
            True
            sage: Ip > I
            False
            sage: E.ideal([x]) < I
            True
            sage: E.ideal([x]) <= I
            True
            sage: I <= E.ideal([x])
            False

            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: p = a + b*c
            sage: IT = E.ideal([p], side='twosided')
            sage: IR = E.ideal([p], side='right')
            sage: IL = E.ideal([p], side='left')
            sage: IR == IL
            False
            sage: IR <= IL
            False
            sage: IR >= IL
            False
            sage: IL.reduce(p * d)
            2*a*d
            sage: IR.reduce(d * p)
            -2*a*d

            sage: IR <= IT
            True
            sage: IL <= IT
            True
            sage: IT <= IL
            False
            sage: IT <= IR
            False
        """
        if not isinstance(other, ExteriorAlgebraIdeal):
            if op == op_EQ:
                return False
            if op == op_NE:
                return True
            return NotImplemented

        if self is other:
            return rich_to_bool(op, 0)

        # comparison for >= and > : swap the arguments
        if op == op_GE:
            return other.__richcmp__(self, op_LE)
        elif op == op_GT:
            return other.__richcmp__(self, op_LT)

        s_gens = {g for g in self.gens() if g}
        o_gens = {g for g in other.gens() if g}

        if self.side() != other.side():
            if other.side() == "right":
                X = {t * f for t in self.ring().basis() for f in s_gens}
                s_gens.update(X)
            elif other.side() == "left":
                X = {f * t for t in self.ring().basis() for f in s_gens}
                s_gens.update(X)

        if set(s_gens) == set(o_gens):
            return rich_to_bool(op, 0)

        contained = all(f in other for f in s_gens)
        if op == op_LE:
            return contained
        if op == op_NE and not contained:
            return True

        if self.side() != other.side():
            if self.side() == "right":
                X = {t * f for t in self.ring().basis() for f in o_gens}
                s_gens.update(X)
            elif self.side() == "left":
                X = {f * t for t in self.ring().basis() for f in o_gens}
                s_gens.update(X)

        contains = all(f in self for f in o_gens)
        if op == op_EQ:
            return contained and contains
        if op == op_NE:
            return not (contained and contains)
         # remaining case <
        return contained and not contains

    def __mul__(self, other):
        """
        Return the product of ``self`` with ``other``.

        .. WARNING::

            If ``self`` is a right ideal and ``other`` is a left ideal,
            this returns a submodule rather than an ideal.

        EXAMPLES::

            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)

            sage: I = E.ideal([a + 1], side='left')
            sage: J = I * I; J
            Left Ideal (2*a + 1, a, b, c, d, a*b, a*c, a*d, 2*a*b*c + b*c, 2*a*b*d + b*d,
                        2*a*c*d + c*d, a*b*c, a*b*d, a*c*d, b*c*d, a*b*c*d)
             of The exterior algebra of rank 4 over Rational Field
            sage: J.groebner_basis()
            (1,)
            sage: I.gen(0)^2
            2*a + 1

            sage: J = E.ideal([b+c])
            sage: I * J
            Twosided Ideal (a*b + a*c + b + c) of The exterior algebra of rank 4 over Rational Field
            sage: J * I
            Left Ideal (-a*b - a*c + b + c) of The exterior algebra of rank 4 over Rational Field

            sage: K = J * I
            sage: K
            Left Ideal (-a*b - a*c + b + c) of The exterior algebra of rank 4 over Rational Field
            sage: E.ideal([J.gen(0) * d * I.gen(0)], side='left') <= K
            True

            sage: J = E.ideal([b + c*d], side='right')
            sage: I * J
            Twosided Ideal (a*c*d + a*b + c*d + b) of The exterior algebra of rank 4 over Rational Field
            sage: X = J * I; X
            Free module generated by {0, 1, 2, 3, 4, 5, 6, 7} over Rational Field
            sage: [X.lift(b) for b in X.basis()]
            [c*d + b, -a*c*d + a*b, b*c, b*d, a*b*c, a*b*d, b*c*d, a*b*c*d]
            sage: p = X.lift(X.basis()[0])
            sage: p
            c*d + b
            sage: a * p  # not a left ideal
            a*c*d + a*b

            sage: I = E.ideal([a + 1], side='right')
            sage: E.ideal([1]) * I
            Twosided Ideal (a + 1) of The exterior algebra of rank 4 over Rational Field
            sage: I * E.ideal([1])
            Right Ideal (a + 1) of The exterior algebra of rank 4 over Rational Field
        """
        if not isinstance(other, ExteriorAlgebraIdeal) or self.ring() != other.ring():
            return super().__mul__(other)

        if self._homogeneous or other._homogeneous or (self.side() == "left" and other.side() == "right"):
            gens = (x * y for x in self.gens() for y in other.gens())
        else:
            gens = (x * t * y for t in self.ring().basis() for x in self.gens() for y in other.gens())
        gens = [z for z in gens if z]

        if self.side() == "right" and other.side() == "left":
            return self.ring().submodule(gens)

        if self.side() == "left" or self.side() == "twosided":
            if other.side() == "right" or other.side() == "twosided":
                return self.ring().ideal(gens, side='twosided')
            return self.ring().ideal(gens, side='left')
        return self.ring().ideal(gens, side='right')

    def groebner_basis(self, term_order=None, reduced=True):
        r"""
        Return the (reduced) Gröbner basis of ``self``.

        INPUT:

        - ``term_order`` -- the term order used to compute the Gröbner basis;
          must be one of the following:

          * ``'neglex'`` -- (default) negative (read right-to-left) lex order
          * ``'degrevlex'`` -- degree reverse lex order
          * ``'deglex'`` -- degree lex order

        - ``reduced`` -- boolean (default: ``True``); whether or not to return
          the reduced Gröbner basis

        EXAMPLES:

        We compute an example::

            sage: E.<a,b,c,d,e> = ExteriorAlgebra(QQ)
            sage: rels = [c*d*e - b*d*e + b*c*e - b*c*d,
            ....:         c*d*e - a*d*e + a*c*e - a*c*d,
            ....:         b*d*e - a*d*e + a*b*e - a*b*d,
            ....:         b*c*e - a*c*e + a*b*e - a*b*c,
            ....:         b*c*d - a*c*d + a*b*d - a*b*c]
            sage: I = E.ideal(rels)
            sage: I.groebner_basis()
            (-a*b*c + a*b*d - a*c*d + b*c*d,
             -a*b*c + a*b*e - a*c*e + b*c*e,
             -a*b*d + a*b*e - a*d*e + b*d*e,
             -a*c*d + a*c*e - a*d*e + c*d*e)

        With different term orders::

            sage: I.groebner_basis("degrevlex")
            (b*c*d - b*c*e + b*d*e - c*d*e,
             a*c*d - a*c*e + a*d*e - c*d*e,
             a*b*d - a*b*e + a*d*e - b*d*e,
             a*b*c - a*b*e + a*c*e - b*c*e)

            sage: I.groebner_basis("deglex")
            (-a*b*c + a*b*d - a*c*d + b*c*d,
             -a*b*c + a*b*e - a*c*e + b*c*e,
             -a*b*d + a*b*e - a*d*e + b*d*e,
             -a*c*d + a*c*e - a*d*e + c*d*e)

        The example above was computed first using M2, which agrees with
        the ``'degrevlex'`` ordering::

            E = QQ[a..e, SkewCommutative => true]
            I = ideal( c*d*e - b*d*e + b*c*e - b*c*d,
                        c*d*e - a*d*e + a*c*e - a*c*d,
                        b*d*e - a*d*e + a*b*e - a*b*d,
                        b*c*e - a*c*e + a*b*e - a*b*c,
                        b*c*d - a*c*d + a*b*d - a*b*c)
            groebnerBasis(I)

            returns:
            o3 = | bcd-bce+bde-cde acd-ace+ade-cde abd-abe+ade-bde abc-abe+ace-bce |

        By default, the Gröbner basis is reduced, but we can get non-reduced
        Gröber bases (which are not unique)::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([x+y*z])
            sage: I.groebner_basis(reduced=False)
            (x*y, x*z, y*z + x, x*y*z)
            sage: I.groebner_basis(reduced=True)
            (x*y, x*z, y*z + x)

        However, if we have already computed a reduced Gröbner basis (with
        a given term order), then we return that::

            sage: I = E.ideal([x+y*z])  # A fresh ideal
            sage: I.groebner_basis()
            (x*y, x*z, y*z + x)
            sage: I.groebner_basis(reduced=False)
            (x*y, x*z, y*z + x)

        TESTS::

            sage: E.<a,b,c,d,e> = ExteriorAlgebra(ZZ)
            sage: I = E.ideal([a+1, b*c+d])
            sage: I.groebner_basis()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented over fields
        """
        if self.ring().base_ring() not in Fields():
            raise NotImplementedError("only implemented over fields")
        if term_order is None:
            if self._groebner_strategy is not None:
                strategy = type(self._groebner_strategy)
            else:
                from sage.algebras.exterior_algebra_groebner import GroebnerStrategyNegLex as strategy
        else:
            if term_order == "neglex":
                from sage.algebras.exterior_algebra_groebner import GroebnerStrategyNegLex as strategy
            elif term_order == "degrevlex":
                from sage.algebras.exterior_algebra_groebner import GroebnerStrategyDegRevLex as strategy
            elif term_order == "deglex":
                from sage.algebras.exterior_algebra_groebner import GroebnerStrategyDegLex as strategy
            else:
                raise ValueError("invalid term order")
        if isinstance(self._groebner_strategy, strategy):
            if self._reduced or not reduced:
                return self._groebner_strategy.groebner_basis
            self._reduced = reduced
            self._groebner_strategy.reduce_computed_gb()
            return self._groebner_strategy.groebner_basis
        self._groebner_strategy = strategy(self)
        self._groebner_strategy.compute_groebner(reduced=reduced)
        self._reduced = reduced
        return self._groebner_strategy.groebner_basis
