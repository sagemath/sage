r"""
Quotients of Lie algebras

AUTHORS:

- Eero Hakavuori (2018-09-02): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                 https://www.gnu.org/licenses/
# ****************************************************************************

from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.algebras.lie_algebras.subalgebra import LieSubalgebra_finite_dimensional_with_basis
from sage.categories.homset import Hom
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.morphism import SetMorphism
from sage.structure.indexed_generators import standardize_names_index_set


class LieQuotient_finite_dimensional_with_basis(LieAlgebraWithStructureCoefficients):
    r"""
    A quotient Lie algebra.

    INPUT:

    - ``I`` -- an ideal or a list of generators of the ideal
    - ``ambient`` -- (optional) the Lie algebra to be quotiented;
      will be deduced from ``I`` if not given
    - ``names`` -- (optional) a string or a list of strings;
      names for the basis elements of the quotient. If ``names`` is a
      string, the basis will be named ``names_1``,...,``names_n``.

    EXAMPLES:

    The Engel Lie algebra as a quotient of the free nilpotent Lie algebra
    of step 3 with 2 generators::

        sage: L = LieAlgebra(QQ, 2, step=3)
        sage: L.inject_variables()
        Defining X_1, X_2, X_12, X_112, X_122
        sage: I = L.ideal(X_122)
        sage: E = L.quotient(I); E
        Lie algebra quotient L/I of dimension 4 over Rational Field where
        L: Free Nilpotent Lie algebra on 5 generators (X_1, X_2, X_12, X_112, X_122) over Rational Field
        I: Ideal (X_122)
        sage: E.category()
        Join of Category of finite dimensional nilpotent lie algebras with basis
        over Rational Field and Category of subquotients of sets
        sage: E.basis().list()
        [X_1, X_2, X_12, X_112]
        sage: E.inject_variables()
        Defining X_1, X_2, X_12, X_112
        sage: X_1.bracket(X_2)
        X_12
        sage: X_1.bracket(X_12)
        X_112
        sage: X_2.bracket(X_12)
        0

    Shorthand for taking a quotient without creating an ideal first::

        sage: E2 = L.quotient(X_122); E2
        Lie algebra quotient L/I of dimension 4 over Rational Field where
        L: Free Nilpotent Lie algebra on 5 generators (X_1, X_2, X_12, X_112, X_122) over Rational Field
        I: Ideal (X_122)
        sage: E is E2
        True

    Custom names for the basis can be given::

        sage: E.<X,Y,Z,W> = L.quotient(X_122)
        sage: E.basis().list()
        [X, Y, Z, W]
        sage: X.bracket(Z)
        W
        sage: Y.bracket(Z)
        0

    The elements can be relabeled linearly by passing a string to the
    ``names`` parameter::

        sage: E = L.quotient(X_122, names='Y')
        sage: E.basis().list()
        [Y_1, Y_2, Y_3, Y_4]
        sage: E.inject_variables()
        Defining Y_1, Y_2, Y_3, Y_4
        sage: Y_1.bracket(Y_3)
        Y_4
        sage: Y_2.bracket(Y_3)
        0

    Conversion from the ambient Lie algebra uses the quotient projection::

        sage: L = LieAlgebra(QQ, 2, step=3)
        sage: L.inject_variables()
        Defining X_1, X_2, X_12, X_112, X_122
        sage: E = L.quotient(X_122, names='Y')
        sage: E(X_1), E(X_2), E(X_12), E(X_112), E(X_122)
        (Y_1, Y_2, Y_3, Y_4, 0)

    A non-stratifiable Lie algebra as a quotient of the free nilpotent Lie
    algebra of step 4 on 2 generators by the relation
    `[X_2, [X_1, X_2]] = [X_1, [X_1, [X_1, X_2]]]`::

        sage: L = LieAlgebra(QQ, 2, step=4)
        sage: X_1, X_2 = L.homogeneous_component_basis(1)
        sage: rel = L[X_2, [X_1, X_2]] - L[X_1, [X_1, [X_1, X_2]]]
        sage: Q = L.quotient(rel, names='Y')
        sage: Q.dimension()
        5
        sage: Q.inject_variables()
        Defining Y_1, Y_2, Y_3, Y_4, Y_5
        sage: lcs = Q.lower_central_series()
        sage: [I.basis().list() for I in lcs]
        [[Y_1, Y_2, Y_3, Y_4, Y_5], [Y_3, Y_4, Y_5], [Y_4, Y_5], [Y_5], []]
        sage: Y_2.bracket(Y_3)
        -Y_5

    Quotients when the base ring is not a field are not implemented::

        sage: L = lie_algebras.Heisenberg(ZZ, 1)
        sage: L.quotient(L.an_element())
        Traceback (most recent call last):
        ...
        NotImplementedError: quotients over non-fields not implemented

    TESTS:

    Verify that iterated quotient constructions work::

        sage: L = LieAlgebra(QQ, 2, step=3)
        sage: quots = [L]
        sage: for k in range(5):
        ....:     L = L.quotient(L.basis().list()[-1])
        ....:     quots.append(L)
        sage: [Q.dimension() for Q in quots]
        [5, 4, 3, 2, 1, 0]
        sage: all(Lp is Ln.ambient() for Lp, Ln in zip(quots,quots[1:]))
        True
        sage: X = quots[-2].an_element()
        sage: lifts = [X]
        sage: quots = list(reversed(quots[1:-1]))
        sage: for Q in quots:
        ....:     X = Q.lift(X)
        ....:     lifts.append(X)
        sage: all(X.parent() is L for X, L in zip(lifts,quots))
        True

    Verify a quotient construction when the basis ordering and indices ordering
    are different, see :trac:`26352`::

        sage: L.<c,b,a> = LieAlgebra(QQ, abelian=True)
        sage: I2 = L.ideal([a+b, a+c], order=sorted)
        sage: I2.basis()
        Family (b + a, c + a)
        sage: Q = L.quotient(I2)
        sage: Q.basis()
        Finite family {'a': a}

    A test suite::

        sage: L.<x,y,z,w,u> = LieAlgebra(QQ, 2, step=3)
        sage: K = L.quotient(z + w)
        sage: K.dimension()
        2
        sage: TestSuite(K).run()
    """

    @staticmethod
    def __classcall_private__(cls, I, ambient=None, names=None,
                              index_set=None, category=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES:

        Specifying the ambient Lie algebra is not necessary::

            sage: from sage.algebras.lie_algebras.quotient import LieQuotient_finite_dimensional_with_basis
            sage: L.<X,Y> = LieAlgebra(QQ, {('X','Y'): {'X': 1}})
            sage: Q1 = LieQuotient_finite_dimensional_with_basis(X, ambient=L)
            sage: Q2 = LieQuotient_finite_dimensional_with_basis(X)
            sage: Q1 is Q2
            True

        Variable names are extracted from the ambient Lie algebra by default::

            sage: Q3 = L.quotient(X, names=['Y'])
            sage: Q1 is Q3
            True
        """
        if not isinstance(I, LieSubalgebra_finite_dimensional_with_basis):
            # assume I is an element or list of elements of some lie algebra
            if ambient is None:
                if not isinstance(I, (list, tuple)):
                    ambient = I.parent()
                else:
                    ambient = I[0].parent()
            I = ambient.ideal(I)
        if ambient is None:
            ambient = I.ambient()

        if not ambient.base_ring().is_field():
            raise NotImplementedError("quotients over non-fields "
                                      "not implemented")

        # extract an index set from a complementary basis to the ideal
        I_supp = [X.leading_support() for X in I.leading_monomials()]
        inv = ambient.basis().inverse_family()
        sorted_indices = [inv[X] for X in ambient.basis()]
        index_set = [i for i in sorted_indices if i not in I_supp]

        if names is None:
            amb_names = dict(zip(sorted_indices, ambient.variable_names()))
            names = [amb_names[i] for i in index_set]
        elif isinstance(names, str):
            if len(index_set) == 1:
                names = [names]
            else:
                names = ['%s_%d' % (names, k + 1)
                         for k in range(len(index_set))]
        names, index_set = standardize_names_index_set(names, index_set)

        cat = LieAlgebras(ambient.base_ring()).FiniteDimensional().WithBasis()
        if ambient in LieAlgebras(ambient.base_ring()).Nilpotent():
            cat = cat.Nilpotent()
        category = cat.Subquotients().or_subcategory(category)

        sup = super(LieQuotient_finite_dimensional_with_basis, cls)
        return sup.__classcall__(cls, I, ambient, names, index_set,
                                 category=category)

    def __init__(self, I, L, names, index_set, category=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: L.<x,y,z> = LieAlgebra(SR, {('x','y'): {'x':1}})
            sage: K = L.quotient(y)
            sage: K.dimension()
            1
            sage: TestSuite(K).run()
        """
        B = L.basis()
        sm = L.module().submodule_with_basis([I.reduce(B[i]).to_vector()
                                              for i in index_set])
        SB = sm.basis()

        # compute and normalize structural coefficients for the quotient
        s_coeff = {}
        for i, ind_i in enumerate(index_set):
            for j in range(i + 1, len(index_set)):
                ind_j = index_set[j]

                brkt = I.reduce(L.bracket(SB[i], SB[j]))
                brktvec = sm.coordinate_vector(brkt.to_vector())
                s_coeff[(ind_i, ind_j)] = dict(zip(index_set, brktvec))
        s_coeff = LieAlgebraWithStructureCoefficients._standardize_s_coeff(
            s_coeff, index_set)

        self._ambient = L
        self._I = I
        self._sm = sm

        LieAlgebraWithStructureCoefficients.__init__(
            self, L.base_ring(), s_coeff, names, index_set, category=category)

        # register the quotient morphism as a conversion
        H = Hom(L, self)
        f = SetMorphism(H, self.retract)
        self.register_conversion(f)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.gl(QQ, 2)
            sage: a,b,c,d = L.basis()
            sage: Q = L.quotient(d); Q
            Lie algebra quotient L/I of dimension 0 over Rational Field where
            L: General linear Lie algebra of rank 2 over Rational Field
            I: Ideal ([0 0]
            [0 1])
        """
        return ("Lie algebra quotient L/I of dimension %s"
                " over %s where\nL: %s\nI: Ideal %s" % (
                    self.dimension(), self.base_ring(),
                    self.ambient(), self._I._repr_short()))

    def _repr_generator(self, i):
        r"""
        Return the string representation of the basis element
        indexed by ``i`` in ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, 2, step=2)
            sage: Q = L.quotient(x + y)
            sage: Q._repr_generator(Q.indices()[0])
            'x'
        """
        ind = self.indices().index(i)
        return self.variable_names()[ind]

    def ambient(self):
        r"""
        Return the ambient Lie algebra of ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, 2, step=2)
            sage: Q = L.quotient(z)
            sage: Q.ambient() == L
            True
        """
        return self._ambient

    def lift(self, X):
        r"""
        Return some preimage of ``X`` under the quotient projection
        into ``self``.

        INPUT:

        - ``X`` -- an element of ``self``

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, 2, step=2)
            sage: Q = L.quotient(x + y)
            sage: Q(y)
            -x
            sage: el = Q.lift(Q(y)); el
            -x
            sage: el.parent()
            Free Nilpotent Lie algebra on 3 generators (x, y, z) over Rational Field
        """
        L = self.ambient()
        B = L.basis()
        return L.sum(ck * B[ik] for ik, ck in X)

    def retract(self, X):
        r"""
        Map ``X`` under the quotient projection to ``self``.

        INPUT:

        - ``X`` -- an element of the ambient Lie algebra

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, step=2)
            sage: L.inject_variables()
            Defining X_1, X_2, X_3, X_12, X_13, X_23
            sage: Q = L.quotient(X_1 + X_2 + X_3)
            sage: Q.retract(X_1), Q.retract(X_2), Q.retract(X_3)
            (X_1, X_2, -X_1 - X_2)
            sage: all(Q.retract(Q.lift(X)) == X for X in Q.basis())
            True
        """
        X_vec = self._I.reduce(X).to_vector()
        return self.from_vector(self._sm.coordinate_vector(X_vec))

    def defining_ideal(self):
        r"""
        Return the ideal generating this quotient Lie algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1)
            sage: p,q,z = L.basis()
            sage: Q = L.quotient(p)
            sage: Q.defining_ideal()
            Ideal (p1) of Heisenberg algebra of rank 1 over Rational Field
        """
        return self._I

    def from_vector(self, v):
        r"""
        Return the element of ``self`` corresponding to the vector ``v``.

        INPUT:

        - ``v`` -- a vector in ``self.module()`` or ``self.ambient().module()``

        EXAMPLES:

        An element from a vector of the intrinsic module::

            sage: L.<X,Y,Z> = LieAlgebra(QQ, 3, abelian=True)
            sage: Q = L.quotient(X + Y + Z)
            sage: Q.dimension()
            2
            sage: el = Q.from_vector([1, 2]); el
            X + 2*Y
            sage: el.parent() == Q
            True

        An element from a vector of the ambient module

            sage: el = Q.from_vector([1, 2, 3]); el
            -2*X - Y
            sage: el.parent() == Q
            True
        """
        if len(v) == self.ambient().dimension():
            return self.retract(self.ambient().from_vector(v))

        sup = super(LieQuotient_finite_dimensional_with_basis, self)
        return sup.from_vector(v)
