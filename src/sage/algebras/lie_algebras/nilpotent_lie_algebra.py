r"""
Nilpotent Lie algebras

AUTHORS:

- Eero Hakavuori (2018-08-16): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.categories.lie_algebras import LieAlgebras
from sage.structure.indexed_generators import standardize_names_index_set
from sage.rings.integer_ring import ZZ


class NilpotentLieAlgebra_dense(LieAlgebraWithStructureCoefficients):
    r"""
    A nilpotent Lie algebra `L` over a base ring.

    INPUT:

    - ``R`` -- the base ring
    - ``s_coeff`` -- a dictionary of structural coefficients
    - ``names`` -- (default:``None``) list of strings to use as names of basis
      elements; if ``None``, the names will be inferred from the structural
      coefficients
    - ``index_set`` -- (default:``None``) list of hashable and comparable
      elements to use for indexing
    - ``step`` -- (optional) an integer; the nilpotency step of the
      Lie algebra if known; otherwise it will be computed when needed
    - ``category`` -- (optional) a subcategory of finite dimensional
      nilpotent Lie algebras with basis

    EXAMPLES:

    The input to a :class:`NilpotentLieAlgebra_dense` should be of the
    same form as to a
    :class:`~sage.algebras.lie_algebras.structure_coefficients.LieAlgebraWithStructureCoefficients`::

        sage: L.<X,Y,Z,W> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
        sage: L
        Nilpotent Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: L[X, Y]
        Z
        sage: L[X, W]
        0

    If the parameter ``names`` is omitted, then the terms appearing in the
    structural coefficients are used as names::

        sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True); L
        Nilpotent Lie algebra on 3 generators (X, Y, Z) over Rational Field

    TESTS::

        sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1},
        ....:                     ('X','Z'): {'W': 1},
        ....:                     ('Y','Z'): {'T': 1}}, nilpotent=True)
        sage: TestSuite(L).run()
    """

    @staticmethod
    def __classcall_private__(cls, R, s_coeff, names=None, index_set=None,
                              category=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES:

        If the variable order is specified, the order of structural
        coefficients does not matter::

            sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra_dense
            sage: L1.<x,y,z> = NilpotentLieAlgebra_dense(QQ, {('x','y'): {'z': 1}})
            sage: L2.<x,y,z> = NilpotentLieAlgebra_dense(QQ, {('y','x'): {'z': -1}})
            sage: L1 is L2
            True

        If the variables are implicitly defined by the structural coefficients,
        the ordering may be different and the Lie algebras will be considered
        different::

            sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra_dense
            sage: L1 = NilpotentLieAlgebra_dense(QQ, {('x','y'): {'z': 1}})
            sage: L2 = NilpotentLieAlgebra_dense(QQ, {('y','x'): {'z': -1}})
            sage: L1
            Nilpotent Lie algebra on 3 generators (x, y, z) over Rational Field
            sage: L2
            Nilpotent Lie algebra on 3 generators (y, x, z) over Rational Field
            sage: L1 is L2
            False

        Constructed using two different methods from :class:`LieAlgebra`
        yields the same Lie algebra::

            sage: sc = {('X','Y'): {'Z': 1}}
            sage: C = LieAlgebras(QQ).Nilpotent().FiniteDimensional().WithBasis()
            sage: L1.<X,Y,Z> = LieAlgebra(QQ, sc, category=C)
            sage: L2 = LieAlgebra(QQ, sc, nilpotent=True, names=['X','Y','Z'])
            sage: L1 is L2
            True
        """
        if not names:
            # extract names from structural coefficients
            names = []
            for (X, Y), d in s_coeff.items():
                if X not in names:
                    names.append(X)
                if Y not in names:
                    names.append(Y)
                for k in d:
                    if k not in names:
                        names.append(k)

        from sage.structure.indexed_generators import standardize_names_index_set
        names, index_set = standardize_names_index_set(names, index_set)
        s_coeff = LieAlgebraWithStructureCoefficients._standardize_s_coeff(
            s_coeff, index_set)

        cat = LieAlgebras(R).FiniteDimensional().WithBasis().Nilpotent()
        category = cat.or_subcategory(category)

        return super(NilpotentLieAlgebra_dense, cls).__classcall__(
            cls, R, s_coeff, names, index_set, category=category, **kwds)

    def __init__(self, R, s_coeff, names, index_set, step=None, **kwds):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L.<X,Y,Z,W> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
            sage: TestSuite(L).run()
        """
        if step is not None:
            self._step = step

        LieAlgebraWithStructureCoefficients.__init__(self, R, s_coeff,
                                                     names, index_set,
                                                     **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L.<X,Y,Z,W> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
            sage: L
            Nilpotent Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        """
        return "Nilpotent %s" % (super(NilpotentLieAlgebra_dense, self)._repr_())


class FreeNilpotentLieAlgebra(NilpotentLieAlgebra_dense):
    r"""
    Return the free nilpotent Lie algebra of step ``s`` with ``r`` generators.

    The free nilpotent Lie algebra `L` of step `s` with `r` generators is
    the quotient of the free Lie algebra on `r` generators by the `(s+1)`-th
    term of the lower central series. That is, the only relations in the
    Lie algebra `L` are anticommutativity, the Jacobi identity, and the
    vanishing of all brackets of length more than `s`.

    INPUT:

    - ``R`` -- the base ring
    - ``r`` -- an integer; the number of generators
    - ``s`` -- an integer; the nilpotency step of the algebra
    - ``names`` -- (optional) a string or a list of strings used to name the
      basis elements; if ``names`` is a string, then names for the basis
      will be autogenerated as determined by the ``naming`` parameter
    - ``naming`` -- (optional) a string; the naming scheme to use for
      the basis; valid values are:

      * ``'index'`` - (default for `r < 10`) the basis elements are
        ``names_w``, where ``w`` are Lyndon words indexing the basis
      * ``'linear'`` - (default for `r \geq 10`) the basis is indexed
        ``names_1``, ..., ``names_n`` in the ordering of the Lyndon basis

    .. NOTE::

        The ``'index'`` naming scheme is not supported if `r \geq 10`
        since it leads to ambiguous names.

    EXAMPLES:

    We compute the free step 4 Lie algebra on 2 generators and
    verify the only non-trivial relation
    `[[X_1,[X_1,X_2]],X_2] = [X_1,[[X_1,X_2],X_2]]`::

        sage: L = LieAlgebra(QQ, 2, step=4)
        sage: L.basis().list()
        [X_1, X_2, X_12, X_112, X_122, X_1112, X_1122, X_1222]
        sage: X_1, X_2 = L.basis().list()[:2]
        sage: L[[X_1, [X_1, X_2]], X_2]
        X_1122
        sage: L[[X_1, [X_1, X_2]], X_2] == L[X_1, [[X_1, X_2], X_2]]
        True

    The linear naming scheme on the same Lie algebra::

        sage: K = LieAlgebra(QQ, 2, step=4, names='Y', naming='linear')
        sage: K.basis().list()
        [Y_1, Y_2, Y_3, Y_4, Y_5, Y_6, Y_7, Y_8]
        sage: K.inject_variables()
        Defining Y_1, Y_2, Y_3, Y_4, Y_5, Y_6, Y_7, Y_8
        sage: Y_2.bracket(Y_3)
        -Y_5
        sage: Y_5.bracket(Y_1)
        -Y_7
        sage: Y_3.bracket(Y_4)
        0

    A fully custom naming scheme on the Heisenberg algebra::

        sage: L = LieAlgebra(ZZ, 2, step=2, names=('X', 'Y', 'Z'))
        sage: a, b, c = L.basis()
        sage: L.basis().list()
        [X, Y, Z]
        sage: a.bracket(b)
        Z

    An equivalent way to define custom names for the basis elements and
    bind them as local variables simultaneously::

        sage: L.<X,Y,Z> = LieAlgebra(ZZ, 2, step=2)
        sage: L.basis().list()
        [X, Y, Z]
        sage: X.bracket(Y)
        Z

    A free nilpotent Lie algebra is a stratified nilpotent Lie algebra::

        sage: L = LieAlgebra(QQ, 3, step=3)
        sage: L.category()
        Category of finite dimensional stratified lie algebras with basis over Rational Field
        sage: L in LieAlgebras(QQ).Nilpotent()
        True

    Being graded means that each basis element has a degree::

        sage: L in LieAlgebras(QQ).Graded()
        True
        sage: L.homogeneous_component_basis(1).list()
        [X_1, X_2, X_3]
        sage: L.homogeneous_component_basis(2).list()
        [X_12, X_13, X_23]
        sage: L.homogeneous_component_basis(3).list()
        [X_112, X_113, X_122, X_123, X_132, X_133, X_223, X_233]

    TESTS:

    Verify all bracket relations in the free nilpotent Lie algebra of step 5
    with 2 generators::

        sage: L = LieAlgebra(QQ, 2, step=5)
        sage: L.inject_variables()
        Defining X_1, X_2, X_12, X_112, X_122, X_1112, X_1122, X_1222,
        X_11112, X_11122, X_11212, X_11222, X_12122, X_12222
        sage: [X_1.bracket(Xk) for Xk in L.basis()]
        [0, X_12, X_112, X_1112, X_1122,
        X_11112, X_11122, X_11222, 0, 0, 0, 0, 0, 0]
        sage: [X_2.bracket(Xk) for Xk in L.basis()]
        [-X_12, 0, -X_122, -X_1122, -X_1222,
        -X_11122 + X_11212, -X_11222 - X_12122, -X_12222, 0, 0, 0, 0, 0, 0]
        sage: [X_12.bracket(Xk) for Xk in L.basis()]
        [-X_112, X_122, 0, -X_11212, X_12122, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_112.bracket(Xk) for Xk in L.basis()]
        [-X_1112, X_1122, X_11212, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_122.bracket(Xk) for Xk in L.basis()]
        [-X_1122, X_1222, -X_12122, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_1112.bracket(Xk) for Xk in L.basis()]
        [-X_11112, X_11122 - X_11212, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_1122.bracket(Xk) for Xk in L.basis()]
        [-X_11122, X_11222 + X_12122, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_1222.bracket(Xk) for Xk in L.basis()]
        [-X_11222, X_12222, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_11112.bracket(Xk) for Xk in L.basis()]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_11122.bracket(Xk) for Xk in L.basis()]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_11212.bracket(Xk) for Xk in L.basis()]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_11222.bracket(Xk) for Xk in L.basis()]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_12122.bracket(Xk) for Xk in L.basis()]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [X_12222.bracket(Xk) for Xk in L.basis()]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    The dimensions of the smallest free nilpotent Lie algebras on
    2 and 3 generators::

        sage: l = [LieAlgebra(QQ, 2, step=k) for k in range(1, 7)]
        sage: [L.dimension() for L in l]
        [2, 3, 5, 8, 14, 23]
        sage: l = [LieAlgebra(QQ, 3, step=k) for k in range(1, 4)]
        sage: [L.dimension() for L in l]
        [3, 6, 14]

    Verify that a free nilpotent Lie algebra of step `>2` with `>10`
    generators can be created, see :trac:`27018` (see also :trac:`27069`)::

        sage: L = LieAlgebra(QQ, 11, step=3)
        sage: L.dimension() == 11 + (11^2-11)/2 + (11^3-11)/3
        True
    """
    @staticmethod
    def __classcall_private__(cls, R, r, s, names=None, naming=None, category=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: L1.<X> = LieAlgebra(ZZ, 2, step=2)
            sage: L2 = LieAlgebra(ZZ, 2, step=2, names=['X'])
            sage: L3 = LieAlgebra(ZZ, 2, step=2)
            sage: L1 is L2 and L2 is L3
            True
        """
        if names is None:
            names = 'X'

        cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        category = cat.Graded().Stratified().or_subcategory(category)

        return super(FreeNilpotentLieAlgebra, cls).__classcall__(
            cls, R, r, s, names=tuple(names), naming=naming,
            category=category, **kwds)

    def __init__(self, R, r, s, names, naming, category, **kwds):
        r"""
        Initialize ``self``

        EXAMPLES::

            sage: L = LieAlgebra(ZZ, 2, step=2)
            sage: TestSuite(L).run()

            sage: L = LieAlgebra(QQ, 4, step=3)
            sage: TestSuite(L).run()  # long time
        """
        if r not in ZZ or r <= 0:
            raise ValueError("number of generators %s is not "
                             "a positive integer" % r)
        if s not in ZZ or s <= 0:
            raise ValueError("step %s is not a positive integer" % s)

        # extract an index set from the Lyndon words of the corresponding
        # free Lie algebra, and store the corresponding elements in a dict
        from sage.algebras.lie_algebras.lie_algebra import LieAlgebra

        free_gen_names = ['F%d' % k for k in range(r)]
        free_gen_names_inv = {val: i + 1 for i, val in enumerate(free_gen_names)}
        L = LieAlgebra(R, free_gen_names).Lyndon()

        basis_by_deg = {d: [] for d in range(1, s + 1)}
        for d in range(1, s + 1):
            for X in L.graded_basis(d):
                # convert brackets of form [X_1, [X_1, X_2]] to words (1,1,2)
                w = tuple(free_gen_names_inv[s]
                          for s in X.leading_support().to_word())
                basis_by_deg[d].append((w, X))

        index_set = [ind for d in basis_by_deg for ind, val in basis_by_deg[d]]

        if len(names) == 1 and len(index_set) > 1:
            if not naming:
                if r >= 10:
                    naming = 'linear'
                else:
                    naming = 'index'
            if naming == 'linear':
                names = ['%s_%d' % (names[0], k + 1)
                         for k in range(len(index_set))]
            elif naming == 'index':
                if r >= 10:
                    raise ValueError("'index' naming scheme not supported for "
                                     "10 or more generators")
                names = ['%s_%s' % (names[0], "".join(str(s) for s in ind))
                         for ind in index_set]
            else:
                raise ValueError("unknown naming scheme %s" % naming)

        # extract structural coefficients from the free Lie algebra
        s_coeff = {}
        for dx in range(1, s + 1):
            # Brackets are only computed when deg(X) + deg(Y) <= s
            # We also require deg(Y) >= deg(X) by the ordering
            for dy in range(dx, s + 1 - dx):
                if dx == dy:
                    for i, val in enumerate(basis_by_deg[dx]):
                        X_ind, X = val
                        for Y_ind, Y in basis_by_deg[dy][i + 1:]:
                            Z = L[X, Y]
                            if not Z.is_zero():
                                s_coeff[(X_ind, Y_ind)] = {W_ind: Z[W.leading_support()]
                                                           for W_ind, W in basis_by_deg[dx + dy]}
                else:
                    for X_ind, X in basis_by_deg[dx]:
                        for Y_ind, Y in basis_by_deg[dy]:
                            Z = L[X, Y]
                            if not Z.is_zero():
                                s_coeff[(X_ind, Y_ind)] = {W_ind: Z[W.leading_support()]
                                                           for W_ind, W in basis_by_deg[dx + dy]}

        names, index_set = standardize_names_index_set(names, index_set)
        s_coeff = LieAlgebraWithStructureCoefficients._standardize_s_coeff(
            s_coeff, index_set)

        NilpotentLieAlgebra_dense.__init__(self, R, s_coeff, names,
                                           index_set, s,
                                           category=category, **kwds)

    def _repr_generator(self, w):
        r"""
        Return the string representation of the basis element
        indexed by the word ``w`` in ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=4)
            sage: L._repr_generator((1, 1, 2, 2))
            'X_1122'
            sage: L = LieAlgebra(QQ, 2, step=4, naming='linear')
            sage: L._repr_generator((1, 1, 2, 2))
            'X_7'
            sage: L.<X,Y,Z> = LieAlgebra(QQ, 2, step=2)
            sage: L._repr_generator((1, 2))
            'Z'
        """
        i = self.indices().index(w)
        return self.variable_names()[i]

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=3)
            sage: L
            Free Nilpotent Lie algebra on 5 generators (X_1, X_2, X_12, X_112, X_122) over Rational Field
        """
        return "Free %s" % (super(FreeNilpotentLieAlgebra, self)._repr_())

