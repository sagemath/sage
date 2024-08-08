# sage.doctest: needs sage.combinat sage.modules
r"""
Graded Lie Conformal Algebras

A (super) Lie conformal algebra `V` is called `H`-graded if there
exists a decomposition `V = \oplus_n V_n` such that the `\lambda`-
bracket is graded of degree `-1`, that is for homogeneous elements
`a \in V_p`, `b \in V_q` with `\lambda`-brackets:

.. MATH::

    [a_\lambda b] = \sum \frac{\lambda^n}{n!} c_n,

we have `c_n \in V_{p+q-n-1}`. This situation arises typically when `V`
has a vector `L \in V` that generates the Virasoro Lie conformal
algebra. Such that for every `a \in V` we have

.. MATH::

    [L_\lambda a] = Ta + \lambda \Delta_a a + O(\lambda^2).

In this situation `V` is graded by the eigenvalues `\Delta_a` of
`L_{(1)}`, the `(1)`-th product with `L`. When the higher order terms
`O(\lambda^2)` vanish we say that `a` is a *primary vector* of
*conformal weight* or degree `\Delta_a`.

.. NOTE::

    Although arbitrary gradings are allowed, many of the constructions
    we implement in these classes work only for positive rational
    gradings.

AUTHORS:

- Reimundo Heluani (2019-08-09): Initial implementation.
"""


#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from .lie_conformal_algebra_with_structure_coefs import \
                                LieConformalAlgebraWithStructureCoefficients


class GradedLieConformalAlgebra(LieConformalAlgebraWithStructureCoefficients):
    r"""
    An H-Graded Lie conformal algebra.

    INPUT:

    - ``R`` -- a commutative ring (default: ``None``); the base
      ring of this Lie conformal algebra. Behaviour is undefined if
      it is not a field of characteristic zero

    - ``s_coeff`` -- dictionary (default: ``None``); as in the
      input of :class:`LieConformalAlgebra`

    - ``names`` -- tuple of strings (default: ``None``); as in the
      input of :class:`LieConformalAlgebra`

    - ``central_elements`` -- tuple of strings (default: ``None``);
      as in the input of :class:`LieConformalAlgebra`

    - ``index_set`` -- enumerated set (default: ``None``); as in the
      input of :class:`LieConformalAlgebra`

    - ``weights`` -- tuple of nonnegative rational numbers
      (default: tuple of ``1``); a list of degrees for this Lie
      conformal algebra.
      This tuple needs to have the same cardinality as
      ``index_set`` or ``names``. Central elements are assumed
      to have weight ``0``.

    - ``category`` -- the category that this Lie conformal algebra
      belongs to

    - ``parity`` -- tuple of ``0`` or ``1`` (default: tuple of
      ``0``); a tuple specifying the parity of each non-central
      generator

    EXAMPLES::

        sage: bosondict = {('a','a'):{1:{('K',0):1}}}
        sage: R = LieConformalAlgebra(QQ,bosondict,names=('a',),central_elements=('K',), weights=(1,))
        sage: R.inject_variables()
        Defining a, K
        sage: a.T(3).degree()
        4
        sage: K.degree()
        0
        sage: R.category()
        Category of H-graded finitely generated Lie conformal algebras with basis over Rational Field
    """
    def __init__(self, R, s_coeff, index_set=None, central_elements=None,
                 category=None, prefix=None, names=None, latex_names=None,
                 parity=None, weights=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: V = lie_conformal_algebras.Virasoro(QQ)
            sage: TestSuite(V).run()
        """
        is_super = kwds.get('super', None)
        default_category = LieConformalAlgebras(R).WithBasis().FinitelyGenerated().Graded()
        if is_super or parity:
            category = default_category.Super().or_subcategory(category)
        else:
            category = default_category.or_subcategory(category)

        LieConformalAlgebraWithStructureCoefficients.__init__(self, R,
            s_coeff, index_set=index_set, central_elements=central_elements,
            category=category, prefix=prefix,
            names=names, latex_names=latex_names, parity=parity, **kwds)

        if weights is None:
            weights = (1,) * (len(self._generators) -
                              len(self.central_elements()))
        if len(weights) != (len(self._generators) -
                            len(self.central_elements())):
            raise ValueError("weights and (non-central) generator lists "
                             "must be of same length")
        self._weights = weights
