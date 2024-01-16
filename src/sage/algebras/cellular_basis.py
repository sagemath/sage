# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Cellular Basis
==============

Cellular algebras are a class of algebras introduced by Graham and Lehrer
[GrLe1996]_. The :class:`CellularBasis` class provides a general framework for
implementing cellular algebras and their cell modules and simple modules.

Let `R` be a commutative ring. A `R`-algebra `A` is a *cellular algebra*
if it has a *cell datum*, which is a tuple `(\Lambda, i, M, C)`,
where `\Lambda` is finite poset with order `\ge`, if `\mu \in \Lambda`
then `T(\mu)` is a finite set and

.. MATH::

    C \colon \coprod_{\mu \in \Lambda} T(\mu) \times T(\mu)
      \longrightarrow A; (\mu,s,t) \mapsto c^\mu_{st}
      \text{ is an injective map}

such that the following holds:

* The set `\{c^\mu_{st}\mid \mu\in\Lambda, s,t\in T(\mu)\}` is a basis of `A`.
* If `a \in A` and `\mu\in\Lambda, s,t \in T(\mu)` then:

  .. MATH::

      a c^\mu_{st} = \sum_{u\in T(\mu)} r_a(s,u) c^\mu_{ut} \pmod{A^{>\mu}},

  where `A^{>\mu}` is spanned by

  .. MATH::

      `\{c^\nu_{ab} | \nu > \mu\text{ and } a,b \in T(\nu)\}`.

  Moreover, the scalar `r_a(s,u)` depends only on `a`, `s` and
  `u` and, in particular, is independent of `t`.

* The map `\iota \colon A \longrightarrow A; c^\mu_{st} \mapsto c^\mu_{ts}`
  is an algebra anti-isomorphism.

A *cellular  basis* for `A` is any basis of the form
`\{c^\mu_{st} \mid \mu \in \Lambda, s,t \in T(\mu) \}`.

Note that the scalars `r_a(s,u) \in R` depend only if `a`, `s` and `u` and,
in particular, they do not depend on `t`. It follows from the definition of
a cell datum that `A^{>\mu}` is a two-sided ideal of `A`. More importantly,
if `\mu \in \Lambda` then the
:class:`~sage.modules.with_basis.cell_module.CellModule` `C^\mu` is the
free `R`-module with basis `\{c^\mu_s \mid \mu \in \Lambda, s \in T(\mu)\}`
and with `A`-action:

.. MATH::

    a c^\mu_{s} = \sum_{u \in T(\mu)} r_a(s,u) c^\mu_{u},

where the scalars `r_a(s,u)` are those appearing in the definition of the cell
datum. It follows from the cellular basis axioms that `C^\mu` comes
equipped with a bilinear form `\langle\ ,\ \rangle` that is determined by:

.. MATH::

    c^\mu_{st} c^\mu_u = \langle c^\mu_{s}, c^\mu_t \rangle c^\mu_u.

The *radical* of `C^\mu` is the `A`-submodule
`\operatorname{rad} C^\mu = \{x \in C^\mu | \langle x,y \rangle = 0 \}`.
Hence, `D^\mu = C^\mu / \operatorname{rad} C^\mu` is also an `A`-module.
It is not difficult to show that `\{ D^\mu \mid D^\mu \neq 0 \}` is a
complete set of pairwise non-isomorphic `A`-modules. Hence, a cell datum for
`A` gives an explicit construction of the irreducible `A`-modules. The module
:meth:`~sage.modules.with_basis.cell_module.CellModule.simple_module()`
`D^\mu` is either zero or absolutely irreducible.

EXAMPLES:

We compute a cellular basis and do some basic computations::

    sage: S = SymmetricGroupAlgebra(QQ, 3)
    sage: C = S.cellular_basis()
    sage: C
    Cellular basis of Symmetric group algebra of order 3
     over Rational Field

.. SEEALSO::

    :class:`~sage.modules.with_basis.cell_module.CellModule`

AUTHOR:

- Travis Scrimshaw (2015-11-5): Initial version

REFERENCES:

- [GrLe1996]_
- [KX1998]_
- [Mat1999]_
- :wikipedia:`Cellular_algebra`
- http://webusers.imj-prg.fr/~bernhard.keller/ictp2006/lecturenotes/xi.pdf
"""

#*****************************************************************************
#       Copyright (C) 2015-2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.algebras import Algebras


class CellularBasis(CombinatorialFreeModule):
    r"""
    The cellular basis of a cellular algebra, in the sense of
    Graham and Lehrer [GrLe1996]_.

    INPUT:

    - ``A`` -- the cellular algebra

    EXAMPLES:

    We compute a cellular basis and do some basic computations::

        sage: S = SymmetricGroupAlgebra(QQ, 3)
        sage: C = S.cellular_basis()
        sage: C
        Cellular basis of Symmetric group algebra of order 3
         over Rational Field
        sage: len(C.basis())
        6
        sage: len(S.basis())
        6
        sage: a,b,c,d,e,f = C.basis()
        sage: f
        C([3], [[1, 2, 3]], [[1, 2, 3]])
        sage: c
        C([2, 1], [[1, 3], [2]], [[1, 2], [3]])
        sage: d
        C([2, 1], [[1, 2], [3]], [[1, 3], [2]])
        sage: f * f
        C([3], [[1, 2, 3]], [[1, 2, 3]])
        sage: f * c
        0
        sage: d * c
        C([2, 1], [[1, 2], [3]], [[1, 2], [3]])
        sage: c * d
        C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
        sage: S(f)
        1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1]
         + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
        sage: S(d)
        1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1]
        sage: B = list(S.basis())
        sage: B[2]
        [2, 1, 3]
        sage: C(B[2])
        -C([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]])
         + C([2, 1], [[1, 2], [3]], [[1, 2], [3]])
         - C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
         + C([3], [[1, 2, 3]], [[1, 2, 3]])
    """
    def __init__(self, A, to_algebra=None, from_algebra=None, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: TestSuite(C).run()
        """
        self._algebra = A
        I = [(la, s, t) for la in A.cell_poset()
             for s in A.cell_module_indices(la)
             for t in A.cell_module_indices(la)]

        # TODO: Use instead A.category().Realizations() so
        #   operations are defined by coercion?
        prefix = kwargs.pop('prefix', 'C')
        cat = Algebras(A.category().base_ring()).FiniteDimensional().WithBasis().Cellular()
        CombinatorialFreeModule.__init__(self, A.base_ring(), I,
                                         prefix=prefix, bracket=False,
                                         category=cat, **kwargs)

        # Register coercions
        if from_algebra is None:
            from_algebra = A._to_cellular_element
        if to_algebra is None:
            to_algebra = A._from_cellular_index
        if from_algebra is not NotImplemented:
            to_cellular = A.module_morphism(from_algebra, codomain=self,
                                            category=cat)
        if to_algebra is NotImplemented:
            from_cellular = ~to_cellular
        else:
            from_cellular = self.module_morphism(to_algebra, codomain=A,
                                                 category=cat)
            if from_algebra is NotImplemented:
                to_cellular = ~from_cellular
        to_cellular.register_as_coercion()
        from_cellular.register_as_coercion()

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: S.cellular_basis()
            Cellular basis of Symmetric group algebra of order 3 over Rational Field
        """
        return "Cellular basis of {}".format(self._algebra)

    def _latex_term(self, x):
        r"""
        Return a latex representation of the term indexed by ``x``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: s = Tableau([[1,2],[3]])
            sage: C._latex_term((Partition([2,1]), s, s))
            'C^{...}_{\\left(...\\right)}'
        """
        from sage.misc.latex import latex
        la = x[0]
        m = (x[1], x[2])
        # m contains "non-LaTeXed" strings, use string representation
        sla = latex(la)
        if sla.find('\\text{\\textt') != -1:
            sla = str(la)
        sm = latex(m)
        if sm.find('\\text{\\textt') != -1:
            sm = str(m)
        return "C^{%s}_{%s}" % (sla, sm)

    def cellular_basis_of(self):
        """
        Return the defining algebra of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.cellular_basis_of() is S
            True
        """
        return self._algebra

    def cell_poset(self):
        """
        Return the cell poset of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.cell_poset()
            Finite poset containing 3 elements
        """
        return self._algebra.cell_poset()

    def cell_module_indices(self, la):
        r"""
        Return the indices of the cell module of ``self``
        indexed by ``la`` .

        This is the finite set `M(\lambda)`.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.cell_module_indices([2,1])
            Standard tableaux of shape [2, 1]
        """
        return self._algebra.cell_module_indices(la)

    def cellular_basis(self):
        """
        Return the cellular basis of ``self``, which is ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.cellular_basis() is C
            True
        """
        return self

    @cached_method
    def one(self):
        """
        Return the element `1` in ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.one()
            C([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]])
             + C([2, 1], [[1, 2], [3]], [[1, 2], [3]])
             + C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
             + C([3], [[1, 2, 3]], [[1, 2, 3]])
        """
        return self(self._algebra.one())

    @cached_method
    def product_on_basis(self, x, y):
        """
        Return the product of basis indices by ``x`` and ``y``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: la = Partition([2,1])
            sage: s = StandardTableau([[1,2],[3]])
            sage: t = StandardTableau([[1,3],[2]])
            sage: C.product_on_basis((la, s, t), (la, s, t))
            0
        """
        A = self._algebra
        return self(A(self.monomial(x)) * A(self.monomial(y)))
